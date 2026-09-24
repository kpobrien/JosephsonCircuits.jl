# Solve benchmark: the hard cases the solvers are judged on. Run in a fresh
# process with one Julia/BLAS thread:
# julia --startup-file=no --project=. benchmark/solve.jl label
# Building the circuits and problems is outside the timers, and the first
# call of each case includes its compilation. Each row ends with what the
# solve reached, since the time of a solve which did not converge measures
# nothing.
using JosephsonCircuits, Statistics
include(joinpath(@__DIR__, "options.jl"))
const label = BENCH_LABEL
# the solves take seconds, so fewer samples than the other scripts take
const SAMPLES = min(BENCH_SAMPLES, 5)

function measure(case, f, input, outcome)
    @nospecialize f input outcome
    firstcall, compile, samples = benchtime(f, input; samples = SAMPLES)
    times = [s.time for s in samples]
    println(join((label, case, firstcall.time, compile, minimum(times),
        median(times), minimum(s.bytes for s in samples),
        outcome(firstcall.value)), '\t'))
    flush(stdout)
    return firstcall.value
end

# A resonantly phase matched JTWPA of `Nj` cells between two 50 ohm ports:
# each cell a junction with its capacitance and its capacitance to ground,
# and a resonator coupled to it which phase matches the pump. The first
# cell and the last node carry half the ground capacitance.
function rpmjtwpa(Nj; Lj = IctoLj(3.4e-6), Cg = 45.0e-15, Cc = 15.0e-15,
        Cr = 2.8153e-12, Lr = 1.70e-10, Cj = 55e-15)
    cell(Cgnd) = Circuit(
        [(:jj, 1, 2, JosephsonJunction(Lj)), (:cj, 1, 2, Capacitor(Cj)),
         (:cg, 1, 0, Capacitor(Cgnd)), (:cc, 1, 3, Capacitor(Cc)),
         (:cr, 3, 0, Capacitor(Cr)), (:lr, 3, 0, Inductor(Lr))];
        pins = [1 => (:jj, 1), 2 => (:jj, 2)])
    first, inner = cell((Cg - Cc)/2), cell(Cg - Cc)
    netlist = Any[(:p1, 1, 0, Port(1; Z0 = 50.0))]
    for i in 1:Nj
        push!(netlist, (Symbol(:cell, i), i, i + 1, i == 1 ? first : inner))
    end
    push!(netlist, (:cend, Nj + 1, 0, Capacitor((Cg - Cc)/2)))
    push!(netlist, (:p2, Nj + 1, 0, Port(2; Z0 = 50.0)))
    return Circuit(netlist)
end

const WP = 2*pi*7.12e9
const WS = 2*pi*6.4e9
converged(sol) = "converged=$(sol.solverinfo.converged)"

println("label\tcase\tfirst_s\tcompile_s\tmin_s\tmedian_s\tbytes\toutcome")

# The pumped line and its gain sweep, the amplifier's everyday solve.
const N512 = benchsize(512)
const WSWEEP = 2*pi*collect(range(4e9, 9e9; length = benchsize(101)))
line = rpmjtwpa(N512)
pumped(c) = hbsolve(WSWEEP, (WP,), [(mode = (1,), port = 1, current = 1.85e-6)],
    (8,), (16,), c; keyedarrays = false)
function gain(sol)
    lin = sol.linearized
    k = argmin(abs.(lin.w .- 2*pi*6e9))
    nm = length(lin.modes)
    s21 = lin.S[lin.signalindex + nm, lin.signalindex, k]
    return "converged=$(sol.nonlinear.solverinfo.converged) gain=$(round(10*log10(abs2(s21)); digits = 2)) dB at $(round(lin.w[k]/(2*pi*1e9); digits = 2)) GHz"
end
measure("rpm$(N512)/hbsolve", pumped, line, gain)

# Two pump tones of equal strength on a shorter line: the pump solve on a
# two dimensional grid of mixing products, with the direct current and
# every order retained.
const N128 = benchsize(128)
twotone(c) = hbnlsolve((WP, WS), (8, 4),
    [(mode = (1, 0), port = 1, current = 1.1e-6),
     (mode = (0, 1), port = 1, current = 1.1e-6)], c;
    dc = true, odd = true, even = true, keyedarrays = false)
measure("rpm$(N128)/hbnlsolve, two pumps", twotone, rpmjtwpa(N128), converged)

# A strong pump on a long line, three quarters of the critical current,
# which the Newton iteration takes many damped steps to reach.
const N1024 = benchsize(1024)
strong(c) = hbnlsolve((WP,), (16,), [(mode = (1,), port = 1, current = 2.5e-6)],
    c; keyedarrays = false)
measure("rpm$(N1024)/hbnlsolve, strong pump", strong, rpmjtwpa(N1024), converged)

# The pumped line in time: a batch of four pump powers, the pump rising
# over 4 ns under a weak signal, solved together.
rise(t, T) = t <= 0 ? 0.0 : t >= T ? 1.0 : sinpi(t/(2*T))^2
function lineproblem(x, ip)
    drive(t) = 2*ip*rise(t, 4e-9)*cospi(2*7.12e9*t) +
        2*2e-9*rise(t, 1e-9)*cospi(2*6.0e9*t)
    return transientproblem(x; sources = [TransientSource(1, drive)])
end
const NSTEPS = benchsize(8000)
const DT = 2.5e-12
problem = lineproblem(line, 1.85e-6)
batch = [lineproblem(problem, 1.85e-6*s) for s in (0.8, 0.9, 1.0, 1.1)]
const REUSE = TransientReuse()
stepped(p) = transientsolve(p, (0.0, NSTEPS*DT); dt = DT, reuse = REUSE)
measure("rpm$(N512)/transientsolve, batch of 4", stepped, batch,
    b -> "finite=$(all(isfinite, b.outgoing))")
