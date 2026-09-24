# Latency benchmark: the time to load the package and to its first solves,
# and how many methods and method instances it carries. The process is the
# measurement, so run it fresh, with one Julia/BLAS thread:
# julia --startup-file=no --project=. benchmark/latency.jl label [--precompile]
# `--precompile` first builds the package's precompiled image anew and
# times that; it takes minutes and writes a new image into the depot.
#
# A method instance is a method specialized for the types of one call,
# inferred and compiled for it; those the precompile workload made are in
# the image, and those a first solve adds were inferred and compiled on
# that call, so their count is how much the solve left to compile.
include(joinpath(@__DIR__, "options.jl"))
const label = BENCH_LABEL
const PRECOMPILE = "--precompile" in ARGS && !BENCH_SMOKE

row(measure, value) = (println(join((label, measure, value), '\t')); flush(stdout))
println("label\tmeasure\tvalue")

if PRECOMPILE
    let pkg = Base.identify_package("JosephsonCircuits")
        row("precompile_s", @elapsed Base.compilecache(pkg))
    end
end
row("load_s", @elapsed @eval using JosephsonCircuits)

# the methods the package defines for its own functions and types, and the
# instances specialized from them so far
function counts(mod::Module)
    seen = Set{Method}()
    ninstances = 0
    for name in names(mod; all = true)
        isdefined(mod, name) || continue
        v = getfield(mod, name)
        (v isa Function || v isa Type) || continue
        for m in methods(v)
            parentmodule(m) === mod || continue
            m in seen && continue
            push!(seen, m)
            ninstances += count(!isnothing, Base.specializations(m))
        end
    end
    return length(seen), ninstances
end
methodsloaded, instancesloaded = counts(JosephsonCircuits)
row("methods", methodsloaded)
row("instances_loaded", instancesloaded)

# a first call: its time, its compile time, which `@timed` reports from
# Julia 1.11 on, and the instances it added
function firstcall(case, f)
    before = last(counts(JosephsonCircuits))
    t = @timed Base.invokelatest(f)
    row("$(case)_s", t.time)
    row("$(case)_compile_s", hasproperty(t, :compile_time) ? t.compile_time : NaN)
    row("$(case)_instances", last(counts(JosephsonCircuits)) - before)
end

# the JPA of the precompile workload, whose solves the image holds
jpa() = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1000e-12)),
    (:cj, 2, 0, Capacitor(1000e-15))])
const IP = 0.00565e-6
const WP = 2*pi*4.75001e9

# its gain sweep in frequency and its pumped response in time
firstsweep() = hbsolve(2*pi*(4.5:0.05:5.0)*1e9, (WP,),
    [(mode = (1,), port = 1, current = IP)], (8,), (16,), jpa())
firsttransient() = transientsolve(
    transientproblem(jpa(); sources = [TransientSource(1,
        t -> 2*IP*min(t/2e-9, 1.0)*cos(WP*t))]),
    (0.0, 10e-9); dt = 2.6e-12)

# a JTWPA of subcircuit cells and two ports pumped by two tones, a
# circuit and a grid the workload does not hold
function firstline()
    cell = Circuit([(:jj, 1, 2, JosephsonJunction(IctoLj(3.4e-6))),
        (:cj, 1, 2, Capacitor(55e-15)), (:cg, 2, 0, Capacitor(45e-15))];
        pins = [1 => (:jj, 1), 2 => (:jj, 2)])
    line = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0));
        [(Symbol(:cell, i), i, i + 1, cell) for i in 1:16];
        (:p2, 17, 0, Port(2; Z0 = 50.0))])
    return hbsolve(2*pi*(5.0:0.5:7.0)*1e9, (2*pi*7.12e9, 2*pi*6.4e9),
        [(mode = (1, 0), port = 1, current = 1.0e-6),
         (mode = (0, 1), port = 1, current = 1.0e-6)], (2, 2), (4, 4), line)
end

firstcall("firstsweep", firstsweep)
firstcall("firsttransient", firsttransient)
firstcall("firstline", firstline)
row("instances", last(counts(JosephsonCircuits)))
