using JosephsonCircuits
using JosephsonCircuits: TransientBatchSolution, TransientStepError
using LinearAlgebra
using SparseArrays
using Random
using Test

# a backend which is not the host, to check that a device batch keeps the
# one chunk its uniform batch already is, without needing a device
struct NotTheHost <: JosephsonCircuits.KernelAbstractions.GPU end

# an adjoint with a sink which captures a buffer, and a weak reference to
# the buffer, for the check that a workspace kept by a reuse holds
# nothing of the call once it returns
@noinline function capturingsink(batch, weights, reuse)
    buffer = zeros(1000)
    sink = (k, c) -> (buffer[1] += sum(c); nothing)
    transientadjoint(batch, weights; sink, reuse)
    return WeakRef(buffer)
end

# the bytes a state extracted from a solution allocates, measured in a
# function of its own: `@allocated` compiles the expression it sits in,
# which at the top level of a file is the whole testset around it
stateallocations(sol) = @allocated transientstate(sol)

# the same component at a value scaled by `r`, for finite differences
rescale(c::Capacitor, r) = Capacitor(r*c.C)
rescale(c::Inductor, r) = Inductor(r*c.L)
rescale(c::Resistor, r) = Resistor(r*c.R)
rescale(c::NonlinearInductor, r) = JosephsonJunction(r*c.L0)
scaledentries(netlist, name, r) =
    [c[1] == name ? (c[1], c[2], c[3], rescale(c[4], r)) : c for c in netlist]

# the fixtures every testset reads, constant so that each testset,
# compiled on its own, reads them without a global lookup
const JC = JosephsonCircuits
const rc = [("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "0", Capacitor(1e-12))]

@testset "the RC response, the rules, the ports and the sources" begin
    prob = transientproblem(Circuit(rc); sources = [TransientSource(1, 1e-6)])
    exact = 50e-6*(1 - exp(-4))
    errors = map((5e-12, 2.5e-12, 1.25e-12)) do dt
        abs(transientsolve(prob, (0.0, 200e-12); dt, method = Trapezoidal()).voltage[1, end] - exact)
    end
    # the trapezoidal rule is second order: the error falls by four per
    # halving
    @test errors[1]/errors[2] ≈ 4 rtol=0.01
    @test errors[2]/errors[3] ≈ 4 rtol=0.01
    coarse = transientsolve(prob, (0.0, 200e-12); dt = 5e-12)
    @test coarse.stats.factorizations == 1
    @test coarse.incident ≈ fill(1e-6*sqrt(50)/2, size(coarse.incident))
    @test coarse.voltage ≈ sqrt(50)*(coarse.incident + coarse.outgoing)
    @test isnothing(coarse.flux)
    decimated = transientsolve(prob, (0.0, 200e-12); dt = 5e-12, saveevery = 7)
    @test decimated.times[[1, end]] == [0.0, 200e-12]
    @test decimated.finalflux == coarse.finalflux
    @test length(decimated.times) == 7
    # backward Euler has a closed form on the RC
    be = transientsolve(prob, (0.0, 200e-12); dt = 5e-12, method = BackwardEuler())
    @test be.voltage[1, end] ≈ 50e-6*(1 - (1 + 0.1)^(-40)) rtol=1e-10
    # a named current source draws its current from the node at its
    # first terminal and delivers it to the node at its second, the
    # opposite sense of a port source, which injects into the port's
    # positive terminal; a named drive replaces the constant value, an
    # unnamed constant source stays
    net = vcat(rc, [("I1", "1", "0", CurrentSource(2e-6))])
    named = transientsolve(transientproblem(Circuit(net); sources = [TransientSource(:I1, 1e-6)]),
        (0.0, 200e-12); dt = 5e-12)
    @test named.voltage ≈ -coarse.voltage
    static = transientsolve(transientproblem(Circuit(net)), (0.0, 200e-12); dt = 5e-12)
    @test static.voltage ≈ -2coarse.voltage
    typed = Circuit([("p", "1", "0", Port(1)), ("c", "1", "0", Capacitor(1e-12))])
    ts = transientsolve(transientproblem(typed; sources = [TransientSource(1, 1e-6)]),
        (0.0, 200e-12); dt = 5e-12)
    @test ts.voltage ≈ coarse.voltage
    # the port axis is in the order of the port numbers, harmonic
    # balance's, whatever order the ports are declared in: a two port
    # declared port 2 first, driven at port 1, reflects and transmits
    # what the linearized solver says S11 and S21 are
    swapped = Circuit([("p2", "2", "0", Port(2; Z0 = 50.0)), ("p1", "1", "0", Port(1; Z0 = 50.0)),
        ("c12", "1", "2", Capacitor(1e-12)), ("c2", "2", "0", Capacitor(2e-12)), ("r1", "1", "0", Resistor(100.0))])
    S = hblinsolve([2pi*1e9], swapped; keyedarrays = false).S
    drive = transientsolve(transientproblem(swapped; sources = [TransientSource(1, t -> 1e-6*sinpi(2e9*t))]),
        (0.0, 40e-9); dt = 1e-12)
    @test maximum(abs, drive.incident[2, :]) <= 1e-12*maximum(abs, drive.incident[1, :])
    window = t -> t < 20e-9 ? 0.0 : sinpi((t - 20e-9)/20e-9)^2
    wave = (q, k) -> transientdemodulate(drive, q, 2pi*1e9; quantity = k, window)
    @test wave(1, :outgoing)/wave(1, :incident) ≈ S[1, 1, 1] rtol=1e-6
    @test wave(2, :outgoing)/wave(1, :incident) ≈ S[2, 1, 1] rtol=1e-6
end

@testset "a weak drive is converged relative to itself" begin
    # a linear circuit's response scales with its drive, so with the
    # absolute tolerance below the drive a picoampere is solved to
    # `rtol` of itself as a microampere is, under every rule
    lc = Circuit([(:p, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(50e-15)),
        (:l, 2, 0, Inductor(1e-9)), (:c, 2, 0, Capacitor(1e-12))])
    tone(a) = t -> a*sinpi(2e9*t)*(t <= 0 ? 0.0 : t >= 1e-9 ? 1.0 : sinpi(t/2e-9)^2)
    strong = transientproblem(lc; sources = [TransientSource(1, tone(1e-6))])
    weak = transientproblem(strong; sources = [TransientSource(1, tone(1e-12))])
    for method in (GaussLegendre(), Trapezoidal(), BackwardEuler())
        s = transientsolve(strong, (0.0, 2e-9); dt = 2e-12, method, atol = 1e-20)
        w = transientsolve(weak, (0.0, 2e-9); dt = 2e-12, method, atol = 1e-20)
        @test maximum(abs, 1e6 .* w.outgoing .- s.outgoing) < 1e-8*maximum(abs, s.outgoing)
    end
    # each row is converged relative to its own terms: beside a node
    # held at a milliampere, which no element joins to the tone's, the
    # picoampere is solved as it is alone, under every rule, and a
    # resistive node whose inductor current the state misstates by a
    # tenth of a picoampere is refused at the start
    held = vcat([(:p, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(50e-15)), (:l, 2, 0, Inductor(1e-9)),
        (:c, 2, 0, Capacitor(1e-12))], [(:lb, 3, 0, Inductor(1e-9)), (:cb, 3, 0, Capacitor(1e-12)),
        (:ib, 0, 3, CurrentSource(1e-3))])
    biased = transientproblem(Circuit(held); sources = [TransientSource(1, tone(1e-12))])
    rest = transientstate(biased; flux = [0.0, 0.0, 1e-12])
    for method in (GaussLegendre(), Trapezoidal(), BackwardEuler())
        w = transientsolve(weak, (0.0, 2e-9); dt = 2e-12, method, atol = 1e-20)
        b = transientsolve(biased, (0.0, 2e-9); dt = 2e-12, method, atol = 1e-20, initialstate = rest)
        @test maximum(abs, b.outgoing .- w.outgoing) < 1e-8*maximum(abs, w.outgoing)
    end
    island = transientproblem(Circuit(vcat(held, [(:ld, 4, 0, Inductor(1e-9)), (:rd, 4, 0, Resistor(50.0))])))
    @test_throws ArgumentError transientsolve(island, (0.0, 1e-11); dt = 1e-12, atol = 1e-20,
        initialstate = transientstate(island; flux = [0.0, 0.0, 1e-12, 1e-22]))
    # and the start is checked against the same tolerances: a resistor
    # driven from the first sample by a current far below the scaled
    # unit is refused from rest
    faint = transientproblem(Circuit(rc[1:1]); sources = [TransientSource(1, 1e-18)])
    @test_throws ArgumentError transientsolve(faint, (0.0, 1e-11); dt = 1e-12, atol = 1e-20)
end

@testset "a running junction's phase does not loosen its step" begin
    # a junction biased past its critical current winds its phase
    # without bound, and nothing in the circuit sees a whole turn of
    # it: started a thousand turns up it is the same circuit, whose
    # steps converge as closely as from zero under every rule
    running = transientproblem(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)),
        ("Lj1", "1", "0", JosephsonJunction(1e-9)), ("C1", "1", "0", Capacitor(0.2e-12))]);
        sources = [TransientSource(1, t -> 0.6e-6*(t <= 0 ? 0.0 : t >= 0.5e-9 ? 1.0 : (1 - cospi(t/0.5e-9))/2))])
    turned = transientstate(running; flux = [2000pi*JC.phi0])
    for method in (Trapezoidal(), BackwardEuler(), GaussLegendre())
        s0 = transientsolve(running, (0.0, 2e-9); dt = 1e-12, method)
        s1 = transientsolve(running, (0.0, 2e-9); dt = 1e-12, method, initialstate = turned)
        @test maximum(abs, s1.voltage .- s0.voltage) < 1e-9*maximum(abs, s0.voltage)
    end
end

@testset "a lossless LC keeps its energy, and coupled inductors their modes" begin
    L, C, V = 1e-9, 1e-12, 1e-6
    circuit = Circuit([("p", "1", "0", Port(1; termination = nothing)),
        ("c", "1", "0", Capacitor(C)), ("l", "1", "0", Inductor(L))])
    prob = transientproblem(circuit)
    period = 2pi*sqrt(L*C)
    # the trapezoidal rule, free of numerical damping, keeps it, and
    # warps its frequency by (2 pi dt/period)^2/12
    sol = transientsolve(prob, (0.0, 3period); dt = period/500,
        initialstate = transientstate(prob; voltage = [V]), record = :states, method = Trapezoidal())
    v = sol.rate[1, :] .* JC.phi0
    flux = sol.flux[1, :] .* JC.phi0
    energy = C .* v .^ 2 ./ 2 .+ flux .^ 2 ./ (2L)
    @test maximum(abs.(energy ./ energy[1] .- 1)) < 2e-10
    @test maximum(abs.(v .- V .* cos.(sol.times ./ sqrt(L*C)))) < 3e-4V
    # mutually coupled inductors go through the augmentation of the
    # modified nodal analysis: two auxiliary currents, and the normal
    # modes of the coupled tanks under the trapezoidal rule
    coupled = Circuit([("C1", "1", "0", Capacitor(C)), ("C2", "2", "0", Capacitor(C)),
        ("L1", "1", "0", Inductor(L)), ("L2", "2", "0", Inductor(L)), ("K1", "L1", "L2", MutualInductor(0.3))])
    p2 = transientproblem(coupled)
    @test p2.Naux == 2 && length(p2) == 4 && isempty(p2.gaugeindices)
    s2 = transientsolve(p2, (0.0, 2period); dt = period/800,
        initialstate = transientstate(p2; voltage = [V, 0]), record = :states, method = Trapezoidal())
    wp, wm = 1/sqrt(L*C*1.3), 1/sqrt(L*C*0.7)
    e1 = V/2 .* (cos.(wp .* s2.times) .+ cos.(wm .* s2.times))
    e2 = V/2 .* (cos.(wp .* s2.times) .- cos.(wm .* s2.times))
    @test maximum(abs.(s2.rate[1, :] .* JC.phi0 .- e1)) < 2e-4V
    @test maximum(abs.(s2.rate[2, :] .* JC.phi0 .- e2)) < 2e-4V
    # a nonzero initial flux on a coupled tank: the auxiliary currents
    # of the state are in the problem's fixed units, so the state is
    # consistent at any step, and the fluxes follow the normal modes
    F = 0.1*JC.phi0
    s3 = transientsolve(p2, (0.0, 2period); dt = period/800,
        initialstate = transientstate(p2; flux = [F, 0]), record = :states, method = Trapezoidal())
    f1 = F/2 .* (cos.(wp .* s3.times) .+ cos.(wm .* s3.times))
    f2 = F/2 .* (cos.(wp .* s3.times) .- cos.(wm .* s3.times))
    @test maximum(abs.(s3.flux[1, :] .* JC.phi0 .- f1)) < 2e-4F
    @test maximum(abs.(s3.flux[2, :] .* JC.phi0 .- f2)) < 2e-4F
    # the auxiliary rates follow the node rates by the same constitutive
    # equations, and a restart at another step continues the solve
    state = transientstate(p2; flux = [F, 0], voltage = [V, 0])
    @test state.flux[3:4] ≈ (p2.Lscale/JC.phi0) .* ([L 0.3L; 0.3L L] \ [F, 0])
    @test state.rate[3:4] ≈ (p2.Lscale/JC.phi0) .* ([L 0.3L; 0.3L L] \ [V, 0])
    half = transientsolve(p2, (0.0, period); dt = period/800, initialstate = state)
    # the end of one solve starts the next, as a state and not as a
    # pair of arrays, whose units the solver could not tell
    continued = transientstate(half)
    @test continued isa TransientState && continued.flux == half.finalflux && continued.rate == half.finalrate
    @test_throws TypeError transientsolve(p2, (period, 2period); dt = period/800, initialstate = (half.finalflux, half.finalrate))
    @test_throws ArgumentError transientsolve([p2, p2], (period, 2period); dt = period/800, initialstate = (half.finalflux, half.finalrate))
    rest = transientsolve(p2, (period, 2period); dt = period/800, initialstate = continued)
    whole = transientsolve(p2, (0.0, 2period); dt = period/800, initialstate = state)
    @test rest.finalflux ≈ whole.finalflux rtol=1e-9
    @test rest.finalrate ≈ whole.finalrate rtol=1e-9
    # the rate along the inductor constraint is read from it wherever
    # the state is reported, under either rule, so the end of a long
    # solve satisfies the constraint and restarts
    for method in (GaussLegendre(), Trapezoidal())
        long = transientsolve(p2, (0.0, 16period); dt = period/800, initialstate = state,
            record = :states, method)
        sysm = JC.transientsystem(p2, period/800, method, JC.CPU(), KLUfactorization())
        @test JC.transientconsistency(sysm, long.finalflux, long.finalrate, 16period).violation < 1e-12
        @test JC.transientconsistency(sysm, long.flux[:, end], long.rate[:, end], 16period).violation < 1e-12
        @test transientsolve(p2, (16period, 17period); dt = period/800, method,
            initialstate = transientstate(long)).stats.steps == 800
    end
    coarse = transientsolve(p2, (period, 2period); dt = period/400,
        initialstate = continued)
    # across a transmission line the state carries the waves over the
    # delay window, so a split solve is the uninterrupted one at the
    # same step, and follows it at another step through the
    # interpolation of the history, on a reactive load and past the
    # delay
    cabled = transientproblem(Circuit([(:p1, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(10e-12)),
        (:line, 1, 2, TransmissionLine(60.0, 0.09)), (:c2, 2, 0, Capacitor(20e-12)), (:p2, 2, 0, Port(2))]);
        sources = [TransientSource(1, t -> 1e-6*sinpi(2e9*t))])
    gl = (; method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
    uninterrupted = transientsolve(cabled, (0.0, 2.4e-9); dt = 2e-12, record = :phases, gl...)
    first = transientsolve(cabled, (0.0, 0.8e-9); dt = 2e-12, gl...)
    restart = transientstate(first)
    @test size(restart.waves, 2) == JC.lineprehistory(cabled, 2e-12) && restart.wavesdt == 2e-12
    second = transientsolve(cabled, (0.8e-9, 2.4e-9); dt = 2e-12, initialstate = restart, gl...)
    after = length(first.times):length(uninterrupted.times)
    @test second.voltage ≈ uninterrupted.voltage[:, after] rtol=1e-8
    @test second.finalrate ≈ uninterrupted.finalrate rtol=1e-8
    finer = transientsolve(cabled, (0.8e-9, 2.4e-9); dt = 1e-12, initialstate = restart, gl...)
    @test finer.voltage[:, 1:2:end] ≈ uninterrupted.voltage[:, after] rtol=1e-4
    # a constant prehistory is not a continuation: it is wrong past
    # the delay too
    constant = TransientState(restart.flux, restart.rate, restart.waves[:, end:end], 0.0, restart.blockstates)
    wrong = transientsolve(cabled, (0.8e-9, 2.4e-9); dt = 2e-12, initialstate = constant, gl...)
    late = 4*150:length(after)
    @test !isapprox(wrong.voltage[:, late], uninterrupted.voltage[:, after][:, late]; rtol = 1e-2)
    # a solve keeps the history it started from, so a segment shorter
    # than the delay still hands the next one a complete window; the
    # segment's span is a whole number of steps, so the grid is the same
    rng = Random.default_rng()
    tmid = first.times[end] + 50*(first.times[2] - first.times[1])
    short = transientsolve(cabled, (first.times[end], tmid); dt = 2e-12, initialstate = restart, gl...)
    @test length(short.times) == 51
    @test size(short.history) == size(restart.waves)
    @test short.history ≈ restart.waves rtol=1e-12
    third = transientsolve(cabled, (tmid, 2.4e-9); dt = 2e-12, initialstate = transientstate(short), gl...)
    @test third.voltage ≈ uninterrupted.voltage[:, length(first.times) + length(short.times) - 1:end] rtol=1e-8
    @test third.finalrate ≈ uninterrupted.finalrate rtol=1e-8
    # a checkpointed solve started from that history replays its first
    # window from it, so its sensitivities are those of the full record
    rstates = transientsolve(cabled, (0.8e-9, 1.8e-9); dt = 2e-12, initialstate = restart, record = :states, gl...)
    rcps = transientsolve(cabled, (0.8e-9, 1.8e-9); dt = 2e-12, initialstate = restart, record = :checkpoints, checkpointevery = 50, gl...)
    @test transientsensitivity(rcps, ["c1"]).outgoing ≈ transientsensitivity(rstates, ["c1"]).outgoing rtol=1e-7
    rweights = randn(rng, 2, length(rstates.times))
    @test transientadjoint(rcps, rweights; components = ["c1"]).sensitivity ≈
        transientadjoint(rstates, rweights; components = ["c1"]).sensitivity rtol=1e-7
    # a line between algebraic terminals, no capacitor at either end:
    # the check of the initial state reads the arriving waves from the
    # history as the steps do, so the restart is accepted and exact
    bare = transientproblem(Circuit([(:p1, 1, 0, Port(1)), (:line, 1, 2, TransmissionLine(60.0, 0.09)), (:p2, 2, 0, Port(2))]);
        sources = [TransientSource(1, t -> 1e-6*sinpi(2e9*t))])
    bwhole = transientsolve(bare, (0.0, 1.6e-9); dt = 2e-12, gl...)
    bfirst = transientsolve(bare, (0.0, 0.8e-9); dt = 2e-12, gl...)
    bsecond = transientsolve(bare, (0.8e-9, 1.6e-9); dt = 2e-12, initialstate = transientstate(bfirst), gl...)
    @test bsecond.voltage ≈ bwhole.voltage[:, length(bfirst.times):end] rtol=1e-8
    # the check reads those waves at the state's own step, so a restart
    # at a coarser one is accepted: the interpolation which resamples
    # the history is an error of the continuation and not an
    # inconsistency of the state
    fast = transientproblem(Circuit([(:p1, 1, 0, Port(1)), (:line, 1, 2, TransmissionLine(60.0, 0.09)), (:p2, 2, 0, Port(2))]);
        sources = [TransientSource(1, t -> 1e-6*sinpi(20e9*t))])
    ffirst = transientsolve(fast, (0.0, 0.8e-9); dt = 2e-12, gl...)
    fstate = transientstate(ffirst)
    fwhole = transientsolve(fast, (0.0, 1.6e-9); dt = 4e-12, gl...)
    fcoarse = transientsolve(fast, (0.8e-9, 1.6e-9); dt = 4e-12, initialstate = fstate, gl...)
    @test fcoarse.voltage[:, end] ≈ fwhole.voltage[:, end] rtol=1e-4
    @test transientsolve(fast, (0.8e-9, 1.6e-9); dt = 1e-11, initialstate = fstate, gl...) isa JC.TransientSolution
    # extracting a state reads the delay window alone, so it costs what
    # the window does for a record of any length
    brief = transientsolve(cabled, (0.0, 0.4e-9); dt = 2e-12, record = :phases, gl...)
    stateallocations(brief), stateallocations(uninterrupted)
    @test size(uninterrupted.linewaves, 2) > 5*size(brief.linewaves, 2)
    @test stateallocations(uninterrupted) < 2*stateallocations(brief)
    # the end of a solve continues whatever it recorded: a rational
    # block's states under a record of the ports or of checkpoints,
    # and the waves of a line short enough that checkpoints keep its
    # history at each of them rather than the record of its waves
    a = 2*50.0/2e-9
    u = [1.0, -1.0]
    ind = RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
    blocked = transientproblem(Circuit([(:p1, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(0.3e-12)), (:b, 1, 2, ind),
        (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2))]); sources = [TransientSource(1, t -> 1e-6*sinpi(2e9*t))])
    shortline = transientproblem(Circuit([(:p1, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(10e-12)),
        (:line, 1, 2, TransmissionLine(60.0, 0.009)), (:c2, 2, 0, Capacitor(20e-12)), (:p2, 2, 0, Port(2))]);
        sources = [TransientSource(1, t -> 1e-6*sinpi(2e9*t))])
    for (prob, record) in ((blocked, :ports), (blocked, :checkpoints), (shortline, :checkpoints))
        onego = transientsolve(prob, (0.0, 1e-9); dt = 2e-12, gl...)
        part = transientsolve(prob, (0.0, 0.5e-9); dt = 2e-12, record, checkpointevery = 50, gl...)
        restof = transientsolve(prob, (0.5e-9, 1e-9); dt = 2e-12, initialstate = transientstate(part), gl...)
        @test restof.voltage ≈ onego.voltage[:, length(part.times):end] rtol=1e-8
    end
    @test coarse.finalflux ≈ whole.finalflux rtol=1e-3
end

@testset "the initial state, the gauge and the unsupported cases" begin
    # a resistor driven from the first sample needs its voltage supplied
    resistor = transientproblem(Circuit(rc[1:1]); sources = [TransientSource(1, 1e-6)])
    @test_throws ArgumentError transientsolve(resistor, (0.0, 1e-9); dt = 1e-12)
    rs = transientsolve(resistor, (0.0, 1e-9); dt = 1e-12,
        initialstate = transientstate(resistor; voltage = [50e-6]))
    @test rs.voltage ≈ fill(50e-6, size(rs.voltage))
    # a singular capacitance matrix without a zero row: two terminated
    # nodes joined by one capacitor constrain the sum of their voltages
    # to the drive, so the zero state is rejected, and the consistent
    # one relaxes without ringing to the resistive division
    pair = transientproblem(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("P2", "2", "0", Port(2; Z0 = 50.0)),
        ("C1", "1", "2", Capacitor(1e-12))]); sources = [TransientSource(1, 1e-6)])
    @test pair.inertialess == [[1, 2]] && isempty(pair.algebraic)
    @test_throws ArgumentError transientsolve(pair, (0.0, 1e-9); dt = 1e-12)
    ps = transientsolve(pair, (0.0, 1e-9); dt = 1e-12,
        initialstate = transientstate(pair; voltage = [25e-6, 25e-6]))
    decay = exp.(-ps.times ./ 100e-12)
    @test ps.voltage[1, :] ≈ 50e-6 .- 25e-6 .* decay rtol=2e-4
    @test ps.voltage[2, :] ≈ 25e-6 .* decay rtol=2e-4 atol=1e-9
    @test maximum(abs.(ps.voltage[1, :] .+ ps.voltage[2, :] .- 50e-6)) < 1e-12
    # a capacitive island reached through an inductor, and a node no
    # capacitor touches between resistors: the same check
    island2 = transientproblem(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(1e-12)),
        ("L2", "2", "0", Inductor(1e-9))]); sources = [TransientSource(1, 1e-6)])
    @test island2.inertialess == [[1, 2]] && isempty(island2.algebraic)
    @test_throws ArgumentError transientsolve(island2, (0.0, 1e-9); dt = 1e-12)
    i2 = transientsolve(island2, (0.0, 1e-9); dt = 1e-12,
        initialstate = transientstate(island2; voltage = [50e-6, 50e-6]))
    @test i2.voltage[1, 1] ≈ 50e-6
    divider = transientproblem(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "0", Capacitor(1e-12)),
        ("R2", "1", "2", Resistor(50.0)), ("R3", "2", "0", Resistor(50.0))]); sources = [TransientSource(1, 1e-6)])
    @test divider.inertialess == [[2]] && isempty(divider.algebraic)
    @test_throws ArgumentError transientsolve(divider, (0.0, 1e-9); dt = 1e-12,
        initialstate = transientstate(divider; voltage = [50e-6, 0.0]))
    d2 = transientsolve(divider, (0.0, 1e-9); dt = 1e-12,
        initialstate = transientstate(divider; voltage = [100e-6/3, 50e-6/3]))
    @test d2.voltage[1, :] ≈ fill(100e-6/3, size(d2.voltage, 2)) rtol=1e-8
    # a drive starting with the record, `t <= 0 ? 0.0 : ...`, into a node
    # without capacitance or resistance: the zero state is consistent,
    # since the integration reads no drive before its start, and the
    # record is the one started earlier from rest, under either rule
    ramped = Circuit([(:p, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(1e-12)), (:l12, 1, 2, Inductor(1e-9)),
        (:l2, 2, 0, Inductor(2e-9)), (:jj, 2, 0, JosephsonJunction(1e-9)), (:ib, 0, 2, CurrentSource(0.0))])
    rise(t) = t <= 0 ? 0.0 : t >= 0.5e-9 ? 1.0 : sinpi(t/1e-9)^2
    rp = transientproblem(ramped; sources = [TransientSource(:ib, t -> 0.1e-6*rise(t))])
    @test rp.algebraic == [[2]]
    for method in (GaussLegendre(), Trapezoidal())
        early = transientsolve(rp, (-0.1e-9, 1e-9); dt = 2e-12, method)
        @test transientsolve(rp, (0.0, 1e-9); dt = 2e-12, method).voltage ≈ early.voltage[:, 51:end] rtol=1e-8
    end
    # a resistor inside a capacitive island carries no current along
    # the island's direction, so the direction is algebraic and its
    # rate is read by the differentiated equation: equal voltages on
    # the two nodes would ring between the inductors, opposite ones
    # are the loop's own mode
    loop = transientproblem(Circuit([("C1", "1", "2", Capacitor(1e-12)), ("R1", "1", "2", Resistor(50.0)),
        ("L1", "1", "0", Inductor(1e-9)), ("L2", "2", "0", Inductor(1e-9))]))
    @test loop.inertialess == [[1, 2]] && loop.algebraic == [[1, 2]]
    @test_throws ArgumentError transientsolve(loop, (0.0, 10e-12); dt = 1e-12,
        initialstate = transientstate(loop; voltage = [1e-6, 1e-6]))
    ls = transientsolve(loop, (0.0, 100e-12); dt = 1e-12, record = :states,
        initialstate = transientstate(loop; voltage = [1e-6, -1e-6]))
    @test maximum(abs.(ls.rate[1, :] .+ ls.rate[2, :])) < 1e-6*maximum(abs.(ls.rate[1, :]))
    @test abs(ls.rate[1, 2] - ls.rate[1, 1]) < 0.05*abs(ls.rate[1, 1])
    # two islands joined by a resistor are one algebraic direction, and
    # none once a resistor reaches ground
    twin = [("C1", "1", "2", Capacitor(1e-12)), ("C2", "3", "4", Capacitor(1e-12)), ("R1", "2", "3", Resistor(50.0)),
        ("L1", "1", "0", Inductor(1e-9)), ("L2", "2", "0", Inductor(1e-9)), ("L3", "3", "0", Inductor(1e-9)), ("L4", "4", "0", Inductor(1e-9))]
    twins = transientproblem(Circuit(twin))
    @test twins.inertialess == [[1, 2], [3, 4]] && twins.algebraic == [[1, 2, 3, 4]]
    @test_throws ArgumentError transientsolve(twins, (0.0, 10e-12); dt = 1e-12,
        initialstate = transientstate(twins; voltage = [1e-6, 1e-6, 1e-6, 1e-6]))
    grounded = transientproblem(Circuit(vcat(twin, [("R2", "4", "0", Resistor(50.0))])))
    @test grounded.inertialess == [[1, 2], [3, 4]] && isempty(grounded.algebraic)
    @test transientsolve(grounded, (0.0, 10e-12); dt = 1e-12,
        initialstate = transientstate(grounded; voltage = [1e-6, 1e-6, 1e-6, 0.0])) isa JC.TransientSolution
    # a node no element connects to ground has a free flux offset, and
    # only such a node gets a gauge row
    island = transientproblem(Circuit([("C1", "1", "0", Capacitor(1e-12)), ("L2", "2", "3", Inductor(1e-9)), ("C3", "2", "3", Capacitor(1e-12))]))
    @test island.gaugeindices == [2]
    @test isempty(transientproblem(Circuit(rc)).gaugeindices)
    # and a circuit with such a node, which harmonic balance refuses,
    # solves
    islanded = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "0", Capacitor(1e-12)),
        ("C2", "2", "3", Capacitor(1e-12)), ("L2", "2", "3", Inductor(1e-9))])
    @test all(isfinite, transientsolve(transientproblem(islanded), (0.0, 1e-10); dt = 1e-12).outgoing)
    # but a net current of the sources into the subnetwork leaves it no
    # path back, and is refused as harmonic balance refuses the
    # subnetwork: a constant source's when the problem is built, a
    # drive's when the drives are evaluated; a source within it drives no
    # net current and stays
    fed = [("P1", "3", "0", Port(1; Z0 = 50.0)), ("C3", "3", "0", Capacitor(1e-12)), ("I1", "3", "1", CurrentSource(0.0)),
        ("C12", "1", "2", Capacitor(1e-12)), ("L12", "1", "2", Inductor(1e-9))]
    quiet = transientproblem(Circuit(fed))
    @test_throws ArgumentError transientsolve(transientproblem(quiet;
        sources = [TransientSource(:I1, t -> 1e-6*sinpi(1e9*t))]), (0.0, 10e-12); dt = 1e-12)
    @test_throws ArgumentError transientproblem(Circuit(vcat(fed[1:2], [("I1", "3", "1", CurrentSource(1e-6))], fed[4:5])))
    # Sources whose net currents into the subnetwork cancel give it the
    # forcing of one source across it, and solve as that source does: two
    # drives of one waveform, into the floating branch and out of it,
    # against a differential drive (review of 2026-10-05, finding 2).
    # Unequal waveforms do not cancel, and are refused when evaluated
    pair = Circuit([("R1", "1", "2", Resistor(50.0)), ("C1", "1", "2", Capacitor(1e-12)),
        ("I1", "0", "1", CurrentSource(0.0)), ("I2", "2", "0", CurrentSource(0.0)),
        ("I3", "2", "1", CurrentSource(0.0))])
    wave(a) = t -> a*sinpi(2e9*t)
    branch(sources) = (sol = transientsolve(transientproblem(pair; sources), (0.0, 1e-9); dt = 1e-12,
        record = :states); sol.rate[1, :] .- sol.rate[2, :])
    @test isapprox(branch([TransientSource(:I1, wave(1e-6)), TransientSource(:I2, wave(1e-6))]),
        branch([TransientSource(:I3, wave(1e-6))]); rtol = 1e-12)
    @test_throws ArgumentError branch([TransientSource(:I1, wave(1e-6)), TransientSource(:I2, wave(1.1e-6))])
    # The balance holds each drive's net into the subnetworks it feeds,
    # two at most, and judges a subnetwork by those alone (review of
    # 2026-10-07, finding 2): along a chain of floating branches joined by
    # sources from ground back to ground it grows with the chain, not with
    # its square, and passes balanced drives and refuses one more current
    # out of the last branch
    joined(m) = transientproblem(Circuit(vcat([("C$(r)", "a$(r)", "b$(r)", Capacitor(1e-12)) for r in 1:m],
            [("I$(r)", r == 0 ? "0" : "b$(r)", r == m ? "0" : "a$(r + 1)", CurrentSource(0.0)) for r in 0:m]));
        sources = [TransientSource(Symbol("I$(r)"), t -> 1e-6) for r in 0:m])
    long = joined(128)
    @test Base.summarysize(long.balance) < 3*Base.summarysize(joined(64).balance)
    @test isnothing(JC.checkbalance(long, ones(129), 0.0))
    @test_throws ArgumentError JC.checkbalance(long, [ones(128); 2.0], 0.0)
    # A tangent's currents into the floating branch are judged as the
    # sources are, and an adjoint's targets by whether their currents can
    # cancel there (follow-up review of 2026-10-05): along the pair
    # together the tangent and the adjoint are those of one source across
    # the branch, under either rule's responses; along I1 alone, whose
    # solve is refused, the tangent is refused, and I1 is refused as an
    # adjoint's target without I2, its derivative alone set by which node
    # of the branch the gauge row sits on
    across = Circuit([("C1", "1", "2", Capacitor(1e-12)), ("P1", "1", "2", Port(1; Z0 = 50.0)),
        ("I1", "0", "1", CurrentSource(0.0)), ("I2", "2", "0", CurrentSource(0.0)),
        ("I3", "2", "1", CurrentSource(0.0))])
    driven = transientproblem(across; sources = [TransientSource(:I1, wave(1e-6)), TransientSource(:I2, wave(1e-6))])
    for method in (GaussLegendre(), Trapezoidal())
        sol = transientsolve(driven, (0.0, 0.2e-9); dt = 1e-12, method, record = :phases)
        d = 1e-7 .* cospi.(3e9 .* sol.times)
        @test isapprox(transienttangent(sol, permutedims([d d]); targets = ["I1", "I2"]).voltage,
            transienttangent(sol, permutedims(d); targets = ["I3"]).voltage; rtol = 1e-12)
        @test_throws ArgumentError transienttangent(sol, permutedims([d zero(d)]); targets = ["I1", "I2"])
        # the staged form in two directions is judged where the rule reads
        # it (review of bundle 8, finding 1): the trapezoidal rule reads
        # the grid alone, and the stages at the last time begin no step,
        # so a net current there perturbs nothing and passes; one at a
        # stage Gauss-Legendre reads, or on the grid in the second
        # direction, is refused
        staged = [j*d[k] for _ in 1:2, s in 1:3, k in eachindex(d), j in 1:2]
        unread = copy(staged)
        method isa GaussLegendre ? (unread[1, 2:3, end, :] .= 0) : (unread[1, 2:3, :, :] .= 0)
        @test isapprox(transienttangent(sol, unread; targets = ["I1", "I2"]).voltage,
            transienttangent(sol, staged[1:1, :, :, :]; targets = ["I3"]).voltage; rtol = 1e-12)
        for at in (method isa GaussLegendre ? ((2, 3, 2), (1, 3, 2)) : ((1, 3, 2),))
            unbalanced = copy(staged)
            unbalanced[1, at...] = 0
            @test_throws ArgumentError transienttangent(sol, unbalanced; targets = ["I1", "I2"])
        end
        w = ones(1, length(sol.times))
        pairsum = sum(transientadjoint(sol, w; quantity = :voltage, targets = ["I1", "I2"]).currents; dims = 1)
        @test isapprox(pairsum, transientadjoint(sol, w; quantity = :voltage, targets = ["I3"]).currents; rtol = 1e-12)
        @test_throws ArgumentError transientadjoint(sol, w; quantity = :voltage, targets = ["I1", 1])
    end
    # The adjoint's targets are the edges of a graph on the floating
    # subnetworks and the grounded rest of the circuit, and a target is
    # refused exactly when its edge is a bridge (review of bundle 9): two
    # floating branches on a path of sources from ground to ground, one
    # source between them. The three return one another's currents, and
    # without the last the first's current has no way back; the source
    # between the branches is refused alone, and taken twice the two are a
    # loop of their own
    path = transientproblem(Circuit([("C1", "1", "2", Capacitor(1e-12)), ("P1", "1", "2", Port(1; Z0 = 50.0)),
        ("C2", "3", "4", Capacitor(1e-12)), ("P2", "3", "4", Port(2; Z0 = 50.0)),
        ("I1", "0", "1", CurrentSource(0.0)), ("I2", "2", "3", CurrentSource(0.0)), ("I3", "4", "0", CurrentSource(0.0))]))
    judge(p, targets) = JC.checkadjointtargets(p, first(JC.targetinjection(p, targets)), targets)
    @test isnothing(judge(path, ["I1", "I2", "I3"]))
    @test_throws ArgumentError judge(path, ["I1", "I2"])
    @test_throws ArgumentError judge(path, ["I2"])
    @test isnothing(judge(path, ["I2", "I2"]))
    # its memory, and the tangent's check's, is a few words per node and
    # target, whatever the number of subnetworks and of the targets that
    # drive no net current: a ring of a hundred floating branches joined by
    # sources through ground, and a thousand ports (review of bundle 8,
    # finding 3)
    ring = transientproblem(Circuit(vcat([("C$(r)", "a$(r)", "b$(r)", Capacitor(1e-12)) for r in 1:100],
        [("P1", "a1", "b1", Port(1; Z0 = 50.0))],
        [("I$(r)", r == 0 ? "0" : "b$(r)", r == 100 ? "0" : "a$(r + 1)", CurrentSource(0.0)) for r in 0:100])))
    many = vcat(["I$(r)" for r in 0:100], fill(1, 1000))
    injection, _ = JC.targetinjection(ring, many)
    judged(p, injection, targets) = @allocated JC.checkadjointtargets(p, injection, targets)
    judged(ring, injection, many)
    @test judged(ring, injection, many) <= 64*(size(injection, 1) + length(many))
    balanced(p, injection, currents) = @allocated JC.checktangentbalance(p, injection, currents, GaussLegendre())
    ringcurrents = JC.tangentcurrents(ones(length(many), 2), length(many), 2, 1)
    balanced(ring, injection, ringcurrents)
    @test balanced(ring, injection, ringcurrents) <= 64*(size(injection, 1) + length(many))
    # the search runs without recursion, along a path and around a cycle
    # longer than a call stack holds
    @test all(>(0), JC.multigraphbridges(permutedims([1:99_999 2:100_000]), 100_000))
    @test all(iszero, JC.multigraphbridges(permutedims([1:100_000 [2:100_000; 1]]), 100_000))
    within = transientproblem(Circuit(vcat(fed, [("I2", "1", "2", CurrentSource(0.0))])))
    @test transientproblem(within; sources = [TransientSource(:I2, t -> 1e-6*sinpi(1e9*t))]) isa JC.TransientProblem
    @test_throws ArgumentError transientproblem(Circuit([("C1", "1", "0", Capacitor(1e-12 + 1e-15im))]))
    @test_throws ArgumentError transientproblem(Circuit([("R1", "1", "0", Resistor(FrequencyDependent(w -> 50.0)))]))
    @test_throws ArgumentError transientproblem(Circuit([("L1", "1", "0", Inductor(0.0))]))
    # a negative capacitance between two nodes whose totals are
    # positive makes the capacitance matrix indefinite, and is refused;
    # a capacitor of no capacitance is none
    @test_throws ArgumentError transientproblem(Circuit([("P1", "1", "0", Port(1)), ("C1", "1", "0", Capacitor(1e-12)),
        ("C2", "2", "0", Capacitor(1e-12)), ("C12", "1", "2", Capacitor(-0.8e-12)), ("L2", "2", "0", Inductor(1e-9))]))
    @test transientproblem(Circuit(vcat(rc, [("C0", "1", "0", Capacitor(0.0))]))) isa JC.TransientProblem
    @test_throws ArgumentError transientproblem(Circuit(rc); sources = [TransientSource(9, 0.0)])
    @test_throws ArgumentError transientproblem(Circuit(rc); sources = [TransientSource("P1/termination", 0.0)])
    @test_throws ArgumentError transientsolve(transientproblem(Circuit(rc)), (0.0, 1e-9); dt = 0.0)
    @test_throws ArgumentError transientsolve(transientproblem(Circuit(rc)), (0.0, 1e-9); dt = 1e-12, record = :phases, saveevery = 2)
    @test_throws ArgumentError transientsolve(transientproblem(Circuit(rc)), (0.0, 1e-9); dt = 1e-12, maxsteps = 5)
end

@testset "a rejected correction is retried from its base point" begin
    # the Newton engine of the step on a scalar junction residual
    # `x + 10 sin(x) - 1` with a kept Jacobian of 1 from a phase of
    # pi/2: the first correction, to 1, is rejected, and the retry
    # must build the fresh Jacobian 11 and the correction at the base
    # point 0, not at the rejected trial, so that one retry converges
    jac = Ref(1.0)
    base = Ref(0.0)
    res(norms, r, y) = (r[1] = y[1] + 10sin(y[1]) - 1; norms[1] = abs(r[1]); nothing)
    refresh = mask -> (jac[] = 1 + 10cos(base[]); nothing)
    x = [0.0]
    baseres = (norms, r, y) -> (base[] = y[1]; res(norms, r, y))
    solve = (c, r) -> (c[1] = r[1]/jac[]; (false, 0))
    work = JC.NewtonWork(JC.CPU(), 1)
    converged, corrections = JC.newtonsolve!(x, [0.0], [0.0], [0.0], [0.0], baseres, res, refresh, solve,
        [1e-12], 15, [false], true, false, work)
    @test converged && work.fresh[1] && work.retries[1] == 1 && work.factorizations[1] == 1
    @test x[1] + 10sin(x[1]) ≈ 1 atol=1e-12
    @test corrections <= 6
    # the stepping rule reaches the same discrete solution whatever
    # the path of its factorizations
    circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "0", Capacitor(1e-15)),
        ("Lj1", "1", "0", JosephsonJunction(1e-10)), ("L1", "1", "0", Inductor(1e-9))])
    prob = transientproblem(circuit; sources = [TransientSource(1, t -> 2e-6*sinpi(2*5e9*t))])
    sol = transientsolve(prob, (0.0, 1e-9); dt = 2e-12)
    tight = transientsolve(prob, (0.0, 1e-9); dt = 2e-12, rtol = 1e-12, atol = 1e-13, iterations = 40)
    @test sol.finalflux ≈ tight.finalflux rtol=1e-6
    @test sol.voltage ≈ tight.voltage rtol=1e-6 atol=1e-12
end

@testset "the step Jacobian, the tangent and the adjoint" begin
    rng = Random.default_rng()
    circuit = Circuit(vcat(rc, [("Lj1", "1", "0", JosephsonJunction(1e-9))]))
    drive(t) = 0.12e-6*sinpi(2*3e9*t) + 0.01e-6*sinpi(2*1.3e9*t)
    perturb(t) = 0.02e-6*sinpi(2*1.7e9*t)
    prob = transientproblem(circuit; sources = [TransientSource(1, drive)])
    # the assembled step Jacobian against a finite difference of the
    # step residual, through the package's plan at one mode
    sys = JC.transientsystem(prob, 2e-12, Trapezoidal(), JC.CPU(), KLUfactorization())
    n = length(prob)
    x, d = randn(rng, n), randn(rng, n)
    phi = zeros(length(sys.lmolj)); jwork = similar(phi)
    mul!(phi, sys.RJ, x)
    JC.stepjacobian!(sys, sys.factorization, phi, nothing)
    J = copy(sys.jacobian)
    r = (y -> (res = zeros(n); JC.stepresidual!(res, sys, zeros(n), y, zeros(n), phi, zeros(n), jwork, zeros(n), ones(n)); res))
    h = 1e-6
    @test J*d ≈ (r(x + h*d) - r(x - h*d))/(2h) rtol=1e-7
    @test J ≈ transpose(J)
    for method in (Trapezoidal(), BackwardEuler())
        sol = transientsolve(prob, (0.0, 1e-9); dt = 2e-12, method, record = :phases, rtol = 1e-12)
        currents = reshape(perturb.(sol.times), 1, :)
        tangent = transienttangent(sol, currents)
        weights = randn(rng, 1, length(sol.times))
        for quantity in (:voltage, :incident, :outgoing)
            adj = transientadjoint(sol, weights; quantity)
            @test sum(weights .* getproperty(tangent, quantity)) ≈ sum(adj.currents .* currents) rtol=1e-10 atol=1e-15
        end
        eps = 1e-3
        plus = transientsolve(transientproblem(circuit; sources = [TransientSource(1, t -> drive(t) + eps*perturb(t))]),
            (0.0, 1e-9); dt = 2e-12, method, rtol = 1e-12)
        minus = transientsolve(transientproblem(circuit; sources = [TransientSource(1, t -> drive(t) - eps*perturb(t))]),
            (0.0, 1e-9); dt = 2e-12, method, rtol = 1e-12)
        @test tangent.voltage ≈ (plus.voltage - minus.voltage)/(2eps) rtol=2e-6
        # the initial state enters too; the RC junction node has a
        # capacitance so any initial state is consistent
        dx0, dv0 = [0.1], [0.05]
        initial = transienttangent(sol, zero(currents); initialstate = (dx0, dv0))
        adj = transientadjoint(sol, weights; quantity = :voltage)
        # the sums cancel by an amount the draw of the weights decides;
        # the tolerance is measured from their terms
        cancellation = sum(abs, weights .* initial.voltage) + abs(adj.initialflux[1]*dx0[1]) + abs(adj.initialrate[1]*dv0[1])
        @test sum(weights .* initial.voltage) ≈ dot(adj.initialflux, dx0) + dot(adj.initialrate, dv0) rtol=1e-10 atol=1e-10*cancellation
    end
    # the factorizations go through the package's machinery: QR, which has
    # no in place refactorization nor an in place solve, agrees with KLU on
    # a junction in series with another on a node without capacitance,
    # whose stiffness refreshes the step's factorization, under either
    # rule, in the solve, the tangent and the adjoint
    series = Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "2", Capacitor(100e-15)),
        ("jj", "2", "3", JosephsonJunction(1e-9)), ("lj", "3", "0", JosephsonJunction(2e-9)),
        ("c2", "2", "0", Capacitor(1e-12))])
    qprob = transientproblem(series; sources = [TransientSource(1,
        t -> 2e-6*sinpi(2*3e9*t)*(t <= 0 ? 0.0 : t >= 0.3e-9 ? 1.0 : sinpi(t/0.6e-9)^2))])
    for method in (GaussLegendre(), Trapezoidal())
        klu = transientsolve(qprob, (0.0, 0.5e-9); dt = 2e-12, method, record = :phases)
        qr = transientsolve(qprob, (0.0, 0.5e-9); dt = 2e-12, method, record = :phases, factorization = QRfactorization())
        @test qr.stats.factorizations > 100
        @test qr.voltage ≈ klu.voltage rtol=1e-8
        currents = reshape([1e-8*sinpi(2*1.1e9*t) for t in klu.times], 1, :)
        weights = reshape([cospi(2*1.7e9*t) for t in klu.times], 1, :)
        @test transienttangent(qr, currents; factorization = QRfactorization()).outgoing ≈
            transienttangent(klu, currents).outgoing rtol=1e-8
        @test transientadjoint(qr, weights; factorization = QRfactorization()).currents ≈
            transientadjoint(klu, weights).currents rtol=1e-8
    end
    # a block factorization has no node blocks of modes to factorize in a
    # time step, and is refused where the step's factorization is chosen
    @test_throws ArgumentError transientsolve(qprob, (0.0, 0.5e-9); dt = 2e-12,
        factorization = BlockFactorization())
    klu = transientsolve(qprob, (0.0, 0.1e-9); dt = 2e-12, record = :phases)
    @test_throws ArgumentError transienttangent(klu,
        reshape([1e-8*sinpi(2*1.1e9*t) for t in klu.times], 1, :);
        factorization = BlockFactorization())
end

@testset "the sensitivity to the component values" begin
    rng = Random.default_rng()
    # one component of each kind the linearized solve differentiates:
    # a series capacitor, a junction, a capacitor and a resistor to
    # ground, and an inductor to a second port
    netlist = [("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)),
        ("Lj1", "2", "0", JosephsonJunction(1e-9)), ("C2", "2", "0", Capacitor(500e-15)), ("R2", "2", "0", Resistor(2000.0)),
        ("L1", "2", "3", Inductor(2e-9)), ("C3", "3", "0", Capacitor(300e-15)), ("P2", "3", "0", Port(2; Z0 = 50.0))]
    circuit = Circuit(netlist)
    names = ["C1", "Lj1", "C2", "R2", "L1"]
    drive(t) = 0.15e-6*sinpi(2*3e9*t) + 0.01e-6*sinpi(2*1.3e9*t)
    sources = [TransientSource(1, drive)]
    prob = transientproblem(circuit; sources)
    scaled(name, r) = transientproblem(Circuit(scaledentries(netlist, name, r)); sources)
    tspan, dt = (0.0, 1e-9), 2e-12
    tol = (; rtol = 1e-12, atol = 1e-13)
    weights = randn(rng, 2, 501)
    for method in (Trapezoidal(), BackwardEuler(), GaussLegendre())
        sol = transientsolve(prob, tspan; dt, method, record = :states, tol...)
        sens = transientsensitivity(sol, names)
        @test size(sens.outgoing) == (2, length(sol.times), length(names))
        # against a central difference of the solve in each value
        eps = 1e-4
        for (c, name) in enumerate(names)
            plus = transientsolve(scaled(name, 1 + eps), tspan; dt, method, tol...)
            minus = transientsolve(scaled(name, 1 - eps), tspan; dt, method, tol...)
            fd = (plus.outgoing .- minus.outgoing) ./ (2eps)
            @test sens.outgoing[:, :, c] ≈ fd rtol=1e-5 atol=1e-7*maximum(abs, fd)
        end
        # the adjoint's derivative of an objective is the weighted tangent
        adj = transientadjoint(sol, weights; components = names)
        @test adj.sensitivity ≈ [sum(weights .* sens.outgoing[:, :, c]) for c in eachindex(names)] rtol=1e-9
        # the junction alone reads the phases
        phases = transientsolve(prob, tspan; dt, method, record = :phases, tol...)
        @test transientsensitivity(phases, ["Lj1"]).outgoing ≈ sens.outgoing[:, :, 2] rtol=1e-9
        @test_throws ArgumentError transientsensitivity(phases, ["C1"])
    end
    # the Gauss-Legendre rule from its checkpoints, replaying the states
    cps = transientsolve(prob, tspan; dt, method = GaussLegendre(), record = :checkpoints, tol...)
    full = transientsolve(prob, tspan; dt, method = GaussLegendre(), record = :states, tol...)
    @test transientsensitivity(cps, names).outgoing ≈ transientsensitivity(full, names).outgoing rtol=1e-7
    @test transientadjoint(cps, weights; components = names).sensitivity ≈
        transientadjoint(full, weights; components = names).sensitivity rtol=1e-7
    # a batch, every condition on one pass
    batch = transientsolve([prob, transientproblem(prob; sources = [TransientSource(1, t -> drive(t)/2)])], tspan;
        dt, record = :states, tol...)
    sb = transientsensitivity(batch, names)
    ab = transientadjoint(batch, weights; components = names)
    for j in 1:2
        @test sb.outgoing[:, :, :, j] ≈ transientsensitivity(batch[j], names).outgoing rtol=1e-10
        @test ab.sensitivity[:, j] ≈ transientadjoint(batch[j], weights; components = names).sensitivity rtol=1e-9
    end
    @test_throws ArgumentError transientsensitivity(full, ["nothere"])
    @test_throws ArgumentError transientsensitivity(full, ["P1"])

    # a port's own termination moves the port's reference impedance
    # and conductance with it, which the port waves read directly
    let name = "p1/termination",
            portcircuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12))])
        psources = [TransientSource(1, t -> 1e-6*sinpi(2e9*t))]
        pprob = transientproblem(portcircuit; sources = psources)
        scaledport(r) = transientproblem(Circuit([(:p1, 1, 0, Port(1; Z0 = r*50.0)), (:c1, 1, 0, Capacitor(1e-12))]); sources = psources)
        for method in (Trapezoidal(), BackwardEuler(), GaussLegendre())
            psol = transientsolve(pprob, (0.0, 1e-9); dt = 2e-12, method, record = :states, tol...)
            ps = transientsensitivity(psol, [name])
            eps = 1e-5
            plus = transientsolve(scaledport(1 + eps), (0.0, 1e-9); dt = 2e-12, method, tol...)
            minus = transientsolve(scaledport(1 - eps), (0.0, 1e-9); dt = 2e-12, method, tol...)
            for q in (:voltage, :incident, :outgoing)
                fd = (getproperty(plus, q) .- getproperty(minus, q)) ./ (2eps)
                @test getproperty(ps, q)[:, :, 1] ≈ fd rtol=1e-5 atol=1e-7*maximum(abs, fd)
            end
            pweights = randn(rng, 1, length(psol.times))
            for q in (:voltage, :incident, :outgoing)
                padj = transientadjoint(psol, pweights; quantity = q, components = [name])
                @test padj.sensitivity ≈ [sum(pweights .* getproperty(ps, q)[:, :, 1])] rtol=1e-9
            end
        end
    end
    # two ports declared out of numerical order with different
    # impedances, each with its own branch and drive: the environments
    # and the traces are listed by port number, and each termination's
    # derivative lands on its own port's row
    twoport(z1, z2) = transientproblem(Circuit([(:p2, 1, 0, Port(2; Z0 = z2)), (:c2, 1, 0, Capacitor(1e-12)),
        (:p1, 2, 0, Port(1; Z0 = z1)), (:c1, 2, 0, Capacitor(2e-12))]);
        sources = [TransientSource(1, t -> 1e-6*sinpi(2e9*t)), TransientSource(2, t -> 0.5e-6*cospi(3e9*t))])
    tp = twoport(50.0, 75.0)
    @test tp.portimpedances == [50.0, 75.0]
    tnames = ["p2/termination", "p1/termination"]
    @test JC.componentperturbation(tp, tnames, JC.CPU(); forcing = true).ports == [2, 1]
    tsol = transientsolve(tp, (0.0, 1e-9); dt = 2e-12, method = GaussLegendre(), record = :states, tol...)
    ts = transientsensitivity(tsol, tnames)
    eps = 1e-5
    for (c, scaled) in enumerate((r -> twoport(50.0, r*75.0), r -> twoport(r*50.0, 75.0)))
        plus = transientsolve(scaled(1 + eps), (0.0, 1e-9); dt = 2e-12, method = GaussLegendre(), tol...)
        minus = transientsolve(scaled(1 - eps), (0.0, 1e-9); dt = 2e-12, method = GaussLegendre(), tol...)
        for q in (:voltage, :incident, :outgoing)
            fd = (getproperty(plus, q) .- getproperty(minus, q)) ./ (2eps)
            @test getproperty(ts, q)[:, :, c] ≈ fd rtol=1e-5 atol=1e-7*maximum(abs, fd)
        end
    end
    # an adjoint holds the entries of the components' derivatives and
    # never the stacked derivatives, whose rows are the state times the
    # components
    adjcp = JC.componentperturbation(tp, tnames, JC.CPU(); forcing = false)
    @test isnothing(adjcp.forcing)
    @test adjcp.entries.G.count <= 4*length(tnames) && adjcp.entries.C.count == 0
    @test size(JC.componentperturbation(tp, tnames, JC.CPU(); forcing = true).forcing.dG) == (2*length(tp), length(tp))
    # a junction's entries are read from its column of the incidence, so
    # an adjoint's perturbation of every junction of a chain allocates in
    # proportion to the chain: linear growth allocates four times as much
    # at four times the junctions, quadratic sixteen
    function perturbationbytes(m)
        chain = transientproblem(Circuit(vcat([("P1", "1", "0", Port(1; Z0 = 50.0))],
            [("J$(i)", "$(i)", "$(i + 1)", JosephsonJunction(100e-12)) for i in 1:m],
            [("C$(i)", "$(i)", "0", Capacitor(30e-15)) for i in 1:m])))
        junctions = ["J$(i)" for i in 1:m]
        run() = JC.componentperturbation(chain, junctions, JC.CPU(); forcing = false)
        run()
        return @allocated run()
    end
    @test perturbationbytes(512) < 6*perturbationbytes(128)

    # a junction on a node without capacitor or resistor to ground, whose
    # endpoint the Gauss-Legendre rule projects onto the constraint: the
    # projection and the reading carry the perturbation too
    L = 1e-9
    floating(r) = Circuit([("p", "1", "0", Port(1; termination = nothing)),
        ("jj", "1", "0", JosephsonJunction(r[1]*L)), ("l", "1", "2", Inductor(r[2]*2e-9)),
        ("c", "2", "0", Capacitor(r[3]*300e-15)), ("r", "2", "0", Resistor(r[4]*200.0)), ("p2", "2", "0", Port(2))])
    w = 2pi*1e9
    fsources = [TransientSource(1, t -> JC.phi0/L*0.2*(1 - cos(w*t)))]
    fprob = transientproblem(floating(ones(4)); sources = fsources)
    @test !isempty(fprob.algebraic)
    fnames = ["jj", "l", "c", "r"]
    fsol = transientsolve(fprob, (0.0, 1e-9); dt = 1e-9/80, method = GaussLegendre(), record = :states, tol...)
    fsens = transientsensitivity(fsol, fnames)
    eps = 1e-5
    for (c, name) in enumerate(fnames)
        rp, rm = ones(4), ones(4)
        rp[c] += eps; rm[c] -= eps
        plus = transientsolve(transientproblem(floating(rp); sources = fsources), (0.0, 1e-9); dt = 1e-9/80, method = GaussLegendre(), tol...)
        minus = transientsolve(transientproblem(floating(rm); sources = fsources), (0.0, 1e-9); dt = 1e-9/80, method = GaussLegendre(), tol...)
        fd = (plus.voltage .- minus.voltage) ./ (2eps)
        @test fsens.voltage[:, :, c] ≈ fd rtol=1e-5 atol=1e-7*maximum(abs, fd)
    end
    fweights = randn(rng, 2, length(fsol.times))
    fadj = transientadjoint(fsol, fweights; quantity = :voltage, components = fnames)
    @test fadj.sensitivity ≈ [sum(fweights .* fsens.voltage[:, :, c]) for c in eachindex(fnames)] rtol=1e-9
    fcps = transientsolve(fprob, (0.0, 1e-9); dt = 1e-9/80, method = GaussLegendre(), record = :checkpoints, tol...)
    @test transientsensitivity(fcps, fnames).voltage ≈ fsens.voltage rtol=1e-7
    @test transientadjoint(fcps, fweights; quantity = :voltage, components = fnames).sensitivity ≈ fadj.sensitivity rtol=1e-7
    # the junction alone, from the phases, through the projection
    fph = transientsolve(fprob, (0.0, 1e-9); dt = 1e-9/80, method = GaussLegendre(), record = :phases, tol...)
    @test transientsensitivity(fph, ["jj"]).voltage ≈ fsens.voltage[:, :, 1] rtol=1e-10
end

@testset "many tones, and their demodulation" begin
    # a linear RC has an analytic transfer function per tone; the
    # transient carries all the tones in one state
    frequencies = [1e9, 2.3e9, 3.7e9, 4.2e9, 5.1e9]
    amplitudes = [1.0, 0.08, 0.06, 0.04, 0.02]*1e-6
    current(t) = sum(a*cospi(2f*t) for (a, f) in zip(amplitudes, frequencies))
    prob = transientproblem(Circuit(rc); sources = [TransientSource(1, current)])
    sol = transientsolve(prob, (0.0, 12e-9); dt = 1e-12)
    # a Hann window after the start's transient, whose leakage from the
    # other tones is below the rule's error
    hann = t -> 2e-9 <= t <= 12e-9 ? sinpi((t - 2e-9)/10e-9)^2 : 0.0
    for (f, a) in zip(frequencies, amplitudes)
        measured = transientdemodulate(sol, 1, 2pi*f; quantity = :voltage, window = hann)
        @test measured ≈ 50a/(1 + 2pi*im*f*50e-12) rtol=1e-6
    end
    @test length(prob) == 1
end

@testset "the Gauss-Legendre rule: fourth order, the responses, the reuse" begin
    # the driven RC against a fine reference: the error falls sixteen
    # fold per halving where the trapezoidal rule's falls four fold
    smooth(t) = 1e-6*sinpi(t/1e-9)^2*sinpi(2*4e9*t)
    prob = transientproblem(Circuit(rc); sources = [TransientSource(1, smooth)])
    ref = transientsolve(prob, (0.0, 1e-9); dt = 1e-9/8192, method = GaussLegendre())
    errors = [maximum(abs.(transientsolve(prob, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre()).voltage[1, :] .-
        ref.voltage[1, 1:8192÷n:end])) for n in (64, 128)]
    @test 14 < errors[1]/errors[2] < 18
    # the lossless LC keeps its energy at twenty samples per period
    L, C, V = 1e-9, 1e-12, 1e-6
    lc = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)),
        ("c", "1", "0", Capacitor(C)), ("l", "1", "0", Inductor(L))]))
    period = 2pi*sqrt(L*C)
    sol = transientsolve(lc, (0.0, 3period); dt = period/20, method = GaussLegendre(),
        initialstate = transientstate(lc; voltage = [V]), record = :states)
    v = sol.rate[1, :] .* JC.phi0
    energy = C .* v .^ 2 ./ 2 .+ (sol.flux[1, :] .* JC.phi0) .^ 2 ./ (2L)
    @test maximum(abs.(energy ./ energy[1] .- 1)) < 1e-12
    @test maximum(abs.(v .- V .* cos.(sol.times ./ sqrt(L*C)))) < 3e-4V
    @test sol.stats.factorizations == 1
    # the coupled tanks through their auxiliary rows
    coupled = transientproblem(Circuit([("C1", "1", "0", Capacitor(C)), ("C2", "2", "0", Capacitor(C)),
        ("L1", "1", "0", Inductor(L)), ("L2", "2", "0", Inductor(L)), ("K1", "L1", "L2", MutualInductor(0.3))]))
    s2 = transientsolve(coupled, (0.0, 2period); dt = period/40, method = GaussLegendre(),
        initialstate = transientstate(coupled; voltage = [V, 0]), record = :states)
    wp, wm = 1/sqrt(L*C*1.3), 1/sqrt(L*C*0.7)
    @test maximum(abs.(s2.rate[1, :] .* JC.phi0 .- V/2 .* (cos.(wp .* s2.times) .+ cos.(wm .* s2.times)))) < 2e-5V
    # the pumped junction against harmonic balance at twenty samples per
    # period, where the trapezoidal rule is off by more than its value;
    # one complex factorization serves the whole solve
    circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)),
        ("Lj1", "2", "0", JosephsonJunction(1e-9)), ("C2", "2", "0", Capacitor(1e-12))])
    fp, ip = 4.75e9, 0.00565e-6
    ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
    pump(t) = 2ip*ramp(t)*cospi(2fp*t)
    pa = transientproblem(circuit; sources = [TransientSource(1, pump)])
    hb = hbnlsolve((2pi*fp,), (10,), [(mode = (1,), port = 1, current = ip)], circuit,
        Dict{Symbol,Float64}(); keyedarrays = false)
    expected = 2im*2pi*fp*JC.phi0*hb.nodeflux[1]
    window(t) = 200e-9 <= t <= 300e-9 ? sinpi((t - 200e-9)/100e-9)^2 : 0.0
    gauss = transientsolve(pa, (0.0, 300e-9); dt = 10e-12, method = GaussLegendre())
    @test transientdemodulate(gauss, 1, 2pi*fp; quantity = :voltage, window) ≈ expected rtol=2e-2
    @test gauss.stats.factorizations == 1
    trap = transientsolve(pa, (0.0, 300e-9); dt = 10e-12, method = Trapezoidal())
    @test !isapprox(transientdemodulate(trap, 1, 2pi*fp; quantity = :voltage, window), expected; rtol = 0.5)
    # the tangent against finite differences and the adjoint against the
    # tangent, on the full stage equations with their two stiffnesses
    loaded = Circuit([("p1", "1", "0", Port(1)), ("p2", "2", "1", Port(2; Z0 = 75.0)),
        ("c1", "1", "0", Capacitor(1e-12)), ("c2", "2", "0", Capacitor(1e-12)),
        ("jj", "2", "1", JosephsonJunction(1e-9)), ("l", "2", "0", Inductor(2e-9))])
    drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(2*3e9*t)
    probe(t, p) = 1e-8*sinpi(2*1.1e9*t + p)
    lp = transientproblem(loaded; sources = [TransientSource(1, drive), TransientSource(2, 2e-8)])
    dt, T = 4e-12, 0.5e-9
    rec = transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12, method = GaussLegendre())
    @test size(rec.phases) == (1, 2, length(rec.times))
    currents = [probe(t, p) for p in 1:2, t in rec.times]
    tg = transienttangent(rec, currents)
    eps = 1e-4
    plus = transientproblem(loaded; sources = [TransientSource(1, t -> drive(t) + eps*probe(t, 1)),
        TransientSource(2, t -> 2e-8 + eps*probe(t, 2))])
    minus = transientproblem(loaded; sources = [TransientSource(1, t -> drive(t) - eps*probe(t, 1)),
        TransientSource(2, t -> 2e-8 - eps*probe(t, 2))])
    sp = transientsolve(plus, (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
    sm = transientsolve(minus, (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
    @test tg.outgoing ≈ (sp.outgoing .- sm.outgoing) ./ (2eps) rtol=1e-6
    weights = [cospi(2*1.7e9*t + p) for p in 1:2, t in rec.times]
    ad = transientadjoint(rec, weights)
    @test sum(weights .* tg.outgoing) ≈ sum(ad.currents .* currents) rtol=1e-10
    dx0, dv0 = [0.3, -0.2], [1e9, 2e9]
    tg0 = transienttangent(rec, zeros(2, length(rec.times)); initialstate = (dx0, dv0))
    @test sum(weights .* tg0.outgoing) ≈ dot(ad.initialflux, dx0) + dot(ad.initialrate, dv0) rtol=1e-9
    # the tangent takes each stage at its own stiffness: a junction on a
    # node without capacitance, pulsed hard at a coarse step, whose solve
    # refreshes at nearly every step, and its tangent matches the solve's
    # differences
    pulsed = Circuit([(:p, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(1e-12)), (:l, 1, 2, Inductor(0.5e-9)),
        (:lj, 2, 0, JosephsonJunction(1e-9))])
    hard(e) = [TransientSource(1, t -> 4e-6*sinpi(8e9*t)*(t <= 0 || t >= 1e-9 ? 0.0 : sinpi(t/1e-9)^2) + e*1e-9*sinpi(t/2e-9)^2)]
    hs = transientsolve(transientproblem(pulsed; sources = hard(0.0)), (0.0, 2e-9); dt = 2e-12, record = :phases)
    @test hs.stats.factorizations > hs.stats.steps/2
    htg = transienttangent(hs, reshape(1e-9 .* sinpi.(hs.times ./ 2e-9) .^ 2, 1, :))
    hfd = (transientsolve(transientproblem(pulsed; sources = hard(1e-3)), (0.0, 2e-9); dt = 2e-12, rtol = 1e-12, iterations = 60).voltage .-
        transientsolve(transientproblem(pulsed; sources = hard(-1e-3)), (0.0, 2e-9); dt = 2e-12, rtol = 1e-12, iterations = 60).voltage) ./ 2e-3
    @test htg.voltage ≈ hfd rtol=1e-5
    # a sink receives every column of the currents once it is final,
    # under both rules, and the adjoint then stores none
    for (s, a) in ((rec, ad), (transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12, method = Trapezoidal()), nothing))
        stored = isnothing(a) ? transientadjoint(s, weights) : a
        received = zeros(size(stored.currents))
        seen = Int[]
        streamed = transientadjoint(s, weights; sink = (k, values) -> (push!(seen, k); received[:, k] .= values[:, 1]; nothing))
        @test isnothing(streamed.currents)
        @test seen == collect(length(s.times):-1:1)
        @test received ≈ stored.currents rtol=1e-12
        @test streamed.initialflux ≈ stored.initialflux
    end
    # the responses need the stage phases, and the reuse carries the
    # complex factorization across the solve and the responses
    bare = transientsolve(lp, (0.0, T); dt, method = GaussLegendre())
    @test_throws ArgumentError transienttangent(bare, currents)
    @test_throws ArgumentError transientadjoint(bare, weights)
    reuse = TransientReuse()
    r1 = transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12, method = GaussLegendre(), reuse)
    @test reuse.system.method isa GaussLegendre
    @test transienttangent(r1, currents; reuse).outgoing ≈ tg.outgoing rtol=1e-10
    @test_throws ArgumentError transientsolve(lp, (0.0, T); dt, method = GaussLegendre(), linearsolver = GMRES())
end

@testset "the algebraic directions of the Gauss-Legendre rule" begin
    # a junction alone across an unterminated port: the flux is the
    # constraint's own solution and the voltage its derivative, which
    # the projected endpoint holds to the Newton tolerance and the
    # reading of the rate from the differentiated constraint holds to
    # roundoff at every step, where the plain rule's rate would not
    # converge
    L, w, a = 1e-9, 2pi*1e9, 0.2
    jj = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)),
        ("jj", "1", "0", JosephsonJunction(L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
    @test jj.inertialess == [[1]] && jj.algebraic == [[1]]
    errors = map((20, 40, 80)) do n
        s = transientsolve(jj, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
        phase = asin.(a .* (1 .- cos.(w .* s.times)))
        voltage = JC.phi0 .* (a*w .* sin.(w .* s.times)) ./ sqrt.(1 .- (a .* (1 .- cos.(w .* s.times))) .^ 2)
        (maximum(abs.(s.flux[1, :] .- phase)), maximum(abs.(s.voltage[1, :] .- voltage))/(JC.phi0*w))
    end
    @test all(e -> e[1] < 1e-10, errors)
    @test all(e -> e[2] < 1e-11, errors)
    # The classification reads the rate system of the equations, not a
    # graph of the ports. An open block on the junction's node adds no
    # conductance: the node stays algebraic and the solution and its
    # order are unchanged (a port taken for a resistor lost the
    # constraint and left a 0.42 error at every step). A short block
    # grounds the node: no direction, and no voltage. A through to a
    # second node with an inductor joins the two into one algebraic
    # direction, the junction and the inductor in parallel; and a
    # rational inductor block on the bare junction's node equals the
    # explicit inductor there.
    jjo = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
        ("open", "1", ScatteringParameters(ones(1, 1)))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
    @test jjo.algebraic == [[1]] && jjo.inertialess == [[1], [2]]
    @test size(jjo.constraints) == (1, 2) && abs(jjo.constraints[1, 1]) > 0.1 && abs(jjo.constraints[1, 2]) > 0.01
    oerrors = map((20, 40, 80)) do n
        s = transientsolve(jjo, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
        voltage = JC.phi0 .* (a*w .* sin.(w .* s.times)) ./ sqrt.(1 .- (a .* (1 .- cos.(w .* s.times))) .^ 2)
        maximum(abs.(s.voltage[1, :] .- voltage))/(JC.phi0*w)
    end
    @test all(<(1e-11), oerrors)
    jjs = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
        ("short", "1", ScatteringParameters(-ones(1, 1)))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
    @test isempty(jjs.algebraic) && jjs.inertialess == [[1], [2]]
    ss = transientsolve(jjs, (0.0, 1e-9); dt = 1e-9/40, method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
    @test maximum(abs, ss.voltage) < 1e-20
    jjt = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
        ("through", "1", "2", ScatteringParameters([0.0 1.0; 1.0 0.0])), ("l", "2", "0", Inductor(L))]);
        sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
    jjp = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
        ("l", "1", "0", Inductor(L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
    @test jjt.algebraic == [[1, 2]] && jjp.algebraic == [[1]]
    st = transientsolve(jjt, (0.0, 1e-9); dt = 1e-9/40, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
    sp = transientsolve(jjp, (0.0, 1e-9); dt = 1e-9/40, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
    @test st.voltage ≈ sp.voltage rtol=1e-8
    @test st.flux[1, :] ≈ sp.flux[1, :] rtol=1e-8
    @test st.flux[2, :] ≈ sp.flux[1, :] rtol=1e-8
    ra = 2*50.0/L
    rind = RationalScattering(fill(-ra, 1, 1), reshape([1.0, -1.0], 1, 2), reshape(-ra .* [1.0, -1.0], 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
    jjr = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
        ("lb", "1", "2", rind), ("l2", "2", "0", Inductor(L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
    jjp2 = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
        ("l", "1", "0", Inductor(2L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
    # the block's ports are open at infinite frequency, so both nodes
    # are algebraic, their constraints carrying the block's resting
    # waves; the block's states are stepped from the stage rates along
    # those directions, so the agreement converges at second order
    @test jjr.algebraic == [[1], [2]]
    rdiffs = map((160, 320, 640)) do n
        sr = transientsolve(jjr, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
        sp2 = transientsolve(jjp2, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
        maximum(abs, sr.voltage .- sp2.voltage)/maximum(abs, sp2.voltage)
    end
    @test 3.5 < rdiffs[1]/rdiffs[2] < 4.5 && 3.5 < rdiffs[2]/rdiffs[3] < 4.5 && rdiffs[3] < 2e-5
    # a driven inductor without capacitance: a linear constraint, met by
    # one correction, and the voltage read to roundoff
    ind = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("l", "1", "0", Inductor(L))]);
        sources = [TransientSource(1, t -> 1e-6*(1 - cos(w*t)))])
    @test ind.algebraic == [[1]]
    verrors = map((20, 40)) do n
        s = transientsolve(ind, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
        maximum(abs.(s.voltage[1, :] .- L*1e-6*w .* sin.(w .* s.times)))/(L*1e-6*w)
    end
    @test all(<(1e-11), verrors)
    # the index one unknowns of the endpoint are read from their
    # equations: a terminated port on an inductor, whose node has no
    # capacitance, keeps fourth order in its voltage, where the
    # rule's endpoint formula alone would give second, and so does a
    # block's port current; the tangent and the adjoint follow the
    # reading
    indport = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l, 1, 2, Inductor(1e-9)), (:c, 2, 0, Capacitor(1e-12)),
        (:jj, 2, 0, JosephsonJunction(1e-9))])
    blockport = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:att, 1, 2, ScatteringParameters([0.0 0.5; 0.5 0.0]; zref = 50.0)),
        (:l, 2, 3, Inductor(1e-9)), (:c, 3, 0, Capacitor(1e-12)), (:jj, 3, 0, JosephsonJunction(1e-9))])
    pulse(t) = t <= 0 ? 0.0 : 1e-6*sinpi(t/1e-9)^2*sinpi(2*4e9*t)
    for (c, resistive, aux) in ((indport, [[1]], Int[]), (blockport, [[1], [2]], [4, 5]))
        p = transientproblem(c; sources = [TransientSource(1, pulse)])
        pr = JC.transientsystem(p, 1e-12, GaussLegendre(), JC.CPU(), KLUfactorization()).projection
        @test pr.readrows == resistive && pr.auxrows == aux && isempty(pr.directions)
        ref = transientsolve(p, (0.0, 1e-9); dt = 1e-9/4096, method = GaussLegendre(), rtol = 1e-13, atol = 1e-15)
        errs = [maximum(abs.(transientsolve(p, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), rtol = 1e-13, atol = 1e-15).voltage[1, :] .-
            ref.voltage[1, 1:4096÷n:end])) for n in (64, 128)]
        @test 12 < errs[1]/errs[2] < 20
    end
    ip = transientproblem(indport; sources = [TransientSource(1, pulse)])
    irec = transientsolve(ip, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre())
    iprobe(t) = 1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + 1)
    ieps = 1e-4
    icurrents = [iprobe(t) for _ in 1:1, t in irec.times]
    itg = transienttangent(irec, icurrents)
    ishift(s) = transientproblem(indport; sources = [TransientSource(1, t -> pulse(t) + s*ieps*iprobe(t))])
    isp = transientsolve(ishift(1), (0.0, 0.5e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre())
    ism = transientsolve(ishift(-1), (0.0, 0.5e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre())
    @test itg.voltage ≈ (isp.voltage .- ism.voltage) ./ (2ieps) rtol=1e-5
    iweights = [cospi(2*1.7e9*t) for _ in 1:1, t in irec.times]
    iad = transientadjoint(irec, iweights; quantity = :voltage)
    @test sum(iweights .* itg.voltage) ≈ sum(iad.currents .* icurrents) rtol=1e-9
    icps = transientsolve(ip, (0.0, 0.5e-9); dt = 2e-12, record = :checkpoints, checkpointevery = 16, rtol = 1e-12, method = GaussLegendre())
    @test transientadjoint(icps, iweights; quantity = :voltage).currents ≈ iad.currents rtol=1e-7
    # the coupled tanks and the lossless LC project nothing
    @test isnothing(JC.transientsystem(transientproblem(Circuit([("C1", "1", "0", Capacitor(1e-12)), ("C2", "2", "0", Capacitor(1e-12)),
        ("L1", "1", "0", Inductor(L)), ("L2", "2", "0", Inductor(L)), ("K1", "L1", "L2", MutualInductor(0.3))])), 1e-12, GaussLegendre(), JC.CPU(),
        KLUfactorization()).projection)
    # the tangent against finite differences and the adjoint against
    # the tangent through the projection, on the record, on the
    # checkpoints and on a batch: a junction from a driven capacitive
    # node to a node with only an inductor
    c = Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "0", Capacitor(1e-12)), ("jj", "1", "2", JosephsonJunction(1e-9)),
        ("l2", "2", "0", Inductor(2e-9)), ("p2", "2", "0", Port(2; termination = nothing))])
    drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(2*3e9*t)
    probe(t, k) = 1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k)
    lp = transientproblem(c; sources = [TransientSource(1, drive), TransientSource(2, t -> 0.0)])
    sys = JC.transientsystem(lp, 2e-12, GaussLegendre(), JC.CPU(), KLUfactorization())
    @test sys.projection.directions == [[2]] && sys.projection.pj == [1]
    dt, T = 2e-12, 0.5e-9
    rec = transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12, method = GaussLegendre())
    @test size(rec.endphases) == (1, length(rec.times))
    currents = [probe(t, k) for k in 1:2, t in rec.times]
    tg = transienttangent(rec, currents)
    eps = 1e-4
    shifted(s) = transientproblem(c; sources = [TransientSource(1, t -> drive(t) + s*eps*probe(t, 1)),
        TransientSource(2, t -> s*eps*probe(t, 2))])
    sp = transientsolve(shifted(1), (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
    sm = transientsolve(shifted(-1), (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
    @test tg.outgoing ≈ (sp.outgoing .- sm.outgoing) ./ (2eps) rtol=1e-5
    @test tg.voltage ≈ (sp.voltage .- sm.voltage) ./ (2eps) rtol=3e-5
    weights = [cospi(2*1.7e9*t + k) for k in 1:2, t in rec.times]
    ad = transientadjoint(rec, weights)
    @test sum(weights .* tg.outgoing) ≈ sum(ad.currents .* currents) rtol=1e-9
    # every direction converges against its own right hand side: a
    # direction sixteen orders weaker than its companion propagates as
    # it does alone
    dirs = zeros(2, length(rec.times), 2)
    dirs[1, :, 1] .= currents[1, :]
    dirs[2, :, 2] .= 1e-16 .* currents[2, :]
    both = transienttangent(rec, dirs)
    @test both.outgoing[:, :, 2] ≈ transienttangent(rec, dirs[:, :, 2]).outgoing rtol=1e-8
    # a current given at the grid and the stage times is read as it
    # is, and agrees with the stencil's reading of the smooth grid
    # current to the stencil's interpolation
    cs = JC.gausscoefficients().c
    staged = [probe(t + (s == 1 ? 0.0 : cs[s - 1]*dt), k) for k in 1:2, s in 1:3, t in rec.times, _ in 1:1]
    @test transienttangent(rec, staged).outgoing ≈ tg.outgoing rtol=1e-5
    @test transienttangent(rec, staged).outgoing ≈ (sp.outgoing .- sm.outgoing) ./ (2eps) rtol=1e-5
    cps = transientsolve(lp, (0.0, T); dt, record = :checkpoints, checkpointevery = 16, rtol = 1e-12, method = GaussLegendre())
    @test transienttangent(cps, currents).outgoing ≈ tg.outgoing rtol=1e-7
    @test transientadjoint(cps, weights).currents ≈ ad.currents rtol=1e-7
    # a reuse keeps the stepper the checkpoints are replayed on, with
    # its factorizations, across the responses of one shape
    reuse = TransientReuse()
    tc1, ac1 = transienttangent(cps, currents; reuse), transientadjoint(cps, weights; reuse)
    kept = (only(reuse.tangent).replay, only(reuse.adjoint).replay)
    @test kept[1] isa JC.GaussStepper && kept[2] isa JC.GaussStepper
    tc2, ac2 = transienttangent(cps, currents; reuse), transientadjoint(cps, weights; reuse)
    @test only(reuse.tangent).replay === kept[1] && only(reuse.adjoint).replay === kept[2]
    @test tc2.outgoing == tc1.outgoing && ac2.currents == ac1.currents
    @test isnothing(JC.replaystepper!(only(reuse.tangent), JC.batchof(rec), reuse.system))
    half = transientproblem(lp; sources = [TransientSource(1, t -> drive(t)/2), TransientSource(2, t -> 0.0)])
    b = transientsolve([lp, half], (0.0, T); dt, record = :phases, rtol = 1e-12)
    @test b.voltage[:, :, 1] == rec.voltage
    @test size(b.endphases) == (1, length(rec.times), 2)
    tb = transienttangent(b, currents)
    @test tb.outgoing[:, :, 1] == tg.outgoing
    @test tb.outgoing[:, :, 2] ≈ transienttangent(b[2], currents).outgoing rtol=1e-9
    @test transientadjoint(b, weights).currents[:, :, 1] == ad.currents
end

# Many subnetworks without capacitance to ground: the classification and
# the projection hold them block by block of their couplings, the build
# and a step growing with the circuit rather than with its square or cube.
# A chain of floating capacitors, each driven through by the same current
# from the sources joining it to its neighbours, has a projected direction
# for every branch, and charges each capacitor to the integral of the
# current, `I0 (1 - cos wt)/(w C)`. The build and a short solve of the
# floating chain, of a junction chain with a lossy block at every node, a
# port current at every node, and of a junction chain without capacitance,
# a projected direction at every inner node, allocate in proportion to the
# chain.
@testset "many inertialess subnetworks" begin
    I0, w, C = 1e-6, 2pi*5e9, 1e-12
    function floatingchain(m)
        c = vcat([("C$(r)", "$(2r - 1)", "$(2r)", Capacitor(C)) for r in 1:m],
            [("I$(r)", r == 0 ? "0" : "$(2r)", r == m ? "0" : "$(2r + 1)", CurrentSource(0.0)) for r in 0:m])
        return Circuit(c), [TransientSource(Symbol("I$(r)"), t -> I0*sin(w*t)) for r in 0:m]
    end
    c, s = floatingchain(16)
    p = transientproblem(c; sources = s)
    @test length(p.algebraic) == 16
    sol = transientsolve(p, (0.0, 200e-12); dt = 1e-12, record = :states)
    charge = I0 .* (1 .- cos.(w .* sol.times)) ./ (w*C)
    # a node's state is its index among the circuit's node names less the ground's
    node = name -> findfirst(==(name), p.circuit.nodenames) - 1
    a, b = [node("$(2r - 1)") for r in 1:16], [node("$(2r)") for r in 1:16]
    voltage = JC.phi0 .* (sol.rate[a, :] .- sol.rate[b, :])
    @test maximum(abs.(voltage .- transpose(charge))) < 1e-6*maximum(charge)
    function blockchain(m)
        c = Any[("P1", "1", "0", Port(1; Z0 = 50.0))]
        for i in 1:m
            push!(c, ("J$(i)", "$(i)", "$(i + 1)", JosephsonJunction(100e-12)), ("C$(i)", "$(i)", "0", Capacitor(30e-15)),
                ("B$(i)", "$(i)", ScatteringParameters(fill(0.5, 1, 1); zref = 50.0)))
        end
        push!(c, ("Cend", "$(m + 1)", "0", Capacitor(30e-15)))
        return Circuit(c), [TransientSource(1, t -> 1e-7*sinpi(1e10*t))]
    end
    function bytes(make, m)
        c, s = make(m)
        run() = transientsolve(transientproblem(c; sources = s), (0.0, 2e-12); dt = 1e-12)
        run()
        return @allocated run()
    end
    for make in (floatingchain, blockchain)
        @test bytes(make, 512) < 3*bytes(make, 256)
    end
    function junctionchain(m)
        c = vcat([("P1", "1", "0", Port(1; Z0 = 50.0)), ("P2", "$(m + 1)", "0", Port(2; Z0 = 50.0))],
            [("J$(i)", "$(i)", "$(i + 1)", JosephsonJunction(100e-12)) for i in 1:m])
        return Circuit(c), [TransientSource(1, t -> 1e-7*sinpi(1e10*t))]
    end
    # growth beyond linear is small beside the linear part at these sizes,
    # so a fourfold step tells them apart: linear growth allocates four
    # times as much, quadratic sixteen
    @test bytes(junctionchain, 2048) < 6*bytes(junctionchain, 512)
    # The projection keeps each condition's factorization of its Jacobian
    # across corrections and steps, a chord on it, refreshed where a
    # correction contracts poorly: a lattice of junctions without
    # capacitance, whose endpoint every step corrects, steps 40 times on
    # one factorization of its projection, where one a correction is one a
    # step. A replay from checkpoints sets the kept factorizations at each
    # checkpoint as the solve did, so each window ends exactly on the
    # checkpoint after it.
    function junctionlattice(m)
        site(i, j) = "n$(i)_$(j)"
        elements = Any[("P1", site(1, 1), "0", Port(1; Z0 = 50.0)), ("P2", site(m, m), "0", Port(2; Z0 = 50.0))]
        for i in 1:m, j in 1:m
            i < m && push!(elements, ("Jv$(i)_$(j)", site(i, j), site(i + 1, j), JosephsonJunction(100e-12)))
            j < m && push!(elements, ("Jh$(i)_$(j)", site(i, j), site(i, j + 1), JosephsonJunction(100e-12)))
        end
        return Circuit(elements)
    end
    lp = transientproblem(junctionlattice(4); sources = [TransientSource(1, t -> 2e-6*sinpi(1e10*t))])
    lsys = JC.transientsystem(lp, 1e-12, GaussLegendre(), JC.CPU(), JC.KLUfactorization())
    st = JC.gaussstepper(lsys, [lp], 1e-9, 1e-10, 15, JC.gaussbatchfactor(lsys, 1))
    JC.setstate!(st, zeros(length(lp), 1), zeros(length(lp), 1), nothing)
    for k in 1:40
        JC.advance!(st, (k - 1)*1e-12, k*1e-12, k)
    end
    @test st.pw.factorizations[1] <= 4
    cps = transientsolve([lp], (0.0, 120e-12); dt = 1e-12, record = :checkpoints, checkpointevery = 20)
    rst = JC.replaystepper(lsys, cps)
    windows = JC.responsewindows(cps, lsys, true; stepper = rst)
    @test all(reverse(1:length(windows) - 1)) do k
        windows[k].replay()
        rst.x == view(cps.checkpoints.flux, :, k + 1, :)
    end
end

# A cascade of blocks through nodes without capacitance is one block of
# the rate system, too large to decompose densely, and is factorized
# sparse. A cascade of lossy throughs transmits the product of their
# transmissions, as one through of that transmission does; the build and
# a short solve of the cascade allocate in proportion to it. On a chain
# with a dependent row and a dependent column, the factorization's null
# spaces and minimum norm solutions are those of the dense decomposition
# of the same block, plain and transposed.
@testset "a cascade of blocks" begin
    through(t) = ScatteringParameters([0.0 t; t 0.0]; zref = 50.0)
    function cascade(ts)
        c = Any[("P1", "1", "0", Port(1; Z0 = 50.0))]
        for (i, t) in enumerate(ts)
            push!(c, ("T$(i)", "$(i)", "$(i + 1)", through(t)))
        end
        k = length(ts) + 1
        push!(c, ("C", "$(k)", "0", Capacitor(1e-12)), ("P2", "$(k)", "0", Port(2; Z0 = 50.0)))
        return transientproblem(Circuit(c); sources = [TransientSource(1, t -> 1e-7*sinpi(1e10*t))])
    end
    waves(p) = transientsolve(p, (0.0, 100e-12); dt = 1e-12).outgoing
    @test isapprox(waves(cascade(fill(0.99, 40))), waves(cascade([0.99^40])); rtol = 1e-10)
    function bytes(m)
        run() = waves(cascade(fill(0.99, m)))
        run()
        return @allocated run()
    end
    @test bytes(512) < 3*bytes(256)
    # the planted chain: the islands are its first 60 unknowns and the
    # block rows the rest, so the rate system is the chain itself
    n, k0 = 120, 60
    rng = Random.Xoshiro(3)
    K = spdiagm(0 => 1 .+ rand(rng, n), 1 => rand(rng, n - 1), -1 => rand(rng, n - 1))
    K[:, 31] = K[:, 30]
    K[90, :] = K[91, :]
    dropzeros!(K)
    G, L = hcat(K[:, 1:k0], spzeros(n, n - k0)), hcat(spzeros(n, k0), K[:, k0 + 1:n])
    Z0 = sparse(1:k0, 1:k0, ones(k0), n, k0)
    Ea = sparse(k0 + 1:n, 1:n - k0, ones(n - k0), n, n - k0)
    dense = JC.ratesystem(G, L, Z0, Ea; maxdense = n)
    factored = JC.ratesystem(G, L, Z0, Ea; maxdense = 16)
    @test isempty(dense.pinv.factors) && length(factored.pinv.factors) == 1
    # the null vectors of every block as columns over the whole system,
    # compared by the projections onto their spans
    function nullspan(rs, right)
        vs = [(v = zeros(n); v[right ? b.cols : b.rows] = (right ? b.rightnull : b.leftnull)[:, j]; v)
            for b in rs.blocks for j in axes(right ? b.rightnull : b.leftnull, 2)]
        B = stack(vs; dims = 2)
        return B*pinv(B)
    end
    for right in (true, false)
        @test nullspan(dense, right) ≈ nullspan(factored, right) atol = 1e-8
    end
    g = randn(rng, n, 3)
    for transposed in (false, true)
        a = JC.ratesolve!(zeros(n, 3), dense.pinv, g, JC.ratework(dense.pinv, 3); transposed)
        b = JC.ratesolve!(zeros(n, 3), factored.pinv, g, JC.ratework(factored.pinv, 3); transposed)
        @test a ≈ b rtol = 1e-8
    end
end

@testset "a batch of drive conditions" begin
    # the amplifier under three pump amplitudes as one solve: each
    # member equals its own solve, the responses run on a member, and
    # the problems must share the circuit and the targets
    circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)),
        ("Lj1", "2", "0", JosephsonJunction(1e-9)), ("C2", "2", "0", Capacitor(1e-12))])
    fp = 4.75e9
    ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
    base = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
    make(ip) = transientproblem(base; sources = [TransientSource(1, let ip = ip; t -> 2ip*ramp(t)*cospi(2fp*t); end)])
    problems = make.([0.002e-6, 0.004e-6, 0.00565e-6])
    @test problems[2].circuit === base.circuit && problems[2].matrices === base.matrices
    batch = transientsolve(problems, (0.0, 5e-9); dt = 5e-12, record = :states)
    @test batch isa TransientBatchSolution && length(batch) == 3
    @test size(batch.voltage) == (1, 1001, 3) && size(batch.phases) == (1, 2, 1001, 3)
    weights = [cospi(2*4.74e9*t) for p in 1:1, t in batch.times]
    for j in 1:3
        single = transientsolve(problems[j], (0.0, 5e-9); dt = 5e-12, method = GaussLegendre(), record = :states)
        member = batch[j]
        @test member.voltage ≈ single.voltage rtol=1e-12
        @test member.flux ≈ single.flux rtol=1e-12
        @test member.phases ≈ single.phases rtol=1e-12
        @test member.finalrate ≈ single.finalrate rtol=1e-12
        @test transientadjoint(member, weights).currents ≈ transientadjoint(single, weights).currents rtol=1e-10
        # a condition taken as a batch of one is a view of the batch's
        # arrays, of the type a range of conditions is
        @test typeof(JosephsonCircuits.batchof(member).phases) === typeof(batch[j:j].phases)
        @test JosephsonCircuits.batchof(member).phases == batch[j:j].phases
    end
    @test batch.stats.factorizations == 1
    # A host batch splits its conditions across the threads of the
    # session and steps the chunks at once. The split is exact: the
    # conditions are independent, so any layout of chunks fills the
    # batch's arrays with the same bits. Only the Newton counter
    # differs, since a joint step iterates until the worst condition of
    # its own chunk has converged, never more than the whole batch.
    let sys = JosephsonCircuits.transientsystem(base, 5e-12, GaussLegendre(),
            JosephsonCircuits.CPU(), JosephsonCircuits.transientfactorization(JosephsonCircuits.CPU())),
        states = [transientstate(q) for q in problems],
        run = chunks -> JosephsonCircuits.gaussbatchintegrate(sys, problems, 0.0, 5e-9, 1000,
            states, 1, :states, 0, 1e-9, 1e-10, 15, nothing; chunks = chunks)
        whole, split = run([1:3]), run([1:1, 2:3])
        for f in (:voltage, :incident, :outgoing, :phases, :flux, :rate, :finalflux, :initialflux)
            @test getfield(whole, f) == getfield(split, f)
        end
        @test whole.stats.factorizations == split.stats.factorizations == 1
        @test split.stats.newtoncorrections <= whole.stats.newtoncorrections
        @test JosephsonCircuits.batchchunks(JosephsonCircuits.CPU(), 1) == [1:1]
        # a device keeps the one chunk its uniform batch already is
        @test JosephsonCircuits.batchchunks(NotTheHost(), 8) == [1:8]
    end
    # A condition refreshes its own factorization: a junction biased
    # hard enough that its frozen operator drifts beside one that
    # never does, stepped as one chunk, as two and each alone, has the
    # same bits every way, forward, in the tangent and in the adjoint,
    # and the weak one factorizes once.
    let circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("Lj1", "1", "0", JosephsonJunction(1e-9)),
            ("C1", "1", "0", Capacitor(0.2e-12))]),
        ramp = t -> t <= 0 ? 0.0 : t >= 0.5e-9 ? 1.0 : (1 - cospi(t/0.5e-9))/2,
        base = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)]),
        pair = [transientproblem(base; sources = [TransientSource(1, t -> ib*ramp(t) + 1e-9*sinpi(2*3e9*t))])
            for ib in (0.2e-6, 0.6e-6)],
        sys = JC.transientsystem(base, 15e-12, GaussLegendre(), JC.CPU(), JC.transientfactorization(JC.CPU())),
        run = (js, chunks) -> JC.gaussbatchintegrate(sys, pair[js], 0.0, 333*15e-12, 333,
            [transientstate(q) for q in pair[js]], 1, :states, 0, 1e-9, 1e-10, 15, nothing; chunks)
        whole, split = run(1:2, [1:2]), run(1:2, [1:1, 2:2])
        alone = [run(j:j, [1:1]) for j in 1:2]
        @test alone[1].stats.factorizations == 1 && alone[2].stats.factorizations > 10
        @test whole.stats.factorizations == split.stats.factorizations == alone[2].stats.factorizations
        member = (a, j) -> selectdim(a, ndims(a), j)
        for j in 1:2, f in (:voltage, :flux, :phases)
            @test member(getfield(whole, f), j) == member(getfield(split, f), j) == member(getfield(alone[j], f), 1)
        end
        nt = length(whole.times)
        (injh, tp) = JC.targetinjection(base, JC.porttargets(base))
        tc = JC.tangentcurrents([1e-9*sinpi(2*2.9e9*t) for _ in 1:1, t in whole.times], 1, nt, 1)
        init = JC.tangentinitial(nothing, length(base), 1, 2, 0, 1, 0)
        wh = JC.adjointweights([cospi(2*3.1e9*t) for _ in 1:1, t in whole.times], 1, nt)
        tangent = chunks -> JC.gaussbatchtangent(whole, tc, injh, tp, init, sys, nothing, nothing, nothing; chunks)
        adjoint = chunks -> JC.gaussbatchadjoint(whole, wh, true, :outgoing, injh, tp, sys, nothing, nothing, nothing,
            nothing; chunks)
        @test tangent([1:2]).outgoing == tangent([1:1, 2:2]).outgoing
        @test adjoint([1:2]).currents == adjoint([1:1, 2:2]).currents
    end
    # so does the projection of its endpoint: two unequal junctions in
    # series around a node without capacitance, one pair driven into
    # the voltage state beside one which is not
    let Ic = JC.phi0/1e-9,
        ramp = t -> t <= 0 ? 0.0 : t >= 0.2e-9 ? 1.0 : (1 - cospi(t/0.2e-9))/2,
        base = transientproblem(Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "0", Capacitor(1e-12)),
            ("J1", "1", "2", JosephsonJunction(1e-9)), ("J2", "2", "0", JosephsonJunction(1.3e-9))]);
            sources = [TransientSource(1, t -> 0.0)]),
        pair = [transientproblem(base; sources = [TransientSource(1, t -> a*Ic*ramp(t) + 0.05Ic*sinpi(2*7e9*t))])
            for a in (0.3, 1.2)],
        sys = JC.transientsystem(base, 1e-12, GaussLegendre(), JC.CPU(), JC.transientfactorization(JC.CPU())),
        run = (js, chunks) -> JC.gaussbatchintegrate(sys, pair[js], 0.0, 1e-9, 1000,
            [transientstate(q) for q in pair[js]], 1, :states, 0, 1e-9, 1e-10, 15, nothing; chunks)
        whole = run(1:2, [1:2])
        @test !isempty(sys.projection.directions)
        @test all(j -> whole.flux[:, :, j] == run(j:j, [1:1]).flux[:, :, 1], 1:2)
    end
    # a step which fails throws the typed error naming the conditions
    # which failed it, the others having converged there: a pump the
    # step cannot follow beside a weak drive, which converges alone
    let base = transientproblem(JC.warmupcircuit(50.0, 100.0e-15, 1000.0e-12, 1000.0e-15);
            sources = [TransientSource(1, t -> 0.0)]),
        h = 1/(2*4.75e9),
        weak = transientproblem(base; sources = [TransientSource(1, t -> 1e-9*cospi(2*4.75e9*t))]),
        strong = transientproblem(base; sources = [TransientSource(1, t -> 1e-6*cospi(2*4.75e9*t))])
        err = try
            transientsolve([weak, strong], (0.0, 40h); dt = h)
        catch e
            e
        end
        @test err isa TransientStepError && err.conditions == [2] && err.cause == :newton
        @test err.step == 2 && err.time ≈ 2h
        @test transientsolve(weak, (0.0, 40h); dt = h).stats.steps == 40
        # stepped as chunks in either order, which a threaded session
        # does, the batch meets the same error
        sys = JC.transientsystem(base, h, GaussLegendre(), JC.CPU(), JC.transientfactorization(JC.CPU()))
        for order in ([weak, strong], [strong, weak])
            cerr = try
                JC.gaussbatchintegrate(sys, order, 0.0, 40h, 40, [transientstate(q) for q in order], 1, :ports, 0,
                    1e-9, 1e-10, 15, nothing; chunks = [1:1, 2:2])
            catch e
                e
            end
            @test cerr isa TransientStepError && cerr.step == 2 && cerr.conditions == [findfirst(==(strong), order)]
        end
    end
    # the chunks step no further than the first step one of them failed
    # at, where a failure of theirs would still be the first: the error
    # is the earliest failure whatever order the chunks run in, and on
    # one thread, where they run in turn, a chunk after a failure stops
    # at its step
    let ran = zeros(Int, 3),
        stepper = (k, ch) -> for step in 1:100
            JC.checkchunks(step)
            ran[k] = step
            k == 1 && step == 30 && throw(TransientStepError(30, 30.0, [1], :newton))
            k == 3 && step == 20 && throw(TransientStepError(20, 20.0, [3], :newton))
        end
        err = try
            JC.runchunks(stepper, [1:1, 2:2, 3:3])
        catch e
            e
        end
        @test err isa TransientStepError && err.step == 20 && err.conditions == [3]
        @test Threads.nthreads() > 1 || ran == [30, 30, 20]
    end
    # the responses of the whole batch on one pass equal the members'
    currents = [1e-8*sinpi(2*4.7e9*t) for p in 1:1, t in batch.times]
    tb = transienttangent(batch, currents)
    ab = transientadjoint(batch, weights)
    @test size(tb.outgoing) == (1, 1001, 3) && size(ab.currents) == (1, 1001, 3)
    for j in 1:3
        @test tb.outgoing[:, :, j] ≈ transienttangent(batch[j], currents).outgoing rtol=1e-10
        @test ab.currents[:, :, j] ≈ transientadjoint(batch[j], weights).currents rtol=1e-9
        @test ab.initialflux[:, j] ≈ transientadjoint(batch[j], weights).initialflux rtol=1e-9
    end
    # a tangent measured by a sink stores no history
    measured = zeros(1, 3)
    sink = (k, v, i, o) -> (measured .+= o; nothing)
    ts = transienttangent(batch, currents; outputsink = sink)
    @test isnothing(ts.outgoing) && isnothing(ts.voltage)
    # the stored tangent runs in chunks across the threads of the
    # session and the one with a sink on one task, each condition
    # refreshing its own stage operator
    @test measured ≈ sum(tb.outgoing; dims = 2)[:, 1, :] rtol=1e-10
    # a reuse keeps the workspaces of the tangent and the adjoint for a
    # response of the same shape, which gives the same results on them;
    # another shape builds its own
    reuse = TransientReuse()
    t1 = transienttangent(batch, currents; reuse)
    a1 = transientadjoint(batch, weights; reuse)
    kept = (reuse.tangent, reuse.adjoint)
    t2 = transienttangent(batch, 2 .* currents; reuse)
    a2 = transientadjoint(batch, 2 .* weights; reuse)
    @test reuse.tangent === kept[1] && reuse.adjoint === kept[2]
    @test t1.outgoing ≈ tb.outgoing rtol=1e-12
    @test t2.outgoing ≈ 2 .* t1.outgoing rtol=1e-12
    @test a1.currents ≈ ab.currents rtol=1e-12
    @test a2.currents ≈ 2 .* a1.currents rtol=1e-12
    transienttangent(batch, cat(currents, currents; dims = 3); reuse)
    @test reuse.tangent !== kept[1]
    # a sensitivity keys its workspace on the components, in order,
    # which every call resolves anew from the names: the same list
    # takes the kept workspace, and another list or order builds its
    # own, on both paths to the components
    s1 = transientsensitivity(batch, ["C1", "Lj1"]; reuse)
    c1 = transientadjoint(batch, weights; components = ["C1", "Lj1"], reuse)
    kept = (reuse.tangent, reuse.adjoint)
    s2 = transientsensitivity(batch, [:C1, :Lj1]; reuse)
    c2 = transientadjoint(batch, weights; components = [:C1, :Lj1], reuse)
    @test reuse.tangent === kept[1] && reuse.adjoint === kept[2]
    @test s2.outgoing ≈ s1.outgoing rtol=1e-12
    @test c2.sensitivity ≈ c1.sensitivity rtol=1e-12
    s3 = transientsensitivity(batch, ["Lj1", "C1"]; reuse)
    c3 = transientadjoint(batch, weights; components = ["Lj1", "C1"], reuse)
    @test reuse.tangent !== kept[1] && reuse.adjoint !== kept[2]
    @test s3.outgoing[:, :, [2, 1], :] ≈ s1.outgoing rtol=1e-12
    @test c3.sensitivity[[2, 1], :] ≈ c1.sensitivity rtol=1e-12
    kept = (reuse.tangent, reuse.adjoint)
    transientsensitivity(batch, ["C1"]; reuse)
    transientadjoint(batch, weights; components = ["C1"], reuse)
    @test reuse.tangent !== kept[1] && reuse.adjoint !== kept[2]
    # On the host a batch's stored responses are split into chunks of
    # conditions across the threads of the session, each chunk on its
    # own workspace, and joined; the conditions are independent, so
    # any split gives the same responses. A reuse keeps one workspace
    # per chunk.
    let sys = reuse.system, p = first(batch.problems), nt = length(batch.times),
        (injh, tp) = JC.targetinjection(p, JC.porttargets(p)),
        tc = JC.tangentcurrents(currents, 1, nt, 1), init = JC.tangentinitial(nothing, length(p), 1, 3, 0, 1, 0),
        wh = JC.adjointweights(weights, 1, nt), cp = JC.componentperturbation(p, ["C1"], JC.CPU(); forcing = false)
        split = TransientReuse()
        whole = JC.gaussbatchtangent(batch, tc, injh, tp, init, sys, nothing, nothing, nothing; chunks = [1:3])
        parts = JC.gaussbatchtangent(batch, tc, injh, tp, init, sys, nothing, nothing, split; chunks = [1:1, 2:3])
        @test length(split.tangent) == 2 && split.tangent[1].N == 1 && split.tangent[2].N == 2
        @test parts.outgoing ≈ whole.outgoing rtol=1e-10
        @test parts.finalflux ≈ whole.finalflux rtol=1e-10
        kept = split.tangent
        JC.gaussbatchtangent(batch, tc, injh, tp, init, sys, nothing, nothing, split; chunks = [1:1, 2:3])
        @test split.tangent === kept
        aw = JC.gaussbatchadjoint(batch, wh, true, :outgoing, injh, tp, sys, nothing, nothing, cp, nothing; chunks = [1:3])
        ap = JC.gaussbatchadjoint(batch, wh, true, :outgoing, injh, tp, sys, nothing, nothing, cp, split; chunks = [1:2, 3:3])
        @test length(split.adjoint) == 2
        @test ap.currents ≈ aw.currents rtol=1e-10
        @test ap.sensitivity ≈ aw.sensitivity rtol=1e-10
        @test ap.initialflux ≈ aw.initialflux rtol=1e-10
        # the outputs of a chunked response are the call's result, each
        # chunk writing its own columns, a later call leaves the result
        # of an earlier one as it was, and the workspaces keep no view of
        # either once the call returns
        again = JC.gaussbatchtangent(batch, tc, injh, tp, init, sys, nothing, nothing, split; chunks = [1:1, 2:3])
        @test again.outgoing == parts.outgoing && again.outgoing !== parts.outgoing
        @test all(w -> isempty(parent(w.voltage)) && isempty(parent(w.incident)) && isempty(parent(w.outgoing)), split.tangent)
        @test all(w -> isempty(parent(w.currents)), split.adjoint)
        # the schedule of the noise: the memory's tiles one after another,
        # whole on a device and split across the threads on the host, so
        # that the accumulators live at once are one tile's
        @test JC.conditionschedule(NotTheHost(), 5, 2) == [(1:2, [1:2]), (3:4, [3:4]), (5:5, [5:5])]
        schedule = JC.conditionschedule(JC.CPU(), 5, 2)
        @test [t[1] for t in schedule] == [1:2, 3:4, 5:5]
        @test all(t -> vcat(t[2]...) == t[1] && length(t[2]) <= Threads.nthreads(), schedule)
        # a workspace kept by a reuse holds nothing of a call once the
        # call returns: the adjoint's sink and what it captured are
        # collectible, and the tangent's host currents are dropped,
        # after a return and after an error alike
        # (in a function of its own, so that the frame which held the
        # sink is gone when the collector looks)
        captured = capturingsink(batch, weights, split)
        GC.gc(); GC.gc()
        @test isnothing(captured.value)
        @test all(w -> isnothing(w.ring.sink) && isnothing(w.ring.feedthrough!), split.adjoint)
        @test all(w -> size(w.currents.grid, 2) == 0, split.tangent)
        @test_throws ErrorException transientadjoint(batch, weights; sink = (k, c) -> error("the sink failed"), reuse = split)
        @test all(w -> isnothing(w.ring.sink), split.adjoint)
        # a call which fails in the reset itself releases as one which
        # fails in a step: an initial perturbation of two conditions
        # handed to the workspace of three
        bad = JC.TangentInitial(zeros(length(p), 1, 2), zeros(length(p), 1, 2), zeros(0, 1, 1, 0), zeros(0, 1, 0), true)
        @test_throws BoundsError JC.gaussbatchtangent(batch, tc, injh, tp, bad, sys, sink, nothing, split; chunks = [1:3])
        @test all(w -> size(w.currents.grid, 2) == 0, split.tangent)
        # the entry refuses a line history or a block state given for
        # other than one or every condition before any workspace sees it
        x0, v0 = zeros(length(p), 1), zeros(length(p), 1)
        @test_throws DimensionMismatch JC.tangentinitial((x0, v0, zeros(2, 4, 1, 2)), length(p), 1, 3, 2, 4, 0)
        @test_throws DimensionMismatch JC.tangentinitial((x0, v0, zeros(2, 4, 1, 3), zeros(5, 1, 2)), length(p), 1, 3, 2, 4, 5)
        @test JC.tangentinitial((x0, v0, zeros(2, 4, 1, 3), zeros(5, 1, 3)), length(p), 1, 3, 2, 4, 5).given
        @test transientadjoint(batch, weights; reuse = split).currents ≈ aw.currents rtol=1e-10
        # the workers' reuses are the caller's children on its system,
        # kept across calls and dropped with the system
        rs = JC.workerreuses(split, sys, 3)
        @test rs[1] === split && length(split.children) == 2 && all(c -> c.system === sys, split.children)
        @test JC.workerreuses(split, sys, 3)[2:3] == rs[2:3] && JC.workerreuses(split, sys, 2)[2] === rs[2]
        @test isnothing(JC.workerreuses(nothing, sys, 2)[1].children[1].children)
    end
    # a reuse serves every problem of one compiled circuit driving the
    # same targets: a solve of another drive of the base, under either
    # rule, keeps the system and steps under its own drives, and a
    # problem of another compiled circuit replaces it
    for method in (GaussLegendre(), Trapezoidal())
        reuse = TransientReuse()
        transientsolve(problems[1], (0.0, 1e-9); dt = 5e-12, method, reuse)
        kept = reuse.system
        shared = transientsolve(problems[2], (0.0, 1e-9); dt = 5e-12, method, reuse)
        @test reuse.system === kept
        @test shared.outgoing == transientsolve(problems[2], (0.0, 1e-9); dt = 5e-12, method).outgoing
        @test shared.outgoing != transientsolve(problems[1], (0.0, 1e-9); dt = 5e-12, method).outgoing
        transientsolve(transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)]), (0.0, 1e-9); dt = 5e-12, method, reuse)
        @test reuse.system !== kept
    end
    @test ts.finalflux ≈ tb.finalflux
    # a record of checkpoints: the phases of every window are replayed
    # from its checkpoint, so the responses equal the full record's to
    # the Newton tolerance and nothing per time is stored
    cps = transientsolve(problems, (0.0, 5e-9); dt = 5e-12, record = :checkpoints, checkpointevery = 128)
    @test isnothing(cps.phases) && isnothing(cps.flux)
    @test cps.checkpoints.every == 128 && size(cps.checkpoints.flux) == (2, 8, 3) && size(cps.checkpoints.increment) == (2, 2, 8, 3)
    tc = transienttangent(cps, currents)
    ac = transientadjoint(cps, weights)
    @test tc.outgoing ≈ tb.outgoing rtol=1e-7
    @test ac.currents ≈ ab.currents rtol=1e-7
    @test ac.initialflux ≈ ab.initialflux rtol=1e-7
    @test transientadjoint(cps[2], weights).currents ≈ ab.currents[:, :, 2] rtol=1e-7
    auto = transientsolve(problems[1], (0.0, 5e-9); dt = 5e-12, method = GaussLegendre(), record = :checkpoints)
    @test auto.checkpoints.every == 32 && size(auto.checkpoints.flux) == (2, 32)
    @test transientadjoint(auto, weights).currents ≈ ab.currents[:, :, 1] rtol=1e-7
    @test_throws ArgumentError transientsolve(problems[1], (0.0, 1e-9); dt = 5e-12, method = Trapezoidal(), record = :checkpoints)
    # by default only the ports are recorded
    ports = transientsolve(problems, (0.0, 1e-9); dt = 5e-12)
    @test isnothing(ports.phases) && isnothing(ports.flux) && size(ports.voltage, 3) == 3
    @test_throws ArgumentError transientadjoint(ports, weights[:, 1:201])
    # a common initial state, a state per member, and the reuse
    state = transientstate(base; voltage = [1e-6, 0.0])
    b1 = transientsolve(problems, (0.0, 1e-9); dt = 5e-12, initialstate = state)
    b2 = transientsolve(problems, (0.0, 1e-9); dt = 5e-12, initialstate = fill(state, 3))
    @test b1.finalflux == b2.finalflux
    reuse = TransientReuse()
    b3 = transientsolve(problems, (0.0, 1e-9); dt = 5e-12, reuse)
    # the kept factorizations are one per chunk of conditions, and the
    # chunks together are the batch
    @test reuse.factor isa Vector && all(f -> f isa JC.GaussBatchFactor, reuse.factor)
    @test sum(f -> f.ncolumns, reuse.factor) == 3
    @test length(reuse.factor) == length(JC.batchchunks(JC.CPU(), 3))
    @test transientsolve(problems, (0.0, 1e-9); dt = 5e-12, reuse).finalflux ≈ b3.finalflux rtol=1e-12
    # each condition's state is checked under its own drive, and the
    # sources a batch leaves constant must agree
    r1 = transientproblem(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0))]); sources = [TransientSource(1, 1e-6)])
    r2 = transientproblem(r1; sources = [TransientSource(1, 2e-6)])
    rstates = [transientstate(r1; voltage = [50e-6]), transientstate(r2; voltage = [100e-6])]
    rb = transientsolve([r1, r2], (0.0, 10e-12); dt = 1e-12, initialstate = rstates)
    @test rb.voltage[1, :, 1] ≈ fill(50e-6, 11) rtol=1e-9
    @test rb.voltage[1, :, 2] ≈ fill(100e-6, 11) rtol=1e-9
    @test_throws ArgumentError transientsolve([r1, r2], (0.0, 10e-12); dt = 1e-12, initialstate = rstates[1])
    @test_throws ArgumentError transientsolve([r1, r2], (0.0, 10e-12); dt = 1e-12, initialstate = reverse(rstates))
    # the chunks of a threaded batch: the states are checked before any
    # chunk steps, and an error inside a chunk's steps is thrown as it
    # is rather than as the failure of the task which met it
    sys = JC.transientsystem(r1, 1e-12, GaussLegendre(), JC.CPU(), KLUfactorization())
    chunked(ps, states) = JC.gaussbatchintegrate(sys, ps, 0.0, 10e-12, 10, states, 1, :ports, 0,
        1e-9, 1e-10, 15, nothing; chunks = [1:1, 2:2])
    @test chunked([r1, r2], rstates).voltage ≈ rb.voltage
    @test_throws ArgumentError chunked([r1, r2], reverse(rstates))
    blowup = transientproblem(r1; sources = [TransientSource(1, t -> t > 5e-12 ? NaN : 1e-6)])
    @test_throws ArgumentError chunked([r1, blowup], rstates)
    two = transientproblem(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "0", Capacitor(1e-12)),
        ("I1", "1", "0", CurrentSource(1e-6)), ("I2", "1", "0", CurrentSource(2e-6))]); sources = [TransientSource(:I1, 0.0)])
    @test_throws ArgumentError transientsolve([two, transientproblem(two; sources = [TransientSource(:I2, 0.0)])],
        (0.0, 10e-12); dt = 1e-12)
    same = transientsolve([two, transientproblem(two; sources = [TransientSource(:I1, 1e-6)])], (0.0, 1e-9); dt = 5e-12)
    @test same.voltage[1, end, 1] ≈ -100e-6 rtol=1e-6
    @test same.voltage[1, end, 2] ≈ -150e-6 rtol=1e-6
    # the contracts
    other = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
    @test_throws ArgumentError transientsolve([problems[1], other], (0.0, 1e-9); dt = 5e-12)
    twodrives = transientproblem(base; sources = [TransientSource(1, t -> 0.0), TransientSource(1, t -> 0.0)])
    @test_throws ArgumentError transientsolve([problems[1], twodrives], (0.0, 1e-9); dt = 5e-12)
    @test_throws ArgumentError transientsolve(problems, (0.0, 1e-9); dt = 5e-12, method = Trapezoidal())
    @test_throws ArgumentError transientproblem(base; sources = [TransientSource(3, 1e-6)])
    @test_throws BoundsError batch[4]
end

@testset "the rate along the algebraic constraints and its tangent" begin
    # an inductive divider: node 2 has neither capacitance nor a
    # source, so the rate along it is read from the differentiated
    # constraint, whose row a perturbed inductor moves with it; the
    # whole final rate against a central difference, under either
    # rule, from a record of the states, from checkpoints and in a
    # batch
    dividernet = [("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "0", Capacitor(1e-12)),
        ("L1", "1", "2", Inductor(1e-9)), ("L2", "2", "0", Inductor(1e-9))]
    divider = Circuit(dividernet)
    sources = [TransientSource(1, t -> 1e-6*sinpi(2*2e9*t))]
    prob = transientproblem(divider; sources)
    scaled(r) = transientproblem(Circuit(scaledentries(dividernet, "L1", r)); sources)
    tspan, dt = (0.0, 153e-12), 1e-12
    eps = 1e-4
    for method in (GaussLegendre(), Trapezoidal())
        fd = (transientsolve(scaled(1 + eps), tspan; dt, method).finalrate .-
            transientsolve(scaled(1 - eps), tspan; dt, method).finalrate) ./ (2eps)
        sens = transientsensitivity(transientsolve(prob, tspan; dt, method, record = :states), ["L1"])
        @test sens.finalrate[:, 1] ≈ fd rtol=1e-6
        method isa GaussLegendre || continue
        sens = transientsensitivity(transientsolve(prob, tspan; dt, record = :checkpoints), ["L1"])
        @test sens.finalrate[:, 1] ≈ fd rtol=1e-6
        sb = transientsensitivity(transientsolve([prob, prob], tspan; dt, record = :states), ["L1"])
        @test sb.finalrate[:, 1, 1] ≈ fd rtol=1e-6
        @test sb.finalrate[:, 1, 2] ≈ fd rtol=1e-6
    end
    # a tangent current into a direction the solve leaves undriven: an
    # unterminated port on an inductor, at rest, with a quadratic
    # current whose rate at the end the reading carries, against a
    # central difference and the analytic rate
    L, T, Iprobe = 1e-9, 100e-12, 1e-8
    rest = Circuit([(:p, 1, 0, Port(1; termination = nothing)), (:l, 1, 0, Inductor(L))])
    probe(t) = Iprobe*(t/T)^2
    shifted(a) = transientproblem(rest; sources = [TransientSource(1, t -> a*probe(t))])
    for method in (GaussLegendre(), Trapezoidal()), record in (:states, :checkpoints)
        record == :checkpoints && !(method isa GaussLegendre) && continue
        s = transientsolve(transientproblem(rest), (0.0, T); dt = 1e-12, record, method)
        tg = transienttangent(s, reshape(probe.(s.times), 1, :))
        fd = (transientsolve(shifted(1e-3), (0.0, T); dt = 1e-12, method).finalrate .-
            transientsolve(shifted(-1e-3), (0.0, T); dt = 1e-12, method).finalrate) ./ 2e-3
        @test tg.finalrate ≈ fd rtol=1e-9
        @test tg.finalrate[1] ≈ 2L*Iprobe/(T*JC.phi0) rtol=1e-9
    end
    # a constraint a junction touches, and one a source drives: the
    # reported rate satisfies the differentiated constraint to
    # roundoff under either rule, where the rules' own rates do not
    ratepart(p, sys, x, v, t) = begin
        Lm, RJ = JC.hostsparse(sys.L), JC.hostsparse(sys.RJ)
        lmolj, hr = Array(sys.lmolj), JC.hostrelations(sys.relations)
        delta = 1e-3*sys.h
        bdot = (JC.hostdrivecurrent(sys, t + delta, p) .- JC.hostdrivecurrent(sys, t - delta, p)) ./ (2delta)
        lv = Lm*v .+ transpose(RJ)*(lmolj .* JC.derivativeat(hr, RJ*x) .* (RJ*v)) .- bdot
        norm(p.constraints*lv, Inf)/max(norm(bdot, Inf), norm(Lm*v, Inf), 1.0)
    end
    Lc, C = 1e-9, 1e-12
    junction = transientproblem(Circuit([("C1", "1", "2", Capacitor(C)), ("L1", "1", "0", Inductor(Lc)), ("L2", "2", "0", Inductor(Lc)), ("Lj1", "1", "0", JosephsonJunction(Lc))]))
    wj = sqrt(3/(Lc*C))
    jstate = transientstate(junction; voltage = [1e-6, -2e-6])
    wd = 0.37*sqrt(2/(Lc*C))
    driven = transientproblem(Circuit([("C1", "1", "2", Capacitor(C)), ("L1", "1", "0", Inductor(Lc)), ("L2", "2", "0", Inductor(Lc)), ("I1", "1", "0", CurrentSource(0.0))]);
        sources = [TransientSource(:I1, t -> 1e-6*(1 - cos(wd*t)))])
    for (p, w, state) in ((junction, wj, jstate), (driven, wd, transientstate(driven))), method in (GaussLegendre(), Trapezoidal())
        @test length(p.algebraic) == 1
        period = 2pi/w
        sol = transientsolve(p, (0.0, 4period); dt = period/800, method, initialstate = state, record = :states)
        sys = JC.transientsystem(p, sol.dt, method, JC.CPU(), KLUfactorization())
        # to the accuracy of this check's own difference of the drive
        @test ratepart(p, sys, sol.finalflux, sol.finalrate, 4period) < 1e-10
        @test ratepart(p, sys, sol.flux[:, end - 1], sol.rate[:, end - 1], sol.times[end - 1]) < 1e-10
    end
    # the tangent's reading on a constraint a junction and a source
    # touch, from rest so that the start stays consistent whatever the
    # inductances: along the inductances, with the junction's stiffness
    # and its curvature in the constraint's row, and along the source's
    # own current, whose rate at the end the reading carries, under
    # either rule. The stage iteration of the Gauss-Legendre tangent
    # stops at a roundoff floor that allows for the cancellation of
    # `C d` along the direction; on the island's two nodes the product
    # itself is a thousandth of its rounding.
    bothnet = [("C1", "1", "2", Capacitor(C)), ("L1", "1", "0", Inductor(Lc)), ("L2", "2", "0", Inductor(Lc)), ("Lj1", "1", "0", JosephsonJunction(Lc)), ("I1", "1", "0", CurrentSource(0.0))]
    both = Circuit(bothnet)
    wb = 0.37*sqrt(3/(Lc*C))
    bdrive = t -> 1e-7*(1 - cos(wb*t))
    bspan, bdt = (0.0, 8pi/wb), 2pi/wb/800
    bdI = t -> 1e-8*(1 - cos(2wb*t))
    for method in (GaussLegendre(), Trapezoidal())
        bsol = transientsolve(transientproblem(both; sources = [TransientSource(:I1, bdrive)]), bspan; dt = bdt,
            method, record = :states)
        bsens = transientsensitivity(bsol, ["Lj1", "L1"])
        for (k, name) in enumerate(["Lj1", "L1"])
            fd = (transientsolve(transientproblem(Circuit(scaledentries(bothnet, name, 1 + 1e-4));
                    sources = [TransientSource(:I1, bdrive)]), bspan; dt = bdt, method).finalrate .-
                transientsolve(transientproblem(Circuit(scaledentries(bothnet, name, 1 - 1e-4));
                    sources = [TransientSource(:I1, bdrive)]), bspan; dt = bdt, method).finalrate) ./ 2e-4
            @test bsens.finalrate[:, k] ≈ fd rtol=1e-5
        end
        btg = transienttangent(bsol, reshape(bdI.(bsol.times), 1, :); targets = [:I1])
        bfd = (transientsolve(transientproblem(both; sources = [TransientSource(:I1, t -> bdrive(t) + 1e-3*bdI(t))]), bspan;
                dt = bdt, method).finalrate .-
            transientsolve(transientproblem(both; sources = [TransientSource(:I1, t -> bdrive(t) - 1e-3*bdI(t))]), bspan;
                dt = bdt, method).finalrate) ./ 2e-3
        @test btg.finalrate ≈ bfd rtol=1e-4
        # every recorded state satisfies the constraint to the roundoff
        # of the terms the projection balances, under either rule, and
        # its differentiated form to the rounding of the difference
        # which reads the drive's rate; so a restart from the end is
        # accepted
        bsys = JC.transientsystem(bsol.problem, bdt, method, JC.CPU(), KLUfactorization())
        @test maximum(JC.transientconsistency(bsys, bsol.flux[:, k], bsol.rate[:, k], bsol.times[k]).violation
            for k in eachindex(bsol.times)) < 1e-10
        @test transientsolve(bsol.problem, (bspan[2], bspan[2] + 8bdt); dt = bdt, method,
            initialstate = transientstate(bsol)).stats.steps == 8
    end
    # two junctions in series on a capacitor free node balance their
    # currents in the constraint, so the residual's floor comes from
    # the sizes of the terms and not from their sum, under either rule
    series = Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "2", Capacitor(100e-15)),
        ("jj", "2", "3", JosephsonJunction(1e-9)), ("lj", "3", "0", NonlinearInductor(1e-9, PolynomialCPR([1.0, 0.25, -1/6]))),
        ("c2", "2", "0", Capacitor(1e-12))])
    sprob = transientproblem(series; sources = [TransientSource(1, t -> 0.3e-6*sinpi(2*3e9*t))])
    for method in (GaussLegendre(), Trapezoidal())
        ssol = transientsolve(sprob, (0.0, 0.1e-9); dt = 2e-12, method, record = :states)
        ssys = JC.transientsystem(sprob, 2e-12, method, JC.CPU(), KLUfactorization())
        @test maximum(JC.transientconsistency(ssys, ssol.flux[:, k], ssol.rate[:, k], ssol.times[k]).violation
            for k in eachindex(ssol.times)) < 1e-10
    end
    # a junction's phase is a difference of node fluxes which may be
    # far larger than it, whose rounding the floor allows for: a
    # series array of identical junctions without capacitance on its
    # inner nodes, which divide the flux evenly, and a pair of unequal
    # junctions driven into the voltage state, whose node fluxes wind
    # up while their currents stay balanced at the node between them
    M = 24
    array = transientproblem(Circuit(vcat([("p1", "1", "0", Port(1)), ("c1", "1", "0", Capacitor(1e-12))],
            [("J$j", string(j), j == M ? "0" : string(j + 1), JosephsonJunction(1e-9/M)) for j in 1:M]));
        sources = [TransientSource(1, t -> 0.3e-6*sinpi(2*3e9*t))])
    order = [findfirst(==(string(j)), array.circuit.nodenames) - 1 for j in 1:M]
    Ic = JC.phi0/1e-9
    pair = transientproblem(Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "0", Capacitor(1e-12)),
            ("J1", "1", "2", JosephsonJunction(1e-9)), ("J2", "2", "0", JosephsonJunction(1.3e-9))]);
        sources = [TransientSource(1, t -> 2.5Ic*(t <= 0 ? 0.0 : t >= 0.2e-9 ? 1.0 : (1 - cospi(t/0.2e-9))/2))])
    for method in (GaussLegendre(), Trapezoidal())
        asol = transientsolve(array, (0.0, 0.2e-9); dt = 1e-12, method)
        @test asol.finalflux[order] ≈ asol.finalflux[order[1]] .* (M:-1:1) ./ M rtol=1e-10
        psol = transientsolve(pair, (0.0, 0.5e-9); dt = 1e-12, method)
        x1, x2 = psol.finalflux[1], psol.finalflux[2]
        @test x2 > 10
        @test sin(x1 - x2) ≈ sin(x2)/1.3 rtol=1e-10
    end
    # a tangent current on a record of two or three samples: its rate
    # at each time from the line or the quadratic through them, under
    # every rule, against the exact voltage `L dI/dt` of a ramp across
    # an unterminated port on an inductor, and the adjoint's
    # contraction against the same reference
    ramp = Circuit([(:p, 1, 0, Port(1; termination = nothing)), (:l, 1, 0, Inductor(1e-9))])
    rprob = transientproblem(ramp)
    T, Iprobe = 100e-12, 1e-8
    for method in (GaussLegendre(), Trapezoidal(), BackwardEuler()), steps in 1:3
        rsol = transientsolve(rprob, (0.0, T); dt = T/steps, method, record = :states)
        di = reshape(Iprobe .* rsol.times ./ T, 1, :)
        initial = ([0.0], [1e-9*Iprobe/(T*JC.phi0)])
        rtg = transienttangent(rsol, di; initialstate = initial)
        @test rtg.voltage ≈ fill(1e-9*Iprobe/T, 1, steps + 1) rtol=1e-8
        @test JC.phi0 .* rtg.finalrate ≈ [1e-9*Iprobe/T] rtol=1e-8
        rw = reshape([1.0 + k for k in 0:steps], 1, :)
        rad = transientadjoint(rsol, rw; quantity = :voltage)
        @test dot(vec(rad.currents), vec(di)) + dot(rad.initialflux, initial[1]) + dot(rad.initialrate, initial[2]) ≈
            sum(rw)*1e-9*Iprobe/T rtol=1e-8
    end
    # a subcircuit the stiffness couples to nothing touched keeps its
    # invariant reading: coupled pairs beside a junction on a winding
    # add invariant directions and leave the projected ones alone
    core = [("Lj1", "1", "0", JosephsonJunction(1e-9)), ("L1", "1", "0", Inductor(1e-9)), ("L2", "2", "0", Inductor(1e-9)), ("K1", "L1", "L2", MutualInductor(0.3)),
        ("C2", "2", "0", Capacitor(1e-12))]
    pairs(m) = Circuit(vcat(core, [c for j in 1:m for c in (("La$j", "$(2j + 1)", "0", Inductor(1e-9)), ("Ca$j", "$(2j + 1)", "0", Capacitor(1e-12)),
        ("Lb$j", "$(2j + 2)", "0", Inductor(1e-9)), ("Cb$j", "$(2j + 2)", "0", Capacitor(1e-12)), ("Kab$j", "La$j", "Lb$j", MutualInductor(0.3)))]))
    for m in (0, 4)
        csys = JC.transientsystem(transientproblem(pairs(m)), 1e-12, GaussLegendre(), JC.CPU(), KLUfactorization())
        @test length(csys.projection.directions) == 3
        @test (isnothing(csys.invariant) ? 0 : size(csys.invariant.Z, 2)) == 2m
    end
    # the bases of the directions and of the constraints are chosen
    # apart, so a coupling's sign can differ from its transpose's: a
    # resistor dangling from the junction's node changes the bases and
    # nothing else, and the partition, the port voltage under every
    # rule and the constraints at every state are those of the
    # circuit without it; nor does the partition depend on the sign of
    # any constraint's row
    hanging(dangling) = Circuit(vcat([(:jj, 1, 0, JosephsonJunction(1e-9)), (:l1, 1, 0, Inductor(1e-9)),
        (:l3, 3, 0, Inductor(1e-9)), (:k, :l1, :l3, MutualInductor(0.3)), (:c3, 3, 0, Capacitor(1e-12)),
        (:p, 1, 0, Port(1; termination = nothing))], dangling ? [(:r12, 1, 2, Resistor(50.0))] : []))
    hsources = [TransientSource(1, t -> 0.2e-6*(1 - cos(2pi*1e9*t)))]
    hprob = transientproblem(hanging(true); sources = hsources)
    hsys = JC.transientsystem(hprob, 50e-12, GaussLegendre(), JC.CPU(), KLUfactorization())
    @test length(hsys.projection.directions) == 3 && isnothing(hsys.invariant)
    for method in (GaussLegendre(), Trapezoidal(), BackwardEuler())
        with = transientsolve(hprob, (0.0, 1e-9); dt = 50e-12, method, rtol = 1e-12, record = :states)
        without = transientsolve(transientproblem(hanging(false); sources = hsources), (0.0, 1e-9);
            dt = 50e-12, method, rtol = 1e-12)
        @test with.voltage ≈ without.voltage rtol=1e-8
        msys = JC.transientsystem(hprob, 50e-12, method, JC.CPU(), KLUfactorization())
        @test maximum(JC.transientconsistency(msys, with.flux[:, k], with.rate[:, k], with.times[k]).violation
            for k in eachindex(with.times)) < 1e-10
    end
    hmats = (JC.hostsparse(hsys.L), JC.hostsparse(hsys.RJ), JC.hostsparse(hsys.injection),
        JC.hostsparse(hsys.lineinjection), JC.hostsparse(hsys.blockscatter))
    hZ, hZt = droptol!(sparse(hprob.directions), 1e-14), droptol!(sparse(hprob.constraints), 1e-14)
    htouched = JC.touchedalgebraic(hZ, hZt, hmats..., JC.statefulcolumns(hprob, hmats[5]))[1]
    @test length(htouched) == 3
    for i in axes(hZt, 1)
        flipped = copy(hZt)
        flipped[i, :] .*= -1
        @test JC.touchedalgebraic(hZ, flipped, hmats..., JC.statefulcolumns(hprob, hmats[5]))[1] == htouched
    end
    # a port on an algebraic direction reads the rate at every time, so
    # its voltage's sensitivity to the junction and the inductor on the
    # direction goes through the reading at every time, under either
    # rule, from the record and from checkpoints, and the adjoint
    # through the reading's transpose
    island(rj, rl) = Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "0", Capacitor(1e-12)),
        ("jj", "1", "2", JosephsonJunction(rj*1e-9)), ("l2", "2", "0", Inductor(rl*2e-9)),
        ("p2", "2", "0", Port(2; termination = nothing))])
    isources = [TransientSource(1, t -> 0.3e-6*sinpi(2*3e9*t)), TransientSource(2, t -> 0.0)]
    ispan, idt = (0.0, 0.5e-9), 2e-12
    iweights = [cospi(2*1.7e9*t + q) for q in 1:2, t in 0:idt:ispan[2]]
    for method in (GaussLegendre(), Trapezoidal())
        isol = transientsolve(transientproblem(island(1.0, 1.0); sources = isources), ispan; dt = idt, method,
            record = :states, rtol = 1e-12)
        isens = transientsensitivity(isol, ["jj", "l2"])
        for (k, scaled) in enumerate((r -> island(r, 1.0), r -> island(1.0, r)))
            fd = (transientsolve(transientproblem(scaled(1 + 1e-4); sources = isources), ispan; dt = idt, method, rtol = 1e-12).voltage .-
                transientsolve(transientproblem(scaled(1 - 1e-4); sources = isources), ispan; dt = idt, method, rtol = 1e-12).voltage) ./ 2e-4
            @test isens.voltage[:, :, k] ≈ fd rtol=1e-5 atol=1e-7*maximum(abs, fd)
        end
        iad = transientadjoint(isol, iweights; quantity = :voltage, components = ["jj", "l2"])
        @test iad.sensitivity ≈ [sum(iweights .* isens.voltage[:, :, k]) for k in 1:2] rtol=1e-8
        method isa GaussLegendre || continue
        icps = transientsolve(transientproblem(island(1.0, 1.0); sources = isources), ispan; dt = idt, method,
            record = :checkpoints, checkpointevery = 16, rtol = 1e-12)
        @test transientsensitivity(icps, ["jj", "l2"]).voltage ≈ isens.voltage rtol=1e-7
        @test transientadjoint(icps, iweights; quantity = :voltage, components = ["jj", "l2"]).sensitivity ≈ iad.sensitivity rtol=1e-7
    end
end

@testset "the stationary limit agrees with harmonic balance" begin
    circuit = Circuit(vcat(rc, [("Lj1", "1", "0", JosephsonJunction(1e-9))]))
    f, ip = 3e9, 0.12e-6
    prob = transientproblem(circuit; sources = [TransientSource(1, t -> ip*cospi(2f*t))])
    sol = transientsolve(prob, (0.0, 12e-9); dt = 0.5e-12, method = Trapezoidal())
    # a factorization is kept while its corrections contract the
    # residual
    @test sol.stats.factorizations < 10
    # the iterative step, on the package's GMRES with the factorization
    # as its preconditioner, gives the same trajectory from one
    # factorization, to the Newton tolerance accumulated over the steps
    iterative = transientsolve(prob, (0.0, 12e-9); dt = 0.5e-12, method = Trapezoidal(), linearsolver = GMRES())
    @test iterative.stats.factorizations == 1
    @test iterative.stats.kryloviterations > 0
    @test iterative.voltage ≈ sol.voltage rtol=1e-4 atol=1e-12
    # a pumped junction whose steps take a second correction: the
    # iterative step converges it as the direct one does
    jpa = Circuit([("p1","1","0",Port(1)), ("c1","1","2",Capacitor(100e-15)),
        ("lj","2","0",JosephsonJunction(1e-9)), ("c2","2","0",Capacitor(1e-12))])
    ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
    strong = transientproblem(jpa; sources = [TransientSource(1,
        t -> 2*4*0.00565e-6*ramp(t)*cospi(2*4.75e9*t))])
    direct = transientsolve(strong, (0.0, 4e-9); dt = 2e-12, method = Trapezoidal())
    @test direct.stats.newtoncorrections > direct.stats.steps
    krylov = transientsolve(strong, (0.0, 4e-9); dt = 2e-12, method = Trapezoidal(), linearsolver = GMRES())
    @test krylov.stats.newtoncorrections > krylov.stats.steps
    @test krylov.voltage ≈ direct.voltage rtol=1e-4 atol=1e-12
    # the Krylov budget is the solver's own: one iteration a solve
    # cannot follow the pump on the kept factorization at every step,
    # so the step refreshes it, and reaches the same trajectory; the
    # kept workspace has the solver's restart length
    budget = transientsolve(strong, (0.0, 4e-9); dt = 2e-12, method = Trapezoidal(),
        linearsolver = GMRES(restart = 1, maxrestarts = 1))
    @test budget.stats.factorizations > krylov.stats.factorizations
    @test budget.voltage ≈ direct.voltage rtol=1e-4 atol=1e-12
    @test_throws ArgumentError transientsolve(prob, (0.0, 1e-9); dt = 1e-12, linearsolver = :gmres)
    # a reuse object carries the system, its factorization and the Krylov
    # workspace to the next solve, tangent and adjoint of the problem
    reuse = TransientReuse()
    a = transientsolve(prob, (0.0, 2e-9); dt = 2e-12, record = :phases, reuse)
    b = transientsolve(prob, (0.0, 2e-9); dt = 2e-12, record = :phases, reuse)
    @test a.voltage == b.voltage
    @test reuse.system === JC.transientsystem(reuse, prob, 2e-12, GaussLegendre(), JC.CPU(), KLUfactorization())
    currents = [1e-8*sinpi(2*1.7e9*t) for _ in 1:1, t in a.times]
    weights = ones(1, length(a.times))
    @test transienttangent(a, currents; reuse).voltage ≈ transienttangent(a, currents).voltage
    @test transientadjoint(a, weights; reuse).currents ≈ transientadjoint(a, weights).currents
    g = transientsolve(prob, (0.0, 2e-9); dt = 2e-12, method = Trapezoidal(), linearsolver = GMRES(restart = 7), reuse)
    @test reuse.workspace isa JC.GMRESWorkspace && size(reuse.workspace.H, 2) == 7
    # a different step replaces the kept system
    transientsolve(prob, (0.0, 2e-9); dt = 1e-12, reuse)
    @test reuse.system.h == 1e-12
    @test_throws ArgumentError transientsolve(prob, (0.0, 1e-9); dt = 1e-12, reuse = 3)
    # harmonic balance stores the conjugacy representatives: a physical
    # cosine of peak Ip has the positive frequency coefficient Ip/2
    hb = hbnlsolve((2pi*f,), (9,), [(mode = (1,), port = 1, current = ip/2)],
        circuit, Dict{Symbol,Float64}(); method = Newton(), keyedarrays = false)
    window(t) = 2e-9 <= t <= 12e-9 ? sinpi((t - 2e-9)/10e-9)^2 : 0.0
    for (k, m) in enumerate(hb.frequencies.modes)
        m[1] > 5 && continue
        measured = transientdemodulate(sol, 1, 2pi*m[1]*f; quantity = :voltage, window)
        expected = 2im*2pi*f*m[1]*JC.phi0*hb.nodeflux[k]
        @test measured ≈ expected rtol=0.003 atol=1e-13
    end
end
