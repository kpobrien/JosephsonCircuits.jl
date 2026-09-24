# The solver paths which divide work between threads, run by compare.jl in
# a process of a stated thread count. The results are written to the path
# given as the first argument, and compare.jl checks that the run with one
# thread and the run with four wrote the same numbers.
using JosephsonCircuits
using LinearAlgebra
using Serialization
using Test

const JC = JosephsonCircuits
BLAS.set_num_threads(1)

@test Threads.nthreads() == parse(Int, ARGS[2])
# the threads are really there: a static loop over one index per thread
# runs on that many distinct threads
ids = zeros(Int, Threads.nthreads())
Threads.@threads :static for i in eachindex(ids)
    ids[i] = Threads.threadid()
end
@test length(unique(ids)) == Threads.nthreads()

@testset "the threaded solver paths" begin
    c = compile(Circuit([(:p1, 1, 0, Port(1)),
        (:cc, 1, 2, Capacitor(100e-15)),
        (:jj, 2, 0, JosephsonJunction(:Lj)),
        (:cj, 2, 0, Capacitor(1e-12)),
        (:loss, 2, 0, Resistor(2000.0)),
        (:p2, 2, 0, Port(2))]))
    defs = Dict(:Lj => 1e-9)
    ws = 2*pi*collect(range(4.3e9, 4.9e9; length = 7))

    # the linearized sweep divides its frequencies between batches: every
    # output is the output of the sweep run in one batch, including the
    # noise of the lossy resistor and the sensitivity of the junction
    batched = hblinsolve(ws, c, defs; nbatches = Threads.nthreads(),
        keyedarrays = false, sensitivitynames = [:jj],
        returnSsensitivity = true, returnSnoise = true)
    serial = hblinsolve(ws, c, defs; nbatches = 1, keyedarrays = false,
        sensitivitynames = [:jj], returnSsensitivity = true,
        returnSnoise = true)
    @test batched.S ≈ serial.S
    @test batched.Ssensitivity ≈ serial.Ssensitivity
    @test batched.Snoise ≈ serial.Snoise
    @test !isempty(batched.Snoise) && norm(batched.Snoise) > 0

    # caches over one compiled circuit solve at the same time: the
    # compilation is shared and every other piece of state is the cache's
    src = [(mode = (1,), port = 1, current = 2e-9)]
    caches = [hbcache((2*pi*4.75e9,), (4,), src, c, defs) for _ in 1:7]
    @test all(cache -> cache.compiled === c, caches)
    fluxes = Vector{Any}(undef, length(caches))
    @sync for j in eachindex(caches)
        Threads.@spawn fluxes[j] =
            copy(hbsolve!(caches[j], (Lj = (1 + 0.01*j)*1e-9,)).nodeflux)
    end
    @test all(cache -> cache.converged, caches)
    @test caches[1].reuse.sys !== caches[2].reuse.sys
    @test caches[1].nm.Cnm !== caches[2].nm.Cnm

    # sensitivities at the same time from one pump solution: what
    # differentiates its operating point evaluates on a system of its own,
    # so the calls neither see each other nor move the point
    nl = hbnlsolve((2*pi*4.75e9,), (4,), src, c, defs;
        returnoperatingpoint = true, keyedarrays = false)
    names = ["jj", "cc", "cj"]
    nm = JC.numericmatrices(c, defs; Nmodes = length(nl.modes))
    idx = [JC.componentindex(c, n) for n in names]
    function derivatives()
        dFr = JC.calcresidualsensitivity(nl.operatingpoint, c, nm, idx)
        S = hblinsolve(ws, c, defs; Nmodulationharmonics = (2,),
            nonlinear = nl, nbatches = 1, keyedarrays = false,
            sensitivitynames = names, returnSsensitivity = true,
            sensitivityresidual = dFr,
            sensitivitymode = :forward).Ssensitivity
        return Matrix(dFr), S
    end
    alone = derivatives()
    together = [fetch(t) for _ in 1:4
        for t in [Threads.@spawn derivatives() for _ in 1:8]]
    @test all(d -> isapprox(d[1], alone[1]; rtol = 1e-12) &&
        isapprox(d[2], alone[2]; rtol = 1e-12), together)

    # the same with a direct current block, a scattering block and a
    # polynomial current-phase relation, whose derivatives also read the
    # block rows of the canonical work and the cached second derivative of
    # the relation, with the forward and the reverse contraction
    # interleaved; no array of the operating point moves
    Z0, R = 50.0, 5.0
    series(w) = (z = R/Z0; [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
    dseries(w) = (z = R/Z0; d = 2/(Z0*(z+2)^2); [d -d; -d d])
    cb = compile(Circuit([:p1 => Port(1; termination = nothing),
            :r1 => Resistor(50.0),
            :cc => ScatteringParameters(series; nports = 2,
                grounded = false, derivatives = (R = dseries,)),
            :jj => NonlinearInductor(1e-9,
                PolynomialCPR([1.0, 0.25, -1/6, 0.0, 1/120])),
            :c2 => Capacitor(1000e-15)],
        [((:p1, 1), (:r1, 1), (:cc, 1, 1)), ((:cc, 2, 1), (:jj, 1), (:c2, 1)),
            ((:cc, 1, 2), (:cc, 2, 2), (:jj, 2), (:c2, 2), (:r1, 2), (:p1, 2),
                Ground)]))
    nlb = hbnlsolve((2*pi*4.75e9,), (4,),
        [(mode = (1,), port = 1, current = 0.00565e-6),
            (mode = (0,), port = 1, current = 1e-8)], cb, Dict{Any,Any}();
        dc = true, odd = true, even = true, method = Newton(),
        returnoperatingpoint = true, keyedarrays = false)
    opb = nlb.operatingpoint
    nmb = JC.numericmatrices(cb, Dict{Any,Any}(); Nmodes = opb.Nmodes)
    namesb = ["jj", "r1", "c2"]
    idxb = [JC.componentindex(cb, n) for n in namesb]
    pairs = JC.designblockjacobian(cb, [:R])
    sig = JC.truncfreqs(JC.calcfreqsdft((2,)); dc = true, odd = true,
        even = true)
    arrays(x) = [copy(getfield(x, f)) for f in fieldnames(typeof(x))
        if getfield(x, f) isa AbstractArray]
    snapshot() = map(arrays, (opb, opb.sys, opb.dc, opb.dc.work))
    held = snapshot()
    function blockderivatives(mode)
        dFr = JC.calcresidualsensitivity(opb, cb, nmb, idxb)
        S = hblinsolve(ws[1:2], cb, Dict{Any,Any}(), sig; nonlinear = nlb,
            nbatches = 2, keyedarrays = false, sensitivitynames = namesb,
            returnSsensitivity = true, sensitivityresidual = dFr,
            sensitivitymode = mode).Ssensitivity
        return Matrix(dFr),
            Matrix(JC.calcblockresidualsensitivity(opb, cb, pairs)), S
    end
    fwd, rev = blockderivatives(:forward), blockderivatives(:reverse)
    @test norm(fwd[2]) > 0 && norm(fwd[3]) > 0
    @test isapprox(fwd[3], rev[3]; rtol = 1e-8)
    interleaved = [fetch(t) for t in [Threads.@spawn blockderivatives(
        isodd(k) ? :forward : :reverse) for k in 1:8]]
    @test all(k -> all(isapprox(a, b; rtol = 1e-12) for (a, b) in
        zip(interleaved[k], isodd(k) ? fwd : rev)), 1:8)
    @test isequal(snapshot(), held)

    # the transient batches its problems over the threads, and the
    # tangent, the adjoint, the noise and the gain share one reuse object
    # between them
    nsteps, T, fp = 128, 0.5e-9, 4.75e9
    base = transientproblem(c, defs; sources = [TransientSource(1, t -> 0.0)])
    problems = [transientproblem(base; sources = [TransientSource(1,
        let ip = j*1e-10; t -> t <= 0 ? 0.0 : 2*ip*sinpi(2*fp*t) end)])
        for j in 1:7]
    chunks = JC.batchchunks(JC.CPU(), length(problems))
    @test length(chunks) == min(Threads.nthreads(), length(problems))
    @test vcat(chunks...) == collect(eachindex(problems))
    sol = transientsolve(problems, (0.0, T*(nsteps-1)/nsteps);
        dt = T/nsteps, method = GaussLegendre(), record = :checkpoints,
        checkpointevery = 16)
    currents = [sinpi(2*4.6e9*t) for _ in 1:2, t in sol.times]
    weights = [cospi(2*4.4e9*t) for _ in 1:2, t in sol.times]
    reuse = TransientReuse()
    tangent = transienttangent(sol, currents; reuse)
    adjoint = transientadjoint(sol, weights; reuse)
    plan = transientquantumplan(sol, sol.times, [2/T]; ports = [2])
    noise = transientnoise(sol, plan; frequencies = [1/T, 2/T, 3/T],
        weights = fill(1/T, 3), inputs = plan, reuse)
    gain = transientgain(sol, plan, plan; reuse)
    # one workspace per chunk, so a threaded run keeps more than one
    @test length([reuse; reuse.children]) >= (Threads.nthreads() > 1 ? 2 : 1)
    @test transientnoise(sol, plan; frequencies = [1/T, 2/T, 3/T],
        weights = fill(1/T, 3), inputs = plan, reuse).covariance ≈
        noise.covariance
    # a batched problem's sensitivities are the single problem's
    for j in (1, 4, 7)
        @test tangent.outgoing[:, :, j] ≈
            transienttangent(sol[j], currents).outgoing rtol = 1e-8
        @test adjoint.currents[:, :, j] ≈
            transientadjoint(sol[j], weights).currents rtol = 1e-8
    end

    # a pump solve whose nonlinear term maps are longer than `hostlooplimit`,
    # which a threaded process applies as KernelAbstractions kernels on the
    # host rather than as plain loops
    chain = Any[(:p1, 1, 0, Port(1))]
    for i in 1:4096
        push!(chain, (Symbol(:lj, i), i, i + 1, JosephsonJunction(100e-12)),
            (Symbol(:c, i), i, 0, Capacitor(40e-15)))
    end
    push!(chain, (:r2, 4097, 0, Resistor(50.0)))
    long = hbnlsolve((2*pi*7e9,), (16,),
        [(mode = (1,), port = 1, current = 1e-7)], Circuit(chain);
        keyedarrays = false)
    @test long.solverinfo.converged

    serialize(ARGS[1], (; S = batched.S, Ssensitivity = batched.Ssensitivity,
        Snoise = batched.Snoise, fluxes, finalflux = sol.finalflux,
        tangent = tangent.outgoing, adjoint = adjoint.currents,
        covariance = noise.covariance, gain, longflux = long.nodeflux))
end
