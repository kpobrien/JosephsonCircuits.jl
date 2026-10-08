using JosephsonCircuits, LinearAlgebra, SparseArrays, Random, Test, Logging
using JosephsonCircuits: BlockDiagonal, Floquet, HarmonicBand

isdefined(Main, :testchaincircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))

# The Newton-Krylov engine over that solver: the forcing terms it picks,
# the parameters it validates, and the reason every solve ends with.
@testset "nlsolvekrylov! forcing terms and parameter validation" begin
    # a small nonlinear system preconditioned by its own exact Jacobian,
    # wrapped in the minimal AbstractPreconditioner interface
    mutable struct ExactP <: JosephsonCircuits.AbstractPreconditioner
        J::Matrix{Float64}
        F::Any
    end
    JosephsonCircuits.updatepreconditioner!(pc::ExactP, x::AbstractVector) =
        (pc.J .= [2x[1] 0.0; 0.0 3x[2]^2]; pc.F = lu(pc.J); pc)
    JosephsonCircuits.applypreconditioner!(z::AbstractVector, pc::ExactP,
        r::AbstractVector) = ldiv!(z, pc.F, r)

    xpt = [1.0, 1.0]
    fj!(F, J, x) = begin
        copyto!(xpt, x)
        isnothing(F) || (F .= [x[1]^2 - 2, x[2]^3 - 3])
        nothing
    end
    jvp!(y, v) = (y .= [2*xpt[1]*v[1], 3*xpt[2]^2*v[2]])

    # every refresh policy and linear solver setting reaches the root
    for (refresh, ls) in ((Always(), GMRES()), (Probe(), GMRES()),
            (Never(), GMRES(restart = 5, maxrestarts = 2)))
        x = [1.0, 1.0]
        F = zeros(2)
        info = JosephsonCircuits.nlsolvekrylov!(fj!, jvp!, F, x,
            ExactP(zeros(2, 2), nothing),
            NewtonKrylov(refresh = refresh, linearsolver = ls); atol = 1e-10)
        @test info.converged
        @test isapprox(x, [sqrt(2), cbrt(3)]; rtol = 1e-6)
        @test length(info.krylov) >= info.iterations
    end

    # a residual whose norm is not finite, from an entry or from finite
    # entries whose norm overflows, scales no relative tolerance: the solve
    # ends at the start, not converged, before any linear solve
    for bad in (Inf, NaN, 1.5e308)
        info = JosephsonCircuits.nlsolvekrylov!((F, J, x) -> (isnothing(F) ||
                fill!(F, bad); nothing), jvp!, zeros(2), [1.0, 1.0],
            ExactP(zeros(2, 2), nothing); rtol = 1e-6)
        @test !info.converged
        @test info.reason === :nonfinite
        @test isempty(info.krylov)
    end

    # A rebuild forced by an escalation is never left to the probe: the
    # escalation drops the factorization, and a probe which applied it
    # would find nothing there. The preconditioner is the Jacobian scaled
    # by two, so its one step reduction is exactly one half at every point
    # and the probe, seeing no change since the last rebuild, always
    # declines to rebuild; the linear solver never converges and returns
    # the half Newton step, so every solve is a failure and the escalation
    # fires while the refresh reason is still `:stale`. Once escalated the
    # preconditioner is exact and the solver converges.
    mutable struct EscalatingP <: JosephsonCircuits.AbstractPreconditioner
        valid::Bool
        escalated::Bool
        refusals::Int
    end
    jac(x) = [2x[1] 0.0; 0.0 3x[2]^2]
    JosephsonCircuits.updatepreconditioner!(pc::EscalatingP, x::AbstractVector) =
        (pc.valid = true; pc)
    JosephsonCircuits.applypreconditioner!(z::AbstractVector, pc::EscalatingP,
        r::AbstractVector) = (pc.valid ||
            error("the preconditioner was applied after its escalation dropped the factorization");
        z .= (pc.escalated ? 1.0 : 0.5) .* (jac(xpt) \ r); z)
    JosephsonCircuits.escalatepreconditioner!(pc::EscalatingP) =
        pc.escalated ? false : pc.refusals < 1 ? (pc.refusals += 1; false) :
            (pc.escalated = true; pc.valid = false; true)
    struct HalfStepLS <: JosephsonCircuits.AbstractHBLinearSolver
        pc::EscalatingP
    end
    function JosephsonCircuits.hblinearsolve!(ls::HalfStepLS, deltax, jvp!, F,
            ws, Mop!; rtol, atol, maxrestarts, oncycle = nothing)
        ex = ls.pc.escalated
        deltax .= (ex ? 1.0 : 0.5) .* (jac(xpt) \ F)
        return (converged = ex, residual = ex ? 0.0 : norm(F), iterations = 1,
            cycles = 1, reason = ex ? :converged : :notconverged,
            precondtime = NaN)
    end
    for refresh in (Probe(), Always())
        x = [1.0, 1.0]
        F = zeros(2)
        pc = EscalatingP(false, false, 0)
        info = JosephsonCircuits.nlsolvekrylov!(fj!, jvp!, F, x, pc,
            NewtonKrylov(refresh = refresh, escalate = true,
                linearsolver = HalfStepLS(pc)); atol = 1e-10)
        @test info.converged
        @test pc.escalated
        @test isapprox(x, [sqrt(2), cbrt(3)]; rtol = 1e-6)
    end

    # a wrapper forwards every hook it does not define to what it wraps,
    # escalation included, so wrapping never turns a hook off
    mutable struct CountingP <: JosephsonCircuits.AbstractPreconditioner
        stalls::Int
        escalated::Bool
    end
    JosephsonCircuits.updatepreconditioner!(pc::CountingP, ::AbstractVector) = pc
    JosephsonCircuits.applypreconditioner!(z::AbstractVector, pc::CountingP,
        r::AbstractVector) = copyto!(z, r)
    JosephsonCircuits.stalled!(pc::CountingP) = (pc.stalls += 1; pc)
    JosephsonCircuits.escalatepreconditioner!(pc::CountingP) = (pc.escalated = true; true)
    JosephsonCircuits.deflationsize(::CountingP) = 7
    struct PlainWrap{P} <: JosephsonCircuits.AbstractWrappedPreconditioner
        inner::P
    end
    JosephsonCircuits.innerpreconditioner(w::PlainWrap) = w.inner
    JosephsonCircuits.updatepreconditioner!(w::PlainWrap, x::AbstractVector) =
        (JosephsonCircuits.updatepreconditioner!(w.inner, x); w)
    JosephsonCircuits.applypreconditioner!(z::AbstractVector, w::PlainWrap,
        r::AbstractVector) = JosephsonCircuits.applypreconditioner!(z, w.inner, r)
    inner = CountingP(0, false)
    w = PlainWrap(PlainWrap(inner))
    JosephsonCircuits.stalled!(w)
    @test inner.stalls == 1
    @test JosephsonCircuits.escalatepreconditioner!(w)
    @test inner.escalated
    @test JosephsonCircuits.deflationsize(w) == 7
    @test !JosephsonCircuits.isexactpreconditioner(w)
    @test JosephsonCircuits.pointmoved!(w) === w

    # under `Never` no linear solve is reported slow, whatever its rate: a
    # linear solver whose residual can grow, as BiCGStab's can, reports a
    # rate above one, which `Always` reports and `Never` does not. The
    # system is linear and the preconditioner exact, so the stagnated solve
    # is replaced by the preconditioner's, which is the root
    mutable struct SlowReportsP <: JosephsonCircuits.AbstractPreconditioner
        A::Matrix{Float64}
        stalls::Int
    end
    JosephsonCircuits.updatepreconditioner!(pc::SlowReportsP, ::AbstractVector) = pc
    JosephsonCircuits.applypreconditioner!(z::AbstractVector, pc::SlowReportsP,
        r::AbstractVector) = ldiv!(z, lu(pc.A), r)
    JosephsonCircuits.stalled!(pc::SlowReportsP) = (pc.stalls += 1; pc)
    struct GrowingResidualLS <: JosephsonCircuits.AbstractHBLinearSolver end
    function JosephsonCircuits.hblinearsolve!(::GrowingResidualLS, deltax,
            jvp!, F, ws, Mop!; rtol, atol, maxrestarts, oncycle = nothing)
        fill!(deltax, 0)
        return (converged = false, residual = 2*norm(F), iterations = 1,
            cycles = 1, reason = :notconverged)
    end
    Alin = [3.0 1.0; 1.0 2.0]
    blin = [1.0, -1.0]
    linres!(F, J, x) = (isnothing(F) || (mul!(F, Alin, x); F .-= blin); nothing)
    linjvp!(y, v) = mul!(y, Alin, v)
    for (refresh, reported) in ((Always(), true), (Never(), false))
        pc = SlowReportsP(Alin, 0)
        info = JosephsonCircuits.nlsolvekrylov!(linres!, linjvp!, zeros(2),
            zeros(2), pc, NewtonKrylov(refresh = refresh, escalate = false,
                linearsolver = GrowingResidualLS()); atol = 1e-10)
        @test info.converged
        @test (pc.stalls > 0) == reported
    end

    # the objects validate their own options
    @test_throws ArgumentError Staged(maxattempts = 0)
    @test_throws ArgumentError Staged(smin = 0.0)
    @test_throws ArgumentError Staged(s0 = 0.1, smin = 0.5)
    @test_throws ArgumentError Staged(interioriterations = 0)
    @test_throws ArgumentError QuasiNewton(anderson = -1)
    @test_throws ArgumentError Newton(factorization = BlockFactorization())
    @test_throws ArgumentError QuasiNewton(factorization = BlockFactorization())
    @test_throws ArgumentError GMRES(restart = 0)
    @test_throws ArgumentError GMRES(maxrestarts = 0)
    @test_throws ArgumentError GMRES(400, 0)
    @test_throws ArgumentError NewtonKrylov(linearsolver = 1)
    @test_throws TypeError NewtonKrylov(refresh = :always)
    @test_throws ArgumentError Floquet(size = 0)
    @test_throws ArgumentError Floquet(Floquet())
    @test_throws ArgumentError HarmonicBand(-1)
    @test_throws ArgumentError HarmonicBand((2, -1))
    @test_throws ArgumentError Staged(grids = [(0,), (8,)])
    # and the loop its tolerances and the forcing sequence it is handed
    @test_throws ArgumentError JosephsonCircuits.nlsolvekrylov!(fj!, jvp!,
        zeros(2), [1.0, 1.0], ExactP(zeros(2, 2), nothing); rtol = NaN)
    @test_throws ArgumentError JosephsonCircuits.nlsolvekrylov!(fj!, jvp!,
        zeros(2), [1.0, 1.0], ExactP(zeros(2, 2), nothing); forcingmax = 1.0)
end

@testset "a single precision solve at the default tolerance" begin
    # the rounding floor the tolerance is raised to is that of the
    # precision the residual is evaluated in, so a single precision solve
    # stops where single precision does, above the default tolerance which
    # a double precision solve reaches, at the double precision solve's
    # point to what single precision holds
    circuit, defs = testjpacircuit()
    wp = (2*pi*4.75001e9,)
    src = [(mode = (1,), port = 1, current = 0.00565e-6)]
    ref = hbnlsolve(wp, (8,), src, circuit, defs; keyedarrays = false,
        atol = 1e-12)
    single = hbnlsolve(wp, (8,), src, circuit, defs; keyedarrays = false,
        method = NewtonKrylov(precision = Float32))
    @test single.solverinfo.converged
    @test single.solverinfo.stages[end].iterations < 10
    @test single.solverinfo.finalresidual > 1e-8
    @test isapprox(single.nodeflux, ref.nodeflux; rtol = 1e-3)
end

# The direct current block is held in Float64 whatever precision the
# periodic solve runs in, because it is small, exactly solved, and the
# worst conditioned part of the problem. A single precision solve is then
# as accurate as single precision allows, and needs a tolerance it can
# meet: the default is absolute and sized for double precision.
@testset "a single precision solve carries the direct current block" begin
    c = Circuit(
        [:p1 => Port(1; Z0 = 50.0), :c1 => Capacitor(1.0e-12)],
        [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]])
    src = [(mode=(0,), port=1, current=1e-6),
           (mode=(1,), port=1, current=1e-6)]
    go(P; kw...) = hbnlsolve((2*pi*5e9,), (4,), src, c, Dict{Any,Any}();
        keyedarrays = false, dc = true, odd = true, even = true,
        method = NewtonKrylov(precision = P), kw...)

    a = go(Float64)
    b = go(Float32; rtol = 1e-6)
    @test a.solverinfo.converged
    @test b.solverinfo.converged
    # the same answer, to what single precision can hold
    @test isapprox(only(a.dcnodevoltage), only(b.dcnodevoltage);
        rtol = 1e-6)
    # and it stops where single precision runs out rather than at the
    # tolerance it was handed: `rtol` would accept a relative residual
    # of 1e-6 and the arithmetic reaches a few times `eps(Float32)`, so
    # the bound sits between them rather than at one ulp, which a
    # residual summed over the modes does not land on exactly
    r = b.solverinfo.finalresidual/b.solverinfo.initialresidual
    @test r <= 4*eps(Float32)
end

@testset "every solve ends with a reason" begin
    circuit, defs = testchaincircuit()
    wp = 2*pi*4.75e9
    src = [(mode=(1,), port=1, current=0.3e-6)]
    ok = hbnlsolve((wp,), (8,), src, circuit, defs)
    @test ok.solverinfo.converged
    @test ok.solverinfo.stages[end].reason == :converged
    for m in (NewtonKrylov(), Newton(), QuasiNewton())
        spent = @test_logs (:warn,) match_mode=:any hbnlsolve(
            (wp,), (8,), src, circuit, defs; iterations = 1, method = m)
        @test !spent.solverinfo.converged
        @test spent.solverinfo.stages[end].reason == :iterations
    end
    # a work budget of two Arnoldi steps, one restart length per Newton
    # step, which a block diagonal preconditioner on a longer chain spends
    # within its first two linear solves
    long, _ = testchaincircuit(12)
    work = @test_logs (:warn,) hbnlsolve((2*pi*8e9,), (8,),
        [(mode=(1,), port=1, current=3.2e-6)], long, defs; iterations = 2,
        method = NewtonKrylov(preconditioner = BlockDiagonal(),
            escalate = false,
            linearsolver = GMRES(restart = 1, maxrestarts = 20)))
    @test !work.solverinfo.converged
    @test work.solverinfo.stages[end].reason === :work
    # a drive far beyond the self oscillation threshold has no operating
    # point the direct solvers can reach: they stop and say so, well within
    # their budgets, rather than spending them
    hard = [(mode=(1,), port=1, current=60e-6)]
    for m in (NewtonKrylov(preconditioner = BlockDiagonal(), escalate = false),
            NewtonKrylov(), Newton(), QuasiNewton())
        r = @test_logs (:warn,) match_mode=:any hbnlsolve(
            (wp,), (8,), hard, circuit, defs; iterations = 400, method = m)
        @test !r.solverinfo.converged
        @test r.solverinfo.stages[end].reason in (:linesearch, :progress, :work)
        @test r.solverinfo.stages[end].iterations < 400
        # the direct loop judges its creep against the steps left too, so
        # it ends within a few stall windows rather than near its budget
        m isa Union{Newton,QuasiNewton} &&
            @test r.solverinfo.stages[end].iterations < 100
    end
end

@testset "a step with no decrease is retried only when a rebuild can change it" begin
    # with no backtrack the first full step of a strong drive raises the
    # merit; the full Jacobian, rebuilt at that point, is exact, so a retry
    # from it would repeat the step, and the solve ends after one search
    circuit, defs = testjpacircuit()
    r = @test_logs (:warn,) match_mode = :any hbnlsolve((2*pi*5e9,), (8,),
        [(mode = (1,), port = 1, current = 2e-6)], circuit, defs;
        method = NewtonKrylov(linesearch = Backtracking(maxbacktracks = 0)))
    st = r.solverinfo.stages[end]
    @test st.reason === :linesearch
    @test count(iszero, st.alpha) == 1
end

@testset "a probed step builds its deflation once" begin
    # the probe measures the preconditioner the last step left, and the
    # deflation is built once a step, by the refresh or at the first
    # application after the point moved, never by both. The probe decides
    # on measured times, so the solve is repeated once compiled
    long, defs = testchaincircuit(12)
    solve() = hbnlsolve((2*pi*8e9,), (8,), [(mode = (1,), port = 1,
        current = 3.2e-6)], long, defs; method = NewtonKrylov(
        preconditioner = Floquet(size = 12, harvest = 4), refresh = Probe(),
        escalate = false))
    solve()
    r = solve()
    st = r.solverinfo.stages[end]
    @test r.solverinfo.converged
    @test last(st.krylov).deflationrebuilds <= length(st.krylov)
end
