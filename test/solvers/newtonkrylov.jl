using JosephsonCircuits, LinearAlgebra, SparseArrays, Random, Test, Logging

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
            NewtonKrylov(refresh = refresh, linearsolver = ls); ftol = 1e-10)
        @test info.converged
        @test isapprox(x, [sqrt(2), cbrt(3)]; rtol = 1e-6)
        @test length(info.krylov) >= info.iterations
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
                linearsolver = HalfStepLS(pc)); ftol = 1e-10)
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

    # the objects validate their own options
    @test_throws ArgumentError Staged(maxattempts = 0)
    @test_throws ArgumentError Staged(smin = 0.0)
    @test_throws ArgumentError Staged(s0 = 0.1, smin = 0.5)
    @test_throws ArgumentError Staged(interioriterations = 0)
    @test_throws ArgumentError Newton(factorization = BlockFactorization())
    @test_throws ArgumentError QuasiNewton(factorization = BlockFactorization())
    @test_throws ArgumentError GMRES(restart = 0)
    @test_throws ArgumentError GMRES(maxrestarts = 0)
    @test_throws ArgumentError NewtonKrylov(linearsolver = 1)
    @test_throws TypeError NewtonKrylov(refresh = :always)
    @test_throws ArgumentError Floquet(size = 0)
    @test_throws ArgumentError Floquet(Floquet())
end

@testset "every solve ends with a reason" begin
    circuit, defs = testchaincircuit()
    wp = 2*pi*4.75e9
    src = [(mode=(1,), port=1, current=0.3e-6)]
    ok = hbnlsolve((wp,), (8,), src, circuit, defs)
    @test ok.solverinfo.converged
    @test ok.solverinfo.stages[end].reason == :converged
    for m in (NewtonKrylov(), Newton(), QuasiNewton())
        spent = @test_logs (:warn, r"did not converge: the Newton iteration budget") match_mode=:any hbnlsolve(
            (wp,), (8,), src, circuit, defs; iterations = 1, method = m)
        @test !spent.solverinfo.converged
        @test spent.solverinfo.stages[end].reason == :iterations
    end
    # a work budget of one Arnoldi step per Newton step, spent on the
    # first. Whether the budget or the iteration count runs out first, or
    # the one step is enough after all, is not fixed, so the warning the
    # solver makes when it gives up is taken rather than asserted; a
    # message at error level would still fail the test
    work = @test_logs min_level=Logging.Error hbnlsolve((wp,), (8,), src, circuit, defs;
        iterations = 2,
        method = NewtonKrylov(linearsolver = GMRES(restart = 1, maxrestarts = 1)))
    @test work.solverinfo.stages[end].reason in (:work, :iterations, :converged)
    # a drive far beyond the self oscillation threshold has no operating
    # point the direct solvers can reach: they stop and say so, well within
    # their budgets, rather than spending them
    hard = [(mode=(1,), port=1, current=60e-6)]
    for m in (NewtonKrylov(preconditioner = BlockDiagonal(), escalate = false),
            NewtonKrylov(), Newton())
        r = @test_logs (:warn, r"did not converge") match_mode=:any hbnlsolve(
            (wp,), (8,), hard, circuit, defs; iterations = 400, method = m)
        @test !r.solverinfo.converged
        @test r.solverinfo.stages[end].reason in (:linesearch, :progress, :work)
        @test r.solverinfo.stages[end].iterations < 400
    end
end
