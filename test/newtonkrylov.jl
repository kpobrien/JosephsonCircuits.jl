using JosephsonCircuits, LinearAlgebra, SparseArrays, Random, Test

isdefined(Main, :testchaincircuit) || include("testcircuits.jl")
Random.seed!(20260905)

@testset verbose=true "krylov" begin

    @testset "argument checking" begin
        @test_throws ArgumentError JosephsonCircuits.GMRESWorkspace(5, 0)
        @test_throws ArgumentError JosephsonCircuits.GMRESWorkspace(-1, 3)
        ws = JosephsonCircuits.GMRESWorkspace(5, 3)
        A = Matrix(1.0I, 5, 5)
        op!(w, v) = mul!(w, A, v)
        @test_throws DimensionMismatch JosephsonCircuits.gmres!(
            zeros(4), op!, zeros(5), ws)
        @test_throws DimensionMismatch JosephsonCircuits.gmres!(
            zeros(6), op!, zeros(6), ws)
        @test_throws ArgumentError JosephsonCircuits.gmres!(
            zeros(5), op!, ones(5), ws; rtol = -1.0)
        @test_throws ArgumentError JosephsonCircuits.gmres!(
            zeros(5), op!, ones(5), ws; maxrestarts = 0)
    end

    @testset "the basis grows to what the iteration uses" begin
        # a workspace for a long restart is born with a few columns and
        # widens as the Arnoldi steps need them, so a solve which takes a
        # handful of steps never touches a basis it would not fill
        n, m = 200, 100
        A = Matrix(Diagonal(1.0 .+ (1:n)./n)) + 0.01*randn(n, n)
        b = ones(n)
        ws = JosephsonCircuits.GMRESWorkspace(b, m)
        @test size(ws.V) == (n, JosephsonCircuits.GMRESINITIALCOLUMNS)
        x = zeros(n)
        out = JosephsonCircuits.gmres!(x, (w, v) -> mul!(w, A, v), b, ws;
            rtol = 1e-12, maxrestarts = 1)
        @test out.converged
        @test out.iterations > JosephsonCircuits.GMRESINITIALCOLUMNS
        # wide enough for what was built, no wider than the restart allows
        @test size(ws.V, 2) >= out.iterations
        @test size(ws.V, 2) <= m + 1
        @test norm(A*x - b) <= 1e-10*norm(b)
        # and the same answer as a basis allocated in full
        full = JosephsonCircuits.GMRESWorkspace(b, m)
        full.V = similar(b, n, m + 1)
        y = zeros(n)
        JosephsonCircuits.gmres!(y, (w, v) -> mul!(w, A, v), b, full;
            rtol = 1e-12, maxrestarts = 1)
        @test x == y
        # growth stops at the restart length: a basis asked for more keeps
        # its width
        @test size(JosephsonCircuits.ensurecolumns!(ws, 10*m), 2) == m + 1
    end

    @testset "zero right hand side" begin
        ws = JosephsonCircuits.GMRESWorkspace(4, 3)
        A = Matrix(2.0I, 4, 4)
        out = JosephsonCircuits.gmres!(ones(4), (w, v) -> mul!(w, A, v), zeros(4), ws)
        @test out.converged
        @test out.iterations == 0
    end

    @testset "dense nonsymmetric systems" begin
        for n in (10, 40)
            A = randn(n, n) + 5n*I     # diagonally dominant, well conditioned
            b = randn(n)
            xref = A \ b
            ws = JosephsonCircuits.GMRESWorkspace(n, n)
            x = zeros(n)
            out = JosephsonCircuits.gmres!(x, (w, v) -> mul!(w, A, v), b, ws;
                rtol = 1e-12, maxrestarts = 4)
            @test out.converged
            @test isapprox(x, xref; rtol = 1e-8)
            @test norm(b - A*x) <= 1e-10*norm(b)
        end
    end

    @testset "exact in at most n steps" begin
        # unpreconditioned GMRES on an n x n system reaches the exact solution
        # in at most n Arnoldi steps in exact arithmetic
        n = 12
        A = randn(n, n) + 3n*I
        b = randn(n)
        ws = JosephsonCircuits.GMRESWorkspace(n, n)
        x = zeros(n)
        out = JosephsonCircuits.gmres!(x, (w, v) -> mul!(w, A, v), b, ws;
            rtol = 1e-13, maxrestarts = 1)
        @test out.iterations <= n
        @test isapprox(x, A \ b; rtol = 1e-7)
    end

    @testset "restarting reaches the same solution" begin
        n = 60
        A = randn(n, n) + 6n*I
        b = randn(n)
        xref = A \ b
        for m in (5, 15, n)
            ws = JosephsonCircuits.GMRESWorkspace(n, m)
            x = zeros(n)
            out = JosephsonCircuits.gmres!(x, (w, v) -> mul!(w, A, v), b, ws;
                rtol = 1e-11, maxrestarts = 50)
            @test out.converged
            @test isapprox(x, xref; rtol = 1e-7)
        end
    end

    @testset "preconditioning reduces iterations" begin
        n = 120
        # a badly scaled operator that GMRES struggles with unpreconditioned
        A = sprandn(n, n, 0.05) + spdiagm(0 => range(1.0, 500.0; length = n))
        b = randn(n)
        Aop!(w, v) = mul!(w, A, v)

        ws1 = JosephsonCircuits.GMRESWorkspace(n, 40)
        x1 = zeros(n)
        out1 = JosephsonCircuits.gmres!(x1, Aop!, b, ws1; rtol = 1e-10,
            maxrestarts = 30)

        # exact preconditioner: one iteration, and the answer is the direct one
        F = lu(A)
        Mop!(z, v) = ldiv!(z, F, v)
        ws2 = JosephsonCircuits.GMRESWorkspace(n, 40)
        x2 = zeros(n)
        out2 = JosephsonCircuits.gmres!(x2, Aop!, b, ws2; Mop! = Mop!,
            rtol = 1e-10, maxrestarts = 30)

        @test out1.converged && out2.converged
        @test out2.iterations < out1.iterations
        @test out2.iterations <= 2
        @test isapprox(x1, x2; rtol = 1e-6)
        @test isapprox(x2, A \ b; rtol = 1e-8)
    end

    @testset "stale preconditioner still converges" begin
        # the operating point moves but the factorization does not, which is
        # how the preconditioner is reused across Newton steps
        n = 80
        A0 = sprandn(n, n, 0.06) + spdiagm(0 => fill(30.0, n))
        F = lu(A0)
        Mop!(z, v) = ldiv!(z, F, v)
        b = randn(n)
        for pert in (0.0, 0.05, 0.25)
            A = A0 + pert*spdiagm(0 => randn(n))
            ws = JosephsonCircuits.GMRESWorkspace(n, 30)
            x = zeros(n)
            out = JosephsonCircuits.gmres!(x, (w, v) -> mul!(w, A, v), b, ws;
                Mop! = Mop!, rtol = 1e-10, maxrestarts = 20)
            @test out.converged
            @test isapprox(x, A \ b; rtol = 1e-6)
        end
    end

    @testset "warm start" begin
        n = 30
        A = randn(n, n) + 4n*I
        b = randn(n)
        xref = A \ b
        ws = JosephsonCircuits.GMRESWorkspace(n, n)
        x = copy(xref)                       # start at the solution
        out = JosephsonCircuits.gmres!(x, (w, v) -> mul!(w, A, v), b, ws;
            rtol = 1e-10, initialzero = false)
        @test out.converged
        @test out.iterations == 0            # nothing to do
        @test isapprox(x, xref; rtol = 1e-10)
    end

    @testset "allocation does not scale with iterations" begin
        n = 50
        A = randn(n, n) + 5n*I
        b = randn(n)
        Aop!(w, v) = mul!(w, A, v)
        # a small fixed overhead remains from the views handed to the BLAS
        # calls; what matters is that it does not grow with the iteration
        # count, so compare a short solve against a long one
        ws1 = JosephsonCircuits.GMRESWorkspace(n, 3)
        x1 = zeros(n)
        JosephsonCircuits.gmres!(x1, Aop!, b, ws1; rtol = 1e-10, maxrestarts = 1)
        short = @allocated JosephsonCircuits.gmres!(x1, Aop!, b, ws1;
            rtol = 1e-10, maxrestarts = 1)

        ws2 = JosephsonCircuits.GMRESWorkspace(n, 40)
        x2 = zeros(n)
        JosephsonCircuits.gmres!(x2, Aop!, b, ws2; rtol = 1e-12, maxrestarts = 20)
        long = @allocated JosephsonCircuits.gmres!(x2, Aop!, b, ws2;
            rtol = 1e-12, maxrestarts = 20)

        @test short <= 1024
        @test long <= 1024
    end

end

@testset "gmres! singular, breakdown and validation" begin
    # A = 0 with a nonzero right hand side is an unhappy breakdown: the
    # Krylov space is invariant on the first step but the residual is not
    # reducible in it. The solver must not report success, must not call it
    # lucky, and must not spend its whole cycle budget rebuilding the same
    # useless space.
    A = zeros(3, 3)
    b = [1.0, 2.0, 3.0]
    x = zeros(3)
    ws = JosephsonCircuits.GMRESWorkspace(3, 3, Float64)
    out = JosephsonCircuits.gmres!(x, (w, v) -> mul!(w, A, v), b, ws;
        rtol = 1e-8, maxrestarts = 5)
    @test !out.converged
    @test out.reason != :converged
    @test out.iterations <= 2          # not one wasted cycle per restart
    @test out.residual ≈ norm(b)
    @test all(isfinite, x)

    # a rank deficient, inconsistent system must not amplify a tiny
    # triangular pivot into a spurious solution component
    A2 = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 0.0]
    b2 = [1.0, 1.0, 1.0]
    x2 = zeros(3)
    ws2 = JosephsonCircuits.GMRESWorkspace(3, 3, Float64)
    out2 = JosephsonCircuits.gmres!(x2, (w, v) -> mul!(w, A2, v), b2, ws2;
        rtol = 1e-8, maxrestarts = 3)
    @test !out2.converged
    @test all(isfinite, x2)
    @test norm(x2) < 10                # not the 1e16 a raw divide produces
    # and it is at least as good as the zero step
    @test norm(b2 - A2*x2) <= norm(b2)

    # non-finite tolerances must be rejected rather than reporting a
    # convergence which did not happen
    A3 = [1.0 0.0; 0.0 1.0]
    b3 = [1.0, 1.0]
    x3 = zeros(2)
    ws3 = JosephsonCircuits.GMRESWorkspace(2, 2, Float64)
    @test_throws ArgumentError JosephsonCircuits.gmres!(x3,
        (w, v) -> mul!(w, A3, v), b3, ws3; rtol = Inf)
    @test_throws ArgumentError JosephsonCircuits.gmres!(x3,
        (w, v) -> mul!(w, A3, v), b3, ws3; atol = NaN)

    # a zero right hand side with initialzero = false must measure the warm
    # start rather than discard it
    A4 = [2.0 0.0; 0.0 3.0]
    b4 = [0.0, 0.0]
    x4 = [7.0, -4.0]
    ws4 = JosephsonCircuits.GMRESWorkspace(2, 2, Float64)
    out4 = JosephsonCircuits.gmres!(x4, (w, v) -> mul!(w, A4, v), b4, ws4;
        rtol = 1e-10, initialzero = false)
    @test out4.converged
    @test norm(A4*x4 - b4) <= 1e-10

    # the reported termination reason and cycle count are present and
    # consistent with convergence
    A5 = Diagonal(exp10.(range(-1, 0, length = 20)))
    b5 = ones(20)
    x5 = zeros(20)
    ws5 = JosephsonCircuits.GMRESWorkspace(20, 20, Float64)
    out5 = JosephsonCircuits.gmres!(x5, (w, v) -> mul!(w, A5, v), b5, ws5;
        rtol = 1e-10, maxrestarts = 3)
    @test out5.converged
    @test out5.reason == :converged
    @test out5.cycles >= 1
    # the explicit residual of the returned solution comes back with it,
    # so a caller which needs `A x` does not take another product
    @test out5.residualvector ≈ b5 - A5*x5
    @test norm(out5.residualvector) == out5.residual

    # the merit function slope read off that residual is the one the
    # product gives: with p = -x, F'Jp = F'(w - F)
    F6 = b5; p6 = -x5
    ϕ0 = JosephsonCircuits.merit(F6)
    Jv = zeros(20)
    op = JosephsonCircuits.asoperator((y, v) -> mul!(y, A5, v), 20)
    fromproduct = JosephsonCircuits.meritslope!(Jv, op, p6, F6, ϕ0, nothing)
    @test fromproduct ≈ dot(F6, A5*p6)
    fromresidual = JosephsonCircuits.meritslope!(Jv, op, p6, F6, ϕ0,
        out5.residualvector)
    @test fromresidual ≈ fromproduct
end

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
    # a work budget of one Arnoldi step per Newton step, spent on the first
    work = hbnlsolve((wp,), (8,), src, circuit, defs; iterations = 2,
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
