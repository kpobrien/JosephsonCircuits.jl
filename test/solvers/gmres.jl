using JosephsonCircuits, LinearAlgebra, SparseArrays, Random, Test, Logging

# The product with a matrix and the solve with a factorization, in the
# form the solver takes them. The tests make their operators here, so the
# solver is compiled once for each type of matrix rather than once for
# each test.
matrixproduct(A) = (w, v) -> mul!(w, A, v)
factorsolve(F) = (z, v) -> ldiv!(z, F, v)

# the inverse of a diagonal, as a preconditioner of the package's interface
struct JacobiP <: JosephsonCircuits.AbstractPreconditioner
    d::Vector{Float64}
end
JosephsonCircuits.applypreconditioner!(z::AbstractVector, p::JacobiP,
    r::AbstractVector) = (z .= r ./ p.d; z)

# The Krylov solver itself: the basis it grows, the systems it is exact
# on, restarting, preconditioning, its allocation, and the breakdowns and
# arguments it refuses.
@testset verbose=true "krylov" begin

    @testset "argument checking" begin
        @test_throws ArgumentError JosephsonCircuits.GMRESWorkspace(5, 0)
        @test_throws ArgumentError JosephsonCircuits.GMRESWorkspace(-1, 3)
        ws = JosephsonCircuits.GMRESWorkspace(5, 3)
        A = Matrix(1.0I, 5, 5)
        op! = matrixproduct(A)
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
        out = JosephsonCircuits.gmres!(x, matrixproduct(A), b, ws;
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
        JosephsonCircuits.gmres!(y, matrixproduct(A), b, full;
            rtol = 1e-12, maxrestarts = 1)
        @test x == y
        # growth stops at the restart length: a basis asked for more keeps
        # its width
        @test size(JosephsonCircuits.ensurecolumns!(ws, 10*m), 2) == m + 1
    end

    @testset "zero right hand side" begin
        ws = JosephsonCircuits.GMRESWorkspace(4, 3)
        A = Matrix(2.0I, 4, 4)
        out = JosephsonCircuits.gmres!(ones(4), matrixproduct(A), zeros(4), ws)
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
            out = JosephsonCircuits.gmres!(x, matrixproduct(A), b, ws;
                rtol = 1e-12, maxrestarts = 4)
            @test out.converged
            @test isapprox(x, xref; rtol = 1e-8)
            @test norm(b - A*x) <= 1e-10*norm(b)
        end
    end

    @testset "exact in at most n steps" begin
        # unpreconditioned GMRES on an n x n system reaches the exact solution
        # in at most n Arnoldi steps in exact arithmetic; the restart length
        # is twice that, so it is the solver which stops, at the breakdown
        # of the n dimensional Krylov space at the latest
        n = 12
        A = randn(n, n) + 3n*I
        b = randn(n)
        ws = JosephsonCircuits.GMRESWorkspace(n, 2n)
        x = zeros(n)
        out = JosephsonCircuits.gmres!(x, matrixproduct(A), b, ws;
            rtol = 1e-13, maxrestarts = 1)
        @test out.converged
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
            out = JosephsonCircuits.gmres!(x, matrixproduct(A), b, ws;
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
        Aop! = matrixproduct(A)

        ws1 = JosephsonCircuits.GMRESWorkspace(n, 40)
        x1 = zeros(n)
        out1 = JosephsonCircuits.gmres!(x1, Aop!, b, ws1; rtol = 1e-10,
            maxrestarts = 30)

        # exact preconditioner: one iteration, and the answer is the direct one
        F = lu(A)
        Mop! = factorsolve(F)
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
        Mop! = factorsolve(F)
        b = randn(n)
        for pert in (0.0, 0.05, 0.25)
            A = A0 + pert*spdiagm(0 => randn(n))
            ws = JosephsonCircuits.GMRESWorkspace(n, 30)
            x = zeros(n)
            out = JosephsonCircuits.gmres!(x, matrixproduct(A), b, ws;
                Mop! = Mop!, rtol = 1e-10, maxrestarts = 20)
            @test out.converged
            @test isapprox(x, A \ b; rtol = 1e-6)
        end
    end

    @testset "allocation does not scale with iterations" begin
        n = 50
        A = randn(n, n) + 5n*I
        b = randn(n)
        Aop! = matrixproduct(A)
        # a small fixed overhead remains from the views handed to the BLAS
        # calls; what matters is that it does not grow with the iteration
        # count, so compare a short solve against a long one
        # (measured inside a function: `@allocated` has Julia compile the
        # whole top-level expression it is written in, here the file's
        # testset)
        allocations(x, ws, rtol, maxrestarts) = @allocated JosephsonCircuits.gmres!(
            x, Aop!, b, ws; rtol, maxrestarts)
        ws1 = JosephsonCircuits.GMRESWorkspace(n, 3)
        x1 = zeros(n)
        JosephsonCircuits.gmres!(x1, Aop!, b, ws1; rtol = 1e-10, maxrestarts = 1)
        short = allocations(x1, ws1, 1e-10, 1)

        ws2 = JosephsonCircuits.GMRESWorkspace(n, 40)
        x2 = zeros(n)
        JosephsonCircuits.gmres!(x2, Aop!, b, ws2; rtol = 1e-12, maxrestarts = 20)
        long = allocations(x2, ws2, 1e-12, 20)

        @test short <= 1024
        @test long <= 1024
    end

    @testset "the erased operator and preconditioner allocate nothing per step" begin
        # the Newton-Krylov iteration hands GMRES its product and its
        # preconditioner erased, and every Arnoldi step calls both through
        # a dynamic dispatch, whose arguments must not be boxed
        n = 50
        A = randn(n, n) + 2*sqrt(n)*I
        b = randn(n)
        Aop = JosephsonCircuits.asoperator(
            JosephsonCircuits.erased(matrixproduct(A)), n)
        M = JosephsonCircuits.erased(JacobiP(diag(A)))
        allocations(x, ws, rtol, maxrestarts) = @allocated JosephsonCircuits.gmres!(
            x, Aop, b, ws; Mop! = M, rtol, maxrestarts)
        ws1 = JosephsonCircuits.GMRESWorkspace(n, 3)
        x1 = zeros(n)
        JosephsonCircuits.gmres!(x1, Aop, b, ws1; Mop! = M, rtol = 1e-12,
            maxrestarts = 1)
        short = allocations(x1, ws1, 1e-12, 1)
        ws2 = JosephsonCircuits.GMRESWorkspace(n, 40)
        x2 = zeros(n)
        out = JosephsonCircuits.gmres!(x2, Aop, b, ws2; Mop! = M,
            rtol = 1e-12, maxrestarts = 4)
        long = allocations(x2, ws2, 1e-12, 4)
        # the long solve takes ten times the steps of the short one
        @test out.iterations >= 30
        @test long <= short
        @test norm(A*x2 - b) <= 1e-10*norm(b)
    end

    @testset "the norm is formed through the inner product in double alone" begin
        # a shorter float keeps the scaled norm, whose squares cannot
        # overflow: 400^2 exceeds the largest Float16
        @test JosephsonCircuits.norm2(Float16[300, 400]) == Float16(500)
        @test JosephsonCircuits.norm2([3.0, 4.0]) == 5.0
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
    out = JosephsonCircuits.gmres!(x, matrixproduct(A), b, ws;
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
    out2 = JosephsonCircuits.gmres!(x2, matrixproduct(A2), b2, ws2;
        rtol = 1e-8, maxrestarts = 3)
    @test !out2.converged
    @test all(isfinite, x2)
    @test norm(x2) < 10                # not the 1e16 a raw divide produces
    # and it is at least as good as the zero step
    @test norm(b2 - A2*x2) <= norm(b2)

    # the solve does not depend on the scale of the operator: an operator
    # scaled by 1e-18 gives the solution of the unscaled one scaled by 1e18
    Ar = Matrix(Diagonal(1.0 .+ (1:20) ./ 20)) .+ 0.05 .* sin.((1:20) .* (1:20)')
    br = ones(20)
    xr1 = zeros(20); xrs = zeros(20)
    wsr = JosephsonCircuits.GMRESWorkspace(20, 20, Float64)
    outr1 = JosephsonCircuits.gmres!(xr1, matrixproduct(Ar), br, wsr;
        rtol = 1e-10, maxrestarts = 3)
    outrs = JosephsonCircuits.gmres!(xrs, matrixproduct(1e-18*Ar), br,
        wsr; rtol = 1e-10, maxrestarts = 3)
    @test outr1.converged && outrs.converged
    @test isapprox(1e-18*xrs, xr1; rtol = 1e-9)

    # a preconditioner which is not fixed makes the recurrence estimate run
    # ahead of the explicit residual, so a cycle ends early and another
    # follows, whose lengths the per cycle callback sees
    calls = Ref(0)
    Mv!(z, v) = (calls[] += 1; z .= v .* (1 + 0.3*sin(calls[])); z)
    Av = Matrix(Diagonal(1.0 .+ (1:40) ./ 40)) .+ 0.05 .* sin.((1:40) .* (1:40)')
    wsv = JosephsonCircuits.GMRESWorkspace(40, 30, Float64)
    lengths = Int[]
    outv = JosephsonCircuits.gmres!(zeros(40), matrixproduct(Av),
        ones(40), wsv; Mop! = Mv!, rtol = 1e-10, maxrestarts = 6,
        oncycle = (ws, j) -> push!(lengths, j))
    @test outv.cycles == length(lengths) >= 2
    @test first(lengths) < 30

    # a non-finite value from the preconditioner or the operator ends the
    # cycle at the step which met it: the solve does not run the cycle out
    # on a spoiled basis, `x` stays where the cycle started, and the
    # residual reported is that point's
    An = Matrix(Diagonal(1.0 .+ (1:40) ./ 40)) .+ 0.05 .* sin.((1:40) .* (1:40)')
    bn = ones(40)
    for poisoned in (:preconditioner, :operator)
        calls = Ref(0)
        Mn!(z, v) = (calls[] += 1; z .= v;
            poisoned === :preconditioner && calls[] == 3 && (z[5] = NaN); z)
        An!(w, v) = (mul!(w, An, v);
            poisoned === :operator && calls[] == 3 && (w[1] = Inf); w)
        xn = zeros(40)
        outn = JosephsonCircuits.gmres!(xn, An!, bn,
            JosephsonCircuits.GMRESWorkspace(40, 30, Float64); Mop! = Mn!,
            rtol = 1e-12, maxrestarts = 4)
        @test !outn.converged
        @test outn.reason == :stagnation
        @test outn.iterations == 3
        @test xn == zeros(40)
        @test outn.residual == norm(bn)
        @test outn.residualvector ≈ bn
    end

    # a cycle of one step against an exact preconditioner takes its
    # correction from the application its step made, so the solve applies
    # the preconditioner once
    Fe = lu(An)
    applied = Ref(0)
    Me!(z, v) = (applied[] += 1; ldiv!(z, Fe, v))
    xe = zeros(40)
    oute = JosephsonCircuits.gmres!(xe, matrixproduct(An), bn,
        JosephsonCircuits.GMRESWorkspace(40, 30, Float64); Mop! = Me!,
        rtol = 1e-12, maxrestarts = 4)
    @test oute.converged && oute.cycles == 1 && oute.iterations == 1
    @test applied[] == 1
    @test xe ≈ An \ bn rtol = 1e-12
    # a longer cycle of rank one does not: its last step left the image of
    # its second basis vector in the workspace. Here the second step is a
    # breakdown of [1 0; 0 0], and the correction is the preconditioned
    # multiple of the first basis vector, x = D*[1, 1]
    D = [1.0, 2.0]
    xd = zeros(2)
    outd = JosephsonCircuits.gmres!(xd, matrixproduct([1.0 0.0; 0.0 0.0]),
        [1.0, 1.0], JosephsonCircuits.GMRESWorkspace(2, 2, Float64);
        Mop! = (z, v) -> (z .= D .* v), rtol = 1e-12, maxrestarts = 1)
    @test outd.cycles == 1 && outd.iterations == 2
    @test xd ≈ [1.0, 2.0]

    # non-finite tolerances must be rejected rather than reporting a
    # convergence which did not happen
    A3 = [1.0 0.0; 0.0 1.0]
    b3 = [1.0, 1.0]
    x3 = zeros(2)
    ws3 = JosephsonCircuits.GMRESWorkspace(2, 2, Float64)
    @test_throws ArgumentError JosephsonCircuits.gmres!(x3,
        matrixproduct(A3), b3, ws3; rtol = Inf)
    @test_throws ArgumentError JosephsonCircuits.gmres!(x3,
        matrixproduct(A3), b3, ws3; atol = NaN)
    # and so must a right hand side whose norm is not finite, from an entry
    # or from finite entries whose norm overflows: its relative tolerance
    # would accept any residual, the zero start's infinite one included
    for bnf in ([Inf, 1.0], fill(1.5e308, 2), [NaN, 1.0])
        xnf = ones(2)
        outnf = JosephsonCircuits.gmres!(xnf, matrixproduct(A3), bnf, ws3)
        @test !outnf.converged
        @test outnf.reason === :nonfinite
        @test iszero(xnf)
    end

    # the reported termination reason and cycle count are present and
    # consistent with convergence
    A5 = Diagonal(exp10.(range(-1, 0, length = 20)))
    b5 = ones(20)
    x5 = zeros(20)
    ws5 = JosephsonCircuits.GMRESWorkspace(20, 20, Float64)
    out5 = JosephsonCircuits.gmres!(x5, matrixproduct(A5), b5, ws5;
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
    op = JosephsonCircuits.asoperator(matrixproduct(A5), 20)
    fromproduct = JosephsonCircuits.meritslope!(Jv, op, p6, F6, ϕ0, nothing)
    @test fromproduct ≈ dot(F6, A5*p6)
    fromresidual = JosephsonCircuits.meritslope!(Jv, op, p6, F6, ϕ0,
        out5.residualvector)
    @test fromresidual ≈ fromproduct
end
