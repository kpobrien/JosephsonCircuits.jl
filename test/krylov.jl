using JosephsonCircuits, LinearAlgebra, SparseArrays, Test

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
