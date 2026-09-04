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

    # the initial forcing term and the upper clamp are separate parameters
    for (rtol0, rtolmax) in ((0.3, 0.9), (0.9, 0.9), (0.01, 0.5))
        x = [1.0, 1.0]
        F = zeros(2)
        info = JosephsonCircuits.nlsolvekrylov!(fj!, jvp!, F, x,
            ExactP(zeros(2, 2), nothing); ftol = 1e-10,
            krylovrtol0 = rtol0, krylovrtolmax = rtolmax)
        @test info.converged
        @test isapprox(x, [sqrt(2), cbrt(3)]; rtol = 1e-6)
        @test length(info.krylov) >= info.iterations
    end

    # the Choice 2 parameters are only defined on their proper ranges
    for bad in ((; krylovgamma = 0.0), (; krylovgamma = 1.5),
                (; krylovalpha = 1.0), (; krylovalpha = 3.0),
                (; krylovrtol0 = 0.99, krylovrtolmax = 0.5),
                (; krylovrefreshrate = 0.0), (; krylovstagnation = 1.5))
        @test_throws ArgumentError JosephsonCircuits.nlsolvekrylov!(fj!, jvp!,
            zeros(2), [1.0, 1.0], ExactP(zeros(2, 2), nothing); bad...)
    end
end

@testset verbose=true "recycling deflation forms (adef1 / adef2)" begin
    using Random
    Random.seed!(20260902)
    JC = JosephsonCircuits

    # a base preconditioner which is a fixed factorization of `P`, wrapped
    # in the minimal AbstractPreconditioner interface
    struct FixedP{F} <: JC.AbstractPreconditioner
        F::F
    end
    JC.updatepreconditioner!(pc::FixedP, ::AbstractVector) = pc
    JC.applypreconditioner!(z::AbstractVector, pc::FixedP, r::AbstractVector) =
        ldiv!(z, pc.F, r)

    # The operator is `A = P*(I + Q3*D*Q3')` with orthonormal `Q3`, so that
    # `A*inv(P)` is similar to `I + Q3*D*Q3'`: its spectrum is one, `n - 3`
    # times, plus `1 + d_i`, with right eigenvectors `P*q_i`. Two of the
    # three are placed near zero, which is what stalls GMRES and what the
    # deflation exists to remove. The corrections of the unknowns along
    # those eigenvectors are `inv(P)*P*q_i = q_i`.
    n = 60
    P = randn(n, n) + 8I
    Q3 = Matrix(qr(randn(n, 3)).Q)[:, 1:3]
    d = [-0.99, -0.95, 2.0]
    A = P*(I + Q3*Diagonal(d)*Q3')
    Aref = Ref(A)
    jvp!(y, v) = mul!(y, Aref[], v)
    Aop = JC.asoperator(jvp!, n)
    base = FixedP(lu(P))
    b = randn(n)

    # the exact candidates for the two eigenvalues near zero: the right
    # eigenvectors of A*inv(P), whose corrections inv(P)*u span Q3[:, 1:2]
    Uexact = Matrix(qr(P*Q3[:, 1:2]).Q)[:, 1:2]
    Xexact = Q3[:, 1:2]

    function build(form; U = Uexact, kmax = 4)
        pc = JC.RecyclingPreconditioner(base, jvp!, n; kmax = kmax,
            kharvest = 2, form = form)
        pc.state.U = copy(U)
        pc.fresh = false
        JC.updatepreconditioner!(pc, zeros(n))
        return pc
    end
    # the preconditioned operator as a dense matrix, column by column
    function densepreconditioned(pc)
        M = zeros(n, n)
        z = zeros(n)
        for j in 1:n
            e = zeros(n); e[j] = 1.0
            JC.applypreconditioner!(z, pc, e)
            M[:, j] = Aref[]*z
        end
        return M
    end

    @testset "construction and the pair invariants" begin
        @test_throws ArgumentError JC.RecyclingPreconditioner(base, jvp!, n;
            form = :bogus)
        @test_throws DimensionMismatch JC.RecyclingPreconditioner(base, jvp!,
            n; state = JC.RecyclingState(zeros(n + 1)))
        pc1 = build(:adef1); pc2 = build(:adef2)
        @test JC.deflationform(pc1) == :adef1
        @test JC.deflationform(pc2) == :adef2
        @test JC.deflationform(base) == :none
        @test JC.deflationsize(pc1) == 2 && JC.deflationsize(pc2) == 2
        @test JC.candidatecount(pc1) == 2 && JC.deflationsize(base) == 0
        @test JC.deflationrebuilds(pc1) == 1
        # a build costs one product per candidate, and nothing else
        @test JC.deflationproducts(pc1) == 2 && JC.deflationproducts(pc2) == 2
        for pc in (pc1, pc2)
            @test isapprox(A*pc.X, pc.C; rtol = 1e-12)
            @test isapprox(pc.C'*pc.C, I; atol = 1e-12)
            # the candidates are kept orthonormal and span what they spanned,
            # and the corrections are their base solves
            @test isapprox(pc.state.U'*pc.state.U, I; atol = 1e-12)
            @test norm((I - Uexact*Uexact')*pc.state.U) < 1e-12
            @test norm((I - Xexact*Xexact')*pc.X) < 1e-10
        end
        # adef1 holds the folded block W and adef2 does not
        @test size(pc1.W, 2) == 2 && size(pc2.W, 2) == 0
        @test isapprox(pc1.W, pc1.X - P\pc1.C; rtol = 1e-12)
        @test pc1.fresh && pc2.fresh
    end

    @testset "a zero Galerkin pairing keeps its direction" begin
        # the review's operator: `J = [0 -e; 1 0]` with the candidate `e2`
        # has `Z'*J*Z = 0` while `norm(J*Z) = e`, so a Galerkin coarse
        # matrix would drop the one direction that matters; the image pair
        # keeps it, and the deflated operator has the eigenvalue one on it
        e = 1e-3
        J2 = [0.0 -e; 1.0 0.0]
        jvp2!(y, v) = mul!(y, J2, v)
        base2 = FixedP(lu(Matrix(1.0I, 2, 2)))
        for form in (:adef1, :adef2)
            pc = JC.RecyclingPreconditioner(base2, jvp2!, 2; kmax = 2,
                kharvest = 1, form = form)
            pc.state.U = [0.0; 1.0;;]
            pc.fresh = false
            JC.updatepreconditioner!(pc, zeros(2))
            @test JC.deflationsize(pc) == 1
            @test isapprox(J2*pc.X, pc.C; atol = 1e-15)
            M = zeros(2, 2)
            for j in 1:2
                ej = zeros(2); ej[j] = 1.0
                M[:, j] = J2*JC.applypreconditioner!(zeros(2), pc, ej)
            end
            @test isapprox(sort(real.(eigvals(M))), [0.0, 1.0]; atol = 1e-12)
            # and the application differs from the base's
            @test !isapprox(JC.applypreconditioner!(zeros(2), pc, [1.0, 0.0]),
                [1.0, 0.0]; atol = 1e-6)
        end
    end

    @testset "the two forms are spectrally equivalent with an exact subspace" begin
        pc1 = build(:adef1); pc2 = build(:adef2)
        M0 = A/P
        M1 = densepreconditioned(pc1)
        M2 = densepreconditioned(pc2)
        ev0 = sort(real.(eigvals(M0)))
        ev1 = sort(real.(eigvals(M1)))
        ev2 = sort(real.(eigvals(M2)))
        # the undeflated operator has the three planted eigenvalues
        @test isapprox(ev0[1:2], [0.01, 0.05]; atol = 1e-8)
        @test isapprox(ev0[end], 3.0; atol = 1e-8)
        # both deflated operators have moved the two near zero to one and
        # left the third alone: spectrum {1 (n-1 times), 3}
        for ev in (ev1, ev2)
            @test isapprox(ev[end], 3.0; atol = 1e-6)
            @test all(isapprox.(ev[1:end-1], 1.0; atol = 1e-6))
        end
        @test isapprox(ev1, ev2; atol = 1e-6)
        # the invariant subspaces differ: adef1 fixes range(C), adef2 fixes
        # the complementary projection
        @test M1*pc1.C ≈ pc1.C
        Pi = I - pc2.C*pc2.C'
        @test Pi*M2 ≈ Pi*M0
    end

    @testset "the fused product equals the standalone form and the operator" begin
        pc2 = build(:adef2)
        v = randn(n)
        zs = JC.applypreconditioner!(zeros(n), pc2, v)
        # the standalone form paid a product; the fused one does not
        @test JC.deflationproducts(pc2) == 3
        z = zeros(n); w = zeros(n)
        tp = JC.preconditionedproduct!(w, z, Aop, pc2, v)
        @test JC.deflationproducts(pc2) == 3
        @test tp isa Float64 && tp >= 0
        @test isapprox(z, zs; rtol = 1e-12)
        @test isapprox(w, A*z; rtol = 1e-12)
        # adef1 goes through the generic path, which is exact by construction
        pc1 = build(:adef1)
        zs1 = JC.applypreconditioner!(zeros(n), pc1, v)
        JC.preconditionedproduct!(w, z, Aop, pc1, v)
        @test z == zs1
        @test isapprox(w, A*z; rtol = 1e-12)
        # a closure preconditioner still works through the same entry point
        Mop!(zz, vv) = ldiv!(zz, base.F, vv)
        JC.preconditionedproduct!(w, z, Aop, Mop!, v)
        @test isapprox(z, P\v; rtol = 1e-12)
        @test isapprox(w, A*z; rtol = 1e-12)
    end

    @testset "a moved operator is refreshed lazily and the product stays exact" begin
        pc2 = build(:adef2)
        A2 = A + 0.2*randn(n, n)
        Aref[] = A2
        try
            v = randn(n)
            z = zeros(n); w = zeros(n)
            # nothing told the preconditioner, so C is stale and the fused
            # product would not be the operator's; the standalone form is
            # exact regardless
            zs = JC.applypreconditioner!(zeros(n), pc2, v)
            JC.pointmoved!(pc2)
            @test !pc2.fresh
            JC.preconditionedproduct!(w, z, Aop, pc2, v)
            @test pc2.fresh
            @test isapprox(A2*pc2.X, pc2.C; rtol = 1e-12)
            @test isapprox(w, A2*z; rtol = 1e-12)
            # the refreshed pair differs from the stale one, so the
            # standalone result taken before the refresh differs from the
            # fused one, while both are exact preconditioner applications
            @test !isapprox(z, zs; rtol = 1e-8)
            @test isapprox(z, JC.applypreconditioner!(zeros(n), pc2, v);
                rtol = 1e-12)
            # adef1 is refreshed the same way
            pc1 = build(:adef1)
            JC.pointmoved!(pc1)
            JC.preconditionedproduct!(w, z, Aop, pc1, v)
            @test isapprox(A2*pc1.X, pc1.C; rtol = 1e-12)
            @test isapprox(w, A2*z; rtol = 1e-12)
            # and forwards through the wrappers
            sp = JC.SizedPreconditioner(pc2, n)
            JC.pointmoved!(sp)
            @test !pc2.fresh
            @test JC.deflationform(sp) == :adef2
            @test JC.deflationsize(sp) == 2
            @test JC.candidatecount(sp) == 2
            @test JC.deflationproducts(sp) == JC.deflationproducts(pc2)
        finally
            Aref[] = A
        end
    end

    @testset "the build resolves the rank and trims by quality" begin
        # two candidates whose images are dependent: the pair has rank one
        # and one candidate is dropped rather than inverted
        Udup = hcat(Uexact[:, 1], Uexact[:, 1])
        pc = JC.RecyclingPreconditioner(base, jvp!, n; kmax = 4,
            kharvest = 2, form = :adef2)
        pc.state.U = Udup
        pc.fresh = false
        JC.updatepreconditioner!(pc, zeros(n))
        @test JC.deflationsize(pc) == 1
        @test JC.candidatecount(pc) == 1
        @test isapprox(A*pc.X, pc.C; rtol = 1e-12)
        # a candidate the operator annihilates falls below the resolution
        # of the build and is dropped, leaving the resolved one active
        Anull = copy(A); Anull[:, 1] .= 0.0
        Aref[] = Anull
        try
            pc = JC.RecyclingPreconditioner(base, jvp!, n; kmax = 4,
                kharvest = 2, form = :adef1)
            # the candidate whose correction is e1 is u = P*e1
            u1 = P[:, 1]/norm(P[:, 1])
            pc.state.U = hcat(u1, Matrix(qr(hcat(u1, Uexact[:, 2])).Q)[:, 2])
            pc.fresh = false
            JC.updatepreconditioner!(pc, zeros(n))
            @test JC.deflationsize(pc) == 1
            @test JC.candidatecount(pc) == 1
            @test isapprox(Anull*pc.X, pc.C; rtol = 1e-12)
        finally
            Aref[] = A
        end
        # more candidates than kmax: the ones the preconditioned operator
        # shrinks most stay. The trim keeps the smallest singular directions
        # of A*inv(P) on the candidate span, which are its two near-zero
        # eigendirections up to the non-orthogonality of the three
        # eigenvectors, a leakage of order 1e-3 here
        U3 = Matrix(qr(P*Q3).Q)[:, 1:3]
        pc = build(:adef2; U = U3, kmax = 2)
        @test JC.deflationsize(pc) == 2 && JC.candidatecount(pc) == 2
        @test norm((I - Uexact*Uexact')*pc.state.U) < 1e-2
        @test norm((I - U3*U3')*pc.state.U) < 1e-10
    end

    @testset "appending candidates never leaves their span" begin
        x = randn(n)
        # duplicates: one direction, not two
        B = JC._orthappend(zeros(n, 0), hcat(x, x))
        @test size(B, 2) == 1
        @test norm(B - x/norm(x)*sign(dot(B, x))) < 1e-12
        # a near duplicate is one direction too
        B = JC._orthappend(zeros(n, 0), hcat(x, x*(1 + 1e-9)))
        @test size(B, 2) == 1
        # against an existing basis, a column already in its span is dropped
        # and an independent one is added orthogonal to it
        U = Matrix(qr(randn(n, 3)).Q)[:, 1:3]
        y = randn(n)
        B = JC._orthappend(U, hcat(U*randn(3), y))
        @test size(B, 2) == 4
        @test isapprox(B'*B, I; atol = 1e-12)
        @test norm((I - hcat(U, y)*pinv(hcat(U, y)))*B) < 1e-10
        # nothing to add
        @test JC._orthappend(U, U*randn(3, 2)) === U
    end

    @testset "gmres! with each form: iteration counts from the structure" begin
        # A*inv(P) = I + rank 3 with four distinct eigenvalues has minimal
        # polynomial degree four, so the base alone needs four Arnoldi
        # steps. Both deflated operators have spectrum {1, 3}, but only one
        # of them is diagonalizable: with `range(C)` the invariant subspace
        # of the two deflated eigenvalues,
        #
        #   adef2:  J*B2 - I = Pi*(J*B - I)      rank one (Pi kills c1, c2)
        #   adef1:  J*B1 - I = (J*B - I)*Pi      rank three, one nonzero eigenvalue
        #
        # where Pi = I - C*C'. So adef2 has minimal polynomial
        # (t - 1)(t - 3) and converges in two steps, while adef1 carries a
        # nilpotent part at eigenvalue one and needs up to four. This is the
        # diagonalizability asymmetry between the two forms that the
        # literature attributes their different robustness to.
        function solve(M)
            ws = JC.GMRESWorkspace(n, 20)
            x = zeros(n)
            out = JC.gmres!(x, jvp!, b, ws; Mop! = M, rtol = 1e-11,
                maxrestarts = 5)
            return x, out
        end
        x0, o0 = solve(base)
        @test o0.converged && o0.iterations <= 4
        @test o0.products == o0.iterations + o0.cycles
        @test isapprox(x0, A\b; rtol = 1e-8)
        bound = Dict(:adef1 => 4, :adef2 => 2)
        for form in (:adef1, :adef2)
            pc = build(form)
            x, o = solve(pc)
            @test o.converged
            @test o.iterations <= bound[form]
            @test isapprox(x, A\b; rtol = 1e-8)
        end
        # the same linear maps through a closure, which cannot fuse: the
        # counts are unchanged, and adef2 pays one extra product per step
        for form in (:adef1, :adef2)
            pc = build(form)
            x, o = solve((z, v) -> JC.applypreconditioner!(z, pc, v))
            @test o.converged && o.iterations <= bound[form]
            @test isapprox(x, A\b; rtol = 1e-8)
        end
    end

    @testset "exact product accounting across restarts" begin
        # a short restart forces several cycles; every product is counted,
        # including the one adef2 pays at each restart correction
        for form in (:adef1, :adef2)
            calls = Ref(0)
            jvpc!(y, v) = (calls[] += 1; mul!(y, A, v))
            pc = JC.RecyclingPreconditioner(base, jvpc!, n; kmax = 4,
                kharvest = 2, form = form)
            pc.state.U = Matrix(qr(randn(n, 2)).Q)[:, 1:2]
            pc.fresh = false
            JC.updatepreconditioner!(pc, zeros(n))
            built = calls[]
            @test built == 2
            ws = JC.GMRESWorkspace(n, 3)
            x = zeros(n)
            out = JC.gmres!(x, jvpc!, b, ws; Mop! = pc, rtol = 1e-8,
                maxrestarts = 30)
            @test out.cycles > 1
            expected = out.iterations + out.cycles +
                (form === :adef2 ? out.cycles : 0)
            @test calls[] - built == expected
            @test out.products + (JC.deflationproducts(pc) - built) == expected
        end
    end

    @testset "nlsolvekrylov! harvests and converges with either form" begin
        # A mildly nonlinear system, F(x) = (L + R)*x + 0.3*x.^3 - c, whose
        # base preconditioner is a fixed factorization of L alone: the rank
        # three R is what the Krylov space has to resolve at every step, so
        # the harvest has something to find and the deflation something to
        # do, while the cubic term keeps Newton's path short.
        Random.seed!(20260903)
        m = 40
        L = randn(m, m) + 8I
        R = 3.0*randn(m, 3)*randn(3, m)/m
        c = 0.5*randn(m)
        xpt = zeros(m)
        fjn!(F, J, x) = begin
            copyto!(xpt, x)
            isnothing(F) || (F .= L*x .+ R*x .+ 0.3 .* x.^3 .- c)
            nothing
        end
        jvpn!(y, v) = (mul!(y, L, v); mul!(y, R, v, true, true);
            y .+= 0.9 .* xpt.^2 .* v; y)
        basen = FixedP(lu(L))
        # reference root from a dense Newton with the exact Jacobian
        xref = zeros(m)
        for _ in 1:50
            Fr = L*xref .+ R*xref .+ 0.3 .* xref.^3 .- c
            xref .-= (L + R + Diagonal(0.9 .* xref.^2))\Fr
        end
        @test norm(L*xref .+ R*xref .+ 0.3 .* xref.^3 .- c) < 1e-12
        for form in (:adef1, :adef2)
            pc = JC.RecyclingPreconditioner(basen, jvpn!, m; kmax = 6,
                kharvest = 2, form = form)
            x = zeros(m); F = zeros(m)
            info = JC.nlsolvekrylov!(fjn!, jvpn!, F, x, pc; ftol = 1e-10,
                krylovescalate = typemax(Int))
            @test info.converged
            @test isapprox(x, xref; rtol = 1e-7)
            @test any(k -> k.deflationsize > 0, info.krylov)
            @test JC.deflationform(pc) == form
            # the records carry the products: the solve's own and the
            # wrapper's running count
            @test all(k -> k.products >= k.iterations + k.cycles, info.krylov)
            @test info.krylov[end].deflationproducts == JC.deflationproducts(pc)
            @test JC.deflationproducts(pc) > 0
        end
        # a cycle shorter than kharvest is not harvested: its singular
        # directions are not a measurement of the operator
        pc = JC.RecyclingPreconditioner(basen, jvpn!, m; kmax = 6,
            kharvest = 4, form = :adef2)
        ws = JC.GMRESWorkspace(m, 10)
        x = zeros(m)
        out = JC.gmres!(x, jvpn!, c, ws; Mop! = pc, rtol = 0.5)
        @test out.iterations >= 1 && out.iterations < 4
        JC.harvest!(pc, ws, out)
        @test JC.candidatecount(pc) == 0
        out = JC.gmres!(x, jvpn!, c, ws; Mop! = pc, rtol = 1e-6)
        @test out.iterations >= 4
        JC.harvest!(pc, ws, out)
        @test JC.candidatecount(pc) == 4
        # the same with the base frozen across the whole Newton path, which
        # is the regime where the pair is refreshed lazily
        for form in (:adef1, :adef2)
            pc = JC.RecyclingPreconditioner(basen, jvpn!, m; kmax = 6,
                kharvest = 2, form = form)
            x = zeros(m); F = zeros(m)
            info = JC.nlsolvekrylov!(fjn!, jvpn!, F, x, pc; ftol = 1e-10,
                krylovescalate = typemax(Int),
                krylovrefreshiterations = typemax(Int),
                krylovrefreshrate = 1.0)
            @test info.converged
            @test isapprox(x, xref; rtol = 1e-7)
            @test info.krylov[end].deflationrebuilds > 0
        end
        # a state handed in is what the solve starts from, and it is the
        # object the solve mutates
        st = JC.RecyclingState(zeros(m))
        pc = JC.RecyclingPreconditioner(basen, jvpn!, m; kmax = 6,
            kharvest = 2, form = :adef2, state = st)
        x = zeros(m); F = zeros(m)
        JC.nlsolvekrylov!(fjn!, jvpn!, F, x, pc; ftol = 1e-10,
            krylovescalate = typemax(Int))
        @test pc.state === st && size(st.U, 2) > 0
        pc2 = JC.RecyclingPreconditioner(basen, jvpn!, m; kmax = 6,
            kharvest = 2, form = :adef2, state = copy(st))
        @test !pc2.fresh
        x = zeros(m); F = zeros(m)
        info = JC.nlsolvekrylov!(fjn!, jvpn!, F, x, pc2; ftol = 1e-10,
            krylovescalate = typemax(Int))
        @test info.converged && info.krylov[1].deflationsize > 0
    end
end
