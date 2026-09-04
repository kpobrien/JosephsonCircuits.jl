using JosephsonCircuits
using LinearAlgebra
using Test
using Random

const JC = JosephsonCircuits

# A dense matrix dressed as a base preconditioner, so the invariants can be
# checked against an operator whose exact inverse is available.
struct DensePC{T} <: JC.AbstractPreconditioner
    B::Matrix{T}
    exact::Bool
end
DensePC(B) = DensePC(B, false)
JC.applypreconditioner!(z, pc::DensePC, r) = (mul!(z, pc.B, r); z)
JC.updatepreconditioner!(pc::DensePC, x) = pc
JC.isexactpreconditioner(pc::DensePC) = pc.exact
JC.escalatepreconditioner!(::DensePC) = false

# `J` with a few directions the base `B0` inverts badly: `B0` is the exact
# inverse on the complement and is wrong by a large factor on `nbad` of the
# singular directions, which is the low-rank-defect situation the
# preconditioner is built for.
function defectsystem(n, nbad; seed = 1, factor = 50.0)
    rng = MersenneTwister(seed)
    Q = Matrix(qr(randn(rng, n, n)).Q)
    d = collect(range(1.0, 2.0; length = n))
    J = Q*Diagonal(d)*Q'
    dinv = 1 ./ d
    # spoil the base's inverse on the leading `nbad` directions
    for i in 1:nbad
        dinv[i] *= factor
    end
    B0 = Q*Diagonal(dinv)*Q'
    return J, B0, Q
end

@testset verbose = true "floquetdeflation" begin

    @testset "construction and argument checking" begin
        J, B0, _ = defectsystem(20, 2)
        b = zeros(20)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, b)
        @test JC.deflationform(pc) === :floquet
        @test JC.deflationsize(pc) == 0
        @test JC.candidatecount(pc) == 0
        @test JC.deflationrebuilds(pc) == 0
        # the integer constructor agrees with the vector one
        pc2 = JC.FloquetPreconditioner(DensePC(B0), jvp!, 20)
        @test JC.deflationsize(pc2) == 0

        @test_throws ArgumentError JC.FloquetPreconditioner(DensePC(B0), jvp!, b;
            kmax = 0)
        @test_throws ArgumentError JC.FloquetPreconditioner(DensePC(B0), jvp!, b;
            kharvest = 0, nritz = 0)
        @test_throws ArgumentError JC.FloquetPreconditioner(DensePC(B0), jvp!, b;
            kmax = 10, kcandidate = 4)
        @test_throws ArgumentError JC.FloquetPreconditioner(DensePC(B0), jvp!, b;
            ranktol = 0.0)
        @test_throws ArgumentError JC.FloquetPreconditioner(DensePC(B0), jvp!, b;
            benefittol = -1.0)
        @test_throws DimensionMismatch JC.seeddeflation!(pc, randn(7, 2))
        # a preconditioner which does not deflate ignores a seed
        @test JC.candidatecount(DensePC(B0)) == 0
        @test_throws DimensionMismatch JC.FloquetPreconditioner(DensePC(B0), jvp!, b;
            state = JC.FloquetState(zeros(7)))
        @test JC.seeddeflation!(DensePC(B0), randn(20, 2)) isa DensePC
    end

    # 33.1: the active residual-image subspace is an exact eigenspace of the
    # right preconditioned operator at eigenvalue one.
    @testset "unit eigenspace invariant J*B*C = C" begin
        n, nbad = 40, 3
        J, B0, Q = defectsystem(n, nbad)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        # seed the bad directions themselves, plus noise
        JC.seeddeflation!(pc, hcat(Q[:, 1:nbad], randn(MersenneTwister(3), n, 2));
            source = :test)
        JC._rebuildfloquet!(pc)
        @test JC.deflationsize(pc) >= nbad

        # C is orthonormal and J*X = C, by construction
        @test pc.C'*pc.C ≈ I atol = 1e-10
        @test J*pc.X ≈ pc.C atol = 1e-9

        for j in axes(pc.C, 2)
            c = Vector(pc.C[:, j])
            z = similar(c)
            JC.applypreconditioner!(z, pc, c)
            @test J*z ≈ c atol = 1e-8
        end
    end

    # 33.2: duplicated candidates must not inflate the active rank.
    @testset "rank-deficient candidates are compressed away" begin
        n = 30
        J, B0, Q = defectsystem(n, 2)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        x = Q[:, 1:2]
        # the same two directions offered five times over, in different
        # scalings and linear combinations
        Xdup = hcat(x, 2.5*x, -x, x*[1.0 1.0; 1.0 -1.0], 1e-3*x)
        JC.seeddeflation!(pc, Xdup; source = :test)
        @test JC.candidatecount(pc) == 10
        JC._rebuildfloquet!(pc)
        @test JC.deflationsize(pc) == 2
        @test J*pc.X ≈ pc.C atol = 1e-9
    end

    # 33.3: a direction the base already inverts must be discarded.
    @testset "base-exact directions are filtered out" begin
        n, nbad = 30, 2
        J, B0, Q = defectsystem(n, nbad)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        # columns nbad+1 onward satisfy B0*J*x = x exactly; only the first
        # nbad are missing from the base
        JC.seeddeflation!(pc, Q[:, 1:6]; source = :test)
        JC._rebuildfloquet!(pc)
        @test JC.deflationsize(pc) == nbad
        @test all(JC.correctionstrengths(pc) .> 1.0)

        # every candidate handled by the base leaves an empty active set
        pc2 = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        JC.seeddeflation!(pc2, Q[:, nbad+1:nbad+4]; source = :test)
        JC._rebuildfloquet!(pc2)
        @test JC.deflationsize(pc2) == 0
        # ... and against an exact base nothing is retained at all
        pcx = JC.FloquetPreconditioner(DensePC(inv(J), true), jvp!, zeros(n))
        JC.seeddeflation!(pcx, Q[:, 1:4]; source = :test)
        JC.updatepreconditioner!(pcx, zeros(n))
        @test JC.deflationsize(pcx) == 0
    end

    @testset "kmax caps the active rank" begin
        n = 40
        J, B0, Q = defectsystem(n, 8)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n); kmax = 3)
        JC.seeddeflation!(pc, Q[:, 1:8]; source = :test)
        JC._rebuildfloquet!(pc)
        @test JC.deflationsize(pc) == 3
        # the compression is orthonormal, so the invariants survive it
        @test pc.C'*pc.C ≈ I atol = 1e-10
        @test J*pc.X ≈ pc.C atol = 1e-9
    end

    @testset "the candidate bank stays bounded" begin
        n = 25
        J, B0, _ = defectsystem(n, 2)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n);
            kmax = 4, kcandidate = 10)
        rng = MersenneTwister(11)
        for _ in 1:20
            JC.seeddeflation!(pc, randn(rng, n, 3); source = :test)
            @test JC.candidatecount(pc) <= 10
            @test length(pc.state.source) == JC.candidatecount(pc)
        end
        # a zero column carries no direction and is not banked
        before = JC.candidatecount(pc)
        JC.seeddeflation!(pc, zeros(n, 2); source = :test)
        @test JC.candidatecount(pc) == before
    end

    # 33.4: harmonic Ritz values must find the small eigenvalues.
    @testset "harmonic Ritz extraction near zero" begin
        # a full Arnoldi factorization of a matrix with two small
        # eigenvalues, built here so H is exactly the projection
        n = 24
        rng = MersenneTwister(5)
        Qm = Matrix(qr(randn(rng, n, n)).Q)
        d = collect(range(1.0, 3.0; length = n))
        d[1] = 1e-3; d[2] = 2e-3
        A = Qm*Diagonal(d)*Qm'
        m = 16
        V = zeros(n, m+1); H = zeros(n, m)
        b = randn(rng, n); V[:, 1] = b/norm(b)
        for j in 1:m
            w = A*V[:, j]
            for i in 1:j
                H[i, j] = dot(V[:, i], w); w -= H[i, j]*V[:, i]
            end
            for i in 1:j   # reorthogonalize
                c = dot(V[:, i], w); H[i, j] += c; w -= c*V[:, i]
            end
            H[j+1, j] = norm(w); V[:, j+1] = w/H[j+1, j]
        end
        Hbar = H[1:m+1, 1:m]
        Y = JC.harmonicritznearzero(Hbar, 2)
        @test size(Y, 2) >= 2
        # the harmonic Ritz directions should live in the span of the two
        # small eigenvectors
        U = Qm[:, 1:2]
        for j in axes(Y, 2)
            q = V[:, 1:m]*Y[:, j]
            q ./= norm(q)
            @test norm(U'*q) > 0.9
        end
        # a singular projected block gives an empty result rather than an
        # error, and the caller falls back on the singular directions
        @test size(JC.harmonicritznearzero(zeros(5, 4), 2), 2) == 0
        @test size(JC.harmonicritznearzero(Hbar, 0), 2) == 0
    end

    # 33.5: on a strongly nonnormal matrix the smallest singular direction
    # and the smallest eigendirection differ; both families are harvested
    # and the residual image picks a usable subspace.
    @testset "nonnormal harvesting keeps both families" begin
        n = 30
        rng = MersenneTwister(7)
        # upper triangular with a large off-diagonal: eigenvectors are far
        # from orthogonal, so eigen and singular information disagree
        A = triu(randn(rng, n, n), 1)*8.0 + Diagonal(range(0.05, 2.0; length = n))
        B0 = Matrix{Float64}(I, n, n)
        jvp!(y, v) = (mul!(y, A, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n);
            kharvest = 2, nritz = 2, kmax = 8)
        # a real Arnoldi factorization of A on a random start
        m = 12
        V = zeros(n, m+1); H = zeros(n, m)
        b = randn(rng, n); V[:, 1] = b/norm(b)
        for j in 1:m
            w = A*V[:, j]
            for i in 1:j
                H[i, j] = dot(V[:, i], w); w -= H[i, j]*V[:, i]
            end
            H[j+1, j] = norm(w); V[:, j+1] = w/H[j+1, j]
        end
        JC._harvestfloquet!(pc, V[:, 1:m], H[1:m+1, 1:m])
        # both families contributed
        @test JC.candidatecount(pc) >= 3
        @test all(s -> s === :gmres_floquet, pc.state.source)
        JC._rebuildfloquet!(pc)
        @test JC.deflationsize(pc) >= 1
        @test A*pc.X ≈ pc.C atol = 1e-7
        @test pc.C'*pc.C ≈ I atol = 1e-9
    end

    @testset "harvesting does not rebuild under a running solve" begin
        n = 30
        J, B0, Q = defectsystem(n, 2)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        JC.seeddeflation!(pc, Q[:, 1:2]; source = :test)
        JC._rebuildfloquet!(pc)
        rebuilds = JC.deflationrebuilds(pc)
        active = copy(pc.W)
        m = 8
        rng = MersenneTwister(13)
        V = Matrix(qr(randn(rng, n, m)).Q)
        H = triu(randn(rng, m+1, m), -1)
        JC._harvestfloquet!(pc, V, H)
        # the bank grew but the applied operator did not change
        @test JC.candidatecount(pc) > 2
        @test JC.deflationrebuilds(pc) == rebuilds
        @test pc.W == active
    end

    # the application is exactly B0 + W*C'
    @testset "application matches the closed form" begin
        n = 30
        J, B0, Q = defectsystem(n, 3)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        JC.seeddeflation!(pc, Q[:, 1:3]; source = :test)
        JC._rebuildfloquet!(pc)
        rng = MersenneTwister(17)
        r = randn(rng, n)
        z = similar(r)
        JC.applypreconditioner!(z, pc, r)
        @test z ≈ B0*r + pc.W*(pc.C'*r) atol = 1e-10
        # and it beats the base on the deflated directions
        @test norm(J*z - r) < norm(J*(B0*r) - r)
    end

    @testset "a moved point is refreshed lazily" begin
        n = 30
        J, B0, Q = defectsystem(n, 2)
        Jmoved = J + 0.01*Q*Diagonal(randn(MersenneTwister(19), n))*Q'
        moved = Ref(false)
        jvp!(y, v) = (mul!(y, moved[] ? Jmoved : J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        JC.seeddeflation!(pc, Q[:, 1:2]; source = :test)
        JC._rebuildfloquet!(pc)
        @test J*pc.X ≈ pc.C atol = 1e-9
        rebuilds = JC.deflationrebuilds(pc)

        moved[] = true
        JC.pointmoved!(pc)
        # stale until something asks for an application
        @test JC.deflationrebuilds(pc) == rebuilds
        z = zeros(n)
        JC.applypreconditioner!(z, pc, randn(MersenneTwister(23), n))
        @test JC.deflationrebuilds(pc) == rebuilds + 1
        # the identity now holds against the *new* Jacobian
        @test Jmoved*pc.X ≈ pc.C atol = 1e-9
        for j in axes(pc.C, 2)
            c = Vector(pc.C[:, j]); zz = similar(c)
            JC.applypreconditioner!(zz, pc, c)
            @test Jmoved*zz ≈ c atol = 1e-8
        end
    end

    @testset "gmres! converges faster with the correction" begin
        n, nbad = 60, 4
        J, B0, Q = defectsystem(n, nbad; seed = 29, factor = 200.0)
        jvp!(y, v) = (mul!(y, J, v); y)
        rng = MersenneTwister(31)
        b = randn(rng, n)

        base = DensePC(B0)
        ws = JC.GMRESWorkspace(b, 40)
        x = zeros(n)
        out0 = JC.gmres!(x, jvp!, b, ws; Mop! = base, rtol = 1e-10,
            maxrestarts = 4)
        @test out0.converged

        pc = JC.FloquetPreconditioner(base, jvp!, zeros(n); kmax = 8)
        JC.seeddeflation!(pc, Q[:, 1:nbad]; source = :test)
        JC._rebuildfloquet!(pc)
        ws2 = JC.GMRESWorkspace(b, 40)
        x2 = zeros(n)
        out1 = JC.gmres!(x2, jvp!, b, ws2; Mop! = pc, rtol = 1e-10,
            maxrestarts = 4)
        @test out1.converged
        @test out1.iterations < out0.iterations
        @test norm(J*x2 - b) <= 1e-8*norm(b)
    end

    # 33.7: a sequence of slowly changing operators, the subspace carried
    # between them, must beat a cold start at the second point.
    @testset "sweep reuse across a changing operator" begin
        n, nbad = 60, 3
        J0, B0, Q = defectsystem(n, nbad; seed = 37, factor = 150.0)
        # a nearby point: the difficult directions move slowly
        P = Q*Diagonal(1 .+ 0.02*randn(MersenneTwister(41), n))*Q'
        J1 = J0*P
        rng = MersenneTwister(43)
        b = randn(rng, n)
        base = DensePC(B0)

        # cold start at the second operator
        jvp1!(y, v) = (mul!(y, J1, v); y)
        cold = JC.FloquetPreconditioner(base, jvp1!, zeros(n); kmax = 8)
        wsc = JC.GMRESWorkspace(b, 40); xc = zeros(n)
        outcold = JC.gmres!(xc, jvp1!, b, wsc; Mop! = cold, rtol = 1e-10,
            maxrestarts = 4)

        # warm: harvest at the first operator, carry the physical vectors
        Jcur = Ref(J0)
        jvp!(y, v) = (mul!(y, Jcur[], v); y)
        warm = JC.FloquetPreconditioner(base, jvp!, zeros(n); kmax = 8,
            kharvest = 4, nritz = 2)
        ws = JC.GMRESWorkspace(b, 40); x = zeros(n)
        out0 = JC.gmres!(x, jvp!, b, ws; Mop! = warm, rtol = 1e-10,
            maxrestarts = 4)
        JC.harvest!(warm, ws, out0)
        @test JC.candidatecount(warm) > 0

        # move to the second operator; the bank is rebuilt against it
        Jcur[] = J1
        JC.pointmoved!(warm)
        ws2 = JC.GMRESWorkspace(b, 40); x2 = zeros(n)
        outwarm = JC.gmres!(x2, jvp!, b, ws2; Mop! = warm, rtol = 1e-10,
            maxrestarts = 4)
        @test outwarm.converged
        @test JC.deflationsize(warm) > 0
        @test J1*warm.X ≈ warm.C atol = 1e-8
        # the inherited directions pay for themselves
        @test outwarm.iterations < outcold.iterations
        # and stale directions are re-evaluated rather than trusted
        @test JC.deflationrebuilds(warm) >= 1
    end

    # 33.6: the difficult direction can appear in an early full cycle while
    # the cycle left in the workspace at the end is short. Harvesting per
    # cycle must not lose it.
    @testset "restart-cycle harvesting keeps early cycles" begin
        n, nbad = 60, 6
        # `defectsystem` leaves `J*B0` with only two distinct eigenvalues,
        # on which GMRES converges in two steps and never restarts. Here the
        # base is imperfect everywhere, mildly, and badly wrong on `nbad`
        # directions, so the preconditioned spectrum is spread and the solve
        # takes several cycles.
        rng = MersenneTwister(53)
        Q = Matrix(qr(randn(rng, n, n)).Q)
        d = collect(range(1.0, 2.0; length = n))
        J = Q*Diagonal(d)*Q'
        g = 1 .+ 0.9*rand(rng, n)
        g[1:nbad] .= range(100.0, 500.0; length = nbad)
        B0 = Q*Diagonal(g ./ d)*Q'
        jvp!(y, v) = (mul!(y, J, v); y)
        b = randn(MersenneTwister(59), n)
        base = DensePC(B0)
        # a restart cycle short enough that the solve takes several of
        # them, which is the situation the per-cycle harvest exists for
        m = 4

        percycle = JC.FloquetPreconditioner(base, jvp!, zeros(n);
            kmax = 8, kharvest = 2, nritz = 2, kcandidate = 64)
        @test JC.usescycleharvest(percycle)
        ws = JC.GMRESWorkspace(b, m); x = zeros(n)
        before = copy(percycle.W)
        out = JC.gmres!(x, jvp!, b, ws; Mop! = percycle, rtol = 1e-10,
            maxrestarts = 10,
            oncycle = (w, j) -> JC.harvestcycle!(percycle, w, j))
        @test out.cycles >= 2
        # the preconditioner did not move under the running solve
        @test percycle.W == before
        @test JC.deflationrebuilds(percycle) == 0

        # against harvesting only the cycle left in the workspace
        finalonly = JC.FloquetPreconditioner(base, jvp!, zeros(n);
            kmax = 8, kharvest = 2, nritz = 2, kcandidate = 64,
            cycleharvest = false)
        @test !JC.usescycleharvest(finalonly)
        ws2 = JC.GMRESWorkspace(b, m); x2 = zeros(n)
        out2 = JC.gmres!(x2, jvp!, b, ws2; Mop! = finalonly, rtol = 1e-10,
            maxrestarts = 10)
        JC.harvest!(finalonly, ws2, out2)
        @test JC.candidatecount(percycle) >
            JC.candidatecount(finalonly)

        # and what the extra cycles found survives the rebuild
        JC.pointmoved!(percycle)
        JC._refreshfloquet!(percycle)
        @test JC.deflationsize(percycle) >= 1
        @test J*percycle.X ≈ percycle.C atol = 1e-7
    end

    @testset "an external seed takes effect at the next solve; a harvest waits" begin
        n = 40
        J, B0, Q = defectsystem(n, 3)
        jvp!(y, v) = (mul!(y, J, v); y)
        pc = JC.FloquetPreconditioner(DensePC(B0), jvp!, zeros(n))
        JC.seeddeflation!(pc, Q[:, 1:2]; source = :test)
        JC._rebuildfloquet!(pc)
        @test JC.deflationsize(pc) == 2
        W = copy(pc.W); rebuilds = JC.deflationrebuilds(pc)
        # a candidate banked by the harvest, mid-solve, leaves the applied
        # operator alone until the point moves
        JC._bankcandidates!(pc, Q[:, 3:3]; source = :test)
        z = zeros(n)
        JC.applypreconditioner!(z, pc, randn(MersenneTwister(67), n))
        @test pc.W == W
        @test JC.deflationrebuilds(pc) == rebuilds
        # an external seed between solves is active at the next application
        JC.seeddeflation!(pc, Q[:, 3:3]; source = :external)
        JC.applypreconditioner!(z, pc, randn(MersenneTwister(71), n))
        @test JC.deflationrebuilds(pc) == rebuilds + 1
        @test JC.deflationsize(pc) == 3
        c = J*Q[:, 3]; c ./= norm(c)
        @test norm(pc.C*(pc.C'*c) - c) < 1e-8
        # a seed is banked unfiltered and marks the blocks stale; the
        # rebuild is what removes what it duplicates
        JC.seeddeflation!(pc, Q[:, 1:1]; source = :external)
        @test JC.candidatecount(pc) == 4
        @test !pc.fresh
        JC.applypreconditioner!(z, pc, randn(MersenneTwister(73), n))
        @test JC.deflationsize(pc) == 3
    end

    @testset "harmonic Ritz through the production path reads the Arnoldi Hessenberg" begin
        # `gmres!` triangularizes `ws.H` in place with its Givens rotations;
        # the harmonic Ritz pencil needs the Arnoldi relation, which only
        # `ws.Harnoldi` still satisfies after a cycle
        n = 24
        rng = MersenneTwister(5)
        Qm = Matrix(qr(randn(rng, n, n)).Q)
        d = collect(range(1.0, 3.0; length = n))
        d[1] = 1e-3; d[2] = 2e-3
        A = Qm*Diagonal(d)*Qm'
        jvp!(y, v) = (mul!(y, A, v); y)
        b = randn(rng, n)
        m = 16
        ws = JC.GMRESWorkspace(b, m)
        x = zeros(n)
        out = JC.gmres!(x, jvp!, b, ws; rtol = 1e-14, maxrestarts = 1)
        j = out.iterations
        @test j == m
        Hr = ws.Harnoldi[1:j+1, 1:j]
        # the Arnoldi relation holds for the raw matrix ...
        @test A*ws.V[:, 1:j] ≈ ws.V[:, 1:j+1]*Hr atol = 1e-10
        @test all(!iszero, [Hr[i+1, i] for i in 1:j])
        # ... and not for the rotated one, whose subdiagonal is gone
        @test all(iszero, [ws.H[i+1, i] for i in 1:j])
        @test svdvals(Hr) ≈ svdvals(ws.H[1:j+1, 1:j])
        Y = JC.harmonicritznearzero(Hr, 2)
        @test size(Y, 2) >= 2
        U = Qm[:, 1:2]
        for c in axes(Y, 2)
            q = ws.V[:, 1:j]*Y[:, c]
            q ./= norm(q)
            @test norm(U'*q) > 0.9
        end
        # and the harvest, given the workspace, finds the same subspace
        pc = JC.FloquetPreconditioner(DensePC(Matrix{Float64}(I, n, n)), jvp!,
            zeros(n); kharvest = 0, nritz = 2)
        JC.harvest!(pc, ws, out)
        @test JC.candidatecount(pc) >= 2
        Xc = pc.state.X
        # the candidates are inv(P)*u = u here, so they should overlap the
        # small eigenspace
        for c in axes(Xc, 2)
            q = Xc[:, c]/norm(Xc[:, c])
            @test norm(U'*q) > 0.9
        end
    end

    # the whole path, on a driven Josephson chain
    function chain(Ncell, Lj)
        c = Tuple{String,String,String,Any}[]
        push!(c, ("P1", "1", "0", 1)); push!(c, ("R1", "1", "0", 50.0))
        for i in 1:Ncell
            push!(c, ("Lj$(i)", "$(i)", "$(i+1)", Lj))
            push!(c, ("C$(i)", "$(i)", "0", 40e-15))
        end
        push!(c, ("C$(Ncell+1)", "$(Ncell+1)", "0", 40e-15))
        push!(c, ("R2", "$(Ncell+1)", "0", 50.0))
        return c
    end

    @testset "the floquet form converges and harvests on a circuit" begin
        w = (2*pi*8e9,); Nh = (6,)
        src = [(mode = (1,), port = 1, current = 2.5e-6)]
        circ = chain(12, 100e-12)
        noesc = (; krylovescalate = typemax(Int))

        ref = JC.hbnlsolve(w, Nh, src, circ, Dict(); method = :newton,
            keyedarrays = false)
        sol = JC.hbnlsolve(w, Nh, src, circ, Dict(); method = :newtonkrylov,
            krylovrecycle = 8, krylovharvest = 3,
            krylovdeflationform = :floquet, krylovkwargs = noesc,
            keyedarrays = false)
        @test sol.solverinfo.converged
        st = sol.solverinfo.stages[1]
        @test !isempty(st.krylov)
        # something was harvested and built into the preconditioner
        @test any(k -> k.deflationrebuilds >= 1, st.krylov)
        @test any(k -> k.deflationsize >= 1, st.krylov)
        # and the answer is the one the direct solve gives
        @test isapprox(sol.nodeflux, ref.nodeflux; rtol = 1e-6,
            atol = 1e-12*maximum(abs, ref.nodeflux))

        @test_throws ArgumentError JC.hbnlsolve(w, Nh, src, circ, Dict();
            method = :newtonkrylov, krylovrecycle = 8,
            krylovdeflationform = :nosuchform, keyedarrays = false)

        # and the same through hbsolve, which forwards the three recycling
        # keywords
        hs = JC.hbsolve(2*pi*8.1e9, w, src, (1,), Nh, circ, Dict();
            krylovrecycle = 8, krylovharvest = 3,
            krylovdeflationform = :floquet, keyedarrays = false)
        @test hs.nonlinear.solverinfo.converged
        @test any(k -> k.deflationsize > 0,
            hs.nonlinear.solverinfo.stages[1].krylov)
    end

    @testset "the subspace is inherited across a cached sweep" begin
        w = (2*pi*8e9,); Nh = (6,)
        src = [(mode = (1,), port = 1, current = 2.5e-6)]
        Ljs = range(100e-12, 102e-12; length = 3)
        noesc = (; krylovescalate = typemax(Int))
        builder = (; Lj) -> chain(12, Lj)

        cache = JC.hbcache(w, Nh, src, builder, (; Lj = Ljs[1]);
            method = :newtonkrylov, krylovrecycle = 8, krylovharvest = 3,
            krylovdeflationform = :floquet, krylovkwargs = noesc)
        first = Int[]
        for Lj in Ljs
            nl = JC.hbsolve!(cache, (; Lj = Lj))
            @test nl.solverinfo.converged
            st = nl.solverinfo.stages[end]
            push!(first, isempty(st.krylov) ? 0 : st.krylov[1].deflationsize)
        end
        # the first point starts cold; every later one starts from the
        # directions its predecessor found
        @test first[1] == 0
        @test all(>(0), first[2:end])
    end
end