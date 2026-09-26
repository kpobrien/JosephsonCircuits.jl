using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# The sparse factorizations behind a Newton step: the guarded factorize
# and solve, the non-conjugating transpose solve, and the ordering chosen
# by predicted fill.
@testset verbose=true "factorizations" begin

    @testset verbose=true "tryfactorize! error" begin

        begin
            factorization = JosephsonCircuits.KLUfactorization()
            J1 = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,2)
            cache = JosephsonCircuits.FactorizationCache()
            JosephsonCircuits.tryfactorize!(cache,factorization,J1)
            J2 = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[0.0, 0.0, 0.0, 0.0],2,2)
            # as of 2023-09-17 1.9.3 and older throws the first error and
            # 1.10.0-beta2 throws the second error
            @test_throws(
                str -> isequal("SingularException(0)",str) || 
                isequal("Unknown KLU error code: 2",str) ||
                isequal("SingularException: matrix is singular; factorization failed. Zero pivot found at index 0",str),
                JosephsonCircuits.tryfactorize!(cache,factorization,J2),
            )
        end

        begin
            cache = JosephsonCircuits.FactorizationCache()
            factorization = JosephsonCircuits.KLUfactorization()
            J3 = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,3)
            @test_throws DimensionMismatch JosephsonCircuits.tryfactorize!(
                cache, factorization, J3)
        end

        begin
            factorization = JosephsonCircuits.KLUfactorization()
            J1 = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,2)
            cache = JosephsonCircuits.FactorizationCache()
            JosephsonCircuits.tryfactorize!(cache,factorization,J1)
            J3 = JosephsonCircuits.sparse([1, 1, 2],[1, 2, 1],[1.3, 0.5, 0.1],2,2)

            @test_throws DimensionMismatch JosephsonCircuits.tryfactorize!(
                cache, factorization, J3)
        end

    end

    @testset verbose=true "tryfactorize! elseif path" begin
        A = [1.0 2.0; 3.0 5.0]
        cache = JosephsonCircuits.FactorizationCache()
        fact = JosephsonCircuits.QRfactorization()
        JosephsonCircuits.tryfactorize!(cache, fact, A)
        @test cache.factorization !== nothing
        A[1,1] = 1.1
        JosephsonCircuits.tryfactorize!(cache, fact, A)
        @test cache.factorization !== nothing
    end

    @testset "a refactorization whose pivots have grown is repivoted" begin
        # Node 4 couples to two nodes of a capacitive mesh and is
        # eliminated first, on its own diagonal, which vanishes at its
        # resonance f0 while the matrix stays well conditioned. A sweep
        # factorizes at 4 GHz and refactorizes at f0 with that pivot; its
        # scattering matrix at f0 must be the one a solve at f0 alone,
        # which pivots afresh, gives.
        w0 = 2pi*5.0e9
        L4, Ca, Cb = 1.0e-9, 30.0e-15, 50.0e-15
        C4 = (1 - 1e-13)/(L4*w0^2) - Ca - Cb
        c = Any[("P1", "3", "0", Port(1; Z0 = 50.0)),
            ("P2", "5", "0", Port(2; Z0 = 50.0))]
        nodes = ["1", "2", "3", "5"]
        k = 0
        for i in 1:4, j in i+1:4
            k += 1
            push!(c, ("Cx$(k)", nodes[i], nodes[j], Capacitor(40.0e-15*(1 + 0.3k))))
        end
        for (i, nd) in enumerate(nodes)
            push!(c, ("Lg$(i)", nd, "0", Inductor(2.0e-9*(1 + 0.2i))))
            push!(c, ("Cg$(i)", nd, "0", Capacitor(100.0e-15*(1 + 0.1i))))
        end
        append!(c, [("Ca", "4", "1", Capacitor(Ca)), ("Cb", "4", "2", Capacitor(Cb)),
            ("L4", "4", "0", Inductor(L4)), ("C4", "4", "0", Capacitor(C4))])
        circuit = Circuit(c)
        kw = (keyedarrays = false, returnQE = false, returnCM = false)
        atw0(S) = selectdim(S, ndims(S), size(S, ndims(S)))
        sweep = hblinsolve([2pi*4.0e9, w0], circuit, Dict{Symbol,Any}(); kw...)
        alone = hblinsolve([w0], circuit, Dict{Symbol,Any}(); kw...)
        @test isapprox(atw0(sweep.S), atw0(alone.S); atol = 1e-12)
        @test_throws ArgumentError JosephsonCircuits.KLUfactorization(pivottol = -1)
    end

    @testset "trysolvetranspose!" begin
        A = sparse([1,2,3,1,2,3], [1,1,2,3,3,3],
            Complex{Float64}[2.0, 1.0, 3.0, im, 4.0, 1.0], 3, 3)
        b = Complex{Float64}[1.0, 2.0, 3.0]
        x = zeros(Complex{Float64}, 3)
        for f in (JosephsonCircuits.KLUfactorization(),
                JosephsonCircuits.LUfactorization())
            cache = JosephsonCircuits.FactorizationCache()
            JosephsonCircuits.tryfactorize!(cache, f, A)
            fill!(x, 0)
            JosephsonCircuits.trysolvetranspose!(x, cache.factorization, b)
            @test isapprox(transpose(A)*x, b, rtol = 1e-10)
            # the non-conjugating transpose, not the adjoint
            @test !isapprox(adjoint(A)*x, b, rtol = 1e-10)
        end
    end
end

@testset verbose=true "kluordered: the ordering chosen by predicted fill" begin
    using SparseArrays, LinearAlgebra, Random
    JC = JosephsonCircuits
    rng = Random.default_rng()

    @testset "symbolicfill matches a Cholesky factorization" begin
        # a random sparse SPD matrix; the fill of L under a permutation is
        # what the elimination tree predicts, exactly
        n = 120
        B = sprand(rng, n, n, 0.03)
        S = B + B' + 4n*I
        for perm in (collect(1:n), randperm(rng, n), reverse(1:n))
            Sp = S[perm, perm]
            Lf = sparse(cholesky(Matrix(Sp)).L)
            fillcount, flops = JC.symbolicfill(S, perm)
            @test fillcount == nnz(Lf)
            @test flops == sum(abs2, [count(!iszero, Lf[j:end, j]) for j in 1:n])
        end
        @test_throws DimensionMismatch JC.symbolicfill(S, 1:n-1)
        @test_throws DimensionMismatch JC.symbolicfill(sprand(rng, 3, 4, 0.5), 1:3)
    end

    @testset "kluordered solves and never predicts worse than AMD" begin
        n = 400
        A = sprand(rng, n, n, 0.01) + 5I
        b = randn(rng, n)
        F = JC.kluordered(A)
        @test norm(A*(F\b) - b) <= 1e-10*norm(b)
        # the refactorization from new values reuses the ordering
        A2 = A + sparse(1.0I, n, n)
        JC.klurefactor!(F, A2, 1e-6)
        @test norm(A2*(F\b) - b) <= 1e-10*norm(b)
        # the chosen permutation is a permutation and its predicted flops
        # are at most AMD's
        S = JC._symmetricpattern(A)
        best = JC._bestordering(A)
        perm = best.perm
        @test isperm(perm)
        # with the fill it predicts, which is the elimination tree's
        @test best.fill == JC.symbolicfill(S, perm)[1]
        common = JC.CHOLMOD.getcommon()
        pamd = Vector{Int64}(undef, n)
        @test JC.LibSuiteSparse.cholmod_l_amd(JC.CHOLMOD.Sparse(S, 1), C_NULL, 0, pamd, common) == 1
        pamd .+= 1
        @test JC.symbolicfill(S, perm)[2] <= JC.symbolicfill(S, pamd)[2]
        # a grid-like pattern, where nested dissection is the better
        # ordering and must be the one taken: a 3d grid Laplacian
        m = 14
        idx(i, j, k) = i + m*(j - 1) + m^2*(k - 1)
        I3 = Int[]; J3 = Int[]
        for i in 1:m, j in 1:m, k in 1:m
            for (di, dj, dk) in ((1, 0, 0), (0, 1, 0), (0, 0, 1))
                i + di <= m && j + dj <= m && k + dk <= m || continue
                push!(I3, idx(i, j, k)); push!(J3, idx(i + di, j + dj, k + dk))
            end
        end
        G = sparse(I3, J3, -1.0, m^3, m^3)
        G = G + G' + 7I
        permg = JC._bestordering(G).perm
        Sg = JC._symmetricpattern(G)
        pamdg = Vector{Int64}(undef, m^3)
        JC.LibSuiteSparse.cholmod_l_amd(JC.CHOLMOD.Sparse(Sg, 1), C_NULL, 0, pamdg, common)
        pamdg .+= 1
        @test JC.symbolicfill(Sg, permg)[2] < JC.symbolicfill(Sg, pamdg)[2]
        # the symmetric pattern of both is that of X + X', whether X is
        # structurally symmetric, as the grid is, or not
        pattern(X) = (X = sparse(X); (SparseArrays.getcolptr(X), rowvals(X)))
        @test pattern(Sg) == pattern(G + G')
        @test pattern(S) == pattern(abs.(A) + abs.(A)')
        Fg = JC.kluordered(G)
        bg = randn(rng, m^3)
        @test norm(G*(Fg\bg) - bg) <= 1e-10*norm(bg)
        # the factorization the package hands out is this one
        @test JC.factorize(JC.KLUfactorization(), G) isa typeof(Fg)
    end

    @testset "a cache keeps the ordering of its pattern, and takes a seeded one" begin
        n = 300
        A = sprand(rng, n, n, 0.02) + 5I
        b = randn(rng, n)
        f = JC.KLUfactorization()
        cache = JC.FactorizationCache()
        JC.tryfactorize!(cache, f, A)
        ordering = cache.ordering
        @test ordering isa JC.FillOrdering
        @test isperm(ordering.perm)
        # a fresh factorization of the same pattern, with new values, takes
        # the ordering the cache holds and solves as a fresh choice does
        A2 = copy(A)
        nonzeros(A2) .*= 1 .+ rand(rng, nnz(A2))
        cache.factorization = nothing
        JC.tryfactorize!(cache, f, A2)
        @test cache.ordering === ordering
        @test cache.factorization\b ≈ JC.kluordered(A2)\b
        # a cache seeded with it factorizes as the one which chose it, and
        # so does one seeded with the bare permutation
        seeded = JC.seedordering!(JC.FactorizationCache(), A2, ordering)
        JC.tryfactorize!(seeded, f, A2)
        @test seeded.ordering === ordering
        @test seeded.factorization.q == cache.factorization.q
        bare = JC.seedordering!(JC.FactorizationCache(), A2, ordering.perm)
        JC.tryfactorize!(bare, f, A2)
        @test bare.factorization.q == cache.factorization.q
        @test_throws ArgumentError JC.seedordering!(JC.FactorizationCache(), A, 1:n-1)
        @test_throws ArgumentError JC.seedordering!(JC.FactorizationCache(), A,
            JC.FillOrdering(collect(1:n-1), n))
        # another pattern gets an ordering of its own
        B = A + sparse(1:n-1, 2:n, 1.0, n, n)
        cache.factorization = nothing
        JC.tryfactorize!(cache, f, B)
        @test cache.ordering == JC.fillordering(f, B)
        @test cache.ordering != ordering
        @test cache.factorization\b ≈ Matrix(B)\b
    end

    @testset "KLU sizes the factors of a given ordering by its fill" begin
        # handed only a permutation, KLU reserves ten times the matrix for
        # each factor; handed the fill the ordering choice predicted, it
        # reserves that. A banded matrix fills little, so the reservation
        # is most of what a fresh factorization allocates
        n = 20000
        A = spdiagm(-2 => fill(-1.0, n - 2), -1 => fill(0.5, n - 1),
            0 => fill(4.0, n), 1 => fill(0.5, n - 1), 2 => fill(-1.0, n - 2))
        o = JC.fillordering(JC.KLUfactorization(), A)
        JC.kluordered(A, o); JC.kluordered(A, o.perm)
        # (measured inside a function: `@allocated` has Julia compile the
        # whole top-level expression it is written in, here the file's
        # testset)
        allocations(A, ordering) =
            @allocated JosephsonCircuits.kluordered(A, ordering)
        sized = allocations(A, o)
        unsized = allocations(A, o.perm)
        @test sized < unsized/2
        b = randn(rng, n)
        @test JC.kluordered(A, o)\b ≈ JC.kluordered(A, o.perm)\b
        @test norm(A*(JC.kluordered(A, o)\b) - b) <= 1e-12*norm(b)
    end
end
