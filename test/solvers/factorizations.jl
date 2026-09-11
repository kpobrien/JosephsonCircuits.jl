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
            @test_throws(
                DimensionMismatch(""),
                JosephsonCircuits.tryfactorize!(cache,factorization,J3),
            )
        end

        begin
            factorization = JosephsonCircuits.KLUfactorization()
            J1 = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,2)
            cache = JosephsonCircuits.FactorizationCache()
            JosephsonCircuits.tryfactorize!(cache,factorization,J1)
            J3 = JosephsonCircuits.sparse([1, 1, 2],[1, 2, 1],[1.3, 0.5, 0.1],2,2)

            @test_throws(
                DimensionMismatch(""),
                JosephsonCircuits.tryfactorize!(cache,factorization,J3),
            )
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
    rng = MersenneTwister(3)

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
        JC.klunzval!(F, A2)
        @test norm(A2*(F\b) - b) <= 1e-10*norm(b)
        # the chosen permutation is a permutation and its predicted flops
        # are at most AMD's
        S = JC._symmetricpattern(A)
        perm = JC._bestordering(A)
        @test isperm(perm)
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
        permg = JC._bestordering(G)
        Sg = JC._symmetricpattern(G)
        pamdg = Vector{Int64}(undef, m^3)
        JC.LibSuiteSparse.cholmod_l_amd(JC.CHOLMOD.Sparse(Sg, 1), C_NULL, 0, pamdg, common)
        pamdg .+= 1
        @test JC.symbolicfill(Sg, permg)[2] <= JC.symbolicfill(Sg, pamdg)[2]
        Fg = JC.kluordered(G)
        bg = randn(rng, m^3)
        @test norm(G*(Fg\bg) - bg) <= 1e-10*norm(bg)
        # the factorization the package hands out is this one
        @test JC.factorize(JC.KLUfactorization(), G) isa typeof(Fg)
    end
end
