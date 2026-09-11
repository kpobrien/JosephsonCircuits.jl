using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# The Newton solver: what it refuses, how it reports a factorization that
# fails under it, and the Anderson correction it accelerates with.
@testset verbose=true "nlsolve" begin

    @testset verbose=true "nlsolve errors" begin

        function fj!(F, J, x)
            if !(F == nothing)
                F[1] = (x[1]+3)*(x[2]^3-7)+18
                F[2] = sin(x[2]*exp(x[1])-1)
            end
            if !(J == nothing)
                J[1, 1] = x[2]^3-7
                J[1, 2] = 3*x[2]^2*(x[1]+3)
                u = exp(x[1])*cos(x[2]*exp(x[1])-1)
                J[2, 1] = x[2]*u
                J[2, 2] = u
            end
            return nothing
        end

        # x = [ 0.1, 1.2]
        # F = [0.0, 0.0]
        # J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2])
        # @test_throws(
        #     DimensionMismatch("Number of columns in C must equal number of columns in B."),
        #     JosephsonCircuits.nlsolve!(fj!, F, J, x)
        # )
        begin
            x = [ 0.1, 1.2]
            F = [0.0, 0.0, 0.0]
            J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2])
            @test_throws(
                DimensionMismatch("First axis of the Jacobian `J` must have the same length as the residual `F`."),
                JosephsonCircuits.nlsolve!(fj!, F, J, x)
            )
        end

        begin
            x = [ 0.1, 1.2, 1.0]
            F = [0.0, 0.0]
            J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2])
            @test_throws(
                DimensionMismatch("Second axis of Jacobian `J` must have the same length as the input `x`."),
                JosephsonCircuits.nlsolve!(fj!, F, J, x)
            )
        end

        begin
            x = [ 0.1, 1.2]
            F = [0.0, 0.0]
            J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,3)
            @test_throws(
                DimensionMismatch("The Jacobian `J` matrix must be square."),
                JosephsonCircuits.nlsolve!(fj!, F, J, x)
            )
        end
    end

    @testset verbose=true "nlsolve klu error" begin

        function fj!(F, J, x)
            if !(F == nothing)
                F[1] = (x[1]+3)*(x[2]^3-7)+18
                F[2] = sin(x[2]*exp(x[1])-1)
            end
            if !(J == nothing)
                J[1, 1] = 0
                J[1, 2] = 0
                J[2, 1] = 0
                J[2, 2] = 0
            end
            return nothing
        end

        x = [ 0.1, 1.2]
        F = [0.0, 0.0]
        J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,2)
        # as of 2025-03-03 1.9.3 and older throws the first error and
        # 1.10.0-beta2 throws the second error
        # 1.13.0-DEV throws the third error
        @test_throws(
            str -> isequal("SingularException(0)",str) || 
            isequal("Unknown KLU error code: 2",str) ||
            isequal("SingularException: matrix is singular; factorization failed. Zero pivot found at index 0",str),
            JosephsonCircuits.nlsolve!(fj!, F, J, x)
        )

    end

    @testset "anderson correction solves the real least squares problem" begin
        # the coefficients are constrained real because the harmonic balance
        # error operator is antilinear, so the reference is least squares on
        # the stacked real and imaginary parts. The correction is computed
        # through a thin QR rather than the Gram matrix, which would square
        # the condition number exactly when the history is near collinear.
        for T in (Float64, ComplexF64)
            for (n, depth, m) in ((40,5,5), (40,5,3), (30,4,4))
                s = JosephsonCircuits.AndersonState(zeros(T,n), depth)
                s.histcount = m
                s.histpos = m == depth ? 1 : m+1
                s.historyready = true
                for k in 1:m
                    s.deltafhistory[:,k] .= randn(T,n)
                    s.deltaxhistory[:,k] .= randn(T,n)
                end
                dx = randn(T,n)
                @test JosephsonCircuits.andersoncorrection!(s, dx)
                A = s.deltafhistory[:,1:m]
                Ar = T <: Complex ? [real(A); imag(A)] : A
                br = T <: Complex ? [real(dx); imag(dx)] : dx
                g = Ar \ br
                ref = zeros(T,n)
                for j in 1:m
                    ref .+= g[j].*s.deltaxhistory[:,j] .+ g[j].*s.deltafhistory[:,j]
                end
                @test isapprox(s.correction, ref; rtol = 1e-8)
            end
        end
    end

    @testset "anderson truncates a rank deficient history" begin
        n = 30
        s = JosephsonCircuits.AndersonState(zeros(n), 4)
        s.histcount = 4; s.histpos = 1; s.historyready = true
        c1 = randn(n); c2 = randn(n)
        s.deltafhistory[:,1] .= c1
        s.deltafhistory[:,2] .= c2
        s.deltafhistory[:,3] .= 2 .*c1 .- 3 .*c2   # exactly dependent
        s.deltafhistory[:,4] .= randn(n)
        for k in 1:4
            s.deltaxhistory[:,k] .= randn(n)
        end
        @test JosephsonCircuits.andersoncorrection!(s, randn(n))
        @test all(isfinite, s.correction)
    end
end
