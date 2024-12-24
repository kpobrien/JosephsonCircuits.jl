using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "nlsolve" begin

    @testset verbose=true "linesearch" begin
        @test(JosephsonCircuits.linesearch(0.0,-0.22,-0.02,0.0) == (1.0, -0.22))
        @test(JosephsonCircuits.linesearch(0.0,0.2,0.0,0.0) == (1.0, 0.2))
        @test(JosephsonCircuits.linesearch(0.0,0.18000000000000002,-0.02,0.1) == (1.0, 0.18000000000000002))
        @test(JosephsonCircuits.linesearch(0.1,0.1,0.0,0.0) == (1.0, 0.1))
        @test(JosephsonCircuits.linesearch(0.0,0.0,-0.2,0.0) == (0.5, -0.05000000000000001))
    end

    @testset verbose=true "linesearch error" begin

        @test_throws(
            ErrorException("NaN in nonlinear solver."),
            JosephsonCircuits.linesearch(0.0,NaN,-0.02,0.0)
        )

    end

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

        x = [ 0.1, 1.2]
        F = [0.0, 0.0, 0.0]
        J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2])
        @test_throws(
            DimensionMismatch("First axis of the Jacobian `J` must have the same length as the residual `F`."),
            JosephsonCircuits.nlsolve!(fj!, F, J, x)
        )

        x = [ 0.1, 1.2, 1.0]
        F = [0.0, 0.0]
        J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2])
        @test_throws(
            DimensionMismatch("Second axis of Jacobian `J` must have the same length as the input `x`."),
            JosephsonCircuits.nlsolve!(fj!, F, J, x)
        )

        x = [ 0.1, 1.2]
        F = [0.0, 0.0]
        J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,3)
        @test_throws(
            DimensionMismatch("The Jacobian `J` matrix must be square."),
            JosephsonCircuits.nlsolve!(fj!, F, J, x)
        )
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
        # as of 2023-09-17 1.9.3 and older throws the first error and
        # 1.10.0-beta2 throws the second error
        @test_throws(
            str -> isequal("SingularException(0)",str) || isequal("Unknown KLU error code: 2",str),
            JosephsonCircuits.nlsolve!(fj!, F, J, x)
        )

    end

    @testset verbose=true "tryfactorize! error" begin

        factorization = JosephsonCircuits.KLUfactorization()
        J1 = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2],2,2)
        cache = JosephsonCircuits.FactorizationCache()
        JosephsonCircuits.tryfactorize!(cache,factorization,J1)
        J2 = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[0.0, 0.0, 0.0, 0.0],2,2)
        # as of 2023-09-17 1.9.3 and older throws the first error and
        # 1.10.0-beta2 throws the second error
        @test_throws(
            str -> isequal("SingularException(0)",str) || isequal("Unknown KLU error code: 2",str),
            JosephsonCircuits.tryfactorize!(cache,factorization,J2)
        )
    end

end