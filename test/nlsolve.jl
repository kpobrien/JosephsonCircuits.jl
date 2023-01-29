using JosephsonCircuits
using Test

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

end