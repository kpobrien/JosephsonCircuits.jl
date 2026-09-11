using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# The line search's trial step: the quadratic model of the merit function
# and the arguments it refuses.
@testset verbose=true "linesearch" begin

    @testset verbose=true "linesearch" begin
        @test(all(isapprox.(
            JosephsonCircuits.quadratic_trial_step(0.0,-0.22,-0.02),
            (1.0, -0.22,true),
        )))
        @test(all(isapprox.(
            JosephsonCircuits.quadratic_trial_step(0.0,0.0,-0.2),
            (0.5, -0.05000000000000001, false),
        )))
        @test(all(isapprox.(
            JosephsonCircuits.quadratic_trial_step(0.1,NaN,-0.02),
            (0.5, 0.09000000000000001, false),
        )))
    end

    @testset verbose=true "linesearch error" begin

        @test_throws(
            ArgumentError("`dϕ0dα` = 0.0 must be finite and negative."),
            JosephsonCircuits.quadratic_trial_step(0.0,0.2,0.0)
        )

        @test_throws(
            ArgumentError("`ϕ0` = NaN must be finite."),
            JosephsonCircuits.quadratic_trial_step(NaN,0.0,-0.02)
        )

    end
end
