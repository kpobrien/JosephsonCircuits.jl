using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# What a solve reports when it stops: the stall projection and the
# message naming the budget that ran out.
@testset "projectedstall and stallmessage" begin
    ps = JosephsonCircuits.projectedstall
    # too short a window is never a stall
    @test !ps([1.0, 0.5, 0.25], 1, 1e-8, 100)
    # a steady rate of one half needs 23 more steps: a stall only when
    # fewer remain
    @test !ps([1.0, 0.5, 0.25, 0.125, 0.0625], 1, 1e-8, 100)
    @test ps([1.0, 0.5, 0.25, 0.125, 0.0625], 1, 1e-8, 10)
    # no decrease at all
    @test ps([1.0, 1.0, 1.0, 1.0, 1.0], 1, 1e-8, 1000)
    @test ps([1.0, 0.999, 1.0, 0.999, 1.0], 1, 1e-8, 1000)
    # an accelerating history is never stopped, however far from the
    # tolerance
    @test !ps([1.0, 0.99, 0.98, 0.9, 0.5], 1, 1e-8, 5)
    # the window starts where the caller says
    @test ps(vcat([1.0, 0.5, 0.25, 0.1, 0.05], ones(5)), 6, 1e-8, 1000)
    @test !ps(vcat(ones(5), [1.0, 0.5, 0.25, 0.125, 0.0625]), 6, 1e-8, 100)
    for r in (:iterations, :work, :linesearch, :progress, :external)
        @test occursin(r"budget|stall|external",
            JosephsonCircuits.stallmessage(r))
    end
    @test occursin("unknown", JosephsonCircuits.stallmessage(:unknown))
end
