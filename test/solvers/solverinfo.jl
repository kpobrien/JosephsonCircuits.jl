using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# What a solve reports when it stops: the stall rule and the message
# naming the budget that ran out.
@testset "residualstalled and stallmessage" begin
    ps = JosephsonCircuits.residualstalled
    W = JosephsonCircuits.STALLHISTORY
    flat(n) = fill(1.0, n)
    geo(n, r) = [r^(k-1) for k in 1:n]
    # a history shorter than the rule needs is never a stall, whatever it
    # does
    @test !ps(flat(W - 1), 1)
    @test ps(flat(W), 1)
    # without a budget a steady geometric decrease is not a stall, however
    # slow
    @test !ps(geo(W, 0.5), 1)
    @test !ps(geo(W, 0.999), 1)
    @test !ps(geo(4*W, 0.99999), 1)
    # a residual which is flat or rising over the history is
    @test ps(flat(W), 1)
    @test ps(geo(W, 1.01), 1)          # steadily rising
    # a rise which is itself slowing counts as improving, not stalled,
    # because the later rate is better than the earlier one
    @test !ps([1.0 + 0.01*k for k in 1:W], 1)
    # a plateau shorter than the history the rule needs is not judged
    plateau = vcat([1.0, 0.05], fill(0.0475, 14))
    @test !ps(plateau, 1)
    # an accelerating history is never stopped
    @test !ps(vcat(flat(W), [0.9, 0.5, 0.1, 0.01]), 1)
    # the history starts where the caller says, and the whole of it is
    # judged: a plateau after a long descent is a stall only from a start
    # inside the plateau, or once it has outlasted the descent
    @test ps(vcat(geo(W, 0.5), flat(W)), W + 1)
    @test !ps(vcat(flat(W), geo(W, 0.5)), W + 1)
    descent = geo(5W, 0.9)
    plateau(n) = fill(descent[end], n)
    @test !ps(vcat(descent, plateau(W)), 1)
    @test ps(vcat(descent, plateau(W)), 5W - 1)
    @test ps(vcat(descent, plateau(20W)), 1)
    # the history the rule needs may be given explicitly
    @test !ps(flat(5), 1, 10)
    @test ps(flat(5), 1, 5)
    # with a tolerance and a budget of further steps, a rate which cannot
    # reach the tolerance within the budget is a stall too, and one which
    # can is not
    @test ps(geo(W, 0.999), 1; atol = 1e-8, remaining = 100)
    @test !ps(geo(W, 0.999), 1; atol = 1e-8, remaining = 100_000)
    @test !ps(geo(W, 0.5), 1; atol = 1e-8, remaining = 100)
    @test ps(flat(W), 1; atol = 1e-8, remaining = 1_000_000)
    # the history length and the acceleration guard come first, budget or
    # not
    @test !ps(geo(W - 1, 0.999), 1; atol = 1e-8, remaining = 1)
    @test !ps(vcat(flat(W), [0.9, 0.5, 0.1, 0.01]), 1; atol = 1e-8, remaining = 1)
    @test ps(geo(5, 0.999), 1, 5; atol = 1e-8, remaining = 1)
    for r in (:iterations, :work, :linesearch, :progress, :external)
        @test occursin(r"budget|stall|external",
            JosephsonCircuits.stallmessage(r))
    end
    @test occursin("unknown", JosephsonCircuits.stallmessage(:unknown))
end
