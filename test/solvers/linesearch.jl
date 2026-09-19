using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test
using Logging

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
        # both safeguards bound the fitted step: the minimizer of this fit
        # is 0.4525
        @test JosephsonCircuits.quadratic_trial_step(0.5, 0.605, -1.0;
            safeguard_high = 0.2)[1] == 0.2
        @test JosephsonCircuits.quadratic_trial_step(0.5, 0.605, -1.0;
            safeguard_low = 0.46, safeguard_high = 0.6)[1] == 0.46
        @test JosephsonCircuits.quadratic_trial_step(0.5, 0.605, -1.0)[1] ≈ 1/2.21
        @test_throws ArgumentError JosephsonCircuits.quadratic_trial_step(
            0.5, 0.605, -1.0; safeguard_low = 0.3, safeguard_high = 0.2)
        # a full step whose residual overflows is halved within the
        # safeguards too, with the linear estimate at that step
        for bad in (Inf, NaN)
            @test JosephsonCircuits.quadratic_trial_step(0.5, bad, -1.0; safeguard_high = 0.2) ==
                (0.2, 0.5 - 0.2, false)
            @test JosephsonCircuits.quadratic_trial_step(0.5, bad, -1.0;
                safeguard_low = 0.6, safeguard_high = 0.7)[1] == 0.6
        end
    end

    @testset verbose=true "the safeguards reach the first proposal" begin
        # a scalar residual whose quadratic fit through the full step
        # proposes 0.4525: with the upper safeguard below that the search
        # tries the full step and then the safeguard, which meets the
        # Armijo condition and is the step
        trials = Float64[]
        f! = (F, x) -> (push!(trials, x[1]); F[1] = 1 - x[1] + 1.1*x[1]^2; nothing)
        x0, deltax, F = [0.0], [1.0], [1.0]
        α, ϕα, accepted, backtracks = JosephsonCircuits.backtracking_linesearch!(
            f!, F, similar(x0), x0, deltax, 0.5, -1.0;
            ls = Backtracking(safeguardhigh = 0.2))
        @test trials == [1.0, 0.2]
        @test α == 0.2 && accepted && backtracks == 1
        @test ϕα ≈ (1 - 0.2 + 1.1*0.04)^2/2
        # the same with the full step's residual overflowing
        empty!(trials)
        g! = (F, x) -> (push!(trials, x[1]); F[1] = x[1] == 1 ? Inf : 1 - x[1] + 1.1*x[1]^2; nothing)
        α, ϕα, accepted, backtracks = JosephsonCircuits.backtracking_linesearch!(
            g!, [1.0], similar(x0), x0, deltax, 0.5, -1.0;
            ls = Backtracking(safeguardhigh = 0.2))
        @test trials == [1.0, 0.2]
        @test α == 0.2 && accepted
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

    @testset verbose=true "the Backtracking option" begin
        # every method interpolates by default, which suits the exact
        # step the direct loops take and the exact preconditioner the
        # Krylov loop picks for one tone
        @test Backtracking().interpolate
        @test Backtracking().safeguardlow == 0.1
        @test NewtonKrylov().linesearch.interpolate
        @test NewtonKrylov().linesearch.safeguardlow == 0.1
        @test Newton().linesearch.interpolate
        @test Newton().linesearch.safeguardlow == 0.1
        @test QuasiNewton().linesearch.interpolate
        b = Backtracking(interpolate = false, safeguardlow = 0.3,
            safeguardhigh = 0.7, c1 = 1e-3, maxbacktracks = 4, maxfailures = 1)
        @test !b.interpolate && b.safeguardlow == 0.3 && b.safeguardhigh == 0.7
        @test b.c1 == 1e-3 && b.maxbacktracks == 4 && b.maxfailures == 1
        @test NewtonKrylov(linesearch = b).linesearch === b
        @test Newton(linesearch = b).linesearch === b
        @test QuasiNewton(linesearch = b).linesearch === b
        # the option survives the escalation rewrite of a method
        @test JosephsonCircuits.withescalation(
            NewtonKrylov(linesearch = b), false).linesearch === b

        # each keyword validates itself, as the other option objects do
        @test_throws ArgumentError Backtracking(safeguardlow = 0.0)
        # the bounds match those the line searches enforce, so a value they
        # would reject is refused at construction
        @test_throws ArgumentError Backtracking(safeguardlow = 0.5)
        @test_throws ArgumentError Backtracking(safeguardlow = 0.3,
            safeguardhigh = 0.3)
        @test_throws ArgumentError Backtracking(c1 = 0.6)
        @test_throws ArgumentError Backtracking(safeguardlow = 0.6,
            safeguardhigh = 0.5)
        @test_throws ArgumentError Backtracking(safeguardhigh = 1.0)
        @test_throws ArgumentError Backtracking(c1 = 0.0)
        @test_throws ArgumentError Backtracking(c1 = 1.0)
        @test_throws ArgumentError Backtracking(maxbacktracks = -1)
        @test_throws ArgumentError Backtracking(maxfailures = 0)
    end

    @testset verbose=true "the option reaches both loops" begin
        # a pumped junction hard enough that the step is shortened, so the
        # two ways of shortening it leave different traces while reaching
        # the same point
        circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(circuit, ("P1","1","0",1)); push!(circuit, ("R1","1","0",:Rleft))
        push!(circuit, ("C1","1","2",:Cc)); push!(circuit, ("Lj1","2","0",:Lj))
        push!(circuit, ("C2","2","0",:Cj))
        defs = Dict{Symbol,Complex{Float64}}(:Lj => 1000e-12, :Cc => 100.0e-15,
            :Cj => 1000e-15, :Rleft => 50.0)
        wp = (2*pi*5e9,)
        src = [(mode = (1,), port = 1, current = 2.0e-6)]
        run(m; kw...) = hbnlsolve(wp, (8,), src, circuit, defs;
            keyedarrays = false, method = m, kw...)
        # halving can only ever accept a power of one half; interpolating
        # is not so constrained, which is the whole of the difference
        halfpower(a) = a == 0 || isapprox(log2(a), round(log2(a)); atol = 1e-12)
        steps(sol) = filter(!isnan, sol.solverinfo.stages[end].alpha)

        # every loop takes the option, and it has to reach the line search
        # rather than merely be accepted by the method: the first step is
        # shortened, to a fitted length or to a power of one half. With a
        # direct current block the direct loop takes its explicit direct
        # current branch.
        for kw in ((;), (; dc = true)),
                M in (NewtonKrylov, Newton, QuasiNewton)
            interp = run(M(linesearch = Backtracking(interpolate = true)); kw...)
            halve = run(M(linesearch = Backtracking(interpolate = false)); kw...)
            @test interp.solverinfo.converged && halve.solverinfo.converged
            @test isapprox(interp.nodeflux, halve.nodeflux; rtol = 1e-6)
            @test all(halfpower, steps(halve))
            @test any(a -> !halfpower(a), steps(interp))
            @test any(<(1), steps(interp))
        end

        # the fields have to reach the line search, not merely be accepted
        # by the method: with no trial budget the step cannot be shortened,
        # so a circuit which needs it shortened stalls on the line search
        # instead of converging
        for m in (Newton(linesearch = Backtracking(maxbacktracks = 0)),
                  NewtonKrylov(linesearch = Backtracking(maxbacktracks = 0)))
            starved = @test_logs (:warn, r"did not converge") match_mode=:any run(m)
            @test !starved.solverinfo.converged
            @test starved.solverinfo.stages[end].reason === :linesearch
        end
    end
end
