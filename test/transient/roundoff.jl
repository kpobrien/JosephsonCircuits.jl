using JosephsonCircuits, LinearAlgebra, Test

# The same passive two-state low-pass block in different state coordinates.
# T commutes with A, and (C/T)*(sI-A)^(-1)*(T*B) is unchanged.
function roundoffproblem(q)
    a = 2pi*5e9
    T = Diagonal([q, 1/q])
    block = RationalScattering(-a*Matrix{Float64}(I, 2, 2),
        T*(a*Matrix{Float64}(I, 2, 2)), 0.8*[0.0 1.0; 1.0 0.0]/T,
        zeros(2, 2); zref = 50.0)
    c = Circuit([(:p1, 1, 0, Port(1)), (:p2, 2, 0, Port(2)), (:b, 1, 2, block)])
    return transientproblem(c; sources = [TransientSource(1, t -> 1e-6*sinpi(2*1e9*t)),
        TransientSource(2, 0.0)])
end

@testset "Gauss-Legendre residual roundoff" begin
    @testset "a trial's floor is adopted only with its point" begin
        # The stale Newton step goes to x=1, where the residual is larger
        # than both the base residual and the trial's finite floor. Neither
        # that rejected floor nor a nonfinite one may accept the base x=0.
        JC = JosephsonCircuits
        for trialbound in (2.0, Inf)
            x, base, jac = [0.0], Ref(0.0), Ref(1.0)
            floors = ([0.0], [0.0])
            function residual!(norms, r, y, floor)
                r[1] = y[1] + 10sin(y[1]) - 1
                norms[1] = abs(r[1])
                floor[1] = abs(y[1]) > 0.5 ? trialbound : 0.0
                return nothing
            end
            baseres = (norms, r, y) -> (base[] = y[1]; residual!(norms, r, y, floors[1]))
            trialres = (norms, r, y) -> residual!(norms, r, y, floors[2])
            refresh = mask -> (jac[] = 1 + 10cos(base[]); nothing)
            solve = (c, r) -> (c[1] = r[1]/jac[]; (false, 0))
            accept = mask -> (mask[1] && (base[] = x[1]); nothing)
            result = JC.newtonsolve!(x, [0.0], [0.0], [0.0], [0.0], baseres, trialres,
                refresh, solve, [1e-12], 15, [false], true, false, JC.NewtonWork(JC.CPU(), 1);
                accept! = accept, roundoff = floors)
            @test result[1]
            @test abs(x[1] + 10sin(x[1]) - 1) <= 1e-12
            @test floors[1][1] == 0.0
        end
    end

    opts = (; dt = 2e-12, method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
    span = (0.0, 0.2e-9)
    ref = transientsolve(roundoffproblem(1.0), span; opts..., record = :states)
    currents = [1e-8*sinpi(2*1e9*t + k) for k in 1:2, t in ref.times]
    weights = [cospi(2*1.7e9*t + k) for k in 1:2, t in ref.times]
    tangent = transienttangent(ref, currents)
    adjoint = transientadjoint(ref, weights)

    @testset "a realization's coordinates do not set the accuracy" begin
        for q in (1e4, 1e6, 1e8)
            p = roundoffproblem(q)
            sol = transientsolve(p, span; opts..., record = :states)
            @test sol.outgoing ≈ ref.outgoing rtol=1e-10
            @test transienttangent(sol, currents).outgoing ≈ tangent.outgoing rtol=1e-10
            @test transientadjoint(sol, weights).currents ≈ adjoint.currents rtol=1e-10
        end
    end

    @testset "checkpoint windows reconstruct their convergence scales" begin
        # Replay in both orders, including one-step windows and a short
        # final window. No unrecorded residual workspace may affect it.
        p = roundoffproblem(1e6)
        for every in (1, 10, 33)
            cp = transientsolve(p, span; opts..., record = :checkpoints, checkpointevery = every)
            @test cp.outgoing ≈ ref.outgoing rtol=1e-10
            @test transienttangent(cp, currents).outgoing ≈ tangent.outgoing rtol=1e-10
            @test transientadjoint(cp, weights).currents ≈ adjoint.currents rtol=1e-10
        end
    end

    @testset "a record of checkpoints replays its solve exactly" begin
        # a driven damped junction in its chaotic regime, which amplifies
        # any departure of a replayed window from the solve's steps: every
        # window, replayed backward by the adjoint and forward by the
        # tangent, ends on the next checkpoint to roundoff
        Lj, C = 1e-9, 1e-12
        wp = 1/sqrt(Lj*C)
        wd = 2wp/3
        chaos = Circuit([("p1", "1", "0", Port(1; termination = nothing)), ("R", "1", "0", Resistor(2/(wp*C))),
            ("C", "1", "0", Capacitor(C)), ("Lj", "1", "0", JosephsonJunction(Lj))])
        ramp(t) = t <= 0 ? 0.0 : t >= 1e-9 ? 1.0 : (1 - cospi(t/1e-9))/2
        drive(t) = 1.5*JosephsonCircuits.phi0/Lj*ramp(t)*cos(wd*t)
        dt = 2pi/wd/100
        sol = transientsolve(transientproblem(chaos; sources = [TransientSource(1, drive)]), (0.0, 6000dt); dt,
            record = :checkpoints, checkpointevery = 200)
        w = zeros(1, length(sol.times))
        w[1, end] = 1.0
        @test all(isfinite, transientadjoint(sol, w; quantity = :voltage).currents)
        @test all(isfinite, transienttangent(sol, w).voltage)
    end

    @testset "a tolerance below roundoff still reaches the discrete solution" begin
        # A requested tolerance far below double precision exercises the
        # floor itself. It must still require corrections.
        for q in (1.0, 1e6)
            sol = transientsolve(roundoffproblem(q), span; dt = opts.dt,
                method = GaussLegendre(), rtol = 0.0, atol = 1e-30)
            @test sol.outgoing ≈ ref.outgoing rtol=1e-10
            @test sol.stats.newtoncorrections >= sol.stats.steps
        end
    end

    @testset "conditions retain their own floor" begin
        c = Circuit([(:p, 1, 0, Port(1)), (:c, 1, 0, Capacitor(1e-12))])
        weak = transientproblem(c; sources = [TransientSource(1, t -> 1e-12*sinpi(2*1e9*t))])
        strong = transientproblem(weak; sources = [TransientSource(1, t -> 1e-3*sinpi(2*1e9*t))])
        kw = (; dt = 1e-12, method = GaussLegendre(), rtol = 1e-15, atol = 1e-20)
        one = transientsolve(weak, (0.0, 0.1e-9); kw...)
        both = transientsolve([strong, weak], (0.0, 0.1e-9); kw...)
        @test both[2].voltage ≈ one.voltage rtol=1e-10
    end

    @testset "a quiet node does not set another's floor" begin
        # a node of 100 nF, apart from the circuit, beside a junction on
        # 100 fF: the step of the circuit stops where it does alone, and the
        # tangent and the adjoint stay transposes
        drive(t) = 0.5e-6*sinpi(2*5e9*t)*(t <= 0 ? 0.0 : t >= 0.5e-9 ? 1.0 : sinpi(t/1e-9)^2)
        core = [(:p, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(100e-15)), (:jj, 1, 0, JosephsonJunction(1e-9)),
            (:l2, 1, 2, Inductor(1e-9)), (:c2, 2, 0, Capacitor(200e-15))]
        quiet = [(:cb, 3, 0, Capacitor(1e-7)), (:rb, 3, 0, Resistor(50.0))]
        kw = (; dt = 1e-12, rtol = 1e-9, atol = 1e-30, record = :phases)
        alone = transientsolve(transientproblem(Circuit(core); sources = [TransientSource(1, drive)]), (0.0, 1e-9); kw...)
        beside = transientsolve(transientproblem(Circuit(vcat(core, quiet)); sources = [TransientSource(1, drive)]),
            (0.0, 1e-9); kw...)
        @test beside.outgoing ≈ alone.outgoing rtol=1e-10
        cur = [1e-9*sinpi(2*4.7e9*t) for _ in 1:1, t in beside.times]
        w = [cospi(2*4.9e9*t) for _ in 1:1, t in beside.times]
        forward = sum(w .* transienttangent(beside, cur).outgoing)
        @test sum(transientadjoint(beside, w).currents .* cur) ≈ forward rtol=1e-9
    end

    @testset "a tangent direction is solved apart from the others" begin
        a = 50e9
        u = [1.0, -1.0]
        block = RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2),
            reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
        c = Circuit([(:p1, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(0.3e-12)), (:blk, 1, 2, block),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2))])
        drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(t/1e-9)^2*sinpi(2*3e9*t)
        p = transientproblem(c; sources = [TransientSource(1, drive), TransientSource(2, 0.0)])
        sol = transientsolve(p, (0.0, 0.5e-9); dt = 1e-11,
            method = GaussLegendre(), rtol = 1e-12, record = :phases)
        cur = [1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k) for k in 1:2, t in sol.times]
        one = transienttangent(sol, cur).outgoing
        for epsilon in (1e-6, 1e-10, 1e-14)
            weak = transienttangent(sol, epsilon .* cur).outgoing
            both = transienttangent(sol, cat(cur, epsilon .* cur; dims = 3)).outgoing
            @test both[:, :, 1] ≈ one rtol=1e-10
            @test weak ./ epsilon ≈ one rtol=1e-10
            @test both[:, :, 2] ≈ weak rtol=1e-10
        end
    end
end
