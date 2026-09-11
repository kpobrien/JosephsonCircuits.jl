using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Random
using Test

# The scattering blocks and the transmission lines in time: constant
# blocks, rational blocks and blocks fitted to scattering data, and
# ideal lines with their delays, on the stepping, the tangent and the
# adjoint. The core of the solver is tested in transient.jl; this file is
# separate so that the two run on different workers.
@testset "scattering blocks and transmission lines in time" begin
    JC = JosephsonCircuits
    rc = [("P1", "1", "0", 1), ("R1", "1", "0", 50.0), ("C1", "1", "0", 1e-12)]

    @testset "constant scattering blocks" begin
        # a real constant block is its hybrid rows on auxiliary port
        # currents, nothing inverted: a through is a wire, an attenuator
        # equals its resistive network at equal and unequal reference
        # impedances, a short and an open are exact, and a circulator, a
        # lossless nonreciprocal block, makes the operator unsymmetric,
        # which the transposed solves of the adjoint must follow
        function pinetwork(S, R)
            Rh = Diagonal(sqrt.(R))
            Y = inv(Rh)*(I - S)*inv(I + S)*inv(Rh)
            return (1/(Y[1, 1] + Y[1, 2]), -1/Y[1, 2], 1/(Y[2, 2] + Y[2, 1]))
        end
        drive(t) = t <= 0 ? 0.0 : 1e-6*sinpi(t/1e-9)^2*sinpi(2*4e9*t)
        T, dt = 1e-9, 2e-12
        run(c; kw...) = transientsolve(transientproblem(c; sources = [TransientSource(1, drive)]), (0.0, T); dt,
            method = GaussLegendre(), rtol = 1e-12, kw...)
        through = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)),
            (:thru, 1, 2, ScatteringParameters([0.0 1.0; 1.0 0.0]; zref = 50.0)),
            (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        wire = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)), (:rw, 1, 2, Resistor(1e-3)),
            (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        pt = transientproblem(through; sources = [TransientSource(1, drive)])
        @test length(pt.blocks) == 1 && pt.Naux == 2 && pt.inertialess == [[3], [4]] && isempty(pt.algebraic)
        st, sw = run(through), run(wire)
        @test st.outgoing ≈ sw.outgoing rtol=1e-4
        @test st.voltage[1, :] ≈ st.voltage[2, :] rtol=1e-12
        @test st.stats.factorizations == 1
        @test_throws ArgumentError transientsolve(pt, (0.0, T); dt)
        for (g, R) in ((0.5, [50.0, 50.0]), (0.7, [50.0, 75.0]))
            Sa = [0.0 g; g 0.0]
            r1, rs, r2 = pinetwork(Sa, R)
            att = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)), (:jj, 1, 0, JosephsonJunction(2e-9)),
                (:att, 1, 2, ScatteringParameters(Sa; zref = R)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = R[2]))])
            res = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)), (:jj, 1, 0, JosephsonJunction(2e-9)),
                (:r1, 1, 0, Resistor(r1)), (:rs, 1, 2, Resistor(rs)), (:r2, 2, 0, Resistor(r2)),
                (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = R[2]))])
            @test run(att).outgoing ≈ run(res).outgoing rtol=1e-12
        end
        short = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)), (:l, 1, 2, Inductor(1e-9)),
            (:sh, 2, ScatteringParameters(fill(-1.0, 1, 1); zref = 50.0))])
        grounded = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)), (:l, 1, 0, Inductor(1e-9))])
        @test run(short).outgoing ≈ run(grounded).outgoing rtol=1e-12
        open = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)),
            (:op, 1, ScatteringParameters(fill(1.0, 1, 1); zref = 50.0))])
        alone = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12))])
        @test run(open).outgoing ≈ run(alone).outgoing rtol=1e-12
        # a block that depends on frequency has no realization
        @test_throws ArgumentError transientproblem(Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:f, 1, 2, ScatteringParameters(w -> [0.0 exp(-im*w*1e-12); exp(-im*w*1e-12) 0.0]; nports = 2, zref = 50.0)),
            (:p2, 2, 0, Port(2; Z0 = 50.0))]))
        # the circulator
        Sc = [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
        circ = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)),
            (:circ, 1, 2, 3, ScatteringParameters(Sc; zref = 50.0)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(1e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0)),
            (:c3, 3, 0, Capacitor(0.2e-12)), (:p3, 3, 0, Port(3; Z0 = 50.0))])
        probe(t, k) = 1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k)
        pc = transientproblem(circ; sources = [TransientSource(1, drive), TransientSource(2, t -> 0.0), TransientSource(3, t -> 0.0)])
        rec = transientsolve(pc, (0.0, T); dt, method = GaussLegendre(), rtol = 1e-12, record = :phases)
        energy = vec(sum(abs2, rec.outgoing; dims = 2))
        @test energy[2] > 5energy[3] > 5energy[1]
        currents = [probe(t, k) for k in 1:3, t in rec.times]
        tg = transienttangent(rec, currents)
        eps = 1e-4
        shifted(s) = transientproblem(circ; sources = [TransientSource(k, let k = k; t -> (k == 1 ? drive(t) : 0.0) + s*eps*probe(t, k); end) for k in 1:3])
        sp = transientsolve(shifted(1), (0.0, T); dt, method = GaussLegendre(), rtol = 1e-12)
        sm = transientsolve(shifted(-1), (0.0, T); dt, method = GaussLegendre(), rtol = 1e-12)
        @test tg.outgoing ≈ (sp.outgoing .- sm.outgoing) ./ (2eps) rtol=1e-6
        weights = [cospi(2*1.7e9*t + k) for k in 1:3, t in rec.times]
        ad = transientadjoint(rec, weights)
        @test sum(weights .* tg.outgoing) ≈ sum(ad.currents .* currents) rtol=1e-9
        cps = transientsolve(pc, (0.0, T); dt, method = GaussLegendre(), rtol = 1e-12, record = :checkpoints, checkpointevery = 25)
        @test transienttangent(cps, currents).outgoing ≈ tg.outgoing rtol=1e-9
        @test transientadjoint(cps, weights).currents ≈ ad.currents rtol=1e-9
        half = transientproblem(pc; sources = [TransientSource(1, t -> drive(t)/2), TransientSource(2, t -> 0.0), TransientSource(3, t -> 0.0)])
        b = transientsolve([pc, half], (0.0, T); dt, rtol = 1e-12, record = :phases)
        @test b.voltage[:, :, 1] == rec.voltage
        @test transientadjoint(b, weights).currents[:, :, 1] ≈ ad.currents rtol=1e-9
        # the JPA behind a 3 dB pad against the linearized solver
        g = 10^(-3/20)
        jpa = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:att, 1, 2, ScatteringParameters([0.0 g; g 0.0]; zref = 50.0)),
            (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)), (:cj, 3, 0, Capacitor(1000e-15))])
        fp, Ip, fs = 4.75001e9, 0.00565e-6, 4.76e9
        hb = hbsolve(2pi*[fs], (2pi*fp,), [(mode = (1,), port = 1, current = Ip)], (8,), (16,), jpa)
        s11 = hb.linearized.S((0,), 1, (0,), 1, 1)
        rise(t) = 1 - 2/(exp(t/10e-9) + exp(-t/10e-9))
        pj = transientproblem(jpa; sources = [TransientSource(1, t -> 2Ip*rise(t)*cospi(2fp*t))])
        steps = 80
        h = 1/(steps*fp)
        pump = transientsolve(pj, (0.0, 76000h); dt = h, method = GaussLegendre(), record = :checkpoints)
        cur = zeros(1, length(pump.times), 2)
        cur[1, :, 1] .= rise.(pump.times) .* cospi.(2fs .* pump.times)
        cur[1, :, 2] .= rise.(pump.times) .* sinpi.(2fs .* pump.times)
        sig = transienttangent(pump, cur)
        last = length(pump.times) - steps + 1:length(pump.times)
        v = sig.voltage[1, last, 1] .+ im .* sig.voltage[1, last, 2]
        S11 = 2*sum(v .* cispi.(-2fs .* pump.times[last]))/steps/50 - 1
        @test S11 ≈ s11 rtol=2e-4
    end


    @testset "rational scattering blocks" begin
        # a passive rational multiport, the form a vector fit delivers: the
        # series inductor between 50 ohm ports as a one state realization
        # is exact at every frequency in the linearized solver and, with
        # the exact elimination of its stage states, in time, where it
        # equals the explicit inductor to roundoff at every step size;
        # unstable and active realizations are refused
        R0, L = 50.0, 2e-9
        a = 2R0/L
        u = [1.0, -1.0]
        block = RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
        @test block.provider isa JC.RationalScatteringProvider && block.nports == 2
        w = 2pi*3e9
        S = zeros(ComplexF64, 2, 2, 1)
        JC.evaluateprovider!(S, block.provider, [w])
        @test S[:, :, 1] ≈ [im*w*L 2R0; 2R0 im*w*L] ./ (im*w*L + 2R0) rtol=1e-12
        @test_throws ArgumentError RationalScattering(fill(a, 1, 1), reshape(u, 1, 2), reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
        @test_throws ArgumentError RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-2a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
        # The passivity test finds the crossings of one by the pencil of
        # the equations and tests the largest singular value between them,
        # not at them: a block active below its crossing, S = 2/(s + 1)
        # with S(0) = 2, is refused, and so is a resonance peaking at
        # 1.05 over three percent of its frequency, which samples miss;
        # one peaking at 0.9 is accepted; a lossless block is certified
        # at more frequencies than its degree, and the pencil is singular
        # for it
        @test_throws ArgumentError RationalScattering(fill(-1.0, 1, 1), ones(1, 1), fill(2.0, 1, 1), zeros(1, 1))
        Ar, Br = [-0.05 1.0; -1.0 -0.05], [0.0; 1.0;;]
        @test_throws ArgumentError RationalScattering(Ar, Br, [1.05*0.1 0.0], zeros(1, 1))
        @test RationalScattering(Ar, Br, [0.9*0.1 0.0], zeros(1, 1)) isa ScatteringParameters
        crossings, _ = JC.passivitycrossings(Ar, Br, [1.05*0.1 0.0], zeros(1, 1))
        @test length(crossings) == 2 && 0.98 < crossings[1] < 1.0 < crossings[2] < 1.02
        @test isempty(JC.passivitycrossings(Ar, Br, [0.9*0.1 0.0], zeros(1, 1))[1])
        # nothing is inferred about a rational block's loss, and a declared
        # Lossless() is validated by the norms of the block and its inverse
        @test !JC.provablylossless(block.provider) && JC.losslessnorms(block.provider)
        @test isnothing(JC.passivitycrossings(block.provider.A, block.provider.B, block.provider.C, block.provider.D)[1])
        @test !JC.losslessnorms(RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-0.5a .* u, 2, 1), zeros(2, 2); zref = 50.0).provider)
        @test_throws ArgumentError RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-0.5a .* u, 2, 1), zeros(2, 2); zref = 50.0, noise = Lossless())
        # The passivity correction solves inequalities, and nonnegative
        # multipliers are only half of what that takes: the correction
        # has to meet the constraints too, since a constraint retired
        # from the working set can still be violated by the correction.
        # Four constraints in four unknowns whose answer keeps the
        # first, third and fourth active: feasibility and
        # complementarity are both checked.
        Aqp = [-1.1906621828062782 0.42012187633737086 -0.6172507381537534 -1.4612716694208;
               -0.23945086068732713 1.0930588528521001 -0.7800231394880406 1.5077262436554324;
               -0.13120004932279417 0.08933336796773754 -0.8860335576592437 0.3404464422612151;
               -1.2936863507016898 1.2907917529108746 -0.010470479326517982 -0.45829556655455883]
        bqp = [-0.7128670438978472, -0.7794876860372566, -1.411668334100207, -0.4487100476703757]
        μqp, okqp = JC.dualactiveset(Symmetric(Aqp*transpose(Aqp)), bqp)
        xqp = -(transpose(Aqp)*μqp)
        @test okqp
        @test all(>=(0), μqp)                              # dual feasibility
        @test maximum(Aqp*xqp .- bqp) <= 1e-9              # primal feasibility
        @test all(k -> μqp[k]*(bqp[k] - (Aqp*xqp)[k]) <= 1e-9, 1:4)   # complementarity
        @test findall(>(1e-10), μqp) == [1, 3, 4]
        @test xqp ≈ [0.3000041255892028, -0.1742802574578624,
                     1.381323342946531, -0.39019323846030324] rtol=1e-7
        # Constraints with parallel gradients and different bounds
        # cannot both hold with equality, and only one of them is active
        # at the answer, so a dependent working set cannot be solved as
        # equalities. `min x'x/2` subject to `2x <= -2` and `x <= -1.5`
        # is the smallest case: the answer is `x = -1.5`, with the first
        # constraint slack.
        depsolve = (Ad, bd) -> begin
            μd, okd = JC.dualactiveset(Symmetric(Ad*transpose(Ad)), bd)
            xd = -(transpose(Ad)*μd)
            (xd, okd, maximum(Ad*xd .- bd))
        end
        for (Ad, bd, want) in (
                (reshape([2.0, 1.0], 2, 1), [-2.0, -1.5], [-1.5]),      # first redundant
                (reshape([2.0, 1.0], 2, 1), [-3.0, -1.0], [-1.5]),      # second redundant
                (reshape([1.0, 1.0], 2, 1), [-1.0, -1.0], [-1.0]),      # exact duplicates
                (reshape([1.0, 2.0, 0.5], 3, 1), [-1.0, -3.0, -0.4], [-1.5]),
                ([1.0 0.0; 2.0 0.0; 0.0 1.0], [-1.0, -3.0, -0.5], [-1.5, -0.5]))
            xd, okd, vd = depsolve(Ad, bd)
            @test okd
            @test vd <= 1e-12                       # feasible
            @test xd ≈ want atol=1e-10              # and the least such
        end
        # and the thresholds are relative to the constraints, not to
        # one: a correction to a block with little loss has them at
        # 1e-10, where a threshold with a unit floor would admit nothing
        for ε in (1.0, 1e-6, 1e-10, 1e-12)
            xs2, ok2, v2 = depsolve(reshape([2.0, 1.0], 2, 1), [-2.0ε, -1.5ε])
            @test ok2
            @test v2 <= 1e-12*ε
            @test only(xs2) ≈ -1.5ε rtol=1e-9
        end
        # opposed parallel rows are infeasible together, and are reported
        # as unsolved rather than as a correction
        _, okx, vx = depsolve(reshape([1.0, -1.0], 2, 1), [-1.0, -1.0])
        @test !okx && vx > 0
        # and it is the least correction, not merely a feasible one: on
        # problems small enough to enumerate, it matches the best over
        # every active set which is feasible at all
        # deterministic rather than seeded, so the case set does not move
        # with the version of the random number stream
        let
            for trial in 1:40
                mq = 2 + trial % 3
                nq = mq + trial % 4
                Aq = [sinpi(0.37*(i + 3j + 5trial)) + 0.4*cospi(0.11*(2i + j + trial))
                      for i in 1:mq, j in 1:nq]
                bq = [-(0.2 + abs(sinpi(0.29*(k + trial)))) for k in 1:mq]
                μq, okq = JC.dualactiveset(Symmetric(Aq*transpose(Aq)), bq)
                xq = -(transpose(Aq)*μq)
                @test okq && maximum(Aq*xq .- bq) <= 1e-8
                best = Inf
                for mask in 1:(2^mq - 1)
                    act = [k for k in 1:mq if (mask >> (k - 1)) & 1 == 1]
                    Aa = Aq[act, :]
                    rank(Aa) < length(act) && continue
                    xc = transpose(Aa)*(Symmetric(Aa*transpose(Aa)) \ bq[act])
                    maximum(Aq*xc .- bq) <= 1e-9 && (best = min(best, norm(xc)))
                end
                @test norm(xq) <= best*(1 + 1e-8)
            end
        end
        # The passivity test is the largest singular value over every
        # frequency by the level set iteration, which finds a peak however
        # narrow: an all pass times a notch or a peak of relative width
        # 1e-4 to 1e-7, deviating from unitarity by the square of the width
        # between samples and by one at its center. The peaks are refused
        # at their true norm of two, the notches accepted, and the notches
        # are not certified lossless, so their loss keeps its noise in the
        # linearized solver: a notch absorbing at 1 rad/s emits the noise
        # of that loss there.
        for epsilon in (1e-4, 1e-6, 1e-7), sgn in (-1.0, 1.0)
            An = [0.0 1.0; -1.0 -2epsilon]
            Cn = [0.0 sgn*2epsilon]
            Am = [An zeros(2, 1); Cn fill(-10.0, 1, 1)]
            Bm = [0.0; 1.0; 1.0;;]
            Cm = [Cn fill(-20.0, 1, 1)]
            Dm = ones(1, 1)
            worst, where = JC.hinfnorm(Am, Bm, Cm, Dm)
            @test !JC.losslessnorms(JC.RationalScatteringProvider(Am, Bm, Cm, Dm))
            if sgn > 0
                @test worst ≈ 2 rtol=1e-6
                @test where ≈ 1 rtol=1e-6
                @test_throws ArgumentError RationalScattering(Am, Bm, Cm, Dm)
            else
                @test worst ≈ 1 rtol=1e-6
                notch = RationalScattering(Am, Bm, Cm, Dm)
                hbn = hblinsolve([1.0], Circuit([(:p, 1, 0, Port(1)), (:block, 1, notch)]); keyedarrays = false, returnCnoise = true)
                @test abs(hbn.S[1, 1, 1]) < 1e-6
                @test real(hbn.Cnoise[1, 1, 1]) ≈ 1 rtol=1e-8
            end
        end
        @test JC.hinfnorm(fill(-1.0, 1, 1), ones(1, 1), fill(2.0, 1, 1), zeros(1, 1))[1] ≈ 2 rtol=1e-8
        # What the norm search establishes, and what it does not. Its
        # first value is a lower bound and its third the level its
        # termination reached, and the question of whether a block is
        # passive to a tolerance can fall between them, so the assessment
        # has three answers and not two.
        #
        # The lower bound is worth having sharp. A resonance peaks near
        # its pole's frequency but not at it, and for a narrow one that
        # difference is larger than the tolerance being tested:
        # `k/(s^2 + 2 z s + 1)` with `z = 3e-4` and `k` set to put the
        # exact peak at `1 + 2e-8` evaluates to 1.0000000087 at the pole
        # frequency itself, and the search around each pole finds the
        # peak.
        let z = 3e-4, kk = (1 + 2e-8)*2z*sqrt(1 - z^2)
            Ah = [0.0 1.0; -1.0 -2z]
            Bh = reshape([0.0, kk], 2, 1)
            Ch = reshape([1.0, 0.0], 1, 2)
            Dh = zeros(1, 1)
            wpk = sqrt(1 - 2z^2)
            peak = abs(only(Dh + Ch*((im*wpk*I - Ah) \ Bh)))
            @test peak ≈ 1 + 2e-8 rtol=1e-9
            lower, _, level = JC.hinfnorm(Ah, Bh, Ch, Dh)
            @test lower ≈ peak rtol=1e-12       # the peak itself, not a probe near it
            # and still a bracket, to the roundoff of evaluating it
            @test lower <= peak*(1 + 1e-12) && peak <= level
            verdict, _, _, _ = JC.passivityassessment(Ah, Bh, Ch, Dh; atol = 1e-8)
            @test verdict === :active
            @test_throws ArgumentError RationalScattering(Ah, Bh, Ch, Dh; zref = 50.0, atol = 1e-8)
            # and the enforcement brings it under
            Ae, Be, Ce, De = @test_logs (:warn,) match_mode = :any JC.enforcepassivity(
                Ah, Bh, Ch, Dh, 2pi .* collect(range(0.01, 1.0; length = 200)))
            @test abs(only(De + Ce*((im*wpk*I - Ae) \ Be))) <= 1 + 1e-8
            # The level is a bound only where the pencil
            # resolves the crossings of it. Here it resolves none, the
            # peak standing over the level by less than roundoff, so the
            # level is the lower bound inflated by twice the tolerance and
            # nothing more.
            @test isempty(JC.passivitycrossings(Ah, Bh, Ch, Dh)[1])
            @test level ≈ lower*(1 + 2e-8) rtol=1e-12
        end
        # A peak need not be near any pole, so the probe grid has to
        # span the magnitudes of the poles the system has rather than
        # fixed decades: `k s/((s + a)(s + b))` with `a = 2.76e-5`,
        # `b = 1` and `k = (1 + 2e-8)(a + b)` peaks at
        # `sqrt(a b) = 5.3e-3`, decades from both poles, which are
        # real.
        let aa = 2.7592661119815163e-5, bb = 1.0
            kk = (1 + 2e-8)*(aa + bb)
            Ar2 = [0.0 1.0; -aa*bb -(aa + bb)]
            Br2 = reshape([0.0, 1.0], 2, 1)
            Cr2 = reshape([0.0, kk], 1, 2)
            Dr2 = zeros(1, 1)
            wpk2 = sqrt(aa*bb)
            peak2 = abs(only(Dr2 + Cr2*((im*wpk2*I - Ar2) \ Br2)))
            @test peak2 ≈ 1 + 2e-8 rtol=1e-9
            @test all(isreal, eigvals(Ar2))          # both poles real
            lower2, _, _ = JC.hinfnorm(Ar2, Br2, Cr2, Dr2)
            @test lower2 ≈ peak2 rtol=1e-9
            @test JC.passivityassessment(Ar2, Br2, Cr2, Dr2; atol = 1e-8)[1] === :active
            @test_throws ArgumentError RationalScattering(Ar2, Br2, Cr2, Dr2;
                zref = 50.0, atol = 1e-8)
            Ae2, Be2, Ce2, De2 = @test_logs (:warn,) match_mode = :any JC.enforcepassivity(
                Ar2, Br2, Cr2, Dr2, exp.(range(log(1e-4), log(1e2); length = 101)))
            @test abs(only(De2 + Ce2*((im*wpk2*I - Ae2) \ Be2))) <= 1 + 1e-8
        end
        # A block whose largest singular value is exactly one, which every
        # lossless block has, cannot be called passive to a tolerance
        # below the search's own: its level stands at `1 + 2 rtol`. The
        # assessment says so rather than picking a side.
        let Al = fill(-1.0, 1, 1), Bl = ones(1, 1), Cl = fill(-2.0, 1, 1), Dl = ones(1, 1)
            # the all pass (s - 1)/(s + 1), of modulus one everywhere
            @test abs(only(Dl + Cl*((im*0.7*I - Al) \ Bl))) ≈ 1 rtol=1e-14
            @test JC.passivityassessment(Al, Bl, Cl, Dl; atol = 1e-8)[1] === :indeterminate
            @test JC.passivityassessment(Al, Bl, Cl, Dl; atol = 1e-6)[1] === :passive
            @test JC.passivityassessment(Al, Bl, 2 .* Cl, Dl; atol = 1e-8)[1] === :active
        end
        @test RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-0.5a .* u, 2, 1), zeros(2, 2); zref = 50.0) isa ScatteringParameters
        drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(t/1e-9)^2*sinpi(2*3e9*t)
        withblock = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)), (:blk, 1, 2, block),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        explicit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)), (:l, 1, 2, Inductor(L)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        pb = transientproblem(withblock; sources = [TransientSource(1, drive)])
        pe = transientproblem(explicit; sources = [TransientSource(1, drive)])
        @test JC.blockstates(pb) == 1 && transientstate(pb).states == [0.0]
        for dt in (4e-12, 1e-12)
            sb = transientsolve(pb, (0.0, 1.5e-9); dt, method = GaussLegendre(), rtol = 1e-12, record = :states)
            se = transientsolve(pe, (0.0, 1.5e-9); dt, method = GaussLegendre(), rtol = 1e-12)
            @test sb.outgoing ≈ se.outgoing rtol=1e-9
            @test sb.stats.factorizations == 1 && size(sb.blockstates) == (1, length(sb.times))
        end
        cps = transientsolve(pb, (0.0, 1.5e-9); dt = 2e-12, method = GaussLegendre(), rtol = 1e-12, record = :checkpoints, checkpointevery = 50)
        @test size(cps.checkpoints.states) == (1, 15)
        @test hblinsolve(2pi*[3e9], withblock; keyedarrays = false).S ≈ hblinsolve(2pi*[3e9], explicit; keyedarrays = false).S rtol=1e-12
        # the tangent and the adjoint through the states equal the explicit
        # inductor's, the adjoint the tangent, the states' cotangent a
        # tangent of the states, and the checkpoints the record
        rprobe(t, k) = 1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k)
        pb2 = transientproblem(withblock; sources = [TransientSource(1, drive), TransientSource(2, t -> 0.0)])
        pe2 = transientproblem(explicit; sources = [TransientSource(1, drive), TransientSource(2, t -> 0.0)])
        rb = transientsolve(pb2, (0.0, 1.5e-9); dt = 2e-12, method = GaussLegendre(), rtol = 1e-12, record = :phases)
        re = transientsolve(pe2, (0.0, 1.5e-9); dt = 2e-12, method = GaussLegendre(), rtol = 1e-12, record = :phases)
        rcurrents = [rprobe(t, k) for k in 1:2, t in rb.times]
        rtb, rte = transienttangent(rb, rcurrents), transienttangent(re, rcurrents)
        @test rtb.outgoing ≈ rte.outgoing rtol=1e-10
        rweights = [cospi(2*1.7e9*t + k) for k in 1:2, t in rb.times]
        rab, rae = transientadjoint(rb, rweights), transientadjoint(re, rweights)
        @test rab.currents ≈ rae.currents rtol=1e-10
        @test sum(rweights .* rtb.outgoing) ≈ sum(rab.currents .* rcurrents) rtol=1e-10
        dz = fill(1e-6, 1, 1)
        rtz = transienttangent(rb, zeros(2, length(rb.times)); initialstate = (zeros(4), zeros(4), zeros(0, 1, 1), dz))
        @test sum(rweights .* rtz.outgoing) ≈ sum(rab.initialstates .* dz[:, 1]) rtol=1e-10
        rcps = transientsolve(pb2, (0.0, 1.5e-9); dt = 2e-12, method = GaussLegendre(), rtol = 1e-12, record = :checkpoints, checkpointevery = 50)
        @test transienttangent(rcps, rcurrents).outgoing ≈ rtb.outgoing rtol=1e-9
        @test transientadjoint(rcps, rweights).currents ≈ rab.currents rtol=1e-9
    end

    @testset "a rational block fitted to scattering data" begin
        # the scattering data of an RLC two-port, tabulated as the
        # linearized solver gives it, fitted by vector fitting at its
        # three poles and at more: the fit is exact, the poles the RLC's,
        # the extra poles dropped, and in time the fitted block equals the
        # explicit RLC around a junction to roundoff; the raw fit is
        # returned on request, and the sampling is checked
        rlc(R) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l1, 1, 2, Inductor(1.5e-9)), (:c, 2, 0, Capacitor(0.6e-12)),
            (:r, 2, 0, R), (:l2, 2, 3, Inductor(1.0e-9)), (:p2, 3, 0, Port(2; Z0 = 50.0))])
        fs = collect(range(0.2e9, 12e9; length = 240))
        hb = hblinsolve(2pi .* fs, rlc(Resistor(120.0)); keyedarrays = false)
        data = ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0)
        fitted = RationalScattering(data, 3)
        # the realization is minimal: the RLC's residues have rank one, so
        # its three poles take three states
        @test fitted.provider isa JC.RationalScatteringProvider && size(fitted.provider.A) == (3, 3)
        # the RLC's poles are the finite eigenvalues of its descriptor
        # pencil in the node voltages and the inductor currents, each of
        # them a pole of the fit once per port
        L1, C1, R1, L2, Z = 1.5e-9, 0.6e-12, 120.0, 1.0e-9, 50.0
        E = zeros(5, 5); E[2, 2] = C1; E[4, 4] = L1; E[5, 5] = L2
        M = zeros(5, 5); M[1, 1] = 1/Z; M[1, 4] = 1; M[2, 2] = 1/R1; M[2, 4] = -1; M[2, 5] = 1; M[3, 3] = 1/Z; M[3, 5] = -1
        M[4, 1] = -1; M[4, 2] = 1; M[5, 2] = -1; M[5, 3] = 1
        order(x) = (imag(x), real(x))
        expected = sort(filter(x -> abs(x) < 1e12, eigvals(-M, E)); by = order)
        poles = sort(eigvals(fitted.provider.A); by = order)
        @test poles ≈ sort(expected; by = order) rtol=1e-6
        fit = zeros(ComplexF64, 2, 2, length(fs))
        JC.evaluateprovider!(fit, fitted.provider, 2pi .* fs)
        @test maximum(abs.(fit .- hb.S)) < 1e-10
        @test (JC.checkpassive(fitted.provider); true)
        more = RationalScattering(data, 8; frequencies = fs)
        @test size(more.provider.A) == (3, 3)
        JC.evaluateprovider!(fit, more.provider, 2pi .* fs)
        @test maximum(abs.(fit .- hb.S)) < 1e-10
        raw = RationalScattering(data, 4; passivity = false)
        @test size(raw.provider.A) == (3, 3)
        @test (JC.checkpassive(raw.provider); true)
        @test_throws ArgumentError RationalScattering(data, 4; frequencies = reverse(fs))
        @test_throws ArgumentError RationalScattering(fitted, 4)
        @test_throws ArgumentError RationalScattering(data, 0)
        # The tolerances of the fit are the caller's, since a block with
        # little loss needs them tighter than the defaults: the enforced
        # margin, the rounds of enforcement, and how much error a dropped
        # pole may cost. They are checked, and the exact fit above is
        # reached whatever they are, since it needs no enforcement.
        @test_throws ArgumentError RationalScattering(data, 3; margin = -1.0)
        @test_throws ArgumentError RationalScattering(data, 3; margin = Inf)
        @test_throws ArgumentError RationalScattering(data, 3; pruneslack = -0.5)
        @test_throws ArgumentError RationalScattering(data, 3; rounds = 0)
        for kw in ((; margin = 1e-9), (; rounds = 40), (; pruneslack = 0.0), (; pruneslack = 1.0))
            tuned = RationalScattering(data, 3; kw...)
            JC.evaluateprovider!(fit, tuned.provider, 2pi .* fs)
            @test maximum(abs.(fit .- hb.S)) < 1e-10
        end
        # a pole is kept unless the fit is as close without it: the
        # permissive threshold of a factor of two would let the error grow
        # by that much at every pole it drops
        @test size(RationalScattering(data, 8; pruneslack = 0.0).provider.A) == (3, 3)
        # A block which states what it does at zero frequency has the fit
        # meet it exactly, which nothing in the data can do when the
        # samples begin above zero. At zero this RLC is a shunt resistor
        # between two ports: its inductors are shorts and its capacitor an
        # open, so `S11 = -y/(2 + y)` and `S21 = 2/(2 + y)` for the
        # resistor's conductance `y` in units of the port's.
        y = 50.0/120.0
        S0 = [(-y/(2 + y)) (2/(2 + y)); (2/(2 + y)) (-y/(2 + y))]
        stated = ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0,
            dcmodel = JC.ScatteringDC(S0))
        told = RationalScattering(stated, 6)
        P = told.provider
        reached = real.(P.D .+ P.C*((0.0*I - P.A) \ P.B))
        # exactly, and through the passivity enforcement, not only out of
        # the residue solve: the condition is carried into the correction
        @test reached ≈ S0 atol=1e-12
        # and the fit is no worse for it, the statement being the truth
        JC.evaluateprovider!(fit, P, 2pi .* fs)
        @test maximum(abs.(fit .- hb.S)) < 1e-8
        # The passivity correction holds a stated value at zero frequency,
        # rather than imposing it in the residue solve and then walking
        # away from it. `S(0) = D + C X0` with `X0 = (-A)^-1 B`, so a
        # perturbation holds it exactly when `dC X0 + dD = 0`, which is
        # linear and is eliminated before the inequalities are solved; the
        # answer is then the least change which holds it and not the least
        # change which nearly does. This resonance needs correcting, and
        # correcting it without the condition moves its value at zero by
        # nine parts in a hundred.
        Ae = [-0.25 1.0; -1.0 -0.25]
        Be = [1.0 0.0; 0.0 1.0]
        Ce = [0.35 0.0; 0.0 0.35]
        De = [0.0 0.0; 0.0 0.0]
        wse = collect(range(0.05, 4.0; length = 160))
        S0e = real.(De .+ Ce*((0.0*I - Ae) \ Be))
        @test first(JC.hinfnorm(Ae, Be, Ce, De)) > 1.3
        @test opnorm(S0e) < 1
        held = @test_logs (:warn,) match_mode = :any JC.enforcepassivity(
            Ae, Be, Ce, De, wse; dc = S0e)
        @test real.(held[4] .+ held[3]*((0.0*I - held[1]) \ held[2])) ≈ S0e atol=1e-12
        @test first(JC.hinfnorm(held...)) <= 1 + 1e-8
        free = @test_logs (:warn,) match_mode = :any JC.enforcepassivity(
            Ae, Be, Ce, De, wse)
        @test first(JC.hinfnorm(free...)) <= 1
        @test maximum(abs, real.(free[4] .+ free[3]*((0.0*I - free[1]) \ free[2])) .- S0e) > 0.05
        # Contracting toward a statement is a different step from
        # contracting toward nothing, and the step is found on the path
        # rather than derived from a bound: the norm along the path is a
        # norm of something affine in the step, so it is convex in it,
        # and the steps which meet a target are an interval containing
        # zero whose end the search finds. The bounds do not stand in for
        # it. `1/M` answers only contraction toward nothing:
        # `0.5 + 0.501 s/(s + 1)` anchored at `S(0) = 0.5` plainly has a
        # step, while repeated `1/M` steps stall short of one. The
        # triangle bound `t <= (1 - a)/(M - a)`, `a` the statement's
        # norm, is sufficient and not necessary: at `a = 1` its
        # numerator vanishes and it says nothing, though a positive step
        # can exist.
        @test JC.contractionstep(τ -> 1.001*τ, 1.0) ≈ 1/1.001 rtol=1e-9
        @test JC.contractionstep(τ -> abs(1 - 2.001τ), 1.0) ≈ 2/2.001 rtol=1e-9
        @test JC.contractionstep(τ -> 0.5, 1.0) == 1        # nothing to do
        @test JC.contractionstep(τ -> 2.0, 1.0) == 0        # no step helps
        # the resolution of the step is the caller's; a coarse one still
        # returns a step which measures under the target
        @test JC.contractionstep(τ -> 1.001*τ, 1.0; tol = 1e-3) ≈ 1/1.001 atol=1e-3
        let Aa = fill(-1.0, 1, 1), Ba = ones(1, 1), Ca = fill(-0.501, 1, 1),
            Da = fill(1.001, 1, 1), S0a = fill(0.5, 1, 1)
            @test first(JC.hinfnorm(Aa, Ba, Ca, Da)) > 1
            got = JC.enforcepassivity(Aa, Ba, Ca, Da, [0.001, 0.01, 0.1]; dc = S0a)
            @test first(JC.hinfnorm(got...)) <= 1 + 1e-8
            @test only(real.(got[4] + got[3]*((0.0*I - got[1]) \ got[2]))) ≈ 0.5 atol=1e-12
            # and an anchor of unit norm is contracted toward, not
            # refused: `-1.001 + 2.001/(s + 1)` is anchored at
            # `S(0) = 1`, and `t = 2/2.001` takes it to
            # `(1 - s)/(1 + s)`, exactly all pass, holding the anchor
            got1 = JC.enforcepassivity(fill(-1.0, 1, 1), ones(1, 1), fill(2.001, 1, 1),
                fill(-1.001, 1, 1), [0.1, 1.0, 10.0]; dc = ones(1, 1))
            @test first(JC.hinfnorm(got1...)) <= 1 + 1e-8
            @test only(real.(got1[4] + got1[3]*((0.0*I - got1[1]) \ got1[2]))) ≈ 1 atol=1e-10
        end
        # The correction is assembled one output port block at a time,
        # which rests on the normal matrix being the same block for
        # every port: restricted to port `i`'s own unknowns the
        # perturbation reads `[X(w)' I]`, and that does not depend on
        # `i`. The property itself is tested.
        let Ab = [-0.25 1.0; -1.0 -0.25], Bb = [1.0 0.0; 0.0 1.0], nzb = 2, nb = 2
            blocks = JC.portblocks(nb, nzb)
            @test length(blocks) == nb
            @test sort(vcat(blocks...)) == 1:(nb*nzb + nb*nb)   # a partition
            for w in (0.3, 1.0, 2.5)
                X = (im*w*I - Ab) \ Bb
                M = zeros(ComplexF64, nb*nb, nb*nzb + nb*nb)
                for i in 1:nb, j in 1:nb
                    row = (i - 1)*nb + j
                    for k in 1:nzb
                        M[row, (k - 1)*nb + i] = X[k, j]
                    end
                    M[row, nb*nzb + (i - 1)*nb + j] = 1.0
                end
                first = real.(transpose(M[:, blocks[1]])*M[:, blocks[1]])
                for i in 2:nb
                    @test real.(transpose(M[:, blocks[i]])*M[:, blocks[i]]) ≈ first atol=1e-12
                end
            end
        end
        # and the correction it produces is feasible and makes the block
        # passive, on a resonance which needs one
        let Ad = [-0.25 1.0; -1.0 -0.25], Bd = [1.0 0.0; 0.0 1.0],
            Cd2 = [0.35 0.0; 0.0 0.35], Dd = [0.0 0.0; 0.0 0.0],
            wsd = collect(range(0.05, 4.0; length = 60))
            @test first(JC.hinfnorm(Ad, Bd, Cd2, Dd)) > 1
            got = @test_logs (:warn,) match_mode = :any JC.enforcepassivity(
                Ad, Bd, Cd2, Dd, wsd)
            @test first(JC.hinfnorm(got...)) <= 1
            S0d = real.(Dd .+ Cd2*((0.0*I - Ad) \ Bd))
            held = @test_logs (:warn,) match_mode = :any JC.enforcepassivity(
                Ad, Bd, Cd2, Dd, wsd; dc = S0d)
            @test first(JC.hinfnorm(held...)) <= 1 + 1e-8
            @test real.(held[4] .+ held[3]*((0.0*I - held[1]) \ held[2])) ≈ S0d atol=1e-12
        end
        # A statement of unit norm, which a through, an open and a short
        # all are, pins the norm of any fit meeting it at one: there is
        # nothing to contract away, and contracting toward anything else
        # would move the statement. Such a fit is accepted at one to the
        # tolerance it is validated against. Where the two cannot both
        # hold, the fit is refused and says so, rather than returning a
        # block which meets neither.
        @test_throws ArgumentError RationalScattering(ScatteringParameters(
            (2pi .* fs, hb.S); nports = 2, zref = 50.0, dcmodel = JC.ThroughDC()), 6)
        # a statement at odds with the poles the data gives is visible at
        # once, in the error it costs in band. That is not the data
        # judging the statement: nothing below the samples is measured,
        # and poles placed below the band would let the same fit meet a
        # value it refuses here. It is only that a fit cannot reach where
        # its poles do not go.
        wrong = ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0,
            dcmodel = JC.ThroughDC())
        pw, rw, Dw = JC.vectorfit(hb.S, 2pi .* fs, 6, 30; pruneslack = 0.05,
            dc = JC.dcscatteringmatrix(JC.ThroughDC(), 2))
        Aw, Bw, Cw = JC.realization(pw, rw, 2)
        JC.evaluateprovider!(fit, JC.RationalScatteringProvider(Aw, Bw, Cw, Dw), 2pi .* fs)
        @test maximum(abs.(fit .- hb.S)) > 1e-3
        # A sample at zero frequency is data like any other, and the fit
        # is scaled by the geometric centre of the band, which zero has
        # no part in: the lowest positive frequency stands in for it, in
        # the scaling and in the spread the starting poles are given. A
        # pole started at zero would sit on top of the sample there,
        # where the basis it belongs to is infinite, and the failure
        # would surface as infinities inside an eigensolver.
        withdc = vcat(0.0, fs)
        Sdc = zeros(ComplexF64, 2, 2, length(withdc))
        Sdc[:, :, 2:end] .= hb.S
        Sdc[:, :, 1] .= hb.S[:, :, 1]
        pdc, rdc, Ddc = JC.vectorfit(Sdc, 2pi .* withdc, 6, 30; pruneslack = 0.05)
        @test all(isfinite, pdc) && all(isfinite, rdc) && all(isfinite, Ddc)
        @test all(p -> real(p) < 0, pdc)
        # and a fit with no positive frequency has no band to be scaled by
        @test_throws ArgumentError JC.vectorfit(Sdc[:, :, 1:2], [0.0, 0.0], 4, 30)
        # Zero frequency is data, and the public fits take it: a sample
        # there is the block's value at DC, and the coordinates the fit
        # works in are normalized by the lowest positive frequency
        # rather than the lowest, so a band beginning at zero still has
        # a scale.
        Az = reshape([-1.0], 1, 1); Bz = reshape([1.0], 1, 1)
        Cz = reshape([0.5], 1, 1); Dz = reshape([0.0], 1, 1)
        wz = 2pi .* collect(range(0.0, 1.0; length = 20))
        Sz = zeros(ComplexF64, 1, 1, length(wz))
        for (k, w) in enumerate(wz)
            Sz[1, 1, k] = only(Dz + Cz*((im*w*I - Az) \ Bz))
        end
        blkz = ScatteringParameters((wz, Sz); nports = 1, zref = 50.0)
        for got in (RationalScattering(blkz, 1), RationalScattering(blkz; tol = 1e-6))
            Sf1 = zeros(ComplexF64, 1, 1, length(wz))
            JC.evaluatescattering!(Sf1, got, wz)
            @test maximum(abs, Sf1 .- Sz) < 1e-12
            @test real(Sf1[1, 1, 1]) ≈ 0.5 atol=1e-12
        end
        # a scan from minpoles needs minpoles + 1 samples, and too few
        # are refused by name rather than reported as a scan of no orders
        @test_throws ArgumentError RationalScattering(blkz; tol = 1e-6, minpoles = 30)
        # A tabulated block is evaluated at the angular frequencies it
        # stores, not at those divided by `2 pi` and multiplied by it
        # again: that round trip puts this table's first sample at
        # 0.19999999999999998, one unit in the last place outside the
        # table's own range.
        let wt = [0.2, 1.0]
            @test 2pi*(wt[1]/(2pi)) != wt[1]
            St = reshape(ComplexF64[0.3/(1 + im*w) for w in wt], 1, 1, :)
            bt = ScatteringParameters((wt, St); nports = 1, zref = 50.0)
            ft = RationalScattering(bt, 1)
            @test size(ft.provider.A, 1) == 1
        end
        # The samples are checked once, by one function: duplicates are
        # refused, since two samples at one frequency make the divided
        # differences of the Loewner pencil singular and leave a seeded
        # pole pair no width.
        @test JC.checkfrequencies([0.0, 1.0, 2.0]) == [0.0, 1.0, 2.0]
        @test_throws ArgumentError JC.checkfrequencies([0.0, 0.0])
        @test_throws ArgumentError JC.checkfrequencies([1.0, 1.0, 2.0])
        @test_throws ArgumentError JC.checkfrequencies([0.0, Inf])
        for bad in ([0.0, 0.1, 0.1, 0.2], [0.0, 0.0, 0.0, 0.0],
                    [-0.1, 0.0, 0.1, 0.2], [0.2, 0.1, 0.3, 0.4])
            @test_throws ArgumentError RationalScattering(blkz, 1; frequencies = bad)
        end
        # `vectorfit` is the fit and not the enforcement. Samples of
        # `1.5 - 1.2/(s + 1)` over a low band have an exact one pole fit
        # whose constant term is 1.5, and the fit returns it: bringing a
        # constant term under one is the enforcement's work, since no
        # perturbation over a band reaches infinite frequency, and it is
        # done for the caller who asks for it rather than to every fit.
        wsc = 2pi .* collect(range(0.01, 0.3; length = 40))
        Sc = zeros(ComplexF64, 1, 1, length(wsc))
        for (k, w) in enumerate(wsc)
            Sc[1, 1, k] = 1.5 - 1.2/(im*w + 1)
        end
        praw, rraw, Draw = JC.vectorfit(Sc, wsc, 1, 30)
        @test only(Draw) ≈ 1.5 rtol=1e-8
        @test only(praw) ≈ -1 rtol=1e-6
        @test only(rraw) ≈ -1.2 rtol=1e-6
        # and when it is asked for, the trigger and the target are the
        # caller's, not a threshold of the fit's own
        _, _, Dc = JC.vectorfit(Sc, wsc, 1, 30; constanttol = 1e-8, constantmargin = 1e-6)
        @test only(Dc) ≈ 1 - 1e-6 rtol=1e-9
        _, _, Dw = JC.vectorfit(Sc, wsc, 1, 30; constanttol = 1e-8, constantmargin = 1e-3)
        @test only(Dw) ≈ 1 - 1e-3 rtol=1e-9
        # a constant term under one is left alone whatever is asked
        @test JC.passiveconstant([0.5 0.0; 0.0 0.25]; tol = 1e-8)[2] == false
        # a value stated at zero frequency must be real, a real rational
        # function being real there
        @test_throws ArgumentError JC.vectorfit(Sc, wsc, 1, 30; dc = fill(0.5 + 0.1im, 1, 1))
        # More iterations of the pole relocation never return a worse
        # fit, because what it returns is the iterate measured to fit
        # best and not the last one, nor the one whose poles moved least
        # between steps: a relocation which is circling can pass its
        # closest fit on a step where the poles happen to be moving
        # quickly. A delay, which no order fits exactly, shows it.
        gsr = 2pi .* collect(range(0.5e9, 12e9; length = 120))
        Sdr = zeros(ComplexF64, 2, 2, length(gsr))
        for (k, w) in enumerate(gsr)
            ee = cis(-w*40e-12)
            Sdr[:, :, k] .= [0 ee; ee 0]
        end
        xsr = gsr ./ sqrt(first(gsr)*last(gsr))
        startr = ComplexF64[]
        for w in range(first(xsr), last(xsr); length = 3)
            push!(startr, complex(-0.01w, w))
            push!(startr, complex(-0.01w, -w))
        end
        errsr = [JC.fiterror(Sdr, xsr, JC.converge(Sdr, xsr, copy(startr), it)) for it in 1:12]
        @test all(k -> errsr[k + 1] <= errsr[k]*(1 + 1e-12), 1:length(errsr) - 1)
        @test minimum(errsr) == errsr[end]
        # The order can be searched for instead of given: the fewest poles
        # which hold the error over the samples under a tolerance, as a
        # fraction of the largest response. This RLC is exact at three
        # poles, so every tolerance it can meet it meets there, and the
        # search returns three however loose or tight the tolerance is.
        auto = RationalScattering(data; tol = 1e-6)
        @test size(auto.provider.A) == (3, 3)
        JC.evaluateprovider!(fit, auto.provider, 2pi .* fs)
        @test maximum(abs.(fit .- hb.S)) < 1e-10
        @test JC.relativefiterror(auto, hb.S, fs) <= 1e-6
        # the tolerance bounds what the caller receives: a looser one may
        # be met by fewer poles but never by a fit which misses it
        for tol in (1e-2, 1e-8)
            got = RationalScattering(data; tol = tol)
            @test JC.relativefiterror(got, hb.S, fs) <= tol
            @test size(got.provider.A, 1) <= 3
        end
        # minpoles is a floor the search does not fit below, and maxpoles
        # a ceiling: a tolerance no order up to it can meet is an error
        # naming the closest fit found rather than a block quietly worse
        # than was asked for
        @test size(RationalScattering(data; tol = 1e-6, minpoles = 6).provider.A, 1) >= 3
        # The search fits every order from `minpoles` up and returns the
        # first which meets the tolerance, so whatever error some order
        # achieves, the search asked for that error meets it and does so
        # with no more states. This is the guarantee, and it needs the
        # scan: more poles do not always fit better, so the orders
        # meeting a tolerance are a window rather than a tail, and a
        # search which skips orders can step over the window entirely.
        @test_logs (:warn,) match_mode = :any for np in 1:6
            reachable = try
                JC.relativefiterror(RationalScattering(data, np), hb.S, fs)
            catch e
                e isa ArgumentError ? Inf : rethrow()
            end
            isfinite(reachable) || continue
            found = RationalScattering(data; tol = reachable, minpoles = 1, maxpoles = 6)
            @test JC.relativefiterror(found, hb.S, fs) <= reachable
            @test size(found.provider.A, 1) <= 2np
        end
        # One evaluator for the relocation's choice of iterate, the
        # pruning, and the acceptance: it measures in the spectral norm
        # the public acceptance uses, and under the same zero frequency
        # condition the final residue solve will impose.
        let ps = [complex(-0.3, 1.0), complex(-0.3, -1.0)]
            xq = collect(range(0.4, 2.5; length = 30))
            Sq = zeros(ComplexF64, 2, 2, length(xq))
            for (k, x) in enumerate(xq)
                Sq[:, :, k] .= [0.2 0.9; 0.9 0.2] ./ (1 + im*x)
            end
            res, Dq = JC.fitresidues(Sq, xq, ps)
            byhand = maximum(eachindex(xq)) do k
                opnorm(Dq .+ sum(res[:, :, q] ./ (im*xq[k] - ps[q]) for q in eachindex(ps)) .-
                       view(Sq, :, :, k))
            end
            @test JC.fiterror(Sq, xq, ps) ≈ byhand
            # and the condition is carried, so the error reported is the
            # error of the fit that will be built
            resd, Dd = JC.fitresidues(Sq, xq, ps; dc = [0.1 0.8; 0.8 0.1])
            byhandd = maximum(eachindex(xq)) do k
                opnorm(Dd .+ sum(resd[:, :, q] ./ (im*xq[k] - ps[q]) for q in eachindex(ps)) .-
                       view(Sq, :, :, k))
            end
            @test JC.fiterror(Sq, xq, ps; dc = [0.1 0.8; 0.8 0.1]) ≈ byhandd
            @test JC.fiterror(Sq, xq, ps; dc = [0.1 0.8; 0.8 0.1]) > JC.fiterror(Sq, xq, ps)
        end
        # a stated value at zero survives the pruning: the poles that are
        # left still reach it, because every candidate was judged with it
        let xr = 2pi .* fs ./ sqrt(2pi*fs[1]*2pi*fs[end]),
            thru = JC.dcscatteringmatrix(JC.ThroughDC(), 2)
            start = ComplexF64[]
            for w in range(xr[1], xr[end]; length = 5)
                push!(start, complex(-0.01w, w)); push!(start, complex(-0.01w, -w))
            end
            settled = JC.converge(hb.S, xr, start, 30; dc = thru)
            kept = JC.prunepoles(hb.S, xr, copy(settled), 30, 0.05; dc = thru)
            @test !isempty(kept)
            resk, Dk = JC.fitresidues(hb.S, xr, kept; dc = thru)
            reached = real.(Dk .+ sum(resk[:, :, q] ./ (0.0 - kept[q]) for q in eachindex(kept)))
            @test reached ≈ thru atol=1e-10
        end
        # and the window is real, not hypothetical: a pure delay is
        # passive and irrational, so no order fits it exactly and the
        # error is not monotone in the order, with orders which fit
        # better than either of their neighbours. A tolerance only such
        # an order meets is found only by a scan. The orders and the
        # tolerances here are measured rather than written down, so this
        # does not depend on where a particular machine's arithmetic
        # puts the noise floor.
        delay = 40e-12
        gs = collect(range(0.5e9, 12e9; length = 200))
        Sdelay = zeros(ComplexF64, 2, 2, length(gs))
        for (k, g) in enumerate(gs)
            e = cis(-2pi*g*delay)
            Sdelay[:, :, k] .= [0 e; e 0]
        end
        delayed = ScatteringParameters((2pi .* gs, Sdelay); nports = 2, zref = 50.0)
        # `pruneslack` is a budget for the pruning and not for each
        # deletion: measured per deletion, a run of them could each
        # spend the whole slack and the fit drift as far from where it
        # started as the number of deletions allowed.
        let xd = 2pi .* gs ./ sqrt(2pi*gs[1]*2pi*gs[end])
            start = ComplexF64[]
            for w in range(xd[1], xd[end]; length = 8)
                push!(start, complex(-0.01w, w)); push!(start, complex(-0.01w, -w))
            end
            settled = JC.converge(Sdelay, xd, start, 30)
            base = JC.fiterror(Sdelay, xd, settled)
            for slack in (0.0, 0.05, 0.5)
                kept = JC.prunepoles(Sdelay, xd, copy(settled), 30, slack)
                @test length(kept) <= length(settled)
                @test JC.fiterror(Sdelay, xd, kept) <= (1 + slack)*base + JC.roundoff(Sdelay)
            end
        end
        reach = np -> try
            JC.relativefiterror(RationalScattering(delayed, np), Sdelay, gs)
        catch e
            e isa ArgumentError ? Inf : rethrow()
        end
        errs = @test_logs (:warn,) match_mode = :any [reach(np) for np in 1:16]
        window = [np for np in 2:15 if errs[np] < min(errs[np-1], errs[np+1])]
        @test !isempty(window)
        @test_logs (:warn,) match_mode = :any for np in window
            # a tolerance between what this order reaches and what the
            # better of its neighbours reaches: only this order meets it
            tol = sqrt(errs[np]*min(errs[np-1], errs[np+1]))
            got = RationalScattering(delayed; tol = tol, minpoles = np - 1, maxpoles = np + 1)
            @test JC.relativefiterror(got, Sdelay, gs) <= tol
        end
        # The degree the samples determine budgets the search rather than
        # walling it. It is the numerical rank of a pencil built along
        # cycling directions from at most four hundred samples, with the
        # constant term contributing to it, so it can be short of what a
        # block needs; where the scan reaches it without meeting the
        # tolerance and the error is still falling, the search goes on.
        # This delay needs six poles and is estimated at four once the
        # noise floor is put high enough, and the search finds the six.
        @test JC.supporteddegree(Sdelay, 2pi .* gs, 1e-2) == 4
        expanded = @test_logs (:warn,) match_mode = :any RationalScattering(
            delayed; tol = 1e-8, noisefloor = 1e-2)
        @test size(expanded.provider.A, 1) ÷ 2 > 4
        @test JC.relativefiterror(expanded, Sdelay, gs) <= 1e-8
        # a `maxpoles` given by the caller is a wall, because the caller
        # made it one
        @test_logs (:warn,) match_mode = :any @test_throws ArgumentError RationalScattering(
            delayed; tol = 1e-8, noisefloor = 1e-2, maxpoles = 4)
        # and the expansion stops rather than running to the sample count:
        # a tolerance nothing reaches is still reported
        @test_logs (:warn,) match_mode = :any @test_throws ArgumentError RationalScattering(
            delayed; tol = 1e-16, maxpoles = 8)
        # The degree the samples determine is a property of the data and
        # not of the unit its frequencies are written in. The two halves
        # of the Loewner pencil do not carry the same units, a divided
        # difference of the response being an inverse frequency and a
        # shifted one dimensionless, so a threshold on the pencil built
        # at the frequencies as they come counts a different number as
        # the unit changes; normalized to the band, the count is the
        # same in every unit.
        @test allequal(JC.supporteddegree(Sdelay, (2pi .* gs) .* c, 1e-10)
                       for c in (1e12, 1e6, 1.0, 1e-6, 1e-12))
        @test allequal(JC.supporteddegree(Sdelay, (2pi .* gs) .* c, 1e-12)
                       for c in (1e9, 1.0, 1e-9))
        @test_throws ArgumentError RationalScattering(data; tol = 1e-16, maxpoles = 8)
        @test_throws ArgumentError RationalScattering(data; tol = 0.0)
        @test_throws ArgumentError RationalScattering(data; tol = Inf)
        @test_throws ArgumentError RationalScattering(data; tol = 1e-3, minpoles = 0)
        @test_throws ArgumentError RationalScattering(data; tol = 1e-3, noisefloor = 0.0)
        @test_throws ArgumentError RationalScattering(data; tol = 1e-3, minpoles = 9, maxpoles = 4)
        # The ceiling on the search comes from the data when it is not
        # given: the rank of the Loewner pencil is the degree the samples
        # determine, which bounds the poles because each takes at least
        # one state. This RLC has three poles of rank one residues, so its
        # degree is at least three and the ceiling admits its fit. A
        # noisefloor so loose that it counts nothing still leaves a floor
        # of one rather than a ceiling of none.
        @test JC.supporteddegree(hb.S, 2pi .* fs, 1e-12) >= 3
        @test JC.supporteddegree(hb.S, 2pi .* fs, 0.5) >= 1
        # the rank is a property of the block, not of how finely it was
        # sampled, so thinning the samples does not change it
        @test JC.supporteddegree(hb.S[:, :, 1:2:end], 2pi .* fs[1:2:end], 1e-12) ==
              JC.supporteddegree(hb.S, 2pi .* fs, 1e-12)
        # a tolerance the data cannot support names the degree it carries
        cannot = try
            RationalScattering(data; tol = 1e-16); ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("determine a degree of about", cannot)
        # and it names why the orders which could not be fitted failed,
        # since that is a different problem from a tolerance too tight
        @test occursin("could not be fitted at all", cannot) ||
              occursin("closest was", cannot)
        # A block may state an active value at zero frequency, as it may
        # be active at any other, so long as it declares its own noise.
        # No fit of one can be made passive, and the search says which of
        # the orders it tried failed and why rather than only that none
        # met the tolerance.
        Nz = zeros(ComplexF64, 2, 2, length(fs))
        for k in eachindex(fs)
            Nz[:, :, k] .= 0.5*Matrix(I, 2, 2)
        end
        active = ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0,
            noise = JC.NoiseCovariance((2pi .* fs, Nz)),
            dcmodel = JC.ScatteringDC([1.4 0.0; 0.0 1.4]))
        @test_throws ArgumentError RationalScattering(active, 4)
        why = try
            RationalScattering(active; tol = 1e-3, minpoles = 2, maxpoles = 6); ""
        catch e
            sprint(showerror, e)
        end
        @test occursin("could not be fitted at all", why)
        # and without a declared noise the same statement is refused by
        # the block itself, before any fitting
        @test_throws ArgumentError ScatteringParameters((2pi .* fs, hb.S); nports = 2,
            zref = 50.0, dcmodel = JC.ScatteringDC([1.4 0.0; 0.0 1.4]))
        # A feedthrough over one which the stated value cannot be
        # contracted toward: both are positive here, so every step along
        # the path moves away from one rather than toward it, and there is
        # no step. The same feedthrough with the opposite sign is repaired
        # exactly, which is the case above.
        @test_throws ArgumentError JC.enforcepassivity(fill(-1.0, 1, 1), ones(1, 1),
            fill(0.1, 1, 1), fill(1.001, 1, 1), [0.1, 1.0]; dc = ones(1, 1))
        # The constant term is the model at infinite frequency, which no
        # perturbation over a band of frequencies can reach, so a fit whose
        # constant term is active cannot be enforced passive at all. It is
        # brought under one before the residues are fitted, and the
        # residues then fitted to what is left.
        @test JC.passiveconstant([0.5 0.0; 0.0 0.25])[2] == false
        @test JC.passiveconstant([0.5 0.0; 0.0 0.25])[1] == [0.5 0.0; 0.0 0.25]
        # a block which is lossless at infinite frequency has a unitary
        # constant term by right, and moving it would put an error into a
        # fit which was exact
        @test JC.passiveconstant([0.0 1.0; 1.0 0.0])[2] == false
        @test JC.passiveconstant(Matrix(1.0I, 3, 3))[2] == false
        # one genuinely above unity is brought down, its singular vectors
        # left alone
        let (Dc, hit) = JC.passiveconstant([4.0 0.0; 0.0 0.5])
            @test hit
            @test opnorm(Dc) <= 1
            @test svdvals(Dc) ≈ [1 - 1e-6, 0.5] rtol=1e-9
            @test svd(Dc).U ≈ svd([4.0 0.0; 0.0 0.5]).U rtol=1e-9
        end
        # the residue fit takes a fixed constant, fitting only the
        # strictly proper part to what is left
        let ws2 = 2pi .* collect(range(1e9, 5e9; length = 20)),
            pol = ComplexF64[-1e10 + 3e10im, -1e10 - 3e10im],
            Sd = zeros(ComplexF64, 1, 1, 20)
            for (k, w) in enumerate(ws2)
                Sd[1,1,k] = 0.3 + 2e10/(im*w - pol[1]) + 2e10/(im*w - pol[2])
            end
            r1, d1 = JC.fitresidues(Sd, ws2, pol)
            @test d1[1,1] ≈ 0.3 rtol=1e-8
            r2, d2 = JC.fitresidues(Sd, ws2, pol; constant = fill(0.1, 1, 1))
            @test d2[1,1] == 0.1
            # the strictly proper part absorbs the difference
            v1 = [d1[1,1] + sum(real(r1[1,1,q]/(im*w - pol[q])) for q in 1:2) for w in ws2]
            v2 = [d2[1,1] + sum(real(r2[1,1,q]/(im*w - pol[q])) for q in 1:2) for w in ws2]
            @test maximum(abs, v1 .- v2) < 0.25
        end

        # Every fit which is returned is passive over the whole imaginary
        # axis and not merely at the samples, though to the accuracy of
        # the norm search and not exactly: it decides on the largest
        # singular value that search evaluated, which is a lower bound,
        # and the level the search ends at stands `2 rtol` above it. The
        # analytic case among the norm tests above pins that gap. A
        # block whose largest singular value exceeds one is active, and a
        # transient built on it grows without bound, so a nearly passive
        # fit is scaled the rest of the way rather than returned as it is.
        # The near lossless case is the one which needs it: there the
        # enforcement's own level of `1 + atol/2` sits above the block's
        # entire dissipation. Whether a given fit is contracted, which
        # warns, or refused, which throws into the catch below, turns on
        # the roundoff of the norm search, so the logs are captured
        # without requiring a warning.
        @test_logs match_mode = :any for (rr, nps) in (
            (Resistor(120.0), (2, 4, 8)), (Resistor(1e9), (2, 4, 8)))
            lossless = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l1, 1, 2, Inductor(2e-9)),
                (:c, 2, 0, Capacitor(0.3e-12)), (:r, 2, 0, rr), (:l2, 2, 3, Inductor(2e-9)),
                (:p2, 3, 0, Port(2; Z0 = 50.0))])
            gs = collect(range(1e9, 12e9; length = 60))
            hbl = hblinsolve(2pi .* gs, lossless; keyedarrays = false)
            dat = ScatteringParameters((2pi .* gs, hbl.S); nports = 2, zref = 50.0)
            for np in nps
                f = try
                    RationalScattering(dat, np)
                catch e
                    # a fit far from passive is refused rather than scaled
                    # into a block which transmits nothing
                    @test e isa ArgumentError
                    continue
                end
                q = f.provider
                lower, _, level = JC.hinfnorm(q.A, q.B, q.C, q.D)
                @test lower <= 1
                # and the level, which is what termination establishes, is
                # over one by no more than the search's own tolerance
                @test level <= 1 + 4e-8
                @test all(real.(eigvals(q.A)) .< 0)
            end
        end
        # A network which is a perfect open at one port and a perfect
        # short at the other at infinite frequency fits with its
        # feedthrough exactly on the unit circle. The stamps snap the
        # roundoff residues of I - S and I + S to the exact zeros the
        # algebra has, so the endpoint's rate system sees zero rows
        # rather than equations of machine epsilon; without the snap the
        # reading amplifies the residual by their inverse at every step
        # and the state overflows within tens of steps. The fitted block
        # in front of a junction is checked against the same circuit as
        # lumped elements.
        embed = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:le, 1, 2, Inductor(1e-9)), (:re, 2, 3, Resistor(0.5)),
            (:ce, 3, 0, Capacitor(1e-12)), (:p2, 3, 0, Port(2; Z0 = 50.0))])
        ghz = collect(range(0.5e9, 10e9; length = 80))
        hbe = hblinsolve(2pi .* ghz, embed; keyedarrays = false)
        fite = @test_logs match_mode = :any RationalScattering(
            ScatteringParameters((2pi .* ghz, hbe.S); nports = 2, zref = 50.0), 2)
        @test maximum(abs.(abs.(diag(fite.provider.D)) .- 1)) < 1e-9
        edrive(t) = 0.1e-6*sin(2pi*4e9*t)*(1 - exp(-t/0.5e-9))
        eblock = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:cp, 1, 0, Capacitor(50e-15)), (:blk, 1, 2, fite),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(1e-12))])
        elump = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:cp, 1, 0, Capacitor(50e-15)),
            (:le, 1, 2, Inductor(1e-9)), (:re, 2, 3, Resistor(0.5)),
            (:ce, 3, 0, Capacitor(1e-12)),
            (:jj, 3, 0, JosephsonJunction(1e-9)), (:cj, 3, 0, Capacitor(1e-12))])
        se = [transientsolve(transientproblem(c; sources = [TransientSource(1, edrive)]),
            (0.0, 2e-9); dt = 1e-12, method = GaussLegendre()) for c in (eblock, elump)]
        @test maximum(abs, se[1].voltage .- se[2].voltage) <
            1e-3*maximum(abs, se[2].voltage)
        # a fit needs its last pole unless a constant reproduces the data:
        # one real pole and one conjugate pair fitted at their order and
        # above keep it, and constant data is refused as a rational block
        onepole = RationalScattering(fill(-1.0, 1, 1), ones(1, 1), fill(0.5, 1, 1), zeros(1, 1))
        for np in (1, 2, 3, 4)
            f1 = RationalScattering(onepole, np; frequencies = range(0.01, 1.0; length = 41))
            @test size(f1.provider.A) == (1, 1)
            @test f1.provider.A[1, 1] ≈ -1.0 rtol=1e-8
        end
        pair = RationalScattering([-0.1 1.0; -1.0 -0.1], [0.0; 1.0;;], [0.05 0.0], zeros(1, 1))
        for np in (2, 3, 5)
            f2 = RationalScattering(pair, np; frequencies = range(0.05, 3.0; length = 121))
            @test size(f2.provider.A) == (2, 2)
            @test sort(eigvals(f2.provider.A); by = imag) ≈ [-0.1 - im, -0.1 + im] rtol=1e-8
        end
        @test_throws ArgumentError RationalScattering(ScatteringParameters(fill(0.3, 1, 1)), 2; frequencies = range(0.1, 1.0; length = 21))
        # a three port, whose nine entries are eliminated one by one in the
        # fit, so nothing the size of every entry's every sample is formed
        tee = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:r1, 1, 4, Resistor(20.0)), (:p2, 2, 0, Port(2; Z0 = 50.0)), (:r2, 2, 4, Resistor(20.0)),
            (:p3, 3, 0, Port(3; Z0 = 50.0)), (:l3, 3, 5, Inductor(2e-9)), (:c3, 5, 4, Capacitor(0.4e-12)), (:r4, 4, 0, Resistor(80.0)), (:c4, 4, 0, Capacitor(0.2e-12))])
        hb3 = hblinsolve(2pi .* fs, tee; keyedarrays = false)
        fit3 = RationalScattering(ScatteringParameters((2pi .* fs, hb3.S); nports = 3, zref = 50.0), 8)
        S3 = zeros(ComplexF64, 3, 3, length(fs))
        JC.evaluateprovider!(S3, fit3.provider, 2pi .* fs)
        @test maximum(abs.(S3 .- hb3.S)) < 1e-10 && size(fit3.provider.A, 1) <= 6
        @test_throws ArgumentError RationalScattering(data, 4; frequencies = fs[1:4])
        drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(t/1e-9)^2*sinpi(2*3e9*t)
        withblock = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:blk, 1, 2, fitted),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.3e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        explicit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:l1, 1, 4, Inductor(1.5e-9)),
            (:c, 4, 0, Capacitor(0.6e-12)), (:r, 4, 0, Resistor(120.0)), (:l2, 4, 2, Inductor(1.0e-9)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.3e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        pb = transientproblem(withblock; sources = [TransientSource(1, drive)])
        pe = transientproblem(explicit; sources = [TransientSource(1, drive)])
        sb = transientsolve(pb, (0.0, 1.5e-9); dt = 2e-12, method = GaussLegendre(), rtol = 1e-12)
        se = transientsolve(pe, (0.0, 1.5e-9); dt = 2e-12, method = GaussLegendre(), rtol = 1e-12)
        @test sb.outgoing ≈ se.outgoing rtol=1e-8
    end

    @testset "ideal transmission lines" begin
        # the line is the method of characteristics: a conductance at its
        # impedance and the far port's wave a delay earlier as a current,
        # nothing inverted, so the half wave resonance that makes its
        # admittance singular in frequency is only the recursion through
        # the history; the scattering parameters of a mismatched line
        # between capacitive loads are the linearized solver's, at that
        # resonance too, an open end echoes with the sign of an open and a
        # short with the sign of a short after two delays, a matched line
        # passes a pulse's energy, and the transmitted pulse converges at
        # the rule's order through the six point history read
        tau = 0.3e-9
        line = TransmissionLine(60.0, tau*3e8; vp = 3e8)
        loaded = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:line, 1, 2, line),
            (:c2, 2, 0, Capacitor(0.4e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        pl = transientproblem(loaded; sources = [TransientSource(1, t -> 0.0)])
        @test length(pl.lines) == 1 && isempty(pl.blocks) && pl.lines[1].delay ≈ tau
        for f in (1.5e9, 1/(2tau), 4.1e9)
            a = 1e-7
            pf = transientproblem(loaded; sources = [TransientSource(1, t -> t <= 0 ? 0.0 : 2a*(1 - exp(-(t/2e-9)^2))*cospi(2f*t))])
            dt = 1/(80f)
            T = 40e-9
            nsteps = round(Int, T/dt)
            sol = transientsolve(pf, (0.0, nsteps*dt); dt, method = GaussLegendre())
            window(t) = t >= T - 8/f ? sinpi((t - (T - 8/f))*f/8)^2 : 0.0
            inc = transientdemodulate(sol, 1, f; quantity = :incident, window)
            hb = hblinsolve(2pi*[f], loaded; keyedarrays = false)
            @test transientdemodulate(sol, 1, f; quantity = :outgoing, window)/inc ≈ hb.S[1, 1, 1] atol=1e-5
            @test transientdemodulate(sol, 2, f; quantity = :outgoing, window)/inc ≈ hb.S[2, 1, 1] atol=1e-5
        end
        pulse(t) = 2e-6*exp(-((t - 0.1e-9)/0.02e-9)^2)
        open = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(50.0, tau*3e8; vp = 3e8)),
            (:p2, 2, 0, Port(2; Z0 = 50.0, termination = nothing))])
        po = transientproblem(open; sources = [TransientSource(1, pulse)])
        @test isempty(po.gaugeindices) && po.inertialess == [[1], [2]] && isempty(po.algebraic)
        so = transientsolve(po, (0.0, 1.2e-9); dt = 2e-12, method = GaussLegendre())
        peak = argmax(abs.(so.outgoing[1, :]))
        @test so.times[peak] ≈ 0.1e-9 + 2tau atol=1e-12
        @test so.outgoing[1, peak] ≈ maximum(so.incident[1, :]) rtol=1e-3
        @test maximum(abs.(so.outgoing[1, so.times .< 0.1e-9 + 2tau - 0.1e-9])) < 1e-15
        short = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 0, TransmissionLine(50.0, tau*3e8; vp = 3e8))])
        ss = transientsolve(transientproblem(short; sources = [TransientSource(1, pulse)]), (0.0, 1.2e-9); dt = 2e-12, method = GaussLegendre())
        speak = argmax(abs.(ss.outgoing[1, :]))
        @test ss.times[speak] ≈ 0.1e-9 + 2tau atol=1e-12
        @test ss.outgoing[1, speak] ≈ -maximum(ss.incident[1, :]) rtol=1e-3
        matched = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(50.0, tau*3e8; vp = 3e8)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        pm = transientproblem(matched; sources = [TransientSource(1, pulse)])
        sm = transientsolve(pm, (0.0, 1e-9); dt = 2e-12, method = GaussLegendre())
        @test sum(abs2, sm.outgoing[2, :]) ≈ sum(abs2, sm.incident[1, :]) rtol=1e-6
        @test maximum(abs.(sm.outgoing[1, :])) < 1e-15
        @test sm.times[argmax(abs.(sm.outgoing[2, :]))] ≈ 0.1e-9 + tau atol=1e-12
        @test size(sm.linewaves) == (2, 501)
        @test_throws ArgumentError transientsolve(pm, (0.0, 1e-9); dt = 0.5e-9, method = GaussLegendre())
        @test_throws ArgumentError transientsolve(pm, (0.0, 1e-9); dt = 2e-12)
        @test_throws ArgumentError transientproblem(Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(50.0, 0.0)),
            (:p2, 2, 0, Port(2; Z0 = 50.0))]))
        # the order of the transmitted pulse through a fractional delay (an
        # integer one reads samples exactly), and the batch and the checkpoints
        fractional = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(50.0, tau*3e8)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        smooth = transientproblem(fractional; sources = [TransientSource(1, t -> 2e-6*exp(-((t - 0.45e-9)/0.06e-9)^2))])
        ref = transientsolve(smooth, (0.0, 5tau); dt = tau/2000, method = GaussLegendre())
        errs = [maximum(abs.(transientsolve(smooth, (0.0, 5tau); dt = tau/m, method = GaussLegendre()).outgoing[2, :] .-
            ref.outgoing[2, 1:(2000 ÷ m):end])) for m in (50, 100)]
        @test errs[1]/errs[2] > 14
        half = transientproblem(pm; sources = [TransientSource(1, t -> pulse(t)/2)])
        b = transientsolve([pm, half], (0.0, 1e-9); dt = 2e-12, record = :phases)
        @test b.outgoing[:, :, 1] == sm.outgoing
        @test b.linewaves[:, :, 2] ≈ sm.linewaves ./ 2 rtol=1e-12
        # with checkpoints the record of the waves gives way to the history
        # before each checkpoint when that is smaller, its columns before
        # the start the initial waves and after it the record's; here the
        # delay exceeds the checkpoint interval, so the record is kept and
        # the replay reads it, and with a short delay the tails are kept
        cps = transientsolve(pm, (0.0, 1e-9); dt = 2e-12, method = GaussLegendre(), record = :checkpoints, checkpointevery = 50)
        npre = JC.lineprehistory(pm, 2e-12)
        @test npre*10 > 501 && size(cps.checkpoints.waves, 2) == 0
        @test cps.linewaves ≈ sm.linewaves rtol=1e-12
        shortline = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(50.0, 10e-12; vp = 1.0)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        psh = transientproblem(shortline; sources = [TransientSource(1, pulse)])
        ssh = transientsolve(psh, (0.0, 1e-9); dt = 2e-12, method = GaussLegendre())
        csh = transientsolve(psh, (0.0, 1e-9); dt = 2e-12, method = GaussLegendre(), record = :checkpoints, checkpointevery = 50)
        npres = JC.lineprehistory(psh, 2e-12)
        @test isnothing(csh.linewaves) && size(csh.checkpoints.waves) == (2, npres, 10)
        for c in 1:10, j in 1:npres
            column = (c - 1)*50 + j
            @test csh.checkpoints.waves[:, j, c] ≈ (column >= npres ? ssh.linewaves[:, column - npres + 1] : ssh.linewaves[:, 1]) rtol=1e-12
        end
        # both replays give the phases record's responses
        for (prob, cp) in ((pm, cps), (psh, csh))
            ph = transientsolve(prob, (0.0, 1e-9); dt = 2e-12, method = GaussLegendre(), record = :phases)
            cur = [1e-8*sinpi(2*1.3e9*t)^2 for k in 1:1, t in ph.times]
            wts = [cospi(2*1.7e9*t) for k in 1:2, t in ph.times]
            @test transienttangent(cp, cur; targets = [1]).outgoing ≈ transienttangent(ph, cur; targets = [1]).outgoing rtol=1e-9
            @test transientadjoint(cp, wts).currents ≈ transientadjoint(ph, wts).currents rtol=1e-9
        end
        # A short line between mismatched loads: a delay of 1.2 steps
        # leaves one accepted sample past every read, so the read is the
        # linear stencil, a contraction like the cubic and the quintic; the
        # pulse rings down instead of growing, as it did through a one sided
        # stencil of gain 2.4. The gain of every stencil, over delays,
        # stages, endpoints and history lengths, is at most one, and the
        # recurrence of a read through the reflection 9/11 is stable.
        short = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(5.0, 1.2e-12; vp = 1.0)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        sp = transientproblem(short; sources = [TransientSource(1, t -> 0 < t < 4e-12 ? 1e-6*sinpi(t/4e-12)^2 : 0.0)])
        ss = transientsolve(sp, (0.0, 200e-12); dt = 1e-12, method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
        early, late = maximum(abs, ss.voltage[:, ss.times .<= 10e-12]), maximum(abs, ss.voltage[:, ss.times .>= 150e-12])
        @test late < 1e-9*early
        worst = 0.0
        for delay in (1.0, 1.05, 1.2, 1.8, 2.5, 3.5, 4.0, 10.5), c in (JC.gausscoefficients().c..., 1.0), accepted in (2, 3, 4, 5, 6, 21)
            first, weights = JC.linestencil(accepted - 1 + c - delay, 0.0, 1.0, accepted)
            @test length(weights) in (1, 2, 4, 6) && first >= 0 && first + length(weights) <= accepted
            worst = max(worst, maximum(abs(sum(weights[k]*cis(-w*(first + k - 1)) for k in eachindex(weights))) for w in range(0, pi; length = 2001)))
        end
        @test worst <= 1 + 1e-12
        @test JC.linestencil(21 - 10.5, 0.0, 1.0, 21)[2] ≈ [3, -25, 150, 150, -25, 3] ./ 256 rtol=1e-12
        radius = 0.0
        for delay in (1.0, 1.2, 1.8, 2.5, 3.5), sign in (-1, 1)
            first, weights = JC.linestencil(21.0 - delay, 0.0, 1.0, 21)
            lag = 21 .- (first .+ (0:length(weights) - 1))
            M = zeros(maximum(lag), maximum(lag))
            for k in eachindex(weights); M[1, lag[k]] += sign*(9/11)*weights[k]; end
            for k in 2:maximum(lag); M[k, k - 1] = 1.0; end
            radius = max(radius, maximum(abs, eigvals(M)))
        end
        @test radius < 1
        # the tangent of a quiet matched line along a drive given at the
        # stages is the linear circuit's forward solve under that drive:
        # the same history and stencils from the start, the prehistory
        # constant before it
        tenhalf = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(50.0, 10.5e-12; vp = 1.0)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        quiet = transientproblem(tenhalf; sources = [TransientSource(1, t -> 0.0)])
        sine(t) = t <= 0 ? 0.0 : 1e-6*sinpi(2*6e9*t)
        driven = transientproblem(quiet; sources = [TransientSource(1, sine)])
        qb = transientsolve(quiet, (0.0, 100e-12); dt = 1e-12, method = GaussLegendre(), record = :phases, rtol = 1e-12, atol = 1e-13)
        qf = transientsolve(driven, (0.0, 100e-12); dt = 1e-12, method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
        staged = zeros(1, 3, length(qb.times), 1)
        for (k, t) in enumerate(qb.times), (i, c) in enumerate((0.0, JC.gausscoefficients().c...))
            staged[1, i, k, 1] = sine(t + c*1e-12)
        end
        qt = transienttangent(qb, staged; targets = [1])
        @test maximum(abs, qf.outgoing[2, :]) > 1e-8
        @test qt.outgoing[:, :, 1] ≈ qf.outgoing rtol=1e-10
        # a batch's conditions each start from their own prehistory: two
        # direct currents on the line, each an equilibrium, stay so in the
        # batch as alone (a repeat of the initial waves once interleaved
        # them)
        dc1 = transientproblem(tenhalf; sources = [TransientSource(1, 1e-6)])
        dc2 = transientproblem(dc1; sources = [TransientSource(1, 2e-6)])
        eq = [transientstate(dc1; voltage = [25e-6, 25e-6], linecurrents = [0.5e-6]), transientstate(dc2; voltage = [50e-6, 50e-6], linecurrents = [1e-6])]
        for rec in (:phases, :checkpoints)
            bdc = transientsolve([dc1, dc2], (0.0, 100e-12); dt = 1e-12, method = GaussLegendre(), record = rec, checkpointevery = 10,
                initialstate = eq, rtol = 1e-12, atol = 1e-13)
            for (j, prob) in enumerate((dc1, dc2))
                alone = transientsolve(prob, (0.0, 100e-12); dt = 1e-12, method = GaussLegendre(), record = rec, checkpointevery = 10,
                    initialstate = eq[j], rtol = 1e-12, atol = 1e-13)
                @test bdc.voltage[:, :, j] ≈ alone.voltage rtol=1e-12
                @test maximum(abs.(alone.voltage .- 25e-6*j)) < 1e-15
            end
        end
        # a state with a direct current on the line, and the waves it carries
        state = transientstate(pm; linecurrents = [1e-6])
        @test state.waves ≈ [sqrt(50)*1e-6/2, -sqrt(50)*1e-6/2]
        # the tangent and the adjoint through the histories: a mismatched
        # line before a pumped junction, the tangent against finite
        # differences, the adjoint against the tangent, the cotangent of
        # the prehistory against a tangent of it, checkpoints and a batch
        cable = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(60.0, 0.09)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.4e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        cdrive(t) = t <= 0 ? 0.0 : 0.2e-6*sinpi(t/1e-9)^2*sinpi(2*3e9*t)
        cprobe(t, k) = 1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k)
        cT, cdt = 1.5e-9, 2e-12
        cp = transientproblem(cable; sources = [TransientSource(1, cdrive), TransientSource(2, t -> 0.0)])
        crec = transientsolve(cp, (0.0, cT); dt = cdt, method = GaussLegendre(), rtol = 1e-12, record = :phases)
        ccurrents = [cprobe(t, k) for k in 1:2, t in crec.times]
        ctg = transienttangent(crec, ccurrents)
        ceps = 1e-4
        cshift(s) = transientproblem(cable; sources = [TransientSource(1, t -> cdrive(t) + s*ceps*cprobe(t, 1)), TransientSource(2, t -> s*ceps*cprobe(t, 2))])
        csp = transientsolve(cshift(1), (0.0, cT); dt = cdt, method = GaussLegendre(), rtol = 1e-12)
        csm = transientsolve(cshift(-1), (0.0, cT); dt = cdt, method = GaussLegendre(), rtol = 1e-12)
        @test ctg.outgoing ≈ (csp.outgoing .- csm.outgoing) ./ (2ceps) rtol=1e-5
        cweights = [cospi(2*1.7e9*t + k) for k in 1:2, t in crec.times]
        cad = transientadjoint(crec, cweights)
        @test sum(cweights .* ctg.outgoing) ≈ sum(cad.currents .* ccurrents) rtol=1e-9
        @test size(cad.initialwaves, 1) == 2 && size(cad.initialwaves, 2) == JC.lineprehistory(cp, cdt)
        dw = zeros(2, size(cad.initialwaves, 2), 1)
        dw[1, :, 1] .= 1e-9
        ctw = transienttangent(crec, zeros(2, length(crec.times)); initialstate = (zeros(2), zeros(2), dw))
        @test sum(cweights .* ctw.outgoing) ≈ sum(cad.initialwaves .* dw[:, :, 1]) rtol=1e-9
        ccps = transientsolve(cp, (0.0, cT); dt = cdt, method = GaussLegendre(), rtol = 1e-12, record = :checkpoints, checkpointevery = 50)
        @test transienttangent(ccps, ccurrents).outgoing ≈ ctg.outgoing rtol=1e-9
        @test transientadjoint(ccps, cweights).currents ≈ cad.currents rtol=1e-9
        chalf = transientproblem(cp; sources = [TransientSource(1, t -> cdrive(t)/2), TransientSource(2, t -> 0.0)])
        cb = transientsolve([cp, chalf], (0.0, cT); dt = cdt, rtol = 1e-12, record = :phases)
        @test transienttangent(cb, ccurrents).outgoing[:, :, 1] == ctg.outgoing
        @test transientadjoint(cb, cweights).currents[:, :, 1] == cad.currents
    end
end
