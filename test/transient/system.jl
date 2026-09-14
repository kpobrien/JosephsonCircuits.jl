using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Random
using Test

# The scattering blocks and the transmission lines in time: constant
# blocks, rational blocks and ideal lines with their delays, on the
# stepping, the tangent and the adjoint. The fit which builds a rational
# block is in circuit/vectorfit.jl and the core of the solver in
# transient/solve.jl; the files are separate so that they run on
# different workers.
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

# A pumped device as a block in time: the fit of its harmonic transfer
# functions to filters whose outputs the step modulates at the harmonics
# of the pump. The references are the tabulated block itself in harmonic
# balance, the device as junctions in time, and the linearized solve of
# the fitted block for the noise, which runs the tangent and the adjoint
# through the modulated coupling.
@testset "a pumped block realized in time" begin
    JC = JosephsonCircuits
    Z0 = 50.0
    jpa = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)), (:cc, 1, 2, Capacitor(100.0e-15)),
        (:jj, 2, 0, JosephsonJunction(1000.0e-12)), (:cj, 2, 0, Capacitor(1000.0e-15))])
    fp, ip = 4.75e9, 0.00565e-6
    ws = 2pi*collect(range(4.0e9, 5.5e9; length = 76))
    sol = hbsolve(ws, (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (4,), (8,), jpa; ftol = 1e-14)
    blk = LinearizedScattering(sol.linearized, 2pi*fp)
    # an unfitted block has no realization in time
    @test_throws ArgumentError transientproblem(Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)), (:b, 1, blk)]))
    fit = RationalScattering(blk, 10)
    @test fit.harmonics == blk.harmonics && fit.phase == 0.0
    # the fit is not lossless, and states the noise its commutator requires
    @test fit.noise isa NoiseCovariance && fit.noise.completed && fit.atol == blk.atol
    @test fit.providers[1] isa JC.RationalScatteringProvider
    @test all(p -> p isa JC.ModulatedRationalProvider, fit.providers[2:end])
    @test all(p -> maximum(abs, p.cosine.D) == 0 && maximum(abs, p.sine.D) == 0, fit.providers[2:end])
    @test_throws ArgumentError RationalScattering(blk, 0)
    # the fit against the tables it came from, at their knots
    for (j, k) in enumerate(blk.harmonics)
        nus = JC.piecewisefrequencies(blk.providers[j])
        Ht = zeros(ComplexF64, 1, 1, length(blk.harmonics), length(nus))
        Hf = similar(Ht)
        JC.evaluateharmonics!(Ht, blk, nus)
        JC.evaluateharmonics!(Hf, fit, nus)
        @test maximum(abs, Hf[1, 1, j, :] .- Ht[1, 1, j, :]) < 1e-4*max(1.0, maximum(abs, Ht[1, 1, j, :]))
    end
    # the fitted block in harmonic balance against the junctions
    c = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)), (:b, 1, fit)])
    hb = hbsolve(ws, (2pi*fp,), [], (4,), (8,), c)
    @test maximum(abs, Array(hb.linearized.S) .- Array(sol.linearized.S)) < 1e-4
    # a weak signal through the block in time, demodulated at the signal
    # and the idler after the block has settled, against the scattering
    # matrix of the fitted block in power waves, and against the
    # junctions in time with their pump ramped on
    fs, Is, dt = 4.6e9, 1e-9, 2.5e-12
    fi = 2fp - fs
    a0 = Is*sqrt(Z0)/2
    tsol = transientsolve(transientproblem(c; sources = [TransientSource(1, t -> Is*sinpi(2fs*t))]), (0.0, 80e-9);
        dt, method = GaussLegendre())
    plan = transientiqplan(tsol.times, [fs, fi]; duration = 4/(fp - fs), ports = [1, 1], stride = 400)
    iq = transientiq(plan, tsol.outgoing)
    k = argmin(abs.(ws .- 2pi*fs))
    Ss = hb.linearized.S((0,), 1, (0,), 1, k)
    Si = hb.linearized.S((-2,), 1, (0,), 1, k)
    @test isapprox(abs(iq[1, end])/a0, abs(Ss); rtol = 1e-3)
    @test isapprox(abs(iq[2, end])/a0, abs(Si)*sqrt(fi/fs); rtol = 1e-2)
    ramp(t) = t <= 0 ? 0.0 : t >= 10e-9 ? 1.0 : (1 - cospi(t/10e-9))/2
    jsol = transientsolve(transientproblem(jpa; sources = [TransientSource(1, t -> 2ip*ramp(t)*cospi(2fp*t) + Is*sinpi(2fs*t))]),
        (0.0, 80e-9); dt, method = GaussLegendre())
    jiq = transientiq(plan, jsol.outgoing)
    @test isapprox(abs(jiq[1, end]), abs(iq[1, end]); rtol = 1e-2)
    @test isapprox(abs(jiq[2, end]), abs(iq[2, end]); rtol = 1e-2)
    @test abs(angle(jiq[2, end]/iq[2, end])) < 0.02
    # the noise of the block with its conversion ramped on, the tangent
    # and the adjoint through the modulated coupling: the gain and the
    # quantum efficiency of the linearized solve, and the two methods
    # agree with each other
    ramped = RationalScattering(LinearizedScattering(sol.linearized, 2pi*fp; envelope = ramp), 10)
    cr = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)), (:b, 1, ramped)])
    fn = 4.7e9
    hbn = hbsolve([2pi*fn], (2pi*fp,), [], (4,), (8,), cr)
    # a settling time which is not a whole number of pump periods, so
    # that the reference of the bath quadratures, the start of the
    # record, is not a zero of the pump's phase
    settle, record = 60.0625e-9, 20e-9
    nsol = transientsolve(transientproblem(cr), (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
    first = round(Int, settle/dt) + 1
    nplan = transientquantumplan(nsol.times[first:end], [fn])
    nfreqs = sort!(abs.([fn + 2m*fp for m in -2:2]))
    noise = transientnoise(nsol, nplan; frequencies = nfreqs, weights = fill(1/record, 5), inputs = nplan, commutationrtol = 3e-3)
    @test noise.diagnostics.passed
    metrics = transientquantumefficiency(noise.gain, noise.covariance; rtol = 3e-3)
    @test isapprox(metrics.gain, abs2(hbn.linearized.S((0,), 1, (0,), 1, 1)); rtol = 1e-4)
    @test isapprox(metrics.QE, hbn.linearized.QE((0,), 1, (0,), 1, 1); rtol = 1e-4)
    forward = transientnoise(nsol, nplan; frequencies = nfreqs, weights = fill(1/record, 5), inputs = nplan, method = :forward, commutationrtol = 3e-3)
    @test forward.covariance ≈ noise.covariance rtol=1e-10
    @test forward.gain ≈ noise.gain rtol=1e-10
    # a device with loss states its noise, which in time is the pair
    # terms of its group: the bath frequencies a harmonic apart, and
    # those summing to one, the signal and the idler, are correlated as
    # the block's harmonic covariances say, and the quantum efficiency is
    # that of the linearized solve of the fitted block, and of the
    # junctions
    lossy = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)), (:cc, 1, 2, Capacitor(100.0e-15)),
        (:jj, 2, 0, JosephsonJunction(1000.0e-12)), (:cj, 2, 0, Capacitor(1000.0e-15)), (:r, 2, 0, Resistor(2.0e4))])
    soll = hbsolve(ws, (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (4,), (8,), lossy; ftol = 1e-14, returnCnoise = true)
    stated = RationalScattering(LinearizedScattering(soll.linearized, 2pi*fp; noise = NoiseCovariance(soll.linearized.Cnoise), envelope = ramp), 10)
    @test stated.noise isa NoiseCovariance && stated.noise.completed && stated.atol == 1e-6 && stated.noise.atol == 1e-8
    cl = Circuit([(:p1, 1, 0, Port(1; Z0 = Z0)), (:b, 1, stated)])
    hbl = hbsolve([2pi*fn], (2pi*fp,), [], (4,), (8,), cl; returnCnoise = true)
    bathsl = transientnoisebaths(transientproblem(cl))
    @test length(bathsl) == 2 && length(bathsl.groups) == 1 && bathsl.channels[2].temperature == 0.0
    lsol = transientsolve(transientproblem(cl), (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
    lplan = transientquantumplan(lsol.times[first:end], [fn])
    lnoise = transientnoise(lsol, lplan; frequencies = nfreqs, weights = fill(1/record, 5), inputs = lplan, commutationrtol = 3e-3)
    @test lnoise.diagnostics.passed
    lmetrics = transientquantumefficiency(lnoise.gain, lnoise.covariance; rtol = 3e-3)
    @test isapprox(lmetrics.gain, abs2(hbl.linearized.S((0,), 1, (0,), 1, 1)); rtol = 1e-4)
    @test isapprox(lmetrics.QE, hbl.linearized.QE((0,), 1, (0,), 1, 1); rtol = 1e-4)
    @test isapprox(lmetrics.QE, hbsolve([2pi*fn], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (4,), (8,), lossy;
        ftol = 1e-14).linearized.QE((0,), 1, (0,), 1, 1); rtol = 1e-4)
end

@testset "the checks and the stage solve of a pumped block in time" begin
    JC = JosephsonCircuits
    Z0 = 50.0
    wp = 2pi*1e9
    zero1 = zeros(ComplexF64, 1, 1)
    # a stable one state block converting by one harmonic, without a fit
    function model(; amp = 0.1, direct = 0.0, noise = Lossless(), envelope = nothing)
        p0 = JC.RationalScatteringProvider(zeros(0, 0), zeros(0, 1), zeros(1, 0), fill(direct, 1, 1))
        pc = JC.RationalScatteringProvider(fill(-wp, 1, 1), fill(wp, 1, 1), fill(amp, 1, 1), zeros(1, 1))
        pz = JC.RationalScatteringProvider(zeros(0, 0), zeros(0, 1), zeros(1, 0), zeros(1, 1))
        return LinearizedScattering([p0, JC.ModulatedRationalProvider(pc, pz)], wp; harmonics = [0, 1], nports = 1, zref = Z0, noise, envelope)
    end
    one(b, rest...) = Circuit([(:p, 1, 0, Port(1; Z0 = Z0)), (:b, 1, b), rest...])
    dt = 2e-11
    # the noise model of a pumped block is checked over the modes its
    # pair terms are read from: a covariance below the commutation
    # relations, and a declared losslessness the block does not have
    for b in (model(noise = NoiseCovariance([zero1, zero1])), model())
        sol = transientsolve(transientproblem(one(b)), (0.0, 40e-9 - dt); dt, method = GaussLegendre(), record = :checkpoints)
        plan = transientquantumplan(sol.times, [0.4e9])
        @test_throws ArgumentError transientnoise(sol, plan; frequencies = [0.4e9, 0.6e9], weights = fill(1/40e-9, 2))
    end
    # the pair terms correlate the bath frequencies, which are not split
    # over tiles: a budget too small for them is refused
    stated = model(noise = NoiseCovariance([fill(10.0, 1, 1), zero1]))
    sol = transientsolve(transientproblem(one(stated)), (0.0, 40e-9 - dt); dt, method = GaussLegendre(), record = :checkpoints)
    plan = transientquantumplan(sol.times, [0.4e9])
    JC.noisememorybudget[] = 1
    try
        @test_throws ArgumentError transientnoise(sol, plan; frequencies = [0.4e9, 0.6e9], weights = fill(1/40e-9, 2))
    finally
        JC.noisememorybudget[] = 0
    end
    # two outputs which share an input are correlated whatever harmonic
    # apart they are: the outputs at 0.4 and 2.4 GHz of a block
    # converting by one harmonic are both fed by the input at 1.4 GHz,
    # so the pair terms are read from the family of the bath frequencies
    # at once, and the whole commutator matrix of the two outputs is
    # canonical and their covariance the family's
    shared = model(amp = 0.5, noise = NoiseCovariance([fill(10.0, 1, 1), zero1]))
    ssol = transientsolve(transientproblem(one(shared)), (0.0, 20e-9 - dt); dt, method = GaussLegendre(), record = :checkpoints)
    splan = transientquantumplan(ssol.times[501:end], [0.4e9, 2.4e9]; ports = [1, 1])
    sfreqs = [0.4e9, 0.6e9, 1.4e9, 2.4e9, 3.4e9]
    sn = transientnoise(ssol, splan; frequencies = sfreqs, weights = fill(1/10e-9, 5))
    @test sn.diagnostics.passed
    # 0.4 and 2.4 GHz are two pump frequencies apart, so one ladder of
    # the bath holds both
    onladder(b, fs, f) = only(filter(L -> any(r -> abs(r - 2pi*f) < 1e-3, L.rows), JC.bathfamily(b, fs)))
    # the bath falls into the ladders of the pump: the modes at 0.4 and
    # 0.6 GHz sum to the pump frequency, so the ladder which holds the
    # positive mode of one holds the negative mode of the other, and the
    # anomalous term of that pair is an entry of its matrix
    fam = sort(JC.bathfamily(shared, sfreqs); by = L -> length(L.positive))
    @test length(fam) == 2
    @test fam[1].positive == [2] && fam[1].negative == [1, 3, 4, 5]
    @test fam[2].positive == [1, 3, 4, 5] && fam[2].negative == [2]
    # the pair terms of a two port block, against the one family over
    # every signed bath frequency they were read from before: a block
    # which converts, with a pump phase and a reference time, over
    # frequencies on independent ladders, a chain a pump apart, a signal
    # and its idler, and the degenerate pair at half the pump, whose
    # difference and whose sum are both a multiple of it
    function densepairs(blk, fs, reference)
        np = blk.nports
        rows, cols, Kd = JC.pumpedfamily(blk, vcat(2pi .* fs, -2pi .* fs); reach = false)
        Sd, Kcd, Vd = JC.pumpednoisematrices(blk, rows, cols, Kd)
        isnothing(Vd) && (Vd = zeros(ComplexF64, size(Kcd)))
        nrd = length(rows)
        at(nu) = findfirst(r -> abs(r - nu) <= 1e-9*(abs(nu) + blk.wp), rows)
        block2(A, ia, ib) = (isnothing(ia) || isnothing(ib)) ? zeros(ComplexF64, np, np) :
            ComplexF64[A[(p - 1)*nrd + ia, (q - 1)*nrd + ib] for p in 1:np, q in 1:np]
        multiple(d) = abs(d - round(d/blk.wp)*blk.wp) <= 1e-6*blk.wp
        quads = (X, Y) -> (real.(X .+ Y) ./ 2, (imag.(X) .- imag.(Y)) ./ 2, .-(imag.(X) .+ imag.(Y)) ./ 2, real.(X .- Y) ./ 2)
        out = Tuple{Int,Int,NTuple{4,Matrix{Float64}},NTuple{4,Matrix{Float64}}}[]
        for (a, fa) in enumerate(fs), (b, fb) in enumerate(fs)
            nua, nub = 2pi*fa, 2pi*fb
            normal, anomalous = multiple(nua - nub), multiple(nua + nub)
            (normal || anomalous) || continue
            ia, ib, ibm = at(nua), at(nub), at(-nub)
            zed = zeros(ComplexF64, np, np)
            N, Kn = normal ? (block2(Vd, ia, ib), block2(Kcd, ia, ib)) : (zed, zed)
            M, Km = anomalous ? (block2(Vd, ia, ibm), block2(Kcd, ia, ibm)) : (zed, zed)
            rn, ra = cis((nua - nub)*reference), cis((nua + nub)*reference)
            push!(out, (a, b, quads(rn .* N, ra .* M), quads(2im .* rn .* Kn, 2im .* ra .* Km)))
        end
        return out
    end
    q0 = JC.RationalScatteringProvider(zeros(0, 0), zeros(0, 2), zeros(2, 0), [0.3 0.1; 0.2 0.4])
    qc = JC.RationalScatteringProvider(fill(-wp, 1, 1), [0.4wp 0.1wp], reshape([0.2, 0.15], 2, 1), zeros(2, 2))
    qz = JC.RationalScatteringProvider(zeros(0, 0), zeros(0, 2), zeros(2, 0), zeros(2, 2))
    pfreqs = [0.3e9, 0.42e9, 0.5e9, 0.7e9, 1.3e9]
    for (noise, stated) in ((NoiseCovariance([[8.0 1.0; 1.0 9.0], [0.5 0.2; 0.2 0.4]]), true),
            (NoiseCovariance([zeros(2, 2), zeros(2, 2)]; completed = true), false))
        blk = LinearizedScattering([q0, JC.ModulatedRationalProvider(qc, qz)], wp; harmonics = [0, 1],
            nports = 2, zref = Z0, phase = 0.73, noise)
        for reference in (0.0, 1.7e-9)
            got = JC.pumpedpairterms(blk, pfreqs, JC.CPU(), reference)
            want = densepairs(blk, pfreqs, reference)
            @test [(t.a, t.b) for t in got] == [(w[1], w[2]) for w in want]
            @test all(zip(got, want)) do (t, w)
                all(isapprox.(t.E, w[3]; rtol = 1e-10, atol = 1e-10)) && all(isapprox.(t.C, w[4]; rtol = 1e-10, atol = 1e-10))
            end
            # the pair at half the pump is one term and not two, and it
            # carries both correlations: the covariance a harmonic apart
            # is its anomalous one, which is why its two quadratures
            # differ where the block states it
            degenerate = only(filter(t -> t.a == 3 && t.b == 3, got))
            if stated
                @test !isapprox(degenerate.E[1], degenerate.E[4]; rtol = 1e-3)
            end
        end
    end
    L4 = onladder(shared, sfreqs, 0.4e9)
    rows = L4.rows
    Sf, Kc, Vf = JC.pumpednoisematrices(shared, rows, L4.cols, L4.K)
    i4, i24 = argmin(abs.(rows .- 2pi*0.4e9)), argmin(abs.(rows .- 2pi*2.4e9))
    @test abs(Kc[i4, i24]) > 0.1
    @test isapprox(sn.covariance[1, 1] + sn.covariance[2, 2], sum(abs2, Sf[i4, :]) + real(Vf[i4, i4]); rtol = 1e-4)
    @test isapprox(norm(sn.covariance[1:2, 3:4]), abs(Vf[i4, i24] + dot(Sf[i24, :], Sf[i4, :]))/sqrt(2); rtol = 1e-4)
    # a covariance completed over the padded ladder of the modes of a
    # solve, of a block far from lossless: the transient emits the
    # completed covariance of the family of its bath frequencies, which
    # is the harmonic balance solve's over its modes whatever modes
    # either keeps, since both restrict one padded completion
    bare = model(amp = 0.5, direct = 1.0, noise = NoiseCovariance([zero1, zero1]; completed = true, padding = 8))
    bsol = transientsolve(transientproblem(one(bare)), (0.0, 20e-9 - dt); dt, method = GaussLegendre(), record = :checkpoints)
    bplan = transientquantumplan(bsol.times[501:end], [0.4e9])
    v2 = Float64[]
    for nm in (2, 4)
        bath = sort!(abs.([0.4e9 + m*1e9 for m in -nm:nm]))
        bn = transientnoise(bsol, bplan; frequencies = bath, weights = fill(1/10e-9, length(bath)))
        @test bn.diagnostics.passed
        Lb = onladder(bare, bath, 0.4e9)
        rows = Lb.rows
        Sf, _, Vf = JC.pumpednoisematrices(bare, rows, Lb.cols, Lb.K)
        i = argmin(abs.(rows .- 2pi*0.4e9))
        @test isapprox(bn.covariance[1, 1] + bn.covariance[2, 2], sum(abs2, Sf[i, :]) + real(Vf[i, i]); rtol = 1e-4)
        hbb = hbsolve([2pi*0.4e9], (wp,), [], (nm,), (4,), one(bare); threewavemixing = true, returnCnoise = true)
        @test isapprox(real(Vf[i, i]), real(hbb.linearized.Cnoise((0,), 1, (0,), 1, 1)); rtol = 1e-10)
        push!(v2, real(Vf[i, i]))
    end
    @test isapprox(v2[1], v2[2]; rtol = 1e-4)
    # the noise of an idler survives export, fitting and the transient:
    # a solve holds its idler at a negative frequency, and the noise
    # there is the noise at the positive frequency a transient measures,
    # transposed, so the fitted block emits what the source does rather
    # than the vacuum its commutator alone would require
    src = RationalScattering(fill(-wp, 1, 1), fill(wp, 1, 1), fill(0.5, 1, 1), zeros(1, 1);
        noise = NoiseCovariance(fill(10.0, 1, 1)))
    sideband = hbsolve(wp .* collect(range(0.35, 0.45; length = 5)), (wp,), [], (1,), (2,), one(src);
        threewavemixing = true, returnCnoise = true).linearized
    fit = RationalScattering(LinearizedScattering(sideband, wp; noise = NoiseCovariance(sideband.Cnoise)), 1; noisetol = 1e-6)
    var = Float64[]
    for blk in (src, fit)
        isol = transientsolve(transientproblem(one(blk)), (0.0, 20e-9 - dt); dt, method = GaussLegendre(), record = :checkpoints)
        iplan = transientquantumplan(isol.times[501:end], [0.6e9])
        inoise = transientnoise(isol, iplan; frequencies = [0.6e9], weights = [1/10e-9])
        @test inoise.diagnostics.passed
        push!(var, tr(inoise.covariance))
    end
    @test isapprox(var[1], var[2]; rtol = 1e-10)
    @test var[1] > 10
    # the stage values of a pumped block belong to the factor of the
    # batch which steps it, so batches on one system do not share them
    p = transientproblem(one(model()))
    sys = JC.transientsystem(p, dt, GaussLegendre(), JC.CPU(), JC.KLUfactorization())
    bf1, bf2 = JC.gaussbatchfactor(sys, 1), JC.gaussbatchfactor(sys, 1)
    r1, r2 = JC.rationalwork(p, JC.CPU(), length(p), 1), JC.rationalwork(p, JC.CPU(), length(p), 1)
    JC.stageweights!(r1, sys, dt, 2dt)
    JC.refreshstageoperator!(r1, sys, bf1)
    v1 = copy(bf1.rationalvals)
    JC.stageweights!(r2, sys, 20dt, 21dt)
    JC.refreshstageoperator!(r2, sys, bf2)
    @test bf1.rationalvals == v1 && bf1.rationalvals != bf2.rationalvals
    @test !(bf1.rationalvals === sys.gauss.rationalvals) && !(bf1.hostvals === sys.gauss.hostvals)
    # the stage correction is rebuilt with the factorization when the
    # junctions of the circuit ask for a fresh one
    cj = one(model(amp = 0.5, noise = NoiseCovariance([fill(10.0, 1, 1), zero1])), (:jj, 1, 0, JosephsonJunction(1e-9)), (:cap, 1, 0, Capacitor(1e-12)))
    pj = transientproblem(cj)
    sysj = JC.transientsystem(pj, dt, GaussLegendre(), JC.CPU(), JC.KLUfactorization())
    bfj = JC.gaussbatchfactor(sysj, 1)
    st = JC.gaussstepper(sysj, [pj], 1e-12, 1e-13, 100, bfj)
    JC.setstate!(st, zeros(length(pj), 1), zeros(length(pj), 1), nothing)
    JC.stageweights!(st.rw, sysj, 0.18e-9, 0.22e-9)
    JC.refreshstageoperator!(st.rw, sysj, st.bf)
    st.refresh!()
    Kold = copy(st.rw.correction.K)
    fill!(st.phi, 1.2)
    st.refresh!()
    @test norm(st.rw.correction.K .- Kold) > 0
    r = randn(size(st.delta))
    a, b = zero(r), zero(r)
    st.solve!(a, r)
    JC.stagecorrection!(st.rw, sysj, st.bf, sysj.gauss.coefficients, st.rc, st.zc, false)
    st.solve!(b, r)
    @test a == b
    # and that circuit in time, a weak signal through the block beside
    # the junction, against the linearized solve of the same circuit
    fs, Is = 0.4e9, 1e-9
    fi = 1e9 - fs
    a0 = Is*sqrt(Z0)/2
    hb = hbsolve([2pi*fs], (wp,), [], (2,), (4,), cj; threewavemixing = true)
    tsol = transientsolve(transientproblem(cj; sources = [TransientSource(1, t -> Is*sinpi(2fs*t))]), (0.0, 400e-9);
        dt, method = GaussLegendre())
    iqplan = transientiqplan(tsol.times, [fs, fi]; duration = 40e-9, ports = [1, 1], stride = 400)
    iq = transientiq(iqplan, tsol.outgoing)
    @test isapprox(abs(iq[1, end])/a0, abs(hb.linearized.S((0,), 1, (0,), 1, 1)); rtol = 1e-3)
    @test isapprox(abs(iq[2, end])/a0, abs(hb.linearized.S((-1,), 1, (0,), 1, 1))*sqrt(fi/fs); rtol = 1e-2)
end
