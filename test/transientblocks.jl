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
