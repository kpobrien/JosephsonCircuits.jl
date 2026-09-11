using JosephsonCircuits
using LinearAlgebra
using Test

# The quantum noise of a transient against harmonic balance: a passive two
# port, whose linearized scattering and vacuum covariance the temporal mode
# noise must reproduce, warm loss against the linearized noise covariance,
# and a pumped Josephson amplifier against the harmonic balance gain and
# quantum efficiency. The forward and the adjoint method must agree with
# each other to roundoff, since they contract the same responses.
@testset "the quantum noise of a transient" begin
    JC = JosephsonCircuits

    @testset "a passive two port against the linearized solver" begin
        c = Circuit([:p1 => Port(1), :p2 => Port(2), :c1 => Capacitor(0.3e-12),
                :c2 => Capacitor(0.5e-12), :loss => Resistor(30.0)],
            [((:p1, 1), (:c1, 1), (:loss, 1)), ((:p2, 1), (:c2, 1), (:loss, 2)),
                ((:p1, 2), (:p2, 2), (:c1, 2), (:c2, 2), Ground)])
        prob = transientproblem(c)
        n, T = 512, 1e-9
        sol = transientsolve(prob, (0.0, T*(n - 1)/n); dt = T/n, record = :phases)
        plan = transientquantumplan(sol.times, [3e9, 3e9]; ports = [1, 2])
        hb = hblinsolve(2pi*[3e9], c; keyedarrays = false, returnSnoise = true, returnCnoise = true)
        S = hb.S[:, :, 1]
        expected = zeros(4, 4)
        for j in 1:2, k in 1:2
            s = S[j, k]
            expected[(2j - 1):2j, (2k - 1):2k] .= [real(s) imag(s); -imag(s) real(s)]
        end
        forward = transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], inputs = plan, method = :forward)
        adjoint = transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], inputs = plan)
        for r in (forward, adjoint)
            @test r.diagnostics.passed
            @test r.covariance ≈ plan.vacuum rtol=1e-5 atol=1e-6
            @test r.commutator ≈ plan.commutator rtol=1e-5 atol=1e-6
            @test r.gain ≈ expected rtol=1e-4 atol=1e-4
            for j in 1:2, k in 1:2
                qe = transientquantumefficiency(r.gain[(2j - 1):2j, (2k - 1):2k], r.covariance[(2j - 1):2j, (2j - 1):2j])
                @test qe.QE ≈ hb.QE[j, k, 1] rtol=3e-4
            end
        end
        @test adjoint.covariance ≈ forward.covariance rtol=1e-10
        @test adjoint.gain ≈ forward.gain rtol=1e-10
        # the pulsed gain of a probe applied inside the window against the
        # periodic gain: a rectangular probe of three cycles in the 1 ns
        # window differs by its edges, and the difference falls as the
        # inverse of the window
        pulsed = transientgain(sol, plan, plan)
        @test pulsed ≈ adjoint.gain rtol=6e-2
        long = transientsolve(prob, (0.0, 4T*(4n - 1)/(4n)); dt = T/n, record = :phases)
        longplan = transientquantumplan(long.times, [3e9, 3e9]; ports = [1, 2])
        @test transientgain(long, longplan, longplan) ≈ expected rtol=2e-2
        @test norm(transientgain(long, longplan, longplan) .- expected) < norm(pulsed .- expected)
        # the probe has the support of its window: a window measured
        # before the input window starts sees nothing of it
        before = transientquantumplan(sol.times[1:64], [1/(64*sol.dt)])
        after = transientquantumplan(sol.times[65:128], [1/(64*sol.dt)])
        @test all(iszero, transientgain(sol, before, after))
        @test !all(iszero, transientgain(sol, after, before))
        # warm loss: the occupation weights the covariance and not the
        # commutator, and the excess is the linearized noise covariance
        hot = Circuit([:p1 => Port(1), :p2 => Port(2), :c1 => Capacitor(0.3e-12),
                :c2 => Capacitor(0.5e-12), :loss => Resistor(30.0; temperature = 0.3)],
            [((:p1, 1), (:c1, 1), (:loss, 1)), ((:p2, 1), (:c2, 1), (:loss, 2)),
                ((:p1, 2), (:p2, 2), (:c1, 2), (:c2, 2), Ground)])
        baths = transientnoisebaths(transientproblem(hot))
        @test [b.temperature for b in baths.channels] == [0.0, 0.0, 0.3]
        warm = transientnoise(sol, plan; frequencies = [3e9], weights = [1/T],
            baths = JC.TransientNoiseBaths(prob, baths.channels))
        hbhot = hblinsolve(2pi*[3e9], hot; keyedarrays = false, returnSnoise = true, returnCnoise = true)
        @test warm.commutator ≈ adjoint.commutator rtol=1e-10
        @test warm.diagnostics.passed
        excess = (hbhot.Cnoise[:, :, 1] - hb.Cnoise[:, :, 1])/2
        expectedexcess = zeros(4, 4)
        for j in 1:2, k in 1:2
            z = excess[j, k]
            expectedexcess[(2j - 1):2j, (2k - 1):2k] .= [real(z) imag(z); -imag(z) real(z)]
        end
        @test warm.covariance - adjoint.covariance ≈ expectedexcess rtol=1e-4 atol=1e-5
        # the contracts: a state whose drift balances a constant drive with
        # a resistor while the junction phase moves is not an equilibrium
        moving = transientproblem([("P1", "1", "0", 1), ("R1", "1", "0", 50.0), ("C1", "1", "0", 1e-12),
            ("Lj1", "1", "0", 1e-9)]; sources = [TransientSource(1, 1e-6)])
        ms = transientsolve(moving, (0.0, T*(n - 1)/n); dt = T/n, record = :phases,
            initialstate = transientstate(moving; voltage = [50e-6]))
        mplan = transientquantumplan(ms.times, [3e9])
        @test_throws ArgumentError transientnoise(ms, mplan; frequencies = [3e9], weights = [1/T])
        @test_throws ArgumentError transientnoise(sol, plan; frequencies = [2e9], weights = [1/T], inputs = plan)
        @test_throws ArgumentError transientnoise(sol, plan; frequencies = [3e9], weights = [2/T], inputs = plan)
        @test_throws ArgumentError transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], method = :other)
        @test_throws ArgumentError transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], baths = transientnoisebaths(transientproblem(hot)))
    end

    @testset "scattering blocks as baths" begin
        # a lossy block emits the noise its loss requires, Bosma's
        # I - S S' as the linearized solver has it: a cold attenuator
        # leaves the vacuum the vacuum, a warm one equals its resistive
        # network at the same temperature and has the excess the
        # linearized solver gives it, a cascade at two temperatures
        # likewise, and a lossless block emits nothing
        function pinetwork(S, R)
            Rh = Diagonal(sqrt.(R))
            Y = inv(Rh)*(I - S)*inv(I + S)*inv(Rh)
            return (1/(Y[1, 1] + Y[1, 2]), -1/Y[1, 2], 1/(Y[2, 2] + Y[2, 1]))
        end
        block(x) = [real(x) imag(x); -imag(x) real(x)]
        quadratures(M) = reduce(vcat, [reduce(hcat, [block(M[j, k]) for k in 1:2]) for j in 1:2])
        g = 0.6
        S = [0.0 g; g 0.0]
        r1, rs, r2 = pinetwork(S, [50.0, 50.0])
        n, T = 512, 1e-9
        mk(att) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)), (:att, 1, 2, att),
            (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        res(temp) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)),
            (:r1, 1, 0, Resistor(r1; temperature = temp)), (:rs, 1, 2, Resistor(rs; temperature = temp)),
            (:r2, 2, 0, Resistor(r2; temperature = temp)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        function noiseof(c; method = :adjoint)
            prob = transientproblem(c)
            sol = transientsolve(prob, (0.0, T*(n - 1)/n); dt = T/n, record = :phases, method = GaussLegendre())
            plan = transientquantumplan(sol.times, [3e9, 3e9]; ports = [1, 2])
            return transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], inputs = plan, method), plan
        end
        cold, plan = noiseof(mk(ScatteringParameters(S; zref = 50.0)))
        baths = transientnoisebaths(transientproblem(mk(ScatteringParameters(S; zref = 50.0))))
        @test length(baths) == 4 && all(b -> iszero(b.resistance) && length(b.rows) == 2, baths.channels[3:4])
        @test cold.diagnostics.passed
        @test cold.covariance ≈ plan.vacuum rtol=1e-6
        @test cold.commutator ≈ plan.commutator rtol=1e-6
        hb = hblinsolve(2pi*[3e9], mk(ScatteringParameters(S; zref = 50.0)); keyedarrays = false, returnCnoise = true)
        @test cold.gain ≈ quadratures(hb.S[:, :, 1]) rtol=1e-6
        forward, _ = noiseof(mk(ScatteringParameters(S; zref = 50.0)); method = :forward)
        @test forward.covariance ≈ cold.covariance rtol=1e-10
        warm, _ = noiseof(mk(ScatteringParameters(S; zref = 50.0, noise = ThermalEquilibrium(0.3))))
        network, _ = noiseof(res(0.3))
        @test warm.covariance ≈ network.covariance rtol=1e-10
        @test warm.commutator ≈ cold.commutator rtol=1e-10
        hbw = hblinsolve(2pi*[3e9], mk(ScatteringParameters(S; zref = 50.0, noise = ThermalEquilibrium(0.3))); keyedarrays = false, returnCnoise = true)
        @test warm.covariance - cold.covariance ≈ quadratures((hbw.Cnoise[:, :, 1] - hb.Cnoise[:, :, 1])/2) rtol=1e-6
        cascade(t1, t2) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)),
            (:a1, 1, 2, ScatteringParameters(S; zref = 50.0, noise = ThermalEquilibrium(t1))), (:c2, 2, 0, Capacitor(0.2e-12)),
            (:a2, 2, 3, ScatteringParameters(S; zref = 50.0, noise = ThermalEquilibrium(t2))), (:c3, 3, 0, Capacitor(0.5e-12)),
            (:p2, 3, 0, Port(2; Z0 = 50.0))])
        cw, _ = noiseof(cascade(0.3, 4.0))
        cc, _ = noiseof(cascade(0.0, 0.0))
        hbc = hblinsolve(2pi*[3e9], cascade(0.3, 4.0); keyedarrays = false, returnCnoise = true)
        hbc0 = hblinsolve(2pi*[3e9], cascade(0.0, 0.0); keyedarrays = false, returnCnoise = true)
        @test cc.covariance ≈ plan.vacuum rtol=1e-6
        @test cw.covariance - cc.covariance ≈ quadratures((hbc.Cnoise[:, :, 1] - hbc0.Cnoise[:, :, 1])/2) rtol=1e-6
        Sc = [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
        circ = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)),
            (:circ, 1, 2, 3, ScatteringParameters(Sc; zref = 50.0)), (:c2, 2, 0, Capacitor(1e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0)),
            (:c3, 3, 0, Capacitor(0.2e-12)), (:p3, 3, 0, Port(3; Z0 = 50.0))])
        @test length(transientnoisebaths(transientproblem(circ))) == 3
        @test_throws ArgumentError transientnoisebaths(transientproblem(mk(ScatteringParameters(S; zref = 50.0, noise = Lossless()))))
        # the pumped amplifier behind a cold 3 dB pad: its gain and its
        # quantum efficiency against the linearized solver
        ga = 10^(-3/20)
        jpa = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:pad, 1, 2, ScatteringParameters([0.0 ga; ga 0.0]; zref = 50.0)),
            (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)), (:cj, 3, 0, Capacitor(1000e-15))])
        fp, fs, ip = 4.75e9, 4.7e9, 0.00565e-6
        jhb = hbsolve([2pi*fs], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (16,), jpa; ftol = 1e-14)
        jgain, jqe = abs2(jhb.linearized.S((0,), 1, (0,), 1, 1)), jhb.linearized.QE((0,), 1, (0,), 1, 1)
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        jprob = transientproblem(jpa; sources = [TransientSource(1, t -> 2ip*ramp(t)*cospi(2fp*t))])
        settle, record, dt = 100e-9, 20e-9, 2.5e-12
        jsol = transientsolve(jprob, (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
        first = round(Int, settle/dt) + 1
        jplan = transientquantumplan(jsol.times[first:end], [fs])
        jfreqs = sort!(abs.([fs + 2k*fp for k in -2:2]))
        jnoise = transientnoise(jsol, jplan; frequencies = jfreqs, weights = fill(1/record, 5), inputs = jplan, commutationrtol = 3e-3)
        @test jnoise.diagnostics.passed
        jmetrics = transientquantumefficiency(jnoise.gain, jnoise.covariance; rtol = 3e-3)
        @test jmetrics.gain ≈ jgain rtol=1e-3
        @test jmetrics.QE ≈ jqe rtol=1e-3
    end

    @testset "transmission lines carry the baths' prehistory" begin
        # a line stores the fluctuations that entered it before the start:
        # the stationary response with the lines' wave phasors gives the
        # tangent its prehistory and the adjoint its initial term, so a
        # cold mismatched line between capacitive loads leaves the vacuum
        # the vacuum and has the linearized solver's gain, and a pumped
        # amplifier behind a 60 ohm cable has its gain and quantum
        # efficiency
        block(x) = [real(x) imag(x); -imag(x) real(x)]
        quadratures(M) = reduce(vcat, [reduce(hcat, [block(M[j, k]) for k in 1:2]) for j in 1:2])
        n, T, tau = 512, 1e-9, 0.3e-9
        c = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)),
            (:line, 1, 2, TransmissionLine(60.0, tau*3e8; vp = 3e8)), (:c2, 2, 0, Capacitor(0.4e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        sol = transientsolve(transientproblem(c), (0.0, T*(n - 1)/n); dt = T/n, record = :phases, method = GaussLegendre())
        plan = transientquantumplan(sol.times, [3e9, 3e9]; ports = [1, 2])
        adjoint = transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], inputs = plan)
        forward = transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], inputs = plan, method = :forward)
        hb = hblinsolve(2pi*[3e9], c; keyedarrays = false)
        @test adjoint.diagnostics.passed
        @test adjoint.covariance ≈ plan.vacuum rtol=1e-6
        @test adjoint.commutator ≈ plan.commutator rtol=1e-6
        @test adjoint.gain ≈ quadratures(hb.S[:, :, 1]) rtol=1e-6
        @test forward.covariance ≈ adjoint.covariance rtol=1e-10
        @test forward.gain ≈ adjoint.gain rtol=1e-10
        jpa = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:cable, 1, 2, TransmissionLine(60.0, tau*3e8; vp = 3e8)),
            (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)), (:cj, 3, 0, Capacitor(1000e-15))])
        fp, fs, ip = 4.75e9, 4.7e9, 0.00565e-6
        jhb = hbsolve([2pi*fs], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (16,), jpa; ftol = 1e-14)
        jgain, jqe = abs2(jhb.linearized.S((0,), 1, (0,), 1, 1)), jhb.linearized.QE((0,), 1, (0,), 1, 1)
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        jprob = transientproblem(jpa; sources = [TransientSource(1, t -> 2ip*ramp(t)*cospi(2fp*t))])
        settle, record, dt = 100e-9, 20e-9, 2.5e-12
        jsol = transientsolve(jprob, (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
        first = round(Int, settle/dt) + 1
        jplan = transientquantumplan(jsol.times[first:end], [fs])
        jfreqs = sort!(abs.([fs + 2k*fp for k in -2:2]))
        jnoise = transientnoise(jsol, jplan; frequencies = jfreqs, weights = fill(1/record, 5), inputs = jplan, commutationrtol = 3e-3)
        @test jnoise.diagnostics.passed
        jmetrics = transientquantumefficiency(jnoise.gain, jnoise.covariance; rtol = 3e-3)
        @test jmetrics.gain ≈ jgain rtol=1e-5
        @test jmetrics.QE ≈ jqe rtol=1e-5
    end

    @testset "rational blocks emit the noise of their loss at every frequency" begin
        # a lossless rational block, the series inductor, has the noise
        # of the explicit inductor and adds no channel of its own; a lossy
        # one carries its ports as correlated channels with Bosma's
        # covariance at each frequency, leaves the vacuum the vacuum when
        # cold, has the linearized solver's gain and, warm, its excess;
        # forward and adjoint agree, and the amplifier behind it has the
        # gain and quantum efficiency of harmonic balance
        block(x) = [real(x) imag(x); -imag(x) real(x)]
        quadratures(M) = reduce(vcat, [reduce(hcat, [block(M[j, k]) for k in 1:2]) for j in 1:2])
        n, T = 512, 1e-9
        R0, L = 50.0, 2e-9
        a = 2R0/L
        u = [1.0, -1.0]
        ind = RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
        al = 2pi*2e9
        lossy(noise) = RationalScattering(-al .* Matrix(1.0I, 2, 2), al .* Matrix(1.0I, 2, 2), 0.8 .* [0.0 1.0; 1.0 0.0], zeros(2, 2); zref = 50.0, noise)
        mk(b) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)), (:b, 1, 2, b), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        explicit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)), (:l, 1, 2, Inductor(L)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        function noiseof(c; method = :adjoint)
            sol = transientsolve(transientproblem(c), (0.0, T*(n - 1)/n); dt = T/n, record = :phases, method = GaussLegendre())
            plan = transientquantumplan(sol.times, [3e9, 3e9]; ports = [1, 2])
            return transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], inputs = plan, method), plan
        end
        ni, plan = noiseof(mk(ind))
        ne, _ = noiseof(explicit)
        @test ni.covariance ≈ ne.covariance rtol=1e-10
        @test ni.gain ≈ ne.gain rtol=1e-10
        @test length(transientnoisebaths(transientproblem(mk(ind)))) == 4
        @test_throws ArgumentError transientnoisebaths(transientproblem(mk(lossy(Lossless()))))
        @test length(transientnoisebaths(transientproblem(mk(RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0, noise = Lossless()))))) == 2
        nl, _ = noiseof(mk(lossy(Passive())))
        hb = hblinsolve(2pi*[3e9], mk(lossy(Passive())); keyedarrays = false, returnCnoise = true)
        @test nl.diagnostics.passed
        @test nl.covariance ≈ plan.vacuum rtol=1e-6
        @test nl.commutator ≈ plan.commutator rtol=1e-6
        @test nl.gain ≈ quadratures(hb.S[:, :, 1]) rtol=1e-6
        nf, _ = noiseof(mk(lossy(Passive())); method = :forward)
        @test nf.covariance ≈ nl.covariance rtol=1e-10
        @test nf.gain ≈ nl.gain rtol=1e-10
        nw, _ = noiseof(mk(lossy(ThermalEquilibrium(0.3))))
        hbw = hblinsolve(2pi*[3e9], mk(lossy(ThermalEquilibrium(0.3))); keyedarrays = false, returnCnoise = true)
        @test nw.covariance - nl.covariance ≈ quadratures((hbw.Cnoise[:, :, 1] - hb.Cnoise[:, :, 1])/2) rtol=1e-6
        # the loss matrix `I - S S'` of that block is a scalar, which hides
        # the sign of the sine quadrature in the group's correction; a
        # feedthrough of opposite signs at the two ports gives a loss
        # matrix with imaginary off-diagonal entries, which does not
        mixed(noise) = RationalScattering(-al .* Matrix(1.0I, 2, 2), al .* Matrix(1.0I, 2, 2), 0.6 .* [0.0 1.0; 1.0 0.0], [0.3 0.0; 0.0 -0.3]; zref = 50.0, noise)
        K = JosephsonCircuits.groupcovariance(transientnoisebaths(transientproblem(mk(mixed(Passive())))), (channels = 3:4, block = 1), 3e9)
        @test abs(imag(K[1, 2])) > 0.1
        nm, _ = noiseof(mk(mixed(Passive())))
        @test nm.diagnostics.passed
        @test nm.covariance ≈ plan.vacuum rtol=1e-6
        nmw, _ = noiseof(mk(mixed(ThermalEquilibrium(0.3))))
        hbm = hblinsolve(2pi*[3e9], mk(mixed(Passive())); keyedarrays = false, returnCnoise = true)
        hbmw = hblinsolve(2pi*[3e9], mk(mixed(ThermalEquilibrium(0.3))); keyedarrays = false, returnCnoise = true)
        @test nmw.covariance - nm.covariance ≈ quadratures((hbmw.Cnoise[:, :, 1] - hbm.Cnoise[:, :, 1])/2) rtol=1e-6
        # Nothing is inferred about a rational block's loss: a notch of
        # relative width 1e-6 at 3 GHz, an all pass elsewhere to a part in
        # 1e12, keeps its channels and, warm, emits the noise of its loss
        # at the notch as the linearized solver has it; its groups are
        # contracted directly with their covariance, so the small loss is
        # not the difference of two large terms; and a declaration of
        # Lossless() on it is refused, while the inductor block's holds
        wn, eps = 2pi*3e9, 1e-6
        Ann = wn .* [0.0 1.0; -1.0 -2eps]
        Cnn = [0.0 -2eps]
        notch(noise) = RationalScattering([Ann zeros(2, 1); wn .* Cnn fill(-10wn, 1, 1)], wn .* [0.0; 1.0; 1.0;;], [Cnn fill(-20.0, 1, 1)], ones(1, 1);
            zref = 50.0, noise)
        one(b) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:b, 1, b)])
        function noiseone(c)
            sol = transientsolve(transientproblem(c), (0.0, T*(n - 1)/n); dt = T/n, record = :phases, method = GaussLegendre())
            plan1 = transientquantumplan(sol.times, [3e9]; ports = [1])
            return transientnoise(sol, plan1; frequencies = [3e9], weights = [1/T], inputs = plan1), plan1
        end
        @test length(transientnoisebaths(transientproblem(one(notch(Passive()))))) == 2
        ncold, plan1 = noiseone(one(notch(Passive())))
        nwarm, _ = noiseone(one(notch(ThermalEquilibrium(0.3))))
        @test ncold.covariance ≈ plan1.vacuum rtol=1e-6
        hcold = hblinsolve(2pi*[3e9], one(notch(Passive())); keyedarrays = false, returnCnoise = true)
        hwarm = hblinsolve(2pi*[3e9], one(notch(ThermalEquilibrium(0.3))); keyedarrays = false, returnCnoise = true)
        Sn = zeros(ComplexF64, 1, 1, 1)
        JosephsonCircuits.evaluateprovider!(Sn, notch(Passive()).provider, [wn])
        @test abs(Sn[1, 1, 1]) < 1e-3
        excess = real(hwarm.Cnoise[1, 1, 1] - hcold.Cnoise[1, 1, 1])/2
        @test excess > 0.1
        @test nwarm.covariance - ncold.covariance ≈ excess .* Matrix(1.0I, 2, 2) rtol=1e-5
        @test_throws ArgumentError notch(Lossless())
        @test RationalScattering(fill(-a, 1, 1), reshape(u, 1, 2), reshape(-a .* u, 2, 1), Matrix(1.0I, 2, 2); zref = 50.0, noise = Lossless()) isa ScatteringParameters
        jpa = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:b, 1, 2, lossy(Passive())), (:cc, 2, 3, Capacitor(100e-15)),
            (:jj, 3, 0, JosephsonJunction(1000e-12)), (:cj, 3, 0, Capacitor(1000e-15))])
        fp, fs, ip = 4.75e9, 4.7e9, 0.00565e-6
        jhb = hbsolve([2pi*fs], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (16,), jpa; ftol = 1e-14)
        jgain, jqe = abs2(jhb.linearized.S((0,), 1, (0,), 1, 1)), jhb.linearized.QE((0,), 1, (0,), 1, 1)
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        jprob = transientproblem(jpa; sources = [TransientSource(1, t -> 2ip*ramp(t)*cospi(2fp*t))])
        settle, record, dt = 100e-9, 20e-9, 2.5e-12
        jsol = transientsolve(jprob, (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
        first = round(Int, settle/dt) + 1
        jplan = transientquantumplan(jsol.times[first:end], [fs])
        jfreqs = sort!(abs.([fs + 2k*fp for k in -2:2]))
        jnoise = transientnoise(jsol, jplan; frequencies = jfreqs, weights = fill(1/record, 5), inputs = jplan, commutationrtol = 3e-3)
        @test jnoise.diagnostics.passed
        jmetrics = transientquantumefficiency(jnoise.gain, jnoise.covariance; rtol = 3e-3)
        @test jmetrics.gain ≈ jgain rtol=1e-5
        @test jmetrics.QE ≈ jqe rtol=1e-5
    end

    @testset "an equilibrium with a direct current on a line or in a block" begin
        # the noise needs a classical equilibrium at the start, and the
        # check reads the whole of it: a matched line carrying a direct
        # current, and a rational block at rest under a bias, are
        # equilibria (the check once left out the lines' forcing and the
        # blocks' resting waves and refused them); a line whose waves are
        # not the bias's, or a block whose states are not at rest, is not
        n, T = 512, 1e-9
        dt = T/n
        line = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line, 1, 2, TransmissionLine(50.0, 10.5e-12; vp = 1.0)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        dc = transientproblem(line; sources = [TransientSource(1, 1e-6)])
        biased = transientstate(dc; voltage = [25e-6, 25e-6], linecurrents = [0.5e-6])
        ds = transientsolve(dc, (0.0, T*(n - 1)/n); dt, record = :phases, method = GaussLegendre(), initialstate = biased, rtol = 1e-12, atol = 1e-13)
        @test maximum(abs.(ds.voltage .- 25e-6)) < 1e-15
        plan = transientquantumplan(ds.times, [3e9, 3e9]; ports = [1, 2])
        nd = transientnoise(ds, plan; frequencies = [3e9], weights = [1/T], inputs = plan)
        @test nd.covariance ≈ plan.vacuum rtol=1e-6
        @test_throws ArgumentError transientnoise(transientsolve(dc, (0.0, T*(n - 1)/n); dt, record = :phases, method = GaussLegendre(),
            initialstate = transientstate(dc; voltage = [25e-6, 25e-6], linecurrents = [0.0])), plan; frequencies = [3e9], weights = [1/T])
        al = 2pi*2e9
        lossy = RationalScattering(-al .* Matrix(1.0I, 2, 2), al .* Matrix(1.0I, 2, 2), 0.6 .* [0.0 1.0; 1.0 0.0], [0.3 0.0; 0.0 -0.3]; zref = 50.0)
        blockc = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:b, 1, 2, lossy), (:c2, 2, 0, Capacitor(0.3e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        bdc = transientproblem(blockc; sources = [TransientSource(1, 1e-6)])
        # the bias's node voltages from a long settle, then the rest state
        settled = transientsolve(bdc, (0.0, 20e-9); dt = 2e-12, method = GaussLegendre())
        vb = settled.voltage[:, end]
        rest = transientstate(bdc; voltage = vb)
        bs = transientsolve(bdc, (0.0, T*(n - 1)/n); dt, record = :phases, method = GaussLegendre(), initialstate = rest, rtol = 1e-12, atol = 1e-13)
        @test maximum(abs.(bs.voltage .- vb)) < 1e-6*maximum(abs, vb)
        nb = transientnoise(bs, plan; frequencies = [3e9], weights = [1/T], inputs = plan)
        @test nb.covariance ≈ plan.vacuum rtol=1e-6
        moved = (rest.x, rest.v, rest.waves, rest.states .* 0.5)
        @test_throws ArgumentError transientnoise(transientsolve(bdc, (0.0, T*(n - 1)/n); dt, record = :phases, method = GaussLegendre(),
            initialstate = moved), plan; frequencies = [3e9], weights = [1/T])
    end

    @testset "a fitted block has the noise of the circuit it was fitted to" begin
        # the RLC two-port's data carries the temperature of its resistor
        # as the block's noise model: fitted, the block leaves the vacuum
        # the vacuum when cold and has the explicit warm resistor's excess
        # in time as in the linearized solver
        block(x) = [real(x) imag(x); -imag(x) real(x)]
        quadratures(M) = reduce(vcat, [reduce(hcat, [block(M[j, k]) for k in 1:2]) for j in 1:2])
        n, T = 512, 1e-9
        rlc(R) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l1, 1, 2, Inductor(1.5e-9)), (:c, 2, 0, Capacitor(0.6e-12)),
            (:r, 2, 0, R), (:l2, 2, 3, Inductor(1.0e-9)), (:p2, 3, 0, Port(2; Z0 = 50.0))])
        fs = collect(range(0.2e9, 12e9; length = 240))
        hb = hblinsolve(2pi .* fs, rlc(Resistor(120.0)); keyedarrays = false)
        fitted(noise) = RationalScattering(ScatteringParameters((2pi .* fs, hb.S); nports = 2, zref = 50.0, noise), 4)
        mk(b) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:b, 1, 2, b), (:c2, 2, 0, Capacitor(0.3e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        explicit(T) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:l1, 1, 4, Inductor(1.5e-9)),
            (:c, 4, 0, Capacitor(0.6e-12)), (:r, 4, 0, Resistor(120.0; temperature = T)), (:l2, 4, 2, Inductor(1.0e-9)),
            (:c2, 2, 0, Capacitor(0.3e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        function noiseof(c)
            sol = transientsolve(transientproblem(c), (0.0, T*(n - 1)/n); dt = T/n, record = :phases, method = GaussLegendre())
            plan = transientquantumplan(sol.times, [3e9, 3e9]; ports = [1, 2])
            return transientnoise(sol, plan; frequencies = [3e9], weights = [1/T], inputs = plan), plan
        end
        nc, plan = noiseof(mk(fitted(Passive())))
        @test nc.diagnostics.passed
        @test nc.covariance ≈ plan.vacuum rtol=1e-6
        @test nc.gain ≈ quadratures(hblinsolve(2pi*[3e9], explicit(0.0); keyedarrays = false).S[:, :, 1]) rtol=1e-6
        nw, _ = noiseof(mk(fitted(ThermalEquilibrium(0.3))))
        ec, _ = noiseof(explicit(0.0))
        ew, _ = noiseof(explicit(0.3))
        @test nw.covariance - nc.covariance ≈ ew.covariance - ec.covariance rtol=1e-6
        hbc = hblinsolve(2pi*[3e9], explicit(0.0); keyedarrays = false, returnCnoise = true)
        hbw = hblinsolve(2pi*[3e9], explicit(0.3); keyedarrays = false, returnCnoise = true)
        @test nw.covariance - nc.covariance ≈ quadratures((hbw.Cnoise[:, :, 1] - hbc.Cnoise[:, :, 1])/2) rtol=1e-6
    end

    @testset "a composed front end: lines, a warm fitted block and a pumped junction" begin
        # the composition the blocks and lines are for: a mismatched line,
        # a lossy block fitted from data at its own temperature, a second
        # mismatched line, and a pumped junction, driven by the pump's
        # rise; the stationary gain and quantum efficiency, with the warm
        # block's excess, are harmonic balance's
        fs = collect(range(0.2e9, 12e9; length = 240))
        rlc = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l1, 1, 2, Inductor(1.5e-9)), (:c, 2, 0, Capacitor(0.6e-12)),
            (:r, 2, 0, Resistor(120.0)), (:l2, 2, 3, Inductor(1.0e-9)), (:p2, 3, 0, Port(2; Z0 = 50.0))])
        data = ScatteringParameters((2pi .* fs, hblinsolve(2pi .* fs, rlc; keyedarrays = false).S); nports = 2, zref = 50.0,
            noise = ThermalEquilibrium(0.3))
        warm = RationalScattering(data, 4)
        front = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line1, 1, 2, TransmissionLine(60.0, 0.09)), (:b, 2, 3, warm),
            (:line2, 3, 4, TransmissionLine(40.0, 0.05)), (:cc, 4, 5, Capacitor(100e-15)),
            (:jj, 5, 0, JosephsonJunction(1000e-12)), (:cj, 5, 0, Capacitor(1000e-15))])
        fp, fsig, ip = 4.75e9, 4.7e9, 0.00565e-6
        hbr = hbsolve([2pi*fsig], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (16,), front; ftol = 1e-14)
        hgain, hqe = abs2(hbr.linearized.S((0,), 1, (0,), 1, 1)), hbr.linearized.QE((0,), 1, (0,), 1, 1)
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        prob = transientproblem(front; sources = [TransientSource(1, t -> 2ip*ramp(t)*cospi(2fp*t))])
        settle, record, dt = 100e-9, 20e-9, 2.5e-12
        sol = transientsolve(prob, (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
        first = round(Int, settle/dt) + 1
        plan = transientquantumplan(sol.times[first:end], [fsig])
        freqs = sort!(abs.([fsig + 2k*fp for k in -2:2]))
        noise = transientnoise(sol, plan; frequencies = freqs, weights = fill(1/record, 5), inputs = plan, commutationrtol = 3e-3)
        @test noise.diagnostics.passed
        metrics = transientquantumefficiency(noise.gain, noise.covariance; rtol = 3e-3)
        @test metrics.gain ≈ hgain rtol=1e-4
        @test metrics.QE ≈ hqe rtol=1e-4
        # the warm block's excess: the same circuit with the block cold
        cold = RationalScattering(ScatteringParameters((2pi .* fs, hblinsolve(2pi .* fs, rlc; keyedarrays = false).S); nports = 2, zref = 50.0), 4)
        frontcold = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:line1, 1, 2, TransmissionLine(60.0, 0.09)), (:b, 2, 3, cold),
            (:line2, 3, 4, TransmissionLine(40.0, 0.05)), (:cc, 4, 5, Capacitor(100e-15)),
            (:jj, 5, 0, JosephsonJunction(1000e-12)), (:cj, 5, 0, Capacitor(1000e-15))])
        hbc = hbsolve([2pi*fsig], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (16,), frontcold; ftol = 1e-14)
        @test hqe < hbc.linearized.QE((0,), 1, (0,), 1, 1)
        solc = transientsolve(transientproblem(frontcold; sources = [TransientSource(1, t -> 2ip*ramp(t)*cospi(2fp*t))]),
            (0.0, settle + record - dt); dt, method = GaussLegendre(), record = :checkpoints)
        noisec = transientnoise(solc, plan; frequencies = freqs, weights = fill(1/record, 5), inputs = plan, commutationrtol = 3e-3)
        metricsc = transientquantumefficiency(noisec.gain, noisec.covariance; rtol = 3e-3)
        @test metricsc.QE ≈ hbc.linearized.QE((0,), 1, (0,), 1, 1) rtol=1e-4
        @test metricsc.QE > metrics.QE
    end
end
