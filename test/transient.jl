using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Random
using Test

@testset "the circuit in time" begin
    JC = JosephsonCircuits
    rc = [("P1", "1", "0", 1), ("R1", "1", "0", 50.0), ("C1", "1", "0", 1e-12)]

    @testset "the RC response, the rules, the ports and the sources" begin
        prob = transientproblem(rc; sources = [TransientSource(1, 1e-6)])
        exact = 50e-6*(1 - exp(-4))
        errors = map((5e-12, 2.5e-12, 1.25e-12)) do dt
            abs(transientsolve(prob, (0.0, 200e-12); dt).voltage[1, end] - exact)
        end
        # second order: the error falls by four per halving
        @test errors[1]/errors[2] ≈ 4 rtol=0.01
        @test errors[2]/errors[3] ≈ 4 rtol=0.01
        coarse = transientsolve(prob, (0.0, 200e-12); dt = 5e-12)
        @test coarse.stats.factorizations == 1
        @test coarse.incident ≈ fill(1e-6*sqrt(50)/2, size(coarse.incident))
        @test coarse.voltage ≈ sqrt(50)*(coarse.incident + coarse.outgoing)
        @test isnothing(coarse.flux)
        decimated = transientsolve(prob, (0.0, 200e-12); dt = 5e-12, saveevery = 7)
        @test decimated.times[[1, end]] == [0.0, 200e-12]
        @test decimated.finalflux == coarse.finalflux
        @test length(decimated.times) == 7
        # backward Euler has a closed form on the RC
        be = transientsolve(prob, (0.0, 200e-12); dt = 5e-12, method = BackwardEuler())
        @test be.voltage[1, end] ≈ 50e-6*(1 - (1 + 0.1)^(-40)) rtol=1e-10
        # a named current source flows out of its first terminal, a port
        # source into the port's positive terminal; a named drive replaces
        # the constant value, an unnamed constant source stays
        net = vcat(rc, [("I1", "1", "0", 2e-6)])
        named = transientsolve(transientproblem(net; sources = [TransientSource(:I1, 1e-6)]),
            (0.0, 200e-12); dt = 5e-12)
        @test named.voltage ≈ -coarse.voltage
        static = transientsolve(transientproblem(net), (0.0, 200e-12); dt = 5e-12)
        @test static.voltage ≈ -2coarse.voltage
        typed = Circuit([("p", "1", "0", Port(1)), ("c", "1", "0", Capacitor(1e-12))])
        ts = transientsolve(transientproblem(typed; sources = [TransientSource(1, 1e-6)]),
            (0.0, 200e-12); dt = 5e-12)
        @test ts.voltage ≈ coarse.voltage
    end

    @testset "a lossless LC keeps its energy, and coupled inductors their modes" begin
        L, C, V = 1e-9, 1e-12, 1e-6
        circuit = Circuit([("p", "1", "0", Port(1; termination = nothing)),
            ("c", "1", "0", Capacitor(C)), ("l", "1", "0", Inductor(L))])
        prob = transientproblem(circuit)
        period = 2pi*sqrt(L*C)
        sol = transientsolve(prob, (0.0, 3period); dt = period/500,
            initialstate = transientstate(prob; voltage = [V]), record = :states)
        v = sol.rate[1, :] .* JC.phi0
        flux = sol.flux[1, :] .* JC.phi0
        energy = C .* v .^ 2 ./ 2 .+ flux .^ 2 ./ (2L)
        @test maximum(abs.(energy ./ energy[1] .- 1)) < 2e-10
        @test maximum(abs.(v .- V .* cos.(sol.times ./ sqrt(L*C)))) < 3e-4V
        # mutually coupled inductors go through the augmentation of the
        # modified nodal analysis: two auxiliary currents, and the normal
        # modes of the coupled tanks
        coupled = [("C1", "1", "0", C), ("C2", "2", "0", C),
            ("L1", "1", "0", L), ("L2", "2", "0", L), ("K1", "L1", "L2", 0.3)]
        p2 = transientproblem(coupled)
        @test p2.Naux == 2 && length(p2) == 4 && isempty(p2.gaugeindices)
        s2 = transientsolve(p2, (0.0, 2period); dt = period/800,
            initialstate = transientstate(p2; voltage = [V, 0]), record = :states)
        wp, wm = 1/sqrt(L*C*1.3), 1/sqrt(L*C*0.7)
        e1 = V/2 .* (cos.(wp .* s2.times) .+ cos.(wm .* s2.times))
        e2 = V/2 .* (cos.(wp .* s2.times) .- cos.(wm .* s2.times))
        @test maximum(abs.(s2.rate[1, :] .* JC.phi0 .- e1)) < 2e-4V
        @test maximum(abs.(s2.rate[2, :] .* JC.phi0 .- e2)) < 2e-4V
        # a nonzero initial flux on a coupled tank: the auxiliary currents
        # of the state are in the problem's fixed units, so the state is
        # consistent at any step, and the fluxes follow the normal modes
        F = 0.1*JC.phi0
        s3 = transientsolve(p2, (0.0, 2period); dt = period/800,
            initialstate = transientstate(p2; flux = [F, 0]), record = :states)
        f1 = F/2 .* (cos.(wp .* s3.times) .+ cos.(wm .* s3.times))
        f2 = F/2 .* (cos.(wp .* s3.times) .- cos.(wm .* s3.times))
        @test maximum(abs.(s3.flux[1, :] .* JC.phi0 .- f1)) < 2e-4F
        @test maximum(abs.(s3.flux[2, :] .* JC.phi0 .- f2)) < 2e-4F
        # the auxiliary rates follow the node rates by the same constitutive
        # equations, and a restart at another step continues the solve
        state = transientstate(p2; flux = [F, 0], voltage = [V, 0])
        @test state.x[3:4] ≈ (p2.Lscale/JC.phi0) .* ([L 0.3L; 0.3L L] \ [F, 0])
        @test state.v[3:4] ≈ (p2.Lscale/JC.phi0) .* ([L 0.3L; 0.3L L] \ [V, 0])
        half = transientsolve(p2, (0.0, period); dt = period/800, initialstate = state)
        rest = transientsolve(p2, (period, 2period); dt = period/800,
            initialstate = (half.finalflux, half.finalrate))
        whole = transientsolve(p2, (0.0, 2period); dt = period/800, initialstate = state)
        @test rest.finalflux ≈ whole.finalflux rtol=1e-9
        @test rest.finalrate ≈ whole.finalrate rtol=1e-9
        coarse = transientsolve(p2, (period, 2period); dt = period/400,
            initialstate = (half.finalflux, half.finalrate))
        @test coarse.finalflux ≈ whole.finalflux rtol=1e-3
    end

    @testset "the initial state, the gauge and the unsupported cases" begin
        # a resistor driven from the first sample needs its voltage supplied
        resistor = transientproblem(rc[1:2]; sources = [TransientSource(1, 1e-6)])
        @test_throws ArgumentError transientsolve(resistor, (0.0, 1e-9); dt = 1e-12)
        rs = transientsolve(resistor, (0.0, 1e-9); dt = 1e-12,
            initialstate = transientstate(resistor; voltage = [50e-6]))
        @test rs.voltage ≈ fill(50e-6, size(rs.voltage))
        # a singular capacitance matrix without a zero row: two terminated
        # nodes joined by one capacitor constrain the sum of their voltages
        # to the drive, so the zero state is rejected, and the consistent
        # one relaxes without ringing to the resistive division
        pair = transientproblem([("P1", "1", "0", 1), ("R1", "1", "0", 50.0), ("P2", "2", "0", 2),
            ("R2", "2", "0", 50.0), ("C1", "1", "2", 1e-12)]; sources = [TransientSource(1, 1e-6)])
        @test pair.inertialess == [[1, 2]] && isempty(pair.algebraic)
        @test_throws ArgumentError transientsolve(pair, (0.0, 1e-9); dt = 1e-12)
        ps = transientsolve(pair, (0.0, 1e-9); dt = 1e-12,
            initialstate = transientstate(pair; voltage = [25e-6, 25e-6]))
        decay = exp.(-ps.times ./ 100e-12)
        @test ps.voltage[1, :] ≈ 50e-6 .- 25e-6 .* decay rtol=2e-4
        @test ps.voltage[2, :] ≈ 25e-6 .* decay rtol=2e-4 atol=1e-9
        @test maximum(abs.(ps.voltage[1, :] .+ ps.voltage[2, :] .- 50e-6)) < 1e-12
        # a capacitive island reached through an inductor, and a node no
        # capacitor touches between resistors: the same check
        island2 = transientproblem([("P1", "1", "0", 1), ("R1", "1", "0", 50.0), ("C1", "1", "2", 1e-12),
            ("L2", "2", "0", 1e-9)]; sources = [TransientSource(1, 1e-6)])
        @test island2.inertialess == [[1, 2]] && isempty(island2.algebraic)
        @test_throws ArgumentError transientsolve(island2, (0.0, 1e-9); dt = 1e-12)
        i2 = transientsolve(island2, (0.0, 1e-9); dt = 1e-12,
            initialstate = transientstate(island2; voltage = [50e-6, 50e-6]))
        @test i2.voltage[1, 1] ≈ 50e-6
        divider = transientproblem([("P1", "1", "0", 1), ("R1", "1", "0", 50.0), ("C1", "1", "0", 1e-12),
            ("R2", "1", "2", 50.0), ("R3", "2", "0", 50.0)]; sources = [TransientSource(1, 1e-6)])
        @test divider.inertialess == [[2]] && isempty(divider.algebraic)
        @test_throws ArgumentError transientsolve(divider, (0.0, 1e-9); dt = 1e-12,
            initialstate = transientstate(divider; voltage = [50e-6, 0.0]))
        d2 = transientsolve(divider, (0.0, 1e-9); dt = 1e-12,
            initialstate = transientstate(divider; voltage = [100e-6/3, 50e-6/3]))
        @test d2.voltage[1, :] ≈ fill(100e-6/3, size(d2.voltage, 2)) rtol=1e-8
        # a resistor inside a capacitive island carries no current along
        # the island's direction, so the direction is algebraic and its
        # rate is read by the differentiated equation: equal voltages on
        # the two nodes would ring between the inductors, opposite ones
        # are the loop's own mode
        loop = transientproblem([("C1", "1", "2", 1e-12), ("R1", "1", "2", 50.0),
            ("L1", "1", "0", 1e-9), ("L2", "2", "0", 1e-9)])
        @test loop.inertialess == [[1, 2]] && loop.algebraic == [[1, 2]]
        @test_throws ArgumentError transientsolve(loop, (0.0, 10e-12); dt = 1e-12,
            initialstate = transientstate(loop; voltage = [1e-6, 1e-6]))
        ls = transientsolve(loop, (0.0, 100e-12); dt = 1e-12, record = :states,
            initialstate = transientstate(loop; voltage = [1e-6, -1e-6]))
        @test maximum(abs.(ls.rate[1, :] .+ ls.rate[2, :])) < 1e-6*maximum(abs.(ls.rate[1, :]))
        @test abs(ls.rate[1, 2] - ls.rate[1, 1]) < 0.05*abs(ls.rate[1, 1])
        # two islands joined by a resistor are one algebraic direction, and
        # none once a resistor reaches ground
        twin = [("C1", "1", "2", 1e-12), ("C2", "3", "4", 1e-12), ("R1", "2", "3", 50.0),
            ("L1", "1", "0", 1e-9), ("L2", "2", "0", 1e-9), ("L3", "3", "0", 1e-9), ("L4", "4", "0", 1e-9)]
        twins = transientproblem(twin)
        @test twins.inertialess == [[1, 2], [3, 4]] && twins.algebraic == [[1, 2, 3, 4]]
        @test_throws ArgumentError transientsolve(twins, (0.0, 10e-12); dt = 1e-12,
            initialstate = transientstate(twins; voltage = [1e-6, 1e-6, 1e-6, 1e-6]))
        grounded = transientproblem(vcat(twin, [("R2", "4", "0", 50.0)]))
        @test grounded.inertialess == [[1, 2], [3, 4]] && isempty(grounded.algebraic)
        @test transientsolve(grounded, (0.0, 10e-12); dt = 1e-12,
            initialstate = transientstate(grounded; voltage = [1e-6, 1e-6, 1e-6, 0.0])) isa JC.TransientSolution
        # a node no element connects to ground has a free flux offset, and
        # only such a node gets a gauge row
        island = transientproblem([("C1", "1", "0", 1e-12), ("L2", "2", "3", 1e-9), ("C3", "2", "3", 1e-12)])
        @test island.gaugeindices == [2]
        @test isempty(transientproblem(rc).gaugeindices)
        @test_throws ArgumentError transientproblem([("C1", "1", "0", 1e-12 + 1e-15im)])
        @test_throws ArgumentError transientproblem([("R1", "1", "0", FrequencyDependent(w -> 50.0))])
        @test_throws ArgumentError transientproblem([("L1", "1", "0", 0.0)])
        @test_throws ArgumentError transientproblem(rc; sources = [TransientSource(9, 0.0)])
        @test_throws ArgumentError transientproblem(rc; sources = [TransientSource("R1", 0.0)])
        @test_throws ArgumentError transientsolve(transientproblem(rc), (0.0, 1e-9); dt = 0.0)
        @test_throws ArgumentError transientsolve(transientproblem(rc), (0.0, 1e-9); dt = 1e-12, record = :phases, saveevery = 2)
        @test_throws ArgumentError transientsolve(transientproblem(rc), (0.0, 1e-9); dt = 1e-12, maxsteps = 5)
    end

    @testset "a rejected correction is retried from its base point" begin
        # the Newton engine of the step on a scalar junction residual
        # `x + 10 sin(x) - 1` with a kept Jacobian of 1 from a phase of
        # pi/2: the first correction, to 1, is rejected, and the retry
        # must build the fresh Jacobian 11 and the correction at the base
        # point 0, not at the rejected trial, so that one retry converges
        jac = Ref(1.0)
        base = Ref(0.0)
        res(norms, r, y) = (r[1] = y[1] + 10sin(y[1]) - 1; norms[1] = abs(r[1]); nothing)
        refresh = () -> (jac[] = 1 + 10cos(base[]); nothing)
        x = [0.0]
        baseres = (norms, r, y) -> (base[] = y[1]; res(norms, r, y))
        solve = (c, r) -> (c[1] = r[1]/jac[]; (false, 0))
        out = JC.newtonsolve!(x, [0.0], [0.0], [0.0], [0.0], baseres, res, refresh, solve,
            [1e-12], 15, false, true, false, JC.NewtonWork(JC.CPU(), 1))
        converged, fresh, corrections, factorizations, retries = out
        @test converged && fresh && retries == 1 && factorizations == 1
        @test x[1] + 10sin(x[1]) ≈ 1 atol=1e-12
        @test corrections <= 6
        # the stepping rule reaches the same discrete solution whatever
        # the path of its factorizations
        circuit = [("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0), ("C1", "1", "0", 1e-15),
            ("Lj1", "1", "0", 1e-10), ("L1", "1", "0", 1e-9)]
        prob = transientproblem(circuit; sources = [TransientSource(1, t -> 2e-6*sinpi(2*5e9*t))])
        sol = transientsolve(prob, (0.0, 1e-9); dt = 2e-12)
        tight = transientsolve(prob, (0.0, 1e-9); dt = 2e-12, rtol = 1e-12, atol = 1e-13, maxiters = 40)
        @test sol.finalflux ≈ tight.finalflux rtol=1e-6
        @test sol.voltage ≈ tight.voltage rtol=1e-6 atol=1e-12
    end

    @testset "the step Jacobian, the tangent and the adjoint" begin
        rng = MersenneTwister(45)
        circuit = vcat(rc, [("Lj1", "1", "0", 1e-9)])
        drive(t) = 0.12e-6*sinpi(2*3e9*t) + 0.01e-6*sinpi(2*1.3e9*t)
        perturb(t) = 0.02e-6*sinpi(2*1.7e9*t)
        prob = transientproblem(circuit; sources = [TransientSource(1, drive)])
        # the assembled step Jacobian against a finite difference of the
        # step residual, through the package's plan at one mode
        sys = JC.transientsystem(prob, 2e-12, Trapezoidal(), JC.CPU(), KLUfactorization())
        n = length(prob)
        x, d = randn(rng, n), randn(rng, n)
        phi = zeros(length(sys.lmolj)); jwork = similar(phi)
        JC.junctionphases!(phi, sys, x)
        JC.stepjacobian!(sys, phi, nothing)
        J = copy(sys.jacobian)
        r = (y -> (res = zeros(n); JC.stepresidual!(res, sys, y, phi, zeros(n), jwork, zeros(n), zeros(n)); res))
        h = 1e-6
        @test J*d ≈ (r(x + h*d) - r(x - h*d))/(2h) rtol=1e-7
        @test J ≈ transpose(J)
        for method in (Trapezoidal(), BackwardEuler())
            sol = transientsolve(prob, (0.0, 1e-9); dt = 2e-12, method, record = :phases, rtol = 1e-12)
            currents = reshape(perturb.(sol.times), 1, :)
            tangent = transienttangent(sol, currents)
            weights = randn(rng, 1, length(sol.times))
            for quantity in (:voltage, :incident, :outgoing)
                adj = transientadjoint(sol, weights; quantity)
                @test sum(weights .* getproperty(tangent, quantity)) ≈ sum(adj.currents .* currents) rtol=1e-10 atol=1e-15
            end
            eps = 1e-3
            plus = transientsolve(transientproblem(circuit; sources = [TransientSource(1, t -> drive(t) + eps*perturb(t))]),
                (0.0, 1e-9); dt = 2e-12, method, rtol = 1e-12)
            minus = transientsolve(transientproblem(circuit; sources = [TransientSource(1, t -> drive(t) - eps*perturb(t))]),
                (0.0, 1e-9); dt = 2e-12, method, rtol = 1e-12)
            @test tangent.voltage ≈ (plus.voltage - minus.voltage)/(2eps) rtol=2e-6
            # the initial state enters too; the RC junction node has a
            # capacitance so any initial state is consistent
            dx0, dv0 = [0.1], [0.05]
            initial = transienttangent(sol, zero(currents); initialstate = (dx0, dv0))
            adj = transientadjoint(sol, weights; quantity = :voltage)
            @test sum(weights .* initial.voltage) ≈ dot(adj.initialflux, dx0) + dot(adj.initialrate, dv0) rtol=1e-10
        end
    end

    @testset "many tones, and their demodulation" begin
        # a linear RC has an analytic transfer function per tone; the
        # transient carries all the tones in one state
        frequencies = [1e9, 2.3e9, 3.7e9, 4.2e9, 5.1e9]
        amplitudes = [1.0, 0.08, 0.06, 0.04, 0.02]*1e-6
        current(t) = sum(a*cospi(2f*t) for (a, f) in zip(amplitudes, frequencies))
        prob = transientproblem(rc; sources = [TransientSource(1, current)])
        sol = transientsolve(prob, (0.0, 12e-9); dt = 1e-12)
        for (f, a) in zip(frequencies, amplitudes)
            measured = transientdemodulate(sol, 1, f; quantity = :voltage, window = t -> t >= 2e-9 ? 1.0 : 0.0)
            @test measured ≈ 50a/(1 + 2pi*im*f*50e-12) rtol=0.008
        end
        @test length(prob) == 1
    end

    @testset "a floating junction and port, several drives" begin
        circuit = Circuit([
            ("p1", "1", "0", Port(1)), ("p2", "2", "1", Port(2; Z0 = 75.0)),
            ("c1", "1", "0", Capacitor(1e-12)), ("c2", "2", "0", Capacitor(1e-12)),
            ("jj", "2", "1", JosephsonJunction(1e-9))])
        drive(t) = 0.1e-6*sinpi(2*3e9*t)
        prob = transientproblem(circuit; sources = [TransientSource(1, drive)])
        sol = transientsolve(prob, (0.0, 1e-9); dt = 2e-12, record = :phases, rtol = 1e-12)
        rng = MersenneTwister(3)
        currents = [1e-8*sinpi(2*1.1e9*t + p) - 1e-8*sinpi(p) for p in 1:2, t in sol.times]
        response = transienttangent(sol, currents)
        weights = randn(rng, size(currents))
        adj = transientadjoint(sol, weights; quantity = :outgoing)
        @test sum(weights .* response.outgoing) ≈ sum(adj.currents .* currents) rtol=1e-10
        eps = 1e-3
        function loaded(sign)
            sources = [TransientSource(1, t -> drive(t) + sign*eps*1e-8*(sinpi(2*1.1e9*t + 1) - sinpi(1))),
                TransientSource(2, t -> sign*eps*1e-8*(sinpi(2*1.1e9*t + 2) - sinpi(2)))]
            transientsolve(transientproblem(circuit; sources), (0.0, 1e-9); dt = 2e-12, rtol = 1e-12)
        end
        plus, minus = loaded(1), loaded(-1)
        @test response.outgoing ≈ (plus.outgoing - minus.outgoing)/(2eps) rtol=2e-6
        @test_throws ArgumentError transientadjoint(plus, weights)
    end

    @testset "the Gauss-Legendre rule: fourth order, the responses, the reuse" begin
        # the driven RC against a fine reference: the error falls sixteen
        # fold per halving where the trapezoidal rule's falls four fold
        rc = [("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0), ("C1", "1", "0", 1e-12)]
        smooth(t) = 1e-6*sinpi(t/1e-9)^2*sinpi(2*4e9*t)
        prob = transientproblem(rc; sources = [TransientSource(1, smooth)])
        ref = transientsolve(prob, (0.0, 1e-9); dt = 1e-9/8192, method = GaussLegendre())
        errors = [maximum(abs.(transientsolve(prob, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre()).voltage[1, :] .-
            ref.voltage[1, 1:8192÷n:end])) for n in (64, 128)]
        @test 14 < errors[1]/errors[2] < 18
        # the lossless LC keeps its energy at twenty samples per period
        L, C, V = 1e-9, 1e-12, 1e-6
        lc = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)),
            ("c", "1", "0", Capacitor(C)), ("l", "1", "0", Inductor(L))]))
        period = 2pi*sqrt(L*C)
        sol = transientsolve(lc, (0.0, 3period); dt = period/20, method = GaussLegendre(),
            initialstate = transientstate(lc; voltage = [V]), record = :states)
        v = sol.rate[1, :] .* JC.phi0
        energy = C .* v .^ 2 ./ 2 .+ (sol.flux[1, :] .* JC.phi0) .^ 2 ./ (2L)
        @test maximum(abs.(energy ./ energy[1] .- 1)) < 1e-12
        @test maximum(abs.(v .- V .* cos.(sol.times ./ sqrt(L*C)))) < 3e-4V
        @test sol.stats.factorizations == 1
        # the coupled tanks through their auxiliary rows
        coupled = transientproblem([("C1", "1", "0", C), ("C2", "2", "0", C),
            ("L1", "1", "0", L), ("L2", "2", "0", L), ("K1", "L1", "L2", 0.3)])
        s2 = transientsolve(coupled, (0.0, 2period); dt = period/40, method = GaussLegendre(),
            initialstate = transientstate(coupled; voltage = [V, 0]), record = :states)
        wp, wm = 1/sqrt(L*C*1.3), 1/sqrt(L*C*0.7)
        @test maximum(abs.(s2.rate[1, :] .* JC.phi0 .- V/2 .* (cos.(wp .* s2.times) .+ cos.(wm .* s2.times)))) < 2e-5V
        # the pumped junction against harmonic balance at twenty samples per
        # period, where the trapezoidal rule is off by more than its value;
        # one complex factorization serves the whole solve
        circuit = [("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0), ("C1", "1", "2", 100e-15),
            ("Lj1", "2", "0", 1e-9), ("C2", "2", "0", 1e-12)]
        fp, ip = 4.75e9, 0.00565e-6
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        pump(t) = 2ip*ramp(t)*cospi(2fp*t)
        pa = transientproblem(circuit; sources = [TransientSource(1, pump)])
        hb = hbnlsolve((2pi*fp,), (10,), [(mode = (1,), port = 1, current = ip)], circuit,
            Dict{Symbol,Float64}(); keyedarrays = false)
        expected = 2im*2pi*fp*JC.phi0*hb.nodeflux[1]
        window(t) = 200e-9 <= t <= 300e-9 ? sinpi((t - 200e-9)/100e-9)^2 : 0.0
        gauss = transientsolve(pa, (0.0, 300e-9); dt = 10e-12, method = GaussLegendre())
        @test transientdemodulate(gauss, 1, fp; quantity = :voltage, window) ≈ expected rtol=2e-2
        @test gauss.stats.factorizations == 1
        trap = transientsolve(pa, (0.0, 300e-9); dt = 10e-12)
        @test !isapprox(transientdemodulate(trap, 1, fp; quantity = :voltage, window), expected; rtol = 0.5)
        # the tangent against finite differences and the adjoint against the
        # tangent, on the full stage equations with their two stiffnesses
        loaded = Circuit([("p1", "1", "0", Port(1)), ("p2", "2", "1", Port(2; Z0 = 75.0)),
            ("c1", "1", "0", Capacitor(1e-12)), ("c2", "2", "0", Capacitor(1e-12)),
            ("jj", "2", "1", JosephsonJunction(1e-9)), ("l", "2", "0", Inductor(2e-9))])
        drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(2*3e9*t)
        probe(t, p) = 1e-8*sinpi(2*1.1e9*t + p)
        lp = transientproblem(loaded; sources = [TransientSource(1, drive), TransientSource(2, 2e-8)])
        dt, T = 4e-12, 0.5e-9
        rec = transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12, method = GaussLegendre())
        @test size(rec.phases) == (1, 2, length(rec.times))
        currents = [probe(t, p) for p in 1:2, t in rec.times]
        tg = transienttangent(rec, currents)
        eps = 1e-4
        plus = transientproblem(loaded; sources = [TransientSource(1, t -> drive(t) + eps*probe(t, 1)),
            TransientSource(2, t -> 2e-8 + eps*probe(t, 2))])
        minus = transientproblem(loaded; sources = [TransientSource(1, t -> drive(t) - eps*probe(t, 1)),
            TransientSource(2, t -> 2e-8 - eps*probe(t, 2))])
        sp = transientsolve(plus, (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
        sm = transientsolve(minus, (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
        @test tg.outgoing ≈ (sp.outgoing .- sm.outgoing) ./ (2eps) rtol=1e-6
        weights = [cospi(2*1.7e9*t + p) for p in 1:2, t in rec.times]
        ad = transientadjoint(rec, weights)
        @test sum(weights .* tg.outgoing) ≈ sum(ad.currents .* currents) rtol=1e-10
        dx0, dv0 = [0.3, -0.2], [1e9, 2e9]
        tg0 = transienttangent(rec, zeros(2, length(rec.times)); initialstate = (dx0, dv0))
        @test sum(weights .* tg0.outgoing) ≈ dot(ad.initialflux, dx0) + dot(ad.initialrate, dv0) rtol=1e-9
        # a sink receives every column of the currents once it is final,
        # under both rules, and the adjoint then stores none
        for (s, a) in ((rec, ad), (transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12), nothing))
            stored = isnothing(a) ? transientadjoint(s, weights) : a
            received = zeros(size(stored.currents))
            seen = Int[]
            streamed = transientadjoint(s, weights; sink = (k, values) -> (push!(seen, k); received[:, k] .= values[:, 1]; nothing))
            @test isnothing(streamed.currents)
            @test seen == collect(length(s.times):-1:1)
            @test received ≈ stored.currents rtol=1e-12
            @test streamed.initialflux ≈ stored.initialflux
        end
        # the responses need the stage phases, and the reuse carries the
        # complex factorization across the solve and the responses
        bare = transientsolve(lp, (0.0, T); dt, method = GaussLegendre())
        @test_throws ArgumentError transienttangent(bare, currents)
        reuse = TransientReuse()
        r1 = transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12, method = GaussLegendre(), reuse)
        @test reuse.system.method isa GaussLegendre
        @test transienttangent(r1, currents; reuse).outgoing ≈ tg.outgoing rtol=1e-10
        @test_throws ArgumentError transientsolve(lp, (0.0, T); dt, method = GaussLegendre(), linearsolver = GMRES())
    end

    @testset "the algebraic directions of the Gauss-Legendre rule" begin
        # a junction alone across an unterminated port: the flux is the
        # constraint's own solution and the voltage its derivative, which
        # the projected endpoint holds to the Newton tolerance and the
        # reconstructed rate recovers at third order, where the plain
        # rule's rate would not converge
        L, w, a = 1e-9, 2pi*1e9, 0.2
        jj = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)),
            ("jj", "1", "0", JosephsonJunction(L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        @test jj.inertialess == [[1]] && jj.algebraic == [[1]]
        errors = map((20, 40, 80)) do n
            s = transientsolve(jj, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
            phase = asin.(a .* (1 .- cos.(w .* s.times)))
            voltage = JC.phi0 .* (a*w .* sin.(w .* s.times)) ./ sqrt.(1 .- (a .* (1 .- cos.(w .* s.times))) .^ 2)
            (maximum(abs.(s.flux[1, :] .- phase)), maximum(abs.(s.voltage[1, :] .- voltage))/(JC.phi0*w))
        end
        @test all(e -> e[1] < 1e-10, errors)
        @test 7 < errors[1][2]/errors[2][2] < 9 && 7 < errors[2][2]/errors[3][2] < 9
        @test errors[3][2] < 2e-6
        # The classification reads the rate system of the equations, not a
        # graph of the ports. An open block on the junction's node adds no
        # conductance: the node stays algebraic and the solution and its
        # order are unchanged (a port taken for a resistor lost the
        # constraint and left a 0.42 error at every step). A short block
        # grounds the node: no direction, and no voltage. A through to a
        # second node with an inductor joins the two into one algebraic
        # direction, the junction and the inductor in parallel; and a
        # rational inductor block on the bare junction's node equals the
        # explicit inductor there.
        jjo = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
            ("open", "1", ScatteringParameters(ones(1, 1)))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        @test jjo.algebraic == [[1]] && jjo.inertialess == [[1], [2]]
        @test size(jjo.constraints) == (1, 2) && abs(jjo.constraints[1, 1]) > 0.1 && abs(jjo.constraints[1, 2]) > 0.01
        oerrors = map((20, 40, 80)) do n
            s = transientsolve(jjo, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
            voltage = JC.phi0 .* (a*w .* sin.(w .* s.times)) ./ sqrt.(1 .- (a .* (1 .- cos.(w .* s.times))) .^ 2)
            maximum(abs.(s.voltage[1, :] .- voltage))/(JC.phi0*w)
        end
        @test all(k -> isapprox(oerrors[k], errors[k][2]; rtol = 1e-6), 1:3)
        jjs = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
            ("short", "1", ScatteringParameters(-ones(1, 1)))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        @test isempty(jjs.algebraic) && jjs.inertialess == [[1], [2]]
        ss = transientsolve(jjs, (0.0, 1e-9); dt = 1e-9/40, method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
        @test maximum(abs, ss.voltage) < 1e-20
        jjt = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
            ("through", "1", "2", ScatteringParameters([0.0 1.0; 1.0 0.0])), ("l", "2", "0", Inductor(L))]);
            sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        jjp = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
            ("l", "1", "0", Inductor(L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        @test jjt.algebraic == [[1, 2]] && jjp.algebraic == [[1]]
        st = transientsolve(jjt, (0.0, 1e-9); dt = 1e-9/40, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
        sp = transientsolve(jjp, (0.0, 1e-9); dt = 1e-9/40, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
        @test st.voltage ≈ sp.voltage rtol=1e-8
        @test st.flux[1, :] ≈ sp.flux[1, :] rtol=1e-8
        @test st.flux[2, :] ≈ sp.flux[1, :] rtol=1e-8
        ra = 2*50.0/L
        rind = RationalScattering(fill(-ra, 1, 1), reshape([1.0, -1.0], 1, 2), reshape(-ra .* [1.0, -1.0], 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
        jjr = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
            ("lb", "1", "2", rind), ("l2", "2", "0", Inductor(L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        jjp2 = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("jj", "1", "0", JosephsonJunction(L)),
            ("l", "1", "0", Inductor(2L))]); sources = [TransientSource(1, t -> JC.phi0/L*a*(1 - cos(w*t)))])
        # the block's ports are open at infinite frequency, so both nodes
        # are algebraic, their constraints carrying the block's resting
        # waves; the block's states are stepped from the stage rates along
        # those directions, so the agreement converges at second order
        @test jjr.algebraic == [[1], [2]]
        rdiffs = map((160, 320, 640)) do n
            sr = transientsolve(jjr, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
            sp2 = transientsolve(jjp2, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), record = :states, rtol = 1e-12, atol = 1e-13)
            maximum(abs, sr.voltage .- sp2.voltage)/maximum(abs, sp2.voltage)
        end
        @test 3.5 < rdiffs[1]/rdiffs[2] < 4.5 && 3.5 < rdiffs[2]/rdiffs[3] < 4.5 && rdiffs[3] < 2e-5
        # a driven inductor without capacitance: a linear constraint, met by
        # one correction, and the same third order in the voltage
        ind = transientproblem(Circuit([("p", "1", "0", Port(1; termination = nothing)), ("l", "1", "0", Inductor(L))]);
            sources = [TransientSource(1, t -> 1e-6*(1 - cos(w*t)))])
        @test ind.algebraic == [[1]]
        verrors = map((20, 40)) do n
            s = transientsolve(ind, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), rtol = 1e-12, atol = 1e-13)
            maximum(abs.(s.voltage[1, :] .- L*1e-6*w .* sin.(w .* s.times)))/(L*1e-6*w)
        end
        @test 7 < verrors[1]/verrors[2] < 9 && verrors[2] < 3e-5
        # the index one unknowns of the endpoint are read from their
        # equations: a terminated port on an inductor, whose node has no
        # capacitance, keeps fourth order in its voltage, where the
        # rule's endpoint formula alone would give second, and so does a
        # block's port current; the tangent and the adjoint follow the
        # reading
        indport = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:l, 1, 2, Inductor(1e-9)), (:c, 2, 0, Capacitor(1e-12)),
            (:jj, 2, 0, JosephsonJunction(1e-9))])
        blockport = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:att, 1, 2, ScatteringParameters([0.0 0.5; 0.5 0.0]; zref = 50.0)),
            (:l, 2, 3, Inductor(1e-9)), (:c, 3, 0, Capacitor(1e-12)), (:jj, 3, 0, JosephsonJunction(1e-9))])
        pulse(t) = t <= 0 ? 0.0 : 1e-6*sinpi(t/1e-9)^2*sinpi(2*4e9*t)
        for (c, resistive, aux) in ((indport, [[1]], Int[]), (blockport, [[1], [2]], [4, 5]))
            p = transientproblem(c; sources = [TransientSource(1, pulse)])
            pr = JC.transientsystem(p, 1e-12, GaussLegendre(), JC.CPU(), KLUfactorization()).gauss.projection
            @test pr.readrows == resistive && pr.auxrows == aux && isempty(pr.directions)
            ref = transientsolve(p, (0.0, 1e-9); dt = 1e-9/4096, method = GaussLegendre(), rtol = 1e-13, atol = 1e-15)
            errs = [maximum(abs.(transientsolve(p, (0.0, 1e-9); dt = 1e-9/n, method = GaussLegendre(), rtol = 1e-13, atol = 1e-15).voltage[1, :] .-
                ref.voltage[1, 1:4096÷n:end])) for n in (64, 128)]
            @test 12 < errs[1]/errs[2] < 20
        end
        ip = transientproblem(indport; sources = [TransientSource(1, pulse)])
        irec = transientsolve(ip, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre())
        iprobe(t) = 1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + 1)
        ieps = 1e-4
        icurrents = [iprobe(t) for _ in 1:1, t in irec.times]
        itg = transienttangent(irec, icurrents)
        ishift(s) = transientproblem(indport; sources = [TransientSource(1, t -> pulse(t) + s*ieps*iprobe(t))])
        isp = transientsolve(ishift(1), (0.0, 0.5e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre())
        ism = transientsolve(ishift(-1), (0.0, 0.5e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre())
        @test itg.voltage ≈ (isp.voltage .- ism.voltage) ./ (2ieps) rtol=1e-5
        iweights = [cospi(2*1.7e9*t) for _ in 1:1, t in irec.times]
        iad = transientadjoint(irec, iweights; quantity = :voltage)
        @test sum(iweights .* itg.voltage) ≈ sum(iad.currents .* icurrents) rtol=1e-9
        icps = transientsolve(ip, (0.0, 0.5e-9); dt = 2e-12, record = :checkpoints, checkpointevery = 16, rtol = 1e-12, method = GaussLegendre())
        @test transientadjoint(icps, iweights; quantity = :voltage).currents ≈ iad.currents rtol=1e-7
        # the coupled tanks and the lossless LC project nothing
        @test isnothing(JC.transientsystem(transientproblem([("C1", "1", "0", 1e-12), ("C2", "2", "0", 1e-12),
            ("L1", "1", "0", L), ("L2", "2", "0", L), ("K1", "L1", "L2", 0.3)]), 1e-12, GaussLegendre(), JC.CPU(),
            KLUfactorization()).gauss.projection)
        # the tangent against finite differences and the adjoint against
        # the tangent through the projection, on the record, on the
        # checkpoints and on a batch: a junction from a driven capacitive
        # node to a node with only an inductor
        c = Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "0", Capacitor(1e-12)), ("jj", "1", "2", JosephsonJunction(1e-9)),
            ("l2", "2", "0", Inductor(2e-9)), ("p2", "2", "0", Port(2; termination = nothing))])
        drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(2*3e9*t)
        probe(t, k) = 1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k)
        lp = transientproblem(c; sources = [TransientSource(1, drive), TransientSource(2, t -> 0.0)])
        sys = JC.transientsystem(lp, 2e-12, GaussLegendre(), JC.CPU(), KLUfactorization())
        @test sys.gauss.projection.directions == [[2]] && sys.gauss.projection.pj == [1]
        dt, T = 2e-12, 0.5e-9
        rec = transientsolve(lp, (0.0, T); dt, record = :phases, rtol = 1e-12, method = GaussLegendre())
        @test size(rec.endphases) == (1, length(rec.times))
        currents = [probe(t, k) for k in 1:2, t in rec.times]
        tg = transienttangent(rec, currents)
        eps = 1e-4
        shifted(s) = transientproblem(c; sources = [TransientSource(1, t -> drive(t) + s*eps*probe(t, 1)),
            TransientSource(2, t -> s*eps*probe(t, 2))])
        sp = transientsolve(shifted(1), (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
        sm = transientsolve(shifted(-1), (0.0, T); dt, rtol = 1e-12, method = GaussLegendre())
        @test tg.outgoing ≈ (sp.outgoing .- sm.outgoing) ./ (2eps) rtol=1e-5
        @test tg.voltage ≈ (sp.voltage .- sm.voltage) ./ (2eps) rtol=3e-5
        weights = [cospi(2*1.7e9*t + k) for k in 1:2, t in rec.times]
        ad = transientadjoint(rec, weights)
        @test sum(weights .* tg.outgoing) ≈ sum(ad.currents .* currents) rtol=1e-9
        # every direction converges against its own right hand side: a
        # direction sixteen orders weaker than its companion propagates as
        # it does alone
        dirs = zeros(2, length(rec.times), 2)
        dirs[1, :, 1] .= currents[1, :]
        dirs[2, :, 2] .= 1e-16 .* currents[2, :]
        both = transienttangent(rec, dirs)
        @test both.outgoing[:, :, 2] ≈ transienttangent(rec, dirs[:, :, 2]).outgoing rtol=1e-8
        # a current given at the grid and the stage times is read as it
        # is, and agrees with the stencil's reading of the smooth grid
        # current to the stencil's interpolation
        cs = JC.gausscoefficients().c
        staged = [probe(t + (s == 1 ? 0.0 : cs[s - 1]*dt), k) for k in 1:2, s in 1:3, t in rec.times, _ in 1:1]
        @test transienttangent(rec, staged).outgoing ≈ tg.outgoing rtol=1e-5
        @test transienttangent(rec, staged).outgoing ≈ (sp.outgoing .- sm.outgoing) ./ (2eps) rtol=1e-5
        cps = transientsolve(lp, (0.0, T); dt, record = :checkpoints, checkpointevery = 16, rtol = 1e-12, method = GaussLegendre())
        @test transienttangent(cps, currents).outgoing ≈ tg.outgoing rtol=1e-7
        @test transientadjoint(cps, weights).currents ≈ ad.currents rtol=1e-7
        half = transientproblem(lp; sources = [TransientSource(1, t -> drive(t)/2), TransientSource(2, t -> 0.0)])
        b = transientsolve([lp, half], (0.0, T); dt, record = :phases, rtol = 1e-12)
        @test b.voltage[:, :, 1] == rec.voltage
        @test size(b.endphases) == (1, length(rec.times), 2)
        tb = transienttangent(b, currents)
        @test tb.outgoing[:, :, 1] == tg.outgoing
        @test tb.outgoing[:, :, 2] ≈ transienttangent(b[2], currents).outgoing rtol=1e-9
        @test transientadjoint(b, weights).currents[:, :, 1] == ad.currents
    end

    @testset "a batch of drive conditions" begin
        # the amplifier under three pump amplitudes as one solve: each
        # member equals its own solve, the responses run on a member, and
        # the problems must share the circuit and the targets
        circuit = [("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0), ("C1", "1", "2", 100e-15),
            ("Lj1", "2", "0", 1e-9), ("C2", "2", "0", 1e-12)]
        fp = 4.75e9
        ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
        base = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
        make(ip) = transientproblem(base; sources = [TransientSource(1, let ip = ip; t -> 2ip*ramp(t)*cospi(2fp*t); end)])
        problems = make.([0.002e-6, 0.004e-6, 0.00565e-6])
        @test problems[2].circuit === base.circuit && problems[2].matrices === base.matrices
        batch = transientsolve(problems, (0.0, 5e-9); dt = 5e-12, record = :states)
        @test batch isa TransientBatchSolution && length(batch) == 3
        @test size(batch.voltage) == (1, 1001, 3) && size(batch.phases) == (1, 2, 1001, 3)
        weights = [cospi(2*4.74e9*t) for p in 1:1, t in batch.times]
        for j in 1:3
            single = transientsolve(problems[j], (0.0, 5e-9); dt = 5e-12, method = GaussLegendre(), record = :states)
            member = batch[j]
            @test member.voltage ≈ single.voltage rtol=1e-12
            @test member.flux ≈ single.flux rtol=1e-12
            @test member.phases ≈ single.phases rtol=1e-12
            @test member.finalrate ≈ single.finalrate rtol=1e-12
            @test transientadjoint(member, weights).currents ≈ transientadjoint(single, weights).currents rtol=1e-10
        end
        @test batch.stats.factorizations == 1
        # the responses of the whole batch on one pass equal the members'
        currents = [1e-8*sinpi(2*4.7e9*t) for p in 1:1, t in batch.times]
        tb = transienttangent(batch, currents)
        ab = transientadjoint(batch, weights)
        @test size(tb.outgoing) == (1, 1001, 3) && size(ab.currents) == (1, 1001, 3)
        for j in 1:3
            @test tb.outgoing[:, :, j] ≈ transienttangent(batch[j], currents).outgoing rtol=1e-10
            @test ab.currents[:, :, j] ≈ transientadjoint(batch[j], weights).currents rtol=1e-9
            @test ab.initialflux[:, j] ≈ transientadjoint(batch[j], weights).initialflux rtol=1e-9
        end
        # a tangent measured by a sink stores no history
        measured = zeros(1, 3)
        sink = (k, v, i, o) -> (measured .+= o; nothing)
        ts = transienttangent(batch, currents; outputsink = sink)
        @test isnothing(ts.outgoing) && isnothing(ts.voltage)
        @test measured ≈ sum(tb.outgoing; dims = 2)[:, 1, :] rtol=1e-12
        @test ts.finalflux ≈ tb.finalflux
        # a record of checkpoints: the phases of every window are replayed
        # from its checkpoint, so the responses equal the full record's to
        # the Newton tolerance and nothing per time is stored
        cps = transientsolve(problems, (0.0, 5e-9); dt = 5e-12, record = :checkpoints, checkpointevery = 128)
        @test isnothing(cps.phases) && isnothing(cps.flux)
        @test cps.checkpoints.every == 128 && size(cps.checkpoints.flux) == (2, 8, 3) && size(cps.checkpoints.increment) == (2, 2, 8, 3)
        tc = transienttangent(cps, currents)
        ac = transientadjoint(cps, weights)
        @test tc.outgoing ≈ tb.outgoing rtol=1e-7
        @test ac.currents ≈ ab.currents rtol=1e-7
        @test ac.initialflux ≈ ab.initialflux rtol=1e-7
        @test transientadjoint(cps[2], weights).currents ≈ ab.currents[:, :, 2] rtol=1e-7
        auto = transientsolve(problems[1], (0.0, 5e-9); dt = 5e-12, method = GaussLegendre(), record = :checkpoints)
        @test auto.checkpoints.every == 32 && size(auto.checkpoints.flux) == (2, 32)
        @test transientadjoint(auto, weights).currents ≈ ab.currents[:, :, 1] rtol=1e-7
        @test_throws ArgumentError transientsolve(problems[1], (0.0, 1e-9); dt = 5e-12, record = :checkpoints)
        # by default only the ports are recorded
        ports = transientsolve(problems, (0.0, 1e-9); dt = 5e-12)
        @test isnothing(ports.phases) && isnothing(ports.flux) && size(ports.voltage, 3) == 3
        @test_throws ArgumentError transientadjoint(ports, weights[:, 1:201])
        # a common initial state, a state per member, and the reuse
        state = transientstate(base; voltage = [1e-6, 0.0])
        b1 = transientsolve(problems, (0.0, 1e-9); dt = 5e-12, initialstate = state)
        b2 = transientsolve(problems, (0.0, 1e-9); dt = 5e-12, initialstate = fill(state, 3))
        @test b1.finalflux == b2.finalflux
        reuse = TransientReuse()
        b3 = transientsolve(problems, (0.0, 1e-9); dt = 5e-12, reuse)
        @test reuse.factor isa JC.GaussBatchFactor && reuse.factor.ncolumns == 3
        @test transientsolve(problems, (0.0, 1e-9); dt = 5e-12, reuse).finalflux ≈ b3.finalflux rtol=1e-12
        # each condition's state is checked under its own drive, and the
        # sources a batch leaves constant must agree
        r1 = transientproblem([("P1", "1", "0", 1), ("R1", "1", "0", 50.0)]; sources = [TransientSource(1, 1e-6)])
        r2 = transientproblem(r1; sources = [TransientSource(1, 2e-6)])
        rstates = [transientstate(r1; voltage = [50e-6]), transientstate(r2; voltage = [100e-6])]
        rb = transientsolve([r1, r2], (0.0, 10e-12); dt = 1e-12, initialstate = rstates)
        @test rb.voltage[1, :, 1] ≈ fill(50e-6, 11) rtol=1e-9
        @test rb.voltage[1, :, 2] ≈ fill(100e-6, 11) rtol=1e-9
        @test_throws ArgumentError transientsolve([r1, r2], (0.0, 10e-12); dt = 1e-12, initialstate = rstates[1])
        @test_throws ArgumentError transientsolve([r1, r2], (0.0, 10e-12); dt = 1e-12, initialstate = reverse(rstates))
        two = transientproblem([("P1", "1", "0", 1), ("R1", "1", "0", 50.0), ("C1", "1", "0", 1e-12),
            ("I1", "1", "0", 1e-6), ("I2", "1", "0", 2e-6)]; sources = [TransientSource(:I1, 0.0)])
        @test_throws ArgumentError transientsolve([two, transientproblem(two; sources = [TransientSource(:I2, 0.0)])],
            (0.0, 10e-12); dt = 1e-12)
        same = transientsolve([two, transientproblem(two; sources = [TransientSource(:I1, 1e-6)])], (0.0, 1e-9); dt = 5e-12)
        @test same.voltage[1, end, 1] ≈ -100e-6 rtol=1e-6
        @test same.voltage[1, end, 2] ≈ -150e-6 rtol=1e-6
        # the contracts
        other = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
        @test_throws ArgumentError transientsolve([problems[1], other], (0.0, 1e-9); dt = 5e-12)
        twodrives = transientproblem(base; sources = [TransientSource(1, t -> 0.0), TransientSource(1, t -> 0.0)])
        @test_throws ArgumentError transientsolve([problems[1], twodrives], (0.0, 1e-9); dt = 5e-12)
        @test_throws ArgumentError transientsolve(problems, (0.0, 1e-9); dt = 5e-12, method = Trapezoidal())
        @test_throws ArgumentError transientproblem(base; sources = [TransientSource(3, 1e-6)])
        @test_throws BoundsError batch[4]
    end

    @testset "the stationary limit agrees with harmonic balance" begin
        circuit = vcat(rc, [("Lj1", "1", "0", 1e-9)])
        f, ip = 3e9, 0.12e-6
        prob = transientproblem(circuit; sources = [TransientSource(1, t -> ip*cospi(2f*t))])
        sol = transientsolve(prob, (0.0, 12e-9); dt = 0.5e-12)
        # a factorization is kept while Newton converges in one correction
        @test sol.stats.factorizations < 10
        # the iterative step, on the package's GMRES with the factorization
        # as its preconditioner, gives the same trajectory from one
        # factorization, to the Newton tolerance accumulated over the steps
        iterative = transientsolve(prob, (0.0, 12e-9); dt = 0.5e-12, linearsolver = GMRES())
        @test iterative.stats.factorizations == 1
        @test iterative.stats.kryloviterations > 0
        @test iterative.voltage ≈ sol.voltage rtol=1e-4 atol=1e-12
        @test_throws ArgumentError transientsolve(prob, (0.0, 1e-9); dt = 1e-12, linearsolver = :gmres)
        # a reuse object carries the system, its factorization and the Krylov
        # workspace to the next solve, tangent and adjoint of the problem
        reuse = TransientReuse()
        a = transientsolve(prob, (0.0, 2e-9); dt = 2e-12, record = :phases, reuse)
        b = transientsolve(prob, (0.0, 2e-9); dt = 2e-12, record = :phases, reuse)
        @test a.voltage == b.voltage
        @test reuse.system === JC.transientsystem(reuse, prob, 2e-12, Trapezoidal(), JC.CPU(), KLUfactorization())
        currents = [1e-8*sinpi(2*1.7e9*t) for _ in 1:1, t in a.times]
        weights = ones(1, length(a.times))
        @test transienttangent(a, currents; reuse).voltage ≈ transienttangent(a, currents).voltage
        @test transientadjoint(a, weights; reuse).currents ≈ transientadjoint(a, weights).currents
        g = transientsolve(prob, (0.0, 2e-9); dt = 2e-12, linearsolver = GMRES(), reuse)
        @test reuse.workspace isa JC.GMRESWorkspace
        # a different step replaces the kept system
        transientsolve(prob, (0.0, 2e-9); dt = 1e-12, reuse)
        @test reuse.system.h == 1e-12
        @test_throws ArgumentError transientsolve(prob, (0.0, 1e-9); dt = 1e-12, reuse = 3)
        # harmonic balance stores the conjugacy representatives: a physical
        # cosine of peak Ip has the positive frequency coefficient Ip/2
        hb = hbnlsolve((2pi*f,), (9,), [(mode = (1,), port = 1, current = ip/2)],
            circuit, Dict{Symbol,Float64}(); method = Newton(), keyedarrays = false)
        window(t) = 2e-9 <= t <= 12e-9 ? sinpi((t - 2e-9)/10e-9)^2 : 0.0
        for (k, m) in enumerate(hb.frequencies.modes)
            m[1] > 5 && continue
            measured = transientdemodulate(sol, 1, m[1]*f; quantity = :voltage, window)
            expected = 2im*2pi*f*m[1]*JC.phi0*hb.nodeflux[k]
            @test measured ≈ expected rtol=0.003 atol=1e-13
        end
    end
end
