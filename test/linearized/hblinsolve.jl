using JosephsonCircuits, LinearAlgebra
using Test

# The linearized sweep: what it assembles at each signal frequency and what
# it accepts.
@testset verbose=true "the linearized sweep" begin

    @testset "a lossy coupled inductor at the negative frequency modes" begin
        # a JPA whose series inductor is lossy and coupled by K to a
        # spectator inductor: as K -> 0 the coupled circuit, where the
        # inductor is an auxiliary current of the modified nodal analysis,
        # must approach the uncoupled one, where it is a node term
        function jpa(L1; K = nothing)
            c = Any[(:P1, 1, 0, Port(1; Z0 = 50.0)), (:L1, 1, 3, Inductor(L1)),
                (:C1, 3, 2, Capacitor(100e-15)),
                (:Lj1, 2, 0, JosephsonJunction(1000e-12)),
                (:C2, 2, 0, Capacitor(1000e-15)),
                (:L2, 4, 0, Inductor(1e-9)), (:C4, 4, 0, Capacitor(100e-15))]
            isnothing(K) || push!(c, (:K1, :L1, :L2, MutualInductor(K)))
            return Circuit(c)
        end
        L1 = 1e-9*(1 - 0.05im)
        ws = 2*pi*(4.5:0.25:5.0)*1e9
        close(a, b; rtol) = norm(a - b) <= rtol*norm(a)

        # the sweep about a pump, whose idlers are at negative frequencies
        wp = (2*pi*4.75001e9,)
        src = [(mode = (1,), port = 1, current = 0.00565e-6)]
        # the noise of a lossy coupled inductance is not modeled, so the
        # scattering parameters are asked for alone
        sonly = (; keyedarrays = false, returnQE = false, returnCM = false,
            returnnbar = false)
        sweep(c) = hbsolve(ws, wp, src, (2,), (8,), c; sonly...).linearized
        uncoupled = sweep(jpa(L1))
        coupled = sweep(jpa(L1; K = 1e-6))
        @test close(uncoupled.S, coupled.S; rtol = 1e-9)

        # without a pump, a signal at -w is the conjugate of one at +w
        for c in (jpa(L1), jpa(L1; K = 0.3))
            p = hblinsolve(collect(ws), c; sonly...)
            n = hblinsolve(-collect(ws), c; sonly...)
            @test close(p.S, conj.(n.S); rtol = 1e-13)
        end
        # and its noise outputs are refused, the default ones included
        @test_throws ArgumentError hblinsolve(collect(ws), jpa(L1; K = 0.3))
        @test_throws ArgumentError hblinsolve(collect(ws), jpa(1e-9; K = 0.3 - 0.01im);
            returnQE = false, returnCM = false, returnSnoise = true)

        # two pumps, whose mode (1,-2) is at a negative frequency, in the
        # pump solve itself
        wp2 = (2*pi*4.75e9, 2*pi*5.1e9)
        src2 = [(mode = (1, 0), port = 1, current = 0.004e-6),
            (mode = (0, 1), port = 1, current = 0.004e-6)]
        pump(c) = hbnlsolve(wp2, (2, 2), src2, c; keyedarrays = false)
        u = pump(jpa(L1))
        k = pump(jpa(L1; K = 1e-6))
        @test any(m -> 4.75*m[1] + 5.1*m[2] < 0, u.modes)
        @test u.solverinfo.converged && k.solverinfo.converged
        # the fluxes of the nodes of the JPA; the spectator's are first
        # order in K
        jpanodes = 1:3*length(u.modes)
        @test close(u.nodeflux[jpanodes], k.nodeflux[jpanodes]; rtol = 1e-8)
    end

    @testset "an error in the sweep is thrown as it was met" begin
        # a block whose data ends inside the sweep throws when the sweep
        # evaluates it there, from however many batches
        data(w) = abs(w) > 2pi*5.2e9 ?
            throw(DomainError(w, "the data ends at 5.2 GHz")) :
            fill(0.1 + 0.0im, 1, 1)
        c = Circuit([("p1", "1", "0", Port(1; Z0 = 50.0)),
            ("c", "1", "2", Capacitor(100e-15)),
            ("b", "2", ScatteringParameters(data; nports = 1)),
            ("L", "2", "0", Inductor(1e-9))])
        for nbatches in (1, 2)
            @test_throws DomainError hblinsolve(2pi*[4.0e9, 6.0e9], c;
                nbatches = nbatches)
        end
    end

    @testset "ports at a temperature" begin
        # each port's termination sends in the thermal field of a matched
        # load at its temperature, nbar + 1/2 at every mode; the outputs
        # are the occupation of the waves leaving the ports and the
        # quantum efficiency with every input in its state
        w = 2pi*3e9
        warm(n, T) = Port(n; termination = MatchedTermination(temperature = T))
        # equilibrium: a lossy passive two port with its terminations and
        # its loss at one temperature emits the thermal state (Bosma)
        twoport(p1, p2, loss) = Circuit([:p1 => p1, :p2 => p2,
                :c1 => Capacitor(0.3e-12), :c2 => Capacitor(0.5e-12), :loss => loss],
            [((:p1, 1), (:c1, 1), (:loss, 1)), ((:p2, 1), (:c2, 1), (:loss, 2)),
                ((:p1, 2), (:p2, 2), (:c1, 2), (:c2, 2), Ground)])
        T = 0.3
        n = thermaloccupation(w, T)
        eq = hblinsolve([w], twoport(warm(1, T), warm(2, T), Resistor(30.0; temperature = T));
            keyedarrays = false, returnVout = true)
        @test eq.porttemperatures == [T, T]
        @test eq.Vout[:, :, 1] ≈ (n + 1/2)*I atol = 1e-12
        @test eq.nbar[:, 1] ≈ [n, n] rtol = 1e-12
        @test eq.QE[:, :, 1] ≈ abs2.(eq.S[:, :, 1])/(2n + 1) rtol = 1e-12
        # a matched attenuator of transmission tau whose loss is at Tl,
        # between ports at T1 and T2: what leaves each port is what the
        # other port sends in, attenuated, and the loss, and never its own
        # termination, which the matched attenuator does not reflect; the
        # ports listed out of their order, which is their numbers'
        tau, Tl, T1, T2 = 0.4, 0.1, 0.05, 4.0
        att = ScatteringParameters([0.0 sqrt(tau); sqrt(tau) 0.0]; zref = 50.0,
            noise = ThermalEquilibrium(Tl))
        pad = Circuit([(:p2, 2, 0, warm(2, T2)), (:att, 1, 2, att), (:p1, 1, 0, warm(1, T1))])
        a = hblinsolve([w], pad; keyedarrays = false)
        @test a.porttemperatures == [T1, T2]
        n1, n2, nl = thermaloccupation(w, T1), thermaloccupation(w, T2), thermaloccupation(w, Tl)
        @test a.nbar[2, 1] ≈ tau*n1 + (1 - tau)*nl rtol = 1e-10
        @test a.nbar[1, 1] ≈ tau*n2 + (1 - tau)*nl rtol = 1e-10
        @test a.QE[2, 1, 1] ≈ tau/(2(tau*(n1 + 1/2) + (1 - tau)*(nl + 1/2))) rtol = 1e-10
        # the output covariance, requested alone, is the inputs' noise
        # scattered and the noise the circuit adds, and its diagonal the
        # occupation plus the vacuum
        v = hblinsolve([w], pad; keyedarrays = false, returnVout = true,
            returnCnoise = true, returnQE = false, returnCM = false,
            returnnbar = false)
        sigma = Diagonal([n1 + 1/2, n2 + 1/2])
        @test v.Vout[:, :, 1] ≈ v.S[:, :, 1]*sigma*v.S[:, :, 1]' + v.Cnoise[:, :, 1] rtol = 1e-12
        @test real(diag(v.Vout[:, :, 1])) .- 1/2 ≈ a.nbar[:, 1] rtol = 1e-12
        # a pumped block built from part of a warm device takes the ports it
        # leaves out into itself: the device's port 2 dropped gives the
        # device with that port's termination as a resistor at its
        # temperature
        wp = 2pi*1e9
        ws = 2pi .* collect(range(2.9e9, 3.1e9; length = 5))
        device = hbsolve(ws, (wp,), [], (1,), (2,),
            twoport(Port(1), warm(2, T2), Resistor(30.0; temperature = Tl));
            threewavemixing = true, returnCnoise = true).linearized
        blk = LinearizedScattering(device, wp; ports = [1],
            noise = NoiseCovariance(device.Cnoise))
        viablock = hbsolve([ws[3]], (wp,), [], (1,), (2,),
            Circuit([(:p1, 1, 0, Port(1)), (:b, 1, blk)]);
            threewavemixing = true, keyedarrays = false).linearized
        lumped = hbsolve([ws[3]], (wp,), [], (1,), (2,),
            twoport(Port(1), Resistor(50.0; temperature = T2), Resistor(30.0; temperature = Tl));
            threewavemixing = true, keyedarrays = false).linearized
        @test viablock.nbar[:, 1] ≈ lumped.nbar[:, 1] rtol = 1e-8
        @test viablock.QE[:, :, 1] ≈ lumped.QE[:, :, 1] rtol = 1e-8
        @test_throws ArgumentError MatchedTermination(temperature = -1.0)
    end
end
