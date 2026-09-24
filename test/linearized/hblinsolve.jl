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
        sonly = (; keyedarrays = false, returnQE = false, returnCM = false)
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
end
