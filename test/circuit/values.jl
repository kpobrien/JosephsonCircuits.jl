using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# A component value given as a closure of the mode frequency: the
# provider the circuit value carries, through a solve.
@testset verbose=true "circuit values" begin

    wp = (2*pi*4.75001*1e9,)
    src = [(mode=(1,), port=1, current=0.00565e-6)]
    ws = 2*pi*(4.5:0.1:5.0)*1e9

    @testset "FrequencyDependent provider" begin
        # a lossy frequency dependent resistor with an arbitrary law (trig
        # inside the closure -- the generality a closed expression set
        # cannot offer)
        law = w -> 50.0*(1 + 0.1*sin(w/2e10)) + im*abs(w)*1e-9
        circuit = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = FrequencyDependent(law))), ("C1", "1", "2", Capacitor(100e-15)), ("Lj1", "2", "0", JosephsonJunction(1000e-12)),
            ("C2", "2", "0", Capacitor(1000e-15))])
        out = hbnlsolve(wp, (8,), src, circuit;
            keyedarrays = false)
        @test out.solverinfo.converged
        o2 = hbsolve(ws, wp, src, (2,), (8,), circuit)
        @test all(isfinite, Array(o2.linearized.S))

        # a CONSTANT law must agree exactly with a plain numeric value
        circa = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = FrequencyDependent(w -> 50.0))), ("C1", "1", "2", Capacitor(100e-15)), ("Lj1", "2", "0", JosephsonJunction(1000e-12)),
            ("C2", "2", "0", Capacitor(1000e-15))])
        circb = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)), ("Lj1", "2", "0", JosephsonJunction(1000e-12)),
            ("C2", "2", "0", Capacitor(1000e-15))])
        Sa = hbsolve(ws, wp, src, (2,), (8,), circa).linearized.S
        Sb = hbsolve(ws, wp, src, (2,), (8,), circb).linearized.S
        @test isapprox(Array(Sa), Array(Sb), rtol = 1e-12)

        # providers combine with numbers through the value arithmetic
        circc = Circuit(Any[
            ("P1", "1", "0", Port(1; Z0 = 2*FrequencyDependent(w -> 25.0))), ("C1", "1", "2", Capacitor(100e-15)), ("Lj1", "2", "0", JosephsonJunction(1000e-12)),
            ("C2", "2", "0", Capacitor(1000e-15))])
        Sc = hbsolve(ws, wp, src, (2,), (8,), circc).linearized.S
        @test isapprox(Array(Sc), Array(Sb), rtol = 1e-12)

        # the linear-only entry point also accepts a circuit alone
        lin = hblinsolve(ws, circb; Nmodulationharmonics = (2,))
        @test all(isfinite, Array(lin.S))
    end

    @testset "a frequency dependent value at the idlers" begin
        # a series resistor and inductor written as one closure of the
        # frequency, as an expression in the frequency, and as two elements:
        # the idlers of the pumped solve lie at negative frequencies, where a
        # real element's impedance is the conjugate of its impedance at the
        # magnitude of the frequency
        R0, L0 = 5.0, 1e-9
        jpa(z...) = Circuit(Any[(:P1, 1, 0, Port(1; Z0 = 50.0)), z...,
            (:C1, 3, 2, Capacitor(100e-15)), (:Lj1, 2, 0, JosephsonJunction(1000e-12)),
            (:C2, 2, 0, Capacitor(1000e-15))])
        S(c) = Array(hbsolve(ws, wp, src, (2,), (8,), c).linearized.S)
        explicit = S(jpa((:R1, 1, 4, Resistor(R0)), (:L1, 4, 3, Inductor(L0))))
        closure = S(jpa((:Z1, 1, 3, Resistor(FrequencyDependent(w -> R0 + im*w*L0)))))
        expression = S(jpa((:Z1, 1, 3, Resistor(R0 + im*L0*FrequencyDependent(identity)))))
        @test isapprox(closure, explicit; rtol = 1e-10)
        @test isapprox(expression, explicit; rtol = 1e-10)
    end

    @testset "a parameter defined under any of its names" begin
        # a value written as a symbol or as a parameter is defined under its
        # symbol, its string or its parameter object alike, and a parameter
        # defined twice with different values is refused
        JosephsonCircuits.@params Cp
        onec(v) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(v))])
        S(c, d) = hblinsolve(ws, c, d; keyedarrays = false).S
        ref = S(onec(1e-12), Dict())
        for (value, key) in ((:Cc, :Cc), (:Cc, "Cc"),
                (:Cc, JosephsonCircuits.CircuitValues.Parameter(:Cc)),
                (Cp, :Cp), (Cp, "Cp"), (Cp, Cp))
            @test S(onec(value), Dict(key => 1e-12)) ≈ ref
        end
        @test_throws ArgumentError S(onec(:Cc), Dict(:Cc => 1e-12, "Cc" => 2e-12))
        # a string value names its parameter as a symbol does, and a name
        # defined as a frequency dependent value serves a parameter as it
        # serves a symbol
        for key in (:Cc, "Cc", JosephsonCircuits.CircuitValues.Parameter(:Cc))
            @test S(onec("Cc"), Dict(key => 1e-12)) ≈ ref
        end
        @test_throws ArgumentError S(onec("Cq"), Dict())
        fd = FrequencyDependent(w -> 1e-12)
        @test S(onec(Cp), Dict(:Cp => fd)) ≈ S(onec(:Cp), Dict(:Cp => fd)) ≈ ref
        # a name may be defined as an expression in parameters the
        # definitions give as numbers, and a value which is neither a
        # number nor symbolic is refused naming its component
        for value in (:Cc, "Cc", JosephsonCircuits.CircuitValues.Parameter(:Cc))
            @test S(onec(value), Dict(:Cc => Cp/2, :Cp => 2e-12)) ≈ ref
        end
        for d in (Dict(:Cc => nothing), Dict(:Cc => "1e-12"), Dict(:Cc => :C0))
            e = try S(onec(:Cc), d); nothing catch err; err end
            @test e isa ArgumentError && occursin("c1", sprint(showerror, e))
        end
        @test_throws ArgumentError S(onec(nothing), Dict())
    end

    @testset "a parameter defined in terms of others" begin
        # a name defined by an expression in other parameters, to any
        # depth, means inside an expression what it means as a whole
        # value; a definition which comes back to itself is refused
        JosephsonCircuits.@params Lj Lj0 Lj1
        onel(v) = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:l1, 1, 0, Inductor(v)), (:c1, 1, 0, Capacitor(1e-12))])
        S(c, d) = hblinsolve(ws, c, d; keyedarrays = false).S
        defs = Dict(Lj => Lj0/2, Lj0 => 2*Lj1, :Lj1 => 1e-9)
        ref = S(onel(2e-9), Dict())
        for v in (2*Lj, Lj + Lj1, :Lj0, Lj0, "Lj0")
            @test S(onel(v), defs) ≈ ref
        end
        @test S(onel(:Lj), defs) ≈ S(onel(1e-9), Dict())
        @test_throws ArgumentError S(onel(Lj), Dict(Lj => Lj0/2, Lj0 => 2*Lj))
    end

    @testset "two junctions on one branch" begin
        # refused naming both, since they cannot be combined into one element
        c = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:Lj1, 1, 0, JosephsonJunction(1e-9)),
            (:Lj2, 1, 0, JosephsonJunction(2e-9)), (:C1, 1, 0, Capacitor(1e-12))])
        e = try numericmatrices(c, Dict()); nothing catch err; err end
        @test e isa ArgumentError
        @test occursin("Lj1", sprint(showerror, e)) && occursin("Lj2", sprint(showerror, e))
    end

    @testset "a loss which vanishes at some frequencies of a sweep" begin
        # a capacitor lossy above 4.8 GHz and lossless below: its noise
        # channel is kept for the whole sweep in either order of the
        # frequencies, carries the loss where there is one, as the same
        # capacitor with a constant loss does, and nothing where there is not
        C0 = 1e-12
        lossyC = C0*(1 - 0.02im)
        onecap(C) = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 2, Capacitor(100e-15)), (:C2, 2, 0, Capacitor(C)),
            (:L2, 2, 0, Inductor(1e-9))])
        stepped = onecap(FrequencyDependent(w -> w > 2*pi*4.8e9 ? lossyC : C0 + 0im))
        qe(c, f) = hblinsolve(2*pi*f, c; keyedarrays = false).QE
        up = qe(stepped, [4.5e9, 5.0e9])
        down = qe(stepped, [5.0e9, 4.5e9])
        @test up[:, :, 2] ≈ qe(onecap(lossyC), [5.0e9])[:, :, 1]
        @test up[:, :, 1] ≈ qe(onecap(C0), [4.5e9])[:, :, 1]
        @test down ≈ up[:, :, [2, 1]]
    end

    @testset "a frequency dependent inductance in a pumped solve" begin
        # a linear inductor written as a closure against the same inductor
        # as a number; the inductance of a junction and the values of a
        # mutual coupling are read once rather than per mode frequency, and
        # a frequency dependent one is refused by name
        jpa(L, Lj, extra...) = Circuit(Any[(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 2, Capacitor(100e-15)), (:Lj1, 2, 3, Lj), (:L1, 3, 0, L),
            (:C2, 2, 0, Capacitor(1000e-15)), extra...])
        coupled(L, K) = jpa(L, JosephsonJunction(1000e-12),
            (:L2, 4, 0, Inductor(10e-12)), (:P2, 4, 0, Port(2; Z0 = 50.0)),
            (:K1, :L1, :L2, MutualInductor(K)))
        S(c) = Array(hbsolve(ws, wp, src, (2,), (8,), c).linearized.S)
        Lj = JosephsonJunction(1000e-12)
        @test isapprox(S(jpa(Inductor(FrequencyDependent(w -> 10e-12)), Lj)),
            S(jpa(Inductor(10e-12), Lj)); rtol = 1e-12)
        for c in (jpa(Inductor(10e-12), JosephsonJunction(FrequencyDependent(w -> 1000e-12))),
                coupled(Inductor(FrequencyDependent(w -> 10e-12)), 0.5),
                coupled(Inductor(10e-12), FrequencyDependent(w -> 0.5)))
            @test_throws ArgumentError hbsolve(ws, wp, src, (2,), (8,), c)
        end
    end
end
