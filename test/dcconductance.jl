using JosephsonCircuits
using Test

# Direct current through resistors. The harmonic balance state is periodic
# node flux, so a voltage is its time derivative and vanishes at zero
# frequency: a resistor was an open circuit at DC and a current source
# driving one could not develop I*R. The missing coordinate is the average
# voltage, which is constant on each static flux component because an
# inductor or a zero-voltage junction is a short there.

@testset verbose=true "direct current through resistors" begin

    ws = (2*pi*5e9,)
    dcsolve(c, sources) = hbnlsolve(ws, (1,), sources, c, Dict{Any,Any}();
        keyedarrays = false, dc = true, odd = true)

    @testset "a resistor obeys Ohm's law" begin
        # no inductor or junction on the driven node: either is a short at
        # DC and would hold it at zero, which is correct physics but would
        # hide the behavior under test. A capacitor is an open circuit at DC
        # and does not disturb the resistive path.
        R = 50.0; Idc = 1.0e-6
        c = Circuit([:p1 => Port(1; Z0 = R), :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2), Ground]])
        for sgn in (1, -1)
            sol = dcsolve(c, [(mode=(0,), port=1, current=sgn*Idc)])
            @test sol.solverinfo.converged
            # ground first and identically zero
            @test sol.dcnodevoltage[1] == 0
            @test isapprox(sol.dcnodevoltage[2], sgn*Idc*R; rtol = 1e-9)
        end
    end

    @testset "a resistive divider divides" begin
        # the port environment in parallel with the series pair
        Rp, Rs, Rl, Idc = 50.0, 50.0, 150.0, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Rp), :rs => Resistor(Rs), :rl => Resistor(Rl),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:rs,1),(:c1,1)], [(:rs,2),(:rl,1)],
             [(:p1,2),(:rl,2),(:c1,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=Idc)])
        @test sol.solverinfo.converged
        Vtop = Idc/(1/Rp + 1/(Rs+Rl))
        v = sort(sol.dcnodevoltage)
        @test isapprox(v[3], Vtop; rtol = 1e-9)
        @test isapprox(v[2], Vtop*Rl/(Rs+Rl); rtol = 1e-9)
    end

    @testset "an inductor shorts the resistor across it" begin
        # at DC the inductor holds both of its nodes at the same average
        # voltage, and reaching ground through one fixes it at zero, so the
        # resistor carries no direct current however hard it is driven
        c = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :l1 => Inductor(1.0e-9),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:l1,1),(:c1,1)], [(:p1,2),(:l1,2),(:c1,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=1.0e-6)])
        @test sol.solverinfo.converged
        @test isnothing(sol.dcnodevoltage) || all(iszero, sol.dcnodevoltage)
    end

    @testset "resistors carry current between floating components" begin
        # two floating nodes, each a capacitor to ground, joined by a
        # resistor and driven in and out. The bridge carries it and develops
        # I*R; the port environments are large enough to divert a part in a
        # ten million, which the tolerance allows for.
        Rb, Rbig, Idc = 100.0, 1.0e9, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Rbig), :p2 => Port(2; Z0 = Rbig),
             :rb => Resistor(Rb), :c1 => Capacitor(1e-12),
             :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:rb,1),(:c1,1)], [(:rb,2),(:p2,1),(:c2,1)],
             [(:p1,2),(:p2,2),(:c1,2),(:c2,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=Idc),
                          (mode=(0,), port=2, current=-Idc)])
        @test sol.solverinfo.converged
        v = sol.dcnodevoltage
        @test isapprox(maximum(v) - minimum(v), Idc*Rb; rtol = 1e-5)
    end

    @testset "no direct current, no voltage" begin
        # the common case: the elimination costs one projection and reports
        # nothing, and the answer is what it always was
        c = Circuit([:p1 => Port(1), :cc => Capacitor(100e-15),
                     :jj => JosephsonJunction(1000e-12),
                     :cj => Capacitor(1000e-15)],
            [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
             [(:p1,2),(:jj,2),(:cj,2), Ground]])
        sol = hbnlsolve((2*pi*4.75001e9,), (8,),
            [(mode=(1,), port=1, current=0.00565e-6)], c, Dict{Any,Any}();
            keyedarrays = false)
        @test sol.solverinfo.converged
        @test isnothing(sol.dcnodevoltage)
    end

    @testset "a floating island fixes only its differences" begin
        # A differential port has both terminals off ground, so the
        # environment it owns sits between them and the pair has no
        # conductance to ground at all. Such an island has no absolute
        # voltage: the solve pins one component and reports differences,
        # which are what is physical. The earlier cases all reach ground
        # through a port environment and so never take this path.
        Z0, Idc = 200.0, 1.0e-6
        c = Circuit(
            [:p1 => Port(1; Z0 = Z0), :ca => Capacitor(1e-12),
             :cb => Capacitor(1e-12)],
            [[(:p1,1),(:ca,1)], [(:p1,2),(:cb,1)],
             [(:ca,2),(:cb,2), Ground]])
        sol = dcsolve(c, [(mode=(0,), port=1, current=Idc)])
        @test sol.solverinfo.converged
        v = sol.dcnodevoltage
        @test isapprox(maximum(v) - minimum(v), Idc*Z0; rtol = 1e-9)
        # the reference is held at zero, so one of them is
        @test count(iszero, v) >= 2   # ground, and the pinned component
    end

    @testset "a conductance with no real DC limit is refused" begin
        # A frequency dependent resistance is evaluated at zero frequency
        # like everywhere else. If what comes back is not a finite real
        # conductance then the component has no direct current behavior to
        # use, and taking a part of it would be inventing one.
        Idc = 1.0e-6
        complexatdc = Circuit(
            [:p1 => Port(1; Z0 = 50.0),
             :r2 => Resistor(FrequencyDependent(w -> 50.0 + 10.0im)),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:r2,1),(:c1,1)], [(:p1,2),(:r2,2),(:c1,2), Ground]])
        @test_throws ArgumentError dcsolve(complexatdc,
            [(mode=(0,), port=1, current=Idc)])

        # and one whose limit is real is fine, and carries its share
        realatdc = Circuit(
            [:p1 => Port(1; Z0 = 50.0),
             :r2 => Resistor(FrequencyDependent(w -> 100.0*(1 + (w/1e11)^2))),
             :c1 => Capacitor(1.0e-12)],
            [[(:p1,1),(:r2,1),(:c1,1)], [(:p1,2),(:r2,2),(:c1,2), Ground]])
        sol = dcsolve(realatdc, [(mode=(0,), port=1, current=Idc)])
        @test sol.solverinfo.converged
        @test isapprox(sol.dcnodevoltage[2], Idc/(1/50.0 + 1/100.0);
            rtol = 1e-9)
    end
end
