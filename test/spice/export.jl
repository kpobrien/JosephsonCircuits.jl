using JosephsonCircuits
using Test

@testset verbose=true "exportnetlist" begin

    @testset "coupled inductor names and scattering blocks" begin
        # the K line names the inductors as their own lines name them,
        # prefix included, so WRSPICE can resolve the coupling
        c = Circuit([(:p1, "1", "0", Port(1)), (:coil1, "1", "0", Inductor(1e-9)),
            (:coil2, "2", "0", Inductor(1e-9)),
            (:k1, :coil1, :coil2, MutualInductor(0.5)),
            (:c2, "2", "0", Capacitor(1e-12)), (:r2, "2", "0", Resistor(50.0))])
        lines = split(JosephsonCircuits.exportnetlist(c, Dict()).netlist, "\n")
        @test any(l -> startswith(l, "Lcoil1 "), lines)
        @test any(l -> startswith(l, "Lcoil2 "), lines)
        @test any(l -> l == "k1 Lcoil1 Lcoil2 0.5", lines)
        # a scattering block has no SPICE element and is refused, not dropped
        blk = Circuit([(:p1, "1", "0", Port(1)),
            (:s, "1", ScatteringParameters(reshape([0.5], 1, 1)))])
        @test_throws JosephsonCircuits.ComponentNotSupportedError JosephsonCircuits.exportnetlist(blk, Dict())
    end

    @testset "sumvalues" begin
        @test_throws(
            ErrorException("unknown component type in sumvalues"),
            JosephsonCircuits.sumvalues(:V, 1.0, 4.0)
        )
    end

    @testset "a junction's shunt capacitance wherever it is listed" begin
        # the capacitors on a junction's branch are one branch of the
        # export's table whatever their order, so the capacitor listed
        # before the junction, after it, or split on both sides of it, gives
        # the same netlist
        before = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 0, Capacitor(1e-12)), (:Lj1, 1, 0, JosephsonJunction(1e-9))])
        after = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:Lj1, 1, 0, JosephsonJunction(1e-9)), (:C1, 1, 0, Capacitor(1e-12))])
        split = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 0, Capacitor(0.5e-12)), (:Lj1, 1, 0, JosephsonJunction(1e-9)),
            (:C2, 1, 0, Capacitor(0.5e-12))])
        reference = JosephsonCircuits.exportnetlist(after).netlist
        @test JosephsonCircuits.exportnetlist(before).netlist == reference
        @test JosephsonCircuits.exportnetlist(split).netlist == reference
        @test JosephsonCircuits.exportnetlist(before; jj = false).netlist ==
            JosephsonCircuits.exportnetlist(split; jj = false).netlist
    end

    @testset "the jj model's conditions" begin
        # two junctions whose critical currents differ by a factor of a
        # hundred and ten either way: WRSPICE's junction model cannot span
        # them
        pair(Lj1, Lj2) = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)),
            ("Lj1", "2", "0", JosephsonJunction(Lj1)), ("Cj1", "2", "0", Capacitor(1e-12)), ("C2", "2", "3", Capacitor(100e-15)),
            ("Lj2", "3", "0", JosephsonJunction(Lj2)), ("Cj2", "3", "0", Capacitor(1e-12))])
        @test_throws ErrorException JosephsonCircuits.exportnetlist(pair(1.0e-9, 100*1.1e-9))
        @test_throws ErrorException JosephsonCircuits.exportnetlist(pair(100*1.1e-9, 1.0e-9))
        # a junction without shunt capacitance
        unshunted = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)),
            ("C1", "1", "2", Capacitor(100e-15)), ("Lj1", "2", "0", JosephsonJunction(1e-9))])
        @test_throws(ErrorException("Cj cannot be zero in the WRSPICE JJ model."),
            JosephsonCircuits.exportnetlist(unshunted))
        # a relation other than the sinusoidal one, a SNAIL's
        snail = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:cc, 1, 2, Capacitor(100e-15)),
            (:nl, 2, 0, NonlinearInductor(1e-9, PolynomialCPR([1.0, 0.3, -1/6]))),
            (:cj, 2, 0, Capacitor(1e-12))])
        @test_throws JosephsonCircuits.ComponentNotSupportedError JosephsonCircuits.exportnetlist(snail)
        # written as their linear inductances the junctions are inductors,
        # which none of the model's conditions apply to
        lines(c) = split(JosephsonCircuits.exportnetlist(c; jj = false).netlist, "\n")
        @test "Lj2 3 0 109999.99999999999p" in lines(pair(1.0e-9, 100*1.1e-9))
        @test "Lj1 2 0 1000.0000000000001p" in lines(unshunted)
        @test "Lnl 2 0 1000.0000000000001p" in lines(snail)
    end

    @testset "parallel elements" begin
        # an inductor a mutual inductor couples keeps its own line beside a
        # parallel one, so the coupling names a line of the netlist
        g = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:la, 1, 0, Inductor(1e-9)),
            (:lb, 1, 0, Inductor(2e-9)), (:lc, 2, 0, Inductor(1e-9)),
            (:k1, :lb, :lc, MutualInductor(0.5)), (:c2, 2, 0, Capacitor(1e-12)),
            (:r2, 2, 0, Resistor(50.0))])
        glines = split(JosephsonCircuits.exportnetlist(g).netlist, "\n")
        @test "la 1 0 1000.0000000000001p" in glines
        @test "lb 1 0 2000.0000000000002p" in glines
        @test "k1 lb lc 0.5" in glines
        # the resistors of one branch combine in parallel: a port's
        # termination and a load across it
        loaded = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:rl, 1, 0, Resistor(150.0)),
            (:c1, 1, 0, Capacitor(1.0e-12))])
        @test "Rp1_termination 1 0 37.5" in
            split(JosephsonCircuits.exportnetlist(loaded).netlist, "\n")
    end

    @testset "the phase nodes are not nets" begin
        # a net named by an integer past the node count: the phase node of
        # the junction is numbered past it
        c = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:cc, 1, 3, Capacitor(100e-15)),
            (:jj, 3, 0, JosephsonJunction(1e-9)), (:cj, 3, 0, Capacitor(1e-12))])
        n = JosephsonCircuits.exportnetlist(c)
        @test !(n.junctions[1].phasenode in compile(c).nodenames)
        @test occursin("Bjj 3 0 $(n.junctions[1].phasenode) jjk ", n.netlist)
    end

    @testset "values SPICE cannot express" begin
        # a lossy capacitor has a complex value, which no SPICE element
        # takes; a value of complex type without an imaginary part is real
        jpa(Cj) = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:cc, 1, 2, Capacitor(100e-15)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(Cj))])
        @test_throws ArgumentError JosephsonCircuits.exportnetlist(jpa(1e-12*(1 - 1e-2im)))
        @test JosephsonCircuits.exportnetlist(jpa(1e-12 + 0im)).netlist ==
            JosephsonCircuits.exportnetlist(jpa(1e-12)).netlist
    end

    @testset "SPICE element names" begin
        JC = JosephsonCircuits
        # SPICE takes an element's type from the first character of its name
        # and does not accept "/" in one, so a hierarchical instance path is
        # not a name it can read
        @test JC.spicename("R1", 'R') == "R1"
        @test JC.spicename("r1", 'R') == "r1"      # SPICE is case insensitive
        @test JC.spicename("foo", 'R') == "Rfoo"
        @test JC.spicename("p1/termination", 'R') == "Rp1_termination"
        @test JC.spicename("Lj1", 'B') == "B1"     # the legacy junction name
        @test JC.spicename("jj", 'B') == "Bjj"

        # a netlist with string names is written under those names, the
        # generated termination under the port's
        named = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)),
                  ("Lj1", "2", "0", JosephsonJunction(1e-9)), ("C2", "2", "0", Capacitor(1e-12))])
        lines = split(JC.exportnetlist(named, Dict{Any,Any}()).netlist, "\n")
        @test any(startswith("RP1_termination "), lines)
        @test any(startswith("C1 "), lines)
        @test any(startswith("B1 "), lines)

        # and a typed circuit's instance paths are written as names SPICE
        # reads
        c = Circuit([:p1 => Port(1; termination = nothing),
                     :foo => Resistor(50.0), :bar => Capacitor(100e-15),
                     :baz => Inductor(1e-9)],
            [[(:p1,1),(:foo,1),(:bar,1),(:baz,1)],
             [(:p1,2),(:foo,2),(:bar,2),(:baz,2), Ground]])
        tlines = split(JC.exportnetlist(c).netlist, "\n")
        @test any(startswith("Rfoo "), tlines)
        @test any(startswith("Cbar "), tlines)
        @test any(startswith("Lbaz "), tlines)
    end

end