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

    @testset "componentdictionaries errors" begin
        begin
            componenttypes = [:P, :I, :R, :L, :K, :K, :L, :C]
            nodeindices = [2 2 2 2 0 0 3 3 3; 1 1 1 1 0 0 1 1 1]
            couplings = [(5, 4, 7), (6, 4, 7)]
            @test_throws(
                DimensionMismatch("Input arrays must have the same length"),
                JosephsonCircuits.componentdictionaries(componenttypes,
                    nodeindices,couplings)
            )
        end

        begin
            componenttypes = [:P, :I, :R, :L, :K, :K, :L, :C, :C]
            nodeindices = [2 2 2 2 0 0 3 3 3; 1 1 1 1 0 0 1 1 1; 1 1 1 1 0 0 1 1 1]
            couplings = [(5, 4, 7), (6, 4, 7)]
            @test_throws(
                DimensionMismatch("The length of the first axis must be 2"),
                JosephsonCircuits.componentdictionaries(componenttypes,
                    nodeindices,couplings)
            )
        end

    end

    @testset "componentdictionaries" begin
        begin
            JosephsonCircuits.@params Ipump Rleft L1 K1 L2 C2 C3
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rleft)))
            push!(circuit,("I1", "1", "0", CurrentSource(Ipump)))
            push!(circuit,("L1", "1", "0", Inductor(L1)))
            push!(circuit,("K1", "L1", "L2", MutualInductor(K1)))
            push!(circuit,("L2", "2", "0", Inductor(L2)))
            push!(circuit,("C2", "2", "0", Capacitor(C2)))
            push!(circuit,("C3", "2", "0", Capacitor(C3)))
            circuit = Circuit(circuit)
            psc = compile(circuit)
            countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.componenttypes,psc.nodeindices,psc.couplings)

            @test isequal(countdict,Dict((:L, 1, 3) => 1, (:K, 4, 6) => 1, (:R, 1, 2) => 1, (:I, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 2, (:L, 1, 2) => 1))
            @test isequal(indexdict,Dict((:C, 1, 3, 1) => 7, (:I, 1, 2, 1) => 3, (:R, 1, 2, 1) => 2, (:L, 1, 3, 1) => 6, (:C, 1, 3, 2) => 8, (:L, 1, 2, 1) => 4, (:P, 1, 2, 1) => 1, (:K, 4, 6, 1) => 5))
        end

        begin
            JosephsonCircuits.@params Ipump Rleft L1 K1 K2 L2 C2 C3
            circuit = Any[]
            push!(circuit,("P1", "1", "0", Port(1; Z0 = Rleft)))
            push!(circuit,("I1", "1", "0", CurrentSource(Ipump)))
            push!(circuit,("L1", "1", "0", Inductor(L1)))
            push!(circuit,("K1", "L1", "L2", MutualInductor(K1)))
            push!(circuit,("K2", "L1", "L2", MutualInductor(K2)))
            push!(circuit,("L2", "2", "0", Inductor(L2)))
            push!(circuit,("C2", "2", "0", Capacitor(C2)))
            push!(circuit,("C3", "2", "0", Capacitor(C3)))
            circuit = Circuit(circuit)
            psc = compile(circuit)
            countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.componenttypes,psc.nodeindices,psc.couplings)

            @test isequal(countdict,Dict((:L, 1, 3) => 1, (:K, 4, 7) => 2, (:R, 1, 2) => 1, (:I, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 2, (:L, 1, 2) => 1))
            @test isequal(indexdict,Dict((:C, 1, 3, 1) => 8, (:I, 1, 2, 1) => 3, (:R, 1, 2, 1) => 2, (:K, 4, 7, 1) => 5, (:K, 4, 7, 2) => 6, (:L, 1, 2, 1) => 4, (:L, 1, 3, 1) => 7, (:P, 1, 2, 1) => 1, (:C, 1, 3, 2) => 9))
        end

    end

    @testset "calcCjIcmean errors" begin
        # two junctions whose critical currents differ by a factor of a
        # hundred either way: WRSPICE's junction model cannot span them
        function tables(Lj1, Lj2)
            circuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)),
                ("Lj1", "2", "0", JosephsonJunction(Lj1)), ("Cj1", "2", "0", Capacitor(1e-12)), ("C2", "2", "3", Capacitor(100e-15)),
                ("Lj2", "3", "0", JosephsonJunction(Lj2)), ("Cj2", "3", "0", Capacitor(1e-12))])
            psc = compile(circuit)
            vvn = JosephsonCircuits.componentvaluestonumber(psc.componentvalues, Dict{Any,Any}())
            countdict, indexdict = JosephsonCircuits.componentdictionaries(
                psc.componenttypes, psc.nodeindices, psc.couplings)
            return (psc.componenttypes, psc.nodeindices, vvn, psc.couplings,
                countdict, indexdict)
        end
        @test_throws ErrorException JosephsonCircuits.calcCjIcmean(tables(1.0e-9, 100*1.1e-9)...)
        @test_throws ErrorException JosephsonCircuits.calcCjIcmean(tables(100*1.1e-9, 1.0e-9)...)

        begin
            JosephsonCircuits.@params R Cc Lj Cj
            circuit = Circuit([
                ("P1", "1", "0", Port(1; Z0 = R)),
                ("C1", "1", "2", Capacitor(Cc)),
                ("Lj1", "2", "0", JosephsonJunction(Lj)),
            #    ("C2", "2", "0", Capacitor(Cj)),
                ])
            circuitdefs = Dict(
                Lj =>1000.0e-12,
                Cc => 100.0e-15,
                Cj => 1000.0e-15,
                R => 50.0)
            psc = compile(circuit)
            vvn = JosephsonCircuits.componentvaluestonumber(psc.componentvalues,circuitdefs)
            countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.componenttypes,psc.nodeindices,psc.couplings)
            @test_throws(
                ErrorException("Cj cannot be zero in the WRSPICE JJ model."),
                JosephsonCircuits.calcCjIcmean(psc.componenttypes, psc.nodeindices,
                    vvn, psc.couplings, countdict, indexdict)
                )
        end
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