using JosephsonCircuits
using LinearAlgebra
using Test

# Tests of the legacy tuple netlist, kept together so the compatibility
# surface is one file on the test side as it is on the source side. When the
# legacy format is deprecated, this file goes with `src/legacyadapter.jl`.
# A test belongs here if it would stop making sense once the tuple netlist
# is gone; a test that merely uses a netlist as a convenient input does not.

@testset verbose=true "legacy tuple netlist" begin

    @testset "parsecomponenttype" begin
        @test_throws(
            ArgumentError("parsecomponenttype() currently only works for two letter components"),
            JosephsonCircuits.parsecomponenttype("BAD1",["Lj","BAD","L","C","K","I","R","P"])
        )

        @test_throws(
            ArgumentError("No matching component found in allowedcomponents."),
            JosephsonCircuits.parsecomponenttype("B1",["Lj","L","C","K","I","R","P"])
        )
    end

    @testset "checkcomponenttypes" begin
        @test_throws(
            ArgumentError("Allowed components parsing check has failed for Lj. This can happen if a two letter long component comes after a one letter component. Please reorder allowedcomponents."),
            JosephsonCircuits.checkcomponenttypes(["L","Lj","C","K","I","R","P"])
        )
    end

    @testset "legacy tuple round trip" begin
        JosephsonCircuits.@params Ipump Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Any}[
            ("P1","1","0",1),
            ("I1","1","0",Ipump),
            ("R1","1","0",Rleft),
            ("C1","1","2",Cc),
            ("Lj1","2","0",Lj),
            ("C2","2","0",Cj),
        ]
        # a netlist of tuples reaches the compiler through the typed circuit,
        # so what is checked is that the adapter reads the name prefixes as
        # component types, keeps netlist order, and passes symbolic values
        # through untouched
        psc = compile(circuit; sorting = :number)
        @test psc.componentnames == ["P1","I1","R1","C1","Lj1","C2"]
        @test psc.componenttypes == [:P, :I, :R, :C, :Lj, :C]
        @test all(isequal.(psc.componentvalues[2:end],
            [Ipump, Rleft, Cc, Lj, Cj]))
        # the port's value is its reference impedance, which a legacy
        # netlist states as the resistor across it
        @test isequal(psc.componentvalues[1], Rleft)

        # the value vector is narrowed from Vector{Any} to the tightest
        # element type the netlist admits. Every entry is a quantity, so
        # that is the parameter type here; it used to be `Real`, the join of
        # a quantity with the port number the port slot carried
        @test eltype(psc.componentvalues) ===
            JosephsonCircuits.CircuitValues.Parameter

        # node "0" is ground and sorts first, so it is node index 1
        @test psc.nodenames == ["0","1","2"]
        @test psc.nodeindices == [2 2 2 2 3 3; 1 1 1 3 1 1]

        # these node names sort the same either way
        @test JosephsonCircuits.comparestruct(psc,
            compile(circuit; sorting = :name))
    end

    @testset "legacy round trip with mutual inductors" begin
        JosephsonCircuits.@params Ipump Rleft L1v L2v C2v
        circuit = Tuple{String,String,String,Any}[
            ("P1","1","0",1),
            ("I1","1","0",Ipump),
            ("R1","1","0",Rleft),
            ("L1","1","0",L1v),
            ("K1","L1","L2",0.9),
            ("L2","2","0",L2v),
            ("C2","2","0",C2v),
        ]
        psc = compile(circuit; sorting = :number)
        @test psc.componenttypes == [:P, :I, :R, :L, :K, :L, :C]
        @test psc.inductors == [4, 6]
        @test psc.mutualinductors == [5]

        # a mutual inductor names the two branches it couples instead of
        # attaching to nodes, so its column of the node index array is empty
        @test psc.mutualinductorbranchnames == ["L1","L2"]
        @test psc.nodeindices[:,5] == [0, 0]
    end

    @testset "legacy round trip with NL and circuitdefs substitution" begin
        circuit = [
            ("P1","1","0",1),
            ("R1","1","0",:Rleft),
            ("NL1","1","0",:nlvalue),
            ("C1","1","0",:Cval),
        ]
        circuitdefs = Dict(:Rleft => 50.0, :nlvalue => 2.0, :Cval => 1e-12)
        c = Circuit(circuit, circuitdefs)
        psc = compile(c; sorting = :number)
        @test psc.componenttypes == [:P, :R, :NL, :C]
        @test psc.componentvalues[2] == 50.0
        @test psc.componentvalues[3] == 2.0
        # without substitution the symbols pass through
        psc2 = compile(Circuit(circuit); sorting = :number)
        @test psc2.componentvalues[2] == :Rleft
    end

    @testset "hbsolve equivalence: legacy tuples vs adapter vs native" begin
        circuit = [
            ("P1","1","0",1),
            ("R1","1","0",50.0),
            ("C1","1","2",100e-15),
            ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15),
        ]
        circuitdefs = Dict{Symbol,Float64}()
        ws = 2*pi*(4.5:0.05:5.0)*1e9
        wp = (2*pi*4.75001e9,)
        sources = [(mode=(1,), port=1, current=0.00565e-6)]

        sol1 = hbsolve(ws, wp, sources, (8,), (8,), circuit, circuitdefs)
        sol2 = hbsolve(ws, wp, sources, (8,), (8,), Circuit(circuit),
            circuitdefs; sorting = :number)

        # the native form states no termination of its own and keeps `r1` as
        # an ordinary device resistor, where the legacy netlist adopts `R1`
        # as the port environment. The two assign the resistor different
        # roles and so count its noise differently, but they are the same
        # electrical circuit, which is what the scattering parameters see
        native = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
             :c1 => Capacitor(100e-15), :jj => JosephsonJunction(1000e-12),
             :c2 => Capacitor(1000e-15)],
            [((:p1, 1), (:r1, 1), (:c1, 1)),
             ((:c1, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:c2, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        sol3 = hbsolve(ws, wp, sources, (8,), (8,), native)

        S1 = sol1.linearized.S((0,), 1, (0,), 1, :)
        S2 = sol2.linearized.S((0,), 1, (0,), 1, :)
        S3 = sol3.linearized.S((0,), 1, (0,), 1, :)
        @test isapprox(S1, S2; rtol = 1e-10)
        @test isapprox(S1, S3; rtol = 1e-6)
        @test maximum(abs2.(S1)) > 1.0 # it is an amplifier
    end

    # A legacy netlist has no way to say which resistor is a port's
    # environment, so two across one port was refused and must stay refused:
    # picking one silently would reinterpret a netlist which used to be an
    # error. The typed format has no such restriction, because a port states
    # its own termination.
    @testset "two resistors across a legacy port are still refused" begin
        @test_throws ArgumentError Circuit([("P1","1","0",1),
            ("R1","1","0",50.0), ("R2","1","0",50.0), ("C1","1","0",1e-12)])

        # one is fine, and becomes the port's environment
        c = Circuit([("P1","1","0",1), ("R1","1","0",50.0)])
        p = only([v for (k,v) in c.components if v isa Port])
        @test p.termination isa JosephsonCircuits.LegacyTermination
        @test p.termination.component == "R1"

        # and the message names the offenders
        err = try
            Circuit([("P1","1","0",1), ("R1","1","0",50.0),
                     ("R2","1","0",25.0)]); nothing
        catch e; e; end
        @test occursin("R1", sprint(showerror, err))
        @test occursin("R2", sprint(showerror, err))
    end
end
