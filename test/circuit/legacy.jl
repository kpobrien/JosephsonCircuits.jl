using JosephsonCircuits
using LinearAlgebra
using Test

# Tests of the deprecated input formats, the tuple netlist and the symbolic
# frequency variable, and of the deprecated entry points, kept together so
# the compatibility surface is one file on the test side as it is on the
# source side: this file goes with `src/circuit/legacy.jl` when they are
# removed. Every use of them warns, and each is captured here.

# the value of a call which warns of the deprecation
macro deprecated(ex)
    return esc(:(@test_logs (:warn,) match_mode = :any $ex))
end

# the matrices of two circuits agree, the elements and where they sit
samematrices(a, b) = a.Cnm == b.Cnm && a.Gnm == b.Gnm &&
    a.invLnm == b.invLnm && a.Ljb == b.Ljb

@testset verbose=true "deprecated tuple netlist" begin

    @testset "every entry point warns once per call" begin
        netlist = [("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        typed = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 2, Capacitor(100e-15)), (:Lj1, 2, 0, JosephsonJunction(1000e-12)),
            (:C2, 2, 0, Capacitor(1000e-15))])
        wp = (2*pi*4.75001e9,); ws = 2*pi*[4.5e9, 5.0e9]
        sources = [(mode=(1,), port=1, current=0.00565e-6)]
        @test (@test_logs (:warn,) Circuit(netlist)) isa Circuit
        nl = @test_logs (:warn,) hbnlsolve(wp, (4,), sources, netlist; keyedarrays = false)
        @test isapprox(nl.nodeflux, hbnlsolve(wp, (4,), sources, typed; keyedarrays = false).nodeflux; rtol = 1e-10)
        hb = @test_logs (:warn,) hbsolve(ws, wp, sources, (2,), (4,), netlist; keyedarrays = false)
        @test isapprox(hb.linearized.S, hbsolve(ws, wp, sources, (2,), (4,), typed; keyedarrays = false).linearized.S; rtol = 1e-10)
        lin = @test_logs (:warn,) hblinsolve(ws, netlist; keyedarrays = false)
        @test isapprox(lin.S, hblinsolve(ws, typed; keyedarrays = false).S; rtol = 1e-10)
        @test samematrices(@test_logs((:warn,), numericmatrices(netlist, Dict())),
            numericmatrices(typed, Dict()))
        @test samematrices(@test_logs((:warn,), symbolicmatrices(netlist)),
            symbolicmatrices(typed))
        @test (@test_logs (:warn,) JosephsonCircuits.exportnetlist(netlist, Dict())).netlist isa String
        io = IOBuffer()
        @test_logs (:warn,) JosephsonCircuits.export_netlist!(io, netlist, Dict())
        back = Tuple{String,String,String,Any}[]
        @test_logs (:warn,) JosephsonCircuits.import_netlist!(io, back)
        @test length(back) == length(netlist)
    end

    @testset "a list of typed components is not a tuple netlist" begin
        # written without the Circuit around it, it is refused as what it
        # is, before any warning about the tuple netlist
        typedlist = [(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 0, Capacitor(1e-12))]
        @test_logs min_level = Base.CoreLogging.Warn @test_throws(
            ArgumentError, hblinsolve(2*pi*[5e9], typedlist))
    end

    @testset "the node order a tuple entry point is given" begin
        # nodes named rather than numbered compile in the order of their
        # names, which the entry points of the tuple format took as
        # `sorting`; the same netlist numbered is the reference
        named = [("P1","in","0",1), ("R1","in","0",50.0),
            ("C1","in","jj",100e-15), ("Lj1","jj","0",1000e-12),
            ("C2","jj","0",1000e-15)]
        numbered = [("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        wp = (2*pi*4.75001e9,); ws = 2*pi*[4.5e9, 5.0e9]
        sources = [(mode=(1,), port=1, current=0.00565e-6)]
        solve(n; kw...) = @deprecated hbsolve(ws, wp, sources, (2,), (4,), n;
            keyedarrays = false, kw...)
        @test isapprox(solve(named; sorting = :name).linearized.S,
            solve(numbered).linearized.S; rtol = 1e-10)
        @test isapprox((@deprecated hbnlsolve(wp, (4,), sources, named;
            sorting = :name, keyedarrays = false)).S,
            (@deprecated hbnlsolve(wp, (4,), sources, numbered;
            keyedarrays = false)).S; rtol = 1e-10)
        @test isapprox((@deprecated hblinsolve(ws, named; sorting = :name,
            keyedarrays = false)).S, (@deprecated hblinsolve(ws, numbered;
            keyedarrays = false)).S; rtol = 1e-10)
        @test samematrices(@deprecated(numericmatrices(named, Dict(); sorting = :name)),
            @deprecated(numericmatrices(numbered, Dict())))
        @test samematrices(@deprecated(symbolicmatrices(named; sorting = :none)),
            @deprecated(symbolicmatrices(numbered)))
    end

    @testset "resistor lookup retains values and branch orientation" begin
        # a port reads the resistor across its nodes whichever way the
        # resistor's terminals are written, and two ports across the same
        # nodes read the same resistor, by name and value
        c = @deprecated Circuit(Any[("P1", "a", "0", 1),
            ("P2", "0", "a", 2), ("R1", "0", "a", :R)])
        ports = [p for (_, p) in c.components if p isa Port]
        @test all(p -> p.Z0 === :R, ports)
        @test all(p -> p.termination.component == "R1", ports)
        @test isconcretetype(eltype(c.connections))

        # a resistor from a node to itself is one resistor across that pair
        loop = @deprecated Circuit([("P1", "0", "0", 1),
            ("R1", "0", "0", 50.0)])
        @test first(loop.components).second.Z0 == 50.0
        @deprecated @test_throws ArgumentError Circuit([("P1", "1", "0", 1)])
        @deprecated @test_throws ArgumentError Circuit([("P1", "1", "0", 1),
            ("R1", "0", "1", 50.0), ("R2", "1", "0", 75.0)])
    end

    @testset "the two netlist forms do not mix" begin
        # an entry ending in a value is a tuple netlist entry, which cannot
        # sit next to a typed one, and a tuple netlist has no interface
        @test_throws ArgumentError Circuit([
            (:p1, 1, 0, Port(1; Z0 = 50.0)), ("C1", "1", "0", 1e-12)])
        @deprecated @test_throws ArgumentError Circuit([("P1","1","0",1),
            ("R1","1","0",50.0)]; pins = [1 => ("P1", 1)])
        # the adapter's port owns the resistor the netlist already contains,
        # and no second termination is generated
        cl = compile(@deprecated Circuit([("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","0",1e-12)]))
        @test cl.componenttypes == [:P, :R, :C]
        @test cl.componentnames[only(cl.ports).environment] == "R1"
        @test length(cl.resistors) == 1
    end

    @testset "the tuple netlist file" begin
        path = tempdir()
        filename = joinpath(path, "JosephsonCircuits-" * string(JosephsonCircuits.UUIDs.uuid1()) * ".net")
        circuit1 = [("P","1","0",1),("R","1","0",50.0)]
        @deprecated JosephsonCircuits.export_netlist(filename, circuit1, Dict())
        circuit2 = @deprecated JosephsonCircuits.import_netlist(filename)
        rm(filename)
        @test isequal(circuit1, circuit2)
        # the importer parses component values into the package's own
        # expression type, so the round trip is checked against that
        R1, = JosephsonCircuits.@params R1
        circuit1 = [("P","1","0",1),("R","1","0",R1)]
        @deprecated JosephsonCircuits.export_netlist(filename, circuit1)
        circuit2 = @deprecated JosephsonCircuits.import_netlist(filename)
        rm(filename)
        @test isequal(circuit1, circuit2)
        # a line with the wrong number of fields
        io = IOBuffer("P 1 1")
        @deprecated @test_throws(
            ErrorException("each line should have component name, node1, node2, component value"),
            JosephsonCircuits.import_netlist!(io, Tuple{String,String,String,Any}[]))
    end

    @testset "the tuple netlist is exported as it was written" begin
        netlist = [("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",100e-15),
            ("Lj1","2","0",1e-9), ("C2","2","0",1e-12)]
        lines = split((@deprecated JosephsonCircuits.exportnetlist(netlist, Dict{Any,Any}())).netlist, "\n")
        @test any(startswith("R1 "), lines)
        @test any(startswith("C1 "), lines)
        @test any(startswith("B1 "), lines)
    end

    @testset "parsecomponenttype" begin
        @test_throws(
            ArgumentError("parsecomponenttype() currently only works for two letter components"),
            JosephsonCircuits.parsecomponenttype("BAD1",["Lj","BAD","L","C","K","I","R","P"])
        )

        @test_throws(
            ArgumentError("No component in allowedcomponents matches the name B1."),
            JosephsonCircuits.parsecomponenttype("B1",["Lj","L","C","K","I","R","P"])
        )
    end

    @testset "checkcomponenttypes" begin
        @test_throws(
            ArgumentError("Allowed components parsing check has failed for Lj. This can happen if a two letter long component comes after a one letter component. Please reorder allowedcomponents."),
            JosephsonCircuits.checkcomponenttypes(["L","Lj","C","K","I","R","P"])
        )
        # the order the tuple netlist reads its prefixes in
        @test JosephsonCircuits.checkcomponenttypes(
            JosephsonCircuits.legacyallowedcomponents)
    end

    @testset "tuple round trip" begin
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
        psc = compile((@deprecated Circuit(circuit)); sorting = :number)
        @test psc.componentnames == ["P1","I1","R1","C1","Lj1","C2"]
        @test psc.componenttypes == [:P, :I, :R, :C, :Lj, :C]
        @test all(isequal.(psc.componentvalues[2:end],
            [Ipump, Rleft, Cc, Lj, Cj]))
        # the port's value is its reference impedance, which a legacy
        # netlist states as the resistor across it
        @test isequal(psc.componentvalues[1], Rleft)

        # the value vector is narrowed from Vector{Any} to the tightest
        # element type the netlist admits. Every entry is a quantity, the
        # port slot carrying its impedance rather than the port number, so
        # that is the parameter type here
        @test eltype(psc.componentvalues) ===
            JosephsonCircuits.CircuitValues.Parameter

        # node "0" is ground and sorts first, so it is node index 1
        @test psc.nodenames == ["0","1","2"]
        @test psc.nodeindices == [2 2 2 2 3 3; 1 1 1 3 1 1]

        # these node names sort the same either way
        @test JosephsonCircuits.comparestruct(psc,
            compile((@deprecated Circuit(circuit)); sorting = :name))
    end

    @testset "round trip with mutual inductors" begin
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
        psc = compile((@deprecated Circuit(circuit)); sorting = :number)
        @test psc.componenttypes == [:P, :I, :R, :L, :K, :L, :C]
        @test psc.inductors == [4, 6]
        @test psc.mutualinductors == [5]

        # a mutual inductor names the two branches it couples instead of
        # attaching to nodes, so its column of the node index array is empty
        @test JosephsonCircuits.coupledinductornames(psc) == ["L1","L2"]
        @test psc.nodeindices[:,5] == [0, 0]
    end

    @testset "round trip with circuitdefs substitution" begin
        circuit = [
            ("P1","1","0",1),
            ("R1","1","0",:Rleft),
            ("L1","1","0",:lvalue),
            ("C1","1","0",:Cval),
        ]
        circuitdefs = Dict(:Rleft => 50.0, :lvalue => 2.0, :Cval => 1e-12)
        c = @deprecated Circuit(circuit, circuitdefs)
        psc = compile(c; sorting = :number)
        @test psc.componenttypes == [:P, :R, :L, :C]
        @test psc.componentvalues[2] == 50.0
        @test psc.componentvalues[3] == 2.0
        # without substitution the symbols pass through
        psc2 = compile((@deprecated Circuit(circuit)); sorting = :number)
        @test psc2.componentvalues[2] == :Rleft
    end

    @testset "hbsolve equivalence: tuples vs native" begin
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

        sol1 = @deprecated hbsolve(ws, wp, sources, (8,), (8,), circuit, circuitdefs)

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
        S3 = sol3.linearized.S((0,), 1, (0,), 1, :)
        @test isapprox(S1, S3; rtol = 1e-6)
        @test maximum(abs2.(S1)) > 1.0 # it is an amplifier
    end

    # A legacy netlist has no way to say which resistor is a port's
    # environment, so two across one port are refused: picking one silently
    # would reinterpret a netlist that is an error. The typed format has no
    # such restriction, because a port states its own termination.
    @testset "two resistors across a tuple netlist port are refused" begin
        @deprecated @test_throws ArgumentError Circuit([("P1","1","0",1),
            ("R1","1","0",50.0), ("R2","1","0",50.0), ("C1","1","0",1e-12)])

        # one is fine, and becomes the port's environment
        c = @deprecated Circuit([("P1","1","0",1), ("R1","1","0",50.0)])
        p = only([v for (k,v) in c.components if v isa Port])
        @test p.termination isa JosephsonCircuits.LegacyTermination
        @test p.termination.component == "R1"

        # and the message names the offenders
        err = @deprecated try
            Circuit([("P1","1","0",1), ("R1","1","0",50.0),
                     ("R2","1","0",25.0)]); nothing
        catch e; e; end
        @test occursin("R1", sprint(showerror, err))
        @test occursin("R2", sprint(showerror, err))
    end
end

@testset verbose=true "deprecated symbolic frequency variable" begin

    wp = (2*pi*4.75001*1e9,)
    ws = 2*pi*(4.5:0.1:5.0)*1e9
    sources = [(mode=(1,),port=1,current=0.00565e-6)]
    law(w) = 50.0*(1 + (w/1e11)^2)
    JosephsonCircuits.@params wsym
    mk(z) = Circuit([("P1", "1", "0", Port(1; Z0 = z)),
        ("C1", "1", "2", Capacitor(100.0e-15)),
        ("Lj1", "2", "0", JosephsonJunction(1000.0e-12)),
        ("C2", "2", "0", Capacitor(1000.0e-15))])

    @testset "a value in the frequency variable is the closure written out" begin
        # the deprecated spelling is rewritten into a FrequencyDependent
        # closure, so it must give the same scattering parameters
        Ssym = (@deprecated hbsolve(ws, wp, sources, (2,), (8,),
            mk(law(wsym)), Dict(); symfreqvar = wsym)).linearized.S
        Sfun = hbsolve(ws, wp, sources, (2,), (8,),
            mk(FrequencyDependent(law))).linearized.S
        @test isapprox(Array(Ssym), Array(Sfun), rtol = 1e-12)
    end

    @testset "a value mixing the variable with a closure resolves both" begin
        # the rewritten closure evaluates the frequency dependent leaves the
        # expression already carries as well as substituting the variable
        both(z) = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)),
            ("R1", "1", "2", Resistor(z)),
            ("P2", "2", "0", Port(2; Z0 = 50.0))])
        wl = [2e9, 3e9]
        Smix = (@deprecated hblinsolve(wl,
            both(50.0*(1 + 0.1*wsym/1e9) + FrequencyDependent(w -> 0.5*w/1e9)),
            Dict(); symfreqvar = wsym, keyedarrays = false)).S
        Sfun = hblinsolve(wl,
            both(FrequencyDependent(w -> 50.0*(1 + 0.1*w/1e9) + 0.5*w/1e9));
            keyedarrays = false).S
        @test isapprox(Array(Smix), Array(Sfun), rtol = 1e-12)
        # and a value named by a symbol whose definition is written in the
        # variable is the closure too
        Snamed = (@deprecated hblinsolve(wl, both(:Rfd),
            Dict(:Rfd => 50.0*(1 + 0.1*wsym/1e9)); symfreqvar = wsym,
            keyedarrays = false)).S
        Sfun = hblinsolve(wl, both(FrequencyDependent(w -> 50.0*(1 + 0.1*w/1e9)));
            keyedarrays = false).S
        @test isapprox(Array(Snamed), Array(Sfun), rtol = 1e-12)
    end

    @testset "every entry point which takes it warns" begin
        # hbsolve's warning is asserted above
        @deprecated JosephsonCircuits.hbnlsolve(wp, (2,), sources,
            mk(law(wsym)), Dict(); symfreqvar = wsym)
        @deprecated hblinsolve(ws, mk(law(wsym)), Dict(); symfreqvar = wsym)
    end

    @testset "the lossy value it resolves is a noise channel" begin
        # the channel is chosen on the resolved value, so a loss written
        # through the frequency variable is found as one written plainly
        C0, tand = 100.0e-15, 1e-3
        lossy(v) = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)),
            ("C1", "1", "2", Capacitor(v)),
            ("Lj1", "2", "0", JosephsonJunction(1000.0e-12)),
            ("C2", "2", "0", Capacitor(500.0e-15))])
        sym = @deprecated hbsolve(ws, wp, sources, (2,), (8,),
            lossy(C0*(1 - im*tand)*(1 + 0*wsym)), Dict();
            symfreqvar = wsym, keyedarrays = false)
        plain = hbsolve(ws, wp, sources, (2,), (8,),
            lossy(C0*(1 - im*tand)), Dict(); keyedarrays = false)
        @test [sym.linearized.componentnames[i]
            for i in sym.linearized.noiseportimpedanceindices] == ["C1"]
        @test isapprox(sym.linearized.QE, plain.linearized.QE, rtol = 1e-12)
    end

    @testset "a value which also needs a definition is named" begin
        # the rewrite leaves a value depending on something else alone, so
        # the undefined parameter is still reported with its component
        JosephsonCircuits.@params Rundef
        c = Circuit([("P1", "1", "0", Port(1; Z0 = Rundef/(1 + wsym*1e-12))),
            ("C1", "1", "0", Capacitor(100.0e-15)),
            ("L1", "1", "0", Inductor(1.0e-9))])
        err = @deprecated try
            JosephsonCircuits.hbnlsolve((2*pi*5e9,), (1,),
                [(mode=(1,),port=1,current=1e-8)], c, Dict();
                symfreqvar = wsym)
            nothing
        catch e; e; end
        @test err isa ArgumentError
        @test occursin("Rundef", sprint(showerror, err))
        @test occursin("P1", sprint(showerror, err))
    end

    @testset "a resolved value table cannot take one" begin
        # the keyword needs the definitions, which a table of already
        # resolved values does not carry
        psc = compile(mk(50.0))
        sf = JosephsonCircuits.truncfreqs(
            JosephsonCircuits.calcfreqsdft((2,)); dc = true, odd = false,
            even = true, maxintermodorder = Inf)
        vvn = JosephsonCircuits.componentvaluestonumber(psc.componentvalues,
            Dict{Any,Any}())
        @test_throws ArgumentError hblinsolve(ws, psc, vvn, sf;
            symfreqvar = wsym)
    end
end

# The deprecated forms warn and give the same numbers as the forms which
# replace them; the wording of the warnings is not pinned.
@testset verbose=true "deprecated" begin

    @testset "connectS is intraconnectS and interconnectS" begin
        Sa = rand(Complex{Float64},3,3)
        Sb = rand(Complex{Float64},3,3)
        Sout1 = zeros(Complex{Float64},1,1)
        Sout2 = zeros(Complex{Float64},4,4)
        a = @test_logs (:warn,) JosephsonCircuits.connectS(Sa,1,2)
        @test a == JosephsonCircuits.intraconnectS(Sa,1,2)
        b = @test_logs (:warn,) JosephsonCircuits.connectS(Sa,Sb,1,2)
        @test b == JosephsonCircuits.interconnectS(Sa,Sb,1,2)
        @test_logs (:warn,) JosephsonCircuits.connectS!(Sout1,Sa,1,2)
        @test Sout1 == JosephsonCircuits.intraconnectS(Sa,1,2)
        @test_logs (:warn,) JosephsonCircuits.connectS!(Sout2,Sa,Sb,1,2)
        @test Sout2 == JosephsonCircuits.interconnectS(Sa,Sb,1,2)
    end

    @testset "X_Y_to_sympletic is X_Y_to_symplectic" begin
        X, Y = JosephsonCircuits.rand_cptp_quadrature_pair(2)
        S = @test_logs (:warn,) JosephsonCircuits.X_Y_to_sympletic_pair(X, Y)
        @test S == JosephsonCircuits.X_Y_to_symplectic_pair(X, Y)
        S = @test_logs (:warn,) JosephsonCircuits.X_Y_to_sympletic_block(X, Y)
        @test S == JosephsonCircuits.X_Y_to_symplectic_block(X, Y)
    end

    @testset "the noise keyword of connectS_initialize" begin
        networks = [("A", rand(Complex{Float64}, 2, 2)),
            ("B", rand(Complex{Float64}, 2, 2))]
        connections = [[("A", 2), ("B", 1)]]
        init = @test_logs (:warn,) JosephsonCircuits.connectS_initialize(
            networks, connections; noise = true)
        @test init == JosephsonCircuits.connectS_initialize(networks,
            connections)
    end

    @testset "the deprecated solver forms and keywords" begin
        circuit = Circuit([(:P1, 1, 0, Port(1; Z0 = :Rleft)),
            (:C1, 1, 2, Capacitor(:Cc)), (:Lj1, 2, 0, JosephsonJunction(:Lj)),
            (:C2, 2, 0, Capacitor(:Cj))])
        circuitdefs = Dict{Symbol,Complex{Float64}}(
            :Lj =>1000.0e-12, :Cc => 100.0e-15, :Cj => 1000.0e-15,
            :Rleft => 50.0)
        ws = 2*pi*[4.5e9]
        wp = (2*pi*4.75001*1e9,)
        Ip = 0.00565e-6
        sources = [(mode=(1,),port=1,current=Ip)]
        ref = hbsolve(ws, wp, sources, (2,), (2,), circuit, circuitdefs;
            keyedarrays = false)
        same(sol) = isapprox(Array(sol.linearized.S), ref.linearized.S;
            rtol = 1e-10) && isapprox(Array(sol.nonlinear.nodeflux),
            ref.nonlinear.nodeflux; rtol = 1e-10)

        # the single pump frequency, integer harmonic count form of hbsolve:
        # its pump count Npumpmodes is the tuple form's (2*Npumpmodes,), so
        # the operating points agree; its signal mode set is the legacy
        # solver's own, so of the linearized outputs the signal to signal
        # entry is compared, with one signal mode either way
        old = @test_logs (:warn,) match_mode = :any hbsolve(ws, wp[1], Ip, 1, 2,
            circuit, circuitdefs, pumpports = [1], keyedarrays = true)
        oldref = hbsolve(ws, wp, sources, (1,), (4,), circuit, circuitdefs)
        @test isapprox(Array(old.nonlinear.nodeflux(outputmode = (1,))),
            Array(oldref.nonlinear.nodeflux(outputmode = (1,))); rtol = 1e-10)
        @test isapprox(old.linearized.S((0,), 1, (0,), 1, 1),
            oldref.linearized.S((0,), 1, (0,), 1, 1); rtol = 1e-6)
        # a compiled circuit keeps the order it was compiled with
        compiled = @test_logs (:warn,) hbsolve(ws, wp[1], Ip, 1, 2,
            compile(circuit), circuitdefs; pumpports = [1], keyedarrays = true)
        @test compiled.nonlinear.nodeflux == old.nonlinear.nodeflux
        # the line search keywords it took are ignored, with a warning of
        # their own besides the form's
        for kw in ((switchofflinesearchtol = 1,), (alphamin = 0.1,))
            sol = @test_logs (:warn,) (:warn,) hbsolve(ws, wp[1], Ip, 1, 2,
                circuit, circuitdefs; pumpports = [1], keyedarrays = true, kw...)
            @test sol.nonlinear.nodeflux == old.nonlinear.nodeflux
        end

        nlref = hbnlsolve(wp, (2,), sources, circuit, circuitdefs;
            keyedarrays = false)
        # whichever method solves
        for method in (NewtonKrylov(), Staged()), kw in
                ((switchofflinesearchtol = 1,), (alphamin = 0.1,),
                (maxharmonics = (2,),))
            sol = @test_logs (:warn,) match_mode = :any hbnlsolve(wp, (2,),
                sources, circuit, circuitdefs; keyedarrays = false,
                method = method, kw...)
            @test isapprox(sol.nodeflux, nlref.nodeflux;
                rtol = method isa Staged ? 1e-6 : 1e-10)
        end
        sol = @test_logs (:warn,) match_mode = :any hbsolve(ws, wp, sources,
            (2,), (2,), circuit, circuitdefs; keyedarrays = false,
            maxpumpharmonics = (2,))
        @test same(sol)
        # the line search keywords warn whichever method solves the pump
        for method in (NewtonKrylov(), Staged()), kw in
                ((switchofflinesearchtol = 1,), (alphamin = 0.1,))
            sol = @test_logs (:warn,) match_mode = :any hbsolve(ws, wp,
                sources, (2,), (2,), circuit, circuitdefs; keyedarrays = false,
                method = method, kw...)
            @test sol.nonlinear.solverinfo.converged
        end

        linref = hblinsolve(ws, circuit, circuitdefs; keyedarrays = false)
        for kw in ((returnZ = true,), (returnZadjoint = true,),
                (returnZsensitivity = true,), (returnZsensitivityadjoint = true,))
            lin = @test_logs (:warn,) match_mode = :any hblinsolve(ws, circuit,
                circuitdefs; keyedarrays = false, kw...)
            @test lin.S == linref.S
        end
    end

    @testset "ftol is atol" begin
        circuit = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 0, Capacitor(1e-12)), (:Lj1, 1, 0, JosephsonJunction(1e-9))])
        wp = (2pi*5e9,)
        src = [(mode = (1,), port = 1, current = 1e-8)]
        old = @test_logs (:warn,) hbnlsolve(wp, (2,), src, circuit; ftol = 1e-10, keyedarrays = false)
        new = hbnlsolve(wp, (2,), src, circuit; atol = 1e-10, keyedarrays = false)
        @test old.nodeflux == new.nodeflux
        ws = 2pi*[4e9]
        olds = @test_logs (:warn,) hbsolve(ws, wp, src, (1,), (2,), circuit; ftol = 1e-10, keyedarrays = false)
        news = hbsolve(ws, wp, src, (1,), (2,), circuit; atol = 1e-10, keyedarrays = false)
        @test olds.linearized.S == news.linearized.S
    end

    @testset "solveS! without the fill reducing ordering" begin
        # the argument list v0.5.4 took, without the ordering
        # solveS_initialize now returns last
        networks = [("S1", rand(ComplexF64, 4, 4, 3)),
            ("S2", rand(ComplexF64, 3, 3, 3))]
        connections = [[("S1", 1), ("S1", 2), ("S1", 3)], [("S1", 4), ("S2", 2)]]
        init = JosephsonCircuits.solveS_initialize(networks, connections)
        new = copy(JosephsonCircuits.solveS!(init...)[1])
        old = @test_logs (:warn,) JosephsonCircuits.solveS!(init[1:end-1]...)
        @test old[1] ≈ new
    end

end
