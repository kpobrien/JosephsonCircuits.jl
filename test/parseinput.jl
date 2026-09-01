using Symbolics
using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "input path" begin

    @testset "hierarchy: repeated subcircuit vs flat legacy" begin
        function laddercell(L, Lj, C)
            return Circuit(
                [:l => Inductor(L), :jj => JosephsonJunction(Lj),
                 :c => Capacitor(C)],
                [((:l, 2), (:jj, 1), (:c, 1)),
                 ((:jj, 2), Ground), ((:c, 2), Ground)],
                Interface(pins = [1 => (:l, 1), 2 => (:l, 2)]),
            )
        end
        L, Lj, C = 100e-12, 400e-12, 50e-15
        cell = laddercell(L, Lj, C)
        Ncells = 3
        components = Any[Symbol("cell", i) => cell for i in 1:Ncells]
        push!(components, :p1 => Port(1; termination = nothing))
        push!(components, :r1 => Resistor(50.0))
        push!(components, :p2 => Port(2; termination = nothing))
        push!(components, :r2 => Resistor(50.0))
        connections = Any[
            ((Symbol("cell", i), 2), (Symbol("cell", i+1), 1))
            for i in 1:Ncells-1]
        push!(connections, ((:p1, 1), (:r1, 1), (:cell1, 1)))
        push!(connections, ((:p2, 1), (:r2, 1), (Symbol("cell", Ncells), 2)))
        push!(connections, ((:p1, 2), (:r1, 2), Ground))
        push!(connections, ((:p2, 2), (:r2, 2), Ground))
        hier = Circuit(components, connections)

        elab = elaborate(hier)
        # shared definitions are deduplicated by identity: one Inductor, one
        # NonlinearInductor, one Capacitor from the shared cell, plus two
        # ports and one shared-value resistor definition... Resistor(50.0)
        # was constructed twice, so it appears twice; the cell internals
        # appear once.
        @test count(d -> d isa Inductor, elab.definitions) == 1
        @test count(d -> d isa NonlinearInductor, elab.definitions) == 1
        @test count(d -> d isa Capacitor, elab.definitions) == 1
        @test JosephsonCircuits.ninstances(elab) == 3*Ncells + 4
        @test "cell1/l" in elab.instancepaths
        @test "cell3/jj" in elab.instancepaths
        # fresh internal nets per instance: each cell has one internal net
        # (between l, jj, c), and they are distinct
        @test length(unique(elab.terminalnets)) == JosephsonCircuits.nnets(elab)

        # flat legacy equivalent
        flat = [
            ("P1","1","0",1),
            ("R1","1","0",50.0),
            ("L1","1","2",L), ("Lj1","2","0",Lj), ("C1","2","0",C),
            ("L2","2","3",L), ("Lj2","3","0",Lj), ("C2","3","0",C),
            ("L3","3","4",L), ("Lj3","4","0",Lj), ("C3","4","0",C),
            ("P2","4","0",1),
            ("R2","4","0",50.0),
        ]
        # the port numbers must match: legacy P2 has value 1 -> fix to 2
        flat[end-1] = ("P2","4","0",2)
        ws = 2*pi*(3.0:0.5:7.0)*1e9
        sol1 = hblinsolve(ws, flat, Dict{Symbol,Float64}())
        sol2 = hblinsolve(ws, hier)
        S21flat = sol1.S((0,), 2, (0,), 1, :)
        S21hier = sol2.S((0,), 2, (0,), 1, :)
        @test isapprox(S21flat, S21hier; rtol = 1e-6)
    end

    @testset "pair syntax and subcircuit ports" begin
        function resonator()
            return Circuit(
                [:l => Inductor(1e-9), :c => Capacitor(1e-12)],
                [((:l, 1), (:c, 1)), ((:l, 2), (:c, 2))],
                Interface(pins = [1 => (:l, 1), 2 => (:l, 2)],
                    ports = [1 => (1, 2)]),
            )
        end
        r1 = resonator()
        # bond two resonator ports with pair syntax; key 1 is both a pin
        # and a port, so PortRef selects the port namespace
        c = Circuit(
            [:a => r1, :b => r1, :p1 => Port(1; termination = nothing), :rt => Resistor(50.0)],
            [PortRef(:a, 1) => PortRef(:b, 1),
             ((:p1, 1), (:rt, 1), (:a, 1)),
             ((:p1, 2), (:rt, 2), (:a, 2), Ground)],
        )
        elab = elaborate(c)
        # port bond merges a/l terminal1 with b/l terminal1, and a's pin 2
        # net with b's pin 2 net
        @test JosephsonCircuits.ninstances(elab) == 6
        psc = compile(elab)
        # ground plus one signal net: the bond merges a's pin 2 net with
        # b's pin 2 net, and the top level ties that net to Ground
        @test psc.Nnodes == 2

        # a pin/port key collision requires PortRef or PinRef
        c2 = Circuit([:a => r1, :b => r1], Any[]; validate = false)
        @test_throws ArgumentError JosephsonCircuits.parsecircuitlevel(
            [:a => r1, :b => r1], [(:a, 1) => (:b, 1)], nothing)
        # explicit PortRef resolves it
        pl = JosephsonCircuits.parsecircuitlevel([:a => r1, :b => r1],
            [PortRef(:a, 1) => PortRef(:b, 1)], nothing)
        @test length(pl.groups) == 2
        # explicit PinRef bonds single pins
        pl2 = JosephsonCircuits.parsecircuitlevel([:a => r1, :b => r1],
            [PinRef(:a, 1) => PinRef(:b, 1)], nothing)
        @test length(pl2.groups) == 1
    end

    @testset "validation diagnostics" begin
        # duplicate identifiers
        @test_throws ArgumentError Circuit(
            [:x => Inductor(1e-9), :x => Capacitor(1e-12)], Any[])
        # dangling endpoint
        @test_throws ArgumentError Circuit([:x => Inductor(1e-9)],
            [((:x, 1), (:y, 1))])
        # invalid terminal
        @test_throws ArgumentError Circuit([:x => Inductor(1e-9)],
            [((:x, 3), Ground)])
        # mutual inductor in connections
        @test_throws ArgumentError Circuit(
            [:l1 => Inductor(1e-9), :l2 => Inductor(1e-9),
             :k => MutualInductor(0.5, :l1, :l2)],
            [((:k, 1), Ground)])
        # mutual inductor partner missing
        @test_throws ArgumentError Circuit(
            [:l1 => Inductor(1e-9), :k => MutualInductor(0.5, :l1, :l2)],
            Any[])
        # mutual inductor partner is not an inductor
        @test_throws ArgumentError Circuit(
            [:l1 => Inductor(1e-9), :c1 => Capacitor(1e-12),
             :k => MutualInductor(0.5, :l1, :c1)],
            Any[])
        # subcircuit without an interface
        inner = Circuit([:l => Inductor(1e-9)], [((:l, 1), Ground)])
        @test_throws ArgumentError Circuit([:sub => inner], Any[])
        # identifiers with the path separator
        @test_throws ArgumentError Circuit(["a/b" => Inductor(1e-9)], Any[])
        # Instance passthrough and override rejection
        c = Circuit([:l => Instance(Inductor(1e-9))], [((:l, 1), Ground)])
        @test JosephsonCircuits.instancedefinition(elaborate(c), 1) isa Inductor
        @test_throws ArgumentError Instance(Inductor(1e-9);
            thermal_bindings = (:body => :t1,))
        # recursion
        comps = Pair{Symbol,Any}[]
        rec = Circuit(comps, Any[],
            Interface(pins = [1 => (:self, 1)]); validate = false)
        push!(comps, :self => rec)
        @test_throws ArgumentError elaborate(rec)
        # net name collision
        @test_throws ArgumentError elaborate(Circuit(
            [:l1 => Inductor(1e-9), :l2 => Inductor(1e-9)],
            [Net("a", ((:l1, 1),)), Net("a", ((:l2, 1),)),
             ((:l1, 2), Ground), ((:l2, 2), Ground)]))
        # names attached to the ground class are ignored
        elab = elaborate(Circuit([:l1 => Inductor(1e-9)],
            [Net("mygnd", ((:l1, 2), Ground)), Net("sig", ((:l1, 1),))]))
        @test elab.netnames[1] == "0"
        @test "sig" in elab.netnames
        @test !("mygnd" in elab.netnames)
        # a single endpoint written where a group is expected gets a
        # targeted hint, in both raw groups and Net entries
        err = try
            Circuit([:l1 => Inductor(1e-9)], [(:l1, 1)])
        catch e
            e
        end
        @test err isa ArgumentError
        @test occursin("one-endpoint group", err.msg)
        @test_throws ArgumentError Circuit([:l1 => Inductor(1e-9)],
            [Net("a", (:l1, 1))])
        # compact display
        cshow = Circuit([:l1 => Inductor(1e-9)],
            [((:l1, 1), (:l1, 2), Ground)])
        @test repr(cshow) == "Circuit(1 components, 1 connections)"
        @test repr(Ground) == "Ground"
        rshow = Circuit([:l => Inductor(1e-9), :c => Capacitor(1e-12)],
            [((:l, 1), (:c, 1)), ((:l, 2), (:c, 2))],
            Interface(pins = [1 => (:l, 1), 2 => (:l, 2)],
                ports = [1 => (1, 2)]))
        @test repr(rshow) == "Circuit(2 components, 2 connections, 2 pins, 1 ports)"
    end

    @testset "ground as a component" begin
        # `Ground()` is the same singleton as `Ground`, so the sentinel and
        # component spellings cannot diverge
        @test Ground() === Ground

        # a declared ground instance referenced through its terminal is
        # equivalent to the bare sentinel, down to identical parse output
        mk(g) = Circuit(
            Any[:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                :cc => Capacitor(100e-15), :jj => JosephsonJunction(1e-9),
                :cj => Capacitor(100e-15),
                (g ? [:gnd => Ground()] : [])...],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:cj, 1)),
             g ? ((:p1, 2), (:r1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)) :
                 ((:p1, 2), (:r1, 2), (:jj, 2), (:cj, 2), Ground)])
        @test JosephsonCircuits.comparestruct(
            compile(mk(false)), compile(mk(true)))

        # the ground instance flattens to no component
        elab = JosephsonCircuits.elaborate(mk(true))
        @test JosephsonCircuits.ninstances(elab) == 5
        @test !any(d -> d isa JosephsonCircuits.GroundType, elab.definitions)

        # pair syntax, PinRef, and both definition spellings
        cpair = Circuit(
            [:r1 => Resistor(50.0), :p1 => Port(1; termination = nothing),
             :g1 => Ground, :g2 => Ground()],
            [((:p1, 1), (:r1, 1)), (:p1, 2) => (:g1, 1),
             (JosephsonCircuits.PinRef(:r1, 2),
              JosephsonCircuits.PinRef(:g2, 1))])
        @test JosephsonCircuits.nnets(JosephsonCircuits.elaborate(cpair)) == 2

        # a subcircuit may ground itself through its own ground instance,
        # and every instance shares the one global reference net
        sub = Circuit(
            [:l => Inductor(1e-9), :c => Capacitor(1e-13), :gnd => Ground()],
            [((:l, 2), (:c, 1)), ((:c, 2), (:gnd, 1))],
            Interface(pins = [1 => (:l, 1)]))
        top = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :s1 => sub, :s2 => sub,
             :gnd => Ground()],
            [((:p1, 1), (:r1, 1), (:s1, 1), (:s2, 1)),
             ((:p1, 2), (:r1, 2), (:gnd, 1))])
        etop = JosephsonCircuits.elaborate(top)
        @test JosephsonCircuits.ninstances(etop) == 6
        capnets = [JosephsonCircuits.instanceterminals(etop, i)[2]
                   for i in 1:JosephsonCircuits.ninstances(etop)
                   if JosephsonCircuits.instancedefinition(etop, i) isa Capacitor]
        @test all(==(1), capnets)

        # a ground instance has exactly one terminal and no ports
        @test_throws ArgumentError Circuit(
            [:r => Resistor(50.0), :gnd => Ground()],
            [((:r, 1), (:gnd, 2)), ((:r, 2),)])
        @test_throws ArgumentError Circuit(
            [:r => Resistor(50.0), :gnd => Ground()],
            [((:r, 1), (:gnd, 1, 1)), ((:r, 2),)])
        # an interface pin may not map to ground through an instance either
        @test_throws ArgumentError Circuit(
            [:l => Inductor(1e-9), :gnd => Ground()],
            [((:l, 2), (:gnd, 1))],
            Interface(pins = [1 => (:l, 1), 2 => (:gnd, 1)]))
        # a mutual inductor may not couple a ground instance
        @test_throws ArgumentError Circuit(
            [:l1 => Inductor(1e-9), :gnd => Ground(),
             :k => MutualInductor(0.5, :l1, :gnd)],
            [((:l1, 1),), ((:l1, 2), (:gnd, 1))])
    end

    @testset "interface through pins/ports keywords" begin
        comps = [:jj => JosephsonJunction(1e-9), :cj => Capacitor(40e-15),
                 :cg => Capacitor(76.6e-15)]
        conns = [((:jj, 1), (:cj, 1), (:cg, 1)), ((:jj, 2), (:cj, 2)),
                 ((:cg, 2), Ground)]
        thepins = [1 => (:jj, 1), 2 => (:jj, 2)]
        cpos = Circuit(comps, conns, Interface(pins = thepins))
        ckw = Circuit(comps, conns; pins = thepins)
        @test ckw.interface isa Interface
        @test isequal(cpos.interface.pins, ckw.interface.pins)

        # identical as a subcircuit
        mkt(sub) = Circuit(
            [:p => Port(1; termination = nothing), :r => Resistor(50.0), :s => sub],
            [((:p, 1), (:r, 1), (:s, 1)), ((:p, 2), (:r, 2), (:s, 2), Ground)])
        e1 = JosephsonCircuits.elaborate(mkt(cpos))
        e2 = JosephsonCircuits.elaborate(mkt(ckw))
        @test e1.terminalnets == e2.terminalnets
        @test e1.netnames == e2.netnames

        # ports keyword
        cports = Circuit([:l => Inductor(1e-9)], [((:l, 1),), ((:l, 2),)];
            pins = [1 => (:l, 1), 2 => (:l, 2)], ports = [1 => (1, 2)])
        @test JosephsonCircuits.hasports(cports)

        # positional interface and keywords are mutually exclusive
        @test_throws ArgumentError Circuit(comps, conns,
            Interface(pins = thepins); pins = thepins)
        # ports require pins, as in Interface
        @test_throws ArgumentError Circuit(
            [:l => Inductor(1e-9)], [((:l, 1),), ((:l, 2),)];
            ports = [1 => (1, 2)])
    end

    @testset "port owns its termination" begin
        # Keyword arguments do not participate in dispatch, so there is one
        # Port constructor and not several. Separate methods for the
        # terminated and unterminated spellings would redefine one another
        # and the survivor would silently answer for both.
        @test Port(1).termination isa JosephsonCircuits.MatchedTermination
        @test Port(1; termination = nothing).termination isa
            JosephsonCircuits.NoPortTermination
        @test length(methods(Port, Tuple{Int})) == 1
        @test repr(Port(1)) == "Port(1; Z0 = 50.0)"
        @test repr(Port(2; Z0 = 1000.0, termination = nothing)) ==
            "Port(2; Z0 = 1000.0, termination = nothing)"

        # a reference impedance must be finite, real and positive once numeric
        @test_throws ArgumentError Port(1; Z0 = 0.0)
        @test_throws ArgumentError Port(1; Z0 = -50.0)
        @test_throws ArgumentError Port(1; Z0 = Inf)
        @test_throws ArgumentError Port(1; Z0 = 50.0 + 1.0im)
        @test_throws ArgumentError Port(1; termination = :matched)

        # a matched port and the explicit resistor it replaces are the same
        # circuit, electrically and numerically
        mk(matched) = Circuit(
            matched ?
              [:p1 => Port(1), :cc => Capacitor(100e-15),
               :jj => JosephsonJunction(1e-9), :cj => Capacitor(1e-12)] :
              [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
               :cc => Capacitor(100e-15), :jj => JosephsonJunction(1e-9),
               :cj => Capacitor(1e-12)],
            matched ?
              [[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
               [(:p1,2),(:jj,2),(:cj,2),Ground]] :
              [[(:p1,1),(:r1,1),(:cc,1)], [(:cc,2),(:jj,1),(:cj,1)],
               [(:p1,2),(:r1,2),(:jj,2),(:cj,2),Ground]])
        nmm, nme = map((true, false)) do matched
            psc = compile(mk(matched))
            numericmatrices(psc, calccircuitgraph(psc), Dict{Any,Any}();
                Nmodes = 1)
        end
        @test nmm.Gnm == nme.Gnm
        @test nmm.Cnm == nme.Cnm
        @test nmm.Lmean == nme.Lmean
        ws = 2*pi*(4.5:0.05:5.0)*1e9
        sol(c) = hbsolve(ws, (2*pi*4.75001e9,),
            [(mode=(1,),port=1,current=0.00565e-6)], (8,), (16,), c)
        @test Array(sol(mk(true)).linearized.S) ==
            Array(sol(mk(false)).linearized.S)

        # Z0 reaches the solver: it is the environment, not a discarded label
        pscz = compile(Circuit(
            [:p1 => Port(1), :p2 => Port(2; Z0 = 1000.0),
             :l => Inductor(1e-9), :jj => JosephsonJunction(1e-9)],
            [[(:p1,1),(:jj,1),(:l,1)], [(:l,2),(:p2,1)],
             [(:p1,2),(:p2,2),(:jj,2),Ground]]))
        rvals = [pscz.componentvalues[i]
                 for (i,t) in enumerate(pscz.componenttypes) if t == :R]
        @test rvals == [50.0, 1000.0]
        @test all(endswith("/termination"),
            [pscz.componentnames[i]
             for (i,t) in enumerate(pscz.componenttypes) if t == :R])

        # a differential port terminates between its own terminals and does
        # not reach for ground
        pscd = compile(Circuit(
            [:p1 => Port(1), :ca => Capacitor(1e-12), :cb => Capacitor(1e-12),
             :l => Inductor(1e-9)],
            [[(:p1,1),(:ca,1),(:l,1)], [(:p1,2),(:cb,1),(:l,2)],
             [(:ca,2),(:cb,2),Ground]]))
        ip = findfirst(==(:P), pscd.componenttypes)
        ir = findfirst(==(:R), pscd.componenttypes)
        @test pscd.nodeindices[:,ir] == pscd.nodeindices[:,ip]
        @test !any(isone, pscd.nodeindices[:,ir])   # node 1 is ground

        # the legacy netlist states its port impedance as a resistor, which
        # the adapter discovers once and records; it must not gain a second
        # termination
        legacy = [("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",100e-15),
                  ("Lj1","2","0",1e-9), ("C2","2","0",1e-12)]
        lc = Circuit(legacy)
        lport = only([v for (k,v) in lc.components if v isa Port])
        @test lport.Z0 == 50.0
        # the netlist's own resistor is the port's environment, named rather
        # than rediscovered later
        @test lport.termination isa JosephsonCircuits.LegacyTermination
        @test count(==(:R), compile(lc).componenttypes) == 1
        @test Array(sol(lc).linearized.S) == Array(sol(mk(true)).linearized.S)
    end

    @testset "compiled circuit" begin
        # The compiler is the only front end, so the parity evidence for it
        # is the rest of this suite continuing to produce the numbers it did
        # against the parsed representation. What is asserted here is the
        # structure that representation could not express.
        JC = JosephsonCircuits
        blk = ScatteringParameters([0.0 1.0; 1.0 0.0])
        cell = Circuit(
            [:jj => JosephsonJunction(1e-9), :cj => Capacitor(4e-14),
             :cg => Capacitor(7e-14), :gnd => Ground()],
            [[(:jj,1),(:cj,1),(:cg,1)], [(:jj,2),(:cj,2)], [(:cg,2),(:gnd,1)]];
            pins = [1 => (:jj,1), 2 => (:jj,2)])
        c = Circuit(
            Any[:p1 => Port(1), :p2 => Port(2; Z0 = 1000.0),
                :cc => Capacitor(1e-13), :s1 => cell, :s2 => cell,
                :l1 => Inductor(1e-9), :l2 => Inductor(2e-9),
                :k => MutualInductor(0.5, :l1, :l2), :blk => blk,
                :i1 => CurrentSource(1e-9), :gnd => Ground()],
            Any[[(:p1,1),(:cc,1)], [(:cc,2),(:s1,1)], [(:s1,2),(:s2,1)],
                [(:s2,2),(:l1,1)], [(:l1,2),(:blk,1)], [(:blk,2),(:l2,1)],
                [(:l2,2),(:p2,1),(:i1,1)],
                [(:p1,2),(:p2,2),(:i1,2),(:gnd,1)]])
        cc = JC.compile(c)

        # the groups partition the flat table exactly once
        grouped = sort(vcat(cc.capacitors, cc.resistors, cc.inductors,
            cc.junctions, cc.nonlinearinductors, cc.currentsources,
            cc.mutualinductors, [p.component for p in cc.ports]))
        @test grouped == collect(1:JC.ncomponents(cc))

        # and each group holds only its own type
        for (g, t) in ((cc.capacitors, :C), (cc.resistors, :R),
                (cc.inductors, :L), (cc.junctions, :Lj),
                (cc.currentsources, :I), (cc.mutualinductors, :K))
            @test all(i -> cc.componenttypes[i] === t, g)
        end

        # a port carries its reference impedance and a direct index to the
        # environment it owns, rather than being matched to a resistor by
        # sharing a branch with it
        @test [p.number for p in cc.ports] == [1, 2]
        @test [p.zref for p in cc.ports] == [50.0, 1000.0]
        for p in cc.ports
            @test p.environment != 0
            @test cc.componenttypes[p.environment] === :R
            @test cc.componentvalues[p.environment] == p.zref
            @test cc.nodeindices[:, p.environment] ==
                [p.positivenode, p.negativenode]
            @test cc.componentnames[p.environment] ==
                cc.componentnames[p.component]*"/termination"
        end

        # a port which declares no termination of its own owns none, whatever
        # the user placed across it. The resistors stay device resistors,
        # any number of them may share the terminals, and the reference
        # impedance is the one the port declared
        cu = JC.compile(Circuit(
            [:p1 => Port(1; Z0 = 50.0, termination = nothing),
             :r1 => Resistor(100.0), :r2 => Resistor(200.0),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:r1,1),(:r2,1),(:c1,1)],
             [(:p1,2),(:r1,2),(:r2,2),(:c1,2),Ground]]))
        @test only(cu.ports).environment == 0
        @test only(cu.ports).zref == 50.0
        nmu = JC.numericmatrices(cu, calccircuitgraph(cu), Dict{Any,Any}();
            Nmodes = 1)
        # normalized to the declared Z0, not to the resistor beside it
        @test nmu.portimpedances == [50.0]
        @test nmu.portenvironmentindices == [0]
        # and both resistors are still internal noise channels
        @test [cu.componentnames[i] for i in nmu.noiseportimpedanceindices] ==
            ["r1", "r2"]
        # A port owning a matched environment beside a device resistor of
        # the same value is loaded twice. That is legal and occasionally
        # meant, but it is much more often a circuit written when a resistor
        # was how a port got its impedance, so it is reported rather than
        # refused or silently accepted.
        dup = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :r1 => Resistor(50.0),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:r1,1),(:c1,1)], [(:p1,2),(:r1,2),(:c1,2),Ground]])
        @test_logs (:warn,) JC.compile(dup)
        # a device resistor of a different value is a device resistor
        notdup = Circuit(
            [:p1 => Port(1; Z0 = 50.0), :r1 => Resistor(100.0),
             :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:r1,1),(:c1,1)], [(:p1,2),(:r1,2),(:c1,2),Ground]])
        @test_logs JC.compile(notdup)

        # with nothing across it the port is simply unloaded
        cnone = JC.compile(Circuit(
            [:p1 => Port(1; termination = nothing), :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2),Ground]]))
        @test only(cnone.ports).environment == 0

        # and the waves are normalized to the declared reference impedance
        # rather than to whatever resistor shares the branch, which is what
        # the geometric rule could not express. A port which owns its
        # matched environment reflects (R - Z0)/(R + Z0) off a load R; an
        # unterminated one is an ideal current source into the same load and
        # gives (2R - Z0)/Z0, and both follow Z0 while the circuit stands
        # still
        onenode(port, R) = hblinsolve([2*pi*5e9], Circuit(
            [:p1 => port, :r => Resistor(R)],
            [[(:p1,1),(:r,1)], [(:p1,2),(:r,2),Ground]]);
            keyedarrays = false).S[1,1,1,1,1]
        @test isapprox(onenode(Port(1; Z0 = 50.0), 100.0), 1/3; atol = 1e-12)
        @test isapprox(onenode(Port(1; Z0 = 100.0), 100.0), 0; atol = 1e-12)
        @test isapprox(onenode(Port(1; Z0 = 50.0, termination = nothing),
            100.0), 3.0; atol = 1e-12)
        @test isapprox(onenode(Port(1; Z0 = 75.0, termination = nothing),
            100.0), (2*100.0 - 75.0)/75.0; atol = 1e-12)

        # the block keeps its whole terminal map, and it is the only place
        # the block appears: it has no entries in the flat table
        b = only(cc.scatteringblocks)
        @test b.definition === blk
        @test length(b.signalnodes) == blk.nports
        @test !any(==(:S), cc.componenttypes)
        for p in 1:blk.nports
            @test b.signalnodes[p] > 0
            @test b.refnodes[p] > 0
        end
        @test b.path == "blk"

        # a legacy netlist compiles too, and its port owns nothing because
        # the netlist already states its termination
        cl = JC.compile(Circuit([("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","0",1e-12)]))
        # a legacy port owns the resistor the netlist already contains, and
        # no second one is generated
        @test cl.componentnames[only(cl.ports).environment] == "R1"
        @test length(cl.resistors) == 1
    end

    @testset "port and noise roles from the compiled circuit" begin
        JC = JosephsonCircuits
        # a legacy netlist names the resistor it already contains as the
        # port's environment, so nothing downstream searches for it
        lc = Circuit([("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","2",100e-15), ("Lj1","2","0",1e-9), ("C2","2","0",1e-12)])
        lport = only([v for (k,v) in lc.components if v isa Port])
        @test lport.termination isa JC.LegacyTermination
        @test lport.termination.component == "R1"
        lcc = JC.compile(lc)
        @test lcc.componentnames[only(lcc.ports).environment] == "R1"
        # and no second termination is generated
        @test count(==(:R), lcc.componenttypes) == 1

        # the three lists a circuit's roles produce, by name so that the
        # assertion does not restate the component ordering
        function roles(c)
            cc = JC.compile(c)
            vals = JC.componentvaluestonumber(cc.componentvalues,
                Dict{Any,Any}())
            idx, numbers = JC.portindicesnumbers(cc)
            nm(v) = [cc.componentnames[i] for i in v]
            return (ports = nm(idx), numbers = numbers,
                environments = nm(JC.portenvironmentindices(cc)),
                noise = nm(JC.noiseindices(cc, vals)))
        end

        # the port's environment is the resistor the netlist already had, and
        # a circuit whose only resistor is that environment has no internal
        # noise channel at all
        @test roles(lc) == (ports = ["P1"], numbers = [1],
            environments = ["R1"], noise = String[])

        # a second resistor is an ordinary device resistor and a capacitor
        # with a nonzero imaginary part is lossy, so both are noise channels
        @test roles(Circuit([("P1","1","0",1), ("R1","1","0",50.0),
            ("R2","1","2",75.0), ("C1","1","2",100e-15 + 1e-18im),
            ("Lj1","2","0",1e-9)])) ==
            (ports = ["P1"], numbers = [1], environments = ["R1"],
             noise = ["R2", "C1"])

        # two typed ports each generate and own their own termination
        @test roles(Circuit(
            Any[:p1 => Port(1), :p2 => Port(2; Z0 = 1000.0),
                :cc => Capacitor(1e-13), :jj => JosephsonJunction(1e-9),
                :gnd => Ground()],
            Any[[(:p1,1),(:cc,1)], [(:cc,2),(:jj,1),(:p2,1)],
                [(:p1,2),(:p2,2),(:jj,2),(:gnd,1)]])) ==
            (ports = ["p1", "p2"], numbers = [1, 2],
             environments = ["p1/termination", "p2/termination"],
             noise = String[])

        # the restriction the geometric rule imposed is gone: a port may
        # share its terminals with an ordinary device resistor, which loads
        # it in parallel and remains an internal noise channel
        cshared = Circuit(
            Any[:p1 => Port(1), :rload => Resistor(100.0),
                :cc => Capacitor(1e-13), :jj => JosephsonJunction(1e-9),
                :gnd => Ground()],
            Any[[(:p1,1),(:rload,1),(:cc,1)], [(:cc,2),(:jj,1)],
                [(:p1,2),(:rload,2),(:jj,2),(:gnd,1)]])
        cc = JC.compile(cshared); b = JC.bind(cc)
        cg = calccircuitgraph(cc)
        nm = JC.assemblematrices(JC.circuitmatrixplan(cc, cg, b; Nmodes = 1), b)
        # the geometric rule refused this circuit, because it could only read
        # a resistor across a port's terminals as that port's environment;
        # now both assemblies take the roles from the ports and agree
        ref = numericmatrices(cc, cg, Dict{Any,Any}(); Nmodes = 1)
        @test nm.portimpedances == ref.portimpedances
        @test nm.portenvironmentindices == ref.portenvironmentindices
        @test nm.noiseportimpedanceindices == ref.noiseportimpedanceindices
        @test cc.componentnames[only(nm.portenvironmentindices)] ==
            "p1/termination"
        @test [cc.componentnames[i] for i in nm.noiseportimpedanceindices] ==
            ["rload"]
        @test nm.Gnm[1,1] == 1/50 + 1/100

        # the checks the geometric path made are still made
        @test_throws ArgumentError JC.portindicesnumbers(JC.compile(Circuit(
            Any[:p1 => Port(1), :p2 => Port(1), :c1 => Capacitor(1e-12),
                :gnd => Ground()],
            Any[[(:p1,1),(:p2,1),(:c1,1)],
                [(:p1,2),(:p2,2),(:c1,2),(:gnd,1)]])))
        @test_throws ArgumentError JC.portindicesnumbers(JC.compile(Circuit(
            Any[:p1 => Port(1), :p2 => Port(2), :c1 => Capacitor(1e-12),
                :gnd => Ground()],
            Any[[(:p1,1),(:p2,1),(:c1,1)],
                [(:p1,2),(:p2,2),(:c1,2),(:gnd,1)]])))
        # a port which owns nothing reports zero rather than failing: it has
        # a reference impedance like every port, it simply loads nothing
        @test JC.portenvironmentindices(JC.compile(
            Circuit([:p1 => Port(1; termination = nothing),
                     :c1 => Capacitor(1e-12)],
                [[(:p1,1),(:c1,1)], [(:p1,2),(:c1,2),Ground]]))) == [0]

        # the per group element types are what the whole table scan produced
        function typesagree(c)
            cc = JC.compile(c)
            vals = JC.componentvaluestonumber(cc.componentvalues,
                Dict{Any,Any}())
            return all(((g, t),) -> JC.grouptype(vals, g, true) ==
                    eltype(JC.calcvaluetype(cc.componenttypes, vals, [t])),
                ((cc.capacitors, :C), (cc.resistors, :R),
                 (cc.inductors, :L), (cc.junctions, :Lj),
                 (cc.nonlinearinductors, :NL), (cc.currentsources, :I),
                 (cc.mutualinductors, :K)))
        end
        @test typesagree(lc)
        @test typesagree(Circuit(
            Any[:p1 => Port(1), :l1 => Inductor(1e-9),
                :c1 => Capacitor(1e-12 + 1im*1e-15), :c2 => Capacitor(2e-12),
                :jj => JosephsonJunction(1e-9), :i1 => CurrentSource(1),
                :l2 => Inductor(2e-9), :k => MutualInductor(0.5, :l1, :l2),
                :gnd => Ground()],
            Any[[(:p1,1),(:l1,1),(:c1,1),(:i1,1)],
                [(:l1,2),(:c2,1),(:jj,1),(:l2,1)],
                [(:l2,2),(:p1,2),(:c1,2),(:c2,2),(:jj,2),(:i1,2),(:gnd,1)]]))
    end

    @testset "vector connection groups" begin
        # groups written as vectors of endpoints are equivalent to tuple
        # groups; with the ground instance spelling the connections are a
        # fully concrete vector of vectors of tuples
        comps = [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0),
                 :cc => Capacitor(100e-15), :jj => JosephsonJunction(1e-9),
                 :cj => Capacitor(40e-15), :gnd => Ground()]
        ctup = Circuit(comps,
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:cj, 1)),
             ((:p1, 2), (:r1, 2), (:jj, 2), (:cj, 2), Ground)])
        vconns = [[(:p1, 1), (:r1, 1), (:cc, 1)],
                  [(:cc, 2), (:jj, 1), (:cj, 1)],
                  [(:p1, 2), (:r1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]]
        @test vconns isa Vector{Vector{Tuple{Symbol,Int}}}
        cvec = Circuit(comps, vconns)
        @test JosephsonCircuits.comparestruct(
            compile(ctup), compile(cvec))
    end

    @testset "scattering block endpoint rules" begin
        S = [0.0 1.0; 1.0 0.0]
        # grounded is the default; a floating block asks for it explicitly
        @test ScatteringParameters(S).grounded
        @test TransmissionLine(50.0, 1e-3).grounded
        blk = ScatteringParameters(S; grounded = false)
        @test blk.nports == 2
        @test JosephsonCircuits.nterminals(blk) == 4
        # floating block: bare (inst, p) in a group is an error
        @test_throws ArgumentError Circuit([:b => blk],
            [((:b, 1), Ground)])
        # explicit terminals work
        cfl = Circuit([:b => blk, :l => Inductor(1e-9)],
            [((:b, 1, 1), (:l, 1)), ((:b, 1, 2), (:l, 2), Ground),
             ((:b, 2, 1), (:b, 2, 2), Ground)])
        @test JosephsonCircuits.ninstances(elaborate(cfl)) == 2

        gblk = ScatteringParameters(S; grounded = true)
        # grounded block: (inst, p) is the signal terminal
        cg = Circuit([:b => gblk, :l => Inductor(1e-9)],
            [((:b, 1), (:l, 1)), ((:l, 2), Ground), ((:b, 2), Ground)])
        elab = elaborate(cg)
        # reference terminals are on the ground net
        bi = findfirst(==("b"), elab.instancepaths)
        terms = JosephsonCircuits.instanceterminals(elab, bi)
        @test terms[2] == 1 && terms[4] == 1
        # explicitly connecting a grounded reference terminal is an error
        @test_throws ArgumentError Circuit([:b => gblk],
            [((:b, 1, 2), Ground)])
        # lowering emits the block once, and nothing in the flat table
        psc = compile(cg)
        @test !any(==(:S), psc.componenttypes)
        b = only(psc.scatteringblocks)
        @test b.definition === gblk
        @test b.path == "b"
        @test length(b.signalnodes) == 2
        # reference terminals lower to the ground node
        @test all(==(1), b.refnodes)
    end

    @testset "voltage sources and ground requirement" begin
        c = Circuit([:v => VoltageSource(1.0), :r => Resistor(50.0)],
            [((:v, 1), (:r, 1)), ((:v, 2), (:r, 2), Ground)])
        @test_throws ComponentNotSupportedError compile(c)
        # no ground
        cng = Circuit([:l => Inductor(1e-9), :c => Capacitor(1e-12)],
            [((:l, 1), (:c, 1)), ((:l, 2), (:c, 2))])
        @test_throws ArgumentError compile(cng)
    end

    @testset "symbolic values through the typed path" begin
        @variables Lj Cc Cj
        circuit = Tuple{String,String,String,Num}[
            ("P1","1","0",1),
            ("R1","1","0",50.0),
            ("C1","1","2",Cc),
            ("Lj1","2","0",Lj),
            ("C2","2","0",Cj),
        ]
        circuitdefs = Dict(Lj => 1000e-12, Cc => 100e-15, Cj => 1000e-15)
        native = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :c1 => Capacitor(Cc),
             :jj => NonlinearInductor(Lj, sin, cos),
             :c2 => Capacitor(Cj)],
            [((:p1, 1), (:r1, 1), (:c1, 1)),
             ((:c1, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:c2, 2), (:r1, 2), (:p1, 2), Ground)],
        )
        ws = 2*pi*(4.5:0.1:5.0)*1e9
        wp = (2*pi*4.75001e9,)
        sources = [(mode=(1,), port=1, current=0.00565e-6)]
        sol1 = hbsolve(ws, wp, sources, (8,), (8,), circuit, circuitdefs)
        sol2 = hbsolve(ws, wp, sources, (8,), (8,), native, circuitdefs)
        @test isapprox(sol1.linearized.S((0,), 1, (0,), 1, :),
            sol2.linearized.S((0,), 1, (0,), 1, :); rtol = 1e-6)
    end

    @testset "calcnodesorting" begin
        @test_throws(
            ArgumentError("Unknown sorting type."),
            JosephsonCircuits.calcnodesorting(["30","11","0","2"];sorting=:test)
        )

        @test_throws(
            ArgumentError("No ground node found in netlist."),
            JosephsonCircuits.calcnodesorting(["30","11","1","2"];sorting=:none)
        )

        @test_throws(
            ArgumentError("No ground node found in netlist."),
            JosephsonCircuits.calcnodesorting(String[];sorting=:none)
        )

        @test_throws(
            ArgumentError("Failed to parse the nodes as integers. Try setting the keyword argument `sorting=:name` or `sorting=:none`."),
            JosephsonCircuits.calcnodesorting(["30","11","0","a"];sorting=:number)
        )
    end

    @testset "calcnodesorting negative node names" begin
        # the returned indices must always be a permutation with ground first.
        # ground's index in the input is not in general its position in the
        # sorting permutation, and using one for the other duplicated an index
        # and dropped another whenever ground did not sort first, which happens
        # as soon as any node name is a negative number.
        function isgroundfirstperm(nodes, indices)
            sort(indices) == collect(1:length(nodes)) && nodes[indices][1] == "0"
        end

        # smallest failing case: ground first in the input, one negative node
        @test JosephsonCircuits.calcnodesorting(["0","-2"];sorting=:number) == [1,2]
        @test JosephsonCircuits.calcnodesorting(["0","-1"];sorting=:name) == [1,2]

        # ground neither first in the input nor first after sorting
        @test JosephsonCircuits.calcnodesorting(["-2","-1","5","0"];sorting=:number) ==
            [4,1,2,3]

        for sorting in (:number, :name, :none)
            for nodes in (["0","-2"], ["0","-1","2"], ["-2","-1","5","0"],
                    ["-1","0","3","-4"], ["0","-3","-1","1","2"], ["3","-1","0"])
                indices = JosephsonCircuits.calcnodesorting(copy(nodes);sorting=sorting)
                @test isgroundfirstperm(nodes, indices)
            end
        end

        # the non-ground nodes must still be in the requested order
        @test ["-2","-1","5","0"][
            JosephsonCircuits.calcnodesorting(["-2","-1","5","0"];sorting=:number)] ==
            ["0","-2","-1","5"]

        # a negative node name must survive parsing rather than being silently
        # renamed to ground, which is what the corrupted permutation did
        @variables R Cc Lj Cj
        circuit = [("P1","0","-1",1),("R1","0","-1",R),("C1","-1","2",Cc),
            ("Lj1","2","0",Lj),("C2","2","0",Cj)]
        psc = JosephsonCircuits.compile(circuit;sorting=:number)
        @test psc.nodenames == ["0","-1","2"]
        @test length(unique(psc.nodenames)) == length(psc.nodenames)
    end
end
