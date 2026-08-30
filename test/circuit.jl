using Symbolics
using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "typed circuit representation" begin

    @testset "legacy tuple round trip" begin
        @variables Ipump Rleft Cc Lj Cj
        circuit = Tuple{String,String,String,Num}[
            ("P1","1","0",1),
            ("I1","1","0",Ipump),
            ("R1","1","0",Rleft),
            ("C1","1","2",Cc),
            ("Lj1","2","0",Lj),
            ("C2","2","0",Cj),
        ]
        psclegacy = parsesortcircuit(circuit; sorting = :number)
        pscnew = parsesortcircuit(Circuit(circuit); sorting = :number)
        @test JosephsonCircuits.comparestruct(psclegacy, pscnew)

        # with name sorting
        @test JosephsonCircuits.comparestruct(
            parsesortcircuit(circuit; sorting = :name),
            parsesortcircuit(Circuit(circuit); sorting = :name))
    end

    @testset "legacy round trip with mutual inductors" begin
        @variables Ipump Rleft L1v L2v C2v
        circuit = Tuple{String,String,String,Num}[
            ("P1","1","0",1),
            ("I1","1","0",Ipump),
            ("R1","1","0",Rleft),
            ("L1","1","0",L1v),
            ("K1","L1","L2",0.9),
            ("L2","2","0",L2v),
            ("C2","2","0",C2v),
        ]
        psclegacy = parsesortcircuit(circuit; sorting = :number)
        pscnew = parsesortcircuit(Circuit(circuit); sorting = :number)
        @test JosephsonCircuits.comparestruct(psclegacy, pscnew)
        @test pscnew.mutualinductorbranchnames == ["L1","L2"]
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
        psc = parsesortcircuit(c; sorting = :number)
        @test psc.componenttypes == [:P, :R, :NL, :C]
        @test psc.componentvalues[2] == 50.0
        @test psc.componentvalues[3] == 2.0
        # without substitution the symbols pass through
        psc2 = parsesortcircuit(Circuit(circuit); sorting = :number)
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

        native = Circuit(
            [:p1 => Port(1), :r1 => Resistor(50.0),
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
        push!(components, :p1 => Port(1))
        push!(components, :r1 => Resistor(50.0))
        push!(components, :p2 => Port(2))
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
            [:a => r1, :b => r1, :p1 => Port(1), :rt => Resistor(50.0)],
            [PortRef(:a, 1) => PortRef(:b, 1),
             ((:p1, 1), (:rt, 1), (:a, 1)),
             ((:p1, 2), (:rt, 2), (:a, 2), Ground)],
        )
        elab = elaborate(c)
        # port bond merges a/l terminal1 with b/l terminal1, and a's pin 2
        # net with b's pin 2 net
        @test JosephsonCircuits.ninstances(elab) == 6
        psc = parsesortcircuit(elab)
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
            Any[:p1 => Port(1), :r1 => Resistor(50.0),
                :cc => Capacitor(100e-15), :jj => JosephsonJunction(1e-9),
                :cj => Capacitor(100e-15),
                (g ? [:gnd => Ground()] : [])...],
            [((:p1, 1), (:r1, 1), (:cc, 1)),
             ((:cc, 2), (:jj, 1), (:cj, 1)),
             g ? ((:p1, 2), (:r1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)) :
                 ((:p1, 2), (:r1, 2), (:jj, 2), (:cj, 2), Ground)])
        @test JosephsonCircuits.comparestruct(
            parsesortcircuit(mk(false)), parsesortcircuit(mk(true)))

        # the ground instance flattens to no component
        elab = JosephsonCircuits.elaborate(mk(true))
        @test JosephsonCircuits.ninstances(elab) == 5
        @test !any(d -> d isa JosephsonCircuits.GroundType, elab.definitions)

        # pair syntax, PinRef, and both definition spellings
        cpair = Circuit(
            [:r1 => Resistor(50.0), :p1 => Port(1),
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
            [:p1 => Port(1), :r1 => Resistor(50.0), :s1 => sub, :s2 => sub,
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
            [:p => Port(1), :r => Resistor(50.0), :s => sub],
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

    @testset "vector connection groups" begin
        # groups written as vectors of endpoints are equivalent to tuple
        # groups; with the ground instance spelling the connections are a
        # fully concrete vector of vectors of tuples
        comps = [:p1 => Port(1), :r1 => Resistor(50.0),
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
            parsesortcircuit(ctup), parsesortcircuit(cvec))
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
        # lowering emits one :S component per port, sharing the block
        psc = parsesortcircuit(cg)
        sindices = findall(==(:S), psc.componenttypes)
        @test length(sindices) == 2
        @test psc.componentnames[sindices] == ["b/port1", "b/port2"]
        @test all(psc.componentvalues[i].block === gblk for i in sindices)
        @test [psc.componentvalues[i].port for i in sindices] == [1, 2]
        # reference terminals lower to the ground node
        @test all(psc.nodeindices[2, i] == 1 for i in sindices)
    end

    @testset "matrix providers and evaluation" begin
        # constant provider with conjugate symmetry
        S = [0.0 0.8im; 0.8im 0.0]
        blk = ScatteringParameters(S)
        ws = [-2*pi*5e9, 2*pi*5e9]
        dest = zeros(Complex{Float64}, 2, 2, 2)
        JosephsonCircuits.evaluatescattering!(dest, blk, ws)
        @test dest[:,:,2] == S
        @test dest[:,:,1] == conj.(S)

        # tabulated provider: interpolation and extrapolation policies
        freqs = [1.0, 2.0, 3.0]
        vals = zeros(Complex{Float64}, 1, 1, 3)
        vals[1,1,:] .= [0.1, 0.3, 0.5]
        tblk = ScatteringParameters((freqs, vals))
        d = zeros(Complex{Float64}, 1, 1, 1)
        JosephsonCircuits.evaluatescattering!(d, tblk, [1.5])
        @test d[1,1,1] ≈ 0.2
        @test_throws ArgumentError JosephsonCircuits.evaluatescattering!(
            d, tblk, [4.0])
        tconst = ScatteringParameters((freqs, vals); extrapolation = :constant)
        JosephsonCircuits.evaluatescattering!(d, tconst, [4.0])
        @test d[1,1,1] ≈ 0.5
        tlin = ScatteringParameters((freqs, vals); extrapolation = :linear)
        JosephsonCircuits.evaluatescattering!(d, tlin, [4.0])
        @test d[1,1,1] ≈ 0.7
        # strictly increasing frequencies required
        @test_throws ArgumentError ScatteringParameters(([2.0, 1.0],
            vals[:,:,1:2]))

        # callable provider requires nports
        f(w) = [0.0 exp(-im*w*1e-12); exp(-im*w*1e-12) 0.0]
        @test_throws ArgumentError ScatteringParameters(f)
        cblk = ScatteringParameters(f; nports = 2)
        d2 = zeros(Complex{Float64}, 2, 2, 1)
        JosephsonCircuits.evaluatescattering!(d2, cblk, [2*pi*5e9])
        @test d2[2,1,1] ≈ exp(-im*2*pi*5e9*1e-12)

        # passivity validation
        @test_throws ArgumentError ScatteringParameters([0.0 2.0; 2.0 0.0])
        active = ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([1.0 0.0; 0.0 1.0]))
        @test active.nports == 2
        # noise covariance must be Hermitian
        @test_throws ArgumentError ScatteringParameters([0.0 2.0; 2.0 0.0];
            noise = NoiseCovariance([1.0 1.0; 0.0 1.0]))
        # thermal equilibrium noise model carries the temperature
        blkT = ScatteringParameters(S; noise = ThermalEquilibrium(20e-3))
        @test blkT.noise.temperature == 20e-3

        # reference impedance handling
        @test ScatteringParameters(S; zref = 30.0).zref == [30.0, 30.0]
        @test ScatteringParameters(S; zref = [50.0, 75.0]).zref == [50.0, 75.0]
        @test_throws DimensionMismatch ScatteringParameters(S; zref = [50.0])
        @test_throws ArgumentError ScatteringParameters(S; zref = -50.0)
    end

    @testset "transmission line" begin
        Z0, len = 50.0, 1e-3
        tl = TransmissionLine(Z0, len)
        @test tl isa ScatteringParameters
        @test tl.zref == [Z0, Z0]
        @test tl.negative_frequency isa Native
        w = 2*pi*5e9
        d = zeros(Complex{Float64}, 2, 2, 2)
        JosephsonCircuits.evaluatescattering!(d, tl, [w, -w])
        delay = len/JosephsonCircuits.speed_of_light
        @test d[2,1,1] ≈ exp(-im*w*delay)
        @test d[1,1,1] == 0
        # native evaluation satisfies the conjugation identity by construction
        @test d[:,:,2] ≈ conj.(d[:,:,1])
    end

    @testset "gaussian channels" begin
        # attenuator at zero temperature: exactly at the CP boundary
        η = 0.5
        ch = GaussianChannel(sqrt(η)*Matrix(1.0I, 2, 2),
            (1-η)/2*Matrix(1.0I, 2, 2); nmodes = 1)
        @test abs(ch.cp_margin) < 1e-10
        @test ch.nmodes == 1
        @test JosephsonCircuits.nterminals(ch) == 2
        # quantum limited phase insensitive amplifier
        G = 4.0
        amp = GaussianChannel(sqrt(G)*Matrix(1.0I, 2, 2),
            (G-1)/2*Matrix(1.0I, 2, 2); nmodes = 1)
        @test abs(amp.cp_margin) < 1e-10
        # a noiseless amplifier is not completely positive
        @test_throws ArgumentError GaussianChannel(
            sqrt(G)*Matrix(1.0I, 2, 2), zeros(2, 2); nmodes = 1)
        # Y must be symmetric
        @test_throws ArgumentError GaussianChannel(Matrix(1.0I, 2, 2),
            [0.5 0.1; -0.1 0.5]; nmodes = 1)
        # odd dimension is rejected
        @test_throws DimensionMismatch GaussianChannel(zeros(3,3),
            zeros(3,3))
        # ideal squeezer: symplectic X, Y = 0 is completely positive
        r = 0.5
        sq = GaussianChannel([exp(r) 0.0; 0.0 exp(-r)], zeros(2,2);
            nmodes = 1)
        @test abs(sq.cp_margin) < 1e-10
        # Bogoliubov conversion: a beamsplitter swap
        X = quadraturetransform([0 1; 1 0], zeros(2, 2))
        @test X == [0 1 0 0; 1 0 0 0; 0 0 0 1; 0 0 1 0]
        # two mode channel from the Bogoliubov form of a two mode squeezer
        A = cosh(r)*Matrix(1.0I, 2, 2)
        B = sinh(r)*[0.0 1.0; 1.0 0.0]
        Xtms = quadraturetransform(A, B)
        tms = GaussianChannel(Xtms, zeros(4, 4); nmodes = 2)
        @test abs(tms.cp_margin) < 1e-10
        # channels embed in circuits and are rejected by the solver bridge
        cg = Circuit([:ch => GaussianChannel(sqrt(η)*Matrix(1.0I, 2, 2),
                (1-η)/2*Matrix(1.0I, 2, 2); nmodes = 1, grounded = true)],
            [((:ch, 1), Ground)])
        @test_throws ComponentNotSupportedError parsesortcircuit(cg)
    end

    @testset "nonlinear inductors and current phase relations" begin
        jj = JosephsonJunction(100e-12)
        @test jj isa NonlinearInductor
        @test JosephsonCircuits.issinusoidal(jj)
        @test JosephsonJunction(Ic = JosephsonCircuits.phi0/100e-12).L0 ≈ 100e-12

        p = PolynomialCPR([1.0, 0.0, -1/6])
        @test p(0.1) ≈ 0.1 - 0.1^3/6
        dp = JosephsonCircuits.cprderivative(p)
        @test dp(0.0) ≈ 1.0
        @test dp(0.2) ≈ 1.0 - 0.2^2/2
        # the linear coefficient must be one
        @test_throws ArgumentError PolynomialCPR([2.0, 0.0])
        @test_throws ArgumentError PolynomialCPR(Float64[])
        # unknown callables require an explicit derivative
        mycpr(x) = x - x^3/6
        @test_throws ArgumentError NonlinearInductor(1e-9, mycpr)
        nl = NonlinearInductor(1e-9, mycpr, x -> 1 - x^2/2)
        @test !JosephsonCircuits.issinusoidal(nl)
        # snail-style direct specification elaborates but does not yet lower
        c = Circuit([:snail => NonlinearInductor(1e-9, p),
                     :r => Resistor(50.0), :p1 => Port(1)],
            [((:p1, 1), (:r, 1), (:snail, 1)),
             ((:p1, 2), (:r, 2), (:snail, 2), Ground)])
        @test JosephsonCircuits.ninstances(elaborate(c)) == 3
        @test_throws ComponentNotSupportedError parsesortcircuit(c)
        # a sinusoidal junction lowers to the legacy Lj component
        c2 = Circuit([:jj => JosephsonJunction(100e-12), :p1 => Port(1),
                      :r => Resistor(50.0)],
            [((:p1, 1), (:r, 1), (:jj, 1)),
             ((:p1, 2), (:r, 2), (:jj, 2), Ground)])
        psc = parsesortcircuit(c2)
        @test psc.componenttypes[findfirst(==("jj"),
            psc.componentnames)] == :Lj
    end

    @testset "voltage sources and ground requirement" begin
        c = Circuit([:v => VoltageSource(1.0), :r => Resistor(50.0)],
            [((:v, 1), (:r, 1)), ((:v, 2), (:r, 2), Ground)])
        @test_throws ComponentNotSupportedError parsesortcircuit(c)
        # no ground
        cng = Circuit([:l => Inductor(1e-9), :c => Capacitor(1e-12)],
            [((:l, 1), (:c, 1)), ((:l, 2), (:c, 2))])
        @test_throws ArgumentError parsesortcircuit(cng)
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
            [:p1 => Port(1), :r1 => Resistor(50.0), :c1 => Capacitor(Cc),
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

    @testset "touchstone loading" begin
        path = joinpath(mktempdir(), "attenuator.s2p")
        open(path, "w") do io
            write(io, "# GHz S MA R 50\n")
            write(io, "1.0 0.0 0.0 0.5 0.0 0.5 0.0 0.0 0.0\n")
            write(io, "2.0 0.0 0.0 0.5 0.0 0.5 0.0 0.0 0.0\n")
        end
        blk = ScatteringParameters(path)
        @test blk.nports == 2
        @test blk.zref == [50.0, 50.0]
        d = zeros(Complex{Float64}, 2, 2, 1)
        JosephsonCircuits.evaluatescattering!(d, blk, [2*pi*1.5e9])
        @test d[2,1,1] ≈ 0.5
        # a conflicting explicit zref is an error
        @test_throws ArgumentError ScatteringParameters(path; zref = 30.0)
    end
end
