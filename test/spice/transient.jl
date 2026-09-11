using JosephsonCircuits
using Test
using XicTools_jll

# The WRspice back end of the transient without the executable: the
# options, the refusals, and the input a run writes. The runs themselves,
# against the package's own stepping, are in wrspicecrosscheck.jl, which
# needs the wrspice executable.
@testset "the transient through WRspice" begin
    JC = JosephsonCircuits

    @testset "the options" begin
        m = WRspice()
        @test isnothing(m.executable)
        @test m.dphimax == 0.01
        @test m.jjaccel
        @test m.maxdata == 2e9
        @test WRspice(dphimax = 0.1).dphimax == 0.1
        @test !WRspice(jjaccel = false).jjaccel
        @test_throws ArgumentError WRspice(dphimax = 0.0)
        @test_throws ArgumentError WRspice(dphimax = 4.0)
        @test_throws ArgumentError WRspice(maxdata = 10.0)
    end

    circuit = Circuit(
        [:P1 => Port(1; Z0 = 50.0),
         :C1 => Capacitor(100e-15),
         :Lj1 => JosephsonJunction(1000e-12),
         :C2 => Capacitor(1000e-15),
         :gnd => Ground()],
        [Net("1", [(:P1, 1), (:C1, 1)]),
         Net("2", [(:C1, 2), (:Lj1, 1), (:C2, 1)]),
         Net("0", [(:P1, 2), (:Lj1, 2), (:C2, 2), (:gnd, 1)])])
    p = transientproblem(circuit;
        sources = [TransientSource(1, t -> 1e-6*sin(2pi*5e9*t))])

    @testset "the refusals" begin
        run(; kwargs...) = transientsolve(p, (0.0, 1e-9); dt = 1e-12,
            method = WRspice(), kwargs...)
        @test_throws ArgumentError run(record = :states)
        @test_throws ArgumentError run(record = :checkpoints)
        @test_throws ArgumentError run(linearsolver = GMRES())
        @test_throws ArgumentError run(factorization = KLUfactorization())
        @test_throws ArgumentError run(reuse = TransientReuse())
        @test_throws ArgumentError run(checkpointevery = 2)
        @test_throws ArgumentError run(
            initialstate = transientstate(p; voltage = fill(1e-6, 2)))
        # the saved grid must be a whole number of saves
        @test_throws ArgumentError run(saveevery = 3)
        # a constant scattering block steps natively but has no SPICE
        # element
        blocked = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:att, 1, 2, ScatteringParameters([0.0 0.5; 0.5 0.0]; zref = 50.0)),
            (:r1, 2, 0, Resistor(50.0)), (:c1, 2, 0, Capacitor(1e-15))])
        pb = transientproblem(blocked)
        @test_throws ArgumentError transientsolve(pb, (0.0, 1e-9);
            dt = 1e-12, method = WRspice())
        # a current phase relation beyond the sinusoidal one
        bent = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:c1, 1, 2, Capacitor(100e-15)),
            (:lj, 2, 0, NonlinearInductor(1000e-12,
                PolynomialCPR([1.0, 0.0, -1/6, 0.0, 1/120]))),
            (:cj, 2, 0, Capacitor(1000e-15))])
        pc = transientproblem(bent)
        @test_throws ArgumentError transientsolve(pc, (0.0, 1e-9);
            dt = 1e-12, method = WRspice())
    end

    @testset "the input of a run" begin
        input, junctions = JC.wrspiceinput(p, WRspice(), 0.0, 1e-12, 10, 1e-12)
        # the netlist with its jj model, the drive as a piecewise linear
        # source into the port's positive node, and the control block
        @test occursin("jjk", input)
        @test occursin("isrcd1 0 1 pwl(0.0 0.0 ", input)
        @test occursin(".tran 1.0e-12 1.0e-11 uic", input)
        @test occursin("set dphimax=0.01", input)
        @test occursin("set jjaccel=1", input)
        @test occursin("set maxdata=2.0e9", input)
        # the run keeps the port terminals and nothing else; ground is
        # not a trace and the junction's phase node comes with :phases
        @test occursin("\nsave v(1)\nrun", input)
        inputph, _ = JC.wrspiceinput(p, WRspice(), 0.0, 1e-12, 10, 1e-12, :phases)
        @test occursin("\nsave v(1) v(3)\nrun", inputph)
        @test length(junctions) == 1
        @test junctions[1].phasenode == "3"
        @test p.circuit.componentnames[junctions[1].index] == "Lj1"
        # the piecewise linear samples are the waveform's on the grid
        @test occursin(string(1e-12, " ", 1e-6*sin(2pi*5e9*1e-12)), input)
        # a waveform given as a number is a constant source, and jjaccel
        # can be left off
        pconst = transientproblem(circuit;
            sources = [TransientSource(1, 1.5e-6)])
        inputc, _ = JC.wrspiceinput(pconst, WRspice(jjaccel = false),
            0.0, 1e-12, 10, 1e-12)
        @test occursin("isrcd1 0 1 1.5e-6", inputc)
        @test !occursin("jjaccel", inputc)
        # a netlist current source no drive replaces is a constant source
        biased = Circuit(
            [:P1 => Port(1; Z0 = 50.0),
             :C1 => Capacitor(100e-15),
             :Lj1 => JosephsonJunction(1000e-12),
             :C2 => Capacitor(1000e-15),
             :I1 => CurrentSource(1e-8),
             :gnd => Ground()],
            [Net("1", [(:P1, 1), (:C1, 1)]),
             Net("2", [(:C1, 2), (:Lj1, 1), (:C2, 1), (:I1, 1)]),
             Net("0", [(:P1, 2), (:Lj1, 2), (:C2, 2), (:I1, 2), (:gnd, 1)])])
        pbias = transientproblem(biased)
        inputb, _ = JC.wrspiceinput(pbias, WRspice(), 0.0, 1e-12, 10, 1e-12)
        @test occursin("isrcc2 0 2 -1.0e-8", inputb) ||
            occursin("isrcc2 0 2 1.0e-8", inputb)
        # a net named as an integer past the node count collides with a
        # phase node
        colliding = [("P1", "1", "0", 1), ("R1", "1", "0", 50.0),
            ("C1", "1", "2", 100e-15), ("Lj1", "2", "0", 1000e-12),
            ("C2", "2", "0", 1000e-15), ("C3", "2", "4", 10e-15),
            ("R2", "4", "0", 1e4)]
        pcol = transientproblem(colliding)
        @test_throws ArgumentError JC.wrspiceinput(pcol, WRspice(),
            0.0, 1e-12, 10, 1e-12)
    end

    @testset "the executable through the extension" begin
        # Loading XicTools_jll registers its wrspice as the default
        # executable, through the package's extension; this file loads it
        # at the top, as the other tests of the wrapper do, so that the
        # extension is in place before any testset runs rather than part
        # way through one. The artifact carries a binary only on the
        # platforms WRSPICE is built for, and on the others the product
        # is not defined, the extension registers nothing, and an
        # executable has to be installed or given as a path, which
        # wrapper.jl tests.
        if isdefined(XicTools_jll, :wrspice)
            @test !isnothing(JC.wrspicedefaultcmd[])
            @test !isnothing(JC.wrspice_cmd())
        else
            @test isnothing(JC.wrspicedefaultcmd[])
        end
    end

    @testset "a transmission line is the lossless line element" begin
        lined = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
            (:line, 1, 2, TransmissionLine(30.0, 0.1e-9*3e8; vp = 3e8)),
            (:r1, 2, 0, Resistor(30.0)), (:c1, 2, 0, Capacitor(1e-15))])
        n = JosephsonCircuits.exportnetlist(lined, Dict())
        @test occursin("Tline 1 0 2 0 z0=30.0 td=1.0e-10", n.netlist)
        pl = transientproblem(lined)
        inputl, _ = JC.wrspiceinput(pl, WRspice(), 0.0, 1e-12, 10, 1e-12)
        @test occursin("Tline 1 0 2 0", inputl)
    end

    @testset "the netlist export names the junctions" begin
        n = JosephsonCircuits.exportnetlist(circuit, Dict())
        @test length(n.junctions) == 1
        @test n.junctions[1].phasenode == "3"
        @test occursin(" 3 jjk ", n.netlist)
        # without the jj model there are no instances to name
        @test isempty(JosephsonCircuits.exportnetlist(circuit; jj = false).junctions)
        # an infinite resistance is an open and writes no line
        open = [("P1", "1", "0", 1), ("R1", "1", "0", 50.0),
            ("C1", "1", "2", 100e-15), ("Lj1", "2", "0", 1000e-12),
            ("C2", "2", "0", 1000e-15), ("R2", "2", "0", Inf)]
        @test !occursin("R2", JosephsonCircuits.exportnetlist(open, Dict()).netlist)
    end
end
