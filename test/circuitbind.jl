using JosephsonCircuits
using LinearAlgebra
using Test

# a function barrier for allocation measurement: a closure defined inside a
# testset captures test locals, and boxing those is counted as the kernel's
refillallocs(nz, seen, p, v, n) = @allocated (for _ in 1:n
    JosephsonCircuits.assemblenodal!(nz, seen, p, v)
end)

@testset verbose=true "binding and assembly" begin

    @testset "bound circuit and nodal assembly" begin
        JC = JosephsonCircuits
        # grounded and floating components, two components in parallel on one
        # node pair, a lossy capacitor, and two ports with different Z0
        c = Circuit(
            Any[:p1 => Port(1), :p2 => Port(2; Z0 = 1000.0),
                :c1 => Capacitor(1e-13), :c2 => Capacitor(2e-13),
                :c3 => Capacitor(3e-13 + 1e-16im), :r9 => Resistor(75.0),
                :r8 => Resistor(120.0), :l1 => Inductor(1e-9),
                :jj => JosephsonJunction(1e-9), :gnd => Ground()],
            Any[[(:p1,1),(:c1,1),(:c2,1),(:r9,1),(:l1,1)],
                [(:l1,2),(:c3,1),(:jj,1),(:r8,1),(:p2,1)],
                [(:p1,2),(:p2,2),(:c1,2),(:c2,2),(:c3,2),(:r9,2),(:r8,2),
                 (:jj,2),(:gnd,1)]])
        cc = JC.compile(c)
        b = JC.bind(cc)

        # values are grouped, concrete, and in the compiled group order
        @test b.capacitors == [cc.componentvalues[i] for i in cc.capacitors]
        @test b.resistors == [cc.componentvalues[i] for i in cc.resistors]
        @test isconcretetype(eltype(b.resistors))
        @test isconcretetype(eltype(b.capacitors))
        # one lossy capacitor makes the capacitances complex and leaves the
        # resistances real
        @test eltype(b.capacitors) <: Complex
        @test eltype(b.resistors) <: Real
        # the port environments carry the reference impedances
        @test [b.values[p.environment] for p in cc.ports] == [50.0, 1000.0]

        # the planned assembly reproduces the coordinate assembly exactly,
        # including the summation order of parallel components
        vvn = JC.componentvaluestonumber(cc.componentvalues, Dict{Any,Any}())
        pC = JC.nodalstampplan(cc, cc.capacitors, cc.Nnodes)
        pG = JC.nodalstampplan(cc, cc.resistors, cc.Nnodes; invert = true)
        same(A, B) = A.colptr == B.colptr && A.rowval == B.rowval &&
            A.nzval == B.nzval
        for Nmodes in (1, 4, 8)
            Cref = JC.calcCn(cc.componenttypes, cc.nodeindices, vvn, Nmodes,
                cc.Nnodes)
            Gref = JC.calcGn(cc.componenttypes, cc.nodeindices, vvn, Nmodes,
                cc.Nnodes)
            @test same(JC.assemblenodal(eltype(Cref), pC, b.capacitors,
                Nmodes), Cref)
            @test same(JC.assemblenodal(eltype(Gref), pG, b.resistors,
                Nmodes), Gref)
        end

        # refilling reuses the pattern: the same plan at new values agrees
        # with a full rebuild at those values
        b2 = JC.bind(JC.compile(Circuit(
            Any[:p1 => Port(1), :p2 => Port(2; Z0 = 1000.0),
                :c1 => Capacitor(5e-13), :c2 => Capacitor(7e-13),
                :c3 => Capacitor(1e-12 + 3e-16im), :r9 => Resistor(75.0),
                :r8 => Resistor(120.0), :l1 => Inductor(1e-9),
                :jj => JosephsonJunction(1e-9), :gnd => Ground()],
            Any[[(:p1,1),(:c1,1),(:c2,1),(:r9,1),(:l1,1)],
                [(:l1,2),(:c3,1),(:jj,1),(:r8,1),(:p2,1)],
                [(:p1,2),(:p2,2),(:c1,2),(:c2,2),(:c3,2),(:r9,2),(:r8,2),
                 (:jj,2),(:gnd,1)]])))
        nz = Vector{ComplexF64}(undef, length(pC.rowval))
        sn = Vector{Bool}(undef, length(pC.rowval))
        JC.assemblenodal!(nz, sn, pC, b2.capacitors)
        @test nz == JC.assemblenodal(ComplexF64, pC, b2.capacitors, 1).nzval
        # a refill allocates nothing
        refillallocs(nz, sn, pC, b2.capacitors, 1)
        @test refillallocs(nz, sn, pC, b2.capacitors, 100) == 0

        # a structural key distinguishes a value change from a topology
        # change: the same circuit at new numbers keeps its key, an infinite
        # inductance does not
        @test JC.structuralkey(b) == JC.structuralkey(b2)
        bopen = JC.bind(JC.compile(Circuit(
            [:p1 => Port(1), :l1 => Inductor(Inf), :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:l1,1),(:c1,1)], [(:p1,2),(:l1,2),(:c1,2),Ground]])))
        bfinite = JC.bind(JC.compile(Circuit(
            [:p1 => Port(1), :l1 => Inductor(1e-9), :c1 => Capacitor(1e-12)],
            [[(:p1,1),(:l1,1),(:c1,1)], [(:p1,2),(:l1,2),(:c1,2),Ground]])))
        @test JC.structuralkey(bopen) != JC.structuralkey(bfinite)
    end

    @testset "planned circuit matrices" begin
        JC = JosephsonCircuits
        same(a, b) = a.colptr == b.colptr && a.rowval == b.rowval &&
            a.nzval == b.nzval
        samev(a, b) = a.n == b.n && a.nzind == b.nzind && a.nzval == b.nzval
        function matricesagree(c; Nmodes = 8)
            cc = JC.compile(c); b = JC.bind(cc)
            cg = calccircuitgraph(cc)
            vvn = JC.componentvaluestonumber(cc.componentvalues,
                Dict{Any,Any}())
            ref = numericmatrices(cc, cg, vvn; Nmodes = Nmodes)
            new = JC.assemblematrices(
                JC.circuitmatrixplan(cc, cg, b; Nmodes = Nmodes), b)
            return same(new.Cnm, ref.Cnm) && same(new.Gnm, ref.Gnm) &&
                same(new.invLnm, ref.invLnm) && same(new.Mb, ref.Mb) &&
                same(new.Rbnm, ref.Rbnm) &&
                samev(new.Lb, ref.Lb) && samev(new.Lbm, ref.Lbm) &&
                samev(new.Ljb, ref.Ljb) && samev(new.Ljbm, ref.Ljbm) &&
                new.Lmean == ref.Lmean &&
                new.portindices == ref.portindices &&
                new.portnumbers == ref.portnumbers &&
                new.portimpedances == ref.portimpedances &&
                new.portenvironmentindices == ref.portenvironmentindices &&
                new.noiseportimpedanceindices == ref.noiseportimpedanceindices
        end

        # a legacy netlist
        @test matricesagree(Circuit([("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","2",100e-15), ("Lj1","2","0",1e-9), ("C2","2","0",1e-12)]))

        # two inductors on one branch combine as a parallel inductance, and a
        # complex capacitance must not make the resistances complex
        @test matricesagree(Circuit(
            Any[:p1 => Port(1), :la => Inductor(2e-9), :lb => Inductor(3e-9),
                :l2 => Inductor(5e-9), :jj => JosephsonJunction(1e-9),
                :c1 => Capacitor(1e-13 + 1e-16im), :gnd => Ground()],
            Any[[(:p1,1),(:la,1),(:lb,1),(:c1,1)],
                [(:la,2),(:lb,2),(:l2,1),(:jj,1)],
                [(:l2,2),(:jj,2),(:p1,2),(:c1,2),(:gnd,1)]]))

        # a node where four inductive branches meet, so the inverse
        # inductance entries take more than two contributions and their
        # summation order is observable
        @test matricesagree(Circuit(
            Any[:p1 => Port(1), :la => Inductor(1e-9), :lb => Inductor(2.5e-9),
                :lc => Inductor(7e-9), :ld => Inductor(0.3e-9),
                :ca => Capacitor(1e-12), :cb => Capacitor(2e-12),
                :cc => Capacitor(3e-12), :gnd => Ground()],
            Any[[(:p1,1),(:la,1),(:lb,1),(:lc,1),(:ld,1)],
                [(:la,2),(:ca,1)], [(:lb,2),(:cb,1)], [(:lc,2),(:cc,1)],
                [(:ld,2),(:p1,2),(:ca,2),(:cb,2),(:cc,2),(:gnd,1)]]))

        # mutually coupled branches are dropped from the inverse inductance
        # matrix and carried as auxiliary MNA currents instead, so no
        # inductance matrix is inverted anywhere in the assembly
        @test matricesagree(Circuit(
            Any[:p1 => Port(1), :l1 => Inductor(1e-9), :l2 => Inductor(2e-9),
                :l3 => Inductor(4e-9), :k => MutualInductor(0.9, :l1, :l2),
                :c1 => Capacitor(1e-12), :gnd => Ground()],
            Any[[(:p1,1),(:l1,1),(:c1,1)], [(:l1,2),(:l2,1)],
                [(:l2,2),(:l3,1)],
                [(:l3,2),(:p1,2),(:c1,2),(:gnd,1)]]))

        # and a circuit with no inductance at all
        @test matricesagree(Circuit(
            [:p1 => Port(1), :c1 => Capacitor(1e-12), :jj => JosephsonJunction(1e-9)],
            [[(:p1,1),(:c1,1),(:jj,1)], [(:p1,2),(:c1,2),(:jj,2),Ground]]))

        # at one mode as well as many
        @test matricesagree(Circuit([("P1","1","0",1), ("R1","1","0",50.0),
            ("L1","1","2",1e-9), ("C2","2","0",1e-12)]); Nmodes = 1)
    end
end
