using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

# a function barrier for allocation measurement: a closure defined inside a
# testset captures test locals, and boxing those is counted as the kernel's
refillallocs(nz, seen, p, v, n) = @allocated (for _ in 1:n
    JosephsonCircuits.assemblenodal!(nz, seen, p, v)
end)

# a compiled circuit bound at its own numeric values
bound(cc) = JosephsonCircuits.bindvalues(cc,
    JosephsonCircuits.componentvaluestonumber(cc.componentvalues, Dict{Any,Any}()))

# The nodal matrix of a group of components written from its definition,
# `f` of each value stamped at its nodes in netlist order, with the modes
# of a node adjacent; ground (node 1) has no row and a component whose
# terminals are one node no stamp.
function coordinatenodal(cc, vvn, group, f, Nmodes)
    n = cc.Nnodes - 1
    T = JosephsonCircuits.grouptype(vvn, group, true)
    I, J, V = Int[], Int[], T[]
    at(node, m) = (node - 1)*Nmodes + m
    for i in group
        n1, n2 = cc.nodeindices[1, i] - 1, cc.nodeindices[2, i] - 1
        n1 == n2 && continue
        v = f(vvn[i])
        for m in 1:Nmodes
            n1 > 0 && (push!(I, at(n1, m)); push!(J, at(n1, m)); push!(V, v))
            n2 > 0 && (push!(I, at(n2, m)); push!(J, at(n2, m)); push!(V, v))
            if n1 > 0 && n2 > 0
                push!(I, at(n1, m)); push!(J, at(n2, m)); push!(V, -v)
                push!(I, at(n2, m)); push!(J, at(n1, m)); push!(V, -v)
            end
        end
    end
    return sparse(I, J, V, Nmodes*n, Nmodes*n)
end

# The mutual inductance matrix written from its definition: K*sqrt(L1*L2)
# between the branches of the two inductors, signed by whether each is
# declared along its branch, which runs from its lower node to its higher.
function coordinatemutual(cc, vvn)
    nb = cc.topology.Nbranches
    I, J, V = Int[], Int[], Float64[]
    along(e) = e[1] < e[2] ? 1 : -1
    for (k, i, j) in cc.couplings
        ei = (cc.nodeindices[1, i], cc.nodeindices[2, i])
        ej = (cc.nodeindices[1, j], cc.nodeindices[2, j])
        bi = cc.topology.edge2indexdict[ei]
        bj = cc.topology.edge2indexdict[ej]
        M = along(ei)*along(ej)*(vvn[k]*sqrt(vvn[i]*vvn[j]))
        append!(I, (bi, bj)); append!(J, (bj, bi)); append!(V, (M, M))
    end
    return sparse(I, J, V, nb, nb)
end

# The circuit matrices written from their definitions: the capacitance, the
# conductance of the resistors and of the ports' terminations, the inverse
# inductance of the inductors on branches no coupling touches, the branch
# inductance of each branch (the reciprocal of its inductors' summed
# reciprocals) and of each junction, repeated per mode, the mutual
# inductances, and the mean of the inductances.
function coordinatematrices(cc, vvn, Nmodes)
    nodes(i) = (cc.nodeindices[1, i], cc.nodeindices[2, i])
    branch(i) = cc.topology.edge2indexdict[nodes(i)]
    nb = cc.topology.Nbranches
    coupled = Set(minmax(nodes(i)...) for (_, i1, i2) in cc.couplings
        for i in (i1, i2))
    uncoupled = [i for i in cc.inductors if !(minmax(nodes(i)...) in coupled)]
    reciprocals = Dict{Int,Float64}()
    for i in cc.inductors
        reciprocals[branch(i)] = get(reciprocals, branch(i), 0.0) + 1/vvn[i]
    end
    Lb = sparsevec(collect(keys(reciprocals)), 1 ./ collect(values(reciprocals)), nb)
    Ljb = sparsevec([branch(i) for i in cc.junctions],
        Float64[vvn[i] for i in cc.junctions], nb)
    repeated(v) = sparsevec([(b - 1)*Nmodes + m for b in v.nzind for m in 1:Nmodes],
        [x for x in v.nzval for m in 1:Nmodes], nb*Nmodes)
    inductances = [vvn[i] for i in vcat(cc.inductors, cc.junctions)]
    return (Cnm = coordinatenodal(cc, vvn, cc.capacitors, identity, Nmodes),
        Gnm = coordinatenodal(cc, vvn, cc.resistors, inv, Nmodes),
        invLnm = coordinatenodal(cc, vvn, uncoupled, inv, Nmodes),
        Lb = Lb, Lbm = repeated(Lb), Ljb = Ljb, Ljbm = repeated(Ljb),
        Mb = coordinatemutual(cc, vvn),
        Lmean = isempty(inductances) ? 0.0 : sum(inductances)/length(inductances))
end

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
        b = bound(cc)

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
        # including the summation order of parallel components, the entries
        # of every component pushed in netlist order
        vvn = JC.componentvaluestonumber(cc.componentvalues, Dict{Any,Any}())
        pC = JC.nodalstampplan(cc, cc.capacitors, cc.Nnodes)
        pG = JC.nodalstampplan(cc, cc.resistors, cc.Nnodes; invert = true)
        same(A, B) = A.colptr == B.colptr && A.rowval == B.rowval &&
            A.nzval == B.nzval
        for Nmodes in (1, 4, 8)
            Cref = coordinatenodal(cc, vvn, cc.capacitors, identity, Nmodes)
            Gref = coordinatenodal(cc, vvn, cc.resistors, inv, Nmodes)
            @test same(JC.assemblenodal(eltype(Cref), pC, b.capacitors,
                Nmodes), Cref)
            @test same(JC.assemblenodal(eltype(Gref), pG, b.resistors,
                Nmodes), Gref)
        end

        # values written as integers assemble in the floating point storage
        # their groups take, so a reciprocal never asks an integer to hold
        # a fraction, and the matrices are those of the same values written
        # as floats
        ints = JC.compile(Circuit([(:p1, 1, 0, Port(1; Z0 = 50)), (:c1, 1, 0, Capacitor(1)),
            (:r1, 1, 0, Resistor(40)), (:l1, 1, 0, Inductor(2))]))
        floats = JC.compile(Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1.0)),
            (:r1, 1, 0, Resistor(40.0)), (:l1, 1, 0, Inductor(2.0))]))
        nmi = numericmatrices(ints, Dict{Any,Any}())
        nmf = numericmatrices(floats, Dict{Any,Any}())
        @test eltype(nmi.Cnm) == eltype(nmi.Gnm) == eltype(nmi.invLnm) == Float64
        @test nmi.invLnm.nzval == [0.5]
        @test same(nmi.Cnm, nmf.Cnm) && same(nmi.Gnm, nmf.Gnm) && same(nmi.invLnm, nmf.invLnm)

        # refilling reuses the pattern: the same plan at new values agrees
        # with a full rebuild at those values
        b2 = bound(JC.compile(Circuit(
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
        # a refill allocates nothing, and still nothing on a circuit an
        # order of magnitude larger: the count is asserted to be zero at
        # each size rather than compared against a figure, which would move
        # with the Julia version while the scaling it protects would not
        refillallocs(nz, sn, pC, b2.capacitors, 1)
        @test refillallocs(nz, sn, pC, b2.capacitors, 100) == 0
        for n in (8, 64)
            comps = Any[:p1 => Port(1), :gnd => Ground()]
            nets = Any[]
            for k in 1:n
                push!(comps, Symbol(:c, k) => Capacitor(1e-13*k),
                    Symbol(:l, k) => Inductor(1e-9), Symbol(:r, k) => Resistor(100.0*k))
            end
            push!(nets, Any[(:p1,1), (:c1,1), (:l1,1), (:r1,1)])
            for k in 1:n-1
                push!(nets, Any[(Symbol(:l,k),2), (Symbol(:c,k+1),1),
                    (Symbol(:l,k+1),1), (Symbol(:r,k+1),1)])
            end
            push!(nets, vcat(Any[(:p1,2), (Symbol(:l,n),2), (:gnd,1)],
                Any[(Symbol(:c,k),2) for k in 1:n], Any[(Symbol(:r,k),2) for k in 1:n]))
            ccn = JC.compile(Circuit(comps, nets))
            bn = bound(ccn)
            pn = JC.nodalstampplan(ccn, ccn.capacitors, ccn.Nnodes)
            nzn = Vector{ComplexF64}(undef, length(pn.rowval))
            snn = Vector{Bool}(undef, length(pn.rowval))
            refillallocs(nzn, snn, pn, bn.capacitors, 1)
            @test refillallocs(nzn, snn, pn, bn.capacitors, 100) == 0
        end
    end

    @testset "planned circuit matrices" begin
        JC = JosephsonCircuits
        same(a, b) = a.colptr == b.colptr && a.rowval == b.rowval &&
            a.nzval == b.nzval
        samev(a, b) = a.n == b.n && a.nzind == b.nzind && a.nzval == b.nzval
        # the planned assembly gives the matrices written from their
        # definitions (see `coordinatematrices`): the same patterns, the
        # capacitance, conductance and junction entries exactly, and the
        # entries formed from reciprocals and sums of inductances to
        # rounding
        function matchesdefinition(c; Nmodes = 8)
            cc = JC.compile(c); b = bound(cc)
            nm = JC.assemblematrices(
                JC.circuitmatrixplan(cc; Nmodes = Nmodes), b)
            ref = coordinatematrices(cc, b.values, Nmodes)
            return same(nm.Cnm, ref.Cnm) && same(nm.Gnm, ref.Gnm) &&
                nm.invLnm.colptr == ref.invLnm.colptr &&
                nm.invLnm.rowval == ref.invLnm.rowval &&
                nm.invLnm.nzval ≈ ref.invLnm.nzval && nm.Mb == ref.Mb &&
                nm.Lb.nzind == ref.Lb.nzind && nm.Lb ≈ ref.Lb &&
                nm.Lbm.nzind == ref.Lbm.nzind && nm.Lbm ≈ ref.Lbm &&
                samev(nm.Ljb, ref.Ljb) && samev(nm.Ljbm, ref.Ljbm) &&
                nm.Lmean ≈ ref.Lmean
        end
        # the same matrices refilled at other values are the matrices a
        # fresh assembly gives at those values, exactly, in the storage
        # they already had
        function refillagrees(c; Nmodes = 8)
            cc = JC.compile(c); b = bound(cc)
            plan = JC.circuitmatrixplan(cc; Nmodes = Nmodes)
            nm = JC.assemblematrices(plan, b)
            vvn = JC.componentvaluestonumber(cc.componentvalues,
                Dict{Any,Any}())
            b2 = JC.bindvalues(cc, [v isa Number ? 1.3*v : v for v in vvn])
            ref = JC.assemblematrices(plan, b2)
            new = JC.assemblematrices!(nm, plan, b2)
            return same(new.Cnm, ref.Cnm) && same(new.Gnm, ref.Gnm) &&
                same(new.invLnm, ref.invLnm) && same(new.Mb, ref.Mb) &&
                samev(new.Lb, ref.Lb) && samev(new.Lbm, ref.Lbm) &&
                samev(new.Ljb, ref.Ljb) && samev(new.Ljbm, ref.Ljbm) &&
                new.Lmean == ref.Lmean &&
                new.portimpedances == ref.portimpedances &&
                new.Cnm === nm.Cnm && new.Gnm === nm.Gnm &&
                new.invLnm === nm.invLnm && new.Lbm === nm.Lbm &&
                new.Ljbm === nm.Ljbm && new.Rbnm === nm.Rbnm
        end

        # a netlist with string names
        @test matchesdefinition(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)), ("Lj1", "2", "0", JosephsonJunction(1e-9)), ("C2", "2", "0", Capacitor(1e-12))]))

        # two inductors on one branch combine as a parallel inductance, and a
        # complex capacitance must not make the resistances complex
        @test matchesdefinition(Circuit(
            Any[:p1 => Port(1), :la => Inductor(2e-9), :lb => Inductor(3e-9),
                :l2 => Inductor(5e-9), :jj => JosephsonJunction(1e-9),
                :c1 => Capacitor(1e-13 + 1e-16im), :gnd => Ground()],
            Any[[(:p1,1),(:la,1),(:lb,1),(:c1,1)],
                [(:la,2),(:lb,2),(:l2,1),(:jj,1)],
                [(:l2,2),(:jj,2),(:p1,2),(:c1,2),(:gnd,1)]]))

        # a node where four inductive branches meet, so the inverse
        # inductance entries take more than two contributions
        @test matchesdefinition(Circuit(
            Any[:p1 => Port(1), :la => Inductor(1e-9), :lb => Inductor(2.5e-9),
                :lc => Inductor(7e-9), :ld => Inductor(0.3e-9),
                :ca => Capacitor(1e-12), :cb => Capacitor(2e-12),
                :cc => Capacitor(3e-12), :gnd => Ground()],
            Any[[(:p1,1),(:la,1),(:lb,1),(:lc,1),(:ld,1)],
                [(:la,2),(:ca,1)], [(:lb,2),(:cb,1)], [(:lc,2),(:cc,1)],
                [(:ld,2),(:p1,2),(:ca,2),(:cb,2),(:cc,2),(:gnd,1)]]))

        @test refillagrees(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("C1", "1", "2", Capacitor(100e-15)), ("Lj1", "2", "0", JosephsonJunction(1e-9)), ("C2", "2", "0", Capacitor(1e-12))]))
        @test refillagrees(Circuit(
            Any[:p1 => Port(1), :la => Inductor(2e-9), :lb => Inductor(3e-9),
                :l2 => Inductor(5e-9), :jj => JosephsonJunction(1e-9),
                :c1 => Capacitor(1e-13 + 1e-16im), :gnd => Ground()],
            Any[[(:p1,1),(:la,1),(:lb,1),(:c1,1)],
                [(:la,2),(:lb,2),(:l2,1),(:jj,1)],
                [(:l2,2),(:jj,2),(:p1,2),(:c1,2),(:gnd,1)]]); Nmodes = 1)
        @test refillagrees(Circuit(
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
        @test matchesdefinition(Circuit(
            Any[:p1 => Port(1), :l1 => Inductor(1e-9), :l2 => Inductor(2e-9),
                :l3 => Inductor(4e-9), :k => MutualInductor(0.9, :l1, :l2),
                :c1 => Capacitor(1e-12), :gnd => Ground()],
            Any[[(:p1,1),(:l1,1),(:c1,1)], [(:l1,2),(:l2,1)],
                [(:l2,2),(:l3,1)],
                [(:l3,2),(:p1,2),(:c1,2),(:gnd,1)]]))

        # and a circuit with no inductance at all
        @test matchesdefinition(Circuit(
            [:p1 => Port(1), :c1 => Capacitor(1e-12), :jj => JosephsonJunction(1e-9)],
            [[(:p1,1),(:c1,1),(:jj,1)], [(:p1,2),(:c1,2),(:jj,2),Ground]]))

        # at one mode as well as many
        @test matchesdefinition(Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("L1", "1", "2", Inductor(1e-9)), ("C2", "2", "0", Capacitor(1e-12))]); Nmodes = 1)
    end
end

@testset "definitions are normalized once and may hold non numbers" begin
    JosephsonCircuits.@params La Lb
    vals = JosephsonCircuits.componentvaluestonumber(
        Any[La, La + Lb, 3.0e-12, :Lc],
        Dict(La => 1e-12, Lb => 2e-12, :Lc => 4e-12, :note => "not a value"))
    @test vals == [1e-12, 3e-12, 3e-12, 4e-12]
    d = JosephsonCircuits.normalizedefinitions(Dict("La" => 1, :Lb => 2.0im,
        :note => "text"))
    @test d == Dict(:La => 1.0 + 0im, :Lb => 2.0im)
end

@testset "the sign of the mutual stamps" begin
    JC = JosephsonCircuits
    # the second inductor declared either way round, so the graph turns one
    # of the two branches against the netlist in one case
    chain(rev) = Circuit([(:p1, 1, 0, Port(1)), (:p2, 3, 0, Port(2)),
        (:l1, 1, 2, Inductor(1e-9)),
        (:l2, (rev ? 3 : 2), (rev ? 2 : 3), Inductor(2e-9)),
        (:k, :l1, :l2, MutualInductor(0.3))])
    # the planned and the in place assembly give the mutual matrix written
    # from its definition, at the plan's own values and at new ones, with
    # the plan built once at the first: K through zero and both signs, and
    # each self inductance moved on its own, so that the mutual inductance
    # is formed from the new values of all three
    for rev in (false, true), Nmodes in (1, 4)
        cc = JC.compile(chain(rev))
        b = bound(cc)
        vvn = b.values
        plan = JC.circuitmatrixplan(cc; Nmodes = Nmodes)
        ki, l1i, l2i = (cc.componentnamedict[n] for n in ("k", "l1", "l2"))
        nm = JC.assemblematrices(plan, b)
        @test nm.Mb == coordinatemutual(cc, vvn)
        for (K, s1, s2) in ((0.3, 1.0, 1.0), (-0.3, 1.0, 1.0),
                (0.0, 1.0, 1.0), (0.7, 2.5, 0.4), (0.3, 2.0, 1.0),
                (0.3, 1.0, 0.5))
            v2 = copy(vvn)
            v2[ki] = K
            v2[l1i] = s1*vvn[l1i]
            v2[l2i] = s2*vvn[l2i]
            b2 = JC.bindvalues(cc, v2)
            ref = coordinatemutual(cc, v2)
            @test JC.assemblematrices(plan, b2).Mb == ref
            @test JC.assemblematrices!(nm, plan, b2).Mb == ref
        end
    end
end

@testset "resolved coupling and per-operation assembly storage" begin
    JC = JosephsonCircuits
    cell = Circuit([(:p,1,0,Port(2)), (:l1,1,2,Inductor(1.0)),
        (:l2,2,0,Inductor(2.0)), (:k,:l1,:l2,MutualInductor(0.2)),
        (:k2,:l1,:l2,MutualInductor(-0.2))]; pins=[1=>(:l1,1)])
    cc = compile(Circuit([:sub => cell], []))
    @test cc.couplings == [(cc.componentnamedict["sub/"*k],
        cc.componentnamedict["sub/l1"],cc.componentnamedict["sub/l2"])
        for k in ("k","k2")]
    # a port's own termination sits in the flat table between the
    # instances, so the flat indices are not the instances'
    @test cc.couplings[1][2] == 3
    b = bound(cc)
    for modes in (1,3)
        plan = JC.circuitmatrixplan(cc;Nmodes=modes)
        nm = JC.assemblematrices(plan,b)
        other = JC.assemblematrices(plan,b)
        work = JC.CircuitMatrixWorkspace(plan,nm)
        otherwork = JC.CircuitMatrixWorkspace(plan,other)
        @test work.mutualvalues !== otherwork.mutualvalues
        @test nm.Mb.nzval !== other.Mb.nzval
        @test all(iszero,nm.Mb.nzval) # opposite couplings cancel, retaining support
        for kval in (0.,0.1,1.)
            v = copy(b.values); v[cc.couplings[2][1]] = kval
            bound = JC.bindvalues(cc,v)
            got = JC.assemblematrices!(nm,plan,bound,work)
            want = JC.assemblematrices(plan,bound)
            for field in (:Cnm,:Gnm,:invLnm,:Mb,:Lb,:Ljb)
                @test getproperty(got,field) == getproperty(want,field)
            end
            @test got.Mb === nm.Mb
            @test all(iszero,other.Mb.nzval)
        end
    end
    # a complex coupling promotes the mutual inductance matrix and leaves
    # the inductances as they are
    v = Any[b.values...]; v[cc.couplings[1][1]] = 0.2+0.01im
    fresh = numericmatrices(cc, v)
    @test eltype(fresh.Lb) == Float64
    @test eltype(fresh.Mb) == ComplexF64
    # a plan is built from its compiled circuit alone: the ports in the
    # order of their numbers, whatever order the netlist gave them
    ports = Circuit([(:p2,2,0,Port(2)),(:p1,1,0,Port(1;termination=nothing))])
    pc = compile(ports); pb = bound(pc); pg = calccircuitgraph(pc)
    pp = JC.circuitmatrixplan(pc)
    @test [p.number for p in pp.ports] == [1,2]
    @test JC.assemblematrices(pp,pb).portenvironmentindices[1] == 0
end

@testset "a self loop branch and the inverse inductance" begin
    JC = JosephsonCircuits
    for modes in (1,3), value in (1.,2.,Inf)
        c = Circuit([(:p,1,0,Port(1)),(:l,1,2,Inductor(1.0)),
            (:self,2,2,Inductor(value)),(:c,2,0,Capacitor(1.0))])
        cc = compile(c)
        b = bound(cc); plan = JC.circuitmatrixplan(cc;Nmodes=modes)
        nm = JC.assemblematrices(plan,b)
        # the sparse product `Rbn' diag(1/L) Rbn` as an independent reference
        D = sparse(nm.Lbm.nzind, nm.Lbm.nzind, 1 ./ nm.Lbm.nzval, length(nm.Lbm), length(nm.Lbm))
        ref = transpose(nm.Rbnm)*D*nm.Rbnm
        @test nm.invLnm == ref
        @test nm.invLnm.colptr == ref.colptr
        @test nm.invLnm.rowval == ref.rowval
        work = JC.CircuitMatrixWorkspace(plan,nm)
        @test JC.assemblematrices!(nm,plan,b,work).invLnm == ref
    end
    JosephsonCircuits.@params L
    cc = compile(Circuit([(:l,1,0,Inductor(L))]))
    nm = numericmatrices(cc, Dict{Symbol,Any}())
    @test JC.componentvaluestonumber(nm.invLnm.nzval,Dict(:L=>2.)) == [0.5]
end
