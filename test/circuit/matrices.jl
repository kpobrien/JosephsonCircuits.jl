using JosephsonCircuits
using Test

@testset verbose=true "the circuit matrices" begin

    @testset "calcMb JJ as first inductor" begin
        Nmodes = 2
        Nbranches = 2
        componenttypes = [:Lj,:K,:L,:C]
        nodeindices = [2 0 3 3; 1 0 1 1]
        componentvalues = [1.0e-9, 0.1, 4.0e-9, 2.0e-12]
        componentnamedict = Dict{Symbol, Int}(:C1 => 4,:L2 => 3,:Lj1 => 1,:K1 => 2)
        edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
        mutualinductorbranchnames = [ :Lj1, :L2]
        Rbn = JosephsonCircuits.SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 2)

        @test_throws(
            ArgumentError("Mutual coupling coefficient K must couple two inductors. Lj1 is not an inductor."),
            JosephsonCircuits.calcMb(componenttypes,nodeindices,componentvalues,componentnamedict,mutualinductorbranchnames,edge2indexdict,Rbn,Nmodes,Nbranches)
        )
    end

    @testset "calcMb JJ as second inductor" begin
        Nmodes = 2
        Nbranches = 2
        componenttypes = [:L,:K,:Lj,:C]
        nodeindices = [2 0 3 3; 1 0 1 1]
        componentvalues = [1.0e-9, 0.1, 4.0e-9, 2.0e-12]
        componentnamedict = Dict{Symbol, Int}(:C1 => 4,:Lj2 => 3,:L1 => 1,:K1 => 2)
        edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
        mutualinductorbranchnames = [ :L1, :Lj2]
        Rbn = JosephsonCircuits.SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 2)

        @test_throws(
            ArgumentError("Mutual coupling coefficient K must couple two inductors. Lj2 is not an inductor."),
            JosephsonCircuits.calcMb(componenttypes,nodeindices,componentvalues,componentnamedict,mutualinductorbranchnames,edge2indexdict,Rbn,Nmodes,Nbranches)
        )
    end

    @testset "calcLmean_inner" begin
        @test_throws(
            DimensionMismatch("componenttypes and componentvalues should have the same length"),
            JosephsonCircuits.calcLmean_inner([:L,:C,:Lj],[10,4,5,1],Float64[])
        )
    end

    @testset "calcnodematrix" begin
        @test_throws(
            DimensionMismatch("nodeindices should have a first dimension size of 2."),
            JosephsonCircuits.calcnodematrix(
                [:R,:R],[2 3;1 1;0 0],[1.0,2.0],Float64[],1,3,:R,false)
        )
        @test_throws(
            DimensionMismatch("componenttypes, nodeindices, and componentvalues should have the same length"),
            JosephsonCircuits.calcnodematrix([:R],[2 3;1 1],[1.0,2.0],
                Float64[],1,3,:R,false)
        )
    end

    @testset "combine" begin

        a = rand()
        b = rand()
        @test(JosephsonCircuits.combine_sum(a,b) == a+b)
        @test(JosephsonCircuits.combine_reciprocal_sum(a,b) == a*b/(a+b))
        @test_throws(
            ArgumentError("Components 1 and 2 cannot be combined to a single element. Please place the two components between different nodes."),
            JosephsonCircuits.combine_error(1,2),
        )
    end

    # the element type the matrix builders assemble in
    @testset "calcvaluetype" begin
        @test_throws(
            DimensionMismatch("componenttypes and componentvalues should have the same length"),
            JosephsonCircuits.calcvaluetype(
                [:C,:R],
                [1,2,3],
                [:R]
            )
        )
    end

    @testset "the sign of a mutual inductance does not follow the node names" begin
        # two inductors in series sharing a node and coupled to each other:
        # walking the chain the series inductance is L1 + L2 + 2M, whatever
        # the nodes are called.  The incidence matrix orients a branch by
        # the spanning tree, so the netlist's own order has to be carried
        # over to it.
        L = 1e-9
        K = 0.5
        Z0 = 50.0
        w = 2pi*1e9
        chain(a, b, c) = Circuit([
            (:p1, a, 0, Port(1; Z0 = Z0)),
            (:p2, c, 0, Port(2; Z0 = Z0)),
            (:l1, a, b, Inductor(L)),
            (:l2, b, c, Inductor(L)),
            (:k, :l1, :l2, MutualInductor(K))])
        function seriesinductance(nodes)
            S21 = hblinsolve([w], chain(nodes...)).S((0,), 2, (0,), 1, 1)
            return 2*Z0*imag(1/S21 - 1)/w
        end
        for nodes in ((1, 2, 3), (1, 3, 2), (2, 1, 3), (3, 2, 1), (3, 1, 2))
            @test isapprox(seriesinductance(nodes), 2*L*(1 + K); rtol = 1e-10)
        end
        # and reversing the terminals of one of them opposes the currents
        opposed = Circuit([
            (:p1, 1, 0, Port(1; Z0 = Z0)),
            (:p2, 3, 0, Port(2; Z0 = Z0)),
            (:l1, 1, 2, Inductor(L)),
            (:l2, 3, 2, Inductor(L)),
            (:k, :l1, :l2, MutualInductor(K))])
        S21 = hblinsolve([w], opposed).S((0,), 2, (0,), 1, 1)
        @test isapprox(2*Z0*imag(1/S21 - 1)/w, 2*L*(1 - K); rtol = 1e-10)
    end

    @testset "the mutual orientation cache" begin
        JC = JosephsonCircuits
        Z0 = 50.0
        # the second inductor declared either way round, so the graph turns
        # one of the two branches against the netlist in one case
        chain(rev) = Circuit([
            (:p1, 1, 0, Port(1; Z0 = Z0)),
            (:p2, 3, 0, Port(2; Z0 = Z0)),
            (:l1, 1, 2, Inductor(1e-9)),
            (:l2, (rev ? 3 : 2), (rev ? 2 : 3), Inductor(2e-9)),
            (:k, :l1, :l2, MutualInductor(0.3))])

        function pieces(c, Nmodes)
            cc = JC.compile(c)
            cg = calccircuitgraph(cc)
            b = JC.bind(cc)
            vvn = JC.componentvaluestonumber(cc.componentvalues,
                Dict{Any,Any}())
            return cc, cg, b, vvn, JC.circuitmatrixplan(cc, cg, b;
                Nmodes = Nmodes)
        end

        # one walk of the incidence matrix names each branch's endpoints the
        # way the matrix itself does
        let (cc, cg, b, vvn, plan) = pieces(chain(false), 1)
            from, to = JC.branchendpoints(cg.Rbn, cg.Nbranches)
            for bi in 1:cg.Nbranches
                from[bi] > 1 && @test cg.Rbn[bi, from[bi]-1] == -1
                to[bi] > 1 && @test cg.Rbn[bi, to[bi]-1] == 1
                @test from[bi] != to[bi]
            end
        end

        # reversing the declared terminals flips the cached product, and it
        # is the sign the matrix is actually filled with
        let (_, _, _, _, p1) = pieces(chain(false), 1),
            (_, _, _, _, p2) = pieces(chain(true), 1)
            @test length(p1.mutualorientations) == 1
            @test p1.mutualorientations[1] == -p2.mutualorientations[1]
        end

        # direct, planned and in place assembly agree on the mutual matrix
        # at the plan's own values and at new ones, with the plan built once
        # at the first: K through zero and both signs, and each self
        # inductance moved on its own, so that the mutual inductance is
        # formed from the new values of all three
        for rev in (false, true), Nmodes in (1, 4)
            cc, cg, b, vvn, plan = pieces(chain(rev), Nmodes)
            key(n) = keytype(cc.componentnamedict) === Symbol ?
                Symbol(n) : string(n)
            ki = cc.componentnamedict[key("k")]
            l1i = cc.componentnamedict[key("l1")]
            l2i = cc.componentnamedict[key("l2")]
            nm = JC.assemblematrices(plan, b)
            @test nm.Mb == numericmatrices(cc, cg, vvn; Nmodes = Nmodes).Mb
            for (K, s1, s2) in ((0.3, 1.0, 1.0), (-0.3, 1.0, 1.0),
                    (0.0, 1.0, 1.0), (0.7, 2.5, 0.4), (0.3, 2.0, 1.0),
                    (0.3, 1.0, 0.5))
                v2 = copy(vvn)
                v2[ki] = K
                v2[l1i] = s1*vvn[l1i]
                v2[l2i] = s2*vvn[l2i]
                b2 = JC.bindvalues(cc, v2)
                ref = numericmatrices(cc, cg, v2; Nmodes = Nmodes).Mb
                @test JC.assemblematrices(plan, b2).Mb == ref
                @test JC.assemblematrices!(nm, plan, b2).Mb == ref
            end
        end
    end
end
