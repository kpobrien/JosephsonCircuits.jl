using Random
using JosephsonCircuits
using Test
import Graphs
import SparseArrays

@testset verbose=true "graphproc" begin
    @testset "calcgraphs" begin
        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(2, 1), (2, 1), (2, 1), (3, 1)], 3; loops = true),
            JosephsonCircuits.CircuitGraph(JosephsonCircuits.CircuitTopology(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 2), 2), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]])),
            )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(4, 3), (3, 6), (5, 3), (3, 7), (2, 4), (6, 8), (2, 5), (8, 7), (2, 8)], 8; loops = true),
            JosephsonCircuits.CircuitGraph(JosephsonCircuits.CircuitTopology(Dict((6, 8) => 8, (7, 8) => 9, (2, 5) => 2, (3, 6) => 6, (8, 6) => 8, (5, 2) => 2, (2, 8) => 3, (6, 3) => 6, (3, 5) => 5, (3, 4) => 4, (5, 3) => 5, (3, 7) => 7, (8, 7) => 9, (2, 4) => 1, (4, 3) => 4, (8, 2) => 3, (7, 3) => 7, (4, 2) => 1), SparseArrays.sparse([1, 2, 3, 4, 5, 6, 7, 1, 4, 2, 5, 6, 8, 7, 9, 3, 8, 9], [1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7], [-1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1], 9, 7), 9), [(2, 4), (2, 5), (2, 8), (3, 4), (3, 6), (3, 7)], [(5, 3), (8, 6), (8, 7)], [(2, 4), (2, 5), (2, 8), (3, 4), (3, 5), (3, 6), (3, 7), (6, 8), (7, 8)], [[3, 4, 2, 5], [6, 3, 4, 2, 8], [7, 3, 4, 2, 8]], [1], Graphs.SimpleGraphs.SimpleGraph{Int64}(9, [Int64[], [4, 5, 8], [4, 5, 6, 7], [2, 3], [2, 3], [3, 8], [3, 8], [2, 6, 7]])),
            )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(2, 1), (2, 1), (3, 1)], 4; loops = true),
            JosephsonCircuits.CircuitGraph(JosephsonCircuits.CircuitTopology(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 3), 2), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]])),
            )
    end

    # the branches the incidence matrix is built from
    @testset "extractbranches" begin
        @test_throws(
            DimensionMismatch("componenttypes must have the same length as the number of node indices"),
            JosephsonCircuits.extractbranches(
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3; 1 1 1 3 1]
            )
        )

        @test_throws(
            DimensionMismatch("the length of the first axis must be 2"),
            JosephsonCircuits.extractbranches(
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3 3; 1 1 1 3 1 1; 0 0 0 0 0 0],
            )
        )
    end
end

# The loop of a closure branch is the unique path through the spanning tree
# between its endpoints, whatever its length: a search bounded in the number
# of edges would report a longer loop as no loop at all, silently, since the
# empty entry is also what a branch with no loop produces.
@testset verbose=true "fundamental loops of any length" begin
    JC = JosephsonCircuits

    # a ring of N inductors is one loop through all N of them
    for N in (3, 5, 10, 11, 12, 40)
        edges = [(i, i == N ? 1 : i+1) for i in 1:N]
        cg = JC.calcgraphs(edges, N; loops = true)
        @test length(cg.lvarray) == 1
        loop = only(cg.lvarray)
        @test length(loop) == N
        @test sort(loop) == collect(1:N)     # every node, once
        # and it is a walk: consecutive vertices are joined by an edge, as
        # are the last and the first
        joined(a, b) = (a, b) in edges || (b, a) in edges
        @test all(joined(loop[i], loop[i+1]) for i in 1:N-1)
        @test joined(loop[end], loop[1])
    end

    # two inductors between one pair of nodes are one edge of a simple
    # graph, so there is no closure branch and no loop to report; the two
    # vertex guard in the walk is there for the case where a closure branch
    # and a tree edge join the same pair for some other reason
    @test isempty(JC.calcgraphs([(1,2), (1,2)], 2; loops = true).lvarray)

    # a tree has no closure branches and so no loops
    @test isempty(JC.calcgraphs([(1,2), (2,3), (3,4)], 4; loops = true).lvarray)

    # the walk itself, on a rooted tree: 3 -> 2 -> 1 -> 4
    parent, depth = JC.rootedtree(
        JC.Graphs.SimpleGraphFromIterator(JC.tuple2edge([(1,2),(2,3),(1,4)])))
    @test JC.treepath(parent, depth, 3, 4) == [3, 2, 1, 4]
    @test JC.treepath(parent, depth, 4, 3) == [4, 1, 2, 3]
    @test JC.treepath(parent, depth, 2, 2) == [2]
end

@testset "the incidence matrix of the branches" begin
    # written from its definition, without the graph library: the distinct
    # branches in ascending order of their endpoints, each row -1 at its
    # lower node and +1 at its higher, ground (node 1) having no column
    # and a self loop no entries, and either orientation of a branch naming
    # its row
    rng = Random.MersenneTwister(190926)
    JC = JosephsonCircuits
    for trial in 1:20
        n = rand(rng, 1:40)
        edges = [(rand(rng,1:n), rand(rng,1:n)) for _ in 1:rand(rng,0:80)]
        t = JC.circuittopology(edges, n)
        branches = sort(unique([minmax(e...) for e in edges]))
        Rbn = zeros(Int, length(branches), n - 1)
        for (i, (a, b)) in enumerate(branches)
            a == b && continue
            a > 1 && (Rbn[i, a-1] = -1)
            b > 1 && (Rbn[i, b-1] = 1)
        end
        row(e) = findfirst(==(minmax(e...)), branches)
        @test t.Nbranches == length(branches) && Matrix(t.Rbn) == Rbn &&
            all(t.edge2indexdict[e] == t.edge2indexdict[reverse(e)] == row(e)
                for e in edges)
    end
end

@testset "components with both terminals on one node" begin
    # such a component carries no current: the response of the circuit, its
    # noise included, is the one without it, and a port so placed is refused
    ws = 2pi*(4.5:0.5:5.5)*1e9
    sol(c) = hblinsolve(ws, Circuit(c); keyedarrays = false, returnSnoise = true)
    # node 2 has no inductive branch, so a self loop there is its only one
    base = [(:p1, 1, 0, Port(1)), (:c1, 1, 2, Capacitor(100e-15)),
        (:c2, 2, 0, Capacitor(1e-12)), (:r2, 2, 0, Resistor(1e4))]
    ref = sol(base)
    # and on node 5, which nothing else touches
    for extra in ((:cx, 0, 0, Capacitor(1e-12)), (:rx, 0, 0, Resistor(50.0)),
            (:rx, 2, 2, Resistor(50.0)), (:lx, 0, 0, Inductor(1e-9)),
            (:lx, 2, 2, Inductor(1e-9)), (:jx, 2, 2, JosephsonJunction(1e-9)),
            (:cx, 5, 5, Capacitor(1e-12)), (:lx, 5, 5, Inductor(1e-9)),
            (:jx, 5, 5, JosephsonJunction(1e-9)))
        s = sol(vcat(base, [extra]))
        @test s.S ≈ ref.S && s.QE ≈ ref.QE && size(s.Snoise) == size(ref.Snoise)
    end
    @test_throws ArgumentError sol([(:p1, 1, 1, Port(1)), (:c1, 1, 0, Capacitor(1e-12))])
    # on a node of its own, a secondary coupled to a primary shorts it, which
    # leaves the primary L1(1 - K^2), and a line whose far port is shorted is
    # a stub of input impedance i Z0 tan(w l/vp)
    L1, L2, K = 1e-9, 2e-9, 0.6
    shorted = sol([(:p1, 1, 0, Port(1)), (:l1, 1, 0, Inductor(L1)),
        (:l2, 3, 3, Inductor(L2)), (:k, :l1, :l2, MutualInductor(K)),
        (:c1, 1, 0, Capacitor(1e-12))])
    reduced = sol([(:p1, 1, 0, Port(1)), (:l1, 1, 0, Inductor(L1*(1 - K^2))),
        (:c1, 1, 0, Capacitor(1e-12))])
    @test shorted.S ≈ reduced.S
    Z0, len = 50.0, 0.004
    stub = sol([(:p1, 1, 0, Port(1)),
        (:line, 1, 0, 2, 2, TransmissionLine(Z0, len; grounded = false))])
    zin = [im*Z0*tan(w*len/JosephsonCircuits.speed_of_light) for w in ws]
    @test vec(stub.S) ≈ (zin .- 50)./(zin .+ 50)
end

@testset "a subnetwork no element connects to ground" begin
    # harmonic balance cannot determine its potential, so both solvers
    # refuse it by name; the transient fixes that potential and solves it
    c = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(1e-12)),
        (:c2, 2, 3, Capacitor(1e-12)), (:l2, 2, 3, Inductor(1e-9))])
    psc = compile(c)
    @test [psc.nodenames[i] for i in only(JosephsonCircuits.isolatedsubnetworks(psc))] ==
        ["2", "3"]
    @test_throws ArgumentError hblinsolve(2pi*[4e9, 5e9], c)
    @test_throws ArgumentError hbnlsolve((2pi*5e9,), (4,),
        [(mode = (1,), port = 1, current = 1e-8)], c)
    s = transientsolve(transientproblem(c), (0.0, 1e-10); dt = 1e-12)
    @test all(isfinite, s.outgoing)
end
