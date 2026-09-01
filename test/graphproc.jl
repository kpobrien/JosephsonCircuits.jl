using JosephsonCircuits
using Test
import Graphs
import SparseArrays

@testset verbose=true "graphproc" begin
    @testset "calcgraphs" begin
        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(2, 1), (2, 1), (2, 1), (3, 1)], 3),
            JosephsonCircuits.CircuitGraph(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 2), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]]), 2),
            )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(4, 3), (3, 6), (5, 3), (3, 7), (2, 4), (6, 8), (2, 5), (8, 7), (2, 8)], 8),
            JosephsonCircuits.CircuitGraph(Dict((6, 8) => 8, (7, 8) => 9, (2, 5) => 2, (3, 6) => 6, (8, 6) => 8, (5, 2) => 2, (2, 8) => 3, (6, 3) => 6, (3, 5) => 5, (3, 4) => 4, (5, 3) => 5, (3, 7) => 7, (8, 7) => 9, (2, 4) => 1, (4, 3) => 4, (8, 2) => 3, (7, 3) => 7, (4, 2) => 1), SparseArrays.sparse([1, 2, 3, 4, 5, 6, 7, 1, 4, 2, 5, 6, 8, 7, 9, 3, 8, 9], [1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7], [-1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1], 9, 7), [(2, 4), (2, 5), (2, 8), (3, 4), (3, 6), (3, 7)], [(5, 3), (8, 6), (8, 7)], [(2, 4), (2, 5), (2, 8), (3, 4), (3, 5), (3, 6), (3, 7), (6, 8), (7, 8)], [[3, 4, 2, 5], [6, 3, 4, 2, 8], [7, 3, 4, 2, 8]], [1], Graphs.SimpleGraphs.SimpleGraph{Int64}(9, [Int64[], [4, 5, 8], [4, 5, 6, 7], [2, 3], [2, 3], [3, 8], [3, 8], [2, 6, 7]]), 9),
            )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(2, 1), (2, 1), (3, 1)],4),
            JosephsonCircuits.CircuitGraph(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 3), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]]), 2),
            )
    end

    @testset "edge2index" begin
        @test isequal(
            JosephsonCircuits.edge2index(JosephsonCircuits.Graphs.path_digraph(4)),
            Dict((1, 2) => 1, (2, 1) => 1, (3, 2) => 2, (3, 4) => 3, (2, 3) => 2, (4, 3) => 3),
            )
    end

    @testset "tuple2edge" begin

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Int}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2)),
            Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Int64}(Graphs.Edge(3 => 4) => 3, Graphs.Edge(2 => 3) => 2, Graphs.Edge(1 => 2) => 1),
            )

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Float64}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2)),
            Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Float64}(Graphs.Edge(3 => 4) => 3.0, Graphs.Edge(2 => 3) => 2.0, Graphs.Edge(1 => 2) => 1.0),
            )

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Complex{Float64}}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2)),
            Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, ComplexF64}(Graphs.Edge(3 => 4) => 3.0 + 0.0im, Graphs.Edge(2 => 3) => 2.0 + 0.0im, Graphs.Edge(1 => 2) => 1.0 + 0.0im),
            )

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Any}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2)),
             Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Any}(Graphs.Edge(3 => 4) => 3, Graphs.Edge(2 => 3) => 2, Graphs.Edge(1 => 2) => 1),
            )

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Int}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3)),
            Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, Int64}((Graphs.Edge(5 => 6), Graphs.Edge(7 => 8),) => 3, (Graphs.Edge(1 => 2), Graphs.Edge(3 => 4),) => 1),
            )

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Float64}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3)),
            Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, Float64}((Graphs.Edge(5 => 6), Graphs.Edge(7 => 8),) => 3.0, (Graphs.Edge(1 => 2), Graphs.Edge(3 => 4),) => 1.0),
            )

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Complex{Float64}}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3)),
            Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, ComplexF64}((Graphs.Edge(5 => 6), Graphs.Edge(7 => 8),) => 3.0 + 0.0im, (Graphs.Edge(1 => 2), Graphs.Edge(3 => 4),) => 1.0 + 0.0im),
            )

        @test isequal(
            JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Any}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3)),
            Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, Any}((Graphs.Edge(5 => 6), Graphs.Edge(7 => 8),) => 3, (Graphs.Edge(1 => 2), Graphs.Edge(3 => 4),) => 1),
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

    @testset "extractbranches!" begin
        @test_throws(
            DimensionMismatch("branchvector should be length zero"),
            JosephsonCircuits.extractbranches!(
                [1],
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3 3; 1 1 1 3 1 1],
            )
        )
    end
end

# The loop of a closure branch is the unique path through the spanning tree
# between its endpoints. It used to be found by adding the branch back and
# enumerating cycles of at most ten edges, so a longer loop was reported as
# no loop at all -- silently, since the empty entry is also what a branch
# with no loop produces.
@testset verbose=true "fundamental loops of any length" begin
    JC = JosephsonCircuits

    # a ring of N inductors is one loop through all N of them
    for N in (3, 5, 10, 11, 12, 40)
        edges = [(i, i == N ? 1 : i+1) for i in 1:N]
        cg = JC.calcgraphs(edges, N)
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
    # graph, so there is no closure branch and no loop to report. That is
    # what the graph did before this too, and the two vertex guard in the
    # walk is there for the case where a closure branch and a tree edge
    # join the same pair for some other reason
    @test isempty(JC.calcgraphs([(1,2), (1,2)], 2).lvarray)

    # a tree has no closure branches and so no loops
    @test isempty(JC.calcgraphs([(1,2), (2,3), (3,4)], 4).lvarray)

    # the walk itself, on a rooted tree: 3 -> 2 -> 1 -> 4
    parent, depth = JC.rootedtree(
        JC.Graphs.SimpleGraphFromIterator(JC.tuple2edge([(1,2),(2,3),(1,4)])))
    @test JC.treepath(parent, depth, 3, 4) == [3, 2, 1, 4]
    @test JC.treepath(parent, depth, 4, 3) == [4, 1, 2, 3]
    @test JC.treepath(parent, depth, 2, 2) == [2]
end
