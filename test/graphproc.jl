using JosephsonCircuits
using Test
import Graphs
import SparseArrays

@testset verbose=true "graphproc" begin
    @testset "storeuniqueloops!" begin
        lvarray = Vector{Int}[];
        @test_throws(
            ErrorException("There should only be one loop associated with each closure branch."),
            JosephsonCircuits.storeuniqueloops!(lvarray,[1, 2, 3],[[1,2,3],[4,5,6]])
        )
    end

    @testset "calcgraphs" begin
        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(2, 1), (2, 1), (2, 1), (3, 1)], 3),
            JosephsonCircuits.CircuitGraph(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 2), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]]), 2),
            )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcgraphs([(4, 3), (3, 6), (5, 3), (3, 7), (2, 4), (6, 8), (2, 5), (8, 7), (2, 8)], 8),
            JosephsonCircuits.CircuitGraph(Dict((6, 8) => 8, (7, 8) => 9, (2, 5) => 2, (3, 6) => 6, (8, 6) => 8, (5, 2) => 2, (2, 8) => 3, (6, 3) => 6, (3, 5) => 5, (3, 4) => 4, (5, 3) => 5, (3, 7) => 7, (8, 7) => 9, (2, 4) => 1, (4, 3) => 4, (8, 2) => 3, (7, 3) => 7, (4, 2) => 1), SparseArrays.sparse([1, 2, 3, 4, 5, 6, 7, 1, 4, 2, 5, 6, 8, 7, 9, 3, 8, 9], [1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7], [-1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1], 9, 7), [(2, 4), (2, 5), (2, 8), (3, 4), (3, 6), (3, 7)], [(5, 3), (8, 6), (8, 7)], [(2, 4), (2, 5), (2, 8), (3, 4), (3, 5), (3, 6), (3, 7), (6, 8), (7, 8)], [[2, 4, 3, 5], [2, 4, 3, 6, 8], [2, 4, 3, 7, 8]], [1], Graphs.SimpleGraphs.SimpleGraph{Int64}(9, [Int64[], [4, 5, 8], [4, 5, 6, 7], [2, 3], [2, 3], [3, 8], [3, 8], [2, 6, 7]]), 9),
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
end