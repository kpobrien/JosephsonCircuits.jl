
"""
    CircuitGraph(edge2indexdict,Rbn,searray,cearray,glearray,lvarray,
        isolatednodes,gl,Nbranches)

A simple structure to hold the circuit graph information.
"""
struct CircuitGraph
    edge2indexdict
    Rbn
    searray
    cearray
    glearray
    lvarray
    isolatednodes
    gl
    Nbranches
end

"""
    calccircuitgraph(parsedsortedcircuit::ParsedSortedCircuit)

Calculate the superconducting spanning tree, incidence matrix, closure branches,
and loops from the parsed and sorted circuit.

See also [`CircuitGraph`](@ref), [`calcgraphs`](@ref), and [`extractbranches`](@ref) 
for more explanation.

# Examples
```jldoctest
@variables Ipump Rleft L1 K1 L2 C2
psc = JosephsonCircuits.ParsedSortedCircuit(
    [2 2 2 2 0 3 3; 1 1 1 1 0 1 1],
    ["0", "1", "2"],
    ["L1", "L2"],
    ["P1", "I1", "R1", "L1", "K1", "L2", "C2"],
    [:P, :I, :R, :L, :K, :L, :C],
    Num[1, Ipump, Rleft, L1, K1, L2, C2],
    Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5),
    3)
cg = JosephsonCircuits.calccircuitgraph(psc)
# output
JosephsonCircuits.CircuitGraph(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), sparse([1, 2], [1, 2], [1, 1], 2, 2), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]]), 2)
```
```jldoctest
@variables Ipump Rleft L Lj Cj
circuit = Tuple{String,String,String,Num}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","2",L))
push!(circuit,("Lj1","2","0",Lj))
push!(circuit,("C2","2","0",Cj))
psc = JosephsonCircuits.parsesortcircuit(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
# output
JosephsonCircuits.CircuitGraph(Dict((3, 2) => 3, (1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1, (2, 3) => 3), sparse([1, 3, 2, 3], [1, 1, 2, 2], [1, -1, 1, 1], 3, 2), [(1, 2), (1, 3)], [(3, 2)], [(1, 2), (1, 3), (2, 3)], [[1, 2, 3]], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(3, [[2, 3], [1, 3], [1, 2]]), 3)
```
"""
function calccircuitgraph(parsedsortedcircuit::ParsedSortedCircuit)

    branchvector = extractbranches(parsedsortedcircuit.typevector,
                                parsedsortedcircuit.nodeindexarraysorted)

    # calculate the graph of inductive components glelist, the
    # superconducting spanning tree selist, and the list of loop
    # indices celist. 
    return calcgraphs(branchvector,parsedsortedcircuit.Nnodes)

end


"""
    calcgraphs(gl,Nnodes)

Calculate the superconducting spanning tree, closure branches, and loops.
Accepts the graph of linear inductors and Josephson junctions. Outputs lists
of edges that can be used to generate graphs.

# Examples
```jldoctest
julia> JosephsonCircuits.calcgraphs([(2, 1), (2, 1), (2, 1), (3, 1)], 3)
JosephsonCircuits.CircuitGraph(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), sparse([1, 2], [1, 2], [1, 1], 2, 2), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]]), 2)

julia> JosephsonCircuits.calcgraphs([(4, 3), (3, 6), (5, 3), (3, 7), (2, 4), (6, 8), (2, 5), (8, 7), (2, 8)], 8)
JosephsonCircuits.CircuitGraph(Dict((6, 8) => 8, (2, 5) => 2, (3, 7) => 7, (6, 3) => 6, (7, 8) => 9, (3, 4) => 4, (7, 3) => 7, (2, 8) => 3, (4, 2) => 1, (8, 6) => 8…), sparse([1, 2, 3, 4, 5, 6, 7, 1, 4, 2, 5, 6, 8, 7, 9, 3, 8, 9], [1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 7], [-1, -1, -1, -1, -1, -1, -1, 1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1], 9, 7), [(2, 4), (2, 5), (2, 8), (3, 4), (3, 6), (3, 7)], [(5, 3), (8, 6), (8, 7)], [(2, 4), (2, 5), (2, 8), (3, 4), (3, 5), (3, 6), (3, 7), (6, 8), (7, 8)], [[2, 4, 3, 5], [2, 4, 3, 6, 8], [2, 4, 3, 7, 8]], [1], Graphs.SimpleGraphs.SimpleGraph{Int64}(9, [Int64[], [4, 5, 8], [4, 5, 6, 7], [2, 3], [2, 3], [3, 8], [3, 8], [2, 6, 7]]), 9)
```
"""
function calcgraphs(Ledgearray::Array{Tuple{Int, Int}, 1},Nnodes::Int)

    gl = Graphs.SimpleGraphFromIterator(tuple2edge(Ledgearray))

    searray = Vector{Tuple{Int, Int}}(undef, 0)
    cearray = Vector{Tuple{Int, Int}}(undef, 0)
    lvarray = Vector{Vector{Int}}(undef, 0)
    glearray = Vector{Tuple{Int, Int}}(undef, 0)
    isolatednodes = Vector{Int}(undef,0)

    # break into groups of connected components.
    # loop over these and construct the spanning trees.
    for v in Graphs.connected_components(gl)

        # list of vertices has length one, this node is isolated and
        # should be removed later
        if length(v) == 1
            push!(isolatednodes,v[1])
        end

        #generate the subgraph
        gli, vmap = Graphs.induced_subgraph(gl,v)

        #superconducting spanning tree
        si = Graphs.SimpleGraph(Graphs.kruskal_mst(gli))
        #si = SimpleGraph(prim_mst(gli))

        #calculate the closure branches
        ci = collect(Graphs.edges(Graphs.difference(gli,si)))
    
        #find the loop indices associated with each closure branch by starting
        #with the superconducting spanning tree (which has no loops), then
        #adding a closure branch and looking for the cycle this creates. 
        for cj in ci
            #push!(cearray,Edge(vmap[dst(cj)],vmap[src(cj)]))
            push!(cearray,(vmap[Graphs.dst(cj)],vmap[Graphs.src(cj)]))

            stmp = copy(si)
    #         stmp = SimpleGraphFromIterator(edges(si))
            Graphs.add_edge!(stmp, Graphs.src(cj), Graphs.dst(cj))
            #l = simplecycles_limited_length(stmp,nv(gl))
            l = Graphs.simplecycles_limited_length(stmp,10)
            ul = unique(x->sort(x),l[length.(l).>2])
            storeuniqueloops!(lvarray,vmap,ul)
        end
        
        #create a directed version of the superconducting spanning tree
        #starting vertex for the tree will always be 1
        if Graphs.ne(si) == 0
            sid = si
        else
            sid = Graphs.SimpleGraph(Graphs.bfs_tree(si,1))
        end
        
        # add the closure branches back to get a directed version of the
        # superconducting graph gl
        glid = copy(sid)
        for cj in ci
            Graphs.add_edge!(glid, Graphs.dst(cj), Graphs.src(cj))
        end
        
        # offset the vertices of the superconducting spanning tree by
        # the vmap.
        for e in Graphs.edges(sid)
            #push!(searray,Edge(vmap[src(e)],vmap[dst(e)]))
            push!(searray,(vmap[Graphs.src(e)],vmap[Graphs.dst(e)]))

        end
        for e in Graphs.edges(glid)
            #push!(glearray,Edge(vmap[src(e)],vmap[dst(e)]))
            push!(glearray,(vmap[Graphs.src(e)],vmap[Graphs.dst(e)]))

        end
    end

    gl2 = Graphs.SimpleDiGraphFromIterator(tuple2edge(glearray))
    # if more vertices were found when parsing all the compoments than there are in
    # the graph of inductive components (gl), then add vertices.this can happen if
    # there are only non-inductive components connected the last nodes in the graph. 
    if Graphs.nv(gl2) < Nnodes
        Graphs.add_vertices!(gl2,Nnodes-Graphs.nv(gl))
    end

    # create a dictionary that maps index to indices
    edge2indexdict = edge2index(gl2)

    # convert the branch inductance matrices 
    # to inverse node inductance matrices. 
    # get rid of the first node (the node to ground)
    Rbn=sparse(transpose(Graphs.incidence_matrix(gl2,oriented=true)))[:,2:end]

    Nbranches = Graphs.ne(gl2)

    return CircuitGraph(edge2indexdict,Rbn,searray,cearray,glearray,lvarray,
        isolatednodes,gl,Nbranches)
end


"""
    storeuniqueloops!(lvarray,vmap,ul)

# Examples
```jldoctest
julia> lvarray = Vector{Int}[];JosephsonCircuits.storeuniqueloops!(lvarray,[1, 2, 3],[[1,2,3]]);lvarray
1-element Vector{Vector{Int64}}:
 [1, 2, 3]

julia> lvarray = Vector{Int}[];JosephsonCircuits.storeuniqueloops!(lvarray,[1, 2, 3],[[1,2,3],[4,5,6]]);lvarray
ERROR: There should only be one loop associated with each closure branch.
```
"""
function storeuniqueloops!(lvarray,vmap,ul)

    if length(ul) > 1
        error("There should only be one loop associated with each closure branch.")
    elseif length(ul) < 1
        # println("Warning: Loop exists but max loop size too small to find vertices.")
        push!(lvarray,Int[])
    else
        # println("closure branch: ",cj,", loop vertices:", first(ul))
        push!(lvarray,vmap[first(ul)])
    end
    return nothing
end

"""
    edge2index(graph::SimpleDiGraph{Int})

Generate a dictionary where the tuple of nodes defining an edge of a graph
is the key and the value is an index. The index gives the order the edge
is found when iterating over the edges of the graph. The same index is used
for both orderings of source and destination nodes on the edge. We don't care
about the ordering of the indices as long as they are sequential and unique.

# Examples
```jldoctest
julia> JosephsonCircuits.edge2index(JosephsonCircuits.Graphs.path_digraph(4))
Dict{Tuple{Int64, Int64}, Int64} with 6 entries:
  (3, 2) => 2
  (1, 2) => 1
  (2, 1) => 1
  (3, 4) => 3
  (4, 3) => 3
  (2, 3) => 2
```
"""
function edge2index(graph::Graphs.SimpleDiGraph{Int})
    edge2indexdict = Dict{Tuple{Int, Int},Int}()
    for (i,e) in enumerate(Graphs.edges(graph))
        edge2indexdict[(Graphs.src(e),Graphs.dst(e))] = i
        edge2indexdict[(Graphs.dst(e),Graphs.src(e))] = i
    end
    return edge2indexdict
end


"""
    tuple2edge(tuplevector::Vector{Tuple{Int, Int}})

Converts a vector of edges specified with tuples of integers to a vector of
Graphs edges.

# Examples
```jldoctest
julia> JosephsonCircuits.tuple2edge([(1,2),(3,4)])
2-element Vector{Graphs.SimpleGraphs.SimpleEdge{Int64}}:
 Edge 1 => 2
 Edge 3 => 4
```
"""
function tuple2edge(tuplevector::Vector{Tuple{Int, Int}})
    edgevector = Vector{Graphs.SimpleGraphs.SimpleEdge{Int}}(undef, 0)

    for i in 1:length(tuplevector)
        push!(edgevector,Graphs.Edge(tuplevector[i][1],tuplevector[i][2]))
    end
    return edgevector
end

"""
    tuple2edge(tuplevector::Vector{Tuple{Int, Int, Int, Int}})

Converts a vector of edges specified with tuples of integers to a vector of
Graphs edges.

# Examples
```jldoctest
julia> JosephsonCircuits.tuple2edge([(1,2,3,4),(5,6,7,8)])
2-element Vector{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}}:
 (Edge 1 => 2, Edge 3 => 4)
 (Edge 5 => 6, Edge 7 => 8)
```
"""
function tuple2edge(tuplevector::Vector{Tuple{Int, Int, Int, Int}})
    edgevector = Vector{
        Tuple{
            Graphs.SimpleGraphs.SimpleEdge{Int},
            Graphs.SimpleGraphs.SimpleEdge{Int}
        }
    }(undef, 0)

    for i in 1:length(tuplevector)
        push!(
            edgevector,
            (
                Graphs.Edge(tuplevector[i][1],tuplevector[i][2]),
                Graphs.Edge(tuplevector[i][3],tuplevector[i][4])
            )
        )
    end

    return edgevector
end

"""
    tuple2edge(tupledict)

Converts a dictionary whose keys are edges specified by tuples of integers to
a dictionary whose keys are Graphs edges. The values associated with each key
are preserved.

# Examples
```jldoctest
julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Int}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2))
Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Int64} with 3 entries:
  Edge 1 => 2 => 1
  Edge 3 => 4 => 3
  Edge 2 => 3 => 2

julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Float64}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2))
Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Float64} with 3 entries:
  Edge 1 => 2 => 1.0
  Edge 3 => 4 => 3.0
  Edge 2 => 3 => 2.0

julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Complex{Float64}}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2))
Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, ComplexF64} with 3 entries:
  Edge 1 => 2 => 1.0+0.0im
  Edge 3 => 4 => 3.0+0.0im
  Edge 2 => 3 => 2.0+0.0im

julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int}, Any}((1, 2) => 1, (3, 4) => 3, (2, 3) => 2))
Dict{Graphs.SimpleGraphs.SimpleEdge{Int64}, Any} with 3 entries:
  Edge 1 => 2 => 1
  Edge 3 => 4 => 3
  Edge 2 => 3 => 2
```
"""
function tuple2edge(tupledict::Dict{Tuple{Int, Int},T}) where T
    edgedict = Dict{Graphs.SimpleGraphs.SimpleEdge{Int},T}()

    for (key,val) in tupledict
        edgedict[Graphs.Edge(key[1],key[2])]=val
    end

    return edgedict
end


"""
    tuple2edge(tupledict)

Converts a dictionary whose keys are edges specified by tuples of integers to
a dictionary whose keys are Graphs edges. The values associated with each key
are preserved.

# Examples
```jldoctest
julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Int}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3))
Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, Int64} with 2 entries:
  (Edge 1 => 2, Edge 3 => 4) => 1
  (Edge 5 => 6, Edge 7 => 8) => 3

julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Float64}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3))
Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, Float64} with 2 entries:
  (Edge 1 => 2, Edge 3 => 4) => 1.0
  (Edge 5 => 6, Edge 7 => 8) => 3.0

julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Complex{Float64}}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3))
Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, ComplexF64} with 2 entries:
  (Edge 1 => 2, Edge 3 => 4) => 1.0+0.0im
  (Edge 5 => 6, Edge 7 => 8) => 3.0+0.0im

julia> JosephsonCircuits.tuple2edge(Dict{Tuple{Int, Int, Int, Int}, Any}((1, 2, 3, 4) => 1, (5, 6, 7, 8) => 3))
Dict{Tuple{Graphs.SimpleGraphs.SimpleEdge{Int64}, Graphs.SimpleGraphs.SimpleEdge{Int64}}, Any} with 2 entries:
  (Edge 1 => 2, Edge 3 => 4) => 1
  (Edge 5 => 6, Edge 7 => 8) => 3
```
"""
function tuple2edge(tupledict::Dict{Tuple{Int, Int, Int, Int},T}) where T
    edgedict = Dict{
        Tuple{
            Graphs.SimpleGraphs.SimpleEdge{Int},
            Graphs.SimpleGraphs.SimpleEdge{Int}
        },
        T
    }()

    for (key,val) in tupledict
        edgedict[(Graphs.Edge(key[1],key[2]),Graphs.Edge(key[3],key[4]))]=val
    end

    return edgedict
end
