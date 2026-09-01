
"""
    CircuitGraph(edge2indexdict, Rbn, searray, cearray, glearray, lvarray,
        isolatednodes, gl, Nbranches)

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
    calccircuitgraph(compiledcircuit::CompiledCircuit)

Calculate the superconducting spanning tree, incidence matrix, closure branches,
and loops from the parsed and sorted circuit.

See also [`CircuitGraph`](@ref), [`calcgraphs`](@ref), and [`extractbranches`](@ref) 
for more explanation.

# Examples
```jldoctest
circuit = Circuit(
    [:p1 => Port(1; Z0 = :Rleft),
     :i1 => CurrentSource(:Ipump),
     :l1 => Inductor(:L),
     :jj => JosephsonJunction(:Lj),
     :cj => Capacitor(:Cj),
     :gnd => Ground()],
    [[(:p1, 1), (:i1, 1), (:l1, 1)],
     [(:l1, 2), (:jj, 1), (:cj, 1)],
     [(:p1, 2), (:i1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
psc = JosephsonCircuits.compile(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
JosephsonCircuits.comparestruct(cg,JosephsonCircuits.CircuitGraph(Dict((3, 2) => 3, (1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1, (2, 3) => 3), JosephsonCircuits.SparseArrays.sparse([1, 3, 2, 3], [1, 1, 2, 2], [1, -1, 1, 1], 3, 2), [(1, 2), (1, 3)], [(3, 2)], [(1, 2), (1, 3), (2, 3)], [[2, 1, 3]], Int64[], JosephsonCircuits.Graphs.SimpleGraphs.SimpleGraph{Int64}(3, [[2, 3], [1, 3], [1, 2]]), 3))
# output
true
```
"""
function calccircuitgraph(compiledcircuit::CompiledCircuit;
        loops::Bool = true)

    branchvector = extractbranches(compiledcircuit.componenttypes,
                                compiledcircuit.nodeindices)

    # calculate the graph of inductive components glelist, the
    # superconducting spanning tree selist, and the list of loop
    # indices celist. 
    return calcgraphs(branchvector, compiledcircuit.Nnodes;
        loops = loops)

end

"""
    calcgraphs(Ledgearray::Array{Tuple{Int, Int}, 1}, Nnodes::Int)

Calculate the superconducting spanning tree, closure branches, and loops.
Accepts the graph of linear inductors and Josephson junctions. Outputs lists
of edges that can be used to generate graphs.
"""
function calcgraphs(Ledgearray::Array{Tuple{Int, Int}, 1}, Nnodes::Int;
        loops::Bool = true)

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
        #si = Graphs.SimpleGraph(prim_mst(gli))

        #calculate the closure branches
        ci = collect(Graphs.edges(Graphs.difference(gli,si)))

        # the spanning tree, rooted once, so that each closure branch costs
        # the length of its own loop rather than a search of the whole tree
        parent, depth = loops ? rootedtree(si) : (Int[], Int[])
    
        #find the loop indices associated with each closure branch by starting
        #with the superconducting spanning tree (which has no loops), then
        #adding a closure branch and looking for the cycle this creates. 
        for cj in ci
            #push!(cearray,Edge(vmap[dst(cj)],vmap[src(cj)]))
            push!(cearray,(vmap[Graphs.dst(cj)],vmap[Graphs.src(cj)]))

            # The loop of a closure branch is the unique path through the
            # spanning tree between its endpoints, closed by the branch
            # itself. A tree has exactly one such path, so it is walked
            # rather than searched for: the enumeration this replaces added
            # the branch back to the tree and looked for cycles of at most
            # ten edges, which silently returned no loop for a longer one.
            #
            # Nothing outside this file reads `lvarray`, so a caller which
            # does not want the loops says so and pays for none of this. The
            # spanning tree and the closure branches are still built,
            # because the incidence matrix is derived from them.
            loops || continue

            cyc = treepath(parent, depth, Graphs.src(cj), Graphs.dst(cj))
            # a closure branch parallel to a tree edge closes a two vertex
            # loop, which the enumeration this replaces also dropped
            push!(lvarray, length(cyc) > 2 ? vmap[cyc] : Int[])
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

    return CircuitGraph(edge2indexdict, Rbn, searray, cearray, glearray,
        lvarray, isolatednodes, gl, Nbranches)
end

"""
    rootedtree(tree)

Root a spanning tree at its first vertex, returning the parent and the depth
of every vertex.

Both are needed to walk between two vertices: the depths bring the two walks
to the same level and the parents carry them up to where they meet.
"""
function rootedtree(tree)
    n = Graphs.nv(tree)
    parent = zeros(Int, n)
    depth = fill(-1, n)
    n == 0 && return parent, depth
    depth[1] = 0
    queue = Int[1]
    head = 1
    while head <= length(queue)
        u = queue[head]; head += 1
        for w in Graphs.neighbors(tree, u)
            depth[w] < 0 || continue
            depth[w] = depth[u] + 1
            parent[w] = u
            push!(queue, w)
        end
    end
    return parent, depth
end

"""
    treepath(parent, depth, u, v)

The vertices of the unique path from `u` to `v` through a rooted tree,
`u` first and `v` last.

Together with the closure branch from `v` back to `u` this is the
fundamental loop of that branch, and there is exactly one of them, which is
why it can be walked instead of searched for.

# Examples
```jldoctest
julia> JosephsonCircuits.treepath([0, 1, 2, 1], [0, 1, 2, 1], 3, 4)
4-element Vector{Int64}:
 3
 2
 1
 4
```
"""
function treepath(parent::Vector{Int}, depth::Vector{Int}, u::Integer,
        v::Integer)
    a, b = Int(u), Int(v)
    (depth[a] < 0 || depth[b] < 0) && return Int[]
    up = Int[a]
    down = Int[b]
    while depth[a] > depth[b]
        a = parent[a]; push!(up, a)
    end
    while depth[b] > depth[a]
        b = parent[b]; push!(down, b)
    end
    while a != b
        a = parent[a]; push!(up, a)
        b = parent[b]; push!(down, b)
    end
    pop!(down)                       # the meeting vertex, already in `up`
    return vcat(up, reverse(down))
end

"""
    edge2index(graph::Graphs.SimpleDiGraph{Int})

Generate a dictionary where the tuple of nodes defining an edge of a graph
is the key and the value is an index. The index gives the order the edge
is found when iterating over the edges of the graph. The same index is used
for both orderings of source and destination nodes on the edge. We don't care
about the ordering of the indices as long as they are sequential and unique.
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
    tuple2edge(tupledict::Dict{Tuple{Int, Int},T})

Converts a dictionary whose keys are edges specified by tuples of integers to
a dictionary whose keys are Graphs edges. The values associated with each key
are preserved.
"""
function tuple2edge(tupledict::Dict{Tuple{Int, Int},T}) where T
    edgedict = Dict{Graphs.SimpleGraphs.SimpleEdge{Int},T}()

    for (key,val) in tupledict
        edgedict[Graphs.Edge(key[1],key[2])]=val
    end

    return edgedict
end

"""
    tuple2edge(tupledict::Dict{Tuple{Int, Int, Int, Int},T})

Converts a dictionary whose keys are edges specified by tuples of integers to
a dictionary whose keys are Graphs edges. The values associated with each key
are preserved.
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

# The branches the incidence matrix is built from, read off the compiled
# tables. Only `calccircuitgraph` above uses them.
"""
    extractbranches(componenttypes::Vector{Symbol},nodeindexarray::Matrix{Int})

Return an array of tuples of pairs of node indices (branches) which we will
use to calculate the incidence matrix.

This will contain duplicates if multiple components are on the same branch. All
checking for duplicate branches will occur in the graph procesing code.

NOTE: the list of component types considered to lie on branches is hardcoded.

# Examples
```jldoctest
julia> JosephsonCircuits.extractbranches([:P,:I,:R,:C,:Lj,:C],[2 2 2 2 3 3; 1 1 1 3 1 1])
3-element Vector{Tuple{Int64, Int64}}:
 (2, 1)
 (2, 1)
 (3, 1)
```
"""
function extractbranches(componenttypes::Vector{Symbol},nodeindexarray::Matrix{Int})

    branchvector = Array{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)},1}(undef,0)
    extractbranches!(branchvector,componenttypes,nodeindexarray)

    return branchvector
end

"""
    extractbranches!(branchvector::Vector,componenttypes::Vector{Symbol},
        nodeindexarray::Matrix{Int})

Append tuples consisting of a pair of node indices (branches) which we will
use to calculate the incidence matrix. Appends the tuples to branchvector.
"""
function extractbranches!(branchvector::Vector,componenttypes::Vector{Symbol},nodeindexarray::Matrix{Int})

    if  length(componenttypes) != size(nodeindexarray,2)
        throw(DimensionMismatch(lazy"componenttypes must have the same length as the number of node indices"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch(lazy"the length of the first axis must be 2"))
    end

    if length(branchvector) != 0
        throw(DimensionMismatch(lazy"branchvector should be length zero"))
    end

    allowedcomponenttypes = [:Lj,:NL,:L,:I,:P,:V]
    for i in eachindex(componenttypes)
        type = componenttypes[i]
        if type in allowedcomponenttypes
            push!(branchvector,(nodeindexarray[1,i],nodeindexarray[2,i]))
        end
    end

    return nothing
end
