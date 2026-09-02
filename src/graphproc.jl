# The graph of a compiled circuit: the branches the incidence matrix is
# built from, a spanning tree of them, and the loops the remaining branches
# close.

"""
    CircuitGraph(edge2indexdict, Rbn, searray, cearray, glearray, lvarray,
        isolatednodes, gl, Nbranches)

The graph of the branch carrying components of a circuit, as computed by
[`calccircuitgraph`](@ref).

# Fields
- `edge2indexdict`: maps a branch `(node1, node2)` in either orientation to
    its branch index, the row of `Rbn` it occupies.
- `Rbn`: the sparse oriented incidence matrix, `Nbranches` by
    `Nnodes - 1`; the ground node column is omitted.
- `searray`: the edges of the spanning tree, as `(node1, node2)` tuples.
- `cearray`: the closure branches, the edges not in the spanning tree.
- `glearray`: all edges, spanning tree first and closure branches after.
- `lvarray`: for each closure branch, the vertices of the loop it closes
    through the spanning tree (empty for a loop of only two vertices).
    Empty altogether when `calccircuitgraph` was called with
    `loops = false`.
- `isolatednodes`: nodes which appear in the graph but have no branch to
    any other node.
- `gl`: the undirected `Graphs.SimpleGraph` of all branches.
- `Nbranches`: the number of branches, `size(Rbn, 1)`.
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
    calccircuitgraph(compiledcircuit::CompiledCircuit; loops = true)

Compute the [`CircuitGraph`](@ref) of a compiled circuit: the incidence
matrix, a spanning tree, the closure branches, and (when `loops = true`)
the loop each closure branch closes.

The graph is built from the branches of the inductive components, the
Josephson junctions, the current and voltage sources and the ports; see
[`extractbranches`](@ref) for the list. Nothing in the solvers reads the
loops, and enumerating them costs a tree walk per closure branch, so a
caller which does not need them passes `loops = false`.

See also [`calcgraphs`](@ref).

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

    return calcgraphs(branchvector, compiledcircuit.Nnodes;
        loops = loops)

end

"""
    calcgraphs(Ledgearray::Array{Tuple{Int, Int}, 1}, Nnodes::Int;
        loops = true)

Build the [`CircuitGraph`](@ref) of the branches `Ledgearray`, given as
`(node1, node2)` tuples over `Nnodes` nodes with ground being node 1.

Each connected component of the branch graph gets a minimum spanning tree
(Kruskal, on unit weights) rooted at its first vertex. The edges not in
the tree are the closure branches, and when `loops = true` the loop of
each closure branch is the unique path between its endpoints through the
tree. The oriented incidence matrix is assembled from the tree edges
followed by the closure branches, with vertices added for any nodes which
carry no branch so that the matrix has `Nnodes - 1` columns.
"""
function calcgraphs(Ledgearray::Array{Tuple{Int, Int}, 1}, Nnodes::Int;
        loops::Bool = true)

    gl = Graphs.SimpleGraphFromIterator(tuple2edge(Ledgearray))

    searray = Vector{Tuple{Int, Int}}(undef, 0)
    cearray = Vector{Tuple{Int, Int}}(undef, 0)
    lvarray = Vector{Vector{Int}}(undef, 0)
    glearray = Vector{Tuple{Int, Int}}(undef, 0)
    isolatednodes = Vector{Int}(undef,0)

    # one spanning tree per connected component of the branch graph
    for v in Graphs.connected_components(gl)

        # a component of one vertex is a node with no branch to anywhere
        if length(v) == 1
            push!(isolatednodes,v[1])
        end

        # the subgraph of this component; `vmap` takes its local vertex
        # numbers back to the circuit's node indices
        gli, vmap = Graphs.induced_subgraph(gl,v)

        # a minimum spanning tree of the component
        si = Graphs.SimpleGraph(Graphs.kruskal_mst(gli))

        # the closure branches: every edge not in the tree
        ci = collect(Graphs.edges(Graphs.difference(gli,si)))

        # root the tree once so that each loop below costs a walk of its
        # own length rather than a search of the whole tree
        parent, depth = loops ? rootedtree(si) : (Int[], Int[])

        for cj in ci
            push!(cearray,(vmap[Graphs.dst(cj)],vmap[Graphs.src(cj)]))

            # The loop of a closure branch is the unique path through the
            # spanning tree between its endpoints, closed by the branch
            # itself. Nothing in the solvers reads the loops, so they are
            # only computed on request; the tree and the closure branches
            # are always built because the incidence matrix needs them.
            loops || continue

            cyc = treepath(parent, depth, Graphs.src(cj), Graphs.dst(cj))
            # a closure branch parallel to a tree edge closes a loop of two
            # vertices, which is recorded as empty
            push!(lvarray, length(cyc) > 2 ? vmap[cyc] : Int[])
        end

        # orient the tree edges away from local vertex 1 by a breadth first
        # search, so that every branch has a definite direction
        if Graphs.ne(si) == 0
            sid = si
        else
            sid = Graphs.SimpleGraph(Graphs.bfs_tree(si,1))
        end

        # add the closure branches back to get the directed graph of every
        # branch in this component
        glid = copy(sid)
        for cj in ci
            Graphs.add_edge!(glid, Graphs.dst(cj), Graphs.src(cj))
        end

        # record the edges in the circuit's node numbering
        for e in Graphs.edges(sid)
            push!(searray,(vmap[Graphs.src(e)],vmap[Graphs.dst(e)]))

        end
        for e in Graphs.edges(glid)
            push!(glearray,(vmap[Graphs.src(e)],vmap[Graphs.dst(e)]))

        end
    end

    gl2 = Graphs.SimpleDiGraphFromIterator(tuple2edge(glearray))
    # A graph built from an edge list has as many vertices as the largest
    # node index in it. Nodes above that, which can happen when only
    # capacitors or resistors touch the highest numbered nodes, are added so
    # that the incidence matrix has a column for every node. (`gl` and `gl2`
    # have the same vertex count, both being built from the same branches.)
    if Graphs.nv(gl2) < Nnodes
        Graphs.add_vertices!(gl2,Nnodes-Graphs.nv(gl))
    end

    edge2indexdict = edge2index(gl2)

    # the oriented incidence matrix, transposed to branches by nodes, with
    # the ground node column (node 1) dropped
    Rbn=sparse(transpose(Graphs.incidence_matrix(gl2,oriented=true)))[:,2:end]

    Nbranches = Graphs.ne(gl2)

    return CircuitGraph(edge2indexdict, Rbn, searray, cearray, glearray,
        lvarray, isolatednodes, gl, Nbranches)
end

"""
    rootedtree(tree)

Root the tree `tree` at vertex 1 by breadth first search and return the
parent and the depth of every vertex (parent 0 and depth 0 for the root;
depth -1 for a vertex unreachable from it).

[`treepath`](@ref) walks between two vertices using both: the depths bring
the two ends to the same level and the parents carry them up to where they
meet.
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

The vertices of the unique path from `u` to `v` through the tree described
by `parent` and `depth` (see [`rootedtree`](@ref)), `u` first and `v` last.
Returns an empty vector when either vertex is unreachable from the root.

Together with a closure branch from `v` back to `u` this path is the
fundamental loop of that branch.

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

A dictionary from the `(src, dst)` tuple of each edge of `graph`, in both
orientations, to the position of that edge in `Graphs.edges(graph)`. The
positions are the branch indices of the incidence matrix.
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

Convert a vector of `(src, dst)` tuples to a vector of `Graphs` edges.

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

Convert a vector of `(src1, dst1, src2, dst2)` tuples to a vector of pairs
of `Graphs` edges.

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

Convert a dictionary keyed by `(src, dst)` tuples to one keyed by `Graphs`
edges, keeping the values.
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

Convert a dictionary keyed by `(src1, dst1, src2, dst2)` tuples to one
keyed by pairs of `Graphs` edges, keeping the values.
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

"""
    extractbranches(componenttypes::Vector{Symbol},nodeindexarray::Matrix{Int})

The `(node1, node2)` branches of the components which define the circuit
graph: inductors (`:L`), Josephson junctions (`:Lj`), legacy nonlinear
elements (`:NL`), current sources (`:I`), ports (`:P`) and voltage sources
(`:V`). Capacitors, resistors and mutual inductors do not create branches.

Components sharing a branch produce duplicate tuples; the graph
construction in [`calcgraphs`](@ref) merges them.

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

Push the branches described in [`extractbranches`](@ref) onto the empty
vector `branchvector`.
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
