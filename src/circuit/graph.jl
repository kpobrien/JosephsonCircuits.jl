# The graph of a compiled circuit: the branches its incidence matrix is
# built from, and the connected components of its nodes.

"""
    nodecomponents(Nnodes::Int, edges)

The connected components of `Nnodes` nodes joined by `edges`, given as
`(node1, node2)` pairs, leaving out the component which contains ground
(node 1). Each component is a sorted vector of node indices and the
components are sorted by their lowest node.

A union-find with path halving, unioning toward the lower node index so
that ground is the root of its own component. Which components are
floating depends on which branches count as an edge, so each caller
supplies its own: the static flux stiffness of the direct current gauge
(`calcstaticfluxcomponents`), the conduction paths of the transient
(`transientfloatingcomponents`), and every element and block port of the
subnetworks harmonic balance refuses (`isolatedsubnetworks`).
"""
function nodecomponents(Nnodes::Int, edges)
    parent = collect(1:Nnodes)
    function findroot(i::Int)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    for (n1, n2) in edges
        a, b = findroot(n1), findroot(n2)
        a == b || (parent[max(a, b)] = min(a, b))
    end
    components = Dict{Int,Vector{Int}}()
    for node in 2:Nnodes
        root = findroot(node)
        root == 1 && continue
        push!(get!(components, root, Int[]), node)
    end
    out = collect(values(components))
    foreach(sort!, out)
    return sort!(out; by = first)
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
    extractbranches(componenttypes::Vector{Symbol},nodeindexarray::Matrix{Int})

The `(node1, node2)` branches of the components which define the circuit
graph: inductors (`:L`), Josephson junctions (`:Lj`), current sources
(`:I`) and ports (`:P`). Capacitors, resistors and mutual inductors do
not create branches.

Components sharing a branch produce duplicate tuples, which
[`circuittopology`](@ref) merges into one branch.

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

    if  length(componenttypes) != size(nodeindexarray,2)
        throw(DimensionMismatch(lazy"componenttypes must have the same length as the number of node indices"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch(lazy"the length of the first axis must be 2"))
    end

    branchvector = Tuple{Int,Int}[]
    for i in eachindex(componenttypes)
        if componenttypes[i] in (:Lj, :L, :I, :P)
            push!(branchvector,(nodeindexarray[1,i],nodeindexarray[2,i]))
        end
    end

    return branchvector
end
