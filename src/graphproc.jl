
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
and loops.

See also [`CircuitGraph`](@ref), [`calcgraphs`](@ref), and [`extractbranches`](@ref) 
for more explanation.

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

Calculate the superconducting spanning tree, closure branches,
and loops. Accepts the graph of linear inductors and Josephson
junctions. Outputs lists of edges that can be used to generate
graphs.

"""
function calcgraphs(Ledgearray::Array{Tuple{Int64, Int64}, 1},Nnodes::Int64)

    gl = Graphs.SimpleGraphFromIterator(tuple2edge(Ledgearray))

    searray = Array{Tuple{Int64, Int64}, 1}(undef, 0)
    cearray = Array{Tuple{Int64, Int64}, 1}(undef, 0)
    lvarray = Array{Array{Int64,1}, 1}(undef, 0)
    glearray = Array{Tuple{Int64, Int64}, 1}(undef, 0)
    isolatednodes = Array{Int64,1}(undef,0)

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
            if length(ul) > 1
                error("There should only be one loop associated with 
                     each closure branch.")
            elseif length(ul) < 1
                #error("No loops found.")
                # println("Warning: Loop exists but max loop size too small to find vertices.")
                push!(lvarray,[])
            else
    #             println("closure branch: ",cj,", loop vertices:", first(ul))
                push!(lvarray,vmap[first(ul)])
            end
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
    # if mode vertices were found when parsing all the compoments than there are in
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
    edge2index(graph::SimpleDiGraph{Int64})

Generate a dictionary that maps an edge to an index.
"""
function edge2index(graph::Graphs.SimpleDiGraph{Int64})
    edge2indexdict = Dict{Tuple{Int64, Int64},Int64}()
    i = 1
    for e in Graphs.edges(graph)
        edge2indexdict[(Graphs.src(e),Graphs.dst(e))] = i
        edge2indexdict[(Graphs.dst(e),Graphs.src(e))] = i
        i+=1
    end
    return edge2indexdict
end


"""
    tuple2edge(tuplearray)

Changes edges specified with tuples of integers to LightGraph edges.
"""
function tuple2edge(tuplearray::Array{Tuple{Int64, Int64}, 1})
    edgearray = Array{Graphs.SimpleGraphs.SimpleEdge{Int64}, 1}(undef, 0)

    for i = 1:length(tuplearray)
        push!(edgearray,Graphs.Edge(tuplearray[i][1],tuplearray[i][2]))
    end
    return edgearray
end

function tuple2edge(tuplearray::Array{Tuple{Int64, Int64, Int64, Int64}, 1})
    edgearray = Array{
        Tuple{
            Graphs.SimpleGraphs.SimpleEdge{Int64},
            Graphs.SimpleGraphs.SimpleEdge{Int64}
        }
        ,1
    }(undef, 0)

    for i = 1:length(tuplearray)
        push!(
            edgearray,
            (
                Graphs.Edge(tuplearray[i][1],tuplearray[i][2]),
                Graphs.Edge(tuplearray[i][3],tuplearray[i][4])
            )
        )
    end

    return edgearray
end

"""
    tuple2edge(tupledict)

Changes dictionary keys specified with tuples of integers to LightGraph edges.
"""
function tuple2edge(tupledict::Dict{Tuple{Int64, Int64},Any})
    edgedict = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64},Any}()

    for (key,val) in tupledict
        edgedict[Graphs.Edge(key[1],key[2])]=val
    end

    return edgedict
end

function tuple2edge(tupledict::Dict{Tuple{Int64, Int64},Complex{Float64}})
    edgedict = Dict{Graphs.SimpleGraphs.SimpleEdge{Int64},Complex{Float64}}()

    for (key,val) in tupledict
        edgedict[Graphs.Edge(key[1],key[2])]=val
    end

    return edgedict
end

function tuple2edge(tupledict::Dict{Tuple{Int64, Int64, Int64, Int64},Any})
    edgedict = Dict{
        Tuple{
            Graphs.SimpleGraphs.SimpleEdge{Int64},
            Graphs.SimpleGraphs.SimpleEdge{Int64}
        },
        Any
    }()

    for (key,val) in tupledict
        edgedict[(Graphs.Edge(key[1],key[2]),Graphs.Edge(key[3],key[4]))]=val
    end

    return edgedict
end

function tuple2edge(tupledict::Dict{Tuple{Int64, Int64, Int64, Int64},Complex{Float64}})
    edgedict = Dict{
        Tuple{
            Graphs.SimpleGraphs.SimpleEdge{Int64},
            Graphs.SimpleGraphs.SimpleEdge{Int64}
        },
        Complex{Float64}
    }()

    for (key,val) in tupledict
        edgedict[(Graphs.Edge(key[1],key[2]),Graphs.Edge(key[3],key[4]))]=val
    end

    return edgedict
end
