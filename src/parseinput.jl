
"""
    ParsedSortedCircuit(nodeindices::Matrix{Int},
        nodenames::Vector{String},
        mutualinductorbranchnames::Vector{String},
        componentnames::Vector{String},
        componenttypes::Vector{Symbol}, componentvalues::Vector,
        componentnamedict::Dict{String, Int}, Nnodes::Int)

A simple structure to hold the parsed and sorted circuit. See also
[`parsesortcircuit`](@ref), [`parsecircuit`](@ref), and [`sortnodes`](@ref)
for more explanation.

# Fields
- `nodeindices::Matrix{Int}`: sorted array of node indices (where the length
    of the first axis is 2).
- `nodenames::Vector{String}`: the sorted unique node names.
- `mutualinductorbranchnames::Vector{String}`: the inductors coupled by the
    mutual inductors.
- `componentnames::Vector{String}`: component names.
- `componenttypes::Vector{Symbol}`: the component (electrical engineering) types.
- `componentvalues::Vector`: the component values.
- `componentnamedict::Dict{String, Int}`: names of the components as keys and
    the index at which the component occurs as the value.
- `Nnodes::Int`: number of nodes including the ground node.

# Examples
```jldoctest
@variables Ipump Rleft L1 K1 L2 C2
println(JosephsonCircuits.ParsedSortedCircuit(
    [2 2 2 2 0 3 3; 1 1 1 1 0 1 1],
    ["0", "1", "2"],
    ["L1", "L2"],
    ["P1", "I1", "R1", "L1", "K1", "L2", "C2"],
    [:P, :I, :R, :L, :K, :L, :C],
    Num[1, Ipump, Rleft, L1, K1, L2, C2],
    Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5),
    3))

# output
JosephsonCircuits.ParsedSortedCircuit([2 2 2 2 0 3 3; 1 1 1 1 0 1 1], ["0", "1", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Num[1, Ipump, Rleft, L1, K1, L2, C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3)
```
"""
struct ParsedSortedCircuit
    nodeindices::Matrix{Int}
    nodenames::Vector{String}
    mutualinductorbranchnames::Vector{String}
    componentnames::Vector{String}
    componenttypes::Vector{Symbol}
    componentvalues::Vector
    componentnamedict::Dict{String, Int}
    Nnodes::Int
end

"""
    parsesortcircuit(circuit; sorting = :name)

Parse and sort the circuit. See [`parsecircuit`](@ref), [`sortnodes`](@ref)
for more explanation.

# Arguments
- `circuit`: vector of tuples each of which contain the component name, the
    first node, the second node, and the component value. The first three must
    be strings.

# Keywords
- `sorting = :name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
- `sorting = :number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
- `sorting = :none`: Don't perform any sorting except to place the ground node
    first. In other words, order the nodes in the order they are found in
    `circuit`.

# Examples
```jldoctest
@variables Ipump Rleft Cc Lj Cj
circuit = Tuple{String,String,String,Num}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc))
push!(circuit,("Lj1","2","0",Lj))
push!(circuit,("C2","2","0",Cj))
println(parsesortcircuit(circuit))

# output
JosephsonCircuits.ParsedSortedCircuit([2 2 2 2 3 3; 1 1 1 3 1 1], ["0", "1", "2"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Num[1, Ipump, Rleft, Cc, Lj, Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3)
```
```jldoctest
@variables Ipump Rleft L1 L2 C2
Kfun(L) = sin(L);@register_symbolic Kfun(L1)
circuit = Tuple{String,String,String,Num}[]
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1)) 
push!(circuit,("K1","L1","L2",Kfun(L1)))
push!(circuit,("L2","2","0",L2)) 
push!(circuit,("C2","2","0",C2))
println(parsesortcircuit(circuit))

# output
JosephsonCircuits.ParsedSortedCircuit([2 2 2 2 0 3 3; 1 1 1 1 0 1 1], ["0", "1", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Num[1, Ipump, Rleft, L1, Kfun(L1), L2, C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3)
```
"""
function parsesortcircuit(circuit; sorting = :number)

    # parse the circuit components
    parsedcircuit = parsecircuit(circuit)

    # sort the nodes
    nodenames,nodeindices = sortnodes(
        parsedcircuit.uniquenodevector,
        parsedcircuit.nodeindexvector,
        sorting=sorting)

    return ParsedSortedCircuit(
        nodeindices,
        nodenames,
        parsedcircuit.mutualinductorbranchnames,
        parsedcircuit.componentnames,
        parsedcircuit.componenttypes,
        parsedcircuit.componentvalues,
        parsedcircuit.componentnamedict,
        parsedcircuit.Nnodes)
end

"""
    ParsedCircuit(nodeindexvector::Vector{Int},
        uniquenodevector::Vector{String},
        mutualinductorbranchnames::Vector{String},
        componentnames::Vector{String}, componenttypes::Vector{Symbol},
        componentvalues::Vector, componentnamedict::Dict{String, Int},
        Nnodes::Int)

A simple structure to hold the parsed circuit including a vector of node
indices, the unique node names, the inductors coupled by the mutual inductors,
the component names, the component types, the values of the components, a
dictionary of the names of the components as keys and the index at which the
component occurs as the value, and dictionaries for the ports and resistors
where the pair of nodes is the key and value is the component value.

See also [`parsecircuit`](@ref).

# Fields
- `nodeindexvector::Vector{Int}`: sorted vector of node indices where the
    two nodes for each component occur as consecutive elements (pairs).
- `uniquenodevector::Vector{String}`: the unique node names.
- `mutualinductorbranchnames::Vector{String}`: the inductors coupled by the
    mutual inductors.
- `componentnames::Vector{String}`: component names.
- `componenttypes::Vector{Symbol}`: the component (electrical engineering)
    types.
- `componentvalues::Vector`: the component values.
- `componentnamedict::Dict{String, Int}`: names of the components as keys and
    the index at which the component occurs as the value.
- `Nnodes::Int`: number of nodes including the ground node.

# Examples
```jldoctest
@variables Ipump Rleft L1 K1 L2 C2
println(JosephsonCircuits.ParsedCircuit(
    [1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2],
    ["1", "0", "2"], ["L1", "L2"],
    ["P1", "I1", "R1", "L1", "K1", "L2", "C2"],
    [:P, :I, :R, :L, :K, :L, :C],
    Num[1, Ipump, Rleft, L1, K1, L2, C2],
    Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5),
    3))

# output
JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2], ["1", "0", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Num[1, Ipump, Rleft, L1, K1, L2, C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3)
```
"""
struct ParsedCircuit
    nodeindexvector::Vector{Int}
    uniquenodevector::Vector{String}
    mutualinductorbranchnames::Vector{String}
    componentnames::Vector{String}
    componenttypes::Vector{Symbol}
    componentvalues::Vector
    componentnamedict::Dict{String, Int}
    Nnodes::Int
end

"""
    parsecircuit(circuit)

Parse `circuit` which is a vector where each element contains a tuple with the
component name, the first node, the second node, and the component value.
Component values can be numbers, symbols, or symbolic variables (including
symbolic functions).

The nodes can be arbitrary strings for SPICE compatibility. Integers are also
supported but are converted internally to strings. The ground node is "0" and
is required. Specifying the type of the vector `circuit` is optional; although,
typically a vector with a type union is preferable to an array of type Any.

# Arguments
- `circuit`: vector of tuples each of which contain the component name, the
    first node, the second node, and the component value.

# Examples
```jldoctest
@variables Ipump Rleft L1 K1 L2 C2
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1))
push!(circuit,("K1","L1","L2",K1))
push!(circuit,("L2","2","0",L2))
push!(circuit,("C2","2","0",C2))
parsecircuit(circuit)

# output
JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2], ["1", "0", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Num[1, Ipump, Rleft, L1, K1, L2, C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3)
```
```jldoctest
@variables Ipump Rleft L1 L2 C2
Kfun(L) = sin(L);@register_symbolic Kfun(L1)
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1))
push!(circuit,("K1","L1","L2",Kfun(L1)))
push!(circuit,("L2","2","0",L2))
push!(circuit,("C2","2","0",C2))
parsecircuit(circuit)

# output
JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2], ["1", "0", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Num[1, Ipump, Rleft, L1, Kfun(L1), L2, C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3)
```
```jldoctest
circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",:Ipump))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("C1","1","2",:Cc))
push!(circuit,("Lj1","2","0",:Lj))
push!(circuit,("C2","2","0",:Cj))
parsecircuit(circuit)

# output
JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 3, 3, 2, 3, 2], ["1", "0", "2"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Union{Int64, Symbol, ComplexF64}[1, :Ipump, :Rleft, :Cc, :Lj, :Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3)
```
```jldoctest
circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
push!(circuit,("P1","One","0",1))
push!(circuit,("I1","One","0",:Ipump))
push!(circuit,("R1","One","0",:Rleft))
push!(circuit,("C1","One","Two",:Cc))
push!(circuit,("Lj1","Two","0",:Lj))
push!(circuit,("C2","Two","0",:Cj))
parsecircuit(circuit)

# output
JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 3, 3, 2, 3, 2], ["One", "0", "Two"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Union{Int64, Symbol, ComplexF64}[1, :Ipump, :Rleft, :Cc, :Lj, :Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3)
```
```jldoctest
circuit = []
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",:Ipump))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("C1","1","2",:Cc))
push!(circuit,("Lj1","2","0",:Lj))
push!(circuit,("C2","2","0",:Cj))
parsecircuit(circuit)

# output
JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 3, 3, 2, 3, 2], ["1", "0", "2"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Any[1, :Ipump, :Rleft, :Cc, :Lj, :Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3)
```
```jldoctest
circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",:Ipump))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("L1","1","0",:L1))
push!(circuit,("K1","L1","L2",:K1))
push!(circuit,("L2","2","0",:L2))
push!(circuit,("C2","2","0",:C2))
parsecircuit(circuit)

# output
JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2], ["1", "0", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Union{Int64, Symbol, ComplexF64}[1, :Ipump, :Rleft, :L1, :K1, :L2, :C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3)
```
"""
function parsecircuit(circuit)

    # the component types we can handle. Note: Lj must come before L otherwise
    # the JJs will be first matched as inductors. 
    # NOTE: voltage sources currently not supported
    allowedcomponents = ["Lj","L","C","K","I","R","P"]

    # check that we can properly parse these
    checkcomponenttypes(allowedcomponents)

    # turn these into symbols
    # allowedsymbols =  Symbol.(allowedcomponents)    
    allowedsymbols = Vector{Symbol}(undef,length(allowedcomponents))
    for (i,c) in enumerate(allowedcomponents)
        allowedsymbols[i]=Symbol(c)
    end

    # empty dictionary for all of the component names. to check if they
    # are unique.
    componentnamedict = Dict{String,Int}()
    sizehint!(componentnamedict,length(circuit))

    # vector to hold the component names
    componentnames = Vector{String}(undef,length(circuit))

    # vector to hold the component type symbols, eg :L, :Lj, etc. 
    componenttypes = Vector{Symbol}(undef,length(circuit))

    # vector to hold the component values. take the type from the circuit array.
    componentvalues = try
        Vector{fieldtype(eltype(circuit),4)}(undef,length(circuit))
    catch
        Vector{Any}(undef,length(circuit))
    end

    # vector to hold the nodes
    nodeindexvector = Vector{Int}(undef,2*length(circuit))

    # vector to hold the inductors which the mutual inductance connects
    mutualinductorbranchnames = Array{String,1}(undef,0)

    # dictionary of unique nodes where the key is the node and the value is the
    # node index
    uniquenodedict = Dict{String,Int}()
    sizehint!(uniquenodedict,length(circuit))

    # a vector where the value is the node and the position in the array
    # is given by the node index
    uniquenodevector = Vector{String}(undef,0)
    sizehint!(uniquenodevector,length(circuit))


    for (i,(name,node1,node2,value)) in enumerate(circuit)

        # find the label, nodes, and the value of this circuit component. each
        # line should have the following format:
        # (name,node1,node2,value or function definition)

        componenttypeindex = parsecomponenttype(name,allowedcomponents)
        componenttype=allowedsymbols[componenttypeindex]

        if haskey(componentnamedict,name)
           throw(ArgumentError("Name \"$(name)\" on line $(i) is not unique."))
        end

        # add the name to the dictionary
        componentnamedict[name] = i

        # store the component name, type, and value
        componentnames[i] = name
        componenttypes[i] = componenttype
        componentvalues[i] = value

        # mutual inductors are treated differently because their "nodes"
        # refer to inductors rather than actual nodes. add zeros to
        # nodeindexvector whenever we find an inductor in order
        # to keep their length the same as the other arrays.

        if componenttype == :K
            push!(mutualinductorbranchnames,node1)
            push!(mutualinductorbranchnames,node2)
            nodeindexvector[2*i-1] = 0
            nodeindexvector[2*i] = 0
        else
            nodeindexvector[2*i-1] = processnode(uniquenodedict,uniquenodevector,node1)
            nodeindexvector[2*i] = processnode(uniquenodedict,uniquenodevector,node2)
        end

    end

    return ParsedCircuit(nodeindexvector,uniquenodevector,mutualinductorbranchnames,
        componentnames, componenttypes, componentvalues, componentnamedict,length(uniquenodevector))
end

"""
    processnode(uniquenodedict::Dict{String, Int},
        uniquenodevector::Vector{String},node::String)

Return the node index when given a node. Add the node string
to the vector `uniquenodevector` and the dictionary `uniquenodedict` with the
node string as the key and the node index (index at which it appears in
`uniquenodevector`) as the value.

# Examples
```jldoctest
uniquenodedict = Dict("10" =>1)
uniquenodevector = ["10"]
println(JosephsonCircuits.processnode(uniquenodedict,uniquenodevector,"15"))
println(uniquenodevector)
println(uniquenodedict)

# output
2
["10", "15"]
Dict("10" => 1, "15" => 2)
```
```jldoctest
uniquenodedict = Dict("10" =>1)
uniquenodevector = ["10"]
println(JosephsonCircuits.processnode(uniquenodedict,uniquenodevector,"10"))
println(uniquenodevector)
println(uniquenodedict)

# output
1
["10"]
Dict("10" => 1)
```
"""
function processnode(uniquenodedict::Dict{String, Int},
    uniquenodevector::Vector{String},node::String)
    if !haskey(uniquenodedict,node)
        # if this is a new node, add to the unique node vector
        push!(uniquenodevector,node)
        
        # use the length plus one so it starts with one
        # and we can use them as indices in an array
        return uniquenodedict[node] = length(uniquenodedict)+1
    else
        return uniquenodedict[node]
    end
end

"""
    processnode(uniquenodedict::Dict{String, Int},
        uniquenodevector::Vector{String},node)

Return the node index when given a node. Add the node string
to the vector `uniquenodevector` and the dictionary `uniquenodedict` with the
node string as the key and the node index (index at which it appears in 
`uniquenodevector`) as the value. If "node" is not a string, make it a string.

# Examples
```jldoctest
uniquenodedict = Dict("10" =>1)
uniquenodevector = ["10"]
println(JosephsonCircuits.processnode(uniquenodedict,uniquenodevector,15))
println(uniquenodevector)
println(uniquenodedict)

# output
2
["10", "15"]
Dict("10" => 1, "15" => 2)
```
```jldoctest
uniquenodedict = Dict("10" =>1)
uniquenodevector = ["10"]
println(JosephsonCircuits.processnode(uniquenodedict,uniquenodevector,:A))
println(uniquenodevector)
println(uniquenodedict)

# output
2
["10", "A"]
Dict("A" => 2, "10" => 1)
```
"""
function processnode(uniquenodedict::Dict{String, Int},
    uniquenodevector::Vector{String},node)
    # if the node isn't a string, turn it into a string
    return processnode(uniquenodedict,uniquenodevector,string(node))
end

"""
    parsecomponenttype(name::String,allowedcomponents::Vector{String})

The first one or two characters of the component name in the string `name`
should match one of the strings in the vector `allowedcomponents`. Return the 
index first of the match found.

NOTE: if a two letter component appears in allowedcomponents after a one 
letter component with the same starting letter this function will match on the
first value.

# Examples
```jldoctest
julia> JosephsonCircuits.parsecomponenttype("L10",["Lj","L","C","K","I","R","P"])
2

julia> [JosephsonCircuits.parsecomponenttype(c,["Lj","L","C","K","I","R","P"]) for c in ["Lj","L","C","K","I","R","P"]]
7-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
 7

julia> JosephsonCircuits.parsecomponenttype("L10",["Lj","L","C","K","I","R","P"])
2
```
"""
function parsecomponenttype(name::String,allowedcomponents::Vector{String})

    # loop over the labels
    @inbounds for j in eachindex(allowedcomponents)
        l=allowedcomponents[j]
        if l[1] == name[1]
            if length(l) == 2
                if length(name) >= 2 && l[2] == name[2]
                    return j
                end
            elseif length(l) == 1
                return j
            else
                throw(ArgumentError("parsecomponenttype() currently only works for two letter components"))
            end
        end
    end
    throw(ArgumentError("No matching component found in allowedcomponents."))
end

"""
    checkcomponenttypes(allowedcomponents::Vector{String})

Check that each element in `allowedcomponents` is found at the correct place.
This will detect the case where a two letter component appears in 
`allowedcomponents` after a one letter component with the same starting letter.
The function parsecomponenttype() will match on the first value and this
function will throw an error.

# Examples
```jldoctest
julia> JosephsonCircuits.checkcomponenttypes(["Lj","L","C","K","I","R","P"])
true
```
"""
function checkcomponenttypes(allowedcomponents::Vector{String})
    for i in eachindex(allowedcomponents)
        if i != parsecomponenttype(allowedcomponents[i],allowedcomponents)
            throw(ArgumentError("Allowed components parsing check has failed for $(allowedcomponents[i]). This can happen if a two letter long component comes after a one letter component. Please reorder allowedcomponents."))
        end
    end
    return true
end

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
        throw(DimensionMismatch("componenttypes must have the same length as the number of node indices"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("the length of the first axis must be 2"))
    end

    if length(branchvector) != 0
        throw(DimensionMismatch("branchvector should be length zero"))
    end

    allowedcomponenttypes = [:Lj,:L,:I,:P,:V]
    for i in eachindex(componenttypes)
        type = componenttypes[i]
        if type in allowedcomponenttypes
            push!(branchvector,(nodeindexarray[1,i],nodeindexarray[2,i]))
        end
    end

    return nothing
end

"""
    findgroundnodeindex(uniquenodevector::Vector{String})

Find the index of the ground node.

# Examples
```jldoctest
julia> JosephsonCircuits.findgroundnodeindex(["1","0","2"])
2

julia> JosephsonCircuits.findgroundnodeindex(["1","2"])
0

julia> JosephsonCircuits.findgroundnodeindex(String[])
0
```
"""
function findgroundnodeindex(uniquenodevector::Vector{String})

    # find the ground node. error if we don't find it.
    # groundnodeindex = 0
    for i in eachindex(uniquenodevector)
        if uniquenodevector[i] == "0"
            # groundnodeindex = i
            return i
            # break
        end
    end

    # if groundnodeindex == 0
    #     throw(ArgumentError("No ground node found in netlist."))
    # end

    # return groundnodeindex
    return 0
end

"""
    calcnodesorting(uniquenodevector::Vector{String};sorting=:number)

Sort the unique node names in `uniquenodevector` according to the specified
sorting scheme, always placing the ground node at the beginning. Return the
indices which sort `uniquenodevector`.

# Keywords
- `sorting = :name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
- `sorting = :number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
- `sorting = :none`: Don't perform any sorting except to place the ground node
    first. In other words, order the nodes in the order they are found in
    `circuit`.

# Examples
```jldoctest
julia> JosephsonCircuits.calcnodesorting(["30","11","0","2"];sorting=:name)
4-element Vector{Int64}:
 3
 2
 4
 1

julia> JosephsonCircuits.calcnodesorting(["30","11","0","2"];sorting=:number)
4-element Vector{Int64}:
 3
 4
 2
 1

julia> JosephsonCircuits.calcnodesorting(["30","11","0","2"];sorting=:none)
4-element Vector{Int64}:
 3
 1
 2
 4
```
"""
function calcnodesorting(uniquenodevector::Vector{String};sorting=:number)

    # vector of indices for the sortperm. if sorting is nothing, use this
    # uniquenodevectorsortindices = ones(Int,length(uniquenodevector))
    uniquenodevectorsortindices = Vector{Int}(undef,length(uniquenodevector))
    # uniquenodevectorsortindices .= 1:length(uniquenodevector)
    for i in eachindex(uniquenodevectorsortindices)
        uniquenodevectorsortindices[i] = i
    end

    # sort according to the desired scheme
    if sorting == :name
        # sort the vector of unique node strings
        sortperm!(uniquenodevectorsortindices,uniquenodevector,initialized=true)

    elseif sorting == :number
        # convert the unique node strings to integers and sort those
        uniquenodevectorints = Vector{Int}(undef,length(uniquenodevector))
        for i in eachindex(uniquenodevectorints)
            parsednode = tryparse(Int,uniquenodevector[i])
            if !isnothing(parsednode)
                uniquenodevectorints[i] = parsednode
            else
                throw(ArgumentError("Failed to parse the nodes as integers. Try setting the keyword argument `sorting=:name` or `sorting=:none`."))
            end
        end
        sortperm!(uniquenodevectorsortindices, uniquenodevectorints, initialized=true)

    elseif sorting == :none
        # don't perform any sorting. keep the nodes in
        # order of first appearance, except move the
        # ground node first later as always
        nothing
    else
        throw(ArgumentError("Unknown sorting type."))
    end

    # find the ground node. error if we don't find it.
    groundnodeindex = findgroundnodeindex(uniquenodevector)

    if groundnodeindex == 0
        throw(ArgumentError("No ground node found in netlist."))
    end

    # if the ground index is not the first after sorting, make it first and
    # increment earlier node indices
    if uniquenodevectorsortindices[1] != groundnodeindex
        for i = 2:groundnodeindex
            j = groundnodeindex-i+2
            uniquenodevectorsortindices[j]= uniquenodevectorsortindices[j-1]
        end
        uniquenodevectorsortindices[1] = groundnodeindex
    end

    return uniquenodevectorsortindices
end

"""
    sortnodes(uniquenodevector::Vector{String},
        nodeindexvector::Vector{Int};sorting=:name)

Sort the unique node names in `uniquenodevector` according to the specified
sorting scheme, always placing the ground node at the beginning.

Return the sorted `uniquenodevector` and `nodeindexvector` (with the vector
reshaped from a vector of length 2*Nnodes into a matrix with dimensions 2 by
Nnodes).

# Keywords
- `sorting = :name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
- `sorting = :number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
- `sorting = :none`: Don't perform any sorting except to place the ground node
    first.

# Examples
```jldoctest
julia> nodenames,nodeindexarray=JosephsonCircuits.sortnodes(["101","0","111","11"],[1,2,1,2,1,2,1,3,3,2,3,2,4,1],sorting=:none);println(nodenames);println(nodeindexarray);
["0", "101", "111", "11"]
[2 2 2 2 3 3 4; 1 1 1 3 1 1 2]

julia> nodenames,nodeindexarray=JosephsonCircuits.sortnodes(["101","0","111","11"],[1,2,1,2,1,2,1,3,3,2,3,2,4,1],sorting=:name);println(nodenames);println(nodeindexarray);
["0", "101", "11", "111"]
[2 2 2 2 4 4 3; 1 1 1 4 1 1 2]

julia> nodenames,nodeindexarray=JosephsonCircuits.sortnodes(["101","0","111","11"],[1,2,1,2,1,2,1,3,3,2,3,2,4,1],sorting=:number);println(nodenames);println(nodeindexarray);
["0", "11", "101", "111"]
[3 3 3 3 4 4 2; 1 1 1 4 1 1 3]

julia> nodenames,nodeindexarray=JosephsonCircuits.sortnodes(["1", "0", "2"],[1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2],sorting=:number);println(nodenames);println(nodeindexarray);
["0", "1", "2"]
[2 2 2 2 0 3 3; 1 1 1 1 0 1 1]
```
"""
function sortnodes(uniquenodevector::Vector{String},
        nodeindexvector::Vector{Int};sorting=:name)

    nodeindices = zeros(eltype(nodeindexvector),2,length(nodeindexvector)÷2)

    uniquenodevectorsortindices=calcnodesorting(uniquenodevector;sorting=sorting)

    nodevectorsortindices = sortperm(uniquenodevectorsortindices)

    # for (i,j) in enumerate(eachindex(nodeindexvector))
    for (i,j) in enumerate(nodeindexvector)
        # if it's a mutual inductor the node index will be zero because the
        # mutual inductor is between two inductors not between two nodes.
        # it not a mutual inductor, assign the sorted node index.
        if j == 0
            nothing
        else
            nodeindices[i] = nodevectorsortindices[j]
        end
    end

    return uniquenodevector[uniquenodevectorsortindices],nodeindices
end

"""
    calcvaluetype(componenttypes::Vector{Symbol},componentvalues::Vector,
        components::Vector{Symbol};checkinverse::Bool=true)

Returns a zero length vector with the (computer science) type which will hold
a set of circuit components of the (electrical engineering) types given in
`components`. This function is not type stable by design, but exists to make
the later function calls type stable.

# Arguments
- `componenttypes::Vector{Symbol}`: the component (electrical engineering) types.
- `componentvalues::Vector`: the component values.
- `components::Vector{Symbol}`: find a (computer science) type which will
    hold the component (electrical engineering) types in this vector.

# Keywords
- `checkinverse = true`: also check the inverse of each element. This is
    useful if the type would be integer but we later want to take the inverse
    and want an array with a type that supports this operation.

# Examples
```jldoctest
julia> JosephsonCircuits.calcvaluetype([:R,:C,:R],[1,2,3],[:R])
Float64[]

julia> JosephsonCircuits.calcvaluetype([:R,:C,:R],[1,2,3+0.0im],[:R])
ComplexF64[]

julia> @variables R1 C1 R2;JosephsonCircuits.calcvaluetype([:R,:C,:R],[R1,C1,R2],[:R])
Num[]

julia> @syms R1 C1 R2;JosephsonCircuits.calcvaluetype([:R,:C,:R],[R1,C1,R2],[:R])
SymbolicUtils.BasicSymbolic{Number}[]
```
"""
function calcvaluetype(componenttypes::Vector{Symbol},componentvalues::Vector,
    components::Vector{Symbol};checkinverse::Bool=true)

    if length(componenttypes) !== length(componentvalues)
         throw(DimensionMismatch("componenttypes and componentvalues should have the same length"))
    end

    # use this to store the types we have seen so we don't call promote_type
    # or take the inverse for the same type more than once.
    typestoredict = Dict{DataType,Nothing}()

    componentsdict = Dict{Symbol,Nothing}()
    sizehint!(componentsdict,length(components))
    for component in components
        componentsdict[component] = nothing
    end

    # find the first one then break the loop so we have to execute the first
    # element logic only once. 
    valuetype = Nothing
    for (i,type) in enumerate(componenttypes)
        if haskey(componentsdict,type)
            valuetype = typeof(componentvalues[i])
            # add the original type to the typestore
            typestoredict[valuetype] = nothing
            if checkinverse
                valuetype = promote_type(typeof(1/componentvalues[i]),valuetype)
            end
            break
        end
    end

    for (i,type) in enumerate(componenttypes)
        if haskey(componentsdict,type)
            # if a different type is found, promote valuetype
            if typeof(componentvalues[i]) != valuetype
                # if it is a type we have seen before, do nothing
                valuetype = promote_type(typeof(componentvalues[i]),valuetype)
                if !haskey(typestoredict,valuetype)
                    typestoredict[valuetype] = nothing
                    if checkinverse
                        valuetype = promote_type(typeof(1/componentvalues[i]),valuetype)
                    end
                end
            end
        end
    end
    return Array{valuetype, 1}(undef, 0)
end

"""
    calcnoiseportimpedanceindices(componenttypes::Vector{Symbol},
        nodeindexarray::Matrix{Int}, mutualinductorbranchnames::Vector,
        componentvalues::Vector)

Find the resistors (not located at a port) or lossy capacitors or lossy
inductors and return their indices.

# Examples
```jldoctest
JosephsonCircuits.calcnoiseportimpedanceindices(
    [:R,:C,:Lj,:C],
    [2 2 3 3; 1 3 1 1],
    [],
    [50,5e-15,1e-12,30e-15])

# output
1-element Vector{Int64}:
 1
```
```jldoctest
JosephsonCircuits.calcnoiseportimpedanceindices(
    [:P,:R,:C,:Lj,:C],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,50,5e-15,1e-12,30e-15])

# output
Int64[]
```
```jldoctest
JosephsonCircuits.calcnoiseportimpedanceindices(
    [:R,:C,:Lj,:C],
    [2 2 3 3; 1 3 1 1],
    [],
    [50,5e-15,1e-12,(30+1im)*1e-15])

# output
2-element Vector{Int64}:
 1
 4
```
```jldoctest
JosephsonCircuits.calcnoiseportimpedanceindices(
    [:R,:C,:L,:C],
    [2 2 3 3; 1 3 1 1],
    [],
    [50,5e-15,(1+1im)*1e-12,30e-15])

# output
2-element Vector{Int64}:
 1
 3
```
"""
function calcnoiseportimpedanceindices(componenttypes::Vector{Symbol},
    nodeindexarray::Matrix{Int},mutualinductorbranchnames::Vector,
    componentvalues::Vector)

    if  length(componenttypes) != size(nodeindexarray,2) || length(componenttypes) != length(componentvalues)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end

    portimpedanceindices = calcportimpedanceindices(componenttypes,nodeindexarray,mutualinductorbranchnames,componentvalues)
    noiseportimpedanceindices = Int[]
    portindices = Int[]

    for i in eachindex(componenttypes)
        type=componenttypes[i]
        if type == :R
            if i in portimpedanceindices
                nothing
            else
                push!(noiseportimpedanceindices,i)
            end
        elseif type == :C && componentvalues[i] isa Complex
            if !iszero(componentvalues[i].im)
                push!(noiseportimpedanceindices,i)
            end
        elseif type == :L && componentvalues[i] isa Complex
            if !iszero(componentvalues[i].im)
                push!(noiseportimpedanceindices,i)
            end
        end
    end

    return noiseportimpedanceindices
end

"""
    calcportindicesnumbers(componenttypes::Vector{Symbol},
        nodeindexarray::Matrix{Int},mutualinductorbranchnames::Vector,
        componentvalues::Vector)

Return vectors containing the indices of the ports and their numbers.

# Examples
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:P,:R,:C,:Lj,:C],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,50,5e-15,1e-12,30e-15])

# output
([1], [1])
```
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:P,:R,:C,:Lj,:P],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,50,5e-15,1e-12,2])

# output
([1, 5], [1, 2])
```
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:P,:R,:C,:Lj,:P],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [2,50,5e-15,1e-12,1])

# output
([5, 1], [1, 2])
```
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:R,:C,:Lj,:C],
    [2 2 3 3; 1 3 1 1],
    [],
    [50,5e-15,1e-12,30e-15])

# output
(Int64[], Int64[])
```
"""
function calcportindicesnumbers(componenttypes::Vector{Symbol},
    nodeindexarray::Matrix{Int},mutualinductorbranchnames::Vector,
    componentvalues::Vector)

    if  length(componenttypes) != size(nodeindexarray,2) || length(componenttypes) != length(componentvalues)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end

    # define vectors to hold the numbers of the ports, the indices in the
    # componenttypes, componentvalues, and nodeindexarray, and the branches.
    # portnumbers = Int[]
    portindices = Int[]
    portbranches = Vector{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)}}(undef,0)

    # find the port indices
    for i in eachindex(componenttypes)
        type=componenttypes[i]
        if type == :P
            key = (nodeindexarray[1,i],nodeindexarray[2,i])
            keyreversed = (nodeindexarray[2,i],nodeindexarray[1,i])
            push!(portindices,i)
            # push!(portnumbers,componentvalues[i])
            push!(portbranches,key)
            push!(portbranches,keyreversed)
        end
    end

    # extract the port numbers. this causes runtime dispatch if the type of
    # componentvalues is any
    portnumbers = zeros(Int,length(portindices))
    for (i,index) in enumerate(portindices)
        portnumbers[i] = componentvalues[index]
    end

    # check that all of the port numbers are unique (two ports should not
    # have the same port number)
    if !allunique(portnumbers)
        throw(ArgumentError("Duplicate ports are not allowed."))
    end

    # check that the port branches are unique (two ports should not be defined
    # on the same branch)
    if !allunique(portbranches)
        throw(ArgumentError("Only one port allowed per branch."))
    end

    # sort by the portnumber
    sp = sortperm(portnumbers)

    return portindices[sp],portnumbers[sp]
end

"""
    calcportimpedanceindices(componenttypes::Vector{Symbol},
        nodeindexarray::Matrix{Int},mutualinductorbranchnames::Vector,
        componentvalues::Vector)

Find the resistors located at a port and return their indices.

# Examples
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:R,:C,:Lj,:C],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,50,5e-15,1e-12,30e-15])

# output
1-element Vector{Int64}:
 2
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:R,:C,:Lj,:C],
    [2 2 3 3; 1 3 1 1],
    [],
    [50,5e-15,1e-12,30e-15])

# output
Int64[]
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:R,:C,:Lj,:C,:P,:R],
    [2 3 2 3 3 3 2; 1 1 3 1 1 1 1],
    [],
    [1,50,5e-15,1e-12,30e-15,2,50.0])

# output
2-element Vector{Int64}:
 7
 2
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:R,:C,:Lj,:C,:P,:R],
    [2 2 2 3 3 3 3; 1 1 3 1 1 1 1],
    [],
    [1,50,5e-15,1e-12,30e-15,2,50.0])

# output
2-element Vector{Int64}:
 2
 7
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:R,:C,:Lj,:C,:P,:R],
    [2 2 2 3 3 3 3; 1 1 3 1 1 1 1],
    [],
    [2,50,5e-15,1e-12,30e-15,1,50.0])

# output
2-element Vector{Int64}:
 7
 2
```
"""
function calcportimpedanceindices(componenttypes::Vector{Symbol},
    nodeindexarray::Matrix{Int},mutualinductorbranchnames::Vector,
    componentvalues::Vector)

    if  length(componenttypes) != size(nodeindexarray,2) || length(componenttypes) != length(componentvalues)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end

    # calculate the indices and numbers of the ports
    portindices,portnumbers = calcportindicesnumbers(componenttypes, nodeindexarray,
        mutualinductorbranchnames, componentvalues)

    # make a dictionary so we can easily lookup the port number from the branch
    portbranchdict = Dict{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)},Int}()
    sizehint!(portbranchdict,2*length(portindices))
    for (i,portindex) in enumerate(portindices)
        key= (nodeindexarray[1,portindex],nodeindexarray[2,portindex])
        keyreversed = (nodeindexarray[2,portindex],nodeindexarray[1,portindex])
        portbranchdict[key] = portnumbers[i]
        portbranchdict[keyreversed] = portnumbers[i]
    end

    resistorbranchdict = Dict{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)},Nothing}()
    resistorindices = Int[]
    resistornumbers = Int[]

    # find the resistor associated with that port
    for i in eachindex(componenttypes)
        type=componenttypes[i]
        if type == :R
            key = (nodeindexarray[1,i],nodeindexarray[2,i])
            keyreversed = (nodeindexarray[2,i],nodeindexarray[1,i])

            if haskey(portbranchdict,key)
                if haskey(resistorbranchdict,key) || haskey(resistorbranchdict,keyreversed)
                    throw(ArgumentError("Only one resistor allowed per port."))
                else
                    push!(resistorindices,i)
                    push!(resistornumbers,portbranchdict[key])
                    resistorbranchdict[key] = nothing
                end
            end
        end
    end

    # throw an error if we haven't found a resistor for every port
    if length(resistorindices) != length(portindices)
        throw(ArgumentError("Ports without resistors detected. Each port must have a resistor to define the impedance."))
    end

    # sort by the portnumber
    sp = sortperm(resistornumbers)

    return resistorindices[sp]
end

"""
    componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)

Convert the array of component values to numbers, if defined in `circuitdefs`. 
This function is not type stable by design because we want the output array 
to use a concrete type if all of the values are evaluated to numbers. 

# Examples
```jldoctest
julia> JosephsonCircuits.componentvaluestonumber([:Lj1,:Lj2],Dict(:Lj1=>1e-12,:Lj2=>2e-12))
2-element Vector{Float64}:
 1.0e-12
 2.0e-12

julia> @variables Lj1 Lj2;JosephsonCircuits.componentvaluestonumber([Lj1,Lj1+Lj2],Dict(Lj1=>1e-12,Lj2=>2e-12))
2-element Vector{Float64}:
 1.0e-12
 3.0e-12
```
```jldoctest
# define a frequency dependent impedance function
Zfun(w,R) = ifelse(w>10,R,100*R);
# create symbolic variables including a two argument function
@variables w R
@register_symbolic Zfun(w,R)
# substitute in numerical values and functions for everything but w
out=JosephsonCircuits.componentvaluestonumber([R,Zfun(w,R)],Dict(R=>50));
println(out)
# evaluate with w = 2
println(JosephsonCircuits.Symbolics.substitute.(out,(Dict(w=>2),)))
# evaluate with w = 11
println(JosephsonCircuits.Symbolics.substitute.(out,(Dict(w=>11),)))

# output
Any[50, Zfun(w, 50)]
[50, 5000]
[50, 50]
```
"""
function componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)
    return [valuetonumber(v,circuitdefs) for v in componentvalues]
    # return map(valuetonumber,componentvalues,Base.Iterators.repeated(circuitdefs,length(componentvalues)))
end

"""
    valuetonumber(value::Symbol,circuitdefs)

If the component value is a symbol, assume it is a dictionary key.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(:Lj1,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0e-12
```
"""
function valuetonumber(value::Symbol,circuitdefs)
    return circuitdefs[value]
end

"""
    valuetonumber(value::String,circuitdefs)

If the component value is a string, assume it is a dictionary key.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber("Lj1",Dict("Lj1"=>1e-12,"Lj2"=>2e-12))
1.0e-12
```
"""
function valuetonumber(value::String,circuitdefs)
    return circuitdefs[value]
end

"""
    valuetonumber(value::Symbolics.Num,circuitdefs)

If the component value is Symbolics.Num, then try substituting in the definition
from `circuitdefs`.

# Examples
```jldoctest
julia> @variables Lj1;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
3.0e-12

julia> @variables Lj1 Lj2;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
4.0e-12
```
"""
function valuetonumber(value::Symbolics.Num,circuitdefs)
    # for Num types unwrap helps speed up
    # their evaluation and evalutes to a number. 
    return Symbolics.substitute(Symbolics.unwrap(value),circuitdefs)
end

"""
    valuetonumber(value::Complex{Symbolics.Num},circuitdefs)

If the component value `value` is Complex{Symbolics.Num}, then try substituting in the
definition from `circuitdefs`.

# Examples
```jldoctest
julia> @variables Lj1::Complex;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
3.0e-12

julia> @variables Lj1::Complex Lj2::Complex;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
ComplexTerm(real(Lj2) + real(Lj1) + im*(imag(Lj2) + imag(Lj1)))
```
"""
function valuetonumber(value::Complex{Symbolics.Num},circuitdefs)
    return Symbolics.unwrap(Symbolics.substitute(value,circuitdefs))
end

"""
    valuetonumber(value::Symbolics.Symbol, circuitdefs)

If the component value `value` has a type Complex{Symbolics.Num}, then try
substituting in the definition from `circuitdefs`.

# Examples
```jldoctest
julia> @syms Lj1;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
3.0e-12

julia> @syms Lj1 Lj2;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
4.0e-12
```
"""
function valuetonumber(value::Symbolics.Symbolic, circuitdefs)
    return Symbolics.substitute(value, circuitdefs)
end

"""
    valuetonumber(value, circuitdefs)

If the component value `value` is a number (or a type we haven't considered,
return it as is.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(1.0,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0
```
"""
function valuetonumber(value, circuitdefs)
    return value
end