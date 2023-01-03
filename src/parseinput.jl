
"""
    ParsedSortedCircuit(nodeindexarraysorted::Matrix{Int},
        uniquenodevectorsorted::Vector{String},
        mutualinductorvector::Vector{String},namevector::Vector{String},
        typevector::Vector{Symbol},valuevector::Vector,
        namedict::Dict{String, Int},Nnodes::Int)

A simple structure to hold the parsed and sorted circuit. See also 
[`parsesortcircuit`](@ref), [`parsecircuit`](@ref), and [`sortnodes`](@ref) 
for more explanation.

# Fields
- `nodeindexarraysorted::Matrix{Int}`: sorted array of node indices (where
    the length of the first axis is 2).
- `uniquenodevectorsorted::Vector{String}`: the sorted unique node names.
- `mutualinductorvector::Vector{String}`: the inductors coupled by the mutual
    inductors.
- `namevector::Vector{String}`: component names.
- `typevector::Vector{Symbol}`: the component (electrical engineering) types 
- `valuevector::Vector`: the component values
- `namedict::Dict{String, Int}`: names of the components as keys and the
    index at which the component occurs as the value.
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
    nodeindexarraysorted::Matrix{Int}
    uniquenodevectorsorted::Vector{String}
    mutualinductorvector::Vector{String}
    namevector::Vector{String}
    typevector::Vector{Symbol}
    valuevector::Vector
    namedict::Dict{String, Int}
    Nnodes::Int
end

"""
    parsesortcircuit(circuit;sorting=:name)

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
function parsesortcircuit(circuit;sorting=:number)

    # parse the circuit components
    parsedcircuit = parsecircuit(circuit)

    # sort the nodes
    uniquenodevectorsorted,nodeindexarraysorted = sortnodes(
        parsedcircuit.uniquenodevector,
        parsedcircuit.nodeindexvector,
        sorting=sorting)

    return ParsedSortedCircuit(
        nodeindexarraysorted,
        uniquenodevectorsorted,
        parsedcircuit.mutualinductorvector,
        parsedcircuit.namevector,
        parsedcircuit.typevector,
        parsedcircuit.valuevector,
        parsedcircuit.namedict,
        parsedcircuit.Nnodes)
end


"""
    ParsedCircuit(nodeindexvector::Vector{Int},
        uniquenodevector::Vector{String},mutualinductorvector::Vector{String},
        namevector::Vector{String},typevector::Vector{Symbol},
        valuevector::Vector,namedict::Dict{String, Int},Nnodes::Int)

A simple structure to hold the parsed circuit including a vector of node indices,
the unique node names, the inductors coupled by the mutual inductors, the
component names, the component types, the values of the components, a
dictionary of the names of the components as keys and the index at which the
component occurs as the value, and dictionaries for the ports and resistors
where the pair of nodes is the key and value is the component value.

See also [`parsecircuit`](@ref).

# Fields
- `nodeindexvector::Vector{Int}`: sorted vector of node indices where the
    two nodes for each component occur as consecutive elements (pairs).
- `uniquenodevector::Vector{String}`: the unique node names.
- `mutualinductorvector::Vector{String}`: the inductors coupled by the mutual
    inductors.
- `namevector::Vector{String}`: component names.
- `typevector::Vector{Symbol}`: the component (electrical engineering) types
- `valuevector::Vector`: the component values
- `namedict::Dict{String, Int}`: names of the components as keys and the
    index at which the component occurs as the value.
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
    mutualinductorvector::Vector{String} 
    namevector::Vector{String}
    typevector::Vector{Symbol}
    valuevector::Vector
    namedict::Dict{String, Int}
    Nnodes::Int
end


"""
    parsecircuit(circuit)

Parse `circuit` which is a vector where each element contains a tuple with the
component name, the first node, the second node, and the component value.
Component values can be numbers, symbols, or symbolic variables (including
symbolic functions). Definitions for component values are stored in
`circuitdefs`.

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
```jldoctest
circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",:Ipump))
push!(circuit,("R1","1","0",:Rleft))
push!(circuit,("L1","1","0",:L1))
push!(circuit,("L1","1","0",:L1))
push!(circuit,("K1","L1","L2",:K1))
push!(circuit,("L2","2","0",:L2))
push!(circuit,("C2","2","0",:C2))
parsecircuit(circuit)

# output
ERROR: ArgumentError: Name "L1" on line 5 is not unique.
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
    namedict = Dict{String,Int}()
    sizehint!(namedict,length(circuit))

    # vector to hold the component names
    namevector = Vector{String}(undef,length(circuit))

    # vector to hold the component type symbols, eg :L, :Lj, etc. 
    typevector = Vector{Symbol}(undef,length(circuit))

    # vector to hold the component values. take the type from the circuit array.
    valuevector = try
        Vector{fieldtype(eltype(circuit),4)}(undef,length(circuit))
    catch
        Vector{Any}(undef,length(circuit))
    end

    # vector to hold the nodes
    nodeindexvector = Vector{Int}(undef,2*length(circuit))

    # vector to hold the inductors which the mutual inductance connects
    mutualinductorvector = Array{String,1}(undef,0)

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

        if haskey(namedict,name)
           throw(ArgumentError("Name \"$(name)\" on line $(i) is not unique."))
        end

        # add the name to the dictionary
        namedict[name] = i

        # store the component name, type, and value
        namevector[i] = name
        typevector[i] = componenttype
        valuevector[i] = value

        # mutual inductors are treated differently because their "nodes"
        # refer to inductors rather than actual nodes. add zeros to
        # nodeindexvector whenever we find an inductor in order
        # to keep their length the same as the other arrays.

        if componenttype == :K
            push!(mutualinductorvector,node1)
            push!(mutualinductorvector,node2)
            nodeindexvector[2*i-1] = 0
            nodeindexvector[2*i] = 0
        else
            nodeindexvector[2*i-1] = processnode(uniquenodedict,uniquenodevector,node1)
            nodeindexvector[2*i] = processnode(uniquenodedict,uniquenodevector,node2)
        end

    end

    return ParsedCircuit(nodeindexvector,uniquenodevector,mutualinductorvector,
        namevector, typevector, valuevector, namedict,length(uniquenodevector))
end


"""
    processnode(uniquenodedict,uniquenodevector,node::String)

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
    processnode(uniquenodedict,uniquenodevector,node)

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

julia> JosephsonCircuits.parsecomponenttype("BAD1",["Lj","BAD","L","C","K","I","R","P"])
ERROR: ArgumentError: parsecomponenttype() currently only works for two letter components

julia> JosephsonCircuits.parsecomponenttype("B1",["Lj","L","C","K","I","R","P"])
ERROR: ArgumentError: No matching component found in allowedcomponents.
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
    checkcomponenttypes(allowedcomponents)

Check that each element in allowedcomponents is found at the correct place.
This will detect the case where a two letter component appears in 
allowedcomponents after a one letter component with the same starting letter.
parsecomponenttype() will match on the first value and this function will
throw an error.

# Examples
```jldoctest
julia> JosephsonCircuits.checkcomponenttypes(["Lj","L","C","K","I","R","P"])
true

julia> JosephsonCircuits.checkcomponenttypes(["L","Lj","C","K","I","R","P"])
ERROR: ArgumentError: Allowed components parsing check has failed for Lj. This can happen if a two letter long component comes after a one letter component. Please reorder allowedcomponents.
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
    extractbranches(typevector::Vector{Symbol},nodeindexarray::Matrix{Int})

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

julia> JosephsonCircuits.extractbranches([:P,:I,:R,:C,:Lj,:C],[2 2 2 2 3; 1 1 1 3 1])
ERROR: DimensionMismatch: typevector must have the same length as the number of node indices

julia> JosephsonCircuits.extractbranches([:P,:I,:R,:C,:Lj,:C],[2 2 2 2 3 3; 1 1 1 3 1 1; 0 0 0 0 0 0])
ERROR: DimensionMismatch: the length of the first axis must be 2

julia> JosephsonCircuits.extractbranches!([1],[:P,:I,:R,:C,:Lj,:C],[2 2 2 2 3 3; 1 1 1 3 1 1])
ERROR: DimensionMismatch: branchvector should be length zero
```
"""
function extractbranches(typevector::Vector{Symbol},nodeindexarray::Matrix{Int})

    branchvector = Array{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)},1}(undef,0)
    extractbranches!(branchvector,typevector,nodeindexarray)

    return branchvector
end

"""
    extractbranches!(branchvector,typevector,nodeindexarray)

Append tuples consisting of a pair of node indices (branches) which we will
use to calculate the incidence matrix.  Appends the tuples to branchvector.
"""
function extractbranches!(branchvector::Vector,typevector::Vector{Symbol},nodeindexarray::Matrix{Int})

    if  length(typevector) != size(nodeindexarray,2)
        throw(DimensionMismatch("typevector must have the same length as the number of node indices"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("the length of the first axis must be 2"))
    end

    if length(branchvector) != 0
        throw(DimensionMismatch("branchvector should be length zero"))
    end

    componenttypes = [:Lj,:L,:I,:P,:V]
    for i in eachindex(typevector)
        type = typevector[i]
        if type in componenttypes
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
ERROR: ArgumentError: No ground node found in netlist.
```
"""
function findgroundnodeindex(uniquenodevector::Vector{String})

    # find the ground node. error if we don't find it.
    groundnodeindex = 0
    for i in eachindex(uniquenodevector)
        if uniquenodevector[i] == "0"
            groundnodeindex = i
            break
        end
    end

    if groundnodeindex == 0
        throw(ArgumentError("No ground node found in netlist."))
    end

    return groundnodeindex
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

julia> JosephsonCircuits.calcnodesorting(["30","11","0","2"];sorting=:test)
ERROR: ArgumentError: Unknown sorting type.

julia> JosephsonCircuits.calcnodesorting(["30","11","1","2"];sorting=:none)
ERROR: ArgumentError: No ground node found in netlist.
```
"""
function calcnodesorting(uniquenodevector::Vector{String};sorting=:number)

    # vector of indices for the sortperm
    uniquenodevectorsortindices = ones(Int,length(uniquenodevector))

    # sort according to the desired scheme
    if sorting == :name
        # sort the vector of unique node strings
        sortperm!(uniquenodevectorsortindices,uniquenodevector)

    elseif sorting == :number
        # convert the unique node strings to integers and sort those
        sortperm!(uniquenodevectorsortindices,parse.(Int,uniquenodevector))

    elseif sorting == :none
        # don't perform any sorting. keep the nodes in
        # order of first appearance, except move the
        # ground node first later as always
        uniquenodevectorsortindices .= 1:length(uniquenodevectorsortindices)

    else
        throw(ArgumentError("Unknown sorting type."))
    end

    # find the ground node. error if we don't find it.
    groundnodeindex = findgroundnodeindex(uniquenodevector)

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
    sortnodes(uniquenodevector,nodeindexvector;sorting=:name)

Sort the unique node names in `uniquenodevector` according to the specified
sorting scheme, always placing the ground node at the beginning.

Return the sorted `uniquenodevector` and `nodeindexvector` (with the vector reshaped
from a vector of length 2*Nnodes into a matrix with dimensions 2 by Nnodes).

# Keywords
- `sorting = :name`: Sort the vector of strings. This always works but leads
    to results like "101" comes before "11".
- `sorting = :number`: Convert the node strings to integer and sort by these
    (this errors if the nodes names cannot be converted to integers).
- `sorting = :none`: Don't perform any sorting except to place the ground node
    first.

# Examples
```jldoctest
julia> uniquenodevectorsorted,nodeindexarray=JosephsonCircuits.sortnodes(["101","0","111","11"],[1,2,1,2,1,2,1,3,3,2,3,2,4,1],sorting=:none);println(uniquenodevectorsorted);println(nodeindexarray);
["0", "101", "111", "11"]
[2 2 2 2 3 3 4; 1 1 1 3 1 1 2]

julia> uniquenodevectorsorted,nodeindexarray=JosephsonCircuits.sortnodes(["101","0","111","11"],[1,2,1,2,1,2,1,3,3,2,3,2,4,1],sorting=:name);println(uniquenodevectorsorted);println(nodeindexarray);
["0", "101", "11", "111"]
[2 2 2 2 4 4 3; 1 1 1 4 1 1 2]

julia> uniquenodevectorsorted,nodeindexarray=JosephsonCircuits.sortnodes(["101","0","111","11"],[1,2,1,2,1,2,1,3,3,2,3,2,4,1],sorting=:number);println(uniquenodevectorsorted);println(nodeindexarray);
["0", "11", "101", "111"]
[3 3 3 3 4 4 2; 1 1 1 4 1 1 3]

julia> uniquenodevectorsorted,nodeindexarray=JosephsonCircuits.sortnodes(["1", "0", "2"],[1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2],sorting=:number);println(uniquenodevectorsorted);println(nodeindexarray);
["0", "1", "2"]
[2 2 2 2 0 3 3; 1 1 1 1 0 1 1]
```
"""
function sortnodes(uniquenodevector::Vector{String},nodeindexvector::Vector{Int};sorting=:name)

    nodeindexarraysorted = zeros(eltype(nodeindexvector),2,length(nodeindexvector)÷2)

    uniquenodevectorsortindices=calcnodesorting(uniquenodevector;sorting=sorting)

    nodevectorsortindices = indexin(1:length(uniquenodevectorsortindices),uniquenodevectorsortindices)

    # for (i,j) in enumerate(eachindex(nodeindexvector))
    for (i,j) in enumerate(nodeindexvector)
        # if it's a mutual inductor the node index will be zero because the
        # mutual inductor is between two inductors not between two nodes.
        # it not a mutual inductor, assign the sorted node index.
        if j == 0
            nothing
        else
            nodeindexarraysorted[i] = nodevectorsortindices[j]
        end
    end

    return uniquenodevector[uniquenodevectorsortindices],nodeindexarraysorted
end


"""
    calcvaluetype(typevector::Vector{Symbol},valuevector::Vector,
        components::Vector{Symbol};checkinverse=true)

Returns a zero length vector with the (computer science) type which will hold
a set of circuit components of the (electrical engineering) types given in
`components`.

# Arguments
- `typevector::Vector{Symbol}`: the component (electrical engineering) types.
- `valuevector::Vector`: the component values.
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
SymbolicUtils.Symbolic{Number}[]

julia> uniquenodevectorsorted,nodeindexarray=JosephsonCircuits.sortnodes(["1", "0", "2"],[1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2],sorting=:number);println(uniquenodevectorsorted);println(nodeindexarray);
["0", "1", "2"]
[2 2 2 2 0 3 3; 1 1 1 1 0 1 1]

julia> JosephsonCircuits.calcvaluetype([:C,:R],[1,2,3],[:R])
ERROR: DimensionMismatch: typevector and valuevector should have the same length
```
"""
function calcvaluetype(typevector::Vector{Symbol},valuevector::Vector,
    components::Vector{Symbol};checkinverse=true)

    if length(typevector) !== length(valuevector)
         throw(DimensionMismatch("typevector and valuevector should have the same length"))
    end

    # use this to store the types we have seen so we don't call promote_type
    # or take the inverse for the same type more than once.
    typestore = Vector{DataType}(undef,0)

    valuetype = Nothing

    # find the first one then break the loop so we have to execute the first
    # element logic only once. 
    for (i,type) in enumerate(typevector)
        if type in components

            valuetype = typeof(valuevector[i])
            # add the original type to the typestore
            push!(typestore,valuetype)
            if checkinverse
                valuetype = promote_type(typeof(1/valuevector[i]),valuetype)
            end
            break
        end
    end

    for (i,type) in enumerate(typevector)
        if type in components

            # if a different type is found, promote valuetype
            if typeof(valuevector[i]) != valuetype
                # if it is a type we have seen before, do nothing
                valuetype = promote_type(typeof(valuevector[i]),valuetype)
                if !(valuetype in typestore)
                    push!(typestore,valuetype)
                    if checkinverse
                        valuetype = promote_type(typeof(1/valuevector[i]),valuetype)
                    end
                end
            end
        end
    end
    return Array{valuetype, 1}(undef, 0)
end

"""
    calcnoiseportimpedanceindices(typevector::Vector{Symbol},nodeindexarray::Matrix{Int},
        mutualinductorvector::Vector,valuevector::Vector)

Find the resistors (not located at a port) or lossy capacitors or lossy inductors 
and return their indices. 

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
```jldoctest
julia> JosephsonCircuits.calcnoiseportimpedanceindices([:R,:C,:Lj],[2 2 3 3; 1 3 1 1],[],[50,5e-15,1e-12,30e-15])
ERROR: DimensionMismatch: Input arrays must have the same length

julia> JosephsonCircuits.calcnoiseportimpedanceindices([:R,:C,:Lj,:C],[2 2 3 3; 1 3 1 1; 0 0 0 0],[],[50,5e-15,1e-12,30e-15])
ERROR: DimensionMismatch: The length of the first axis must be 2
```
"""
function calcnoiseportimpedanceindices(typevector::Vector{Symbol},
    nodeindexarray::Matrix{Int},mutualinductorvector::Vector,
    valuevector::Vector)

    if  length(typevector) != size(nodeindexarray,2) || length(typevector) != length(valuevector)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end

    portimpedanceindices = calcportimpedanceindices(typevector,nodeindexarray,mutualinductorvector,valuevector)
    out = Int[]
    portindices = Int[]

    for i in eachindex(typevector)
        type=typevector[i]
        if type == :R
            if i in portimpedanceindices
                nothing
            else
                push!(out,i)
            end
        elseif type == :C && valuevector[i] isa Complex
            if !iszero(valuevector[i].im)
                push!(out,i)
            end
        elseif type == :L && valuevector[i] isa Complex
            if !iszero(valuevector[i].im)
                push!(out,i)
            end
        end
    end

    return out
end

"""
    calcportindicesnumbers(typevector::Vector{Symbol},
        nodeindexarray::Matrix{Int},mutualinductorvector::Vector,
        valuevector::Vector)

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
    [:R,:C,:Lj,:C],
    [2 2 3 3; 1 3 1 1],
    [],
    [50,5e-15,1e-12,30e-15])

# output
(Int64[], Int64[])
```
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:P,:R,:C,:Lj,:C],
    [2 2 2 3 3; 1 1 3 1 1; 0 0 0 0 0],
    [],
    [1,50,5e-15,1e-12,30e-15])

# output
ERROR: DimensionMismatch: The length of the first axis must be 2
```
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:P,:R,:C,:Lj,:C],
    [2 2 2 3; 1 1 3 1],
    [],
    [1,50,5e-15,1e-12,30e-15])

# output
ERROR: DimensionMismatch: Input arrays must have the same length
```
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:P,:P,:C,:Lj,:C],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,2,5e-15,1e-12,30e-15])

# output
ERROR: Only one port allowed per branch.
```
```jldoctest
JosephsonCircuits.calcportindicesnumbers(
    [:P,:R,:C,:Lj,:P],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,50,5e-15,1e-12,1])

# output
ERROR: Duplicate ports are not allowed.
```
"""
function calcportindicesnumbers(typevector::Vector{Symbol},
    nodeindexarray::Matrix{Int},mutualinductorvector::Vector,
    valuevector::Vector)

    if  length(typevector) != size(nodeindexarray,2) || length(typevector) != length(valuevector)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end


    portarray = Vector{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)}}(undef,0)
    # portnumbers = Vector{eltype(valuevector)}(undef,0)
    portnumbers = Vector{Int}(undef,0)
    portindices = Int[]

    # find the ports
    for i in eachindex(typevector)
        type=typevector[i]
        if type == :P
            key= (nodeindexarray[1,i],nodeindexarray[2,i])

            if valuevector[i] in portnumbers
                error("Duplicate ports are not allowed.")
            elseif key in portarray
                error("Only one port allowed per branch.")
            else
                push!(portindices,i)
                push!(portnumbers,valuevector[i])
                push!(portarray,key)
            end
        end
    end

    # sort by the portnumber
    sp = sortperm(portnumbers)

    return portindices[sp],portnumbers[sp]
end


"""
    calcportimpedanceindices(typevector::Vector{Symbol},
        nodeindexarray::Matrix{Int},mutualinductorvector::Vector,
        valuevector::Vector)

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
    [:P,:R,:C,:Lj,:C],
    [2 2 2 3 3; 1 1 3 1 1; 0 0 0 0 0],
    [],
    [1,50,5e-15,1e-12,30e-15])

# output
ERROR: DimensionMismatch: The length of the first axis must be 2
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:R,:C,:Lj,:C],
    [2 2 2 3; 1 1 3 1],
    [],
    [1,50,5e-15,1e-12,30e-15])

# output
ERROR: DimensionMismatch: Input arrays must have the same length
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:P,:C,:Lj,:C],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,2,5e-15,1e-12,30e-15])

# output
ERROR: Only one port allowed per branch.
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:R,:C,:Lj,:P],
    [2 2 2 3 3; 1 1 3 1 1],
    [],
    [1,50,5e-15,1e-12,1])

# output
ERROR: Duplicate ports are not allowed.
```
```jldoctest
JosephsonCircuits.calcportimpedanceindices(
    [:P,:R,:R,:Lj,:C],
    [2 2 2 2 3; 1 1 1 3 1],
    [],
    [1,50.0,50.0,1e-12,30e-15])

# output
ERROR: Only one resistor allowed per port.
```
"""
function calcportimpedanceindices(typevector::Vector{Symbol},
    nodeindexarray::Matrix{Int},mutualinductorvector::Vector,
    valuevector::Vector)

    if  length(typevector) != size(nodeindexarray,2) || length(typevector) != length(valuevector)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end

    portarray = Vector{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)}}(undef,0)
    # portnumbers = Vector{eltype(valuevector)}(undef,0)
    portnumbers = Vector{Int}(undef,0)
    portindices = Int[]

    # find the ports
    for i in eachindex(typevector)
        type=typevector[i]
        if type == :P
            key= (nodeindexarray[1,i],nodeindexarray[2,i])

            if valuevector[i] in portnumbers
                error("Duplicate ports are not allowed.")
            elseif key in portarray
                error("Only one port allowed per branch.")
            else
                push!(portindices,i)
                push!(portnumbers,valuevector[i])
                push!(portarray,key)
            end
        end
    end

    resistorarray = Vector{Tuple{eltype(nodeindexarray),eltype(nodeindexarray)}}(undef,0)
    resistorindices = Int[]

    # find the resistor associated with that port
    for i in eachindex(typevector)
        type=typevector[i]
        if type == :R
            key= (nodeindexarray[1,i],nodeindexarray[2,i])
            if key in portarray
                if key in resistorarray
                    error("Only one resistor allowed per port.")
                else
                    push!(resistorindices,i)
                    push!(resistorarray,key)
                end
            end
        end
    end

    # sort by the portnumber
    sp = sortperm(portnumbers)

    return resistorindices[sp]
end


"""
    valuevectortonumber(valuevector,circuitdefs)

Convert the array of component values to numbers, if defined in `circuitdefs`. 
This function is not type stable by design because we want the output array 
to use a concrete type if all of the values are evaluated to numbers. 

# Examples
```jldoctest
julia> JosephsonCircuits.valuevectortonumber([:Lj1,:Lj2],Dict(:Lj1=>1e-12,:Lj2=>2e-12))
2-element Vector{Float64}:
 1.0e-12
 2.0e-12

julia> @variables Lj1 Lj2;JosephsonCircuits.valuevectortonumber([Lj1,Lj1+Lj2],Dict(Lj1=>1e-12,Lj2=>2e-12))
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
out=JosephsonCircuits.valuevectortonumber([R,Zfun(w,R)],Dict(R=>50));
println(out)
# evaluate with w = 2
println(JosephsonCircuits.substitute.(out,(Dict(w=>2),)))
# evaluate with w = 11
println(JosephsonCircuits.substitute.(out,(Dict(w=>11),)))

# output
Any[50, Zfun(w, 50)]
[50, 5000]
[50, 50]
```
"""
function valuevectortonumber(valuevector::Vector,circuitdefs::Dict)
    return [valuetonumber(v,circuitdefs) for v in valuevector]
    # return map(valuetonumber,valuevector,Base.Iterators.repeated(circuitdefs,length(valuevector)))
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

If the component value is Complex{Symbolics.Num}, then try substituting in the
definition from `circuitdefs`. 

NOTE: Below fails because Symbolics.jl doesn't have good support for complex 
numbers. 
@variables Lj1::Complex Lj2::Complex;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
4.0e-12

# Examples
```jldoctest
julia> @variables Lj1::Complex;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
3.0e-12
```
"""
function valuetonumber(value::Complex{Symbolics.Num},circuitdefs)
    return Symbolics.unwrap(Symbolics.substitute(value,circuitdefs))
end


"""
    valuetonumber(value::Symbolics.Symbol,circuitdefs)

If the component value is Complex{Symbolics.Num}, then try substituting in the
definition from `circuitdefs`. 

# Examples
```jldoctest
julia> @syms Lj1;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
3.0e-12

julia> @syms Lj1 Lj2;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
4.0e-12
```
"""
function valuetonumber(value::Symbolics.Symbolic,circuitdefs)
    return Symbolics.substitute(value,circuitdefs)
end


"""
    valuetonumber(value,circuitdefs)

If the component value is a number (or a type we haven't considered, return it
as is. 

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(1.0,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0
```
"""
function valuetonumber(value,circuitdefs)
    return value
end