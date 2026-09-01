# Everything here exists only to keep the legacy tuple netlist working. The
# typed `Circuit` is the input format; a netlist of `(name, node1, node2,
# value)` tuples is adapted into one here and then takes exactly the same
# path. Nothing outside this file should grow a dependency on it, so that
# when the legacy format is deprecated the file can be deleted whole.
#
# What that costs today: `LegacyNL`, which has no typed counterpart because
# component behavior is never inferred from an instance name outside this
# file; the name prefix table and the two functions that read it; and the
# port impedance convention where a port takes the resistor beside it.

# unwrap a wrapped symbolic value; the Symbolics extension adds the `Num`
# method
unwrapvalue(value) = value


# === legacy tuple netlist -> Circuit ===

"""
    LegacyNL(value)

An internal two terminal component holding the raw value of a legacy "NL"
netlist entry verbatim. The legacy parser accepts "NL" components but the
solver on this branch does not define their value semantics, so the adapter
preserves the value for exact round trips. Use
[`NonlinearInductor`](@ref) for nonlinear elements in the typed
representation.
"""
struct LegacyNL{T} <: AbstractComponent
    value::T
end
nterminals(::LegacyNL) = 2
# the only component type the legacy netlist has that the typed circuit
# does not, so it adds the one lowering method the compiler lacks
lowercomponent(def::LegacyNL, path) = :NL, def.value

const legacyallowedcomponents = ["Lj","NL","L","C","K","I","R","P"]

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
julia> JosephsonCircuits.parsecomponenttype("L10",["Lj","NL","L","C","K","I","R","P"])
3

julia> [JosephsonCircuits.parsecomponenttype(c,["Lj","NL","L","C","K","I","R","P"]) for c in ["Lj","NL","L","C","K","I","R","P"]]
8-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
 7
 8
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
                throw(ArgumentError(lazy"parsecomponenttype() currently only works for two letter components"))
            end
        end
    end
    throw(ArgumentError(lazy"No matching component found in allowedcomponents."))
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
julia> JosephsonCircuits.checkcomponenttypes(["Lj","NL","L","C","K","I","R","P"])
true
```
"""
function checkcomponenttypes(allowedcomponents::Vector{String})
    for i in eachindex(allowedcomponents)
        if i != parsecomponenttype(allowedcomponents[i],allowedcomponents)
            throw(ArgumentError(lazy"Allowed components parsing check has failed for $(allowedcomponents[i]). This can happen if a two letter long component comes after a one letter component. Please reorder allowedcomponents."))
        end
    end
    return true
end

function legacycomponent(typesymbol::Symbol, name, node1, node2, value)
    if typesymbol == :L
        return Inductor(value)
    elseif typesymbol == :C
        return Capacitor(value)
    elseif typesymbol == :R
        return Resistor(value)
    elseif typesymbol == :Lj
        return NonlinearInductor(value, sin, cos)
    elseif typesymbol == :NL
        return LegacyNL(value)
    elseif typesymbol == :I
        return CurrentSource(value)
    elseif typesymbol == :P
        return Port(legacyportnumber(name, value); termination = nothing)
    elseif typesymbol == :K
        return MutualInductor(value, string(node1), string(node2))
    else
        throw(ArgumentError(lazy"Unknown legacy component type $(typesymbol) for $(name)."))
    end
end

function legacyportnumber(name, value)
    v = unwrapvalue(value)
    if v isa Integer
        return Int(v)
    elseif v isa Real && isinteger(v)
        return Int(v)
    elseif v isa Complex && isreal(v) && isinteger(real(v))
        return Int(real(v))
    else
        throw(ArgumentError(lazy"The port $(name) has the value $(value), which cannot be interpreted as an integer port number."))
    end
end

"""
    Circuit(netlist::AbstractVector{<:Tuple})
    Circuit(netlist::AbstractVector{<:Tuple}, circuitdefs)

Construct a typed [`Circuit`](@ref) from a legacy tuple netlist, where each
entry is `(name, node1, node2, value)` and the component type is inferred
from the name prefix: Lj, NL, L, C, K, I, R, and P. Component name
prefixes are interpreted only inside this
adapter; typed component models never infer behavior from instance names.

Node labels become net names, so lowering the result back with
[`compile`](@ref) reproduces the legacy parse exactly. When
`circuitdefs` is supplied, symbolic values are substituted with
[`valuetonumber`](@ref) during conversion; otherwise values pass through
verbatim and `circuitdefs` may be supplied to the analysis as usual.

# Examples
```jldoctest
julia> Circuit([("P1","1","0",1),("R1","1","0",50.0),("C1","1","0",1e-12)]) isa Circuit
true
```
"""
function Circuit(netlist::AbstractVector)
    # The element type is not restricted to a tuple: a netlist assembled by
    # pushing onto a Vector{Any}, which the original parser accepts, is a
    # netlist. legacycircuit checks each entry and says what is wrong with it.
    return legacycircuit(netlist, nothing)
end

function Circuit(netlist::AbstractVector, circuitdefs::AbstractDict)
    return legacycircuit(netlist, circuitdefs)
end

# A legacy netlist states no port reference impedance: it is the resistor the
# user placed across the port, which the legacy solver finds by looking for a
# resistor on the port's branch. That search is the one piece of geometric
# discovery which stays, because legacy syntax carries no role marker, and
# doing it here once means the typed circuit downstream carries an explicit
# reference impedance like any other. The port itself is constructed
# unterminated, so the netlist's own resistor remains its only load.
function legacyportimpedances!(components, netlist)
    resistorat = Dict{Tuple{String,String},Any}()  # node pair -> (value, name)
    for (i, entry) in enumerate(netlist)
        components[i].second isa Resistor || continue
        n1, n2 = string(entry[2]), string(entry[3])
        haskey(resistorat, (n1, n2)) && continue
        r = (value = components[i].second.R, name = components[i].first)
        resistorat[(n1, n2)] = r
        resistorat[(n2, n1)] = r
    end
    for (i, entry) in enumerate(netlist)
        p = components[i].second
        p isa Port || continue
        R = get(resistorat, (string(entry[2]), string(entry[3])), nothing)
        isnothing(R) && continue
        components[i] = components[i].first =>
            Port(p.number; Z0 = R.value,
                termination = LegacyTermination(R.name))
    end
    return components
end

function legacycircuit(netlist, circuitdefs)
    checkcomponenttypes(legacyallowedcomponents)
    components = Vector{Pair{String,Any}}(undef, length(netlist))
    nodeorder = String[]
    nodegroups = Dict{String,Vector{Any}}()
    for (i, entry) in enumerate(netlist)
        if length(entry) != 4
            throw(ArgumentError(lazy"The netlist entry $(entry) on line $(i) must be a tuple of (name, node1, node2, value)."))
        end
        name, node1, node2, value = entry
        if !(name isa AbstractString)
            throw(ArgumentError(lazy"The component name $(name) on line $(i) must be a string."))
        end
        if occursin('/', name)
            throw(ArgumentError(lazy"The component name $(name) on line $(i) contains the reserved hierarchical path separator \"/\"."))
        end
        typeindex = parsecomponenttype(String(name), legacyallowedcomponents)
        typesymbol = Symbol(legacyallowedcomponents[typeindex])
        if !isnothing(circuitdefs)
            value = valuetonumber(value, circuitdefs)
        end
        components[i] = String(name) => legacycomponent(typesymbol, name,
            node1, node2, value)
        if typesymbol != :K
            for (t, node) in enumerate((node1, node2))
                label = string(node)
                if occursin('/', label)
                    throw(ArgumentError(lazy"The node label $(label) on line $(i) contains the reserved hierarchical path separator \"/\"."))
                end
                group = get!(() -> (push!(nodeorder, label); Any[]),
                    nodegroups, label)
                push!(group, (String(name), t))
            end
        end
    end
    legacyportimpedances!(components, netlist)
    connections = Vector{Any}(undef, length(nodeorder))
    for (i, label) in enumerate(nodeorder)
        group = nodegroups[label]
        if label == "0"
            push!(group, Ground)
        end
        # vector endpoint groups: a tuple of thousands of endpoints (a large
        # ground net) would create a fresh NTuple type and force
        # recompilation against it
        connections[i] = Net(label, group)
    end
    return Circuit(components, connections, nothing)
end
