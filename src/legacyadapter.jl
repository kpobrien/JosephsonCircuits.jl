# The legacy tuple netlist.
#
# The typed `Circuit` is the input format of the package. A netlist of
# `(name, node1, node2, value)` tuples, the original input format, is
# converted into a `Circuit` here and then follows the same path as one.
# Everything specific to the tuple format lives in this file, so that the
# format can be removed by deleting the file: the `LegacyNL` component, the
# name prefix table with the two functions that read it, and the
# convention that a port's reference impedance is the resistor placed
# across it.

# Unwrap a wrapped symbolic value to whatever it holds. The Symbolics
# extension adds the method for `Num`; everything else is already unwrapped.
unwrapvalue(value) = value


# === legacy tuple netlist -> Circuit ===

"""
    LegacyNL(value)

An internal two terminal component holding the value of a legacy `"NL"`
netlist entry verbatim. The tuple format accepts `NL` components and they
are carried through compilation as the `:NL` type, but the solvers do not
define what their value means. In the typed representation use
[`NonlinearInductor`](@ref).
"""
struct LegacyNL{T} <: AbstractComponent
    value::T
end
nterminals(::LegacyNL) = 2
# the one component the tuple format has and the typed circuit does not,
# with the lowering method the compiler would otherwise lack
lowercomponent(def::LegacyNL, path) = :NL, def.value

# The component type prefixes of the tuple format. Two letter prefixes must
# come before one letter prefixes with the same first letter; see
# `checkcomponenttypes`.
const legacyallowedcomponents = ["Lj","NL","L","C","K","I","R","P"]

"""
    parsecomponenttype(name::String,allowedcomponents::Vector{String})

The index in `allowedcomponents` of the one or two letter prefix which
matches the start of the component name `name`. Prefixes are tried in
order and the first match wins, so a two letter prefix listed after a one
letter prefix with the same first letter can never match;
[`checkcomponenttypes`](@ref) detects that ordering mistake.

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

Check that [`parsecomponenttype`](@ref) maps each prefix in
`allowedcomponents` back to its own index, and throw an `ArgumentError`
otherwise. This fails when a two letter prefix is listed after a one letter
prefix with the same first letter, which would shadow it.

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

# the typed component model of one tuple netlist entry
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

# the port number of a `P` entry, whose value must be an integer however it
# is written (an Int, a whole Float64, a whole real Complex)
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
    Circuit(netlist::AbstractVector)
    Circuit(netlist::AbstractVector, circuitdefs)

Construct a typed [`Circuit`](@ref) from a legacy tuple netlist. Each entry
is `(name, node1, node2, value)` and the component type is taken from the
prefix of `name`: `Lj` (Josephson junction), `NL`, `L`, `C`, `K` (mutual
inductor, whose "nodes" are the two inductor names), `I`, `R`, and `P`
(port, whose value is the port number). Only this adapter reads a name
prefix; typed component models never infer behavior from an instance name.

Node labels become net names, so [`compile`](@ref) of the result gives the
same tables the tuple netlist always produced. A port's reference
impedance is the value of the single resistor placed across it, which the
adapter records as the port's [`LegacyTermination`](@ref). When
`circuitdefs` is given, values are resolved with [`valuetonumber`](@ref)
during conversion; otherwise they pass through unchanged and
`circuitdefs` is given to the analysis as usual.

# Examples
```jldoctest
julia> Circuit([("P1","1","0",1),("R1","1","0",50.0),("C1","1","0",1e-12)]) isa Circuit
true
```
"""
function Circuit(netlist::AbstractVector)
    # The element type is not restricted to `Tuple`: a netlist built by
    # pushing onto a `Vector{Any}` is accepted, and `legacycircuit` checks
    # each entry and reports what is wrong with it.
    return legacycircuit(netlist, nothing)
end

function Circuit(netlist::AbstractVector, circuitdefs::AbstractDict)
    return legacycircuit(netlist, circuitdefs)
end

# A tuple netlist has no syntax for a port's reference impedance: by
# convention it is the one resistor placed across the port. This finds that
# resistor for every port, once, and rewrites the port with its value as
# `Z0` and the resistor as its `LegacyTermination`, so that downstream
# nothing needs to look for a resistor on a port's branch. The port is
# constructed without a matched termination of its own, so the netlist's
# resistor remains its only load.
function legacyportimpedances!(components, netlist)
    # Collect every resistor on each branch rather than the first one: a
    # port with two resistors across it has always been an error in this
    # format, and picking one silently would change the meaning of such a
    # netlist.
    resistorsat = Dict{Tuple{String,String},Vector{Any}}()
    for (i, entry) in enumerate(netlist)
        components[i].second isa Resistor || continue
        n1, n2 = string(entry[2]), string(entry[3])
        r = (value = components[i].second.R, name = components[i].first)
        push!(get!(Vector{Any}, resistorsat, (n1, n2)), r)
        n1 == n2 || push!(get!(Vector{Any}, resistorsat, (n2, n1)), r)
    end
    for (i, entry) in enumerate(netlist)
        p = components[i].second
        p isa Port || continue
        rs = get(resistorsat, (string(entry[2]), string(entry[3])), nothing)
        # a port with no resistor across it has no reference impedance
        if isnothing(rs)
            throw(ArgumentError(lazy"Ports without resistors detected. Each port must have a resistor to define the impedance. Port $(p.number) has none; place a resistor across it, or write the circuit in the typed format, where a port states its own reference impedance."))
        end
        if length(rs) > 1
            names = join([r.name for r in rs], ", ")
            throw(ArgumentError(lazy"Only one resistor allowed per port. Port $(p.number) has $(length(rs)) resistors across it ($(names)), and a legacy netlist has no way to say which one is its environment. Give the port a single resistor, or write the circuit in the typed format, where a port states its own termination and any number of device resistors may share its terminals."))
        end
        R = only(rs)
        components[i] = components[i].first =>
            Port(p.number; Z0 = R.value,
                termination = LegacyTermination(R.name))
    end
    return components
end

# Convert the tuple netlist to components and connection groups. Each
# distinct node label becomes one `Net` holding every terminal on it, in
# order of first appearance, with `Ground` appended to net "0".
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
        # the endpoints are a vector, not a tuple: a large ground net as an
        # `NTuple` of thousands of endpoints would be a new type to compile
        # against
        connections[i] = Net(label, group)
    end
    return Circuit(components, connections, nothing)
end
