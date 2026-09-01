# The input path, in the order a circuit travels it.
#
#   Circuit          what the user writes: components, connections, and the
#                    interface a subcircuit presents to its parent. Parsed
#                    one level at a time by `parsecircuitlevel`.
#   ElaboratedCircuit  the hierarchy flattened: every instance given a path,
#                    every net a single wire, by `elaborate`.
#   CompiledCircuit  flat integer indexed tables plus the groups the
#                    assembly reads, by `compile`.
#
# The components themselves are in `circuitmodel.jl`, which is a catalogue
# of what each one is rather than a step on this path. How a component
# value becomes a number is in `circuitvalue.jl`. The legacy tuple netlist
# is adapted into a `Circuit` in `legacyadapter.jl` and joins the path at
# the top.


# === the circuit a user writes ===


"""
    GroundType

The singleton type of [`Ground`](@ref).
"""
struct GroundType end

"""
    Ground

The distinguished global electrical reference. `Ground` may appear as an
endpoint in any connection group, and as the negative pin of an interface
port. The ground net is always named "0".

For uniformity with ordinary components, ground may also be declared in the
components list and referenced through its single terminal:

```julia
Circuit([:r1 => Resistor(50.0), :gnd => Ground()],
    [[(:r1, 1)], [(:r1, 2), (:gnd, 1)]])
```

`Ground()` and `Ground` are the same object, so both spellings work in both
positions. A ground instance is sugar for the reference net rather than a
device: every reference to its terminal resolves to the global ground net,
however many instances are declared and at whatever level of the hierarchy,
and it contributes no flattened component.
"""
const Ground = GroundType()

# `Ground()` constructs like a component model but is the same singleton, so
# the sentinel and component spellings cannot diverge.
(::GroundType)() = Ground

Base.show(io::IO, ::GroundType) = print(io, "Ground")

"""
    PortRef(instance, key)

An explicit reference to the bundled two terminal port `key` of the
component instance `instance`, for use in pair connections when a bare
`(instance, key)` tuple would be ambiguous between a scalar pin and a port.
"""
struct PortRef{I,K}
    instance::I
    key::K
end

"""
    PinRef(instance, key)

An explicit reference to the scalar pin or terminal `key` of the component
instance `instance`, for use when a bare `(instance, key)` tuple would be
ambiguous between a scalar pin and a port.
"""
struct PinRef{I,K}
    instance::I
    key::K
end

"""
    Net(name, endpoints)

A named connection group: all `endpoints` belong to the same electrical net,
which is given the name `name` for diagnostics and outputs. Unnamed nets are
named automatically. Names attached to the ground net are ignored; the
ground net is always named "0".

# Examples
```julia
Net(:bias, ((:source, 1), (:device, 2)))
```
"""
struct Net{N,E}
    name::N
    endpoints::E
end

"""
    Instance(definition)

An explicit instance wrapper around a component definition. `:id => model`
and `:id => Instance(model)` are equivalent. Keyword overrides (parameters,
thermal bindings) are reserved for future use and currently raise an error.
"""
struct Instance{D}
    definition::D
    function Instance(definition; kwargs...)
        if !isempty(kwargs)
            throw(ArgumentError(lazy"Instance overrides are not yet supported; got keywords $(keys(kwargs)). Remove the keywords or construct a separate definition."))
        end
        return new{typeof(definition)}(definition)
    end
end

"""
    Interface(; pins, ports = nothing)

The interface of a hierarchical circuit, exposing internal endpoints as
scalar pins and optionally grouping pins into oriented two terminal wave
port views.

`pins` maps external keys (integers or symbols) to internal scalar
endpoints, analogous to the pin list of a SPICE `.subckt`:

```julia
pins = [1 => (:jj1, 1), 2 => (:jj2, 2), 3 => (:cap, 2)]
```

`ports` optionally maps external port keys to `(positive, negative)` pairs
of pin keys, where the negative entry may be [`Ground`](@ref):

```julia
ports = [1 => (1, 3), 2 => (2, 3)]
```

The parent circuit physically binds pins; connecting port to port with pair
syntax is shorthand which expands to the pin connections.

The explicit call is optional: `Circuit(components, connections;
pins = ..., ports = ...)` constructs the same interface through keywords.
"""
struct Interface{P,W}
    pins::P
    ports::W
end
Interface(; pins, ports = nothing) = Interface(pins, ports)

"""
    Circuit(components, connections, interface = nothing;
        pins = nothing, ports = nothing, validate = true)

The public typed circuit representation.

- `components` associates unique instance identifiers (symbols, strings, or
  integers) with component models, typically as a vector of pairs.
- `connections` describes which component endpoints are electrically
  connected, as groups (tuples of endpoints on one net), pairs (port to port
  bonds), and [`Net`](@ref) entries.
- `interface` optionally exposes pins and ports so that the circuit can be
  used as a component inside another circuit.

The `pins` and `ports` keywords are sugar for the positional interface:
`Circuit(components, connections; pins = ..., ports = ...)` is
`Circuit(components, connections, Interface(pins = ..., ports = ...))`.
Give the interface one way or the other, not both.

The constructor validates identifiers, endpoint references, connector
namespaces, and the interface, so that errors point at the construction
site, but stores the collections exactly as given: no data is copied, and a
thousand instances of one subcircuit hold a thousand references to the same
object. Use [`elaborate`](@ref) to flatten the hierarchy.

# Endpoint grammar

| written | meaning |
|---|---|
| `(:inst, k)` | scalar terminal or pin `k` in a group; port `k` in a pair |
| `(:inst, p, t)` | terminal `t` (1 signal, 2 reference) of port `p` |
| `Ground` | the global reference net "0" |
| `(:gnd, 1)` with `:gnd => Ground()` | the same reference net, component style |
| `PortRef`/`PinRef` | explicit namespace selection |

In a group every endpoint is scalar. In a pair `a => b`, endpoints resolve
in the port namespace of components which expose ports, and the pair
expands to signal-to-signal and reference-to-reference groups; components
without ports fall back to scalar endpoints, making a pair of scalar
endpoints sugar for a two element group. A key which exists both as a pin
and as a port of a subcircuit is an error in a pair and requires `PortRef`
or `PinRef`.

Connection groups may be written as tuples or vectors of endpoints.
Vectors are recommended for large or generated groups: a vector is one
type whatever its length, where every distinct tuple shape is a separate
type for the compiler to specialize on.

# Examples
```julia
circuit = Circuit(
    [:l1 => Inductor(1e-9), :c1 => Capacitor(100e-15), :p1 => Port(1)],
    [[(:p1, 1), (:l1, 1)],
     [(:l1, 2), (:c1, 1)],
     [(:c1, 2), (:p1, 2), Ground]],
)
```
"""
struct Circuit{C,K,I} <: AbstractComponent
    components::C
    connections::K
    interface::I
    function Circuit(components, connections, interface = nothing;
            pins = nothing, ports = nothing, validate::Bool = true)
        if !isnothing(pins) || !isnothing(ports)
            if !isnothing(interface)
                throw(ArgumentError("Give the interface either positionally or through the pins/ports keywords, not both."))
            end
            if isnothing(pins)
                throw(ArgumentError("The ports keyword requires pins: ports group pin keys into two terminal views. Pass pins as well."))
            end
            interface = Interface(pins, ports)
        end
        if validate
            parsecircuitlevel(components, connections, interface)
        end
        return new{typeof(components),typeof(connections),
            typeof(interface)}(components, connections, interface)
    end
end

countelements(x) = Base.IteratorSize(x) isa Union{Base.HasLength,Base.HasShape} ?
    length(x) : count(Returns(true), x)

# compact display: the stored collections (and their types) can be large,
# especially for generated circuits, and nested subcircuits would otherwise
# print recursively
function Base.show(io::IO, c::Circuit)
    print(io, "Circuit(", countelements(c.components), " components, ",
        countelements(c.connections), " connections")
    if !isnothing(c.interface)
        print(io, ", ", countelements(c.interface.pins), " pins")
        if !isnothing(c.interface.ports)
            print(io, ", ", countelements(c.interface.ports), " ports")
        end
    end
    print(io, ")")
end

# === connector protocol ===

"""
    nterminals(component)

The number of scalar electrical terminals of a component model. Two
terminal lumped components have 2; a [`ScatteringParameters`](@ref) has two per
port; a hierarchical [`Circuit`](@ref) has one per interface pin; a
[`MutualInductor`](@ref) has none because it couples branches, not nets; a
[`Ground`](@ref) instance has one, which is the reference net itself.
"""
nterminals(::GroundType) = 1
nterminals(::Inductor) = 2
nterminals(::Capacitor) = 2
nterminals(::Resistor) = 2
nterminals(::CurrentSource) = 2
nterminals(::VoltageSource) = 2
nterminals(::Port) = 2
nterminals(::NonlinearInductor) = 2
nterminals(::MutualInductor) = 0
nterminals(c::ScatteringParameters) = 2*c.nports
nterminals(c::GaussianChannel) = 2*c.nmodes
function nterminals(c::Circuit)
    if isnothing(c.interface)
        throw(ArgumentError("A Circuit used as a component must have an Interface."))
    end
    return length(c.interface.pins)
end
nterminals(c) = throw(ArgumentError(lazy"$(typeof(c)) is not a known component model."))

"""
    hasports(component)

Whether the component exposes bundled two terminal port views addressable
in pair connections.
"""
hasports(c::ScatteringParameters) = true
hasports(c::GaussianChannel) = true
hasports(c::Circuit) = !isnothing(c.interface) && !isnothing(c.interface.ports)
hasports(c) = false

componentnports(c::ScatteringParameters) = c.nports
componentnports(c::GaussianChannel) = c.nmodes

isgrounded(c::ScatteringParameters) = c.grounded
isgrounded(c::GaussianChannel) = c.grounded

# === component table ===

struct ComponentTable
    ids::Vector{Any}
    defs::Vector{Any}
    index::Dict{Any,Int}
end

function componenttable(components)
    ids = Vector{Any}()
    defs = Vector{Any}()
    index = Dict{Any,Int}()
    for entry in components
        if !(entry isa Pair)
            throw(ArgumentError(lazy"Each element of components must be a Pair of an identifier and a component model; got $(typeof(entry))."))
        end
        id = entry.first
        def = entry.second
        if !(id isa Symbol || id isa AbstractString || id isa Integer)
            throw(ArgumentError(lazy"Instance identifiers must be symbols, strings, or integers; got $(typeof(id)) for $(id)."))
        end
        if occursin('/', string(id))
            throw(ArgumentError(lazy"Instance identifier $(id) contains the reserved hierarchical path separator \"/\"."))
        end
        if def isa Instance
            def = def.definition
        end
        if def isa Circuit && isnothing(def.interface)
            throw(ArgumentError(lazy"The Circuit used as instance $(id) has no Interface. A subcircuit must expose pins through an Interface."))
        end
        if haskey(index, id)
            throw(ArgumentError(lazy"Instance identifier $(id) is not unique."))
        end
        push!(ids, id)
        push!(defs, def)
        index[id] = length(ids)
    end
    return ComponentTable(ids, defs, index)
end

function instanceindex(table::ComponentTable, id, context::AbstractString)
    i = get(table.index, id, 0)
    if i == 0
        throw(ArgumentError(lazy"The endpoint $(context) references the instance $(id), which does not exist in this circuit."))
    end
    return i
end

# === scalar endpoint resolution ===

# Resolve a scalar key on a definition to a terminal index in 1:nterminals.
function scalarterminal(def, id, k)
    n = nterminals(def)
    if !(k isa Integer)
        throw(ArgumentError(lazy"Terminal keys of $(id) must be integers; got $(k)."))
    end
    if !(1 <= k <= n)
        throw(ArgumentError(lazy"The instance $(id) has terminals 1:$(n); got terminal $(k)."))
    end
    return Int(k)
end

function scalarterminal(def::MutualInductor, id, k)
    throw(ArgumentError(lazy"The mutual inductor $(id) couples two inductor branches and has no terminals; it must not appear in connections."))
end

function scalarterminal(def::Union{ScatteringParameters,GaussianChannel}, id, k)
    if isgrounded(def)
        if !(k isa Integer) || !(1 <= k <= componentnports(def))
            throw(ArgumentError(lazy"The grounded multiport $(id) has ports 1:$(componentnports(def)); got $(k)."))
        end
        return 2*(Int(k)-1) + 1 # signal terminal of port k
    end
    throw(ArgumentError(lazy"Port $(k) of the multiport $(id) is a two terminal port; write ($(repr(id)), $(k), 1) or ($(repr(id)), $(k), 2) for a single terminal, bond port to port with pair syntax, or construct the block with grounded = true to tie all reference terminals to Ground."))
end

function scalarterminal(def::Circuit, id, k)
    pins = def.interface.pins
    for (i, pin) in enumerate(pins)
        if isequal(pin.first, k)
            return i
        end
    end
    throw(ArgumentError(lazy"The subcircuit $(id) has no pin $(k). Its pins are $(collect(p.first for p in pins))."))
end

# Resolve a (port, terminal) scalar address.
function portterminal(def, id, p, t)
    throw(ArgumentError(lazy"The instance $(id) has no ports; address its terminals as ($(repr(id)), terminal)."))
end

function portterminal(def::Union{ScatteringParameters,GaussianChannel}, id, p, t)
    np = componentnports(def)
    if !(p isa Integer) || !(1 <= p <= np)
        throw(ArgumentError(lazy"The multiport $(id) has ports 1:$(np); got port $(p)."))
    end
    if !(t isa Integer) || !(1 <= t <= 2)
        throw(ArgumentError(lazy"Port terminals are 1 (signal) and 2 (reference); got $(t) on port $(p) of $(id)."))
    end
    if isgrounded(def) && t == 2
        throw(ArgumentError(lazy"The reference terminal of port $(p) of $(id) is auto-tied to Ground by grounded = true; remove this connection or set grounded = false."))
    end
    return 2*(Int(p)-1) + Int(t)
end

# === port view resolution ===

# Return the (positive, negative) scalar sides of a bundled port view.
# Each side is either (instanceindex, terminal) or Ground.
function portview(def, id, p)
    throw(ArgumentError(lazy"The instance $(id) exposes no ports, so $(p) cannot be used as a port in a pair connection."))
end

function portview(def::Union{ScatteringParameters,GaussianChannel}, id, p)
    np = componentnports(def)
    if !(p isa Integer) || !(1 <= p <= np)
        throw(ArgumentError(lazy"The multiport $(id) has ports 1:$(np); got port $(p)."))
    end
    return (2*(Int(p)-1) + 1, 2*(Int(p)-1) + 2)
end

function portview(def::Circuit, id, p)
    ports = def.interface.ports
    if isnothing(ports)
        throw(ArgumentError(lazy"The subcircuit $(id) exposes no ports, so $(p) cannot be used as a port in a pair connection."))
    end
    for port in ports
        if isequal(port.first, p)
            pk, nk = port.second
        return (scalarterminal(def, id, pk),
                nk === Ground ? Ground : scalarterminal(def, id, nk))
        end
    end
    throw(ArgumentError(lazy"The subcircuit $(id) has no port $(p). Its ports are $(collect(w.first for w in ports))."))
end

# key collision check for bare 2-tuples in pair context
function pincollision(def, k)
    return false
end
function pincollision(def::Circuit, k)
    haspin = any(pin -> isequal(pin.first, k), def.interface.pins)
    hasport = !isnothing(def.interface.ports) &&
        any(w -> isequal(w.first, k), def.interface.ports)
    return haspin && hasport
end

# === connection normalization ===

# The normalized form of one circuit level:
struct ParsedLevel
    table::ComponentTable
    # (name or nothing, endpoints as (instanceindex, terminal), containsground)
    groups::Vector{Tuple{Union{Nothing,String},Vector{Tuple{Int,Int}},Bool}}
    # per interface pin, in interface order: (instanceindex, terminal)
    pinendpoints::Vector{Tuple{Int,Int}}
    pinkeys::Vector{Any}
    # terminals auto-tied to ground by grounded multiports: (instanceindex, terminal)
    groundties::Vector{Tuple{Int,Int}}
    # mutual inductors: (kindex, inductor1index, inductor2index)
    mutuals::Vector{NTuple{3,Int}}
end

# Resolve one scalar endpoint written in group context. Returns
# (instanceindex, terminal) or Ground.
function resolvescalar(table::ComponentTable, ep, context::AbstractString)
    if ep === Ground || ep isa GroundType
        return Ground
    elseif ep isa PinRef
        i = instanceindex(table, ep.instance, context)
        t = scalarterminal(table.defs[i], ep.instance, ep.key)
        # a declared ground instance is the reference net, not a device:
        # its terminal is the ground sentinel, so every downstream rule
        # (net naming, interface pin restrictions) applies uniformly
        return table.defs[i] isa GroundType ? Ground : (i, t)
    elseif ep isa PortRef
        throw(ArgumentError(lazy"A PortRef is a two terminal port view and cannot be a member of a scalar connection group ($(context)). Use the port in a pair connection or address its terminals individually."))
    elseif ep isa Tuple && length(ep) == 2
        i = instanceindex(table, ep[1], context)
        t = scalarterminal(table.defs[i], ep[1], ep[2])
        return table.defs[i] isa GroundType ? Ground : (i, t)
    elseif ep isa Tuple && length(ep) == 3
        i = instanceindex(table, ep[1], context)
        return (i, portterminal(table.defs[i], ep[1], ep[2], ep[3]))
    else
        throw(ArgumentError(lazy"Unrecognized endpoint $(ep) in $(context). Endpoints are (instance, terminal), (instance, port, terminal), Ground, PinRef, or PortRef."))
    end
end

# Resolve one side of a pair connection. Returns either
# (:scalar, endpoint) or (:port, positive, negative) where endpoints are
# (instanceindex, terminal) or Ground.
function resolvepairside(table::ComponentTable, ep, context::AbstractString)
    if ep isa PortRef
        i = instanceindex(table, ep.instance, context)
        pos, neg = portview(table.defs[i], ep.instance, ep.key)
        return (:port, attach(i, pos), attach(i, neg))
    elseif ep isa Tuple && length(ep) == 2
        i = instanceindex(table, ep[1], context)
        def = table.defs[i]
        if hasports(def)
            if pincollision(def, ep[2])
                throw(ArgumentError(lazy"The key $(ep[2]) of the subcircuit $(ep[1]) exists both as a pin and as a port, which is ambiguous in a pair connection. Use PortRef($(repr(ep[1])), $(repr(ep[2]))) or PinRef($(repr(ep[1])), $(repr(ep[2])))."))
            end
            if def isa Circuit && !any(w -> isequal(w.first, ep[2]), def.interface.ports)
                # the key is only a pin: scalar fallback
                return (:scalar, resolvescalar(table, ep, context))
            end
            pos, neg = portview(def, ep[1], ep[2])
            return (:port, attach(i, pos), attach(i, neg))
        else
            return (:scalar, resolvescalar(table, ep, context))
        end
    else
        return (:scalar, resolvescalar(table, ep, context))
    end
end

attach(i::Int, t::Int) = (i, t)
attach(i::Int, g::GroundType) = Ground

# Normalize the user connection collection into scalar groups.
function parseconnections(table::ComponentTable, connections)
    groups = Vector{Tuple{Union{Nothing,String},Vector{Tuple{Int,Int}},Bool}}()
    for (ci, entry) in enumerate(connections)
        context = lazy"connection $(ci)"
        if entry isa Net
            name = string(entry.name)
            pushgroup!(groups, name, table, entry.endpoints, context)
        elseif entry isa Pair
            a = resolvepairside(table, entry.first, context)
            b = resolvepairside(table, entry.second, context)
            if a[1] == :scalar && b[1] == :scalar
                addgroup!(groups, nothing, Any[a[2], b[2]])
            elseif a[1] == :port && b[1] == :port
                addgroup!(groups, nothing, Any[a[2], b[2]])
                addgroup!(groups, nothing, Any[a[3], b[3]])
            else
                scalarside = a[1] == :scalar ? entry.first : entry.second
                portside = a[1] == :port ? entry.first : entry.second
                throw(ArgumentError(lazy"The pair connection $(ci) bonds the two terminal port $(portside) to the scalar endpoint $(scalarside), which have different arities. Bond two ports, or connect scalar terminals individually."))
            end
        elseif entry isa Union{Tuple,AbstractVector}
            pushgroup!(groups, nothing, table, entry, context)
        else
            throw(ArgumentError(lazy"Unrecognized connection entry $(entry). Connections are endpoint groups (tuples), pairs (port bonds), or Net entries."))
        end
    end
    return groups
end

function pushgroup!(groups, name, table::ComponentTable, endpoints, context)
    if length(endpoints) < 1
        throw(ArgumentError(lazy"The connection group $(context) is empty."))
    end
    # a common mistake is a single endpoint where a group of endpoints is
    # expected, e.g. ((:l1, 1)) — which is just (:l1, 1) — instead of
    # ((:l1, 1),); diagnose it instead of complaining about the elements
    if endpointlike(table, endpoints)
        ep = Tuple(endpoints)
        throw(ArgumentError(lazy"The connection group $(context) is the single endpoint $(ep) rather than a collection of endpoints. Wrap it as ($(ep),) or [$(ep)] to make a one-endpoint group."))
    end
    resolved = Any[resolvescalar(table, ep, context) for ep in endpoints]
    addgroup!(groups, name, resolved)
    return nothing
end

# whether a would-be group of endpoints is itself shaped like one endpoint:
# (instance, terminal) or (instance, port, terminal) with a known instance
function endpointlike(table::ComponentTable, endpoints)
    if (endpoints isa Tuple || endpoints isa AbstractVector) &&
            2 <= length(endpoints) <= 3
        id = first(endpoints)
        if (id isa Symbol || id isa AbstractString || id isa Integer) &&
                haskey(table.index, id) &&
                all(k isa Integer for k in Iterators.drop(endpoints, 1))
            return true
        end
    end
    return false
end

function addgroup!(groups, name, resolved::Vector)
    endpoints = Tuple{Int,Int}[]
    hasground = false
    for r in resolved
        if r === Ground
            hasground = true
        else
            push!(endpoints, r)
        end
    end
    push!(groups, (name, endpoints, hasground))
    return nothing
end

# === interface parsing ===

function parseinterface(table::ComponentTable, interface)
    if isnothing(interface)
        return Tuple{Int,Int}[], Any[]
    end
    if !(interface isa Interface)
        throw(ArgumentError(lazy"The interface must be an Interface or nothing; got $(typeof(interface))."))
    end
    pinendpoints = Tuple{Int,Int}[]
    pinkeys = Any[]
    for pin in interface.pins
        if !(pin isa Pair)
            throw(ArgumentError(lazy"Each interface pin must be a Pair of a key and an endpoint; got $(typeof(pin))."))
        end
        k = pin.first
        if any(isequal(k), pinkeys)
            throw(ArgumentError(lazy"The interface pin key $(k) is not unique."))
        end
        r = resolvescalar(table, pin.second, lazy"interface pin $(k)")
        if r === Ground
            throw(ArgumentError(lazy"The interface pin $(k) maps to Ground. Pins must map to component terminals; use Ground directly in the parent or as the negative side of an interface port."))
        end
        push!(pinkeys, k)
        push!(pinendpoints, r)
    end
    if isempty(pinendpoints)
        throw(ArgumentError("An Interface must expose at least one pin."))
    end
    if !isnothing(interface.ports)
        portkeys = Any[]
        for port in interface.ports
            if !(port isa Pair)
                throw(ArgumentError(lazy"Each interface port must be a Pair of a key and a (positive, negative) tuple of pin keys; got $(typeof(port))."))
            end
            k = port.first
            if any(isequal(k), portkeys)
                throw(ArgumentError(lazy"The interface port key $(k) is not unique."))
            end
            push!(portkeys, k)
            v = port.second
            if !(v isa Tuple && length(v) == 2)
                throw(ArgumentError(lazy"The interface port $(k) must map to a (positive, negative) tuple of pin keys; got $(v)."))
            end
            pk, nk = v
            if !any(isequal(pk), pinkeys)
                throw(ArgumentError(lazy"The positive side of interface port $(k) is $(pk), which is not a pin key."))
            end
            if !(nk === Ground) && !any(isequal(nk), pinkeys)
                throw(ArgumentError(lazy"The negative side of interface port $(k) is $(nk), which is neither a pin key nor Ground."))
            end
        end
    end
    return pinendpoints, pinkeys
end

# === mutual inductors and grounded ties ===

function parsemutuals(table::ComponentTable)
    mutuals = NTuple{3,Int}[]
    for (i, def) in enumerate(table.defs)
        if def isa MutualInductor
            i1 = get(table.index, def.inductor1, 0)
            i2 = get(table.index, def.inductor2, 0)
            if i1 == 0 || i2 == 0
                missingid = i1 == 0 ? def.inductor1 : def.inductor2
                throw(ArgumentError(lazy"The mutual inductor $(table.ids[i]) couples $(def.inductor1) and $(def.inductor2), but $(missingid) does not exist in this circuit."))
            end
            for j in (i1, i2)
                if !(table.defs[j] isa Inductor)
                    throw(ArgumentError(lazy"The mutual inductor $(table.ids[i]) couples $(table.ids[j]), which is a $(typeof(table.defs[j])); mutual inductors couple Inductor instances."))
                end
            end
            push!(mutuals, (i, i1, i2))
        end
    end
    return mutuals
end

function parsegroundties(table::ComponentTable)
    ties = Tuple{Int,Int}[]
    for (i, def) in enumerate(table.defs)
        if (def isa ScatteringParameters || def isa GaussianChannel) && isgrounded(def)
            for p in 1:componentnports(def)
                push!(ties, (i, 2*(p-1) + 2))
            end
        end
    end
    return ties
end

# === one level parse (also the constructor validation) ===

function parsecircuitlevel(components, connections, interface)
    table = componenttable(components)
    groups = parseconnections(table, connections)
    pinendpoints, pinkeys = parseinterface(table, interface)
    groundties = parsegroundties(table)
    mutuals = parsemutuals(table)
    return ParsedLevel(table, groups, pinendpoints, pinkeys, groundties,
        mutuals)
end

parsecircuitlevel(c::Circuit) = parsecircuitlevel(c.components, c.connections,
    c.interface)


# === node names and node ordering ===
#
# Both are used by `compile` below, which interns every terminal it
# visits and then renumbers the result.

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
function calcnodesorting(uniquenodevector::Vector{String};
    sorting::Symbol = :number)

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
                throw(ArgumentError(lazy"Failed to parse the nodes as integers. Try setting the keyword argument `sorting=:name` or `sorting=:none`."))
            end
        end
        sortperm!(uniquenodevectorsortindices, uniquenodevectorints, initialized=true)

    elseif sorting == :none
        # don't perform any sorting. keep the nodes in
        # order of first appearance, except move the
        # ground node first later as always
        nothing
    else
        throw(ArgumentError(lazy"Unknown sorting type."))
    end

    # find the ground node. error if we don't find it.
    groundnodeindex = findgroundnodeindex(uniquenodevector)

    if groundnodeindex == 0
        throw(ArgumentError(lazy"No ground node found in netlist."))
    end

    # if the ground index is not the first after sorting, move it to the front
    # and shift the nodes which sorted before it back by one
    if uniquenodevectorsortindices[1] != groundnodeindex
        groundpos = findfirst(==(groundnodeindex), uniquenodevectorsortindices)
        for j = groundpos:-1:2
            uniquenodevectorsortindices[j] = uniquenodevectorsortindices[j-1]
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
        nodeindexvector::Vector{Int};sorting::Symbol = :name)
    order = calcnodesorting(uniquenodevector; sorting = sorting)
    return sortnodes(uniquenodevector, nodeindexvector, order)[1:2]
end

"""
    noderenumbering(uniquenodevector, order)

Where each node ends up once the nodes are sorted: `renumber[j]` is the new
index of the node which was `j`.

Anything holding a node index from before the sort is renumbered with this,
which is how something that is not a component -- the terminals of a
scattering block -- can be carried through the sort without being given a
row in the component table to ride on.
"""
noderenumbering(order::Vector{Int}) = sortperm(order)

"""
    sortnodes(uniquenodevector, nodeindexvector, order)

Sort the nodes by a precomputed ordering, returning the sorted names, the
component node indices, and the renumbering the sort applied.
"""
function sortnodes(uniquenodevector::Vector{String},
        nodeindexvector::Vector{Int}, order::Vector{Int})

    nodeindices = zeros(eltype(nodeindexvector),2,length(nodeindexvector)÷2)

    nodevectorsortindices = noderenumbering(order)

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

    return uniquenodevector[order], nodeindices, nodevectorsortindices
end
