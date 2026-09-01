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


# === flattening the hierarchy ===


# === a small union-find over wire indices ===

mutable struct WireForest
    parent::Vector{Int}
    size::Vector{Int}
end
WireForest() = WireForest(Int[], Int[])

function newwire!(f::WireForest)
    push!(f.parent, length(f.parent) + 1)
    push!(f.size, 1)
    return length(f.parent)
end

function findwire(f::WireForest, i::Int)
    while f.parent[i] != i
        f.parent[i] = f.parent[f.parent[i]] # path halving
        i = f.parent[i]
    end
    return i
end

function unionwires!(f::WireForest, a::Int, b::Int)
    ra = findwire(f, a)
    rb = findwire(f, b)
    if ra == rb
        return ra
    end
    if f.size[ra] < f.size[rb]
        ra, rb = rb, ra
    end
    f.parent[rb] = ra
    f.size[ra] += f.size[rb]
    return ra
end

# === elaborated circuit ===

"""
    ElaboratedCircuit(definitions, definitionof, instancepaths,
        terminaloffsets, terminalnets, netnames, couplings)

The flattened result of [`elaborate`](@ref): hierarchy resolved, definitions
deduplicated by identity, and nets assigned dense integer indices. This
structure is immutable after construction and safe to share across threads.

# Fields
- `definitions::Vector{Any}`: the unique component definitions, deduplicated
    by object identity, so that a thousand instances of one shared
    definition store its data once.
- `definitionof::Vector{Int}`: for each flattened primitive instance, the
    index of its definition in `definitions`.
- `instancepaths::Vector{String}`: the hierarchical path of each instance,
    such as "cell37/cap", with "/" as the separator.
- `terminaloffsets::Vector{Int}`: offsets into `terminalnets` in CSR layout;
    the terminals of instance `i` are
    `terminalnets[terminaloffsets[i]:terminaloffsets[i+1]-1]`.
- `terminalnets::Vector{Int}`: the net index of every instance terminal.
    Net 1 is the ground net.
- `netnames::Vector{String}`: the net names; `netnames[1] == "0"` is the
    ground net. User supplied `Net` names win over automatic names; nested
    names are hierarchical, such as "cell37/net2".
- `couplings::Vector{NTuple{3,Int}}`: for each mutual inductor, the
    flattened instance indices `(mutualinductor, inductor1, inductor2)`.
"""
struct ElaboratedCircuit
    definitions::Vector{Any}
    definitionof::Vector{Int}
    instancepaths::Vector{String}
    terminaloffsets::Vector{Int}
    terminalnets::Vector{Int}
    netnames::Vector{String}
    couplings::Vector{NTuple{3,Int}}
end

"""
    ninstances(elab::ElaboratedCircuit)

The number of flattened primitive instances.
"""
ninstances(elab::ElaboratedCircuit) = length(elab.definitionof)

"""
    nnets(elab::ElaboratedCircuit)

The number of nets including the ground net.
"""
nnets(elab::ElaboratedCircuit) = length(elab.netnames)

"""
    instanceterminals(elab::ElaboratedCircuit, i::Integer)

A view of the net indices of the terminals of flattened instance `i`.
"""
function instanceterminals(elab::ElaboratedCircuit, i::Integer)
    return view(elab.terminalnets,
        elab.terminaloffsets[i]:(elab.terminaloffsets[i+1]-1))
end

"""
    instancedefinition(elab::ElaboratedCircuit, i::Integer)

The component definition of flattened instance `i`.
"""
instancedefinition(elab::ElaboratedCircuit, i::Integer) =
    elab.definitions[elab.definitionof[i]]

# === flattening state ===

mutable struct FlattenState
    wires::WireForest
    groundwire::Int
    definitions::Vector{Any}
    defindex::IdDict{Any,Int}
    definitionof::Vector{Int}
    instancepaths::Vector{String}
    terminaloffsets::Vector{Int}
    terminalwires::Vector{Int}
    couplings::Vector{NTuple{3,Int}}
    # net name candidates: wire => (depth, isauto, sequence, name or path)
    usernames::Vector{Tuple{Int,Int,String,Int}}   # (depth, seq, qualifiedname, wire)
    # The automatic name of a net comes from the shallowest, earliest
    # terminal on it. That was one (depth, seq, path, wire) tuple per
    # terminal, walked again afterwards; it is now the depth and the level
    # path of each instance, which its terminals share, and the pass which
    # already walks every terminal to find its net does the selection as it
    # goes. The terminal's own index orders the candidates, because the
    # sequence counter only ever increased along that same walk.
    levelpaths::Vector{String}
    autodepth::Vector{Int}                         # per instance
    autopathid::Vector{Int}                        # per instance
    parsedcache::IdDict{Any,ParsedLevel}
    active::IdDict{Any,Nothing}
    maxdepth::Int
    seq::Int
end

function FlattenState(maxdepth::Int)
    wires = WireForest()
    ground = newwire!(wires)
    return FlattenState(wires, ground, Any[], IdDict{Any,Int}(), Int[],
        String[], Int[1], Int[], NTuple{3,Int}[],
        Tuple{Int,Int,String,Int}[], String[], Int[], Int[],
        IdDict{Any,ParsedLevel}(), IdDict{Any,Nothing}(), maxdepth, 0)
end

nextseq!(st::FlattenState) = (st.seq += 1; st.seq)

function definitionindex!(st::FlattenState, def)
    i = get(st.defindex, def, 0)
    if i == 0
        push!(st.definitions, def)
        i = length(st.definitions)
        st.defindex[def] = i
    end
    return i
end

joinpath_(path::String, id) = isempty(path) ? string(id) : path * "/" * string(id)

"""
    elaborate(circuit::Circuit; maxdepth = 64)

Recursively flatten the hierarchy of `circuit` into an
[`ElaboratedCircuit`](@ref). Elaboration:

1. assigns each primitive instance a stable hierarchical path such as
   "cell37/cap";
2. substitutes subcircuit interface pins with parent nets and allocates
   fresh internal nets for every subcircuit instance;
3. deduplicates component definitions by object identity, so shared
   definitions and their data appear once;
4. parses and validates each unique circuit definition once, however many
   times it is instantiated;
5. resolves mutual inductor couplings to flattened instance indices;
6. rejects recursive circuit definitions.

The structural analysis of a repeated subcircuit definition is performed
once and reused for every instance; per instance work is proportional to
the instance's own size.
"""
function elaborate(circuit::Circuit; maxdepth::Integer = 64)
    st = FlattenState(Int(maxdepth))
    flattencircuit!(st, circuit, "", 0)
    return finishelaboration(st)
end

function flattencircuit!(st::FlattenState, c::Circuit, path::String,
        depth::Int)
    if haskey(st.active, c)
        location = isempty(path) ? "the top level" : path
        throw(ArgumentError(lazy"The circuit definition at $(location) contains itself, directly or indirectly. Recursive circuit definitions are not allowed."))
    end
    if depth > st.maxdepth
        throw(ArgumentError(lazy"The circuit hierarchy at $(path) exceeds the maximum depth $(st.maxdepth). Pass a larger maxdepth to elaborate if this is intentional."))
    end
    st.active[c] = nothing

    pd = get!(() -> parsecircuitlevel(c), st.parsedcache, c)
    table = pd.table
    push!(st.levelpaths, path)
    pathid = length(st.levelpaths)

    # wires of each local instance's terminals (for subcircuits: pin wires)
    instwires = Vector{Vector{Int}}(undef, length(table.ids))
    # flattened global index of each local primitive instance (0 for subcircuits)
    localglobal = zeros(Int, length(table.ids))

    for (i, def) in enumerate(table.defs)
        id = table.ids[i]
        if def isa Circuit
            instwires[i] = flattencircuit!(st, def, joinpath_(path, id),
                depth + 1)
        elseif def isa GroundType
            # a declared ground instance is the reference net, not a device:
            # it flattens to no instance, and every reference to its
            # terminal already resolved to the ground sentinel at parse
            # time, so its wire list is never consulted
            instwires[i] = Int[]
        else
            n = nterminals(def)
            w = Vector{Int}(undef, n)
            for t in 1:n
                w[t] = newwire!(st.wires)
                push!(st.terminalwires, w[t])
            end
            instwires[i] = w
            push!(st.autodepth, depth)
            push!(st.autopathid, pathid)
            push!(st.definitionof, definitionindex!(st, def))
            push!(st.instancepaths, joinpath_(path, id))
            push!(st.terminaloffsets, st.terminaloffsets[end] + n)
            localglobal[i] = length(st.definitionof)
        end
    end

    # grounded multiport reference ties
    for (i, t) in pd.groundties
        unionwires!(st.wires, instwires[i][t], st.groundwire)
    end

    # connection groups
    for (name, endpoints, hasground) in pd.groups
        first = hasground ? st.groundwire :
            instwires[endpoints[1][1]][endpoints[1][2]]
        for (i, t) in endpoints
            unionwires!(st.wires, first, instwires[i][t])
        end
        if !isnothing(name)
            push!(st.usernames, (depth, nextseq!(st), joinpath_(path, name),
                hasground ? st.groundwire : first))
        end
    end

    # mutual inductor couplings resolved to flattened indices
    for (k, i1, i2) in pd.mutuals
        push!(st.couplings, (localglobal[k], localglobal[i1],
            localglobal[i2]))
    end

    # interface pin wires for the parent
    pinwires = Vector{Int}(undef, length(pd.pinendpoints))
    for (j, (i, t)) in enumerate(pd.pinendpoints)
        pinwires[j] = instwires[i][t]
    end

    delete!(st.active, c)
    return pinwires
end

function finishelaboration(st::FlattenState)
    # compress the union-find into dense net indices; net 1 is ground
    netofroot = Dict{Int,Int}()
    netofroot[findwire(st.wires, st.groundwire)] = 1
    terminalnets = Vector{Int}(undef, length(st.terminalwires))
    # the shallowest, earliest terminal on each net, which names it when no
    # user name does. Found on this walk rather than on a second one over a
    # parallel record of every terminal.
    autopath = Dict{Int,Tuple{Int,Int,String}}()
    nnets = 1
    inst = 1
    ninst = length(st.autodepth)
    for (i, w) in enumerate(st.terminalwires)
        r = findwire(st.wires, w)
        n = get(netofroot, r, 0)
        if n == 0
            nnets += 1
            n = nnets
            netofroot[r] = n
        end
        terminalnets[i] = n
        # which instance this terminal belongs to, walked alongside
        while inst < ninst && i >= st.terminaloffsets[inst+1]
            inst += 1
        end
        if n != 1 && inst <= ninst
            d = st.autodepth[inst]
            best = get(autopath, n, (typemax(Int), typemax(Int), ""))
            if (d, i) < (best[1], best[2])
                autopath[n] = (d, i, st.levelpaths[st.autopathid[inst]])
            end
        end
    end

    # assign names: user names win, shallowest then earliest; otherwise an
    # automatic hierarchical name from the shallowest, earliest level which
    # touches the net.
    username = Dict{Int,Tuple{Int,Int,String}}()
    for (depth, seq, name, w) in st.usernames
        n = get(netofroot, findwire(st.wires, w), 0)
        (n == 0 || n == 1) && continue # dropped or ground (always "0")
        best = get(username, n, (typemax(Int), typemax(Int), ""))
        if (depth, seq) < (best[1], best[2])
            username[n] = (depth, seq, name)
        end
    end
    netnames = Vector{String}(undef, nnets)
    netnames[1] = "0"
    autocounter = Dict{String,Int}()
    for n in 2:nnets
        if haskey(username, n)
            netnames[n] = username[n][3]
        else
            path = autopath[n][3]
            k = get(autocounter, path, 0) + 1
            autocounter[path] = k
            netnames[n] = joinpath_(path, "net$(k)")
        end
    end

    # net names must be unique
    seen = Dict{String,Int}()
    for (n, name) in enumerate(netnames)
        prev = get(seen, name, 0)
        if prev != 0
            throw(ArgumentError(lazy"The net name \"$(name)\" is used for two distinct nets. Rename one of the Net entries."))
        end
        seen[name] = n
    end

    return ElaboratedCircuit(st.definitions, st.definitionof,
        st.instancepaths, st.terminaloffsets, terminalnets, netnames,
        st.couplings)
end


# === lowering to the compiled tables ===

# Compiling an elaborated typed circuit.
#
# `ElaboratedCircuit` is the last heterogeneous representation: its
# definitions are arbitrary component objects and reading it means asking
# what each one is. Everything below this file works from grouped, integer
# indexed tables instead, so the question is asked once, here.
#
# The compiled form is a flat component table plus groups indexing into it.
# The flat table is what the matrix builders, the netlist export and the
# sensitivities read, component by component. The groups are what the
# assembly plans, the ports and the scattering stamps read, and they are
# what makes those O(components in the group) rather than a scan.
#
# A port is compiled with its own reference impedance and a direct index to
# the environment it owns. Nothing here or below looks for a resistor which
# happens to share a port's branch.

"""
    CompiledPort

An analysis port with its reference impedance and the environment it owns.

`environment` is the flat component index of the port's own termination, or
`0` when the port owns none. It is recorded when the port is compiled rather
than discovered later by looking for a resistor on the port's branch, so a
port may share its terminals with any number of ordinary device resistors.

`zref` is the reference impedance as written, which may still be symbolic;
binding resolves it to a finite positive real number.
"""
struct CompiledPort
    number::Int
    positivenode::Int
    negativenode::Int
    zref::Any
    environment::Int
    component::Int
end

"""
    CompiledScatteringBlock

A multiport scattering block with its whole terminal map intact.

The block keeps its identity through compilation: `signalnodes[p]` and
`refnodes[p]` are the nodes of port `p`, and `components[p]` is the flat
index of the entry the parsed view exposes for that port. The parsed
representation can only hold two terminal components, so it sees one entry
per port and has to reassemble the block afterwards by object identity;
nothing reading the compiled form needs to.
"""
struct CompiledScatteringBlock
    definition::Any
    signalnodes::Vector{Int}
    refnodes::Vector{Int}
    components::Vector{Int}
    path::String
end

"""
    CompiledCircuit

An elaborated circuit lowered to a flat component table with typed groups
indexing into it.

# Fields

The flat table, in the order the components were elaborated:

- `componentnames`: hierarchical instance paths.
- `componenttypes`: the lowered type symbol of each component.
- `componentvalues`: the value of each component, as written.
- `nodeindices`: the two node indices of each component, ground being `1`.
- `componenttemperatures`: the temperature of each component which states one.
- `mutualinductorbranchnames`: the coupled inductor names, two per `:K`.
- `nodenames`, `Nnodes`: the sorted node names and their count.
- `componentnamedict`: name to flat index.

The groups, each a vector of flat indices:

- `capacitors`, `resistors`, `inductors`, `junctions`, `nonlinearinductors`,
  `currentsources`, `mutualinductors`.

and the two which keep their own structure:

- `ports::Vector{CompiledPort}`, ordered as elaborated.
- `scatteringblocks::Vector{CompiledScatteringBlock}`, one per block instance.

See [`compile`](@ref).
"""
struct CompiledCircuit
    nodenames::Vector{String}
    nodeindices::Matrix{Int}
    Nnodes::Int
    componentnames::Vector{String}
    componenttypes::Vector{Symbol}
    componentvalues::Vector
    componentnamedict::Dict{String,Int}
    componenttemperatures::Dict{Int,Float64}
    mutualinductorbranchnames::Vector{String}
    capacitors::Vector{Int}
    resistors::Vector{Int}
    inductors::Vector{Int}
    junctions::Vector{Int}
    nonlinearinductors::Vector{Int}
    currentsources::Vector{Int}
    mutualinductors::Vector{Int}
    ports::Vector{CompiledPort}
    scatteringblocks::Vector{CompiledScatteringBlock}
end

function Base.show(io::IO, c::CompiledCircuit)
    print(io, "CompiledCircuit(", length(c.componenttypes), " components, ",
        c.Nnodes, " nodes, ", length(c.ports), " ports")
    isempty(c.scatteringblocks) ||
        print(io, ", ", length(c.scatteringblocks), " scattering blocks")
    print(io, ")")
end

"""
    ncomponents(c::CompiledCircuit)

The number of entries in the flat component table.
"""
ncomponents(c::CompiledCircuit) = length(c.componenttypes)

# A legacy netlist names its nodes with integers and has always been sorted
# by their numeric value, which the entry points below preserve. Anything
# else is a typed circuit, whose hierarchical net names are not integers.
defaultsorting(circuit) = circuit isa AbstractVector ? :number : :name

# === lowering an elaborated circuit to component entries ===
#
# The `isa` chain in `compile` handles the common components inline;
# everything rarer, and every component the solver does not support,
# falls through to `lowercomponent` so the diagnostics stay in one place.
# The legacy netlist's untyped nonlinear element adds its own method in
# `legacyadapter.jl`.

function lowercomponent(def::Inductor, path) 
    return :L, def.L
end
lowercomponent(def::Capacitor, path) = :C, def.C
lowercomponent(def::Resistor, path) = :R, def.R
lowercomponent(def::CurrentSource, path) = :I, def.I
lowercomponent(def::Port, path) = :P, def.number
lowercomponent(def::MutualInductor, path) = :K, def.K
function lowercomponent(def::NonlinearInductor, path)
    if issinusoidal(def)
        return :Lj, def.L0
    end
    throw(ComponentNotSupportedError(lazy"the NonlinearInductor at $(path) has a non-sinusoidal current-phase relation, which the solver does not yet support. It parsed, validated, and elaborated successfully. Currently solvable nonlinear elements have the sinusoidal relation of JosephsonJunction."))
end
function lowercomponent(def::VoltageSource, path)
    throw(ComponentNotSupportedError(lazy"the VoltageSource at $(path) is not supported by the solver, which matches the legacy parser (voltage sources are not currently supported)."))
end
function lowercomponent(def::GaussianChannel, path)
    throw(ComponentNotSupportedError(lazy"the GaussianChannel at $(path) is not yet supported by the harmonic balance solvers. It parsed, validated, and elaborated successfully; solver support for Gaussian channels is planned. Currently solvable components: Inductor, Capacitor, Resistor, JosephsonJunction, MutualInductor, CurrentSource, and Port."))
end
function lowercomponent(def, path)
    throw(ComponentNotSupportedError(lazy"the component $(typeof(def)) at $(path) is not supported by the solver."))
end

# narrow a Vector{Any} to a concrete element type when possible, which is
# what makes a fully numeric circuit's values a Vector{Float64} rather than
# a vector of boxes.
function tightenvalues(values::Vector{Any})
    return map(identity, values)
end

# the temperature a component states, or `nothing`. Only the components which
# can dissipate carry one; a scattering block states its temperature through
# its noise model instead (see [`ThermalEquilibrium`](@ref)).
componenttemperature(def::Resistor) = def.temperature
componenttemperature(def::Capacitor) = def.temperature
componenttemperature(def::Inductor) = def.temperature
componenttemperature(def) = nothing

"""
    compile(elab::ElaboratedCircuit; sorting = :name)
    compile(circuit::Circuit; sorting = :name)

Lower an elaborated typed circuit to a [`CompiledCircuit`](@ref).

Component and node ordering are exactly those of the parsed representation
this replaces, so no downstream result moves.

Only components the solvers support can be lowered; a
[`GaussianChannel`](@ref), a non-sinusoidal [`NonlinearInductor`](@ref), or a
[`ScatteringParameters`](@ref) with an arbitrary noise covariance raises a
[`ComponentNotSupportedError`](@ref) naming the instance. The default
`sorting` is `:name` because automatic hierarchical net names are not
integers; see [`sortnodes`](@ref).
"""
function compile(elab::ElaboratedCircuit; sorting::Symbol = :name)

    N = ninstances(elab)
    componentnames = String[]
    componenttypes = Symbol[]
    componentvalues = Any[]
    nodeindexvector = Int[]
    componenttemperatures = Dict{Int,Float64}()
    mutualinductorbranchnames = String[]
    sizehint!(componentnames, N)
    sizehint!(componenttypes, N)
    sizehint!(componentvalues, N)
    sizehint!(nodeindexvector, 2*N)

    uniquenodedict = Dict{String,Int}()
    uniquenodevector = String[]

    # the two coupled inductor names of each mutual inductor, by instance
    couplingnames = Dict{Int,Tuple{String,String}}()
    for (k, i1, i2) in elab.couplings
        couplingnames[k] = (elab.instancepaths[i1], elab.instancepaths[i2])
    end

    ports = CompiledPort[]
    scatteringblocks = CompiledScatteringBlock[]
    legacyenvironments = Pair{Int,String}[]

    for i in 1:N
        def = instancedefinition(elab, i)
        path = elab.instancepaths[i]

        # A scattering block is one component with two terminals per port.
        # The parsed table holds two terminal components only, so the block
        # is spread over one entry per port there; the compiled block keeps
        # the whole terminal map and the entries it owns.
        if def isa ScatteringParameters
            if def.noise isa NoiseCovariance
                throw(ComponentNotSupportedError(lazy"the ScatteringParameters at $(path) has an arbitrary noise covariance, which the harmonic balance solvers do not yet support. Use Passive, Lossless, or ThermalEquilibrium."))
            end
            terminals = instanceterminals(elab, i)
            n = def.nports
            signalnodes = Vector{Int}(undef, n)
            refnodes = Vector{Int}(undef, n)
            entries = Vector{Int}(undef, n)
            for p in 1:n
                push!(componentnames, path * "/port" * string(p))
                push!(componenttypes, :S)
                push!(componentvalues, ScatteringStamp(def, p))
                if def.noise isa ThermalEquilibrium
                    componenttemperatures[length(componentnames)] =
                        Float64(def.noise.temperature)
                end
                s = processnode(uniquenodedict, uniquenodevector,
                    elab.netnames[terminals[2*p-1]])
                r = processnode(uniquenodedict, uniquenodevector,
                    elab.netnames[terminals[2*p]])
                push!(nodeindexvector, s)
                push!(nodeindexvector, r)
                signalnodes[p] = s
                refnodes[p] = r
                entries[p] = length(componentnames)
            end
            push!(scatteringblocks, CompiledScatteringBlock(def, signalnodes,
                refnodes, entries, path))
            continue
        end

        # An `isa` chain rather than dispatch on `lowercomponent`, ordered
        # by how many of each a large circuit holds. Selecting an operation
        # by dynamic dispatch is expensive in the abstract: over 8192
        # heterogeneous components whose type changes every element, it
        # costs 114 us at two types and 420 us at eight against 2 us for a
        # branch chain, the step between four and eight being the callsite
        # method cache giving out.
        #
        # It is not expensive here, and this chain is not why the compiler
        # is fast: on a circuit of 8402 instances spanning seven component
        # types, lowering by dispatch and lowering by this chain both take
        # 2.45 ms, because the work around them -- building the instance
        # path, interning the nodes -- dwarfs the selection. The chain is
        # kept because it costs nothing and does not depend on the component
        # set staying small, not because it was measured to help.
        #
        # Where the selection genuinely is in a loop, it is not made at all:
        # the assembly reads one homogeneous group at a time and never asks
        # what a component is. Everything rarer, and every unsupported
        # component, falls through to `lowercomponent`, which keeps the
        # diagnostics in one place.
        typesymbol, value = if def isa Capacitor
            (:C, def.C)
        elseif def isa NonlinearInductor
            issinusoidal(def) ? (:Lj, def.L0) : lowercomponent(def, path)
        elseif def isa Inductor
            (:L, def.L)
        elseif def isa Resistor
            (:R, def.R)
        elseif def isa Port
            (:P, def.number)
        elseif def isa CurrentSource
            (:I, def.I)
        elseif def isa MutualInductor
            (:K, def.K)
        else
            lowercomponent(def, path)
        end
        push!(componentnames, path)
        push!(componenttypes, typesymbol)
        push!(componentvalues, value)
        marker = length(componentnames)
        # a component which states its temperature keeps it; the rest take
        # the one the analysis is run at
        t = componenttemperature(def)
        isnothing(t) || (componenttemperatures[marker] = t)

        if typesymbol == :K
            l1, l2 = couplingnames[i]
            push!(mutualinductorbranchnames, l1)
            push!(mutualinductorbranchnames, l2)
            push!(nodeindexvector, 0)
            push!(nodeindexvector, 0)
            continue
        end

        terminals = instanceterminals(elab, i)
        if length(terminals) != 2
            throw(ComponentNotSupportedError(lazy"the component $(typeof(def)) at $(path) has $(length(terminals)) terminals; the solver supports two terminal components."))
        end
        n1 = processnode(uniquenodedict, uniquenodevector,
            elab.netnames[terminals[1]])
        n2 = processnode(uniquenodedict, uniquenodevector,
            elab.netnames[terminals[2]])
        push!(nodeindexvector, n1)
        push!(nodeindexvector, n2)

        if def isa Port
            # A matched port owns its external environment. It is emitted
            # here as an ordinary resistor across the port terminals, which
            # the conductance stamping, the solver scale and the parsed view
            # all consume unchanged. What makes it the port's own is the
            # index recorded on the port, not its position in the circuit,
            # so any further resistor on these terminals is an ordinary
            # device resistor.
            #
            # The generated name cannot collide: "<path>/termination" would
            # require the instance at "<path>" to be a subcircuit containing
            # an instance named "termination", and it is a Port.
            environment = 0
            if def.termination isa MatchedTermination
                push!(componentnames, path * "/termination")
                push!(componenttypes, :R)
                push!(componentvalues, def.Z0)
                push!(nodeindexvector, n1)
                push!(nodeindexvector, n2)
                environment = length(componentnames)
            elseif def.termination isa LegacyTermination
                # the environment already exists; it is named relative to the
                # port, so a legacy circuit instanced as a subcircuit still
                # resolves
                k = findlast('/', path)
                prefix = isnothing(k) ? "" : path[1:k]
                push!(legacyenvironments,
                    length(ports) + 1 => prefix*string(def.termination.component))
            end
            push!(ports, CompiledPort(def.number, n1, n2, def.Z0,
                environment, marker))
        end
    end

    if !haskey(uniquenodedict, "0")
        throw(ArgumentError("The circuit has no connection to Ground. Connect at least one endpoint to Ground; the ground net is required by the solver."))
    end

    nodenames, nodeindices = sortnodes(uniquenodevector, nodeindexvector;
        sorting = sorting)

    componentnamedict = Dict{String,Int}()
    sizehint!(componentnamedict, length(componentnames))
    for (i, name) in enumerate(componentnames)
        componentnamedict[name] = i
    end

    # a legacy port owns a resistor which already existed, so its index is
    # resolved once the whole table is known
    for (k, name) in legacyenvironments
        i = get(componentnamedict, name, 0)
        if iszero(i) || componenttypes[i] !== :R
            throw(ArgumentError(lazy"The port $(componentnames[ports[k].component]) names $(name) as its termination, which is not a resistor in this circuit."))
        end
        ports[k] = CompiledPort(ports[k].number, ports[k].positivenode,
            ports[k].negativenode, ports[k].zref, i, ports[k].component)
    end

    # sortnodes renumbers the nodes, so the port and block node indices
    # recorded above are in the pre-sort numbering and are refreshed from the
    # sorted table rather than translated.
    ports = [CompiledPort(p.number, nodeindices[1, p.component],
        nodeindices[2, p.component], p.zref, p.environment, p.component)
        for p in ports]

    # A port declared with `termination = nothing` is terminated by a
    # resistor the user placed across it, which is the one case where the
    # environment is not declared anywhere and has to be found. It is found
    # once, here, and recorded on the port like any other, so nothing below
    # this point searches for it. Two resistors across such a port leave the
    # reference impedance ambiguous, which is the error the geometric path
    # raised and still raises.
    for (k, p) in enumerate(ports)
        iszero(p.environment) || continue
        found = 0
        for i in eachindex(componenttypes)
            componenttypes[i] === :R || continue
            n1, n2 = nodeindices[1, i], nodeindices[2, i]
            ((n1 == p.positivenode && n2 == p.negativenode) ||
             (n1 == p.negativenode && n2 == p.positivenode)) || continue
            if !iszero(found)
                throw(ArgumentError(lazy"Only one resistor allowed per port."))
            end
            found = i
        end
        ports[k] = CompiledPort(p.number, p.positivenode, p.negativenode,
            p.zref, found, p.component)
    end
    for b in scatteringblocks
        for (p, k) in enumerate(b.components)
            b.signalnodes[p] = nodeindices[1, k]
            b.refnodes[p] = nodeindices[2, k]
        end
    end

    group(t) = [i for (i, s) in enumerate(componenttypes) if s === t]

    return CompiledCircuit(nodenames, nodeindices, length(uniquenodevector),
        componentnames, componenttypes, tightenvalues(componentvalues),
        componentnamedict, componenttemperatures, mutualinductorbranchnames,
        group(:C), group(:R), group(:L), group(:Lj), group(:NL), group(:I),
        group(:K), ports, scatteringblocks)
end

compile(circuit::Circuit; sorting::Symbol = :name) =
    compile(elaborate(circuit); sorting = sorting)

# a netlist of tuples becomes a typed circuit first
compile(netlist::AbstractVector; sorting::Symbol = :name) =
    compile(Circuit(netlist); sorting = sorting)

# and a circuit which is already compiled is itself, so the solver entry
# points accept one
compile(c::CompiledCircuit; sorting::Symbol = :name) = c

# === port and noise roles, read from the compiled circuit ===
#
# The parsed representation carries no role, so it recovers these by
# geometry: it looks for a resistor sharing a port's branch and calls that
# the port impedance, which forces the rule that a port may have exactly one
# resistor across it. A compiled port states which component is its own
# environment, so the same lists come out of the structure, an ordinary
# device resistor may share the port's terminals, and a port with no
# environment at all is representable.

"""
    portindicesnumbers(c::CompiledCircuit)

The flat indices and numbers of the ports, ordered by port number.
"""
function portindicesnumbers(c::CompiledCircuit)
    numbers = [p.number for p in c.ports]
    if !allunique(numbers)
        throw(ArgumentError(lazy"Duplicate ports are not allowed."))
    end
    branches = Tuple{Int,Int}[]
    for p in c.ports
        push!(branches, (p.positivenode, p.negativenode))
        push!(branches, (p.negativenode, p.positivenode))
    end
    if !allunique(branches)
        throw(ArgumentError(lazy"Only one port allowed per branch."))
    end
    sp = sortperm(numbers)
    return [p.component for p in c.ports][sp], numbers[sp]
end

"""
    portenvironmentindices(c::CompiledCircuit)

The flat index of each port's own environment, ordered by port number.

Every port must own one for the wave normalization the solvers currently
perform. A port declared with `termination = nothing` and no resistor of its
own has none, which is an error here for the same reason it was before.
"""
function portenvironmentindices(c::CompiledCircuit)
    sp = sortperm([p.number for p in c.ports])
    idx = [c.ports[i].environment for i in sp]
    for (k, i) in enumerate(idx)
        if iszero(i)
            throw(ArgumentError(lazy"The port $(c.componentnames[c.ports[sp[k]].component]) owns no termination and no resistor was found for it. Give the port a matched termination, or place a resistor across it."))
        end
    end
    return idx
end

"""
    noiseindices(c::CompiledCircuit, values)

The flat indices of the internal dissipative components: every resistor which
is not a port's own environment, and every capacitor or inductor with a
nonzero imaginary part.

A port environment is an external bath rather than an internal channel, so it
is excluded by role. Any other resistor across the port's terminals is an
ordinary device resistor and is included.
"""
function noiseindices(c::CompiledCircuit, values)
    owned = Set(p.environment for p in c.ports if !iszero(p.environment))
    out = Int[]
    for i in eachindex(c.componenttypes)
        t = c.componenttypes[i]
        if t === :R
            i in owned || push!(out, i)
        elseif (t === :C || t === :L) && values[i] isa Complex &&
                !iszero(values[i].im)
            push!(out, i)
        end
    end
    return out
end
