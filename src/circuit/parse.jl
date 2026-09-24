# The typed circuit a user writes, and the parse of one level of it.
#
# A circuit passes through three representations:
#
#   Circuit            what the user writes: components, connections, and
#                      the interface a subcircuit presents to its parent.
#                      Validated one level at a time by `parsecircuitlevel`
#                      (this file).
#   ElaboratedCircuit  the hierarchy flattened, by `elaborate`
#                      (circuit/compile.jl).
#   CompiledCircuit    flat integer indexed tables, by `compile`
#                      (circuit/compile.jl).
#
# The component models are in circuit/components.jl, how a component value
# becomes a number is in circuit/values.jl, and the deprecated tuple netlist
# is converted to a `Circuit` in circuit/legacy.jl.


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

# `Ground()` reads like a component constructor but returns the same
# singleton, so the two spellings cannot diverge.
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
        circuit = new{typeof(components),typeof(connections),
            typeof(interface)}(components, connections, interface)
        validate && parsecircuitlevel(circuit)
        return circuit
    end
end

# the length of a collection which knows it, for a capacity hint, and
# zero for an iterator which would have to be walked to count it
knownlength(x) = Base.IteratorSize(typeof(x)) isa Union{Base.HasLength,Base.HasShape} ?
    length(x) : 0

countelements(x) = Base.IteratorSize(x) isa Union{Base.HasLength,Base.HasShape} ?
    length(x) : count(Returns(true), x)

# === the netlist form ===

"""
    Circuit(netlist::AbstractVector; pins = nothing, ports = nothing)

Construct a [`Circuit`](@ref) from a netlist: a vector of entries
`(name, nodes..., component)`, one per component instance, each listing the
node of every terminal in terminal order, as a SPICE netlist does.

- `name` is the instance identifier, a `Symbol` or a string.
- The nodes are integers, strings or symbols. Node `0` (or `"0"`) is
  ground. Every entry naming a node joins its net, and the nets carry the
  node names.
- `component` is a typed component model, and the entry lists one node
  per terminal in terminal order: two for a lumped element or a port; for
  a [`ScatteringParameters`](@ref) block the signal and reference
  terminal of each port in turn, or the signal terminal alone when the
  block is `grounded`; for a subcircuit [`Circuit`](@ref) one per
  interface pin, in the order the pins were declared.
- A [`MutualInductor`](@ref) has no terminals. Its entry names the two
  inductors it couples in place of nodes, and the component is written
  `MutualInductor(K)`.

`pins` and `ports` declare an interface as they do for the connection-group
form, so a netlist can define a subcircuit.

The netlist is the connection-group form with the groups written out by
node, and builds exactly that: the components in netlist order and one
named connection group per node. What the connection-group form expresses
beyond a node list, such as the bundled port views of pair connections, is
written in that form.

A vector of entries carrying values under type-prefixed names,
`("C1", "1", "0", 1e-12)`, is the deprecated tuple netlist, which this
method reads, with a warning, through the adapter in circuit/legacy.jl.

# Examples
```jldoctest
julia> circuit = Circuit([
           (:p1, 1, 0, Port(1; Z0 = 50.0)),
           (:cc, 1, 2, Capacitor(100e-15)),
           (:jj, 2, 0, JosephsonJunction(1000e-12)),
           (:cj, 2, 0, Capacitor(1000e-15))]);

julia> compile(circuit).nodenames
3-element Vector{String}:
 "0"
 "1"
 "2"
```
"""
function Circuit(netlist::AbstractVector; pins = nothing, ports = nothing)
    typed = count(isnetlistentry, netlist)
    if typed == length(netlist) && typed > 0
        return netlistcircuit(netlist; pins = pins, ports = ports)
    elseif typed > 0
        throw(ArgumentError("The netlist mixes entries whose last element is a typed component with entries which are not; a netlist is one or the other."))
    end
    # entries ending in values rather than components are the deprecated
    # tuple netlist, which circuit/legacy.jl reads
    return legacycircuit(netlist, nothing; pins = pins, ports = ports)
end

# an entry of the netlist form ends in a component model; a tuple netlist
# entry ends in a value
isnetlistentry(e) = e isa Tuple && length(e) >= 2 &&
    (last(e) isa AbstractComponent || last(e) isa GroundType)

function netlistcircuit(netlist::AbstractVector; pins = nothing,
        ports = nothing)
    # one container type whatever the mixture and order of the component
    # models and the identifiers, so that a circuit is one type
    components = Vector{Pair{Any,Any}}(undef, length(netlist))
    nodeorder = String[]
    groups = Dict{String,Vector{Any}}()
    for (i, entry) in enumerate(netlist)
        name = first(entry)
        def = last(entry)
        nnodes = length(entry) - 2
        if !(name isa Union{Symbol,AbstractString})
            throw(ArgumentError(lazy"The netlist entry $(i) names its component $(name), which is not a Symbol or a string."))
        end
        if def isa MutualInductor
            # the entry names the inductors where the others name nodes
            if nnodes != 2
                throw(ArgumentError(lazy"The mutual inductor $(name) must name the two inductors it couples, as (name, inductor1, inductor2, MutualInductor(K)); its entry has $(nnodes) names."))
            end
            if isnothing(def.inductor1) && isnothing(def.inductor2)
                def = MutualInductor(def.K, entry[2], entry[3])
            elseif !(isequal(def.inductor1, entry[2]) &&
                    isequal(def.inductor2, entry[3]))
                throw(ArgumentError(lazy"The mutual inductor $(name) names $(entry[2]) and $(entry[3]) in its entry but $(def.inductor1) and $(def.inductor2) in its component; write MutualInductor(K) and let the entry name them."))
            end
            components[i] = Pair{Any,Any}(name, def)
            continue
        end
        nexpected = netlistterminalcount(def, name)
        if nnodes != nexpected
            throw(ArgumentError(lazy"The component $(name) has $(nexpected) terminals but its netlist entry lists $(nnodes) nodes; an entry is (name, nodes..., component) with one node per terminal."))
        end
        components[i] = Pair{Any,Any}(name, def)
        appendnetlistnodes!(groups, nodeorder, entry, def, name)
    end
    connections = Vector{Net{String,Vector{Any}}}(undef, length(nodeorder))
    for (i, label) in enumerate(nodeorder)
        group = groups[label]
        label == "0" && push!(group, Ground)
        connections[i] = Net(label, group)
    end
    return Circuit(components, connections; pins = pins, ports = ports)
end

# The endpoints of a netlist entry go straight into their node groups,
# one dispatch per component, without a vector of endpoints or a slice of
# the entry in between.
netlistterminalcount(def, name) = nterminals(def)
function netlistterminalcount(def::Union{ScatteringParameters,LinearizedScattering,GaussianChannel}, name)
    return isgrounded(def) ? componentnports(def) : 2*componentnports(def)
end
function netlistterminalcount(def::Circuit, name)
    if isnothing(def.interface)
        throw(ArgumentError(lazy"The subcircuit $(name) has no interface; give it pins to use it as a component."))
    end
    return length(def.interface.pins)
end

function appendnetlistnode!(groups, nodeorder, node, endpoint, name)
    label = string(node)
    if occursin('/', label)
        throw(ArgumentError(lazy"The node $(label) of $(name) contains the reserved hierarchical path separator \"/\"."))
    end
    group = get!(() -> (push!(nodeorder, label); Any[]), groups, label)
    push!(group, endpoint)
    return nothing
end

function appendnetlistnodes!(groups, nodeorder, entry, def, name)
    for t in 1:length(entry)-2
        appendnetlistnode!(groups, nodeorder, entry[t+1], (name, t), name)
    end
    return nothing
end
function appendnetlistnodes!(groups, nodeorder, entry,
        def::Union{ScatteringParameters,LinearizedScattering,GaussianChannel}, name)
    if isgrounded(def)
        for p in 1:componentnports(def)
            appendnetlistnode!(groups, nodeorder, entry[p+1], (name, p), name)
        end
    else
        for p in 1:componentnports(def), t in 1:2
            appendnetlistnode!(groups, nodeorder, entry[2*(p-1)+t+1],
                (name, p, t), name)
        end
    end
    return nothing
end
function appendnetlistnodes!(groups, nodeorder, entry, def::Circuit, name)
    for (t, pin) in enumerate(def.interface.pins)
        appendnetlistnode!(groups, nodeorder, entry[t+1], (name, pin.first), name)
    end
    return nothing
end

# A compact display: the stored collections can be large for generated
# circuits, and a nested subcircuit would otherwise print recursively.
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
port; a [`GaussianChannel`](@ref) has two per mode; a hierarchical
[`Circuit`](@ref) has one per interface pin; a
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
nterminals(c::LinearizedScattering) = 2*c.nports
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
hasports(c::LinearizedScattering) = true
hasports(c::GaussianChannel) = true
hasports(c::Circuit) = !isnothing(c.interface) && !isnothing(c.interface.ports)
hasports(c) = false

"""
    componentnports(c)

The number of ports of a multiport component: the ports of a
[`ScatteringParameters`](@ref), the modes of a [`GaussianChannel`](@ref).
Part of the connector protocol with [`nterminals`](@ref) and
[`hasports`](@ref).
"""
componentnports(c::ScatteringParameters) = c.nports
componentnports(c::LinearizedScattering) = c.nports
componentnports(c::GaussianChannel) = c.nmodes

"""
    isgrounded(c)

Whether a multiport component's second terminals are all tied to ground,
so that only its first terminals connect (`grounded = true` at
construction).
"""
isgrounded(c::ScatteringParameters) = c.grounded
isgrounded(c::LinearizedScattering) = c.grounded
isgrounded(c::GaussianChannel) = c.grounded

# === interface key lookup ===

# The keys of an interface, pins and ports, looked up in dictionaries
# whatever the interface's size; the declarations are kept in their order
# for the pin numbering and the diagnostics, and a port maps to its pair
# of terminals, a negative terminal of zero being Ground.
struct InterfaceIndex
    pins::Vector{Pair{Any,Any}}
    pinmap::Dict{Any,Int}
    portkeys::Vector{Any}
    portmap::Dict{Any,Tuple{Int,Int}}
end

function InterfaceIndex(interface)
    if !(interface isa Interface)
        throw(ArgumentError(lazy"The interface must be an Interface or nothing; got $(typeof(interface))."))
    end
    pins = Pair{Any,Any}[]
    pinmap = Dict{Any,Int}()
    n = knownlength(interface.pins)
    sizehint!(pins, n)
    sizehint!(pinmap, n)
    for pin in interface.pins
        if !(pin isa Pair)
            throw(ArgumentError(lazy"Each interface pin must be a Pair of a key and an endpoint; got $(typeof(pin))."))
        end
        k = pin.first
        if haskey(pinmap, k)
            throw(ArgumentError(lazy"The interface pin key $(k) is not unique."))
        end
        push!(pins, Pair{Any,Any}(k, pin.second))
        pinmap[k] = length(pins)
    end
    isempty(pins) && throw(ArgumentError("An Interface must expose at least one pin."))
    portkeys = Any[]
    portmap = Dict{Any,Tuple{Int,Int}}()
    if !isnothing(interface.ports)
        n = knownlength(interface.ports)
        sizehint!(portkeys, n)
        sizehint!(portmap, n)
        for port in interface.ports
            if !(port isa Pair)
                throw(ArgumentError(lazy"Each interface port must be a Pair of a key and a (positive, negative) tuple of pin keys; got $(typeof(port))."))
            end
            k = port.first
            if haskey(portmap, k)
                throw(ArgumentError(lazy"The interface port key $(k) is not unique."))
            end
            v = port.second
            if !(v isa Tuple && length(v) == 2)
                throw(ArgumentError(lazy"The interface port $(k) must map to a (positive, negative) tuple of pin keys; got $(v)."))
            end
            pk, nk = v
            pos = get(pinmap, pk, 0)
            if pos == 0
                throw(ArgumentError(lazy"The positive side of interface port $(k) is $(pk), which is not a pin key."))
            end
            neg = nk === Ground ? 0 : get(pinmap, nk, 0)
            if nk !== Ground && neg == 0
                throw(ArgumentError(lazy"The negative side of interface port $(k) is $(nk), which is neither a pin key nor Ground."))
            end
            push!(portkeys, k)
            portmap[k] = (pos, neg)
        end
    end
    return InterfaceIndex(pins, pinmap, portkeys, portmap)
end

# === component table ===

# The instances of one circuit level: their identifiers, their definitions
# (with `Instance` wrappers removed), and a lookup from identifier to
# position.
# The identifiers and the component definitions are held as they came,
# heterogeneous, so that the parsed representation encodes neither the
# topology nor the mixture of models in its type; the connectivity is
# resolved to integer tuples in `ParsedLevel` below.
struct ComponentTable
    ids::Vector{Any}
    defs::Vector{Any}
    index::Dict{Any,Int}
    # the interface indexes of the subcircuits, keyed by the subcircuit
    # object, built within one parse or elaboration and dropped with it,
    # so that an edit of a user's collection is seen by the next
    interfacecache::IdDict{Any,InterfaceIndex}
end

function componenttable(components, interfacecache = IdDict{Any,InterfaceIndex}())
    ids = Vector{Any}()
    defs = Vector{Any}()
    index = Dict{Any,Int}()
    names = Set{String}()
    n = knownlength(components)
    sizehint!(ids, n)
    sizehint!(defs, n)
    sizehint!(index, n)
    sizehint!(names, n)
    for entry in components
        if !(entry isa Pair)
            throw(ArgumentError(lazy"Each element of components must be a Pair of an identifier and a component model; got $(typeof(entry))."))
        end
        id = entry.first
        def = entry.second
        if !(id isa Symbol || id isa AbstractString || id isa Integer)
            throw(ArgumentError(lazy"Instance identifiers must be symbols, strings, or integers; got $(typeof(id)) for $(id)."))
        end
        name = string(id)
        if occursin('/', name)
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
        # identifiers become path strings, so `:R1` and `"R1"` are one
        # component further down and are refused here where both are visible
        if name in names
            throw(ArgumentError(lazy"Instance identifier $(id) is not unique: another identifier of a different type has the same name $(name)."))
        end
        push!(names, name)
        push!(ids, id)
        push!(defs, def)
        index[id] = length(ids)
    end
    return ComponentTable(ids, defs, index, interfacecache)
end

function interfaceindex(table::ComponentTable, circuit::Circuit)
    return get!(() -> InterfaceIndex(circuit.interface), table.interfacecache, circuit)
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

function scalarterminal(def::Union{ScatteringParameters,LinearizedScattering,GaussianChannel}, id, k)
    if isgrounded(def)
        if !(k isa Integer) || !(1 <= k <= componentnports(def))
            throw(ArgumentError(lazy"The grounded multiport $(id) has ports 1:$(componentnports(def)); got $(k)."))
        end
        return 2*(Int(k)-1) + 1 # signal terminal of port k
    end
    throw(ArgumentError(lazy"Port $(k) of the multiport $(id) is a two terminal port; write ($(repr(id)), $(k), 1) or ($(repr(id)), $(k), 2) for a single terminal, bond port to port with pair syntax, or construct the block with grounded = true to tie all reference terminals to Ground."))
end

function scalarterminal(table::ComponentTable, i::Int, id, k)
    def = table.defs[i]
    def isa Circuit || return scalarterminal(def, id, k)
    index = interfaceindex(table, def)
    t = get(index.pinmap, k, 0)
    t != 0 && return t
    throw(ArgumentError(lazy"The subcircuit $(id) has no pin $(k). Its pins are $(first.(index.pins))."))
end

# Resolve a (port, terminal) scalar address.
function portterminal(def, id, p, t)
    throw(ArgumentError(lazy"The instance $(id) has no ports; address its terminals as ($(repr(id)), terminal)."))
end

function portterminal(def::Union{ScatteringParameters,LinearizedScattering,GaussianChannel}, id, p, t)
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

function portview(def::Union{ScatteringParameters,LinearizedScattering,GaussianChannel}, id, p)
    np = componentnports(def)
    if !(p isa Integer) || !(1 <= p <= np)
        throw(ArgumentError(lazy"The multiport $(id) has ports 1:$(np); got port $(p)."))
    end
    return (2*(Int(p)-1) + 1, 2*(Int(p)-1) + 2)
end

function portview(table::ComponentTable, i::Int, id, p)
    def = table.defs[i]
    def isa Circuit || return portview(def, id, p)
    if isnothing(def.interface.ports)
        throw(ArgumentError(lazy"The subcircuit $(id) exposes no ports, so $(p) cannot be used as a port in a pair connection."))
    end
    index = interfaceindex(table, def)
    pos, neg = get(index.portmap, p, (0, 0))
    if pos != 0
        return (pos, neg == 0 ? Ground : neg)
    end
    throw(ArgumentError(lazy"The subcircuit $(id) has no port $(p). Its ports are $(index.portkeys)."))
end

# === connection normalization ===

# one electrical net as the parser hands it to the elaborator: its name,
# its endpoints as integer instance and terminal pairs, and whether it is
# tied to ground
struct ParsedGroup
    name::Union{Nothing,String}
    endpoints::Vector{Tuple{Int,Int}}
    hasground::Bool
end

"""
    ParsedLevel

The normalized form of one circuit level, what the elaboration reads: the
component table, the connection groups, the interface pins, the ground
ties and the mutual inductors, each with names resolved to instance
indices and terminals. Built by [`parsecircuitlevel`](@ref).
"""
struct ParsedLevel
    table::ComponentTable
    groups::Vector{ParsedGroup}
    # per interface pin, in interface order: (instanceindex, terminal)
    pinendpoints::Vector{Tuple{Int,Int}}
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
        t = scalarterminal(table, i, ep.instance, ep.key)
        # a declared ground instance is the reference net, not a device:
        # its terminal is the ground sentinel, so every downstream rule
        # (net naming, interface pin restrictions) applies uniformly
        return table.defs[i] isa GroundType ? Ground : (i, t)
    elseif ep isa PortRef
        throw(ArgumentError(lazy"A PortRef is a two terminal port view and cannot be a member of a scalar connection group ($(context)). Use the port in a pair connection or address its terminals individually."))
    elseif ep isa Tuple && length(ep) == 2
        i = instanceindex(table, ep[1], context)
        t = scalarterminal(table, i, ep[1], ep[2])
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
        pos, neg = portview(table, i, ep.instance, ep.key)
        return (:port, attach(i, pos), attach(i, neg))
    elseif ep isa Tuple && length(ep) == 2
        i = instanceindex(table, ep[1], context)
        def = table.defs[i]
        if hasports(def)
            index = def isa Circuit ? interfaceindex(table, def) : nothing
            haspin = !isnothing(index) && haskey(index.pinmap, ep[2])
            hasport = isnothing(index) || haskey(index.portmap, ep[2])
            if haspin && hasport
                throw(ArgumentError(lazy"The key $(ep[2]) of the subcircuit $(ep[1]) exists both as a pin and as a port, which is ambiguous in a pair connection. Use PortRef($(repr(ep[1])), $(repr(ep[2]))) or PinRef($(repr(ep[1])), $(repr(ep[2])))."))
            end
            if !hasport
                # the key is only a pin: scalar fallback
                return (:scalar, resolvescalar(table, ep, context))
            end
            pos, neg = portview(table, i, ep[1], ep[2])
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
    groups = ParsedGroup[]
    sizehint!(groups, knownlength(connections))
    for (ci, entry) in enumerate(connections)
        context = lazy"connection $(ci)"
        if entry isa Net
            name = string(entry.name)
            pushgroup!(groups, name, table, entry.endpoints, context)
        elseif entry isa Pair
            a = resolvepairside(table, entry.first, context)
            b = resolvepairside(table, entry.second, context)
            if a[1] == :scalar && b[1] == :scalar
                addgroup!(groups, nothing, (a[2], b[2]))
            elseif a[1] == :port && b[1] == :port
                addgroup!(groups, nothing, (a[2], b[2]))
                addgroup!(groups, nothing, (a[3], b[3]))
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
    # the endpoints resolved straight into the integer group, without a
    # vector of boxed endpoints in between
    resolved = (resolvescalar(table, ep, context) for ep in endpoints)
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

function addgroup!(groups, name, resolved)
    endpoints = Tuple{Int,Int}[]
    sizehint!(endpoints, knownlength(resolved))
    hasground = false
    for r in resolved
        if r === Ground
            hasground = true
        else
            push!(endpoints, r)
        end
    end
    push!(groups, ParsedGroup(name, endpoints, hasground))
    return nothing
end

# === interface parsing ===

function parseinterface(table::ComponentTable, circuit::Circuit)
    isnothing(circuit.interface) && return Tuple{Int,Int}[]
    index = interfaceindex(table, circuit)
    pinendpoints = Vector{Tuple{Int,Int}}(undef, length(index.pins))
    for (i, (k, target)) in enumerate(index.pins)
        r = resolvescalar(table, target, lazy"interface pin $(k)")
        if r === Ground
            throw(ArgumentError(lazy"The interface pin $(k) maps to Ground. Pins must map to component terminals; use Ground directly in the parent or as the negative side of an interface port."))
        end
        pinendpoints[i] = r
    end
    return pinendpoints
end

# === mutual inductors and grounded ties ===

function parsemutuals(table::ComponentTable)
    mutuals = NTuple{3,Int}[]
    for (i, def) in enumerate(table.defs)
        if def isa MutualInductor
            if isnothing(def.inductor1) || isnothing(def.inductor2)
                throw(ArgumentError(lazy"The mutual inductor $(table.ids[i]) names no inductors. Write MutualInductor(K, inductor1, inductor2) in the connection-group form; MutualInductor(K) alone is for the netlist form, where the entry names the inductors."))
            end
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
        if (def isa ScatteringParameters || def isa LinearizedScattering ||
                def isa GaussianChannel) && isgrounded(def)
            for p in 1:componentnports(def)
                push!(ties, (i, 2*(p-1) + 2))
            end
        end
    end
    return ties
end

# === one level parse (also the constructor validation) ===

"""
    parsecircuitlevel(components, connections, interface)
    parsecircuitlevel(c::Circuit)

Parse and validate one level of a circuit description into a
[`ParsedLevel`](@ref): the component table, the connection groups, the
interface pins, the ground ties and the mutual inductors. The
[`Circuit`](@ref) constructor runs it to validate its arguments, and the
elaboration runs it on every level it flattens.
"""
function parsecircuitlevel(c::Circuit,
        interfacecache = IdDict{Any,InterfaceIndex}())
    table = componenttable(c.components, interfacecache)
    groups = parseconnections(table, c.connections)
    pinendpoints = parseinterface(table, c)
    groundties = parsegroundties(table)
    mutuals = parsemutuals(table)
    return ParsedLevel(table, groups, pinendpoints, groundties, mutuals)
end

parsecircuitlevel(components, connections, interface) =
    parsecircuitlevel(Circuit(components, connections, interface; validate = false))
