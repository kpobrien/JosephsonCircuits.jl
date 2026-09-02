# From a typed `Circuit` to the tables the matrix builders read.
#
# `parseinput.jl` ends with a `Circuit` whose instances may themselves be
# circuits. `elaborate` flattens that hierarchy into an `ElaboratedCircuit`,
# in which every primitive instance has a path and every net is a single
# integer, and `compile` lowers the result to a `CompiledCircuit`: a flat
# table of two terminal components, groups of indices into it by component
# kind, and the ports and scattering blocks as their own records.

# === flattening the hierarchy ===

# A union-find over wire indices. Every terminal of every primitive
# instance starts on its own wire, and each connection group unions the
# wires of its endpoints; the roots at the end are the nets.

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

The flattened result of [`elaborate`](@ref): the hierarchy resolved to a
list of primitive instances, definitions deduplicated by identity, and
nets numbered densely with the ground net first.

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
    # the names the user gave to nets, as (depth, sequence, qualified
    # name, wire); the shallowest and then earliest one wins for a net
    usernames::Vector{Tuple{Int,Int,String,Int}}
    # What an unnamed net is named from: the hierarchy path of the level
    # containing the shallowest, earliest terminal on it. `levelpaths` holds
    # the path of each level visited, and `autodepth` and `autopathid` give
    # the depth and level of each primitive instance, which its terminals
    # share. The terminal index orders candidates at equal depth.
    levelpaths::Vector{String}
    autodepth::Vector{Int}                         # per instance
    autopathid::Vector{Int}                        # per instance
    # each distinct circuit definition is parsed once
    parsedcache::IdDict{Any,ParsedLevel}
    # the definitions on the current path from the root, to detect recursion
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

# the index of `def` in the deduplicated definition list, adding it if new
function definitionindex!(st::FlattenState, def)
    i = get(st.defindex, def, 0)
    if i == 0
        push!(st.definitions, def)
        i = length(st.definitions)
        st.defindex[def] = i
    end
    return i
end

# hierarchical paths use "/" as the separator; the top level has the empty path
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

A repeated subcircuit definition is parsed once and the parse reused for
every instance, so the work per instance is proportional to the instance's
own size. `maxdepth` bounds the nesting depth.
"""
function elaborate(circuit::Circuit; maxdepth::Integer = 64)
    st = FlattenState(Int(maxdepth))
    flattencircuit!(st, circuit, "", 0)
    return finishelaboration(st)
end

# Flatten one level of the hierarchy at `path`, returning the wires of its
# interface pins so the parent can connect them.
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

    # the wires of each local instance's terminals; for a subcircuit, the
    # wires of its interface pins
    instwires = Vector{Vector{Int}}(undef, length(table.ids))
    # the flattened index of each local primitive instance; 0 for a
    # subcircuit or a ground instance, which contribute none
    localglobal = zeros(Int, length(table.ids))

    for (i, def) in enumerate(table.defs)
        id = table.ids[i]
        if def isa Circuit
            instwires[i] = flattencircuit!(st, def, joinpath_(path, id),
                depth + 1)
        elseif def isa GroundType
            # a declared ground instance is a spelling of the reference net,
            # not a device: it produces no instance, and every reference to
            # its terminal was already resolved to `Ground` by the parser,
            # so its (empty) wire list is never read
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

    # the reference terminals of grounded multiport blocks are tied to ground
    for (i, t) in pd.groundties
        unionwires!(st.wires, instwires[i][t], st.groundwire)
    end

    # each connection group unions the wires of its endpoints, and records
    # its name as a candidate name for the resulting net
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

    # mutual inductor couplings, resolved to flattened instance indices
    for (k, i1, i2) in pd.mutuals
        push!(st.couplings, (localglobal[k], localglobal[i1],
            localglobal[i2]))
    end

    # the wires of the interface pins, for the parent to connect
    pinwires = Vector{Int}(undef, length(pd.pinendpoints))
    for (j, (i, t)) in enumerate(pd.pinendpoints)
        pinwires[j] = instwires[i][t]
    end

    delete!(st.active, c)
    return pinwires
end

# Number the nets and name them, and assemble the `ElaboratedCircuit`.
function finishelaboration(st::FlattenState)
    # number the union-find roots densely, in order of first appearance
    # over the terminals; the ground net is 1
    netofroot = Dict{Int,Int}()
    netofroot[findwire(st.wires, st.groundwire)] = 1
    terminalnets = Vector{Int}(undef, length(st.terminalwires))
    # for each net, the (depth, terminal index, level path) of the
    # shallowest and then earliest terminal on it, which names the net when
    # no user name does
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
        # advance to the instance which owns terminal `i`
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

    # A net named by the user takes the shallowest, then earliest, of its
    # user names. Any other net is named "<level path>/net<k>" from the
    # shallowest, earliest level which touches it, with `k` counting the
    # automatically named nets of that level.
    username = Dict{Int,Tuple{Int,Int,String}}()
    for (depth, seq, name, w) in st.usernames
        n = get(netofroot, findwire(st.wires, w), 0)
        (n == 0 || n == 1) && continue # a net with no terminal, or ground
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
#
# `ElaboratedCircuit` is the last representation whose definitions are
# arbitrary component objects. Everything downstream works from the
# `CompiledCircuit`: a flat table of two terminal components which the
# matrix builders, the netlist export and the sensitivities walk entry by
# entry, and index groups by component kind which the assembly plans and
# the ports read without scanning the table.
#
# Two things are not entries in the table. A port's entry holds its
# reference impedance (so that the table's element type is a quantity, not
# a label); the port number, its nodes, and the index of the termination it
# owns are on a `CompiledPort`. A scattering block is not a two terminal
# component and has no table entry at all; it is a
# `CompiledScatteringBlock`.

"""
    CompiledPort

An analysis port: its number, its two nodes, its reference impedance and
the termination it owns.

`environment` is the flat table index of the port's own termination, or
`0` when the port owns none. It is recorded here when the port is
compiled, so nothing downstream needs to look for a resistor on the port's
branch, and a port may share its terminals with any number of ordinary
device resistors.

`zref` is the reference impedance as written, which may still be symbolic;
binding resolves it to a finite positive real number. `component` is the
port's own index in the flat table.
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

A multiport [`ScatteringParameters`](@ref) instance after compilation:
its `definition`, its instance `path`, and the node of the signal and
reference terminal of each port, `signalnodes[p]` and `refnodes[p]`. A
block has no entries in the flat component table.
"""
struct CompiledScatteringBlock
    definition::Any
    signalnodes::Vector{Int}
    refnodes::Vector{Int}
    path::String
end

"""
    CompiledCircuit

An elaborated circuit lowered to a flat table of two terminal components,
with index groups by component kind.

# Fields

The flat table, in elaboration order:

- `componentnames`: the hierarchical instance path of each entry. A matched
    port's own termination is the entry named `"<port path>/termination"`.
- `componenttypes`: the type symbol of each entry: `:C`, `:R`, `:L`, `:Lj`
    (a sinusoidal [`NonlinearInductor`](@ref)), `:NL` (a legacy nonlinear
    element), `:I`, `:K` (a mutual inductor) or `:P` (a port).
- `componentvalues`: the value of each entry as written; the reference
    impedance for a port.
- `nodeindices`: a 2 by `ncomponents` matrix of the node indices of each
    entry, ground being node 1; both zero for a mutual inductor.
- `componenttemperatures`: the temperature of each entry which states one,
    keyed by flat index.
- `mutualinductorbranchnames`: the names of the coupled inductors, two per
    `:K` entry in order.
- `nodenames`, `Nnodes`: the node names in sorted order (ground first) and
    their count.
- `componentnamedict`: component name to flat index.

The groups, each a vector of flat indices in table order:

- `capacitors`, `resistors`, `inductors`, `junctions` (`:Lj`),
  `nonlinearinductors` (`:NL`), `currentsources`, `mutualinductors`.

The records which keep their own structure:

- `ports::Vector{CompiledPort}`, in elaboration order.
- `scatteringblocks::Vector{CompiledScatteringBlock}`, one per block
  instance.

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

# A tuple netlist names its nodes with integers and is sorted by their
# numeric value; a typed circuit has hierarchical net names, which are not
# integers, and is sorted by name.
defaultsorting(circuit) = circuit isa AbstractVector ? :number : :name

# === lowering one component to a table entry ===
#
# `lowercomponent` returns the `(typesymbol, value)` of a component, or
# throws `ComponentNotSupportedError` for one the solvers cannot use. The
# `isa` chain in `compile` handles the common components inline and falls
# through to `lowercomponent` for the rest, so the diagnostics are in one
# place. The legacy `NL` element adds its own method in legacyadapter.jl.

function lowercomponent(def::Inductor, path)
    return :L, def.L
end
lowercomponent(def::Capacitor, path) = :C, def.C
lowercomponent(def::Resistor, path) = :R, def.R
lowercomponent(def::CurrentSource, path) = :I, def.I
lowercomponent(def::Port, path) = :P, def.Z0
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

# Narrow a `Vector{Any}` to the element type its contents allow, so that a
# fully numeric circuit gets a `Vector{Float64}` of values rather than a
# vector of boxed numbers.
function tightenvalues(values::Vector{Any})
    return map(identity, values)
end

# The temperature a component states, or `nothing`. Only the lumped
# components which can dissipate carry one; a scattering block states its
# temperature through its noise model (see `ThermalEquilibrium`).
componenttemperature(def::Resistor) = def.temperature
componenttemperature(def::Capacitor) = def.temperature
componenttemperature(def::Inductor) = def.temperature
componenttemperature(def) = nothing

"""
    compile(elab::ElaboratedCircuit; sorting = :name)
    compile(circuit::Circuit; sorting = :name)
    compile(netlist::AbstractVector; sorting = :name)

Lower a circuit to a [`CompiledCircuit`](@ref). A [`Circuit`](@ref) is
elaborated first, and a legacy tuple netlist is converted to a `Circuit`
first; a `CompiledCircuit` is returned unchanged.

Components appear in the table in elaboration order, with a matched port's
own termination emitted as a resistor entry directly after the port. Nodes
are numbered by [`calcnodesorting`](@ref) with ground first; the default
`sorting = :name` sorts the net names as strings, since hierarchical net
names are not integers (the tuple netlist entry points of the solvers
default to `:number` instead).

Only components the solvers support can be lowered: a
[`GaussianChannel`](@ref), a [`VoltageSource`](@ref), a non-sinusoidal
[`NonlinearInductor`](@ref), a [`ScatteringParameters`](@ref) with a
[`NoiseCovariance`](@ref) noise model, or any component with other than
two terminals throws a [`ComponentNotSupportedError`](@ref) naming the
instance. A circuit with no connection to [`Ground`](@ref) throws an
`ArgumentError`.
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

    # the names of the two inductors coupled by each mutual inductor
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

        # A scattering block has two terminals per port and is not a two
        # terminal component, so it gets no table entry: it is compiled as
        # one `CompiledScatteringBlock` holding the nodes of every port.
        if def isa ScatteringParameters
            if def.noise isa NoiseCovariance
                throw(ComponentNotSupportedError(lazy"the ScatteringParameters at $(path) has an arbitrary noise covariance, which the harmonic balance solvers do not yet support. Use Passive, Lossless, or ThermalEquilibrium."))
            end
            terminals = instanceterminals(elab, i)
            n = def.nports
            signalnodes = Vector{Int}(undef, n)
            refnodes = Vector{Int}(undef, n)
            for p in 1:n
                signalnodes[p] = processnode(uniquenodedict, uniquenodevector,
                    elab.netnames[terminals[2*p-1]])
                refnodes[p] = processnode(uniquenodedict, uniquenodevector,
                    elab.netnames[terminals[2*p]])
            end
            push!(scatteringblocks, CompiledScatteringBlock(def, signalnodes,
                refnodes, path))
            continue
        end

        # The common components are lowered by an `isa` chain, ordered by
        # how many of each a large circuit typically holds, and everything
        # else falls through to `lowercomponent`. On a heterogeneous vector
        # a branch chain avoids dynamic dispatch, although the interning of
        # node names below dominates this loop either way.
        typesymbol, value = if def isa Capacitor
            (:C, def.C)
        elseif def isa NonlinearInductor
            issinusoidal(def) ? (:Lj, def.L0) : lowercomponent(def, path)
        elseif def isa Inductor
            (:L, def.L)
        elseif def isa Resistor
            (:R, def.R)
        elseif def isa Port
            # a port's table value is its reference impedance, a quantity of
            # the same kind as the other values, so that the table's element
            # type stays concrete; the port number is on the `CompiledPort`
            (:P, def.Z0)
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
        # a component which states its own temperature records it; the rest
        # take the temperature the analysis is run at
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
            # A matched port owns its external environment, emitted here as
            # an ordinary resistor entry across the port's nodes, which the
            # conductance stamping and the solver scale consume like any
            # other resistor. What marks it as the port's own is the index
            # recorded on the `CompiledPort`, so any further resistor on the
            # same nodes is a device resistor.
            #
            # The generated name "<path>/termination" cannot collide with an
            # instance: that would require the instance at "<path>" to be a
            # subcircuit containing an instance named "termination", and it
            # is a Port.
            environment = 0
            if def.termination isa MatchedTermination
                push!(componentnames, path * "/termination")
                push!(componenttypes, :R)
                push!(componentvalues, def.Z0)
                push!(nodeindexvector, n1)
                push!(nodeindexvector, n2)
                environment = length(componentnames)
            elseif def.termination isa LegacyTermination
                # the termination is a resistor which already exists in the
                # table. Its name is relative to the port's level, so a
                # legacy circuit instanced as a subcircuit still resolves;
                # the index is looked up once the table is complete.
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

    nodenames, nodeindices, renumber = sortnodes(uniquenodevector,
        nodeindexvector, calcnodesorting(uniquenodevector; sorting = sorting))

    componentnamedict = Dict{String,Int}()
    sizehint!(componentnamedict, length(componentnames))
    for (i, name) in enumerate(componentnames)
        componentnamedict[name] = i
    end

    # resolve the legacy terminations now that every name is in the table
    for (k, name) in legacyenvironments
        i = get(componentnamedict, name, 0)
        if iszero(i) || componenttypes[i] !== :R
            throw(ArgumentError(lazy"The port $(componentnames[ports[k].component]) names $(name) as its termination, which is not a resistor in this circuit."))
        end
        ports[k] = CompiledPort(ports[k].number, ports[k].positivenode,
            ports[k].negativenode, ports[k].zref, i, ports[k].component)
    end

    # `sortnodes` renumbered the nodes. The port nodes recorded above are in
    # the pre-sort numbering and are re-read from the sorted table.
    ports = [CompiledPort(p.number, nodeindices[1, p.component],
        nodeindices[2, p.component], p.zref, p.environment, p.component)
        for p in ports]

    warnduplicatematchedload(ports, componentnames, componenttypes,
        componentvalues, nodeindices)

    # the block nodes are also pre-sort and are translated with the
    # renumbering, since blocks have no table entry to re-read them from
    for b in scatteringblocks
        for p in eachindex(b.signalnodes)
            b.signalnodes[p] = renumber[b.signalnodes[p]]
            b.refnodes[p] = renumber[b.refnodes[p]]
        end
    end

    group(t) = [i for (i, s) in enumerate(componenttypes) if s === t]

    return CompiledCircuit(nodenames, nodeindices, length(uniquenodevector),
        componentnames, componenttypes, tightenvalues(componentvalues),
        componentnamedict, componenttemperatures, mutualinductorbranchnames,
        group(:C), group(:R), group(:L), group(:Lj), group(:NL), group(:I),
        group(:K), ports, scatteringblocks)
end

"""
    warnduplicatematchedload(ports, componentnames, componenttypes,
        componentvalues, nodeindices)

Warn when a matched port has a device resistor of exactly its own reference
impedance across the same two nodes.

Such a circuit is legal, and is what a user who wants two loads means, so
it is not refused. It is far more often a circuit written in the tuple
netlist style, where the resistor across a port *was* its termination, and
now carries two loads instead of one. Only an exact match of the reference
impedance is reported, because that is what makes the resistor a likely
duplicate rather than a device.
"""
function warnduplicatematchedload(ports, componentnames, componenttypes,
        componentvalues, nodeindices)
    isempty(ports) && return nothing
    for p in ports
        iszero(p.environment) && continue
        z = componentvalues[p.environment]
        z isa Number || continue
        for i in eachindex(componenttypes)
            componenttypes[i] === :R || continue
            i == p.environment && continue
            n1, n2 = nodeindices[1, i], nodeindices[2, i]
            ((n1 == p.positivenode && n2 == p.negativenode) ||
             (n1 == p.negativenode && n2 == p.positivenode)) || continue
            v = componentvalues[i]
            (v isa Number && v == z) || continue
            @warn "This port owns a matched environment of its own and a device resistor of the same value sits across the same terminals, so the port is loaded twice. If the resistor was written as the port's termination, which is how a port was terminated before a port could own one, either delete it or write the port as `termination = nothing` to keep it as the only load. If two loads are intended, this is correct and the warning can be ignored." port=p.number resistor=componentnames[i] value=z
        end
    end
    return nothing
end

compile(circuit::Circuit; sorting::Symbol = :name) =
    compile(elaborate(circuit); sorting = sorting)

# a tuple netlist becomes a typed circuit first
compile(netlist::AbstractVector; sorting::Symbol = :name) =
    compile(Circuit(netlist); sorting = sorting)

# a compiled circuit is returned as is, so the solver entry points accept one
compile(c::CompiledCircuit; sorting::Symbol = :name) = c

# === port and noise roles, read from the compiled circuit ===
#
# A compiled port states its own reference impedance and, through
# `environment`, which table entry realizes it. The functions below read
# those records; none of them looks for a resistor on a port's branch. The
# only place a role is recovered from circuit geometry is the legacy
# adapter, where the tuple format has no way to state one.

"""
    scatteringblockindex(c::CompiledCircuit, name)

The position in `c.scatteringblocks` of the block whose instance path is
`name`, or zero when there is none. The spelling `"<path>/port1"` is also
accepted for compatibility with names the design sensitivities once
produced.
"""
function scatteringblockindex(c::CompiledCircuit, name)
    s = String(name)
    for (k, b) in enumerate(c.scatteringblocks)
        (b.path == s || b.path*"/port1" == s) && return k
    end
    return 0
end

"""
    portindicesnumbers(c::CompiledCircuit)

The flat table indices and the numbers of the ports, both ordered by port
number. Throws an `ArgumentError` for duplicate port numbers or two ports
on the same branch.
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

The flat table index of each port's own termination, ordered by port
number, or zero for a port which owns none.

A port owns a termination when it was written with the default
[`MatchedTermination`](@ref), or when the legacy adapter recorded the
resistor a tuple netlist placed across it. A port written with
`termination = nothing` owns none. The index identifies a role rather than
a value: it says which entry realizes the port's reference impedance, which
the noise classification needs (a port termination is an external bath,
not an internal noise channel) and the sensitivities need (perturbing that
entry also moves the wave normalization). The impedance itself comes from
[`portreferenceimpedances`](@ref), which is defined for every port.
"""
portenvironmentindices(c::CompiledCircuit) =
    [c.ports[i].environment
     for i in sortperm([p.number for p in c.ports])]

"""
    portreferenceimpedances(c::CompiledCircuit, values)

The reference impedance of each port, ordered by port number, which the
incoming and outgoing waves are normalized to.

For a port which owns a termination the impedance is read from that entry
of the resolved value table `values`, so that a swept or symbolic
impedance resolves the same way any component value does; the entry holds
the port's own `Z0` in any case. For a port which owns none it is the
port's `zref` as written.
"""
portreferenceimpedances(c::CompiledCircuit, values) =
    [iszero(c.ports[i].environment) ? c.ports[i].zref :
        values[c.ports[i].environment]
     for i in sortperm([p.number for p in c.ports])]

"""
    noiseindices(c::CompiledCircuit, values)

The flat table indices of the internal dissipative components, which are
the noise channels of the linearized analysis: every resistor which is not
a port's own termination, and every capacitor or inductor whose resolved
value in `values` has a nonzero imaginary part.

A port termination is an external bath rather than an internal channel and
is excluded by its role; any other resistor across a port's nodes is an
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
