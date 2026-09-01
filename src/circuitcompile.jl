# Lowering a typed circuit to the tables the assembly reads.
#
# The second half of the input path. `parseinput.jl` ends with a `Circuit`
# whose instances may themselves be circuits; this file flattens that
# hierarchy into an `ElaboratedCircuit`, in which every instance has a path
# and every net is a single wire, and then lowers it to a
# `CompiledCircuit`, the flat integer indexed tables plus the groups the
# assembly reads.
#
# The two steps share a file because `compile` dispatches on what
# `elaborate` returns, and separating them would put a type and its only
# consumer either side of a file boundary for no reader's benefit.

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
# The compiled form is a flat component table plus groups indexing into it,
# and it is honest to say that the table is still where much of the solver
# reads from rather than that it has been retired. The matrix builders, the
# netlist export and the sensitivities walk it component by component. The
# groups are what the assembly plans and the ports read, and they are what
# makes those O(components in the group) rather than a scan.
#
# What the table no longer holds is anything which is not a component. A
# port's entry carries its reference impedance and not its number, so the
# table's element type is the type of the quantities in it; and a scattering
# block, which is not a two terminal component, has no entry there and lives
# in `scatteringblocks` alone.
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
`refnodes[p]` are the nodes of port `p`, and the whole block is one object.
The flat component table can only hold two terminal components, so a block
had to appear there as one entry per port and be reassembled afterwards;
nothing does that now, and the block has no entries in that table at all.
"""
struct CompiledScatteringBlock
    definition::Any
    signalnodes::Vector{Int}
    refnodes::Vector{Int}
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
            # A block is one thing and is compiled as one. Its ports used to
            # be emitted as one flat component each, every one carrying the
            # whole block definition, and the block was then reconstructed
            # downstream by scanning for runs of them -- a contract that the
            # ports appear consecutively with port 1 first, stated nowhere
            # and checked nowhere. The block below already holds everything
            # those entries did, so they are not written.
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
            # the reference impedance, not the port number. A port has one
            # quantity and the number is not it: a label in the value table
            # makes the table's element type the join of a label and a
            # quantity, which for the common circuit is `Number` rather than
            # `ComplexF64`, so every read of every component value is boxed.
            # The number is on the `CompiledPort`, which is where a number
            # that identifies rather than measures belongs.
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

    nodenames, nodeindices, renumber = sortnodes(uniquenodevector,
        nodeindexvector, calcnodesorting(uniquenodevector; sorting = sorting))

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

    warnduplicatematchedload(ports, componentnames, componenttypes,
        componentvalues, nodeindices)

    # the block terminals are node indices from before the sort, and are
    # renumbered like every other one
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

Warn when a matched port has a device resistor of its own value across it.

Before ports owned their environments, a port was terminated by a resistor
the user placed across it, and the two together were one 50 ohm load.
Written unchanged now, the same circuit has two: the port's own matched
environment and the resistor beside it. That is a legal circuit and is
sometimes what is wanted, so it is not refused; but it is much more often a
circuit which was written under the old rule and whose loading has quietly
halved, which is not the kind of change that announces itself in an answer.

Only an exact match of the reference impedance is reported, because that is
what makes it a likely duplicate rather than a device.
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

# a netlist of tuples becomes a typed circuit first
compile(netlist::AbstractVector; sorting::Symbol = :name) =
    compile(Circuit(netlist); sorting = sorting)

# and a circuit which is already compiled is itself, so the solver entry
# points accept one
compile(c::CompiledCircuit; sorting::Symbol = :name) = c

# === port and noise roles, read from the compiled circuit ===
#
# The legacy netlist carries no role, so the adapter recovers one by
# geometry: it looks for a resistor sharing a port's branch and calls that
# the port impedance, which is why the legacy language allows a port exactly
# one resistor across it. That is the only place geometry decides a role. A
# compiled port states its own reference impedance and which component, if
# any, realizes it, so the lists below come out of the structure: a port may
# share its terminals with any number of ordinary device resistors, none of
# which it adopts, and a port which loads the circuit with nothing at all is
# an ordinary thing to write.

"""
    scatteringblockindex(c::CompiledCircuit, name)

The ordinal of the scattering block named `name`, or zero.

A block is named by its instance path, which is the name the user wrote. It
used to be named by the flat component its first port lowered to, so a two
port block at `x` was addressed as `x/port1`; that spelling is still
accepted, because it is what the design sensitivities produced.
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

The flat index of each port's own environment, ordered by port number, or
zero for a port which owns none.

A port owns an environment when it was declared with one, either as the
matched source and load of `termination = nothing`'s complement or as the
resistor a legacy netlist placed across it. `termination = nothing` means
there is no such component, and zero says so. Nothing is searched for here:
a resistor a user placed across a port is that user's device resistor, and
the only place a role is recovered from geometry is the legacy adapter,
where the input language has no way to state one.

The index names a role, not a number. It says which component a port's
reference impedance is realized by, for the two things which need to know
that: the noise classification, which counts an external bath differently
from an internal channel, and the sensitivities, where differentiating that
component also moves the wave normalization. The impedance itself comes from
[`portreferenceimpedances`](@ref), which is defined for every port.
"""
portenvironmentindices(c::CompiledCircuit) =
    [c.ports[i].environment
     for i in sortperm([p.number for p in c.ports])]

"""
    portreferenceimpedances(c::CompiledCircuit, values)

The reference impedance of each port, ordered by port number, resolved
against a flat value table.

This is the impedance the incoming and outgoing waves are normalized to. It
is the port's declared `Z0`, which is what the reference impedance means;
when the port owns an environment the value is read from the table instead,
so that a swept or symbolic impedance resolves the same way a component
value does and the two cannot disagree. They describe the same quantity: a
matched environment is generated with the port's own `Z0`, and a legacy
port's `Z0` is the value of the resistor it adopted.
"""
portreferenceimpedances(c::CompiledCircuit, values) =
    [iszero(c.ports[i].environment) ? c.ports[i].zref :
        values[c.ports[i].environment]
     for i in sortperm([p.number for p in c.ports])]

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
