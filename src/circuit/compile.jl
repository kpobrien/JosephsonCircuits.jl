# From a typed `Circuit` to the tables the matrix builders read.
#
# `circuit/parse.jl` ends with a `Circuit` whose instances may themselves be
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
    # the names the user gave to nets, as (depth, qualified name, wire);
    # the shallowest and then earliest entry wins for a net
    usernames::Vector{Tuple{Int,String,Int}}
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
    interfacecache::IdDict{Any,InterfaceIndex}
    # the definitions on the current path from the root, to detect recursion
    active::IdDict{Any,Nothing}
    maxdepth::Int
end

function FlattenState(maxdepth::Int)
    wires = WireForest()
    ground = newwire!(wires)
    return FlattenState(wires, ground, Any[], IdDict{Any,Int}(), Int[],
        String[], Int[1], Int[], NTuple{3,Int}[],
        Tuple{Int,String,Int}[], String[], Int[], Int[],
        IdDict{Any,ParsedLevel}(), IdDict{Any,InterfaceIndex}(),
        IdDict{Any,Nothing}(), maxdepth)
end

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

# Flatten one level of the hierarchy at `path`, appending the wires of its
# interface pins to `pinwires` so the parent can connect them.
function flattencircuit!(st::FlattenState, c::Circuit, path::String,
        depth::Int, pinwires::Vector{Int} = Int[])
    if haskey(st.active, c)
        location = isempty(path) ? "the top level" : path
        throw(ArgumentError(lazy"The circuit definition at $(location) contains itself, directly or indirectly. Recursive circuit definitions are not allowed."))
    end
    if depth > st.maxdepth
        throw(ArgumentError(lazy"The circuit hierarchy at $(path) exceeds the maximum depth $(st.maxdepth). Pass a larger maxdepth to elaborate if this is intentional."))
    end
    st.active[c] = nothing

    pd = get!(() -> parsecircuitlevel(c, st.interfacecache), st.parsedcache, c)
    flattenlevel!(st, pd, path, depth, pinwires)
    delete!(st.active, c)
    return pinwires
end

# The walk of one level reads the parsed connectivity alone, so it is
# compiled once, whatever the container and interface types of the
# circuits it came from.
@noinline function flattenlevel!(st::FlattenState, pd::ParsedLevel,
        path::String, depth::Int, pinwires::Vector{Int})
    table = pd.table
    push!(st.levelpaths, path)
    pathid = length(st.levelpaths)

    # the wires of each local instance's terminals; for a subcircuit, the
    # wires of its interface pins
    # the wires of every instance's terminals in one array with an offset
    # per instance, the children's interface wires appended as they are
    # flattened, rather than a two-element array per lumped component
    instwires = Int[]
    sizehint!(instwires, 2*length(table.ids))
    instoffsets = Vector{Int}(undef, length(table.ids))
    wireat((i, t)) = instwires[instoffsets[i] + t - 1]
    # the flattened index of each local primitive instance; 0 for a
    # subcircuit or a ground instance, which contribute none
    localglobal = zeros(Int, length(table.ids))

    for (i, def) in enumerate(table.defs)
        id = table.ids[i]
        instoffsets[i] = length(instwires) + 1
        if def isa Circuit
            flattencircuit!(st, def, joinpath_(path, id),
                depth + 1, instwires)
        elseif def isa GroundType
            # a declared ground instance is a spelling of the reference net,
            # not a device: it produces no instance, and every reference to
            # its terminal was already resolved to `Ground` by the parser,
            # so its empty range of local wires is never read
        else
            n = nterminals(def)
            for t in 1:n
                w = newwire!(st.wires)
                push!(instwires, w)
                push!(st.terminalwires, w)
            end
            push!(st.autodepth, depth)
            push!(st.autopathid, pathid)
            push!(st.definitionof, definitionindex!(st, def))
            push!(st.instancepaths, joinpath_(path, id))
            push!(st.terminaloffsets, st.terminaloffsets[end] + n)
            localglobal[i] = length(st.definitionof)
        end
    end

    # the reference terminals of grounded multiport blocks are tied to ground
    for endpoint in pd.groundties
        unionwires!(st.wires, wireat(endpoint), st.groundwire)
    end

    # each connection group unions the wires of its endpoints, and records
    # its name as a candidate name for the resulting net
    for group in pd.groups
        first = group.hasground ? st.groundwire : wireat(group.endpoints[1])
        for endpoint in group.endpoints
            unionwires!(st.wires, first, wireat(endpoint))
        end
        if !isnothing(group.name)
            push!(st.usernames, (depth, joinpath_(path, group.name),
                first))
        end
    end

    # mutual inductor couplings, resolved to flattened instance indices
    for (k, i1, i2) in pd.mutuals
        push!(st.couplings, (localglobal[k], localglobal[i1],
            localglobal[i2]))
    end

    # the wires of the interface pins, for the parent to connect
    for endpoint in pd.pinendpoints
        push!(pinwires, wireat(endpoint))
    end

    return pinwires
end

# Number the nets and name them, and assemble the `ElaboratedCircuit`.
function finishelaboration(st::FlattenState)
    # number the union-find roots densely, in order of first appearance
    # over the terminals; the ground net is 1. The wires are numbered
    # densely, so the net of a root is an array over them; the naming
    # metadata below is per net.
    netofroot = zeros(Int, length(st.wires.parent))
    netofroot[findwire(st.wires, st.groundwire)] = 1
    terminalnets = Vector{Int}(undef, length(st.terminalwires))
    # the depth and the level path of the terminal which names each net,
    # the first at the shallowest depth, the terminals being visited in
    # order
    unset = (typemax(Int), 0)
    autopath = Tuple{Int,Int}[unset]
    nnets = 1
    inst = 1
    ninst = length(st.autodepth)
    for (i, w) in enumerate(st.terminalwires)
        r = findwire(st.wires, w)
        n = netofroot[r]
        if n == 0
            nnets += 1
            n = nnets
            netofroot[r] = n
            push!(autopath, unset)
        end
        terminalnets[i] = n
        # advance to the instance which owns terminal `i`
        while inst < ninst && i >= st.terminaloffsets[inst+1]
            inst += 1
        end
        if n != 1 && inst <= ninst
            d = st.autodepth[inst]
            if d < autopath[n][1]
                autopath[n] = (d, st.autopathid[inst])
            end
        end
    end

    # A net named by the user takes the shallowest, then earliest, of its
    # user names. Any other net is named "<level path>/net<k>" from the
    # shallowest, earliest level which touches it, with `k` counting the
    # automatically named nets of that level.
    # the record which names each net, by index rather than by a copy of
    # its name
    username = zeros(Int, nnets)
    for (j, (depth, name, w)) in enumerate(st.usernames)
        n = netofroot[findwire(st.wires, w)]
        (n == 0 || n == 1) && continue # a net with no terminal, or ground
        winner = username[n]
        if winner == 0 || depth < st.usernames[winner][1]
            username[n] = j
        end
    end
    netnames = Vector{String}(undef, nnets)
    netnames[1] = "0"
    # a counter of automatic names per level, whose instance path is unique
    autocounter = zeros(Int, length(st.levelpaths))
    for n in 2:nnets
        winner = username[n]
        if winner != 0
            netnames[n] = st.usernames[winner][2]
        else
            pathid = autopath[n][2]
            autocounter[pathid] += 1
            netnames[n] = joinpath_(st.levelpaths[pathid], "net$(autocounter[pathid])")
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

An analysis port and the environment it owns.

`environment` is the flat table index of the port's own termination, or
`0` when the port owns none. It is recorded here when the port is
compiled, so nothing downstream needs to look for a resistor on the port's
branch, and a port may share its terminals with any number of ordinary
device resistors.

The reference impedance is not stored here: it is the value of the port's
own entry in the flat component table, at `component`, so it is bound like
every other value and read from the bound table by
[`portreferenceimpedances`](@ref).
"""
struct CompiledPort
    number::Int
    positivenode::Int
    negativenode::Int
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
    CircuitTopology(edge2indexdict, Rbn, Nbranches)

The branches of a compiled circuit, which is what the solvers read of its
graph.

# Fields
- `edge2indexdict`: maps a branch `(node1, node2)` in either orientation to
    its branch index, the row of `Rbn` it occupies.
- `Rbn`: the sparse oriented incidence matrix, `Nbranches` by
    `Nnodes - 1`; the ground node column is omitted.
- `Nbranches`: the number of branches, `size(Rbn, 1)`.

[`compile`](@ref) builds one for every circuit, see
[`circuittopology`](@ref). The spanning tree, the loops and the isolated
nodes are diagnostics of the same branches and are computed on request by
[`calccircuitgraph`](@ref).
"""
struct CircuitTopology
    edge2indexdict::Dict{Tuple{Int,Int},Int}
    Rbn::SparseMatrixCSC{Int,Int}
    Nbranches::Int
end

"""
    circuittopology(componenttypes, nodeindices, Nnodes)
    circuittopology(branchvector, Nnodes)

The [`CircuitTopology`](@ref) of the branch carrying components listed by
`componenttypes` and `nodeindices`, or of the branches `branchvector`
given as `(node1, node2)` tuples, over `Nnodes` nodes with ground being
node 1; see [`extractbranches`](@ref) for which components make a branch.

The branches are the edges of the undirected graph of those endpoints, in
ascending order of their endpoints, oriented from the lower to the higher
node index. A component whose two terminals are the same node is a branch
with no incidence entries: it carries its value and couples no node.
"""
circuittopology(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    Nnodes::Int) = circuittopology(
        extractbranches(componenttypes, nodeindices), Nnodes)

function circuittopology(branchvector::Vector{Tuple{Int,Int}}, Nnodes::Int)
    gl = Graphs.SimpleGraphFromIterator(tuple2edge(branchvector))
    edge2indexdict = Dict{Tuple{Int,Int},Int}()
    I = Int[]; J = Int[]; V = Int[]
    Nbranches = Graphs.ne(gl)
    sizehint!(edge2indexdict, 2Nbranches)
    sizehint!(I, 2Nbranches); sizehint!(J, 2Nbranches); sizehint!(V, 2Nbranches)
    for (i, edge) in enumerate(Graphs.edges(gl))
        a, b = Graphs.src(edge), Graphs.dst(edge)
        edge2indexdict[(a,b)] = i
        edge2indexdict[(b,a)] = i
        a == b && continue
        if a > 1
            push!(I, i); push!(J, a-1); push!(V, -1)
        end
        if b > 1
            push!(I, i); push!(J, b-1); push!(V, 1)
        end
    end
    return CircuitTopology(edge2indexdict,
        sparse(I, J, V, Nbranches, Nnodes-1), Nbranches)
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
    (a sinusoidal [`NonlinearInductor`](@ref)), `:I`, `:K` (a mutual
    inductor) or `:P` (a port).
- `componentvalues`: the value of each entry as written; the reference
    impedance for a port.
- `nodeindices`: a 2 by `ncomponents` matrix of the node indices of each
    entry, ground being node 1; both zero for a mutual inductor.
- `junctioncprs`: the [`PolynomialCPR`](@ref) of each `:Lj` entry whose
    current-phase relation is not the sinusoidal Josephson one. Empty for
    every circuit which does not ask for another, and the solvers then
    evaluate `sin` and `cos`.
- `componenttemperatures`: the temperature of each entry which states one,
    keyed by flat index.
- `couplings`: resolved `(coupling, inductor1, inductor2)` flat indices,
    in coupling table order; the names of the coupled inductors are
    `componentnames` at those indices, see
    [`coupledinductornames`](@ref).
- `nodenames`, `Nnodes`: the node names in sorted order (ground first) and
    their count.
- `componentnamedict`: component name to flat index.
- `topology`: the branches and the oriented incidence matrix, see
    [`CircuitTopology`](@ref).

The groups, each a vector of flat indices in table order:

- `capacitors`, `resistors`, `inductors`, `junctions` (`:Lj`),
  `currentsources`, `mutualinductors`.

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
    junctioncprs::Dict{Int,PolynomialCPR{Float64}}
    capacitors::Vector{Int}
    resistors::Vector{Int}
    inductors::Vector{Int}
    junctions::Vector{Int}
    currentsources::Vector{Int}
    mutualinductors::Vector{Int}
    ports::Vector{CompiledPort}
    scatteringblocks::Vector{CompiledScatteringBlock}
    # (coupling, first inductor, second inductor), in flat table order
    couplings::Vector{NTuple{3,Int}}
    topology::CircuitTopology
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

# === lowering one component to a table entry ===
#
# `lowercomponent` returns the `(typesymbol, value)` of a component the
# `isa` chain in `compile` does not handle inline, or throws
# `ComponentNotSupportedError` for one the solvers cannot use, so the
# diagnostics are in one place. The lumped elements are lowered by the
# chain alone; a second table for them would be a second place to get one
# wrong.

function lowercomponent(def::VoltageSource, path)
    throw(ComponentNotSupportedError(lazy"the VoltageSource at $(path) is not supported by the solvers."))
end
function lowercomponent(def::GaussianChannel, path)
    throw(ComponentNotSupportedError(lazy"the GaussianChannel at $(path) is not yet supported by the harmonic balance solvers. It parsed, validated, and elaborated successfully; solver support for Gaussian channels is planned. Currently solvable components: Inductor, Capacitor, Resistor, JosephsonJunction and the other NonlinearInductors, MutualInductor, CurrentSource, Port, and the scattering blocks (ScatteringParameters, TransmissionLine, RationalScattering and LinearizedScattering)."))
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

# === node ordering ===
#
# `compile` numbers the nets in the order their terminals are met, sorts
# that list with `calcnodesorting` and renumbers every recorded node index
# with `sortnodes`.

"""
    findgroundnodeindex(uniquenodevector::Vector{String})

The index of the ground node `"0"` in `uniquenodevector`, or `0` if there
is none.

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
    for i in eachindex(uniquenodevector)
        if uniquenodevector[i] == "0"
            return i
        end
    end
    return 0
end

"""
    calcnodesorting(uniquenodevector::Vector{String};sorting=:number)

The permutation which sorts the node names in `uniquenodevector` according
to `sorting`, with the ground node `"0"` moved to the front in every case.
Throws an `ArgumentError` if there is no ground node.

# Keywords
- `sorting = :number`: parse the names as integers and sort numerically.
    Throws an `ArgumentError` if a name is not an integer.
- `sorting = :name`: sort the names as strings, so that `"101"` sorts
    before `"11"`.
- `sorting = :none`: keep the names in order of first appearance, apart
    from moving ground to the front.

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

    # the identity permutation, which `:none` keeps
    uniquenodevectorsortindices = Vector{Int}(undef,length(uniquenodevector))
    for i in eachindex(uniquenodevectorsortindices)
        uniquenodevectorsortindices[i] = i
    end

    if sorting == :name
        sortperm!(uniquenodevectorsortindices,uniquenodevector,initialized=true)

    elseif sorting == :number
        uniquenodevectorints = Vector{Int}(undef,length(uniquenodevector))
        for i in eachindex(uniquenodevectorints)
            parsednode = tryparse(Int,uniquenodevector[i])
            if !isnothing(parsednode)
                uniquenodevectorints[i] = parsednode
            else
                throw(ArgumentError(lazy"The node $(repr(uniquenodevector[i])) is not an integer. Name the nodes with integers, or set the keyword argument `sorting=:name` or `sorting=:none`."))
            end
        end
        sortperm!(uniquenodevectorsortindices, uniquenodevectorints, initialized=true)

    elseif sorting == :none
        nothing
    else
        throw(ArgumentError(lazy"Unknown sorting $(repr(sorting)); use :number, :name or :none."))
    end

    groundnodeindex = findgroundnodeindex(uniquenodevector)

    if groundnodeindex == 0
        throw(ArgumentError("The circuit has no connection to Ground. Connect at least one endpoint to Ground; the ground net is required by the solver."))
    end

    # move ground to the front, shifting the nodes which sorted before it
    # back by one
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
    noderenumbering(order)

The renumbering induced by the sorting permutation `order` returned by
[`calcnodesorting`](@ref): `renumber[j]` is the new index of the node whose
old index was `j`. `compile` uses it to renumber the node indices of
scattering blocks, which have no component table entry to be re-read from.
"""
noderenumbering(order::Vector{Int}) = invperm(order)

"""
    sortnodes(uniquenodevector, nodeindexvector, order)

Apply the precomputed sorting permutation `order` (see
[`calcnodesorting`](@ref)), returning the sorted names, the renumbered
component node indices as a 2 by `Ncomponents` matrix, and the
renumbering itself (see [`noderenumbering`](@ref)).
"""
function sortnodes(uniquenodevector::Vector{String},
        nodeindexvector::Vector{Int}, order::Vector{Int})

    nodeindices = zeros(eltype(nodeindexvector),2,length(nodeindexvector)÷2)

    nodevectorsortindices = noderenumbering(order)

    for (i,j) in enumerate(nodeindexvector)
        # a mutual inductor couples two inductors rather than two nodes, so
        # its node indices are zero and stay zero
        if j == 0
            nothing
        else
            nodeindices[i] = nodevectorsortindices[j]
        end
    end

    return uniquenodevector[order], nodeindices, nodevectorsortindices
end

"""
    shortednets(elab::ElaboratedCircuit)

The net each net of `elab` is compiled as: the ground net (1) for a net
whose every terminal belongs to a component, or a port of a scattering
block, shorted across that net, and the net itself otherwise.

A component whose two terminals are one node carries no current and
couples no node, so a net which only such components touch is a node with
no equation of its own: every frequency domain solve would meet it as an
empty row. On the ground net the components are the self loops they are
at any other node, and the node is gone.
"""
function shortednets(elab::ElaboratedCircuit)
    shorted = falses(nnets(elab))   # a shorted pair's terminal is on the net
    other = falses(nnets(elab))     # so is some other terminal
    for i in 1:ninstances(elab)
        def = instancedefinition(elab, i)
        terminals = instanceterminals(elab, i)
        # the terminals come in pairs for a two terminal component and for
        # each port of a scattering block; any other component is compiled
        # into a refusal, and its terminals keep their nets
        if length(terminals) == 2 ||
                def isa ScatteringParameters || def isa LinearizedScattering
            for k in 1:2:length(terminals)
                a, b = terminals[k], terminals[k+1]
                if a == b
                    shorted[a] = true
                else
                    other[a] = other[b] = true
                end
            end
        else
            for n in terminals
                other[n] = true
            end
        end
    end
    return [n != 1 && shorted[n] && !other[n] ? 1 : n for n in 1:nnets(elab)]
end

"""
    compile(elab::ElaboratedCircuit; sorting = :name)
    compile(circuit::Circuit; sorting = :name)
    compile(c::CompiledCircuit; sorting = :name)

Lower a circuit to a [`CompiledCircuit`](@ref). A [`Circuit`](@ref) is
elaborated first; a `CompiledCircuit` is returned unchanged, and asking it
for a node order other than the one it carries is an error.

Components appear in the table in elaboration order, with a matched port's
own termination emitted as a resistor entry directly after the port. A net
which only components shorted across it touch is the ground net (see
[`shortednets`](@ref)). Nodes are numbered by [`calcnodesorting`](@ref)
with ground first; the default `sorting = :name` sorts the net names as
strings, since hierarchical net names are not integers, and `:number`
sorts integer node names by value.

Only components the solvers support can be lowered: a
[`GaussianChannel`](@ref), a [`VoltageSource`](@ref), or a component with
other than two terminals which is not a scattering block throws a
[`ComponentNotSupportedError`](@ref) naming the instance. A circuit with no connection to [`Ground`](@ref) throws an
`ArgumentError`.
"""
function compile(elab::ElaboratedCircuit; sorting::Symbol = :name)

    N = ninstances(elab)
    componentnames = String[]
    componenttypes = Symbol[]
    componentvalues = Any[]
    nodeindexvector = Int[]
    componenttemperatures = Dict{Int,Float64}()
    junctioncprs = Dict{Int,PolynomialCPR{Float64}}()
    sizehint!(componentnames, N)
    sizehint!(componenttypes, N)
    sizehint!(componentvalues, N)
    sizehint!(nodeindexvector, 2*N)

    # the nets are numbered in the order their terminals are met, by the
    # dense net ids the elaboration resolved, a net which only components
    # shorted across it touch being the ground net
    compiledas = shortednets(elab)
    netnumber = zeros(Int, nnets(elab))
    uniquenodevector = String[]
    for net in elab.terminalnets
        net = compiledas[net]
        if netnumber[net] == 0
            push!(uniquenodevector, elab.netnames[net])
            netnumber[net] = length(uniquenodevector)
        end
    end
    for net in eachindex(compiledas)
        netnumber[net] = netnumber[compiledas[net]]
    end

    ports = CompiledPort[]
    scatteringblocks = CompiledScatteringBlock[]
    namedenvironments = Pair{Int,String}[]
    hascouplings = !isempty(elab.couplings)
    compiledindex = zeros(Int, hascouplings ? N : 0)

    for i in 1:N
        def = instancedefinition(elab, i)
        path = elab.instancepaths[i]

        # A scattering block has two terminals per port and is not a two
        # terminal component, so it gets no table entry: it is compiled as
        # one `CompiledScatteringBlock` holding the nodes of every port; a
        # pumped block is one too, its mode coupling the stamp's affair.
        if def isa ScatteringParameters || def isa LinearizedScattering
            terminals = instanceterminals(elab, i)
            n = def.nports
            signalnodes = Vector{Int}(undef, n)
            refnodes = Vector{Int}(undef, n)
            for p in 1:n
                signalnodes[p] = netnumber[terminals[2*p-1]]
                refnodes[p] = netnumber[terminals[2*p]]
            end
            push!(scatteringblocks, CompiledScatteringBlock(def, signalnodes,
                refnodes, path))
            continue
        end

        # The common components are lowered by an `isa` chain, ordered by
        # how many of each a large circuit typically holds, and everything
        # else falls through to `lowercomponent`; on a heterogeneous vector
        # the chain of branches is one dispatch where a method per model
        # would be one per component
        typesymbol, value = if def isa Capacitor
            (:C, def.C)
        elseif def isa NonlinearInductor
            # a junction whatever its relation: the branch it makes, the
            # matrices it enters and the small signal inductance `L0` are
            # the same, and only the pointwise relation differs, which is
            # recorded below
            (:Lj, def.L0)
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
        hascouplings && (compiledindex[i] = marker)
        # a component which states its own temperature records it; the rest
        # take the temperature the analysis is run at
        t = componenttemperature(def)
        isnothing(t) || (componenttemperatures[marker] = t)
        # a junction whose relation is not the sinusoidal one records it;
        # the entry is otherwise a junction like any other
        if def isa NonlinearInductor
            cpr = junctioncpr(def, path)
            isnothing(cpr) || (junctioncprs[marker] = cpr)
        end

        if typesymbol == :K
            push!(nodeindexvector, 0)
            push!(nodeindexvector, 0)
            continue
        end

        terminals = instanceterminals(elab, i)
        if length(terminals) != 2
            throw(ComponentNotSupportedError(lazy"the component $(typeof(def)) at $(path) has $(length(terminals)) terminals; the solver supports two terminal components."))
        end
        n1 = netnumber[terminals[1]]
        n2 = netnumber[terminals[2]]
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
            elseif !isnothing(namedtermination(def.termination))
                # the termination names a resistor which already exists in
                # the table (the resistor a tuple netlist places across a
                # port, at the top level); the index is looked up once the
                # table is complete
                push!(namedenvironments, length(ports) + 1 =>
                    string(namedtermination(def.termination)))
            end
            push!(ports, CompiledPort(def.number, n1, n2, environment,
                marker))
        end
    end

    nodenames, nodeindices, renumber = sortnodes(uniquenodevector,
        nodeindexvector, calcnodesorting(uniquenodevector; sorting = sorting))

    componentnamedict = Dict{String,Int}()
    sizehint!(componentnamedict, length(componentnames))
    for (i, name) in enumerate(componentnames)
        # the table is keyed by name: a name written twice (or a Symbol
        # and a String which print alike) would be one component here and
        # two in the netlist
        haskey(componentnamedict, name) && throw(ArgumentError(
            lazy"The component name $(name) appears more than once in the circuit; component names must be unique."))
        componentnamedict[name] = i
    end

    # resolve the named terminations now that every name is in the table
    for (k, name) in namedenvironments
        i = get(componentnamedict, name, 0)
        if iszero(i) || componenttypes[i] !== :R
            throw(ArgumentError(lazy"The port $(componentnames[ports[k].component]) names $(name) as its termination, which is not a resistor in this circuit."))
        end
        ports[k] = CompiledPort(ports[k].number, ports[k].positivenode,
            ports[k].negativenode, i, ports[k].component)
    end

    # `sortnodes` renumbered the nodes. The port nodes recorded above are in
    # the pre-sort numbering and are re-read from the sorted table.
    ports = [CompiledPort(p.number, nodeindices[1, p.component],
        nodeindices[2, p.component], p.environment, p.component)
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

    couplings = sort!([(compiledindex[k], compiledindex[i], compiledindex[j])
        for (k, i, j) in elab.couplings]; by = first)
    group(t) = [i for (i, s) in enumerate(componenttypes) if s === t]

    Nnodes = length(uniquenodevector)
    return CompiledCircuit(nodenames, nodeindices, Nnodes,
        componentnames, componenttypes, tightenvalues(componentvalues),
        componentnamedict, componenttemperatures, junctioncprs,
        group(:C), group(:R), group(:L), group(:Lj), group(:I),
        group(:K), ports, scatteringblocks, couplings,
        circuittopology(componenttypes, nodeindices, Nnodes))
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
    any(p -> !iszero(p.environment), ports) || return nothing
    resistors = Dict{Tuple{Int,Int},Vector{Int}}()
    for i in eachindex(componenttypes)
        componenttypes[i] === :R || continue
        componentvalues[i] isa Number || continue
        key = minmax(nodeindices[1,i], nodeindices[2,i])
        push!(get!(() -> Int[], resistors, key), i)
    end
    for p in ports
        iszero(p.environment) && continue
        z = componentvalues[p.environment]
        z isa Number || continue
        key = minmax(p.positivenode, p.negativenode)
        candidates = get(resistors, key, nothing)
        isnothing(candidates) && continue
        for i in candidates
            i == p.environment && continue
            v = componentvalues[i]
            # `Number` admits a symbolic value, whose equality is another
            # symbolic value rather than a Bool, so the comparison is asked
            # for a definite `true`: a warning about two loads of the same
            # value cannot be made about values which are not yet numbers,
            # and a circuit carrying one is compiled rather than refused
            ((v isa Number && z isa Number) && (v == z) === true) || continue
            @warn "This port owns a matched environment of its own and a device resistor of the same value sits across the same terminals, so the port is loaded twice. If the resistor was written as the port's termination, either delete it or write the port as `termination = nothing` to keep it as the only load. If two loads are intended, this is correct and the warning can be ignored." port=p.number resistor=componentnames[i] value=z
        end
    end
    return nothing
end

compile(circuit::Circuit; sorting::Symbol = :name) =
    compile(elaborate(circuit); sorting = sorting)

# a compiled circuit is returned as it is, so the solver entry points
# accept one; it carries the node order it was compiled with, and asking
# for another one here would be silently ignored
function compile(c::CompiledCircuit; sorting::Symbol = :name)
    sorting === :name || throw(ArgumentError(
        lazy"the nodes of a compiled circuit are already ordered; pass the Circuit to compile it with sorting = $(repr(sorting))."))
    return c
end

# what the entry points of the analyses take as a circuit: anything
# `compile` accepts
const CompilableCircuit = Union{Circuit,ElaboratedCircuit,CompiledCircuit}

# === port and noise roles, read from the compiled circuit ===
#
# A compiled port states its own reference impedance and, through
# `environment`, which table entry realizes it. The functions below read
# those records; none of them looks for a resistor on a port's branch.

"""
    scatteringblockindex(c::CompiledCircuit, name)

The position in `c.scatteringblocks` of the block whose instance path is
`name`, a string, a symbol or an integer instance id, or zero when there
is none. A block is also named by the
`"<path>/port1"` spelling its stamp carries, which is what the solver
messages and the sensitivity labels print.
"""
function scatteringblockindex(c::CompiledCircuit, name)
    s = string(name)
    bare = endswith(s, "/port1") ? chop(s; tail = 6) : s
    for (k, b) in enumerate(c.scatteringblocks)
        (b.path == s || b.path == bare) && return k
    end
    return 0
end

"""
    coupledinductornames(c::CompiledCircuit)

The names of the two inductors each mutual inductor couples, two per `:K`
entry in the order the couplings are listed. Read from the `couplings` the
compilation resolved, which is where the pairing lives; the names are for
printing and for the output objects.
"""
coupledinductornames(c::CompiledCircuit) =
    String[c.componentnames[i] for (_, i1, i2) in c.couplings for i in (i1, i2)]

"""
    componentindex(c::CompiledCircuit, name)

The flat table index of the component named `name`, a `Symbol` or a
`String`. A name which is not in the circuit is an `ArgumentError` naming
it. This is the one place a name becomes an index, so every entry point
refuses an unknown one the same way.
"""
function componentindex(c::CompiledCircuit, name)
    idx = get(c.componentnamedict, string(name), 0)
    iszero(idx) && throw(ArgumentError(
        lazy"The component $(name) is not in this circuit."))
    return idx
end

"""
    orderedports(c::CompiledCircuit)

The ports of a compiled circuit ordered by port number. Throws an
`ArgumentError` for duplicate port numbers, a port with both terminals on
one node, or two ports on the same branch. Every port list the assembly
reads is built from this one, so the indices, the numbers, the
environments and the reference impedances are in one order.
"""
function orderedports(c::CompiledCircuit)
    name(p) = c.componentnames[p.component]
    numbers = [p.number for p in c.ports]
    if !allunique(numbers)
        n = first(k for k in numbers if count(==(k), numbers) > 1)
        names = join([name(p) for p in c.ports if p.number == n], ", ")
        throw(ArgumentError(lazy"The port number $(n) is given to more than one port: $(names)."))
    end
    pairs = Dict{Tuple{Int,Int},Int}()
    for (k, p) in enumerate(c.ports)
        p.positivenode == p.negativenode && throw(ArgumentError(
            lazy"The port $(name(p)) has both terminals on one node."))
        pair = minmax(p.positivenode, p.negativenode)
        other = get(pairs, pair, 0)
        iszero(other) || throw(ArgumentError(lazy"The ports $(name(c.ports[other])) and $(name(p)) are across the same pair of nodes; only one port is allowed between two nodes."))
        pairs[pair] = k
    end
    sp = sortperm(numbers)
    return c.ports[sp]
end

"""
    portreferenceimpedances(ports::Vector{CompiledPort}, values)

The reference impedance of each port of [`orderedports`](@ref), read from
a bound flat value table.

A port's `environment` field carries the other half of its role: the flat
index of the termination it owns, or zero for a port which owns none. That
index says which entry realizes the impedance, which the noise
classification needs (a port termination is an external bath, not an
internal noise channel) and the sensitivities need (perturbing that entry
also moves the wave normalization).

This is the impedance the incoming and outgoing waves are normalized to. It
is the port's declared `Z0`, which is the value of the port's own entry in
the table, and it is read from there for every port, whatever the port owns:
a symbolic or swept impedance then resolves the same way a component value
does, and an unterminated port resolves the same way a matched one. The two
cannot disagree with an environment either, because a matched environment
is generated with the port's own `Z0`, and a port whose termination names a
resistor has that resistor's value as its `Z0`.

A bound value must be a finite positive real number. The constructor lets a
symbol or a deferred value through so that it can be bound; this is where
what it bound to is checked, before any matrix or wave is built from it.
"""
function portreferenceimpedances(ports::Vector{CompiledPort}, values)
    z = [values[p.component] for p in ports]
    for (k, p) in enumerate(ports)
        checkportimpedance(z[k], p.number)
    end
    return z
end

"""
    noiseindices(c::CompiledCircuit, values,
        candidates = noisecandidates(c))

The flat table indices of the internal dissipative components, which are
the noise channels of the linearized analysis: every resistor which is not
a port's own termination, and every capacitor or inductor whose resolved
value in `values` has a nonzero imaginary part.

A port termination is an external bath rather than an internal channel and
is excluded by its role; any other resistor across a port's nodes is an
ordinary device resistor and is included. `candidates` are the components
examined, those which can be noise channels whatever their values, which a
plan computes once and passes in.
"""
function noiseindices(c::CompiledCircuit, values, candidates = noisecandidates(c))
    return [i for i in candidates if c.componenttypes[i] === :R ||
        (values[i] isa Complex && !iszero(values[i].im))]
end

"""
    isolatedsubnetworks(c::CompiledCircuit)

The subnetworks which no element connects to ground: the connected
components of the nodes joined by every capacitor, resistor, inductor,
junction and scattering block port, whatever their values, leaving out
the one which holds ground. Each is a sorted vector of node indices.
"""
function isolatedsubnetworks(c::CompiledCircuit)
    edges = Tuple{Int,Int}[]
    for k in eachindex(c.componenttypes)
        c.componenttypes[k] in (:C, :R, :L, :Lj) || continue
        push!(edges, (c.nodeindices[1, k], c.nodeindices[2, k]))
    end
    for cb in c.scatteringblocks, q in eachindex(cb.signalnodes)
        push!(edges, (cb.signalnodes[q], cb.refnodes[q]))
    end
    return nodecomponents(c.Nnodes, edges)
end

# Harmonic balance writes one equation per node and mode, the current law,
# and the currents of an isolated subnetwork sum to zero whatever its
# potential, so its equations leave that potential free at every mode and
# the system is singular. The solvers refuse such a circuit by name; the
# transient fixes the potential by a gauge row (see
# `transientfloatingcomponents`).
function checkisolatedsubnetworks(c::CompiledCircuit)
    islands = isolatedsubnetworks(c)
    isempty(islands) && return nothing
    names = join(["(" * join(c.nodenames[island], ", ") * ")"
        for island in islands], ", ")
    throw(ArgumentError(lazy"the nodes $names form a subnetwork which no element connects to ground, whose potential harmonic balance cannot determine; connect it to ground (a resistor or a capacitor will do), or solve the circuit in time with `transientsolve`, which fixes that potential."))
end

# the components which can be noise channels whatever their values: the
# resistors which are not a port's own environment, and the capacitors
# and inductors, which are when their value has an imaginary part; a plan
# holds them so that an assembly reads the values of these alone. A
# component whose terminals are one node carries no current and is none.
function noisecandidates(c::CompiledCircuit)
    owned = Set(p.environment for p in c.ports if !iszero(p.environment))
    return [i for (i, t) in enumerate(c.componenttypes)
        if ((t === :R && !(i in owned)) || t === :C || t === :L) &&
            c.nodeindices[1, i] != c.nodeindices[2, i]]
end
