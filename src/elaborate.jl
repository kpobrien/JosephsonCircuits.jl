
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
    autopaths::Vector{Tuple{Int,Int,String,Int}}   # (depth, seq, path, wire)
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
        Tuple{Int,Int,String,Int}[], Tuple{Int,Int,String,Int}[],
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

    # wires of each local instance's terminals (for subcircuits: pin wires)
    instwires = Vector{Vector{Int}}(undef, length(table.ids))
    # flattened global index of each local primitive instance (0 for subcircuits)
    localglobal = zeros(Int, length(table.ids))

    for (i, def) in enumerate(table.defs)
        id = table.ids[i]
        if def isa Circuit
            instwires[i] = flattencircuit!(st, def, joinpath_(path, id),
                depth + 1)
        else
            n = nterminals(def)
            w = Vector{Int}(undef, n)
            for t in 1:n
                w[t] = newwire!(st.wires)
                push!(st.terminalwires, w[t])
                push!(st.autopaths, (depth, nextseq!(st), path, w[t]))
            end
            instwires[i] = w
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
    nnets = 1
    for (i, w) in enumerate(st.terminalwires)
        r = findwire(st.wires, w)
        n = get(netofroot, r, 0)
        if n == 0
            nnets += 1
            n = nnets
            netofroot[r] = n
        end
        terminalnets[i] = n
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
    autopath = Dict{Int,Tuple{Int,Int,String}}()
    for (depth, seq, path, w) in st.autopaths
        n = get(netofroot, findwire(st.wires, w), 0)
        (n == 0 || n == 1) && continue
        best = get(autopath, n, (typemax(Int), typemax(Int), ""))
        if (depth, seq) < (best[1], best[2])
            autopath[n] = (depth, seq, path)
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
