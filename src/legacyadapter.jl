
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

const legacyallowedcomponents = ["Lj","NL","L","C","K","I","R","P"]

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
        return Port{typeof(value)}(legacyportnumber(name, value), 50.0, value)
    elseif typesymbol == :K
        return MutualInductor(value, string(node1), string(node2))
    else
        throw(ArgumentError(lazy"Unknown legacy component type $(typesymbol) for $(name)."))
    end
end

function legacyportnumber(name, value)
    v = value isa Symbolics.Num ? Symbolics.value(value) : value
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
from the name prefix exactly as in [`parsecircuit`](@ref): Lj, NL, L, C, K,
I, R, and P. Component name prefixes are interpreted only inside this
adapter; typed component models never infer behavior from instance names.

Node labels become net names, so lowering the result back with
[`parsesortcircuit`](@ref) reproduces the legacy parse exactly. When
`circuitdefs` is supplied, symbolic values are substituted with
[`valuetonumber`](@ref) during conversion; otherwise values pass through
verbatim and `circuitdefs` may be supplied to the analysis as usual.

# Examples
```jldoctest
julia> Circuit([("P1","1","0",1),("R1","1","0",50.0),("C1","1","0",1e-12)]) isa Circuit
true
```
"""
function Circuit(netlist::AbstractVector{<:Tuple})
    return legacycircuit(netlist, nothing)
end

function Circuit(netlist::AbstractVector{<:Tuple}, circuitdefs::AbstractDict)
    return legacycircuit(netlist, circuitdefs)
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

# === lowering: ElaboratedCircuit -> ParsedSortedCircuit ===

function lowercomponent(def::Inductor, path) 
    return :L, def.L
end
lowercomponent(def::Capacitor, path) = :C, def.C
lowercomponent(def::Resistor, path) = :R, def.R
lowercomponent(def::CurrentSource, path) = :I, def.I
lowercomponent(def::Port, path) = :P, def.value
lowercomponent(def::MutualInductor, path) = :K, def.K
lowercomponent(def::LegacyNL, path) = :NL, def.value
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

"""
    parsesortcircuit(elab::ElaboratedCircuit; sorting = :name)
    parsesortcircuit(circuit::Circuit; sorting = :name)

Lower an elaborated typed circuit to the [`ParsedSortedCircuit`](@ref)
consumed by the existing solvers. Hierarchical instance paths become
component names and net names become node names, so results keep the
keyed array behavior with hierarchical keys such as "cell37/net2"; circuits
constructed by the legacy tuple adapter round trip exactly.

Only components supported by the current solvers can be lowered. Each port
of a passive [`ScatteringBlock`](@ref) lowers to one `:S` component whose
two nodes are the signal and reference terminals of the port and whose
value is a [`ScatteringStamp`](@ref) sharing the block definition; the
solvers stamp the multiport admittance of the block at every mode
frequency. A [`GaussianChannel`](@ref), a non-sinusoidal
[`NonlinearInductor`](@ref), or a `ScatteringBlock` with an arbitrary noise
covariance raises a [`ComponentNotSupportedError`](@ref) naming the
instance. The default `sorting` is `:name` because automatic
hierarchical net names are not integers; see [`sortnodes`](@ref).
"""
function parsesortcircuit(elab::ElaboratedCircuit; sorting::Symbol = :name)
    N = ninstances(elab)
    componentnames = String[]
    componenttypes = Symbol[]
    componentvalues = Any[]
    nodeindexvector = Int[]
    sizehint!(componentnames, N)
    sizehint!(componenttypes, N)
    sizehint!(componentvalues, N)
    sizehint!(nodeindexvector, 2*N)
    mutualinductorbranchnames = String[]
    uniquenodedict = Dict{String,Int}()
    uniquenodevector = String[]

    componenttemperatures = Dict{Int,Float64}()
    couplingnames = Dict{Int,Tuple{String,String}}()
    for (k, i1, i2) in elab.couplings
        couplingnames[k] = (elab.instancepaths[i1], elab.instancepaths[i2])
    end

    for i in 1:N
        def = instancedefinition(elab, i)
        path = elab.instancepaths[i]
        if def isa ScatteringBlock
            # each port lowers to one :S component whose two nodes are the
            # signal and reference terminals of the port, carrying the
            # shared block definition and the port index; the solvers
            # reassemble the multiport admittance coupling from the shared
            # identity (see ScatteringStampSystem).
            if def.noise isa NoiseCovariance
                throw(ComponentNotSupportedError(lazy"the ScatteringBlock at $(path) has an arbitrary noise covariance, which permits active blocks; the solver supports passive scattering blocks (noise = Passive() or ThermalEquilibrium). It parsed, validated, and elaborated successfully."))
            end
            terminals = instanceterminals(elab, i)
            for p in 1:def.nports
                push!(componentnames, path * "/port" * string(p))
                push!(componenttypes, :S)
                push!(componentvalues, ScatteringStamp(def, p))
                if def.noise isa ThermalEquilibrium
                    componenttemperatures[length(componentnames)] =
                        Float64(def.noise.temperature)
                end
                push!(nodeindexvector, processnode(uniquenodedict,
                    uniquenodevector, elab.netnames[terminals[2*p-1]]))
                push!(nodeindexvector, processnode(uniquenodedict,
                    uniquenodevector, elab.netnames[terminals[2*p]]))
            end
            continue
        end
        typesymbol, value = lowercomponent(def, path)
        push!(componentnames, path)
        push!(componenttypes, typesymbol)
        push!(componentvalues, value)
        # a component which states its temperature keeps it; the rest take
        # the one the analysis is run at
        t = componenttemperature(def)
        isnothing(t) || (componenttemperatures[length(componentnames)] = t)
        if typesymbol == :K
            l1, l2 = couplingnames[i]
            push!(mutualinductorbranchnames, l1)
            push!(mutualinductorbranchnames, l2)
            push!(nodeindexvector, 0)
            push!(nodeindexvector, 0)
        else
            terminals = instanceterminals(elab, i)
            if length(terminals) != 2
                throw(ComponentNotSupportedError(lazy"the component $(typeof(def)) at $(path) has $(length(terminals)) terminals; the solver supports two terminal components."))
            end
            for t in 1:2
                push!(nodeindexvector, processnode(uniquenodedict,
                    uniquenodevector, elab.netnames[terminals[t]]))
            end
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

    return ParsedSortedCircuit(nodeindices, nodenames,
        mutualinductorbranchnames, componentnames, componenttypes,
        tightenvalues(componentvalues), componentnamedict,
        length(uniquenodevector), componenttemperatures)
end

# the temperature a component states, or `nothing`. Only the components which
# can dissipate carry one; a scattering block states its temperature through
# its noise model instead (see [`ThermalEquilibrium`](@ref)).
componenttemperature(def::Resistor) = def.temperature
componenttemperature(def::Capacitor) = def.temperature
componenttemperature(def::Inductor) = def.temperature
componenttemperature(def) = nothing

function parsesortcircuit(circuit::Circuit; sorting::Symbol = :name)
    return parsesortcircuit(elaborate(circuit); sorting = sorting)
end

# Identity passthrough so that the existing solver entry points, which call
# parsesortcircuit on their circuit argument, accept a pre-parsed circuit.
function parsesortcircuit(psc::ParsedSortedCircuit; sorting::Symbol = :name)
    return psc
end

# narrow a Vector{Any} to a concrete element type when possible, matching
# the value vector types produced by parsecircuit from typed tuple vectors.
function tightenvalues(values::Vector{Any})
    return map(identity, values)
end

# === solver entry points for the typed representation ===

"""
    hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        circuit::Circuit, circuitdefs = Dict{Symbol,Number}();
        sorting = :name, keyword arguments...)

Harmonic balance solution of a typed [`Circuit`](@ref). The circuit is
elaborated and lowered with [`parsesortcircuit`](@ref) and solved with the
existing solver; all keyword arguments of the legacy method are supported.
`circuitdefs` is only needed when component values are symbolic. The
default `sorting` is `:name` because hierarchical net names are not
integers.
"""
function hbsolve(ws, wp::NTuple{N,Number}, sources::Vector,
        Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
        circuit::Circuit, circuitdefs::AbstractDict = Dict{Symbol,Number}();
        sorting::Symbol = :name, kwargs...) where {N,M}
    return hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        parsesortcircuit(circuit; sorting = sorting), circuitdefs; kwargs...)
end

"""
    hbnlsolve(w, Nharmonics, sources, circuit::Circuit,
        circuitdefs = Dict{Symbol,Number}(); sorting = :name,
        keyword arguments...)

Nonlinear harmonic balance solution of a typed [`Circuit`](@ref); see
[`hbsolve`](@ref).
"""
function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        circuit::Circuit, circuitdefs::AbstractDict = Dict{Symbol,Number}();
        sorting::Symbol = :name, kwargs...) where N
    return hbnlsolve(w, Nharmonics, sources,
        parsesortcircuit(circuit; sorting = sorting), circuitdefs; kwargs...)
end

"""
    hblinsolve(w, circuit::Circuit, circuitdefs = Dict{Symbol,Number}();
        sorting = :name, keyword arguments...)

Linearized harmonic balance solution of a typed [`Circuit`](@ref); see
[`hbsolve`](@ref).
"""
function hblinsolve(w, circuit::Circuit,
        circuitdefs::AbstractDict = Dict{Symbol,Number}();
        sorting::Symbol = :name, kwargs...)
    return hblinsolve(w, parsesortcircuit(circuit; sorting = sorting),
        circuitdefs; kwargs...)
end
