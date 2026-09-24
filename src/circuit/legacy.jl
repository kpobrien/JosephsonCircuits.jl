# The deprecated input formats and entry points.
#
# The typed `Circuit` whose frequency dependent values are closures of the
# frequency is the input format of the package. Two older ways of writing
# a circuit are deprecated, and each is converted here, with a warning,
# into one written the current way, which then follows the same path. The
# deprecated entry points, at the end of the file, warn and forward to
# their replacements. Everything specific to them lives in this file, so
# that they are removed by deleting it and test/circuit/legacy.jl.
#
# A netlist of `(name, node1, node2, value)` tuples, the original input
# format: the name prefix table with the two functions that read it, the
# convention that a port's reference impedance is the resistor placed
# across it and the port termination which records that resistor, the
# tuple forms of the entry points which took a tuple netlist in v0.5.4,
# which sort the nodes by number as the format always did, and the tuple
# netlist file reader and writer.
#
# A value written as an expression in a parameter named by the
# `symfreqvar` keyword of the solvers, rewritten as a
# `FrequencyDependent` closure.
#
# The keywords of the nonlinear solvers which v0.5.4 had and the solvers
# no longer read: `ftol`, the absolute residual tolerance under its old
# name, the line search settings `switchofflinesearchtol` and `alphamin`,
# and `maxharmonics` and `maxpumpharmonics`, whose role the retained
# harmonics took.
#
# Outside this file, the `Circuit(netlist)` constructor in
# circuit/parse.jl hands a netlist whose entries end in values rather than
# components to `legacycircuit`, the Symbolics extension unwraps a `Num`
# port number for `legacyportnumber`, the three solver entry points
# which hold a compiled circuit beside its definitions hand a
# `symfreqvar` to `frequencydependentcircuit`, and `hbnlsolve` and
# `hbsolve` hand their deprecated keywords to `deprecatedsolverkeywords`.

# Unwrap a wrapped symbolic value to whatever it holds. The Symbolics
# extension adds the method for `Num`; everything else is already unwrapped.
unwrapvalue(value) = value

"""
    LegacyTermination(component)

The port owned environment of a tuple netlist: a resistor the netlist
already contains, named by its instance identifier.

Internal. A tuple netlist states a port's impedance by placing a resistor
across it and carries no role marker, so the adapter finds that resistor once
and records which one it is. Everything downstream then reads the port's
environment from the port, exactly as for a native matched port, and nothing
searches for a resistor sharing a port's branch.
"""
struct LegacyTermination{I} <: AbstractPortTermination
    component::I
end
namedtermination(t::LegacyTermination) = t.component
showtermination(io::IO, t::LegacyTermination) =
    print(io, ", termination = LegacyTermination(", repr(t.component), ")")

# === the tuple netlist -> Circuit ===

# The component type prefixes of the tuple format. Two letter prefixes must
# come before one letter prefixes with the same first letter; see
# `checkcomponenttypes`.
const legacyallowedcomponents = ["Lj","L","C","K","I","R","P"]

"""
    parsecomponenttype(name::String,allowedcomponents::Vector{String})

The index in `allowedcomponents` of the one or two letter prefix which
matches the start of the component name `name`. Prefixes are tried in
order and the first match wins, so a two letter prefix listed after a one
letter prefix with the same first letter can never match;
[`checkcomponenttypes`](@ref) detects that ordering mistake.

# Examples
```jldoctest
julia> JosephsonCircuits.parsecomponenttype("L10",["Lj","L","C","K","I","R","P"])
2

julia> [JosephsonCircuits.parsecomponenttype(c,["Lj","L","C","K","I","R","P"]) for c in ["Lj","L","C","K","I","R","P"]]
7-element Vector{Int64}:
 1
 2
 3
 4
 5
 6
 7
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
    throw(ArgumentError(lazy"No component in allowedcomponents matches the name $(name)."))
end

"""
    checkcomponenttypes(allowedcomponents::Vector{String})

Check that [`parsecomponenttype`](@ref) maps each prefix in
`allowedcomponents` back to its own index, and throw an `ArgumentError`
otherwise. This fails when a two letter prefix is listed after a one letter
prefix with the same first letter, which would shadow it.

# Examples
```jldoctest
julia> JosephsonCircuits.checkcomponenttypes(["Lj","L","C","K","I","R","P"])
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
        return NonlinearInductor(value, sin)
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
    Circuit(netlist::AbstractVector, circuitdefs::AbstractDict)

Construct a typed [`Circuit`](@ref) from a tuple netlist, which is
deprecated: this warns and will be removed in a future release. Each
entry is `(name, node1, node2, value)` and the component type is taken
from the prefix of `name`: `Lj` (Josephson junction), `L`, `C`, `K`
(mutual inductor, whose "nodes" are the two inductor names), `I`, `R`, and
`P` (port, whose value is the port number). Only this adapter reads a name
prefix; typed component models never infer behavior from an instance name.
The one argument `Circuit(netlist)` reads a tuple netlist the same way
when its entries end in values rather than typed components.

Node labels become net names, so [`compile`](@ref) of the result gives the
same tables the tuple netlist always produced. A port's reference
impedance is the value of the single resistor placed across it, which the
adapter records as the port's [`LegacyTermination`](@ref). When
`circuitdefs` is given, values are resolved with [`valuetonumber`](@ref)
during conversion; otherwise they pass through unchanged and
`circuitdefs` is given to the analysis as usual.

# Examples
```julia
julia> Circuit([("P1","1","0",1),("R1","1","0",50.0),("C1","1","0",1e-12)]) isa Circuit
true
```
"""
function Circuit(netlist::AbstractVector, circuitdefs::AbstractDict)
    # The element type is not restricted to `Tuple`: a netlist built by
    # pushing onto a `Vector{Any}` is accepted, and `legacycircuit` checks
    # each entry and reports what is wrong with it.
    if any(isnetlistentry, netlist)
        throw(ArgumentError("A netlist of typed components takes no circuitdefs here; parameterized values are resolved by the analysis, which takes the definitions as usual."))
    end
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
    # the resistors across each pair of nodes, by index into the component
    # table rather than by a copy of their names and values, which may be
    # symbolic
    resistorsat = Dict{Tuple{String,String},Vector{Int}}()
    for (i, entry) in enumerate(netlist)
        components[i].second isa Resistor || continue
        n1, n2 = string(entry[2]), string(entry[3])
        push!(get!(Vector{Int}, resistorsat, (n1, n2)), i)
        n1 == n2 || push!(get!(Vector{Int}, resistorsat, (n2, n1)), i)
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
            names = join((components[j].first for j in rs), ", ")
            throw(ArgumentError(lazy"Only one resistor allowed per port. Port $(p.number) has $(length(rs)) resistors across it ($(names)), and a legacy netlist has no way to say which one is its environment. Give the port a single resistor, or write the circuit in the typed format, where a port states its own termination and any number of device resistors may share its terminals."))
        end
        resistor = components[only(rs)]
        components[i] = Pair{String,Any}(components[i].first,
            Port(p.number; Z0 = resistor.second.R,
                termination = LegacyTermination(resistor.first)))
    end
    return components
end

const tuplenetlistmessage = "The netlist of (name, node1, node2, value) tuples is deprecated and will be removed in a future release. Write the circuit as a Circuit of typed components, `Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(100e-15)), ...])`, where a port states its own reference impedance and needs no resistor across it; see the Circuit docstring."

# Convert the tuple netlist to components and connection groups, with the
# deprecation warning attributed to the entry point `caller` the netlist
# was given to. Each distinct node label becomes one `Net` holding every
# terminal on it, in order of first appearance, with `Ground` appended to
# net "0".
function legacycircuit(netlist, circuitdefs; pins = nothing, ports = nothing,
        caller::Symbol = :Circuit)
    # a list of typed components handed over without the Circuit around it,
    # or components without their connections, is not a tuple netlist, and
    # is said to be so before a warning about one
    if any(e -> e isa Pair && last(e) isa AbstractComponent, netlist)
        throw(ArgumentError("The list holds name => component pairs, which are the components of a circuit: give their connections as well, as Circuit(components, connections)."))
    end
    if any(e -> e isa Tuple && last(e) isa AbstractComponent, netlist)
        throw(ArgumentError("The netlist holds typed components, which are passed as a circuit: write Circuit(netlist)."))
    end
    Base.depwarn(tuplenetlistmessage, caller; force = true)
    if !isnothing(pins) || !isnothing(ports)
        throw(ArgumentError("A tuple netlist has no interface; give pins and ports to a netlist of typed components or to the connection-group form."))
    end
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
        components[i] = Pair{String,Any}(String(name),
            legacycomponent(typesymbol, name, node1, node2, value))
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
    connections = Vector{Net{String,Vector{Any}}}(undef, length(nodeorder))
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

# === the symbolic frequency variable, deprecated ===
#
# The deprecated way of writing a frequency dependent value is an
# expression in a parameter named by the `symfreqvar` keyword of the
# solvers. `FrequencyDependent` states the same law as a closure of the
# frequency, which needs neither a free parameter nor a keyword, so a
# circuit given a symbolic frequency variable is rewritten here into one
# whose dependent values are closures, and then follows the same path as
# one written that way.

const symfreqvarmessage = "The `symfreqvar` keyword is deprecated and will be removed in a future release. Write a frequency dependent value as a closure of the frequency, `Capacitor(FrequencyDependent(w -> C0*(1 + im*w/wc)))`, which needs neither a symbolic variable nor a keyword; see the FrequencyDependent docstring."

"""
    frequencydependentcircuit(psc, circuitdefs, symfreqvar, caller)

A compiled circuit whose values depending on the symbolic frequency
variable `symfreqvar` are rewritten as [`FrequencyDependent`](@ref)
closures resolving them at `circuitdefs` and at the frequency, and a
deprecation warning attributed to `caller`.

Every other value is left as it was, including one which still depends on
a parameter the definitions do not give, so that the undefined parameter
is still reported by [`checkcomponentvaluesdefined`](@ref) naming its
component rather than failing later inside a closure.
"""
function frequencydependentcircuit(psc::CompiledCircuit, circuitdefs,
        symfreqvar, caller::Symbol)
    Base.depwarn(symfreqvarmessage, caller; force = true)
    # the definitions resolved first, once for the table, so that only the
    # frequency is left for a closure to give, and a value named by a symbol
    # whose definition is written in the variable is found
    partials = componentvaluestonumber(psc.componentvalues,
        definitiontable(circuitdefs))
    values = Any[frequencydependentvalue(v, partial, symfreqvar)
        for (v, partial) in zip(psc.componentvalues, partials)]
    return CompiledCircuit(psc.nodenames, psc.nodeindices, psc.Nnodes,
        psc.componentnames, psc.componenttypes, values,
        psc.componentnamedict, psc.componenttemperatures, psc.junctioncprs,
        psc.capacitors, psc.resistors,
        psc.inductors, psc.junctions, psc.currentsources,
        psc.mutualinductors, psc.ports, psc.scatteringblocks,
        psc.couplings, psc.topology)
end

# the value `value`, resolved at the definitions to `partial`, as a closure
# of the frequency when `partial` depends on the variable and on nothing
# else
function frequencydependentvalue(value, partial, symfreqvar)
    checkissymbolic(partial) || return value
    variables = circuitvariables(partial)
    any(v -> isequal(v, symfreqvar), variables) || return value
    all(v -> isequal(v, symfreqvar), variables) || return value
    # the frequency reaches the value two ways: as the symbolic variable it
    # substitutes, and as the argument of any frequency dependent leaf the
    # same expression already carries, which `substitutefreq` evaluates.
    return FrequencyDependent(
        w -> substitutefreq(
            valuetonumber(partial, Dict{Any,Any}(symfreqvar => w)), w))
end

# === the deprecated keywords of the nonlinear solvers ===

"""
    deprecatedsolverkeywords(caller::Symbol, atol; ftol = nothing,
        switchofflinesearchtol = nothing, alphamin = nothing,
        maxharmonics = nothing, maxpumpharmonics = nothing)

The absolute residual tolerance of a solve given the deprecated keywords
of the nonlinear solvers: `ftol` when it is given, which is `atol` under
its old name, and `atol` otherwise. Each deprecated keyword given warns,
attributed to `caller`; `switchofflinesearchtol` and `alphamin`, which the
line search no longer has, and `maxharmonics` of `hbnlsolve` and
`maxpumpharmonics` of `hbsolve`, whose role the retained harmonics took,
are ignored.
"""
function deprecatedsolverkeywords(caller::Symbol, atol; ftol = nothing,
        switchofflinesearchtol = nothing, alphamin = nothing,
        maxharmonics = nothing, maxpumpharmonics = nothing)
    if !isnothing(ftol)
        Base.depwarn("The `ftol` kwarg is deprecated: the absolute residual tolerance is `atol` in every solver of the package. Please use `atol` to avoid errors in future versions.", caller; force = true)
        atol = ftol
    end
    unused(name) = Base.depwarn(lazy"The `$(name)` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions.", caller; force = true)
    isnothing(switchofflinesearchtol) || unused(:switchofflinesearchtol)
    isnothing(alphamin) || unused(:alphamin)
    retained(name, by) = Base.depwarn(lazy"The `$(name)` kwarg is deprecated and no longer used. `$(by)` is the retained set of modes and `Nevaluationharmonics` the grid on which the nonlinearity is sampled. Please remove it to avoid errors in future versions.", caller; force = true)
    isnothing(maxharmonics) || retained(:maxharmonics, :Nharmonics)
    isnothing(maxpumpharmonics) || retained(:maxpumpharmonics, :Npumpharmonics)
    return atol
end

# === the tuple forms of the entry points ===
#
# Each entry point of v0.5.4 which took a tuple netlist converts it, which
# warns, compiles it with the nodes sorted by number as the tuple format
# always did, and forwards the compiled circuit to the typed form. The
# solvers and the matrix builders of the tuple format took the order as
# their keyword `sorting`, which is now `compile`'s, and still take it
# here.

# the netlist compiled the way the tuple format ordered its nodes
legacycompiled(netlist, caller::Symbol; sorting::Symbol = :number) =
    compile(legacycircuit(netlist, nothing; caller = caller); sorting = sorting)

function hbsolve(ws, wp::NTuple{N,Number}, sources,
        Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
        netlist::AbstractVector, circuitdefs::AbstractDict = Dict{Symbol,Any}();
        sorting::Symbol = :number, kwargs...) where {N,M}
    return hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        legacycompiled(netlist, :hbsolve; sorting = sorting), circuitdefs;
        kwargs...)
end

function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        netlist::AbstractVector, circuitdefs::AbstractDict = Dict{Symbol,Any}();
        sorting::Symbol = :number, kwargs...) where {N}
    return hbnlsolve(w, Nharmonics, sources,
        legacycompiled(netlist, :hbnlsolve; sorting = sorting), circuitdefs;
        kwargs...)
end

function hblinsolve(w, netlist::AbstractVector,
        circuitdefs::AbstractDict = Dict{Symbol,Any}();
        sorting::Symbol = :number, kwargs...)
    return hblinsolve(w, legacycompiled(netlist, :hblinsolve; sorting = sorting),
        circuitdefs; kwargs...)
end

function numericmatrices(netlist::AbstractVector, circuitdefs::AbstractDict;
        Nmodes::Int = 1, sorting::Symbol = :number)
    return numericmatrices(legacycompiled(netlist, :numericmatrices;
        sorting = sorting), circuitdefs; Nmodes = Nmodes)
end

function symbolicmatrices(netlist::AbstractVector; Nmodes::Int = 1,
        sorting::Symbol = :number)
    return symbolicmatrices(legacycompiled(netlist, :symbolicmatrices;
        sorting = sorting); Nmodes = Nmodes)
end

function exportnetlist(netlist::AbstractVector, circuitdefs::Dict;
        port::Int = 1, jj::Bool = true)
    return exportnetlist(legacycompiled(netlist, :exportnetlist), circuitdefs;
        port = port, jj = jj)
end

# === the tuple netlist file ===
#
# A text file with one `name node1 node2 value` entry per line, read into
# a tuple netlist and written from one.

"""
    export_netlist(filename, circuit, circuitdefs)

Export the netlist in `circuit` to the file with name and path `filename`.
"""
function export_netlist(filename, circuit, circuitdefs)
    open(filename, "w") do io
        export_netlist!(io, circuit, circuitdefs)
    end
    return nothing
end

"""
    export_netlist(filename, circuit)

Export the netlist in `circuit` to the file with name and path `filename`.
"""
function export_netlist(filename, circuit)
    return export_netlist(filename, circuit, Dict())
end

"""
    export_netlist!(io::IO, circuit, circuitdefs)

Export the netlist in `circuit` to the IOBuffer or IOStream `io`.

# Examples
```julia
julia> io = IOBuffer();JosephsonCircuits.export_netlist!(io, [("P","1","0",1),("R","1","0",50.0)],Dict());println(String(take!(io)))
P 1 0 1
R 1 0 50.0
```
"""
function export_netlist!(io::IO, circuit::AbstractVector, circuitdefs::Dict)
    Base.depwarn(tuplenetlistmessage, :export_netlist; force = true)
    for i in eachindex(circuit)
        c = circuit[i]
        for j in eachindex(c)
            cj = c[j]
            if j > 1
                write(io," ")
            end
            write(io,string(substitutedefs(cj,circuitdefs)))
        end
        write(io,"\n")
    end
end

# A field of a netlist line with the definitions substituted, for writing:
# a symbolic value resolved as far as the definitions go, and anything
# else (a name, a node label, a number) as it is written, since a name
# field is not a value to look up.
substitutedefs(value, circuitdefs) =
    checkissymbolic(value) ? valuetonumber(value, circuitdefs) : value

"""
    import_netlist(filename)

Import the netlist from the file with name and path `filename` and return
it as a vector of `(name, node1, node2, value)` tuples. The value field is
`Any`: a number for a literal and a `CircuitValue` for an expression.
"""
function import_netlist(filename)
    # the value field is a number when the netlist holds a literal and a
    # `CircuitValue` when it holds an expression, so the tuple is
    # heterogeneous. Pass your own vector to `import_netlist!` to pin a
    # narrower element type.
    circuit = Tuple{String,String,String,Any}[]
    open(filename, "r") do io
        import_netlist!(io, circuit)
    end
    return circuit
end

"""
    import_netlist!(io::IO, circuit)

Import the netlist from the IOBuffer or IOStream `io` to the vector of tuples
`circuit`.

# Examples
```julia
julia> io = IOBuffer();circuit1=[("P","1","0",1),("R","1","0",50.0)];JosephsonCircuits.export_netlist!(io,circuit1,Dict());circuit2 = Tuple{String,String,String,Any}[];JosephsonCircuits.import_netlist!(io,circuit2);circuit2
2-element Vector{Tuple{String, String, String, Any}}:
 ("P", "1", "0", 1.0)
 ("R", "1", "0", 50.0)
```
"""
function import_netlist!(io::IO, circuit::AbstractVector)
    Base.depwarn(tuplenetlistmessage, :import_netlist; force = true)
    seekstart(io)
    for line in eachline(io)
        split_line = split(strip(line),r"\s+")
        if length(split_line) != 4
            error(lazy"each line should have component name, node1, node2, component value")
        end
        value = parsecomponentvalue(split_line[4])
        push!(circuit,(split_line[1],split_line[2],split_line[3],value))
    end
    return nothing
end



# Reading a component value back out of a netlist line. Only
# `import_netlist!` above uses it.
"""
    parsecomponentvalue(s::AbstractString)

Parse a SPICE netlist component value into a number or a `CircuitValue`.
Nothing is evaluated into a module, so a netlist cannot introduce arbitrary
code.
"""
function parsecomponentvalue(s::AbstractString)
    v = tryparse(Float64, s)
    isnothing(v) || return v
    return CircuitValues.fromexpr(Meta.parse(s))
end

# === the deprecated entry points ===
#
# Each one warns once and forwards to its replacement.

# `connectS(Sa, k, l)` and `connectS(Sa, Sb, k, l)` were split into
# `intraconnectS` and `interconnectS` so that each can also take noise
# covariance matrices: with one name, a second matrix argument would be
# ambiguous between a second scattering matrix and a noise covariance.
function connectS(Sa::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}
    Base.depwarn(lazy"`connectS(Sa::AbstractArray, k::Int, l::Int)` is deprecated, use `intraconnectS(Sa, k, l)` instead.", :connectS; force=true)
    return intraconnectS(Sa, k, l; nbatches = nbatches)
end

function connectS(Sa::AbstractArray{T,N}, Sb::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}
    Base.depwarn(lazy"`connectS(Sa::AbstractArray, Sb::AbstractArray, k::Int, l::Int)` is deprecated, use `interconnectS(Sa, Sb, k, l)` instead.", :connectS; force=true)
    return interconnectS(Sa, Sb, k, l; nbatches = nbatches)
end

function connectS!(Sout, Sa, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())
    Base.depwarn(lazy"`connectS!(Sout, Sa, k::Int, l::Int)` is deprecated, use `intraconnectS!(Sout, Sa, k, l)` instead.", :connectS!; force=true)
    return intraconnectS!(Sout, Sa, k, l; nbatches = nbatches)
end

function connectS!(Sout, Sa, Sb, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())
    Base.depwarn(lazy"`connectS!(Sout, Sa, Sb, k::Int, l::Int)` is deprecated, use `interconnectS!(Sout, Sa, Sb, k, l)` instead.", :connectS!; force=true)
    return interconnectS!(Sout, Sa, Sb, k, l; nbatches = nbatches)
end

# `X_Y_to_sympletic_pair` and `X_Y_to_sympletic_block` were misspelled.
function X_Y_to_sympletic_pair(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    Base.depwarn("`X_Y_to_sympletic_pair` is deprecated, use `X_Y_to_symplectic_pair` instead.", :X_Y_to_sympletic_pair; force=true)
    return X_Y_to_symplectic_pair(X, Y)
end

function X_Y_to_sympletic_block(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real})
    Base.depwarn("`X_Y_to_sympletic_block` is deprecated, use `X_Y_to_symplectic_block` instead.", :X_Y_to_sympletic_block; force=true)
    return X_Y_to_symplectic_block(X, Y)
end


#     hbsolve(ws, wp, Ip, Nsignalmodes::Int, Npumpmodes::Int, circuit,
#         circuitdefs; pumpports = [1], keyword arguments...)
#
# The original `hbsolve` signature: a single pump at the scalar frequency
# `wp`, applied to `pumpports` with the currents `Ip`, four wave mixing
# only, and mode counts given as integers rather than tuples. It is
# translated into a call of the current solvers and warns that it is
# deprecated. A tuple netlist is compiled with its nodes sorted as
# `sorting` says, and any other circuit as `compile` compiles it; the line
# search keywords are handed to `hbnlsolve`, which warns that they are
# ignored. Note that the signal modes of the result are ordered as
# `hblinsolve` orders them (the signal at index 1, the rest as listed in
# `modes`), which need not match the order the original solver used.
function hbsolve(ws, wp, Ip, Nsignalmodes::Int, Npumpmodes::Int, circuit,
    circuitdefs; pumpports = [1], iterations = 1000, ftol = 1e-8,
    switchofflinesearchtol = nothing, alphamin = nothing,
    symfreqvar = nothing,
    nbatches = Base.Threads.nthreads(), sorting = :number,
    returnS::Bool = true, returnSnoise::Bool = false, returnQE::Bool = true,
    returnCM::Bool = true, returnnodeflux::Bool = false,
    returnvoltage::Bool = false, returnnodefluxadjoint::Bool = false,
    returnvoltageadjoint::Bool = false, keyedarrays::Bool = false,
    sensitivitynames::AbstractVector = String[],
    returnSsensitivity::Bool = false, returnZ = nothing,
    returnZadjoint = nothing, returnZsensitivity = nothing,
    returnZsensitivityadjoint = nothing,
    factorization = nothing)

    Base.depwarn("""
    This form of hbsolve, with a single pump frequency and integer harmonic
    counts, is deprecated: it calls the harmonic balance solvers hbnlsolve
    and hblinsolve, which take any number of pump tones and ports, with the
    syntax of the legacy solver, which supported four wave mixing of one
    strong tone only. Please switch to hbsolve(ws, (wp,), sources,
    (Nmodulationharmonics,), (Npumpharmonics,), circuit, circuitdefs).
        """, :hbsolve; force=true)

    # the single pump as a one element frequency tuple
    w = (wp,)
    Nharmonics = (2*Npumpmodes,)

    # one source per pump port, all at the pump frequency
    length(pumpports) == length(Ip) || throw(ArgumentError(
        lazy"there are $(length(pumpports)) pump ports and $(length(Ip)) pump currents; give one current per port."))
    sources = [(mode = (1,), port = pumpports[i], current = Ip[i])
        for i in eachindex(pumpports)]

    # the pump harmonics: odd harmonics only, which is four wave mixing
    freq, indices = pumpmodeset((wp,), Nharmonics, Nharmonics;
        dc = false, odd = true, even = false)

    Nmodes = length(freq.modes)

    psc = circuit isa AbstractVector ?
        legacycompiled(circuit, :hbsolve; sorting = sorting) : compile(circuit)
    # the deprecated symbolic frequency variable, in circuit/legacy.jl
    isnothing(symfreqvar) || (psc = frequencydependentcircuit(psc,
        circuitdefs, symfreqvar, :hbsolve))
    nm=numericmatrices(psc, circuitdefs, Nmodes = Nmodes)

    # `factorization = nothing` leaves the preconditioner to `Automatic`; a
    # factorization pins the block diagonal, the only preconditioner this
    # form builds with one
    nonlinear = hbnlsolve(w, sources, freq, indices, psc, nm;
        iterations = iterations, atol = ftol,
        switchofflinesearchtol = switchofflinesearchtol, alphamin = alphamin,
        keyedarrays = keyedarrays,
        sensitivitynames = sensitivitynames,
        method = NewtonKrylov(preconditioner = isnothing(factorization) ?
            Automatic() : BlockDiagonal(factorization = factorization)))

    # the signal modes: the signal and the even pump harmonics on either
    # side of it
    signalfreq =truncfreqs(
        calcfreqsdft((Nsignalmodes,)),
        dc=true,odd=false,even=true,maxintermodorder=Inf,
    )

    # the original solver kept one mode fewer when Nsignalmodes is even;
    # drop the highest one to match
    if mod(Nsignalmodes,2) == 0 && Nsignalmodes > 0
        signalfreq = JosephsonCircuits.removefreqs(
            signalfreq,
            [(Nsignalmodes,)],
        )
    end

    linearized = hblinsolve(ws, psc, circuitdefs, signalfreq;
        nonlinear = nonlinear, nbatches = nbatches,
        returnS = returnS, returnSnoise = returnSnoise, returnQE = returnQE,
        returnCM = returnCM, returnnodeflux = returnnodeflux,
        returnnodefluxadjoint = returnnodefluxadjoint,
        returnvoltage = returnvoltage,
        returnvoltageadjoint = returnvoltageadjoint,
        keyedarrays = keyedarrays, sensitivitynames = sensitivitynames,
        returnSsensitivity = returnSsensitivity, returnZ = returnZ,
        returnZadjoint = returnZadjoint, returnZsensitivity = returnZsensitivity,
        returnZsensitivityadjoint = returnZsensitivityadjoint,
        factorization = factorization)

    return HB(nonlinear, linearized)
end

# `solveS!` as v0.5.4 took it, without the fill reducing ordering that
# `solveS_initialize` now chooses once and returns last in its tuple: the
# ordering is chosen afresh by each batch's factorization
function solveS!(Se, Si, Ce, Ci, portse, portsi, gammaii, See, Sei, Sie, Sii,
    See_indices, Sei_indices, Sie_indices, Sii_indices, gammaii_indexmap,
    Sii_indexmap, scattering_parameters, noise_covariances, nbatches,
    factorization, internal_ports, noise)
    Base.depwarn("`solveS!` takes the fill reducing ordering `solveS_initialize` returns as its last argument; call it as `solveS!(init...)` with the whole tuple `solveS_initialize` returns to avoid errors in future versions.", :solveS!; force = true)
    return solveS!(Se, Si, Ce, Ci, portse, portsi, gammaii, See, Sei, Sie,
        Sii, See_indices, Sei_indices, Sie_indices, Sii_indices,
        gammaii_indexmap, Sii_indexmap, scattering_parameters,
        noise_covariances, nbatches, factorization, internal_ports, noise,
        nothing)
end
