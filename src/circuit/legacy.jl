# The deprecated input formats.
#
# The typed `Circuit` whose frequency dependent values are closures of the
# frequency is the input format of the package. Two older ways of writing
# a circuit are deprecated, and each is converted here, with a warning,
# into one written the current way, which then follows the same path.
# Everything specific to them lives in this file, so that they are removed
# by deleting it and test/circuit/legacy.jl.
#
# A netlist of `(name, node1, node2, value)` tuples, the original input
# format: the name prefix table with the two functions that read it, the
# convention that a port's reference impedance is the resistor placed
# across it and the port termination which records that resistor, the
# tuple forms of every entry point, which sort the nodes by number as the
# format always did, and the tuple netlist file reader and writer.
#
# A value written as an expression in a parameter named by the
# `symfreqvar` keyword of the solvers, rewritten as a
# `FrequencyDependent` closure.
#
# Outside this file, the `Circuit(netlist)` constructor in
# circuit/parse.jl hands a netlist whose entries end in values rather than
# components to `legacycircuit`, the Symbolics extension unwraps a `Num`
# port number for `legacyportnumber`, and the three solver entry points
# which hold a compiled circuit beside its definitions hand a
# `symfreqvar` to `frequencydependentcircuit`.

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
        return NonlinearInductor(value, sin, cos)
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
    Base.depwarn(tuplenetlistmessage, caller; force = true)
    if !isnothing(pins) || !isnothing(ports)
        throw(ArgumentError("A tuple netlist has no interface; give pins and ports to a netlist of typed components or to the connection-group form."))
    end
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
    values = Any[frequencydependentvalue(v, circuitdefs, symfreqvar)
        for v in psc.componentvalues]
    return CompiledCircuit(psc.nodenames, psc.nodeindices, psc.Nnodes,
        psc.componentnames, psc.componenttypes, values,
        psc.componentnamedict, psc.componenttemperatures, psc.junctioncprs,
        psc.capacitors, psc.resistors,
        psc.inductors, psc.junctions, psc.currentsources,
        psc.mutualinductors, psc.ports, psc.scatteringblocks,
        psc.couplings, psc.topology)
end

function frequencydependentvalue(value, circuitdefs, symfreqvar)
    checkissymbolic(value) || return value
    any(v -> isequal(v, symfreqvar), circuitvariables(value)) || return value
    # the definitions resolved once, so that only the frequency is left for
    # the closure to give
    partial = valuetonumber(value, circuitdefs)
    if checkissymbolic(partial) &&
            any(v -> !isequal(v, symfreqvar), circuitvariables(partial))
        return value
    end
    # the frequency reaches the value two ways: as the symbolic variable it
    # substitutes, and as the argument of any frequency dependent leaf the
    # same expression already carries, which `substitutefreq` evaluates.
    return FrequencyDependent(
        w -> substitutefreq(
            valuetonumber(partial, Dict{Any,Any}(symfreqvar => w)), w))
end

# === the tuple forms of the entry points ===
#
# Each converts the netlist, which warns, compiles it with the nodes
# sorted by number as the tuple format always did, and forwards the
# compiled circuit to the typed form.

compile(netlist::AbstractVector; sorting::Symbol = :number) =
    compile(legacycircuit(netlist, nothing; caller = :compile); sorting = sorting)

# the netlist compiled the way the tuple format ordered its nodes
legacycompiled(netlist, caller::Symbol) =
    compile(legacycircuit(netlist, nothing; caller = caller); sorting = :number)

function hbsolve(ws, wp::NTuple{N,Number}, sources,
        Nmodulationharmonics::NTuple{M,Int}, Npumpharmonics::NTuple{N,Int},
        netlist::AbstractVector, circuitdefs::AbstractDict = Dict{Symbol,Any}();
        kwargs...) where {N,M}
    return hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        legacycompiled(netlist, :hbsolve), circuitdefs; kwargs...)
end

function hbnlsolve(w::NTuple{N,Number}, Nharmonics::NTuple{N,Int}, sources,
        netlist::AbstractVector, circuitdefs::AbstractDict = Dict{Symbol,Any}();
        kwargs...) where {N}
    return hbnlsolve(w, Nharmonics, sources,
        legacycompiled(netlist, :hbnlsolve), circuitdefs; kwargs...)
end

function hblinsolve(w, netlist::AbstractVector,
        circuitdefs::AbstractDict = Dict{Symbol,Any}(); kwargs...)
    return hblinsolve(w, legacycompiled(netlist, :hblinsolve), circuitdefs;
        kwargs...)
end

function transientproblem(netlist::AbstractVector,
        circuitdefs::AbstractDict = Dict{Symbol,Any}(); sources = ())
    return transientproblem(legacycompiled(netlist, :transientproblem),
        circuitdefs; sources = sources)
end

function numericmatrices(netlist::AbstractVector, circuitdefs::Dict;
        Nmodes::Int = 1)
    return numericmatrices(legacycompiled(netlist, :numericmatrices),
        circuitdefs; Nmodes = Nmodes)
end

function symbolicmatrices(netlist::AbstractVector; Nmodes::Int = 1)
    return symbolicmatrices(legacycompiled(netlist, :symbolicmatrices);
        Nmodes = Nmodes)
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
        value = try
            parse(Float64,split_line[4])
        catch
            # https://docs.sciml.ai/Symbolics/stable/manual/parsing/
            # Symbolics.parse_expr_to_symbolic(Meta.parse(split_line[4]),Main)
            parsecomponentvalue(split_line[4])
        end
        push!(circuit,(split_line[1],split_line[2],split_line[3],value))
    end
    return nothing
end



# export_netlist("test1.net", circuit,circuitdefs)

# Reading a component value back out of a netlist line. Only
# `import_netlist!` above uses it.
"""
    parsecomponentvalue(s::AbstractString)

Parse a SPICE netlist component value into a number or a `CircuitValue`.
Replaces `Symbolics.parse_expr_to_symbolic`; unlike it, this does not
evaluate into a module, so a netlist cannot introduce arbitrary code.
"""
function parsecomponentvalue(s::AbstractString)
    v = tryparse(Float64, s)
    isnothing(v) || return v
    return CircuitValues.fromexpr(Meta.parse(s))
end
