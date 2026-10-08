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
# format: the name prefix table with the function that reads it, the
# convention that a port's reference impedance is the resistor placed
# across it and the port termination which records that resistor, the
# tuple forms of the entry points which took a tuple netlist in v0.5.4,
# which sort the nodes by number as the format always did, and the tuple
# netlist file reader, with its parser of a value written as an
# expression, and writer.
#
# A value written as an expression in a parameter named by the
# `symfreqvar` keyword of the solvers, rewritten as a
# `FrequencyDependent` closure.
#
# The keywords of the nonlinear solvers which v0.5.4 had and the solvers
# no longer read: `ftol`, the absolute residual tolerance under its old
# name, the line search settings `switchofflinesearchtol` and `alphamin`,
# and `maxharmonics` and `maxpumpharmonics`, whose role the retained
# harmonics took; the `returnZ` keywords of the linearized solvers, whose
# impedances they no longer compute; and the `noise` keyword of
# `connectS_initialize`.
#
# Outside this file: the `Circuit(netlist)` constructor in
# circuit/parse.jl hands a netlist whose entries end in values rather than
# components to `legacycircuit`; `compile` (circuit/compile.jl) finds the
# resistor a port's termination names through `namedtermination`, whose
# method for every other termination, naming none, is in
# circuit/components.jl; the Symbolics extension unwraps a `Num` port
# number for `legacyportnumber`; the three solver entry points which hold
# a compiled circuit beside its definitions hand a `symfreqvar` to
# `frequencydependentcircuit`; `hbnlsolve` and `hbsolve` hand their
# deprecated keywords to `deprecatedsolverkeywords`, and `hbnlsolve` its
# `factorization` to `deprecatedfactorization`; `hbsolve` and
# `hblinsolve` hand the `returnZ` keywords to `removedimpedancekeywords`;
# `hbnlsolve` and `hbsolve` collect what these find and warn once with it
# (`warndeprecations`); and `connectS_initialize`
# (networks/connections.jl) hands its `noise` to `deprecatedconnectnoise`.

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

# The component type a tuple netlist entry's name gives by its prefix:
# the first of the prefixes which begins the name, so that `Lj`, a
# junction, is read before `L`, an inductor.
const legacyprefixes = (:Lj, :L, :C, :K, :I, :R, :P)
function legacycomponenttype(name::AbstractString)
    for t in legacyprefixes
        startswith(name, String(t)) && return t
    end
    throw(ArgumentError(lazy"The component name $(name) begins with none of the prefixes of the tuple netlist, Lj, L, C, K, I, R and P, which give its type."))
end

# the typed component model of a tuple netlist entry of the type `t`; a
# mutual inductor's entry names the inductors it couples in place of its
# nodes, as the netlist form of a `Circuit` does
function legacycomponent(t::Symbol, name, value)
    t === :L && return Inductor(value)
    t === :C && return Capacitor(value)
    t === :R && return Resistor(value)
    t === :Lj && return NonlinearInductor(value, sin)
    t === :I && return CurrentSource(value)
    t === :P && return Port(legacyportnumber(name, value); termination = nothing)
    return MutualInductor(value)
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

# Convert the tuple netlist to a `Circuit`, its deprecation going `to` the
# entry point the netlist was given to or to the messages of the call (see
# `deprecate!`): each entry becomes the entry of the netlist form of a
# `Circuit` with its typed component, its nodes labelled by their text, a
# port with the resistor across it as its termination, and the netlist
# form groups the nodes into nets, in order of first appearance, with
# `Ground` appended to net "0".
function legacycircuit(netlist, circuitdefs; pins = nothing, ports = nothing,
        to = :Circuit)
    # a list of typed components handed over without the Circuit around it,
    # or components without their connections, is not a tuple netlist, and
    # is said to be so before a warning about one
    if any(e -> e isa Pair && last(e) isa AbstractComponent, netlist)
        throw(ArgumentError("The list holds name => component pairs, which are the components of a circuit: give their connections as well, as Circuit(components, connections)."))
    end
    if any(e -> e isa Tuple && last(e) isa AbstractComponent, netlist)
        throw(ArgumentError("The netlist holds typed components, which are passed as a circuit: write Circuit(netlist)."))
    end
    deprecate!(to, tuplenetlistmessage)
    if !isnothing(pins) || !isnothing(ports)
        throw(ArgumentError("A tuple netlist has no interface; give pins and ports to a netlist of typed components or to the connection-group form."))
    end
    components = Vector{Pair{String,Any}}(undef, length(netlist))
    for (i, entry) in enumerate(netlist)
        if length(entry) != 4
            throw(ArgumentError(lazy"The netlist entry $(entry) on line $(i) must be a tuple of (name, node1, node2, value)."))
        end
        name, _, _, value = entry
        if !(name isa AbstractString)
            throw(ArgumentError(lazy"The component name $(name) on line $(i) must be a string."))
        end
        if !isnothing(circuitdefs)
            value = valuetonumber(value, circuitdefs)
        end
        components[i] = Pair{String,Any}(String(name),
            legacycomponent(legacycomponenttype(name), name, value))
    end
    legacyportimpedances!(components, netlist)
    return netlistcircuit([(c.first, string(entry[2]), string(entry[3]), c.second)
        for (c, entry) in zip(components, netlist)])
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
    frequencydependentcircuit(psc, circuitdefs, symfreqvar, to)

A compiled circuit whose values depending on the symbolic frequency
variable `symfreqvar` are rewritten as [`FrequencyDependent`](@ref)
closures resolving them at `circuitdefs` and at the frequency, its
deprecation going `to` the entry point or the messages of the call (see
[`warndeprecations`](@ref)).

Every other value is left as it was, including one which still depends on
a parameter the definitions do not give, so that the undefined parameter
is still reported by [`checkcomponentvaluesdefined`](@ref) naming its
component rather than failing later inside a closure.
"""
function frequencydependentcircuit(psc::CompiledCircuit, circuitdefs,
        symfreqvar, to)
    deprecate!(to, symfreqvarmessage)
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

# === the deprecations of one call ===
#
# `Base.depwarn` shows a warning once per calling frame and function name,
# so of two deprecations one call of an entry point meets, the second
# would not be shown. A deprecation therefore goes `to` one of two
# places: the name of the entry point, which warns at once, when the call
# can meet no other, or the messages of a call which can meet several,
# which warns once with all of them (`warndeprecations`).

deprecate!(caller::Symbol, message) =
    (Base.depwarn(message, caller; force = true); nothing)
deprecate!(messages::Vector{String}, message) =
    (push!(messages, message); nothing)

"""
    warndeprecations(messages::Vector{String}, caller::Symbol)

Warn once, attributed to `caller`, with `messages`, the deprecations one
call of `caller` met; nothing when there are none.
"""
function warndeprecations(messages::Vector{String}, caller::Symbol)
    isempty(messages) ||
        Base.depwarn(join(messages, "\n"), caller; force = true)
    return nothing
end

# === the deprecated keywords of the nonlinear solvers ===

"""
    deprecatedsolverkeywords(to, atol; ftol = nothing,
        switchofflinesearchtol = nothing, alphamin = nothing,
        maxharmonics = nothing, maxpumpharmonics = nothing)

The absolute residual tolerance of a solve given the deprecated keywords
of the nonlinear solvers: `ftol` when it is given, which is `atol` under
its old name, and `atol` otherwise. The deprecation of each keyword given
goes `to` the entry point or the messages of the call (see
[`warndeprecations`](@ref)); `switchofflinesearchtol` and `alphamin`,
which the line search no longer has, and `maxharmonics` of `hbnlsolve`
and `maxpumpharmonics` of `hbsolve`, whose role the retained harmonics
took, are ignored.
"""
function deprecatedsolverkeywords(to, atol; ftol = nothing,
        switchofflinesearchtol = nothing, alphamin = nothing,
        maxharmonics = nothing, maxpumpharmonics = nothing)
    if !isnothing(ftol)
        deprecate!(to, "The `ftol` kwarg is deprecated: the absolute residual tolerance is `atol` in every solver of the package. Please use `atol` to avoid errors in future versions.")
        atol = ftol
    end
    unused(name) = deprecate!(to, lazy"The `$(name)` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions.")
    isnothing(switchofflinesearchtol) || unused(:switchofflinesearchtol)
    isnothing(alphamin) || unused(:alphamin)
    retained(name, by) = deprecate!(to, lazy"The `$(name)` kwarg is deprecated and no longer used. `$(by)` is the retained set of modes and `Nevaluationharmonics` the grid on which the nonlinearity is sampled. Please remove it to avoid errors in future versions.")
    isnothing(maxharmonics) || retained(:maxharmonics, :Nharmonics)
    isnothing(maxpumpharmonics) || retained(:maxpumpharmonics, :Npumpharmonics)
    return atol
end

"""
    deprecatedfactorization(to, method, factorization)

The method of a solve given `factorization`, the keyword with which
`hbnlsolve` took the sparse factorization of its Newton iteration in
v0.5.4. The factorization is now an option of the method: given beside
the default, a [`NewtonKrylov`](@ref), its deprecation goes `to` the
entry point or the messages of the call (see [`warndeprecations`](@ref)),
and the solve takes `Newton(factorization = factorization)`, as v0.5.4
did; given beside another method it is refused.
"""
deprecatedfactorization(to, method, ::Nothing) = method
function deprecatedfactorization(to, method, factorization)
    method isa NewtonKrylov || throw(ArgumentError(lazy"`factorization` is an option of the method, and `method` = $(method) is given as well; pass the factorization to the method alone."))
    deprecate!(to, "The `factorization` kwarg is deprecated: the factorization is an option of the method. Please pass `method = Newton(factorization = ...)` to avoid errors in future versions.")
    return Newton(; factorization)
end

# === the tuple forms of the entry points ===
#
# Each entry point of v0.5.4 which took a tuple netlist converts it, which
# warns, compiles it with the nodes sorted by number as the tuple format
# always did, and forwards the compiled circuit to the typed form. The
# solvers and the matrix builders of the tuple format took the order as
# their keyword `sorting`, which is now `compile`'s, and still take it
# here.

# the netlist compiled the way the tuple format ordered its nodes, its
# deprecation going `to` the entry point or the messages of the call. The
# conversion is called through an inference barrier: a solver called with
# a circuit inference cannot see, such as a global read in a function,
# infers its tuple form too, and would infer the conversion of a netlist
# of unknown entries with it
legacycompiled(netlist, to; sorting::Symbol = :number) =
    Base.inferencebarrier(compilenetlist)(netlist, to, sorting)::CompiledCircuit
compilenetlist(netlist, to, sorting) =
    compile(legacycircuit(netlist, nothing; to = to); sorting = sorting)

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
    return expressionvalue(Meta.parse(s))
end

# A value written in a netlist file arrives as a parsed Julia `Expr`, which
# is converted to a `CircuitValue` accepting only the closed operator set of
# circuit/values.jl and the constants `im` and `pi`, so that no Symbolics
# parser is needed.
const expressionconstants = Dict{Symbol,Any}(:im => im, :pi => pi)
const expressionbinary = Dict{Symbol,Function}(:+ => +, :- => -, :* => *, :/ => /, :^ => ^)
const expressionunary = Dict{Symbol,Function}(:- => -, :inv => inv, :sqrt => sqrt,
    :exp => exp, :log => log, :conj => conj, :real => real, :imag => imag)
expressionvalue(x::Number) = CircuitValues.Constant(x)
function expressionvalue(s::Symbol)
    haskey(expressionconstants, s) && return CircuitValues.Constant(expressionconstants[s])
    return CircuitValues.Parameter(s)
end
function expressionvalue(e::Expr)
    e.head === :call || error("unsupported expression head $(e.head)")
    op = e.args[1]; args = map(expressionvalue, e.args[2:end])
    length(args) == 1 && haskey(expressionunary, op) && return expressionunary[op](args[1])
    haskey(expressionbinary, op) && return reduce(expressionbinary[op], args)
    error("unsupported operator $(op) in a component value")
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

# The `noise` keyword `connectS_initialize` took in v0.5.4, which has no
# effect: `connectS!` computes the noise covariances of the networks given
# without them when it is called with `noise = true`. Given, it warns.
function deprecatedconnectnoise(noise)
    isnothing(noise) || Base.depwarn("The `noise` kwarg of `connectS_initialize` is deprecated and has no effect: `connectS!` computes the noise covariances of the networks given without them, from their scattering parameters at the time, when it is called with `noise = true`. Please remove it to avoid errors in future versions.", :connectS_initialize; force=true)
    return nothing
end

# `X_Y_to_sympletic_pair` and `X_Y_to_sympletic_block` were misspelled.
function X_Y_to_sympletic_pair(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real}; kwargs...)
    Base.depwarn("`X_Y_to_sympletic_pair` is deprecated, use `X_Y_to_symplectic_pair` instead.", :X_Y_to_sympletic_pair; force=true)
    return X_Y_to_symplectic_pair(X, Y; kwargs...)
end

function X_Y_to_sympletic_block(X::AbstractMatrix{<:Real}, Y::AbstractMatrix{<:Real}; kwargs...)
    Base.depwarn("`X_Y_to_sympletic_block` is deprecated, use `X_Y_to_symplectic_block` instead.", :X_Y_to_sympletic_block; force=true)
    return X_Y_to_symplectic_block(X, Y; kwargs...)
end


# the deprecation of the single pump form of `hbsolve`, below
const singlepumpmessage = "This form of hbsolve, with a single pump frequency and integer harmonic counts, is deprecated: it calls the harmonic balance solvers hbnlsolve and hblinsolve, which take any number of pump tones and ports, with the syntax of the legacy solver, which supported four wave mixing of one strong tone only. Please switch to hbsolve(ws, (wp,), sources, (Nmodulationharmonics,), (Npumpharmonics,), circuit, circuitdefs)."

#     hbsolve(ws, wp, Ip, Nsignalmodes::Int, Npumpmodes::Int, circuit,
#         circuitdefs; pumpports = [1], keyword arguments...)
#
# The original `hbsolve` signature: a single pump at the scalar frequency
# `wp`, applied to `pumpports` with the currents `Ip`, four wave mixing
# only, and mode counts given as integers rather than tuples. It is
# translated into a call of the current solvers and warns that it is
# deprecated. A tuple netlist is compiled with its nodes sorted as
# `sorting` says, and any other circuit as `compile` compiles it; the line
# search keywords and the impedance outputs are ignored. One warning names
# the form and every other deprecation the call meets. Note that the
# signal modes of the result are ordered as `hblinsolve` orders them (the
# signal at index 1, the rest as listed in `modes`), which need not match
# the order the original solver used.
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

    # the deprecations of this call, warned together below
    deprecations = String[singlepumpmessage]

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
        legacycompiled(circuit, deprecations; sorting = sorting) :
        compile(circuit)
    isnothing(symfreqvar) || (psc = frequencydependentcircuit(psc,
        circuitdefs, symfreqvar, deprecations))
    deprecatedsolverkeywords(deprecations, ftol; switchofflinesearchtol,
        alphamin)
    removedimpedancekeywords(deprecations; returnZ, returnZadjoint,
        returnZsensitivity, returnZsensitivityadjoint)
    warndeprecations(deprecations, :hbsolve)
    nm=numericmatrices(psc, circuitdefs, Nmodes = Nmodes)

    # `factorization = nothing` leaves the preconditioner to `Automatic`; a
    # factorization pins the block diagonal, the only preconditioner this
    # form builds with one
    nonlinear = hbnlsolve(w, sources, freq, indices, psc, nm;
        iterations = iterations, atol = ftol,
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
        returnSsensitivity = returnSsensitivity,
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

# === the impedance outputs of the linearized solve, removed ===

"""
    removedimpedancekeywords(to; returnZ = nothing, returnZadjoint = nothing,
        returnZsensitivity = nothing, returnZsensitivityadjoint = nothing)

A deprecation, going `to` the entry point or the messages of the call
(see [`warndeprecations`](@ref)), when any of the impedance outputs
v0.5.4's `hbsolve` and `hblinsolve` took is given: no output replaces
them, and they are otherwise ignored. `hbsolve`, its single pump form
and `hblinsolve` hand them here.
"""
function removedimpedancekeywords(to; returnZ = nothing,
        returnZadjoint = nothing, returnZsensitivity = nothing,
        returnZsensitivityadjoint = nothing)
    all(isnothing, (returnZ, returnZadjoint, returnZsensitivity,
        returnZsensitivityadjoint)) && return nothing
    deprecate!(to, "The `returnZ`, `returnZadjoint`, `returnZsensitivity`, and `returnZsensitivityadjoint` kwargs have been removed. Please compute them from scattering parameters matrices.")
    return nothing
end

# === the Xyce harmonic balance reader, deprecated ===
#
# `spice_hb_load` reads the frequency domain output of a Xyce harmonic
# balance simulation, a `.HB.FD.prn` file, and is deprecated, since the
# package neither writes nor runs Xyce. It returns a named tuple with the
# value of each output variable at each frequency in `data`, one row per
# variable in the order of the header, the name of each row in
# `variables`, the frequencies in `f`, the index column in `index`, and
# the column names in `header`. A variable printed as the two columns
# `Re(name)` and `Im(name)` is the row `name`; any other column is a row
# of its own, its real part.
function spice_hb_load(filename)
    Base.depwarn("`spice_hb_load` is deprecated: the package neither writes nor runs Xyce. It still reads a Xyce `.HB.FD.prn` file, and will be removed in a future version.", :spice_hb_load; force = true)

    data = Float64[]
    header = SubString{String}[]

    open(filename, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            if startswith(s, "Index")
                append!(header, split(s, r"\s+"))
            elseif s == "End of Xyce(TM) Simulation"
                break
            else
                append!(data, parse.(Float64, split(s, r"\s+")))
            end
        end
    end

    values = reshape(data, length(header), :)
    index = values[1, :]
    f = values[2, :]

    # the columns of each variable, its real and imaginary parts by name
    variables = String[]
    recolumn = Int[]
    imcolumn = Int[]
    for c in 3:length(header)
        m = match(r"^(Re|Im)\((.*)\)$", header[c])
        name = isnothing(m) ? String(header[c]) : String(m.captures[2])
        k = findfirst(==(name), variables)
        if isnothing(k)
            push!(variables, name)
            push!(recolumn, 0)
            push!(imcolumn, 0)
            k = length(variables)
        end
        if !isnothing(m) && m.captures[1] == "Im"
            imcolumn[k] = c
        else
            recolumn[k] = c
        end
    end

    data1 = zeros(Complex{Float64}, length(variables), size(values, 2))
    for k in eachindex(variables)
        recolumn[k] > 0 && (data1[k, :] .+= view(values, recolumn[k], :))
        imcolumn[k] > 0 && (data1[k, :] .+= im .* view(values, imcolumn[k], :))
    end

    return (data=data1, f=f, index=index, header=header, variables=variables)
end
