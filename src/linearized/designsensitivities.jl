# =========================================================================
# Design parameter sensitivities of a typed circuit.
#
# A circuit's component values are written in terms of design parameters
# (`Capacitor(:Cc)`, `JosephsonJunction(phi0/Ic)` with `Ic` a parameter of
# `@params`), and the definitions give every parameter a number. The
# derivative of a component value with respect to a parameter is then the
# derivative of the expression it was written as, exact, and the
# derivative of the scattering parameters follows by the chain rule through
# the component sensitivities:
#
#     dS/dp_j = sum_k (dS/dv_k) (dv_k/dp_j),
#
# reverse mode through the solve, once, via the adjoint component
# sensitivities, and the exact direction of every dependent value carried
# into the contraction.
# =========================================================================

"""
    designderivative(value, name::Symbol, definitions)

The derivative of a component value, as it is written, with respect to the
parameter `name` it is written in, at `definitions`, as a number: zero for
a number, one for the parameter itself written as a symbol or a string,
and for a [`CircuitValue`](@ref) expression the derivative of the
expression with the definitions substituted (see
[`resolvedefinitions`](@ref)). A frequency dependent leaf depends on no
parameter. A derivative which does not come to a number, one which
depends on the frequency or on a parameter the definitions do not give, is
refused. This differentiates the value alone: when `name` is itself
defined in terms of other parameters, [`DesignChain`](@ref) carries the
derivative on through its definition. The Symbolics extension adds the
method for a `Num`.
"""
designderivative(::Number, ::Symbol, definitions) = zero(ComplexF64)
designderivative(v::Union{Symbol,AbstractString}, name::Symbol, definitions) =
    definitionname(v) === name ? one(ComplexF64) : zero(ComplexF64)
function designderivative(v::CircuitValue, name::Symbol, definitions)
    d = valuetonumber(CircuitValues.derivative(v, name), definitions)
    isnumeric(d) || throw(ArgumentError(lazy"the derivative of the value $(v) with respect to $(name) is $(d) at these definitions, which is not a number: a design parameter must move a value in a direction which depends neither on the frequency nor on an undefined parameter."))
    return ComplexF64(d)
end

"""
    DesignChain(definitions, names)

The derivatives of the parameters of `definitions` with respect to the
design parameters `names`, by the chain rule through the definitions:
a parameter's derivative is one with respect to itself, and a parameter
defined by a value written in other parameters adds the derivative of
that value in each of them ([`designderivative`](@ref)) times that one's
own. So a component written `Capacitor(:Cd)` with `:Cd => 2*C0` has the
derivative 2 with respect to `C0`, and a selected parameter which is itself
defined in terms of others is differentiated as the definitions would move
it, its definition replaced by its value. The derivatives are evaluated at
the definitions resolved once ([`resolvedefinitions`](@ref)), the
resolution the values themselves take, kept sparse, as pairs of a
parameter's index in `names` and a nonzero derivative, and memoized by
name, so a parameter shared by many values is differentiated once.
"""
struct DesignChain
    definitions::Dict{Symbol,Any}
    resolved::ResolvedDefinitions
    selected::Dict{Symbol,Int}
    gradients::Dict{Symbol,Vector{Pair{Int,ComplexF64}}}
end
function DesignChain(definitions, names::AbstractVector{Symbol})
    byname = definitionsbyname(definitions)
    return DesignChain(byname, resolvedefinitions(byname),
        Dict(n => j for (j, n) in enumerate(names)),
        Dict{Symbol,Vector{Pair{Int,ComplexF64}}}())
end

# the derivative of the parameter `name`: one with respect to itself when
# it is selected, and the derivative of its definition when that is written
# in other parameters
function namegradient!(chain::DesignChain, name::Symbol)
    g = get(chain.gradients, name, nothing)
    isnothing(g) || return g
    v = get(chain.definitions, name, nothing)
    g = checkissymbolic(v) ? valuegradient!(chain, v) :
        Pair{Int,ComplexF64}[]
    j = get(chain.selected, name, 0)
    iszero(j) || (g = vcat(j => one(ComplexF64), g))
    chain.gradients[name] = g
    return g
end

# the derivative of a value as it is written: its derivative in each
# parameter it names times that parameter's own, summed by design parameter
function valuegradient!(chain::DesignChain, v)
    g = Pair{Int,ComplexF64}[]
    for n in valuenames(v)
        gn = namegradient!(chain, n)
        isempty(gn) && continue
        d = designderivative(v, n, chain.resolved)
        iszero(d) && continue
        for (j, x) in gn
            k = findfirst(p -> first(p) == j, g)
            if isnothing(k)
                push!(g, j => d*x)
            else
                g[k] = j => last(g[k]) + d*x
            end
        end
    end
    return filter!(p -> !iszero(last(p)), g)
end

# The value of the component `i` at the definitions, as a complex number,
# and its derivative with respect to the selected parameters, or `nothing`
# and an empty derivative for a frequency dependent value which no selected
# parameter moves. A value which depends on a parameter the definitions do
# not give is refused, as is a frequency dependent value which a selected
# parameter moves: the sensitivities rescale a component's own stamp,
# which needs its value to be a number.
function componentgradient!(chain::DesignChain, psc::CompiledCircuit,
        i::Integer)
    written = psc.componentvalues[i]
    value = resolvevalue(written, chain.resolved)
    if !isnumeric(value)
        (checkissymbolic(value) && isempty(circuitvariables(value))) ||
            throw(ArgumentError(lazy"the component $(psc.componentnames[i]) has the value $(value) at these definitions, which is not a number; design sensitivities need every parameter defined."))
        isempty(valuegradient!(chain, written)) || throw(ArgumentError(
            lazy"the component $(psc.componentnames[i]) has the frequency dependent value $(value), which a selected design parameter moves; a design parameter may move only components whose values are numbers."))
        return nothing, Pair{Int,ComplexF64}[]
    end
    return ComplexF64(value), valuegradient!(chain, written)
end

# the analytic derivatives of a block's scattering matrix, as
# `name => provider` pairs: only a `ScatteringParameters` block states any
blockderivatives(d::ScatteringParameters) = d.derivatives
blockderivatives(d) = Pair{Symbol,Any}[]

# the names of the design parameters: the ones given, as names or as
# definition keys, or by default every parameter the definitions give a
# number and every parameter a scattering block states a derivative for,
# in sorted order
function designparameters(parameters, definitions, psc::CompiledCircuit)
    if isnothing(parameters)
        names = Symbol[]
        for (n, v) in definitionsbyname(definitions)
            isnumeric(v) && push!(names, n)
        end
        for b in psc.scatteringblocks
            append!(names, first.(blockderivatives(b.definition)))
        end
        return sort!(unique!(names))
    end
    names = Symbol[]
    for q in parameters
        n = definitionname(q)
        isnothing(n) && throw(ArgumentError(lazy"$(q) is not the name of a design parameter."))
        push!(names, n)
    end
    return names
end

"""
    designdirections(psc::CompiledCircuit, definitions, names)

The nonzero derivatives of the component values of `psc` with respect to
the design parameters `names`, from one sparse enumeration: each value is
differentiated once, through the definitions by the chain rule
([`DesignChain`](@ref)), and only the derivatives which are not zero are
kept, so the cost grows with the dependences rather than with the
components times the parameters. Returns `(components, values, entries)`:
the flat indices of the components a design parameter may move, in
compiled order, their values at the definitions, and a triplet
`(k, j, dv/dp)` per nonzero derivative of the `k`th of them with respect
to parameter `j`, by parameter and then by component. The components are
those of [`designjacobian`](@ref), which materializes its dense output
from these.
"""
function designdirections(psc::CompiledCircuit, definitions,
        names::AbstractVector{Symbol})
    chain = DesignChain(definitions, names)
    # a port's reference impedance is the value of the termination it owns,
    # a component of its own, or of the port itself when it owns none
    terminated = Set(p.component for p in psc.ports if !iszero(p.environment))
    components = Int[]
    values = ComplexF64[]
    entries = Tuple{Int,Int,ComplexF64}[]
    for i in eachindex(psc.componentvalues)
        i in terminated && continue
        v, g = componentgradient!(chain, psc, i)
        isnothing(v) && continue
        push!(components, i)
        push!(values, v)
        for (j, x) in g
            push!(entries, (length(components), j, x))
        end
    end
    sort!(entries; by = e -> (e[2], e[1]))
    return (; components, values, entries)
end

"""
    designjacobian(circuit, circuitdefs; parameters = nothing)

The Jacobian of the component values of a typed [`Circuit`](@ref), or of
its compiled circuit, with respect to its design parameters, exact, at the
definitions `circuitdefs`: the entry `(k, j)` is
`d(value of component k)/d(parameter j)`. A definition is keyed by the
parameter's symbol, its string, its parameter of [`@params`](@ref) or,
with Symbolics, its `Num`, and is one of

- a number, real or complex: an independent parameter, whose derivative
  with respect to itself is one;
- a value written in other parameters, an expression of them or a
  Symbolics expression, to any depth (`:Cc => Cj/10`): its derivative
  follows the definition by the chain rule ([`DesignChain`](@ref)), so
  `Capacitor(:Cc)` has the derivative `1/10` with respect to `Cj`;
- a frequency dependent value ([`FrequencyDependent`](@ref)), which no
  parameter moves.

A component value is a number, a name, an expression in names (`Cj/4`) or
a frequency dependent value, and its derivative is the expression's,
evaluated at the definitions, complex wherever a parameter moves a value
along a complex direction (a loss tangent rotates it).

Returns `(names, values, J)`: the component names in compiled order, their
values at the definitions, and the complex `length(names)` by
`length(parameters)` Jacobian, dense, materialized from the sparse
enumeration of the nonzero derivatives ([`designdirections`](@ref)).
`parameters` are names or definition keys; by default every parameter the
definitions give a number and every parameter a
[`ScatteringParameters`](@ref) block states a derivative for, sorted by
name. A parameter defined in terms of others may be selected by name, and
is differentiated as replacing its definition by its value would move it,
the parameters it is written in held. A port's slot holds its reference
impedance: a port which owns a termination is not among the components,
its reference impedance being the value of that termination, which is; a
port which owns none (`termination = nothing`) is, and a parameter moves
only the normalization of its waves through it. A component whose value
is frequency dependent is not among them either, and one a selected
parameter moves is refused: the sensitivities rescale a component's own
stamp, which needs its value to be a number. A scattering block has no
scalar value to differentiate; its dependence is carried separately by
[`designblockjacobian`](@ref).
"""
function designjacobian(circuit::CompilableCircuit, circuitdefs::AbstractDict;
        parameters = nothing)
    psc = compile(circuit)
    definitions = definitiontable(circuitdefs)
    names = designparameters(parameters, definitions, psc)
    d = designdirections(psc, definitions, names)
    J = zeros(ComplexF64, length(d.components), length(names))
    for (k, j, x) in d.entries
        J[k, j] = x
    end
    return String[psc.componentnames[i] for i in d.components], d.values, J
end

"""
    designblockjacobian(circuit, parameters)

The scattering block dependence of a circuit on its design parameters: for
each [`ScatteringParameters`](@ref) block of the compiled circuit and each
parameter its `derivatives` names, a derivative block whose scattering
matrix is the block's stated `dS/dp`. A block which states no derivative
for a parameter does not depend on it. Returns a vector of
`(blockpath, parameterindex, derivativeblock)`, by block and, within a
block, in the order of `parameters`. Each block's derivatives are read
once, against a table of the parameters' positions.
"""
function designblockjacobian(circuit::CompilableCircuit, parameters)
    psc = compile(circuit)
    positions = Dict{Symbol,Vector{Int}}()
    for (j, name) in enumerate(parameters)
        push!(get!(Vector{Int}, positions, name), j)
    end
    out = Tuple{String,Int,Any}[]
    found = Tuple{Int,Any}[]
    for b in psc.scatteringblocks
        empty!(found)
        for (name, provider) in blockderivatives(b.definition),
                j in get(positions, name, ())
            push!(found, (j, provider))
        end
        for (j, provider) in sort!(found; by = first)
            # a derivative is not a passive scattering matrix and is not
            # checked as one: the positional constructor, called with the
            # provider the block holds untyped, builds a block of its
            # concrete type
            push!(out, (b.path, j, ScatteringParameters(provider,
                b.definition.nports, b.definition.zref, b.definition.grounded,
                b.definition.noise, b.definition.negative_frequency)))
        end
    end
    return out
end

"""
    designsensitivities(circuit, circuitdefs, ws, wp, sources,
        Nmodulationharmonics, Npumpharmonics; parameters = nothing,
        keyedarrays = true, kwargs...)

The derivative of the scattering parameters with respect to the design
parameters of a typed [`Circuit`](@ref) whose component values are
written in terms of them, by the chain rule through the component
sensitivities:

    dS/dp_j = sum_k (dS/dv_k) (dv_k/dp_j).

`circuitdefs` defines the parameters in the forms
[`designjacobian`](@ref) lists: by numbers, by values written in other
parameters, to any depth, or by frequency dependent values. The
components which depend on the selected parameters, and the exact
derivative of each, come from the expressions the values and the
definitions were written as, by the chain rule, in one sparse enumeration
of the nonzero derivatives ([`designdirections`](@ref)), so a derived
value like `Capacitor(C/4)` contributes its factor of one quarter without
being declared, as does `Capacitor(:Cd)` with `:Cd => C/4`; a
[`ScatteringParameters`](@ref) block
depends on a parameter through the analytic derivative its `derivatives`
states ([`designblockjacobian`](@ref)). A port's reference impedance
normalizes its waves, so a parameter it depends on moves that
normalization too, through the termination the port owns or, for a port
which owns none, through the port itself. The solve runs once, with the
adjoint sensitivities of [`hbsolve`](@ref) carrying the exact direction
`dv_k/dp_j` of each dependent component into the contraction
(`sensitivityoperatingpoint = true`, so the shift of the pump operating
point is included).

Carrying the direction, rather than scaling a relative derivative, is
what makes the result exact for complex component values: a parameter
which rotates a value in the complex plane, a loss tangent, say, has
`dv/dp` not parallel to `v`, which no single relative derivative can
represent. All the components a parameter touches merge into one
contraction, so a design variable shared across a long line costs one
contraction rather than one per cell.

`parameters` are the names, or the definition keys, of the parameters to
differentiate with respect to; by default every parameter the definitions
give a number and every parameter a block states a derivative for, sorted
by name. Returns
`(out, dSdp)`: the full [`hbsolve`](@ref) output, and a keyed array
`dS/dp` with the axes of the scattering sensitivity and a `parameter`
axis in place of the `component` axis; with `keyedarrays = false`, the
plain array of the solve's scattering sensitivity, `(output, input,
parameter, frequency)`, the port and the mode of each side flattened as in
the plain `S`, its component axis the parameters. Additional keyword
arguments are forwarded to `hbsolve`, as `keyedarrays` is.

# Extended help

A gradient based optimizer wants a closure from a parameter vector to a
value and a derivative, evaluated at the same point in sequence, so memoize
one solve for both. For the gain in dB,
`G_k = 20*log10(abs(S[out, in, k]))`, the chain rule is
`dG/dp = (20/log(10))*real(conj(S)*dSdp)/abs2(S)`:

```julia
circuit = Circuit([(:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(:Cc)),
    (:jj, 2, 0, JosephsonJunction(:Lj)), (:cj, 2, 0, Capacitor(1000e-15))])
mutable struct Objective; lastp; lastr; end
const obj = Objective(nothing, nothing)
function solveat(pvec)
    if obj.lastp != pvec
        # the parameters in the order of pvec, not sorted by name
        obj.lastr = designsensitivities(circuit, Dict(:Lj => pvec[1], :Cc => pvec[2]),
            ws, wp, sources, (2,), (8,); parameters = (:Lj, :Cc))
        obj.lastp = copy(pvec)
    end
    return obj.lastr
end
value(pvec) = [20*log10(abs(s))
    for s in solveat(pvec).out.linearized.S((0,),1,(0,),1,:)]
function jacobian(pvec)
    r = solveat(pvec)
    S = r.out.linearized.S((0,),1,(0,),1,:)
    # the parameter axis precedes the frequency axis of dSdp
    d = permutedims(r.dSdp((0,),1,(0,),1,:,:))
    return (20/log(10)).*real.(conj.(S).*d)./abs2.(S)
end
```
"""
function designsensitivities(circuit::CompilableCircuit, circuitdefs::AbstractDict,
        ws, wp, sources, Nmodulationharmonics, Npumpharmonics;
        parameters = nothing, keyedarrays::Bool = true, kwargs...)
    psc = compile(circuit)
    definitions = definitiontable(circuitdefs)
    names = designparameters(parameters, definitions, psc)
    d = designdirections(psc, definitions, names)
    blockpairs = designblockjacobian(psc, names)
    # One pair per (component, parameter) dependence, carrying the exact
    # direction of the component value under the parameter as the rescale
    # alpha = (dv_k/dp_j)/v_k. The solver applies alpha to the component's
    # stamp before the negative frequency conjugation, which is what makes
    # this exact for complex component values.
    pairs = Tuple{String,Int,Complex{Float64}}[]
    for (k, j, x) in d.entries
        name = psc.componentnames[d.components[k]]
        iszero(d.values[k]) && throw(ArgumentError(
            lazy"the component $(name) depends on the parameter $(names[j]) but has the value zero at this point, so its stamp carries no direction to rescale."))
        push!(pairs, (name, j, x/d.values[k]))
    end
    isempty(pairs) && isempty(blockpairs) && throw(ArgumentError(
        "no component value depends on the selected design parameters."))
    out = hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        psc, definitions;
        sensitivitypairs = pairs,
        sensitivityblockpairs = blockpairs,
        nsensitivityparameters = length(names),
        sensitivitylabels = String.(names),
        returnSsensitivity = true,
        sensitivityoperatingpoint = true, keyedarrays = keyedarrays,
        kwargs...)
    # Ssensitivity already carries one slot per design parameter: plain, it
    # is the derivative as it stands, and keyed it is rewrapped with the
    # documented axis names
    keyedarrays || return (out = out, dSdp = out.linearized.Ssensitivity)
    Ss = Array(out.linearized.Ssensitivity)
    dSdpout = AxisKeys.KeyedArray(ComplexF64.(Ss),
        outputmode = out.linearized.modes,
        outputport = collect(out.linearized.portnumbers),
        inputmode = out.linearized.modes,
        inputport = collect(out.linearized.portnumbers),
        parameter = names,
        freqindex = 1:size(Ss, 6))
    return (out = out, dSdp = dSdpout)
end
