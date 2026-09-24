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

The derivative of a component value with respect to the design parameter
`name`, at `definitions`, as a number: zero for a number, one for the
parameter itself written as a symbol or a string, and for a
[`CircuitValue`](@ref) expression the derivative of the expression
evaluated at the definitions. A value which does not resolve to a number
is refused, as is a frequency dependent one. The Symbolics extension adds
the method for a `Num`.
"""
designderivative(::Number, ::Symbol, definitions) = zero(ComplexF64)
designderivative(v::Union{Symbol,AbstractString}, name::Symbol, definitions) =
    Symbol(v) === name ? one(ComplexF64) : zero(ComplexF64)
function designderivative(v::CircuitValue, name::Symbol, definitions)
    d = valuetonumber(CircuitValues.derivative(v, name), definitions)
    d isa Number || throw(ArgumentError(lazy"the derivative of the value $(v) with respect to $(name) is $(d) at these definitions, which is not a number; design sensitivities need every parameter defined."))
    return ComplexF64(d)
end
designderivative(v, ::Symbol, definitions) =
    throw(ArgumentError(lazy"the value $(v) cannot be differentiated with respect to a design parameter; write component values as numbers, parameters or expressions in parameters (frequency dependent values are not supported)."))

# the names of the design parameters: the ones given, as names or as
# definition keys, or by default every defined parameter and every
# parameter a scattering block states a derivative for, in sorted order
function designparameters(parameters, definitions, psc::CompiledCircuit)
    if isnothing(parameters)
        names = Symbol[]
        for k in keys(definitions)
            n = definitionname(k)
            isnothing(n) || push!(names, n)
        end
        for b in psc.scatteringblocks
            append!(names, keys(b.definition.derivatives))
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
    designjacobian(circuit, circuitdefs; parameters = nothing)

The Jacobian of the component values of a typed [`Circuit`](@ref), or of
its compiled circuit, with respect to its design parameters: the values
are written in terms of parameters (symbols, or the parameters of
[`@params`](@ref)) and `circuitdefs` gives every parameter a number; the
entry `(k, j)` is `d(value of component k)/d(parameter j)` at the
definitions, exact (see [`designderivative`](@ref)).

Returns `(names, values, J)`: the component names in compiled order, their
values at the definitions, and the complex `length(names)` by
`length(parameters)` Jacobian. `parameters` are names or definition keys;
by default every defined parameter and every parameter a
[`ScatteringParameters`](@ref) block states a derivative for, sorted by
name. Analysis ports are not among the components: a port's slot holds
its reference impedance, which appears under the termination generated for
it. A scattering block has no scalar value to differentiate; its dependence
is carried separately by [`designblockjacobian`](@ref).
"""
function designjacobian(circuit::CompilableCircuit, circuitdefs::AbstractDict;
        parameters = nothing)
    psc = compile(circuit)
    definitions = definitiontable(circuitdefs)
    names = designparameters(parameters, definitions, psc)
    vvn = componentvaluestonumber(psc.componentvalues, definitions)
    keep = [i for i in eachindex(vvn) if psc.componenttypes[i] !== :P]
    for i in keep
        vvn[i] isa Number || throw(ArgumentError(lazy"the component $(psc.componentnames[i]) has the value $(vvn[i]) at these definitions, which is not a number; design sensitivities need every parameter defined and no frequency dependent value."))
    end
    v0 = ComplexF64[vvn[i] for i in keep]
    J = zeros(ComplexF64, length(keep), length(names))
    # the definitions an expression substitutes are normalized once for the
    # whole Jacobian, as `componentvaluestonumber` normalizes them once for
    # the whole value table; every other representation of a value takes
    # the definitions as they were given.
    normalized = normalizedefinitions(definitions)
    for (j, name) in enumerate(names), (k, i) in enumerate(keep)
        value = psc.componentvalues[i]
        J[k, j] = value isa CircuitValue ?
            designderivative(value, name, normalized) :
            designderivative(value, name, definitions)
    end
    return String[psc.componentnames[i] for i in keep], v0, J
end

"""
    designblockjacobian(circuit, parameters)

The scattering block dependence of a circuit on its design parameters: for
each [`ScatteringParameters`](@ref) block of the compiled circuit and each
parameter its `derivatives` names, a derivative block whose scattering
matrix is the block's stated `dS/dp`. A block which states no derivative
for a parameter does not depend on it. Returns a vector of
`(blockpath, parameterindex, derivativeblock)`.
"""
function designblockjacobian(circuit::CompilableCircuit, parameters)
    psc = compile(circuit)
    out = Tuple{String,Int,Any}[]
    for b in psc.scatteringblocks, (j, name) in enumerate(parameters)
        d = b.definition.derivatives
        haskey(d, name) || continue
        # a derivative is not a passive scattering matrix and is not
        # checked as one: the positional constructor
        push!(out, (b.path, j, ScatteringParameters(d[name], b.definition.nports,
            b.definition.zref, b.definition.grounded, b.definition.noise,
            b.definition.negative_frequency)))
    end
    return out
end

"""
    designsensitivities(circuit, circuitdefs, ws, wp, sources,
        Nmodulationharmonics, Npumpharmonics; parameters = nothing,
        kwargs...)

The derivative of the scattering parameters with respect to the design
parameters of a typed [`Circuit`](@ref) whose component values are
written in terms of them, by the chain rule through the component
sensitivities:

    dS/dp_j = sum_k (dS/dv_k) (dv_k/dp_j).

`circuitdefs` gives every parameter a number. The components which depend
on the selected parameters, and the exact derivative of each, come from
the expressions the values were written as ([`designjacobian`](@ref)), so
a derived value like `Capacitor(C/4)` contributes its factor of one
quarter without being declared; a [`ScatteringParameters`](@ref) block
depends on a parameter through the analytic derivative its `derivatives`
states ([`designblockjacobian`](@ref)). The solve runs once, with the
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
differentiate with respect to; by default every defined parameter and
every parameter a block states a derivative for, sorted by name. Returns
`(out, dSdp)`: the full [`hbsolve`](@ref) output, and a keyed array
`dS/dp` with the axes of the scattering sensitivity and a `parameter`
axis in place of the `component` axis. Additional keyword arguments are
forwarded to `hbsolve`.

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
        obj.lastr = designsensitivities(circuit, Dict(:Lj => pvec[1], :Cc => pvec[2]),
            ws, wp, sources, (2,), (8,))
        obj.lastp = copy(pvec)
    end
    return obj.lastr
end
value(pvec) = [20*log10(abs(s))
    for s in solveat(pvec).out.linearized.S((0,),2,(0,),1,:)]
function jacobian(pvec)
    r = solveat(pvec)
    S = r.out.linearized.S((0,),2,(0,),1,:)
    # the parameter axis precedes the frequency axis of dSdp
    d = permutedims(r.dSdp((0,),2,(0,),1,:,:))
    return (20/log(10)).*real.(conj.(S).*d)./abs2.(S)
end
```
"""
function designsensitivities(circuit::CompilableCircuit, circuitdefs::AbstractDict,
        ws, wp, sources, Nmodulationharmonics, Npumpharmonics;
        parameters = nothing, kwargs...)
    psc = compile(circuit)
    definitions = definitiontable(circuitdefs)
    names = designparameters(parameters, definitions, psc)
    componentnames, v0, J = designjacobian(psc, definitions; parameters = names)
    blockpairs = designblockjacobian(psc, names)
    # One pair per (component, parameter) dependence, carrying the exact
    # direction of the component value under the parameter as the rescale
    # alpha = (dv_k/dp_j)/v_k. The solver applies alpha to the component's
    # stamp before the negative frequency conjugation, which is what makes
    # this exact for complex component values.
    pairs = Tuple{String,Int,Complex{Float64}}[]
    for j in eachindex(names), k in eachindex(componentnames)
        iszero(J[k, j]) && continue
        iszero(v0[k]) && throw(ArgumentError(
            lazy"the component $(componentnames[k]) depends on the parameter $(names[j]) but has the value zero at this point, so its stamp carries no direction to rescale."))
        push!(pairs, (componentnames[k], j, J[k, j]/v0[k]))
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
        sensitivityoperatingpoint = true, kwargs...)
    # Ssensitivity already carries one slot per design parameter; rewrap
    # with the documented axis names.
    Ss = Array(out.linearized.Ssensitivity)
    ndims(Ss) == 6 || throw(DimensionMismatch(
        lazy"unexpected Ssensitivity layout with $(ndims(Ss)) axes."))
    dSdpout = AxisKeys.KeyedArray(ComplexF64.(Ss),
        outputmode = out.linearized.modes,
        outputport = collect(out.linearized.portnumbers),
        inputmode = out.linearized.modes,
        inputport = collect(out.linearized.portnumbers),
        parameter = names,
        freqindex = 1:size(Ss, 6))
    return (out = out, dSdp = dSdpout)
end
