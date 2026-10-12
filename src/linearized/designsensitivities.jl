# =========================================================================
# Design parameter sensitivities of a typed circuit.
#
# A circuit's component values are written as names (`Capacitor(:Cc)`),
# and a Julia function of the design parameters gives the value of each
# name: `values(p)`, a NamedTuple. The derivative of the values with
# respect to the parameters is that function's Jacobian, by ForwardDiff
# (an extension, loaded with `using ForwardDiff`) or as the caller states
# it, and the derivative of the scattering parameters
# follows by the chain rule through the component sensitivities:
#
#     dS/dp_j = sum_k (dS/dv_k) (dv_k/dp_j),
#
# reverse mode through the solve, once, via the adjoint component
# sensitivities, and the exact direction of every dependent value carried
# into the contraction.
# =========================================================================

# The values of the circuit's names at the design parameters `p`, as the
# function `f` gives them, a NamedTuple of numbers. A frequency dependent
# value is no value a design parameter may move: it is written into the
# circuit as the component's value instead.
function designvalues(f, p::NamedTuple)
    all(x -> x isa Real && isfinite(x), p) || throw(ArgumentError(
        lazy"the design parameters $(p) must be finite real numbers."))
    v = f(p)
    v isa NamedTuple || throw(ArgumentError(lazy"the values function returns a $(typeof(v)); it returns a NamedTuple of the values of the circuit's names."))
    for (name, x) in pairs(v)
        x isa Number || throw(ArgumentError(lazy"the values function gives $(name) the value $(x); a value a design parameter may move is a number, and a frequency dependent value, which none may move, is written into the circuit as the component's value."))
    end
    return v
end

# The derivatives of the values `v = f(p)` with respect to the design
# parameters, complex, a row per name of `v` and a column per parameter of
# `p`, in their orders: `jacobian(p)` when it is given, and otherwise
# ForwardDiff's (`automaticjacobian`).
function designvaluejacobian(f, p::NamedTuple, v::NamedTuple{V},
        jacobian) where {V}
    J = isnothing(jacobian) ? automaticjacobian(f, p, V) : jacobian(p)
    (J isa AbstractMatrix && size(J) == (length(v), length(p))) ||
        throw(DimensionMismatch(lazy"the Jacobian of the values has a row per value, $(length(v)), and a column per design parameter, $(length(p)); `jacobian(p)` returned a $(summary(J))."))
    all(isfinite, J) || throw(ArgumentError(
        "the Jacobian of the values is not finite."))
    return J
end

# ForwardDiff's Jacobian of the real and imaginary parts of the values of
# `f`, the names `names`, as real functions of the real parameters `p`,
# which the ForwardDiff extension computes: ForwardDiff is loaded for it.
function automaticjacobian(f, p::NamedTuple, names)
    ext = Base.get_extension(@__MODULE__, :JosephsonCircuitsForwardDiffExt)
    isnothing(ext) && throw(ArgumentError("the derivatives of the values function are ForwardDiff's, which is loaded for them (`using ForwardDiff`), or are given as `jacobian = p -> J`, a row per value and a column per design parameter."))
    return ext.valuejacobian(f, p, names)
end

"""
    designdirections(psc::CompiledCircuit, v::NamedTuple, J::AbstractMatrix)

The nonzero derivatives of the component values of `psc` with respect to
the design parameters, from the values `v` of the names its components
are written in and their Jacobian `J`, a row per name of `v` and a column
per parameter: a component written as a name takes that name's row, and
one written as a number depends on no parameter. A component whose value
is frequency dependent is not among them, and neither is a port which owns
a termination, whose reference impedance is the value of that
termination, which is. Returns `(components, values, entries)`: the flat
indices of the components, in compiled order, their values, and a triplet
`(k, j, dv/dp)` per nonzero derivative of the `k`th of them with respect
to parameter `j`, by parameter and then by component, so that the cost
grows with the dependences rather than with the components times the
parameters.
"""
function designdirections(psc::CompiledCircuit, v::NamedTuple,
        J::AbstractMatrix)
    row = Dict(name => r for (r, name) in enumerate(keys(v)))
    terminated = Set(p.component for p in psc.ports if !iszero(p.environment))
    components = Int[]
    values = ComplexF64[]
    entries = Tuple{Int,Int,ComplexF64}[]
    for (i, written) in enumerate(psc.componentvalues)
        i in terminated && continue
        name = definitionname(written)
        if isnothing(name)
            written isa Number && !checkissymbolic(written) || continue
            push!(components, i)
            push!(values, ComplexF64(written))
            continue
        end
        r = get(row, name, 0)
        iszero(r) && throw(ArgumentError(lazy"the component $(psc.componentnames[i]) is written as $(written), a name the values function does not give."))
        push!(components, i)
        push!(values, ComplexF64(v[r]))
        for j in axes(J, 2)
            x = J[r, j]
            iszero(x) || push!(entries, (length(components), j, x))
        end
    end
    sort!(entries; by = e -> (e[2], e[1]))
    return (; components, values, entries)
end

# The drive of a design sensitivity, its tone frequencies `wp` and the
# currents of its `sources`, may be written as names of the values `v`,
# which the solve takes; the derivatives carry the components and the
# blocks alone, so a drive which a design parameter moves, through its row
# of the Jacobian `J`, is refused rather than left out of the derivative.
function checkdesigndrive(v::NamedTuple{V}, J::AbstractMatrix, wp, sources,
        parameters) where {V}
    # the row of the name a drive is written as when a parameter moves it,
    # else zero
    moving(written) = (r = something(findfirst(==(definitionname(written)), V), 0);
        !iszero(r) && any(!iszero, view(J, r, :)) ? r : 0)
    for (i, w) in enumerate(wp)
        r = moving(w)
        iszero(r) || throwdesigndrive(lazy"the frequency of tone $(i)", w, J, r,
            parameters)
    end
    for source in sources
        r = moving(source[:current])
        iszero(r) || throwdesigndrive(
            lazy"the current of the source at port $(source[:port]) and mode $(source[:mode])",
            source[:current], J, r, parameters)
    end
    return nothing
end

# the refusal of a drive which the parameters move, built only when it is
# thrown
@noinline function throwdesigndrive(what, written, J, r, parameters)
    moved = join((parameters[j] for j in axes(J, 2) if !iszero(J[r, j])), ", ")
    throw(ArgumentError(lazy"$(what) is written as $(written), which the design parameters $(moved) move; the derivatives carry the motion of the component values and not of the drive, so write it as a number or as a name no parameter moves."))
end

"""
    designjacobian(circuit, values, p::NamedTuple; jacobian = nothing)

The Jacobian of the component values of a typed [`Circuit`](@ref), or of
its compiled circuit, with respect to the design parameters `p`, a
NamedTuple of real numbers. `values(p)` returns the values of the names
the circuit's components are written in as a NamedTuple of numbers, the
definitions a solve of the circuit takes; a derived value is computed
there, in Julia. The derivatives of the values are `jacobian(p)` when it
is given, a matrix with a row per name of `values(p)` and a column per
parameter of `p`, in their orders, and otherwise those of ForwardDiff
through `values`, which is loaded for them (`using ForwardDiff`), whose
real and imaginary parts are differentiated as real functions of the real
parameters, so a value a parameter rotates in the complex plane is exact.
A frequency dependent value is no value of `values(p)`: no design
parameter may move it, and it is written into the circuit as the
component's value.

Returns `(names, values, J)`: the component names in compiled order, their
values, and the complex `length(names)` by `length(p)` Jacobian, dense,
materialized from the sparse enumeration of the nonzero derivatives
([`designdirections`](@ref)). A port's slot holds its reference
impedance: a port which owns a termination is not among the components,
its reference impedance being the value of that termination, which is; a
port which owns none (`termination = nothing`) is, and a parameter moves
only the normalization of its waves through it. A component whose value is
frequency dependent is not among them either. A scattering block has no
scalar value to differentiate; its dependence is carried separately by
[`designblockjacobian`](@ref).
"""
function designjacobian(circuit::CompilableCircuit, values, p::NamedTuple;
        jacobian = nothing)
    psc = compile(circuit)
    v = designvalues(values, p)
    d = designdirections(psc, v, designvaluejacobian(values, p, v, jacobian))
    J = zeros(ComplexF64, length(d.components), length(p))
    for (k, j, x) in d.entries
        J[k, j] = x
    end
    return String[psc.componentnames[i] for i in d.components], d.values, J
end

# the analytic derivatives of a block's scattering matrix, as
# `name => provider` pairs: only a `ScatteringParameters` block states any
blockderivatives(d::ScatteringParameters) = d.derivatives
blockderivatives(d) = Pair{Symbol,Any}[]

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
    designsensitivities(circuit, values, p::NamedTuple, ws, wp, sources,
        Nmodulationharmonics, Npumpharmonics; jacobian = nothing,
        keyedarrays = true, kwargs...)

The derivative of the scattering parameters with respect to the design
parameters `p`, a NamedTuple of real numbers, of a typed
[`Circuit`](@ref) whose component values are written as names, which the
function `values` gives, `values(p)` a NamedTuple of numbers, by the chain
rule through the component sensitivities:

    dS/dp_j = sum_k (dS/dv_k) (dv_k/dp_j).

The derivatives of the values are those of [`designjacobian`](@ref), by
ForwardDiff through `values`, loaded for them (`using ForwardDiff`), or
from `jacobian(p)`, so a derived value
like `C/4` contributes its factor of one quarter without being declared;
a [`ScatteringParameters`](@ref) block depends on a parameter of `p`
through the analytic derivative its `derivatives` states for that
parameter's name ([`designblockjacobian`](@ref)). A port's reference
impedance normalizes its waves, so a parameter it depends on moves that
normalization too, through the termination the port owns or, for a port
which owns none, through the port itself. The solve runs once, at the
definitions `values(p)`, with the adjoint sensitivities of
[`hbsolve`](@ref) carrying the exact direction `dv_k/dp_j` of each
dependent component into the contraction (`sensitivityoperatingpoint =
true`, so the shift of the pump operating point is included). A tone
frequency or a source current may be written as a name `values(p)` gives,
which the solve takes, but the derivative carries the component values
alone: a drive which a design parameter moves is an `ArgumentError`
rather than a derivative which leaves its motion out.

Carrying the direction, rather than scaling a relative derivative, is
what makes the result exact for complex component values: a parameter
which rotates a value in the complex plane, a loss tangent, say, has
`dv/dp` not parallel to `v`, which no single relative derivative can
represent. All the components a parameter touches merge into one
contraction, so a design variable shared across a long line costs one
contraction rather than one per cell.

Returns `(out, dSdp)`: the full [`hbsolve`](@ref) output, and a keyed
array `dS/dp` with the axes of the scattering sensitivity and a
`parameter` axis, the names of `p`, in place of the `component` axis;
with `keyedarrays = false`, the plain array of the solve's scattering
sensitivity, `(output, input, parameter, frequency)`, the port and the
mode of each side flattened as in the plain `S`, its component axis the
parameters. Additional keyword arguments are forwarded to `hbsolve`, as
`keyedarrays` is.

# Extended help

A gradient based optimizer wants a closure from a parameter vector to a
value and a derivative, evaluated at the same point in sequence, so memoize
one solve for both. For the gain in dB,
`G_k = 20*log10(abs(S[out, in, k]))`, the chain rule is
`dG/dp = (20/log(10))*real(conj(S)*dSdp)/abs2(S)`:

```julia
circuit = Circuit([(:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(:Cc)),
    (:jj, 2, 0, JosephsonJunction(:Lj)), (:cj, 2, 0, Capacitor(1000e-15))])
values(p) = (Lj = p.Lj, Cc = p.Cc)
mutable struct Objective; lastp; lastr; end
const obj = Objective(nothing, nothing)
function solveat(pvec)
    if obj.lastp != pvec
        obj.lastr = designsensitivities(circuit, values,
            (Lj = pvec[1], Cc = pvec[2]), ws, wp, sources, (2,), (8,))
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
function designsensitivities(circuit::CompilableCircuit, values,
        p::NamedTuple, ws, wp, sources, Nmodulationharmonics, Npumpharmonics;
        jacobian = nothing, keyedarrays::Bool = true, kwargs...)
    psc = compile(circuit)
    v = designvalues(values, p)
    J = designvaluejacobian(values, p, v, jacobian)
    names = collect(keys(p))
    checkdesigndrive(v, J, wp, sources, names)
    d = designdirections(psc, v, J)
    blockpairs = designblockjacobian(psc, names)
    # One pair per (component, parameter) dependence, carrying the exact
    # direction of the component value under the parameter as the rescale
    # alpha = (dv_k/dp_j)/v_k. The solver applies alpha to the component's
    # stamp before the negative frequency conjugation, which is what makes
    # this exact for complex component values.
    componentpairs = Tuple{String,Int,Complex{Float64}}[]
    for (k, j, x) in d.entries
        name = psc.componentnames[d.components[k]]
        iszero(d.values[k]) && throw(ArgumentError(
            lazy"the component $(name) depends on the parameter $(names[j]) but has the value zero at this point, so its stamp carries no direction to rescale."))
        push!(componentpairs, (name, j, x/d.values[k]))
    end
    isempty(componentpairs) && isempty(blockpairs) && throw(ArgumentError(
        "no component value depends on the design parameters."))
    out = hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        psc, Dict{Any,Any}(pairs(v));
        sensitivitypairs = componentpairs,
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
