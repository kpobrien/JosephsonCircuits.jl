# =========================================================================
# Design-parameter sensitivities through circuit builders.
#
# A parameterized circuit is an ordinary Julia function (a builder) from
# real design parameters to a numeric circuit, so the derivative of the
# component values with respect to the design parameters is the Jacobian
# of ordinary Julia code, and the derivative of the scattering parameters
# follows by the chain rule through the existing component-sensitivity
# machinery:
#
#     dS/dp_j = sum_k (dS/dv_k) (dv_k/dp_j).
#
# The factorization is deliberate: reverse mode where it is expensive
# (through the solve, once, via the adjoint component sensitivities) and
# forward differences where they are cheap (through the builder, twice per
# parameter, microseconds each).
# =========================================================================

"""
    designparse(circuit)

Parse a builder output -- a legacy tuple netlist or a typed
[`Circuit`](@ref) -- into `(names, values)` in the solver's sorted
component order, so the names returned to the caller are the names the
solver uses. A [`ScatteringParameters`](@ref) port appears with its
`ScatteringStamp` as the value; every other value must be numeric.

Analysis ports are not among them. A port occupies a slot in the flat table
whose value is the port *number*, which is a label and not a quantity, so
differencing it is meaningless: at best it is an axis which is always zero,
and at worst a builder whose port numbering moved with a parameter would
report a derivative for it and the sensitivity would try to rescale a stamp
which does not exist.

A port's reference impedance is a design quantity, and it appears here under
the component which realizes it: the resistor a legacy netlist placed across
the port, or the environment the lowering generated for a matched one, named
after the port. That component is the one whose value the matrices read, so
it is the one a sensitivity has to name.
"""
function designparse(circuit)
    psc = compile(circuit)
    keep = [i for i in eachindex(psc.componentvalues)
            if psc.componenttypes[i] !== :P]
    names = String[psc.componentnames[i] for i in keep]
    vals = Vector{Any}(undef, length(keep))
    for (k, i) in enumerate(keep)
        v = psc.componentvalues[i]
        if v isa Number
            vals[k] = v
        else
            throw(ArgumentError(lazy"the component $(psc.componentnames[i]) has the non-numeric value $(v); design sensitivities require a fully numeric builder output (frequency dependent and symbolic values are not supported)."))
        end
    end
    # a scattering block is one entry, named by its instance path and
    # carrying its definition. It used to be one entry per port, each
    # carrying the whole definition, so a two port block appeared twice.
    for b in psc.scatteringblocks
        push!(names, b.path)
        push!(vals, b.definition)
    end
    return names, vals
end

"""
    designjacobian(builder, p::NamedTuple; parameters = keys(p),
        delta = 1e-6)

The Jacobian of the component values with respect to the design
parameters: `builder(; p...)` returns a numeric circuit, and the entry
`(k, j)` is `d(value of component k)/d(parameter j)`, by central finite
differences with relative step `delta`.

Returns `(names, values, J)`: the component names in circuit order, their
values at `p`, and the complex `length(names)` by `length(parameters)`
Jacobian. The builder must keep the circuit topology fixed as the
parameters vary; a parameter-dependent component *list* has no derivative
and is reported as an error.
"""
function designjacobian(builder, p::NamedTuple;
        parameters = keys(p), delta = 1e-6)
    names, raw0 = designparse(builder(; p...))
    v0 = ComplexF64[v isa Number ? ComplexF64(v) : 0 for v in raw0]
    J = zeros(ComplexF64, length(names), length(parameters))
    # a ScatteringParameters has no scalar value to difference; its dependence
    # is a derivative of a whole S(w), detected and carried separately by
    # `designblockjacobian`. Blocks appear here as zero rows.
    numeric = [v isa Number for v in raw0]
    for (j, q) in enumerate(parameters)
        x = float(real(getproperty(p, q)))
        # a RELATIVE step: physical parameters span many decades (a
        # critical current is ~1e-7), so any absolute floor can exceed the
        # parameter itself and step the difference through zero
        h = iszero(x) ? delta : delta*abs(x)
        namesp, rawp = designparse(
            builder(; merge(p, NamedTuple{(q,)}((x + h,)))...))
        namesm, rawm = designparse(
            builder(; merge(p, NamedTuple{(q,)}((x - h,)))...))
        namesp == names && namesm == names ||
            throw(ArgumentError(lazy"the builder changed the circuit topology when the parameter $(q) was perturbed; design sensitivities require a fixed component list."))
        for k in eachindex(names)
            numeric[k] || continue
            J[k, j] = (ComplexF64(rawp[k]) - ComplexF64(rawm[k]))/(2h)
        end
    end
    return names, v0, J
end

"""
    designblockjacobian(builder, p::NamedTuple; parameters = keys(p),
        delta = 1e-6)

The scattering block dependence of a circuit builder: for each design
parameter and each [`ScatteringParameters`](@ref) whose definition the builder
rebuilds when that parameter is perturbed, a derivative block whose
scattering matrix is `dS/dp`: the block's own analytic derivative when its
`derivatives` named tuple has an entry for the parameter, otherwise
central finite differences through the block's provider, evaluated lazily
at whatever frequencies the solver requests.

A block is treated as parameter independent when the perturbed builder
calls return *the identical block object* (`===`): hoist a fixed block --
measured Touchstone data, say -- out of the builder and it costs nothing,
while a block constructed inside the builder is differentiated. Returns a
vector of `(firstportname, parameterindex, derivativeblock)`.
"""
function designblockjacobian(builder, p::NamedTuple;
        parameters = keys(p), delta = 1e-6)
    names0, raw0 = designparse(builder(; p...))
    # the scattering blocks, one entry each
    firstport = [i for i in eachindex(raw0)
                 if raw0[i] isa ScatteringParameters]
    out = Tuple{String,Int,Any}[]
    isempty(firstport) && return out
    for (j, q) in enumerate(parameters)
        x = float(real(getproperty(p, q)))
        h = iszero(x) ? delta : delta*abs(x)
        _, rawp = designparse(
            builder(; merge(p, NamedTuple{(q,)}((x + h,)))...))
        _, rawm = designparse(
            builder(; merge(p, NamedTuple{(q,)}((x - h,)))...))
        for i in firstport
            bp = rawp[i]
            bm = rawm[i]
            b0 = raw0[i]
            (bp === b0 && bm === b0) && continue
            bp.nports == b0.nports && bm.nports == b0.nports ||
                throw(ArgumentError(lazy"the builder changed the port count of the scattering block at $(names0[i]) when the parameter $(q) was perturbed."))
            dblock = if haskey(b0.derivatives, q)
                # the analytic escape hatch: the block carries its own
                # dS/dtheta for this parameter
                ScatteringParameters(b0.derivatives[q], b0.nports, b0.zref,
                    b0.grounded, b0.noise, b0.negative_frequency)
            else
                fdscatteringderivative(b0, bp, bm, 2h)
            end
            push!(out, (names0[i], j, dblock))
        end
    end
    return out
end

"""
    fdscatteringderivative(b0, bp, bm, twoh)

A [`ScatteringParameters`](@ref) whose scattering matrix is the central finite
difference `(S_plus(w) - S_minus(w))/twoh` of the two perturbed blocks,
carrying the reference impedances and the negative frequency convention of
the unperturbed block, so the solver treats the derivative with exactly
the conventions of the value. Constructed directly rather than through the
public constructor, because a derivative is not a passive scattering
matrix and must not be checked as one.
"""
function fdscatteringderivative(b0, bp, bm, twoh)
    n = b0.nports
    f = function (w)
        Sp = Array{Complex{Float64},3}(undef, n, n, 1)
        Sm = Array{Complex{Float64},3}(undef, n, n, 1)
        evaluateprovider!(Sp, bp.provider, [w])
        evaluateprovider!(Sm, bm.provider, [w])
        return (Sp[:, :, 1] .- Sm[:, :, 1])./twoh
    end
    return ScatteringParameters(CallableMatrixProvider(f, n, :matrix), n,
        b0.zref, b0.grounded, b0.noise, b0.negative_frequency)
end

"""
    designsensitivities(builder, p::NamedTuple, ws, wp, sources,
        Nmodulationharmonics, Npumpharmonics; parameters = keys(p),
        delta = 1e-6, circuitdefs = Dict{Any,Any}(), kwargs...)

The derivative of the scattering parameters with respect to the design
parameters of a circuit builder, by the chain rule through the component
sensitivities:

    dS/dp_j = sum_k (dS/dv_k) (dv_k/dp_j).

`builder(; p...)` must return a fully numeric circuit. The components
which depend on the selected parameters are found automatically from the
builder Jacobian, so a derived value like `L2 = L/2` contributes its
factor of one half without being declared. The solve runs once, with the
adjoint sensitivities of [`hbsolve`](@ref) carrying the exact direction
`dv_k/dp_j` of each dependent component into the contraction
(`sensitivityoperatingpoint = true`, so the shift of the pump operating
point is included); the builder is evaluated twice per parameter, which is
negligible next to the solve.

Carrying the direction, rather than scaling a relative derivative, is
what makes the result exact for complex component values: a parameter
which rotates a value in the complex plane -- a loss tangent, say -- has
`dv/dp` not parallel to `v`, which no single relative derivative can
represent. All the components a parameter touches merge into one
contraction, so a design variable shared across a long line costs one
contraction rather than one per cell.

Returns `(out, dSdp)`: the full [`hbsolve`](@ref) output, and a keyed
array `dS/dp` with the axes of the scattering sensitivity and a
`parameter` axis in place of the `component` axis. Additional keyword
arguments are forwarded to `hbsolve`.

# Extended help
A gradient based optimizer wants a closure from a parameter vector to a
value and a derivative, evaluated at the same point in sequence, so memoize
one solve for both. For the gain in dB,
`G_k = 20*log10(abs(S[out, in, k]))`, the chain rule is
`dG/dp = (20/log(10))*real(conj(S)*dSdp)/abs2(S)`:

```julia
mutable struct Objective; lastp; lastr; end
const obj = Objective(nothing, nothing)
function solveat(pvec)
    p = (Lj = pvec[1], Cc = pvec[2])
    if obj.lastp != pvec
        obj.lastr = designsensitivities(make, p, ws, wp, sources, (2,), (8,))
        obj.lastp = copy(pvec)
    end
    return obj.lastr
end
value(pvec) = [20*log10(abs(s))
    for s in solveat(pvec).out.linearized.S((0,),2,(0,),1,:)]
function jacobian(pvec)
    r = solveat(pvec)
    S = r.out.linearized.S((0,),2,(0,),1,:)
    d = r.dSdp((0,),2,(0,),1,:,:)
    return (20/log(10)).*real.(conj.(S).*d)./abs2.(S)
end
```
"""
function designsensitivities(builder, p::NamedTuple, ws, wp, sources,
        Nmodulationharmonics, Npumpharmonics;
        parameters = keys(p), delta = 1e-6,
        circuitdefs = Dict{Any,Any}(), kwargs...)

    names, v0, J = designjacobian(builder, p;
        parameters = parameters, delta = delta)
    blockpairs = designblockjacobian(builder, p;
        parameters = parameters, delta = delta)

    # One pair per (component, parameter) dependence, carrying the exact
    # direction of the component value under the parameter as the rescale
    # alpha = (dv_k/dp_j)/v_k. The solver applies alpha to the component's
    # stamp *before* the negative frequency conjugation, which is what
    # makes this exact for complex component values: a parameter which
    # rotates a value in the complex plane (a loss tangent, say) has a
    # direction which is not parallel to the value, and no single relative
    # derivative can represent it.
    pairs = Tuple{String,Int,Complex{Float64}}[]
    for j in eachindex(parameters), k in eachindex(names)
        iszero(J[k, j]) && continue
        iszero(v0[k]) && throw(ArgumentError(
            lazy"the component $(names[k]) depends on the parameter $(parameters[j]) but has the value zero at this point, so its stamp carries no direction to rescale."))
        push!(pairs, (names[k], j, J[k, j]/v0[k]))
    end
    isempty(pairs) && isempty(blockpairs) && throw(ArgumentError(
        "no component value depends on the selected design parameters."))

    out = hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics,
        builder(; p...), circuitdefs;
        sensitivitypairs = pairs,
        sensitivityblockpairs = blockpairs,
        nsensitivityparameters = length(parameters),
        sensitivitylabels = String.(Symbol.(collect(parameters))),
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
        parameter = collect(Symbol.(parameters)),
        freqindex = 1:size(Ss, 6))
    return (out = out, dSdp = dSdpout)
end
