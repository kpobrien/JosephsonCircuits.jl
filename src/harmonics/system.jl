# The harmonic balance system at a point: the residual, the Jacobian-vector
# and Hessian-vector products through the transforms, the assembled real and
# complex Jacobians, and rebinding to new component values.

"""
    HBSystem

Everything needed to evaluate the harmonic balance nonlinear system and its
derivatives at a point: the linear term matrices, the Josephson junction
data, the frequency domain packing maps and Fourier transform plans, the
real representation layout, optional precomputed Jacobian assembly plans,
and preallocated workspaces.

Set the evaluation point with [`setpoint!`](@ref), then evaluate any of:

- [`residual!`](@ref): the residual `F(x) = B(sin.(A*x)) + K*x - b`,
- [`jacobianvectorproduct!`](@ref): the exact matrix-free Jacobian-vector
  product `J(x)*v = B(cos.(A*x) .* (A*v)) + K*v`,
- [`hessianvectorproduct!`](@ref): the exact matrix-free second directional
  derivative `H(x)[v, w] = B(-sin.(A*x) .* (A*v) .* (A*w))`,
- [`jacobian!`](@ref): the assembled complex (holomorphic) or exact real
  Jacobian via the precomputed plans,

where `A` is the linear map from the unknowns to the time domain branch
fluxes on the Josephson junctions, `B` the linear map from a time domain
signal back to the node vector, and `K` the frequency dependent linear
terms. Each entry point has a complex representation method and an
equivalent real representation method, dispatched on the element type of the
output (and direction) vectors, so the same object serves both the
[`QuasiNewton`](@ref) and [`Newton`](@ref) methods of [`hbnlsolve`](@ref) as
well as matrix-free solvers.

The time domain branch fluxes and the pointwise sine and cosine at the
current point are cached, so repeated products at the same point (eg. the
many Jacobian-vector products of a Krylov solve) cost only two Fourier
transforms and the linear term each.

The fields are intentionally loosely typed; all performance critical loops
are behind function barriers which specialize on the concrete argument
types. The workspaces of the residual and the matrix-free products are
parameterized on their array types rather than fixed to `Array`, so they can
live on whichever KernelAbstractions backend the system was built for.
"""
struct HBSystem{TR,TinvL,TG,TC,TWm,TK,Tb,TLjb,TLjbm,TLm,TIP,TFP,TRJ,TCJ,TNP,TVC,TVR,TAC,TAR,TAM,TAB}
    # linear term matrices and source vector (complex representation, scaled
    # by Lscale, conjugated for negative frequency modes with conjnegfreq!)
    Rbnm::TR
    invLnm::TinvL
    Gnm::TG
    Cnm::TC
    wmodesm::TWm
    wmodes2m::TWm
    # the frequency dependent linear terms collapsed into a single matrix,
    # K = invLnm + im*Gnm*wmodesm - Cnm*wmodes2m. The mode frequencies are
    # fixed for the lifetime of the system, so this is formed once here
    # rather than at every residual and Jacobian-vector product. The
    # individual matrices are retained because the Jacobian assembly plans
    # scatter them separately, with the current frequencies, at assembly
    # time.
    Knm::TK
    # the source vector, on the backend: the complex representation residual
    # subtracts it from a vector which lives there, so it cannot stay on the
    # host the way it was handed in. bnmr below is the same vector in the
    # real representation and moves with it.
    bnm::Tb
    # Josephson junction data
    Ljb::TLjb
    Ljbm::TLjbm
    Lscale::TLm
    # frequency domain packing and Fourier transform plans
    freqindexmap::Vector{Int}
    conjsourceindices::Vector{Int}
    conjtargetindices::Vector{Int}
    irfftplan::TIP
    rfftplan::TFP
    # real mode layout
    modelayout::ModeLayout{Int}
    # optional precomputed Jacobian assembly plans
    realjacobianplan::TRJ
    complexjacobianplan::TCJ
    # the precomputed, device-generic plan for the two linear maps which
    # surround the pointwise time domain nonlinearity, used by the real
    # representation entry points
    nonlineartermplan::TNP
    # the source vector in the equivalent real representation, formed once
    # so the residual subtracts it as a plain vector operation
    bnmr::TVR
    # the current point (set with setpoint!), in both representations
    x::TVC
    xr::TVR
    # cached time domain branch fluxes and pointwise nonlinearities at the
    # current point
    phitd::TAR
    sintd::TAR
    costd::TAR
    sincurrent::Base.RefValue{Bool}
    coscurrent::Base.RefValue{Bool}
    # the negative of the second derivative of the relation at the point,
    # which the Hessian and the derivative of the linearized system with
    # respect to the operating point are written in. The sinusoidal
    # relation has `-f'' = sin`, so that case reads `sintd` and this array
    # is empty; a polynomial one fills it.
    negsecondtd::TAR
    negsecondcurrent::Base.RefValue{Bool}
    # the third derivative of the relation at the point, `-cos` for the
    # Josephson one, which the trilinear form of the problem interface
    # takes; empty for the Josephson relation, whose caller reads `-costd`
    thirdtd::TAR
    thirdcurrent::Base.RefValue{Bool}
    # the current-phase relation of every junction. The empty table is
    # every junction sinusoidal, which is every circuit that does not ask
    # for another, and the evaluations then take the `sin` and `cos` they
    # always did. Held as a table rather than as `nothing` so that the type
    # of a system, and with it everything compiled for it, does not depend
    # on what the circuit holds.
    relations::JunctionRelations{TAM,TAB}
    # workspaces: `phimatrix` is the frequency domain scratch of every
    # residual and product
    phimatrix::TAC
    dirtd::TAR
    dirtd2::TAR
    worktd::TAR
    # the Fourier coefficients of cos(phi(t)) at the current point, which
    # the Jacobian assemblies read, in a buffer of their own so that no
    # product between an update and an assembly can overwrite them;
    # `cosfdcurrent` says whether they are those of the current point
    cosfd::TAC
    cosfdcurrent::Base.RefValue{Bool}
end

"""
    HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm, Ljb,
        Ljbm, Lscale, Nbranches, freqindexmap, conjsourceindices, conjtargetindices,
        phimatrix, phimatrixtd, irfftplan, rfftplan, modelayout,
        realjacobianplan, complexjacobianplan, backend = CPU();
        realbackward = true, relations = nothing)

Construct an [`HBSystem`](@ref) from the ingredients assembled by
[`hbnlsolve`](@ref), allocating the workspaces. `phimatrix` and
`phimatrixtd` are adopted as the frequency domain and one of the time domain
workspaces. `realjacobianplan` and `complexjacobianplan` may be `nothing`,
in which case the corresponding [`jacobian!`](@ref) method is unavailable.
`realbackward = false` skips building the real representation of the
linear term in the nonlinear term plan, which only a solve in the real
representation applies; see [`plannonlinearterm`](@ref).
`backend` is the KernelAbstractions backend on which the index mapped
kernels run and on which the plan and the workspaces of the residual and the
matrix-free products are allocated; `CPU()` is the default and the reference.
`phimatrix` and `phimatrixtd` are adopted as given, so pass them already on
the backend, and the time domain workspaces derived from them with `similar`
follow.
`relations` are the [`JunctionRelations`](@ref) of the junctions, or
`nothing` when every one of them is the sinusoidal Josephson relation,
which is the case the evaluations take as the plain `sin` and `cos` they
always did. They are moved to the backend and to the working precision
here, and the cache of the second derivative is allocated only when there
is one to hold.
"""
function HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
    Ljb, Ljbm, Lscale, Nbranches, freqindexmap, conjsourceindices, conjtargetindices,
    phimatrix, phimatrixtd, irfftplan, rfftplan, modelayout,
    realjacobianplan, complexjacobianplan, backend = CPU();
    realbackward::Bool = true, relations = nothing)

    n = size(Rbnm, 2)
    # the working precision of everything the residual and the matrix-free
    # products touch is the one of the frequency domain array handed in, so
    # allocating `phimatrix` as Complex{Float32} carries the whole evaluation
    # path to Float32. The host side linear term matrices are left as they
    # were: they are only read by the assembly plan builders, which run on
    # the host at plan time and write into a matrix of their own element
    # type; the assembly itself runs on the backend.
    TF = real(eltype(phimatrix))
    # collapse the frequency dependent linear terms into a single matrix. the
    # mode frequencies do not change for the lifetime of the system, so the
    # three sparse products and the two diagonal scalings are done once here
    # instead of at every evaluation.
    Knm = linearterm(invLnm, Gnm, Cnm, wmodesm, wmodes2m, TF)
    # the index maps of the two linear maps which surround the pointwise time
    # domain nonlinearity, built once here. with realbackward = false the real
    # form of the linear term is left out, for a caller which only applies the
    # complex representation entry points.
    nonlineartermplan = plannonlinearterm(Rbnm, Ljb, Lscale, Nbranches,
        freqindexmap, conjsourceindices, conjtargetindices, phimatrix, Knm,
        modelayout, backend; realbackward = realbackward)
    # both representations of the source vector move to the backend. the real
    # one already did; the complex one did not, so `residual!` on a complex
    # vector broadcast a device array against a host one. `complex_to_real`
    # below reads the host argument, which `tobackend` does not modify.
    return HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, Knm,
        tobackend(backend, convert(Vector{Complex{TF}}, bnm)),
        Ljb, Ljbm, Lscale, freqindexmap, conjsourceindices,
        conjtargetindices, irfftplan, rfftplan, modelayout,
        realjacobianplan, complexjacobianplan, nonlineartermplan,
        tobackend(backend, convert(Vector{TF},
            complex_to_real(bnm, modelayout.isreal))),
        tobackend(backend, zeros(Complex{TF}, n)),
        tobackend(backend, zeros(TF, modelayout.rdim)),
        similar(phimatrixtd), similar(phimatrixtd), similar(phimatrixtd),
        Ref(false), Ref(false),
        # the second derivative is cached only where it is not the sine
        isnothing(relations) ? similar(phimatrixtd, ntuple(_ -> 0,
            ndims(phimatrixtd))) : similar(phimatrixtd),
        Ref(false),
        isnothing(relations) ? similar(phimatrixtd, ntuple(_ -> 0,
            ndims(phimatrixtd))) : similar(phimatrixtd),
        Ref(false), torelations(relations, phimatrixtd),
        phimatrix, similar(phimatrixtd), similar(phimatrixtd), phimatrixtd,
        similar(phimatrix), Ref(false))
end

# the structure of two sparse matrices, compared exactly: a value which
# moved under a fixed structure is a refresh, and one which moved the
# structure is a new system
function samestructure(A::SparseMatrixCSC, B::SparseMatrixCSC)
    return size(A) == size(B) &&
        SparseArrays.getcolptr(A) == SparseArrays.getcolptr(B) &&
        rowvals(A) == rowvals(B)
end

"""
    rebind!(sys::HBSystem, invLnm, Gnm, Cnm, bnm, Ljb, Ljbm, Lscale;
        maps = nothing, realjacobianplan = sys.realjacobianplan,
        complexjacobianplan = sys.complexjacobianplan)

The same system at new component values: the linear term matrices, the
source, the junction inductances and the scale are replaced by the new
ones, in place where the arrays live and by a new struct sharing them where
a scalar does. The transforms, the index maps, the kernels and every
workspace stay, and nothing the size of the plan is allocated; with `maps`
given the only allocations are the two temporaries of the source vector.

The structure must not have moved: a sparse pattern which differs from the
one the system was built on is refused, because the plans are built on the
pattern and a value which moves it is a new circuit and not a new point.

The assembly plans of the Jacobians hold the gathered linear term and, in
their [`JunctionStructure`](@ref), `Lscale/Lj`; a plan the system holds is
refreshed to the new values here, and a plan passed as `realjacobianplan`
or `complexjacobianplan`, built for the new values by a solve which needs
a Jacobian the system did not hold, is installed instead.
"""
function rebind!(sys::HBSystem, invLnm, Gnm, Cnm, bnm, Ljb, Ljbm, Lscale;
        maps::Union{Nothing,ValueMaps} = nothing,
        realjacobianplan = sys.realjacobianplan,
        complexjacobianplan = sys.complexjacobianplan)
    for (old, new, what) in ((sys.invLnm, invLnm, "inverse inductance"),
            (sys.Gnm, Gnm, "conductance"), (sys.Cnm, Cnm, "capacitance"))
        samestructure(old, new) || throw(ArgumentError(
            lazy"the $(what) matrix changed its sparse structure between points, which a reused system cannot follow; build a new one."))
        copyto!(nonzeros(old), nonzeros(new))
    end
    (Ljb.nzind == sys.Ljb.nzind && Ljbm.nzind == sys.Ljbm.nzind) ||
        throw(ArgumentError("the junctions changed between points, which a reused system cannot follow; build a new one."))
    copyto!(sys.Ljb.nzval, Ljb.nzval)
    copyto!(sys.Ljbm.nzval, Ljbm.nzval)
    TF = real(eltype(sys.phimatrix))
    if isnothing(maps)
        Knm = linearterm(sys.invLnm, sys.Gnm, sys.Cnm, sys.wmodesm,
            sys.wmodes2m, TF)
        samestructure(Knm, sys.Knm) || throw(ArgumentError(
            "the linear term changed its sparse structure between points, which a reused system cannot follow; build a new one."))
        copyto!(nonzeros(sys.Knm), nonzeros(Knm))
        refreshvalues!(sys.nonlineartermplan, sys.Rbnm, sys.Ljb, Lscale,
            sys.Knm, sys.modelayout, sys.freqindexmap)
    else
        linearterm!(sys.Knm, maps, sys.invLnm, sys.Gnm, sys.Cnm, sys.wmodesm,
            sys.wmodes2m)
        refreshvalues!(sys.nonlineartermplan, maps, sys.Knm, sys.Ljb, Lscale)
    end
    # the plans the system keeps follow the values; every plan's junction
    # structure is refreshed, harmlessly twice when two plans share one
    if !isnothing(realjacobianplan) && realjacobianplan === sys.realjacobianplan
        refreshvalues!(realjacobianplan, sys.invLnm, sys.Gnm, sys.Cnm,
            sys.wmodesm, sys.wmodes2m, sys.Ljb, Lscale)
    end
    if !isnothing(complexjacobianplan) && complexjacobianplan === sys.complexjacobianplan
        refreshvalues!(complexjacobianplan, sys.invLnm, sys.Gnm, sys.Cnm,
            sys.wmodesm, sys.wmodes2m, sys.Ljb, Lscale)
    end
    copyto!(sys.bnm, convert(Vector{Complex{TF}}, bnm))
    copyto!(sys.bnmr, convert(Vector{TF},
        complex_to_real(bnm, sys.modelayout.isreal)))
    # the point is stale: it was set against the old values
    sys.sincurrent[] = false
    sys.coscurrent[] = false
    sys.negsecondcurrent[] = false
    sys.thirdcurrent[] = false
    sys.cosfdcurrent[] = false
    return HBSystem(sys.Rbnm, sys.invLnm, sys.Gnm, sys.Cnm, sys.wmodesm,
        sys.wmodes2m, sys.Knm, sys.bnm, sys.Ljb, sys.Ljbm, Lscale,
        sys.freqindexmap, sys.conjsourceindices, sys.conjtargetindices,
        sys.irfftplan, sys.rfftplan, sys.modelayout, realjacobianplan,
        complexjacobianplan, sys.nonlineartermplan, sys.bnmr, sys.x,
        sys.xr, sys.phitd, sys.sintd, sys.costd, sys.sincurrent,
        sys.coscurrent, sys.negsecondtd, sys.negsecondcurrent,
        sys.thirdtd, sys.thirdcurrent, sys.relations, sys.phimatrix, sys.dirtd, sys.dirtd2, sys.worktd,
        sys.cosfd, sys.cosfdcurrent)
end

"""
    valuemaps(sys::HBSystem)

The [`ValueMaps`](@ref) of a system, for rebinding it without rebuilding
the plan, or `nothing` when a conversion has no fixed map.
"""
valuemaps(sys::HBSystem) = valuemaps(sys.nonlineartermplan, sys.Knm,
    sys.invLnm, sys.Gnm, sys.Cnm, sys.Rbnm, sys.Ljb, sys.modelayout,
    sys.freqindexmap)

"""
    linearterm!(Knm, maps::ValueMaps, invLnm, Gnm, Cnm, wmodesm, wmodes2m)

[`linearterm`](@ref) in place, through the maps: the same three terms in the
same order, so the values are the ones the allocating form gives.
"""
function linearterm!(Knm::SparseMatrixCSC, maps::ValueMaps, invLnm, Gnm,
        Cnm, wmodesm::Diagonal, wmodes2m::Diagonal)
    knz = nonzeros(Knm)
    lnz, gnz, cnz = nonzeros(invLnm), nonzeros(Gnm), nonzeros(Cnm)
    wm, wm2 = wmodesm.diag, wmodes2m.diag
    ml, mg, mc = maps.ml, maps.mg, maps.mc
    T = eltype(knz)
    @inbounds for j in axes(Knm, 2), t in nzrange(Knm, j)
        v = zero(T)
        iszero(ml[t]) || (v += lnz[ml[t]])
        iszero(mg[t]) || (v += im*(gnz[mg[t]]*wm[j]))
        iszero(mc[t]) || (v -= cnz[mc[t]]*wm2[j])
        knz[t] = v
    end
    return Knm
end

"""
    linearterm(invLnm, Gnm, Cnm, wmodesm::Diagonal, wmodes2m::Diagonal,
        T = Float64)

Collapse the frequency dependent linear terms of the harmonic balance system
into the single sparse matrix of element type `Complex{T}`

    K = invLnm + im*Gnm*wmodesm - Cnm*wmodes2m

so that applying them is one gather over the entries of an output row rather
than three sparse matrix-vector products and two diagonal scalings through
two temporaries. The mode frequency diagonals are fixed for the lifetime of
an [`HBSystem`](@ref), so this is formed once when the system is constructed
and handed to [`plannonlinearterm`](@ref), which stores it transposed in both
representations.
"""
function linearterm(invLnm, Gnm, Cnm, wmodesm::Diagonal, wmodes2m::Diagonal,
    ::Type{T} = Float64) where {T<:AbstractFloat}
    K = invLnm + im*(Gnm*wmodesm) - Cnm*wmodes2m
    return SparseMatrixCSC{Complex{T},Int}(K)
end

"""
    applyifft!(td::AbstractArray{T}, fd::AbstractArray{Complex{T}}, irfftplan)

The first half of [`applynl!`](@ref): the physical time domain signal `td`
from the frequency domain coefficients `fd`, with the normalization
convention of `applynl!`, so that linear operations can be interleaved
with the pointwise time domain nonlinearities. [`applyfft!`](@ref) is the
other half, and `applynl!` is their composition with a pointwise function
between them. A `nothing` plan (a system with no junctions) is the identity.

NOTE: `applyifft!` may overwrite `fd`.
"""
function applyifft!(td::AbstractArray{T}, fd::AbstractArray{Complex{T}},
    irfftplan) where T
    mul!(td, irfftplan, fd)
    normalization = prod(size(td)[1:end-1])
    # broadcasting keeps this device generic
    td .*= normalization
    return td
end

# a system with no Josephson junctions has nothing to transform, and
# `plan_applynl` gives it no plan; both arrays are empty, so this is the
# identity rather than a special case of the transform
applyifft!(td::AbstractArray{T}, ::AbstractArray{Complex{T}},
    ::Nothing) where T = td

"""
    applyfft!(fd::AbstractArray{Complex{T}}, td::AbstractArray{T}, rfftplan)

The second half of [`applynl!`](@ref): the frequency domain coefficients
`fd` of the time domain signal `td`, the inverse of [`applyifft!`](@ref)
with the same normalization convention. A `nothing` plan is the identity.

The real to complex transform of an out of place plan reads `td` and
leaves it as it was (FFTW and cuFFT alike; only the complex to real
direction of [`applyifft!`](@ref) may destroy its input), so a caller may
hand in a cached time domain array directly.
"""
function applyfft!(fd::AbstractArray{Complex{T}}, td::AbstractArray{T},
    rfftplan) where T
    mul!(fd, rfftplan, td)
    invnormalization = 1/prod(size(td)[1:end-1])
    # broadcasting keeps this device generic
    fd .*= invnormalization
    return fd
end

# the counterpart of the `applyifft!` method above, for a system with no
# junctions and so no plan
applyfft!(fd::AbstractArray{Complex{T}}, ::AbstractArray{T},
    ::Nothing) where T = fd

"""
    plan_applyffttranspose(fd::Array{Complex{T}}, td::Array{T})

Create the complex transform plan for [`applyffttranspose!`](@ref), the
transpose of [`applyfft!`](@ref) on the same grid. A work array the size of
the time domain grid is allocated for planning and discarded; the caller
passes its own `padded` at apply time. The plan is created with
`FFTW.UNALIGNED` so it can be executed against per-thread work arrays
allocated elsewhere.
"""
function plan_applyffttranspose(fd::Array{Complex{T}}, td::Array{T}) where T
    padded = zeros(Complex{T}, size(td))
    fftplan = FFTW.plan_fft(padded, 1:ndims(td)-1;
        flags = FFTW.ESTIMATE | FFTW.UNALIGNED, timelimit = Inf)
    return fftplan
end

"""
    applyffttranspose!(alpha::Array{Complex{T}}, P::Array{Complex{T}},
        padded::Array{Complex{T}}, fftplan)

Apply the transpose of [`applyfft!`](@ref): given a covector `P` on the
stored frequency domain coefficients, compute the covector `alpha` on the
time domain samples such that `sum(alpha .* td) == sum(P .* applyfft(td))`
for every real time domain array `td`, exactly, including the truncated
first dimension of the real transform.

This is a forward transform again, not a new kernel: the discrete Fourier
matrix is symmetric, so the transpose of the transform which produces the
stored coefficients is the full complex transform of those coefficients
zero padded along the first dimension (the only dimension the real
transform truncates; [`applyfft!`](@ref) stores only its non negative
harmonics), with the same `1/prod(Nt)` normalization. `padded` is a work
array the size of the time domain grid and `alpha` may not alias it. The
last dimension, the Josephson junction index, is not transformed, exactly
as in [`applyfft!`](@ref).

Used by the reverse order sensitivity contraction to transpose the map from
the time domain samples to the Fourier coefficients of `cos(phi(t))`,
through the same transform plans and normalization which define the forward
map, for any number of tones.
"""
function applyffttranspose!(alpha::Array{Complex{T}}, P::Array{Complex{T}},
    padded::Array{Complex{T}}, fftplan) where T
    fill!(padded, 0)
    padded[CartesianIndices(P)] .= P
    mul!(alpha, fftplan, padded)
    invnormalization = 1/prod(size(padded)[1:end-1])
    @inbounds for i in eachindex(alpha)
        alpha[i] = alpha[i]*invnormalization
    end
    return alpha
end

# out .= f.(src), specialized behind a function barrier
function _applypointwise!(out::AbstractArray{T}, f::F,
    src::AbstractArray{T}) where {T,F}
    out .= f.(src)
    return out
end

# Ensure the cached pointwise current-phase relation, or its derivative,
# of the time domain branch fluxes at the current point is up to date.
#
# A circuit whose junctions are all sinusoidal, which is every circuit that
# does not ask for another relation, takes the same single broadcast of
# `sin` or `cos` over the whole array that it always did: the table is
# empty, the branch is one comparison outside the loop, and nothing about
# the junction path changes. A circuit which does ask evaluates every
# junction's polynomial by Horner along the junction axis.
function _ensuresin!(sys::HBSystem)
    if !sys.sincurrent[]
        r = sys.relations
        if allsinusoidal(r)
            _applypointwise!(sys.sintd, sin, sys.phitd)
        else
            applyrelationlast!(sys.sintd, sys.phitd, r.value, r.sinusoidal,
                r.anysinusoidal, sin)
        end
        sys.sincurrent[] = true
    end
    return sys
end

function _ensurecos!(sys::HBSystem)
    if !sys.coscurrent[]
        r = sys.relations
        if allsinusoidal(r)
            _applypointwise!(sys.costd, cos, sys.phitd)
        else
            applyrelationlast!(sys.costd, sys.phitd, r.derivative,
                r.sinusoidal, r.anysinusoidal, cos)
        end
        sys.coscurrent[] = true
    end
    return sys
end

# The array holding the negative of the second derivative of the relation
# at the point, which is what the Hessian and the derivative of the
# linearized system multiply. The sinusoidal relation has `-f'' = sin`, so
# that case returns the cached sine and allocates nothing.
function _negsecond!(sys::HBSystem)
    r = sys.relations
    if allsinusoidal(r)
        _ensuresin!(sys)
        return sys.sintd
    end
    if !sys.negsecondcurrent[]
        applyrelationlast!(sys.negsecondtd, sys.phitd, r.negsecond,
            r.sinusoidal, r.anysinusoidal, sin)
        sys.negsecondcurrent[] = true
    end
    return sys.negsecondtd
end

# The array holding the third derivative of the relation at the point, for
# a circuit which has a polynomial one. The Josephson relation has
# `f''' = -cos`, which its caller reads from the cached cosine instead, so
# this is never called for it.
function _third!(sys::HBSystem)
    r = sys.relations
    if !sys.thirdcurrent[]
        applyrelationlast!(sys.thirdtd, sys.phitd, r.third, r.sinusoidal,
            r.anysinusoidal, x -> -cos(x))
        sys.thirdcurrent[] = true
    end
    return sys.thirdtd
end

"""
    setpoint!(sys::HBSystem, x::AbstractVector)

Set the point at which [`residual!`](@ref), [`jacobianvectorproduct!`](@ref),
[`hessianvectorproduct!`](@ref) and [`jacobian!`](@ref) evaluate the
harmonic balance nonlinear system, and cache the time domain branch fluxes
there. Accepts the complex vector of node fluxes or the equivalent real
representation, dispatched on the element type. Returns `sys`.
"""
function setpoint!(sys::HBSystem, x::AbstractVector{<:Complex})
    copyto!(sys.x, x)
    applycomplextoreal!(sys.xr, sys.nonlineartermplan, sys.x)
    return _setpoint!(sys)
end

function setpoint!(sys::HBSystem, xr::AbstractVector{<:Real})
    copyto!(sys.xr, xr)
    applyrealtocomplex!(sys.x, sys.nonlineartermplan, sys.xr)
    # the time domain branch fluxes come from the forward map of the plan,
    # which reads the real representation directly
    applyforwardterm!(sys.phimatrix, sys.nonlineartermplan, sys.xr)
    applyifft!(sys.phitd, sys.phimatrix, sys.irfftplan)
    sys.sincurrent[] = false
    sys.coscurrent[] = false
    sys.negsecondcurrent[] = false
    sys.thirdcurrent[] = false
    sys.cosfdcurrent[] = false
    return sys
end

function _setpoint!(sys::HBSystem)
    applyforwardterm!(sys.phimatrix, sys.nonlineartermplan, sys.x)
    applyifft!(sys.phitd, sys.phimatrix, sys.irfftplan)
    sys.sincurrent[] = false
    sys.coscurrent[] = false
    sys.negsecondcurrent[] = false
    sys.thirdcurrent[] = false
    sys.cosfdcurrent[] = false
    return sys
end

"""
    residual!(F::AbstractVector, sys::HBSystem)

Evaluate the residual of the harmonic balance nonlinear system,
`F = B(sin.(A*x)) + K*x - b`, at the point set with [`setpoint!`](@ref), in
place. Dispatches on the element type of `F`: a complex vector receives the
complex representation and a real vector the equivalent real representation.
"""
function residual!(F::AbstractVector{<:Complex}, sys::HBSystem)
    _ensuresin!(sys)
    applyfft!(sys.phimatrix, sys.sintd, sys.rfftplan)
    applybackwardterm!(F, sys.nonlineartermplan, sys.phimatrix, sys.x)
    @. F -= sys.bnm
    return F
end

function residual!(Fr::AbstractVector{<:Real}, sys::HBSystem)
    _ensuresin!(sys)
    applyfft!(sys.phimatrix, sys.sintd, sys.rfftplan)
    # the backward map of the plan writes the node vector in the real
    # representation and adds the linear term in the same pass
    applybackwardterm!(Fr, sys.nonlineartermplan, sys.phimatrix, sys.xr)
    Fr .-= sys.bnmr
    return Fr
end

"""
    jacobianvectorproduct!(Jv::AbstractVector, sys::HBSystem,
        v::AbstractVector)

Evaluate the exact matrix-free Jacobian-vector product of the harmonic
balance nonlinear system, `J(x)*v = B(cos.(A*x) .* (A*v)) + K*v`, at the
point set with [`setpoint!`](@ref), in place. Dispatches on the element
types: complex vectors receive the complex representation and real vectors
the equivalent real representation, for which the product equals the
assembled real Jacobian ([`jacobian!`](@ref)) applied to `vr` up to floating
point roundoff, including the self-conjugate (eg. DC) modes. Each product
costs two Fourier transforms and the linear term; the time domain cosine is
cached across products at the same point. Suitable as the operator for
Krylov methods.
"""
function jacobianvectorproduct!(Jv::AbstractVector{<:Complex},
    sys::HBSystem, v::AbstractVector{<:Complex})
    _ensurecos!(sys)
    plan = sys.nonlineartermplan
    applyforwardterm!(sys.phimatrix, plan, v)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, sys.costd, sys.dirtd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    applybackwardterm!(Jv, plan, sys.phimatrix, v)
    return Jv
end

function jacobianvectorproduct!(Jvr::AbstractVector{<:Real},
    sys::HBSystem, vr::AbstractVector{<:Real})
    _ensurecos!(sys)
    plan = sys.nonlineartermplan
    applyforwardterm!(sys.phimatrix, plan, vr)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, sys.costd, sys.dirtd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    applybackwardterm!(Jvr, plan, sys.phimatrix, vr)
    return Jvr
end

# out .= a .* b and out .= -a .* b .* c behind function barriers
function _multiplyintowork!(out::AbstractArray{T}, a::AbstractArray{T},
    b::AbstractArray{T}) where T
    out .= a .* b
    return out
end

function _multiplyintowork!(out::AbstractArray{T}, a::AbstractArray{T},
    b::AbstractArray{T}, c::AbstractArray{T}) where T
    out .= .-a .* b .* c
    return out
end

"""
    hessianvectorproduct!(Hvw::AbstractVector, sys::HBSystem,
        v::AbstractVector, w::AbstractVector)

Evaluate the exact matrix-free second directional derivative of the harmonic
balance nonlinear system, `H(x)[v, w] = B(-sin.(A*x) .* (A*v) .* (A*w))`, at
the point set with [`setpoint!`](@ref), in place. The frequency dependent
linear terms are linear in `x` so they do not contribute. Dispatches on the
element types: complex vectors receive the complex representation and real
vectors the equivalent real representation. The product is symmetric in `v`
and `w`. Useful for continuation and bifurcation tracking methods which
require directional second derivatives.
"""
function hessianvectorproduct!(Hvw::AbstractVector{<:Complex},
    sys::HBSystem, v::AbstractVector{<:Complex},
    w::AbstractVector{<:Complex})
    negsecond = _negsecond!(sys)
    plan = sys.nonlineartermplan
    applyforwardterm!(sys.phimatrix, plan, v)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    applyforwardterm!(sys.phimatrix, plan, w)
    applyifft!(sys.dirtd2, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, negsecond, sys.dirtd, sys.dirtd2)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    # the linear terms are linear in x so they do not contribute
    applybackwardterm!(Hvw, plan, sys.phimatrix, v; addlinearterm = false)
    return Hvw
end

function hessianvectorproduct!(Hvwr::AbstractVector{<:Real},
    sys::HBSystem, vr::AbstractVector{<:Real}, wr::AbstractVector{<:Real})
    negsecond = _negsecond!(sys)
    plan = sys.nonlineartermplan
    applyforwardterm!(sys.phimatrix, plan, vr)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    applyforwardterm!(sys.phimatrix, plan, wr)
    applyifft!(sys.dirtd2, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, negsecond, sys.dirtd, sys.dirtd2)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    # the linear terms are linear in x so they do not contribute
    applybackwardterm!(Hvwr, plan, sys.phimatrix, vr; addlinearterm = false)
    return Hvwr
end

"""
    cosdirectionalderivative!(dcos::Array, sys::HBSystem, v::AbstractVector)

Evaluate the directional derivative, along `v`, of the Fourier coefficients
of `cos(phi_b(t))` of the Josephson junction branch fluxes at the point set
with [`setpoint!`](@ref), in place. Those coefficients parameterize the pump
modulation of the linearized harmonic balance system (see
[`HBLinearizedSystem`](@ref)), so this is what propagates a shift of the pump
operating point into the linearized system matrix. The derivative has the
coefficients of `-sin.(A*x).*(A*v)`, computed on the same time grid and with
the same normalization as the residual, so no separate transform convention
is introduced. Shares the cached time domain sine with
[`hessianvectorproduct!`](@ref), of which this is the first half.

`dcos` is a host array, and so therefore is `sys`: the only caller is the
sensitivity contraction, and an [`HBOperatingPoint`](@ref) holds a host
system whichever backend solved for it.
"""
function cosdirectionalderivative!(dcos::Array, sys::HBSystem,
    v::AbstractVector{<:Complex})
    negsecond = _negsecond!(sys)
    applyforwardterm!(sys.phimatrix, sys.nonlineartermplan, v)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, negsecond, sys.dirtd)
    applyfft!(dcos, sys.worktd, sys.rfftplan)
    @inbounds for i in eachindex(dcos)
        dcos[i] = -dcos[i]
    end
    return dcos
end

"""
    tohost(x::AbstractArray)

The array on the host, for the host loops which cannot read a device array.
Returns `x` itself when it is already an `Array`, so nothing about the CPU
path changes; anything else is copied back with `Array`.
"""
tohost(x::Array) = x
tohost(x::AbstractArray) = Array(x)

"""
    jacobian!(Jx::SparseMatrixCSC{<:Complex}, sys::HBSystem)
    jacobian!(Jr::SparseMatrixCSC{<:Real}, sys::HBSystem)
    jacobian!(A::DeviceValuedSparseMatrix, sys::HBSystem)
    jacobian!(A::DeviceValuedSparseMatrix{<:Complex}, plan, sys::HBSystem)

Assemble the Jacobian of the harmonic balance nonlinear system at the point
set with [`setpoint!`](@ref), in place, using the precomputed plans.
Dispatches on the element type of the matrix: a complex matrix receives the complex
holomorphic Jacobian (an approximation to the exact Jacobian, used by the
[`QuasiNewton`](@ref) method, via [`assemblecomplexjacobian!`](@ref)) and a
real matrix the exact Jacobian of the equivalent real system (used by the
[`Newton`](@ref) method, via [`assemblerealjacobian!`](@ref)). A
[`DeviceValuedSparseMatrix`](@ref) receives the same on its backend. The
corresponding plan must have been provided when the [`HBSystem`](@ref) was
constructed; the method taking a plan explicitly is documented below.
"""
function jacobian!(Jx::SparseMatrixCSC{<:Complex}, sys::HBSystem)
    isnothing(sys.complexjacobianplan) && throw(ArgumentError(
        "no complex Jacobian plan was provided to this HBSystem."))
    assemblecomplexjacobian!(nonzeros(Jx), sys.complexjacobianplan,
        tohost(cosphimatrix(sys)))
    return Jx
end

function jacobian!(Jr::SparseMatrixCSC{<:Real}, sys::HBSystem)
    isnothing(sys.realjacobianplan) && throw(ArgumentError(
        "no real Jacobian plan was provided to this HBSystem."))
    return jacobian!(Jr, sys.realjacobianplan, sys)
end

"""
    jacobian!(Jr::SparseMatrixCSC, plan::StructureRealJacobianPlan,
        sys::HBSystem)

Assemble the real Jacobian described by `plan` at the point set with
[`setpoint!`](@ref), in place. This is the method [`jacobian!`](@ref)`(Jr, sys)`
delegates to with the plan stored in `sys`; taking the plan explicitly lets a
*different* plan, built over the same system, be assembled from the same Fourier
coefficients of `cos(phi(t))` and the same linear term matrices.
"""
function jacobian!(Jr::SparseMatrixCSC{<:Real},
    plan::StructureRealJacobianPlan, sys::HBSystem)
    # the assembly writes a host sparse matrix, so on a device backend the
    # coefficients it reads come back first. This is the one host round trip
    # left on this path, and it is the exact Jacobian of the direct solve
    # rather than anything in the Krylov iteration.
    assemblerealjacobian!(nonzeros(Jr), plan, tohost(cosphimatrix(sys)))
    return Jr
end

"""
    jacobian!(A::DeviceValuedSparseMatrix, plan::StructureRealJacobianPlan,
        sys::HBSystem)

Assemble the real Jacobian into the device values of `A`. Nothing crosses to
the host: the coefficients are already on the backend and the assembly runs
there.
"""
function jacobian!(A::DeviceValuedSparseMatrix,
    plan::StructureRealJacobianPlan, sys::HBSystem)
    assemblerealjacobian!(A.nzval, plan, cosphimatrix(sys))
    return A
end

function jacobian!(A::DeviceValuedSparseMatrix, sys::HBSystem)
    isnothing(sys.realjacobianplan) && throw(ArgumentError(
        "no real Jacobian plan was provided to this HBSystem."))
    return jacobian!(A, sys.realjacobianplan, sys)
end

function jacobian!(A::DeviceValuedSparseMatrix{<:Complex},
    plan::StructureComplexJacobianPlan, sys::HBSystem)
    assemblecomplexjacobian!(A.nzval, plan, cosphimatrix(sys))
    return A
end

function jacobian!(A::DeviceValuedSparseMatrix{<:Complex}, sys::HBSystem)
    isnothing(sys.complexjacobianplan) && throw(ArgumentError(
        "no complex Jacobian plan was provided to this HBSystem."))
    return jacobian!(A, sys.complexjacobianplan, sys)
end

"""
    cosphimatrix(sys::HBSystem)

The Fourier coefficients of `cos(phi(t))` at the current point, from which
the Jacobians are assembled: transformed once per point into a buffer of
the system's own, which no residual or product touches, and returned as is
thereafter. Read it, do not write it.
"""
function cosphimatrix(sys::HBSystem)
    if !sys.cosfdcurrent[]
        _ensurecos!(sys)
        applyfft!(sys.cosfd, sys.costd, sys.rfftplan)
        sys.cosfdcurrent[] = true
    end
    return sys.cosfd
end
