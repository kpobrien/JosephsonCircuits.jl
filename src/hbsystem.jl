
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
:quasinewton and :newton methods of [`hbnlsolve`](@ref) as well as
matrix-free solvers.

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
struct HBSystem{TR,TinvL,TG,TC,TWm,TK,Tb,TLjb,TLjbm,TLm,TIP,TFP,TRJ,TCJ,TNP,TVC,TVR,TAC,TAR}
    # linear term matrices and source vector (complex representation, scaled
    # by Lmean, conjugated for negative frequency modes with conjnegfreq!)
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
    Lmean::TLm
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
    # workspaces
    phimatrix::TAC
    dirtd::TAR
    dirtd2::TAR
    worktd::TAR
end

"""
    HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm, Ljb,
        Ljbm, Lmean, Nbranches, freqindexmap, conjsourceindices, conjtargetindices,
        phimatrix, phimatrixtd, irfftplan, rfftplan, modelayout,
        realjacobianplan, complexjacobianplan, backend = CPU())

Construct an [`HBSystem`](@ref) from the ingredients assembled by
[`hbnlsolve`](@ref), allocating the workspaces. `phimatrix` and
`phimatrixtd` are adopted as the frequency domain and one of the time domain
workspaces. `realjacobianplan` and `complexjacobianplan` may be `nothing`,
in which case the corresponding [`jacobian!`](@ref) method is unavailable.
`backend` is the KernelAbstractions backend on which the index mapped
kernels run and on which the plan and the workspaces of the residual and the
matrix-free products are allocated; `CPU()` is the default and the reference.
`phimatrix` and `phimatrixtd` are adopted as given, so pass them already on
the backend, and the time domain workspaces derived from them with `similar`
follow.
"""
function HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
    Ljb, Ljbm, Lmean, Nbranches, freqindexmap, conjsourceindices, conjtargetindices,
    phimatrix, phimatrixtd, irfftplan, rfftplan, modelayout,
    realjacobianplan, complexjacobianplan, backend = CPU();
    realbackward::Bool = true)

    n = size(Rbnm, 2)
    # the working precision of everything the residual and the matrix-free
    # products touch is the one of the frequency domain array handed in, so
    # allocating `phimatrix` as Complex{Float32} carries the whole evaluation
    # path to Float32. The host side linear term matrices are left as they
    # were: they are only read by the assembled Jacobian plans, which run on
    # the host and write into a matrix of their own element type.
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
    nonlineartermplan = plannonlinearterm(Rbnm, Ljb, Lmean, Nbranches,
        freqindexmap, conjsourceindices, conjtargetindices, phimatrix, Knm,
        modelayout, backend; realbackward = realbackward)
    # both representations of the source vector move to the backend. the real
    # one already did; the complex one did not, so `residual!` on a complex
    # vector broadcast a device array against a host one. `complex_to_real`
    # below reads the host argument, which `tobackend` does not modify.
    return HBSystem(Rbnm, invLnm, Gnm, Cnm, wmodesm, wmodes2m, Knm,
        tobackend(backend, convert(Vector{Complex{TF}}, bnm)),
        Ljb, Ljbm, Lmean, freqindexmap, conjsourceindices,
        conjtargetindices, irfftplan, rfftplan, modelayout,
        realjacobianplan, complexjacobianplan, nonlineartermplan,
        tobackend(backend, convert(Vector{TF},
            complex_to_real(bnm, modelayout.isreal))),
        tobackend(backend, zeros(Complex{TF}, n)),
        tobackend(backend, zeros(TF, modelayout.rdim)),
        similar(phimatrixtd), similar(phimatrixtd), similar(phimatrixtd),
        Ref(false), Ref(false),
        phimatrix, similar(phimatrixtd), similar(phimatrixtd), phimatrixtd)
end


"""
    linearterm(invLnm, Gnm, Cnm, wmodesm::Diagonal, wmodes2m::Diagonal)

Collapse the frequency dependent linear terms of the harmonic balance system
into the single sparse matrix

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

applyfft! and applyifft! and the two halves of applynl!, so linear operations
can be interleaved with the pointwise time domain nonlinearities. applyifft!
computes the physical time domain signal from the frequency domain
coefficients and applyfft! the inverse, with the same normalization convention
as applynl!, of which applynl! is the composition with a pointwise function
between the halves.

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

applyfft! and applyifft! and the two halves of applynl!, so linear operations
can be interleaved with the pointwise time domain nonlinearities. applyifft!
computes the physical time domain signal from the frequency domain
coefficients and applyfft! the inverse, with the same normalization convention
as applynl!, of which applynl! is the composition with a pointwise function
between the halves.

NOTE: `applyifft!` may overwrite `td`.
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

Create the padded work array and the complex transform plan for
[`applyffttranspose!`](@ref), the transpose of [`applyfft!`](@ref) on the
same grid. The plan is created with `FFTW.UNALIGNED` so it can be executed
against per-thread work arrays allocated elsewhere.
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

# ensure the cached pointwise sine or cosine of the time domain branch
# fluxes at the current point are up to date.
function _ensuresin!(sys::HBSystem)
    if !sys.sincurrent[]
        _applypointwise!(sys.sintd, sin, sys.phitd)
        sys.sincurrent[] = true
    end
    return sys
end

function _ensurecos!(sys::HBSystem)
    if !sys.coscurrent[]
        _applypointwise!(sys.costd, cos, sys.phitd)
        sys.coscurrent[] = true
    end
    return sys
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
    return sys
end

function _setpoint!(sys::HBSystem)
    applyforwardterm!(sys.phimatrix, sys.nonlineartermplan, sys.x)
    applyifft!(sys.phitd, sys.phimatrix, sys.irfftplan)
    sys.sincurrent[] = false
    sys.coscurrent[] = false
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
    copyto!(sys.worktd, sys.sintd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    applybackwardterm!(F, sys.nonlineartermplan, sys.phimatrix, sys.x)
    @. F -= sys.bnm
    return F
end

function residual!(Fr::AbstractVector{<:Real}, sys::HBSystem)
    _ensuresin!(sys)
    copyto!(sys.worktd, sys.sintd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
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
    _ensuresin!(sys)
    plan = sys.nonlineartermplan
    applyforwardterm!(sys.phimatrix, plan, v)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    applyforwardterm!(sys.phimatrix, plan, w)
    applyifft!(sys.dirtd2, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, sys.sintd, sys.dirtd, sys.dirtd2)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    # the linear terms are linear in x so they do not contribute
    applybackwardterm!(Hvw, plan, sys.phimatrix, v; addlinearterm = false)
    return Hvw
end

function hessianvectorproduct!(Hvwr::AbstractVector{<:Real},
    sys::HBSystem, vr::AbstractVector{<:Real}, wr::AbstractVector{<:Real})
    _ensuresin!(sys)
    plan = sys.nonlineartermplan
    applyforwardterm!(sys.phimatrix, plan, vr)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    applyforwardterm!(sys.phimatrix, plan, wr)
    applyifft!(sys.dirtd2, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, sys.sintd, sys.dirtd, sys.dirtd2)
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
    _ensuresin!(sys)
    applyforwardterm!(sys.phimatrix, sys.nonlineartermplan, v)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    _multiplyintowork!(sys.worktd, sys.sintd, sys.dirtd)
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
    jacobian!(J::SparseMatrixCSC, sys::HBSystem)

Assemble the Jacobian of the harmonic balance nonlinear system at the point
set with [`setpoint!`](@ref), in place, using the precomputed plans.
Dispatches on the element type of `J`: a complex matrix receives the complex
holomorphic Jacobian (an approximation to the exact Jacobian, used by the
:quasinewton method, via [`assemblecomplexjacobian!`](@ref)) and a real
matrix the exact Jacobian of the equivalent real system (used by the :newton
method, via [`assemblerealjacobian!`](@ref)). The corresponding plan must
have been provided when the [`HBSystem`](@ref) was constructed.
"""
function jacobian!(Jx::SparseMatrixCSC{<:Complex}, sys::HBSystem)
    isnothing(sys.complexjacobianplan) && throw(ArgumentError(
        "no complex Jacobian plan was provided to this HBSystem."))
    _updatecosphimatrix!(sys)
    assemblecomplexjacobian!(Jx, sys.complexjacobianplan, tohost(sys.phimatrix),
        sys.invLnm, sys.Gnm, sys.Cnm, sys.wmodesm, sys.wmodes2m)
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
    _updatecosphimatrix!(sys)
    # the assembly writes a host sparse matrix, so on a device backend the
    # coefficients it reads come back first. This is the one host round trip
    # left on this path, and it is the exact Jacobian of the direct solve
    # rather than anything in the Krylov iteration.
    assemblerealjacobian!(nonzeros(Jr), plan, tohost(sys.phimatrix))
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
    _updatecosphimatrix!(sys)
    assemblerealjacobian!(A.nzval, plan, sys.phimatrix)
    return A
end

function jacobian!(A::DeviceValuedSparseMatrix, sys::HBSystem)
    isnothing(sys.realjacobianplan) && throw(ArgumentError(
        "no real Jacobian plan was provided to this HBSystem."))
    return jacobian!(A, sys.realjacobianplan, sys)
end

function jacobian!(A::DeviceValuedSparseMatrix{<:Complex},
    plan::StructureComplexJacobianPlan, sys::HBSystem)
    _updatecosphimatrix!(sys)
    assemblecomplexjacobian!(A.nzval, plan, sys.phimatrix)
    return A
end

function jacobian!(A::DeviceValuedSparseMatrix{<:Complex}, sys::HBSystem)
    isnothing(sys.complexjacobianplan) && throw(ArgumentError(
        "no complex Jacobian plan was provided to this HBSystem."))
    return jacobian!(A, sys.complexjacobianplan, sys)
end

# update sys.phimatrix with the Fourier coefficients of cos.(A*x) at the
# current point, from which the Jacobians are assembled.
function _updatecosphimatrix!(sys::HBSystem)
    _ensurecos!(sys)
    copyto!(sys.worktd, sys.costd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    return sys
end

"""
    HBLinearizedSystem

Everything needed to assemble the linearized harmonic balance system matrix
on the signal frequency grid at any signal frequency, built with the same
plan machinery ([`plancomplexjacobian`](@ref) and
[`addjosephsonterm!`](@ref)) used for the Jacobians of the nonlinear system
by [`HBSystem`](@ref). [`hblinsolve`](@ref) constructs one of these and its
per-frequency solves are expressed through it, so any external use of this
object is exercising exactly the production linearized solver path.

Use [`assemblesystemmatrix!`](@ref) to assemble

    A(ws) = AoLjnm + invLnm + im*Gnm*W - Cnm*W^2

(with the negative frequency mode conjugations and symbolic frequency
substitutions of the linearized solver) into a matrix sharing the sparsity
structure of the `Asparse` field, at the signal frequency `ws`, for either
the pump modulation `AoLjnm` or its complex conjugate (the adjoint system of
the noise and quantum efficiency calculations). Operator products and
adjoint products are then sparse matrix-vector products with the assembled
matrix, `mul!(y, A, v)` and `mul!(y, A', v)`, suitable for iterative solves
and sensitivity adjoints.

Unlike [`HBSystem`](@ref), no time domain matrix-free product is provided,
deliberately: a Fourier transform based product computes a cyclic
convolution on the pump grid, whereas the assembled matrix uses the explicit
truncation of [`hbmatind`](@ref) (the zeros of `Amatrixindices`), dropping
the couplings which fall outside the grid instead of wrapping them. The
assembled matrix defines the linearized solver, so products must match it
exactly; since the pump modulation contribution is precomputed, assembling
at a frequency costs about the same as one matrix-vector product would, and
every product thereafter is a plain sparse matrix-vector product.
"""
struct HBLinearizedSystem{TinvL,TG,TC,TAG}
    # sparsity structure of the system matrix, with the pump modulation
    # contribution in its values. per-frequency assembly operates on copies
    # sharing this structure, so this object can be shared across threads.
    Asparse
    # the plan which built Asparse, holding the scatter lists for the pump
    # modulation contribution
    complexjacobianplan
    # the pump modulation contribution AoLjnm = Rbnm'*AoLjbm*Rbnm and its
    # complex conjugate, assembled once as nonzero value vectors aligned
    # with the sparsity structure of Asparse
    AoLjnmnzval::Vector{Complex{Float64}}
    AoLjnmconjnzval::Vector{Complex{Float64}}
    # the frequency dependent linear term matrices (possibly containing
    # symbolic frequency variables), their index maps into Asparse, and
    # the indices of their symbolic entries
    invLnm::TinvL
    Gnm::TG
    Cnm::TC
    invLnmindexmap::Vector{Int}
    Gnmindexmap::Vector{Int}
    Cnmindexmap::Vector{Int}
    invLnmfreqsubstindices
    Gnmfreqsubstindices
    Cnmfreqsubstindices
    # the modified nodal analysis augmentation of the promoted resistors
    # (see calcAmnasplit): the frequency independent entries, the
    # conductances of the constitutive equations (added per frequency as
    # im*AmnaG*wmodesm), their index maps, and the indices of any symbolic
    # entries of AmnaG. empty matrices when no resistors are promoted.
    Amna0::SparseMatrixCSC{Complex{Float64},Int}
    AmnaG::TAG
    Amna0indexmap::Vector{Int}
    AmnaGindexmap::Vector{Int}
    AmnaGfreqsubstindices::Vector{Int}
    # the symbolic frequency variable, or nothing
    symfreqvar
    # the pump mode frequency offsets of the signal modes, and the numbers
    # of modes and nodes, for computing the mode frequency matrices
    wpumpmodes::Vector{Float64}
    Nmodes::Int
    Nnodes::Int
    # the multiport admittance contribution of the scattering block
    # components (see ScatteringStampSystem), or nothing
    scattering
end

"""
    HBLinearizedSystem(Amatrixindices::Matrix, Ljb::SparseVector,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        phimatrix::Array, invLnmcopy::SparseMatrixCSC,
        Gnmcopy::SparseMatrixCSC, Cnmcopy::SparseMatrixCSC, invLnm, Gnm,
        Cnm, invLnmfreqsubstindices, Gnmfreqsubstindices,
        Cnmfreqsubstindices, symfreqvar, wpumpmodes, Nnodes::Integer)

Construct an [`HBLinearizedSystem`](@ref) from the signal frequency grid
index matrix `Amatrixindices` (see [`hbmatind`](@ref)), the Josephson
junction data, the Fourier coefficients of `cos(phi(t))` of the pump in
`phimatrix`, the numeric copies of the linear term matrices (which define
the sparsity structure) and the possibly symbolic originals with their
symbolic entry indices. Builds the sparsity structure and plan with
[`plancomplexjacobian`](@ref) and assembles the pump modulation contribution
and its conjugate with [`addjosephsonterm!`](@ref).
"""
function HBLinearizedSystem(Amatrixindices::Matrix, Ljb::SparseVector,
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    phimatrix::Array, invLnmcopy::SparseMatrixCSC,
    Gnmcopy::SparseMatrixCSC, Cnmcopy::SparseMatrixCSC, invLnm, Gnm, Cnm,
    invLnmfreqsubstindices, Gnmfreqsubstindices, Cnmfreqsubstindices,
    Amna0::SparseMatrixCSC, AmnaG::SparseMatrixCSC,
    symfreqvar, wpumpmodes, Nnodes::Integer; scattering = nothing)

    # the sparsity structure must contain the modified nodal analysis
    # augmentation as well, so merge its entries into the numeric copies
    # used only for the structure. the source index stride is the size of
    # the pump frequency grid, the leading dimensions of phimatrix.
    invLnmpattern = spaddkeepzeros(spaddkeepzeros(invLnmcopy, Amna0), AmnaG)
    if !isnothing(scattering)
        # the scattering block contribution shares the structure
        invLnmpattern = spaddkeepzeros(invLnmpattern, scattering.pattern)
    end
    Asparse, complexjacobianplan = plancomplexjacobian(Amatrixindices, Ljb,
        1, Rbnm, Nmodes, Nbranches, prod(size(phimatrix)[1:end-1]),
        invLnmpattern, Gnmcopy, Cnmcopy)

    # the plan's linear index maps were built for the merged pattern
    # matrix, so build the runtime index maps for the individual matrices.
    invLnmindexmap = sparseaddmap(Asparse, invLnmcopy)
    Gnmindexmap = sparseaddmap(Asparse, Gnmcopy)
    Cnmindexmap = sparseaddmap(Asparse, Cnmcopy)
    Amna0indexmap = sparseaddmap(Asparse, Amna0)
    AmnaGindexmap = sparseaddmap(Asparse, AmnaG)
    AmnaGfreqsubstindices = symbolicindices(AmnaG)
    if !isnothing(scattering)
        setscatteringindexmap!(scattering, Asparse)
    end

    # assemble the pump modulation contribution and its complex conjugate
    # (for the adjoint system) once, so resetting the system matrix at each
    # frequency is a single copy.
    AoLjnmnzval = zeros(Complex{Float64}, nnz(Asparse))
    addjosephsonterm!(AoLjnmnzval, complexjacobianplan, phimatrix)
    AoLjnmconjnzval = conj.(AoLjnmnzval)

    # also place the pump modulation contribution in Asparse so it holds
    # something reasonable to factorize.
    copyto!(Asparse.nzval, AoLjnmnzval)

    return HBLinearizedSystem(Asparse, complexjacobianplan, AoLjnmnzval,
        AoLjnmconjnzval, invLnm, Gnm, Cnm, invLnmindexmap, Gnmindexmap,
        Cnmindexmap, invLnmfreqsubstindices, Gnmfreqsubstindices,
        Cnmfreqsubstindices, Amna0, AmnaG, Amna0indexmap, AmnaGindexmap,
        AmnaGfreqsubstindices, symfreqvar, wpumpmodes, Nmodes, Nnodes,
        scattering)
end

"""
    assemblesystemmatrix!(A::SparseMatrixCSC, lsys::HBLinearizedSystem,
        wmodes::AbstractVector; conjugatepump::Bool = false)
    assemblesystemmatrix!(A::SparseMatrixCSC, lsys::HBLinearizedSystem,
        ws::Number; conjugatepump::Bool = false)

Assemble the linearized harmonic balance system matrix into `A`, which must
share the sparsity structure of `lsys.Asparse`, either from the mode
frequency vector `wmodes` or at the signal frequency `ws` (from which the
mode frequencies are computed as `wmodes = ws .+ lsys.wpumpmodes`). The
frequency scaling, the negative frequency conjugation, and any symbolic
frequency substitution are applied per column from the mode index, without
materializing system sized diagonals. With `conjugatepump = true` the complex
conjugate of the pump modulation contribution is used, which for a circuit
without scattering blocks is a similarity transformation of the transposed
system, as below. The
negative frequency mode entries of the linear term matrices are conjugated
and any symbolic frequency variables substituted, exactly as in the
per-frequency loop of [`hblinsolve`](@ref), which calls this function.
Returns `A`.

The conjugated pump system is a diagonal similarity transformation of the
transposed forward system,

    A(conjugate pump) = D*transpose(A(pump))*inv(D),

with `D` diagonal, equal to one on every node flux row and, on the auxiliary
branch current rows of the modified nodal analysis augmentation, equal to the
assembled conductance entry `im*w_m/R_r` of the constitutive equation of the
corresponding promoted resistor and mode. (The promoted resistances are
constant and real by construction, see [`mnaportresistorindices`](@ref), so
the negative frequency conjugation of `sparseaddconjsubst!`, which acts on
the stored conductance and not on the frequency factor, is trivial here.)
Nothing else contributes, so long as the circuit has no scattering blocks: the
auxiliary rows of the promoted coupled inductors are already symmetric (see
[`calcAmnaind`](@ref)); the linear term matrices are symmetric and mode
diagonal, so the column indexed frequency scaling and conjugation of
`sparseaddconjsubst!` are symmetric under transposition; and the transpose of
the pump modulation contribution exchanges the mode pair, mapping each
difference harmonic to its complex conjugate, which is the same as
conjugating the pump. The hybrid rows of a [`ScatteringBlock`](@ref) do break
it: their constant Kirchhoff couplings and frequency dependent constitutive
entries exchange under transposition, and no diagonal `D` undoes that, so with
blocks the two are different matrices.

The adjoint solutions the noise, quantum efficiency, commutation relation and
adjoint node output calculations read are in every case the solutions of the
*transposed* system: by the adjoint identity, the response at an output port
to a source anywhere in the circuit is that source contracted against the
transposed solution driven at the port. [`hblinsolve`](@ref) obtains them with
[`trysolvetranspose!`](@ref) on the factorization of the forward system, which
costs a pair of triangular solves rather than an assembly and a factorization
at every signal frequency. Where the similarity holds, the conjugated pump
assembly is the independent construction that equivalence is tested against.
"""
function assemblesystemmatrix!(A::SparseMatrixCSC,
    lsys::HBLinearizedSystem, wmodes::AbstractVector;
    conjugatepump::Bool = false)

    # the pump modulation contribution, precomputed
    if conjugatepump
        copyto!(A.nzval, lsys.AoLjnmconjnzval)
    else
        copyto!(A.nzval, lsys.AoLjnmnzval)
    end

    # take the complex conjugate of the negative frequency terms in
    # the capacitance and conductance matrices. substitute in the symbolic
    # frequency variable if present. the frequency scaling of each term is
    # the per column mode frequency raised to the given power.
    sparseaddconjsubst!(A, -1, lsys.Cnm, lsys.Cnmindexmap, wmodes, 2,
        lsys.Cnmfreqsubstindices, lsys.symfreqvar)
    sparseaddconjsubst!(A, im, lsys.Gnm, lsys.Gnmindexmap, wmodes, 1,
        lsys.Gnmfreqsubstindices, lsys.symfreqvar)
    sparseaddconjsubst!(A, 1, lsys.invLnm, lsys.invLnmindexmap, wmodes, 0,
        lsys.invLnmfreqsubstindices, lsys.symfreqvar)

    # the modified nodal analysis augmentation of the promoted resistors:
    # the frequency independent entries and the conductances of the
    # constitutive equations, the latter with the same per-mode frequency
    # scaling and negative frequency conjugation as the conductance matrix.
    sparseadd!(A, 1, lsys.Amna0, lsys.Amna0indexmap)
    sparseaddconjsubst!(A, im, lsys.AmnaG, lsys.AmnaGindexmap, wmodes, 1,
        lsys.AmnaGfreqsubstindices, lsys.symfreqvar)

    # the multiport admittance contribution of the scattering block
    # components: sign * im * w_m * Y[p,q](w_m) per contribution, the
    # frequency dependent generalization of the conductance term
    if !isnothing(lsys.scattering)
        assemblescattering!(A, lsys.scattering, wmodes)
    end

    return A
end

function assemblesystemmatrix!(A::SparseMatrixCSC,
    lsys::HBLinearizedSystem, ws::Number; conjugatepump::Bool = false)
    return assemblesystemmatrix!(A, lsys, ws .+ lsys.wpumpmodes;
        conjugatepump = conjugatepump)
end
