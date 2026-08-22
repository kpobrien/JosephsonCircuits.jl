
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
types.
"""
struct HBSystem
    # linear term matrices and source vector (complex representation, scaled
    # by Lmean, conjugated for negative frequency modes with conjnegfreq!)
    Rbnm
    Rbnmt
    invLnm
    Gnm
    Cnm
    wmodesm
    wmodes2m
    bnm
    # Josephson junction data
    Ljb
    Ljbm
    Lmean
    # frequency domain packing and Fourier transform plans
    freqindexmap
    conjsourceindices
    conjtargetindices
    irfftplan
    rfftplan
    # real representation layout
    modelayout
    # optional precomputed Jacobian assembly plans
    realjacobianplan
    complexjacobianplan
    # the current point (set with setpoint!)
    x::Vector{Complex{Float64}}
    # cached time domain branch fluxes and pointwise nonlinearities at the
    # current point, with validity flags
    phitd
    sintd
    costd
    sincurrent::Base.RefValue{Bool}
    coscurrent::Base.RefValue{Bool}
    # workspaces
    phimatrix
    dirtd
    dirtd2
    worktd
    AoLjbmvector::Vector{Complex{Float64}}
    workn1::Vector{Complex{Float64}}
    workn2::Vector{Complex{Float64}}
    workv::Vector{Complex{Float64}}
    workw::Vector{Complex{Float64}}
    workF::Vector{Complex{Float64}}
end

"""
    HBSystem(Rbnm, Rbnmt, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm, Ljb,
        Ljbm, Lmean, freqindexmap, conjsourceindices, conjtargetindices,
        phimatrix, phimatrixtd, irfftplan, rfftplan, modelayout,
        realjacobianplan, complexjacobianplan)

Construct an [`HBSystem`](@ref) from the ingredients assembled by
[`hbnlsolve`](@ref), allocating the workspaces. `phimatrix` and
`phimatrixtd` are adopted as the frequency domain and one of the time domain
workspaces. `realjacobianplan` and `complexjacobianplan` may be `nothing`,
in which case the corresponding [`jacobian!`](@ref) method is unavailable.
"""
function HBSystem(Rbnm, Rbnmt, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
    Ljb, Ljbm, Lmean, freqindexmap, conjsourceindices, conjtargetindices,
    phimatrix, phimatrixtd, irfftplan, rfftplan, modelayout,
    realjacobianplan, complexjacobianplan)

    n = size(Rbnm, 2)
    return HBSystem(Rbnm, Rbnmt, invLnm, Gnm, Cnm, wmodesm, wmodes2m, bnm,
        Ljb, Ljbm, Lmean, freqindexmap, conjsourceindices,
        conjtargetindices, irfftplan, rfftplan, modelayout,
        realjacobianplan, complexjacobianplan,
        zeros(Complex{Float64}, n),
        similar(phimatrixtd), similar(phimatrixtd), similar(phimatrixtd),
        Ref(false), Ref(false),
        phimatrix, similar(phimatrixtd), similar(phimatrixtd), phimatrixtd,
        zeros(Complex{Float64}, size(Rbnm, 1)),
        zeros(Complex{Float64}, n), zeros(Complex{Float64}, n),
        zeros(Complex{Float64}, n), zeros(Complex{Float64}, n),
        zeros(Complex{Float64}, n))
end

# the two halves of applynl!, so linear operations can be interleaved with
# the pointwise time domain nonlinearities. applyifft! computes the physical
# time domain signal from the frequency domain coefficients and applyfft!
# the inverse, with the same normalization convention as applynl!, of which
# applynl! is the composition with a pointwise function between the halves.
# NOTE: the multidimensional real Fourier transform may destroy its input,
# so only pass workspaces to applyfft!.
function applyifft!(td::Array{T}, fd::Array{Complex{T}}, irfftplan) where T
    mul!(td, irfftplan, fd)
    normalization = prod(size(td)[1:end-1])
    @inbounds for i in eachindex(td)
        td[i] = td[i]*normalization
    end
    return td
end

function applyfft!(fd::Array{Complex{T}}, td::Array{T}, rfftplan) where T
    mul!(fd, rfftplan, td)
    invnormalization = 1/prod(size(td)[1:end-1])
    @inbounds for i in eachindex(fd)
        fd[i] = fd[i]*invnormalization
    end
    return fd
end

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
function _applypointwise!(out::Array{T}, f::F, src::Array{T}) where {T,F}
    @inbounds for i in eachindex(out)
        out[i] = f(src[i])
    end
    return out
end

# compute the physical time domain branch fluxes on the Josephson junctions
# for the complex vector z (a point or a direction) into td, using
# sys.phimatrix as the frequency domain workspace.
function _branchtimedomain!(td, sys::HBSystem, z::AbstractVector)
    phib = sys.Rbnm*z
    phivectortomatrix!(phib[sys.Ljbm.nzind], sys.phimatrix,
        sys.freqindexmap, sys.conjsourceindices, sys.conjtargetindices,
        length(sys.Ljb.nzval))
    applyifft!(td, sys.phimatrix, sys.irfftplan)
    return td
end

# pack the frequency domain coefficients in sys.phimatrix into the branch
# vector, scale by Lmean/Lj, and apply the transpose incidence matrix:
# out .= Rbnmt*(Lmean/Lj .* pack(phimatrix)). this is the linear map B.
function _nonlinearterm!(out::AbstractVector, sys::HBSystem)
    fill!(sys.AoLjbmvector, 0)
    AoLjbmvectorview = view(sys.AoLjbmvector, sys.Ljbm.nzind)
    phimatrixtovector!(AoLjbmvectorview, sys.phimatrix, sys.freqindexmap,
        sys.conjsourceindices, sys.conjtargetindices,
        length(sys.Ljb.nzval))
    _scalebranchvector!(AoLjbmvectorview, sys.Lmean, sys.Ljbm)
    mul!(out, sys.Rbnmt, sys.AoLjbmvector)
    return out
end

function _scalebranchvector!(v, Lmean, Ljbm::SparseVector)
    for i in eachindex(v)
        v[i] = v[i]*(Lmean/Ljbm.nzval[i])
    end
    return v
end

# out .+= invLnm*z + im*Gnm*wmodesm*z - Cnm*wmodes2m*z, the frequency
# dependent linear terms K applied to z, without allocations.
function _addlinearterm!(out::AbstractVector, sys::HBSystem,
    z::AbstractVector)
    mul!(sys.workn1, sys.wmodesm, z)
    mul!(sys.workn2, sys.Gnm, sys.workn1)
    @. out += im*sys.workn2
    mul!(sys.workn1, sys.wmodes2m, z)
    mul!(sys.workn2, sys.Cnm, sys.workn1)
    @. out -= sys.workn2
    mul!(sys.workn2, sys.invLnm, z)
    @. out += sys.workn2
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
    return _setpoint!(sys)
end

function setpoint!(sys::HBSystem, xr::AbstractVector{<:Real})
    real_to_complex!(sys.x, xr, sys.modelayout.isreal)
    return _setpoint!(sys)
end

function _setpoint!(sys::HBSystem)
    _branchtimedomain!(sys.phitd, sys, sys.x)
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
    _nonlinearterm!(F, sys)
    _addlinearterm!(F, sys, sys.x)
    @. F -= sys.bnm
    return F
end

function residual!(Fr::AbstractVector{<:Real}, sys::HBSystem)
    residual!(sys.workF, sys)
    complex_to_real!(Fr, sys.workF, sys.modelayout.isreal)
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
    _branchtimedomain!(sys.dirtd, sys, v)
    _multiplyintowork!(sys.worktd, sys.costd, sys.dirtd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    _nonlinearterm!(Jv, sys)
    _addlinearterm!(Jv, sys, v)
    return Jv
end

function jacobianvectorproduct!(Jvr::AbstractVector{<:Real},
    sys::HBSystem, vr::AbstractVector{<:Real})
    real_to_complex!(sys.workv, vr, sys.modelayout.isreal)
    jacobianvectorproduct!(sys.workF, sys, sys.workv)
    complex_to_real!(Jvr, sys.workF, sys.modelayout.isreal)
    return Jvr
end

# out .= a .* b and out .= -a .* b .* c behind function barriers
function _multiplyintowork!(out::Array{T}, a::Array{T},
    b::Array{T}) where T
    @inbounds for i in eachindex(out)
        out[i] = a[i]*b[i]
    end
    return out
end

function _multiplyintowork!(out::Array{T}, a::Array{T}, b::Array{T},
    c::Array{T}) where T
    @inbounds for i in eachindex(out)
        out[i] = -a[i]*b[i]*c[i]
    end
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
    _branchtimedomain!(sys.dirtd, sys, v)
    _branchtimedomain!(sys.dirtd2, sys, w)
    _multiplyintowork!(sys.worktd, sys.sintd, sys.dirtd, sys.dirtd2)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    _nonlinearterm!(Hvw, sys)
    return Hvw
end

function hessianvectorproduct!(Hvwr::AbstractVector{<:Real},
    sys::HBSystem, vr::AbstractVector{<:Real}, wr::AbstractVector{<:Real})
    real_to_complex!(sys.workv, vr, sys.modelayout.isreal)
    real_to_complex!(sys.workw, wr, sys.modelayout.isreal)
    hessianvectorproduct!(sys.workF, sys, sys.workv, sys.workw)
    complex_to_real!(Hvwr, sys.workF, sys.modelayout.isreal)
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
"""
function cosdirectionalderivative!(dcos::Array, sys::HBSystem,
    v::AbstractVector{<:Complex})
    _ensuresin!(sys)
    _branchtimedomain!(sys.dirtd, sys, v)
    _multiplyintowork!(sys.worktd, sys.sintd, sys.dirtd)
    applyfft!(dcos, sys.worktd, sys.rfftplan)
    @inbounds for i in eachindex(dcos)
        dcos[i] = -dcos[i]
    end
    return dcos
end

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
    assemblecomplexjacobian!(Jx, sys.complexjacobianplan, sys.phimatrix,
        sys.invLnm, sys.Gnm, sys.Cnm, sys.wmodesm, sys.wmodes2m)
    return Jx
end

function jacobian!(Jr::SparseMatrixCSC{<:Real}, sys::HBSystem)
    isnothing(sys.realjacobianplan) && throw(ArgumentError(
        "no real Jacobian plan was provided to this HBSystem."))
    _updatecosphimatrix!(sys)
    assemblerealjacobian!(Jr, sys.realjacobianplan, sys.phimatrix,
        sys.invLnm, sys.Gnm, sys.Cnm, sys.wmodesm, sys.wmodes2m)
    return Jr
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
    symfreqvar, wpumpmodes, Nnodes::Integer)

    # the sparsity structure must contain the modified nodal analysis
    # augmentation as well, so merge its entries into the numeric copies
    # used only for the structure. the source index stride is the size of
    # the pump frequency grid, the leading dimensions of phimatrix.
    invLnmpattern = spaddkeepzeros(spaddkeepzeros(invLnmcopy, Amna0), AmnaG)
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
        AmnaGfreqsubstindices, symfreqvar, wpumpmodes, Nmodes, Nnodes)
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
conjugate of the pump modulation contribution is used, which is the adjoint
system solved for the noise and quantum efficiency calculations. The
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
Nothing else contributes: the
auxiliary rows of the promoted coupled inductors are already symmetric (see
[`calcAmnaind`](@ref)); the linear term matrices are symmetric and mode
diagonal, so the column indexed frequency scaling and conjugation of
`sparseaddconjsubst!` are symmetric under transposition; and the transpose of
the pump modulation contribution exchanges the mode pair, mapping each
difference harmonic to its complex conjugate, which is the same as
conjugating the pump.

Because the source terms of [`hblinsolve`](@ref) are supported on the node
flux rows, `inv(D)` acts trivially on them and the solutions of the two
systems agree exactly in the node flux rows, which are the only rows the
noise, quantum efficiency, commutation relation, and adjoint node flux,
voltage and impedance outputs read. [`hblinsolve`](@ref) therefore obtains
its adjoint solutions with [`trysolvetranspose!`](@ref) on the factorization
of the forward system, rather than assembling and factorizing this matrix
with `conjugatepump = true` at every signal frequency. The conjugated pump
assembly is retained as the independent construction that equivalence is
tested against.
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

    return A
end

function assemblesystemmatrix!(A::SparseMatrixCSC,
    lsys::HBLinearizedSystem, ws::Number; conjugatepump::Bool = false)
    return assemblesystemmatrix!(A, lsys, ws .+ lsys.wpumpmodes;
        conjugatepump = conjugatepump)
end
