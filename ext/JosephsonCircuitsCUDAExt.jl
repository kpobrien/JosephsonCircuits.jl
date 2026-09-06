"""
    JosephsonCircuitsCUDAExt

Package extension loaded with `using CUDA`. It supplies three things: the
real transform plans, through CUFFT, which the residual and the matrix-free
Jacobian-vector and Hessian-vector products of
[`JosephsonCircuits.HBSystem`](@ref) need on a CUDA device; the device's
free memory, which the automatic choices of preconditioner and linearized
factorization are sized against; and the batched dense primitives
`batchedinverse!` and `batchedmul!`, cuBLAS `getrf`/`getri` strided
batched and `gemm_strided_batched!`, which every dense operation of a
`SparseBlockFactorization` is one call to.

Everything else on that path is device generic: the linear maps around the
pointwise time domain nonlinearity are KernelAbstractions kernels of
`NonlinearTermPlan`, the pointwise and normalization steps are broadcasts,
the transform is a batched real transform over all but the last dimension,
which is the layout CUFFT batches over, and the Jacobian assembly is a
KernelAbstractions kernel per stored entry. For the direct solves a device
needs either a device sparse factorization (see `CUDSSFactorization` and
the CUDSS extension) or a `BlockFactorization`, whose dense kernels run on
the primitives here.

"""
module JosephsonCircuitsCUDAExt

using CUDA
using CUDA.CUFFT
using KernelAbstractions
import LinearAlgebra
using LinearAlgebra: lu!, ldiv!
import JosephsonCircuits: fftplans, freememory, batchedinverse!, batchedmul!,
    blockidentity!

# Real transform plans on the device with the same dimensions, direction
# and normalization convention as the FFTW plans of the CPU backend: the
# transform runs over all but the last dimension and the caller applies the
# `prod(size(td)[1:end-1])` normalization.
function fftplans(fd::AbstractArray{Complex{T}}, td::AbstractArray{T},
    stepsperperiod::Int, backend::CUDABackend) where T
    dims = 1:length(size(fd))-1
    irfftplan = CUFFT.plan_irfft(fd, stepsperperiod, dims)
    rfftplan = CUFFT.plan_rfft(td, dims)
    return irfftplan, rfftplan
end


freememory(::CUDABackend) = Int(CUDA.free_memory())

# the batched dense primitives of the block factorization of the
# linearized system: one cuBLAS call over the batch of systems
function batchedinverse!(Dinv::CuArray{T,3}, D::CuArray{T,3},
    F::CuArray{T,3}, backend::CUDABackend) where {T}
    if size(D, 3) == 1
        # a batch of one (the preconditioner's clusters, whose supernodes
        # are amalgamated to hundreds of rows): the dense LU of cuSOLVER,
        # which the batched routines below, made for many small blocks,
        # are far slower than at that size
        n = size(D, 1)
        Fm = reshape(F, n, n)
        copyto!(Fm, reshape(D, n, n))
        LU = lu!(Fm)
        blockidentity!(Dinv, backend)
        ldiv!(LU, reshape(Dinv, n, n))
        return Dinv
    end
    copyto!(F, D)
    pivots, info = CUDA.CUBLAS.getrf_strided_batched!(F, true)
    # a zero pivot in any system of the batch is a singular diagonal block;
    # the host path throws the same from `lu!`, and `tryfactorize!` then
    # refactorizes afresh rather than solve with garbage factors
    k = findfirst(!=(0), Array(info))
    isnothing(k) || throw(LinearAlgebra.SingularException(Int(k)))
    CUDA.CUBLAS.getri_strided_batched!(F, Dinv, pivots)
    return Dinv
end
function batchedmul!(C::AbstractArray{T,3}, A::AbstractArray{T,3},
    B::AbstractArray{T,3}, alpha, beta, tA::Bool, tB::Bool,
    ::CUDABackend) where {T}
    CUDA.CUBLAS.gemm_strided_batched!(tA ? 'T' : 'N', tB ? 'T' : 'N', T(alpha),
        A, B, T(beta), C)
    return C
end

end # module
