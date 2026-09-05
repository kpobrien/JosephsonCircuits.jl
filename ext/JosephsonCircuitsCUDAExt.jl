"""
    JosephsonCircuitsCUDAExt

Package extension loaded with `using CUDA`. It supplies the real transform
plans, through CUFFT, which the residual and the matrix-free
Jacobian-vector and Hessian-vector products of
[`JosephsonCircuits.HBSystem`](@ref) need on a CUDA device.

Everything else on that path is device generic: the linear maps around the
pointwise time domain nonlinearity are KernelAbstractions kernels of
`NonlinearTermPlan`, the pointwise and normalization steps are broadcasts,
and the transform is a batched real transform over all but the last
dimension, which is the layout CUFFT batches over. The Jacobian assembly
is written as host loops and is not addressed here; a device factorization
(see `CUDSSFactorization` and the CUDSS extension) is needed for the direct
solve path on a device.

"""
module JosephsonCircuitsCUDAExt

using CUDA
using CUDA.CUFFT
using KernelAbstractions
import JosephsonCircuits: fftplans, freememory, batchedinverse!, batchedmul!

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
    F::CuArray{T,3}, ::CUDABackend) where {T}
    copyto!(F, D)
    pivots, info = CUDA.CUBLAS.getrf_strided_batched!(F, true)
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
