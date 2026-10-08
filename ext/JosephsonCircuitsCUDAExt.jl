"""
    JosephsonCircuitsCUDAExt

Package extension loaded with `using CUDA`. It supplies four things: the
real transform plans, through CUFFT, which the residual and the matrix-free
Jacobian-vector and Hessian-vector products of
[`JosephsonCircuits.HBSystem`](@ref) need on a CUDA device; the device's
free memory, which the automatic choices of preconditioner and linearized
factorization are sized against; the batched dense primitives
`batchedinverse!` and `batchedmul!`, cuBLAS `getrf`/`getri` strided
batched and `gemm_strided_batched!`, which every dense operation of a
`SparseBlockFactorization` is one call to; and the period map's dense
eigensolve, cuSOLVER's `geev` with the left vectors it needs checked
(see `mapspectrum!`).

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
    blockidentity!, transientiqfftplans, mapspectrum!
using JosephsonCircuits: eigenvector

# Real transform plans on the device with the same dimensions, direction
# and normalization convention as the FFTW plans of the CPU backend: the
# transform runs over all but the last dimension, the inverse is the
# unnormalized backward transform, and the caller scales the forward one by
# `1/prod(size(td)[1:end-1])`.
function fftplans(fd::AbstractArray{Complex{T}}, td::AbstractArray{T},
    stepsperperiod::Int, backend::CUDABackend) where T
    dims = 1:length(size(fd))-1
    irfftplan = CUFFT.plan_brfft(fd, stepsperperiod, dims)
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


# the complex in place plans of the transient's windowed I/Q measurement
transientiqfftplans(work, ::CUDABackend) = (CUFFT.plan_fft!(work), CUFFT.plan_bfft!(work))

# The multipliers of the balanced period map `M` on the device (see
# `JosephsonCircuits.mapspectrum!`), which it leaves as it is: cuSOLVER's
# geev of every right vector, and the left vectors of the chosen ones from
# the inverse of the right ones, by one factorization: the row of a real
# multiplier's column is its left vector, and the rows of a complex pair's
# two columns hold the first's left vector as LAPACK holds it, twice over.
# The inverse holds a left vector to working accuracy only where the right
# vectors are well conditioned together, which a defective multiplier
# anywhere in the map undoes, so each is taken where its backward error,
# its residual in the map's transpose over the map's 1-norm and its own
# norm, is within `tolerance`, `n` eps, the bound of a backward stable
# vector; where the right vectors are singular or a left vector's error
# exceeds it, `vectors` gives nothing, and the host's take their place.
function mapspectrum!(M::Matrix{Float64}, ilo::Int, ihi::Int, ::CUDABackend;
        tolerance::Float64 = size(M, 1)*eps())
    n = size(M, 1)
    B = CuArray(M)
    W, _, V = CUDA.CUSOLVER.Xgeev!('N', 'V', copy(B))
    values = Array(W)
    wi = imag.(values)
    norm1 = LinearAlgebra.opnorm(M, 1)
    factors = Ref{Any}(nothing)
    function vectors(ks)
        firsts = sort!(unique!([wi[k] < 0 ? k - 1 : k for k in ks]))
        at, columns, column, c = Pair{Int,Int}[], Int[], zeros(Int, n), 1
        for k in firsts
            push!(at, k => c); push!(columns, k); column[k] = c
            wi[k] == 0 || (push!(at, k + 1 => c + 1); push!(columns, k + 1); column[k + 1] = c + 1)
            c += wi[k] == 0 ? 1 : 2
        end
        if isnothing(factors[])
            F, ipiv, info = CUDA.CUSOLVER.getrf!(copy(V))
            info == 0 || return nothing
            factors[] = (F, ipiv)
        end
        F, ipiv = factors[]
        E = zeros(n, length(columns))
        for (j, k) in enumerate(columns)
            E[k, j] = 1.0
        end
        L = CUDA.CUSOLVER.getrs!('T', F, ipiv, CuArray(E))
        left, image = Array(L), Array(transpose(B)*L)
        for k in firsts
            u = eigenvector(left, k, column, wi)
            LinearAlgebra.norm(eigenvector(image, k, column, wi) .- conj(values[k]) .* u) <=
                tolerance*norm1*LinearAlgebra.norm(u) || return nothing
        end
        return Array(V[:, columns]), left, at
    end
    return values, vectors
end

end # module
