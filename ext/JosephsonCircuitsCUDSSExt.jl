module JosephsonCircuitsCUDSSExt

using JosephsonCircuits
using CUDA
using CUDA.CUSPARSE
using CUDSS
using SparseArrays
using LinearAlgebra

# Replacing the values behind a solver was `cudss_set(solver, A)` in CUDSS
# 0.4 and became `cudss_update` in later versions. Bind whichever this CUDSS
# provides, so the extension works across both.
const setmatrix! = isdefined(CUDSS, :cudss_update) ? CUDSS.cudss_update :
    CUDSS.cudss_set

import JosephsonCircuits: _cudss_factorize, _cudss_factorize!,
    blockdiagonallayout, gatherblocks!, scatterblocks!, blockmatrix,
    BatchedBlockLayout, myldiv!, tobackend, cscvaluepermutation

# ---------------------------------------------------------------------------
# unbatched: one sparse system on the device, analysis reused across steps
# ---------------------------------------------------------------------------

mutable struct CUDSSSolve{TS,TM,TV}
    solver::TS
    A::TM              # device copy of the matrix, values overwritten in place
    x::TV
    b::TV
    # the device matrix is CSR and the host one CSC, so its stored values are
    # a permutation of theirs. the pattern never changes, so the permutation
    # is computed once here and the values are reordered on the host, where
    # they already are, into one contiguous transfer.
    perm::Vector{Int}
    vals::Vector{Float64}
end

# The batched path needs the block assignment, which a bare matrix does not
# carry, so it is selected by building the factorization from a layout
# (`cudssbatchedfactorization`) rather than by a flag here. `CUDSSFactorization`
# always takes the unbatched path; `batchedfactorization` below picks between
# them and is what callers should use.
function _cudss_factorize(A::SparseMatrixCSC{Float64,<:Integer};
    batched::Bool = true, kwargs...)
    n = size(A, 1)
    Agpu = CuSparseMatrixCSR(A)
    x = CUDA.zeros(Float64, n)
    b = CUDA.zeros(Float64, n)
    solver = CudssSolver(Agpu, "G", 'F')
    cudss("analysis", solver, x, b)
    cudss("factorization", solver, x, b)
    perm = cscvaluepermutation(A)
    return CUDSSSolve(solver, Agpu, x, b, perm,
        Vector{Float64}(undef, length(perm)))
end

function _cudss_factorize!(F::CUDSSSolve, A::SparseMatrixCSC; kwargs...)
    # the pattern is unchanged, so only the values move and only the numeric
    # phase is redone. they are reordered into the row major order of the
    # device matrix first: copying `nonzeros(A)` straight across would put
    # column major values into a row major array.
    nz = nonzeros(A)
    @inbounds for k in eachindex(F.perm)
        F.vals[k] = nz[F.perm[k]]
    end
    copyto!(F.A.nzVal, F.vals)
    setmatrix!(F.solver, F.A)
    cudss("refactorization", F.solver, F.x, F.b)
    return F
end

function myldiv!(x::AbstractVector, F::CUDSSSolve, b::AbstractVector)
    copyto!(F.b, b)
    cudss("solve", F.solver, F.x, F.b)
    copyto!(x, F.x)
    return x
end

# ---------------------------------------------------------------------------
# batched: the mode block diagonal as a uniform batch
# ---------------------------------------------------------------------------

mutable struct CUDSSBatchedSolve{TS,TL,TM,TV,TB}
    solver::TS
    # the layout with its slot map on the device, so the per solve gather and
    # scatter are kernels rather than scalar indexing of a device vector
    layout::TL
    blocks::Vector{TM}     # device blocks, one per mode, shared pattern
    xs::Vector{TV}
    bs::Vector{TV}
    rhs::TB                # device, blocksize x nblocks
    sol::TB
    # as in the unbatched case: the device blocks are CSR and the host ones
    # CSC. every block shares one pattern, so one permutation serves the batch.
    perm::Vector{Int}
    vals::Vector{Float64}
end

"""
    cudssbatchedfactorization(layout::BatchedBlockLayout; kwargs...)

A [`Factorization`](@ref) which treats the matrix as the uniform batch
described by `layout`: one symbolic analysis for the shared pattern, then a
batched numeric refactorization per Newton step. Build `layout` once with
[`blockdiagonallayout`](@ref); if it returns `nothing` the matrix is not a
uniform batch and [`CUDSSFactorization`](@ref) should be used instead.
"""
function JosephsonCircuits.cudssbatchedfactorization(layout::BatchedBlockLayout;
    kwargs...)
    return JosephsonCircuits.Factorization(
        (A; kws...) -> _batched_factorize(A, layout; kws...),
        (F, A; kws...) -> _batched_factorize!(F, A; kws...),
        kwargs)
end

function _batched_factorize(A::SparseMatrixCSC, layout::BatchedBlockLayout;
    kwargs...)

    nb, bs = layout.nblocks, layout.blocksize
    blocks = [CuSparseMatrixCSR(blockmatrix(A, layout, b)) for b in 1:nb]
    xs = [CUDA.zeros(Float64, bs) for _ in 1:nb]
    bsv = [CUDA.zeros(Float64, bs) for _ in 1:nb]
    solver = CudssBatchedSolver(blocks, "G", 'F')
    cudss("analysis", solver, xs, bsv)
    cudss("factorization", solver, xs, bsv)
    # the slot map moves to the device once; nzindex stays on the host, where
    # it indexes the host sparse matrix during a numeric refactorization
    devlayout = tobackend(CUDABackend(), layout)
    perm = cscvaluepermutation(blockmatrix(A, layout, 1))
    return CUDSSBatchedSolve(solver, devlayout, blocks, xs, bsv,
        CUDA.zeros(Float64, bs, nb), CUDA.zeros(Float64, bs, nb),
        perm, Vector{Float64}(undef, length(perm)))
end

function _batched_factorize!(F::CUDSSBatchedSolve, A::SparseMatrixCSC;
    kwargs...)
    # one permuted gather on the host, then values only to the device: the
    # shared pattern and the analysis are untouched
    vals = Matrix{Float64}(undef, size(F.layout.nzindex))
    # the host method, gathering from the host sparse matrix through nzindex
    gatherblocks!(vals, A, F.layout)
    for b in 1:F.layout.nblocks
        # reorder each block's values into the row major order of its device
        # copy before the transfer, as the unbatched path does
        @inbounds for k in eachindex(F.perm)
            F.vals[k] = vals[F.perm[k], b]
        end
        copyto!(F.blocks[b].nzVal, F.vals)
    end
    setmatrix!(F.solver, F.blocks)
    cudss("refactorization", F.solver, F.xs, F.bs)
    return F
end

# `b` and `x` are the vectors of the Krylov iteration, which are allocated
# like its right hand side and therefore live on the device. The gather and
# the scatter are kernels on that backend, and every copy below is device to
# device, so a preconditioner solve makes no host round trip.
function myldiv!(x::AbstractVector, F::CUDSSBatchedSolve, b::AbstractVector)
    gatherblocks!(F.rhs, b, F.layout)
    for k in 1:F.layout.nblocks
        copyto!(F.bs[k], view(F.rhs, :, k))
    end
    cudss("solve", F.solver, F.xs, F.bs)
    for k in 1:F.layout.nblocks
        copyto!(view(F.sol, :, k), F.xs[k])
    end
    scatterblocks!(x, F.sol, F.layout)
    return x
end

end # module
