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
    _cudss_sweep, _cudss_sweepsolve!,
    blockdiagonallayout, gatherblocks!, scatterblocks!, gathervalues!,
    blockpattern, blockrowmajorindex, blockmatrix, BatchedBlockLayout,
    myldiv!, tobackend, cscvaluepermutation, rowpointer, columnindices

# ---------------------------------------------------------------------------
# binding the caller's vectors to a cuDSS descriptor
# ---------------------------------------------------------------------------

# cuDSS reads the right hand side and writes the solution through a descriptor
# which carries a device pointer, and `cudss(phase, solver, x::CuVector,
# b::CuVector)` builds a fresh pair of descriptors on every call. Building them
# once and rebinding them to the caller's vectors removes both the construction
# and the staging copies that existed only because the descriptor was tied to
# an owned buffer.
#
# Rebinding needs the vector to be a device vector of the right element type
# and length. The right hand sides of the Krylov iteration are columns of the
# Arnoldi basis; CUDA.jl represents a contiguous view of a `CuMatrix` as a
# `CuArray`, so they qualify, but a caller with anything else still gets a
# correct answer through the staging path below.
bindable(v, ::Type{T}, n::Integer) where {T} = v isa CuVector{T} && length(v) == n

# ---------------------------------------------------------------------------
# unbatched: one sparse system on the device, analysis reused across steps
# ---------------------------------------------------------------------------

mutable struct CUDSSSolve{TS,TM,TV,TD,Tv<:Union{AbstractFloat,Complex}}
    solver::TS
    A::TM              # device copy of the matrix, values overwritten in place
    x::TV
    b::TV
    # the solution and right hand side descriptors, created once and rebound
    # to the caller's vectors on each solve
    xdesc::TD
    bdesc::TD
    # the device matrix is CSR and the host one CSC, so its stored values are
    # a permutation of theirs. the pattern never changes, so the permutation
    # is computed once here and the values are reordered on the host, where
    # they already are, into one contiguous transfer.
    perm::Vector{Int}
    # the staging buffer for that reorder, in the element type of the matrix:
    # cuDSS factorizes in single as well as in double precision, and the
    # working precision of the solve is what the caller handed in.
    vals::Vector{Tv}
end

# The batched path needs the block assignment, which a bare matrix does not
# carry, so it is selected by building the factorization from a layout
# (`cudssbatchedfactorization`) rather than by a flag here. `CUDSSFactorization`
# always takes the unbatched path; `batchedfactorization` picks between them
# from the sparsity structure and is what callers should use.
function _cudss_factorize(A::SparseMatrixCSC{Tv,<:Integer};
    kwargs...) where {Tv<:AbstractFloat}
    n = size(A, 1)
    Agpu = CuSparseMatrixCSR(A)
    F = newsolver(Agpu, Tv, n, cscvaluepermutation(A))
    return F
end

# ---------------------------------------------------------------------------
# values already on the device, in the row major order the device matrix
# stores: nothing crosses to the host at all
# ---------------------------------------------------------------------------

# The permutation and the host staging buffer of the SparseMatrixCSC methods
# exist to reorder column major host values into row major device ones. A
# DeviceValuedSparseMatrix already carries its values in that order, on the
# device, so both are empty here and a refactorization is one device to device
# copy. The host `pattern` is read once, for its sparsity structure, and its
# stored values are never read.
function _cudss_factorize(A::JosephsonCircuits.DeviceValuedSparseMatrix{Tv};
    kwargs...) where {Tv<:Union{AbstractFloat,Complex}}
    n = size(A, 1)
    # the structure is already a row pointer and a column index array, so the
    # device matrix is built from it directly: no conversion, and nothing of
    # the values crosses to the host
    nzval = CUDA.similar(A.nzval)
    copyto!(nzval, A.nzval)
    Agpu = CuSparseMatrixCSR{Tv,Int32}(
        CuVector(convert(Vector{Int32}, rowpointer(A))),
        CuVector(convert(Vector{Int32}, columnindices(A))), nzval, size(A))
    return newsolver(Agpu, Tv, n)
end

# the analysis and the first numeric factorization, shared by both entry
# points. The descriptors start bound to the owned buffers, which is what the
# analysis and every later refactorization use.
function newsolver(Agpu, ::Type{Tv}, n::Integer,
    perm::Vector{Int} = Int[]) where {Tv}
    x = CUDA.zeros(Tv, n)
    b = CUDA.zeros(Tv, n)
    solver = CudssSolver(Agpu, "G", 'F')
    xdesc = CudssMatrix(x)
    bdesc = CudssMatrix(b)
    cudss("analysis", solver, xdesc, bdesc)
    cudss("factorization", solver, xdesc, bdesc)
    return CUDSSSolve(solver, Agpu, x, b, xdesc, bdesc, perm,
        Vector{Tv}(undef, length(perm)))
end

function _cudss_factorize!(F::CUDSSSolve,
    A::JosephsonCircuits.DeviceValuedSparseMatrix; kwargs...)
    copyto!(F.A.nzVal, A.nzval)
    return refactorize!(F)
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
    return refactorize!(F)
end

# the refactorization phase reads whatever the descriptors point at, and the
# last solve left them bound to the caller's vectors, which may since have
# been freed. Rebind them to the owned buffers first.
function refactorize!(F::CUDSSSolve)
    setmatrix!(F.solver, F.A)
    CUDSS.cudss_update(F.xdesc, F.x)
    CUDSS.cudss_update(F.bdesc, F.b)
    cudss("refactorization", F.solver, F.xdesc, F.bdesc)
    return F
end

function myldiv!(x::AbstractVector, F::CUDSSSolve, b::AbstractVector)
    T, n = eltype(F.x), length(F.x)
    if bindable(x, T, n) && bindable(b, T, n)
        CUDSS.cudss_update(F.xdesc, x)
        CUDSS.cudss_update(F.bdesc, b)
        cudss("solve", F.solver, F.xdesc, F.bdesc)
    else
        copyto!(F.b, b)
        CUDSS.cudss_update(F.xdesc, F.x)
        CUDSS.cudss_update(F.bdesc, F.b)
        cudss("solve", F.solver, F.xdesc, F.bdesc)
        copyto!(x, F.x)
    end
    return x
end

# ---------------------------------------------------------------------------
# batched: the mode block diagonal as a uniform batch
# ---------------------------------------------------------------------------

# cuDSS calls a set of systems which share one sparsity pattern a *uniform
# batch*, and represents it as one row pointer, one column index array, and a
# `blocknnz` by `nblocks` matrix of values; the right hand sides and solutions
# are likewise one `blocksize` by `nblocks` matrix each. Every array here is
# therefore a single contiguous device buffer, which is what removes the per
# block staging copies an array-of-arrays interface would need: the gather
# writes the batch's right hand sides directly and the assembly writes the
# batch's values directly.
mutable struct CUDSSBatchedSolve{TS,TL,TM,TV,TD,TI,Tv<:AbstractFloat}
    solver::TS
    # the layout with its slot map on the device, so the per solve gather and
    # scatter are kernels rather than scalar indexing of a device vector
    layout::TL
    rowPtr::TV             # the shared pattern, one copy for the whole batch
    colVal::TV
    nzVal::TM              # blocknnz x nblocks
    rhs::TM                # blocksize x nblocks
    sol::TM
    xdesc::TD              # bound to sol and rhs once, at construction
    bdesc::TD
    # a stored entry of the batch, in the row major order the device wants,
    # back to where its value is found. `src` indexes the device value array
    # of a DeviceValuedSparseMatrix; `hostsrc` indexes `nonzeros` of a host
    # SparseMatrixCSC. Only one of the two is ever populated.
    src::TI
    hostsrc::Matrix{Int}
    hostvals::Matrix{Tv}
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

# the shared pattern as compressed sparse rows, which is compressed sparse
# columns of its transpose. Built on the host from the layout rather than by
# converting a device matrix, so nothing depends on how a conversion treats a
# stored zero.
function sharedrowmajorpattern(layout::BatchedBlockLayout)
    blkT = sparse(transpose(blockpattern(layout)))
    return CuVector{Cint}(SparseArrays.getcolptr(blkT)),
        CuVector{Cint}(rowvals(blkT))
end

# the analysis, the first numeric factorization, and the descriptors, shared by
# the host valued and device valued entry points. `fill!` is the caller's job:
# `nzVal` already holds the values of the first operating point.
function newbatchedsolver(layout::BatchedBlockLayout, rowPtr, colVal, nzVal,
    ::Type{Tv}) where {Tv}
    nb, bs = layout.nblocks, layout.blocksize
    solver = CudssSolver(rowPtr, colVal, nzVal, "G", 'F')
    cudss_set(solver, "ubatch_size", nb)
    rhs = CUDA.zeros(Tv, bs, nb)
    sol = CUDA.zeros(Tv, bs, nb)
    # a batched dense descriptor is created for the shape and then pointed at
    # the buffers, which never move, so this is the only binding either needs
    xdesc = CudssMatrix(Tv, bs; nbatch = nb)
    bdesc = CudssMatrix(Tv, bs; nbatch = nb)
    CUDSS.cudss_update(xdesc, sol)
    CUDSS.cudss_update(bdesc, rhs)
    cudss("analysis", solver, xdesc, bdesc)
    cudss("factorization", solver, xdesc, bdesc)
    # the slot map moves to the device once, because the right hand side
    # gather and the solution scatter are kernels on the Krylov vectors
    devlayout = tobackend(CUDABackend(), layout)
    return solver, devlayout, rhs, sol, xdesc, bdesc
end

function _batched_factorize(A::SparseMatrixCSC{Tv}, layout::BatchedBlockLayout;
    kwargs...) where {Tv}

    rowPtr, colVal = sharedrowmajorpattern(layout)
    hostsrc = blockrowmajorindex(layout)
    hostvals = Matrix{Tv}(undef, size(hostsrc))
    nz = nonzeros(A)
    @inbounds for k in eachindex(hostsrc)
        hostvals[k] = nz[hostsrc[k]]
    end
    nzVal = CuMatrix{Tv}(hostvals)
    solver, devlayout, rhs, sol, xdesc, bdesc =
        newbatchedsolver(layout, rowPtr, colVal, nzVal, Tv)
    return CUDSSBatchedSolve(solver, devlayout, rowPtr, colVal, nzVal, rhs,
        sol, xdesc, bdesc, nothing, hostsrc, hostvals)
end

# The values are already on the device, in the row major order of the *full*
# matrix. Each stored entry of the batch is one of them, so the whole numeric
# update is a single indexed copy on the device: the composition of "which
# entry of the full matrix is this block slot" with "where does the device
# value array keep that entry" is a fixed map, built once here.
function _batched_factorize(A::JosephsonCircuits.DeviceValuedSparseMatrix{Tv},
    layout::BatchedBlockLayout; kwargs...) where {Tv<:AbstractFloat}

    rowPtr, colVal = sharedrowmajorpattern(layout)
    hostsrc = blockrowmajorindex(layout)
    # `nonzeros(A)` is in the row major order of the matrix, and the structure
    # is stored as the transpose, whose row major order is the matrix's column
    # major order. So transposing the stored structure's value order takes a
    # column major index straight to the device slot holding it.
    csrslot = cscvaluepermutation(A.patterntranspose)
    src = CuMatrix{Cint}(Cint[csrslot[i] for i in hostsrc])
    nzVal = CuMatrix{Tv}(undef, size(src))
    gathervalues!(nzVal, A.nzval, src)
    solver, devlayout, rhs, sol, xdesc, bdesc =
        newbatchedsolver(layout, rowPtr, colVal, nzVal, Tv)
    return CUDSSBatchedSolve(solver, devlayout, rowPtr, colVal, nzVal, rhs,
        sol, xdesc, bdesc, src, Matrix{Int}(undef, 0, 0),
        Matrix{Tv}(undef, 0, 0))
end

function _batched_factorize!(F::CUDSSBatchedSolve,
    A::JosephsonCircuits.DeviceValuedSparseMatrix; kwargs...)
    gathervalues!(F.nzVal, A.nzval, F.src)
    return refactorize!(F)
end

function _batched_factorize!(F::CUDSSBatchedSolve, A::SparseMatrixCSC;
    kwargs...)
    # one permuted gather on the host, then values only to the device: the
    # shared pattern and the analysis are untouched
    nz = nonzeros(A)
    @inbounds for k in eachindex(F.hostsrc)
        F.hostvals[k] = nz[F.hostsrc[k]]
    end
    copyto!(F.nzVal, F.hostvals)
    return refactorize!(F)
end

# `nzVal` is overwritten in place, so the solver's matrix descriptor already
# points at the new values and only the numeric phase is redone.
function refactorize!(F::CUDSSBatchedSolve)
    cudss("refactorization", F.solver, F.xdesc, F.bdesc)
    return F
end

# `b` and `x` are the vectors of the Krylov iteration, which are allocated
# like its right hand side and therefore live on the device. The gather and
# the scatter are kernels on that backend, and the batch's right hand sides
# and solutions are the contiguous buffers the descriptors are bound to, so a
# preconditioner solve is two kernels and a cuDSS call with no copies at all.
function myldiv!(x::AbstractVector, F::CUDSSBatchedSolve, b::AbstractVector)
    gatherblocks!(F.rhs, b, F.layout)
    cudss("solve", F.solver, F.xdesc, F.bdesc)
    scatterblocks!(x, F.sol, F.layout)
    return x
end

# ---------------------------------------------------------------------------
# a frequency sweep as a uniform batch
# ---------------------------------------------------------------------------

# The linearized system matrices of a frequency sweep differ only in their
# stored values, so the batch shares one row pointer, one column index array
# and one symbolic analysis, and holds an `nnz` by `nbatch` matrix of values.
# The right hand sides are the same for every frequency but cuDSS wants one
# per system, so they are replicated into an `n` by `nrhs` by `nbatch` array;
# that costs a few tens of megabytes and saves a descriptor rebind per system.
mutable struct CUDSSSweep{T,INT,TS,TD}
    solver::TS
    rowptr::CuVector{INT}
    colind::CuVector{INT}
    nzval::CuMatrix{T}
    X::CuArray{T,3}
    B::CuArray{T,3}
    xdesc::TD
    bdesc::TD
    nbatch::Int
end

function _cudss_sweep(rowptr::CuVector{INT}, colind::CuVector{INT},
    nzval::CuMatrix{T}, X::CuArray{T,3}, B::CuArray{T,3};
    kwargs...) where {T<:Union{AbstractFloat,Complex},INT}

    n = length(rowptr) - 1
    nrhs, nbatch = size(X, 2), size(X, 3)
    size(nzval, 2) == nbatch || throw(DimensionMismatch(
        "the value matrix and the solution array must agree on the batch size."))
    solver = CudssSolver(rowptr, colind, vec(nzval), "G", 'F')
    for (k, v) in kwargs
        cudss_set(solver, string(k), v)
    end
    # a batch of one is not a batch: cuDSS rejects the batched descriptors,
    # so bind the single system's matrices directly
    xdesc, bdesc = if nbatch > 1
        cudss_set(solver, "ubatch_size", nbatch)
        xd = CudssMatrix(T, n, nrhs; nbatch = nbatch)
        bd = CudssMatrix(T, n, nrhs; nbatch = nbatch)
        CUDSS.cudss_update(xd, X); CUDSS.cudss_update(bd, B)
        xd, bd
    else
        CudssMatrix(reshape(X, n, nrhs)), CudssMatrix(reshape(B, n, nrhs))
    end
    cudss("analysis", solver, xdesc, bdesc)
    cudss("factorization", solver, xdesc, bdesc)
    return CUDSSSweep{T,INT,typeof(solver),typeof(xdesc)}(
        solver, rowptr, colind, nzval, X, B, xdesc, bdesc, nbatch)
end

# refactorize the whole batch against whatever values `S.nzval` now holds and
# solve every system for every right hand side, reusing the one analysis
function _cudss_sweepsolve!(S::CUDSSSweep)
    setmatrix!(S.solver, S.rowptr, S.colind, vec(S.nzval))
    cudss("refactorization", S.solver, S.xdesc, S.bdesc)
    cudss("solve", S.solver, S.xdesc, S.bdesc)
    CUDA.synchronize()
    return S.X
end

end # module
