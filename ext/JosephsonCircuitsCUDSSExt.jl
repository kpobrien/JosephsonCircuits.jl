"""
    JosephsonCircuitsCUDSSExt

Package extension loaded with `using CUDSS`. It supplies the cuDSS
factorizations of the device path: `_cudss_factorize` and
`_cudss_factorize!` for one sparse system, and `_cudss_sweep` and
`_cudss_sweepsolve!` for the uniform batch of a frequency sweep. cuDSS 0.8
has no transposed solve, which is why the sweep's adjoint direction is a
second factorization there; a `BlockFactorization` on the device does not
need this extension.
"""
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

# one sparse system: the analysis, the first numeric factorization and the
# descriptors bound to owned buffers
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
    # An ideal short scattering block has S = -1, so its constitutive
    # coefficient C = R^(1/2)(I + S) is exactly zero: the block's auxiliary
    # port current appears with a zero diagonal, a pure constraint row.
    # Host KLU pivots through that structure; cuDSS with its defaults
    # produces a factorization whose preconditioned residual grows by
    # orders of magnitude, which stalls the Newton-Krylov solve. Perturbing
    # zero pivots and cleaning up with two iterative refinement steps makes
    # the device factorization follow the host path exactly (same iteration
    # count, same final residual). The perturbation is below the O(1)
    # scaled Jacobian entries, so well pivoted systems are unaffected.
    CUDSS.cudss_set(solver, "pivot_epsilon", 1e-8)
    CUDSS.cudss_set(solver, "ir_n_steps", 2)
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
