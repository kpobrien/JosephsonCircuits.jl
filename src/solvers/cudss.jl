
"""
    cscvaluepermutation(A::SparseMatrixCSC)

The permutation `p` for which `nonzeros(A)[p]` is the stored-value order of the
*compressed sparse row* form of `A`.

A device sparse matrix is CSR, and building one from a host `SparseMatrixCSC`
converts the layout, which permutes the stored values. A numeric
refactorization which reuses the symbolic analysis must therefore push new
values through the same permutation rather than copying `nonzeros` straight
into the device array: the pattern is unchanged, but the order is not.

The permutation depends only on the sparsity pattern, so it is computed once
alongside the analysis and reused for every subsequent refactorization.
"""
function cscvaluepermutation(A::SparseMatrixCSC)
    # a matrix with the same pattern whose values are the positions of
    # `nonzeros(A)`; transposing it to CSC of transpose(A), which is CSR of A,
    # carries each position to its CSR slot
    positions = SparseMatrixCSC(size(A, 1), size(A, 2),
        copy(SparseArrays.getcolptr(A)), copy(rowvals(A)),
        collect(1:nnz(A)))
    return nonzeros(sparse(transpose(positions)))
end

"""
    CUDSSFactorization(; kwargs...)

An [`AbstractFactorization`](@ref) backed by NVIDIA's cuDSS direct sparse
solver, for use on a GPU.

Requires `CUDSS.jl` and `CUDA.jl` to be loaded; without them the returned
factorization raises an informative error when used, so the constructor itself
is always available and a script can select it unconditionally. Nothing else in
the solver changes: `AbstractFactorization` already separates `factorize` (the
symbolic analysis and the first numeric factorization) from `refactorize!` (a
numeric refactorization reusing the analysis), which is exactly the split cuDSS
wants.

That split is what makes this worth doing. The symbolic analysis depends only
on the circuit topology and the retained mode coupling, neither of which
changes across Newton steps, while the values change at every step; the
analysis is almost all of the cost of a factorization and is paid once.

cuDSS is handed the whole matrix, block diagonal or not: it discovers the
independent blocks of a mode block diagonal itself and works on them
together, and measured faster that way than when handed the blocks as a
uniform batch.
"""
struct CUDSSFactorization{K} <: AbstractFactorization
    kwargs::K
end
CUDSSFactorization(; kwargs...) = CUDSSFactorization(kwargs)
factorize(f::CUDSSFactorization, A) = _cudss_factorize(A; f.kwargs...)
refactorize!(f::CUDSSFactorization, F, A) = _cudss_factorize!(F, A; f.kwargs...)

# Overridden by the CUDSS extension. The error names the missing packages
# rather than failing with a MethodError somewhere inside the solver.
function _cudss_factorize(A; kwargs...)
    throw(ArgumentError(
        "CUDSSFactorization requires CUDSS.jl and CUDA.jl to be loaded. Run `using CUDA, CUDSS` before calling the solver, or use the default KLUfactorization()."))
end

function _cudss_factorize!(F, A; kwargs...)
    throw(ArgumentError(
        "CUDSSFactorization requires CUDSS.jl and CUDA.jl to be loaded."))
end

"""
    uniformbatchlimit(nrhs::Integer)

The largest uniform batch of systems which cuDSS solves correctly with `nrhs`
right hand sides each.

!!! warning "This works around a wrong answer, not a failure"
    cuDSS 0.7 (through CUDSS.jl 0.8.0) returns silently wrong solutions from a
    uniform batch of sixteen or more systems once each has six or more right
    hand sides. Every system of the batch comes back wrong, by order one, while
    `cudss_get(solver, "info")` reports success and `"lu_nnz"` is unchanged, so
    nothing downstream can detect it. A batch of fifteen is correct with
    twelve right hand sides and a batch of sixteen is wrong by order one;
    with one or two right hand sides batches of well over a hundred are
    correct. The cap costs nothing on this path, since the speedup of
    batching a frequency sweep through cuDSS saturates by about a dozen
    systems; it applies only to the cuDSS batch, a
    [`SparseBlockFactorization`](@ref) sweep sizes its batch by memory
    instead ([`blocksystembytes`](@ref)) and profits from batches well past
    this. Re-check against newer cuDSS releases before raising it.
"""
uniformbatchlimit(nrhs::Integer) = nrhs >= 6 ? 15 : 128

# Overridden by the CUDSS extension: a uniform batch of systems sharing one
# sparsity pattern, analyzed once and then refactorized and solved as a batch.
function _cudss_sweep(rowptr, colind, nzval, X, B; kwargs...)
    throw(ArgumentError(
        "solving a frequency sweep on a device requires CUDSS.jl and CUDA.jl to be loaded. Run `using CUDA, CUDSS` before calling the solver, or leave `backend` at its default."))
end

function _cudss_sweepsolve!(S)
    throw(ArgumentError(
        "solving a frequency sweep on a device requires CUDSS.jl and CUDA.jl to be loaded."))
end
