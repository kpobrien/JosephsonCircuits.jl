
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
    CUDSSFactorization(; batched = false, kwargs...)

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

`batched = true` offers the mode block diagonal to cuDSS as a *uniform batch*
instead: one analysis for the pattern every mode block shares, then a batched
numeric refactorization. See [`BatchedBlockLayout`](@ref) for the structure and
[`blockdiagonallayout`](@ref) for the conditions, which
[`batchedfactorization`](@ref) checks; when they do not hold, the unbatched
path is used anyway.

!!! note "Why the batch is off by default"
    With the cuDSS versions tested it loses. cuDSS does not need to be told
    that a block diagonal matrix is a batch: given the whole matrix it
    discovers the independent blocks and works on them together, while given
    the blocks as a uniform batch it works through them far more slowly per
    solve and per refactorization, with a cost which jumps beyond about a
    dozen blocks. The batch only wins the symbolic analysis, which is paid
    once. It is kept because it is cheap to re-measure against a later
    cuDSS, and because the saving it does find is in the reordering, which
    the unbatched path could have too by analyzing one block and handing
    the replicated permutation to the whole matrix as `"user_perm"`.
"""
struct CUDSSFactorization{K} <: AbstractFactorization
    batched::Bool
    kwargs::K
end
CUDSSFactorization(; batched::Bool = false, kwargs...) =
    CUDSSFactorization(batched, kwargs)
factorize(f::CUDSSFactorization, A) = _cudss_factorize(A; f.kwargs...)
refactorize!(f::CUDSSFactorization, F, A) = _cudss_factorize!(F, A; f.kwargs...)

"""
    CUDSSBatchedFactorization(layout::BatchedBlockLayout; kwargs...)

The batched form of [`CUDSSFactorization`](@ref): the matrix treated as the
uniform batch of diagonal blocks described by `layout`, one symbolic
analysis for the shared pattern, then a batched numeric refactorization per
Newton step. Selected by [`batchedfactorization`](@ref) when the matrix is
such a batch and the factorization asked for it.
"""
struct CUDSSBatchedFactorization{L,K} <: AbstractFactorization
    layout::L
    kwargs::K
end
CUDSSBatchedFactorization(layout; kwargs...) =
    CUDSSBatchedFactorization(layout, kwargs)
factorize(f::CUDSSBatchedFactorization, A) =
    _batched_factorize(A, f.layout; f.kwargs...)
refactorize!(f::CUDSSBatchedFactorization, F, A) =
    _batched_factorize!(F, A; f.kwargs...)
batchedform(f::CUDSSFactorization, layout) =
    f.batched ? CUDSSBatchedFactorization(layout, f.kwargs) : f

# Overridden by the CUDSS extension. The error names the missing packages
# rather than failing with a MethodError somewhere inside the solver.
function _cudss_factorize(A; kwargs...)
    throw(ArgumentError(
        "CUDSSFactorization requires CUDSS.jl and CUDA.jl to be loaded. Run `using CUDA, CUDSS` before calling the solver, or use the default KLUfactorization()."))
end
function _batched_factorize(A, layout; kwargs...)
    throw(ArgumentError(
        "CUDSSFactorization requires CUDSS.jl and CUDA.jl to be loaded. Run `using CUDA, CUDSS` before calling the solver, or use the default KLUfactorization()."))
end
function _batched_factorize!(F, A; kwargs...)
    throw(ArgumentError(
        "CUDSSFactorization requires CUDSS.jl and CUDA.jl to be loaded."))
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

"""
    batchedfactorization(A::SparseMatrixCSC, blockof, nblocks,
        gpu, cpu = KLUfactorization())

Pick the factorization to use for `A`, preferring a batched GPU one when the
matrix really is a uniform batch of diagonal blocks.

The three-way choice is the whole point of the architecture:

1. `A` is block diagonal with uniform blocks, and a GPU factorization was
   supplied: use [`CUDSSBatchedFactorization`](@ref) on the shared pattern, so
   one symbolic analysis serves every block and every Newton step.
2. `A` is not a uniform batch (mode coupling is retained after an escalation,
   or the truncation contains a self-conjugate mode such as dc): use `gpu`
   unbatched.
3. No GPU factorization was supplied, or the packages are not loaded: use
   `cpu`, which is the ordinary KLU path and the default everywhere.

Case 3 is not a degraded mode. It is what the solver does today, and the GPU
path is an opt-in on top of it.
"""
function batchedfactorization(A::SparseMatrixCSC, blockof, nblocks,
    gpu::Union{AbstractFactorization,Nothing},
    cpu::AbstractFactorization = KLUfactorization())

    isnothing(gpu) && return cpu
    # a factorization with no batched form is used unchanged, and asking for
    # the layout would be a waste: it is a walk over every stored entry.
    (gpu isa CUDSSFactorization && gpu.batched) || return gpu
    layout = blockdiagonallayout(A, blockof, nblocks)
    isnothing(layout) && return gpu
    return batchedform(gpu, layout)
end
