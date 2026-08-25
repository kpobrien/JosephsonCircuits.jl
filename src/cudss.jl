
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

A [`Factorization`](@ref) backed by NVIDIA's cuDSS direct sparse solver, for
use on a GPU.

Requires `CUDSS.jl` and `CUDA.jl` to be loaded; without them the returned
factorization raises an informative error when used, so the constructor itself
is always available and a script can select it unconditionally. Nothing else in
the solver changes: `Factorization` already separates `factorize` (the symbolic
analysis and the first numeric factorization) from `factorize!` (a numeric
refactorization reusing the analysis), which is exactly the split cuDSS wants.

That split is what makes this worth doing. The symbolic analysis depends only
on the circuit topology and the retained mode coupling, neither of which
changes across Newton steps, while the values change at every step. On a two
tone travelling wave amplifier the mode block diagonal preconditioner takes
about 72 ms to analyze and factorize from scratch and about 0.44 ms to
refactorize, so the analysis is roughly 99% of the cost and is paid once.

`batched = true` offers the mode block diagonal to cuDSS as a *uniform batch*
instead: one analysis for the pattern every mode block shares, then a batched
numeric refactorization. See [`BatchedBlockLayout`](@ref) for the structure and
[`blockdiagonallayout`](@ref) for the conditions, which
[`batchedfactorization`](@ref) checks; when they do not hold, the unbatched
path is used anyway.

!!! note "Why the batch is off by default"
    It was measured and it loses. cuDSS does not need to be told that a block
    diagonal matrix is a batch: given the whole matrix it discovers the
    independent blocks and works on them together, and given the blocks
    separately it works through them one at a time. On a two tone line with 30
    mode blocks of 4504 rows, one solve against the whole block diagonal is
    165 us against 5920 us for the batch, and a refactorization is 0.44 ms
    against 10.7 ms; the batched cost grows linearly in the number of blocks
    up to twelve and then jumps by a factor of forty and stays flat, which is
    a fallback inside the library rather than a load. The batch does win the
    symbolic analysis, 12 ms against 72 ms, because it analyzes one block
    rather than all thirty, but that is paid once and the per step losses
    swamp it: end to end the batch turned a 3.16 s solve into 5.63 s.

    Kept, off, because it is cheap to re-measure against a later cuDSS, and
    because the analysis result points somewhere useful: the saving is in the
    reordering, which the unbatched path could have too by analyzing one block
    and handing the replicated permutation to the whole matrix as
    `"user_perm"`.
"""
function CUDSSFactorization(; batched::Bool = false, kwargs...)
    return Factorization(
        (A; kws...) -> _cudss_factorize(A; kws...),
        (F, A; kws...) -> _cudss_factorize!(F, A; kws...),
        kwargs,
        # the batched form, for a caller which has established that its matrix
        # is a uniform batch. `nothing` disables the batched path entirely.
        batched ? (layout -> cudssbatchedfactorization(layout; kwargs...)) :
            nothing)
end

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
    nothing downstream can detect it. Measured on a 15372 row complex matrix:
    a batch of fifteen is correct to 5e-12 with twelve right hand sides and a
    batch of sixteen is wrong by a relative 1.4; with one or two right hand
    sides batches of 128 and 20 are correct.

    The cap costs nothing. Batching a frequency sweep saturates by about
    twelve systems anyway, so the limit sits above the useful range: on that
    same matrix a sweep took 0.0768 s at a batch of eight, 0.0723 s at twelve
    and 0.0730 s at fifteen.
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
        gpu::Factorization, cpu::Factorization = KLUfactorization())

Pick the factorization to use for `A`, preferring a batched GPU one when the
matrix really is a uniform batch of diagonal blocks.

The three-way choice is the whole point of the architecture:

1. `A` is block diagonal with uniform blocks, and a GPU factorization was
   supplied: use [`cudssbatchedfactorization`](@ref) on the shared pattern, so
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
    gpu::Union{Factorization,Nothing},
    cpu::Factorization = KLUfactorization())

    isnothing(gpu) && return cpu
    # a factorization with no batched form is used unchanged, and asking for
    # the layout would be a waste: it is a walk over every stored entry.
    isnothing(gpu.batched) && return gpu
    layout = blockdiagonallayout(A, blockof, nblocks)
    isnothing(layout) && return gpu
    return gpu.batched(layout)
end
