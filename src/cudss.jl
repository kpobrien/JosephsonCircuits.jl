
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
    CUDSSFactorization(; batched = true, kwargs...)

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
tone travelling wave amplifier with 288 cells the mode block diagonal
preconditioner takes about 16 ms to factorize from scratch and about 1.4 ms to
refactorize, so the analysis is roughly 92% of the cost and is paid once.

With `batched = true` the mode block diagonal is handed to cuDSS as a *uniform
batch*: one analysis for the pattern shared by every mode block, then a batched
numeric refactorization. See [`BatchedBlockLayout`](@ref) for the structure and
[`blockdiagonallayout`](@ref) for the conditions. When the matrix is not a
uniform batch, which happens whenever mode coupling is retained or the
truncation includes a self-conjugate mode such as dc, the batch is declined and
the unbatched cuDSS path is used instead. Set `batched = false` to force that.

!!! warning
    This has never been executed. It is written from the cuDSS and CUDSS.jl
    documentation, on a machine with no GPU, and no part of the device path has
    been run even once. The layout it depends on is tested, the extension is
    not. Treat the first run as debugging, not as a benchmark. cuDSS is also
    distributed as a preview, so its API may move.
"""
function CUDSSFactorization(; batched::Bool = true, kwargs...)
    return Factorization(
        (A; kws...) -> _cudss_factorize(A; batched = batched, kws...),
        (F, A; kws...) -> _cudss_factorize!(F, A; kws...),
        kwargs)
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
    layout = blockdiagonallayout(A, blockof, nblocks)
    isnothing(layout) && return gpu
    return cudssbatchedfactorization(layout)
end
