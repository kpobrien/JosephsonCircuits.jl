
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
struct CUDSSFactorization <: AbstractFactorization
    kwargs::NamedTuple
end
CUDSSFactorization(; kwargs...) = CUDSSFactorization(NamedTuple(kwargs))
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

The largest uniform batch of systems to hand cuDSS, whatever the number of
right hand sides `nrhs`. Fifteen, for two independent reasons, one of
correctness and one of speed.

!!! warning "A wrong answer above fifteen systems with six or more right hand sides"
    cuDSS 0.7 and 0.8 (through CUDSS.jl 0.8.0) return silently wrong
    solutions from a uniform batch of sixteen or more systems once each has
    six or more right hand sides. Every system of the batch comes back
    wrong, by order one, while `cudss_get(solver, "info")` reports success
    and `"lu_nnz"` is unchanged, so nothing downstream can detect it. A
    batch of fifteen is correct with twelve right hand sides and a batch of
    sixteen is wrong by order one.

!!! warning "A step in the cost at sixteen systems, at every right hand side count"
    cuDSS 0.8 takes about eight times as long per refactorization and solve
    for a batch of sixteen as for a batch of fifteen, and then the same
    time for every batch from sixteen to sixty four. On a 600 by 600 sparse
    system the cost per refactorization and solve was 2.0 ms at fifteen
    systems and 16.6 ms at sixteen, and the ratio was 8.2, 8.4, 8.4 and 7.8
    at one, two, four and eight right hand sides. Because the cost above
    the step does not grow with the batch, splitting into chunks of fifteen
    always wins: sixty four systems as five chunks is about 9.5 ms against
    16.2 ms as one batch. This is not documented by NVIDIA and does not
    appear to have been reported.

The cap costs nothing on this path, since the speedup of batching a
frequency sweep through cuDSS saturates by about a dozen systems. It
applies only to the cuDSS batch: a [`SparseBlockFactorization`](@ref)
sweep sizes its batch by memory instead ([`blocksystembytes`](@ref)) and
profits from batches well past this. The `nrhs` argument is kept because
the first bound depends on it and the second does not, so a cuDSS release
which fixes one can be accommodated without changing the callers. Re-check
both against newer releases before raising it.
"""
uniformbatchlimit(nrhs::Integer) = 15

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
function _cudss_sweeprefactorize!(S)
    throw(ArgumentError("a batched solve on a device requires CUDSS.jl and CUDA.jl to be loaded."))
end
function _cudss_sweepapply!(S, X, B)
    throw(ArgumentError("a batched solve on a device requires CUDSS.jl and CUDA.jl to be loaded."))
end
