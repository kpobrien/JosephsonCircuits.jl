
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
    CUDSSFactorization(; precision = nothing, kwargs...)

An [`AbstractFactorization`](@ref) backed by NVIDIA's cuDSS direct sparse
solver, for use on a GPU.

`precision` is the floating point type the factors are held and the
triangular solves run in, `nothing` for the iteration's own. It is a
setting of a [`NewtonKrylov`](@ref) preconditioner: the factors of a
preconditioner need only make the Krylov solve converge, so `Float32`
halves their memory and runs them at a device's single precision rate,
the preconditioner assembling its matrix in that precision and
converting the residual and the correction around each solve. Such a
preconditioner is not exact, and a Krylov solve it fails is escalated to
the same factors in the iteration's precision
([`escalatepreconditioner!`](@ref)). A direct solve factorizes in the
precision of its iteration and refuses a factorization asking for
another. cuDSS and [`BlockFactorization`](@ref) honour a precision; KLU
and UMFPACK factorize in double whatever they are handed.

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
together, which is faster than handing it the blocks as a uniform batch.
"""
struct CUDSSFactorization <: AbstractFactorization
    precision::Union{Nothing,Type{<:AbstractFloat}}
    kwargs::NamedTuple
end
CUDSSFactorization(; precision::Union{Nothing,Type{<:AbstractFloat}} = nothing,
    kwargs...) = CUDSSFactorization(precision, NamedTuple(kwargs))
factorizationprecision(f::CUDSSFactorization) = f.precision
withprecision(f::CUDSSFactorization, ::Type{T}) where {T<:AbstractFloat} =
    CUDSSFactorization(T, f.kwargs)
# a factorization asking for a precision is handed a matrix built in it,
# which only a preconditioner does; any other matrix is refused rather
# than factorized in its own precision as if nothing had been asked
function factorize(f::CUDSSFactorization, A)
    isnothing(f.precision) || real(eltype(A)) === f.precision ||
        throw(ArgumentError(
            lazy"a CUDSSFactorization with `precision` = $(f.precision) was handed a matrix of $(eltype(A)); the precision is a setting of a NewtonKrylov preconditioner, and a direct solve factorizes in the precision of its iteration."))
    return _cudss_factorize(A; f.kwargs...)
end
refactorize!(f::CUDSSFactorization, F, A) = _cudss_factorize!(F, A; f.kwargs...)
solverkwargs(f::CUDSSFactorization) = f.kwargs

# Overridden by the CUDSS extension. The error names the missing packages
# rather than failing with a MethodError somewhere inside the solver.
function _cudss_factorize(A; kwargs...)
    throw(ArgumentError(
        "CUDSSFactorization requires CUDSS.jl and CUDA.jl to be loaded. Run `using CUDA, CUDSS` before calling the solver, or use the default KLUfactorization()."))
end

# Defined by the CUDSS extension: the refactorization of a factorization
# `_cudss_factorize` made, which only the extension makes.
function _cudss_factorize! end

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
    cuDSS 0.8 takes several times as long per refactorization and solve
    for a batch of sixteen as for a batch of fifteen, at every number of
    right hand sides, and then the same time for every batch from sixteen
    to sixty four. Because the cost above the step does not grow with the
    batch, splitting into chunks of fifteen always wins.

The cap costs nothing on this path, since the speedup of batching a
frequency sweep through cuDSS saturates by about a dozen systems. It
applies only to the cuDSS batch, which a frequency sweep also holds to its
memory budget ([`cudssbatchlimit`](@ref)): a
[`SparseBlockFactorization`](@ref) sweep sizes its batch by memory alone
([`blocksystembytes`](@ref)) and profits from batches well past this. The
`nrhs` argument is kept because the first bound depends on it and the
second does not, so a cuDSS release which fixes one can be accommodated
without changing the callers. Re-check both against newer releases before
raising it.
"""
uniformbatchlimit(nrhs::Integer) = 15

"""
    cudssbatchlimit(persystem::Integer, nrhs::Integer, budget::Integer)

The systems of a uniform cuDSS batch: as many as `budget` bytes hold at
`persystem` bytes a system, at least one, and no more than
[`uniformbatchlimit`](@ref) allows.
"""
cudssbatchlimit(persystem::Integer, nrhs::Integer, budget::Integer) =
    clamp(budget ÷ max(persystem, 1), 1, uniformbatchlimit(nrhs))

"""
    cudsssystembytes(A::SparseMatrixCSC, T, nrhs::Integer, backend;
        kwargs...)

The device memory cuDSS takes for one system of the pattern of `A`, with
`nrhs` right hand sides of element type `T`, by its own estimate after an
analysis of the pattern: the peak of its factorization, which a uniform
batch takes once for each system. `kwargs` are the options of the
factorization ([`solverkwargs`](@ref)), on which the estimate depends.
"""
function cudsssystembytes(A::SparseMatrixCSC, ::Type{T}, nrhs::Integer,
        backend; kwargs...) where {T}
    # the compressed rows of `A`, as the sweep hands it to cuDSS, returned
    # to the pool with the estimate
    At = sparse(transpose(A))
    pattern = (tobackend(backend, SparseArrays.getcolptr(At)),
        tobackend(backend, rowvals(At)), tobackend(backend, ones(T, nnz(At))))
    bytes = _cudss_systembytes(pattern..., nrhs; kwargs...)
    foreach(releasearray!, pattern)
    return bytes
end

# Overridden by the CUDSS extension: cuDSS's estimate, after its analysis,
# of the device memory one system of a pattern takes.
function _cudss_systembytes(rowptr, colind, nzval, nrhs; kwargs...)
    throw(ArgumentError(
        "estimating cuDSS's memory requires CUDSS.jl and CUDA.jl to be loaded. Run `using CUDA, CUDSS` before calling the solver, or leave `backend` at its default."))
end

# Overridden by the CUDSS extension: a uniform batch of systems sharing one
# sparsity pattern, analyzed once and then refactorized and solved as a batch.
function _cudss_sweep(rowptr, colind, nzval, X, B; kwargs...)
    throw(ArgumentError(
        "solving a frequency sweep on a device requires CUDSS.jl and CUDA.jl to be loaded. Run `using CUDA, CUDSS` before calling the solver, or leave `backend` at its default."))
end

# Defined by the CUDSS extension, for a batch `_cudss_sweep` made, which
# only the extension makes: the refactorization and solve of the batch,
# and the two apart for a time stepper which solves many times per
# refactorization.
function _cudss_sweepsolve! end
function _cudss_sweeprefactorize! end
function _cudss_sweepapply! end

# Overridden by the CUDSS extension: the factorization data of a finished
# sweep, destroyed at once rather than when the collector finds it.
_cudss_release!(S) = nothing
