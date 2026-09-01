

"""
    BatchedBlockLayout

The description of a mode block diagonal matrix as a *uniform batch*: a set of
diagonal blocks which are all the same size and share one sparsity pattern,
differing only in their values.

This is the structure the mode block diagonal preconditioner has whenever
every mode contributes the same number of real slots. Each block is the
circuit's node graph at one mode frequency, so all blocks have identical
structure and differ only through the frequency on the diagonal and the shared
pump-averaged nonlinear inductance. A batched direct solver wants exactly this:
one symbolic analysis for the shared pattern, then a numeric refactorization
per block with no further analysis.

The layout is computed once, from the sparsity pattern alone, and is valid for
every subsequent operating point because the pattern does not change with the
Newton iterate.

# Fields
- `nblocks`, `blocksize`: the batch shape.
- `slots`: `blocksize` by `nblocks`, mapping a block-local row to its row of
    the full matrix. Moved to the backend by [`tobackend`](@ref), because the
    right hand side gather and the solution scatter run once per Krylov
    iteration on whatever device the vectors live on.
- `nzindex`: `blocknnz` by `nblocks`, mapping a block-local stored entry to its
    index in `nonzeros` of the full matrix, so gathering the batch is a permuted
    copy with no searching. Host resident, because it is read only alongside
    the host sparse matrix it indexes, once per numeric refactorization.
- `colptr`, `rowval`: the shared CSC pattern of one block.

Returns `nothing` from [`blockdiagonallayout`](@ref) when the blocks are not
uniform, which is not an error: a circuit with a dc or Nyquist mode has
self-conjugate modes that collapse to a single real slot, so those blocks are
genuinely a different size and the caller must fall back to an unbatched path.
"""
struct BatchedBlockLayout{TS<:AbstractMatrix{Int}}
    nblocks::Int
    blocksize::Int
    slots::TS
    nzindex::Matrix{Int}
    colptr::Vector{Int}
    rowval::Vector{Int}
end

"""
    blockdiagonallayout(A::SparseMatrixCSC, blockof::AbstractVector{<:Integer},
        nblocks::Integer)

Describe `A` as a uniform batch of diagonal blocks, where `blockof[i]` gives
the block each row and column `i` belongs to. Returns a
[`BatchedBlockLayout`](@ref), or `nothing` if `A` is not block diagonal with
respect to `blockof`, if the blocks differ in size, or if they do not share a
sparsity pattern.

Returning `nothing` rather than throwing is deliberate: the batched path is an
optimization, and every condition it requires is a property of the circuit
rather than a mistake by the caller.
"""
function blockdiagonallayout(A::SparseMatrixCSC,
    blockof::AbstractVector{<:Integer}, nblocks::Integer)

    n = size(A, 1)
    (size(A, 2) == n && length(blockof) == n) || return nothing
    nblocks >= 1 || return nothing

    # group the slots of each block, in increasing order of global index
    counts = zeros(Int, nblocks)
    for i in 1:n
        1 <= blockof[i] <= nblocks || return nothing
        counts[blockof[i]] += 1
    end
    blocksize = counts[1]
    all(==(blocksize), counts) || return nothing
    blocksize == 0 && return nothing

    slots = Matrix{Int}(undef, blocksize, nblocks)
    filled = zeros(Int, nblocks)
    local2global = Vector{Int}(undef, n)   # global slot -> block-local index
    for i in 1:n
        b = blockof[i]
        filled[b] += 1
        slots[filled[b], b] = i
        local2global[i] = filled[b]
    end

    # walk the stored entries once, binning them by block. an entry which
    # couples two different blocks means A is not block diagonal here.
    rows = rowvals(A)
    localrow = [Int[] for _ in 1:nblocks]
    localcol = [Int[] for _ in 1:nblocks]
    localnz = [Int[] for _ in 1:nblocks]
    for j in 1:n
        b = blockof[j]
        for k in nzrange(A, j)
            blockof[rows[k]] == b || return nothing
            push!(localrow[b], local2global[rows[k]])
            push!(localcol[b], local2global[j])
            push!(localnz[b], k)
        end
    end

    blocknnz = length(localnz[1])
    all(b -> length(localnz[b]) == blocknnz, 1:nblocks) || return nothing

    # the shared pattern, taken from the first block, and the CSC column
    # pointers implied by it. entries arrive in column major order already,
    # because the outer loop above walks columns of A in order.
    colptr = zeros(Int, blocksize + 1)
    colptr[1] = 1
    for c in localcol[1]
        colptr[c+1] += 1
    end
    cumsum!(colptr, colptr)
    rowval = localrow[1]

    # every block must present the same pattern in the same order
    for b in 2:nblocks
        localrow[b] == rowval || return nothing
        localcol[b] == localcol[1] || return nothing
    end

    nzindex = Matrix{Int}(undef, blocknnz, nblocks)
    for b in 1:nblocks, k in 1:blocknnz
        nzindex[k, b] = localnz[b][k]
    end

    return BatchedBlockLayout(Int(nblocks), blocksize, slots, nzindex,
        colptr, rowval)
end

"""
    gatherblocks!(values::AbstractMatrix, A::SparseMatrixCSC,
        layout::BatchedBlockLayout)

Copy the stored entries of `A` into `values`, one column per block, in the
shared pattern order of `layout`. A permuted copy: no searching, and the same
permutation for every Newton step.
"""
function gatherblocks!(values::AbstractMatrix, A::SparseMatrixCSC,
    layout::BatchedBlockLayout)
    nz = nonzeros(A)
    nzindex = layout.nzindex
    size(values) == size(nzindex) || throw(DimensionMismatch(
        lazy"`values` is $(size(values)) but the layout needs $(size(nzindex))."))
    @inbounds for b in axes(nzindex, 2), k in axes(nzindex, 1)
        values[k, b] = nz[nzindex[k, b]]
    end
    return values
end

# one work item per block slot. `slots` lists every row of the full matrix
# exactly once, so the gather and the scatter are both permutations: no slot
# is read or written twice and neither needs an atomic.
@kernel function gatherblockskernel!(rhs, @Const(v), @Const(slots))
    k = @index(Global)
    @inbounds rhs[k] = v[slots[k]]
end

@kernel function scatterblockskernel!(v, @Const(sol), @Const(slots))
    k = @index(Global)
    @inbounds v[slots[k]] = sol[k]
end

"""
    gathervalues!(dest::AbstractArray, src::AbstractVector,
        index::AbstractArray)

`dest[k] = src[index[k]]` for every `k`, as a KernelAbstractions kernel on the
backend of `dest`, which `src` and `index` must share.

This is the one operation the batched path needs on the hot path, in two
guises: gathering a Krylov vector into per block right hand sides, where
`index` is `layout.slots`, and gathering the assembled values of the full
matrix into the batch's value array, where it is
[`blockrowmajorindex`](@ref). Both are permutations or injections, so no
element is written twice and no atomic is needed.
"""
function gathervalues!(dest::AbstractArray, src::AbstractVector,
    index::AbstractArray)
    size(dest) == size(index) || throw(DimensionMismatch(
        lazy"`dest` is $(size(dest)) but the index is $(size(index))."))
    backend = KernelAbstractions.get_backend(dest)
    kernel! = gatherblockskernel!(backend)
    kernel!(dest, src, index; ndrange = length(index))
    KernelAbstractions.synchronize(backend)
    return dest
end

"""
    scattervalues!(dest::AbstractVector, src::AbstractVector,
        index::AbstractArray)

`dest[index[k]] = src[k]` for every `k`, as a KernelAbstractions kernel on
the backend of `src`. The inverse of [`gathervalues!`](@ref) when `index` is
a permutation; `index` must not repeat, or the writes race.
"""
function scattervalues!(dest::AbstractVector, src::AbstractVector,
        index::AbstractArray)
    length(src) == length(index) || throw(DimensionMismatch(
        lazy"`src` has length $(length(src)) but the index has $(length(index))."))
    backend = KernelAbstractions.get_backend(src)
    kernel! = scatterblockskernel!(backend)
    kernel!(dest, src, index; ndrange = length(index))
    KernelAbstractions.synchronize(backend)
    return dest
end

"""
    gatherblocks!(rhs::AbstractMatrix, v::AbstractVector,
        layout::BatchedBlockLayout)

Scatter the vector `v` into per block right hand sides, one column per block.

Runs as a KernelAbstractions kernel on the backend of `rhs`, so it works
wherever the vectors of the Krylov iteration live. `layout.slots` must be on
that backend too; see [`tobackend`](@ref). This is applied once per
preconditioner solve, which is once per Krylov iteration, so it is on the
hot path of a device solve rather than in its setup.
"""
function gatherblocks!(rhs::AbstractMatrix, v::AbstractVector,
    layout::BatchedBlockLayout)
    return gathervalues!(rhs, v, layout.slots)
end

"""
    scatterblocks!(v::AbstractVector, sol::AbstractMatrix,
        layout::BatchedBlockLayout)

The inverse of [`gatherblocks!`](@ref) for vectors: collect per block solutions
back into a full length vector. Runs on the backend of `sol`.
"""
function scatterblocks!(v::AbstractVector, sol::AbstractMatrix,
    layout::BatchedBlockLayout)
    slots = layout.slots
    size(sol) == size(slots) || throw(DimensionMismatch(
        lazy"`sol` is $(size(sol)) but the layout needs $(size(slots))."))
    backend = KernelAbstractions.get_backend(sol)
    kernel! = scatterblockskernel!(backend)
    kernel!(v, sol, slots; ndrange = length(slots))
    KernelAbstractions.synchronize(backend)
    return v
end

"""
    tobackend(backend, layout::BatchedBlockLayout)

Move the parts of `layout` which the per solve gather and scatter read to the
given KernelAbstractions backend, leaving the rest on the host. Only `slots`
moves: `nzindex`, `colptr` and `rowval` are read alongside the host sparse
matrix, in the setup and in the numeric refactorization, and never by a
kernel.
"""
function tobackend(backend::Backend, layout::BatchedBlockLayout)
    return BatchedBlockLayout(layout.nblocks, layout.blocksize,
        tobackend(backend, layout.slots), layout.nzindex, layout.colptr,
        layout.rowval)
end

"""
    blockpattern(layout::BatchedBlockLayout)

The sparsity pattern every block of the batch shares, as a sparse matrix whose
stored values are all zero. The batch is uniform precisely because there is
only one of these, so a batched solver analyzes this matrix once and reuses the
result for every block and every Newton step.
"""
function blockpattern(layout::BatchedBlockLayout)
    return SparseMatrixCSC(layout.blocksize, layout.blocksize,
        copy(layout.colptr), copy(layout.rowval),
        zeros(length(layout.rowval)))
end

"""
    blockrowmajorindex(layout::BatchedBlockLayout)

`blocknnz` by `nblocks`, mapping a stored entry of one block *in row major
order* to its index in `nonzeros` of the full matrix.

`layout.nzindex` is the same map for column major order, which is how a host
`SparseMatrixCSC` stores its values. A device sparse matrix is compressed by
rows, so a batch handed to a device solver wants this one instead. Composing
the two permutations once, here, is what lets the per Newton step gather be a
single indexed copy rather than a gather followed by a permutation.
"""
function blockrowmajorindex(layout::BatchedBlockLayout)
    return layout.nzindex[cscvaluepermutation(blockpattern(layout)), :]
end

"""
    blockmatrix(A::SparseMatrixCSC, layout::BatchedBlockLayout, b::Integer)

The `b`th diagonal block of `A` as its own sparse matrix, built from the
shared pattern. Used to hand one representative block to a symbolic analysis,
and by the tests to check the layout against ordinary submatrix extraction.
"""
function blockmatrix(A::SparseMatrixCSC, layout::BatchedBlockLayout,
    b::Integer)
    1 <= b <= layout.nblocks || throw(BoundsError(layout, b))
    nz = nonzeros(A)
    vals = [nz[layout.nzindex[k, b]] for k in axes(layout.nzindex, 1)]
    return SparseMatrixCSC(layout.blocksize, layout.blocksize,
        copy(layout.colptr), copy(layout.rowval), vals)
end
