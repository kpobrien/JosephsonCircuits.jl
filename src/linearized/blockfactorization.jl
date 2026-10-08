# A block factorization of a sparse matrix, batched over the frequencies of a
# sweep: the linearized system's direct solve with dense node blocks (the
# `BlockFactorization` option of `hblinsolve`), sharing the symbolic structure
# of the preconditioner's clusters in solvers/blockclusters.jl.

# ---------------------------------------------------------------------------
# A block factorization of a sparse matrix, batched: the linearized
# system's direct solve with dense node blocks, the same symbolic structure
# as the preconditioner's clusters, filled from the stored values of the
# assembled matrix, and holding a batch of systems with one pattern (the
# frequencies of a device sweep) so that every dense operation is one
# batched call over the batch dimension

# B[dst[i], k] = nzval[src[i], k], `dst` the linear index within one of the
# batch's blocks B[:, :, k], nzval (entries x batch)
@kernel function blockfillkernel!(B, @Const(nzval), @Const(src), @Const(dst))
    gid = @index(Global)
    @inbounds begin
        m = length(src)
        i = (gid - 1) % m + 1
        k = (gid - 1) ÷ m + 1
        # `dst` indexes one slice, the batch's slices follow each other
        B[dst[i] + (k - 1)*(size(B, 1)*size(B, 2))] = nzval[src[i], k]
    end
end
# the equilibration of a single precision factorization: the matrix is
# scaled symmetrically by the inverse square roots of its diagonal
# magnitudes, `d[i, k] = 1/sqrt(|A_k[i, i]|)`, so that entries spanning
# many orders of magnitude (inverse inductances against capacitances
# times squared frequencies) become of order one before single precision
# arithmetic sees them; the solves scale the right-hand side and the
# solution back by the same vector
@kernel function blockscalekernel!(d, @Const(nzval), @Const(diagidx))
    gid = @index(Global)
    @inbounds begin
        n = length(diagidx)
        i = (gid - 1) % n + 1
        k = (gid - 1) ÷ n + 1
        p = Int(diagidx[i])
        a = p == 0 ? zero(real(eltype(nzval))) : abs(nzval[p, k])
        d[i, k] = a > 0 ? one(a)/sqrt(a) : one(a)
    end
end
# the same, times d[rows[i], k]*d[cols[i], k]
@kernel function blockfillscaledkernel!(B, @Const(nzval), @Const(src),
        @Const(dst), @Const(rows), @Const(cols), @Const(d))
    gid = @index(Global)
    @inbounds begin
        m = length(src)
        i = (gid - 1) % m + 1
        k = (gid - 1) ÷ m + 1
        B[dst[i] + (k - 1)*(size(B, 1)*size(B, 2))] =
            nzval[src[i], k]*(d[rows[i], k]*d[cols[i], k])
    end
end
# Z[i, j, k] = R[idx[i], j(, k)]*d[idx[i], k]
@kernel function blockgatherscaledkernel!(Z, @Const(R), @Const(idx), @Const(d))
    gid = @index(Global)
    @inbounds begin
        m = length(idx); W = size(Z, 2)
        i = (gid - 1) % m + 1
        j = ((gid - 1) ÷ m) % W + 1
        k = (gid - 1) ÷ (m*W) + 1
        r = ndims(R) == 2 ? R[idx[i], j] : R[idx[i], j, k]
        Z[i, j, k] = r*d[idx[i], k]
    end
end
# X[idx[i], j, k] = Z[i, j, k]*d[idx[i], k]
@kernel function blockscatterscaledkernel!(X, @Const(Z), @Const(idx), @Const(d))
    gid = @index(Global)
    @inbounds begin
        m = length(idx); W = size(Z, 2)
        i = (gid - 1) % m + 1
        j = ((gid - 1) ÷ m) % W + 1
        k = (gid - 1) ÷ (m*W) + 1
        X[idx[i], j, k] = Z[i, j, k]*d[idx[i], k]
    end
end
"""
    SparseBlockFactorization

The [`BlockFactorization`](@ref) of a sparse matrix whose unknowns come
in node blocks: the matrix of the linearized solve, `Nmodes` unknowns per
circuit node (and per auxiliary variable of the modified nodal analysis),
whose sparse LU has the block structure of the matrix itself. Built by
[`factorize`](@ref) from the pattern: the node graph is read off the
pattern and ordered by KLU, one node per supernode
([`blocksymbolic`](@ref)), and every stored entry of the matrix is mapped
once to its place in a diagonal block or a panel, so a refactorization is
one scatter of the stored values and one block LU.

The factorization holds `nb` systems with the one pattern, the
frequencies of a batch of a device sweep (one on the host): every block
is an array `(rows, columns, nb)` and every dense operation of the
factorization and the solves is one batched call over the batch
([`batchedinverse!`](@ref), [`batchedmul!`](@ref)), which is what fills a
device with the many small blocks of a chain. Solves take a right-hand
side shared by the batch or one per system, a GEMM per block, and the
transposed system from the same factors; factors in single precision
refine against the double residual formed from the matrix's own blocks
kept in double ([`blockresidual!`](@ref)), at close to twice the memory,
since the originals cost as much as the single precision factors.

# Fields
- `lu`: the [`BlockLU`](@ref), the factors and the Schur schedule shared
    with the preconditioner's clusters.
- `fills`: per block, the stored entries which land in it.
- `original`: the matrix's own blocks in its precision when refining.
- `scale`, `diagidx`: the equilibration of single precision factors, the
    symmetric diagonal scaling by the inverse square roots of the diagonal
    magnitudes, which brings a linearized matrix's entries (inverse
    inductances against capacitances times squared frequencies) to order
    one before single precision arithmetic sees them; `nothing` in double.
- `A`: the sparse pattern matrix; `refine` and `refinesteps`: whether the
    solves refine against the double residual and the most steps they
    take; `backend`; `work`: the work arrays of the last right-hand side
    width.
"""
mutable struct SparseBlockFactorization{T,A3,VI}
    const lu::BlockLU{T,A3,VI}
    # (slot, stored entry indices, linear indices in the block, rows,
    # columns): slot `P` is `D[P]`, `N + P` is `U[P]`, `2N + P` is `L[P]`
    const fills::Vector{Tuple{Int,VI,VI,VI,VI}}
    const original::Any
    A::SparseMatrixCSC
    const refine::Bool
    const refinesteps::Int
    # single precision factors are of the equilibrated matrix: `scale` is
    # `(n, nb)` in the factors' real precision, `diagidx` the stored
    # position of each diagonal entry (0 when not stored)
    const scale::Any
    const diagidx::Any
    const backend
    work::Any
end

"""
    blocknodegraph(A::SparseMatrixCSC, blocksize::Integer)

The node graph of a matrix whose unknowns come in contiguous blocks of
`blocksize` (a trailing shorter block allowed): the slot lists of the
nodes and the symmetric adjacency read off the pattern.
"""
function blocknodegraph(A::SparseMatrixCSC, blocksize::Integer)
    n = size(A, 1)
    size(A, 2) == n || throw(DimensionMismatch("the matrix must be square."))
    blocksize >= 1 || throw(ArgumentError("`blocksize` must be positive."))
    nnodes = cld(n, blocksize)
    noderows = [collect((a - 1)*blocksize + 1:min(a*blocksize, n)) for a in 1:nnodes]
    nodeof(s) = (s - 1) ÷ blocksize + 1
    adj = [Set{Int}() for _ in 1:nnodes]
    rows = rowvals(A)
    for j in 1:n
        b = nodeof(j)
        for k in nzrange(A, j)
            a = nodeof(rows[k])
            a == b && continue
            push!(adj[a], b); push!(adj[b], a)
        end
    end
    return noderows, [sort!(collect(x)) for x in adj]
end

"""
    blocksymbolic(A::SparseMatrixCSC, blocksize::Integer)

The symbolic structure ([`clustersymbolic`](@ref)) of the block
factorization of the linearized solve: the node graph of `A` with
`blocksize` unknowns per node ([`blocknodegraph`](@ref)), ordered by KLU,
one node per supernode. The factorization and the memory estimates that
choose it and size a device sweep's batches take it from here. Pivoting
stays within a supernode's diagonal block, and a merged block, which the
separators of a meshed circuit would form, can be ill conditioned at a
frequency where the matrix is not, which costs the solution more digits
than node by node elimination does; the preconditioner's clusters merge
([`amalgamate`](@ref)), since their accuracy costs only iterations.
"""
function blocksymbolic(A::SparseMatrixCSC, blocksize::Integer)
    noderows, adj = blocknodegraph(A, blocksize)
    return clustersymbolic(noderows, adj, klunodeorder(adj); maxrows = blocksize)
end

"""
    blocksystembytes(::Type{T}, sym; refine = false, TA = T)

The bytes one system of a [`SparseBlockFactorization`](@ref) in precision
`T` holds, from the symbolic structure: blocks, panels, inverses, scratch,
and the original blocks in `TA` when refining. What sizes the batch of a
device sweep.
"""
function blocksystembytes(::Type{T}, sym; refine::Bool = false,
    TA::Type = T) where {T}
    bytes = 0; scratch = Set{Tuple{Int,Int}}()
    for P in 1:sym.N
        n = length(sym.range[P]); m = length(sym.rowidxh[P])
        bytes += (2*n*n + 2*m*n)*sizeof(T)
        refine && (bytes += (n*n + 2*m*n)*sizeof(TA))
        push!(scratch, (n, n)); m > 0 && (push!(scratch, (m, n)); push!(scratch, (m, m)))
    end
    return bytes + sum(a*b for (a, b) in scratch; init = 0)*sizeof(T)
end

# the most steps the refinement of single precision block factors against
# the double residual takes by default: each step multiplies the error by
# about the accuracy of the factors, so a moderately conditioned system
# reaches double precision within a few, and a system stops refining at
# the first step that fails to reduce its residual enough
# (`refinedsolve!`)
const BLOCKREFINESTEPS = 6

"""
    blockanalysis(f::BlockFactorization, A::SparseMatrixCSC; blocksize,
        backend = CPU(), nb = 1, refine = BLOCKREFINESTEPS)

The [`SparseBlockFactorization`](@ref) of the pattern of `A` before its
first numeric factorization: the symbolic analysis with node blocks of
`blocksize` unknowns (the linearized solve passes its mode count), one
node per supernode ([`blocksymbolic`](@ref)), and the blocks of `nb`
systems of the pattern allocated on `backend` with the destinations of
the stored entries in them.
[`fillandfactorize!`](@ref) factorizes it from the values of a batch, as
the device sweep does for each of its batches, and [`factorize`](@ref)
from `A`'s. Precision from `f`; `Float32` factors of a double matrix
refine against the double residual for at most `refine` steps
(`BLOCKREFINESTEPS`, six, by default, until the residual stops halving,
see [`refinedsolve!`](@ref)); `refine = 0` leaves the
solutions single precision solutions of the double system, computed
entirely in single precision.
"""
function blockanalysis(f::BlockFactorization, A::SparseMatrixCSC;
    blocksize::Integer, backend = CPU(), nb::Integer = 1,
    refine::Integer = BLOCKREFINESTEPS)
    n = size(A, 1)
    # the stored entries and their block positions are Int32 on the backend
    nnz(A) < typemax(Int32) || throw(ArgumentError(
        lazy"the block factorization indexes the $(nnz(A)) stored entries in Int32."))
    Tf = something(f.precision, real(eltype(A)))
    T = eltype(A) <: Complex ? Complex{Tf} : Tf
    TA = eltype(A)
    refinesteps = Tf === Float32 && real(TA) !== Float32 ? Int(refine) : 0
    refine = refinesteps > 0
    sym = blocksymbolic(A, blocksize)
    (; N, perm, range, rowsnodes, rowidxh, paneloff, snode, nodepos, nrows) = sym
    dI = x -> tobackend(backend, Vector{Int32}(x))
    lu = blocklu(T, sym, backend; nb = nb)
    A3 = eltype(lu.D); VI = typeof(lu.perm)
    D, L, U = lu.D, lu.L, lu.U
    alloc = (S, r, c) -> KernelAbstractions.zeros(backend, S, r, c, nb)
    original = refine ? (D = [alloc(TA, size(D[P], 1), size(D[P], 2)) for P in 1:N],
        L = [alloc(TA, size(L[P], 1), size(L[P], 2)) for P in 1:N],
        U = [alloc(TA, size(U[P], 1), size(U[P], 2)) for P in 1:N]) : nothing
    # every stored entry's destination: its supernode pair decides the
    # diagonal block or panel, its position within them the linear index
    nodeof(s) = (s - 1) ÷ blocksize + 1
    pos = invperm(perm)
    srcs = [Int32[] for _ in 1:3N]; dsts = [Int32[] for _ in 1:3N]
    rows = rowvals(A)
    for j in 1:n
        pj = pos[j]; Q = snode[nodeof(j)]
        for k in nzrange(A, j)
            i = rows[k]; pi = pos[i]; P = snode[nodeof(i)]
            if P == Q
                slot = P
                r = pi - first(range[P]) + 1; c = pj - first(range[P]) + 1
                lin = r + (c - 1)*length(range[P])
            elseif P < Q
                slot = N + P                # U[P]: rows of P, panel columns
                r = pi - first(range[P]) + 1
                c = paneloff[P][nodeof(j)] + (pj - nodepos[nodeof(j)]) + 1
                lin = r + (c - 1)*length(range[P])
            else
                slot = 2N + Q               # L[Q]: panel rows, columns of Q
                r = paneloff[Q][nodeof(i)] + (pi - nodepos[nodeof(i)]) + 1
                c = pj - first(range[Q]) + 1
                lin = r + (c - 1)*size(L[Q], 1)
            end
            push!(srcs[slot], k); push!(dsts[slot], lin)
        end
    end
    scaled = Tf === Float32
    fills = Tuple{Int,VI,VI,VI,VI}[]
    colof = Int32[]
    for j in 1:n, _ in nzrange(A, j)
        push!(colof, j)
    end
    for slot in 1:3N
        isempty(srcs[slot]) && continue
        src = srcs[slot]
        push!(fills, (slot, dI(src), dI(dsts[slot]),
            scaled ? dI(rows[src]) : dI(Int32[]),
            scaled ? dI(colof[src]) : dI(Int32[])))
    end
    scale = scaled ? KernelAbstractions.zeros(backend, Tf, n, nb) : nothing
    diagidx = if scaled
        di = zeros(Int32, n)
        for j in 1:n, k in nzrange(A, j)
            rows[k] == j && (di[j] = k)
        end
        dI(di)
    else
        nothing
    end
    return SparseBlockFactorization{T,A3,VI}(lu, fills,
        original, A, refine, refinesteps, scale, diagidx, backend, nothing)
end

"""
    factorize(f::BlockFactorization, A::SparseMatrixCSC; blocksize,
        backend = CPU(), nb = 1, refine = BLOCKREFINESTEPS)

The [`SparseBlockFactorization`](@ref) of `A`, its analysis
([`blockanalysis`](@ref), whose keywords these are) factorized from `A`'s
values in every slot of its `nb` systems.
"""
factorize(f::BlockFactorization, A::SparseMatrixCSC; kwargs...) =
    fillandfactorize!(blockanalysis(f, A; kwargs...), A)

# refactorize from a matrix with the pattern, in every slot of the batch;
# one system on the host is filled from the matrix's stored values
# themselves, without a copy
function fillandfactorize!(F::SparseBlockFactorization, A::SparseMatrixCSC)
    size(A) == size(F.A) && SparseArrays.getcolptr(A) == SparseArrays.getcolptr(F.A) &&
        rowvals(A) == rowvals(F.A) || throw(DimensionMismatch(
        "the matrix does not have the sparsity pattern the factorization was built for."))
    vals = F.lu.nb == 1 && F.backend isa CPU ? reshape(nonzeros(A), :, 1) :
        tobackend(F.backend, repeat(nonzeros(A), 1, F.lu.nb))
    fillandfactorize!(F, vals)
    F.A = A
    return F
end

"""
    fillandfactorize!(F::SparseBlockFactorization, vals::AbstractMatrix)

Refactorize the batch from the stored values `vals`, one column per
system in the order of `nonzeros` of the pattern, on the backend: the
device sweep assembles the values of a batch of frequencies there.
"""
function fillandfactorize!(F::SparseBlockFactorization{T},
    vals::AbstractMatrix) where {T}
    size(vals) == (nnz(F.A), F.lu.nb) || throw(DimensionMismatch(
        lazy"the values must be $(nnz(F.A)) entries by $(F.lu.nb) systems."))
    # the fields held untyped are read behind this call
    fillblocks!(F.lu, F.fills, vals, F.original, F.scale, F.diagidx,
        F.backend)
    blocklu!(F.lu, F.backend)
    return F
end

# The blocks of the batch filled from the stored values `vals`: every
# block from zero, since the positions the pattern does not cover hold the
# last factorization's Schur updates; the originals, unscaled in the
# matrix's precision, when refining; and the blocks to factorize, of the
# equilibrated matrix when `scale` is given (single precision factors)
function fillblocks!(lu::BlockLU, fills, vals::AbstractMatrix, original,
    scale, diagidx, backend)
    N = lu.N
    # slot `P` is `D[P]`, `N + P` is `U[P]`, `2N + P` is `L[P]`
    block(bs, slot) = slot <= N ? bs.D[slot] :
        slot <= 2N ? bs.U[slot - N] : bs.L[slot - 2N]
    for bs in (lu.D, lu.L, lu.U), Bk in bs
        fill!(Bk, zero(eltype(Bk)))
    end
    if !isnothing(original)
        for bs in (original.D, original.L, original.U), Bk in bs
            fill!(Bk, zero(eltype(Bk)))
        end
        for (slot, src, dst, _, _) in fills
            fillblock!(block(original, slot), vals, src, dst, backend)
        end
    end
    if isnothing(scale)
        for (slot, src, dst, _, _) in fills
            fillblock!(block(lu, slot), vals, src, dst, backend)
        end
    else
        equilibrationscale!(scale, vals, diagidx, backend)
        for (slot, src, dst, rows, cols) in fills
            fillblock!(block(lu, slot), vals, src, dst, backend, rows, cols,
                scale)
        end
    end
    KernelAbstractions.synchronize(backend)
    return lu
end
refactorize!(::BlockFactorization, F::SparseBlockFactorization, A::SparseMatrixCSC) =
    fillandfactorize!(F, A)

# The fill of one block from the stored values of the batch, `B[dst[i], k]
# = vals[src[i], k]`, times the equilibration of its row and column when
# `d` is given; and the equilibration itself. The kernels above on a
# device, plain loops on the host, where a launch per block costs more
# than the entries it moves (`hostloop`).
function fillblock!(B::AbstractArray{<:Any,3}, vals::AbstractMatrix,
    src::AbstractVector, dst::AbstractVector, backend, rows = nothing,
    cols = nothing, d = nothing)
    nb = size(B, 3)
    if hostloop(backend, length(src)*nb)
        # `dst` indexes one slice, the batch's slices follow each other
        slice = size(B, 1)*size(B, 2)
        @inbounds for k in 1:nb, i in eachindex(src)
            v = vals[src[i], k]
            B[dst[i] + (k - 1)*slice] = isnothing(d) ? v :
                v*(d[rows[i], k]*d[cols[i], k])
        end
    # the block itself to the kernels, not a reshape of it, which on a
    # device would be a second array holding the block's memory until the
    # collector finds it (see `releasesweep!`)
    elseif isnothing(d)
        blockfillkernel!(backend, 256)(B, vals, src, dst;
            ndrange = length(src)*nb)
    else
        blockfillscaledkernel!(backend, 256)(B, vals, src, dst, rows, cols,
            d; ndrange = length(src)*nb)
    end
    return B
end
function equilibrationscale!(d::AbstractMatrix, vals::AbstractMatrix,
    diagidx::AbstractVector, backend)
    if hostloop(backend, length(d))
        @inbounds for k in axes(d, 2), i in axes(d, 1)
            p = Int(diagidx[i])
            a = p == 0 ? zero(real(eltype(vals))) : abs(vals[p, k])
            d[i, k] = a > 0 ? one(a)/sqrt(a) : one(a)
        end
    else
        blockscalekernel!(backend, 256)(d, vals, diagidx; ndrange = length(d))
    end
    return d
end

# the work arrays for a right-hand side of `W` columns: those of the
# substitutions, and when refining those of the residual and the residual
# itself, the correction and the residual norm and progress of each system
function blockwork(F::SparseBlockFactorization{T}, W::Integer) where {T}
    w = F.work
    if isnothing(w) || w.W != W
        maxpanel = max(maximum(length.(F.lu.rowidx); init = 0), 1)
        alloc = (S, r, c) -> KernelAbstractions.zeros(F.backend, S, r, c, F.lu.nb)
        TA = eltype(F.A)
        nb = F.lu.nb
        w = (; W = Int(W), Z = alloc(T, F.lu.n, W), Y = alloc(T, F.lu.n, W),
            P = alloc(T, maxpanel, W),
            Zo = F.refine ? alloc(TA, F.lu.n, W) : nothing,
            Yo = F.refine ? alloc(TA, F.lu.n, W) : nothing,
            Po = F.refine ? alloc(TA, maxpanel, W) : nothing,
            R = F.refine ? alloc(TA, F.lu.n, W) : nothing,
            dX = F.refine ? alloc(TA, F.lu.n, W) : nothing,
            rnorm = zeros(real(TA), nb), active = zeros(Bool, nb))
        F.work = w
    end
    return w
end

"""
    blocksolve!(X, F::SparseBlockFactorization, B; transposed = false)

Overwrite `X` (`n x W x nb`) with the solutions of `A_k X_k = B` (or
`transpose(A_k) X_k = B`) for the batch, `B` (`n x W`) shared by the
batch or (`n x W x nb`) one per system: every
operation a batched dense product over a whole block, the right-hand
side's columns and the batch. Forward substitution through the scaled
panels, back substitution through the panels and the inverses; for the
transposed system the same factors read the other way round, `(D ⊕ U)ᵀ`
first as a lower block triangular solve with the transposed inverses,
then `(I + L)ᵀ` backward.
"""
function blocksolve!(X::AbstractArray{<:Any,3}, F::SparseBlockFactorization,
    B::AbstractArray; transposed::Bool = false)
    # the fields held untyped are read behind this call
    return blocksolve!(X, F.lu, B, blockwork(F, size(B, 2)), F.scale,
        F.backend, transposed)
end
function blocksolve!(X::AbstractArray{<:Any,3}, lu::BlockLU, B::AbstractArray,
    w::NamedTuple, scale, backend, transposed::Bool)
    Z, Y = w.Z, w.Y
    if isnothing(scale)
        gatherrows!(Z, B, lu.perm, backend)
    else
        gatherscaled!(Z, B, lu.perm, scale, backend)
    end
    substitute!(Y, lu, Z, w.P, backend; transposed)
    if isnothing(scale)
        scatterrows!(X, Y, lu.perm, backend)
    else
        scatterscaled!(X, Y, lu.perm, scale, backend)
    end
    KernelAbstractions.synchronize(backend)
    return X
end

# the gather of an equilibrated right-hand side and the scatter of its
# solution back to the unscaled unknowns: the kernels above on a device,
# plain loops on the host (`hostloop`)
function gatherscaled!(Z::AbstractArray{<:Any,3}, R::AbstractArray,
    idx::AbstractVector, d::AbstractMatrix, backend)
    if hostloop(backend, length(Z))
        @inbounds for k in axes(Z, 3), j in axes(Z, 2), i in axes(Z, 1)
            r = ndims(R) == 2 ? R[idx[i], j] : R[idx[i], j, k]
            Z[i, j, k] = r*d[idx[i], k]
        end
    else
        blockgatherscaledkernel!(backend, 256)(Z, R, idx, d;
            ndrange = length(Z))
    end
    return Z
end
function scatterscaled!(X::AbstractArray{<:Any,3}, Z::AbstractArray{<:Any,3},
    idx::AbstractVector, d::AbstractMatrix, backend)
    if hostloop(backend, length(Z))
        @inbounds for k in axes(Z, 3), j in axes(Z, 2), i in axes(Z, 1)
            X[idx[i], j, k] = Z[i, j, k]*d[idx[i], k]
        end
    else
        blockscatterscaledkernel!(backend, 256)(X, Z, idx, d;
            ndrange = length(Z))
    end
    return X
end

"""
    blockresidual!(R, lu::BlockLU, o, X, B, w, backend, transposed::Bool)

`R_k = B - A_k X_k` for the batch, or `B - transpose(A_k) X_k` with
`transposed`, from the matrix's own blocks `o` kept for the refinement, a
batched dense product per block, in the matrix's precision, on the work
arrays `w` of the right-hand side's width.
"""
function blockresidual!(R::AbstractArray{<:Any,3}, lu::BlockLU, o::NamedTuple,
    X::AbstractArray{<:Any,3}, B::AbstractMatrix, w::NamedTuple, backend,
    transposed::Bool)
    TA = eltype(B)
    Z, Y, Pn = w.Zo, w.Yo, w.Po
    gatherrows!(Z, X, lu.perm, backend)
    gatherrows!(Y, B, lu.perm, backend)
    for P in 1:lu.N
        zP = view(Z, lu.range[P], :, :); yP = view(Y, lu.range[P], :, :)
        m = length(lu.rowidx[P])
        t = view(Pn, 1:m, :, :)
        batchedmul!(yP, o.D[P], zP, -one(TA), one(TA), transposed, backend)
        m == 0 && continue
        gatherrows!(t, Z, lu.rowidx[P], backend)
        if !transposed
            batchedmul!(yP, o.U[P], t, -one(TA), one(TA), false, backend)
            batchedmul!(t, o.L[P], zP, one(TA), zero(TA), false, backend)
        else
            batchedmul!(yP, o.L[P], t, -one(TA), one(TA), true, backend)
            batchedmul!(t, o.U[P], zP, one(TA), zero(TA), true, backend)
        end
        scattersubrows!(Y, t, lu.rowidx[P], backend)
    end
    scatterrows!(R, Y, lu.perm, backend)
    KernelAbstractions.synchronize(backend)
    return R
end

"""
    refinedsolve!(X, F::SparseBlockFactorization, B; transposed = false,
        rtol = 4*eps(real(eltype(B))), contraction = 1/2)

The batched solve, refined against the residual when the factors are in
single precision: `X += F \\ (B - A X)` for at most `refinesteps`
steps, each system until its residual is within `rtol*norm(B)`, the
rounding of the double residual, or until a step fails to bring its
residual below `contraction` times the last. Exact factors solve once.
"""
function refinedsolve!(X::AbstractArray{<:Any,3}, F::SparseBlockFactorization,
    B::AbstractMatrix; transposed::Bool = false,
    rtol::Real = 4*eps(real(eltype(B))), contraction::Real = 1/2)
    blocksolve!(X, F, B; transposed)
    F.refine || return X
    # the fields held untyped are read behind this call
    return refinedsolve!(X, F.lu, F.original, B, blockwork(F, size(B, 2)),
        F.scale, F.backend, F.refinesteps, transposed, rtol, contraction)
end
# the refinement of a solution `X` of the batch
function refinedsolve!(X::AbstractArray{<:Any,3}, lu::BlockLU, original,
    B::AbstractMatrix, w::NamedTuple, scale, backend, refinesteps::Integer,
    transposed::Bool, rtol::Real, contraction::Real)
    R, dX, rnorm, active = w.R, w.dX, w.rnorm, w.active
    blockresidual!(R, lu, original, X, B, w, backend, transposed)
    # each system of the batch is judged on its own residual: one which
    # stagnates stops its own corrections and keeps its best iterate, and
    # the others go on, so one difficult system neither ends the refinement
    # of the rest nor keeps converged systems iterating.
    nb = size(X, 3)
    floor = rtol*norm(B)
    for k in 1:nb
        rnorm[k] = norm(view(R, :, :, k))
        active[k] = rnorm[k] > floor
    end
    for _ in 1:refinesteps
        any(active) || break
        blocksolve!(dX, lu, R, w, scale, backend, transposed)
        for k in 1:nb
            active[k] || continue
            view(X, :, :, k) .+= view(dX, :, :, k)
        end
        blockresidual!(R, lu, original, X, B, w, backend, transposed)
        for k in 1:nb
            active[k] || continue
            rnew = norm(view(R, :, :, k))
            if rnew < contraction*rnorm[k]
                rnorm[k] = rnew
                rnew <= floor && (active[k] = false)
            else
                # no useful progress: keep the previous iterate when the
                # correction made it worse, and stop correcting this one
                rnew > rnorm[k] && (view(X, :, :, k) .-= view(dX, :, :, k))
                active[k] = false
            end
        end
    end
    return X
end

# the host form: one system with matrix right-hand sides, which is what
# the linearized sweep solves; the direct Newton solvers, which solve
# against vectors, refuse this factorization
function myldiv!(X::AbstractMatrix, F::SparseBlockFactorization, B::AbstractMatrix)
    refinedsolve!(reshape(X, size(X, 1), size(X, 2), 1), F, B)
    return X
end
# the transposed system against the same factors, what
# `trysolvetranspose!` solves with
struct TransposedSparseBlockFactorization{F}
    parent::F
end
LinearAlgebra.transpose(F::SparseBlockFactorization) =
    TransposedSparseBlockFactorization(F)
function myldiv!(X::AbstractMatrix, Ft::TransposedSparseBlockFactorization, B::AbstractMatrix)
    refinedsolve!(reshape(X, size(X, 1), size(X, 2), 1), Ft.parent, B; transposed = true)
    return X
end

"""
    linearizedfactorization(A::SparseMatrixCSC, Nmodes::Integer,
        ntones::Integer, backend; nbatches = 1,
        budget = memorybudget(backend))

The factorization of the linearized solve when none is given. One tone
takes the backend's sparse factorization, [`KLUfactorization`](@ref) on
the host and [`CUDSSFactorization`](@ref) on a device. Two or more tones
take [`BlockFactorization`](@ref) in double when the factors of one
system ([`blocksystembytes`](@ref)), times the `nbatches` host batches,
which hold one each, fit in `budget`, and the sparse factorization
otherwise. Several tones make the node blocks large and dense, which the
block factorization eliminates with dense LAPACK and BLAS-3 kernels where
a scalar factorization works one entry at a time. It pivots within a
node block only and stops at a singular one with a `SingularException`,
where [`hblinsolve`](@ref) solves its sweep again with the sparse
factorization, which pivots across the whole matrix; KLU also keeps the
sparsity within the node blocks, which the resonators of a line have,
and takes less memory.
"""
function linearizedfactorization(A::SparseMatrixCSC, Nmodes::Integer,
    ntones::Integer, backend; nbatches::Integer = 1,
    budget::Integer = memorybudget(backend))
    sparsefactorization = defaultfactorization(backend)
    ntones < 2 && return sparsefactorization
    bytes = blocksystembytes(Complex{Float64}, blocksymbolic(A, Nmodes))
    systems = backend isa CPU ? max(nbatches, 1) : 1
    return bytes*systems <= budget ? BlockFactorization() : sparsefactorization
end
