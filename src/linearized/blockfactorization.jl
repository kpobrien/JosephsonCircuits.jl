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

# B[dst[i], k] = nzval[src[i], k], B and nzval as (entries x batch)
@kernel function blockfillkernel!(B, @Const(nzval), @Const(src), @Const(dst))
    gid = @index(Global)
    @inbounds begin
        m = length(src)
        i = (gid - 1) % m + 1
        k = (gid - 1) ÷ m + 1
        B[dst[i], k] = nzval[src[i], k]
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
# B[dst[i], k] = nzval[src[i], k]*d[rows[i], k]*d[cols[i], k]
@kernel function blockfillscaledkernel!(B, @Const(nzval), @Const(src),
        @Const(dst), @Const(rows), @Const(cols), @Const(d))
    gid = @index(Global)
    @inbounds begin
        m = length(src)
        i = (gid - 1) % m + 1
        k = (gid - 1) ÷ m + 1
        B[dst[i], k] = nzval[src[i], k]*(d[rows[i], k]*d[cols[i], k])
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
pattern ([`blocknodegraph`](@ref)), ordered by KLU (`klunodeorder`),
amalgamated into supernodes ([`clustersymbolic`](@ref)), and every
stored entry of the matrix is mapped once to its place in a diagonal
block or a panel, so a refactorization is one scatter of the stored
values and one block LU.

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
    with the preconditioner's clusters; `blocksize`: the node block size.
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
    const blocksize::Int
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

"""
    factorize(f::BlockFactorization, A::SparseMatrixCSC; blocksize,
        backend = CPU(), nb = 1, target, refine = 6)

The [`SparseBlockFactorization`](@ref) of `A` with node blocks of
`blocksize` unknowns (the linearized solve passes its mode count), on
`backend`, holding `nb` systems of the pattern, filled from `A`'s values
in every slot. Precision from `f`; `Float32` factors of a double matrix
refine against the double residual for at most `refine` steps (six by
default, until the residual stops halving); `refine = 0` leaves the
solutions single precision solutions of the double system, computed
entirely in single precision. `target` is the amalgamation target in rows: on the host
`BLOCKTARGETROWS`, where LAPACK's LU of a supernode and the BLAS-3
panel products want large blocks; on a device the block size, no
amalgamation, since the batched LU and inverse of the diagonal blocks
(`getrf`/`getri` batched) are efficient only for small blocks and the
batch supplies the parallelism amalgamation gave a single system.
"""
function factorize(f::BlockFactorization, A::SparseMatrixCSC;
    blocksize::Integer, backend = CPU(), nb::Integer = 1,
    target::Integer = backend isa CPU ? BLOCKTARGETROWS : blocksize,
    refine::Integer = 6)
    n = size(A, 1)
    # the stored entries and their block positions are Int32 on the backend
    nnz(A) < typemax(Int32) || throw(ArgumentError(
        lazy"the block factorization indexes the $(nnz(A)) stored entries in Int32."))
    Tf = something(f.precision, real(eltype(A)))
    T = eltype(A) <: Complex ? Complex{Tf} : Tf
    TA = eltype(A)
    refinesteps = Tf === Float32 && real(TA) !== Float32 ? Int(refine) : 0
    refine = refinesteps > 0
    noderows, adj = blocknodegraph(A, blocksize)
    order = klunodeorder(adj)
    sym = clustersymbolic(noderows, adj, order; target)
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
    F = SparseBlockFactorization{T,A3,VI}(lu, Int(blocksize), fills,
        original, A, refine, refinesteps, scale, diagidx, backend, nothing)
    vals = tobackend(backend, repeat(nonzeros(A), 1, nb))
    return fillandfactorize!(F, vals)
end

# refactorize from a matrix with the pattern, in every slot of the batch
function fillandfactorize!(F::SparseBlockFactorization, A::SparseMatrixCSC)
    size(A) == size(F.A) && SparseArrays.getcolptr(A) == SparseArrays.getcolptr(F.A) &&
        rowvals(A) == rowvals(F.A) || throw(DimensionMismatch(
        "the matrix does not have the sparsity pattern the factorization was built for."))
    fillandfactorize!(F, tobackend(F.backend, repeat(nonzeros(A), 1, F.lu.nb)))
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
    backend = F.backend
    N = F.lu.N
    target(slot, o) = slot <= N ? (o ? F.original.D : F.lu.D)[slot] :
        slot <= 2N ? (o ? F.original.U : F.lu.U)[slot - N] :
        (o ? F.original.L : F.lu.L)[slot - 2N]
    # every block starts from zero: the fill positions the pattern does not
    # cover hold the last factorization's Schur updates
    for bs in (F.lu.D, F.lu.L, F.lu.U), Bk in bs
        fill!(Bk, zero(eltype(Bk)))
    end
    if F.refine
        # the originals, unscaled, in the matrix's precision
        for bs in (F.original.D, F.original.L, F.original.U), Bk in bs
            fill!(Bk, zero(eltype(Bk)))
        end
        kern = blockfillkernel!(backend, 256)
        for (slot, src, dst, _, _) in F.fills
            kern(reshape(target(slot, true), :, F.lu.nb), vals, src, dst;
                ndrange = length(src)*F.lu.nb)
        end
    end
    if isnothing(F.scale)
        kern = blockfillkernel!(backend, 256)
        for (slot, src, dst, _, _) in F.fills
            kern(reshape(target(slot, false), :, F.lu.nb), vals, src, dst;
                ndrange = length(src)*F.lu.nb)
        end
    else
        # single precision factors of the equilibrated matrix
        sk = blockscalekernel!(backend, 256)
        sk(F.scale, vals, F.diagidx; ndrange = F.lu.n*F.lu.nb)
        kern = blockfillscaledkernel!(backend, 256)
        for (slot, src, dst, rows, cols) in F.fills
            kern(reshape(target(slot, false), :, F.lu.nb), vals, src, dst, rows,
                cols, F.scale; ndrange = length(src)*F.lu.nb)
        end
    end
    KernelAbstractions.synchronize(backend)
    blocklu!(F.lu, backend)
    return F
end
refactorize!(::BlockFactorization, F::SparseBlockFactorization, A::SparseMatrixCSC) =
    fillandfactorize!(F, A)
refactorize!(::BlockFactorization, F::SparseBlockFactorization, vals::AbstractMatrix) =
    fillandfactorize!(F, vals)

# the work arrays for a right-hand side of `W` columns
function blockwork(F::SparseBlockFactorization{T}, W::Integer) where {T}
    w = F.work
    if isnothing(w) || w.W != W
        maxpanel = max(maximum(length.(F.lu.rowidx); init = 0), 1)
        alloc = (S, r, c) -> KernelAbstractions.zeros(F.backend, S, r, c, F.lu.nb)
        TA = eltype(F.A)
        w = (; W = Int(W), Z = alloc(T, F.lu.n, W), Y = alloc(T, F.lu.n, W),
            P = alloc(T, maxpanel, W),
            Zo = F.refine ? alloc(TA, F.lu.n, W) : nothing,
            Yo = F.refine ? alloc(TA, F.lu.n, W) : nothing,
            Po = F.refine ? alloc(TA, maxpanel, W) : nothing)
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
function blocksolve!(X::AbstractArray{<:Any,3}, F::SparseBlockFactorization{T},
    B::AbstractArray; transposed::Bool = false) where {T}
    backend = F.backend; nb = F.lu.nb
    W = size(B, 2)
    w = blockwork(F, W)
    Z, Y = w.Z, w.Y
    gather = blockgatherrowskernel!(backend, 256)
    scatter = blockscatterrowskernel!(backend, 256)
    scattersub = blockscattersubrowskernel!(backend, 256)
    if isnothing(F.scale)
        gather(Z, B, F.lu.perm; ndrange = F.lu.n*W*nb)
    else
        gathers = blockgatherscaledkernel!(backend, 256)
        gathers(Z, B, F.lu.perm, F.scale; ndrange = F.lu.n*W*nb)
    end
    substitute!(Y, F.lu, Z, w.P, backend; transposed)
    if isnothing(F.scale)
        scatter(X, Y, F.lu.perm; ndrange = F.lu.n*W*nb)
    else
        scatters = blockscatterscaledkernel!(backend, 256)
        scatters(X, Y, F.lu.perm, F.scale; ndrange = F.lu.n*W*nb)
    end
    KernelAbstractions.synchronize(backend)
    return X
end

"""
    blockresidual!(R, F::SparseBlockFactorization, X, B; transposed = false)

`R_k = B - A_k X_k` for the batch, or `B - transpose(A_k) X_k` with
`transposed`, from the matrix's own blocks kept for the refinement, a
batched dense product per block, in the matrix's precision.
"""
function blockresidual!(R::AbstractArray{<:Any,3}, F::SparseBlockFactorization,
    X::AbstractArray{<:Any,3}, B::AbstractMatrix; transposed::Bool = false)
    backend = F.backend; o = F.original; nb = F.lu.nb
    W = size(B, 2); TA = eltype(B)
    w = blockwork(F, W)
    Z, Y, Pn = w.Zo, w.Yo, w.Po
    gather = blockgatherrowskernel!(backend, 256)
    scatter = blockscatterrowskernel!(backend, 256)
    scattersub = blockscattersubrowskernel!(backend, 256)
    gather(Z, X, F.lu.perm; ndrange = F.lu.n*W*nb)
    gather(Y, B, F.lu.perm; ndrange = F.lu.n*W*nb)
    for P in 1:F.lu.N
        zP = view(Z, F.lu.range[P], :, :); yP = view(Y, F.lu.range[P], :, :)
        m = length(F.lu.rowidx[P])
        t = view(Pn, 1:m, :, :)
        batchedmul!(yP, o.D[P], zP, -one(TA), one(TA), transposed, false, backend)
        m == 0 && continue
        gather(t, Z, F.lu.rowidx[P]; ndrange = m*W*nb)
        if !transposed
            batchedmul!(yP, o.U[P], t, -one(TA), one(TA), false, false, backend)
            batchedmul!(t, o.L[P], zP, one(TA), zero(TA), false, false, backend)
        else
            batchedmul!(yP, o.L[P], t, -one(TA), one(TA), true, false, backend)
            batchedmul!(t, o.U[P], zP, one(TA), zero(TA), true, false, backend)
        end
        scattersub(Y, t, F.lu.rowidx[P]; ndrange = m*W*nb)
    end
    scatter(R, Y, F.lu.perm; ndrange = F.lu.n*W*nb)
    KernelAbstractions.synchronize(backend)
    return R
end

"""
    refinedsolve!(X, F::SparseBlockFactorization, B; transposed = false)

The batched solve, refined against the residual when the factors are in
single precision: `X += F \\ (B - A X)` while each step lowers the
residual, at most `refinesteps` steps. Exact factors solve once.
"""
function refinedsolve!(X::AbstractArray{<:Any,3}, F::SparseBlockFactorization,
    B::AbstractMatrix; transposed::Bool = false)
    blocksolve!(X, F, B; transposed)
    F.refine || return X
    R = similar(X); dX = similar(X)
    blockresidual!(R, F, X, B; transposed)
    # each system of the batch is judged on its own residual: one which
    # stagnates stops its own corrections and keeps its best iterate, and
    # the others go on. A batch-wide norm let one difficult system end the
    # refinement of the rest, or one large residual keep converged systems
    # iterating.
    nb = size(X, 3)
    rnorm = [norm(view(R, :, :, k)) for k in 1:nb]
    floor = 4*eps(real(eltype(B)))*norm(B)
    active = [rnorm[k] > floor for k in 1:nb]
    for _ in 1:F.refinesteps
        any(active) || break
        blocksolve!(dX, F, R; transposed)
        for k in 1:nb
            active[k] || continue
            view(X, :, :, k) .+= view(dX, :, :, k)
        end
        blockresidual!(R, F, X, B; transposed)
        for k in 1:nb
            active[k] || continue
            rnew = norm(view(R, :, :, k))
            if rnew < rnorm[k]/2
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
        budget = freememory(backend) ÷ 2)

The factorization of the linearized solve when none is given, by the
number of tones and the memory, the rule [`Automatic`](@ref) applies to
the nonlinear solve. One tone keeps the backend's sparse factorization
(KLU on the host, cuDSS on a device): its node blocks are small and the
sparse factorizations are as fast as or faster than the block one on
them (measured on the README's JTWPA examples: equal on the host, cuDSS
3.5x faster on a device). Two or more tones take
[`BlockFactorization`](@ref) in double when the factors of one system
([`blocksystembytes`](@ref)), times the host batches, fit in `budget`:
2 to 3.5x faster than KLU on the host and at parity with cuDSS on an
RTX 4090, where the double rate bounds both; otherwise the sparse
factorization.
"""
function linearizedfactorization(A::SparseMatrixCSC, Nmodes::Integer,
    ntones::Integer, backend; nbatches::Integer = 1,
    budget::Integer = freememory(backend) ÷ 2)
    sparsefactorization = backend isa CPU ? KLUfactorization() :
        CUDSSFactorization()
    ntones >= 2 || return sparsefactorization
    noderows, adj = blocknodegraph(A, Nmodes)
    sym = clustersymbolic(noderows, adj, klunodeorder(adj);
        target = backend isa CPU ? BLOCKTARGETROWS : Nmodes)
    bytes = blocksystembytes(Complex{Float64}, sym)
    systems = backend isa CPU ? max(nbatches, 1) : 1
    return bytes*systems <= budget ? BlockFactorization() : sparsefactorization
end
