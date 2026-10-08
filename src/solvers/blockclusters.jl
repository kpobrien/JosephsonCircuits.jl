# The dense block factorization of a mode coupling preconditioner over the
# circuit graph: the `BlockFactorization` option, the symbolic elimination of
# the node graph, one cluster of modes per factorization, and the memory
# estimate `Automatic` chooses by.

"""
    BlockFactorization(singletons = nothing; precision = nothing, refine = true)

The [`AbstractFactorization`](@ref) of a [`ModeCouplingPreconditioner`](@ref)
by dense blocks over the circuit graph rather than by a scalar sparse solver,
and, handed a sparse matrix with a block size, of the linearized solve.

The harmonic balance Jacobian has two structures: its sparsity follows the
circuit graph, and every nonlinear connection carries a dense coupling
between harmonic modes. A scalar sparse solver sees only their product and
rediscovers the circuit's block structure, entry by entry, in its symbolic
analysis. This factorization eliminates the circuit graph instead, treating
every circuit node as a supernode whose block holds the real-layout slots of
all retained modes at that node: the nodes are ordered by KLU's analysis of
the circuit-node graph (block triangular form, then a fill reducing
ordering), merged along the elimination tree into larger supernodes only
where the merged block holds no entry the factors would not hold node by
node (the separators of a meshed circuit; a chain stays node by node, see
[`amalgamate`](@ref)), and eliminated with pivoted dense LU on the
diagonal blocks, dense products for the Schur updates and dense
matrix-vector products for the solves. The blocks are assembled on the
backend straight from the Fourier coefficients with
[`realstructureentry`](@ref), the same per entry value the sparse assembly
uses, so no sparse Jacobian is ever formed. The structure comes from the
circuit graph alone, whatever the circuit is.

The same specification serves two operators. Handed a
[`BlockJacobian`](@ref) it factorizes the preconditioner's mode clusters
over the circuit graph, as above; handed a `SparseMatrixCSC` with a
`blocksize` (see [`factorize`](@ref)) it builds a
[`SparseBlockFactorization`](@ref), the direct solve of the linearized
system in dense node blocks, batched over the frequencies of a device
sweep, with `precision` then meaning the precision of the factors of a
double matrix (equilibrated and refined when single). That solve merges
no nodes: a merged supernode of a meshed circuit can be ill conditioned
at a frequency where the matrix is not, and pivoting within it then
costs the solution digits ([`blocksymbolic`](@ref)), where in the
preconditioner it costs only iterations.

The linearized solve takes it at two or more tones, on the host and on a
device, when its factors fit the memory budget
([`linearizedfactorization`](@ref)): several tones make the node blocks
large and dense, and BLAS-3 kernels eliminate them, and the separators a
junction lattice fills, faster than a scalar factorization does. A scalar
factorization, [`KLUfactorization`](@ref) on the host, pivots across the
whole matrix where this one pivots within a supernode, and keeps the
sparsity within the node blocks which this one stores dense.

The coupling set of the preconditioner is honored at the level of its
*clusters*: the retained coupling graph of the modes is split into its
connected components, every component of two or more modes becomes one
block factorization over the circuit graph restricted to those modes'
slots, and the modes left single are solved by the mode block diagonal.
[`FullJacobian`](@ref) is therefore one factorization of the complete
Jacobian, an exact solve; a mask made of complete clusters (as
[`Clusters`](@ref) produces) is one factorization per cluster;
and a coupling set which is not a union of complete clusters, a [`HarmonicBand`](@ref)
say, is factorized on its closure, which keeps at least every coupling the
set asked for.

Factor storage grows as the square of the number of retained slots per
node and the arithmetic as its cube, which is what bounds it.

`singletons` is the sparse [`AbstractFactorization`](@ref) of the block diagonal
of the modes left single, the backend's default (KLU on the host, cuDSS on
a device) when `nothing`. `precision` is the floating point type of the
blocks, that of the preconditioner when `nothing`: `Float32` is the mixed
precision form, which halves the storage and runs at a device's single
precision rate while the iteration stays in double precision.

The block factorization does not pivot across its supernodes. A supernode
whose diagonal block is singular, which happens when the part of the
circuit eliminated up to it resonates at a mode frequency with the rest
held at zero, stops it with a `SingularException`, though the matrix is
not singular. The mode coupling preconditioner then falls back to the
backend's sparse factorization (see [`refactorize!`](@ref)); the
linearized solve does too when it chose the block factorization itself,
solving its sweep again with `KLUfactorization()` or
`CUDSSFactorization()`, and throws when it was given one.

`refine` concerns the linearized solve: single precision factors of a
double system refine their solutions against the double residual to double
accuracy by default; `refine = false` leaves them single precision
solutions computed entirely in single precision, equilibrated, with no
refinement, for the cases where single precision scattering parameters
suffice, at the speed of a device whose single precision rate far exceeds
its double one; the accuracy then falls with the conditioning of the
system, furthest on strongly resonant multi-tone lines.
The outputs are returned in double either way.
"""
struct BlockFactorization <: AbstractFactorization
    singletons::Union{Nothing,AbstractFactorization}
    precision::Union{Nothing,Type{<:AbstractFloat}}
    refine::Bool
end
function BlockFactorization(singletons::Union{Nothing,AbstractFactorization} = nothing;
    precision::Union{Nothing,Type{<:AbstractFloat}} = nothing,
    refine::Bool = true)
    singletons isa BlockFactorization && throw(ArgumentError(
        "the singleton modes are factorized by a sparse factorization, not a block one."))
    return BlockFactorization(singletons, precision, refine)
end
withprecision(f::BlockFactorization, ::Type{T}) where {T<:AbstractFloat} =
    BlockFactorization(f.singletons, T, f.refine)

# the sparse factorization of the singleton modes' block diagonal: the one
# given, or the backend's default
singletonfactorization(f::BlockFactorization, backend) =
    something(f.singletons, defaultfactorization(backend))

"""
    circuitnodegraph(pairptr, pairrow, invLnm::SparseMatrixCSC,
        Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, Nmodes::Integer,
        nnodes::Integer)

The adjacency lists of the circuit-node graph: two nodes are adjacent when a
junction pair table entry or a stored entry of the linear term matrices
couples them. This is the graph the block factorization eliminates.
"""
function circuitnodegraph(pairptr, pairrow, invLnm::SparseMatrixCSC,
    Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, Nmodes::Integer,
    nnodes::Integer)
    adj = [Set{Int}() for _ in 1:nnodes]
    for n2 in 1:nnodes, k in Int(pairptr[n2]):Int(pairptr[n2+1])-1
        n1 = Int(pairrow[k])
        n1 != n2 && (push!(adj[n1], n2); push!(adj[n2], n1))
    end
    for M in (invLnm, Gnm, Cnm)
        rv = rowvals(M)
        for j in 1:size(M, 2), k in nzrange(M, j)
            a = (rv[k] - 1) ÷ Nmodes + 1
            b = (j - 1) ÷ Nmodes + 1
            a != b && (push!(adj[a], b); push!(adj[b], a))
        end
    end
    return [sort!(collect(s)) for s in adj]
end

"""
    klunodeorder(adj::AbstractVector{<:AbstractVector{<:Integer}})

The elimination order of the circuit nodes from KLU's symbolic analysis of
the node graph: its block triangular form permutation followed by the fill
reducing ordering within the blocks. On a chain this walks the chain; on a
meshed circuit it fills less than a bandwidth reducing order and gives a
far shallower elimination tree.
"""
function klunodeorder(adj::AbstractVector{<:AbstractVector{<:Integer}})
    N = length(adj)
    N == 0 && return Int[]
    I_ = Int[]; J_ = Int[]
    for a in 1:N
        push!(I_, a); push!(J_, a)
        for b in adj[a]
            push!(I_, a); push!(J_, b)
        end
    end
    A = sparse(I_, J_, ones(length(I_)), N, N)
    K = KLU.KLUFactorization(A)
    KLU.klu_analyze!(K)
    # the column permutation, 0-based in KLU, read from the symbolic
    # analysis `K` owns, which its finalizer frees
    return GC.@preserve K unsafe_wrap(Array, K.symbolic.Q, N) .+ 1
end

"""
    eliminationtree(adj, order::AbstractVector{<:Integer})

The symbolic elimination of the node graph `adj` under `order`, returned as
`(parent, post, later)`. `parent[a]` is the node eliminated first among
those coupled to `a` after `a`, following the fill, or zero for a root;
`post` is the postorder of this elimination tree, the order the supernodes
are eliminated in; and `later[a]` lists the nodes after `a` that `a` is
coupled to once the fill is included. Every postorder of the tree
eliminates with the same fill, and every coupling of the filled graph joins
a node to one of its ancestors, so `later` holds in the postorder as well.
"""
function eliminationtree(adj, order::AbstractVector{<:Integer})
    N = length(adj)
    rank = invperm(order)
    nbr = [Set{Int}(adj[a]) for a in 1:N]
    parent = zeros(Int, N)
    later = [Int[] for _ in 1:N]
    for a in order
        lat = [b for b in nbr[a] if rank[b] > rank[a]]
        later[a] = lat
        isempty(lat) && continue
        parent[a] = lat[argmin(rank[lat])]
        for b in lat, c in lat
            b != c && push!(nbr[b], c)
        end
    end
    children = [Int[] for _ in 1:N]
    for a in order
        parent[a] == 0 || push!(children[parent[a]], a)
    end
    # children in elimination order, so that the postorder is a refinement
    # of `order`
    post = Int[]
    for r in order
        parent[r] == 0 || continue
        stack = [(r, 1)]
        while !isempty(stack)
            v, k = stack[end]
            if k <= length(children[v])
                stack[end] = (v, k + 1)
                push!(stack, (children[v][k], 1))
            else
                pop!(stack)
                push!(post, v)
            end
        end
    end
    return parent, post, later
end

"""
    amalgamate(parent, post, later, nrows; maxrows = typemax(Int))

Merge the elimination tree into its fundamental supernodes: walking the
postorder, a node joins the supernode of its child when it is that child's
parent, the child is its only child, and the child's structure after fill,
`later[child]`, is the node and the node's own structure `later[node]`.
The rest of a child's structure always lies within its parent's, so the
test is `length(later[child]) == length(later[node]) + 1`. The merged
diagonal block and panel then hold no entry the node by node factors would
not hold: a merge takes no zero into the factors, whatever the circuit, so
a chain keeps its block bidiagonal factors while the separators of a
meshed circuit, dense after fill, are eliminated as one block. `maxrows`
bounds the rows of a supernode; the linearized solve, which merges
nothing, passes the rows of one node ([`blocksymbolic`](@ref)). Each
supernode is returned as its list of nodes, in elimination order.
"""
function amalgamate(parent, post, later, nrows; maxrows::Integer = typemax(Int))
    N = length(parent)
    nchildren = zeros(Int, N)
    for a in 1:N
        parent[a] == 0 || (nchildren[parent[a]] += 1)
    end
    supernodes = Vector{Int}[]
    current = Int[]
    rows = 0
    for (k, a) in enumerate(post)
        push!(current, a)
        rows += nrows[a]
        nested = k < length(post) && parent[a] == post[k+1] &&
            nchildren[post[k+1]] == 1 &&
            length(later[a]) == length(later[post[k+1]]) + 1 &&
            rows + nrows[post[k+1]] <= maxrows
        if !nested
            push!(supernodes, current)
            current = Int[]
            rows = 0
        end
    end
    isempty(current) || push!(supernodes, current)
    return supernodes
end

# the dense block assembly: one work item per entry of one block whose rows
# and columns are given natural slots of the real layout
@kernel function blockassemblykernel!(B, @Const(rslots), @Const(cslots),
        Nmodes, Nfreq, @Const(ami), @Const(amc), @Const(pairptr),
        @Const(pairrow), @Const(pairjunc), @Const(paircoef), @Const(lmolj),
        @Const(rlinv), @Const(rlptr), @Const(phimatrix), @Const(lcolptr),
        @Const(lrowval), @Const(lnzval), @Const(gcolptr), @Const(growval),
        @Const(gnzval), @Const(wm), @Const(ccolptr), @Const(crowval),
        @Const(cnzval), @Const(wm2))
    gid = @index(Global)
    T = eltype(B)
    @inbounds begin
        m = size(B, 1)
        r = (gid - 1) % m + 1
        c = (gid - 1) ÷ m + 1
        rri = Int(rslots[r])
        rci = Int(cslots[c])
        acc = realstructureentry(T, rri, rci, Nmodes, Nfreq, ami, amc,
            pairptr, pairrow, pairjunc, paircoef, lmolj, rlinv, rlptr, phimatrix)
        ci = Int(rlinv[rri]); dr = rri - Int(rlptr[ci])
        cj = Int(rlinv[rci]); dc = rci - Int(rlptr[cj])
        lin = realblockterm(sparselookup(lcolptr, lrowval, lnzval, ci, cj),
            dr, dc)
        lin += realblockterm((im * wm[cj]) *
            sparselookup(gcolptr, growval, gnzval, ci, cj), dr, dc)
        lin += realblockterm((-1 * wm2[cj]) *
            sparselookup(ccolptr, crowval, cnzval, ci, cj), dr, dc)
        B[r, c] = acc + T(lin)
    end
end
# === the block LU shared by the preconditioner's clusters and the sweep ===
#
# One symbolic structure ([`clustersymbolic`](@ref)) serves two factorizations:
# the preconditioner's, one cluster of modes of the Jacobian assembled into
# the blocks, and the linearized sweep's, a batch of system matrices with
# one pattern filled from their stored values. Both are the same
# right-looking block LU over supernodes and the same substitutions, so the
# factors, the Schur schedule and those operations live here, batched: every
# block is an array `(rows, columns, nb)` and every dense operation is one
# batched call ([`batchedinverse!`](@ref), [`batchedmul!`](@ref)), a batch
# of one for a cluster. What differs, how the values get into the blocks and
# what wraps the solve, stays with each caller.

# Z[i, j, k] = R[idx[i], j, k]; R may be a matrix, shared by every k
@kernel function blockgatherrowskernel!(Z, @Const(R), @Const(idx))
    gid = @index(Global)
    @inbounds begin
        m = length(idx); W = size(Z, 2)
        i = (gid - 1) % m + 1
        j = ((gid - 1) ÷ m) % W + 1
        k = (gid - 1) ÷ (m*W) + 1
        Z[i, j, k] = ndims(R) == 2 ? R[idx[i], j] : R[idx[i], j, k]
    end
end
# X[idx[i], j, k] = Z[i, j, k]
@kernel function blockscatterrowskernel!(X, @Const(Z), @Const(idx))
    gid = @index(Global)
    @inbounds begin
        m = length(idx); W = size(Z, 2)
        i = (gid - 1) % m + 1
        j = ((gid - 1) ÷ m) % W + 1
        k = (gid - 1) ÷ (m*W) + 1
        X[idx[i], j, k] = Z[i, j, k]
    end
end
# X[idx[i], j, k] -= Z[i, j, k]
@kernel function blockscattersubrowskernel!(X, @Const(Z), @Const(idx))
    gid = @index(Global)
    @inbounds begin
        m = length(idx); W = size(Z, 2)
        i = (gid - 1) % m + 1
        j = ((gid - 1) ÷ m) % W + 1
        k = (gid - 1) ÷ (m*W) + 1
        X[idx[i], j, k] -= Z[i, j, k]
    end
end

# The gathers and scatters of the substitutions: the kernels above on a
# device, plain loops on the host, where a launch costs more than the rows
# it moves and allocates at every call (`hostloop`). A gather's source may
# be a matrix, the right-hand side every system of the batch shares.
function gatherrows!(Z::AbstractArray{<:Any,3}, R::AbstractArray,
    idx::AbstractVector, backend)
    if hostloop(backend, length(Z))
        if ndims(R) == 2
            @inbounds for k in axes(Z, 3), j in axes(Z, 2), i in axes(Z, 1)
                Z[i, j, k] = R[idx[i], j]
            end
        else
            @inbounds for k in axes(Z, 3), j in axes(Z, 2), i in axes(Z, 1)
                Z[i, j, k] = R[idx[i], j, k]
            end
        end
    else
        blockgatherrowskernel!(backend, 256)(Z, R, idx; ndrange = length(Z))
    end
    return Z
end
function scatterrows!(X::AbstractArray{<:Any,3}, Z::AbstractArray{<:Any,3},
    idx::AbstractVector, backend)
    if hostloop(backend, length(Z))
        @inbounds for k in axes(Z, 3), j in axes(Z, 2), i in axes(Z, 1)
            X[idx[i], j, k] = Z[i, j, k]
        end
    else
        blockscatterrowskernel!(backend, 256)(X, Z, idx; ndrange = length(Z))
    end
    return X
end
function scattersubrows!(X::AbstractArray{<:Any,3}, Z::AbstractArray{<:Any,3},
    idx::AbstractVector, backend)
    if hostloop(backend, length(Z))
        @inbounds for k in axes(Z, 3), j in axes(Z, 2), i in axes(Z, 1)
            X[idx[i], j, k] -= Z[i, j, k]
        end
    else
        blockscattersubrowskernel!(backend, 256)(X, Z, idx; ndrange = length(Z))
    end
    return X
end

# a batched identity: one work item per row of every slice
@kernel function blockidentitykernel!(A)
    gid = @index(Global)
    @inbounds begin
        n = size(A, 1)
        i = (gid - 1) % n + 1
        k = (gid - 1) ÷ n + 1
        A[i, i, k] = one(eltype(A))
    end
end

# `A[:, :, k] = I` for every slice, on the backend, without scalar indexing
function blockidentity!(A::AbstractArray{T,3}, backend) where {T}
    fill!(A, zero(T))
    blockidentitykernel!(backend, 256)(A; ndrange = size(A, 1)*size(A, 3))
    return A
end

# target[rowmap[i], colmap[j], k] -= W[i0 + i, j0 + j, k]
@kernel function blockschurkernel!(Tm, @Const(rowmap), @Const(colmap),
        @Const(W), i0, j0)
    gid = @index(Global)
    @inbounds begin
        m = length(rowmap); c = length(colmap)
        i = (gid - 1) % m + 1
        j = ((gid - 1) ÷ m) % c + 1
        k = (gid - 1) ÷ (m*c) + 1
        Tm[rowmap[i], colmap[j], k] -= W[i0 + i, j0 + j, k]
    end
end

"""
    batchedinverse!(Dinv, D, F, backend, pivots, status)

For each `k`, `Dinv[:, :, k]` becomes the inverse of `D[:, :, k]`, with
`F` scratch of the same size: a loop of pivoted dense LU solves on the
host, which pivot into `pivots`, at least as long as a block; one batched
call on a device (the CUDA extension), which keeps its pivots there. A
singular block throws a `SingularException` on the host, and on a device
sets the entry of its system in `status`, the device vector
[`blocklu!`](@ref) reads once after its last supernode (a batch of one
block, which cuSOLVER factorizes, throws at once).
"""
batchedinverse!(Dinv, D, F, backend, pivots, status) =
    batchedinverse!(Dinv, D, F, backend, status)
function batchedinverse!(Dinv::AbstractArray{T,3}, D::AbstractArray{T,3},
    F::AbstractArray{T,3}, ::CPU, pivots::Vector{BlasInt}, status) where {T}
    ipiv = view(pivots, 1:size(D, 1))
    for k in axes(D, 3)
        Fk = view(F, :, :, k)
        copyto!(Fk, view(D, :, :, k))
        Dk = view(Dinv, :, :, k)
        fill!(Dk, zero(T))
        for i in axes(Dk, 1); Dk[i, i] = one(T); end
        # LAPACK takes the pivot indices from Julia 1.11 on; before, and
        # for a type LAPACK does not have, the factorization allocates them
        @static if VERSION >= v"1.11"
            if T <: LinearAlgebra.BlasFloat
                info = last(LAPACK.getrf!(Fk, ipiv))
                info > 0 && throw(SingularException(info))
                LAPACK.getrs!('N', Fk, ipiv, Dk)
            else
                ldiv!(lu!(Fk), Dk)
            end
        else
            ldiv!(lu!(Fk), Dk)
        end
    end
    return Dinv
end

"""
    batchedmul!(C, A, B, alpha, beta, tA::Bool, backend)

`C[:, :, k] = alpha*op(A[:, :, k])*B[:, :, k] + beta*C[:, :, k]` for every
`k`, `op` the transpose when `tA` is set: a loop of `mul!` on the host, one
strided batched GEMM on a device (the CUDA extension).
"""
function batchedmul!(C::AbstractArray{T,3}, A::AbstractArray{T,3},
    B::AbstractArray{T,3}, alpha, beta, tA::Bool, ::CPU) where {T}
    for k in axes(C, 3)
        Ak = view(A, :, :, k); Bk = view(B, :, :, k); Ck = view(C, :, :, k)
        # a branch per transposition, so that each `mul!` is resolved
        # where it is compiled
        if tA
            mul!(Ck, transpose(Ak), Bk, alpha, beta)
        else
            mul!(Ck, Ak, Bk, alpha, beta)
        end
    end
    return C
end

# one Schur update: the product of a supernode's panels lands, through the
# row and column maps, at the offsets `i0`, `j0` of the diagonal block or a
# panel of a later supernode, in every slice of the batch
struct SchurTask{A3,VI}
    target::A3
    rowmap::VI
    colmap::VI
    i0::Int
    j0::Int
end

"""
    BlockLU

The factors of a batched block LU over a symbolic structure: the supernodes'
diagonal blocks, panels and explicit inverses, each `(rows, columns, nb)`,
the Schur updates between them and the work blocks. Allocated by
[`blocklu`](@ref) from a [`clustersymbolic`](@ref) structure, factorized in
place by [`blocklu!`](@ref) once its blocks hold values, and applied by
[`substitute!`](@ref).

# Fields
- `N`, `nb`, `n`: supernodes, batch, order.
- `perm`, `range`, `rowidx`: the elimination order (position in the
    factorization to natural slot) and each supernode's positions and panel
    rows.
- `D`, `L`, `U`, `Dinv`: the blocks; after `blocklu!`, `L` is scaled by the
    inverse and `Dinv` holds it.
- `tasks`, `scratch`: the Schur updates of each supernode and the work
    blocks by size.
- `pivots`: the pivot indices of the host's dense LU of a diagonal block.
"""
struct BlockLU{T,A3,VI}
    N::Int
    nb::Int
    n::Int
    perm::VI
    range::Vector{UnitRange{Int}}
    rowidx::Vector{VI}
    D::Vector{A3}
    L::Vector{A3}
    U::Vector{A3}
    Dinv::Vector{A3}
    tasks::Vector{Vector{SchurTask{A3,VI}}}
    scratch::Dict{Tuple{Int,Int},A3}
    pivots::Vector{BlasInt}
end

"""
    blocklu(::Type{T}, sym, backend; nb = 1)

Allocate the [`BlockLU`](@ref) of the symbolic structure `sym` in precision
`T` on `backend`, `nb` systems deep, with its Schur schedule: a supernode's
panel product lands, through index maps, in the diagonal block or a panel of
each later supernode it reaches.
"""
function blocklu(::Type{T}, sym, backend; nb::Integer = 1) where {T}
    (; N, perm, range, rowsnodes, rowidxh, paneloff, snode, nodepos, nrows) = sym
    # the positions and slots are Int32 on the backend
    length(perm) < typemax(Int32) || throw(ArgumentError(
        lazy"the block factorization indexes its $(length(perm)) unknowns in Int32."))
    dI = x -> tobackend(backend, Vector{Int32}(x))
    alloc = (r, c) -> KernelAbstractions.zeros(backend, T, r, c, nb)
    D = [alloc(length(range[P]), length(range[P])) for P in 1:N]
    L = [alloc(length(rowidxh[P]), length(range[P])) for P in 1:N]
    U = [alloc(length(range[P]), length(rowidxh[P])) for P in 1:N]
    Dinv = [alloc(length(range[P]), length(range[P])) for P in 1:N]
    A3 = eltype(D); VI = typeof(dI(Int32[]))
    loc(Q, us) = reduce(vcat, [collect(nodepos[u] - first(range[Q]) + 1 :
        nodepos[u] - first(range[Q]) + nrows[u]) for u in us]; init = Int[])
    pan(Q, us) = reduce(vcat, [collect(paneloff[Q][u] + 1 :
        paneloff[Q][u] + nrows[u]) for u in us]; init = Int[])
    tasks = Vector{Vector{SchurTask{A3,VI}}}(undef, N)
    for P in 1:N
        R = rowsnodes[P]
        tk = SchurTask{A3,VI}[]
        # the panel nodes grouped by supernode; they are contiguous in the panel
        groups = Vector{Tuple{Int,UnitRange{Int}}}()
        i = 1
        while i <= length(R)
            j = i
            while j < length(R) && snode[R[j+1]] == snode[R[i]]
                j += 1
            end
            push!(groups, (snode[R[i]], i:j))
            i = j + 1
        end
        for (Qa, ga) in groups, (Qb, gb) in groups
            i0 = paneloff[P][R[first(ga)]]
            j0 = paneloff[P][R[first(gb)]]
            if Qa == Qb
                tgt = D[Qa]; rowmap = loc(Qa, R[ga]); colmap = loc(Qb, R[gb])
            elseif Qa < Qb
                tgt = U[Qa]; rowmap = loc(Qa, R[ga]); colmap = pan(Qa, R[gb])
            else
                tgt = L[Qb]; rowmap = pan(Qb, R[ga]); colmap = loc(Qb, R[gb])
            end
            push!(tk, SchurTask{A3,VI}(tgt, dI(rowmap), dI(colmap), i0, j0))
        end
        tasks[P] = tk
    end
    scratch = Dict{Tuple{Int,Int},A3}()
    for P in 1:N
        nP = length(range[P]); m = length(rowidxh[P])
        haskey(scratch, (nP, nP)) || (scratch[(nP, nP)] = alloc(nP, nP))
        if m > 0
            haskey(scratch, (m, nP)) || (scratch[(m, nP)] = alloc(m, nP))
            haskey(scratch, (m, m)) || (scratch[(m, m)] = alloc(m, m))
        end
    end
    return BlockLU{T,A3,VI}(N, Int(nb), length(perm), dI(perm), range,
        [dI(r) for r in rowidxh], D, L, U, Dinv, tasks, scratch,
        Vector{BlasInt}(undef, maximum(length, range; init = 0)))
end

"""
    blocklu!(lu::BlockLU, backend)

The right-looking block LU of the blocks `lu` holds, in place, over the
batch: for each supernode in order, the pivoted dense LU of its diagonal
block and the explicit inverse from it, the panel below scaled by that
inverse, and the product of the scaled panel with the panel to the right
subtracted from the later blocks it reaches. After this the panels and
inverses are the factors. A singular diagonal block throws a
`SingularException`, on a device once the last supernode is done (see
[`batchedinverse!`](@ref)).
"""
function blocklu!(lu::BlockLU{T}, backend) where {T}
    # on a device the systems with a singular diagonal block, gathered
    # where the blocks are factorized and read once, after the last
    # supernode, rather than once per supernode
    status = backend isa CPU ? nothing :
        KernelAbstractions.zeros(backend, Cint, lu.nb)
    for P in 1:lu.N
        Dp = lu.D[P]; nP = size(Dp, 1)
        X = lu.Dinv[P]
        batchedinverse!(X, Dp, lu.scratch[(nP, nP)], backend, lu.pivots,
            status)
        m = size(lu.L[P], 1)
        m == 0 && continue
        Lp = lu.L[P]
        tmp = lu.scratch[(m, nP)]
        batchedmul!(tmp, Lp, X, one(T), zero(T), false, backend)
        copyto!(Lp, tmp)
        Wm = lu.scratch[(m, m)]
        batchedmul!(Wm, Lp, lu.U[P], one(T), zero(T), false, backend)
        for t in lu.tasks[P]
            schurupdate!(t, Wm, lu.nb, backend)
        end
    end
    KernelAbstractions.synchronize(backend)
    if !isnothing(status)
        k = findfirst(!=(0), Array(status))
        isnothing(k) || throw(SingularException(k))
    end
    return lu
end

# one Schur update, `t.target[rowmap[i], colmap[j], k] -= W[i0 + i, j0 + j,
# k]`: the kernel on a device, a plain loop on the host (`hostloop`)
function schurupdate!(t::SchurTask, W::AbstractArray{<:Any,3}, nb::Integer,
    backend)
    m = length(t.rowmap); c = length(t.colmap)
    if hostloop(backend, m*c*nb)
        Tm = t.target
        @inbounds for k in 1:nb, j in 1:c, i in 1:m
            Tm[t.rowmap[i], t.colmap[j], k] -= W[t.i0 + i, t.j0 + j, k]
        end
    else
        blockschurkernel!(backend, 256)(t.target, t.rowmap, t.colmap, W,
            t.i0, t.j0; ndrange = m*c*nb)
    end
    return t.target
end

"""
    substitute!(Y, lu::BlockLU, Z, Pw, backend; transposed = false)

The substitutions through the factors: overwrite `Y` (`n x W x nb`, in
factorization order) with the solution against `Z` (the right-hand side in
factorization order, overwritten), with `Pw` a panel work array of at least
`(maxpanel, W, nb)`. Forward substitution through the scaled panels, back
substitution through the panels and the inverses; for the transposed system
the same factors read the other way round, `(D ⊕ U)ᵀ` first as a lower
block triangular solve with the transposed inverses, then `(I + L)ᵀ`
backward. The caller gathers into `Z` and scatters out of `Y` by `perm`.
"""
function substitute!(Y::AbstractArray{<:Any,3}, lu::BlockLU{T},
    Z::AbstractArray{<:Any,3}, Pw::AbstractArray{<:Any,3}, backend;
    transposed::Bool = false) where {T}
    if !transposed
        for P in 1:lu.N
            m = length(lu.rowidx[P]); m == 0 && continue
            t = view(Pw, 1:m, :, :)
            batchedmul!(t, lu.L[P], view(Z, lu.range[P], :, :), one(T), zero(T),
                false, backend)
            scattersubrows!(Z, t, lu.rowidx[P], backend)
        end
        for P in lu.N:-1:1
            zP = view(Z, lu.range[P], :, :)
            m = length(lu.rowidx[P])
            if m > 0
                t = view(Pw, 1:m, :, :)
                gatherrows!(t, Y, lu.rowidx[P], backend)
                batchedmul!(zP, lu.U[P], t, -one(T), one(T), false, backend)
            end
            batchedmul!(view(Y, lu.range[P], :, :), lu.Dinv[P], zP, one(T),
                zero(T), false, backend)
        end
    else
        for P in 1:lu.N
            zP = view(Z, lu.range[P], :, :)
            yP = view(Y, lu.range[P], :, :)
            batchedmul!(yP, lu.Dinv[P], zP, one(T), zero(T), true, backend)
            m = length(lu.rowidx[P]); m == 0 && continue
            t = view(Pw, 1:m, :, :)
            batchedmul!(t, lu.U[P], yP, one(T), zero(T), true, backend)
            scattersubrows!(Z, t, lu.rowidx[P], backend)
        end
        for P in lu.N:-1:1
            yP = view(Y, lu.range[P], :, :)
            m = length(lu.rowidx[P]); m == 0 && continue
            t = view(Pw, 1:m, :, :)
            gatherrows!(t, Y, lu.rowidx[P], backend)
            batchedmul!(yP, lu.L[P], t, -one(T), one(T), true, backend)
        end
    end
    return Y
end

# === the preconditioner's clusters ===

"""
    ClusterBlocks

The block factorization of one cluster of modes over the circuit graph: the
[`BlockLU`](@ref) of the supernodes, the natural slots its blocks are
assembled from, and the work arrays of a solve. Built by
[`clusterblocks`](@ref), assembled by [`assembleblocks!`](@ref), factorized
by [`blocklu!`](@ref) and applied by [`clustersolve!`](@ref).
"""
struct ClusterBlocks{T,A3,VI}
    modes::Vector{Int}
    lu::BlockLU{T,A3,VI}
    colslots::Vector{VI}           # natural slots of a supernode's columns
    rowslots::Vector{VI}           # natural slots of its panel rows
    z::A3                          # the gathered right-hand side, (n, 1, 1)
    w::A3                          # the solution in factorization order
    tmp::A3                        # a panel's worth
end

"""
    clustersymbolic(modes, adj, order, Nmodes::Integer, layout::ModeLayout;
        tree = eliminationtree(adj, order))

The symbolic block structure of one cluster, on the host and without
allocating any factor storage: the fundamental supernodes
([`amalgamate`](@ref)) of the elimination tree of the circuit-node graph
`adj` under the node `order`, restricted to the real-layout slots of
`modes`, the positions of every supernode, its panel rows after fill, and
the offsets the Schur updates scatter through. [`clusterblocks`](@ref)
allocates from it and [`blockfactorbytes`](@ref) sizes it. The elimination
`tree` ([`eliminationtree`](@ref)) depends on the graph and the order
alone, so a caller structuring several clusters of one graph computes it
once.
"""
function clustersymbolic(modes, adj, order, Nmodes::Integer,
    layout::ModeLayout; tree = eliminationtree(adj, order))
    nnodes = length(adj)
    # the slots of this cluster's modes at each node, in the real layout
    noderows = [Int[] for _ in 1:nnodes]
    for a in 1:nnodes, k in modes
        c = (a - 1)*Nmodes + k
        append!(noderows[a], Int(layout.ptr[c]):Int(layout.ptr[c+1])-1)
    end
    return clustersymbolic(noderows, adj, order; tree)
end

# the core: `noderows[a]` are the slots node `a` contributes, in order
function clustersymbolic(noderows::Vector{Vector{Int}}, adj, order;
    maxrows::Integer = typemax(Int), tree = eliminationtree(adj, order))
    nnodes = length(adj)
    nrows = length.(noderows)
    parent, post, later = tree
    nodes = amalgamate(parent, post, later, nrows; maxrows)
    N = length(nodes)
    snode = zeros(Int, nnodes); nodepos = zeros(Int, nnodes)
    noderank = zeros(Int, nnodes)
    perm = Int[]; range = UnitRange{Int}[]
    for P in 1:N
        lo = length(perm) + 1
        for a in nodes[P]
            snode[a] = P
            nodepos[a] = length(perm) + 1
            noderank[a] = length(perm)
            append!(perm, noderows[a])
        end
        push!(range, lo:length(perm))
    end
    # the panel rows of a supernode: the later nodes its nodes couple to,
    # fill included, outside the supernode
    rowsnodes = Vector{Int}[]
    for P in 1:N
        R = Set{Int}()
        for a in nodes[P], b in later[a]
            snode[b] > P && push!(R, b)
        end
        push!(rowsnodes, sort!(collect(R); by = b -> noderank[b]))
    end
    paneloff = [Dict{Int,Int}() for _ in 1:N]
    rowidxh = Vector{Vector{Int}}(undef, N)
    for P in 1:N
        idx = Int[]
        for b in rowsnodes[P]
            paneloff[P][b] = length(idx)
            append!(idx, nodepos[b]:nodepos[b]+nrows[b]-1)
        end
        rowidxh[P] = idx
    end
    return (; nodes, N, perm, range, rowsnodes, rowidxh, paneloff, snode,
        nodepos, nrows)
end

"""
    blockfactorbytes(::Type{T}, keep::AbstractMatrix{Bool}, adj, order,
        Nmodes::Integer, layout::ModeLayout)

The bytes a [`BlockFactorization`](@ref) in precision `T` of the coupling
mask `keep` will hold on the backend: the diagonal blocks, their inverses,
the panels, the scratch of the largest blocks and the work vectors, over
every cluster of the mask. Exact for the
floating point storage the factorization allocates (the `Int32` index
maps are not counted), from the symbolic analysis alone, so a caller can
decide whether the factors fit before building anything.
"""
function blockfactorbytes(::Type{T}, keep::AbstractMatrix{Bool}, adj, order,
    Nmodes::Integer, layout::ModeLayout) where {T}
    bytes = 0
    clusters = modeclusters(keep)
    isempty(clusters) && return 0
    tree = eliminationtree(adj, order)
    for modes in clusters
        sym = clustersymbolic(modes, adj, order, Nmodes, layout; tree)
        scratch = Set{Tuple{Int,Int}}()
        for P in 1:sym.N
            n = length(sym.range[P]); m = length(sym.rowidxh[P])
            bytes += 2*n*n + 2*m*n          # D, Dinv, L, U
            push!(scratch, (n, n))
            m > 0 && (push!(scratch, (m, n)); push!(scratch, (m, m)))
        end
        bytes += sum(a*b for (a, b) in scratch; init = 0)
        bytes += 2*length(sym.perm) + max(maximum(length, sym.rowidxh; init = 0), 1)
    end
    return bytes*sizeof(T)
end

"""
    clusterblocks(::Type{T}, modes, adj, order, Nmodes::Integer,
        layout::ModeLayout, backend; tree = eliminationtree(adj, order))

The [`ClusterBlocks`](@ref) of one cluster, its storage allocated on
`backend` in precision `T`: the block LU, a batch of one, of the symbolic
block structure [`clustersymbolic`](@ref) finds for `modes` on the
circuit-node graph `adj` under the node `order`, with the natural slots its
blocks are assembled from. `tree` is the elimination of the graph, shared
by the clusters of one graph.
"""
function clusterblocks(::Type{T}, modes, adj, order, Nmodes::Integer,
    layout::ModeLayout, backend; tree = eliminationtree(adj, order)) where {T}
    sym = clustersymbolic(modes, adj, order, Nmodes, layout; tree)
    (; N, perm, range, rowidxh) = sym
    lu = blocklu(T, sym, backend)
    dI = x -> tobackend(backend, Vector{Int32}(x))
    nc = length(perm)
    maxpanel = maximum(length.(rowidxh); init = 0)
    dV = n -> KernelAbstractions.zeros(backend, T, n, 1, 1)
    A3 = eltype(lu.D); VI = typeof(lu.perm)
    return ClusterBlocks{T,A3,VI}(collect(modes), lu,
        [dI(perm[range[P]]) for P in 1:N], [dI(perm[rowidxh[P]]) for P in 1:N],
        dV(nc), dV(nc), dV(max(maxpanel, 1)))
end

"""
    BlockStructure

What a [`ModeCouplingPreconditioner`](@ref) holds in place of a sparse
structure when its factorization is a [`BlockFactorization`](@ref): the
block factorizations of the mode clusters, the mode block diagonal for the
modes left single (`nothing` when there are none), the assembly ingredients
on the backend, refreshed when the system is rebound
([`refreshvalues!`](@ref)), work vectors in the factorization's precision,
and the backend.
"""
mutable struct BlockStructure{T,VT}
    const clusters::Vector
    singletons                     # a ModeCouplingPreconditioner of the block diagonal, or nothing
    const ingredients::NamedTuple
    const backend
    const rT::VT
    const xT::VT
end

# the ingredient arrays of the assembly on the backend, refreshed by
# `refreshvalues!` when the system is rebound
function blockingredients(::Type{T}, sys, junctions::JunctionStructure{T},
    layout::ModeLayout, backend) where {T}
    dI = x -> tobackend(backend, Vector{Int32}(x))
    return (; junctions = junctions,
        rlinv = dI(collect(layout.inv)), rlptr = dI(collect(layout.ptr)),
        lcolptr = dI(SparseArrays.getcolptr(sys.invLnm)),
        lrowval = dI(rowvals(sys.invLnm)),
        lnzval = tobackend(backend, copy(nonzeros(sys.invLnm))),
        gcolptr = dI(SparseArrays.getcolptr(sys.Gnm)),
        growval = dI(rowvals(sys.Gnm)),
        gnzval = tobackend(backend, copy(nonzeros(sys.Gnm))),
        wm = tobackend(backend, copy(sys.wmodesm.diag)),
        ccolptr = dI(SparseArrays.getcolptr(sys.Cnm)),
        crowval = dI(rowvals(sys.Cnm)),
        cnzval = tobackend(backend, copy(nonzeros(sys.Cnm))),
        wm2 = tobackend(backend, copy(sys.wmodes2m.diag)))
end

"""
    blockstructure(::Type{T}, sys, Amatrixindices::Matrix,
        Amatrixconjindices::Matrix, keep::AbstractMatrix{Bool},
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, layout::ModeLayout, singletons)

The [`BlockStructure`](@ref) for the coupling mask `keep`: one
[`ClusterBlocks`](@ref) per connected component of two or more modes of
the retained coupling graph, over the circuit-node graph ordered by
[`klunodeorder`](@ref). `singletons` is the block diagonal preconditioner
for the remaining modes, or `nothing`.
"""
function blockstructure(::Type{T}, sys, Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, keep::AbstractMatrix{Bool},
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, layout::ModeLayout, singletons) where {T}
    backend = sys.nonlineartermplan.backend
    nnodes = layout.dim ÷ Nmodes
    junctions = junctionstructure(T, Amatrixindices, Amatrixconjindices,
        sys.Ljb, sys.Lscale, Rbnm, Nmodes, Nbranches, Nfreq, backend)
    adj = circuitnodegraph(junctions.hostpairptr, junctions.hostpairrow,
        sys.invLnm, sys.Gnm, sys.Cnm, Nmodes, nnodes)
    order = klunodeorder(adj)
    clusters = Any[]
    mc = modeclusters(keep)
    if !isempty(mc)
        tree = eliminationtree(adj, order)
        for modes in mc
            push!(clusters, clusterblocks(T, modes, adj, order, Nmodes,
                layout, backend; tree))
        end
    end
    ing = blockingredients(T, sys, junctions, layout, backend)
    n = layout.rdim
    rT = KernelAbstractions.zeros(backend, T, n)
    return BlockStructure{T,typeof(rT)}(clusters, singletons, ing, backend,
        rT, similar(rT))
end

"""
    modeclusters(keep::AbstractMatrix{Bool})

The connected components of two or more modes of the coupling graph whose
edges are the off-diagonal `true` entries of `keep`, read symmetrically;
each as its sorted mode list.
"""
function modeclusters(keep::AbstractMatrix{Bool})
    N = size(keep, 1)
    parent = collect(1:N)
    find(i) = (while parent[i] != i; parent[i] = parent[parent[i]]; i = parent[i]; end; i)
    for i in 1:N, j in i+1:N
        if keep[i, j] || keep[j, i]
            a = find(i); b = find(j)
            a != b && (parent[a] = b)
        end
    end
    groups = Dict{Int,Vector{Int}}()
    for i in 1:N
        push!(get!(groups, find(i), Int[]), i)
    end
    return sort!([sort!(g) for g in values(groups) if length(g) >= 2];
        by = first)
end

# the modes of `keep` which belong to no cluster
function singletonmodes(keep::AbstractMatrix{Bool})
    N = size(keep, 1)
    inc = falses(N)
    for g in modeclusters(keep), m in g
        inc[m] = true
    end
    return [m for m in 1:N if !inc[m]]
end

"""
    refreshvalues!(S::BlockStructure, sys::HBSystem)

The assembly ingredients of a block structure refreshed from a system
rebound to new component values: the linear term matrices, the junction
coefficients and the pair table. The structure and the ordering are kept.
"""
function refreshvalues!(S::BlockStructure{T}, sys::HBSystem) where {T}
    g = S.ingredients
    refreshvalues!(g.junctions, sys.Ljb, sys.Lscale)
    copyto!(g.lnzval, nonzeros(sys.invLnm))
    copyto!(g.gnzval, nonzeros(sys.Gnm))
    copyto!(g.cnzval, nonzeros(sys.Cnm))
    copyto!(g.wm, sys.wmodesm.diag)
    copyto!(g.wm2, sys.wmodes2m.diag)
    isnothing(S.singletons) || rebind!(S.singletons, sys)
    return S
end

"""
    assembleblocks!(C::ClusterBlocks, S::BlockStructure, phimatrix)

Assemble the diagonal blocks and panels of one cluster from the Fourier
coefficients `phimatrix` of `cos(phi(t))` and the linear terms.
"""
function assembleblocks!(C::ClusterBlocks{T}, S::BlockStructure,
    phimatrix) where {T}
    g = S.ingredients; js = g.junctions
    kern = blockassemblykernel!(S.backend, 256)
    launch(B, rs, cs) = kern(B, rs, cs, js.nmodes, js.nfreq, js.ami, js.amc,
        js.pairptr, js.pairrow, js.pairjunc, js.paircoef, js.lmolj, g.rlinv,
        g.rlptr, phimatrix, g.lcolptr, g.lrowval, g.lnzval, g.gcolptr,
        g.growval, g.gnzval, g.wm, g.ccolptr, g.crowval, g.cnzval, g.wm2;
        ndrange = length(B))
    lu = C.lu
    slice(B) = reshape(B, size(B, 1), size(B, 2))
    for P in 1:lu.N
        launch(slice(lu.D[P]), C.colslots[P], C.colslots[P])
        if length(lu.rowidx[P]) > 0
            launch(slice(lu.L[P]), C.rowslots[P], C.colslots[P])
            launch(slice(lu.U[P]), C.colslots[P], C.rowslots[P])
        end
    end
    KernelAbstractions.synchronize(S.backend)
    return C
end

"""
    clustersolve!(x::AbstractVector, C::ClusterBlocks, r::AbstractVector,
        backend)

Overwrite the cluster's slots of `x` with the solution of the cluster's
factorized operator against the cluster's slots of `r`; the other slots of
`x` are untouched. The gather by the elimination order, the substitutions
of [`substitute!`](@ref) and the scatter back, a batch of one.
"""
function clustersolve!(x::AbstractVector, C::ClusterBlocks{T},
    r::AbstractVector, backend) where {T}
    lu = C.lu
    n = lu.n
    if hostloop(backend, n)
        # the gather and the scatter as loops over the vectors themselves
        perm = lu.perm
        @inbounds for i in 1:n
            C.z[i, 1, 1] = r[perm[i]]
        end
        substitute!(C.w, lu, C.z, C.tmp, backend)
        @inbounds for i in 1:n
            x[perm[i]] = C.w[i, 1, 1]
        end
        return x
    end
    gather = blockgatherrowskernel!(backend, 256)
    scatter = blockscatterrowskernel!(backend, 256)
    gather(C.z, reshape(r, :, 1), lu.perm; ndrange = n)
    substitute!(C.w, lu, C.z, C.tmp, backend)
    scatter(reshape(x, :, 1, 1), C.w, lu.perm; ndrange = n)
    KernelAbstractions.synchronize(backend)
    return x
end

"""
    BlockJacobian

The operator a [`BlockFactorization`](@ref) factorizes: a
[`BlockStructure`](@ref) and the Fourier coefficients of `cos(phi(t))` at
the current point. Handed to [`tryfactorize!`](@ref) in place of a sparse
matrix.
"""
struct BlockJacobian{S,P}
    structure::S
    phimatrix::P
end

# assemble every cluster and factorize it, and refactorize the block
# diagonal of the singleton modes; the structure is the factorization.
# Every factorization is into the structure's own storage, so there is no
# separate refactorization, and one which fails is not repeated.
function factorize(::BlockFactorization, A::BlockJacobian)
    S = A.structure
    for C in S.clusters
        assembleblocks!(C, S, A.phimatrix)
        blocklu!(C.lu, S.backend)
    end
    isnothing(S.singletons) || refactorize!(S.singletons)
    return S
end

# the solve: the block diagonal on every slot, then each cluster's exact
# solve on its own slots; in the factorization's precision
function myldiv!(x::AbstractVector, F::BlockStructure, b::AbstractVector)
    if isnothing(F.singletons)
        fill!(x, zero(eltype(x)))
    else
        applypreconditioner!(x, F.singletons, b)
    end
    isempty(F.clusters) && return x
    F.rT .= b
    F.xT .= x
    for C in F.clusters
        clustersolve!(F.xT, C, F.rT, F.backend)
    end
    x .= F.xT
    return x
end

"""
    freememory(backend)

The free memory of `backend` in bytes: the host's for `CPU()`; on a CUDA
backend (defined by the CUDA extension) the device's, with what CUDA.jl's
memory pool holds without an array in it, which an allocation takes
first. What [`Automatic`](@ref), [`linearizedfactorization`](@ref) and the
device sweep's batch size their choices against.
"""
freememory(::CPU) = Int(Sys.free_memory())
function freememory(backend)
    throw(ArgumentError(
        "the free memory of this backend is unknown; load CUDA.jl for a CUDA device."))
end

"""
    releasearray!(x)

Return the memory of the device array `x` to its backend's memory pool at
once, rather than when the collector finds the array; `x` must not be used
again. A no-op on the host; the CUDA extension frees a CUDA array.
"""
releasearray!(x) = nothing

# the share of a backend's free memory the factors and batches a solve
# sizes for itself may take, the rest left to its other arrays and to
# whatever else the machine runs
const FREEMEMORYSHARE = 1//2

"""
    memorybudget(backend)

The bytes the factors and batches a solve sizes for itself may take on
`backend`: `FREEMEMORYSHARE`, half, of its free memory
([`freememory`](@ref)). What
[`Automatic`](@ref) ([`resolveautomatic`](@ref)), the growth of a
preconditioner's coupling set ([`escalatepreconditioner!`](@ref)),
[`linearizedfactorization`](@ref) and the batch of the device sweep
([`devicesolutions`](@ref)) are held to.
"""
memorybudget(backend) = floor(Int, FREEMEMORYSHARE*freememory(backend))

"""
    circuitorder(sys, Rbnm::SparseMatrixCSC, Nmodes::Integer,
        Nbranches::Integer, layout::ModeLayout)

The circuit-node graph of the system and KLU's elimination order of it,
the two symbolic ingredients a block factorization and its memory
prediction share.
"""
function circuitorder(sys, Rbnm::SparseMatrixCSC, Nmodes::Integer,
    Nbranches::Integer, layout::ModeLayout)
    nnodes = layout.dim ÷ Nmodes
    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
    pairptr, pairrow, _, _ = junctionpairtable(Int32, Float32, sys.Ljb,
        nodesandsigns, nnodes)
    adj = circuitnodegraph(pairptr, pairrow, sys.invLnm, sys.Gnm, sys.Cnm,
        Nmodes, nnodes)
    return adj, klunodeorder(adj)
end

"""
    sparsefactorbytes(P::SparseMatrixCSC, ::Type{T}, ordering)

The bytes a sparse LU of the pattern `P` in precision `T` would hold,
each entry with its index, from the symbolic analysis alone: the entries
of `L` and `U` under `ordering`, the fill reducing ordering a KLU
factorization of `P` takes ([`fillordering`](@ref)), which are twice the
fill that ordering predicts, `fill` entries each with the diagonal counted
in both, or under the natural order, which bounds them, when `ordering` is
`nothing`. What escalation to a larger coupling set is budgeted against on
any backend; a device factorization orders differently, but the fill of
the same pattern is of the same size.
"""
function sparsefactorbytes(P::SparseMatrixCSC, ::Type{T},
    ordering::Union{Nothing,FillOrdering}) where {T}
    # without an ordering, which only a pattern of one column or a failure
    # of both orderings leaves, the natural order's fill bounds it
    fillcount = isnothing(ordering) ?
        first(symbolicfill(_symmetricpattern(P), 1:size(P, 1))) : ordering.fill
    return 2*fillcount*(sizeof(T) + sizeof(Int))
end
