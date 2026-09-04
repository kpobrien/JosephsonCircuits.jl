"""
    BlockFactorization(singletons = nothing; precision = nothing)

The [`Factorization`](@ref) of a [`ModeCouplingPreconditioner`](@ref) by
dense blocks over the circuit graph rather than by a scalar sparse solver.

The harmonic balance Jacobian has two structures: its sparsity follows the
circuit graph, and every nonlinear connection carries a dense coupling
between harmonic modes. A scalar sparse solver sees only their product and
rediscovers the circuit's block structure, entry by entry, in its symbolic
analysis. This factorization eliminates the circuit graph instead, treating
every circuit node as a supernode whose block holds the real-layout slots of
all retained modes at that node: the nodes are ordered by KLU's analysis of
the circuit-node graph (block triangular form, then a fill reducing
ordering), amalgamated along the elimination tree into supernodes of a few
hundred rows, and eliminated with pivoted dense LU on the diagonal blocks,
dense products for the Schur updates and dense matrix-vector products for
the solves. The blocks are assembled on the backend straight from the
Fourier coefficients with [`realstructureentry`](@ref), the same per entry
value the sparse assembly uses, so no sparse Jacobian is ever formed. The
structure comes from the circuit graph alone, whatever the circuit is.

The coupling set of the preconditioner is honored at the level of its
*clusters*: the retained coupling graph of the modes is split into its
connected components, every component of two or more modes becomes one
block factorization over the circuit graph restricted to those modes'
slots, and the modes left single are solved by the mode block diagonal as
before. `couplingmodes = :all` is therefore one factorization of the
complete Jacobian, an exact solve; a mask made of complete clusters (as the
spectral rule of `:clusters` produces) is one factorization per cluster;
and a coupling set which is not a union of complete clusters, `:band => p`
say, is factorized on its closure, which keeps at least every coupling the
set asked for.

Measured on a three tone line of 128 junctions with 527 modes on a GPU, the
single precision factorization of the complete Jacobian takes 0.4 s and
2.8 GB against 2.6 s and 10.5 GB for cuDSS on the measured band, converges
the linear solves in three Arnoldi steps against 85, and solves the
nonlinear problem in 3.2 s against 38.5 s; clusters halve the memory again
at 18 steps. Factor storage grows as the square of the number of retained
slots per node and the arithmetic as its cube, which is what bounds it.

`singletons` is the sparse [`Factorization`](@ref) of the block diagonal
of the modes left single, the backend's default (KLU on the host, cuDSS on
a device) when `nothing`. `precision` is the floating point type of the
blocks, that of the preconditioner when `nothing`: `Float32` is the mixed
precision form, which halves the storage and is two to four times faster
on a GPU at three outer steps per solve while the iteration stays in
double precision.
"""
function BlockFactorization(singletons::Union{Nothing,Factorization} = nothing;
    precision::Union{Nothing,Type{<:AbstractFloat}} = nothing)
    return Factorization(blockfactorize, blockrefactorize!,
        (; singletons, precision))
end

isblockfactorization(f::Factorization) = f.factorize === blockfactorize

# the sparse factorization of the singleton modes' block diagonal: the one
# given, or the backend's default as `hbnlsolve` picks it
function singletonfactorization(f::Factorization, backend)
    isnothing(f.kwargs.singletons) || return f.kwargs.singletons
    return backend isa CPU ? KLUfactorization() : CUDSSFactorization()
end

# Supernodes of a few hundred rows keep the dense kernels busy without
# spending fill on merged zeros: on the measured circuits the factorization
# time is flat between 256 and 1024 rows and the solve is six times cheaper
# than with single node blocks. A kernel size, not a numerical parameter.
const BLOCKTARGETROWS = 512

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
meshed circuit it halves the fill of a bandwidth reducing order and gives an
elimination tree many times shallower.
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
    sym = K.symbolic
    # the column permutation, 0-based in KLU
    q = unsafe_wrap(Array, sym.Q, N) .+ 1
    return copy(q)
end

"""
    eliminationtree(adj, order::AbstractVector{<:Integer})

The elimination tree of the node graph `adj` under `order`: `parent[a]` is
the node eliminated first among those coupled to `a` after `a`, following the
fill, or zero for a root. Returned with the postorder of the tree, which is
the order the supernodes are eliminated in.
"""
function eliminationtree(adj, order::AbstractVector{<:Integer})
    N = length(adj)
    rank = invperm(order)
    nbr = [Set{Int}(adj[a]) for a in 1:N]
    parent = zeros(Int, N)
    for a in order
        lat = [b for b in nbr[a] if rank[b] > rank[a]]
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
    return parent, post
end

"""
    amalgamate(parent, post, nrows, target::Integer)

Merge chains of the elimination tree into supernodes: walking the postorder,
a node joins the supernode of its child when it is that child's parent and
the child its only child, until the supernode holds `target` rows. Each
supernode is returned as its list of nodes, in elimination order.
"""
function amalgamate(parent, post, nrows, target::Integer)
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
        chain = k < length(post) && parent[a] == post[k+1] &&
            nchildren[post[k+1]] == 1 && rows < target
        if !chain
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
            pairptr, pairrow, pairjunc, paircoef, lmolj, rlinv, rlptr, rlinv,
            rlptr, phimatrix)
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
@kernel function blockgatherkernel!(z, @Const(r), @Const(idx))
    i = @index(Global)
    @inbounds z[i] = r[idx[i]]
end
@kernel function blockscatterkernel!(x, @Const(z), @Const(idx))
    i = @index(Global)
    @inbounds x[idx[i]] = z[i]
end
@kernel function blockscattersubkernel!(x, @Const(z), @Const(idx))
    i = @index(Global)
    @inbounds x[idx[i]] -= z[i]
end
# target[rowmap[i], colmap[j]] -= W[i0 + i, j0 + j]
@kernel function blockschurkernel!(Tm, @Const(rowmap), @Const(colmap),
        @Const(W), i0, j0)
    gid = @index(Global)
    @inbounds begin
        m = length(rowmap)
        i = (gid - 1) % m + 1
        j = (gid - 1) ÷ m + 1
        Tm[rowmap[i], colmap[j]] -= W[i0 + i, j0 + j]
    end
end

# one Schur update: the product of a supernode's panels lands, through
# index maps, in the diagonal block or a panel of a later supernode
struct BlockSchurTask{MT,VI}
    target::MT
    rowmap::VI
    colmap::VI
    i0::Int
    j0::Int
end

"""
    ClusterBlocks

The block factorization of one cluster of modes over the circuit graph: the
supernodes, their diagonal blocks, panels and explicit inverses, and the
scattered Schur updates between them. Built by [`clusterblocks`](@ref) and
factorized in place by [`factorizeblocks!`](@ref).
"""
struct ClusterBlocks{T,VT,VI,MT}
    modes::Vector{Int}
    N::Int
    perm::VI                       # position in the cluster -> natural slot
    range::Vector{UnitRange{Int}}  # positions of each supernode
    rowidx::Vector{VI}             # positions of each supernode's panel rows
    colslots::Vector{VI}           # natural slots of a supernode's columns
    rowslots::Vector{VI}           # natural slots of its panel rows
    D::Vector{MT}
    L::Vector{MT}                  # panel rows x columns; after the factorization, scaled by the inverse
    U::Vector{MT}                  # columns x panel rows
    Dinv::Vector{MT}
    tasks::Vector{Vector{BlockSchurTask{MT,VI}}}
    scratch::Dict{Tuple{Int,Int},MT}
    eyes::Dict{Int,MT}
    z::VT
    w::VT
    tmp::VT
end

# the bytes of the factors
function factorbytes(C::ClusterBlocks{T}) where {T}
    return (sum(length, C.D) + sum(length, C.L) + sum(length, C.U) +
        sum(length, C.Dinv))*sizeof(T)
end

"""
    clusterblocks(::Type{T}, modes, adj, order, Nmodes::Integer,
        layout::ModeLayout, backend; target = BLOCKTARGETROWS)

The symbolic block structure of one cluster: the supernodes of the
amalgamated elimination tree of the circuit-node graph `adj` under the node
`order`, restricted to the real-layout slots of `modes`; the panels from a
symbolic elimination on the node graph; and the index maps of the Schur
updates. Storage is allocated on `backend` in precision `T`.
"""
function clusterblocks(::Type{T}, modes, adj, order, Nmodes::Integer,
    layout::ModeLayout, backend; target = BLOCKTARGETROWS) where {T}
    nnodes = length(adj)
    # the slots of this cluster's modes at each node, in the real layout
    noderows = [Int[] for _ in 1:nnodes]
    for a in 1:nnodes, k in modes
        c = (a - 1)*Nmodes + k
        append!(noderows[a], Int(layout.ptr[c]):Int(layout.ptr[c+1])-1)
    end
    nrows = length.(noderows)
    parent, post = eliminationtree(adj, order)
    nodes = amalgamate(parent, post, nrows, target)
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
    # symbolic elimination on the node graph in the final order gives each
    # supernode the later nodes it couples to, fill included
    nbr = [Set{Int}(adj[a]) for a in 1:nnodes]
    for P in 1:N, a in nodes[P]
        lat = [b for b in nbr[a] if noderank[b] > noderank[a]]
        for b in lat, c in lat
            b != c && push!(nbr[b], c)
        end
    end
    rowsnodes = Vector{Int}[]
    for P in 1:N
        R = Set{Int}()
        for a in nodes[P], b in nbr[a]
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
    dI = x -> tobackend(backend, Vector{Int32}(x))
    dM = (r, c) -> KernelAbstractions.zeros(backend, T, r, c)
    D = [dM(length(range[P]), length(range[P])) for P in 1:N]
    L = [dM(length(rowidxh[P]), length(range[P])) for P in 1:N]
    U = [dM(length(range[P]), length(rowidxh[P])) for P in 1:N]
    Dinv = [dM(length(range[P]), length(range[P])) for P in 1:N]
    MT = eltype(D); VI = typeof(dI(Int32[]))
    tasks = Vector{Vector{BlockSchurTask{MT,VI}}}(undef, N)
    loc(Q, us) = reduce(vcat, [collect(nodepos[u] - first(range[Q]) + 1 :
        nodepos[u] - first(range[Q]) + nrows[u]) for u in us]; init = Int[])
    pan(Q, us) = reduce(vcat, [collect(paneloff[Q][u] + 1 :
        paneloff[Q][u] + nrows[u]) for u in us]; init = Int[])
    for P in 1:N
        R = rowsnodes[P]
        tk = BlockSchurTask{MT,VI}[]
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
            push!(tk, BlockSchurTask{MT,VI}(tgt, dI(rowmap), dI(colmap), i0, j0))
        end
        tasks[P] = tk
    end
    scratch = Dict{Tuple{Int,Int},MT}()
    eyes = Dict{Int,MT}()
    for P in 1:N
        n = length(range[P]); m = length(rowidxh[P])
        haskey(scratch, (n, n)) || (scratch[(n, n)] = dM(n, n))
        haskey(eyes, n) || (eyes[n] = tobackend(backend, Matrix{T}(I, n, n)))
        if m > 0
            haskey(scratch, (m, n)) || (scratch[(m, n)] = dM(m, n))
            haskey(scratch, (m, m)) || (scratch[(m, m)] = dM(m, m))
        end
    end
    nc = length(perm)
    maxpanel = maximum(length.(rowidxh); init = 0)
    dV = n -> KernelAbstractions.zeros(backend, T, n)
    return ClusterBlocks{T,typeof(dV(0)),VI,MT}(collect(modes), N, dI(perm),
        range, [dI(r) for r in rowidxh], [dI(perm[range[P]]) for P in 1:N],
        [dI(perm[rowidxh[P]]) for P in 1:N], D, L, U, Dinv, tasks, scratch,
        eyes, dV(nc), dV(nc), dV(max(maxpanel, 1)))
end

"""
    BlockStructure

What a [`ModeCouplingPreconditioner`](@ref) holds in place of a sparse
structure when its factorization is a [`BlockFactorization`](@ref): the
block factorizations of the mode clusters, the mode block diagonal for the
modes left single (`nothing` when there are none), the assembly ingredients
on the backend, and work vectors in the factorization's precision.
"""
mutable struct BlockStructure{T,VT}
    const clusters::Vector
    singletons                     # a ModeCouplingPreconditioner of the block diagonal, or nothing
    const ingredients::NamedTuple
    const backend
    const rT::VT
    const xT::VT
    const Nmodes::Int
    const Rbnm::SparseMatrixCSC
    const Nbranches::Int
end

# the ingredient arrays of the assembly on the backend, refreshed by
# `refreshvalues!` when the system is rebound
function blockingredients(::Type{T}, sys, Amatrixindices, Amatrixconjindices,
    pairptr, pairrow, pairjunc, paircoef, lmolj, layout::ModeLayout,
    Nmodes::Integer, Nfreq::Integer, backend) where {T}
    dI = x -> tobackend(backend, Vector{Int32}(x))
    dT = x -> tobackend(backend, Vector{T}(x))
    return (; Nmodes = Int(Nmodes), Nfreq = Int(Nfreq),
        ami = tobackend(backend, Matrix{Int32}(Amatrixindices)),
        amc = tobackend(backend, Matrix{Int32}(Amatrixconjindices)),
        pairptr = dI(pairptr), pairrow = dI(pairrow), pairjunc = dI(pairjunc),
        paircoef = dT(paircoef), lmolj = dT(lmolj),
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
    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
    pairptr, pairrow, pairjunc, paircoef = junctionpairtable(Int32, T,
        sys.Ljb, nodesandsigns, nnodes)
    lmolj = T[T(sys.Lmean / sys.Ljb.nzval[i]) for i in eachindex(sys.Ljb.nzval)]
    adj = circuitnodegraph(pairptr, pairrow, sys.invLnm, sys.Gnm, sys.Cnm,
        Nmodes, nnodes)
    order = klunodeorder(adj)
    clusters = Any[]
    for modes in modeclusters(keep)
        push!(clusters, clusterblocks(T, modes, adj, order, Nmodes, layout,
            backend))
    end
    ing = blockingredients(T, sys, Amatrixindices, Amatrixconjindices,
        pairptr, pairrow, pairjunc, paircoef, lmolj, layout, Nmodes, Nfreq,
        backend)
    n = layout.rdim
    rT = KernelAbstractions.zeros(backend, T, n)
    return BlockStructure{T,typeof(rT)}(clusters, singletons, ing, backend,
        rT, similar(rT), Int(Nmodes), Rbnm, Int(Nbranches))
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
    nnodes = (length(g.rlptr) - 1) ÷ S.Nmodes
    nodesandsigns = branchnodesandsigns(S.Rbnm, S.Nmodes, S.Nbranches)
    pairptr, pairrow, pairjunc, paircoef = junctionpairtable(Int32, T,
        sys.Ljb, nodesandsigns, nnodes)
    length(paircoef) == length(g.paircoef) || throw(ArgumentError(
        "the junctions changed between points, which a refreshed block structure cannot follow; build a new one."))
    copyto!(g.paircoef, Vector{T}(paircoef))
    copyto!(g.lmolj, T[T(sys.Lmean / sys.Ljb.nzval[i]) for i in eachindex(sys.Ljb.nzval)])
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
    g = S.ingredients
    kern = blockassemblykernel!(S.backend, 256)
    launch(B, rs, cs) = kern(B, rs, cs, g.Nmodes, g.Nfreq, g.ami, g.amc,
        g.pairptr, g.pairrow, g.pairjunc, g.paircoef, g.lmolj, g.rlinv,
        g.rlptr, phimatrix, g.lcolptr, g.lrowval, g.lnzval, g.gcolptr,
        g.growval, g.gnzval, g.wm, g.ccolptr, g.crowval, g.cnzval, g.wm2;
        ndrange = length(B))
    for P in 1:C.N
        launch(C.D[P], C.colslots[P], C.colslots[P])
        if length(C.rowidx[P]) > 0
            launch(C.L[P], C.rowslots[P], C.colslots[P])
            launch(C.U[P], C.colslots[P], C.rowslots[P])
        end
    end
    KernelAbstractions.synchronize(S.backend)
    return C
end

"""
    factorizeblocks!(C::ClusterBlocks, backend)

The right-looking block LU of an assembled cluster: for each supernode in
order, the pivoted dense LU of its diagonal block and the explicit inverse
from it, the panel below scaled by that inverse, and the product of the
scaled panel with the panel to the right subtracted from the later blocks
it reaches. After this the panels and inverses are the factors.
"""
function factorizeblocks!(C::ClusterBlocks{T}, backend) where {T}
    scat = blockschurkernel!(backend, 256)
    for P in 1:C.N
        Dm = C.D[P]
        n = size(Dm, 1)
        F = C.scratch[(n, n)]
        copyto!(F, Dm)
        Fl = lu!(F)
        # the inverse, as the solve of the factors against the identity
        X = C.Dinv[P]
        copyto!(X, C.eyes[n])
        ldiv!(Fl, X)
        m = size(C.L[P], 1)
        m == 0 && continue
        Lp = C.L[P]
        tmp = C.scratch[(m, n)]
        mul!(tmp, Lp, X)
        copyto!(Lp, tmp)
        W = C.scratch[(m, m)]
        mul!(W, Lp, C.U[P])
        for t in C.tasks[P]
            scat(t.target, t.rowmap, t.colmap, W, t.i0, t.j0;
                ndrange = length(t.rowmap)*length(t.colmap))
        end
    end
    KernelAbstractions.synchronize(backend)
    return C
end

"""
    clustersolve!(x::AbstractVector, C::ClusterBlocks, r::AbstractVector,
        backend)

Overwrite the cluster's slots of `x` with the solution of the cluster's
factorized operator against the cluster's slots of `r`; the other slots of
`x` are untouched. Forward substitution through the scaled panels, then
back substitution through the panels and the inverses, all dense
matrix-vector products.
"""
function clustersolve!(x::AbstractVector, C::ClusterBlocks{T},
    r::AbstractVector, backend) where {T}
    gather = blockgatherkernel!(backend, 256)
    scatter = blockscatterkernel!(backend, 256)
    scattersub = blockscattersubkernel!(backend, 256)
    z = C.z; w = C.w
    gather(z, r, C.perm; ndrange = length(z))
    for P in 1:C.N
        m = length(C.rowidx[P])
        m == 0 && continue
        t = view(C.tmp, 1:m)
        mul!(t, C.L[P], view(z, C.range[P]))
        scattersub(z, t, C.rowidx[P]; ndrange = m)
    end
    for P in C.N:-1:1
        zP = view(z, C.range[P])
        m = length(C.rowidx[P])
        if m > 0
            t = view(C.tmp, 1:m)
            gather(t, w, C.rowidx[P]; ndrange = m)
            mul!(zP, C.U[P], t, -one(T), one(T))
        end
        mul!(view(w, C.range[P]), C.Dinv[P], zP)
    end
    scatter(x, w, C.perm; ndrange = length(w))
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
# diagonal of the singleton modes; the structure is the factorization
function blockfactorize(A::BlockJacobian; kwargs...)
    S = A.structure
    for C in S.clusters
        assembleblocks!(C, S, A.phimatrix)
        factorizeblocks!(C, S.backend)
    end
    isnothing(S.singletons) || refactorize!(S.singletons)
    return S
end
blockrefactorize!(F::BlockStructure, A::BlockJacobian; kwargs...) =
    blockfactorize(A)

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
