"""
    BlockFactorization(singletons = nothing; precision = nothing)

The [`AbstractFactorization`](@ref) of a [`ModeCouplingPreconditioner`](@ref)
by dense blocks over the circuit graph rather than by a scalar sparse solver.

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
before. [`FullJacobian`](@ref) is therefore one factorization of the complete
Jacobian, an exact solve; a mask made of complete clusters (as
[`Clusters`](@ref) produces) is one factorization per cluster;
and a coupling set which is not a union of complete clusters, a [`HarmonicBand`](@ref)
say, is factorized on its closure, which keeps at least every coupling the
set asked for.

Measured on a three tone line of 128 junctions with 527 modes on a GPU, the
single precision factorization of the complete Jacobian takes 0.4 s and
2.8 GB against 2.6 s and 10.5 GB for cuDSS on the measured band, converges
the linear solves in three Arnoldi steps against 85, and solves the
nonlinear problem in 3.2 s against 38.5 s; clusters halve the memory again
at 18 steps. Factor storage grows as the square of the number of retained
slots per node and the arithmetic as its cube, which is what bounds it.

`singletons` is the sparse [`AbstractFactorization`](@ref) of the block diagonal
of the modes left single, the backend's default (KLU on the host, cuDSS on
a device) when `nothing`. `precision` is the floating point type of the
blocks, that of the preconditioner when `nothing`: `Float32` is the mixed
precision form, which halves the storage and is two to four times faster
on a GPU at three outer steps per solve while the iteration stays in
double precision.
"""
struct BlockFactorization{S,P} <: AbstractFactorization
    singletons::S
    precision::P
end
function BlockFactorization(singletons::Union{Nothing,AbstractFactorization} = nothing;
    precision::Union{Nothing,Type{<:AbstractFloat}} = nothing)
    singletons isa BlockFactorization && throw(ArgumentError(
        "the singleton modes are factorized by a sparse factorization, not a block one."))
    return BlockFactorization(singletons, precision)
end

# the sparse factorization of the singleton modes' block diagonal: the one
# given, or the backend's default as `hbnlsolve` picks it
function singletonfactorization(f::BlockFactorization, backend)
    isnothing(f.singletons) || return f.singletons
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
    clustersymbolic(modes, adj, order, Nmodes::Integer, layout::ModeLayout;
        target = BLOCKTARGETROWS)

The symbolic block structure of one cluster, on the host and without
allocating any factor storage: the supernodes of the amalgamated
elimination tree of the circuit-node graph `adj` under the node `order`,
restricted to the real-layout slots of `modes`, the positions of every
supernode, its panel rows after fill, and the offsets the Schur updates
scatter through. [`clusterblocks`](@ref) allocates from it and
[`blockfactorbytes`](@ref) sizes it.
"""
function clustersymbolic(modes, adj, order, Nmodes::Integer,
    layout::ModeLayout; target = BLOCKTARGETROWS)
    nnodes = length(adj)
    # the slots of this cluster's modes at each node, in the real layout
    noderows = [Int[] for _ in 1:nnodes]
    for a in 1:nnodes, k in modes
        c = (a - 1)*Nmodes + k
        append!(noderows[a], Int(layout.ptr[c]):Int(layout.ptr[c+1])-1)
    end
    return clustersymbolic(noderows, adj, order; target)
end

# the core: `noderows[a]` are the slots node `a` contributes, in order
function clustersymbolic(noderows::Vector{Vector{Int}}, adj, order;
    target = BLOCKTARGETROWS)
    nnodes = length(adj)
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
    return (; nodes, N, perm, range, rowsnodes, rowidxh, paneloff, snode,
        nodepos, nrows)
end

"""
    blockfactorbytes(::Type{T}, keep::AbstractMatrix{Bool}, adj, order,
        Nmodes::Integer, layout::ModeLayout)

The bytes a [`BlockFactorization`](@ref) in precision `T` of the coupling
mask `keep` will hold on the backend: the diagonal blocks, their inverses,
the panels, the scratch of the largest blocks and the identity of each
block size, over every cluster of the mask. Exact for the storage the
factorization allocates, from the symbolic analysis alone, so a caller can
decide whether the factors fit before building anything.
"""
function blockfactorbytes(::Type{T}, keep::AbstractMatrix{Bool}, adj, order,
    Nmodes::Integer, layout::ModeLayout) where {T}
    bytes = 0
    for modes in modeclusters(keep)
        sym = clustersymbolic(modes, adj, order, Nmodes, layout)
        scratch = Set{Tuple{Int,Int}}(); eyes = Set{Int}()
        for P in 1:sym.N
            n = length(sym.range[P]); m = length(sym.rowidxh[P])
            bytes += 2*n*n + 2*m*n          # D, Dinv, L, U
            push!(scratch, (n, n)); push!(eyes, n)
            m > 0 && (push!(scratch, (m, n)); push!(scratch, (m, m)))
        end
        bytes += sum(a*b for (a, b) in scratch; init = 0) + sum(n*n for n in eyes; init = 0)
        bytes += 2*length(sym.perm) + max(maximum(length, sym.rowidxh; init = 0), 1)
    end
    return bytes*sizeof(T)
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
    sym = clustersymbolic(modes, adj, order, Nmodes, layout; target)
    return clusterblocks(T, modes, sym, backend)
end

# the core: allocate the blocks and the Schur tasks of a symbolic structure
function clusterblocks(::Type{T}, modes, sym::NamedTuple, backend) where {T}
    (; nodes, N, perm, range, rowsnodes, rowidxh, paneloff, snode, nodepos,
        nrows) = sym
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
function factorize(::BlockFactorization, A::BlockJacobian)
    S = A.structure
    for C in S.clusters
        assembleblocks!(C, S, A.phimatrix)
        factorizeblocks!(C, S.backend)
    end
    isnothing(S.singletons) || refactorize!(S.singletons)
    return S
end
refactorize!(f::BlockFactorization, F::BlockStructure, A::BlockJacobian) =
    factorize(f, A)

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

The free memory of `backend` in bytes: the host's for `CPU()`, the
device's on a CUDA backend (defined by the CUDA extension). What
[`Automatic`](@ref) sizes its choice of factorization against.
"""
freememory(::CPU) = Int(Sys.free_memory())
function freememory(backend)
    throw(ArgumentError(
        "the free memory of this backend is unknown; load CUDA.jl for a CUDA device."))
end

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
    sparsefactorbytes(P::SparseMatrixCSC, ::Type{T})

The bytes a sparse LU of the pattern `P` in precision `T` would hold,
from KLU's symbolic analysis of the pattern (its block triangular form
and fill-reducing order) and nothing numeric: the entries of L and U and
of the off-diagonal blocks, each with its index. What escalation to a
larger coupling set is budgeted against on any backend; a device
factorization orders differently, but the fill of the same pattern is of
the same size.
"""
function sparsefactorbytes(P::SparseMatrixCSC, ::Type{T}) where {T}
    A = SparseMatrixCSC(size(P, 1), size(P, 2), Vector{Int64}(P.colptr),
        Vector{Int64}(P.rowval), ones(Float64, nnz(P)))
    K = KLU.KLUFactorization(A)
    KLU.klu_analyze!(K)
    sym = K.symbolic
    entries = sym.lnz + sym.unz + sym.nzoff
    return round(Int, entries*(sizeof(T) + sizeof(Int)))
end

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
# target[rowmap[i], colmap[j], k] -= W[i0 + i, j0 + j, k]
@kernel function blockschurbatchkernel!(Tm, @Const(rowmap), @Const(colmap),
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
    batchedinverse!(Dinv, D, F, backend)

For each `k`, `Dinv[:, :, k]` becomes the inverse of `D[:, :, k]`, with
`F` scratch of the same size: a loop of pivoted dense LU solves on the
host, one batched call on a device (the CUDA extension).
"""
function batchedinverse!(Dinv::AbstractArray{T,3}, D::AbstractArray{T,3},
    F::AbstractArray{T,3}, ::CPU) where {T}
    for k in axes(D, 3)
        Fk = view(F, :, :, k)
        copyto!(Fk, view(D, :, :, k))
        LU = lu!(Fk)
        Dk = view(Dinv, :, :, k)
        fill!(Dk, zero(T))
        for i in axes(Dk, 1); Dk[i, i] = one(T); end
        ldiv!(LU, Dk)
    end
    return Dinv
end

"""
    batchedmul!(C, A, B, alpha, beta, tA::Bool, tB::Bool, backend)

`C[:, :, k] = alpha*op(A[:, :, k])*op(B[:, :, k]) + beta*C[:, :, k]` for
every `k`, `op` the transpose when the flag is set: a loop of `mul!` on
the host, one strided batched GEMM on a device (the CUDA extension).
"""
function batchedmul!(C::AbstractArray{T,3}, A::AbstractArray{T,3},
    B::AbstractArray{T,3}, alpha, beta, tA::Bool, tB::Bool, ::CPU) where {T}
    for k in axes(C, 3)
        Ak = view(A, :, :, k); Bk = view(B, :, :, k)
        mul!(view(C, :, :, k), tA ? transpose(Ak) : Ak, tB ? transpose(Bk) : Bk,
            alpha, beta)
    end
    return C
end

struct BatchSchurTask{A3,VI}
    target::A3
    rowmap::VI
    colmap::VI
    i0::Int
    j0::Int
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
device with the many small blocks of a chain. Solves take a matrix
right-hand side shared by the batch, a GEMM per block, and the transposed
system from the same factors; factors in single precision refine against
the double residual formed from the matrix's own blocks kept in double
([`blockresidual!`](@ref)), at one and a half times the memory.

# Fields
- `N`, `nb`, `n`, `blocksize`: supernodes, batch, order, node block size.
- `perm`, `range`, `rowidx`: the elimination order and each supernode's
    positions and panel rows.
- `D`, `L`, `U`, `Dinv`: the blocks, panels and inverses, `(., ., nb)`.
- `tasks`, `scratch`: the Schur updates and the work blocks.
- `fills`: per block, the stored entries which land in it.
- `original`: the matrix's own blocks in its precision when refining.
- `scale`, `diagidx`: the equilibration of single precision factors, the
    symmetric diagonal scaling by the inverse square roots of the diagonal
    magnitudes, which brings a linearized matrix's entries (inverse
    inductances against capacitances times squared frequencies) to order
    one before single precision arithmetic sees them; `nothing` in double.
- `A`: the sparse pattern matrix; `refine`; `backend`; `work`: the work
    arrays of the last right-hand side width.
"""
mutable struct SparseBlockFactorization{T,A3,VI}
    const N::Int
    const nb::Int
    const n::Int
    const blocksize::Int
    const perm::VI
    const range::Vector{UnitRange{Int}}
    const rowidx::Vector{VI}
    const D::Vector{A3}
    const L::Vector{A3}
    const U::Vector{A3}
    const Dinv::Vector{A3}
    const tasks::Vector{Vector{BatchSchurTask{A3,VI}}}
    const scratch::Dict{Tuple{Int,Int},A3}
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
        backend = CPU(), nb = 1, target, refine = true)

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
    alloc = (S, r, c) -> KernelAbstractions.zeros(backend, S, r, c, nb)
    D = [alloc(T, length(range[P]), length(range[P])) for P in 1:N]
    L = [alloc(T, length(rowidxh[P]), length(range[P])) for P in 1:N]
    U = [alloc(T, length(range[P]), length(rowidxh[P])) for P in 1:N]
    Dinv = [alloc(T, length(range[P]), length(range[P])) for P in 1:N]
    A3 = eltype(D); VI = typeof(dI(Int32[]))
    original = refine ? (D = [alloc(TA, size(D[P], 1), size(D[P], 2)) for P in 1:N],
        L = [alloc(TA, size(L[P], 1), size(L[P], 2)) for P in 1:N],
        U = [alloc(TA, size(U[P], 1), size(U[P], 2)) for P in 1:N]) : nothing
    # the Schur updates: a supernode's panel product lands, through index
    # maps, in the diagonal block or a panel of each later supernode it
    # reaches (the same construction as the preconditioner's clusters)
    loc(Q, us) = reduce(vcat, [collect(nodepos[u] - first(range[Q]) + 1 :
        nodepos[u] - first(range[Q]) + nrows[u]) for u in us]; init = Int[])
    pan(Q, us) = reduce(vcat, [collect(paneloff[Q][u] + 1 :
        paneloff[Q][u] + nrows[u]) for u in us]; init = Int[])
    tasks = Vector{Vector{BatchSchurTask{A3,VI}}}(undef, N)
    for P in 1:N
        R = rowsnodes[P]
        tk = BatchSchurTask{A3,VI}[]
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
            push!(tk, BatchSchurTask{A3,VI}(tgt, dI(rowmap), dI(colmap), i0, j0))
        end
        tasks[P] = tk
    end
    scratch = Dict{Tuple{Int,Int},A3}()
    for P in 1:N
        nP = length(range[P]); m = length(rowidxh[P])
        haskey(scratch, (nP, nP)) || (scratch[(nP, nP)] = alloc(T, nP, nP))
        if m > 0
            haskey(scratch, (m, nP)) || (scratch[(m, nP)] = alloc(T, m, nP))
            haskey(scratch, (m, m)) || (scratch[(m, m)] = alloc(T, m, m))
        end
    end
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
    F = SparseBlockFactorization{T,A3,VI}(N, Int(nb), n, Int(blocksize),
        dI(perm), range, [dI(r) for r in rowidxh], D, L, U, Dinv, tasks,
        scratch, fills, original, A, refine, refinesteps, scale, diagidx,
        backend, nothing)
    vals = tobackend(backend, repeat(nonzeros(A), 1, nb))
    return fillandfactorize!(F, vals)
end

# refactorize from a matrix with the pattern, in every slot of the batch
function fillandfactorize!(F::SparseBlockFactorization, A::SparseMatrixCSC)
    nnz(A) == nnz(F.A) && size(A) == size(F.A) || throw(DimensionMismatch(
        "the matrix does not have the pattern the factorization was built for."))
    fillandfactorize!(F, tobackend(F.backend, repeat(nonzeros(A), 1, F.nb)))
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
    size(vals) == (nnz(F.A), F.nb) || throw(DimensionMismatch(
        lazy"the values must be $(nnz(F.A)) entries by $(F.nb) systems."))
    backend = F.backend
    N = F.N
    target(slot, o) = slot <= N ? (o ? F.original.D : F.D)[slot] :
        slot <= 2N ? (o ? F.original.U : F.U)[slot - N] :
        (o ? F.original.L : F.L)[slot - 2N]
    # every block starts from zero: the fill positions the pattern does not
    # cover hold the last factorization's Schur updates
    for bs in (F.D, F.L, F.U), Bk in bs
        fill!(Bk, zero(eltype(Bk)))
    end
    if F.refine
        # the originals, unscaled, in the matrix's precision
        for bs in (F.original.D, F.original.L, F.original.U), Bk in bs
            fill!(Bk, zero(eltype(Bk)))
        end
        kern = blockfillkernel!(backend, 256)
        for (slot, src, dst, _, _) in F.fills
            kern(reshape(target(slot, true), :, F.nb), vals, src, dst;
                ndrange = length(src)*F.nb)
        end
    end
    if isnothing(F.scale)
        kern = blockfillkernel!(backend, 256)
        for (slot, src, dst, _, _) in F.fills
            kern(reshape(target(slot, false), :, F.nb), vals, src, dst;
                ndrange = length(src)*F.nb)
        end
    else
        # single precision factors of the equilibrated matrix
        sk = blockscalekernel!(backend, 256)
        sk(F.scale, vals, F.diagidx; ndrange = F.n*F.nb)
        kern = blockfillscaledkernel!(backend, 256)
        for (slot, src, dst, rows, cols) in F.fills
            kern(reshape(target(slot, false), :, F.nb), vals, src, dst, rows,
                cols, F.scale; ndrange = length(src)*F.nb)
        end
    end
    KernelAbstractions.synchronize(backend)
    # the right-looking block LU over the batch: the inverse of each
    # diagonal block, the panel below scaled by it, the panel product
    # subtracted from the later blocks it reaches
    scat = blockschurbatchkernel!(backend, 256)
    for P in 1:F.N
        Dp = F.D[P]; nP = size(Dp, 1)
        X = F.Dinv[P]
        batchedinverse!(X, Dp, F.scratch[(nP, nP)], backend)
        m = size(F.L[P], 1)
        m == 0 && continue
        Lp = F.L[P]
        tmp = F.scratch[(m, nP)]
        batchedmul!(tmp, Lp, X, one(T), zero(T), false, false, backend)
        copyto!(Lp, tmp)
        Wm = F.scratch[(m, m)]
        batchedmul!(Wm, Lp, F.U[P], one(T), zero(T), false, false, backend)
        for t in F.tasks[P]
            scat(t.target, t.rowmap, t.colmap, Wm, t.i0, t.j0;
                ndrange = length(t.rowmap)*length(t.colmap)*F.nb)
        end
    end
    KernelAbstractions.synchronize(backend)
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
        maxpanel = max(maximum(length.(F.rowidx); init = 0), 1)
        alloc = (S, r, c) -> KernelAbstractions.zeros(F.backend, S, r, c, F.nb)
        TA = eltype(F.A)
        w = (; W = Int(W), Z = alloc(T, F.n, W), Y = alloc(T, F.n, W),
            P = alloc(T, maxpanel, W),
            Zo = F.refine ? alloc(TA, F.n, W) : nothing,
            Yo = F.refine ? alloc(TA, F.n, W) : nothing,
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
    backend = F.backend; nb = F.nb
    W = size(B, 2)
    w = blockwork(F, W)
    Z, Y = w.Z, w.Y
    gather = blockgatherrowskernel!(backend, 256)
    scatter = blockscatterrowskernel!(backend, 256)
    scattersub = blockscattersubrowskernel!(backend, 256)
    if isnothing(F.scale)
        gather(Z, B, F.perm; ndrange = F.n*W*nb)
    else
        gathers = blockgatherscaledkernel!(backend, 256)
        gathers(Z, B, F.perm, F.scale; ndrange = F.n*W*nb)
    end
    if !transposed
        for P in 1:F.N
            m = length(F.rowidx[P]); m == 0 && continue
            t = view(w.P, 1:m, :, :)
            batchedmul!(t, F.L[P], view(Z, F.range[P], :, :), one(T), zero(T),
                false, false, backend)
            scattersub(Z, t, F.rowidx[P]; ndrange = m*W*nb)
        end
        for P in F.N:-1:1
            zP = view(Z, F.range[P], :, :)
            m = length(F.rowidx[P])
            if m > 0
                t = view(w.P, 1:m, :, :)
                gather(t, Y, F.rowidx[P]; ndrange = m*W*nb)
                batchedmul!(zP, F.U[P], t, -one(T), one(T), false, false, backend)
            end
            batchedmul!(view(Y, F.range[P], :, :), F.Dinv[P], zP, one(T),
                zero(T), false, false, backend)
        end
    else
        for P in 1:F.N
            zP = view(Z, F.range[P], :, :)
            yP = view(Y, F.range[P], :, :)
            batchedmul!(yP, F.Dinv[P], zP, one(T), zero(T), true, false, backend)
            m = length(F.rowidx[P]); m == 0 && continue
            t = view(w.P, 1:m, :, :)
            batchedmul!(t, F.U[P], yP, one(T), zero(T), true, false, backend)
            scattersub(Z, t, F.rowidx[P]; ndrange = m*W*nb)
        end
        for P in F.N:-1:1
            yP = view(Y, F.range[P], :, :)
            m = length(F.rowidx[P]); m == 0 && continue
            t = view(w.P, 1:m, :, :)
            gather(t, Y, F.rowidx[P]; ndrange = m*W*nb)
            batchedmul!(yP, F.L[P], t, -one(T), one(T), true, false, backend)
        end
    end
    if isnothing(F.scale)
        scatter(X, Y, F.perm; ndrange = F.n*W*nb)
    else
        scatters = blockscatterscaledkernel!(backend, 256)
        scatters(X, Y, F.perm, F.scale; ndrange = F.n*W*nb)
    end
    KernelAbstractions.synchronize(backend)
    return X
end

"""
    blockresidual!(R, F::SparseBlockFactorization, X, B; transposed = false)

`R_k = B - A_k X_k` for the batch from the matrix's own blocks kept for
the refinement, a batched dense product per block, in the matrix's
precision.
"""
function blockresidual!(R::AbstractArray{<:Any,3}, F::SparseBlockFactorization,
    X::AbstractArray{<:Any,3}, B::AbstractMatrix; transposed::Bool = false)
    backend = F.backend; o = F.original; nb = F.nb
    W = size(B, 2); TA = eltype(B)
    w = blockwork(F, W)
    Z, Y, Pn = w.Zo, w.Yo, w.Po
    gather = blockgatherrowskernel!(backend, 256)
    scatter = blockscatterrowskernel!(backend, 256)
    scattersub = blockscattersubrowskernel!(backend, 256)
    gather(Z, X, F.perm; ndrange = F.n*W*nb)
    gather(Y, B, F.perm; ndrange = F.n*W*nb)
    for P in 1:F.N
        zP = view(Z, F.range[P], :, :); yP = view(Y, F.range[P], :, :)
        m = length(F.rowidx[P])
        t = view(Pn, 1:m, :, :)
        batchedmul!(yP, o.D[P], zP, -one(TA), one(TA), transposed, false, backend)
        m == 0 && continue
        gather(t, Z, F.rowidx[P]; ndrange = m*W*nb)
        if !transposed
            batchedmul!(yP, o.U[P], t, -one(TA), one(TA), false, false, backend)
            batchedmul!(t, o.L[P], zP, one(TA), zero(TA), false, false, backend)
        else
            batchedmul!(yP, o.L[P], t, -one(TA), one(TA), true, false, backend)
            batchedmul!(t, o.U[P], zP, one(TA), zero(TA), true, false, backend)
        end
        scattersub(Y, t, F.rowidx[P]; ndrange = m*W*nb)
    end
    scatter(R, Y, F.perm; ndrange = F.n*W*nb)
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
    rnorm = norm(R); bnorm = norm(B)*sqrt(F.nb)
    for _ in 1:F.refinesteps
        rnorm <= 4*eps(real(eltype(B)))*bnorm && break
        blocksolve!(dX, F, R; transposed)
        X .+= dX
        blockresidual!(R, F, X, B; transposed)
        rnew = norm(R)
        rnew < rnorm/2 || break
        rnorm = rnew
    end
    return X
end

# the host forms: one system, matrix or vector right-hand sides
function myldiv!(X::AbstractMatrix, F::SparseBlockFactorization, B::AbstractMatrix)
    refinedsolve!(reshape(X, size(X, 1), size(X, 2), 1), F, B)
    return X
end
function myldiv!(x::AbstractVector, F::SparseBlockFactorization, b::AbstractVector)
    refinedsolve!(reshape(x, :, 1, 1), F, reshape(b, :, 1))
    return x
end
struct TransposedSparseBlockFactorization{F}
    parent::F
end
LinearAlgebra.transpose(F::SparseBlockFactorization) =
    TransposedSparseBlockFactorization(F)
function myldiv!(X::AbstractMatrix, Ft::TransposedSparseBlockFactorization, B::AbstractMatrix)
    refinedsolve!(reshape(X, size(X, 1), size(X, 2), 1), Ft.parent, B; transposed = true)
    return X
end
function myldiv!(x::AbstractVector, Ft::TransposedSparseBlockFactorization, b::AbstractVector)
    refinedsolve!(reshape(x, :, 1, 1), Ft.parent, reshape(b, :, 1); transposed = true)
    return x
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
