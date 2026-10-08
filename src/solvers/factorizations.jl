# The sparse factorizations a direct solve or a preconditioner is built on:
# the KLU ordering chosen by predicted fill, the solve and transposed solve
# over any factorization, and the cache which keeps a symbolic analysis
# between refactorizations.

"""
    symbolicfill(S::SparseMatrixCSC, perm::AbstractVector{<:Integer})

The number of nonzeros of the Cholesky factor of the symmetric pattern `S`
under the symmetric permutation `perm` (`perm[k]` is the original index of
the `k`th pivot), and the flop count of that factorization, both from the
elimination tree without forming anything: `(fill, flops)`.

Only the pattern of `S` is read, and only its structural symmetry matters;
for the pattern of an unsymmetric `A` use that of `A + A'`. The fill of an
LU factorization with the same permutation on both sides is about twice
this and its flops about the same, which is enough to rank two orderings.
Cost `O(fill)`: the column counts are accumulated by walking the row
subtrees of the elimination tree, one step per nonzero of the factor.
"""
function symbolicfill(S::SparseMatrixCSC, perm::AbstractVector{<:Integer})
    n = size(S, 1)
    n == size(S, 2) || throw(DimensionMismatch("the pattern must be square."))
    length(perm) == n || throw(DimensionMismatch(
        lazy"the permutation has length $(length(perm)) but the pattern is $(n) by $(n)."))
    iperm = invperm(perm)
    rows = rowvals(S)
    # the permuted upper pattern, column by column: the row indices below
    # the diagonal of PAP' in the original storage are collected per pivot
    # column, since the elimination tree wants the entries above each
    # pivot and the walk below wants the entries left of each row, and both
    # come from the same list
    counts = zeros(Int, n)
    for j in 1:n, k in nzrange(S, j)
        i = rows[k]
        pi, pj = iperm[i], iperm[j]
        pi < pj && (counts[pj] += 1)
    end
    ptr = Vector{Int}(undef, n + 1)
    ptr[1] = 1
    for j in 1:n
        ptr[j+1] = ptr[j] + counts[j]
    end
    above = Vector{Int}(undef, ptr[end] - 1)
    fill!(counts, 0)
    for j in 1:n, k in nzrange(S, j)
        i = rows[k]
        pi, pj = iperm[i], iperm[j]
        if pi < pj
            above[ptr[pj] + counts[pj]] = pi
            counts[pj] += 1
        end
    end
    # Liu's elimination tree with path compression
    parent = zeros(Int, n)
    ancestor = zeros(Int, n)
    for j in 1:n
        for k in ptr[j]:ptr[j+1]-1
            i = above[k]
            while i != 0 && i < j
                inext = ancestor[i]
                ancestor[i] = j
                if inext == 0
                    parent[i] = j
                    break
                end
                i = inext
            end
        end
    end
    # column counts by row subtrees: row i of L is the union of the paths
    # from each entry left of the diagonal up to i; each node on a path
    # not yet marked for this row is one nonzero of L
    colcount = ones(Int, n)
    mark = zeros(Int, n)
    fillcount = 0
    flops = 0.0
    for i in 1:n
        mark[i] = i
        for k in ptr[i]:ptr[i+1]-1
            j = above[k]
            while j != 0 && mark[j] != i
                mark[j] = i
                colcount[j] += 1
                j = parent[j]
            end
        end
    end
    for j in 1:n
        fillcount += colcount[j]
        flops += float(colcount[j])^2
    end
    return fillcount, flops
end

"""
    FillOrdering(perm, fill)
    FillOrdering(perm, fill, rows)

A fill reducing ordering of a sparsity pattern as [`fillordering`](@ref)
chooses it: the permutation `perm`, `perm[k]` the original index of the
column of the `k`th pivot, `rows` that of its row, `perm` itself unless
the diagonal of the pattern has structural zeros, and `fill`, the entries
of one triangular factor of the pattern under it, diagonal included, as
[`symbolicfill`](@ref) predicts them. [`kluordered`](@ref) sizes the
first allocation of a fresh factorization from `fill`, and
[`sparsefactorbytes`](@ref) the memory of the factors. The choice depends
only on the pattern, so one serves every factorization of it: a
[`FactorizationCache`](@ref) holds the one it was handed or chose
([`seedordering!`](@ref)).
"""
struct FillOrdering
    perm::Vector{Int}
    fill::Int
    rows::Vector{Int}
end
FillOrdering(perm, fill) = FillOrdering(perm, fill, perm)

# the permutation of an ordering handed to a factorization
orderingpermutation(o::FillOrdering) = o.perm

"""
    kluordered(A::SparseMatrixCSC; kwargs...)
    kluordered(A::SparseMatrixCSC, ordering; kwargs...)

`KLU.klu(A)` with its fill reducing ordering chosen by measurement
([`fillordering`](@ref)), or handed in as `ordering`: a
[`FillOrdering`](@ref), or `nothing` for KLU's own. KLU's
own default, AMD on the pattern of `A + A'`, is the right ordering for
most circuit matrices and a pathological one for some of the mode-coupling
patterns the preconditioners of this package factorize: the harmonic band
of a two-tone line is a mode lattice crossed with the spatial chain, a
grid-like graph, on which minimum degree fills many times more than nested
dissection. Nested dissection is not uniformly better either: on the full
Jacobian of the same line it fills more than AMD.

So both permutations are computed, AMD and METIS nested dissection, each
through the CHOLMOD library that ships with Julia, the flops of the
factorization each would need are predicted from the elimination tree
([`symbolicfill`](@ref)), and the cheaper one is handed to KLU as a given
ordering; a harmonic balance sweep has METIS order the graph of its nodes
rather than of its rows ([`fillordering`](@ref)'s `blocksize`). KLU takes
a given ordering without the maximum transversal it finds for its own, so
where the diagonal of `A` has structural zeros, as the Jacobian of a
circuit with scattering blocks can have, a maximum transversal matches its
rows to columns first: the pattern so permuted is ordered, its fill
predicted, and KLU is handed the row and the column permutations, a
zero-free diagonal for its pivots. Should either ordering fail, KLU's
default is used. Everything before the numeric factorization is symbolic
and depends only on the sparsity pattern: the numeric refactorizations of
the pattern reuse the whole analysis, and a [`FactorizationCache`](@ref)
keeps the ordering for every fresh factorization of it.

KLU has no estimate of the fill of an ordering it is handed, and reserves
ten times the entries of `A` for each factor, which it trims once the
factorization is done. A `FillOrdering` carries the fill its choice
predicted, which is the size of each factor under it without pivoting, so
KLU is asked for that, with the margin it gives its own AMD estimate, and
grows a factor which pivoting fills beyond it. `kwargs` are `check` and
`allowsingular` of `KLU.klu`.
"""
function kluordered(A::SparseMatrixCSC{Tv,Ti},
    ordering::Union{Nothing,FillOrdering} = fillordering(KLUfactorization(), A);
    check::Bool = true, allowsingular::Bool = false) where {Tv,Ti}
    isnothing(ordering) && return KLU.klu(A; check = check, allowsingular = allowsingular)
    perm = orderingpermutation(ordering)
    nzval = Tv <: Complex ? convert(Vector{ComplexF64}, A.nzval) :
        convert(Vector{Float64}, A.nzval)
    K = KLU.KLUFactorization(size(A, 1), A.colptr .- one(Ti), A.rowval .- one(Ti), nzval)
    if nnz(A) > 0
        # KLU reserves `initmem*nnz(A) + n` for each factor of a given
        # ordering; `initmem_amd` is the margin it puts on its own estimate
        K.common.initmem = K.common.initmem_amd*ordering.fill/nnz(A)
    end
    KLU.klu_analyze!(K, Ti.(ordering.rows .- 1), Ti.(perm .- 1); check = check)
    return KLU.klu_factor!(K; check = check, allowsingular = allowsingular)
end

"""
    fillordering(factorization::AbstractFactorization, A; blocksize = 1)

The fill reducing ordering `factorization` chooses for the sparsity
pattern of `A`: for a [`KLUfactorization`](@ref) the better of AMD and
METIS nested dissection by predicted flops ([`kluordered`](@ref)), as a
[`FillOrdering`](@ref) carrying its predicted fill, or `nothing` for KLU's
own when neither can be formed or `A` is not square; `nothing` for any
other factorization, which orders a matrix itself. `blocksize` is the
number of rows of each node of a pattern whose rows come in contiguous
node blocks, as a harmonic balance system's modes do: METIS then orders
the graph of the nodes, and the default of one, or a size which does not
divide the rows, has it order the rows. The ordering depends only on the
pattern, so one choice serves every factorization of it: a
[`FactorizationCache`](@ref) keeps the one it chose, and one chosen
elsewhere can be handed to it ([`seedordering!`](@ref)).
"""
fillordering(::AbstractFactorization, A; blocksize::Integer = 1) = nothing
function fillordering(::KLUfactorization, A::SparseMatrixCSC;
    blocksize::Integer = 1)
    # a matrix which is not square is left to `KLU.klu` to refuse
    size(A, 1) == size(A, 2) || return nothing
    # CHOLMOD reports an ordering it cannot form, its memory exhausted
    # included, as a `CHOLMODException`, and KLU's own ordering is taken
    # then; anything else, an interrupt or an error of the code, propagates
    return try
        _bestordering(A; blocksize)
    catch e
        e isa CHOLMOD.CHOLMODException || rethrow()
        nothing
    end
end

# the symmetric pattern of `A`, as CHOLMOD wants it: `A + A'` with unit
# values and 64 bit indices. The Jacobians of this package are
# structurally symmetric, and such a pattern is its own: it is then `A`'s
# index arrays under unit values, shared when they are 64 bit already,
# since nothing here writes them.
function _symmetricpattern(A::SparseMatrixCSC)
    P = SparseMatrixCSC{Float64,Int64}(size(A, 1), size(A, 2),
        convert(Vector{Int64}, SparseArrays.getcolptr(A)),
        convert(Vector{Int64}, rowvals(A)), ones(nnz(A)))
    issymmetric(P) && return P
    return P + sparse(transpose(P))
end

# AMD and METIS nested dissection on the symmetric pattern, the one with
# the smaller predicted flop count with its fill; `nothing` if neither
# could be formed. With `blocksize` rows a node, a size which divides the
# rows, METIS orders the graph of the nodes, each node's rows taken in
# turn: on the graph of the rows its time and memory grow faster than the
# pattern where the node blocks are banded, rows it cannot merge
function _bestordering(A::SparseMatrixCSC; blocksize::Integer = 1)
    n = size(A, 1)
    n <= 1 && return nothing
    # where the diagonal has structural zeros, the pattern with its columns
    # matched to the rows is ordered, and the pivots' columns are the
    # matched ones
    cols = _transversal(A)
    S = _symmetricpattern(isnothing(cols) ? A : A[:, cols])
    common = CHOLMOD.getcommon()
    Sc = CHOLMOD.Sparse(S, 1)
    blocked = blocksize > 1 && n % blocksize == 0
    Gc = blocked ? CHOLMOD.Sparse(_nodegraph(S, blocksize), 1) : Sc
    best = nothing
    bestflops = Inf
    for order in (:amd, :metis)
        perm = Vector{Int64}(undef, order === :metis && blocked ? n ÷ blocksize : n)
        ok = if order === :amd
            LibSuiteSparse.cholmod_l_amd(Sc, C_NULL, 0, perm, common)
        else
            LibSuiteSparse.cholmod_l_metis(Gc, C_NULL, 0, true, perm, common)
        end
        ok == 1 || continue
        perm .+= 1
        order === :metis && blocked && (perm = _rowordering(perm, blocksize))
        fillcount, flops = symbolicfill(S, perm)
        if flops < bestflops
            best = isnothing(cols) ? FillOrdering(perm, fillcount) :
                FillOrdering(cols[perm], fillcount, perm)
            bestflops = flops
        end
    end
    return best
end

# The graph of the nodes of the symmetric pattern `S` whose rows come in
# contiguous blocks of `blocksize`, which divides them: an entry between
# two nodes where a row of one has an entry in a row of the other, read
# off `S` in one pass over its entries.
function _nodegraph(S::SparseMatrixCSC, blocksize::Integer)
    rows = rowvals(S)
    nnodes = size(S, 1) ÷ blocksize
    colptr = Vector{Int64}(undef, nnodes + 1)
    rowval = Int64[]
    mark = zeros(Int, nnodes)
    colptr[1] = 1
    for b in 1:nnodes
        for j in (b - 1)*blocksize + 1:b*blocksize, k in nzrange(S, j)
            a = (rows[k] - 1) ÷ blocksize + 1
            mark[a] == b && continue
            mark[a] = b
            push!(rowval, a)
        end
        colptr[b + 1] = length(rowval) + 1
        sort!(view(rowval, colptr[b]:length(rowval)))
    end
    return SparseMatrixCSC(nnodes, nnodes, colptr, rowval, ones(length(rowval)))
end

# the ordering of the rows from the ordering `perm` of their nodes of
# `blocksize` rows each, each node's rows in turn
_rowordering(perm::Vector{Int64}, blocksize::Integer) =
    vec([(c - 1)*blocksize + a for a in 1:blocksize, c in perm])

# The columns of `A` matched to its rows by a maximum transversal (BTF's,
# as KLU finds it for its own ordering), `cols[i]` the column on the
# diagonal of row `i`, so that `A[:, cols]` has a zero-free diagonal as far
# as the structural rank allows, rows a singular pattern leaves unmatched
# taking the columns left over in order; `nothing` where the diagonal of
# `A` has no structural zero. `maxwork` bounds the search as BTF bounds it,
# in multiples of the entries of `A`, and zero, KLU's default, leaves it
# unbounded.
function _transversal(A::SparseMatrixCSC; maxwork::Real = 0.0)
    n = size(A, 2)
    rows = rowvals(A)
    all(j -> insorted(j, view(rows, nzrange(A, j))), 1:n) && return nothing
    Ap = convert(Vector{Int64}, SparseArrays.getcolptr(A)) .- 1
    Ai = convert(Vector{Int64}, rows) .- 1
    match = Vector{Int64}(undef, n)
    # from the BTF library Julia ships beside KLU's, called by name as
    # KLU.jl calls its own, whose wrapper of BTF names no library
    ccall((:btf_l_maxtrans, :libbtf), Int64, (Int64, Int64, Ptr{Int64}, Ptr{Int64}, Cdouble,
        Ptr{Cdouble}, Ptr{Int64}, Ptr{Int64}), n, n, Ap, Ai, Float64(maxwork), Ref(0.0), match,
        Vector{Int64}(undef, 5n))
    cols = match .+ 1
    matched = falses(n)
    for c in cols
        c > 0 && (matched[c] = true)
    end
    left = findall(!, matched)
    k = 0
    for i in 1:n
        cols[i] > 0 && continue
        k += 1
        cols[i] = left[k]
    end
    return cols
end

"""
    klupivotgrowth(F::KLU.KLUFactorization)

The largest pivot of the KLU factorization `F` in magnitude. KLU factorizes
the matrix with each row scaled to a largest entry of one, so this bounds
the growth of the elimination from below, and it is what a pivot which
has become small inflates: eliminating with a pivot of size `d` adds
entries of size `1/d` to the rows it updates, and to their diagonal when
the pattern is structurally symmetric, as circuit matrices are. `O(n)`,
read from the diagonal of `U` KLU keeps.
"""
function klupivotgrowth(F::KLU.KLUFactorization{Tv}) where {Tv}
    growth = 0.0
    GC.@preserve F begin
        udiag = unsafe_wrap(Array, Ptr{Tv}(F.numeric.Udiag), F.n)
        for u in udiag
            growth = max(growth, abs(u))
        end
    end
    return growth
end

"""
    klurefactor!(F::KLU.KLUFactorization, A::SparseMatrixCSC, pivottol;
        kwargs...)

Refactorize `F` from the values of `A`, whose pattern is the one `F` was
analyzed for, with the pivot sequence of its last factorization
(`KLU.klu!`), and factorize the same values again with fresh partial
pivoting (`KLU.klu_factor!`, which reuses the symbolic analysis) when that
sequence is unstable for them: when a pivot is exactly zero, or when the
pivot growth ([`klupivotgrowth`](@ref)) exceeds `1/pivottol`. `kwargs` are
those of `KLU.klu!`. Returns `F`.
"""
function klurefactor!(F::KLU.KLUFactorization, A::SparseMatrixCSC,
    pivottol::Real; kwargs...)
    stable = try
        # the values straight into the factorization, without the structure
        # check `KLU.klu!` makes of a sparse matrix: the pattern is fixed
        KLU.klu!(F, nonzeros(A); kwargs...)
        iszero(pivottol) || klupivotgrowth(F) <= inv(pivottol)
    catch e
        e isa SingularException || rethrow()
        false
    end
    stable || KLU.klu_factor!(F; kwargs...)
    return F
end

# fallback method to handle everything but QR adjoint
function myldiv!(x,F,b)
    return ldiv!(x,F,b)
end

# workaround for lack of QR adjoint factorization
# modification of solution proposed here
# and https://github.com/JuliaSparse/SparseArrays.jl/issues/656
# explanation: https://www.netlib.org/lapack/lug/node41.html
# and based in part on ldiv! in spqr.jl to handle the rank deficient case
# https://github.com/JuliaSparse/SparseArrays.jl/blob/main/src/solvers/spqr.jl
function myldiv!(x::StridedVecOrMat{<:Number},Fadj::LinearAlgebra.AdjointFactorization{<:Number, <:SparseArrays.SPQR.QRSparse},
                b::StridedVecOrMat{<:Number})
    F = parent(Fadj)
    m, n = size(F)
    pcol = F.pcol
    prow = F.prow
    rnk = rank(F)

    nrhs = size(b, 2)

    # define a workspace
    W = Matrix{eltype(b)}(undef, m, nrhs)

    # c[1:rnk] = first rnk permuted columns of b with the rest of c zero.
    @inbounds for j in 1:nrhs
        for i in 1:rnk
            W[i, j] = b[pcol[i], j]
        end
        for i in (rnk+1):m
            W[i, j] = zero(eltype(W))
        end
    end

    # solve R11^H * y[1:rnk] = c[1:rnk]
    if rnk > 0
        R11 = F.R[1:rnk, 1:rnk]
        ldiv!(UpperTriangular(R11)', view(W, 1:rnk, :))
    end

    # W = Q*y
    lmul!(F.Q, W)

    # get the solution
    @inbounds for j in 1:nrhs
        for i in 1:m
            x[prow[i], j] = W[i, j]
        end
    end

    return x
end

"""
    FactorizationCache(factorization = nothing)

A mutable holder for a factorization object, so that
[`tryfactorize!`](@ref) can refactorize into it across calls, and for the
fill reducing ordering of the sparsity pattern it factorizes, so that a
fresh factorization of the same pattern takes that ordering rather than
choosing again ([`fillordering`](@ref)). Starts empty (`nothing`) when
constructed without an argument. The ordering is chosen by the first fresh
factorization of a pattern, or handed in by [`seedordering!`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.FactorizationCache(JosephsonCircuits.KLU.klu(JosephsonCircuits.sparse([1, 2], [1, 2], [1/2, 1/2], 2, 2)));

```
"""
mutable struct FactorizationCache
    factorization
    # the fill reducing ordering, a `FillOrdering`, and the pattern it was
    # chosen for, the `colptr` and `rowval` of the matrix; `nothing` until
    # one is chosen
    ordering
    pattern
end

FactorizationCache(factorization = nothing) =
    FactorizationCache(factorization, nothing, nothing)

"""
    seedordering!(cache::FactorizationCache, A::SparseMatrixCSC, ordering)

Hand `cache` the fill reducing ordering `ordering` of the sparsity pattern
of `A`, so that its fresh factorizations of that pattern take it instead
of choosing one: how the caches of several workers factorizing one
pattern, or the successive solves of one pattern, share one choice.
`ordering` is what [`fillordering`](@ref) returns, a
[`FillOrdering`](@ref) or `nothing`. The ordering a cache holds,
`cache.ordering`, is one to hand to another. Returns `cache`.
"""
function seedordering!(cache::FactorizationCache, A::SparseMatrixCSC,
    ordering::Union{Nothing,FillOrdering})
    isnothing(ordering) || (length(orderingpermutation(ordering)) == size(A, 2) &&
        isperm(orderingpermutation(ordering)) && length(ordering.rows) == size(A, 1) &&
        isperm(ordering.rows)) ||
        throw(ArgumentError(
            "an ordering is a permutation of the columns and of the rows of the pattern it is seeded for."))
    cache.ordering = ordering
    cache.pattern = (SparseArrays.getcolptr(A), rowvals(A))
    return cache
end

# whether the ordering `cache` holds was chosen for the pattern of `A`: the
# pattern arrays are compared by identity, and by value for a copy
function orderedfor(cache::FactorizationCache, A::SparseMatrixCSC)
    isnothing(cache.pattern) && return false
    colptr, rowval = cache.pattern
    Acolptr, Arowval = SparseArrays.getcolptr(A), rowvals(A)
    return (colptr === Acolptr || colptr == Acolptr) &&
        (rowval === Arowval || rowval == Arowval)
end

"""
    tryfactorize!(cache::FactorizationCache,
        factorization::AbstractFactorization, A; kwargs...)

Factorize `A`, a matrix or a [`BlockJacobian`](@ref), with the method
`factorization` and store the result in `cache`. When the cache already
holds a factorization and the method supports refactorization, its symbolic
analysis is reused; a `SingularException` during that refactorization falls
back to a fresh factorization, since reusing the symbolic analysis can
fail numerically where a fresh one succeeds. (KLU checks the pivots it
reuses itself and repivots in place; see [`klurefactor!`](@ref).)
A fresh factorization by a method which is handed its fill reducing
ordering (KLU) takes the ordering the cache holds for the pattern of `A`,
choosing and keeping one when it holds none. `kwargs` are forwarded to
`factorize` (the block size of a [`BlockFactorization`](@ref) of a sparse
matrix).
"""
function tryfactorize!(cache::FactorizationCache,
    factorization::AbstractFactorization, A; kwargs...)

    if !isnothing(cache.factorization)
        refreshed = try
            # the sparsity structure is unchanged, so refactorize in place;
            # a method without in place refactorization (QR) returns nothing
            refactorize!(factorization, cache.factorization, A)
        catch e
            # reusing the symbolic analysis can fail numerically; factorize
            # afresh
            isa(e, SingularException) || rethrow()
            nothing
        end
        isnothing(refreshed) || return cache
    end
    cache.factorization = freshfactorization!(cache, factorization, A;
        kwargs...)
    return cache
end

# a fresh factorization of `A`, with the ordering the cache holds for its
# pattern when the method is handed one
freshfactorization!(cache::FactorizationCache,
    factorization::AbstractFactorization, A; kwargs...) =
    factorize(factorization, A; kwargs...)
function freshfactorization!(cache::FactorizationCache,
    factorization::KLUfactorization, A::SparseMatrixCSC)
    orderedfor(cache, A) ||
        seedordering!(cache, A, fillordering(factorization, A))
    return kluordered(A, cache.ordering; factorization.kwargs...)
end


"""
    trysolve!(x,factorization,b)

Solve the linear system factorized by `factorization` for the right hand
side `b` into `x` with `ldiv!`, and with `\\` for a factorization which has
no `ldiv!` for these arguments: after a `MethodError` of `ldiv!` itself or
an `ArgumentError`. Any other error is rethrown.
"""
function trysolve!(x,factorization,b)
    try
        myldiv!(x,factorization,b)
    catch e
        if (e isa MethodError && e.f === ldiv!) || e isa ArgumentError
            x .= factorization \ b
        else
            rethrow()
        end
    end
    return x
end

"""
    trysolvetranspose!(x,factorization,b)

Solve the transposed linear system `transpose(A)*x = b` using an existing
factorization of `A`, without refactorizing. The non-conjugating transpose is
used, not the adjoint. Sparse LU factorizations support this directly with a
pair of triangular solves against the stored factors (`klu_tsolve` for KLU), so
an adjoint solve costs a solve rather than a factorization. As in
[`trysolve!`](@ref), fall back to `\\\\` for factorizations which do not support
`ldiv!` with a transposed factorization.

Used by [`hblinsolve`](@ref) to obtain the solutions of the transposed
linearized system, which are the adjoint solutions required by the noise,
quantum efficiency, and sensitivity calculations.
"""
function trysolvetranspose!(x,factorization,b)
    return trysolve!(x,transpose(factorization),b)
end
