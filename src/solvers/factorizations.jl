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
    kluordered(A::SparseMatrixCSC; kwargs...)
    kluordered(A::SparseMatrixCSC, ordering; kwargs...)

`KLU.klu(A)` with its fill reducing ordering chosen by measurement
([`fillordering`](@ref)), or handed in as `ordering`, a permutation of
the columns of `A` or `nothing` for KLU's own. KLU's
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
ordering. Should either ordering fail, KLU's default is used. Everything
before the numeric factorization is symbolic and depends only on the
sparsity pattern: the numeric refactorizations of the pattern reuse the
whole analysis, and a [`FactorizationCache`](@ref) keeps the ordering for
every fresh factorization of it. `kwargs` are `check` and `allowsingular`
of `KLU.klu`.
"""
function kluordered(A::SparseMatrixCSC{Tv,Ti},
    ordering = fillordering(KLUfactorization(), A); check::Bool = true,
    allowsingular::Bool = false) where {Tv,Ti}
    isnothing(ordering) && return KLU.klu(A; check = check, allowsingular = allowsingular)
    nzval = Tv <: Complex ? convert(Vector{ComplexF64}, A.nzval) :
        convert(Vector{Float64}, A.nzval)
    K = KLU.KLUFactorization(size(A, 1), A.colptr .- one(Ti), A.rowval .- one(Ti), nzval)
    p = Ti.(ordering .- 1)
    KLU.klu_analyze!(K, p, copy(p); check = check)
    return KLU.klu_factor!(K; check = check, allowsingular = allowsingular)
end

"""
    fillordering(factorization::AbstractFactorization, A)

The fill reducing ordering `factorization` chooses for the sparsity pattern
of `A`: for a [`KLUfactorization`](@ref) the better of AMD and METIS
nested dissection by predicted flops ([`kluordered`](@ref)), or `nothing`
for KLU's own when neither can be formed or `A` is not square; `nothing`
for any other factorization, which orders a matrix itself. The ordering
depends only on the pattern, so one choice serves every factorization of
it: a [`FactorizationCache`](@ref) keeps the one it chose, and one chosen
elsewhere can be handed to it ([`seedordering!`](@ref)).
"""
fillordering(::AbstractFactorization, A) = nothing
function fillordering(::KLUfactorization, A::SparseMatrixCSC)
    # a matrix which is not square is left to `KLU.klu` to refuse
    size(A, 1) == size(A, 2) || return nothing
    return try _bestordering(A) catch; nothing end
end

# the symmetric pattern of `A`, as CHOLMOD wants it: `A + A'` with unit
# values, 64 bit indices, stored as its upper triangle
function _symmetricpattern(A::SparseMatrixCSC)
    ones_ = SparseMatrixCSC(A.m, A.n, copy(A.colptr), copy(A.rowval), ones(nnz(A)))
    S = SparseMatrixCSC{Float64,Int64}(ones_ + sparse(transpose(ones_)))
    return S
end

# AMD and METIS nested dissection on the symmetric pattern, the one with
# the smaller predicted flop count; `nothing` if neither could be formed
function _bestordering(A::SparseMatrixCSC)
    n = size(A, 1)
    n <= 1 && return nothing
    S = _symmetricpattern(A)
    common = CHOLMOD.getcommon()
    Sc = CHOLMOD.Sparse(S, 1)
    best = nothing
    bestflops = Inf
    for order in (:amd, :metis)
        perm = Vector{Int64}(undef, n)
        ok = if order === :amd
            LibSuiteSparse.cholmod_l_amd(Sc, C_NULL, 0, perm, common)
        else
            LibSuiteSparse.cholmod_l_metis(Sc, C_NULL, 0, true, perm, common)
        end
        ok == 1 || continue
        perm .+= 1
        _, flops = symbolicfill(S, perm)
        if flops < bestflops
            best = perm
            bestflops = flops
        end
    end
    return best
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
    # the fill reducing ordering and the pattern it was chosen for, the
    # `colptr` and `rowval` of the matrix; `nothing` until one is chosen
    ordering
    pattern
end

FactorizationCache(factorization = nothing) =
    FactorizationCache(factorization, nothing, nothing)

"""
    seedordering!(cache::FactorizationCache, A::SparseMatrixCSC, ordering)

Hand `cache` the fill reducing ordering `ordering` of the sparsity pattern
of `A`, as [`fillordering`](@ref) chose it, so that its fresh
factorizations of that pattern take it instead of choosing one: how the
caches of several workers factorizing one pattern share one choice.
Returns `cache`.
"""
function seedordering!(cache::FactorizationCache, A::SparseMatrixCSC, ordering)
    isnothing(ordering) || (length(ordering) == size(A, 2) && isperm(ordering)) ||
        throw(ArgumentError(
            "an ordering is a permutation of the columns of the pattern it is seeded for."))
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
