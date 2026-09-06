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

`KLU.klu(A)` with its fill reducing ordering chosen by measurement. KLU's
own default, AMD on the pattern of `A + A'`, is the right ordering for
most circuit matrices and a pathological one for some of the mode-coupling
patterns the preconditioners of this package factorize: the harmonic band
of a two-tone line is a mode lattice crossed with the spatial chain, a
grid-like graph, and on the bandwidth-one pattern of a 128-junction line
AMD produced 23 million fill entries and a 20 s factorization where METIS
nested dissection gave 6 million and 0.4 s. Nested dissection is not
uniformly better either: on the full Jacobian of the same line it fills
60% more than AMD and factorizes in twice the time.

So both permutations are computed, AMD and METIS nested dissection, each
through the CHOLMOD library that ships with Julia, the flops of the
factorization each would need are predicted from the elimination tree
([`symbolicfill`](@ref)), and the cheaper one is handed to KLU as a given
ordering. Everything before the numeric factorization is symbolic and
costs a few tenths of a second on a matrix of a million nonzeros, once per
sparsity pattern; the numeric refactorizations of the same pattern reuse
the choice. Should either ordering fail, KLU's default is used.
"""
function kluordered(A::SparseMatrixCSC{Tv,Ti}; check::Bool = true,
    allowsingular::Bool = false) where {Tv,Ti}
    n = size(A, 1)
    # a matrix which is not square is left to `KLU.klu` to refuse
    perm = n == size(A, 2) ? (try _bestordering(A) catch; nothing end) : nothing
    isnothing(perm) && return KLU.klu(A; check = check, allowsingular = allowsingular)
    nzval = Tv <: Complex ? convert(Vector{ComplexF64}, A.nzval) :
        convert(Vector{Float64}, A.nzval)
    K = KLU.KLUFactorization(n, A.colptr .- one(Ti), A.rowval .- one(Ti), nzval)
    p = Ti.(perm .- 1)
    KLU.klu_analyze!(K, p, copy(p); check = check)
    return KLU.klu_factor!(K; check = check, allowsingular = allowsingular)
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

# Refactorize from the nonzero values directly, skipping the sparse matrix
# structure check `KLU.klu!` would do; the pattern is fixed by construction.
function klunzval!(F,A;kwargs...)
    return KLU.klu!(F,A.nzval;kwargs...)
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
    FactorizationCache(factorization)

A mutable holder for a factorization object, so that
[`tryfactorize!`](@ref) can refactorize into it across calls. Starts
empty (`nothing`) when constructed without an argument.

# Examples
```jldoctest
julia> JosephsonCircuits.FactorizationCache(JosephsonCircuits.KLU.klu(JosephsonCircuits.sparse([1, 2], [1, 2], [1/2, 1/2], 2, 2)));

```
"""
mutable struct FactorizationCache
    factorization
end

function FactorizationCache()
    return FactorizationCache(nothing)
end

"""
    tryfactorize!(cache::FactorizationCache,
        factorization::AbstractFactorization, A; kwargs...)

Factorize `A`, a matrix or a [`BlockJacobian`](@ref), with the method
`factorization` and store the result in `cache`. When the cache already
holds a factorization and the method supports refactorization, its symbolic
analysis is reused; a `SingularException` during that refactorization falls
back to a fresh factorization, since reusing the symbolic analysis
occasionally fails numerically where a fresh one succeeds.
`kwargs` are forwarded to `factorize` (the block size of a
[`BlockFactorization`](@ref) of a sparse matrix).
"""
function tryfactorize!(cache::FactorizationCache,
    factorization::AbstractFactorization, A; kwargs...)

    if isnothing(cache.factorization)
        cache.factorization = factorize(factorization, A; kwargs...)
        return cache
    end
    refreshed = try
        # the sparsity structure is unchanged, so refactorize in place; a
        # method without in place refactorization (QR) returns nothing
        refactorize!(factorization, cache.factorization, A)
    catch e
        # reusing the symbolic analysis occasionally fails numerically;
        # factorize afresh
        isa(e, SingularException) || rethrow()
        nothing
    end
    isnothing(refreshed) && (cache.factorization = factorize(factorization, A;
        kwargs...))
    return cache
end


"""
    trysolve!(x,factorization,b)

First try to solve a linear system using ldiv! then if it errors, use \\. The
motivation for this function is some factorizations such as `qr` with sparse
matrices don't support ldiv!. 
"""
function trysolve!(x,factorization,b)
    try
        myldiv!(x,factorization,b)
    catch e
        if e isa MethodError || e isa ArgumentError
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
