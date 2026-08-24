
"""
    storednzindex(A::SparseMatrixCSC, i::Integer, j::Integer)

Return the index into `nonzeros(A)` of the stored entry at position `(i, j)`,
or throw an error if there is no stored entry there. Uses a binary search
within the column, like `SparseArrays` indexing.
"""
function storednzindex(A::SparseMatrixCSC, i::Integer, j::Integer)
    rows = rowvals(A)
    colptr = SparseArrays.getcolptr(A)
    lo, hi = Int(colptr[j]), Int(colptr[j+1]) - 1
    k = searchsortedfirst(view(rows, lo:hi), i) + lo - 1
    if k > hi || rows[k] != i
        throw(ArgumentError(lazy"no stored entry at ($(i), $(j)); the sparsity structure of the destination does not contain the source."))
    end
    return k
end

"""
    branchnodesandsigns(Rbnm::SparseMatrixCSC, Nmodes::Integer,
        Nbranches::Integer)

Recover, for each branch, the list of (node, sign) pairs from the incidence
matrix `Rbnm = diagrepeat(Rbn, Nmodes)` which converts node fluxes to branch
fluxes. Verifies that `Rbnm` has the expected mode-diagonal structure (each
entry connects a branch-mode to a node-mode of the same mode index, with the
same value for every mode) and throws an error otherwise.
"""
function branchnodesandsigns(Rbnm::SparseMatrixCSC, Nmodes::Integer,
    Nbranches::Integer)

    size(Rbnm, 1) == Nbranches * Nmodes || throw(DimensionMismatch(
        lazy"Rbnm has $(size(Rbnm,1)) rows, expected Nbranches*Nmodes = $(Nbranches*Nmodes)."))
    size(Rbnm, 2) % Nmodes == 0 || throw(DimensionMismatch(
        lazy"the number of columns of Rbnm is not a multiple of Nmodes."))

    nodesandsigns = [Tuple{Int,eltype(Rbnm)}[] for _ in 1:Nbranches]
    rows = rowvals(Rbnm)
    vals = nonzeros(Rbnm)
    @inbounds for j in axes(Rbnm, 2)
        n = (j - 1) ÷ Nmodes + 1
        mcol = (j - 1) % Nmodes + 1
        for k in nzrange(Rbnm, j)
            i = rows[k]
            b = (i - 1) ÷ Nmodes + 1
            mrow = (i - 1) % Nmodes + 1
            mrow == mcol || throw(ArgumentError(
                "Rbnm does not have the expected mode-diagonal (diagrepeat) structure."))
            if mcol == 1
                push!(nodesandsigns[b], (n, vals[k]))
            else
                # verify every mode carries the same incidence entry
                any(ns -> ns[1] == n && ns[2] == vals[k], nodesandsigns[b]) ||
                    throw(ArgumentError(
                        "Rbnm does not have the expected mode-diagonal (diagrepeat) structure."))
            end
        end
    end
    return nodesandsigns
end

"""
    segmentscatterbydest(dest::Vector{Ti}, src::Vector{Ti}, coef::Vector{T},
        nbins::Integer)

Group the parallel vectors `(dest, src, coef)`, whose destinations lie in
`1:nbins`, into contiguous segments by destination, returning
`(seg, src, coef)` where the contributions to destination `q` are
`seg[q]:seg[q+1]-1`. This is the compressed sparse row form of the scatter,
so the assembly it describes is a gather: one work item owns destination `q`
and reduces its own segment, with no accumulation into a shared slot and
therefore no atomics.

A stable counting sort by destination, `O(length + nbins)`, which produces
the segment pointer as the prefix sum it needs anyway. The destination
vector is not returned because a gather over `seg` never reads it, which is
why this form is smaller than the scatter it replaces: there are typically
several contributions per destination, so `dest` is the longest array in the
plan.
"""
function segmentscatterbydest(dest::Vector{Ti}, src::Vector{Ti},
    coef::Vector{T}, nbins::Integer) where {Ti<:Integer,T}

    count = zeros(Ti, nbins + 1)
    @inbounds for d in dest
        count[d+1] += 1
    end
    # prefix sum: count[q] is the first output position of destination q, and
    # count[nbins+1] is one past the last, which is exactly the segment
    # pointer the gather needs
    count[1] = 1
    @inbounds for d in 1:nbins
        count[d+1] += count[d]
    end
    isempty(dest) && return count, src, coef
    seg = copy(count)
    src2 = similar(src)
    coef2 = similar(coef)
    # the permutation destroys count, which is why the pointer is copied out
    # of it first
    @inbounds for k in eachindex(dest)
        d = dest[k]
        p = count[d]
        count[d] = p + 1
        src2[p] = src[k]
        coef2[p] = coef[k]
    end
    return seg, src2, coef2
end

"""
    jjnodeadjacency(Ljb::SparseVector, nodesandsigns, Nnodes::Integer)

For each node, return the sorted, deduplicated list of nodes which share at
least one Josephson junction branch with it (including itself), from the
per-branch (node, sign) lists in `nodesandsigns` (see
[`branchnodesandsigns`](@ref)) restricted to the Josephson branches in `Ljb`.
Used to enumerate the columns of the Jacobian sparsity structure directly in
compressed sparse column form.
"""
function jjnodeadjacency(Ljb::SparseVector, nodesandsigns, Nnodes::Integer)
    adjacency = [Int[] for _ in 1:Nnodes]
    for i in eachindex(Ljb.nzind)
        ns = nodesandsigns[Ljb.nzind[i]]
        for (n2, _) in ns, (n1, _) in ns
            push!(adjacency[n2], n1)
        end
    end
    for a in adjacency
        sort!(a)
        unique!(a)
    end
    return adjacency
end

"""
    activemoderows(Nmodes::Integer, Amatrixindices::Matrix,
        Bmatrixindices::Union{Matrix,Nothing} = nothing)

For each column mode, return the sorted list of row modes for which the
frequency domain index matrix `Amatrixindices` (or, if provided,
`Bmatrixindices`) has a nonzero entry, ie. the row modes which contribute to
the Jacobian sparsity structure in that column mode.
"""
function activemoderows(Nmodes::Integer, Amatrixindices::Matrix,
    Bmatrixindices::Union{Matrix,Nothing} = nothing)
    return [[m1 for m1 in 1:Nmodes if !iszero(Amatrixindices[m1, m2]) ||
        (!isnothing(Bmatrixindices) && !iszero(Bmatrixindices[m1, m2]))]
        for m2 in 1:Nmodes]
end

"""
    complexjacobianpattern(n::Integer, Nmodes::Integer, adjacency, activem1,
        linearmatrices)

Build the sparsity structure (colptr, rowval) of the complex Jacobian of the
harmonic balance nonlinear system directly in compressed sparse column form:
the union of the Josephson contributions, whose rows in the column of node
`n2` and mode `m2` are the modes `activem1[m2]` (see
[`activemoderows`](@ref)) of the nodes `adjacency[n2]` (see
[`jjnodeadjacency`](@ref)), and of the sparsity structures of the linear term
matrices. Numerical zeros are kept as stored entries so the structure does
not change between assemblies. This reproduces the structure `sparse` would
produce from the corresponding list of coordinates, or `spaddkeepzeros`
applied to `Rbnm'*AoLjbm*Rbnm` and the linear term matrices, without forming
either.
"""
function complexjacobianpattern(n::Integer, Nmodes::Integer, adjacency,
    activem1, linearmatrices)

    for M in linearmatrices
        size(M) == (n, n) || throw(DimensionMismatch(
            lazy"the linear term matrices must be $(n) x $(n)."))
    end

    # count and collect the rows of each column: mark the candidate rows in
    # a workspace to deduplicate, collect, and sort. the workspace marks are
    # cleared by iterating the collected rows so the cost is proportional to
    # the number of entries, not the dimension.
    seen = fill(false, n)
    colrows = Int[]
    colptr = Vector{Int}(undef, n + 1)
    colptr[1] = 1
    rowval = Int[]
    for cj in 1:n
        n2 = (cj - 1) ÷ Nmodes + 1
        m2 = (cj - 1) % Nmodes + 1
        empty!(colrows)
        # Josephson contributions
        for n1 in adjacency[n2], m1 in activem1[m2]
            ci = (n1 - 1) * Nmodes + m1
            if !seen[ci]
                seen[ci] = true
                push!(colrows, ci)
            end
        end
        # linear contributions
        for M in linearmatrices
            rows = rowvals(M)
            for k in nzrange(M, cj)
                ci = Int(rows[k])
                if !seen[ci]
                    seen[ci] = true
                    push!(colrows, ci)
                end
            end
        end
        sort!(colrows)
        append!(rowval, colrows)
        for ci in colrows
            seen[ci] = false
        end
        colptr[cj+1] = colptr[cj] + length(colrows)
    end
    return colptr, rowval
end

"""
    ComplexJacobianPlan{Ti<:Integer,T<:Real}

Fully precomputed assembly plan for constructing the complex (holomorphic)
Jacobian `Jx` of the harmonic balance nonlinear system directly, without
forming the Josephson branch matrix `AoLjbm` and multiplying by the incidence
matrices each iteration.

The nonlinear (Josephson) part of the Jacobian is a fixed linear map from the
Fourier coefficients of `cos(phi(t))` stored in `phimatrix` into slots of
`nonzeros(Jx)`. That map, which folds together the scatter of the mode coupling
coefficients into the Josephson branch matrix `AoLjbm` (following
`Amatrixindices`, with negative entries denoting complex conjugation and
zeros denoting dropped couplings, scaled by `Lmean/Lj`) and the incidence
matrix triple product `Rbnm'*AoLjbm*Rbnm` ([`spmatmul!`](@ref)), is stored
as two segmented gather lists ([`segmentscatterbydest`](@ref)), one for the
plain and one for the complex conjugated coefficients, so that nonzero `q`
of `Jx` is

    sum(coef[t]*phimatrix[src[t]] for t in seg[q]:seg[q+1]-1) +
    sum(ccoef[t]*conj(phimatrix[csrc[t]]) for t in cseg[q]:cseg[q+1]-1)

One work item owns one nonzero and reduces its own segments, so the assembly
needs no accumulation into shared slots and no atomics.

The frequency dependent linear terms are not baked into the plan. Instead the
plan stores [`sparseaddmap`](@ref) index maps for `invLnm`, `Gnm` and `Cnm`
so they can be scattered at assembly time with the current mode frequencies
using [`sparseadd!`](@ref). This keeps the plan valid when the frequencies
change between assemblies.

See also [`plancomplexjacobian`](@ref), [`assemblecomplexjacobian!`](@ref)
and, for the analogous plan for the exact real Jacobian,
[`RealJacobianPlan`](@ref).
"""
struct ComplexJacobianPlan{Ti<:Integer,T<:Real}
    # nonlinear (Josephson) part: segmented gather from phimatrix
    seg::Vector{Ti}
    src::Vector{Ti}
    coef::Vector{T}
    # as above but with the source coefficient complex conjugated
    cseg::Vector{Ti}
    csrc::Vector{Ti}
    ccoef::Vector{T}
    # frequency dependent linear part: index maps into nonzeros(Jx)
    invLnmindexmap::Vector{Int}
    Gnmindexmap::Vector{Int}
    Cnmindexmap::Vector{Int}
end

# function barrier for filling the scatter lists of plancomplexjacobian with
# a concrete index type Ti, then counting-sorting them by destination. See
# plancomplexjacobian for the meaning of the arguments.
function fillcomplexscatterlists(::Type{Ti}, ::Type{T}, Jx::SparseMatrixCSC,
    Amatrixindices::Matrix, Ljb::SparseVector, Lmean, nodesandsigns,
    Nmodes::Integer, Nfreq::Integer, nplain::Integer,
    nconj::Integer) where {Ti<:Integer,T<:Real}

    dest = Vector{Ti}(undef, nplain)
    src = Vector{Ti}(undef, nplain)
    coef = Vector{T}(undef, nplain)
    cdest = Vector{Ti}(undef, nconj)
    csrc = Vector{Ti}(undef, nconj)
    ccoef = Vector{T}(undef, nconj)
    kplain = 1
    kconj = 1

    @inbounds for i in eachindex(Ljb.nzval)
        b = Ljb.nzind[i]
        ns = nodesandsigns[b]
        LmoLj = T(Lmean / Ljb.nzval[i])
        for m2 in 1:Nmodes, m1 in 1:Nmodes
            ind = Amatrixindices[m1, m2]
            iszero(ind) && continue
            fsrc = abs(ind) + Nfreq * (i - 1)
            for (n2, s2) in ns, (n1, s1) in ns
                ci = (n1 - 1) * Nmodes + m1
                cj = (n2 - 1) * Nmodes + m2
                k = storednzindex(Jx, ci, cj)
                if ind > 0
                    dest[kplain] = k
                    src[kplain] = fsrc
                    coef[kplain] = T(s1 * s2) * LmoLj
                    kplain += 1
                else
                    cdest[kconj] = k
                    csrc[kconj] = fsrc
                    ccoef[kconj] = T(s1 * s2) * LmoLj
                    kconj += 1
                end
            end
        end
    end

    (kplain == nplain + 1 && kconj == nconj + 1) || throw(ErrorException(
        "internal error: scatter entry count mismatch in plancomplexjacobian."))

    # group the contributions into contiguous segments by destination, which
    # turns the assembly into a gather. the destinations are indices into
    # nonzeros(Jx) so a stable counting sort is O(entries + nnz(Jx)), and it
    # produces the segment pointer as its prefix sum.
    seg, src, coef = segmentscatterbydest(dest, src, coef, nnz(Jx))
    cseg, csrc, ccoef = segmentscatterbydest(cdest, csrc, ccoef, nnz(Jx))

    return seg, src, coef, cseg, csrc, ccoef
end

"""
    plancomplexjacobian(Amatrixindices::Matrix, Ljb::SparseVector, Lmean,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC)

Build the complex Jacobian sparse matrix `Jx` (with the same sparsity
structure `spaddkeepzeros` applied to `Rbnm'*AoLjbm*Rbnm` and the linear term
matrices would produce, including stored numerical zeros) and a
[`ComplexJacobianPlan`](@ref) for assembling it directly with
[`assemblecomplexjacobian!`](@ref).

The plan folds together, at build time, the map from the Fourier coefficients
of `cos(phi(t))` to the Josephson branch matrix `AoLjbm` (`Amatrixindices`,
with negative entries denoting complex conjugation and zeros denoting
dropped couplings) and the incidence matrix triple product
`Rbnm'*AoLjbm*Rbnm`. The frequency dependent linear terms `invLnm`, `Gnm` and
`Cnm` are stored as [`sparseaddmap`](@ref) index maps and scattered at
assembly time so the mode frequencies may change between assemblies.

Returns the tuple `(Jx, plan)`.
"""
function plancomplexjacobian(Amatrixindices::Matrix, Ljb::SparseVector,
    Lmean, Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC)

    size(Amatrixindices) == (Nmodes, Nmodes) || throw(DimensionMismatch(
        lazy"Amatrixindices must be Nmodes x Nmodes."))
    isreal(Lmean) || throw(ArgumentError(
        "plancomplexjacobian requires a real Lmean."))
    isempty(Ljb.nzval) || all(isreal, Ljb.nzval) || throw(ArgumentError(
        "plancomplexjacobian requires real Josephson inductances."))

    # a circuit with no Josephson junctions has an Ljb with element type
    # Nothing and an empty nzval, in which case only the linear terms
    # contribute to the Jacobian.
    T = if isempty(Ljb.nzval) || eltype(Ljb) === Nothing
        real(float(typeof(Lmean)))
    else
        real(promote_type(typeof(Lmean), real(eltype(Ljb))))
    end
    Tc = Complex{float(T)}

    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)

    # count the runtime scatter entries so the lists can be preallocated
    nplain = 0
    nconj = 0
    for i in eachindex(Ljb.nzval)
        npairs = length(nodesandsigns[Ljb.nzind[i]])^2
        for m2 in 1:Nmodes, m1 in 1:Nmodes
            ind = Amatrixindices[m1, m2]
            iszero(ind) && continue
            if ind > 0
                nplain += npairs
            else
                nconj += npairs
            end
        end
    end

    # Build the sparsity structure (the union of the Josephson contributions
    # and the linear terms, keeping numerical zeros as stored entries so the
    # structure does not change between assemblies) directly in compressed
    # sparse column form. This reproduces the structure produced by
    # spaddkeepzeros on Rbnm'*AoLjbm*Rbnm and the linear term matrices.
    n = size(Rbnm, 2)
    adjacency = jjnodeadjacency(Ljb, nodesandsigns, n ÷ Nmodes)
    activem1 = activemoderows(Nmodes, Amatrixindices)
    colptr, rowval = complexjacobianpattern(n, Nmodes, adjacency, activem1,
        (invLnm, Gnm, Cnm))
    Jx = SparseMatrixCSC(n, n, colptr, rowval, zeros(Tc, length(rowval)))

    # Second pass: emit the runtime scatter entries for the nonlinear part,
    # with destinations resolved to indices into nonzeros(Jx). use 32 bit
    # indices for the scatter lists when the destination and source ranges
    # permit, halving the memory traffic of the plan during both construction
    # and assembly. the branch on the index type is resolved before calling
    # the function barrier fillcomplexscatterlists so the inner loops are
    # type stable.
    srcmax = Nfreq * length(Ljb.nzval)
    Ti = (nnz(Jx) <= typemax(Int32) && srcmax <= typemax(Int32)) ? Int32 : Int

    seg, src, coef, cseg, csrc, ccoef = fillcomplexscatterlists(Ti, T, Jx,
        Amatrixindices, Ljb, Lmean, nodesandsigns, Nmodes, Nfreq, nplain,
        nconj)

    # index maps for the frequency dependent linear terms
    invLnmindexmap = sparseaddmap(Jx, invLnm)
    Gnmindexmap = sparseaddmap(Jx, Gnm)
    Cnmindexmap = sparseaddmap(Jx, Cnm)

    plan = ComplexJacobianPlan(seg, src, coef, cseg, csrc, ccoef,
        invLnmindexmap, Gnmindexmap, Cnmindexmap)

    return Jx, plan
end

"""
    assemblecomplexjacobian!(Jx::SparseMatrixCSC, plan::ComplexJacobianPlan,
        phimatrix::Array, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC, wmodesm::Diagonal, wmodes2m::Diagonal)

Assemble the complex (holomorphic) Jacobian `Jx` in place from the Fourier
coefficients of `cos(phi(t))` in `phimatrix` and the linear term matrices,
using a [`ComplexJacobianPlan`](@ref) from [`plancomplexjacobian`](@ref).
Computes

    Jx = Rbnm'*AoLjbm*Rbnm + invLnm + im*Gnm*wmodesm - Cnm*wmodes2m

with no branch matrix update and no sparse matrix multiplications, just a
flat scatter loop for the nonlinear part and index mapped additions
([`sparseadd!`](@ref)) for the linear part.
"""
function assemblecomplexjacobian!(Jx::SparseMatrixCSC,
    plan::ComplexJacobianPlan, phimatrix::Array, invLnm::SparseMatrixCSC,
    Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, wmodesm::Diagonal,
    wmodes2m::Diagonal)

    Jxnz = nonzeros(Jx)

    # nonlinear (Josephson) part: a segmented gather from phimatrix, which
    # writes every nonzero, so no zeroing pass is needed
    addjosephsonterm!(Jxnz, plan, phimatrix)

    # frequency dependent linear part
    sparseadd!(Jx, invLnm, plan.invLnmindexmap)
    sparseadd!(Jx, im, Gnm, wmodesm, plan.Gnmindexmap)
    sparseadd!(Jx, -1, Cnm, wmodes2m, plan.Cnmindexmap)

    return Jx
end

"""
    addjosephsonterm!(nzval::AbstractVector, plan::ComplexJacobianPlan,
        phimatrix::Array, conjugate::Bool = false)

Write the Josephson (nonlinear) contribution `Rbnm'*AoLjbm*Rbnm` into the
nonzero value vector `nzval`, which must be aligned with the sparsity
structure the [`ComplexJacobianPlan`](@ref) was built for, by gathering the
Fourier coefficients of `cos(phi(t))` in `phimatrix` through the plan. Every
entry of `nzval` is written, so the destination need not be zeroed first.
With `conjugate = true` the elementwise complex conjugate of the contribution
is written instead, as needed for the adjoint of the linearized harmonic
balance system. Used by [`assemblecomplexjacobian!`](@ref) for the Jacobian of the
nonlinear system and by [`hblinsolve`](@ref) for the pump modulation term of
the linearized system, so the two constructions share the same machinery.
"""
function addjosephsonterm!(nzval::AbstractVector, plan::ComplexJacobianPlan,
    phimatrix::Array, conjugate::Bool = false)

    seg, src, coef = plan.seg, plan.src, plan.coef
    cseg, csrc, ccoef = plan.cseg, plan.csrc, plan.ccoef
    # conj(sum) = sum of conjugated terms: with conjugate = true the plain
    # coefficients are conjugated and the conjugated ones are not
    T = eltype(nzval)
    @inbounds for q in eachindex(nzval)
        acc = zero(T)
        for t in seg[q]:seg[q+1]-1
            v = phimatrix[src[t]]
            acc += coef[t] * (conjugate ? conj(v) : v)
        end
        for t in cseg[q]:cseg[q+1]-1
            v = phimatrix[csrc[t]]
            acc += ccoef[t] * (conjugate ? v : conj(v))
        end
        nzval[q] = acc
    end
    return nzval
end
