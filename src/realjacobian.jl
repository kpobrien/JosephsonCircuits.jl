
"""
    RealHolomorphicIndexMap{Ti<:Integer}

Precomputed index map for scattering a complex sparse matrix `M`, treated as a
holomorphic (non-conjugate) term of the Jacobian, directly into the nonzero
vector of a real Jacobian `Jr` whose axes are described by [`ModeLayout`](@ref)s.

For each nonzero `k` of `M` at complex position `(i, j)`, the real form
occupies a `wr` x `wc` block of `Jr` where `wr` and `wc` are 1 for
self-conjugate (real) modes and 2 otherwise. The four fields store the indices
into `nonzeros(Jr)` of the block slots, with 0 marking slots which do not
exist because the row and/or column mode is self-conjugate:

    k11[k] : (r0,   c0  ) receives +real(v)
    k21[k] : (r0+1, c0  ) receives +imag(v), 0 if the row mode is real
    k12[k] : (r0,   c0+1) receives -imag(v), 0 if the column mode is real
    k22[k] : (r0+1, c0+1) receives +real(v), 0 unless both modes are complex

where `v` is the (scaled) value of the nonzero of `M`. This reproduces the
holomorphic `[re -im; im re]` block structure of the equivalent real form
of a complex linear operator (see [`complex_to_real!`](@ref)).

See [`realsparseaddmap`](@ref) and [`realsparseadd!`](@ref).
"""
struct RealHolomorphicIndexMap{Ti<:Integer}
    k11::Vector{Ti}
    k21::Vector{Ti}
    k12::Vector{Ti}
    k22::Vector{Ti}
end

"""
    RealJacobianPlan{Ti<:Integer,T<:Real}

Fully precomputed assembly plan for constructing the real Jacobian of the
harmonic balance nonlinear system directly, without forming the complex
Jacobians `Jx` (holomorphic part) and `Jxconj` (non-holomorphic part) and
converting them to the equivalent real form (the `[re -im; im re]` blocks
of the holomorphic part plus the `[re im; im -re]` blocks of the
non-holomorphic part, with the self-conjugate rows and columns collapsed
to their real representatives).

The nonlinear (Josephson) part of the Jacobian is a fixed linear map from the
real and imaginary parts of the Fourier coefficients of `cos(phi(t))` stored
in `phimatrix` into slots of `nonzeros(Jr)`. That map, which folds together
the scatter of the mode coupling coefficients into the Josephson branch
matrices `AoLjbm` and `AoLjbmconj` (following the difference and sum mode
coupling index matrices), the incidence matrix triple products
`Rbnm'*AoLjbm*Rbnm`, and the complex-to-real block conversion, is stored
as two flat scatter lists:

    nonzeros(Jr)[redest[k]] += recoef[k]*real(phimatrix[resrc[k]])
    nonzeros(Jr)[imdest[k]] += imcoef[k]*imag(phimatrix[imsrc[k]])

with any complex conjugation folded into the sign of `imcoef`. Both the
holomorphic and non-holomorphic (conjugate) Josephson contributions appear in
these lists; the non-holomorphic contribution to self-conjugate (eg. DC)
columns is dropped at plan time, reproducing the `realcolscale_b=0` behavior
of the reference complex-to-real construction. See
[`ComplexJacobianPlan`](@ref) for the analogous plan for the complex
(holomorphic) Jacobian used by the :quasinewton method.

The frequency dependent linear terms are not baked into the plan. Instead the
plan stores [`RealHolomorphicIndexMap`](@ref)s for `invLnm`, `Gnm` and `Cnm`
so they can be scattered at assembly time with the current mode frequencies,
mirroring [`sparseadd!`](@ref) in the complex code path. This keeps the plan
valid when the frequencies change between assemblies.
"""
struct RealJacobianPlan{Ti<:Integer,T<:Real}
    # nonlinear (Josephson) part: scatter from phimatrix
    redest::Vector{Ti}
    resrc::Vector{Ti}
    recoef::Vector{T}
    imdest::Vector{Ti}
    imsrc::Vector{Ti}
    imcoef::Vector{T}
    # frequency dependent linear part: index maps into nonzeros(Jr)
    invLnmindexmap::RealHolomorphicIndexMap{Ti}
    Gnmindexmap::RealHolomorphicIndexMap{Ti}
    Cnmindexmap::RealHolomorphicIndexMap{Ti}
end

"""
    realsparseaddmap(Jr::SparseMatrixCSC, M::SparseMatrixCSC, rl::ModeLayout,
        cl::ModeLayout)

Compute a [`RealHolomorphicIndexMap`](@ref) for scattering the complex sparse
matrix `M` into the real sparse matrix `Jr`, whose row and column axes are the
real forms of the axes of `M` under the layouts `rl` and `cl`. `Ti` selects
the integer type of the stored indices. Every real block slot of every nonzero
of `M` must be a stored entry of `Jr`. This is the real-destination analogue
of [`sparseaddmap`](@ref).
"""
function realsparseaddmap(Jr::SparseMatrixCSC, M::SparseMatrixCSC,
    rl::ModeLayout, cl::ModeLayout, ::Type{Ti} = eltype(rowvals(Jr))) where {Ti<:Integer}

    size(M, 1) == rl.dim || throw(DimensionMismatch(
        lazy"M has $(size(M,1)) rows but the row layout describes $(rl.dim)."))
    size(M, 2) == cl.dim || throw(DimensionMismatch(
        lazy"M has $(size(M,2)) columns but the column layout describes $(cl.dim)."))
    size(Jr, 1) == rl.rdim || throw(DimensionMismatch(
        lazy"Jr has $(size(Jr,1)) rows but the row layout describes $(rl.rdim)."))
    size(Jr, 2) == cl.rdim || throw(DimensionMismatch(
        lazy"Jr has $(size(Jr,2)) columns but the column layout describes $(cl.rdim)."))

    k11 = Vector{Ti}(undef, nnz(M))
    k21 = Vector{Ti}(undef, nnz(M))
    k12 = Vector{Ti}(undef, nnz(M))
    k22 = Vector{Ti}(undef, nnz(M))

    rows = rowvals(M)
    rptr, cptr = rl.ptr, cl.ptr
    @inbounds for j in axes(M, 2)
        c0 = Int(cptr[j])
        wc = Int(cptr[j+1]) - c0
        for k in nzrange(M, j)
            i = rows[k]
            r0 = Int(rptr[i])
            wr = Int(rptr[i+1]) - r0
            k11[k] = storednzindex(Jr, r0, c0)
            k21[k] = wr == 2 ? storednzindex(Jr, r0 + 1, c0) : zero(Ti)
            k12[k] = wc == 2 ? storednzindex(Jr, r0, c0 + 1) : zero(Ti)
            k22[k] = (wr == 2 && wc == 2) ? storednzindex(Jr, r0 + 1, c0 + 1) : zero(Ti)
        end
    end
    return RealHolomorphicIndexMap(k11, k21, k12, k22)
end

"""
    realsparseadd!(Jr::SparseMatrixCSC{T}, c::Number, M::SparseMatrixCSC,
        indexmap::RealHolomorphicIndexMap)
    realsparseadd!(Jr::SparseMatrixCSC{T}, c::Number, M::SparseMatrixCSC,
        D::Diagonal, indexmap::RealHolomorphicIndexMap)

Add the real form of the holomorphic term `c*M` (or `c*M*D` with `D` a
diagonal matrix on the complex column axis, eg. a matrix of mode frequencies)
to the real sparse matrix `Jr` in place using a precomputed
[`RealHolomorphicIndexMap`](@ref). This is the real-destination analogue of
[`sparseadd!`](@ref) and reproduces the holomorphic `[re -im; im re]` block
structure of the equivalent real form.
"""
function realsparseadd!(Jr::SparseMatrixCSC{T}, c::Number, M::SparseMatrixCSC,
    indexmap::RealHolomorphicIndexMap) where {T<:Real}
    nnz(M) == length(indexmap.k11) || throw(DimensionMismatch(
        lazy"The indexmap length $(length(indexmap.k11)) must equal nnz(M) = $(nnz(M))."))
    Jrnz = nonzeros(Jr)
    Mnz = nonzeros(M)
    k11, k21, k12, k22 = indexmap.k11, indexmap.k21, indexmap.k12, indexmap.k22
    @inbounds for k in eachindex(Mnz)
        v = c * Mnz[k]
        re, im_ = real(v), imag(v)
        Jrnz[k11[k]] += re
        iszero(k21[k]) || (Jrnz[k21[k]] += im_)
        iszero(k12[k]) || (Jrnz[k12[k]] -= im_)
        iszero(k22[k]) || (Jrnz[k22[k]] += re)
    end
    return Jr
end

function realsparseadd!(Jr::SparseMatrixCSC{T}, c::Number, M::SparseMatrixCSC,
    D::Diagonal, indexmap::RealHolomorphicIndexMap) where {T<:Real}
    nnz(M) == length(indexmap.k11) || throw(DimensionMismatch(
        lazy"The indexmap length $(length(indexmap.k11)) must equal nnz(M) = $(nnz(M))."))
    size(M, 2) == size(D, 1) || throw(DimensionMismatch(
        lazy"M and D must have compatible sizes."))
    Jrnz = nonzeros(Jr)
    Mnz = nonzeros(M)
    d = D.diag
    k11, k21, k12, k22 = indexmap.k11, indexmap.k21, indexmap.k12, indexmap.k22
    @inbounds for j in axes(M, 2)
        cd = c * d[j]
        for k in nzrange(M, j)
            v = cd * Mnz[k]
            re, im_ = real(v), imag(v)
            Jrnz[k11[k]] += re
            iszero(k21[k]) || (Jrnz[k21[k]] += im_)
            iszero(k12[k]) || (Jrnz[k12[k]] -= im_)
            iszero(k22[k]) || (Jrnz[k22[k]] += re)
        end
    end
    return Jr
end

# number of (re, im) runtime scatter entries produced by one complex-valued
# contribution to a real block of widths (wr, wc). `conjpart` selects the
# non-holomorphic (Jxconj-like) contribution, which is dropped for
# self-conjugate (eg. DC) columns (realcolscale_b = 0 in the
# historical complex-to-real conversion it reproduces).
@inline function countjjentries(wr, wc, conjpart::Bool)
    if conjpart && wc == 1
        return 0, 0
    end
    nre = 1 + ((wr == 2 && wc == 2) ? 1 : 0)
    nim = ((wr == 2) ? 1 : 0) + ((wc == 2) ? 1 : 0)
    return nre, nim
end

# write the runtime scatter entries for one complex-valued contribution
# `coef*phimatrix[src]` (conjugated first if `conjflag`) to the complex
# Jacobian position with real block anchored at (r0, c0) and widths (wr, wc),
# at cursor positions `kre` and `kim` of the preallocated scatter lists.
# `conjpart` selects the non-holomorphic (Jxconj-like) block pattern. Returns
# the updated cursors.
@inline function writejjentries!(redest, resrc, recoef, imdest, imsrc, imcoef,
    kre::Int, kim::Int, Jr::SparseMatrixCSC, r0, c0, wr, wc, src, coef,
    conjflag::Bool, conjpart::Bool)

    # sign applied to imag(phimatrix[src]) from complex conjugation
    s = conjflag ? -one(coef) : one(coef)

    if conjpart && wc == 1
        # the non-holomorphic contribution to self-conjugate (eg. DC) columns
        # is dropped (the non-holomorphic contribution to a self-conjugate
        # column is folded into the holomorphic one, matching the historical
        # convention).
        return kre, kim
    end

    # (r0, c0) += coef*real(v)
    @inbounds begin
        redest[kre] = storednzindex(Jr, r0, c0)
        resrc[kre] = src
        recoef[kre] = coef
        kre += 1
        if wr == 2
            # (r0+1, c0) += coef*imag(v)
            imdest[kim] = storednzindex(Jr, r0 + 1, c0)
            imsrc[kim] = src
            imcoef[kim] = s * coef
            kim += 1
        end
        if wc == 2
            # holomorphic: (r0, c0+1) -= coef*imag(v)
            # non-holomorphic: (r0, c0+1) += coef*imag(v)
            imdest[kim] = storednzindex(Jr, r0, c0 + 1)
            imsrc[kim] = src
            imcoef[kim] = (conjpart ? s : -s) * coef
            kim += 1
            if wr == 2
                # holomorphic: (r0+1, c0+1) += coef*real(v)
                # non-holomorphic: (r0+1, c0+1) -= coef*real(v)
                redest[kre] = storednzindex(Jr, r0 + 1, c0 + 1)
                resrc[kre] = src
                recoef[kre] = (conjpart ? -one(coef) : one(coef)) * coef
                kre += 1
            end
        end
    end
    return kre, kim
end

"""
    planrealjacobian(Amatrixindices::Matrix, Amatrixconjindices::Matrix,
        Ljb::SparseVector, Lmean, Rbnm::SparseMatrixCSC, Nmodes::Integer,
        Nbranches::Integer, Nfreq::Integer, invLnm::SparseMatrixCSC,
        Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, rl::ModeLayout,
        cl::ModeLayout)

Build the real Jacobian sparse matrix `Jr` (with the same sparsity structure
the complex-to-real block conversion would produce from `Jx` and `Jxconj`
(the equivalent real form of `x -> Jx*x + Jxconj*conj(x)`), including
stored numerical zeros) and a [`RealJacobianPlan`](@ref) for assembling it
directly with [`assemblerealjacobian!`](@ref).

The plan folds together, at build time, the maps from the Fourier coefficients
of `cos(phi(t))` to the Josephson branch matrices `AoLjbm` and `AoLjbmconj`
(`Amatrixindices` and `Amatrixconjindices`, with negative entries denoting
complex conjugation and zeros denoting dropped couplings), the incidence matrix
triple products `Rbnm'*AoLjbm*Rbnm` and `Rbnm'*AoLjbmconj*Rbnm`, and the
complex-to-real conversion including the special handling of self-conjugate
(eg. DC) modes. The frequency dependent linear terms `invLnm`, `Gnm` and `Cnm`
are stored as [`RealHolomorphicIndexMap`](@ref)s and scattered at assembly
time so the mode frequencies may change between assemblies.

Returns the tuple `(Jr, plan)`.
"""
function planrealjacobian(Amatrixindices::Matrix, Amatrixconjindices::Matrix,
    Ljb::SparseVector, Lmean, Rbnm::SparseMatrixCSC, Nmodes::Integer,
    Nbranches::Integer, Nfreq::Integer, invLnm::SparseMatrixCSC,
    Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, rl::ModeLayout,
    cl::ModeLayout)

    size(Amatrixindices) == (Nmodes, Nmodes) || throw(DimensionMismatch(
        lazy"Amatrixindices must be Nmodes x Nmodes."))
    size(Amatrixconjindices) == (Nmodes, Nmodes) || throw(DimensionMismatch(
        lazy"Amatrixconjindices must be Nmodes x Nmodes."))
    isreal(Lmean) || throw(ArgumentError(
        "planrealjacobian requires a real Lmean."))
    isempty(Ljb.nzval) || all(isreal, Ljb.nzval) || throw(ArgumentError(
        "planrealjacobian requires real Josephson inductances."))

    # a circuit with no Josephson junctions has an Ljb with element type
    # Nothing and an empty nzval, in which case only the linear terms
    # contribute to the Jacobian.
    T = if isempty(Ljb.nzval) || eltype(Ljb) === Nothing
        real(float(typeof(Lmean)))
    else
        real(promote_type(typeof(Lmean), real(eltype(Ljb))))
    end

    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)

    # The real block widths depend only on the mode index, not the node,
    # because the layout mask repeats every Nmodes entries. Precompute them.
    wrow = [Int(rl.ptr[m+1]) - Int(rl.ptr[m]) for m in 1:Nmodes]
    wcol = [Int(cl.ptr[m+1]) - Int(cl.ptr[m]) for m in 1:Nmodes]

    # Build the sparsity structure of the underlying complex Jacobian (the
    # union of the holomorphic and non-holomorphic Josephson contributions
    # and the linear terms) directly in compressed sparse column form, then
    # expand each complex entry to its real block, keeping numerical zeros as
    # stored entries so the structure does not change between assemblies.
    # This reproduces the structure the complex-to-real conversion of the
    # complex Jacobians would produce.
    n = size(Rbnm, 2)
    adjacency = jjnodeadjacency(Ljb, nodesandsigns, n ÷ Nmodes)
    activem1 = activemoderows(Nmodes, Amatrixindices, Amatrixconjindices)
    ccolptr, crowval = complexjacobianpattern(n, Nmodes, adjacency, activem1,
        (invLnm, Gnm, Cnm))

    # expand the complex structure to the real structure: the real block of
    # the complex entry at (ci, cj) occupies the wr consecutive real rows of
    # ci in each of the wc consecutive real columns of cj. complex rows are
    # sorted within each column and the layout pointer is increasing, so the
    # expanded real rows are sorted as well.
    nrealnz = 0
    for cj in 1:n
        wc = Int(cl.ptr[cj+1]) - Int(cl.ptr[cj])
        rowwidthsum = 0
        for k in ccolptr[cj]:ccolptr[cj+1]-1
            ci = crowval[k]
            rowwidthsum += Int(rl.ptr[ci+1]) - Int(rl.ptr[ci])
        end
        nrealnz += wc * rowwidthsum
    end
    rcolptr = Vector{Int}(undef, cl.rdim + 1)
    rrowval = Vector{Int}(undef, nrealnz)
    rcolptr[1] = 1
    krow = 1
    @inbounds for cj in 1:n
        c0 = Int(cl.ptr[cj])
        wc = Int(cl.ptr[cj+1]) - c0
        for dc in 0:wc-1
            for k in ccolptr[cj]:ccolptr[cj+1]-1
                ci = crowval[k]
                r0 = Int(rl.ptr[ci])
                wr = Int(rl.ptr[ci+1]) - r0
                for dr in 0:wr-1
                    rrowval[krow] = r0 + dr
                    krow += 1
                end
            end
            rcolptr[c0+dc+1] = krow
        end
    end
    krow == nrealnz + 1 || throw(ErrorException(
        "internal error: real structure count mismatch in planrealjacobian."))
    Jr = SparseMatrixCSC(rl.rdim, cl.rdim, rcolptr, rrowval,
        zeros(T, nrealnz))

    # Second pass: emit the runtime scatter entries for the nonlinear part,
    # with destinations resolved to indices into nonzeros(Jr).

    # count the entries so the scatter lists can be preallocated
    nre = 0
    nim = 0
    for i in eachindex(Ljb.nzval)
        npairs = length(nodesandsigns[Ljb.nzind[i]])^2
        for m2 in 1:Nmodes, m1 in 1:Nmodes
            if !iszero(Amatrixindices[m1, m2])
                nre1, nim1 = countjjentries(wrow[m1], wcol[m2], false)
                nre += nre1 * npairs
                nim += nim1 * npairs
            end
            if !iszero(Amatrixconjindices[m1, m2])
                nre1, nim1 = countjjentries(wrow[m1], wcol[m2], true)
                nre += nre1 * npairs
                nim += nim1 * npairs
            end
        end
    end

    # use 32 bit indices for the scatter lists when the destination and
    # source ranges permit, halving the memory traffic of the plan during
    # both construction and assembly. the branch on the index type is
    # resolved before calling the function barrier fillscatterlists so the
    # inner loops are type stable.
    srcmax = Nfreq * length(Ljb.nzval)
    Ti = (nnz(Jr) <= typemax(Int32) && srcmax <= typemax(Int32)) ? Int32 : Int

    redest, resrc, recoef, imdest, imsrc, imcoef = fillscatterlists(Ti, T,
        Jr, Amatrixindices, Amatrixconjindices, Ljb, Lmean, nodesandsigns,
        rl, cl, Nmodes, Nfreq, nre, nim)

    # index maps for the frequency dependent linear terms
    invLnmindexmap = realsparseaddmap(Jr, invLnm, rl, cl, Ti)
    Gnmindexmap = realsparseaddmap(Jr, Gnm, rl, cl, Ti)
    Cnmindexmap = realsparseaddmap(Jr, Cnm, rl, cl, Ti)

    plan = RealJacobianPlan(redest, resrc, recoef, imdest, imsrc, imcoef,
        invLnmindexmap, Gnmindexmap, Cnmindexmap)

    return Jr, plan
end

# function barrier for filling the scatter lists of planrealjacobian with a
# concrete index type Ti, then counting-sorting them by destination. See
# planrealjacobian for the meaning of the arguments.
function fillscatterlists(::Type{Ti}, ::Type{T}, Jr::SparseMatrixCSC,
    Amatrixindices::Matrix, Amatrixconjindices::Matrix, Ljb::SparseVector,
    Lmean, nodesandsigns, rl::ModeLayout, cl::ModeLayout, Nmodes::Integer,
    Nfreq::Integer, nre::Integer, nim::Integer) where {Ti<:Integer,T<:Real}

    redest = Vector{Ti}(undef, nre)
    resrc = Vector{Ti}(undef, nre)
    recoef = Vector{T}(undef, nre)
    imdest = Vector{Ti}(undef, nim)
    imsrc = Vector{Ti}(undef, nim)
    imcoef = Vector{T}(undef, nim)
    kre = 1
    kim = 1

    for i in eachindex(Ljb.nzval)
        b = Ljb.nzind[i]
        ns = nodesandsigns[b]
        LmoLj = T(Lmean / Ljb.nzval[i])
        for m2 in 1:Nmodes, m1 in 1:Nmodes
            ind = Amatrixindices[m1, m2]
            indconj = Amatrixconjindices[m1, m2]
            iszero(ind) && iszero(indconj) && continue
            for (n2, s2) in ns, (n1, s1) in ns
                ci = (n1 - 1) * Nmodes + m1
                cj = (n2 - 1) * Nmodes + m2
                r0 = Int(rl.ptr[ci]); wr = Int(rl.ptr[ci+1]) - r0
                c0 = Int(cl.ptr[cj]); wc = Int(cl.ptr[cj+1]) - c0
                coef = T(s1 * s2) * LmoLj
                if !iszero(ind)
                    src = abs(ind) + Nfreq * (i - 1)
                    kre, kim = writejjentries!(redest, resrc, recoef, imdest,
                        imsrc, imcoef, kre, kim, Jr, r0, c0, wr, wc, src,
                        coef, ind < 0, false)
                end
                if !iszero(indconj)
                    src = abs(indconj) + Nfreq * (i - 1)
                    kre, kim = writejjentries!(redest, resrc, recoef, imdest,
                        imsrc, imcoef, kre, kim, Jr, r0, c0, wr, wc, src,
                        coef, indconj < 0, true)
                end
            end
        end
    end

    (kre == nre + 1 && kim == nim + 1) || throw(ErrorException(
        "internal error: scatter entry count mismatch in planrealjacobian."))

    # sort the scatter lists by destination for memory locality when
    # assembling. the destinations are indices into nonzeros(Jr) so a stable
    # counting sort is O(entries + nnz(Jr)).
    redest, resrc, recoef = sortscatterbydest(redest, resrc, recoef, nnz(Jr))
    imdest, imsrc, imcoef = sortscatterbydest(imdest, imsrc, imcoef, nnz(Jr))

    return redest, resrc, recoef, imdest, imsrc, imcoef
end

"""
    assemblerealjacobian!(Jr::SparseMatrixCSC, plan::RealJacobianPlan,
        phimatrix::Array, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC, wmodesm::Diagonal, wmodes2m::Diagonal)

Assemble the real Jacobian `Jr` in place from the Fourier coefficients of
`cos(phi(t))` in `phimatrix` and the linear term matrices, using a
[`RealJacobianPlan`](@ref) from [`planrealjacobian`](@ref). Computes

    Jr = real(AoLjnm + AoLjnmconj-part + invLnm + im*Gnm*wmodesm - Cnm*wmodes2m)

in the real representation, equivalent to assembling the complex Jacobians
`Jx = Rbnm'*AoLjbm*Rbnm + invLnm + im*Gnm*wmodesm - Cnm*wmodes2m` and
`Jxconj = Rbnm'*AoLjbmconj*Rbnm` and converting the pair to the equivalent
real form of `x -> Jx*x + Jxconj*conj(x)` (dropping the non-holomorphic
contribution to the self-conjugate columns, which the residual absorbs into
the holomorphic part), but with no complex intermediate matrices, no sparse
matrix multiplications, and a single flat scatter loop for the nonlinear
part.
"""
function assemblerealjacobian!(Jr::SparseMatrixCSC, plan::RealJacobianPlan,
    phimatrix::Array, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC, wmodesm::Diagonal, wmodes2m::Diagonal)

    Jrnz = nonzeros(Jr)
    fill!(Jrnz, 0)

    # nonlinear (Josephson) part: direct scatter from phimatrix
    redest, resrc, recoef = plan.redest, plan.resrc, plan.recoef
    @inbounds for k in eachindex(redest)
        Jrnz[redest[k]] += recoef[k] * real(phimatrix[resrc[k]])
    end
    imdest, imsrc, imcoef = plan.imdest, plan.imsrc, plan.imcoef
    @inbounds for k in eachindex(imdest)
        Jrnz[imdest[k]] += imcoef[k] * imag(phimatrix[imsrc[k]])
    end

    # frequency dependent linear part
    realsparseadd!(Jr, 1, invLnm, plan.invLnmindexmap)
    realsparseadd!(Jr, im, Gnm, wmodesm, plan.Gnmindexmap)
    realsparseadd!(Jr, -1, Cnm, wmodes2m, plan.Cnmindexmap)

    return Jr
end
