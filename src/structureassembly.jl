# Assembling the real Jacobian from the circuit's structure, with no plan.
#
# The Josephson part of the Jacobian is the incidence triple product
# `Rbnm' * AoLjbm * Rbnm`, where `AoLjbm` is block diagonal in branches: one
# dense mode block per junction. [`planrealjacobian`](@ref) precomputes that
# product as a segmented gather, one entry per contribution, which on a two
# tone line is 57.6 million entries costing 0.43 s to build and 1.0 GB to read
# at every assembly.
#
# The triple product itself is tiny. A junction touches two nodes, so it
# deposits its mode block at four ordered node pairs with signs: 7996 entries
# for a circuit with 2000 junctions. What the plan stores is the *forward* map,
# contribution to destination. Reading it backwards -- destination to
# contributions -- needs only that table, because an output entry names the
# node pair and the mode pair it belongs to, and the junctions incident on a
# node pair are what contribute to it.

"""
    junctionpairtable(::Type{Ti}, ::Type{T}, Ljb::SparseVector, nodesandsigns,
        Nnodes::Integer)

The junctions incident on each ordered node pair, as a compressed sparse
column structure over the second node: `ptr`, the row `n1`, the junction, and
the product of the two incidence signs.

This is the whole of the incidence triple product, and it has one entry per
(junction, node, node) triple: four per junction for a two terminal one.

The entries of a pair are ordered by junction index, because that is the order
the assembly must sum them in to agree with [`fillscatterlists`](@ref), whose
outer loop is over junctions.
"""
function junctionpairtable(::Type{Ti}, ::Type{T}, Ljb::SparseVector,
    nodesandsigns, Nnodes::Integer) where {Ti<:Integer,T<:Real}

    I = Int[]; J = Int[]; junc = Int[]; sgn = T[]
    for i in eachindex(Ljb.nzval)
        ns = nodesandsigns[Ljb.nzind[i]]
        for (n2, s2) in ns, (n1, s1) in ns
            push!(J, n2); push!(I, n1); push!(junc, i); push!(sgn, T(s1 * s2))
        end
    end
    p = sortperm(collect(zip(J, I, junc)))
    ptr = zeros(Ti, Nnodes + 1)
    ptr[1] = one(Ti)
    for j in J
        ptr[j+1] += one(Ti)
    end
    cumsum!(ptr, ptr)
    return ptr, convert(Vector{Ti}, I[p]), convert(Vector{Ti}, junc[p]), sgn[p]
end

"""
    StructureRealJacobianPlan

Everything the structure aware assembly reads, on a backend. There is no
segmented gather here and nothing proportional to the number of contributions:
the largest arrays are the sparsity structure and the precomputed linear term,
both of which are one entry per *stored entry* of the Jacobian rather than one
per contribution.

# Fields
- `colptr`, `rowval`: the stored structure. The Jacobian is held transposed,
    so a column of it is a row of the Jacobian.
- `lin`: the constant frequency dependent linear contribution.
- `pairptr`, `pairrow`, `pairjunc`, `paircoef`: the incidence triple product,
    from [`junctionpairtable`](@ref).
- `lmolj`: `Lmean/Lj` per junction.
- `perrow`: whether the assembly runs one work item per stored row rather than
    per stored entry, which is the better trade on a host.
- `transposed`: whether the stored structure is the Jacobian's transpose, which
    it is on a device and is not for a matrix meant to be factorized directly.
- `ami`, `amc`: the mode coupling index matrices.
- `rlinv`, `rlptr`, `clinv`, `clptr`: the mode layout and its inverse, which
    turn a stored entry back into a (node, mode) pair.
"""
struct StructureRealJacobianPlan{Ti<:Integer,T<:Real,VI,VT,MI,K,B}
    colptr::VI
    rowval::VI
    lin::VT
    pairptr::VI
    pairrow::VI
    pairjunc::VI
    paircoef::VT
    lmolj::VT
    perrow::Bool
    transposed::Bool
    ami::MI
    amc::MI
    rlinv::VI
    rlptr::VI
    clinv::VI
    clptr::VI
    assemble!::K
    backend::B
    nmodes::Int
    nfreq::Int
    n::Int
end

"""
    realstructureentry(::Type{T}, rri, rci, Nmodes, Nfreq, ami, amc, pairptr,
        pairrow, pairjunc, paircoef, lmolj, rlinv, rlptr, clinv, clptr,
        phimatrix)

The Josephson contribution to the stored entry of the real Jacobian at row
`rri` and column `rci`, which is what both real assembly kernels compute and
the only thing they share.

The entry is decoded to the (node, mode) pair of its row and of its column,
the junctions incident on that node pair are looked up, and their
contributions summed. `(r0, c0)` and `(r0+1, c0+1)` are the real part entries
and `(r0+1, c0)` and `(r0, c0+1)` the imaginary part ones, so a stored entry
belongs to exactly one of the two and only one kind of contribution can reach
it.
"""
@inline function realstructureentry(::Type{T}, rri, rci, Nmodes, Nfreq, ami,
        amc, pairptr, pairrow, pairjunc, paircoef, lmolj, rlinv, rlptr, clinv,
        clptr, phimatrix) where {T}
    @inbounds begin
        ci = Int(rlinv[rri]); r0 = Int(rlptr[ci])
        dr = rri - r0; wr = Int(rlptr[ci+1]) - r0
        n1 = (ci - 1) ÷ Nmodes + 1; m1 = (ci - 1) % Nmodes + 1
        cj = Int(clinv[rci]); c0 = Int(clptr[cj])
        dc = rci - c0; wc = Int(clptr[cj+1]) - c0
        n2 = (cj - 1) ÷ Nmodes + 1; m2 = (cj - 1) % Nmodes + 1

        acc = zero(T)
        ind = Int(ami[m1, m2]); indconj = Int(amc[m1, m2])
        if !(ind == 0 && indconj == 0)
            lo2 = Int(pairptr[n2]); hi2 = Int(pairptr[n2+1]) - 1
            isre = dr == dc
            live = isre ? (dr == 0 || (wr == 2 && wc == 2)) :
                          (dr == 1 ? wr == 2 : wc == 2)
            if live
                for k in lo2:hi2
                    if Int(pairrow[k]) == n1
                        b = Int(pairjunc[k])
                        coef = paircoef[k] * lmolj[b]
                        for which in 1:2
                            v = which == 1 ? ind : indconj
                            if v != 0
                                conjpart = which == 2
                                if !(conjpart && wc == 1)
                                    w = phimatrix[abs(v) + Nfreq * (b - 1)]
                                    if isre
                                        acc += (dr == 1 && conjpart ? -one(T) : one(T)) *
                                            coef * real(w)
                                    else
                                        s = v < 0 ? -one(T) : one(T)
                                        acc += (dr == 1 ? s :
                                                (conjpart ? s : -s)) * coef * imag(w)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
        return acc
    end
end

"""
    structureassemblykernel!(nzval, colptr, rowval, lin, phimatrix, ...)

Assemble one stored entry of the real Jacobian per work item, from the
circuit's structure.

The work item decodes the entry it owns into the node pair and mode pair it
belongs to, looks up the junctions incident on that node pair, and sums their
contributions. The two passes are the real part contributions and then the
imaginary part ones, and the linear term is added last, because that is the
order [`assemblerealjacobian!`](@ref) sums its two segment lists and its
`lin` in, and floating point addition is not associative.
"""
@kernel function structureassemblykernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(lin), @Const(phimatrix), @Const(pairptr), @Const(pairrow),
        @Const(pairjunc), @Const(paircoef), @Const(lmolj), @Const(ami),
        @Const(amc), @Const(rlinv), @Const(rlptr), @Const(clinv),
        @Const(clptr), Nmodes, Nfreq, transposed)

    gid = @index(Global)
    T = eltype(nzval)
    @inbounds begin
        q = gid
        # the stored column this entry is in, which is a row of the Jacobian.
        # A precomputed row index was measured against this and made no
        # difference, so the search stays and the four bytes per entry do not.
        lo = 1
        hi = length(colptr) - 1
        while lo < hi
            mid = (lo + hi + 1) >>> 1
            if colptr[mid] <= q
                lo = mid
            else
                hi = mid - 1
            end
        end
        # the stored matrix is the transpose on a device, so its column is the
        # Jacobian's row; for a normally oriented one it is the other way
        rri = transposed ? lo : Int(rowval[q])
        rci = transposed ? Int(rowval[q]) : lo
        nzval[q] = realstructureentry(T, rri, rci, Nmodes, Nfreq, ami, amc,
            pairptr, pairrow, pairjunc, paircoef, lmolj, rlinv, rlptr, clinv,
            clptr, phimatrix) + lin[q]
    end
end

"""
    structureassemblerowkernel!(nzval, colptr, rowval, lin, phimatrix, ...)

As [`structureassemblykernel!`](@ref), one work item per stored *row* rather
than per stored entry.

The two differ in what they amortize against what they lose. Per entry, the
row has to be found by a binary search and the row side of the decode is
redone for every entry; per row, both are paid once for the whole row, but the
entries a work item writes are contiguous rather than interleaved with its
neighbours', which costs coalescing on a device and nothing on a host.
"""
@kernel function structureassemblerowkernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(lin), @Const(phimatrix), @Const(pairptr), @Const(pairrow),
        @Const(pairjunc), @Const(paircoef), @Const(lmolj), @Const(ami),
        @Const(amc), @Const(rlinv), @Const(rlptr), @Const(clinv),
        @Const(clptr), Nmodes, Nfreq, transposed)

    gid = @index(Global)
    T = eltype(nzval)
    @inbounds begin
        rr = gid
        for q in Int(colptr[rr]):Int(colptr[rr+1])-1
            # the stored matrix is the transpose on a device, so its column is
            # the Jacobian's row; for a normally oriented one it is reversed
            rri = transposed ? rr : Int(rowval[q])
            rci = transposed ? Int(rowval[q]) : rr
            nzval[q] = realstructureentry(T, rri, rci, Nmodes, Nfreq, ami, amc,
                pairptr, pairrow, pairjunc, paircoef, lmolj, rlinv, rlptr,
                clinv, clptr, phimatrix) + lin[q]
        end
    end
end

"""
    planstructurerealjacobian(Jt::SparseMatrixCSC, Amatrixindices,
        Amatrixconjindices, Ljb, Lmean, nodesandsigns, invLnm, Gnm, Cnm,
        wmodesm, wmodes2m, rl::ModeLayout, cl::ModeLayout, Nmodes, Nfreq,
        backend)

Build a [`StructureRealJacobianPlan`](@ref) for the transposed real Jacobian
`Jt`. Nothing here is proportional to the number of contributions, and nothing
needs a [`RealJacobianPlan`](@ref): the constant linear term is gathered on
`backend` by [`linearcontributionkernel!`](@ref) rather than scattered through
index maps on the host.
"""
function planstructurerealjacobian(Jt, ::Type{T},
    Amatrixindices::Matrix, Amatrixconjindices::Matrix,
    Ljb::SparseVector, Lmean, nodesandsigns, invLnm::SparseMatrixCSC,
    Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, wmodesm::Diagonal,
    wmodes2m::Diagonal, rl::ModeLayout, cl::ModeLayout, Nmodes::Integer,
    Nfreq::Integer, backend; transposed::Bool = true) where {T<:Real}

    n = nnz(Jt)
    Ti = n <= typemax(Int32) ? Int32 : Int
    nnodes = cl.dim ÷ Nmodes
    pairptr, pairrow, pairjunc, paircoef =
        junctionpairtable(Ti, T, Ljb, nodesandsigns, nnodes)
    lmolj = T[T(Lmean / Ljb.nzval[i]) for i in eachindex(Ljb.nzval)]

    d = x -> tobackend(backend, convert(Vector{Ti}, x))
    dt = x -> tobackend(backend, convert(Vector{T}, x))
    ami = tobackend(backend, convert(Matrix{Ti}, Amatrixindices))
    amc = tobackend(backend, convert(Matrix{Ti}, Amatrixconjindices))

    # the structure, adopted when it is already on the backend
    dcolptr, drowval = if Jt isa DeviceSparsePattern
        Jt.colptr, Jt.rowval
    else
        d(SparseArrays.getcolptr(Jt)), d(rowvals(Jt))
    end

    # the constant linear term, gathered on the backend
    drlinv = d(collect(rl.inv)); drlptr = d(collect(rl.ptr))
    dclinv = d(collect(cl.inv)); dclptr = d(collect(cl.ptr))
    lin = KernelAbstractions.allocate(backend, T, n)
    linearcontributionkernel!(backend, 64)(lin, dcolptr, drowval,
        drlinv, drlptr, dclinv, dclptr,
        d(SparseArrays.getcolptr(invLnm)), d(rowvals(invLnm)),
        tobackend(backend, nonzeros(invLnm)),
        d(SparseArrays.getcolptr(Gnm)), d(rowvals(Gnm)),
        tobackend(backend, nonzeros(Gnm)), tobackend(backend, collect(wmodesm.diag)),
        d(SparseArrays.getcolptr(Cnm)), d(rowvals(Cnm)),
        tobackend(backend, nonzeros(Cnm)), tobackend(backend, collect(wmodes2m.diag)),
        transposed; ndrange = n)
    KernelAbstractions.synchronize(backend)

    # one work item per entry on a device, where coalescing pays for the
    # repeated decode, and one per row on a host, where it does not
    perrow = backend isa CPU
    assemble! = perrow ? structureassemblerowkernel!(backend, 64) :
        structureassemblykernel!(backend, 64)
    return StructureRealJacobianPlan{Ti,T,typeof(d(pairptr)),typeof(dt(lmolj)),
        typeof(ami),typeof(assemble!),typeof(backend)}(
        dcolptr, drowval, lin,
        d(pairptr), d(pairrow), d(pairjunc), dt(paircoef), dt(lmolj), perrow,
        transposed, ami, amc, drlinv, drlptr, dclinv, dclptr, assemble!, backend,
        Int(Nmodes), Int(Nfreq), n)
end

"""
    assemblerealjacobian!(nzval::AbstractVector,
        plan::StructureRealJacobianPlan, phimatrix::AbstractArray)

Assemble the stored values of the real Jacobian into `nzval` from the Fourier
coefficients of `cos(phi(t))`, using the circuit's structure rather than a
precomputed gather.
"""
function assemblerealjacobian!(nzval::AbstractVector,
    plan::StructureRealJacobianPlan, phimatrix::AbstractArray)
    length(nzval) == plan.n || throw(DimensionMismatch(
        lazy"`nzval` has length $(length(nzval)) but the plan assembles $(plan.n) entries."))
    plan.assemble!(nzval, plan.colptr, plan.rowval, plan.lin, phimatrix,
        plan.pairptr, plan.pairrow, plan.pairjunc, plan.paircoef, plan.lmolj,
        plan.ami, plan.amc, plan.rlinv, plan.rlptr, plan.clinv, plan.clptr,
        plan.nmodes, plan.nfreq, plan.transposed;
        ndrange = plan.perrow ? length(plan.colptr) - 1 : plan.n)
    KernelAbstractions.synchronize(plan.backend)
    return nzval
end

# the value of a sparse matrix at (i, j), or zero, by binary search in its
# column. The linear term matrices are small and their columns short, so this
# is a handful of cached reads.
@inline function sparselookup(colptr, rowval, nzval, i, j)
    lo = Int(colptr[j]); hi = Int(colptr[j+1]) - 1
    z = zero(eltype(nzval))
    lo > hi && return z
    @inbounds while lo < hi
        mid = (lo + hi) >>> 1
        if Int(rowval[mid]) < i
            lo = mid + 1
        else
            hi = mid
        end
    end
    @inbounds return Int(rowval[lo]) == i ? nzval[lo] : z
end

# the real block entry a complex value contributes at offset (dr, dc), which is
# what `realsparseadd!` scatters through its index map
@inline function realblockterm(v, dr, dc)
    if dc == 0
        return dr == 0 ? real(v) : imag(v)
    else
        return dr == 0 ? -imag(v) : real(v)
    end
end

"""
    linearcontributionkernel!(lin, colptr, rowval, ...)

The constant frequency dependent contribution to each stored entry of the
Jacobian: `invLnm + im*Gnm*wmodesm - Cnm*wmodes2m`, in the real
representation.

[`realsparseadd!`](@ref) computes the same thing by scattering each stored
entry of the three matrices through a precomputed index map, which costs a map
per matrix, a pass over a Jacobian sized array per matrix, and, on a device, an
upload of the result. Each stored entry of the Jacobian takes at most one
entry from each of the three, so it can be gathered instead: decode the entry
into the complex position it belongs to and look that position up. The three
are summed in the order `realsparseadd!` is called in, because floating point
addition is not associative.
"""
@kernel function linearcontributionkernel!(lin, @Const(colptr), @Const(rowval),
        @Const(rlinv), @Const(rlptr), @Const(clinv), @Const(clptr),
        @Const(lcolptr), @Const(lrowval), @Const(lnzval),
        @Const(gcolptr), @Const(growval), @Const(gnzval), @Const(wm),
        @Const(ccolptr), @Const(crowval), @Const(cnzval), @Const(wm2),
        transposed)

    gid = @index(Global)
    T = eltype(lin)
    @inbounds begin
        q = gid
        lo = 1
        hi = length(colptr) - 1
        while lo < hi
            mid = (lo + hi + 1) >>> 1
            if colptr[mid] <= q
                lo = mid
            else
                hi = mid - 1
            end
        end
        rri = transposed ? lo : Int(rowval[q])
        rci = transposed ? Int(rowval[q]) : lo
        ci = Int(rlinv[rri]); dr = rri - Int(rlptr[ci])
        cj = Int(clinv[rci]); dc = rci - Int(clptr[cj])

        acc = realblockterm(sparselookup(lcolptr, lrowval, lnzval, ci, cj), dr, dc)
        acc += realblockterm((im * wm[cj]) *
            sparselookup(gcolptr, growval, gnzval, ci, cj), dr, dc)
        acc += realblockterm((-1 * wm2[cj]) *
            sparselookup(ccolptr, crowval, cnzval, ci, cj), dr, dc)
        lin[q] = T(acc)
    end
end

# ---------------------------------------------------------------------------
# the complex Jacobian, from the same structure
# ---------------------------------------------------------------------------
#
# The holomorphic Jacobian is the same incidence triple product without the
# real block expansion, so the same table serves it. It is in fact simpler: a
# stored entry names one mode pair, so `Amatrixindices` at that pair is one
# number, and every contribution to the entry is therefore either conjugated
# or not. There is no equivalent of the real path's two segment lists.

"""
    StructureComplexJosephsonPlan{Ti,T,VI,VT,VR,MI,K,B}

The Josephson contribution to the complex Jacobian, as a linear map from the
Fourier coefficients of `cos(phi(t))` to the stored entries of a matrix with a
given structure.

This replaces the two segmented gathers of [`ComplexJacobianPlan`](@ref), which
stored one entry per contribution. It is used both to assemble the Jacobian and
on its own, as the map applied to other coefficient arrays by the linearized
solve, which is why it is a plan for the Josephson term rather than for the
whole Jacobian.
"""
struct StructureComplexJosephsonPlan{Ti<:Integer,VI,VT,MI,K,B}
    colptr::VI
    rowval::VI
    pairptr::VI
    pairrow::VI
    pairjunc::VI
    paircoef::VT
    lmolj::VT
    perrow::Bool
    transposed::Bool
    ami::MI
    assemble!::K
    backend::B
    nmodes::Int
    nfreq::Int
    n::Int
end

# `conj(sum)` is the sum of conjugated terms, so `conjugate` flips which of the
# two kinds of coupling is conjugated rather than conjugating afterwards
@inline needsconj(ind::Int, conjugate::Bool) = (ind > 0) == conjugate

@inline function josephsonentry(::Type{T}, ci, cj, Nmodes, Nfreq, ami,
        pairptr, pairrow, pairjunc, paircoef, lmolj, phimatrix,
        conjugate) where {T}
    @inbounds begin
        n1 = (ci - 1) ÷ Nmodes + 1; m1 = (ci - 1) % Nmodes + 1
        n2 = (cj - 1) ÷ Nmodes + 1; m2 = (cj - 1) % Nmodes + 1
        ind = Int(ami[m1, m2])
        acc = zero(T)
        ind == 0 && return acc
        for k in Int(pairptr[n2]):Int(pairptr[n2+1])-1
            if Int(pairrow[k]) == n1
                b = Int(pairjunc[k])
                v = phimatrix[abs(ind) + Nfreq * (b - 1)]
                acc += (paircoef[k] * lmolj[b]) *
                    (needsconj(ind, conjugate) ? conj(v) : v)
            end
        end
        return acc
    end
end

@kernel function complexjosephsonkernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(phimatrix), @Const(pairptr), @Const(pairrow), @Const(pairjunc),
        @Const(paircoef), @Const(lmolj), @Const(ami), Nmodes, Nfreq, conjugate,
        transposed)
    gid = @index(Global)
    T = eltype(nzval)
    @inbounds begin
        q = gid
        lo = 1; hi = length(colptr) - 1
        while lo < hi
            mid = (lo + hi + 1) >>> 1
            if colptr[mid] <= q
                lo = mid
            else
                hi = mid - 1
            end
        end
        # a compressed column of the stored structure is a column of the
        # Jacobian, or a row of it when the transpose is what is stored
        ri = transposed ? lo : Int(rowval[q])
        ci = transposed ? Int(rowval[q]) : lo
        nzval[q] = josephsonentry(T, ri, ci, Nmodes, Nfreq, ami,
            pairptr, pairrow, pairjunc, paircoef, lmolj, phimatrix, conjugate)
    end
end

@kernel function complexjosephsonrowkernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(phimatrix), @Const(pairptr), @Const(pairrow), @Const(pairjunc),
        @Const(paircoef), @Const(lmolj), @Const(ami), Nmodes, Nfreq, conjugate,
        transposed)
    gid = @index(Global)
    T = eltype(nzval)
    @inbounds for q in Int(colptr[gid]):Int(colptr[gid+1])-1
        ri = transposed ? gid : Int(rowval[q])
        ci = transposed ? Int(rowval[q]) : gid
        nzval[q] = josephsonentry(T, ri, ci, Nmodes, Nfreq, ami,
            pairptr, pairrow, pairjunc, paircoef, lmolj, phimatrix, conjugate)
    end
end

"""
    planstructurecomplexjosephson(Jx::SparseMatrixCSC, ::Type{T},
        Amatrixindices, Ljb, Lmean, nodesandsigns, Nmodes, Nfreq, backend)

Build a [`StructureComplexJosephsonPlan`](@ref) for the structure of `Jx`.
"""
function planstructurecomplexjosephson(Jx::SparseMatrixCSC, ::Type{T},
    Amatrixindices::Matrix, Ljb::SparseVector, Lmean, nodesandsigns,
    Nmodes::Integer, Nfreq::Integer, backend;
    transposed::Bool = false) where {T<:Real}

    n = nnz(Jx)
    Ti = n <= typemax(Int32) ? Int32 : Int
    nnodes = size(Jx, 2) ÷ Nmodes
    pairptr, pairrow, pairjunc, paircoef =
        junctionpairtable(Ti, T, Ljb, nodesandsigns, nnodes)
    lmolj = T[T(Lmean / Ljb.nzval[i]) for i in eachindex(Ljb.nzval)]
    d = x -> tobackend(backend, convert(Vector{Ti}, x))
    dt = x -> tobackend(backend, convert(Vector{T}, x))
    perrow = backend isa CPU
    assemble! = perrow ? complexjosephsonrowkernel!(backend, 64) :
        complexjosephsonkernel!(backend, 64)
    ami = tobackend(backend, convert(Matrix{Ti}, Amatrixindices))
    return StructureComplexJosephsonPlan{Ti,typeof(d(pairptr)),typeof(dt(lmolj)),
        typeof(ami),typeof(assemble!),typeof(backend)}(
        d(SparseArrays.getcolptr(Jx)), d(rowvals(Jx)), d(pairptr), d(pairrow),
        d(pairjunc), dt(paircoef), dt(lmolj), perrow, transposed, ami,
        assemble!, backend, Int(Nmodes), Int(Nfreq), n)
end

"""
    addjosephsonterm!(nzval::AbstractVector,
        plan::StructureComplexJosephsonPlan, phimatrix, conjugate::Bool = false)

Write the Josephson contribution into `nzval`, overwriting it. With
`conjugate` the map is the one whose source coefficients are conjugated, which
is what the adjoint of the linearized system needs.
"""
function addjosephsonterm!(nzval::AbstractVector,
    plan::StructureComplexJosephsonPlan, phimatrix, conjugate::Bool = false)
    length(nzval) == plan.n || throw(DimensionMismatch(
        lazy"`nzval` has length $(length(nzval)) but the plan assembles $(plan.n) entries."))
    plan.assemble!(nzval, plan.colptr, plan.rowval, phimatrix, plan.pairptr,
        plan.pairrow, plan.pairjunc, plan.paircoef, plan.lmolj, plan.ami,
        plan.nmodes, plan.nfreq, conjugate, plan.transposed;
        ndrange = plan.perrow ? length(plan.colptr) - 1 : plan.n)
    KernelAbstractions.synchronize(plan.backend)
    return nzval
end

"""
    josephsonadjoint!(P, Q, plan::StructureComplexJosephsonPlan,
        w::AbstractVector)

Apply the transpose of the Josephson map: accumulate `w`, which is indexed by
stored entry, back onto the Fourier coefficients it was gathered from.

The plain and conjugated contributions land in `P` and `Q` respectively,
because the caller treats the two halves differently downstream. Which of the
two a stored entry belongs to is decided by the sign of its mode coupling
index, and a stored entry names one mode pair, so an entry contributes to one
of them and never both.

This runs on the host: it is a scatter with collisions, it is used only by the
sensitivity calculation, and that calculation is a host loop throughout.
"""
function josephsonadjoint!(P, Q, plan::StructureComplexJosephsonPlan,
    w::AbstractVector)

    length(w) == plan.n || throw(DimensionMismatch(
        lazy"`w` has length $(length(w)) but the plan has $(plan.n) stored entries."))
    plan.transposed && throw(ArgumentError(
        "josephsonadjoint! reads the plan's stored structure as the Jacobian's own, not as its transpose."))
    colptr, rowval = plan.colptr, plan.rowval
    Nmodes, Nfreq = plan.nmodes, plan.nfreq
    @inbounds for cj in 1:(length(colptr)-1)
        n2 = (cj - 1) ÷ Nmodes + 1
        m2 = (cj - 1) % Nmodes + 1
        for q in Int(colptr[cj]):Int(colptr[cj+1])-1
            ci = Int(rowval[q])
            n1 = (ci - 1) ÷ Nmodes + 1
            m1 = (ci - 1) % Nmodes + 1
            ind = Int(plan.ami[m1, m2])
            ind == 0 && continue
            dst = ind > 0 ? P : Q
            wq = w[q]
            for k in Int(plan.pairptr[n2]):Int(plan.pairptr[n2+1])-1
                if Int(plan.pairrow[k]) == n1
                    b = Int(plan.pairjunc[k])
                    dst[abs(ind) + Nfreq * (b - 1)] +=
                        (plan.paircoef[k] * plan.lmolj[b]) * wq
                end
            end
        end
    end
    return P, Q
end

"""
    complexlinearcontributionkernel!(lin, colptr, rowval, ...)

The constant frequency dependent contribution to each stored entry of the
complex Jacobian: `invLnm + im*Gnm*wmodesm - Cnm*wmodes2m`.

The real path's counterpart splits each complex value across a real block; here
the stored entry is the complex value itself, so the three are simply summed,
in the order `assemblecomplexjacobian!` adds them.
"""
@kernel function complexlinearcontributionkernel!(lin, @Const(colptr),
        @Const(rowval), @Const(lcolptr), @Const(lrowval), @Const(lnzval),
        @Const(gcolptr), @Const(growval), @Const(gnzval), @Const(wm),
        @Const(ccolptr), @Const(crowval), @Const(cnzval), @Const(wm2),
        transposed)
    gid = @index(Global)
    @inbounds begin
        q = gid
        lo = 1; hi = length(colptr) - 1
        while lo < hi
            mid = (lo + hi + 1) >>> 1
            if colptr[mid] <= q
                lo = mid
            else
                hi = mid - 1
            end
        end
        cj = transposed ? Int(rowval[q]) : lo
        ci = transposed ? lo : Int(rowval[q])
        acc = sparselookup(lcolptr, lrowval, lnzval, ci, cj)
        acc += (im * wm[cj]) * sparselookup(gcolptr, growval, gnzval, ci, cj)
        acc += (-1 * wm2[cj]) * sparselookup(ccolptr, crowval, cnzval, ci, cj)
        lin[q] = acc
    end
end

"""
    StructureComplexJacobianPlan{TJ,VT}

The complex Jacobian on a backend: the Josephson map of
[`StructureComplexJosephsonPlan`](@ref) and the constant linear term gathered
once, which is what [`assemblecomplexjacobian!`](@ref) adds to it.
"""
struct StructureComplexJacobianPlan{TJ,VT}
    josephson::TJ
    lin::VT
end

"""
    planstructurecomplexjacobian(Jx::SparseMatrixCSC, ::Type{T},
        Amatrixindices, Ljb, Lmean, nodesandsigns, invLnm, Gnm, Cnm, wmodesm,
        wmodes2m, Nmodes, Nfreq, backend)

Build a [`StructureComplexJacobianPlan`](@ref) for the structure of `Jx`.
"""
function planstructurecomplexjacobian(Jx::SparseMatrixCSC, ::Type{T},
    Amatrixindices::Matrix, Ljb::SparseVector, Lmean, nodesandsigns,
    invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC,
    wmodesm::Diagonal, wmodes2m::Diagonal, Nmodes::Integer, Nfreq::Integer,
    backend; transposed::Bool = false) where {T<:Real}

    josephson = planstructurecomplexjosephson(Jx, T, Amatrixindices, Ljb,
        Lmean, nodesandsigns, Nmodes, Nfreq, backend;
        transposed = transposed)
    n = nnz(Jx)
    Ti = eltype(josephson.pairrow)
    d = x -> tobackend(backend, convert(Vector{Ti}, x))
    lin = KernelAbstractions.allocate(backend, Complex{T}, n)
    complexlinearcontributionkernel!(backend, 64)(lin, josephson.colptr,
        josephson.rowval,
        d(SparseArrays.getcolptr(invLnm)), d(rowvals(invLnm)),
        tobackend(backend, nonzeros(invLnm)),
        d(SparseArrays.getcolptr(Gnm)), d(rowvals(Gnm)),
        tobackend(backend, nonzeros(Gnm)),
        tobackend(backend, collect(wmodesm.diag)),
        d(SparseArrays.getcolptr(Cnm)), d(rowvals(Cnm)),
        tobackend(backend, nonzeros(Cnm)),
        tobackend(backend, collect(wmodes2m.diag)), transposed; ndrange = n)
    KernelAbstractions.synchronize(backend)
    return StructureComplexJacobianPlan{typeof(josephson),typeof(lin)}(
        josephson, lin)
end

"""
    assemblecomplexjacobian!(nzval::AbstractVector,
        plan::StructureComplexJacobianPlan, phimatrix)

Assemble the stored values of the complex Jacobian: the Josephson map applied
to the Fourier coefficients, plus the constant linear term.
"""
function assemblecomplexjacobian!(nzval::AbstractVector,
    plan::StructureComplexJacobianPlan, phimatrix)
    addjosephsonterm!(nzval, plan.josephson, phimatrix)
    nzval .+= plan.lin
    return nzval
end
