# Assembling the real Jacobian from the circuit's structure.
#
# The Josephson part of the Jacobian is the incidence triple product
# `Rbnm' * AoLjbm * Rbnm`, where `AoLjbm` is block diagonal in branches: one
# dense mode block per junction. Precomputing that product as a gather with
# one entry per contribution is what an assembly plan would ordinarily do,
# and on a multi tone line that gather is larger than the Jacobian it
# fills, and is read at every assembly.
#
# The triple product itself is tiny. A junction touches two nodes, so it
# deposits its mode block at four ordered node pairs with signs, a table of
# a few entries per junction (see `junctionpairtable`). That table is
# enough on its own, because it can be read backwards, from destination to
# contributions: an output entry names the node pair and the mode pair it
# belongs to, and the junctions incident on that node pair are exactly what
# contribute to it. So the assembly here stores the structure and no gather
# at all.

"""
    junctionpairtable(::Type{Ti}, ::Type{T}, Ljb::SparseVector, nodesandsigns,
        Nnodes::Integer)

The junctions incident on each ordered node pair, as a compressed sparse
column structure over the second node: `ptr`, the row `n1`, the junction, and
the product of the two incidence signs.

This is the whole of the incidence triple product, and it has one entry per
(junction, node, node) triple: four per junction for a two terminal one.

The entries of a pair are ordered by junction index, which fixes the order the
assembly sums them in. That order is part of the result, because floating
point addition is not associative.
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
    JunctionStructure

What every assembly of the Josephson term reads about the junctions, built
once per system and precision and shared by the plans which assemble from
it: the incidence triple product as a per node pair table, the coefficient
`Lscale/Lj` of each junction, and the mode coupling index matrices, on the
backend. The values which move with the component values, `lmolj`, are
refreshed in place by [`refreshvalues!`](@ref), so every plan holding the
structure sees the new values at once; the table and the index matrices
are fixed by the topology and the mode grid.

# Fields
- `pairptr`, `pairrow`, `pairjunc`, `paircoef`: the incidence triple product
    per ordered node pair, from [`junctionpairtable`](@ref), on the backend;
    `hostpairptr`, `hostpairrow`: the pointer and row of the table on the
    host, for the node graph of a block factorization.
- `lmolj`: `Lscale/Lj` per junction, on the backend.
- `ami`, `amc`: the mode coupling index matrices, on the backend.
- `nodesandsigns`: the (node, sign) pairs of each branch, on the host.
- `nmodes`, `nfreq`: the mode count and the frequency grid stride into
    `phimatrix`.
"""
struct JunctionStructure{T<:Real,VI,VT,MI}
    pairptr::VI
    pairrow::VI
    pairjunc::VI
    paircoef::VT
    lmolj::VT
    ami::MI
    amc::MI
    hostpairptr::Vector{Int32}
    hostpairrow::Vector{Int32}
    nodesandsigns::Vector{Vector{Tuple{Int,Float64}}}
    nmodes::Int
    nfreq::Int
end

"""
    junctionstructure(::Type{T}, Amatrixindices::Matrix,
        Amatrixconjindices::Matrix, Ljb::SparseVector, Lscale,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, backend)

Build the [`JunctionStructure`](@ref) of a system in precision `T` on
`backend`. The pair table is indexed in `Int32`, which the node and junction
counts never exceed; a plan's own structure keeps whatever index type its
size needs.
"""
function junctionstructure(::Type{T}, Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, Ljb::SparseVector, Lscale,
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, backend) where {T<:Real}
    size(Amatrixindices) == (Nmodes, Nmodes) || throw(DimensionMismatch(
        lazy"Amatrixindices must be Nmodes x Nmodes."))
    size(Amatrixconjindices) == (Nmodes, Nmodes) || throw(DimensionMismatch(
        lazy"Amatrixconjindices must be Nmodes x Nmodes."))
    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
    nnodes = size(Rbnm, 2) ÷ Nmodes
    pairptr, pairrow, pairjunc, paircoef =
        junctionpairtable(Int32, T, Ljb, nodesandsigns, nnodes)
    lmolj = junctioncoefficients(T, Ljb, Lscale)
    d = x -> tobackend(backend, Vector{Int32}(x))
    dt = x -> tobackend(backend, Vector{T}(x))
    ami = tobackend(backend, Matrix{Int32}(Amatrixindices))
    amc = tobackend(backend, Matrix{Int32}(Amatrixconjindices))
    ns = [Tuple{Int,Float64}[(Int(n), Float64(s)) for (n, s) in b] for b in nodesandsigns]
    # each array is moved to the backend once, and the structure's type read
    # from the moved arrays
    ptrd = d(pairptr)
    lmoljd = dt(lmolj)
    return JunctionStructure{T,typeof(ptrd),typeof(lmoljd),typeof(ami)}(ptrd,
        d(pairrow), d(pairjunc), dt(paircoef), lmoljd, ami, amc,
        Vector{Int32}(pairptr), Vector{Int32}(pairrow), ns, Int(Nmodes),
        Int(Nfreq))
end

"""
    refreshvalues!(js::JunctionStructure, Ljb::SparseVector, Lscale)

Rewrite the junction coefficients `Lscale/Lj` of a structure for new
component values; the table and the index matrices are what they were. The
junctions must be the same ones.
"""
function refreshvalues!(js::JunctionStructure{T}, Ljb::SparseVector,
        Lscale) where {T}
    length(Ljb.nzval) == length(js.lmolj) || throw(ArgumentError(
        "the junctions changed between points, which a refreshed structure cannot follow; build a new one."))
    copyto!(js.lmolj, junctioncoefficients(T, Ljb, Lscale))
    return js
end

# the coefficient `Lscale/Lj` of each junction in precision `T`, the
# quotient rounded once, which every plan holding one and every refresh of
# it writes
junctioncoefficients(::Type{T}, Ljb::SparseVector, Lscale) where {T} =
    T[T(Lscale/Lj) for Lj in Ljb.nzval]

"""
    StructureRealJacobianPlan

Everything the structure aware assembly reads, on a backend. There is no
segmented gather here and nothing proportional to the number of contributions:
the largest arrays are the sparsity structure and the precomputed linear term,
both of which are one entry per *stored entry* of the Jacobian rather than one
per contribution.

# Fields
- `colptr`, `rowval`: the stored structure, which is the Jacobian's
    transpose when `transposed`, so that a column of it is a row of the
    Jacobian.
- `lin`: the constant frequency dependent linear contribution of each
    stored entry.
- `linear`: the [`LinearTermGather`](@ref) which refreshes `lin`.
- `junctions`: the [`JunctionStructure`](@ref), the incidence triple
    product, the junction coefficients and the mode coupling index
    matrices, shared with every other plan of the system.
- `transposed`: whether the stored structure is the Jacobian's transpose, which
    it is on a device and is not for a matrix meant to be factorized directly.
- `linv`, `lptr`: the mode layout of the rows and the columns, the square
    Jacobian having one, and its inverse, which turn a stored entry back
    into a (node, mode) pair.
- `slots`: on a host, for a structure in the natural orientation, the
    decode of every real index of the layout ([`realslot`](@ref)), from
    which the assembly of a stored column reads the rows of its entries;
    `nothing` on a device and for a transposed structure, whose kernel
    decodes both sides of each entry.
- `assemble!`, `backend`: the compiled kernel, sized at plan time, and the
    KernelAbstractions backend it runs on.
- `n`: the number of stored entries, checked by
    [`assemblerealjacobian!`](@ref).
"""
struct StructureRealJacobianPlan{Ti<:Integer,T<:Real,VI,VT,LG,JS<:JunctionStructure{T},VS,K,B}
    colptr::VI
    rowval::VI
    lin::VT
    linear::LG
    junctions::JS
    transposed::Bool
    linv::VI
    lptr::VI
    slots::VS
    assemble!::K
    backend::B
    n::Int
end

"""
    realstructureentry(::Type{T}, rri, rci, Nmodes, Nfreq, ami, amc, pairptr,
        pairrow, pairjunc, paircoef, lmolj, linv, lptr, phimatrix)

The Josephson contribution to the stored entry of the real Jacobian at row
`rri` and column `rci`, both decoded here ([`realslot`](@ref)): what the
per entry kernel and the block preconditioner's `blockassemblykernel!`
compute for each entry. The per column assembly of a host decodes its stored
column once instead, and both sum the contribution with
[`realjosephsonentry`](@ref).
"""
@inline function realstructureentry(::Type{T}, rri, rci, Nmodes, Nfreq, ami,
        amc, pairptr, pairrow, pairjunc, paircoef, lmolj, linv, lptr,
        phimatrix) where {T}
    return realjosephsonentry(T, realslot(rri, linv, lptr, Nmodes),
        realslot(rci, linv, lptr, Nmodes), Nfreq, ami, amc, pairptr, pairrow,
        pairjunc, paircoef, lmolj, phimatrix)
end

"""
    realslot(r, linv, lptr, Nmodes)

The decode of the real index `r` of a layout with inverse `linv` and
pointer `lptr`: the offset of `r` within the real block of its complex
index, the width of that block, and the node and the mode of the complex
index.
"""
@inline function realslot(r, linv, lptr, Nmodes)
    @inbounds begin
        ci = Int(linv[r]); r0 = Int(lptr[ci])
        return (r - r0, Int(lptr[ci+1]) - r0, (ci - 1) ÷ Nmodes + 1,
            (ci - 1) % Nmodes + 1)
    end
end

# the decode of every real index of `layout` in `Int32`, which a host
# assembly reads rather than decoding each entry
realslottable(layout::ModeLayout, Nmodes::Integer) =
    [map(Int32, realslot(r, layout.inv, layout.ptr, Nmodes))
        for r in 1:layout.rdim]

"""
    realjosephsonentry(::Type{T}, rowslot, colslot, Nfreq, ami, amc, pairptr,
        pairrow, pairjunc, paircoef, lmolj, phimatrix)

The Josephson contribution to the stored entry of the real Jacobian whose row
and column decode to `rowslot` and `colslot` ([`realslot`](@ref)), which is
what every real assembly computes: the two kernels and the host loop here and
the block preconditioner's `blockassemblykernel!`.

The junctions incident on the node pair of the entry are looked up, and
their contributions summed. `(r0, c0)` and `(r0+1, c0+1)` are the real part
entries and `(r0+1, c0)` and `(r0, c0+1)` the imaginary part ones, so a
stored entry belongs to exactly one of the two and only one kind of
contribution can reach it.
"""
@inline function realjosephsonentry(::Type{T}, rowslot, colslot, Nfreq, ami,
        amc, pairptr, pairrow, pairjunc, paircoef, lmolj, phimatrix) where {T}
    @inbounds begin
        dr, wr, n1, m1 = map(Int, rowslot)
        dc, wc, n2, m2 = map(Int, colslot)

        acc = zero(T)
        ind = Int(ami[m1, m2]); indconj = Int(amc[m1, m2])
        if !(ind == 0 && indconj == 0)
            lo2 = Int(pairptr[n2]); hi2 = Int(pairptr[n2+1]) - 1
            isre = dr == dc
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
        return acc
    end
end

"""
    structureassemblykernel!(nzval, colptr, rowval, lin, phimatrix, ...)

Assemble one stored entry of the real Jacobian per work item, from the
circuit's structure.

The work item decodes the entry it owns into the node pair and mode pair it
belongs to, looks up the junctions incident on that node pair, and sums their
contributions: per incident junction the difference frequency coupling
`ami` and its conjugate partner `amc`, and the linear term last. Floating
point addition is not associative, so that order is part of the result.
"""
@kernel function structureassemblykernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(lin), @Const(phimatrix), @Const(pairptr), @Const(pairrow),
        @Const(pairjunc), @Const(paircoef), @Const(lmolj), @Const(ami),
        @Const(amc), @Const(linv), @Const(lptr), Nmodes, Nfreq, transposed)

    gid = @index(Global)
    T = eltype(nzval)
    @inbounds begin
        q = gid
        # the stored column this entry is in, which is a row of the Jacobian
        lo = storedcolumn(colptr, q)
        # the stored matrix is the transpose on a device, so its column is the
        # Jacobian's row; for a normally oriented one it is the other way
        rri = transposed ? lo : Int(rowval[q])
        rci = transposed ? Int(rowval[q]) : lo
        nzval[q] = realstructureentry(T, rri, rci, Nmodes, Nfreq, ami, amc,
            pairptr, pairrow, pairjunc, paircoef, lmolj, linv, lptr,
            phimatrix) + lin[q]
    end
end

"""
    structureassemblerowkernel!(nzval, colptr, rowval, lin, phimatrix, ...)

As [`structureassemblykernel!`](@ref), one work item per stored column
rather than per stored entry, which is how a host assembles a structure in
the natural orientation ([`realjacobiancolumnitem!`](@ref)).

The two differ in what they amortize against what they lose. Per entry, the
column has to be found by a binary search and both sides of the entry are
decoded; per column, the search is paid once, the column is decoded once,
and the rows are read from a table, but the entries a work item writes are
contiguous rather than interleaved with its neighbours', which costs
coalescing on a device and nothing on a host.
"""
@kernel function structureassemblerowkernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(lin), @Const(phimatrix), @Const(pairptr), @Const(pairrow),
        @Const(pairjunc), @Const(paircoef), @Const(lmolj), @Const(ami),
        @Const(amc), @Const(slots), Nfreq)
    gid = @index(Global)
    realjacobiancolumnitem!(nzval, gid, colptr, rowval, lin,
        phimatrix, pairptr, pairrow, pairjunc, paircoef, lmolj, ami, amc,
        slots, Nfreq)
end

"""
    realjacobiancolumnitem!(nzval, rr, colptr, rowval, lin, phimatrix,
        pairptr, pairrow, pairjunc, paircoef, lmolj, ami, amc, slots, Nfreq)

Assemble the stored column `rr` of a host plan in the natural orientation,
which the row kernel and the host loop share. The column is decoded once,
and the row of each entry read from `slots`, the decode of every real index
(`realslottable`); the entry is then what [`realstructureentry`](@ref)
gives, plus the linear term.
"""
@inline function realjacobiancolumnitem!(nzval, rr, colptr, rowval, lin,
        phimatrix, pairptr, pairrow, pairjunc, paircoef, lmolj, ami, amc,
        slots, Nfreq)
    T = eltype(nzval)
    @inbounds begin
        column = slots[rr]
        for q in Int(colptr[rr]):Int(colptr[rr+1])-1
            nzval[q] = realjosephsonentry(T, slots[Int(rowval[q])], column,
                Nfreq, ami, amc, pairptr, pairrow, pairjunc, paircoef,
                lmolj, phimatrix) + lin[q]
        end
    end
    return nothing
end

"""
    planstructurerealjacobian(Jt, T::Type{<:Real}, junctions::JunctionStructure,
        invLnm, Gnm, Cnm, wmodesm, wmodes2m, layout::ModeLayout, backend;
        transposed = true)

Build a [`StructureRealJacobianPlan`](@ref) with values of type `T` for
the real Jacobian whose structure is `Jt`, its rows and columns both in the
real representation `layout`, stored transposed (as a device factorization
wants it) when `transposed = true` and in the natural orientation
otherwise. Nothing here is proportional to the number of contributions:
what is stored is [`junctionpairtable`](@ref), whose size is set by the
circuit rather than by the mode count, and the constant linear term, which
is gathered on `backend` at the entries it reaches
([`LinearTermGather`](@ref)) and is zero at the others.
"""
function planstructurerealjacobian(Jt, ::Type{T},
    junctions::JunctionStructure{T}, invLnm::SparseMatrixCSC,
    Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, wmodesm::Diagonal,
    wmodes2m::Diagonal, layout::ModeLayout, backend;
    transposed::Bool = true) where {T<:Real}

    n = nnz(Jt)
    Ti = n < typemax(Int32) ? Int32 : Int
    d = x -> tobackend(backend, convert(Vector{Ti}, x))

    # the structure, adopted when it is already on the backend
    dcolptr, drowval = if Jt isa DeviceSparsePattern
        Jt.colptr, Jt.rowval
    else
        d(SparseArrays.getcolptr(Jt)), d(rowvals(Jt))
    end

    # the constant linear term, gathered on the backend
    dlinv = d(collect(layout.inv)); dlptr = d(collect(layout.ptr))
    linear = lineargather(Ti, dcolptr, drowval, invLnm, Gnm, Cnm, wmodesm,
        wmodes2m, layout, dlptr, transposed, backend)
    lin = gatherlinear!(fill!(KernelAbstractions.allocate(backend, T, n),
        zero(T)), linear, backend)

    # one work item per entry on a device, where coalescing pays for the
    # repeated decode, and one per stored column on a host, where it does
    # not and the decode of the rows is read from a table; a transposed
    # structure, the device's orientation, takes the per entry kernel
    # wherever it is
    slots = backend isa CPU && !transposed ?
        realslottable(layout, junctions.nmodes) : nothing
    assemble! = isnothing(slots) ? structureassemblykernel!(backend, 64) :
        structureassemblerowkernel!(backend, 64)
    return StructureRealJacobianPlan{Ti,T,typeof(dcolptr),typeof(lin),
        typeof(linear),typeof(junctions),typeof(slots),typeof(assemble!),
        typeof(backend)}(dcolptr, drowval, lin, linear, junctions, transposed,
        dlinv, dlptr, slots, assemble!, backend, n)
end

"""
    assemblerealjacobian!(nzval::AbstractVector,
        plan::StructureRealJacobianPlan, phimatrix::AbstractArray;
        synchronize = true)

Assemble the stored values of the real Jacobian into `nzval` from the Fourier
coefficients of `cos(phi(t))`, using the circuit's structure rather than a
precomputed gather. On a device the assembly kernel is synchronized before
returning unless `synchronize = false`, for a caller which orders the work
on the stream itself, as the transient step does.
"""
function assemblerealjacobian!(nzval::AbstractVector,
    plan::StructureRealJacobianPlan, phimatrix::AbstractArray; synchronize::Bool = true)
    length(nzval) == plan.n || throw(DimensionMismatch(
        lazy"`nzval` has length $(length(nzval)) but the plan assembles $(plan.n) entries."))
    js = plan.junctions
    colptr, rowval, lin, slots = plan.colptr, plan.rowval, plan.lin, plan.slots
    if isnothing(slots)
        # a device or a transposed structure: one work item per entry, which
        # decodes both sides
        plan.assemble!(nzval, colptr, rowval, lin, phimatrix, js.pairptr,
            js.pairrow, js.pairjunc, js.paircoef, js.lmolj, js.ami, js.amc,
            plan.linv, plan.lptr, js.nmodes, js.nfreq, plan.transposed;
            ndrange = plan.n)
    elseif hostloop(plan.backend, plan.n)
        # the per column assembly of the row kernel as a plain loop
        for rr in 1:length(colptr)-1
            realjacobiancolumnitem!(nzval, rr, colptr, rowval, lin,
                phimatrix, js.pairptr, js.pairrow, js.pairjunc, js.paircoef,
                js.lmolj, js.ami, js.amc, slots, js.nfreq)
        end
        return nzval
    else
        plan.assemble!(nzval, colptr, rowval, lin, phimatrix, js.pairptr,
            js.pairrow, js.pairjunc, js.paircoef, js.lmolj, js.ami, js.amc,
            slots, js.nfreq; ndrange = length(colptr) - 1)
    end
    # a caller assembling many matrices in a row synchronizes once after
    synchronize && KernelAbstractions.synchronize(plan.backend)
    return nzval
end

# Whether `n` work items on `backend` run as a plain loop rather than as a
# kernel. Launching a kernel on the CPU backend costs tens of microseconds of
# task scheduling, more than a small map or assembly costs, and these run once
# per Krylov product or time step; a large launch with threads to use keeps the
# threaded kernel.
hostloop(backend, n) = backend isa CPU &&
    (Threads.nthreads() == 1 || n <= hostlooplimit)

# the item count below which a CPU launch runs as a plain loop even with
# threads available: the launch costs more than the items
const hostlooplimit = 1 << 16

# the stored column of entry `q` of a compressed structure, the last `j`
# with `colptr[j] <= q`, by binary search: a stored column index per entry
# would cost memory and save no time
@inline function storedcolumn(colptr, q)
    lo = 1; hi = length(colptr) - 1
    @inbounds while lo < hi
        mid = (lo + hi + 1) >>> 1
        if colptr[mid] <= q
            lo = mid
        else
            hi = mid - 1
        end
    end
    return lo
end

# the stored entry at row `i` of column `j` of a compressed structure, or
# zero where there is none, by binary search in the column
@inline function storedposition(colptr, rowval, i, j)
    @inbounds begin
        lo = Int(colptr[j]); hi = Int(colptr[j+1]) - 1
        lo > hi && return 0
        while lo < hi
            mid = (lo + hi) >>> 1
            if Int(rowval[mid]) < i
                lo = mid + 1
            else
                hi = mid
            end
        end
        return Int(rowval[lo]) == i ? lo : 0
    end
end

# the value of a sparse matrix at (i, j), or zero. The linear term matrices
# are small and their columns short, so this is a handful of cached reads.
@inline function sparselookup(colptr, rowval, nzval, i, j)
    q = storedposition(colptr, rowval, i, j)
    return iszero(q) ? zero(eltype(nzval)) : @inbounds(nzval[q])
end

"""
    linearentry(lcolptr, lrowval, lnzval, gcolptr, growval, gnzval, wm,
        ccolptr, crowval, cnzval, wm2, ci, cj, dr, dc)
    linearentry(lcolptr, lrowval, lnzval, gcolptr, growval, gnzval, wm,
        ccolptr, crowval, cnzval, wm2, ci, cj)

The constant linear term `invLnm + im*Gnm*wmodesm - Cnm*wmodes2m` at the
complex position `(ci, cj)`: the entry at offset `(dr, dc)` of its real block
([`realblockterm`](@ref)), or without an offset the complex value. The
three matrices, compressed by columns, are looked up at the position, the
last two scaled by the mode frequency diagonals `wm` and `wm2` at the
column, and summed in the order `invLnm`, `Gnm`, `Cnm`: floating point
addition is not associative, so the order is part of the result.
"""
@inline function linearentry(lcolptr, lrowval, lnzval, gcolptr, growval,
        gnzval, wm, ccolptr, crowval, cnzval, wm2, ci, cj, dr, dc)
    @inbounds begin
        acc = realblockterm(sparselookup(lcolptr, lrowval, lnzval, ci, cj),
            dr, dc)
        acc += realblockterm((im * wm[cj]) *
            sparselookup(gcolptr, growval, gnzval, ci, cj), dr, dc)
        acc += realblockterm((-1 * wm2[cj]) *
            sparselookup(ccolptr, crowval, cnzval, ci, cj), dr, dc)
        return acc
    end
end

@inline function linearentry(lcolptr, lrowval, lnzval, gcolptr, growval,
        gnzval, wm, ccolptr, crowval, cnzval, wm2, ci, cj)
    @inbounds begin
        acc = sparselookup(lcolptr, lrowval, lnzval, ci, cj)
        acc += (im * wm[cj]) * sparselookup(gcolptr, growval, gnzval, ci, cj)
        acc += (-1 * wm2[cj]) * sparselookup(ccolptr, crowval, cnzval, ci, cj)
        return acc
    end
end

"""
    LinearTermGather{VI,NT}

The constant linear term `invLnm + im*Gnm*wmodesm - Cnm*wmodes2m` of an
assembly plan: the stored entries of the plan it reaches, and the three
matrices it is read from, on the plan's backend.

A stored entry takes something from the linear term only if its complex
position is stored in one of the three matrices, which the structure fixes
whatever the values, and most entries of a Jacobian couple modes through the
junctions alone. So the entries the linear term reaches are found once, when
the plan is built, and [`refreshlinear!`](@ref) rewrites those and nothing
else for new component values: every other entry takes zero, whatever the
values.

# Fields
- `pos`: the stored entry of the plan each reached entry is, or zero where
    the plan's structure does not store it.
- `row`, `col`: the complex position it belongs to.
- `part`: for a real plan, its offset `(dr, dc)` in the real block of that
    position, as `2dr + dc`; empty for a complex plan.
- `inputs`: the three matrices compressed by columns, each with the mode
    frequency diagonal it multiplies, named and ordered as
    [`linearentry`](@ref) takes them. The values are the plan's own copies,
    which a refresh overwrites.
"""
struct LinearTermGather{VI,NT}
    pos::VI
    row::VI
    col::VI
    part::VI
    inputs::NT
end

"""
    lineargather(::Type{Ti}, colptr, rowval, invLnm, Gnm, Cnm, wmodesm,
        wmodes2m, layout, lptr, transposed::Bool, backend)

The [`LinearTermGather`](@ref) of a plan whose structure is `colptr`,
`rowval` on `backend`, the Jacobian's transpose when `transposed`. Every
position stored in one of the three matrices is expanded to the entries of
its real block through `layout`, whose pointer `lptr` is on the backend, or
kept whole for a complex plan (`layout = lptr = nothing`), and located in
the structure. The walk is over the entries of the three matrices, so it
costs their size and not the plan's.
"""
function lineargather(::Type{Ti}, colptr, rowval, invLnm::SparseMatrixCSC,
        Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, wmodesm::Diagonal,
        wmodes2m::Diagonal, layout, lptr, transposed::Bool,
        backend) where {Ti}
    width(i) = isnothing(layout) ? 1 : Int(layout.ptr[i+1] - layout.ptr[i])
    # the entries, counted and then listed by column, block column, row and
    # block row, so that on a host a refresh writes the entries of an
    # untransposed structure in the order they are stored
    rows = Int[]
    n = 0
    for cj in axes(invLnm, 2)
        unionrows!(rows, invLnm, Gnm, Cnm, cj)
        n += width(cj) * sum(width, rows; init = 0)
    end
    row = Vector{Ti}(undef, n); col = Vector{Ti}(undef, n)
    part = Vector{Ti}(undef, isnothing(layout) ? 0 : n)
    k = 0
    for cj in axes(invLnm, 2)
        unionrows!(rows, invLnm, Gnm, Cnm, cj)
        for dc in 0:width(cj)-1, ci in rows, dr in 0:width(ci)-1
            k += 1
            row[k] = ci; col[k] = cj
            isnothing(layout) || (part[k] = 2dr + dc)
        end
    end
    d = x -> tobackend(backend, x)
    drow, dcol, dpart = d(row), d(col), d(part)
    pos = KernelAbstractions.allocate(backend, Ti, n)
    if hostloop(backend, n)
        for k in 1:n
            storedpositionitem!(pos, colptr, rowval, drow, dcol, dpart, lptr,
                transposed, k)
        end
    elseif n > 0
        storedpositionkernel!(backend, 64)(pos, colptr, rowval, drow, dcol,
            dpart, lptr, transposed; ndrange = n)
        KernelAbstractions.synchronize(backend)
    end
    # the values are the gather's own copies, which a refresh overwrites,
    # held as complex and the frequencies as real doubles whatever the
    # matrices hold, so that a gather has one type whichever of its
    # matrices are real or complex
    di = x -> tobackend(backend, convert(Vector{Ti}, x))
    dv = x -> tobackend(backend, Vector{ComplexF64}(x))
    df = x -> tobackend(backend, Vector{Float64}(x))
    inputs = (lcolptr = di(SparseArrays.getcolptr(invLnm)),
        lrowval = di(rowvals(invLnm)), lnzval = dv(nonzeros(invLnm)),
        gcolptr = di(SparseArrays.getcolptr(Gnm)), growval = di(rowvals(Gnm)),
        gnzval = dv(nonzeros(Gnm)), wm = df(wmodesm.diag),
        ccolptr = di(SparseArrays.getcolptr(Cnm)), crowval = di(rowvals(Cnm)),
        cnzval = dv(nonzeros(Cnm)), wm2 = df(wmodes2m.diag))
    return LinearTermGather{typeof(pos),typeof(inputs)}(pos, drow, dcol,
        dpart, inputs)
end

# the rows stored in column `j` of any of three sparse matrices, ascending
# and each once, into `rows`
function unionrows!(rows, A, B, C, j)
    empty!(rows)
    for M in (A, B, C)
        append!(rows, view(rowvals(M), nzrange(M, j)))
    end
    return unique!(sort!(rows))
end

# the stored entry which entry `k` of a gather is, or zero where the
# structure stores none: its real row and column from its complex position
# and block offset (the complex position itself for a complex plan, with
# `lptr === nothing`), searched for in the stored column, which is the
# Jacobian's row when `transposed`
@inline function storedpositionitem!(pos, colptr, rowval, row, col, part,
        lptr, transposed, k)
    @inbounds begin
        ci = Int(row[k]); cj = Int(col[k])
        if isnothing(lptr)
            r, c = ci, cj
        else
            p = Int(part[k])
            r = Int(lptr[ci]) + (p >> 1); c = Int(lptr[cj]) + (p & 1)
        end
        pos[k] = transposed ? storedposition(colptr, rowval, c, r) :
            storedposition(colptr, rowval, r, c)
    end
    return nothing
end

@kernel function storedpositionkernel!(pos, @Const(colptr), @Const(rowval),
        @Const(row), @Const(col), @Const(part), lptr, transposed)
    gid = @index(Global)
    storedpositionitem!(pos, colptr, rowval, row, col, part, lptr,
        transposed, gid)
end

"""
    refreshlinear!(lin, g::LinearTermGather, invLnm, Gnm, Cnm, wmodesm,
        wmodes2m, backend)

Write the linear term at the values of `invLnm`, `Gnm` and `Cnm` into the
entries of `lin` it reaches: the new values are copied into the gather's
own, and [`linearentry`](@ref) is evaluated at each reached entry, a plain
loop on the host and a kernel on a device. The matrices must have the
structure the gather was built on.
"""
function refreshlinear!(lin, g::LinearTermGather, invLnm, Gnm, Cnm,
        wmodesm::Diagonal, wmodes2m::Diagonal, backend)
    x = g.inputs
    (length(x.lnzval) == nnz(invLnm) && length(x.gnzval) == nnz(Gnm) &&
        length(x.cnzval) == nnz(Cnm)) || throw(DimensionMismatch(
        "the linear term matrices do not have the structure the plan was built on; build a new plan."))
    copyto!(x.lnzval, nonzeros(invLnm))
    copyto!(x.gnzval, nonzeros(Gnm))
    copyto!(x.cnzval, nonzeros(Cnm))
    copyto!(x.wm, wmodesm.diag)
    copyto!(x.wm2, wmodes2m.diag)
    return gatherlinear!(lin, g, backend)
end

# the evaluation of `refreshlinear!`, on the values the gather holds
function gatherlinear!(lin, g::LinearTermGather, backend)
    x = g.inputs
    args = (g.pos, g.row, g.col, g.part, x.lcolptr, x.lrowval, x.lnzval,
        x.gcolptr, x.growval, x.gnzval, x.wm, x.ccolptr, x.crowval, x.cnzval,
        x.wm2)
    n = length(g.pos)
    if hostloop(backend, n)
        for k in 1:n
            linearitem!(lin, k, args...)
        end
    elseif n > 0
        linearkernel!(backend, 64)(lin, args...; ndrange = n)
        KernelAbstractions.synchronize(backend)
    end
    return lin
end

# the reached entry `k` of a gather: the real block entry of its position
# for a real `lin`, the complex value for a complex one
@inline function linearitem!(lin, k, pos, row, col, part, lcolptr, lrowval,
        lnzval, gcolptr, growval, gnzval, wm, ccolptr, crowval, cnzval, wm2)
    @inbounds begin
        q = Int(pos[k])
        iszero(q) && return nothing
        ci = Int(row[k]); cj = Int(col[k])
        if eltype(lin) <: Complex
            lin[q] = linearentry(lcolptr, lrowval, lnzval, gcolptr, growval,
                gnzval, wm, ccolptr, crowval, cnzval, wm2, ci, cj)
        else
            p = Int(part[k])
            lin[q] = eltype(lin)(linearentry(lcolptr, lrowval, lnzval,
                gcolptr, growval, gnzval, wm, ccolptr, crowval, cnzval, wm2,
                ci, cj, p >> 1, p & 1))
        end
    end
    return nothing
end

@kernel function linearkernel!(lin, @Const(pos), @Const(row), @Const(col),
        @Const(part), @Const(lcolptr), @Const(lrowval), @Const(lnzval),
        @Const(gcolptr), @Const(growval), @Const(gnzval), @Const(wm),
        @Const(ccolptr), @Const(crowval), @Const(cnzval), @Const(wm2))
    gid = @index(Global)
    linearitem!(lin, gid, pos, row, col, part, lcolptr, lrowval,
        lnzval, gcolptr, growval, gnzval, wm, ccolptr, crowval, cnzval, wm2)
end

"""
    refreshvalues!(plan::StructureRealJacobianPlan, invLnm, Gnm, Cnm,
        wmodesm, wmodes2m, Ljb, Lscale)

Rewrite the value arrays of an assembly plan for new component values under
the same structure: the constant linear term at the entries it reaches
([`refreshlinear!`](@ref)), and `Lscale/Lj` per junction. The incidence
products and the structure are what they were.
"""
function refreshvalues!(plan::StructureRealJacobianPlan, invLnm, Gnm, Cnm,
        wmodesm, wmodes2m, Ljb::SparseVector, Lscale)
    refreshlinear!(plan.lin, plan.linear, invLnm, Gnm, Cnm, wmodesm,
        wmodes2m, plan.backend)
    refreshvalues!(plan.junctions, Ljb, Lscale)
    return plan
end

# ---------------------------------------------------------------------------
# the complex Jacobian, from the same structure
# ---------------------------------------------------------------------------
#
# The holomorphic Jacobian is the same incidence triple product without the
# real block expansion, so the same table serves it. It is in fact simpler: a
# stored entry names one mode pair, so `Amatrixindices` at that pair is one
# number, and every contribution to the entry is therefore either conjugated
# or not.

"""
    StructureComplexJosephsonPlan{Ti,VI,JS,K,B}

The Josephson contribution to the complex Jacobian, as a linear map from the
Fourier coefficients of `cos(phi(t))` to the stored entries of a matrix with a
given structure.

Each stored entry is computed from the circuit's structure rather than
gathered from one entry per contribution. It is used both to assemble the
Jacobian and on its own, as the map applied to other coefficient arrays by
the linearized solve, which is why it is a plan for the Josephson term
rather than for the whole Jacobian.

# Fields
- `colptr`, `rowval`: the stored structure, the transpose when `transposed`.
- `junctions`: the [`JunctionStructure`](@ref), of which the pair table,
    the junction coefficients and `ami` are read.
- `transposed`: as in [`StructureRealJacobianPlan`](@ref).
- `assemble!`, `backend`: the compiled kernel, one work item per stored
    column for a host structure in the natural orientation and one per
    stored entry otherwise, and its backend.
- `n`: the number of stored entries.
"""
struct StructureComplexJosephsonPlan{Ti<:Integer,VI,JS<:JunctionStructure,K,B}
    colptr::VI
    rowval::VI
    junctions::JS
    transposed::Bool
    assemble!::K
    backend::B
    n::Int
end

# the Josephson contribution to the stored entry of the complex Jacobian at
# row `ci` and column `cj`, the complex counterpart of `realstructureentry`:
# one lookup of the coupling per incident junction, no conjugate partner,
# the Fourier coefficient conjugated where the coupling index is negative
@inline function josephsonentry(::Type{T}, ci, cj, Nmodes, Nfreq, ami,
        pairptr, pairrow, pairjunc, paircoef, lmolj, phimatrix) where {T}
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
                acc += (paircoef[k] * lmolj[b]) * (ind < 0 ? conj(v) : v)
            end
        end
        return acc
    end
end

# one work item per stored entry of the complex Jacobian, the counterpart of
# `structureassemblykernel!`
@kernel function complexjosephsonkernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(phimatrix), @Const(pairptr), @Const(pairrow), @Const(pairjunc),
        @Const(paircoef), @Const(lmolj), @Const(ami), Nmodes, Nfreq,
        transposed)
    gid = @index(Global)
    T = eltype(nzval)
    @inbounds begin
        q = gid
        lo = storedcolumn(colptr, q)
        # a compressed column of the stored structure is a column of the
        # Jacobian, or a row of it when the transpose is what is stored
        ri = transposed ? lo : Int(rowval[q])
        ci = transposed ? Int(rowval[q]) : lo
        nzval[q] = josephsonentry(T, ri, ci, Nmodes, Nfreq, ami,
            pairptr, pairrow, pairjunc, paircoef, lmolj, phimatrix)
    end
end

# one work item per stored column of the complex Jacobian, the counterpart
# of `structureassemblerowkernel!`
@kernel function complexjosephsonrowkernel!(nzval, @Const(colptr), @Const(rowval),
        @Const(phimatrix), @Const(pairptr), @Const(pairrow), @Const(pairjunc),
        @Const(paircoef), @Const(lmolj), @Const(ami), Nmodes, Nfreq)
    gid = @index(Global)
    complexjosephsoncolumnitem!(nzval, gid, colptr, rowval,
        phimatrix, pairptr, pairrow, pairjunc, paircoef, lmolj, ami, Nmodes,
        Nfreq)
end

# the stored column `gid` of a host plan in the natural orientation, which
# the row kernel and the host loop share: a column of the stored structure
# is a column of the Jacobian
@inline function complexjosephsoncolumnitem!(nzval, gid, colptr, rowval,
        phimatrix, pairptr, pairrow, pairjunc, paircoef, lmolj, ami, Nmodes,
        Nfreq)
    T = eltype(nzval)
    @inbounds for q in Int(colptr[gid]):Int(colptr[gid+1])-1
        nzval[q] = josephsonentry(T, Int(rowval[q]), gid, Nmodes, Nfreq, ami,
            pairptr, pairrow, pairjunc, paircoef, lmolj, phimatrix)
    end
    return nothing
end

"""
    planstructurecomplexjosephson(Jx::SparseMatrixCSC,
        junctions::JunctionStructure, backend; transposed = false)

Build a [`StructureComplexJosephsonPlan`](@ref) for the structure of `Jx`.
"""
function planstructurecomplexjosephson(Jx::SparseMatrixCSC,
    junctions::JunctionStructure, backend; transposed::Bool = false)

    n = nnz(Jx)
    Ti = n < typemax(Int32) ? Int32 : Int
    d = x -> tobackend(backend, convert(Vector{Ti}, x))
    # one work item per stored column on a host, and per stored entry on a
    # device and for a transposed structure, as the real plan assembles
    assemble! = backend isa CPU && !transposed ?
        complexjosephsonrowkernel!(backend, 64) :
        complexjosephsonkernel!(backend, 64)
    colptr = d(SparseArrays.getcolptr(Jx))
    return StructureComplexJosephsonPlan{Ti,typeof(colptr),typeof(junctions),
        typeof(assemble!),typeof(backend)}(colptr, d(rowvals(Jx)), junctions,
        transposed, assemble!, backend, n)
end

"""
    addjosephsonterm!(nzval::AbstractVector,
        plan::StructureComplexJosephsonPlan, phimatrix)

Write the Josephson contribution into `nzval`, overwriting it.
"""
function addjosephsonterm!(nzval::AbstractVector,
    plan::StructureComplexJosephsonPlan, phimatrix)
    length(nzval) == plan.n || throw(DimensionMismatch(
        lazy"`nzval` has length $(length(nzval)) but the plan assembles $(plan.n) entries."))
    js = plan.junctions
    args = (plan.colptr, plan.rowval, phimatrix, js.pairptr, js.pairrow,
        js.pairjunc, js.paircoef, js.lmolj, js.ami, js.nmodes, js.nfreq)
    if !(plan.backend isa CPU) || plan.transposed
        # a device or a transposed structure: one work item per entry
        plan.assemble!(nzval, args..., plan.transposed; ndrange = plan.n)
    elseif hostloop(plan.backend, plan.n)
        # the per column assembly of the row kernel as a plain loop
        for gid in 1:length(plan.colptr)-1
            complexjosephsoncolumnitem!(nzval, gid, args...)
        end
        return nzval
    else
        plan.assemble!(nzval, args...; ndrange = length(plan.colptr) - 1)
    end
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
    js = plan.junctions
    Nmodes, Nfreq = js.nmodes, js.nfreq
    @inbounds for cj in 1:(length(colptr)-1)
        n2 = (cj - 1) ÷ Nmodes + 1
        m2 = (cj - 1) % Nmodes + 1
        for q in Int(colptr[cj]):Int(colptr[cj+1])-1
            ci = Int(rowval[q])
            n1 = (ci - 1) ÷ Nmodes + 1
            m1 = (ci - 1) % Nmodes + 1
            ind = Int(js.ami[m1, m2])
            ind == 0 && continue
            dst = ind > 0 ? P : Q
            wq = w[q]
            for k in Int(js.pairptr[n2]):Int(js.pairptr[n2+1])-1
                if Int(js.pairrow[k]) == n1
                    b = Int(js.pairjunc[k])
                    dst[abs(ind) + Nfreq * (b - 1)] +=
                        (js.paircoef[k] * js.lmolj[b]) * wq
                end
            end
        end
    end
    return P, Q
end

"""
    StructureComplexJacobianPlan{TJ,VT,LG}

The complex Jacobian on a backend: the Josephson map of
[`StructureComplexJosephsonPlan`](@ref) and the constant linear term of each
stored entry, gathered once where it reaches ([`LinearTermGather`](@ref)),
which is what [`assemblecomplexjacobian!`](@ref) adds to it.
"""
struct StructureComplexJacobianPlan{TJ,VT,LG}
    josephson::TJ
    lin::VT
    linear::LG
end

"""
    planstructurecomplexjacobian(Jx::SparseMatrixCSC,
        junctions::JunctionStructure, invLnm, Gnm, Cnm, wmodesm, wmodes2m,
        backend; transposed = false, josephson = nothing)

Build a [`StructureComplexJacobianPlan`](@ref) for the structure of `Jx`: the
Josephson map over `junctions`, built here unless one built for the same
structure, backend and orientation is handed in as `josephson` (the plan
[`plancomplexjacobian`](@ref) returns), and the constant linear term at the
mode frequencies `wmodesm`, gathered once on `backend`.
"""
function planstructurecomplexjacobian(Jx::SparseMatrixCSC,
    junctions::JunctionStructure{T}, invLnm::SparseMatrixCSC,
    Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, wmodesm::Diagonal,
    wmodes2m::Diagonal, backend; transposed::Bool = false,
    josephson = nothing) where {T<:Real}

    if isnothing(josephson)
        josephson = planstructurecomplexjosephson(Jx, junctions, backend;
            transposed = transposed)
    end
    linear = lineargather(eltype(josephson.colptr), josephson.colptr,
        josephson.rowval, invLnm, Gnm, Cnm, wmodesm, wmodes2m, nothing,
        nothing, josephson.transposed, backend)
    lin = gatherlinear!(fill!(KernelAbstractions.allocate(backend,
        Complex{T}, nnz(Jx)), zero(Complex{T})), linear, backend)
    return StructureComplexJacobianPlan{typeof(josephson),typeof(lin),
        typeof(linear)}(josephson, lin, linear)
end

"""
    refreshvalues!(plan::StructureComplexJacobianPlan, invLnm, Gnm, Cnm,
        wmodesm, wmodes2m, Ljb, Lscale)

Rewrite the gathered linear term and the junction coefficients of a complex
plan for new component values under the same structure.
"""
function refreshvalues!(plan::StructureComplexJacobianPlan, invLnm, Gnm, Cnm,
        wmodesm, wmodes2m, Ljb::SparseVector, Lscale)
    refreshlinear!(plan.lin, plan.linear, invLnm, Gnm, Cnm, wmodesm,
        wmodes2m, plan.josephson.backend)
    refreshvalues!(plan.josephson.junctions, Ljb, Lscale)
    return plan
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
