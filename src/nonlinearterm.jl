
"""
    NonlinearTermPlan{Ti<:Integer,T<:Real,VI,VT,VF,KF,KB,B}

Fully precomputed, device-generic plan for the two linear maps which surround
the pointwise time domain nonlinearity of the harmonic balance system:

- the map `A` from the unknowns to the frequency domain coefficients of the
  Josephson junction branch fluxes, which folds together the real to complex
  conversion of the equivalent real representation, the incidence matrix
  product `Rbnm*z`, the gather onto the Josephson branches, and the packing
  into the array the inverse real transform consumes
  ([`phivectortomatrix!`](@ref)), and
- the map `B` from the frequency domain coefficients back to the node vector,
  which folds together the unpacking ([`phimatrixtovector!`](@ref)), the
  `Lmean/Lj` branch scaling, the transposed incidence matrix product, the
  frequency dependent linear term `K`, and the complex to real conversion.

Both are expressed, in the complex representation and in the equivalent real
one, as KernelAbstractions kernels in which **one work item owns one output
slot and reads only precomputed index maps**, so neither contains a
scatter, an atomic, or a write conflict, and both run unchanged on a device
once the arrays are moved there. The maps are fixed at the sparsity structure
and mode frequencies of the system, so they are built once by
[`plannonlinearterm`](@ref). The topology of the maps does not depend on the
representation, so the two share the incidence entries and the Josephson
segment lists and differ only in how a source is addressed and in which form
of the linear term matrix they gather.

Only the pointwise time domain function between the two maps distinguishes the
entry points of the system, so one plan serves all of them:

| entry point                   | between the transforms          | linear term |
|:------------------------------|:--------------------------------|:------------|
| [`residual!`](@ref)           | `sin(phitd)`                    | `+ K*x - b` |
| [`jacobianvectorproduct!`](@ref) | `cos(phitd) .* dirtd`        | `+ K*v`     |
| [`hessianvectorproduct!`](@ref)  | `-sin(phitd).*dirtd.*dirtd2` | none        |

# Fields

The forward map stores, per slot of the frequency domain array, the at most
two entries of the corresponding row of the incidence matrix. A row of a
branch incidence matrix has at most two entries because a branch is an edge
and one of its two nodes may be the ground node, so no segment structure and
no inner loop is needed:

- `n1`, `n2`: the real slot holding the real part of the source node flux, or
  zero when the entry is absent.
- `cn1`, `cn2`: the same entries addressed as complex indices, for the complex
  representation.
- `s1`, `s2`: the corresponding incidence matrix entries.
- `flags`: [`FCONJ`](@ref) when the slot is a conjugate symmetry target and
  [`FWIDE`](@ref) when the source mode is not self conjugate, so its imaginary
  part occupies the next real slot.

The backward map stores, per complex output index, a segment of contributions
gathered from the frequency domain array, and the linear term in compressed
sparse row form:

- `bptr`, `bsrc`, `bcoef`: the segment of frequency domain slots of the
  Josephson branches incident to that node, with the incidence matrix entry
  times `Lmean/Lj` as the coefficient.
- `kptr`, `kidx`, `kcoef`: the real form of the linear term matrix `K`, stored
  transposed so applying it is a gather over the entries of one output row.
  Empty when the plan was built with `realbackward = false`, which
  [`hasrealbackward`](@ref) reports.
- `cptr`, `cidx`, `ccoef`: the same for the complex form of `K`.
- `lptr`, `lwide`: the first real slot of each complex index and whether that
  mode is not self conjugate.
"""
struct NonlinearTermPlan{Ti<:Integer,T<:Real,VI,VT,VC,VF,KF,KB,KFC,KBC,KRC,KCR,B}
    # The forward map, one entry per slot of the frequency domain array.
    # The incidence entries `s1`, `s2` and the flags are shared by the two
    # representations; only the addressing of the source differs, so the
    # real representation stores the real slot holding the real part
    # (`n1`, `n2`) and the complex one the complex index (`cn1`, `cn2`).
    n1::VI
    s1::VT
    n2::VI
    s2::VT
    cn1::VI
    cn2::VI
    flags::VF
    # backward map, the Josephson contribution
    bptr::VI
    bsrc::VI
    bcoef::VT
    # the backward map of the linear term, in compressed sparse row form,
    # in the real representation and in the complex one
    kptr::VI
    kidx::VI
    kcoef::VT
    cptr::VI
    cidx::VI
    ccoef::VC
    # the real representation layout
    lptr::VI
    lwide::VF
    # row pointers whose ranges are all empty, which switch the linear term
    # off without a second kernel
    kptrzero::VI
    cptrzero::VI
    # cached kernel objects and the backend they were built for
    forward!::KF
    backward!::KB
    forwardcomplex!::KFC
    backwardcomplex!::KBC
    realtocomplex!::KRC
    complextoreal!::KCR
    backend::B
    # sizes, so applying the plan does not query the arrays
    nslots::Int
    ncomplex::Int
end

"""
    FCONJ

Flag bit marking a slot of the frequency domain array which is a conjugate
symmetry target, so the value gathered for it is conjugated. See
[`NonlinearTermPlan`](@ref).
"""
const FCONJ = Int32(1)

"""
    FWIDE

Flag bit marking a slot of the frequency domain array whose source mode is not
self conjugate, so the imaginary part of the source occupies the real slot
after the one addressed. See [`NonlinearTermPlan`](@ref).
"""
const FWIDE = Int32(2)

"""
    forwardtermkernel!(phimatrix, n1, s1, n2, s2, flags, xr)

Evaluate the map from the equivalent real representation of the unknowns to
the frequency domain coefficients of the Josephson junction branch fluxes, one
work item per slot of `phimatrix`. Each work item gathers the at most two node
fluxes of its branch through the incidence matrix entries held in the plan,
assembles the complex coefficient from the real slots of the real
representation, and conjugates it when the slot is a conjugate symmetry
target. Every slot is written exactly once, so the destination needs no
zeroing and the stores are contiguous. See [`NonlinearTermPlan`](@ref).
"""
@kernel function forwardtermkernel!(phimatrix, @Const(n1), @Const(s1),
        @Const(n2), @Const(s2), @Const(flags), @Const(xr))
    q = @index(Global)
    T = eltype(xr)
    @inbounds begin
        f = flags[q]
        wide = (f & FWIDE) != Int32(0)
        p1 = n1[q]
        p2 = n2[q]
        re = zero(T)
        im_ = zero(T)
        if !iszero(p1)
            c = s1[q]
            re += c*xr[p1]
            wide && (im_ += c*xr[p1+1])
        end
        if !iszero(p2)
            c = s2[q]
            re += c*xr[p2]
            wide && (im_ += c*xr[p2+1])
        end
        phimatrix[q] = (f & FCONJ) != Int32(0) ? Complex(re, -im_) :
            Complex(re, im_)
    end
end

"""
    forwardtermkernelcomplex!(phimatrix, cn1, s1, cn2, s2, flags, xc)

The complex representation counterpart of [`forwardtermkernel!`](@ref), which
gathers the same at most two node fluxes through the same incidence entries
but addresses them as complex indices, so no real slot arithmetic is needed
and the [`FWIDE`](@ref) flag is not consulted. See
[`NonlinearTermPlan`](@ref).
"""
@kernel function forwardtermkernelcomplex!(phimatrix, @Const(cn1), @Const(s1),
        @Const(cn2), @Const(s2), @Const(flags), @Const(xc))
    q = @index(Global)
    @inbounds begin
        p1 = cn1[q]
        p2 = cn2[q]
        v = zero(eltype(phimatrix))
        if !iszero(p1)
            v += s1[q]*xc[p1]
        end
        if !iszero(p2)
            v += s2[q]*xc[p2]
        end
        phimatrix[q] = (flags[q] & FCONJ) != Int32(0) ? conj(v) : v
    end
end

"""
    backwardtermkernelcomplex!(out, bptr, bsrc, bcoef, phimatrix, kptr, kidx,
        kcoef, xc)

The complex representation counterpart of [`backwardtermkernel!`](@ref). Each
work item owns one complex output entry, gathers the same segment of
Josephson contributions and one row of the complex form of the linear term
matrix, so the real slot bookkeeping of the real representation is not
needed. Passing the `cptrzero` field of the plan as `kptr` switches the linear
term off. See [`NonlinearTermPlan`](@ref).
"""
@kernel function backwardtermkernelcomplex!(out, @Const(bptr), @Const(bsrc),
        @Const(bcoef), @Const(phimatrix), @Const(kptr), @Const(kidx),
        @Const(kcoef), @Const(xc))
    k = @index(Global)
    @inbounds begin
        acc = zero(eltype(out))
        for t in bptr[k]:bptr[k+1]-1
            acc += bcoef[t]*phimatrix[bsrc[t]]
        end
        for t in kptr[k]:kptr[k+1]-1
            acc += kcoef[t]*xc[kidx[t]]
        end
        out[k] = acc
    end
end

"""
    backwardtermkernel!(out, bptr, bsrc, bcoef, phimatrix, kptr, kidx, kcoef,
        xr, lptr, lwide)

Evaluate the map from the frequency domain coefficients back to the node
vector in the equivalent real representation, one work item per complex output
index. Each work item gathers the coefficients of the Josephson branches
incident to its node, scaled by the incidence matrix entry and `Lmean/Lj`,
then gathers the one or two rows of the real form of the linear term matrix
applied to `xr`, and writes the one or two real slots it owns. Those slots are
disjoint across work items. Passing the `kptrzero` field of the plan as `kptr`
switches the linear term off. See [`NonlinearTermPlan`](@ref).
"""
@kernel function backwardtermkernel!(out, @Const(bptr), @Const(bsrc),
        @Const(bcoef), @Const(phimatrix), @Const(kptr), @Const(kidx),
        @Const(kcoef), @Const(xr), @Const(lptr), @Const(lwide))
    k = @index(Global)
    T = eltype(out)
    @inbounds begin
        p = lptr[k]
        accre = zero(T)
        accim = zero(T)
        for t in bptr[k]:bptr[k+1]-1
            v = phimatrix[bsrc[t]]
            c = bcoef[t]
            accre += c*real(v)
            accim += c*imag(v)
        end
        for t in kptr[p]:kptr[p+1]-1
            accre += kcoef[t]*xr[kidx[t]]
        end
        out[p] = accre
        if lwide[k] == Int32(1)
            for t in kptr[p+1]:kptr[p+2]-1
                accim += kcoef[t]*xr[kidx[t]]
            end
            out[p+1] = accim
        end
    end
end

"""
    realtocomplexkernel!(xc, xr, lptr, lwide)

Expand the equivalent real representation of a node vector into the complex
one, one work item per complex entry, reading the real slot layout the plan
already carries. The imaginary part of a self conjugate mode is written as an
explicit zero, so the result does not depend on what was in `xc` beforehand.
Each work item owns its own entry, so this is conflict free. The scalar
[`real_to_complex!`](@ref) it replaces on this path advances a serial cursor
and reads a `BitVector`, neither of which a device can do. See
[`NonlinearTermPlan`](@ref).
"""
@kernel function realtocomplexkernel!(xc, @Const(xr), @Const(lptr),
        @Const(lwide))
    k = @index(Global)
    T = eltype(xr)
    @inbounds begin
        p = lptr[k]
        xc[k] = lwide[k] == Int32(1) ? Complex(xr[p], xr[p+1]) :
            Complex(xr[p], zero(T))
    end
end

"""
    complextorealkernel!(xr, xc, lptr, lwide)

Contract the complex representation of a node vector into the equivalent real
one, one work item per complex entry. The imaginary part of a self conjugate
mode is never read. The one or two real slots each work item writes are
disjoint from every other's, so this is conflict free. The inverse of
[`realtocomplexkernel!`](@ref). See [`NonlinearTermPlan`](@ref).
"""
@kernel function complextorealkernel!(xr, @Const(xc), @Const(lptr),
        @Const(lwide))
    k = @index(Global)
    @inbounds begin
        p = lptr[k]
        v = xc[k]
        xr[p] = real(v)
        if lwide[k] == Int32(1)
            xr[p+1] = imag(v)
        end
    end
end

"""
    plannonlinearterm(Rbnm::SparseMatrixCSC, Ljb::SparseVector, Lmean,
        Nbranches::Integer, freqindexmap::Vector{Int},
        conjsourceindices::Vector{Int}, conjtargetindices::Vector{Int},
        phimatrix::AbstractArray, Knm::SparseMatrixCSC, layout::ModeLayout,
        backend = CPU(); realbackward = true)

Build the [`NonlinearTermPlan`](@ref) for the harmonic balance system
described by the incidence matrix `Rbnm`, the Josephson inductances `Ljb`,
the frequency domain packing maps, the frequency domain array `phimatrix`
which fixes the slot layout, the collapsed linear term matrix `Knm` (see
[`linearterm`](@ref)) and the real representation `layout`.

The per branch node and sign lists of the forward map come from
[`branchnodesandsigns`](@ref), the same function the assembled Jacobian plans
use, which also verifies that `Rbnm` has the expected mode-diagonal
`diagrepeat` structure. The backward map reads the columns of `Rbnm`
directly, because the entries of column `k` are exactly the branch rows
contributing to output `k`, which is what makes the transposed product a
gather.

Throws an `ArgumentError` if a branch touches more than two nodes, which the
forward map's flat two entry form assumes.

`realbackward = false` leaves out the real representation of the linear
term, the largest allocation of the plan, which only the backward map of
the real representation reads; a solve which stays in the complex
representation (`method = :quasinewton`) never applies it.
"""
function plannonlinearterm(Rbnm::SparseMatrixCSC, Ljb::SparseVector, Lmean,
    Nbranches::Integer, freqindexmap::Vector{Int},
    conjsourceindices::Vector{Int}, conjtargetindices::Vector{Int},
    phimatrix::AbstractArray, Knm::SparseMatrixCSC, layout::ModeLayout,
    backend = CPU(); realbackward::Bool = true)

    NJJ = length(Ljb.nzval)
    Nmodes = length(freqindexmap)
    Nmatrix = prod(size(phimatrix)[1:end-1])
    nc = layout.dim
    nslots = length(phimatrix)
    # 32 bit indices whenever the ranges permit, halving the memory traffic
    # of the plan, as `plancomplexjacobian` does
    Ti = (nslots <= typemax(Int32) && nc <= typemax(Int32) &&
        nnz(Knm) <= typemax(Int32) && size(Rbnm, 1) <= typemax(Int32)) ?
        Int32 : Int
    # the precision of the plan follows the frequency domain array it is
    # built around: a `Complex{Float32}` phimatrix gives `Float32` incidence
    # entries, branch coefficients and linear term
    T = real(eltype(phimatrix))

    # invert the frequency index map, so that a conjugate symmetry target
    # can be traced back to the entry of the branch flux vector which feeds
    # its source and the forward map is one kernel rather than a scatter
    # followed by a conjugating pass
    invfreqindexmap = Dict{Int,Int}(freqindexmap[j] => j for j in 1:Nmodes)
    slotmode = zeros(Int, nslots)
    slotflag = zeros(Int32, nslots)
    for i in 1:NJJ
        off = (i-1)*Nmatrix
        for j in 1:Nmodes
            slotmode[freqindexmap[j]+off] = j
        end
        for j in eachindex(conjtargetindices)
            haskey(invfreqindexmap, conjsourceindices[j]) || throw(
                ArgumentError(lazy"the conjugate source index $(conjsourceindices[j]) is not in the frequency index map."))
            slotmode[conjtargetindices[j]+off] = invfreqindexmap[conjsourceindices[j]]
            slotflag[conjtargetindices[j]+off] = FCONJ
        end
    end

    n1 = zeros(Ti, nslots)
    n2 = zeros(Ti, nslots)
    cn1 = zeros(Ti, nslots)
    cn2 = zeros(Ti, nslots)
    s1 = zeros(T, nslots)
    s2 = zeros(T, nslots)
    flags = zeros(Int32, nslots)
    # the incidence matrix is mode diagonal, so a branch has the same node
    # and sign list for every mode and the list is shared with the assembled
    # Jacobian plans
    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
    @inbounds for q in 1:nslots
        j = slotmode[q]
        iszero(j) && continue
        i = (q-1) ÷ Nmatrix + 1
        ns = nodesandsigns[Ljb.nzind[i]]
        length(ns) <= 2 || throw(ArgumentError(
            lazy"a branch touches $(length(ns)) nodes; the forward map assumes at most two."))
        f = slotflag[q]
        for (c, (node, sgn)) in enumerate(ns)
            k = (node-1)*Nmodes + j
            layout.w[k] || (f |= FWIDE)
            if isone(c)
                n1[q] = layout.ptr[k]
                cn1[q] = k
                s1[q] = T(sgn)
            else
                n2[q] = layout.ptr[k]
                cn2[q] = k
                s2[q] = T(sgn)
            end
        end
        flags[q] = f
    end

    # The backward map: for each complex output index, the frequency domain
    # slots of the Josephson branches incident on it. The entries of column
    # k of `Rbnm` are the branch rows which contribute to output k, so the
    # transposed product is a gather and needs no scatter or atomic.
    bptr, bsrc, bcoef = backwardjosephsonmap(Ti, T, Rbnm, Ljb, Lmean, Nmodes,
        Nmatrix, freqindexmap, nc)

    # The linear term in both representations, each transposed so that
    # applying it is a gather over the entries of one output row. Only the
    # real representation backward map reads the real form, so with
    # `realbackward = false` it is not built and its three arrays are left
    # empty; converting `K` to the real layout is the largest allocation of
    # the plan, and a solve in the complex representation never applies it.
    kptr, kidx, kcoef = if realbackward
        reallinearmap(Ti, T, Knm, layout)
    else
        (Ti[], Ti[], T[])
    end
    cptr, cidx, ccoef = complexlinearmap(Ti, T, Knm)

    lptr = convert(Vector{Ti}, layout.ptr[1:nc])
    lwide = Int32[layout.w[k] ? 0 : 1 for k in 1:nc]
    kptrzero = ones(Ti, length(kptr))
    cptrzero = ones(Ti, length(cptr))

    # the work sizes are baked into the kernel objects, since they are fixed
    # by the sparsity structure; partitioning the index space on every
    # launch would allocate
    groupsize = 64
    forward! = forwardtermkernel!(backend, groupsize, nslots)
    backward! = backwardtermkernel!(backend, groupsize, nc)
    forwardcomplex! = forwardtermkernelcomplex!(backend, groupsize, nslots)
    backwardcomplex! = backwardtermkernelcomplex!(backend, groupsize, nc)
    realtocomplex! = realtocomplexkernel!(backend, groupsize, nc)
    complextoreal! = complextorealkernel!(backend, groupsize, nc)

    # the index computation above is sequential, so the maps are built on
    # the host and then moved to the backend; on `CPU()` a dense array is
    # adopted without a copy (see `tobackend`)
    d = x -> tobackend(backend, x)
    dn1, ds1, dn2, ds2 = d(n1), d(s1), d(n2), d(s2)
    dcn1, dcn2, dflags = d(cn1), d(cn2), d(flags)
    dbptr, dbsrc, dbcoef = d(bptr), d(bsrc), d(bcoef)
    dkptr, dkidx, dkcoef = d(kptr), d(kidx), d(kcoef)
    dcptr, dcidx, dccoef = d(cptr), d(cidx), d(ccoef)
    dlptr, dlwide = d(lptr), d(lwide)
    dkptrzero, dcptrzero = d(kptrzero), d(cptrzero)

    return NonlinearTermPlan{Ti,T,typeof(dn1),typeof(ds1),typeof(dccoef),
        typeof(dflags), typeof(forward!), typeof(backward!),
        typeof(forwardcomplex!), typeof(backwardcomplex!),
        typeof(realtocomplex!), typeof(complextoreal!), typeof(backend)}(
        dn1, ds1, dn2, ds2, dcn1, dcn2, dflags, dbptr, dbsrc, dbcoef,
        dkptr, dkidx, dkcoef, dcptr, dcidx, dccoef, dlptr, dlwide,
        dkptrzero, dcptrzero, forward!, backward!, forwardcomplex!,
        backwardcomplex!, realtocomplex!, complextoreal!, backend,
        nslots, nc)
end

# The three parts of the backward map which carry values rather than
# structure: the Josephson coefficients `Lmean/Lj` times the incidence entry,
# and the linear term in each representation. They are built here for the
# plan and again by `refreshvalues!` when the values move under a fixed
# structure, so the two cannot disagree.
function backwardjosephsonmap(::Type{Ti}, ::Type{T}, Rbnm::SparseMatrixCSC,
        Ljb::SparseVector, Lmean, Nmodes::Integer, Nmatrix::Integer,
        freqindexmap::Vector{Int}, nc::Integer) where {Ti,T}
    NJJ = length(Ljb.nzval)
    jjofbranch = Dict{Int,Int}(Ljb.nzind[i] => i for i in 1:NJJ)
    bptr = Vector{Ti}(undef, nc+1)
    bptr[1] = 1
    bsrc = Ti[]
    bcoef = T[]
    rowsrbnm = rowvals(Rbnm)
    valsrbnm = nonzeros(Rbnm)
    @inbounds for k in 1:nc
        for t in nzrange(Rbnm, k)
            r = rowsrbnm[t]
            i = get(jjofbranch, (r-1) ÷ Nmodes + 1, 0)
            iszero(i) && continue
            j = (r-1) % Nmodes + 1
            push!(bsrc, freqindexmap[j] + (i-1)*Nmatrix)
            push!(bcoef, T(valsrbnm[t])*T(Lmean/Ljb.nzval[i]))
        end
        bptr[k+1] = length(bsrc)+1
    end
    return bptr, bsrc, bcoef
end

function reallinearmap(::Type{Ti}, ::Type{T}, Knm::SparseMatrixCSC,
        layout::ModeLayout) where {Ti,T}
    KrT = sparse(transpose(complex_to_real(Knm, layout, layout)))
    return (convert(Vector{Ti}, SparseArrays.getcolptr(KrT)),
        convert(Vector{Ti}, rowvals(KrT)),
        convert(Vector{T}, nonzeros(KrT)))
end

function complexlinearmap(::Type{Ti}, ::Type{T},
        Knm::SparseMatrixCSC) where {Ti,T}
    KnmT = sparse(transpose(Knm))
    return (convert(Vector{Ti}, SparseArrays.getcolptr(KnmT)),
        convert(Vector{Ti}, rowvals(KnmT)),
        convert(Vector{Complex{T}}, nonzeros(KnmT)))
end

"""
    ValueMaps

Where each value array of a [`NonlinearTermPlan`](@ref) and the linear term
it was built from read their numbers, under a fixed structure, so that new
component values can be written into them without rebuilding anything.

# Fields
- `ml`, `mg`, `mc`: for each stored entry of `Knm`, the stored entry of
  `invLnm`, `Gnm` and `Cnm` it takes, or zero.
- `kmap`, `kre`, `kim`: for each entry of the real form of the linear term,
  the entry of `Knm` it comes from and the coefficient of its real and its
  imaginary part, which are `0` or `±1`.
- `cmap`: for each entry of the complex form, the entry of `Knm`.
- `bjunc`, `bsgn`: for each Josephson coefficient, its junction and the
  incidence entry, so the coefficient is `bsgn*Lmean/Lj`.

Found by [`valuemaps`](@ref), which probes the conversions with tagged
values rather than restating their rules.
"""
struct ValueMaps{VI,VT,VC}
    ml::Vector{Int}
    mg::Vector{Int}
    mc::Vector{Int}
    kmap::VI
    kre::VT
    kim::VT
    cmap::VI
    bjunc::VI
    bsgn::VT
    knz::VC            # the entries of `Knm`, where the plan lives
end

# the entry a probe tagged, read back through a linear map: zero where the
# entry depends on nothing, else the tag with the sign the map gave it
_tagof(v) = (t = round(Int, abs(v)); (t, iszero(t) ? 0.0 : sign(v)))

"""
    valuemaps(plan::NonlinearTermPlan, Knm, invLnm, Gnm, Cnm, Rbnm, Ljb,
        layout, freqindexmap)

The [`ValueMaps`](@ref) of a plan, found by probing.

Each conversion is linear and each of its output entries reads one input
entry, so tagging the inputs with their own indices and reading the outputs
gives the map, and does so for whatever rule the conversion follows. A
conversion whose structure under the tags differs from the plan's -- one
which drops a zero, say -- has no fixed map and `nothing` is returned; the
refresh then rebuilds the arrays instead.
"""
function valuemaps(plan::NonlinearTermPlan{Ti,T}, Knm::SparseMatrixCSC,
        invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC,
        Rbnm::SparseMatrixCSC, Ljb::SparseVector, layout::ModeLayout,
        freqindexmap::Vector{Int}) where {Ti,T}
    # The three linear term matrices into `Knm`. Every entry of every
    # matrix carries a value, so the structure is the union the real values
    # produce, and the one whose map is read carries its tags on an axis of
    # its own: the inductance and the capacitance contribute on the real
    # axis and the conductance, through the `im`, on the imaginary one, so
    # the other two are given values which keep off the axis being read and
    # cannot cancel each other where they overlap, since the sum drops an
    # exact zero. Unit frequency diagonals, so a zero frequency does not
    # erase a tag.
    n = size(Knm, 1)
    ones_ = Diagonal(ones(n))
    tagged(A, v) = SparseMatrixCSC(size(A)..., copy(SparseArrays.getcolptr(A)),
        copy(rowvals(A)), v)
    tags(A) = Complex{Float64}.(1:nnz(A))
    fillc(A, v) = fill(Complex{Float64}(v), nnz(A))
    # The probe's structure is the union of the three, which holds every
    # entry of `Knm` and possibly more: a term which vanishes at the zero
    # frequency column is dropped from `Knm` by the sum and kept by the
    # probe. So the probe is read onto the structure of `Knm`, entry by
    # entry, and a missing one means there is no map.
    function probe(zl, zg, zc, part)
        K = linearterm(zl, zg, zc, ones_, ones_, Float64)
        size(K) == size(Knm) || return nothing
        m = zeros(Int, nnz(Knm))
        kcol, krow, knz = SparseArrays.getcolptr(K), rowvals(K), nonzeros(K)
        for j in axes(Knm, 2), t in nzrange(Knm, j)
            i = rowvals(Knm)[t]
            r = kcol[j]:(kcol[j+1] - 1)
            q = searchsortedfirst(view(krow, r), i)
            (q <= length(r) && krow[r[q]] == i) || return nothing
            m[t] = round(Int, part(knz[r[q]]))
        end
        return m
    end
    ml = probe(tagged(invLnm, tags(invLnm)), tagged(Gnm, fillc(Gnm, 1)),
        tagged(Cnm, fillc(Cnm, 2im)), real)
    mg = probe(tagged(invLnm, fillc(invLnm, 1)), tagged(Gnm, tags(Gnm)),
        tagged(Cnm, fillc(Cnm, 2)), imag)
    mc = probe(tagged(invLnm, fillc(invLnm, im)), tagged(Gnm, fillc(Gnm, 1)),
        tagged(Cnm, tags(Cnm)), v -> -real(v))
    (isnothing(ml) || isnothing(mg) || isnothing(mc)) && return nothing

    # `Knm` into its real form. Each real entry is one part of one entry
    # of `Knm` with a sign, so tagging the real parts with `k` and the
    # imaginary parts with `nz + k` tells the two apart by size.
    nz = nnz(Knm)
    kmap = zeros(Ti, length(plan.kidx)); kre = zeros(T, length(plan.kidx))
    kim = zeros(T, length(plan.kidx))
    if hasrealbackward(plan)
        both = tagged(Knm, Complex{Float64}[k + im*(nz + k) for k in 1:nz])
        p1, i1, v1 = reallinearmap(Ti, Float64, both, layout)
        (i1 == Array(plan.kidx) && p1 == Array(plan.kptr)) || return nothing
        for t in eachindex(i1)
            a, sa = _tagof(v1[t])
            iszero(a) && continue
            if a <= nz
                kmap[t] = a; kre[t] = sa
            else
                kmap[t] = a - nz; kim[t] = sa
            end
        end
    end
    # and into its complex form, a permutation
    pc, ic, vc = complexlinearmap(Ti, Float64, tagged(Knm, tags(Knm)))
    (ic == Array(plan.cidx) && pc == Array(plan.cptr)) || return nothing
    cmap = Ti[round(Ti, real(v)) for v in vc]

    # the Josephson coefficients: the junction and the incidence entry of
    # each, read from a map built with unit inductances
    nc = plan.ncomplex
    Nmodes = length(freqindexmap)
    Nmatrix = plan.nslots ÷ max(length(Ljb.nzval), 1)
    unit = SparseVector(length(Ljb), copy(Ljb.nzind), ones(Float64, length(Ljb.nzval)))
    _, bsrc, bsgn = backwardjosephsonmap(Ti, Float64, Rbnm, unit, 1.0, Nmodes,
        Nmatrix, freqindexmap, nc)
    Array(plan.bsrc) == bsrc || return nothing
    bjunc = Ti[(Int(src) - 1) ÷ Nmatrix + 1 for src in bsrc]

    d = x -> tobackend(plan.backend, x)
    return ValueMaps(ml, mg, mc, d(kmap), d(convert(Vector{T}, kre)),
        d(convert(Vector{T}, kim)), d(cmap), d(bjunc),
        d(convert(Vector{T}, bsgn)),
        KernelAbstractions.allocate(plan.backend, Complex{T}, nz))
end

# the three refreshes, one work item per output entry
@kernel function refreshrealkernel!(kcoef, @Const(knz), @Const(kmap),
        @Const(kre), @Const(kim))
    t = @index(Global)
    @inbounds begin
        v = knz[kmap[t]]
        kcoef[t] = kre[t]*real(v) + kim[t]*imag(v)
    end
end

@kernel function refreshcomplexkernel!(ccoef, @Const(knz), @Const(cmap))
    t = @index(Global)
    @inbounds ccoef[t] = knz[cmap[t]]
end

@kernel function refreshjosephsonkernel!(bcoef, @Const(ljnz), @Const(bjunc),
        @Const(bsgn), lmean)
    t = @index(Global)
    @inbounds bcoef[t] = bsgn[t]*(lmean/ljnz[bjunc[t]])
end

"""
    refreshvalues!(plan::NonlinearTermPlan, maps::ValueMaps, Knm, Ljb, Lmean)

Rewrite the value arrays of a plan through its [`ValueMaps`](@ref): three
kernels over the arrays where they live, and nothing allocated.
"""
function refreshvalues!(plan::NonlinearTermPlan{Ti,T}, maps::ValueMaps,
        Knm::SparseMatrixCSC, Ljb::SparseVector, Lmean) where {Ti,T}
    backend = plan.backend
    copyto!(maps.knz, nonzeros(Knm))
    ljnz = tobackend(backend, convert(Vector{T}, real.(Ljb.nzval)))
    if hasrealbackward(plan)
        refreshrealkernel!(backend, 64)(plan.kcoef, maps.knz, maps.kmap,
            maps.kre, maps.kim; ndrange = length(plan.kcoef))
    end
    refreshcomplexkernel!(backend, 64)(plan.ccoef, maps.knz, maps.cmap;
        ndrange = length(plan.ccoef))
    isempty(plan.bcoef) || refreshjosephsonkernel!(backend, 64)(plan.bcoef,
        ljnz, maps.bjunc, maps.bsgn, T(Lmean); ndrange = length(plan.bcoef))
    KernelAbstractions.synchronize(backend)
    return plan
end

"""
    refreshvalues!(plan::NonlinearTermPlan, Rbnm, Ljb, Lmean, Knm, layout,
        freqindexmap)

Rewrite the value arrays of a plan for new component values under the same
structure: the Josephson coefficients `Lmean/Lj`, and the linear term in
both representations. Everything else in the plan is structure, which the
new values must share; the maps are rebuilt and their structure checked
against the plan's before the values are copied in, on whichever backend
the plan lives on.

This is what makes a system reusable across the points of a sweep: the
transforms, the index maps and the kernels stay, and the numbers move.
"""
function refreshvalues!(plan::NonlinearTermPlan{Ti,T}, Rbnm::SparseMatrixCSC,
        Ljb::SparseVector, Lmean, Knm::SparseMatrixCSC, layout::ModeLayout,
        freqindexmap::Vector{Int}) where {Ti,T}
    nc = plan.ncomplex
    Nmodes = length(freqindexmap)
    Nmatrix = plan.nslots ÷ max(length(Ljb.nzval), 1)
    bptr, bsrc, bcoef = backwardjosephsonmap(Ti, T, Rbnm, Ljb, Lmean,
        Nmodes, Nmatrix, freqindexmap, nc)
    length(bcoef) == length(plan.bcoef) || throw(ArgumentError(
        "the Josephson structure moved between points, which a refreshed plan cannot follow; build a new system."))
    copyto!(plan.bcoef, bcoef)
    if hasrealbackward(plan)
        kptr, kidx, kcoef = reallinearmap(Ti, T, Knm, layout)
        (length(kcoef) == length(plan.kcoef) && kidx == Array(plan.kidx)) ||
            throw(ArgumentError("the linear term's structure moved between points, which a refreshed plan cannot follow; build a new system."))
        copyto!(plan.kcoef, kcoef)
    end
    cptr, cidx, ccoef = complexlinearmap(Ti, T, Knm)
    (length(ccoef) == length(plan.ccoef) && cidx == Array(plan.cidx)) ||
        throw(ArgumentError("the linear term's structure moved between points, which a refreshed plan cannot follow; build a new system."))
    copyto!(plan.ccoef, ccoef)
    return plan
end

"""
    tobackend(backend, v::AbstractArray)

Move a host array to the given KernelAbstractions backend.

On a device backend this allocates there and copies. On `CPU()` the vector is
already the type and in the memory the plan wants, so it is adopted as is
rather than duplicated: `allocate` would return an ordinary `Vector` and the
copy would leave the original as garbage, which for a plan of a large circuit
is a transient copy of every index map. Only a dense `Array` is adopted, so
any other `AbstractArray` (a view, a range) is still materialized.

The caller must therefore hand over arrays it does not retain, which every
call site does: they are all freshly built.
"""
function tobackend(backend::Backend, v::AbstractArray)
    d = KernelAbstractions.allocate(backend, eltype(v), size(v)...)
    copyto!(d, v)
    return d
end

tobackend(::CPU, v::Array) = v

"""
    applyrealtocomplex!(xc::AbstractVector, plan::NonlinearTermPlan,
        xr::AbstractVector)
    applycomplextoreal!(xr::AbstractVector, plan::NonlinearTermPlan,
        xc::AbstractVector)

Convert a node vector between the complex representation and the equivalent
real one through the layout the plan carries, as a conflict free kernel on
the plan's backend. These are the device capable counterparts of
[`real_to_complex!`](@ref) and [`complex_to_real!`](@ref) for the two
representations of an [`HBSystem`](@ref) point, which advance a serial cursor
and read a `BitVector`.
"""
function applyrealtocomplex!(xc::AbstractVector, plan::NonlinearTermPlan,
    xr::AbstractVector)
    plan.realtocomplex!(xc, xr, plan.lptr, plan.lwide)
    KernelAbstractions.synchronize(plan.backend)
    return xc
end

function applycomplextoreal!(xr::AbstractVector, plan::NonlinearTermPlan,
    xc::AbstractVector)
    plan.complextoreal!(xr, xc, plan.lptr, plan.lwide)
    KernelAbstractions.synchronize(plan.backend)
    return xr
end

"""
    applyforwardterm!(phimatrix::AbstractArray, plan::NonlinearTermPlan,
        z::AbstractVector)

Apply the forward map of a [`NonlinearTermPlan`](@ref), writing the frequency
domain coefficients of the Josephson junction branch fluxes for the point or
direction `z`. Dispatches on the element type: a real vector is taken in the
equivalent real representation and a complex vector in the complex one.
"""
function applyforwardterm!(phimatrix::AbstractArray,
    plan::NonlinearTermPlan, xr::AbstractVector{<:Real})
    plan.forward!(phimatrix, plan.n1, plan.s1, plan.n2, plan.s2, plan.flags,
        xr)
    KernelAbstractions.synchronize(plan.backend)
    return phimatrix
end

function applyforwardterm!(phimatrix::AbstractArray,
    plan::NonlinearTermPlan, xc::AbstractVector{<:Complex})
    plan.forwardcomplex!(phimatrix, plan.cn1, plan.s1, plan.cn2, plan.s2,
        plan.flags, xc)
    KernelAbstractions.synchronize(plan.backend)
    return phimatrix
end

"""
    hasrealbackward(plan::NonlinearTermPlan)

Whether `plan` carries the real form of the linear term, and so whether the
real representation backward map can be applied. False for a plan built with
`realbackward = false`. The complex representation is always available.
"""
hasrealbackward(plan::NonlinearTermPlan) = !isempty(plan.kptr)

"""
    applybackwardterm!(out::AbstractVector, plan::NonlinearTermPlan,
        phimatrix::AbstractArray, z::AbstractVector;
        addlinearterm::Bool = true)

Apply the backward map of a [`NonlinearTermPlan`](@ref), writing the node
vector. Dispatches on the element type of `out`: a real vector receives the
equivalent real representation and a complex vector the complex one, and `z`
must be in the same representation. With `addlinearterm = true` the frequency
dependent linear term applied to `z` is added in the same pass; with
`addlinearterm = false` only the Josephson contribution is written, as the
second directional derivative requires, and `z` is not read.
"""
function applybackwardterm!(out::AbstractVector{<:Real},
    plan::NonlinearTermPlan, phimatrix::AbstractArray,
    xr::AbstractVector{<:Real}; addlinearterm::Bool = true)
    hasrealbackward(plan) || throw(ArgumentError(
        "this plan was built with realbackward = false, so the real form of "*
        "the linear term is not available and the real representation "*
        "backward map cannot be applied. Build the system with "*
        "realbackward = true to use the real representation entry points."))
    kptr = addlinearterm ? plan.kptr : plan.kptrzero
    plan.backward!(out, plan.bptr, plan.bsrc, plan.bcoef, phimatrix, kptr,
        plan.kidx, plan.kcoef, xr, plan.lptr, plan.lwide)
    KernelAbstractions.synchronize(plan.backend)
    return out
end

function applybackwardterm!(out::AbstractVector{<:Complex},
    plan::NonlinearTermPlan, phimatrix::AbstractArray,
    xc::AbstractVector{<:Complex}; addlinearterm::Bool = true)
    kptr = addlinearterm ? plan.cptr : plan.cptrzero
    plan.backwardcomplex!(out, plan.bptr, plan.bsrc, plan.bcoef, phimatrix,
        kptr, plan.cidx, plan.ccoef, xc)
    KernelAbstractions.synchronize(plan.backend)
    return out
end
