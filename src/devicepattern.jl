# Building the real sparsity structure on a backend.
#
# `expandrealpattern` turns the complex pattern into the real one by walking
# it and emitting each complex entry's block of real entries. On a two tone
# line that is 21.7 million row indices written by a four deep host loop,
# 0.070 s, and it runs inside the rebuild at escalation with the device idle.
# It is also embarrassingly parallel: the number of real entries each complex
# index contributes is fixed by the layouts, so a prefix sum over those counts
# gives every output column a known place to write.

"""
    DeviceSparsePattern{Ti,V}

The sparsity structure of a sparse matrix, compressed by columns, resident on
a backend. Carries no values: it describes where the stored entries are, which
is all the device path needs of the Jacobian's structure.

The Jacobian is held transposed on a device, so this is usually the pattern of
the transpose and its columns are the Jacobian's rows.
"""
struct DeviceSparsePattern{Ti<:Integer,V<:AbstractVector{Ti}}
    colptr::V
    rowval::V
    m::Int
    n::Int
end

Base.size(p::DeviceSparsePattern) = (p.m, p.n)
Base.size(p::DeviceSparsePattern, d::Integer) = size(p)[d]
SparseArrays.nnz(p::DeviceSparsePattern) = length(p.rowval)

# the number of real entries one complex index contributes to each of its real
# columns: the sum of the widths of the complex indices in its list
@kernel function outerwidthsumkernel!(widthsum, @Const(outerptr),
        @Const(outerrowval), @Const(innerptr))
    o = @index(Global)
    @inbounds begin
        s = 0
        for k in outerptr[o]:outerptr[o+1]-1
            i = Int(outerrowval[k])
            s += Int(innerptr[i+1]) - Int(innerptr[i])
        end
        widthsum[o] = s
    end
end

# one work item per real column, writing that column's run of row indices in
# the order the host loop emits them
@kernel function expandpatternkernel!(rowval, @Const(colptr), @Const(outerinv),
        @Const(outerptr), @Const(outerrowval), @Const(innerptr))
    c = @index(Global)
    @inbounds begin
        o = Int(outerinv[c])
        p = Int(colptr[c])
        for k in outerptr[o]:outerptr[o+1]-1
            i = Int(outerrowval[k])
            i0 = Int(innerptr[i])
            for e in 0:(Int(innerptr[i+1]) - i0 - 1)
                rowval[p] = i0 + e
                p += 1
            end
        end
    end
end

"""
    deviceexpandrealpattern(outerptr, outerrowval, outerlayout::ModeLayout,
        innerlayout::ModeLayout, n::Integer, backend)

As [`expandrealpattern`](@ref), built on `backend` and returned there as a
[`DeviceSparsePattern`](@ref).

The number of real entries a complex index contributes is the same for every
one of its real columns, so the column pointer is a gather of per index counts
followed by a prefix sum, and the row indices are then one work item per
column with nothing to coordinate.
"""
function deviceexpandrealpattern(outerptr, outerrowval,
    outerlayout::ModeLayout, innerlayout::ModeLayout, n::Integer, backend)

    d = x -> tobackend(backend, convert(Vector{Int}, x))
    # a pattern built on the backend already is adopted rather than copied
    adopt = x -> x isa Vector ? d(x) : x
    douterptr, douterrowval = adopt(outerptr), adopt(outerrowval)
    dinnerptr = d(collect(innerlayout.ptr))
    douterinv = d(collect(outerlayout.inv))

    widthsum = KernelAbstractions.allocate(backend, Int, n)
    outerwidthsumkernel!(backend, 64)(widthsum, douterptr, douterrowval,
        dinnerptr; ndrange = n)
    KernelAbstractions.synchronize(backend)

    # every real column of complex index `o` holds `widthsum[o]` entries
    rdim = outerlayout.rdim
    wide = KernelAbstractions.allocate(backend, Int, rdim + 1)
    copyto!(view(wide, 2:rdim+1), widthsum[douterinv])
    copyto!(view(wide, 1:1), [1])
    cumsum!(wide, wide)
    nz = Int(Array(view(wide, rdim+1:rdim+1))[1]) - 1

    # the index type is chosen from the count, as the host plan chooses it:
    # narrow indices halve what the assembly and the factorization read
    Ti = nz <= typemax(Int32) ? Int32 : Int
    colptr = Ti === Int ? wide : convert(typeof(similar(wide, Ti)), wide)
    rowval = KernelAbstractions.allocate(backend, Ti, nz)
    expandpatternkernel!(backend, 64)(rowval, colptr, douterinv, douterptr,
        douterrowval, dinnerptr; ndrange = rdim)
    KernelAbstractions.synchronize(backend)
    return DeviceSparsePattern{Ti,typeof(rowval)}(colptr, rowval,
        innerlayout.rdim, rdim)
end

# ---------------------------------------------------------------------------
# the complex pattern, on the backend
# ---------------------------------------------------------------------------
#
# `complexjacobianpattern` collects each column's candidate rows, deduplicates
# them through a workspace the size of the matrix, sorts, and appends. The
# workspace is what makes it serial, and it exists only because the host
# version cannot afford per column state.
#
# It is not needed. The candidates arrive already sorted, in four sequences:
# the Josephson product `adjacency[n2] x activem1[m2]`, whose entries
# `(n1-1)*Nmodes + m1` ascend with `n1` and then `m1` when both factors do,
# and the row lists of the three linear term matrices in that column. So a
# column is a four way merge with deduplication, which needs four cursors and
# no workspace at all.

# the `t`th entry of the Josephson product for a column, or `typemax` past its
# end, enumerated so that the sequence ascends
@inline function jjcandidate(t, na, nm, adjval, a0, amval, m0, Nmodes)
    t >= na * nm && return typemax(Int)
    @inbounds n1 = Int(adjval[a0 + t ÷ nm])
    @inbounds m1 = Int(amval[m0 + t % nm])
    return (n1 - 1) * Nmodes + m1
end

@inline function sortedcandidate(t, ptr, rowval, j)
    @inbounds lo = Int(ptr[j]); hi = Int(ptr[j+1]) - 1
    lo + t > hi && return typemax(Int)
    @inbounds return Int(rowval[lo+t])
end

# the merge itself, shared by the counting and the filling pass so the two
# cannot describe different patterns. `out === nothing` counts only, and the
# branch on it is resolved when the kernel is compiled.
@inline function mergecolumn!(out, p0, cj, Nmodes, adjptr, adjval, amptr,
        amval, lptr, lrow, gptr, grow, cptr, crow)
    @inbounds begin
        n2 = (cj - 1) ÷ Nmodes + 1
        m2 = (cj - 1) % Nmodes + 1
        a0 = Int(adjptr[n2]); na = Int(adjptr[n2+1]) - a0
        m0 = Int(amptr[m2]); nm = Int(amptr[m2+1]) - m0
        t0 = 0; t1 = 0; t2 = 0; t3 = 0
        count = 0
        prev = 0
        while true
            v0 = jjcandidate(t0, na, nm, adjval, a0, amval, m0, Nmodes)
            v1 = sortedcandidate(t1, lptr, lrow, cj)
            v2 = sortedcandidate(t2, gptr, grow, cj)
            v3 = sortedcandidate(t3, cptr, crow, cj)
            v = min(min(v0, v1), min(v2, v3))
            v == typemax(Int) && break
            v0 == v && (t0 += 1)
            v1 == v && (t1 += 1)
            v2 == v && (t2 += 1)
            v3 == v && (t3 += 1)
            if v != prev
                if !isnothing(out)
                    out[p0+count] = v
                end
                count += 1
                prev = v
            end
        end
        return count
    end
end

@kernel function complexpatterncountkernel!(counts, Nmodes, @Const(adjptr),
        @Const(adjval), @Const(amptr), @Const(amval), @Const(lptr),
        @Const(lrow), @Const(gptr), @Const(grow), @Const(cptr), @Const(crow))
    cj = @index(Global)
    @inbounds counts[cj] = mergecolumn!(nothing, 1, cj, Nmodes, adjptr, adjval,
        amptr, amval, lptr, lrow, gptr, grow, cptr, crow)
end

@kernel function complexpatternfillkernel!(rowval, @Const(colptr), Nmodes,
        @Const(adjptr), @Const(adjval), @Const(amptr), @Const(amval),
        @Const(lptr), @Const(lrow), @Const(gptr), @Const(grow), @Const(cptr),
        @Const(crow))
    cj = @index(Global)
    @inbounds mergecolumn!(rowval, Int(colptr[cj]), cj, Nmodes, adjptr, adjval,
        amptr, amval, lptr, lrow, gptr, grow, cptr, crow)
end

# a vector of vectors as a pointer and a value array
function flattenlists(::Type{Ti}, lists) where {Ti<:Integer}
    ptr = Vector{Ti}(undef, length(lists) + 1); ptr[1] = one(Ti)
    val = Ti[]
    for (i, l) in enumerate(lists)
        append!(val, l)
        ptr[i+1] = ptr[i] + length(l)
    end
    return ptr, val
end

"""
    devicecomplexjacobianpattern(n, Nmodes, adjacency, activem1,
        linearmatrices, backend)

As [`complexjacobianpattern`](@ref), built on `backend` and returned there as
a [`DeviceSparsePattern`](@ref).
"""
function devicecomplexjacobianpattern(n::Integer, Nmodes::Integer, adjacency,
    activem1, linearmatrices, backend)

    Ti = Int
    d = x -> tobackend(backend, convert(Vector{Ti}, x))
    adjptr, adjval = flattenlists(Ti, adjacency)
    amptr, amval = flattenlists(Ti, activem1)
    L, G, C = linearmatrices
    args = (d(adjptr), d(adjval), d(amptr), d(amval),
        d(SparseArrays.getcolptr(L)), d(rowvals(L)),
        d(SparseArrays.getcolptr(G)), d(rowvals(G)),
        d(SparseArrays.getcolptr(C)), d(rowvals(C)))

    colptr = KernelAbstractions.allocate(backend, Ti, n + 1)
    complexpatterncountkernel!(backend, 64)(view(colptr, 2:n+1), Int(Nmodes),
        args...; ndrange = n)
    KernelAbstractions.synchronize(backend)
    copyto!(view(colptr, 1:1), [one(Ti)])
    cumsum!(colptr, colptr)
    nz = Int(Array(view(colptr, n+1:n+1))[1]) - 1

    rowval = KernelAbstractions.allocate(backend, Ti, nz)
    complexpatternfillkernel!(backend, 64)(rowval, colptr, Int(Nmodes),
        args...; ndrange = n)
    KernelAbstractions.synchronize(backend)
    return DeviceSparsePattern{Ti,typeof(rowval)}(colptr, rowval, Int(n), Int(n))
end

"""
    transposepattern(p::DeviceSparsePattern, backend)

The transpose of a sparsity structure, on the same backend.

Compressed sparse columns of the transpose are compressed sparse rows of the
original, so this is a stable counting sort of the stored entries by their row:
[`segmentbydest!`](@ref) does exactly that, and its stability is what leaves
the column indices ascending within each row, because a column major traversal
meets them in that order.
"""
function transposepattern(p::DeviceSparsePattern{Ti}, backend) where {Ti}
    nz = nnz(p)
    colptr = KernelAbstractions.allocate(backend, Ti, p.m + 1)
    perm = KernelAbstractions.allocate(backend, Ti, nz)
    segmentbydest!(colptr, perm, p.rowval, backend)
    # the row of the transpose is the column of the original, which is the
    # column each stored entry sits in
    colof = KernelAbstractions.allocate(backend, Ti, nz)
    columnofkernel!(backend, 64)(colof, p.colptr; ndrange = p.n)
    KernelAbstractions.synchronize(backend)
    return DeviceSparsePattern{Ti,typeof(colof)}(colptr, colof[perm], p.n, p.m)
end

# the column each stored entry belongs to, one work item per column
@kernel function columnofkernel!(colof, @Const(colptr))
    c = @index(Global)
    @inbounds for k in Int(colptr[c]):Int(colptr[c+1])-1
        colof[k] = c
    end
end

"""
    realjacobianstructure(Amatrixindices::Matrix, Amatrixconjindices::Matrix,
        Ljb::SparseVector, Rbnm::SparseMatrixCSC, Nmodes::Integer,
        Nbranches::Integer, invLnm, Gnm, Cnm, rl::ModeLayout, cl::ModeLayout,
        ::Type{T} = Float64; transposed = false, backend = CPU())

The sparsity structure of the real Jacobian, and the branch incidence the
assembly needs with it: `(P, nodesandsigns)`.

`P` is a `SparseMatrixCSC` with zero values on a host and a
[`DeviceSparsePattern`](@ref) on a backend, which is what a device
factorization can be built from without a conversion. `transposed` returns the
structure of the transpose, whose stored order is the row major order of the
Jacobian.

Nothing sized by the circuit is built on the host on a backend: the complex
pattern is built and transposed there and expanded there, and the only host
inputs are the node adjacency and the active mode rows.
"""
function realjacobianstructure(Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, Ljb::SparseVector, Rbnm::SparseMatrixCSC,
    Nmodes::Integer, Nbranches::Integer, invLnm, Gnm, Cnm, rl::ModeLayout,
    cl::ModeLayout, ::Type{T} = Float64; transposed::Bool = false,
    backend = CPU()) where {T<:Real}

    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
    n = size(Rbnm, 2)
    C = devicecomplexjacobianpattern(n, Nmodes,
        jjnodeadjacency(Ljb, nodesandsigns, n ÷ Nmodes),
        activemoderows(Nmodes, Amatrixindices, Amatrixconjindices),
        (invLnm, Gnm, Cnm), backend)
    Cx = transposed ? transposepattern(C, backend) : C
    Pd = deviceexpandrealpattern(Cx.colptr, Cx.rowval, rl, cl, n, backend)
    P = backend isa CPU ?
        SparseMatrixCSC(Pd.m, Pd.n, Array(Pd.colptr), Array(Pd.rowval),
            zeros(T, nnz(Pd))) : Pd
    return P, nodesandsigns
end

"""
    segmentbydest!(seg, perm, dest, backend)

Group `dest` into contiguous segments by value, writing the segment pointer
into `seg` and, into `perm`, the positions of `dest` ordered by segment.

This is the device form of [`segmentscatterbydest`](@ref), and it is a
permutation rather than a copy: the caller gathers whatever payload it has
through `perm`, which keeps this independent of how many arrays travel with
the destinations.

The host version is a stable counting sort, and stability is what makes the
assembly reproducible: it fixes the order in which the contributions to one
entry are added, and floating point addition is not associative. A device
histogram with an atomic cursor is not stable, so the order is restored
afterwards by sorting each segment. That is cheap because the segments are
tiny -- 28.8 million contributions over 21.7 million entries is an average of
1.33 -- and it is why this is not done with a general sort, which measured
0.315 s on 28.8 million keys, as much as the host counting sort it would
replace.
"""
function segmentbydest!(seg::AbstractVector, perm::AbstractVector,
    dest::AbstractVector, backend)

    nbins = length(seg) - 1
    length(perm) == length(dest) || throw(DimensionMismatch(
        lazy"`perm` has length $(length(perm)) but `dest` has $(length(dest))."))

    # count, then turn the counts into the segment pointer
    fill!(seg, zero(eltype(seg)))
    counthistogram!(backend)(seg, dest; ndrange = length(dest))
    KernelAbstractions.synchronize(backend)
    exclusiveprefix!(seg)

    # scatter the positions into their segments, then put each segment back
    # into the order the positions were emitted in
    cursor = copy(seg)
    scatterpositions!(backend)(perm, cursor, dest; ndrange = length(dest))
    KernelAbstractions.synchronize(backend)
    sortsegments!(backend)(perm, seg; ndrange = nbins)
    KernelAbstractions.synchronize(backend)
    return seg, perm
end

# one work item per contribution, counting into the bin one past its
# destination so that the prefix sum below lands the pointer in place
@kernel function counthistogram!(seg, @Const(dest))
    k = @index(Global)
    @inbounds Atomix.@atomic seg[dest[k]+1] += one(eltype(seg))
end

# `seg` arrives holding the count of each destination at one past it; leave it
# holding the first position of each. Runs on the host: it is one pass over the
# stored entries and it is serial by nature, but it is also the only part of
# this that is, and it is a hundredth of what it replaces.
function exclusiveprefix!(seg::AbstractVector)
    s = Array(seg)
    s[1] = one(eltype(s))
    @inbounds for k in 2:length(s)
        s[k] += s[k-1]
    end
    copyto!(seg, s)
    return seg
end

@kernel function scatterpositions!(perm, cursor, @Const(dest))
    k = @index(Global)
    @inbounds begin
        d = dest[k]
        # the atomic add returns the value after the increment, so the slot
        # this work item claimed is one below it
        one_ = one(eltype(cursor))
        slot = (Atomix.@atomic cursor[d] += one_) - one_
        perm[slot] = k
    end
end

# an insertion sort per segment, restoring the emission order the atomic
# cursor above scrambled. The segments average between one and two entries, so
# this is a handful of comparisons for the whole matrix.
@kernel function sortsegments!(perm, @Const(seg))
    q = @index(Global)
    @inbounds begin
        lo = seg[q]
        hi = seg[q+1] - 1
        for a in lo+1:hi
            v = perm[a]
            b = a - 1
            while b >= lo && perm[b] > v
                perm[b+1] = perm[b]
                b -= 1
            end
            perm[b+1] = v
        end
    end
end


"""
    DeviceValuedSparseMatrix{Tv,Ti,V} <: AbstractMatrix{Tv}

A sparse matrix whose *structure* lives on the host and whose *values* live on
a device, in the stored order of a compressed sparse row matrix built from that
structure.

The structure is held as the **transpose**, because the compressed columns of
the transpose are the compressed rows of this matrix: `patterntranspose` is
therefore a row pointer and a column index array already, and a device sparse
matrix can be built from it with no conversion at all. That is also the form
[`planrealjacobian`](@ref) produces directly when asked for it with
`transposed = true`, which is why nothing here permutes anything.

`patterntranspose` carries the sparsity structure only; its stored values are
whatever they were last set to and must not be read. Use [`nonzeros`](@ref) to
reach the device values, and [`rowpointer`](@ref) and
[`columnindices`](@ref) to reach the structure.

A factorization which does not have a method for this type will fail on it
rather than silently reading the pattern, which is deliberate: the values in
`patterntranspose` are stale by construction.
"""
struct DeviceValuedSparseMatrix{Tv,P,V<:AbstractVector{Tv}} <: AbstractMatrix{Tv}
    patterntranspose::P
    nzval::V
end

"""
    rowpointer(A::DeviceValuedSparseMatrix)

The compressed sparse row pointer of `A`, which is the column pointer of the
transpose it stores.
"""
rowpointer(A::DeviceValuedSparseMatrix{<:Any,<:SparseMatrixCSC}) =
    SparseArrays.getcolptr(A.patterntranspose)
rowpointer(A::DeviceValuedSparseMatrix{<:Any,<:DeviceSparsePattern}) =
    A.patterntranspose.colptr

"""
    columnindices(A::DeviceValuedSparseMatrix)

The compressed sparse row column indices of `A`, which are the row indices of
the transpose it stores.
"""
columnindices(A::DeviceValuedSparseMatrix{<:Any,<:SparseMatrixCSC}) =
    rowvals(A.patterntranspose)
columnindices(A::DeviceValuedSparseMatrix{<:Any,<:DeviceSparsePattern}) =
    A.patterntranspose.rowval

Base.size(A::DeviceValuedSparseMatrix) = reverse(size(A.patterntranspose))
Base.size(A::DeviceValuedSparseMatrix, d::Integer) = size(A)[d]
SparseArrays.nnz(A::DeviceValuedSparseMatrix) = nnz(A.patterntranspose)
SparseArrays.nonzeros(A::DeviceValuedSparseMatrix) = A.nzval

# reading an entry would read the stale host values, so refuse rather than
# return something wrong
function Base.getindex(A::DeviceValuedSparseMatrix, ::Integer, ::Integer)
    throw(ArgumentError("A DeviceValuedSparseMatrix holds its values on a device; its host structure carries the sparsity pattern only and its stored values are stale. Index the device vector returned by `nonzeros` instead."))
end

# one work item per row of the matrix, which the stored transpose makes a
# compressed column of the structure
@kernel function devicecsrmulkernel!(w, @Const(colptr), @Const(rowval),
        @Const(nzval), @Const(x))
    i = @index(Global)
    T = eltype(w)
    @inbounds begin
        acc = zero(T)
        for k in Int(colptr[i]):Int(colptr[i+1])-1
            acc += nzval[k] * x[Int(rowval[k])]
        end
        w[i] = acc
    end
end

"""
    mul!(w::AbstractVector{<:BlasFloat}, A::DeviceValuedSparseMatrix,
        x::AbstractVector)

`w = A*x` on the backend `A`'s values live on.

`A` stores its structure as the transpose, so a compressed column of that
structure is a row of `A` and the product is a row wise contraction with no
atomics.

The destination is restricted to the element types a device backend actually
stores, which is also what keeps this from being ambiguous with the `mul!`
methods packages define for their own scalar types.
"""
function LinearAlgebra.mul!(w::AbstractVector{<:LinearAlgebra.BlasFloat},
    A::DeviceValuedSparseMatrix, x::AbstractVector)
    length(w) == size(A, 1) || throw(DimensionMismatch(
        lazy"`w` has length $(length(w)) but `A` has $(size(A,1)) rows."))
    length(x) == size(A, 2) || throw(DimensionMismatch(
        lazy"`x` has length $(length(x)) but `A` has $(size(A,2)) columns."))
    backend = KernelAbstractions.get_backend(A.nzval)
    devicecsrmulkernel!(backend, 64)(w, rowpointer(A), columnindices(A),
        A.nzval, x; ndrange = size(A, 1))
    KernelAbstractions.synchronize(backend)
    return w
end

"""
    dot(y::AbstractVector{<:BlasFloat}, A::DeviceValuedSparseMatrix,
        x::AbstractVector)

`y' * A * x`, the model slope the line search of [`nlsolve!`](@ref) needs.
"""
function LinearAlgebra.dot(y::AbstractVector{<:LinearAlgebra.BlasFloat},
    A::DeviceValuedSparseMatrix, x::AbstractVector)
    w = similar(y, size(A, 1))
    mul!(w, A, x)
    return dot(y, w)
end
