

"""
    ModeLayout(isreal::AbstractVector{Bool}, dim::Integer, ::Type{Ti}=Int)
    ModeLayout(realindices, nmodes::Integer, dim::Integer, ::Type{Ti}=Int)

Layout of one axis of length `dim` (a complex dimension), built from a length-
`nmodes` mask of which modes are real. `dim` must be an integer multiple of
`nmodes`. Complex index `i` owns real slots `ptr[i]:ptr[i+1]-1`; `inv` maps a
real slot back to its complex index; `isfirst` marks the slots that start a
mode; `rdim` is the resulting real dimension.

`w` is a bit per index recording whether that mode is real.
"""
struct ModeLayout{Ti<:Integer}
    nmodes::Int
    dim::Int
    rdim::Int
    nreal::Int
    isreal::BitVector
    ptr::Vector{Ti}       # length dim+1
    inv::Vector{Ti}       # length rdim
    isfirst::Vector{Bool} # length rdim
    w::BitVector          # length dim, true where the mode is real
end

@inline _rowwidth(w::BitVector, i::Integer) = @inbounds 2 - w[i]
@inline _width(L::ModeLayout, i::Integer)   = @inbounds Int(L.ptr[i+1] - L.ptr[i])

function ModeLayout(isreal::AbstractVector{Bool}, dim::Integer, ::Type{Ti} = Int) where {Ti<:Integer}
    nmodes = length(isreal)
    nmodes >= 1 || throw(ArgumentError("nmodes must be at least 1"))
    dim >= 0 || throw(ArgumentError("dim must be nonnegative"))
    dim % nmodes == 0 || throw(DimensionMismatch(
        "dimension $dim is not an integer multiple of nmodes = $nmodes"))
    nreal = count(isreal)
    rdim  = (dim ÷ nmodes) * (2nmodes - nreal)

    ptr     = Vector{Ti}(undef, dim + 1)
    inv     = Vector{Ti}(undef, rdim)
    isfirst = Vector{Bool}(undef, rdim)
    w       = falses(dim)
    p = 1
    @inbounds for i in 1:dim
        r = isreal[(i - 1) % nmodes + 1]
        ptr[i] = p
        isfirst[p] = true
        inv[p] = i
        w[i] = r
        if !r
            inv[p+1] = i;  isfirst[p+1] = false
            p += 2
        else
            p += 1
        end
    end
    @inbounds ptr[dim+1] = p
    return ModeLayout(nmodes, dim, rdim, nreal, BitVector(isreal), ptr, inv, isfirst, w)
end

function ModeLayout(realindices, nmodes::Integer, dim::Integer, ::Type{Ti} = Int) where {Ti<:Integer}
    mask = falses(nmodes)
    for r in realindices
        1 <= r <= nmodes || throw(ArgumentError("real index $r outside 1:$nmodes"))
        mask[r] && throw(ArgumentError("duplicate real index $r"))
        mask[r] = true
    end
    return ModeLayout(mask, dim, Ti)
end

# per-entry scale factor; `cf` is loop-invariant per column
@inline _colfac(cs::T, wc) where {T} = wc == 1 ? cs : one(T)
@inline _rowfac(rs::T, wr) where {T} = ifelse(wr == 1, rs, one(T))

#  vectors and dense matrices

"""
    realdim(dim, isreal) -> Int

Length of the real form of a `dim`-long complex axis under mask `isreal`.
"""
function realdim(dim::Integer, isreal::AbstractVector{Bool})
    nm = length(isreal)
    nm >= 1 || throw(ArgumentError("mask must not be empty"))
    dim % nm == 0 || throw(DimensionMismatch(
        "dimension $dim is not an integer multiple of nmodes = $nm"))
    return (dim ÷ nm) * (2nm - count(isreal))
end

"""
    complexdim(rdim, isreal) -> Int

Inverse of [`realdim`](@ref).
"""
function complexdim(rdim::Integer, isreal::AbstractVector{Bool})
    nm = length(isreal)
    per = 2nm - count(isreal)
    rdim % per == 0 || throw(DimensionMismatch(
        "real dimension $rdim is not an integer multiple of $per"))
    return (rdim ÷ per) * nm
end

@inline _next(t, n) = ifelse(t == n, 1, t + 1)

# `rmask` describes the row axis (the output vector) and `cmask` the column axis
# (the input vector); they are the dense analogue of the `rl`/`cl` layout pair.
# When the two axes share a mode pattern pass it once. 
complex_to_real(A::AbstractMatrix{<:Complex}, mask::AbstractVector{Bool}; kw...) =
    complex_to_real(A, mask, mask; kw...)
complex_to_real!(Ar::AbstractMatrix{<:Real}, A::AbstractMatrix{<:Complex},
         mask::AbstractVector{Bool}; kw...) = complex_to_real!(Ar, A, mask, mask; kw...)
real_to_complex(Ar::AbstractMatrix{<:Real}, mask::AbstractVector{Bool}; kw...) =
    real_to_complex(Ar, mask, mask; kw...)
real_to_complex!(A::AbstractMatrix{<:Complex}, Ar::AbstractMatrix{<:Real},
            mask::AbstractVector{Bool}; kw...) = real_to_complex!(A, Ar, mask, mask; kw...)

#  vectors

"""
    complex_to_real!(xr, xc, isreal; conj_input=false, realscale=1) -> xr

Complex -> real. The imaginary part of a real mode is never read. `realscale`
multiplies the real modes; `conj_input=true` converts `conj(xc)` instead.
"""
function complex_to_real!(xr::AbstractVector{T}, xc::AbstractVector{Complex{T}},
               isreal::AbstractVector{Bool}; conj_input::Bool = false,
               realscale = one(T)) where {T<:Real}
    length(xr) == realdim(length(xc), isreal) || throw(DimensionMismatch(
        "destination has length $(length(xr)), expected $(realdim(length(xc), isreal))"))
    nm = length(isreal)
    s  = T(realscale)
    qs = conj_input ? -one(T) : one(T)
    p, t = 1, 1
    @inbounds for i in eachindex(xc)
        a = xc[i]
        if isreal[t]
            xr[p] = s * real(a)
            p += 1
        else
            xr[p]   = real(a)
            xr[p+1] = qs * imag(a)
            p += 2
        end
        t = _next(t, nm)
    end
    return xr
end

"""
    real_to_complex!(xc, xr, isreal; conj_input=false, realscale=1) -> xc

Real -> complex. The imaginary part of a real mode is written as an explicit
zero, so the result does not depend on what was in `xc` beforehand.
"""
function real_to_complex!(xc::AbstractVector{Complex{T}}, xr::AbstractVector{T},
                 isreal::AbstractVector{Bool}; conj_input::Bool = false,
                 realscale = one(T)) where {T<:Real}
    length(xr) == realdim(length(xc), isreal) || throw(DimensionMismatch(
        "source has length $(length(xr)), expected $(realdim(length(xc), isreal))"))
    nm = length(isreal)
    s  = T(realscale)
    qs = conj_input ? -one(T) : one(T)
    p, t = 1, 1
    @inbounds for i in eachindex(xc)
        if isreal[t]
            xc[i] = Complex(s * xr[p], zero(T))
            p += 1
        else
            xc[i] = Complex(xr[p], qs * xr[p+1])
            p += 2
        end
        t = _next(t, nm)
    end
    return xc
end

complex_to_real(xc::AbstractVector{Complex{T}}, isreal::AbstractVector{Bool}; kw...) where {T} =
    complex_to_real!(Vector{T}(undef, realdim(length(xc), isreal)), xc, isreal; kw...)
real_to_complex(xr::AbstractVector{T}, isreal::AbstractVector{Bool}; kw...) where {T} =
    real_to_complex!(Vector{Complex{T}}(undef, complexdim(length(xr), isreal)), xr, isreal; kw...)

#  complex -> real

"""
    complex_to_real!(Ar, A, rmask, cmask; conj_input=false,
             realrowscale=1, realcolscale=1) -> Ar

In-place real form of `x -> A*x`, or of `x -> A*conj(x)` with `conj_input=true`,
so that `complex_to_real(A, rm, cm) * complex_to_real(x, cm) == complex_to_real(A*x, rm)`.
"""
function complex_to_real!(Ar::AbstractMatrix{T}, A::AbstractMatrix{Complex{T}},
                  rmask::AbstractVector{Bool}, cmask::AbstractVector{Bool};
                  conj_input::Bool = false, realrowscale = one(T),
                  realcolscale = one(T)) where {T<:Real}
    m, n = size(A)
    _checkshape(Ar, m, n, rmask, cmask)
    nr, nc = length(rmask), length(cmask)
    rs, cs = T(realrowscale), T(realcolscale)
    c0, tc = 1, 1
    @inbounds for j in 1:n
        creal = cmask[tc]
        cf = creal ? cs : one(T)
        r0, tr = 1, 1
        if creal
            for i in 1:m
                wr = 2 - rmask[tr]
                a  = A[i,j] * (cf * ifelse(wr == 1, rs, one(T)))
                Ar[r0, c0]      = real(a)
                Ar[r0+wr-1, c0] = ifelse(wr == 2, imag(a), real(a))
                r0 += wr
                tr = _next(tr, nr)
            end
            c0 += 1
        else
            for i in 1:m
                wr = 2 - rmask[tr]
                a  = A[i,j] * (cf * ifelse(wr == 1, rs, one(T)))
                d  = conj_input ? -a : a
                Ar[r0, c0]        = real(a)
                Ar[r0+wr-1, c0]   = ifelse(wr == 2, imag(a), real(a))
                Ar[r0, c0+1]      = -imag(d)
                Ar[r0+wr-1, c0+1] = ifelse(wr == 2, real(d), -imag(d))
                r0 += wr
                tr = _next(tr, nr)
            end
            c0 += 2
        end
        tc = _next(tc, nc)
    end
    return Ar
end

"""
    complex_to_real(A, rmask, cmask; conj_input=false, realrowscale=1, realcolscale=1)
        -> Matrix{T}
"""
complex_to_real(A::AbstractMatrix{Complex{T}}, rmask::AbstractVector{Bool},
        cmask::AbstractVector{Bool}; kw...) where {T<:Real} =
    complex_to_real!(Matrix{T}(undef, realdim(size(A,1), rmask), realdim(size(A,2), cmask)),
             A, rmask, cmask; kw...)



#  real -> complex

"""
    real_to_complex!(A, Ar, rmask, cmask; conj_input=false) -> A

Inverse of `complex_to_real!` (unit scales only). For a complex row mode both parts sit
in the first column of the block pair; for a real row mode with a complex column
mode the imaginary part lives in the second column, with a sign set by
`conj_input`. A 1x1 block has no stored `q` and comes back with an exact zero
imaginary part.
"""
function real_to_complex!(A::AbstractMatrix{Complex{T}}, Ar::AbstractMatrix{T},
                     rmask::AbstractVector{Bool}, cmask::AbstractVector{Bool};
                     conj_input::Bool = false) where {T<:Real}
    m, n = size(A)
    _checkshape(Ar, m, n, rmask, cmask)
    nr, nc = length(rmask), length(cmask)
    qs = conj_input ? one(T) : -one(T)
    c0, tc = 1, 1
    @inbounds for j in 1:n
        creal = cmask[tc]
        r0, tr = 1, 1
        if creal
            for i in 1:m
                wr = 2 - rmask[tr]
                A[i,j] = Complex(Ar[r0, c0],
                                 ifelse(wr == 2, Ar[r0+wr-1, c0], zero(T)))
                r0 += wr
                tr = _next(tr, nr)
            end
            c0 += 1
        else
            for i in 1:m
                wr = 2 - rmask[tr]
                A[i,j] = Complex(Ar[r0, c0],
                                 ifelse(wr == 2, Ar[r0+wr-1, c0], qs * Ar[r0, c0+1]))
                r0 += wr
                tr = _next(tr, nr)
            end
            c0 += 2
        end
        tc = _next(tc, nc)
    end
    return A
end

"""
    real_to_complex(Ar, rmask, cmask; conj_input=false) -> Matrix{Complex{T}}
"""
real_to_complex(Ar::AbstractMatrix{T}, rmask::AbstractVector{Bool},
           cmask::AbstractVector{Bool}; kw...) where {T<:Real} =
    real_to_complex!(Matrix{Complex{T}}(undef, complexdim(size(Ar,1), rmask),
                                   complexdim(size(Ar,2), cmask)),
                Ar, rmask, cmask; kw...)

@inline function _checkshape(Ar, m, n, rmask, cmask)
    want = (realdim(m, rmask), realdim(n, cmask))
    size(Ar) == want || throw(DimensionMismatch(
        "real matrix is $(size(Ar)), expected $want"))
    nothing
end

#  sparse matrices

_layouts(mask::AbstractVector{Bool}, m::Integer, n::Integer) =
    m == n ? ((L = ModeLayout(mask, m)); (L, L)) :
             (ModeLayout(mask, m), ModeLayout(mask, n))

function complex_to_real(A::SparseMatrixCSC{Complex{T},Ti}, mask::AbstractVector{Bool},
                 ::Type{Tj} = Ti; kw...) where {T<:Real,Ti,Tj<:Integer}
    rl, cl = _layouts(mask, size(A)...)
    return complex_to_real(A, rl, cl, Tj; kw...)
end

function complex_to_real!(Ar::SparseMatrixCSC{T}, A::SparseMatrixCSC{Complex{T}},
                  mask::AbstractVector{Bool}; kw...) where {T<:Real}
    rl, cl = _layouts(mask, size(A)...)
    return complex_to_real!(Ar, A, rl, cl; kw...)
end


function real_to_complex(Ar::SparseMatrixCSC{T,Ti}, mask::AbstractVector{Bool},
                    ::Type{Tj} = Ti; kw...) where {T<:Real,Ti,Tj<:Integer}
    rl, cl = _layouts(mask, complexdim(size(Ar,1), mask), complexdim(size(Ar,2), mask))
    return real_to_complex(Ar, rl, cl, Tj; kw...)
end

function real_to_complex!(A::SparseMatrixCSC{Complex{T}}, Ar::SparseMatrixCSC{T},
                     mask::AbstractVector{Bool}; kw...) where {T<:Real}
    rl, cl = _layouts(mask, size(A)...)
    return real_to_complex!(A, Ar, rl, cl; kw...)
end

function is_complex_to_real_pattern(Ar::SparseMatrixCSC, A::SparseMatrixCSC, mask::AbstractVector{Bool})
    rl, cl = _layouts(mask, size(A)...)
    return is_complex_to_real_pattern(Ar, A, rl, cl)
end

#  complex -> real

"""
    complex_to_real(A, rowlayout, collayout, ::Type{Tj}=Ti;
            conj_input=false, realrowscale=1, realcolscale=1) -> SparseMatrixCSC{T,Tj}

Real form of `x -> A*x`, or of `x -> A*conj(x)` with `conj_input=true`. `Tj`
selects the index type of the result; `Int32` halves `rowval` and is worth using
whenever the dimensions fit. Runs one O(nnz) pass to size the output, then one
to fill it.
"""
function complex_to_real(A::SparseMatrixCSC{Complex{T},Ti}, rl::ModeLayout, cl::ModeLayout,
                 ::Type{Tj} = Ti; conj_input::Bool = false,
                 realrowscale = one(T), realcolscale = one(T)) where {T<:Real,Ti,Tj<:Integer}
    _checkdims(A, rl, cl)
    n = size(A, 2)
    Ap, Ai, Av = SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)
    rptr, cptr = rl.ptr, cl.ptr
    rs, cs = T(realrowscale), T(realcolscale)

    total = 0
    @inbounds for j in 1:n
        S = 0
        for idx in Ap[j]:Ap[j+1]-1
            i = Ai[idx]
            S += rptr[i+1] - rptr[i]
        end
        total += (cptr[j+1] - cptr[j]) * S
    end

    colptr = Vector{Tj}(undef, cl.rdim + 1)
    rowval = Vector{Tj}(undef, total)
    nzval  = Vector{T}(undef, total)
    colptr[1] = 1
    k = 1
    @inbounds for j in 1:n
        c0 = cptr[j]
        wc = cptr[j+1] - c0
        cf = _colfac(cs, wc)
        k0 = k
        for idx in Ap[j]:Ap[j+1]-1
            i  = Ai[idx]
            r0 = rptr[i]
            wr = rptr[i+1] - r0
            a  = Av[idx] * (cf * _rowfac(rs, wr))
            rowval[k]      = r0
            rowval[k+wr-1] = ifelse(wr == 2, r0 + one(r0), r0)
            nzval[k]       = real(a)
            nzval[k+wr-1]  = ifelse(wr == 2, imag(a), real(a))
            k += wr
        end
        S = k - k0
        colptr[c0+1] = k
        if wc == 2
            copyto!(rowval, k, rowval, k0, S)
            for idx in Ap[j]:Ap[j+1]-1
                i  = Ai[idx]
                wr = rptr[i+1] - rptr[i]
                d  = _second(Av[idx] * (cf * _rowfac(rs, wr)), conj_input)
                nzval[k]      = -imag(d)
                nzval[k+wr-1] = ifelse(wr == 2, real(d), -imag(d))
                k += wr
            end
            colptr[c0+2] = k
        end
    end
    return SparseMatrixCSC{T,Tj}(rl.rdim, cl.rdim, colptr, rowval, nzval)
end

# value whose (-imag, real) pair is the second column of the block
@inline _second(a::Complex, conj_input::Bool) = conj_input ? -a : a

"""
    complex_to_real!(Ar, A, rowlayout, collayout; conj_input=false,
             realrowscale=1, realcolscale=1) -> Ar

In-place `complex_to_real`. `Ar` must already have exactly the pattern `complex_to_real` would
produce; only the shapes are verified. Allocation-free. The pattern does not
depend on `conj_input` or on the scales, so one `Ar` can be refilled any way.
"""
function complex_to_real!(Ar::SparseMatrixCSC{T}, A::SparseMatrixCSC{Complex{T}},
                  rl::ModeLayout, cl::ModeLayout; conj_input::Bool = false,
                  realrowscale = one(T), realcolscale = one(T)) where {T<:Real}
    _checkdims(A, rl, cl)
    _checkdest(Ar, rl, cl)
    Ap, Ai, Av = SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)
    Rp, Rv = SparseArrays.getcolptr(Ar), nonzeros(Ar)
    rw, cptr = rl.w, cl.ptr
    rs, cs = T(realrowscale), T(realcolscale)
    @inbounds for j in 1:size(A, 2)
        c0 = cptr[j]
        wc = cptr[j+1] - c0
        cf = _colfac(cs, wc)
        k  = Rp[c0]
        k2 = wc == 2 ? Rp[c0+1] : k
        for idx in Ap[j]:Ap[j+1]-1
            wr = _rowwidth(rw, Ai[idx])
            a  = Av[idx] * (cf * _rowfac(rs, wr))
            Rv[k]      = real(a)
            Rv[k+wr-1] = ifelse(wr == 2, imag(a), real(a))
            if wc == 2
                d = _second(a, conj_input)
                Rv[k2]      = -imag(d)
                Rv[k2+wr-1] = ifelse(wr == 2, real(d), -imag(d))
                k2 += wr
            end
            k += wr
        end
    end
    return Ar
end


#  real -> complex

"""
    real_to_complex(Ar, rowlayout, collayout, ::Type{Tj}=Ti; conj_input=false)

Inverse of `complex_to_real` (unit scales only). For a complex row mode both parts sit
in the first column of the block pair; for a real row mode with a complex column
mode the imaginary part lives in the second column, with a sign that depends on
`conj_input`, and is read from there using the fact that paired columns share
their `rowval`. A 1x1 block has no stored `q` and comes back with an exact zero
imaginary part.
"""
function real_to_complex(Ar::SparseMatrixCSC{T,Ti}, rl::ModeLayout, cl::ModeLayout,
                    ::Type{Tj} = Ti; conj_input::Bool = false) where {T<:Real,Ti,Tj<:Integer}
    _checkdest(Ar, rl, cl)
    Rp, Ri, Rv = SparseArrays.getcolptr(Ar), rowvals(Ar), nonzeros(Ar)
    cptr, rinv, isfirst = cl.ptr, rl.inv, rl.isfirst
    qs = conj_input ? one(T) : -one(T)

    # Scan every slot and count the ones that start a mode
    total = 0
    @inbounds for j in 1:cl.dim
        c0 = cptr[j]
        for k in Rp[c0]:Rp[c0+1]-1
            total += isfirst[Ri[k]]
        end
    end

    colptr = Vector{Tj}(undef, cl.dim + 1)
    rowval = Vector{Tj}(undef, total)
    nzval  = Vector{Complex{T}}(undef, total)
    colptr[1] = 1
    t = 0
    @inbounds for j in 1:cl.dim
        c0 = cptr[j]
        wc = cptr[j+1] - c0
        k2 = wc == 2 ? Rp[c0+1] : Rp[c0]
        # the end of the paired column, used to bound the lockstep cursor `k2`
        k2end = wc == 2 ? Rp[c0+2] : Rp[c0+1]
        for k in Rp[c0]:Rp[c0+1]-1
            r = Ri[k]
            if wc == 2
                # the two real columns of a complex column carry the same row
                # pattern, so the partner entry must sit at the matching offset
                # of the second column. without this check a shorter or
                # differently patterned second column makes `Rv[k2]` read past
                # the column, or past the end of `Rv` entirely, which
                # `@inbounds` does not catch.
                (k2 < k2end && Ri[k2] == r) || throw(ArgumentError(
                    lazy"the two real columns of complex column $(j) do not share a row pattern; `Ar` is not a pattern produced by `complex_to_real`."))
            end
            if isfirst[r]
                t += 1
                rowval[t] = rinv[r]
                nzval[t]  = Complex(Rv[k], wc == 2 ? qs * Rv[k2] : zero(T))
            else
                # the imaginary part of a complex row mode is folded into the
                # entry written for that mode's first real slot, so that slot
                # must already be stored in this column. if it is not, `t`
                # either still refers to a previous mode or is zero, and
                # writing here would silently corrupt an unrelated entry or
                # write out of bounds (`@inbounds` elides the check).
                (t >= 1 && rowval[t] == rinv[r]) || throw(ArgumentError(
                    lazy"the second real slot of the mode at row $(rinv[r]) of column $(j) is stored but the first is not; `Ar` is not a pattern produced by `complex_to_real`."))
                nzval[t] = Complex(real(nzval[t]), Rv[k])
            end
            k2 += 1
        end
        # the paired column must not carry entries the first one lacks
        wc == 2 && k2 != k2end && throw(ArgumentError(
            lazy"the two real columns of complex column $(j) do not share a row pattern; `Ar` is not a pattern produced by `complex_to_real`."))
        colptr[j+1] = t + 1
    end
    return SparseMatrixCSC{Complex{T},Tj}(rl.dim, cl.dim, colptr, rowval, nzval)
end

"""
    real_to_complex!(A, Ar, rowlayout, collayout; conj_input=false) -> A

In-place `real_to_complex`. `A` must have exactly the pattern `Ar` was built from;
only the shapes are verified. Pass the same `conj_input` used to build `Ar`.
"""
function real_to_complex!(A::SparseMatrixCSC{Complex{T}}, Ar::SparseMatrixCSC{T},
                     rl::ModeLayout, cl::ModeLayout; conj_input::Bool = false) where {T<:Real}
    _checkdims(A, rl, cl)
    _checkdest(Ar, rl, cl)
    Ap, Ai, Av = SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)
    Rp, Rv = SparseArrays.getcolptr(Ar), nonzeros(Ar)
    rw, cptr = rl.w, cl.ptr
    qs = conj_input ? one(T) : -one(T)
    @inbounds for j in 1:size(A, 2)
        c0 = cptr[j]
        k  = Rp[c0]
        if cptr[j+1] - c0 == 2
            k2 = Rp[c0+1]
            for idx in Ap[j]:Ap[j+1]-1
                wr = _rowwidth(rw, Ai[idx])
                Av[idx] = Complex(Rv[k], ifelse(wr == 2, Rv[k+wr-1], qs * Rv[k2]))
                k += wr;  k2 += wr
            end
        else
            for idx in Ap[j]:Ap[j+1]-1
                wr = _rowwidth(rw, Ai[idx])
                Av[idx] = Complex(Rv[k], ifelse(wr == 2, Rv[k+wr-1], zero(T)))
                k += wr
            end
        end
    end
    return A
end

# -------------------------------------------------------------------
#  checks
# -------------------------------------------------------------------

@inline function _checkdims(A, rl::ModeLayout, cl::ModeLayout)
    rl.dim == size(A, 1) || throw(DimensionMismatch(
        "row layout is for dimension $(rl.dim), matrix has $(size(A,1))"))
    cl.dim == size(A, 2) || throw(DimensionMismatch(
        "column layout is for dimension $(cl.dim), matrix has $(size(A,2))"))
    nothing
end

@inline function _checkdest(Ar, rl::ModeLayout, cl::ModeLayout)
    size(Ar) == (rl.rdim, cl.rdim) || throw(DimensionMismatch(
        "real matrix is $(size(Ar)), expected $((rl.rdim, cl.rdim))"))
    nothing
end

"""
    is_complex_to_real_pattern(Ar, A, rowlayout, collayout) -> Bool

Returns true if `Ar` has exactly the size and pattern 
`complex_to_real(A, rl, cl)` produces. Intended for use in setup.
"""
function is_complex_to_real_pattern(Ar::SparseMatrixCSC, A::SparseMatrixCSC,
                            rl::ModeLayout, cl::ModeLayout)
    (rl.dim == size(A,1) && cl.dim == size(A,2)) || return false
    size(Ar) == (rl.rdim, cl.rdim) || return false
    Ap, Ai = SparseArrays.getcolptr(A), rowvals(A)
    Rp, Ri = SparseArrays.getcolptr(Ar), rowvals(Ar)
    rptr, cptr = rl.ptr, cl.ptr
    @inbounds for j in 1:size(A,2)
        c0 = cptr[j]
        for h in 0:(cptr[j+1] - c0 - 1)
            k = Rp[c0+h]
            for idx in Ap[j]:Ap[j+1]-1
                i  = Ai[idx]
                r0 = rptr[i];  wr = rptr[i+1] - r0
                (k + wr - 1 < Rp[c0+h+1] && Ri[k] == r0) || return false
                wr == 2 && Ri[k+1] != r0 + 1 && return false
                k += wr
            end
            k == Rp[c0+h+1] || return false
        end
    end
    return true
end