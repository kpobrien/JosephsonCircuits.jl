# The scattering blocks on a backend: their stamps into the batched system
# matrices and the evaluation of their tabulated or entry-form data on the
# device, one work item per contribution and frequency.

# ---------------------------------------------------------------------------
# the scattering block stamps on a backend
# ---------------------------------------------------------------------------

# add one value per (contribution, frequency of the batch) into the stored
# values of that frequency's system matrix
@kernel function scatteringstampkernel!(nzval, @Const(index), @Const(values),
        ncontrib, nzA)
    gid = @index(Global)
    @inbounds begin
        q = gid - 1
        c = q % ncontrib + 1
        f = q ÷ ncontrib + 1
        nzval[Int(index[c]) + (f-1)*nzA] += values[c, f]
    end
end

"""
    sweepdestinations(A::SparseMatrixCSC, Aindex, adjoint::Bool)

The destination of each scattering contribution in the stored order a
[`FrequencySweepPlan`](@ref) uses, given its destinations in the stored order
of `A`.

A plan built for the transpose stores its values in `A`'s own order, so the
indices pass through; one built for the matrix itself stores them permuted
into compressed sparse row order, so they must be permuted the same way.
"""
function sweepdestinations(A::SparseMatrixCSC, Aindex::AbstractVector,
    adjoint::Bool)
    adjoint && return collect(Int, Aindex)
    # cscvaluepermutation carries a compressed sparse row slot to the stored
    # entry it holds; this is the other direction
    tocsr = invperm(cscvaluepermutation(A))
    return Int[tocsr[j] for j in Aindex]
end

"""
    DeviceScatteringStamps

The scattering block contribution to a batch of system matrices on a backend.

[`assemblesweep!`](@ref) covers the terms which are a constant quadratic in
the signal frequency, which the scattering blocks are not: each is evaluated
through its own provider, an arbitrary callable or an interpolation of
tabulated data, so its values have to be computed on the host. What is left
is a gather-add into the stored entries, which is this.

The values serve both directions of a sweep: a block's contribution does not
depend on the direction the system is assembled in, only on where in that
system's stored order each scalar lands. So they are computed and sent once,
and the adjoint direction gets a second view of them with its own
destinations (see [`transposedestinations`](@ref)).

The values are computed on the host, because a callable provider is an
arbitrary Julia function. A block whose data is tabulated or constant needs
no host at all; see [`plandeviceproviders`](@ref).
"""
struct DeviceScatteringStamps{VI,M,H,B}
    ssys::Any
    index::VI
    values::M               # ncontrib by nbatch, on the backend
    host::H
    wmodes::Vector{Float64}
    work::Any               # a ScatteringWorkspace, reused across batches
    backend::B
    nzA::Int
end

"""
    plandevicescattering(ssys, Aindexcsr, nzA, nbatch, backend, Nmodes)

Build a [`DeviceScatteringStamps`](@ref) from a scattering stamp system whose
destination indices have already been mapped into the stored order the batch
uses, or `nothing` when `ssys` is `nothing`.
"""
function plandevicescattering(ssys, Aindexcsr::AbstractVector{Int},
    nzA::Integer, nbatch::Integer, backend, Nmodes::Integer)

    isnothing(ssys) && return nothing
    n = length(Aindexcsr)
    Ti = nzA*nbatch < typemax(Int32) ? Int32 : Int
    values = KernelAbstractions.allocate(backend, Complex{Float64}, n, nbatch)
    host = Matrix{Complex{Float64}}(undef, n, nbatch)
    return DeviceScatteringStamps(ssys,
        tobackend(backend, convert(Vector{Ti}, Aindexcsr)), values, host,
        zeros(Float64, Nmodes), ScatteringWorkspace(), backend, Int(nzA))
end

"""
    transposedestinations(st::DeviceScatteringStamps, Aindexcsr, backend)

A view of the same scattering stamp values with different destinations, for
the transposed (adjoint) system: the contribution of a block does not depend
on the direction the system is assembled in, only on where in the stored
order of that system each scalar lands. Sharing `values` is what lets the
providers be evaluated once for both directions.
"""
function transposedestinations(st::DeviceScatteringStamps,
    Aindexcsr::AbstractVector{Int}, backend)
    Ti = eltype(st.index)
    return DeviceScatteringStamps(st.ssys,
        tobackend(backend, convert(Vector{Ti}, Aindexcsr)), st.values,
        st.host, st.wmodes, st.work, backend, st.nzA)
end

# compute the values of the batch beginning at `lo` into host buffer `buf`.
# Pure host work, and the only thing which touches the providers.
function computescatteringvalues!(st::DeviceScatteringStamps, buf, w,
    lo::Integer, k::Integer, wpumpmodes)

    nbatch = size(buf, 2)
    wmodes = st.wmodes
    @inbounds for j in 1:k
        for m in eachindex(wmodes)
            wmodes[m] = w[lo + j - 1] + wpumpmodes[m]
        end
        scatteringvalues!(view(buf, :, j), st.ssys, wmodes, st.work)
    end
    # a short final batch: repeat the last real frequency, as the assembly
    # does, so nothing reads uninitialized memory
    @inbounds for j in (k+1):nbatch
        copyto!(view(buf, :, j), view(buf, :, k))
    end
    return buf
end

"""
    stagescatteringstamps!(st::DeviceScatteringStamps, w, lo, k, wpumpmodes)

Compute the scattering values of the batch beginning at `lo` on the host and
send them to the backend.
"""
function stagescatteringstamps!(st::DeviceScatteringStamps, w, lo::Integer,
    k::Integer, wpumpmodes)
    computescatteringvalues!(st, st.host, w, lo, k, wpumpmodes)
    copyto!(st.values, st.host)
    return st
end

"""
    applyscatteringstamps!(nzval::AbstractMatrix, st::DeviceScatteringStamps)

Add the staged scattering values into the stored values `nzval` of a batch of
system matrices. [`stagescatteringstamps!`](@ref) must have run for the batch
`nzval` was assembled at.
"""
function applyscatteringstamps!(nzval::AbstractMatrix,
    st::DeviceScatteringStamps)
    scatteringstampkernel!(st.backend, 64)(nzval, st.index, st.values,
        size(st.values, 1), st.nzA; ndrange = length(st.values))
    KernelAbstractions.synchronize(st.backend)
    return nzval
end

# ---------------------------------------------------------------------------
# evaluating the scattering blocks on the backend
# ---------------------------------------------------------------------------

# the interpolated scattering parameter S[p,q] of one block at frequency w,
# from its table. A block whose data is a single matrix is a one point table,
# so it needs no separate path.
@inline function tableentry(freqs, vals, foff, nf, voff, n, p, q, w,
        extrapcode)
    @inbounds begin
        s = n*n
        if nf == 1
            return vals[voff + (q-1)*n + p]
        end
        f1 = freqs[foff+1]
        fn = freqs[foff+nf]
        # which segment, by the same branchless search the assembly kernels
        # use over a column pointer
        lo = 1; hi = nf - 1
        if w <= f1
            lo = 1
        elseif w >= fn
            lo = nf - 1
        else
            while lo < hi
                mid = (lo + hi + 1) >>> 1
                if freqs[foff+mid] <= w
                    lo = mid
                else
                    hi = mid - 1
                end
            end
        end
        a = freqs[foff+lo]
        b = freqs[foff+lo+1]
        t = (w - a)/(b - a)
        # `:constant` holds the end value beyond the ends; `:linear` lets the
        # end segment continue, which is what an unclamped t already does
        if extrapcode == Int8(1)
            t = t < 0 ? zero(t) : (t > 1 ? one(t) : t)
        end
        k1 = voff + (lo-1)*s + (q-1)*n + p
        return vals[k1]*(1-t) + vals[k1+s]*t
    end
end

# the value one contribution adds, from the scattering entry it reads. The
# two kernels below differ only in where that entry comes from: a table or a
# callable.
@inline function hybridcontribution(::Type{T}, S, wm, samepq::Bool, zrefq,
        coeffc, sgnc, scale) where {T}
    Bpq = zero(T); Cpq = zero(T)
    if iszero(wm)
        # the zero frequency rows are i = 0
        Cpq = samepq ? one(T) : zero(T)
    else
        # the impedance of the column port q, whose voltage or current the
        # entry multiplies
        r2 = sqrt(zrefq)
        rinv2 = one(r2)/r2
        Bpq = -rinv2*S
        Cpq = r2*S
        if samepq
            Bpq += rinv2
            Cpq += r2
        end
    end
    return coeffc == 1 ? sgnc*(im*wm*scale)*Bpq : -Cpq
end

# one work item per (contribution, frequency of the batch): interpolate the
# block, form the one hybrid coefficient the contribution needs, and write the
# value the gather-add will add
@kernel function deviceproviderkernel!(values, @Const(modeindex),
        @Const(blockindex), @Const(pindex), @Const(qindex), @Const(coeff),
        @Const(sgn), @Const(nports), @Const(zrefoff), @Const(zref),
        @Const(freqoff), @Const(nfreq), @Const(freqs), @Const(valoff),
        @Const(vals), @Const(conjsym), @Const(extrapcode), @Const(wpump),
        @Const(ws), scale, ncontrib)
    gid = @index(Global)
    @inbounds begin
        g = gid - 1
        c = g % ncontrib + 1
        f = g ÷ ncontrib + 1
        m = Int(modeindex[c]); bi = Int(blockindex[c])
        p = Int(pindex[c]); q = Int(qindex[c])
        wm = ws[f] + wpump[m]
        n = Int(nports[bi])
        T = eltype(values)
        S = zero(T)
        if !iszero(wm)
            isconj = conjsym[bi] != 0
            wq = isconj ? abs(wm) : wm
            S = tableentry(freqs, vals, Int(freqoff[bi]), Int(nfreq[bi]),
                Int(valoff[bi]), n, p, q, wq, extrapcode[bi])
            if isconj && wm < 0
                S = conj(S)
            end
        end
        values[c, f] = hybridcontribution(T, S, wm, p == q,
            zref[Int(zrefoff[bi]) + q], coeff[c], sgn[c], scale)
    end
end

# the same, for blocks whose scattering parameters come from a callable of the
# `:entry` form. The callables live in a device array indexed by block: they
# capture only numbers, so they are `isbits`, and blocks built from one helper
# share a type, which is what lets them be stored and called this way.
@kernel function deviceentrykernel!(values, @Const(modeindex),
        @Const(blockindex), @Const(pindex), @Const(qindex), @Const(coeff),
        @Const(sgn), @Const(zrefoff), @Const(zref), @Const(funcs),
        @Const(conjsym), @Const(wpump), @Const(ws), scale, ncontrib)
    gid = @index(Global)
    @inbounds begin
        g = gid - 1
        c = g % ncontrib + 1
        f = g ÷ ncontrib + 1
        m = Int(modeindex[c]); bi = Int(blockindex[c])
        p = Int(pindex[c]); q = Int(qindex[c])
        wm = ws[f] + wpump[m]
        T = eltype(values)
        S = zero(T)
        if !iszero(wm)
            isconj = conjsym[bi] != 0
            wq = isconj ? abs(wm) : wm
            S = T(funcs[bi](p, q, wq))
            if isconj && wm < 0
                S = conj(S)
            end
        end
        values[c, f] = hybridcontribution(T, S, wm, p == q,
            zref[Int(zrefoff[bi]) + q], coeff[c], sgn[c], scale)
    end
end

"""
    DeviceProviders

The scattering blocks of a stamp system on a backend, as flat tables.

The values a batch of system matrices needs are otherwise computed on the
host, because a callable provider is an arbitrary Julia function. A block
whose data is tabulated or constant is not: its evaluation is a search and an
interpolation, which is the shape a kernel wants. When every block of a
circuit is one of those, the whole evaluation moves to the backend and the
host does nothing per frequency at all.

A constant block is stored as a one point table, so it needs no separate path.
Which frequency a block is evaluated at still depends on its negative
frequency rule, and a range which must not be extrapolated is checked on the
host before each batch, because a kernel cannot raise.
"""
struct DeviceProviders{VI,VZ,VR,VC,VF,B}
    nports::VI
    zrefoff::VI
    zref::VR
    freqoff::VI
    nfreq::VI
    freqs::VR
    valoff::VI
    vals::VC
    # the callables of the `:entry` form, or `nothing` when the blocks are
    # tabulated; exactly one of `vals` and `funcs` is used
    funcs::VF
    conjsym::VZ
    extrapcode::VZ
    modeindex::VI
    blockindex::VI
    pindex::VI
    qindex::VI
    coeff::VZ
    sgn::VZ
    wpump::VR
    ws::VR
    wshost::Vector{Float64}
    scale::Float64
    ncontrib::Int
    ranges::Vector{Tuple{Float64,Float64}}
    strict::Vector{Bool}
    conjhost::Vector{Bool}
    names::Vector{String}
    ssys::Any
    backend::B
end

"""
    candeviceevaluate(ssys)

Whether every block of a stamp system has data a kernel can evaluate, which
is tabulated or constant. A callable provider is an arbitrary Julia function
and stays on the host.
"""
function candeviceevaluate(ssys)
    isnothing(ssys) && return false
    isempty(ssys.blocks) && return false
    first = ssys.blocks[1].block.provider
    if first isa CallableMatrixProvider
        # a callable can be called from a kernel only in the `:entry` form,
        # and only if the closures can live in a device array, which needs
        # them to be `isbits` and all of one type
        first.form === :entry || return false
        isbits(first.f) || return false
        for sb in ssys.blocks
            p = sb.block.provider
            p isa CallableMatrixProvider || return false
            p.form === :entry || return false
            typeof(p.f) === typeof(first.f) || return false
        end
        return true
    end
    for sb in ssys.blocks
        p = sb.block.provider
        (p isa TabulatedMatrixProvider || p isa ConstantMatrixProvider) ||
            return false
        p isa TabulatedMatrixProvider && p.interpolation != :linear &&
            return false
    end
    return true
end

# the callables of an entry-wise stamp system, or `nothing` when its blocks
# are tabulated
function entrycallables(ssys)
    p1 = ssys.blocks[1].block.provider
    (p1 isa CallableMatrixProvider && p1.form === :entry) || return nothing
    return [sb.block.provider.f for sb in ssys.blocks]
end

# The stamp system holds its blocks untyped, so a loop which reads a
# provider's arrays element by element dispatches on every element. These
# barriers give the compiler the concrete types once.
function copyfrequencies!(dest::Vector{Float64}, off::Int,
    src::AbstractVector)
    @inbounds for k in eachindex(src)
        dest[off + k] = Float64(src[k])
    end
    return nothing
end

function copytable!(dest::Vector{Complex{Float64}}, off::Int,
    src::AbstractArray{<:Number,3}, n::Int, nf::Int)
    @inbounds for k in 1:nf, j in 1:n, i in 1:n
        dest[off + (k-1)*n*n + (j-1)*n + i] = src[i, j, k]
    end
    return nothing
end

function copymatrix!(dest::Vector{Complex{Float64}}, off::Int,
    src::AbstractMatrix{<:Number}, n::Int)
    @inbounds for j in 1:n, i in 1:n
        dest[off + (j-1)*n + i] = src[i, j]
    end
    return nothing
end

"""
    plandeviceproviders(ssys, nbatch, backend, wpumpmodes, scale)

Build a [`DeviceProviders`](@ref) for a stamp system whose blocks can all be
evaluated by a kernel, or `nothing` otherwise.
"""
function plandeviceproviders(ssys, nbatch::Integer, backend, wpumpmodes,
    scale::Real)

    candeviceevaluate(ssys) || return nothing
    nb = length(ssys.blocks)
    nports = Int32[sb.block.nports for sb in ssys.blocks]
    callables = entrycallables(ssys)
    istable = isnothing(callables)
    # size everything first and fill it in place: a line whose every cell is
    # its own block has hundreds of thousands of table points, and growing
    # the flat arrays a block at a time would be slow
    ntot = 0; vtot = 0; ztot = 0
    for sb in ssys.blocks
        p = sb.block.provider
        n = sb.block.nports
        ztot += n
        istable || continue
        if p isa ConstantMatrixProvider
            ntot += 1; vtot += n*n
        else
            nf = length(p.frequencies)
            ntot += nf; vtot += n*n*nf
        end
    end
    zrefoff = Vector{Int32}(undef, nb); zref = Vector{Float64}(undef, ztot)
    freqoff = Vector{Int32}(undef, nb); nfreq = Vector{Int32}(undef, nb)
    freqs = Vector{Float64}(undef, ntot)
    valoff = Vector{Int32}(undef, nb)
    vals = Vector{Complex{Float64}}(undef, vtot)
    conjsym = Vector{Int8}(undef, nb); extrapcode = Vector{Int8}(undef, nb)
    ranges = Vector{Tuple{Float64,Float64}}(undef, nb)
    strict = Vector{Bool}(undef, nb)
    conjhost = Vector{Bool}(undef, nb)
    names = String[sb.name for sb in ssys.blocks]
    zi = 0; fi = 0; vi = 0
    for (bi, sb) in enumerate(ssys.blocks)
        blk = sb.block
        n = blk.nports
        zrefoff[bi] = zi
        @inbounds for p in 1:n
            zref[zi + p] = real(blk.zref[p])
        end
        zi += n
        prov = blk.provider
        freqoff[bi] = fi
        valoff[bi] = vi
        if !istable
            # an entry-wise callable has no table, no range, and nothing to
            # extrapolate
            nfreq[bi] = 0
            extrapcode[bi] = Int8(1)
            ranges[bi] = (-Inf, Inf)
            strict[bi] = false
        elseif prov isa ConstantMatrixProvider
            # a single matrix is a one point table
            nfreq[bi] = 1
            freqs[fi + 1] = 0.0
            fi += 1
            copymatrix!(vals, vi, prov.A, n)
            vi += n*n
            extrapcode[bi] = Int8(1)
            ranges[bi] = (0.0, 0.0)
            strict[bi] = false
        else
            nf = length(prov.frequencies)
            nfreq[bi] = nf
            copyfrequencies!(freqs, fi, prov.frequencies)
            fi += nf
            copytable!(vals, vi, prov.values, n, nf)
            vi += n*n*nf
            extrapcode[bi] = prov.extrapolation == :linear ? Int8(2) : Int8(1)
            ranges[bi] = (Float64(prov.frequencies[1]),
                Float64(prov.frequencies[nf]))
            strict[bi] = prov.extrapolation == :error && nf > 1
        end
        conjhost[bi] = !(blk.negative_frequency isa Native)
        conjsym[bi] = conjhost[bi] ? Int8(1) : Int8(0)
    end

    d(x) = tobackend(backend, x)
    funcs = istable ? nothing : tobackend(backend, callables)
    ncontrib = length(ssys.Aindex)
    return DeviceProviders(d(nports), d(zrefoff), d(zref), d(freqoff),
        d(nfreq), d(freqs), d(valoff), d(vals), funcs, d(conjsym),
        d(extrapcode),
        d(Int32.(ssys.modeindex)), d(Int32.(ssys.blockindex)),
        d(Int32.(ssys.pindex)), d(Int32.(ssys.qindex)), d(Int8.(ssys.coeff)),
        d(Int8.(ssys.sign)), d(collect(Float64, wpumpmodes)),
        d(zeros(Float64, nbatch)), zeros(Float64, nbatch), Float64(scale),
        ncontrib, ranges, strict, conjhost, names, ssys, backend)
end

# a kernel cannot raise, so a block whose data must not be extrapolated has
# its range checked here, before the batch runs
function checkdeviceranges(dp::DeviceProviders, k::Integer)
    any(dp.strict) || return nothing
    lonat = Inf; hinat = -Inf; loabs = Inf; hiabs = -Inf
    nm = length(dp.wpump)
    wp = Array(dp.wpump)
    @inbounds for j in 1:k
        for m in 1:nm
            w = dp.wshost[j] + wp[m]
            iszero(w) && continue
            lonat = min(lonat, w); hinat = max(hinat, w)
            a = abs(w)
            loabs = min(loabs, a); hiabs = max(hiabs, a)
        end
    end
    for bi in eachindex(dp.strict)
        dp.strict[bi] || continue
        lo, hi = dp.ranges[bi]
        l, h = dp.conjhost[bi] ? (loabs, hiabs) : (lonat, hinat)
        if l < lo || h > hi
            throw(ArgumentError(lazy"The scattering block at $(dp.names[bi]) is evaluated over [$(l), $(h)] rad/s but its tabulated range is [$(lo), $(hi)] rad/s. Extrapolation of tabulated data is opt-in: pass extrapolation = :constant or :linear if extrapolation is intended."))
        end
    end
    return nothing
end

"""
    stagedeviceproviders!(values, dp::DeviceProviders, w, lo, k)

Compute the scattering values of the batch of frequencies beginning at `lo`
directly on the backend, into `values`.

The host part of this is the batch's own signal frequencies and, for a block
whose data must not be extrapolated, a check of the range it will be
evaluated over. Nothing per block or per contribution crosses to the host.
"""
function stagedeviceproviders!(values::AbstractMatrix, dp::DeviceProviders,
    w, lo::Integer, k::Integer)

    nbatch = length(dp.wshost)
    @inbounds for j in 1:nbatch
        dp.wshost[j] = w[min(lo + j - 1, lo + k - 1)]
    end
    checkdeviceranges(dp, k)
    copyto!(dp.ws, dp.wshost)
    if isnothing(dp.funcs)
        deviceproviderkernel!(dp.backend, 64)(values, dp.modeindex,
            dp.blockindex, dp.pindex, dp.qindex, dp.coeff, dp.sgn, dp.nports,
            dp.zrefoff, dp.zref, dp.freqoff, dp.nfreq, dp.freqs, dp.valoff,
            dp.vals, dp.conjsym, dp.extrapcode, dp.wpump, dp.ws, dp.scale,
            dp.ncontrib; ndrange = length(values))
    else
        deviceentrykernel!(dp.backend, 64)(values, dp.modeindex,
            dp.blockindex, dp.pindex, dp.qindex, dp.coeff, dp.sgn,
            dp.zrefoff, dp.zref, dp.funcs, dp.conjsym, dp.wpump, dp.ws,
            dp.scale, dp.ncontrib; ndrange = length(values))
    end
    KernelAbstractions.synchronize(dp.backend)
    return values
end
