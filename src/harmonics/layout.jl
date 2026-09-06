

"""
    ModeLayout(isreal::AbstractVector{Bool}, dim::Integer, ::Type{Ti}=Int)
    ModeLayout(realindices, nmodes::Integer, dim::Integer, ::Type{Ti}=Int)

Layout of one axis of length `dim` (a complex dimension), built from a length-
`nmodes` mask of which modes are real. `dim` must be an integer multiple of
`nmodes`. Complex index `i` owns real slots `ptr[i]:ptr[i+1]-1`; `inv` maps a
real slot back to its complex index; `isfirst` marks the slots that start a
mode; `rdim` is the resulting real dimension.

`w` is a bit per index recording whether that mode is real.

# Fields
- `nmodes`, `dim`, `rdim`: the mode count, the complex dimension and the real
    dimension.
- `nreal`, `isreal`: the number of real modes and the mask over the modes.
- `ptr`, `inv`, `isfirst`: the slot ranges, the inverse map and the mode
    starts, as above.
- `w`: the bit per complex index, whether its mode is real.
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

"""
    realblockterm(v, dr, dc)

The entry at offset `(dr, dc)` of the real block a complex value `v`
contributes: the two by two block `[real(v) -imag(v); imag(v) real(v)]`
which multiplication by `v` is in the `(real, imag)` coordinates of a
complex mode, of which a self conjugate (real) mode keeps only the first
row, the first column, or both. The one rule behind every real form in the
package: the sparse [`complex_to_real`](@ref), the real Jacobian assembly
and the block factorization values.
"""
@inline function realblockterm(v, dr, dc)
    if dc == 0
        return dr == 0 ? real(v) : imag(v)
    else
        return dr == 0 ? -imag(v) : real(v)
    end
end
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

# The conversions between the complex node vector and its real
# representation, in the forms the solvers use: a vector either way, and a
# sparse matrix into its real form under a row and a column layout. Every
# self conjugate mode contributes one real entry and every conjugate pair
# two, `(real, imag)`, so a complex entry `a` of a matrix expands to the
# block `[re -im; im re]` where both modes are pairs and to a single row or
# column of it where one is self conjugate. The wider family (dense forms,
# the mask entry points, a conjugated input, a scaling of the real modes)
# lives in the tests as the reference these are checked against.

#  vectors

"""
    complex_to_real!(xr, xc, isreal) -> xr

Complex -> real. The imaginary part of a real mode is never read.
"""
function complex_to_real!(xr::AbstractVector{T}, xc::AbstractVector{Complex{T}},
               isreal::AbstractVector{Bool}) where {T<:Real}
    length(xr) == realdim(length(xc), isreal) || throw(DimensionMismatch(
        "destination has length $(length(xr)), expected $(realdim(length(xc), isreal))"))
    nm = length(isreal)
    p, t = 1, 1
    @inbounds for i in eachindex(xc)
        a = xc[i]
        if isreal[t]
            xr[p] = real(a)
            p += 1
        else
            xr[p]   = real(a)
            xr[p+1] = imag(a)
            p += 2
        end
        t = _next(t, nm)
    end
    return xr
end

"""
    real_to_complex!(xc, xr, isreal) -> xc

Real -> complex. The imaginary part of a real mode is written as an explicit
zero, so the result does not depend on what was in `xc` beforehand.
"""
function real_to_complex!(xc::AbstractVector{Complex{T}}, xr::AbstractVector{T},
                 isreal::AbstractVector{Bool}) where {T<:Real}
    length(xr) == realdim(length(xc), isreal) || throw(DimensionMismatch(
        "source has length $(length(xr)), expected $(realdim(length(xc), isreal))"))
    nm = length(isreal)
    p, t = 1, 1
    @inbounds for i in eachindex(xc)
        if isreal[t]
            xc[i] = Complex(xr[p], zero(T))
            p += 1
        else
            xc[i] = Complex(xr[p], xr[p+1])
            p += 2
        end
        t = _next(t, nm)
    end
    return xc
end

# the allocating vector forms of the conversions above, for a mask
complex_to_real(xc::AbstractVector{Complex{T}}, isreal::AbstractVector{Bool}) where {T} =
    complex_to_real!(Vector{T}(undef, realdim(length(xc), isreal)), xc, isreal)
real_to_complex(xr::AbstractVector{T}, isreal::AbstractVector{Bool}) where {T} =
    real_to_complex!(Vector{Complex{T}}(undef, complexdim(length(xr), isreal)), xr, isreal)

#  sparse matrices

"""
    complex_to_real(A, rowlayout, collayout, ::Type{Tj}=Ti) -> SparseMatrixCSC{T,Tj}

Real form of `x -> A*x` under the row layout of its output and the column
layout of its input, so that `complex_to_real(A, rl, cl) * complex_to_real(x, cl.isreal)
== complex_to_real(A*x, rl.isreal)`. `Tj` selects the index type of the
result; `Int32` halves `rowval` and is worth using whenever the dimensions
fit. Runs one O(nnz) pass to size the output, then one to fill it.
"""
function complex_to_real(A::SparseMatrixCSC{Complex{T},Ti}, rl::ModeLayout, cl::ModeLayout,
                 ::Type{Tj} = Ti) where {T<:Real,Ti,Tj<:Integer}
    _checkdims(A, rl, cl)
    n = size(A, 2)
    Ap, Ai, Av = SparseArrays.getcolptr(A), rowvals(A), nonzeros(A)
    rptr, cptr = rl.ptr, cl.ptr

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
        k0 = k
        for idx in Ap[j]:Ap[j+1]-1
            i  = Ai[idx]
            r0 = rptr[i]
            wr = rptr[i+1] - r0
            a  = Av[idx]
            rowval[k]      = r0
            rowval[k+wr-1] = ifelse(wr == 2, r0 + one(r0), r0)
            nzval[k]       = realblockterm(a, 0, 0)
            nzval[k+wr-1]  = realblockterm(a, wr - 1, 0)
            k += wr
        end
        S = k - k0
        colptr[c0+1] = k
        if wc == 2
            # the second column of each block, with the rows of the first
            copyto!(rowval, k, rowval, k0, S)
            for idx in Ap[j]:Ap[j+1]-1
                i  = Ai[idx]
                wr = rptr[i+1] - rptr[i]
                d  = Av[idx]
                nzval[k]      = realblockterm(d, 0, 1)
                nzval[k+wr-1] = realblockterm(d, wr - 1, 1)
                k += wr
            end
            colptr[c0+2] = k
        end
    end
    return SparseMatrixCSC{T,Tj}(rl.rdim, cl.rdim, colptr, rowval, nzval)
end

@inline function _checkdims(A, rl::ModeLayout, cl::ModeLayout)
    rl.dim == size(A, 1) || throw(DimensionMismatch(
        "row layout is for dimension $(rl.dim), matrix has $(size(A,1))"))
    cl.dim == size(A, 2) || throw(DimensionMismatch(
        "column layout is for dimension $(cl.dim), matrix has $(size(A,2))"))
    nothing
end

# =====================================================================
# The canonical state layout `[internal | vdc]` and the workspaces of a
# canonical evaluation, which append the direct current unknowns to the
# internal real state described above.

# The canonical state layout: the internal real state as the solver
# evaluates it, followed by the explicit direct current coordinates.
#
# The internal real state is node major and mode minor: every node
# contributes one entry per self conjugate mode and two per conjugate pair,
# in mode order, so the zero frequency entry of a node sits between that
# node's alternating current entries. That is the layout `ModeLayout`
# describes and every kernel in the solve is written against.
#
# The canonical state is that state, untouched, with the explicit average
# voltages appended:
#
#     [ internal | vdc ]
#
# The direct current block -- the zero frequency flux of each node and of
# each auxiliary unknown, and the explicit voltages -- is then not
# contiguous: its flux entries sit where the internal layout put them and
# are named by an index, `dcpos`. That is the whole cost of leaving the
# internal state alone, and it is paid on a window of a few hundred entries
# rather than on the state. The residual, the Jacobian vector product and
# the preconditioner run on the internal part of the canonical vector
# directly, through a view, and only the window is gathered, worked on and
# scattered back.
#
# The earlier layout grouped the state by role, `[phiac | phidc | vdc]`, to
# make the block contiguous, and bracketed every residual, product and
# preconditioner application with a permutation of the whole state to get
# there. Measured on a 2048 cell line with a direct current bias those
# passes were under one percent of the solve on the host, so this is not a
# change made for time there; the copies were overhead with nothing behind
# them, and on a device each was a kernel launch and a synchronization.

"""
    CompositeLayout

The canonical state layout, `[internal | vdc]`, and where the zero frequency
entries of the internal state are.

# Fields
- `rdim`: the length of the internal real state, which is the first block.
- `ndc`: the number of zero frequency entries in it: one per node and one
  per auxiliary unknown, in internal order.
- `nvdc`: the length of the explicit voltage block, which follows. Zero
  until a circuit injects direct current.
- `dcpos`: for the `k`th zero frequency entry, its position in the internal
  state.

The direct current window is the `ndc + nvdc` entries the block reads and
writes: the zero frequency entries first, then the voltages. Window entry
`k` sits at canonical position [`windowindex`](@ref)`(L, k)`.

See [`compositelayout`](@ref).
"""
struct CompositeLayout
    rdim::Int
    ndc::Int
    nvdc::Int
    dcpos::Vector{Int}
end

"""
    canonicaldim(L::CompositeLayout)

The length of the canonical state.
"""
canonicaldim(L::CompositeLayout) = L.rdim + L.nvdc

"""
    nwindow(L::CompositeLayout)

The length of the direct current window.
"""
nwindow(L::CompositeLayout) = L.ndc + L.nvdc

"""
    isinternal(L::CompositeLayout)

Whether the canonical state is the internal one, which it is exactly while
there are no explicit direct current coordinates.
"""
isinternal(L::CompositeLayout) = iszero(L.nvdc)

"""
    windowindex(L::CompositeLayout, k)

The canonical position of window entry `k`: a zero frequency entry of the
internal state for `k <= L.ndc`, an explicit voltage after that.
"""
windowindex(L::CompositeLayout, k::Integer) =
    k <= L.ndc ? L.dcpos[k] : L.rdim + (k - L.ndc)

"""
    windowindices(L::CompositeLayout)

The canonical positions of the whole window, in window order.
"""
windowindices(L::CompositeLayout) =
    vcat(L.dcpos, collect((L.rdim + 1):(L.rdim + L.nvdc)))

# the canonical index range of the explicit voltage block
voltagerange(L::CompositeLayout) = (L.rdim + 1):(L.rdim + L.nvdc)

"""
    internalpart(u::AbstractVector, L::CompositeLayout)

The internal block of a canonical vector `u`, as a view: the first
`L.rdim` entries, which the harmonic system reads and writes in place.
"""
internalpart(u::AbstractVector, L::CompositeLayout) = view(u, 1:L.rdim)

"""
    compositelayout(ml::ModeLayout, isdc::AbstractVector{Bool};
        nvdc = 0)

Build the canonical layout for a state whose internal layout is `ml`, where
`isdc[t]` marks mode `t` as the zero frequency one.

A zero frequency mode is self conjugate, so it occupies one internal entry
rather than two; `isdc` must therefore imply `ml.isreal`. The converse does
not hold, since a Nyquist mode is self conjugate without being direct
current, which is why the split is made on `isdc` and not on `ml.isreal`.
"""
function compositelayout(ml::ModeLayout, isdc::AbstractVector{Bool};
        nvdc::Int = 0)
    length(isdc) == ml.nmodes || throw(DimensionMismatch(
        lazy"`isdc` has length $(length(isdc)) but the layout has $(ml.nmodes) modes."))
    for t in eachindex(isdc)
        if isdc[t] && !ml.isreal[t]
            throw(ArgumentError(lazy"mode $(t) is marked zero frequency but is not self conjugate; a zero frequency component of a real signal has no imaginary part."))
        end
    end
    nvdc >= 0 || throw(ArgumentError("the voltage block length must be nonnegative."))

    nmodes = ml.nmodes
    nper = ml.dim ÷ nmodes            # entries per mode across the state
    ndc = nper * count(isdc)

    dcpos = Vector{Int}(undef, ndc)
    d = 1                             # cursor into the window
    p, t = 1, 1                       # cursor into the internal real state
    @inbounds for _ in 1:ml.dim
        w = ml.isreal[t] ? 1 : 2
        if isdc[t]
            for k in 0:w-1
                dcpos[d] = p + k
                d += 1
            end
        end
        p += w
        t = _next(t, nmodes)
    end
    d == ndc + 1 || error("internal error: found $(d-1) of $(ndc) direct current entries.")

    return CompositeLayout(ml.rdim, ndc, nvdc, dcpos)
end

"""
    compositelayout(ml::ModeLayout, modes::AbstractVector{<:Tuple}; kwargs...)

Build the canonical layout from the mode index tuples, taking the zero
frequency mode to be the one whose indices are all zero.
"""
function compositelayout(ml::ModeLayout, modes::AbstractVector{<:Tuple};
        kwargs...)
    return compositelayout(ml, [all(iszero, m) for m in modes]; kwargs...)
end

"""
    gathercanonical!(u, rint, L::CompositeLayout)

Write the internal real state `rint` into the internal block of the
canonical state `u`, a copy.

The `vdc` block has no internal counterpart and is left untouched, so a
caller which keeps explicit direct current coordinates in `u` does not lose
them here. A solve does not call this: it hands the harmonic system the
internal block of `u` through [`internalpart`](@ref) and no copy is made.
The interfaces which want a state of their own do.
"""
function gathercanonical!(u::AbstractVector, rint::AbstractVector,
        L::CompositeLayout)
    length(u) == canonicaldim(L) || throw(DimensionMismatch(
        lazy"the canonical state has length $(length(u)) but the layout needs $(canonicaldim(L))."))
    length(rint) == L.rdim || throw(DimensionMismatch(
        lazy"the internal state has length $(length(rint)) but the layout needs $(L.rdim)."))
    copyto!(u, 1, rint, 1, L.rdim)
    return u
end

"""
    scattercanonical!(rint, u, L::CompositeLayout)

Write the internal block of the canonical state `u` into `rint`, a copy.
The inverse of [`gathercanonical!`](@ref) on that block.
"""
function scattercanonical!(rint::AbstractVector, u::AbstractVector,
        L::CompositeLayout)
    length(u) == canonicaldim(L) || throw(DimensionMismatch(
        lazy"the canonical state has length $(length(u)) but the layout needs $(canonicaldim(L))."))
    length(rint) == L.rdim || throw(DimensionMismatch(
        lazy"the internal state has length $(length(rint)) but the layout needs $(L.rdim)."))
    copyto!(rint, 1, u, 1, L.rdim)
    return rint
end

# The window is gathered and scattered by index: on the host a plain loop,
# on a device the kernel. It is a few hundred entries, so what this costs is
# the launch and not the copy.
_onhost(x) = KernelAbstractions.get_backend(x) isa CPU

"""
    gathervalues!(dest::AbstractArray, src::AbstractVector,
        index::AbstractArray)

`dest[k] = src[index[k]]` for every `k`, as a KernelAbstractions kernel on the
backend of `dest`, which `src` and `index` must share. The device side of
the permuted copies of the canonical layout (`_gatherperm!`): `index` is an
injection, so no element is written twice and no atomic is needed.
"""
function gathervalues!(dest::AbstractArray, src::AbstractVector,
    index::AbstractArray)
    size(dest) == size(index) || throw(DimensionMismatch(
        lazy"`dest` is $(size(dest)) but the index is $(size(index))."))
    backend = KernelAbstractions.get_backend(dest)
    kernel! = gatherkernel!(backend)
    kernel!(dest, src, index; ndrange = length(index))
    KernelAbstractions.synchronize(backend)
    return dest
end

"""
    scattervalues!(dest::AbstractVector, src::AbstractVector,
        index::AbstractArray)

`dest[index[k]] = src[k]` for every `k`, as a KernelAbstractions kernel on
the backend of `src`. The inverse of [`gathervalues!`](@ref) when `index` is
a permutation; `index` must not repeat, or the writes race.
"""
function scattervalues!(dest::AbstractVector, src::AbstractVector,
        index::AbstractArray)
    length(src) == length(index) || throw(DimensionMismatch(
        lazy"`src` has length $(length(src)) but the index has $(length(index))."))
    backend = KernelAbstractions.get_backend(src)
    kernel! = scatterkernel!(backend)
    kernel!(dest, src, index; ndrange = length(index))
    KernelAbstractions.synchronize(backend)
    return dest
end

# one work item per entry of the index; an injective index means no entry
# is read or written twice and neither kernel needs an atomic
@kernel function gatherkernel!(dest, @Const(src), @Const(index))
    k = @index(Global)
    @inbounds dest[k] = src[index[k]]
end

@kernel function scatterkernel!(dest, @Const(src), @Const(index))
    k = @index(Global)
    @inbounds dest[index[k]] = src[k]
end


function _gatherperm!(dest, src, index)
    if _onhost(src)
        @inbounds for k in eachindex(index)
            dest[k] = src[index[k]]
        end
        return dest
    end
    gathervalues!(view(dest, 1:length(index)), src, index)
    return dest
end

function _scatterperm!(dest, src, index)
    if _onhost(dest)
        @inbounds for k in eachindex(index)
            dest[index[k]] = src[k]
        end
        return dest
    end
    scattervalues!(dest, view(src, 1:length(index)), index)
    return dest
end

# the layout is a host object: the index a device needs is carried by the
# work built on it, on the backend the state is on
tobackend(::Backend, L::CompositeLayout) = L

# The direct current block is held in `Float64` whatever precision the
# periodic solve runs in. It is small, dense and solved exactly, and its
# conditioning is the worst in the problem (scaled conductances of order
# 1e-11 against injected currents of order 1e8), so a rank decision or a
# factorization taken in a lower precision would be the least accurate
# part of the answer. A single precision solve still converges to single
# precision, and the block is a few hundred numbers beside a state of tens
# of thousands, so this costs nothing.

"""
    DCPinning

The references a singular direct current subsystem needs: which redundant
equations to give up, and which coordinate each one fixes at zero.

A direction the descriptor does not determine may be pinned only when
nothing outside the descriptor can see it. Write the rest of the harmonic
residual's dependence on the direct current unknowns as `H`, the zero
frequency nodal currents they drive; a null direction `N` is a gauge exactly
when `H N = 0`.

The common case is a floating island's average voltage. A conductance
island with no path to ground has zero row sums, so raising every voltage in
it by the same amount drives no current anywhere, and only differences were
ever physical. Choosing a reference there changes nothing.

An ideal through in parallel with an inductor is the opposite. Its free
direction is the division of current between the two ideal branches, which
cancels in the transport row because both terminals lie in one static flux
component, and does not cancel at the nodes: it injects `+d` at one and `-d`
at the other, moving the inductor current, the static flux across it and,
through a junction, the nonlinear operating point. `H N` is nonzero, and the
circuit is refused rather than given one of infinitely many answers.

A reference is written as `y_c = 0` for a chosen coordinate rather than as a
minimum norm condition on the whole direction. The subsystem mixes volts and
amperes, so a minimum norm row is not invariant under a change of units,
while fixing one coordinate is; and the row it produces has a single entry.

# Fields
- `rows`: positions in the subsystem whose equation is replaced.
- `cols`: the coordinate each replaced row fixes at zero.
"""
struct DCPinning
    rows::Vector{Int}      # positions in the subsystem to replace
    cols::Vector{Int}      # the coordinate each replaced row fixes at zero
end

"""
    CanonicalWork

The workspaces a canonical evaluation needs, and the direct current block
when it is explicit.

# Fields
- `layout`: the [`CompositeLayout`](@ref) the canonical state is written in.
- `xint`, `Fint`: an internal state and an internal residual, for the
    interfaces which need one of their own; the solve reads and writes the
    internal block of the canonical vector in place.
- `transport`: the transport rows, or `nothing` when there is no explicit
    block.
- `blockrows`: the scattering blocks' own zero frequency rows, or `nothing`.
- `dwork`: the resistor current the coupling drives into the nodes.
- `nnodaldc`: where the zero frequency entries split. The first are the zero
    frequency flux of each node, which the transport coupling drives; any
    after them belong to auxiliary branch currents, which it does not.
- `pinning`: the reference rows, or `nothing`. See [`DCPinning`](@ref).
- `dcindex`, `dclocal`: the canonical positions of the direct current
    subsystem's unknowns, and the same positions local to the window.
- `window`: the canonical positions of the window, on the backend the state
    is on, for the gather and the scatter.
- `Fwindow`, `uwindow`: the window itself, on that backend.
- `update`: the block in its matrix form, where the state lives, or
    `nothing` on the host and when there is no explicit block, where the
    scalar walk of `addtransportwindow!` is already the cheaper of the two.
"""
struct CanonicalWork{T,V<:AbstractVector{T},L<:CompositeLayout,TR,BR,I}
    layout::L
    xint::V
    Fint::V
    transport::TR
    blockrows::BR
    dwork::Vector{Float64}
    nnodaldc::Int
    pinning::Union{Nothing,DCPinning}
    dcindex::Vector{Int}
    dclocal::Vector{Int}
    window::I
    Fwindow::V
    uwindow::V
    update::Any
end

function CanonicalWork(L::CompositeLayout, proto::AbstractVector{T};
        transport = nothing, blockrows = nothing,
        nnodaldc::Int = L.ndc) where {T}
    nd = isnothing(transport) ? 0 : size(transport.coupling, 1)
    if !isnothing(transport)
        nvoltages(transport) == L.nvdc || throw(DimensionMismatch(
            lazy"the transport rows carry $(nvoltages(transport)) voltages but the layout has $(L.nvdc)."))
        nd <= nnodaldc || throw(DimensionMismatch(
            lazy"the coupling drives $(nd) nodal rows but only $(nnodaldc) of the direct current block are nodal."))
    end
    nw = nwindow(L)
    bk = KernelAbstractions.get_backend(proto)
    window = tobackend(bk, windowindices(L))
    w = CanonicalWork(L, similar(proto, L.rdim), similar(proto, L.rdim),
        transport, blockrows, zeros(Float64, nd), nnodaldc, nothing, Int[],
        Int[], window, similar(proto, nw), similar(proto, nw), nothing)
    isnothing(transport) && return w
    # the pinning is read off the subsystem this work describes, so it is
    # found once here and then carried
    full = CanonicalWork(L, w.xint, w.Fint, transport, blockrows, w.dwork,
        nnodaldc, dcpinning(w), dcsubsystemindices(w), dcsubsystemlocal(w),
        window, w.Fwindow, w.uwindow, nothing)
    # the matrix form, on the backend the state is on. On the host the
    # scalar walk is already the cheaper of the two at these sizes, so it is
    # built only where the state is not host resident.
    _onhost(proto) && return full
    up = dcupdate(full)
    isnothing(up) && return full
    dev = DCUpdate(tobackend(bk, up.keep), tobackend(bk, up.rowptr),
        tobackend(bk, up.colval), tobackend(bk, up.nzval),
        tobackend(bk, up.cresidual))
    return CanonicalWork(L, full.xint, full.Fint, transport, blockrows,
        full.dwork, nnodaldc, full.pinning, full.dcindex, full.dclocal,
        window, full.Fwindow, full.uwindow, dev)
end
