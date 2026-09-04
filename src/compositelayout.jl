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
    compositelayout(ml::ModeLayout, modes::AbstractVector; kwargs...)

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
- `update`: the block in its matrix form, where the state lives.
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

# The explicit direct current block's contribution: the resistor current the
# average voltages drive into the zero frequency nodal rows, and the
# transport rows themselves.
#
# The sign is the elimination's, read backwards. The residual is
# `A x - bnm`, and eliminating replaces `bnm` by `bnm - G0 P v`, so carrying
# `v` instead means adding `G0 P v` here and leaving the applied source
# alone. Doing both would count the resistor current twice.
#
# Everything the block reads and writes is the window: the zero frequency
# entries, scattered through the internal state, and the explicit voltages
# after it. The window is gathered by index into a buffer of its own, worked
# on with local indices, and scattered back, which is what lets the
# arithmetic run on a device without scalar indexing; it is a few hundred
# numbers where the state is tens of thousands.
function addtransport!(Fc::AbstractVector, work::CanonicalWork,
        u::AbstractVector; residual::Bool = true)
    isnothing(work.transport) && return Fc
    Fw, uw = work.Fwindow, work.uwindow
    _gatherperm!(Fw, Fc, work.window)
    _gatherperm!(uw, u, work.window)
    if isnothing(work.update)
        addtransportwindow!(Fw, uw, work; residual)
    else
        # in place, where the state is: three array operations and no copy
        applydcupdate!(Fw, uw, work.update; residual)
    end
    _scatterperm!(Fc, Fw, work.window)
    return Fc
end

# `Fw` and `uw` are the window; every index below is local to it.
function addtransportwindow!(Fw::AbstractVector, uw::AbstractVector,
        work::CanonicalWork; residual::Bool = true)
    t = work.transport
    L = work.layout
    vr = (L.ndc + 1):(L.ndc + L.nvdc)
    v = view(uw, vr)
    d = work.dwork
    transportcurrent!(d, t, v)
    @inbounds for k in eachindex(d)
        Fw[k] += d[k]
    end
    if residual
        transportresidual!(view(Fw, vr), t, v)
    else
        mul!(view(Fw, vr), t.Y, v)
    end
    # the blocks' own zero frequency rows, which replace `i = 0`, and the
    # currents they exchange across a component boundary
    br = work.blockrows
    if !isnothing(br)
        addblockdc!(Fw, br, uw, v)
        addblocktransport!(view(Fw, vr), br, uw)
    end
    # the rows spent on the directions nothing determines, written over what
    # was there: the equation they replaced was the redundant one
    pn = work.pinning
    if !isnothing(pn)
        idx = work.dclocal
        @inbounds for j in eachindex(pn.rows)
            Fw[idx[pn.rows[j]]] = uw[idx[pn.cols[j]]]
        end
    end
    return Fw
end

# =====================================================================
# Evaluating in canonical coordinates.
#
# The residual, the Jacobian vector product and the preconditioner are all
# written against the internal layout, and the canonical vector holds that
# layout in its first block. So each is handed a view of that block and
# nothing is copied; what is added is the direct current block, on the
# window.

"""
    canonicalresidual(fjreal!, work::CanonicalWork)

Wrap an internal coordinate residual and Jacobian closure so it takes and
returns canonical vectors.
"""
function canonicalresidual(fjreal!, work::CanonicalWork)
    L = work.layout
    return function (Fc, Jr, uc)
        fjreal!(isnothing(Fc) ? nothing : internalpart(Fc, L), Jr,
            internalpart(uc, L))
        isnothing(Fc) || addtransport!(Fc, work, uc)
        return nothing
    end
end

"""
    canonicaljvp(jvpreal!, work::CanonicalWork)

Wrap an internal coordinate Jacobian vector product so it takes and returns
canonical vectors.
"""
function canonicaljvp(jvpreal!, work::CanonicalWork)
    L = work.layout
    return function (Jvc, vc)
        jvpreal!(internalpart(Jvc, L), internalpart(vc, L))
        # the transport rows are affine, so their product drops the constant
        addtransport!(Jvc, work, vc; residual = false)
        return Jvc
    end
end

"""
    CanonicalPreconditioner(inner, work::CanonicalWork)

An internal coordinate preconditioner presented in canonical coordinates.

The mode blocks the inner preconditioner is built from include the zero
frequency mode, so it is applied where it was built and the canonical
vectors are transformed around it rather than the preconditioner being
rederived. That is exact, and it keeps the two paths taking the same
iterations.
"""
struct CanonicalPreconditioner{P,W<:CanonicalWork,F,D} <: AbstractPreconditioner
    inner::P
    work::W
    Yfact::F          # the direct current block, factorized; `nothing` without one
    dcindices::Vector{Int}
    dcwork::Vector{Float64}
    device::D         # the same solve, resident where the state is
end

function CanonicalPreconditioner(inner, work::CanonicalWork, F)
    idx = isnothing(F) ? Int[] : dcsubsystemindices(work)
    dev = (isnothing(F) || _onhost(work.xint) ||
           length(idx) > DCDEVICESOLVEMAX) ? nothing :
        DCFactorization(F, idx, KernelAbstractions.get_backend(work.xint))
    return CanonicalPreconditioner(inner, work, F, idx,
        zeros(Float64, length(idx)), dev)
end

function CanonicalPreconditioner(inner, work::CanonicalWork)
    t = work.transport
    isnothing(t) && return CanonicalPreconditioner(inner, work, nothing)
    # The direct current subsystem is small, dense and constant, with one
    # unknown per static flux component and one per block port current, and
    # its rows see no periodic state, so it is factorized once and solved
    # exactly. It has to be solved jointly: the transport rows carry the
    # block currents and the block rows carry the average voltages, and
    # solving only the transport half leaves the coupling to the Krylov
    # iteration, which stalls it.
    F = lu(dcsubsystem(work))
    return CanonicalPreconditioner(inner, work, F)
end

# =====================================================================
# The direct current solve, where the state is.
#
# The subsystem is constant, so its exact solve is a fixed linear map, but
# that map cannot be stored as a matrix: the subsystem is ill conditioned
# (a circuit with a ten million to one impedance ratio gives it a condition
# number of that order) and its unknowns are large in the solver's scaled
# units, so an explicit inverse leaves an absolute error the solve cannot
# drive to zero. It stays a linear solve, and what moves to the device is
# the factorization rather than the answer: the same factors, permutation
# and substitution order as the host path, in a kernel with one work item,
# so the two agree exactly and the window does not cross the bus per
# application.
#
# One work item because the substitutions are sequential, so the kernel's
# cost grows with the subsystem while the host copies it replaces cost a
# fixed latency. It wins for a handful of unknowns and loses badly beyond
# that. A handful is the common case, since there is one unknown per
# floating static flux component and one per scattering block port current;
# a larger subsystem keeps the host path, which is correct and merely
# copies.

# the largest subsystem the sequential device solve is launched for; the
# measured crossover against the host copies
const DCDEVICESOLVEMAX = 8

"""
    DCFactorization

The factorization of the direct current subsystem and the indices it acts
on, resident on a backend.

`factors` is the packed unit lower and upper triangle `lu` produced, `perm`
its row permutation, and `index` the canonical position of each subsystem
coordinate.
"""
struct DCFactorization{V,I}
    index::I
    perm::I
    factors::V
    work::V
    n::Int
end

function DCFactorization(F::LinearAlgebra.LU, index::Vector{Int}, backend)
    n = size(F.factors, 1)
    return DCFactorization(tobackend(backend, index),
        tobackend(backend, Vector{Int}(F.p)),
        tobackend(backend, Vector{Float64}(vec(F.factors))),
        tobackend(backend, zeros(Float64, n)), n)
end

# One work item: the substitutions are sequential and the system is a
# handful of unknowns, so what is being avoided is the bus and not the
# arithmetic.
@kernel function dcsolvekernel!(z, @Const(r), @Const(index),
        @Const(perm), @Const(factors), b, n)
    @index(Global)
    @inbounds begin
        for k in 1:n
            b[k] = r[index[perm[k]]]
        end
        for k in 1:n, i in (k+1):n            # L y = P b, unit lower
            b[i] -= factors[i + (k-1)*n]*b[k]
        end
        for k in n:-1:1                       # U x = y
            for i in (k+1):n
                b[k] -= factors[k + (i-1)*n]*b[i]
            end
            b[k] /= factors[k + (k-1)*n]
        end
        for k in 1:n
            z[index[k]] = b[k]
        end
    end
end

"""
    applydcsolve!(z, r, d::DCFactorization)

Solve the direct current subsystem for the coordinates of `r` it owns and
write the answer into the same coordinates of `z`, in place and on the
backend the factorization lives on.
"""
function applydcsolve!(z::AbstractVector, r::AbstractVector,
        d::DCFactorization)
    backend = KernelAbstractions.get_backend(z)
    kernel! = dcsolvekernel!(backend)
    kernel!(z, r, d.index, d.perm, d.factors, d.work, d.n; ndrange = 1)
    KernelAbstractions.synchronize(backend)
    return z
end

"""
    dcsubsystem(work::CanonicalWork)

The dense direct current block `[v; i]` of the canonical Jacobian: the
transport rows and the blocks' zero frequency rows, in the order
[`dcsubsystemindices`](@ref) gives.
"""
function dcsubsystem(work::CanonicalWork)
    t, br = work.transport, work.blockrows
    nc = nvoltages(t)
    idx = dcsubsystemindices(work)
    nb = length(idx) - nc
    A = zeros(Float64, nc + nb, nc + nb)
    A[1:nc, 1:nc] .= t.Y
    # the reference rows are written last and unconditionally: the transport
    # rows no longer carry one, so a circuit with no blocks still needs them.
    # A block current is named by its window position, as the block rows
    # name it.
    slots = dcsubsystemlocal(work)
    local_ = Dict(slots[nc+k] => nc + k for k in 1:nb)
    if !isnothing(br)
        # the block currents each component exchanges across its boundary
        for (c, ci, sgn) in br.transportterms
            A[c, local_[ci]] += sgn
        end
        # and the blocks' own rows: B0 (scale dv) - C0 i
        for (b, d) in enumerate(br.descriptors)
            ci = br.currentindex[b]
            sc, rc = br.signalcomponent[b], br.refcomponent[b]
            for p in eachindex(ci)
                row = local_[ci[p]]
                for q in eachindex(ci)
                    A[row, local_[ci[q]]] -= d.C0[p,q]
                    iszero(sc[q]) || (A[row, sc[q]] += d.B0[p,q]*br.scale)
                    iszero(rc[q]) || (A[row, rc[q]] -= d.B0[p,q]*br.scale)
                end
            end
        end
    end
    pn = work.pinning
    if !isnothing(pn)
        for j in eachindex(pn.rows)
            A[pn.rows[j], :] .= 0.0
            A[pn.rows[j], pn.cols[j]] = 1.0
        end
    end
    return A
end

"""
    dcsubsystemindices(work::CanonicalWork)

The canonical positions the direct current subsystem occupies: the explicit
voltage block, then every block port's zero frequency current.
"""
function dcsubsystemindices(work::CanonicalWork)
    L = work.layout
    return [windowindex(L, k) for k in dcsubsystemlocal(work)]
end

"""
    dcsubsystemlocal(work::CanonicalWork)

The same positions, local to the window.
"""
function dcsubsystemlocal(work::CanonicalWork)
    L, br = work.layout, work.blockrows
    idx = collect(L.ndc .+ (1:L.nvdc))
    isnothing(br) || for v in br.currentindex, i in v
        push!(idx, i)
    end
    return idx
end

function updatepreconditioner!(pc::CanonicalPreconditioner, u::AbstractVector)
    updatepreconditioner!(pc.inner, internalpart(u, pc.work.layout))
    return pc
end

# Block diagonal, and that is enough. The canonical Jacobian is
#
#     [ Jpp  Jpd ]    p: every flux entry, the zero frequency ones included
#     [  0   Jdd ]    d: the average voltages and the block port currents
#
# since the transport rows and the block relations see no periodic state,
# while the average voltages drive the resistor current `G0 P v` into the
# nodal zero frequency rows and a block current enters its two terminals'
# rows. Solving `d` exactly, the flux block with the inner preconditioner,
# and dropping `Jpd` leaves an error whose range is the columns of `Jpd`,
# one per direct current unknown, and a Krylov iteration removes an error
# of that rank in as many steps. The exact triangular form, `d` first and
# `Jpd d` taken off the nodal rows before the inner solve, was measured to
# save no iterations on junction chains with a resistive bias network or
# on a bridge between nearly open ports, and with an exact inner
# preconditioner both take one step per solve; it cost a pass over the
# nodal rows per application, so it is not here.
function applypreconditioner!(z::AbstractVector, pc::CanonicalPreconditioner,
        r::AbstractVector)
    L = pc.work.layout
    applypreconditioner!(internalpart(z, L), pc.inner, internalpart(r, L))
    return _solvedcblock!(z, pc, r)
end

# the direct current coordinates of `z`, solved exactly from `r`
function _solvedcblock!(z::AbstractVector, pc::CanonicalPreconditioner,
        r::AbstractVector)
    isnothing(pc.Yfact) && return z
    # where the state is not on the host, the same factors and the same
    # substitutions run there; the subsystem used to cross the bus three
    # times for this and that cost more than the residual it preconditions
    if !isnothing(pc.device)
        applydcsolve!(z, r, pc.device)
        return z
    end
    idx = pc.dcindices
    b = pc.dcwork
    @inbounds for k in eachindex(idx)
        b[k] = r[idx[k]]
    end
    ldiv!(pc.Yfact, b)
    # overwritten, not added: this solve is exact on these coordinates
    # and the inner preconditioner's guess at them is not
    @inbounds for k in eachindex(idx)
        z[idx[k]] = b[k]
    end
    return z
end

# The rest of the preconditioner interface is forwarded to the inner
# preconditioner rather than inherited: the defaults are inert (no
# escalation, no deflation), so a wrapper which did not forward them would
# silently turn escalation off, and on a problem which needs the full
# Jacobian that is the difference between converging and stalling.
escalatepreconditioner!(pc::CanonicalPreconditioner) =
    escalatepreconditioner!(pc.inner)
deflationsize(pc::CanonicalPreconditioner) = deflationsize(pc.inner)
deflationrebuilds(pc::CanonicalPreconditioner) = deflationrebuilds(pc.inner)
deflationform(pc::CanonicalPreconditioner) = deflationform(pc.inner)
pointmoved!(pc::CanonicalPreconditioner) = (pointmoved!(pc.inner); pc)
stalled!(pc::CanonicalPreconditioner) = (stalled!(pc.inner); pc)

candidatecount(pc::CanonicalPreconditioner) = candidatecount(pc.inner)
deflationproducts(pc::CanonicalPreconditioner) = deflationproducts(pc.inner)
harvest!(pc::CanonicalPreconditioner, ws::GMRESWorkspace, out::NamedTuple) =
    (harvest!(pc.inner, ws, out); pc)

# Exact when the inner is: the error the wrapper then leaves, the dropped
# `Jpd` block, has rank at most the number of direct current unknowns and
# costs that many Arnoldi steps per solve, which is fewer than the products
# and exact solves a deflation of it would cost at every Newton step. A
# recycling wrapper outside this one therefore clears its pair after an
# escalation rather than learning `Jpd`.
isexactpreconditioner(pc::CanonicalPreconditioner) =
    isexactpreconditioner(pc.inner)

"""
    dcsubsystemrhs(work::CanonicalWork)

The constant side of the direct current subsystem: the injected current on
the transport rows, zero on the blocks' own rows.
"""
function dcsubsystemrhs(work::CanonicalWork)
    t = work.transport
    n = length(dcsubsystemindices(work))
    b = zeros(Float64, n)
    b[1:nvoltages(t)] .= t.j
    return b
end

# The direct current subsystem mixes volts and amperes, and its rows are
# Kirchhoff sums in one place and constitutive relations in another, so its
# entries carry the circuit's impedance scale. A rank decision on the raw
# matrix would then depend on that scale: the same circuit written at a
# different impedance could be called singular or not. Scaling the rows
# and the columns to unit infinity norm first (two passes) removes the
# units from a question whose answer is a structural fact about the
# circuit.
function equilibrate(A::AbstractMatrix)
    dr = ones(Float64, size(A, 1))
    dc = ones(Float64, size(A, 2))
    B = copy(A)
    for _ in 1:2
        for i in axes(B, 1)
            m = maximum(abs, view(B, i, :); init = 0.0)
            iszero(m) && continue
            dr[i] /= m
            view(B, i, :) ./= m
        end
        for j in axes(B, 2)
            m = maximum(abs, view(B, :, j); init = 0.0)
            iszero(m) && continue
            dc[j] /= m
            view(B, :, j) ./= m
        end
    end
    return B, dr, dc
end

"""
    dccoupling(work::CanonicalWork)

`H`: the zero frequency nodal current each direct current unknown drives,
with the unknowns in the order [`dcsubsystemindices`](@ref) gives.

This is everything outside the direct current subsystem which can see those
unknowns. The average voltages reach the rest of the residual only as the
resistor current `G0 P v` added to the zero frequency nodal rows, and a
block port current only as the `+1` and `-1` it contributes at its two
terminals. Nothing else in the harmonic residual reads them, so a direction
with `H N = 0` is invisible to the whole problem and not only to the
subsystem, which is the condition a reference has to meet.

The transport rows are the component sums of these same nodal rows, which is
why `H` and the subsystem cannot disagree about a sign: `Y = P' (G0 P)` and
a block current enters its signal component's row with the sign it enters
its signal node's row.
"""
function dccoupling(work::CanonicalWork)
    t, br = work.transport, work.blockrows
    nc = nvoltages(t)
    n = size(t.coupling, 1)
    nb = length(dcsubsystemindices(work)) - nc
    I, J, V = Int[], Int[], Float64[]
    C = t.coupling
    for j in axes(C, 2), k in nzrange(C, j)
        push!(I, C.rowval[k]); push!(J, j); push!(V, C.nzval[k])
    end
    if !isnothing(br)
        c = nc
        for (b, d) in enumerate(br.descriptors)
            for p in eachindex(br.currentindex[b])
                c += 1
                # node 1 is ground and has no row
                d.signalnodes[p] > 1 &&
                    (push!(I, d.signalnodes[p]-1); push!(J, c); push!(V, 1.0))
                d.refnodes[p] > 1 &&
                    (push!(I, d.refnodes[p]-1); push!(J, c); push!(V, -1.0))
            end
        end
    end
    return sparse(I, J, V, n, nc + nb)
end

# the name of the block each block current coordinate belongs to, for
# messages; the voltage coordinates come first and have none
function dccoordinatenames(work::CanonicalWork)
    t, br = work.transport, work.blockrows
    names = fill("", nvoltages(t))
    isnothing(br) && return names
    for (b, d) in enumerate(br.descriptors)
        for _ in eachindex(br.currentindex[b])
            push!(names, d.name)
        end
    end
    return names
end

"""
    dcpinning(work::CanonicalWork)

Return the [`DCPinning`](@ref) a singular direct current subsystem needs,
`nothing` when it is nonsingular, or throw when it has no solution or an
undetermined direction the rest of the circuit can see.

The subsystem passed in must be the complete unreferenced descriptor: every
transport row and every block relation, with no reference chosen. Choosing
one earlier, from the resistors alone, can discard a row a block has made
necessary, and no later check can recover it.
"""
function dcpinning(work::CanonicalWork, nodenames = String[])
    isnothing(work.transport) && return nothing
    A = dcsubsystem(work)
    isempty(A) && return nothing
    B, dr, dc = equilibrate(A)
    F = svd(B)
    tol = maximum(F.S; init = 0.0) * maximum(size(B)) * eps()
    k = count(<=(tol), F.S)
    k == 0 && return nothing

    # the left null space: directions in which the equations say nothing, so
    # a constant side with a component along one of them cannot be met
    Y = F.U[:, end-k+1:end]
    b = dcsubsystemrhs(work) .* dr
    if norm(Y'b) > max(tol, eps()) * max(1.0, norm(b))
        throw(ArgumentError(lazy"No direct current solution exists: direct current is injected into a subnetwork which has no path carrying it away. The zero frequency mode is the average voltage, so a subnetwork whose average voltage is unconstrained cannot absorb a net current; give it a path to ground, or drive it differentially."))
    end

    # Which of the undetermined directions are gauges: a direction is one
    # only if the rest of the residual cannot see it, `H N = 0`. Anything
    # else is a physical quantity the circuit leaves undetermined, which is
    # refused rather than pinned to one of infinitely many answers.
    Nhat = F.V[:, end-k+1:end]
    H = dccoupling(work) * Diagonal(dc)
    G = Matrix(H * Nhat)
    hs = maximum(abs, H; init = 0.0)
    gtol = max(hs, 1.0) * maximum(size(G)) * sqrt(eps())
    if any(>(gtol), svdvals(G))
        names = dccoordinatenames(work)
        w = vec(maximum(abs, Nhat; dims = 2))
        seen = String[]
        for c in eachindex(names)
            isempty(names[c]) && continue
            w[c] <= sqrt(eps()) && continue
            names[c] in seen || push!(seen, names[c])
        end
        # the joined list is built outside the message: a comma inside a
        # `lazy` interpolation ends the interpolated expression
        involved = join(seen, ", ")
        throw(ArgumentError(lazy"The direct current network leaves a branch current undetermined, and that current is visible to the rest of the circuit: changing it moves the current at a node, and with it the static flux of any inductor or junction in parallel. The blocks whose zero frequency currents are involved are $(involved). An ideal short or through in parallel with an inductive branch has no unique direct current solution; give the block a finite series impedance at zero frequency, or give the parallel branch one, so that the division is determined."))
    end

    # Which equations to give up, and which coordinate each one fixes.
    # Pivoting on the left null space picks rows which carry the redundancy,
    # so what is left still spans the row space; pivoting on the null space
    # picks coordinates the directions move, so the references are
    # independent. Both are done in the equilibrated coordinates.
    rows = sort!(qr(Y', ColumnNorm()).p[1:k])
    cols = sort!(qr(Nhat', ColumnNorm()).p[1:k])
    pn = DCPinning(rows, cols)

    # the references have to leave a nonsingular system, which the two
    # pivoted choices give but do not guarantee jointly
    Ap = copy(A)
    for j in 1:k
        Ap[rows[j], :] .= 0.0
        Ap[rows[j], cols[j]] = 1.0
    end
    Bp, _, _ = equilibrate(Ap)
    if minimum(svdvals(Bp)) <= maximum(size(Bp)) * eps()
        error("the direct current references left a singular subsystem, which they must not: this is a bug in `dcpinning`, not a property of the circuit.")
    end
    return pn
end


# The canonical Jacobian, as a plan.
#
# Every term the direct current block contributes is constant: the transport
# rows, the coupling `G0 P`, the block pencils `B0` and `C0`, the boundary
# currents and the reference rows do not depend on the periodic state. What
# moves between Newton iterations is the internal Jacobian, and only its
# values, since its pattern and the permutation which reorders it are
# fixed. So the whole assembly is a fixed pattern, a fixed scatter of the
# internal values into it, and a fixed list of constant additions, found
# once.
#
# Without the plan the rebuild costs several times the evaluation of the
# internal Jacobian and the factorization of the result together, most of
# the linear algebra of a Newton step spent rediscovering a pattern which
# has not moved.

"""
    CanonicalJacobianPlan

The canonical Jacobian's pattern and the fixed arithmetic which fills it.

# Fields
- `J`: the pattern, and the buffer the values are written into.
- `source`: for each stored entry of the internal Jacobian, where it lands,
  or zero when it lands in a row a reference replaces.
- `fixedindex`, `fixedvalue`: the direct current block's constant entries,
  including the single one of each reference row.

See [`canonicaljacobianplan`](@ref) and [`canonicaljacobian!`](@ref).
"""
struct CanonicalJacobianPlan
    J::SparseMatrixCSC{Float64,Int}
    source::Vector{Int}
    fixedindex::Vector{Int}
    fixedvalue::Vector{Float64}
end

# the position of (i, j) in the value array of `S`, whose row indices are
# sorted within a column
function nzposition(S::SparseMatrixCSC, i::Integer, j::Integer)
    r = nzrange(S, j)
    k = searchsortedfirst(view(S.rowval, r), i)
    (k <= length(r) && S.rowval[r[k]] == i) ||
        error("the canonical Jacobian's pattern is missing an entry it was built to hold, which is a bug in `canonicaljacobianplan`.")
    return r[k]
end

# The direct current block's constant entries, in the canonical numbering of
# the whole matrix. This is the same arithmetic the residual does, read as a
# matrix; it is written once here rather than once per iteration.
function dcjacobianentries(work::CanonicalWork)
    L = work.layout
    n = L.rdim
    dcpos = L.dcpos
    I, J, V = Int[], Int[], Float64[]
    t = work.transport
    isnothing(t) && return I, J, V

    # the resistor current the voltages drive into the zero frequency nodal
    # rows, and the transport rows themselves
    C = t.coupling
    for j in axes(C, 2), k in nzrange(C, j)
        push!(I, dcpos[C.rowval[k]]); push!(J, n + j); push!(V, C.nzval[k])
    end
    for j in axes(t.Y, 2), i in axes(t.Y, 1)
        iszero(t.Y[i,j]) && continue
        push!(I, n + i); push!(J, n + j); push!(V, t.Y[i,j])
    end

    br = work.blockrows
    if !isnothing(br)
        # the block currents each component exchanges across its boundary
        for (c, ci, sgn) in br.transportterms
            push!(I, n + c); push!(J, dcpos[ci]); push!(V, float(sgn))
        end
        # and each block's own row, replacing the `-i` already in the stamp
        for (b, d) in enumerate(br.descriptors)
            ci = br.currentindex[b]
            sc, rc = br.signalcomponent[b], br.refcomponent[b]
            for p in eachindex(ci)
                rp = dcpos[ci[p]]
                push!(I, rp); push!(J, rp); push!(V, 1.0)
                for q in eachindex(ci)
                    push!(I, rp); push!(J, dcpos[ci[q]]); push!(V, -d.C0[p,q])
                    w = d.B0[p,q]*br.scale
                    iszero(sc[q]) ||
                        (push!(I, rp); push!(J, n + sc[q]); push!(V, w))
                    iszero(rc[q]) ||
                        (push!(I, rp); push!(J, n + rc[q]); push!(V, -w))
                end
            end
        end
    end
    return I, J, V
end

"""
    canonicaljacobianplan(Jint::SparseMatrixCSC, work::CanonicalWork)

Build the [`CanonicalJacobianPlan`](@ref) for an internal Jacobian pattern.

Only the pattern of `Jint` is read, not its values, so the plan is valid for
every point the solve visits. It stops being valid if that pattern moves,
which it must not.
"""
function canonicaljacobianplan(Jint::SparseMatrixCSC, work::CanonicalWork)
    L = work.layout
    n = L.rdim
    N = n + L.nvdc

    # where each stored entry of the internal Jacobian lands: where it is,
    # since the internal state is the first block of the canonical one
    di = Vector{Int}(undef, nnz(Jint))
    dj = Vector{Int}(undef, nnz(Jint))
    for col in axes(Jint, 2), k in nzrange(Jint, col)
        di[k] = Jint.rowval[k]
        dj[k] = col
    end

    fi, fj, fv = dcjacobianentries(work)

    # a reference row is a replacement, so nothing else in it survives
    pn = work.pinning
    dead = Set{Int}()
    refi, refj = Int[], Int[]
    if !isnothing(pn)
        idx = work.dcindex
        for j in eachindex(pn.rows)
            push!(dead, idx[pn.rows[j]])
            push!(refi, idx[pn.rows[j]])
            push!(refj, idx[pn.cols[j]])
        end
    end

    live(i) = !(i in dead)
    keepd = [k for k in eachindex(di) if live(di[k])]
    keepf = [k for k in eachindex(fi) if live(fi[k])]
    I = vcat(di[keepd], fi[keepf], refi)
    Jc = vcat(dj[keepd], fj[keepf], refj)
    J = sparse(I, Jc, ones(Float64, length(I)), N, N)

    source = zeros(Int, length(di))
    for k in keepd
        source[k] = nzposition(J, di[k], dj[k])
    end
    fixedindex = Int[nzposition(J, fi[k], fj[k]) for k in keepf]
    fixedvalue = Float64[fv[k] for k in keepf]
    for k in eachindex(refi)
        push!(fixedindex, nzposition(J, refi[k], refj[k]))
        push!(fixedvalue, 1.0)
    end
    return CanonicalJacobianPlan(J, source, fixedindex, fixedvalue)
end

"""
    canonicaljacobian!(plan::CanonicalJacobianPlan, Jint::SparseMatrixCSC)

Fill the plan's matrix from an internal Jacobian and return it.

The flux block is `Jint` as it is, since the internal state is the first
block of the canonical one. What is added is the explicit direct
current block and its couplings, in the same places the residual adds them:
the resistor current the average voltages drive into the zero frequency
nodal rows, the transport rows and the block currents they carry across a
component boundary, and each block's own zero frequency row in place of the
`i = 0` the stamp wrote. A reference row is a replacement rather than an
addition, which is why the internal entries landing in it are dropped when
the plan is built rather than zeroed here.

The matrix-free product is the reference this is checked against: for every
unit vector the two must agree exactly, which is a sharper test than any
finite difference of the residual would be.
"""
function canonicaljacobian!(plan::CanonicalJacobianPlan,
        Jint::SparseMatrixCSC)
    length(plan.source) == nnz(Jint) || throw(DimensionMismatch(
        lazy"the internal Jacobian has $(nnz(Jint)) stored entries but the plan was built for $(length(plan.source)); its pattern must not move between iterations."))
    nz = plan.J.nzval
    fill!(nz, 0.0)
    src = plan.source
    v = Jint.nzval
    @inbounds for k in eachindex(src)
        d = src[k]
        iszero(d) || (nz[d] += v[k])
    end
    @inbounds for k in eachindex(plan.fixedindex)
        nz[plan.fixedindex[k]] += plan.fixedvalue[k]
    end
    return plan.J
end

"""
    canonicaljacobian(Jint::SparseMatrixCSC, work::CanonicalWork)

The Jacobian in canonical coordinates, assembled from the internal one.

This builds a plan and applies it, which is what a caller wanting one matrix
at one point should do. A solve builds the plan once instead; see
[`canonicaljacobianplan`](@ref).
"""
canonicaljacobian(Jint::SparseMatrixCSC, work::CanonicalWork) =
    canonicaljacobian!(canonicaljacobianplan(Jint, work), Jint)

"""
    canonicalfj(fjreal!, work::CanonicalWork, Jint, plan)

Wrap an internal coordinate residual and Jacobian closure for a direct solve
method: the residual in canonical coordinates, and the canonical Jacobian
filled through `plan`.

The matrix is filled rather than replaced, because the solver factorizes the
one it was handed. Its pattern is fixed -- the internal pattern under a
permutation plus the direct current block, none of which moves between
iterations -- and the plan is built from it once.
"""
function canonicalfj(fjreal!, work::CanonicalWork, Jint,
        plan::CanonicalJacobianPlan)
    L = work.layout
    return function (Fc, Jout, uc)
        fjreal!(isnothing(Fc) ? nothing : internalpart(Fc, L),
            isnothing(Jout) ? nothing : Jint, internalpart(uc, L))
        isnothing(Fc) || addtransport!(Fc, work, uc)
        if !isnothing(Jout)
            canonicaljacobian!(plan, Jint)
            Jout === plan.J || copyto!(Jout.nzval, plan.J.nzval)
        end
        return nothing
    end
end

# =====================================================================
# The direct current block as a matrix.
#
# Every term of `addtransportwindow!` is linear in the window, and some rows
# are added to while others are overwritten, so the whole update is
#
#     Fw <- keep .* Fw + M uw + c
#
# with `keep` zero on the overwritten rows. In that form it is three array
# operations rather than a walk over scattered indices, so it can run where
# the state lives instead of being copied to the host and back.
#
# `M` and `c` are read off the scalar implementation by probing it one
# basis vector at a time rather than assembled a second time by hand, which
# is O(window) work once and makes the two forms agree by construction.

"""
    DCUpdate

The direct current block's contribution to the residual as `keep .* Fw +
M*uw + c`, in whatever array type the state uses.

`cresidual` carries the injected current; the Jacobian vector product uses
the same `M` and `keep` with no constant, the two differing only by that.
"""
struct DCUpdate{V,I}
    keep::V
    rowptr::I          # the matrix by rows, so one work item owns one row
    colval::I
    nzval::V
    cresidual::V
end

# One work item per row of the window: its own entry, kept or not, plus its
# row of the matrix, plus the constant when this is a residual.
@kernel function dcupdatekernel!(Fw, @Const(keep), @Const(rowptr),
        @Const(colval), @Const(nzval), @Const(uw), @Const(c), alpha)
    i = @index(Global)
    @inbounds begin
        acc = keep[i]*Fw[i] + alpha*c[i]
        for k in rowptr[i]:(rowptr[i+1] - 1)
            acc += nzval[k]*uw[colval[k]]
        end
        Fw[i] = acc
    end
end

"""
    dcupdate(work::CanonicalWork)

Build the [`DCUpdate`](@ref) by probing `addtransportwindow!`, so the matrix
form and the scalar form agree by construction.
"""
function dcupdate(work::CanonicalWork)
    L = work.layout
    nw = L.ndc + L.nvdc
    nw == 0 && return nothing
    Fw = zeros(Float64, nw); uw = zeros(Float64, nw)

    # `keep` is one where the row is added to and zero where it is written
    fill!(Fw, 1.0); fill!(uw, 0.0)
    addtransportwindow!(Fw, uw, work; residual = false)
    keep = copy(Fw)

    # the constant, from a zero point
    fill!(Fw, 0.0)
    addtransportwindow!(Fw, uw, work; residual = true)
    c = copy(Fw)

    # and the linear part, column by column, against the product form which
    # carries no constant
    I, J, V = Int[], Int[], Float64[]
    for k in 1:nw
        fill!(Fw, 0.0); fill!(uw, 0.0); uw[k] = 1.0
        addtransportwindow!(Fw, uw, work; residual = false)
        for i in 1:nw
            iszero(Fw[i]) && continue
            push!(I, i); push!(J, k); push!(V, Fw[i])
        end
    end
    # by rows: the transpose of a column major sparse matrix is the row
    # major form of the original, which is what the kernel walks
    Mt = sparse(J, I, V, nw, nw)
    return DCUpdate(keep, Mt.colptr, Mt.rowval, Mt.nzval, c)
end

"""
    applydcupdate!(Fw, uw, up::DCUpdate; residual = true)

Apply the direct current block in its matrix form.
"""
function applydcupdate!(Fw::AbstractVector, uw::AbstractVector,
        up::DCUpdate; residual::Bool = true)
    backend = KernelAbstractions.get_backend(Fw)
    kernel! = dcupdatekernel!(backend)
    kernel!(Fw, up.keep, up.rowptr, up.colval, up.nzval, uw, up.cresidual,
        residual ? one(eltype(Fw)) : zero(eltype(Fw));
        ndrange = length(Fw))
    KernelAbstractions.synchronize(backend)
    return Fw
end
