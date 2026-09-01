# The canonical state layout, and the transform between it and the layout
# the solver evaluates in.
#
# The internal real state is node major and mode minor: every node
# contributes one entry per self conjugate mode and two per conjugate pair,
# in mode order, so the zero frequency entry of a node sits between that
# node's alternating current entries. That is the layout `ModeLayout`
# describes and every kernel in the solve is written against.
#
# The canonical state groups by role instead:
#
#     [ phiac | phidc | vdc | idc ]
#
# `phiac` is every nonzero frequency flux entry, `phidc` the zero frequency
# flux of each node, and `vdc` and `idc` the explicit direct current node
# voltages and branch currents that stage 5 adds. Grouping by role is what
# makes the direct current block a contiguous subsystem with its own rows,
# which is the whole point: the transport rows can then be assembled,
# factorized and solved apart from the alternating current problem.
#
# This file is the transitional half of that change. `vdc` and `idc` are
# empty here, so the canonical state is a permutation of the internal one
# and nothing about the physics moves. What it buys is the layout itself,
# and the measurement of what bracketing every residual and every Jacobian
# vector product with a scatter and a gather costs, before anything harder
# is built on top of it.

"""
    CompositeLayout

The canonical state layout, `[phiac | phidc | vdc | idc]`, together with the
map into the internal layout the solver evaluates in.

# Fields
- `nac`, `ndc`, `nvdc`, `nidc`: the length of each block. `nvdc` and `nidc`
  are zero until explicit direct current coordinates exist.
- `rdim`: the length of the internal real state.
- `perm`: for each canonical index, the internal index holding it. Length
  `nac + ndc`; the `vdc` and `idc` blocks have no internal counterpart.

While `nvdc` and `nidc` are zero, `perm` is a permutation of `1:rdim`, so
[`scattercanonical!`](@ref) and [`gathercanonical!`](@ref) are exact
inverses and the canonical Jacobian is `P J Pᵀ` with `P` orthogonal. A
Krylov method is invariant under an orthogonal change of basis, so the two
paths take the same iterations on the same problem, which is what the tests
assert.

See [`compositelayout`](@ref).
"""
struct CompositeLayout{Ti<:Integer,V<:AbstractVector{Ti}}
    nac::Int
    ndc::Int
    nvdc::Int
    nidc::Int
    rdim::Int
    perm::V
end

"""
    canonicaldim(L::CompositeLayout)

The length of the canonical state.
"""
canonicaldim(L::CompositeLayout) = L.nac + L.ndc + L.nvdc + L.nidc

"""
    ispermutation(L::CompositeLayout)

Whether the canonical state is a permutation of the internal one, which it
is exactly while there are no explicit direct current coordinates.
"""
ispermutation(L::CompositeLayout) =
    iszero(L.nvdc) && iszero(L.nidc) && length(L.perm) == L.rdim

"""
    compositelayout(ml::ModeLayout, isdc::AbstractVector{Bool};
        nvdc = 0, nidc = 0)

Build the canonical layout for a state whose internal layout is `ml`, where
`isdc[t]` marks mode `t` as the zero frequency one.

A zero frequency mode is self conjugate, so it occupies one internal entry
rather than two; `isdc` must therefore imply `ml.isreal`. The converse does
not hold, since a Nyquist mode is self conjugate without being direct
current, which is why the split is made on `isdc` and not on `ml.isreal`.
"""
function compositelayout(ml::ModeLayout, isdc::AbstractVector{Bool};
        nvdc::Int = 0, nidc::Int = 0)
    length(isdc) == ml.nmodes || throw(DimensionMismatch(
        lazy"`isdc` has length $(length(isdc)) but the layout has $(ml.nmodes) modes."))
    for t in eachindex(isdc)
        if isdc[t] && !ml.isreal[t]
            throw(ArgumentError(lazy"mode $(t) is marked zero frequency but is not self conjugate; a zero frequency component of a real signal has no imaginary part."))
        end
    end
    (nvdc >= 0 && nidc >= 0) || throw(ArgumentError("block lengths must be nonnegative."))

    nmodes = ml.nmodes
    # count first so both blocks can be filled in one pass
    nper = ml.dim ÷ nmodes            # entries per mode across the state
    ndcper = count(isdc)
    ndc = nper * ndcper
    nac = ml.rdim - ndc

    perm = Vector{Int}(undef, nac + ndc)
    a, d = 1, nac + 1                 # cursors into the two blocks
    p, t = 1, 1                       # cursor into the internal real state
    @inbounds for _ in 1:ml.dim
        w = ml.isreal[t] ? 1 : 2
        if isdc[t]
            for k in 0:w-1
                perm[d] = p + k
                d += 1
            end
        else
            for k in 0:w-1
                perm[a] = p + k
                a += 1
            end
        end
        p += w
        t = _next(t, nmodes)
    end
    a == nac + 1 || error("internal error: filled $(a-1) of $(nac) alternating current entries.")
    d == nac + ndc + 1 || error("internal error: filled $(d-nac-1) of $(ndc) direct current entries.")

    return CompositeLayout(nac, ndc, nvdc, nidc, ml.rdim, perm)
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

Write the canonical form of the internal real state `rint` into `u`.

The `vdc` and `idc` blocks have no internal counterpart and are left
untouched, so a caller which keeps explicit direct current coordinates in
`u` does not lose them here.
"""
function gathercanonical!(u::AbstractVector, rint::AbstractVector,
        L::CompositeLayout)
    length(u) == canonicaldim(L) || throw(DimensionMismatch(
        lazy"the canonical state has length $(length(u)) but the layout needs $(canonicaldim(L))."))
    length(rint) == L.rdim || throw(DimensionMismatch(
        lazy"the internal state has length $(length(rint)) but the layout needs $(L.rdim)."))
    _gatherperm!(u, rint, L.perm)
    return u
end

# On the host a plain loop, on a device the kernel. A kernel launch and the
# synchronize after it cost more than the copy itself at these sizes: going
# through the device path on the CPU measured 26 us against 3.6 us for the
# loop on a 12240 entry state, which is the difference between this being a
# fifth of a Jacobian vector product and a twentieth of one.
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

"""
    scattercanonical!(rint, u, L::CompositeLayout)

Write the internal real form of the canonical state `u` into `rint`.

Every internal entry the layout covers is written, so `rint` needs no
initialization. While the layout is a permutation that is all of them; once
explicit direct current coordinates exist it will not be, and the entries
they replace are zeroed here rather than left stale.

The zero frequency flux entries are written from the `phidc` block, never
zeroed. They carry the static phase across a junction, and dropping them
would silently change `sin(phi)`.
"""
function scattercanonical!(rint::AbstractVector, u::AbstractVector,
        L::CompositeLayout)
    length(u) == canonicaldim(L) || throw(DimensionMismatch(
        lazy"the canonical state has length $(length(u)) but the layout needs $(canonicaldim(L))."))
    length(rint) == L.rdim || throw(DimensionMismatch(
        lazy"the internal state has length $(length(rint)) but the layout needs $(L.rdim)."))
    ispermutation(L) || fill!(rint, zero(eltype(rint)))
    _scatterperm!(rint, u, L.perm)
    return rint
end

"""
    tobackend(backend, L::CompositeLayout)

Move the index the scatter and the gather read to the given backend.
"""
function tobackend(backend::Backend, L::CompositeLayout)
    return CompositeLayout(L.nac, L.ndc, L.nvdc, L.nidc, L.rdim,
        tobackend(backend, L.perm))
end

"""
    DCPinning

The directions of the direct current subsystem which nothing determines,
and the rows to spend on fixing them.

A free port current can be a genuine gauge: two ideal shorts between one
pair of nodes divide the current between them in a way no equation fixes,
while every node voltage is determined. Refusing that would be refusing a
circuit which has an answer. Pinning it chooses one of the equivalent
current divisions, exactly as a floating island's reference voltage is
already chosen, and nothing observable moves.

What cannot be pinned is an inconsistent block -- direct current injected
into a subnetwork with nowhere to go -- which is singular *and* has no
solution. The two are told apart by whether the constant side lies in the
range of the matrix, not by inspecting the circuit.
"""
struct DCPinning
    rows::Vector{Int}      # positions in the subsystem to replace
    N::Matrix{Float64}     # the null directions, one per column
end

# The workspaces a canonical evaluation needs: one internal state to scatter
# into, one internal residual (or product) to gather out of, and the
# transport rows when the direct current block is explicit.
#
# `nnodaldc` splits the `phidc` block. Its first entries are the zero
# frequency flux of each node, which the transport coupling drives; any
# after them belong to auxiliary branch currents, which it does not.
struct CanonicalWork{T,V<:AbstractVector{T},L<:CompositeLayout,TR,BR}
    layout::L
    xint::V
    Fint::V
    transport::TR
    blockrows::BR
    dwork::Vector{Float64}
    nnodaldc::Int
    pinning::Union{Nothing,DCPinning}
    dcindex::Vector{Int}
    dclocal::Vector{Int}            # `dcindex`, local to the window
    blocklocal::Dict{Int,Int}       # canonical -> window, for the blocks
    Fwindow::Vector{Float64}        # device path only: the window on the host
    uwindow::Vector{Float64}
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
    nw = L.ndc + L.nvdc
    w = CanonicalWork(L, similar(proto, L.rdim), similar(proto, L.rdim),
        transport, blockrows, zeros(Float64, nd), nnodaldc, nothing, Int[],
        Int[], Dict{Int,Int}(), zeros(Float64, nw), zeros(Float64, nw))
    isnothing(transport) && return w
    # the pinning is read off the subsystem this work describes, so it is
    # found once here and then carried
    idx = dcsubsystemindices(w)
    shift = L.nac
    blocklocal = Dict{Int,Int}()
    if !isnothing(blockrows)
        for v in blockrows.currentindex, i in v
            blocklocal[i] = i - shift
        end
    end
    return CanonicalWork(L, w.xint, w.Fint, transport, blockrows, w.dwork,
        nnodaldc, dcpinning(w), idx, idx .- shift, blocklocal,
        w.Fwindow, w.uwindow)
end

# the canonical index range of the explicit voltage block
voltagerange(L::CompositeLayout) =
    (L.nac + L.ndc + 1):(L.nac + L.ndc + L.nvdc)

# The explicit direct current block's contribution: the resistor current the
# average voltages drive into the zero frequency nodal rows, and the
# transport rows themselves.
#
# The sign is the elimination's, read backwards. The residual is
# `A x - bnm`, and eliminating replaces `bnm` by `bnm - G0 P v`, so carrying
# `v` instead means adding `G0 P v` here and leaving the applied source
# alone. Doing both would count the resistor current twice.
# Everything the direct current block reads and writes lies in one
# contiguous window of the canonical vector: the zero frequency block and
# the explicit voltages, which are adjacent. So the arithmetic is done on
# that window with local indices, which is what lets it run on a device
# without scalar indexing -- the window is copied to the host, worked on,
# and copied back, and it is a few hundred numbers where the state is tens
# of thousands. On the host it is a view and there is no copy at all.
dcwindow(L::CompositeLayout) = (L.nac + 1):(L.nac + L.ndc + L.nvdc)

function addtransport!(Fc::AbstractVector, work::CanonicalWork,
        u::AbstractVector; residual::Bool = true)
    isnothing(work.transport) && return Fc
    win = dcwindow(work.layout)
    if _onhost(Fc)
        addtransportwindow!(view(Fc, win), view(u, win), work; residual)
        return Fc
    end
    Fw, uw = work.Fwindow, work.uwindow
    copyto!(Fw, view(Fc, win))
    copyto!(uw, view(u, win))
    addtransportwindow!(Fw, uw, work; residual)
    copyto!(view(Fc, win), Fw)
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
        mul!(view(Fw, vr), t.Ytr, v)
    end
    # the blocks' own zero frequency rows, which replace `i = 0`, and the
    # currents they exchange across a component boundary
    br = work.blockrows
    if !isnothing(br)
        addblockdc!(Fw, br, uw, v, work.blocklocal)
        addblocktransport!(view(Fw, vr), br, uw, t.pinned, work.blocklocal)
    end
    # the rows spent on the directions nothing determines, written over what
    # was there: the equation they replaced was the redundant one
    pn = work.pinning
    if !isnothing(pn)
        idx = work.dclocal
        @inbounds for j in axes(pn.N, 2)
            acc = zero(eltype(Fw))
            for k in eachindex(idx)
                acc += pn.N[k, j] * uw[idx[k]]
            end
            Fw[idx[pn.rows[j]]] = acc
        end
    end
    return Fw
end

# =====================================================================
# Evaluating in canonical coordinates.
#
# The residual, the Jacobian vector product and the preconditioner are all
# written against the internal layout. Rather than rewrite them, the
# canonical path brackets each with a scatter in and a gather out, which is
# what makes the two paths comparable: while the layout is a permutation
# the canonical operator is `P J Pᵀ` with `P` orthogonal, so a Krylov
# method takes the same iterations and reaches the same point, and any
# difference between the paths is a bug rather than a change of problem.
#
# The cost of that bracketing is the thing worth knowing before the rest of
# stage 5 is built on it, which is why it is measured rather than assumed.

"""
    canonicalresidual(fjreal!, work::CanonicalWork)

Wrap an internal coordinate residual and Jacobian closure so it takes and
returns canonical vectors.

`fjreal!` writes the canonical real representation of the evaluation point
back into its state argument, so the wrapper gathers that back out too;
otherwise the solver's point would drift from the system's.
"""
function canonicalresidual(fjreal!, work::CanonicalWork)
    L, xint, Fint = work.layout, work.xint, work.Fint
    return function (Fc, Jr, uc)
        scattercanonical!(xint, uc, L)
        fjreal!(isnothing(Fc) ? nothing : Fint, Jr, xint)
        if !isnothing(Fc)
            gathercanonical!(Fc, Fint, L)
            addtransport!(Fc, work, uc)
        end
        # writes the flux blocks only, so the explicit voltages the solver
        # is carrying in `uc` survive
        gathercanonical!(uc, xint, L)
        return nothing
    end
end

"""
    canonicaljvp(jvpreal!, work::CanonicalWork)

Wrap an internal coordinate Jacobian vector product so it takes and returns
canonical vectors.
"""
function canonicaljvp(jvpreal!, work::CanonicalWork)
    L, vint, Jvint = work.layout, work.xint, work.Fint
    return function (Jvc, vc)
        scattercanonical!(vint, vc, L)
        jvpreal!(Jvint, vint)
        gathercanonical!(Jvc, Jvint, L)
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
struct CanonicalPreconditioner{P,W<:CanonicalWork,F} <: AbstractPreconditioner
    inner::P
    work::W
    Yfact::F          # the direct current block, factorized; `nothing` without one
    dcindices::Vector{Int}
    dcwork::Vector{Float64}
end

function CanonicalPreconditioner(inner, work::CanonicalWork, F)
    idx = isnothing(F) ? Int[] : dcsubsystemindices(work)
    return CanonicalPreconditioner(inner, work, F, idx,
        zeros(Float64, length(idx)))
end

function CanonicalPreconditioner(inner, work::CanonicalWork)
    t = work.transport
    isnothing(t) && return CanonicalPreconditioner(inner, work, nothing)
    # The direct current subsystem is small, dense and constant: one unknown
    # per static flux component and one per block port current, and its rows
    # see no periodic state. Factorizing it once and solving it exactly
    # beats any approximation, and it has to be solved jointly -- the
    # transport rows carry the block currents and the block rows carry the
    # average voltages, so solving only the transport half leaves the
    # coupling to the Krylov iteration, which is what stalls it.
    F = lu(dcsubsystem(work))
    return CanonicalPreconditioner(inner, work, F)
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
    A[1:nc, 1:nc] .= t.Ytr
    isnothing(br) && return A
    local_ = Dict(idx[nc+k] => nc + k for k in 1:nb)
    # the block currents each component exchanges across its boundary
    for (c, ci, sgn) in br.transportterms
        (c in t.pinned) && continue
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
    pn = work.pinning
    if !isnothing(pn)
        for j in axes(pn.N, 2)
            A[pn.rows[j], :] .= view(pn.N, :, j)
        end
    end
    return A
end

"""
    dcsubsystemindices(work::CanonicalWork)

The canonical indices the direct current subsystem occupies: the explicit
voltage block, then every block port's zero frequency current.
"""
function dcsubsystemindices(work::CanonicalWork)
    L, br = work.layout, work.blockrows
    idx = collect(voltagerange(L))
    isnothing(br) || for v in br.currentindex, i in v
        push!(idx, i)
    end
    return idx
end

function updatepreconditioner!(pc::CanonicalPreconditioner, u::AbstractVector)
    w = pc.work
    scattercanonical!(w.xint, u, w.layout)
    updatepreconditioner!(pc.inner, w.xint)
    return pc
end

function applypreconditioner!(z::AbstractVector, pc::CanonicalPreconditioner,
        r::AbstractVector)
    w = pc.work
    L = w.layout
    scattercanonical!(w.xint, r, L)
    applypreconditioner!(w.Fint, pc.inner, w.xint)
    gathercanonical!(z, w.Fint, L)
    # Block diagonal, not block triangular: the coupling `Jpv` is dropped
    # here. The exact triangular solve is 4D's, and needs the flux block
    # solved first; as a preconditioner this is already the right shape,
    # since the block it adds is solved exactly.
    if !isnothing(pc.Yfact)
        # through the window, for the same reason `addtransport!` is: these
        # indices are scattered, and a device vector cannot be indexed one
        # element at a time
        win = dcwindow(L)
        rw, zw = w.uwindow, w.Fwindow
        copyto!(rw, view(r, win))
        copyto!(zw, view(z, win))
        loc = w.dclocal
        b = pc.dcwork
        @inbounds for k in eachindex(loc)
            b[k] = rw[loc[k]]
        end
        ldiv!(pc.Yfact, b)
        # overwritten, not added: this solve is exact on these coordinates
        # and the inner preconditioner's guess at them is not
        @inbounds for k in eachindex(loc)
            zw[loc[k]] = b[k]
        end
        copyto!(view(z, win), zw)
    end
    return z
end

# The rest of the preconditioner interface is delegated rather than
# inherited. Its defaults are inert -- no escalation, no deflation, an
# unknown level -- so a wrapper which does not forward them silently turns
# escalation off, and on a case hard enough to need the full Jacobian that
# is the difference between converging and stalling with an O(1e-2)
# residual. It does not look like a layout bug when it happens.
escalatepreconditioner!(pc::CanonicalPreconditioner) =
    escalatepreconditioner!(pc.inner)
deflationsize(pc::CanonicalPreconditioner) = deflationsize(pc.inner)
deflationrebuilds(pc::CanonicalPreconditioner) = deflationrebuilds(pc.inner)
escalationrefusals(pc::CanonicalPreconditioner) = escalationrefusals(pc.inner)
preconditionerlevel(pc::CanonicalPreconditioner) =
    preconditionerlevel(pc.inner)

"""
    dcsubsystemrhs(work::CanonicalWork)

The constant side of the direct current subsystem: the injected current on
the transport rows, zero on the blocks' own rows.
"""
function dcsubsystemrhs(work::CanonicalWork)
    t = work.transport
    n = length(dcsubsystemindices(work))
    b = zeros(Float64, n)
    b[1:nvoltages(t)] .= t.jtr
    return b
end

"""
    dcpinning(work::CanonicalWork)

Return the [`DCPinning`](@ref) a singular direct current subsystem needs,
`nothing` when it is nonsingular, or throw when it has no solution at all.
"""
function dcpinning(work::CanonicalWork, nodenames = String[])
    isnothing(work.transport) && return nothing
    A = dcsubsystem(work)
    isempty(A) && return nothing
    F = svd(A)
    tol = maximum(F.S; init = 0.0) * maximum(size(A)) * eps()
    k = count(<=(tol), F.S)
    k == 0 && return nothing

    # the left null space: directions in which the equations say nothing, so
    # a constant side with a component along one of them cannot be met
    Y = F.U[:, end-k+1:end]
    b = dcsubsystemrhs(work)
    if norm(Y'b) > max(tol, eps()) * max(1.0, norm(b))
        throw(ArgumentError(lazy"No direct current solution exists: direct current is injected into a subnetwork which has no path carrying it away. The zero frequency mode is the average voltage, so a subnetwork whose average voltage is unconstrained cannot absorb a net current; give it a path to ground, or drive it differentially."))
    end

    # spend the redundant rows on the free directions. The pivoted
    # factorization picks rows which actually carry the redundancy, so what
    # is left of the original system still spans its row space.
    rows = sort!(qr(Y', ColumnNorm()).p[1:k])
    return DCPinning(rows, Matrix(F.V[:, end-k+1:end]))
end

"""
    checkdcsubsystem(work::CanonicalWork, nodenames)

Refuse a circuit whose explicit direct current block has no solution.

The block is `[v; i]` against the transport rows and the blocks' zero
frequency rows. It is singular when a current the blocks leave undetermined
is not determined by anything else either -- a scattering block whose free
current direction does not cross a static flux component's boundary, so that
it cancels in every transport row and no equation sees it. That is §20's
`H ker(M) = 0` condition, asked of the matrix rather than of a taxonomy of
block types, so a short which does reach across a boundary passes and one
which does not is named.
"""
function checkdcsubsystem(work::CanonicalWork, nodenames = String[])
    isnothing(work.transport) && return nothing
    A = dcsubsystem(work)
    isempty(A) && return nothing
    n = size(A, 1)
    # a relative tolerance, which is what `rank` defaults to: the direct
    # current conductances are of order 1e-11 in these units, so an absolute
    # floor anywhere near `sqrt(eps())` calls every circuit singular
    if rank(A) < n
        br = work.blockrows
        free = isnothing(br) ? String[] :
            [d.name for d in br.descriptors if d.freecurrents > 0]
        detail = isempty(free) ?
            "no scattering block leaves a current undetermined, so this is a floating subnetwork with nowhere for its direct current to go" :
            lazy"the block(s) $(join(free, \", \")) leave a port current undetermined, and nothing else determines it: a free current direction which stays inside one static flux component cancels in that component's transport row"
        throw(ArgumentError(lazy"No direct current solution exists: the explicit direct current block is singular. $(detail)."))
    end
    return nothing
end

"""
    canonicaljacobian(Jint::SparseMatrixCSC, work::CanonicalWork)

The Jacobian in canonical coordinates, assembled from the internal one.

The flux block is `Jint` under the layout's permutation, which is a
symmetric reordering and nothing more. What is added is the explicit direct
current block and its couplings, in the same places the residual adds them:
the resistor current the average voltages drive into the zero frequency
nodal rows, the transport rows and the block currents they carry across a
component boundary, and each block's own zero frequency row in place of the
`i = 0` the stamp wrote. A pinned row is overwritten rather than added to,
for the same reason it is in the residual.

The matrix-free product is the reference this is checked against: for every
unit vector the two must agree exactly, which is a sharper test than any
finite difference of the residual would be.
"""
function canonicaljacobian(Jint::SparseMatrixCSC, work::CanonicalWork)
    L = work.layout
    p = L.perm
    n = length(p)
    nv = L.nvdc
    A = Jint[p, p]
    iszero(nv) && isnothing(work.transport) && return A

    t = work.transport
    br = work.blockrows
    # rows and columns of the direct current block, in canonical numbering
    I1, J1, V1 = Int[], Int[], Float64[]      # additions to the flux block
    I2, J2, V2 = Int[], Int[], Float64[]      # flux rows, voltage columns
    I3, J3, V3 = Int[], Int[], Float64[]      # voltage rows, flux columns
    I4, J4, V4 = Int[], Int[], Float64[]      # voltage rows, voltage columns

    # the resistor current the voltages drive into the zero frequency nodal
    # rows, and the transport rows themselves
    C = t.coupling
    for j in axes(C, 2), k in nzrange(C, j)
        push!(I2, L.nac + C.rowval[k]); push!(J2, j); push!(V2, C.nzval[k])
    end
    for j in axes(t.Ytr, 2), i in axes(t.Ytr, 1)
        iszero(t.Ytr[i,j]) && continue
        push!(I4, i); push!(J4, j); push!(V4, t.Ytr[i,j])
    end

    if !isnothing(br)
        # the block currents each component exchanges across its boundary
        for (c, ci, sgn) in br.transportterms
            (c in t.pinned) && continue
            push!(I3, c); push!(J3, ci); push!(V3, float(sgn))
        end
        # and each block's own row, replacing the `-i` already in `A`
        for (b, d) in enumerate(br.descriptors)
            ci = br.currentindex[b]
            sc, rc = br.signalcomponent[b], br.refcomponent[b]
            for p_ in eachindex(ci)
                push!(I1, ci[p_]); push!(J1, ci[p_]); push!(V1, 1.0)
                for q in eachindex(ci)
                    push!(I1, ci[p_]); push!(J1, ci[q])
                    push!(V1, -d.C0[p_,q])
                    w = d.B0[p_,q]*br.scale
                    iszero(sc[q]) ||
                        (push!(I2, ci[p_]); push!(J2, sc[q]); push!(V2, w))
                    iszero(rc[q]) ||
                        (push!(I2, ci[p_]); push!(J2, rc[q]); push!(V2, -w))
                end
            end
        end
    end

    A = A + sparse(I1, J1, V1, n, n)
    Bpv = sparse(I2, J2, V2, n, nv)
    Bvp = sparse(I3, J3, V3, nv, n)
    Yv  = sparse(I4, J4, V4, nv, nv)
    J = [A Bpv; Bvp Yv]

    # the pinned rows are replacements, so they are written last
    pn = work.pinning
    if !isnothing(pn)
        idx = work.dcindex
        for j in axes(pn.N, 2)
            r = idx[pn.rows[j]]
            J[r, :] .= 0.0
            for k in eachindex(idx)
                iszero(pn.N[k,j]) || (J[r, idx[k]] = pn.N[k,j])
            end
        end
        dropzeros!(J)
    end
    return J
end

"""
    canonicalfj(fjreal!, work::CanonicalWork, Jint, Jc)

Wrap an internal coordinate residual and Jacobian closure for a direct solve
method: the residual in canonical coordinates, and the canonical Jacobian
assembled into `Jc`.

`Jc` is filled rather than replaced, because the solver factorizes the
matrix it was handed. Its pattern is fixed -- the internal pattern under a
permutation plus the direct current block, none of which moves between
iterations -- and it is checked once rather than assumed.
"""
function canonicalfj(fjreal!, work::CanonicalWork, Jint, Jc)
    L = work.layout
    xint, Fint = work.xint, work.Fint
    return function (Fc, Jout, uc)
        scattercanonical!(xint, uc, L)
        fjreal!(isnothing(Fc) ? nothing : Fint,
            isnothing(Jout) ? nothing : Jint, xint)
        if !isnothing(Fc)
            gathercanonical!(Fc, Fint, L)
            addtransport!(Fc, work, uc)
        end
        if !isnothing(Jout)
            Jnew = canonicaljacobian(Jint, work)
            (Jnew.colptr == Jout.colptr && Jnew.rowval == Jout.rowval) ||
                error("the canonical Jacobian's pattern moved between iterations, which it must not: this is a bug in `canonicaljacobian`, not a property of the circuit.")
            copyto!(Jout.nzval, Jnew.nzval)
        end
        gathercanonical!(uc, xint, L)
        return nothing
    end
end
