# =====================================================================
# Transposed gather maps for the vector-Jacobian product.
#
# Device-clean transposed maps for the vjp.
#
# The forward maps of a NonlinearTermPlan are gathers: one work item owns
# one output slot and reads only precomputed index maps, so there is no
# scatter, no atomic and no write conflict, and they run unchanged on a
# device. Their naive transposes are scatters and lose that property.
#
# This builds the transposed *gather* maps once, at plan time, so the
# transposed products keep the same one-item-per-output structure. The
# build is a counting sort over the same index arrays and costs nothing at
# run time.
#
# The transposes of the two transforms need no new kernel at all. Writing
# F = applyfft! and G = applyifft! (which is the unnormalized inverse, so
# G = inv(F) on the packing F produces), the adjoints with respect to the
# real inner product are
#
#     adjoint(G)(u) = N * s .* F(u)
#     adjoint(F)(P) = (1/N) * G(P ./ s)
#
# where s is the conjugate multiplicity of a stored bin (see gtscale) and
# N = prod(Nt). The second follows from the first: on the subspace F maps
# into, adjoint(F) = inv(adjoint(G)), and a covector's components on the
# imaginary parts of the self-conjugate bins are never seen, because F
# writes zero there. Composing the two through the diagonal cosine, the
# factors of N cancel exactly:
#
#     adjoint(G) . D . adjoint(F) (P) = s .* F( cos .* G(P ./ s) )
#
# so the vjp needs only the existing forward transform plans, one division
# by s folded into the first kernel and one multiplication by s folded
# into the coefficients of the second. The host-only applyffttranspose!,
# with its zero padded work array and its separate complex transform plan,
# is not needed.
# =====================================================================

"""
    NonlinearTermTransposePlan

The transposed gather maps of a [`NonlinearTermPlan`](@ref), and the two
kernels which apply them. Built once by
[`plannonlineartermtranspose`](@ref). See [`hbvjp!`](@ref).

# Fields
- `tbptr`, `tbnode`, `tbcoef`: per frequency domain slot, the segment of
  complex node indices which gather that slot in the backward map, with
  the same coefficients. The transpose of the Josephson contribution.
- `tfptr`, `tfslot`, `tfcoef`, `tfimag`: per real slot of the unknowns,
  the segment of frequency domain slots which read it in the forward map.
  `tfcoef` carries the incidence entry with the `FCONJ` sign and the
  conjugate multiplicity `s` already folded in; `tfimag` says whether the
  entry reads the imaginary rather than the real part of the slot.
- `ktptr`, `ktrow`, `ktcoef`: the transpose of the real form of the linear
  term, gathered by output index.
- `gtscale`: the conjugate multiplicity of each frequency domain slot, one
  for a bin whose conjugate partner is also stored (the zero and Nyquist
  frequencies of the first dimension, the only one the real transform
  truncates) and two otherwise.
"""
struct NonlinearTermTransposePlan{VI,VT,VF,KB,KF,B}
    tbptr::VI
    tbnode::VI
    tbcoef::VT
    tfptr::VI
    tfslot::VI
    tfcoef::VT
    tfimag::VF
    ktptr::VI
    ktrow::VI
    ktcoef::VT
    gtscale::VT
    backwardtranspose!::KB
    forwardtranspose!::KF
    backend::B
    nslots::Int
    rdim::Int
end

"""
    conjugatemultiplicity(fd, td)

The number of harmonics of the real time domain signal each stored bin of
the frequency domain array represents: one when the bin's conjugate
partner is also stored, two otherwise.

The real transform truncates only the first dimension, so the partner of
a bin is stored exactly when its first index is the zero or the Nyquist
frequency. This is the factor which distinguishes the adjoint of
[`applyifft!`](@ref) from [`applyfft!`](@ref), and getting it wrong gives
a vjp which is correct on self-conjugate modes and wrong by a factor of
two on every other one.
"""
function conjugatemultiplicity(fd::AbstractArray, td::AbstractArray)
    s = ones(Float64, length(fd))
    isempty(fd) && return s
    Nw1 = size(fd, 1)
    hasnyquist = iseven(size(td, 1))
    @inbounds for (q, I) in enumerate(CartesianIndices(fd))
        i1 = I[1]
        s[q] = (i1 == 1 || (hasnyquist && i1 == Nw1)) ? 1.0 : 2.0
    end
    return s
end

# counting sort of (key, payload...) triples into a compressed gather map
function _buildgather(keys::Vector{Int}, nkeys::Integer, Ti::Type)
    counts = zeros(Int, nkeys)
    @inbounds for k in keys
        counts[k] += 1
    end
    ptr = Vector{Ti}(undef, nkeys+1)
    ptr[1] = 1
    @inbounds for i in 1:nkeys
        ptr[i+1] = ptr[i] + counts[i]
    end
    cursor = Vector{Int}(undef, nkeys)
    @inbounds for i in 1:nkeys
        cursor[i] = ptr[i]
    end
    perm = Vector{Int}(undef, length(keys))
    @inbounds for (e, k) in enumerate(keys)
        perm[cursor[k]] = e
        cursor[k] += 1
    end
    return ptr, perm
end

"""
    plannonlineartermtranspose(plan, modelayout, fd, td; backend = CPU())

Build the [`NonlinearTermTransposePlan`](@ref) of a
[`NonlinearTermPlan`](@ref).
"""
function plannonlineartermtranspose(plan::NonlinearTermPlan, modelayout,
        fd::AbstractArray, td::AbstractArray; backend = CPU())
    Ti = Int32
    nslots = plan.nslots
    ncomplex = plan.ncomplex
    rdim = modelayout.rdim
    gtscale = conjugatemultiplicity(fd, td)

    # ---- transpose of the Josephson part of the backward map ----
    # forward: node k gathers slots bsrc[t]. transpose: slot q gathers the
    # nodes k which named it.
    bkeys = Int[]; bnodes = Int[]; bcoefs = Float64[]
    @inbounds for k in 1:ncomplex, t in plan.bptr[k]:plan.bptr[k+1]-1
        push!(bkeys, Int(plan.bsrc[t]))
        push!(bnodes, k)
        push!(bcoefs, plan.bcoef[t])
    end
    tbptr, bperm = _buildgather(bkeys, nslots, Ti)
    tbnode = Ti[bnodes[e] for e in bperm]
    tbcoef = Float64[bcoefs[e] for e in bperm]

    # ---- transpose of the forward map ----
    # forward: slot q reads real slot n1[q] (and n1[q]+1 when the mode is
    # not self conjugate), likewise n2. transpose: real slot p gathers the
    # slots which read it. The FCONJ sign and the conjugate multiplicity
    # are folded into the coefficient here so the kernel is a plain gather.
    fkeys = Int[]; fslots = Int[]; fcoefs = Float64[]; fimag = Int32[]
    @inbounds for q in 1:nslots
        f = plan.flags[q]
        wide = (f & FWIDE) != Int32(0)
        sgn = (f & FCONJ) != Int32(0) ? -1.0 : 1.0
        sq = gtscale[q]
        for (nv, sv) in ((plan.n1, plan.s1), (plan.n2, plan.s2))
            p = nv[q]
            iszero(p) && continue
            c = sv[q]
            push!(fkeys, Int(p)); push!(fslots, q)
            push!(fcoefs, c*sq);  push!(fimag, Int32(0))
            if wide
                push!(fkeys, Int(p)+1); push!(fslots, q)
                push!(fcoefs, sgn*c*sq); push!(fimag, Int32(1))
            end
        end
    end
    tfptr, fperm = _buildgather(fkeys, rdim, Ti)
    tfslot = Ti[fslots[e] for e in fperm]
    tfcoef = Float64[fcoefs[e] for e in fperm]
    tfimag = Int32[fimag[e] for e in fperm]

    # ---- transpose of the linear term ----
    kkeys = Int[]; krows = Int[]; kcoefs = Float64[]
    if !isempty(plan.kptr)
        @inbounds for p in 1:(length(plan.kptr)-1), t in plan.kptr[p]:plan.kptr[p+1]-1
            push!(kkeys, Int(plan.kidx[t])); push!(krows, p)
            push!(kcoefs, plan.kcoef[t])
        end
    end
    ktptr, kperm = _buildgather(kkeys, rdim, Ti)
    ktrow = Ti[krows[e] for e in kperm]
    ktcoef = Float64[kcoefs[e] for e in kperm]

    bk = backwardjosephsontransposekernel!(backend, 64, max(nslots, 1))
    fk = forwardtransposekernel!(backend, 64, max(rdim, 1))

    return NonlinearTermTransposePlan(
        tobackend(backend, tbptr), tobackend(backend, tbnode),
        tobackend(backend, tbcoef),
        tobackend(backend, tfptr), tobackend(backend, tfslot),
        tobackend(backend, tfcoef), tobackend(backend, tfimag),
        tobackend(backend, ktptr), tobackend(backend, ktrow),
        tobackend(backend, ktcoef), tobackend(backend, gtscale),
        bk, fk, backend, nslots, rdim)
end

"""
    backwardjosephsontransposekernel!(P, tbptr, tbnode, tbcoef, w, lptr,
        lwide, gtscale)

Transpose of the Josephson contribution of the backward map, one work item
per frequency domain slot. Each item gathers the nodes which read its slot
and writes its own entry of `P`, so the writes are disjoint and no atomic
is needed.

Divides by the conjugate multiplicity in the same pass, which is the first
half of the adjoint of [`applyfft!`](@ref).

The imaginary part of a self-conjugate bin is deliberately *not* zeroed
here. The pairing never sees it, so zeroing looks harmless, and in one
dimension it is. In more than one it is wrong: the real part has to be
taken after the transform along the remaining dimensions, not before, and
the two differ as soon as the first-dimension zero-frequency slice has
more than one element. The complex-to-real transform already discards
exactly the right component, so the correct action here is none.
"""
@kernel function backwardjosephsontransposekernel!(P, @Const(tbptr),
        @Const(tbnode), @Const(tbcoef), @Const(w), @Const(lptr),
        @Const(lwide), @Const(gtscale))
    q = @index(Global)
    T = eltype(w)
    @inbounds begin
        accre = zero(T); accim = zero(T)
        for t in tbptr[q]:tbptr[q+1]-1
            k = tbnode[t]
            c = tbcoef[t]
            p = lptr[k]
            accre += c*w[p]
            if lwide[k] == Int32(1)
                accim += c*w[p+1]
            end
        end
        s = gtscale[q]
        P[q] = Complex(accre/s, accim/s)
    end
end

"""
    forwardtransposekernel!(out, tfptr, tfslot, tfcoef, tfimag, Q, ktptr,
        ktrow, ktcoef, w)

Transpose of the forward map plus the transposed linear term, one work item
per real slot of the unknowns. Each item gathers the frequency domain slots
which read it and the entries of the transposed linear term which land on
it, and writes the single slot it owns.
"""
@kernel function forwardtransposekernel!(out, @Const(tfptr), @Const(tfslot),
        @Const(tfcoef), @Const(tfimag), @Const(Q), @Const(ktptr),
        @Const(ktrow), @Const(ktcoef), @Const(w))
    p = @index(Global)
    T = eltype(out)
    @inbounds begin
        acc = zero(T)
        for t in tfptr[p]:tfptr[p+1]-1
            v = Q[tfslot[t]]
            acc += tfcoef[t]*(tfimag[t] == Int32(1) ? imag(v) : real(v))
        end
        for t in ktptr[p]:ktptr[p+1]-1
            acc += ktcoef[t]*w[ktrow[t]]
        end
        out[p] = acc
    end
end

"""
    applybackwardjosephsontranspose!(P, tplan, plan, w)

Apply the transposed Josephson backward map, scaled by the inverse
conjugate multiplicity.
"""
function applybackwardjosephsontranspose!(P::AbstractArray,
        tplan::NonlinearTermTransposePlan, plan::NonlinearTermPlan,
        w::AbstractVector)
    tplan.backwardtranspose!(P, tplan.tbptr, tplan.tbnode, tplan.tbcoef, w,
        plan.lptr, plan.lwide, tplan.gtscale)
    KernelAbstractions.synchronize(tplan.backend)
    return P
end

"""
    applyforwardtranspose!(out, tplan, Q, w)

Apply the transposed forward map and add the transposed linear term.
Overwrites `out`.
"""
function applyforwardtranspose!(out::AbstractVector,
        tplan::NonlinearTermTransposePlan, Q::AbstractArray,
        w::AbstractVector)
    tplan.forwardtranspose!(out, tplan.tfptr, tplan.tfslot, tplan.tfcoef,
        tplan.tfimag, Q, tplan.ktptr, tplan.ktrow, tplan.ktcoef, w)
    KernelAbstractions.synchronize(tplan.backend)
    return out
end
