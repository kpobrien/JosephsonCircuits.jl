# The canonical operators: the residual, the Jacobian vector product and the
# preconditioner of the internal system presented in canonical coordinates
# (the layout of harmonics/layout.jl), with the direct current block of
# harmonics/directcurrent.jl added on the window; and the canonical Jacobian
# as a plan, a fixed scatter of the internal Jacobian into a fixed pattern.

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
frequency mode, so it is applied where it was built, on the internal block
of the canonical vector through [`internalpart`](@ref), with no copy and no
permutation, rather than the preconditioner being rederived. That is
exact, and it keeps the two paths taking the same iterations.
"""
struct CanonicalPreconditioner{P,W<:CanonicalWork,F,D} <: AbstractWrappedPreconditioner
    inner::P
    work::W
    Yfact::F          # the direct current subsystem, factorized
    dcwork::Vector{Float64}
    device::D         # the same solve, resident where the state is
end

function CanonicalPreconditioner(inner, work::CanonicalWork)
    # The direct current subsystem is sparse and constant, with one unknown
    # per static flux component and one per block port current, and its
    # rows see no periodic state, so it is factorized once and solved
    # exactly. It has to be solved jointly: the transport rows carry the
    # block currents and the block rows carry the average voltages, and
    # solving only the transport half leaves the coupling to the Krylov
    # iteration, which stalls it. Where the state is not on the host, it is
    # solved there, at the positions found when the work was built (see
    # `DCFactorization`).
    A = dcsubsystem(work)
    F = kluordered(A)
    device = _onhost(work.xint) ? nothing : DCFactorization(A, F,
        work.dcindex, KernelAbstractions.get_backend(work.xint))
    return CanonicalPreconditioner(inner, work, F,
        zeros(Float64, length(work.dcindex)), device)
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
# rows. Solving `d` exactly, the flux block with an exact inner
# preconditioner, and dropping `Jpd` leaves the preconditioned operator the
# identity plus a part, from `Jpd`, whose square is zero: a Krylov
# iteration takes at most two steps per solve, whatever the number of
# direct current unknowns, and one when the residual has no direct current
# part. The exact triangular form, `d` first and `Jpd d` taken off the
# nodal rows before the inner solve, would save that one step at the cost
# of a pass over the nodal rows per application, so it is not here.
function applypreconditioner!(z::AbstractVector, pc::CanonicalPreconditioner,
        r::AbstractVector)
    L = pc.work.layout
    applypreconditioner!(internalpart(z, L), pc.inner, internalpart(r, L))
    return _solvedcblock!(z, pc, r)
end

# the direct current coordinates of `z`, solved exactly from `r`
function _solvedcblock!(z::AbstractVector, pc::CanonicalPreconditioner,
        r::AbstractVector)
    # where the state is not on the host, the same factors solve the
    # subsystem there (see `applydcsolve!`)
    if !isnothing(pc.device)
        applydcsolve!(z, r, pc.device)
        return z
    end
    idx = pc.work.dcindex
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

# The rest of the preconditioner interface forwards to the inner
# preconditioner through `AbstractWrappedPreconditioner`, escalation
# included: the defaults are inert, and a wrapper which did not forward
# them would silently turn escalation off. The inner preconditioner reads
# the internal block of the canonical vectors the harvest hooks carry, so
# forwarding them unchanged is right: the deflation which reads them wraps
# this one from outside, in canonical coordinates.
#
# `isexactpreconditioner` forwards as well, so this is exact when the inner
# is: what the wrapper then leaves, the dropped `Jpd` block, costs at most
# one more Arnoldi step per solve, when the residual has a direct current
# part, which is fewer than the products and exact solves a deflation of
# it would cost at every Newton step. A recycling wrapper outside this one
# therefore clears its pair after an escalation rather than learning `Jpd`.
innerpreconditioner(pc::CanonicalPreconditioner) = pc.inner

# The canonical Jacobian, as a plan.
#
# Every term the direct current block contributes is constant: the transport
# rows, the coupling `G0 P`, the block pencils `B0` and `C0`, the boundary
# currents and the reference rows do not depend on the periodic state. What
# moves between Newton iterations is the internal Jacobian, and only its
# values, since its pattern is fixed. So the whole assembly is a fixed
# pattern, a fixed scatter of the internal values into it, and a fixed list
# of constant additions, found once.

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

"""
    canonicaljacobianplan(Jint::SparseMatrixCSC, work::CanonicalWork)

Build the [`CanonicalJacobianPlan`](@ref) for an internal Jacobian pattern.

Only the pattern of `Jint` is read, not its values, so the plan is valid for
every point the solve visits. It stops being valid if that pattern moves,
which it must not. The direct current block's entries are those of its
matrix form, [`dcupdate`](@ref), which is read off the residual, so the
Jacobian and the residual agree by construction.
"""
function canonicaljacobianplan(Jint::SparseMatrixCSC, work::CanonicalWork)
    L = work.layout
    up = dcupdate(work)
    window = windowindices(L)

    # a row the block writes over, a transport row or a reference, keeps
    # nothing of the internal Jacobian
    written = Set(window[k] for k in eachindex(up.keep) if iszero(up.keep[k]))

    # where each stored entry of the internal Jacobian lands: where it is,
    # since the internal state is the first block of the canonical one
    di = rowvals(Jint)
    dj = Vector{Int}(undef, nnz(Jint))
    for col in axes(Jint, 2), k in nzrange(Jint, col)
        dj[k] = col
    end
    keepd = [k for k in eachindex(di) if !(di[k] in written)]

    # the block's constant entries: its matrix, stored by rows on the
    # window, in the canonical numbering
    fi, fj = Int[], Int[]
    for i in eachindex(window), k in up.rowptr[i]:(up.rowptr[i+1] - 1)
        push!(fi, window[i]); push!(fj, window[up.colval[k]])
    end

    N = canonicaldim(L)
    J = sparse(vcat(di[keepd], fi), vcat(dj[keepd], fj),
        ones(Float64, length(keepd) + length(fi)), N, N)
    source = zeros(Int, length(di))
    for k in keepd
        source[k] = nzposition(J, di[k], dj[k])
    end
    fixedindex = Int[nzposition(J, fi[k], fj[k]) for k in eachindex(fi)]
    return CanonicalJacobianPlan(J, source, fixedindex, up.nzval)
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
    canonicalfj(fjreal!, work::CanonicalWork, Jint, plan)

Wrap an internal coordinate residual and Jacobian closure for a direct solve
method: the residual in canonical coordinates, and the canonical Jacobian
filled through `plan`.

The matrix is filled rather than replaced, because the solver factorizes the
one it was handed. Its pattern is fixed -- the internal pattern plus the
direct current block, neither of which moves between iterations -- and the
plan is built from it once.
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
