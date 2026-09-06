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

# The rest of the preconditioner interface forwards to the inner
# preconditioner through `AbstractWrappedPreconditioner`, escalation
# included: the defaults are inert, and a wrapper which did not forward
# them would silently turn escalation off. The inner preconditioner reads
# the internal block of the canonical vectors the harvest hooks carry, so
# forwarding them unchanged is right: the deflation which reads them wraps
# this one from outside, in canonical coordinates.
innerpreconditioner(pc::CanonicalPreconditioner) = pc.inner

# Exact when the inner is: the error the wrapper then leaves, the dropped
# `Jpd` block, has rank at most the number of direct current unknowns and
# costs that many Arnoldi steps per solve, which is fewer than the products
# and exact solves a deflation of it would cost at every Newton step. A
# recycling wrapper outside this one therefore clears its pair after an
# escalation rather than learning `Jpd`.
isexactpreconditioner(pc::CanonicalPreconditioner) =
    isexactpreconditioner(pc.inner)

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
