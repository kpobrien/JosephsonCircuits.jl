# =========================================================================
# Sensitivities of the scattering parameters to component values.
# =========================================================================

"""
    HBOperatingPoint(sys, x, jacobian, modelayout, Nnodal, Lmean, wmodes,
        Amna, mnaindices, coupledbranches, Nmodes, Nnodes)

The converged pump operating point of [`hbnlsolve`](@ref) together with
everything needed to propagate a component perturbation through it: the
[`HBSystem`](@ref) evaluation object, the converged augmented state, the
exact Jacobian of the equivalent real system assembled there, and the
scaled matrices and layout of the augmented system.

Requested with `returnoperatingpoint = true`. The Jacobian is the exact
Jacobian of the equivalent real system, assembled with
[`assemblerealjacobian!`](@ref), rather than the complex holomorphic
Jacobian of the `:quasinewton` method, which is only an approximation: the
harmonic balance residual is not complex differentiable, so the implicit
function theorem does not hold with the holomorphic Jacobian, while in the
real representation it applies directly.
"""
struct HBOperatingPoint
    # the HBSystem evaluation object and the ModeLayout; untyped because
    # hbsolve.jl is included before hbsystem.jl and realcomplexconv.jl
    sys
    x::Vector{Complex{Float64}}
    jacobian::SparseMatrixCSC{Float64,Int}
    modelayout
    Nnodal::Int
    Lmean::Complex{Float64}
    wmodes::Vector{Float64}
    Amna::SparseMatrixCSC
    mnaindices::Vector{Int}
    coupledbranches::Vector{Int}
    Nmodes::Int
    Nnodes::Int
    # The explicit direct current block, when the circuit has one: the
    # canonical work, the canonical state the solve actually converged, and
    # the canonical Jacobian there. `jacobian` above stays the harmonic one,
    # because everything which differentiates the harmonic system alone
    # reads it; what differentiates the whole system reads these.
    dc
end

# an operating point of a circuit with no direct current block, which is
# every one taken before the block existed
HBOperatingPoint(sys, x, jacobian, modelayout, Nnodal, Lmean, wmodes, Amna,
    mnaindices, coupledbranches, Nmodes, Nnodes) =
    HBOperatingPoint(sys, x, jacobian, modelayout, Nnodal, Lmean, wmodes,
        Amna, mnaindices, coupledbranches, Nmodes, Nnodes, nothing)

"""
    DCOperatingPoint

The explicit direct current block at a converged point.

# Fields
- `work`: the [`CanonicalWork`](@ref) carrying the layout, the transport
  rows and the blocks' zero frequency rows.
- `u`: the converged canonical state, `[phiac | phidc | vdc]`.
- `jacobian`: the canonical Jacobian there, which is the one the implicit
  function theorem applies to when the block is active.
- `plan`: the [`CanonicalJacobianPlan`](@ref) that filled it.
"""
struct DCOperatingPoint{W,P}
    work::W
    u::Vector{Float64}
    jacobian::SparseMatrixCSC{Float64,Int}
    plan::P
end

"""
    dcresidualsensitivity(dc::DCOperatingPoint, psc, nm, scale,
        sensitivityindices, alphas)

The direct current block's own contribution to the residual sensitivity, in
canonical coordinates.

The canonical residual is `D G F(S u) + M u + s c`, so its derivative with
respect to a component value is the harmonic one gathered and masked, plus
`(dM/dr) u`. Only the conductance depends on a lumped component value, and
it enters twice: as the transport rows `Y = P' G0 P` and as the coupling
`G0 P` into the zero frequency nodal rows. A capacitor, an inductor and a
junction are open circuits, a short and a short at zero frequency, none of
which carries a conductance, so they contribute nothing here -- which is
correct and is why the harmonic rows alone were right until a resistor
carried direct current.

The perturbation is relative, as it is everywhere in this file: `G0` is
proportional to `1/R`, so a relative change in `R` scales the whole stamp by
`-1`.
"""
function dcresidualsensitivity(dc::DCOperatingPoint, psc::CompiledCircuit,
        nm::CircuitMatrices, scale, sensitivityindices,
        alphas::AbstractVector = ones(Complex{Float64},
            length(sensitivityindices)))
    w = dc.work
    L = w.layout
    plan = w.transport.plan
    n = canonicaldim(L)
    # the average voltage of every node, ground dropped, from the component
    # voltages the solve carries
    nodev = plan.lift * dc.u[voltagerange(L)]
    # the solver scale, which the whole system's rows are multiplied by and
    # which is not the mean inductance: it is `Z0/w0` (see
    # [`calcsolverscale`](@ref)), and using the wrong one leaves the block's
    # own rows and their derivative scaled differently
    gscale = real(scale)
    I, J, V = Int[], Int[], Float64[]
    for (k, idx) in enumerate(sensitivityindices)
        psc.componenttypes[idx] === :R || continue
        r = nm.vvn[idx]
        (r isa Number && isfinite(r) && !iszero(r)) || continue
        # the conductance the assembly gives this resistor at zero
        # frequency, which is its own value scaled the way every other
        # entry of the system is; a relative perturbation scales it by -1
        g = -real(alphas[k])*gscale/real(r)
        n1, n2 = psc.nodeindices[1, idx], psc.nodeindices[2, idx]
        v1 = n1 > 1 ? nodev[n1-1] : zero(eltype(nodev))
        v2 = n2 > 1 ? nodev[n2-1] : zero(eltype(nodev))
        cur = g*(v1 - v2)
        iszero(cur) && continue
        # the current it drives into the zero frequency nodal rows
        n1 > 1 && (push!(I, L.nac + n1 - 1); push!(J, k); push!(V, cur))
        n2 > 1 && (push!(I, L.nac + n2 - 1); push!(J, k); push!(V, -cur))
        # and into the transport rows, which are the component sums of
        # those; a resistor inside one component cancels there, exactly as
        # an inductor branch does
        c1 = plan.componentof[n1]
        c2 = plan.componentof[n2]
        c1 == c2 && continue
        iszero(c1) || (push!(I, L.nac + L.ndc + c1); push!(J, k); push!(V, cur))
        iszero(c2) || (push!(I, L.nac + L.ndc + c2); push!(J, k); push!(V, -cur))
    end
    return sparse(I, J, V, n, length(sensitivityindices))
end

"""
    sensitivityjacobian(op::HBOperatingPoint)

The Jacobian the implicit function theorem applies to at `op`: the canonical
one when an explicit direct current block is active, and the harmonic one
otherwise.

The forward and the reverse contraction both solve against this, forward
through it and reverse through its transpose, so naming it once is what
keeps the two orders solving the same system.
"""
sensitivityjacobian(op::HBOperatingPoint) =
    isnothing(op.dc) ? op.jacobian : op.dc.jacobian

"""
    sensitivitydim(op::HBOperatingPoint)

The dimension of the space the residual derivatives and the adjoint
covectors live in: canonical with a direct current block, and the real
representation of the augmented harmonic state without one.
"""
sensitivitydim(op::HBOperatingPoint) =
    isnothing(op.dc) ? size(op.jacobian, 1) : canonicaldim(op.dc.work.layout)

"""
    componentlookups(mnaindices, coupledbranches, Ljb)

Constant time lookups for [`componentstamp`](@ref), built once per stamp
table rather than searched per component: the ordinal of a promoted
resistor within `mnaindices`, the set of mutually coupled branches, and the
ordinal of a junction branch within `Ljb.nzind`. Without these the
classification repeats linear searches per component, which becomes
quadratic over a large sensitivity set.
"""
function componentlookups(mnaindices, coupledbranches, Ljb)
    return (
        coupled = Set(coupledbranches),
        junctionordinal = Dict(b => j for (j, b) in enumerate(Ljb.nzind)))
end

"""
    componentstamp(idx::Integer, psc::CompiledCircuit, cg::CircuitGraph,
        nm::CircuitMatrices, lookups, Nmodes::Integer, Nnodes::Integer)

Classify the component at index `idx` for sensitivity analysis and build its
raw one-component matrix, without any solver scaling, negative frequency
conjugation, or padding, which the callers apply for their own grids. The
component matrices are built with the same functions which build the system
matrices, [`calcCn`](@ref), [`calcGn`](@ref), [`calcLb`](@ref) and
[`calcinvLn`](@ref), applied to the single component, so the node and mode
conventions agree by construction. Returns one of

- `(:C, M)`: the component's capacitance matrix,
- `(:G, M)`: the component's conductance matrix, for a resistor which is not
    promoted by the modified nodal analysis formulation,
- `(:Lj, j)`: the ordinal `j` of a Josephson junction within the junction
    branch vector `nm.Ljb`,
- `(:invL, M)`: the component's inverse inductance matrix.

This is the single definition of which components are supported: `:C`, `:L`,
`:R` and `:Lj` with numeric values. Mutually coupled inductors and
components with symbolic (frequency dependent) values throw, with the same
message from both the fixed operating point stamps
([`calcsensitivitystamps`](@ref)) and the residual derivatives
([`calcresidualsensitivity`](@ref)). `lookups` are the constant time
tables of [`componentlookups`](@ref).
"""
function componentstamp(idx::Integer, psc::CompiledCircuit,
    cg::CircuitGraph, nm::CircuitMatrices, lookups,
    Nmodes::Integer, Nnodes::Integer)

    componenttypes = psc.componenttypes
    nodeindices = psc.nodeindices
    vvn = nm.vvn
    componenttype = componenttypes[idx]
    value = vvn[idx]
    if !(value isa Number)
        throw(ArgumentError(lazy"Sensitivities require a numeric component value, but the value of $(psc.componentnames[idx]) is $(value). Components with symbolic frequency dependent values are not supported."))
    end
    n1 = nodeindices[1, idx]
    n2 = nodeindices[2, idx]
    if componenttype == :C
        return (:C, calcCn(componenttypes[[idx]], nodeindices[:,[idx]],
            vvn[[idx]], Nmodes, Nnodes))
    elseif componenttype == :R
        return (:G, calcGn(componenttypes[[idx]], nodeindices[:,[idx]],
            vvn[[idx]], Nmodes, Nnodes))
    elseif componenttype == :L
        b = cg.edge2indexdict[(n1, n2)]
        if b in lookups.coupled
            throw(ArgumentError(lazy"Sensitivities are not supported for the mutually coupled inductor $(psc.componentnames[idx])."))
        end
        Lb = calcLb(componenttypes[[idx]], nodeindices[:,[idx]],
            vvn[[idx]], cg.edge2indexdict, 1, cg.Nbranches)
        return (:invL, calcinvLn(Lb, cg.Rbn, Nmodes))
    elseif componenttype == :Lj
        b = cg.edge2indexdict[(n1, n2)]
        j = get(lookups.junctionordinal, b, nothing)
        if isnothing(j)
            throw(ArgumentError(lazy"The Josephson junction $(psc.componentnames[idx]) was not found in the branch inductance vector."))
        end
        return (:Lj, j)
    else
        throw(ArgumentError(lazy"Sensitivities are only supported for C, L, R, and Lj components, not $(componenttype), the type of $(psc.componentnames[idx])."))
    end
end

"""
    calcresidualsensitivity(op::HBOperatingPoint, psc, cg, nm,
        sensitivityindices)

Calculate the derivative of the harmonic balance residual with respect to a
relative (logarithmic) perturbation of each component value, at the
operating point. Combined with the implicit function theorem applied to
`F(x, r) = 0` in the equivalent real representation,

    dx/dr = -inv(J)*(dF/dr),

with `J` the exact real Jacobian retained by [`hbnlsolve`](@ref), this gives
the derivative of the operating point itself
(see [`calcnodefluxsensitivity`](@ref)). Returns a sparse matrix whose
columns are `dF/dr` for each component, in the real representation of the
augmented residual: each component touches only its own rows (its nodes and
modes, a promoted resistor's constitutive rows, or a junction branch's
Kirchhoff rows), so the storage scales with the touched entries rather than
with `Nstate*Ncomponents`.

The residual is affine in `C`, `1/R` and `1/L`, so those parameter
derivatives are that component's own contribution to the linear term applied
to the converged state, with a sign, built from the shared classification of
[`componentstamp`](@ref) with the same nondimensionalization and negative
frequency conjugation the solver applies. The Josephson junction term is
the residual's own sine contribution restricted to that junction.

Note that the auxiliary branch currents of the modified nodal analysis
formulation are scaled by the solver scale (see [`calcsolverscale`](@ref)),
which itself depends on the port impedances, so for a port resistor the
auxiliary rows of the returned derivative differ from a finite difference of
a re-solve by that change of normalization. The node flux rows, which are
the physical quantity and the only rows the linearized system depends on,
are unaffected.
"""
function calcresidualsensitivity(op::HBOperatingPoint,
    psc::CompiledCircuit, cg::CircuitGraph, nm::CircuitMatrices,
    sensitivityindices,
    alphas::AbstractVector = ones(Complex{Float64},
        length(sensitivityindices)))

    if isnothing(op.jacobian)
        throw(ArgumentError("The operating point does not contain a Jacobian."))
    end
    # evaluate at the operating point regardless of what the shared
    # evaluation object was last used for: the Josephson term below reads
    # the cached time domain branch fluxes.
    setpoint!(op.sys, op.x)
    Ntot = length(op.x)
    Nmodes = op.Nmodes
    Nnodes = op.Nnodes
    wmodes = op.wmodes

    # the component's own matrix, nondimensionalized, negative frequency
    # conjugated and padded exactly as hbnlsolve does with the full matrices
    function scaledpadded(M, alpha = one(Complex{Float64}))
        Ms = SparseMatrixCSC{Complex{Float64},Int}(copy(M))
        # the design parameter rescale is a direction in component value
        # space, so it multiplies the stored value before the negative
        # frequency conjugation, exactly as in `reparameterize`
        isone(alpha) || rmul!(Ms, alpha)
        conjnegfreq!(Ms, wmodes)
        rmul!(Ms, op.Lmean)
        return mnapadto(Ms, Ntot)
    end

    # The residual derivatives are sparse by nature: each component touches
    # only its own rows of the state (its nodes and modes, one promoted
    # resistor's constitutive rows, or one junction branch's Kirchhoff
    # rows), so they are accumulated as triplets of the real representation
    # and returned as a sparse matrix. The dense alternative peaks at
    # O(Nstate*Ncomponents) memory, which is exactly the many component
    # regime the reverse contraction order exists for. Duplicate triplets
    # (a component matrix with several entries in one row) are summed by
    # sparse, and the real representation is linear, so per entry
    # realification distributes over the sums.
    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)
    # the offset of complex entry r in the real representation (one real
    # row for the self conjugate modes, a real and an imaginary row
    # otherwise), matching complex_to_real! with the default scale.
    realindexmap = zeros(Int, Ntot)
    k = 1
    for j in eachindex(realindexmap)
        realindexmap[j] = k
        k += isrealmode[(j-1) % nmd + 1] ? 1 : 2
    end
    Nreal = realdim(Ntot, isrealmode)
    Ir = Int[]; Jc = Int[]; Vr = Float64[]
    function pushentry!(r, comp, v)
        kr = realindexmap[r]
        push!(Ir, kr); push!(Jc, comp); push!(Vr, real(v))
        if !isrealmode[(r-1) % nmd + 1]
            push!(Ir, kr+1); push!(Jc, comp); push!(Vr, imag(v))
        end
        return nothing
    end

    # work vector and branch support for the Josephson terms, which are
    # evaluated through the full residual machinery and then read off on
    # the Kirchhoff rows of that junction's branch only.
    cwork = Complex{Float64}[]
    branchsupport = Vector{Tuple{Int,Int}}[]

    x = op.x
    Nm = length(wmodes)
    lookups = componentlookups(op.mnaindices, op.coupledbranches,
        op.sys.Ljb)
    for (comp, idx) in enumerate(sensitivityindices)
        kind, info = componentstamp(idx, psc, cg, nm, lookups,
            Nmodes, Nnodes)
        if kind == :C || kind == :G || kind == :invL
            # dF_comp = c * Ms * Diagonal(w.^power) * x, accumulated per
            # stored entry of the component's own (tiny) matrix.
            Ms = scaledpadded(info, alphas[comp])
            c, power = kind == :C ? (-1.0 + 0im, 2) :
                kind == :G ? (0.0 - 1im, 1) : (-1.0 + 0im, 0)
            rows = rowvals(Ms)
            vals = nonzeros(Ms)
            for col in axes(Ms, 2)
                xc = x[col]
                iszero(xc) && continue
                w = wmodes[(col-1) % Nm + 1]
                scale = power == 0 ? c*xc : power == 1 ? c*w*xc : c*w^2*xc
                for pp in nzrange(Ms, col)
                    pushentry!(rows[pp], comp, vals[pp]*scale)
                end
            end
        else # :Lj
            if isempty(cwork)
                cwork = zeros(Complex{Float64}, Ntot)
                branchsupport = branchnodesandsigns(op.sys.Rbnm, Nmodes,
                    size(op.sys.Rbnm, 1) ÷ Nmodes)
            end
            residualjosephsonterm!(cwork, op.sys, info)
            branch = op.sys.Ljb.nzind[info]
            for (node, _) in branchsupport[branch]
                for m in 1:Nmodes
                    r = (node-1)*Nmodes + m
                    pushentry!(r, comp, -cwork[r]*alphas[comp])
                end
            end
        end
    end
    harmonic = sparse(Ir, Jc, Vr, Nreal, length(sensitivityindices))
    isnothing(op.dc) && return harmonic

    # In canonical coordinates the residual is `D G F(S u) + M u`, so its
    # parameter derivative is this gathered, masked by the rows the block
    # replaces rather than adds to, plus the block's own dependence on the
    # component values. The gather is a permutation of the flux rows and
    # writes nothing into the voltage rows, which is where the block's part
    # lands.
    L = op.dc.work.layout
    N = canonicaldim(L)
    dccols = dcresidualsensitivity(op.dc, psc, nm, op.Lmean,
        sensitivityindices, alphas)
    Ic, Jc2, Vc = Int[], Int[], Float64[]
    keep = dckeep(op.dc.work)
    col = zeros(Float64, Nreal)
    gathered = zeros(Float64, N)
    for k in axes(harmonic, 2)
        fill!(col, 0.0)
        col .= view(harmonic, :, k)
        fill!(gathered, 0.0)
        gathercanonical!(gathered, col, L)
        for i in eachindex(gathered)
            v = gathered[i]*keep[i]
            iszero(v) && continue
            push!(Ic, i); push!(Jc2, k); push!(Vc, v)
        end
    end
    return sparse(Ic, Jc2, Vc, N, length(sensitivityindices)) + dccols
end

"""
    dckeep(work::CanonicalWork)

The diagonal which is zero on the rows the direct current block replaces
rather than adds to, over the whole canonical vector.

Read off the block by probing it, so it cannot disagree with the residual it
describes.
"""
function dckeep(work::CanonicalWork)
    L = work.layout
    keep = ones(Float64, canonicaldim(L))
    up = dcupdate(work)
    isnothing(up) && return keep
    copyto!(view(keep, dcwindow(L)), Array(up.keep))
    return keep
end

"""
    dcvoltagesensitivity(op::HBOperatingPoint, dFr::AbstractMatrix;
        factorization = KLUfactorization())

The derivative of the average node voltages with respect to a relative
perturbation of each component value, in volts, indexed by node with ground
dropped as [`hbnlsolve`](@ref) reports them.

The same solve as [`calcnodefluxsensitivity`](@ref) and the other half of
its answer: that returns the node fluxes, which are the periodic part, and
this returns the average voltages, which are the direct current part and
are unknowns of the canonical system rather than of the harmonic one.
`nothing` for a circuit with no direct current block.
"""
function dcvoltagesensitivity(op::HBOperatingPoint, dFr::AbstractMatrix;
        factorization = KLUfactorization())
    isnothing(op.dc) && return nothing
    L = op.dc.work.layout
    lift = op.dc.work.transport.plan.lift
    cache = FactorizationCache()
    tryfactorize!(cache, factorization, op.dc.jacobian)
    rhs = zeros(Float64, canonicaldim(L))
    duc = zeros(Float64, canonicaldim(L))
    dv = zeros(Float64, size(lift, 1), size(dFr, 2))
    for k in axes(dFr, 2)
        rhs .= view(dFr, :, k)
        trysolve!(duc, cache.factorization, rhs)
        rmul!(duc, -1)
        dv[:, k] .= phi0 .* (lift*view(duc, voltagerange(L)))
    end
    return dv
end

"""
    calcnodefluxsensitivity(op::HBOperatingPoint, dFr::AbstractMatrix;
        factorization = KLUfactorization())

Solve `dx/dr = -inv(J)*(dF/dr)` for the residual sensitivities `dFr` of
[`calcresidualsensitivity`](@ref), with one factorization of the exact real
Jacobian of the operating point. Returns a matrix whose columns are `dx/dr`
for each component, in the complex representation of the augmented state.
"""
function calcnodefluxsensitivity(op::HBOperatingPoint, dFr::AbstractMatrix;
    factorization = KLUfactorization())

    Ntot = length(op.x)
    cache = FactorizationCache()
    if !isnothing(op.dc)
        # the implicit function theorem applies to the canonical system, so
        # the solve is there and the answer is scattered back: the average
        # voltages are unknowns of that system and not of this one, and the
        # node fluxes are what the caller asked for
        L = op.dc.work.layout
        tryfactorize!(cache, factorization, op.dc.jacobian)
        rhs = zeros(Float64, canonicaldim(L))
        duc = zeros(Float64, canonicaldim(L))
        dxr = zeros(Float64, L.rdim)
        dx = zeros(Complex{Float64}, Ntot, size(dFr, 2))
        for k in axes(dx, 2)
            rhs .= view(dFr, :, k)
            trysolve!(duc, cache.factorization, rhs)
            rmul!(duc, -1)
            scattercanonical!(dxr, duc, L)
            real_to_complex!(view(dx,:,k), dxr, op.modelayout.isreal)
        end
        return dx
    end
    tryfactorize!(cache, factorization, op.jacobian)
    # dFr is sparse from construction; the factorization wants a dense
    # right hand side, so densify one column at a time rather than the
    # whole state by component matrix.
    rhs = zeros(Float64, size(dFr, 1))
    dxr = zeros(Float64, size(dFr, 1))
    dx = zeros(Complex{Float64}, Ntot, size(dFr, 2))
    for k in axes(dx, 2)
        rhs .= view(dFr, :, k)
        trysolve!(dxr, cache.factorization, rhs)
        rmul!(dxr, -1)
        real_to_complex!(view(dx,:,k), dxr, op.modelayout.isreal)
    end
    return dx
end

# the Josephson contribution of the residual, restricted to junction j: the
# node vector of the Fourier coefficients of sin(phi_b(t))/Lj of that
# junction alone, which is the derivative of the residual with respect to
# minus the logarithm of that junction inductance.
function residualjosephsonterm!(out, sys, j::Integer)
    _ensuresin!(sys)
    copyto!(sys.worktd, sys.sintd)
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    for i in axes(sys.phimatrix, ndims(sys.phimatrix))
        if i != j
            selectdim(sys.phimatrix, ndims(sys.phimatrix), i) .= 0
        end
    end
    # the Josephson contribution alone, without the linear term
    applybackwardterm!(out, sys.nonlineartermplan, sys.phimatrix, sys.x;
        addlinearterm = false)
    return out
end

# the entries of M in the given row and column ranges, with every other entry
# dropped, keeping the size of M.
function maskrowscolumns(M::SparseMatrixCSC, firstrow::Integer,
    lastrow::Integer, firstcol::Integer, lastcol::Integer)
    I, J, V = findnz(M)
    keep = (firstrow .<= I .<= lastrow) .& (firstcol .<= J .<= lastcol)
    return SparseMatrixCSC{Complex{Float64},Int}(sparse(I[keep], J[keep],
        Complex{Float64}.(V[keep]), size(M,1), size(M,2)))
end




"""
    ReverseSensitivity(op, dFr, T, Em, nzrow, nzcol, realindexmap,
        branchnodes)

Everything the reverse mode contraction of [`calcSsensitivityreverse!`](@ref)
needs, precomputed once and shared read only across the signal frequencies.
"""
struct ReverseSensitivity
    op::HBOperatingPoint
    # the residual derivatives, sparse: each component touches only its own
    # nodes, so the inner product per component is over a handful of entries
    # rather than over the whole state.
    dFr::SparseMatrixCSC{Float64,Int}
    T::Matrix{Float64}
    # the transpose of the frequency to time domain map is a forward
    # transform again (the discrete Fourier matrix is symmetric), executed
    # with this plan by applyffttranspose! against per-thread work arrays.
    fftplan
    nzrow::Vector{Int}
    nzcol::Vector{Int}
    realindexmap::Vector{Int}
    branchnodes
    # the output slot of each residual derivative column: the design
    # parameter it belongs to, or its own index in the relative form
    slots::Vector{Int}
end

"""
    ReverseSensitivityBuffers(rev::ReverseSensitivity)

The mutable work arrays of one invocation of
[`calcSsensitivityreverse!`](@ref), allocated once per batch of signal
frequencies rather than at every frequency. Each thread of
[`hblinsolve`](@ref) owns its own set; the [`ReverseSensitivity`](@ref)
itself is shared read only.
"""
struct ReverseSensitivityBuffers
    P::Vector{Complex{Float64}}
    Q::Vector{Complex{Float64}}
    # the output functional covectors of a chunk of output pairs, and their
    # solutions through the transposed pump Jacobian, batched so the sparse
    # solver amortizes its per-call overhead over many right hand sides.
    G::Matrix{Complex{Float64}}
    Psi::Matrix{Complex{Float64}}
    eta::Matrix{Complex{Float64}}
    c::Matrix{Complex{Float64}}
    # the zero padded input and the single output grid of the transposed
    # transform: the holomorphic and antiholomorphic halves are folded into
    # eta one after the other, so their transforms need not coexist.
    padded::Array{Complex{Float64}}
    tgrid::Array{Complex{Float64}}
    # the covector as the junction branches produce it, in the real
    # representation of the harmonic state. With a direct current block the
    # columns of `G` are canonical and this is gathered into one of them;
    # without, it is copied straight across.
    gint::Vector{Complex{Float64}}
end

# the byte budget of the two nr x chunk complex right hand side and
# solution buffers of one frequency batch of the reverse contraction
# (together they cost 32*nr*chunk bytes, so at this budget a pump system of
# nr = 8000 real unknowns gets a chunk of about 128 columns, and a system a
# hundred times larger degrades gracefully toward single column solves).
# batching amortizes the per-call overhead of the sparse triangular solves,
# measured about 3x for tens of right hand sides, and the budget rather
# than a fixed column count keeps the per batch memory bounded for large
# pump systems; nbatches controls how many batches (and therefore budgets)
# are active at once.
const REVERSESENSITIVITYCHUNKBYTES = 32*2^20

function ReverseSensitivityBuffers(rev::ReverseSensitivity, NPM::Integer)
    sys = rev.op.sys
    NLj = size(sys.phimatrix)[end]
    NF = length(sys.phimatrix)
    nr = sensitivitydim(rev.op)
    chunk = clamp(REVERSESENSITIVITYCHUNKBYTES ÷ (32*nr), 1, NPM^2)
    return ReverseSensitivityBuffers(
        zeros(Complex{Float64}, NF), zeros(Complex{Float64}, NF),
        zeros(Complex{Float64}, nr, chunk),
        zeros(Complex{Float64}, nr, chunk),
        zeros(Complex{Float64}, size(rev.T, 1), NLj),
        zeros(Complex{Float64}, size(rev.T, 2), NLj),
        zeros(Complex{Float64}, size(sys.phitd)),
        zeros(Complex{Float64}, size(sys.phitd)),
        zeros(Complex{Float64}, size(rev.op.jacobian, 1)))
end

"""
    calcbranchtimedomainmap(sys, Nmodes, NLj)

The matrix of the map from the branch fluxes of one Josephson junction to its
physical time domain branch flux, with the real and the imaginary part of
each mode as separate columns. The map is the same for every junction,
because the packing and the inverse transform act on each junction
independently, so it is built once by transforming unit branch fluxes with
[`phivectortomatrix!`](@ref) and `applyifft!`, which keeps the conjugate mode
bookkeeping inside the functions which define it. This is the transpose of
the linear map from the unknowns to the time domain branch fluxes, restricted
to one junction.
"""
function calcbranchtimedomainmap(sys, Nmodes::Integer, NLj::Integer)
    branch = sys.Ljb.nzind[1]
    Nbranches = size(sys.Rbnm, 1) ÷ Nmodes
    T = zeros(Float64, length(sys.phitd) ÷ NLj, 2*Nmodes)
    fd = similar(sys.phimatrix)
    td = similar(sys.phitd)
    tdflat = reshape(td, :, NLj)
    bv = zeros(Complex{Float64}, Nbranches*Nmodes)
    for m in 1:Nmodes
        for (q, part) in enumerate((one(Complex{Float64}), im))
            fill!(bv, 0)
            bv[(branch-1)*Nmodes+m] = part
            fill!(fd, 0)
            phivectortomatrix!(bv[sys.Ljbm.nzind], fd, sys.freqindexmap,
                sys.conjsourceindices, sys.conjtargetindices, NLj)
            applyifft!(td, fd, sys.irfftplan)
            T[:, 2*(m-1)+q] .= view(tdflat, :, 1)
        end
    end
    return T
end

"""
    ReverseSensitivity(op::HBOperatingPoint, lsys, dFr)

Precompute the reverse mode contraction data: the branch flux map
([`calcbranchtimedomainmap`](@ref)), the transform of the pump harmonic grid,
the row and the column of each nonzero of the linearized system matrix, and
the offset of each entry of the augmented state in its real representation.
"""
function ReverseSensitivity(op::HBOperatingPoint, lsys, dFr,
        slots::Vector{Int} = collect(1:size(dFr, 2)))
    # the operating point, and therefore the branch flux map and the
    # incidence lists, live on the pump mode grid, not the signal mode grid
    Nmodes = op.Nmodes
    sys = op.sys
    # the per-frequency contraction reads the cached time domain sine of the
    # branch fluxes, so pin it to the operating point here, in this serial
    # constructor: the threads of hblinsolve only read it, and updating a
    # shared cache from them would be a race. without this the cache would
    # hold whatever point the shared evaluation object was last used at.
    setpoint!(sys, op.x)
    _ensuresin!(sys)
    NLj = size(sys.phimatrix)[end]
    A = lsys.Asparse
    nzrow = zeros(Int, nnz(A))
    nzcol = zeros(Int, nnz(A))
    rows = rowvals(A)
    for j in axes(A, 2)
        for p in nzrange(A, j)
            nzrow[p] = rows[p]
            nzcol[p] = j
        end
    end
    # the offset of complex entry j in the real representation. note that
    # isreal is indexed by mode, not by entry of the augmented state.
    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)
    realindexmap = zeros(Int, length(op.x))
    k = 1
    for j in eachindex(realindexmap)
        realindexmap[j] = k
        k += isrealmode[(j-1) % nmd + 1] ? 1 : 2
    end
    Nbranches = size(sys.Rbnm, 1) ÷ Nmodes
    return ReverseSensitivity(op, SparseMatrixCSC{Float64,Int}(dFr),
        calcbranchtimedomainmap(sys, Nmodes, NLj),
        plan_applyffttranspose(sys.phimatrix, sys.phitd), nzrow, nzcol,
        realindexmap, branchnodesandsigns(sys.Rbnm, Nmodes, Nbranches),
        slots)
end

"""
    calcSsensitivityreverse!(Ssensitivity, rev::ReverseSensitivity, lsys,
        phin, phinadjoint, gamma, beta, cache,
        bufs::ReverseSensitivityBuffers)

Add the contribution of the shift of the pump operating point to the
scattering parameter sensitivities, contracting in the reverse order so that
the cost per component is a sparse inner product rather than a product
against a matrix which is dense on the sparsity structure of the linearized
system.

For each pair of output and input port modes `(a,b)`, the transpose of the
Josephson scatter of [`addjosephsonterm!`](@ref) gives

    transpose(lam_a)*dAop_k*phi_b
        = sum_s P[s]*dcos_k[s] + Q[s]*conj(dcos_k[s]),

with `P` and `Q` accumulated over the scatter lists of the plan. Since
`dcos_k` is the directional derivative of the Fourier coefficients of
`cos(phi_b(t))` along the operating point shift
([`cosdirectionalderivative!`](@ref)), that is a linear functional of the
shift, and with

    alpha = Em*P,  gam = conj(Em)*Q,  eta = -sin(phi_b(t)).*(alpha + gam)

its covector is the transpose of the branch flux map applied to `eta`.
Finally the implicit function theorem gives `dx_k = -inv(J)*dF_k`, so pushing
that covector through the transposed Jacobian once per output pair leaves a
sparse inner product with `dF_k` for each component. The cost per signal
frequency is `(Nports*Nmodes)^2` transposed solves, independent of the number
of components, instead of one product against the full sparsity structure per
component. The solves are batched into multi right hand side calls whose
column count is set by the `REVERSESENSITIVITYCHUNKBYTES` budget, which
amortizes the per-call overhead of the sparse triangular solves while
keeping the per batch work matrices memory bounded.

The transform of the pump harmonic grid is applied one dimension at a time,
so this supports any number of pump tones.
"""
function calcSsensitivityreverse!(Ssensitivity, rev::ReverseSensitivity,
    lsys, phin, phinadjoint, gamma, beta, cache,
    bufs::ReverseSensitivityBuffers)

    op = rev.op
    # the pump mode count: the operating point shift lives on the pump grid
    Nmodes = op.Nmodes
    sys = op.sys
    plan = lsys.complexjacobianplan
    NPM = size(phin, 2)
    NLj = size(sys.phimatrix)[end]
    Ncomponents = size(rev.dFr, 2)
    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)

    # the cached time domain sine of the pump branch fluxes, pinned to the
    # operating point by the ReverseSensitivity constructor and read only
    # here.
    sintd = reshape(sys.sintd, :, NLj)
    P = bufs.P
    Q = bufs.Q
    G = bufs.G
    Psi = bufs.Psi
    eta = bufs.eta
    c = bufs.c
    padded = bufs.padded
    tgrid = bufs.tgrid
    wsize = size(sys.phimatrix)
    dFrows = rowvals(rev.dFr)
    dFvals = nonzeros(rev.dFr)

    # The output pairs are independent, so their solves through the
    # transposed pump Jacobian are batched: the covectors of a chunk of
    # pairs are accumulated as the columns of G and pushed through the
    # factorization in one multi right hand side call, which amortizes the
    # per-call overhead of the sparse triangular solves over the chunk.
    wcov = zeros(eltype(P), plan.josephson.n)
    pairs = vec(CartesianIndices((NPM, NPM)))
    for chunk in Iterators.partition(eachindex(pairs), size(G, 2))
        for (col, pi) in enumerate(chunk)
            a, b = Tuple(pairs[pi])
            # the transpose of the Josephson scatter
            fill!(P, 0)
            fill!(Q, 0)
            # the covector per destination, then the transpose of the
            # Josephson map applied to it. The plain and conjugated halves
            # separate by the sign of the mode coupling index.
            @inbounds for p in eachindex(wcov)
                wcov[p] = phinadjoint[rev.nzrow[p],a]*phin[rev.nzcol[p],b]
            end
            josephsonadjoint!(P, Q, plan.josephson, wcov)

            # the transpose of the cos directional derivative: a forward
            # transform of the zero padded coefficients through the plan of
            # the operating point grid (the transform matrix is symmetric),
            # and the antiholomorphic half by conjugating around the same
            # transform. the two halves are folded into eta one after the
            # other through the same output grid.
            applyffttranspose!(tgrid, reshape(P, wsize), padded, rev.fftplan)
            tf = reshape(tgrid, :, NLj)
            @inbounds for i in eachindex(eta)
                eta[i] = -sintd[i]*tf[i]
            end
            @inbounds for i in eachindex(Q)
                Q[i] = conj(Q[i])
            end
            applyffttranspose!(tgrid, reshape(Q, wsize), padded, rev.fftplan)
            @inbounds for i in eachindex(eta)
                eta[i] -= sintd[i]*conj(tf[i])
            end
            mul!(c, transpose(rev.T), eta)

            # the transpose of the branch flux map, into the real
            # representation of the augmented state. The outputs depend on
            # the state through the harmonic coordinates, so this is where
            # the covector is built whether or not there is a direct current
            # block; with one it is carried into the canonical coordinates
            # afterwards, which is the transpose of the scatter the residual
            # applies going the other way.
            g = bufs.gint
            fill!(g, 0)
            @inbounds for (jj, branch) in enumerate(sys.Ljb.nzind)
                for m in 1:Nmodes
                    cre = c[2*(m-1)+1, jj]
                    cim = c[2*(m-1)+2, jj]
                    holo = (cre - im*cim)/2
                    anti = (cre + im*cim)/2
                    for (node, sgn) in rev.branchnodes[branch]
                        j = (node-1)*Nmodes + m
                        k = rev.realindexmap[j]
                        g[k] += sgn*(holo + anti)
                        if !isrealmode[(j-1) % nmd + 1]
                            g[k+1] += sgn*im*(holo - anti)
                        end
                    end
                end
            end
            gcol = view(G, :, col)
            if isnothing(op.dc)
                copyto!(gcol, g)
            else
                # the gather leaves the voltage rows alone, and they are
                # zero: no output reads an average voltage directly, only
                # through the state the block moves
                fill!(gcol, 0)
                gathercanonical!(gcol, g, op.dc.work.layout)
            end
        end

        # through the transposed Jacobian once for the whole chunk
        ncols = length(chunk)
        Gc = view(G, :, 1:ncols)
        Psic = view(Psi, :, 1:ncols)
        trysolvetranspose!(Psic, cache.factorization, Gc)

        # the sparse inner products with the residual derivatives
        for (col, pi) in enumerate(chunk)
            a, b = Tuple(pairs[pi])
            @inbounds for k in 1:Ncomponents
                acc = zero(Complex{Float64})
                for r in nzrange(rev.dFr, k)
                    acc += Psi[dFrows[r], col]*dFvals[r]
                end
                Ssensitivity[a,b,rev.slots[k]] += gamma[a]*beta[b]*acc
            end
        end
    end
    return nothing
end

"""
    calcoperatingpointstamps(op::HBOperatingPoint, lsys, dx)

Calculate the contribution of a shift of the pump operating point to the
derivative of the linearized system matrix, for each column of `dx`. The
linearized system matrix depends on the operating point only through the
Fourier coefficients of `cos(phi_b(t))` of the Josephson junction branch
fluxes, so the contribution is the directional derivative of those
coefficients along the operating point shift
([`cosdirectionalderivative!`](@ref)) scattered into the system matrix
through the same plan which assembles it ([`addjosephsonterm!`](@ref)), so
the mode coupling and its truncation agree exactly. Returns a vector of
nonzero value vectors aligned with the sparsity structure of the system
matrix, which are frequency independent.
"""
function calcoperatingpointstamps(op::HBOperatingPoint, lsys, dx)
    setpoint!(op.sys, op.x)
    dcos = similar(op.sys.phimatrix)
    stamps = Vector{Vector{Complex{Float64}}}(undef, size(dx,2))
    for k in axes(dx, 2)
        cosdirectionalderivative!(dcos, op.sys,
            Vector{Complex{Float64}}(view(dx,:,k)))
        nzval = zeros(Complex{Float64}, nnz(lsys.Asparse))
        addjosephsonterm!(nzval, lsys.complexjacobianplan, dcos)
        stamps[k] = nzval
    end
    return stamps
end

"""
    SensitivityStamp(kind, M, indexmap, freqsubstindices, nzval, portindex)

The derivative of the linearized harmonic balance system matrix with respect
to a relative (logarithmic) perturbation of one component value, `p -> r*p`
evaluated at `r = 1`. The system matrix is affine in `C`, `1/R`, `1/L` and
`1/Lj`, so the derivative is that component's own contribution to the system
matrix, with a sign: positive for the capacitance and negative for the
quantities which enter inversely.

`kind` selects how the stamp is assembled at each signal frequency, mirroring
[`assemblesystemmatrix!`](@ref):

- `:C`: `-M*wmodes2m`, the component's capacitance matrix,
- `:G`: `-im*M*wmodesm`, the component's conductance matrix, either the nodal
    stamp or, for a resistor promoted by the modified nodal analysis
    formulation, the conductances of its constitutive equations,
- `:invL`: `-M`, the component's inverse inductance matrix,
- `:Lj`: the constant nonzero values `nzval`, the negative of the pump
    modulated Josephson contribution of that junction alone.

`indexmap` maps the nonzeros of `M` into the nonzeros of the system matrix
and `freqsubstindices` locates any symbolic entries. `portindex` is the
index of the port whose impedance this component is, or zero, which selects
the additional wave normalization term of [`calcSsensitivity!`](@ref).

`parameter` is the design parameter this stamp contributes to, or zero for
the legacy relative form where every component is its own output slot.
`portscale` is `(dZport/dtheta)/Zport` for a port impedance stamp, one for
the relative form.
"""
struct SensitivityStamp
    kind::Symbol
    rows::Vector{Int}
    cols::Vector{Int}
    vals::Vector{Complex{Float64}}
    portindex::Int
    parameter::Int
    portscale::Complex{Float64}
end

# the relative form is the special case dc/dtheta = c, one output slot per
# component
SensitivityStamp(kind, rows, cols, vals, portindex) =
    SensitivityStamp(kind, rows, cols, vals, portindex, 0,
        one(Complex{Float64}))

"""
    reparameterize(stamp::SensitivityStamp, alpha, parameter)

The same stamp expressed as the derivative with respect to a real design
parameter rather than a relative perturbation of the component value.

A component's contribution to the system matrix is linear in its value `c`,
and the negative frequency conjugation is applied to the *stored* value at
contraction time ([`sensitivitystampvalue`](@ref)), so

    d/dtheta modevalue(c, w) = modevalue(dc/dtheta, w),

and the stamp for `theta` is the stamp for `c` rescaled by
`alpha = (dc/dtheta)/c`. Nothing else changes: the sparsity, the kind and
the frequency scaling are all properties of the component, not of the
parameterization.

This is what makes the real parameter form correct for complex component
values, where the relative form is not. A relative perturbation moves `c`
along itself, so `dS/dr` is a single complex number which cannot resolve
the two real directions of a complex `c`; carrying `dc/dtheta` as the
direction resolves them, because a real `theta` which rotates `c` in the
complex plane produces a `dc/dtheta` which is not parallel to `c`.
"""
function reparameterize(stamp::SensitivityStamp, alpha::Number,
        parameter::Integer)
    return SensitivityStamp(stamp.kind, stamp.rows, stamp.cols,
        stamp.vals .* alpha, stamp.portindex, parameter,
        Complex{Float64}(alpha))
end

# =====================================================================
# Scattering block design parameter sensitivities.
#
# A block's contribution to the system matrix is the hybrid stamp of
# B = R^(-1/2)(I - S) and C = R^(1/2)(I + S), affine in S, and the
# assembly is linear in B and C, so the derivative of the stamped values
# with respect to a parameter is the assembly run with S replaced by
# dS/dtheta minus the assembly with S replaced by zero: the subtraction
# removes the identity parts of B and C and cancels the zero frequency
# rows, which do not depend on S. The stamp values depend on the signal
# frequency through S itself, so unlike the lumped stamps they are
# rebuilt at every frequency -- by each worker into its own buffers,
# because the workers of hblinsolve run concurrently.
# =====================================================================

"""
    zeroscatteringblock(b::ScatteringParameters)

A block with the same ports, reference impedances and conventions whose
scattering matrix is identically zero. Constructed directly rather than
through the public constructor: it needs no passivity check.
"""
zeroscatteringblock(b) = ScatteringParameters(
    ConstantMatrixProvider(zeros(Complex{Float64}, b.nports, b.nports)),
    b.nports, b.zref, b.grounded, b.noise, b.negative_frequency)

"""
    derivativestampsystems(ssys, targetdef, dblock)

The pair of provider-swapped stamp systems whose value difference is the
derivative of the block contribution with respect to a parameter of the
block `targetdef`: the first evaluates `dblock` (the block whose S is
dS/dtheta) in place of the target and zero for every other block, the
second zero for all. Everything else -- the sparsity, the contribution
lists, the auxiliary layout -- is shared with `ssys` by construction.
"""
function derivativestampsystems(ssys, targetdef, dblock)
    swap(pick) = ScatteringStampSystem(
        [StampedScatteringBlock(pick(sb.block), sb.signalnodes, sb.refnodes,
             sb.auxbase, sb.name) for sb in ssys.blocks],
        ssys.kcl, ssys.pattern, ssys.patternindex, ssys.Aindex,
        ssys.blockindex, ssys.pindex, ssys.qindex, ssys.coeff, ssys.sign,
        ssys.modeindex, ssys.Nmodes, ssys.Nauxports, ssys.scale)
    dsys = swap(b -> b === targetdef ? dblock : zeroscatteringblock(b))
    zsys = swap(zeroscatteringblock)
    return dsys, zsys
end

"""
    blockstampvals!(vals, dsys, zsys, wmodes, workd, workz, dbuf, zbuf,
        patternindex)

Overwrite `vals` (in the nonzero order of the stamp pattern) with the
derivative of the block contribution at the signed mode frequencies
`wmodes`: the values of the dS system minus the values of the zero system,
scattered through `patternindex`.
"""
function blockstampvals!(vals, dsys, zsys, wmodes, workd, workz, dbuf,
        zbuf, patternindex)
    scatteringvalues!(dbuf, dsys, wmodes, workd)
    scatteringvalues!(zbuf, zsys, wmodes, workz)
    fill!(vals, 0)
    @inbounds for c in eachindex(patternindex)
        vals[patternindex[c]] += dbuf[c] - zbuf[c]
    end
    return vals
end

"""
    WorkerBlockSensitivity

One worker's private state for the scattering block sensitivity stamps:
its own copy of the stamp vector (the lumped stamps are shared read only,
the `:S` stamps' values are private because they are rebuilt at every
signal frequency), the provider-swapped systems, and the evaluation
scratch.
"""
struct WorkerBlockSensitivity
    stamps::Vector{SensitivityStamp}
    # (index into stamps, dsys, zsys), one per block parameter pair
    entries::Vector{Tuple{Int,Any,Any}}
    patternindex::Vector{Int}
    workd::ScatteringWorkspace
    workz::ScatteringWorkspace
    dbuf::Vector{Complex{Float64}}
    zbuf::Vector{Complex{Float64}}
end

function WorkerBlockSensitivity(stamps::Vector{SensitivityStamp},
        blockentries, patternindex)
    private = SensitivityStamp[st.kind == :S ?
        SensitivityStamp(st.kind, st.rows, st.cols, copy(st.vals),
            st.portindex, st.parameter, st.portscale) : st
        for st in stamps]
    m = length(patternindex)
    return WorkerBlockSensitivity(private, collect(blockentries),
        collect(Int, patternindex), ScatteringWorkspace(),
        ScatteringWorkspace(), Vector{Complex{Float64}}(undef, m),
        Vector{Complex{Float64}}(undef, m))
end

"""
    refreshblockstamps!(wbs::WorkerBlockSensitivity, wmodes)

Rebuild this worker's `:S` stamp values at the signed mode frequencies
`wmodes`.
"""
function refreshblockstamps!(wbs::WorkerBlockSensitivity, wmodes)
    for (idx, dsys, zsys) in wbs.entries
        blockstampvals!(wbs.stamps[idx].vals, dsys, zsys, wmodes,
            wbs.workd, wbs.workz, wbs.dbuf, wbs.zbuf, wbs.patternindex)
    end
    return wbs
end

"""
    blocksensitivitystamp(ssys, parameter)

The [`SensitivityStamp`](@ref) skeleton of a scattering block parameter:
the triplet positions of the frequency dependent block pattern with zero
values, which each worker fills at each signal frequency
([`refreshblockstamps!`](@ref)). Kind `:S` carries no further frequency
scaling, exactly as `:Lj` does.
"""
function blocksensitivitystamp(ssys, parameter::Integer)
    rows, cols, _ = findnz(ssys.pattern)
    return SensitivityStamp(:S, rows, cols,
        zeros(Complex{Float64}, length(rows)), 0, Int(parameter),
        one(Complex{Float64}))
end

"""
    calcblockresidualsensitivity(op::HBOperatingPoint, psc, blockpairs)

The derivative of the harmonic balance residual with respect to each
scattering block design parameter of `blockpairs = [(firstportname,
parameterindex, derivativeblock)]`, at the operating point, in the real
representation of the augmented residual (one column per pair, matching
the columns [`calcresidualsensitivity`](@ref) produces for lumped pairs).

A block enters the nonlinear system as a constant matrix folded into the
inverse inductance augmentation, so it moves the pump operating point and
its residual derivative is `+dA_block/dtheta` applied to the converged
state: the linear term enters the residual with a plus sign, and unlike
the lumped `:invL` case -- whose minus encodes `d(1/L)/d(ln L) = -1/L`,
the derivative of the component law, not a residual sign -- the block
derivative is computed directly. The derivative
matrix comes from the same affine subtraction as the linearized stamps,
on a stamp system rebuilt for the pump mode grid with the geometry the
operating point already determines: the augmented dimension is the state
length and the block auxiliary variables occupy its tail.
"""
function calcblockresidualsensitivity(op::HBOperatingPoint,
        psc::CompiledCircuit, blockpairs::AbstractVector)
    Ntot = length(op.x)
    Nmodes = op.Nmodes
    Nsc = countscatteringports(psc)*Nmodes
    pumpssys = scatteringstampsystem(psc.scatteringblocks, Nmodes;
        auxoffset = Ntot - Nsc, Ntotal = Ntot, scale = real(op.Lmean))
    isnothing(pumpssys) && throw(ArgumentError(
        "the circuit has no scattering blocks to take a block sensitivity of"))

    # the offset of complex entry r in the real representation, matching
    # complex_to_real! with the default scale (see calcresidualsensitivity)
    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)
    realindexmap = zeros(Int, Ntot)
    k = 1
    for j in eachindex(realindexmap)
        realindexmap[j] = k
        k += isrealmode[(j-1) % nmd + 1] ? 1 : 2
    end
    Nreal = realdim(Ntot, isrealmode)

    dFr = zeros(Float64, Nreal, length(blockpairs))
    work = ScatteringWorkspace()
    workz = ScatteringWorkspace()
    m = length(pumpssys.patternindex)
    dbuf = Vector{Complex{Float64}}(undef, m)
    zbuf = Vector{Complex{Float64}}(undef, m)
    dA = copy(pumpssys.pattern)
    for (col, bp) in enumerate(blockpairs)
        bi = scatteringblockindex(psc, bp[1])
        iszero(bi) && throw(ArgumentError(lazy"the block pair names $(bp[1]), which is not a scattering block of this circuit"))
        targetdef = psc.scatteringblocks[bi].definition
        dsys, zsys = derivativestampsystems(pumpssys, targetdef, bp[3])
        blockstampvals!(nonzeros(dA), dsys, zsys, op.wmodes, work, workz,
            dbuf, zbuf, pumpssys.patternindex)
        dAx = dA*op.x
        for r in eachindex(dAx)
            iszero(dAx[r]) && continue
            v = dAx[r]
            kr = realindexmap[r]
            dFr[kr, col] += real(v)
            if !isrealmode[(r-1) % nmd + 1]
                dFr[kr+1, col] += imag(v)
            end
        end
    end
    return SparseMatrixCSC{Float64,Int}(sparse(dFr))
end

"""
    parametergrouping(pairs, componentindices, componenttypes, portordinal)

The grouping of the sensitivity pairs `(componentname, parameterindex,
alpha)` into merged stamps: pairs which share a stamp kind, a design
parameter and a port ordinal merge into one contraction.
`componentindices[i]` is the parsed circuit index of pair `i`. Returns
`(grouping, slots)`, the pair indices of each group and the design
parameter each group accumulates into.

Computed from the parsed circuit alone, before any stamp exists, because
the residual derivative columns and the operating point solves must be
merged with the same grouping as the stamps, and both are needed earlier
than the stamps are built.
"""
function parametergrouping(pairs, componentindices, componenttypes,
        portordinal)
    kindof(t) = t == :C ? :C : t == :R ? :G : t == :L ? :invL :
        t == :Lj ? :Lj : t
    groups = Dict{Tuple{Symbol,Int,Int},Int}()
    grouping = Vector{Int}[]
    slots = Int[]
    for i in eachindex(pairs)
        ci = componentindices[i]
        pj = pairs[i][2]
        key = (kindof(componenttypes[ci]), pj, get(portordinal, ci, 0))
        g = get!(groups, key) do
            push!(grouping, Int[])
            push!(slots, pj)
            length(grouping)
        end
        push!(grouping[g], i)
    end
    return grouping, slots
end

"""
    mergestamps(stamps, grouping)

Concatenate the stamps of each group of `grouping` into one. A design
parameter typically touches many components -- a single junction
inductance across a two thousand cell line -- and the contraction cost is
per stamp, so merging turns one contraction per component into one per
parameter (and kind). The grouping comes from
[`parametergrouping`](@ref), so it is the same one applied to the residual
derivative columns.
"""
function mergestamps(stamps::AbstractVector{SensitivityStamp}, grouping)
    out = SensitivityStamp[]
    for idx in grouping
        if length(idx) == 1
            push!(out, stamps[idx[1]])
        else
            push!(out, SensitivityStamp(stamps[idx[1]].kind,
                vcat((stamps[i].rows for i in idx)...),
                vcat((stamps[i].cols for i in idx)...),
                vcat((stamps[i].vals for i in idx)...),
                stamps[idx[1]].portindex, stamps[idx[1]].parameter,
                stamps[idx[1]].portscale))
        end
    end
    return out
end
"""
    calcsensitivitystamps(sensitivityindices, psc, cg, nm, lsys, phimatrix,
        mnaindices, coupledbranches, Nnodalmna, Nmodes, Nnodes)

Build the [`SensitivityStamp`](@ref) of each component in
`sensitivityindices`. The classification and the raw one-component matrices
come from [`componentstamp`](@ref), which is shared with the residual
derivatives of [`calcresidualsensitivity`](@ref), so the two grids cannot
disagree on which components are supported or how they are built. The
Josephson junction stamp is the pump modulated contribution of that
junction alone, obtained by scattering the Fourier coefficients of
`cos(phi(t))` of that junction through the same plan
([`addjosephsonterm!`](@ref)) which assembles the system matrix, so the mode
coupling and its truncation agree exactly.
"""
function calcsensitivitystamps(sensitivityindices, psc::CompiledCircuit,
    cg::CircuitGraph, nm::CircuitMatrices, lsys,
    phimatrix, mnaindices, coupledbranches, Nnodalmna, Nmodes, Nnodes)

    Ntot = size(lsys.Asparse, 1)
    stamps = Vector{SensitivityStamp}(undef, length(sensitivityindices))
    lookups = componentlookups(mnaindices, coupledbranches, nm.Ljb)
    # a sensitivity taken with respect to a port's own environment also
    # moves the wave normalization, so those components are recognized by
    # their role; a port which owns no environment contributes none
    portordinal = Dict(idx => p
        for (p, idx) in enumerate(nm.portenvironmentindices) if !iszero(idx))

    for (k, idx) in enumerate(sensitivityindices)
        portindex = get(portordinal, idx, 0)
        kind, info = componentstamp(idx, psc, cg, nm, lookups,
            Nmodes, Nnodes)
        if kind == :C
            stamps[k] = tripletstamp(:C, mnapadto(info, Ntot), portindex)
        elseif kind == :G
            stamps[k] = tripletstamp(:G, mnapadto(info, Ntot), portindex)
        elseif kind == :invL
            stamps[k] = tripletstamp(:invL, mnapadto(info, Ntot), portindex)
        else # :Lj
            # the contribution of this junction alone, by zeroing the Fourier
            # coefficients of every other junction before the scatter.
            onejunction = zero(phimatrix)
            selectdim(onejunction, ndims(onejunction), info) .=
                selectdim(phimatrix, ndims(phimatrix), info)
            nzval = zeros(Complex{Float64}, nnz(lsys.Asparse))
            addjosephsonterm!(nzval, lsys.complexjacobianplan, onejunction)
            rmul!(nzval, -1)
            stamps[k] = tripletstamp(:Lj,
                SparseMatrixCSC(size(lsys.Asparse,1), size(lsys.Asparse,2),
                    copy(lsys.Asparse.colptr), copy(lsys.Asparse.rowval),
                    nzval), portindex)
        end
    end
    return stamps
end

# pad a matrix with empty rows and columns for the auxiliary variables of the
# modified nodal analysis formulation, and convert to the element type of the
# system matrix.
function mnapadto(M::SparseMatrixCSC, Ntot::Integer)
    padded = size(M,1) == Ntot ? M : mnapad(M, Ntot - size(M,1))
    return SparseMatrixCSC{Complex{Float64},Int}(padded)
end

# convert a component matrix to the compact triplet form of a
# SensitivityStamp, dropping structural zeros. The stamps of the individual
# components are extremely sparse compared with the system matrix (a
# capacitor to ground touches one node), so the contraction is driven by
# these entries rather than by the sparsity structure of the system matrix.
function tripletstamp(kind::Symbol, M::SparseMatrixCSC, portindex::Integer)
    I, J, V = findnz(M)
    keep = .!iszero.(V)
    return SensitivityStamp(kind, I[keep], J[keep],
        Complex{Float64}.(V[keep]), portindex)
end

"""
    sensitivitystampvalue(stamp::SensitivityStamp, t::Integer, wmodes,
        Nmodes)

The value of entry `t` of the derivative of the linearized harmonic balance
system matrix with respect to a relative perturbation of the component of
`stamp`, at the mode frequencies `wmodes`. Applies the same per mode
frequency scaling and negative frequency mode conjugation as
[`assemblesystemmatrix!`](@ref), which are indexed by the column.
"""
@inline function sensitivitystampvalue(stamp::SensitivityStamp, t::Integer,
    wmodes, Nmodes)

    (stamp.kind == :Lj || stamp.kind == :S) && return stamp.vals[t]
    m = (stamp.cols[t] - 1) % Nmodes + 1
    w = wmodes[m]
    v = modevalue(stamp.vals[t], w)
    if stamp.kind == :C
        return -v*w^2
    elseif stamp.kind == :G
        return -im*v*w
    else
        return -v
    end
end

"""
    calcsensitivityscaling!(gamma, beta, inputwave, bnm,
        portindices, portimpedances, componenttypes, nodeindices,
        wmodes, Nmodes)

Calculate the scalars which convert the derivative of the node fluxes into
the derivative of the scattering parameters. Writing the output wave of
[`calcinputoutput!`](@ref) as a linear functional of the node fluxes, the
part which depends on them is
`(1/2)*kval*(1 + conj(Z)/Z)*im*w_n*(phi_n1 - phi_n2)`, and the scattering
parameters divide by the input wave, so

    dS[(j,n),(i,m)] = gamma[(j,n)]*beta[(i,m)]
                      *(dphi[node1,(j,n)] - dphi[node2,(j,n)])

with `gamma = (1/2)*kval*(1 + conj(Z)/Z)*im*w_n/s_{(j,n)}` and
`beta = 1/inputwave[(i,m),(i,m)]`. The node flux difference is contracted
with the adjoint solution by [`calcSsensitivity!`](@ref), which is why the
`im*w_n` factor of the port voltage is folded into `gamma` here: it is
exactly the factor relating the adjoint source vector to the source vector
of the forward problem.

`s_{(j,n)}` is the source current of the port's own unit drive
([`calcsourcecurrent`](@ref) on the diagonal), `±1` depending on whether the
canonical orientation of the port branch in the incidence matrix agrees with
the node order of the port component. The adjoint solution is the solve
against the source columns of `bnm`, which carry the canonical branch
orientation, while the output functional differences the node fluxes in the
component node order, so their ratio enters the contraction. Without it the
sensitivities of any output at a port written with its nodes in the opposite
order of the branch orientation would have the wrong sign, even though the
scattering parameters themselves, which use `s` consistently in both the
input and the output waves, would be correct.
"""
function calcsensitivityscaling!(gamma, beta, inputwave, bnm,
    portindices, portimpedances, componenttypes, nodeindices,
    wmodes, Nmodes)

    for i in eachindex(portindices)
        for j in 1:Nmodes
            row = (i-1)*Nmodes + j
            portimpedance = calcimpedance(portimpedances[i],
                componenttypes[portindices[i]], wmodes[j], nothing)
            kval = portwavescale(portimpedance, wmodes[j])
            # the orientation of the port branch relative to the node order
            # of the port component, from the source current of the port's
            # own unit drive.
            sourcecurrent = calcsourcecurrent(
                nodeindices[1,portindices[i]],
                nodeindices[2,portindices[i]], bnm, Nmodes, j, row)
            gamma[row] = iszero(sourcecurrent) ? 0 :
                1/2*kval*(1 + conj(portimpedance)/portimpedance)*
                im*wmodes[j]/sourcecurrent
            beta[row] = iszero(inputwave[row,row]) ? 0 :
                1/inputwave[row,row]
        end
    end
    return nothing
end

"""
    calcSsensitivity!(Ssensitivity, stamps, dA, dAphin, phin, phinadjoint,
        S, gamma, beta, wmodes, Nmodes, symfreqvar)

Calculate the derivative of the scattering matrix with respect to a relative
(logarithmic) perturbation of each component value, `p -> r*p` evaluated at
`r = 1`, at the pump operating point, with the adjoint method. Overwrites
`Ssensitivity`.

Differentiating the linearized system `A*phi = b`, whose source terms do not
depend on any component value, gives `dphi = -inv(A)*dA*phi`, so with the
adjoint solutions `lam` of the transposed system driven by the output
functionals of [`calcsensitivityscaling!`](@ref),

    dS[(j,n),(i,m)] = -gamma[(j,n)]*beta[(i,m)]
                      *transpose(lam[:,(j,n)])*dA*phi[:,(i,m)].

The adjoint source vectors are the source vectors of the forward problem
scaled by `im*w_n`, which is folded into `gamma`, so `phinadjoint`, the
solution of the transposed system already computed for the noise and quantum
efficiency calculations, is used directly.

When the perturbed component is itself a port impedance the wave
normalization of [`calcinputoutput!`](@ref) moves as well, contributing the
additional closed form term `-(1/2)*(P*(S+I) + (S+I)*P)` with `P` the
projector onto that port, which is exact for constant real port impedances.
"""
function calcSsensitivity!(Ssensitivity, stamps, dAop, dA, dAphin, phin,
    phinadjoint, S, gamma, beta, contraction, wmodes, Nmodes, symfreqvar)

    NPM = size(phin, 2)
    fill!(Ssensitivity, 0)
    for (k, stamp) in enumerate(stamps)
        # a stamp declares which output slot it accumulates into: the design
        # parameter it belongs to, or its own index in the legacy relative
        # form where every component is its own parameter
        slot = stamp.parameter == 0 ? k : stamp.parameter
        # the component's own contribution, driven by the entries of its
        # stamp rather than by the sparsity structure of the system matrix,
        # of which the stamp of a single component is a tiny part.
        fill!(contraction, 0)
        @inbounds for t in eachindex(stamp.rows)
            i = stamp.rows[t]
            j = stamp.cols[t]
            v = sensitivitystampvalue(stamp, t, wmodes, Nmodes)
            iszero(v) && continue
            for b in 1:NPM
                vphi = v*phin[j,b]
                iszero(vphi) && continue
                for a in 1:NPM
                    contraction[a,b] += phinadjoint[i,a]*vphi
                end
            end
        end

        # the contribution of the shift of the pump operating point, which is
        # frequency independent but dense on the sparsity structure of the
        # system matrix, so it goes through a sparse matrix vector product.
        if !isempty(dAop)
            copyto!(nonzeros(dA), dAop[k])
            mul!(dAphin, dA, phin)
            mul!(contraction, transpose(phinadjoint), dAphin, 1, 1)
        end

        for b in 1:NPM
            for a in 1:NPM
                Ssensitivity[a,b,slot] += -gamma[a]*beta[b]*contraction[a,b]
            end
        end

        # the wave normalization term when the component is a port impedance
        if stamp.portindex > 0
            p = stamp.portindex
            for b in 1:NPM
                for a in 1:NPM
                    sI = S[a,b] + (a == b ? 1 : 0)
                    correction = zero(Complex{Float64})
                    if (a-1) ÷ Nmodes + 1 == p
                        correction += sI/2
                    end
                    if (b-1) ÷ Nmodes + 1 == p
                        correction += sI/2
                    end
                    Ssensitivity[a,b,slot] -= correction*stamp.portscale
                end
            end
        end
    end
    return nothing
end


