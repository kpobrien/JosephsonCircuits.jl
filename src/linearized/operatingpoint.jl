# Sensitivities of the scattering parameters to component values, by the
# adjoint method, at a fixed pump operating point or including the shift of
# the operating point through the implicit function theorem.

"""
    HBOperatingPoint(sys, x, jacobian, modelayout, Nnodal, Lscale, wmodes,
        Amna, coupledbranches, Nmodes, Nnodes[, dc])

The converged pump operating point of [`hbnlsolve`](@ref) together with
everything needed to propagate a component perturbation through it: the
[`HBSystem`](@ref) evaluation object, the converged augmented state, the
exact Jacobian of the equivalent real system assembled there, and the
scaled matrices and layout of the augmented system.

Requested with `returnoperatingpoint = true`. The Jacobian is the exact
Jacobian of the equivalent real system, assembled with
[`assemblerealjacobian!`](@ref), rather than the complex holomorphic
Jacobian of the [`QuasiNewton`](@ref) method, which is only an approximation: the
harmonic balance residual is not complex differentiable, so the implicit
function theorem does not hold with the holomorphic Jacobian, while in the
real representation it applies directly.
"""
struct HBOperatingPoint
    # the HBSystem evaluation object and the ModeLayout, untyped: the
    # system's type depends on the backend and the layout's on its index type
    sys
    x::Vector{Complex{Float64}}
    jacobian::SparseMatrixCSC{Float64,Int}
    modelayout
    Nnodal::Int
    Lscale::Complex{Float64}
    wmodes::Vector{Float64}
    Amna::SparseMatrixCSC
    coupledbranches::Vector{Int}
    Nmodes::Int
    Nnodes::Int
    # The explicit direct current block, when the circuit has one: the
    # canonical work, the canonical state the solve converged to, and the
    # canonical Jacobian there. `jacobian` above is the harmonic one, read
    # by everything which differentiates the harmonic system alone; what
    # differentiates the whole system reads these.
    dc
end

# an operating point of a circuit with no direct current block
HBOperatingPoint(sys, x, jacobian, modelayout, Nnodal, Lscale, wmodes, Amna,
    coupledbranches, Nmodes, Nnodes) =
    HBOperatingPoint(sys, x, jacobian, modelayout, Nnodal, Lscale, wmodes,
        Amna, coupledbranches, Nmodes, Nnodes, nothing)

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

The perturbation is relative by default: `G0` is proportional to `1/R`,
so a relative change in `R` scales the whole stamp by `-1`; `alphas`
carries an absolute direction `dv/dp` instead, of which only the real part
reaches the direct current conductance rows.
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
    # the solver scale `Z0/w0` (see `calcsolverscale`), which every row of
    # the system is multiplied by; the block's own rows and their derivative
    # must be scaled by the same one
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
        n1 > 1 && (push!(I, L.dcpos[n1 - 1]); push!(J, k); push!(V, cur))
        n2 > 1 && (push!(I, L.dcpos[n2 - 1]); push!(J, k); push!(V, -cur))
        # and into the transport rows, which are the component sums of
        # those; a resistor inside one component cancels there, exactly as
        # an inductor branch does
        c1 = plan.componentof[n1]
        c2 = plan.componentof[n2]
        c1 == c2 && continue
        iszero(c1) || (push!(I, L.rdim + c1); push!(J, k); push!(V, cur))
        iszero(c2) || (push!(I, L.rdim + c2); push!(J, k); push!(V, -cur))
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
    componentlookups(coupledbranches, Ljb)

Constant time lookups for [`componentstamp`](@ref), built once per stamp
table rather than searched per component: the set of mutually coupled
branches and the ordinal of a junction branch within `Ljb.nzind`. Without
these the classification repeats linear searches per component, which
becomes quadratic over a large sensitivity set.
"""
function componentlookups(coupledbranches, Ljb)
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
- `(:G, M)`: the component's conductance matrix,
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
        sensitivityindices, alphas = ones(Complex{Float64}, length(sensitivityindices)))

Calculate the derivative of the harmonic balance residual with respect to a
relative (logarithmic) perturbation of each component value, at the
operating point. `alphas` scales the direction of each column: the derivative
of column `k` is with respect to `alphas[k]*r` applied to component
`sensitivityindices[k]`, which is how a design parameter's direction
`alpha = (dv/dp)/v` is folded into the residual derivative of a component
(see [`designsensitivities`](@ref)). Combined with the implicit function theorem applied to
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
        rmul!(Ms, op.Lscale)
        return mnapadto(Ms, Ntot)
    end

    # The residual derivatives are sparse: each component touches only its
    # own rows of the state (its nodes and modes, or one junction branch's
    # Kirchhoff rows), so they are accumulated as triplets of the real
    # representation and returned as a sparse matrix, where a dense array
    # would cost O(Nstate*Ncomponents), the many component regime the
    # reverse contraction order exists for. Duplicate triplets (a component
    # matrix with several entries in one row) are summed by `sparse`, which
    # is right because the real representation is linear.
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
    lookups = componentlookups(op.coupledbranches, op.sys.Ljb)
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
    # component values. The gather copies the flux rows and writes nothing
    # into the voltage rows, which is where the block's part lands.
    L = op.dc.work.layout
    N = canonicaldim(L)
    dccols = dcresidualsensitivity(op.dc, psc, nm, op.Lscale,
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
    keep[windowindices(L)] .= Array(up.keep)
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
[`calcresidualsensitivity`](@ref), with one factorization of
[`sensitivityjacobian`](@ref) at the operating point: the canonical
Jacobian when a direct current block is active, the harmonic one
otherwise. Returns a matrix whose columns are `dx/dr`
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
    applyfft!(sys.phimatrix, sys.sintd, sys.rfftplan)
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

"""
    ReverseSensitivity(op, dFr, T, fftplan, nzrow, nzcol, realindexmap,
        branchnodes, slots)

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
    ReverseSensitivityBuffers(rev::ReverseSensitivity, NPM::Integer)

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

# The byte budget of the two `nr` by `chunk` complex right hand side and
# solution buffers of one frequency batch of the reverse contraction.
# Together they cost `32*nr*chunk` bytes, so a pump system of 8000 real
# unknowns gets a chunk of about 128 columns and a system a hundred times
# larger degrades toward single column solves. Batching amortizes the per
# call overhead of the sparse triangular solves, and a byte budget rather
# than a fixed column count keeps the memory per batch bounded; one set of
# these buffers is active per `hblinsolve` batch.
"""
    REVERSESENSITIVITYCHUNKBYTES

The byte budget of one chunk of the reverse contraction's dense buffers,
which sets the column count [`calcSsensitivityreverse!`](@ref) works in.
"""
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
    ReverseSensitivity(op::HBOperatingPoint, lsys, dFr,
        slots = collect(1:size(dFr, 2)))

Precompute the reverse mode contraction data: the branch flux map
([`calcbranchtimedomainmap`](@ref)), the transform of the pump harmonic grid,
the row and the column of each nonzero of the linearized system matrix, and
the offset of each entry of the augmented state in its real representation.
`dFr` holds the residual derivative columns and `slots[k]` is the output
slot of `Ssensitivity` column `k` accumulates into, so that several columns
belonging to one design parameter can share a slot.
"""
function ReverseSensitivity(op::HBOperatingPoint, lsys, dFr,
        slots::Vector{Int} = collect(1:size(dFr, 2)))
    # the operating point, and therefore the branch flux map and the
    # incidence lists, live on the pump mode grid, not the signal mode grid
    Nmodes = op.Nmodes
    sys = op.sys
    # the per frequency contraction reads the cached time domain sine of the
    # branch fluxes, so it is pinned to the operating point here, in this
    # serial constructor: the threads of hblinsolve only read the shared
    # evaluation object, and updating it from them would be a race
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

    alpha = T(P),  gam = conj(T(Q)),  eta = -sin(phi_b(t)).*(alpha + gam)

with `T` the transposed transform ([`applyffttranspose!`](@ref) through
`rev.fftplan`), and

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
    wcov = zeros(eltype(P), plan.n)
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
            josephsonadjoint!(P, Q, plan, wcov)

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
