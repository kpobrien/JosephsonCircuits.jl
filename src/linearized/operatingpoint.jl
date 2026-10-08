# Sensitivities of the scattering parameters to component values, by the
# adjoint method, at a fixed pump operating point or including the shift of
# the operating point through the implicit function theorem.

"""
    HBOperatingPoint(sys, x, jacobian, modelayout, Lscale, wmodes,
        coupledbranches, Nmodes, dc)

The converged pump operating point of [`hbnlsolve`](@ref) together with
everything needed to propagate a component perturbation through it: the
[`HBSystem`](@ref) evaluation object, the converged augmented state, the
exact Jacobian of the equivalent real system assembled there, and the
scale and layout of the augmented system.

Requested with `returnoperatingpoint = true`. The Jacobian is the exact
Jacobian of the equivalent real system, assembled with
[`assemblerealjacobian!`](@ref), rather than the complex holomorphic
Jacobian of the [`QuasiNewton`](@ref) method, which is only an approximation: the
harmonic balance residual is not complex differentiable, so the implicit
function theorem does not hold with the holomorphic Jacobian, while in the
real representation it applies directly.
"""
struct HBOperatingPoint
    # the HBSystem evaluation object, untyped: its type depends on the
    # backend
    sys
    x::Vector{Complex{Float64}}
    jacobian::SparseMatrixCSC{Float64,Int}
    modelayout::ModeLayout
    Lscale::Complex{Float64}
    wmodes::Vector{Float64}
    coupledbranches::Vector{Int}
    Nmodes::Int
    # The explicit direct current block when the circuit has one, and
    # nothing otherwise: the canonical work, the canonical state the solve
    # converged to, and the canonical Jacobian there. `jacobian` above is the
    # harmonic one, read by everything which differentiates the harmonic
    # system alone; what differentiates the whole system reads these.
    dc
end

"""
    pointsystem(op::HBOperatingPoint)

A system at the operating point with a point and workspaces of its own,
sharing the rest of `op.sys` ([`workspacetwin`](@ref)). What differentiates
an operating point evaluates on one of these rather than on `op.sys`, since
an evaluation writes the point and the workspaces of its system: the
operating point is left as it was built, and one operating point serves
concurrent calls.
"""
function pointsystem(op::HBOperatingPoint)
    sys = workspacetwin(op.sys)
    setpoint!(sys, op.x)
    return sys
end

"""
    DCOperatingPoint

The explicit direct current block at a converged point.

# Fields
- `work`: the [`CanonicalWork`](@ref) carrying the layout, the transport
  rows and the blocks' zero frequency rows.
- `u`: the converged canonical state, `[internal | vdc]`: the internal
  real state as the solver evaluates it, followed by the explicit average
  voltages (see [`CompositeLayout`](@ref)).
- `jacobian`: the canonical Jacobian there, which is the one the implicit
  function theorem applies to when the block is active.
- `keep`: the rows the block adds to rather than replaces
  ([`dckeep`](@ref)), which mask the harmonic part of a residual derivative.
"""
struct DCOperatingPoint{W}
    work::W
    u::Vector{Float64}
    jacobian::SparseMatrixCSC{Float64,Int}
    keep::Vector{Float64}
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
which carries a conductance, so they contribute nothing here, and their
derivative is the harmonic rows alone.

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
    componentstamp(idx::Integer, psc::CompiledCircuit, nm::CircuitMatrices,
        lookups, Nmodes::Integer)

Classify the component at index `idx` for sensitivity analysis and build its
raw one-component stamp, without any solver scaling, negative frequency
conjugation, or padding, which the callers apply for their own grids.
Returns a named tuple `(kind, rows, cols, vals, junction)`:

- `kind = :C`, `:G` or `:invL`: `rows`, `cols` and `vals` are the entries
    of the component's capacitance, conductance or inverse inductance
    matrix over the node flux unknowns, ground dropped and the mode
    fastest, in the order a compressed sparse column matrix stores them,
    and `junction` is zero;
- `kind = :Lj`: there are no entries, and `junction` is the ordinal of the
    Josephson junction within the junction branch vector `nm.Ljb`.

The entries are the ones the stamp plans of the system matrices
([`nodalstampplan`](@ref) and [`inverseinductanceplan`](@ref)) give a
single component, written from its own two terminals: its value on the
diagonal of each of its nodes and its negative between them, a
conductance being the reciprocal of the resistance, formed after the
sign as the plan forms it, and an inverse inductance that of the
inductance. So a stamp costs a constant, one to four entries per mode,
whatever the size of the circuit.

This is the single definition of which components are supported: `:C`, `:L`,
`:R` and `:Lj` with numeric values. Mutually coupled inductors and
components with symbolic (frequency dependent) values throw, with the same
message from both the fixed operating point stamps
([`calcsensitivitystamps`](@ref)) and the residual derivatives
([`calcresidualsensitivity`](@ref)). `lookups` are the constant time
tables of [`componentlookups`](@ref).
"""
function componentstamp(idx::Integer, psc::CompiledCircuit,
    nm::CircuitMatrices, lookups, Nmodes::Integer)

    topology = psc.topology
    componenttypes = psc.componenttypes
    nodeindices = psc.nodeindices
    vvn = nm.vvn
    componenttype = componenttypes[idx]
    value = vvn[idx]
    if !(value isa Number)
        throw(ArgumentError(lazy"Sensitivities require a numeric component value, but the value of $(psc.componentnames[idx]) is $(value). Frequency dependent values are not supported."))
    end
    n1 = nodeindices[1, idx]
    n2 = nodeindices[2, idx]
    # the storage type the value's group would assemble in: floating point
    # for a plain number however it was written, so that the reciprocal of
    # an integer value has somewhere to go
    T = grouptype(vvn, (idx,))
    v = convert(T, value)
    if componenttype == :C
        rows, cols, vals = twoterminalstamp(n1, n2, v, -v, Nmodes)
        return (kind = :C, rows, cols, vals, junction = 0)
    elseif componenttype == :R
        rows, cols, vals = twoterminalstamp(n1, n2, 1/v, 1/(-v), Nmodes)
        return (kind = :G, rows, cols, vals, junction = 0)
    elseif componenttype == :L
        b = topology.edge2indexdict[(n1, n2)]
        if b in lookups.coupled
            throw(ArgumentError(lazy"Sensitivities are not supported for the mutually coupled inductor $(psc.componentnames[idx])."))
        end
        y = 1/v
        rows, cols, vals = twoterminalstamp(n1, n2, y, -y, Nmodes)
        return (kind = :invL, rows, cols, vals, junction = 0)
    elseif componenttype == :Lj
        b = topology.edge2indexdict[(n1, n2)]
        j = get(lookups.junctionordinal, b, nothing)
        if isnothing(j)
            throw(ArgumentError(lazy"The Josephson junction $(psc.componentnames[idx]) was not found in the branch inductance vector."))
        end
        return (kind = :Lj, rows = Int[], cols = Int[],
            vals = Complex{Float64}[], junction = j)
    else
        throw(ArgumentError(lazy"Sensitivities are only supported for C, L, R, and Lj components, not $(componenttype), the type of $(psc.componentnames[idx])."))
    end
end

# The entries of a two terminal component between the nodes `n1` and `n2`
# (`1` being ground) whose stamp is `d` on the diagonal of each of its
# nodes and `o` between them, repeated for each of `Nmodes` modes, as the
# rows, the columns and the values over the node flux unknowns with ground
# dropped and the mode fastest, in the order a compressed sparse column
# matrix stores them: the columns of the lower node and then of the higher,
# each with its rows in order. A grounded component has one entry per mode,
# a floating one four, and one whose terminals are one node, which carries
# no current, none.
function twoterminalstamp(n1::Integer, n2::Integer, d, o, Nmodes::Integer)
    rows = Int[]; cols = Int[]; vals = Complex{Float64}[]
    n1 == n2 && return rows, cols, vals
    dc = Complex{Float64}(d)
    if n1 == 1 || n2 == 1
        a = max(n1, n2) - 1
        sizehint!(rows, Nmodes); sizehint!(cols, Nmodes); sizehint!(vals, Nmodes)
        for m in 1:Nmodes
            k = (a-1)*Nmodes + m
            push!(rows, k); push!(cols, k); push!(vals, dc)
        end
    else
        oc = Complex{Float64}(o)
        lo, hi = minmax(n1, n2) .- 1
        sizehint!(rows, 4Nmodes); sizehint!(cols, 4Nmodes); sizehint!(vals, 4Nmodes)
        for (c, first, second) in ((lo, dc, oc), (hi, oc, dc))
            for m in 1:Nmodes
                col = (c-1)*Nmodes + m
                push!(rows, (lo-1)*Nmodes + m); push!(cols, col); push!(vals, first)
                push!(rows, (hi-1)*Nmodes + m); push!(cols, col); push!(vals, second)
            end
        end
    end
    return rows, cols, vals
end

# the offset of each complex entry of the operating point in its real
# representation, one real row for a self conjugate mode and a real and an
# imaginary row otherwise: the first slot of each entry in the state's
# layout, which spans the whole augmented state
realoffsets(op::HBOperatingPoint) =
    Vector{Int}(op.modelayout.ptr[1:length(op.x)])

"""
    calcresidualsensitivity(op::HBOperatingPoint, psc, nm,
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
modes, or a junction branch's Kirchhoff rows), so the storage scales with
the touched entries rather than with `Nstate*Ncomponents`.

The residual is affine in `C`, `1/R` and `1/L`, so those parameter
derivatives are that component's own contribution to the linear term applied
to the converged state, with a sign, built from the shared classification of
[`componentstamp`](@ref) with the same nondimensionalization and negative
frequency conjugation the solver applies. The Josephson junction term is
the residual's own sine contribution restricted to that junction.

Note that the auxiliary unknowns of the modified nodal analysis
formulation, the currents of the mutually coupled inductors and of the
scattering block ports, are scaled by the solver scale (see
[`calcsolverscale`](@ref)), which itself depends on the port impedances, so
for a port's termination the auxiliary rows of the returned derivative
differ from a finite difference of a re-solve by that change of
normalization. The node flux rows, which are the physical quantity and the
only rows the linearized system depends on, are unaffected.
"""
function calcresidualsensitivity(op::HBOperatingPoint,
    psc::CompiledCircuit, nm::CircuitMatrices,
    sensitivityindices,
    alphas::AbstractVector = ones(Complex{Float64},
        length(sensitivityindices)))

    Ntot = length(op.x)
    Nmodes = op.Nmodes
    lookups = componentlookups(op.coupledbranches, op.sys.Ljb)
    stamps = [sensitivitystamp(idx, psc, nm, lookups, Nmodes)
        for idx in sensitivityindices]
    # the Josephson terms of the junctions asked for, on their own rows,
    # evaluated on a system of its own at the operating point, whose cached
    # time domain branch fluxes they read
    js = [s.junction for s in stamps if s.kind == :Lj]
    ljterms = isempty(js) ?
        Dict{Int,Vector{Tuple{Int,Complex{Float64}}}}() :
        josephsonterms(pointsystem(op), js, Nmodes)

    # The residual derivatives are sparse: each component touches only its
    # own rows of the state (its nodes and modes, or one junction branch's
    # Kirchhoff rows), so they are accumulated as triplets of the real
    # representation and returned as a sparse matrix, where a dense array
    # would cost O(Nstate*Ncomponents), the many component regime the
    # reverse contraction order exists for. Duplicate triplets (a component
    # with several entries in one row) are summed by `sparse`, which is
    # right because the real representation is linear.
    isrealmode = op.modelayout.isreal
    Ir = Int[]; Jc = Int[]; Vr = Float64[]
    residualentries!(Ir, Jc, Vr, stamps, ljterms, alphas, op.x, op.wmodes,
        op.Lscale, realoffsets(op), isrealmode)
    harmonic = sparse(Ir, Jc, Vr, realdim(Ntot, isrealmode),
        length(sensitivityindices))
    isnothing(op.dc) && return harmonic

    # In canonical coordinates the residual is `D G F(S u) + M u`, so its
    # parameter derivative is this gathered and masked by the rows the
    # block replaces rather than adds to, plus the block's own dependence
    # on the component values, which lands in the voltage rows the gather
    # leaves alone.
    dccols = dcresidualsensitivity(op.dc, psc, nm, op.Lscale,
        sensitivityindices, alphas)
    return canonicalresidual(harmonic, op.dc.keep,
        canonicaldim(op.dc.work.layout)) + dccols
end

# The triplets of the residual derivative columns, in the real
# representation of the augmented state, for the component stamps of
# `componentstamp` and the Josephson terms of `josephsonterms`: the
# component's own contribution to the linear term, `c*vals*w^power`
# applied to the state `x` column by column, the value directed by
# `alphas`, conjugated at the negative frequency modes and
# nondimensionalized by `Lscale` as the solver does with the full
# matrices, and the negative of a junction's term.
function residualentries!(Ir, Jc, Vr, stamps, ljterms, alphas, x, wmodes,
        Lscale, realindexmap, isrealmode)
    Nm = length(wmodes)
    nmd = length(isrealmode)
    function pushentry!(r, comp, v)
        kr = realindexmap[r]
        push!(Ir, kr); push!(Jc, comp); push!(Vr, real(v))
        if !isrealmode[(r-1) % nmd + 1]
            push!(Ir, kr+1); push!(Jc, comp); push!(Vr, imag(v))
        end
        return nothing
    end
    for (comp, s) in enumerate(stamps)
        alpha = alphas[comp]
        if s.kind == :Lj
            for (r, v) in ljterms[s.junction]
                pushentry!(r, comp, -v*alpha)
            end
            continue
        end
        c, power = s.kind == :C ? (-1.0 + 0im, 2) :
            s.kind == :G ? (0.0 - 1im, 1) : (-1.0 + 0im, 0)
        for t in eachindex(s.vals)
            col = s.cols[t]
            xc = x[col]
            iszero(xc) && continue
            w = wmodes[(col-1) % Nm + 1]
            scale = power == 0 ? c*xc : power == 1 ? c*w*xc : c*w^2*xc
            # the design parameter rescale is a direction in component
            # value space, so it multiplies the stored value before the
            # negative frequency conjugation, exactly as in `reparameterize`
            v = s.vals[t]
            isone(alpha) || (v *= alpha)
            pushentry!(s.rows[t], comp, modevalue(v, w)*Lscale*scale)
        end
    end
    return nothing
end

"""
    canonicalresidual(dF::SparseMatrixCSC, keep::Vector{Float64}, N)

The harmonic residual derivative columns `dF` in the `N` canonical
coordinates of the direct current block: the gather copies the harmonic
rows, which lead the canonical vector, and the rows the block replaces
rather than adds to are masked out by `keep` ([`dckeep`](@ref)). Linear in
the stored entries of `dF`; the block's own dependence on the parameters
is added by the caller.
"""
function canonicalresidual(dF::SparseMatrixCSC, keep::Vector{Float64},
        N::Integer)
    I = Int[]; J = Int[]; V = Float64[]
    rows = rowvals(dF)
    vals = nonzeros(dF)
    for k in axes(dF, 2), p in nzrange(dF, k)
        v = vals[p]*keep[rows[p]]
        iszero(v) && continue
        push!(I, rows[p]); push!(J, k); push!(V, v)
    end
    return sparse(I, J, V, N, size(dF, 2))
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

# The junctions of `js` in groups of which no two share a node, `nodes(j)`
# giving the nodes of junction `j`, greedily: a junction joins the first
# group none of whose junctions meets it. A junction's contribution to the
# residual, or to the linearized system matrix, reaches only the rows and
# columns of its own nodes, so the junctions of a group are evaluated
# together and each reads its own part; the number of groups is set by how
# many junctions meet at a node, not by how many there are.
function nodedisjointgroups(js, nodes)
    groups = Vector{Int}[]
    nodegroups = Dict{Int,BitSet}()
    for j in unique(js)
        taken = BitSet()
        for n in nodes(j)
            union!(taken, get(nodegroups, n, BitSet()))
        end
        g = 1
        while g in taken
            g += 1
        end
        g > length(groups) && push!(groups, Int[])
        push!(groups[g], j)
        for n in nodes(j)
            push!(get!(BitSet, nodegroups, n), g)
        end
    end
    return groups
end

# The Josephson contribution of the residual of each junction of `js`
# alone: the node vector of the Fourier coefficients of sin(phi_b(t))/Lj of
# that junction, which is the derivative of the residual with respect to
# minus the logarithm of its inductance, as the pairs of row and value on
# the Kirchhoff rows of its branch, the only rows it reaches. The transform
# is taken once and one backward map of each group of junctions which
# share no node (`nodedisjointgroups`) gives the terms of all of them.
function josephsonterms(sys, js, Nmodes)
    terms = Dict{Int,Vector{Tuple{Int,Complex{Float64}}}}()
    isempty(js) && return terms
    support = branchnodesandsigns(sys.Rbnm, Nmodes, size(sys.Rbnm, 1) ÷ Nmodes)
    nodes(j) = (node for (node, _) in support[sys.Ljb.nzind[j]])
    groups = nodedisjointgroups(js, nodes)
    _ensuresin!(sys)
    applyfft!(sys.phimatrix, sys.sintd, sys.rfftplan)
    sinfd = copy(sys.phimatrix)
    jdim = ndims(sinfd)
    out = zeros(Complex{Float64}, length(sys.x))
    for group in groups
        fill!(sys.phimatrix, 0)
        for j in group
            selectdim(sys.phimatrix, jdim, j) .= selectdim(sinfd, jdim, j)
        end
        # the Josephson contribution alone, without the linear term
        applybackwardterm!(out, sys.nonlineartermplan, sys.phimatrix, sys.x;
            addlinearterm = false)
        for j in group
            terms[j] = [(r, out[r]) for n in nodes(j)
                for r in (n-1)*Nmodes+1:n*Nmodes]
        end
    end
    return terms
end

"""
    ReverseSensitivity(op, sys, dFr, T, fftplan, nzrow, nzcol,
        realindexmap, branchnodes, slots)

Everything the reverse mode contraction of [`calcSsensitivityreverse!`](@ref)
needs, precomputed once and shared read only across the signal frequencies.
"""
struct ReverseSensitivity
    op::HBOperatingPoint
    # a system of its own at the operating point, whose cached second
    # derivative of the relation the contraction reads at every frequency
    sys
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
itself is shared read only. The time domain arrays of the transposed
transform have one dimension per pump tone and one per junction, and
their type is the parameter, so that the contraction reads them
concretely.
"""
struct ReverseSensitivityBuffers{A<:Array{Complex{Float64}}}
    P::Vector{Complex{Float64}}
    Q::Vector{Complex{Float64}}
    # the output functional covectors of a chunk of output pairs, and their
    # solutions through the transposed pump Jacobian, batched so the sparse
    # solver amortizes its per-call overhead over many right hand sides.
    # The covectors are complex and the Jacobian real, so a chunk of `n`
    # pairs is held as `2n` real columns, the real parts of its covectors
    # followed by their imaginary parts, which one real solve takes as
    # they are, where a complex right hand side would be split into real
    # copies at every solve.
    G::Matrix{Float64}
    Psi::Matrix{Float64}
    eta::Matrix{Complex{Float64}}
    c::Matrix{Complex{Float64}}
    # the zero padded input and the single output grid of the transposed
    # transform: the holomorphic and antiholomorphic halves are folded into
    # eta one after the other, so their transforms need not coexist.
    padded::A
    tgrid::A
    # the covector of one output pair per stored entry of the Josephson plan
    wcov::Vector{Complex{Float64}}
end

# The byte budget of the right hand side and solution buffers of one
# frequency batch of the reverse contraction, `nr` by `2*chunk` real
# each for a chunk of `chunk` output pairs. Together they cost
# `32*nr*chunk` bytes, so a pump system of 8000 real unknowns gets a
# chunk of about 128 output pairs and a system a hundred times
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
    sys = rev.sys
    NLj = size(sys.phimatrix)[end]
    NF = length(sys.phimatrix)
    nr = sensitivitydim(rev.op)
    chunk = clamp(REVERSESENSITIVITYCHUNKBYTES ÷ (32*nr), 1, NPM^2)
    return ReverseSensitivityBuffers(
        zeros(Complex{Float64}, NF), zeros(Complex{Float64}, NF),
        zeros(Float64, nr, 2*chunk),
        zeros(Float64, nr, 2*chunk),
        zeros(Complex{Float64}, size(rev.T, 1), NLj),
        zeros(Complex{Float64}, size(rev.T, 2), NLj),
        zeros(Complex{Float64}, size(sys.phitd)),
        zeros(Complex{Float64}, size(sys.phitd)),
        zeros(Complex{Float64}, length(rev.nzrow)))
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
    ReverseSensitivity(op::HBOperatingPoint, lsys, dFr, slots)

Precompute the reverse mode contraction data: the branch flux map
([`calcbranchtimedomainmap`](@ref)), the transform of the pump harmonic grid,
the row and the column of each nonzero of the linearized system matrix, and
the offset of each entry of the augmented state in its real representation.
`dFr` holds the residual derivative columns and `slots[k]` is the output
slot of `Ssensitivity` column `k` accumulates into, so that several columns
belonging to one design parameter can share a slot.
"""
function ReverseSensitivity(op::HBOperatingPoint, lsys, dFr,
        slots::Vector{Int})
    # the operating point, and therefore the branch flux map and the
    # incidence lists, live on the pump mode grid, not the signal mode grid
    Nmodes = op.Nmodes
    # the per frequency contraction reads the cached negative of the
    # second derivative of the relation at the branch fluxes, `sin` for the
    # Josephson one, so it is pinned to the operating point here, in this
    # serial constructor, on a system of its own: the threads of hblinsolve
    # only read it, and updating it from them would be a race
    sys = pointsystem(op)
    _negsecond!(sys)
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
    realindexmap = realoffsets(op)
    Nbranches = size(sys.Rbnm, 1) ÷ Nmodes
    return ReverseSensitivity(op, sys, SparseMatrixCSC{Float64,Int}(dFr),
        calcbranchtimedomainmap(sys, Nmodes, NLj),
        plan_applyffttranspose(sys.phitd), nzrow, nzcol,
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
`rev.fftplan`), and its covector is the transpose of the branch flux map applied to `eta`.
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

    # The operating point holds its system, its layout and its direct current
    # block untyped, and so does the plan of the transform, so what the loops
    # read of them is read here once and the loops run behind a function
    # barrier, compiled for the concrete types.
    op = rev.op
    sys = rev.sys
    NLj = size(sys.phimatrix)[end]
    # the cached negative of the second derivative of the relation at the
    # pump branch fluxes, `sin` for the Josephson one, pinned to the
    # operating point by the ReverseSensitivity constructor and read only
    # here.
    sintd = reshape(_negsecond!(sys), :, NLj)
    return reversecontraction!(Ssensitivity, rev.dFr, rev.T, rev.fftplan,
        rev.nzrow, rev.nzcol, rev.realindexmap, rev.branchnodes, rev.slots,
        lsys.complexjacobianplan, op.Nmodes, op.modelayout.isreal, sintd,
        size(sys.phimatrix), sys.Ljb.nzind, phin, phinadjoint, gamma, beta,
        cache, bufs)
end

# the loops of `calcSsensitivityreverse!`, with every argument concrete
function reversecontraction!(Ssensitivity, dFr, T, fftplan, nzrow, nzcol,
    realindexmap, branchnodes, slots, plan, Nmodes, isrealmode, sintd,
    wsize, junctionbranches, phin, phinadjoint, gamma, beta, cache,
    bufs::ReverseSensitivityBuffers)

    NPM = size(phin, 2)
    NLj = wsize[end]
    Ncomponents = size(dFr, 2)
    nmd = length(isrealmode)
    P = bufs.P
    Q = bufs.Q
    G = bufs.G
    Psi = bufs.Psi
    eta = bufs.eta
    c = bufs.c
    padded = bufs.padded
    tgrid = bufs.tgrid
    dFrows = rowvals(dFr)
    dFvals = nonzeros(dFr)

    # The output pairs are independent, so their solves through the
    # transposed pump Jacobian are batched: the covectors of a chunk of
    # pairs are accumulated as the columns of G, their real parts and then
    # their imaginary parts, and pushed through the factorization in one
    # multi right hand side call, which amortizes the per-call overhead of
    # the sparse triangular solves over the chunk.
    wcov = bufs.wcov
    pairs = vec(CartesianIndices((NPM, NPM)))
    for chunk in Iterators.partition(eachindex(pairs), size(G, 2) ÷ 2)
        ncols = length(chunk)
        for (col, pi) in enumerate(chunk)
            a, b = Tuple(pairs[pi])
            # the transpose of the Josephson scatter
            fill!(P, 0)
            fill!(Q, 0)
            # the covector per destination, then the transpose of the
            # Josephson map applied to it. The plain and conjugated halves
            # separate by the sign of the mode coupling index.
            @inbounds for p in eachindex(wcov)
                wcov[p] = phinadjoint[nzrow[p],a]*phin[nzcol[p],b]
            end
            josephsonadjoint!(P, Q, plan, wcov)

            # the transpose of the cos directional derivative: a forward
            # transform of the zero padded coefficients through the plan of
            # the operating point grid (the transform matrix is symmetric),
            # and the antiholomorphic half by conjugating around the same
            # transform. the two halves are folded into eta one after the
            # other through the same output grid.
            applyffttranspose!(tgrid, reshape(P, wsize), padded, fftplan)
            tf = reshape(tgrid, :, NLj)
            @inbounds for i in eachindex(eta)
                eta[i] = -sintd[i]*tf[i]
            end
            @inbounds for i in eachindex(Q)
                Q[i] = conj(Q[i])
            end
            applyffttranspose!(tgrid, reshape(Q, wsize), padded, fftplan)
            @inbounds for i in eachindex(eta)
                eta[i] -= sintd[i]*conj(tf[i])
            end
            mul!(c, transpose(T), eta)

            # the transpose of the branch flux map, into the real
            # representation of the augmented state, its real part in
            # column `col` and its imaginary part in column `ncols + col`.
            # The outputs depend on the state through the harmonic
            # coordinates, so this is where the covector is built whether or
            # not there is a direct current block. The canonical coordinates
            # of one lead with the harmonic ones, so the covector is already
            # carried into them, which is the transpose of the scatter the
            # residual applies going the other way; the voltage rows stay
            # zero, since no output reads an average voltage directly, only
            # through the state the block moves.
            gre = view(G, :, col)
            gim = view(G, :, ncols + col)
            fill!(gre, 0)
            fill!(gim, 0)
            @inbounds for (jj, branch) in enumerate(junctionbranches)
                for m in 1:Nmodes
                    # the covector of the real and of the imaginary part of
                    # the mode's branch flux
                    cre = c[2*(m-1)+1, jj]
                    cim = c[2*(m-1)+2, jj]
                    for (node, sgn) in branchnodes[branch]
                        j = (node-1)*Nmodes + m
                        k = realindexmap[j]
                        gre[k] += sgn*real(cre)
                        gim[k] += sgn*imag(cre)
                        if !isrealmode[(j-1) % nmd + 1]
                            gre[k+1] += sgn*real(cim)
                            gim[k+1] += sgn*imag(cim)
                        end
                    end
                end
            end
        end

        # through the transposed Jacobian once for the whole chunk
        Gc = view(G, :, 1:2*ncols)
        Psic = view(Psi, :, 1:2*ncols)
        trysolvetranspose!(Psic, cache.factorization, Gc)

        # the sparse inner products with the residual derivatives
        for (col, pi) in enumerate(chunk)
            a, b = Tuple(pairs[pi])
            @inbounds for k in 1:Ncomponents
                accre = 0.0
                accim = 0.0
                for r in nzrange(dFr, k)
                    row = dFrows[r]
                    v = dFvals[r]
                    accre += Psi[row, col]*v
                    accim += Psi[row, ncols + col]*v
                end
                Ssensitivity[a,b,slots[k]] +=
                    gamma[a]*beta[b]*complex(accre, accim)
            end
        end
    end
    return nothing
end
