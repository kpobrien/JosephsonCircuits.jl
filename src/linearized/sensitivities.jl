# The scattering parameter sensitivities of the linearized solve: the fixed
# point stamps of each component and block, the merged parameter stamps, the
# scaling by the input waves and the contraction against the forward and
# adjoint solutions.

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
    # the directional derivatives write the workspaces of the system they
    # evaluate on, so they run on one of their own
    sys = pointsystem(op)
    dcos = similar(sys.phimatrix)
    stamps = Vector{Vector{Complex{Float64}}}(undef, size(dx,2))
    for k in axes(dx, 2)
        cosdirectionalderivative!(dcos, sys,
            Vector{Complex{Float64}}(view(dx,:,k)))
        nzval = zeros(Complex{Float64}, nnz(lsys.Asparse))
        addjosephsonterm!(nzval, lsys.complexjacobianplan, dcos)
        stamps[k] = nzval
    end
    return stamps
end

"""
    SensitivityStamp(kind, rows, cols, vals, portindex, parameter, portscale)
    SensitivityStamp(kind, rows, cols, vals, portindex)

The derivative of the linearized harmonic balance system matrix with respect
to a relative (logarithmic) perturbation of one component value, `p -> r*p`
evaluated at `r = 1`. The system matrix is affine in `C`, `1/R`, `1/L` and
`1/Lj`, so the derivative is that component's own contribution to the system
matrix, with a sign: positive for the capacitance and negative for the
quantities which enter inversely.

`kind` selects how the stamp is assembled at each signal frequency, mirroring
[`assemblesystemmatrix!`](@ref):

- `:C`: `-vals*wmodes2m`, the component's capacitance matrix,
- `:G`: `-im*vals*wmodesm`, the component's conductance matrix,
- `:invL`: `-vals`, the component's inverse inductance matrix,
- `:Lj`: the constant values, the negative of the pump modulated Josephson
    contribution of that junction alone,
- `:S`: the derivative stamp of a scattering block
    ([`blocksensitivitystamp`](@ref)).

The stamp is the triplet `(rows, cols, vals)` of the component's own
contribution, built by [`calcsensitivitystamps`](@ref), and each value is
scaled by the mode frequency of its column at assembly. `portindex` is the
index of the port whose impedance this component is, or zero, which selects
the additional wave normalization term of [`calcSsensitivity!`](@ref).

`parameter` is the design parameter this stamp contributes to, or zero for
the relative form, where every component is its own output slot.
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

# Scattering block design parameter sensitivities.
#
# A block's contribution to the system matrix is the hybrid stamp of
# B = R^(-1/2)(I - S) and C = R^(1/2)(I + S), affine in S, and the
# assembly is linear in B and C, so the derivative of the stamped values
# with respect to a parameter is the assembly of the derivatives of B and
# C, -R^(-1/2) dS/dtheta and R^(1/2) dS/dtheta, with the zero frequency
# rows, which do not depend on S, zero. It is stamped as it stands, with
# no identity part beside it to lose its digits to however small it is.
# The stamp values depend on the signal frequency through S itself, so
# unlike the lumped stamps they are rebuilt at every frequency, by each
# worker into its own buffers, because the workers of hblinsolve run
# concurrently.

"""
    checkblockpair(block, derivative, name)

Refuse a block pair, for the block instance `name`, whose block or
derivative is not a [`ScatteringParameters`](@ref): the derivative is
stamped in the block's place, and only an ordinary block states one
(see [`designblockjacobian`](@ref)); a pumped block
([`LinearizedScattering`](@ref)) has none.
"""
function checkblockpair(block, derivative, name)
    (block isa ScatteringParameters && derivative isa ScatteringParameters) ||
        throw(ArgumentError(lazy"The block pair of $(name) pairs a $(nameof(typeof(block))) with a derivative given as a $(nameof(typeof(derivative))); a block sensitivity is taken of a ScatteringParameters block, with its derivative given as one."))
    return nothing
end

"""
    targetstampsystems(ssys, target::Integer, dblock)

The stamp system of the derivative of the block contribution with
respect to a parameter of the block instance at ordinal `target` (its
position among the compiled scattering blocks, which the stamped blocks
follow), as `(dsys, position, rows, cols)`. The contributions are affine
in the block's S and no other instance's depend on the parameter, so
`dsys` holds that instance's contributions alone, with `dblock`, the
block whose S is dS/dtheta, in its place (see [`blockstampvals!`](@ref)).
`position` places each contribution among the entries of the stamp,
which are the pattern entries at `rows` and `cols` those contributions
reach. The instance is selected by its ordinal and not by its
definition, because two instances may share one definition object and a
parameter belongs to one of them.
"""
function targetstampsystems(ssys, target::Integer, dblock)
    sb = ssys.blocks[target]
    checkblockpair(sb.block, dblock, sb.name)
    cs = findall(==(target), ssys.blockindex)
    dsys = ScatteringStampSystem(
        [StampedScatteringBlock(dblock, sb.signalnodes, sb.refnodes,
            sb.auxbase, sb.name)],
        ssys.kcl, ssys.pattern, ssys.patternindex[cs], ssys.Aindex[cs],
        ones(Int32, length(cs)), ssys.pindex[cs], ssys.qindex[cs],
        ssys.coeff[cs], ssys.sign[cs], ssys.modeindex[cs],
        ssys.inmodeindex[cs], zeros(Int32, length(cs)), Int[],
        Matrix{Int}[], ssys.modeoffsets, ssys.Nmodes, ssys.Nauxports,
        ssys.scale, ssys.iscale, [1])
    # the pattern entries the contributions reach, in pattern order, and
    # their rows and columns
    entries = sort!(unique(ssys.patternindex[cs]))
    position = [searchsortedfirst(entries, k) for k in ssys.patternindex[cs]]
    colptr = SparseArrays.getcolptr(ssys.pattern)
    rows = rowvals(ssys.pattern)[entries]
    cols = [searchsortedlast(colptr, k) for k in entries]
    return dsys, position, rows, cols
end

"""
    blockstampvals!(vals, dsys, wmodes, work, buf, position)

Add to `vals` (the entries of a block parameter's stamp) the derivative
of the block contribution at the signed mode frequencies `wmodes`: the
values of the derivative system `dsys` (see
[`targetstampsystems`](@ref)), whose block's data is dS/dtheta, with its
hybrid coefficients the derivatives of `B` and `C` (see
[`evaluatehybrid!`](@ref)), scattered through `position`. The instances
of one parameter share a stamp, each through positions of its own.
`buf` is scratch at least as long as the system's contributions.
"""
function blockstampvals!(vals, dsys, wmodes, work, buf, position)
    d = view(buf, 1:length(position))
    scatteringvalues!(d, dsys, wmodes, work; derivative = true)
    @inbounds for c in eachindex(position)
        vals[position[c]] += d[c]
    end
    return vals
end

"""
    WorkerBlockSensitivity

One worker's private state for the scattering block sensitivity stamps:
its own copy of the stamp vector (the lumped stamps are shared read only,
the `:S` stamps' values are private because they are rebuilt at every
signal frequency), the derivative systems, and the evaluation scratch.
"""
struct WorkerBlockSensitivity
    stamps::Vector{SensitivityStamp}
    # (index into stamps, the derivative system, the stamp position of each
    # of the target's contributions), one per block parameter pair
    entries::Vector{Tuple{Int,Any,Vector{Int}}}
    work::ScatteringWorkspace
    buf::Vector{Complex{Float64}}
end

function WorkerBlockSensitivity(stamps::Vector{SensitivityStamp},
        blockentries)
    private = SensitivityStamp[st.kind == :S ?
        SensitivityStamp(st.kind, st.rows, st.cols, copy(st.vals),
            st.portindex, st.parameter, st.portscale) : st
        for st in stamps]
    m = maximum(e -> length(e[3]), blockentries; init = 0)
    return WorkerBlockSensitivity(private, collect(blockentries),
        ScatteringWorkspace(), Vector{Complex{Float64}}(undef, m))
end

"""
    refreshblockstamps!(wbs::WorkerBlockSensitivity, wmodes)

Rebuild this worker's `:S` stamp values at the signed mode frequencies
`wmodes`.
"""
function refreshblockstamps!(wbs::WorkerBlockSensitivity, wmodes)
    for (idx, _, _) in wbs.entries
        fill!(wbs.stamps[idx].vals, 0)
    end
    for (idx, dsys, position) in wbs.entries
        blockstampvals!(wbs.stamps[idx].vals, dsys, wmodes, wbs.work,
            wbs.buf, position)
    end
    return wbs
end

"""
    blocksensitivitystamp(rows, cols, parameter)

The [`SensitivityStamp`](@ref) skeleton of a scattering block parameter:
the entries at `rows` and `cols` its instance's contributions reach (see
[`targetstampsystems`](@ref)) with zero values, which each worker fills at
each signal frequency ([`refreshblockstamps!`](@ref)). Kind `:S` carries
no further frequency scaling, exactly as `:Lj` does.
"""
function blocksensitivitystamp(rows, cols, parameter::Integer)
    return SensitivityStamp(:S, rows, cols,
        zeros(Complex{Float64}, length(rows)), 0, Int(parameter),
        one(Complex{Float64}))
end

"""
    calcblockresidualsensitivity(op::HBOperatingPoint, psc, blockpairs)

The derivative of the harmonic balance residual with respect to each
scattering block design parameter of `blockpairs = [(blockpath,
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
matrix is stamped as the linearized stamps are (see
[`blockstampvals!`](@ref)), on a stamp system rebuilt for the pump mode grid with the geometry the
operating point already determines: the augmented dimension is the state
length and the block auxiliary variables occupy its tail.

With an explicit direct current block the columns are in canonical
coordinates, as the lumped ones are (see [`calcresidualsensitivity`](@ref)):
the harmonic column gathered and masked by the rows the block writes over,
plus the derivative of the block's own zero frequency rows
`B0 (scale dv) - C0 i` (see [`addblockdc!`](@ref)), which depend on the
parameter through `S(0)` when the block reads its zero frequency behavior
from its data rather than stating a `dcmodel`.
"""
function calcblockresidualsensitivity(op::HBOperatingPoint,
        psc::CompiledCircuit, blockpairs::AbstractVector)
    Ntot = length(op.x)
    Nmodes = op.Nmodes
    Nsc = countscatteringports(psc)*Nmodes
    pumpssys = scatteringstampsystem(psc.scatteringblocks, Nmodes;
        auxoffset = Ntot - Nsc, Ntotal = Ntot, scale = real(op.Lscale),
        modeoffsets = op.wmodes)
    isnothing(pumpssys) && throw(ArgumentError(
        "the circuit has no scattering blocks to take a block sensitivity of"))

    isrealmode = op.modelayout.isreal
    nmd = length(isrealmode)
    realindexmap = realoffsets(op)
    Nreal = realdim(Ntot, isrealmode)

    # the columns are gathered as triplets: a pair's derivative reaches the
    # rows of its own instance alone
    Ir, Jr, Vr = Int[], Int[], Float64[]
    work = ScatteringWorkspace()
    buf = Vector{Complex{Float64}}(undef, length(pumpssys.patternindex))
    dAx = zeros(Complex{Float64}, Ntot)
    for (col, bp) in enumerate(blockpairs)
        bi = scatteringblockindex(psc, bp[1])
        iszero(bi) && throw(ArgumentError(lazy"the block pair names $(bp[1]), which is not a scattering block of this circuit"))
        dsys, position, rows, cols = targetstampsystems(pumpssys, bi, bp[3])
        vals = zeros(Complex{Float64}, length(rows))
        blockstampvals!(vals, dsys, op.wmodes, work, buf, position)
        # the derivative matrix, of the target's entries alone, applied to
        # the state, over the rows those entries reach
        @inbounds for k in eachindex(vals)
            dAx[rows[k]] += vals[k]*op.x[cols[k]]
        end
        for r in unique(rows)
            v = dAx[r]
            dAx[r] = 0
            iszero(v) && continue
            kr = realindexmap[r]
            push!(Ir, kr); push!(Jr, col); push!(Vr, real(v))
            if !isrealmode[(r-1) % nmd + 1]
                push!(Ir, kr+1); push!(Jr, col); push!(Vr, imag(v))
            end
        end
    end
    dFr = sparse(Ir, Jr, Vr, Nreal, length(blockpairs))
    isnothing(op.dc) && return dFr

    L = op.dc.work.layout
    N = canonicaldim(L)
    keep = op.dc.keep
    dFc = canonicalresidual(dFr, keep, N)
    br = op.dc.work.blockrows
    isnothing(br) && return dFc
    Ic, Jc, Vc = Int[], Int[], Float64[]
    window = windowindices(L)
    u = op.dc.u
    v = view(u, voltagerange(L))
    for (col, bp) in enumerate(blockpairs)
        sb = pumpssys.blocks[scatteringblockindex(psc, bp[1])]
        dcmodelof(sb.block) isa ScatteringLimit || continue
        b = findfirst(d -> d.auxbase == sb.auxbase, br.descriptors)
        isnothing(b) && continue
        idx = br.currentindex[b]
        n = length(idx)
        dS0 = Array{Complex{Float64},3}(undef, n, n, 1)
        evaluatescattering!(dS0, bp[3], [0.0])
        r2 = sqrt.(float.(sb.block.zref))
        sc, rc = br.signalcomponent[b], br.refcomponent[b]
        # the derivative of B0 (scale dv) - C0 i, with
        # B0 = (I - S(0)) R^(-1/2) and C0 = (I + S(0)) R^(1/2)
        for p in 1:n
            row = window[idx[p]]
            acc = 0.0
            for q in 1:n
                dS = real(dS0[p, q, 1])
                acc -= dS/r2[q]*br.scale*(_vof(v, sc[q]) - _vof(v, rc[q]))
                acc -= dS*r2[q]*u[window[idx[q]]]
            end
            acc *= keep[row]
            iszero(acc) || (push!(Ic, row); push!(Jc, col); push!(Vc, acc))
        end
    end
    return dFc + sparse(Ic, Jc, Vc, N, length(blockpairs))
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
    calcsensitivitystamps(sensitivityindices, psc, nm, lsys, phimatrix,
        coupledbranches, Nmodes)

Build the [`SensitivityStamp`](@ref) of each component in
`sensitivityindices`. The classification and the entries of a capacitor,
an inductor or a resistor come from [`componentstamp`](@ref), which is
shared with the residual derivatives of [`calcresidualsensitivity`](@ref),
so the two grids cannot disagree on which components are supported or how
they are built; the node rows lead the augmented state, so the entries
need no padding. The Josephson junction stamp is the pump modulated
contribution of that junction alone, obtained by scattering the Fourier
coefficients of `cos(phi(t))` of that junction through the same plan
([`addjosephsonterm!`](@ref)) which assembles the system matrix, so the mode
coupling and its truncation agree exactly.
"""
function calcsensitivitystamps(sensitivityindices, psc::CompiledCircuit,
    nm::CircuitMatrices, lsys, phimatrix, coupledbranches, Nmodes)

    stamps = Vector{SensitivityStamp}(undef, length(sensitivityindices))
    lookups = componentlookups(coupledbranches, nm.Ljb)
    # a sensitivity taken with respect to a port's own environment also
    # moves the wave normalization, so those components are recognized by
    # their role; a port which owns no environment contributes none
    portordinal = Dict(idx => p
        for (p, idx) in enumerate(nm.portenvironmentindices) if !iszero(idx))

    kinds = [componentstamp(idx, psc, nm, lookups, Nmodes)
        for idx in sensitivityindices]
    ljstamps = junctionstamps(lsys, phimatrix, nm.Ljb,
        [s.junction for s in kinds if s.kind == :Lj], Nmodes)
    for (k, idx) in enumerate(sensitivityindices)
        portindex = get(portordinal, idx, 0)
        s = kinds[k]
        stamps[k] = if s.kind == :Lj
            rows, cols, vals = ljstamps[s.junction]
            SensitivityStamp(:Lj, copy(rows), copy(cols), copy(vals),
                portindex)
        else
            # the entries of a zero value contribute nothing
            keep = .!iszero.(s.vals)
            SensitivityStamp(s.kind, s.rows[keep], s.cols[keep],
                s.vals[keep], portindex)
        end
    end
    return stamps
end

# The stamp of each junction of `js`, as the rows, columns and values of the
# nonzero entries of the negative of the pump modulated Josephson term of
# that junction alone, in the order of the system matrix's storage. A
# junction's term fills only the entries whose row and column are both of
# its own nodes, so one assembly of each group of junctions which share no
# node (`nodedisjointgroups`) gives the stamps of all of them, each read on
# its own columns.
function junctionstamps(lsys, phimatrix, Ljb, js, Nmodes)
    stamps = Dict{Int,Tuple{Vector{Int},Vector{Int},Vector{Complex{Float64}}}}()
    isempty(js) && return stamps
    plan = lsys.complexjacobianplan
    A = lsys.Asparse
    support = plan.junctions.nodesandsigns
    nodes(j) = (node for (node, _) in support[Ljb.nzind[j]])
    groupphi = zero(phimatrix)
    jdim = ndims(phimatrix)
    nzval = zeros(Complex{Float64}, nnz(A))
    for group in nodedisjointgroups(js, nodes)
        fill!(groupphi, 0)
        for j in group
            selectdim(groupphi, jdim, j) .= selectdim(phimatrix, jdim, j)
        end
        addjosephsonterm!(nzval, plan, groupphi)
        for j in group
            I = Int[]; J = Int[]; V = Complex{Float64}[]
            for c in sort!([(n-1)*Nmodes + m for n in nodes(j) for m in 1:Nmodes])
                for t in nzrange(A, c)
                    v = nzval[t]
                    iszero(v) && continue
                    push!(I, rowvals(A)[t]); push!(J, c); push!(V, -v)
                end
            end
            stamps[j] = (I, J, V)
        end
    end
    return stamps
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
    calcsensitivityscaling!(gamma, beta, inputwave, sourcecurrents,
        portindices, portimpedances, componenttypes, wmodes, Nmodes)

Calculate the scalars which convert the derivative of the node fluxes into
the derivative of the scattering parameters. Writing the output wave of
[`calcoutputwaves!`](@ref) as a linear functional of the node fluxes, the
part which depends on them is
`(1/2)*kval*(1 + conj(Z)/Z)*im*w_n*(phi_n1 - phi_n2)`, and the scattering
parameters divide by the input wave, so

    dS[(j,n),(i,m)] = gamma[(j,n)]*beta[(i,m)]
                      *(dphi[node1,(j,n)] - dphi[node2,(j,n)])

with `gamma = (1/2)*kval*(1 + conj(Z)/Z)*im*w_n/s_{(j,n)}` and
`beta = 1/inputwave[(i,m)]`. The node flux difference is contracted
with the adjoint solution by [`calcSsensitivity!`](@ref), which is why the
`im*w_n` factor of the port voltage is folded into `gamma` here: it is
exactly the factor relating the adjoint source vector to the source vector
of the forward problem.

`s_{(j,n)}` is the source current of the port's own unit drive,
`sourcecurrents` of [`portsourcecurrents`](@ref), `±1` depending on whether the
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
function calcsensitivityscaling!(gamma, beta, inputwave, sourcecurrents,
    portindices, portimpedances, componenttypes, wmodes, Nmodes)

    for i in eachindex(portindices)
        for j in 1:Nmodes
            row = (i-1)*Nmodes + j
            portimpedance = calcimpedance(portimpedances[i],
                componenttypes[portindices[i]], wmodes[j])
            kval = portwavescale(portimpedance, wmodes[j])
            # the orientation of the port branch relative to the node order
            # of the port component
            sourcecurrent = sourcecurrents[row]
            gamma[row] = iszero(sourcecurrent) ? 0 :
                1/2*kval*(1 + conj(portimpedance)/portimpedance)*
                im*wmodes[j]/sourcecurrent
            beta[row] = iszero(inputwave[row]) ? 0 : 1/inputwave[row]
        end
    end
    return nothing
end

"""
    calcSsensitivity!(Ssensitivity, stamps, dAop, dA, dAphin, phin,
        phinadjoint, S, gamma, beta, contraction, wmodes, Nmodes)

Calculate the derivative of the scattering matrix with respect to a relative
(logarithmic) perturbation of each component value, `p -> r*p` evaluated at
`r = 1`, at the pump operating point, with the adjoint method. Overwrites
`Ssensitivity`. `stamps` are the fixed operating point stamps of each
component (see [`SensitivityStamp`](@ref)); `dAop` holds, per component,
the values of the operating point contribution to `dA` in the sparsity
structure of the linearized system matrix, or is empty when the operating
point is held fixed; `dA`, `dAphin` and `contraction` are scratch of the
size of the system matrix, the solution, and the output port mode pairs.

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
normalization of [`calcinputwaves!`](@ref) moves as well, contributing the
additional closed form term `-(portscale/2)*(P*(S+I) + (S+I)*P)` with `P`
the projector onto that port and `portscale` one in the relative form (the
design parameter form carries `(dZport/dp)/Zport`), which is exact for
constant real port impedances.
"""
function calcSsensitivity!(Ssensitivity, stamps, dAop, dA, dAphin, phin,
    phinadjoint, S, gamma, beta, contraction, wmodes, Nmodes)

    NPM = size(phin, 2)
    fill!(Ssensitivity, 0)
    for (k, stamp) in enumerate(stamps)
        # a stamp declares which output slot it accumulates into: the design
        # parameter it belongs to, or its own index in the relative form,
        # where every component is its own parameter
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
