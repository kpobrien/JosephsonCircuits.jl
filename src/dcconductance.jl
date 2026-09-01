# Direct current through resistors.
#
# The harmonic balance state is periodic node flux, so a voltage is its time
# derivative and appears as `V = i*w*phi0*phi`. At zero frequency that is
# identically zero, and the linear term
#
#     K = invLnm + im*Gnm*wmodesm - Cnm*wmodes2m
#
# carries no conductance in its direct current rows: a resistor is an open
# circuit at DC. That is right for a capacitor and wrong for a resistor, so a
# current source driving a resistor could not develop I*R.
#
# The missing coordinate is the linear in time part of the flux, equivalently
# the average node voltage:
#
#     Phi(t) = phi0*phitilde(t) + P*v*t
#
# with `phitilde` periodic and `v` the average voltages. A finite inductor or
# a zero-voltage Josephson junction requires zero average voltage across its
# branch, so `v` is constant on each connected component of the finite L/Lj
# graph -- the components `calcstaticfluxcomponents` already returns.
#
# Summing the DC Kirchhoff equation over the nodes of a component cancels
# every inductor and junction branch current, because both of its terminals
# lie inside the component and the current appears with both signs. That
# holds whatever produced the current, rectification of the drive by the
# junction nonlinearity included, so what remains is linear in the average
# voltages alone:
#
#     Y v = j,   Y = P'G0P,   j = P'i_source,0
#
# For linear conductances and prescribed current sources this is independent
# of the nonlinear state, so it is solved once, before Newton, and its
# resistor current folded into the zero frequency source. That is an exact
# elimination rather than an approximation, and it changes no mode layout, no
# transform plan, no residual or Jacobian kernel, and no backend path.
#
# The conductance is read from the assembled `Gnm` rather than from the
# resistors, so that it keeps working when a port environment becomes a
# boundary stamp with no component behind it.

"""
    DCConductancePlan

The topology of the eliminated direct current voltage block.

Built from the circuit and its assembled conductance; holds nothing which
depends on the sources.

# Fields
- `modeindex`, `dcrows`: the zero frequency mode and its nodal rows.
- `components`: the nodes of each floating static flux component.
- `componentof`: the component index of each node, zero for a node whose
    average voltage is fixed at zero by a path to ground through inductance.
- `lift`: `P`, mapping a component voltage to its nodes.
- `conductance`: `G0`, the direct current conductance in the solver's scaled
    units.
- `reduced`: `Y = P'G0P`, the conductance seen between components.
- `islands`: the resistively connected groups of components, and whether
    each reaches ground.
"""
struct DCConductancePlan{Tv,Ti}
    modeindex::Int
    dcrows::Vector{Int}
    components::Vector{Vector{Int}}
    componentof::Vector{Int}
    lift::SparseMatrixCSC{Tv,Ti}
    conductance::SparseMatrixCSC{Tv,Ti}
    reduced::Matrix{Tv}
    islands::Vector{Tuple{Vector{Int},Bool}}
end

"""
    DCConductanceSolution

The solved average voltages and the direct current they carry.

`nodevoltage` is in volts, indexed by node with ground first and identically
zero. On a floating island only voltage differences are physical; the
island's lowest numbered component is held at zero.
"""
struct DCConductanceSolution{T}
    active::Bool
    componentvoltage::Vector{T}
    nodevoltage::Vector{T}
    scaledcurrent::Vector{T}
end

"""
    dcconductanceplan(floatingcomponents, Gnm, wmodes, Nmodes, Nnodes)

Build the [`DCConductancePlan`](@ref), or `nothing` when the circuit needs
none: no zero frequency mode, no floating static component, or no
conductance joining two distinct components, in which case no resistor can
carry direct current between components and the existing answer is already
right.
"""
function dcconductanceplan(floatingcomponents::Vector{Vector{Int}},
        Gnm::SparseMatrixCSC, wmodes::AbstractVector, Nmodes::Integer,
        Nnodes::Integer)

    m0 = findfirst(iszero, wmodes)
    isnothing(m0) && return nothing
    isempty(floatingcomponents) && return nothing

    n = Nnodes - 1
    dcrows = collect(Int(m0):Nmodes:n*Nmodes)
    G0c = Gnm[dcrows, dcrows]
    # A conductance carries direct current only if it is real: a complex one
    # at zero frequency has no steady state meaning, and would come from a
    # frequency dependent law whose limit at DC is not a conductance. Refuse
    # it rather than take a part of it.
    for (k, g) in enumerate(nonzeros(G0c))
        if !isfinite(g) || !iszero(imag(g))
            throw(ArgumentError(lazy"The zero frequency conductance of this circuit has the entry $(g), which is not finite and real. A component whose conductance at zero frequency is complex or unbounded has no direct current behavior for the solver to use; give it a real finite value at DC or remove the direct current drive."))
        end
    end
    G0 = SparseMatrixCSC(size(G0c)..., G0c.colptr, G0c.rowval,
        real.(nonzeros(G0c)))

    # node -> component, zero for ground and for anything held at zero by an
    # inductive path to it
    componentof = zeros(Int, Nnodes)
    for (k, nodes) in enumerate(floatingcomponents), p in nodes
        componentof[p] = k
    end
    nc = length(floatingcomponents)

    # The conductance seen between components, Y = P'G0P. A conductance to
    # ground appears in G0 only as a diagonal entry, because the ground node
    # has no row of its own; a floating one contributes +g to two diagonals
    # and -g to two off diagonals. So the row sum of G0 at a node is exactly
    # its conductance to ground, which is what tells a grounded component
    # from a floating one.
    Y = zeros(eltype(G0), nc, nc)
    rows = rowvals(G0)
    vals = nonzeros(G0)
    any(!iszero, vals) || return nothing
    for col in 1:n
        b = componentof[col+1]
        for r in nzrange(G0, col)
            a = componentof[rows[r]+1]
            (a > 0 && b > 0) && (Y[a,b] += vals[r])
        end
    end

    lift = sparse(
        [p-1 for c in floatingcomponents for p in c if p > 1],
        [k for (k, c) in enumerate(floatingcomponents) for p in c if p > 1],
        ones(eltype(G0), count(p -> p > 1, reduce(vcat, floatingcomponents))),
        n, nc)


    return DCConductancePlan(Int(m0), dcrows, [Int.(c) for c in
        floatingcomponents], componentof, lift, G0, Y,
        dcislands(G0, componentof, n, nc))
end

# The resistively connected groups of components, and whether each reaches
# ground. Current injected into a grounded island leaves through ground, so
# no compatibility condition applies to it; a floating island supports a
# bounded solution only when its net injection is zero.
function dcislands(G0::SparseMatrixCSC, componentof::Vector{Int}, n::Int,
        nc::Int)
    parent = collect(0:nc)
    findroot(i) = (while parent[i+1] != i; i = parent[i+1]; end; i)
    rows = rowvals(G0)
    vals = nonzeros(G0)
    rowsum = zeros(eltype(G0), n)
    for col in 1:n
        b = componentof[col+1]
        for r in nzrange(G0, col)
            iszero(vals[r]) && continue
            row = rows[r]
            rowsum[row] += vals[r]
            a = componentof[row+1]
            ra, rb = findroot(a), findroot(b)
            ra == rb || (parent[ra+1] = rb)
        end
    end
    # a node whose conductances do not cancel in its row sum has one to
    # ground, which joins its component to the ground island
    for p in 1:n
        iszero(rowsum[p]) && continue
        a = componentof[p+1]
        a == 0 && continue
        ra, rg = findroot(a), findroot(0)
        ra == rg || (parent[ra+1] = rg)
    end
    groups = Dict{Int,Vector{Int}}()
    for c in 1:nc
        push!(get!(Vector{Int}, groups, findroot(c)), c)
    end
    gr = findroot(0)
    return [(sort(v), r == gr) for (r, v) in sort(collect(groups); by = first)]
end

"""
    solvedcconductance(plan, bnm, Nmodes, Lmean)

Solve `Y v = j` for the average component voltages given the applied source.

A grounded island solves directly. A floating island holds its lowest
numbered component at zero and drops that equation, which is the redundant
one; only its voltage differences are physical. An island whose net injected
current is not zero has no bounded solution, and the compatibility check
which follows reports it.

Returns an inactive solution when no direct current is injected anywhere,
which is the common case and costs one projection.
"""
function solvedcconductance(plan::DCConductancePlan, bnm::AbstractVector,
        Nmodes::Integer, Lmean)
    nc = size(plan.reduced, 1)
    T = float(real(eltype(plan.reduced)))

    # the direct current injected into each component, in the solver's units
    j = zeros(T, nc)
    for (k, nodes) in enumerate(plan.components), p in nodes
        p > 1 || continue
        j[k] += real(bnm[(p-2)*Nmodes + plan.modeindex])
    end

    n = size(plan.lift, 1)
    if all(iszero, j)
        return DCConductanceSolution(false, zeros(T, nc),
            zeros(T, n+1), zeros(T, n))
    end

    Yr = plan.reduced
    v = zeros(T, nc)
    for (c, grounded) in plan.islands
        if grounded
            v[c] = Yr[c, c] \ j[c]
        else
            keep = c[2:end]
            isempty(keep) && continue
            v[keep] = Yr[keep, keep] \ j[keep]
        end
    end

    nv = zeros(T, n+1)
    for (k, nodes) in enumerate(plan.components), p in nodes
        p > 1 && (nv[p] = v[k])
    end
    d = plan.conductance * (plan.lift * v)
    # The scaled coordinate is v = V/phi0: the source carries a factor
    # Lmean/phi0 and the conductance a factor Lmean, which cancel, so the
    # physical voltage is phi0 times the solved coordinate. The current d is
    # left in the source's units, which is what it is subtracted from.
    return DCConductanceSolution(true, v, phi0.*nv, d)
end

"""
    applydcconductance(bnm, plan, solution, Nmodes)

Subtract the direct resistor current from the zero frequency rows of the
source.

This is the elimination: the residual the periodic flux solver then forms is
the physical augmented Kirchhoff equation, exactly and not approximately.
The caller keeps the original source, because the port and scattering
outputs need the current which was applied rather than the corrected one.
"""
function applydcconductance(bnm::AbstractVector, plan::DCConductancePlan,
        sol::DCConductanceSolution, Nmodes::Integer)
    sol.active || return bnm
    out = copy(bnm)
    for p in 2:(size(plan.lift, 1) + 1)
        out[(p-2)*Nmodes + plan.modeindex] -= sol.scaledcurrent[p-1]
    end
    return out
end
