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
    applydcconductance(bnm, plan, solution, Nmodes)

Subtract the direct current the resistors and blocks carry from the zero
frequency rows of the source.

The solve does not use this: it carries the average voltages as unknowns and
the currents reach the nodes through the coupling. What needs it is the
Kirchhoff validation, which reconstructs its residual from the system alone
and so knows nothing of that coupling; correcting the source here, and
adding the same current back to the residual, hands it a consistent pair.
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

# =====================================================================
# The same equation as explicit rows.
#
# Everything above eliminates the average voltages before Newton: it solves
# `Y v = j` once, folds the resistor current `G0 P v` into the zero
# frequency source, and the periodic solve never sees `v`. That is exact
# while the direct current devices are linear conductances and the sources
# are prescribed, which is the whole of stage 3.
#
# It cannot go further. A scattering block's direct current relation is a
# pencil between its port voltages and its port currents, so the current is
# a genuine unknown and there is nothing to eliminate; a short or an ideal
# through has a free current direction and no determined value at all.
# Reaching those needs `v` carried as an unknown with its equation as a row,
# which is what this builds.
#
# The row is the same component sum. Adding the zero frequency Kirchhoff
# equations over the nodes of one static flux component cancels every
# inductor and junction branch current, because both terminals of such a
# branch lie inside the component and the current enters with both signs --
# rectified mixing products included, since the cancellation is a statement
# about the topology and not about what produced the current. So the row
# sees no periodic state, and the Jacobian is block triangular:
#
#     [ Jpp  Jpv ]      Jpv = G0 P on the zero frequency nodal rows
#     [  0   Y   ]
#
# which is why an explicit `v` costs a small constant solve and does not
# make the nonlinear problem harder. Solving that triangular system by
# substitution is exactly the elimination above, which is the sense in
# which the two paths must agree and the reason the tests can demand it.
#
# A floating island fixes no absolute voltage, so `Y` is singular on it and
# only differences are physical. The elimination drops the redundant
# equation and holds the island's lowest numbered component at zero; the
# explicit form does the same, by replacing that component's row with
# `v = 0`. Keeping the singularity instead would move it into Newton, where
# it is worse, and for a linear direct current device there is nothing to be
# gained by it.

"""
    TransportRows

The transport rows `Ytr v = jtr` of the explicit direct current block,
together with the coupling into the zero frequency nodal rows.

# Fields
- `plan`: the topology, shared with the elimination.
- `Ytr`: `Y` with each floating island's pinned component replaced by the
  row `v = 0`, so it is nonsingular.
- `jtr`: the injected current, zero on a pinned row.
- `coupling`: `G0 P`, the resistor current each component's voltage drives
  into the nodes, indexed by node with ground dropped.
- `pinned`: the component held at zero on each floating island.

See [`transportrows`](@ref) and [`DCConductancePlan`](@ref).
"""
struct TransportRows{T}
    plan::DCConductancePlan
    Ytr::Matrix{T}
    jtr::Vector{T}
    coupling::SparseMatrixCSC{T,Int}
    pinned::Vector{Int}
end

nvoltages(t::TransportRows) = length(t.jtr)

"""
    transportrows(plan::DCConductancePlan, bnm, Nmodes)

Build the [`TransportRows`](@ref) for a source `bnm`.

`bnm` must be the applied source, not the one the elimination corrects:
the resistor current appears here as the coupling term, and taking it from
a corrected source would count it twice.
"""
function transportrows(plan::DCConductancePlan, bnm::AbstractVector,
        Nmodes::Integer)
    nc = size(plan.reduced, 1)
    T = float(real(eltype(plan.reduced)))

    j = zeros(T, nc)
    for (k, nodes) in enumerate(plan.components), p in nodes
        p > 1 || continue
        j[k] += real(bnm[(p-2)*Nmodes + plan.modeindex])
    end

    Ytr = Matrix{T}(plan.reduced)
    pinned = Int[]
    for (c, grounded) in plan.islands
        grounded && continue
        isempty(c) && continue
        k = first(c)                 # the elimination pins the same one
        push!(pinned, k)
        Ytr[k, :] .= zero(T)
        Ytr[k, k] = one(T)
        j[k] = zero(T)
    end

    coupling = SparseMatrixCSC{T,Int}(plan.conductance * plan.lift)
    return TransportRows(plan, Ytr, j, coupling, pinned)
end

"""
    transportresidual!(Fv, t::TransportRows, v)

The transport row residual `Ytr v - jtr`, in place.
"""
function transportresidual!(Fv::AbstractVector, t::TransportRows,
        v::AbstractVector)
    mul!(Fv, t.Ytr, v)
    Fv .-= t.jtr
    return Fv
end

"""
    transportcurrent!(d, t::TransportRows, v)

The direct resistor current `G0 P v` each node carries, in place, indexed by
node with ground dropped. This is the coupling `Jpv` applied to `v`, and it
is the quantity the elimination subtracts from the source.
"""
function transportcurrent!(d::AbstractVector, t::TransportRows,
        v::AbstractVector)
    mul!(d, t.coupling, v)
    return d
end

"""
    dcsolutionfrom(plan::DCConductancePlan, v::AbstractVector)

Package explicitly solved component voltages as a
[`DCConductanceSolution`](@ref), so the two paths report their answer the
same way.
"""
function dcsolutionfrom(plan::DCConductancePlan, v::AbstractVector)
    T = float(real(eltype(v)))
    n = size(plan.lift, 1)
    nv = zeros(T, n + 1)
    for (k, nodes) in enumerate(plan.components), p in nodes
        p > 1 && (nv[p] = v[k])
    end
    d = plan.conductance * (plan.lift * v)
    return DCConductanceSolution(any(!iszero, v), collect(T, v),
        phi0 .* nv, d)
end

"""
    dcinjected(plan::DCConductancePlan, bnm::AbstractVector, Nmodes)

Whether any direct current is injected into a static flux component.

When none is, every average voltage and every block port current is zero and
the explicit block has nothing to find: the `i = 0` rows the scattering
stamp writes are then the right answer rather than a simplification, and the
whole direct current apparatus is skipped.
"""
function dcinjected(plan::DCConductancePlan, bnm::AbstractVector,
        Nmodes::Integer)
    for nodes in plan.components, p in nodes
        p > 1 || continue
        iszero(real(bnm[(p-2)*Nmodes + plan.modeindex])) || return true
    end
    return false
end
