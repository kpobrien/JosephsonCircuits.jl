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
"""
struct DCConductancePlan{Tv,Ti}
    modeindex::Int
    dcrows::Vector{Int}
    components::Vector{Vector{Int}}
    componentof::Vector{Int}
    lift::SparseMatrixCSC{Tv,Ti}
    conductance::SparseMatrixCSC{Tv,Ti}
    reduced::Matrix{Tv}
end

"""
    DCConductanceSolution

The solved average voltages and the direct current they carry.

`nodevoltage` is in volts, indexed by node with ground first and identically
zero. On a floating island only voltage differences are physical; one
component of it is held at zero as a reference.

A solution exists exactly when the explicit direct current block was built,
which is what the caller tests by asking whether there is one. Whether the
voltages it found happen to be zero is a fact about the circuit and not
about whether it has a direct current model: a node shorted to ground sits
at zero volts, and that is an answer.
"""
struct DCConductanceSolution{T}
    nodevoltage::Vector{T}
    scaledcurrent::Vector{T}
end

"""
    dcconductanceplan(floatingcomponents, Gnm, wmodes, Nmodes, Nnodes)

Build the [`DCConductancePlan`](@ref), or `nothing` when the circuit needs
none: no zero frequency mode, or no floating static component, in which case
there is no average voltage to solve for.

Having no conductance at all is not one of those cases. A circuit whose only
direct current devices are scattering blocks has an empty `G0` and still
needs the block rows, and gating on `G0` left it with the artificial `i = 0`
rows the stamp writes.
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
    #
    # Against a tolerance and not an exact zero. These entries can come from
    # evaluating a frequency dependent law at zero, and an expression whose
    # imaginary part cancels analytically leaves roundoff when it is
    # evaluated. The scale is the largest conductance in the matrix, so the
    # decision does not depend on the units the circuit is written in.
    gscale = maximum(abs, nonzeros(G0c); init = 0.0)
    gtol = sqrt(eps(Float64))
    for g in nonzeros(G0c)
        isfinite(g) ||
            throw(ArgumentError(lazy"The zero frequency conductance of this circuit has the entry $(g), which is not finite. A component whose conductance at zero frequency is unbounded has no direct current behavior for the solver to use; give it a finite value at DC or remove the direct current drive."))
        abs(imag(g)) <= gtol*max(abs(real(g)), gscale) ||
            throw(ArgumentError(lazy"The zero frequency conductance of this circuit has the entry $(g), whose imaginary part is not roundoff. A component whose conductance at zero frequency is complex has no direct current behavior for the solver to use; give it a real value at DC or remove the direct current drive."))
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
        floatingcomponents], componentof, lift, G0, Y)
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
# only differences are physical. Something has to choose a reference, but
# nothing here does: these rows are the physics, and which of them is
# redundant is not decided by the resistors alone. A block joining an island
# to the rest of the circuit can make a row the resistors called redundant
# necessary, and a reference chosen before the blocks were assembled would
# have thrown that row away. The reference is chosen once the whole
# descriptor exists, in `dcpinning`.

"""
    TransportRows

The transport rows `Y v = j` of the explicit direct current block,
together with the coupling into the zero frequency nodal rows.

# Fields
- `plan`: the topology, shared with the elimination.
- `Y`: `P'G0P`, unreferenced. It is singular on a floating island, which is
  correct: these rows state the physics and say nothing about where the
  potential is measured from.
- `j`: the injected current.
- `coupling`: `G0 P`, the resistor current each component's voltage drives
  into the nodes, indexed by node with ground dropped.

See [`transportrows`](@ref) and [`DCConductancePlan`](@ref).
"""
struct TransportRows{T}
    plan::DCConductancePlan
    Y::Matrix{T}
    j::Vector{T}
    coupling::SparseMatrixCSC{T,Int}
end

nvoltages(t::TransportRows) = length(t.j)

"""
    transportrows(plan::DCConductancePlan, bnm, Nmodes)

Build the [`TransportRows`](@ref) for a source `bnm`.

`bnm` must be the applied source, not the one the elimination corrects:
the resistor current appears here as the coupling term, and taking it from
a corrected source would count it twice.
"""
function transportrows(plan::DCConductancePlan, bnm::AbstractVector,
        Nmodes::Integer)
    T = float(real(eltype(plan.reduced)))
    j = dcsourcecurrent(plan, bnm, Nmodes)
    coupling = SparseMatrixCSC{T,Int}(plan.conductance * plan.lift)
    return TransportRows(plan, Matrix{T}(plan.reduced), j, coupling)
end

"""
    transportresidual!(Fv, t::TransportRows, v)

The transport row residual `Y v - j`, in place.
"""
function transportresidual!(Fv::AbstractVector, t::TransportRows,
        v::AbstractVector)
    mul!(Fv, t.Y, v)
    Fv .-= t.j
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
    return DCConductanceSolution(phi0 .* nv, d)
end

"""
    dcsourcecurrent(plan::DCConductancePlan, bnm, Nmodes)

The direct current injected into each static flux component, validated.

The zero frequency mode is self conjugate, so its source coefficient is real
by construction; an imaginary part there is not a small direct current but a
sign that the source assembly or the mode layout is wrong, and taking the
real part of it would carry that error silently into the answer. It is
refused instead, against a tolerance relative to the largest coefficient of
the zero frequency source, so the decision does not depend on the units the
circuit is written in.
"""
function dcsourcecurrent(plan::DCConductancePlan, bnm::AbstractVector,
        Nmodes::Integer)
    T = float(real(eltype(plan.reduced)))
    nc = size(plan.reduced, 1)
    # the scale the imaginary parts are judged against, from the source
    # itself rather than from a fixed number
    scale = zero(T)
    for nodes in plan.components, p in nodes
        p > 1 || continue
        b = bnm[(p-2)*Nmodes + plan.modeindex]
        isfinite(b) || throw(ArgumentError(lazy"The zero frequency source at node $(p) is $(b), which is not finite. A direct current drive has to be a finite real current."))
        scale = max(scale, abs(b))
    end
    tol = sqrt(eps(T))
    j = zeros(T, nc)
    for (k, nodes) in enumerate(plan.components), p in nodes
        p > 1 || continue
        b = bnm[(p-2)*Nmodes + plan.modeindex]
        if abs(imag(b)) > tol*max(abs(real(b)), scale)
            throw(ArgumentError(lazy"The zero frequency source at node $(p) is $(b), which has an imaginary part. The zero frequency mode is self conjugate, so its source is real; an imaginary part there is a source or mode layout error rather than a current, and is refused rather than discarded."))
        end
        j[k] += real(b)
    end
    return j
end

"""
    dcinjected(plan::DCConductancePlan, bnm::AbstractVector, Nmodes)

Whether any direct current is injected into a static flux component.

When none is, every average voltage and every block port current is zero and
the explicit block has nothing to find: the `i = 0` rows the scattering
stamp writes are then the right answer rather than a simplification, and the
whole direct current apparatus is skipped.

The test is an exact zero and not a tolerance. A drive is either declared at
the zero frequency mode or it is not: `calcsources` writes the coefficient
of a mode which no source names as exactly zero, so this asks a structural
question and gets a structural answer. A small but nonzero direct current is
a direct current, and is solved for.
"""
dcinjected(plan::DCConductancePlan, bnm::AbstractVector, Nmodes::Integer) =
    any(!iszero, dcsourcecurrent(plan, bnm, Nmodes))

"""
    checkjunctiondc(sintd::AbstractArray, junctionbranches, branchnames;
        atol = 1e-2)

Warn about a Josephson junction carrying nearly its critical current at zero
frequency.

The static flux partition treats a junction of finite inductance as a short
at zero frequency, which puts its two terminals in one component and is what
lets the transport rows be the component sum of the nodal equations. That is
true of a junction in the zero voltage state and false of one which is
running, and the difference is whether the junction can carry the direct
current asked of it: the branch current is `Ic*sin(phi)`, so its zero
frequency part is `Ic` times the time average of `sin(phi)`, and no zero
voltage state exists once that average would have to exceed `Ic`.

The solver cannot report the failure itself. `sin` is bounded, so the
average it finds is always a fraction of one, and a circuit which has no
periodic solution converges to the nearest thing which is one rather than
announcing that it does not exist. What can be reported is the approach: a
junction whose direct current is within a percent of its critical current is
at the edge of the partition's assumption, and a result there should be
checked against a run at a lower drive.

This is a heuristic about the operating point and not a proof of dynamic
stability, which harmonic balance does not decide.
"""
function checkjunctiondc(sintd::AbstractArray, junctionbranches,
        branchnames::Dict{Int,Vector{String}}; atol::Real = 1e-2)
    # the last axis is the junction branch and the ones before it are the
    # time grid, which has one axis per tone; the zero frequency Fourier
    # coefficient is the average over all of them
    nt = div(length(sintd), max(size(sintd)[end], 1))
    (iszero(nt) || isempty(junctionbranches)) && return nothing
    flat = reshape(sintd, nt, :)
    for k in eachindex(junctionbranches)
        f = sum(view(flat, :, k))/nt
        abs(f) < 1 - atol && continue
        b = junctionbranches[k]
        who = join(get(branchnames, b, ["branch $(b)"]), ", ")
        @warn "A Josephson junction is carrying nearly its critical current at zero frequency, so the assumption that it is in the zero voltage state -- which is what lets the static flux partition treat it as a short at direct current -- is marginal here. A junction which cannot carry the direct current asked of it has no zero voltage state at all, and the solver cannot tell you so: the branch current is Ic*sin(phi), which is bounded, so it converges to the nearest periodic thing instead. Check this result against one at a lower drive." junction=who fractionofcritical=abs(f)
    end
    return nothing
end
