# Direct current in the harmonic balance system.
#
# The zero frequency mode of the periodic state is the static node flux, and
# a voltage is its time derivative, so the periodic state alone carries no
# average voltage: a resistor is an open circuit there and a scattering
# block sees `i = 0`. This file is everything the solver does about that,
# in the order it happens: the gauge which pins the static flux of every
# floating inductive subnetwork; the average voltages carried as unknowns
# beside the periodic state, with the resistor conductance as their
# equations and the resistor current as their coupling into the nodal rows;
# the zero frequency pencil of each scattering block, whose port currents
# the solver already carries; the direct current subsystem those rows form,
# its classification (which redundant equations to give up, and when a free
# direction is refused), its exact factorization where the state lives; and
# the block as a matrix, so that its contribution to the residual and the
# Jacobian-vector product is three array operations on the window.
#
# The canonical state layout these rows are written against, `[internal |
# vdc]`, is in harmonics/layout.jl; the canonical operators and the
# preconditioner wrapper which apply them are in harmonics/canonical.jl.

# =====================================================================
# The gauge at zero frequency: the static flux of a floating inductive
# subnetwork is fixed to a reference, one node per subnetwork.

"""
    calcstaticfluxcomponents(componenttypes::Vector{Symbol},
        nodeindices::Matrix{Int}, vvn::Vector, Nnodes::Int)

Return the connected components of the static flux-stiffness graph which do
not contain the ground node, as a vector of vectors of "one indexed" node
indices (so 1 is the ground node and does not appear in the output). The graph
has the circuit nodes as vertices and an edge for every linear inductor and
every Josephson junction whose value provides static (zero frequency) flux
stiffness: finite numeric values contribute an edge, while a non-finite
numeric value does not. Symbolic values are assumed to provide finite, nonzero
static stiffness; a symbolic element whose zero-frequency stiffness vanishes
or diverges will break this function. Resistors and capacitors contribute no
edges because they provide no static flux stiffness in the node flux basis.
Mutual inductances also contribute no edges because they couple branch fluxes
without providing a galvanic connection. With the coupled branches
represented by auxiliary branch current variables (see
[`calcAmnaind`](@ref)), a singular coupling matrix (a perfectly coupled
pair, `|k| = 1`) does not add flux null directions: the constitutive
equations constrain the branch fluxes rather than freeing them, so this
graph classification of the flux gauge freedom remains sound. Any
degeneracies which exist at perfect coupling live in the branch current
space (a branch current combination which the surrounding circuit leaves
physically undetermined) and produce a singular system caught at
factorization. Mutual coupling between inductors sharing a single branch
is rejected up front (see [`mnacoupledbranches`](@ref)).

The DC flux of each returned "floating" component is defined only up to a
common shift (a gauge degree of freedom). The modified nodal analysis
formulation adds one gauge fixing equation per floating component and
zero-frequency mode, see [`calcdcgaugeindices`](@ref). A net direct current
injected into a floating component is carried by the explicit average
voltage block when the component has a conductive path to somewhere that
can absorb it; only a component with no such path has no solution, which
`dcpinning` refuses (see `harmonics/directcurrent.jl`).

# Examples
```jldoctest
julia> JosephsonCircuits.calcstaticfluxcomponents([:P,:R,:C,:Lj,:C],[2 2 2 3 4;1 1 3 4 1],[1,50.0,1e-13,1e-9,1e-12],4)
2-element Vector{Vector{Int64}}:
 [2]
 [3, 4]

julia> JosephsonCircuits.calcstaticfluxcomponents([:P,:R,:L],[2 2 2;1 1 1],[1,50.0,1e-9],2)
Vector{Int64}[]

julia> JosephsonCircuits.calcstaticfluxcomponents([:P,:R,:L],[2 2 2;1 1 1],[1,50.0,Inf],2)
1-element Vector{Vector{Int64}}:
 [2]
```
"""
function calcstaticfluxcomponents(componenttypes::Vector{Symbol},
    nodeindices::Matrix{Int}, vvn::Vector, Nnodes::Int)


    # a union-find over the nodes (ground is node 1), with path halving
    parent = collect(1:Nnodes)
    function findroot(i::Int)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end

    # add an edge for every component which provides static flux stiffness
    for i in eachindex(componenttypes)
        if componenttypes[i] == :L || componenttypes[i] == :Lj
            # an infinite numeric inductance provides no static stiffness.
            # it is an open circuit at zero frequency. zero and NaN values
            # are rejected by checkstaticstiffnessvalues before this
            # function is called. symbolic values are assumed to provide
            # finite, nonzero static stiffness.
            v = vvn[i]
            if !checkissymbolic(v) && v isa Number && !isfinite(abs(v))
                continue
            end
            r1 = findroot(nodeindices[1, i])
            r2 = findroot(nodeindices[2, i])
            if r1 != r2
                parent[r1] = r2
            end
        end
    end

    # collect the components which do not contain the ground node
    groundroot = findroot(1)
    components = Dict{Int,Vector{Int}}()
    for p in 2:Nnodes
        r = findroot(p)
        if r != groundroot
            push!(get!(components, r, Int[]), p)
        end
    end

    floatingcomponents = collect(values(components))
    for component in floatingcomponents
        sort!(component)
    end
    sort!(floatingcomponents, by = first)

    return floatingcomponents
end

"""
    checkstaticstiffnessvalues(componenttypes::Vector{Symbol}, vvn::Vector)

Check that every linear inductor and Josephson junction has a finite, nonzero
numeric value (or a symbolic value), and throw an `ArgumentError` otherwise.
Symbolic values are accepted under the documented assumption that their
zero-frequency stiffness is finite and nonzero (the analysis will fail
otherwise), see [`calcstaticfluxcomponents`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.checkstaticstiffnessvalues([:P,:R,:L],[1,50.0,1e-9])

```
"""
function checkstaticstiffnessvalues(componenttypes::Vector{Symbol},
    vvn::Vector)
    for i in eachindex(componenttypes)
        if componenttypes[i] == :L || componenttypes[i] == :Lj
            v = vvn[i]
            if !checkissymbolic(v) && v isa Number
                if iszero(v)
                    throw(ArgumentError("A zero value for an inductor or "*
                        "Josephson junction is not supported: it would "*
                        "stamp an infinite inverse inductance. An ideal "*
                        "short is not a supported circuit element; connect "*
                        "the nodes directly or use a small finite "*
                        "inductance."))
                elseif isinf(abs(v))
                    throw(ArgumentError("An infinite value for an inductor "*
                        "or Josephson junction is not supported: it would "*
                        "make the scaling inductance Lscale infinite. An "*
                        "open circuit is represented by omitting the "*
                        "branch."))
                elseif isnan(abs(v))
                    throw(ArgumentError("An inductor or Josephson junction "*
                        "has a NaN value."))
                end
            end
        end
    end
    return nothing
end

"""
    isnumericallyzero(value, terms)

Return `true` if `value`, the floating point result of a linear combination of
the given `terms`, is approximately zero up to roundoff error.

# Examples
```jldoctest
julia> JosephsonCircuits.isnumericallyzero(2*pi*1.0 + (2*pi*(5e9-1) - 2*pi*5e9), (2*pi*1.0, -2*pi*5e9, 2*pi*(5e9-1)))
true

julia> JosephsonCircuits.isnumericallyzero(2*pi*1.0, (2*pi*1.0, 0.0, 0.0))
false
```
"""
function isnumericallyzero(value, terms)
    isfinite(value) || throw(ArgumentError("The combined mode frequency "*
        "must be finite."))
    all(isfinite, terms) || throw(ArgumentError("Every contributing "*
        "frequency term must be finite."))
    return isnumericallyzero(value, sum(abs, terms), length(terms))
end

# the core of the criterion from a precomputed magnitude sum and term
# count, so callers checking many combinations against the same terms (the
# per signal frequency, per mode validation of hblinsolve) can precompute
# the scale once instead of materializing a term vector per combination.
# this is the single definition of the tolerance.
function isnumericallyzero(value, scale::Real, nterms::Integer)
    n = max(1, 2*nterms - 1)
    epsT = eps(one(float(real(scale))))
    gamma = n*epsT/(1 - n*epsT)
    return abs(value) <= 4*gamma*scale
end

"""
    calcdcgaugeindices(floatingcomponents::Vector{Vector{Int}},
        wmodes::Vector, Nmodes::Int)

Return the indices of the node flux variables to which a gauge fixing
equation will be added. For each floating component of the static
flux-stiffness graph from [`calcstaticfluxcomponents`](@ref) and each mode
with zero frequency (DC), one index is returned, corresponding to the
lowest-numbered node of the component (the reference node). The DC flux of
a floating component only enters the equations through differences of node
fluxes within the component, so exactly one constraint per component and
zero-frequency mode removes the gauge degree of freedom without
overconstraining the system. The DC node flux reported as zero depends on the
node ordering (the `sorting` keyword of the solvers): the reference is the
lowest-numbered node of each component after sorting. All physical quantities
are gauge independent and unaffected by this choice.

# Examples
```jldoctest
julia> JosephsonCircuits.calcdcgaugeindices([[2],[3,4]],[0.0,2pi*4e9],2)
2-element Vector{Int64}:
 1
 3

julia> JosephsonCircuits.calcdcgaugeindices([[2],[3,4]],[2pi*4e9,2pi*8e9],2)
Int64[]
```
"""
function calcdcgaugeindices(floatingcomponents::Vector{Vector{Int}},
    wmodes::Vector, Nmodes::Int)

    gaugeindices = Int[]
    for component in floatingcomponents
        # the reference node of the component. nodes are "one indexed" so
        # subtract 2 to find the position in the matrices, which exclude the
        # ground node.
        p = first(component)
        for m in 1:Nmodes
            if iszero(wmodes[m])
                push!(gaugeindices, (p-2)*Nmodes + m)
            end
        end
    end
    sort!(gaugeindices)

    return gaugeindices
end

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
# a zero voltage Josephson junction requires zero average voltage across its
# branch, so `v` is constant on each connected component of the finite L/Lj
# graph, the components `calcstaticfluxcomponents` returns.
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
# The average voltages are carried as unknowns beside the periodic state,
# with these rows as their equations and the resistor current `G0 P v` as
# their coupling into the zero frequency nodal rows; a scattering block adds
# its own zero frequency relation between its port voltages and the port
# currents the solver already carries. That block of unknowns is small, its
# rows see no periodic state, and it is built, classified and solved apart
# from the periodic problem (see the sections below).
#
# It is classified whenever a zero frequency mode exists, and solved only
# when it has something to find. With no direct current injected every
# average voltage is zero, and so is every block current a nonsingular
# relation determines from them, so the periodic system as stamped -- a
# resistor open at zero frequency, a block with `i = 0` -- is exact and is
# solved alone. What the classification decides is whether that is so: a
# short or an ideal through in parallel with an inductive path leaves a
# current undetermined whether or not anything drives it, and is refused
# either way rather than given the zero it happens to start at.
#
# The conductance is read from the assembled `Gnm` rather than from the
# resistors, so that it keeps working when a port environment becomes a
# boundary stamp with no component behind it.

"""
    DCConductancePlan

The topology of the direct current voltage block.

Built from the circuit and its assembled conductance; holds nothing which
depends on the sources. It exists whenever there is a zero frequency mode,
with no components when every node is held at zero by an inductive path to
ground, so that the block relations and the classification have the nodal
rows and the conductance to work from.

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

A solution exists exactly when the analysis has a zero frequency mode. When
no direct current is injected it is the zero every average voltage sits at,
which is an answer and not the absence of one: a node shorted to ground
sits at zero volts. Only an analysis with no zero frequency mode has no
average voltage to report.

# Fields
- `nodevoltage`: volts, indexed by node with ground first and identically
    zero.
- `scaledcurrent`: `G0 P v`, the direct current each node carries in the
    solver's scaled units, with ground dropped; what the Kirchhoff current
    law validation and [`applydcconductance`](@ref) read.
"""
struct DCConductanceSolution{T}
    nodevoltage::Vector{T}
    scaledcurrent::Vector{T}
end

"""
    dcconductanceplan(floatingcomponents, Gnm, wmodes, Nmodes, Nnodes)

Build the [`DCConductancePlan`](@ref), or `nothing` when there is no zero
frequency mode and so no average voltage at all.

Neither an empty conductance nor an empty set of floating components is a
reason to build none. A circuit whose only direct current devices are
scattering blocks has an empty `G0` and still needs the block rows; a
circuit whose every node is held at zero by an inductive path to ground
has no voltage to solve for and still has block currents to classify, and
gating on either left such a block with the artificial `i = 0` row the
stamp writes and nothing looking at whether that row is right.

Throws an `ArgumentError` when the zero frequency conductance has a
non-finite or non-real entry, which has no direct current meaning.
"""
function dcconductanceplan(floatingcomponents::Vector{Vector{Int}},
        Gnm::SparseMatrixCSC, wmodes::AbstractVector, Nmodes::Integer,
        Nnodes::Integer)

    m0 = findfirst(iszero, wmodes)
    isnothing(m0) && return nothing

    n = Nnodes - 1
    dcrows = collect(Int(m0):Nmodes:n*Nmodes)
    G0c = Gnm[dcrows, dcrows]
    # A conductance carries direct current only if it is real: a complex one
    # at zero frequency has no steady state meaning, and would come from a
    # frequency dependent law whose limit at DC is not a conductance, so it
    # is refused rather than partly used. The test is against a tolerance
    # rather than an exact zero, because an expression whose imaginary part
    # cancels analytically leaves roundoff when evaluated; the scale is the
    # largest conductance in the matrix, so the decision does not depend on
    # the circuit's units.
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
        ones(eltype(G0), count(p -> p > 1,
            reduce(vcat, floatingcomponents; init = Int[]))),
        n, nc)


    return DCConductancePlan(Int(m0), dcrows, [Int.(c) for c in
        floatingcomponents], componentof, lift, G0, Y)
end

"""
    applydcconductance(bnm, plan, sol, Nmodes)

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
# The equation as explicit rows.
#
# Eliminating the average voltages before Newton -- solving `Y v = j` once
# and folding the resistor current `G0 P v` into the zero frequency source
# -- is exact while the direct current devices are linear conductances and
# the sources are prescribed, and it is how this began. It cannot go
# further. A scattering block's direct current relation is a pencil between
# its port voltages and its port currents, so the current is a genuine
# unknown and there is nothing to eliminate; a short or an ideal through
# has a free current direction and no determined value at all. Reaching
# those needs `v` carried as an unknown with its equation as a row, which
# is what this builds, and it is the only path now.
#
# The row is the component sum. Adding the zero frequency Kirchhoff
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
# substitution is exactly the elimination, which is the sense in which the
# explicit rows and a hand elimination must agree and the reason the tests
# can demand it.
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
- `plan`: the topology.
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

"""
    nvoltages(t::TransportRows)

The number of explicit average voltages, one per floating static flux
component.
"""
nvoltages(t::TransportRows) = length(t.j)

"""
    transportrows(plan::DCConductancePlan, bnm, Nmodes)

Build the [`TransportRows`](@ref) for a source `bnm`.

`bnm` must be the applied source, not one corrected for the resistor
current: that current appears here as the coupling term, and taking it
from a corrected source would count it twice.
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
is the quantity a hand elimination would subtract from the source.
"""
function transportcurrent!(d::AbstractVector, t::TransportRows,
        v::AbstractVector)
    mul!(d, t.coupling, v)
    return d
end

"""
    dcsolutionfrom(plan::DCConductancePlan, v::AbstractVector)

Package the solved component voltages as a
[`DCConductanceSolution`](@ref), the one form in which the direct current
answer is reported.
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

When none is, every average voltage is zero, and so is every block port
current a nonsingular relation determines from them, so the explicit block
has nothing to find: the `i = 0` rows the scattering stamp writes are then
the right answer rather than a simplification, and the block is classified
but not carried. Whether the relations are nonsingular in a way the circuit
can see is what the classification decides, and it does not depend on this.

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

# =====================================================================
# Scattering blocks at zero frequency: the constitutive pencil of each
# block replaces its `i = 0` row when the circuit injects direct current.

# A scattering block at zero frequency.
#
# `evaluatehybrid!` writes `B = 0`, `C = I` at zero frequency, the equation
# `i = 0`: in the harmonic system alone every block is an open circuit at
# direct current, since the periodic state carries no average voltage for
# a block to respond to. When the circuit injects direct current the
# explicit block (above) carries the average port voltages,
# and each scattering block's zero frequency row is replaced by the same
# pencil it satisfies at every other frequency,
#
#     B(0) V - C(0) i = 0,   B(0) = R^(-1/2)(I - S(0)),
#                            C(0) = R^(1/2)(I + S(0)),
#
# with `V` the average port voltage and `i` the average port current, the
# auxiliary unknown the `i = 0` row would otherwise pin to zero. This adds
# no coordinate; it removes a constraint.
#
# Nothing is inverted, which is the point of the pencil form. A resistive
# block has `I + S(0)` invertible and a determined current. An ideal short
# has `S(0) = -1`, so `C(0) = 0` and the row reads `B(0) V = 0`: the
# voltage is constrained and the current is free, which is why the current
# has to be an unknown rather than a value to be computed.

"""
    DCBlockDescriptor

The zero frequency constitutive pencil of one scattering block.

# Fields
- `B0`, `C0`: `R^(-1/2)(I - S(0))` and `R^(1/2)(I + S(0))`, so the block's
  rows read `B0 V - C0 i = 0`.
- `signalnodes`, `refnodes`: the terminals of each port.
- `auxbase`: the auxiliary index base, as in [`StampedScatteringBlock`](@ref).
- `freecurrents`: the dimension of the null space of `C0`, the number of
  port current directions the block leaves undetermined. Zero for a block
  whose current is fixed by its voltages.
- `name`: the block's first port, for messages.
"""
struct DCBlockDescriptor
    B0::Matrix{Float64}
    C0::Matrix{Float64}
    signalnodes::Vector{Int}
    refnodes::Vector{Int}
    auxbase::Int
    freecurrents::Int
    name::String
end

nports(d::DCBlockDescriptor) = length(d.signalnodes)

# a block without a `dcmodel` field evaluates its own data at zero; every
# `ScatteringParameters` has the field, so this is a defensive fallback
dcmodelof(blk) = hasproperty(blk, :dcmodel) ? blk.dcmodel : ScatteringLimit()

"""
    dclimit(sb::StampedScatteringBlock, n, atol)

The zero frequency scattering matrix read from a block's own data, or
`nothing` and the reason the data has none: `:range` when it is tabulated
and does not reach zero, `:nonfinite` when the value there is unbounded,
`:complex` when it has no real limit.

Read and not decided: what to do about a block with no zero frequency data
depends on whether direct current is asked of it, which is the caller's
question. See [`dcblockdescriptor`](@ref).
"""
function dclimit(sb::StampedScatteringBlock, n::Integer, atol::Real)
    blk = sb.block
    pr = blk.provider
    if pr isa TabulatedMatrixProvider && pr.extrapolation == :error &&
            !(pr.frequencies[1] <= 0.0 <= pr.frequencies[end])
        return nothing, :range
    end
    S = Array{Complex{Float64},3}(undef, n, n, 1)
    evaluatescattering!(S, blk, [0.0])
    S0 = @view S[:,:,1]
    all(isfinite, S0) || return nothing, :nonfinite
    m = maximum(abs∘imag, S0)
    m <= atol*max(1, maximum(abs, S0)) || return nothing, :complex
    return Matrix{Float64}(real.(S0)), :ok
end

"""
    dcblockdescriptor(sb::StampedScatteringBlock; atol = 1e-10,
        required = true)

Evaluate the zero frequency pencil of a stamped block.

The matrix comes from the block's [`AbstractDCModel`](@ref): by default from
its own scattering data evaluated at zero, and from the stated model when
the block declared one.

Evaluated or stated, `S(0)` must be real and finite. A complex zero
frequency scattering matrix describes a block with no direct current limit,
which is refused here rather than resolved by a convention, on the same
principle as a complex direct current conductance.

A block which states no model and whose data has no zero frequency value is
refused when `required`, which is the case when direct current is injected
and the block would have to carry its share. When none is, the block is
open at zero frequency, which is the `i = 0` row the stamp already writes,
and `nothing` is returned so the caller leaves that row alone: a block
which cannot say what it does at direct current, and is asked to carry
none, carries none. Tabulated data which starts above zero is the common
case, and a circuit driven at its pump alone is not refused for it.
"""
function dcblockdescriptor(sb::StampedScatteringBlock; atol::Real = 1e-10,
        required::Bool = true)
    blk = sb.block
    n = blk.nports
    model = dcmodelof(blk)
    S0 = if model isa ScatteringLimit
        S0, why = dclimit(sb, n, atol)
        if isnothing(S0)
            required || return nothing
            if why === :range
                throw(ArgumentError(lazy"the scattering block $(sb.name) is tabulated over a frequency range which does not reach zero, so its direct current behavior cannot be read from it, and this circuit injects direct current. State the zero frequency limit with the dcmodel keyword: OpenDC(), ShortDC(), ThroughDC() or ScatteringDC(S0), or pass extrapolation = :constant to hold the lowest tabulated value down to zero."))
            elseif why === :nonfinite
                throw(ArgumentError(lazy"the scattering block $(sb.name) returned a non-finite S(0), so its direct current behavior cannot be read from it. A block whose limit exists but is not evaluable at zero -- a series capacitance written as 1/(im*w*C), whose limit is the open circuit -- has to state that limit with the dcmodel keyword: OpenDC(), ShortDC(), ThroughDC() or ScatteringDC(S0)."))
            else
                throw(ArgumentError(lazy"the scattering block $(sb.name) has a complex S(0); a block with no real zero frequency limit has no direct current behavior to stamp. State the limit with the dcmodel keyword if the block has one."))
            end
        end
        S0
    else
        # checked for size and passivity when the block was constructed
        dcscatteringmatrix(model, n)
    end

    r2 = sqrt.(float.(blk.zref))
    B0 = Matrix{Float64}(undef, n, n)
    C0 = Matrix{Float64}(undef, n, n)
    @inbounds for q in 1:n, p in 1:n
        d = p == q ? 1.0 : 0.0
        B0[p,q] = (d - S0[p,q]) / r2[q]
        C0[p,q] = (d + S0[p,q]) * r2[q]
    end
    # the current directions the block does not determine: an ideal short
    # has C0 = 0 and leaves all of them free
    free = n - rank(C0; atol = atol*max(1, maximum(abs, C0)))
    return DCBlockDescriptor(B0, C0, sb.signalnodes, sb.refnodes, sb.auxbase,
        free, sb.name)
end

"""
    DCBlockRows

The zero frequency rows of every scattering block in a circuit, ready to be
applied to a canonical state.

# Fields
- `descriptors`: one [`DCBlockDescriptor`](@ref) per block.
- `currentindex`: for each block, the window position of each port's zero
  frequency current: its slot among the zero frequency entries of the state,
  one per node and then one per auxiliary unknown.
- `signalcomponent`, `refcomponent`: for each block, the static flux
  component of each port terminal, or zero where the average voltage is held
  at zero by a path to ground through inductance.
- `scale`: the solver scale (see [`calcsolverscale`](@ref)), which carries
  the average voltage into the units the stamp's rows are written in.
- `transportterms`: `(component, window position of the current, sign)` for
  each block current which crosses a static flux component boundary.

The units are the stamp's throughout. At a nonzero frequency the block's row
is `B (im w scale phi) - C i`, so the stamp's voltage is `scale/phi0` times
the physical one; the explicit direct current coordinate is `v = V/phi0`, so
the same voltage is `scale * v` and the row is `B0 (scale dv) - C0 i` with
the current the solver already carries.
"""
struct DCBlockRows{T}
    descriptors::Vector{DCBlockDescriptor}
    currentindex::Vector{Vector{Int}}
    signalcomponent::Vector{Vector{Int}}
    refcomponent::Vector{Vector{Int}}
    scale::T
    # (component, window position of the current, sign): the block currents
    # which cross a static flux component's boundary and so survive its
    # transport row's sum. A port with both terminals inside one component cancels,
    # exactly as an inductor branch does.
    transportterms::Vector{Tuple{Int,Int,Int}}
end

Base.isempty(r::DCBlockRows) = isempty(r.descriptors)

"""
    freecurrents(r::DCBlockRows)

The total number of port current directions the blocks leave undetermined.
Nonzero means some current is fixed by node level Kirchhoff rather than by
the block, which is the case a short or an ideal through presents.
"""
freecurrents(r::DCBlockRows) = sum(d -> d.freecurrents, r.descriptors;
    init = 0)

"""
    dcblockrows(blocks, componentof, Nmodes, modeindex, nnodaldc, scale;
        required = true)

Build the [`DCBlockRows`](@ref) for the stamped blocks of a circuit.

`nnodaldc` is the point in the zero frequency block where the nodal entries
end and the auxiliary ones begin. A block with no zero frequency data is
refused when `required` and left open otherwise; see
[`dcblockdescriptor`](@ref).
"""
function dcblockrows(blocks::AbstractVector, componentof::Vector{Int},
        Nmodes::Integer, modeindex::Integer, nnodaldc::Integer, scale;
        required::Bool = true)
    descriptors = DCBlockDescriptor[]
    currentindex = Vector{Int}[]
    signalcomponent = Vector{Int}[]
    refcomponent = Vector{Int}[]
    for sb in blocks
        d = dcblockdescriptor(sb; required = required)
        isnothing(d) && continue
        n = nports(d)
        idx = Vector{Int}(undef, n)
        for p in 1:n
            # the complex state index of this port's zero frequency current,
            # and the slot it belongs to; the window holds one entry per slot
            # in the same order
            s = sb.auxbase + (p-1)*Nmodes + modeindex
            idx[p] = (s - 1) ÷ Nmodes + 1
        end
        push!(descriptors, d)
        push!(currentindex, idx)
        push!(signalcomponent, [componentof[k] for k in d.signalnodes])
        push!(refcomponent, [componentof[k] for k in d.refnodes])
    end
    # The transport row of a component is the sum of its nodes' zero
    # frequency Kirchhoff equations, in which a port current enters at its
    # signal node and leaves at its reference node. Both inside one
    # component and it cancels; otherwise it is a current the component
    # exchanges with the rest of the circuit, and the row has to carry it or
    # the block is invisible to the average voltages.
    terms = Tuple{Int,Int,Int}[]
    for b in eachindex(descriptors)
        for p in eachindex(currentindex[b])
            sc, rc = signalcomponent[b][p], refcomponent[b][p]
            sc == rc && continue
            iszero(sc) || push!(terms, (sc, currentindex[b][p], 1))
            iszero(rc) || push!(terms, (rc, currentindex[b][p], -1))
        end
    end
    return DCBlockRows(descriptors, currentindex, signalcomponent,
        refcomponent, scale, terms)
end

"""
    addblocktransport!(Fv, r::DCBlockRows, u)

Add the block currents which cross a component boundary to that component's
transport row, reading each current from `u` at the position
`r.currentindex` gives it.

Every one of them, unconditionally. A reference row is chosen after this,
from the assembled descriptor, so there is no row here which is known in
advance to be redundant and no current which may be dropped on the grounds
that it lands in one.
"""
function addblocktransport!(Fv::AbstractVector, r::DCBlockRows,
        u::AbstractVector)
    @inbounds for (c, idx, sgn) in r.transportterms
        Fv[c] += sgn * u[idx]
    end
    return Fv
end

# the average voltage of a node in the explicit coordinate, zero where a
# path to ground through inductance holds it there
@inline _vof(v, c) = iszero(c) ? zero(eltype(v)) : @inbounds v[c]

"""
    addblockdc!(Fc, r::DCBlockRows, u, v)

Replace each block's zero frequency row in `Fc`, which holds `-i`, by the
block's own relation `B0 (scale dv) - C0 i`, with the block currents read
from `u` at the positions `r.currentindex` gives them and the
average voltages from `v`.

The correction is added rather than written, because the row also carries
the Kirchhoff coupling of the current into the node equations, which is
unchanged and must survive.
"""
function addblockdc!(Fc::AbstractVector, r::DCBlockRows, u::AbstractVector,
        v::AbstractVector)
    for (b, d) in enumerate(r.descriptors)
        idx = r.currentindex[b]
        n = length(idx)
        sc, rc = r.signalcomponent[b], r.refcomponent[b]
        @inbounds for p in 1:n
            # B0 (scale dv) - C0 i, less the -i already in the row
            acc = u[idx[p]]
            for q in 1:n
                acc -= d.C0[p,q] * u[idx[q]]
                acc += d.B0[p,q] * r.scale * (_vof(v, sc[q]) - _vof(v, rc[q]))
            end
            Fc[idx[p]] += acc
        end
    end
    return Fc
end

# =====================================================================
# The block's contribution to a canonical residual or product: the
# transport rows and the coupling, on the window of the state.

# The explicit direct current block's contribution: the resistor current the
# average voltages drive into the zero frequency nodal rows, and the
# transport rows themselves.
#
# The sign is the elimination's, read backwards. The residual is
# `A x - bnm`, and eliminating replaces `bnm` by `bnm - G0 P v`, so carrying
# `v` instead means adding `G0 P v` here and leaving the applied source
# alone. Doing both would count the resistor current twice.
#
# Everything the block reads and writes is the window: the zero frequency
# entries, scattered through the internal state, and the explicit voltages
# after it. The window is gathered by index into a buffer of its own, worked
# on with local indices, and scattered back, which is what lets the
# arithmetic run on a device without scalar indexing; it is a few hundred
# numbers where the state is tens of thousands.
function addtransport!(Fc::AbstractVector, work::CanonicalWork,
        u::AbstractVector; residual::Bool = true)
    isnothing(work.transport) && return Fc
    Fw, uw = work.Fwindow, work.uwindow
    _gatherperm!(Fw, Fc, work.window)
    _gatherperm!(uw, u, work.window)
    if isnothing(work.update)
        addtransportwindow!(Fw, uw, work; residual)
    else
        # in place, where the state is: three array operations and no copy
        applydcupdate!(Fw, uw, work.update; residual)
    end
    _scatterperm!(Fc, Fw, work.window)
    return Fc
end

# `Fw` and `uw` are the window; every index below is local to it.
function addtransportwindow!(Fw::AbstractVector, uw::AbstractVector,
        work::CanonicalWork; residual::Bool = true)
    t = work.transport
    L = work.layout
    vr = (L.ndc + 1):(L.ndc + L.nvdc)
    v = view(uw, vr)
    d = work.dwork
    transportcurrent!(d, t, v)
    @inbounds for k in eachindex(d)
        Fw[k] += d[k]
    end
    if residual
        transportresidual!(view(Fw, vr), t, v)
    else
        mul!(view(Fw, vr), t.Y, v)
    end
    # the blocks' own zero frequency rows, which replace `i = 0`, and the
    # currents they exchange across a component boundary
    br = work.blockrows
    if !isnothing(br)
        addblockdc!(Fw, br, uw, v)
        addblocktransport!(view(Fw, vr), br, uw)
    end
    # the rows spent on the directions nothing determines, written over what
    # was there: the equation they replaced was the redundant one
    pn = work.pinning
    if !isnothing(pn)
        idx = work.dclocal
        @inbounds for j in eachindex(pn.rows)
            Fw[idx[pn.rows[j]]] = uw[idx[pn.cols[j]]]
        end
    end
    return Fw
end

# =====================================================================
# Evaluating in canonical coordinates.
#
# The residual, the Jacobian vector product and the preconditioner are all
# written against the internal layout, and the canonical vector holds that
# layout in its first block. So each is handed a view of that block and
# nothing is copied; what is added is the direct current block, on the
# window.

# =====================================================================
# The direct current solve, where the state is.
#
# The subsystem is constant, so its exact solve is a fixed linear map, but
# that map cannot be stored as a matrix: the subsystem is ill conditioned
# (a circuit with a ten million to one impedance ratio gives it a condition
# number of that order) and its unknowns are large in the solver's scaled
# units, so an explicit inverse leaves an absolute error the solve cannot
# drive to zero. It stays a linear solve, and what moves to the device is
# the factorization rather than the answer: the same factors, permutation
# and substitution order as the host path, in a kernel with one work item,
# so the two agree exactly and the window does not cross the bus per
# application.
#
# One work item because the substitutions are sequential, so the kernel's
# cost grows with the subsystem while the host copies it replaces cost a
# fixed latency. It wins for a handful of unknowns and loses badly beyond
# that. A handful is the common case, since there is one unknown per
# floating static flux component and one per scattering block port current;
# a larger subsystem keeps the host path, which is correct and merely
# copies.

# the largest subsystem the sequential device solve is launched for; the
# measured crossover against the host copies
const DCDEVICESOLVEMAX = 8

"""
    DCFactorization

The factorization of the direct current subsystem and the indices it acts
on, resident on a backend.

# Fields
- `factors`: the packed unit lower and upper triangle `lu` produced.
- `perm`: its row permutation.
- `index`: the canonical position of each subsystem coordinate.
- `work`: a device scratch vector of length `n`, the permuted right hand
    side and then the solution.
- `n`: the subsystem size.
"""
struct DCFactorization{V,I}
    index::I
    perm::I
    factors::V
    work::V
    n::Int
end

function DCFactorization(F::LinearAlgebra.LU, index::Vector{Int}, backend)
    n = size(F.factors, 1)
    return DCFactorization(tobackend(backend, index),
        tobackend(backend, Vector{Int}(F.p)),
        tobackend(backend, Vector{Float64}(vec(F.factors))),
        tobackend(backend, zeros(Float64, n)), n)
end

# One work item: the substitutions are sequential and the system is a
# handful of unknowns, so what is being avoided is the bus and not the
# arithmetic.
@kernel function dcsolvekernel!(z, @Const(r), @Const(index),
        @Const(perm), @Const(factors), b, n)
    @index(Global)
    @inbounds begin
        for k in 1:n
            b[k] = r[index[perm[k]]]
        end
        for k in 1:n, i in (k+1):n            # L y = P b, unit lower
            b[i] -= factors[i + (k-1)*n]*b[k]
        end
        for k in n:-1:1                       # U x = y
            for i in (k+1):n
                b[k] -= factors[k + (i-1)*n]*b[i]
            end
            b[k] /= factors[k + (k-1)*n]
        end
        for k in 1:n
            z[index[k]] = b[k]
        end
    end
end

"""
    applydcsolve!(z, r, d::DCFactorization)

Solve the direct current subsystem for the coordinates of `r` it owns and
write the answer into the same coordinates of `z`, in place and on the
backend the factorization lives on.
"""
function applydcsolve!(z::AbstractVector, r::AbstractVector,
        d::DCFactorization)
    backend = KernelAbstractions.get_backend(z)
    kernel! = dcsolvekernel!(backend)
    kernel!(z, r, d.index, d.perm, d.factors, d.work, d.n; ndrange = 1)
    KernelAbstractions.synchronize(backend)
    return z
end

"""
    dcsubsystem(work::CanonicalWork)

The dense direct current block `[v; i]` of the canonical Jacobian: the
transport rows and the blocks' zero frequency rows, in the order
[`dcsubsystemindices`](@ref) gives.
"""
function dcsubsystem(work::CanonicalWork)
    t, br = work.transport, work.blockrows
    nc = nvoltages(t)
    idx = dcsubsystemindices(work)
    nb = length(idx) - nc
    A = zeros(Float64, nc + nb, nc + nb)
    A[1:nc, 1:nc] .= t.Y
    # the reference rows are written last and unconditionally: the transport
    # rows no longer carry one, so a circuit with no blocks still needs them.
    # A block current is named by its window position, as the block rows
    # name it.
    slots = dcsubsystemlocal(work)
    local_ = Dict(slots[nc+k] => nc + k for k in 1:nb)
    if !isnothing(br)
        # the block currents each component exchanges across its boundary
        for (c, ci, sgn) in br.transportterms
            A[c, local_[ci]] += sgn
        end
        # and the blocks' own rows: B0 (scale dv) - C0 i
        for (b, d) in enumerate(br.descriptors)
            ci = br.currentindex[b]
            sc, rc = br.signalcomponent[b], br.refcomponent[b]
            for p in eachindex(ci)
                row = local_[ci[p]]
                for q in eachindex(ci)
                    A[row, local_[ci[q]]] -= d.C0[p,q]
                    iszero(sc[q]) || (A[row, sc[q]] += d.B0[p,q]*br.scale)
                    iszero(rc[q]) || (A[row, rc[q]] -= d.B0[p,q]*br.scale)
                end
            end
        end
    end
    pn = work.pinning
    if !isnothing(pn)
        for j in eachindex(pn.rows)
            A[pn.rows[j], :] .= 0.0
            A[pn.rows[j], pn.cols[j]] = 1.0
        end
    end
    return A
end

"""
    dcsubsystemindices(work::CanonicalWork)

The canonical positions the direct current subsystem occupies: the explicit
voltage block, then every block port's zero frequency current.
"""
function dcsubsystemindices(work::CanonicalWork)
    L = work.layout
    return [windowindex(L, k) for k in dcsubsystemlocal(work)]
end

"""
    dcsubsystemlocal(work::CanonicalWork)

The same positions, local to the window.
"""
function dcsubsystemlocal(work::CanonicalWork)
    L, br = work.layout, work.blockrows
    idx = collect(L.ndc .+ (1:L.nvdc))
    isnothing(br) || for v in br.currentindex, i in v
        push!(idx, i)
    end
    return idx
end

# =====================================================================
# Classifying the direct current subsystem: its constant side, its
# equilibration, the coupling into the nodal rows, and the references a
# singular subsystem needs.

"""
    dcsubsystemrhs(work::CanonicalWork)

The constant side of the direct current subsystem: the injected current on
the transport rows, zero on the blocks' own rows.
"""
function dcsubsystemrhs(work::CanonicalWork)
    t = work.transport
    n = length(dcsubsystemindices(work))
    b = zeros(Float64, n)
    b[1:nvoltages(t)] .= t.j
    return b
end

# The direct current subsystem mixes volts and amperes, and its rows are
# Kirchhoff sums in one place and constitutive relations in another, so its
# entries carry the circuit's impedance scale. A rank decision on the raw
# matrix would then depend on that scale: the same circuit written at a
# different impedance could be called singular or not. Scaling the rows
# and the columns to unit infinity norm first (two passes) removes the
# units from a question whose answer is a structural fact about the
# circuit.
function equilibrate(A::AbstractMatrix)
    dr = ones(Float64, size(A, 1))
    dc = ones(Float64, size(A, 2))
    B = copy(A)
    for _ in 1:2
        for i in axes(B, 1)
            m = maximum(abs, view(B, i, :); init = 0.0)
            iszero(m) && continue
            dr[i] /= m
            view(B, i, :) ./= m
        end
        for j in axes(B, 2)
            m = maximum(abs, view(B, :, j); init = 0.0)
            iszero(m) && continue
            dc[j] /= m
            view(B, :, j) ./= m
        end
    end
    return B, dr, dc
end

"""
    dccoupling(work::CanonicalWork)

`H`: the zero frequency nodal current each direct current unknown drives,
with the unknowns in the order [`dcsubsystemindices`](@ref) gives.

This is everything outside the direct current subsystem which can see those
unknowns. The average voltages reach the rest of the residual only as the
resistor current `G0 P v` added to the zero frequency nodal rows, and a
block port current only as the `+1` and `-1` it contributes at its two
terminals. Nothing else in the harmonic residual reads them, so a direction
with `H N = 0` is invisible to the whole problem and not only to the
subsystem, which is the condition a reference has to meet.

The transport rows are the component sums of these same nodal rows, which is
why `H` and the subsystem cannot disagree about a sign: `Y = P' (G0 P)` and
a block current enters its signal component's row with the sign it enters
its signal node's row.
"""
function dccoupling(work::CanonicalWork)
    t, br = work.transport, work.blockrows
    nc = nvoltages(t)
    n = size(t.coupling, 1)
    nb = length(dcsubsystemindices(work)) - nc
    I, J, V = Int[], Int[], Float64[]
    C = t.coupling
    for j in axes(C, 2), k in nzrange(C, j)
        push!(I, C.rowval[k]); push!(J, j); push!(V, C.nzval[k])
    end
    if !isnothing(br)
        c = nc
        for (b, d) in enumerate(br.descriptors)
            for p in eachindex(br.currentindex[b])
                c += 1
                # node 1 is ground and has no row
                d.signalnodes[p] > 1 &&
                    (push!(I, d.signalnodes[p]-1); push!(J, c); push!(V, 1.0))
                d.refnodes[p] > 1 &&
                    (push!(I, d.refnodes[p]-1); push!(J, c); push!(V, -1.0))
            end
        end
    end
    return sparse(I, J, V, n, nc + nb)
end

# the name of the block each block current coordinate belongs to, for
# messages; the voltage coordinates come first and have none
function dccoordinatenames(work::CanonicalWork)
    t, br = work.transport, work.blockrows
    names = fill("", nvoltages(t))
    isnothing(br) && return names
    for (b, d) in enumerate(br.descriptors)
        for _ in eachindex(br.currentindex[b])
            push!(names, d.name)
        end
    end
    return names
end

"""
    dcpinning(work::CanonicalWork)

Return the [`DCPinning`](@ref) a singular direct current subsystem needs,
`nothing` when it is nonsingular, or throw when it has no solution or an
undetermined direction the rest of the circuit can see.

The subsystem this reads, `dcsubsystem(work)`, must be the complete
unreferenced descriptor: every transport row and every block relation,
with no reference chosen, which is why this is called once at
[`CanonicalWork`](@ref) construction, before any reference exists. Choosing
one earlier, from the resistors alone, can discard a row a block has made
necessary, and no later check can recover it.
"""
function dcpinning(work::CanonicalWork)
    isnothing(work.transport) && return nothing
    A = dcsubsystem(work)
    isempty(A) && return nothing
    B, dr, dc = equilibrate(A)
    F = svd(B)
    tol = maximum(F.S; init = 0.0) * maximum(size(B)) * eps()
    k = count(<=(tol), F.S)
    k == 0 && return nothing

    # the left null space: directions in which the equations say nothing, so
    # a constant side with a component along one of them cannot be met
    Y = F.U[:, end-k+1:end]
    b = dcsubsystemrhs(work) .* dr
    if norm(Y'b) > max(tol, eps()) * max(1.0, norm(b))
        throw(ArgumentError(lazy"No direct current solution exists: direct current is injected into a subnetwork which has no path carrying it away. The zero frequency mode is the average voltage, so a subnetwork whose average voltage is unconstrained cannot absorb a net current; give it a path to ground, or drive it differentially."))
    end

    # Which of the undetermined directions are gauges: a direction is one
    # only if the rest of the residual cannot see it, `H N = 0`. Anything
    # else is a physical quantity the circuit leaves undetermined, which is
    # refused rather than pinned to one of infinitely many answers.
    Nhat = F.V[:, end-k+1:end]
    H = dccoupling(work) * Diagonal(dc)
    G = Matrix(H * Nhat)
    hs = maximum(abs, H; init = 0.0)
    gtol = max(hs, 1.0) * maximum(size(G)) * sqrt(eps())
    if any(>(gtol), svdvals(G))
        names = dccoordinatenames(work)
        w = vec(maximum(abs, Nhat; dims = 2))
        seen = String[]
        for c in eachindex(names)
            isempty(names[c]) && continue
            w[c] <= sqrt(eps()) && continue
            names[c] in seen || push!(seen, names[c])
        end
        # the joined list is built outside the message: a comma inside a
        # `lazy` interpolation ends the interpolated expression
        involved = join(seen, ", ")
        throw(ArgumentError(lazy"The direct current network leaves a branch current undetermined, and that current is visible to the rest of the circuit: changing it moves the current at a node, and with it the static flux of any inductor or junction in parallel. The blocks whose zero frequency currents are involved are $(involved). An ideal short or through in parallel with an inductive branch has no unique direct current solution; give the block a finite series impedance at zero frequency, or give the parallel branch one, so that the division is determined."))
    end

    # Which equations to give up, and which coordinate each one fixes.
    # Pivoting on the left null space picks rows which carry the redundancy,
    # so what is left still spans the row space; pivoting on the null space
    # picks coordinates the directions move, so the references are
    # independent. Both are done in the equilibrated coordinates.
    rows = sort!(qr(Y', ColumnNorm()).p[1:k])
    cols = sort!(qr(Nhat', ColumnNorm()).p[1:k])
    pn = DCPinning(rows, cols)

    # the references have to leave a nonsingular system, which the two
    # pivoted choices give but do not guarantee jointly
    Ap = copy(A)
    for j in 1:k
        Ap[rows[j], :] .= 0.0
        Ap[rows[j], cols[j]] = 1.0
    end
    Bp, _, _ = equilibrate(Ap)
    if minimum(svdvals(Bp)) <= maximum(size(Bp)) * eps()
        error("the direct current references left a singular subsystem, which they must not: this is a bug in `dcpinning`, not a property of the circuit.")
    end
    return pn
end

# =====================================================================
# The direct current block as a matrix.
#
# Every term of `addtransportwindow!` is linear in the window, and some rows
# are added to while others are overwritten, so the whole update is
#
#     Fw <- keep .* Fw + M uw + c
#
# with `keep` zero on the overwritten rows. In that form it is three array
# operations rather than a walk over scattered indices, so it can run where
# the state lives instead of being copied to the host and back.
#
# `M` and `c` are read off the scalar implementation by probing it one
# basis vector at a time rather than assembled a second time by hand, which
# is O(window) work once and makes the two forms agree by construction.

"""
    DCUpdate

The direct current block's contribution to the residual as `keep .* Fw +
M*uw + c`, in whatever array type the state uses.

`cresidual` carries the injected current; the Jacobian vector product uses
the same `M` and `keep` with no constant, the two differing only by that.
"""
struct DCUpdate{V,I}
    keep::V
    rowptr::I          # the matrix by rows, so one work item owns one row
    colval::I
    nzval::V
    cresidual::V
end

# One work item per row of the window: its own entry, kept or not, plus its
# row of the matrix, plus the constant when this is a residual.
@kernel function dcupdatekernel!(Fw, @Const(keep), @Const(rowptr),
        @Const(colval), @Const(nzval), @Const(uw), @Const(c), alpha)
    i = @index(Global)
    @inbounds begin
        acc = keep[i]*Fw[i] + alpha*c[i]
        for k in rowptr[i]:(rowptr[i+1] - 1)
            acc += nzval[k]*uw[colval[k]]
        end
        Fw[i] = acc
    end
end

"""
    dcupdate(work::CanonicalWork)

Build the [`DCUpdate`](@ref) by probing `addtransportwindow!`, so the matrix
form and the scalar form agree by construction.
"""
function dcupdate(work::CanonicalWork)
    L = work.layout
    nw = L.ndc + L.nvdc
    nw == 0 && return nothing
    Fw = zeros(Float64, nw); uw = zeros(Float64, nw)

    # `keep` is one where the row is added to and zero where it is written
    fill!(Fw, 1.0); fill!(uw, 0.0)
    addtransportwindow!(Fw, uw, work; residual = false)
    keep = copy(Fw)

    # the constant, from a zero point
    fill!(Fw, 0.0)
    addtransportwindow!(Fw, uw, work; residual = true)
    c = copy(Fw)

    # and the linear part, column by column, against the product form which
    # carries no constant
    I, J, V = Int[], Int[], Float64[]
    for k in 1:nw
        fill!(Fw, 0.0); fill!(uw, 0.0); uw[k] = 1.0
        addtransportwindow!(Fw, uw, work; residual = false)
        for i in 1:nw
            iszero(Fw[i]) && continue
            push!(I, i); push!(J, k); push!(V, Fw[i])
        end
    end
    # by rows: the transpose of a column major sparse matrix is the row
    # major form of the original, which is what the kernel walks
    Mt = sparse(J, I, V, nw, nw)
    return DCUpdate(keep, Mt.colptr, Mt.rowval, Mt.nzval, c)
end

"""
    applydcupdate!(Fw, uw, up::DCUpdate; residual = true)

Apply the direct current block in its matrix form.
"""
function applydcupdate!(Fw::AbstractVector, uw::AbstractVector,
        up::DCUpdate; residual::Bool = true)
    backend = KernelAbstractions.get_backend(Fw)
    kernel! = dcupdatekernel!(backend)
    kernel!(Fw, up.keep, up.rowptr, up.colval, up.nzval, uw, up.cresidual,
        residual ? one(eltype(Fw)) : zero(eltype(Fw));
        ndrange = length(Fw))
    KernelAbstractions.synchronize(backend)
    return Fw
end
