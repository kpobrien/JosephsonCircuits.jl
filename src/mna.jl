
"""
    ismnaresistance(value)

Return `true` if `value` is a constant, real, finite, nonzero resistance
suitable for promotion to an auxiliary branch current variable in the
modified nodal analysis (MNA) formulation used by [`hbnlsolve`](@ref) and
[`hblinsolve`](@ref). For example, real numbers and complex numbers with zero
imaginary part are accepted. Symbolic values, values with nonzero imaginary
part, zeros, and non-finite values return `false` and will remain as node
conductances.

# Examples
```jldoctest
julia> JosephsonCircuits.ismnaresistance(50.0)
true

julia> JosephsonCircuits.ismnaresistance(50.0+0.0im)
true

julia> JosephsonCircuits.ismnaresistance(50.0+1.0im)
false

julia> JosephsonCircuits.ismnaresistance(0.0)
false
```
"""
function ismnaresistance(value)
    if checkissymbolic(value)
        # symbolic values, such as frequency dependent resistances, are
        # not promoted. note that Symbolics.Num is a subtype of Real, so
        # this check must come first.
        return false
    elseif value isa Real
        return isfinite(value) && !iszero(value)
    elseif value isa Complex
        return iszero(imag(value)) && isfinite(real(value)) &&
            !iszero(real(value))
    else
        return false
    end
end

"""
    mnaresistance(value)

Return the real resistance of a value accepted by [`ismnaresistance`](@ref).
"""
mnaresistance(value::Real) = value
mnaresistance(value::Complex) = real(value)

"""
    mnaportresistorindices(componenttypes::Vector{Symbol},
        nodeindices::Matrix{Int}, mutualinductorbranchnames::Vector,
        vvn::Vector)

Return the sorted indices of the resistors which are promoted to auxiliary
branch current variables for MNA: the port resistors (resistors sharing a
branch with a port, see [`calcportimpedanceindices`](@ref)) whose value is
accepted by [`ismnaresistance`](@ref). All other resistors remain as node
conductances. Eliminating the auxiliary variable of a promoted resistor
recovers its conductance stamp exactly, so the two forms are algebraically
equivalent wherever the nodal form is well posed, and leaving the interior
resistors in the conductance matrix keeps the augmented system small. The
port resistors are promoted so their branch currents are explicit solver
unknowns; the hybrid form of their constitutive equations is the stamp
which generalizes to elements a conductance cannot express, such as
voltage sources and scattering parameter blocks.

# Examples
```jldoctest
julia> JosephsonCircuits.mnaportresistorindices([:P,:R,:C,:Lj,:R],[2 2 2 3 3;1 1 3 1 1],[],[1,50.0,1e-12,1e-9,25.0])
1-element Vector{Int64}:
 2

julia> JosephsonCircuits.mnaportresistorindices([:P,:R,:C,:Lj],[2 2 2 3;1 1 3 1],[],Any[1,50.0+0.0im,1e-12,1e-9])
1-element Vector{Int64}:
 2

julia> JosephsonCircuits.mnaportresistorindices([:R,:C,:Lj],[2 2 3;1 3 1],[],[50.0,1e-12,1e-9])
Int64[]
```
"""
function mnaportresistorindices(componenttypes::Vector{Symbol},
    nodeindices::Matrix{Int}, mutualinductorbranchnames::Vector,
    vvn::Vector)
    portimpedanceindices = calcportimpedanceindices(componenttypes,
        nodeindices, mutualinductorbranchnames, vvn)
    mnaindices = Int[]
    for i in portimpedanceindices
        if ismnaresistance(vvn[i])
            push!(mnaindices, i)
        end
    end
    return sort!(mnaindices)
end

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
zero-frequency mode, see [`calcdcgaugeindices`](@ref). The net direct current
injected into each floating component must be zero for a periodic solution to
exist; see [`checkdcsourcecompatibility`](@ref).

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


    # union-find over the "one indexed" nodes, with path halving
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
Symbolic values are accepted under the documented assuming that their
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
                        "make the scaling inductance Lmean infinite. An "*
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
    scale = sum(abs, terms)
    n = max(1, 2*length(terms) - 1)
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

"""
    checkdcsourcecompatibility(floatingcomponents::Vector{Vector{Int}},
        bnm::Vector, wmodes::Vector, Nmodes::Int, nodenames::Vector{String})

Check that the net direct current injected into each floating component of
the static flux-stiffness graph is zero, and throw an `ArgumentError`
otherwise. At zero frequency, capacitor currents and resistor currents
vanish in the periodic steady state (the DC mode is the zero Fourier
coefficient of a periodic flux, so the DC voltage `d(phi)/dt` is zero), and
inductor and junction branches internal to the component transport current
only between its nodes. Consequently the sum of the source terms over the
nodes of a floating component must vanish for a solution to exist. Without
this check, the gauge fixing equation would absorb an incompatible net
current into the arbitrary flux reference and report an apparently converged
solution which violates the original Kirchhoff current law equations.

The comparison is relative: the magnitude of the sum is compared against a
summation-roundoff bound, `max(1024, 4*length(component))*eps` times the
sum of the magnitudes of the contributing source terms.
"""
function checkdcsourcecompatibility(floatingcomponents::Vector{Vector{Int}},
    bnm::Vector, wmodes::Vector, Nmodes::Int, nodenames::Vector{String})

    for component in floatingcomponents
        for m in 1:Nmodes
            if iszero(wmodes[m])
                s = zero(eltype(bnm))
                nrm = zero(real(eltype(bnm)))
                for p in component
                    v = bnm[(p-2)*Nmodes + m]
                    s += v
                    nrm += abs(v)
                end
                epsT = eps(one(real(eltype(bnm))))
                # a non-finite sum (from a non-finite source term) is
                # certainly incompatible; note Inf > Inf and NaN > tol are
                # both false, so this must be tested explicitly.
                if !isfinite(abs(s)) || abs(s) > max(1024, 4*length(component))*epsT*nrm
                    names = join([nodenames[p] for p in component], ", ")
                    throw(ArgumentError("No periodic solution exists: the "*
                        "floating inductive/Josephson subnetwork containing "*
                        "node(s) $(names) has nonzero net direct current "*
                        "injection. With dc = true the zero frequency mode "*
                        "is the zero Fourier coefficient of a periodic "*
                        "flux, so capacitors and resistors carry no DC "*
                        "current and the net DC current into a subnetwork "*
                        "with no inductive path to ground must be zero."))
                end
            end
        end
    end

    return nothing
end

"""
    calcAmna(mnaindices::Vector{Int}, nodeindices::Matrix{Int}, vvn::Vector,
        gaugeindices::Vector{Int}, wmodes::Vector, Nmodes::Int, Nnodes::Int,
        Lmean)

Calculate the constant sparse matrix which augments the harmonic balance
system with the modified nodal analysis (MNA) equations. For each resistor in
`mnaindices` and each mode `m` with frequency `w_m` we add an auxiliary
variable and the constitutive equation

`im*w_m*(Lmean/R)*(phi_n1 - phi_n2) - u = 0,`

which is the frequency domain form of Ohm's law `i = (dPhi/dt)/R` in the
scaled node flux basis. Because the node flux unknowns of the nonlinear solver
are normalized by the reduced flux quantum, `phi = Phi/phi0`, and the
Kirchhoff current law rows are scaled by `Lmean`, the auxiliary variable is
the correspondingly scaled branch current

`u = Lmean*i/phi0,`

where `i` is the physical branch current flowing from the first to the second
node of the resistor. The augmented state layout of the nonlinear solver is
therefore `x = [Phi/phi0; Lmean*i/phi0]` with the mode index fastest. (The
linearized solver uses a different, unscaled layout; see
[`calcAmnasplit`](@ref).) The auxiliary variable `u` enters the Kirchhoff
current law equations of the two nodes with coefficients `+1` and `-1`.
Eliminating `u` recovers exactly the conductance matrix contribution
`im*Gnm*wmodesm`, so the MNA formulation is algebraically equivalent to the
nodal formulation wherever the latter is well posed. At zero frequency the
constitutive equation reduces to `u = 0`, consistent with a periodic steady
state in which resistors carry no direct current.

A gauge fixing equation (a one on the diagonal) is added for each index in
`gaugeindices`, one per floating component of the static flux-stiffness
graph and zero-frequency mode; see [`calcdcgaugeindices`](@ref). Because
the Kirchhoff current law equations of a floating component are consistent
but redundant at DC whenever the compatibility condition of
[`checkdcsourcecompatibility`](@ref) holds, this rank-one term renders the
system nonsingular while the reference node flux is driven to exactly zero
and all original equations remain satisfied.

The layout of the augmented state vector is the `(Nnodes-1)*Nmodes` node
fluxes followed by the auxiliary variables, with the mode index fastest, in
the order of `mnaindices`. Because these equations are linear, this matrix
is constant during the nonlinear solve and is its own contribution to the
Jacobian.
"""
function calcAmna(mnaindices::Vector{Int}, nodeindices::Matrix{Int},
    vvn::Vector, gaugeindices::Vector{Int}, wmodes::Vector, Nmodes::Int,
    Nnodes::Int, Lmean)

    Nnodal = (Nnodes-1)*Nmodes
    Naux = length(mnaindices)*Nmodes
    Ntot = Nnodal + Naux

    I = Int[]
    J = Int[]
    V = Complex{Float64}[]

    for (r, ci) in enumerate(mnaindices)
        # nodeindices is "one indexed" so 1 is the ground node
        p1 = nodeindices[1, ci] - 1
        p2 = nodeindices[2, ci] - 1
        g = 1/mnaresistance(vvn[ci])
        for m in 1:Nmodes
            auxindex = Nnodal + (r-1)*Nmodes + m
            # Kirchhoff current law contributions of the branch current
            if p1 > 0
                push!(I, (p1-1)*Nmodes + m); push!(J, auxindex); push!(V, 1)
            end
            if p2 > 0
                push!(I, (p2-1)*Nmodes + m); push!(J, auxindex); push!(V, -1)
            end
            # constitutive equation. the signed mode frequency is used, which
            # is consistent with taking the complex conjugate of the entries
            # for modes with negative frequencies as conjnegfreq! does for
            # the conductance matrix.
            if p1 > 0
                push!(I, auxindex); push!(J, (p1-1)*Nmodes + m)
                push!(V, im*wmodes[m]*Lmean*g)
            end
            if p2 > 0
                push!(I, auxindex); push!(J, (p2-1)*Nmodes + m)
                push!(V, -im*wmodes[m]*Lmean*g)
            end
            push!(I, auxindex); push!(J, auxindex); push!(V, -1)
        end
    end

    # gauge fixing equations, one per floating component of the static
    # flux-stiffness graph and zero-frequency mode
    for c in gaugeindices
        push!(I, c); push!(J, c); push!(V, 1)
    end

    return sparse(I, J, V, Ntot, Ntot)
end

"""
    mnasubtractpromoted(Gnm::SparseMatrixCSC, Gnmp::SparseMatrixCSC)

Return a copy of the conductance matrix `Gnm` with the stamp `Gnmp` of the
promoted resistors subtracted from `Gnmp`. Every stored entry of `Gnmp` is
also a stored entry of `Gnm` because the promoted resistors are a subset of
the resistors stamped into `Gnm`, so the sparsity structure is unchanged.
Unlike the generic sparse subtraction, this works when `Gnm` has a
non-concrete element type, which occurs when the circuit also contains
symbolic (frequency dependent) resistors.
"""
function mnasubtractpromoted(Gnm::SparseMatrixCSC, Gnmp::SparseMatrixCSC)
    G = copy(Gnm)
    rows = rowvals(G)
    vals = nonzeros(G)
    I, J, V = findnz(Gnmp)
    for k in eachindex(I)
        # locate the entry within the existing structure of G and assert
        # its presence rather than allowing a silent structural insertion
        found = false
        for p in nzrange(G, J[k])
            if rows[p] == I[k]
                vals[p] -= V[k]
                found = true
                break
            end
        end
        if !found
            throw(ArgumentError("The stamp of a promoted resistor at "*
                "position ($(I[k]), $(J[k])) is not present in the "*
                "conductance matrix structure. This indicates an "*
                "inconsistency between the resistor stamps and their "*
                "promoted subset."))
        end
    end
    return G
end

"""
    mnapad(A::SparseMatrixCSC, Naux::Int)

Pad the sparse matrix `A` with `Naux` empty rows and columns, returning a
square matrix suitable for the augmented modified nodal analysis system.
"""
function mnapad(A::SparseMatrixCSC, Naux::Int)
    I, J, V = findnz(A)
    return sparse(I, J, V, size(A,1)+Naux, size(A,2)+Naux)
end

"""
    mnainitialaux!(x::AbstractVector, Amna::SparseMatrixCSC, Nnodal::Int)

Initialize the auxiliary current variables `x[Nnodal+1:end]` consistently
with the node fluxes `x[1:Nnodal]` by adding the residual of each
constitutive row to its auxiliary variable, which zeros those rows exactly.
Together with [`mnagaugenormalize!`](@ref), which must be applied first
when gauge fixing equations are present, this makes the augmented residual
at the initial value equal to the nodal Kirchhoff current law residual at
the corresponding node fluxes, so line searches, residual stopping tests,
and Anderson acceleration see the same initial merit as a purely nodal
formulation would when an initial guess `x0` is supplied.
"""
function mnainitialaux!(x::AbstractVector, Amna::SparseMatrixCSC,
    Nnodal::Int)
    t = Amna*x
    for i in Nnodal+1:length(x)
        # each constitutive row is c'*phi - u, so adding the row residual to
        # u makes the row exactly zero for any starting u.
        x[i] += t[i]
    end
    return x
end

"""
    mnagaugenormalize!(x::AbstractVector,
        floatingcomponents::Vector{Vector{Int}}, wmodes::Vector,
        Nmodes::Int)

Transform an initial guess into the gauge selected by the gauge fixing
equations by subtracting, for each floating component of the static
flux-stiffness graph and each zero-frequency mode, the flux of the
component's reference node from the fluxes of all nodes of the component.
A common shift of the DC fluxes of a floating component leaves every
branch flux, and therefore every physical circuit quantity and every
Kirchhoff current law residual, unchanged - it is exactly the gauge degree
of freedom - but it does enter the gauge fixing rows. Normalizing the
guess makes physically equivalent initial values produce identical
augmented residuals and makes the reference-node gauge rows exactly zero
at the initial value.
"""
function mnagaugenormalize!(x::AbstractVector,
    floatingcomponents::Vector{Vector{Int}}, wmodes::Vector, Nmodes::Int)
    for component in floatingcomponents
        pref = first(component)
        for m in 1:Nmodes
            if iszero(wmodes[m])
                offset = x[(pref-2)*Nmodes + m]
                if !iszero(offset)
                    for p in component
                        x[(p-2)*Nmodes + m] -= offset
                    end
                end
            end
        end
    end
    return x
end

"""
    mnaungaugedkcl(F::AbstractVector, x::AbstractVector,
        gaugeindices::Vector{Int}, Nnodal::Int)

Reconstruct the residuals of the original, ungauged Kirchhoff current law
equations from the augmented residual `F` and the state `x`. The gauge
fixing equations add `x[g]` to the augmented residual of each gauge row
`g`, so the physical residual is `F[g] - x[g]` at the gauge rows and
`F[i]` elsewhere in the node block. This is the quantity which must be
small for the reported solution to satisfy the original circuit
equations: a gauge equation can otherwise absorb an incompatibility (for
example a net direct current injected into a floating subnetwork which
slipped past [`checkdcsourcecompatibility`](@ref)) into the arbitrary
flux reference while the augmented residual converges to zero.
"""
function mnaungaugedkcl(F::AbstractVector, x::AbstractVector,
    gaugeindices::Vector{Int}, Nnodal::Int)
    Fkcl = Vector(view(F, 1:Nnodal))
    for g in gaugeindices
        Fkcl[g] -= x[g]
    end
    return Fkcl
end

"""
    mnavalidatekcl(F::AbstractVector, x::AbstractVector,
        gaugeindices::Vector{Int}, Nnodal::Int, bnm::AbstractVector, ftol)

Validate the original, ungauged Kirchhoff current law equations at a
converged solution by reconstructing their residuals with
[`mnaungaugedkcl`](@ref) and comparing their infinity norm against a
block-relative infinity-norm tolerance,

`10*ftol*(1 + norm(bnm[1:Nnodal], Inf)),`

so both sides have the same per-row interpretation and the accepted error
in any one equation does not grow with the number of driven rows. The
tolerance is deliberately independent of the achieved augmented residual
(which would be circular) and of the auxiliary current entries of the
state, which are not Kirchhoff current law quantities. (With the solver
inductance scale of [`calcsolverscale`](@ref) the auxiliary entries are of
order one; under the earlier mean-inductance scale they reached ~1e9 in
inductor free circuits, which motivated this exclusion.) A non-finite reconstructed norm, or a non-finite source scale
(which would make the tolerance infinite and accept anything), fails the
validation. Returns `(ok, normkcl, kcltol)` so a diagnostic can report
the achieved residual against the applied tolerance.
"""
function mnavalidatekcl(F::AbstractVector, x::AbstractVector,
    gaugeindices::Vector{Int}, Nnodal::Int, bnm::AbstractVector, ftol)
    Fkcl = mnaungaugedkcl(F, x, gaugeindices, Nnodal)
    normkcl = norm(Fkcl, Inf)
    T = real(eltype(Fkcl))
    sourcescale = norm(view(bnm, 1:Nnodal), Inf)
    kcltol = 10*ftol*(one(T) + sourcescale)
    ok = isfinite(normkcl) && isfinite(sourcescale) && normkcl <= kcltol
    return ok, normkcl, kcltol
end

"""
    calcAmnasplit(mnaindices::Vector{Int}, nodeindices::Matrix{Int},
        vvn::Vector, Nmodes::Int, Nnodes::Int)

Like [`calcAmna`](@ref) but for the linearized solver `hblinsolve`, where
the mode frequencies change at every signal frequency. Returns two
matrices: `Amna0` containing the frequency independent entries (the
Kirchhoff current law couplings of the auxiliary branch currents and the
minus identity of the constitutive equations), and `AmnaG` containing the
conductances of the constitutive equations, to be added per signal
frequency as `im*AmnaG*wmodesm` with `sparseaddconjsubst!`, which applies
the per-mode frequencies and the complex conjugation for negative frequency
modes exactly as for the conductance matrix `Gnm`. Eliminating the
auxiliary variables then recovers the (possibly conjugated) nodal
conductance stamp of each promoted resistor exactly. The auxiliary
variables are branch currents in the same normalization as the node fluxes
of the linearized solver (no `Lmean` or `phi0` scaling is applied, unlike
the nonlinear solver; see [`calcAmna`](@ref)), and the conductances are
computed as
`1/vvn[i]`, the same expression used by `calcGn`, so the subtraction of the
promoted resistors from `Gnm` cancels exactly. No gauge fixing equations
are needed because `hblinsolve` does not permit modes at (numerically)
zero total frequency; see [`isnumericallyzero`](@ref).
"""
function calcAmnasplit(mnaindices::Vector{Int}, nodeindices::Matrix{Int},
    vvn::Vector, Nmodes::Int, Nnodes::Int)

    Nnodal = (Nnodes-1)*Nmodes
    Ntot = Nnodal + length(mnaindices)*Nmodes
    I0 = Int[]; J0 = Int[]; V0 = Complex{Float64}[]
    IG = Int[]; JG = Int[]; VG = Complex{Float64}[]
    for (r, ci) in enumerate(mnaindices)
        # nodeindices is "one indexed" so 1 is the ground node
        p1 = nodeindices[1, ci] - 1
        p2 = nodeindices[2, ci] - 1
        g = 1/vvn[ci]
        for m in 1:Nmodes
            auxindex = Nnodal + (r-1)*Nmodes + m
            if p1 > 0
                push!(I0, (p1-1)*Nmodes + m); push!(J0, auxindex); push!(V0, 1)
                push!(IG, auxindex); push!(JG, (p1-1)*Nmodes + m); push!(VG, g)
            end
            if p2 > 0
                push!(I0, (p2-1)*Nmodes + m); push!(J0, auxindex); push!(V0, -1)
                push!(IG, auxindex); push!(JG, (p2-1)*Nmodes + m); push!(VG, -g)
            end
            push!(I0, auxindex); push!(J0, auxindex); push!(V0, -1)
        end
    end
    return (sparse(I0, J0, V0, Ntot, Ntot), sparse(IG, JG, VG, Ntot, Ntot))
end

"""
    calcsolverscale(w, componenttypes::Vector{Symbol}, vvn::Vector,
        portimpedanceindices::Vector{Int}, Lmean)

Calculate the inductance scale used to nondimensionalize the nonlinear
harmonic balance system: the Kirchhoff current law rows are multiplied by
this scale (divided by the reduced flux quantum), the Josephson terms enter
as ratios of this scale to the junction inductances, and the auxiliary
variables of the modified nodal analysis formulation are branch currents in
units of the corresponding natural current scale. The scale is

`Lscale = Z0/w0`

with `Z0` the geometric mean of the constant real port impedance resistors
(falling back to the geometric mean of all constant real resistors and then
to 50 ohms when none are present) and `w0` the geometric mean of the
absolute values of the nonzero drive frequencies. With this choice the
natural current unit is `phi0*w0/Z0`, the entries of the scaled system are
dimensionless and of order one for circuits driven near their characteristic
impedance and frequency, the auxiliary branch currents have magnitudes
comparable to the node fluxes (in particular in circuits without inductors,
where the previous mean-inductance scale degenerated to one henry and
produced auxiliary values of order 1e9), and the residual tolerance `ftol`
becomes independent of the unit system of the problem. Because the scale
multiplies rows only, and the auxiliary variables are internal, the returned
node fluxes and all physical quantities are unchanged in exact arithmetic.

When every drive frequency is zero the mean inductance `Lmean` is returned
instead, since no frequency scale is available.
"""
function calcsolverscale(w, componenttypes::Vector{Symbol}, vvn::Vector,
    portimpedanceindices::Vector{Int}, Lmean)

    # geometric mean of the constant real port impedances, falling back to
    # all constant real resistors, then to 50 ohms.
    logsum = 0.0
    n = 0
    for i in portimpedanceindices
        if ismnaresistance(vvn[i])
            logsum += log(abs(mnaresistance(vvn[i])))
            n += 1
        end
    end
    if n == 0
        for i in eachindex(componenttypes)
            if componenttypes[i] == :R && ismnaresistance(vvn[i])
                logsum += log(abs(mnaresistance(vvn[i])))
                n += 1
            end
        end
    end
    Z0 = n == 0 ? 50.0 : exp(logsum/n)

    # geometric mean of the absolute values of the nonzero drive
    # frequencies. non-finite frequencies are rejected by the solvers
    # before this function is called.
    logsum = 0.0
    n = 0
    for wi in w
        wr = abs(float(real(wi)))
        if !iszero(wr)
            logsum += log(wr)
            n += 1
        end
    end
    if n == 0
        return Lmean
    end
    w0 = exp(logsum/n)

    return Z0/w0
end

"""
    mnacoupledbranches(Mb::SparseMatrixCSC)

Return the sorted branch indices which participate in mutual inductive
coupling, the union of the row and column supports of the branch mutual
inductance matrix `Mb`. These branches are assigned auxiliary branch current
variables by the modified nodal analysis formulation instead of being
eliminated through the inverse of the branch inductance matrix, so the system
matrix entries remain bounded as the coupling coefficient approaches one,
where the inverse inductance entries of the nodal formulation diverge as
`1/(1-k^2)`.
"""
function mnacoupledbranches(Mb::SparseMatrixCSC)
    I, J, V = findnz(Mb)
    for k in eachindex(I)
        if I[k] == J[k] && !iszero(V[k])
            throw(ArgumentError("Mutual coupling between inductors which "*
                "share the same branch (the same pair of nodes) is not "*
                "supported: the parallel inductors are combined into a "*
                "single branch inductance before the coupling is applied, "*
                "which silently misrepresents the coupled pair. Route one "*
                "of the coupled inductors through an intermediate node so "*
                "the two inductors occupy distinct branches."))
        end
    end
    return sort(unique(vcat(I, J)))
end

"""
    mnadropbranches(Lb::SparseVector, branches::Vector{Int})

Return a copy of the branch inductance vector `Lb` with the entries at
`branches` removed, so the nodal inverse inductance matrix can be computed
from the remaining, uncoupled inductors only.
"""
function mnadropbranches(Lb::SparseVector, branches::Vector{Int})
    keep = [i for i in eachindex(Lb.nzind) if !(Lb.nzind[i] in branches)]
    return SparseArrays.sparsevec(Lb.nzind[keep], Lb.nzval[keep], length(Lb))
end

"""
    calcAmnaind(coupledbranches::Vector{Int}, Lb::SparseVector,
        Mb::SparseMatrixCSC, Rbn::SparseMatrixCSC, Nmodes::Int,
        auxoffset::Int, Ntot::Int, Lscale)

Calculate the constant sparse matrix which augments the harmonic balance
system with auxiliary branch current variables for the mutually coupled
inductors. For each coupled branch `b` and mode `m` an auxiliary variable `u`
is added at index `auxoffset + (r-1)*Nmodes + m` (with `r` the position of
`b` in `coupledbranches`) together with the branch flux constitutive equation

`sum_p Rbn[b,p]*phi_p - sum_k (L[b,k]/Lscale)*u_k = 0,`

where `L[b,k]` is the branch inductance matrix (the branch self inductances
on the diagonal and the mutual inductances `Mb` off the diagonal), and the
auxiliary variable enters the Kirchhoff current law equation of each node
`p` of the branch with coefficient `Rbn[b,p]`. In the scaled units of the
nonlinear solver the auxiliary variable is `u = Lscale*i/phi0` with `i` the
physical branch current in the orientation of the incidence matrix; the
linearized solver uses `Lscale = 1` and unscaled branch currents. All
entries are real and frequency independent. Eliminating the auxiliary
variables recovers exactly the coupled part of the nodal inverse inductance
stamp, `Lscale*Rbn'*inv(L)*Rbn` (unit tested as a Schur complement
identity), so the formulation is algebraically equivalent to the nodal one
wherever the branch inductance matrix is invertible. Its entries remain
bounded as the coupling coefficient approaches one, and unlike the nodal
formulation it remains well posed at perfect coupling (|k| = 1) whenever
the surrounding circuit determines the branch currents. Coupling between
inductors sharing a single branch is rejected with an informative error
(see [`mnacoupledbranches`](@ref)); a coupling matrix which leaves some
branch current combination physically undetermined would produce a
singular system caught at factorization.
"""
function calcAmnaind(coupledbranches::Vector{Int}, Lb::SparseVector,
    Mb::SparseMatrixCSC, Rbn::SparseMatrixCSC, Nmodes::Int,
    auxoffset::Int, Ntot::Int, Lscale)

    I = Int[]
    J = Int[]
    V = Complex{Float64}[]
    branchposition = Dict(b => r for (r, b) in enumerate(coupledbranches))
    # transpose so the nodes of each branch are a column
    Rnb = sparse(transpose(Rbn))
    Rnbrows = rowvals(Rnb)
    Rnbvals = nonzeros(Rnb)
    Mbrows = rowvals(Mb)
    Mbvals = nonzeros(Mb)
    for (r, b) in enumerate(coupledbranches)
        if iszero(Lb[b])
            throw(ArgumentError("The mutually coupled branch $(b) has no "*
                "self inductance. Mutual coupling requires the coupled "*
                "inductors to have finite, nonzero inductances."))
        end
        for m in 1:Nmodes
            aux = auxoffset + (r-1)*Nmodes + m
            # the Kirchhoff current law couplings of the branch current and
            # the branch flux entries of the constitutive equation
            for ptr in nzrange(Rnb, b)
                p = Rnbrows[ptr]
                s = Rnbvals[ptr]
                push!(I, (p-1)*Nmodes + m); push!(J, aux); push!(V, s)
                push!(I, aux); push!(J, (p-1)*Nmodes + m); push!(V, s)
            end
            # the branch inductance matrix row, scaled
            push!(I, aux); push!(J, aux); push!(V, -Lb[b]/Lscale)
            for ptr in nzrange(Mb, b)
                k = Mbrows[ptr]
                if k != b
                    kr = get(branchposition, k, 0)
                    if kr == 0
                        throw(ArgumentError("The mutual inductance between "*
                            "branches $(b) and $(k) references a branch "*
                            "outside the coupled set. This indicates an "*
                            "inconsistency in the mutual inductance matrix."))
                    end
                    kaux = auxoffset + (kr-1)*Nmodes + m
                    push!(I, aux); push!(J, kaux)
                    push!(V, -Mbvals[ptr]/Lscale)
                end
            end
        end
    end
    return sparse(I, J, V, Ntot, Ntot)
end

"""
    mnainitialauxind!(x::AbstractVector, coupledbranches::Vector{Int},
        Lb::SparseVector, Mb::SparseMatrixCSC, Rbn::SparseMatrixCSC,
        Nmodes::Int, auxoffset::Int, Lscale)

Initialize the auxiliary branch current variables of the mutually coupled
inductors consistently with the node fluxes in `x`, by solving the small
dense branch inductance system `(L/Lscale)*u = Rbn*phi` over the coupled
branches, which zeros their constitutive rows exactly. If the branch
inductance matrix is singular (a perfectly coupled pair, `|k| = 1`) the
auxiliary variables are left unchanged: the full system can still be well
posed and solvable in that case, because the constitutive equations use
the un-inverted branch inductance matrix, so only this warm start
refinement is skipped.
"""
function mnainitialauxind!(x::AbstractVector, coupledbranches::Vector{Int},
    Lb::SparseVector, Mb::SparseMatrixCSC, Rbn::SparseMatrixCSC,
    Nmodes::Int, auxoffset::Int, Lscale)

    nb = length(coupledbranches)
    nb == 0 && return x
    # the dense branch inductance matrix over the coupled branches, scaled
    L = zeros(Complex{Float64}, nb, nb)
    for (r, b) in enumerate(coupledbranches)
        L[r, r] = Lb[b]/Lscale
        for ptr in nzrange(Mb, b)
            k = rowvals(Mb)[ptr]
            if k != b
                kr = findfirst(==(k), coupledbranches)
                L[kr, r] = nonzeros(Mb)[ptr]/Lscale
            end
        end
    end
    # the branch fluxes of the coupled branches for each mode
    Rnb = sparse(transpose(Rbn))
    for m in 1:Nmodes
        phib = zeros(Complex{Float64}, nb)
        for (r, b) in enumerate(coupledbranches)
            for ptr in nzrange(Rnb, b)
                p = rowvals(Rnb)[ptr]
                phib[r] += nonzeros(Rnb)[ptr]*x[(p-1)*Nmodes + m]
            end
        end
        u = try
            L \ phib
        catch
            nothing
        end
        if !isnothing(u) && all(isfinite, u)
            for r in 1:nb
                x[auxoffset + (r-1)*Nmodes + m] = u[r]
            end
        end
    end
    return x
end

"""
    mnainitialauxall!(x::AbstractVector, Amna::SparseMatrixCSC, Nnodal::Int,
        Nauxr::Int, coupledbranches::Vector{Int}, Lb::SparseVector,
        Mb::SparseMatrixCSC, Rbn::SparseMatrixCSC, Nmodes::Int, Lscale)

Initialize all auxiliary variables consistently with the node fluxes in
`x`: the resistor auxiliary currents with [`mnainitialaux!`](@ref)
restricted to the resistor rows (whose constitutive equations have a unit
self coefficient and no cross coupling) and the coupled inductor auxiliary
currents with [`mnainitialauxind!`](@ref).
"""
function mnainitialauxall!(x::AbstractVector, Amna::SparseMatrixCSC,
    Nnodal::Int, Nauxr::Int, coupledbranches::Vector{Int},
    Lb::SparseVector, Mb::SparseMatrixCSC, Rbn::SparseMatrixCSC,
    Nmodes::Int, Lscale)
    # each resistor constitutive row is c'*phi - u, so adding the row
    # residual to u makes the row exactly zero for any starting u. this is
    # not applied to the inductor rows, whose self coefficients are not
    # unity and which couple other auxiliary variables.
    t = Amna*x
    for i in Nnodal+1:Nnodal+Nauxr
        x[i] += t[i]
    end
    mnainitialauxind!(x, coupledbranches, Lb, Mb, Rbn, Nmodes,
        Nnodal + Nauxr, Lscale)
    return x
end

"""
    checkcoupledbranchinductors(componentnames::Vector,
        componenttypes::Vector{Symbol}, nodeindices::Matrix,
        edge2indexdict::Dict, Mb::SparseMatrixCSC)

Check that no branch which participates in mutual inductive coupling hosts
more than one inductor. Inductors which share a branch (the same pair of
nodes) are combined into a single branch inductance by reciprocal sum
before the mutual coupling is applied - exact for uncoupled parallel
inductors, which share the same branch flux, but a misrepresentation as
soon as any inductor on the branch is mutually coupled: the correct
effective mutual coupling of the merged branch differs from the stamped
one (for an uncoupled `L2` sharing a branch with a coupled `L1`, by the
current division factor `L2/(L1+L2)`). An informative `ArgumentError`
naming the inductors is thrown, directing the user to route the inductors
through intermediate nodes so each mutually coupled inductor occupies its
own branch. Called by [`numericmatrices`](@ref) so both solvers and direct
users of the circuit matrices are protected; see also
[`mnacoupledbranches`](@ref), which rejects coupling between two inductors
on the same branch (a diagonal mutual inductance entry).
"""
function checkcoupledbranchinductors(componentnames::Vector,
    componenttypes::Vector{Symbol}, nodeindices::Matrix,
    edge2indexdict::Dict, Mb::SparseMatrixCSC)

    I, J, V = findnz(Mb)
    coupled = Set{Int}()
    for k in eachindex(I)
        if !iszero(V[k])
            push!(coupled, I[k])
            push!(coupled, J[k])
        end
    end
    isempty(coupled) && return nothing
    counts = Dict{Int,Int}()
    for i in eachindex(componenttypes)
        if componenttypes[i] == :L
            b = edge2indexdict[(nodeindices[1,i], nodeindices[2,i])]
            if b in coupled
                counts[b] = get(counts, b, 0) + 1
            end
        end
    end
    for (b, c) in counts
        if c > 1
            offending = [componentnames[i] for i in eachindex(componenttypes)
                if componenttypes[i] == :L &&
                edge2indexdict[(nodeindices[1,i], nodeindices[2,i])] == b]
            throw(ArgumentError("The inductors "*join(offending, ", ")*
                " share the same branch (the same pair of nodes), and at "*
                "least one inductor on that branch participates in mutual "*
                "coupling. Inductors on a shared branch are combined into "*
                "a single branch inductance before the mutual coupling is "*
                "applied, which misrepresents the coupled system. Route "*
                "the inductors through intermediate nodes so each "*
                "mutually coupled inductor occupies its own branch."))
        end
    end
    return nothing
end
