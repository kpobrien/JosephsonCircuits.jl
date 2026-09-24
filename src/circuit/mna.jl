"""
    ismnaresistance(value)

Return `true` if `value` is a constant, real, finite, nonzero resistance,
which is what the port impedance scale of [`calcsolverscale`](@ref) is
taken from. Real numbers and complex numbers with zero imaginary part are
accepted; symbolic values, values with nonzero imaginary part, zeros, and
non-finite values return `false`.

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
        # symbolic values, such as frequency dependent resistances, are not
        # promoted; `Symbolics.Num` is a subtype of `Real`, so this test must
        # come first
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
    calcAmna(gaugeindices::Vector{Int}, Ntot::Int)

The constant gauge fixing rows of the augmented harmonic balance system of
[`hbnlsolve`](@ref): a one on the diagonal for each index in
`gaugeindices`, one per floating component of the static flux-stiffness
graph and zero-frequency mode (see [`calcdcgaugeindices`](@ref)), in an
`Ntot` square sparse matrix. Because the Kirchhoff current law equations of
a floating component are consistent but redundant at DC whenever the direct
current subsystem has a solution, which `dcpinning` checks, this rank-one
term renders the system nonsingular while the reference node flux is driven
to exactly zero and all original equations remain satisfied.

The augmented state is the `(Nnodes-1)*Nmodes` node fluxes followed by the
auxiliary variables of the mutually coupled inductor branches
([`calcAmnaind`](@ref)) and of the scattering block port currents, with the
mode index fastest. These equations are linear, so the matrix is constant
during the nonlinear solve and is its own contribution to the Jacobian.
"""
function calcAmna(gaugeindices::Vector{Int}, Ntot::Int)
    return sparse(gaugeindices, gaugeindices, ones(Complex{Float64},
        length(gaugeindices)), Ntot, Ntot)
end

"""
    mnapad(A::SparseMatrixCSC, Naux::Int)

Pad the sparse matrix `A` with `Naux` empty rows and columns, returning a
square matrix suitable for the augmented modified nodal analysis system.
"""
function mnapad(A::SparseMatrixCSC, Naux::Int)
    colptr = vcat(SparseArrays.getcolptr(A),
        fill(SparseArrays.getcolptr(A)[end], Naux))
    return SparseMatrixCSC(size(A, 1) + Naux, size(A, 2) + Naux,
        colptr, copy(rowvals(A)), copy(nonzeros(A)))
end

"""
    PaddedLinearTerm

The frequency dependent linear term of the augmented system, padded to the
auxiliary variables with the modified nodal analysis augmentation folded
in, together with the index maps which refill it from a new assembly of the
unpadded matrices. A sweep over component values builds this once and then
moves values only: the padding, the union with the augmentation and the
incidence matrix's empty columns are structure, which the values do not
change (see [`hbcache`](@ref)).

# Fields
- `Rbnm`: the incidence matrix with its empty auxiliary columns.
- `invLnm`, `Gnm`, `Cnm`: the padded matrices, `invLnm` with `Amna` added.
- `Amna`: the augmentation, whose value dependent entries are the coupled
    inductor rows of [`calcAmnaind`](@ref) and the scattering block term
    of [`scatteringlinearterm`](@ref), whose node flux columns carry the
    solver scale.
- `pind`: the positions of the coupled inductor rows' entries in
    `nonzeros(Amna)`.
- `pblock`: the positions of the scattering block term's entries in
    `nonzeros(Amna)`.
- `scale`: the solver scale the block term is stamped at.
- `pinv`, `pamna`: the positions of the unpadded inverse inductance
    entries and of `Amna`'s entries in `nonzeros(invLnm)`.
- `bnm`: the drive in the padded node basis.
- `wmodesm`, `wmodes2m`: the mode frequency diagonals over the padded
    system.
- `stampedblocks`: the scattering blocks as stamped, which the direct
    current path reads.
"""
struct PaddedLinearTerm{TR,TM,TA,TD,TS}
    Rbnm::TR
    invLnm::TM
    Gnm::TM
    Cnm::TM
    Amna::TA
    pind::Vector{Int}
    pblock::Vector{Int}
    scale::Base.RefValue{Float64}
    pinv::Vector{Int}
    pamna::Vector{Int}
    bnm::Vector{Complex{Float64}}
    wmodesm::TD
    wmodes2m::TD
    stampedblocks::TS
end

"""
    PaddedLinearTerm(Rbnm, invLnm, Gnm, Cnm, Amna, AmnaL, blockterm, scale,
        invLnm0, bnm, wmodesm, wmodes2m, stampedblocks)

Record a padded linear term as [`hbnlsolve`](@ref) assembled it, with the
unpadded inverse inductance `invLnm0`, the coupled inductor rows `AmnaL`
and the scattering block term `blockterm` (or `nothing` when the circuit
has no blocks) at the solver scale `scale` it was assembled from, so that
[`refill!`](@ref) can reproduce the assembly at new values. Returns
`nothing` when the coupled inductor rows or the block term share an entry
with the rest of the augmentation, which a refill by overwrite could not
reproduce exactly.
"""
function PaddedLinearTerm(Rbnm, invLnm, Gnm, Cnm, Amna, AmnaL, blockterm,
        scale, invLnm0, bnm, wmodesm, wmodes2m, stampedblocks)
    pind = sparseaddmap(Amna, AmnaL)
    pblock = isnothing(blockterm) ? Int[] : sparseaddmap(Amna, blockterm)
    # the coupled inductor rows and the block term are overwritten, not
    # added, on refill, which is exact only where nothing else is stored
    # under them
    all(k -> nonzeros(Amna)[pind[k]] == nonzeros(AmnaL)[k], eachindex(pind)) ||
        return nothing
    isnothing(blockterm) || (isdisjoint(pind, pblock) &&
        all(k -> nonzeros(Amna)[pblock[k]] == nonzeros(blockterm)[k],
            eachindex(pblock))) || return nothing
    pinv = sparseaddmap(invLnm,
        mnapad(invLnm0, size(invLnm, 1) - size(invLnm0, 1)))
    pamna = sparseaddmap(invLnm, Amna)
    return PaddedLinearTerm(Rbnm, invLnm, Gnm, Cnm, Amna, pind, pblock,
        Ref(Float64(scale)), pinv, pamna, bnm, wmodesm, wmodes2m,
        stampedblocks)
end

# the unpadded matrix has the padded one's structure in its leading block,
# entry for entry
function paddedstructure(A::SparseMatrixCSC, P::SparseMatrixCSC)
    n = size(A, 2)
    return nnz(A) == nnz(P) && rowvals(A) == rowvals(P) &&
        view(SparseArrays.getcolptr(A), 1:n+1) ==
        view(SparseArrays.getcolptr(P), 1:n+1)
end

"""
    refill!(lin::PaddedLinearTerm, invLnm, Gnm, Cnm, AmnaL, bbm, Rbnm0;
        blockterm = nothing, scale = lin.scale[])

Move the values of a new assembly of the unpadded matrices, of the coupled
inductor rows `AmnaL` (or `nothing` when there are none) and of the branch
drive `bbm` into the padded linear term, in place, and with them the
scattering block term `blockterm` stamped at the solver scale `scale`,
which is needed only when the scale moved. The result is entry for entry
what [`hbnlsolve`](@ref) assembles from the same inputs.
"""
function refill!(lin::PaddedLinearTerm, invLnm::SparseMatrixCSC,
        Gnm::SparseMatrixCSC, Cnm::SparseMatrixCSC, AmnaL, bbm, Rbnm0;
        blockterm::Union{Nothing,SparseMatrixCSC} = nothing,
        scale::Real = lin.scale[])
    for (A, P, what) in ((Gnm, lin.Gnm, "conductance"),
            (Cnm, lin.Cnm, "capacitance"), (invLnm, lin.invLnm, "inverse inductance"))
        (what == "inverse inductance" ? length(lin.pinv) == nnz(A) :
            paddedstructure(A, P)) || throw(ArgumentError(
            lazy"the $(what) matrix changed its sparse structure between points, which a reused linear term cannot follow; build a new one."))
    end
    copyto!(nonzeros(lin.Gnm), nonzeros(Gnm))
    copyto!(nonzeros(lin.Cnm), nonzeros(Cnm))
    if !isnothing(AmnaL)
        nnz(AmnaL) == length(lin.pind) || throw(ArgumentError(
            "the coupled inductor rows changed their structure between points; build a new linear term."))
        va = nonzeros(lin.Amna); vl = nonzeros(AmnaL)
        @inbounds for k in eachindex(lin.pind)
            va[lin.pind[k]] = vl[k]
        end
    end
    if !isnothing(blockterm)
        nnz(blockterm) == length(lin.pblock) || throw(ArgumentError(
            "the scattering block term changed its structure between points; build a new linear term."))
        va = nonzeros(lin.Amna); vb = nonzeros(blockterm)
        @inbounds for k in eachindex(lin.pblock)
            va[lin.pblock[k]] = vb[k]
        end
        lin.scale[] = scale
    end
    v = nonzeros(lin.invLnm)
    fill!(v, 0)
    va = nonzeros(lin.Amna)
    @inbounds for k in eachindex(lin.pamna)
        v[lin.pamna[k]] = va[k]
    end
    vi = nonzeros(invLnm)
    @inbounds for k in eachindex(lin.pinv)
        v[lin.pinv[k]] += vi[k]
    end
    Nnodal = size(Rbnm0, 2)
    mul!(view(lin.bnm, 1:Nnodal), transpose(Rbnm0), bbm)
    fill!(view(lin.bnm, Nnodal+1:length(lin.bnm)), 0)
    return lin
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
slipped past the direct current subsystem's solvability check) into the arbitrary
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
        gaugeindices::Vector{Int}, Nnodal::Int, bnm::AbstractVector, atol)

Validate the original, ungauged Kirchhoff current law equations at a
converged solution by reconstructing their residuals with
[`mnaungaugedkcl`](@ref) and comparing their infinity norm against a
block-relative infinity-norm tolerance,

`10*atol*(1 + norm(bnm[1:Nnodal], Inf)),`

so both sides have the same per-row interpretation and the accepted error
in any one equation does not grow with the number of driven rows. The
tolerance is deliberately independent of the achieved augmented residual
(which would be circular) and of the auxiliary current entries of the
state, which are not Kirchhoff current law quantities. A non-finite reconstructed norm, or a non-finite source scale
(which would make the tolerance infinite and accept anything), fails the
validation. Returns `(ok, normkcl, kcltol)` so a diagnostic can report
the achieved residual against the applied tolerance.
"""
function mnavalidatekcl(F::AbstractVector, x::AbstractVector,
    gaugeindices::Vector{Int}, Nnodal::Int, bnm::AbstractVector, atol)
    Fkcl = mnaungaugedkcl(F, x, gaugeindices, Nnodal)
    normkcl = norm(Fkcl, Inf)
    T = real(eltype(Fkcl))
    sourcescale = norm(view(bnm, 1:Nnodal), Inf)
    kcltol = 10*atol*(one(T) + sourcescale)
    ok = isfinite(normkcl) && isfinite(sourcescale) && normkcl <= kcltol
    return ok, normkcl, kcltol
end

"""
    auxcurrentscale(Lscale)

The scale of the auxiliary port current unknowns of the scattering blocks in
the linearized system, one over the solver inductance scale (see
[`calcsolverscale`](@ref)).

A node unknown is a flux and an auxiliary unknown is a current, and flux over
current is an inductance, so this writes the two in the same units, which
puts the auxiliary columns of the system on the scale of the nodal ones.
The nonlinear solve carries the same relation on its node fluxes, since
those are its Newton unknowns; the linearized solve carries it on the port
currents, which nothing outside the solve reads. A circuit with no scale to
take leaves the unknowns as they are.
"""
function auxcurrentscale(Lscale)
    L = abs(float(real(Lscale)))
    return (isfinite(L) && L > 0) ? 1/L : 1.0
end

"""
    calcsolverscale(w, componenttypes::Vector{Symbol}, vvn::Vector,
        portimpedances::Vector, Lscale)

Calculate the inductance scale used to nondimensionalize the nonlinear
harmonic balance system: the Kirchhoff current law rows are multiplied by
this scale (divided by the reduced flux quantum), the Josephson terms enter
as ratios of this scale to the junction inductances, and the auxiliary
variables of the modified nodal analysis formulation are branch currents in
units of the corresponding natural current scale. The scale is

`Lscale = Z0/w0`

with `Z0` the geometric mean of the constant real port reference impedances
(falling back to the geometric mean of all constant real resistors and then
to 50 ohms when none are present) and `w0` the geometric mean of the
absolute values of the nonzero drive frequencies. With this choice the
natural current unit is `phi0*w0/Z0`, the entries of the scaled system are
dimensionless and of order one for circuits driven near their characteristic
impedance and frequency, the auxiliary branch currents have magnitudes
comparable to the node fluxes, in circuits without inductors too, and the
residual tolerance `atol` becomes independent of the unit system of the
problem. Because the scale
multiplies rows only, and the auxiliary variables are internal, the returned
node fluxes and all physical quantities are unchanged in exact arithmetic.

When every drive frequency is zero the mean inductance `Lscale` is returned
instead, since no frequency scale is available.
"""
function calcsolverscale(w, componenttypes::Vector{Symbol}, vvn::Vector,
    portimpedances::Vector, Lscale)

    # geometric mean of the constant real port impedances, falling back to
    # all constant real resistors, then to 50 ohms.
    logsum = 0.0
    n = 0
    for z in portimpedances
        if ismnaresistance(z)
            logsum += log(abs(mnaresistance(z)))
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
        return Lscale
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
    checkcoupleddiagonal(Mb)
    I, J, _ = findnz(Mb)
    return sort(unique(vcat(I, J)))
end

function checkcoupleddiagonal(Mb::SparseMatrixCSC)
    nnz(Mb) == 0 && return nothing
    for j in axes(Mb, 2), k in nzrange(Mb, j)
        if rowvals(Mb)[k] == j && !iszero(nonzeros(Mb)[k])
            throw(ArgumentError("Mutual coupling between inductors which "*
                "share the same branch (the same pair of nodes) is not "*
                "supported: the parallel inductors are combined into a "*
                "single branch inductance before the coupling is applied, "*
                "which silently misrepresents the coupled pair. Route one "*
                "of the coupled inductors through an intermediate node so "*
                "the two inductors occupy distinct branches."))
        end
    end
    return nothing
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
inductors consistently with the node fluxes in `x`, by solving the branch
inductance system `(L/Lscale)*u = Rbn*phi` over the coupled branches, with
one sparse factorization for every mode, which zeros their constitutive
rows exactly. If the branch
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
    # the branch inductance matrix over the coupled branches, scaled: the
    # diagonal and the couplings, sparse, factorized once for every mode.
    # The coupled branches are sorted, so a branch's position is a search.
    I = Int[]; J = Int[]; V = Complex{Float64}[]
    for (r, b) in enumerate(coupledbranches)
        push!(I, r); push!(J, r); push!(V, Lb[b]/Lscale)
        for ptr in nzrange(Mb, b)
            k = rowvals(Mb)[ptr]
            k == b && continue
            push!(I, searchsortedfirst(coupledbranches, k)); push!(J, r)
            push!(V, nonzeros(Mb)[ptr]/Lscale)
        end
    end
    # a singular matrix, a unit coupling, leaves the currents as they are
    F = try
        lu(sparse(I, J, V, nb, nb))
    catch
        return x
    end
    # the branch fluxes of the coupled branches, a column per mode
    Rnb = sparse(transpose(Rbn))
    phib = zeros(Complex{Float64}, nb, Nmodes)
    for m in 1:Nmodes, (r, b) in enumerate(coupledbranches)
        for ptr in nzrange(Rnb, b)
            p = rowvals(Rnb)[ptr]
            phib[r, m] += nonzeros(Rnb)[ptr]*x[(p-1)*Nmodes + m]
        end
    end
    u = F \ phib
    for m in 1:Nmodes
        all(isfinite, view(u, :, m)) || continue
        for r in 1:nb
            x[auxoffset + (r-1)*Nmodes + m] = u[r, m]
        end
    end
    return x
end

# the error for inductors which share the branch of a coupling, raised
# where the mutual stamps are checked
function throwsharedinductors(offending)
    throw(ArgumentError("The inductors "*join(offending, ", ")*
        " share the same branch (the same pair of nodes), and at "*
        "least one inductor on that branch participates in mutual "*
        "coupling. Inductors on a shared branch are combined into "*
        "a single branch inductance before the mutual coupling is "*
        "applied, which misrepresents the coupled system. Route "*
        "the inductors through intermediate nodes so each "*
        "mutually coupled inductor occupies its own branch."))
end
