
"""
    modecouplingmask(Nmodes::Integer, couplingmodes)

Return the `Nmodes` x `Nmodes` boolean matrix of mode coupling blocks retained
by the preconditioner of [`ModeCouplingPreconditioner`](@ref): the block
coupling column mode `m2` into row mode `m1` is kept when `m2` is one of the
`couplingmodes` or when `m1 == m2`.

In the block partition into the retained set `S = couplingmodes` and the shell
`E` (everything else) this is

    [A_SS   0  ]
    [A_ES   D_E]

which is block lower triangular: the soft block with all of its internal
coupling, the mode diagonal blocks of the stiff shell, and the one way
coupling which carries the soft correction into the stiff equations. It is
therefore the multiplicative block Gauss-Seidel sweep (solve the soft block,
then each stiff mode against a residual updated by what the soft block just
did) expressed as a matrix, and a sparse LU of it *is* that sweep: the block
triangular form permutation of KLU discovers the structure and factorizes only
the diagonal blocks, so the factorization costs one factorization of `A_SS`
plus one small factorization per stiff mode, and the forward substitution of
the solve applies `A_ES`.

The two extremes are the endpoints of the accuracy/cost trade: `couplingmodes`
containing every mode keeps every coupling and gives the full Jacobian (an
exact preconditioner and a direct solve), while an empty `couplingmodes` keeps
only the mode diagonal and gives the block diagonal (one small factorization
per mode and no coupling at all).

The upper triangle `A_SE` is the part which is dropped. That choice, rather
than the transposed one, is what makes the stiffness of `E` useful: the error
of the preconditioned operator carries a factor of the stiff shell's inverse,
which is small precisely because those modes are stiff. See
[`selectcouplingmodes`](@ref) for how the set is chosen, and
[`modesoftness`](@ref) for the measurement it is chosen from.

# Examples
```jldoctest
julia> JosephsonCircuits.modecouplingmask(3, [2])
3×3 Matrix{Bool}:
 1  1  0
 0  1  0
 0  1  1

julia> JosephsonCircuits.modecouplingmask(3, Int[])
3×3 Matrix{Bool}:
 1  0  0
 0  1  0
 0  0  1
```
"""
function modecouplingmask(Nmodes::Integer, couplingmodes)
    Nmodes >= 1 || throw(ArgumentError(
        lazy"`Nmodes` = $(Nmodes) must be at least 1."))
    issoft = falses(Nmodes)
    for k in couplingmodes
        1 <= k <= Nmodes || throw(ArgumentError(
            lazy"the soft mode index $(k) is outside 1:$(Nmodes)."))
        issoft[k] = true
    end
    return [issoft[m2] || m1 == m2 for m1 in 1:Nmodes, m2 in 1:Nmodes]
end

"""
    restrictmodecoupling(Amatrixindices::Matrix, keep::AbstractMatrix{Bool})

Return a copy of a mode coupling index matrix (`Amatrixindices` or
`Amatrixconjindices`, see [`hbmatind`](@ref)) with the couplings not selected
by `keep` set to zero.

This is the whole of the coarse frequency grid transformation. A zero entry of
these matrices already means "this coupling falls outside the retained grid and
is dropped", which [`planrealjacobian`](@ref) honors when it builds the
sparsity structure (through [`activemoderows`](@ref)) and when it emits the
scatter lists. Zeroing entries therefore produces a Jacobian assembly plan for
the mode restricted operator directly: the restricted Jacobian is *materialized*
with the restricted structure rather than assembled in full and masked
afterwards, so the saving is realized in the sparsity, the fill in, and the
factorization, not just in the arithmetic.

The frequency dependent linear terms `invLnm`, `Gnm` and `Cnm` (including the
modified nodal analysis augmentation and the gauge fixing equations) are mode
diagonal by construction, so they are unaffected by the restriction and are
kept exactly. All of the mode coupling of the Jacobian lives in these two index
matrices.
"""
function restrictmodecoupling(Amatrixindices::Matrix,
    keep::AbstractMatrix{Bool})
    size(Amatrixindices) == size(keep) || throw(DimensionMismatch(
        lazy"the index matrix is $(size(Amatrixindices)) but the mask is $(size(keep))."))
    return [keep[i, j] ? Amatrixindices[i, j] : zero(eltype(Amatrixindices))
        for i in axes(Amatrixindices, 1), j in axes(Amatrixindices, 2)]
end

"""
    modeslotindex(layout::ModeLayout)

The vector, of length `layout.rdim`, giving the mode index of each slot of the
real representation. This is the block assignment
[`blockdiagonallayout`](@ref) needs to see the mode block diagonal
preconditioner as a uniform batch.
"""
function modeslotindex(layout::ModeLayout)
    nmodes = layout.nmodes
    return Int[(Int(layout.inv[j]) - 1) % nmodes + 1 for j in 1:layout.rdim]
end

"""
    ModeCouplingPreconditioner

A preconditioner for the matrix-free Newton-Krylov solve of the harmonic
balance system: the Jacobian materialized with its mode coupling restricted to
a selected set of modes and reduced to the mode diagonal everywhere else (see
[`modecouplingmask`](@ref)), then factorized.

The restriction is applied *while the Jacobian is materialized*, not afterwards:
the assembly plan is built from mode coupling index matrices whose dropped
couplings have been zeroed by [`restrictmodecoupling`](@ref), so the restricted
operator is assembled directly into its own, sparser, structure by the same
[`assemblerealjacobian!`](@ref) which assembles the full Jacobian, from the same
Fourier coefficients of `cos(phi(t))` and the same linear term matrices. No
submatrices are extracted and no sweep is written by hand; the block triangular
structure is what makes the sparse LU cheap, and the sparse triangular solve is
what applies the Gauss-Seidel sweep.

The default is the empty coupling set, that is, the mode block diagonal. Its
factorization is a batch of small independent per mode factorizations rather
than one large sparse factorization, which is what makes it scale. On a
strongly pumped device the block diagonal alone stalls, and
[`escalatepreconditioner!`](@ref) is the safety net.

# Fields
- `P`: the restricted Jacobian, in the real representation.
- `plan`: its [`RealJacobianPlan`](@ref).
- `sys`: the [`HBSystem`](@ref) the Jacobian is assembled from.
- `cache`: the [`FactorizationCache`](@ref) holding the factorization of `P`.
- `factorization`: the [`Factorization`](@ref) used to factorize `P`.
- `couplingmodes`: the modes whose coupling is retained in full.
- `updates`: the number of times the factorization has been rebuilt.
- `escalations`: the number of times the coupling set has been grown.
"""
mutable struct ModeCouplingPreconditioner{TS,TB} <: AbstractPreconditioner
    P::SparseMatrixCSC{Float64,Int}
    plan::RealJacobianPlan
    const sys::TS
    const cache::FactorizationCache
    const factorization::Factorization
    # rebuilds (P, plan) for a coupling set, closing over the plan ingredients
    # so the set can be grown after construction
    const build::TB
    const Nmodes::Int
    couplingmodes::Vector{Int}
    updates::Int
    escalations::Int
end

"""
    ModeCouplingPreconditioner(sys::HBSystem, Amatrixindices::Matrix,
        Amatrixconjindices::Matrix, Ljb::SparseVector, Lmean,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC, layout::ModeLayout;
        couplingmodes = :none, factorization = KLUfactorization())

Build a [`ModeCouplingPreconditioner`](@ref) from the same ingredients
[`planrealjacobian`](@ref) takes. `couplingmodes` may be

- `:none` (the default): the mode block diagonal,
- `:all`: the full Jacobian, which is an exact preconditioner,
- an `AbstractVector{<:Integer}`: exactly these mode indices.
"""
function ModeCouplingPreconditioner(sys, Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, Ljb::SparseVector, Lmean,
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC, layout::ModeLayout; couplingmodes = :none,
    factorization = KLUfactorization())

    build(S) = planrealjacobian(
        restrictmodecoupling(Amatrixindices, modecouplingmask(Nmodes, S)),
        restrictmodecoupling(Amatrixconjindices, modecouplingmask(Nmodes, S)),
        Ljb, Lmean, Rbnm, Nmodes, Nbranches, Nfreq, invLnm, Gnm, Cnm,
        layout, layout)

    S = if couplingmodes === :none
        Int[]
    elseif couplingmodes === :all
        collect(1:Nmodes)
    elseif couplingmodes isa AbstractVector{<:Integer}
        sort!(unique(Vector{Int}(couplingmodes)))
    else
        throw(ArgumentError(
            lazy"`couplingmodes` = $(couplingmodes) must be `:none`, `:all`, or a vector of mode indices."))
    end

    P, plan = build(S)
    return ModeCouplingPreconditioner(P, plan, sys, FactorizationCache(),
        factorization, build, Int(Nmodes), S, 0, 0)
end

"""
    escalatepreconditioner!(pc::ModeCouplingPreconditioner)

Grow the coupling set to every mode and rebuild, returning `true` if it grew
and `false` if it was already the full set.

This is the safety net of the block diagonal. A block diagonal preconditioner
is cheap, but on a strongly pumped device it leaves the preconditioned operator
nearly singular in a few directions, which stalls GMRES: the residual
polynomial is pinned at `p(0) = 1` and so cannot be made small near the origin.

Escalation goes straight to the full Jacobian rather than growing the set
gradually, because there is no reliable way to know in advance *which* modes
carry those directions. Criteria based on the linear response of the circuit,
on the mode frequencies, and on the mode diagonal blocks were each measured and
none generalized across devices: the deficiency is a specific direction inside
a mode's subspace rather than a property of the mode, and any per mode score
averages it away. The full set is exact, so the method is never less robust
than a direct solve, only faster when the block diagonal suffices. In practice
this fires once or twice on a strongly pumped line and not at all otherwise.

See [`RecyclingPreconditioner`](@ref) for the alternative, which absorbs the
same deficiency by measuring the directions rather than enlarging the
factorization.
"""
function escalatepreconditioner!(pc::ModeCouplingPreconditioner)
    length(pc.couplingmodes) >= pc.Nmodes && return false
    S = collect(1:pc.Nmodes)
    pc.P, pc.plan = pc.build(S)
    pc.couplingmodes = S
    pc.escalations += 1
    # the cache holds a symbolic factorization of the previous sparsity
    # structure, which no longer applies
    pc.cache.factorization = nothing
    return true
end

function updatepreconditioner!(pc::ModeCouplingPreconditioner,
    x::AbstractVector)
    setpoint!(pc.sys, x)
    jacobian!(pc.P, pc.plan, pc.sys)
    tryfactorize!(pc.cache, pc.factorization, pc.P)
    pc.updates += 1
    return pc
end

function applypreconditioner!(z::AbstractVector,
    pc::ModeCouplingPreconditioner, r::AbstractVector)
    return trysolve!(z, pc.cache.factorization, r)
end
