
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
which is small precisely because those modes are stiff. Which modes go in `S`
is the caller's, not chosen here: see the `couplingmodes` argument of
[`ModeCouplingPreconditioner`](@ref).

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
is dropped", which [`planstructurerealjacobian`](@ref) honors when it builds
the sparsity structure (through [`activemoderows`](@ref)) and when it assembles
the values. Zeroing entries therefore produces a Jacobian assembly plan for the
mode restricted operator directly: the restricted Jacobian is *materialized*
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
- `P`: the restricted Jacobian's sparsity structure. A `SparseMatrixCSC` on
    the host, and on a backend a [`DeviceSparsePattern`](@ref) of the
    transpose, which is the row major structure a device sparse matrix wants.
- `sys`: the [`HBSystem`](@ref) the Jacobian is assembled from.
- `cache`: the [`FactorizationCache`](@ref) holding the factorization of `P`.
- `basefactorization`: the [`Factorization`](@ref) the caller asked for.
- `factorization`: the one actually used for the current `P`, which is the
    batched form of `basefactorization` while `P` is the mode block diagonal
    and that form exists. It is reselected whenever `P` is rebuilt, because
    escalation destroys the batch structure.
- `couplingmodes`: the modes whose coupling is retained in full, or a mask of
    the retained couplings when one was given directly.
- `updates`: the number of times the factorization has been rebuilt.
- `escalations`: the number of times the coupling set has been grown.
"""
mutable struct ModeCouplingPreconditioner{TS,TB} <: AbstractPreconditioner
    # loosely typed, like the device half below: on a backend `P` is a
    # `DeviceSparsePattern` rather than a host sparse matrix
    P
    const sys::TS
    const cache::FactorizationCache
    const basefactorization::Factorization
    factorization::Factorization
    # rebuilds (P, plan) for a coupling set, closing over the plan ingredients
    # so the set can be grown after construction
    const build::TB
    const Nmodes::Int
    couplingmodes
    updates::Int
    escalations::Int
    # the assembly and the values it writes. Loosely typed, as the rest of
    # this struct is, because they are rebuilt on escalation and the concrete
    # index type may change with the size of the escalated Jacobian. On a host
    # `nzval` aliases the stored values of `P`.
    deviceplan
    nzval
end

"""
    ModeCouplingPreconditioner(sys::HBSystem, Amatrixindices::Matrix,
        Amatrixconjindices::Matrix, Ljb::SparseVector, Lmean,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC, layout::ModeLayout;
        couplingmodes = :none, factorization = KLUfactorization())

Build a [`ModeCouplingPreconditioner`](@ref) from the same ingredients
[`planstructurerealjacobian`](@ref) takes. `couplingmodes` may be

- `:none` (the default): the mode block diagonal,
- `:all`: the full Jacobian, which is an exact preconditioner,
- an `AbstractVector{<:Integer}`: exactly these mode indices.
"""
function ModeCouplingPreconditioner(sys, Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, Ljb::SparseVector, Lmean,
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC, layout::ModeLayout; couplingmodes = :none,
    factorization = KLUfactorization(),
    precision::Union{Nothing,Type{<:AbstractFloat}} = nothing)

    # On a device backend the Jacobian is built transposed, because its
    # stored order is then the row major order a device sparse matrix and a
    # device direct solver want. Producing it that way costs nothing and
    # removes the permutation, the index the assembly kernel would go through
    # to apply it, and the loss of coalescing that indirection caused.
    # On a device backend there is no segmented gather at all: the assembly
    # reads the circuit's structure directly, so what is built here is the
    # sparsity structure and the linear term index maps, both cheap. The
    # Jacobian is
    # built transposed at the same time, because its stored order is then the
    # row major order a device sparse matrix and a device direct solver want.
    # One path, on every backend. The two differ only in orientation and in
    # what the factorization is handed: a device sparse matrix is compressed by
    # rows, so there the Jacobian is built transposed and its values live on the
    # backend, while a host factorization wants the matrix itself and the
    # assembly writes straight into its stored values.
    function build(S)
        keep = S isa AbstractMatrix{Bool} ? S : modecouplingmask(Nmodes, S)
        Ami = restrictmodecoupling(Amatrixindices, keep)
        Amc = restrictmodecoupling(Amatrixconjindices, keep)
        backend = sys.nonlineartermplan.backend
        transposed = !(backend isa CPU)
        Tv = isnothing(precision) ?
            real(promote_type(typeof(Lmean),
                isempty(Ljb.nzval) ? typeof(Lmean) : real(eltype(Ljb)))) :
            precision

        P, nodesandsigns = realjacobianstructure(Ami, Amc, Ljb, Rbnm, Nmodes,
            Nbranches, invLnm, Gnm, Cnm, layout, layout, Tv;
            transposed = transposed, backend = backend)

        # the linear term ingredients come from the system rather than from the
        # constructor arguments, so the precomputed contribution is by
        # construction the one the assembly would have scattered
        dp = planstructurerealjacobian(P, Tv, Ami, Amc, Ljb, Lmean,
            nodesandsigns, sys.invLnm, sys.Gnm, sys.Cnm, sys.wmodesm,
            sys.wmodes2m, layout, layout, Nmodes, Nfreq, backend;
            transposed = transposed)
        # on a host the assembly writes into the matrix that will be factorized
        nzval = transposed ? tobackend(backend, zeros(Tv, nnz(P))) : nonzeros(P)
        return P, dp, nzval
    end

    S = if couplingmodes === :none
        Int[]
    elseif couplingmodes === :all
        collect(1:Nmodes)
    elseif couplingmodes isa AbstractVector{<:Integer}
        sort!(unique(Vector{Int}(couplingmodes)))
    elseif couplingmodes isa AbstractMatrix{Bool}
        size(couplingmodes) == (Nmodes, Nmodes) || throw(DimensionMismatch(
            lazy"a `couplingmodes` mask is $(size(couplingmodes)) but there are $(Nmodes) modes."))
        couplingmodes
    else
        throw(ArgumentError(
            lazy"`couplingmodes` = $(couplingmodes) must be `:none`, `:all`, a vector of mode indices, or an Nmodes x Nmodes Bool mask."))
    end

    P, dp, nzval = build(S)
    return ModeCouplingPreconditioner(P, sys, FactorizationCache(),
        factorization, selectfactorization(factorization, P, S, layout, Nmodes),
        build, Int(Nmodes), S, 0, 0, dp, nzval)
end

# The factorization to use for a freshly built `P`. The mode block diagonal is
# a uniform batch of one block per mode whenever every mode contributes the
# same number of real slots, and a batched direct solver then analyzes one
# block instead of the whole block diagonal and factorizes the batch in
# parallel. Only the empty coupling set can be such a batch: retaining any
# coupling couples modes by construction, and escalation goes straight to the
# full Jacobian, so the two cases skip `blockdiagonallayout` rather than pay a
# walk over every stored entry to be told what is already known.
function selectfactorization(factorization::Factorization, P, S,
    layout::ModeLayout, Nmodes::Integer)
    (S isa Vector{Int} && isempty(S)) || return factorization
    isnothing(factorization.batched) && return factorization
    # discovering the block layout is a walk over the stored entries of a host
    # matrix, and on a backend the structure is not one. The batched path is
    # off by default anyway, having measured slower than handing cuDSS the
    # whole block diagonal.
    P isa SparseMatrixCSC || return factorization
    return batchedfactorization(P, modeslotindex(layout), Nmodes,
        factorization, factorization)
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
    pc.couplingmodes isa Vector{Int} &&
        length(pc.couplingmodes) >= pc.Nmodes && return false
    pc.couplingmodes isa AbstractMatrix{Bool} && all(pc.couplingmodes) && return false
    S = collect(1:pc.Nmodes)
    pc.P, pc.deviceplan, pc.nzval = pc.build(S)
    # the full Jacobian couples every mode, so whatever batch structure the
    # block diagonal had is gone and the caller's own factorization applies
    pc.factorization = pc.basefactorization
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
    # the Fourier coefficients are already on the backend, the assembly runs
    # there, and it writes the order the factorization stores, so on a device
    # nothing crosses to the host: no copy of phimatrix, no serial assembly
    # loop, and no permutation of the values on their way back down.
    _updatecosphimatrix!(pc.sys)
    assemblerealjacobian!(pc.nzval, pc.deviceplan, pc.sys.phimatrix)
    A = pc.P isa SparseMatrixCSC ? pc.P :
        DeviceValuedSparseMatrix(pc.P, pc.nzval)
    tryfactorize!(pc.cache, pc.factorization, A)
    pc.updates += 1
    return pc
end

function applypreconditioner!(z::AbstractVector,
    pc::ModeCouplingPreconditioner, r::AbstractVector)
    return trysolve!(z, pc.cache.factorization, r)
end
