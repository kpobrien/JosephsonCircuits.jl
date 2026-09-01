
"""
    cosphibandwidths(sys::HBSystem, Amatrixindices::Matrix,
        Amatrixmodes::AbstractMatrix, Nfreq::Integer, Nbranches::Integer;
        tol = 1e-2)

The per-tone harmonic bandwidth of the Josephson coupling at the current point,
as the tuple of the largest offset in each tone which carries a coefficient of
`cos(phi(t))` above `tol` relative to the largest.

This is the bandwidth [`modebandmask`](@ref) should be given, measured rather
than guessed, and it costs nothing beyond the Fourier coefficients
[`updatepreconditioner!`](@ref) already computes: the assembly reads
`phimatrix[abs(ind) + Nfreq*(b-1)]` for `ind = Amatrixindices[m1,m2]`, and this
reads the same entries and records how far from the diagonal they stay
significant.

The measurement is per tone, and that anisotropy is the point. Jacobi-Anger
gives the coefficient at a multi-tone offset `m` as a product over the tones,
`chat_m ~ prod_k J_{m_k}(delta_k)`, with `delta_k` the phase amplitude tone `k`
contributes. Its support is therefore a *rectangle* whose sides are set by each
`delta_k` separately, not a ball: a strongly pumped tone earns a wide bandwidth
and a weak one collapses to zero, because `J_n(delta)` falls off super
geometrically once `n` exceeds `delta`. Keeping a ball instead spends fill on
offsets in the weak tones which carry nothing, and the difference is large. On a
two tone chain with a 5 percent second drive and 74 modes, the measured
rectangle `(2,0)` converges where the ball of any radius does not, at four
percent of the full Jacobian's stored entries.

A tone whose grid is short enough that offsets alias will report a large
bandwidth. That is not a failure: the aliased coupling is really there, and a
bandwidth which saturates that tone's grid costs nothing extra because the grid
had no room to be truncated in the first place.
"""
function cosphibandwidths(sys, Amatrixindices::Matrix,
    Amatrixmodes::AbstractMatrix, Nfreq::Integer = 0, Nbranches::Integer = 0;
    tol = 1e-2, budget = 0.25)

    _updatecosphimatrix!(sys)
    cm = tohost(sys.phimatrix)
    Ntones = length(first(Amatrixmodes))
    # The coefficient array's own shape, not the caller's Nfreq: it is
    # `(frequency dims..., branch)` and multi-tone leaves it multidimensional,
    # so the frequency stride is everything but the last dimension. This is the
    # linear index the assembly kernel itself uses,
    # `phimatrix[abs(ind) + nfreq*(b-1)]`.
    nb = size(cm)[end]
    nfreq = length(cm) ÷ nb

    mag(ind) = begin
        a = abs(ind)
        (1 <= a <= nfreq) || return 0.0
        m = 0.0
        for b in 1:nb
            m = max(m, abs(cm[a + nfreq*(b - 1)]))
        end
        m
    end

    # marginal weight of each per-tone offset magnitude, and the total
    maxoff = zeros(Int, Ntones)
    for o in Amatrixmodes, k in 1:Ntones
        maxoff[k] = max(maxoff[k], abs(o[k]))
    end
    w = [zeros(Float64, maxoff[k] + 1) for k in 1:Ntones]
    total = 0.0
    for i in axes(Amatrixindices, 1), j in axes(Amatrixindices, 2)
        ind = Int(Amatrixindices[i, j])
        ind == 0 && continue
        m = mag(ind)
        m == 0.0 && continue
        total += m
        o = Amatrixmodes[i, j]
        for k in 1:Ntones
            w[k][abs(o[k]) + 1] += m
        end
    end
    total > 0 || return ntuple(_ -> 0, Ntones)

    # smallest per-tone bandwidth whose discarded marginal tail is within `tol`.
    # This sets the *shape*: a tone which contributes little phase amplitude has
    # a tail that is already negligible at zero and collapses to a diagonal,
    # which is the anisotropy the Jacobi-Anger product predicts.
    p = zeros(Int, Ntones)
    for k in 1:Ntones
        pk = maxoff[k]
        for cand in 0:maxoff[k]
            tail = sum(@view w[k][cand+2:end])
            if tail <= tol*total
                pk = cand
                break
            end
        end
        p[k] = pk
    end

    # `tol` alone chooses the bandwidth which makes the preconditioner accurate,
    # and under strong drive that is the whole Jacobian: the coefficients decay,
    # but from a diagonal which is no longer dominant, so the tail stays above
    # any useful tolerance out to the edge of the grid. A preconditioner does
    # not need to be accurate, it needs to be worth its fill, so the *size* is
    # then capped by a budget on the fraction of coupling blocks retained. The
    # widest tone is trimmed first, which preserves the shape the tail chose.
    keptfraction(q) = count(modebandmask(Amatrixmodes,
        NTuple{Ntones,Int}(q))) / length(Amatrixmodes)
    while keptfraction(p) > budget && any(>(0), p)
        k = argmax(p)
        p[k] -= 1
    end

    # Snap each bandwidth down onto an offset the grid actually realizes. A
    # tone whose harmonics are all odd realizes only even offsets, so a
    # bandwidth of one there keeps exactly what zero keeps; reporting one would
    # claim a coupling the mask does not contain, and would make the escalation
    # step of `escalatepreconditioner!` a no-op for a whole increment.
    realized = [sort!(unique(abs(o[k]) for o in Amatrixmodes)) for k in 1:Ntones]
    for k in 1:Ntones
        idx = searchsortedlast(realized[k], p[k])
        p[k] = idx >= 1 ? realized[k][idx] : 0
    end
    return NTuple{Ntones,Int}(p)
end

"""
    modebandmask(Amatrixmodes::AbstractMatrix, p)

Return the `Nmodes` x `Nmodes` boolean matrix of mode coupling blocks retained
by a *bandwidth* restriction on the harmonic offset: the block coupling column
mode `m2` into row mode `m1` is kept when the offset between them,
`Amatrixmodes[m1, m2] = modes[m1] .- modes[m2]`, is within `p`.

This is the restriction the structure of the Jacobian asks for. Multiplication
by `cos(phi(t))` is a convolution in the harmonic index, so the nonlinear part
of the Jacobian is block Toeplitz in the offset `m1 - m2`: every block shares
the junction incidence sparsity and they differ only by the Fourier coefficient
of `cos(phi(t))` at that offset. Those coefficients are Bessel-like in the
junction phase amplitude and fall off quickly once the offset exceeds it, so
truncating by offset keeps the large blocks and drops the small ones.

Truncating by *column*, as [`modecouplingmask`](@ref) does, cuts across the
Toeplitz structure instead: a retained column keeps one large block and a whole
column of small ones, and drops large blocks elsewhere. Measured at equal fill
the difference is not marginal. On an eight mode, eight junction chain driven to
`max|phi| = 1.9` rad, a bandwidth of one converges the linear solves of a Newton
path in 118 GMRES iterations while a two column selection with the same number
of stored nonzeros fails to converge in 1051.

`p` may be

- an `Integer`: the number of offset *shells* retained beyond the diagonal,
  where a shell is one distinct value of the total intermodulation order
  `sum(abs, offset)` that the retained mode set actually realizes, or
- an `NTuple{N,Integer}`: an absolute per-tone bound, `abs.(offset) .<= p`, for
  a grid whose coupling is anisotropic, eg. a strong pump and a weak second
  tone.

Counting shells rather than bounding the raw offset matters, because the offsets
a mode set realizes are not the integers. With the usual odd-harmonic-only
truncation every mode difference is even, so `sum(abs, offset) <= 1` retains
exactly what `<= 0` does and a raw bound of one would silently be the block
diagonal. A shell count is also the quantity the coupling decays in: the Fourier
coefficients of `cos(phi(t))` live on the harmonic lattice of `phi`, and the
shells are the lattice distances that lattice actually has.

`p = 0` is the mode block diagonal and a `p` large enough to cover the grid is
the full Jacobian, so the bandwidth is a graded ladder between the two rather
than the jump [`escalatepreconditioner!`](@ref) had to make when the only
alternatives were a column set and the whole operator.

# Examples
```jldoctest
julia> modes = [(-1,), (0,), (1,)];

julia> A = [modes[i] .- modes[j] for i in 1:3, j in 1:3];

julia> JosephsonCircuits.modebandmask(A, 0)
3×3 Matrix{Bool}:
 1  0  0
 0  1  0
 0  0  1

julia> JosephsonCircuits.modebandmask(A, 1)
3×3 Matrix{Bool}:
 1  1  0
 1  1  1
 0  1  1
```
"""
function modebandmask(Amatrixmodes::AbstractMatrix, p)
    return [_withinband(Amatrixmodes[i, j], p, _shells(Amatrixmodes, p))
        for i in axes(Amatrixmodes, 1), j in axes(Amatrixmodes, 2)]
end

# The distinct total intermodulation orders the retained mode set realizes,
# sorted. Only needed for an Integer bandwidth.
_shells(Amatrixmodes, ::NTuple) = nothing
function _shells(Amatrixmodes::AbstractMatrix, ::Integer)
    return sort!(unique(sum(abs, Amatrixmodes[i, j])
        for i in axes(Amatrixmodes, 1), j in axes(Amatrixmodes, 2)))
end
_shells(Amatrixmodes, p) = throw(ArgumentError(
    lazy"the bandwidth `p` = $(p) must be an Integer or a tuple of Integers."))

# `p` shells of offset beyond the diagonal
function _withinband(offset, p::Integer, shells)
    p >= 0 || throw(ArgumentError(lazy"the bandwidth `p` = $(p) must be nonnegative."))
    cutoff = shells[min(p + 1, length(shells))]
    return sum(abs, offset) <= cutoff
end
_withinband(offset, p::Integer) = sum(abs, offset) <= p
# per tone bound
function _withinband(offset, p::NTuple{N,Integer}, shells) where N
    length(offset) == N || throw(DimensionMismatch(
        lazy"the offset has $(length(offset)) tones but the bandwidth has $(N)."))
    return all(abs(offset[k]) <= p[k] for k in 1:N)
end

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
- `Amatrixmodes`: the harmonic offset `modes[m1] .- modes[m2]` of every mode
    pair, or `nothing`. Needed by `:band` and by its escalation.
- `couplingmodes`: the modes whose coupling is retained in full, `:band => p`
    for a bandwidth restriction on the harmonic offset (see
    [`modebandmask`](@ref)), or a mask of the retained couplings when one was
    given directly.
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
    # the harmonic offset of every mode pair, `modes[m1] .- modes[m2]`, kept so
    # a bandwidth restriction can be rebuilt and stepped after construction.
    # `nothing` when the caller did not supply it, which disables `:band`.
    const Amatrixmodes
    # ingredients of the `:auto` bandwidth measurement, `nothing` when the
    # caller did not ask for it
    const autoindices
    const autotol
    const autoNfreq::Int
    const autoNbranches::Int
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
- an `AbstractVector{<:Integer}`: exactly these mode indices,
- `:band => p`: every coupling whose harmonic offset is within `p`, which
    requires `Amatrixmodes` and is usually the best of these (see
    [`modebandmask`](@ref)),
- `:auto` or `:auto => tol`: a `:band` whose per-tone bandwidth is measured from
    the Fourier coefficients of `cos(phi(t))` at every point and widened when
    the drive demands it (see [`cosphibandwidths`](@ref)). Starts at the block
    diagonal, so nothing is paid until it is needed.
"""
function ModeCouplingPreconditioner(sys, Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, Ljb::SparseVector, Lmean,
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC, layout::ModeLayout; couplingmodes = :none,
    factorization = KLUfactorization(),
    precision::Union{Nothing,Type{<:AbstractFloat}} = nothing,
    Amatrixmodes = nothing)

    # `:auto` measures the per tone bandwidth at every point and rebuilds when
    # it grows. It starts from the block diagonal, so the first factorization is
    # the cheap one and the structure widens only as the drive demands.
    autotol = if couplingmodes === :auto
        1e-2
    elseif couplingmodes isa Pair && couplingmodes.first === :auto
        couplingmodes.second
    else
        nothing
    end
    if !isnothing(autotol)
        isnothing(Amatrixmodes) && throw(ArgumentError(
            "`couplingmodes = :auto` needs the mode offsets; pass `Amatrixmodes`."))
        couplingmodes = :band => ntuple(_ -> 0, length(first(Amatrixmodes)))
    end

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
        keep = if S isa AbstractMatrix{Bool}
            S
        elseif S isa Pair && S.first === :band
            isnothing(Amatrixmodes) && throw(ArgumentError(
                "`couplingmodes = :band => p` needs the mode offsets; pass `Amatrixmodes`."))
            modebandmask(Amatrixmodes, S.second)
        else
            modecouplingmask(Nmodes, S)
        end
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
    elseif couplingmodes isa Pair && couplingmodes.first === :band
        couplingmodes
    elseif couplingmodes isa AbstractVector{<:Integer}
        sort!(unique(Vector{Int}(couplingmodes)))
    elseif couplingmodes isa AbstractMatrix{Bool}
        size(couplingmodes) == (Nmodes, Nmodes) || throw(DimensionMismatch(
            lazy"a `couplingmodes` mask is $(size(couplingmodes)) but there are $(Nmodes) modes."))
        couplingmodes
    else
        throw(ArgumentError(
            lazy"`couplingmodes` = $(couplingmodes) must be `:none`, `:all`, `:band => p`, a vector of mode indices, or an Nmodes x Nmodes Bool mask."))
    end

    P, dp, nzval = build(S)
    return ModeCouplingPreconditioner(P, sys, FactorizationCache(),
        factorization, selectfactorization(factorization, P, S, layout, Nmodes),
        build, Int(Nmodes), Amatrixmodes, isnothing(autotol) ? nothing :
        Amatrixindices, autotol, Int(Nfreq), Int(Nbranches), S, 0, 0, dp,
        nzval)
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
    # a bandwidth escalates by one offset at a time. The jump straight to the
    # full Jacobian below exists because no per mode score identified *which*
    # modes carried the deficiency; a bandwidth is not a choice of modes, so
    # the graded step is available and is what is taken here. It stops when the
    # band covers the grid, at which point it is already the full Jacobian.
    S = if pc.couplingmodes isa Pair && pc.couplingmodes.first === :band
        p = pc.couplingmodes.second
        next = p isa Integer ? p + 1 : map(x -> x + 1, p)
        all(modebandmask(pc.Amatrixmodes, next)) ? collect(1:pc.Nmodes) :
            (:band => next)
    else
        collect(1:pc.Nmodes)
    end
    pc.P, pc.deviceplan, pc.nzval = pc.build(S)
    # any retained coupling couples modes, so whatever batch structure the
    # block diagonal had is gone and the caller's own factorization applies
    pc.factorization = pc.basefactorization
    pc.couplingmodes = S
    pc.escalations += 1
    # the cache holds a symbolic factorization of the previous sparsity
    # structure, which no longer applies
    pc.cache.factorization = nothing
    return true
end

# the coupling set is the full one: every mode retained, or a mask with no
# zero. `:band` never reports exact even when wide, because escalation
# replaces a covering band by the full set before it stops growing.
function isexactpreconditioner(pc::ModeCouplingPreconditioner)
    pc.couplingmodes isa Vector{Int} &&
        length(pc.couplingmodes) >= pc.Nmodes && return true
    pc.couplingmodes isa AbstractMatrix{Bool} && all(pc.couplingmodes) &&
        return true
    pc.couplingmodes === :all && return true
    return false
end

function updatepreconditioner!(pc::ModeCouplingPreconditioner,
    x::AbstractVector)
    setpoint!(pc.sys, x)
    # `:auto` remeasures the bandwidth here, where the Fourier coefficients are
    # about to be computed anyway, and rebuilds only when it has grown. The
    # bandwidth is never reduced: the structure would have to be rebuilt and
    # refactorized symbolically to save fill the solve has already paid for, and
    # the drive generally widens the coupling as the Newton iteration proceeds.
    if !isnothing(pc.autotol)
        want = cosphibandwidths(pc.sys, pc.autoindices, pc.Amatrixmodes;
            tol = pc.autotol)
        have = pc.couplingmodes isa Pair && pc.couplingmodes.second isa NTuple ?
            pc.couplingmodes.second : ntuple(_ -> 0, length(want))
        grown = map(max, want, have)
        if grown != have
            pc.couplingmodes = :band => grown
            pc.P, pc.deviceplan, pc.nzval = pc.build(pc.couplingmodes)
            pc.factorization = pc.basefactorization
            pc.cache.factorization = nothing
        end
    end
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
