
"""
    cosphibandwidths(sys::HBSystem, Amatrixindices::Matrix,
        Amatrixmodes::AbstractMatrix, Nfreq::Integer = 0,
        Nbranches::Integer = 0; tol = 1e-2, budget = 0.25)

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
    # The frequency stride comes from the coefficient array's own shape,
    # `(frequency dims..., branch)`, which is multidimensional for several
    # tones: everything but the last dimension. This is the linear index the
    # assembly kernel uses, `phimatrix[abs(ind) + nfreq*(b-1)]`.
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

    # The smallest per tone bandwidth whose discarded tail is within `tol`.
    # This sets the shape of the band: a tone which contributes little phase
    # amplitude has a tail which is negligible already at zero offset and
    # collapses to a diagonal, the anisotropy the Jacobi-Anger expansion
    # predicts.
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

    # `tol` alone chooses the bandwidth which makes the preconditioner
    # accurate, and under strong drive that is the whole Jacobian: the
    # coefficients decay, but from a diagonal which is not dominant, so the
    # tail stays above any useful tolerance out to the edge of the grid. A
    # preconditioner does not need to be accurate but to be worth its fill,
    # so the size is then capped by a budget on the fraction of coupling
    # blocks retained, trimming the widest tone first to preserve the shape.
    keptfraction(q) = count(modebandmask(Amatrixmodes,
        NTuple{Ntones,Int}(q))) / length(Amatrixmodes)
    while keptfraction(p) > budget && any(>(0), p)
        k = argmax(p)
        p[k] -= 1
    end

    # Snap each bandwidth down to an offset the grid realizes. A tone whose
    # harmonics are all odd realizes only even offsets, so a bandwidth of one
    # keeps exactly what zero keeps; reporting one would claim a coupling the
    # mask does not contain and make one escalation step a no-op.
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
    ClusterProbe

The state of a [`Clusters`](@ref) request: the block diagonal preconditioner
the probe solves with, the probe vector and its per mode norms, the probed
strength matrix `W`, the couplings taken so far and the current cluster
mask, and the flag [`stalled!`](@ref) raises to remeasure. See
[`probecouplings!`](@ref) and [`spectralclusters`](@ref).
"""
mutable struct ClusterProbe
    const bd
    const r
    const slot
    const rnorm::Vector{Float64}
    const W::Matrix{Float64}
    const v
    const y
    const zz
    const w
    const pairs::Set{Tuple{Int,Int}}   # the couplings taken into clusters so far
    mask::Matrix{Bool}
    reprobe::Bool                      # set by `stalled!`
    probes::Int
    radius::Float64
end

@kernel function modenormkernel!(w, @Const(zz), @Const(slot))
    i = @index(Global)
    @inbounds Atomix.@atomic w[slot[i]] += abs2(zz[i])
end

"""
    probecouplings!(pr::ClusterProbe, sys::HBSystem, Nmodes::Integer)

Measure the block norms of the block Jacobi iteration matrix at the current
point, one Jacobian product and one block diagonal solve per mode: the
column `j` of `pr.W` holds, for every mode, the norm of what the block
diagonal solve of the product of the Jacobian with a random vector on mode
`j`'s slots leaves on that mode, relative to the input. This is the
coupling strength the cluster rule sorts by, and it includes the gain of
the receiving mode's block, which a coefficient of `cos(phi(t))` alone
does not: measured on a two tone line, a coefficient proxy ranked inert
near-dc modes first and the probe the difference ladder.
"""
function probecouplings!(pr::ClusterProbe, sys::HBSystem, Nmodes::Integer)
    refactorize!(pr.bd)
    backend = sys.nonlineartermplan.backend
    accumulate! = modenormkernel!(backend, 256)
    for j in 1:Nmodes
        pr.v .= pr.r .* (pr.slot .== j) ./ pr.rnorm[j]
        jacobianvectorproduct!(pr.y, sys, pr.v)
        applypreconditioner!(pr.zz, pr.bd, pr.y)
        fill!(pr.w, 0)
        accumulate!(pr.w, pr.zz, pr.slot; ndrange = length(pr.zz))
        KernelAbstractions.synchronize(backend)
        col = sqrt.(tohost(pr.w))
        col[j] = 0
        pr.W[:, j] .= col
    end
    pr.probes += 1
    return pr.W
end

"""
    spectralclusters(W::AbstractMatrix)

Cluster the modes by coupling strength: starting from every mode alone,
merge the two clusters of the strongest coupling not yet inside a cluster,
in decreasing strength `W[i,j] + W[j,i]`, until the couplings left between
clusters have block Jacobi spectral radius below one. Returns the cluster
index of every mode, the radius before and after, and the couplings taken.

The value one is not a threshold to tune: for a nonnegative comparison
matrix of block norms, a spectral radius below one is the condition under
which the block Jacobi iteration on the omitted couplings contracts, so
the rule keeps exactly enough coupling inside the clusters for what is
left outside to be correctable. It finds collective chains that no pairwise
threshold does: on a two tone line every single coupling of the difference
ladder is weak but the ladder as a whole is not, and the rule closes it.
"""
function spectralclusters(W::AbstractMatrix)
    N = size(W, 1)
    pairs = [(W[i, j] + W[j, i], i, j) for i in 1:N for j in i+1:N
        if W[i, j] > 0 || W[j, i] > 0]
    sort!(pairs; by = p -> -p[1])
    parent = collect(1:N)
    function find(i)
        while parent[i] != i
            parent[i] = parent[parent[i]]
            i = parent[i]
        end
        return i
    end
    nzs = [(i, j) for i in 1:N for j in 1:N if i != j && W[i, j] > 0]
    # the spectral radius of the couplings between clusters, by power
    # iteration on the nonnegative matrix
    function radius()
        v = ones(N)
        nrm = 0.0
        for _ in 1:80
            u = zeros(N)
            for (i, j) in nzs
                find(i) == find(j) || (u[i] += W[i, j]*v[j])
            end
            nrm = norm(u)
            nrm == 0 && return 0.0
            v = u ./ nrm
        end
        return nrm
    end
    k = 0
    r0 = radius()
    r = r0
    chunk = max(1, length(pairs) ÷ 200)
    while r >= 1 && k < length(pairs)
        for _ in 1:chunk
            k < length(pairs) || break
            k += 1
            a = find(pairs[k][2]); b = find(pairs[k][3])
            a != b && (parent[a] = b)
        end
        r = radius()
    end
    roots = [find(i) for i in 1:N]
    ids = Dict(x => m for (m, x) in enumerate(unique(roots)))
    return [ids[x] for x in roots], r0, r, [(p[2], p[3]) for p in pairs[1:k]]
end

# the coupling mask of a clustering: every pair inside a cluster, and the
# diagonal
function clustermask(cids::AbstractVector{<:Integer})
    N = length(cids)
    return [cids[i] == cids[j] for i in 1:N, j in 1:N]
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

The default is [`BlockDiagonal`](@ref), the mode block diagonal. Its
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
- `coupling`: the coupling set currently factorized, a
    [`BlockDiagonal`](@ref), [`FullJacobian`](@ref), [`HarmonicBand`](@ref),
    [`CoupledModes`](@ref) or [`CouplingMask`](@ref); a [`MeasuredBand`](@ref)
    or [`Clusters`](@ref) request is held as the band or mask it has grown
    to so far.
- `updates`: the number of times the factorization has been rebuilt.
- `escalations`: the number of times the coupling set has been grown.
"""
mutable struct ModeCouplingPreconditioner{TS,TB} <: AbstractPreconditioner
    # untyped, because on a backend `P` is a `DeviceSparsePattern` rather
    # than a host sparse matrix
    P
    sys::TS      # replaced when the system is rebound to new values
    const cache::FactorizationCache
    const basefactorization::AbstractFactorization
    factorization::AbstractFactorization
    # rebuilds `(P, plan, nzval)` for a coupling set, closing over the plan
    # ingredients so the set can be grown after construction
    const build::TB
    const Nmodes::Int
    # the harmonic offset of every mode pair, `modes[m1] .- modes[m2]`, kept so
    # a bandwidth restriction can be rebuilt and stepped after construction.
    # `nothing` when the caller did not supply it, which disables `:band`.
    const Amatrixmodes
    # ingredients of the `MeasuredBand` bandwidth measurement, `nothing`
    # when the caller did not ask for it
    const autoindices
    const autotol
    const autobudget
    const autoNfreq::Int
    const autoNbranches::Int
    # the coupling set currently factorized: a `BlockDiagonal`,
    # `FullJacobian`, `HarmonicBand`, `CoupledModes` or `CouplingMask`
    coupling::AbstractModeCoupling
    updates::Int
    escalations::Int
    # the assembly plan and the values it writes, untyped because they are
    # rebuilt on escalation with a possibly different index type; on a host
    # `nzval` aliases the stored values of `P`
    deviceplan
    nzval
    # the state of a `Clusters` request, `nothing` otherwise
    const clusterprobe
end

"""
    ModeCouplingPreconditioner(sys::HBSystem, Amatrixindices::Matrix,
        Amatrixconjindices::Matrix, Ljb::SparseVector, Lmean,
        Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
        Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
        Cnm::SparseMatrixCSC, layout::ModeLayout;
        spec::AbstractModeCoupling = BlockDiagonal(), precision = nothing,
        Amatrixmodes = nothing)

Build a [`ModeCouplingPreconditioner`](@ref) from the same ingredients
[`planstructurerealjacobian`](@ref) takes. `spec` is the member of the mode
coupling family to build ([`BlockDiagonal`](@ref), [`FullJacobian`](@ref),
[`HarmonicBand`](@ref), [`MeasuredBand`](@ref), [`Clusters`](@ref),
[`CoupledModes`](@ref) or [`CouplingMask`](@ref)), carrying the
factorization it is built with, the backend's default (KLU on the host,
cuDSS on a device) when it carries none. A [`BlockFactorization`](@ref)
eliminates the circuit graph with dense blocks over the clusters of the
coupling set instead of factorizing a sparse matrix. `precision` is the
floating point type of the factorization, `nothing` for that of the
system. `Amatrixmodes` is the harmonic offset of every mode pair, needed by
`HarmonicBand` and `MeasuredBand`.
"""
function ModeCouplingPreconditioner(sys, Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, Ljb::SparseVector, Lmean,
    Rbnm::SparseMatrixCSC, Nmodes::Integer, Nbranches::Integer,
    Nfreq::Integer, invLnm::SparseMatrixCSC, Gnm::SparseMatrixCSC,
    Cnm::SparseMatrixCSC, layout::ModeLayout;
    spec::AbstractModeCoupling = BlockDiagonal(),
    precision::Union{Nothing,Type{<:AbstractFloat}} = nothing,
    Amatrixmodes = nothing)

    backend = sys.nonlineartermplan.backend
    spec isa Automatic && (spec = resolveautomatic(sys, Rbnm, Nmodes,
        Nbranches, layout, Amatrixmodes, backend))
    factorization = something(spec.factorization, defaultfactorization(backend))
    # A measured band remeasures the per tone bandwidth at every point and
    # rebuilds when it grows, starting from the block diagonal so that the
    # first factorization is the cheap one.
    autotol = spec isa MeasuredBand ? spec.tol : nothing
    autobudget = spec isa MeasuredBand ? spec.budget : nothing
    if spec isa MeasuredBand || spec isa HarmonicBand
        isnothing(Amatrixmodes) && throw(ArgumentError(
            "a `HarmonicBand` or `MeasuredBand` needs the mode offsets; pass `Amatrixmodes`."))
    end
    # `Clusters` starts from the block diagonal (an empty mask) and grows
    # it by probing; the probe solves with its own block diagonal, sparse
    # whatever the base factorization is
    clusterprobe = if spec isa Clusters
        bdfactorization = factorization isa BlockFactorization ?
            singletonfactorization(factorization, backend) : factorization
        bd = ModeCouplingPreconditioner(sys, Amatrixindices,
            Amatrixconjindices, Ljb, Lmean, Rbnm, Nmodes, Nbranches, Nfreq,
            invLnm, Gnm, Cnm, layout; spec = BlockDiagonal(bdfactorization),
            precision = precision)
        slot = modeslotindex(layout)
        n = layout.rdim
        # a fixed pseudo-random probe vector, so the clusters are
        # reproducible, without a random number dependency
        rh = [2*(hash(i) % 1_000_003)/1_000_003 - 1 for i in 1:n]
        rnorm = zeros(Nmodes)
        for i in 1:n
            rnorm[slot[i]] += abs2(rh[i])
        end
        rnorm .= sqrt.(rnorm)
        dV = x -> tobackend(backend, x)
        ClusterProbe(bd, dV(rh), tobackend(backend, Int32.(slot)), rnorm,
            zeros(Nmodes, Nmodes), dV(zeros(n)), dV(zeros(n)), dV(zeros(n)),
            dV(zeros(Nmodes)), Set{Tuple{Int,Int}}(),
            Matrix{Bool}(I, Nmodes, Nmodes), false, 0, NaN)
    else
        nothing
    end
    # what is factorized first
    coupling = if spec isa MeasuredBand
        HarmonicBand(ntuple(_ -> 0, length(first(Amatrixmodes))), factorization)
    elseif spec isa Clusters
        CouplingMask(clusterprobe.mask, factorization)
    elseif spec isa CoupledModes
        all(k -> 1 <= k <= Nmodes, spec.indices) || throw(ArgumentError(
            lazy"the coupled mode indices $(spec.indices) are outside 1:$(Nmodes)."))
        withfactorization(spec, factorization)
    elseif spec isa CouplingMask
        size(spec.mask) == (Nmodes, Nmodes) || throw(DimensionMismatch(
            lazy"a coupling mask is $(size(spec.mask)) but there are $(Nmodes) modes."))
        withfactorization(spec, factorization)
    else
        withfactorization(spec, factorization)
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
    # everything with a value in it is read from the system, so that a
    # system rebound to new values is what a rebuild sees; the constructor
    # arguments are the same objects at construction
    function build(S::AbstractModeCoupling, sys = sys)
        Ljb, Lmean = sys.Ljb, sys.Lmean
        invLnm, Gnm, Cnm = sys.invLnm, sys.Gnm, sys.Cnm
        keep = couplingmask(S, Nmodes, Amatrixmodes)
        backend = sys.nonlineartermplan.backend
        Tv = isnothing(precision) ?
            real(promote_type(typeof(Lmean),
                isempty(Ljb.nzval) ? typeof(Lmean) : real(eltype(Ljb)))) :
            precision
        # a block factorization eliminates the circuit graph with dense
        # blocks over the clusters of the coupling set, and the modes left
        # single by the block diagonal, which is this same preconditioner
        # with the empty coupling set and the backend's sparse factorization
        if factorization isa BlockFactorization
            singletons = if isempty(singletonmodes(keep))
                nothing
            else
                ModeCouplingPreconditioner(sys, Amatrixindices,
                    Amatrixconjindices, Ljb, Lmean, Rbnm, Nmodes, Nbranches,
                    Nfreq, invLnm, Gnm, Cnm, layout;
                    spec = BlockDiagonal(singletonfactorization(factorization,
                        backend)), precision = precision)
            end
            Tb = something(factorization.precision, Tv)
            return blockstructure(Tb, sys, Amatrixindices, Amatrixconjindices,
                keep, Rbnm, Nmodes, Nbranches, Nfreq, layout, singletons),
                nothing, nothing
        end
        Ami = restrictmodecoupling(Amatrixindices, keep)
        Amc = restrictmodecoupling(Amatrixconjindices, keep)
        transposed = !(backend isa CPU)

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
        # on a host the assembly writes into the matrix which is factorized
        nzval = transposed ? tobackend(backend, zeros(Tv, nnz(P))) : nonzeros(P)
        return P, dp, nzval
    end

    P, dp, nzval = build(coupling)
    return ModeCouplingPreconditioner(P, sys, FactorizationCache(),
        factorization, selectfactorization(factorization, P, coupling, layout,
        Nmodes), build, Int(Nmodes), Amatrixmodes,
        isnothing(autotol) ? nothing : Amatrixindices, autotol, autobudget,
        Int(Nfreq), Int(Nbranches), coupling, 0, 0, dp, nzval, clusterprobe)
end

# the backend's default sparse factorization: KLU on the host, cuDSS on a
# device, where a host factorization could not be applied to the device
# vectors of the Krylov iteration
defaultfactorization(backend) = backend isa CPU ? KLUfactorization() :
    CUDSSFactorization()

"""
    resolveautomatic(sys, Rbnm, Nmodes, Nbranches, layout, Amatrixmodes,
        backend; budget = freememory(backend) ÷ 2)

The member of the mode coupling family an [`Automatic`](@ref) request
stands for on this problem: [`BlockDiagonal`](@ref) for one tone;
otherwise [`FullJacobian`](@ref) with single precision block factors when
[`blockfactorbytes`](@ref) is within `budget`, and [`MeasuredBand`](@ref)
when it is not.
"""
function resolveautomatic(sys, Rbnm::SparseMatrixCSC, Nmodes::Integer,
    Nbranches::Integer, layout::ModeLayout, Amatrixmodes, backend;
    budget::Integer = freememory(backend) ÷ 2)
    ntones = isnothing(Amatrixmodes) ? 1 : length(first(Amatrixmodes))
    ntones == 1 && return BlockDiagonal()
    nnodes = layout.dim ÷ Nmodes
    nodesandsigns = branchnodesandsigns(Rbnm, Nmodes, Nbranches)
    pairptr, pairrow, _, _ = junctionpairtable(Int32, Float32, sys.Ljb,
        nodesandsigns, nnodes)
    adj = circuitnodegraph(pairptr, pairrow, sys.invLnm, sys.Gnm, sys.Cnm,
        Nmodes, nnodes)
    order = klunodeorder(adj)
    bytes = blockfactorbytes(Float32, modecouplingmask(Nmodes, 1:Nmodes), adj,
        order, Nmodes, layout)
    bytes <= budget && return FullJacobian(BlockFactorization(;
        precision = Float32))
    isnothing(Amatrixmodes) && return BlockDiagonal()
    return MeasuredBand()
end

"""
    couplingmask(S::AbstractModeCoupling, Nmodes::Integer, Amatrixmodes)

The `Nmodes` x `Nmodes` mask of the mode couplings the set `S` retains.
"""
couplingmask(::BlockDiagonal, Nmodes::Integer, Amatrixmodes) =
    modecouplingmask(Nmodes, Int[])
couplingmask(::FullJacobian, Nmodes::Integer, Amatrixmodes) =
    modecouplingmask(Nmodes, 1:Nmodes)
couplingmask(S::HarmonicBand, Nmodes::Integer, Amatrixmodes) =
    modebandmask(Amatrixmodes, S.p)
couplingmask(S::CoupledModes, Nmodes::Integer, Amatrixmodes) =
    modecouplingmask(Nmodes, S.indices)
couplingmask(S::CouplingMask, Nmodes::Integer, Amatrixmodes) = S.mask

# The factorization to use for a freshly built `P`. The mode block diagonal is
# a uniform batch of one block per mode whenever every mode contributes the
# same number of real slots, and a batched direct solver then analyzes one
# block instead of the whole block diagonal and factorizes the batch in
# parallel. Only the empty coupling set can be such a batch: retaining any
# coupling couples modes by construction, and escalation goes straight to the
# full Jacobian, so the two cases skip `blockdiagonallayout` rather than pay a
# walk over every stored entry to be told what is already known.
function selectfactorization(factorization::AbstractFactorization, P, S,
    layout::ModeLayout, Nmodes::Integer)
    S isa BlockDiagonal || return factorization
    (factorization isa CUDSSFactorization && factorization.batched) ||
        return factorization
    # discovering the block layout is a walk over the stored entries of a
    # host matrix, and on a backend the structure is not one; the batched
    # path is off by default anyway, being slower than handing cuDSS the
    # whole block diagonal
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
    isexactpreconditioner(pc) && return false
    # A band escalates by one offset per tone at a time, and becomes the
    # full mode set once it covers the grid. Every other coupling set jumps
    # straight to the full Jacobian, since nothing identifies which modes
    # carried the deficiency.
    f = pc.coupling.factorization
    S = if pc.coupling isa HarmonicBand
        p = pc.coupling.p
        next = p isa Integer ? p + 1 : map(x -> x + 1, p)
        all(modebandmask(pc.Amatrixmodes, next)) ? FullJacobian(f) :
            HarmonicBand(next, f)
    else
        FullJacobian(f)
    end
    pc.P, pc.deviceplan, pc.nzval = pc.build(S, pc.sys)
    # any retained coupling couples modes, so whatever batch structure the
    # block diagonal had is gone and the caller's own factorization applies
    pc.factorization = pc.basefactorization
    pc.coupling = S
    pc.escalations += 1
    # the cached symbolic factorization is of the previous structure
    pc.cache.factorization = nothing
    return true
end

# Whether the coupling set is the full one: every mode retained, or a mask
# with no zero. A `:band` never reports exact, because escalation replaces
# a band which covers the grid by the full set.
function isexactpreconditioner(pc::ModeCouplingPreconditioner)
    S = pc.coupling
    S isa FullJacobian && return true
    S isa CoupledModes && length(S.indices) >= pc.Nmodes && return true
    S isa CouplingMask && all(S.mask) && return true
    return false
end

function updatepreconditioner!(pc::ModeCouplingPreconditioner,
    x::AbstractVector)
    setpoint!(pc.sys, x)
    # a measured band remeasures the bandwidth here, where the Fourier coefficients
    # are computed anyway, and rebuilds only when it has grown. It is never
    # reduced: that would cost a symbolic refactorization to save fill
    # already paid for, and the drive generally widens the coupling as the
    # iteration proceeds.
    if !isnothing(pc.autotol) && pc.coupling isa HarmonicBand
        want = cosphibandwidths(pc.sys, pc.autoindices, pc.Amatrixmodes;
            tol = pc.autotol, budget = pc.autobudget)
        have = pc.coupling.p isa NTuple ? pc.coupling.p :
            ntuple(_ -> 0, length(want))
        grown = map(max, want, have)
        if grown != have
            pc.coupling = HarmonicBand(grown, pc.coupling.factorization)
            pc.P, pc.deviceplan, pc.nzval = pc.build(pc.coupling, pc.sys)
            pc.factorization = pc.basefactorization
            pc.cache.factorization = nothing
        end
    end
    # `Clusters` probes at the first point and whenever the driver reports
    # a slow linear solve (`stalled!`), the sign that the coupling has
    # outgrown the clusters; the mask only grows, so the structure is
    # rebuilt only when a probe adds to it
    pr = pc.clusterprobe
    if !isnothing(pr) && pc.coupling isa CouplingMask
        if pr.probes == 0 || pr.reprobe
            pr.reprobe = false
            probecouplings!(pr, pc.sys, pc.Nmodes)
            _, _, pr.radius, taken = spectralclusters(pr.W)
            # the couplings taken accumulate across probes and the clusters
            # are their components: monotone, and two probes which agree on
            # the families keep them apart
            union!(pr.pairs, taken)
            edges = Matrix{Bool}(I, pc.Nmodes, pc.Nmodes)
            for (i, j) in pr.pairs
                edges[i, j] = edges[j, i] = true
            end
            grown = copy(edges)
            for g in modeclusters(edges), i in g, j in g
                grown[i, j] = true
            end
            if grown != pr.mask
                pr.mask = grown
                pc.coupling = CouplingMask(grown, pc.coupling.factorization)
                pc.P, pc.deviceplan, pc.nzval = pc.build(pc.coupling, pc.sys)
                pc.factorization = pc.basefactorization
                pc.cache.factorization = nothing
            end
        end
    end
    refactorize!(pc)
    return pc
end

function stalled!(pc::ModeCouplingPreconditioner)
    isnothing(pc.clusterprobe) || (pc.clusterprobe.reprobe = true)
    return pc
end

"""
    refactorize!(pc::ModeCouplingPreconditioner)

Assemble the restricted Jacobian at the system's current point and
factorize it: what [`updatepreconditioner!`](@ref) does after setting the
point and choosing the coupling set.
"""
function refactorize!(pc::ModeCouplingPreconditioner)
    # the Fourier coefficients are on the backend, the assembly runs there
    # and writes the order the factorization stores, so on a device nothing
    # crosses to the host
    _updatecosphimatrix!(pc.sys)
    A = if pc.P isa BlockStructure
        # the block factorization assembles its own blocks from the same
        # coefficients
        BlockJacobian(pc.P, pc.sys.phimatrix)
    else
        assemblerealjacobian!(pc.nzval, pc.deviceplan, pc.sys.phimatrix)
        pc.P isa SparseMatrixCSC ? pc.P :
            DeviceValuedSparseMatrix(pc.P, pc.nzval)
    end
    tryfactorize!(pc.cache, pc.factorization, A)
    pc.updates += 1
    return pc
end

function applypreconditioner!(z::AbstractVector,
    pc::ModeCouplingPreconditioner, r::AbstractVector)
    return trysolve!(z, pc.cache.factorization, r)
end

"""
    rebind!(pc::ModeCouplingPreconditioner, sys::HBSystem)

The same preconditioner over a system rebound to new component values: the
constant linear contribution of its assembly plan and the junction
coefficients are refreshed from the system, and the structure, the
coupling set and the factorization's symbolic analysis are kept. The next
[`updatepreconditioner!`](@ref) refactorizes the numbers.
"""
function rebind!(pc::ModeCouplingPreconditioner, sys::HBSystem)
    pc.sys = sys
    isnothing(pc.clusterprobe) || rebind!(pc.clusterprobe.bd, sys)
    if pc.P isa BlockStructure
        refreshvalues!(pc.P, sys)
    else
        refreshvalues!(pc.deviceplan, sys.invLnm, sys.Gnm, sys.Cnm,
            sys.wmodesm, sys.wmodes2m, sys.Ljb, sys.Lmean)
    end
    return pc
end
