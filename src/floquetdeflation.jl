# =========================================================================
# Residual-image A-DEF1: low-rank deflation of the directions the base
# preconditioner misses, carried as *physical* correction vectors.
# =========================================================================

"""
    seeddeflation!(pc, Xnew::AbstractMatrix; source = :external)

Offer the columns of `Xnew` to `pc` as candidate physical correction
vectors, between two linear solves, and return `pc`. The default does
nothing, which is correct for any preconditioner which does not deflate.

This is the injection point for an external candidate source: a
continuation secant, a Newton step, or a physically constructed Floquet
mode. Candidates are *not* filtered here. They are concatenated with the
bank and compressed together against their exact residual images at the
next rebuild, which is what lets several sources rediscover the same
channel without any of them having to know about the others. The seed
marks the active blocks stale, so the next application rebuilds them and
the seeded directions take part in the very next solve; for that reason it
must not be called while a GMRES solve is running, where the
preconditioner has to stay fixed. The harvest banks its candidates through
the internal `_bankcandidates!`, which does not mark anything stale.
"""
seeddeflation!(pc::AbstractPreconditioner, ::AbstractMatrix; kwargs...) = pc

"""
    FloquetState(b::AbstractVector)

The part of a [`FloquetPreconditioner`](@ref) that outlives one solve: the
candidate physical correction vectors `X`, an `n` by `k` matrix allocated
like `b`, and the provenance of each column. A candidate is a direction of
the unknowns; its image under the current Jacobian is recomputed at every
rebuild, so it means the same thing after the base has been rebuilt or
rebound to a new operating point. [`HBReuse`](@ref) carries one across a
cached sweep and replaces it only after a converged solve.
"""
mutable struct FloquetState{TM<:AbstractMatrix}
    X::TM
    source::Vector{Symbol}
end

FloquetState(b::AbstractVector) = FloquetState(similar(b, length(b), 0), Symbol[])
Base.copy(s::FloquetState) = FloquetState(copy(s.X), copy(s.source))

"""
    FloquetPreconditioner(inner::AbstractPreconditioner, jvp!, b::AbstractVector;
        kmax = 20, kharvest = 4, nritz = 0, escalateafter = 1,
        ranktol = eps(T)^(3/8), benefittol = 1e-6, kcandidate = 3*kmax,
        cycleharvest = true, state = FloquetState(b))

Augments a base preconditioner with a low-rank correction for the few global
channels it represents badly, in the residual-image A-DEF1 form.

The hypothesis this is built on is not that the harmonic balance Jacobian
`J` is low rank, which it is not, but that the *defect* of a cheap base
preconditioner,

    D = I - J*inv(P),

has low effective rank: that a mode block diagonal already solves almost
every direction of a multi-tone problem, and what it leaves is a handful of
global spectral-spatial channels — an amplifier's signal-idler pair, a
dominant conversion channel, a weakly damped Floquet mode. Those are the
directions this deflates.

# The residual-image form

Candidate correction vectors `X0` are held in the *physical* state space,
not in a residual or Krylov coordinate system. At each rebuild their exact
images `Y0 = J*X0` are formed and factorized, `Y0 = U*S*V'`. Keeping the
numerically significant rank `r`,

    C = U[:, 1:r],    X = X0*V[:, 1:r]*inv(S[1:r])

so that

    J*X = C    and    C'C = I

hold by construction. With `W = X - inv(P)*C` the preconditioner is

    inv(Pdef) = inv(P) + W*C'

applied as `z = inv(P)*r + W*(C'*r)`: one base solve and two dense mat-vecs
with an `n` by `r` matrix, and no Jacobian product. Because `J*X = C`
exactly,

    J*inv(Pdef)*C = J*(inv(P)*C + W) = J*X = C,

so `range(C)` is an exact eigenspace of the right preconditioned operator at
eigenvalue one and GMRES never has to resolve those directions again.

# Why not the Galerkin pairing

The `:adef1` and `:adef2` forms of [`RecyclingPreconditioner`](@ref) build
their coarse correction from `G = Z'*J*Z` and a truncated pseudo inverse of
it. That pairing is the natural one for a symmetric positive definite
system. This Jacobian is none of those things: it is nonsymmetric,
strongly nonnormal, and near a parametric threshold nearly singular. A
perfectly good correction direction can then have `J*z` large while
`z'*J*z` is nearly zero, which makes `G` ill conditioned for a reason that
has nothing to do with the quality of `z`. Orthonormalizing the image `J*X`
instead of the pairing `X'*J*X` removes that failure mode, needs no
projected inverse at application time, and measures the subspace in the
Euclidean norm GMRES actually minimizes.

# Candidates

Candidates arrive through [`seeddeflation!`](@ref) and are tagged by
source. [`harvest!`](@ref) contributes up to two families from the Arnoldi
factorization of each solve: the smallest singular directions of the
rectangular Hessenberg, robust on a nonnormal operator, and, when `nritz`
is positive, the harmonic Ritz directions nearest zero
([`harmonicritznearzero`](@ref)), which target the near-singular
eigendirections of a device close to threshold. Both are
mapped through the current preconditioner into physical coordinates before
being stored, so that a candidate stays meaningful after the Jacobian
moves. A continuation secant or an externally constructed Floquet mode is
injected the same way.

No source is trusted. The residual-image factorization decides the rank,
which removes candidates that rediscovered the same channel, and each
surviving direction is then tested for whether the base preconditioner
needs help with it at all,

    eta = norm(x - inv(P)*J*x)/norm(x),

with directions below `benefittol` dropped: physically interesting but
numerically easy channels cost nothing to carry and are not carried. For an
amplifier this is the asymmetry between the amplified quadrature, which
needs a large correction, and its deamplified partner, which often does
not.

That test is applied in the eigenbasis of `W'W` and nowhere else, which is
not a detail. The rank-revealing factorization pins down `range(C)` but not
a basis of it, and the basis it happens to return is arbitrary whenever the
retained singular values are close together — which is the normal case
here, not a degenerate one, because equalizing candidates by their image
norms is what makes those singular values close. `eta` is not invariant
under a rotation of the subspace, so measured in an arbitrary basis a
correction of true rank two can present as six columns of comparable `eta`,
each holding a mixture of what the base misses and what it handles, and
nothing is droppable. The eigenvectors of `W'W` diagonalize the correction,
so they are the one basis in which the split is visible; they also order
the directions by how much correction each carries, which is the order the
cap at `kmax` wants.

# Cost

Per application: one base solve and two dense mat-vecs, the same as
`:adef1`. Per rebuild: `k` Jacobian products for the candidate images and
`r` base solves for `inv(P)*C`, against `2k` base solves and `k` products
for `:adef1`. A rebuild happens when the base is rebuilt and, lazily, when
the point moves under a frozen base ([`pointmoved!`](@ref)).

# Arguments

- `kmax = 20`: the largest active rank retained. When more directions
  survive both filters the update is compressed to the `kmax` principal
  directions of `W'W`, that is, those carrying the most correction.
- `kharvest = 4`: smallest singular directions taken per harvest.
- `nritz = 0`: harmonic Ritz directions nearest zero taken per harvest.
  A complex Ritz pair contributes its real and imaginary parts as two real
  candidates spanning the same invariant subspace. Off by default: on a
  128-junction two-tone line the singular directions alone take 1590
  Arnoldi steps (1732 Jacobian products in all), and adding one or two
  Ritz directions per harvest takes 1955 and 1943 (2125 and 2119); over a
  seven-point cached sweep 10541 against 12623. The Ritz directions of a
  strongly nonnormal operator carry large corrections without being the
  directions GMRES stalls on, and the trim by correction strength then
  keeps them over the singular ones.
- `kcandidate = 3*kmax`: the candidate bank's capacity. The active
  directions of the last rebuild are always kept; older candidates are
  dropped first. A rebuild resets the bank to the active set, so the
  capacity only matters between rebuilds, and anything at or above
  `kmax + kharvest + 2*nritz` behaves the same; below that a harvest
  displaces its own newest candidates (measured: an escalation returns on
  the 128-junction line at `kcandidate = kmax`).
- `ranktol = eps(T)^(3/8)`: relative singular value threshold for the
  residual-image rank. The candidates' images are equalized to unit norm,
  so a small singular value of their block means two images nearly
  parallel, a channel two sources rediscovered, not a small image; the
  threshold is the angle below which they count as one. The rank is read
  off the Gram matrix of the images, whose roundoff is `k*eps` on the
  squared singular values, so the threshold sits above `sqrt(k*eps)` with
  margin (`eps^(3/8)` is `1.4e-6` in double precision).
- `benefittol`: the `eta` below which a direction is judged already handled
  by the base.
- `escalateafter = 1`: as [`RecyclingPreconditioner`](@ref).
- `cycleharvest = true`: harvest at the end of every restart cycle
  ([`harvestcycle!`](@ref)) rather than only from the cycle left in the
  workspace when the solve returns.
- `state`: a [`FloquetState`](@ref) to start from, the candidates of a
  previous solve of a nearby system; it is mutated by the harvests of this
  solve.

The intended base is the mode block diagonal, and everything here is either
a dense level 3 kernel on a small block or a device gemv, with the `k` by
`k` factorizations on the host, so the pairing ports to a GPU for the same
reason that one does. As for [`RecyclingPreconditioner`](@ref), the
three-argument application reads `r` after writing `z`, so the two must not
alias.
"""
mutable struct FloquetPreconditioner{TI,TJ,T<:AbstractFloat,TM<:AbstractMatrix{T},TV<:AbstractVector{T}} <: AbstractPreconditioner
    const inner::TI
    const jvp!::TJ
    # the candidate physical correction vectors, `n` by `k`, with their
    # provenance: the bank several sources feed and the rebuild compresses,
    # and what persists across solves
    state::FloquetState{TM}
    # the active basis: `J*X = C`, `C'C = I`, `W = X - inv(P)*C`. System
    # sized and allocated like the vectors of the system, so on a device
    # backend the correction is a pair of device gemvs.
    X::TM
    C::TM
    W::TM
    coeff::TV
    # the correction strength `eta` of each *active* column at the last
    # rebuild, kept for diagnostics
    strength::Vector{T}
    const kmax::Int
    const kharvest::Int
    const nritz::Int
    const kcandidate::Int
    const escalateafter::Int
    const ranktol::T
    const benefittol::T
    # whether to read every restart cycle rather than only the one left in
    # the workspace when the solve returns
    const cycleharvest::Bool
    escalationrequests::Int
    rebuilds::Int
    # Jacobian products taken for the candidate images
    products::Int
    # whether the active blocks are those of the current Jacobian, the
    # current base and the current candidate bank: cleared by `pointmoved!`
    # and by a seed, restored by `_rebuildfloquet!`
    fresh::Bool
end

function FloquetPreconditioner(inner::AbstractPreconditioner, jvp!,
    b::AbstractVector{T}; kmax::Integer = 20, kharvest::Integer = 4,
    nritz::Integer = 0, kcandidate::Integer = 3*kmax,
    escalateafter::Integer = 1, ranktol::Real = eps(T)^(3//8),
    benefittol::Real = 1e-6,
    cycleharvest::Bool = true, state = nothing) where {T<:AbstractFloat}
    kmax >= 1 || throw(ArgumentError(lazy"`kmax` = $(kmax) must be at least 1."))
    kharvest >= 0 || throw(ArgumentError(
        lazy"`kharvest` = $(kharvest) must be nonnegative."))
    nritz >= 0 || throw(ArgumentError(
        lazy"`nritz` = $(nritz) must be nonnegative."))
    kharvest + nritz >= 1 || throw(ArgumentError(
        "`kharvest` and `nritz` cannot both be zero; the harvest would "*
        "produce no candidates."))
    kcandidate >= kmax || throw(ArgumentError(
        lazy"`kcandidate` = $(kcandidate) must be at least `kmax` = $(kmax)."))
    ranktol > 0 || throw(ArgumentError(
        lazy"`ranktol` = $(ranktol) must be positive."))
    benefittol >= 0 || throw(ArgumentError(
        lazy"`benefittol` = $(benefittol) must be nonnegative."))
    n = length(b)
    st = isnothing(state) ? FloquetState(b) : state
    size(st.X, 1) == n || throw(DimensionMismatch(
        lazy"the deflation state is for dimension $(size(st.X, 1)) but the system has $(n)."))
    length(st.source) == size(st.X, 2) || throw(ArgumentError(
        "the deflation state has one provenance tag per candidate column."))
    return FloquetPreconditioner(inner, jvp!, st,
        similar(b, n, 0), similar(b, n, 0), similar(b, n, 0), similar(b, 0),
        T[], Int(kmax), Int(kharvest), Int(nritz),
        Int(kcandidate), Int(escalateafter), T(ranktol), T(benefittol),
        cycleharvest, 0, 0, 0, size(st.X, 2) == 0)
end

function FloquetPreconditioner(inner::AbstractPreconditioner, jvp!,
    n::Integer; T::Type{<:AbstractFloat} = Float64, kwargs...)
    return FloquetPreconditioner(inner, jvp!, Vector{T}(undef, n); kwargs...)
end

deflationform(::FloquetPreconditioner) = :floquet
# the *active* rank, which is what enters the cost of an application; the
# candidate bank is reported separately
deflationsize(pc::FloquetPreconditioner) = size(pc.C, 2)
deflationrebuilds(pc::FloquetPreconditioner) = pc.rebuilds
candidatecount(pc::FloquetPreconditioner) = size(pc.state.X, 2)
deflationproducts(pc::FloquetPreconditioner) = pc.products

"""
    correctionstrengths(pc::FloquetPreconditioner)

The correction strength `eta = norm(x - inv(P)*J*x)/norm(x)` of each active
direction at the last rebuild. `eta` near zero means the base preconditioner
already handles that channel, `eta` of order one that a substantial part of
it is missing from the base. Diagnostic only.
"""
correctionstrengths(pc::FloquetPreconditioner) = pc.strength

pointmoved!(pc::FloquetPreconditioner) = (pc.fresh = false; pc)
stalled!(pc::FloquetPreconditioner) = (stalled!(pc.inner); pc)

# as `RecyclingPreconditioner`: withholding an escalation is only defensible
# while there is a subspace which might cover the deficiency instead
function escalatepreconditioner!(pc::FloquetPreconditioner)
    if size(pc.C, 2) == 0
        pc.escalationrequests = 0
        return escalatepreconditioner!(pc.inner)
    end
    pc.escalationrequests += 1
    if pc.escalationrequests < pc.escalateafter
        return false
    end
    pc.escalationrequests = 0
    return escalatepreconditioner!(pc.inner)
end

function updatepreconditioner!(pc::FloquetPreconditioner, x::AbstractVector)
    updatepreconditioner!(pc.inner, x)
    # Against an exact base the correction is identically zero: C = J*X, so
    # W = X - inv(J)*J*X = 0. Rebuilding it anyway would cost k products and
    # r solves against the full Jacobian, which after an escalation is most
    # of the Newton step. The candidates are kept in case the base is later
    # rebuilt restricted again.
    if isexactpreconditioner(pc.inner)
        _clearfloquet!(pc)
        return pc
    end
    _rebuildfloquet!(pc)
    return pc
end

function _clearfloquet!(pc::FloquetPreconditioner{TI,TJ,T}) where {TI,TJ,T}
    Xc = pc.state.X
    n = size(Xc, 1)
    pc.X = similar(Xc, n, 0)
    pc.C = similar(Xc, n, 0)
    pc.W = similar(Xc, n, 0)
    pc.coeff = similar(Xc, 0)
    pc.strength = T[]
    pc.fresh = true
    return pc
end

function seeddeflation!(pc::FloquetPreconditioner,
    Xnew::AbstractMatrix; source::Symbol = :external)
    before = candidatecount(pc)
    _bankcandidates!(pc, Xnew; source = source)
    candidatecount(pc) == before || (pc.fresh = false)
    return pc
end

# Append candidates to the bank without touching the active blocks. A
# candidate can arrive in the middle of a running GMRES, from a per-cycle
# harvest, and GMRES requires a fixed preconditioner: rebuilding on the
# next application would change the operator under a half-built Arnoldi
# factorization and invalidate it. Banked candidates therefore take effect
# at the next rebuild, which is triggered by `pointmoved!` at the next
# Newton step, by a base refresh, or by an external `seeddeflation!` --
# never from inside a solve.
function _bankcandidates!(pc::FloquetPreconditioner{TI,TJ,T},
    Xnew::AbstractMatrix; source::Symbol = :external) where {TI,TJ,T}
    size(Xnew, 2) == 0 && return pc
    size(Xnew, 1) == size(pc.state.X, 1) || throw(DimensionMismatch(
        lazy"candidates have length $(size(Xnew, 1)) but the state is of "*
        lazy"dimension $(size(pc.state.X, 1))."))
    # A candidate carries a direction, not a magnitude, and the residual
    # image factorization below weights a column by the size of its image.
    # Columns are therefore normalized on the way in, and one which is
    # numerically zero is dropped rather than turned into noise.
    B = copy(Xnew)
    kept = Int[]
    for j in axes(B, 2)
        nj = norm(view(B, :, j))
        if isfinite(nj) && nj > 0
            view(B, :, j) ./= nj
            push!(kept, j)
        end
    end
    isempty(kept) && return pc
    B = B[:, kept]
    st = pc.state
    st.X = size(st.X, 2) > 0 ? hcat(st.X, B) : B
    append!(st.source, fill(source, size(B, 2)))
    _trimcandidates!(pc)
    return pc
end

# Keep the bank bounded. The active directions of the last rebuild are the
# proven ones and sit at the front (`_rebuildfloquet!` writes them there),
# so the oldest *unproven* candidates are the ones dropped.
function _trimcandidates!(pc::FloquetPreconditioner)
    st = pc.state
    k = size(st.X, 2)
    k <= pc.kcandidate && return pc
    nactive = min(size(pc.X, 2), k)
    ndrop = k - pc.kcandidate
    # never drop an active direction; if the active set alone overruns the
    # bank the newest candidates go instead
    idx = if ndrop <= k - nactive
        vcat(1:nactive, (nactive+ndrop+1):k)
    else
        1:pc.kcandidate
    end
    st.X = st.X[:, idx]
    st.source = st.source[idx]
    return pc
end

"""
    _rebuildfloquet!(pc::FloquetPreconditioner)

Rebuild the active basis from the candidate bank at the current point:
the exact images `Y0 = J*X0`, the rank-revealing factorization which makes
`J*X = C` with `C'C = I`, the base solves `inv(P)*C`, the correction
strength filter, and the compression to `kmax`.

Costs `k` Jacobian products and `r` base solves for `k` candidates and
active rank `r`.
"""
function _rebuildfloquet!(pc::FloquetPreconditioner{TI,TJ,T}) where {TI,TJ,T}
    X0 = pc.state.X
    n, k = size(X0)
    k == 0 && return _clearfloquet!(pc)

    # Exact residual images. The product is handed contiguous vectors
    # rather than column views, which a device direct solver cannot bind a
    # descriptor to.
    cin = similar(X0, n)
    cout = similar(X0, n)
    Y0 = similar(X0)
    for j in 1:k
        copyto!(cin, view(X0, :, j))
        pc.jvp!(cout, cin)
        copyto!(view(Y0, :, j), cout)
    end
    pc.products += k

    # Equalize the candidates by the size of their images before the
    # factorization. Without this the rank test measures how large a
    # candidate happens to be rather than how independent its image is, and
    # a carried-forward direction (whose image is already of norm one)
    # would be weighted against a freshly harvested one on a scale set by
    # the base preconditioner. Scaling a column of `X0` and of `Y0` by the
    # same factor preserves `J*X0 = Y0`, so the identity the whole
    # construction rests on is untouched. A candidate whose image vanishes
    # is in the null space of `J` and cannot be given a residual image at
    # all, so it is dropped here.
    X0s = copy(X0)
    ynorm = [norm(view(Y0, :, j)) for j in 1:k]
    ymax = maximum(ynorm; init = zero(T))
    live = [j for j in 1:k if isfinite(ynorm[j]) && ynorm[j] > eps(T)*ymax]
    isempty(live) && return _clearfloquet!(pc)
    if length(live) < k
        X0s = X0s[:, live]; Y0 = Y0[:, live]; ynorm = ynorm[live]
    end
    for (j, s) in enumerate(ynorm)
        view(X0s, :, j) ./= s
        view(Y0, :, j) ./= s
    end

    # The rank-revealing step, through the `k` by `k` Gram matrix of the
    # images: `Y0'*Y0 = V*S^2*V'` is formed where `Y0` lives and only the
    # small matrix crosses to the host, where it is factorized; the `k` by
    # `r` transform goes back and the two large products stay on the
    # device. The Gram form squares the condition number, so its
    # resolution is about `sqrt(k*eps)` on the singular values, below the
    # default `ranktol`; a cleanup pass below restores `C'C = I` to working
    # precision.
    E = eigen(Symmetric(Array(Y0'*Y0)))
    sv = sqrt.(max.(E.values, zero(T)))
    smax = maximum(sv; init = zero(T))
    smax > 0 || return _clearfloquet!(pc)
    keep = [i for i in eachindex(sv) if sv[i] > pc.ranktol*smax]
    r = length(keep)
    r >= 1 || return _clearfloquet!(pc)
    Tsmall = E.vectors[:, keep]*Diagonal(inv.(sv[keep]))
    Tdev = _hostbuilt(X0, Tsmall)
    X = X0s*Tdev
    C = Y0*Tdev
    Rc = cholesky(Symmetric(Array(C'*C)), check = false)
    if issuccess(Rc)
        Rinv = _hostbuilt(X0, Matrix(inv(Rc.U)))
        C = C*Rinv
        X = X*Rinv
    end

    # `W = X - inv(P)*C` folds the base solve of `C` into the update, so an
    # application needs no solve beyond the base solve of the residual
    # itself.
    BC = similar(C)
    for j in 1:r
        copyto!(cin, view(C, :, j))
        applypreconditioner!(cout, pc.inner, cin)
        copyto!(view(BC, :, j), cout)
    end
    W = X .- BC

    # Rotate into the principal directions of `W'W` *before* measuring
    # anything. The rank-revealing factorization above fixes `range(C)` but
    # not the basis of it: equalizing the candidates by their image norms
    # tends to make the columns of `Y0` nearly orthonormal, and the singular
    # values of a nearly orthonormal block are all nearly equal, which
    # leaves the right singular vectors — and so the columns of `X`, `C` and
    # `W` — an essentially arbitrary rotation of the subspace. Correction
    # strength is not invariant under that rotation. Measured in the
    # arbitrary basis, a subspace whose correction genuinely has rank two
    # can show six columns of comparable `eta`, none of them droppable,
    # because each holds a mixture of the directions the base misses and the
    # directions it handles. The eigenvectors of `W'W` are the basis in
    # which the correction is diagonal, so they are where the split between
    # the two is visible at all. The transform is orthonormal, so `J*X = C`
    # and `C'C = I` both survive it.
    G = Hermitian(Array(W'*W))
    E = eigen(G)
    ord = sortperm(E.values; rev = true)
    R = _hostbuilt(X0, E.vectors[:, ord])
    X = X*R; C = C*R; W = W*R

    # The correction strength test, now in the basis where it means
    # something: `eta` is the fraction of the direction the base
    # preconditioner fails to reconstruct from its own image. A direction
    # the base already inverts accurately contributes nothing but cost,
    # whatever its physical interest. The columns are ordered by correction
    # magnitude, so the survivors are a prefix and the cap at `kmax` keeps
    # the directions carrying the most correction.
    strength = [norm(view(W, :, j))/
        max(norm(view(X, :, j)), floatmin(T)) for j in 1:r]
    keep = [j for j in 1:r if strength[j] > pc.benefittol]
    if isempty(keep)
        _clearfloquet!(pc)
        # the candidates were all handled by the base; keep them, since a
        # base rebuilt restricted differently may not handle them again
        pc.rebuilds += 1
        return pc
    end
    length(keep) > pc.kmax && (keep = keep[1:pc.kmax])
    if length(keep) < r
        X = X[:, keep]; C = C[:, keep]; W = W[:, keep]
        strength = strength[keep]
        r = length(keep)
    end

    pc.X = X; pc.C = C; pc.W = W
    pc.coeff = similar(X0, r)
    pc.strength = strength
    # Carry the active corrections forward as the seed of the next point's
    # bank. They are the directions which survived both filters, and at a
    # neighboring operating point they are the best available guess at the
    # same physical channels.
    pc.state.X = X
    pc.state.source = fill(:active, r)
    pc.rebuilds += 1
    pc.fresh = true
    return pc
end

# Bring the active blocks up to date when the point has moved or a
# candidate has been seeded, but the base has not been rebuilt. There is
# only one build path here: unlike the Galerkin forms, nothing is carried
# between rebuilds, because `inv(P)*U` is never formed in the first place.
function _refreshfloquet!(pc::FloquetPreconditioner)
    pc.fresh && return pc
    if isexactpreconditioner(pc.inner)
        pc.fresh = true
        return pc
    end
    return _rebuildfloquet!(pc)
end

function applypreconditioner!(z::AbstractVector, pc::FloquetPreconditioner,
    r::AbstractVector)
    _refreshfloquet!(pc)
    applypreconditioner!(z, pc.inner, r)
    size(pc.C, 2) > 0 || return z
    mul!(pc.coeff, transpose(pc.C), r)
    mul!(z, pc.W, pc.coeff, true, true)
    return z
end

"""
    harmonicritznearzero(Hbar::AbstractMatrix, nkeep::Integer)

The harmonic Ritz vectors of an Arnoldi factorization whose Ritz values lie
nearest the origin, returned as the columns of a *real* `m` by `p` matrix of
coefficients in the Arnoldi basis, `p <= 2*nkeep`.

For `A*V = V_{m+1}*Hbar` a harmonic Ritz pair `(theta, y)` makes the
eigenpair residual `A*V*y - theta*V*y` orthogonal to `A*K_m`, which is the
generalized eigenproblem

    Hbar'*Hbar*y = theta*H'*y

with `H` the square leading block. That is the pencil solved here, on the
host, rather than the equivalent `H + abs2(h)*inv(H')*e*e'`, which applies
the inverse of `H` exactly where it is worst conditioned, near the
directions of interest; a singular `H` gives infinite Ritz values, which
are simply not selected. Harmonic Ritz values approximate the eigenvalues
of `A` nearest zero far better than the ordinary Ritz values do, whose
restarted GMRES cycle approximates the *outer* spectrum, the part deflation
has no use for.

`Hbar` must be the *Arnoldi* Hessenberg (`GMRESWorkspace.Harnoldi`), not
the least squares matrix the Givens rotations leave in `GMRESWorkspace.H`:
the singular values of the two agree, the pencil does not. A real matrix
has complex eigenpairs in conjugate pairs, and both the real and the
imaginary part of such an eigenvector lie in the real invariant subspace
the pair spans, so each is returned as its own real column; the
duplication is harmless because the residual-image factorization removes
whatever is redundant.
"""
function harmonicritznearzero(Hbar::AbstractMatrix{T},
    nkeep::Integer) where {T<:AbstractFloat}
    m = size(Hbar, 2)
    (nkeep >= 1 && m >= 1) || return zeros(T, m, 0)
    Hb = Matrix(Hbar)
    H = Hb[1:m, 1:m]
    E = try
        eigen(Hb'*Hb, Matrix(transpose(H)))
    catch
        return zeros(T, m, 0)
    end
    vals, vecs = E.values, E.vectors
    finite = [i for i in eachindex(vals) if isfinite(vals[i])]
    isempty(finite) && return zeros(T, m, 0)
    order = sort(finite; by = i -> abs(vals[i]))
    Y = Matrix{T}(undef, m, 0)
    for i in order
        size(Y, 2) >= 2*nkeep && break
        v = view(vecs, :, i)
        all(isfinite, v) || continue
        vr = real.(v)
        norm(vr) > 0 && (Y = hcat(Y, vr))
        vi = imag.(v)
        # the imaginary part is the second real direction of a conjugate
        # pair, and is exactly zero for a real eigenvalue
        norm(vi) > sqrt(eps(T))*norm(vr) && (Y = hcat(Y, vi))
        count(j -> abs(vals[j]) <= abs(vals[i]), finite) >= nkeep &&
            size(Y, 2) >= nkeep && break
    end
    return Y
end

"""
    _harvestfloquet!(pc::FloquetPreconditioner, Vj, Hj)

Harvest candidate correction vectors from an Arnoldi factorization with
basis `Vj` (`n` by `j`) and rectangular Hessenberg `Hj` (`j + 1` by `j`,
host resident) and return `pc`.

Two families are taken, following the principle that on a strongly
nonnormal operator neither eigenvalue nor singular value information is
reliably the better target:

- the `kharvest` right singular vectors of `Hj` with the smallest singular
  values, the directions this Krylov space found the preconditioned
  operator shrinks most;
- the `nritz` harmonic Ritz directions nearest zero
  ([`harmonicritznearzero`](@ref)), which approximate its near-singular
  eigendirections.

Both are Krylov-space directions. Each is lifted through the Arnoldi basis
and then mapped through the *current complete* preconditioner, `x = B*u`,
into the physical state space before being stored, because that is the
space in which a direction still means something once the Jacobian has
moved: at the next point the exact image `J*x` is recomputed and the coarse
space rebuilt around it, whereas a stored residual direction would have to
be reinterpreted against an operator it was never measured on.

Nothing is rebuilt here. The preconditioner must not change under a running
GMRES, so a harvest only appends to the candidate bank and marks the active
blocks stale; the rebuild happens at the next Newton step.
"""
function _harvestfloquet!(pc::FloquetPreconditioner{TI,TJ,T}, Vj::AbstractMatrix,
    Hj::AbstractMatrix) where {TI,TJ,T}
    # nothing to deflate against an exact base; see `updatepreconditioner!`
    isexactpreconditioner(pc.inner) && return pc
    j = size(Vj, 2)
    j >= 1 || return pc
    Hh = Matrix(Hj)
    Y = Matrix{T}(undef, j, 0)

    if pc.kharvest > 0
        F = svd(Hh)
        nv = size(F.V, 2)
        p = min(pc.kharvest, nv)
        p > 0 && (Y = hcat(Y, F.V[:, nv-p+1:nv]))
    end
    if pc.nritz > 0
        Y = hcat(Y, harmonicritznearzero(Hh, pc.nritz))
    end
    size(Y, 2) == 0 && return pc

    # lift into the state space, then into correction coordinates
    U = Vj*_hostbuilt(Vj, Y)
    n = size(U, 1)
    X = similar(U)
    cin = similar(U, n)
    cout = similar(U, n)
    for c in axes(U, 2)
        copyto!(cin, view(U, :, c))
        applypreconditioner!(cout, pc, cin)
        copyto!(view(X, :, c), cout)
    end
    return _bankcandidates!(pc, X; source = :gmres_floquet)
end

function harvest!(pc::FloquetPreconditioner, ws::GMRESWorkspace,
    out::NamedTuple)
    j = harvestdimension(ws, out)
    j >= 1 || return pc
    return _harvestfloquet!(pc, view(ws.V, :, 1:j),
        Matrix(view(ws.Harnoldi, 1:j+1, 1:j)))
end

usescycleharvest(pc::FloquetPreconditioner) = pc.cycleharvest

function harvestcycle!(pc::FloquetPreconditioner, ws::GMRESWorkspace,
    j::Integer)
    j >= 1 || return pc
    # the *Arnoldi* Hessenberg: `ws.H` has been triangularized by the Givens
    # rotations by now, which the singular directions survive and the
    # harmonic Ritz pencil does not
    return _harvestfloquet!(pc, view(ws.V, :, 1:j),
        Matrix(view(ws.Harnoldi, 1:j+1, 1:j)))
end
