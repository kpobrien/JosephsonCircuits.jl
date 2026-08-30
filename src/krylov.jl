
"""
    FunctionOperator(f!, n)

Wraps an in-place product `f!(y, v)` as a `mul!`-able operator of dimension
`n`, so the Krylov machinery can be written against `mul!` alone while
still accepting a bare closure.
"""
struct FunctionOperator{F}
    f!::F
    n::Int
end
Base.size(A::FunctionOperator) = (A.n, A.n)
Base.size(A::FunctionOperator, i::Integer) = A.n
Base.eltype(::FunctionOperator) = Float64
LinearAlgebra.mul!(y::AbstractVector, A::FunctionOperator, v::AbstractVector) =
    (A.f!(y, v); y)

"""
    asoperator(A, n)

`A` itself when it is already `mul!`-able, and a [`FunctionOperator`](@ref)
when it is a bare in-place product.

The Krylov solvers of this package apply the Jacobian through `mul!`, which
is the interface every external linear algebra package in the ecosystem
consumes. Normalizing here means a `JacobianOperator` can be handed
straight through with no wrapper, while the older closure form keeps
working.
"""
asoperator(f::Function, n::Integer) = FunctionOperator(f, Int(n))
asoperator(A, ::Integer) = A



"""
    AbstractPreconditioner

Supertype of the preconditioners of [`nlsolvekrylov!`](@ref). A preconditioner
`pc` approximates the Jacobian of the nonlinear system by something which can be
inverted cheaply, and must implement two methods:

- `updatepreconditioner!(pc, x)`: rebuild `pc` at the point `x`,
- `applypreconditioner!(z, pc, r)`: overwrite `z` with `inv(P)*r`.

The Jacobian itself is never required, only its approximation, so the nonlinear
solve stays matrix-free. See [`ModeCouplingPreconditioner`](@ref) for the harmonic balance Jacobian with
its mode coupling restricted, and [`RecyclingPreconditioner`](@ref) for that
wrapped in a recycled deflation subspace.
"""
abstract type AbstractPreconditioner end

"""
    updatepreconditioner!(pc::AbstractPreconditioner, x::AbstractVector)

Rebuild the preconditioner `pc` at the point `x` and return `pc`.
"""
function updatepreconditioner! end

"""
    applypreconditioner!(z::AbstractVector, pc::AbstractPreconditioner,
        r::AbstractVector)

Overwrite `z` with the action of the inverse of the preconditioner `pc` on `r`
and return `z`.
"""
function applypreconditioner! end

"""
    KrylovSolveInfo

Diagnostics for a single linear solve of [`nlsolvekrylov!`](@ref).

There is one record per GMRES call, not per Newton step. A Newton step can
contain up to three: the first solve, a refresh-and-retry when the
preconditioner has drifted, and a rescue solve when the resulting direction is
not a descent direction. Recording per solve rather than per step is what makes
the retry structure visible; a one-entry-per-step vector cannot express it, and
the extra solves are exactly where the expensive failures hide.

# Fields
- `iteration`: the outer Newton step this solve belongs to. Repeated values
    mark a retry or a rescue.
- `role`: `:step` for the solve which produces the step, `:retry` for the
    refresh-and-retry after a solve that did not reach tolerance, `:rescue`
    for the solve after a non-descent direction.
- `normF`: the nonlinear residual norm at the point of the solve, the right
    hand side norm.
- `forcing`: the requested relative tolerance, the Eisenstat-Walker forcing
    term for this step.
- `residualratio`: the *achieved* explicit linear residual ratio. The pair
    (`forcing`, `residualratio`) is the diagnosis of an oversolving or
    undersolving forcing sequence: a loose `forcing` satisfied in one
    iteration at a `residualratio` near one is a solve that did almost
    nothing while reporting success.
- `iterations`, `cycles`, `reason`: as returned by [`gmres!`](@ref).
- `refreshed`: whether the preconditioner was rebuilt immediately before this
    solve.
- `escalated`: whether [`escalatepreconditioner!`](@ref) grew the
    preconditioner after this solve.
- `stagnated`: whether the step was discarded and replaced by the
    preconditioner solve.
- `slope`: `dot(F, J*deltax)/normF^2`, a scale free measure of direction
    quality from an exact matrix-free product. `-1` is the Newton direction;
    a small magnitude is a weak descent direction, which is what a loose
    solve against a stale preconditioner produces. `NaN` when not computed.
- `alpha`, `backtracks`, `armijo`: the linesearch outcome for the step this
    solve produced. `armijo` is false when the step was merely the best
    decreasing trial rather than an Armijo accepted one.
- `time`: seconds since the start of the nonlinear solve, so that the
    residual history can be plotted against wall time with the refreshes and
    escalations marked.
"""
struct KrylovSolveInfo
    iteration::Int
    role::Symbol
    normF::Float64
    forcing::Float64
    residualratio::Float64
    iterations::Int
    cycles::Int
    reason::Symbol
    refreshed::Bool
    escalated::Bool
    stagnated::Bool
    slope::Float64
    alpha::Float64
    backtracks::Int
    armijo::Bool
    time::Float64
    # The preconditioner at the time of this solve. `escalated` alone cannot
    # distinguish an escalation which was never requested from one which was
    # requested and refused, and the difference between those two is the
    # whole of the behaviour below.
    escalationrequested::Bool
    # the width of the deflation subspace and how many times it has been
    # rebuilt: the only visible sign of what recycling costs
    deflationsize::Int
    deflationrebuilds::Int
    # wall time applying the preconditioner. A configuration which takes more
    # iterations in less time is running against a weaker and cheaper
    # preconditioner, which an iteration count alone cannot distinguish from
    # one which is simply worse.
    precondtime::Float64
end

function Base.show(io::IO, ::MIME"text/plain", k::KrylovSolveInfo)
    print(io, "KrylovSolveInfo(", k.iteration, " ", k.role,
        ": |F|=", round(k.normF, sigdigits = 3),
        " eta=", round(k.forcing, sigdigits = 2),
        " achieved=", round(k.residualratio, sigdigits = 2),
        " its=", k.iterations, "/", k.cycles, " ", k.reason,
        k.refreshed ? " refreshed" : "",
        k.escalated ? " escalated" :
            (k.escalationrequested ? " escalation-refused" : ""),
        k.deflationsize > 0 ? " deflation=$(k.deflationsize)" : "",
        k.stagnated ? " stagnated" : "",
        " slope=", round(k.slope, sigdigits = 2),
        " alpha=", round(k.alpha, sigdigits = 2),
        k.armijo ? "" : " (non-Armijo)", ")")
end

"""
    krylovtotaliterations(records)

Total Arnoldi steps across every linear solve in a vector of
[`KrylovSolveInfo`](@ref), including the retries and rescues.
"""
krylovtotaliterations(records) = sum(k -> k.iterations, records; init = 0)


# the step outcome is only known after the linesearch, so the record of the
# solve which produced it is completed in place rather than deferred
function _withstep(k::KrylovSolveInfo, slope, alpha, backtracks, armijo)
    return KrylovSolveInfo(k.iteration, k.role, k.normF, k.forcing,
        k.residualratio, k.iterations, k.cycles, k.reason, k.refreshed,
        k.escalated, k.stagnated, slope, alpha, backtracks, armijo, k.time,
        k.escalationrequested, k.deflationsize, k.deflationrebuilds,
        k.precondtime)
end

function _withescalated(k::KrylovSolveInfo)
    return KrylovSolveInfo(k.iteration, k.role, k.normF, k.forcing,
        k.residualratio, k.iterations, k.cycles, k.reason, k.refreshed,
        true, k.stagnated, k.slope, k.alpha, k.backtracks, k.armijo, k.time,
        k.escalationrequested, k.deflationsize, k.deflationrebuilds,
        k.precondtime)
end

# an escalation which the preconditioner refused: the request happened, the
# strengthening did not
function _withescalationrefused(k::KrylovSolveInfo)
    return KrylovSolveInfo(k.iteration, k.role, k.normF, k.forcing,
        k.residualratio, k.iterations, k.cycles, k.reason, k.refreshed,
        k.escalated, k.stagnated, k.slope, k.alpha, k.backtracks, k.armijo,
        k.time, true, k.deflationsize, k.deflationrebuilds, k.precondtime)
end


"""
    LinearAlgebra.ldiv!(z, pc::AbstractPreconditioner, r)
    LinearAlgebra.ldiv!(pc::AbstractPreconditioner, r)
    pc \\ r

Apply the inverse of the preconditioner, forwarding to
[`applypreconditioner!`](@ref).

`ldiv!` is the de facto interface every external Krylov package consumes:
`Pl`/`Pr` in LinearSolve.jl and IterativeSolvers.jl, `M`/`N` in Krylov.jl,
`Pl` in BifurcationKit's linear solvers. The domain knowledge of this
solver lives in the preconditioner, so defining `ldiv!` is what makes
[`ModeCouplingPreconditioner`](@ref) and [`RecyclingPreconditioner`](@ref)
reusable outside the package without an adapter.
"""
LinearAlgebra.ldiv!(z::AbstractVector, pc::AbstractPreconditioner,
    r::AbstractVector) = applypreconditioner!(z, pc, r)

function LinearAlgebra.ldiv!(pc::AbstractPreconditioner, r::AbstractVector)
    z = applypreconditioner!(similar(r), pc, r)
    copyto!(r, z)
    return r
end

Base.:\(pc::AbstractPreconditioner, r::AbstractVector) =
    applypreconditioner!(similar(r), pc, r)

"""
    LinearAlgebra.mul!(z, pc::AbstractPreconditioner, r)

Apply the inverse of the preconditioner, spelled as a multiplication.

The two conventions in the ecosystem disagree. LinearSolve.jl and
IterativeSolvers.jl take a preconditioner and *divide* by it, so they call
`ldiv!`. Krylov.jl takes `M` and `N` to be operators which already
represent the inverse and *multiplies*, so it calls `mul!`. Supporting both
is what lets the same preconditioner object be handed to either without an
adapter, which is the whole point of shipping it as an object.
"""
LinearAlgebra.mul!(z::AbstractVector, pc::AbstractPreconditioner,
    r::AbstractVector) = applypreconditioner!(z, pc, r)


"""
    escalatepreconditioner!(pc::AbstractPreconditioner)

Make the preconditioner `pc` a better approximation of the Jacobian, at greater
cost, and return `true`; return `false` when it cannot be improved further.
Called by [`nlsolvekrylov!`](@ref) when GMRES stagnates repeatedly, which is the
symptom of a preconditioner too crude for the problem. The default method
returns `false`, which is correct for any preconditioner that is already exact
or has no cheaper/costlier settings.
"""
escalatepreconditioner!(::AbstractPreconditioner) = false

"""
    deflationsize(pc)

The number of vectors in the deflation subspace, zero for a preconditioner
which does not recycle.
"""
deflationsize(::AbstractPreconditioner) = 0

"""
    deflationrebuilds(pc)

How many times the deflation subspace has been rebuilt.
"""
deflationrebuilds(::AbstractPreconditioner) = 0

"""
    escalationrefusals(pc)

How many escalation requests were absorbed without being passed on.
"""
escalationrefusals(::AbstractPreconditioner) = 0

"""
    preconditionerlevel(pc)

How many modes the preconditioner couples, as a measure of its strength, or
`-1` when that is not meaningful. Reported so an escalation can be checked
to have taken effect rather than merely to have been requested.
"""
preconditionerlevel(::AbstractPreconditioner) = -1

"""
    GMRESWorkspace{T<:AbstractFloat}

Preallocated storage for [`gmres!`](@ref) with a restart length of `m` on a
system of dimension `n`. Holds the `n x (m+1)` Arnoldi basis `V`, the
`(m+1) x m` Hessenberg matrix `H`, the Givens rotations `cs` and `sn` which
reduce it, the least squares right hand side `s`, its solution `y`, and three
length `n` work vectors.

The dominant cost is `V`, which is `n*(m+1)` numbers, so `m` trades memory and
orthogonalization work against restart frequency.
"""
struct GMRESWorkspace{T<:AbstractFloat,TV<:AbstractVector{T},TM<:AbstractMatrix{T}}
    # system sized, and therefore allocated like the right hand side, so on a
    # device backend they live on the device
    V::TM
    w::TV
    z::TV
    u::TV
    # the projected quantities, which are `m` x `m` at most. These stay on the
    # host: the Givens rotations and the back substitution index them entry by
    # entry, which is exactly what a device array must not be asked to do, and
    # they are small enough that keeping them host resident costs one small
    # transfer per Arnoldi step rather than a kernel launch per entry.
    H::Matrix{T}
    cs::Vector{T}
    sn::Vector{T}
    s::Vector{T}
    y::Vector{T}
    # device resident staging buffers for the two block Gram-Schmidt
    # projections, so the coefficient vectors are formed where `V` lives and
    # only the finished column of `H` crosses to the host
    hd::TV
    cd::TV
end

function GMRESWorkspace(n::Integer, m::Integer, ::Type{T} = Float64) where {T<:AbstractFloat}
    return GMRESWorkspace(Vector{T}(undef, n), m)
end

"""
    GMRESWorkspace(b::AbstractVector{T}, m::Integer)

Build a workspace for a restart length of `m` on a system whose right hand
side is `b`. The system sized arrays are allocated with `similar(b)`, so they
live wherever `b` does and the iteration runs on that device; the projected
`m` x `m` quantities are host arrays regardless, because the Givens rotations
and the back substitution index them entry by entry.
"""
function GMRESWorkspace(b::AbstractVector{T}, m::Integer) where {T<:AbstractFloat}
    n = length(b)
    n >= 0 || throw(ArgumentError(lazy"`n` = $(n) must be nonnegative."))
    m >= 1 || throw(ArgumentError(lazy"the restart length `m` = $(m) must be at least 1."))
    return GMRESWorkspace{T,typeof(similar(b)),typeof(similar(b, n, m+1))}(
        similar(b, n, m + 1), similar(b), similar(b), similar(b),
        zeros(T, m + 1, m),
        Vector{T}(undef, m), Vector{T}(undef, m),
        Vector{T}(undef, m + 1), Vector{T}(undef, m),
        similar(b, m), similar(b, m))
end


"""
    harvest!(pc::AbstractPreconditioner, ws::GMRESWorkspace, out::NamedTuple)

Give the preconditioner `pc` the Arnoldi factorization a solve just built, so
it can extract information for the *next* solve. `out` is the named tuple
returned by [`gmres!`](@ref). Called by [`nlsolvekrylov!`](@ref) after every
GMRES call. The default does nothing, which is correct for any preconditioner
that does not recycle.

Only the *last* restart cycle is still present in the workspace, so
implementations must derive the usable Arnoldi dimension from
`out.iterations` and `out.cycles` rather than from `out.iterations` alone.
"""
harvest!(pc::AbstractPreconditioner, ::GMRESWorkspace, ::NamedTuple) = pc

"""
    isexactpreconditioner(pc::AbstractPreconditioner)

Whether `pc` currently applies the exact Jacobian, so that a deflation or
composition layered on top of it can contribute nothing. `false` for any
preconditioner that does not say otherwise.
"""
isexactpreconditioner(::AbstractPreconditioner) = false

"""
    RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!, n::Integer;
        kmax = 40, kharvest = 8)

Augments any base preconditioner with a recycled deflation subspace, carried
across Newton steps.

The base preconditioner `inner` leaves a handful of directions on which the
preconditioned operator `J*inv(P)` is nearly singular. Those directions are
what stalls GMRES, because its residual polynomial is pinned at `p(0) = 1` and
so cannot be small near the origin. Rather than enlarging `inner` until it
happens to cover them, this wrapper *measures* them: after each solve
[`harvest!`](@ref) extracts the directions the Arnoldi factorization found the
operator shrinks most, and the next solve deflates them explicitly.

Given a subspace `U` and `Z = inv(P)*U`, `Y = J*Z`, `G = Z'Y`, the deflated
preconditioner is

    inv(Pdef)*v = inv(P)*v + W*(inv(G)*(Z'*v)),   W = Z - inv(P)*Y

which makes `range(Y)` an exact eigenspace of `J*inv(Pdef)` at eigenvalue one.
Folding `inv(P)*Y` into `W` at build time keeps the application to one base
solve plus two dense mat-vecs with an `n` by `k` matrix, rather than three.

The subspace is *accumulated* rather than replaced: each harvest is appended
and re-orthogonalized, and the result is trimmed back to `kmax` by keeping the
eigenvectors of `Y'Y` with the smallest eigenvalues, that is, the subspace the
operator shrinks most. Accumulating measurably beats rediscovering the
subspace at every step, because the deficient directions move slowly while the
Newton iterate does not.

The intended base is the mode block diagonal
([`ModeCouplingPreconditioner`](@ref) with `couplingmodes = :none`). That
combination is the scalable one: the base is a batch of small independent
factorizations, one per mode, and the deflation is dense level 2 BLAS, so
neither part needs a large sparse factorization. It is also, for the same
reason, the part of this solver that ports to a GPU.

`escalatepreconditioner!` is delegated to `inner`, so the usual safety net
still applies when recycling alone is not enough.

!!! note
    The pairing this is built for is `couplingmodes = :none` plus recycling.
    Escalation of the base is deliberately reluctant underneath this wrapper;
    see `escalateafter`.
"""
mutable struct RecyclingPreconditioner{TI,TJ,T<:AbstractFloat,TM<:AbstractMatrix{T},TV<:AbstractVector{T}} <: AbstractPreconditioner
    const inner::TI
    const jvp!::TJ
    # system sized, and therefore allocated like the vectors of the system, so
    # on a device backend the deflation is a pair of device gemvs
    U::TM
    Z::TM
    W::TM
    # the k x k projected inverse and its two work vectors, which are small.
    # `Ginv` lives with the basis rather than on the host: the apply
    # multiplies it against `a` and `b`, which are device resident, so keeping
    # it on the host would mix the two in a `mul!`. Only the dense
    # factorizations that *build* it run on the host; see
    # `_rebuilddeflation!`.
    Ginv::TM
    a::TV
    b::TV
    const kmax::Int
    const kharvest::Int
    const escalateafter::Int
    escalationrequests::Int
    rebuilds::Int
    # escalation requests absorbed without being passed to `inner`
    refusals::Int
end

function RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!,
    n::Integer; kmax::Integer = 20, kharvest::Integer = 8,
    escalateafter::Integer = 1, T::Type{<:AbstractFloat} = Float64)
    return RecyclingPreconditioner(inner, jvp!, Vector{T}(undef, n);
        kmax = kmax, kharvest = kharvest, escalateafter = escalateafter)
end

"""
    RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!,
        b::AbstractVector; kmax = 20, kharvest = 8, escalateafter = 1)

Build the deflation wrapper on a system whose vectors are like `b`. The
deflation subspace and the two blocks derived from it are allocated with
`similar(b)`, so they live wherever `b` does and the application is a base
solve plus two device gemvs. The `k` x `k` projected inverse stays on the
host, where it is formed by a dense factorization of a matrix no larger than
`kmax`.
"""
function RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!,
    b::AbstractVector{T}; kmax::Integer = 20, kharvest::Integer = 8,
    escalateafter::Integer = 1) where {T<:AbstractFloat}
    kmax >= 1 || throw(ArgumentError(lazy"`kmax` = $(kmax) must be at least 1."))
    1 <= kharvest <= kmax || throw(ArgumentError(
        lazy"`kharvest` = $(kharvest) must satisfy 1 <= kharvest <= kmax = $(kmax)."))
    n = length(b)
    return RecyclingPreconditioner(inner, jvp!, similar(b, n, 0),
        similar(b, n, 0), similar(b, n, 0), similar(b, 0, 0),
        similar(b, 0), similar(b, 0),
        Int(kmax), Int(kharvest), Int(escalateafter), 0, 0, 0)
end

# Escalating the base and recycling are two answers to the same problem, and
# `escalateafter` lets recycling defer the base's escalation to absorb the
# deficiency itself. The default is 1 -- no deferral. It used to be 3, from a
# measurement in the restart-30 era where a failed linear solve cost at most
# ~120 iterations and a premature escalation cost more than the rest of the
# solve (3.25 s against 1.36 s at 128 cells). At the current restart of 400
# with 4 cycles, one deferred escalation costs up to 1600 Arnoldi iterations,
# and the deferral only pays if the deflation improves between failures.
# Measured on a strongly driven 64 junction RPM line it does not: three
# consecutive deferred solves achieved residual ratios of 0.38, 0.84 and 0.99
# -- monotonically worse -- and the deferral tripled the total linear
# iteration count (2306 to 6384) before escalating anyway. When the deflation
# genuinely absorbs the deficiency (128 junctions, kmax = 40), every solve
# converges and no escalation request ever arrives, so no deferral is needed
# for recycling to win.
function escalatepreconditioner!(pc::RecyclingPreconditioner)
    # Withholding an escalation is only defensible when there is a deflation
    # subspace which might cover the deficiency instead. While `U` is empty
    # the recycling has nothing to absorb it with, so the throttle delays a
    # real remedy for a mechanism which is not running.
    #
    # The delay is expensive. On a two tone travelling wave amplifier the
    # first escalation takes the mode coupling preconditioner from block
    # diagonal to fully coupled, after which every Newton step converges in
    # one Krylov iteration. Absorbing the first two requests left five Newton
    # steps solving against the unescalated preconditioner, reaching relative
    # residuals of 0.46, 0.48, 0.79, 0.89 and 0.93 and spending 324
    # iterations to do it.
    if size(pc.U, 2) == 0
        pc.escalationrequests = 0
        return escalatepreconditioner!(pc.inner)
    end
    pc.escalationrequests += 1
    if pc.escalationrequests < pc.escalateafter
        pc.refusals += 1
        return false
    end
    pc.escalationrequests = 0
    return escalatepreconditioner!(pc.inner)
end

deflationsize(pc::RecyclingPreconditioner) = size(pc.U, 2)
deflationrebuilds(pc::RecyclingPreconditioner) = pc.rebuilds
escalationrefusals(pc::RecyclingPreconditioner) = pc.refusals
preconditionerlevel(pc::RecyclingPreconditioner) = preconditionerlevel(pc.inner)

function updatepreconditioner!(pc::RecyclingPreconditioner, x::AbstractVector)
    updatepreconditioner!(pc.inner, x)
    # Against an exact inner the correction is identically zero: Y = J*P^-1*U
    # = U, so W = Z - P^-1*Y = 0. Rebuilding it anyway costs 2k inner solves
    # and k Jacobian products per Newton step for nothing, which after an
    # escalation to the full Jacobian -- where every linear solve takes one
    # iteration -- is the dominant cost of the whole step. The subspace `U` is
    # kept, in case the base is ever rebuilt restricted again.
    if isexactpreconditioner(pc.inner)
        n = size(pc.U, 1)
        pc.Z = similar(pc.U, n, 0); pc.W = similar(pc.U, n, 0)
        pc.Ginv = similar(pc.U, 0, 0)
        pc.a = similar(pc.U, 0); pc.b = similar(pc.U, 0)
        return pc
    end
    _rebuilddeflation!(pc)
    return pc
end

# A small dense matrix built on the host, moved to wherever `proto` lives.
# `eigen`, `svd` and `cholesky` are scalar indexed dense kernels, so the k x k
# projected quantities are brought across, factorized on the host and the
# result sent back. This is the same split `GMRESWorkspace` makes for the
# Hessenberg and the Givens rotations, and for the same reason: what crosses
# is k x k once per rebuild, not anything of the system dimension.
_hostbuilt(proto::AbstractMatrix, A::AbstractMatrix) =
    copyto!(similar(proto, size(A)...), A)

# rebuild Z, W and inv(G) at the current point. Y = J*Z is recomputed rather
# than carried over from the previous step's Arnoldi: the free but stale
# version was measured to roughly double the Krylov iteration count, which
# costs more than the mat-vecs it saves.
function _rebuilddeflation!(pc::RecyclingPreconditioner{TI,TJ,T}) where {TI,TJ,T}
    k = size(pc.U, 2)
    n = size(pc.U, 1)
    if k == 0
        pc.Z = similar(pc.U, n, 0); pc.W = similar(pc.U, n, 0)
        pc.Ginv = similar(pc.U, 0, 0)
        pc.a = similar(pc.U, 0); pc.b = similar(pc.U, 0)
        return pc
    end
    Z = similar(pc.U)
    Y = similar(pc.U)
    # the base solve and the product are handed contiguous vectors rather than
    # column views: a device direct solver wants a vector it can bind a
    # descriptor to, and a strided view falls through to a scalar kernel
    cin = similar(pc.U, n)
    cout = similar(pc.U, n)
    for j in 1:k
        copyto!(cin, view(pc.U, :, j))
        applypreconditioner!(cout, pc.inner, cin)
        copyto!(view(Z, :, j), cout)
    end
    for j in 1:k
        copyto!(cin, view(Z, :, j))
        pc.jvp!(cout, cin)
        copyto!(view(Y, :, j), cout)
    end
    if k > pc.kmax
        # keep the subspace the operator shrinks most
        C = _hostbuilt(pc.U,
            eigen(Symmetric(Array(Y'*Y))).vectors[:, 1:pc.kmax])
        pc.U = pc.U*C; Z = Z*C; Y = Y*C; k = pc.kmax
    end
    F = svd(Array(Z'*Y))
    tol = eps(T)^(3//4)*maximum(F.S; init = zero(T))
    pc.Ginv = _hostbuilt(pc.U,
        F.V*Diagonal([s > tol ? inv(s) : zero(T) for s in F.S])*F.U')
    B = similar(Y)
    for j in 1:k
        copyto!(cin, view(Y, :, j))
        applypreconditioner!(cout, pc.inner, cin)
        copyto!(view(B, :, j), cout)
    end
    pc.Z = Z
    pc.W = Z .- B
    pc.a = similar(pc.U, k); pc.b = similar(pc.U, k)
    pc.rebuilds += 1
    return pc
end

function applypreconditioner!(z::AbstractVector, pc::RecyclingPreconditioner,
    r::AbstractVector)
    applypreconditioner!(z, pc.inner, r)
    if size(pc.Z, 2) > 0
        mul!(pc.a, transpose(pc.Z), r)
        mul!(pc.b, pc.Ginv, pc.a)
        mul!(z, pc.W, pc.b, true, true)
    end
    return z
end


# Orthonormalize the columns of `A` in place by CholeskyQR2, returning the
# number of columns retained or -1 if the factorization broke down.
#
# CholeskyQR2 is `G = A'A; R = chol(G); A = A/R`, done twice. Everything is
# level 3 BLAS (a syrk and a trsm per pass) and there are no column-by-column
# dependencies, which is what makes it the orthogonalization of choice on a
# GPU; Householder QR is sequential in the columns and, worse here, forming
# `Matrix(qr(A).Q)` materializes a dense n by k matrix on every harvest.
# The price is a squared condition number, so the Cholesky can fail on an
# ill conditioned block. That is not an edge case for a recycling subspace,
# whose whole purpose is to keep rediscovering the same few directions, so
# the caller projects and drops dependent columns first and falls back to a
# Householder QR when this returns -1.
function _choleskyqr2!(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    k = size(A, 2)
    k == 0 && return 0
    for _ in 1:2
        # The gram matrix is formed where `A` lives and factorized on the
        # host. The k x k triangular factor is inverted there and applied as a
        # gemm rather than sent back for a triangular solve, because a
        # `rdiv!` against a device resident `UpperTriangular` falls through to
        # a generic scalar kernel. Inverting is the less stable of the two,
        # which costs nothing here: the second CholeskyQR pass is what fixes
        # the squared condition number either way, and a block too ill
        # conditioned for that already returns -1 to the Householder fallback.
        F = cholesky(Symmetric(Array(A'*A)), check = false)
        issuccess(F) || return -1
        Rinv = _hostbuilt(A, Matrix(inv(F.U)))
        copyto!(A, A*Rinv)
    end
    return k
end

# Append `Unew` to the orthonormal basis `U`, keeping the result orthonormal.
# The projection is block classical Gram-Schmidt with one reorthogonalization,
# again all level 3 BLAS. Columns whose norm collapses under the projection
# were already in the span of `U` and are dropped rather than orthogonalized
# into noise.
function _orthappend(U::AbstractMatrix{T},
    Unew::AbstractMatrix{T}) where {T<:AbstractFloat}
    before = [norm(view(Unew, :, j)) for j in axes(Unew, 2)]
    if size(U, 2) > 0
        for _ in 1:2
            mul!(Unew, U, U'*Unew, -one(T), one(T))
        end
    end
    # A column is dropped when the projection removed essentially all of it,
    # measured against its own norm *before* projection. Comparing against
    # the largest surviving norm instead would retain a block that is pure
    # roundoff, which is exactly what a subspace that keeps rediscovering the
    # same directions produces.
    kept = [j for j in axes(Unew, 2)
        if norm(view(Unew, :, j)) > sqrt(eps(T))*before[j]]
    isempty(kept) && return U
    B = Unew[:, kept]
    if _choleskyqr2!(B) < 0
        # the Householder fallback is host only, and is reached only when the
        # block was too ill conditioned for CholeskyQR2
        B = _hostbuilt(U, Matrix(qr(Array(B)).Q)[:, 1:length(kept)])
    end
    return size(U, 2) > 0 ? hcat(U, B) : B
end

function harvest!(pc::RecyclingPreconditioner{TI,TJ,T}, ws::GMRESWorkspace,
    out::NamedTuple) where {TI,TJ,T}
    # nothing to deflate against an exact inner; see updatepreconditioner!
    isexactpreconditioner(pc.inner) && return pc
    # the workspace holds only the final cycle, so the usable Arnoldi
    # dimension is what that cycle built, not the total across restarts
    m = size(ws.H, 2)
    j = Int(out.iterations) - (Int(out.cycles) - 1)*m
    (1 <= j <= m && j >= pc.kharvest) || return pc
    # the directions this Krylov space found the operator shrinks most are the
    # right singular vectors of the rectangular Hessenberg with the smallest
    # singular values, lifted through the Arnoldi basis
    F = svd(Matrix(view(ws.H, 1:j+1, 1:j)))
    # `ws.H` is host resident, so the right singular vectors are built there
    # and moved to the basis before lifting them through it
    C = _hostbuilt(ws.V, F.V[:, size(F.V, 2)-pc.kharvest+1:size(F.V, 2)])
    Unew = view(ws.V, :, 1:j)*C
    U = _orthappend(pc.U, Unew)
    # guard only; the rebuild trims to kmax by quality on the next update
    pc.U = size(U, 2) > pc.kmax + pc.kharvest ?
        U[:, size(U, 2)-(pc.kmax + pc.kharvest)+1:end] : U
    return pc
end

"""
    norm2(v::AbstractVector)

The Euclidean norm of `v`, formed through the inner product.

`LinearAlgebra.norm` scales its argument before squaring so that an entry
cannot overflow or underflow on its way to the sum. That guard is not free
on a device: cuBLAS routes a Float64 vector to a scaled `nrm2` kernel which
runs at a small fraction of memory bandwidth, measured at 22.9 us against
2.7 us for a `dot` product of the same 51,200 element vector on an RTX 4090,
and it is the largest single device kernel of a double precision solve. In
single precision cuBLAS selects a different kernel and the gap is gone.

The substitution is made only where it pays and only where it is safe, which
is the same place: `Float64`.

- In double precision cuBLAS routes `norm` to a scaled `nrm2` kernel that runs
  at well under a tenth of memory bandwidth, measured at 55.6 us against
  14.5 us for `sqrt(dot(v, v))` on a 51,200 element device vector, and it was
  the largest single device kernel of a double precision solve.
- In single precision cuBLAS selects a different kernel and the two are within
  20% of each other, so there is nothing to win. Single precision is also
  where the exponent range is narrowest and the guard is worth the most.

So `Float32` keeps `norm` and everything else goes through the inner product,
and the change is confined to the precision where the trade is favorable in
both directions. Anything not covered by the `AbstractFloat` method -- complex
vectors, other element types -- falls back to `norm` as well.

The vectors this is applied to are carried at the scale of the harmonic
balance residual of a nondimensionalized system, which runs from a few units
down to the solver tolerance. Squaring that stays far inside the double
precision exponent: overflow would need an entry above 1e154 and underflow to
zero an entry below 1e-162.

This is deliberately not exported and not used outside the Krylov solver.
Anything whose scale is not controlled should keep using `norm`.
"""
norm2(v::AbstractVector{<:AbstractFloat}) = sqrt(dot(v, v))
norm2(v::AbstractVector{Float32}) = norm(v)
norm2(v::AbstractVector) = norm(v)

"""
    gmres_orthogonalize!(w, V, H, hd, c, j)

Orthogonalize `w` against the first `j` Arnoldi basis vectors (the columns of
`V`) by block classical Gram-Schmidt with one reorthogonalization (CGS2),
accumulating the coefficients into column `j` of the Hessenberg matrix `H`.
`c` is scratch of length at least `j`. Both passes accumulate into the same
entries of `H`, so `H` remains the exact projection. Writes the subdiagonal
`H[j+1, j]` and returns `(hsub, normw0)`: the norm of the orthogonalized `w`
and its norm on entry, the pair the caller compares to detect a breakdown.

CGS2 is chosen over modified Gram-Schmidt with a DGKS test for its shape
rather than its accuracy, which is equivalent: both are orthogonal to machine
precision. Each pass here is two level 2 BLAS calls over the whole basis,
where the modified form is `j` dependent pairs of a dot product and an axpy,
each dot having to complete before the axpy that follows it. The second pass
is unconditional, which costs what the DGKS path costs whenever it does
reorthogonalize and removes a branch on a freshly computed scalar. Neither the
coefficient vector nor the branch has to reach the host, which is what makes
this form usable on a device.
"""
function gmres_orthogonalize!(w::AbstractVector{T}, V::AbstractMatrix{T},
    H::AbstractMatrix{T}, hd::AbstractVector{T}, c::AbstractVector{T},
    j::Integer) where {T<:AbstractFloat}
    normw0 = norm2(w)
    if j > 0
        Vj = view(V, :, 1:j)
        # the projections are formed in buffers allocated like `V`, so on a
        # device they stay there; only the finished column crosses to `H`
        hj = view(hd, 1:j)
        cj = view(c, 1:j)
        # Block classical Gram-Schmidt, two unconditional passes. CGS2 is as
        # accurate as modified Gram-Schmidt with a DGKS test, to machine
        # precision, and it is a much better shape for a GPU: each pass is two
        # level 2 BLAS calls over the whole basis rather than `j` dependent
        # pairs of a dot product and an axpy. In the modified form every dot
        # has to finish before the axpy that follows it, which on a device
        # means `j` synchronizations per Arnoldi step; here the coefficient
        # vector never has to reach the host. Making the second pass
        # unconditional also removes a branch on a device-resident scalar,
        # which would force a synchronization of its own. On the CPU the two
        # passes cost what the DGKS path costs when it does reorthogonalize.
        mul!(hj, transpose(Vj), w)
        mul!(w, Vj, hj, -one(T), one(T))
        mul!(cj, transpose(Vj), w)
        mul!(w, Vj, cj, -one(T), one(T))
        hj .+= cj
        # `H` is deliberately host resident. copyto! between a host view and a
        # device view falls back to scalar indexing of the device array, so the
        # column crosses through its contiguous linear range instead: a dense
        # `H` stores rows 1:j of column j contiguously from (j-1)*size(H,1)+1.
        copyto!(H, (j-1)*size(H,1)+1, hd, 1, j)
    end
    hsub = norm2(w)
    H[j+1, j] = hsub
    return hsub, normw0
end

"""
    gmres_givens(a, b)

The Givens rotation `(c, s, r)` with `c*a + s*b = r` and `-s*a + c*b = 0`,
computed through `hypot` so it cannot overflow, with the identity rotation
returned for the zero input.
"""
function gmres_givens(a::T, b::T) where {T<:AbstractFloat}
    r = hypot(a, b)
    iszero(r) && return one(T), zero(T), zero(T)
    return a/r, b/r, r
end

"""
    gmres_applyrotations!(H, cs, sn, s, j)

Reduce column `j` of the Hessenberg matrix `H` to upper triangular form:
apply the `j-1` previous Givens rotations to the new column, compute and
store the rotation which annihilates the new subdiagonal `H[j+1, j]`, and
apply it to the least squares right hand side `s`. After this the magnitude
of `s[j+1]` is the residual norm of the least squares problem, which with
right preconditioning is the true residual norm of the original system.
Allocation free.
"""
function gmres_applyrotations!(H::AbstractMatrix{T}, cs::AbstractVector{T},
    sn::AbstractVector{T}, s::AbstractVector{T}, j::Integer) where {T<:AbstractFloat}
    for i in 1:j-1
        tmp       =  cs[i]*H[i, j] + sn[i]*H[i+1, j]
        H[i+1, j] = -sn[i]*H[i, j] + cs[i]*H[i+1, j]
        H[i, j]   = tmp
    end
    cs[j], sn[j], r = gmres_givens(H[j, j], H[j+1, j])
    H[j, j]   = r
    H[j+1, j] = zero(T)
    s[j+1] = -sn[j]*s[j]
    s[j]   =  cs[j]*s[j]
    return abs(s[j+1])
end

"""
    gmres_correction!(x, ws::GMRESWorkspace, j, Mop!)

Solve the reduced `j x j` triangular least squares problem by back
substitution, assemble the correction `u = V[:, 1:j]*y` in the Krylov basis,
undo the right preconditioning once with `Mop!` (or not at all when
`Mop! === nothing`), and add the result to `x` in place. A zero diagonal
entry, which can only arise from an exact breakdown, contributes a zero
coefficient rather than a division by zero. Allocation free.
"""
function gmres_correction!(x::AbstractVector{T}, ws::GMRESWorkspace{T},
    j::Integer, Mop!) where {T<:AbstractFloat}
    H, s, y, V, u, z = ws.H, ws.s, ws.y, ws.V, ws.u, ws.z

    # The projected problem can be rank deficient, on a singular or
    # inconsistent system or after a near breakdown. Back substituting
    # through a tiny pivot then amplifies roundoff without bound: a pivot of
    # order eps relative to the rest of the triangle produces a coefficient
    # of order 1/eps, and the returned "solution" is dominated by a spurious
    # basis direction. (Testing exactly for a zero pivot does not help, since
    # the damaging case is a pivot which is small but nonzero.) Instead the
    # numerically well determined leading part of the triangle is solved and
    # the remaining coefficients are set to zero, which is the minimizer over
    # the subspace the cycle actually resolved.
    rank = j
    dmax = zero(T)
    for i in 1:j
        dmax = max(dmax, abs(H[i, i]))
    end
    pivottol = eps(T)*max(dmax, one(T))*j
    for i in 1:j
        if abs(H[i, i]) <= pivottol
            rank = i - 1
            break
        end
    end
    for i in rank+1:j
        y[i] = zero(T)
    end
    for i in rank:-1:1
        acc = s[i]
        for k in i+1:rank
            acc -= H[i, k]*y[k]
        end
        y[i] = acc/H[i, i]
    end
    j = rank
    fill!(u, zero(T))
    for i in 1:j
        axpy!(y[i], view(V, :, i), u)
    end
    if isnothing(Mop!)
        @. x += u
    else
        Mop!(z, u)
        @. x += z
    end
    return x
end

"""
    gmres!(x, Aop!, b, ws::GMRESWorkspace; Mop! = nothing, rtol = 1e-6,
        atol = 0.0, maxrestarts = 10, initialzero = true)

Solve `A*x = b` with restarted GMRES, where `mul!(w, Aop, v)` computes `w = A*v` and
the optional `Mop!(z, v)` applies a preconditioner `z = M \\ v`. The matrix `A`
is never formed; only its action is required, which is what makes this usable
with the matrix-free [`jacobianvectorproduct!`](@ref).

Preconditioning is applied on the right, solving `A*inv(M)*u = b` and then
`x = inv(M)*u`. Right preconditioning keeps the recurrence's residual estimate
equal to the true residual of the original system, so the stopping test is on
`norm(b - A*x)` and does not depend on the quality of `M`. Because `M` is held
fixed across a solve, the preconditioner is applied once per Arnoldi step and
once more per restart, rather than being stored for every basis vector.

The Arnoldi basis is built by modified Gram-Schmidt with a conditional second
pass ([`gmres_orthogonalize!`](@ref)). A subdiagonal which collapses relative
to the vector it came from is a (lucky) breakdown: the Krylov space is
invariant, the reduced least squares solution is exact, and the cycle ends
there rather than continuing with a spurious basis vector. The residual is
recomputed explicitly at every restart so restarts cannot drift from the
recurrence estimate.

Converges when `norm(b - A*x) <= max(rtol*norm(b), atol)`. Returns the named
tuple `(iterations, residual, converged, cycles, reason)`, where `iterations`
counts Arnoldi steps across all cycles, `cycles` the number of restart cycles
begun, and `reason` is one of `:converged`, `:breakdown` (an unhappy
breakdown: the Krylov space went invariant without the residual coming down),
`:stagnation` (a cycle failed to reduce the explicit residual, or produced a
non-finite one), or `:iterationlimit`.

`iterations` is *not* the total number of `Aop!` calls: each cycle costs one
further application for the explicit residual recomputation, and a warm start
costs one at the outset. `maxrestarts` bounds the number of cycles including
the first, so the Arnoldi work is capped at `maxrestarts*m` steps.

Allocation free after the workspace is built, apart from whatever `Aop!` and
`Mop!` themselves allocate.
"""
function gmres!(x::AbstractVector{T}, Aop_, b::AbstractVector{T},
    ws::GMRESWorkspace{T}; Mop! = nothing, rtol = 1e-6, atol = 0.0,
    maxrestarts::Integer = 10, initialzero::Bool = true) where {T<:AbstractFloat}

    n = length(b)
    # a bare in-place product is accepted alongside any `mul!`-able operator
    Aop = asoperator(Aop_, n)
    precondtime = 0.0
    length(x) == n || throw(DimensionMismatch(
        lazy"`x` has length $(length(x)) but `b` has length $(n)."))
    size(ws.V, 1) == n || throw(DimensionMismatch(
        lazy"the workspace is for dimension $(size(ws.V,1)) but `b` has length $(n)."))
    (rtol >= 0 && isfinite(rtol)) || throw(ArgumentError(
        lazy"`rtol` = $(rtol) must be nonnegative and finite."))
    (atol >= 0 && isfinite(atol)) || throw(ArgumentError(
        lazy"`atol` = $(atol) must be nonnegative and finite."))
    maxrestarts >= 1 || throw(ArgumentError(
        lazy"`maxrestarts` = $(maxrestarts) must be at least 1."))

    m = size(ws.H, 2)
    V, H, cs, sn, s = ws.V, ws.H, ws.cs, ws.sn, ws.s
    w, z = ws.w, ws.z

    bnorm = norm(b)
    tol = max(rtol*bnorm, atol)

    # a zero right hand side has the zero solution; return it rather than
    # dividing by a zero residual norm below
    if iszero(bnorm)
        # zero is a solution, but so is any point already in the null space
        # of A, so a warm start is measured rather than discarded
        if initialzero
            fill!(x, zero(T))
            return (iterations = 0, residual = zero(T), converged = true,
                reason = :converged)
        end
        mul!(w, Aop, x)
        resnorm = norm(w)
        if resnorm <= atol
            return (iterations = 0, residual = resnorm, converged = true,
                reason = :converged)
        end
    end

    # initial residual w = b - A*x
    if initialzero
        fill!(x, zero(T))
        copyto!(w, b)
    else
        mul!(w, Aop, x)
        @. w = b - w
    end
    resnorm = norm(w)

    totaliterations = 0
    cycles = 0
    unhappy = false
    stagnated = false
    for _ in 1:maxrestarts
        resnorm <= tol && break

        cycles += 1
        beta = resnorm
        @views V[:, 1] .= w ./ beta
        fill!(s, zero(T))
        s[1] = beta

        j = 0
        while j < m
            j += 1
            # the Arnoldi step on A*inv(M)
            if isnothing(Mop!)
                @views copyto!(z, V[:, j])
            else
                tpc = time()
                @views Mop!(z, V[:, j])
                precondtime += time() - tpc
            end
            mul!(w, Aop, z)

            hsub, normw0 = gmres_orthogonalize!(w, V, H, ws.hd, ws.cd, j)
            resnorm = gmres_applyrotations!(H, cs, sn, s, j)
            totaliterations += 1

            # A subdiagonal which collapsed relative to the incoming vector
            # means the Krylov space is invariant and there is no valid next
            # basis vector to normalize. This is not automatically a *lucky*
            # breakdown: that additionally requires the projected problem to
            # be compatible, so that the reduced solution really is the
            # solution. On a singular or inconsistent system the space is
            # invariant while the residual is not reducible in it at all, an
            # unhappy breakdown, and continuing to restart only rebuilds the
            # same useless space. The two are separated by whether the
            # recurrence residual actually came down.
            breakdown = hsub <= eps(T)*normw0
            if breakdown
                unhappy = resnorm > tol && resnorm > (1 - sqrt(eps(T)))*beta
                break
            end
            resnorm <= tol && break

            @views V[:, j+1] .= w ./ hsub
        end

        gmres_correction!(x, ws, j, Mop!)

        # recompute the residual explicitly for the next cycle so restarts
        # cannot drift from the recurrence estimate
        mul!(w, Aop, x)
        @. w = b - w
        previous = beta
        resnorm = norm(w)

        # a cycle which did not reduce the explicit residual will not do
        # better on a rebuild of the same space, so stop rather than burn the
        # remaining cycle budget on it
        if !isfinite(resnorm)
            stagnated = true
            break
        end
        if unhappy || resnorm >= (1 - sqrt(eps(T)))*previous
            stagnated = true
            break
        end
    end

    converged = resnorm <= tol
    reason = if converged
        :converged
    elseif unhappy
        :breakdown
    elseif stagnated
        :stagnation
    else
        :iterationlimit
    end
    return (iterations = totaliterations, residual = resnorm,
        converged = converged, cycles = cycles, reason = reason,
        precondtime = precondtime)
end

abstract type AbstractHBLinearSolver end

"""
    InternalGMRES()

The restarted GMRES of this package, with Givens rotations, the recycling
subspace harvest and the preconditioner escalation the solver was built
around. The default, and the only solver supporting deflation recycling,
because `harvest!` reads the Arnoldi basis out of the internal workspace.
"""
struct InternalGMRES <: AbstractHBLinearSolver end

"""
    KrylovJL(method::Symbol = :gmres; kwargs...)

Solve the Newton step with Krylov.jl: `:gmres`, `:fgmres`, `:bicgstab`,
`:dqgmres` or any other solver taking an operator and a right hand side.
Requires Krylov.jl to be loaded; the method lives in the package extension.

Only the linear solve changes. The forcing term, the line search, the
preconditioner escalation and the stagnation handling are untouched, and
the mode coupling preconditioner is passed through unchanged because it is
applied by `ldiv!`, which is what Krylov.jl's `N` argument consumes.
Deflation recycling is unavailable, since it depends on the internal
workspace.
"""
struct KrylovJL{K} <: AbstractHBLinearSolver
    method::Symbol
    kwargs::K
end
KrylovJL(method::Symbol = :gmres; kwargs...) = KrylovJL(method, kwargs)

"""
    hblinearsolve!(ls, deltax, jvp!, F, ws, Mop!; rtol, atol, maxrestarts)

Solve for the Newton step and return the output named tuple `gmres!`
produces: `converged`, `residual`, `iterations`, `cycles`, `reason`.
`jvp!(y, v)` applies the Jacobian and `Mop!(z, r)` the preconditioner, both
in place.

This is the one part of the Newton-Krylov loop with nothing harmonic
balance specific about it: an operator, a right hand side, a preconditioner
and a tolerance. Putting it behind an interface lets an external Krylov
library be substituted without touching anything else.
"""
function hblinearsolve!(::InternalGMRES, deltax, jvp!, F, ws, Mop!;
        rtol, atol, maxrestarts)
    return gmres!(deltax, jvp!, F, ws; Mop! = Mop!, rtol = rtol,
        atol = atol, maxrestarts = maxrestarts)
end

hblinearsolve!(ls::AbstractHBLinearSolver, args...; kwargs...) =
    throw(ArgumentError(lazy"no hblinearsolve! method for $(typeof(ls)); "*
        "KrylovJL requires Krylov.jl to be loaded."))

"""
    supportsrecycling(ls)

Whether the solver exposes an Arnoldi basis for [`harvest!`](@ref).
"""
supportsrecycling(::AbstractHBLinearSolver) = false
supportsrecycling(::InternalGMRES) = true

"""
    nlsolvekrylov!(fj!, jvp!, F, x, pc::AbstractPreconditioner; kwargs...)

Inexact (Newton-Krylov) solver for a real system: the Newton step is taken
from [`gmres!`](@ref) on the exact matrix-free product `jvp!(y, v)` rather
than from a factorization of an assembled Jacobian. `fj!(F, J, x)` evaluates
the residual and Jacobian as in [`nlsolve!`](@ref), accepting `nothing` for
either. `x` is updated in place and `F` holds the residual at the returned
`x`.

The Jacobian enters only through `jvp!` and through the right preconditioner
`pc` (see [`AbstractPreconditioner`](@ref)), so a preconditioner much cheaper
than the Jacobian makes each Newton step much cheaper than a direct one.
[`ModeCouplingPreconditioner`](@ref) is the harmonic balance Jacobian with its
mode coupling restricted, and is exact when every mode is retained.

`pc` is rebuilt lazily: when the GMRES iteration count says it has
gone stale (`krylovrefreshiterations`), when GMRES fails, when a step is not
a descent direction, or when a linesearch cannot find any decrease. The
linear tolerance follows the Eisenstat-Walker choice 2 forcing sequence
`krylovgamma*(|F_k|/|F_{k-1}|)^krylovalpha` clamped to
`[krylovrtolmin, krylovrtolmax]`, with an absolute floor of `ftol/10` so late
solves are not pushed below the nonlinear tolerance. Because the assembled
Jacobian can be stale, the linesearch slope is always taken from an exact
matrix-free product, and a non-descent direction falls back to the exact
Newton step through a fresh factorization before the iteration is declared
stalled.

The globalization is the plain damped-Newton path of [`nlsolve!`](@ref):
the interpolated [`backtracking_linesearch!`](@ref), which on Armijo failure
still takes the best decreasing trial, with consecutive failures counted
against `maxbacktrackfailures` and a no-decrease step retried once from a
fresh preconditioner before stopping. There is deliberately no Anderson
acceleration here: the Krylov steps are near-exact Newton steps, and the
solver is kept simple.

# Keywords
- `iterations = 1000`: maximum number of Newton iterations.
- `ftol = 1e-8`: convergence when `norm(F) <= ftol`.
- `label = ""`: label for the returned `IterationInfo`.
- `c1 = 1e-4`: Armijo sufficient-decrease constant, in (0, 1/2).
- `safeguard_low = 0.1`, `safeguard_high = 0.5`: backtracking step clamp as
  fractions of the previous trial.
- `maxbacktracks = 10`: trial-point budget per line search.
- `maxbacktrackfailures = 2`: consecutive-failure stall threshold.
- `krylovrestart = 400`: GMRES restart length. A long cycle is the tuned
  default: the block diagonal preconditioner leaves a few near-null
  directions a short Krylov space cannot resolve, the linear solve stalls,
  and escalation to the full Jacobian is what rescues it -- a long cycle
  attacks the stall directly and the basis (`krylovrestart + 1` vectors of
  the system dimension) is cheap next to the sparse factorization
  escalation would build.
- `krylovmaxrestarts = 4`: GMRES restart budget per solve.
- `krylovrefreshiterations = 1`: GMRES iteration count above which the
  preconditioner is considered stale and refreshed before the next step;
  the tuned default refreshes eagerly.
- `krylovrtolmin = 1e-10`, `krylovrtolmax = 0.9`, `krylovrtol0 = 0.3`:
  clamp and initial value of the forcing sequence.
- `krylovgamma = 0.9`, `krylovalpha = (1 + sqrt(5))/2`: Eisenstat-Walker
  forcing parameters.

These are the values `hbnlsolve` runs with: the signature is the single
source of truth for the tuned defaults, so a direct caller gets the same
solver the package uses.

Returns an [`IterationInfo`](@ref) with the same per-iteration diagnostics
as [`nlsolve!`](@ref); the `andersonaccepted` record is always false.
"""
function nlsolvekrylov!(fj!::Function, jvp!, F::AbstractVector{T},
    x::AbstractVector{T}, pc::AbstractPreconditioner; iterations = 1000, ftol = 1e-8,
    linearsolver::AbstractHBLinearSolver = InternalGMRES(),
    label = "",
    c1 = 1e-4, safeguard_low = 0.1, safeguard_high = 0.5,
    maxbacktracks::Integer = 10, maxbacktrackfailures::Integer = 2,
    krylovrestart::Integer = 400, krylovmaxrestarts::Integer = 4,
    krylovrefreshiterations::Integer = 1, krylovrtolmin = 1e-10,
    krylovrtolmax = 0.9, krylovrtol0 = 0.3, krylovgamma = 0.9,
    krylovalpha = (1 + sqrt(5))/2,
    krylovstagnation = 0.9, krylovescalate::Integer = 1,
    krylovrefreshrate = 0.5) where {T<:AbstractFloat}

    length(F) == length(x) || throw(DimensionMismatch(
        lazy"The residual `F` has length $(length(F)) but the point `x` has length $(length(x))."))

    # a bare in-place product is normalized to a `mul!`-able operator once,
    # so the loop below and the pluggable linear solver see one interface
    jvp = asoperator(jvp!, length(x))

    # validate every option before the first residual evaluation, with the
    # same bounds and rationale as nlsolve!
    iterations >= 0 || throw(ArgumentError(
        lazy"`iterations` = $(iterations) must be nonnegative."))
    ftol >= 0 || throw(ArgumentError(lazy"`ftol` = $(ftol) must be nonnegative."))
    0 < c1 < 1//2 || throw(ArgumentError(
        lazy"`c1` = $(c1) must be in (0, 1/2) for the Newton merit function."))
    0 < safeguard_low < 1//2 || throw(ArgumentError(
        lazy"`safeguard_low` = $(safeguard_low) must be in (0, 1/2)."))
    safeguard_low < safeguard_high < 1 || throw(ArgumentError(
        lazy"`safeguard_high` = $(safeguard_high) must satisfy `safeguard_low < safeguard_high < 1`."))
    maxbacktracks >= 0 || throw(ArgumentError(
        lazy"`maxbacktracks` = $(maxbacktracks) must be nonnegative."))
    maxbacktrackfailures >= 1 || throw(ArgumentError(
        lazy"`maxbacktrackfailures` = $(maxbacktrackfailures) must be positive."))
    krylovrestart >= 1 || throw(ArgumentError(
        lazy"`krylovrestart` = $(krylovrestart) must be at least 1."))
    krylovmaxrestarts >= 1 || throw(ArgumentError(
        lazy"`krylovmaxrestarts` = $(krylovmaxrestarts) must be at least 1."))
    krylovrefreshiterations >= 1 || throw(ArgumentError(
        lazy"`krylovrefreshiterations` = $(krylovrefreshiterations) must be at least 1."))
    0 < krylovrtolmin <= krylovrtolmax < 1 || throw(ArgumentError(
        lazy"`krylovrtolmin` = $(krylovrtolmin) and `krylovrtolmax` = $(krylovrtolmax) must satisfy 0 < krylovrtolmin <= krylovrtolmax < 1."))
    krylovrtolmin <= krylovrtol0 <= krylovrtolmax || throw(ArgumentError(
        lazy"`krylovrtol0` = $(krylovrtol0) must satisfy `krylovrtolmin` <= `krylovrtol0` <= `krylovrtolmax`."))
    # the Choice 2 parameters are only defined on these ranges; outside them
    # the forcing sequence is not a forcing sequence
    (isfinite(krylovgamma) && 0 < krylovgamma <= 1) || throw(ArgumentError(
        lazy"`krylovgamma` = $(krylovgamma) must be finite and in (0, 1]."))
    (isfinite(krylovalpha) && 1 < krylovalpha <= 2) || throw(ArgumentError(
        lazy"`krylovalpha` = $(krylovalpha) must be finite and in (1, 2]."))
    0 < krylovstagnation <= 1 || throw(ArgumentError(
        lazy"`krylovstagnation` = $(krylovstagnation) must be in (0, 1]."))
    krylovescalate >= 1 || throw(ArgumentError(
        lazy"`krylovescalate` = $(krylovescalate) must be at least 1."))
    0 < krylovrefreshrate <= 1 || throw(ArgumentError(
        lazy"`krylovrefreshrate` = $(krylovrefreshrate) must be in (0, 1]."))

    ws = GMRESWorkspace(x, min(krylovrestart, length(x)))
    deltax = similar(x)
    xcandidate = similar(x)
    Jv = similar(x)
    Fbest = similar(F)

    # absolute floor for the linear solves: once the linear residual is below
    # the nonlinear tolerance, further accuracy cannot help the Newton
    # iteration, and demanding it makes late GMRES solves "fail"
    gmresatol = real(T)(ftol)/10

    ### diagnostic info
    krylovrecord = KrylovSolveInfo[]
    backtrackrecord = Int[]
    tstart = time()
    alphas = real(T)[]
    normF = real(T)[]
    converged = false
    backtrackfailures = 0
    refresh = true
    refreshedforstall = false
    # consecutive linear solves which failed to reach the forcing tolerance,
    # which trigger an escalation of the preconditioner
    linearfailures = 0
    # Slope-aware safeguard on the forcing sequence. The Eisenstat-Walker
    # term reads only the residual reduction, which makes a crawl
    # self-perpetuating: a loose solve gives a weak direction, the line
    # search accepts a short step, the residual barely falls, and the loose
    # forcing is reproduced -- measured as tens of outer iterations at
    # alpha ~ 0.2-0.4 with normalized slopes of -0.5 to -0.7 on a line of
    # scattering block capacitors. The inexact Newton bound
    # |slope| >= 1 - eta ties the direction quality to the forcing
    # directly, so when a step comes back weak (a backtracked line search
    # or a shallow slope) the cap tightens by a factor of four, and it
    # relaxes by a factor of two once full steps resume. Well behaved
    # problems take full steps at slope ~ -1 from the start and never feel
    # the cap.
    forcingcap = krylovrtolmax

    Mop!(zv, vv) = applypreconditioner!(zv, pc, vv)
    # residual-only adapter for the linesearch, which never needs the
    # Jacobian and therefore does not accept the combined fj! interface
    residual!(Fv, xv) = fj!(Fv, nothing, xv)

    # rebuild the preconditioner at the current point. a preconditioner is
    # free to move the evaluation point of the matrix-free products while
    # rebuilding, so it is resynchronized afterwards
    function refreshpreconditioner!()
        updatepreconditioner!(pc, x)
        fj!(nothing, nothing, x)
        return nothing
    end

    # the residual norm at the initial point; every later entry of normF is
    # pushed immediately after a step is accepted, so convergence is decided
    # on each fresh residual and no preconditioner is ever assembled at a
    # final point
    residual!(F, x)
    push!(normF, norm(F))
    normF[end] <= ftol && (converged = true)

    for n in 1:iterations
        converged && break

        # the matrix-free product reads the evaluation point held by the
        # caller, and the linesearch leaves it at the last trial point
        # rather than the accepted one, so resynchronize it
        fj!(nothing, nothing, x)
        justrefreshed = refresh
        if refresh
            refreshpreconditioner!()
            refresh = false
        end

        # Eisenstat-Walker choice 2 forcing term from the last accepted
        # step, at its clamp maximum before any step has been taken
        forcing = if length(normF) >= 2 && normF[end-1] > 0
            clamp(min(krylovgamma*(normF[end]/normF[end-1])^krylovalpha,
                    forcingcap),
                krylovrtolmin, krylovrtolmax)
        else
            # the *initial* forcing term, which is a separate quantity from
            # the upper clamp: seeding the safeguard at krylovrtolmax would
            # let it walk down from there over several outer steps
            # (0.9, 0.76, 0.58, 0.37, 0.18 for the default gamma and alpha),
            # so that several successive linear solves terminate after very
            # little residual reduction
            krylovrtol0
        end

        out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
            rtol = forcing, atol = gmresatol,
            maxrestarts = krylovmaxrestarts)
        supportsrecycling(linearsolver) && harvest!(pc, ws, out)
        # one record per GMRES call; the step outcome is filled in later
        function record!(o, role, refreshedbefore, stag)
            push!(krylovrecord, KrylovSolveInfo(n, role, normF[end], forcing,
                normF[end] > 0 ? o.residual/normF[end] : NaN, o.iterations,
                o.cycles, o.reason, refreshedbefore, false, stag, NaN, NaN,
                0, false, time() - tstart, false, deflationsize(pc),
                deflationrebuilds(pc), get(o, :precondtime, NaN)))
            return nothing
        end
        stagnated = !out.converged &&
            out.residual > krylovstagnation*normF[end]
        record!(out, :step, justrefreshed, stagnated)
        # `!justrefreshed` matters: a retry only makes sense against a
        # preconditioner which had drifted. If it was rebuilt at this very
        # point immediately before the solve, rebuilding it again reproduces
        # the same operator and the retry reproduces the same failure at the
        # cost of a second full Krylov cycle.
        if !out.converged && !stagnated && !justrefreshed
            # GMRES made real progress but did not reach the tolerance,
            # which a preconditioner that has drifted can cause. Rebuild it
            # and retry once. A *stagnated* solve is not retried: stagnation
            # means the Krylov space contributed nothing, which is a
            # preconditioner too crude for the problem rather than one that
            # is merely stale.
            refreshpreconditioner!()
            out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
                rtol = forcing, atol = gmresatol,
                maxrestarts = krylovmaxrestarts)
            supportsrecycling(linearsolver) && harvest!(pc, ws, out)
            stagnated = !out.converged &&
                out.residual > krylovstagnation*normF[end]
            record!(out, :retry, true, stagnated)
        end
        # A GMRES which ran out of iterations is not automatically a failure
        # to be undone: it still returns the step which minimizes the linear
        # residual over the Krylov space it did build. A stagnated solve
        # produced nothing usable: its Krylov space is no better than the
        # preconditioner solve, which is the better step and, for an exact
        # preconditioner, the Newton step.
        if stagnated
            applypreconditioner!(deltax, pc, F)
        end
        # Escalation is triggered by repeated *non-convergence*, not by
        # stagnation. These are different symptoms: stagnation is a Krylov
        # space which contributed nothing, while the common failure of a
        # merely inadequate preconditioner is a solve which makes real
        # progress and still cannot reach the forcing tolerance within its
        # budget. Only the latter repeats, so keying escalation to stagnation
        # alone would leave it never firing on the problems it exists to
        # rescue.
        if out.converged
            linearfailures = 0
        else
            linearfailures += 1
            if linearfailures >= krylovescalate
                if escalatepreconditioner!(pc)
                    refresh = true
                    krylovrecord[end] = _withescalated(krylovrecord[end])
                else
                    krylovrecord[end] = _withescalationrefused(krylovrecord[end])
                end
                linearfailures = 0
            end
        end
        # Staleness is judged by how fast the linear residual came down per
        # Arnoldi step, not by the raw step count alone. The count on its own
        # is confounded by the requested accuracy: with a loose forcing term a
        # badly stale preconditioner can satisfy the tolerance in a single
        # iteration after removing only a small fraction of the linear
        # residual, and a rule keyed to the count reads that as health.
        if out.iterations >= krylovrefreshiterations
            refresh = true
        elseif out.iterations > 0 && normF[end] > 0
            rate = (out.residual/normF[end])^(1/out.iterations)
            isfinite(rate) && rate > krylovrefreshrate && (refresh = true)
        end
        rmul!(deltax, -1)

        # the merit function and its slope along deltax. the assembled
        # Jacobian is stale, so the slope comes from an exact matrix-free
        # product rather than the model claim dot(F, J, deltax) used by
        # nlsolve!
        ϕ0 = merit(F)
        mul!(Jv, jvp, deltax)
        dϕ0dα = real(dot(F, Jv))
        if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
            # not a descent direction: rebuild the preconditioner and solve
            # again. For an exact preconditioner GMRES then returns the
            # Newton step, which is a descent direction up to roundoff, so
            # this reproduces the exact-Newton rescue; for an inexact one it
            # is the best step in the Krylov space of a fresh operator,
            # which is the strongest direction available without assembling
            # the Jacobian. If that is not a descent direction either, the
            # steepest descent direction -J'F of the merit function is the
            # last resort before declaring the iteration stalled.
            refreshpreconditioner!()
            out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
                rtol = forcing, atol = gmresatol,
                maxrestarts = krylovmaxrestarts)
            record!(out, :rescue, true, false)
            rmul!(deltax, -1)
            mul!(Jv, jvp, deltax)
            dϕ0dα = real(dot(F, Jv))
            if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
                break
            end
        end

        # interpolated backtracking linesearch shared with nlsolve!: on
        # Armijo failure it returns the best decreasing trial (alpha > 0)
        # with F and xcandidate restored there, or alpha == 0 when no trial
        # decreased the merit at all
        alpha1, ϕα, accepted, backtracks = backtracking_linesearch!(
            residual!, F, xcandidate, x, deltax, ϕ0, dϕ0dα;
            c1 = c1, safeguard_low = safeguard_low,
            safeguard_high = safeguard_high,
            maxbacktracks = maxbacktracks, Fbest = Fbest)
        push!(alphas, alpha1)
        push!(backtrackrecord, backtracks)
        if !isempty(krylovrecord) && krylovrecord[end].iteration == n
            krylovrecord[end] = _withstep(krylovrecord[end],
                normF[end] > 0 ? dϕ0dα/normF[end]^2 : NaN, alpha1,
                backtracks, accepted)
        end

        if iszero(alpha1)
            # no decrease anywhere along the direction; F again holds the
            # residual at the unchanged x (the linesearch restore contract).
            # retry once from a fresh preconditioner, then give up
            refreshedforstall && break
            refreshedforstall = true
            refresh = true
            continue
        end
        refreshedforstall = false
        # the safeguard state advances only for a nonlinear step which was
        # actually taken. a retry at the same point is not a new outer
        # iteration and must not walk the forcing sequence forward
        forcingprev = forcing
        # step quality read by the slope-aware forcing cap: the normalized
        # slope of the merit along the step, and whether the line search
        # took the full step
        stepslope = normF[end] > 0 ? dϕ0dα/normF[end]^2 : -one(ϕ0)
        if alpha1 < one(alpha1) || stepslope > -0.9
            # the cap floor stays well above the forcing clamp minimum: a
            # near-exact solve costs a long Krylov cycle per step, and the
            # direction quality it buys beyond this point is marginal
            forcingcap = max(1e-2, forcingcap/4)
        else
            forcingcap = min(krylovrtolmax, 2*forcingcap)
        end

        # accept the trial point: F already holds the residual there (the
        # linesearch postcondition), and convergence is decided on it now,
        # before any preconditioner work
        copyto!(x, xcandidate)
        push!(normF, norm(F))
        if normF[end] <= ftol
            converged = true
            break
        end

        # count consecutive linesearch failures where the best decreasing
        # trial was taken instead of an Armijo-accepted step; an accepted
        # step resets the count, repeated failures are a stall
        if accepted
            backtrackfailures = 0
        else
            backtrackfailures += 1
            backtrackfailures >= maxbacktrackfailures && break
        end
    end

    return IterationInfo(label, NaN, 0.0, converged, length(alphas), normF,
        alphas, backtrackrecord, fill(false, length(alphas)), krylovrecord)
end

