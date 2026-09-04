
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
- `deflationsize`, `deflationrebuilds`, `precondtime`: the active rank of
    the recycled deflation, how many times it has been built, and the wall
    time spent applying the preconditioner in this solve.
- `products`, `deflationproducts`: the exact operator products this linear
    solve took, and the running count of those the recycling wrapper took
    for its builds and standalone applications. The cost of a solve is in
    these, not in `iterations` alone: every restart cycle recomputes the
    residual, and `:adef2` reapplies its product at the restart correction.
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
    # whether an escalation was requested for this solve; with `escalated`
    # this distinguishes an escalation never requested from one requested
    # and refused
    escalationrequested::Bool
    # the width of the deflation subspace and how many times it has been
    # rebuilt, which is what recycling costs
    deflationsize::Int
    deflationrebuilds::Int
    # wall time applying the preconditioner, which tells a weaker and cheaper
    # preconditioner from a worse one when the iteration counts look alike
    precondtime::Float64
    # exact operator products: those the linear solve took (Arnoldi steps,
    # the residual recomputed at every restart, a warm start) and, in the
    # running count `deflationproducts`, those the recycling wrapper took
    # for its builds and standalone `:adef2` applications; Arnoldi steps
    # alone understate the cost of a form or a restart
    products::Int
    deflationproducts::Int
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
        " products=", k.products,
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


# the step outcome is only known after the linesearch, so the record is
# completed then
function _withstep(k::KrylovSolveInfo, slope, alpha, backtracks, armijo)
    return KrylovSolveInfo(k.iteration, k.role, k.normF, k.forcing,
        k.residualratio, k.iterations, k.cycles, k.reason, k.refreshed,
        k.escalated, k.stagnated, slope, alpha, backtracks, armijo, k.time,
        k.escalationrequested, k.deflationsize, k.deflationrebuilds,
        k.precondtime, k.products, k.deflationproducts)
end

function _withescalated(k::KrylovSolveInfo)
    return KrylovSolveInfo(k.iteration, k.role, k.normF, k.forcing,
        k.residualratio, k.iterations, k.cycles, k.reason, k.refreshed,
        true, k.stagnated, k.slope, k.alpha, k.backtracks, k.armijo, k.time,
        k.escalationrequested, k.deflationsize, k.deflationrebuilds,
        k.precondtime, k.products, k.deflationproducts)
end

# an escalation which the preconditioner refused: the request happened, the
# strengthening did not
function _withescalationrefused(k::KrylovSolveInfo)
    return KrylovSolveInfo(k.iteration, k.role, k.normF, k.forcing,
        k.residualratio, k.iterations, k.cycles, k.reason, k.refreshed,
        k.escalated, k.stagnated, k.slope, k.alpha, k.backtracks, k.armijo,
        k.time, true, k.deflationsize, k.deflationrebuilds, k.precondtime,
        k.products, k.deflationproducts)
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
    GMRESWorkspace{T<:AbstractFloat}

Preallocated storage for [`gmres!`](@ref) with a restart length of `m` on a
system of dimension `n`. Holds the `n x (m+1)` Arnoldi basis `V`, the
`(m+1) x m` Hessenberg matrix `H` as the Givens rotations leave it, the raw
Arnoldi Hessenberg `Harnoldi` beside it, the Givens rotations `cs` and `sn` which
reduce it, the least squares right hand side `s`, its solution `y`, and three
length `n` work vectors.

The dominant cost is `V`, which is `n*(m+1)` numbers, so `m` trades memory and
orthogonalization work against restart frequency. It is not paid up front:
the basis is allocated with a few columns and grows, by doubling, to what
the iteration uses, up to `m + 1`. A restart length long enough for the
hardest solve is then free on the easy ones, where a warm started Newton
step takes a handful of Arnoldi steps, and the cost of a solve is no longer
dominated by touching a basis it never fills. See [`ensurecolumns!`](@ref).
"""
mutable struct GMRESWorkspace{T<:AbstractFloat,TV<:AbstractVector{T},TM<:AbstractMatrix{T}}
    # system sized, and therefore allocated like the right hand side, so on a
    # device backend they live on the device. `V` holds the columns built so
    # far and is replaced by a wider one as the iteration needs them; `H`
    # is host resident and small, and is allocated in full.
    V::TM
    w::TV
    z::TV
    u::TV
    # the projected quantities, at most `m` by `m`, kept on the host: the
    # Givens rotations and the back substitution index them entry by entry,
    # which a device array must not be asked to do, and they are small
    # enough that the cost is one small transfer per Arnoldi step
    H::Matrix{T}
    # the Arnoldi Hessenberg before the Givens rotations, column by column
    # as each is finished. `H` itself is the triangularized least squares
    # matrix once the rotations have been applied; a harvest that needs the
    # Arnoldi relation `A*V[:, 1:j] = V[:, 1:j+1]*Harnoldi[1:j+1, 1:j]`
    # (the harmonic Ritz pencil) reads this one. The singular values of the
    # two agree, since a left orthogonal transformation preserves them.
    Harnoldi::Matrix{T}
    cs::Vector{T}
    sn::Vector{T}
    s::Vector{T}
    y::Vector{T}
    # staging buffers for the block Gram-Schmidt projections, allocated like
    # `V`, so that only the finished column of `H` crosses to the host
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
    cols = min(m + 1, GMRESINITIALCOLUMNS)
    return GMRESWorkspace{T,typeof(similar(b)),typeof(similar(b, n, cols))}(
        similar(b, n, cols), similar(b), similar(b), similar(b),
        zeros(T, m + 1, m), zeros(T, m + 1, m),
        Vector{T}(undef, m), Vector{T}(undef, m),
        Vector{T}(undef, m + 1), Vector{T}(undef, m),
        similar(b, m), similar(b, m))
end

# the columns a basis is born with. Sixteen covers the Arnoldi steps of a
# well preconditioned or warm started Newton step without a copy, and a
# solve which needs more doubles its way there; at most log2(m/16) copies of
# what was built, which is a fraction of the orthogonalization that built it
const GMRESINITIALCOLUMNS = 16

"""
    ensurecolumns!(ws::GMRESWorkspace, k)

Return the basis `ws.V` with at least `k` columns, replacing it with a wider
one when it has fewer. The columns built so far are copied across; the new
ones are uninitialized, as the whole basis was before. The width doubles
and is capped at `size(ws.H, 1)`, which is `m + 1`, so a basis never grows
past the restart length.
"""
function ensurecolumns!(ws::GMRESWorkspace, k::Integer)
    V = ws.V
    have = size(V, 2)
    k <= have && return V
    cols = min(max(k, 2*have), size(ws.H, 1))
    Vnew = similar(V, size(V, 1), cols)
    copyto!(view(Vnew, :, 1:have), V)
    ws.V = Vnew
    return Vnew
end


"""
    KrylovVectors(x, F, m)

The system sized vectors of one Newton-Krylov solve: the GMRES workspace
for a restart length of `m`, and the step, the trial point, the product
and the best residual, allocated like `x` and `F`. Handed back to
[`nlsolvekrylov!`](@ref) through its `workspace` argument they are reused
across solves of one system, which a sweep over component values is.
"""
struct KrylovVectors{V,VF,W}
    ws::W
    deltax::V
    xcandidate::V
    Jv::V
    Fbest::VF
end

KrylovVectors(x::AbstractVector, F::AbstractVector, m::Integer) =
    KrylovVectors(GMRESWorkspace(x, m), similar(x), similar(x), similar(x),
        similar(F))

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
    harvestcycle!(pc::AbstractPreconditioner, ws::GMRESWorkspace, j::Integer)

Give the preconditioner the Arnoldi factorization of the restart cycle which
has just ended, `j` vectors of it, before [`gmres!`](@ref) overwrites the
workspace with the next cycle. The default does nothing.

This is the per-cycle counterpart of [`harvest!`](@ref), which sees only the
cycle left in the workspace when the solve returns. A preconditioner opts
into it through [`usescycleharvest`](@ref), and one which does is *not*
harvested again afterwards.
"""
harvestcycle!(pc::AbstractPreconditioner, ::GMRESWorkspace, ::Integer) = pc

"""
    usescycleharvest(pc::AbstractPreconditioner)

Whether `pc` wants [`harvestcycle!`](@ref) at the end of every restart cycle
instead of [`harvest!`](@ref) once the solve is over. `false` by default, so
that a preconditioner which harvests only the final cycle keeps doing
exactly that.
"""
usescycleharvest(::AbstractPreconditioner) = false

"""
    isexactpreconditioner(pc::AbstractPreconditioner)

Whether `pc` currently applies the exact Jacobian, so that a deflation or
composition layered on top of it can contribute nothing. `false` for any
preconditioner that does not say otherwise.
"""
isexactpreconditioner(::AbstractPreconditioner) = false

"""
    RecyclingState(b::AbstractVector)

The part of a [`RecyclingPreconditioner`](@ref) that outlives one solve:
the candidate directions `U`, an `n` by `kc` matrix with orthonormal
columns allocated like `b`. A candidate is a direction of the residual on
which the preconditioned operator `J*inv(P)` was found to be small, the
metric the harvest measures in; the correction it stands for is
`inv(P)*u`, formed at the next build against whatever base is then in
place. [`HBReuse`](@ref) carries one across a cached sweep, where a nearby
operating point rebinds the base to nearly the same operator; the state is
replaced only after a converged solve, so a failed point cannot seed the
next one.
"""
mutable struct RecyclingState{TM<:AbstractMatrix}
    U::TM
end

RecyclingState(b::AbstractVector) = RecyclingState(similar(b, length(b), 0))
Base.copy(s::RecyclingState) = RecyclingState(copy(s.U))

"""
    RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!, n::Integer;
        kmax = 20, kharvest = 8, escalateafter = 1, form = :adef1)

Augments any base preconditioner with a recycled deflation subspace, carried
across Newton steps and, through a [`RecyclingState`](@ref), across the
solves of a sweep.

The base preconditioner `inner` leaves a handful of directions on which the
preconditioned operator `J*inv(P)` is nearly singular. Those directions are
what stalls GMRES, because its residual polynomial is pinned at `p(0) = 1` and
so cannot be small near the origin. Rather than enlarging `inner` until it
happens to cover them, this wrapper *measures* them: after each solve
[`harvest!`](@ref) extracts the directions the Arnoldi factorization found the
operator shrinks most, maps them to corrections of the unknowns, and the next
solve deflates them explicitly.

The coarse space is the residual-image pair of the candidates: with `U` the
candidate directions, `Z = inv(P)*U` their corrections, one Jacobian product
per column and one small singular value decomposition give

    X = Z*V/S,    C = J*X,    C'*C = I,

after which `Q = X*C'` satisfies `J*Q = C*C'`, an orthogonal projector in the
norm GMRES minimizes. This is the recycling space of GCRO-DR (Parks, de
Sturler, Mackey, Johnson and Maiti, SIAM J. Sci. Comput. 28, 2006). It is
preferred to the Galerkin pairing `Z'*J*Z` because a nonsymmetric Jacobian
can pair a useful direction with a zero there (`J = [0 -e; 1 0]`, `Z = e2`
gives `Z'*J*Z = 0` while `norm(J*Z) = e`), leaving the correction silently
absent; the image pair cannot, and it needs no projected inverse and no
tolerance for one. The deflated preconditioner is one of two forms, selected
by `form`, which differ in whether the projection is applied to the *input*
or to the *output* of the base solve:

- `:adef1` (the default): `inv(Pdef)*v = inv(P)*(I - J*Q)*v + Q*v`, applied as
  `inv(P)*v + W*(C'*v)` with `W = X - inv(P)*C` folded in at build time, so
  an application is one base solve plus two dense mat-vecs with an `n` by
  `k` matrix. `range(C)` is an exact eigenspace of `J*inv(Pdef)` at
  eigenvalue one. The build costs `k` Jacobian products and `2k` base
  solves.
- `:adef2`: `inv(Pdef)*v = (I - Q*J)*inv(P)*v + Q*v`, applied
  as `u = inv(P)*v; u + X*(C'*(v - J*u))`. This needs `J*u`, which on its
  own costs a Jacobian product per application; inside [`gmres!`](@ref) the
  product is fused with the Arnoldi step ([`preconditionedproduct!`](@ref)),
  where `J*inv(Pdef)*v = J*u + C*(C'*(v - J*u))` is formed from the one
  product the step takes anyway, so an Arnoldi step costs one base solve,
  one Jacobian product and three dense mat-vecs. The build costs `k`
  Jacobian products and `k` base solves (none of the latter when the base
  has not changed since the last build). Outside the fused step (the restart
  correction of `gmres!`, a stagnated step, an external Krylov solver) the
  standalone form pays its product, which is counted in
  [`deflationproducts`](@ref).

With an exact `Q` the two preconditioned operators are block triangular in a
basis adapted to the projector and have the same spectrum; they differ in
sensitivity to an inexact `Q`, which here means a subspace harvested at an
earlier Newton iterate or an earlier sweep point. In the deflation taxonomy
of Tang, Nabben, Vuik and Erlangga (J. Sci. Comput. 39, 2009) these are the
A-DEF1 and A-DEF2 variants. Which is the robust one depends on the coarse
space it is paired with: on identical candidates from a 128-junction
two-tone line, the image pair here takes 1891 Arnoldi steps as A-DEF1 and
2849 as A-DEF2, while the Galerkin pairing took 1679 as A-DEF2 and 3058 as
A-DEF1. So A-DEF1 is the default on this space, 13% behind the best
Galerkin form on that case and without its failure mode.

The pair is rebuilt at every Newton step from the current Jacobian, whether
or not the base was rebuilt ([`pointmoved!`](@ref)); with a frozen base the
rebuild is what turns newly harvested candidates into an active correction,
and the base solves of the candidates it already holds are kept. The
subspace is *accumulated* rather than replaced: each harvest is appended to
the candidates and re-orthonormalized, and the rebuild trims the candidates
to `kmax` by keeping the directions the preconditioned operator shrinks
most, the smallest singular values of `J*inv(P)*U`. That is the metric of
the harvest, and it is the one that matters: measured under `J` alone the
images of the candidates span the many orders of magnitude of the harmonic
scales, and a trim in that metric keeps the wrong directions. Accumulating
measurably beats rediscovering the subspace at every step, because the
deficient directions move slowly while the Newton iterate does not.
`deflationsize` reports the active rank of the pair, which can be smaller
than the number of candidates when a candidate's image falls below the
numerical resolution of the build.

What this can and cannot do, measured on two-tone lines with an explicit
direct current block: a deficiency of rank comparable to `kmax` is removed
(the 128-junction line above); the high-rank deficiency the block diagonal
leaves on longer lines and stronger drives is not, 20 to 40 directions
change nothing there and the escalation does the work, so the wrapper is a
net cost. On a base strong enough to converge in ten Arnoldi steps
(`krylovcouplingmodes = :auto`) the rebuild is the dominant cost and the
wrapper loses. It is off by default (`krylovrecycle = 0`) for those reasons,
and is the natural fit for a sweep at a size where the escalation is what
it replaces.

The intended base is the mode block diagonal
([`ModeCouplingPreconditioner`](@ref) with `couplingmodes = :none`): a batch
of small independent factorizations, one per mode, plus dense level 2 BLAS
for the deflation, so neither part needs a large sparse factorization, and
both port to a GPU.

`escalatepreconditioner!` is delegated to `inner`, so the usual safety net
still applies when recycling alone is not enough; `escalateafter` lets the
wrapper defer a request while it holds an active subspace which might absorb
the deficiency instead (the default is no deferral).

!!! note "`:adef2` and the Jacobian product"
    The fused product relies on `C` being the product of the *current*
    Jacobian with `X`. Inside [`nlsolvekrylov!`](@ref) that is guaranteed by
    [`pointmoved!`](@ref); a caller driving [`gmres!`](@ref) directly must
    call it whenever the operator changes. Through an external Krylov solver
    the standalone form is used, which is exact and costs one extra Jacobian
    product per application.

!!! note "Aliasing"
    `applypreconditioner!(z, pc, r)` reads `r` after writing `z`, so the two
    must not alias, the same contract as `mul!`. The two-argument `ldiv!`
    allocates its result and is safe.
"""
mutable struct RecyclingPreconditioner{TI,TJ,T<:AbstractFloat,TM<:AbstractMatrix{T},TV<:AbstractVector{T}} <: AbstractPreconditioner
    const inner::TI
    const jvp!::TJ
    # the candidates, which persist across solves, and their base solves
    # `Z = inv(P)*U`, column for column, kept while the base is unchanged
    state::RecyclingState{TM}
    Z::TM
    # the active pair, `J*X = C` with orthonormal `C`, and for `:adef1` the
    # folded block `W = X - inv(P)*C`; system sized and allocated like the
    # vectors of the system, so on a device backend the deflation is a pair
    # of device gemvs
    X::TM
    C::TM
    W::TM
    # a k sized work vector for the coefficients, and a system sized one for
    # the `:adef2` correction `v - J*u`
    a::TV
    t::TV
    const kmax::Int
    const kharvest::Int
    const escalateafter::Int
    const form::Symbol
    escalationrequests::Int
    rebuilds::Int
    # Jacobian products this wrapper took itself: `k` per rebuild and one
    # per standalone `:adef2` application
    products::Int
    # whether the pair is that of the current Jacobian and the current
    # candidates: cleared by `pointmoved!` and by a harvest, restored by a
    # rebuild
    fresh::Bool
end

function RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!,
    n::Integer; kmax::Integer = 20, kharvest::Integer = 8,
    escalateafter::Integer = 1, form::Symbol = :adef1,
    T::Type{<:AbstractFloat} = Float64, state = nothing)
    return RecyclingPreconditioner(inner, jvp!, Vector{T}(undef, n);
        kmax = kmax, kharvest = kharvest, escalateafter = escalateafter,
        form = form, state = state)
end

"""
    RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!,
        b::AbstractVector; kmax = 20, kharvest = 8, escalateafter = 1,
        form = :adef1, state = RecyclingState(b))

Build the deflation wrapper on a system whose vectors are like `b`. The
candidates and the active pair are allocated with `similar(b)`, so they live
wherever `b` does and an application is a base solve plus a few device
gemvs; the `k` by `k` factorizations of the build run on the host. `form` is
`:adef1` or `:adef2`; see the type's documentation. `state` is a
[`RecyclingState`](@ref) to start from, the candidates of a previous solve
of a nearby system; it is mutated by the harvests of this solve.
"""
function RecyclingPreconditioner(inner::AbstractPreconditioner, jvp!,
    b::AbstractVector{T}; kmax::Integer = 20, kharvest::Integer = 8,
    escalateafter::Integer = 1, form::Symbol = :adef1,
    state = nothing) where {T<:AbstractFloat}
    kmax >= 1 || throw(ArgumentError(lazy"`kmax` = $(kmax) must be at least 1."))
    1 <= kharvest <= kmax || throw(ArgumentError(
        lazy"`kharvest` = $(kharvest) must satisfy 1 <= kharvest <= kmax = $(kmax)."))
    form in (:adef1, :adef2) || throw(ArgumentError(
        lazy"`form` = $(form) must be `:adef1` or `:adef2`."))
    n = length(b)
    st = isnothing(state) ? RecyclingState(b) : state
    size(st.U, 1) == n || throw(DimensionMismatch(
        lazy"the recycling state is for dimension $(size(st.U, 1)) but the system has $(n)."))
    return RecyclingPreconditioner(inner, jvp!, st, similar(b, n, 0),
        similar(b, n, 0), similar(b, n, 0), similar(b, n, 0), similar(b, 0),
        similar(b), Int(kmax), Int(kharvest), Int(escalateafter), form, 0, 0,
        0, size(st.U, 2) == 0)
end

"""
    deflationform(pc)

The deflation form of a preconditioner, `:adef1` or `:adef2` for a
[`RecyclingPreconditioner`](@ref) and `:none` for anything else.
"""
deflationform(::AbstractPreconditioner) = :none
deflationform(pc::RecyclingPreconditioner) = pc.form

"""
    pointmoved!(pc::AbstractPreconditioner)

Tell the preconditioner that the operator it approximates has changed
without it having been rebuilt, and return `pc`. The default does nothing.
A [`RecyclingPreconditioner`](@ref) marks its pair stale, so that the next
application rebuilds it from the current Jacobian: for `:adef2` the fused
product needs that to be exact, and for either form a stale pair is a
weaker preconditioner. Called by [`nlsolvekrylov!`](@ref) at every Newton
step; wrappers forward it.
"""
pointmoved!(pc::AbstractPreconditioner) = pc
pointmoved!(pc::RecyclingPreconditioner) = (pc.fresh = false; pc)

"""
    stalled!(pc::AbstractPreconditioner)

Tell the preconditioner that the last linear solve reduced its residual
slowly, by less than `krylovrefreshrate` per Arnoldi step, and return `pc`.
The default does nothing. A [`ModeCouplingPreconditioner`](@ref) with
`couplingmodes = :clusters` takes it as the sign that the coupling has
outgrown its clusters and remeasures them at the next update. Called by
[`nlsolvekrylov!`](@ref) after every linear solve; wrappers forward it.
"""
stalled!(pc::AbstractPreconditioner) = pc
stalled!(pc::RecyclingPreconditioner) = (stalled!(pc.inner); pc)

# Escalating the base and recycling are two answers to the same problem,
# and `escalateafter` lets recycling defer the base's escalation in the
# hope that the deflation absorbs the deficiency. The default is 1, no
# deferral: with a long restart cycle a deferred escalation costs up to a
# full linear solve budget of Arnoldi iterations, and it only pays if the
# deflation improves between failures, which in practice it does not; when
# the deflation does absorb the deficiency no escalation is requested at
# all.
function escalatepreconditioner!(pc::RecyclingPreconditioner)
    # Withholding an escalation is only defensible when there is an active
    # deflation which might cover the deficiency instead; candidates whose
    # pair has not been built, or fell below the build's resolution, are
    # not that
    if deflationsize(pc) == 0
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

# The diagnostics of a RecyclingPreconditioner, forwarded by the wrappers
# and zero for any other preconditioner.
"""
    deflationsize(pc)

The active rank of the deflation pair of a
[`RecyclingPreconditioner`](@ref): how many directions the applied
correction spans. Zero for a preconditioner which does not recycle.
"""
deflationsize(::AbstractPreconditioner) = 0
deflationsize(pc::RecyclingPreconditioner) = size(pc.C, 2)

"""
    candidatecount(pc)

The number of candidate directions a [`RecyclingPreconditioner`](@ref)
holds: at least [`deflationsize`](@ref), and more when a candidate's image
fell below the resolution of the last build or a harvest has added
candidates the next build has not yet seen.
"""
candidatecount(::AbstractPreconditioner) = 0
candidatecount(pc::RecyclingPreconditioner) = size(pc.state.U, 2)

"""
    deflationrebuilds(pc)

How many times the deflation pair of a [`RecyclingPreconditioner`](@ref)
has been built.
"""
deflationrebuilds(::AbstractPreconditioner) = 0
deflationrebuilds(pc::RecyclingPreconditioner) = pc.rebuilds

"""
    deflationproducts(pc)

The number of Jacobian products a [`RecyclingPreconditioner`](@ref) has
taken itself: one per candidate at every build, and one per standalone
`:adef2` application (the restart correction of [`gmres!`](@ref), a
stagnated step, an external Krylov solver). Together with `products` in
[`KrylovSolveInfo`](@ref) this is the exact cost of a solve.
"""
deflationproducts(::AbstractPreconditioner) = 0
deflationproducts(pc::RecyclingPreconditioner) = pc.products

function updatepreconditioner!(pc::RecyclingPreconditioner, x::AbstractVector)
    updatepreconditioner!(pc.inner, x)
    # Against an exact inner preconditioner the correction is identically
    # zero: C = J*inv(J)*... = the images are the candidates' Newton
    # corrections and W = X - inv(P)*C = 0 for `:adef1`, v - J*inv(P)*v = 0
    # for `:adef2`. Building it anyway would cost k Jacobian products (and k
    # exact solves) per Newton step. The candidates are kept in case the
    # base is rebuilt restricted again.
    # the base has changed, so its solves of the candidates are stale
    pc.Z = similar(pc.state.U, size(pc.state.U, 1), 0)
    if isexactpreconditioner(pc.inner)
        _clearactive!(pc)
        return pc
    end
    _builddeflation!(pc)
    return pc
end

function _clearactive!(pc::RecyclingPreconditioner)
    n = size(pc.state.U, 1)
    pc.X = similar(pc.state.U, n, 0); pc.C = similar(pc.state.U, n, 0)
    pc.W = similar(pc.state.U, n, 0)
    pc.a = similar(pc.state.U, 0)
    pc.fresh = true
    return pc
end

# A small dense matrix built on the host and moved to wherever `proto`
# lives. `eigen`, `svd` and `cholesky` are scalar indexed dense kernels, so
# the k by k projected quantities are factorized on the host and the
# result sent back, the same split `GMRESWorkspace` makes for the
# Hessenberg and the Givens rotations.
_hostbuilt(proto::AbstractMatrix, A::AbstractMatrix) =
    copyto!(similar(proto, size(A)...), A)
_hostbuilt(proto::SubArray, A::AbstractMatrix) = _hostbuilt(parent(proto), A)

# apply an in-place map to the columns of `A` into `B`, through contiguous
# vectors rather than column views, which a device direct solver cannot
# bind a descriptor to
function _mapcolumns!(f!, B::AbstractMatrix, A::AbstractMatrix, cin, cout)
    for j in axes(A, 2)
        copyto!(cin, view(A, :, j))
        f!(cout, cin)
        copyto!(view(B, :, j), cout)
    end
    return B
end

"""
    _builddeflation!(pc::RecyclingPreconditioner)

Build the active pair for the current candidates at the current point. The
base solves `Z = inv(P)*U` are kept from the last build for the candidates
it already held, since the base is rebuilt through `updatepreconditioner!`,
which empties them; the new candidates cost one base solve each. The
images `Y = J*Z` cost one product per candidate. Their singular value
decomposition is taken through the `k` by `k` Gram matrix on the host,
`Y'*Y = V*S^2*V'`, and does three things at once: it drops the directions
whose image is below the resolution of the Gram form (`S < sqrt(eps)*Smax`,
where a squared condition number stops resolving; such a direction would
otherwise be inverted into an unbounded correction), it trims to `kmax` by
keeping the *smallest* remaining singular values, the directions the
preconditioned operator shrinks most and the ones GMRES cannot resolve on
its own, and it normalizes the rest, `X <- Z*V/S` and `C <- Y*V/S`, so that
`J*X = C` with `C'*C = I` up to the Gram error, which a second Gram pass
on `C` removes. The candidates and their base solves are replaced by their
trimmed span, `U*V` and `Z*V`. For `:adef1` the folded block
`W = X - inv(P)*C` costs `k` further base solves.
"""
function _builddeflation!(pc::RecyclingPreconditioner{TI,TJ,T}) where {TI,TJ,T}
    U = pc.state.U
    kc = size(U, 2)
    n = size(U, 1)
    kc == 0 && return _clearactive!(pc)
    cin = similar(U, n)
    cout = pc.t
    solve! = (o, i) -> applypreconditioner!(o, pc.inner, i)
    kz = min(size(pc.Z, 2), kc)
    Z = similar(U)
    kz > 0 && copyto!(view(Z, :, 1:kz), view(pc.Z, :, 1:kz))
    kz < kc && _mapcolumns!(solve!, view(Z, :, kz+1:kc), view(U, :, kz+1:kc),
        cin, cout)
    Y = similar(U)
    _mapcolumns!(pc.jvp!, Y, Z, cin, cout)
    pc.products += kc
    F = eigen(Symmetric(Array(Y'*Y)))
    s = sqrt.(max.(F.values, zero(T)))
    smax = maximum(s; init = zero(T))
    floor = sqrt(eps(T))*smax
    # ascending, so the smallest resolved singular values come first
    keep = [i for i in eachindex(s) if s[i] > floor]
    length(keep) > pc.kmax && (keep = keep[1:pc.kmax])
    if isempty(keep)
        pc.state.U = similar(U, n, 0)
        pc.Z = similar(U, n, 0)
        return _clearactive!(pc)
    end
    V = F.vectors[:, keep]
    Vd = _hostbuilt(U, V)
    T1 = _hostbuilt(U, V*Diagonal(inv.(s[keep])))
    X = Z*T1
    C = Y*T1
    # the second pass: `C` is orthonormal to the Gram error, `eps` times the
    # squared condition of the images, which the floor above bounds; one
    # Cholesky of `C'*C` finishes it, applied to `X` as well so that
    # `J*X = C` survives
    R = cholesky(Symmetric(Array(C'*C)), check = false)
    if issuccess(R)
        Rinv = _hostbuilt(U, Matrix(inv(R.U)))
        C = C*Rinv
        X = X*Rinv
    end
    pc.state.U = U*Vd
    pc.Z = Z*Vd
    pc.X = X
    pc.C = C
    k = size(C, 2)
    if pc.form === :adef1
        # fold the base solve of C into the correction, k solves now for
        # one fewer gemv per application
        B = similar(C)
        _mapcolumns!(solve!, B, C, cin, cout)
        pc.W = X .- B
    else
        pc.W = similar(U, n, 0)
    end
    pc.a = similar(U, k)
    pc.rebuilds += 1
    pc.fresh = true
    return pc
end

"""
    _refreshproducts!(pc::RecyclingPreconditioner)

Bring the pair up to date with the current Jacobian and candidates when
the point has moved ([`pointmoved!`](@ref)) or a harvest has added
candidates since the last build, whether or not the base was rebuilt. A
no-op when the pair is fresh. Without this a subspace would only ever be
applied after a base refresh, and a base frozen across Newton steps
(`krylovrefreshiterations` large) would leave the deflation inert.
"""
function _refreshproducts!(pc::RecyclingPreconditioner)
    pc.fresh && return pc
    if isexactpreconditioner(pc.inner)
        # the correction is identically zero against an exact base
        _clearactive!(pc)
        return pc
    end
    return _builddeflation!(pc)
end

function applypreconditioner!(z::AbstractVector, pc::RecyclingPreconditioner,
    r::AbstractVector)
    _refreshproducts!(pc)
    if pc.form === :adef1 && size(pc.C, 2) > 0
        mul!(pc.a, transpose(pc.C), r)
    end
    applypreconditioner!(z, pc.inner, r)
    size(pc.C, 2) > 0 || return z
    if pc.form === :adef1
        mul!(z, pc.W, pc.a, true, true)
    else
        # the standalone form of `:adef2` pays the Jacobian product itself
        pc.jvp!(pc.t, z)
        pc.products += 1
        pc.t .= r .- pc.t
        mul!(pc.a, transpose(pc.C), pc.t)
        mul!(z, pc.X, pc.a, true, true)
    end
    return z
end

# apply a preconditioner handed to `gmres!` either as an
# `AbstractPreconditioner` or as a bare in-place closure `M(z, v)`
_applyprecond!(z::AbstractVector, M::AbstractPreconditioner,
    v::AbstractVector) = applypreconditioner!(z, M, v)
_applyprecond!(z::AbstractVector, M, v::AbstractVector) = (M(z, v); z)

"""
    preconditionedproduct!(w, z, Aop, M, v)

One right preconditioned Arnoldi step: overwrite `z` with `inv(M)*v` and `w`
with `A*z`, and return the seconds spent inside the preconditioner. The
generic method applies `M` and then the operator. A
[`RecyclingPreconditioner`](@ref) of form `:adef2` fuses the two, because
its correction needs `J*inv(P)*v`, which is the product the step takes
anyway: with `u = inv(P)*v`, `w0 = J*u` and `c = C'*(v - w0)`,

    z = u + X*c,    w = J*z = w0 + C*c

exactly, since `J*X = C`. The timing excludes the Jacobian product so that
`precondtime` in [`KrylovSolveInfo`](@ref) means the same thing for every
form.
"""
function preconditionedproduct!(w::AbstractVector, z::AbstractVector, Aop,
    M, v::AbstractVector)
    return _unfusedproduct!(w, z, Aop, M, v)
end

function _unfusedproduct!(w, z, Aop, M, v)
    tpc = time()
    _applyprecond!(z, M, v)
    tpc = time() - tpc
    mul!(w, Aop, z)
    return tpc
end

function preconditionedproduct!(w::AbstractVector, z::AbstractVector, Aop,
    pc::RecyclingPreconditioner, v::AbstractVector)
    t0 = time()
    _refreshproducts!(pc)
    tpc = time() - t0
    if pc.form !== :adef2 || size(pc.C, 2) == 0
        return tpc + _unfusedproduct!(w, z, Aop, pc, v)
    end
    t0 = time()
    applypreconditioner!(z, pc.inner, v)
    tpc += time() - t0
    mul!(w, Aop, z)
    t0 = time()
    pc.t .= v .- w
    mul!(pc.a, transpose(pc.C), pc.t)
    mul!(z, pc.X, pc.a, true, true)
    mul!(w, pc.C, pc.a, true, true)
    return tpc + (time() - t0)
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
# ill conditioned block; the caller resolves the rank first (`_orthappend`)
# and this is the cleanup pass.
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
        # the squared condition number either way.
        F = cholesky(Symmetric(Array(A'*A)), check = false)
        issuccess(F) || return -1
        Rinv = _hostbuilt(A, Matrix(inv(F.U)))
        copyto!(A, A*Rinv)
    end
    return k
end

# Append `Xnew` to the orthonormal basis `X`, keeping the result orthonormal
# and never adding a direction outside `span(X, Xnew)`. The projection
# against `X` is block classical Gram-Schmidt with one reorthogonalization,
# level 3 BLAS. The rank of what remains is then revealed on the host from
# the `k` by `k` Gram matrix: an eigenvalue below `eps` times the largest
# incoming column's squared norm (or the block's own largest eigenvalue) is
# a direction the projection removed, a candidate already in the span or a
# duplicate of another, and is dropped rather than orthogonalized into
# noise. A full orthogonal factor of a rank deficient block would instead
# manufacture arbitrary orthogonal-complement columns after the true range
# is exhausted, which is what a Householder fallback used to do here.
function _orthappend(X::AbstractMatrix{T},
    Xnew::AbstractMatrix{T}) where {T<:AbstractFloat}
    size(Xnew, 2) == 0 && return X
    scale = maximum(Array(vec(sum(abs2, Xnew; dims = 1))); init = zero(T))
    if size(X, 2) > 0
        for _ in 1:2
            mul!(Xnew, X, X'*Xnew, -one(T), one(T))
        end
    end
    F = eigen(Symmetric(Array(Xnew'*Xnew)))
    lmax = max(maximum(F.values; init = zero(T)), scale)
    keep = [i for i in eachindex(F.values) if F.values[i] > eps(T)*lmax]
    isempty(keep) && return X
    R = F.vectors[:, keep]*Diagonal(inv.(sqrt.(F.values[keep])))
    B = Xnew*_hostbuilt(Xnew, R)
    # the rank is settled; one CholeskyQR2 pass cleans up the Gram error
    _choleskyqr2!(B)
    return size(X, 2) > 0 ? hcat(X, B) : B
end

"""
    harvestdimension(ws::GMRESWorkspace, out::NamedTuple)

The number of Arnoldi vectors of the *last* restart cycle still present in
the workspace, which is the usable dimension for a harvest; the workspace
holds only that cycle, so this is derived from `out.iterations` and
`out.cycles` rather than from the iteration count alone. Zero when the
cycle is empty or overran.
"""
function harvestdimension(ws::GMRESWorkspace, out::NamedTuple)
    m = size(ws.H, 2)
    j = Int(out.iterations) - (Int(out.cycles) - 1)*m
    return 1 <= j <= m ? j : 0
end

# The directions this Krylov space found the operator shrinks most are the
# right singular vectors of the rectangular Hessenberg with the smallest
# singular values, lifted through the Arnoldi basis. A cycle shorter than
# `kharvest` is not harvested: the singular values of a Krylov space of two
# or three vectors say little about the operator, and taking them anyway
# (so that a short warm-started solve still contributes) was measured to
# cost 21% more Arnoldi steps on a 128-junction two-tone line, since the
# quality trim of the next build cannot tell such a direction from a good
# one. They are directions `u` of the residual on which `J*inv(Pdef)` is
# small. The
# candidate stored is `u` for `:adef2` and `u - C*(C'*u)` for `:adef1`: the
# base part of `inv(Pdef)*u` is `inv(P)` applied to that vector, and the
# rest of it lies in `span(X)`, which the append would remove anyway. No
# base solve and no Jacobian product here; the build pays them.
function harvest!(pc::RecyclingPreconditioner{TI,TJ,T}, ws::GMRESWorkspace,
    out::NamedTuple) where {TI,TJ,T}
    # nothing to deflate against an exact inner; see updatepreconditioner!
    isexactpreconditioner(pc.inner) && return pc
    j = harvestdimension(ws, out)
    j >= pc.kharvest || return pc
    nh = pc.kharvest
    Vj = view(ws.V, :, 1:j)
    F = svd(Matrix(view(ws.Harnoldi, 1:j+1, 1:j)))
    # the Hessenberg is host resident, so the right singular vectors are
    # built there and moved to the basis before lifting them through it
    Cs = _hostbuilt(Vj, F.V[:, size(F.V, 2)-nh+1:size(F.V, 2)])
    Unew = Vj*Cs
    if pc.form === :adef1 && size(pc.C, 2) > 0
        mul!(Unew, pc.C, pc.C'*Unew, -one(T), one(T))
    end
    before = size(pc.state.U, 2)
    pc.state.U = _orthappend(pc.state.U, Unew)
    # new candidates have no base solve and no image yet: the next build
    # makes them
    size(pc.state.U, 2) == before || (pc.fresh = false)
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
        _applyprecond!(z, Mop!, u)
        @. x += z
    end
    return x
end

"""
    gmres!(x, Aop!, b, ws::GMRESWorkspace; Mop! = nothing, rtol = 1e-6,
        atol = 0.0, maxrestarts = 10, initialzero = true, oncycle = nothing)

Solve `A*x = b` with restarted GMRES, where `mul!(w, Aop, v)` computes `w = A*v` and
the optional `Mop!` applies a preconditioner `z = M \\ v`, either as a bare
in-place closure `Mop!(z, v)` or as an [`AbstractPreconditioner`](@ref), which
is applied through [`applypreconditioner!`](@ref) and may fuse its application
with the operator product ([`preconditionedproduct!`](@ref)). The matrix `A`
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
costs one at the outset; `products` in the returned tuple is that total, not
counting products a preconditioner takes inside its own application. `maxrestarts` bounds the number of cycles including
the first, so the Arnoldi work is capped at `maxrestarts*m` steps.
`oncycle(ws, j)`, when given, is called at the end of every cycle with the
workspace still holding that cycle's `j` Arnoldi vectors, for a caller
which harvests from each cycle ([`harvestcycle!`](@ref)); it must only read
the workspace.

Allocation free after the workspace is built, apart from whatever `Aop!` and
`Mop!` themselves allocate.
"""
function gmres!(x::AbstractVector{T}, Aop_, b::AbstractVector{T},
    ws::GMRESWorkspace{T}; Mop! = nothing, rtol = 1e-6, atol = 0.0,
    maxrestarts::Integer = 10, initialzero::Bool = true,
    oncycle = nothing) where {T<:AbstractFloat}

    n = length(b)
    # a bare in-place product is accepted alongside any `mul!`-able operator
    Aop = asoperator(Aop_, n)
    precondtime = 0.0
    products = 0
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
                reason = :converged, residualvector = nothing, products = 0)
        end
        mul!(w, Aop, x)
        products += 1
        resnorm = norm(w)
        if resnorm <= atol
            return (iterations = 0, residual = resnorm, converged = true,
                reason = :converged, residualvector = nothing,
                products = products)
        end
    end

    # initial residual w = b - A*x
    if initialzero
        fill!(x, zero(T))
        copyto!(w, b)
    else
        mul!(w, Aop, x)
        products += 1
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
                mul!(w, Aop, z)
            else
                precondtime += preconditionedproduct!(w, z, Aop, Mop!,
                    view(V, :, j))
            end

            hsub, normw0 = gmres_orthogonalize!(w, V, H, ws.hd, ws.cd, j)
            # the finished column of the Arnoldi Hessenberg, before the
            # rotations below triangularize it in place
            copyto!(view(ws.Harnoldi, 1:j+1, j), view(H, 1:j+1, j))
            resnorm = gmres_applyrotations!(H, cs, sn, s, j)
            totaliterations += 1
            products += 1

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

            V = ensurecolumns!(ws, j + 1)
            @views V[:, j+1] .= w ./ hsub
        end

        gmres_correction!(x, ws, j, Mop!)

        # The Arnoldi factorization of this cycle is about to be overwritten
        # by the next one, so anything that wants to read it has to read it
        # here. A caller which recycles uses this to harvest from *every*
        # cycle rather than only the one left in the workspace at the end:
        # the difficult directions often show up in an early full cycle,
        # while the cycle that finally converges can be a few steps long and
        # carry almost nothing. The callback may only *read* the workspace;
        # the preconditioner is fixed for the duration of a solve, so a
        # harvest appends candidates and nothing is rebuilt until the solve
        # is over.
        isnothing(oncycle) || oncycle(ws, j)

        # recompute the residual explicitly for the next cycle so restarts
        # cannot drift from the recurrence estimate
        mul!(w, Aop, x)
        products += 1
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
    # `w` is the explicit residual `b - A x` of the returned `x`: it is
    # recomputed after every cycle and never estimated, so a caller which
    # needs `A x` has it without another product
    return (iterations = totaliterations, residual = resnorm,
        converged = converged, cycles = cycles, reason = reason,
        precondtime = precondtime, residualvector = w, products = products)
end

"""
    meritslope!(Jv, jvp, p, F, ϕ0, w)

The slope `real(F' J p)` of the merit function `ϕ = F'F/2` along the step
`p`, where `p = -Δ` for a linear solve `J Δ ≈ F` which left the explicit
residual `w = F - J Δ`.

Then `J p = w - F` and the slope is `real(F'w) - 2ϕ0`: one inner product,
with the Jacobian vector product the solve already paid for. Without a
valid residual (`w === nothing`), the product is taken.
"""
function meritslope!(Jv, jvp, p, F, ϕ0, w)
    if isnothing(w)
        mul!(Jv, jvp, p)
        return real(dot(F, Jv))
    end
    return real(dot(F, w)) - 2ϕ0
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
    hblinearsolve!(ls, deltax, jvp!, F, ws, Mop!; rtol, atol, maxrestarts,
        oncycle = nothing)

Solve for the Newton step and return the output named tuple `gmres!`
produces: `converged`, `residual`, `iterations`, `cycles`, `reason`.
`jvp!(y, v)` applies the Jacobian, in place, and `Mop!` is the
preconditioner, an [`AbstractPreconditioner`](@ref) or a closure
`Mop!(z, r)`. `oncycle` is the per-cycle callback of [`gmres!`](@ref),
which a solver that does not restart may ignore.

This is the one part of the Newton-Krylov loop with nothing harmonic
balance specific about it: an operator, a right hand side, a preconditioner
and a tolerance. Putting it behind an interface lets an external Krylov
library be substituted without touching anything else.
"""
function hblinearsolve!(::InternalGMRES, deltax, jvp!, F, ws, Mop!;
        rtol, atol, maxrestarts, oncycle = nothing)
    return gmres!(deltax, jvp!, F, ws; Mop! = Mop!, rtol = rtol,
        atol = atol, maxrestarts = maxrestarts, oncycle = oncycle)
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
    nlsolvekrylov!(fj!, jvp!, F, x, pc::AbstractPreconditioner;
        iterations = 1000, ftol = 1e-8, rtol = 0.0,
        linearsolver = InternalGMRES(), workspace = nothing, label = "",
        c1 = 1e-4, safeguard_low = 0.1, safeguard_high = 0.5, maxbacktracks = 10,
        maxbacktrackfailures = 2, krylovrestart = 400,
        krylovmaxrestarts = 4, krylovrefreshiterations = 1,
        krylovrtolmin = 1e-10, krylovrtolmax = 0.9, krylovrtol0 = 0.3,
        krylovgamma = 0.9, krylovalpha = (1 + sqrt(5))/2,
        krylovstagnation = 0.9, krylovescalate = 1, krylovrefreshrate = 0.5,
        krylovrefresh = :count)

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
- `iterations = 1000`: the maximum number of Newton iterations.
- `ftol = 1e-8`: converged when `norm(F) <= ftol`.
- `rtol = 0.0`: an additional relative test, `norm(F) <= rtol*norm(F0)`
    with `F0` the initial residual, satisfied when either holds. A
    residual whose terms are of size `s` cannot be driven below about
    `eps*s` however exact the step, so an absolute tolerance is a statement
    about the problem's units; a relative one is not.
- `linearsolver = InternalGMRES()`: the linear solver of the Newton step,
    [`InternalGMRES`](@ref) or [`KrylovJL`](@ref).
- `workspace = nothing`: a `Ref` holding the [`KrylovVectors`](@ref) of a
    previous solve of the same system, or holding `nothing`, in which case
    the vectors are allocated and stored into it for the next solve. With
    no `Ref` at all they are allocated and dropped.
- `label = ""`: the label of the returned `IterationInfo`.
- `c1 = 1e-4`: the Armijo sufficient decrease constant, in (0, 1/2).
- `safeguard_low = 0.1`, `safeguard_high = 0.5`: the backtracking step
    clamp, as fractions of the previous trial. The line search here halves
    rather than interpolates (see [`backtracking_linesearch!`](@ref)), so
    only `safeguard_high` acts.
- `maxbacktracks = 10`: the trial point budget per line search.
- `maxbacktrackfailures = 2`: the number of consecutive line search
    failures after which the iteration is declared stalled.
- `krylovrestart = 400`: the GMRES restart length. A long cycle is the
    default because a restricted preconditioner leaves a few directions a
    short Krylov space cannot resolve, and a restart discards the progress
    on them; the basis of `krylovrestart + 1` vectors is cheap next to the
    sparse factorization an escalation would build.
- `krylovmaxrestarts = 4`: the GMRES restart budget per solve.
- `krylovrefreshiterations = 1`: the GMRES iteration count above which the
    preconditioner is considered stale and rebuilt before the next step;
    the default rebuilds eagerly.
- `krylovrefreshrate = 0.5`: the preconditioner is also rebuilt when the
    mean residual ratio per Arnoldi step of the last linear solve,
    `(residual/norm(F))^(1/iterations)`, exceeded this value; a better
    staleness signal than the iteration count alone when the forcing term
    is loose.
- `krylovrefresh = :count`: how a rebuild the two rules above ask for is
    decided. `:count` rebuilds. `:probe` first applies the stale
    preconditioner once to the residual and takes one product, which
    measures the one-step reduction `rho = |J P^-1 F - F|/|F|`; the same
    measurement on the fresh preconditioner (`rho_fresh`, with the Arnoldi
    count `k_fresh` of its solve) calibrates the prediction
    `k = k_fresh log(rho_fresh)/log(rho)` of the stale solve's Arnoldi
    count, and the rebuild is skipped when `k` steps at the measured cost
    of a step are cheaper than the measured rebuild plus a fresh solve.
    Everything is measured, so the rule adapts to the device and the
    factorization; it pays when a rebuild is expensive next to a solve,
    as with a [`BlockFactorization`](@ref) of three tones, where it saved
    a fifth to a third of the time. A rebuild forced by a failed, stalled
    or non-descent solve is never skipped.
- `krylovrtolmin = 1e-10`, `krylovrtolmax = 0.9`, `krylovrtol0 = 0.3`: the
    clamp and the initial value of the forcing sequence.
- `krylovgamma = 0.9`, `krylovalpha = (1 + sqrt(5))/2`: the
    Eisenstat-Walker forcing parameters.
- `krylovstagnation = 0.9`: a GMRES solve which did not reduce the linear
    residual below this fraction of the residual norm is treated as
    stagnated, and the preconditioner solve is taken as the step instead.
- `krylovescalate = 1`: the number of consecutive linear solves which make
    progress but fail to reach their tolerance after which the
    preconditioner is escalated (see [`escalatepreconditioner!`](@ref));
    `typemax(Int)` disables escalation.

These defaults are the ones `hbnlsolve` runs with; `krylovkwargs` there
overrides any of them.

Returns an [`IterationInfo`](@ref) with the same per-iteration diagnostics
as [`nlsolve!`](@ref); the `andersonaccepted` record is always false.
"""
function nlsolvekrylov!(fj!::Function, jvp!, F::AbstractVector{T},
    x::AbstractVector{T}, pc::AbstractPreconditioner; iterations = 1000, ftol = 1e-8,
    linearsolver::AbstractHBLinearSolver = InternalGMRES(),
    label = "",
    c1 = 1e-4, safeguard_low = 0.1, safeguard_high = 0.5,
    maxbacktracks::Integer = 10, maxbacktrackfailures::Integer = 2,
    workspace::Union{Nothing,Base.RefValue} = nothing,
    krylovrestart::Integer = 400, krylovmaxrestarts::Integer = 4,
    krylovrefreshiterations::Integer = 1, krylovrtolmin = 1e-10,
    krylovrtolmax = 0.9, krylovrtol0 = 0.3, krylovgamma = 0.9,
    krylovalpha = (1 + sqrt(5))/2,
    krylovstagnation = 0.9, krylovescalate::Integer = 1,
    krylovrefreshrate = 0.5, krylovrefresh::Symbol = :count,
    rtol = 0.0) where {T<:AbstractFloat}

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
    krylovrefresh in (:count, :probe) || throw(ArgumentError(
        lazy"`krylovrefresh` = $(krylovrefresh) must be `:count` or `:probe`."))

    m = min(krylovrestart, length(x))
    kv = if isnothing(workspace) || isnothing(workspace[])
        KrylovVectors(x, F, m)
    else
        workspace[]
    end
    (length(kv.deltax) == length(x) && size(kv.ws.H, 2) == m &&
        length(kv.Fbest) == length(F)) || throw(ArgumentError(
        "the Krylov workspace handed in is for a different system or restart length; hand in a `Ref` to `nothing` to allocate one."))
    isnothing(workspace) || (workspace[] = kv)
    ws, deltax, xcandidate = kv.ws, kv.deltax, kv.xcandidate
    Jv, Fbest = kv.Jv, kv.Fbest

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
    # why the next rebuild was asked for: `:forced` by a failure of the
    # last solve or the start, `:stale` by the count and rate rules, which
    # is the only kind `krylovrefresh = :probe` may skip
    refreshreason = :forced
    # the measurements of the probe rule: the last rebuild's time, the
    # one-step reduction and Arnoldi count of the solve right after it, and
    # the time per Arnoldi step of the last solve
    tfactor = 0.0
    rhofresh = NaN
    kfresh = 0
    tstep = 0.0
    # consecutive linear solves which failed to reach the forcing tolerance,
    # which trigger an escalation of the preconditioner
    linearfailures = 0
    # The forcing sequence is Eisenstat-Walker choice 2 clamped to
    # [krylovrtolmin, krylovrtolmax] and nothing else. A slope-aware cap
    # which tightened the clamp by a factor of four after every damped step
    # used to sit here. It read a short step as a weak direction; on a long
    # pumped line every step is short because the residual is nonlinear
    # along a full-slope Newton direction, and the near-exact solve the cap
    # then demanded returned a longer step in exactly the direction the
    # line search had to damp: 1033 Arnoldi steps against 177 for the same
    # 24 Newton steps on a 2048 cell line. The inexact Newton theory needs
    # only eta < 1 and a sufficient-decrease line search, which is what
    # remains.

    # the preconditioner object itself is handed to the linear solve, so a
    # form which fuses its application with the operator product can
    # (`preconditionedproduct!`); a plain closure would hide that
    Mop! = pc
    # Where the recycled subspace is read out of the Arnoldi factorization.
    # A preconditioner which harvests per cycle gets the callback and is not
    # harvested again after the solve; one which harvests only the cycle
    # left in the workspace gets the call afterwards. Either way nothing is
    # *rebuilt* inside a solve: the preconditioner has to stay fixed for
    # GMRES, so a harvest only banks candidates and the rebuild waits for
    # the `pointmoved!` of the next Newton step.
    recycles = supportsrecycling(linearsolver)
    percycle = recycles && usescycleharvest(pc)
    oncycle = percycle ? (wsc, j) -> harvestcycle!(pc, wsc, j) : nothing
    harvestafter!(out) =
        (recycles && !percycle && harvest!(pc, ws, out); nothing)
    # residual-only adapter for the linesearch, which never needs the
    # Jacobian and therefore does not accept the combined fj! interface
    residual!(Fv, xv) = fj!(Fv, nothing, xv)

    # rebuild the preconditioner at the current point. a preconditioner is
    # free to move the evaluation point of the matrix-free products while
    # rebuilding, so it is resynchronized afterwards
    # the one-step reduction of the preconditioned residual: one solve of
    # the residual and one product, into the scratch the solve overwrites
    function onestepreduction()
        nF = norm(F)
        nF > 0 || return 0.0
        applypreconditioner!(deltax, pc, F)
        mul!(Jv, jvp, deltax)
        Jv .-= F
        return norm(Jv)/nF
    end
    function refreshpreconditioner!()
        t0 = time()
        updatepreconditioner!(pc, x)
        fj!(nothing, nothing, x)
        tfactor = time() - t0
        # the fresh reduction calibrates the probe of later steps
        krylovrefresh === :probe && (rhofresh = onestepreduction())
        return nothing
    end

    # the residual norm at the initial point; every later entry of normF is
    # pushed immediately after a step is accepted, so convergence is decided
    # on each fresh residual and no preconditioner is ever assembled at a
    # final point
    residual!(F, x)
    push!(normF, norm(F))
    # `ftol` is absolute; `rtol` adds a relative test beside it, and with
    # the default `rtol = 0` the tolerance is exactly `ftol`
    ftol = max(ftol, rtol*normF[1])
    normF[end] <= ftol && (converged = true)

    for n in 1:iterations
        converged && break

        # the matrix-free product reads the evaluation point held by the
        # caller, and the linesearch leaves it at the last trial point
        # rather than the accepted one, so resynchronize it
        fj!(nothing, nothing, x)
        # the Jacobian the preconditioner's deflation was measured against
        # is that of the previous step; a refresh below rebuilds it, and a
        # form which can refresh cheaply without the base does so lazily
        pointmoved!(pc)
        if refresh && refreshreason === :stale && krylovrefresh === :probe &&
                kfresh > 0 && isfinite(rhofresh) && rhofresh > 0
            rho = onestepreduction()
            kpred = rho >= 1 ? Inf : rho <= 0 ? 0.0 :
                kfresh*log(rhofresh)/log(rho)
            refresh = kpred*tstep > tfactor + kfresh*tstep
        end
        justrefreshed = refresh
        if refresh
            refreshpreconditioner!()
            refresh = false
            refreshreason = :forced
        end

        # Eisenstat-Walker choice 2 forcing term from the last accepted
        # step, at its clamp maximum before any step has been taken
        forcing = if length(normF) >= 2 && normF[end-1] > 0
            clamp(krylovgamma*(normF[end]/normF[end-1])^krylovalpha,
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

        tsolve = time()
        out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
            rtol = forcing, atol = gmresatol,
            maxrestarts = krylovmaxrestarts, oncycle = oncycle)
        tsolve = time() - tsolve
        tstep = tsolve/max(out.iterations, 1)
        justrefreshed && (kfresh = max(out.iterations, 1))
        harvestafter!(out)
        # one record per GMRES call; the step outcome is filled in later
        function record!(o, role, refreshedbefore, stag)
            push!(krylovrecord, KrylovSolveInfo(n, role, normF[end], forcing,
                normF[end] > 0 ? o.residual/normF[end] : NaN, o.iterations,
                o.cycles, o.reason, refreshedbefore, false, stag, NaN, NaN,
                0, false, time() - tstart, false, deflationsize(pc),
                deflationrebuilds(pc), get(o, :precondtime, NaN),
                get(o, :products, 0), deflationproducts(pc)))
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
                maxrestarts = krylovmaxrestarts, oncycle = oncycle)
            harvestafter!(out)
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
        if out.iterations > 0 && normF[end] > 0
            rate = (out.residual/normF[end])^(1/out.iterations)
            # a slow solve is what a preconditioner which measures its own
            # structure needs to hear, whether or not it is refreshed
            isfinite(rate) && rate > krylovrefreshrate && stalled!(pc)
        end
        if out.iterations >= krylovrefreshiterations
            refresh || (refreshreason = :stale)
            refresh = true
        elseif out.iterations > 0 && normF[end] > 0
            rate = (out.residual/normF[end])^(1/out.iterations)
            if isfinite(rate) && rate > krylovrefreshrate
                refresh || (refreshreason = :stale)
                refresh = true
            end
        end
        rmul!(deltax, -1)

        # the merit function and its slope along deltax. The assembled
        # Jacobian is stale, so the slope comes from the exact product,
        # which the linear solve already took: its explicit final residual
        # is `F - J Δ`, so the slope is an inner product away. A stagnated
        # solve replaced `deltax` by the preconditioner solve, whose product
        # nothing took, and a linear solver which reports no residual
        # vector leaves it to the product as well.
        ϕ0 = merit(F)
        dϕ0dα = meritslope!(Jv, jvp, deltax, F, ϕ0,
            stagnated ? nothing : get(out, :residualvector, nothing))
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
                maxrestarts = krylovmaxrestarts, oncycle = oncycle)
            # the rescue is often the most informative solve of the step
            harvestafter!(out)
            record!(out, :rescue, true, false)
            rmul!(deltax, -1)
            dϕ0dα = meritslope!(Jv, jvp, deltax, F, ϕ0,
                get(out, :residualvector, nothing))
            if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
                break
            end
        end

        # interpolated backtracking linesearch shared with nlsolve!: on
        # Armijo failure it returns the best decreasing trial (alpha > 0)
        # with F and xcandidate restored there, or alpha == 0 when no trial
        # decreased the merit at all
        # halving rather than interpolating: a trial here costs two
        # transforms, one Arnoldi step of the direction it tests, and the
        # interpolated first backtrack is the floor step whenever the full
        # step overshoots, which on a long pumped line it does at every
        # step (see backtracking_linesearch!)
        alpha1, ϕα, accepted, backtracks = backtracking_linesearch!(
            residual!, F, xcandidate, x, deltax, ϕ0, dϕ0dα;
            c1 = c1, safeguard_low = safeguard_low,
            safeguard_high = safeguard_high,
            maxbacktracks = maxbacktracks, Fbest = Fbest,
            interpolate = false)
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

