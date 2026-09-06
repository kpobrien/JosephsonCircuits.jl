# The preconditioner interface of the Newton-Krylov solve: the operator wrapper
# handed to a linear solver, the abstract preconditioner and every hook a wrapper
# may implement, and the record of one linear solve.

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
its mode coupling restricted, and [`FloquetPreconditioner`](@ref) for that
wrapped in a deflation subspace.

The rest of the interface is hooks the driver calls with an inert default,
which a preconditioner implements when it has something to say:

| hook | default | who implements it |
|---|---|---|
| `pointmoved!(pc)` | nothing | a deflation, whose pair goes stale |
| `stalled!(pc)` | nothing | `Clusters`, which remeasures |
| `escalatepreconditioner!(pc)` | `false` | mode coupling, a deflation deferring it |
| `harvest!(pc, ws, out)` | nothing | a deflation |
| `usescycleharvest(pc)`, `harvestcycle!(pc, ws, j)` | `false`, nothing | a deflation reading every cycle |
| `seeddeflation!(pc, X; ...)` | nothing | a deflation taking physical candidates |
| `isexactpreconditioner(pc)` | `false` | an exact factorization |
| `deflationsize`, `candidatecount`, `deflationrebuilds`, `deflationproducts` | `0`, `0`, `0`, `0` | a deflation |

A preconditioner which wraps another subtypes
[`AbstractWrappedPreconditioner`](@ref) and defines
[`innerpreconditioner`](@ref); every hook it does not implement itself is then
forwarded to the inner one, so that wrapping never silently turns a hook
off.
"""
abstract type AbstractPreconditioner end

"""
    AbstractWrappedPreconditioner

A preconditioner which wraps another and presents it in different
coordinates, at a different size or under a correction. Subtypes define
[`innerpreconditioner`](@ref) and their own `updatepreconditioner!` and
`applypreconditioner!`; every other hook of the interface forwards to the
inner preconditioner unless the wrapper defines it, so escalation,
deflation and its diagnostics reach the inner one through any number of
wrappers. A wrapper which changes coordinates must define the hooks whose
arguments carry vectors of its own coordinates (`harvest!`,
`harvestcycle!`, `seeddeflation!`) if the inner preconditioner reads them
in its own.
"""
abstract type AbstractWrappedPreconditioner <: AbstractPreconditioner end

"""
    innerpreconditioner(pc::AbstractWrappedPreconditioner)

The preconditioner `pc` wraps.
"""
function innerpreconditioner end

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
- `escalationrequested`: whether an escalation was requested after this
    solve; with `escalated` false, that is an escalation the preconditioner
    refused because the grown factors would not fit its memory budget.
- `deflationsize`, `deflationrebuilds`, `precondtime`: the active rank of
    the recycled deflation, how many times it has been built, and the wall
    time spent applying the preconditioner in this solve.
- `products`, `deflationproducts`: the exact operator products this linear
    solve took, and the running count of those the deflation wrapper took
    for its builds. The cost of a solve is in these, not in `iterations`
    alone: every restart cycle recomputes the residual.
"""
Base.@kwdef struct KrylovSolveInfo
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
    # running count `deflationproducts`, those the deflation wrapper took
    # for its builds; Arnoldi steps alone understate the cost of a restart
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


# the step outcome is only known after the linesearch, so the record is
# completed then
"""
    with(k::KrylovSolveInfo; fields...)
    with(k::IterationInfo; fields...)

A copy of the record with the named fields replaced: the step outcome
filled in after the line search, an escalation marked after the solve, the
drive fraction of a stage.
"""
function with(k::Union{KrylovSolveInfo,IterationInfo}; kwargs...)
    T = typeof(k)
    for name in keys(kwargs)
        hasfield(T, name) || throw(ArgumentError(
            lazy"$(T) has no field $(name)."))
    end
    return T((haskey(kwargs, f) ? kwargs[f] : getfield(k, f)
        for f in fieldnames(T))...)
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
[`ModeCouplingPreconditioner`](@ref) and [`FloquetPreconditioner`](@ref)
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
Called by [`nlsolvekrylov!`](@ref) after repeated linear solves which fail to
reach the forcing tolerance, the symptom of a preconditioner too crude for
the problem (stagnation alone is deliberately not the trigger). The default method
returns `false`, which is correct for any preconditioner that is already exact
or has no cheaper/costlier settings. A preconditioner which can grow must
also refuse when the grown factors would not fit in memory: the driver
records the refusal and carries on with what it has rather than let a
rescue exhaust the machine.
"""
escalatepreconditioner!(::AbstractPreconditioner) = false

"""
    pointmoved!(pc::AbstractPreconditioner)

Tell the preconditioner that the operator it approximates has changed
without it having been rebuilt, and return `pc`. The default does nothing.
A [`FloquetPreconditioner`](@ref) marks its image pair stale, so that the
next application rebuilds it from the current Jacobian. Called by
[`nlsolvekrylov!`](@ref) at every Newton step; wrappers forward it.
"""
pointmoved!(pc::AbstractPreconditioner) = pc

"""
    stalled!(pc::AbstractPreconditioner)

Tell the preconditioner that the last linear solve reduced its residual
slowly, by a factor worse than 0.5 per Arnoldi step (the report is off
under [`Never`](@ref)), and return `pc`.
The default does nothing. A [`ModeCouplingPreconditioner`](@ref) with
[`Clusters`](@ref) takes it as the sign that the coupling has
outgrown its clusters and remeasures them at the next update. Called by
[`nlsolvekrylov!`](@ref) after every linear solve; wrappers forward it.
"""
stalled!(pc::AbstractPreconditioner) = pc

"""
    deflationsize(pc)

The number of directions the active deflation pair of a
[`FloquetPreconditioner`](@ref) spans. Zero for a preconditioner which does
not deflate.
"""
deflationsize(::AbstractPreconditioner) = 0

"""
    candidatecount(pc)

The number of candidate directions a [`FloquetPreconditioner`](@ref)
holds in its bank: at least [`deflationsize`](@ref), and more when a
candidate was left out of the last build or a harvest has added candidates
the next build has not yet seen.
"""
candidatecount(::AbstractPreconditioner) = 0

"""
    deflationrebuilds(pc)

How many times the image pair of a [`FloquetPreconditioner`](@ref) has been
built.
"""
deflationrebuilds(::AbstractPreconditioner) = 0

"""
    deflationproducts(pc)

The number of Jacobian products a [`FloquetPreconditioner`](@ref) has taken
itself, one per candidate at every build. Together with `products` in
[`KrylovSolveInfo`](@ref) this is the exact cost of a solve.
"""
deflationproducts(::AbstractPreconditioner) = 0

# the hooks of a wrapper forward to what it wraps; the two which take a
# `GMRESWorkspace` are in gmres.jl, after the type is defined
pointmoved!(pc::AbstractWrappedPreconditioner) =
    (pointmoved!(innerpreconditioner(pc)); pc)
stalled!(pc::AbstractWrappedPreconditioner) =
    (stalled!(innerpreconditioner(pc)); pc)
escalatepreconditioner!(pc::AbstractWrappedPreconditioner) =
    escalatepreconditioner!(innerpreconditioner(pc))
usescycleharvest(pc::AbstractWrappedPreconditioner) =
    usescycleharvest(innerpreconditioner(pc))
seeddeflation!(pc::AbstractWrappedPreconditioner, X::AbstractMatrix; kwargs...) =
    (seeddeflation!(innerpreconditioner(pc), X; kwargs...); pc)
isexactpreconditioner(pc::AbstractWrappedPreconditioner) =
    isexactpreconditioner(innerpreconditioner(pc))
deflationsize(pc::AbstractWrappedPreconditioner) =
    deflationsize(innerpreconditioner(pc))
candidatecount(pc::AbstractWrappedPreconditioner) =
    candidatecount(innerpreconditioner(pc))
deflationrebuilds(pc::AbstractWrappedPreconditioner) =
    deflationrebuilds(innerpreconditioner(pc))
deflationproducts(pc::AbstractWrappedPreconditioner) =
    deflationproducts(innerpreconditioner(pc))
