# The per-stage diagnostic records every solver returns, and the stall
# diagnostics (`stallmessage`, `residualstalled`) they are read with.

"""
    AbstractStageInfo

The supertype of every per-stage diagnostic record stored in the `stages`
vector of a `SolverInfo`. Each solver contributes its own concrete record
type instead of adding fields to a shared struct: [`IterationInfo`](@ref)
for the Newton family (with the Krylov linear-solve records where the
solver is `nlsolvekrylov!`), `StagedStageInfo` for the source continuation
driver, and future methods add their own. Every record follows three field
conventions so generic reporting works across methods -- `label::String`,
`converged::Bool` and `iterations::Int` -- and everything else belongs to
the record type itself.
"""
abstract type AbstractStageInfo end

"""
    IterationInfo(label, parameter, regularization, converged, iterations,
        normresidual, alpha, backtracks, andersonaccepted)

Diagnostics recorded for a call of [`nlsolve!`](@ref).

# Fields
- `label`: the solver stage this invocation belongs to.
- `parameter`: the continuation parameter of the stage (the source scale
    or the damping coefficient, depending on the stage), or NaN.
- `regularization`: the diagonal regularization of the Jacobian, if any.
- `converged`: whether the iterations converged.
- `iterations`: the number of Newton iterations performed.
- `normresidual`: the norm of the residual at the start of each iteration.
- `alpha`: the accepted step size for each iteration, or NaN for
    iterations where an Anderson extrapolation was accepted instead of a
    Newton step.
- `backtracks`: the number of linesearch backtracks for each iteration.
- `andersonaccepted`: whether an Anderson extrapolation was accepted for
    each iteration.
- `krylov`: a `KrylovSolveInfo` record for every linear solve performed by
    [`nlsolvekrylov!`](@ref), with one entry per GMRES call rather than per
    Newton step, so that retries and rescues are visible. Empty for the direct
    solvers, which take each step from a factorization.
- `reason`: why the iteration ended. `:converged`; `:iterations` when the
    Newton step budget was spent; `:work` when the Krylov work budget was
    spent (`nlsolvekrylov!` only); `:linesearch` when the line search found
    no sufficient decrease along the Newton direction, once with no decrease
    at all (twice in `nlsolvekrylov!`, which retries the first from a
    rebuilt preconditioner, and which also reports a direction that is not
    a descent direction after its exact rescue here) or twice in a row with
    a decrease short of the Armijo condition, which is a stall; `:progress`
    when the residual stopped coming down and its rate is not improving,
    or in `nlsolvekrylov!` comes down too slowly to reach the tolerance
    within the remaining budget ([`residualstalled`](@ref); that loop
    first takes one recovery, a rebuilt preconditioner and exact Newton
    steps, and reports the stall only if it persists); `:external` for a
    failed [`ExternalSolver`](@ref). [`stallmessage`](@ref) spells each
    out.
"""
struct IterationInfo <: AbstractStageInfo
    label::String
    parameter::Float64
    regularization::Float64
    converged::Bool
    iterations::Int
    normresidual::Vector{Float64}
    alpha::Vector{Float64}
    backtracks::Vector{Int}
    andersonaccepted::Vector{Bool}
    krylov::Vector
    reason::Symbol
end
"""
    NewtonTrace{T}

The record and the acceptance rules shared by the two Newton loops,
[`nlsolve!`](@ref) and [`nlsolvekrylov!`](@ref): the residual norm history,
the step lengths and backtrack counts of every trial, the count of
consecutive line searches which returned the best decreasing trial rather
than an Armijo step, the tolerance, and why the iteration ended. The loops
differ in how they compute a direction and in what they do about a stall;
what a step is, when the iteration has converged and when repeated short
steps are a stall, they share here, through [`tracestart!`](@ref),
[`tracetrial!`](@ref), [`tracestep!`](@ref) and [`tracestalled`](@ref).
The history vectors are the ones the [`IterationInfo`](@ref) of the solve
reports, so a loop may read `normresidual` directly for its forcing terms
and records.

# Fields
- `normresidual`: the residual norm at the starting point and after every
    step taken.
- `alpha`, `backtracks`, `andersonaccepted`: per trial, the step length,
    the trial evaluations after the first, and whether the step lies on
    the accelerated path.
- `atol`: the tolerance in force, the absolute one or the relative one
    times the initial norm, whichever is larger.
- `maxbacktrackfailures`, `backtrackfailures`: the stall threshold and the
    consecutive count against it.
- `converged`, `reason`: the outcome.
"""
mutable struct NewtonTrace{T<:Real}
    const normresidual::Vector{T}
    const alpha::Vector{T}
    const backtracks::Vector{Int}
    const andersonaccepted::Vector{Bool}
    atol::T
    const maxbacktrackfailures::Int
    backtrackfailures::Int
    converged::Bool
    reason::Symbol
end

function NewtonTrace{T}(maxbacktrackfailures::Integer) where {T<:Real}
    return NewtonTrace{T}(T[], T[], Int[], Bool[], zero(T),
        Int(maxbacktrackfailures), 0, false, :iterations)
end

"""
    tracestart!(tr::NewtonTrace, F, atol, rtol)

Begin (or, on a restart, begin again) the record at a point whose residual
`F` holds: the history is emptied, the tolerance fixed at `atol` or
`rtol*norm(F)`, whichever is larger, and convergence decided on the
residual before any Jacobian work. Returns whether it has converged.
"""
function tracestart!(tr::NewtonTrace{T}, F, atol, rtol) where {T}
    empty!(tr.normresidual); empty!(tr.alpha)
    empty!(tr.backtracks); empty!(tr.andersonaccepted)
    tr.backtrackfailures = 0
    tr.converged = false
    tr.reason = :iterations
    push!(tr.normresidual, norm(F))
    tr.atol = max(T(atol), T(rtol)*tr.normresidual[1])
    return traceconverged!(tr)
end

# whether the last residual meets the tolerance, recorded as the outcome
function traceconverged!(tr::NewtonTrace)
    if tr.normresidual[end] <= tr.atol
        tr.converged = true
        tr.reason = :converged
    end
    return tr.converged
end

"""
    tracetrial!(tr::NewtonTrace, alpha, backtracks, anderson::Bool = false)

Record the outcome of one line search: the step length it returned (zero
when no trial decreased the merit, `NaN` for an accelerated candidate
accepted outright), the trial evaluations after its first, and whether
the step lies on the accelerated path. Recorded before it is known
whether the step is taken, since a zero step is not.
"""
function tracetrial!(tr::NewtonTrace, alpha, backtracks, anderson::Bool = false)
    push!(tr.alpha, alpha)
    push!(tr.backtracks, backtracks)
    push!(tr.andersonaccepted, anderson)
    return nothing
end

"""
    tracestep!(tr::NewtonTrace, F, accepted::Bool)

Record the step taken to a point whose residual `F` holds: its norm is
appended, convergence is decided on it, and the consecutive line searches
which returned the best decreasing trial rather than an Armijo accepted
step are counted, an accepted step resetting the count so that an isolated
failure in an otherwise recovering solve is ignored. Returns `:converged`,
`:linesearch` when the count reached `maxbacktrackfailures`, which is a
stall recorded as the reason, or `:continue`.
"""
function tracestep!(tr::NewtonTrace, F, accepted::Bool)
    push!(tr.normresidual, norm(F))
    traceconverged!(tr) && return :converged
    if accepted
        tr.backtrackfailures = 0
    else
        tr.backtrackfailures += 1
        if tr.backtrackfailures >= tr.maxbacktrackfailures
            tr.reason = :linesearch
            return :linesearch
        end
    end
    return :continue
end

"""
    tracestalled(tr::NewtonTrace, start::Integer; remaining = nothing)

[`residualstalled`](@ref) on the recorded history from `start`, against
the tolerance in force and `remaining` further steps when a budget is
given.
"""
tracestalled(tr::NewtonTrace, start::Integer; remaining = nothing) =
    residualstalled(tr.normresidual, start; atol = tr.atol, remaining)

"""
    IterationInfo(tr::NewtonTrace, label, krylov = [])

The record of a solve from its trace, with the Krylov records of
[`nlsolvekrylov!`](@ref) when there are any.
"""
function IterationInfo(tr::NewtonTrace, label, krylov = [])
    return IterationInfo(label, NaN, 0.0, tr.converged, length(tr.alpha),
        tr.normresidual, tr.alpha, tr.backtracks, tr.andersonaccepted,
        krylov, tr.reason)
end

"""
    stallmessage(reason::Symbol)

The sentence behind a `reason` of an [`IterationInfo`](@ref), for the
warning a solve which did not converge issues; a reason it does not know
(`:converged`, `:unspecified`) is reported as `reason <name>`.
"""
function stallmessage(reason::Symbol)
    reason === :iterations && return "the Newton iteration budget was spent"
    reason === :work && return "the Krylov work budget (`iterations` restart lengths of Arnoldi steps) was spent"
    reason === :linesearch && return "the line search found no sufficient decrease along the Newton direction (a stall)"
    reason === :progress && return "the residual stopped coming down, or comes down too slowly for the remaining budget, and its rate is not improving (a stall; the recovery did not help)"
    reason === :external && return "the external solver reported failure"
    return "reason $(reason)"
end

# the residual history the stall rule needs before it judges: long enough
# for the plateaus a solve crosses on its way into a Newton basin to end
# within it
const STALLHISTORY = 20

"""
    residualstalled(normF::AbstractVector, start::Integer,
        history::Integer = STALLHISTORY; atol = nothing, remaining = nothing)

Whether the residual history `normF[start:end]` has stalled. The whole
history from `start` is judged, once it is at least `history` points
long: it is split in half and a geometric rate taken over each. A later
rate better than the earlier one is an iteration accelerating into a
Newton basin, never a stall. Otherwise the residual has stalled when it
is flat or rising, and, given the tolerance `atol` and a budget of
`remaining` further steps, when the steps its later rate projects to the
tolerance, `log(atol/normF)/log(rate)`, exceed the budget. A plateau
after a long descent is judged against the whole descent, so it is a
stall only once it is long enough to bring the later rate to one; a
caller which wants a fresh judgement moves `start` forward, as the Krylov
loop does after its recovery.

The projection is what ends a solve whose line search keeps finding a
decrease: along a descent direction a short enough step always meets the
Armijo condition, and with an inexact direction the residual can creep
for the whole iteration budget without its rate reaching one. Near a
plateau the projection diverges, so only a loop which gives a first
stall a recovery over a fresh history takes it, as the Krylov loop does.
The direct loop judges the rate alone: a direct solve whose residual
keeps coming down, however slowly, runs to its iteration budget.
"""
function residualstalled(normF::AbstractVector, start::Integer,
    history::Integer = STALLHISTORY; atol = nothing, remaining = nothing)
    m = length(normF)
    npts = m - start + 1
    npts >= history || return false
    mid = start + (npts - 1) ÷ 2
    k1 = mid - start
    k2 = m - mid
    r1 = (normF[mid]/normF[start])^(1/k1)
    r2 = (normF[end]/normF[mid])^(1/k2)
    # still accelerating: the second half is coming down faster than the
    # first, so whatever it is doing it is not stalled
    r2 < r1 && return false
    # a residual which is not coming down at all is a stall
    r2 >= 1 && return true
    # and without a budget nothing else is; with one, so is a rate which
    # cannot reach the tolerance within it
    isnothing(remaining) && return false
    return log(atol/normF[end])/log(r2) > remaining
end
