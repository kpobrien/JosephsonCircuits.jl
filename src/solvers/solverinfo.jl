# The per-stage diagnostic records every solver returns, and the stall
# diagnostics (`stallmessage`, `projectedstall`) they are read with.

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
    when the residual history projected no convergence within the remaining
    budget without accelerating (`nlsolvekrylov!` first takes one recovery,
    a rebuilt preconditioner and exact Newton steps, and reports it only if
    the stall persists); `:external` for a failed [`ExternalSolver`](@ref).
    [`stallmessage`](@ref) spells each out.
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
- `ftol`: the tolerance in force, the absolute one or the relative one
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
    ftol::T
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
    tracestart!(tr::NewtonTrace, F, ftol, rtol)

Begin (or, on a restart, begin again) the record at a point whose residual
`F` holds: the history is emptied, the tolerance fixed at `ftol` or
`rtol*norm(F)`, whichever is larger, and convergence decided on the
residual before any Jacobian work. Returns whether it has converged.
"""
function tracestart!(tr::NewtonTrace{T}, F, ftol, rtol) where {T}
    empty!(tr.normresidual); empty!(tr.alpha)
    empty!(tr.backtracks); empty!(tr.andersonaccepted)
    tr.backtrackfailures = 0
    tr.converged = false
    tr.reason = :iterations
    push!(tr.normresidual, norm(F))
    tr.ftol = max(T(ftol), T(rtol)*tr.normresidual[1])
    return traceconverged!(tr)
end

# whether the last residual meets the tolerance, recorded as the outcome
function traceconverged!(tr::NewtonTrace)
    if tr.normresidual[end] <= tr.ftol
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
    tracestalled(tr::NewtonTrace, start::Integer, remaining::Integer)

[`projectedstall`](@ref) on the recorded history from `start`, against the
tolerance in force and `remaining` further steps.
"""
tracestalled(tr::NewtonTrace, start::Integer, remaining::Integer) =
    projectedstall(tr.normresidual, start, tr.ftol, remaining)

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
    reason === :progress && return "the residual reduction rate projects no convergence within the remaining budget (a stall; the recovery did not help)"
    reason === :external && return "the external solver reported failure"
    return "reason $(reason)"
end

"""
    projectedstall(normF::AbstractVector, start::Integer, ftol::Real,
        remaining::Integer)

Whether the residual history `normF[start:end]` says the iteration will
not reach `ftol` in `remaining` further steps. The window is split in
half: the geometric reduction rate over the later half projects the steps
still needed, and the verdict is a stall when they exceed `remaining` and
the later rate is no better than the earlier one, so that an iteration
which is accelerating into a Newton basin is never stopped, only one
whose slow progress is steady or worsening. A window shorter than four
steps is never a stall. No constant enters beyond the halving; the
budget the projection is measured against is the caller's own.
"""
function projectedstall(normF::AbstractVector, start::Integer, ftol::Real,
    remaining::Integer)
    m = length(normF)
    npts = m - start + 1
    npts >= 5 || return false
    mid = start + (npts - 1) ÷ 2
    k1 = mid - start
    k2 = m - mid
    r1 = (normF[mid]/normF[start])^(1/k1)
    r2 = (normF[end]/normF[mid])^(1/k2)
    r2 < r1 && return false
    r2 >= 1 && return true
    return log(ftol/normF[end])/log(r2) > remaining
end
