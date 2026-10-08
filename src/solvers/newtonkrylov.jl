# The Newton-Krylov driver: matrix-free Newton steps from GMRES over the
# Jacobian-vector product, with a preconditioner refreshed by policy, a work
# budget, and a reason for every way a solve ends.

"""
    nlsolvekrylov!(fj!, jvp!, F, x, pc::AbstractPreconditioner,
        method::NewtonKrylov = NewtonKrylov(); iterations = 1000,
        atol = 1e-8, rtol = 0.0, workspace = nothing, forcingmin = 1e-10,
        forcingmax = 0.9, forcingstart = 0.3, forcinggamma = 0.9,
        forcingalpha = (1 + sqrt(5))/2, stagnation = 0.9, slowrate = 0.5,
        linearfloor = 0.1)

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

The linear solver of the Newton step, the refresh policy and the
escalation are the `linearsolver`, `refresh` and `escalate` of `method`, a
[`NewtonKrylov`](@ref) whose `preconditioner` and `precision` are not read
here: the preconditioner is `pc`, built by the caller.

`pc` is rebuilt according to `refresh`: before every step for
[`Always`](@ref), by the measured rule of [`Probe`](@ref), and for
[`Never`](@ref) only when forced. A rebuild is forced regardless of the
policy when a solve makes progress but misses its tolerance, when a step is
not a descent direction, when the line search finds no decrease and the
rebuild can change the step (below), and after a successful escalation.
The linear tolerance follows the Eisenstat-Walker choice 2 forcing sequence
`forcinggamma*(|F_k|/|F_{k-1}|)^forcingalpha`, clamped to
`[forcingmin, forcingmax]` and started at `forcingstart`, with an absolute
floor of `linearfloor` times the tolerance in force,
`max(atol, rtol*norm(F0))`, so late solves are not pushed below the
nonlinear tolerance. Because the preconditioner can be stale, the
linesearch slope is always taken from an exact matrix-free product, and a
direction which is not a descent direction is solved for again, at the
same forcing, from a rebuilt preconditioner before the iteration is
declared stalled; for an exact preconditioner that is the exact Newton
step.

The globalization is the plain damped-Newton path of [`nlsolve!`](@ref):
the [`backtracking_linesearch!`](@ref) of [`nlsolve!`](@ref) with the
method's [`Backtracking`](@ref), which on Armijo failure still takes the
best decreasing trial, with consecutive failures counted against its
`maxfailures`, and a step with no decrease retried once from a rebuilt
preconditioner before stopping when the rebuild can change it: when the
preconditioner the step came from was not rebuilt at its point, or has
changed since in a way a rebuild takes in (an escalation granted, a slow
solve reported, on which a [`Clusters`](@ref) request remeasures,
candidates a deflation banked, see [`hasnewcandidates`](@ref)). Against a
preconditioner rebuilt there and unchanged the retry would repeat the step,
and the solve stops at once. There is deliberately no Anderson
acceleration here: the Krylov steps are near-exact Newton steps, and the
solver is kept simple.

# Keywords
- `iterations = 1000`: the maximum number of Newton iterations.
- `atol = 1e-8`: converged when `norm(F) <= atol`.
- `rtol = 0.0`: an additional relative test, `norm(F) <= rtol*norm(F0)`
    with `F0` the initial residual, satisfied when either holds. A
    residual whose terms are of size `s` cannot be driven below about
    `eps*s` however exact the step, so an absolute tolerance is a statement
    about the problem's units; a relative one is not.
- `workspace = nothing`: a `Ref` holding the [`KrylovVectors`](@ref) of a
    previous solve of the same system, or holding `nothing`, in which case
    the vectors are allocated and stored into it for the next solve. With
    no `Ref` at all they are allocated and dropped.
- `forcingmin = 1e-10`, `forcingmax = 0.9`, `forcingstart = 0.3`,
    `forcinggamma = 0.9`, `forcingalpha = (1 + sqrt(5))/2`: the clamp, the
    first term and the parameters of the forcing sequence, on the ranges
    Eisenstat-Walker choice 2 is defined on: `0 < forcingmin <=
    forcingstart <= forcingmax < 1`, `0 < forcinggamma <= 1` and
    `1 < forcingalpha <= 2`.
- `stagnation = 0.9`: a linear solve which does not bring the linear
    residual below this fraction of the residual norm is stagnated, and
    the preconditioner solve is taken as the step.
- `slowrate = 0.5`: a linear solve whose residual was multiplied, on
    average, by more than this factor per Arnoldi step is reported to the
    preconditioner as slow ([`stalled!`](@ref)); none is under
    [`Never`](@ref).
- `linearfloor = 0.1`: the absolute floor of the linear solves, as a
    fraction of the nonlinear tolerance in force.

And, read off `method`: `linearsolver`, the linear solver of the Newton
step, a [`GMRES`](@ref) or a [`KrylovJL`](@ref); `refresh`, when the
preconditioner is rebuilt, [`Always`](@ref) before every step, by the
measured rule of [`Probe`](@ref), or [`Never`](@ref) except when forced
(either way a solve which made progress but missed its tolerance, a
non-descent direction, a line search with no decrease which a rebuild can
change, and a successful escalation force a rebuild); and `escalate`,
whether a preconditioner whose linear solve fails to reach its tolerance
is escalated (see [`escalatepreconditioner!`](@ref)), within the memory
the grown factors are predicted to take; a refused escalation is recorded
(`escalationrequested` in the Krylov record) and the solve carries on.

The forcing sequence, the thresholds of a stagnated and of a slow linear
solve and the floor of the linear solves are keywords of this function
rather than options of the method; [`Never`](@ref) turns off the
slow-solve report and the count rule. The line search is the method's
[`Backtracking`](@ref), interpolating by default.
Two budgets bound the work: `iterations` Newton steps, and `iterations`
restart lengths of Arnoldi steps in total, so that a preconditioner which
runs every linear solve to its limit cannot turn the step budget into
hours. A residual which has stopped coming down, or comes down too slowly
to reach the tolerance within the remaining budget, and whose rate is not
improving ([`residualstalled`](@ref)) gets one recovery, a rebuilt
preconditioner and exact Newton steps from then on, and ends the solve if
it persists.

These are the settings `hbnlsolve` runs with; a caller changes them through
the [`NewtonKrylov`](@ref) method object (`preconditioner`, `linearsolver`,
`refresh`, `escalate`, `linesearch`, `precision`) and through
`hbnlsolve`'s own `iterations`, `atol` and `rtol`.

Returns an [`IterationInfo`](@ref) with the same per-iteration diagnostics
as [`nlsolve!`](@ref) (the `andersonaccepted` record is always false) and a
`reason` of `:converged`, `:iterations`, `:work` (the Arnoldi budget was
spent), `:linesearch` (no decrease, after the retry when one is made, or
a direction which is not a descent direction after the exact rescue),
`:progress`, or `:nonfinite` (the residual norm at the initial point is
not finite, and no step is taken).
"""
function nlsolvekrylov!(fj!::Function, jvp!::Function, F::AbstractVector{T},
    x::AbstractVector{T}, pc::AbstractPreconditioner,
    method::NewtonKrylov = NewtonKrylov(); iterations = 1000, atol = 1e-8,
    rtol = 0.0, workspace::Union{Nothing,Base.RefValue} = nothing,
    kwargs...) where {T<:AbstractFloat}

    length(F) == length(x) || throw(DimensionMismatch(
        lazy"The residual `F` has length $(length(F)) but the point `x` has length $(length(x))."))

    # a bare in-place product is normalized to a `mul!`-able operator once,
    # so the loop below and the pluggable linear solver see one interface.
    # It, the residual and the preconditioner are erased (see `erased`), so
    # that the iteration, the line search and the linear solve are compiled
    # once per vector type rather than once per system
    jvp = asoperator(erased(jvp!), length(x))

    # validate every option before the first residual evaluation; the line
    # search validated its own when it was built
    iterations >= 0 || throw(ArgumentError(
        lazy"`iterations` = $(iterations) must be nonnegative."))
    atol >= 0 || throw(ArgumentError(lazy"`atol` = $(atol) must be nonnegative."))
    0 <= rtol < Inf || throw(ArgumentError(
        lazy"`rtol` = $(rtol) must be finite and nonnegative."))
    # the restart length, which `GMRES` validated when it was built
    m = min(restartlength(method.linearsolver), length(x))
    kv = if isnothing(workspace) || isnothing(workspace[])
        KrylovVectors(x, F, m)
    else
        workspace[]
    end
    (length(kv.deltax) == length(x) && size(kv.ws.H, 2) == m &&
        length(kv.Fbest) == length(F)) || throw(ArgumentError(
        "the Krylov workspace handed in is for a different system or restart length; hand in a `Ref` to `nothing` to allocate one."))
    isnothing(workspace) || (workspace[] = kv)
    # the iteration behind a function barrier, so that it is compiled for
    # the concrete type of a workspace which a reuse holds untyped
    return _nlsolvekrylov!(erased(fj!), jvp, F, x, erased(pc), method, kv;
        iterations = iterations, atol = atol, rtol = rtol, kwargs...)
end

function _nlsolvekrylov!(fj!::Function, jvp, F::AbstractVector{T},
    x::AbstractVector{T}, pc::AbstractPreconditioner, method::NewtonKrylov,
    kv::KrylovVectors; iterations, atol, rtol, forcingmin = 1e-10,
    forcingmax = 0.9, forcingstart = 0.3, forcinggamma = 0.9,
    forcingalpha = (1 + sqrt(5))/2, stagnation = 0.9, slowrate = 0.5,
    linearfloor = 0.1) where {T<:AbstractFloat}

    # the ranges the forcing sequence and the thresholds are defined on
    0 < forcingmin <= forcingstart <= forcingmax < 1 || throw(ArgumentError(
        lazy"the forcing terms must satisfy 0 < `forcingmin` = $(forcingmin) <= `forcingstart` = $(forcingstart) <= `forcingmax` = $(forcingmax) < 1."))
    (0 < forcinggamma <= 1 && 1 < forcingalpha <= 2) || throw(ArgumentError(
        lazy"Eisenstat-Walker choice 2 takes 0 < `forcinggamma` <= 1 and 1 < `forcingalpha` <= 2, not $(forcinggamma) and $(forcingalpha)."))
    (0 < stagnation <= 1 && 0 < slowrate < 1 && 0 <= linearfloor <= 1) ||
        throw(ArgumentError(
            lazy"`stagnation` = $(stagnation) must be in (0, 1], `slowrate` = $(slowrate) in (0, 1) and `linearfloor` = $(linearfloor) in [0, 1]."))

    linearsolver = method.linearsolver
    krylovmaxrestarts = maxrestarts(linearsolver)
    # the refresh policy, read off the method: under `Probe` a measured
    # probe decides whether a rebuild the count asks for pays, and under
    # `Never` the count asks for none and no slow solve is reported
    probing = method.refresh isa Probe
    frozen = method.refresh isa Never
    escalate = method.escalate
    linesearch = method.linesearch

    ws, deltax, xcandidate = kv.ws, kv.deltax, kv.xcandidate
    Jv, Fbest = kv.Jv, kv.Fbest

    ### diagnostic info
    krylovrecord = KrylovSolveInfo[]
    tstart = time()
    # the record and the acceptance rules shared with `nlsolve!`
    tr = NewtonTrace{real(T)}(linesearch.maxfailures)
    normF = tr.normresidual
    refresh = true
    refreshedforstall = false
    # the work budget in Arnoldi steps: `iterations` restart lengths, so
    # that the average Newton step may spend one GMRES cycle and a
    # preconditioner which runs every solve to its limit cannot turn the
    # step budget into hours
    work = 0
    workbudget = iterations*restartlength(linearsolver)
    # the progress rule (`residualstalled`) judges the residual history
    # from `progressstart` against the budget left; its one recovery
    # rebuilds the preconditioner and takes exact Newton steps from then
    # on, ruling out inexact directions before a stall is declared
    progressstart = 1
    exactforcing = false
    # why the next rebuild was asked for: `:forced` by a failure of the
    # last solve or the start, `:stale` by the count and rate rules, which
    # is the only kind `Probe` may skip
    refreshreason = :forced
    # the measurements of the probe rule: the last rebuild's time, the
    # one-step reduction and Arnoldi count of the solve right after it, and
    # the time per Arnoldi step of the last solve
    tfactor = 0.0
    rhofresh = NaN
    kfresh = 0
    tstep = 0.0
    # The forcing sequence is Eisenstat-Walker choice 2 clamped to
    # [forcingmin, forcingmax] and nothing else. A cap which
    # tightened the clamp after a damped step would read a short step as a
    # weak direction; on a long pumped line every step is short because
    # the residual is nonlinear along a full Newton direction, and the near
    # exact solve such a cap demands returns a longer step in exactly the
    # direction the line search has to damp. The inexact Newton theory
    # needs only eta < 1 and a sufficient decrease line search.

    # the preconditioner object itself is handed to the linear solve, which
    # applies it through `applypreconditioner!`
    Mop! = pc
    # A preconditioner which recycles reads the Arnoldi factorization of
    # every restart cycle, when the linear solver exposes it. Nothing is
    # *rebuilt* inside a solve: the preconditioner has to stay fixed for
    # GMRES, so a harvest only banks candidates and the rebuild waits for
    # the next Newton step, its refresh or, without one, its `pointmoved!`.
    oncycle = supportsrecycling(linearsolver) && usescycleharvest(pc) ?
        (wsc, j) -> harvestcycle!(pc, wsc, j) : nothing
    # residual-only adapter for the linesearch, which never needs the
    # Jacobian and therefore does not accept the combined fj! interface
    residual!(Fv, xv) = fj!(Fv, nothing, xv)

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
    # rebuild the preconditioner at the current point. a preconditioner is
    # free to move the evaluation point of the matrix-free products while
    # rebuilding, so it is resynchronized afterwards. Returns the time the
    # rebuild took and, under the probe rule, the one-step reduction of the
    # fresh preconditioner, which calibrate the probe of later steps; the
    # caller assigns them, so no variable is shared with the closure. An
    # exact preconditioner's reduction is roundoff, against which any stale
    # one predicts a rebuild, so it is not measured, and without the
    # calibration no probe is taken: the rebuild is made as under `Always`
    function refreshpreconditioner!()
        t0 = time()
        updatepreconditioner!(pc, x)
        fj!(nothing, nothing, x)
        t = time() - t0
        probe = probing && !isexactpreconditioner(pc)
        return t, probe ? onestepreduction() : NaN
    end

    # the residual norm at the initial point; every later entry of normF is
    # pushed immediately after a step is accepted, so convergence is decided
    # on each fresh residual and no preconditioner is ever assembled at a
    # final point. `atol` is absolute; `rtol` adds a relative test beside
    # it, and with the default `rtol = 0` the tolerance is exactly `atol`.
    # A start which has converged, or whose residual norm is not finite,
    # ends the solve before any preconditioner work
    residual!(F, x)
    tracestart!(tr, F, atol, rtol) && return IterationInfo(tr, krylovrecord)
    # absolute floor for the linear solves, from the tolerance in force, the
    # relative one included: once the linear residual is below the
    # nonlinear tolerance, further accuracy cannot help the Newton
    # iteration, and demanding it makes late GMRES solves "fail"
    gmresatol = linearfloor*tr.atol

    for n in 1:iterations
        # the matrix-free product reads the evaluation point held by the
        # caller, and the linesearch leaves it at the last trial point
        # rather than the accepted one, so resynchronize it
        fj!(nothing, nothing, x)
        # the probe measures the preconditioner the last step left, with the
        # deflation built against that step's Jacobian
        if refresh && refreshreason === :stale && probing &&
                kfresh > 0 && isfinite(rhofresh) && rhofresh > 0
            rho = onestepreduction()
            kpred = rho >= 1 ? Inf : rho <= 0 ? 0.0 :
                kfresh*log(rhofresh)/log(rho)
            refresh = kpred*tstep > tfactor + kfresh*tstep
        end
        justrefreshed = refresh
        if refresh
            tfactor, rhofresh = refreshpreconditioner!()
            refresh = false
            refreshreason = :forced
        else
            # the deflation was built against an earlier step's Jacobian: a
            # form which can rebuild it without the base does so at its next
            # application. A refresh rebuilds the deflation itself, so it is
            # told only when there is none, and a step builds it once
            pointmoved!(pc)
        end
        # whether the preconditioner this step's direction comes from was
        # rebuilt at this point, and whether it has changed since in a way
        # a rebuild here takes in; they decide the retry of a step whose
        # line search finds no decrease
        rebuilthere = justrefreshed
        changed = false

        # Eisenstat-Walker choice 2 forcing term from the last accepted
        # step, and its first term before any step has been taken
        forcing = if length(normF) >= 2 && normF[end-1] > 0
            clamp(forcinggamma*(normF[end]/normF[end-1])^forcingalpha,
                forcingmin, forcingmax)
        else
            # the *initial* forcing term, which is a separate quantity from
            # the upper clamp: seeding the safeguard at forcingmax would let
            # it walk down from there over several outer steps (0.9, 0.76,
            # 0.58, 0.37, 0.18 for the default gamma and alpha), so that
            # several successive linear solves terminate after very little
            # residual reduction
            forcingstart
        end

        exactforcing && (forcing = forcingmin)
        tsolve = time()
        out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
            rtol = forcing, atol = gmresatol,
            maxrestarts = krylovmaxrestarts, oncycle = oncycle)
        tsolve = time() - tsolve
        work += out.iterations
        tstep = tsolve/max(out.iterations, 1)
        justrefreshed && (kfresh = max(out.iterations, 1))
        # a residual which is not finite counts as stagnated
        stagnated = !out.converged &&
            !(out.residual <= stagnation*normF[end])
        push!(krylovrecord, krylovsolverecord(out, n, :step, normF[end],
            forcing, justrefreshed, stagnated, tstart, pc))
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
            tfactor, rhofresh = refreshpreconditioner!()
            rebuilthere = true
            out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
                rtol = forcing, atol = gmresatol,
                maxrestarts = krylovmaxrestarts, oncycle = oncycle)
            # the fresh preconditioner's reduction and Arnoldi count, which
            # calibrate the probe, are this solve's
            kfresh = max(out.iterations, 1)
            work += out.iterations
            stagnated = !out.converged &&
                !(out.residual <= stagnation*normF[end])
            push!(krylovrecord, krylovsolverecord(out, n, :retry, normF[end],
                forcing, true, stagnated, tstart, pc))
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
        # Escalation is triggered by every solve which does not converge,
        # not by stagnation alone. These are different symptoms: stagnation
        # is a Krylov space which contributed nothing, while the common
        # failure of a merely inadequate preconditioner is a solve which
        # makes real progress and still cannot reach the forcing tolerance
        # within its budget, so keying escalation to stagnation alone would
        # leave it never firing on the problems it exists to rescue.
        if escalate && !out.converged
            if escalatepreconditioner!(pc)
                # the escalation dropped the factorization: the rebuild
                # is mandatory, not a probe's economic decision
                refresh = true
                refreshreason = :forced
                changed = true
                krylovrecord[end] = with(krylovrecord[end]; escalated = true)
            else
                krylovrecord[end] = with(krylovrecord[end]; escalationrequested = true)
            end
        end
        # Staleness is judged by how fast the linear residual came down per
        # Arnoldi step, not by the raw step count alone. The count on its own
        # is confounded by the requested accuracy: with a loose forcing term a
        # badly stale preconditioner can satisfy the tolerance in a single
        # iteration after removing only a small fraction of the linear
        # residual, and a rule keyed to the count reads that as health.
        # Under `Never` no solve is reported, whatever its rate, which
        # exceeds one for a linear solver whose residual can grow.
        if !frozen && out.iterations > 0 && normF[end] > 0
            rate = (out.residual/normF[end])^(1/out.iterations)
            # a slow solve is what a preconditioner which measures its own
            # structure needs to hear, whether or not it is refreshed, and
            # one may remeasure at its next rebuild
            if isfinite(rate) && rate > slowrate
                stalled!(pc)
                changed = true
            end
        end
        # the count rule: under `Always` and `Probe` one Arnoldi step is
        # enough to ask for a rebuild (which the probe may then decline),
        # under `Never` no count is
        if !frozen && out.iterations >= 1
            refresh || (refreshreason = :stale)
            refresh = true
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
            # iteration has stalled.
            tfactor, rhofresh = refreshpreconditioner!()
            rebuilthere = true
            changed = false
            out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
                rtol = forcing, atol = gmresatol,
                maxrestarts = krylovmaxrestarts, oncycle = oncycle)
            kfresh = max(out.iterations, 1)
            # the rescue is often the most informative solve of the step,
            # and its Arnoldi steps count against the work budget like any
            # other solve's
            push!(krylovrecord, krylovsolverecord(out, n, :rescue,
                normF[end], forcing, true, false, tstart, pc))
            work += out.iterations
            rmul!(deltax, -1)
            dϕ0dα = meritslope!(Jv, jvp, deltax, F, ϕ0,
                get(out, :residualvector, nothing))
            if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
                tr.reason = :linesearch
                break
            end
        end

        # the backtracking linesearch shared with nlsolve!, fitted or halved
        # as the method's `Backtracking` says: on Armijo failure it returns
        # the best decreasing trial (alpha > 0) with F and xcandidate
        # restored there, or alpha == 0 when no trial decreased the merit
        # at all
        alpha1, ϕα, accepted, backtracks = backtracking_linesearch!(
            residual!, F, xcandidate, x, deltax, ϕ0, dϕ0dα;
            ls = linesearch, Fbest = Fbest)
        tracetrial!(tr, alpha1, backtracks)
        if !isempty(krylovrecord) && krylovrecord[end].iteration == n
            krylovrecord[end] = with(krylovrecord[end];
                slope = normF[end] > 0 ? dϕ0dα/normF[end]^2 : NaN,
                alpha = alpha1, backtracks = backtracks, armijo = accepted)
        end

        if iszero(alpha1)
            # no decrease anywhere along the direction; F again holds the
            # residual at the unchanged x (the linesearch restore contract).
            # A rebuild here gives another direction only when the
            # preconditioner this one came from was not rebuilt at this
            # point, or has changed since in a way a rebuild takes in: an
            # escalation granted, a slow solve reported, on which a cluster
            # request remeasures, or candidates a deflation banked. Then the
            # step is retried once from a rebuilt preconditioner; otherwise
            # a retry would repeat it exactly, and the solve ends, as it
            # does when the retry finds no decrease either
            if refreshedforstall ||
                    (rebuilthere && !changed && !hasnewcandidates(pc))
                tr.reason = :linesearch
                break
            end
            refreshedforstall = true
            refresh = true
            refreshreason = :forced
            continue
        end
        refreshedforstall = false

        # accept the trial point: F already holds the residual there (the
        # linesearch postcondition), and convergence is decided on it now,
        # before any preconditioner work; repeated line searches which took
        # the best decreasing trial instead of an Armijo accepted step are a
        # stall
        copyto!(x, xcandidate)
        tracestep!(tr, F, accepted) === :continue || break
        if work >= workbudget
            tr.reason = :work
            break
        end
        # Armijo accepts a step which reduces the merit by a hair, so a
        # hopeless iteration can satisfy it for its whole budget; the
        # progress rule stops it once the residual has stopped coming
        # down, or comes down too slowly for the budget left, after one
        # recovery
        if tracestalled(tr, progressstart; remaining = iterations - n)
            if exactforcing
                tr.reason = :progress
                break
            end
            exactforcing = true
            progressstart = length(normF)
            refresh = true
            refreshreason = :forced
        end
    end

    return IterationInfo(tr, krylovrecord)
end

# the record of a linear solve `o` of Newton step `n` at the residual norm
# `normF`; the step's outcome is filled in once the line search is done
function krylovsolverecord(o, n, role, normF, forcing, refreshed, stagnated,
        tstart, pc)
    return KrylovSolveInfo(; iteration = n, role = role, normF = normF,
        forcing = forcing,
        residualratio = normF > 0 ? o.residual/normF : NaN,
        iterations = o.iterations, cycles = o.cycles, reason = o.reason,
        refreshed = refreshed, escalated = false, stagnated = stagnated,
        slope = NaN, alpha = NaN, backtracks = 0, armijo = false,
        time = time() - tstart, escalationrequested = false,
        deflationsize = deflationsize(pc),
        deflationrebuilds = deflationrebuilds(pc),
        precondtime = get(o, :precondtime, NaN),
        products = get(o, :products, 0),
        deflationproducts = deflationproducts(pc))
end
