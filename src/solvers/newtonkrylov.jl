# The Newton-Krylov driver: matrix-free Newton steps from GMRES over the
# Jacobian-vector product, with a preconditioner refreshed by policy, a work
# budget, and a reason for every way a solve ends.

"""
    nlsolvekrylov!(fj!, jvp!, F, x, pc::AbstractPreconditioner,
        method::NewtonKrylov = NewtonKrylov(); iterations = 1000,
        ftol = 1e-8, rtol = 0.0, workspace = nothing, label = "")

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
not a descent direction, when the line search finds no decrease, and after
a successful escalation. The linear tolerance follows the Eisenstat-Walker choice 2 forcing sequence
`krylovgamma*(|F_k|/|F_{k-1}|)^krylovalpha` clamped to
`[krylovrtolmin, krylovrtolmax]`, with an absolute floor of `ftol/10` so late
solves are not pushed below the nonlinear tolerance. Because the assembled
Jacobian can be stale, the linesearch slope is always taken from an exact
matrix-free product, and a non-descent direction falls back to the exact
Newton step through a fresh factorization before the iteration is declared
stalled.

The globalization is the plain damped-Newton path of [`nlsolve!`](@ref):
the [`backtracking_linesearch!`](@ref) of [`nlsolve!`](@ref) run in halving
mode (`interpolate = false`), which on Armijo failure
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
- `workspace = nothing`: a `Ref` holding the [`KrylovVectors`](@ref) of a
    previous solve of the same system, or holding `nothing`, in which case
    the vectors are allocated and stored into it for the next solve. With
    no `Ref` at all they are allocated and dropped.
- `label = ""`: the label of the returned `IterationInfo`.

And, read off `method`: `linearsolver`, the linear solver of the Newton
step, a [`GMRES`](@ref) or a [`KrylovJL`](@ref); `refresh`, when the
preconditioner is rebuilt, [`Always`](@ref) before every step, by the
measured rule of [`Probe`](@ref), or [`Never`](@ref) except when forced
(either way a solve which made progress but missed its tolerance, a
non-descent direction, a line search with no decrease and a successful
escalation force a rebuild); and `escalate`, whether a preconditioner
which makes progress but fails to reach its tolerance is escalated (see
[`escalatepreconditioner!`](@ref)) rather than tried again, within the
memory the grown factors are predicted to take; a refused escalation is
recorded (`escalationrequested` in the Krylov record) and the solve
carries on.

The forcing sequence is Eisenstat-Walker choice 2 with `gamma = 0.9` and
`alpha = (1 + sqrt(5))/2`, clamped to `[1e-10, 0.9]` and started at 0.3;
the line search is Armijo backtracking with constant 1e-4, halving with
safeguards 0.1 and 0.5, at most ten trials and two consecutive failures;
a solve which does not bring the linear residual below 0.9 of the residual
norm is treated as stagnated and the preconditioner solve taken as the
step; and a solve whose residual came down by less than 0.5 per Arnoldi
step is reported to the preconditioner as slow ([`stalled!`](@ref); off
under [`Never`](@ref), which also disables the count rule). These are
fixed: none has been changed in any measured case, and each was set by the
inexact Newton theory or by a measurement recorded beside it. Two budgets
bound the work: `iterations` Newton steps, and `iterations` restart
lengths of Arnoldi steps in total, so that a preconditioner which runs
every linear solve to its limit cannot turn the step budget into hours.
A residual history which projects no convergence within the remaining
budget without accelerating ([`projectedstall`](@ref)) gets one recovery,
a rebuilt preconditioner and exact Newton steps from then on, and ends the
solve if it persists.

These are the settings `hbnlsolve` runs with; a caller changes them through
the [`NewtonKrylov`](@ref) method object (`preconditioner`, `linearsolver`,
`refresh`, `escalate`, `precision`) and through `hbnlsolve`'s own
`iterations`, `ftol` and `rtol`.

Returns an [`IterationInfo`](@ref) with the same per-iteration diagnostics
as [`nlsolve!`](@ref) (the `andersonaccepted` record is always false) and a
`reason` of `:converged`, `:iterations`, `:work` (the Arnoldi budget was
spent), `:linesearch` (no decrease twice, or a direction which is not a
descent direction after the exact rescue), or `:progress`.
"""
function nlsolvekrylov!(fj!::Function, jvp!, F::AbstractVector{T},
    x::AbstractVector{T}, pc::AbstractPreconditioner,
    method::NewtonKrylov = NewtonKrylov(); iterations = 1000, ftol = 1e-8,
    rtol = 0.0, workspace::Union{Nothing,Base.RefValue} = nothing,
    label = "") where {T<:AbstractFloat}

    linearsolver = method.linearsolver
    refresh = method.refresh
    escalate = method.escalate

    # The fixed constants of the iteration, under the names the loop below
    # uses. The forcing sequence is Eisenstat-Walker choice 2 and its
    # parameters are only defined on these ranges; the line search halves
    # rather than interpolates, so only the upper safeguard acts; the
    # stagnation and slow solve thresholds are those the docstring
    # describes.
    krylovrestart = restartlength(linearsolver)
    krylovmaxrestarts = maxrestarts(linearsolver)
    krylovrefreshiterations = refresh isa Never ? typemax(Int) : 1
    krylovrefreshrate = refresh isa Never ? 1.0 : 0.5
    krylovrefresh = refresh isa Probe ? :probe : :count
    krylovrtolmin = 1e-10
    krylovrtolmax = 0.9
    krylovrtol0 = 0.3
    krylovgamma = 0.9
    krylovalpha = (1 + sqrt(5))/2
    krylovstagnation = 0.9
    krylovescalate = escalate ? 1 : typemax(Int)
    c1 = 1e-4
    safeguard_low = 0.1
    safeguard_high = 0.5
    maxbacktracks = 10
    maxbacktrackfailures = 2

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
    tstart = time()
    # the record and the acceptance rules shared with `nlsolve!`
    tr = NewtonTrace{real(T)}(maxbacktrackfailures)
    normF = tr.normresidual
    refresh = true
    refreshedforstall = false
    # the work budget in Arnoldi steps: `iterations` restart lengths, so
    # that the average Newton step may spend one GMRES cycle and a
    # preconditioner which runs every solve to its limit cannot turn the
    # step budget into hours
    work = 0
    workbudget = iterations*krylovrestart
    # the progress projection (`projectedstall`) measures the residual
    # history from `progressstart`; its one recovery rebuilds the
    # preconditioner and takes exact Newton steps from then on, ruling out
    # inexact directions before a stall is declared
    progressstart = 1
    exactforcing = false
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
    # [krylovrtolmin, krylovrtolmax] and nothing else. A cap which
    # tightened the clamp after a damped step would read a short step as a
    # weak direction; on a long pumped line every step is short because
    # the residual is nonlinear along a full Newton direction, and the near
    # exact solve such a cap demands returns a longer step in exactly the
    # direction the line search has to damp. The inexact Newton theory
    # needs only eta < 1 and a sufficient decrease line search.

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
    # final point. `ftol` is absolute; `rtol` adds a relative test beside
    # it, and with the default `rtol = 0` the tolerance is exactly `ftol`
    residual!(F, x)
    tracestart!(tr, F, ftol, rtol)

    for n in 1:iterations
        tr.converged && break

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

        exactforcing && (forcing = krylovrtolmin)
        tsolve = time()
        out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
            rtol = forcing, atol = gmresatol,
            maxrestarts = krylovmaxrestarts, oncycle = oncycle)
        tsolve = time() - tsolve
        work += out.iterations
        tstep = tsolve/max(out.iterations, 1)
        justrefreshed && (kfresh = max(out.iterations, 1))
        harvestafter!(out)
        # one record per GMRES call; the step outcome is filled in later
        function record!(o, role, refreshedbefore, stag)
            push!(krylovrecord, KrylovSolveInfo(; iteration = n, role = role,
                normF = normF[end], forcing = forcing,
                residualratio = normF[end] > 0 ? o.residual/normF[end] : NaN,
                iterations = o.iterations, cycles = o.cycles,
                reason = o.reason, refreshed = refreshedbefore,
                escalated = false, stagnated = stag, slope = NaN, alpha = NaN,
                backtracks = 0, armijo = false, time = time() - tstart,
                escalationrequested = false, deflationsize = deflationsize(pc),
                deflationrebuilds = deflationrebuilds(pc),
                precondtime = get(o, :precondtime, NaN),
                products = get(o, :products, 0),
                deflationproducts = deflationproducts(pc)))
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
            work += out.iterations
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
                    # the escalation dropped the factorization: the rebuild
                    # is mandatory, not a probe's economic decision
                    refresh = true
                    refreshreason = :forced
                    krylovrecord[end] = with(krylovrecord[end]; escalated = true)
                else
                    krylovrecord[end] = with(krylovrecord[end]; escalationrequested = true)
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
        # the count rule: under `Always` and `Probe` one Arnoldi step is
        # enough to ask for a rebuild (which the probe may then decline),
        # under `Never` no count is
        if out.iterations >= krylovrefreshiterations
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
            # steepest descent direction -J'F of the merit function is the
            # last resort before declaring the iteration stalled.
            refreshpreconditioner!()
            out = hblinearsolve!(linearsolver, deltax, jvp, F, ws, Mop!;
                rtol = forcing, atol = gmresatol,
                maxrestarts = krylovmaxrestarts, oncycle = oncycle)
            # the rescue is often the most informative solve of the step,
            # and its Arnoldi steps count against the work budget like any
            # other solve's
            harvestafter!(out)
            record!(out, :rescue, true, false)
            work += out.iterations
            rmul!(deltax, -1)
            dϕ0dα = meritslope!(Jv, jvp, deltax, F, ϕ0,
                get(out, :residualvector, nothing))
            if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
                tr.reason = :linesearch
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
        tracetrial!(tr, alpha1, backtracks)
        if !isempty(krylovrecord) && krylovrecord[end].iteration == n
            krylovrecord[end] = with(krylovrecord[end];
                slope = normF[end] > 0 ? dϕ0dα/normF[end]^2 : NaN,
                alpha = alpha1, backtracks = backtracks, armijo = accepted)
        end

        if iszero(alpha1)
            # no decrease anywhere along the direction; F again holds the
            # residual at the unchanged x (the linesearch restore contract).
            # retry once from a fresh preconditioner, then give up
            if refreshedforstall
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
        # projection stops it once its own rate says the budget cannot
        # suffice, after one recovery
        if tracestalled(tr, progressstart, iterations - n)
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

    return IterationInfo(tr, label, krylovrecord)
end
