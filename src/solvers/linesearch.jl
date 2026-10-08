# The backtracking line search shared by the Newton and Newton-Krylov loops:
# quadratic and cubic trial steps with safeguards, the trial point and its
# evaluation, and the Armijo acceptance.

"""
    quadratic_trial_step(ϕ0, ϕ1, dϕ0dα, ls::Backtracking)

Return a tuple `(αtrial, ϕtrial, measured)`: the step `αtrial` which
minimizes a quadratic fitted to the merit function `ϕ(α) = f(xₖ + α pₖ)` on
`[0, 1]`, the merit `ϕtrial` there, and `measured`, whether `ϕtrial` was
evaluated rather than estimated from the fit. The fit takes the merit at
`α = 0` and at `α = 1` and its derivative `dϕ(α)/dα|α = 0`. A full step
`α = 1` which satisfies the Armijo sufficient-decrease condition
`ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0` is returned without a fit.

Based on Nocedal and Wright, chapter 3 section 5.

# Arguments
- `ϕ0`: `ϕ(0)`, the value of the merit function at `α = 0`, finite.
- `ϕ1`: `ϕ(1)`, the value of the merit function at `α = 1`.
- `dϕ0dα`: `dϕ(α)/dα|α=0`, the derivative of the merit function with
    respect to `α` at `α = 0`, finite and negative.
- `ls`: the [`Backtracking`](@ref) which gives the constant `c1` of the
    Armijo condition and the safeguards which bound the proposed step. The
    lower one, `safeguardlow`, protects against large (eg. order of magnitude)
    reductions in the step size without an additional function evaluation,
    which would occur outside of this function; the fitted minimizer of a
    full step which fails the Armijo condition is below `1/(2(1 - c1))`, so
    the upper one, `safeguardhigh`, acts only when it is set below that.

The callers check `ϕ0` and `dϕ0dα` before any search, and `ls` checked its
settings when it was built, so neither is checked here.

# Returns
- `αtrial`: the step the quadratic fit predicts to minimize the merit
    function, within the safeguards, or the full step.
- `ϕtrial`: the merit at `αtrial`, measured or predicted. A predicted one
    is the fit's: the line search evaluates the trial point and tests the
    Armijo condition there before accepting the step.
- `measured`: `true` when the full step is returned, having satisfied the
    Armijo condition, so that `ϕtrial` is its measured merit, and `false`
    when `ϕtrial` is an estimate.
"""
function quadratic_trial_step(ϕ0, ϕ1, dϕ0dα, ls::Backtracking)
    T = float(promote_type(typeof(ϕ0),typeof(ϕ1),typeof(dϕ0dα)))
    ϕ0, ϕ1, dϕ0dα = T(ϕ0), T(ϕ1), T(dϕ0dα)
    safeguard_low = T(ls.safeguardlow)
    safeguard_high = T(ls.safeguardhigh)
    c1 = T(ls.c1)

    # the Armijo sufficient decrease condition, under which the full step
    # is returned
    if isfinite(ϕ1) && ϕ1 <= muladd(c1,dϕ0dα,ϕ0)
        return one(T), ϕ1, true
    end

    if !isfinite(ϕ1)
        # the residual at the full step overflowed, so the full step is too
        # large: halve the step, within the safeguards, and estimate the
        # function value from a linear fit
        αtrial = clamp(one(T)/2, safeguard_low, safeguard_high)
        return αtrial, muladd(αtrial, dϕ0dα, ϕ0), false
    end

    # coefficients of the quadratic equation ϕ(α) = a α² + b α + c to
    # interpolate ϕ(α) vs α. dϕ0dα is negative so be careful to subtract the
    # two ϕ's first to minimize loss of precision.
    a = (ϕ1-ϕ0)-dϕ0dα
    b = dϕ0dα
    c = ϕ0

    # compute the fitted value of alpha and phi, clamped to the safeguards:
    # a fitted step outside them is replaced by the bound and the fitted
    # function value is returned at that step. The minimizer is below
    # 1/(2*(1-c1)) ≈ 0.5 for c1 = 1e-4, so an upper safeguard of one half
    # only catches the floating point error which pushes it above.
    αtrial = clamp(-(b/2)/a, safeguard_low, safeguard_high)
    ϕtrial = muladd(αtrial, muladd(a, αtrial, b), c)

    return αtrial, ϕtrial, false
end

"""
    cubic_trial_step(α0, α1, ϕ0, ϕα0, ϕα1, dϕ0dα, ls::Backtracking)

Return a tuple `(αtrial, ϕtrial, measured)`: the step `αtrial` which
minimizes a cubic fitted to the merit function `ϕ(α) = f(xₖ + α pₖ)` on
`[0, α1]`, the merit `ϕtrial` there, and `measured`, whether `ϕtrial` was
evaluated rather than estimated from the fit. The fit takes the merit at
`α = 0`, `α = α0` and `α = α1` and its derivative `dϕ(α)/dα|α = 0`. A latest
trial `α1` which satisfies the α-scaled Armijo sufficient-decrease
condition `ϕ(α1) <= ϕ(0) + α1 c1 (dϕ(α)/dα|α = 0)` is returned without a
fit.

Based on Nocedal and Wright, chapter 3 section 5.

# Arguments
- `α0`: the trial before the latest, `0 < α1 < α0 <= 1`.
- `α1`: the latest trial.
- `ϕ0`: the value of the merit function at `α = 0`, finite.
- `ϕα0`: `ϕ(α0)`, the value of the merit function at `α = α0`, finite.
- `ϕα1`: `ϕ(α1)`, the value of the merit function at `α = α1`.
- `dϕ0dα`: `dϕ(α)/dα|α=0`, the derivative of the merit function with
    respect to `α` at `α = 0`, finite and negative.
- `ls`: the [`Backtracking`](@ref) which gives the constant `c1` of the
    α-scaled Armijo condition `ϕ(α) <= ϕ(0) + α c1 (dϕ(α)/dα|α = 0)` and
    the safeguards which bound the proposed step, relative to the latest
    trial `α1`: `safeguardlow` the smallest, which bounds the cut a single
    fit makes before the caller evaluates the merit at the step it
    proposes, and `safeguardhigh` the largest, so that every backtrack
    cuts the step by at least this factor. A trial `α1` which satisfies
    the Armijo condition is returned unclamped.

The search makes its trials in this order and hands only finite merits as
`ϕα0`, the callers check `ϕ0` and `dϕ0dα` before any search, and `ls`
checked its settings when it was built, so none of these is checked here.

# Returns
- `αtrial`: the step the cubic fit predicts to minimize the merit
    function, within the safeguards.
- `ϕtrial`: the merit at `αtrial`, measured or predicted. A predicted one
    is the fit's: the line search evaluates the trial point and tests the
    Armijo condition there before accepting the step.
- `measured`: `true` when the latest trial `α1` itself is returned, having
    satisfied the Armijo condition, so that `ϕtrial` is its measured merit,
    and `false` when `ϕtrial` is the fit's estimate.
"""
function cubic_trial_step(α0, α1, ϕ0, ϕα0, ϕα1, dϕ0dα, ls::Backtracking)
    T = float(promote_type(typeof(α0),typeof(α1),typeof(ϕ0),typeof(ϕα0),
        typeof(ϕα1),typeof(dϕ0dα)))
    α0, α1 = T(α0), T(α1)
    ϕ0, ϕα0, ϕα1, dϕ0dα = T(ϕ0), T(ϕα0), T(ϕα1), T(dϕ0dα)
    safeguard_low = T(ls.safeguardlow)
    safeguard_high = T(ls.safeguardhigh)
    c1 = T(ls.c1)

    # the Armijo sufficient decrease condition at α1, which is returned
    # when it holds
    if isfinite(ϕα1) && ϕα1 <= muladd(c1*α1,dϕ0dα,ϕ0)
        return α1, ϕα1, true
    end

    fallback = safeguard_high*α1
    if !isfinite(ϕα1)
        # the residual at the full step overflowed, so the full step is too
        # large. propose safeguard_high*α1. Estimate the function value from
        # a linear fit.
        return fallback, muladd(fallback,dϕ0dα,ϕ0), false
    end

    # scaled coefficients of the cubic equation ϕ(α) = a α³ + b α² + c α + d
    # to interpolate ϕ(α) vs α.
    r0 = muladd(-α0,dϕ0dα,ϕα0 - ϕ0)/(α0*α0)
    r1 = muladd(-α1,dϕ0dα,ϕα1 - ϕ0)/(α1*α1)
    delta = α1-α0
    a = (r1-r0)/delta
    b = (α1*r0-α0*r1)/delta
    disc = muladd(b,b,-3*a*dϕ0dα)
    if disc < zero(T) || !isfinite(disc)
        α = fallback
    else
        s = sqrt(disc)
        if b >= zero(T)
            α = -dϕ0dα/(b+s)
        else
            α = ((s-b)/3)/a
        end

        if !isfinite(α)
            α = fallback
        end
    end

    # clamp the returned trial step to be in the range
    # [α1*safeguard_low,α1*safeguard_high]
    αtrial = clamp(α,α1*safeguard_low,α1*safeguard_high)
    ϕtrial = muladd(αtrial,muladd(αtrial,muladd(αtrial,a,b),dϕ0dα),ϕ0)
    return αtrial, ϕtrial, false
end

"""
    linesearchtrialpoint!(xcandidate, x, α, deltax, beta, correction)

Overwrite `xcandidate` with the line search trial points on the curvilinear
path: `xcandidate = x + α*deltax - beta*α²*correction`, or the straight path
`x + α*deltax` when `correction == nothing` or `beta == 0`. Overwrites and
returns `xcandidate`.
"""
function linesearchtrialpoint!(xcandidate::AbstractVector, x::AbstractVector,
    α::Real, deltax::AbstractVector, beta::Real,
    correction::Union{Nothing, AbstractVector})
    if iszero(α)
        # if α is zero just copy `x` into `xcandidate`.
        copyto!(xcandidate, x)
    elseif isnothing(correction) || iszero(beta)
        @. xcandidate = x + α*deltax
    else
        @. xcandidate = x + α*deltax - (beta*α*α)*correction
    end
    return xcandidate
end

"""
    merit(F)

The line search merit function `ϕ = 0.5*||F||²` = `real(0.5*dot(F, F))`. The
dot product conjugates the first argument so this can be used for real or
complex vectors.
"""
merit(F::AbstractVector) = real(dot(F, F)/2)

"""
    linesearchevaluate!(f!, F, xcandidate, x, α, deltax, beta, correction)

Evaluate the line search merit function at the trial step `α`: generate the
trial point with [`linesearchtrialpoint!`](@ref), evaluate the residual there
with `f!(F, xcandidate)`, and return [`merit`](@ref)`(F)`. `xcandidate`
and `F` contain the trial point and its residual when this function finishes.
"""
function linesearchevaluate!(f!, F::AbstractVector,
    xcandidate::AbstractVector, x::AbstractVector, α::Real,
    deltax::AbstractVector, beta::Real,
    correction::Union{Nothing, AbstractVector})
    linesearchtrialpoint!(xcandidate, x, α, deltax, beta, correction)
    f!(F, xcandidate)
    return merit(F)
end

"""
    backtracking_linesearch!(f!, F, xcandidate, x0, deltax, ϕ0, dϕ0dα;
        ls = Backtracking(), correction = nothing, beta = 1.0,
        Fbest = copy(F), ϕfullstep = nothing)

Backtracking line search on the curvilinear trial path:

    `x(α) = x0 + α*deltax - beta*α²*correction`

with objective `ϕ(α) = 0.5*||F(x(α))||²`, following Nocedal & Wright section
3.5 with the addition of a curvilinear path. The `α²` term is a correction
to the approximate Jacobian, which improves convergence on strongly driven
three wave mixing problems. Its `α²` scaling makes it vanish at `α = 0` with
its derivative, so it changes neither the merit function nor its slope at
the starting point, and the Armijo condition is that of the straight path.
`correction` is the Anderson acceleration (Anderson mixing) correction of
[`nlsolve!`](@ref). When `correction == nothing` or `beta == 0` the path is
the straight path `x + α*deltax`.

`F` holds the residual at `x0`, from which `ϕ0` was computed, and `Fbest`
is scratch: the search keeps there the residual of its best trial,
starting from the residual at `x0`, and restores it into `F` when it ends
without an accepted step. When `ϕfullstep` is given, `F` and `xcandidate`
hold the full step instead, and `Fbest` must hold the residual at `x0`.

[`quadratic_trial_step`](@ref) performs a quadratic interpolation on the
full-step data to estimate the trial step `α` at which the minimum of the
merit function occurs. The full step data consists of the merit function value
`ϕ0` and derivative `dϕ0dα` at the starting point `α=0` and the merit function
value `ϕfullstep` at the full step `α=1`, evaluated first when the caller
does not provide it. A full step which [`quadratic_trial_step`](@ref)
returns as `measured` has passed the Armijo sufficient-decrease condition
`ϕα <= ϕ0 + c1*α*dϕ0dα` and is returned.

Otherwise trial evaluations alternate with the cubic fits of
[`cubic_trial_step`](@ref) until a trial satisfies the condition, which is
returned, or `maxbacktracks` trials have been made, when the best of them
is.
The Armijo constant `c1`, the safeguards of the fits, the trial budget
`maxbacktracks` and whether the fits are made at all are the fields of
`ls`, a [`Backtracking`](@ref). Without interpolation neither fit is made:
the first backtrack is to `α = 1/2` and every later one multiplies `α` by
`safeguardhigh`, with the Armijo test at each trial.

This function always leaves `xcandidate == x(α)` and `F` holds the residual
there (for `α == 0` that is the residual at `x0`).

Returns `(α, ϕα, accepted, backtracks)`:
- `accepted == true`: `α` satisfies the α-scaled Armijo condition
  `ϕα <= ϕ0 + c1*α*dϕ0dα`, and `ϕα` is its measured objective value.
- `accepted == false`: maxbacktracks was reached. `α` is the best
  (lowest measured ϕ) trial found, which may still be a useful step; if no
  trial produced any decrease at all, `α == 0` and `ϕα == ϕ0`.
- `backtracks` counts the trial evaluations after the full step, so the
  total number of `f!` residual evaluations is exactly `backtracks + 1`:
  the failure path restores the best trial's residual from a copy saved
  when it was measured (`Fbest`, a caller-suppliable scratch buffer whose
  contents are clobbered), never by re-evaluating.
"""
function backtracking_linesearch!(f!, F::AbstractVector,
    xcandidate::AbstractVector, x0::AbstractVector, deltax::AbstractVector,
    ϕ0::Real, dϕ0dα::Real; ls::Backtracking = Backtracking(),
    correction::Union{Nothing, AbstractVector} = nothing, beta::Real = 1.0,
    Fbest::AbstractVector = copy(F),
    ϕfullstep::Union{Nothing, Real} = nothing)

    if !isfinite(beta) || beta < zero(beta)
        throw(ArgumentError(lazy"`beta` = $(beta) must be finite and nonnegative."))
    end
    # the settings of the search, validated when `ls` was built
    c1 = ls.c1
    safeguard_high = ls.safeguardhigh
    maxbacktracks = ls.maxbacktracks
    interpolate = ls.interpolate

    # First take a full step, unless the merit function value ϕfullstep is
    # already provided (with F and xcandidate left at that full step), in
    # which case the evaluation is not repeated. linesearchevaluate! returns
    # the merit function value and overwrites F and xcandidate.
    # quadratic_trial_step will later check the Armijo condition at α = 1.
    ϕ1 = if isnothing(ϕfullstep)
        # F holds the residual at x0. save it so the no-decrease
        # failure path can restore it by copy (a caller passing ϕfullstep
        # must pre-populate Fbest with the residual at x0 itself)
        copyto!(Fbest, F)
        linesearchevaluate!(f!, F, xcandidate, x0, one(ϕ0), deltax, beta,
            correction)
    else
        ϕfullstep
    end

    # the quadratic fit through the full step proposes the first trial, and
    # an accepted full step returns at once; otherwise the cubic fit
    # validates each proposal before making the next. The fit is only
    # informative where the merit is still close to quadratic: when the
    # full step overshoots badly the fitted minimizer falls to the floor
    # with nothing measured in between, which is what halving avoids for a
    # caller whose trials are cheap next to the direction they test.
    α, ϕpred, accepted = if interpolate
        quadratic_trial_step(ϕ0, ϕ1, dϕ0dα, ls)
    else
        armijo = isfinite(ϕ1) && ϕ1 <= muladd(c1, dϕ0dα, ϕ0)
        (armijo ? one(ϕ0) : one(ϕ0)/2, ϕ1, armijo)
    end
    if accepted
        return α, ϕ1, true, 0
    end

    # check if the merit function value at the full step ϕ1 is finite. if it
    # is finite, store as (αprev, ϕprev) and use them in the cubic fit. update
    # and store bestα, bestϕ, Fbest to return in case of backtracking failure.
    havefinite = isfinite(ϕ1)
    αprev, ϕprev = one(α), ϕ1
    bestα, bestϕ = zero(α), ϕ0
    if havefinite && ϕ1 < bestϕ
        bestα, bestϕ = one(α), ϕ1
        copyto!(Fbest, F)
    end

    # run the backtracking loop
    backtracks = 0
    # F is currently at the full step, α=1
    αeval = one(α)
    while backtracks < maxbacktracks
        backtracks += 1
        # evaluate the merit function for the proposed trial step
        ϕα = linesearchevaluate!(f!, F, xcandidate, x0, α, deltax, beta,
            correction)
        αeval = α
        if isfinite(ϕα) && ϕα < bestϕ
            bestα, bestϕ = α, ϕα
            copyto!(Fbest, F)
        end
        # if the merit function is finite at the trial step, try the cubic fit
        if havefinite && interpolate
            # cubic fit through (0, ϕ0, dϕ0dα) and the two most recent
            # trials. first performs the Armijo test at α.
            αnext, ϕpred, accepted = cubic_trial_step(αprev, α, ϕ0, ϕprev,
                ϕα, dϕ0dα, ls)
            # an accepted trial ends the search; otherwise the next pass
            # of the loop tests the proposed step αnext
            if accepted
                return α, ϕα, true, backtracks
            end
        else
            # no finite older trial yet (every longer step overflowed), or
            # the caller asked for no interpolation: test Armijo directly
            # and back off geometrically.
            if isfinite(ϕα) && ϕα <= muladd(c1*α, dϕ0dα, ϕ0)
                return α, ϕα, true, backtracks
            end
            αnext = safeguard_high*α
        end
        # slide the window over finite measurements only, so the cubic
        # never receives a non-finite ϕprev (which it would reject).
        if isfinite(ϕα)
            αprev, ϕprev = α, ϕα
            havefinite = true
        end
        α = αnext
    end

    # with the trials spent the search has failed: the best trial's
    # residual is restored from its saved copy and its trial point
    # recomputed, unless the best trial is the last one evaluated, which F
    # and xcandidate hold already
    if bestα != αeval
        copyto!(F, Fbest)
        linesearchtrialpoint!(xcandidate, x0, bestα, deltax, beta, correction)
    end
    return bestα, bestϕ, false, backtracks
end
