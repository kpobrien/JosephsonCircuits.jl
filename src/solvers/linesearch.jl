# The backtracking line search shared by the Newton and Newton-Krylov loops:
# quadratic and cubic trial steps with safeguards, the trial point and its
# evaluation, and the Armijo acceptance.

"""
    quadratic_trial_step(ϕ0, ϕ1, dϕ0dα; c1 = 1e-4, safeguard = 0.1)

Return a tuple `(αfit, ϕfit, measured)` with the proposed step `αfit`, the
estimated merit function value `ϕfit`, and `measured` a boolean indicating if
the function value is based on an evaluation of the merit function (vs an 
estimate) that minimizes a quadratic function fitted to `ϕ(α) = f(xₖ + α pₖ)`
in the range `[0, 1]`. The fitting process uses the merit function values at
`α = 0`, `α = 1`, and the derivative at the first point `dϕ(α)/dα|α = 0`.
If the full step `α = ϕα1` satisfies the Armijo sufficient-decrease condition
`ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0`, then the full step is returned without
fitting. By default `c1 = 1e-4`.

Based on Nocedal and Wright, chapter 3 section 5.

# Arguments
-`ϕ0`: `ϕ(0)`, the value of the merit function at `α = 0`.
-`ϕ1`: `ϕ(1)`, the value of the merit function at `α = 1`.
-`dϕ0dα`:  `dϕ(α)/dα|α=0`, the derivative of the merit function with respect
    to `α` at `α = 0`.

# Keywords
- `c1 = 1e-4`: the constant in the Armijo sufficient-decrease check which
     is typically (heuristicaly) set to be 1e-4,
    `ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0`.
- `safeguard = 0.1`: the smallest value we allow the step to take. This
    protects against large (eg. order of magnitude) reductions in the step
    size without an additional function evaluation, which would occur outside
    of this function.

# Returns
- `αtrial`: `αtrial` is the trial step predicted to minimize the merit
    function based on quadratic interpolation.
- `ϕtrial`: `ϕtrial` is either the predicted or measured value of the merit
    function at the trial step above. If `measured = false`, then the
    linesearch function needs to evaluate the trial point to verify that
    Armijo sufficient-decrease condition is satisfied before accepting the
    step.
- `measured`: `true` if the returned `ϕtrial` has been measured and `false` if it
    is an estimate value based on a fit.
"""
function quadratic_trial_step(ϕ0, ϕ1, dϕ0dα; c1 = 1e-4, safeguard = 0.1)
    T = float(promote_type(typeof(ϕ0),typeof(ϕ1),typeof(dϕ0dα)))
    ϕ0, ϕ1, dϕ0dα = T(ϕ0), T(ϕ1), T(dϕ0dα)
    safeguard = T(safeguard)
    c1 = T(c1)

    # check that safeguard is in (0,0.5)
    if !(zero(T) < safeguard < one(T)/2)
        throw(ArgumentError(lazy"`safeguard` = $(safeguard) must be in (0,0.5)."))
    end

    # check that c1 is in (0,1)
    if !(zero(T) < c1 < one(T))
        throw(ArgumentError(lazy"`c1` = $(c1) must be in (0,1)."))
    end

    # if the function at alpha=0 is finite there isn't much we can do since
    # ϕ0 is required for the algorithm.
    if !isfinite(ϕ0)
        throw(ArgumentError(lazy"`ϕ0` = $(ϕ0) must be finite."))
    end

    # check that the slope is negative.
    if !isfinite(dϕ0dα) || dϕ0dα >= zero(T)
        throw(ArgumentError(lazy"`dϕ0dα` = $(dϕ0dα) must be finite and negative."))
    end

    # check the Armijo sufficient decrease condition.
    # if satified, return the full step.
    if isfinite(ϕ1) && ϕ1 <= muladd(c1,dϕ0dα,ϕ0)
        return one(T), ϕ1, true
    end

    if !isfinite(ϕ1)
        # the residual at the full step overflowed, so the full step is too
        # large. halve the step. Estimate the function value from a linear
        # fit.
        return one(T)/2, muladd(one(T)/2,dϕ0dα,ϕ0), false
    end

    # coefficients of the quadratic equation ϕ(α) = a α² + b α + c to
    # interpolate ϕ(α) vs α. dϕ0dα is negative so be careful to subtract the
    # two ϕ's first to minimize loss of precision.
    a = (ϕ1-ϕ0)-dϕ0dα
    b = dϕ0dα
    c = ϕ0

    # compute the fitted value of alpha and phi. clamp it such that if the
    # fitted step is below the minimum step, take the minimum step and return
    # the fitted function value at that step. clamp to 0.5 on the upper side
    # in case floating point errors push it above 1/(2*(1-c1)) ≈ 0.5 for
    # c1 = 1e-4.
    αtrial = clamp(-(b/2)/a, safeguard, one(T) / 2)
    ϕtrial = muladd(αtrial, muladd(a, αtrial, b), c)

    return αtrial, ϕtrial, false
end

"""
    cubic_trial_step(α0, α1, ϕ0, ϕα0, ϕα1, dϕ0dα; c1 = 1e-4,
        safeguard_low = 0.1, safeguard_high = 0.5)

Return a tuple `(αfit, ϕfit, measured)` with the proposed step `αfit`, the
estimated merit function value `ϕfit`, and `measured` a boolean indicating if
the function value is based on an evaluation of the merit function (vs an 
estimate) that minimizes a cubic function fitted to `ϕ(α) = f(xₖ + α pₖ)` in
the range `[0, α1]`. The fitting process uses the merit function values at
`α = 0`, `α = α0`, `α = α1`, and the derivative at the first point
`dϕ(α)/dα|α = 0`. If the full step `α = α1` satisfies the α-scaled Armijo
sufficient-decrease condition `ϕ(α1) <= ϕ(0) + α1 c1 (dϕ(α)/dα|α = 0)`, then
the full step `α1` is returned without fitting. By default `c1 = 1e-4`.

Based on Nocedal and Wright, chapter 3 section 5.

# Arguments
-`α0`: the previous trial step `α = α0`.
-`α1`: the proposed full trial step `α = α1`.
-`ϕ0`: the value of the function at `α = 0`.
-`ϕα0`: `ϕ(α0)`, the value of the merit function at `α = α0`.
-`ϕα1`: `ϕ(α1)`, the value of the merit function at `α = α1`.
-`dϕ0dα`: `dϕ(α)/dα|α=0`, the derivative of the merit function with respect to
    `α` at `α = 0`.

# Keywords
- `c1 = 1e-4`: the constant in the α-scaled Armijo sufficient-decrease
    condition which is typically (heuristicaly) set to be 1e-4,
    `ϕ(α) <= ϕ(0) + α c1 (dϕ(α)/dα|α = 0)`.
- `safeguard_low = 0.1`: the smallest value we allow the step to take relative
    to the full step `α1`. This protects against large (eg. order of magnitude)
    reductions in the step size without an additional function evaluation,
    which would occur outside of this function.
- `safeguard_high = 0.5`: the largest value we allow the step to take relative
    to the full step `α1`. This forces the linesearch to at least reduce the
    step size by a factor of two for every backtrack. Note that if the full
    step `α=α1` satisfies the Armijo sufficient-decrease condition, then it is
    returned without any clamping.

# Returns
- `αtrial`: the trial step predicted to minimize the merit
    function based on cubic interpolation.
- `ϕtrial`: either the predicted or measured value of the merit
    function at the trial step above. If `measured = false`, then the
    linesearch function needs to evaluate the trial point to verify that
    Armijo sufficient-decrease condition is satisfied before accepting the
    step.
- `measured`: `true` if the returned `ϕtrial` has been measured (only happens
    when `α=ϕα0`) and `false` if it is an estimate value based on a fit.
"""
function cubic_trial_step(α0, α1, ϕ0, ϕα0, ϕα1, dϕ0dα; c1 = 1e-4,
    safeguard_low = 0.1, safeguard_high = 0.5)
    T = float(promote_type(typeof(α0),typeof(α1),typeof(ϕ0),typeof(ϕα0),
        typeof(ϕα1),typeof(dϕ0dα)))
    α0, α1 = T(α0), T(α1)
    ϕ0, ϕα0, ϕα1, dϕ0dα = T(ϕ0), T(ϕα0), T(ϕα1), T(dϕ0dα)
    safeguard_low = T(safeguard_low)
    safeguard_high = T(safeguard_high)
    c1 = T(c1)

    # check safeguard_low
    if !(zero(T) < safeguard_low < safeguard_high)
        throw(ArgumentError(lazy"`safeguard_low` = $(safeguard_low) must satisfy `0 < safeguard_low < safeguard_high`."))
    end

    # check safeguard_high
    if !(safeguard_high < one(T))
        throw(ArgumentError(lazy"`safeguard_high` = $(safeguard_high) must satisfy `safeguard_low < safeguard_high < 1`."))
    end

    # check that c1 is in (0,1)
    if !(zero(T) < c1 < one(T))
        throw(ArgumentError(lazy"`c1` = $(c1) must be in (0,1)."))
    end

    # if the function at alpha=0 is finite there isn't much we can do since
    # ϕ0 is required for the algorithm.
    if !isfinite(ϕ0)
        throw(ArgumentError(lazy"`ϕ0` = $(ϕ0) must be finite."))
    end

    # finite older trial step function value ϕα0 is required for the cubic
    # fit. we should screen out non-finite values before this function. 
    if !isfinite(ϕα0)
        throw(ArgumentError(lazy"`ϕα0` = $(ϕα0) must be finite."))
    end

    # check that the slope is negative.
    if !isfinite(dϕ0dα) || dϕ0dα >= zero(T)
        throw(ArgumentError(lazy"`dϕ0dα` = $(dϕ0dα) must be finite and negative."))
    end

    # the steps must be ordered as described below. this isn't strictly
    # necessary for the fitting, but since the intended use is of this
    # function is for a backtracking linesearch, we will enforce that the
    # most recent trial step is smaller than the previous one.
    if !(α1 < α0 <= one(T))
        throw(ArgumentError(lazy"`α0` = $(α0) must satisfy `α1 < α0 <= 1`."))
    end
    if !(zero(T) < α1)
        throw(ArgumentError(lazy"`α1` = $(α1) must satisfy `0 < α1 < α0`."))
    end

    # check the Armijo sufficient decrease condition for α1.
    # if satified, return the proposed step α1.
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
        c1 = 1e-4, safeguard_low = 0.1, safeguard_high = 0.5,
        maxbacktracks = 10, correction = nothing, beta = 1.0,
        Fbest = copy(F), ϕfullstep = nothing, interpolate = true)

Backtracking line search on the curvilinear trial path:

    `x(α) = x0 + α*deltax - beta*α²*correction`

with objective `ϕ(α) = 0.5*||F(x(α))||²`, following Nocedal & Wright section
3.5 with the addition of a curvilinear path. The `α²`term is a correction to
the approximate Jacobian which improves convergence particularly for strongly
driven 3WM problems. The `α²` scaling enables it to turn off at `α=0` and not
change the merit function or its derivatve at the starting point (so the 
definition of the Armijo condition is not changed). We currently compute
`correction` using Anderson acceleration (Anderson mixing). When
`correction == nothing` or `beta == 0` the path is the straight path
`x + α*deltax`.

When calling [`backtracking_linesearch!`](@ref), `F` should either hold the
residual at `x` (the same residual from which `ϕ0` was computed) or that
residual should be provided with the kwarg `Fbest`. The initial value of
`Fbest` is saved and restored if the linesearch fails to find a point
satisfying the sufficient-decrease condition or the residuals of the best
found point if the line search is successful.

[`quadratic_trial_step`](@ref) performs a quadratic interpolation on the
full-step data to estimate the trial step `α` at which the minimum of the
merit function occurs. The full step data consists of the merit function value
`ϕ0` and derivative `dϕ0dα` at the starting point `α=0` and the merit function
value `ϕfullstep` at the full step `α=1`. If `ϕfullstep` is not provided by
the user, then it is computed before calling [`quadratic_trial_step`](@ref)).
If [`quadratic_trial_step`](@ref) returns a full step with the `measured`
Boolean set to true, then we know it has already passed the Armijo
sufficient-decrease condition `ϕα <= ϕ0 + c1*α*dϕ0dα` and can be used as the
step. Return this step `α` and exit the function.

Otherwise loop over proposed trial step evaluations and cubic interpolations
with [`cubic_trial_step`](@ref). Once a successful trial step is identified
return that or return the best identified once `maxbacktracks` is reached.
With `interpolate = false` neither fit is made: the first backtrack is to
`α = 1/2` and every later one multiplies `α` by `safeguard_high`, with the
Armijo test at each trial. That is the choice for a caller whose residual
evaluations are cheap next to the direction they test.

This function always leaves`xcandidate == x(α)` and `F` holds the residual
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
    ϕ0::Real, dϕ0dα::Real; c1 = 1e-4, safeguard_low = 0.1,
    safeguard_high = 0.5, maxbacktracks::Integer = 10,
    correction::Union{Nothing, AbstractVector} = nothing, beta::Real = 1.0,
    Fbest::AbstractVector = copy(F),
    ϕfullstep::Union{Nothing, Real} = nothing, interpolate::Bool = true)

    if maxbacktracks < 0
        throw(ArgumentError(lazy"`maxbacktracks` = $(maxbacktracks) must be nonnegative."))
    end
    if !isfinite(beta) || beta < zero(beta)
        throw(ArgumentError(lazy"`beta` = $(beta) must be finite and nonnegative."))
    end

    # First take a full step, unless the merit function value ϕfullstep is
    # already provided (with F and xcandidate left at that full step), in
    # which case the evaluation is not repeated. linesearchevaluate! returns
    # the merit function value and overwrites F and xcandidate.
    # quadratic_trial_step will later validate ϕ0 and dϕ0dα and check
    # the Armijo condition at α = 1.
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

    # run the quadratic trial step function, if it returns the accepted full
    # step return that, we're done. otherwise, we will have to validate the
    # proposed step. we will use the cubic trial step function for that
    # validation since it checks the proposed step before fitting.
    #
    # The fit is only informative when the full step lands where the merit
    # is still close to quadratic. When the full step overshoots badly the
    # fitted minimizer falls below the floor and the floor is what gets
    # taken, with nothing measured in between: on a long pumped line the
    # merit along the Newton direction is minimized near 0.2 and the floor
    # step of 0.1 is accepted sixty times in a row. A caller whose trials
    # are cheap next to its directions therefore halves instead
    # (`interpolate = false`): two evaluations reach 0.25, and a trial costs
    # the Newton-Krylov path about as much as one Arnoldi step.
    α, ϕpred, accepted = if interpolate
        quadratic_trial_step(ϕ0, ϕ1, dϕ0dα; c1 = c1, safeguard = safeguard_low)
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
                ϕα, dϕ0dα; c1 = c1, safeguard_low = safeguard_low,
                safeguard_high = safeguard_high)
            # if the trial point α was accepted, then the linesearch is
            # successful and we can return that trial step. otherwise, test
            # the next proposed step αnext in the next iteration of the loop.
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

    # if the backtrack counter backtracks reaches maxbacktracks then the
    # linesearch has failed and we need to restore the best trial's residual
    # from the saved copy and recompute its trial point, xcandidate. if the
    # bestα is from the last evaluation then there is nothign we need to do.
    if bestα != αeval
        copyto!(F, Fbest)
        linesearchtrialpoint!(xcandidate, x0, bestα, deltax, beta, correction)
    end
    return bestα, bestϕ, false, backtracks
end
