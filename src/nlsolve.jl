

"""
    Factorization(outofplace,inplace,kwargs)

A structure to hold the factorizations and their keyword arguments.

```
"""
struct Factorization
    factorize
    factorize!
    kwargs
end

function KLUfactorization(;kwargs...)
    return Factorization(KLU.klu,KLU.klu!,kwargs)
end

function LUfactorization(;kwargs...)
    return Factorization(lu,lu!,kwargs)
end

function QRfactorization(;kwargs...)
    return Factorization(qr,nothing,kwargs)
end

# fallback method to handle everything but QR adjoint
function myldiv!(x,F,b)
    return ldiv!(x,F,b)
end

# workaround for lack of QR adjoint factorization
# modification of solution proposed here
# and https://github.com/JuliaSparse/SparseArrays.jl/issues/656
# explanation: https://www.netlib.org/lapack/lug/node41.html
# and based in part on ldiv! in spqr.jl to handle the rank deficient case
# https://github.com/JuliaSparse/SparseArrays.jl/blob/main/src/solvers/spqr.jl
function myldiv!(x::StridedVecOrMat{<:Number},Fadj::LinearAlgebra.AdjointFactorization{<:Number, <:SparseArrays.SPQR.QRSparse},
                b::StridedVecOrMat{<:Number})
    F = parent(Fadj)
    m, n = size(F)
    pcol = F.pcol
    prow = F.prow
    rnk = rank(F)

    nrhs = size(b, 2)

    # define a workspace
    W = Matrix{eltype(b)}(undef, m, nrhs)

    # c[1:rnk] = first rnk permuted columns of b with the rest of c zero.
    @inbounds for j in 1:nrhs
        for i in 1:rnk
            W[i, j] = b[pcol[i], j]
        end
        for i in (rnk+1):m
            W[i, j] = zero(eltype(W))
        end
    end

    # solve R11^H * y[1:rnk] = c[1:rnk]
    if rnk > 0
        R11 = F.R[1:rnk, 1:rnk]
        ldiv!(UpperTriangular(R11)', view(W, 1:rnk, :))
    end

    # W = Q*y
    lmul!(F.Q, W)

    # get the solution
    @inbounds for j in 1:nrhs
        for i in 1:m
            x[prow[i], j] = W[i, j]
        end
    end

    return x
end

"""
    FactorizationCache(factorization)

A cache for the factorization object.

# Examples
```jldoctest
julia> JosephsonCircuits.FactorizationCache(JosephsonCircuits.KLU.klu(JosephsonCircuits.sparse([1, 2], [1, 2], [1/2, 1/2], 2, 2)));

```
"""
mutable struct FactorizationCache
    factorization
end

function FactorizationCache()
    return FactorizationCache(nothing)
end

"""
    tryfactorize!(cache::FactorizationCache,
        factorization::Factorization, A::AbstractArray)

Factorize the matrix `A` using the factorization from `factorization` and
store the result in `cache`. Attempt to reuse the symbolic factorization. Redo the
symbolic factorization if we get a SingularException.

"""
function tryfactorize!(cache::FactorizationCache,
    factorization::Factorization, A::AbstractMatrix)

    # if the factorization cache is empty then generate a factorizaion
    if isnothing(cache.factorization)
        cache.factorization = factorization.factorize(A;
            factorization.kwargs...)
    # otherwise, try to update the factorization, falling back to generating
    # a new one if that fails
    elseif isnothing(factorization.factorize!)
        cache.factorization = factorization.factorize(A;
            factorization.kwargs...)
    else
        try
            # update the factorization. the sparsity structure does 
            # not change so we can reuse the factorization object.
            factorization.factorize!(cache.factorization, A;
                factorization.kwargs...)
        catch e
            if isa(e, SingularException)
                # reusing the symbolic factorization can sometimes
                # lead to numerical problems. if the first linear
                # solve fails try factoring and solving again
                cache.factorization = factorization.factorize(A;
                    factorization.kwargs...)
            else
                throw(e)
            end
        end
    end
    return cache
end


"""
    trysolve!(x,factorization,b)

First try to solve a linear system using ldiv! then if it errors, use \\. The
motivation for this function is some factorizations such as `qr` with sparse
matrices don't support ldiv!. 
"""
function trysolve!(x,factorization,b)
    try
        myldiv!(x,factorization,b)
    catch e
        if e isa MethodError || e isa ArgumentError
            x .= factorization \ b
        else
            rethrow()
        end
    end
    return x
end

"""
    linesearch(ϕ0, ϕ1, dϕdα, αmin; c1 = 1e-4)

Return a tuple containing the proposed step and estimated function value
(αfit,ϕfit) that minimizes a quadratic function fitted to
`ϕ(α) = f(xₖ + α pₖ)` in the range `[αmin, 1]`. The fitting processes uses the
values at `α = 0`, `α = 1`, and the derivative at the first point
`dϕ(α)/dα|α = 0`. If the full step `α = 1` satisfies the Armijo
sufficient-decrease condition `ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0`, then the
full step is returned without fitting. By default `c1 = 1e-4`.

Quadratic linesearch based on Nocedal and Wright, chapter 3 section 5.

# Arguments
-`ϕ0`: `ϕ(0)`, the value of the function at `α = 0`.
-`ϕ1`: `ϕ(1)`, the value of the function at `α = 1`.
-`dϕdα0`:  `dϕ(α)/dα|α=0`, the derivative of the function with respect to `α`
    at `α = 0`.
-`αmin`: the minimum value of `α` to return.
"""
function linesearch(ϕ0, ϕ1, dϕ0dα, αmin; c1 = 1e-4)
    αtrial, ϕtrial, measured = quadratic_trial_step(ϕ0, ϕ1, dϕ0dα; c1 = c1, safeguard = αmin)
    return αtrial, ϕtrial
end

"""
    quadratic_trial_step(ϕ0, ϕ1, dϕ0dα; c1 = 1e-4, safeguard = 0.1)

Return a tuple containing the proposed step and the estimated function value
(αfit,ϕfit) that minimizes a quadratic function fitted to
`ϕ(α) = f(xₖ + α pₖ)` in the range `[0, 1]`. The fitting process uses the
values at `α = 0`, `α = 1`, and the derivative at the first point
`dϕ(α)/dα|α = 0`. If the full step `α = 1` satisfies the Armijo
sufficient-decrease condition `ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0`, then the
full step is returned without fitting. By default `c1 = 1e-4`.

Quadratic linesearch based on Nocedal and Wright, chapter 3 section 5.

# Arguments
-`ϕ0`: `ϕ(0)`, the value of the function at `α = 0`.
-`ϕ1`: `ϕ(1)`, the value of the function at `α = 1`.
-`dϕ0dα`:  `dϕ(α)/dα|α=0`, the derivative of the function with respect to `α`
    at `α = 0`.

# Keywords
- `c1 = 1e-4`: the constant in the Armijo sufficient-decrease check which
     is typically (heuristicaly) defined to be 1e-4,
    `ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0`.
- `safeguard = 0.1`: the smallest value we allow the step to take. This
    protects against large (eg. order of magnitude) reductions in the step
    size without an additional function evaluation, which would occur outside
    of this function.

# Returns
- `αtrial`: `αtrial` is the trial step predicted to minimize the merit function based
    on quadratic interpolation.
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

Return a tuple containing the proposed step and the estimated function value
(αfit,ϕfit) that minimizes a cubic function fitted to
`ϕ(α) = f(xₖ + α pₖ)` in the range `[0, ϕkm1]`. The fitting process uses the
values at `α = 0`, `α = αkm1`, `α = αk`, and the derivative at the first point
`dϕ(α)/dα|α = 0`. If the full step `α = αkm1` satisfies the Armijo
sufficient-decrease condition `ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0`, then the
full step is returned without fitting. By default `c1 = 1e-4`.

Quadratic linesearch based on Nocedal and Wright, chapter 3 section 5.

# Arguments
-`α0`: the step `α = α0`.
-`α1`: the step `α = α1`.
-`ϕ0`: the value of the function at `α = 0`.
-`ϕα0`: `ϕ(α0)`, the value of the merit function at `α = α0`.
-`ϕα1`: `ϕ(α1)`, the value of the function at `α = α1`.
-`dϕ0dα`: `dϕ(α)/dα|α=0`, the derivative of the function with respect to `α`
    at `α = 0`.

# Keywords
- `c1 = 1e-4`: the constant in the Armijo sufficient-decrease check which
     is typically (heuristicaly) defined to be 1e-4,
    `ϕ(1) <= ϕ(0) + c1 dϕ(α)/dα|α = 0`.
- `safeguard = 0.1`: the smallest value we allow the step to take. This
    protects against large (eg. order of magnitude) reductions in the step
    size without an additional function evaluation, which would occur outside
    of this function.

# Returns
- `αtrial`: `αtrial` is the trial step predicted to minimize the merit
    function based on cubic interpolation.
- `ϕtrial`: `ϕtrial` is either the predicted or measured value of the merit
    function at the trial step above. If `measured = false`, then the
    linesearch function needs to evaluate the trial point to verify that
    Armijo sufficient-decrease condition is satisfied before accepting the
    step.
- `measured`: `true` if the returned `ϕtrial` has been measured and `false` if
    it is an estimate value based on a fit.
"""
function cubic_trial_step(α0, α1, ϕ0, ϕα0, ϕα1, dϕ0dα; c1 = 1e-4,
    safeguard_low = 0.1, safeguard_high = 0.5)
    T = float(promote_type(typeof(α0),typeof(α1),typeof(ϕ0),typeof(ϕα0),typeof(ϕα1),typeof(dϕ0dα)))
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
"""
struct IterationInfo
    label::String
    parameter::Float64
    regularization::Float64
    converged::Bool
    iterations::Int
    normresidual::Vector{Float64}
    alpha::Vector{Float64}
    backtracks::Vector{Int}
    andersonaccepted::Vector{Bool}
end

"""
    nlsolve!(fj!::Function, F::AbstractVector{T}, J::AbstractArray{T},
        x::Vector{T}; iterations=1000, ftol=1e-8, switchofflinesearchtol = 1e-5,
        alphamin = 1e-4, factorization = KLUfactorization(),
        andersondepth::Integer = 5)

A simple nonlinear solver using Newton's method with linesearch based on
Nocedal and Wright, chapter 3 section 5, combined with Anderson acceleration.

This solver attempts to find x such that f(x) == 0, where f is a
nonlinear function with Jacobian J.

A few points to note:
(1) It uses KLU factorization, so only works on sparse matrices.
(2) The Jacobian J cannot change sparsity structure.
(3) This function attempts to reuse the symbolic factorization which can
    sometimes result in a SingularException, which we catch, then create a
    new factorization object.

Returns `true` if the iterations converged and `false` otherwise.

# Arguments
- `fj!`: a function to compute a vector-valued objective function and
its Jacobian.
- `F`: vector for holding intermediate results.
- `J`: sparse matrix with with the desired sparsity structure of the
  Jacobian.
- `x`: initial values for x.

# Keywords
- `iterations = 1000`: the maximum number of iterations.
- `ftol = 1e-8`: the tolerance on the residual norm below which the
  iterations are considered converged.
- `switchofflinesearchtol = 1e-5`: the relative residual norm below which
  the linesearch is switched off and full Newton steps are taken.
- `alphamin = 1e-4`: the minimum linesearch step size.
- `factorization = KLUfactorization()`: the matrix factorization to use
  for the linear solves.
- `andersondepth::Integer = 5`: the depth of the Anderson acceleration of the
  Newton fixed point iteration, the maximum number of previous iterates
  used for the extrapolation. Values less than one disable the
  acceleration.
- `stalliterations::Integer = 0`: when positive, give up when the best
  residual norm of the most recent half of a window of this many
  iterations has improved by less than ten percent over the best of the
  half window before it, so a solve converging too slowly to finish
  within any reasonable budget hands control back to any fallbacks
  quickly instead of exhausting the iteration budget. A solve converging
  at a constant rate r is abandoned when r^(stalliterations/2) > 0.9.
  Off by default.

# Examples
```jldoctest
function fj!(F, J, x)
    if !isnothing(F)
        F[1] = (x[1]+3)*(x[2]^3-7)+18
        F[2] = sin(x[2]*exp(x[1])-1)
    end
    if !isnothing(J)
        J[1, 1] = x[2]^3-7
        J[1, 2] = 3*x[2]^2*(x[1]+3)
        u = exp(x[1])*cos(x[2]*exp(x[1])-1)
        J[2, 1] = x[2]*u
        J[2, 2] = u
    end
    return nothing
end
x = [ 0.1, 1.2]
F = [0.0, 0.0]
J = JosephsonCircuits.sparse([1, 1, 2, 2],[1, 2, 1, 2],[1.3, 0.5, 0.1, 1.2])
JosephsonCircuits.nlsolve!(fj!, F, J, x)
isapprox([0.0,1.0],x)

# output
true
```
"""
function nlsolve!(fj!::Function, F::AbstractVector{T}, J::AbstractArray{T},
    x::Vector{T}; iterations::Integer = 1000, ftol = 1e-8,
    switchofflinesearchtol = 1e-5, alphamin = 1e-4,
    factorization = KLUfactorization(), label = "quasi-Newton",
    andersondepth::Integer = 5, stalliterations::Integer = 0) where T

    if size(J,1) != size(J,2)
        throw(DimensionMismatch(lazy"The Jacobian `J` matrix must be square."))
    end

    if size(J,2) != length(x)
        throw(DimensionMismatch(lazy"Second axis of Jacobian `J` must have the same length as the input `x`."))
    end

    if size(J,1) != length(F)
        throw(DimensionMismatch(lazy"First axis of the Jacobian `J` must have the same length as the residual `F`."))
    end

    cache = FactorizationCache()
    tryfactorize!(cache,factorization,J)

    deltax = copy(x)
    # reusable buffer for trial points, to avoid allocating a new vector
    # at every trial evaluation
    xtrial = copy(x)

    # working arrays and history for Anderson acceleration. the Newton
    # update deltax = -J \\ F is the residual of the fixed point map
    # G(x) = x + deltax, so the history of iterates and updates is
    # available at no extra cost. the extrapolation coefficients are
    # constrained to be real because for the harmonic balance quasi-Newton
    # map the error operator is antilinear (it involves complex
    # conjugation), so complex coefficients cannot cancel the error modes;
    # real coefficients correspond to Anderson acceleration of the
    # equivalent real system.
    if andersondepth > 0
        xanderson = copy(x)
        fanderson = copy(deltax)
        Fanderson = copy(F)
        xcandidate = copy(x)
        andersonready = false
        deltaxhistory = Vector{typeof(x)}(undef, 0)
        deltafhistory = Vector{typeof(x)}(undef, 0)
        # buffers for the small least squares problem for the
        # extrapolation coefficients
        gram = Matrix{Float64}(undef, andersondepth, andersondepth)
        rhs = Vector{Float64}(undef, andersondepth)
    end

    # diagnostic info
    normF = Float64[]
    alphas = Float64[]
    backtrackrecord = Int[]
    andersonrecord = Bool[]


    # perform Newton's method with linesearch based on Nocedal and Wright
    # chapter 3 section 5.
    converged = false
    Ffresh = false
    updatenorm = Inf
    for _ in 1:iterations

        if Ffresh
            # the residual F was evaluated at the current x at the end of
            # the previous iteration, so only the Jacobian needs updating.
            fj!(nothing, J, x)
        else
            # update the residual function and the Jacobian
            fj!(F, J, x)
        end

        push!(normF, norm(F))

        # when stalliterations is positive, give up when the best
        # residual norm of the most recent half window has improved by
        # less than ten percent over the best of the half window before
        # it, so a solve converging too slowly to finish within any
        # reasonable budget hands control back to any fallbacks (such as
        # the retreat of a continuation) quickly instead of exhausting
        # the iteration budget. comparing the minima of two half windows
        # is robust both to the non-monotonic residual of the Anderson
        # accelerated iteration and to slow monotone convergence, for
        # which comparing the current residual against the window
        # minimum would always report a stall. a solve converging at a
        # constant rate r is abandoned when r^(stalliterations/2) > 0.9.
        if stalliterations > 0 && length(normF) > stalliterations
            half = stalliterations ÷ 2
            recentmin = minimum(@view normF[end-half+1:end])
            oldermin = minimum(@view normF[end-stalliterations+1:end-half])
            if recentmin > 0.9*oldermin
                break
            end
        end

        # factor the Jacobian
        tryfactorize!(cache,factorization,J)

        # solve the linear system
        trysolve!(deltax, cache.factorization, F)

        # multiply deltax by -1
        rmul!(deltax, -1)

        # calculate the objective function and the derivative of the objective
        # with respect to the scalar variable alpha which parameterizes the
        # path between the old x and the new x. 
        # Note: the dot product takes the complex conjugate of the first vector
        f = real(0.5*dot(F, F))
        dfdalpha = real(dot(F, J, deltax))

        # evaluate the function at the trial point
        @. xtrial = x + deltax
        fj!(F, nothing, xtrial)

        fp = real(0.5*dot(F,F))

        # calculate the step size based on the last point, the trial point, and
        # derivative at the first point.
        alpha1, _ = linesearch(f,fp,dfdalpha,alphamin)

        # switch to newton once the norm is small enough
        normx = norm(x)
        if normx > 0 && sqrt(fp)/normx <= switchofflinesearchtol && sqrt(f)/normx <= switchofflinesearchtol
            alpha1 = 1.0
        end

        # Anderson acceleration: extrapolate through the history of
        # iterates and fixed point residuals, and accept the extrapolated
        # iterate only if it reduces the norm of the residual, so the
        # accelerated iteration is never less robust than the damped
        # iteration.
        andersonaccepted = false
        if andersondepth > 0
            if andersonready
                # reuse the oldest history vectors when the history is full,
                # to avoid allocating new ones at every iteration
                if length(deltaxhistory) >= andersondepth
                    deltaxnew = popfirst!(deltaxhistory)
                    deltafnew = popfirst!(deltafhistory)
                else
                    deltaxnew = similar(x)
                    deltafnew = similar(x)
                end
                @. deltaxnew = x - xanderson
                @. deltafnew = deltax - fanderson
                push!(deltaxhistory, deltaxnew)
                push!(deltafhistory, deltafnew)
            end
            copyto!(xanderson, x)
            copyto!(fanderson, deltax)
            andersonready = true
            m = length(deltafhistory)
            if m > 0
                # solve the small least squares problem for the real
                # extrapolation coefficients with regularized normal
                # equations
                for i in 1:m
                    for j in i:m
                        gram[i, j] = real(dot(deltafhistory[i], deltafhistory[j]))
                        gram[j, i] = gram[i, j]
                    end
                    rhs[i] = real(dot(deltafhistory[i], deltax))
                end
                gramscale = 0.0
                for i in 1:m
                    gramscale = max(gramscale, gram[i, i])
                end
                for i in 1:m
                    gram[i, i] += 1e-12*gramscale
                end
                gamma = try
                    view(gram, 1:m, 1:m) \ view(rhs, 1:m)
                catch
                    nothing
                end
                if !isnothing(gamma) && all(isfinite, gamma)
                    @. xcandidate = x + deltax
                    for j in 1:m
                        axpy!(-gamma[j], deltaxhistory[j], xcandidate)
                        axpy!(-gamma[j], deltafhistory[j], xcandidate)
                    end
                    fj!(Fanderson, nothing, xcandidate)
                    normFanderson = norm(Fanderson)
                    if isfinite(normFanderson) && normFanderson < 0.9*sqrt(2*f)
                        andersonupdate = 0.0
                        for i in eachindex(x)
                            andersonupdate += abs2(xcandidate[i] - x[i])
                        end
                        andersonupdate = sqrt(andersonupdate)
                        copyto!(x, xcandidate)
                        copyto!(F, Fanderson)
                        andersonaccepted = true
                        # the size of the accepted extrapolation, used by
                        # the common convergence test below
                        updatenorm = andersonupdate
                        push!(alphas, NaN)
                        push!(backtrackrecord, 0)
                        push!(andersonrecord, true)
                    end
                end
            end
        end

        # update x, verifying sufficient decrease of the objective at the
        # chosen step and backtracking when the quadratic interpolation
        # fails, so that the iteration never takes a step which increases
        # the residual when a smaller acceptable step exists. the
        # directional derivative of the objective along the Newton update
        # is dfdalpha, so the Armijo sufficient decrease condition is
        # falpha <= f + c*alpha*dfdalpha with a small constant c.
        if !andersonaccepted
            if alpha1 == 1.0
                # the residual at the full step is already in F
                falpha = fp
            else
                @. xtrial = x + deltax*alpha1
                fj!(F, nothing, xtrial)
                falpha = real(0.5*dot(F,F))
            end
            backtracks = 0
            while !(isfinite(falpha) && falpha <= f + 1e-4*alpha1*dfdalpha) &&
                    alpha1 > alphamin && backtracks < 30
                alpha1 = max(alpha1/2, alphamin)
                @. xtrial = x + deltax*alpha1
                fj!(F, nothing, xtrial)
                falpha = real(0.5*dot(F,F))
                backtracks += 1
            end
            push!(alphas, alpha1)
            push!(backtrackrecord, backtracks)
            push!(andersonrecord, false)
            @. x += deltax*alpha1
            # the size of the accepted update, used by the common
            # convergence test below
            updatenorm = norm(deltax)*abs(alpha1)
        end
        # every path leaves F evaluated at the updated x
        Ffresh = true

        if norm(F,Inf) <= ftol || ( norm(x) > 0 && norm(F)/norm(x) < ftol && updatenorm <= sqrt(ftol)*norm(x))
            # terminate iterations if the infinity norm is less than ftol, or
            # if the relative norm is less than ftol and the Newton update is
            # small relative to the solution. check that norm(x) is greater
            # than zero to avoid divide by zero errors. the relative norm
            # test alone can pass spuriously when the norm of the solution
            # diverges while the residual stays bounded, so also require that
            # the Newton updates are small (ie we are near the solution).
            converged = true
            break
        end
    end

    return IterationInfo(label, NaN, 0.0, converged, length(alphas),
        normF, alphas, backtrackrecord, andersonrecord)
end
