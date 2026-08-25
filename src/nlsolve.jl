
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

The nine argument constructor leaves `krylov` empty, so the direct solvers
construct an `IterationInfo` without mentioning a field which does not apply
to them.
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
    krylov::Vector
end

function IterationInfo(label, parameter, regularization, converged,
    iterations, normresidual, alpha, backtracks, andersonaccepted)
    return IterationInfo(label, parameter, regularization, converged,
        iterations, normresidual, alpha, backtracks, andersonaccepted, [])
end

"""
    AndersonState(x::AbstractVector, depth::Integer)

Preallocated state for Anderson acceleration of the Newton fixed point
iteration `G(x) = x + deltax`: the difference history, the previous
iterate/update pair it is built from, the assembled correction vector, and the
buffers for the extrapolation least squares problem.

The history lives in fixed n×depth matrices with circular column indexing.
[`andersonhistory!`](@ref) overwrites one column per recorded step and updates
the indices. All access to the history columns is through age-ordered indices,
so results are independent of the physical column layout.
"""
mutable struct AndersonState{T, RT, V, M}
    depth::Int
    historyready::Bool
    histcount::Int
    histpos::Int
    # the vectors and the history follow the iterate onto whatever backend it
    # lives on; the Gram matrix and the coefficients are depth by depth and
    # depth long, and stay on the host, where the back substitution that reads
    # them element by element belongs
    xprev::V
    deltaxprev::V
    correction::V
    deltaxhistory::M
    deltafhistory::M
    agecols::Vector{Int}
    qrQ::M
    gram::Matrix{RT}
    gammabuf::Vector{RT}
end

function AndersonState(x::AbstractVector{T}, depth::Integer) where T
    if depth < 1
        throw(ArgumentError(lazy"`depth` = $(depth) must be positive."))
    end
    n = length(x)
    RT = real(T)
    return AndersonState{T, RT, typeof(similar(x, n)),
        typeof(similar(x, n, depth))}(depth, false, 0, 1,
        copy(x), copy(x), similar(x, n),
        similar(x, n, depth), similar(x, n, depth),
        Vector{Int}(undef, depth),
        similar(x, n, depth),
        Matrix{RT}(undef, depth, depth), Vector{RT}(undef, depth))
end


"""
    Factorization(outofplace,inplace,kwargs)
    Factorization(outofplace,inplace,kwargs,batched)

A structure to hold the factorizations and their keyword arguments.

`batched`, when it is not `nothing`, is a function of a
[`BatchedBlockLayout`](@ref) returning the form of this factorization which
treats the matrix as a uniform batch of diagonal blocks: one symbolic analysis
for the pattern every block shares, then a batched numeric factorization. It is
consulted by [`batchedfactorization`](@ref), which is the only place that knows
whether a particular matrix really is such a batch. A factorization with no
batched form leaves it `nothing` and is used unchanged.

```
"""
struct Factorization
    factorize
    factorize!
    kwargs
    batched
end

Factorization(factorize, factorize!, kwargs) =
    Factorization(factorize, factorize!, kwargs, nothing)

function KLUfactorization(;kwargs...)
    # return Factorization(KLU.klu,KLU.klu!,kwargs)
    return Factorization(KLU.klu,klunzval!,kwargs)
end

# this directly uses the nzvals without checking the sparse matrix so should
# be slightly faster.
function klunzval!(F,A;kwargs...)
    return KLU.klu!(F,A.nzval;kwargs...)
end

"""
    cudssbatchedfactorization(layout::BatchedBlockLayout; kwargs...)

A batched cuDSS [`Factorization`](@ref). Defined by the CUDSS extension; see
[`CUDSSFactorization`](@ref). Without CUDSS.jl and CUDA.jl loaded this raises
an informative error.
"""
function cudssbatchedfactorization(layout; kwargs...)
    throw(ArgumentError(
        "cudssbatchedfactorization requires CUDSS.jl and CUDA.jl to be loaded."))
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
    trysolvetranspose!(x,factorization,b)

Solve the transposed linear system `transpose(A)*x = b` using an existing
factorization of `A`, without refactorizing. The non-conjugating transpose is
used, not the adjoint. Sparse LU factorizations support this directly with a
pair of triangular solves against the stored factors (`klu_tsolve` for KLU), so
an adjoint solve costs a solve rather than a factorization. As in
[`trysolve!`](@ref), fall back to `\\\\` for factorizations which do not support
`ldiv!` with a transposed factorization.

Used by [`hblinsolve`](@ref) to obtain the solutions of the transposed
linearized system, which are the adjoint solutions required by the noise,
quantum efficiency, and sensitivity calculations.
"""
function trysolvetranspose!(x,factorization,b)
    return trysolve!(x,transpose(factorization),b)
end

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
the range `[0, ϕα0]`. The fitting process uses the merit function values at
`α = 0`, `α = ϕα0`, `α = ϕα1`, and the derivative at the first point
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
- `safeguard = 0.1`: the smallest value we allow the step to take relative to
    the full step `α1`. This protects against large (eg. order of magnitude)
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
        Fbest = copy(F), ϕfullstep = nothing)

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
    ϕfullstep::Union{Nothing, Real} = nothing)

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
    α, ϕpred, accepted = quadratic_trial_step(ϕ0, ϕ1, dϕ0dα;
        c1 = c1, safeguard = safeguard_low)
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
        if havefinite
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
            # no finite older trial yet (every longer step overflowed):
            # test Armijo directly and back off geometrically. I could
            # probably use the quadratic trial step function for this.
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


"""
    andersonhistory!(s::AndersonState, x, deltax)

Record one step into the Anderson history. Once a previous iterate/update pair
exists, write the differences `x - xprev` and `deltax - deltaxprev` into the
oldest history column in place, then refresh the previous pair, so each stored
difference spans exactly one step. The first call after construction (or
after a full reset) only establishes the pair.
"""
function andersonhistory!(s::AndersonState, x::AbstractVector,
    deltax::AbstractVector)
    if s.historyready
        @views @. s.deltaxhistory[:, s.histpos] = x - s.xprev
        @views @. s.deltafhistory[:, s.histpos] = deltax - s.deltaxprev
        s.histpos = s.histpos == s.depth ? 1 : s.histpos + 1
        s.histcount = min(s.histcount + 1, s.depth)
    end
    copyto!(s.xprev, x)
    copyto!(s.deltaxprev, deltax)
    s.historyready = true
    return s
end

"""
    andersoncorrection!(s::AndersonState, deltax) -> Bool

Assemble the Type-II Anderson correction `cₖ = (Sₖ + Yₖ)γₖ` into
`s.correction` from the current history and the Newton update `deltax`. Return
`true` when a usable correction was produced and `false` otherwise, such as
when the history is empty or the coefficient solve fails  or produces
non-finite values, in which case `s.correction` must not be used.

The real extrapolation coefficients γ minimize `||ΔF*γ - deltax||` via
ridge-regularized normal equations (with a default `ridge = 1e-12` relative to
the largest Gram diagonal), solved by LU with partial pivoting in place in the
preallocated buffers. All history access is through age-ordered column indices
(oldest first). The coefficients are constrained real because for the harmonic
balance quasi-Newton map the error operator is antilinear (involves complex
conjugation), so complex coefficients cannot cancel the error modes; real
coefficients correspond to Anderson acceleration of the equivalent real
system.
"""
function andersoncorrection!(s::AndersonState{T}, deltax::AbstractVector;
    rtol = eps(real(T))^(3//4)) where T
    s.histcount > 0 || return false
    m = s.histcount
    # age-ordered column indices (oldest first): when the buffer is full
    # the oldest column is the next overwrite slot, and before that
    # columns were written sequentially from 1
    for k in 1:m
        s.agecols[k] = s.histcount < s.depth ? k :
            (s.histpos + k - 1 > s.depth ?
                s.histpos + k - 1 - s.depth : s.histpos + k - 1)
    end
    # Thin QR of the age-ordered history ΔF by modified Gram-Schmidt, in the
    # real inner product real(dot(a, b)), which is the Euclidean inner product
    # of the stacked real and imaginary parts.
    #
    # A column whose norm collapses after orthogonalization is linearly
    # dependent on the older ones and carries no new direction. Rather than
    # regularizing it back to invertibility, truncate the history there and
    # keep the well conditioned leading columns; `rtol` is relative to the
    # largest norm seen, so the test is scale free.
    Q = s.qrQ
    R = s.gram
    rank = 0
    maxnorm = zero(real(T))
    for j in 1:m
        qj = view(Q, :, j)
        copyto!(qj, view(s.deltafhistory, :, s.agecols[j]))
        # two passes of modified Gram-Schmidt: one pass can lose orthogonality
        # on the near-dependent histories this is meant to handle. the
        # coefficients from both passes accumulate into the same R entry
        for i in 1:rank
            R[i, j] = zero(real(T))
        end
        for _ in 1:2
            for i in 1:rank
                qi = view(Q, :, i)
                c = real(dot(qi, qj))
                R[i, j] += c
                axpy!(-c, qi, qj)
            end
        end
        nrm = norm(qj)
        maxnorm = max(maxnorm, nrm)
        if nrm <= rtol*maxnorm || !isfinite(nrm)
            break
        end
        rank += 1
        R[rank, j] = nrm
        rmul!(qj, inv(nrm))
    end
    rank > 0 || return false
    m = rank

    # the least squares right hand side in the orthonormal basis
    for i in 1:m
        s.gammabuf[i] = real(dot(view(Q, :, i), deltax))
    end

    # back substitute R*gamma = Q'*deltax
    gamma = view(s.gammabuf, 1:m)
    for i in m:-1:1
        acc = gamma[i]
        for k in i+1:m
            acc -= R[i, k]*gamma[k]
        end
        if iszero(R[i, i]) || !isfinite(acc)
            return false
        end
        gamma[i] = acc/R[i, i]
    end
    if all(isfinite, gamma)
        fill!(s.correction, zero(T))
        for j in 1:m
            axpy!(gamma[j], view(s.deltaxhistory, :, s.agecols[j]), s.correction)
            axpy!(gamma[j], view(s.deltafhistory, :, s.agecols[j]), s.correction)
        end
        # finite coefficients do not guarantee a finite correction: near-
        # collinear history can overflow the accumulation
        if all(isfinite, s.correction)
            return true
        end
        return false
    end
    return false
end

"""
    andersonrestart!(s::AndersonState)

Discard the entire Anderson history (but keep the preallocated storage).
This can be used after, for example, linesearch failures which suggest the
local Anderson model was not useful. The previous iterate/update pair is kept,
so the acceleration resumes from the next accepted step's difference.
"""
function andersonrestart!(s::AndersonState)
    s.histcount = 0
    s.histpos = 1
    return s
end

"""
    dualsearch!(f!, F, xcandidate, x, deltax, ϕ0, dϕ0dα, ϕcand,
        correction, betak, Fbest, Fspare, Fatx; c1, safeguard_low,
        safeguard_high, maxbacktracks)



Search selection for a rejected Anderson candidate: Both the curvilinear
search `x + α*deltax - betak*α²*correction` and the plain damped-Newton search
are run, and the better accepted point is taken. `F` and `xcandidate` should
hold the candidate's residual and trial point (the curved path's α = 1),
`ϕcand` the merit, and `Fbest` the residual at `x`.

The motivation is an iteration costs a Jacobian evaluation, factorization, and
a linear solve which for a typical device like a TWPA take an order of
magnitude more time than a residual evaluation, so extra residual evaluations
are worth it if it results in a better path (which they appear to).

If both searches fail, the better best-effort point is returned and
`accepted` is false. `Fatx` preserves the residual at `x` across the
searches (which clobber `Fbest` with best-trial copies), keeping every
restore contract valid.

Returns `(α, ϕα, accepted, backtracks, usedcorrection)`, where
`usedcorrection` reports whether the returned point lies on the curved
path, and `backtracks` counts ALL trial evaluations after the candidate
(the curved search's backtracks plus the plain search's full step and
backtracks), so the iteration's residual evaluations are exactly
`1 + backtracks`.
"""
function dualsearch!(f!, F::AbstractVector,
    xcandidate::AbstractVector, x::AbstractVector,
    deltax::AbstractVector, ϕ0::Real, dϕ0dα::Real, ϕcand::Real,
    correction::AbstractVector, betak::Real, Fbest::AbstractVector,
    Fspare::AbstractVector, Fatx::AbstractVector;
    c1 = 1e-4, safeguard_low = 0.1, safeguard_high = 0.5,
    maxbacktracks::Integer = 10, curvedpriority::Bool = false)

    # preserve F(x): the searches overwrite Fbest with best trials
    copyto!(Fatx, Fbest)

    # curved search first (it is the seeded one)
    αc, ϕc, accc, btc = backtracking_linesearch!(
        f!, F, xcandidate, x, deltax, ϕ0, dϕ0dα;
        c1 = c1, safeguard_low = safeguard_low,
        safeguard_high = safeguard_high, maxbacktracks = maxbacktracks,
        correction = correction, beta = betak,
        Fbest = Fbest, ϕfullstep = ϕcand)
    # Accept a curvilinear Anderson step that passed the Armijo test at the
    # full step without comparing to the linear step. If damping (α < 1) is
    # required, that suggests doubt over whether this is the better path, so
    # we will compare with the linear step. If `curvedpriority` (set true if
    # both fail) then any curved acceptance is returned and the linear search
    # is only run if no curved step is accepted.
    if accc && (isone(αc) || curvedpriority)
        return αc, ϕc, true, btc, true
    end
    # store the result and run the plain search
    copyto!(Fspare, F)
    copyto!(Fbest, Fatx)
    αp, ϕp, accp, btp = backtracking_linesearch!(
        f!, F, xcandidate, x, deltax, ϕ0, dϕ0dα;
        c1 = c1, safeguard_low = safeguard_low,
        safeguard_high = safeguard_high, maxbacktracks = maxbacktracks,
        Fbest = Fbest)
    backtracks = btc + 1 + btp

    # Select the point to use. If only curvilinear fulfills the Armijo
    # sufficient-decrease condition then go with that.
    curvedwins = if (accc && !accp)
        true
    # If only linear passes go with that.
    elseif (!accc && accp)
        false
    # Otherwise, compare them, with linear as a fallback if the merit function
    # values are identical.
    else
        if isfinite(ϕc) && !(isfinite(ϕp) && ϕp <= ϕc)
            true
        else
            false
        end
    end
    # curvedwins = (accc && !accp) ||
        # (accc == accp && isfinite(ϕc) && !(isfinite(ϕp) && ϕp <= ϕc))
    if curvedwins
        copyto!(F, Fspare)
        linesearchtrialpoint!(xcandidate, x, αc, deltax, betak, correction)
        return αc, ϕc, accc, backtracks, true
    end
    return αp, ϕp, accp, backtracks, false
end

"""
    solveonbackend!(fj!, F, J, x, backend; kwargs...)

Solve the nonlinear system on the backend its Jacobian and factorization live
on, returning the iteration information and writing the converged state back
into the host vectors `F` and `x`.

The state has to go where the Jacobian is, so that the assembly, the linear
solve and the line search all stay on one side. On `CPU()` `tobackend` adopts
the caller's vectors and the copies back are between an array and itself, so
this is exactly [`nlsolve!`](@ref).
"""
function solveonbackend!(fj!::Function, F::AbstractVector, J,
    x::AbstractVector, backend; kwargs...)

    xb = tobackend(backend, x)
    Fb = tobackend(backend, F)
    info = nlsolve!(fj!, Fb, J, xb; kwargs...)
    copyto!(x, tohost(xb))
    copyto!(F, tohost(Fb))
    return info
end


"""
    nlsolve!(fj!, F, J, x; kwargs...)

Newton's method with Anderson acceleration, designed for quasi-Newton problems
with an approximate Jacobian. `fj!(F, J, x)` must write the residual into `F`
when `F !== nothing` and the Jacobian into `J` when `J != nothing`. `x` is
updated in place. `F` holds the residual at the returned `x` after `nlsolve!`
returns.

Each iteration solves `J*pₖ = -F` for the Newton step and forms the
Type-II Anderson correction `cₖ = (Sₖ + Yₖ)γₖ` from a depth-
`andersondepth` history of iterate and update differences (the Newton
fixed point map `G(x) = x + pₖ` has fixed point residual `pₖ`, which is
stored). The step is then chosen by measurement:

1. The Anderson candidate `x + pₖ - andersonbeta*cₖ` is evaluated, and
accepted if its residual norm improves on the current one by at least
`andersonacceptfactor`. An accepted candidate records `alpha = NaN`.

2. Otherwise both line searches are run and compared. The curvilinear path
`x + α*pₖ - andersonbeta*α²*cₖ`, starting with the value from (1), and the 
linear path and the point which results in the lowest merit function is used.
A curvilinear trial point that fulfills the Armijo condition at the full step
is accepted without comparison. If both line searches fail, the solver sets 
curved priority and any curved acceptance is taken with the linear path only
taken if the curved path results in non-finite values.See `dualsearch!`.
3. If there is no usable correction (empty history, failed coefficient
solve), only the linear path line search runs.

# Keywords
- `iterations = 1000`: maximum number of Newton iterations.
- `ftol = 1e-8`: convergence when `norm(F) <= ftol`.
- `factorization = KLUfactorization()`: factorization method.
- `label = ""`: label for the returned `IterationInfo`.
- `c1 = 1e-4`: Armijo sufficient-decrease constant, in (0, 1/2); the
  upper bound keeps the full Newton step acceptable near a root.
- `safeguard_low = 0.1`, `safeguard_high = 0.5`: backtracking step
  clamp as fractions of the previous trial.
- `maxbacktracks = 10`: trial-point budget per line search.
- `maxbacktrackfailures = 2`: consecutive-failure stall threshold.
- `andersondepth = 5`: history depth; `0` disables the acceleration.
  Hard problems may need more depth (5 was required on the hardest
  examples).
- `andersonbeta = 1.0`: correction strength.
- `andersonacceptfactor = 0.9`: candidate accept threshold, in (0, 1);
  smaller is stricter.

Returns an `IterationInfo` with per-iteration diagnostics: `normF`
(residual norms, starting at the initial point), `alphas` (step
lengths; `NaN` marks an accepted candidate), `backtrackrecord` (trial
evaluations after each iteration's first), and `andersonrecord` (true
when the taken step lies on the curved path).
"""
function nlsolve!(fj!::Function, F::AbstractVector{T}, J::AbstractArray{T},
    x::AbstractVector{T}; iterations = 1000, ftol = 1e-8,
    factorization = KLUfactorization(), label = "",
    c1 = 1e-4, safeguard_low = 0.1, safeguard_high = 0.5,
    maxbacktracks::Integer = 10, maxbacktrackfailures::Integer = 2,
    andersondepth::Integer = 5, andersonbeta = 1.0,
    andersonacceptfactor = 0.9) where T

    if size(J, 1) != size(J, 2)
        throw(DimensionMismatch(lazy"The Jacobian `J` matrix must be square."))
    end

    if size(J, 2) != length(x)
        throw(DimensionMismatch(lazy"Second axis of Jacobian `J` must have the same length as the input `x`."))
    end

    if size(J, 1) != length(F)
        throw(DimensionMismatch(lazy"First axis of the Jacobian `J` must have the same length as the residual `F`."))
    end

    cache = FactorizationCache()

    deltax = copy(x)
    xcandidate = copy(x)
    Fbest = similar(F)
    Fspare = similar(F)
    Fatx = similar(F)
    xinitial = copy(x)

    # state for Anderson acceleration of the Newton fixed point map
    # G(x) = x + deltax, whose fixed point residual is the Newton update
    # deltax itself, so the history of iterates and updates is available
    # at no extra cost. see AndersonState, andersonhistory!,
    # andersoncorrection!, and andersonrestart!
    anderson = andersondepth > 0 ? AndersonState(x, andersondepth) : nothing

    ### diagnostic info
    backtrackrecord = Int[]
    andersonrecord = Bool[]
    alphas = real(T)[]
    normF = real(T)[]
    converged = false
    backtrackfailures = 0

    # validate every option before the first residual evaluation
    if iterations < 0
        throw(ArgumentError(lazy"`iterations` = $(iterations) must be nonnegative."))
    end
    if !(ftol >= 0)
        throw(ArgumentError(lazy"`ftol` = $(ftol) must be nonnegative."))
    end
    # for the exact Newton step dϕ0 = -2ϕ0, so the full-step Armijo bound
    # is (1 - 2c1)ϕ0: any c1 >= 1/2 makes full-step acceptance impossible
    # for a nonnegative merit, destroying the fast local convergence
    if !(0 < c1 < 1//2)
        throw(ArgumentError(lazy"`c1` = $(c1) must be in (0, 1/2) for the Newton merit function."))
    end
    if !(0 < safeguard_low < 1//2)
        throw(ArgumentError(lazy"`safeguard_low` = $(safeguard_low) must be in (0, 1/2)."))
    end
    if !(safeguard_low < safeguard_high < 1)
        throw(ArgumentError(lazy"`safeguard_high` = $(safeguard_high) must satisfy `safeguard_low < safeguard_high < 1`."))
    end
    if maxbacktracks < 0
        throw(ArgumentError(lazy"`maxbacktracks` = $(maxbacktracks) must be nonnegative."))
    end
    if maxbacktrackfailures < 1
        throw(ArgumentError(lazy"`maxbacktrackfailures` = $(maxbacktrackfailures) must be positive."))
    end
    if andersondepth < 0
        throw(ArgumentError(lazy"`andersondepth` = $(andersondepth) must be nonnegative."))
    end
    if !(0 <= andersonbeta) || !isfinite(andersonbeta)
        throw(ArgumentError(lazy"`andersonbeta` = $(andersonbeta) must be finite and nonnegative."))
    end
    if !(0 < andersonacceptfactor < 1)
        throw(ArgumentError(lazy"`andersonacceptfactor` = $(andersonacceptfactor) must be in (0, 1)."))
    end

    # residual-only adapter for the linesearch, which never needs the
    # Jacobian and therefore does not accept the combined fj! interface
    residual!(Fv, xv) = fj!(Fv, nothing, xv)


    # run the fast comparison first. on a stall (the failure-counter or
    # no-decrease exit) restart once from the initial point with curved
    # priority on a fresh trajectory.
    curvedpriority = false
    stalled = false

    for attempt in 1:2
        if attempt == 2
            # retry by resetting to initial values and setting curved priority
            # the motivation for this is some problems may require the
            # corrections from Anderson acceleration to converge.
            # @warn string(lazy"Second attempt: restarting with curved priority.")
            copyto!(x, xinitial)
            empty!(normF); empty!(alphas)
            empty!(backtrackrecord); empty!(andersonrecord)
            if !isnothing(anderson)
                andersonrestart!(anderson)
                anderson.historyready = false
            end
            curvedpriority = true
            stalled = false
            converged = false
            backtrackfailures = 0
        end

        # evaluate the residual at the (initial or restarted) point and decide
        # convergence before touching the Jacobian. the residual norm at the
        # initial point; every later entry of normF is pushed immediately
        # after a step is accepted, so convergence is decided on each fresh
        # residual before the Jacobian is refreshed and no Jacobian is ever
        # evaluated at a final point
        residual!(F, x)
        push!(normF, norm(F))
        if normF[end] <= ftol
            converged = true
        else
            # only a point from which a step will be taken needs a Jacobian.
            fj!(nothing, J, x)
            tryfactorize!(cache, factorization, J)
        end

        # perform Newton's method with linesearch based on Nocedal and Wright
        # chapter 3 section 5.
        for n in 1:iterations
            converged && break

            # F and x are consistent here, and cache.factorization matches
            # the J from which deltax will be computed

            # solve the linear system
            trysolve!(deltax, cache.factorization, F)

            # multiply deltax by -1
            rmul!(deltax, -1)

            # record this step into the Anderson history. Every actual step
            # is stored (the history samples the true trajectory); whether
            # a correction is used is decided per iteration by measurement.
            if !isnothing(anderson)
                andersonhistory!(anderson, x, deltax)
            end

            # calculate the objective function and the derivative of the
            # objective with respect to the scalar variable alpha which
            # parameterizes the path between the old x and the new x.
            # Note: the dot product takes the complex conjugate of the first
            # vector
            ϕ0 = real(0.5*dot(F, F))
            # the model slope Re(F'*J*deltax), exactly -||F||² for the J
            # that produced deltax; for a quasi-Newton J this is the model's
            # claim, not the true directional derivative
            dϕ0dα = real(dot(F, J, deltax))

            # check before the Armijo tests below, which would otherwise run
            # with an invalid slope (the trial-step helpers validate too, but
            # only on their paths).
            #
            # a non-finite merit, or a search direction which is not a descent
            # direction, is a numerical outcome of this solve rather than a
            # caller error: `deltax` comes from a factorization of `J`, and for
            # `method = :quasinewton` that `J` is only an approximation, so it
            # can propose a direction the true merit does not decrease along.
            # treat it as a stall, which stops this attempt and lets the robust
            # retry below run on a fresh trajectory, instead of throwing out of
            # the middle of the solve and discarding a usable best point.
            if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
                stalled = true
                break
            end

            # whether this iteration escalated the search policy (set below
            # on a both-fail event); the escalating failure is charged to
            # the old policy in the failure counting at the bottom
            escalated = false

            # Anderson candidate x + pₖ - βcₖ: evaluate its residual
            # and accept it outright if the residual norm improves by at
            # least `andersonacceptfactor`. an accepted candidate
            # records alpha = NaN and andersonrecord = true.
            candidateaccepted = false
            havecorrection = !isnothing(anderson) &&
                andersoncorrection!(anderson, deltax)
            if havecorrection
                betak = andersonbeta
                # save F(x) before the candidate evaluation so the line
                # search restore paths stay valid on rejection
                copyto!(Fbest, F)
                ϕcand = linesearchevaluate!(residual!, F, xcandidate, x,
                    one(ϕ0), deltax, betak, anderson.correction)
                if isfinite(ϕcand) &&
                        sqrt(2*ϕcand) < andersonacceptfactor*normF[end]
                    alpha1 = real(T)(NaN)
                    ϕα, accepted, backtracks = ϕcand, true, 0
                    candidateaccepted = true
                end
            end
            if candidateaccepted
                usecorrection = true
            elseif havecorrection
                # race both searches and take the measured better point:
                # see dualsearch! for the contract and rationale.
                # usecorrection reports whether the taken step lies on
                # the curved path.
                alpha1, ϕα, accepted, backtracks, usecorrection =
                    dualsearch!(residual!, F, xcandidate, x,
                        deltax, ϕ0, dϕ0dα, ϕcand, anderson.correction,
                        betak, Fbest, Fspare, Fatx;
                        c1 = c1, safeguard_low = safeguard_low,
                        safeguard_high = safeguard_high,
                        maxbacktracks = maxbacktracks,
                        curvedpriority = curvedpriority)
                if !accepted && !curvedpriority
                    # both searches failed: permanently set to curved priority,
                    # the option observed to converge when this occurs, and
                    # reset backtrackfailures
                    curvedpriority = true
                    # after the backtrackfailures+=1 below this will be set
                    # to 0.
                    backtrackfailures = -1
                end
            else
                # no usable correction: line search along pₖ
                usecorrection = false
                alpha1, ϕα, accepted, backtracks = backtracking_linesearch!(
                    residual!, F, xcandidate, x, deltax, ϕ0, dϕ0dα;
                    c1 = c1, safeguard_low = safeguard_low,
                    safeguard_high = safeguard_high,
                    maxbacktracks = maxbacktracks, Fbest = Fbest)
            end
            push!(andersonrecord, usecorrection)

            push!(alphas, alpha1)
            push!(backtrackrecord, backtracks)

            if iszero(alpha1)
                # the linesearch exhausted its budget without finding any
                # decrease along the Newton direction; F again holds the
                # residual at the unchanged x, so no progress is possible and
                # we stop, reporting non-convergence (this is a stall, so we
                # will attempt the robust retry).
                stalled = true
                break
            end

            # accept the trial point: F already holds the residual there (the
            # linesearch postcondition), and convergence is decided on it now,
            # before any Jacobian evaluation.
            copyto!(x, xcandidate)
            push!(normF, norm(F))
            if normF[end] <= ftol
                converged = true
                break
            end

            # count consecutive linesearch failures where the best decreasing
            # trial step is returned instead of an Armijo-accepted step.
            # an accepted step resets the count, so an isolated failure in an
            # otherwise recovering solve is ignored. repeated failures suggest
            # a stall. the iterations cannot achieve sufficient decrease,
            # so we will try the robust retry then give up, reporting the
            # failure to converge.
            if accepted
                backtrackfailures = 0
            elseif escalated
                # the escalating failure is charged to the policy it ended;
                # the curved-priority policy starts with a clean budget
                backtrackfailures = 0
            else
                backtrackfailures += 1
                if backtrackfailures >= maxbacktrackfailures
                    # the (decreasing) step was already taken above, so the
                    # reported state is the best point found; stop without
                    # refreshing the Jacobian, which would not be used (a
                    # stall will get one robust retry).
                    stalled = true
                    break
                end
            end

            # refresh and refactor the Jacobian at the new x; convergence was
            # already decided above, so no Jacobian is evaluated at the final
            # point
            fj!(nothing, J, x)
            tryfactorize!(cache, factorization, J)
        end

        # retry only on a stall (not on budget exhaustion), only once, and
        # only if the fast attempt actually used the acceleration (otherwise
        # the robust policy could not differ)
        if converged || !stalled || isnothing(anderson) || !any(andersonrecord)
            break
        end
    end

    return IterationInfo(label, NaN, 0.0, converged, length(alphas),
        normF, alphas, backtrackrecord, andersonrecord)
end

