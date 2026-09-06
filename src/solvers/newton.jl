# Newton and quasi-Newton with an assembled Jacobian and a sparse
# factorization, with Anderson acceleration of the fixed point iteration.

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

# `factorize` of the host sparse factorizations; the types are defined in
# solvers/options.jl. KLU's fill reducing ordering is chosen by measurement
# (`kluordered`) and its refactorization reuses the symbolic analysis.
factorize(f::KLUfactorization, A) = kluordered(A; f.kwargs...)
refactorize!(f::KLUfactorization, F, A) = klunzval!(F, A; f.kwargs...)
factorize(f::LUfactorization, A) = lu(A; f.kwargs...)
refactorize!(f::LUfactorization, F, A) = lu!(F, A; f.kwargs...)
factorize(f::QRfactorization, A) = qr(A; f.kwargs...)

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
    andersoncorrection!(s::AndersonState, deltax; rtol = eps(T)^(3//4)) -> Bool

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
        safeguard_high, maxbacktracks, curvedpriority = false)

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
    nlsolve!(fj!, F, J, x; iterations = 1000, ftol = 1e-8, rtol = 0.0,
        factorization = KLUfactorization(), label = "", c1 = 1e-4,
        safeguard_low = 0.1, safeguard_high = 0.5, maxbacktracks = 10,
        maxbacktrackfailures = 2, andersondepth = 5, andersonbeta = 1.0,
        andersonacceptfactor = 0.9)

Newton's method with a line search and Anderson acceleration, suited to
quasi-Newton problems with an approximate Jacobian. `fj!(F, J, x)` must
write the residual into `F` when `F !== nothing` and the Jacobian into `J`
when `J !== nothing`. `x` is updated in place and holds the solution on
return; `F` holds the residual there.

Each iteration solves `J*pₖ = -F` for the Newton step and forms the
Type-II Anderson correction `cₖ = (Sₖ + Yₖ)γₖ` from a depth-
`andersondepth` history of iterate and update differences (the Newton
fixed point map `G(x) = x + pₖ` has fixed point residual `pₖ`, which is
stored). The step is then chosen by measurement:

1. The Anderson candidate `x + pₖ - andersonbeta*cₖ` is evaluated, and
accepted if its residual norm improves on the current one by at least
`andersonacceptfactor`. An accepted candidate records `alpha = NaN`.

2. Otherwise both line searches are run and compared: the curvilinear path
`x + α*pₖ - andersonbeta*α²*cₖ`, starting from the value of (1), and the
linear path `x + α*pₖ`; the point with the lower merit function is taken.
A curvilinear trial point which satisfies the Armijo condition at the full
step is accepted without comparison. After both line searches have failed
the solver gives the curved path priority, taking the linear path only
when the curved one produces non-finite values. See [`dualsearch!`](@ref).

3. Without a usable correction (an empty history, or a failed coefficient
solve) only the linear line search runs.

# Keywords
- `iterations = 1000`: the maximum number of Newton iterations.
- `ftol = 1e-8`: converged when `norm(F) <= ftol`.
- `rtol = 0.0`: a relative tolerance; the effective tolerance is
    `max(ftol, rtol*norm(F0))` with `F0` the initial residual.
- `factorization = KLUfactorization()`: the sparse factorization of `J`.
- `label = ""`: label for the returned `IterationInfo`.
- `c1 = 1e-4`: Armijo sufficient-decrease constant, in (0, 1/2); the
  upper bound keeps the full Newton step acceptable near a root.
- `safeguard_low = 0.1`, `safeguard_high = 0.5`: backtracking step
  clamp as fractions of the previous trial.
- `maxbacktracks = 10`: trial-point budget per line search.
- `maxbacktrackfailures = 2`: consecutive-failure stall threshold.
- `andersondepth = 5`: the history depth; `0` disables the acceleration.
- `andersonbeta = 1.0`: correction strength.
- `andersonacceptfactor = 0.9`: candidate accept threshold, in (0, 1);
  smaller is stricter.

Returns an [`IterationInfo`](@ref) with per-iteration diagnostics:
`normresidual` (residual norms, starting at the initial point), `alpha`
(step lengths; `NaN` marks an accepted candidate), `backtracks` (trial
evaluations after each iteration's first), `andersonaccepted` (true when
the taken step lies on the curved path), and `reason`, why the iteration
ended: `:converged`, `:iterations`, `:linesearch` (no decrease at all, or
two consecutive steps short of the Armijo condition), or `:progress` (the
residual history projects no convergence within the remaining budget,
[`projectedstall`](@ref)); see [`stallmessage`](@ref).
"""
function nlsolve!(fj!::Function, F::AbstractVector{T}, J::AbstractArray{T},
    x::AbstractVector{T}; iterations = 1000, ftol = 1e-8, rtol = 0.0,
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

    # the record and the acceptance rules shared with `nlsolvekrylov!`
    tr = NewtonTrace{real(T)}(maxbacktrackfailures)
    normF = tr.normresidual

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
            if !isnothing(anderson)
                andersonrestart!(anderson)
                anderson.historyready = false
            end
            curvedpriority = true
            stalled = false
        end

        # evaluate the residual at the (initial or restarted) point and decide
        # convergence before touching the Jacobian. every later entry of the
        # history is pushed immediately after a step is accepted, so
        # convergence is decided on each fresh residual before the Jacobian
        # is refreshed and no Jacobian is ever evaluated at a final point.
        # `rtol` adds the relative test beside the absolute one, satisfied
        # when either holds; at `rtol = 0` the tolerance is exactly `ftol`
        # and nothing already measured moves. See `nlsolvekrylov!`.
        residual!(F, x)
        if !tracestart!(tr, F, ftol, rtol)
            # only a point from which a step will be taken needs a Jacobian.
            fj!(nothing, J, x)
            tryfactorize!(cache, factorization, J)
        end

        # perform Newton's method with linesearch based on Nocedal and Wright
        # chapter 3 section 5.
        for n in 1:iterations
            tr.converged && break

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
                tr.reason = :linesearch
                stalled = true
                break
            end

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
                    # after the count below this will be set to 0.
                    tr.backtrackfailures = -1
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
            tracetrial!(tr, alpha1, backtracks, usecorrection)

            if iszero(alpha1)
                # the linesearch exhausted its budget without finding any
                # decrease along the Newton direction; F again holds the
                # residual at the unchanged x, so no progress is possible and
                # we stop, reporting non-convergence (this is a stall, so we
                # will attempt the robust retry).
                tr.reason = :linesearch
                stalled = true
                break
            end

            # accept the trial point: F already holds the residual there (the
            # linesearch postcondition), and convergence is decided on it now,
            # before any Jacobian evaluation. Repeated line searches which
            # return the best decreasing trial instead of an Armijo accepted
            # step are a stall: the iterations cannot achieve sufficient
            # decrease, so we will try the robust retry then give up. The
            # (decreasing) step was already taken, so the reported state is
            # the best point found, and the Jacobian is not refreshed since
            # it would not be used.
            copyto!(x, xcandidate)
            outcome = tracestep!(tr, F, accepted)
            outcome === :converged && break
            if outcome === :linesearch
                stalled = true
                break
            end

            if tracestalled(tr, 1, iterations - n)
                stalled = true
                tr.reason = :progress
                break
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
        if tr.converged || !stalled || isnothing(anderson) ||
                !any(tr.andersonaccepted)
            break
        end
    end

    return IterationInfo(tr, label)
end
