
"""
    GMRESWorkspace{T<:AbstractFloat}

Preallocated storage for [`gmres!`](@ref) with a restart length of `m` on a
system of dimension `n`. Holds the `n x (m+1)` Arnoldi basis `V`, the
`(m+1) x m` Hessenberg matrix `H`, the Givens rotations `cs` and `sn` which
reduce it, the least squares right hand side `s`, its solution `y`, and three
length `n` work vectors.

The dominant cost is `V`, which is `n*(m+1)` numbers, so `m` trades memory and
orthogonalization work against restart frequency.
"""
struct GMRESWorkspace{T<:AbstractFloat}
    V::Matrix{T}
    H::Matrix{T}
    cs::Vector{T}
    sn::Vector{T}
    s::Vector{T}
    y::Vector{T}
    w::Vector{T}
    z::Vector{T}
    u::Vector{T}
end

function GMRESWorkspace(n::Integer, m::Integer, ::Type{T} = Float64) where {T<:AbstractFloat}
    n >= 0 || throw(ArgumentError(lazy"`n` = $(n) must be nonnegative."))
    m >= 1 || throw(ArgumentError(lazy"the restart length `m` = $(m) must be at least 1."))
    return GMRESWorkspace{T}(
        Matrix{T}(undef, n, m + 1), zeros(T, m + 1, m),
        Vector{T}(undef, m), Vector{T}(undef, m),
        Vector{T}(undef, m + 1), Vector{T}(undef, m),
        Vector{T}(undef, n), Vector{T}(undef, n), Vector{T}(undef, n))
end

"""
    gmres_orthogonalize!(w, V, H, j)

Orthogonalize `w` against the first `j` Arnoldi basis vectors (the columns of
`V`) by modified Gram-Schmidt, accumulating the coefficients into column `j`
of the Hessenberg matrix `H`. A second pass runs when the first one removes
more than `1/sqrt(2)` of the norm of `w` (the DGKS criterion), which is when
a single pass is known to be able to lose orthogonality to roundoff; the
coefficients of both passes accumulate into the same entries of `H`, so `H`
remains the exact projection. Writes the subdiagonal `H[j+1, j]` and returns
`(hsub, normw0)`: the norm of the orthogonalized `w` and its norm on entry,
whose ratio is the breakdown test performed by the caller. Allocation free.
"""
function gmres_orthogonalize!(w::AbstractVector{T}, V::AbstractMatrix{T},
    H::AbstractMatrix{T}, j::Integer) where {T<:AbstractFloat}
    normw0 = norm(w)
    for i in 1:j
        vi = view(V, :, i)
        c = dot(vi, w)
        H[i, j] = c
        axpy!(-c, vi, w)
    end
    hsub = norm(w)
    if hsub < normw0/sqrt(T(2))
        for i in 1:j
            vi = view(V, :, i)
            c = dot(vi, w)
            H[i, j] += c
            axpy!(-c, vi, w)
        end
        hsub = norm(w)
    end
    H[j+1, j] = hsub
    return hsub, normw0
end

"""
    gmres_givens(a, b)

The Givens rotation `(c, s, r)` with `c*a + s*b = r` and `-s*a + c*b = 0`,
computed through `hypot` so it cannot overflow, with the identity rotation
returned for the zero input.
"""
function gmres_givens(a::T, b::T) where {T<:AbstractFloat}
    r = hypot(a, b)
    iszero(r) && return one(T), zero(T), zero(T)
    return a/r, b/r, r
end

"""
    gmres_applyrotations!(H, cs, sn, s, j)

Reduce column `j` of the Hessenberg matrix `H` to upper triangular form:
apply the `j-1` previous Givens rotations to the new column, compute and
store the rotation which annihilates the new subdiagonal `H[j+1, j]`, and
apply it to the least squares right hand side `s`. After this the magnitude
of `s[j+1]` is the residual norm of the least squares problem, which with
right preconditioning is the true residual norm of the original system.
Allocation free.
"""
function gmres_applyrotations!(H::AbstractMatrix{T}, cs::AbstractVector{T},
    sn::AbstractVector{T}, s::AbstractVector{T}, j::Integer) where {T<:AbstractFloat}
    for i in 1:j-1
        tmp       =  cs[i]*H[i, j] + sn[i]*H[i+1, j]
        H[i+1, j] = -sn[i]*H[i, j] + cs[i]*H[i+1, j]
        H[i, j]   = tmp
    end
    cs[j], sn[j], r = gmres_givens(H[j, j], H[j+1, j])
    H[j, j]   = r
    H[j+1, j] = zero(T)
    s[j+1] = -sn[j]*s[j]
    s[j]   =  cs[j]*s[j]
    return abs(s[j+1])
end

"""
    gmres_correction!(x, ws::GMRESWorkspace, j, Mop!)

Solve the reduced `j x j` triangular least squares problem by back
substitution, assemble the correction `u = V[:, 1:j]*y` in the Krylov basis,
undo the right preconditioning once with `Mop!` (or not at all when
`Mop! === nothing`), and add the result to `x` in place. A zero diagonal
entry, which can only arise from an exact breakdown, contributes a zero
coefficient rather than a division by zero. Allocation free.
"""
function gmres_correction!(x::AbstractVector{T}, ws::GMRESWorkspace{T},
    j::Integer, Mop!) where {T<:AbstractFloat}
    H, s, y, V, u, z = ws.H, ws.s, ws.y, ws.V, ws.u, ws.z
    for i in j:-1:1
        acc = s[i]
        for k in i+1:j
            acc -= H[i, k]*y[k]
        end
        y[i] = iszero(H[i, i]) ? zero(T) : acc/H[i, i]
    end
    fill!(u, zero(T))
    for i in 1:j
        axpy!(y[i], view(V, :, i), u)
    end
    if isnothing(Mop!)
        @. x += u
    else
        Mop!(z, u)
        @. x += z
    end
    return x
end

"""
    gmres!(x, Aop!, b, ws::GMRESWorkspace; Mop! = nothing, rtol = 1e-6,
        atol = 0.0, maxrestarts = 10, initialzero = true)

Solve `A*x = b` with restarted GMRES, where `Aop!(w, v)` computes `w = A*v` and
the optional `Mop!(z, v)` applies a preconditioner `z = M \\ v`. The matrix `A`
is never formed; only its action is required, which is what makes this usable
with the matrix-free [`jacobianvectorproduct!`](@ref).

Preconditioning is applied on the right, solving `A*inv(M)*u = b` and then
`x = inv(M)*u`. Right preconditioning keeps the recurrence's residual estimate
equal to the true residual of the original system, so the stopping test is on
`norm(b - A*x)` and does not depend on the quality of `M`. Because `M` is held
fixed across a solve, the preconditioner is applied once per Arnoldi step and
once more per restart, rather than being stored for every basis vector.

The Arnoldi basis is built by modified Gram-Schmidt with a conditional second
pass ([`gmres_orthogonalize!`](@ref)). A subdiagonal which collapses relative
to the vector it came from is a (lucky) breakdown: the Krylov space is
invariant, the reduced least squares solution is exact, and the cycle ends
there rather than continuing with a spurious basis vector. The residual is
recomputed explicitly at every restart so restarts cannot drift from the
recurrence estimate.

Converges when `norm(b - A*x) <= max(rtol*norm(b), atol)`. Returns the named
tuple `(iterations, residual, converged)`, where `iterations` counts Arnoldi
steps (and therefore calls to `Aop!`) across all restarts.

Allocation free after the workspace is built, apart from whatever `Aop!` and
`Mop!` themselves allocate.
"""
function gmres!(x::AbstractVector{T}, Aop!, b::AbstractVector{T},
    ws::GMRESWorkspace{T}; Mop! = nothing, rtol = 1e-6, atol = 0.0,
    maxrestarts::Integer = 10, initialzero::Bool = true) where {T<:AbstractFloat}

    n = length(b)
    length(x) == n || throw(DimensionMismatch(
        lazy"`x` has length $(length(x)) but `b` has length $(n)."))
    size(ws.V, 1) == n || throw(DimensionMismatch(
        lazy"the workspace is for dimension $(size(ws.V,1)) but `b` has length $(n)."))
    rtol >= 0 || throw(ArgumentError(lazy"`rtol` = $(rtol) must be nonnegative."))
    atol >= 0 || throw(ArgumentError(lazy"`atol` = $(atol) must be nonnegative."))
    maxrestarts >= 1 || throw(ArgumentError(
        lazy"`maxrestarts` = $(maxrestarts) must be at least 1."))

    m = size(ws.H, 2)
    V, H, cs, sn, s = ws.V, ws.H, ws.cs, ws.sn, ws.s
    w, z = ws.w, ws.z

    bnorm = norm(b)
    tol = max(rtol*bnorm, atol)

    # a zero right hand side has the zero solution; return it rather than
    # dividing by a zero residual norm below
    if iszero(bnorm)
        fill!(x, zero(T))
        return (iterations = 0, residual = zero(T), converged = true)
    end

    # initial residual w = b - A*x
    if initialzero
        fill!(x, zero(T))
        copyto!(w, b)
    else
        Aop!(w, x)
        @. w = b - w
    end
    resnorm = norm(w)

    totaliterations = 0
    for _ in 1:maxrestarts
        resnorm <= tol && break

        beta = resnorm
        @views V[:, 1] .= w ./ beta
        fill!(s, zero(T))
        s[1] = beta

        j = 0
        while j < m
            j += 1
            # the Arnoldi step on A*inv(M)
            if isnothing(Mop!)
                @views copyto!(z, V[:, j])
            else
                @views Mop!(z, V[:, j])
            end
            Aop!(w, z)

            hsub, normw0 = gmres_orthogonalize!(w, V, H, j)
            resnorm = gmres_applyrotations!(H, cs, sn, s, j)
            totaliterations += 1

            # a subdiagonal which collapsed relative to the incoming vector
            # is a breakdown: the Krylov space is invariant, the reduced
            # solution below is exact, and there is no valid next basis
            # vector to normalize
            breakdown = hsub <= eps(T)*normw0
            (breakdown || resnorm <= tol) && break

            @views V[:, j+1] .= w ./ hsub
        end

        gmres_correction!(x, ws, j, Mop!)

        # recompute the residual explicitly for the next cycle so restarts
        # cannot drift from the recurrence estimate
        Aop!(w, x)
        @. w = b - w
        resnorm = norm(w)
    end

    return (iterations = totaliterations, residual = resnorm,
        converged = resnorm <= tol)
end

"""
    nlsolvekrylov!(fj!, jvp!, F, J, x; kwargs...)

Inexact (Newton-Krylov) solver for a real system: the Newton step is taken
from [`gmres!`](@ref) on the exact matrix-free product `jvp!(y, v)` rather
than from a factorization of an assembled Jacobian. `fj!(F, J, x)` evaluates
the residual and Jacobian as in [`nlsolve!`](@ref), accepting `nothing` for
either. `x` is updated in place and `F` holds the residual at the returned
`x`.

The assembled Jacobian `J` serves only as the (right) preconditioner. It is
reassembled and refactored lazily: when the GMRES iteration count says it has
gone stale (`krylovrefreshiterations`), when GMRES fails, when a step is not
a descent direction, or when a linesearch cannot find any decrease. The
linear tolerance follows the Eisenstat-Walker choice 2 forcing sequence
`krylovgamma*(|F_k|/|F_{k-1}|)^krylovalpha` clamped to
`[krylovrtolmin, krylovrtolmax]`, with an absolute floor of `ftol/10` so late
solves are not pushed below the nonlinear tolerance. Because the assembled
Jacobian can be stale, the linesearch slope is always taken from an exact
matrix-free product, and a non-descent direction falls back to the exact
Newton step through a fresh factorization before the iteration is declared
stalled.

The globalization is the plain damped-Newton path of [`nlsolve!`](@ref):
the interpolated [`backtracking_linesearch!`](@ref), which on Armijo failure
still takes the best decreasing trial, with consecutive failures counted
against `maxbacktrackfailures` and a no-decrease step retried once from a
fresh preconditioner before stopping. There is deliberately no Anderson
acceleration here: the Krylov steps are near-exact Newton steps, and the
solver is kept simple.

# Keywords
- `iterations = 1000`: maximum number of Newton iterations.
- `ftol = 1e-8`: convergence when `norm(F) <= ftol`.
- `factorization = KLUfactorization()`: preconditioner factorization method.
- `label = ""`: label for the returned `IterationInfo`.
- `c1 = 1e-4`: Armijo sufficient-decrease constant, in (0, 1/2).
- `safeguard_low = 0.1`, `safeguard_high = 0.5`: backtracking step clamp as
  fractions of the previous trial.
- `maxbacktracks = 10`: trial-point budget per line search.
- `maxbacktrackfailures = 2`: consecutive-failure stall threshold.
- `krylovrestart = 30`: GMRES restart length.
- `krylovmaxrestarts = 4`: GMRES restart budget per solve.
- `krylovrefreshiterations = 10`: GMRES iteration count above which the
  preconditioner is considered stale and refreshed before the next step.
- `krylovrtolmin = 1e-10`, `krylovrtolmax = 1e-2`: clamp of the forcing
  sequence.
- `krylovgamma = 0.9`, `krylovalpha = (1 + sqrt(5))/2`: Eisenstat-Walker
  forcing parameters.

Returns an [`IterationInfo`](@ref) with the same per-iteration diagnostics
as [`nlsolve!`](@ref); the `andersonaccepted` record is always false.
"""
function nlsolvekrylov!(fj!::Function, jvp!::Function, F::AbstractVector{T},
    J::AbstractArray{T}, x::Vector{T}; iterations = 1000, ftol = 1e-8,
    factorization = KLUfactorization(), label = "",
    c1 = 1e-4, safeguard_low = 0.1, safeguard_high = 0.5,
    maxbacktracks::Integer = 10, maxbacktrackfailures::Integer = 2,
    krylovrestart::Integer = 30, krylovmaxrestarts::Integer = 4,
    krylovrefreshiterations::Integer = 10, krylovrtolmin = 1e-10,
    krylovrtolmax = 1e-2, krylovgamma = 0.9,
    krylovalpha = (1 + sqrt(5))/2) where {T<:AbstractFloat}

    size(J, 1) == size(J, 2) || throw(DimensionMismatch(
        lazy"The Jacobian `J` matrix must be square."))
    size(J, 2) == length(x) || throw(DimensionMismatch(
        lazy"Second axis of Jacobian `J` must have the same length as the input `x`."))
    size(J, 1) == length(F) || throw(DimensionMismatch(
        lazy"First axis of the Jacobian `J` must have the same length as the residual `F`."))

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
    krylovrefreshiterations >= 1 || throw(ArgumentError(
        lazy"`krylovrefreshiterations` = $(krylovrefreshiterations) must be at least 1."))
    0 < krylovrtolmin <= krylovrtolmax < 1 || throw(ArgumentError(
        lazy"`krylovrtolmin` = $(krylovrtolmin) and `krylovrtolmax` = $(krylovrtolmax) must satisfy 0 < krylovrtolmin <= krylovrtolmax < 1."))

    cache = FactorizationCache()
    ws = GMRESWorkspace(length(x), min(krylovrestart, length(x)), T)
    deltax = similar(x)
    xcandidate = similar(x)
    Jv = similar(x)
    Fbest = similar(F)

    # absolute floor for the linear solves: once the linear residual is below
    # the nonlinear tolerance, further accuracy cannot help the Newton
    # iteration, and demanding it makes late GMRES solves "fail"
    gmresatol = real(T)(ftol)/10

    ### diagnostic info
    backtrackrecord = Int[]
    alphas = real(T)[]
    normF = real(T)[]
    converged = false
    backtrackfailures = 0
    refresh = true
    refreshedforstall = false

    Mop!(zv, vv) = trysolve!(zv, cache.factorization, vv)
    # residual-only adapter for the linesearch, which never needs the
    # Jacobian and therefore does not accept the combined fj! interface
    residual!(Fv, xv) = fj!(Fv, nothing, xv)

    # refresh and refactor the preconditioner at the point currently held by
    # the evaluation object
    function refreshpreconditioner!()
        fj!(nothing, J, x)
        tryfactorize!(cache, factorization, J)
        return nothing
    end

    # the residual norm at the initial point; every later entry of normF is
    # pushed immediately after a step is accepted, so convergence is decided
    # on each fresh residual and no preconditioner is ever assembled at a
    # final point
    residual!(F, x)
    push!(normF, norm(F))
    normF[end] <= ftol && (converged = true)

    for n in 1:iterations
        converged && break

        # the matrix-free product reads the evaluation point held by the
        # caller, and the linesearch leaves it at the last trial point
        # rather than the accepted one, so resynchronize it
        fj!(nothing, nothing, x)
        if refresh
            refreshpreconditioner!()
            refresh = false
        end

        # Eisenstat-Walker choice 2 forcing term from the last accepted
        # step, at its clamp maximum before any step has been taken
        forcing = if length(normF) >= 2 && normF[end-1] > 0
            clamp(krylovgamma*(normF[end]/normF[end-1])^krylovalpha,
                krylovrtolmin, krylovrtolmax)
        else
            krylovrtolmax
        end

        out = gmres!(deltax, jvp!, F, ws; Mop! = Mop!, rtol = forcing,
            atol = gmresatol, maxrestarts = krylovmaxrestarts)
        if !out.converged
            # the preconditioner has drifted far enough that GMRES no
            # longer converges. rebuild it and retry once, falling back
            # to the direct solve through the fresh factorization
            refreshpreconditioner!()
            out = gmres!(deltax, jvp!, F, ws; Mop! = Mop!,
                rtol = forcing, atol = gmresatol,
                maxrestarts = krylovmaxrestarts)
            out.converged || trysolve!(deltax, cache.factorization, F)
        end
        if out.iterations > krylovrefreshiterations
            refresh = true
        end
        rmul!(deltax, -1)

        # the merit function and its slope along deltax. the assembled
        # Jacobian is stale, so the slope comes from an exact matrix-free
        # product rather than the model claim dot(F, J, deltax) used by
        # nlsolve!
        ϕ0 = merit(F)
        jvp!(Jv, deltax)
        dϕ0dα = real(dot(F, Jv))
        if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
            # not a descent direction: fall back to the exact Newton step
            # through a fresh factorization, which is guaranteed descent up
            # to roundoff; if even that fails, stop
            refreshpreconditioner!()
            trysolve!(deltax, cache.factorization, F)
            rmul!(deltax, -1)
            jvp!(Jv, deltax)
            dϕ0dα = real(dot(F, Jv))
            if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
                break
            end
        end

        # interpolated backtracking linesearch shared with nlsolve!: on
        # Armijo failure it returns the best decreasing trial (alpha > 0)
        # with F and xcandidate restored there, or alpha == 0 when no trial
        # decreased the merit at all
        alpha1, ϕα, accepted, backtracks = backtracking_linesearch!(
            residual!, F, xcandidate, x, deltax, ϕ0, dϕ0dα;
            c1 = c1, safeguard_low = safeguard_low,
            safeguard_high = safeguard_high,
            maxbacktracks = maxbacktracks, Fbest = Fbest)
        push!(alphas, alpha1)
        push!(backtrackrecord, backtracks)

        if iszero(alpha1)
            # no decrease anywhere along the direction; F again holds the
            # residual at the unchanged x (the linesearch restore contract).
            # retry once from a fresh preconditioner, then give up
            refreshedforstall && break
            refreshedforstall = true
            refresh = true
            continue
        end
        refreshedforstall = false

        # accept the trial point: F already holds the residual there (the
        # linesearch postcondition), and convergence is decided on it now,
        # before any preconditioner work
        copyto!(x, xcandidate)
        push!(normF, norm(F))
        if normF[end] <= ftol
            converged = true
            break
        end

        # count consecutive linesearch failures where the best decreasing
        # trial was taken instead of an Armijo-accepted step; an accepted
        # step resets the count, repeated failures are a stall
        if accepted
            backtrackfailures = 0
        else
            backtrackfailures += 1
            backtrackfailures >= maxbacktrackfailures && break
        end
    end

    return IterationInfo(label, NaN, 0.0, converged, length(alphas), normF,
        alphas, backtrackrecord, fill(false, length(alphas)))
end
