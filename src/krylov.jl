
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
    V, H, cs, sn, s, y = ws.V, ws.H, ws.cs, ws.sn, ws.s, ws.y
    w, z, u = ws.w, ws.z, ws.u

    bnorm = norm(b)
    tol = max(rtol * bnorm, atol)

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

            # modified Gram-Schmidt, which is stable enough at the restart
            # lengths used here
            for i in 1:j
                @views H[i, j] = dot(V[:, i], w)
                @views axpy!(-H[i, j], V[:, i], w)
            end
            H[j+1, j] = norm(w)

            # a zero subdiagonal is a lucky breakdown: the Krylov space is
            # invariant and the least squares solution below is exact
            if H[j+1, j] > eps(T) * beta
                @views V[:, j+1] .= w ./ H[j+1, j]
            else
                @views fill!(V[:, j+1], zero(T))
            end

            # apply the previous rotations to the new column
            for i in 1:j-1
                tmp       =  cs[i] * H[i, j] + sn[i] * H[i+1, j]
                H[i+1, j] = -sn[i] * H[i, j] + cs[i] * H[i+1, j]
                H[i, j]   = tmp
            end

            # the rotation which annihilates the new subdiagonal
            r = hypot(H[j, j], H[j+1, j])
            if iszero(r)
                cs[j] = one(T)
                sn[j] = zero(T)
            else
                cs[j] = H[j, j] / r
                sn[j] = H[j+1, j] / r
            end
            H[j, j]   = cs[j] * H[j, j] + sn[j] * H[j+1, j]
            H[j+1, j] = zero(T)
            s[j+1] = -sn[j] * s[j]
            s[j]   =  cs[j] * s[j]

            totaliterations += 1
            # with right preconditioning this is the true residual norm
            resnorm = abs(s[j+1])
            resnorm <= tol && break
        end

        # back substitute the triangular least squares problem
        for i in j:-1:1
            acc = s[i]
            for k in i+1:j
                acc -= H[i, k] * y[k]
            end
            y[i] = iszero(H[i, i]) ? zero(T) : acc / H[i, i]
        end

        # u = V[:,1:j]*y, then undo the right preconditioning once
        fill!(u, zero(T))
        for i in 1:j
            @views axpy!(y[i], V[:, i], u)
        end
        if isnothing(Mop!)
            @. x += u
        else
            Mop!(z, u)
            @. x += z
        end

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
    nlsolvekrylov!(fj!, jvp!, F, J, x; iterations = 1000, ftol = 1e-8,
        factorization = KLUfactorization(), label = "", c1 = 1e-4,
        backtrackfactor = 0.5, maxbacktracks = 10, krylovrestart = 30,
        krylovmaxrestarts = 4, krylovrefreshiterations = 10,
        krylovrtolmin = 1e-10, krylovrtolmax = 1e-2, krylovgamma = 0.9,
        krylovalpha = (1 + sqrt(5))/2)

Inexact Newton solver for a real system, taking the Newton step from
[`gmres!`](@ref) on the matrix-free product `jvp!(y, v)` rather than from a
factorization of an assembled Jacobian. `fj!(F, J, x)` evaluates the residual
and Jacobian as in [`nlsolve!`](@ref), accepting `nothing` for either.

The assembled Jacobian `J` is used only as a preconditioner, and is
reassembled and refactored when the Krylov iteration count says it has gone
stale, when GMRES fails, or when the linesearch cannot make progress. The
linear tolerance follows an Eisenstat-Walker forcing sequence between
`krylovrtolmin` and `krylovrtolmax`, so early steps are not solved to an
accuracy the nonlinear iteration cannot use. The linesearch is plain Armijo
backtracking, and the slope comes from an exact matrix-free product because
the assembled Jacobian is stale.

Returns an [`IterationInfo`](@ref).
"""
function nlsolvekrylov!(fj!::Function, jvp!::Function, F::AbstractVector{T},
    J::AbstractArray{T}, x::Vector{T}; iterations = 1000, ftol = 1e-8,
    factorization = KLUfactorization(), label = "", c1 = 1e-4,
    backtrackfactor = 0.5, maxbacktracks::Integer = 10,
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
    iterations >= 0 || throw(ArgumentError(
        lazy"`iterations` = $(iterations) must be nonnegative."))
    ftol >= 0 || throw(ArgumentError(lazy"`ftol` = $(ftol) must be nonnegative."))
    0 < c1 < 1//2 || throw(ArgumentError(
        lazy"`c1` = $(c1) must be in (0, 1/2) for the Newton merit function."))
    0 < backtrackfactor < 1 || throw(ArgumentError(
        lazy"`backtrackfactor` = $(backtrackfactor) must be in (0, 1)."))
    maxbacktracks >= 0 || throw(ArgumentError(
        lazy"`maxbacktracks` = $(maxbacktracks) must be nonnegative."))
    krylovrestart >= 1 || throw(ArgumentError(
        lazy"`krylovrestart` = $(krylovrestart) must be at least 1."))
    krylovrefreshiterations >= 1 || throw(ArgumentError(
        lazy"`krylovrefreshiterations` = $(krylovrefreshiterations) must be at least 1."))
    0 < krylovrtolmin <= krylovrtolmax < 1 || throw(ArgumentError(
        lazy"`krylovrtolmin` = $(krylovrtolmin) and `krylovrtolmax` = $(krylovrtolmax) must satisfy 0 < krylovrtolmin <= krylovrtolmax < 1."))

    cache = FactorizationCache()
    ws = GMRESWorkspace(length(x), min(krylovrestart, length(x)), T)
    deltax = similar(x)
    xcandidate = similar(x)
    Jv = similar(x)

    normF = real(T)[]
    alphas = real(T)[]
    backtrackrecord = Int[]
    krylovrecord = Int[]
    converged = false
    refresh = true

    Mop!(zv, vv) = trysolve!(zv, cache.factorization, vv)

    fj!(F, nothing, x)
    push!(normF, norm(F))
    normF[end] <= ftol && (converged = true)

    for _ in 1:iterations
        converged && break

        # the matrix-free product reads the evaluation point held by the
        # caller, and the linesearch leaves it at the last trial point rather
        # than the accepted one, so resynchronize it
        fj!(nothing, nothing, x)
        if refresh
            fj!(nothing, J, x)
            tryfactorize!(cache, factorization, J)
            refresh = false
        end

        forcing = if length(normF) >= 2 && normF[end-1] > 0
            clamp(krylovgamma*(normF[end]/normF[end-1])^krylovalpha,
                krylovrtolmin, krylovrtolmax)
        else
            krylovrtolmax
        end

        out = gmres!(deltax, jvp!, F, ws; Mop! = Mop!, rtol = forcing,
            maxrestarts = krylovmaxrestarts)
        if !out.converged
            # the preconditioner has drifted far enough that GMRES no longer
            # converges. rebuild it here and retry once, falling back to the
            # direct solve through the fresh factorization
            fj!(nothing, J, x)
            tryfactorize!(cache, factorization, J)
            out = gmres!(deltax, jvp!, F, ws; Mop! = Mop!, rtol = forcing,
                maxrestarts = krylovmaxrestarts)
            out.converged || trysolve!(deltax, cache.factorization, F)
        elseif out.iterations > krylovrefreshiterations
            refresh = true
        end
        push!(krylovrecord, out.iterations)
        rmul!(deltax, -1)

        ϕ0 = real(T)(0.5)*normF[end]^2
        # the assembled Jacobian is stale, so take the slope from an exact
        # matrix-free product
        jvp!(Jv, deltax)
        dϕ0dα = real(dot(F, Jv))
        if !isfinite(ϕ0) || !isfinite(dϕ0dα) || dϕ0dα >= zero(dϕ0dα)
            # not a descent direction: refresh the preconditioner and, if it
            # was already fresh, stop
            refresh && break
            refresh = true
            continue
        end

        # Armijo backtracking
        alpha = one(real(T))
        accepted = false
        backtracks = 0
        for _ in 0:maxbacktracks
            @. xcandidate = x + alpha*deltax
            fj!(F, nothing, xcandidate)
            ϕ = real(T)(0.5)*norm(F)^2
            if isfinite(ϕ) && ϕ <= ϕ0 + c1*alpha*dϕ0dα
                accepted = true
                break
            end
            alpha *= backtrackfactor
            backtracks += 1
        end

        if !accepted
            # no sufficient decrease along this direction. a stale
            # preconditioner is the likely cause, so refresh once and retry
            fj!(F, nothing, x)
            if refresh
                break
            end
            refresh = true
            continue
        end

        copyto!(x, xcandidate)
        push!(alphas, alpha)
        push!(backtrackrecord, backtracks)
        push!(normF, norm(F))
        if normF[end] <= ftol
            converged = true
            break
        end
    end

    return IterationInfo(label, NaN, 0.0, converged, length(alphas), normF,
        alphas, backtrackrecord, fill(false, length(alphas)))
end
