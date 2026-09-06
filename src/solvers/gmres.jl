# GMRES with restarts over a reusable workspace, the Krylov vectors a Newton loop
# keeps between steps, and the linear solver objects (`GMRES`, `KrylovJL`) the
# Newton-Krylov driver dispatches on.

"""
    GMRESWorkspace{T<:AbstractFloat}

Preallocated storage for [`gmres!`](@ref) with a restart length of `m` on a
system of dimension `n`. Holds the `n x (m+1)` Arnoldi basis `V`, the
`(m+1) x m` Hessenberg matrix `H` as the Givens rotations leave it, the raw
Arnoldi Hessenberg `Harnoldi` beside it, the Givens rotations `cs` and `sn` which
reduce it, the least squares right hand side `s`, its solution `y`, three
length `n` work vectors, and the two length `m` staging buffers `hd` and
`cd` of the block Gram-Schmidt projection, allocated like `V`.

The dominant cost is `V`, which is `n*(m+1)` numbers, so `m` trades memory and
orthogonalization work against restart frequency. It is not paid up front:
the basis is allocated with a few columns and grows, by doubling, to what
the iteration uses, up to `m + 1`. A restart length long enough for the
hardest solve is then free on the easy ones, where a warm started Newton
step takes a handful of Arnoldi steps, and the cost of a solve is no longer
dominated by touching a basis it never fills. See [`ensurecolumns!`](@ref).
"""
mutable struct GMRESWorkspace{T<:AbstractFloat,TV<:AbstractVector{T},TM<:AbstractMatrix{T}}
    # system sized, and therefore allocated like the right hand side, so on a
    # device backend they live on the device. `V` holds the columns built so
    # far and is replaced by a wider one as the iteration needs them; `H`
    # is host resident and small, and is allocated in full.
    V::TM
    w::TV
    z::TV
    u::TV
    # the projected quantities, at most `m` by `m`, kept on the host: the
    # Givens rotations and the back substitution index them entry by entry,
    # which a device array must not be asked to do, and they are small
    # enough that the cost is one small transfer per Arnoldi step
    H::Matrix{T}
    # the Arnoldi Hessenberg before the Givens rotations, column by column
    # as each is finished. `H` itself is the triangularized least squares
    # matrix once the rotations have been applied; a harvest that needs the
    # Arnoldi relation `A*V[:, 1:j] = V[:, 1:j+1]*Harnoldi[1:j+1, 1:j]`
    # (the harmonic Ritz pencil) reads this one. The singular values of the
    # two agree, since a left orthogonal transformation preserves them.
    Harnoldi::Matrix{T}
    cs::Vector{T}
    sn::Vector{T}
    s::Vector{T}
    y::Vector{T}
    # staging buffers for the block Gram-Schmidt projections, allocated like
    # `V`, so that only the finished column of `H` crosses to the host
    hd::TV
    cd::TV
end

function GMRESWorkspace(n::Integer, m::Integer, ::Type{T} = Float64) where {T<:AbstractFloat}
    return GMRESWorkspace(Vector{T}(undef, n), m)
end

"""
    GMRESWorkspace(b::AbstractVector{T}, m::Integer)

Build a workspace for a restart length of `m` on a system whose right hand
side is `b`. The system sized arrays are allocated with `similar(b)`, so they
live wherever `b` does and the iteration runs on that device; the projected
`m` x `m` quantities are host arrays regardless, because the Givens rotations
and the back substitution index them entry by entry.
"""
function GMRESWorkspace(b::AbstractVector{T}, m::Integer) where {T<:AbstractFloat}
    n = length(b)
    n >= 0 || throw(ArgumentError(lazy"`n` = $(n) must be nonnegative."))
    m >= 1 || throw(ArgumentError(lazy"the restart length `m` = $(m) must be at least 1."))
    cols = min(m + 1, GMRESINITIALCOLUMNS)
    return GMRESWorkspace{T,typeof(similar(b)),typeof(similar(b, n, cols))}(
        similar(b, n, cols), similar(b), similar(b), similar(b),
        zeros(T, m + 1, m), zeros(T, m + 1, m),
        Vector{T}(undef, m), Vector{T}(undef, m),
        Vector{T}(undef, m + 1), Vector{T}(undef, m),
        similar(b, m), similar(b, m))
end

# the columns a basis is born with. Sixteen covers the Arnoldi steps of a
# well preconditioned or warm started Newton step without a copy, and a
# solve which needs more doubles its way there; at most log2(m/16) copies of
# what was built, which is a fraction of the orthogonalization that built it
const GMRESINITIALCOLUMNS = 16

"""
    ensurecolumns!(ws::GMRESWorkspace, k)

Return the basis `ws.V` with at least `k` columns, replacing it with a wider
one when it has fewer. The columns built so far are copied across; the new
ones are uninitialized, as the whole basis was before. The width doubles
and is capped at `size(ws.H, 1)`, which is `m + 1`, so a basis never grows
past the restart length.
"""
function ensurecolumns!(ws::GMRESWorkspace, k::Integer)
    V = ws.V
    have = size(V, 2)
    k <= have && return V
    cols = min(max(k, 2*have), size(ws.H, 1))
    Vnew = similar(V, size(V, 1), cols)
    copyto!(view(Vnew, :, 1:have), V)
    ws.V = Vnew
    return Vnew
end


"""
    KrylovVectors(x, F, m)

The system sized vectors of one Newton-Krylov solve: the GMRES workspace
for a restart length of `m`, and the step, the trial point, the product
and the best residual, allocated like `x` and `F`. Handed back to
[`nlsolvekrylov!`](@ref) through its `workspace` argument they are reused
across solves of one system, which a sweep over component values is.
"""
struct KrylovVectors{V,VF,W}
    ws::W
    deltax::V
    xcandidate::V
    Jv::V
    Fbest::VF
end

KrylovVectors(x::AbstractVector, F::AbstractVector, m::Integer) =
    KrylovVectors(GMRESWorkspace(x, m), similar(x), similar(x), similar(x),
        similar(F))

"""
    harvest!(pc::AbstractPreconditioner, ws::GMRESWorkspace, out::NamedTuple)

Give the preconditioner `pc` the Arnoldi factorization a solve just built, so
it can extract information for the *next* solve. `out` is the named tuple
returned by [`gmres!`](@ref). Called by [`nlsolvekrylov!`](@ref) after every
GMRES call. The default does nothing, which is correct for any preconditioner
that does not recycle.

Only the *last* restart cycle is still present in the workspace, so
implementations must derive the usable Arnoldi dimension from
`out.iterations` and `out.cycles` rather than from `out.iterations` alone.
"""
harvest!(pc::AbstractPreconditioner, ::GMRESWorkspace, ::NamedTuple) = pc
harvest!(pc::AbstractWrappedPreconditioner, ws::GMRESWorkspace, out::NamedTuple) =
    (harvest!(innerpreconditioner(pc), ws, out); pc)

"""
    harvestcycle!(pc::AbstractPreconditioner, ws::GMRESWorkspace, j::Integer)

Give the preconditioner the Arnoldi factorization of the restart cycle which
has just ended, `j` vectors of it, before [`gmres!`](@ref) overwrites the
workspace with the next cycle. The default does nothing.

This is the per-cycle counterpart of [`harvest!`](@ref), which sees only the
cycle left in the workspace when the solve returns. A preconditioner opts
into it through [`usescycleharvest`](@ref), and one which does is *not*
harvested again afterwards.
"""
harvestcycle!(pc::AbstractPreconditioner, ::GMRESWorkspace, ::Integer) = pc
harvestcycle!(pc::AbstractWrappedPreconditioner, ws::GMRESWorkspace, j::Integer) =
    (harvestcycle!(innerpreconditioner(pc), ws, j); pc)

"""
    usescycleharvest(pc::AbstractPreconditioner)

Whether `pc` wants [`harvestcycle!`](@ref) at the end of every restart cycle
instead of [`harvest!`](@ref) once the solve is over. `false` by default, so
that a preconditioner which harvests only the final cycle keeps doing
exactly that.
"""
usescycleharvest(::AbstractPreconditioner) = false

"""
    isexactpreconditioner(pc::AbstractPreconditioner)

Whether `pc` currently applies the exact Jacobian, so that a deflation or
composition layered on top of it can contribute nothing. `false` for any
preconditioner that does not say otherwise.
"""
isexactpreconditioner(::AbstractPreconditioner) = false

"""
    gmres_orthogonalize!(w, V, H, hd, c, j)

Orthogonalize `w` against the first `j` Arnoldi basis vectors (the columns of
`V`) by block classical Gram-Schmidt with one reorthogonalization (CGS2),
accumulating the coefficients into column `j` of the Hessenberg matrix `H`.
`c` is scratch of length at least `j`. Both passes accumulate into the same
entries of `H`, so `H` remains the exact projection. Writes the subdiagonal
`H[j+1, j]` and returns `(hsub, normw0)`: the norm of the orthogonalized `w`
and its norm on entry, the pair the caller compares to detect a breakdown.

CGS2 is chosen over modified Gram-Schmidt with a DGKS test for its shape
rather than its accuracy, which is equivalent: both are orthogonal to machine
precision. Each pass here is two level 2 BLAS calls over the whole basis,
where the modified form is `j` dependent pairs of a dot product and an axpy,
each dot having to complete before the axpy that follows it. The second pass
is unconditional, which costs what the DGKS path costs whenever it does
reorthogonalize and removes a branch on a freshly computed scalar. Neither the
coefficient vector nor the branch has to reach the host, which is what makes
this form usable on a device.
"""
function gmres_orthogonalize!(w::AbstractVector{T}, V::AbstractMatrix{T},
    H::AbstractMatrix{T}, hd::AbstractVector{T}, c::AbstractVector{T},
    j::Integer) where {T<:AbstractFloat}
    normw0 = norm2(w)
    if j > 0
        Vj = view(V, :, 1:j)
        # the projections are formed in buffers allocated like `V`, so on a
        # device they stay there; only the finished column crosses to `H`
        hj = view(hd, 1:j)
        cj = view(c, 1:j)
        # Block classical Gram-Schmidt, two unconditional passes. CGS2 is as
        # accurate as modified Gram-Schmidt with a DGKS test, to machine
        # precision, and it is a much better shape for a GPU: each pass is two
        # level 2 BLAS calls over the whole basis rather than `j` dependent
        # pairs of a dot product and an axpy. In the modified form every dot
        # has to finish before the axpy that follows it, which on a device
        # means `j` synchronizations per Arnoldi step; here the coefficient
        # vector never has to reach the host. Making the second pass
        # unconditional also removes a branch on a device-resident scalar,
        # which would force a synchronization of its own. On the CPU the two
        # passes cost what the DGKS path costs when it does reorthogonalize.
        mul!(hj, transpose(Vj), w)
        mul!(w, Vj, hj, -one(T), one(T))
        mul!(cj, transpose(Vj), w)
        mul!(w, Vj, cj, -one(T), one(T))
        hj .+= cj
        # `H` is deliberately host resident. copyto! between a host view and a
        # device view falls back to scalar indexing of the device array, so the
        # column crosses through its contiguous linear range instead: a dense
        # `H` stores rows 1:j of column j contiguously from (j-1)*size(H,1)+1.
        copyto!(H, (j-1)*size(H,1)+1, hd, 1, j)
    end
    hsub = norm2(w)
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

    # The projected problem can be rank deficient, on a singular or
    # inconsistent system or after a near breakdown. Back substituting
    # through a tiny pivot then amplifies roundoff without bound: a pivot of
    # order eps relative to the rest of the triangle produces a coefficient
    # of order 1/eps, and the returned "solution" is dominated by a spurious
    # basis direction. (Testing exactly for a zero pivot does not help, since
    # the damaging case is a pivot which is small but nonzero.) Instead the
    # numerically well determined leading part of the triangle is solved and
    # the remaining coefficients are set to zero, which is the minimizer over
    # the subspace the cycle actually resolved.
    rank = j
    dmax = zero(T)
    for i in 1:j
        dmax = max(dmax, abs(H[i, i]))
    end
    pivottol = eps(T)*max(dmax, one(T))*j
    for i in 1:j
        if abs(H[i, i]) <= pivottol
            rank = i - 1
            break
        end
    end
    for i in rank+1:j
        y[i] = zero(T)
    end
    for i in rank:-1:1
        acc = s[i]
        for k in i+1:rank
            acc -= H[i, k]*y[k]
        end
        y[i] = acc/H[i, i]
    end
    j = rank
    fill!(u, zero(T))
    for i in 1:j
        axpy!(y[i], view(V, :, i), u)
    end
    if isnothing(Mop!)
        @. x += u
    else
        _applyprecond!(z, Mop!, u)
        @. x += z
    end
    return x
end

_applyprecond!(z::AbstractVector, M::AbstractPreconditioner,
    v::AbstractVector) = applypreconditioner!(z, M, v)
_applyprecond!(z::AbstractVector, M, v::AbstractVector) = (M(z, v); z)

"""
    preconditionedproduct!(w, z, Aop, M, v)

One right preconditioned Arnoldi step: overwrite `z` with `inv(M)*v` and `w`
with `A*z`, and return the seconds spent inside the preconditioner. The timing
excludes the Jacobian product so that `precondtime` in
[`KrylovSolveInfo`](@ref) means the same thing whatever the preconditioner.
(A fused form which folded the deflation's correction into this product
existed and was retired: it measured more Arnoldi steps than the image
pair form of [`FloquetPreconditioner`](@ref) on every case tried.)
"""
function preconditionedproduct!(w::AbstractVector, z::AbstractVector, Aop,
    M, v::AbstractVector)
    return _unfusedproduct!(w, z, Aop, M, v)
end

function _unfusedproduct!(w, z, Aop, M, v)
    tpc = time()
    _applyprecond!(z, M, v)
    tpc = time() - tpc
    mul!(w, Aop, z)
    return tpc
end

"""
    harvestdimension(ws::GMRESWorkspace, out::NamedTuple)

The number of Arnoldi vectors of the *last* restart cycle still present in
the workspace, which is the usable dimension for a harvest; the workspace
holds only that cycle, so this is derived from `out.iterations` and
`out.cycles` rather than from the iteration count alone. Zero when the
cycle is empty or overran.
"""
function harvestdimension(ws::GMRESWorkspace, out::NamedTuple)
    m = size(ws.H, 2)
    j = Int(out.iterations) - (Int(out.cycles) - 1)*m
    return 1 <= j <= m ? j : 0
end

"""
    norm2(v::AbstractVector)

The Euclidean norm of `v`, formed through the inner product.

`LinearAlgebra.norm` scales its argument before squaring so that an entry
cannot overflow or underflow on its way to the sum. That guard is not free
on a device: cuBLAS routes a Float64 vector to a scaled `nrm2` kernel which
runs at a small fraction of memory bandwidth, measured at 22.9 us against
2.7 us for a `dot` product of the same 51,200 element vector on an RTX 4090,
and it is the largest single device kernel of a double precision solve. In
single precision cuBLAS selects a different kernel and the gap is gone.

The substitution is made only where it pays and only where it is safe, which
is the same place: `Float64`.

- In double precision cuBLAS routes `norm` to a scaled `nrm2` kernel that runs
  at well under a tenth of memory bandwidth, measured at 55.6 us against
  14.5 us for `sqrt(dot(v, v))` on a 51,200 element device vector, and it was
  the largest single device kernel of a double precision solve.
- In single precision cuBLAS selects a different kernel and the two are within
  20% of each other, so there is nothing to win. Single precision is also
  where the exponent range is narrowest and the guard is worth the most.

So `Float32` keeps `norm` and everything else goes through the inner product,
and the change is confined to the precision where the trade is favorable in
both directions. Anything not covered by the `AbstractFloat` method -- complex
vectors, other element types -- falls back to `norm` as well.

The vectors this is applied to are carried at the scale of the harmonic
balance residual of a nondimensionalized system, which runs from a few units
down to the solver tolerance. Squaring that stays far inside the double
precision exponent: overflow would need an entry above 1e154 and underflow to
zero an entry below 1e-162.

This is deliberately not exported and not used outside the Krylov solver.
Anything whose scale is not controlled should keep using `norm`.
"""
norm2(v::AbstractVector{<:AbstractFloat}) = sqrt(dot(v, v))
norm2(v::AbstractVector{Float32}) = norm(v)
norm2(v::AbstractVector) = norm(v)

"""
    gmres!(x, Aop!, b, ws::GMRESWorkspace; Mop! = nothing, rtol = 1e-6,
        atol = 0.0, maxrestarts = 10, initialzero = true, oncycle = nothing)

Solve `A*x = b` with restarted GMRES, where `mul!(w, Aop, v)` computes `w = A*v` and
the optional `Mop!` applies a preconditioner `z = M \\ v`, either as a bare
in-place closure `Mop!(z, v)` or as an [`AbstractPreconditioner`](@ref), which
is applied through [`applypreconditioner!`](@ref) and may fuse its application
with the operator product ([`preconditionedproduct!`](@ref)). The matrix `A`
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
tuple `(iterations, residual, converged, cycles, reason, precondtime,
residualvector, products)`, where `iterations` counts Arnoldi steps across
all cycles, `cycles` the number of restart cycles begun, `reason` is one of
`:converged`, `:breakdown` (an unhappy breakdown: the Krylov space went
invariant without the residual coming down), `:stagnation` (a cycle failed
to reduce the explicit residual, or produced a non-finite one), or
`:iterationlimit`, `precondtime` the seconds spent applying the
preconditioner, and `residualvector` the final residual `b - A*x` when it
was formed explicitly (`nothing` otherwise; the caller reads it with `get`,
as it does `precondtime` and `products`).

`iterations` is *not* the total number of `Aop!` calls: each cycle costs one
further application for the explicit residual recomputation, and a warm start
costs one at the outset; `products` in the returned tuple is that total, not
counting products a preconditioner takes inside its own application. `maxrestarts` bounds the number of cycles including
the first, so the Arnoldi work is capped at `maxrestarts*m` steps.
`oncycle(ws, j)`, when given, is called at the end of every cycle with the
workspace still holding that cycle's `j` Arnoldi vectors, for a caller
which harvests from each cycle ([`harvestcycle!`](@ref)); it must only read
the workspace.

Allocation free after the workspace is built, apart from whatever `Aop!` and
`Mop!` themselves allocate.
"""
function gmres!(x::AbstractVector{T}, Aop_, b::AbstractVector{T},
    ws::GMRESWorkspace{T}; Mop! = nothing, rtol = 1e-6, atol = 0.0,
    maxrestarts::Integer = 10, initialzero::Bool = true,
    oncycle = nothing) where {T<:AbstractFloat}

    n = length(b)
    # a bare in-place product is accepted alongside any `mul!`-able operator
    Aop = asoperator(Aop_, n)
    precondtime = 0.0
    products = 0
    length(x) == n || throw(DimensionMismatch(
        lazy"`x` has length $(length(x)) but `b` has length $(n)."))
    size(ws.V, 1) == n || throw(DimensionMismatch(
        lazy"the workspace is for dimension $(size(ws.V,1)) but `b` has length $(n)."))
    (rtol >= 0 && isfinite(rtol)) || throw(ArgumentError(
        lazy"`rtol` = $(rtol) must be nonnegative and finite."))
    (atol >= 0 && isfinite(atol)) || throw(ArgumentError(
        lazy"`atol` = $(atol) must be nonnegative and finite."))
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
        # zero is a solution, but so is any point already in the null space
        # of A, so a warm start is measured rather than discarded
        if initialzero
            fill!(x, zero(T))
            return (iterations = 0, residual = zero(T), converged = true,
                cycles = 0, reason = :converged, precondtime = 0.0,
                residualvector = nothing, products = 0)
        end
        mul!(w, Aop, x)
        products += 1
        resnorm = norm(w)
        if resnorm <= atol
            return (iterations = 0, residual = resnorm, converged = true,
                cycles = 0, reason = :converged, precondtime = 0.0,
                residualvector = nothing, products = products)
        end
    end

    # initial residual w = b - A*x
    if initialzero
        fill!(x, zero(T))
        copyto!(w, b)
    else
        mul!(w, Aop, x)
        products += 1
        @. w = b - w
    end
    resnorm = norm(w)

    totaliterations = 0
    cycles = 0
    unhappy = false
    stagnated = false
    for _ in 1:maxrestarts
        resnorm <= tol && break

        cycles += 1
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
                mul!(w, Aop, z)
            else
                precondtime += preconditionedproduct!(w, z, Aop, Mop!,
                    view(V, :, j))
            end

            hsub, normw0 = gmres_orthogonalize!(w, V, H, ws.hd, ws.cd, j)
            # the finished column of the Arnoldi Hessenberg, before the
            # rotations below triangularize it in place
            copyto!(view(ws.Harnoldi, 1:j+1, j), view(H, 1:j+1, j))
            resnorm = gmres_applyrotations!(H, cs, sn, s, j)
            totaliterations += 1
            products += 1

            # A subdiagonal which collapsed relative to the incoming vector
            # means the Krylov space is invariant and there is no valid next
            # basis vector to normalize. This is not automatically a *lucky*
            # breakdown: that additionally requires the projected problem to
            # be compatible, so that the reduced solution really is the
            # solution. On a singular or inconsistent system the space is
            # invariant while the residual is not reducible in it at all, an
            # unhappy breakdown, and continuing to restart only rebuilds the
            # same useless space. The two are separated by whether the
            # recurrence residual actually came down.
            breakdown = hsub <= eps(T)*normw0
            if breakdown
                unhappy = resnorm > tol && resnorm > (1 - sqrt(eps(T)))*beta
                break
            end
            resnorm <= tol && break

            V = ensurecolumns!(ws, j + 1)
            @views V[:, j+1] .= w ./ hsub
        end

        gmres_correction!(x, ws, j, Mop!)

        # The Arnoldi factorization of this cycle is about to be overwritten
        # by the next one, so anything that wants to read it has to read it
        # here. A caller which recycles uses this to harvest from *every*
        # cycle rather than only the one left in the workspace at the end:
        # the difficult directions often show up in an early full cycle,
        # while the cycle that finally converges can be a few steps long and
        # carry almost nothing. The callback may only *read* the workspace;
        # the preconditioner is fixed for the duration of a solve, so a
        # harvest appends candidates and nothing is rebuilt until the solve
        # is over.
        isnothing(oncycle) || oncycle(ws, j)

        # recompute the residual explicitly for the next cycle so restarts
        # cannot drift from the recurrence estimate
        mul!(w, Aop, x)
        products += 1
        @. w = b - w
        previous = beta
        resnorm = norm(w)

        # a cycle which did not reduce the explicit residual will not do
        # better on a rebuild of the same space, so stop rather than burn the
        # remaining cycle budget on it
        if !isfinite(resnorm)
            stagnated = true
            break
        end
        if unhappy || resnorm >= (1 - sqrt(eps(T)))*previous
            stagnated = true
            break
        end
    end

    converged = resnorm <= tol
    reason = if converged
        :converged
    elseif unhappy
        :breakdown
    elseif stagnated
        :stagnation
    else
        :iterationlimit
    end
    # `w` is the explicit residual `b - A x` of the returned `x`: it is
    # recomputed after every cycle and never estimated, so a caller which
    # needs `A x` has it without another product
    return (iterations = totaliterations, residual = resnorm,
        converged = converged, cycles = cycles, reason = reason,
        precondtime = precondtime, residualvector = w, products = products)
end

# `AbstractHBLinearSolver` is declared in solvers/options.jl, before the
# `NewtonKrylov` option which holds one.

"""
    GMRES(; restart = 400, maxrestarts = 4)

The restarted GMRES of this package, with Givens rotations, the recycling
subspace harvest and the preconditioner escalation the solver was built
around. The default, and the only solver supporting deflation recycling,
because `harvest!` reads the Arnoldi basis out of the internal workspace.
`restart` is the cycle length and `maxrestarts` the restart budget per
solve. A long cycle is the default because a restricted preconditioner
leaves a few directions a short Krylov space cannot resolve, and a restart
discards the progress on them; the basis of `restart + 1` vectors is cheap
next to the sparse factorization an escalation would build.
"""
struct GMRES <: AbstractHBLinearSolver
    restart::Int
    maxrestarts::Int
end
function GMRES(; restart::Integer = 400, maxrestarts::Integer = 4)
    restart >= 1 || throw(ArgumentError(lazy"`restart` = $(restart) must be at least 1."))
    maxrestarts >= 1 || throw(ArgumentError(
        lazy"`maxrestarts` = $(maxrestarts) must be at least 1."))
    return GMRES(Int(restart), Int(maxrestarts))
end
# the Krylov workspace of an external solver is sized like the default
restartlength(ls::GMRES) = ls.restart
restartlength(::AbstractHBLinearSolver) = 400
maxrestarts(ls::GMRES) = ls.maxrestarts
maxrestarts(::AbstractHBLinearSolver) = 4

"""
    KrylovJL(method::Symbol = :gmres; kwargs...)

Solve the Newton step with Krylov.jl: `:gmres`, `:fgmres`, `:bicgstab`,
`:dqgmres` or any other solver taking an operator and a right hand side.
Requires Krylov.jl to be loaded; the method lives in the package extension.

Only the linear solve changes. The forcing term, the line search, the
preconditioner escalation and the stagnation handling are untouched, and
the mode coupling preconditioner is passed through unchanged because it is
applied by `mul!`, which is what Krylov.jl's `N` argument consumes.
Deflation recycling is unavailable, since it depends on the internal
workspace.
"""
struct KrylovJL <: AbstractHBLinearSolver
    method::Symbol
    kwargs::NamedTuple
end
KrylovJL(method::Symbol = :gmres; kwargs...) = KrylovJL(method, NamedTuple(kwargs))

"""
    hblinearsolve!(ls, deltax, jvp!, F, ws, Mop!; rtol, atol, maxrestarts,
        oncycle = nothing)

Solve for the Newton step and return the output named tuple `gmres!`
produces. `converged`, `residual`, `iterations`, `cycles` and `reason` are
required; `residualvector` (the explicit final residual, which the line
search slope reads, at the cost of one extra Jacobian product when it is
missing), `precondtime` and `products` are read with `get` and may be
omitted by an external solver.
`jvp!(y, v)` applies the Jacobian, in place, and `Mop!` is the
preconditioner, an [`AbstractPreconditioner`](@ref) or a closure
`Mop!(z, r)`. `oncycle` is the per-cycle callback of [`gmres!`](@ref),
which a solver that does not restart may ignore.

This is the one part of the Newton-Krylov loop with nothing harmonic
balance specific about it: an operator, a right hand side, a preconditioner
and a tolerance. Putting it behind an interface lets an external Krylov
library be substituted without touching anything else.
"""
function hblinearsolve!(::GMRES, deltax, jvp!, F, ws, Mop!;
        rtol, atol, maxrestarts, oncycle = nothing)
    return gmres!(deltax, jvp!, F, ws; Mop! = Mop!, rtol = rtol,
        atol = atol, maxrestarts = maxrestarts, oncycle = oncycle)
end

hblinearsolve!(ls::AbstractHBLinearSolver, args...; kwargs...) =
    throw(ArgumentError(lazy"no hblinearsolve! method for $(typeof(ls)); "*
        "KrylovJL requires Krylov.jl to be loaded."))

"""
    supportsrecycling(ls)

Whether the solver exposes an Arnoldi basis for [`harvest!`](@ref).
"""
supportsrecycling(::AbstractHBLinearSolver) = false
supportsrecycling(::GMRES) = true
