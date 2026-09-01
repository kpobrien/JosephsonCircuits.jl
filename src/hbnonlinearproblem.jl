# =====================================================================
# The harmonic balance nonlinear system as a standalone, evaluable object.
#
# `hbnlsolve` builds the residual, the matrix-free Jacobian-vector product
# and the mode layout, uses them once, and throws them away. Exposing them
# lets an external nonlinear solver drive the same system: the residual and
# the products are what any of them ask for.
#
# The equivalent *real* representation is the vector space to hand out. The
# harmonic balance residual is not complex differentiable, so the implicit
# function theorem does not hold with the holomorphic Jacobian, while in
# the real representation it applies directly. Every external solver which
# relies on a Jacobian -- Newton, trust region, pseudo-transient -- needs
# the real form.
# =====================================================================

"""
    HBNonlinearProblem

The harmonic balance nonlinear system in the equivalent real
representation, as an evaluable object: `F(u)`, the exact matrix-free
`J(u)*v`, and the assembled real Jacobian when one was built.

Built by [`hbnonlinearproblem`](@ref). See [`hbresidual!`](@ref) and
[`hbjvp!`](@ref).

# Fields
- `sys`: the [`HBSystem`](@ref) evaluation object.
- `modelayout`: the real representation layout.
- `u0`: the initial value, in the real representation.
- `jacobian`: the assembled real Jacobian, or `nothing` when built with
  `assemblejacobian = false`, which is what a matrix-free solver wants: on
  a multi-tone problem the real Jacobian plan is the single largest object
  in the solve.
"""
# =====================================================================
# The direct current augmentation.
#
# With an explicit direct current block the system being solved is not the
# `HBSystem`: its scattering rows still say `i = 0` and its resistors still
# carry no direct current. The average voltages, the transport rows and the
# blocks' own zero frequency rows live in canonical coordinates around it.
#
# So there are two candidate definitions of the same operating point, and
# handing out the harmonic one is how a caller ends up differentiating a
# problem which is not the one that was solved. The augmentation is carried
# on the problem rather than in a second problem type: one object, one
# interface, and an external solver written against it does not have to know
# whether the circuit has direct current in it.
#
# Everything the augmentation adds is affine. Writing `G` for the gather
# into canonical coordinates and `S = G'` for the scatter back, the residual
# is
#
#     R(u) = D (G F(S u)) + M u + s c
#
# with `M` the block's constant matrix, `c` its drive dependent constant,
# `s` the drive scale, and `D` the diagonal which is zero on the rows the
# block replaces rather than adds to -- the transport rows, which the
# harmonic residual knows nothing about, and the reference rows. Every
# derivative follows from that: the Jacobian is `D G J S + M`, the
# transposed product is `S' J' G' D` plus `M'`, and the second and third
# derivatives lose `M` entirely because it is linear.

"""
    DCAugmentation

The explicit direct current block of an [`HBNonlinearProblem`](@ref).

# Fields
- `work`: the [`CanonicalWork`](@ref) holding the layout, the transport rows
  and the blocks' zero frequency rows.
- `jplan`: the [`CanonicalJacobianPlan`](@ref), or `nothing` when no
  Jacobian was assembled.
- `jint`: the internal Jacobian the plan reads, or `nothing`.
- `keep`: the diagonal `D` over the direct current window: one where a row
  is added to and zero where the block replaces it.
- `constant`: the block's drive dependent constant, at the drive the problem
  was built with.
- `dcmatrix`: the transpose of the block's constant matrix `M`, restricted
  to the window, which is what the transposed product needs.
- `scale`: the drive scale, kept in step with [`setdrive!`](@ref).
- `dwork`, `Fwork`, `zwork`: workspaces for the transposed product.
"""
struct DCAugmentation{W,JP,JI}
    work::W
    jplan::JP
    jint::JI
    keep::Vector{Float64}
    constant::Vector{Float64}
    dcmatrix::SparseMatrixCSC{Float64,Int}
    scale::Base.RefValue{Float64}
    dwork::Vector{Float64}
    Fwork::Vector{Float64}
    zwork::Vector{Float64}
end

function DCAugmentation(work, jint)
    L = work.layout
    N = canonicaldim(L)
    nw = L.ndc + L.nvdc
    # `keep`, the constant and the matrix are read off the scalar
    # implementation by probing it, exactly as the device form is, so the
    # three cannot disagree with the residual they describe
    up = dcupdate(work)
    isnothing(up) && throw(ArgumentError("a direct current augmentation needs a direct current block"))
    Mt = SparseMatrixCSC(nw, nw, copy(up.rowptr), copy(up.colval),
        copy(up.nzval))
    keep = ones(Float64, N)
    constant = zeros(Float64, N)
    win = dcwindow(L)
    copyto!(view(keep, win), up.keep)
    copyto!(view(constant, win), up.cresidual)
    jplan = isnothing(jint) ? nothing : canonicaljacobianplan(jint, work)
    return DCAugmentation(work, jplan, jint, keep, constant, Mt,
        Ref(1.0), zeros(Float64, nw), zeros(Float64, L.rdim),
        zeros(Float64, N))
end

"""
    HBNonlinearProblem

The harmonic balance system as a nonlinear problem: the residual, its exact
derivatives and the pieces a preconditioner reads, with the unknowns a plain
real vector so that a solver written elsewhere can drive it.

Build one with [`hbnonlinearproblem`](@ref). `length(p)` is the number of
unknowns; [`hbresidual!`](@ref), [`hbjvp!`](@ref), [`hbjacobian!`](@ref) and
[`hbvjp!`](@ref) evaluate at any point, and [`setdrive!`](@ref) scales the
drive for continuation.

When the circuit injects direct current the unknowns are the canonical state
rather than the harmonic one, and `augmentation` carries the block which
makes the difference. [`isaugmented`](@ref) says which, and every entry
point above routes through the block when it is there, so a caller does not
have to know. See [`DCAugmentation`](@ref).

# Fields
- `sys`: the harmonic balance system the residual is evaluated on.
- `modelayout`: the layout relating the real unknowns to the complex modes.
- `u0`: the point the problem was built at.
- `jacobian`: the assembled Jacobian, or `nothing` when none was asked for.
- `parts`: the pieces a preconditioner needs.
- `bnm0`: the drive as built, which `setdrive!` scales.
- `tplan`: the transposed gather maps the vector-Jacobian product walks.
- `Pwork`, `Qwork`, `betawork`, `dirtd3`: transform workspaces.
- `augmentation`: the direct current block, or `nothing`.
"""
struct HBNonlinearProblem{S,ML,J,X,B,TP,FD,TD,A}
    sys::S
    modelayout::ML
    u0::Vector{Float64}
    jacobian::J
    parts::X          # the pieces a preconditioner needs
    bnm0::B           # the drive as built, for setdrive!
    tplan::TP         # transposed gather maps for the vjp
    Pwork::FD; Qwork::FD; betawork::TD; dirtd3::TD
    augmentation::A   # the direct current block, or `nothing`
end

function HBNonlinearProblem(sys, ml, u0, J, parts; augmentation = nothing)
    fd = sys.phimatrix; td = sys.phitd
    # The transpose plan is built here rather than on first use so the field
    # is concretely typed. A `Ref{Any}` filled lazily makes every call into
    # the transposed product a dynamic dispatch, which allocates on each
    # Jacobian application in the innermost Krylov loop; the build itself is
    # a counting sort over the plan's own index arrays and costs nothing.
    tp = plannonlineartermtranspose(sys.nonlineartermplan, ml,
        fd, td; backend = sys.nonlineartermplan.backend)
    return HBNonlinearProblem(sys, ml, u0, J, parts, copy(sys.bnm),
        tp, Ref(similar(fd)), Ref(similar(fd)), Ref(similar(td)),
        Ref(similar(td)), augmentation)
end

"""
    isaugmented(p::HBNonlinearProblem)

Whether `p` carries an explicit direct current block, so that its unknowns
are the canonical state rather than the harmonic one.
"""
isaugmented(p::HBNonlinearProblem) = !isnothing(p.augmentation)

# the harmonic residual and its derivatives, evaluated at the internal point
# a canonical one scatters to
function _setcanonical!(p::HBNonlinearProblem, u::AbstractVector)
    a = p.augmentation
    scattercanonical!(a.work.xint, u, a.work.layout)
    setpoint!(p.sys, a.work.xint)
    return a
end

Base.length(p::HBNonlinearProblem) = length(p.u0)

"""
    hbnonlinearproblem(w, Nharmonics, sources, circuit, circuitdefs; kwargs...)

Build the harmonic balance system as an [`HBNonlinearProblem`](@ref)
without solving it.
"""
function hbnonlinearproblem(w, Nharmonics, sources, circuit, circuitdefs;
        assemblejacobian::Bool = true, kwargs...)
    d = hbnlsolve(w, Nharmonics, sources, circuit, circuitdefs;
        returnsystem = true, assemblejacobian = assemblejacobian, kwargs...)
    d.dcexplicit || return HBNonlinearProblem(d.sys, d.modelayout,
        copy(d.xr), d.Jr, d)
    # with an explicit direct current block the unknowns are the canonical
    # state and the Jacobian is the canonical one, so the problem handed out
    # is the system which was solved rather than the harmonic part of it
    a = DCAugmentation(d.canonicalwork, d.Jr)
    L = a.work.layout
    u0 = zeros(Float64, canonicaldim(L))
    gathercanonical!(u0, d.xr, L)
    return HBNonlinearProblem(d.sys, d.modelayout, u0,
        isnothing(a.jplan) ? nothing : a.jplan.J, d; augmentation = a)
end

"""
    hbresidual!(F, prob::HBNonlinearProblem, u)

The harmonic balance residual at `u`, in place.
"""
function hbresidual!(F::AbstractVector{<:Real}, p::HBNonlinearProblem,
        u::AbstractVector{<:Real})
    if !isaugmented(p)
        setpoint!(p.sys, u)
        residual!(F, p.sys)
        return F
    end
    a = _setcanonical!(p, u)
    w = a.work
    residual!(w.Fint, p.sys)
    gathercanonical!(F, w.Fint, w.layout)
    # the block's rows, which replace where the harmonic system has nothing
    # and add where it has something. `residual = false` leaves out the
    # drive dependent constant, which is added back at the current scale so
    # that `setdrive!` moves the direct current injection with the rest
    addtransport!(F, w, u; residual = false)
    @. F += a.scale[]*a.constant
    return F
end

"""
    hbjvp!(Jv, prob::HBNonlinearProblem, u, v)

The exact matrix-free Jacobian-vector product `J(u)*v`. Two transforms plus
the linear term, with no Jacobian assembled.

This convenience form re-sets the evaluation point on every call, which
costs one extra forward transform. Inside a Krylov loop use a
[`JacobianOperator`](@ref), which sets the point once at construction and
whose `mul!` pays only the product itself.
"""
function hbjvp!(Jv::AbstractVector{<:Real}, p::HBNonlinearProblem,
        u::AbstractVector{<:Real}, v::AbstractVector{<:Real})
    if !isaugmented(p)
        setpoint!(p.sys, u)
        jacobianvectorproduct!(Jv, p.sys, v)
        return Jv
    end
    a = _setcanonical!(p, u)
    return _canonicaljvp!(Jv, p, v)
end

# the product at the point already set, which is what a Krylov loop wants
function _canonicaljvp!(Jv::AbstractVector, p::HBNonlinearProblem,
        v::AbstractVector)
    a = p.augmentation
    w = a.work
    L = w.layout
    scattercanonical!(w.xint, v, L)
    jacobianvectorproduct!(w.Fint, p.sys, w.xint)
    gathercanonical!(Jv, w.Fint, L)
    addtransport!(Jv, w, v; residual = false)
    return Jv
end

"""
    hbjacobian!(J, prob::HBNonlinearProblem, u)

Assemble the exact real Jacobian at `u`.
"""
function hbjacobian!(J, p::HBNonlinearProblem, u::AbstractVector{<:Real})
    if !isaugmented(p)
        setpoint!(p.sys, u)
        jacobian!(J, p.sys)
        return J
    end
    a = _setcanonical!(p, u)
    isnothing(a.jplan) && throw(ArgumentError("this problem assembled no Jacobian; build it with assemblejacobian = true"))
    jacobian!(a.jint, p.sys)
    canonicaljacobian!(a.jplan, a.jint)
    J === a.jplan.J || copyto!(J.nzval, a.jplan.J.nzval)
    return J
end

# =====================================================================
# The Jacobian and the preconditioner as operators.
#
# `mul!` and `ldiv!` are the closest thing to a universal interface for
# linear algebra in Julia: Krylov.jl, IterativeSolvers.jl, KrylovKit.jl,
# LinearSolve.jl, LinearMaps.jl and BifurcationKit all consume them, and
# they cost no dependency because they live in `Base`/`LinearAlgebra`.
#
# Handing out these objects, rather than the closures they wrap, is what
# lets an external solver be used without a package extension: the user
# writes `Krylov.gmres(J, -F; N = P)` in their own script.
# =====================================================================

"""
    JacobianOperator(prob::HBNonlinearProblem, u)

The Jacobian of `prob` at `u` as a linear operator.

Implements `size`, `eltype`, `LinearAlgebra.mul!`, and the same for its
`adjoint` and `transpose` (which are matrix free through the transposed
gather maps; see [`hbvjp!`](@ref)). Construction sets the system's
evaluation point ONCE -- one forward transform -- and `mul!` never touches
it again, so the products inside a Krylov loop pay exactly two transforms
each and nothing more. After moving `u`, construct a new operator: that is
the point update, and it is the same cost the internal solver pays once
per Newton step.

```julia
J = JacobianOperator(prob, u)
Krylov.gmres(J, -F; N = preconditioner(prob, u), atol = 0.0)
```

!!! warning "Set the absolute tolerance to zero inside a Newton loop"
    `rtol` in a Krylov solver is relative to the norm of the right hand
    side, and in a Newton loop that right hand side is the residual being
    driven to zero. Any absolute floor eventually exceeds it, at which
    point the linear solver correctly reports success for a system it never
    touched, the Newton step is zero, and the iteration stagnates with
    nothing reporting a failure.

    Krylov.jl defaults to `atol = sqrt(eps())`, about 1.5e-8, which is
    sensible standalone and wrong here. This package's own `gmres!`
    defaults to `atol = 0.0`, and `nlsolvekrylov!` passes `ftol/10`, tying
    the floor to the nonlinear tolerance rather than to machine epsilon.
"""
struct JacobianOperator{P<:HBNonlinearProblem,U<:AbstractVector}
    prob::P
    u::U
end
function JacobianOperator(prob::HBNonlinearProblem, u)
    isaugmented(prob) ? _setcanonical!(prob, u) : setpoint!(prob.sys, u)
    return JacobianOperator{typeof(prob),typeof(u)}(prob, u)
end

Base.size(J::JacobianOperator) = (length(J.u), length(J.u))
Base.size(J::JacobianOperator, i::Integer) = size(J)[i]
Base.eltype(::JacobianOperator) = Float64
Base.axes(J::JacobianOperator) = (Base.OneTo(size(J,1)), Base.OneTo(size(J,2)))
Base.axes(J::JacobianOperator, i::Integer) = axes(J)[i]

# the point was set at construction; the product is two transforms plus the
# linear term, with no setpoint on the hot path
LinearAlgebra.mul!(y::AbstractVector, J::JacobianOperator, v::AbstractVector) =
    isaugmented(J.prob) ? _canonicaljvp!(y, J.prob, v) :
        jacobianvectorproduct!(y, J.prob.sys, v)

function LinearAlgebra.mul!(y::AbstractVector, J::JacobianOperator,
        v::AbstractVector, alpha::Number, beta::Number)
    tmp = mul!(similar(y), J, v)
    @. y = alpha*tmp + beta*y
    return y
end

# `JacobianOperator` is deliberately not an `AbstractMatrix` -- it has no
# `getindex` -- so the lazy wrappers are constructed explicitly.
Base.adjoint(J::JacobianOperator) = LinearAlgebra.Adjoint(J)
Base.transpose(J::JacobianOperator) = LinearAlgebra.Transpose(J)
# `AbstractVector` is too wide: AxisKeys and NamedDims define `*` for a
# generic `AbstractMatrix` against their own vector wrappers, and the two
# signatures are mutually unreachable. The operator only ever multiplies
# plain strided vectors, so narrowing removes the ambiguity without losing
# anything.
Base.:*(J::JacobianOperator, v::StridedVector) =
    mul!(similar(v, size(J,1)), J, v)

# The transposed and adjoint products are matrix free through the
# transposed gather maps -- no assembled Jacobian is required. The output
# eltype is narrowed to `Float64` because the `Adjoint`/`Transpose`
# wrappers are `AbstractMatrix`, and a generic `y::AbstractVector` is
# ambiguous against packages which define `mul!` for their own element
# types over all of `AbstractVecOrMat` (MutableArithmetics does).
for W in (:(LinearAlgebra.Adjoint), :(LinearAlgebra.Transpose))
    @eval begin
        function LinearAlgebra.mul!(y::AbstractVector{Float64},
                Jt::$W{<:Any,<:JacobianOperator}, w::AbstractVector)
            J = parent(Jt)
            return hbvjp!(y, J.prob, J.u, w)
        end
        Base.:*(Jt::$W{<:Any,<:JacobianOperator}, w::StridedVector) =
            mul!(similar(w, Float64, size(Jt,1)), Jt, w)
    end
end

"""
    jacobianprototype(prob::HBNonlinearProblem)

A copy of the assembled real Jacobian's sparsity pattern, for solvers which
want a `jac_prototype`, or `nothing`.
"""
jacobianprototype(p::HBNonlinearProblem) =
    isnothing(p.jacobian) ? nothing : copy(p.jacobian)


"""
    preconditioner(prob::HBNonlinearProblem, u; couplingmodes = :none,
        factorization = KLUfactorization())

The mode coupling preconditioner of `prob`, updated at `u`.

Applied by `ldiv!` and by `mul!`, so it can be handed straight to any
external Krylov solver: `Krylov.gmres(J, -F; N = preconditioner(prob, u),
atol = 0.0)`. See [`JacobianOperator`](@ref) for why `atol = 0.0`.
"""
function preconditioner(p::HBNonlinearProblem, u::AbstractVector;
        couplingmodes = :none, factorization = KLUfactorization(),
        precision = Float64)
    d = p.parts
    pc = ModeCouplingPreconditioner(p.sys, d.Amatrixindicesaliased,
        d.Amatrixconjindices, d.Ljb, d.Lmean, d.Rbnm, d.Nmodes, d.Nbranches,
        d.Nfreq, d.invLnm, d.Gnm, d.Cnm, p.modelayout;
        couplingmodes = couplingmodes, factorization = factorization,
        precision = precision, Amatrixmodes = d.Amatrixmodes)
    if isaugmented(p)
        # the same wrapper the internal solve uses: the mode coupling
        # preconditioner in the canonical coordinates, with the direct
        # current subsystem factorized and solved exactly
        cp = CanonicalPreconditioner(pc, p.augmentation.work)
        updatepreconditioner!(cp, u)
        return SizedPreconditioner(cp, length(u))
    end
    setpoint!(p.sys, u)
    updatepreconditioner!(pc, u)
    return SizedPreconditioner(pc, length(u))
end

"""
    SizedPreconditioner(pc, n)

A preconditioner carrying its dimension, so it satisfies the `size`/
`eltype`/`mul!` contract external Krylov solvers check.
"""
struct SizedPreconditioner{P} <: AbstractPreconditioner
    pc::P
    n::Int
end
Base.size(p::SizedPreconditioner) = (p.n, p.n)
Base.size(p::SizedPreconditioner, i::Integer) = p.n
Base.eltype(::SizedPreconditioner) = Float64

# subtyping `AbstractPreconditioner` inherits the whole interface: the two
# and three argument `ldiv!`, `\\` and `mul!`. IterativeSolvers.jl calls the
# two argument in-place form on a view, KrylovKit and Krylov.jl call `mul!`,
# LinearSolve calls the three argument `ldiv!`; all of them now work.
applypreconditioner!(z::AbstractVector, p::SizedPreconditioner,
    r::AbstractVector) = applypreconditioner!(z, p.pc, r)
updatepreconditioner!(p::SizedPreconditioner, u::AbstractVector) =
    (updatepreconditioner!(p.pc, u); p)

# =====================================================================
# The nonlinear solver behind the operating point.
#
# `hbnlsolve` fuses three phases: build the system, solve it, and convert
# the result back. Only the middle one is replaceable, and only it has
# nothing harmonic balance specific about it: given a residual, a
# Jacobian-vector product and a preconditioner, find the root.
#
# The solver is a value rather than a symbol so it carries its own options,
# and so an external one can be passed in without the package naming it.
# =====================================================================

"""
    AbstractHBNonlinearSolver

The solver used for the harmonic balance operating point. See
[`NewtonKrylov`](@ref), [`Newton`](@ref), [`QuasiNewton`](@ref) and
[`ExternalSolver`](@ref).
"""
abstract type AbstractHBNonlinearSolver end

"""
    NewtonKrylov(; kwargs...)

Jacobian-free Newton-Krylov with the mode coupling preconditioner: the
default, and what `method = :newtonkrylov` selects. Keywords are the
`krylov*` options of [`hbnlsolve`](@ref), including `linearsolver`.
"""
struct NewtonKrylov{K} <: AbstractHBNonlinearSolver
    kwargs::K
end
NewtonKrylov(; kwargs...) = NewtonKrylov(kwargs)

"""
    Newton()

Newton's method on the equivalent real system with the exact assembled real
Jacobian. `method = :newton`.
"""
struct Newton <: AbstractHBNonlinearSolver end

"""
    QuasiNewton(; andersondepth = 5)

The holomorphic Jacobian approximation with Anderson acceleration.
`method = :quasinewton`. The harmonic balance residual is not complex
differentiable, so this Jacobian is an approximation.
"""
struct QuasiNewton <: AbstractHBNonlinearSolver
    andersondepth::Int
end
QuasiNewton(; andersondepth::Integer = 5) = QuasiNewton(andersondepth)

"""
    ExternalSolver(f)

Solve the operating point with a caller supplied root finder.

`f(prob, u0)` receives an [`HBNonlinearProblem`](@ref) and the initial
value in the real representation, and returns `(u, converged)`. Everything
it needs is on `prob`: [`hbresidual!`](@ref), [`hbjvp!`](@ref),
[`JacobianOperator`](@ref) and [`preconditioner`](@ref).

This is the plug point for a solver the package does not know about. A
NonlinearSolve.jl algorithm, a hand written continuation stepper or a
homotopy all go here without an extension.

The assembled real Jacobian is available on the problem unless
`assemblejacobian = false` was passed, which is what a matrix-free solver
wants: on a multi-tone problem that plan is the largest object in the
solve.

```julia
ExternalSolver() do prob, u0
    u = copy(u0); F = similar(u)
    hbresidual!(F, prob, u)
    for k in 1:40
        J = JacobianOperator(prob, u)
        P = preconditioner(prob, u)
        d, st = Krylov.gmres(J, -F; N = P, rtol = 1e-10, atol = 0.0)
        st.solved || return (u, false)
        u .+= d; hbresidual!(F, prob, u)
    end
    return (u, norm(F) <= tol)
end
```

!!! warning "Absolute tolerances stall Newton"
    Note `atol = 0.0`. Krylov.jl defaults to `atol = sqrt(eps())`, about
    1.5e-8, and stops as soon as the linear residual falls below it. Once
    the Newton residual is smaller than that -- which is the whole point of
    the last few Newton steps -- every linear solve returns immediately
    having done zero iterations, reports success, and hands back a zero
    step. Newton then stagnates while nothing reports a failure.

    Measured on a JPA with the default `atol`: 40 Newton iterations, final
    residual 3.2e-10, never converged. With `atol = 0.0`: 7 Newton
    iterations, residual 7.4e-17. Any external Krylov solver used inside a
    Newton loop wants its absolute tolerance set to zero and its stopping
    left to the relative one.
"""
struct ExternalSolver{F} <: AbstractHBNonlinearSolver
    f::F
end
ExternalSolver(f::Function) = ExternalSolver{typeof(f)}(f)

"""
    solvermethod(solver)

The `method` symbol of the built-in solvers, for the dispatch inside
[`hbnlsolve`](@ref).
"""
solvermethod(::NewtonKrylov) = :newtonkrylov
solvermethod(::Newton) = :newton
solvermethod(::QuasiNewton) = :quasinewton
solvermethod(::ExternalSolver) = :external
solvermethod(m::Symbol) = m

solverkwargs(s::NewtonKrylov) = s.kwargs
solverkwargs(::AbstractHBNonlinearSolver) = (;)
solverkwargs(::Symbol) = (;)

"""
    setdrive!(prob::HBNonlinearProblem, scale)

Scale the drive of `prob` by `scale`, in place, and return the problem.

The residual is `F(u) = B(sin(A*u)) + K*u - b`, and the drive enters only
through `b`. Scaling it is therefore the one parameter which can be varied
without rebuilding anything: the sparsity, the plans and the preconditioner
structure are all untouched.

This is what makes continuation in pump power possible. `scale = 0` is the
undriven problem, whose solution is zero and whose Jacobian is the linear
circuit; stepping up to `scale = 1` and carrying the converged state
forward walks onto the driven branch, which is the reliable way to reach an
operating point a cold solve cannot find.

The scale is relative to the drive the problem was built with, so
`setdrive!(prob, 1)` restores it.
"""
function setdrive!(p::HBNonlinearProblem, scale::Real)
    isnothing(p.bnm0) && throw(ArgumentError(
        "this problem did not record its drive; build it with "*
        "hbnonlinearproblem"))
    @. p.sys.bnm = scale*p.bnm0
    complex_to_real!(p.sys.bnmr, p.sys.bnm, p.modelayout.isreal)
    # the direct current injection is part of the same drive and moves with
    # it; the block's constant is the injection, so scaling it is enough
    isaugmented(p) && (p.augmentation.scale[] = scale)
    return p
end

"""
    drivenresidual!(F, prob, u, scale)

The residual at `u` with the drive scaled by `scale`: the `F(u, p)` form a
continuation library asks for, with `p` the drive scale.
"""
function drivenresidual!(F::AbstractVector{<:Real}, p::HBNonlinearProblem,
        u::AbstractVector{<:Real}, scale::Real)
    setdrive!(p, scale)
    return hbresidual!(F, p, u)
end

# =====================================================================
# The remaining ingredients for matrix-free continuation.
#
# A continuation library asks for more than the residual and the Jacobian
# product. Codim-2 continuation and the minimally augmented formulations
# need the transposed product; normal forms need the second and third
# directional derivatives; the tangent predictor needs the derivative with
# respect to the continuation parameter.
#
# All of them are available without assembling anything. The second and
# third derivatives in particular are nearly free: the plan's forward and
# backward maps are shared, and only the pointwise time domain function
# between the two transforms changes.
# =====================================================================

"""
    transposeplan!(prob::HBNonlinearProblem)

The [`NonlinearTermTransposePlan`](@ref) of `prob`.
"""
transposeplan!(p::HBNonlinearProblem) = p.tplan

"""
    hbvjp!(out, prob, u, w)

The transposed product `transpose(J(u))*w`, matrix free.

Needed by BifurcationKit for minimally augmented codim-2 continuation, fold
and Hopf continuation and left eigenvectors, and accepted by NonlinearSolve
as `vjp`.

The adjoints of the two transforms compose so that the normalizations
cancel, leaving the conjugate multiplicity divided out in one kernel and
multiplied back in through the coefficients of the other; see
[`NonlinearTermTransposePlan`](@ref).
"""
function hbvjp!(out::AbstractVector{<:Real}, p::HBNonlinearProblem,
        u::AbstractVector{<:Real}, w::AbstractVector{<:Real})
    if !isaugmented(p)
        setpoint!(p.sys, u)
        _ensurecos!(p.sys)
        # the transpose plan and the work buffers are concretely typed
        # fields built at construction, so nothing here dispatches
        # dynamically or allocates
        return _hbvjp!(out, p.sys, p.tplan, p.Pwork[], p.Qwork[],
            p.betawork[], w)
    end
    a = _setcanonical!(p, u)
    _ensurecos!(p.sys)
    L = a.work.layout
    # The canonical Jacobian is `D G J S + M`, so its transpose is
    # `S' J' G' D + M'`: mask the direction by `D` first, take the harmonic
    # transposed product through the scatter and the gather, then add the
    # block's own transpose. `M` lives entirely in the window, and so does
    # its transpose.
    z = a.zwork
    @. z = a.keep*w
    scattercanonical!(a.work.xint, z, L)
    _hbvjp!(a.Fwork, p.sys, p.tplan, p.Pwork[], p.Qwork[], p.betawork[],
        a.work.xint)
    fill!(out, 0.0)
    gathercanonical!(out, a.Fwork, L)
    win = dcwindow(L)
    mul!(view(out, win), a.dcmatrix, view(w, win), 1.0, 1.0)
    return out
end

function _hbvjp!(out, sys, tp, P, Q, beta, w)
    applybackwardjosephsontranspose!(P, tp, sys.nonlineartermplan, w)
    if !isnothing(sys.rfftplan)
        applyifft!(beta, P, sys.irfftplan)
        beta .*= sys.costd
        applyfft!(Q, beta, sys.rfftplan)
    else
        fill!(Q, 0)
    end
    applyforwardtranspose!(out, tp, Q, w)
    return out
end

"""
    hbd2F!(out, prob, u, v, w)

The exact second directional derivative
`H(u)[v,w] = B(-sin(A*u).*(A*v).*(A*w))`.
"""
function hbd2F!(out::AbstractVector{<:Real}, p::HBNonlinearProblem,
        u::AbstractVector{<:Real}, v::AbstractVector{<:Real},
        w::AbstractVector{<:Real})
    if !isaugmented(p)
        setpoint!(p.sys, u)
        hessianvectorproduct!(out, p.sys, v, w)
        return out
    end
    # the block is affine, so it contributes nothing beyond the first
    # derivative; what remains is the harmonic second derivative in the
    # canonical coordinates, with the replaced rows masked out
    a = _setcanonical!(p, u)
    L = a.work.layout
    scattercanonical!(a.work.xint, v, L)
    scattercanonical!(a.Fwork, w, L)
    hessianvectorproduct!(a.work.Fint, p.sys, a.work.xint, a.Fwork)
    fill!(out, 0.0)
    gathercanonical!(out, a.work.Fint, L)
    @. out *= a.keep
    return out
end

"""
    hbd3F!(out, prob, u, v, w, z)

The exact third directional derivative
`d3F(u)[v,w,z] = B(-cos(A*u).*(A*v).*(A*w).*(A*z))`.

Supplying this exactly rather than by finite differences is what makes a
continuation library's normal form computation -- cusp, Bogdanov-Takens,
Bautin -- accurate.
"""
function hbd3F!(out::AbstractVector{<:Real}, p::HBNonlinearProblem,
        u::AbstractVector{<:Real}, v::AbstractVector{<:Real},
        w::AbstractVector{<:Real}, z::AbstractVector{<:Real})
    if isaugmented(p)
        # affine again, so only the harmonic term survives; it is taken in
        # the internal coordinates and carried back
        a = _setcanonical!(p, u)
        L = a.work.layout
        vi = similar(a.Fwork); wi = similar(a.Fwork); zi = similar(a.Fwork)
        scattercanonical!(vi, v, L)
        scattercanonical!(wi, w, L)
        scattercanonical!(zi, z, L)
        _hbd3F!(a.work.Fint, p, vi, wi, zi)
        fill!(out, 0.0)
        gathercanonical!(out, a.work.Fint, L)
        @. out *= a.keep
        return out
    end
    setpoint!(p.sys, u)
    return _hbd3F!(out, p, v, w, z)
end

function _hbd3F!(out::AbstractVector{<:Real}, p::HBNonlinearProblem,
        v::AbstractVector{<:Real}, w::AbstractVector{<:Real},
        z::AbstractVector{<:Real})
    sys = p.sys
    _ensurecos!(sys)
    plan = sys.nonlineartermplan
    applyforwardterm!(sys.phimatrix, plan, v)
    applyifft!(sys.dirtd, sys.phimatrix, sys.irfftplan)
    applyforwardterm!(sys.phimatrix, plan, w)
    applyifft!(sys.dirtd2, sys.phimatrix, sys.irfftplan)
    applyforwardterm!(sys.phimatrix, plan, z)
    applyifft!(p.dirtd3[], sys.phimatrix, sys.irfftplan)
    sys.worktd .= .-sys.costd .* sys.dirtd .* sys.dirtd2 .* p.dirtd3[]
    applyfft!(sys.phimatrix, sys.worktd, sys.rfftplan)
    applybackwardterm!(out, plan, sys.phimatrix, v; addlinearterm = false)
    return out
end

"""
    hbdFdp!(out, prob)

The derivative of the residual with respect to the drive scale, which is
`-b`: the residual is `B(sin(A*u)) + K*u - scale*b`, so the parameter
derivative is exact and constant, and a continuation tangent never needs a
finite difference in the parameter.
"""
function hbdFdp!(out::AbstractVector{<:Real}, p::HBNonlinearProblem)
    isnothing(p.bnm0) && throw(ArgumentError("this problem recorded no drive"))
    b = similar(p.sys.bnm)
    copyto!(b, p.bnm0)
    if !isaugmented(p)
        br = zeros(Float64, length(out))
        complex_to_real!(br, b, p.modelayout.isreal)
        @. out = -br
        return out
    end
    # the alternating current drive through the gather, with the replaced
    # rows masked out, plus the block's own constant, which is the direct
    # current injection and carries its own sign
    a = p.augmentation
    L = a.work.layout
    br = zeros(Float64, L.rdim)
    complex_to_real!(br, b, p.modelayout.isreal)
    @. br = -br
    fill!(out, 0.0)
    gathercanonical!(out, br, L)
    @. out = a.keep*out + a.constant
    return out
end
