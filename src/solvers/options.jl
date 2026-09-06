# The vocabulary of the solvers: the objects a caller composes to say how
# an operating point is solved. Every option is a keyword of the object it
# configures and of nothing else; the objects validate their own keywords.
# A method holds a preconditioner, a linear solver and a refresh policy; a
# preconditioner holds the factorization it is built with; a deflation is a
# preconditioner wrapping another. The runtime objects these specify
# (`ModeCouplingPreconditioner`, `FloquetPreconditioner`, ...) keep their
# suffixed names and are built by the solver from these values.

# ---------------------------------------------------------------- factorizations

"""
    AbstractFactorization

A sparse factorization method. A concrete method `f` implements
`factorize(f, A)`, which computes a factorization of `A` (the symbolic
analysis and the first numeric factorization), and `refactorize!(f, F, A)`,
which refactorizes `A` into the existing `F` reusing its symbolic analysis
and returns `F`, or `nothing` when the method has no in place
refactorization.

The methods are [`KLUfactorization`](@ref), [`LUfactorization`](@ref),
[`QRfactorization`](@ref), [`CUDSSFactorization`](@ref) and
[`BlockFactorization`](@ref).
"""
abstract type AbstractFactorization end

refactorize!(::AbstractFactorization, F, A) = nothing

"""
    KLUfactorization(; kwargs...)

The [`AbstractFactorization`](@ref) using KLU.jl, a sparse LU factorization
suited to circuit matrices. This is the default on the host. `kwargs` are
passed to `KLU.klu`. The fill reducing ordering is chosen by
[`kluordered`](@ref) rather than left at KLU's default.
"""
struct KLUfactorization <: AbstractFactorization
    kwargs::NamedTuple
end
KLUfactorization(; kwargs...) = KLUfactorization(NamedTuple(kwargs))

"""
    LUfactorization(; kwargs...)

The [`AbstractFactorization`](@ref) using the UMFPACK sparse LU
factorization `LinearAlgebra.lu`, with `kwargs` passed to it.
"""
struct LUfactorization <: AbstractFactorization
    kwargs::NamedTuple
end
LUfactorization(; kwargs...) = LUfactorization(NamedTuple(kwargs))

"""
    QRfactorization(; kwargs...)

The [`AbstractFactorization`](@ref) using the SPQR sparse QR factorization
`LinearAlgebra.qr`, with `kwargs` passed to it. QR does not support
refactorization in place, so each call factorizes from scratch.
"""
struct QRfactorization <: AbstractFactorization
    kwargs::NamedTuple
end
QRfactorization(; kwargs...) = QRfactorization(NamedTuple(kwargs))

# `CUDSSFactorization` lives in solvers/cudss.jl and `BlockFactorization` in
# solvers/blockclusters.jl, beside their methods.
#
# None of the option objects below has a type parameter. They are
# configuration, read a few times per solve, and a parameter on the
# factorization or the preconditioner they carry would make every solver
# body that takes them a new specialization per way of configuring them.

# ---------------------------------------------------------------- preconditioners

"""
    AbstractPreconditionerSpec

How the Newton-Krylov solver preconditions its linear solves: one of the
mode coupling family, [`BlockDiagonal`](@ref), [`FullJacobian`](@ref),
[`HarmonicBand`](@ref), [`MeasuredBand`](@ref), [`Clusters`](@ref),
[`CoupledModes`](@ref) and [`CouplingMask`](@ref), each built with a
factorization, [`Automatic`](@ref), which picks among them by the problem
and the memory, or [`Floquet`](@ref), a deflation wrapping one of them. The
solver builds
the runtime preconditioner from the value.
"""
abstract type AbstractPreconditionerSpec end

"""
    AbstractModeCoupling

The mode coupling family of [`AbstractPreconditionerSpec`](@ref): the
Jacobian with its mode coupling restricted to a selected set and reduced to
the mode diagonal elsewhere, factorized; see
[`ModeCouplingPreconditioner`](@ref). Each member carries the
[`AbstractFactorization`](@ref) it is built with in its `factorization`
field, `nothing` for the backend's default (KLU on the host, cuDSS on a
device).
"""
abstract type AbstractModeCoupling <: AbstractPreconditionerSpec end

const MaybeFactorization = Union{Nothing,AbstractFactorization}

"""
    BlockDiagonal(; factorization = nothing)

The mode block diagonal, what [`Automatic`](@ref) picks for one tone: one
small independent factorization per mode, and no coupling. Cheap, and
sufficient for one tone and for weak drives; on a strongly pumped device
it stalls and is grown to the full Jacobian by escalation.
"""
struct BlockDiagonal <: AbstractModeCoupling
    factorization::MaybeFactorization
end
BlockDiagonal(; factorization::MaybeFactorization = nothing) =
    BlockDiagonal(factorization)

"""
    FullJacobian(; factorization = nothing)

Every mode coupling: the full Jacobian, an exact preconditioner and a
direct solve. With a [`BlockFactorization`](@ref) this is the dense block
factorization over the circuit graph, the fastest measured method on three
or more tones.
"""
struct FullJacobian <: AbstractModeCoupling
    factorization::MaybeFactorization
end
FullJacobian(; factorization::MaybeFactorization = nothing) =
    FullJacobian(factorization)

"""
    HarmonicBand(p; factorization = nothing)

The couplings whose harmonic offset is within `p`, an `Integer` number of
offset shells or a per tone tuple of bounds; see [`modebandmask`](@ref).
Grown by one offset per tone on escalation.
"""
struct HarmonicBand <: AbstractModeCoupling
    p::Union{Integer,Tuple{Vararg{Integer}}}
    factorization::MaybeFactorization
end
HarmonicBand(p; factorization::MaybeFactorization = nothing) = HarmonicBand(p, factorization)

"""
    MeasuredBand(; tol = 1e-2, budget = 0.25, factorization = nothing)

A [`HarmonicBand`](@ref) whose per tone width is measured from the Fourier
coefficients of `cos(phi(t))` at every point and widened when the drive
demands it, starting from the block diagonal; see
[`cosphibandwidths`](@ref) for `tol` and `budget`. For two strong tones
this is the setting that matters.
"""
struct MeasuredBand <: AbstractModeCoupling
    tol::Float64
    budget::Float64
    factorization::MaybeFactorization
end
function MeasuredBand(; tol::Real = 1e-2, budget::Real = 0.25,
    factorization::MaybeFactorization = nothing)
    0 < tol < 1 || throw(ArgumentError(lazy"`tol` = $(tol) must be in (0, 1)."))
    0 < budget <= 1 || throw(ArgumentError(
        lazy"`budget` = $(budget) must be in (0, 1]."))
    return MeasuredBand(Float64(tol), Float64(budget), factorization)
end

"""
    Clusters(; factorization = nothing)

Clusters of modes measured from the operator: the block Jacobi coupling
strengths are probed and modes merged in decreasing strength until the
couplings left between clusters are contractive; see
[`spectralclusters`](@ref). Probed at the first point and again whenever
the solver reports a slow linear solve; the clusters only grow within a
solve. With a [`BlockFactorization`](@ref) each cluster is one dense block
factorization over the circuit graph, which halves the memory of
[`FullJacobian`](@ref) on three tones.
"""
struct Clusters <: AbstractModeCoupling
    factorization::MaybeFactorization
end
Clusters(; factorization::MaybeFactorization = nothing) = Clusters(factorization)

"""
    Automatic()

The fastest robust preconditioner measured for the problem which fits in
memory, chosen when the preconditioner is built, by the number of tones and
the memory the factors would take:

- one tone: [`BlockDiagonal`](@ref) with the backend's sparse
  factorization, which the escalation to the full Jacobian backs up;
- two or more tones: [`FullJacobian`](@ref) with a single precision
  [`BlockFactorization`](@ref) when its factors, sized exactly from the
  symbolic analysis by [`blockfactorbytes`](@ref), take at most half the
  free memory of the backend ([`freememory`](@ref)), leaving the rest to
  the system, the Krylov basis and the products; otherwise
  [`MeasuredBand`](@ref) with the backend's sparse factorization.

Measured on a GPU, the full block factors solve three tones at 128
junctions in 4.3 s against 38.9 s for the measured band, at a quarter of
its memory, and match the band on two tones; on the mild 1024-junction
two-tone line they take 4.8 s against 7.9 s. The choice does not look at
the mixing order of the tones or the circuit's topology, so that it
generalizes; nothing about it is tuned beyond the memory margin.
"""
struct Automatic <: AbstractModeCoupling
    factorization::Nothing
end
Automatic() = Automatic(nothing)

"""
    CoupledModes(indices; factorization = nothing)

Exactly these modes coupled in full, the rest on the mode diagonal; see
[`modecouplingmask`](@ref).
"""
struct CoupledModes <: AbstractModeCoupling
    indices::Vector{Int}
    factorization::MaybeFactorization
end
CoupledModes(indices::AbstractVector{<:Integer}; factorization::MaybeFactorization = nothing) =
    CoupledModes(sort!(unique(Vector{Int}(indices))), factorization)

"""
    CouplingMask(mask; factorization = nothing)

The couplings selected by an `Nmodes` by `Nmodes` `Bool` matrix, the block
coupling column mode `m2` into row mode `m1` kept where `mask[m1, m2]`.
"""
struct CouplingMask <: AbstractModeCoupling
    mask::Matrix{Bool}
    factorization::MaybeFactorization
end
function CouplingMask(mask::AbstractMatrix{Bool}; factorization::MaybeFactorization = nothing)
    size(mask, 1) == size(mask, 2) || throw(ArgumentError(
        lazy"a coupling mask must be square, not $(size(mask))."))
    return CouplingMask(Matrix{Bool}(mask), factorization)
end

"""
    Floquet(inner = BlockDiagonal(); size = 20, harvest = 4, ritz = 0,
        candidates = 3*size, ranktol = nothing, benefittol = 1e-6,
        cycleharvest = true)

The preconditioner `inner` wrapped in a [`FloquetPreconditioner`](@ref):
the residual-image deflation with physical candidates. `harvest` is the
number of singular directions per harvest, `ritz` the harmonic Ritz
directions on top of it, `candidates` the size of the candidate bank,
`ranktol` the rank tolerance of the residual image (`nothing` for the
precision's default) and `benefittol` the predicted improvement below
which a candidate is not built in; `cycleharvest` harvests every GMRES
cycle rather than the last.
"""
struct Floquet <: AbstractPreconditionerSpec
    inner::AbstractPreconditionerSpec
    size::Int
    harvest::Int
    ritz::Int
    candidates::Int
    ranktol::Union{Nothing,Float64}
    benefittol::Float64
    cycleharvest::Bool
end
function Floquet(inner::AbstractPreconditionerSpec = BlockDiagonal();
    size::Integer = 20, harvest::Integer = 4, ritz::Integer = 0,
    candidates::Integer = 3*size, ranktol::Union{Nothing,Real} = nothing,
    benefittol::Real = 1e-6, cycleharvest::Bool = true)
    # the same checks the runtime constructor makes, so a bad value is
    # refused where it is written rather than deep inside the solve
    size >= 1 || throw(ArgumentError(lazy"`size` = $(size) must be at least 1."))
    harvest >= 0 || throw(ArgumentError(
        lazy"`harvest` = $(harvest) must be nonnegative."))
    ritz >= 0 || throw(ArgumentError(lazy"`ritz` = $(ritz) must be nonnegative."))
    harvest + ritz >= 1 || throw(ArgumentError(
        "`harvest` and `ritz` cannot both be zero; the harvest would produce no candidates."))
    candidates >= size || throw(ArgumentError(
        lazy"`candidates` = $(candidates) must be at least `size` = $(size)."))
    isnothing(ranktol) || ranktol > 0 || throw(ArgumentError(
        lazy"`ranktol` = $(ranktol) must be positive."))
    benefittol >= 0 || throw(ArgumentError(
        lazy"`benefittol` = $(benefittol) must be nonnegative."))
    inner isa AbstractModeCoupling || throw(ArgumentError(
        "a deflation wraps a mode coupling preconditioner, not another deflation."))
    return Floquet(inner, Int(size), Int(harvest), Int(ritz), Int(candidates),
        ranktol, Float64(benefittol), cycleharvest)
end

"""
    withfactorization(s, f)

The preconditioner spec `s` with its factorization replaced by `f` where it
had none; a deflation applies this to what it wraps. An [`Automatic`](@ref)
is returned unchanged: it carries no factorization, and the member it
resolves to takes the backend's default (`resolveautomatic`).
"""
withfactorization(s::BlockDiagonal, f) =
    isnothing(s.factorization) ? BlockDiagonal(f) : s
withfactorization(s::FullJacobian, f) =
    isnothing(s.factorization) ? FullJacobian(f) : s
withfactorization(s::HarmonicBand, f) =
    isnothing(s.factorization) ? HarmonicBand(s.p, f) : s
withfactorization(s::MeasuredBand, f) =
    isnothing(s.factorization) ? MeasuredBand(s.tol, s.budget, f) : s
withfactorization(s::Clusters, f) =
    isnothing(s.factorization) ? Clusters(f) : s
withfactorization(s::CoupledModes, f) =
    isnothing(s.factorization) ? CoupledModes(s.indices, f) : s
withfactorization(s::CouplingMask, f) =
    isnothing(s.factorization) ? CouplingMask(s.mask, f) : s
withfactorization(s::Automatic, f) = s
withfactorization(s::Floquet, f) = Floquet(withfactorization(s.inner, f),
    s.size, s.harvest, s.ritz, s.candidates, s.ranktol, s.benefittol,
    s.cycleharvest)

# ---------------------------------------------------------------- the refresh policy

"""
    Always()

Rebuild the preconditioner before every Newton step: the default, and
the policy whose solve path depends on nothing measured, so that the same
problem solves the same way every time. [`Probe`](@ref) measures whether
a rebuild pays and skips the ones which do not.
"""
struct Always end

"""
    Probe()

Decide each rebuild by measurement: one application of the stale
preconditioner to the residual and one product give the one-step
reduction `rho = |J P^-1 F - F|/|F|`; the same measurement on the fresh
preconditioner, with the Arnoldi count `k_fresh` of its solve, calibrates
the prediction `k = k_fresh log(rho_fresh)/log(rho)` of the stale solve's
Arnoldi count, and the rebuild is skipped when `k` steps at the measured
cost of a step are cheaper than the measured rebuild plus a fresh solve.
Everything is measured, so the rule adapts to the device and the
factorization; it pays when a rebuild is expensive next to a solve, as
with a [`BlockFactorization`](@ref) of three tones, where it saved a
fifth to a third of the time. A rebuild forced by a failed, stalled or
non-descent solve is never skipped.
Because the decision rests on measured times, the path a solve takes,
and the answer within the tolerance, can differ between two runs of the
same problem; the default [`Always`](@ref) is reproducible.
"""
struct Probe end

"""
    Never()

Rebuild the preconditioner only when it is forced: a linear solve which
makes progress but misses its tolerance, a direction which is not a
descent direction, a line search which finds no decrease, or a successful
escalation. A stagnated solve is not retried and does not rebuild; its
step is replaced by the preconditioner solve. The slow-solve report to the
preconditioner ([`stalled!`](@ref)) is off as well, so a [`Clusters`](@ref)
preconditioner never remeasures under this policy. A frozen
preconditioner, for when a deflation ([`Floquet`](@ref)) is to carry the
solve across the Newton path against a base built once.
"""
struct Never end

const AbstractRefresh = Union{Always,Probe,Never}

"""
    meritslope!(Jv, jvp, p, F, ϕ0, w)

The slope `real(F' J p)` of the merit function `ϕ = F'F/2` along the step
`p`, where `p = -Δ` for a linear solve `J Δ ≈ F` which left the explicit
residual `w = F - J Δ`.

Then `J p = w - F` and the slope is `real(F'w) - 2ϕ0`: one inner product,
with the Jacobian vector product the solve already paid for. Without a
valid residual (`w === nothing`), the product is taken.
"""
function meritslope!(Jv, jvp, p, F, ϕ0, w)
    if isnothing(w)
        mul!(Jv, jvp, p)
        return real(dot(F, Jv))
    end
    return real(dot(F, w)) - 2ϕ0
end

abstract type AbstractHBLinearSolver end

# ---------------------------------------------------------------- the methods

"""
    AbstractHBNonlinearSolver

A method of solving the operating point, the `method` of
[`hbnlsolve`](@ref) and [`hbsolve`](@ref): [`NewtonKrylov`](@ref),
[`Newton`](@ref), [`QuasiNewton`](@ref), [`Staged`](@ref) or
[`ExternalSolver`](@ref).
"""
abstract type AbstractHBNonlinearSolver end

"""
    NewtonKrylov(; preconditioner = Automatic(), linearsolver = GMRES(),
        refresh = Always(), escalate = true, precision = Float64)

Jacobian-free Newton-Krylov with the mode coupling preconditioner: the
default. `preconditioner` is an [`AbstractPreconditionerSpec`](@ref); the
default [`Automatic`](@ref) picks the fastest one measured that fits in
memory. `linearsolver` is a [`GMRES`](@ref) or a [`KrylovJL`](@ref)
solver, `refresh` [`Always`](@ref) (the default), [`Probe`](@ref) (which
rebuilds the preconditioner only when a measured probe says a rebuild
pays, and is faster by a fifth to a third on the hard cases, at the price
of a solve path which depends on measured times and so can differ
between two runs) or [`Never`](@ref). `escalate` allows a preconditioner
which fails to reach its tolerance to be grown (a band by one offset per
tone, any other set to the full Jacobian; see
[`escalatepreconditioner!`](@ref)), within the memory the grown factors
are predicted to take; a refused escalation is recorded and the solve
carries on.

The solve ends promptly when it cannot succeed and says why, in the
`reason` of its [`IterationInfo`](@ref): `:iterations` when the Newton
steps are spent; `:work` when the Arnoldi steps exceed `iterations`
restart lengths, so that a preconditioner which runs every linear solve to
its limit cannot turn the step budget into hours; `:linesearch` when no
sufficient decrease can be found (a step with no decrease at all is
retried once from a rebuilt preconditioner and ends the solve if it fails
again, as does a direction which is still not a descent direction after
the exact rescue, or two consecutive steps short of the Armijo
condition); `:progress` when the residual history
projects no convergence within the remaining budget and is not
accelerating, after one recovery which rebuilds the preconditioner and
takes exact Newton steps from then on ([`projectedstall`](@ref)). A stall
outside the Newton basin is the continuation problem [`Staged`](@ref)
exists for. `precision` is the floating point type
of the iteration: the system on the backend, the Krylov vectors, and the
factors of a sparse preconditioner; a single precision solve needs a
relative tolerance `rtol` it can meet.

The forcing sequence (Eisenstat-Walker choice 2 clamped to `[1e-10, 0.9]`,
starting at 0.3), the line search (Armijo with constant 1e-4, halving
with safeguards 0.1 and 0.5, ten trials, two failures) and the stagnation
threshold (a solve which does not bring the linear residual below 0.9 of
the residual norm) are fixed; see [`nlsolvekrylov!`](@ref).
"""
struct NewtonKrylov{T<:AbstractFloat} <: AbstractHBNonlinearSolver
    # the precision is the one parameter kept: it sets the element types of
    # the system's arrays, so a solve at each precision is a different
    # specialization whatever the option carries
    preconditioner::AbstractPreconditionerSpec
    linearsolver::AbstractHBLinearSolver
    refresh::AbstractRefresh
    escalate::Bool
    precision::Type{T}
end
function NewtonKrylov(; preconditioner::AbstractPreconditionerSpec = Automatic(),
    linearsolver = GMRES(), refresh::AbstractRefresh = Always(),
    escalate::Bool = true, precision::Type{<:AbstractFloat} = Float64)
    linearsolver isa AbstractHBLinearSolver || throw(ArgumentError(
        lazy"`linearsolver` = $(linearsolver) must be a `GMRES()` or a `KrylovJL` solver."))
    return NewtonKrylov(preconditioner, linearsolver, refresh, escalate,
        precision)
end

"""
    Newton(; factorization = nothing)

Newton's method on the equivalent real system with the exact assembled
real Jacobian, factorized by `factorization` (the host's KLU when
`nothing`).
"""
struct Newton <: AbstractHBNonlinearSolver
    factorization::MaybeFactorization
end
function Newton(; factorization::MaybeFactorization = nothing)
    checkdirectfactorization(factorization, "Newton")
    return Newton(factorization)
end

# a direct solve factorizes the assembled sparse Jacobian; the block
# factorization is a preconditioner's (through `NewtonKrylov`) or the
# linearized sweep's, and would fail deep inside the first factorization
# for want of a block size
function checkdirectfactorization(f, method)
    f isa BlockFactorization && throw(ArgumentError(
        lazy"$(method) factorizes the assembled sparse Jacobian and takes a sparse factorization (KLUfactorization(), LUfactorization(), CUDSSFactorization()); BlockFactorization() is the factorization of a NewtonKrylov preconditioner or of the linearized solve."))
    return f
end

"""
    QuasiNewton(; anderson = 5, factorization = nothing)

The holomorphic Jacobian approximation with Anderson acceleration of depth
`anderson` (the maximum number of previous iterates used for the
extrapolation; less than one disables it). The harmonic balance residual is
not complex differentiable, so this Jacobian is an approximation.

!!! warning "The zero frequency flux is complex"
    This method solves for a complex flux at every mode, the zero frequency
    mode included, so the imaginary part of a node's zero frequency flux is
    not held at zero and the converged value can carry a spurious imaginary
    part. The real part is the direct current flux; read only that, and
    prefer [`Newton`](@ref) or [`NewtonKrylov`](@ref), which solve the real
    system, for a circuit with a direct current bias.
"""
struct QuasiNewton <: AbstractHBNonlinearSolver
    anderson::Int
    factorization::MaybeFactorization
end
function QuasiNewton(; anderson::Integer = 5, factorization::MaybeFactorization = nothing)
    checkdirectfactorization(factorization, "QuasiNewton")
    return QuasiNewton(Int(anderson), factorization)
end

"""
    Staged(; grids = nothing, s0 = 0.5, smin = 0.02, interiorftol = 1e-7,
        interioriterations = 60, inner = NewtonKrylov(),
        interiorescalation = false, maxattempts = 60, verbose = false)

Source continuation on an adaptively grown harmonic grid, with `inner`
solving every stage; see [`stagedhbnlsolve`](@ref) for the keywords.
`grids` is the ladder of retained harmonic caps, `nothing` for the default
ladder of the problem's `Nharmonics`.
"""
struct Staged <: AbstractHBNonlinearSolver
    grids::Union{Nothing,AbstractVector}
    s0::Float64
    smin::Float64
    interiorftol::Float64
    interioriterations::Int
    inner::AbstractHBNonlinearSolver
    interiorescalation::Bool
    maxattempts::Int
    verbose::Bool
end
function Staged(; grids = nothing, s0::Real = 0.5, smin::Real = 0.02,
    interiorftol::Real = 1e-7, interioriterations::Integer = 60,
    inner::AbstractHBNonlinearSolver = NewtonKrylov(),
    interiorescalation::Bool = false, maxattempts::Integer = 60,
    verbose::Bool = false)
    inner isa Staged && throw(ArgumentError("`inner` must be a non-staged method."))
    0 < s0 <= 1 || throw(ArgumentError(lazy"`s0` = $(s0) must be in (0, 1]."))
    0 < smin <= s0 || throw(ArgumentError(
        lazy"`smin` = $(smin) must be in (0, `s0` = $(s0)]."))
    interiorftol > 0 || throw(ArgumentError(
        lazy"`interiorftol` = $(interiorftol) must be positive."))
    interioriterations >= 1 || throw(ArgumentError(
        lazy"`interioriterations` = $(interioriterations) must be at least 1."))
    maxattempts >= 1 || throw(ArgumentError(
        lazy"`maxattempts` = $(maxattempts) must be at least 1."))
    return Staged(grids, Float64(s0), Float64(smin), Float64(interiorftol),
        Int(interioriterations), inner, interiorescalation, Int(maxattempts),
        verbose)
end

"""
    solverprecision(m::AbstractHBNonlinearSolver)

The floating point type the method iterates in: `Float64` for every method
but [`NewtonKrylov`](@ref), whose `precision` it is, and the inner method's
for [`Staged`](@ref).
"""
solverprecision(m::NewtonKrylov) = m.precision
solverprecision(m::Staged) = solverprecision(m.inner)
solverprecision(::AbstractHBNonlinearSolver) = Float64

"""
    withescalation(m::AbstractHBNonlinearSolver, flag::Bool)

The method with its `escalate` set to `flag`, for the interior stages of
[`Staged`](@ref); methods without escalation are returned unchanged.
"""
withescalation(m::NewtonKrylov, flag::Bool) = NewtonKrylov(m.preconditioner,
    m.linearsolver, m.refresh, flag, m.precision)
withescalation(m::AbstractHBNonlinearSolver, ::Bool) = m

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

# === the canonical forms of the solver inputs ===
#
# The entry points accept their inputs in every way a user may write them
# (a frequency as an integer, a source as a named tuple of whatever number
# types, the definitions as a dictionary of any key and value types) and
# convert them once here, so that the solves below are compiled for one form
# of each rather than once per way of writing them: a solver body of a
# thousand lines compiled for `Tuple{Int64}` and again for `Tuple{Float64}`
# costs seconds each time and computes the same numbers.

"""
    SourceTuple{N}

The canonical form of a source: a named tuple `(mode, port, current)` with
the mode as `N` integers, one harmonic index per tone, the port as an
integer and the current as a complex number. See [`sourcetable`](@ref).
"""
const SourceTuple{N} = @NamedTuple{mode::NTuple{N,Int}, port::Int, current::ComplexF64}

"""
    tonefrequencies(w)

The tone frequencies `w` as a tuple of `Float64`. Each must be a real
number, in radians per second; the solve checks that they are finite.
"""
function tonefrequencies(w::NTuple{N,Number}) where {N}
    all(x -> x isa Real, w) || throw(ArgumentError(lazy"the tone frequencies $(w) must be real numbers, in radians per second."))
    return map(Float64, w)
end

"""
    sweepfrequencies(ws)

The signal frequencies `ws` as a `Vector{Float64}`: a real number, or any
iterable of real numbers (a vector, a range), in radians per second.
"""
function sweepfrequencies(ws)
    ws isa Real && return Float64[ws]
    ws isa Number && throw(ArgumentError(lazy"the signal frequency $(ws) must be a real number, in radians per second."))
    out = Float64[]
    for w in ws
        w isa Real || throw(ArgumentError(lazy"the signal frequency $(w) must be a real number, in radians per second."))
        push!(out, w)
    end
    return out
end

"""
    sourcetable(sources, w)

The sources as a `Vector{SourceTuple{N}}` for the `N` tones of `w`.
`sources` is any iterable of named tuples with the fields `mode`, `port`
and `current`; a source whose mode does not have one integer per tone,
whose port is not an integer or whose current is not a number is an
`ArgumentError`.
"""
function sourcetable(sources, ::NTuple{N,Number}) where {N}
    table = SourceTuple{N}[]
    for s in sources
        (s isa NamedTuple && hasfield(typeof(s), :mode) &&
            hasfield(typeof(s), :port) && hasfield(typeof(s), :current)) ||
            throw(ArgumentError(lazy"the source $(s) must be a named tuple with the fields mode, port and current, such as (mode = (1,), port = 1, current = 1e-6)."))
        (s.mode isa Tuple && length(s.mode) == N &&
            all(x -> x isa Integer, s.mode)) ||
            throw(ArgumentError(lazy"the source $(s) has the mode $(s.mode); a mode is a tuple of $(N) integers, one harmonic index per tone."))
        s.port isa Integer || throw(ArgumentError(lazy"the source $(s) has the port $(s.port), which must be an integer."))
        s.current isa Number || throw(ArgumentError(lazy"the source $(s) has the current $(s.current), which must be a number."))
        push!(table, (mode = NTuple{N,Int}(s.mode), port = Int(s.port),
            current = ComplexF64(s.current)))
    end
    return table
end

"""
    definitiontable(circuitdefs)

The component definitions as a `Dict{Any,Any}`, whatever the key and value
types of the dictionary given.
"""
definitiontable(circuitdefs::Dict{Any,Any}) = circuitdefs
definitiontable(circuitdefs::AbstractDict) = Dict{Any,Any}(circuitdefs)

"""
    initialguess(x0)

The initial guess `x0` of a nonlinear solve as a `Vector{ComplexF64}`,
empty when there is none (`nothing`); a node flux matrix of a previous
solve, keyed or plain, is flattened in the solver's own layout.
"""
initialguess(::Nothing) = ComplexF64[]
# a keyed node flux matrix of a previous solve is accepted as it is, in
# the layout the solver keeps it (modes fastest)
initialguess(x0::AbstractArray) =
    Vector{ComplexF64}(vec(AxisKeys.keyless_unname(x0)))

"""
    sensitivitypairtable(pairs)
    sensitivityblockpairtable(pairs)

The `(name, parameter, direction)` sensitivity pairs of a component or of a
scattering block as vectors of one tuple type, the direction of a
component pair as a complex number.
"""
sensitivitypairtable(pairs) =
    Tuple{String,Int,ComplexF64}[(String(t[1]), Int(t[2]), ComplexF64(t[3]))
        for t in pairs]
sensitivityblockpairtable(pairs) =
    Tuple{String,Int,Any}[(String(t[1]), Int(t[2]), t[3]) for t in pairs]
