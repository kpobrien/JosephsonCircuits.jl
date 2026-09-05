# The vocabulary of the solvers: the objects a caller composes to say how
# an operating point is solved. Every option is a keyword of the object it
# configures and of nothing else; the objects validate their own keywords.
# A method holds a preconditioner, a linear solver and a refresh policy; a
# preconditioner holds the factorization it is built with; a deflation is a
# preconditioner wrapping another. The runtime objects these specify
# (`ModeCouplingPreconditioner`, `RecyclingPreconditioner`, ...) keep their
# suffixed names and are built by the solver from these values.

# ---------------------------------------------------------------- factorizations

"""
    AbstractFactorization

A sparse factorization method. A concrete method `f` implements
`factorize(f, A)`, which computes a factorization of `A` (the symbolic
analysis and the first numeric factorization), and `refactorize!(f, F, A)`,
which refactorizes `A` into the existing `F` reusing its symbolic analysis
and returns `F`, or `nothing` when the method has no in place
refactorization. `batchedform(f, layout)` is the form of `f` which treats a
matrix as the uniform batch of diagonal blocks described by a
[`BatchedBlockLayout`](@ref), or `f` itself.

The methods are [`KLUfactorization`](@ref), [`LUfactorization`](@ref),
[`QRfactorization`](@ref), [`CUDSSFactorization`](@ref) and
[`BlockFactorization`](@ref).
"""
abstract type AbstractFactorization end

refactorize!(::AbstractFactorization, F, A) = nothing
batchedform(f::AbstractFactorization, layout) = f

"""
    KLUfactorization(; kwargs...)

The [`AbstractFactorization`](@ref) using KLU.jl, a sparse LU factorization
suited to circuit matrices. This is the default on the host. `kwargs` are
passed to `KLU.klu`. The fill reducing ordering is chosen by
[`kluordered`](@ref) rather than left at KLU's default.
"""
struct KLUfactorization{K} <: AbstractFactorization
    kwargs::K
end
KLUfactorization(; kwargs...) = KLUfactorization(kwargs)

"""
    LUfactorization(; kwargs...)

The [`AbstractFactorization`](@ref) using the UMFPACK sparse LU
factorization `LinearAlgebra.lu`, with `kwargs` passed to it.
"""
struct LUfactorization{K} <: AbstractFactorization
    kwargs::K
end
LUfactorization(; kwargs...) = LUfactorization(kwargs)

"""
    QRfactorization(; kwargs...)

The [`AbstractFactorization`](@ref) using the SPQR sparse QR factorization
`LinearAlgebra.qr`, with `kwargs` passed to it. QR does not support
refactorization in place, so each call factorizes from scratch.
"""
struct QRfactorization{K} <: AbstractFactorization
    kwargs::K
end
QRfactorization(; kwargs...) = QRfactorization(kwargs)

# `CUDSSFactorization` lives in cudss.jl and `BlockFactorization` in
# blockfactorization.jl, beside their methods.

# ---------------------------------------------------------------- preconditioners

"""
    AbstractPreconditionerSpec

How the Newton-Krylov solver preconditions its linear solves: one of the
mode coupling family, [`BlockDiagonal`](@ref), [`FullJacobian`](@ref),
[`HarmonicBand`](@ref), [`MeasuredBand`](@ref), [`Clusters`](@ref),
[`CoupledModes`](@ref) and [`CouplingMask`](@ref), each built with a
factorization, [`Automatic`](@ref), which picks among them by the problem
and the memory, or a deflation wrapping
one of them, [`Recycling`](@ref) or [`Floquet`](@ref). The solver builds
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
struct BlockDiagonal{F<:MaybeFactorization} <: AbstractModeCoupling
    factorization::F
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
struct FullJacobian{F<:MaybeFactorization} <: AbstractModeCoupling
    factorization::F
end
FullJacobian(; factorization::MaybeFactorization = nothing) =
    FullJacobian(factorization)

"""
    HarmonicBand(p; factorization = nothing)

The couplings whose harmonic offset is within `p`, an `Integer` number of
offset shells or a per tone tuple of bounds; see [`modebandmask`](@ref).
Grown by one offset per tone on escalation.
"""
struct HarmonicBand{P,F<:MaybeFactorization} <: AbstractModeCoupling
    p::P
    factorization::F
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
struct MeasuredBand{F<:MaybeFactorization} <: AbstractModeCoupling
    tol::Float64
    budget::Float64
    factorization::F
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
struct Clusters{F<:MaybeFactorization} <: AbstractModeCoupling
    factorization::F
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
struct CoupledModes{F<:MaybeFactorization} <: AbstractModeCoupling
    indices::Vector{Int}
    factorization::F
end
CoupledModes(indices::AbstractVector{<:Integer}; factorization::MaybeFactorization = nothing) =
    CoupledModes(sort!(unique(Vector{Int}(indices))), factorization)

"""
    CouplingMask(mask; factorization = nothing)

The couplings selected by an `Nmodes` by `Nmodes` `Bool` matrix, the block
coupling column mode `m2` into row mode `m1` kept where `mask[m1, m2]`.
"""
struct CouplingMask{F<:MaybeFactorization} <: AbstractModeCoupling
    mask::Matrix{Bool}
    factorization::F
end
function CouplingMask(mask::AbstractMatrix{Bool}; factorization::MaybeFactorization = nothing)
    size(mask, 1) == size(mask, 2) || throw(ArgumentError(
        lazy"a coupling mask must be square, not $(size(mask))."))
    return CouplingMask(Matrix{Bool}(mask), factorization)
end

"""
    Recycling(inner = BlockDiagonal(); size = 20, harvest = 8, form = :adef1)

The preconditioner `inner` wrapped in a [`RecyclingPreconditioner`](@ref):
a deflation subspace of at most `size` vectors kept across Newton steps
and, through the reuse object of a cached sweep, across its points, with
`harvest` vectors taken from each linear solve; `form` is `:adef1` or
`:adef2`. This pays when the base leaves a deficiency of rank comparable
to `size`.
"""
struct Recycling{I<:AbstractPreconditionerSpec} <: AbstractPreconditionerSpec
    inner::I
    size::Int
    harvest::Int
    form::Symbol
end
function Recycling(inner::AbstractPreconditionerSpec = BlockDiagonal();
    size::Integer = 20, harvest::Integer = 8, form::Symbol = :adef1)
    size >= 1 || throw(ArgumentError(lazy"`size` = $(size) must be at least 1."))
    1 <= harvest <= size || throw(ArgumentError(
        lazy"`harvest` = $(harvest) must be between 1 and `size` = $(size)."))
    form in (:adef1, :adef2) || throw(ArgumentError(
        lazy"`form` = $(form) must be `:adef1` or `:adef2`."))
    inner isa AbstractModeCoupling || throw(ArgumentError(
        "a deflation wraps a mode coupling preconditioner, not another deflation."))
    return Recycling(inner, Int(size), Int(harvest), form)
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
struct Floquet{I<:AbstractPreconditionerSpec,R} <: AbstractPreconditionerSpec
    inner::I
    size::Int
    harvest::Int
    ritz::Int
    candidates::Int
    ranktol::R
    benefittol::Float64
    cycleharvest::Bool
end
function Floquet(inner::AbstractPreconditionerSpec = BlockDiagonal();
    size::Integer = 20, harvest::Integer = 4, ritz::Integer = 0,
    candidates::Integer = 3*size, ranktol::Union{Nothing,Real} = nothing,
    benefittol::Real = 1e-6, cycleharvest::Bool = true)
    size >= 1 || throw(ArgumentError(lazy"`size` = $(size) must be at least 1."))
    inner isa AbstractModeCoupling || throw(ArgumentError(
        "a deflation wraps a mode coupling preconditioner, not another deflation."))
    return Floquet(inner, Int(size), Int(harvest), Int(ritz), Int(candidates),
        ranktol, Float64(benefittol), cycleharvest)
end

"""
    withfactorization(s, f)

The preconditioner spec `s` with its factorization replaced by `f` where it
had none; a deflation applies this to what it wraps.
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
withfactorization(s::Recycling, f) =
    Recycling(withfactorization(s.inner, f), s.size, s.harvest, s.form)
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

Rebuild the preconditioner only when a linear solve fails its tolerance,
stagnates or produces no descent direction: a frozen preconditioner, for
when a deflation ([`Recycling`](@ref), [`Floquet`](@ref)) is to carry the
solve across the Newton path against a base built once.
"""
struct Never end

const AbstractRefresh = Union{Always,Probe,Never}

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
which fails to reach its tolerance to be grown to the full Jacobian (see
[`escalatepreconditioner!`](@ref)), within the memory the grown factors
are predicted to take.

The solve ends promptly when it cannot succeed and says why, in the
`reason` of its [`IterationInfo`](@ref): `:iterations` when the Newton
steps are spent; `:work` when the Arnoldi steps exceed `iterations`
restart lengths, so that a preconditioner which runs every linear solve to
its limit cannot turn the step budget into hours; `:linesearch` when no
sufficient decrease can be found, once with none at all or twice in a row
short of the Armijo condition; `:progress` when the residual history
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
struct NewtonKrylov{P<:AbstractPreconditionerSpec,L,R<:AbstractRefresh,T} <: AbstractHBNonlinearSolver
    preconditioner::P
    linearsolver::L
    refresh::R
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
struct Newton{F<:MaybeFactorization} <: AbstractHBNonlinearSolver
    factorization::F
end
Newton(; factorization::MaybeFactorization = nothing) = Newton(factorization)

"""
    QuasiNewton(; anderson = 5, factorization = nothing)

The holomorphic Jacobian approximation with Anderson acceleration of depth
`anderson` (the maximum number of previous iterates used for the
extrapolation; less than one disables it). The harmonic balance residual is
not complex differentiable, so this Jacobian is an approximation.
"""
struct QuasiNewton{F<:MaybeFactorization} <: AbstractHBNonlinearSolver
    anderson::Int
    factorization::F
end
QuasiNewton(; anderson::Integer = 5, factorization::MaybeFactorization = nothing) =
    QuasiNewton(Int(anderson), factorization)

"""
    Staged(; grids = nothing, s0 = 0.5, smin = 0.02, interiorftol = 1e-7,
        interioriterations = 60, inner = NewtonKrylov(),
        interiorescalation = false, maxattempts = 60, verbose = false)

Source continuation on an adaptively grown harmonic grid, with `inner`
solving every stage; see [`stagedhbnlsolve`](@ref) for the keywords.
`grids` is the ladder of retained harmonic caps, `nothing` for the default
ladder of the problem's `Nharmonics`.
"""
struct Staged{G,I<:AbstractHBNonlinearSolver} <: AbstractHBNonlinearSolver
    grids::G
    s0::Float64
    smin::Float64
    interiorftol::Float64
    interioriterations::Int
    inner::I
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
    return Staged(grids, Float64(s0), Float64(smin), Float64(interiorftol),
        Int(interioriterations), inner, interiorescalation, Int(maxattempts),
        verbose)
end

# the precision of the iteration a method runs in
solverprecision(m::NewtonKrylov) = m.precision
solverprecision(m::Staged) = solverprecision(m.inner)
solverprecision(::AbstractHBNonlinearSolver) = Float64

# a method with its escalation switched, for the interior stages of the
# staged solver; methods without escalation are returned unchanged
withescalation(m::NewtonKrylov, flag::Bool) = NewtonKrylov(m.preconditioner,
    m.linearsolver, m.refresh, flag, m.precision)
withescalation(m::AbstractHBNonlinearSolver, ::Bool) = m
