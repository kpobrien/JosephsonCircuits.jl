# The same circuit under many drive conditions, stepped as one system:
# the states are matrices with the conditions as columns, every product
# and residual takes all columns in one call, the junction stiffness and
# the factorization are one per condition, KLU on the CPU and the uniform
# cuDSS batch on a device, and the Newton engine accepts and refreshes
# per column. The typical use of the transient is one circuit under many
# pumps or signals, so this is where a device earns its throughput: the
# launches of a step serve every condition. A single problem is a batch
# of one, on the same path.

"""
    transientproblem(problem::TransientProblem; sources)

The compiled circuit of `problem` under other `sources`, without
compiling again: the same matrices and augmentation, new bound drives.
The problems of one batch must drive the same targets in the same order,
which this gives when `sources` differ only in their waveforms.
"""
Base.@nospecializeinfer function transientproblem(p::TransientProblem; @nospecialize(sources))
    psc = p.circuit
    drives, injection, constantcurrent, balance = bindsources(psc, p.matrices.vvn, p.ports, p.portpositive, p.portnegative, length(p),
        p.floatingcomponents, sources)
    return TransientProblem(p; injection, drives, constantcurrent, balance)
end

# The problems of a batch share everything but their waveforms: the
# compiled circuit, the injection of the drives and the ports they drive,
# and the constant current of the sources no drive replaced, which is
# what the binding of the sources amounts to. A system depends on the
# problem through these alone, so a system built for one such problem
# serves the other, and a reuse keeps its system across them.
sharedcircuit(p::TransientProblem, q::TransientProblem) = q.circuit === p.circuit && q.matrices === p.matrices
function shareddrives(p::TransientProblem, q::TransientProblem)
    return q.injection == p.injection && length(q.drives) == length(p.drives) &&
        all(k -> q.drives[k].portindex == p.drives[k].portindex, eachindex(p.drives)) &&
        q.constantcurrent == p.constantcurrent
end
sharedsystem(p::TransientProblem, q::TransientProblem) = p === q || (sharedcircuit(p, q) && shareddrives(p, q))
function batchcompatible(problems)
    isempty(problems) && throw(ArgumentError("a batch needs at least one problem."))
    p = first(problems)
    for q in problems
        sharedcircuit(p, q) || throw(ArgumentError(
            "the problems of a batch must share one compiled circuit: build them with transientproblem(problem; sources)."))
        shareddrives(p, q) || throw(ArgumentError(
            "the problems of a batch must drive the same targets in the same order and leave the same sources constant; only the waveforms may differ."))
    end
    return p
end

"""
    TransientBatchSolution

Result of [`transientsolve`](@ref) on a vector of problems: the arrays of
a [`TransientSolution`](@ref) with the conditions as the trailing
dimension, `voltage[port, time, condition]`, `flux[state, time, condition]`,
`finalflux[state, condition]`, `phases[junction, stage, time, condition]`,
`stages[state, stage, time, condition]`; `problems` holds the batch. Indexing,
`solution[j]`, is the ordinary solution of condition `j`, a view of the
batch's arrays, on which the demodulation, the tangent, the adjoint and
the noise run as on any solution.

`stats` counts the factorizations and the retries of the condition
which needed the most, since each condition refreshes its own
factorization, and the Newton corrections of the chunk which iterated
longest: the conditions of a chunk are corrected together until the
last of them has converged, so that count falls as a host batch split
across threads makes its chunks smaller, while the results do not
change.
"""
struct TransientBatchSolution
    problems::Vector{TransientProblem}
    method::AbstractTransientIntegrator
    dt::Float64
    times::Vector{Float64}
    # every array untyped, as in a `TransientSolution`, so that a batch is
    # one type whatever its backend and its record level, a range of its
    # conditions, a view of its arrays, is that type, and the responses
    # compile once
    voltage::Any
    incident::Any
    outgoing::Any
    phases::Any
    endphases::Any
    endrates::Any
    linewaves::Any
    history::Any
    flux::Any
    rate::Any
    stages::Any
    checkpoints::Any
    finalwaves::Any
    finalstates::Any
    initialflux::Any
    initialrate::Any
    finalflux::Any
    finalrate::Any
    blockstates::Any
    initialwaves::Any
    initialstates::Any
    stats::NamedTuple
end
Base.length(b::TransientBatchSolution) = length(b.problems)

# The arrays a solution and a batch hold alike, in the order of their
# fields, the checkpoints as their named tuple; a solution is built from a
# named tuple of them, so that no array lands in another's place.
const SOLUTIONARRAYS = (:voltage, :incident, :outgoing, :phases, :endphases, :endrates, :linewaves, :history, :flux,
    :rate, :stages, :checkpoints, :finalwaves, :finalstates, :initialflux, :initialrate, :finalflux, :finalrate,
    :blockstates, :initialwaves, :initialstates)
TransientSolution(problem, method, dt, times, arrays::NamedTuple, stats) =
    TransientSolution(problem, method, dt, times, values(arrays[SOLUTIONARRAYS])..., stats)
TransientBatchSolution(problems, method, dt, times, arrays::NamedTuple, stats) =
    TransientBatchSolution(problems, method, dt, times, values(arrays[SOLUTIONARRAYS])..., stats)

# the arrays of a solution or a batch, each mapped by `f`, and the
# checkpoints' each (see `mapcheckpoints`), by name
solutionarrays(f, s) = NamedTuple{SOLUTIONARRAYS}(map(k -> k === :checkpoints ? mapcheckpoints(f, s.checkpoints) :
    f(getfield(s, k)), SOLUTIONARRAYS))

# a range of conditions is a batch of views, so that the noise can tile
# the conditions of a batch within its memory
function Base.getindex(b::TransientBatchSolution, js::AbstractVector{<:Integer})
    all(j -> 1 <= j <= length(b), js) || throw(BoundsError(b, js))
    slice = a -> isnothing(a) ? nothing : selectdim(a, ndims(a), js)
    return TransientBatchSolution(b.problems[js], b.method, b.dt, b.times, solutionarrays(slice, b), b.stats)
end
# a condition is a solution of views, its first and last states vectors of
# their own
function Base.getindex(b::TransientBatchSolution, j::Integer)
    1 <= j <= length(b) || throw(BoundsError(b, j))
    slice = a -> isnothing(a) ? nothing : selectdim(a, ndims(a), j)
    states = (; initialflux = copy(view(b.initialflux, :, j)), initialrate = copy(view(b.initialrate, :, j)),
        finalflux = copy(view(b.finalflux, :, j)), finalrate = copy(view(b.finalrate, :, j)))
    return TransientSolution(b.problems[j], b.method, b.dt, b.times, merge(solutionarrays(slice, b), states), b.stats)
end

# The factorizations of a batch of Gauss-Legendre stage matrices, one per
# condition on one pattern: on the CPU a complex matrix and a KLU
# factorization per condition, on a device one complex value matrix with
# a column per condition and the uniform cuDSS batch over it, in chunks of
# at most the batch size cuDSS solves correctly. `solve!` applies every
# factorization to its column of a right hand side matrix. A device
# chunk's sweep holds the solution and right hand side arrays of its
# conditions, a column each, allocated when it is built; the host has no
# chunks.
struct GaussBatchFactor{J, F, S, I}
    ncolumns::Int
    jacobians::J
    factors::F
    scratch::S
    chunks::Vector{UnitRange{Int}}
    rowptr::I
    colind::I
end

function gaussbatchfactor(sys::TransientSystem, ncolumns::Int)
    g = sys.gauss
    backend = sys.backend
    nnzj = nnz(sys.jacobian)
    if backend isa CPU
        pattern = g.cjacobian
        jacobians = [SparseMatrixCSC(size(pattern)..., SparseArrays.getcolptr(pattern), rowvals(pattern),
            zeros(ComplexF64, nnzj)) for _ in 1:ncolumns]
        # the assembly scratch belongs to the factor, not to the system,
        # so that several batches of one circuit may step at once
        return GaussBatchFactor(ncolumns, jacobians, Vector{Any}(nothing, ncolumns),
            zeros(Float64, nnzj), UnitRange{Int}[], nothing, nothing)
    end
    limit = uniformbatchlimit(1)
    chunks = [first:min(first + limit - 1, ncolumns) for first in 1:limit:ncolumns]
    nzval = KernelAbstractions.zeros(backend, ComplexF64, nnzj, ncolumns)
    A = g.cjacobian
    rowptr = tobackend(backend, convert(Vector{Int32}, rowpointer(A)))
    colind = tobackend(backend, convert(Vector{Int32}, columnindices(A)))
    return GaussBatchFactor(ncolumns, nzval, Vector{Any}(nothing, length(chunks)), nothing, chunks, rowptr, colind)
end

# the solution and right hand side arrays of a device sweep over the
# conditions `chunk`, which it keeps
sweeparrays(bf::GaussBatchFactor, sys::TransientSystem, chunk) =
    ntuple(_ -> KernelAbstractions.zeros(sys.backend, ComplexF64, size(sys.jacobian, 1), 1, length(chunk)), 2)

# The stage matrices at the stage phases of the columns of `columns`, a
# host mask, or of every column with `nothing`, assembled per column by
# the plan from the mean cosine of the two stages, and factorized, or
# refactorized on the analyses of the first time, and factorized anew by a
# method without an in place refactorization (QR). On the host only the
# columns asked for are factorized. A device batch refactorizes a whole
# chunk of the uniform batch at once, so the other columns of a chunk
# with one asked for keep the values of their last assembly and are
# refactorized to the factors they had; a chunk with none is left alone.
# A column's factor therefore depends on its own refreshes alone. With
# `fresh` the factorizations are new ones, their pivots chosen at these
# values rather than kept from the first; on the host a new one takes the
# fill reducing ordering the system chose for the pattern.
function gaussbatchjacobian!(bf::GaussBatchFactor, sys::TransientSystem, phi, cosphi, dwork, columns = nothing;
        fresh::Bool = false)
    isnothing(columns) || any(columns) || return bf
    g = sys.gauss
    r = sys.relations
    if allsinusoidal(r)
        cosphi .= (cos.(stage(phi, 1)) .+ cos.(stage(phi, 2))) ./ 2
    else
        derivativeinto!(dwork, r, stage(phi, 1))
        cosphi .= dwork
        derivativeinto!(dwork, r, stage(phi, 2))
        cosphi .= (cosphi .+ dwork) ./ 2
    end
    ncolumns = size(cosphi, 2)
    asked = j -> isnothing(columns) || columns[j]
    if sys.backend isa CPU
        for j in 1:ncolumns
            asked(j) || continue
            A = bf.jacobians[j]
            assemblerealjacobian!(bf.scratch, sys.plan, view(cosphi, :, j))
            nonzeros(A) .= bf.scratch .+ im .* g.imvals .+ g.rationalvals
            F = bf.factors[j]
            refreshed = (fresh || isnothing(F)) ? nothing : refactorize!(sys.factorization, F, A)
            bf.factors[j] = isnothing(refreshed) ? freshfactorization!(g.ordering, sys.factorization, A) : refreshed
        end
        return bf
    end
    nzval = bf.jacobians
    # every asked condition's assembly launched, then one synchronization
    for j in 1:ncolumns
        asked(j) || continue
        assemblerealjacobian!(nonzeros(sys.jacobian), sys.plan, view(cosphi, :, j); synchronize = false)
        view(nzval, :, j) .= nonzeros(sys.jacobian) .+ im .* g.imvals .+ g.rationalvals
    end
    KernelAbstractions.synchronize(sys.backend)
    rowptr, colind = bf.rowptr, bf.colind
    for (c, chunk) in enumerate(bf.chunks)
        any(asked, chunk) || continue
        if fresh || isnothing(bf.factors[c])
            bf.factors[c] = _cudss_sweep(rowptr, colind, nzval[:, chunk], sweeparrays(bf, sys, chunk)...;
                sys.factorization.kwargs...)
        else
            copyto!(bf.factors[c].nzval, view(nzval, :, chunk))
            _cudss_sweeprefactorize!(bf.factors[c])
        end
    end
    return bf
end

# the solve of the columns `first:first + nrhs - 1` of `R` into `Z` by a
# condition's factor, which the batch holds untyped: the views are taken
# past the dynamic call, where their types are known, so that the call
# boxes none of them
function solvecolumns!(Z, factor, R, first::Int, nrhs::Int)
    block = first:first + nrhs - 1
    matrixsolve!(view(Z, :, block), factor, view(R, :, block))
    return nothing
end

# the solve of every condition's block of right hand sides: the columns
# of `R` are `nrhs` per condition, condition after condition
function gaussbatchsolve!(Z::AbstractMatrix, bf::GaussBatchFactor, R::AbstractMatrix)
    n = size(R, 1)
    nrhs = size(R, 2) ÷ bf.ncolumns
    if isempty(bf.chunks)
        for j in 1:bf.ncolumns
            solvecolumns!(Z, bf.factors[j], R, (j - 1)*nrhs + 1, nrhs)
        end
        return Z
    end
    R3, Z3 = reshape(R, n, nrhs, :), reshape(Z, n, nrhs, :)
    for (c, chunk) in enumerate(bf.chunks)
        S = bf.factors[c]
        S.B .= view(R3, :, :, chunk)
        _cudss_sweepapply!(S, S.X, S.B)
        view(Z3, :, :, chunk) .= S.X
    end
    return Z
end

# The stage transform on a batch: the correction `[J*]^{-1} r` of a real
# stage pair `r` through the complex solve of every condition, `r̃ = Tinv r`,
# `z = S^{-1} r̃_1` and `d = T (z, conj z)`, which is real.
function gaussbatchtransform!(d, r, bf::GaussBatchFactor, gc::GaussCoefficients, rc, zc)
    rc .= gc.tinv11 .* stage(r, 1) .+ gc.tinv12 .* stage(r, 2)
    gaussbatchsolve!(zc, bf, rc)
    stage(d, 1) .= 2 .* real.(gc.t11 .* zc)
    stage(d, 2) .= 2 .* real.(gc.t21 .* zc)
    return d
end

# the junction Jacobian of every condition applied to its block of
# directions: `d` has `ncolumns` columns per condition, `phi` the phases
# of one stage as `(nj, N)`, and `work` is `(nj, size(d, 2))`
function batchjunctionproduct!(y, sys::TransientSystem, phi, d, work, dwork)
    nj, N = size(phi)
    stepmul!(work, sys.RJ, d)
    # the columns per condition, written out rather than left to a colon:
    # a circuit with no junctions has none of the rows the colon divides
    # by, which the long term support release will not reshape
    ncolumns = div(size(work, 2), N)
    r = sys.relations
    if allsinusoidal(r) && work isa Array
        # on the host a loop, without the reshaped arrays
        lmolj = sys.lmolj
        @inbounds for c in 1:N, col in 1:ncolumns, k in 1:nj
            work[k, (c - 1)*ncolumns + col] *= lmolj[k]*cos(phi[k, c])
        end
    elseif allsinusoidal(r)
        reshape(work, nj, ncolumns, N) .*= sys.lmolj .* cos.(reshape(phi, nj, 1, N))
    else
        derivativeinto!(dwork, r, phi)
        reshape(work, nj, ncolumns, N) .*= sys.lmolj .* reshape(dwork, nj, 1, N)
    end
    stepmul!(y, sys.RJt, work)
    return y
end

# the work of the rational blocks, declared ahead of the residual, which
# takes it or nothing unspecialized
abstract type AbstractRationalWork end

# The factorizations of a batch's two stages' matrices (see `StagePlan`),
# which the tangent and the adjoint solve at each step, one per condition:
# on the host a matrix on the plan's pattern and its factorization per
# condition, on a device one value matrix with a column per condition and
# the uniform cuDSS batch over it, in chunks of at most the batch size
# cuDSS solves correctly, of the matrix for the tangent and of its
# transpose for the adjoint, which a device does not solve transposed;
# `gather` assembles the orientation the factorization reads, and a device
# chunk's sweep holds the solution and right hand side arrays of its
# `nrhs` columns per condition. A response's first step pivots afresh at
# its own values, on the analysis the factorization holds, and its later
# steps refactorize, so that a response depends on its record alone and
# not on what the workspace solved before; a device's refactorization
# reuses its analysis alone. A constant matrix, without junctions or
# converted outputs, is factorized at the first step alone. `stiffness`
# holds the junctions' stiffness at both stages, `(2 nj, N)`, `derivative`
# one stage's derivative of their relations, `(nj, N)`, `weights` the
# blocks' output terms' weights at both stages, and `X` and `Y` on the
# host a condition's stacked right hand sides and their solution.
struct StageFactor{G, D, W, V}
    ncolumns::Int
    nrhs::Int
    adjoint::Bool
    gather::G
    matrices::Vector{SparseMatrixCSC{Float64, Int}}
    values::V
    factors::Vector{Any}
    chunks::Vector{UnitRange{Int}}
    stiffness::D
    derivative::D
    weights::W
    X::Matrix{Float64}
    Y::Matrix{Float64}
    fresh::Base.RefValue{Bool}
end

function stagefactor(sys::TransientSystem, ncolumns::Int, nrhs::Int; adjoint::Bool)
    plan = sys.gauss.stages
    backend = sys.backend
    n2, nnz2 = 2plan.n, nnz(plan.pattern)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    stiffness, derivative, weights = allocate(2plan.nj, ncolumns), allocate(plan.nj, ncolumns), allocate(2plan.nterms)
    if backend isa CPU
        colptr, rowval = SparseArrays.getcolptr(plan.pattern), rowvals(plan.pattern)
        matrices = [SparseMatrixCSC(n2, n2, colptr, rowval, zeros(nnz2)) for _ in 1:ncolumns]
        return StageFactor(ncolumns, nrhs, adjoint, plan.natural, matrices, zeros(0, 0), Vector{Any}(nothing, ncolumns),
            UnitRange{Int}[], stiffness, derivative, weights, zeros(n2, nrhs), zeros(n2, nrhs), Ref(true))
    end
    limit = uniformbatchlimit(nrhs)
    chunks = [first:min(first + limit - 1, ncolumns) for first in 1:limit:ncolumns]
    return StageFactor(ncolumns, nrhs, adjoint, adjoint ? plan.natural : plan.transposed, SparseMatrixCSC{Float64, Int}[],
        allocate(nnz2, ncolumns), Vector{Any}(nothing, length(chunks)), chunks, stiffness, derivative, weights, zeros(0, 0),
        zeros(0, 0), Ref(true))
end

# the solution and right hand side arrays of a device sweep over the
# conditions `chunk`, which it keeps
stagesweeparrays(sf::StageFactor, sys::TransientSystem, chunk) =
    ntuple(_ -> KernelAbstractions.zeros(sys.backend, Float64, 2sys.gauss.stages.n, sf.nrhs, length(chunk)), 2)

# the value of stored entry `e` of condition `j`'s matrix from a gather
# (see `StageGather`), the junctions' stiffness and the output terms'
# weights, summed in the gather's order, on the host and on a device alike
@inline function stagevalue(base, jptr, jrow, jcoef, stiffness, bptr, brow, bval, weights, e, j)
    @inbounds begin
        v = base[e]
        for c in Int(jptr[e]):Int(jptr[e + 1]) - 1
            v += jcoef[c]*stiffness[Int(jrow[c]), j]
        end
        for c in Int(bptr[e]):Int(bptr[e + 1]) - 1
            v += bval[c]*weights[Int(brow[c])]
        end
    end
    return v
end

# one work item per stored entry and condition
@kernel function stageassemblykernel!(values, @Const(base), @Const(jptr), @Const(jrow), @Const(jcoef), @Const(stiffness),
        @Const(bptr), @Const(brow), @Const(bval), @Const(weights))
    e, j = @index(Global, NTuple)
    @inbounds values[e, j] = stagevalue(base, jptr, jrow, jcoef, stiffness, bptr, brow, bval, weights, e, j)
end

# Every condition's matrix at its recorded stage phases `phi`, `(nj, N,
# 2)`, with the blocks' output terms at the step's weights `weights`, or
# nothing without blocks, factorized: pivoted afresh at a response's
# first step, refactorized at the steps after, and factorized anew where
# the method cannot do either in place.
function stagefactor!(sf::StageFactor, sys::TransientSystem, phi, weights)
    plan = sys.gauss.stages
    fresh = sf.fresh[]
    sf.fresh[] = false
    plan.constant && !fresh && return sf
    nj = plan.nj
    # each stage's stiffness, its relation's derivative taken into a whole
    # array, which a device's relation reads by index, and scaled into the
    # stage's rows
    for i in 1:2
        derivativeinto!(sf.derivative, sys.relations, view(phi, :, :, i))
        view(sf.stiffness, (i - 1)*nj + 1:i*nj, :) .= sys.lmolj .* sf.derivative
    end
    isnothing(weights) || copyto!(sf.weights, weights)
    g = sf.gather
    if sys.backend isa CPU
        for j in 1:sf.ncolumns
            A = sf.matrices[j]
            vals = nonzeros(A)
            for e in eachindex(vals)
                vals[e] = stagevalue(g.base, g.jptr, g.jrow, g.jcoef, sf.stiffness, g.bptr, g.brow, g.bval, sf.weights, e, j)
            end
            F = sf.factors[j]
            refreshed = isnothing(F) ? nothing : fresh ? pivotafresh!(F, sys.factorization, A) :
                refactorize!(sys.factorization, F, A)
            sf.factors[j] = isnothing(refreshed) ? freshfactorization!(plan.ordering, sys.factorization, A) : refreshed
        end
        return sf
    end
    stageassemblykernel!(sys.backend, 64)(sf.values, g.base, g.jptr, g.jrow, g.jcoef, sf.stiffness, g.bptr, g.brow, g.bval,
        sf.weights; ndrange = size(sf.values))
    KernelAbstractions.synchronize(sys.backend)
    for (c, chunk) in enumerate(sf.chunks)
        if isnothing(sf.factors[c])
            sf.factors[c] = _cudss_sweep(g.colptr, g.rowval, sf.values[:, chunk], stagesweeparrays(sf, sys, chunk)...;
                sys.factorization.kwargs...)
        else
            copyto!(sf.factors[c].nzval, view(sf.values, :, chunk))
            _cudss_sweeprefactorize!(sf.factors[c])
        end
    end
    return sf
end

# a factorization pivoted afresh at the values of `A`, on the pattern it
# was analyzed for: KLU's numeric factorization on the analysis it holds,
# or nothing for a method without one, which factorizes anew
function pivotafresh!(F::KLU.KLUFactorization, factorization::KLUfactorization, A::SparseMatrixCSC)
    F.nzval = nonzeros(A)
    return KLU.klu_factor!(F; factorization.kwargs...)
end
pivotafresh!(F, factorization, A) = nothing

# Every condition's directions solved on its factorization, or for the
# adjoint on the transposed matrix: `r` and `d` are `(n, nrhs N, 2)`, the
# directions of a condition contiguous, stacked by stage for the solve.
function stagesolve!(d, r, sf::StageFactor, sys::TransientSystem)
    n, nrhs = size(r, 1), sf.nrhs
    if sys.backend isa CPU
        X, Y = sf.X, sf.Y
        for j in 1:sf.ncolumns
            cols = (j - 1)*nrhs + 1:j*nrhs
            view(X, 1:n, :) .= view(r, :, cols, 1)
            view(X, n + 1:2n, :) .= view(r, :, cols, 2)
            stagesolvecolumns!(Y, sf.factors[j], X, sf.adjoint)
            view(d, :, cols, 1) .= view(Y, 1:n, :)
            view(d, :, cols, 2) .= view(Y, n + 1:2n, :)
        end
        return d
    end
    r4, d4 = reshape(r, n, nrhs, sf.ncolumns, 2), reshape(d, n, nrhs, sf.ncolumns, 2)
    for (c, chunk) in enumerate(sf.chunks)
        S = sf.factors[c]
        view(S.B, 1:n, :, :) .= view(r4, :, :, chunk, 1)
        view(S.B, n + 1:2n, :, :) .= view(r4, :, :, chunk, 2)
        _cudss_sweepapply!(S, S.X, S.B)
        view(d4, :, :, chunk, 1) .= view(S.X, 1:n, :, :)
        view(d4, :, :, chunk, 2) .= view(S.X, n + 1:2n, :, :)
    end
    return d
end

# a condition's stacked right hand sides `X` solved into `Y` on its
# factorization, which the batch holds untyped, or on its transpose,
# through the package's solves, which fall back to `\` for a method
# without an in place one (QR)
function stagesolvecolumns!(Y, factor, X, transposed::Bool)
    transposed ? trysolvetranspose!(Y, factor, X) : trysolve!(Y, factor, X)
    return nothing
end

# The rounding of the products a stage residual sums, row by row, per
# unit of roundoff, into `bound`, `(n, m, 2)` stage by stage: the
# capacitance and the conductance on both stages' increments `d` at the
# largest of the rule's weights, `a2` and `a1`, and the stiffness on the
# stage's own; `dsum` and `prod` are `(n, m)` work.
function stagerounding!(bound, sys::TransientSystem, a2, a1, d, dsum, prod)
    dsum .= abs.(stage(d, 1)) .+ abs.(stage(d, 2))
    stepmul!(prod, sys.Cabs, dsum)
    stage(bound, 1) .= a2 .* prod
    stepmul!(prod, sys.Gabs, dsum)
    stage(bound, 1) .+= a1 .* prod
    stage(bound, 2) .= stage(bound, 1)
    for i in 1:2
        dsum .= abs.(stage(d, i))
        stepmul!(prod, sys.Labs, dsum)
        stage(bound, i) .+= prod
    end
    return bound
end

# The rounding of the junction phases `RJ v` of each stage, carried to
# the rows by the junctions' stiffness at the stage phases `phi`,
# `(nj, N, 2)`: `|RJ'| (lmolj |relation'(phi)| (|RJ| |v|))`, added to
# `bound`, with `v` the stage values of the flux in a step, the columns of
# a condition contiguous; `vabs` and `prod` are `(n, m)` work, `jb`
# `(nj, m)` and `jd` the work of a polynomial relation's derivative.
function junctionrounding!(bound, sys::TransientSystem, phi, v, vabs, jb, jd, prod)
    nj, N = size(phi, 1), size(phi, 2)
    nj == 0 && return bound
    ncol = div(size(jb, 2), N)
    for i in 1:2
        vabs .= abs.(stage(v, i))
        stepmul!(jb, sys.RJabs, vabs)
        if allsinusoidal(sys.relations)
            reshape(jb, nj, ncol, N) .*= sys.lmolj .* reshape(abs.(cos.(view(phi, :, :, i))), nj, 1, N)
        else
            derivativeinto!(jd, sys.relations, view(phi, :, :, i))
            reshape(jb, nj, ncol, N) .*= sys.lmolj .* reshape(abs.(jd), nj, 1, N)
        end
        stepmul!(prod, sys.RJtabs, jb)
        stage(bound, i) .+= prod
    end
    return bound
end

# the largest residual of each column of `res` over its row's tolerance
# `tol` in a row where it exceeds that row's rounding bound, zero where
# every row is within its bound, on the host
residualexcess(res, bound, tol) = vec(Array(maximum(ifelse.(abs.(res) .> bound, abs.(res) ./ tol, 0.0); dims = (1, 3))))

# the steepest slope of the junctions' relations at the stage phases,
# `(nj, N, 2)`: one for the sinusoidal relation, and the largest
# derivative otherwise, `dwork` being `(nj, N)` work
function relationslope(sys::TransientSystem, phi, dwork)
    r = sys.relations
    allsinusoidal(r) && return 1.0
    slope = 0.0
    for i in 1:2
        derivativeinto!(dwork, r, view(phi, :, :, i))
        slope = max(slope, maximum(abs, dwork; init = 0.0))
    end
    return slope
end

recordedsolution(b::TransientBatchSolution) = recordedsolution(b[1])

# The currents of a tangent in one form, whatever the caller gave, on the
# host: the grid values `(target, time, direction)`, or one zero column
# `(target, 1, direction)` read at every time along the components alone,
# and the stage values `(target, stage, time, direction)` of a current
# given at the grid and the stages, empty otherwise. Converted once at
# the entry, so that a response compiles for one form of its currents
# rather than for every rank a caller may give. Where each direction is
# one waveform at one target, as a bath's quadrature is, `targets` holds
# the target of every direction and the grid and stage values hold the
# waveforms alone, one row, so that they grow with the directions rather
# than with the targets times the directions; `targets` is empty where
# the values are every target's.
struct TangentCurrents
    grid::Array{Float64, 3}
    stages::Array{Float64, 4}
    staged::Bool
    constant::Bool
    # the caller gave a matrix, one direction, and the outputs carry no
    # direction dimension
    single::Bool
    targets::Vector{Int}
end
directions(c::TangentCurrents) = size(c.grid, 3)
# the column of the grid values read at recorded time `k`
gridcolumn(c::TangentCurrents, k::Int) = c.constant ? 1 : k
function tangentcurrents(currents, nq::Int, nt::Int, ndir::Int)
    isnothing(currents) && return TangentCurrents(zeros(nq, 1, ndir), zeros(nq, 2, 0, ndir), false, true, false, Int[])
    staged = ndims(currents) == 4
    grid = reshape(Float64.(collect(staged ? selectdim(currents, 2, 1) : currents)), nq, nt, ndir)
    stages = staged ? reshape(Float64.(collect(selectdim(currents, 2, 2:3))), nq, 2, nt, ndir) : zeros(nq, 2, 0, ndir)
    return TangentCurrents(grid, stages, staged, false, ndims(currents) == 2, Int[])
end
# currents already in this form, the noise's, which `currentshape` checked
tangentcurrents(currents::TangentCurrents, nq::Int, nt::Int, ndir::Int) = currents
# the currents of a call released, their arrays of no times
releasedcurrents(c::TangentCurrents) = TangentCurrents(zeros(size(c.grid, 1), 0, size(c.grid, 3)),
    zeros(size(c.stages, 1), 2, 0, size(c.stages, 4)), c.staged, c.constant, c.single, c.targets)
# The currents of every target at grid column `k` on the host,
# `(targets, directions)`: a view of the grid values, or, where each
# direction drives one target, `out` holding its waveform there and
# zeros elsewhere.
function gridcurrents!(out, c::TangentCurrents, k::Int)
    isempty(c.targets) && return view(c.grid, :, k, :)
    fill!(out, 0)
    for (d, q) in enumerate(c.targets)
        out[q, d] = c.grid[1, k, d]
    end
    return out
end
# the selector of the targets of currents given one waveform per
# direction, one where a direction drives a target and zero elsewhere,
# `(targets, directions)`
function targetselector(targets::Vector{Int}, nq::Int)
    selector = zeros(nq, length(targets))
    for (d, q) in enumerate(targets)
        selector[q, d] = 1.0
    end
    return selector
end
# `out .=` the currents of every target from the values `v` a workspace
# holds at a time or a stage, and `out .+= weight .*` them, on the
# backend: the values themselves, or, through the selector of the
# targets, each direction's waveform at its target and zeros elsewhere
settargetcurrents!(out, ::Nothing, v) = (out .= v)
settargetcurrents!(out, selector, v) = (out .= ifelse.(selector .> 0, v, 0.0))
addtargetcurrents!(out, ::Nothing, v, weight) = (out .+= weight .* v)
addtargetcurrents!(out, selector, v, weight) = (out .+= weight .* ifelse.(selector .> 0, v, 0.0))

# The initial perturbation of a tangent in one form, on the host: the
# flux and the rate `(state, direction, condition)`, the line waves
# `(port, prehistory, direction, condition)` and the block states
# `(state, direction, condition)`, with one condition where the
# perturbation is shared by the conditions, and none where a member was
# not given.
struct TangentInitial
    flux::Array{Float64, 3}
    rate::Array{Float64, 3}
    waves::Array{Float64, 4}
    states::Array{Float64, 3}
    given::Bool
end
# a member of the perturbation for condition `j`, the one condition where
# it is shared, as a host array without the condition dimension
conditionslice(a::Array{Float64, 3}, j::Int) = collect(selectdim(a, 3, size(a, 3) == 1 ? 1 : j))
conditionslice(a::Array{Float64, 4}, j::Int) = collect(selectdim(a, 4, size(a, 4) == 1 ? 1 : j))
function tangentinitial(initialstate, n::Int, ndir::Int, N::Int, nl2::Int, npre::Int, nz::Int)
    isnothing(initialstate) && return TangentInitial(zeros(n, ndir, 0), zeros(n, ndir, 0), zeros(nl2, npre, ndir, 0),
        zeros(nz, ndir, 0), false)
    px, pv = initialstate
    (size(px, 1) == n && size(pv, 1) == n) || throw(DimensionMismatch(lazy"the state has $(n) entries."))
    percondition = ndims(px) == 3 && size(px, 3) > 1
    !percondition || (size(px, 3) == N && size(pv, 3) == N) || throw(DimensionMismatch(
        lazy"give one initial perturbation for every condition, or one per condition ($(N)) as the trailing dimension."))
    Nc = percondition ? N : 1
    (length(px) == n*ndir*Nc && length(pv) == n*ndir*Nc) || throw(DimensionMismatch(
        lazy"the initial perturbation has one column per direction ($(ndir))."))
    flux = reshape(Float64.(collect(px)), n, ndir, Nc)
    rate = reshape(Float64.(collect(pv)), n, ndir, Nc)
    waves = zeros(nl2, npre, ndir, 0)
    if length(initialstate) >= 3 && nl2 > 0
        pw0 = initialstate[3]
        (size(pw0, 1) == nl2 && size(pw0, 2) == npre) || throw(DimensionMismatch(
            lazy"the initial line waves need $(nl2) rows and $(npre) prehistory columns at this delay and step."))
        Nw = ndims(pw0) == 4 ? size(pw0, 4) : 1
        (Nw == 1 || Nw == N) && length(pw0) == nl2*npre*ndir*Nw || throw(DimensionMismatch(
            lazy"give one initial line history for every condition, or one per condition ($(N)) as the trailing dimension, with one column per direction ($(ndir))."))
        waves = reshape(Float64.(collect(pw0)), nl2, npre, ndir, Nw)
    end
    states = zeros(nz, ndir, 0)
    if length(initialstate) >= 4 && nz > 0
        pz = initialstate[4]
        size(pz, 1) == nz || throw(DimensionMismatch(lazy"the initial block states need $(nz) rows."))
        Nz = ndims(pz) == 3 ? size(pz, 3) : 1
        (Nz == 1 || Nz == N) && length(pz) == nz*ndir*Nz || throw(DimensionMismatch(
            lazy"give one initial block state for every condition, or one per condition ($(N)) as the trailing dimension, with one column per direction ($(ndir))."))
        states = reshape(Float64.(collect(pz)), nz, ndir, Nz)
    end
    return TangentInitial(flux, rate, waves, states, true)
end

# the scaled injection of a unit current at each target, and the port
# each target is, zero for a component, converted once at the entry of a
# response for whatever names the targets
function targetinjection(p::TransientProblem, targets)
    return (p.Lscale/phi0) .* transientinjection(p, targets), Vector{Int}(targetports(p, targets))
end

# the weights of an adjoint as `(port, time, objective)` on the host,
# checked
function adjointweights(weights, np::Int, nt::Int)
    ndims(weights) in (2, 3) && size(weights, 1) == np && size(weights, 2) == nt || throw(DimensionMismatch(
        lazy"weights must have one row per port ($(np)), one column per recorded time ($(nt)) and optionally a third dimension of objectives."))
    all(isfinite, weights) || throw(ArgumentError("the weights must be finite."))
    return reshape(Float64.(collect(weights)), np, nt, :)
end

# an array of a solution, untyped there, as a host matrix of a known type
# for the setup of a response
hostmatrix(a, n::Int, N::Int) = reshape(Float64.(collect(a)), n, N)::Matrix{Float64}

"""
    StageCorrection

The exact stage solve of a circuit with a pumped block. The frozen
operator of a step carries every block's unconverted response alone, as
the system assembles it once; a pumped block's converted outputs, with
the weights of each stage's time, act on the block's port rows alone, so
the true stage operator is the frozen one plus a correction of rank
twice the block's ports, `J = J* + U V'`, with `U` the scatter onto the
block's rows at each stage and, at stage `i`, `V' = -sum_j delta_ij W_ij`:
the outputs `W_ij` of each converted term `j` on the stacked states from
the stage unknowns, weighted by the term's weight `delta_ij` at the
stage's time. It is taken exactly by the Woodbury identity: a solve
`c = J*^(-1) r` is corrected to `c - K M^(-1) V' c` with `K = J*^(-1) U`
and `M = I + V' K` per column of the batch. `K` and the products
`W_ij K` depend on the factorization alone and are built when it is
refreshed; `M` is combined from them at each step's weights. The
frozen operator never changes with the modulation, so a pumped block
refactorizes only when the junctions ask. The tangent and the adjoint
solve the true stage operator directly (see `StagePlan`).
"""
mutable struct StageCorrection{K, R, Y}
    # the scatter onto the modulated rows, `U`, as the columns of an
    # `(n, np)` array
    Ucols::Y
    # the solves at the last refresh, `K = J*^(-1) U`, `(n, m, 2, r)`
    K::K
    # the products with the solves per column, `W_ij K` as `(np, r, m,
    # nw)`, on the host, and the small matrices `M` factorized per column
    Q::Array{Float64, 4}
    M::Vector{LU{Float64,Matrix{Float64},Vector{Int}}}
    # scratch: a right hand side and a solution of the stage transform,
    # a product on the block's rows, `(np, m)`, its host copy, the
    # products of both stages on the host, `(r, m)`, and the multipliers
    # of the columns, `(m, r)`, on both
    rhs::R
    sol::R
    y::Y
    ys::Matrix{Float64}
    yh::Matrix{Float64}
    z::Y
    zh::Matrix{Float64}
    Mh::Matrix{Float64}
end

mutable struct RationalWork{M, A} <: AbstractRationalWork
    z::M
    dstack::M
    xstack::M
    tstack::M
    zwork::M
    zwork2::M
    ystack::M
    ywork::M
    swork::M
    source::A
    # the entrywise bound of the blocks' contribution, `(n, N, 2)`, and
    # the magnitudes of the stacked states it is built from
    sbound::A
    yabs::M
    weights::Matrix{Float64}
    endweights::Vector{Float64}
    # the exact stage solve of a pumped block, built at a refresh of the
    # factorization and combined at each step, or nothing
    correction::Union{Nothing,StageCorrection}
end

function rationalwork(p::TransientProblem, backend, n::Int, N::Int)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    nz = blockstates(p)
    nterms = 1 + sum(b -> length(b.modulations), p.blocks; init = 0)
    return RationalWork(allocate(nz, N), allocate(2n, N), allocate(2n, N), allocate(2n, N),
        allocate(nz, N), allocate(nz, N), allocate(2nz, N), allocate(2nz, N), allocate(n, N), allocate(n, N, 2),
        allocate(n, N, 2), allocate(2nz, N),
        ones(nterms, 2), ones(nterms), nothing)
end

# The rational blocks' work of a stepper, or nothing, in a cell typed by
# the union. A closure's type carries the types of the values it
# captures, so a closure capturing the work itself is one type for a
# circuit with blocks and another for one without, and the stepper
# holding it with them; capturing the cell, the closures and the stepper
# are one type for both, and the presence of the work is a branch.
struct RationalWorkCell{M, A}
    work::Union{Nothing, RationalWork{M, A}}
    # the only constructor, and it takes the parameters, which nothing
    # leaves unbound
    RationalWorkCell{M, A}(work) where {M, A} = new{M, A}(work)
end

# The work of a perturbation of the components, which the tangent and
# the adjoint workspaces hold within their own types; built here for
# the workspaces, filled and read in sensitivity.jl.
# The quantities a step's residual reads, over `N` conditions: the
# acceleration, the rate, the flux and the relation values, their
# derivative's forcing when a tangent asks for one, and the state and the
# stage increments of a Gauss-Legendre step. Typed by the backend's
# matrix and stage array, so that a workspace holds it as a union with
# nothing within its own types.
struct PerturbationWork{M, A}
    F::Union{Nothing, M}
    work::Union{Nothing, M}
    a::M
    w::M
    X::M
    f::M
    fprev::M
    x::M
    v::M
    delta::A
end
function perturbationwork(cp, sys::TransientSystem, N::Int; forcing::Bool)
    allocate = (dims...) -> KernelAbstractions.zeros(sys.backend, Float64, dims...)
    nc, n = length(cp.names), cp.n
    nj = length(sys.lmolj)
    return PerturbationWork(forcing ? allocate(nc*n, N) : nothing, forcing ? allocate(nc*n, N) : nothing,
        allocate(n, N), allocate(n, N), allocate(n, N), allocate(nj, N), allocate(nj, N),
        allocate(n, N), allocate(n, N), allocate(n, N, 2))
end

# the work of the contraction of every kind of entries over `N`
# conditions and `nobj` objectives: the gathered state and multipliers,
# their product and the sum into the components
struct EntryContraction{M}
    GA::M
    GM::M
    prod::M
    out::M
end
struct ContractionWork{M}
    C::EntryContraction{M}
    G::EntryContraction{M}
    L::EntryContraction{M}
    J::EntryContraction{M}
end
function contractionwork(cp, entries, N::Int, nobj::Int, backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    nc = length(cp.names)
    one = e -> EntryContraction(allocate(e.count, N), allocate(e.count, nobj*N), allocate(e.count, nobj), allocate(nc, nobj))
    return ContractionWork(one(entries.C), one(entries.G), one(entries.L), one(entries.J))
end

# The work of the endpoint of a Gauss-Legendre step, whose projection and
# reading are host computations: the state on the backend and on the
# host, the junction phases and relation values, the host junction
# incidence and relation table, the reading's rows `Q`, and for a tangent
# the forcing `(n, ncomponents N)` with its stacked buffer, for an adjoint
# the contraction work.
struct EndpointWork{M}
    F::Union{Nothing, Matrix{Float64}}
    work::Union{Nothing, Matrix{Float64}}
    cw::Union{Nothing, ContractionWork{Matrix{Float64}}}
    x::M
    v::M
    xh::Matrix{Float64}
    vh::Matrix{Float64}
    phi::Matrix{Float64}
    f::Matrix{Float64}
    RJ::SparseMatrixCSC{Float64, Int}
    relations::JunctionRelations{Matrix{Float64}, Vector{Int}}
    Q::SparseMatrixCSC{Float64, Int}
end
function endpointwork(cp, sys::TransientSystem, p::TransientProblem, N::Int, nobj::Int; forcing::Bool)
    n, nj, nc = cp.n, length(sys.lmolj), length(cp.names)
    pr = sys.projection
    RJ = p.RJ
    relations = isnothing(p.relations) ? emptyrelations(zeros(0)) : hostrelations(p.relations)
    Q = endpointinjection(pr, sparse(1.0I, n, n))
    return EndpointWork(forcing ? zeros(n, nc*N) : nothing, forcing ? zeros(nc*n, N) : nothing,
        forcing ? nothing : contractionwork(cp, cp.hostentries, N, nobj, CPU()),
        KernelAbstractions.zeros(sys.backend, Float64, n, N), KernelAbstractions.zeros(sys.backend, Float64, n, N),
        zeros(n, N), zeros(n, N), zeros(nj, N), zeros(nj, N), RJ, relations, Q)
end

# The workspace of the Gauss-Legendre tangent of a batch: every array the
# tangent steps with and every piece it reads, sized by the system, the
# conditions, the directions and the targets, built once by
# `gausstangentwork` and driven step by step by `gaussbatchtangent`, as
# the solve's stepper is. A `TransientReuse` keeps one, and a later
# tangent of the same shape takes it over with its factorizations and
# their analyses. The perturbation is untyped, dispatched on at the
# stages alone; its work, the rational blocks' work, the projection's,
# the reading's and the endpoint injections are unions within the
# backend's types, as in the stepper.
mutable struct GaussTangentWork{S, J, M, A, A4, F, LT, CO, O}
    sys::S
    # the shape: the conditions, the directions, the targets with their
    # scaled injection and their ports, the recorded times, the form of
    # the currents, whether the outputs are stored, and the perturbation
    N::Int
    ndir::Int
    nq::Int
    nt::Int
    staged::Bool
    constant::Bool
    storing::Bool
    injh::SparseMatrixCSC{Float64, Int}
    tp::Vector{Int}
    perturbation::Any
    withstates::Bool
    injection::J
    # the currents on the backend, at the grid and at the stages, every
    # target's, or each direction's waveform alone with the selector of
    # the targets, nothing otherwise (see `TangentCurrents`); the call's
    # currents on the host, and the host work of every target's at a time
    # where each direction drives one
    dI::A
    dS::A4
    selector::Union{Nothing, M}
    currents::TangentCurrents
    hcurrents::Matrix{Float64}
    # the tangent's state, and the stage right hand sides and unknowns
    dx::M
    dv::M
    dxnew::M
    lx::M
    cdv::M
    work::M
    r::A
    d::A
    # the lines: the ring of the perturbation of the wave leaving each
    # port, the currents the arriving waves force and the rates across the
    # ports, the rate of the forced currents a reading takes, the
    # stencils and the tables, and the start of the prehistory
    dwaves::A
    ddlinevalues::M
    dlinerates::M
    dforcedrate::M
    linework::M
    stencil::Matrix{Float64}
    dstencil::M
    linetables::LT
    tpre::Float64
    # the rational blocks: the perturbation of their states and the stage
    # values their sources are read at
    rwt::Union{Nothing, RationalWork{M, A}}
    zerostages::Union{Nothing, A}
    fullstages::Union{Nothing, A}
    # the junctions at the stages
    jwork::M
    phi::A
    dwork::M
    # the currents at a stage, or at a time for the ports' direct term,
    # and their injection
    stagecurrent::M
    stageinjection::M
    injectionall::M
    # the port outputs: views of the call's result on the columns of this
    # workspace's conditions, handed over by the reset, of no times where
    # a sink receives them; the coefficients of the three quantities and
    # the work of a time
    voltage::O
    incident::O
    outgoing::O
    coefficients::CO
    portwork::M
    portmap::J
    directwork::M
    directall::M
    outwork::Vector{M}
    reading::Union{Nothing, OutputReading{M, SparseArrays.UMFPACK.UmfpackLU{Float64, Int}}}
    # the factorizations of every condition's two stages' matrix
    stages::F
    # the projection of the endpoint: its work, the injections along its
    # rows and the perturbation's endpoint work; the perturbation's stage
    # work, and the direct term of the port waves along the components
    # which are port terminations, read once per time
    pw::Union{Nothing, ProjectionWork{M, Matrix{Float64}, M}}
    Zci::Union{Nothing, SparseMatrixCSC{Float64, Int}}
    Qinji::Union{Nothing, SparseMatrixCSC{Float64, Int}}
    pwork::Union{Nothing, PerturbationWork{M, A}}
    pend::Union{Nothing, EndpointWork{M}}
    dterm::Any
    # the stepper a record of checkpoints is replayed on, with its
    # factorizations, or nothing until one asks
    replay::Any
end

function gausstangentwork(sys::TransientSystem, N::Int, ndir::Int, nq::Int, nt::Int, injh::SparseMatrixCSC{Float64, Int},
        tp::Vector{Int}, staged::Bool, constant::Bool, storing::Bool, reads::Bool, perturbation, targeted::Bool)
    p = sys.problem
    backend = sys.backend
    n, np = length(p), length(p.portimpedances)
    m = ndir*N
    h = sys.h
    # the components' forcing of the stages and, where the endpoint is
    # projected, of the endpoint, from the recorded or replayed states,
    # and the direct term of the port waves along a port's own termination
    pwork = isnothing(perturbation) ? nothing : perturbationwork(perturbation, sys, N; forcing = true)
    pend = (isnothing(perturbation) || isnothing(sys.projection)) ? nothing :
        endpointwork(perturbation, sys, p, N, ndir; forcing = true)
    dterm = isnothing(perturbation) ? nothing : directterm(perturbation, sys, N, ndir)
    withstates = !isnothing(perturbation) && perturbation.states
    injection = devicesparse(injh, backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    rows = targeted ? 1 : nq
    dI = allocate(rows, constant ? 1 : nt, ndir)
    dS = allocate(rows, 2, staged ? nt : 0, ndir)
    selector = targeted ? allocate(nq, ndir) : nothing
    currents = TangentCurrents(zeros(rows, 0, ndir), zeros(rows, 2, 0, ndir), staged, constant, false, Int[])
    dx, dv, dxnew, lx, cdv, work = [allocate(n, m) for _ in 1:6]
    r, d = allocate(n, m, 2), allocate(n, m, 2)
    nl2 = 2length(p.lines)
    npre = lineprehistory(p, h)
    nring = linering(npre)
    dwaves = allocate(nl2, nring, m)
    ddlinevalues, dlinerates, dforcedrate, linework = allocate(nl2, m), allocate(nl2, m), allocate(nl2, m), allocate(n, m)
    stencil, dstencil = zeros(8, nl2), allocate(8, nl2)
    rwt = isnothing(sys.gauss.coupling) ? nothing : rationalwork(p, backend, n, m)
    zerostages = isnothing(rwt) ? nothing : allocate(n, m, 2)
    fullstages = isnothing(rwt) ? nothing : allocate(n, m, 2)
    nj = length(sys.lmolj)
    jwork, phi = allocate(nj, m), allocate(nj, N, 2)
    # the derivative of the relation at one stage, for a circuit which has
    # a polynomial one; the Josephson relation is broadcast in place
    dwork = allsinusoidal(sys.relations) ? allocate(nj, 0) : allocate(nj, N)
    stagecurrent = allocate(nq, ndir)
    stageinjection, injectionall = allocate(n, ndir), allocate(n, m)
    # the outputs are the reset's, views of the result of a call
    voltage, incident, outgoing = [outputview(allocate(np, 0, m), 1:m) for _ in 1:3]
    coefficients = (outputcoefficients(sys, :voltage), outputcoefficients(sys, :incident), outputcoefficients(sys, :outgoing))
    portwork = allocate(np, m)
    portmap = devicesparse(targetportmap(tp, np), backend)
    directwork, directall = allocate(np, ndir), allocate(np, m)
    outwork = [allocate(np, m) for _ in 1:3]
    reading = reads ? outputreading(sys, backend, n, N, m) : nothing
    pr = sys.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, m)
    Zci = isnothing(pr) ? nothing : pr.Zch*injh
    # the direction's current along the rows of the endpoint reading
    Qinji = isnothing(pr) ? nothing : endpointinjection(pr, injh)
    return GaussTangentWork(sys, N, ndir, nq, nt, staged, constant, storing, injh, tp, perturbation, withstates, injection,
        dI, dS, selector, currents, zeros(targeted ? nq : 0, targeted ? ndir : 0), dx, dv, dxnew, lx, cdv, work, r, d,
        dwaves, ddlinevalues, dlinerates, dforcedrate, linework, stencil, dstencil, linetables(p, backend), 0.0,
        rwt, zerostages, fullstages, jwork, phi, dwork, stagecurrent, stageinjection, injectionall,
        voltage, incident, outgoing, coefficients, portwork, portmap, directwork, directall, outwork, reading,
        stagefactor(sys, N, ndir; adjoint = false), pw, Zci, Qinji, pwork, pend, dterm, nothing)
end

# whether a kept tangent workspace has the shape asked for, the reading,
# the perturbation and the form of the currents
function sameshape(w::GaussTangentWork, sys::TransientSystem, N::Int, ndir::Int, nq::Int, nt::Int, injh, tp, staged::Bool,
        constant::Bool, storing::Bool, reads::Bool, perturbation, targeted::Bool)
    return w.sys === sys && w.N == N && w.ndir == ndir && w.nq == nq && w.nt == nt && w.staged == staged &&
        w.constant == constant && w.storing == storing && !isnothing(w.reading) == reads && w.injh == injh && w.tp == tp &&
        sameperturbation(w.perturbation, perturbation) && isnothing(w.selector) == !targeted
end

# The workspaces of the chunks of a response, one per chunk: the reuse's
# kept ones when every chunk's has the shape asked for, or new ones for
# the chunks whose kept one has not, kept in a new vector. A reuse holds
# the workspaces of the tangent and of the adjoint each as a vector by
# chunk, one chunk being the whole batch on a device or on one thread.
# `shape(N)` tells whether a kept workspace fits a chunk of `N`
# conditions, and `build(N)` builds one.
function chunkworkspaces(@nospecialize(reuse), field::Symbol, chunks, @nospecialize(shape), @nospecialize(build))
    kept = isnothing(reuse) ? nothing : getfield(reuse, field)
    fits = kept isa Vector && length(kept) == length(chunks) ? [shape(kept[c], length(ch)) for (c, ch) in enumerate(chunks)] : fill(false, length(chunks))
    all(fits) && return kept
    ws = Vector{Any}(undef, length(chunks))
    for (c, ch) in enumerate(chunks)
        ws[c] = fits[c] ? kept[c] : build(length(ch))
    end
    isnothing(reuse) || setfield!(reuse, field, ws)
    return ws
end

# the workspace at the start of a tangent: the currents and the initial
# perturbation in, the state and the histories cleared, fresh output
# arrays where they are handed over, and the stage factorizations to
# pivot afresh at the first step
function tangentreset!(w::GaussTangentWork, sol::TransientBatchSolution, currents::TangentCurrents, initial::TangentInitial,
        outputs)
    sys = w.sys
    p = sys.problem
    backend = sys.backend
    n, N, ndir = length(p), w.N, w.ndir
    copyto!(w.dI, currents.grid)
    copyto!(w.dS, currents.stages)
    isnothing(w.selector) || copyto!(w.selector, targetselector(currents.targets, w.nq))
    w.currents = currents
    fill!(w.dx, 0)
    fill!(w.dv, 0)
    fill!(w.dwaves, 0)
    isnothing(w.rwt) || fill!(w.rwt.z, 0)
    nl2 = 2length(p.lines)
    if initial.given
        if size(initial.states, 3) > 0 && !isnothing(w.rwt)
            for j in 1:N
                copyto!(view(w.rwt.z, :, (j - 1)*ndir + 1:j*ndir), conditionslice(initial.states, j))
            end
        end
        if size(initial.waves, 4) > 0 && nl2 > 0
            for j in 1:N
                fillhistory!(view(w.dwaves, :, :, (j - 1)*ndir + 1:j*ndir), tobackend(backend, conditionslice(initial.waves, j)), 1)
            end
        end
        for j in 1:N
            copyto!(view(w.dx, :, (j - 1)*ndir + 1:j*ndir), conditionslice(initial.flux, j))
            copyto!(view(w.dv, :, (j - 1)*ndir + 1:j*ndir), conditionslice(initial.rate, j))
        end
    end
    w.tpre = linestart(first(sol.times), lineprehistory(p, sys.h), sys.h)
    w.voltage, w.incident, w.outgoing = outputs
    w.stages.fresh[] = true
    return w
end

# the columns `cols` of a result array `(rows, times, columns)`, the form
# every workspace holds its outputs in, one chunk or many
outputview(a, cols::UnitRange{Int}) = view(a, :, :, cols)

# the currents the arriving waves force at time `t`, read from the
# history accepted through column `accepted`
function tangentlineread!(w::GaussTangentWork, t, accepted)
    far, readscale, _ = w.linetables
    linestencil!(w.stencil, w.sys.problem.lines, t, w.tpre, w.sys.h, accepted)
    copyto!(w.dstencil, w.stencil)
    wavegather!(w.ddlinevalues, w.dwaves, w.dstencil, far, readscale, w.sys.backend)
    return nothing
end

# the rate of the currents the arriving waves force at time `t`, by the
# central difference the solve's reading takes of the lines' (see
# `readstepper!`), from the history accepted through column `accepted`,
# into the workspace's `dforcedrate`
function tangentlinerate!(w::GaussTangentWork, t, accepted)
    delta = ratedelta(w.sys)
    tangentlineread!(w, t + delta, accepted)
    copyto!(w.dforcedrate, w.ddlinevalues)
    tangentlineread!(w, t - delta, accepted)
    w.dforcedrate .= (w.dforcedrate .- w.ddlinevalues) ./ (2delta)
    return w.dforcedrate
end

# the projected junctions' phases and read rates at time `k` from a
# response's window into the reading work `o`, and with `withstates` the
# state there into the perturbation's work `pwork`, the window holding
# the rates as the solve read them, the lines' forcing included
function windowreading!(o, pwork, withstates::Bool, k::Int, window)
    if !isnothing(o.rn.pw)
        window.endphases!(o.erdev, k)
        copyto!(o.hphi, o.erdev)
        window.endrates!(o.erdev, k)
        copyto!(o.er, o.erdev)
    end
    withstates && window.state!(pwork.x, pwork.v, k)
    return nothing
end

# the port outputs of time `k`, stored, or handed to the sink as the three
# port by column matrices, and the state handed to the state sink, with
# the waves leaving the line ports at the time, the ring's newest; where
# a port or the state sink reads a rate along an algebraic direction,
# from the tangent rate read as the solve's is, the reading linearized at
# the projected junctions' recorded phases and read rates, with the
# tangent currents' rate and the components' perturbation of the
# constraints
function tangentoutputs!(w::GaussTangentWork, sol::TransientBatchSolution, k::Int, window, @nospecialize(outputsink),
        @nospecialize(statesink))
    sys = w.sys
    np, ndir, N = size(w.portwork, 1), w.ndir, w.N
    reading, pwork = w.reading, w.pwork
    if !isnothing(reading)
        windowreading!(reading, pwork, w.withstates, k, window)
        linerate = size(w.dwaves, 1) == 0 ? nothing : tangentlinerate!(w, sol.times[k], lineprehistory(sys.problem, sys.h) + k - 1)
        readoutputs!(reading, sys, k, w.dv, w.dx, w.currents, w.injh, w.withstates ? Array(pwork.x) : zeros(0, 0),
            w.withstates ? Array(pwork.v) : zeros(0, 0), w.perturbation, sys.backend, linerate)
    end
    stepmul!(w.portwork, sys.ports, isnothing(reading) ? w.dv : reading.dvread)
    w.portwork .*= phi0
    current = view(w.dI, :, w.constant ? 1 : k, :)
    if isnothing(w.selector)
        stepmul!(w.directwork, w.portmap, current)
    else
        settargetcurrents!(w.stagecurrent, w.selector, current)
        stepmul!(w.directwork, w.portmap, w.stagecurrent)
    end
    reshape(w.directall, np, ndir, N) .= reshape(w.directwork, np, ndir, 1)
    if w.storing
        tangentwrite!((view(w.voltage, :, k, :), view(w.incident, :, k, :), view(w.outgoing, :, k, :)), w, sol, k)
    else
        tangentwrite!((w.outwork[1], w.outwork[2], w.outwork[3]), w, sol, k)
        outputsink(k, w.outwork[1], w.outwork[2], w.outwork[3])
    end
    isnothing(statesink) || statesink(k, w.dx, isnothing(reading) ? w.dv : reading.dvread, isnothing(w.rwt) ? nothing : w.rwt.z,
        size(w.dwaves, 1) == 0 ? nothing :
            view(w.dwaves, :, ringslot(lineprehistory(sys.problem, sys.h) + k - 1, size(w.dwaves, 2)), :))
    return nothing
end
# the three outputs of time `k` written into `outs`, the stored ones or the
# sink's, each form on its own specialization
function tangentwrite!(outs, w::GaussTangentWork, sol::TransientBatchSolution, k::Int)
    for s in 1:3
        cv, cd = w.coefficients[s]
        outs[s] .= cv .* w.portwork .+ cd .* w.directall
    end
    isnothing(w.dterm) || adddirectterm!(outs, w.dterm, sol, sol.problems, k)
    return nothing
end

# the projection of the endpoint of step `k`, linearized about the
# recorded endpoint: `Z' (K dx - inj dI)` per column, solved with the small
# Jacobians of the condition; then the linearized reading of the index one
# unknowns, whose junction term is the stiffness at the endpoint times the
# flux direction's phase and whose drive is the direction's current at the
# endpoint time along the rows
function tangentproject!(w::GaussTangentWork, sol::TransientBatchSolution, k::Int, window)
    sys, pr, pw, pend, rwt = w.sys, w.sys.projection, w.pw, w.pend, w.rwt
    N, ndir = w.N, w.ndir
    m = ndir*N
    nl2 = size(w.ddlinevalues, 1)
    kcol = w.constant ? 1 : k + 1
    window.endphases!(pw.phip, k + 1)
    copyto!(pw.hphi, pw.phip)
    isnothing(pend) || (endpointquantities!(pend, w.perturbation, window, pw.hphi, pr.pj, k + 1); endpointforcing!(pend, w.perturbation))
    if !isempty(pr.directions)
        constrainttangent!(pw, pr, w.dxnew)
        dIz = repeat(w.Zci*gridcurrents!(w.hcurrents, w.currents, kcol), 1, N)
        nl2 > 0 && (dIz .+= pr.Zcline*Array(w.ddlinevalues))
        isnothing(rwt) || (dIz .+= pr.Zcblock*restingwaves(sys, rwt.z, sol.times[k + 1]))
        isnothing(pend) || (dIz .+= pr.Zch*pend.F)
        pw.g .-= dIz
        w.dxnew .+= constraintcorrection!(pw, pr)
        pr.cubic && projectrate!(w.dv, pr, pw, sys.gauss.coefficients, sys.h, stage(w.d, 1), stage(w.d, 2), w.dxnew, w.dx)
    end
    if !isempty(pr.readrows) || !isempty(pr.auxrows)
        junction = if isempty(pr.pj)
            zeros(0, m)
        else
            stepmul!(pw.phim, pr.RJp, w.dxnew)
            hphim = Array(pw.phim)
            reshape(pr.lmoljp .* reshape(derivativeat(pr.relationsp, pw.hphi), :, 1, N) .*
                reshape(hphim, :, ndir, N), :, m)
        end
        dIq = repeat(w.Qinji*gridcurrents!(w.hcurrents, w.currents, kcol), 1, N)
        nl2 > 0 && (dIq .+= pr.Qline*Array(w.ddlinevalues))
        isnothing(rwt) || (dIq .+= pr.Qblock*restingwaves(sys, rwt.z, sol.times[k + 1]))
        isnothing(pend) || (dIq .+= pend.Q*pend.F)
        endpointread!(w.dv, w.dxnew, pr, pw, junction, dIq)
    end
    return nothing
end

# the tangent's step from time `k` to `k + 1` on the recorded stage phases
# of the window: the linearized stages,
#   r_i = -(L + J_i') dx_n + (A^{-1} 1)_i C dv_n / h + db_i,
# with the current at the stage time read off the grid, the same for every
# condition, the components' forcing of the stage and the lines' and the
# blocks' terms, solved on the two stages' matrix at the recorded phases
# and the step's weights, factorized at the step
function tangentstep!(w::GaussTangentWork, sol::TransientBatchSolution, k::Int, window)
    sys, rwt, pwork = w.sys, w.rwt, w.pwork
    gc = sys.gauss.coefficients
    h = sys.h
    N, ndir, nt = w.N, w.ndir, w.nt
    n = size(w.dx, 1)
    m = ndir*N
    nl2 = size(w.ddlinevalues, 1)
    npre = lineprehistory(sys.problem, h)
    _, _, sqrtz = w.linetables
    window.phases!(w.phi, k)
    w.withstates && stagestates!(pwork, window, k)
    stepmul!(w.lx, sys.L, w.dx)
    stepmul!(w.cdv, sys.C, w.dv)
    for i in 1:2
        batchjunctionproduct!(w.work, sys, view(w.phi, :, :, i), w.dx, w.jwork, w.dwork)
        if w.staged
            settargetcurrents!(w.stagecurrent, w.selector, view(w.dS, :, i, k, :))
        else
            indices, weights = gaussstencil(k, nt, gc.c[i])
            fill!(w.stagecurrent, 0)
            for (idx, wgt) in zip(indices, weights)
                addtargetcurrents!(w.stagecurrent, w.selector, view(w.dI, :, w.constant ? 1 : idx, :), wgt)
            end
        end
        stepmul!(w.stageinjection, w.injection, w.stagecurrent)
        ri = stage(w.r, i)
        reshape(w.injectionall, n, ndir, N) .= reshape(w.stageinjection, n, ndir, 1)
        ri .= w.injectionall .+ (gc.ainvone[i]/h) .* w.cdv .- w.lx .- w.work
        if nl2 > 0
            tangentlineread!(w, sol.times[k] + gc.c[i]*h, npre + k - 1)
            stepmul!(w.linework, sys.lineinjection, w.ddlinevalues)
            ri .+= w.linework
        end
        if !isnothing(w.perturbation)
            stagequantities!(pwork, w.perturbation, sys, gc, view(w.phi, :, :, i), i)
            perturbationforcing!(pwork.F, w.perturbation, pwork, pwork.work)
            ri .+= reshape(pwork.F, n, m)
        end
    end
    if !isnothing(rwt)
        # the blocks' states and the state's part of the stage currents
        # carry part of the right hand side, at the step's weights
        stageweights!(rwt, sys, sol.times[k] + gc.c[1]*h, sol.times[k] + gc.c[2]*h)
        endweights!(rwt, sys, sol.times[k + 1])
        stage(w.fullstages, 1) .= w.dx
        stage(w.fullstages, 2) .= w.dx
        rationalsources!(rwt, sys, w.zerostages, w.fullstages)
        w.r .+= rwt.source
    end
    stagefactor!(w.stages, sys, w.phi, isnothing(rwt) ? nothing : rwt.weights)
    stagesolve!(w.d, w.r, w.stages, sys)
    if !isnothing(rwt)
        stage(w.fullstages, 1) .= w.dx .+ stage(w.d, 1)
        stage(w.fullstages, 2) .= w.dx .+ stage(w.d, 2)
        rationalstates!(rwt, sys, w.d, w.fullstages)
    end
    d1, d2 = stage(w.d, 1), stage(w.d, 2)
    w.dxnew .= w.dx .+ gc.ex[1] .* d1 .+ gc.ex[2] .* d2
    w.dv .+= (gc.ev[1]/h) .* d1 .+ (gc.ev[2]/h) .* d2
    nl2 > 0 && tangentlineread!(w, sol.times[k + 1], npre + k - 1)
    isnothing(sys.projection) || tangentproject!(w, sol, k, window)
    if nl2 > 0
        # the waves leaving the ports at the endpoint into the history
        stepmul!(w.dlinerates, sys.linegather, w.dv)
        view(w.dwaves, :, ringslot(npre + k, size(w.dwaves, 2)), :) .= phi0 .* w.dlinerates ./ sqrtz .- w.ddlinevalues .* sqrtz ./ 2
    end
    copyto!(w.dx, w.dxnew)
    return nothing
end

# The tangent of a Gauss-Legendre solve or batch, over `ndir` directions
# for every condition, the currents shared by the conditions and the
# initial perturbation shared or one per condition as a trailing
# dimension; `transienttangent` checked the arguments. The stored
# outputs are allocated once per call, and on the host a batch is split
# into chunks of conditions across the threads of the session, each
# chunk on its own workspace writing its own columns of them, as the
# solve steps its chunks into views of its arrays; a sink receives the
# columns of every condition at each time, so a call with one runs on
# one task. The workspaces are the reuse's, which holds them untyped,
# or new ones, and the tangent is invoked dynamically on each (see
# transientsolve), so that its steps are compiled for the workspace's
# type alone.
function gaussbatchtangent(sol::TransientBatchSolution, currents::TangentCurrents, injh::SparseMatrixCSC{Float64, Int},
        tp::Vector{Int}, initial::TangentInitial, sys::TransientSystem, @nospecialize(outputsink), perturbation,
        @nospecialize(reuse); @nospecialize(statesink = nothing),
        chunks = isnothing(outputsink) && isnothing(statesink) ? batchchunks(sys.backend, length(sol.problems)) : [1:length(sol.problems)])
    N, ndir, nq, nt = length(sol.problems), directions(currents), size(injh, 2), length(sol.times)
    storing = isnothing(outputsink)
    # the rate along the algebraic directions is read where a port reads
    # it, and for a state sink wherever the circuit has such directions
    reads = sys.portsread || (!isnothing(statesink) && !(isnothing(sys.invariant) && isnothing(sys.projection)))
    targeted = !isempty(currents.targets)
    shape = (w, Nc) -> w isa GaussTangentWork && sameshape(w, sys, Nc, ndir, nq, nt, injh, tp, currents.staged, currents.constant,
        storing, reads, perturbation, targeted)
    build = Nc -> gausstangentwork(sys, Nc, ndir, nq, nt, injh, tp, currents.staged, currents.constant, storing, reads, perturbation,
        targeted)
    ws = chunkworkspaces(reuse, :tangent, chunks, shape, build)
    # the result of the call, fresh, so that the results of earlier calls
    # stay what they were
    np = length(sys.problem.portimpedances)
    voltage, incident, outgoing = [KernelAbstractions.zeros(sys.backend, Float64, np, storing ? nt : 0, ndir*N) for _ in 1:3]
    columns = ch -> (first(ch) - 1)*ndir + 1:last(ch)*ndir
    outputs = ch -> (outputview(voltage, columns(ch)), outputview(incident, columns(ch)), outputview(outgoing, columns(ch)))
    parts = runchunks(chunks) do c, ch
        length(chunks) == 1 ? Base.invokelatest(gaussbatchtangent, ws[1], sol, currents, initial, outputsink, outputs(1:N), statesink) :
            Base.invokelatest(gaussbatchtangent, ws[c], sol[ch], currents, conditionslice(initial, ch), outputsink, outputs(ch), statesink)
    end
    ends = joinconditions(parts)
    shape = a -> begin
        b = reshape(a, np, size(a, 2), ndir, N)
        currents.single ? reshape(b, np, size(a, 2), N) : b
    end
    return (; voltage = storing ? shape(voltage) : nothing, incident = storing ? shape(incident) : nothing,
        outgoing = storing ? shape(outgoing) : nothing, ends.finalflux, ends.finalrate, ends.finalwaves, ends.finalstates)
end

# the members of an initial perturbation for the conditions `ch`, the
# shared one and the absent one as they are
function conditionslice(initial::TangentInitial, ch)
    slice = a -> size(a, ndims(a)) <= 1 ? a : collect(selectdim(a, ndims(a), ch))
    return TangentInitial(slice(initial.flux), slice(initial.rate), slice(initial.waves), slice(initial.states), initial.given)
end

# the endpoint results of the chunks of a response joined along the
# conditions, field by field, a field absent from every chunk absent
function joinconditions(parts)
    length(parts) == 1 && return only(parts)
    return NamedTuple(k => (isnothing(getfield(first(parts), k)) ? nothing :
        cat((getfield(part, k) for part in parts)...; dims = ndims(getfield(first(parts), k)))) for k in keys(first(parts)))
end

# the tangent on its workspace: the reset onto the outputs of its
# conditions, the outputs of the first time, the steps of every window
# with the outputs of their endpoints, and the endpoint results shaped
# by direction and condition
function gaussbatchtangent(w::GaussTangentWork, sol::TransientBatchSolution, currents::TangentCurrents,
        initial::TangentInitial, @nospecialize(outputsink), outputs, @nospecialize(statesink))
    sys = w.sys
    N, ndir = w.N, w.ndir
    try
        tangentreset!(w, sol, currents, initial, outputs)
        windows = responsewindows(sol, sys, true; withstates = w.withstates, stepper = replaystepper!(w, sol, sys))
        for (i, window) in enumerate(windows)
            isnothing(window.replay) || window.replay()
            i == 1 && tangentoutputs!(w, sol, 1, window, outputsink, statesink)
            for k in window.steps
                checkchunks(k)
                tangentstep!(w, sol, k, window)
                tangentoutputs!(w, sol, k + 1, window, outputsink, statesink)
            end
        end
        KernelAbstractions.synchronize(sys.backend)
        single = currents.single
        shape = a -> begin
            b = reshape(a, size(a)[1:end-1]..., ndir, N)
            single ? reshape(b, size(b)[1:end-2]..., N) : b
        end
        # the tangent's rate along the algebraic directions read as the
        # solve's is, the constraints' rows perturbed as its directions
        # perturb them
        finalrate = copy(w.dv)
        nl2 = size(w.dwaves, 1)
        if !(isnothing(sys.invariant) && isnothing(sys.projection))
            linerate = nl2 == 0 ? nothing : tangentlinerate!(w, last(sol.times), lineprehistory(sys.problem, sys.h) + w.nt - 1)
            readtangentrate!(finalrate, w.dv, w.dx, sys, sol, N, ndir, w.perturbation, w.currents, w.injh, sys.backend, linerate)
        end
        # the waves leaving the line ports over the delay window before
        # the end, the history a tangent of the next record starts from
        finalwaves = nl2 == 0 ? nothing : shape(readhistory!(KernelAbstractions.zeros(sys.backend, Float64, nl2,
            lineprehistory(sys.problem, sys.h), ndir*N), w.dwaves, w.nt))
        return (; finalflux = shape(copy(w.dx)), finalrate = shape(finalrate), finalwaves,
            finalstates = isnothing(w.rwt) ? nothing : shape(copy(w.rwt.z)))
    finally
        tangentrelease!(w)
    end
end

# A workspace kept by a reuse holds its own arrays and factorizations
# between calls, and of a call only, in its stepper, the grid and the
# problems of the solution it last replayed, which it reads and never
# writes: what a call hands over and a workspace only reads or writes
# while it runs, the host currents of a tangent, the views of the
# result it writes the outputs or the currents into, and the sink of an
# adjoint with whatever it captured, the noise's accumulator among them,
# is released when the call returns, on an error as well, in the reset
# as in the steps, so that its lifetime is the call's.
function tangentrelease!(w::GaussTangentWork)
    w.currents = releasedcurrents(w.currents)
    w.voltage = w.incident = w.outgoing = emptyoutput(w.voltage)
    return w
end

# an empty output of the kind a workspace holds a view of the result as
emptyoutput(a) = outputview(similar(parent(a), size(a, 1), 0, size(a, 3)), 1:size(a, 3))

# The columns of an adjoint's currents at the recorded times are final a
# few steps after a step first touches them, the trapezoidal step touching
# its own two endpoints, the Gauss-Legendre stencil up to four grid
# points ahead, and the reading of a rate along an algebraic direction
# the four grid points behind a time: a ring of eight columns holds the
# pending ones, and a column is handed to the sink, with the direct
# feedthrough of a port's current into its own wave added, once no
# remaining step touches it, in decreasing time. A sink that stores every
# column is the default; the noise sinks each column into its bath
# contraction and stores none.
const RINGSLOTS = 8
# the sink and the feedthrough are held untyped, so that a ring is one
# type whatever sink a call gives; each is one dynamic call per emitted
# column
struct CurrentRing{A, M}
    slots::A
    column::M
    sink::Any
    feedthrough!::Any
end
function CurrentRing(backend, nq, nobj, sink, feedthrough!)
    slots = KernelAbstractions.zeros(backend, Float64, nq, nobj, RINGSLOTS)
    column = KernelAbstractions.zeros(backend, Float64, nq, nobj)
    return CurrentRing(slots, column, sink, feedthrough!)
end
ringadd!(ring::CurrentRing, j, values, weight) = (view(ring.slots, :, :, mod1(j, RINGSLOTS)) .+= weight .* values; nothing)

# the map of the targets' currents onto the ports they are, `(port,
# target)` on the host, `tp` holding each target's port, zero for a
# component
targetportmap(tp::Vector{Int}, np::Int) = sparse([q for q in tp if q > 0], [k for (k, q) in enumerate(tp) if q > 0],
    ones(count(>(0), tp)), np, length(tp))

# The pieces of an adjoint's currents every rule shares: the targets'
# transposed injection and port map on the backend; the feedthrough of a
# target's current into its own port wave, `cd` of the weights `w` on the
# port, added to the column of time `k` of the currents of `N`
# conditions, with `directwork` and `objectivework` its work over the
# ports and the targets; and the sink storing every column into
# `currents`.
adjointmaps(injh::SparseMatrixCSC, tp::Vector{Int}, np::Int, backend) =
    (devicesparse(sparse(transpose(injh)), backend), devicesparse(sparse(transpose(targetportmap(tp, np))), backend))
function currentfeedthrough(w, cd, portmapt, directwork, objectivework, N::Int)
    nq, nobj = size(objectivework)
    return (column, k) -> begin
        directwork .= cd .* view(w, :, k, :)
        stepmul!(objectivework, portmapt, directwork)
        if N == 1
            column .+= objectivework
        else
            reshape(column, nq, nobj, N) .+= reshape(objectivework, nq, nobj, 1)
        end
        nothing
    end
end
storingsink(currents) = (k, values) -> (copyto!(view(currents, :, k, :), values); nothing)
function ringemit!(ring::CurrentRing, j)
    slot = view(ring.slots, :, :, mod1(j, RINGSLOTS))
    ring.column .= slot
    ring.feedthrough!(ring.column, j)
    ring.sink(j, ring.column)
    fill!(slot, 0)
    return nothing
end

# The workspace of the Gauss-Legendre adjoint of a batch, the counterpart
# of `GaussTangentWork`: every array the adjoint steps with and every
# piece it reads, sized by the system, the conditions, the objectives and
# the targets, built once by `gaussadjointwork` and driven step by step by
# `gaussbatchadjoint`; a `TransientReuse` keeps one for a later adjoint
# of the same shape. The perturbation is untyped, dispatched on at the
# stages alone; its work is typed as the tangent's is.
mutable struct GaussAdjointWork{S, J, M, A, F, LT, V, O}
    sys::S
    # the shape: the conditions, the objectives, the targets with their
    # scaled transposed injection and their ports, the recorded times, the
    # quantity, whether the currents are stored, and the perturbation
    N::Int
    nobj::Int
    nq::Int
    nt::Int
    quantity::Symbol
    storing::Bool
    injh::SparseMatrixCSC{Float64, Int}
    tp::Vector{Int}
    perturbation::Any
    withstates::Bool
    injectiont::J
    portmapt::J
    # the weights on the backend and the coefficients of the quantity
    w::A
    cv::V
    cd::V
    # the feedthrough of a target's current into its port wave, the
    # currents, a view of the call's result on the columns of this
    # workspace's conditions handed over by the reset, of no times where
    # a sink receives them, and the ring of the pending current columns,
    # rebuilt on the sink of a call
    directwork::M
    targetwork::M
    objectivework::M
    currents::O
    ring::CurrentRing{A, M}
    # the cotangents of the state, and the stages' cotangents and
    # multipliers
    xbar::M
    vbar::M
    work::M
    lmu::M
    cmu::M
    wst::A
    mu::A
    # the junctions at the stages
    jwork::M
    phi::A
    dwork::M
    # the objective's weights on the ports and their state
    portwork::M
    objectivestate::M
    transposing::Union{Nothing, OutputTranspose{M, SparseArrays.UMFPACK.UmfpackLU{Float64, Int}}}
    # the final state of every condition, on the backend
    xf::M
    wf::M
    # the factorizations of every condition's two stages' matrix, which
    # solve its transpose
    stages::F
    # the lines: the ring of the cotangent of the wave leaving each port,
    # the forced currents' and the rates' cotangents, the stencils and the
    # tables, and the start of the prehistory
    abar::A
    forcedbar::M
    ratebar::M
    linework::M
    stencil::Matrix{Float64}
    dstencil::M
    linetables::LT
    tpre::Float64
    # the rational blocks: the cotangent of their states, and the part of
    # the flux's cotangent the stage values carry
    rwa::Union{Nothing, RationalWork{M, A}}
    xextra::Union{Nothing, M}
    # the projection transposed, with the rate along the directions where
    # a block is on one
    pw::Union{Nothing, ProjectionWork{M, Matrix{Float64}, M}}
    vrecbar::Union{Nothing, M}
    # the perturbation: the stage work, the contraction work, the endpoint
    # work and the sensitivities on the backend and on the host
    pwork::Union{Nothing, PerturbationWork{M, A}}
    pcwork::Union{Nothing, ContractionWork{M}}
    pend::Union{Nothing, EndpointWork{M}}
    sensd::Union{Nothing, A}
    sensh::Union{Nothing, Array{Float64, 3}}
    # the stepper a record of checkpoints is replayed on, with its
    # factorizations, or nothing until one asks
    replay::Any
end

function gaussadjointwork(sys::TransientSystem, N::Int, nobj::Int, nq::Int, nt::Int, quantity::Symbol,
        injh::SparseMatrixCSC{Float64, Int}, tp::Vector{Int}, storing::Bool, perturbation)
    p = sys.problem
    backend = sys.backend
    n, np = length(p), length(p.portimpedances)
    m = nobj*N
    h = sys.h
    # the components' forcing of the stages and of the endpoint against
    # the multipliers: the stages' contraction on the backend, the
    # endpoint's on the host where the projection works
    nc = isnothing(perturbation) ? 0 : length(perturbation.names)
    pwork = isnothing(perturbation) ? nothing : perturbationwork(perturbation, sys, N; forcing = false)
    pcwork = isnothing(perturbation) ? nothing : contractionwork(perturbation, perturbation.entries, N, nobj, backend)
    pend = (isnothing(perturbation) || isnothing(sys.projection)) ? nothing :
        endpointwork(perturbation, sys, p, N, nobj; forcing = false)
    withstates = !isnothing(perturbation) && perturbation.states
    sensd = isnothing(perturbation) ? nothing : KernelAbstractions.zeros(backend, Float64, nc, nobj, N)
    cv, cd = outputcoefficients(sys, quantity)
    injectiont, portmapt = adjointmaps(injh, tp, np, backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    w = allocate(np, nt, nobj)
    directwork, targetwork, objectivework = allocate(np, nobj), allocate(nq, m), allocate(nq, nobj)
    # the currents are the reset's, a view of the result of a call
    currents = outputview(allocate(nq, 0, m), 1:m)
    ringslots, ringcolumn = allocate(nq, m, RINGSLOTS), allocate(nq, m)
    xbar, vbar, work, lmu, cmu = [allocate(n, m) for _ in 1:5]
    wst, mu = allocate(n, m, 2), allocate(n, m, 2)
    nj = length(sys.lmolj)
    jwork, phi = allocate(nj, m), allocate(nj, N, 2)
    dwork = allsinusoidal(sys.relations) ? allocate(nj, 0) : allocate(nj, N)
    portwork = allocate(np, nobj)
    objectivestate = allocate(n, nobj)
    transposing = sys.portsread ? outputtranspose(sys, backend, n, N, m) : nothing
    xf, wf = allocate(n, N), allocate(n, N)
    nl2 = 2length(p.lines)
    npre = lineprehistory(p, h)
    nring = linering(npre)
    abar = allocate(nl2, nring, m)
    forcedbar, ratebar, linework = allocate(nl2, m), allocate(nl2, m), allocate(n, m)
    stencil, dstencil = zeros(8, nl2), allocate(8, nl2)
    rwa = isnothing(sys.gauss.coupling) ? nothing : rationalwork(p, backend, n, m)
    xextra = isnothing(rwa) ? nothing : allocate(n, m)
    pr = sys.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, m)
    vrecbar = !isnothing(pr) && pr.cubic ? allocate(n, m) : nothing
    return GaussAdjointWork(sys, N, nobj, nq, nt, quantity, storing, injh, tp, perturbation, withstates, injectiont, portmapt,
        w, cv, cd, directwork, targetwork, objectivework, currents, CurrentRing(ringslots, ringcolumn, nothing, nothing),
        xbar, vbar, work, lmu, cmu, wst, mu, jwork, phi, dwork,
        portwork, objectivestate, transposing, xf, wf, stagefactor(sys, N, nobj; adjoint = true),
        abar, forcedbar, ratebar, linework, stencil, dstencil, linetables(p, backend), 0.0, rwa, xextra, pw, vrecbar,
        pwork, pcwork, pend, sensd, nothing, nothing)
end

# whether a kept adjoint workspace has the shape asked for and the
# perturbation
function sameshape(w::GaussAdjointWork, sys::TransientSystem, N::Int, nobj::Int, nq::Int, nt::Int, quantity::Symbol, injh, tp,
        storing::Bool, perturbation)
    return w.sys === sys && w.N == N && w.nobj == nobj && w.nq == nq && w.nt == nt && w.quantity == quantity &&
        w.storing == storing && w.injh == injh && w.tp == tp && sameperturbation(w.perturbation, perturbation)
end

# the workspace at the start of an adjoint: the weights and the final
# states of the conditions in, the cotangents and the histories cleared,
# the ring on the sink of the call with fresh current columns where they
# are handed over, the direct term of the objective along the port
# terminations, and the stage factorizations to pivot afresh at the first
# step
function adjointreset!(w::GaussAdjointWork, sol::TransientBatchSolution, wh::Array{Float64, 3}, @nospecialize(sink), currents)
    sys = w.sys
    p = sys.problem
    n, N = length(p), w.N
    copyto!(w.w, wh)
    for a in (w.xbar, w.vbar, w.abar, w.ring.slots)
        fill!(a, 0)
    end
    isnothing(w.rwa) || fill!(w.rwa.z, 0)
    isnothing(w.sensd) || fill!(w.sensd, 0)
    w.sensh = isnothing(w.perturbation) ? nothing :
        adddirectsensitivity!(zeros(length(w.perturbation.names), w.nobj, N), w.perturbation, sys, sol, sol.problems, wh, w.quantity)
    copyto!(w.xf, hostmatrix(sol.finalflux, n, N))
    copyto!(w.wf, hostmatrix(sol.finalrate, n, N))
    w.tpre = linestart(first(sol.times), lineprehistory(p, sys.h), sys.h)
    w.currents = currents
    w.ring = CurrentRing(w.ring.slots, w.ring.column, w.storing ? storingsink(currents) : sink,
        currentfeedthrough(w.w, w.cd, w.portmapt, w.directwork, w.objectivework, N))
    w.stages.fresh[] = true
    return w
end

# the ring's callbacks released when the call returns (see
# tangentrelease!): the sink of the call and what it captured go with it
function adjointrelease!(w::GaussAdjointWork)
    w.ring = CurrentRing(w.ring.slots, w.ring.column, nothing, nothing)
    w.currents = emptyoutput(w.currents)
    return w
end

# the cotangent of a forced current at time `t` scattered onto the far
# port's samples it read, `dq = 2 forced / sqrt(Z)`, with `sign`
function adjointlinescatter!(w::GaussAdjointWork, t, accepted, sign)
    far, readscale, _ = w.linetables
    linestencil!(w.stencil, w.sys.problem.lines, t, w.tpre, w.sys.h, accepted)
    copyto!(w.dstencil, w.stencil)
    wavescatter!(w.abar, w.forcedbar, w.dstencil, far, readscale, sign, w.sys.backend)
    return nothing
end

# the forced currents' cotangent from a multiplier on the equations
function adjointforced!(w::GaussAdjointWork, mu)
    stepmul!(w.forcedbar, w.sys.linegather, mu)
    w.forcedbar .*= w.sys.Lscale/phi0
    return nothing
end

# the projected junctions' phases and read rates at time `k`, and the
# states where the perturbation reads them: at the end from the final
# state, before it from the window, each with its rate as the solve read
# it, the lines' forcing included
function adjointreading!(w::GaussAdjointWork, k::Int, window)
    o, sys, pwork = w.transposing, w.sys, w.pwork
    prj, pwn = sys.projection, o.rn.pw
    if k == w.nt
        if !isnothing(pwn)
            stepmul!(pwn.phip, prj.RJp, w.xf)
            copyto!(o.hphi, pwn.phip)
            stepmul!(pwn.phip, prj.RJp, w.wf)
            copyto!(o.er, pwn.phip)
        end
        copyto!(o.wk, w.wf)
        w.withstates && (copyto!(pwork.x, w.xf); copyto!(pwork.v, w.wf))
    else
        windowreading!(o, pwork, w.withstates, k, window)
        w.withstates && copyto!(o.wk, pwork.v)
    end
    return nothing
end

# The adjoint of an output at time `k`: the rate receives
# `phi0 P' (cv .* w)`. Where a port reads a rate along an algebraic
# direction, through the transpose of the reading: the rate before the
# reading, the flux through the curvature of the projected junctions'
# stiffness, the currents whose rate the reading read, at the points of
# the stencil, the waves whose forced currents' rate it read, at the far
# ports' samples of the central difference (see `tangentlinerate!`), and
# the components whose perturbation of the constraints it read, on the
# read rate at the time.
function adjointoutput!(w::GaussAdjointWork, sol::TransientBatchSolution, k::Int, window)
    sys = w.sys
    n, nobj, N = size(w.xbar, 1), w.nobj, w.N
    w.portwork .= w.cv .* view(w.w, :, k, :)
    stepmul!(w.objectivestate, sys.portst, w.portwork)
    o = w.transposing
    if isnothing(o)
        reshape(w.vbar, n, nobj, N) .+= phi0 .* reshape(w.objectivestate, n, nobj, 1)
        return nothing
    end
    reshape(o.rbar, n, nobj, N) .= phi0 .* reshape(w.objectivestate, n, nobj, 1)
    adjointreading!(w, k, window)
    readoutputstranspose!(o, sys, k, w.nt, w.vbar, w.xbar, w.ring, w.injectiont, w.targetwork, w.perturbation, w.pcwork,
        w.sensd, w.withstates ? w.pwork.x : o.wk, sys.backend)
    if size(w.forcedbar, 1) > 0
        # the forcing's cotangent `fbar` carried to the forced currents'
        # rate, `-lineinjection' fbar`, and by the difference to the
        # currents at `t + delta` and `t - delta`
        delta = ratedelta(sys)
        accepted = lineprehistory(sys.problem, sys.h) + k - 1
        adjointforced!(w, o.fbar)
        w.forcedbar ./= 2delta
        adjointlinescatter!(w, sol.times[k] + delta, accepted, -1.0)
        adjointlinescatter!(w, sol.times[k] - delta, accepted, 1.0)
    end
    return nothing
end

# the projection transposed, before the stages of step `k`: where a block
# is on a direction the rate along the directions first, the part of the
# objective's rate the cubic carries to the stages and to the state, then
# the endpoint, whose multiplier `gamma` reads the current at time
# `k + 1` and moves the objective by `K Z gamma`
function adjointproject!(w::GaussAdjointWork, sol::TransientBatchSolution, k::Int, window)
    sys, pr, pw, pend, rwa = w.sys, w.sys.projection, w.pw, w.pend, w.rwa
    N, nobj = w.N, w.nobj
    nl2 = size(w.forcedbar, 1)
    npre = lineprehistory(sys.problem, sys.h)
    ew = sys.gauss.coefficients.endrate
    h = sys.h
    window.endphases!(pw.phip, k + 1)
    copyto!(pw.hphi, pw.phip)
    isnothing(pend) || endpointquantities!(pend, w.perturbation, window, pw.hphi, pr.pj, k + 1)
    # the reading of the index one unknowns transposed, last in the step
    # so first here: its row vector moves the current at the endpoint time
    cosp = tobackend(sys.backend, derivativeat(pr.relationsp, pw.hphi))
    if endpointreadtranspose!(w.vbar, w.xbar, pr, pw, sys, cosp, nobj, N)
        stepmul!(w.targetwork, w.injectiont, pw.work)
        ringadd!(w.ring, k + 1, w.targetwork, -1.0)
        isnothing(pend) || endpointcontract!(w.sensh, pend, w.perturbation, Array(pw.work), -1.0)
        if nl2 > 0
            adjointforced!(w, pw.work)
            adjointlinescatter!(w, sol.times[k + 1], npre + k - 1, -1.0)
        end
        # the resting waves the reading saw came from the states
        isnothing(rwa) || restingwavesbar!(rwa, sys, pw.work, -1.0)
    end
    isempty(pr.directions) && return nothing
    if pr.cubic
        # the rate along the directions transposed: `v += Z R (vrec - v)`
        stepmul!(pw.zwork, pr.Zt, w.vbar)
        stepmul!(w.vrecbar, pr.Zratet, pw.zwork)
        w.vbar .-= w.vrecbar
        w.xbar .+= (ew[4]/h) .* w.vrecbar
    end
    # the Newton step transposed (see `constrainttranspose!`): the
    # multiplier from the flux along the directions moves the flux before
    # the projection, and its image along the constraints' rows, which it
    # leaves in the work, reaches the current, the lines' forced currents
    # and the resting waves
    constrainttranspose!(w.xbar, w.lmu, pw, pr, sys, sys.backend; cosp)
    stepmul!(w.targetwork, w.injectiont, pw.work)
    ringadd!(w.ring, k + 1, w.targetwork, 1.0)
    isnothing(pend) || endpointcontract!(w.sensh, pend, w.perturbation, Array(pw.work), 1.0)
    if nl2 > 0
        adjointforced!(w, pw.work)
        adjointlinescatter!(w, sol.times[k + 1], npre + k - 1, 1.0)
    end
    isnothing(rwa) || restingwavesbar!(rwa, sys, pw.work, 1.0)
    return nothing
end

# the adjoint's step from time `k + 1` back to `k` on the recorded stage
# phases of the window: the cotangent of the wave that left each port at
# the endpoint, the projection transposed, the stage multipliers on the
# transposed two stages' matrix at the recorded phases and the step's
# weights, factorized at the step, and the multipliers' pull on the
# currents the stages read, the components, the lines and the state
function adjointstep!(w::GaussAdjointWork, sol::TransientBatchSolution, k::Int, window, @nospecialize(stagesink))
    sys, rwa, pwork = w.sys, w.rwa, w.pwork
    gc = sys.gauss.coefficients
    h = sys.h
    nt = w.nt
    nl2 = size(w.forcedbar, 1)
    npre = lineprehistory(sys.problem, h)
    _, _, sqrtz = w.linetables
    cubic = !isnothing(sys.projection) && sys.projection.cubic
    ew = gc.endrate
    window.phases!(w.phi, k)
    w.withstates && stagestates!(pwork, window, k)
    if !isnothing(rwa)
        stageweights!(rwa, sys, sol.times[k] + gc.c[1]*h, sol.times[k] + gc.c[2]*h)
        endweights!(rwa, sys, sol.times[k + 1])
    end
    if nl2 > 0
        # the wave that left each port at the endpoint: its cotangent goes
        # to the rate across the port and, with the opposite sign, to the
        # far port's samples the arriving wave read
        slot = ringslot(npre + k, size(w.abar, 2))
        w.ratebar .= phi0 .* view(w.abar, :, slot, :) ./ sqrtz
        w.forcedbar .= .- view(w.abar, :, slot, :) .* sqrtz ./ 2
        view(w.abar, :, slot, :) .= 0
        stepmul!(w.linework, sys.linescatter, w.ratebar)
        w.vbar .+= w.linework
        adjointlinescatter!(w, sol.times[k + 1], npre + k - 1, 1.0)
    end
    isnothing(sys.projection) || adjointproject!(w, sol, k, window)
    for i in 1:2
        stage(w.wst, i) .= gc.ex[i] .* w.xbar .+ (gc.ev[i]/h) .* w.vbar
        cubic && (stage(w.wst, i) .+= (ew[i + 1]/h) .* w.vrecbar)
    end
    if !isnothing(rwa)
        # the states after the step depend on the stage unknowns through
        # the update: their cotangent reaches the stages
        fill!(w.xextra, 0)
        statesbartostages!(rwa, sys, w.wst, w.xextra)
    end
    stagefactor!(w.stages, sys, w.phi, isnothing(rwa) ? nothing : rwa.weights)
    stagesolve!(w.mu, w.wst, w.stages, sys)
    # the states before the step: through the update, and through the
    # reflected waves the multipliers weigh, whose incident waves' value
    # part reaches the state's flux as well
    isnothing(rwa) || statesbarstep!(rwa, sys, w.mu, w.xextra)
    for i in 1:2
        mui = stage(w.mu, i)
        if !isnothing(w.perturbation)
            stagequantities!(pwork, w.perturbation, sys, gc, view(w.phi, :, :, i), i)
            contractperturbation!(w.sensd, w.perturbation, w.perturbation.entries, w.pcwork, pwork, mui, 1.0)
        end
        stepmul!(w.targetwork, w.injectiont, mui)
        if isnothing(stagesink)
            indices, wgts = gaussstencil(k, nt, gc.c[i])
            for (idx, wgt) in zip(indices, wgts)
                ringadd!(w.ring, idx, w.targetwork, wgt)
            end
        else
            stagesink(k, i, w.targetwork)
        end
        if nl2 > 0
            adjointforced!(w, mui)
            adjointlinescatter!(w, sol.times[k] + gc.c[i]*h, npre + k - 1, 1.0)
        end
        stepmul!(w.lmu, sys.Lt, mui)
        batchjunctionproduct!(w.work, sys, view(w.phi, :, :, i), mui, w.jwork, w.dwork)
        w.xbar .-= w.lmu .+ w.work
        stepmul!(w.cmu, sys.C, mui)
        w.vbar .+= (gc.ainvone[i]/h) .* w.cmu
    end
    cubic && (w.xbar .-= (ew[4]/h) .* w.vrecbar)
    isnothing(rwa) || (w.xbar .+= w.xextra)
    return nothing
end

# The adjoint of a Gauss-Legendre solve or batch, the weights shared by
# the conditions, the objectives of a condition contiguous in the columns
# the sinks receive; with a stage sink the stage multipliers go to it as
# they are, and the grid columns hold what the grid reads, the
# feedthrough and the projection; `transientadjoint` checked the
# arguments. The stored currents are allocated once per call, and on
# the host a batch is split into chunks across the threads of the
# session writing their own columns of them, as the tangent is; a call
# with a sink runs on one task. The workspaces are the reuse's or new
# ones, and the adjoint is invoked dynamically on each, as the tangent
# is.
function gaussbatchadjoint(sol::TransientBatchSolution, wh::Array{Float64, 3}, single::Bool, quantity::Symbol,
        injh::SparseMatrixCSC{Float64, Int}, tp::Vector{Int}, sys::TransientSystem, @nospecialize(sink),
        @nospecialize(stagesink), perturbation, @nospecialize(reuse);
        chunks = (isnothing(sink) && isnothing(stagesink)) ? batchchunks(sys.backend, length(sol.problems)) : [1:length(sol.problems)])
    N, nobj, nq, nt = length(sol.problems), size(wh, 3), size(injh, 2), length(sol.times)
    storing = isnothing(sink)
    shape = (w, Nc) -> w isa GaussAdjointWork && sameshape(w, sys, Nc, nobj, nq, nt, quantity, injh, tp, storing, perturbation)
    build = Nc -> gaussadjointwork(sys, Nc, nobj, nq, nt, quantity, injh, tp, storing, perturbation)
    ws = chunkworkspaces(reuse, :adjoint, chunks, shape, build)
    # the result of the call, fresh, so that the results of earlier calls
    # stay what they were
    currents = KernelAbstractions.zeros(sys.backend, Float64, nq, storing ? nt : 0, nobj*N)
    columns = ch -> (first(ch) - 1)*nobj + 1:last(ch)*nobj
    parts = runchunks(chunks) do c, ch
        Base.invokelatest(gaussbatchadjoint, ws[c], length(chunks) == 1 ? sol : sol[ch], wh, single, sink, stagesink,
            outputview(currents, columns(ch)))
    end
    ends = joinconditions(parts)
    shape = a -> (b = reshape(a, nq, size(a, 2), nobj, N); single ? reshape(b, nq, size(a, 2), N) : b)
    return (; currents = storing ? shape(currents) : nothing, ends.initialflux, ends.initialrate, ends.initialwaves,
        ends.initialstates, ends.sensitivity)
end

# the adjoint on its workspace: the reset onto the currents of its
# conditions, the objective at the final time, the steps of every window
# backwards with the objective at their start and the emission of the
# current columns no step touches, the columns of the first times, and
# the initial results shaped by objective and condition
function gaussbatchadjoint(w::GaussAdjointWork, sol::TransientBatchSolution, wh::Array{Float64, 3}, single::Bool,
        @nospecialize(sink), @nospecialize(stagesink), currents)
    sys = w.sys
    N, nobj, nt = w.N, w.nobj, w.nt
    try
        adjointreset!(w, sol, wh, sink, currents)
        windows = responsewindows(sol, sys, false; withstates = w.withstates, stepper = replaystepper!(w, sol, sys))
        adjointoutput!(w, sol, nt, first(windows))
        for window in windows
            isnothing(window.replay) || window.replay()
            for k in reverse(window.steps)
                checkchunks(k)
                adjointstep!(w, sol, k, window, stagesink)
                adjointoutput!(w, sol, k, window)
                k + 3 <= nt && ringemit!(w.ring, k + 3)
            end
        end
        for j in min(3, nt):-1:1
            ringemit!(w.ring, j)
        end
        KernelAbstractions.synchronize(sys.backend)
        m = nobj*N
        shape = a -> begin
            isnothing(a) && return nothing
            b = reshape(a, size(a)[1:end-1]..., nobj, N)
            single ? reshape(b, size(b)[1:end-2]..., N) : b
        end
        nl2 = size(w.forcedbar, 1)
        npre = lineprehistory(sys.problem, sys.h)
        initialwaves = nl2 > 0 ? shape(readhistory!(KernelAbstractions.zeros(sys.backend, Float64, nl2, npre, m), w.abar, 1)) : nothing
        initialstates = isnothing(w.rwa) ? nothing : shape(copy(w.rwa.z))
        sensitivity = isnothing(w.perturbation) ? nothing :
            shape(reshape(w.sensd .+ tobackend(sys.backend, w.sensh), length(w.perturbation.names), m))
        return (; initialflux = shape(copy(w.xbar)), initialrate = shape(copy(w.vbar)), initialwaves, initialstates, sensitivity)
    finally
        adjointrelease!(w)
    end
end

"""
    transienttangent(batch::TransientBatchSolution, currents; targets, initialstate)

The tangent of every condition of a batch along the same currents, all
conditions on one pass: the arrays of the single tangent with the
conditions as the trailing dimension. The initial perturbation is one
tuple for every condition, or one per condition with the conditions as a
trailing dimension of its arrays. On the host the conditions are
split across the threads of the session when the outputs are stored,
as the solve splits them; an `outputsink` or a `statesink` receives the
columns of every condition at each time, so a call with one runs on one
thread.
"""
Base.@nospecializeinfer function transienttangent(b::TransientBatchSolution,
        @nospecialize(currents::Union{Nothing,AbstractArray{<:Real},TangentCurrents});
        targets = porttargets(first(b.problems)), initialstate = nothing, factorization = nothing, reuse = nothing,
        outputsink = nothing, statesink = nothing, perturbation = nothing)
    # compiled once whatever the arguments are: they are brought to the
    # forms the steps take here, the noise's baths' quadratures given in
    # that form already, and the steps are invoked dynamically
    @nospecialize targets initialstate factorization reuse outputsink statesink perturbation
    p = first(b.problems)
    sys = responsesystem(b, p, factorization, reuse)
    nq, nt, N = length(targets), length(b.times), length(b)
    ndir = tangentdirections(currents, perturbation, nq, nt)
    injection, ports = targetinjection(p, targets)
    hostcurrents = tangentcurrents(currents, nq, nt, ndir)
    checktangentbalance(p, injection, hostcurrents, b.method)
    initial = tangentinitial(initialstate, length(p), ndir, N, 2length(p.lines), lineprehistory(p, b.dt), blockstates(p))
    # invoked dynamically on the untyped kept system (see transientsolve);
    # a record of another rule is a solution's, a batch of one (see
    # `batchof`), stepped by that rule's own tangent
    b.method isa GaussLegendre || return map(addcondition, Base.invokelatest(steptangent, unbatch(b), hostcurrents,
        injection, ports, initial, sys, sys.factorization, keptfactor(reuse), outputsink, statesink, perturbation, reuse))
    return Base.invokelatest(gaussbatchtangent, b, hostcurrents, injection, ports, initial, sys,
        outputsink, perturbation, reuse; statesink)
end

"""
    transientadjoint(batch::TransientBatchSolution, weights; quantity, targets, sink)

The adjoint of every condition of a batch for the same weights, all
conditions on one pass: the arrays of the single adjoint with the
conditions as the trailing dimension, and a sink's columns holding the
objectives of a condition contiguously, condition after condition. On
the host the conditions are split across the threads of the session
when the currents are stored, as the solve splits them; a call with a
`sink` or a `stagesink` runs on one thread.
"""
Base.@nospecializeinfer function transientadjoint(b::TransientBatchSolution, @nospecialize(weights::AbstractArray{<:Real});
        quantity::Symbol = :outgoing, targets = porttargets(first(b.problems)), factorization = nothing,
        reuse = nothing, sink = nothing, stagesink = nothing, components = String[])
    # compiled once whatever the arguments are (see transienttangent)
    @nospecialize targets factorization reuse sink stagesink components
    p = first(b.problems)
    sys = responsesystem(b, p, factorization, reuse)
    wh = adjointweights(weights, length(p.portimpedances), length(b.times))
    perturbation = isempty(components) ? nothing : componentperturbation(p, components, sys.backend; forcing = false)
    isnothing(perturbation) || recordedstates(b, perturbation)
    injection, ports = targetinjection(p, targets)
    checkadjointtargets(p, injection, targets)
    # invoked dynamically on the untyped kept system (see transientsolve);
    # a record of another rule is a solution's, stepped by that rule's own
    # adjoint, which has no stages to hand a stage sink (see transienttangent)
    b.method isa GaussLegendre || return map(addcondition, Base.invokelatest(stepadjoint, unbatch(b), wh,
        ndims(weights) == 2, quantity, injection, ports, sys, sys.factorization, keptfactor(reuse), sink, perturbation, reuse))
    return Base.invokelatest(gaussbatchadjoint, b, wh, ndims(weights) == 2, quantity, injection, ports, sys, sink,
        stagesink, perturbation, reuse)
end

# the stage residual of a batch on `(n, N, 2)` stage
# increments with the states as the columns of `x`, and its norm per
# column weighted by the tolerance `rowtol` of each row, from the
# columns' maxima in the host's `maxima`, reduced on a device in
# `devicemaxima` (see `columnmaxima!`)
function gaussbatchresidual!(norms, residual, sys::TransientSystem, gc::GaussCoefficients, delta, x, lx, X,
        phi, junction, jwork, cwork, gwork, rhs, roundoff, maxima, devicemaxima, blockcol, rowtol, dwork,
        @nospecialize(rw::Union{Nothing, AbstractRationalWork}))
    h = sys.h
    X .= x .+ delta
    # the rational blocks' reflected waves at the stages, from their states
    # and the incident waves, as sources on their rows
    isnothing(rw) || rationalsources!(rw, sys, delta, X)
    stepmul!(stage(cwork, 1), sys.C, stage(delta, 1)); stepmul!(stage(cwork, 2), sys.C, stage(delta, 2))
    stepmul!(stage(gwork, 1), sys.G, stage(delta, 1)); stepmul!(stage(gwork, 2), sys.G, stage(delta, 2))
    stepmul!(stage(phi, 1), sys.RJ, stage(X, 1)); stepmul!(stage(phi, 2), sys.RJ, stage(X, 2))
    if allsinusoidal(sys.relations)
        jwork .= sys.lmolj .* sin.(phi)
    else
        relationinto!(jwork, sys.relations, phi)
        jwork .= sys.lmolj .* jwork
    end
    stepmul!(stage(junction, 1), sys.RJt, stage(jwork, 1)); stepmul!(stage(junction, 2), sys.RJt, stage(jwork, 2))
    stepmul!(stage(residual, 1), sys.L, stage(delta, 1)); stepmul!(stage(residual, 2), sys.L, stage(delta, 2))
    c1, c2 = stage(cwork, 1), stage(cwork, 2)
    g1, g2 = stage(gwork, 1), stage(gwork, 2)
    for i in 1:2
        stage(residual, i) .+= lx .+ stage(junction, i) .- stage(rhs, i) .+
            (entry(gc.ainv2, i, 1)/h^2) .* c1 .+ (entry(gc.ainv2, i, 2)/h^2) .* c2 .+
            (entry(gc.ainv, i, 1)/h) .* g1 .+ (entry(gc.ainv, i, 2)/h) .* g2
        isnothing(rw) || (stage(residual, i) .-= stage(rw.source, i))
    end
    # the weighted norm per column over both stages, and the floor, on the
    # products' work, which the residual no longer needs
    columnmaxima!(maxima, devicemaxima, residual, rowtol, delta, x, X, jwork, rhs, sys.backend)
    for j in eachindex(norms)
        norms[j] = maxima[1, j]
    end
    gaussresidualroundoff!(roundoff, norms, maxima, sys, gc, residual, delta, x, X, phi, jwork, rhs, rowtol, cwork, gwork,
        blockcol, dwork, rw)
    return nothing
end

# The largest magnitudes of each column of a stage residual `res`, `(n,
# N, 2)`, and of what its roundoff floor reads, the rows of `maxima`,
# `(COLUMNMAXIMA, N)` on the host: the residual weighted by the
# tolerance `tol` of its row, `max_i |r_i|/tol_i`, which is at most one
# where every row is within its own tolerance; the largest residual of a
# row over its tolerance, zero where none is; and the largest increment
# `delta`, state `x`, junction current `jcur`, stage value `X` and right
# hand side `rhs`. A loop on the host; on a device one kernel into
# `devicemaxima` and one copy.
const COLUMNMAXIMA = 7
# the work items of a column's reduction on a device
const COLUMNGROUP = 256

function columnmaxima!(maxima::Matrix{Float64}, devicemaxima, res, tol, delta, x, X, jcur, rhs, backend)
    columnmaximakernel!(backend, COLUMNGROUP)(devicemaxima, res, tol, delta, x, X, jcur, rhs;
        ndrange = COLUMNGROUP*size(res, 2))
    copyto!(maxima, devicemaxima)
    return maxima
end

function columnmaxima!(maxima::Matrix{Float64}, devicemaxima, res, tol, delta, x, X, jcur, rhs, ::CPU)
    @inbounds for col in axes(res, 2)
        weighted, over, dm, xm, jm, Xm, rm = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
        for row in axes(res, 1)
            xm = max(xm, abs(x[row, col]))
            for i in 1:2
                r, t = abs(res[row, col, i]), tol[row, col, i]
                weighted = max(weighted, r/t)
                over = max(over, ifelse(r > t, r, 0.0))
                dm = max(dm, abs(delta[row, col, i]))
                Xm = max(Xm, abs(X[row, col, i]))
                rm = max(rm, abs(rhs[row, col, i]))
            end
        end
        for i in 1:2, k in axes(jcur, 1)
            jm = max(jm, abs(jcur[k, col, i]))
        end
        maxima[1, col], maxima[2, col], maxima[3, col], maxima[4, col] = weighted, over, dm, xm
        maxima[5, col], maxima[6, col], maxima[7, col] = jm, Xm, rm
    end
    return maxima
end

# one workgroup of `COLUMNGROUP` work items per column: each takes the
# maxima of its share of the rows, and the group halves them into the
# column's
@kernel function columnmaximakernel!(maxima, @Const(res), @Const(tol), @Const(delta), @Const(x), @Const(X), @Const(jcur),
        @Const(rhs))
    col = @index(Group, Linear)
    item = @index(Local, Linear)
    partial = @localmem Float64 (COLUMNGROUP, COLUMNMAXIMA)
    weighted, over, dm, xm, jm, Xm, rm = 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
    @inbounds begin
        row = item
        while row <= size(res, 1)
            xm = max(xm, abs(x[row, col]))
            for i in 1:2
                r, t = abs(res[row, col, i]), tol[row, col, i]
                weighted = max(weighted, r/t)
                over = max(over, ifelse(r > t, r, 0.0))
                dm = max(dm, abs(delta[row, col, i]))
                Xm = max(Xm, abs(X[row, col, i]))
                rm = max(rm, abs(rhs[row, col, i]))
            end
            row += COLUMNGROUP
        end
        k = item
        while k <= size(jcur, 1)
            for i in 1:2
                jm = max(jm, abs(jcur[k, col, i]))
            end
            k += COLUMNGROUP
        end
        partial[item, 1], partial[item, 2], partial[item, 3], partial[item, 4] = weighted, over, dm, xm
        partial[item, 5], partial[item, 6], partial[item, 7] = jm, Xm, rm
    end
    @synchronize
    stride = COLUMNGROUP ÷ 2
    while stride > 0
        if item <= stride
            @inbounds for q in 1:COLUMNMAXIMA
                partial[item, q] = max(partial[item, q], partial[item + stride, q])
            end
        end
        @synchronize
        stride ÷= 2
    end
    if item == 1
        @inbounds for q in 1:COLUMNMAXIMA
            maxima[q, col] = partial[1, q]
        end
    end
end

# The roundoff floor of each column's residual just formed, in the units
# of its weighted norm, which Newton stops on where a row cannot reach
# its tolerance: the column's norm where every row over its tolerance is
# at its own rounding, and zero otherwise. It is decided in two stages.
# The cheap one bounds the rounding of every term by a row sum of its
# matrix times the largest entry it multiplies, an upper bound of any
# row's rounding from a few reductions of the column: the capacitance,
# conductance and stiffness on the increments, the stiffness on the
# state, the junction currents and the rounding of their phases, whose
# difference of node fluxes may be far larger than the phase, the right
# hand side, and the blocks' reflected waves. A column with a row over
# its tolerance by more than that has not converged, and one without a
# row over its tolerance has. For a column in between, the rounding is
# bounded row by row from each row's own terms (see
# `stepresidualexcess!`), and the column is at its floor when every
# row's residual lies within its tolerance or that row's own rounding. A
# quiet node with a large capacitance, or a block's state far from
# another's, then sets the threshold of its own rows and not of the
# others. Base and trial evaluations keep separate floors, as they keep
# separate norms.
function gaussresidualroundoff!(roundoff, norms, maxima, sys::TransientSystem, gc::GaussCoefficients, residual, delta, x, X,
        phi, jwork, rhs, rowtol, bound, work, blockcol, dwork, @nospecialize(rw::Union{Nothing, AbstractRationalWork}))
    cs, gs, ls, js, jr = sys.rowsums
    h = sys.h
    a2, a1 = maximum(abs, gc.ainv2)/h^2, maximum(abs, gc.ainv)/h
    isnothing(rw) ? fill!(blockcol, 0.0) : stackedscatterbound!(blockcol, rw, sys.gauss.coupling)
    # the cheap bound from the column maxima (see `columnmaxima!`), per
    # unit of the increments, the state, the junction currents, the stage
    # values' phases and the right hand side, with the blocks' bound
    rounding, jphase = 2a2*cs + 2a1*gs + ls, js*relationslope(sys, phi, dwork)
    for j in eachindex(roundoff)
        roundoff[j] = 8eps(Float64)*(rounding*maxima[3, j] + ls*maxima[4, j] + jr*maxima[5, j] + jphase*maxima[6, j] +
            maxima[7, j]) + blockcol[j]
    end
    between = j -> norms[j] > 1 && maxima[2, j] <= roundoff[j]
    if !any(between, eachindex(norms))
        fill!(roundoff, 0.0)
        return nothing
    end
    excess = stepresidualexcess!(bound, sys, a2, a1, residual, delta, x, X, phi, jwork, rhs, work, dwork, rw, rowtol)
    for j in eachindex(norms)
        roundoff[j] = between(j) && excess[j] <= 1 ? norms[j] : 0.0
    end
    return nothing
end

# The cheap bound of the rounding of the blocks' reflected waves per
# column into the host vector `out`: at each stage the largest row sum of
# the output terms' `|S|` at the stage's weights times the largest of the
# stage's stacked states, the larger of the two, in as many units of
# roundoff as the states it sums. Typed on the work, so that the callers,
# which hold it unspecialized, reach it through one dispatch.
function stackedscatterbound!(out, rw::RationalWork, cp::RationalCoupling)
    nz = cp.nstates
    y = rw.ystack
    s1 = sum(j -> abs(rw.weights[j, 1])*cp.SCrows[j], eachindex(cp.SCrows))
    s2 = sum(j -> abs(rw.weights[j, 2])*cp.SCrows[j], eachindex(cp.SCrows))
    weight = nz*eps(Float64)
    if y isa Array
        @inbounds for col in axes(y, 2)
            y1, y2 = 0.0, 0.0
            for k in 1:nz
                y1 = max(y1, abs(y[k, col]))
                y2 = max(y2, abs(y[nz + k, col]))
            end
            out[col] = weight*max(s1*y1, s2*y2)
        end
    else
        out .= weight .* max.(s1 .* vec(Array(maximum(abs, view(y, 1:nz, :); dims = 1, init = 0.0))),
            s2 .* vec(Array(maximum(abs, view(y, nz + 1:2nz, :); dims = 1, init = 0.0))))
    end
    return out
end

# The excess of a step's residual over its rounding, row by row (see
# `residualexcess`), each row's bound from its own terms: the matrices of
# `stagerounding!` on the increments, the stiffness on the state, the
# junction currents `|RJ'| |lmolj relation(phi)|` from `jwork`, which is
# then overwritten, the rounding of the junction phases of the stage
# values `X`, the right hand side, and the blocks' `|S| |y|`, weighted by
# the tolerance `rowtol` of each row. `bound` and `work` are `(n, N, 2)`
# work.
function stepresidualexcess!(bound, sys::TransientSystem, a2, a1, residual, delta, x, X, phi, jwork, rhs, work, dwork,
        @nospecialize(rw::Union{Nothing, AbstractRationalWork}), rowtol)
    dsum, prod = stage(work, 1), stage(work, 2)
    stagerounding!(bound, sys, a2, a1, delta, dsum, prod)
    dsum .= abs.(x)
    stepmul!(prod, sys.Labs, dsum)
    stage(bound, 1) .+= prod
    stage(bound, 2) .+= prod
    if size(jwork, 1) > 0
        jwork .= abs.(jwork)
        for i in 1:2
            stepmul!(prod, sys.RJtabs, stage(jwork, i))
            stage(bound, i) .+= prod
        end
        junctionrounding!(bound, sys, phi, X, dsum, stage(jwork, 1), dwork, prod)
    end
    bound .= 8eps(Float64) .* (bound .+ abs.(rhs))
    if !isnothing(rw)
        reflectedwavesbound!(rw, sys.gauss.coupling)
        bound .+= sys.gauss.coupling.nstates*eps(Float64) .* rw.sbound
    end
    return residualexcess(residual, bound, rowtol)
end


"""
    stageweights!(rw::RationalWork, sys::TransientSystem, t1, t2)

Set the weights of the output terms of the rational blocks at the two
stage times `t1` and `t2` of the step about to be formed, which the
reflected waves at the stages read; the weight of the unconverted
response is one and those of the modulated outputs of the pumped blocks
their modulation at the time. A circuit without a pumped block has only
the first, and its weights never change.
"""
function stageweights!(rw::RationalWork, sys::TransientSystem, t1, t2)
    size(rw.weights, 1) == 1 && return rw
    cp = sys.gauss.coupling
    modulationweights!(view(rw.weights, :, 1), sys.problem, cp, t1)
    modulationweights!(view(rw.weights, :, 2), sys.problem, cp, t2)
    return rw
end

# the correction of a work, built on first use
function stagecorrectionwork!(rw::RationalWork, sys::TransientSystem)
    isnothing(rw.correction) || return rw.correction
    cp = sys.gauss.coupling
    backend = sys.backend
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    n = size(cp.Phost, 2) ÷ 2
    ports = cp.modulated
    np, nw = length(ports), length(cp.W)
    r = 2np
    m = size(rw.z, 2)
    Ucols = allocate(n, np)
    copyto!(Ucols, Matrix(cp.Shost[:, ports]))
    sc = StageCorrection(Ucols, allocate(n, m, 2, r), zeros(np, r, m, nw), LU{Float64,Matrix{Float64},Vector{Int}}[],
        allocate(n, m, 2), allocate(n, m, 2), allocate(np, m), zeros(np, m), zeros(r, m), allocate(m, r), zeros(m, r),
        zeros(r, r))
    rw.correction = sc
    return sc
end

# The pieces of the correction which depend on the factorization, built
# on the current one after a refresh: the solves `K = J*^(-1) U` with the
# products `W_ij K`; then the step's `M`.
function stagefactors!(rw::RationalWork, sys::TransientSystem, bf::GaussBatchFactor, gc::GaussCoefficients, rc, zc)
    sys.gauss.pumped || return nothing
    cp = sys.gauss.coupling
    sc = stagecorrectionwork!(rw, sys)
    np, nw = length(cp.modulated), length(cp.W)
    for col in 1:2np
        i, q = (col - 1) ÷ np + 1, (col - 1) % np + 1
        fill!(sc.rhs, 0)
        view(sc.rhs, :, :, i) .= view(sc.Ucols, :, q)
        gaussbatchtransform!(sc.sol, sc.rhs, bf, gc, rc, zc)
        view(sc.K, :, :, :, col) .= sc.sol
        stack!(rw.dstack, sc.sol)
        for w in 1:nw
            stepmul!(sc.y, cp.W[w], rw.dstack)
            copyto!(sc.ys, sc.y)
            view(sc.Q, :, col, :, w) .= sc.ys
        end
    end
    stagemultipliers!(rw, sys)
    return nothing
end

# The pieces of the correction at the step's weights, from those of the
# factorization: `M` per column. The term `w` of `W` is stage `i` and
# converted output `j`, the frozen operator's weight of which is zero, so
# its weight is the step's own.
function stagemultipliers!(rw::RationalWork, sys::TransientSystem)
    sys.gauss.pumped || return nothing
    sc = rw.correction
    isnothing(sc) && return nothing
    cp = sys.gauss.coupling
    np, nw = length(cp.modulated), length(cp.W)
    nmod = nw ÷ 2
    r = 2np
    weight = w -> rw.weights[(w - 1) % nmod + 2, (w - 1) ÷ nmod + 1]
    empty!(sc.M)
    for c in axes(sc.Q, 3)
        Mc = sc.Mh
        fill!(Mc, 0)
        for k in 1:r
            Mc[k, k] = 1
        end
        for w in 1:nw
            d = weight(w)
            iszero(d) && continue
            i = (w - 1) ÷ nmod + 1
            view(Mc, (i - 1)*np + 1:i*np, :) .-= d .* view(sc.Q, :, :, c, w)
        end
        push!(sc.M, lu(Mc))
    end
    return nothing
end

# the products on the block's rows at both stages of a stage pair `c`,
# `V' c`, into the host copy `sc.yh` as `(r, m)` stage major
function stagerows!(sc::StageCorrection, rw::RationalWork, cp::RationalCoupling, c)
    np, nw = length(cp.modulated), length(cp.W)
    nmod = nw ÷ 2
    stack!(rw.dstack, c)
    fill!(sc.yh, 0)
    for w in 1:nw
        i = (w - 1) ÷ nmod + 1
        d = rw.weights[(w - 1) % nmod + 2, i]
        iszero(d) && continue
        stepmul!(sc.y, cp.W[w], rw.dstack)
        copyto!(sc.ys, sc.y)
        view(sc.yh, (i - 1)*np + 1:i*np, :) .-= d .* sc.ys
    end
    return sc.yh
end

# the correction of a solution `c` of the frozen operator, `(n, m, 2)`,
# to that of the true one
function stagecorrect!(c, rw::RationalWork, sys::TransientSystem)
    sys.gauss.pumped || return c
    sc = rw.correction
    isnothing(sc) && return c
    cp = sys.gauss.coupling
    r = 2length(cp.modulated)
    yh = stagerows!(sc, rw, cp, c)
    for cc in axes(sc.zh, 1)
        view(sc.zh, cc, :) .= sc.M[cc] \ view(yh, :, cc)
    end
    copyto!(sc.z, sc.zh)
    for col in 1:r
        c .-= view(sc.K, :, :, :, col) .* reshape(view(sc.z, :, col), 1, :, 1)
    end
    return c
end

# the same at the endpoint of a step, which its reading and the resting
# waves see
function endweights!(rw::RationalWork, sys::TransientSystem, t)
    size(rw.weights, 1) == 1 && return rw
    modulationweights!(rw.endweights, sys.problem, sys.gauss.coupling, t)
    return rw
end

# the reflected waves at both stages from the stacked states `rw.ystack`,
# scattered onto the blocks' rows into `rw.source`: the terms weighted by
# the stage's weights
function reflectedwaves!(rw::RationalWork, cp::RationalCoupling)
    nz = cp.nstates
    for i in 1:2
        yi = view(rw.ystack, (i - 1)*nz + 1:i*nz, :)
        si = stage(rw.source, i)
        fill!(si, 0)
        for j in eachindex(cp.SC)
            w = rw.weights[j, i]
            iszero(w) && continue
            stepmul!(rw.swork, cp.SC[j], yi)
            si .+= w .* rw.swork
        end
    end
    return rw
end

# The rounding of `reflectedwaves!` per unit of roundoff: the same sum
# with every entry in magnitude, so a scatter row reaches only the states
# it multiplies. `|S| |y|` rather than `norm(S) norm(y)`, which for a fit
# whose poles span decades overstates the sum by as many.
function reflectedwavesbound!(rw::RationalWork, cp::RationalCoupling)
    nz = cp.nstates
    rw.yabs .= abs.(rw.ystack)
    for i in 1:2
        yi = view(rw.yabs, (i - 1)*nz + 1:i*nz, :)
        si = stage(rw.sbound, i)
        fill!(si, 0)
        for j in eachindex(cp.SCabs)
            w = abs(rw.weights[j, i])
            iszero(w) && continue
            stepmul!(rw.swork, cp.SCabs[j], yi)
            si .+= w .* rw.swork
        end
    end
    return rw
end

# the transpose: multipliers `mu` on the blocks' rows at both stages
# carried to the stacked states, into `rw.ystack`
function reflectedwavesbar!(rw::RationalWork, cp::RationalCoupling, mu)
    nz = cp.nstates
    for i in 1:2
        yi = view(rw.ystack, (i - 1)*nz + 1:i*nz, :)
        fill!(yi, 0)
        mui = stage(mu, i)
        for j in eachindex(cp.SCt)
            w = rw.weights[j, i]
            iszero(w) && continue
            stepmul!(rw.zwork2, cp.SCt[j], mui)
            yi .+= w .* rw.zwork2
        end
    end
    return rw
end

# the two stages of a `(n, N, 2)` array stacked into `(2n, N)`, and the
# two halves of a stacked array summed onto a `(n, N)` array
function stack!(s, a)
    n = size(a, 1)
    view(s, 1:n, :) .= stage(a, 1)
    view(s, n + 1:2n, :) .= stage(a, 2)
    return s
end
function stacksum!(y, s)
    n = size(y, 1)
    y .+= view(s, 1:n, :) .+ view(s, n + 1:2n, :)
    return y
end

# the reflected waves of the rational parts at the stages scattered onto
# the blocks' rows, from the stage increments, the stage values and the
# states: the stacked states `P_d delta + P_x X + P_z z` through the
# output terms at the stages' weights
function rationalsources!(rw::RationalWork, sys::TransientSystem, delta, X)
    cp = sys.gauss.coupling
    stack!(rw.dstack, delta)
    stack!(rw.xstack, X)
    stepmul!(rw.ystack, cp.Pd, rw.dstack)
    stepmul!(rw.ywork, cp.Px, rw.xstack)
    rw.ystack .+= rw.ywork
    stepmul!(rw.ywork, cp.Pz, rw.z)
    rw.ystack .+= rw.ywork
    reflectedwaves!(rw, cp)
    return rw
end

# the states of the rational blocks at the end of a step from the states
# at its start, the stage increments and the stage values
function rationalstates!(rw::RationalWork, sys::TransientSystem, delta, X)
    cp = sys.gauss.coupling
    stack!(rw.dstack, delta)
    stack!(rw.xstack, X)
    stepmul!(rw.zwork, cp.Ez, rw.z)
    stepmul!(rw.zwork2, cp.Ed, rw.dstack)
    rw.zwork .+= rw.zwork2
    stepmul!(rw.zwork2, cp.Ex, rw.xstack)
    rw.zwork .+= rw.zwork2
    copyto!(rw.z, rw.zwork)
    return rw
end

# The adjoint's states. The cotangent of the states after a step reaches
# the stage unknowns through the increment and the value parts of the
# update, and the state's flux through the value part, as the value at
# a stage is the flux plus the increment; then the multipliers of the
# step carry the reflected waves they weigh to the states before the
# step, and the value part of those waves' incident waves to the flux as
# well; and the multiplier of the endpoint reading carries the resting
# waves it saw to the states.
function statesbartostages!(rw::RationalWork, sys::TransientSystem, wst, xextra)
    cp = sys.gauss.coupling
    n = size(xextra, 1)
    stepmul!(rw.tstack, cp.Edt, rw.z)
    stage(wst, 1) .+= view(rw.tstack, 1:n, :)
    stage(wst, 2) .+= view(rw.tstack, n + 1:2n, :)
    stepmul!(rw.tstack, cp.Ext, rw.z)
    stage(wst, 1) .+= view(rw.tstack, 1:n, :)
    stage(wst, 2) .+= view(rw.tstack, n + 1:2n, :)
    stacksum!(xextra, rw.tstack)
    return rw
end
function statesbarstep!(rw::RationalWork, sys::TransientSystem, mu, xextra)
    cp = sys.gauss.coupling
    reflectedwavesbar!(rw, cp, mu)
    stepmul!(rw.zwork, cp.Ezt, rw.z)
    stepmul!(rw.zwork2, cp.Pzt, rw.ystack)
    rw.z .= rw.zwork .+ rw.zwork2
    stepmul!(rw.tstack, cp.Pxt, rw.ystack)
    stacksum!(xextra, rw.tstack)
    return rw
end
# the cotangent of the resting waves the endpoint reading saw, at the
# endpoint's weights, to the states
function restingwavesbar!(rw::RationalWork, sys::TransientSystem, work, sign)
    cp = sys.gauss.coupling
    for j in eachindex(cp.SCt)
        w = rw.endweights[j]
        iszero(w) && continue
        stepmul!(rw.zwork2, cp.SCt[j], work)
        rw.z .+= (sign*w) .* rw.zwork2
    end
    return rw
end

# the reflected waves of the rational parts at rest at their states at
# the time `t`, per port on the host, the source the endpoint reading and
# the consistency check see on the blocks' rows
function restingwaves(sys::TransientSystem, z, t)
    cp = sys.gauss.coupling
    w = modulationweights!(ones(length(cp.terms)), sys.problem, cp, t)
    zh = Array(z)
    out = w[1] .* (cp.Cblk[1]*zh)
    for j in 2:length(cp.terms)
        iszero(w[j]) || (out .+= w[j] .* (cp.Cblk[j]*zh))
    end
    return out
end

# the drives of every condition at a time, as the scaled node currents in
# the columns of `b`, and the drive values; one condition's in vectors,
# the trapezoidal and backward Euler rules', whose problem is the solve's
# and which the system holds only in its injection and its scale
function batchdrivecurrent!(b, sys::TransientSystem, problems, values, hostvalues, t)
    drivevalues!(hostvalues, problems, t)
    values === hostvalues || copyto!(values, hostvalues)
    stepmul!(b, sys.injection, values)
    b .+= sys.constant
    return b
end

# the drive current of a stepper at a time, the lines' forced currents
# from their histories included
function batchdrivecurrent!(b, st, t)
    batchdrivecurrent!(b, st.sys, st.problems, st.values, st.hostvalues, t)
    if !isempty(st.sys.problem.lines)
        linevalues!(st, t)
        stepmul!(st.lwork, st.sys.lineinjection, st.linevalues)
        b .+= st.lwork
    end
    return b
end

# The wave arriving at each line port at time `t`, the far port's wave a
# delay earlier read from the history, and the current it forces,
# `2 q / sqrt(Z)`, into the values on the backend. The history is a ring
# of the wave leaving each port at the grid times, `column c` at
# `tpre + (c - 1) h` in `slot mod1(c, nring)`, the first `npre` columns
# the prehistory before and at the start, the sampled history the
# initial state carries in a solve, constant where that history is one
# column, and a perturbation in a response, and the ring long enough
# for the longest delay and the stencil's reach; the columns are
# accepted through `st.accepted`. A read interpolates with the centered
# Lagrange stencil of as many samples as the history holds on both sides
# of the query's interval, up to three on each side, so the quintic
# where the delay leaves three accepted samples past the query, the
# cubic with two and the linear with one. Every one of these is a
# contraction, its magnitude at most one at every frequency, so a wave
# gains nothing on a round trip through any delay of at least a step,
# and a passive line stays passive; a one sided stencil is not, its
# gain reaching 2.4 near a delay of a step, and would let a mismatched
# short line grow without bound. The quintic keeps the rule's order; a
# delay under four steps reads at lower order. A query before the
# history reads its first column. The stencil of every port is the same
# for every condition, so it is set on the host and the read is one
# gather over the ports and the conditions.
function linevalues!(st, t)
    linestencil!(st.stencil, st.sys.problem.lines, t, linestart(st), st.sys.h, st.accepted)
    copyto!(st.dstencil, st.stencil)
    wavegather!(st.linevalues, st.waves, st.dstencil, st.far, st.readscale, st.sys.backend)
    return st.linevalues
end

# the samples of prehistory a response keeps before the start at step
# `h`: enough for every read before the start, the longest delay and the
# stencil's reach, one for a circuit without lines; the reads of one
# line's ports reach the last of them, as many as its own delay and the
# reach take; and the ring that holds a window of them
lineprehistory(line::TransientLine, h) = ceil(Int, line.delay/h) + 6
lineprehistory(p::TransientProblem, h) = isempty(p.lines) ? 1 : maximum(line -> lineprehistory(line, h), p.lines)
# the first column of the prehistory the reads of each line port reach,
# two ports per line
function lineprehistorystarts(p::TransientProblem, h)
    npre = lineprehistory(p, h)
    return [npre - lineprehistory(line, h) + 1 for line in p.lines for _ in 1:2]
end
linering(npre::Int) = npre + 2
# the time of the first column of a history whose column `npre` is at `t0`
linestart(t0, npre, h) = t0 - (npre - 1)*h
linestart(st) = linestart(first(st.times), st.npre, st.sys.h)

# the stencils of every port's read at time `t`, the far port's wave a
# delay earlier, in the columns of `stencil`: the column before the
# first sample read, the number of samples, and their six weights
function linestencil!(stencil, lines, t, tpre, h, accepted)
    for (l, line) in enumerate(lines), q in 1:2
        row = 2(l - 1) + q
        first, nst = linestencil!(view(stencil, 3:8, row), t - line.delay, tpre, h, accepted)
        stencil[1, row] = first
        stencil[2, row] = nst
    end
    return stencil
end

# the stencil of a read at time `s` of a history whose first column is at
# time `tpre` with `accepted` columns: the query's interval between two
# samples, and as many samples on each side of it as the history holds,
# up to three, so the stencil is centered; the weights into `weights`,
# the column before the first read and the count returned; a query before
# the history reads its first column, and one past it its last interval
function linestencil!(weights, s, tpre, h, accepted)
    fill!(weights, 0)
    if s <= tpre || accepted < 2
        weights[1] = 1.0
        return 0, 1
    end
    u = (s - tpre)/h
    # the interval, its left sample `i` counted from zero
    i = min(floor(Int, u), accepted - 2)
    half = min(i + 1, accepted - 1 - i, 3)
    x = u - (i - half + 1)
    nst = 2half
    for k in 1:nst
        w = 1.0
        for m in 0:nst - 1
            m == k - 1 && continue
            w *= (x - m)/(k - 1 - m)
        end
        weights[k] = w
    end
    return i - half + 1, nst
end
# the read of every port's far wave through its stencil, scaled, over
# the ports and the conditions, and its transpose scattering a value
# onto the samples it read, with a sign, accumulating; the columns in
# the slots of the ring. Each port and condition is one item, whose body
# the kernel and the host loop share.
@inline function wavegatheritem!(values, waves, stencil, far, scale, nring, row, j)
    @inbounds begin
        f = far[row]
        first = Int(stencil[1, row])
        nst = Int(stencil[2, row])
        total = 0.0
        for k in 1:nst
            total += stencil[2 + k, row]*waves[f, mod1(first + k, nring), j]
        end
        values[row, j] = scale[row]*total
    end
    return nothing
end
@inline function wavescatteritem!(wavesbar, values, stencil, far, scale, sign, nring, row, j)
    @inbounds begin
        f = far[row]
        first = Int(stencil[1, row])
        nst = Int(stencil[2, row])
        v = sign*scale[row]*values[row, j]
        for k in 1:nst
            wavesbar[f, mod1(first + k, nring), j] += stencil[2 + k, row]*v
        end
    end
    return nothing
end
@kernel function wavegatherkernel!(values, @Const(waves), @Const(stencil), @Const(far), @Const(scale), nring)
    row, j = @index(Global, NTuple)
    wavegatheritem!(values, waves, stencil, far, scale, nring, row, j)
end
function wavegather!(values, waves, stencil, far, scale, backend)
    isempty(values) && return values
    wavegatherkernel!(backend, (64, 1))(values, waves, stencil, far, scale, size(waves, 2); ndrange = size(values))
    return values
end
# on the host a loop: a kernel launch there costs more than the read
function wavegather!(values, waves, stencil, far, scale, ::CPU)
    nring = size(waves, 2)
    for j in axes(values, 2), row in axes(values, 1)
        wavegatheritem!(values, waves, stencil, far, scale, nring, row, j)
    end
    return values
end
@kernel function wavescatterkernel!(wavesbar, @Const(values), @Const(stencil), @Const(far), @Const(scale), sign, nring)
    row, j = @index(Global, NTuple)
    wavescatteritem!(wavesbar, values, stencil, far, scale, sign, nring, row, j)
end
function wavescatter!(wavesbar, values, stencil, far, scale, sign, backend)
    isempty(values) && return wavesbar
    wavescatterkernel!(backend, (64, 1))(wavesbar, values, stencil, far, scale, Float64(sign), size(wavesbar, 2); ndrange = size(values))
    return wavesbar
end
function wavescatter!(wavesbar, values, stencil, far, scale, sign, ::CPU)
    nring = size(wavesbar, 2)
    for j in axes(values, 2), row in axes(values, 1)
        wavescatteritem!(wavesbar, values, stencil, far, scale, sign, nring, row, j)
    end
    return wavesbar
end

# the slot of a column of the ring, and the columns `c1:c2` of the ring
# set from the columns of `tail`, or read into them
ringslot(c, nring) = mod1(c, nring)
# in at most two contiguous segments, where the columns wrap
function ringsegments(nring, c1, count)
    count <= nring || throw(ArgumentError("the history is longer than the ring."))
    s1 = ringslot(c1, nring)
    s1 + count - 1 <= nring && return [(s1:s1 + count - 1, 1:count)]
    first = nring - s1 + 1
    return [(s1:nring, 1:first), (1:count - first, first + 1:count)]
end
function fillhistory!(ring, tail, c1)
    for (slots, columns) in ringsegments(size(ring, 2), c1, size(tail, 2))
        view(ring, :, slots, :) .= view(tail, :, columns, :)
    end
    return ring
end
function readhistory!(tail, ring, c1)
    for (slots, columns) in ringsegments(size(ring, 2), c1, size(tail, 2))
        view(tail, :, columns, :) .= view(ring, :, slots, :)
    end
    return tail
end

# the tables of the line ports on the backend: each port's far port, the
# read's scale `2/sqrt(Z)`, and `sqrt(Z)`
function linetables(p::TransientProblem, backend)
    far = [2(l - 1) + (3 - q) for l in eachindex(p.lines) for q in 1:2]
    sqrtz = [sqrt(line.Z) for line in p.lines for q in 1:2]
    return tobackend(backend, far), tobackend(backend, 2 ./ sqrtz), tobackend(backend, sqrtz)
end

# the waves leaving the line ports at the endpoint just accepted, from
# the rates across the ports and the waves that arrived, into column
# `c` of the history
function linewaves!(st, c)
    stepmul!(st.linerates, st.sys.linegather, st.v)
    view(st.waves, :, ringslot(c, size(st.waves, 2)), :) .= phi0 .* st.linerates ./ st.sqrtz .- st.linevalues .* st.sqrtz ./ 2
    return st.waves
end

# the drive values of a stepper's conditions at a time on its backend,
# which the port waves of a saved time read
function stepperdrives!(st, t)
    drivevalues!(st.hostvalues, st.problems, t)
    copyto!(st.values, st.hostvalues)
    return st.values
end

# the currents of every condition's drives at a time, on the host, and
# the balance of the floating subnetworks a condition's drives feed, read
# in the same pass over the conditions
function drivevalues!(hostvalues, problems, t)
    balanced = true
    @inbounds for (j, p) in enumerate(problems)
        for (k, d) in enumerate(p.drives)
            hostvalues[k, j] = d.current(t)
        end
        balanced &= isempty(p.balance.islands)
    end
    all(isfinite, hostvalues) || throw(ArgumentError(lazy"a source returned a nonfinite current at t = $(t) s."))
    if !balanced
        for (j, p) in enumerate(problems)
            checkbalance(p, view(hostvalues, :, j), t)
        end
    end
    return hostvalues
end

# The stepper of a batch: the buffers of a Gauss-Legendre step, the
# batch's factorizations, the Newton engine's closures, which hold the
# residual's own work, and the state, so that the solve advances it step
# by step and a response replays a window from a checkpoint with the same
# code.
mutable struct GaussStepper{S, P, A, M, C, F, W, B, R, T, E, U, K, HW, FI, NW}
    sys::S
    problems::P
    N::Int
    x::M
    v::M
    X::A
    delta::A
    lastdelta::A
    residual::A
    trial::A
    trialresidual::A
    correction::A
    rhs::A
    xnew::M
    cv::M
    lx::M
    b1::M
    b2::M
    phi::A
    # the derivative of the relation at one stage, for a circuit which has
    # a polynomial one, and empty for the Josephson relation, which is
    # broadcast in place
    dwork::M
    cosphi::C
    rc::C
    zc::C
    hostvalues::Matrix{Float64}
    values::F
    portwork::M
    bf::B
    baseresidual!::R
    trialresidual!::T
    refresh!::E
    solve!::U
    accept!::K
    # the projection's work, which keeps each condition's factorization
    # of its Jacobian across the steps, or nothing, and the rational
    # blocks' work, or nothing, as unions within the stepper's array types
    # for the same reason the stage's are
    pw::Union{Nothing, ProjectionWork{M, Matrix{Float64}, M}}
    # the lines: the ring of the waves leaving their ports on the
    # backend, `(port, slot, condition)`, the prehistory columns, the
    # grid times, the columns accepted so far, the currents the arriving
    # waves force, the stencils of a read on the host and the backend,
    # the ports' tables, and the rates across the ports
    waves::HW
    npre::Int
    times::Vector{Float64}
    accepted::Int
    linevalues::M
    stencil::Matrix{Float64}
    dstencil::M
    far::FI
    readscale::W
    sqrtz::W
    linerates::M
    lwork::M
    rw::Union{Nothing, RationalWork{M, A}}
    newtonwork::NW
    # the tolerance of every row of every stage, set by the step, which the
    # residual's norm is weighted by, and that of the weighted norm of each
    # column, one
    rowtol::A
    tolerance::Vector{Float64}
    roundoff::Tuple{Vector{Float64}, Vector{Float64}}
    rtol::Float64
    atol::Float64
    iterations::Int
    # per condition, whether the last step refreshed its factorization,
    # which the next then refreshes before its first correction, the
    # refresh the step began with included (see `newtonsolve!`), and which
    # `setstate!` clears
    stalefailed::Vector{Bool}
    # the indices of the conditions in the batch, which a failed step
    # names
    conditions::Vector{Int}
    # the corrections, and the factorizations of every condition at once
    # (those of each condition's own refreshes are the Newton work's)
    corrections::Int
    factorizations::Int
end

function gaussstepper(sys::TransientSystem, problems, rtol, atol, iterations, bf; conditions = 1:length(problems))
    p = sys.problem
    backend = sys.backend
    gc = sys.gauss.coefficients
    N = length(problems)
    n, np, nd = length(p), length(p.portimpedances), length(p.drives)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    x, v = allocate(n, N), allocate(n, N)
    X, delta, lastdelta, residual, trial, trialresidual, correction, rhs = [allocate(n, N, 2) for _ in 1:8]
    junction, trialjunction, cwork, gwork = [allocate(n, N, 2) for _ in 1:4]
    xnew, cv, lx, b1, b2 = [allocate(n, N) for _ in 1:5]
    nj = length(sys.lmolj)
    phi, trialphi, jwork = [allocate(nj, N, 2) for _ in 1:3]
    dwork = allsinusoidal(sys.relations) ? allocate(nj, 0) : allocate(nj, N)
    cosphi = KernelAbstractions.zeros(backend, ComplexF64, nj, N)
    rc, zc = [KernelAbstractions.zeros(backend, ComplexF64, n, N) for _ in 1:2]
    roundoff = (zeros(N), zeros(N))
    hostvalues = zeros(nd, N)
    values = tobackend(backend, zeros(nd, N))
    portwork = allocate(np, N)
    rw = isnothing(sys.gauss.coupling) ? nothing : rationalwork(p, backend, n, N)
    cell = RationalWorkCell{typeof(x), typeof(X)}(rw)
    # the tolerance of each row, set by the step, the columns' maxima of
    # a residual on the host, reduced first on a device, and the blocks'
    # part of the floor on the host
    rowtol = allocate(n, N, 2)
    maxima, devicemaxima = zeros(COLUMNMAXIMA, N), backend isa CPU ? nothing : allocate(COLUMNMAXIMA, N)
    blockcol = zeros(N)
    baseresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, phi, junction, jwork, cwork, gwork, rhs,
        roundoff[1], maxima, devicemaxima, blockcol, rowtol, dwork, cell.work)
    trialresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, trialphi, trialjunction, jwork, cwork,
        gwork, rhs, roundoff[2], maxima, devicemaxima, blockcol, rowtol, dwork, cell.work)
    # a new factorization of the frozen operator of the columns which
    # asked, and the stage correction of a pumped block rebuilt on it
    refresh! = columns -> (gaussbatchjacobian!(bf, sys, phi, cosphi, dwork, columns);
        isnothing(cell.work) || stagefactors!(cell.work, sys, bf, gc, rc, zc); nothing)
    solve! = (c, r) -> (gaussbatchtransform!(c, r, bf, gc, rc, zc); isnothing(cell.work) || stagecorrect!(c, cell.work, sys); (false, 0))
    accept! = mask -> (maskcolumns!(phi, trialphi, mask); maskcolumns!(junction, trialjunction, mask); nothing)
    pr = sys.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, N; kept = true)
    nl = 2length(p.lines)
    far, readscale, sqrtz = linetables(p, backend)
    npre = lineprehistory(p, sys.h)
    return GaussStepper(sys, problems, N, x, v, X, delta, lastdelta, residual, trial, trialresidual, correction, rhs,
        xnew, cv, lx, b1, b2, phi, dwork, cosphi, rc, zc,
        hostvalues, values, portwork, bf, baseresidual!, trialresidual!, refresh!, solve!, accept!, pw,
        allocate(nl, linering(npre), N), npre, Float64[], 0, allocate(nl, N), zeros(8, nl), allocate(8, nl), far, readscale, sqrtz,
        allocate(nl, N), allocate(n, N), rw,
        NewtonWork(backend, N), rowtol, ones(N), roundoff, Float64(rtol), Float64(atol), Int(iterations), fill(false, N),
        collect(Int, conditions), 0, 0)
end

# the stepper given its grid and the `npre` columns of history before
# the start of step `k + 1`, the columns `k + 1` to `npre + k`, the last
# at the step's start, accepted through them
function sethistory!(st::GaussStepper, times, tail, k)
    st.times = times
    fillhistory!(st.waves, tail, k + 1)
    st.accepted = st.npre + k
    return st
end

# the rate of a stepper's state at time `t` of a record starting at `t0`
# read along the algebraic directions into `dst`, with the work `rw` of
# `ratereadwork`: the drives' rate by the difference of `delta` (see
# `ratestencil`), the lines' forced currents' rate by a central
# difference of their histories, which hold the waves before the start
function readstepper!(dst, st::GaussStepper, @nospecialize(rw::Union{Nothing, RateReadWork}), t, t0, delta)
    sys = st.sys
    if !isnothing(sys.projection)
        linerate = if !isempty(sys.problem.lines)
            lr = Array(linevalues!(st, t + delta))
            lr .= (lr .- Array(linevalues!(st, t - delta))) ./ (2delta)
        else
            nothing
        end
        drivedotz!(rw.bdotz, sys.projection, st.problems, t, t0, delta, rw.hv1, rw.hv2, linerate)
    end
    return readrate!(dst, st.v, st.x, sys, rw)
end

# the phases of the projected junctions at a state, for the record
function projectedphases!(st::GaussStepper, x)
    pr = st.sys.projection
    isempty(pr.pj) || stepmul!(st.pw.phip, pr.RJp, x)
    return st.pw.phip
end

# The drive of the constraints' rows at a time, `gb`, and the magnitudes
# of its terms, `gbabs`, which set the residual's floor: from the sources'
# values `drive` and their magnitudes `driveabs`, the currents the lines
# force, `linedrive`, and the blocks' resting waves, `resting`, either of
# them nothing where the circuit has none. The projection of an endpoint
# under every rule, the Gauss-Legendre stepper's and the trapezoidal and
# backward Euler step's, and the start's (see `projectedstate`) take it
# alike, compiled once whichever of the two terms a circuit has.
function constraintdrive!(gb, gbabs, pr::ConstraintProjection, drive, driveabs, @nospecialize(linedrive),
        @nospecialize(resting))
    mul!(gb, pr.Zcinj, drive)
    gb .+= pr.Zcconstant
    mul!(gbabs, pr.Zcinjabs, driveabs)
    gbabs .+= abs.(pr.Zcconstant)
    if !isnothing(linedrive)
        mul!(gb, pr.Zcline, linedrive, 1.0, 1.0)
        mul!(gbabs, pr.Zclineabs, abs.(linedrive), 1.0, 1.0)
    end
    if !isnothing(resting)
        mul!(gb, pr.Zcblock, resting, 1.0, 1.0)
        mul!(gbabs, pr.Zcblockabs, abs.(resting), 1.0, 1.0)
    end
    return gb, gbabs
end

# the drive of the index one rows at a time, `gq`, from the same terms
function readdrive!(gq, pr::ConstraintProjection, drive, @nospecialize(linedrive), @nospecialize(resting))
    mul!(gq, pr.Qinj, drive)
    gq .+= pr.Qconst
    isnothing(linedrive) || mul!(gq, pr.Qline, linedrive, 1.0, 1.0)
    isnothing(resting) || mul!(gq, pr.Qblock, resting, 1.0, 1.0)
    return gq
end

# The projection of every condition's endpoint onto the algebraic
# constraints along the projected directions (see `projectendpoint!`):
# Newton on the coefficients `alpha` of `Z` per condition, the residual
# `Z' (L x + J(x) - b(t))` and its Jacobians `Z' (L + J'(x)) Z` on the
# host from two small products on the backend, to the roundoff of the
# terms the constraint balances, the predictor moved with the state. The
# step's tolerance is set in the stage equations' units and is far
# looser for the constraint, so the drift of the endpoints would
# accumulate below it. The rate along `Z` is left to the reading at the
# read-outs, except where a block is on a direction, when it is the
# derivative of the cubic through the state, the stages and the projected
# endpoint. Leaves the projected junctions' phases at the endpoint in the
# work buffer.
function gaussproject!(st::GaussStepper, t, step)
    sys, pw = st.sys, st.pw
    pr = sys.projection
    gc = sys.gauss.coefficients
    h = sys.h
    drivevalues!(st.hostvalues, st.problems, t)
    isempty(sys.problem.lines) || linevalues!(st, t)
    linedrive = isempty(sys.problem.lines) ? nothing : Array(st.linevalues)
    resting = isnothing(st.rw) ? nothing : restingwaves(sys, st.rw.z, t)
    if !isempty(pr.directions)
        # the drive along the constraints' rows, and the magnitudes of its
        # terms for the residual's floor
        pw.hostabs .= abs.(st.hostvalues)
        gb, gbabs = constraintdrive!(pw.gb, pw.gbabs, pr, st.hostvalues, pw.hostabs, linedrive, resting)
        # the predictor moves with the state
        moved! = w -> (stage(st.lastdelta, 1) .-= w; stage(st.lastdelta, 2) .-= w; nothing)
        projectendpoint!(st.xnew, pw, pr, gb, gbabs, st.iterations, moved!) ||
            throw(TransientStepError(step, t, st.conditions[findall(pw.open)], :projection))
        pr.cubic && projectrate!(st.v, pr, pw, gc, h, stage(st.delta, 1), stage(st.delta, 2), st.xnew, st.x)
    end
    # the index one unknowns of the endpoint from their equations
    if !isempty(pr.readrows) || !isempty(pr.auxrows)
        projectedphases!(st, st.xnew)
        copyto!(pw.hphi, pw.phip)
        gq = readdrive!(pw.gqdrive, pr, st.hostvalues, linedrive, resting)
        relationinto!(pw.current, pr.relationsp, pw.hphi)
        pw.current .*= pr.lmoljp
        endpointread!(st.v, st.xnew, pr, pw, pw.current, gq)
        projectedphases!(st, st.xnew)
    end
    return st
end

# A state made consistent with the circuit's algebraic equations at the
# time `t` as the Gauss-Legendre stepper makes every endpoint (see
# `gaussproject!`): the flux projected onto the constraints of the
# directions a junction or a source touches, the index one unknowns read
# from their equations, and the rate along every algebraic direction read
# from the differentiated constraints as at the start of a record, at
# most `iterations` corrections of the projection. The lines force the
# currents the state's history says arrive at `t`, with the rate of their
# central difference, as the stepper reads them at an endpoint; the
# history is the state's. A harmonic balance orbit meets the constraints
# along a junction's node without capacitance only to the truncation of
# its harmonics, and the transient's check of its start holds them to the
# solve's tolerance.
function projectedstate(state::TransientState, sys::TransientSystem, t::Real; iterations::Int = 15)
    p = sys.problem
    n = length(p)
    x, v = copy(state.flux), copy(state.rate)
    X, V = reshape(x, n, 1), reshape(v, n, 1)
    pr = sys.projection
    delta = ratedelta(sys)
    # the currents the lines force at `t`, the state's history ending
    # there, and their rate
    lines = !isempty(p.lines)
    forced, forcedrate = stateforcing(state, p, delta)
    linedrive = lines ? reshape(forced, :, 1) : nothing
    linerate = lines ? reshape(forcedrate, :, 1) : nothing
    if !isnothing(pr)
        pw = projectionwork(pr, CPU(), n, 1, 1)
        drive = drivevalues!(zeros(length(p.drives), 1), [p], t)
        resting = isnothing(sys.gauss.coupling) ? nothing :
            restingwaves(sys, reshape(initialblockstates(state, p), :, 1), t)
        if !isempty(pr.directions)
            gb, gbabs = constraintdrive!(pw.gb, pw.gbabs, pr, drive, abs.(drive), linedrive, resting)
            projectendpoint!(X, pw, pr, gb, gbabs, iterations) || throw(ArgumentError(
                lazy"the state cannot be projected onto the algebraic constraints at t = $(t) s in $(iterations) corrections."))
        end
        if !isempty(pr.readrows) || !isempty(pr.auxrows)
            isempty(pr.pj) || (stepmul!(pw.phip, pr.RJp, X); copyto!(pw.hphi, pw.phip))
            gq = readdrive!(pw.gqdrive, pr, drive, linedrive, resting)
            relationinto!(pw.current, pr.relationsp, pw.hphi)
            pw.current .*= pr.lmoljp
            endpointread!(V, X, pr, pw, pw.current, gq)
        end
    end
    rw = ratereadwork(sys, CPU(), n, 1, 1)
    isnothing(pr) || drivedotz!(rw.bdotz, pr, [p], t, t, delta, rw.hv1, rw.hv2, linerate)
    readrate!(V, V, X, sys, rw)
    return TransientState(x, v, state.waves, state.wavesdt, state.blockstates)
end

# the stepper set to a state and the stage increments to start from, its
# stage matrices and its projection's Jacobians assembled at that state,
# factorized anew with `fresh`
function setstate!(st::GaussStepper, x, v, lastdelta; fresh::Bool = false)
    copyto!(st.x, x)
    copyto!(st.v, v)
    isnothing(lastdelta) ? fill!(st.lastdelta, 0) : copyto!(st.lastdelta, lastdelta)
    st.X .= st.x
    stepmul!(stage(st.phi, 1), st.sys.RJ, stage(st.X, 1)); stepmul!(stage(st.phi, 2), st.sys.RJ, stage(st.X, 2))
    gaussbatchjacobian!(st.bf, st.sys, st.phi, st.cosphi, st.dwork; fresh)
    isnothing(st.rw) || stagefactors!(st.rw, st.sys, st.bf, st.sys.gauss.coefficients, st.rc, st.zc)
    isnothing(st.pw) || refreshprojection!(st.pw, st.sys.projection, st.x; fresh)
    st.factorizations += 1
    fill!(st.stalefailed, false)
    return st
end

# one step of every condition from `tprev` to `t`: the drives at the stage
# times, the Newton solve of the stage increments from the predictor, the
# endpoint, and the next predictor; the phases of the stages are left in
# `st.phi`
function advance!(st::GaussStepper, tprev, t, step)
    sys, gc = st.sys, st.sys.gauss.coefficients
    h = sys.h
    st.accepted = st.npre + step - 1
    if !isnothing(st.rw)
        stageweights!(st.rw, sys, tprev + gc.c[1]*h, tprev + gc.c[2]*h)
        endweights!(st.rw, sys, t)
        # a pumped block's stage correction at this step's weights
        stagemultipliers!(st.rw, sys)
    end
    batchdrivecurrent!(st.b1, st, tprev + gc.c[1]*h)
    batchdrivecurrent!(st.b2, st, tprev + gc.c[2]*h)
    stepmul!(st.cv, sys.C, st.v)
    stage(st.rhs, 1) .= st.b1 .+ (gc.ainvone[1]/h) .* st.cv
    stage(st.rhs, 2) .= st.b2 .+ (gc.ainvone[2]/h) .* st.cv
    stepmul!(st.lx, sys.L, st.x)
    # Each row's tolerance relative to the magnitudes of its own terms at
    # the step's start, with `atol` alone absolute: its right hand side,
    # the capacitive rate term and the stiffness on the state, each summed
    # in magnitude, so that a row whose terms cancel, a node coupled to a
    # moving one or a coupled inductor's current, is held to them rather
    # than to their small net. A weak drive is then converged relative to
    # itself, whatever a bias held on another row carries, and a zero
    # right hand side stops on the roundoff floor. The stage drives' work
    # is free once the right hand side holds them.
    st.xnew .= abs.(st.v)
    stepmul!(st.b1, sys.Cabs, st.xnew)
    st.xnew .= abs.(st.x)
    stepmul!(st.b2, sys.Labs, st.xnew)
    for i in 1:2
        stage(st.rowtol, i) .= st.atol .+ st.rtol .* max.(abs.(stage(st.rhs, i)), (abs(gc.ainvone[i])/h) .* st.b1, st.b2)
    end
    copyto!(st.delta, st.lastdelta)
    converged, ncorr, _, _ = newtonsolve!(st.delta, st.correction, st.trial,
        st.residual, st.trialresidual, st.baseresidual!, st.trialresidual!, st.refresh!, st.solve!,
        st.tolerance, st.iterations, st.stalefailed, length(sys.lmolj) > 0, false, st.newtonwork;
        accept! = st.accept!, roundoff = st.roundoff)
    st.corrections += ncorr
    converged || throw(TransientStepError(step, t, st.conditions[unconverged(st.newtonwork, st.tolerance, st.roundoff[1])], :newton))
    st.stalefailed .= st.newtonwork.fresh
    # the rational blocks' states at the end of the step
    isnothing(st.rw) || rationalstates!(st.rw, sys, st.delta, st.X)
    d1, d2 = stage(st.delta, 1), stage(st.delta, 2)
    st.xnew .= st.x .+ gc.ex[1] .* d1 .+ gc.ex[2] .* d2
    st.v .+= (gc.ev[1]/h) .* d1 .+ (gc.ev[2]/h) .* d2
    # the next step's stages start on this step's collocation polynomial
    # extrapolated
    stage(st.lastdelta, 1) .= entry(gc.predict, 1, 1) .* d1 .+ entry(gc.predict, 1, 2) .* d2
    stage(st.lastdelta, 2) .= entry(gc.predict, 2, 1) .* d1 .+ entry(gc.predict, 2, 2) .* d2
    # the endpoint onto the algebraic constraints, and the waves leaving
    # the line ports into the history
    isnothing(st.pw) || gaussproject!(st, t, step)
    if !isempty(sys.problem.lines)
        isnothing(st.pw) && linevalues!(st, t)
        linewaves!(st, st.npre + step)
        st.accepted = st.npre + step
    end
    copyto!(st.x, st.xnew)
    return st
end

# The conditions each task of a host batch steps. The conditions of a
# batch are independent of one another, so a host splits them across the
# threads of the session and steps the chunks at once: the whole step
# parallelizes that way, not only the assembly, the factorization and the
# solve, and the arrays of the batch are filled in place through views, so
# nothing is copied to put the batch back together. A device keeps the one
# chunk it has always had: its batch is already the parallelism, and
# driving one device from several tasks at once is not how the package
# uses it.
function batchchunks(backend, N::Integer)
    (backend isa CPU && N > 1) || return [1:N]
    nc = min(Base.Threads.nthreads(), N)
    nc > 1 || return [1:N]
    return [round(Int, (k - 1)*N/nc) + 1:round(Int, k*N/nc) for k in 1:nc]
end

# The arrays a batch of `N` conditions fills, allocated before the
# conditions are split. Every one carries the condition last, so the
# arrays of a chunk are views along that axis and the batch's solution is
# the whole set. Every array is typed by its rank, and one the record
# does not keep is there with no saved times, so that the outputs of
# every record level are one type, the views of them a chunk steps into
# another, and the stepping compiles once for each with every save
# typed; the solution is built with `nothing` for the arrays not kept
# (see `recordedarray`). The checkpoints are the state and the stage
# predictor every `every` steps, `every` zero without any.
struct GaussCheckpoints{A3, A4}
    every::Int
    flux::A3
    rate::A3
    increment::A4
    states::A3
    waves::A4
end
struct GaussBatchOutputs{A2, A3, A4}
    voltage::A3
    incident::A3
    outgoing::A3
    phases::A4
    endphases::A3
    endrates::A3
    linewaves::A3
    history::A3
    flux::A3
    rate::A3
    stages::A4
    blockstates::A3
    checkpoints::GaussCheckpoints{A3, A4}
    # the waves leaving the line ports over the delay window before the
    # end, and the block states at the end
    finalwaves::A3
    finalstates::A2
    initialflux::A2
    initialrate::A2
    finalflux::A2
    finalrate::A2
    # the initial waves and block states on the host, written once
    initialwaves::Any
    initialstates::Any
end
# whether an array of the outputs is kept by the record: it then has its
# saved times, the dimension before the conditions
recorded(a) = size(a, ndims(a) - 1) > 0
recordedarray(a) = recorded(a) ? a : nothing

function gaussbatchoutputs(sys::TransientSystem, N, nsteps, saveevery, record, checkpointevery)
    savephases, savestates, savecheckpoints = recordlevel(record, saveevery)
    p = sys.problem
    backend = sys.backend
    n, np, nj = length(p), length(p.portimpedances), length(sys.lmolj)
    nl = 2length(p.lines)
    npre = lineprehistory(p, sys.h)
    nzs = blockstates(p)
    pr = sys.projection
    npj = isnothing(pr) ? 0 : length(pr.pj)
    nsaved = cld(nsteps, saveevery) + 1
    # the checkpoints, and the history of the line waves before each of
    # them when that history would not outweigh the record of the waves
    K = checkpointevery > 0 ? Int(checkpointevery) : max(1, round(Int, sqrt(nsteps)))
    nc = savecheckpoints ? cld(nsteps, K) : 0
    tails = savecheckpoints && npre*nc <= nsteps + 1
    zeroed = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    voltage, incident, outgoing = [KernelAbstractions.allocate(backend, Float64, np, nsaved, N) for _ in 1:3]
    kept = keep -> keep ? nsaved : 0
    checkpoints = GaussCheckpoints(savecheckpoints ? K : 0, zeroed(n, nc, N), zeroed(n, nc, N), zeroed(n, 2, nc, N),
        zeroed(nzs, nc, N), zeroed(nl, nl > 0 && tails ? npre : 0, tails ? nc : 0, N))
    return GaussBatchOutputs(voltage, incident, outgoing,
        zeroed(nj, 2, kept(savephases), N),
        zeroed(npj, kept(savephases && npj > 0), N),
        zeroed(npj, kept(savephases && npj > 0), N),
        zeroed(nl, nl > 0 && (savephases || (savecheckpoints && !tails)) ? nsteps + 1 : 0, N),
        zeroed(nl, nl > 0 ? npre : 0, N),
        KernelAbstractions.allocate(backend, Float64, n, kept(savestates), N),
        KernelAbstractions.allocate(backend, Float64, n, kept(savestates), N),
        zeroed(n, 2, kept(savestates), N),
        zeroed(nzs, kept(savestates && nzs > 0), N),
        checkpoints, zeroed(nl, nl > 0 ? npre : 0, N), zeroed(nzs, N),
        zeroed(n, N), zeroed(n, N), zeroed(n, N), zeroed(n, N), zeros(nl, N), zeros(nzs, N))
end

# the arrays of the conditions `ch`, views of the batch's along the
# condition axis, which a chunk fills as though they were its own
function chunkoutputs(out::GaussBatchOutputs, ch)
    slice = a -> selectdim(a, ndims(a), ch)
    checkpoints = GaussCheckpoints(values(mapcheckpoints(slice, out.checkpoints))...)
    return GaussBatchOutputs(slice(out.voltage), slice(out.incident), slice(out.outgoing),
        slice(out.phases), slice(out.endphases), slice(out.endrates), slice(out.linewaves), slice(out.history),
        slice(out.flux), slice(out.rate), slice(out.stages), slice(out.blockstates), checkpoints,
        slice(out.finalwaves), slice(out.finalstates), slice(out.initialflux), slice(out.initialrate), slice(out.finalflux), slice(out.finalrate),
        slice(out.initialwaves), slice(out.initialstates))
end

# the checkpoints of a batch's solution: the record's named tuple, or
# nothing without any
recordedcheckpoints(cp::GaussCheckpoints) = cp.every > 0 ? mapcheckpoints(identity, cp) : nothing

# the checkpoints of a record, a named tuple or the outputs', with `f`
# applied to each of their arrays, as a named tuple, or nothing without
# any
mapcheckpoints(f, cp) = isnothing(cp) ? nothing : (; every = cp.every, flux = f(cp.flux), rate = f(cp.rate),
    increment = f(cp.increment), states = f(cp.states), waves = f(cp.waves))

# The factorizations of every chunk, kept between the solves of one
# batch: one per chunk, taken over when the chunks are laid out as the
# kept ones are, since a factor holds the analyses of its conditions and
# an ordering is not cheap to find again.
function batchfactors(sys::TransientSystem, chunks, reuse)
    kept = isnothing(reuse) ? nothing : reuse.factor
    if kept isa Vector && length(kept) == length(chunks) && all(k -> kept[k] isa GaussBatchFactor &&
            kept[k].ncolumns == length(chunks[k]), eachindex(chunks))
        return kept
    end
    return [gaussbatchfactor(sys, length(ch)) for ch in chunks]
end

# The statistics of a batch stepped in chunks. Every chunk walks the same
# grid; the factorizations and the retries of a chunk are those of its
# condition which needed the most, and its corrections those of its
# joint iterations, so across chunks each is the largest; summing them
# would count one grid several times over.
function mergebatchstats(each)
    s = first(each)
    return (; steps = s.steps, newtoncorrections = maximum(c -> c.newtoncorrections, each),
        factorizations = maximum(c -> c.factorizations, each), retries = maximum(c -> c.retries, each),
        kryloviterations = 0, rtol = s.rtol, atol = s.atol, iterations = s.iterations)
end

# The integration of a batch under the Gauss-Legendre rule: the arrays of
# every condition allocated once, the conditions split into chunks, and
# each chunk stepped into its own views of them. One chunk is the whole
# batch on one task, which is what a device and a single threaded session
# do.
function gaussbatchintegrate(sys::TransientSystem, problems, t0, tf, nsteps, initialstates,
        saveevery, record, checkpointevery, rtol, atol, iterations, reuse;
        chunks = batchchunks(sys.backend, length(problems)))
    h = sys.h
    times = [k == nsteps ? tf : t0 + k*h for k in 0:nsteps]
    savedtimes = [times[1]; [times[step + 1] for step in 1:nsteps if step % saveevery == 0 || step == nsteps]]
    out = gaussbatchoutputs(sys, length(problems), nsteps, saveevery, record, checkpointevery)
    init = batchinitial(sys, problems, initialstates, t0, rtol, atol)
    factors = batchfactors(sys, chunks, reuse)
    # one chunk steps into the same views of the arrays as the chunks of
    # a threaded batch do, so that the stepping compiles for one form of
    # its outputs
    # the chunks share the system, which they read and do not write, and
    # nothing else: their steppers, factors and arrays are their own
    each = runchunks(chunks) do k, ch
        gaussbatchrun!(chunkoutputs(out, ch), sys, problems[ch], chunkinitial(init, ch), times, saveevery, rtol, atol,
            iterations, factors[k], ch)
    end
    stats = mergebatchstats(each)
    isnothing(reuse) || (reuse.factor = factors)
    arrays = (; out.voltage, out.incident, out.outgoing, phases = recordedarray(out.phases),
        endphases = recordedarray(out.endphases), endrates = recordedarray(out.endrates),
        linewaves = recordedarray(out.linewaves), history = recordedarray(out.history), flux = recordedarray(out.flux),
        rate = recordedarray(out.rate), stages = recordedarray(out.stages), checkpoints = recordedcheckpoints(out.checkpoints),
        finalwaves = isempty(sys.problem.lines) ? nothing : out.finalwaves,
        finalstates = isnothing(sys.gauss.coupling) ? nothing : out.finalstates, out.initialflux, out.initialrate,
        out.finalflux, out.finalrate, blockstates = recordedarray(out.blockstates), out.initialwaves, out.initialstates)
    return TransientBatchSolution(problems, sys.method, h, savedtimes, arrays, stats)
end

# The initial states of every condition of a batch, read on the host and
# checked against the algebraic equations before any chunk steps, so that
# a state a condition refuses is refused from the caller's task before a
# thread has begun stepping another: the fluxes and the rates, the
# history of the line waves before the start of every condition and its
# last column, the waves at the start, and the block states. The
# currents the lines force at the start and their rates, read from each
# state's own history at its own step as a stepper reads them (see
# `stateforcing`), are what the check reads.
function batchinitial(sys::TransientSystem, problems, initialstates, t0, rtol, atol)
    p = sys.problem
    N = length(problems)
    n = length(p)
    h = sys.h
    nl = 2length(p.lines)
    npre = lineprehistory(p, h)
    x0, v0, w0, tailh = zeros(n, N), zeros(n, N), zeros(nl, N), zeros(nl, npre, N)
    q0, qdot = zeros(nl, N), zeros(nl, N)
    for (j, state) in enumerate(initialstates)
        xj, vj = state.flux, state.rate
        (length(xj) == n && length(vj) == n) || throw(DimensionMismatch(lazy"the state has $(n) entries; use transientstate."))
        (all(isfinite, xj) && all(isfinite, vj)) || throw(ArgumentError("the initial state must be finite."))
        x0[:, j] .= xj
        v0[:, j] .= vj
        wj = initialwaves(state, p, h, npre)
        all(isfinite, wj) || throw(ArgumentError("the initial line waves must be finite."))
        tailh[:, :, j] .= wj
        w0[:, j] .= view(wj, :, npre)
        if nl > 0
            forced, forcedrate = stateforcing(state, p, ratedelta(sys))
            q0[:, j] .= forced
            qdot[:, j] .= forcedrate
        end
    end
    z0 = zeros(blockstates(p), N)
    withblocks = !isnothing(sys.gauss.coupling)
    if withblocks
        for (j, state) in enumerate(initialstates)
            zj = initialblockstates(state, p)
            length(zj) == size(z0, 1) && all(isfinite, zj) || throw(DimensionMismatch(lazy"the state needs $(size(z0, 1)) finite block states; use transientstate."))
            z0[:, j] .= zj
        end
    end
    # the check of the algebraic equations at the start, per condition
    # under its own drives, the lines forcing what their history says
    # arrives at the start
    resting = withblocks ? restingwaves(sys, z0, t0) : nothing
    for j in 1:N
        c = transientconsistency(sys, view(x0, :, j), view(v0, :, j), t0, problems[j],
            q0[:, j], isnothing(resting) ? nothing : resting[:, j], qdot[:, j])
        consistentstate(c, atol, rtol, h) || throw(ArgumentError(
            lazy"the initial state of condition $(j) violates the algebraic equations of the circuit along a direction without capacitance to ground (a node no capacitor touches, a capacitive island, a coupled inductor or gauge row); supply a consistent transientstate, or start the drive from an equilibrium."))
    end
    return (; x0, v0, w0, tailh, z0)
end

# the initial arrays of the conditions `ch` of a batch, copied so that a
# chunk moves plain host arrays to its backend
chunkinitial(init, ch) = (; x0 = init.x0[:, ch], v0 = init.v0[:, ch], w0 = init.w0[:, ch],
    tailh = init.tailh[:, :, ch], z0 = init.z0[:, ch])

# the error the chunks' tasks met, out of the failure of the tasks the
# threaded batch reports it as: the failed steps of the chunks as the
# first of them, with the conditions of every chunk which failed there,
# or else the first other error
function chunkerror(err::CompositeException)
    isempty(err.exceptions) && return err
    errs = [e for e in (chunkerror(e) for e in err.exceptions) if !(e isa ChunkAborted)]
    isempty(errs) && return err
    all(e -> e isa TransientStepError, errs) || return first(e for e in errs if !(e isa TransientStepError))
    step = minimum(e -> e.step, errs)
    at = [e for e in errs if e.step == step]
    return TransientStepError(step, first(at).time, sort!(reduce(vcat, [e.conditions for e in at])), first(at).cause)
end
chunkerror(err::TaskFailedException) = chunkerror(err.task.result)
chunkerror(err) = err

# The chunks of one `runchunks` share a flag, the step past which none
# goes on, kept in each chunk's task storage so that the stepping loops
# read it without an argument (`checkchunks`). A chunk whose step fails
# lowers it to that step: the others step on up to it, where a failure
# of theirs would still be the first, and stop there rather than at the
# end, so the error is the one the whole batch would meet. Any other
# error lowers it to zero, and a chunk stopped by it throws
# `ChunkAborted`, which the chunks' error leaves out.
struct ChunkAborted <: Exception end
const CHUNKFLAG = :transientchunkflag

# the stop of a chunk at a step past the flag its task holds, if it holds
# one; a step counts forward in a solve and a tangent and backward in an
# adjoint, whose errors all lower the flag to zero
function checkchunks(step::Integer)
    flag = get(task_local_storage(), CHUNKFLAG, nothing)
    flag isa Threads.Atomic{Int} && step > flag[] && throw(ChunkAborted())
    return nothing
end

# `f(k, chunk)` for every chunk of conditions, returned as a vector: on
# the calling task for one chunk, and otherwise each chunk on a task of
# its own across the threads of the session, which stop after a failure
# as the flag above says, a failure thrown as the chunks' error (see
# `chunkerror`). `f` is called once a chunk, and compiled for its own
# type; this is compiled once whatever closure it is given.
function runchunks(@nospecialize(f), chunks)
    length(chunks) == 1 && return Any[f(1, only(chunks))]
    each = Vector{Any}(undef, length(chunks))
    flag = Threads.Atomic{Int}(typemax(Int))
    try
        Base.Threads.@sync for (k, ch) in enumerate(chunks)
            Base.Threads.@spawn begin
                task_local_storage(CHUNKFLAG, flag)
                try
                    each[k] = f(k, ch)
                catch err
                    err isa ChunkAborted || Threads.atomic_min!(flag, err isa TransientStepError ? err.step : 0)
                    rethrow()
                end
            end
        end
    catch err
        throw(chunkerror(err))
    end
    return each
end

# One chunk of a batch stepped along `times` into the arrays `out`, on
# the factorizations `bf`, from the initial arrays `init` of its
# conditions, whose indices in the batch are `conditions`: every step
# advanced and saved. The record is whichever arrays of `out` are there to
# be filled. Returns the chunk's statistics.
function gaussbatchrun!(out, sys::TransientSystem, problems, init, times, saveevery,
        rtol, atol, iterations, bf::GaussBatchFactor, conditions)
    p = sys.problem
    backend = sys.backend
    N = length(problems)
    n = length(p)
    t0, nsteps = first(times), length(times) - 1
    savephases, savestates = recorded(out.phases), recorded(out.flux)
    saveendphases, saveendrates, saveblockstates = recorded(out.endphases), recorded(out.endrates), recorded(out.blockstates)
    savelinewaves, savehistory = recorded(out.linewaves), recorded(out.history)
    checkpoints = out.checkpoints
    savecheckpoints = checkpoints.every > 0
    K = checkpoints.every
    tails = savecheckpoints && size(checkpoints.waves, 2) > 0
    st = gaussstepper(sys, problems, rtol, atol, iterations, bf; conditions)
    nl = 2length(p.lines)
    npre = st.npre
    # the rate along the algebraic directions, read from the
    # differentiated constraints wherever the state is reported, the
    # lines' forced currents' rate by the same difference the drives' is
    rw = ratereadwork(sys, backend, n, N, N)
    delta = ratedelta(sys)
    vread = KernelAbstractions.zeros(backend, Float64, n, N)
    reading = savestates || sys.portsread || saveendrates
    readout!(dst, t) = readstepper!(dst, st, rw, t, t0, delta)
    # the outputs of a saved time: the port waves from the read rate where
    # a port reads one along an algebraic direction, the rate across the
    # projected junctions, and the states
    saveoutputs!(saved, t) = begin
        reading && readout!(vread, t)
        portwaves!(view(out.voltage, :, saved, :), view(out.incident, :, saved, :),
            view(out.outgoing, :, saved, :), sys, sys.portsread ? vread : st.v, st.values, st.portwork)
        if savephases
            view(out.phases, :, 1, saved, :) .= stage(st.phi, 1)
            view(out.phases, :, 2, saved, :) .= stage(st.phi, 2)
        end
        saveendphases && copyto!(view(out.endphases, :, saved, :), st.pw.phip)
        if saveendrates
            stepmul!(rw.pw.phip, sys.projection.RJp, vread)
            copyto!(view(out.endrates, :, saved, :), rw.pw.phip)
        end
        if savestates
            copyto!(view(out.flux, :, saved, :), st.x)
            copyto!(view(out.rate, :, saved, :), vread)
            view(out.stages, :, 1, saved, :) .= stage(st.delta, 1)
            view(out.stages, :, 2, saved, :) .= stage(st.delta, 2)
            saveblockstates && (view(out.blockstates, :, saved, :) .= st.rw.z)
        end
        nothing
    end
    isnothing(st.rw) || copyto!(st.rw.z, init.z0)
    copyto!(out.initialwaves, init.w0)
    copyto!(out.initialstates, init.z0)
    # With checkpoints, the factorizations begin anew at the initial
    # state, whatever a kept factor held, and every window of steps
    # starts on its checkpoint as the replay starts it (see
    # `responsewindows`), so the replay retraces the solve's steps exactly
    setstate!(st, init.x0, init.v0, nothing; fresh = savecheckpoints)
    # the history of the line waves before the start; with a record of
    # the phases or the states the waves leaving every port at every step,
    # which the responses read, and with checkpoints the history before
    # each of them, or those waves where that history would outweigh them
    tail = tobackend(backend, init.tailh)
    sethistory!(st, times, tail, 0)
    savehistory && copyto!(out.history, tail)
    savelinewaves && (view(out.linewaves, :, 1, :) .= tobackend(backend, init.w0))
    saveendphases && projectedphases!(st, st.x)
    stepperdrives!(st, t0)
    saveoutputs!(1, t0)
    copyto!(out.initialflux, st.x)
    copyto!(out.initialrate, st.v)
    saved = 1
    tprev = t0
    for step in 1:nsteps
        checkchunks(step)
        t = times[step + 1]
        if savecheckpoints && (step - 1) % K == 0
            c = (step - 1) ÷ K + 1
            copyto!(view(checkpoints.flux, :, c, :), st.x)
            copyto!(view(checkpoints.rate, :, c, :), st.v)
            isnothing(st.rw) || (view(checkpoints.states, :, c, :) .= st.rw.z)
            view(checkpoints.increment, :, 1, c, :) .= stage(st.lastdelta, 1)
            view(checkpoints.increment, :, 2, c, :) .= stage(st.lastdelta, 2)
            nl > 0 && tails && readhistory!(view(checkpoints.waves, :, :, c, :), st.waves, step)
            setstate!(st, st.x, st.v, st.lastdelta)
        end
        advance!(st, tprev, t, step)
        tprev = t
        savelinewaves && (view(out.linewaves, :, step + 1, :) .=
            view(st.waves, :, ringslot(npre + step, size(st.waves, 2)), :))
        if step % saveevery == 0 || step == nsteps
            saved += 1
            stepperdrives!(st, t)
            saveoutputs!(saved, t)
        end
    end
    KernelAbstractions.synchronize(backend)
    copyto!(out.finalflux, st.x)
    readout!(out.finalrate, tprev)
    nl > 0 && readhistory!(out.finalwaves, st.waves, nsteps + 1)
    isnothing(st.rw) || copyto!(out.finalstates, st.rw.z)
    nw = st.newtonwork
    return (; steps = nsteps, newtoncorrections = st.corrections, factorizations = st.factorizations + maximum(nw.factorizations),
        retries = maximum(nw.retries), kryloviterations = 0, rtol = Float64(rtol), atol = Float64(atol), iterations = Int(iterations))
end

# A window of steps a response walks: its steps, `phases!(buffer, k)`
# filling the `(nj, N, 2)` buffer with the phases of the step from time
# `k` to `k + 1`, `endphases!(buffer, k)` and `endrates!(buffer, k)`
# filling the `(npj, N)` buffer with the projected junctions' phases and
# read rates at time `k`, the replay to run before the window, if any,
# and `state!(x, v, k)` and `stages!(buffer, k)` filling the states, with
# their rates read as the solve read them, and the stage increments where
# a response reads them. The readers are held
# untyped, so that a window of a phase record and a window of
# checkpoints, with or without the states, are one type and a response
# compiles once for every kind; a reader is one dynamic call per step.
struct ResponseWindow
    steps::UnitRange{Int}
    phases!::Any
    endphases!::Any
    endrates!::Any
    replay::Any
    state!::Any
    stages!::Any
end

# The windows of steps a response walks and the phases of each. With a
# phase record, one window of every step reading the record; with
# checkpoints, one window per checkpoint, replayed from its state by the
# stepper into a buffer of the window's stage phases, walked forward by
# the tangent and backward by the adjoint, on the stepper given, a
# workspace's kept one, or one built here. The return type is declared,
# as the call is dynamic on the solution's untyped record and the
# responses step on the windows.
function responsewindows(sol::TransientBatchSolution, sys::TransientSystem, forward::Bool;
        withstates::Bool = false, stepper = nothing)::Vector{ResponseWindow}
    nt = length(sol.times)
    pr = sys.projection
    npj = isnothing(pr) ? 0 : length(pr.pj)
    if !isnothing(sol.phases)
        phases, endphases = sol.phases, sol.endphases
        phases! = (buffer, k) -> begin
            view(buffer, :, :, 1) .= view(phases, :, 1, k + 1, :)
            view(buffer, :, :, 2) .= view(phases, :, 2, k + 1, :)
            nothing
        end
        endphases! = (buffer, k) -> (npj > 0 && copyto!(buffer, view(endphases, :, k, :)); nothing)
        endrates = sol.endrates
        endrates! = (buffer, k) -> (npj > 0 && copyto!(buffer, view(endrates, :, k, :)); nothing)
        # the states and the stage increments of the record, when asked for
        flux, rate, stages = sol.flux, sol.rate, sol.stages
        withstates && isnothing(stages) && throw(ArgumentError("the record holds no states; solve with record = :states or :checkpoints."))
        state! = !withstates ? nothing : (x, v, k) -> (copyto!(x, view(flux, :, k, :)); copyto!(v, view(rate, :, k, :)); nothing)
        stages! = !withstates ? nothing : (buffer, k) -> begin
            view(buffer, :, :, 1) .= view(stages, :, 1, k + 1, :)
            view(buffer, :, :, 2) .= view(stages, :, 2, k + 1, :)
            nothing
        end
        return [ResponseWindow(1:nt - 1, phases!, endphases!, endrates!, nothing, state!, stages!)]
    end
    cps = sol.checkpoints
    K, nc = cps.every, size(cps.flux, 2)
    N = length(sol.problems)
    nj = length(sys.lmolj)
    st = isnothing(stepper) ? replaystepper(sys, sol) : stepper
    # the factorizations begun anew at the first checkpoint, the initial
    # state, as the solve began them, so that each window's are the ones
    # the solve had at its checkpoint whatever window was replayed before
    setstate!(st, view(cps.flux, :, 1, :), view(cps.rate, :, 1, :), nothing; fresh = true)
    nl = 2length(sys.problem.lines)
    npre = st.npre
    record = sol.linewaves
    buffer = KernelAbstractions.zeros(sys.backend, Float64, nj, 2, K + 1, N)
    ebuffer = KernelAbstractions.zeros(sys.backend, Float64, npj, K + 2, N)
    # the rate across the projected junctions at every time of the
    # window, read as the solve read it
    rbuffer = KernelAbstractions.zeros(sys.backend, Float64, npj, K + 2, N)
    # the states of the window, their rates read as the solve read them,
    # and the stage increments of its steps, for the responses which read
    # them
    n = length(sys.problem)
    ns = withstates ? n : 0
    xbuffer = KernelAbstractions.zeros(sys.backend, Float64, ns, K + 2, N)
    vbuffer = KernelAbstractions.zeros(sys.backend, Float64, ns, K + 2, N)
    dbuffer = KernelAbstractions.zeros(sys.backend, Float64, ns, 2, K + 1, N)
    # the work and the buffers the replay captures are there whether or
    # not it reads them, empty then, so that its closure is one type
    rw = ratereadwork(sys, sys.backend, length(sys.problem), N, N)
    vread = KernelAbstractions.zeros(sys.backend, Float64, length(sys.problem), N)
    delta, t0 = ratedelta(sys), first(sol.times)
    reads = npj > 0 || withstates
    readrates! = (column, t) -> begin
        readstepper!(vread, st, rw, t, t0, delta)
        if npj > 0
            stepmul!(rw.pw.phip, pr.RJp, vread)
            copyto!(view(rbuffer, :, column, :), rw.pw.phip)
        end
        withstates && copyto!(view(vbuffer, :, column, :), vread)
        nothing
    end
    windows = map(1:nc) do c
        kstart = (c - 1)*K + 1
        kend = min(c*K, nt - 1)
        replay = () -> begin
            lastdelta = similar(st.lastdelta)
            stage(lastdelta, 1) .= view(cps.increment, :, 1, c, :)
            stage(lastdelta, 2) .= view(cps.increment, :, 2, c, :)
            setstate!(st, view(cps.flux, :, c, :), view(cps.rate, :, c, :), lastdelta)
            isnothing(st.rw) || (st.rw.z .= view(cps.states, :, c, :))
            if nl > 0 && size(cps.waves, 2) > 0
                sethistory!(st, sol.times, view(cps.waves, :, :, c, :), kstart - 1)
            elseif nl > 0
                # the history columns `kstart` to `npre + kstart - 1` from the
                # record, whose column `j` is history column `npre + j - 1`,
                # and before it the history the solve started from, whose
                # column `npre` the record's first repeats
                st.times = sol.times
                before = max(npre - kstart, 0)
                if before > 0
                    fillhistory!(st.waves, view(sol.history, :, kstart:npre - 1, :), kstart)
                end
                fillhistory!(st.waves, view(record, :, max(kstart - npre + 1, 1):kstart, :), kstart + before)
                st.accepted = npre + kstart - 1
            end
            npj > 0 && copyto!(view(ebuffer, :, 1, :), projectedphases!(st, st.x))
            reads && readrates!(1, sol.times[kstart])
            withstates && copyto!(view(xbuffer, :, 1, :), st.x)
            for k in kstart:kend
                advance!(st, sol.times[k], sol.times[k + 1], k)
                view(buffer, :, 1, k - kstart + 2, :) .= stage(st.phi, 1)
                view(buffer, :, 2, k - kstart + 2, :) .= stage(st.phi, 2)
                npj > 0 && copyto!(view(ebuffer, :, k - kstart + 2, :), st.pw.phip)
                reads && readrates!(k - kstart + 2, sol.times[k + 1])
                if withstates
                    view(dbuffer, :, 1, k - kstart + 1, :) .= stage(st.delta, 1)
                    view(dbuffer, :, 2, k - kstart + 1, :) .= stage(st.delta, 2)
                    copyto!(view(xbuffer, :, k - kstart + 2, :), st.x)
                end
            end
            # the replayed window ends where the next checkpoint was taken,
            # to roundoff, as it retraces the solve's steps
            if c < nc
                next = view(cps.flux, :, c + 1, :)
                gap = maximum(abs, st.x .- next)
                gap <= 8eps(Float64)*maximum(abs, next) || error(
                    lazy"the replay of a window of checkpoints does not reach the next checkpoint, by $(gap); the record was not reproduced.")
            end
            nothing
        end
        phases! = (b, k) -> begin
            view(b, :, :, 1) .= view(buffer, :, 1, k - kstart + 2, :)
            view(b, :, :, 2) .= view(buffer, :, 2, k - kstart + 2, :)
            nothing
        end
        endphases! = (b, k) -> (npj > 0 && copyto!(b, view(ebuffer, :, k - kstart + 1, :)); nothing)
        endrates! = (b, k) -> (npj > 0 && copyto!(b, view(rbuffer, :, k - kstart + 1, :)); nothing)
        state! = !withstates ? nothing : (x, v, k) -> begin
            copyto!(x, view(xbuffer, :, k - kstart + 1, :))
            copyto!(v, view(vbuffer, :, k - kstart + 1, :))
            nothing
        end
        stages! = !withstates ? nothing : (b, k) -> begin
            view(b, :, :, 1) .= view(dbuffer, :, 1, k - kstart + 1, :)
            view(b, :, :, 2) .= view(dbuffer, :, 2, k - kstart + 1, :)
            nothing
        end
        ResponseWindow(kstart:kend, phases!, endphases!, endrates!, replay, state!, stages!)
    end
    return forward ? windows : reverse(windows)
end

# The stepper a record of checkpoints is replayed on: the batch's
# stepper on its own factorizations, at the solve's tolerances and under
# the record's problems. A workspace keeps its own across the responses
# of its shape, the factorizations with their analyses along with it,
# and hands it over set to the record in hand.
function replaystepper(sys::TransientSystem, sol::TransientBatchSolution)
    return gaussstepper(sys, sol.problems, sol.stats.rtol, sol.stats.atol, sol.stats.iterations,
        gaussbatchfactor(sys, length(sol.problems)))
end
function replaystepper!(w, sol::TransientBatchSolution, sys::TransientSystem)
    isnothing(sol.checkpoints) && return nothing
    st = w.replay
    if !(st isa GaussStepper) || st.sys !== sys || st.N != length(sol.problems)
        st = replaystepper(sys, sol)
        w.replay = st
    end
    st.problems = sol.problems
    st.rtol, st.atol, st.iterations = sol.stats.rtol, sol.stats.atol, sol.stats.iterations
    return st
end

# an array of a solution with a trailing dimension of one: a condition
# of a batch is a view of the batch's array at its index, and gets that
# index back as a range, a view of the same array again rather than a
# reshape of the view, so that the readers of a response walk a strided
# view as they do for a range of conditions
function expand(a)
    isnothing(a) && return nothing
    if a isa SubArray && last(parentindices(a)) isa Integer
        j = last(parentindices(a))
        return view(parent(a), Base.front(parentindices(a))..., j:j)
    end
    return reshape(a, size(a)..., 1)
end

# an ordinary solution as a batch of one, its arrays with a trailing
# dimension of one: the responses' entries take the batch form, and a
# solution through this
batchof(sol::TransientSolution) =
    TransientBatchSolution([sol.problem], sol.method, sol.dt, sol.times, solutionarrays(expand, sol), sol.stats)
# the arrays of a response of a batch of one without the trailing
# dimension of one, as the response of the solution, and the reverse
dropcondition(a) = isnothing(a) ? nothing : reshape(a, size(a)[1:end-1]...)
addcondition(a) = isnothing(a) ? nothing : reshape(a, size(a)..., 1)

# a batch of one is an ordinary solution
function unbatch(b::TransientBatchSolution)
    squeeze = a -> isnothing(a) ? nothing : reshape(a, size(a)[1:end-1]...)
    return TransientSolution(b.problems[1], b.method, b.dt, b.times, solutionarrays(squeeze, b), b.stats)
end
