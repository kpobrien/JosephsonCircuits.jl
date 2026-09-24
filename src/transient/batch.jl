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
function transientproblem(p::TransientProblem; sources)
    psc = p.circuit
    drives, injection, constantcurrent = bindsources(psc, p.matrices.vvn, p.ports, p.portpositive, p.portnegative, length(p), sources)
    return TransientProblem(psc, p.matrices, p.Nnodal, p.Naux, p.Lscale, p.coupledbranches,
        p.floatingcomponents, p.gaugeindices, p.inertialess, p.algebraic, p.directions, p.constraints,
        injection, drives, constantcurrent, p.ports, p.portpositive, p.portnegative, p.portimpedances, p.portconductances, p.blocks, p.lines,
        p.relations, p.C, p.G, p.L, p.lineE, p.RJ, p.lmolj)
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
# a range of conditions is a batch of views, so that the noise can tile
# the conditions of a batch within its memory
function Base.getindex(b::TransientBatchSolution, js::AbstractVector{<:Integer})
    all(j -> 1 <= j <= length(b), js) || throw(BoundsError(b, js))
    slice = a -> isnothing(a) ? nothing : selectdim(a, ndims(a), js)
    cp = mapcheckpoints(slice, b.checkpoints)
    return TransientBatchSolution(b.problems[js], b.method, b.dt, b.times, slice(b.voltage), slice(b.incident),
        slice(b.outgoing), slice(b.phases), slice(b.endphases), slice(b.endrates), slice(b.linewaves), slice(b.history), slice(b.flux), slice(b.rate), slice(b.stages), cp,
        slice(b.finalwaves), slice(b.finalstates), slice(b.initialflux),
        slice(b.initialrate), slice(b.finalflux), slice(b.finalrate), slice(b.blockstates), slice(b.initialwaves), slice(b.initialstates), b.stats)
end
function Base.getindex(b::TransientBatchSolution, j::Integer)
    1 <= j <= length(b) || throw(BoundsError(b, j))
    slice = a -> isnothing(a) ? nothing : selectdim(a, ndims(a), j)
    cp = mapcheckpoints(slice, b.checkpoints)
    return TransientSolution(b.problems[j], b.method, b.dt, b.times, slice(b.voltage), slice(b.incident),
        slice(b.outgoing), slice(b.phases), slice(b.endphases), slice(b.endrates), slice(b.linewaves), slice(b.history), slice(b.flux), slice(b.rate), slice(b.stages), cp,
        slice(b.finalwaves), slice(b.finalstates), copy(view(b.initialflux, :, j)),
        copy(view(b.initialrate, :, j)), copy(view(b.finalflux, :, j)), copy(view(b.finalrate, :, j)), slice(b.blockstates),
        slice(b.initialwaves), slice(b.initialstates), b.stats)
end

# The factorizations of a batch of Gauss-Legendre stage matrices, one per
# condition on one pattern: on the CPU a complex matrix and a KLU
# factorization per condition, on a device one complex value matrix with
# a column per condition and the uniform cuDSS batch over it, in chunks of
# at most the batch size cuDSS solves correctly. `solve!` applies every
# factorization to its column of a right hand side matrix.
# On a device the transpose of an unsymmetric operator is a second
# factorization on the transposed pattern, built from the same values
# through the permutation `tperm` when an adjoint first asks for it and
# refreshed with the forward factor thereafter; on the host KLU solves
# the transpose from the one factorization. A device chunk's sweep holds
# the solution and right hand side arrays of its `nrhs` columns per
# condition, allocated when it is built; the host has no chunks.
struct GaussBatchFactor{J, F, S, I, IT}
    ncolumns::Int
    nrhs::Int
    jacobians::J
    factors::F
    scratch::S
    chunks::Vector{UnitRange{Int}}
    rowptr::I
    colind::I
    symmetric::Bool
    tfactors::Vector{Any}
    trowptr::IT
    tcolind::IT
    tperm::Any
    tnzval::Any
    tstale::Base.RefValue{Bool}
end

function gaussbatchfactor(sys::TransientSystem, ncolumns::Int; nrhs::Int = 1)
    g = sys.gauss
    backend = sys.backend
    nnzj = nnz(sys.jacobian)
    if backend isa CPU
        pattern = g.cjacobian
        jacobians = [SparseMatrixCSC(size(pattern)..., SparseArrays.getcolptr(pattern), rowvals(pattern),
            zeros(ComplexF64, nnzj)) for _ in 1:ncolumns]
        # the assembly scratch belongs to the factor, not to the system,
        # so that several batches of one circuit may step at once
        return GaussBatchFactor(ncolumns, nrhs, jacobians, Vector{Any}(nothing, ncolumns),
            zeros(Float64, nnzj), UnitRange{Int}[],
            nothing, nothing, sys.symmetric, Any[], nothing, nothing, nothing, nothing, Ref(true))
    end
    n = size(sys.jacobian, 1)
    limit = uniformbatchlimit(nrhs)
    chunks = [first:min(first + limit - 1, ncolumns) for first in 1:limit:ncolumns]
    nzval = KernelAbstractions.zeros(backend, ComplexF64, nnzj, ncolumns)
    A = g.cjacobian
    hostrowptr, hostcolind = convert(Vector{Int}, rowpointer(A)), convert(Vector{Int}, columnindices(A))
    rowptr = tobackend(backend, convert(Vector{Int32}, hostrowptr))
    colind = tobackend(backend, convert(Vector{Int32}, hostcolind))
    trowptr, tcolind, tperm, tnzval = nothing, nothing, nothing, nothing
    if !sys.symmetric
        # the row structure of the transpose is the column structure of the
        # matrix, and the values of the transpose in its row order are the
        # matrix's values gathered by the permutation
        Kt = SparseMatrixCSC(n, n, hostrowptr, hostcolind, collect(1:nnzj))
        Kcsc = sparse(transpose(Kt))
        trowptr = tobackend(backend, convert(Vector{Int32}, SparseArrays.getcolptr(Kcsc)))
        tcolind = tobackend(backend, convert(Vector{Int32}, rowvals(Kcsc)))
        tperm = tobackend(backend, nonzeros(Kcsc))
        tnzval = KernelAbstractions.zeros(backend, ComplexF64, nnzj, ncolumns)
    end
    return GaussBatchFactor(ncolumns, nrhs, nzval, Vector{Any}(nothing, length(chunks)), nothing, chunks, rowptr, colind,
        sys.symmetric, Vector{Any}(nothing, length(chunks)), trowptr, tcolind, tperm, tnzval, Ref(true))
end

# the solution and right hand side arrays of a device sweep over the
# conditions `chunk`, which it keeps
sweeparrays(bf::GaussBatchFactor, sys::TransientSystem, chunk) =
    ntuple(_ -> KernelAbstractions.zeros(sys.backend, ComplexF64, size(sys.jacobian, 1), bf.nrhs, length(chunk)), 2)

# the transposed factors of a device batch from its current values, built
# or refreshed when an adjoint asks and the forward factor has been
# refreshed since
function gaussbatchtransposed!(bf::GaussBatchFactor, sys::TransientSystem)
    bf.tstale[] || return bf
    bf.tnzval .= view(bf.jacobians, bf.tperm, :)
    for (c, chunk) in enumerate(bf.chunks)
        if isnothing(bf.tfactors[c])
            bf.tfactors[c] = _cudss_sweep(bf.trowptr, bf.tcolind, bf.tnzval[:, chunk], sweeparrays(bf, sys, chunk)...;
                sys.factorization.kwargs...)
        else
            copyto!(bf.tfactors[c].nzval, view(bf.tnzval, :, chunk))
            _cudss_sweeprefactorize!(bf.tfactors[c])
        end
    end
    bf.tstale[] = false
    return bf
end

# The stage matrices at the stage phases of the columns of `columns`, a
# host mask, or of every column with `nothing`, assembled per column by
# the plan from the mean cosine of the two stages, and factorized, or
# refactorized on the analyses of the first time. On the host only the
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
            bf.factors[j] = (fresh || isnothing(bf.factors[j])) ? freshfactorization!(g.ordering, sys.factorization, A) :
                refactorize!(sys.factorization, bf.factors[j], A)
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
    bf.tstale[] = true
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

# the solve of every condition's block of right hand sides, or of the
# transposed operator's: the columns of `R` are `nrhs` per condition,
# condition after condition
function gaussbatchsolve!(Z::AbstractMatrix, bf::GaussBatchFactor, R::AbstractMatrix, transposed::Bool = false, sys = nothing)
    n = size(R, 1)
    nrhs = size(R, 2) ÷ bf.ncolumns
    transposed = transposed && !bf.symmetric
    if isempty(bf.chunks)
        for j in 1:bf.ncolumns
            solvecolumns!(Z, transposed ? transpose(bf.factors[j]) : bf.factors[j], R, (j - 1)*nrhs + 1, nrhs)
        end
        return Z
    end
    transposed && gaussbatchtransposed!(bf, sys)
    R3, Z3 = reshape(R, n, nrhs, :), reshape(Z, n, nrhs, :)
    for (c, chunk) in enumerate(bf.chunks)
        S = transposed ? bf.tfactors[c] : bf.factors[c]
        S.B .= view(R3, :, :, chunk)
        _cudss_sweepapply!(S, S.X, S.B)
        view(Z3, :, :, chunk) .= S.X
    end
    return Z
end

# The stage transform on a batch: the correction `[J*]^{-1} r` of a real
# stage pair `r` through the complex solve of every condition, `r̃ = Tinv r`,
# `z = S^{-1} r̃_1` and `d = T (z, conj z)`, which is real; with `transposed`
# the transpose `[J*]^{-T} r`, with the roles of `T` and `Tinv` exchanged
# and the transposed solve of `S`, which is `S` itself when symmetric.
function gaussbatchtransform!(d, r, bf::GaussBatchFactor, gc::GaussCoefficients, rc, zc, transposed::Bool = false, sys = nothing)
    if transposed
        rc .= gc.t11 .* stage(r, 1) .+ gc.t21 .* stage(r, 2)
    else
        rc .= gc.tinv11 .* stage(r, 1) .+ gc.tinv12 .* stage(r, 2)
    end
    gaussbatchsolve!(zc, bf, rc, transposed, sys)
    if transposed
        stage(d, 1) .= 2 .* real.(gc.tinv11 .* zc)
        stage(d, 2) .= 2 .* real.(gc.tinv12 .* zc)
    else
        stage(d, 1) .= 2 .* real.(gc.t11 .* zc)
        stage(d, 2) .= 2 .* real.(gc.t21 .* zc)
    end
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

# the work of the rational blocks, declared ahead of the stage solve and
# the residual, which take it or nothing unspecialized
abstract type AbstractRationalWork end

# The host vectors of a stage solve over `m` columns of `N` conditions,
# the mask of the columns still iterating on the backend as an `(m, 1)`
# array, and per condition whether its next correction's contraction is
# measured and whether its operator is refreshed, which the tangent and
# the adjoint keep so that their steps allocate none of them.
struct StageSolveWork{M}
    scale::Vector{Float64}
    previous::Vector{Float64}
    current::Vector{Float64}
    floor::Vector{Float64}
    blockfloor::Vector{Float64}
    active::Vector{Bool}
    between::Vector{Bool}
    activehost::Vector{Float64}
    activemask::M
    measuring::Vector{Bool}
    refresh::Vector{Bool}
end
stagesolvework(backend, m::Int, N::Int) = StageSolveWork([zeros(m) for _ in 1:5]..., fill(false, m), fill(false, m), zeros(m),
    KernelAbstractions.zeros(backend, Float64, m, 1), fill(false, N), fill(false, N))

# The exact solve of the stage equations of every condition linearized at
# its recorded stage phases, `J d = r` with the true stage Jacobian whose
# two stiffness blocks differ, or its transpose: a fixed point iteration
# on the frozen complex operators, `d += [J*]^{-1} (r - J d)`, contracting
# as the simplified Newton of the step does. `d` and `r` are `(n, ndir*N, 2)`
# with the directions of a condition contiguous; `phi` is `(nj, N, 2)`.
# A column iterates until it has converged and no further, so that its
# solution does not depend on the columns beside it, and a response is
# the same for a condition whatever batch or chunk it is stepped in.
# The contraction of each condition's first correction, the worst over
# its directions, measures how well its frozen operator fits this step's
# two stiffnesses. Where it leaves more than `stale` of the residual on an
# operator not assembled at this step's phases, those not in `fresh`, a
# vector per condition the solve adds to, the operator is refreshed at
# them, with `cosphi` the work of their mean stiffness and a pumped
# block's correction rebuilt on it, and the iteration goes on from where
# it is, as the step's Newton refreshes when it drifts, so a response
# solves wherever the forward solve did. Returns in `contraction` that of
# the first correction on the operator each condition ends on, for the
# caller to decide whether to refresh it at the next step.
function gaussbatchstagesolve!(contraction, d, r, sys::TransientSystem, gc::GaussCoefficients, bf::GaussBatchFactor, phi,
        transposed::Bool, res, dd, cwork, gwork, lwork, jwork, dwork, rc, zc,
        @nospecialize(rw::Union{Nothing, AbstractRationalWork}), sw::StageSolveWork, cosphi, fresh::Vector{Bool};
        rtol, iterations, stale = 0.25)
    h = sys.h
    fill!(d, 0)
    # the convergence per column: each direction of each condition against
    # its own right hand side, so a weak direction is not left at the
    # tolerance of a strong one, down to the roundoff floor of the terms
    # the residual sums
    scale, previous, current, floor, blockfloor = sw.scale, sw.previous, sw.current, sw.floor, sw.blockfloor
    active, between, measuring, refresh = sw.active, sw.between, sw.measuring, sw.refresh
    columnmax!(scale, r)
    scale .= max.(scale, floatmin(Float64))
    cmax, gmax = maximum(abs, gc.ainv2)/h^2, maximum(abs, gc.ainv)/h
    # The rounding of the terms the residual sums, which the iteration
    # stops on, bounded in two stages (see `gaussresidualroundoff!`): per
    # unit of the increment, the row or column sums of each matrix,
    # whichever is larger, since the adjoint multiplies by the
    # transposes, at the rule's weights, and the junction stamp at the
    # steepest slope of the relations at the stages, an upper bound of
    # any row's rounding; and for a column whose residual lies within that
    # bound and above its tolerance, the bound of every row from its own
    # terms. A product's own size would understate it along an algebraic
    # direction, where `C d` cancels. The eight units of roundoff allow
    # for the order of the summation; they are not a bound for rows of
    # every length.
    cs, gs, ls, js, _ = sys.rowsums
    rounding = 2cmax*cs + 2gmax*gs + ls + js*relationslope(sys, phi, dwork)
    copyto!(previous, scale)
    fill!(contraction, 0.0)
    fill!(blockfloor, 0.0)
    fill!(measuring, true)
    ndir = length(scale) ÷ size(phi, 2)
    nonlinear = length(sys.lmolj) > 0
    for iteration in 1:iterations
        for i in 1:2
            di = stage(d, i)
            stepmul!(stage(cwork, i), sys.C, di)
            stepmul!(stage(gwork, i), transposed ? sys.Gt : sys.G, di)
            stepmul!(stage(lwork, i), transposed ? sys.Lt : sys.L, di)
            batchjunctionproduct!(stage(res, i), sys, view(phi, :, :, i), di, jwork, dwork)
        end
        # the rational blocks' coupling of the stages, or its transpose,
        # and the rounding of the sum that carries it: as many units of
        # roundoff as the states it sums, which the matrices' eight do
        # not cover, `|S| |y|` row by row, which keeps every scatter entry
        # with the state it multiplies where a norm of each would pair the
        # largest row of one with the largest entry of the other, and the
        # poles of a fit span decades, and bounded by a row sum of `|S|`
        # times the largest state for the cheap bound. The bounds are per
        # column, so a strong direction does not lift a weak one's
        # threshold and the tangent and the adjoint stay transposes.
        if !isnothing(rw)
            transposed ? rationalsourcestranspose!(rw, sys, d) : rationalsources!(rw, sys, d, d; withstate = false)
            res .-= rw.source
            transposed || stackedscatterbound!(blockfloor, rw, sys.gauss.coupling)
        end
        columnmax!(floor, d)
        floor .= 8eps(Float64) .* (scale .+ rounding .* floor) .+ blockfloor
        for i in 1:2
            ri = stage(res, i)
            ri .+= stage(lwork, i)
            for l in 1:2
                a = transposed ? entry(gc.ainv, l, i) : entry(gc.ainv, i, l)
                b = transposed ? entry(gc.ainv2, l, i) : entry(gc.ainv2, i, l)
                ri .+= (b/h^2) .* stage(cwork, l) .+ (a/h) .* stage(gwork, l)
            end
        end
        res .= r .- res
        columnmax!(current, res)
        active .= current .> max.(rtol .* scale, floor)
        # the columns within the cheap bound and above their tolerance,
        # read row by row: at their floor when every row's residual lies
        # within the tolerance or that row's own rounding
        between .= .!active .& (current .> rtol .* scale)
        if any(between)
            excess = stageresidualexcess!(cwork, sys, cmax, gmax, res, r, d, phi, transposed, gwork, jwork, dwork, rw)
            active .|= between .& (excess .> rtol .* scale)
        end
        any(active) || return contraction
        if iteration >= 2 && any(measuring)
            # the contraction of the first correction on each condition's
            # operator, and the refresh of those that drifted
            for j in eachindex(measuring)
                measuring[j] || continue
                cols = (j - 1)*ndir + 1:j*ndir
                contraction[j] = maximum(k -> current[k]/previous[k], cols)
                refresh[j] = nonlinear && !fresh[j] && contraction[j] > stale && any(view(active, cols))
            end
            measuring .= refresh
            if any(refresh)
                gaussbatchjacobian!(bf, sys, phi, cosphi, dwork, refresh)
                isnothing(rw) || stagefactors!(rw, sys, bf, gc, rc, zc, transposed)
                fresh .|= refresh
            end
        end
        copyto!(previous, current)
        gaussbatchtransform!(dd, res, bf, gc, rc, zc, transposed, sys)
        isnothing(rw) || stagecorrect!(dd, rw, sys)
        # the columns which have converged are left as they are
        if !all(active)
            sw.activehost .= active
            copyto!(sw.activemask, sw.activehost)
            dd .*= reshape(sw.activemask, 1, :, 1)
        end
        d .+= dd
    end
    error("the stage equations of a Gauss-Legendre step did not converge in the tangent or the adjoint within stageiterations on an operator assembled at the step's own phases: the two stage stiffnesses differ too much from their mean; reduce dt or raise stageiterations of GaussLegendre.")
end

# The rounding of the products a stage residual sums, row by row, per
# unit of roundoff, into `bound`, `(n, m, 2)` stage by stage: the
# capacitance and the conductance on both stages' increments `d` at the
# largest of the rule's weights, `a2` and `a1`, and the stiffness on the
# stage's own, their transposes for an adjoint; `dsum` and `prod` are
# `(n, m)` work.
function stagerounding!(bound, sys::TransientSystem, a2, a1, d, transposed::Bool, dsum, prod)
    dsum .= abs.(stage(d, 1)) .+ abs.(stage(d, 2))
    stepmul!(prod, sys.Cabs, dsum)
    stage(bound, 1) .= a2 .* prod
    stepmul!(prod, transposed ? sys.Gtabs : sys.Gabs, dsum)
    stage(bound, 1) .+= a1 .* prod
    stage(bound, 2) .= stage(bound, 1)
    for i in 1:2
        dsum .= abs.(stage(d, i))
        stepmul!(prod, transposed ? sys.Ltabs : sys.Labs, dsum)
        stage(bound, i) .+= prod
    end
    return bound
end

# The rounding of the junction phases `RJ v` of each stage, carried to
# the rows by the junctions' stiffness at the stage phases `phi`,
# `(nj, N, 2)`: `|RJ'| (lmolj |relation'(phi)| (|RJ| |v|))`, added to
# `bound`, with `v` the stage values of the flux in a step or the
# increments of a stage solve, the columns of a condition contiguous;
# `vabs` and `prod` are `(n, m)` work, `jb` `(nj, m)` and `jd` the work of
# a polynomial relation's derivative.
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

# the largest residual of each column of `res` in a row where it exceeds
# that row's rounding bound, zero where every row is within its bound,
# on the host, or with `tol` the largest such residual over its row's
# tolerance
residualexcess(res, bound) = vec(Array(maximum(ifelse.(abs.(res) .> bound, abs.(res), 0.0); dims = (1, 3))))
residualexcess(res, bound, tol) = vec(Array(maximum(ifelse.(abs.(res) .> bound, abs.(res) ./ tol, 0.0); dims = (1, 3))))

# The excess of the stage solve's residual `res` over its rounding, row
# by row (see `residualexcess`), each row's bound from its own terms: the
# right hand side `r`, the matrices of `stagerounding!`, the junctions'
# stiffness on the rounding of the increments' phases and the blocks'
# `|S| |y|`. `bound` and `work` are `(n, m, 2)` work.
function stageresidualexcess!(bound, sys::TransientSystem, a2, a1, res, r, d, phi, transposed::Bool, work, jwork, dwork,
        @nospecialize(rw::Union{Nothing, AbstractRationalWork}))
    dsum, prod = stage(work, 1), stage(work, 2)
    stagerounding!(bound, sys, a2, a1, d, transposed, dsum, prod)
    junctionrounding!(bound, sys, phi, d, dsum, jwork, dwork, prod)
    bound .= 8eps(Float64) .* (bound .+ abs.(r))
    if !isnothing(rw) && !transposed
        reflectedwavesbound!(rw, sys.gauss.coupling)
        bound .+= sys.gauss.coupling.nstates*eps(Float64) .* rw.sbound
    end
    return residualexcess(res, bound)
end

# The largest magnitude of each column of `a`, over its rows and its
# stages or other trailing dimension, into the host vector `out`: a loop
# on the host, where a reduction over dimensions would allocate its
# result at every call, and a reduction on a device.
function columnmax!(out::Vector{Float64}, a::AbstractArray{<:Any,3})
    if a isa Array
        @inbounds for col in axes(a, 2)
            m = 0.0
            for i in axes(a, 3), row in axes(a, 1)
                m = max(m, abs(a[row, col, i]))
            end
            out[col] = m
        end
        return out
    end
    copyto!(out, vec(maximum(abs, a; dims = (1, 3), init = 0.0)))
    return out
end

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
# rather than for every rank a caller may give.
struct TangentCurrents
    grid::Array{Float64, 3}
    stages::Array{Float64, 4}
    staged::Bool
    constant::Bool
    # the caller gave a matrix, one direction, and the outputs carry no
    # direction dimension
    single::Bool
end
directions(c::TangentCurrents) = size(c.grid, 3)
# the column of the grid values read at recorded time `k`
gridcolumn(c::TangentCurrents, k::Int) = c.constant ? 1 : k
function tangentcurrents(currents, nq::Int, nt::Int, ndir::Int)
    isnothing(currents) && return TangentCurrents(zeros(nq, 1, ndir), zeros(nq, 2, 0, ndir), false, true, false)
    staged = ndims(currents) == 4
    grid = reshape(Float64.(collect(staged ? selectdim(currents, 2, 1) : currents)), nq, nt, ndir)
    stages = staged ? reshape(Float64.(collect(selectdim(currents, 2, 2:3))), nq, 2, nt, ndir) : zeros(nq, 2, 0, ndir)
    return TangentCurrents(grid, stages, staged, false, ndims(currents) == 2)
end

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
transposed operator, which the adjoint solves, has `U` and `V`
exchanged: `J*^(-T) W_ij'` and `U' J*^(-T) W_ij'` are built at a
refresh, and `K = J*^(-T) V` and `M = I + U' K` combined at each step.
The frozen operator never changes with the modulation, so a pumped
block refactorizes only when the junctions ask.
"""
mutable struct StageCorrection{S, K, R, Y}
    # the scatter onto the modulated rows, `U`, as the columns of an
    # `(n, np)` array, its transpose, and the rows of every `W_ij` on the
    # backend as the columns of a `(2n, np nw)` array, for the right hand
    # sides of the transposed solves
    Ucols::Y
    St::S
    Wcols::Y
    # the solves at the last refresh: `K = J*^(-1) U`, `(n, m, 2, r)`, or
    # for the transposed operator `J*^(-T) W_ij'` for every term, stage
    # and port, `(n, m, 2, np nw)`, of which `K` is combined at each step
    K::K
    X::K
    # the products with the solves per column, `W_ij K` as `(np, r, m, nw)`
    # or `U' J*^(-T) W_ij'` as `(r, np, m, nw)`, on the host, and the
    # small matrices `M` factorized per column
    Q::Array{Float64, 4}
    M::Vector{LU{Float64,Matrix{Float64},Vector{Int}}}
    transposed::Bool
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
    sstack::M
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
    return RationalWork(allocate(nz, N), allocate(2n, N), allocate(2n, N), allocate(2n, N), allocate(2n, N),
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
mutable struct GaussTangentWork{S, J, M, A, A4, C, B, LT, CO, O}
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
    # the currents on the backend, at the grid and at the stages, and the
    # grid values on the host
    dI::A
    dS::A4
    dIh::Array{Float64, 3}
    # the tangent's state, the stage unknowns and the stage work
    dx::M
    dv::M
    dxnew::M
    lx::M
    cdv::M
    work::M
    r::A
    d::A
    res::A
    dd::A
    cwork::A
    gwork::A
    lwork::A
    # the lines: the ring of the perturbation of the wave leaving each
    # port, the currents the arriving waves force and their rates, the
    # stencils and the tables, and the start of the prehistory
    dwaves::A
    ddlinevalues::M
    dlinerates::M
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
    cosphi::C
    rc::C
    zc::C
    # the currents at a stage and their injection
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
    # the initial state of every condition, on the backend
    xs0::M
    vs0::M
    # the stage factorizations, whether each condition's stage operator
    # is stale, the contraction of its last stage solve, and the stage
    # solve's host vectors
    bf::B
    stale::Vector{Bool}
    contraction::Vector{Float64}
    sw::StageSolveWork{M}
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
        tp::Vector{Int}, staged::Bool, constant::Bool, storing::Bool, perturbation)
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
    dI = allocate(nq, constant ? 1 : nt, ndir)
    dS = allocate(nq, 2, staged ? nt : 0, ndir)
    dx, dv, dxnew, lx, cdv, work = [allocate(n, m) for _ in 1:6]
    r, d, res, dd, cwork, gwork, lwork = [allocate(n, m, 2) for _ in 1:7]
    nl2 = 2length(p.lines)
    npre = lineprehistory(p, h)
    nring = linering(npre)
    dwaves = allocate(nl2, nring, m)
    ddlinevalues, dlinerates, linework = allocate(nl2, m), allocate(nl2, m), allocate(n, m)
    stencil, dstencil = zeros(8, nl2), allocate(8, nl2)
    rwt = isnothing(sys.gauss.coupling) ? nothing : rationalwork(p, backend, n, m)
    zerostages = isnothing(rwt) ? nothing : allocate(n, m, 2)
    fullstages = isnothing(rwt) ? nothing : allocate(n, m, 2)
    nj = length(sys.lmolj)
    jwork, phi = allocate(nj, m), allocate(nj, N, 2)
    # the derivative of the relation at one stage, for a circuit which has
    # a polynomial one; the Josephson relation is broadcast in place
    dwork = allsinusoidal(sys.relations) ? allocate(nj, 0) : allocate(nj, N)
    cosphi = KernelAbstractions.zeros(backend, ComplexF64, nj, N)
    rc, zc = [KernelAbstractions.zeros(backend, ComplexF64, n, m) for _ in 1:2]
    stagecurrent = allocate(nq, ndir)
    stageinjection, injectionall = allocate(n, ndir), allocate(n, m)
    # the outputs are the reset's, views of the result of a call
    voltage, incident, outgoing = [outputview(allocate(np, 0, m), 1:m) for _ in 1:3]
    coefficients = (outputcoefficients(sys, :voltage), outputcoefficients(sys, :incident), outputcoefficients(sys, :outgoing))
    portwork = allocate(np, m)
    portmap = devicesparse(sparse([q for q in tp if q > 0], [k for (k, q) in enumerate(tp) if q > 0],
        ones(count(>(0), tp)), np, nq), backend)
    directwork, directall = allocate(np, ndir), allocate(np, m)
    outwork = [allocate(np, m) for _ in 1:3]
    reading = sys.portsread ? outputreading(sys, backend, n, N, m) : nothing
    xs0, vs0 = allocate(n, N), allocate(n, N)
    bf = gaussbatchfactor(sys, N; nrhs = ndir)
    pr = sys.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, m)
    Zci = isnothing(pr) ? nothing : pr.Zch*injh
    # the direction's current along the rows of the endpoint reading
    Qinji = isnothing(pr) ? nothing : endpointinjection(pr, injh)
    return GaussTangentWork(sys, N, ndir, nq, nt, staged, constant, storing, injh, tp, perturbation, withstates, injection,
        dI, dS, zeros(nq, 1, ndir), dx, dv, dxnew, lx, cdv, work, r, d, res, dd, cwork, gwork, lwork,
        dwaves, ddlinevalues, dlinerates, linework, stencil, dstencil, linetables(p, backend), 0.0,
        rwt, zerostages, fullstages, jwork, phi, dwork, cosphi, rc, zc, stagecurrent, stageinjection, injectionall,
        voltage, incident, outgoing, coefficients, portwork, portmap, directwork, directall, outwork, reading, xs0, vs0,
        bf, fill(true, N), zeros(N), stagesolvework(backend, m, N), pw, Zci, Qinji, pwork, pend, dterm, nothing)
end

# whether a kept tangent workspace has the shape asked for and the
# perturbation
function sameshape(w::GaussTangentWork, sys::TransientSystem, N::Int, ndir::Int, nq::Int, nt::Int, injh, tp, staged::Bool,
        constant::Bool, storing::Bool, perturbation)
    return w.sys === sys && w.N == N && w.ndir == ndir && w.nq == nq && w.nt == nt && w.staged == staged &&
        w.constant == constant && w.storing == storing && w.injh == injh && w.tp == tp &&
        sameperturbation(w.perturbation, perturbation)
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

# the workspace at the start of a tangent: the currents, the initial
# perturbation and the initial states of the conditions in, the state
# and the histories cleared, fresh output arrays where they are handed
# over, and the stage operator to be assembled at the first step
function tangentreset!(w::GaussTangentWork, sol::TransientBatchSolution, currents::TangentCurrents, initial::TangentInitial,
        outputs)
    sys = w.sys
    p = sys.problem
    backend = sys.backend
    n, N, ndir = length(p), w.N, w.ndir
    copyto!(w.dI, currents.grid)
    copyto!(w.dS, currents.stages)
    w.dIh = currents.grid
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
    copyto!(w.xs0, hostmatrix(sol.initialflux, n, N))
    copyto!(w.vs0, hostmatrix(sol.initialrate, n, N))
    w.tpre = linestart(first(sol.times), lineprehistory(p, sys.h), sys.h)
    w.voltage, w.incident, w.outgoing = outputs
    fill!(w.stale, true)
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

# the projected junctions' phases and read rates at time `k`, and the
# states where the perturbation reads them, their rate read: at the start
# from the initial state, later from the window
function tangentreading!(w::GaussTangentWork, sol::TransientBatchSolution, k::Int, window)
    o, sys, pwork = w.reading, w.sys, w.pwork
    prj, pwn = sys.projection, o.rn.pw
    delta = ratedelta(sys)
    if k == 1
        isnothing(pwn) || (stepmul!(pwn.phip, prj.RJp, w.xs0); copyto!(o.hphi, pwn.phip))
        isnothing(prj) || drivedotz!(o.rn.bdotz, prj, sol.problems, sol.times[1], delta, o.rn.hv1, o.rn.hv2)
        readrate!(o.wk, w.vs0, w.xs0, sys, o.rn)
        isnothing(pwn) || (stepmul!(pwn.phip, prj.RJp, o.wk); copyto!(o.er, pwn.phip))
        w.withstates && (copyto!(pwork.x, w.xs0); copyto!(pwork.v, o.wk))
    else
        if !isnothing(pwn)
            window.endphases!(o.erdev, k)
            copyto!(o.hphi, o.erdev)
            window.endrates!(o.erdev, k)
            copyto!(o.er, o.erdev)
        end
        if w.withstates
            window.state!(pwork.x, pwork.v, k)
            isnothing(prj) || drivedotz!(o.rn.bdotz, prj, sol.problems, sol.times[k], delta, o.rn.hv1, o.rn.hv2)
            readrate!(o.wk, pwork.v, pwork.x, sys, o.rn)
            copyto!(pwork.v, o.wk)
        end
    end
    return nothing
end

# the port outputs of time `k`, stored, or handed to the sink as the three
# port by column matrices; where a port reads a rate along an algebraic
# direction, from the tangent rate read as the solve's is, the reading
# linearized at the projected junctions' recorded phases and read rates,
# with the tangent currents' rate and the components' perturbation of the
# constraints
function tangentoutputs!(w::GaussTangentWork, sol::TransientBatchSolution, k::Int, window, @nospecialize(outputsink))
    sys = w.sys
    np, ndir, N = size(w.portwork, 1), w.ndir, w.N
    reading, pwork = w.reading, w.pwork
    if !isnothing(reading)
        tangentreading!(w, sol, k, window)
        readoutputs!(reading, sys, k, w.dv, w.dx, w.dIh, w.injh, w.withstates ? Array(pwork.x) : zeros(0, 0),
            w.withstates ? Array(pwork.v) : zeros(0, 0), w.perturbation, sys.backend)
    end
    stepmul!(w.portwork, sys.ports, isnothing(reading) ? w.dv : reading.dvread)
    w.portwork .*= phi0
    stepmul!(w.directwork, w.portmap, view(w.dI, :, w.constant ? 1 : k, :))
    reshape(w.directall, np, ndir, N) .= reshape(w.directwork, np, ndir, 1)
    if w.storing
        tangentwrite!((view(w.voltage, :, k, :), view(w.incident, :, k, :), view(w.outgoing, :, k, :)), w, sol, k)
    else
        tangentwrite!((w.outwork[1], w.outwork[2], w.outwork[3]), w, sol, k)
        outputsink(k, w.outwork[1], w.outwork[2], w.outwork[3])
    end
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
        dIz = repeat(w.Zci*view(w.dIh, :, kcol, :), 1, N)
        nl2 > 0 && (dIz .+= pr.Zcline*Array(w.ddlinevalues))
        isnothing(rwt) || (dIz .+= pr.Zcblock*restingwaves(sys, rwt.z, sol.times[k + 1]))
        isnothing(pend) || (dIz .+= pr.Zch*pend.F)
        pw.g .-= dIz
        projectionsolve!(pw, pr, false)
        pw.alpha .*= -1
        copyto!(pw.dalpha, pw.alpha)
        stepmul!(pw.work, pr.Z, pw.dalpha)
        w.dxnew .+= pw.work
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
        dIq = repeat(w.Qinji*view(w.dIh, :, kcol, :), 1, N)
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
# blocks' terms, solved on the stage operator, refreshed at the recorded
# phases before the solve when the last step's first correction did not
# contract, and within it when this step's does not
function tangentstep!(w::GaussTangentWork, sol::TransientBatchSolution, k::Int, window)
    sys, rwt, pwork, bf = w.sys, w.rwt, w.pwork, w.bf
    gc = sys.gauss.coefficients
    h = sys.h
    N, ndir, nt = w.N, w.ndir, w.nt
    n = size(w.dx, 1)
    m = ndir*N
    nj = length(sys.lmolj)
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
            w.stagecurrent .= view(w.dS, :, i, k, :)
        else
            indices, weights = gaussstencil(k, nt, gc.c[i])
            fill!(w.stagecurrent, 0)
            for (idx, wgt) in zip(indices, weights)
                w.stagecurrent .+= wgt .* view(w.dI, :, w.constant ? 1 : idx, :)
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
    # the stale operators refreshed, and a pumped block's correction at
    # this step's weights, on the refreshed factorizations if any
    if any(w.stale)
        gaussbatchjacobian!(bf, sys, w.phi, w.cosphi, w.dwork, w.stale)
        isnothing(rwt) || stagefactors!(rwt, sys, bf, gc, w.rc, w.zc, false)
    else
        isnothing(rwt) || stagemultipliers!(rwt, sys)
    end
    # the stale operators, just refreshed, are the fresh ones of the solve
    contraction = gaussbatchstagesolve!(w.contraction, w.d, w.r, sys, gc, bf, w.phi, false, w.res, w.dd, w.cwork, w.gwork, w.lwork,
        w.jwork, w.dwork, w.rc, w.zc, rwt, w.sw, w.cosphi, w.stale; rtol = sys.gauss.stagertol, iterations = sys.gauss.stageiterations)
    w.stale .= (nj > 0) .& (contraction .> 0.25)
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
        @nospecialize(reuse); chunks = isnothing(outputsink) ? batchchunks(sys.backend, length(sol.problems)) : [1:length(sol.problems)])
    N, ndir, nq, nt = length(sol.problems), directions(currents), size(injh, 2), length(sol.times)
    storing = isnothing(outputsink)
    shape = (w, Nc) -> w isa GaussTangentWork && sameshape(w, sys, Nc, ndir, nq, nt, injh, tp, currents.staged, currents.constant, storing, perturbation)
    build = Nc -> gausstangentwork(sys, Nc, ndir, nq, nt, injh, tp, currents.staged, currents.constant, storing, perturbation)
    ws = chunkworkspaces(reuse, :tangent, chunks, shape, build)
    # the result of the call, fresh, so that the results of earlier calls
    # stay what they were
    np = length(sys.problem.portimpedances)
    voltage, incident, outgoing = [KernelAbstractions.zeros(sys.backend, Float64, np, storing ? nt : 0, ndir*N) for _ in 1:3]
    columns = ch -> (first(ch) - 1)*ndir + 1:last(ch)*ndir
    outputs = ch -> (outputview(voltage, columns(ch)), outputview(incident, columns(ch)), outputview(outgoing, columns(ch)))
    parts = runchunks(chunks) do c, ch
        length(chunks) == 1 ? Base.invokelatest(gaussbatchtangent, ws[1], sol, currents, initial, outputsink, outputs(1:N)) :
            Base.invokelatest(gaussbatchtangent, ws[c], sol[ch], currents, conditionslice(initial, ch), outputsink, outputs(ch))
    end
    ends = joinconditions(parts)
    shape = a -> begin
        b = reshape(a, np, size(a, 2), ndir, N)
        currents.single ? reshape(b, np, size(a, 2), N) : b
    end
    return (; voltage = storing ? shape(voltage) : nothing, incident = storing ? shape(incident) : nothing,
        outgoing = storing ? shape(outgoing) : nothing, ends.finalflux, ends.finalrate)
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
        initial::TangentInitial, @nospecialize(outputsink), outputs)
    sys = w.sys
    N, ndir = w.N, w.ndir
    try
        tangentreset!(w, sol, currents, initial, outputs)
        windows = responsewindows(sol, sol.problems, sys, true; withstates = w.withstates, stepper = replaystepper!(w, sol, sys))
        tangentoutputs!(w, sol, 1, first(windows), outputsink)
        for window in windows
            isnothing(window.replay) || window.replay()
            for k in window.steps
                checkchunks(k)
                tangentstep!(w, sol, k, window)
                tangentoutputs!(w, sol, k + 1, window, outputsink)
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
        (isnothing(sys.invariant) && isnothing(sys.projection)) ||
            readtangentrate!(finalrate, w.dv, w.dx, sys, sol, N, ndir, w.perturbation, w.dIh, w.injh, sys.backend)
        return (; finalflux = shape(copy(w.dx)), finalrate = shape(finalrate))
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
    w.dIh = zeros(size(w.dIh, 1), 0, size(w.dIh, 3))
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
mutable struct GaussAdjointWork{S, J, M, A, C, B, LT, V, O}
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
    # the cotangents of the state, the stage multipliers and the stage work
    xbar::M
    vbar::M
    work::M
    lmu::M
    cmu::M
    wst::A
    mu::A
    res::A
    dd::A
    cwork::A
    gwork::A
    lwork::A
    # the junctions at the stages
    jwork::M
    phi::A
    dwork::M
    cosphi::C
    rc::C
    zc::C
    # the objective's weights on the ports and their state
    portwork::M
    objectivestate::M
    transposing::Union{Nothing, OutputTranspose{M, SparseArrays.UMFPACK.UmfpackLU{Float64, Int}}}
    # the final state of every condition, on the backend
    xf::M
    wf::M
    # the stage factorizations, whether each condition's stage operator
    # is stale, the contraction of its last stage solve, and the stage
    # solve's host vectors
    bf::B
    stale::Vector{Bool}
    contraction::Vector{Float64}
    sw::StageSolveWork{M}
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
    injectiont = devicesparse(sparse(transpose(injh)), backend)
    portmapt = devicesparse(sparse([k for (k, q) in enumerate(tp) if q > 0], [q for q in tp if q > 0],
        ones(count(>(0), tp)), nq, np), backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    w = allocate(np, nt, nobj)
    directwork, targetwork, objectivework = allocate(np, nobj), allocate(nq, m), allocate(nq, nobj)
    # the currents are the reset's, a view of the result of a call
    currents = outputview(allocate(nq, 0, m), 1:m)
    ringslots, ringcolumn = allocate(nq, m, RINGSLOTS), allocate(nq, m)
    xbar, vbar, work, lmu, cmu = [allocate(n, m) for _ in 1:5]
    wst, mu, res, dd, cwork, gwork, lwork = [allocate(n, m, 2) for _ in 1:7]
    nj = length(sys.lmolj)
    jwork, phi = allocate(nj, m), allocate(nj, N, 2)
    dwork = allsinusoidal(sys.relations) ? allocate(nj, 0) : allocate(nj, N)
    cosphi = KernelAbstractions.zeros(backend, ComplexF64, nj, N)
    rc, zc = [KernelAbstractions.zeros(backend, ComplexF64, n, m) for _ in 1:2]
    portwork = allocate(np, nobj)
    objectivestate = allocate(n, nobj)
    transposing = sys.portsread ? outputtranspose(sys, backend, n, N, m) : nothing
    xf, wf = allocate(n, N), allocate(n, N)
    bf = gaussbatchfactor(sys, N; nrhs = nobj)
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
        xbar, vbar, work, lmu, cmu, wst, mu, res, dd, cwork, gwork, lwork, jwork, phi, dwork, cosphi, rc, zc,
        portwork, objectivestate, transposing, xf, wf, bf, fill(true, N), zeros(N), stagesolvework(backend, m, N),
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
# terminations, and the stage operator to be assembled at the first step
function adjointreset!(w::GaussAdjointWork, sol::TransientBatchSolution, wh::Array{Float64, 3}, @nospecialize(sink), currents)
    sys = w.sys
    p = sys.problem
    n, N, nq = length(p), w.N, w.nq
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
    weights, directwork, objectivework, portmapt, cd = w.w, w.directwork, w.objectivework, w.portmapt, w.cd
    nobj = w.nobj
    feedthrough! = (column, k) -> begin
        directwork .= cd .* view(weights, :, k, :)
        stepmul!(objectivework, portmapt, directwork)
        reshape(column, nq, nobj, N) .+= reshape(objectivework, nq, nobj, 1)
        nothing
    end
    store = w.storing ? (k, values) -> (copyto!(view(currents, :, k, :), values); nothing) : sink
    w.ring = CurrentRing(w.ring.slots, w.ring.column, store, feedthrough!)
    fill!(w.stale, true)
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
# states where the perturbation reads them, their rate read: at the end
# from the final state, before it from the window
function adjointreading!(w::GaussAdjointWork, sol::TransientBatchSolution, k::Int, window)
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
        if !isnothing(pwn)
            window.endphases!(o.erdev, k)
            copyto!(o.hphi, o.erdev)
            window.endrates!(o.erdev, k)
            copyto!(o.er, o.erdev)
        end
        if w.withstates
            window.state!(pwork.x, pwork.v, k)
            isnothing(prj) || drivedotz!(o.rn.bdotz, prj, sol.problems, sol.times[k], ratedelta(sys), o.rn.hv1, o.rn.hv2)
            readrate!(o.wk, pwork.v, pwork.x, sys, o.rn)
            copyto!(pwork.v, o.wk)
        end
    end
    return nothing
end

# The adjoint of an output at time `k`: the rate receives
# `phi0 P' (cv .* w)`. Where a port reads a rate along an algebraic
# direction, through the transpose of the reading: the rate before the
# reading, the flux through the curvature of the projected junctions'
# stiffness, the currents whose rate the reading read, at the points of
# the stencil, and the components whose perturbation of the constraints
# it read, on the read rate at the time.
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
    adjointreading!(w, sol, k, window)
    readoutputstranspose!(o, sys, k, w.nt, w.vbar, w.xbar, w.ring, w.injectiont, w.targetwork, w.perturbation, w.pcwork,
        w.sensd, w.withstates ? w.pwork.x : o.wk, sys.backend)
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
    # the Newton step transposed: the multiplier from the flux along the
    # directions, carried by the constraints' rows to the flux, the
    # current, the lines' forced currents and the resting waves
    stepmul!(pw.zwork, pr.Zt, w.xbar)
    copyto!(pw.g, pw.zwork)
    projectionsolve!(pw, pr, true)
    copyto!(pw.dalpha, pw.alpha)
    stepmul!(pw.work, pr.Zl, pw.dalpha)
    stepmul!(w.targetwork, w.injectiont, pw.work)
    ringadd!(w.ring, k + 1, w.targetwork, 1.0)
    isnothing(pend) || endpointcontract!(w.sensh, pend, w.perturbation, Array(pw.work), 1.0)
    if nl2 > 0
        adjointforced!(w, pw.work)
        adjointlinescatter!(w, sol.times[k + 1], npre + k - 1, 1.0)
    end
    isnothing(rwa) || restingwavesbar!(rwa, sys, pw.work, 1.0)
    stepmul!(w.lmu, sys.Lt, pw.work)
    w.xbar .-= w.lmu
    if !isempty(pr.pj)
        stepmul!(pw.phim, pr.RJp, pw.work)
        reshape(pw.phim, :, nobj, N) .*= pr.lmoljpdev .* reshape(cosp, :, 1, N)
        stepmul!(pw.work2, pr.RJpt, pw.phim)
        w.xbar .-= pw.work2
    end
    return nothing
end

# the adjoint's step from time `k + 1` back to `k` on the recorded stage
# phases of the window: the cotangent of the wave that left each port at
# the endpoint, the projection transposed, the stage multipliers on the
# transposed stage operator, refreshed at the recorded phases before the
# solve when the last step's first correction did not contract, and
# within it when this step's does not, and the multipliers' pull on the
# currents the stages read, the components, the lines and the state
function adjointstep!(w::GaussAdjointWork, sol::TransientBatchSolution, k::Int, window, @nospecialize(stagesink))
    sys, rwa, pwork, bf = w.sys, w.rwa, w.pwork, w.bf
    gc = sys.gauss.coefficients
    h = sys.h
    nt = w.nt
    nj = length(sys.lmolj)
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
    if any(w.stale)
        gaussbatchjacobian!(bf, sys, w.phi, w.cosphi, w.dwork, w.stale)
        isnothing(rwa) || stagefactors!(rwa, sys, bf, gc, w.rc, w.zc, true)
    else
        isnothing(rwa) || stagemultipliers!(rwa, sys)
    end
    # the stale operators, just refreshed, are the fresh ones of the solve
    contraction = gaussbatchstagesolve!(w.contraction, w.mu, w.wst, sys, gc, bf, w.phi, true, w.res, w.dd, w.cwork, w.gwork, w.lwork,
        w.jwork, w.dwork, w.rc, w.zc, rwa, w.sw, w.cosphi, w.stale; rtol = sys.gauss.stagertol, iterations = sys.gauss.stageiterations)
    w.stale .= (nj > 0) .& (contraction .> 0.25)
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
        windows = responsewindows(sol, sol.problems, sys, false; withstates = w.withstates, stepper = replaystepper!(w, sol, sys))
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
as the solve splits them; an `outputsink` receives the columns of every
condition at each time, so a call with one runs on one thread.
"""
function transienttangent(b::TransientBatchSolution, currents::Union{Nothing,AbstractArray{<:Real}};
        targets = porttargets(first(b.problems)), initialstate = nothing, factorization = nothing, reuse = nothing,
        outputsink = nothing, perturbation = nothing)
    # a batch of one under another rule is its solution's own tangent
    b.method isa GaussLegendre || return map(addcondition, transienttangent(singlecondition(b), currents;
        targets, initialstate, factorization, reuse, outputsink, perturbation))
    recordedsolution(b)
    p = first(b.problems)
    backend = KernelAbstractions.get_backend(b.finalflux)
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, b.dt, b.method, backend, fact)
    nq, nt, N = length(targets), length(b.times), length(b)
    ndir = tangentdirections(currents, perturbation, nq, nt)
    injection, ports = targetinjection(p, targets)
    initial = tangentinitial(initialstate, length(p), ndir, N, 2length(p.lines), lineprehistory(p, b.dt), blockstates(p))
    # invoked dynamically on the untyped kept system (see transientsolve)
    return Base.invokelatest(gaussbatchtangent, b, tangentcurrents(currents, nq, nt, ndir), injection, ports, initial, sys,
        outputsink, perturbation, reuse)
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
function transientadjoint(b::TransientBatchSolution, weights::AbstractArray{<:Real};
        quantity::Symbol = :outgoing, targets = porttargets(first(b.problems)), factorization = nothing,
        reuse = nothing, sink = nothing, stagesink = nothing, components = String[])
    # a batch of one under another rule is its solution's own adjoint
    b.method isa GaussLegendre || return map(addcondition, transientadjoint(singlecondition(b), weights;
        quantity, targets, factorization, reuse, sink, stagesink, components))
    recordedsolution(b)
    p = first(b.problems)
    backend = KernelAbstractions.get_backend(b.finalflux)
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, b.dt, b.method, backend, fact)
    wh = adjointweights(weights, length(p.portimpedances), length(b.times))
    perturbation = isempty(components) ? nothing : componentperturbation(p, components, backend; forcing = false)
    isnothing(perturbation) || recordedstates(b, perturbation)
    injection, ports = targetinjection(p, targets)
    # invoked dynamically on the untyped kept system (see transientsolve)
    return Base.invokelatest(gaussbatchadjoint, b, wh, ndims(weights) == 2, quantity, injection, ports, sys, sink,
        stagesink, perturbation, reuse)
end

# the stage residual of a batch on `(n, N, 2)` stage
# increments with the states as the columns of `x`, and its norm per
# column weighted by the tolerance `rowtol` of each row (see
# `weightedcolumnmax!`), `over` the host work of the rows over theirs
function gaussbatchresidual!(norms, residual, sys::TransientSystem, gc::GaussCoefficients, delta, x, lx, X,
        phi, junction, jwork, cwork, gwork, rhs, roundoff, colfloor, blockcol, rowtol, over, dwork,
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
    # the weighted norm per column over both stages
    weightedcolumnmax!(norms, over, residual, rowtol)
    # the floor, on the products' work, which the residual no longer needs
    gaussresidualroundoff!(roundoff, norms, over, sys, gc, residual, delta, x, X, phi, jwork, rhs, rowtol, cwork, gwork, colfloor,
        blockcol, dwork, rw)
    return nothing
end

# The residual `res` of each column weighted by the tolerance `tol` of
# each of its rows, `max_i |r_i|/tol_i`, into `norms`, which is at most
# one where every row is within its own tolerance, and the largest
# residual of a row over its tolerance, zero where none is, into `over`,
# both host vectors: a loop on the host, and reductions on a device.
function weightedcolumnmax!(norms::Vector{Float64}, over::Vector{Float64}, res::AbstractArray{<:Any,3}, tol)
    if res isa Array && tol isa Array
        @inbounds for col in axes(res, 2)
            m, e = 0.0, 0.0
            for i in axes(res, 3), row in axes(res, 1)
                r, t = abs(res[row, col, i]), tol[row, col, i]
                m = max(m, r/t)
                e = max(e, ifelse(r > t, r, 0.0))
            end
            norms[col], over[col] = m, e
        end
        return nothing
    end
    copyto!(norms, vec(mapreduce((r, t) -> abs(r)/t, max, res, tol; dims = (1, 3), init = 0.0)))
    copyto!(over, vec(mapreduce((r, t) -> abs(r) > t ? abs(r) : 0.0, max, res, tol; dims = (1, 3), init = 0.0)))
    return nothing
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
function gaussresidualroundoff!(roundoff, norms, over, sys::TransientSystem, gc::GaussCoefficients, residual, delta, x, X,
        phi, jwork, rhs, rowtol, bound, work, colfloor, blockcol, dwork, @nospecialize(rw::Union{Nothing, AbstractRationalWork}))
    cs, gs, ls, js, jr = sys.rowsums
    h = sys.h
    a2, a1 = maximum(abs, gc.ainv2)/h^2, maximum(abs, gc.ainv)/h
    isnothing(rw) ? fill!(blockcol, 0.0) : stackedscatterbound!(blockcol, rw, sys.gauss.coupling)
    gaussroundoffcolumns!(roundoff, 2a2*cs + 2a1*gs + ls, ls, jr, js*relationslope(sys, phi, dwork), delta, x, X, jwork,
        rhs, blockcol, colfloor, sys.backend)
    between = j -> norms[j] > 1 && over[j] <= roundoff[j]
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
    stagerounding!(bound, sys, a2, a1, delta, false, dsum, prod)
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

# The cheap bound per column, on the backend by reductions, into
# `roundoff` on the host through `colfloor`: `rounding` per unit of the
# increments, `lrows` of the state, `jrows` of the junction currents
# `jcur`, `jphase` of the stage values' phases, and the blocks' bound
# `blockcol` on the host (see `stackedscatterbound!`).
function gaussroundoffcolumns!(roundoff, rounding, lrows, jrows, jphase, delta, x, X, jcur, rhs, blockcol, colfloor, backend)
    colmax = a -> vec(maximum(abs, a; dims = (1, 3), init = 0.0))
    colfloor .= 8eps(Float64) .* (rounding .* colmax(delta) .+ lrows .* vec(maximum(abs, x; dims = 1, init = 0.0)) .+
        jrows .* colmax(jcur) .+ jphase .* colmax(X) .+ colmax(rhs))
    copyto!(roundoff, colfloor)
    roundoff .+= blockcol
    return nothing
end

# on the host a loop over each column, without the temporary arrays of
# the reductions, at every residual evaluation
function gaussroundoffcolumns!(roundoff, rounding, lrows, jrows, jphase, delta, x, X, jcur, rhs, blockcol, colfloor, ::CPU)
    @inbounds for col in axes(delta, 2)
        dm, xm, Xm, rm, jm = 0.0, 0.0, 0.0, 0.0, 0.0
        for row in axes(delta, 1)
            xm = max(xm, abs(x[row, col]))
            for i in 1:2
                dm = max(dm, abs(delta[row, col, i]))
                Xm = max(Xm, abs(X[row, col, i]))
                rm = max(rm, abs(rhs[row, col, i]))
            end
        end
        for i in 1:2, k in axes(jcur, 1)
            jm = max(jm, abs(jcur[k, col, i]))
        end
        roundoff[col] = 8eps(Float64)*(rounding*dm + lrows*xm + jrows*jm + jphase*Xm + rm) + blockcol[col]
    end
    return nothing
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

# the correction of a work, built on first use for the operator asked
function stagecorrectionwork!(rw::RationalWork, sys::TransientSystem, transposed::Bool)
    sc = rw.correction
    !isnothing(sc) && sc.transposed == transposed && return sc
    cp = sys.gauss.coupling
    backend = sys.backend
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    n = size(cp.Phost, 2) ÷ 2
    ports = cp.modulated
    np, nw = length(ports), length(cp.W)
    r = 2np
    m = size(rw.z, 2)
    St = devicesparse(sparse(transpose(cp.Shost[:, ports])), backend)
    Ucols, Wcols = allocate(n, np), allocate(2n, np*nw)
    copyto!(Ucols, Matrix(cp.Shost[:, ports]))
    transposed && copyto!(Wcols, reduce(hcat, [Matrix(transpose(W)) for W in cp.Whost]; init = zeros(2n, 0)))
    sc = StageCorrection(Ucols, St, Wcols, allocate(n, m, 2, r), allocate(n, m, 2, transposed ? np*nw : 0),
        transposed ? zeros(r, np, m, nw) : zeros(np, r, m, nw), LU{Float64,Matrix{Float64},Vector{Int}}[], transposed,
        allocate(n, m, 2), allocate(n, m, 2), allocate(np, m), zeros(np, m), zeros(r, m), allocate(m, r), zeros(m, r),
        zeros(r, r))
    rw.correction = sc
    return sc
end

# The pieces of the correction which depend on the factorization, built
# on the current one after a refresh: the solves `K = J*^(-1) U` with the
# products `W_ij K`, or for the transposed operator the solves
# `J*^(-T) W_ij'` with the products `U'` of them; then the step's `M`.
function stagefactors!(rw::RationalWork, sys::TransientSystem, bf::GaussBatchFactor, gc::GaussCoefficients, rc, zc,
        transposed::Bool)
    sys.gauss.pumped || return nothing
    cp = sys.gauss.coupling
    sc = stagecorrectionwork!(rw, sys, transposed)
    n = size(sc.rhs, 1)
    np, nw = length(cp.modulated), length(cp.W)
    if transposed
        for w in 1:nw, q in 1:np
            col = (w - 1)*np + q
            view(sc.rhs, :, :, 1) .= view(sc.Wcols, 1:n, col)
            view(sc.rhs, :, :, 2) .= view(sc.Wcols, n + 1:2n, col)
            gaussbatchtransform!(sc.sol, sc.rhs, bf, gc, rc, zc, true, sys)
            view(sc.X, :, :, :, col) .= sc.sol
            for i in 1:2
                stepmul!(sc.y, sc.St, stage(sc.sol, i))
                copyto!(sc.ys, sc.y)
                view(sc.Q, (i - 1)*np + 1:i*np, q, :, w) .= sc.ys
            end
        end
    else
        for col in 1:2np
            i, q = (col - 1) ÷ np + 1, (col - 1) % np + 1
            fill!(sc.rhs, 0)
            view(sc.rhs, :, :, i) .= view(sc.Ucols, :, q)
            gaussbatchtransform!(sc.sol, sc.rhs, bf, gc, rc, zc, false, sys)
            view(sc.K, :, :, :, col) .= sc.sol
            stack!(rw.dstack, sc.sol)
            for w in 1:nw
                stepmul!(sc.y, cp.W[w], rw.dstack)
                copyto!(sc.ys, sc.y)
                view(sc.Q, :, col, :, w) .= sc.ys
            end
        end
    end
    stagemultipliers!(rw, sys)
    return nothing
end

# The pieces of the correction at the step's weights, from those of the
# factorization: `M` per column, and for the transposed operator `K`.
# The term `w` of `W` is stage `i` and converted output `j`, the frozen
# operator's weight of which is zero, so its weight is the step's own.
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
            if sc.transposed
                view(Mc, :, (i - 1)*np + 1:i*np) .-= d .* view(sc.Q, :, :, c, w)
            else
                view(Mc, (i - 1)*np + 1:i*np, :) .-= d .* view(sc.Q, :, :, c, w)
            end
        end
        push!(sc.M, lu(Mc))
    end
    if sc.transposed
        fill!(sc.K, 0)
        for w in 1:nw
            d = weight(w)
            iszero(d) && continue
            i = (w - 1) ÷ nmod + 1
            view(sc.K, :, :, :, (i - 1)*np + 1:i*np) .-= d .* view(sc.X, :, :, :, (w - 1)*np + 1:w*np)
        end
    end
    return nothing
end

# the products on the block's rows at both stages of a stage pair `c`,
# `V' c`, or `U' c` for the transposed operator, into the host copy
# `sc.yh` as `(r, m)` stage major
function stagerows!(sc::StageCorrection, rw::RationalWork, cp::RationalCoupling, c)
    np, nw = length(cp.modulated), length(cp.W)
    nmod = nw ÷ 2
    if sc.transposed
        for i in 1:2
            stepmul!(sc.y, sc.St, stage(c, i))
            copyto!(sc.ys, sc.y)
            view(sc.yh, (i - 1)*np + 1:i*np, :) .= sc.ys
        end
        return sc.yh
    end
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

# the two stages of a `(n, N, 2)` array stacked into `(2n, N)`, back,
# and the two halves of a stacked array summed onto a `(n, N)` array
function stack!(s, a)
    n = size(a, 1)
    view(s, 1:n, :) .= stage(a, 1)
    view(s, n + 1:2n, :) .= stage(a, 2)
    return s
end
function unstack!(a, s)
    n = size(a, 1)
    stage(a, 1) .= view(s, 1:n, :)
    stage(a, 2) .= view(s, n + 1:2n, :)
    return a
end
function stacksum!(y, s)
    n = size(y, 1)
    y .+= view(s, 1:n, :) .+ view(s, n + 1:2n, :)
    return y
end

# the reflected waves of the rational parts at the stages scattered onto
# the blocks' rows, from the stage increments, the stage values and, with
# the state, the states: the stacked states `P_d delta + P_x X + P_z z`
# through the output terms at the stages' weights
function rationalsources!(rw::RationalWork, sys::TransientSystem, delta, X; withstate::Bool = true)
    cp = sys.gauss.coupling
    stack!(rw.dstack, delta)
    stack!(rw.xstack, X)
    stepmul!(rw.ystack, cp.Pd, rw.dstack)
    stepmul!(rw.ywork, cp.Px, rw.xstack)
    rw.ystack .+= rw.ywork
    if withstate
        stepmul!(rw.ywork, cp.Pz, rw.z)
        rw.ystack .+= rw.ywork
    end
    reflectedwaves!(rw, cp)
    return rw
end

# the transpose of the stages' rational coupling: the multipliers `mu`
# on the rows carried through the output terms to the stacked states and
# on to the stage unknowns, into `rw.source`
function rationalsourcestranspose!(rw::RationalWork, sys::TransientSystem, mu)
    cp = sys.gauss.coupling
    reflectedwavesbar!(rw, cp, mu)
    stepmul!(rw.sstack, cp.Pdt, rw.ystack)
    stepmul!(rw.tstack, cp.Pxt, rw.ystack)
    rw.sstack .+= rw.tstack
    unstack!(rw.source, rw.sstack)
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
# stencil's reach, one for a circuit without lines; and the ring that
# holds a window of them
lineprehistory(p::TransientProblem, h) = isempty(p.lines) ? 1 : ceil(Int, maximum(l -> l.delay, p.lines)/h) + 6
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

# the currents of every condition's drives at a time, on the host
function drivevalues!(hostvalues, problems, t)
    @inbounds for (j, p) in enumerate(problems), (k, d) in enumerate(p.drives)
        hostvalues[k, j] = d.current(t)
    end
    all(isfinite, hostvalues) || throw(ArgumentError(lazy"a source returned a nonfinite current at t = $(t) s."))
    return hostvalues
end

# The stepper of a batch: every buffer of a Gauss-Legendre step, the
# batch's factorizations, the Newton engine's closures and the state, so
# that the solve advances it step by step and a response replays a window
# from a checkpoint with the same code.
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
    junction::A
    cwork::A
    gwork::A
    xnew::M
    cv::M
    lx::M
    b1::M
    b2::M
    phi::A
    jwork::A
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
    # the projection's work, or nothing, and the rational blocks' work, or
    # nothing, as unions within the stepper's array types for the same
    # reason the stage's are
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
    colfloor = allocate(N)
    roundoff = (zeros(N), zeros(N))
    hostvalues = zeros(nd, N)
    values = tobackend(backend, zeros(nd, N))
    portwork = allocate(np, N)
    rw = isnothing(sys.gauss.coupling) ? nothing : rationalwork(p, backend, n, N)
    cell = RationalWorkCell{typeof(x), typeof(X)}(rw)
    # the tolerance of each row, set by the step, the host work of the
    # rows over theirs, and the blocks' part of the floor on the host
    rowtol = allocate(n, N, 2)
    over, blockcol = zeros(N), zeros(N)
    baseresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, phi, junction, jwork, cwork, gwork, rhs,
        roundoff[1], colfloor, blockcol, rowtol, over, dwork, cell.work)
    trialresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, trialphi, trialjunction, jwork, cwork,
        gwork, rhs, roundoff[2], colfloor, blockcol, rowtol, over, dwork, cell.work)
    # a new factorization of the frozen operator of the columns which
    # asked, and the stage correction of a pumped block rebuilt on it
    refresh! = columns -> (gaussbatchjacobian!(bf, sys, phi, cosphi, dwork, columns);
        isnothing(cell.work) || stagefactors!(cell.work, sys, bf, gc, rc, zc, false); nothing)
    solve! = (c, r) -> (gaussbatchtransform!(c, r, bf, gc, rc, zc); isnothing(cell.work) || stagecorrect!(c, cell.work, sys); (false, 0))
    accept! = mask -> (maskcolumns!(phi, trialphi, mask); maskcolumns!(junction, trialjunction, mask); nothing)
    pr = sys.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, N)
    nl = 2length(p.lines)
    far, readscale, sqrtz = linetables(p, backend)
    npre = lineprehistory(p, sys.h)
    return GaussStepper(sys, problems, N, x, v, X, delta, lastdelta, residual, trial, trialresidual, correction, rhs,
        junction, cwork, gwork, xnew, cv, lx, b1, b2, phi, jwork, dwork, cosphi, rc, zc,
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

# the rate of a stepper's state at time `t` read along the algebraic
# directions into `dst`, with the work `rw` of `ratereadwork`: the
# drives' rate by the central difference of `delta`, the lines' forced
# currents' rate by the same difference of their histories
function readstepper!(dst, st::GaussStepper, @nospecialize(rw::Union{Nothing, RateReadWork}), t, delta)
    sys = st.sys
    if !isnothing(sys.projection)
        linerate = if !isempty(sys.problem.lines)
            lr = Array(linevalues!(st, t + delta))
            lr .= (lr .- Array(linevalues!(st, t - delta))) ./ (2delta)
        else
            nothing
        end
        drivedotz!(rw.bdotz, sys.projection, st.problems, t, delta, rw.hv1, rw.hv2, linerate)
    end
    return readrate!(dst, st.v, st.x, sys, rw)
end

# the phases of the projected junctions at a state, for the record
function projectedphases!(st::GaussStepper, x)
    pr = st.sys.projection
    isempty(pr.pj) || stepmul!(st.pw.phip, pr.RJp, x)
    return st.pw.phip
end

# The projection of every condition's endpoint onto the algebraic
# constraints along the projected directions (see `projectendpoint!`):
# Newton on the coefficients `alpha` of `Z` per condition, the residual
# `Z' (L x + J(x) - b(t))` and its Jacobians `Z' (L + J'(x)) Z` on the
# host from two small products on the backend, to the roundoff of the
# terms the constraint balances, the predictor moved with the state. The
# step's tolerance is set in the stage equations' units and is far
# looser for the constraint, so the drift of the endpoints would
# accumulate below it. The rate along `Z`
# is left to the reading at the read-outs, except where a block is on a
# direction, when it is the derivative of the cubic through the state,
# the stages and the projected endpoint. Leaves the projected junctions'
# phases at the endpoint in the work buffer.
function gaussproject!(st::GaussStepper, t, step)
    sys, pw = st.sys, st.pw
    pr = sys.projection
    gc = sys.gauss.coefficients
    h = sys.h
    drivevalues!(st.hostvalues, st.problems, t)
    isempty(sys.problem.lines) || linevalues!(st, t)
    linedrive = isempty(sys.problem.lines) ? nothing : Array(st.linevalues)
    if !isempty(pr.directions)
        # the drive along the constraints' rows, and the magnitudes of its
        # terms for the residual's floor
        gb, gbabs = pw.gb, pw.gbabs
        mul!(gb, pr.Zcinj, st.hostvalues)
        gb .+= pr.Zcconstant
        pw.hostabs .= abs.(st.hostvalues)
        mul!(gbabs, pr.Zcinjabs, pw.hostabs)
        gbabs .+= abs.(pr.Zcconstant)
        if !isnothing(linedrive)
            mul!(gb, pr.Zcline, linedrive, 1.0, 1.0)
            mul!(gbabs, pr.Zclineabs, abs.(linedrive), 1.0, 1.0)
        end
        if !isnothing(st.rw)
            waves = restingwaves(sys, st.rw.z, t)
            mul!(gb, pr.Zcblock, waves, 1.0, 1.0)
            mul!(gbabs, pr.Zcblockabs, abs.(waves), 1.0, 1.0)
        end
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
        gq = pw.gqdrive
        mul!(gq, pr.Qinj, st.hostvalues)
        gq .+= pr.Qconst
        isnothing(linedrive) || mul!(gq, pr.Qline, linedrive, 1.0, 1.0)
        isnothing(st.rw) || mul!(gq, pr.Qblock, restingwaves(sys, st.rw.z, t), 1.0, 1.0)
        relationinto!(pw.current, pr.relationsp, pw.hphi)
        pw.current .*= pr.lmoljp
        endpointread!(st.v, st.xnew, pr, pw, pw.current, gq)
        projectedphases!(st, st.xnew)
    end
    return st
end

# the stepper set to a state and the stage increments to start from, its
# stage matrices assembled at that state, factorized anew with `fresh`
function setstate!(st::GaussStepper, x, v, lastdelta; fresh::Bool = false)
    copyto!(st.x, x)
    copyto!(st.v, v)
    isnothing(lastdelta) ? fill!(st.lastdelta, 0) : copyto!(st.lastdelta, lastdelta)
    st.X .= st.x
    stepmul!(stage(st.phi, 1), st.sys.RJ, stage(st.X, 1)); stepmul!(stage(st.phi, 2), st.sys.RJ, stage(st.X, 2))
    gaussbatchjacobian!(st.bf, st.sys, st.phi, st.cosphi, st.dwork; fresh)
    isnothing(st.rw) || stagefactors!(st.rw, st.sys, st.bf, st.sys.gauss.coefficients, st.rc, st.zc, false)
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
        simplified = true, accept! = st.accept!, roundoff = st.roundoff)
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
    return TransientBatchSolution(problems, sys.method, h, savedtimes, out.voltage, out.incident, out.outgoing,
        recordedarray(out.phases), recordedarray(out.endphases), recordedarray(out.endrates), recordedarray(out.linewaves),
        recordedarray(out.history), recordedarray(out.flux), recordedarray(out.rate), recordedarray(out.stages),
        recordedcheckpoints(out.checkpoints), isempty(sys.problem.lines) ? nothing : out.finalwaves,
        isnothing(sys.gauss.coupling) ? nothing : out.finalstates, out.initialflux, out.initialrate, out.finalflux, out.finalrate,
        recordedarray(out.blockstates), out.initialwaves, out.initialstates, stats)
end

# The initial states of every condition of a batch, read on the host and
# checked against the algebraic equations before any chunk steps, so that
# a state a condition refuses is refused from the caller's task before a
# thread has begun stepping another: the fluxes and the rates, the
# history of the line waves before the start of every condition and its
# last column, the waves at the start, and the block states. The
# currents the lines force at the start and their rates, read from each
# state's own history at its own step, are what the check reads, by the
# difference it reads the rate of the drives with.
function batchinitial(sys::TransientSystem, problems, initialstates, t0, rtol, atol)
    p = sys.problem
    N = length(problems)
    n = length(p)
    h = sys.h
    nl = 2length(p.lines)
    npre = lineprehistory(p, h)
    x0, v0, w0, tailh = zeros(n, N), zeros(n, N), zeros(nl, N), zeros(nl, npre, N)
    q0, qdot = zeros(nl, N), zeros(nl, N)
    delta = 1e-3*h
    for (j, state) in enumerate(initialstates)
        xj, vj = state.flux, state.rate
        (length(xj) == n && length(vj) == n) || throw(DimensionMismatch(lazy"the state has $(n) entries; use transientstate."))
        (all(isfinite, xj) && all(isfinite, vj)) || throw(ArgumentError("the initial state must be finite."))
        x0[:, j] .= Float64.(collect(xj))
        v0[:, j] .= Float64.(collect(vj))
        wj = initialwaves(state, p, h, npre)
        all(isfinite, wj) || throw(ArgumentError("the initial line waves must be finite."))
        tailh[:, :, j] .= wj
        w0[:, j] .= view(wj, :, npre)
        if nl > 0
            q0[:, j] .= lineforcing(p, initialarrivals(state, p))
            qdot[:, j] .= (lineforcing(p, initialarrivals(state, p, delta)) .-
                lineforcing(p, initialarrivals(state, p, -delta))) ./ (2delta)
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
# `chunkerror`)
function runchunks(f, chunks)
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
# conditions: every step advanced and saved. The record is whichever
# arrays of `out` are there to be filled. Returns the chunk's statistics.
function gaussbatchrun!(out, sys::TransientSystem, problems, init, times, saveevery,
        rtol, atol, iterations, bf::GaussBatchFactor, conditions = 1:length(problems))
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
    readout!(dst, t) = readstepper!(dst, st, rw, t, delta)
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
            copyto!(view(out.rate, :, saved, :), saved == 1 ? st.v : vread)
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
# and `state!(x, v, k)` and `stages!(buffer, k)` filling the states and
# the stage increments where a response reads them. The readers are held
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
function responsewindows(sol::TransientBatchSolution, problems, sys::TransientSystem, forward::Bool;
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
    N = length(problems)
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
    rw = npj > 0 ? ratereadwork(sys, sys.backend, length(sys.problem), N, N) : nothing
    vread = npj > 0 ? KernelAbstractions.zeros(sys.backend, Float64, length(sys.problem), N) : nothing
    delta = ratedelta(sys)
    readrates! = (column, t) -> begin
        readstepper!(vread, st, rw, t, delta)
        stepmul!(rw.pw.phip, pr.RJp, vread)
        copyto!(view(rbuffer, :, column, :), rw.pw.phip)
        nothing
    end
    # the states of the window and the stage increments of its steps, for
    # the responses which read them
    n = length(sys.problem)
    xbuffer = withstates ? KernelAbstractions.zeros(sys.backend, Float64, n, K + 2, N) : nothing
    vbuffer = withstates ? KernelAbstractions.zeros(sys.backend, Float64, n, K + 2, N) : nothing
    dbuffer = withstates ? KernelAbstractions.zeros(sys.backend, Float64, n, 2, K + 1, N) : nothing
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
            npj > 0 && readrates!(1, sol.times[kstart])
            if withstates
                copyto!(view(xbuffer, :, 1, :), st.x)
                copyto!(view(vbuffer, :, 1, :), st.v)
            end
            for k in kstart:kend
                advance!(st, sol.times[k], sol.times[k + 1], k)
                view(buffer, :, 1, k - kstart + 2, :) .= stage(st.phi, 1)
                view(buffer, :, 2, k - kstart + 2, :) .= stage(st.phi, 2)
                npj > 0 && copyto!(view(ebuffer, :, k - kstart + 2, :), st.pw.phip)
                npj > 0 && readrates!(k - kstart + 2, sol.times[k + 1])
                if withstates
                    view(dbuffer, :, 1, k - kstart + 1, :) .= stage(st.delta, 1)
                    view(dbuffer, :, 2, k - kstart + 1, :) .= stage(st.delta, 2)
                    copyto!(view(xbuffer, :, k - kstart + 2, :), st.x)
                    copyto!(view(vbuffer, :, k - kstart + 2, :), st.v)
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
# dimension of one: the Gauss-Legendre responses run on the batch form
# alone, and take a solution through this
function batchof(sol::TransientSolution)
    cps = mapcheckpoints(expand, sol.checkpoints)
    return TransientBatchSolution([sol.problem], sol.method, sol.dt, sol.times, expand(sol.voltage), expand(sol.incident),
        expand(sol.outgoing), expand(sol.phases), expand(sol.endphases), expand(sol.endrates), expand(sol.linewaves),
        expand(sol.history), expand(sol.flux), expand(sol.rate), expand(sol.stages), cps, expand(sol.finalwaves),
        expand(sol.finalstates), expand(sol.initialflux),
        expand(sol.initialrate), expand(sol.finalflux), expand(sol.finalrate), expand(sol.blockstates),
        expand(sol.initialwaves), expand(sol.initialstates), sol.stats)
end
# the arrays of a response of a batch of one without the trailing
# dimension of one, as the response of the solution, and the reverse
dropcondition(a) = isnothing(a) ? nothing : reshape(a, size(a)[1:end-1]...)
addcondition(a) = isnothing(a) ? nothing : reshape(a, size(a)..., 1)
# the solution of a batch of one, which is the only batch another rule
# than Gauss-Legendre makes (see batchof)
function singlecondition(b::TransientBatchSolution)
    length(b) == 1 || throw(ArgumentError("a batch of conditions steps under GaussLegendre()."))
    return unbatch(b)
end

# a batch of one is an ordinary solution
function unbatch(b::TransientBatchSolution)
    squeeze = a -> isnothing(a) ? nothing : reshape(a, size(a)[1:end-1]...)
    cp = mapcheckpoints(squeeze, b.checkpoints)
    return TransientSolution(b.problems[1], b.method, b.dt, b.times, squeeze(b.voltage), squeeze(b.incident),
        squeeze(b.outgoing), squeeze(b.phases), squeeze(b.endphases), squeeze(b.endrates), squeeze(b.linewaves), squeeze(b.history), squeeze(b.flux), squeeze(b.rate), squeeze(b.stages), cp,
        squeeze(b.finalwaves), squeeze(b.finalstates), vec(b.initialflux),
        vec(b.initialrate), vec(b.finalflux), vec(b.finalrate), squeeze(b.blockstates), squeeze(b.initialwaves), squeeze(b.initialstates), b.stats)
end
