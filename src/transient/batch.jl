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
    drives = TransientDrive[]
    rows, cols, vals = Int[], Int[], Float64[]
    replaced = Set{Int}()
    for (k, source) in enumerate(sources)
        source isa TransientSource || throw(ArgumentError("sources must contain TransientSource objects."))
        if source.target isa Int
            q = findfirst(port -> port.number == source.target, psc.ports)
            isnothing(q) && throw(ArgumentError(lazy"there is no port $(source.target)."))
            n1, n2 = p.portpositive[q], p.portnegative[q]
            push!(drives, TransientDrive(q, source.current))
        else
            c = get(psc.componentnamedict, source.target, 0)
            (c > 0 && psc.componenttypes[c] == :I) || throw(ArgumentError(
                lazy"$(source.target) does not name a CurrentSource of the circuit."))
            push!(replaced, c)
            n1, n2 = psc.nodeindices[2, c] - 1, psc.nodeindices[1, c] - 1
            push!(drives, TransientDrive(0, source.current))
        end
        n1 > 0 && (push!(rows, n1); push!(cols, k); push!(vals, 1.0))
        n2 > 0 && (push!(rows, n2); push!(cols, k); push!(vals, -1.0))
    end
    injection = sparse(rows, cols, vals, length(p), length(drives))
    constantcurrent = zeros(length(p))
    vvn = p.matrices.vvn
    for c in psc.currentsources
        c in replaced && continue
        n1, n2 = psc.nodeindices[2, c] - 1, psc.nodeindices[1, c] - 1
        n1 > 0 && (constantcurrent[n1] += vvn[c])
        n2 > 0 && (constantcurrent[n2] -= vvn[c])
    end
    return TransientProblem(psc, p.matrices, p.Nnodal, p.Naux, p.Lscale, p.coupledbranches,
        p.floatingcomponents, p.gaugeindices, p.inertialess, p.algebraic, p.directions, p.constraints,
        injection, drives, constantcurrent, p.portpositive, p.portnegative, p.portimpedances, p.portconductances, p.blocks, p.lines,
        p.relations)
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

`stats` counts the work of one chunk of conditions: a host batch split
across threads reports the counters of the chunk which worked hardest,
as a batch on one chunk reports those of the condition which converged
worst. The count of Newton corrections therefore falls as the chunks get
smaller, while the results do not change.
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
    cp = isnothing(b.checkpoints) ? nothing : (; every = b.checkpoints.every, flux = slice(b.checkpoints.flux),
        rate = slice(b.checkpoints.rate), increment = slice(b.checkpoints.increment), states = slice(b.checkpoints.states),
        waves = slice(b.checkpoints.waves))
    return TransientBatchSolution(b.problems[js], b.method, b.dt, b.times, slice(b.voltage), slice(b.incident),
        slice(b.outgoing), slice(b.phases), slice(b.endphases), slice(b.endrates), slice(b.linewaves), slice(b.history), slice(b.flux), slice(b.rate), slice(b.stages), cp, slice(b.initialflux),
        slice(b.initialrate), slice(b.finalflux), slice(b.finalrate), slice(b.blockstates), slice(b.initialwaves), slice(b.initialstates), b.stats)
end
function Base.getindex(b::TransientBatchSolution, j::Integer)
    1 <= j <= length(b) || throw(BoundsError(b, j))
    slice = a -> isnothing(a) ? nothing : selectdim(a, ndims(a), j)
    cp = isnothing(b.checkpoints) ? nothing : (; every = b.checkpoints.every, flux = slice(b.checkpoints.flux),
        rate = slice(b.checkpoints.rate), increment = slice(b.checkpoints.increment), states = slice(b.checkpoints.states),
        waves = slice(b.checkpoints.waves))
    return TransientSolution(b.problems[j], b.method, b.dt, b.times, slice(b.voltage), slice(b.incident),
        slice(b.outgoing), slice(b.phases), slice(b.endphases), slice(b.endrates), slice(b.linewaves), slice(b.history), slice(b.flux), slice(b.rate), slice(b.stages), cp, copy(view(b.initialflux, :, j)),
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
# the transpose from the one factorization.
struct GaussBatchFactor{J, F, X, B, S, I, IT, R}
    ncolumns::Int
    jacobians::J
    factors::F
    X::X
    B::B
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
    # the entries of the rational blocks in the stage operator, on the
    # backend and on the host: those of the system when they never
    # change, and this factor's own for a pumped block, whose stage
    # operator is refreshed every step by the stepper which owns the
    # factor (see refreshstageoperator!), so that batches stepping at
    # once on one system do not write over each other
    rationalvals::R
    hostvals::Vector{ComplexF64}
end

function gaussbatchfactor(sys::TransientSystem, ncolumns::Int; nrhs::Int = 1)
    g = sys.gauss
    backend = sys.backend
    nnzj = nnz(sys.jacobian)
    rationalvals = g.pumped ? copy(g.rationalvals) : g.rationalvals
    hostvals = g.pumped ? copy(g.hostvals) : g.hostvals
    if backend isa CPU
        pattern = g.cjacobian
        jacobians = [SparseMatrixCSC(size(pattern)..., SparseArrays.getcolptr(pattern), rowvals(pattern),
            zeros(ComplexF64, nnzj)) for _ in 1:ncolumns]
        # the assembly scratch belongs to the factor, not to the system,
        # so that several batches of one circuit may step at once
        return GaussBatchFactor(ncolumns, jacobians, Vector{Any}(nothing, ncolumns), nothing, nothing,
            zeros(Float64, nnzj), UnitRange{Int}[],
            nothing, nothing, sys.symmetric, Any[], nothing, nothing, nothing, nothing, Ref(true), rationalvals, hostvals)
    end
    n = size(sys.jacobian, 1)
    limit = uniformbatchlimit(nrhs)
    chunks = [first:min(first + limit - 1, ncolumns) for first in 1:limit:ncolumns]
    nzval = KernelAbstractions.zeros(backend, ComplexF64, nnzj, ncolumns)
    X = KernelAbstractions.zeros(backend, ComplexF64, n, nrhs, ncolumns)
    B = KernelAbstractions.zeros(backend, ComplexF64, n, nrhs, ncolumns)
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
    return GaussBatchFactor(ncolumns, nzval, Vector{Any}(nothing, length(chunks)), X, B, nothing, chunks, rowptr, colind,
        sys.symmetric, Vector{Any}(nothing, length(chunks)), trowptr, tcolind, tperm, tnzval, Ref(true), rationalvals, hostvals)
end

# the transposed factors of a device batch from its current values, built
# or refreshed when an adjoint asks and the forward factor has been
# refreshed since
function gaussbatchtransposed!(bf::GaussBatchFactor, sys::TransientSystem)
    bf.tstale[] || return bf
    bf.tnzval .= view(bf.jacobians, bf.tperm, :)
    for (c, chunk) in enumerate(bf.chunks)
        if isnothing(bf.tfactors[c])
            bf.tfactors[c] = _cudss_sweep(bf.trowptr, bf.tcolind, bf.tnzval[:, chunk], bf.X[:, :, chunk], bf.B[:, :, chunk];
                sys.factorization.kwargs...)
        else
            copyto!(bf.tfactors[c].nzval, view(bf.tnzval, :, chunk))
            _cudss_sweeprefactorize!(bf.tfactors[c])
        end
    end
    bf.tstale[] = false
    return bf
end

# the stage matrices at the stage phases of every column, assembled per
# column by the plan from the mean cosine of the two stages, and
# factorized, or refactorized on the analyses of the first time
function gaussbatchjacobian!(bf::GaussBatchFactor, sys::TransientSystem, phi, cosphi, dwork)
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
    if sys.backend isa CPU
        for j in 1:ncolumns
            A = bf.jacobians[j]
            assemblerealjacobian!(bf.scratch, sys.plan, view(cosphi, :, j))
            nonzeros(A) .= bf.scratch .+ im .* g.imvals .+ bf.rationalvals
            bf.factors[j] = isnothing(bf.factors[j]) ? factorize(sys.factorization, A) :
                refactorize!(sys.factorization, bf.factors[j], A)
        end
        return bf
    end
    nzval = bf.jacobians
    # every condition's assembly launched, then one synchronization
    for j in 1:ncolumns
        assemblerealjacobian!(nonzeros(sys.jacobian), sys.plan, view(cosphi, :, j); synchronize = false)
        view(nzval, :, j) .= nonzeros(sys.jacobian) .+ im .* g.imvals .+ bf.rationalvals
    end
    KernelAbstractions.synchronize(sys.backend)
    bf.tstale[] = true
    rowptr, colind = bf.rowptr, bf.colind
    for (c, chunk) in enumerate(bf.chunks)
        if isnothing(bf.factors[c])
            bf.factors[c] = _cudss_sweep(rowptr, colind, nzval[:, chunk], bf.X[:, :, chunk], bf.B[:, :, chunk];
                sys.factorization.kwargs...)
        else
            copyto!(bf.factors[c].nzval, view(nzval, :, chunk))
            _cudss_sweeprefactorize!(bf.factors[c])
        end
    end
    return bf
end

# the solve of every condition's block of right hand sides, or of the
# transposed operator's: the columns of `R` are `nrhs` per condition,
# condition after condition
function gaussbatchsolve!(Z::AbstractMatrix, bf::GaussBatchFactor, R::AbstractMatrix, transposed::Bool = false, sys = nothing)
    n = size(R, 1)
    nrhs = size(R, 2) ÷ bf.ncolumns
    transposed = transposed && !bf.symmetric
    if isnothing(bf.X)
        for j in 1:bf.ncolumns
            block = (j - 1)*nrhs + 1:j*nrhs
            matrixsolve!(view(Z, :, block), transposed ? transpose(bf.factors[j]) : bf.factors[j], view(R, :, block))
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
    if allsinusoidal(r)
        reshape(work, nj, ncolumns, N) .*= sys.lmolj .* reshape(cos.(phi), nj, 1, N)
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

# The exact solve of the stage equations of every condition linearized at
# its recorded stage phases, `J d = r` with the true stage Jacobian whose
# two stiffness blocks differ, or its transpose: a fixed point iteration
# on the frozen complex operators, `d += [J*]^{-1} (r - J d)`, contracting
# as the simplified Newton of the step does. `d` and `r` are `(n, ndir*N, 2)`
# with the directions of a condition contiguous; `phi` is `(nj, N, 2)`.
# A column iterates until it has converged and no further, so that its
# solution does not depend on the columns beside it, and a response is
# the same for a condition whatever batch or chunk it is stepped in.
function gaussbatchstagesolve!(d, r, sys::TransientSystem, gc::GaussCoefficients, bf::GaussBatchFactor, phi,
        transposed::Bool, res, dd, cwork, gwork, lwork, jwork, dwork, rc, zc,
        @nospecialize(rw::Union{Nothing, AbstractRationalWork}); rtol = 1e-12, iterations = 100)
    h = sys.h
    fill!(d, 0)
    # the convergence per column: each direction of each condition against
    # its own right hand side, so a weak direction is not left at the
    # tolerance of a strong one, down to the roundoff floor of the terms
    # the residual sums
    colmax = a -> vec(Array(maximum(abs, a; dims = (1, 3))))
    scale = max.(colmax(r), floatmin(Float64))
    cmax, gmax = maximum(abs, gc.ainv2)/h^2, maximum(abs, gc.ainv)/h
    # the rounding of the terms the residual sums, per unit of the
    # increment: the row or column sums of each matrix, whichever is
    # larger, since the adjoint multiplies by the transposes, at the
    # rule's weights, and the junction stamp at the steepest slope of the
    # relations at the stages. A product's own size would understate it
    # along an algebraic direction, where `C d` cancels. The eight units
    # of roundoff allow for the order of the summation; they are not a
    # bound for rows of every length.
    cs, gs, ls, js = sys.rowsums
    rounding = 2cmax*cs + 2gmax*gs + ls + js*relationslope(sys, phi, dwork)
    # the contraction of the first correction, the measure of how well the
    # frozen operator fits this step's two stiffnesses, returned for the
    # caller to decide whether to refresh it
    first = copy(scale)
    contraction = 0.0
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
        # not cover. `|S| |y|` keeps every scatter entry with the state
        # it multiplies, where a norm of each would pair the largest row
        # of one with the largest entry of the other whatever rows they
        # are in, and the poles of a fit span decades; the bound is per
        # column, so a strong direction does not lift a weak one's
        # threshold and the tangent and the adjoint stay transposes.
        blockfloor = 0.0
        if !isnothing(rw)
            transposed ? rationalsourcestranspose!(rw, sys, d) : rationalsources!(rw, sys, d, d; withstate = false)
            res .-= rw.source
            if !transposed
                reflectedwavesbound!(rw, sys.gauss.coupling)
                blockfloor = sys.gauss.coupling.nstates*eps(Float64) .* colmax(rw.sbound)
            end
        end
        floor = 8eps(Float64) .* (scale .+ rounding .* colmax(d) .+ colmax(res)) .+ blockfloor
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
        current = colmax(res)
        iteration == 2 && (contraction = maximum(current ./ first))
        active = current .> max.(rtol .* scale, floor)
        any(active) || return contraction
        gaussbatchtransform!(dd, res, bf, gc, rc, zc, transposed, sys)
        isnothing(rw) || stagecorrect!(dd, rw, sys)
        # the columns which have converged are left as they are
        all(active) || (dd .*= reshape(tobackend(sys.backend, Float64.(active)), 1, :, 1))
        d .+= dd
    end
    error("the stage equations of a Gauss-Legendre step did not converge in the tangent or the adjoint: the two stage stiffnesses differ too much from their mean; reduce dt.")
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
operator of a step carries the block's converted coupling at the mean
of the two stages' weights (see [`refreshstageoperator!`](@ref)); what
is left, the difference of each stage's weights from the mean, acts on
the block's port rows alone, so the true stage operator is the frozen
one plus a correction of rank twice the block's ports,
`J = J* + U V'`, with `U` the scatter onto the block's rows at each
stage and `V'` the difference of the weights times the output terms on
the stacked states from the stage unknowns. It is taken exactly by the
Woodbury identity: `K = J*^(-1) U` by one solve per column of `U` on
the step's factorization, `M = I + V' K` per column of the batch, and a
solve `c = J*^(-1) r` is corrected to `c - K M^(-1) V' c`; the
transposed operator, which the adjoint solves, has `U` and `V`
exchanged. A converted coupling as large as the unconverted one, an
amplifier's, would not converge on the frozen operator alone.
"""
mutable struct StageCorrection{S, K, R, Y, Z}
    # `V'` per stage on the backend, `(nports, 2n)`, its rows as the
    # columns of a `(2n, r)` array, and the transpose of the scatter
    # with its columns as an `(n, nports)` array
    Vt::Vector{S}
    Vrows::Z
    St::S
    Ucols::Y
    # `K` on the backend, `(n, m, 2, r)`, the small matrices `M`
    # factorized per column on the host, and whether they are of the
    # transposed operator
    K::K
    M::Vector{LU{Float64,Matrix{Float64},Vector{Int}}}
    transposed::Bool
    # scratch: a right hand side and a solution of the stage transform,
    # a product on the block's rows, `(nports, m)`, its host copy, and
    # the multipliers of the columns, `(m, r)`, on both
    rhs::R
    sol::R
    y::Y
    ys::Matrix{Float64}
    yh::Matrix{Float64}
    z::Y
    zh::Matrix{Float64}
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
    # the exact stage solve of a pumped block, built per step, or nothing
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
    relations::JunctionRelations{Matrix{Float64}, Vector{Bool}}
    Q::SparseMatrixCSC{Float64, Int}
end
function endpointwork(cp, sys::TransientSystem, p::TransientProblem, N::Int, nobj::Int; forcing::Bool)
    n, nj, nc = cp.n, length(sys.lmolj), length(cp.names)
    pr = sys.projection
    RJ, _ = junctionincidence(p)
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
    # the stage factorizations and whether the stage operator is stale
    bf::B
    stale::Bool
    # the projection of the endpoint: its work, the injections along its
    # rows and the perturbation's endpoint work; the perturbation's stage
    # work, and the direct term of the port waves along the components
    # which are port terminations, read once per time
    pw::Union{Nothing, ProjectionWork{M, Matrix{Float64}, M}}
    Zti::Union{Nothing, SparseMatrixCSC{Float64, Int}}
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
    rwt = isempty(sys.gauss.rational) ? nothing : rationalwork(p, backend, n, m)
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
    Zti = isnothing(pr) ? nothing : pr.Zth*injh
    # the direction's current along the rows of the endpoint reading
    Qinji = isnothing(pr) ? nothing : endpointinjection(pr, injh)
    return GaussTangentWork(sys, N, ndir, nq, nt, staged, constant, storing, injh, tp, perturbation, withstates, injection,
        dI, dS, zeros(nq, 1, ndir), dx, dv, dxnew, lx, cdv, work, r, d, res, dd, cwork, gwork, lwork,
        dwaves, ddlinevalues, dlinerates, linework, stencil, dstencil, linetables(p, backend), 0.0,
        rwt, zerostages, fullstages, jwork, phi, dwork, cosphi, rc, zc, stagecurrent, stageinjection, injectionall,
        voltage, incident, outgoing, coefficients, portwork, portmap, directwork, directall, outwork, reading, xs0, vs0,
        bf, true, pw, Zti, Qinji, pwork, pend, dterm, nothing)
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
    n, np, N, ndir, nt = length(p), length(p.portimpedances), w.N, w.ndir, w.nt
    m = ndir*N
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
    w.stale = true
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
    outs = w.storing ? (view(w.voltage, :, k, :), view(w.incident, :, k, :), view(w.outgoing, :, k, :)) : w.outwork
    for s in 1:3
        cv, cd = w.coefficients[s]
        outs[s] .= cv .* w.portwork .+ cd .* w.directall
    end
    isnothing(w.dterm) || adddirectterm!(outs, w.dterm, sol, sol.problems, k)
    w.storing || outputsink(k, w.outwork[1], w.outwork[2], w.outwork[3])
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
        Ms = projectionmatrices(pr, pw.hphi)
        constrainttangent!(pw, pr, w.dxnew)
        dIz = repeat(w.Zti*view(w.dIh, :, kcol, :), 1, N)
        nl2 > 0 && (dIz .+= pr.Ztline*Array(w.ddlinevalues))
        isnothing(rwt) || (dIz .+= pr.Ztblock*restingwaves(sys, rwt.z, sol.times[k + 1]))
        isnothing(pend) || (dIz .+= pr.Zth*pend.F)
        pw.g .-= dIz
        for j in 1:N
            cols = (j - 1)*ndir + 1:j*ndir
            view(pw.alpha, :, cols) .= -(Ms[j] \ view(pw.g, :, cols))
        end
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
# phases when its first correction did not contract as the step's own did
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
        refreshstageoperator!(rwt, sys, bf) && (w.stale = true)
        stage(w.fullstages, 1) .= w.dx
        stage(w.fullstages, 2) .= w.dx
        rationalsources!(rwt, sys, w.zerostages, w.fullstages)
        w.r .+= rwt.source
    end
    w.stale && gaussbatchjacobian!(bf, sys, w.phi, w.cosphi, w.dwork)
    isnothing(rwt) || stagecorrection!(rwt, sys, bf, gc, w.rc, w.zc, false)
    contraction = gaussbatchstagesolve!(w.d, w.r, sys, gc, bf, w.phi, false, w.res, w.dd, w.cwork, w.gwork, w.lwork,
        w.jwork, w.dwork, w.rc, w.zc, rwt)
    w.stale = nj > 0 && contraction > 0.25
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
    parts = Vector{Any}(undef, length(chunks))
    if length(chunks) == 1
        parts[1] = Base.invokelatest(gaussbatchtangent, ws[1], sol, currents, initial, outputsink, outputs(1:N))
    else
        try
            Base.Threads.@sync for (c, ch) in enumerate(chunks)
                Base.Threads.@spawn parts[c] = Base.invokelatest(gaussbatchtangent, ws[c], sol[ch], currents,
                    conditionslice(initial, ch), outputsink, outputs(ch))
            end
        catch err
            throw(chunkerror(err))
        end
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
# between calls, and of a call only the views of its result: what a
# call hands over and a workspace only reads while it runs, the host
# currents of a tangent and the sink of an adjoint with whatever it
# captured, the noise's accumulator among them, is released when the
# call returns, on an error as well, in the reset as in the steps, so
# that its lifetime is the call's.
function tangentrelease!(w::GaussTangentWork)
    w.dIh = zeros(size(w.dIh, 1), 0, size(w.dIh, 3))
    return w
end

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
    # the stage factorizations and whether the stage operator is stale
    bf::B
    stale::Bool
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
    rwa = isempty(sys.gauss.rational) ? nothing : rationalwork(p, backend, n, m)
    xextra = isnothing(rwa) ? nothing : allocate(n, m)
    pr = sys.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, m)
    vrecbar = !isnothing(pr) && pr.cubic ? allocate(n, m) : nothing
    return GaussAdjointWork(sys, N, nobj, nq, nt, quantity, storing, injh, tp, perturbation, withstates, injectiont, portmapt,
        w, cv, cd, directwork, targetwork, objectivework, currents, CurrentRing(ringslots, ringcolumn, nothing, nothing),
        xbar, vbar, work, lmu, cmu, wst, mu, res, dd, cwork, gwork, lwork, jwork, phi, dwork, cosphi, rc, zc,
        portwork, objectivestate, transposing, xf, wf, bf, true,
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
    backend = sys.backend
    n, N, nq, nt, m = length(p), w.N, w.nq, w.nt, w.nobj*w.N
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
    w.stale = true
    return w
end

# the ring's callbacks released when the call returns (see
# tangentrelease!): the sink of the call and what it captured go with it
function adjointrelease!(w::GaussAdjointWork)
    w.ring = CurrentRing(w.ring.slots, w.ring.column, nothing, nothing)
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
    Ms = projectionmatrices(pr, pw.hphi)
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
    for j in 1:N
        cols = (j - 1)*nobj + 1:j*nobj
        view(pw.alpha, :, cols) .= transpose(Ms[j]) \ view(pw.g, :, cols)
    end
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
# transposed stage operator, refreshed at the recorded phases when its
# first correction did not contract, and the multipliers' pull on the
# currents the stages read, the components, the lines and the state
function adjointstep!(w::GaussAdjointWork, sol::TransientBatchSolution, k::Int, window, @nospecialize(stagesink))
    sys, rwa, pwork, bf = w.sys, w.rwa, w.pwork, w.bf
    gc = sys.gauss.coefficients
    h = sys.h
    nt, N, nobj = w.nt, w.N, w.nobj
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
        refreshstageoperator!(rwa, sys, bf) && (w.stale = true)
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
    w.stale && gaussbatchjacobian!(bf, sys, w.phi, w.cosphi, w.dwork)
    isnothing(rwa) || stagecorrection!(rwa, sys, bf, gc, w.rc, w.zc, true)
    contraction = gaussbatchstagesolve!(w.mu, w.wst, sys, gc, bf, w.phi, true, w.res, w.dd, w.cwork, w.gwork, w.lwork,
        w.jwork, w.dwork, w.rc, w.zc, rwa)
    w.stale = nj > 0 && contraction > 0.25
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
    parts = Vector{Any}(undef, length(chunks))
    if length(chunks) == 1
        parts[1] = Base.invokelatest(gaussbatchadjoint, ws[1], sol, wh, single, sink, stagesink, outputview(currents, columns(1:N)))
    else
        try
            Base.Threads.@sync for (c, ch) in enumerate(chunks)
                Base.Threads.@spawn parts[c] = Base.invokelatest(gaussbatchadjoint, ws[c], sol[ch], wh, single, sink, stagesink,
                    outputview(currents, columns(ch)))
            end
        catch err
            throw(chunkerror(err))
        end
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
pair for every condition, or one per condition with the conditions as a
trailing dimension of the pair's arrays. On the host the conditions are
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
# increments with the states as the columns of `x`, the norm per column
function gaussbatchresidual!(norms, residual, sys::TransientSystem, gc::GaussCoefficients, delta, x, lx, X,
        phi, junction, jwork, cwork, gwork, rhs, colnorm, roundoff, colfloor,
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
    # the norm per column over both stages
    colnorm .= max.(dropdims(maximum(abs, stage(residual, 1); dims = 1); dims = 1),
        dropdims(maximum(abs, stage(residual, 2); dims = 1); dims = 1))
    copyto!(norms, colnorm)
    gaussresidualroundoff!(roundoff, sys, gc, delta, lx, junction, rhs, colfloor, rw)
    return nothing
end

# A floor for each condition at the point whose residual was just formed.
# The matrix norms cover cancellation within a product; the other terms
# are already in residual units. The blocks use |S| |y| so unrelated
# state coordinates cannot inflate one another's bound. Base and trial
# evaluations keep separate floors, just as they keep separate norms.
function gaussresidualroundoff!(roundoff, sys::TransientSystem, gc::GaussCoefficients,
        delta, lx, junction, rhs, colfloor, @nospecialize(rw::Union{Nothing, AbstractRationalWork}))
    cs, gs, ls, _ = sys.rowsums
    h = sys.h
    rounding = 2maximum(abs, gc.ainv2)*cs/h^2 + 2maximum(abs, gc.ainv)*gs/h + ls
    isnothing(rw) || reflectedwavesbound!(rw, sys.gauss.coupling)
    blockbound = isnothing(rw) ? nothing : rw.sbound
    blockweight = isnothing(rw) ? 0.0 : sys.gauss.coupling.nstates*eps(Float64)
    gaussroundoffcolumns!(roundoff, rounding, delta, lx, junction, rhs,
        blockbound, blockweight, colfloor, sys.backend)
    return nothing
end

function gaussroundoffcolumns!(roundoff, rounding, delta, lx, junction, rhs,
        blockbound, blockweight, colfloor, backend)
    colmax = a -> vec(maximum(abs, a; dims = (1, 3)))
    colfloor .= 8eps(Float64) .* (rounding .* colmax(delta) .+
        vec(maximum(abs, lx; dims = 1)) .+ colmax(junction) .+ colmax(rhs))
    isnothing(blockbound) || (colfloor .+= blockweight .* colmax(blockbound))
    copyto!(roundoff, colfloor)
    return nothing
end

# A host reduction writes the column bounds directly without temporary
# reduction arrays in every residual evaluation. Devices keep the batched
# reductions above and copy only the bounds to the Newton engine.
function gaussroundoffcolumns!(roundoff, rounding, delta, lx, junction, rhs,
        blockbound, blockweight, colfloor, ::CPU)
    @inbounds for col in axes(delta, 2)
        dm, lm, jm, rm, bm = 0.0, 0.0, 0.0, 0.0, 0.0
        for row in axes(delta, 1)
            lm = max(lm, abs(lx[row, col]))
            for i in 1:2
                dm = max(dm, abs(delta[row, col, i]))
                jm = max(jm, abs(junction[row, col, i]))
                rm = max(rm, abs(rhs[row, col, i]))
                isnothing(blockbound) || (bm = max(bm, abs(blockbound[row, col, i])))
            end
        end
        roundoff[col] = 8eps(Float64)*(rounding*dm + lm + jm + rm) + blockweight*bm
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

"""
    refreshstageoperator!(rw::RationalWork, sys::TransientSystem,
        bf::GaussBatchFactor)

Refresh the entries of the rational blocks in the stage operator of the
factor `bf` at the step's weights, the mean of the two stages' (see
[`rationalvalues!`](@ref)), for a circuit with a pumped block; nothing
otherwise. The values are the factor's own, so the system is read only
to the batches stepping on it at once. The caller factorizes the
operator again after it, and rebuilds the stage correction on the new
factorization.
"""
function refreshstageoperator!(rw::RationalWork, sys::TransientSystem, bf::GaussBatchFactor)
    g = sys.gauss
    g.pumped || return false
    weights = [(rw.weights[j, 1] + rw.weights[j, 2])/2 for j in axes(rw.weights, 1)]
    rationalvalues!(bf.hostvals, sys.problem, g.coefficients.mu/sys.h, sys.Lscale,
        g.stageterms.pattern, g.transposedvals, g.stageterms.terms, weights)
    copyto!(bf.rationalvals, bf.hostvals)
    return true
end

# the correction of the current step: `V'` from the stage weights, then
# `K` and `M` on the current factorization
function stagecorrection!(rw::RationalWork, sys::TransientSystem, bf::GaussBatchFactor, gc::GaussCoefficients, rc, zc, transposed::Bool)
    g = sys.gauss
    g.pumped || return nothing
    cp = g.coupling
    n = size(cp.Phost, 2) ÷ 2
    nz = cp.nstates
    # the correction acts on the port rows with a modulated output alone
    ports = cp.modulated
    nports = length(ports)
    r = 2nports
    m = size(rw.z, 2)
    backend = sys.backend
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    if isnothing(rw.correction)
        St = devicesparse(sparse(transpose(cp.Shost[:, ports])), backend)
        Ucols = allocate(n, nports)
        copyto!(Ucols, Matrix(cp.Shost[:, ports]))
        rw.correction = StageCorrection(typeof(St)[], allocate(2n, r), St, Ucols, allocate(n, m, 2, r),
            LU{Float64,Matrix{Float64},Vector{Int}}[], transposed, allocate(n, m, 2), allocate(n, m, 2),
            allocate(nports, m), zeros(nports, m), zeros(r, m), allocate(m, r), zeros(m, r))
    end
    sc = rw.correction
    sc.transposed = transposed
    # V' per stage: the difference of the weights times the output terms
    # on the stacked states from the stage unknowns
    Vrows = zeros(2n, r)
    empty!(sc.Vt)
    for i in 1:2
        mean = (rw.weights[:, 1] .+ rw.weights[:, 2]) ./ 2
        delta = rw.weights[:, i] .- mean
        Pi = cp.Phost[(i - 1)*nz + 1:i*nz, :]
        V = spzeros(nports, 2n)
        for j in eachindex(cp.Cblk)
            iszero(delta[j]) && continue
            V = V - delta[j] .* (cp.Cblk[j]*Pi)[ports, :]
        end
        push!(sc.Vt, devicesparse(sparse(V), backend))
        Vrows[:, (i - 1)*nports + 1:i*nports] .= Matrix(transpose(V))
    end
    copyto!(sc.Vrows, Vrows)
    # the columns of U, or of V for the transposed operator, solved on
    # the factorization, each replicated over the columns of the batch
    Mh = zeros(r, r, m)
    for col in 1:r
        i, q = (col - 1) ÷ nports + 1, (col - 1) % nports + 1
        fill!(sc.rhs, 0)
        if transposed
            view(sc.rhs, :, :, 1) .= view(sc.Vrows, 1:n, col)
            view(sc.rhs, :, :, 2) .= view(sc.Vrows, n + 1:2n, col)
        else
            view(sc.rhs, :, :, i) .= view(sc.Ucols, :, q)
        end
        gaussbatchtransform!(sc.sol, sc.rhs, bf, gc, rc, zc, transposed, sys)
        view(sc.K, :, :, :, col) .= sc.sol
        # the rows of M from this column: V' K, or U' K transposed
        Mh[:, col, :] .= stagerows!(sc, rw, sc.sol, nports)
    end
    empty!(sc.M)
    for c in 1:m
        Mc = Matrix{Float64}(I, r, r) .+ view(Mh, :, :, c)
        push!(sc.M, lu(Mc))
    end
    return nothing
end

# the products on the block's rows at both stages of a stage pair `c`,
# `V' c`, or `U' c` for the transposed operator, into the host copy
# `sc.yh` as `(r, m)` stage major
function stagerows!(sc::StageCorrection, rw::RationalWork, c, nports)
    sc.transposed || stack!(rw.dstack, c)
    for i in 1:2
        if sc.transposed
            stepmul!(sc.y, sc.St, stage(c, i))
        else
            stepmul!(sc.y, sc.Vt[i], rw.dstack)
        end
        copyto!(sc.ys, sc.y)
        view(sc.yh, (i - 1)*nports + 1:i*nports, :) .= sc.ys
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
    nports = length(cp.modulated)
    r = 2nports
    m = size(sc.zh, 1)
    yh = stagerows!(sc, rw, c, nports)
    for cc in 1:m
        sc.zh[cc, :] .= sc.M[cc] \ view(yh, :, cc)
    end
    copyto!(sc.z, sc.zh)
    for col in 1:r
        c .-= view(sc.K, :, :, :, col) .* reshape(view(sc.z, :, col), 1, m, 1)
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
    nz = cp.nstates
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
# the columns of `b`
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
# the same as the column before the first read and the weights, for
# the tests of the stencil
function linestencil(s, tpre, h, accepted)
    weights = zeros(6)
    first, nst = linestencil!(weights, s, tpre, h, accepted)
    return first, weights[1:nst]
end

# the read of every port's far wave through its stencil, scaled, over
# the ports and the conditions, and its transpose scattering a value
# onto the samples it read, with a sign, accumulating; the columns in
# the slots of the ring
@kernel function wavegatherkernel!(values, @Const(waves), @Const(stencil), @Const(far), @Const(scale), nring)
    row, j = @index(Global, NTuple)
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
end
function wavegather!(values, waves, stencil, far, scale, backend)
    isempty(values) && return values
    wavegatherkernel!(backend, (64, 1))(values, waves, stencil, far, scale, size(waves, 2); ndrange = size(values))
    return values
end
# on the host a loop: a kernel launch there costs more than the read
function wavegather!(values, waves, stencil, far, scale, ::CPU)
    nring = size(waves, 2)
    @inbounds for j in axes(values, 2), row in axes(values, 1)
        f = far[row]
        first = Int(stencil[1, row])
        nst = Int(stencil[2, row])
        total = 0.0
        for k in 1:nst
            total += stencil[2 + k, row]*waves[f, mod1(first + k, nring), j]
        end
        values[row, j] = scale[row]*total
    end
    return values
end
@kernel function wavescatterkernel!(wavesbar, @Const(values), @Const(stencil), @Const(far), @Const(scale), sign, nring)
    row, j = @index(Global, NTuple)
    @inbounds begin
        f = far[row]
        first = Int(stencil[1, row])
        nst = Int(stencil[2, row])
        v = sign*scale[row]*values[row, j]
        for k in 1:nst
            wavesbar[f, mod1(first + k, nring), j] += stencil[2 + k, row]*v
        end
    end
end
function wavescatter!(wavesbar, values, stencil, far, scale, sign, backend)
    isempty(values) && return wavesbar
    wavescatterkernel!(backend, (64, 1))(wavesbar, values, stencil, far, scale, Float64(sign), size(wavesbar, 2); ndrange = size(values))
    return wavesbar
end
function wavescatter!(wavesbar, values, stencil, far, scale, sign, ::CPU)
    nring = size(wavesbar, 2)
    @inbounds for j in axes(values, 2), row in axes(values, 1)
        f = far[row]
        first = Int(stencil[1, row])
        nst = Int(stencil[2, row])
        v = sign*scale[row]*values[row, j]
        for k in 1:nst
            wavesbar[f, mod1(first + k, nring), j] += stencil[2 + k, row]*v
        end
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
    bend::M
    phi::A
    jwork::A
    # the derivative of the relation at one stage, for a circuit which has
    # a polynomial one, and empty for the Josephson relation, which is
    # broadcast in place
    dwork::M
    cosphi::C
    rc::C
    zc::C
    colnorm::W
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
    tolerance::Vector{Float64}
    roundoff::Tuple{Vector{Float64}, Vector{Float64}}
    rtol::Float64
    atol::Float64
    iterations::Int
    stalefailed::Bool
    corrections::Int
    factorizations::Int
    retries::Int
end

function gaussstepper(sys::TransientSystem, problems, rtol, atol, iterations, bf)
    p = sys.problem
    backend = sys.backend
    gc = sys.gauss.coefficients
    N = length(problems)
    n, np, nd = length(p), length(p.portimpedances), length(p.drives)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    x, v = allocate(n, N), allocate(n, N)
    X, delta, lastdelta, residual, trial, trialresidual, correction, rhs = [allocate(n, N, 2) for _ in 1:8]
    junction, trialjunction, cwork, gwork = [allocate(n, N, 2) for _ in 1:4]
    xnew, cv, lx, b1, b2, bend = [allocate(n, N) for _ in 1:6]
    nj = length(sys.lmolj)
    phi, trialphi, jwork = [allocate(nj, N, 2) for _ in 1:3]
    dwork = allsinusoidal(sys.relations) ? allocate(nj, 0) : allocate(nj, N)
    cosphi = KernelAbstractions.zeros(backend, ComplexF64, nj, N)
    rc, zc = [KernelAbstractions.zeros(backend, ComplexF64, n, N) for _ in 1:2]
    colnorm = allocate(N)
    colfloor = allocate(N)
    roundoff = (zeros(N), zeros(N))
    hostvalues = zeros(nd, N)
    values = tobackend(backend, zeros(nd, N))
    portwork = allocate(np, N)
    rw = isempty(sys.gauss.rational) ? nothing : rationalwork(p, backend, n, N)
    cell = RationalWorkCell{typeof(x), typeof(X)}(rw)
    baseresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, phi, junction, jwork, cwork, gwork, rhs, colnorm, roundoff[1], colfloor, cell.work)
    trialresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, trialphi, trialjunction, jwork, cwork, gwork, rhs, colnorm, roundoff[2], colfloor, cell.work)
    # a new factorization of the frozen operator, and the stage
    # correction of a pumped block rebuilt on it, whichever asked
    refresh! = () -> (gaussbatchjacobian!(bf, sys, phi, cosphi, dwork);
        isnothing(cell.work) || stagecorrection!(cell.work, sys, bf, gc, rc, zc, false); nothing)
    solve! = (c, r) -> (gaussbatchtransform!(c, r, bf, gc, rc, zc); isnothing(cell.work) || stagecorrect!(c, cell.work, sys); (false, 0))
    accept! = mask -> (maskcolumns!(phi, trialphi, mask); maskcolumns!(junction, trialjunction, mask); nothing)
    pr = sys.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, N)
    nl = 2length(p.lines)
    far, readscale, sqrtz = linetables(p, backend)
    npre = lineprehistory(p, sys.h)
    return GaussStepper(sys, problems, N, x, v, X, delta, lastdelta, residual, trial, trialresidual, correction, rhs,
        junction, cwork, gwork, xnew, cv, lx, b1, b2, bend, phi, jwork, dwork, cosphi, rc, zc,
        colnorm, hostvalues, values, portwork, bf, baseresidual!, trialresidual!, refresh!, solve!, accept!, pw,
        allocate(nl, linering(npre), N), npre, Float64[], 0, allocate(nl, N), zeros(8, nl), allocate(8, nl), far, readscale, sqrtz,
        allocate(nl, N), allocate(n, N), rw,
        NewtonWork(backend, N), zeros(N), roundoff, Float64(rtol), Float64(atol), Int(iterations), false, 0, 0, 0)
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
function gaussproject!(st::GaussStepper, t)
    sys, pw = st.sys, st.pw
    pr = sys.projection
    gc = sys.gauss.coefficients
    h = sys.h
    N = st.N
    drivevalues!(st.hostvalues, st.problems, t)
    isempty(sys.problem.lines) || linevalues!(st, t)
    linedrive = isempty(sys.problem.lines) ? nothing : Array(st.linevalues)
    if !isempty(pr.directions)
        # the drive along the constraints' rows, and the magnitudes of its
        # terms for the residual's floor
        gb = pr.Ztinj*st.hostvalues .+ pr.Ztconstant
        gbabs = pr.Ztinjabs*abs.(st.hostvalues) .+ abs.(pr.Ztconstant)
        if !isnothing(linedrive)
            gb .+= pr.Ztline*linedrive
            gbabs .+= pr.Ztlineabs*abs.(linedrive)
        end
        if !isnothing(st.rw)
            waves = restingwaves(sys, st.rw.z, t)
            gb .+= pr.Ztblock*waves
            gbabs .+= pr.Ztblockabs*abs.(waves)
        end
        # the predictor moves with the state
        moved! = w -> (stage(st.lastdelta, 1) .-= w; stage(st.lastdelta, 2) .-= w; nothing)
        projectendpoint!(st.xnew, pw, pr, gb, gbabs, st.iterations, moved!) ||
            error(lazy"the projection onto the algebraic constraints at t = $(t) s did not converge; reduce dt or check the state.")
        pr.cubic && projectrate!(st.v, pr, pw, gc, h, stage(st.delta, 1), stage(st.delta, 2), st.xnew, st.x)
    end
    # the index one unknowns of the endpoint from their equations
    if !isempty(pr.readrows) || !isempty(pr.auxrows)
        projectedphases!(st, st.xnew)
        copyto!(pw.hphi, pw.phip)
        gq = pr.Qinj*st.hostvalues .+ pr.Qconst
        isnothing(linedrive) || (gq .+= pr.Qline*linedrive)
        isnothing(st.rw) || (gq .+= pr.Qblock*restingwaves(sys, st.rw.z, t))
        endpointread!(st.v, st.xnew, pr, pw,
            pr.lmoljp .* relationat(pr.relationsp, pw.hphi), gq)
        projectedphases!(st, st.xnew)
    end
    return st
end

# the stepper set to a state and the stage increments to start from, its
# stage matrices assembled at that state
function setstate!(st::GaussStepper, x, v, lastdelta)
    copyto!(st.x, x)
    copyto!(st.v, v)
    isnothing(lastdelta) ? fill!(st.lastdelta, 0) : copyto!(st.lastdelta, lastdelta)
    st.X .= st.x
    stepmul!(stage(st.phi, 1), st.sys.RJ, stage(st.X, 1)); stepmul!(stage(st.phi, 2), st.sys.RJ, stage(st.X, 2))
    gaussbatchjacobian!(st.bf, st.sys, st.phi, st.cosphi, st.dwork)
    st.factorizations += 1
    st.stalefailed = false
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
        # a pumped block's converted coupling in the operator, at this
        # step's weights
        if refreshstageoperator!(st.rw, sys, st.bf)
            st.refresh!()
            st.factorizations += 1
        end
    end
    batchdrivecurrent!(st.b1, st, tprev + gc.c[1]*h)
    batchdrivecurrent!(st.b2, st, tprev + gc.c[2]*h)
    stepmul!(st.cv, sys.C, st.v)
    stage(st.rhs, 1) .= st.b1 .+ (gc.ainvone[1]/h) .* st.cv
    stage(st.rhs, 2) .= st.b2 .+ (gc.ainvone[2]/h) .* st.cv
    stepmul!(st.lx, sys.L, st.x)
    st.colnorm .= max.(dropdims(maximum(abs, stage(st.rhs, 1); dims = 1); dims = 1),
        dropdims(maximum(abs, stage(st.rhs, 2); dims = 1); dims = 1),
        dropdims(maximum(abs, st.lx; dims = 1); dims = 1), 1.0)
    copyto!(st.tolerance, st.colnorm)
    st.tolerance .= st.atol .+ st.rtol .* st.tolerance
    copyto!(st.delta, st.lastdelta)
    converged, fresh, ncorr, nfact, nretry, _, _ = newtonsolve!(st.delta, st.correction, st.trial,
        st.residual, st.trialresidual, st.baseresidual!, st.trialresidual!, st.refresh!, st.solve!,
        st.tolerance, st.iterations, st.stalefailed, length(sys.lmolj) > 0, false, st.newtonwork;
        simplified = true, accept! = st.accept!, roundoff = st.roundoff)
    st.corrections += ncorr
    st.factorizations += nfact
    st.retries += nretry
    converged || error(lazy"the Newton solve of step $(step) at t = $(t) s did not converge on every condition; reduce dt or check the initial states and the circuit.")
    st.stalefailed = fresh
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
    isnothing(st.pw) || gaussproject!(st, t)
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
        zeroed(nl, nl > 0 && (!savecheckpoints || !tails) ? nsteps + 1 : 0, N),
        zeroed(nl, nl > 0 ? npre : 0, N),
        KernelAbstractions.allocate(backend, Float64, n, kept(savestates), N),
        KernelAbstractions.allocate(backend, Float64, n, kept(savestates), N),
        zeroed(n, 2, kept(savestates), N),
        zeroed(nzs, kept(savestates && nzs > 0), N),
        checkpoints, zeroed(n, N), zeroed(n, N), zeroed(n, N), zeroed(n, N), zeros(nl, N), zeros(nzs, N))
end

# the arrays of the conditions `ch`, views of the batch's along the
# condition axis, which a chunk fills as though they were its own
function chunkoutputs(out::GaussBatchOutputs, ch)
    slice = a -> selectdim(a, ndims(a), ch)
    cp = out.checkpoints
    checkpoints = GaussCheckpoints(cp.every, slice(cp.flux), slice(cp.rate), slice(cp.increment), slice(cp.states),
        slice(cp.waves))
    return GaussBatchOutputs(slice(out.voltage), slice(out.incident), slice(out.outgoing),
        slice(out.phases), slice(out.endphases), slice(out.endrates), slice(out.linewaves), slice(out.history),
        slice(out.flux), slice(out.rate), slice(out.stages), slice(out.blockstates), checkpoints,
        slice(out.initialflux), slice(out.initialrate), slice(out.finalflux), slice(out.finalrate),
        slice(out.initialwaves), slice(out.initialstates))
end

# the checkpoints of a batch's solution: the record's named tuple, or
# nothing without any
function recordedcheckpoints(cp::GaussCheckpoints)
    cp.every > 0 || return nothing
    return (; every = cp.every, flux = cp.flux, rate = cp.rate, increment = cp.increment, states = cp.states,
        waves = cp.waves)
end

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
# grid, and within a chunk the counters are those of the condition which
# converged worst, so across chunks they are those of the worst chunk;
# summing them would count one grid several times over.
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
    stats = if length(chunks) == 1
        ch = only(chunks)
        gaussbatchrun!(chunkoutputs(out, ch), sys, problems[ch], chunkinitial(init, ch), times, saveevery, rtol, atol,
            iterations, factors[1])
    else
        each = Vector{Any}(undef, length(chunks))
        # the chunks share the system, which they read and do not write,
        # and nothing else: their steppers, factors and arrays are their own
        try
            Base.Threads.@sync for (k, ch) in enumerate(chunks)
                Base.Threads.@spawn each[k] = gaussbatchrun!(chunkoutputs(out, ch), sys, problems[ch],
                    chunkinitial(init, ch), times, saveevery, rtol, atol, iterations, factors[k])
            end
        catch err
            throw(chunkerror(err))
        end
        mergebatchstats(each)
    end
    isnothing(reuse) || (reuse.factor = factors)
    return TransientBatchSolution(problems, sys.method, h, savedtimes, out.voltage, out.incident, out.outgoing,
        recordedarray(out.phases), recordedarray(out.endphases), recordedarray(out.endrates), recordedarray(out.linewaves),
        recordedarray(out.history), recordedarray(out.flux), recordedarray(out.rate), recordedarray(out.stages),
        recordedcheckpoints(out.checkpoints), out.initialflux, out.initialrate, out.finalflux, out.finalrate,
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
    withblocks = !isempty(sys.gauss.rational)
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
        violation, _ = transientconsistency(sys, view(x0, :, j), view(v0, :, j), t0, problems[j],
            q0[:, j], isnothing(resting) ? nothing : resting[:, j], qdot[:, j])
        violation <= atol + rtol || throw(ArgumentError(
            lazy"the initial state of condition $(j) violates the algebraic equations of the circuit along a direction without capacitance to ground (a node no capacitor touches, a capacitive island, a coupled inductor or gauge row); supply a consistent transientstate, or start the drive from an equilibrium."))
    end
    return (; x0, v0, w0, tailh, z0)
end

# the initial arrays of the conditions `ch` of a batch, copied so that a
# chunk moves plain host arrays to its backend
chunkinitial(init, ch) = (; x0 = init.x0[:, ch], v0 = init.v0[:, ch], w0 = init.w0[:, ch],
    tailh = init.tailh[:, :, ch], z0 = init.z0[:, ch])

# the error a chunk's task met, out of the failure of the task the
# threaded batch reports it as
chunkerror(err::CompositeException) = isempty(err.exceptions) ? err : chunkerror(first(err.exceptions))
chunkerror(err::TaskFailedException) = chunkerror(err.task.result)
chunkerror(err) = err

# One chunk of a batch stepped along `times` into the arrays `out`, on
# the factorizations `bf`, from the initial arrays `init` of its
# conditions: every step advanced and saved. The record is whichever
# arrays of `out` are there to be filled. Returns the chunk's statistics.
function gaussbatchrun!(out, sys::TransientSystem, problems, init, times, saveevery,
        rtol, atol, iterations, bf::GaussBatchFactor)
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
    st = gaussstepper(sys, problems, rtol, atol, iterations, bf)
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
    setstate!(st, init.x0, init.v0, nothing)
    # the history of the line waves before the start; the record of the
    # waves leaving every port at every step, or with checkpoints the
    # history before each of them
    tail = tobackend(backend, init.tailh)
    sethistory!(st, times, tail, 0)
    savehistory && copyto!(out.history, tail)
    savelinewaves && (view(out.linewaves, :, 1, :) .= tobackend(backend, init.w0))
    saveendphases && projectedphases!(st, st.x)
    batchdrivecurrent!(st.bend, st, t0)
    saveoutputs!(1, t0)
    copyto!(out.initialflux, st.x)
    copyto!(out.initialrate, st.v)
    saved = 1
    tprev = t0
    for step in 1:nsteps
        t = times[step + 1]
        if savecheckpoints && (step - 1) % K == 0
            c = (step - 1) ÷ K + 1
            copyto!(view(checkpoints.flux, :, c, :), st.x)
            copyto!(view(checkpoints.rate, :, c, :), st.v)
            isnothing(st.rw) || (view(checkpoints.states, :, c, :) .= st.rw.z)
            view(checkpoints.increment, :, 1, c, :) .= stage(st.lastdelta, 1)
            view(checkpoints.increment, :, 2, c, :) .= stage(st.lastdelta, 2)
            nl > 0 && tails && readhistory!(view(checkpoints.waves, :, :, c, :), st.waves, step)
        end
        advance!(st, tprev, t, step)
        tprev = t
        savelinewaves && (view(out.linewaves, :, step + 1, :) .=
            view(st.waves, :, ringslot(npre + step, size(st.waves, 2)), :))
        if step % saveevery == 0 || step == nsteps
            saved += 1
            batchdrivecurrent!(st.bend, st, t)
            saveoutputs!(saved, t)
        end
    end
    KernelAbstractions.synchronize(backend)
    copyto!(out.finalflux, st.x)
    readout!(out.finalrate, tprev)
    return (; steps = nsteps, newtoncorrections = st.corrections, factorizations = st.factorizations,
        retries = st.retries, kryloviterations = 0, rtol = Float64(rtol), atol = Float64(atol), iterations = Int(iterations))
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
            # the replayed window ends where the next checkpoint was taken
            if c < nc
                next = view(cps.flux, :, c + 1, :)
                gap = maximum(abs, st.x .- next)
                gap <= 1e3*(sol.stats.atol + sol.stats.rtol*maximum(abs, next)) || error(
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
    cp = sol.checkpoints
    cps = isnothing(cp) ? nothing : (; every = cp.every, flux = expand(cp.flux), rate = expand(cp.rate),
        increment = expand(cp.increment), states = expand(cp.states), waves = expand(cp.waves))
    return TransientBatchSolution([sol.problem], sol.method, sol.dt, sol.times, expand(sol.voltage), expand(sol.incident),
        expand(sol.outgoing), expand(sol.phases), expand(sol.endphases), expand(sol.endrates), expand(sol.linewaves),
        expand(sol.history), expand(sol.flux), expand(sol.rate), expand(sol.stages), cps, expand(sol.initialflux),
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
    cp = isnothing(b.checkpoints) ? nothing : (; every = b.checkpoints.every, flux = squeeze(b.checkpoints.flux),
        rate = squeeze(b.checkpoints.rate), increment = squeeze(b.checkpoints.increment), states = squeeze(b.checkpoints.states),
        waves = squeeze(b.checkpoints.waves))
    return TransientSolution(b.problems[1], b.method, b.dt, b.times, squeeze(b.voltage), squeeze(b.incident),
        squeeze(b.outgoing), squeeze(b.phases), squeeze(b.endphases), squeeze(b.endrates), squeeze(b.linewaves), squeeze(b.history), squeeze(b.flux), squeeze(b.rate), squeeze(b.stages), cp, vec(b.initialflux),
        vec(b.initialrate), vec(b.finalflux), vec(b.finalrate), squeeze(b.blockstates), squeeze(b.initialwaves), squeeze(b.initialstates), b.stats)
end
