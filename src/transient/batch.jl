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
    return TransientProblem(psc, p.graph, p.matrices, p.Nnodal, p.Naux, p.Lscale, p.coupledbranches,
        p.floatingcomponents, p.gaugeindices, p.inertialess, p.algebraic, p.directions, p.constraints, p.rateextraction,
        injection, drives, constantcurrent, p.portpositive, p.portnegative, p.portimpedances, p.portconductances, p.blocks, p.lines)
end

# The problems of a batch share everything but their waveforms: the
# compiled circuit, the injection of the drives and the ports they drive,
# and the constant current of the sources no drive replaced, which is
# what the binding of the sources amounts to.
function batchcompatible(problems)
    isempty(problems) && throw(ArgumentError("a batch needs at least one problem."))
    p = first(problems)
    for q in problems
        (q.circuit === p.circuit && q.matrices === p.matrices) || throw(ArgumentError(
            "the problems of a batch must share one compiled circuit: build them with transientproblem(problem; sources)."))
        (q.injection == p.injection && length(q.drives) == length(p.drives) &&
            all(k -> q.drives[k].portindex == p.drives[k].portindex, eachindex(p.drives)) &&
            q.constantcurrent == p.constantcurrent) || throw(ArgumentError(
            "the problems of a batch must drive the same targets in the same order and leave the same sources constant; only the waveforms may differ."))
    end
    return p
end

"""
    TransientBatchSolution

Result of [`transientsolve`](@ref) on a vector of problems: the arrays of
a [`TransientSolution`](@ref) with the conditions as the trailing
dimension, `voltage[port, time, condition]`, `flux[state, time, condition]`,
`finalflux[state, condition]`, `phases[junction, stage, time, condition]`; `problems` holds the batch. Indexing,
`solution[j]`, is the ordinary solution of condition `j`, a view of the
batch's arrays, on which the demodulation, the tangent, the adjoint and
the noise run as on any solution.
"""
struct TransientBatchSolution{P, M, V}
    problems::P
    method::AbstractTransientIntegrator
    dt::Float64
    times::Vector{Float64}
    voltage::M
    incident::M
    outgoing::M
    phases::Any
    endphases::Any
    linewaves::Any
    # the flux and rate records untyped, as the phases are, so that the
    # record level does not make a new solution type and the responses
    # compile once
    flux::Any
    rate::Any
    checkpoints::Any
    initialflux::V
    initialrate::V
    finalflux::V
    finalrate::V
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
        slice(b.outgoing), slice(b.phases), slice(b.endphases), slice(b.linewaves), slice(b.flux), slice(b.rate), cp, slice(b.initialflux),
        slice(b.initialrate), slice(b.finalflux), slice(b.finalrate), slice(b.blockstates), slice(b.initialwaves), slice(b.initialstates), b.stats)
end
function Base.getindex(b::TransientBatchSolution, j::Integer)
    1 <= j <= length(b) || throw(BoundsError(b, j))
    slice = a -> isnothing(a) ? nothing : selectdim(a, ndims(a), j)
    cp = isnothing(b.checkpoints) ? nothing : (; every = b.checkpoints.every, flux = slice(b.checkpoints.flux),
        rate = slice(b.checkpoints.rate), increment = slice(b.checkpoints.increment), states = slice(b.checkpoints.states),
        waves = slice(b.checkpoints.waves))
    return TransientSolution(b.problems[j], b.method, b.dt, b.times, slice(b.voltage), slice(b.incident),
        slice(b.outgoing), slice(b.phases), slice(b.endphases), slice(b.linewaves), slice(b.flux), slice(b.rate), cp, copy(view(b.initialflux, :, j)),
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
struct GaussBatchFactor{J, F, X, B, I, IT}
    ncolumns::Int
    jacobians::J
    factors::F
    X::X
    B::B
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
        return GaussBatchFactor(ncolumns, jacobians, Vector{Any}(nothing, ncolumns), nothing, nothing, UnitRange{Int}[],
            nothing, nothing, sys.symmetric, Any[], nothing, nothing, nothing, nothing, Ref(true))
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
    return GaussBatchFactor(ncolumns, nzval, Vector{Any}(nothing, length(chunks)), X, B, chunks, rowptr, colind,
        sys.symmetric, Vector{Any}(nothing, length(chunks)), trowptr, tcolind, tperm, tnzval, Ref(true))
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
function gaussbatchjacobian!(bf::GaussBatchFactor, sys::TransientSystem, phi, cosphi)
    g = sys.gauss
    cosphi .= (cos.(stage(phi, 1)) .+ cos.(stage(phi, 2))) ./ 2
    ncolumns = size(cosphi, 2)
    if sys.backend isa CPU
        for j in 1:ncolumns
            A = bf.jacobians[j]
            assemblerealjacobian!(nonzeros(sys.jacobian), sys.plan, view(cosphi, :, j))
            nonzeros(A) .= nonzeros(sys.jacobian) .+ im .* g.imvals .+ g.rationalvals
            bf.factors[j] = isnothing(bf.factors[j]) ? factorize(sys.factorization, A) :
                refactorize!(sys.factorization, bf.factors[j], A)
        end
        return bf
    end
    nzval = bf.jacobians
    # every condition's assembly launched, then one synchronization
    for j in 1:ncolumns
        assemblerealjacobian!(nonzeros(sys.jacobian), sys.plan, view(cosphi, :, j); synchronize = false)
        view(nzval, :, j) .= nonzeros(sys.jacobian) .+ im .* g.imvals .+ g.rationalvals
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
function batchjunctionproduct!(y, sys::TransientSystem, phi, d, work)
    nj, N = size(phi)
    stepmul!(work, sys.RJ, d)
    reshape(work, nj, :, N) .*= sys.lmolj .* reshape(cos.(phi), nj, 1, N)
    stepmul!(y, sys.RJt, work)
    return y
end

# The exact solve of the stage equations of every condition linearized at
# its recorded stage phases, `J d = r` with the true stage Jacobian whose
# two stiffness blocks differ, or its transpose: a fixed point iteration
# on the frozen complex operators, `d += [J*]^{-1} (r - J d)`, contracting
# as the simplified Newton of the step does. `d` and `r` are `(n, ndir*N, 2)`
# with the directions of a condition contiguous; `phi` is `(nj, N, 2)`.
function gaussbatchstagesolve!(d, r, sys::TransientSystem, gc::GaussCoefficients, bf::GaussBatchFactor, phi,
        transposed::Bool, res, dd, cwork, gwork, lwork, jwork, rc, zc; rtol = 1e-12, maxiters = 100, rw = nothing)
    h = sys.h
    fill!(d, 0)
    # the convergence per column: each direction of each condition against
    # its own right hand side, so a weak direction is not left at the
    # tolerance of a strong one, down to the roundoff floor of the terms
    # the residual sums, the unit roundoff times their magnitudes
    colmax = a -> vec(Array(maximum(abs, a; dims = (1, 3))))
    scale = max.(colmax(r), floatmin(Float64))
    cmax, gmax = maximum(abs, gc.ainv2)/h^2, maximum(abs, gc.ainv)/h
    # the contraction of the first correction, the measure of how well the
    # frozen operator fits this step's two stiffnesses, returned for the
    # caller to decide whether to refresh it
    first = copy(scale)
    contraction = 0.0
    for iteration in 1:maxiters
        for i in 1:2
            di = stage(d, i)
            stepmul!(stage(cwork, i), sys.C, di)
            stepmul!(stage(gwork, i), transposed ? sys.Gt : sys.G, di)
            stepmul!(stage(lwork, i), transposed ? sys.Lt : sys.L, di)
            batchjunctionproduct!(stage(res, i), sys, view(phi, :, :, i), di, jwork)
        end
        # the rational blocks' coupling of the stages, or its transpose
        if !isnothing(rw)
            transposed ? rationalsourcestranspose!(rw, sys, d) : rationalsources!(rw, sys, d, d; withstate = false)
            res .-= rw.source
        end
        floor = 8eps(Float64) .* (scale .+ 2cmax .* colmax(cwork) .+ 2gmax .* colmax(gwork) .+ colmax(lwork) .+ colmax(res))
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
        all(current .<= max.(rtol .* scale, floor)) && return contraction
        gaussbatchtransform!(dd, res, bf, gc, rc, zc, transposed, sys)
        d .+= dd
    end
    error("the stage equations of a Gauss-Legendre step did not converge in the tangent or the adjoint: the two stage stiffnesses differ too much from their mean; reduce dt.")
end

# a solution as a batch of one for the responses: its problems, its
# phases as `(nj, 2, nt, N)`, its initial states as `(n, N)`, and its
# endpoint phases as `(npj, nt, N)`
function batchview(sol::TransientSolution)
    phases = isnothing(sol.phases) || ndims(sol.phases) != 3 ? sol.phases : reshape(sol.phases, size(sol.phases)..., 1)
    endphases = isnothing(sol.endphases) ? nothing : reshape(sol.endphases, size(sol.endphases)..., 1)
    linewaves = isnothing(sol.linewaves) ? nothing : reshape(sol.linewaves, size(sol.linewaves)..., 1)
    w0 = isnothing(sol.initialwaves) ? zeros(2length(sol.problem.lines), 1) : reshape(sol.initialwaves, :, 1)
    z0 = isnothing(sol.initialstates) ? zeros(blockstates(sol.problem), 1) : reshape(sol.initialstates, :, 1)
    return [sol.problem], phases, reshape(sol.initialflux, :, 1), reshape(sol.initialrate, :, 1), endphases, linewaves, w0, z0
end
batchview(b::TransientBatchSolution) = b.problems, b.phases, b.initialflux, b.initialrate, b.endphases, b.linewaves, b.initialwaves, b.initialstates
recordedsolution(b::TransientBatchSolution) = recordedsolution(b[1])

# the tangent of a Gauss-Legendre solve or batch, over `ndir` directions
# for every condition, the currents shared by the conditions and the
# initial perturbation shared or one per condition as a trailing
# dimension; `transienttangent` checked the arguments
function gaussbatchtangent(sol, currents, targets, initialstate, sys::TransientSystem, outputsink = nothing)
    problems, phases, x0s, v0s = batchview(sol)
    N = length(problems)
    p = first(problems)
    backend = sys.backend
    gc = sys.gauss.coefficients
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq, ndir = length(targets), (ndims(currents) == 2 ? 1 : size(currents, ndims(currents)))
    m = ndir*N
    h = sys.h
    injection = devicesparse((sys.Lscale/phi0) .* transientinjection(p, targets), backend)
    # the currents on the grid, read at the stage times through the
    # stencil, or given at the grid and the two stage times of each step
    staged = ndims(currents) == 4
    grid = Float64.(collect(staged ? selectdim(currents, 2, 1) : currents))
    dI = tobackend(backend, reshape(grid, nq, nt, ndir))
    dS = staged ? tobackend(backend, reshape(Float64.(collect(selectdim(currents, 2, 2:3))), nq, 2, nt, ndir)) : nothing
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    dx, dv, dxnew, lx, cdv, work = [allocate(n, m) for _ in 1:6]
    r, d, res, dd, cwork, gwork, lwork = [allocate(n, m, 2) for _ in 1:7]
    # the lines: the ring of the perturbation of the wave leaving each
    # port, with a prehistory long enough for every read before the
    # start, from the initial state's third member or zero, the same
    # history and stencils as the solve's so that the tangent is the
    # derivative of the solve
    lines = p.lines
    nl2 = 2length(lines)
    npre = lineprehistory(p, h)
    nring = linering(npre)
    t0 = first(sol.times)
    tpre = linestart(t0, npre, h)
    dwaves = allocate(nl2, nring, m)
    ddlinevalues, dlinerates, linework = allocate(nl2, m), allocate(nl2, m), allocate(n, m)
    stencil, dstencil = zeros(8, nl2), allocate(8, nl2)
    far, readscale, sqrtz = linetables(p, backend)
    # the currents the arriving waves force at time `t`, read from the
    # history accepted through column `accepted`
    function lineread!(t, accepted)
        linestencil!(stencil, lines, t, tpre, h, accepted)
        copyto!(dstencil, stencil)
        wavegather!(ddlinevalues, dwaves, dstencil, far, readscale, backend)
        nothing
    end
    # the rational blocks: the perturbation of their states, from the
    # initial state's fourth member or zero
    rwt = isempty(sys.gauss.rational) ? nothing : rationalwork(p, backend, n, m)
    zerostages = isnothing(rwt) ? nothing : allocate(n, m, 2)
    fullstages = isnothing(rwt) ? nothing : allocate(n, m, 2)
    if !isnothing(initialstate)
        px, pv = initialstate
        if length(initialstate) >= 4 && !isnothing(rwt)
            pz = initialstate[4]
            size(pz, 1) == size(rwt.z, 1) || throw(DimensionMismatch(lazy"the initial block states need $(size(rwt.z, 1)) rows."))
            for j in 1:N
                pzj = ndims(pz) == 3 ? view(pz, :, :, j) : pz
                copyto!(view(rwt.z, :, (j - 1)*ndir + 1:j*ndir), reshape(Float64.(collect(pzj)), size(rwt.z, 1), :))
            end
        end
        if length(initialstate) >= 3 && nl2 > 0
            pw0 = initialstate[3]
            (size(pw0, 1) == nl2 && size(pw0, 2) == npre) || throw(DimensionMismatch(
                lazy"the initial line waves need $(nl2) rows and $(npre) prehistory columns at this delay and step."))
            for j in 1:N
                pwj = ndims(pw0) == 4 ? view(pw0, :, :, :, j) : pw0
                fillhistory!(view(dwaves, :, :, (j - 1)*ndir + 1:j*ndir), tobackend(backend, reshape(Float64.(collect(pwj)), nl2, npre, :)), 1)
            end
        end
        (size(px, 1) == n && size(pv, 1) == n) || throw(DimensionMismatch(lazy"the state has $(n) entries."))
        percondition = ndims(px) == 3 && size(px, 3) > 1
        !percondition || (size(px, 3) == N && size(pv, 3) == N) || throw(DimensionMismatch(
            lazy"give one initial perturbation for every condition, or one per condition ($(N)) as the trailing dimension."))
        for j in 1:N
            pxj, pvj = percondition ? (view(px, :, :, j), view(pv, :, :, j)) : (px, pv)
            copyto!(view(dx, :, (j - 1)*ndir + 1:j*ndir), reshape(Float64.(collect(pxj)), n, :))
            copyto!(view(dv, :, (j - 1)*ndir + 1:j*ndir), reshape(Float64.(collect(pvj)), n, :))
        end
    end
    nj = length(sys.lmolj)
    jwork, phi = allocate(nj, m), allocate(nj, N, 2)
    cosphi = KernelAbstractions.zeros(backend, ComplexF64, nj, N)
    rc, zc = [KernelAbstractions.zeros(backend, ComplexF64, n, m) for _ in 1:2]
    stagecurrent = allocate(nq, ndir)
    stageinjection, injectionall = allocate(n, ndir), allocate(n, m)
    voltage, incident, outgoing = isnothing(outputsink) ? [allocate(np, nt, m) for _ in 1:3] : (nothing, nothing, nothing)
    coefficients = Dict(q => outputcoefficients(sys, q) for q in (:voltage, :incident, :outgoing))
    portwork = allocate(np, m)
    tp = targetports(p, targets)
    portmap = devicesparse(sparse([q for q in tp if q > 0], [k for (k, q) in enumerate(tp) if q > 0],
        ones(count(>(0), tp)), np, nq), backend)
    directwork, directall = allocate(np, ndir), allocate(np, m)
    # the port outputs of a time, stored, or handed to the sink as the
    # three port by column matrices
    outwork = isnothing(outputsink) ? nothing : [allocate(np, m) for _ in 1:3]
    function outputs!(k)
        stepmul!(portwork, sys.ports, dv)
        portwork .*= phi0
        stepmul!(directwork, portmap, view(dI, :, k, :))
        reshape(directall, np, ndir, N) .= reshape(directwork, np, ndir, 1)
        for (s, q) in enumerate((:voltage, :incident, :outgoing))
            cv, cd = coefficients[q]
            out = isnothing(outputsink) ? view((voltage, incident, outgoing)[s], :, k, :) : outwork[s]
            out .= cv .* portwork .+ cd .* directall
        end
        isnothing(outputsink) || outputsink(k, outwork[1], outwork[2], outwork[3])
        nothing
    end
    outputs!(1)
    bf = gaussbatchfactor(sys, N; nrhs = ndir)
    # the stage operator is kept across steps while its first correction
    # contracts as the step's own does, and refreshed at the recorded
    # phases when it does not; a linear circuit assembles it once
    stale = true
    # the projection of the endpoint, linearized about the recorded
    # endpoint: `Z' (K dx - inj dI)` per column, solved with the small
    # Jacobians of the condition, then the rate along the directions
    pr = sys.gauss.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, m)
    dIh = isnothing(pr) ? nothing : Array(reshape(grid, nq, nt, ndir))
    Zti = isnothing(pr) ? nothing : pr.Zth*((sys.Lscale/phi0) .* transientinjection(p, targets))
    # the direction's current along the rows of the endpoint reading
    Qinji = isnothing(pr) ? nothing : endpointinjection(pr, (sys.Lscale/phi0) .* transientinjection(p, targets))
    function project!(k, endat!)
        endat!(pw.phip, k + 1)
        copyto!(pw.hphi, pw.phip)
        if !isempty(pr.directions)
            Ms = projectionmatrices(pr, pw.hphi)
            stepmul!(pw.zwork, pr.ZtL, dxnew)
            copyto!(pw.g, pw.zwork)
            if !isempty(pr.pj)
                stepmul!(pw.phim, pr.RJp, dxnew)
                hphim = Array(pw.phim)
                for j in 1:N
                    cols = (j - 1)*ndir + 1:j*ndir
                    view(pw.g, :, cols) .+= transpose(pr.RJZl)*(pr.lmoljp .* cos.(view(pw.hphi, :, j)) .* view(hphim, :, cols))
                end
            end
            dIz = repeat(Zti*view(dIh, :, k + 1, :), 1, N)
            nl2 > 0 && (dIz .+= pr.Ztline*Array(ddlinevalues))
            isnothing(rwt) || (dIz .+= pr.Ztblock*restingwaves(sys, rwt.z))
            pw.g .-= dIz
            for j in 1:N
                cols = (j - 1)*ndir + 1:j*ndir
                view(pw.alpha, :, cols) .= -(Ms[j] \ view(pw.g, :, cols))
            end
            copyto!(pw.dalpha, pw.alpha)
            stepmul!(pw.work, pr.Z, pw.dalpha)
            dxnew .+= pw.work
            projectrate!(dv, pr, pw, gc, h, stage(d, 1), stage(d, 2), dxnew, dx)
        end
        if !isempty(pr.readrows) || !isempty(pr.auxrows)
            # the linearized reading: the junction term is the stiffness at
            # the endpoint times the flux direction's phase, the drive the
            # direction's current at the endpoint time along the rows
            junction = if isempty(pr.pj)
                zeros(0, m)
            else
                stepmul!(pw.phim, pr.RJp, dxnew)
                hphim = Array(pw.phim)
                reshape(pr.lmoljp .* reshape(cos.(pw.hphi), :, 1, N) .* reshape(hphim, :, ndir, N), :, m)
            end
            dIq = repeat(Qinji*view(dIh, :, k + 1, :), 1, N)
            nl2 > 0 && (dIq .+= pr.Qline*Array(ddlinevalues))
            isnothing(rwt) || (dIq .+= pr.Qblock*restingwaves(sys, rwt.z))
            endpointread!(dv, dxnew, pr, pw, junction, dIq)
        end
        nothing
    end
    for window in responsewindows(sol, problems, sys, true)
        krange, phaseat!, endat! = window.steps, window.phases!, window.endphases!
        isnothing(window.replay) || window.replay()
      for k in krange
        # the recorded stage phases of every condition, stage last
        phaseat!(phi, k)
        # the right hand side of the linearized stages,
        #   r_i = -(L + J_i') dx_n + (A^{-1} 1)_i C dv_n / h + db_i,
        # the current at the stage time read off the grid, the same for
        # every condition
        stepmul!(lx, sys.L, dx)
        stepmul!(cdv, sys.C, dv)
        for i in 1:2
            batchjunctionproduct!(work, sys, view(phi, :, :, i), dx, jwork)
            if staged
                stagecurrent .= view(dS, :, i, k, :)
            else
                indices, weights = gaussstencil(k, nt, gc.c[i])
                fill!(stagecurrent, 0)
                for (idx, wgt) in zip(indices, weights)
                    stagecurrent .+= wgt .* view(dI, :, idx, :)
                end
            end
            stepmul!(stageinjection, injection, stagecurrent)
            ri = stage(r, i)
            reshape(injectionall, n, ndir, N) .= reshape(stageinjection, n, ndir, 1)
            ri .= injectionall .+ (gc.ainvone[i]/h) .* cdv .- lx .- work
            if nl2 > 0
                lineread!(sol.times[k] + gc.c[i]*h, npre + k - 1)
                stepmul!(linework, sys.lineinjection, ddlinevalues)
                ri .+= linework
            end
        end
        if !isnothing(rwt)
            # the blocks' states and the state's part of the stage currents
            # carry part of the right hand side
            stage(fullstages, 1) .= dx
            stage(fullstages, 2) .= dx
            rationalsources!(rwt, sys, zerostages, fullstages)
            r .+= rwt.source
        end
        stale && gaussbatchjacobian!(bf, sys, phi, cosphi)
        contraction = gaussbatchstagesolve!(d, r, sys, gc, bf, phi, false, res, dd, cwork, gwork, lwork, jwork, rc, zc; rw = rwt)
        stale = nj > 0 && contraction > 0.25
        if !isnothing(rwt)
            stage(fullstages, 1) .= dx .+ stage(d, 1)
            stage(fullstages, 2) .= dx .+ stage(d, 2)
            rationalstates!(rwt, sys, d, fullstages)
        end
        d1, d2 = stage(d, 1), stage(d, 2)
        dxnew .= dx .+ gc.ex[1] .* d1 .+ gc.ex[2] .* d2
        dv .+= (gc.ev[1]/h) .* d1 .+ (gc.ev[2]/h) .* d2
        nl2 > 0 && lineread!(sol.times[k + 1], npre + k - 1)
        isnothing(pr) || project!(k, endat!)
        if nl2 > 0
            # the waves leaving the ports at the endpoint into the history
            stepmul!(dlinerates, sys.linegather, dv)
            view(dwaves, :, ringslot(npre + k, nring), :) .= phi0 .* dlinerates ./ sqrtz .- ddlinevalues .* sqrtz ./ 2
        end
        copyto!(dx, dxnew)
        outputs!(k + 1)
      end
    end
    KernelAbstractions.synchronize(backend)
    single = ndims(currents) == 2
    shape = a -> begin
        isnothing(a) && return nothing
        b = reshape(a, size(a)[1:end-1]..., ndir, N)
        sol isa TransientSolution ? (single ? reshape(b, size(b)[1:end-2]...) : reshape(b, size(b)[1:end-1]...)) :
            (single ? reshape(b, size(b)[1:end-2]..., N) : b)
    end
    return (; voltage = shape(voltage), incident = shape(incident), outgoing = shape(outgoing),
        finalflux = shape(copy(dx)), finalrate = shape(copy(dv)))
end

# the adjoint of a Gauss-Legendre solve or batch, the weights shared by
# the conditions, the objectives of a condition contiguous in the
# columns the sinks receive; with a stage sink the stage multipliers go
# to it as they are, and the grid columns hold what the grid reads, the
# feedthrough and the projection; `transientadjoint` checked the
# arguments
function gaussbatchadjoint(sol, weights, quantity::Symbol, targets, sys::TransientSystem, sink, stagesink = nothing)
    problems, phases, x0s, v0s = batchview(sol)
    N = length(problems)
    p = first(problems)
    backend = sys.backend
    gc = sys.gauss.coefficients
    n, np, nt = length(p), length(p.portimpedances), length(sol.times)
    nq, nobj = length(targets), (ndims(weights) == 3 ? size(weights, 3) : 1)
    m = nobj*N
    h = sys.h
    w = tobackend(backend, reshape(Float64.(collect(weights)), np, nt, nobj))
    cv, cd = outputcoefficients(sys, quantity)
    injectiont = devicesparse(sparse(transpose((sys.Lscale/phi0) .* transientinjection(p, targets))), backend)
    tp = targetports(p, targets)
    portmapt = devicesparse(sparse([k for (k, q) in enumerate(tp) if q > 0], [q for q in tp if q > 0],
        ones(count(>(0), tp)), nq, np), backend)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    directwork, targetwork, objectivework = allocate(np, nobj), allocate(nq, m), allocate(nq, nobj)
    feedthrough! = (column, k) -> begin
        directwork .= cd .* view(w, :, k, :)
        stepmul!(objectivework, portmapt, directwork)
        reshape(column, nq, nobj, N) .+= reshape(objectivework, nq, nobj, 1)
        nothing
    end
    currents = isnothing(sink) ? allocate(nq, nt, m) : nothing
    store = isnothing(sink) ? (k, values) -> (copyto!(view(currents, :, k, :), values); nothing) : sink
    ring = CurrentRing(backend, nq, m, store, feedthrough!)
    xbar, vbar, work, lmu, cmu = [allocate(n, m) for _ in 1:5]
    wst, mu, res, dd, cwork, gwork, lwork = [allocate(n, m, 2) for _ in 1:7]
    nj = length(sys.lmolj)
    jwork, phi = allocate(nj, m), allocate(nj, N, 2)
    cosphi = KernelAbstractions.zeros(backend, ComplexF64, nj, N)
    rc, zc = [KernelAbstractions.zeros(backend, ComplexF64, n, m) for _ in 1:2]
    portwork = allocate(np, nobj)
    objectivestate = allocate(n, nobj)
    function output!(k)
        portwork .= cv .* view(w, :, k, :)
        stepmul!(objectivestate, sys.portst, portwork)
        reshape(vbar, n, nobj, N) .+= phi0 .* reshape(objectivestate, n, nobj, 1)
    end
    output!(nt)
    bf = gaussbatchfactor(sys, N; nrhs = nobj)
    stale = true
    # the lines: the ring of the cotangent of the wave leaving each port,
    # which every read scatters onto and each endpoint hands to the rates
    # and to the waves it read, then clears for the column that will
    # reuse its slot; a read at a step touches only columns before the
    # step's endpoint, so a column is complete when its step is reached,
    # and at the end the ring holds the prehistory's cotangent
    lines = p.lines
    nl2 = 2length(lines)
    npre = lineprehistory(p, h)
    nring = linering(npre)
    t0 = first(sol.times)
    tpre = linestart(t0, npre, h)
    abar = allocate(nl2, nring, m)
    forcedbar, ratebar, linework = allocate(nl2, m), allocate(nl2, m), allocate(n, m)
    stencil, dstencil = zeros(8, nl2), allocate(8, nl2)
    far, readscale, sqrtz = linetables(p, backend)
    # the rational blocks: the cotangent of their states, and the part of
    # the flux's cotangent the stage values carry
    rwa = isempty(sys.gauss.rational) ? nothing : rationalwork(p, backend, n, m)
    xextra = isnothing(rwa) ? nothing : allocate(n, m)
    # the cotangent of a forced current at time `t` scattered onto the
    # far port's samples it read, `dq = 2 forced / sqrt(Z)`, with `sign`
    function linescatter!(t, accepted, sign)
        linestencil!(stencil, lines, t, tpre, h, accepted)
        copyto!(dstencil, stencil)
        wavescatter!(abar, forcedbar, dstencil, far, readscale, sign, backend)
        nothing
    end
    # the forced currents' cotangent from a multiplier on the equations
    function forced!(mu)
        stepmul!(forcedbar, sys.linegather, mu)
        forcedbar .*= sys.Lscale/phi0
        nothing
    end
    # the projection transposed: the rate along the directions first, the
    # part of the objective's rate the cubic carries to the stages and to
    # the state, then the endpoint, whose multiplier `gamma` reads the
    # current at time `k + 1` and moves the objective by `K Z gamma`
    pr = sys.gauss.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, m)
    vrecbar = isnothing(pr) ? nothing : allocate(n, m)
    er = gc.endrate
    function projectbefore!(k, endat!)
        endat!(pw.phip, k + 1)
        copyto!(pw.hphi, pw.phip)
        # the reading of the index one unknowns transposed, last in the
        # step so first here: its row vector moves the current at the
        # endpoint time
        cosp = tobackend(backend, cos.(pw.hphi))
        if endpointreadtranspose!(vbar, xbar, pr, pw, sys, cosp, nobj, N)
            stepmul!(targetwork, injectiont, pw.work)
            ringadd!(ring, k + 1, targetwork, -1.0)
            if nl2 > 0
                forced!(pw.work)
                linescatter!(sol.times[k + 1], npre + k - 1, -1.0)
            end
            # the resting waves the reading saw came from the states
            isnothing(rwa) || restingwavesbar!(rwa, sys, pw.work, -1.0)
        end
        isempty(pr.directions) && return nothing
        Ms = projectionmatrices(pr, pw.hphi)
        # the rate along the directions transposed: `v += Z R (vrec - v)`
        stepmul!(pw.zwork, pr.Zt, vbar)
        stepmul!(vrecbar, pr.Zratet, pw.zwork)
        vbar .-= vrecbar
        xbar .+= (er[4]/h) .* vrecbar
        # the Newton step transposed: the multiplier from the flux along
        # the directions, carried by the constraints' rows to the flux,
        # the current, the lines' forced currents and the resting waves
        stepmul!(pw.zwork, pr.Zt, xbar)
        copyto!(pw.g, pw.zwork)
        for j in 1:N
            cols = (j - 1)*nobj + 1:j*nobj
            view(pw.alpha, :, cols) .= transpose(Ms[j]) \ view(pw.g, :, cols)
        end
        copyto!(pw.dalpha, pw.alpha)
        stepmul!(pw.work, pr.Zl, pw.dalpha)
        stepmul!(targetwork, injectiont, pw.work)
        ringadd!(ring, k + 1, targetwork, 1.0)
        if nl2 > 0
            forced!(pw.work)
            linescatter!(sol.times[k + 1], npre + k - 1, 1.0)
        end
        isnothing(rwa) || restingwavesbar!(rwa, sys, pw.work, 1.0)
        stepmul!(lmu, sys.Lt, pw.work)
        xbar .-= lmu
        if !isempty(pr.pj)
            stepmul!(pw.phim, pr.RJp, pw.work)
            reshape(pw.phim, :, nobj, N) .*= pr.lmoljpdev .* reshape(cosp, :, 1, N)
            stepmul!(pw.work2, pr.RJpt, pw.phim)
            xbar .-= pw.work2
        end
        nothing
    end
    for window in responsewindows(sol, problems, sys, false)
        krange, phaseat!, endat! = window.steps, window.phases!, window.endphases!
        isnothing(window.replay) || window.replay()
      for k in reverse(krange)
        phaseat!(phi, k)
        if nl2 > 0
            # the wave that left each port at the endpoint: its cotangent
            # goes to the rate across the port and, with the opposite
            # sign, to the far port's samples the arriving wave read
            slot = ringslot(npre + k, nring)
            ratebar .= phi0 .* view(abar, :, slot, :) ./ sqrtz
            forcedbar .= .- view(abar, :, slot, :) .* sqrtz ./ 2
            view(abar, :, slot, :) .= 0
            stepmul!(linework, sys.linescatter, ratebar)
            vbar .+= linework
            linescatter!(sol.times[k + 1], npre + k - 1, 1.0)
        end
        isnothing(pr) || projectbefore!(k, endat!)
        for i in 1:2
            stage(wst, i) .= gc.ex[i] .* xbar .+ (gc.ev[i]/h) .* vbar
            (isnothing(pr) || isempty(pr.directions)) || (stage(wst, i) .+= (er[i + 1]/h) .* vrecbar)
        end
        if !isnothing(rwa)
            # the states after the step depend on the stage unknowns
            # through the update: their cotangent reaches the stages
            fill!(xextra, 0)
            statesbartostages!(rwa, sys, wst, xextra)
        end
        stale && gaussbatchjacobian!(bf, sys, phi, cosphi)
        contraction = gaussbatchstagesolve!(mu, wst, sys, gc, bf, phi, true, res, dd, cwork, gwork, lwork, jwork, rc, zc; rw = rwa)
        stale = nj > 0 && contraction > 0.25
        # the states before the step: through the update, and through
        # the reflected waves the multipliers weigh, whose incident
        # waves' value part reaches the state's flux as well
        isnothing(rwa) || statesbarstep!(rwa, sys, mu, xextra)
        for i in 1:2
            mui = stage(mu, i)
            stepmul!(targetwork, injectiont, mui)
            if isnothing(stagesink)
                indices, wgts = gaussstencil(k, nt, gc.c[i])
                for (idx, wgt) in zip(indices, wgts)
                    ringadd!(ring, idx, targetwork, wgt)
                end
            else
                stagesink(k, i, targetwork)
            end
            if nl2 > 0
                forced!(mui)
                linescatter!(sol.times[k] + gc.c[i]*h, npre + k - 1, 1.0)
            end
            stepmul!(lmu, sys.Lt, mui)
            batchjunctionproduct!(work, sys, view(phi, :, :, i), mui, jwork)
            xbar .-= lmu .+ work
            stepmul!(cmu, sys.C, mui)
            vbar .+= (gc.ainvone[i]/h) .* cmu
        end
        (isnothing(pr) || isempty(pr.directions)) || (xbar .-= (er[4]/h) .* vrecbar)
        isnothing(rwa) || (xbar .+= xextra)
        output!(k)
        k + 3 <= nt && ringemit!(ring, k + 3)
      end
    end
    for j in min(3, nt):-1:1
        ringemit!(ring, j)
    end
    KernelAbstractions.synchronize(backend)
    shape = a -> begin
        isnothing(a) && return nothing
        b = reshape(a, size(a)[1:end-1]..., nobj, N)
        sol isa TransientSolution ? (ndims(weights) == 2 ? reshape(b, size(b)[1:end-2]...) : reshape(b, size(b)[1:end-1]...)) :
            (ndims(weights) == 2 ? reshape(b, size(b)[1:end-2]..., N) : b)
    end
    initialwaves = nl2 > 0 ? shape(readhistory!(allocate(nl2, npre, m), abar, 1)) : nothing
    initialstates = isnothing(rwa) ? nothing : shape(copy(rwa.z))
    return (; currents = shape(currents), initialflux = shape(xbar), initialrate = shape(vbar), initialwaves, initialstates)
end

"""
    transienttangent(batch::TransientBatchSolution, currents; targets, initialstate)

The tangent of every condition of a batch along the same currents, all
conditions on one pass: the arrays of the single tangent with the
conditions as the trailing dimension. The initial perturbation is one
pair for every condition, or one per condition with the conditions as a
trailing dimension of the pair's arrays.
"""
function transienttangent(b::TransientBatchSolution, currents::AbstractArray{<:Real};
        targets = porttargets(first(b.problems)), initialstate = nothing, factorization = nothing, reuse = nothing,
        outputsink = nothing)
    recordedsolution(b)
    p = first(b.problems)
    backend = KernelAbstractions.get_backend(b.finalflux)
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, b.dt, b.method, backend, fact)
    nq, nt = length(targets), length(b.times)
    currentshape(currents, nq, nt)
    return gaussbatchtangent(b, currents, targets, initialstate, sys, outputsink)
end

"""
    transientadjoint(batch::TransientBatchSolution, weights; quantity, targets, sink)

The adjoint of every condition of a batch for the same weights, all
conditions on one pass: the arrays of the single adjoint with the
conditions as the trailing dimension, and a sink's columns holding the
objectives of a condition contiguously, condition after condition.
"""
function transientadjoint(b::TransientBatchSolution, weights::AbstractArray{<:Real};
        quantity::Symbol = :outgoing, targets = porttargets(first(b.problems)), factorization = nothing,
        reuse = nothing, sink = nothing, stagesink = nothing)
    recordedsolution(b)
    p = first(b.problems)
    backend = KernelAbstractions.get_backend(b.finalflux)
    fact = isnothing(factorization) ? transientfactorization(backend) : factorization
    sys = transientsystem(reuse, p, b.dt, b.method, backend, fact)
    np, nt = length(p.portimpedances), length(b.times)
    ndims(weights) in (2, 3) && size(weights, 1) == np && size(weights, 2) == nt || throw(DimensionMismatch(
        lazy"weights must have one row per port ($(np)), one column per recorded time ($(nt)) and optionally a third dimension of objectives."))
    all(isfinite, weights) || throw(ArgumentError("the weights must be finite."))
    return gaussbatchadjoint(b, weights, quantity, targets, sys, sink, stagesink)
end

# the stage residual of a batch on `(n, N, 2)` stage
# increments with the states as the columns of `x`, the norm per column
function gaussbatchresidual!(norms, residual, sys::TransientSystem, gc::GaussCoefficients, delta, x, lx, X,
        phi, junction, jwork, cwork, gwork, rhs, colnorm, rw = nothing)
    h = sys.h
    X .= x .+ delta
    # the rational blocks' reflected waves at the stages, from their states
    # and the incident waves, as sources on their rows
    isnothing(rw) || rationalsources!(rw, sys, delta, X)
    stepmul!(stage(cwork, 1), sys.C, stage(delta, 1)); stepmul!(stage(cwork, 2), sys.C, stage(delta, 2))
    stepmul!(stage(gwork, 1), sys.G, stage(delta, 1)); stepmul!(stage(gwork, 2), sys.G, stage(delta, 2))
    stepmul!(stage(phi, 1), sys.RJ, stage(X, 1)); stepmul!(stage(phi, 2), sys.RJ, stage(X, 2))
    jwork .= sys.lmolj .* sin.(phi)
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
    return nothing
end

# The work of the rational blocks over `N` conditions, all on the
# backend: the states, the stage stacked increments, values and sources,
# `[stage 1; stage 2]`, on which the grouped operators of the coupling
# act, and the sources per stage the residual reads.
struct RationalWork{M, A}
    z::M
    dstack::M
    xstack::M
    sstack::M
    tstack::M
    zwork::M
    zwork2::M
    source::A
end

function rationalwork(p::TransientProblem, backend, n::Int, N::Int)
    allocate = (dims...) -> KernelAbstractions.zeros(backend, Float64, dims...)
    nz = blockstates(p)
    return RationalWork(allocate(nz, N), allocate(2n, N), allocate(2n, N), allocate(2n, N), allocate(2n, N),
        allocate(nz, N), allocate(nz, N), allocate(n, N, 2))
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
# the state, the states: `M_d delta + M_x X + M_s z` on the stacked stages
function rationalsources!(rw::RationalWork, sys::TransientSystem, delta, X; withstate::Bool = true)
    cp = sys.gauss.coupling
    stack!(rw.dstack, delta)
    stack!(rw.xstack, X)
    stepmul!(rw.sstack, cp.Md, rw.dstack)
    stepmul!(rw.tstack, cp.Mx, rw.xstack)
    rw.sstack .+= rw.tstack
    if withstate
        stepmul!(rw.tstack, cp.Ms, rw.z)
        rw.sstack .+= rw.tstack
    end
    unstack!(rw.source, rw.sstack)
    return rw
end

# the reflected waves at the stages from the states alone, `M_s z`, the
# part of a linearized step's right hand side its states carry
function rationalstatesource!(rw::RationalWork, sys::TransientSystem, z)
    stepmul!(rw.sstack, sys.gauss.coupling.Ms, z)
    unstack!(rw.source, rw.sstack)
    return rw
end

# the transpose of the stages' rational coupling: the multipliers `mu`
# on the rows carried to the stage unknowns, into `rw.source`
function rationalsourcestranspose!(rw::RationalWork, sys::TransientSystem, mu)
    cp = sys.gauss.coupling
    stack!(rw.dstack, mu)
    stepmul!(rw.sstack, cp.Mdt, rw.dstack)
    stepmul!(rw.tstack, cp.Mxt, rw.dstack)
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
    stack!(rw.dstack, mu)
    stepmul!(rw.zwork, cp.Ezt, rw.z)
    stepmul!(rw.zwork2, cp.Mst, rw.dstack)
    rw.z .= rw.zwork .+ rw.zwork2
    stepmul!(rw.tstack, cp.Mxt, rw.dstack)
    stacksum!(xextra, rw.tstack)
    return rw
end
function restingwavesbar!(rw::RationalWork, sys::TransientSystem, work, sign)
    stepmul!(rw.zwork, sys.gauss.coupling.Rt, work)
    rw.z .+= sign .* rw.zwork
    return rw
end

# the reflected waves of the rational parts at rest at their states,
# `C z` per port on the host, the source the endpoint reading and the
# consistency check see on the blocks' rows
restingwaves(sys::TransientSystem, z) = sys.gauss.coupling.Cblk*Array(z)

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
# the prehistory before and at the start, constant at the initial waves
# in a solve and a perturbation in a response, and the ring long enough
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
    trialjunction::A
    cwork::A
    gwork::A
    xnew::M
    cv::M
    lx::M
    b1::M
    b2::M
    bend::M
    phi::A
    trialphi::A
    jwork::A
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
    rtol::Float64
    atol::Float64
    maxiters::Int
    stalefailed::Bool
    corrections::Int
    factorizations::Int
    retries::Int
end

function gaussstepper(sys::TransientSystem, problems, rtol, atol, maxiters, bf)
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
    cosphi = KernelAbstractions.zeros(backend, ComplexF64, nj, N)
    rc, zc = [KernelAbstractions.zeros(backend, ComplexF64, n, N) for _ in 1:2]
    colnorm = allocate(N)
    hostvalues = zeros(nd, N)
    values = tobackend(backend, zeros(nd, N))
    portwork = allocate(np, N)
    rw = isempty(sys.gauss.rational) ? nothing : rationalwork(p, backend, n, N)
    baseresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, phi, junction, jwork, cwork, gwork, rhs, colnorm, rw)
    trialresidual! = (norms, r, D) -> gaussbatchresidual!(norms, r, sys, gc, D, x, lx, X, trialphi, trialjunction, jwork, cwork, gwork, rhs, colnorm, rw)
    refresh! = () -> (gaussbatchjacobian!(bf, sys, phi, cosphi); nothing)
    solve! = (c, r) -> (gaussbatchtransform!(c, r, bf, gc, rc, zc); (false, 0))
    accept! = mask -> (maskcolumns!(phi, trialphi, mask); maskcolumns!(junction, trialjunction, mask); nothing)
    pr = sys.gauss.projection
    pw = isnothing(pr) ? nothing : projectionwork(pr, backend, n, N, N)
    nl = 2length(p.lines)
    far, readscale, sqrtz = linetables(p, backend)
    npre = lineprehistory(p, sys.h)
    return GaussStepper(sys, problems, N, x, v, X, delta, lastdelta, residual, trial, trialresidual, correction, rhs,
        junction, trialjunction, cwork, gwork, xnew, cv, lx, b1, b2, bend, phi, trialphi, jwork, cosphi, rc, zc,
        colnorm, hostvalues, values, portwork, bf, baseresidual!, trialresidual!, refresh!, solve!, accept!, pw,
        allocate(nl, linering(npre), N), npre, Float64[], 0, allocate(nl, N), zeros(8, nl), allocate(8, nl), far, readscale, sqrtz,
        allocate(nl, N), allocate(n, N), rw,
        NewtonWork(backend, N), zeros(N), Float64(rtol), Float64(atol), Int(maxiters), false, 0, 0, 0)
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

# the phases of the projected junctions at a state, for the record
function projectedphases!(st::GaussStepper, x)
    pr = st.sys.gauss.projection
    isempty(pr.pj) || stepmul!(st.pw.phip, pr.RJp, x)
    return st.pw.phip
end

# The projection of every condition's endpoint onto the algebraic
# constraints along the projected directions, and the rate along them.
# Newton on the coefficients `alpha` of `Z` per condition, the residual
# `Z' (L x + J(x) - b(t))` and its Jacobians `Z' (L + J'(x)) Z` on the host
# from two small products on the backend, to the step's tolerance, the
# predictor moved with the state; a linear constraint is met by one
# correction. Then the rate along `Z` is replaced by the derivative at
# the endpoint of the cubic through the state, the stages and the
# projected endpoint, all of which satisfy the constraint. Leaves the
# projected junctions' phases at the endpoint in the work buffer.
function gaussproject!(st::GaussStepper, t)
    sys, pw = st.sys, st.pw
    pr = sys.gauss.projection
    gc = sys.gauss.coefficients
    h = sys.h
    N = st.N
    drivevalues!(st.hostvalues, st.problems, t)
    isempty(sys.problem.lines) || linevalues!(st, t)
    linedrive = isempty(sys.problem.lines) ? nothing : Array(st.linevalues)
    if !isempty(pr.directions)
        gb = pr.Ztinj*st.hostvalues .+ pr.Ztconstant
        isnothing(linedrive) || (gb .+= pr.Ztline*linedrive)
        isnothing(st.rw) || (gb .+= pr.Ztblock*restingwaves(sys, st.rw.z))
        converged = false
        for iteration in 1:st.maxiters + 1
            projectedphases!(st, st.xnew)
            copyto!(pw.hphi, pw.phip)
            stepmul!(pw.zwork, pr.ZtL, st.xnew)
            copyto!(pw.g, pw.zwork)
            pw.g .+= transpose(pr.RJZl)*(pr.lmoljp .* sin.(pw.hphi)) .- gb
            converged = all(j -> maximum(abs, view(pw.g, :, j)) <= st.tolerance[j], 1:N)
            (converged || iteration > st.maxiters) && break
            for (j, M) in enumerate(projectionmatrices(pr, pw.hphi))
                pw.alpha[:, j] .= -(M \ view(pw.g, :, j))
            end
            copyto!(pw.dalpha, pw.alpha)
            stepmul!(pw.work, pr.Z, pw.dalpha)
            st.xnew .+= pw.work
            stage(st.lastdelta, 1) .-= pw.work
            stage(st.lastdelta, 2) .-= pw.work
            isempty(pr.pj) && (converged = true; break)
        end
        converged || error(lazy"the projection onto the algebraic constraints at t = $(t) s did not converge; reduce dt or check the state.")
        projectrate!(st.v, pr, pw, gc, h, stage(st.delta, 1), stage(st.delta, 2), st.xnew, st.x)
    end
    # the index one unknowns of the endpoint from their equations
    if !isempty(pr.readrows) || !isempty(pr.auxrows)
        projectedphases!(st, st.xnew)
        copyto!(pw.hphi, pw.phip)
        gq = pr.Qinj*st.hostvalues .+ pr.Qconst
        isnothing(linedrive) || (gq .+= pr.Qline*linedrive)
        isnothing(st.rw) || (gq .+= pr.Qblock*restingwaves(sys, st.rw.z))
        endpointread!(st.v, st.xnew, pr, pw, pr.lmoljp .* sin.(pw.hphi), gq)
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
    gaussbatchjacobian!(st.bf, st.sys, st.phi, st.cosphi)
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
        st.tolerance, st.maxiters, st.stalefailed, length(sys.lmolj) > 0, false, st.newtonwork;
        simplified = true, accept! = st.accept!)
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

# the integration of a batch under the Gauss-Legendre rule
function gaussbatchintegrate(sys::TransientSystem, problems, t0, tf, nsteps, initialstates,
        saveevery, record, checkpointevery, rtol, atol, maxiters, reuse)
    savephases, savestates, savecheckpoints = recordlevel(record, saveevery)
    p = sys.problem
    backend = sys.backend
    N = length(problems)
    n, np = length(p), length(p.portimpedances)
    nj = length(sys.lmolj)
    h = sys.h
    bf = if !isnothing(reuse) && reuse.factor isa GaussBatchFactor && reuse.factor.ncolumns == N
        reuse.factor
    else
        gaussbatchfactor(sys, N)
    end
    st = gaussstepper(sys, problems, rtol, atol, maxiters, bf)
    nl = 2length(p.lines)
    x0, v0, w0 = zeros(n, N), zeros(n, N), zeros(nl, N)
    for (j, state) in enumerate(initialstates)
        xj, vj = state
        (length(xj) == n && length(vj) == n) || throw(DimensionMismatch(lazy"the state has $(n) entries; use transientstate."))
        (all(isfinite, xj) && all(isfinite, vj)) || throw(ArgumentError("the initial state must be finite."))
        x0[:, j] .= Float64.(collect(xj))
        v0[:, j] .= Float64.(collect(vj))
        wj = initialwaves(state, p)
        length(wj) == nl && all(isfinite, wj) || throw(DimensionMismatch(lazy"the state needs $(nl) finite line waves; use transientstate."))
        w0[:, j] .= wj
    end
    z0 = zeros(blockstates(p), N)
    if !isnothing(st.rw)
        for (j, state) in enumerate(initialstates)
            zj = initialblockstates(state, p)
            length(zj) == size(z0, 1) && all(isfinite, zj) || throw(DimensionMismatch(lazy"the state needs $(size(z0, 1)) finite block states; use transientstate."))
            z0[:, j] .= zj
        end
        copyto!(st.rw.z, z0)
    end
    setstate!(st, x0, v0, nothing)
    # the grid, and the history of the line waves before the start,
    # constant at the initial waves; the record of the waves leaving
    # every port at every step, or with checkpoints the history before
    # each of them
    times = [k == nsteps ? tf : t0 + k*h for k in 0:nsteps]
    npre = st.npre
    tail = KernelAbstractions.zeros(backend, Float64, nl, npre, N)
    tail .= reshape(tobackend(backend, w0), nl, 1, N)
    sethistory!(st, times, tail, 0)
    # the record of the waves is kept without checkpoints, and with them
    # when the history before every checkpoint would outweigh it, as it
    # does once the longest delay exceeds the checkpoint interval
    K = checkpointevery > 0 ? Int(checkpointevery) : max(1, round(Int, sqrt(nsteps)))
    nc = savecheckpoints ? cld(nsteps, K) : 0
    tails = savecheckpoints && npre*nc <= nsteps + 1
    waves = nl > 0 && (!savecheckpoints || !tails) ? KernelAbstractions.zeros(backend, Float64, nl, nsteps + 1, N) : nothing
    isnothing(waves) || (view(waves, :, 1, :) .= tobackend(backend, w0))
    # the check of the algebraic equations at the start, per condition
    # under its own drives, the lines carrying their initial waves
    resting = isnothing(st.rw) ? nothing : restingwaves(sys, st.rw.z)
    for j in 1:N
        violation, _ = transientconsistency(sys, view(st.x, :, j), view(st.v, :, j), t0, problems[j],
            lineforcing(p, arrivingwaves(p, w0[:, j])), isnothing(resting) ? nothing : resting[:, j])
        violation <= atol + rtol || throw(ArgumentError(
            lazy"the initial state of condition $(j) violates the algebraic equations of the circuit along a direction without capacitance to ground (a node no capacitor touches, a capacitive island, a coupled inductor or gauge row); supply a consistent transientstate, or start the drive from an equilibrium."))
    end
    # the saved outputs
    nsaved = cld(nsteps, saveevery) + 1
    savedtimes = Vector{Float64}(undef, nsaved)
    voltage, incident, outgoing = [KernelAbstractions.allocate(backend, Float64, np, nsaved, N) for _ in 1:3]
    phases = savephases ? KernelAbstractions.zeros(backend, Float64, nj, 2, nsaved, N) : nothing
    pr = sys.gauss.projection
    npj = isnothing(pr) ? 0 : length(pr.pj)
    endphases = savephases && npj > 0 ? KernelAbstractions.zeros(backend, Float64, npj, nsaved, N) : nothing
    isnothing(endphases) || copyto!(view(endphases, :, 1, :), projectedphases!(st, st.x))
    flux = savestates ? KernelAbstractions.allocate(backend, Float64, n, nsaved, N) : nothing
    rate = savestates ? KernelAbstractions.allocate(backend, Float64, n, nsaved, N) : nothing
    nzs = blockstates(p)
    states = savestates && nzs > 0 ? KernelAbstractions.zeros(backend, Float64, nzs, nsaved, N) : nothing
    isnothing(states) || (view(states, :, 1, :) .= st.rw.z)
    # the checkpoints: the state and the stage predictor every K steps,
    # from which a response replays the steps between, and the history of
    # the line waves before each, unless the record is kept instead
    checkpoints = savecheckpoints ? (; every = K, flux = KernelAbstractions.zeros(backend, Float64, n, nc, N),
        rate = KernelAbstractions.zeros(backend, Float64, n, nc, N),
        increment = KernelAbstractions.zeros(backend, Float64, n, 2, nc, N),
        states = KernelAbstractions.zeros(backend, Float64, nzs, nc, N),
        waves = KernelAbstractions.zeros(backend, Float64, nl, nl > 0 && tails ? npre : 0, tails ? nc : 0, N)) : nothing
    savedtimes[1] = t0
    batchdrivecurrent!(st.bend, st, t0)
    portwaves!(view(voltage, :, 1, :), view(incident, :, 1, :), view(outgoing, :, 1, :), sys, st.v, st.values, st.portwork)
    if savestates
        copyto!(view(flux, :, 1, :), st.x)
        copyto!(view(rate, :, 1, :), st.v)
    end
    initialflux, initialrate = copy(st.x), copy(st.v)
    saved = 1
    tprev = t0
    for step in 1:nsteps
        t = step == nsteps ? tf : t0 + step*h
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
        isnothing(waves) || (view(waves, :, step + 1, :) .= view(st.waves, :, ringslot(npre + step, size(st.waves, 2)), :))
        if step % saveevery == 0 || step == nsteps
            saved += 1
            savedtimes[saved] = t
            batchdrivecurrent!(st.bend, st, t)
            portwaves!(view(voltage, :, saved, :), view(incident, :, saved, :), view(outgoing, :, saved, :), sys, st.v, st.values, st.portwork)
            if savephases
                view(phases, :, 1, saved, :) .= stage(st.phi, 1)
                view(phases, :, 2, saved, :) .= stage(st.phi, 2)
            end
            isnothing(endphases) || copyto!(view(endphases, :, saved, :), st.pw.phip)
            if savestates
                copyto!(view(flux, :, saved, :), st.x)
                copyto!(view(rate, :, saved, :), st.v)
                isnothing(states) || (view(states, :, saved, :) .= st.rw.z)
            end
        end
    end
    KernelAbstractions.synchronize(backend)
    isnothing(reuse) || (reuse.factor = bf)
    return TransientBatchSolution(problems, sys.method, h, savedtimes, voltage, incident, outgoing,
        phases, endphases, waves, flux, rate, checkpoints, initialflux, initialrate, copy(st.x), copy(st.v),
        states, w0, z0,
        (; steps = nsteps, newtoncorrections = st.corrections, factorizations = st.factorizations,
            retries = st.retries, kryloviterations = 0, rtol = Float64(rtol), atol = Float64(atol), maxiters = Int(maxiters)))
end

# The windows of steps a response walks and the phases of each. With a
# phase record, one window of every step reading the record; with
# checkpoints, one window per checkpoint, replayed from its state by the
# stepper into a buffer of the window's stage phases, walked forward by
# the tangent and backward by the adjoint. Each window names its steps,
# `phases!(buffer, k)` filling the `(nj, N, 2)` buffer with the phases of
# the step from time `k` to `k + 1`, `endphases!(buffer, k)` filling the
# `(npj, N)` buffer with the projected junctions' phases at time `k`, and
# the replay to run before the window, if any.
function responsewindows(sol, problems, sys::TransientSystem, forward::Bool)
    nt = length(sol.times)
    pr = sys.gauss.projection
    npj = isnothing(pr) ? 0 : length(pr.pj)
    if !isnothing(sol.phases)
        phases, endphases = batchview(sol)[2], batchview(sol)[5]
        phases! = (buffer, k) -> begin
            view(buffer, :, :, 1) .= view(phases, :, 1, k + 1, :)
            view(buffer, :, :, 2) .= view(phases, :, 2, k + 1, :)
            nothing
        end
        endphases! = (buffer, k) -> (npj > 0 && copyto!(buffer, view(endphases, :, k, :)); nothing)
        return [(; steps = 1:nt - 1, phases!, endphases!, replay = nothing)]
    end
    cp = sol.checkpoints
    cps = sol isa TransientSolution ? (; every = cp.every, flux = reshape(cp.flux, size(cp.flux)..., 1),
        rate = reshape(cp.rate, size(cp.rate)..., 1), increment = reshape(cp.increment, size(cp.increment)..., 1),
        states = reshape(cp.states, size(cp.states)..., 1), waves = reshape(cp.waves, size(cp.waves)..., 1)) : cp
    K, nc = cps.every, size(cps.flux, 2)
    N = length(problems)
    nj = length(sys.lmolj)
    st = gaussstepper(sys, problems, sol.stats.rtol, sol.stats.atol, sol.stats.maxiters, gaussbatchfactor(sys, N))
    nl = 2length(sys.problem.lines)
    npre = st.npre
    record = batchview(sol)[6]
    buffer = KernelAbstractions.zeros(sys.backend, Float64, nj, 2, K + 1, N)
    ebuffer = KernelAbstractions.zeros(sys.backend, Float64, npj, K + 2, N)
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
                # and the constant initial waves before it
                st.times = sol.times
                before = max(npre - kstart, 0)
                if before > 0
                    initial = KernelAbstractions.zeros(sys.backend, Float64, nl, before, N)
                    initial .= view(record, :, 1:1, :)
                    fillhistory!(st.waves, initial, kstart)
                end
                fillhistory!(st.waves, view(record, :, max(kstart - npre + 1, 1):kstart, :), kstart + before)
                st.accepted = npre + kstart - 1
            end
            npj > 0 && copyto!(view(ebuffer, :, 1, :), projectedphases!(st, st.x))
            for k in kstart:kend
                advance!(st, sol.times[k], sol.times[k + 1], k)
                view(buffer, :, 1, k - kstart + 2, :) .= stage(st.phi, 1)
                view(buffer, :, 2, k - kstart + 2, :) .= stage(st.phi, 2)
                npj > 0 && copyto!(view(ebuffer, :, k - kstart + 2, :), st.pw.phip)
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
        (; steps = kstart:kend, phases!, endphases!, replay)
    end
    return forward ? windows : reverse(windows)
end

# a batch of one is an ordinary solution
function unbatch(b::TransientBatchSolution)
    squeeze = a -> isnothing(a) ? nothing : reshape(a, size(a)[1:end-1]...)
    cp = isnothing(b.checkpoints) ? nothing : (; every = b.checkpoints.every, flux = squeeze(b.checkpoints.flux),
        rate = squeeze(b.checkpoints.rate), increment = squeeze(b.checkpoints.increment), states = squeeze(b.checkpoints.states),
        waves = squeeze(b.checkpoints.waves))
    return TransientSolution(b.problems[1], b.method, b.dt, b.times, squeeze(b.voltage), squeeze(b.incident),
        squeeze(b.outgoing), squeeze(b.phases), squeeze(b.endphases), squeeze(b.linewaves), squeeze(b.flux), squeeze(b.rate), cp, vec(b.initialflux),
        vec(b.initialrate), vec(b.finalflux), vec(b.finalrate), squeeze(b.blockstates), squeeze(b.initialwaves), squeeze(b.initialstates), b.stats)
end
