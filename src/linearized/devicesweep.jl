# The linearized sweep on a backend: the system matrices of a batch of
# frequencies assembled and solved together, the forward and adjoint
# solutions gathered at the ports, and the outputs reduced in place.

"""
    FrequencySweepPlan

Everything needed to assemble the linearized system matrix of
[`hblinsolve`](@ref) at many signal frequencies at once, on a backend.

The per-frequency assembly of [`assemblesystemmatrix!`](@ref) is

    A = AoLjnm + invLnm + im*Gnm*w - Cnm*w^2 + Amna0

with the frequency of each stored entry taken from its column's mode, and
with the stored value of every frequency dependent term conjugated where that
mode frequency is negative (see [`modevalue`](@ref)). Each stored entry is
therefore an independent quadratic in its own mode frequency, and the whole
assembly collapses to four constant coefficient vectors and one kernel:

    A[q] = cst[q] + sel(kinvL[q]) + im*wm*sel(kG[q]) - wm^2*sel(kC[q])

where `wm = ws + wpumpmodes[mode of q]` and `sel` conjugates when `wm < 0`.
The terms which share a frequency power and a conjugation rule are summed
into one coefficient at build time, `AoLjnm` with `Amna0`. Conjugation
distributes over that sum, so this is exact.

Because the coefficients do not depend on the signal frequency, one kernel
fills the stored values of a whole batch of frequencies, which is what a
uniform batch wants, whether cuDSS's or the batched block factorization's:
the batch shares one sparsity pattern and one symbolic analysis, and only
the values differ.

# Fields
- `colof`: the matrix column of each stored entry, from which its mode
    follows, so the kernel needs no search. For the forward structure this is
    also the structure's own column index array; for the transposed one it is
    not, so it is carried separately.
- `cst`, `kinvL`, `kG`, `kC`: the four coefficient vectors, in the stored
    order of whichever structure the plan was built for (compressed sparse
    row of the matrix for the forward plan, the matrix's own stored order
    for the adjoint one).
- `wpump`: the pump mode frequency offsets of the signal modes.
- `assemble!`, `backend`: the compiled assembly kernel and the backend it
    was compiled for.
- `nmodes`, `nnz`: the mode count the kernel takes a column's mode from,
    and the number of stored entries, which sizes the value matrix.
"""
struct FrequencySweepPlan{VI,VC,VR,K,B}
    colof::VI
    cst::VC
    kinvL::VC
    kG::VC
    kC::VC
    wpump::VR
    assemble!::K
    backend::B
    nmodes::Int
    nnz::Int
end

# one work item per (stored entry, frequency) pair of the batch
@kernel function sweepassemblykernel!(nzval, @Const(colof), @Const(cst),
        @Const(kinvL), @Const(kG), @Const(kC), @Const(wpump), @Const(ws),
        Nmodes, nnzA)
    gid = @index(Global)
    @inbounds begin
        q = gid
        e = (q - 1) % nnzA + 1
        b = (q - 1) ÷ nnzA + 1
        # a compressed column of the stored structure is a row of the system
        # matrix, so the stored row index is the matrix's column
        m = (Int(colof[e]) - 1) % Nmodes + 1
        wm = ws[b] + wpump[m]
        neg = wm < 0
        vL = neg ? conj(kinvL[e]) : kinvL[e]
        vG = neg ? conj(kG[e]) : kG[e]
        vC = neg ? conj(kC[e]) : kC[e]
        nzval[q] = cst[e] + vL + (im*wm)*vG - (wm*wm)*vC
    end
end

"""
    cansweepondevice(lsys::HBLinearizedSystem)

Whether the linearized system's per-frequency assembly can be reduced to the
constant coefficients of a [`FrequencySweepPlan`](@ref).

It cannot when a component value depends on the symbolic frequency variable:
then the stored values themselves change with the frequency and there is no
constant quadratic to precompute.
"""
function cansweepondevice(lsys::HBLinearizedSystem)
    return isnothing(lsys.symfreqvar) && !lsys.symbolicvalues
end

# scatter the stored values of `As` into the slots of the system matrix its
# index map names, accumulating, which is what the per-frequency assembly
# does before applying the frequency factor
function scattercoefficient!(v::Vector{Complex{Float64}}, As, indexmap)
    nz = nonzeros(As)
    length(indexmap) == length(nz) || throw(DimensionMismatch(
        lazy"the index map has length $(length(indexmap)) but the matrix has $(length(nz)) stored entries."))
    @inbounds for j in eachindex(nz)
        v[indexmap[j]] += convert(Complex{Float64}, nz[j])
    end
    return v
end

# the column of each stored entry of a compressed sparse column matrix
function storedcolumns(A::SparseMatrixCSC)
    colof = Vector{Int}(undef, nnz(A))
    @inbounds for j in axes(A, 2)
        for q in nzrange(A, j)
            colof[q] = j
        end
    end
    return colof
end

"""
    planfrequencysweep(lsys::HBLinearizedSystem, backend;
        adjoint::Bool = false)

Build a [`FrequencySweepPlan`](@ref) for `lsys` on `backend`, together with the
compressed sparse row structure a device direct solver factorizes, as
`(plan, rowptr, colind)`.

With `adjoint` the structure and the coefficients describe the transpose of
the system matrix, whose solutions are the adjoint ones the noise, quantum
efficiency and sensitivity calculations need. That transpose is free
to form: compressed sparse row of the transpose is compressed sparse column of
the matrix, which is how the host holds it, so the adjoint plan is the same
coefficients in their original order against the original structure.

With cuDSS the adjoint costs a second symbolic analysis and a second
numeric factorization per batch, because cuDSS 0.8 has no transposed solve
(its `"solve_mode"` is documented in the header as unsupported), so the
factors of the forward system cannot be reused the way the host reuses
them with [`trysolvetranspose!`](@ref). The block path asks for the
adjoint plan for a different reason: its coefficients are in the matrix's
own stored order, which is the order a [`SparseBlockFactorization`](@ref)
fills its blocks from, and it costs no second factorization, since the
block factors solve both directions.
"""
function planfrequencysweep(lsys::HBLinearizedSystem, backend;
    adjoint::Bool = false)
    cansweepondevice(lsys) || throw(ArgumentError(
        "the linearized system's component values depend on the symbolic frequency variable, so its assembly is not a constant quadratic in the signal frequency."))
    A = lsys.Asparse
    n = size(A, 1)
    nz = nnz(A)

    # the four coefficients, in the stored order of A
    cst = copy(lsys.AoLjnmnzval)
    scattercoefficient!(cst, lsys.Amna0, lsys.Amna0indexmap)
    kinvL = zeros(Complex{Float64}, nz)
    scattercoefficient!(kinvL, lsys.invLnm, lsys.invLnmindexmap)
    kG = zeros(Complex{Float64}, nz)
    scattercoefficient!(kG, lsys.Gnm, lsys.Gnmindexmap)
    kC = zeros(Complex{Float64}, nz)
    scattercoefficient!(kC, lsys.Cnm, lsys.Cnmindexmap)

    # ... in the stored order of whichever structure is handed to the solver,
    # and with the matrix column of each stored entry, from which the kernel
    # takes its mode
    perm, rowptrhost, colindhost, colofhost = if adjoint
        collect(1:nz), SparseArrays.getcolptr(A), rowvals(A), storedcolumns(A)
    else
        p = cscvaluepermutation(A)
        At = sparse(transpose(A))
        p, SparseArrays.getcolptr(At), rowvals(At), rowvals(At)
    end
    d = x -> tobackend(backend, x[perm])
    colof = tobackend(backend, colofhost)
    assemble! = sweepassemblykernel!(backend, 64)
    plan = FrequencySweepPlan{typeof(colof),typeof(d(cst)),
        typeof(tobackend(backend, lsys.wpumpmodes)),typeof(assemble!),
        typeof(backend)}(
        colof, d(cst), d(kinvL), d(kG), d(kC),
        tobackend(backend, lsys.wpumpmodes), assemble!, backend,
        Int(lsys.Nmodes), nz)
    return plan, tobackend(backend, rowptrhost),
        tobackend(backend, colindhost)
end

"""
    assemblesweep!(nzval::AbstractMatrix, plan::FrequencySweepPlan,
        ws::AbstractVector)

Assemble the stored values of the linearized system matrix at each signal
frequency of `ws` into the corresponding column of `nzval`, which must have
one row per stored entry and one column per frequency.

Every stored value is written, so `nzval` need not be zeroed first.
"""
function assemblesweep!(nzval::AbstractMatrix, plan::FrequencySweepPlan,
    ws::AbstractVector)

    size(nzval, 1) == plan.nnz || throw(DimensionMismatch(
        lazy"`nzval` has $(size(nzval,1)) rows but the plan assembles $(plan.nnz) stored entries."))
    size(nzval, 2) == length(ws) || throw(DimensionMismatch(
        lazy"`nzval` has $(size(nzval,2)) columns but there are $(length(ws)) frequencies."))
    plan.assemble!(nzval, plan.colof, plan.cst, plan.kinvL, plan.kG, plan.kC,
        plan.wpump, ws, plan.nmodes, plan.nnz; ndrange = length(nzval))
    KernelAbstractions.synchronize(plan.backend)
    return nzval
end

"""
    portsolutionrows(nodeindices, portindices, Nmodes::Integer)

The rows of a solution of the linearized system which the scattering
parameter calculation reads: for each port, each of its two nodes which is
not ground, and each mode.

[`calcinputoutput!`](@ref) reads a solution only through
[`calcportvoltage`](@ref), which touches these rows and no others. On a
backend the solutions are produced there, so gathering these rows and
copying back only them replaces a transfer of the whole solution, which is
the same size as the state of the whole circuit, with one the size of the
scattering matrix.
"""
function portsolutionrows(nodeindices, portindices, Nmodes::Integer)
    rows = Int[]
    for p in portindices
        for t in 1:2
            key = nodeindices[t, p]
            key == 1 && continue        # ground carries no variable
            for j in 1:Nmodes
                push!(rows, (key - 2)*Nmodes + j)
            end
        end
    end
    sort!(rows)
    unique!(rows)
    return rows
end

# gather the named rows of every solution of every system in the batch
@kernel function gatherrowskernel!(out, @Const(X), @Const(rows), nrows, nrhs)
    gid = @index(Global)
    @inbounds begin
        q = gid - 1
        r = q % nrows + 1
        k = (q ÷ nrows) % nrhs + 1
        b = q ÷ (nrows*nrhs) + 1
        out[r, k, b] = X[Int(rows[r]), k, b]
    end
end

"""
    gatherportrows!(out, X, rows, backend)

Gather the rows named by `rows` from every right hand side of every system of
the batch `X` into `out`, with a kernel on `backend`. See
[`portsolutionrows`](@ref).
"""
function gatherportrows!(out::AbstractArray{<:Any,3}, X::AbstractArray{<:Any,3},
    rows::AbstractVector, backend)

    nrows, nrhs, nb = size(out)
    nrows == length(rows) || throw(DimensionMismatch(
        lazy"`out` has $(nrows) rows but $(length(rows)) were named."))
    (size(X, 2), size(X, 3)) == (nrhs, nb) || throw(DimensionMismatch(
        "`out` and `X` must agree in the number of right hand sides and the batch size."))
    gatherrowskernel!(backend, 64)(out, X, rows, nrows, nrhs;
        ndrange = length(out))
    KernelAbstractions.synchronize(backend)
    return out
end

"""
    needsadjointsolve(arrays::LinearizedArrays,
        noiseportimpedanceindices, noiseplan = nothing)

Whether the transposed (adjoint) linearized system must be solved at each
signal frequency: for the scattering parameter sensitivities always, and
otherwise when a consumer of the adjoint solution (the noise scattering
parameters, the quantum efficiency, the commutation relations, or the adjoint
node outputs) is requested together with a source of it. The dissipative
scattering blocks of a [`ScatteringNoisePlan`](@ref) are such a source, as
the lumped noise ports are.

This does not depend on the frequency, and both the host loop and the device
sweep test it, the latter to decide whether to allocate and solve the
adjoint direction at all (with cuDSS a whole second factorization, with a
block factorization a second solve against the same factors).
"""
function needsadjointsolve(arrays::LinearizedArrays,
    noiseportimpedanceindices, noiseplan = nothing)
    isempty(arrays.Ssensitivity) || return true
    hassource = !isempty(noiseportimpedanceindices) ||
        !isnothing(noiseplan) ||
        !isempty(arrays.nodefluxadjoint) || !isempty(arrays.voltageadjoint)
    hasconsumer = !isempty(arrays.Snoise) || !isempty(arrays.QE) ||
        !isempty(arrays.CM) || !isempty(arrays.Cnoise) ||
        !isempty(arrays.nodefluxadjoint) || !isempty(arrays.voltageadjoint)
    return hassource && hasconsumer
end

"""
    DeviceSweep

The state of a linearized sweep computed on a device a batch of
frequencies at a time, built by [`devicesolutions`](@ref) and driven by
four verbs: [`solvebatch!`](@ref) solves the batch beginning at a
frequency index and stages its solutions on the host,
[`forwardsolution!`](@ref) and [`adjointsolution!`](@ref) fill the
solution of any frequency of the staged batch, and
[`adjointdevice`](@ref) hands out the adjoint solution of a frequency
where it was computed, for the noise scattering parameters.

# Fields
- `backend`, `nb`, `F`: the backend, the batch size and the number of
    frequencies.
- `w`, `wpumpmodes`, `wshost`, `wsdev`: the signal frequencies, the pump
    mode offsets, and the frequencies of the current batch on the host and
    on the device.
- `fwd`, `adj`: the batch of each direction, `(plan, rowptr, colind,
    nzval, X, B)` on the sparse device factorization path or `(X,)` on
    the block path, and `nothing` when no adjoint was asked for.
- `fstage`, `astage`: the staging of each direction, `(full, rows, rowsd,
    gathered, host)`.
- `blocks`: the block path's `(plan, nzval, B, F, X, Xadj)`, or `nothing`.
- `scatstamps`, `scatstampsadjoint`, `providers`: the scattering block
    stamps of each direction and the providers they are evaluated through.
- `sweeps`: the cuDSS sweep of each direction, made on the first batch.
- `batchlo`: the first signal frequency of the batch currently staged.
"""
struct DeviceSweep{TB,TFw,TAd,TFs,TAs,TBl,TSt,TSa,TPr}
    backend::TB
    nb::Int
    F::Int
    n::Int
    nrhs::Int
    w::Vector{Float64}
    wpumpmodes::Vector{Float64}
    wshost::Vector{Float64}
    wsdev::Any
    fwd::TFw
    adj::TAd
    fstage::TFs
    astage::TAs
    blocks::TBl
    scatstamps::TSt
    scatstampsadjoint::TSa
    providers::TPr
    sweeps::Vector{Any}
    batchlo::Base.RefValue{Int}
end

"""
    devicesolutions(lsys::HBLinearizedSystem, bnm, w, backend, forward,
        adjoint = nothing; factorization = nothing, refine = true)

The [`DeviceSweep`](@ref) which computes the solutions of the linearized
system on `backend` a batch of frequencies at a time.
`factorization` is the linearized solve's factorization: a
[`BlockFactorization`](@ref) takes the batched block path below, anything
else the cuDSS one. `refine` asks single precision block factors to refine
against the double residual; `false` is the fully single precision sweep.

The sweep is driven by [`solvebatch!`](@ref), which solves the batch of
frequencies beginning at an index and stages its solutions on the host;
[`forwardsolution!`](@ref) and [`adjointsolution!`](@ref) then fill the
solution of any frequency of that batch, and [`adjointdevice`](@ref) hands
out the adjoint solution of a frequency where it was computed, for the noise
scattering parameters, which are formed there rather than brought back (see
[`devicenoise`](@ref)). The `providers` field is the
[`DeviceProviders`](@ref) the scattering blocks are evaluated through, or
`nothing` when there are none or they cannot be evaluated on `backend`; the
noise channels of the dissipative blocks are read through the same one.

The split is what keeps the host work parallel. `solvebatch!` is the only verb
which touches the device; once it has returned, its batch's frequencies can be
post-processed by as many workers as the host path uses, because reading a
solution touches staged host memory and nothing else.

The systems of a batch share one sparsity pattern, so cuDSS analyzes it once
and then refactorizes and solves the whole batch together, from the values of
one [`assemblesweep!`](@ref). With a sparse device factorization the batch
size is capped by [`uniformbatchlimit`](@ref); with a
[`BlockFactorization`](@ref) the cap does not apply and the batch is sized by
[`blocksystembytes`](@ref), the factors, the originals when refining, the
solutions of both directions and the value matrix of one system, against
half the backend's free memory and the length of the sweep.

`forward` and `adjoint` each describe what a direction needs, as a named tuple
`(full, rows)`. With `full` the whole solution is copied back, which the node
flux, voltage and sensitivity outputs need; otherwise only `rows` are gathered
and returned, which for the scattering parameters is a handful of port rows
out of the whole circuit (see [`portsolutionrows`](@ref)).

With a sparse device factorization the adjoint direction is a second uniform
batch over the transposed system: cuDSS 0.8 cannot solve against the
transpose of a factorization, so unlike the host, which gets its adjoint
solutions from the forward factors with [`trysolvetranspose!`](@ref), it
pays a second analysis and a second factorization per batch. A
[`BlockFactorization`](@ref) solves both directions from the one
factorization, reading the same factors the other way round
([`blocksolve!`](@ref) with `transposed = true`), and assembles one batch
of values in the stored column order for both.
"""
function devicesolutions(lsys::HBLinearizedSystem, bnm, w, backend, forward,
    adjoint = nothing; factorization = nothing, refine::Bool = true)

    T = Complex{Float64}
    n = size(lsys.Asparse, 1)
    nrhs = size(bnm, 2)
    nzA = nnz(lsys.Asparse)
    F = length(w)
    # a block factorization solves both directions from one factorization
    # per frequency, filled from the values assembled in the stored order;
    # its batch is sized by the memory of the factors and the solutions of
    # one system, within half the device's free memory. The cuDSS batch cap
    # does not apply to it; the memory and the sweep do
    usesblocks = factorization isa BlockFactorization
    nb = if usesblocks
        noderows, adjn = blocknodegraph(lsys.Asparse, lsys.Nmodes)
        blocksym = clustersymbolic(noderows, adjn, klunodeorder(adjn);
            target = lsys.Nmodes)
        Tf = something(factorization.precision, Float64)
        # the solutions of each direction, the shared right-hand side, the
        # values, and when refining the residual and correction of a solve
        refining = refine && Tf === Float32
        persystem = blocksystembytes(Complex{Tf}, blocksym; refine = refining,
            TA = T) + ((isnothing(adjoint) ? 2 : 3) + (refining ? 2 : 0))*
            n*nrhs*sizeof(T) + nzA*sizeof(T)
        clamp(freememory(backend) ÷ 2 ÷ persystem, 1, F)
    else
        min(uniformbatchlimit(nrhs), F)
    end

    # the right hand sides do not depend on the frequency, and are the same for
    # both directions, but cuDSS wants one set per system of the batch
    bhost = convert(Matrix{T}, bnm)
    hasscattering = !isnothing(lsys.scattering)

    # One set of scattering values serves both directions: a contribution
    # does not depend on the direction the system is assembled in, only its
    # destination does, so the adjoint's view shares the values and differs
    # only in where they land.
    scatstamps = plandevicescattering(lsys.scattering,
        hasscattering ? sweepdestinations(lsys.Asparse,
            lsys.scattering.Aindex, usesblocks) : Int[],
        nzA, nb, backend, lsys.Nmodes)
    scatstampsadjoint = if hasscattering && !isnothing(adjoint) && !usesblocks
        transposedestinations(scatstamps, sweepdestinations(lsys.Asparse,
            lsys.scattering.Aindex, true), backend)
    else
        scatstamps
    end
    # when every block's data is tabulated or constant the values are
    # computed on the backend and the host does nothing per frequency; a
    # callable provider is an arbitrary Julia function, so those stay on the
    # host
    scatproviders = plandeviceproviders(lsys.scattering, nb, backend,
        lsys.wpumpmodes, isnothing(lsys.scattering) ? 1.0 :
            lsys.scattering.scale)

    fstage = devicestage(forward, n, nrhs, nb, backend)
    astage = isnothing(adjoint) ? nothing :
        devicestage(adjoint, n, nrhs, nb, backend)
    # the block path: the values of a batch in the stored (column) order,
    # one shared right-hand side, the solutions of each direction, and the
    # factorization built on the backend from the pattern
    blocks = if usesblocks
        plan, _, _ = planfrequencysweep(lsys, backend; adjoint = true)
        Fb = factorize(factorization, lsys.Asparse;
            blocksize = lsys.Nmodes, backend = backend, nb = nb,
            refine = refine ? 6 : 0)
        (plan = plan,
            nzval = KernelAbstractions.allocate(backend, T, nzA, nb),
            B = tobackend(backend, bhost), F = Fb,
            X = KernelAbstractions.allocate(backend, T, n, nrhs, nb),
            Xadj = isnothing(adjoint) ? nothing :
                KernelAbstractions.allocate(backend, T, n, nrhs, nb))
    else
        nothing
    end
    fwd = usesblocks ? (X = blocks.X,) :
        devicebatch(lsys, backend, false, bhost, nb)
    adj = isnothing(adjoint) ? nothing :
        usesblocks ? (X = blocks.Xadj,) :
        devicebatch(lsys, backend, true, bhost, nb)

    wshost = zeros(Float64, nb)
    wsdev = tobackend(backend, wshost)

    return DeviceSweep(backend, nb, F, n, nrhs, collect(Float64, w),
        collect(Float64, lsys.wpumpmodes), wshost, wsdev, fwd, adj, fstage,
        astage, blocks, scatstamps, scatstampsadjoint, scatproviders,
        Any[nothing, nothing], Ref(0))
end

# the uniform batch of one direction on the sparse device factorization
# path: its sweep plan, pattern, values, solutions and right hand sides
function devicebatch(lsys, backend, isadjoint::Bool, bhost::Matrix, nb::Int)
    T = eltype(bhost)
    n = size(lsys.Asparse, 1)
    nrhs = size(bhost, 2)
    plan, rowptr, colind = planfrequencysweep(lsys, backend;
        adjoint = isadjoint)
    nzval = KernelAbstractions.allocate(backend, T, nnz(lsys.Asparse), nb)
    X = KernelAbstractions.allocate(backend, T, n, nrhs, nb)
    B = KernelAbstractions.allocate(backend, T, n, nrhs, nb)
    fill!(X, zero(T))
    bd = tobackend(backend, bhost)
    for k in 1:nb
        copyto!(view(B, :, :, k), bd)
    end
    return (plan = plan, rowptr = rowptr, colind = colind, nzval = nzval,
        X = X, B = B)
end

# The staging for one direction. A whole batch is brought back at once and
# read from host memory afterwards, so reading a solution does not touch
# the device and any number of workers may do it at once.
#
# Gathering only the named rows pays while they are a small part of the
# solution; past that the whole solution is staged instead. A circuit whose
# loss is spread along the line reaches this, since its noise ports touch
# almost every node.
function devicestage(spec, n::Int, nrhs::Int, nb::Int, backend)
    T = Complex{Float64}
    full = spec.full || 4*length(spec.rows) >= n
    rows = full ? Int[] : spec.rows
    return (full = full, rows = rows,
        rowsd = tobackend(backend, rows),
        gathered = KernelAbstractions.allocate(backend, T,
            max(length(rows), 1), nrhs, nb),
        host = Array{T}(undef, full ? n : length(rows), nrhs, nb))
end

# the frequencies of the batch beginning at `lo`, padded with the last real
# one so every system of a short final batch is well posed
function loadfrequencies!(ds::DeviceSweep, lo::Int)
    hi = min(lo + ds.nb - 1, ds.F)
    k = hi - lo + 1
    @inbounds for j in 1:ds.nb
        ds.wshost[j] = ds.w[j <= k ? lo + j - 1 : hi]
    end
    copyto!(ds.wsdev, ds.wshost)
    return hi
end

# the staging of a solved direction: nothing to bring back when the
# direction's solutions are read only where they were computed, as the
# noise scattering parameters read the adjoint ones
function stagebatch!(ds::DeviceSweep, b, stage)
    if stage.full
        copyto!(stage.host, b.X)
    elseif !isempty(stage.rows)
        gatherportrows!(stage.gathered, b.X, stage.rowsd, ds.backend)
        copyto!(stage.host, stage.gathered)
    end
    return nothing
end

# the block path: the batch of frequencies is factorized from its assembled
# values in one batched block LU and solved in both directions from the one
# factorization
function runblockbatch!(ds::DeviceSweep, stamps)
    bl = ds.blocks
    assemblesweep!(bl.nzval, bl.plan, ds.wsdev)
    isnothing(stamps) || applyscatteringstamps!(bl.nzval, stamps)
    fillandfactorize!(bl.F, bl.nzval)
    refinedsolve!(bl.X, bl.F, bl.B)
    isnothing(bl.Xadj) || refinedsolve!(bl.Xadj, bl.F, bl.B; transposed = true)
    stagebatch!(ds, ds.fwd, ds.fstage)
    isnothing(ds.adj) || stagebatch!(ds, ds.adj, ds.astage)
    return nothing
end

# the sparse device factorization path of one direction: the batch is
# assembled, analyzed once on the first batch, refactorized and solved
function runbatch!(ds::DeviceSweep, slot::Int, b, stage, stamps)
    assemblesweep!(b.nzval, b.plan, ds.wsdev)
    isnothing(stamps) || applyscatteringstamps!(b.nzval, stamps)
    if isnothing(ds.sweeps[slot])
        ds.sweeps[slot] = _cudss_sweep(b.rowptr, b.colind, b.nzval, b.X, b.B)
    end
    _cudss_sweepsolve!(ds.sweeps[slot])
    stagebatch!(ds, b, stage)
    return nothing
end

"""
    solvebatch!(ds::DeviceSweep, lo::Integer)

Solve the batch of frequencies beginning at index `lo` on the device and
stage its solutions on the host. This is the only verb which touches the
device; once it has returned, the batch's frequencies can be read by as
many workers as the host path uses through [`forwardsolution!`](@ref) and
[`adjointsolution!`](@ref), because reading a solution touches staged host
memory and nothing else.
"""
function solvebatch!(ds::DeviceSweep, lo::Integer)
    hi = loadfrequencies!(ds, Int(lo))
    k = hi - lo + 1
    if !isnothing(ds.scatstamps)
        if isnothing(ds.providers)
            stagescatteringstamps!(ds.scatstamps, ds.w, lo, k, ds.wpumpmodes)
        else
            stagedeviceproviders!(ds.scatstamps.values, ds.providers, ds.w,
                lo, k)
        end
    end
    if !isnothing(ds.blocks)
        runblockbatch!(ds, ds.scatstamps)
    else
        runbatch!(ds, 1, ds.fwd, ds.fstage, ds.scatstamps)
        isnothing(ds.adj) || runbatch!(ds, 2, ds.adj, ds.astage,
            ds.scatstampsadjoint)
    end
    ds.batchlo[] = lo
    return nothing
end

# the slot of frequency `i` in the staged batch
function batchslot(ds::DeviceSweep, i::Integer)
    s = i - ds.batchlo[] + 1
    (1 <= s <= ds.nb) || throw(ArgumentError(
        lazy"frequency $(i) is not in the batch beginning at $(ds.batchlo[])."))
    return s
end

function fillsolution!(phin, ds::DeviceSweep, stage, i::Integer)
    s = batchslot(ds, i)
    if stage.full
        copyto!(phin, view(stage.host, :, :, s))
    else
        # down the rows for each right hand side in turn, which reads the
        # staged block contiguously
        @inbounds for k in axes(phin, 2)
            for (r, row) in enumerate(stage.rows)
                phin[row, k] = stage.host[r, k, s]
            end
        end
    end
    return phin
end

"""
    forwardsolution!(phin, ds::DeviceSweep, i::Integer)
    adjointsolution!(phin, ds::DeviceSweep, i::Integer)

Fill `phin` with the forward or adjoint solution of frequency `i` of the
batch staged by the last [`solvebatch!`](@ref): the whole solution when the
direction was staged in full, otherwise only its gathered rows. The form
[`hblinsolve_inner!`](@ref) takes as `presolved` and `presolvedadjoint`.
"""
forwardsolution!(phin, ds::DeviceSweep, i::Integer) =
    fillsolution!(phin, ds, ds.fstage, i)

"""
    adjointsolution!(phin, ds::DeviceSweep, i::Integer)

The adjoint counterpart of [`forwardsolution!`](@ref).
"""
adjointsolution!(phin, ds::DeviceSweep, i::Integer) =
    fillsolution!(phin, ds, ds.astage, i)

"""
    adjointdevice(ds::DeviceSweep, i::Integer)

The adjoint solution of frequency `i` of the staged batch where it was
computed, a view into the device batch, for the noise scattering parameters
which are formed there rather than brought back (see
[`devicenoise`](@ref)). A read of the solved batch, so several workers may
hold different frequencies of it at once.
"""
adjointdevice(ds::DeviceSweep, i::Integer) =
    view(ds.adj.X, :, :, batchslot(ds, i))
