
"""
    FrequencySweepPlan

Everything needed to assemble the linearized system matrix of
[`hblinsolve`](@ref) at many signal frequencies at once, on a backend.

The per-frequency assembly of [`assemblesystemmatrix!`](@ref) is

    A = AoLjnm + invLnm + im*Gnm*w - Cnm*w^2 + Amna0 + im*AmnaG*w

with the frequency of each stored entry taken from its column's mode, and
with the stored value of every frequency dependent term conjugated where that
mode frequency is negative (see [`modevalue`](@ref)). Each stored entry is
therefore an independent quadratic in its own mode frequency, and the whole
assembly collapses to four constant coefficient vectors and one kernel:

    A[q] = cst[q] + sel(kinvL[q]) + im*wm*sel(kG[q]) - wm^2*sel(kC[q])

where `wm = ws + wpumpmodes[mode of q]` and `sel` conjugates when `wm < 0`.
The terms which share a frequency power and a conjugation rule are summed
into one coefficient at build time: `AoLjnm` with `Amna0`, and `Gnm` with
`AmnaG`. Conjugation distributes over that sum, so this is exact.

Because the coefficients do not depend on the signal frequency, one kernel
fills the stored values of a whole batch of frequencies, which is what cuDSS
wants for a uniform batch: the batch shares one sparsity pattern and one
symbolic analysis, and only the values differ.

# Fields
- `colof`: the matrix column of each stored entry, from which its mode
    follows, so the kernel needs no search. For the forward structure this is
    also the structure's own column index array; for the transposed one it is
    not, so it is carried separately.
- `cst`, `kinvL`, `kG`, `kC`: the four coefficient vectors, in the stored
    order of the transposed structure.
- `wpump`: the pump mode frequency offsets of the signal modes.
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
    isnothing(lsys.symfreqvar) || return false
    return isempty(lsys.Cnmfreqsubstindices) &&
        isempty(lsys.Gnmfreqsubstindices) &&
        isempty(lsys.invLnmfreqsubstindices) &&
        isempty(lsys.AmnaGfreqsubstindices)
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

The adjoint costs a second symbolic analysis and a second numeric
factorization per batch, because cuDSS 0.8 has no transposed solve (its
`"solve_mode"` is documented in the header as unsupported), so the factors of
the forward system cannot be reused the way the host reuses them with
[`trysolvetranspose!`](@ref).
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
    scattercoefficient!(kG, lsys.AmnaG, lsys.AmnaGindexmap)
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
sweep test it, the latter because it has no factorization to solve the
transposed system against.
"""
function needsadjointsolve(arrays::LinearizedArrays,
    noiseportimpedanceindices, noiseplan = nothing)
    isempty(arrays.Ssensitivity) || return true
    hassource = !isempty(noiseportimpedanceindices) ||
        !isnothing(noiseplan) ||
        !isempty(arrays.nodefluxadjoint) || !isempty(arrays.voltageadjoint)
    hasconsumer = !isempty(arrays.Snoise) || !isempty(arrays.QE) ||
        !isempty(arrays.CM) || !isempty(arrays.nodefluxadjoint) ||
        !isempty(arrays.voltageadjoint)
    return hassource && hasconsumer
end

"""
    devicesolutions(lsys::HBLinearizedSystem, bnm, w, backend, forward,
        adjoint = nothing)

Callbacks which fill the solutions of the linearized system at a signal
frequency, computing them on `backend` a batch of frequencies at a time.

Returns
`(batchsize, solvebatch!, providers, forward, adjoint, adjointdevice)`.
`solvebatch!(lo)` solves the batch of frequencies beginning at index `lo` and
stages its solutions on the host; `forward` and `adjoint` then fill the
solution of any frequency of that batch, in the form
[`hblinsolve_inner!`](@ref) takes as `presolved` and `presolvedadjoint`, with
`adjoint` `nothing` when none was asked for. `adjointdevice` hands out the
adjoint solution of a frequency where it was computed, for the noise
scattering parameters, which are formed there rather than brought back (see
[`devicenoise`](@ref)). `providers` is the [`DeviceProviders`](@ref) the
scattering blocks are evaluated through, or `nothing` when there are none or
they cannot be evaluated on `backend`; the noise channels of the dissipative
blocks are read through the same one.

The split is what keeps the host work parallel. `solvebatch!` is the only part
which touches the device; once it has returned, its batch's frequencies can be
post-processed by as many workers as the host path uses, because reading a
solution touches staged host memory and nothing else.

The systems of a batch share one sparsity pattern, so cuDSS analyzes it once
and then refactorizes and solves the whole batch together, from the values of
one [`assemblesweep!`](@ref). The batch size is capped by
[`uniformbatchlimit`](@ref).

`forward` and `adjoint` each describe what a direction needs, as a named tuple
`(full, rows)`. With `full` the whole solution is copied back, which the node
flux, voltage and sensitivity outputs need; otherwise only `rows` are gathered
and returned, which for the scattering parameters is a handful of port rows
out of the whole circuit (see [`portsolutionrows`](@ref)).

The adjoint direction is a second uniform batch, over the transposed system.
cuDSS 0.8 cannot solve against the transpose of a factorization, so unlike the
host, which gets its adjoint solutions from the forward factors with
[`trysolvetranspose!`](@ref), this pays a second analysis and a second
factorization per batch.
"""
function devicesolutions(lsys::HBLinearizedSystem, bnm, w, backend, forward,
    adjoint = nothing)

    T = Complex{Float64}
    n = size(lsys.Asparse, 1)
    nrhs = size(bnm, 2)
    nzA = nnz(lsys.Asparse)
    F = length(w)
    nb = min(uniformbatchlimit(nrhs), F)

    # the right hand sides do not depend on the frequency, and are the same for
    # both directions, but cuDSS wants one set per system of the batch
    bhost = convert(Matrix{T}, bnm)
    hasscattering = !isnothing(lsys.scattering)
    function newbatch(isadjoint)
        plan, rowptr, colind = planfrequencysweep(lsys, backend;
            adjoint = isadjoint)
        nzval = KernelAbstractions.allocate(backend, T, nzA, nb)
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

    # One set of scattering values serves both directions: a contribution
    # does not depend on the direction the system is assembled in, only its
    # destination does, so the adjoint's view shares the values and differs
    # only in where they land.
    scatstamps = plandevicescattering(lsys.scattering,
        hasscattering ? sweepdestinations(lsys.Asparse,
            lsys.scattering.Aindex, false) : Int[],
        nzA, nb, backend, lsys.Nmodes)
    scatstampsadjoint = if hasscattering && !isnothing(adjoint)
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

    # The staging for one direction. A whole batch is brought back at once
    # and read from host memory afterwards, so reading a solution does not
    # touch the device and any number of workers may do it at once.
    #
    # Gathering only the named rows pays while they are a small part of the
    # solution; past that the whole solution is staged instead. A circuit
    # whose loss is spread along the line reaches this, since its noise
    # ports touch almost every node.
    function newstage(spec)
        full = spec.full || 4*length(spec.rows) >= n
        rows = full ? Int[] : spec.rows
        return (full = full, rows = rows,
            rowsd = tobackend(backend, rows),
            gathered = KernelAbstractions.allocate(backend, T,
                max(length(rows), 1), nrhs, nb),
            host = Array{T}(undef, full ? n : length(rows), nrhs, nb))
    end

    fwd = newbatch(false)
    fstage = newstage(forward)
    adj = isnothing(adjoint) ? nothing : newbatch(true)
    astage = isnothing(adjoint) ? nothing : newstage(adjoint)

    wshost = zeros(Float64, nb)
    wsdev = tobackend(backend, wshost)
    sweeps = Any[nothing, nothing]

    # the frequencies of the batch beginning at `lo`, padded with the last real
    # one so every system of a short final batch is well posed
    function loadfrequencies!(lo)
        hi = min(lo + nb - 1, F)
        k = hi - lo + 1
        @inbounds for j in 1:nb
            wshost[j] = w[j <= k ? lo + j - 1 : hi]
        end
        copyto!(wsdev, wshost)
        return hi
    end

    function runbatch!(slot, b, stage, stamps)
        assemblesweep!(b.nzval, b.plan, wsdev)
        isnothing(stamps) || applyscatteringstamps!(b.nzval, stamps)
        if isnothing(sweeps[slot])
            sweeps[slot] = _cudss_sweep(b.rowptr, b.colind, b.nzval, b.X, b.B)
        end
        _cudss_sweepsolve!(sweeps[slot])
        # nothing to bring back when the direction's solutions are read only
        # where they were computed, as the noise scattering parameters read
        # the adjoint ones
        if stage.full
            copyto!(stage.host, b.X)
        elseif !isempty(stage.rows)
            gatherportrows!(stage.gathered, b.X, stage.rowsd, backend)
            copyto!(stage.host, stage.gathered)
        end
        return nothing
    end

    # the first signal frequency of the batch currently staged
    batchlo = Ref(0)
    solvebatch! = function(lo)
        hi = loadfrequencies!(lo)
        k = hi - lo + 1
        if !isnothing(scatstamps)
            if isnothing(scatproviders)
                stagescatteringstamps!(scatstamps, w, lo, k, lsys.wpumpmodes)
            else
                stagedeviceproviders!(scatstamps.values, scatproviders, w,
                    lo, k)
            end
        end
        runbatch!(1, fwd, fstage, scatstamps)
        isnothing(adj) || runbatch!(2, adj, astage, scatstampsadjoint)
        batchlo[] = lo
        return nothing
    end

    function fillsolution!(phin, stage, i)
        s = i - batchlo[] + 1
        (1 <= s <= nb) || throw(ArgumentError(
            lazy"frequency $(i) is not in the batch beginning at $(batchlo[])."))
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

    # `adjointdevice` is a read of the solved batch, so several workers may
    # hold different frequencies of it at once.
    return (batchsize = nb,
        solvebatch! = solvebatch!,
        providers = scatproviders,
        forward = (i, phin) -> fillsolution!(phin, fstage, i),
        adjoint = isnothing(adj) ? nothing :
            (i, phin) -> fillsolution!(phin, astage, i),
        adjointdevice = isnothing(adj) ? nothing :
            i -> view(adj.X, :, :, i - batchlo[] + 1))
end

# ---------------------------------------------------------------------------
# the noise scattering parameters on the backend
# ---------------------------------------------------------------------------

"""
    noiseoutputwavekernel!

The output power waves at the noise ports, from the adjoint solution, one work
item per (noise port, mode, right hand side).

This is [`calcinputoutputnoise!`](@ref) restricted to its output waves, which
is all of it that reads a solution: the input waves come from the source terms
alone and are computed on the host.
"""
@kernel function noiseoutputwavekernel!(out, @Const(phin), @Const(node1),
        @Const(node2), @Const(values), @Const(codes), @Const(wmodes),
        Nmodes, Nnoise, nrhs)
    gid = @index(Global)
    @inbounds begin
        q = gid - 1
        j = q % Nmodes + 1
        i = (q ÷ Nmodes) % Nnoise + 1
        k = q ÷ (Nmodes*Nnoise) + 1
        w = wmodes[j]
        z = impedance(values[i], codes[i], w)
        kval = portwavescale(z, w)
        k1 = Int(node1[i]); k2 = Int(node2[i])
        v = if k1 == 1
            -phin[(k2-2)*Nmodes+j, k]
        elseif k2 == 1
            phin[(k1-2)*Nmodes+j, k]
        else
            phin[(k1-2)*Nmodes+j, k] - phin[(k2-2)*Nmodes+j, k]
        end
        v *= im*w
        # no source current at a noise port, so the whole port current is the
        # one the port voltage drives through the port impedance
        current = -v/z
        out[(i-1)*Nmodes+j, k] = (kval*(v - conj(z)*current))/2
    end
end

# the sign of the mode frequency of each row of the noise scattering matrix,
# which is what the commutation relations weight that row by
@kernel function modesignkernel!(signs, @Const(wmodes), Nmodes)
    gid = @index(Global)
    @inbounds signs[gid] = sign(wmodes[(gid - 1) % Nmodes + 1])
end


"""
    DeviceNoisePlan

The noise ports of a circuit, on a backend, for
[`noiseoutputwavekernel!`](@ref).
"""
struct DeviceNoisePlan{VI,VC,VZ,B}
    node1::VI
    node2::VI
    values::VC
    codes::VZ
    nmodes::Int
    nnoise::Int
    backend::B
end

"""
    plandevicenoise(nodeindices, componenttypes, noiseportimpedanceindices,
        noiseportimpedances, Nmodes, backend)

Build a [`DeviceNoisePlan`](@ref) for the noise ports of a circuit.
"""
function plandevicenoise(nodeindices, componenttypes,
    noiseportimpedanceindices, noiseportimpedances, Nmodes::Integer, backend)

    nn = length(noiseportimpedanceindices)
    node1 = Int[nodeindices[1, p] for p in noiseportimpedanceindices]
    node2 = Int[nodeindices[2, p] for p in noiseportimpedanceindices]
    codes = Int32[impedancecode(componenttypes[p])
        for p in noiseportimpedanceindices]
    values = Complex{Float64}[convert(Complex{Float64}, z)
        for z in noiseportimpedances]
    return DeviceNoisePlan(tobackend(backend, node1),
        tobackend(backend, node2), tobackend(backend, values),
        tobackend(backend, codes), Int(Nmodes), nn, backend)
end

"""
    devicenoise(plan::DeviceNoisePlan, blockplan, providers,
        adjointsolution, nrhs, wpumpmodes, w, keepmatrix)

A callback which computes the noise scattering parameters of a signal
frequency from the adjoint solutions on the backend.

`temperatures` is one temperature per noise channel, or `nothing`; the waves
of a warm channel are scaled where they are computed, so nothing downstream
of this knows about temperature.

`blockplan` is the [`DeviceBlockNoisePlan`](@ref) of the dissipative
scattering blocks, or `nothing` when there are none; their channels follow
the noise ports of the lumped components in the rows, as they do on the
host, and are evaluated through the same [`DeviceProviders`](@ref)
`providers` the stamps are.

Returned in the form [`hblinsolve_inner!`](@ref) takes as `presolvednoise`:
called with the frequency index, its input waves and the destination for the
noise scattering matrix, it returns the [`NoiseReduction`](@ref) the quantum
efficiency and the commutation relations read.

The adjoint solution never leaves the backend, and neither does the noise
scattering matrix unless `keepmatrix`: on a line with loss spread along it
that matrix has a row per noise port mode and is the largest thing in the
sweep, while what is read of it is one number per port mode.

`Snoise = noiseoutputwave/inputwave` is formed as a product with the inverse
of the input waves, which is a dense matrix the size of the scattering matrix
and so is inverted on the host.
"""
function devicenoise(plan::DeviceNoisePlan, blockplan, providers,
    adjointsolution, nrhs::Integer, wpumpmodes, w, keepmatrix::Bool,
    temperatures = nothing)

    backend = plan.backend
    T = Complex{Float64}
    lumpedrows = plan.nnoise*plan.nmodes
    nrows = lumpedrows +
        (isnothing(blockplan) ? 0 : blockplan.nchannels*plan.nmodes)
    out = KernelAbstractions.allocate(backend, T, nrows, nrhs)
    Snoise = KernelAbstractions.allocate(backend, T, nrows, nrhs)
    invinput = KernelAbstractions.allocate(backend, T, nrhs, nrhs)
    # the two reductions go through the backend's own `sum!`, a tree over
    # the noise index; summing each column in one work item would be a
    # dozen threads each walking every noise port in turn
    absq = KernelAbstractions.allocate(backend, Float64, nrows, nrhs)
    signs = KernelAbstractions.allocate(backend, Float64, nrows, 1)
    denomd = KernelAbstractions.allocate(backend, Float64, 1, nrhs)
    signedd = KernelAbstractions.allocate(backend, Float64, 1, nrhs)
    wmodesd = KernelAbstractions.allocate(backend, Float64, plan.nmodes)
    denom = zeros(Float64, nrhs)
    signed = zeros(Float64, nrhs)
    wmodes = zeros(Float64, plan.nmodes)
    reduction = NoiseReduction(denom, signed)
    # The occupation of each channel mode, when any channel is warm. It
    # depends on the mode frequencies, so it is rebuilt per signal frequency
    # and sent; one entry per channel mode, which is small beside the
    # waves. It enters the sum the quantum efficiency reads and not the one
    # the commutation relations read, which is why the latter stay at one.
    warm = !isnothing(temperatures) && !all(iszero, temperatures)
    occupationhost = warm ? zeros(Float64, nrows) : Float64[]
    occupationd = warm ? KernelAbstractions.allocate(backend, Float64, nrows) :
        KernelAbstractions.allocate(backend, Float64, 0)
    # the added noise covariance, when it is asked for, is formed here rather
    # than by bringing the noise scattering matrix home: it is one product of
    # the matrix with itself and comes back as a port mode square.
    Cwork = KernelAbstractions.allocate(backend, T, nrows, nrhs)
    Cdev = KernelAbstractions.allocate(backend, T, nrhs, nrhs)
    Chost = Matrix{T}(undef, nrhs, nrhs)

    return function(i, inputwave, Snoiseview, Cnoiseview = nothing)
        @inbounds for m in eachindex(wmodes)
            wmodes[m] = w[i] + wpumpmodes[m]
        end
        copyto!(wmodesd, wmodes)
        if lumpedrows > 0
            noiseoutputwavekernel!(backend, 64)(out, adjointsolution(i),
                plan.node1, plan.node2, plan.values, plan.codes, wmodesd,
                plan.nmodes, plan.nnoise, nrhs;
                ndrange = lumpedrows*nrhs)
        end
        if !isnothing(blockplan)
            # the factor of each block's covariance at each mode frequency,
            # then the contraction of the adjoint solution against it
            if isnothing(providers.funcs)
                blocknoisefactorkernel!(backend, 64)(blockplan.factors,
                    blockplan.blockindex, blockplan.factoroff,
                    providers.nports, providers.freqoff, providers.nfreq,
                    providers.freqs, providers.valoff, providers.vals,
                    providers.conjsym, providers.extrapcode, wmodesd,
                    plan.nmodes, blockplan.nentries;
                    ndrange = blockplan.nentries*plan.nmodes)
            else
                blocknoiseentryfactorkernel!(backend, 64)(blockplan.factors,
                    blockplan.blockindex, blockplan.factoroff,
                    providers.nports, providers.funcs, providers.conjsym,
                    wmodesd, plan.nmodes, blockplan.nentries;
                    ndrange = blockplan.nentries*plan.nmodes)
            end
            KernelAbstractions.synchronize(backend)
            blocknoisecontractkernel!(backend, 64)(out, adjointsolution(i),
                blockplan.factors, blockplan.blockindex, blockplan.factoroff,
                blockplan.auxbase, providers.nports, blockplan.channelentry,
                blockplan.channellocal, wmodesd, plan.nmodes,
                blockplan.nchannels, lumpedrows, nrhs;
                ndrange = blockplan.nchannels*plan.nmodes*nrhs)
        end
        KernelAbstractions.synchronize(backend)
        copyto!(invinput, inv(Matrix{T}(inputwave)))
        mul!(Snoise, out, invinput)
        modesignkernel!(backend, 64)(signs, wmodesd, plan.nmodes;
            ndrange = nrows)
        KernelAbstractions.synchronize(backend)
        # the same two passes and two reductions as at zero temperature, with
        # the occupation folded into the one the quantum efficiency reads
        if warm
            noiseoccupation!(occupationhost, temperatures, wmodes,
                plan.nmodes)
            copyto!(occupationd, occupationhost)
            absq .= abs2.(Snoise) .* occupationd
        else
            absq .= abs2.(Snoise)
        end
        sum!(denomd, absq)
        absq .= abs2.(Snoise) .* signs
        sum!(signedd, absq)
        KernelAbstractions.synchronize(backend)
        copyto!(denom, vec(Array(denomd)))
        copyto!(signed, vec(Array(signedd)))
        keepmatrix && isempty(Snoiseview) == false && copyto!(Snoiseview,
            Array(Snoise))
        if !isnothing(Cnoiseview)
            # Cnoise[i,j] = sum_c occupation[c] Snoise[c,i] conj(Snoise[c,j])
            if warm
                Cwork .= occupationd .* conj.(Snoise)
            else
                Cwork .= conj.(Snoise)
            end
            mul!(Cdev, transpose(Snoise), Cwork)
            KernelAbstractions.synchronize(backend)
            copyto!(Chost, Cdev)
            copyto!(Cnoiseview, Chost)
        end
        return reduction
    end
end

"""
    blocknoisefactorkernel!

The factor `L` of the vacuum noise covariance `I - S S'` of each dissipative
scattering block at each mode frequency, one work item per (block, mode),
from tabulated or constant scattering data.

The factorization is [`psdcholesky!`](@ref), which the host path runs too, so
the two agree on the channels themselves and not only on the sums over them.
Its covariance is built from `2 n^3` interpolations rather than caching the
`n^2` scattering parameters, because a work item has nowhere to cache them
and a block has few ports.
"""
@kernel function blocknoisefactorkernel!(L, @Const(blockindex),
        @Const(factoroff), @Const(nports), @Const(freqoff), @Const(nfreq),
        @Const(freqs), @Const(valoff), @Const(vals), @Const(conjsym),
        @Const(extrapcode), @Const(wmodes), Nmodes, nentries)
    gid = @index(Global)
    @inbounds begin
        g = gid - 1
        m = g % Nmodes + 1
        e = g ÷ Nmodes + 1
        bi = Int(blockindex[e])
        n = Int(nports[bi])
        off = Int(factoroff[e]) + (m-1)*n*n
        T = eltype(L)
        wm = wmodes[m]
        if iszero(wm)
            # the wave normalization is singular at zero frequency, where the
            # lumped noise ports are zero too
            for j in 1:n*n
                L[off + j] = zero(T)
            end
        else
            isconj = conjsym[bi] != 0
            wq = isconj ? abs(wm) : wm
            neg = isconj && wm < 0
            fo = Int(freqoff[bi]); nf = Int(nfreq[bi]); vo = Int(valoff[bi])
            ec = extrapcode[bi]
            for c in 1:n
                for p in c:n
                    acc = p == c ? one(T) : zero(T)
                    for l in 1:n
                        spl = tableentry(freqs, vals, fo, nf, vo, n, p, l,
                            wq, ec)
                        scl = tableentry(freqs, vals, fo, nf, vo, n, c, l,
                            wq, ec)
                        if neg
                            spl = conj(spl); scl = conj(scl)
                        end
                        acc -= spl*conj(scl)
                    end
                    L[off + (c-1)*n + p] = acc
                end
            end
            psdcholesky!(L, off, n)
        end
    end
end

# the same, for blocks whose scattering parameters come from a callable of the
# `:entry` form
@kernel function blocknoiseentryfactorkernel!(L, @Const(blockindex),
        @Const(factoroff), @Const(nports), @Const(funcs), @Const(conjsym),
        @Const(wmodes), Nmodes, nentries)
    gid = @index(Global)
    @inbounds begin
        g = gid - 1
        m = g % Nmodes + 1
        e = g ÷ Nmodes + 1
        bi = Int(blockindex[e])
        n = Int(nports[bi])
        off = Int(factoroff[e]) + (m-1)*n*n
        T = eltype(L)
        wm = wmodes[m]
        if iszero(wm)
            for j in 1:n*n
                L[off + j] = zero(T)
            end
        else
            isconj = conjsym[bi] != 0
            wq = isconj ? abs(wm) : wm
            neg = isconj && wm < 0
            f = funcs[bi]
            for c in 1:n
                for p in c:n
                    acc = p == c ? one(T) : zero(T)
                    for l in 1:n
                        spl = T(f(p, l, wq))
                        scl = T(f(c, l, wq))
                        if neg
                            spl = conj(spl); scl = conj(scl)
                        end
                        acc -= spl*conj(scl)
                    end
                    L[off + (c-1)*n + p] = acc
                end
            end
            psdcholesky!(L, off, n)
        end
    end
end

"""
    blocknoisecontractkernel!

The noise output waves of the scattering block channels, one work item per
(channel, mode, right hand side).

The noise a block adds is a source in its auxiliary port current rows, so by
the adjoint identity its contribution at an output port is that source
contracted against those same rows of the adjoint solution:
`sqrt(abs(w)) sum_p L[p,c] i[p]`. Nothing else of the solution is read, which
is why the adjoint solutions of a circuit with dissipative blocks need not
come back from the backend at all.
"""
@kernel function blocknoisecontractkernel!(out, @Const(phin), @Const(L),
        @Const(blockindex), @Const(factoroff), @Const(auxbase),
        @Const(nports), @Const(channelentry), @Const(channellocal),
        @Const(wmodes), Nmodes, nchannels, rowoffset, nrhs)
    gid = @index(Global)
    @inbounds begin
        g = gid - 1
        m = g % Nmodes + 1
        ch = (g ÷ Nmodes) % nchannels + 1
        k = g ÷ (Nmodes*nchannels) + 1
        e = Int(channelentry[ch]); c = Int(channellocal[ch])
        bi = Int(blockindex[e])
        n = Int(nports[bi])
        off = Int(factoroff[e]) + (m-1)*n*n
        ab = Int(auxbase[e])
        acc = zero(eltype(out))
        for p in 1:n
            acc += L[off + (c-1)*n + p]*phin[ab + (p-1)*Nmodes + m, k]
        end
        out[rowoffset + (ch-1)*Nmodes + m, k] = sqrt(abs(wmodes[m]))*acc
    end
end

"""
    DeviceBlockNoisePlan

The vacuum noise channels of the dissipative scattering blocks, on a backend.

Where the host path reads the auxiliary port current rows of an adjoint
solution it has brought back, this describes the same channels as flat
tables, so [`blocknoisefactorkernel!`](@ref) and
[`blocknoisecontractkernel!`](@ref) can form them where the solution already
is. The scattering data itself comes from the [`DeviceProviders`](@ref) the
stamps are evaluated through, so this holds only what those do not: which
block each entry is, where its auxiliary rows and its factor live, and which
entry each channel belongs to.
"""
struct DeviceBlockNoisePlan{VI,VC,B}
    blockindex::VI
    factoroff::VI
    auxbase::VI
    channelentry::VI
    channellocal::VI
    factors::VC
    nentries::Int
    nchannels::Int
    nmodes::Int
    backend::B
end

"""
    withfactors(bp::DeviceBlockNoisePlan)

A copy of the plan with scratch of its own for the covariance factors,
sharing everything else, so that several workers can form the noise channels
of different signal frequencies at once. The factors are written by
[`blocknoisefactorkernel!`](@ref) and read by
[`blocknoisecontractkernel!`](@ref) at one frequency, so they cannot be
shared; the tables which say where each block is can be.
"""
function withfactors(bp::DeviceBlockNoisePlan)
    return DeviceBlockNoisePlan(bp.blockindex, bp.factoroff, bp.auxbase,
        bp.channelentry, bp.channellocal, similar(bp.factors), bp.nentries,
        bp.nchannels, bp.nmodes, bp.backend)
end

"""
    plandeviceblocknoise(ssys, noiseplan, Nmodes, backend)

Build a [`DeviceBlockNoisePlan`](@ref) from a
[`ScatteringNoisePlan`](@ref), or `nothing` when there is none.
"""
plandeviceblocknoise(ssys, ::Nothing, Nmodes, backend) = nothing
function plandeviceblocknoise(ssys, noiseplan::ScatteringNoisePlan,
    Nmodes::Integer, backend)

    ne = length(noiseplan.blockindices)
    blockindex = Int32[bi for bi in noiseplan.blockindices]
    auxbase = Int32[ssys.blocks[bi].auxbase for bi in noiseplan.blockindices]
    factoroff = Vector{Int32}(undef, ne)
    channelentry = Vector{Int32}(undef, noiseplan.Nchannels)
    channellocal = Vector{Int32}(undef, noiseplan.Nchannels)
    off = 0
    for (e, bi) in enumerate(noiseplan.blockindices)
        n = ssys.blocks[bi].block.nports
        factoroff[e] = off
        off += n*n*Nmodes
        for c in 1:n
            channelentry[noiseplan.channelbase[e] + c] = e
            channellocal[noiseplan.channelbase[e] + c] = c
        end
    end
    factors = KernelAbstractions.allocate(backend, Complex{Float64}, off)
    return DeviceBlockNoisePlan(tobackend(backend, blockindex),
        tobackend(backend, factoroff), tobackend(backend, auxbase),
        tobackend(backend, channelentry), tobackend(backend, channellocal),
        factors, ne, noiseplan.Nchannels, Int(Nmodes), backend)
end

# ---------------------------------------------------------------------------
# the scattering block stamps on a backend
# ---------------------------------------------------------------------------

# add one value per (contribution, frequency of the batch) into the stored
# values of that frequency's system matrix
@kernel function scatteringstampkernel!(nzval, @Const(index), @Const(values),
        ncontrib, nzA)
    gid = @index(Global)
    @inbounds begin
        q = gid - 1
        c = q % ncontrib + 1
        f = q ÷ ncontrib + 1
        nzval[Int(index[c]) + (f-1)*nzA] += values[c, f]
    end
end

"""
    sweepdestinations(A::SparseMatrixCSC, Aindex, adjoint::Bool)

The destination of each scattering contribution in the stored order a
[`FrequencySweepPlan`](@ref) uses, given its destinations in the stored order
of `A`.

A plan built for the transpose stores its values in `A`'s own order, so the
indices pass through; one built for the matrix itself stores them permuted
into compressed sparse row order, so they must be permuted the same way.
"""
function sweepdestinations(A::SparseMatrixCSC, Aindex::AbstractVector,
    adjoint::Bool)
    adjoint && return collect(Int, Aindex)
    # cscvaluepermutation carries a compressed sparse row slot to the stored
    # entry it holds; this is the other direction
    tocsr = invperm(cscvaluepermutation(A))
    return Int[tocsr[j] for j in Aindex]
end

"""
    DeviceScatteringStamps

The scattering block contribution to a batch of system matrices on a backend.

[`assemblesweep!`](@ref) covers the terms which are a constant quadratic in
the signal frequency, which the scattering blocks are not: each is evaluated
through its own provider, an arbitrary callable or an interpolation of
tabulated data, so its values have to be computed on the host. What is left
is a gather-add into the stored entries, which is this.

The values serve both directions of a sweep: a block's contribution does not
depend on the direction the system is assembled in, only on where in that
system's stored order each scalar lands. So they are computed and sent once,
and the adjoint direction gets a second view of them with its own
destinations (see [`transposedestinations`](@ref)).

The values are computed on the host, because a callable provider is an
arbitrary Julia function. A block whose data is tabulated or constant needs
no host at all; see [`plandeviceproviders`](@ref).
"""
struct DeviceScatteringStamps{VI,M,H,B}
    ssys::Any
    index::VI
    values::M               # ncontrib by nbatch, on the backend
    host::H
    wmodes::Vector{Float64}
    work::Any               # a ScatteringWorkspace, reused across batches
    backend::B
    nzA::Int
end

"""
    plandevicescattering(ssys, Aindexcsr, nzA, nbatch, backend, Nmodes)

Build a [`DeviceScatteringStamps`](@ref) from a scattering stamp system whose
destination indices have already been mapped into the stored order the batch
uses, or `nothing` when `ssys` is `nothing`.
"""
function plandevicescattering(ssys, Aindexcsr::AbstractVector{Int},
    nzA::Integer, nbatch::Integer, backend, Nmodes::Integer)

    isnothing(ssys) && return nothing
    n = length(Aindexcsr)
    Ti = nzA*nbatch <= typemax(Int32) ? Int32 : Int
    values = KernelAbstractions.allocate(backend, Complex{Float64}, n, nbatch)
    host = Matrix{Complex{Float64}}(undef, n, nbatch)
    return DeviceScatteringStamps(ssys,
        tobackend(backend, convert(Vector{Ti}, Aindexcsr)), values, host,
        zeros(Float64, Nmodes), ScatteringWorkspace(), backend, Int(nzA))
end

"""
    transposedestinations(st::DeviceScatteringStamps, Aindexcsr, backend)

A view of the same scattering stamp values with different destinations, for
the transposed (adjoint) system: the contribution of a block does not depend
on the direction the system is assembled in, only on where in the stored
order of that system each scalar lands. Sharing `values` is what lets the
providers be evaluated once for both directions.
"""
function transposedestinations(st::DeviceScatteringStamps,
    Aindexcsr::AbstractVector{Int}, backend)
    Ti = eltype(st.index)
    return DeviceScatteringStamps(st.ssys,
        tobackend(backend, convert(Vector{Ti}, Aindexcsr)), st.values,
        st.host, st.wmodes, st.work, backend, st.nzA)
end

# compute the values of the batch beginning at `lo` into host buffer `buf`.
# Pure host work, and the only thing which touches the providers.
function computescatteringvalues!(st::DeviceScatteringStamps, buf, w,
    lo::Integer, k::Integer, wpumpmodes)

    nbatch = size(buf, 2)
    wmodes = st.wmodes
    @inbounds for j in 1:k
        for m in eachindex(wmodes)
            wmodes[m] = w[lo + j - 1] + wpumpmodes[m]
        end
        scatteringvalues!(view(buf, :, j), st.ssys, wmodes, st.work)
    end
    # a short final batch: repeat the last real frequency, as the assembly
    # does, so nothing reads uninitialized memory
    @inbounds for j in (k+1):nbatch
        copyto!(view(buf, :, j), view(buf, :, k))
    end
    return buf
end

"""
    stagescatteringstamps!(st::DeviceScatteringStamps, w, lo, k, wpumpmodes)

Compute the scattering values of the batch beginning at `lo` on the host and
send them to the backend.
"""
function stagescatteringstamps!(st::DeviceScatteringStamps, w, lo::Integer,
    k::Integer, wpumpmodes)
    computescatteringvalues!(st, st.host, w, lo, k, wpumpmodes)
    copyto!(st.values, st.host)
    return st
end

"""
    applyscatteringstamps!(nzval::AbstractMatrix, st::DeviceScatteringStamps)

Add the staged scattering values into the stored values `nzval` of a batch of
system matrices. [`stagescatteringstamps!`](@ref) must have run for the batch
`nzval` was assembled at.
"""
function applyscatteringstamps!(nzval::AbstractMatrix,
    st::DeviceScatteringStamps)
    scatteringstampkernel!(st.backend, 64)(nzval, st.index, st.values,
        size(st.values, 1), st.nzA; ndrange = length(st.values))
    KernelAbstractions.synchronize(st.backend)
    return nzval
end

# ---------------------------------------------------------------------------
# evaluating the scattering blocks on the backend
# ---------------------------------------------------------------------------

# the interpolated scattering parameter S[p,q] of one block at frequency w,
# from its table. A block whose data is a single matrix is a one point table,
# so it needs no separate path.
@inline function tableentry(freqs, vals, foff, nf, voff, n, p, q, w,
        extrapcode)
    @inbounds begin
        s = n*n
        if nf == 1
            return vals[voff + (q-1)*n + p]
        end
        f1 = freqs[foff+1]
        fn = freqs[foff+nf]
        # which segment, by the same branchless search the assembly kernels
        # use over a column pointer
        lo = 1; hi = nf - 1
        if w <= f1
            lo = 1
        elseif w >= fn
            lo = nf - 1
        else
            while lo < hi
                mid = (lo + hi + 1) >>> 1
                if freqs[foff+mid] <= w
                    lo = mid
                else
                    hi = mid - 1
                end
            end
        end
        a = freqs[foff+lo]
        b = freqs[foff+lo+1]
        t = (w - a)/(b - a)
        # `:constant` holds the end value beyond the ends; `:linear` lets the
        # end segment continue, which is what an unclamped t already does
        if extrapcode == Int8(1)
            t = t < 0 ? zero(t) : (t > 1 ? one(t) : t)
        end
        k1 = voff + (lo-1)*s + (q-1)*n + p
        return vals[k1]*(1-t) + vals[k1+s]*t
    end
end

# the value one contribution adds, from the scattering entry it reads. The
# two kernels below differ only in where that entry comes from: a table or a
# callable.
@inline function hybridcontribution(::Type{T}, S, wm, samepq::Bool, zrefp,
        coeffc, sgnc, scale) where {T}
    Bpq = zero(T); Cpq = zero(T)
    if iszero(wm)
        # the zero frequency rows are i = 0
        Cpq = samepq ? one(T) : zero(T)
    else
        r2 = sqrt(zrefp)
        rinv2 = one(r2)/r2
        Bpq = -rinv2*S
        Cpq = r2*S
        if samepq
            Bpq += rinv2
            Cpq += r2
        end
    end
    return coeffc == 1 ? sgnc*(im*wm*scale)*Bpq : -Cpq
end

# one work item per (contribution, frequency of the batch): interpolate the
# block, form the one hybrid coefficient the contribution needs, and write the
# value the gather-add will add
@kernel function deviceproviderkernel!(values, @Const(modeindex),
        @Const(blockindex), @Const(pindex), @Const(qindex), @Const(coeff),
        @Const(sgn), @Const(nports), @Const(zrefoff), @Const(zref),
        @Const(freqoff), @Const(nfreq), @Const(freqs), @Const(valoff),
        @Const(vals), @Const(conjsym), @Const(extrapcode), @Const(wpump),
        @Const(ws), scale, ncontrib)
    gid = @index(Global)
    @inbounds begin
        g = gid - 1
        c = g % ncontrib + 1
        f = g ÷ ncontrib + 1
        m = Int(modeindex[c]); bi = Int(blockindex[c])
        p = Int(pindex[c]); q = Int(qindex[c])
        wm = ws[f] + wpump[m]
        n = Int(nports[bi])
        T = eltype(values)
        S = zero(T)
        if !iszero(wm)
            isconj = conjsym[bi] != 0
            wq = isconj ? abs(wm) : wm
            S = tableentry(freqs, vals, Int(freqoff[bi]), Int(nfreq[bi]),
                Int(valoff[bi]), n, p, q, wq, extrapcode[bi])
            if isconj && wm < 0
                S = conj(S)
            end
        end
        values[c, f] = hybridcontribution(T, S, wm, p == q,
            zref[Int(zrefoff[bi]) + p], coeff[c], sgn[c], scale)
    end
end

# the same, for blocks whose scattering parameters come from a callable of the
# `:entry` form. The callables live in a device array indexed by block: they
# capture only numbers, so they are `isbits`, and blocks built from one helper
# share a type, which is what lets them be stored and called this way.
@kernel function deviceentrykernel!(values, @Const(modeindex),
        @Const(blockindex), @Const(pindex), @Const(qindex), @Const(coeff),
        @Const(sgn), @Const(zrefoff), @Const(zref), @Const(funcs),
        @Const(conjsym), @Const(wpump), @Const(ws), scale, ncontrib)
    gid = @index(Global)
    @inbounds begin
        g = gid - 1
        c = g % ncontrib + 1
        f = g ÷ ncontrib + 1
        m = Int(modeindex[c]); bi = Int(blockindex[c])
        p = Int(pindex[c]); q = Int(qindex[c])
        wm = ws[f] + wpump[m]
        T = eltype(values)
        S = zero(T)
        if !iszero(wm)
            isconj = conjsym[bi] != 0
            wq = isconj ? abs(wm) : wm
            S = T(funcs[bi](p, q, wq))
            if isconj && wm < 0
                S = conj(S)
            end
        end
        values[c, f] = hybridcontribution(T, S, wm, p == q,
            zref[Int(zrefoff[bi]) + p], coeff[c], sgn[c], scale)
    end
end

"""
    DeviceProviders

The scattering blocks of a stamp system on a backend, as flat tables.

The values a batch of system matrices needs are otherwise computed on the
host, because a callable provider is an arbitrary Julia function. A block
whose data is tabulated or constant is not: its evaluation is a search and an
interpolation, which is the shape a kernel wants. When every block of a
circuit is one of those, the whole evaluation moves to the backend and the
host does nothing per frequency at all.

A constant block is stored as a one point table, so it needs no separate path.
Which frequency a block is evaluated at still depends on its negative
frequency rule, and a range which must not be extrapolated is checked on the
host before each batch, because a kernel cannot raise.
"""
struct DeviceProviders{VI,VZ,VR,VC,VF,B}
    nports::VI
    zrefoff::VI
    zref::VR
    freqoff::VI
    nfreq::VI
    freqs::VR
    valoff::VI
    vals::VC
    # the callables of the `:entry` form, or `nothing` when the blocks are
    # tabulated; exactly one of `vals` and `funcs` is used
    funcs::VF
    conjsym::VZ
    extrapcode::VZ
    modeindex::VI
    blockindex::VI
    pindex::VI
    qindex::VI
    coeff::VZ
    sgn::VZ
    wpump::VR
    ws::VR
    wshost::Vector{Float64}
    scale::Float64
    ncontrib::Int
    ranges::Vector{Tuple{Float64,Float64}}
    strict::Vector{Bool}
    conjhost::Vector{Bool}
    names::Vector{String}
    ssys::Any
    backend::B
end

"""
    candeviceevaluate(ssys)

Whether every block of a stamp system has data a kernel can evaluate, which
is tabulated or constant. A callable provider is an arbitrary Julia function
and stays on the host.
"""
function candeviceevaluate(ssys)
    isnothing(ssys) && return false
    isempty(ssys.blocks) && return false
    first = ssys.blocks[1].block.provider
    if first isa CallableMatrixProvider
        # a callable can be called from a kernel only in the `:entry` form,
        # and only if the closures can live in a device array, which needs
        # them to be `isbits` and all of one type
        first.form === :entry || return false
        isbits(first.f) || return false
        for sb in ssys.blocks
            p = sb.block.provider
            p isa CallableMatrixProvider || return false
            p.form === :entry || return false
            typeof(p.f) === typeof(first.f) || return false
        end
        return true
    end
    for sb in ssys.blocks
        p = sb.block.provider
        (p isa TabulatedMatrixProvider || p isa ConstantMatrixProvider) ||
            return false
        p isa TabulatedMatrixProvider && p.interpolation != :linear &&
            return false
    end
    return true
end

# the callables of an entry-wise stamp system, or `nothing` when its blocks
# are tabulated
function entrycallables(ssys)
    p1 = ssys.blocks[1].block.provider
    (p1 isa CallableMatrixProvider && p1.form === :entry) || return nothing
    return [sb.block.provider.f for sb in ssys.blocks]
end

# The stamp system holds its blocks untyped, so a loop which reads a
# provider's arrays element by element dispatches on every element. These
# barriers give the compiler the concrete types once.
function copyfrequencies!(dest::Vector{Float64}, off::Int,
    src::AbstractVector)
    @inbounds for k in eachindex(src)
        dest[off + k] = Float64(src[k])
    end
    return nothing
end

function copytable!(dest::Vector{Complex{Float64}}, off::Int,
    src::AbstractArray{<:Number,3}, n::Int, nf::Int)
    @inbounds for k in 1:nf, j in 1:n, i in 1:n
        dest[off + (k-1)*n*n + (j-1)*n + i] = src[i, j, k]
    end
    return nothing
end

function copymatrix!(dest::Vector{Complex{Float64}}, off::Int,
    src::AbstractMatrix{<:Number}, n::Int)
    @inbounds for j in 1:n, i in 1:n
        dest[off + (j-1)*n + i] = src[i, j]
    end
    return nothing
end

"""
    plandeviceproviders(ssys, nbatch, backend, wpumpmodes, scale)

Build a [`DeviceProviders`](@ref) for a stamp system whose blocks can all be
evaluated by a kernel, or `nothing` otherwise.
"""
function plandeviceproviders(ssys, nbatch::Integer, backend, wpumpmodes,
    scale::Real)

    candeviceevaluate(ssys) || return nothing
    nb = length(ssys.blocks)
    nports = Int32[sb.block.nports for sb in ssys.blocks]
    callables = entrycallables(ssys)
    istable = isnothing(callables)
    # size everything first and fill it in place: a line whose every cell is
    # its own block has hundreds of thousands of table points, and growing
    # the flat arrays a block at a time would be slow
    ntot = 0; vtot = 0; ztot = 0
    for sb in ssys.blocks
        p = sb.block.provider
        n = sb.block.nports
        ztot += n
        istable || continue
        if p isa ConstantMatrixProvider
            ntot += 1; vtot += n*n
        else
            nf = length(p.frequencies)
            ntot += nf; vtot += n*n*nf
        end
    end
    zrefoff = Vector{Int32}(undef, nb); zref = Vector{Float64}(undef, ztot)
    freqoff = Vector{Int32}(undef, nb); nfreq = Vector{Int32}(undef, nb)
    freqs = Vector{Float64}(undef, ntot)
    valoff = Vector{Int32}(undef, nb)
    vals = Vector{Complex{Float64}}(undef, vtot)
    conjsym = Vector{Int8}(undef, nb); extrapcode = Vector{Int8}(undef, nb)
    ranges = Vector{Tuple{Float64,Float64}}(undef, nb)
    strict = Vector{Bool}(undef, nb)
    conjhost = Vector{Bool}(undef, nb)
    names = String[sb.name for sb in ssys.blocks]
    zi = 0; fi = 0; vi = 0
    for (bi, sb) in enumerate(ssys.blocks)
        blk = sb.block
        n = blk.nports
        zrefoff[bi] = zi
        @inbounds for p in 1:n
            zref[zi + p] = real(blk.zref[p])
        end
        zi += n
        prov = blk.provider
        freqoff[bi] = fi
        valoff[bi] = vi
        if !istable
            # an entry-wise callable has no table, no range, and nothing to
            # extrapolate
            nfreq[bi] = 0
            extrapcode[bi] = Int8(1)
            ranges[bi] = (-Inf, Inf)
            strict[bi] = false
        elseif prov isa ConstantMatrixProvider
            # a single matrix is a one point table
            nfreq[bi] = 1
            freqs[fi + 1] = 0.0
            fi += 1
            copymatrix!(vals, vi, prov.A, n)
            vi += n*n
            extrapcode[bi] = Int8(1)
            ranges[bi] = (0.0, 0.0)
            strict[bi] = false
        else
            nf = length(prov.frequencies)
            nfreq[bi] = nf
            copyfrequencies!(freqs, fi, prov.frequencies)
            fi += nf
            copytable!(vals, vi, prov.values, n, nf)
            vi += n*n*nf
            extrapcode[bi] = prov.extrapolation == :linear ? Int8(2) : Int8(1)
            ranges[bi] = (Float64(prov.frequencies[1]),
                Float64(prov.frequencies[nf]))
            strict[bi] = prov.extrapolation == :error && nf > 1
        end
        conjhost[bi] = !(blk.negative_frequency isa Native)
        conjsym[bi] = conjhost[bi] ? Int8(1) : Int8(0)
    end

    d(x) = tobackend(backend, x)
    funcs = istable ? nothing : tobackend(backend, callables)
    ncontrib = length(ssys.Aindex)
    return DeviceProviders(d(nports), d(zrefoff), d(zref), d(freqoff),
        d(nfreq), d(freqs), d(valoff), d(vals), funcs, d(conjsym),
        d(extrapcode),
        d(Int32.(ssys.modeindex)), d(Int32.(ssys.blockindex)),
        d(Int32.(ssys.pindex)), d(Int32.(ssys.qindex)), d(Int8.(ssys.coeff)),
        d(Int8.(ssys.sign)), d(collect(Float64, wpumpmodes)),
        d(zeros(Float64, nbatch)), zeros(Float64, nbatch), Float64(scale),
        ncontrib, ranges, strict, conjhost, names, ssys, backend)
end

# a kernel cannot raise, so a block whose data must not be extrapolated has
# its range checked here, before the batch runs
function checkdeviceranges(dp::DeviceProviders, k::Integer)
    any(dp.strict) || return nothing
    lonat = Inf; hinat = -Inf; loabs = Inf; hiabs = -Inf
    @inbounds for j in 1:k, m in eachindex(dp.wshost)
        nothing
    end
    nm = length(dp.wpump)
    wp = Array(dp.wpump)
    @inbounds for j in 1:k
        for m in 1:nm
            w = dp.wshost[j] + wp[m]
            iszero(w) && continue
            lonat = min(lonat, w); hinat = max(hinat, w)
            a = abs(w)
            loabs = min(loabs, a); hiabs = max(hiabs, a)
        end
    end
    for bi in eachindex(dp.strict)
        dp.strict[bi] || continue
        lo, hi = dp.ranges[bi]
        l, h = dp.conjhost[bi] ? (loabs, hiabs) : (lonat, hinat)
        if l < lo || h > hi
            throw(ArgumentError(lazy"The scattering block at $(dp.names[bi]) is evaluated over [$(l), $(h)] rad/s but its tabulated range is [$(lo), $(hi)] rad/s. Extrapolation of tabulated data is opt-in: pass extrapolation = :constant or :linear if extrapolation is intended."))
        end
    end
    return nothing
end

"""
    stagedeviceproviders!(values, dp::DeviceProviders, w, lo, k)

Compute the scattering values of the batch of frequencies beginning at `lo`
directly on the backend, into `values`.

The host part of this is the batch's own signal frequencies and, for a block
whose data must not be extrapolated, a check of the range it will be
evaluated over. Nothing per block or per contribution crosses to the host.
"""
function stagedeviceproviders!(values::AbstractMatrix, dp::DeviceProviders,
    w, lo::Integer, k::Integer)

    nbatch = length(dp.wshost)
    @inbounds for j in 1:nbatch
        dp.wshost[j] = w[min(lo + j - 1, lo + k - 1)]
    end
    checkdeviceranges(dp, k)
    copyto!(dp.ws, dp.wshost)
    if isnothing(dp.funcs)
        deviceproviderkernel!(dp.backend, 64)(values, dp.modeindex,
            dp.blockindex, dp.pindex, dp.qindex, dp.coeff, dp.sgn, dp.nports,
            dp.zrefoff, dp.zref, dp.freqoff, dp.nfreq, dp.freqs, dp.valoff,
            dp.vals, dp.conjsym, dp.extrapcode, dp.wpump, dp.ws, dp.scale,
            dp.ncontrib; ndrange = length(values))
    else
        deviceentrykernel!(dp.backend, 64)(values, dp.modeindex,
            dp.blockindex, dp.pindex, dp.qindex, dp.coeff, dp.sgn,
            dp.zrefoff, dp.zref, dp.funcs, dp.conjsym, dp.wpump, dp.ws,
            dp.scale, dp.ncontrib; ndrange = length(values))
    end
    KernelAbstractions.synchronize(dp.backend)
    return values
end

