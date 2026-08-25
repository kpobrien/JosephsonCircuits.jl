
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
the system matrix instead, whose solutions are the adjoint ones the noise,
quantum efficiency and sensitivity calculations need. That transpose is free
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
    gatherportrows!(out, X, rows)

Gather the rows named by `rows` from every right hand side of every system of
the batch `X`, into `out`. See [`portsolutionrows`](@ref).
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
    needsadjointsolve(arrays::LinearizedArrays, noiseportimpedanceindices)

Whether the transposed (adjoint) linearized system must be solved at each
signal frequency: for the scattering parameter sensitivities always, and
otherwise when a consumer of the adjoint solution (the noise scattering
parameters, the quantum efficiency, the commutation relations, or the adjoint
node outputs) is requested together with a source of it.

This does not depend on the frequency, and both the host loop and the device
sweep test it, the latter because it has no factorization to solve the
transposed system against.
"""
function needsadjointsolve(arrays::LinearizedArrays, noiseportimpedanceindices)
    isempty(arrays.Ssensitivity) || return true
    hassource = !isempty(noiseportimpedanceindices) ||
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

Returns `(batchsize, solvebatch!, forward, adjoint, adjointdevice)`.
`solvebatch!(lo)` solves the batch of frequencies beginning at index `lo` and
stages its solutions on the host; `forward` and `adjoint` then fill the
solution of any frequency of that batch, in the form
[`hblinsolve_inner!`](@ref) takes as `presolved` and `presolvedadjoint`, with
`adjoint` `nothing` when none was asked for.

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
    function newbatch(adj)
        plan, rowptr, colind = planfrequencysweep(lsys, backend; adjoint = adj)
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

    # The staging for one direction. A whole batch is brought back at once and
    # read from host memory afterwards, so that reading a solution touches the
    # device not at all and any number of workers may do it at once.
    #
    # Gathering pays only while the named rows are a small part of the
    # solution. Past that the copy back is nearly the same size either way and
    # the gather is just work, so the whole solution is staged instead. A
    # circuit whose loss is spread along the line reaches this: its noise
    # ports touch almost every node, so the adjoint solutions the noise
    # scattering parameters read are almost the entire solution.
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

    function runbatch!(slot, b, stage)
        assemblesweep!(b.nzval, b.plan, wsdev)
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
        loadfrequencies!(lo)
        runbatch!(1, fwd, fstage)
        isnothing(adj) || runbatch!(2, adj, astage)
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

    # `adjointdevice` hands out the adjoint solution where it was computed,
    # for the noise scattering parameters, which are reduced there rather than
    # brought back (see [`devicenoise`](@ref)). It is a read of the solved
    # batch, so several workers may hold different frequencies of it at once.
    return (batchsize = nb,
        solvebatch! = solvebatch!,
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
    devicenoise(plan::DeviceNoisePlan, adjointsolution, nrhs, wpumpmodes, w,
        keepmatrix)

A callback which computes the noise scattering parameters of a signal
frequency from the adjoint solutions on the backend.

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
function devicenoise(plan::DeviceNoisePlan, adjointsolution, nrhs::Integer,
    wpumpmodes, w, keepmatrix::Bool)

    backend = plan.backend
    T = Complex{Float64}
    nrows = plan.nnoise*plan.nmodes
    out = KernelAbstractions.allocate(backend, T, nrows, nrhs)
    Snoise = KernelAbstractions.allocate(backend, T, nrows, nrhs)
    invinput = KernelAbstractions.allocate(backend, T, nrhs, nrhs)
    # the two reductions go through the backend's own `sum!`, which is a tree
    # over the noise index. Summing each column in one work item instead was
    # measured and it loses badly: it is one thread per port mode, a dozen of
    # them, each walking every noise port in turn.
    absq = KernelAbstractions.allocate(backend, Float64, nrows, nrhs)
    signs = KernelAbstractions.allocate(backend, Float64, nrows, 1)
    denomd = KernelAbstractions.allocate(backend, Float64, 1, nrhs)
    signedd = KernelAbstractions.allocate(backend, Float64, 1, nrhs)
    wmodesd = KernelAbstractions.allocate(backend, Float64, plan.nmodes)
    denom = zeros(Float64, nrhs)
    signed = zeros(Float64, nrhs)
    wmodes = zeros(Float64, plan.nmodes)
    reduction = NoiseReduction(denom, signed)

    return function(i, inputwave, Snoiseview)
        @inbounds for m in eachindex(wmodes)
            wmodes[m] = w[i] + wpumpmodes[m]
        end
        copyto!(wmodesd, wmodes)
        noiseoutputwavekernel!(backend, 64)(out, adjointsolution(i),
            plan.node1, plan.node2, plan.values, plan.codes, wmodesd,
            plan.nmodes, plan.nnoise, nrhs;
            ndrange = nrows*nrhs)
        KernelAbstractions.synchronize(backend)
        copyto!(invinput, inv(Matrix{T}(inputwave)))
        mul!(Snoise, out, invinput)
        absq .= abs2.(Snoise)
        sum!(denomd, absq)
        modesignkernel!(backend, 64)(signs, wmodesd, plan.nmodes;
            ndrange = nrows)
        KernelAbstractions.synchronize(backend)
        absq .*= signs
        sum!(signedd, absq)
        KernelAbstractions.synchronize(backend)
        copyto!(denom, vec(Array(denomd)))
        copyto!(signed, vec(Array(signedd)))
        keepmatrix && isempty(Snoiseview) == false && copyto!(Snoiseview,
            Array(Snoise))
        return reduction
    end
end
