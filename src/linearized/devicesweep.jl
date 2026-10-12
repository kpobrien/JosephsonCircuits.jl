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
with the stored value of every term conjugated where that mode frequency is
negative (see [`modevalue`](@ref)). Each stored entry is therefore an
independent quadratic in its own mode frequency, and the whole assembly
collapses to four constant coefficient vectors and one kernel:

    A[q] = cst[q] + sel(kinvL[q]) + im*wm*sel(kG[q]) - wm^2*sel(kC[q])

where `wm = ws + wpumpmodes[mode of q]` and `sel` conjugates when `wm < 0`.
The terms which share a frequency power and a conjugation rule are summed
into one coefficient at build time: `Amna0`, the coupled inductances, is
conjugated at a negative mode frequency as `invLnm` is, and goes with it,
while `AoLjnm`, the pump's junction term, is evaluated at the signed modes
and is not conjugated. Conjugation distributes over a sum, so this is
exact.

Because the coefficients do not depend on the signal frequency, one kernel
fills the stored values of a whole batch of frequencies, which is what a
uniform batch wants, whether cuDSS's or the batched block factorization's:
the batch shares one sparsity pattern and one symbolic analysis, and only
the values differ. The matrices hold the components whose values depend on
the frequency open, so the coefficients leave them out; their stamps are
added to the batch by [`applyfrequencydependent!`](@ref).

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
`(plan, rowptr, colind, colindhost)`. The column index array comes back on
the host as well, which is where the grouping of
[`planequilibration`](@ref) is built.

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
    A = lsys.Asparse
    n = size(A, 1)
    nz = nnz(A)

    # the four coefficients, in the stored order of A
    cst = copy(lsys.AoLjnmnzval)
    kinvL = zeros(Complex{Float64}, nz)
    scattercoefficient!(kinvL, lsys.invLnm, lsys.invLnmindexmap)
    scattercoefficient!(kinvL, lsys.Amna0, lsys.Amna0indexmap)
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
    # each array is moved to the backend once, and the plan's type read
    # from the moved arrays
    cstd = d(cst)
    wpump = tobackend(backend, lsys.wpumpmodes)
    plan = FrequencySweepPlan{typeof(colof),typeof(cstd),typeof(wpump),
        typeof(assemble!),typeof(backend)}(
        colof, cstd, d(kinvL), d(kG), d(kC), wpump, assemble!, backend,
        Int(lsys.Nmodes), nz)
    return plan, tobackend(backend, rowptrhost),
        tobackend(backend, colindhost), colindhost
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

# ---------------------------------------------------------------------------
# the frequency dependent components on a backend
# ---------------------------------------------------------------------------

# Add the frequency dependent stamps of one stored entry, for one frequency
# of the batch, into that frequency's system matrix: each stamp is the
# value of its function at its mode frequency, read from the batch's
# table, or the inverse of that value, with its sign, conjugated at a
# negative mode frequency and multiplied by its frequency factor, `-w^2`
# for a capacitance and `im*w` for a conductance, as
# `addfrequencydependent!` adds them on the host. The stamps are grouped by
# destination, as the scattering stamps are, so each address is written by
# one work item and the sum needs no atomic.
@kernel function frequencydependentkernel!(nzval, @Const(dest),
        @Const(order), @Const(segptr), @Const(mode), @Const(tablerow),
        @Const(kind), @Const(sign), @Const(table), @Const(wpump), @Const(ws),
        ndest, nzA)
    gid = @index(Global)
    @inbounds begin
        q = gid - 1
        k = q % ndest + 1
        f = q ÷ ndest + 1
        acc = zero(eltype(table))
        for s in Int(segptr[k]):(Int(segptr[k+1]) - 1)
            p = Int(order[s])
            wm = ws[f] + wpump[Int(mode[p])]
            c = kind[p]
            v = modevalue(termvalue(c, sign[p], table[Int(tablerow[p]), f]),
                wm)
            acc += c == TERM_CAPACITANCE ? -(wm*wm)*v :
                c == TERM_CONDUCTANCE ? (im*wm)*v : v
        end
        nzval[Int(dest[k]) + (f-1)*nzA] += acc
    end
end

"""
    DeviceFrequencyDependentStamps

The stamps of the frequency dependent components of a linearized system
(see [`FrequencyDependentStamps`](@ref)) on a backend, for a batch of
signal frequencies.

[`assemblesweep!`](@ref) covers the terms which are a constant quadratic in
the signal frequency, and the matrices it reads hold these components
open. Their values are those of their functions at the mode frequencies,
arbitrary Julia callables, so the functions are evaluated on the host, each
distinct one once a mode frequency, into the batch's `table`, a column per
frequency and a row per function and mode, which is sent; a kernel then
forms each stamp from its row of the table and adds it
([`applyfrequencydependent!`](@ref)).

The table serves both directions of a sweep: the value of a stamp does not
depend on the direction the system is assembled in, only its destination
does, so the adjoint direction gets a second view with destinations of its
own (see [`transposedestinations`](@ref)).
"""
struct DeviceFrequencyDependentStamps{VI,VK,M,B}
    providers::Vector{Any}
    # the stamps grouped by destination; see `groupdestinations`
    dest::VI
    order::VI
    segptr::VI
    # of each stamp: its mode, the row of its function's value at that mode
    # in a column of the table, its kind and its sign
    mode::VI
    tablerow::VI
    kind::VK
    sign::VK
    table::M                # nproviders*Nmodes by nbatch, on the backend
    host::Matrix{Complex{Float64}}
    # one frequency's values, a row per function and a column per mode, and
    # its mode frequencies
    values::Matrix{Complex{Float64}}
    wmodes::Vector{Float64}
    backend::B
    nzA::Int
end

"""
    plandevicefrequencydependent(fd, dest, nzA, nbatch, backend, Nmodes)

Build a [`DeviceFrequencyDependentStamps`](@ref) from the stamps `fd` of a
linearized system, with their destinations `dest` in the stored order the
batch uses (see [`sweepdestinations`](@ref)), or `nothing` when `fd` is
`nothing`.
"""
function plandevicefrequencydependent(fd, dest::AbstractVector{Int},
    nzA::Integer, nbatch::Integer, backend, Nmodes::Integer)
    isnothing(fd) && return nothing
    np = length(fd.providers)
    Ti = nzA*nbatch < typemax(Int32) ? Int32 : Int
    order, groups, segptr = groupdestinations(dest)
    # copies, which the backend adopts on the host
    d = x -> tobackend(backend, Vector{Ti}(x))
    return DeviceFrequencyDependentStamps(fd.providers, d(groups), d(order),
        d(segptr), d(fd.mode), d(fd.provider .+ np .* (fd.mode .- 1)),
        tobackend(backend, copy(fd.kind)), tobackend(backend, copy(fd.sign)),
        KernelAbstractions.allocate(backend, Complex{Float64}, np*Nmodes,
            nbatch),
        Matrix{Complex{Float64}}(undef, np*Nmodes, nbatch),
        Matrix{Complex{Float64}}(undef, np, Nmodes), zeros(Float64, Nmodes),
        backend, Int(nzA))
end

# the same stamps and table with the destinations of the transposed
# (adjoint) system
function transposedestinations(st::DeviceFrequencyDependentStamps,
    dest::AbstractVector{Int}, backend)
    Ti = eltype(st.dest)
    order, groups, segptr = groupdestinations(dest)
    d = x -> tobackend(backend, Vector{Ti}(x))
    return DeviceFrequencyDependentStamps(st.providers, d(groups), d(order),
        d(segptr), st.mode, st.tablerow, st.kind, st.sign, st.table, st.host,
        st.values, st.wmodes, backend, st.nzA)
end

"""
    stagefrequencydependent!(st::DeviceFrequencyDependentStamps, w, lo, k,
        wpumpmodes)

Evaluate the functions of the frequency dependent components on the host
at the mode frequencies of the batch of `k` signal frequencies of `w`
beginning at `lo`, each distinct function once a mode frequency, and send
the table to the backend. A short final batch repeats its last frequency,
as its assembly does.
"""
function stagefrequencydependent!(st::DeviceFrequencyDependentStamps, w,
    lo::Integer, k::Integer, wpumpmodes)
    n = length(st.values)
    @inbounds for j in 1:size(st.host, 2)
        if j <= k
            for m in eachindex(st.wmodes)
                st.wmodes[m] = w[lo + j - 1] + wpumpmodes[m]
            end
            providervalues!(st.values, st.providers, st.wmodes)
            copyto!(st.host, (j - 1)*n + 1, st.values, 1, n)
        else
            copyto!(st.host, (j - 1)*n + 1, st.host, (k - 1)*n + 1, n)
        end
    end
    copyto!(st.table, st.host)
    return st
end

"""
    applyfrequencydependent!(nzval::AbstractMatrix,
        st::DeviceFrequencyDependentStamps, wpump, ws)

Add the stamps of the frequency dependent components into the stored
values `nzval` of a batch of system matrices at the signal frequencies
`ws`, from the table [`stagefrequencydependent!`](@ref) staged for the
batch; `wpump` holds the pump mode offsets of the signal modes.
"""
function applyfrequencydependent!(nzval::AbstractMatrix,
    st::DeviceFrequencyDependentStamps, wpump, ws)
    ndest = length(st.dest)
    frequencydependentkernel!(st.backend, 64)(nzval, st.dest, st.order,
        st.segptr, st.mode, st.tablerow, st.kind, st.sign, st.table, wpump,
        ws, ndest, st.nzA; ndrange = ndest*size(st.table, 2))
    KernelAbstractions.synchronize(st.backend)
    return nzval
end

"""
    SweepEquilibration

The row and column scaling of the system matrices of one direction of a
batch.

cuDSS chooses its pivots within the blocks its ordering fixes, by the size
of the candidates, and perturbs a pivot below an absolute threshold (see
the extension's `configure!`), so the units the equations and the unknowns
are written in decide which candidates it can tell apart and which it
perturbs. The linearized system is assembled in physical units: a node's
row holds admittances, while the constitutive row of a coupled inductor
branch holds the incidence of its node fluxes and, on its diagonal, the
inductance itself, many orders smaller. So each system of a batch is
scaled here: every row is divided by its largest entry, then every column
of the result by its largest, which makes the pivot candidates of every
row and column comparable; the right-hand side of each system is divided
by the same row scales, and its solution by the column scales on the way
out. cuDSS's matching, which permutes and scales the matrix ahead of the
factorization, is not supported for a uniform batch.

# Fields
- `rowptr`: the row pointer of the structure the solver is handed, whose
    rows are contiguous in the stored values.
- `order`, `segptr`: the stored entries of that structure grouped by its
    column, so one work item owns a column and no two write the same scale.
- `rowscale`, `scale`: the scale of each row and of each column of each
    system of the batch.
"""
struct SweepEquilibration{VR,VI,VS}
    rowptr::VR
    order::VI
    segptr::VI
    rowscale::VS
    scale::VS
end

# one work item per (row, system): the largest entry of the row, then the
# row divided by it
@kernel function rowequilibratekernel!(nzval, rowscale, @Const(rowptr), n)
    gid = @index(Global)
    @inbounds begin
        i = (gid - 1) % n + 1
        k = (gid - 1) ÷ n + 1
        lo = Int(rowptr[i]); hi = Int(rowptr[i+1]) - 1
        m = zero(real(eltype(nzval)))
        for q in lo:hi
            m = max(m, abs(nzval[q, k]))
        end
        s = m > 0 ? m : one(m)
        rowscale[i, k] = s
        for q in lo:hi
            nzval[q, k] /= s
        end
    end
end

# one work item per right-hand side entry of each system: the unscaled
# right-hand side divided by the system's row scale
@kernel function scalerighthandsidekernel!(B, @Const(B0), @Const(rowscale),
        n, nrhs)
    gid = @index(Global)
    @inbounds begin
        i = (gid - 1) % n + 1
        r = ((gid - 1) ÷ n) % nrhs + 1
        k = (gid - 1) ÷ (n*nrhs) + 1
        B[i, r, k] = B0[i, r]/rowscale[i, k]
    end
end

# one work item per (column, system): the largest entry of the column, then
# the column divided by it
@kernel function columnequilibratekernel!(nzval, scale, @Const(order),
        @Const(segptr), n)
    gid = @index(Global)
    @inbounds begin
        j = (gid - 1) % n + 1
        k = (gid - 1) ÷ n + 1
        lo = Int(segptr[j]); hi = Int(segptr[j+1]) - 1
        m = zero(real(eltype(nzval)))
        for q in lo:hi
            m = max(m, abs(nzval[Int(order[q]), k]))
        end
        s = m > 0 ? m : one(m)
        scale[j, k] = s
        for q in lo:hi
            nzval[Int(order[q]), k] /= s
        end
    end
end

# one work item per solution entry, undoing the change of variable
@kernel function unscalesolutionkernel!(X, @Const(scale), n, nrhs)
    gid = @index(Global)
    @inbounds begin
        i = (gid - 1) % n + 1
        r = ((gid - 1) ÷ n) % nrhs + 1
        k = (gid - 1) ÷ (n*nrhs) + 1
        X[i, r, k] /= scale[i, k]
    end
end

# the stored entries of each column of the structure handed to the solver,
# grouped by column: a counting sort of the column index array
function groupstoredcolumns(colind::AbstractVector{<:Integer}, n::Integer)
    segptr = zeros(Int32, n + 1)
    @inbounds for c in colind
        segptr[c + 1] += 1
    end
    segptr[1] = 1
    @inbounds for j in 1:n
        segptr[j+1] += segptr[j]
    end
    next = copy(segptr)
    order = Vector{Int32}(undef, length(colind))
    @inbounds for k in eachindex(colind)
        c = colind[k]
        order[next[c]] = k
        next[c] += 1
    end
    return order, segptr
end

"""
    planequilibration(rowptr, colind::AbstractVector{<:Integer}, n, nb,
        backend)

The [`SweepEquilibration`](@ref) of a batch of `nb` systems of order `n`
handed to the solver as the structure `rowptr` (on the backend) and
`colind` (on the host).
"""
function planequilibration(rowptr, colind::AbstractVector{<:Integer},
    n::Integer, nb::Integer, backend)
    order, segptr = groupstoredcolumns(colind, n)
    return SweepEquilibration(rowptr, tobackend(backend, order),
        tobackend(backend, segptr),
        KernelAbstractions.allocate(backend, Float64, Int(n), Int(nb)),
        KernelAbstractions.allocate(backend, Float64, Int(n), Int(nb)))
end

"""
    equilibrate!(nzval::AbstractMatrix, B, B0, eq::SweepEquilibration,
        backend)

Divide each row of each system of the batch by its largest entry, then
each column of the result by its largest, recording the scales in `eq`,
and write into `B` the right-hand side `B0` of every system divided by
that system's row scales.
"""
function equilibrate!(nzval::AbstractMatrix, B, B0, eq::SweepEquilibration,
    backend)
    n, nb = size(eq.scale)
    rowequilibratekernel!(backend, 64)(nzval, eq.rowscale, eq.rowptr, n;
        ndrange = n*nb)
    columnequilibratekernel!(backend, 64)(nzval, eq.scale, eq.order,
        eq.segptr, n; ndrange = n*nb)
    nrhs = size(B0, 2)
    scalerighthandsidekernel!(backend, 64)(B, B0, eq.rowscale, n, nrhs;
        ndrange = n*nrhs*nb)
    KernelAbstractions.synchronize(backend)
    return nzval
end

"""
    unscalesolution!(X::AbstractArray{<:Any,3}, eq::SweepEquilibration,
        backend)

Undo on the solutions of a batch the column scaling [`equilibrate!`](@ref)
applied to its matrices; the row scaling changes the equations and not
the unknowns, so it has nothing to undo.
"""
function unscalesolution!(X::AbstractArray{<:Any,3}, eq::SweepEquilibration,
    backend)
    n, nrhs, nb = size(X)
    unscalesolutionkernel!(backend, 64)(X, eq.scale, n, nrhs;
        ndrange = n*nrhs*nb)
    KernelAbstractions.synchronize(backend)
    return X
end

"""
    portsolutionrows(nodeindices, portindices, Nmodes::Integer)

The rows of a solution of the linearized system which the scattering
parameter calculation reads: for each port, each of its two nodes which is
not ground, and each mode.

[`calcoutputwaves!`](@ref) reads a solution only through
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
        noiseportimpedanceindices, noiseplan)

Whether the transposed (adjoint) linearized system must be solved at each
signal frequency: for the scattering parameter sensitivities always, and
otherwise when a consumer of the adjoint solution (the noise scattering
parameters, the noise covariances, the quantum efficiency, the occupations,
the commutation relations, or the adjoint node outputs) is requested
together with a source of it. The dissipative scattering blocks of the
[`ScatteringNoisePlan`](@ref) `noiseplan`, `nothing` without any, are such
a source, as the lumped noise ports are.

This does not depend on the frequency, and both the host loop and the device
sweep test it, the latter to decide whether to allocate and solve the
adjoint direction at all (with cuDSS a whole second factorization, with a
block factorization a second solve against the same factors).
"""
function needsadjointsolve(arrays::LinearizedArrays,
    noiseportimpedanceindices, noiseplan)
    isempty(arrays.Ssensitivity) || return true
    hassource = !isempty(noiseportimpedanceindices) ||
        !isnothing(noiseplan) ||
        !isempty(arrays.nodefluxadjoint) || !isempty(arrays.voltageadjoint)
    hasconsumer = readsnoise(arrays) ||
        !isempty(arrays.nodefluxadjoint) || !isempty(arrays.voltageadjoint)
    return hassource && hasconsumer
end

# whether an output of `arrays` reads the noise channels: the noise
# scattering parameters, the noise covariances, the quantum efficiency,
# the occupations or the commutation relations
readsnoise(arrays::LinearizedArrays) = !isempty(arrays.Snoise) ||
    !isempty(arrays.QE) || !isempty(arrays.CM) || !isempty(arrays.Cnoise) ||
    !isempty(arrays.nbar) || !isempty(arrays.Vout)

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
    nzval, X, B, B0, equil)` on the sparse device factorization path, `B0`
    the unscaled right-hand side each batch scales into `B`, or `(X,)` on
    the block path, and `nothing` when no adjoint was asked for.
- `fstage`, `astage`: the staging of each direction, `(full, rows, rowsd,
    gathered, host)`.
- `blocks`: the block path's `(plan, nzval, B, F, X, Xadj)`, or `nothing`.
- `scatstamps`, `scatstampsadjoint`, `providers`: the scattering block
    stamps of each direction and the providers they are evaluated through.
- `fdstamps`, `fdstampsadjoint`: the stamps of the frequency dependent
    components of each direction ([`DeviceFrequencyDependentStamps`](@ref)),
    or `nothing` when there are none.
- `sweeps`: the cuDSS sweep of each direction, made on the first batch.
- `solverkwargs`: the factorization's options, forwarded to each sweep as
    it is made so a caller can override the sweep's own defaults.
- `batchlo`: the first signal frequency of the batch currently staged.
"""
struct DeviceSweep{TB,TFw,TAd,TFs,TAs,TBl,TSt,TSa,TPr,TFd,TFa,TK<:NamedTuple}
    backend::TB
    nb::Int
    F::Int
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
    fdstamps::TFd
    fdstampsadjoint::TFa
    sweeps::Vector{Any}
    solverkwargs::TK
    batchlo::Base.RefValue{Int}
end

"""
    devicesolutions(lsys::HBLinearizedSystem, bnm, w, backend, forward,
        adjoint = nothing; factorization = nothing, refine = true,
        budget = memorybudget(backend))

The [`DeviceSweep`](@ref) which computes the solutions of the linearized
system on `backend` a batch of frequencies at a time.
`factorization` is the linearized solve's factorization: a
[`BlockFactorization`](@ref) takes the batched block path below, anything
else the cuDSS one. `refine` asks single precision block factors to refine
against the double residual; `false` is the fully single precision sweep.
`budget` is the device memory the batch may take.

The sweep is driven by [`solvebatch!`](@ref), which solves the batch of
frequencies beginning at an index and stages its solutions on the host;
[`forwardsolution!`](@ref) and [`adjointsolution!`](@ref) then fill the
solution of any frequency of that batch, and [`adjointdevice`](@ref) hands
out the adjoint solution of a frequency where it was computed, for the noise
scattering parameters, which are formed there rather than brought back (see
[`devicenoise`](@ref)). The `providers` field is the
[`DeviceProviders`](@ref) the scattering blocks are evaluated through, or
`nothing` when there are none or they cannot be evaluated on `backend`; the
noise channels of the dissipative blocks are read through the same one. The
components whose values depend on the frequency are stamped from a table
of their functions' values at the batch's modes, which the host fills for
each batch ([`DeviceFrequencyDependentStamps`](@ref)).

The split is what keeps the host work parallel. `solvebatch!` is the only verb
which touches the device; once it has returned, its batch's frequencies can be
post-processed by as many workers as the host path uses, because reading a
solution touches staged host memory and nothing else.

The systems of a batch share one sparsity pattern, so cuDSS, or the block
factorization ([`blockanalysis`](@ref)), analyzes it once and then factorizes
and solves each batch together, from the values of one
[`assemblesweep!`](@ref). With a sparse device factorization the batch
size is capped by [`uniformbatchlimit`](@ref) and by `budget`
([`memorybudget`](@ref) by default) at cuDSS's own estimate of the batch's
factorization, from the analysis of its forward direction, with the values,
solutions and right hand sides of each direction
([`cudssbatchlimit`](@ref)): the largest batch those arrays allow is
analyzed and kept when the estimate fits beside them, and a smaller one
made when it does not. With a [`BlockFactorization`](@ref) the cap does
not apply and the batch is sized by [`blocksystembytes`](@ref), the
factors, the originals when refining, the solutions of both directions and
the value matrix of one system, against `budget` and the length of the
sweep.

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
    adjoint = nothing; factorization = nothing, refine::Bool = true,
    budget::Integer = memorybudget(backend))

    T = Complex{Float64}
    n = size(lsys.Asparse, 1)
    nrhs = size(bnm, 2)
    nzA = nnz(lsys.Asparse)
    F = length(w)
    # the right hand sides do not depend on the frequency, and are the same for
    # both directions, but cuDSS wants one set per system of the batch
    bhost = convert(Matrix{T}, bnm)
    # a block factorization solves both directions from one factorization
    # per frequency, filled from the values assembled in the stored order;
    # its batch is sized by the memory of the factors and the solutions of
    # one system, within the device's memory budget. The cuDSS batch cap
    # does not apply to it; the memory and the sweep do. A cuDSS batch is
    # held to the same budget, at cuDSS's own estimate of its factorization
    # and the values, solutions and right hand sides each system holds, in
    # each direction. The estimate needs an analysis, which the batch size
    # must precede, so the largest batch the cap, the sweep and those arrays
    # allow is made and its forward direction analyzed; the first batch
    # factorizes that analysis. When the estimate does not fit beside the
    # arrays, a smaller batch is made instead, analyzed at its first batch.
    usesblocks = factorization isa BlockFactorization
    fwdbatch, fwdsweep = nothing, nothing
    nb = if usesblocks
        blocksym = blocksymbolic(lsys.Asparse, lsys.Nmodes)
        Tf = something(factorization.precision, Float64)
        # the solutions of each direction, the shared right-hand side, the
        # values, and when refining the residual and correction of a solve
        refining = refine && Tf === Float32
        persystem = blocksystembytes(Complex{Tf}, blocksym; refine = refining,
            TA = T) + ((isnothing(adjoint) ? 2 : 3) + (refining ? 2 : 0))*
            n*nrhs*sizeof(T) + nzA*sizeof(T)
        clamp(budget ÷ persystem, 1, F)
    else
        directions = isnothing(adjoint) ? 1 : 2
        arrays = nzA*sizeof(T) + 2n*nrhs*sizeof(T)
        largest = min(cudssbatchlimit(directions*arrays, nrhs, budget), F)
        fwdbatch = devicebatch(lsys, backend, false, bhost, largest)
        fwdsweep = _cudss_sweep(fwdbatch.rowptr, fwdbatch.colind,
            fwdbatch.nzval, fwdbatch.X, fwdbatch.B;
            solverkwargs(factorization)...)
        persystem = directions*(cld(_cudss_sweepbytes(fwdsweep), largest) +
            arrays)
        fits = min(cudssbatchlimit(persystem, nrhs, budget), largest)
        if fits < largest
            _cudss_release!(fwdsweep)
            releasearrays!(fwdbatch)
            fwdbatch = devicebatch(lsys, backend, false, bhost, fits)
            fwdsweep = nothing
        end
        fits
    end

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
    # the frequency dependent components: one table of their functions'
    # values serves both directions, whose stamps differ only in where
    # they land
    fd = lsys.frequencydependent
    fdstamps = plandevicefrequencydependent(fd, isnothing(fd) ? Int[] :
        sweepdestinations(lsys.Asparse, fd.Aindex, usesblocks), nzA, nb,
        backend, lsys.Nmodes)
    fdstampsadjoint = if !isnothing(fd) && !isnothing(adjoint) && !usesblocks
        transposedestinations(fdstamps, sweepdestinations(lsys.Asparse,
            fd.Aindex, true), backend)
    else
        fdstamps
    end

    fstage = devicestage(forward, n, nrhs, nb, backend)
    astage = isnothing(adjoint) ? nothing :
        devicestage(adjoint, n, nrhs, nb, backend)
    # the block path: the values of a batch in the stored (column) order,
    # one shared right-hand side, the solutions of each direction, and the
    # factorization's analysis of the pattern on the backend, which each
    # batch factorizes from its own values
    blocks = if usesblocks
        plan, _, _ = planfrequencysweep(lsys, backend; adjoint = true)
        Fb = blockanalysis(factorization, lsys.Asparse;
            blocksize = lsys.Nmodes, backend = backend, nb = nb,
            refine = refine ? BLOCKREFINESTEPS : 0)
        (plan = plan,
            nzval = KernelAbstractions.allocate(backend, T, nzA, nb),
            B = tobackend(backend, bhost), F = Fb,
            X = KernelAbstractions.allocate(backend, T, n, nrhs, nb),
            Xadj = isnothing(adjoint) ? nothing :
                KernelAbstractions.allocate(backend, T, n, nrhs, nb))
    else
        nothing
    end
    fwd = usesblocks ? (X = blocks.X,) : fwdbatch
    adj = isnothing(adjoint) ? nothing :
        usesblocks ? (X = blocks.Xadj,) :
        devicebatch(lsys, backend, true, bhost, nb)

    wshost = zeros(Float64, nb)
    wsdev = tobackend(backend, wshost)

    return DeviceSweep(backend, nb, F, collect(Float64, w),
        collect(Float64, lsys.wpumpmodes), wshost, wsdev, fwd, adj, fstage,
        astage, blocks, scatstamps, scatstampsadjoint, scatproviders,
        fdstamps, fdstampsadjoint, Any[fwdsweep, nothing],
        solverkwargs(factorization), Ref(0))
end

"""
    releasesweep!(ds::DeviceSweep)

Return the device memory of a finished sweep at once rather than when the
collector finds it ([`releasearray!`](@ref)): the factors, values,
solutions and right-hand sides of its batches, the staging of its
solutions, and cuDSS's factorizations. The next solve or sweep sizes
itself against the free memory ([`freememory`](@ref)), which memory the
collector has not yet found would shrink. `ds` is not used afterwards.
"""
function releasesweep!(ds::DeviceSweep)
    for part in (ds.fwd, ds.adj, ds.blocks, ds.fstage, ds.astage)
        releasearrays!(part)
    end
    for S in ds.sweeps
        isnothing(S) || _cudss_release!(S)
    end
    return nothing
end

# the arrays of a sweep's batches, through its named tuples, its assembly
# plans and equilibrations, the block factorization and its blocks
releasearrays!(x) = nothing
releasearrays!(x::AbstractArray{<:Number}) = releasearray!(x)
releasearrays!(x::Union{Tuple,NamedTuple,AbstractVector{<:AbstractArray}}) =
    foreach(releasearrays!, x)
releasearrays!(x::AbstractDict) = foreach(releasearrays!, values(x))
releasearrays!(lu::BlockLU) =
    foreach(releasearrays!, (lu.D, lu.L, lu.U, lu.Dinv, lu.scratch))
releasearrays!(F::SparseBlockFactorization) =
    foreach(releasearrays!, (F.lu, F.original, F.work, F.scale))
releasearrays!(p::FrequencySweepPlan) =
    foreach(releasearrays!, (p.colof, p.cst, p.kinvL, p.kG, p.kC, p.wpump))
releasearrays!(e::SweepEquilibration) =
    foreach(releasearrays!, (e.rowptr, e.order, e.segptr, e.rowscale, e.scale))

# the uniform batch of one direction on the sparse device factorization
# path: its sweep plan, pattern, values, solutions and right hand sides
function devicebatch(lsys, backend, isadjoint::Bool, bhost::Matrix, nb::Int)
    T = eltype(bhost)
    n = size(lsys.Asparse, 1)
    nrhs = size(bhost, 2)
    plan, rowptr, colind, colindhost = planfrequencysweep(lsys, backend;
        adjoint = isadjoint)
    nzval = KernelAbstractions.allocate(backend, T, nnz(lsys.Asparse), nb)
    X = KernelAbstractions.allocate(backend, T, n, nrhs, nb)
    B = KernelAbstractions.allocate(backend, T, n, nrhs, nb)
    fill!(X, zero(T))
    # the unscaled right-hand side, from which each batch writes every
    # system's own scaled one into `B`
    B0 = tobackend(backend, bhost)
    return (plan = plan, rowptr = rowptr, colind = colind, nzval = nzval,
        X = X, B = B, B0 = B0,
        equil = planequilibration(rowptr, colindhost, n, nb, backend))
end

# The staging for one direction. A whole batch is brought back at once and
# read from host memory afterwards, so reading a solution does not touch
# the device and any number of workers may do it at once.
#
# Gathering only the named rows pays while they are a small part of the
# solution, below `gatherfraction` of its rows; past that the whole
# solution is staged instead. A circuit whose loss is spread along the line
# reaches this, since its noise ports touch almost every node.
function devicestage(spec, n::Int, nrhs::Int, nb::Int, backend;
    gatherfraction::Real = 1/4)
    T = Complex{Float64}
    full = spec.full || length(spec.rows) >= gatherfraction*n
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
function runblockbatch!(ds::DeviceSweep, stamps, fdstamps)
    bl = ds.blocks
    assemblesweep!(bl.nzval, bl.plan, ds.wsdev)
    isnothing(stamps) || applyscatteringstamps!(bl.nzval, stamps)
    isnothing(fdstamps) || applyfrequencydependent!(bl.nzval, fdstamps,
        bl.plan.wpump, ds.wsdev)
    fillandfactorize!(bl.F, bl.nzval)
    refinedsolve!(bl.X, bl.F, bl.B)
    isnothing(bl.Xadj) || refinedsolve!(bl.Xadj, bl.F, bl.B; transposed = true)
    stagebatch!(ds, ds.fwd, ds.fstage)
    isnothing(ds.adj) || stagebatch!(ds, ds.adj, ds.astage)
    return nothing
end

# the sparse device factorization path of one direction: the batch is
# assembled, equilibrated, analyzed on the first batch unless the sizing
# analyzed it, factorized and solved, and the solutions are returned to the
# unscaled unknowns
function runbatch!(ds::DeviceSweep, slot::Int, b, stage, stamps, fdstamps)
    assemblesweep!(b.nzval, b.plan, ds.wsdev)
    isnothing(stamps) || applyscatteringstamps!(b.nzval, stamps)
    isnothing(fdstamps) || applyfrequencydependent!(b.nzval, fdstamps,
        b.plan.wpump, ds.wsdev)
    equilibrate!(b.nzval, b.B, b.B0, b.equil, ds.backend)
    if isnothing(ds.sweeps[slot])
        ds.sweeps[slot] = _cudss_sweep(b.rowptr, b.colind, b.nzval, b.X, b.B;
            ds.solverkwargs...)
    end
    _cudss_sweepsolve!(ds.sweeps[slot])
    unscalesolution!(b.X, b.equil, ds.backend)
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
    isnothing(ds.fdstamps) || stagefrequencydependent!(ds.fdstamps, ds.w, lo,
        k, ds.wpumpmodes)
    if !isnothing(ds.blocks)
        runblockbatch!(ds, ds.scatstamps, ds.fdstamps)
    else
        runbatch!(ds, 1, ds.fwd, ds.fstage, ds.scatstamps, ds.fdstamps)
        isnothing(ds.adj) || runbatch!(ds, 2, ds.adj, ds.astage,
            ds.scatstampsadjoint, ds.fdstampsadjoint)
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
    stagedtable!(table, ds::DeviceSweep, i::Integer)

Fill `table` with the values of the frequency dependent functions at the
modes of frequency `i` of the batch staged by the last
[`solvebatch!`](@ref), a row per function and a column per mode, as
[`stagefrequencydependent!`](@ref) evaluated them for its assembly. The
form [`hblinsolve_inner!`](@ref) takes as `presolvedtable`.
"""
function stagedtable!(table::Matrix{ComplexF64}, ds::DeviceSweep, i::Integer)
    n = length(table)
    copyto!(table, 1, ds.fdstamps.host, (batchslot(ds, i) - 1)*n + 1, n)
    return table
end

"""
    adjointdevice(ds::DeviceSweep, i::Integer)

The adjoint solution of frequency `i` of the staged batch where it was
computed, a view into the device batch, for the noise scattering parameters
which are formed there rather than brought back (see
[`devicenoise`](@ref)). A read of the solved batch, so several workers may
hold different frequencies of it at once.
"""
function adjointdevice(ds::DeviceSweep, i::Integer)
    # the view built as such: `view` makes a contiguous view of a GPU array
    # an array of its own, which shares the batch's memory until the
    # collector finds it and so keeps a finished sweep's solutions from the
    # pool (see `releasesweep!`)
    X = ds.adj.X
    return SubArray(X, (Base.Slice(axes(X, 1)), Base.Slice(axes(X, 2)),
        batchslot(ds, i)))
end
