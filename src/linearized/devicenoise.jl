# The noise scattering parameters, covariance, quantum efficiency and
# commutation relations formed on the backend from the adjoint solutions of
# a device sweep.

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
