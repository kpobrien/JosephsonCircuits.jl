
# === lowering value for scattering blocks ===

"""
    ScatteringStamp(block, port)

The value carried by one port of a [`ScatteringParameters`](@ref) lowered into a
[`ParsedSortedCircuit`](@ref): the shared block definition and the port
index. Each port of an n-port block lowers to one `:S` component whose two
nodes are the signal and reference terminals of that port, so the two node
per component layout of the parsed circuit is preserved; all ports of one
block share the block object by identity, which is how the solver
reassembles the full multiport coupling.
"""
struct ScatteringStamp{B}
    block::B
    port::Int
end

# === hybrid wave to modified nodal analysis coefficients ===

"""
    HybridWorkspace()

Reusable scratch for [`evaluatehybrid!`](@ref).

Every call needs the same five small buffers: which frequencies are nonzero,
those frequencies, the scattering parameters there, and the two square roots
of the reference impedances. A line whose every cell is its own block
evaluates thousands of blocks at every signal frequency, so allocating those
per call dominated the evaluation. They are reallocated only when a block
needs a larger one, which for a circuit whose blocks are all the same size is
once.

A workspace is mutable and is written by every call, so it belongs to one
task. [`scatteringvalues!`](@ref) makes one per call, which is what keeps it
safe to call from several threads at once.
"""
mutable struct HybridWorkspace
    nonzeroindices::Vector{Int}
    wnz::Vector{Float64}
    S::Array{Complex{Float64},3}
    rinv2::Vector{Float64}
    r2::Vector{Float64}
    absws::Vector{Float64}
end

HybridWorkspace() = HybridWorkspace(Int[], Float64[],
    Array{Complex{Float64},3}(undef, 0, 0, 0), Float64[], Float64[],
    Float64[])

"""
    evaluatehybrid!(B, C, block::ScatteringParameters, ws::AbstractVector,
        [work::HybridWorkspace])

The hybrid (wave to modified nodal analysis) coefficients `B` and `C` of a
scattering block at the signed frequencies `ws`.

The constitutive equation of the block is `B(w) phi' = C(w) i` with
`B = R^(-1/2)(I-S)` and `C = R^(1/2)(I+S)`; nothing is inverted, so a block
whose `I+S` is singular somewhere stamps exactly. At exactly zero frequency
the rows become `i = 0`, consistent with resistors in the node flux basis.

        C::AbstractArray{Complex{Float64},3}, block::ScatteringParameters,
        ws::AbstractVector)

Evaluate the coefficient matrices of the hybrid (wave to modified nodal
analysis) constitutive equations of `block` at the signed angular
frequencies `ws`:

    B(w) v - C(w) i = 0,  B = R^(-1/2) (I - S(w)),  C = R^(1/2) (I + S(w)),

with `R` the diagonal of the reference impedances, `v` the port voltages
and `i` the port currents (which the solvers carry as auxiliary variables).
This representation exists for every scattering matrix: unlike the
admittance `Y = R^(-1/2)(I-S)(I+S)^(-1)R^(-1/2)`, nothing is inverted, so
blocks whose `I+S` is singular somewhere (an ideal short `S = -1`, or a
lossless through line at each of its half wavelength resonances, where
`det(I+S) = 1-exp(-2*im*theta) = 0`) are stamped exactly.

The scattering parameters are evaluated with the negative frequency rule of
the block through [`evaluatescattering!`](@ref) at the native reference
impedances; no renormalization of the data is performed. At exactly zero
frequency the rows are replaced by `i = 0` (`B = 0`, `C = I`): the node
flux basis carries no DC voltage, and direct currents flow only through the
inductive branches of the static flux stiffness graph, so a scattering
block carries no direct current, consistent with the treatment of
resistors.

Pass a [`HybridWorkspace`](@ref) to reuse its scratch across calls.
"""
function evaluatehybrid!(B::AbstractArray{Complex{Float64},3},
    C::AbstractArray{Complex{Float64},3}, block::ScatteringParameters,
    ws::AbstractVector)
    return evaluatehybrid!(B, C, block, ws, HybridWorkspace())
end

function evaluatehybrid!(B::AbstractArray{Complex{Float64},3},
    C::AbstractArray{Complex{Float64},3}, block::ScatteringParameters,
    ws::AbstractVector, work::HybridWorkspace)

    n = block.nports
    if size(B) != (n, n, length(ws)) || size(C) != (n, n, length(ws))
        throw(DimensionMismatch(lazy"The destination arrays have sizes $(size(B)) and $(size(C)) but ($(n), $(n), $(length(ws))) is required."))
    end

    # evaluate the scattering parameters only at the nonzero frequencies
    empty!(work.nonzeroindices)
    empty!(work.wnz)
    for i in eachindex(ws)
        if !iszero(ws[i])
            push!(work.nonzeroindices, i)
            push!(work.wnz, ws[i])
        end
    end
    nonzeroindices = work.nonzeroindices
    k = length(nonzeroindices)
    if size(work.S) != (n, n, k)
        work.S = Array{Complex{Float64},3}(undef, n, n, k)
    end
    S = work.S
    evaluatescattering!(S, block, work.wnz, work.absws)

    resize!(work.rinv2, n)
    resize!(work.r2, n)
    rinv2 = work.rinv2
    r2 = work.r2
    @inbounds for p in 1:n
        r2[p] = sqrt(block.zref[p])
        rinv2[p] = 1/r2[p]
    end
    fill!(B, zero(Complex{Float64}))
    fill!(C, zero(Complex{Float64}))
    # zero frequency rows: i = 0
    for i in eachindex(ws)
        if iszero(ws[i])
            for p in 1:n
                C[p,p,i] = one(Complex{Float64})
            end
        end
    end
    for (kk, i) in enumerate(nonzeroindices)
        for q in 1:n
            for p in 1:n
                B[p,q,i] = -rinv2[p]*S[p,q,kk]
                C[p,q,i] = r2[p]*S[p,q,kk]
            end
        end
        for p in 1:n
            B[p,p,i] += rinv2[p]
            C[p,p,i] += r2[p]
        end
    end
    return B, C
end


# === the stamp system ===

# one scattering block with its port terminal nodes in the parsed circuit
# and the position of its auxiliary port current variables
struct StampedScatteringBlock
    block::Any            # the shared ScatteringParameters definition
    signalnodes::Vector{Int}
    refnodes::Vector{Int}
    auxbase::Int          # aux index of port p mode m: auxbase+(p-1)*Nmodes+m
    name::String          # component name of the first port, for messages
end

"""
    ScatteringStampSystem

The contribution of the [`ScatteringParameters`](@ref) components of a parsed
circuit to the harmonic balance system matrix, as hybrid (wave to modified
nodal analysis) stamps: one auxiliary port current variable per port and
mode, the constant Kirchhoff current law couplings `kcl` of those currents
into the node equations, and the constitutive equations

    im*w_m*scale*B(w_m) phi - C(w_m) i = 0

per port and mode, whose frequency dependent coefficient entries are
described by the sparsity `pattern` plus, per scalar contribution, the
block, port pair, sign, mode, destination index, and which coefficient
(`B`, stamped with the `im*w_m*scale` factor of a voltage in the node flux
basis, mirroring the constitutive equations of the promoted port resistors,
or `C`, stamped as `-C`). The contribution is diagonal in mode space
because a linear time invariant multiport cannot convert frequencies; all
mode coupling in the system comes from the junction pump modulation term.
This representation exists for every scattering matrix (see
[`evaluatehybrid!`](@ref)); no admittance conversion is performed.

The values side is a pure gather-add over precomputable per (block, mode)
coefficients, so the same structure extends to a GPU kernel; the current
implementation assembles on the CPU.
"""
struct ScatteringStampSystem
    blocks::Vector{StampedScatteringBlock}
    # the constant Kirchhoff current law couplings of the auxiliary port
    # currents into the node equations
    kcl::SparseMatrixCSC{Complex{Float64},Int}
    # the frequency dependent constitutive entries
    pattern::SparseMatrixCSC{Complex{Float64},Int}
    # per scalar contribution: destination nonzero index in the pattern and
    # in the target system matrix (set by setscatteringindexmap!)
    patternindex::Vector{Int}
    Aindex::Vector{Int}
    blockindex::Vector{Int32}
    pindex::Vector{Int32}
    qindex::Vector{Int32}
    # 1 for a B entry (times sign*im*w_m*scale), 2 for a C entry (times -1)
    coeff::Vector{Int8}
    sign::Vector{Int8}
    modeindex::Vector{Int32}
    Nmodes::Int
    Nauxports::Int
    scale::Float64
end

"""
    countscatteringports(psc::ParsedSortedCircuit)

The total number of scattering block ports (`:S` components) of the parsed
circuit, which is the number of auxiliary port current variables per mode.
"""
countscatteringports(psc::ParsedSortedCircuit) = count(==(:S),
    psc.componenttypes)

"""
    hasscattering(psc::ParsedSortedCircuit)

Whether the parsed circuit contains scattering block components.
"""
hasscattering(psc::ParsedSortedCircuit) = any(==(:S), psc.componenttypes)

"""
    scatteringstampsystem(psc::ParsedSortedCircuit, Nmodes::Integer;
        auxoffset::Integer, Ntotal::Integer, scale::Real = 1.0)

Collect the `:S` components of the parsed circuit into a
[`ScatteringStampSystem`](@ref), or return `nothing` when there are none.
Ports are grouped into blocks by the identity of the shared
[`ScatteringParameters`](@ref) definition, so instances sharing a definition
share its data; every port of each block must appear exactly once. The
auxiliary port current variables occupy the `countscatteringports(psc) *
Nmodes` indices starting after `auxoffset` of the `Ntotal` dimensional
system, and `scale` multiplies the stamped `B` entries (`Lmean` in the
scaled nonlinear solver, one in the linearized solver), mirroring the
constitutive equations of the promoted port resistors.
"""
function scatteringstampsystem(psc::ParsedSortedCircuit, Nmodes::Integer;
    auxoffset::Integer, Ntotal::Integer, scale::Real = 1.0)

    blockofdef = IdDict{Any,Int}()
    blocks = StampedScatteringBlock[]
    auxbase = auxoffset
    for (i, type) in enumerate(psc.componenttypes)
        type == :S || continue
        stamp = psc.componentvalues[i]
        if !(stamp isa ScatteringStamp)
            throw(ArgumentError(lazy"The component $(psc.componentnames[i]) has type :S but its value is a $(typeof(stamp)) rather than a ScatteringStamp."))
        end
        block = stamp.block
        bi = get!(blockofdef, block) do
            n = block.nports
            push!(blocks, StampedScatteringBlock(block, zeros(Int, n),
                zeros(Int, n), auxbase, psc.componentnames[i]))
            auxbase += n*Nmodes
            length(blocks)
        end
        sb = blocks[bi]
        p = stamp.port
        if !(1 <= p <= sb.block.nports)
            throw(ArgumentError(lazy"The component $(psc.componentnames[i]) references port $(p) of a $(sb.block.nports) port scattering block."))
        end
        if sb.signalnodes[p] != 0
            throw(ArgumentError(lazy"Port $(p) of the scattering block at $(sb.name) appears more than once in the parsed circuit."))
        end
        sb.signalnodes[p] = psc.nodeindices[1,i]
        sb.refnodes[p] = psc.nodeindices[2,i]
    end

    if isempty(blocks)
        return nothing
    end
    for sb in blocks
        missingports = [p for p in 1:sb.block.nports if sb.signalnodes[p] == 0]
        if !isempty(missingports)
            throw(ArgumentError(lazy"The ports $(missingports) of the scattering block at $(sb.name) are missing from the parsed circuit; every port of a block must lower to one :S component."))
        end
    end

    # the constant Kirchhoff current law couplings: the port current of
    # port p enters the signal node and leaves the reference node, exactly
    # as the branch currents of the promoted port resistors do
    kclrows = Int[]
    kclcols = Int[]
    kclvals = Complex{Float64}[]
    # the frequency dependent constitutive entries, with their metadata.
    # sparse() sums duplicate positions, so contributions map onto the
    # pattern afterwards.
    rows = Int[]
    cols = Int[]
    blockindex = Int32[]
    pindex = Int32[]
    qindex = Int32[]
    coeff = Int8[]
    sign = Int8[]
    modeindex = Int32[]
    for (bi, sb) in enumerate(blocks)
        n = sb.block.nports
        for p in 1:n
            for m in 1:Nmodes
                auxp = sb.auxbase + (p-1)*Nmodes + m
                if sb.signalnodes[p] != 1
                    noderow = (sb.signalnodes[p]-2)*Nmodes + m
                    push!(kclrows, noderow); push!(kclcols, auxp)
                    push!(kclvals, 1)
                end
                if sb.refnodes[p] != 1
                    noderow = (sb.refnodes[p]-2)*Nmodes + m
                    push!(kclrows, noderow); push!(kclcols, auxp)
                    push!(kclvals, -1)
                end
                for q in 1:n
                    auxq = sb.auxbase + (q-1)*Nmodes + m
                    # -C[p,q] on the auxiliary current columns
                    push!(rows, auxp); push!(cols, auxq)
                    push!(blockindex, bi); push!(pindex, p); push!(qindex, q)
                    push!(coeff, 2); push!(sign, 1); push!(modeindex, m)
                    # im*w*scale*B[p,q] times the node flux of the signal
                    # terminal minus that of the reference terminal
                    for (node, s) in ((sb.signalnodes[q], Int8(1)),
                                      (sb.refnodes[q], Int8(-1)))
                        node == 1 && continue # grounded terminal: phi = 0
                        push!(rows, auxp)
                        push!(cols, (node-2)*Nmodes + m)
                        push!(blockindex, bi); push!(pindex, p)
                        push!(qindex, q); push!(coeff, 1); push!(sign, s)
                        push!(modeindex, m)
                    end
                end
            end
        end
    end

    kcl = sparse(kclrows, kclcols, kclvals, Ntotal, Ntotal)
    pattern = sparse(rows, cols, ones(Complex{Float64}, length(rows)),
        Ntotal, Ntotal)

    # map each contribution to its nonzero index in the pattern
    patternindex = Vector{Int}(undef, length(rows))
    for c in eachindex(rows)
        col = cols[c]
        r = pattern.colptr[col]:(pattern.colptr[col+1]-1)
        k = searchsortedfirst(view(pattern.rowval, r), rows[c])
        patternindex[c] = first(r) + k - 1
    end

    return ScatteringStampSystem(blocks, kcl, pattern, patternindex,
        copy(patternindex), blockindex, pindex, qindex, coeff, sign,
        modeindex, Int(Nmodes), auxbase - auxoffset, Float64(scale))
end

"""
    setscatteringindexmap!(ssys::ScatteringStampSystem,
        A::SparseMatrixCSC)

Point the destination indices of the frequency dependent contributions at
the nonzero values of the system matrix `A`, whose sparsity structure must
contain the pattern of `ssys` (ensure this by merging `ssys.pattern` into
the structure before calling, as [`HBLinearizedSystem`](@ref) does; the
constant `ssys.kcl` couplings are folded into the constant augmentation
matrix instead).
"""
function setscatteringindexmap!(ssys::ScatteringStampSystem,
    A::SparseMatrixCSC)
    indexmap = sparseaddmap(A, ssys.pattern)
    for c in eachindex(ssys.patternindex)
        ssys.Aindex[c] = indexmap[ssys.patternindex[c]]
    end
    return ssys
end

"""
    ScatteringWorkspace()

Reusable scratch for [`scatteringvalues!`](@ref): a [`HybridWorkspace`](@ref)
and the coefficient arrays of every block.

Those arrays have to be live at once, because the contributions are ordered
by their destination in the system matrix rather than by block, so they
cannot be one buffer reused block by block. They can be reused from one
frequency to the next, which is what this is for: a sweep over a line whose
every cell is its own block otherwise allocates two arrays per block at every
signal frequency.

A workspace is written by every call, so it belongs to one task. A caller
which sweeps on several threads gives each one its own.
"""
mutable struct ScatteringWorkspace
    hybrid::HybridWorkspace
    Bs::Vector{Array{Complex{Float64},3}}
    Cs::Vector{Array{Complex{Float64},3}}
end

ScatteringWorkspace() = ScatteringWorkspace(HybridWorkspace(),
    Array{Complex{Float64},3}[], Array{Complex{Float64},3}[])

# grow the caches to the shapes this stamp system and mode count need,
# reallocating only what has the wrong shape
function preparecoefficients!(work::ScatteringWorkspace,
    ssys::ScatteringStampSystem, nw::Integer)
    nb = length(ssys.blocks)
    if length(work.Bs) != nb
        resize!(work.Bs, nb)
        resize!(work.Cs, nb)
        for bi in 1:nb
            n = ssys.blocks[bi].block.nports
            work.Bs[bi] = Array{Complex{Float64},3}(undef, n, n, nw)
            work.Cs[bi] = Array{Complex{Float64},3}(undef, n, n, nw)
        end
    elseif nb > 0 && size(work.Bs[1], 3) != nw
        for bi in 1:nb
            n = ssys.blocks[bi].block.nports
            work.Bs[bi] = Array{Complex{Float64},3}(undef, n, n, nw)
            work.Cs[bi] = Array{Complex{Float64},3}(undef, n, n, nw)
        end
    end
    return work
end

"""
    assemblescattering!(A::SparseMatrixCSC, ssys::ScatteringStampSystem,
        wmodes::AbstractVector)

Add the frequency dependent constitutive entries of the scattering blocks
at the signed mode frequencies `wmodes` into the values of the system
matrix `A`: `sign*im*w_m*scale*B[p,q](w_m)` on the node flux columns and
`-C[p,q](w_m)` on the auxiliary current columns, with the coefficients
evaluated per block by [`evaluatehybrid!`](@ref) (which applies the
negative frequency rule of each block, so the negative mode entries carry
the complex conjugate data exactly as the conjugation of the conductance
matrix does for resistors). The destination indices must have been set with
[`setscatteringindexmap!`](@ref). Thread safe: evaluation buffers are local
to the call and blocks are read only.
"""
function assemblescattering!(A::SparseMatrixCSC,
    ssys::ScatteringStampSystem, wmodes::AbstractVector,
    work::ScatteringWorkspace = ScatteringWorkspace())

    values = Vector{Complex{Float64}}(undef, length(ssys.Aindex))
    scatteringvalues!(values, ssys, wmodes, work)
    nzval = A.nzval
    @inbounds for c in eachindex(ssys.Aindex)
        nzval[ssys.Aindex[c]] += values[c]
    end
    return A
end

"""
    scatteringvalues!(values::AbstractVector, ssys::ScatteringStampSystem,
        wmodes::AbstractVector, [work::ScatteringWorkspace])

The value each scalar contribution of the scattering blocks adds to the
system matrix at the signed mode frequencies `wmodes`, in the order of
`ssys.Aindex`.

This is the half of [`assemblescattering!`](@ref) which has to run on the
host, because it evaluates each block through its provider, which may be an
arbitrary callable or an interpolation of tabulated data. What is left is a
gather-add of these values into the stored entries, which is where the two
backends part company: the host adds them here, and a device adds them with
a kernel (see [`DeviceScatteringStamps`](@ref)). Splitting it this way is
what lets both backends stamp identical values.
"""
function scatteringvalues!(values::AbstractVector,
    ssys::ScatteringStampSystem, wmodes::AbstractVector)
    return scatteringvalues!(values, ssys, wmodes, ScatteringWorkspace())
end

function scatteringvalues!(values::AbstractVector,
    ssys::ScatteringStampSystem, wmodes::AbstractVector,
    work::ScatteringWorkspace)

    if length(wmodes) != ssys.Nmodes
        throw(DimensionMismatch(lazy"scatteringvalues! received $(length(wmodes)) mode frequencies but the stamp system was built for $(ssys.Nmodes) modes."))
    end
    if length(values) != length(ssys.Aindex)
        throw(DimensionMismatch(lazy"`values` has length $(length(values)) but the stamp system has $(length(ssys.Aindex)) contributions."))
    end

    preparecoefficients!(work, ssys, length(wmodes))
    Bs = work.Bs
    Cs = work.Cs
    for (bi, sb) in enumerate(ssys.blocks)
        evaluatehybrid!(Bs[bi], Cs[bi], sb.block, wmodes, work.hybrid)
    end

    @inbounds for c in eachindex(ssys.Aindex)
        m = ssys.modeindex[c]
        bi = ssys.blockindex[c]
        p = ssys.pindex[c]
        q = ssys.qindex[c]
        values[c] = if ssys.coeff[c] == 1
            ssys.sign[c] * (im*wmodes[m]*ssys.scale) * Bs[bi][p, q, m]
        else
            -Cs[bi][p, q, m]
        end
    end
    return values
end

"""
    scatteringlinearterm(psc::ParsedSortedCircuit, wmodes::AbstractVector,
        Nmodes::Integer; auxoffset::Integer, Ntotal::Integer,
        scale::Real = 1.0)

The constant sparse matrix of the scattering block contribution at the
fixed mode frequencies `wmodes` (the constitutive entries plus the
Kirchhoff current law couplings of the auxiliary port currents), or
`nothing` when the circuit has no scattering blocks. Used by the nonlinear
(pump) solver, where the mode frequencies do not change: the contribution
is folded into the frequency independent linear term alongside the
augmentation matrix of the promoted resistors, so the residual, Jacobian,
and solver machinery operate on the augmented system unchanged.
"""
function scatteringlinearterm(psc::ParsedSortedCircuit,
    wmodes::AbstractVector, Nmodes::Integer; auxoffset::Integer,
    Ntotal::Integer, scale::Real = 1.0)

    ssys = scatteringstampsystem(psc, Nmodes; auxoffset = auxoffset,
        Ntotal = Ntotal, scale = scale)
    if isnothing(ssys)
        return nothing
    end
    Snm = copy(ssys.pattern)
    fill!(Snm.nzval, zero(Complex{Float64}))
    # the pattern indices are the destination indices for the pattern itself
    copyto!(ssys.Aindex, ssys.patternindex)
    assemblescattering!(Snm, ssys, wmodes)
    return spaddkeepzeros(Snm, ssys.kcl)
end

# === the vacuum noise of dissipative blocks ===

"""
    ScatteringNoisePlan

The vacuum noise channels of the dissipative [`ScatteringParameters`](@ref)
components of a circuit: which blocks of a [`ScatteringStampSystem`](@ref)
carry noise and where their channels sit in the rows of the noise
scattering matrix.

A block which absorbs must add noise, or its output would violate the
commutation relations. In the wave domain its constitutive equation is
`b = S a + n` with an added noise wave whose vacuum covariance is
`I - S S'`, and in the hybrid stamp

    im*w_m*scale*B(w_m) phi - C(w_m) i = 2 n

that noise is a source in the auxiliary port current rows. A block with `n`
ports therefore carries `n` noise channels, one per eigenvector of
`I - S S'`; those of a lossless block are identically zero, so only blocks
which are not [`provablylossless`](@ref) are given channels.

The channels of the blocks follow the noise ports of the dissipative
lumped components in the rows of `Snoise`, with the same
channel-major-mode-minor ordering.
"""
struct ScatteringNoisePlan
    # index into the blocks of the stamp system, and the 0 based channel
    # offset of each, in channel (not row) units
    blockindices::Vector{Int}
    channelbase::Vector{Int}
    Nchannels::Int
    Nmodes::Int
end

"""
    planscatteringnoise(ssys)

The [`ScatteringNoisePlan`](@ref) of a [`ScatteringStampSystem`](@ref), or
`nothing` when the circuit has no scattering blocks or all of them are
[`provablylossless`](@ref).

The noise models which carry a covariance the analysis derives from the
scattering data are supported: [`Passive`](@ref) and
[`ThermalEquilibrium`](@ref), which differ only in the temperature the
channels are at, and [`Lossless`](@ref), which declares there are none. A
[`NoiseCovariance`](@ref) gives the covariance outright and is not lowered
into a circuit.
"""
planscatteringnoise(::Nothing) = nothing
function planscatteringnoise(ssys::ScatteringStampSystem)
    blockindices = Int[]
    channelbase = Int[]
    nch = 0
    for (bi, sb) in enumerate(ssys.blocks)
        checknoisemodel(sb.block, sb.name)
        # a block asserted lossless, or shown to be, adds nothing
        sb.block.noise isa Lossless && continue
        provablylossless(sb.block) && continue
        push!(blockindices, bi)
        push!(channelbase, nch)
        nch += sb.block.nports
    end
    isempty(blockindices) && return nothing
    return ScatteringNoisePlan(blockindices, channelbase, nch, ssys.Nmodes)
end

# the noise models other than Passive() are part of the component API but
# not yet of the noise calculation; silently ignoring one would return a
# quantum efficiency computed from a different block than the user asked
# for, so they are an error at the point the noise is planned rather than a
# warning
function checknoisemodel(block::ScatteringParameters, name)
    noise = block.noise
    if noise isa Passive || noise isa Lossless ||
            noise isa ThermalEquilibrium
        return nothing
    end
    throw(ArgumentError(lazy"The scattering block at $(name) has the noise model $(noise), which the noise and quantum efficiency calculations do not yet support. Use Passive(), Lossless(), or ThermalEquilibrium(T), or request no noise outputs."))
end
checknoisemodel(block, name) = nothing

"""
    scatteringnoisenames(plan::ScatteringNoisePlan,
        ssys::ScatteringStampSystem)

The name of each noise channel of the plan, for labelling the rows of a
keyed noise scattering matrix. A block with more than one port has one
channel per port, distinguished by a channel number, because the
eigenvectors of `I - S S'` mix the ports and no single port owns a channel.
"""
function scatteringnoisenames(plan::ScatteringNoisePlan,
    ssys::ScatteringStampSystem)
    names = Vector{String}(undef, plan.Nchannels)
    for (e, bi) in enumerate(plan.blockindices)
        sb = ssys.blocks[bi]
        n = sb.block.nports
        for c in 1:n
            names[plan.channelbase[e]+c] = n == 1 ? sb.name :
                string(sb.name, "#", c)
        end
    end
    return names
end

"""
    noisecovariance!(L, off, n, S, soff)

The vacuum noise covariance `I - S S'` of an `n` port block, into the length
`n*n` column major block of `L` at `off`, from the scattering matrix in the
same layout at `soff` in `S`.

Only the lower triangle is written, which is all
[`psdcholesky!`](@ref) reads: the covariance is Hermitian.
"""
@inline function noisecovariance!(L, off::Integer, n::Integer, S, soff::Integer)
    T = eltype(L)
    @inbounds for c in 1:n
        for p in c:n
            acc = p == c ? one(T) : zero(T)
            for l in 1:n
                acc -= S[soff + (l-1)*n + p]*conj(S[soff + (l-1)*n + c])
            end
            L[off + (c-1)*n + p] = acc
        end
    end
    return nothing
end

"""
    psdcholesky!(L, off, n)

Overwrite the length `n*n` column major block of `L` at `off`, whose lower
triangle holds a Hermitian positive semidefinite matrix `V`, with a lower
triangular factor satisfying `L L' = V`.

The noise covariance of a block which is transparent in some direction, such
as a series element, is singular, and rounding can make that of a lossless
block slightly indefinite, so this is not a plain Cholesky factorization: a
pivot which is not positive is taken as zero, which for a positive
semidefinite matrix means its whole column is zero and there is nothing to
divide by. The strict upper triangle held the covariance and is zeroed, so
the result is the factor and nothing else.

Any factor of the covariance describes the same noise: the channels are
defined only up to a unitary mixing among them, and the quantum efficiency
and the commutation relations read only sums over them. A triangular factor
is chosen over an eigendecomposition because it is a few lines of arithmetic
with no iteration, so the host and a kernel can share it and agree to the
last bit.
"""
@inline function psdcholesky!(L, off::Integer, n::Integer)
    @inbounds for c in 1:n
        d = real(L[off + (c-1)*n + c])
        for l in 1:c-1
            d -= abs2(L[off + (l-1)*n + c])
        end
        if d > 0
            r = sqrt(d)
            L[off + (c-1)*n + c] = r
            for p in c+1:n
                acc = L[off + (c-1)*n + p]
                for l in 1:c-1
                    acc -= L[off + (l-1)*n + p]*conj(L[off + (l-1)*n + c])
                end
                L[off + (c-1)*n + p] = acc/r
            end
        else
            for p in c:n
                L[off + (c-1)*n + p] = 0
            end
        end
    end
    @inbounds for c in 2:n
        for p in 1:c-1
            L[off + (c-1)*n + p] = 0
        end
    end
    return nothing
end

"""
    checklosslessblocks(ssys::ScatteringStampSystem, w, wpumpmodes;
        atol = 1e-6, nsamples = 32)

Check the blocks of `ssys` which declare [`Lossless`](@ref) against the
frequencies the sweep will solve at, and throw if one of them dissipates
there.

A block whose data is stored is held to the declaration when it is
constructed. A callable cannot be, which is what the declaration is for, so
this is the one check available: evaluate it at up to `nsamples` of the
signal frequencies, at every pump mode, and see. Sampling can show that a
block dissipates and can never show that it does not, so this catches a
declaration which is wrong at a frequency it looked at and makes no promise
about the rest. It is bounded rather than exhaustive because an exhaustive
pass would cost a share of what the declaration saves.
"""
function checklosslessblocks(ssys::Union{Nothing,ScatteringStampSystem}, w,
    wpumpmodes; atol::Real = 1e-6, nsamples::Integer = 32)

    isnothing(ssys) && return nothing
    any(sb -> sb.block.noise isa Lossless &&
        !provablylossless(sb.block), ssys.blocks) || return nothing
    nw = length(w)
    step = max(1, cld(nw, nsamples))
    ws = Float64[]
    for i in 1:step:nw
        for m in wpumpmodes
            push!(ws, w[i] + m)
        end
    end
    isempty(ws) && return nothing
    work = HybridWorkspace()
    for sb in ssys.blocks
        (sb.block.noise isa Lossless && !provablylossless(sb.block)) || continue
        checkoneblocklossless(sb.block, sb.name, ws, atol, work)
    end
    return nothing
end

# behind a function barrier: the block is stored untyped
function checkoneblocklossless(block::ScatteringParameters, name, ws, atol, work)
    n = block.nports
    S = Array{Complex{Float64},3}(undef, n, n, length(ws))
    evaluatescattering!(S, block, ws, work.absws)
    worst = 0.0
    worstw = 0.0
    for k in eachindex(ws)
        iszero(ws[k]) && continue
        d = unitaritydeviation(view(S, :, :, k))
        if d > worst
            worst = d
            worstw = ws[k]
        end
    end
    if worst > atol
        throw(ArgumentError(lazy"the scattering block at $(name) declares noise = Lossless(), but at $(worstw) rad/s the largest absolute entry of I - S*S' is $(worst). A block which dissipates must carry the noise its loss requires; use the default Passive() noise model."))
    end
    return nothing
end

"""
    ScatteringNoiseWorkspace()

The scratch of [`scatteringnoisewaves!`](@ref): the scattering parameters
at the mode frequencies, the noise covariance and its factor, the port
Norton currents, and the buffer of unsigned frequencies. Reused across
frequencies by one worker.
"""
mutable struct ScatteringNoiseWorkspace
    S::Array{Complex{Float64},3}
    # the covariance and its factor, flat and column major, so that the
    # factorization is the one the kernel runs
    L::Vector{Complex{Float64}}
    absws::Vector{Float64}
end

function ScatteringNoiseWorkspace()
    return ScatteringNoiseWorkspace(
        Array{Complex{Float64},3}(undef, 0, 0, 0),
        Complex{Float64}[], Float64[])
end

# the node flux of a terminal in the adjoint solution, which is zero at
# ground because ground carries no row
@inline function terminalflux(phi::AbstractMatrix, node::Integer,
    m::Integer, Nmodes::Integer, k::Integer)
    node == 1 && return zero(eltype(phi))
    return @inbounds phi[(node-2)*Nmodes + m, k]
end

"""
    scatteringnoisewaves!(noiseoutputwave, plan::ScatteringNoisePlan,
        ssys::ScatteringStampSystem, phiadj, wmodes, rowoffset,
        work = ScatteringNoiseWorkspace())

Write the noise output waves of the scattering block channels of `plan`
into the rows of `noiseoutputwave` after `rowoffset`, from the adjoint
solution `phiadj` at the mode frequencies `wmodes`.

The noise wave `n` of a block enters its constitutive equation as a source
in the auxiliary port current rows, so by the adjoint identity its
contribution to the output is that source contracted against those same rows
of the adjoint solution, weighted by the factor `L` of the vacuum covariance
`L L' = I - S S'`:

    noiseoutputwave[channel c] = sqrt(abs(w)) sum_p L[p,c] i[p]

`phiadj` must be the solution of the *transposed* system, which is what
[`hblinsolve`](@ref) solves for its adjoint. The conjugated pump system,
which is the same matrix for a circuit without blocks, is not: it agrees with
the transposed system in the node flux rows and not in the auxiliary port
current rows, and on a non reciprocal block the two differ by the direction
the block transmits in, so contracting it would give a block which emits its
noise backwards.

The channel of a mode whose frequency is zero is zero, matching the wave
normalization of the lumped noise ports, which is singular there.
"""
function scatteringnoisewaves!(noiseoutputwave::AbstractMatrix,
    plan::ScatteringNoisePlan, ssys::ScatteringStampSystem,
    phiadj::AbstractMatrix, wmodes::AbstractVector, rowoffset::Integer,
    work::ScatteringNoiseWorkspace = ScatteringNoiseWorkspace())

    Nmodes = plan.Nmodes
    if length(wmodes) != Nmodes
        throw(DimensionMismatch(lazy"scatteringnoisewaves! received $(length(wmodes)) mode frequencies but the plan was built for $(Nmodes) modes."))
    end
    if size(noiseoutputwave,1) < rowoffset + plan.Nchannels*Nmodes
        throw(DimensionMismatch(lazy"`noiseoutputwave` has $(size(noiseoutputwave,1)) rows but the scattering block channels need $(rowoffset + plan.Nchannels*Nmodes)."))
    end
    for (e, bi) in enumerate(plan.blockindices)
        sb = ssys.blocks[bi]
        # the block is stored untyped, so the per port and per mode work
        # goes behind a function barrier
        blocknoisewaves!(noiseoutputwave, sb.block, sb, wmodes, phiadj,
            rowoffset + plan.channelbase[e]*Nmodes, Nmodes, work)
    end
    return noiseoutputwave
end

function blocknoisewaves!(noiseoutputwave::AbstractMatrix,
    block::ScatteringParameters, sb::StampedScatteringBlock,
    wmodes::AbstractVector, phiadj::AbstractMatrix, rowoffset::Integer,
    Nmodes::Integer, work::ScatteringNoiseWorkspace)

    n = block.nports
    nrhs = size(phiadj, 2)
    if size(work.S) != (n, n, Nmodes)
        work.S = Array{Complex{Float64},3}(undef, n, n, Nmodes)
        resize!(work.L, n*n)
    end
    S = work.S
    L = work.L
    evaluatescattering!(S, block, wmodes, work.absws)
    for m in 1:Nmodes
        if iszero(wmodes[m])
            for c in 1:n
                for k in 1:nrhs
                    noiseoutputwave[rowoffset + (c-1)*Nmodes + m, k] = 0
                end
            end
            continue
        end
        # the vacuum covariance of the added noise wave and its factor,
        # `L L' = I - S S'`, in the flat layout the kernel of the device
        # path uses so that the two compute the same channels
        noisecovariance!(L, 0, n, S, (m-1)*n*n)
        psdcholesky!(L, 0, n)
        # the noise wave enters the constitutive equation as the source
        # `2n` of the auxiliary port current rows, in the sqrt(power)
        # normalization of the hybrid stamp; `sqrt(abs(w))` carries it to
        # the sqrt(photons/second) normalization the outputs are in
        kw = sqrt(abs(wmodes[m]))
        @inbounds for k in 1:nrhs
            for c in 1:n
                acc = zero(Complex{Float64})
                for p in 1:n
                    acc += L[(c-1)*n + p]*
                        phiadj[sb.auxbase + (p-1)*Nmodes + m, k]
                end
                noiseoutputwave[rowoffset + (c-1)*Nmodes + m, k] = kw*acc
            end
        end
    end
    return noiseoutputwave
end

"""
    noisechanneltemperatures(psc::ParsedSortedCircuit,
        noiseportimpedanceindices, noiseplan, ssys, temperature)

The temperature of each row of the noise scattering matrix, in the order
[`noisechannelnames`](@ref) gives them.

`temperature` is the analysis default, which every dissipative element takes
unless it states one of its own. A lumped component states it as
`Resistor(R; temperature = T)` and a [`ScatteringParameters`](@ref) as
`noise = ThermalEquilibrium(T)`, both of which are recorded by
[`parsesortcircuit`](@ref) as it lowers the circuit. Only the typed circuit
format carries them; a netlist of tuples states none and everything in it
takes the default.

A block with a [`NoiseCovariance`](@ref) is left alone entirely: its
covariance is given rather than derived, and scaling something the caller
supplied by a temperature would be modifying an answer they already know.
"""
function noisechanneltemperatures(psc, noiseportimpedanceindices, noiseplan,
    ssys, temperature)

    stated = psc.componenttemperatures
    ts = Float64[get(stated, i, Float64(temperature))
        for i in noiseportimpedanceindices]
    isnothing(noiseplan) && return ts
    # a block's ports lower to consecutive :S components, and its noise model
    # was recorded against the first of them
    index = psc.componentnamedict
    for (e, bi) in enumerate(noiseplan.blockindices)
        sb = ssys.blocks[bi]
        i = get(index, sb.name, 0)
        t = iszero(i) ? Float64(temperature) :
            get(stated, i, Float64(temperature))
        for _ in 1:sb.block.nports
            push!(ts, t)
        end
    end
    return ts
end

"""
    noisechannelnames(componentnames, noiseportimpedanceindices, noiseplan,
        ssys)

The name of each row of the noise scattering matrix: the dissipative lumped
components first, then the channels of the dissipative scattering blocks.
"""
function noisechannelnames(componentnames, noiseportimpedanceindices,
    noiseplan, ssys)
    names = String[String(componentnames[i]) for i in noiseportimpedanceindices]
    isnothing(noiseplan) && return names
    return vcat(names, scatteringnoisenames(noiseplan, ssys))
end
