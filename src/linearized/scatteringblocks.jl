# === hybrid wave to modified nodal analysis coefficients ===

"""
    HybridWorkspace()

Reusable scratch for [`evaluatehybrid!`](@ref).

Every call needs the same six small buffers: which frequencies are nonzero,
those frequencies, the scattering parameters there, the two square roots of
the reference impedances, and the unsigned frequencies the providers are
evaluated at. A line whose every cell is its own block
evaluates thousands of blocks at every signal frequency, so allocating those
per call dominated the evaluation. They are reallocated only when a block
needs a larger one, which for a circuit whose blocks are all the same size is
once.

A workspace is mutable and is written by every call, so it belongs to one
task. [`scatteringvalues!`](@ref) makes one per call when it is not given a
workspace, which is what keeps it safe to call from several threads at
once; a caller which passes one gives each thread its own.
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
        work::HybridWorkspace)

Evaluate the coefficient matrices of the hybrid (wave to modified nodal
analysis) constitutive equations of `block` at the signed angular
frequencies `ws`:

    B(w) v - C(w) i = 0,  B = (I - S(w)) R^(-1/2),  C = (I + S(w)) R^(1/2),

with `R` the diagonal of the reference impedances, `v` the port voltages
and `i` the port currents (which the solvers carry as auxiliary variables).
The diagonal multiplies on the right: with the power waves
`a = (R^(-1/2) v + R^(1/2) i)/2` and `b = (R^(-1/2) v - R^(1/2) i)/2` and
`b = S a`, the entry `(p, q)` of `I - S` scales with the impedance of port
`q`, the port whose voltage it multiplies, so a block whose ports have
different reference impedances is not symmetric between the two sides.
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

`work` is the [`HybridWorkspace`](@ref) whose scratch is reused across
calls.
"""
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
                B[p,q,i] = -rinv2[q]*S[p,q,kk]
                C[p,q,i] = r2[q]*S[p,q,kk]
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

"""
    StampedScatteringBlock

One scattering block instance, with the nodes its ports attach to and the
position of the auxiliary port current variables which carry its currents.

`auxbase` locates them: port `p` at mode `m` is the state index
`auxbase + (p-1)*Nmodes + m`. That is the same current the zero frequency
row treats as an unknown, which is how a block becomes visible at direct
current; see [`DCBlockRows`](@ref).
"""
struct StampedScatteringBlock
    block::Any            # the shared ScatteringParameters definition
    signalnodes::Vector{Int}
    refnodes::Vector{Int}
    auxbase::Int          # aux index of port p mode m: auxbase+(p-1)*Nmodes+m
    name::String          # the instance path with a "/port1" suffix, for
                          # messages and `scatteringblockindex`
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
coefficients, so the values are computed here and a device adds them with
a kernel; see [`DeviceScatteringStamps`](@ref).

# Fields
- `blocks`: the [`StampedScatteringBlock`](@ref)s.
- `kcl`: the constant Kirchhoff current law couplings of the auxiliary port
    currents, folded into the modified nodal analysis augmentation.
- `pattern`: the sparsity pattern of the constitutive rows.
- `patternindex`, `Aindex`: where each stored entry of `pattern` lands in
    the system matrix, and the stored position it is added at.
- `blockindex`, `pindex`, `qindex`, `modeindex`: for each contribution, the
    block, the port pair and the mode it reads its coefficient from.
- `coeff`, `sign`: the per contribution coefficient buffer and sign.
- `Nmodes`, `Nauxports`: the mode count and the number of auxiliary port
    current unknowns.
- `scale`: the solver scale the rows are written in.
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
    # the input mode of each contribution, the output mode's own except
    # for a pumped block, and which pumped block a contribution belongs
    # to (its ordinal among `pumped`, zero for a block which does not
    # couple modes)
    inmodeindex::Vector{Int32}
    coupled::Vector{Int32}
    # the pumped blocks, as indices into `blocks`, the harmonic each
    # couples each pair of modes by (see [`pumpedharmonics`](@ref)), and
    # the mode frequency offsets those were read from
    pumped::Vector{Int}
    pumpedk::Vector{Matrix{Int}}
    modeoffsets::Vector{Float64}
    Nmodes::Int
    Nauxports::Int
    scale::Float64
end

"""
    countscatteringports(psc::CompiledCircuit)

The total number of scattering block ports of the circuit, which is the
number of auxiliary port current variables per mode.
"""
countscatteringports(psc::CompiledCircuit) =
    sum(b -> b.definition.nports, psc.scatteringblocks; init = 0)

"""
    scatteringstampsystem(blocks::Vector{StampedScatteringBlock}, Nmodes,
        Ntotal, scale; modeoffsets = nothing)

The positional inner form of [`scatteringstampsystem`](@ref), everything
past the point where the blocks and their terminals are known: the
Kirchhoff current law couplings, the constitutive pattern and the
contribution tables of a [`ScatteringStampSystem`](@ref). Unlike the
keyword form it never returns `nothing`.

`modeoffsets` are the frequency offsets of the modes from the signal, in
radians per second. A [`LinearizedScattering`](@ref) block couples the modes
whose offsets differ by a harmonic of its pump, so with one present the
offsets are required, and the contributions between such pairs are
entered on top of the diagonal ones; a block which does not convert
takes no notice of them.
"""
function scatteringstampsystem(blocks::Vector{StampedScatteringBlock},
    Nmodes::Integer, Ntotal::Integer, scale::Real; modeoffsets = nothing)

    # how many auxiliary port currents the blocks occupy in total: each has
    # one per port per mode, laid out consecutively from its own base
    naux = sum(sb -> sb.block.nports*Nmodes, blocks)

    # the pumped blocks and the harmonic by which each couples each pair
    # of modes
    pumped = Int[]
    pumpedk = Matrix{Int}[]
    offsets = isnothing(modeoffsets) ? zeros(Float64, Nmodes) :
        Float64.(collect(modeoffsets))
    length(offsets) == Nmodes || throw(DimensionMismatch(lazy"$(length(offsets)) mode offsets were given for $(Nmodes) modes."))
    for (bi, sb) in enumerate(blocks)
        sb.block isa LinearizedScattering || continue
        isnothing(modeoffsets) && throw(ArgumentError(lazy"the pumped block at $(sb.name) couples the modes of the solve, whose frequency offsets this stamp was not given."))
        push!(pumped, bi)
        push!(pumpedk, pumpedharmonics(sb.block, offsets))
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
    inmodeindex = Int32[]
    coupled = Int32[]
    # one contribution: the row of output port p at mode m, and the
    # columns of input port q at mode nn, its auxiliary current and the
    # node fluxes of its two terminals
    function contribute!(bi, sb, p, m, q, nn, j)
        n = sb.block.nports
        auxp = sb.auxbase + (p-1)*Nmodes + m
        auxq = sb.auxbase + (q-1)*Nmodes + nn
        # -C[p,q] on the auxiliary current columns
        push!(rows, auxp); push!(cols, auxq)
        push!(blockindex, bi); push!(pindex, p); push!(qindex, q)
        push!(coeff, 2); push!(sign, 1); push!(modeindex, m)
        push!(inmodeindex, nn); push!(coupled, j)
        # im*w*scale*B[p,q] times the node flux of the signal
        # terminal minus that of the reference terminal
        for (node, s) in ((sb.signalnodes[q], Int8(1)),
                          (sb.refnodes[q], Int8(-1)))
            node == 1 && continue # grounded terminal: phi = 0
            push!(rows, auxp)
            push!(cols, (node-2)*Nmodes + nn)
            push!(blockindex, bi); push!(pindex, p)
            push!(qindex, q); push!(coeff, 1); push!(sign, s)
            push!(modeindex, m); push!(inmodeindex, nn); push!(coupled, j)
        end
        return nothing
    end
    for (bi, sb) in enumerate(blocks)
        n = sb.block.nports
        j = something(findfirst(==(bi), pumped), 0)
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
                    contribute!(bi, sb, p, m, q, m, j)
                end
                j == 0 && continue
                # a pumped block: the columns of every other mode a
                # harmonic of its pump away
                K = pumpedk[j]
                for nn in 1:Nmodes
                    nn == m && continue
                    K[m, nn] == typemin(Int) && continue
                    for q in 1:n
                        contribute!(bi, sb, p, m, q, nn, j)
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
        modeindex, inmodeindex, coupled, pumped, pumpedk, offsets,
        Int(Nmodes), naux, Float64(scale))
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
    # the coefficient arrays of the pumped blocks, over (output port,
    # input port, output mode, input mode), and the harmonic transfer
    # functions at the mode frequencies and at their negatives
    Bp::Vector{Array{Complex{Float64},4}}
    Cp::Vector{Array{Complex{Float64},4}}
    Hpos::Array{Complex{Float64},4}
    Hneg::Array{Complex{Float64},4}
    negws::Vector{Float64}
    # the value of every contribution at one frequency
    values::Vector{Complex{Float64}}
end

ScatteringWorkspace() = ScatteringWorkspace(HybridWorkspace(),
    Array{Complex{Float64},3}[], Array{Complex{Float64},3}[],
    Array{Complex{Float64},4}[], Array{Complex{Float64},4}[],
    Array{Complex{Float64},4}(undef, 0, 0, 0, 0),
    Array{Complex{Float64},4}(undef, 0, 0, 0, 0), Float64[],
    Complex{Float64}[])

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
    np = length(ssys.pumped)
    if length(work.Bp) != np || (np > 0 && size(work.Bp[1], 3) != nw)
        resize!(work.Bp, np)
        resize!(work.Cp, np)
        for j in 1:np
            n = ssys.blocks[ssys.pumped[j]].block.nports
            work.Bp[j] = Array{Complex{Float64},4}(undef, n, n, nw, nw)
            work.Cp[j] = Array{Complex{Float64},4}(undef, n, n, nw, nw)
        end
    end
    return work
end

# The entry of a pumped block's conversion from port `q` at the input
# frequency of index `nn` to port `p` through the harmonic `k`, whose
# magnitude is the block's harmonic of index `j`: `H_k` at the input
# frequency for `k >= 0`, and for `k < 0` the conjugate of `H_{-k}` at the
# negative of the input frequency, since `H_{-k}(nu) = conj(H_k(-nu))`.
# `Hin` holds the harmonics at the input frequencies and `Hnegin` at
# their negatives, as `(p, q, j, nn)`.
conversionentry(Hin, Hnegin, k, j, p, q, nn) = k >= 0 ? Hin[p, q, j, nn] : conj(Hnegin[p, q, j, nn])

"""
    checkpumpconjugates(blocks::Vector{StampedScatteringBlock}, wmodes)

Refuse a pumped block of `blocks` which converts, with an entry above its
tolerance, the conjugate of one retained mode of a pump solve into
another, the two frequencies of `wmodes` summing to a harmonic of its
pump: the pump solve's linear term is complex linear in the retained
modes and does not carry a term in their conjugates, so with a source
the operating point would be wrong. Called for every solve with a
source, whether its operator is new or reused.
"""
function checkpumpconjugates(blocks::Vector{StampedScatteringBlock}, wmodes::AbstractVector)
    for sb in blocks
        block = sb.block
        block isa LinearizedScattering || continue
        n, nk, nm = block.nports, length(block.harmonics), length(wmodes)
        # the conjugate of a retained mode is the input, at the negative of
        # its frequency, so the entries are read as the stamp reads them
        # with the input frequencies negated
        Hin = Array{Complex{Float64},4}(undef, n, n, nk, nm)
        Hnegin = Array{Complex{Float64},4}(undef, n, n, nk, nm)
        evaluateharmonics!(Hin, block, -wmodes)
        evaluateharmonics!(Hnegin, block, wmodes)
        # an entry below the tolerance of the block's data is zero
        tol = block.atol*max(1.0, maximum(abs, Hin), maximum(abs, Hnegin))
        for m in 1:nm, nn in 1:nm
            (iszero(wmodes[m]) || iszero(wmodes[nn])) && continue
            d = wmodes[m] + wmodes[nn]
            k = round(Int, d/block.wp)
            (abs(d - k*block.wp) <= 1e-6*block.wp && abs(k) in block.harmonics) || continue
            j = findfirst(==(abs(k)), block.harmonics)
            worst = maximum(abs(conversionentry(Hin, Hnegin, k, j, p, q, nn)) for p in 1:n, q in 1:n)
            worst <= tol && continue
            throw(ArgumentError(lazy"the pumped block at $(sb.name) converts the conjugate of the mode at $(wmodes[nn]) rad/s into the mode at $(wmodes[m]) rad/s through its harmonic $(k), with an entry of $(worst) against the $(tol) its atol takes as zero, which the pump solve, complex linear in its retained modes, does not carry: solve the circuit with no source, its pump being the block's, give the pumped device as its junctions, or raise the block's atol if an entry of that size is to be taken as zero."))
        end
    end
    return nothing
end

"""
    evaluatehybridpumped!(B, C, block::LinearizedScattering, ws, K,
        work::ScatteringWorkspace)

The coefficient matrices of the hybrid constitutive equations of a pumped
block over the modes at the signed frequencies `ws`, as
[`evaluatehybrid!`](@ref) gives them for a block which does not convert,
with the modes coupled: `B[p, q, m, n]` and `C[p, q, m, n]` are the
entries of `(I - S) R^(-1/2)` and `(I + S) R^(1/2)` from input port `q`
at mode `n` to output port `p` at mode `m`, where the multi-mode
scattering matrix is `S[(p, m), (q, n)] = H_k[p, q](ws[n])` for the
harmonic `k = K[m, n]` (see [`pumpedharmonics`](@ref)), the negative
harmonics by `H_{-k}(nu) = conj(H_k(-nu))`, and zero where the modes are
not coupled. The rows and columns of a mode at zero frequency carry only
the identity, the `i = 0` row of the stamp.
"""
function evaluatehybridpumped!(B::AbstractArray{Complex{Float64},4},
    C::AbstractArray{Complex{Float64},4}, block::LinearizedScattering,
    ws::AbstractVector, K::AbstractMatrix{Int}, work::ScatteringWorkspace)

    n = block.nports
    nw = length(ws)
    nk = length(block.harmonics)
    if size(work.Hpos) != (n, n, nk, nw)
        work.Hpos = Array{Complex{Float64},4}(undef, n, n, nk, nw)
        work.Hneg = Array{Complex{Float64},4}(undef, n, n, nk, nw)
    end
    resize!(work.negws, nw)
    @inbounds for i in 1:nw
        work.negws[i] = -ws[i]
    end
    Hpos = evaluateharmonics!(work.Hpos, block, ws)
    Hneg = evaluateharmonics!(work.Hneg, block, work.negws)
    fill!(B, zero(Complex{Float64}))
    fill!(C, zero(Complex{Float64}))
    @inbounds for m in 1:nw
        if iszero(ws[m])
            for p in 1:n
                C[p, p, m, m] = one(Complex{Float64})
            end
            continue
        end
        for p in 1:n
            B[p, p, m, m] = 1/sqrt(block.zref[p])
            C[p, p, m, m] = sqrt(block.zref[p])
        end
        for nn in 1:nw
            iszero(ws[nn]) && continue
            k = K[m, nn]
            k == typemin(Int) && continue
            j = findfirst(==(abs(k)), block.harmonics)
            for q in 1:n
                rinv2 = 1/sqrt(block.zref[q])
                r2 = sqrt(block.zref[q])
                for p in 1:n
                    s = conversionentry(Hpos, Hneg, k, j, p, q, nn)
                    B[p, q, m, nn] -= rinv2*s
                    C[p, q, m, nn] += r2*s
                end
            end
        end
    end
    return B, C
end

"""
    assemblescattering!(A::SparseMatrixCSC, ssys::ScatteringStampSystem,
        wmodes::AbstractVector, work::ScatteringWorkspace = ScatteringWorkspace())

Add the frequency dependent constitutive entries of the scattering blocks
at the signed mode frequencies `wmodes` into the values of the system
matrix `A`: `sign*im*w_m*scale*B[p,q](w_m)` on the node flux columns and
`-C[p,q](w_m)` on the auxiliary current columns, with the coefficients
evaluated per block by [`evaluatehybrid!`](@ref) (which applies the
negative frequency rule of each block, so the negative mode entries carry
the complex conjugate data exactly as the conjugation of the conductance
matrix does for resistors). The destination indices must have been set with
[`setscatteringindexmap!`](@ref). Thread safe when each thread has its own
`work` and its own `A`: blocks are read only.
"""
function assemblescattering!(A::SparseMatrixCSC,
    ssys::ScatteringStampSystem, wmodes::AbstractVector,
    work::ScatteringWorkspace = ScatteringWorkspace())

    values = work.values
    length(values) == length(ssys.Aindex) ||
        resize!(values, length(ssys.Aindex))
    scatteringvalues!(values, ssys, wmodes, work)
    nzval = A.nzval
    @inbounds for c in eachindex(ssys.Aindex)
        nzval[ssys.Aindex[c]] += values[c]
    end
    return A
end

"""
    scatteringvalues!(values::AbstractVector, ssys::ScatteringStampSystem,
        wmodes::AbstractVector, work::ScatteringWorkspace)

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
    Bp = work.Bp
    Cp = work.Cp
    for (bi, sb) in enumerate(ssys.blocks)
        sb.block isa LinearizedScattering && continue
        evaluatehybrid!(Bs[bi], Cs[bi], sb.block, wmodes, work.hybrid)
    end
    for (j, bi) in enumerate(ssys.pumped)
        evaluatehybridpumped!(Bp[j], Cp[j], ssys.blocks[bi].block, wmodes,
            ssys.pumpedk[j], work)
    end

    @inbounds for c in eachindex(ssys.Aindex)
        m = ssys.modeindex[c]
        bi = ssys.blockindex[c]
        p = ssys.pindex[c]
        q = ssys.qindex[c]
        j = ssys.coupled[c]
        if j == 0
            values[c] = if ssys.coeff[c] == 1
                ssys.sign[c] * (im*wmodes[m]*ssys.scale) * Bs[bi][p, q, m]
            else
                -Cs[bi][p, q, m]
            end
        else
            # a pumped block: the column's own mode frequency scales the
            # node flux, whichever mode the row is
            nn = ssys.inmodeindex[c]
            values[c] = if ssys.coeff[c] == 1
                ssys.sign[c] * (im*wmodes[nn]*ssys.scale) * Bp[j][p, q, m, nn]
            else
                -Cp[j][p, q, m, nn]
            end
        end
    end
    return values
end

"""
    scatteringlinearterm(psc::CompiledCircuit, wmodes::AbstractVector,
        Nmodes::Integer; auxoffset::Integer, Ntotal::Integer,
        scale::Real = 1.0, blocks = nothing)

The constant sparse matrix of the scattering block contribution at the
fixed mode frequencies `wmodes` (the constitutive entries plus the
Kirchhoff current law couplings of the auxiliary port currents), or
`nothing` when the circuit has no scattering blocks. Used by the nonlinear
(pump) solver, where the mode frequencies do not change: the contribution
is folded into the frequency independent linear term alongside the
augmentation matrix of the promoted resistors, so the residual, Jacobian,
and solver machinery operate on the augmented system unchanged. `blocks`
stamps the given compiled blocks instead of `psc.scatteringblocks`.
"""
function scatteringlinearterm(psc::CompiledCircuit,
    wmodes::AbstractVector, Nmodes::Integer; auxoffset::Integer,
    Ntotal::Integer, scale::Real = 1.0, blocks = nothing)

    # the caller's blocks when it has them, and the circuit's own otherwise
    ssys = scatteringstampsystem(
        isnothing(blocks) ? psc.scatteringblocks : blocks,
        Nmodes; auxoffset = auxoffset, Ntotal = Ntotal, scale = scale,
        modeoffsets = wmodes)
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
ports therefore carries `n` noise channels, one per column of the
triangular factor of `I - S S'` ([`psdcholesky!`](@ref)); those of a
lossless block are identically zero, so only blocks
which are not [`provablylossless`](@ref) are given channels.

A block which states its noise with a [`NoiseCovariance`](@ref) `V`,
which is how an active block declares it, carries `2n` channels: `n`
from the factor of `(V + K)/2` with `K = I - S S'`, which emit like modes
in their vacuum, and `n` from the factor of `(V - K)/2`, which emit like
the conjugates of modes and enter the commutation relations with the
opposite sign (see [`noisechannelsigns`](@ref)). Their sum is `V`, their
difference `K`, so the block adds the noise it states and its output
obeys the commutation relations, and neither kind carries a temperature.

The channels of the blocks follow the noise ports of the dissipative
lumped components in the rows of `Snoise`, with the same
channel-major-mode-minor ordering.
"""
struct ScatteringNoisePlan
    # index into the blocks of the stamp system, the 0 based channel
    # offset of each, in channel (not row) units, and how many channels
    # each carries
    blockindices::Vector{Int}
    channelbase::Vector{Int}
    channelcounts::Vector{Int}
    Nchannels::Int
    Nmodes::Int
end

"""
    statednoise(block)

Whether a block states its noise with a [`NoiseCovariance`](@ref), and so
carries the channels of both kinds rather than those of its loss.
"""
statednoise(block::ScatteringParameters) = block.noise isa NoiseCovariance
statednoise(block::LinearizedScattering) = block.noise isa NoiseCovariance
statednoise(block) = false
"""
    planscatteringnoise(ssys)

The [`ScatteringNoisePlan`](@ref) of a [`ScatteringStampSystem`](@ref), or
`nothing` when the circuit has no scattering blocks or every block either
declares [`Lossless`](@ref) or is [`provablylossless`](@ref).

[`Passive`](@ref) and [`ThermalEquilibrium`](@ref) carry the covariance
of the block's loss, differing only in the temperature the channels are
at, and give a block one channel per port; [`Lossless`](@ref) declares
there are none; and a [`NoiseCovariance`](@ref) states the covariance
outright and gives a block two channels per port.
"""
planscatteringnoise(::Nothing) = nothing
function planscatteringnoise(ssys::ScatteringStampSystem)
    blockindices = Int[]
    channelbase = Int[]
    channelcounts = Int[]
    nch = 0
    for (bi, sb) in enumerate(ssys.blocks)
        checknoisemodel(sb.block, sb.name)
        # a block asserted lossless, or shown to be, adds nothing; a block
        # which states its noise adds it whatever its loss
        if !statednoise(sb.block)
            sb.block.noise isa Lossless && continue
            provablylossless(sb.block) && continue
        end
        push!(blockindices, bi)
        push!(channelbase, nch)
        count = statednoise(sb.block) ? 2*sb.block.nports : sb.block.nports
        push!(channelcounts, count)
        nch += count
    end
    isempty(blockindices) && return nothing
    return ScatteringNoisePlan(blockindices, channelbase, channelcounts, nch,
        ssys.Nmodes)
end

# The noise models the noise calculation supports. Silently ignoring an
# unsupported one would return a quantum efficiency computed for a
# different block than the user asked for, so it is an error here rather
# than a warning.
function checknoisemodel(block::ScatteringParameters, name)
    noise = block.noise
    if noise isa Passive || noise isa Lossless ||
            noise isa ThermalEquilibrium || noise isa NoiseCovariance
        return nothing
    end
    throw(ArgumentError(lazy"The scattering block at $(name) has the noise model $(noise), which the noise and quantum efficiency calculations do not support. Use Passive(), Lossless(), ThermalEquilibrium(T) or NoiseCovariance(V), or request no noise outputs."))
end
function checknoisemodel(block::LinearizedScattering, name)
    (block.noise isa Lossless || block.noise isa NoiseCovariance) && return nothing
    throw(ArgumentError(lazy"The pumped scattering block at $(name) has the noise model $(block.noise); a pumped block is lossless and emits nothing, or states the covariance its solve reports."))
end
checknoisemodel(block, name) = nothing

"""
    scatteringnoisenames(plan::ScatteringNoisePlan,
        ssys::ScatteringStampSystem)

The name of each noise channel of the plan, for labelling the rows of a
keyed noise scattering matrix. A block with more than one port has one
channel per port, distinguished by a channel number, because the columns
of the factor of `I - S S'` mix the ports and no single port owns a channel.
The channels of the conjugate kind of a block which states its noise
follow its channels of the first kind, marked with a prime.
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
        if plan.channelcounts[e] == 2n
            for c in 1:n
                names[plan.channelbase[e]+n+c] = n == 1 ? string(sb.name, "'") :
                    string(sb.name, "#", c, "'")
            end
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
    checklosslessblocks(ssys::Union{Nothing,ScatteringStampSystem}, w,
        wpumpmodes; atol = 1e-6, nsamples = 32)

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
    # a pumped block's declaration is checked over the modes of the solve
    # by checkpumpedblockmodels, its losslessness with its conversion
    unpumped(sb) = !(sb.block isa LinearizedScattering)
    any(sb -> unpumped(sb) && sb.block.noise isa Lossless &&
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
        (unpumped(sb) && sb.block.noise isa Lossless && !provablylossless(sb.block)) || continue
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
at the mode frequencies, the noise covariance and its factor in one flat
array, and the buffer of unsigned frequencies. Reused across frequencies
by one worker.
"""
mutable struct ScatteringNoiseWorkspace
    S::Array{Complex{Float64},3}
    # the covariance and its factor, flat and column major, so that the
    # factorization is the one the kernel runs
    L::Vector{Complex{Float64}}
    absws::Vector{Float64}
    # for a block which states its noise: the stated covariance at the
    # mode frequencies, and the covariance and factor of its channels of
    # the conjugate kind
    V::Array{Complex{Float64},3}
    M::Vector{Complex{Float64}}
end

function ScatteringNoiseWorkspace()
    return ScatteringNoiseWorkspace(
        Array{Complex{Float64},3}(undef, 0, 0, 0),
        Complex{Float64}[], Float64[],
        Array{Complex{Float64},3}(undef, 0, 0, 0), Complex{Float64}[])
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

A block which states its noise with a [`NoiseCovariance`](@ref) `V` has
its channels of the first kind weighted by the factor of `(V + K)/2` and,
in the rows after them, its channels of the conjugate kind by the factor
of `(V - K)/2`, with `K = I - S S'`, which [`checkblocknoisemodels`](@ref)
has admitted at every mode frequency of the sweep, or which a completed
covariance meets by construction.

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
        # goes behind a function barrier; a pumped block's channels span
        # its modes and read which harmonic couples each pair
        if sb.block isa LinearizedScattering
            j = findfirst(==(bi), ssys.pumped)
            pumpedblocknoisewaves!(noiseoutputwave, sb.block, sb, wmodes,
                phiadj, rowoffset + plan.channelbase[e]*Nmodes, Nmodes,
                ssys.pumpedk[j])
        else
            blocknoisewaves!(noiseoutputwave, sb.block, sb, wmodes, phiadj,
                rowoffset + plan.channelbase[e]*Nmodes, Nmodes, work)
        end
    end
    return noiseoutputwave
end

"""
    pumpednoisematrices(block::LinearizedScattering, wmodes, K)
    pumpednoisematrices(block::LinearizedScattering, rows, cols, K)

The multi-mode matrices of a pumped block at the signed mode frequencies
`wmodes`, over the index `(p - 1)*Nmodes + m` of port `p` at mode `m`,
in the units of the solver's outputs, waves of photons per second: its
scattering matrix `S`, from its harmonic transfer functions on power
waves through the harmonic `K[m, n]` coupling each pair of modes (see
[`pumpedharmonics`](@ref)) and the square root of the frequency ratio;
the commutator `J - S J S'` of the noise its output needs, with `J` the
signs of the mode frequencies; and, for a block which states its noise,
the covariance `V` from its harmonic covariances, `V_{-k}(nu)` being
`V_k(nu - k wp)'`, or `nothing`. The rows and columns of a mode at zero
frequency are left zero, the wave normalization being singular there.
With the output modes `rows` and the input modes `cols` given apart,
`K[m, n]` coupling `rows[m]` to `cols[n]`, `S` is the block from the
inputs to the outputs, the commutator `J_rows - S J_cols S'` is that of
the outputs with every input which reaches them, and `V` is over the
outputs: what validates a set of outputs without truncating the inputs
they are fed by (see [`pumpedfamily`](@ref)). The covariance of a
block whose noise model is completed is the completed one, from
[`completedcovariance`](@ref), unless `complete = false` asks for the
covariance as stated.
"""
pumpednoisematrices(block::LinearizedScattering, wmodes::AbstractVector, K::AbstractMatrix{Int}; complete::Bool = true) =
    pumpednoisematrices(block, wmodes, wmodes, K; complete = complete)
function pumpednoisematrices(block::LinearizedScattering, rows::AbstractVector, cols::AbstractVector, K::AbstractMatrix{Int};
        complete::Bool = true)
    n = block.nports
    nr, nc = length(rows), length(cols)
    size(K) == (nr, nc) || throw(DimensionMismatch(lazy"the harmonic map has size $(size(K)) for $(nr) output and $(nc) input modes."))
    nk = length(block.harmonics)
    Hp = Array{Complex{Float64},4}(undef, n, n, nk, nc)
    Hn = Array{Complex{Float64},4}(undef, n, n, nk, nc)
    evaluatecoveredharmonics!(Hp, block, cols)
    evaluatecoveredharmonics!(Hn, block, -cols)
    S = zeros(Complex{Float64}, n*nr, n*nc)
    idxr = (p, m) -> (p-1)*nr + m
    idxc = (q, nn) -> (q-1)*nc + nn
    for m in 1:nr, nn in 1:nc
        (iszero(rows[m]) || iszero(cols[nn])) && continue
        k = K[m, nn]
        k == typemin(Int) && continue
        j = findfirst(==(abs(k)), block.harmonics)
        photon = sqrt(abs(cols[nn])/abs(rows[m]))
        for q in 1:n, p in 1:n
            h = conversionentry(Hp, Hn, k, j, p, q, nn)
            S[idxr(p, m), idxc(q, nn)] = photon*h
        end
    end
    Jr = Diagonal([iszero(rows[m]) ? 0.0 : sign(rows[m]) for p in 1:n for m in 1:nr])
    Jc = Diagonal([iszero(cols[nn]) ? 0.0 : sign(cols[nn]) for q in 1:n for nn in 1:nc])
    Kc = Matrix(Jr) - S*Jc*S'
    V = statedcovariance(block, rows)
    complete && !isnothing(V) && block.noise.completed && (V = completedcovariance(block, rows, cols))
    return S, Kc, V
end

# the covariance a block states over the output modes `rows`, laid out
# over the modes and the ports as `pumpednoisematrices` lays out its
# multi-mode matrices, and `nothing` for a block which states none: the
# harmonic covariance of the difference of two modes where they are a
# harmonic apart, the negative harmonics from
# `V_{-k}(nu) = V_k(nu - k wp)'`, which makes the matrix Hermitian
# wherever each `V_0` is
function statedcovariance(block::LinearizedScattering, rows::AbstractVector)
    block.noise isa NoiseCovariance || return nothing
    n = block.nports
    nr = length(rows)
    nk = length(block.harmonics)
    Vp = Array{Complex{Float64},4}(undef, n, n, nk, nr)
    evaluatecoveredharmonics!(Vp, block, rows; covariance = true)
    V = zeros(Complex{Float64}, n*nr, n*nr)
    for m in 1:nr, nn in 1:nr
        (iszero(rows[m]) || iszero(rows[nn])) && continue
        d = rows[m] - rows[nn]
        k = round(Int, d/block.wp)
        (abs(d - k*block.wp) <= 1e-6*block.wp && abs(k) in block.harmonics) || continue
        j = findfirst(==(abs(k)), block.harmonics)
        for q in 1:n, p in 1:n
            V[(p-1)*nr + m, (q-1)*nr + nn] = k >= 0 ? Vp[p, q, j, nn] : conj(Vp[q, p, j, m])
        end
    end
    return V
end

"""
    completedcovariance(block::LinearizedScattering, rows, cols)

The covariance of a pumped block whose noise model is completed, over
the output modes `rows` fed by the input modes `cols` of a solve: the
covariance the block states, completed to the commutation relations
(see [`completecovariance`](@ref)) over the ladder of the rows padded
by the noise model's `padding` multiples of the pump frequency on
either side, with every input which feeds it, and restricted to the
rows, so that the block's noise is one model whatever modes a solve
keeps, to the precision of the padding; plus the vacuum of the inputs
which feed the rows and are not among `cols`, `S_c S_c'` for the
scattering `S_c` from them, which a solve without those modes traces
out, so that the covariance meets the commutator of the rows over the
inputs the solve has.

Every row a solve asks for is completed, including one the block's data
does not reach, which scatters nothing and so absorbs everything and
carries the vacuum its commutator requires, as the stamp of the block
takes it. The block couples nothing between frequencies which are not a
multiple of the pump apart, so the rows are completed one ladder at a
time (see `pumpladders`), the negative part of a matrix of independent
blocks being the negative part of each; a block which does not convert
couples nothing along a ladder either and needs no padding around a row.
"""
function completedcovariance(block::LinearizedScattering, rows::AbstractVector, cols::AbstractVector)
    n = block.nports
    nr = length(rows)
    out = zeros(Complex{Float64}, n*nr, n*nr)
    nr == 0 && return out
    for g in pumpladders(block.wp, rows)
        ng = length(g)
        Vg = completedladder(block, Float64[rows[i] for i in g], cols)
        for mm in 1:ng, m in 1:ng, q in 1:n, p in 1:n
            out[(p-1)*nr + g[m], (q-1)*nr + g[mm]] = Vg[(p-1)*ng + m, (q-1)*ng + mm]
        end
    end
    return out
end

# the indices of the frequencies `nus` grouped into the ladders of the
# pump, the sets whose members are a multiple of `wp` apart, which is
# the only way a pumped block couples two frequencies, in the order the
# groups first appear. Two frequencies are on one ladder when their
# residues modulo the pump agree, so they are grouped by sorted residue
# rather than compared with one another, which a bath of many ladders
# would make quadratic; the residues lie in a half open interval of one
# pump frequency whose two ends are the same ladder, which is why the
# outermost groups can still be one
function pumpladders(wp::Real, nus::AbstractVector)
    tol = 1e-6*wp
    res = Float64[nu - round(nu/wp)*wp for nu in nus]
    groups = Vector{Int}[]
    for i in sortperm(res)
        if !isempty(groups) && abs(res[i] - res[first(groups[end])]) <= tol
            push!(groups[end], i)
        else
            push!(groups, [i])
        end
    end
    if length(groups) > 1 && abs(abs(res[first(groups[end])] - res[first(groups[1])]) - wp) <= tol
        append!(groups[1], pop!(groups))
    end
    foreach(sort!, groups)
    return sort!(groups; by = first)
end

# the modes the completion of one ladder of rows reads: the rows padded
# by the noise model's multiples of the pump frequency on either side,
# with frequencies within roundoff of one another taken as one, and no
# padding for a block which does not convert, since it couples nothing
# along a ladder. This is the whole domain the negative part is taken
# over, so it is the domain a covariance is checked on
function paddedladder(block::LinearizedScattering, rows::AbstractVector)
    wp = block.wp
    pad = maximum(block.harmonics) == 0 ? 0 : block.noise.padding
    out = Float64[]
    for w in sort(Float64[r + j*wp for r in rows for j in -pad:pad])
        (isempty(out) || abs(w - out[end]) > 1e-9*(abs(w) + wp)) && push!(out, w)
    end
    return out
end

# the completed covariance of one ladder of rows, the body of
# completedcovariance over a set of frequencies the block couples
function completedladder(block::LinearizedScattering, rows::Vector{Float64}, cols::AbstractVector)
    n = block.nports
    wp = block.wp
    near(a, b) = abs(a - b) <= 1e-9*(abs(a) + wp)
    nr = length(rows)
    prows, pcols, pK = pumpedfamily(block, paddedladder(block, rows); reach = false, keep = true)
    S, Kc, V = pumpednoisematrices(block, prows, pcols, pK; complete = false)
    Vc = completecovariance(V, Kc)
    np, npc = length(prows), length(pcols)
    # every row is kept, so each is one of the padded rows
    index = Int[findfirst(r -> near(r, x), prows)::Int for x in rows]
    out = zeros(Complex{Float64}, n*nr, n*nr)
    # the inputs of the rows the solve lacks, whose vacuum is traced out
    absent = [nn for nn in 1:npc if !any(c -> near(c, pcols[nn]), cols)]
    for mm in 1:nr, m in 1:nr
        for q in 1:n, p in 1:n
            v = Vc[(p-1)*np + index[m], (q-1)*np + index[mm]]
            for nn in absent, pp in 1:n
                v += S[(p-1)*np + index[m], (pp-1)*npc + nn]*conj(S[(q-1)*np + index[mm], (pp-1)*npc + nn])
            end
            out[(p-1)*nr + m, (q-1)*nr + mm] = v
        end
    end
    return out
end

# the negative part of a Hermitian matrix, `-lambda v v'` summed over its
# negative eigenvalues, which is positive semidefinite
function negativepart(A::AbstractMatrix)
    e = eigen(Hermitian(Matrix{Complex{Float64}}(A)))
    return e.vectors*Diagonal([l < 0 ? -l : 0.0 for l in e.values])*e.vectors'
end

"""
    completecovariance(V, K)

The covariance `V` completed to the commutation relations of the
commutator `K`: `V + neg(V - K) + neg(V + K)` with `neg` the negative
part, which makes `V - K` and `V + K` positive semidefinite, so that
the channels formed from them exist (see [`NoiseCovariance`](@ref));
each addition is positive semidefinite and vanishes where `V` already
meets the relations. For `V = 0` it is `|K|`, the least total noise a
Gaussian channel with the map of `K` can add for its output to obey
the commutation relations, the least in the trace, which is the
`Ymin` of the quantum optics functions in the basis of the modes; for
a stated `V` the addition is sufficient and not in general the least.
The negative part of a matrix is not that of its parts, so a pumped
block is completed over a padded ladder (see
[`completedcovariance`](@ref)).
"""
completecovariance(V::AbstractMatrix, K::AbstractMatrix) = V + negativepart(V - K) + negativepart(V + K)

# a factor `F F' = A` of a Hermitian positive semidefinite matrix by its
# eigendecomposition: every positive eigenvalue is kept, however small,
# since a block with little loss emits little noise and not none, and a
# negative one is roundoff or a violation the check before the solve
# admitted, and is taken as zero
function psdfactor(A::AbstractMatrix)
    e = eigen(Hermitian(Matrix{Complex{Float64}}(A)))
    return e.vectors*Diagonal([l > 0 ? sqrt(l) : 0.0 for l in e.values])
end

"""
    checkpumpedblock(block::LinearizedScattering, wmodes, K, name, at)
    checkpumpedblock(block::LinearizedScattering, rows, cols, K, name, at)

Check the noise model of a pumped block over the signed mode frequencies
`wmodes` coupled by the harmonics `K`, or over the output modes `rows`
fed by the input modes `cols` (see [`pumpednoisematrices`](@ref)):
a declared [`Lossless`](@ref) requires `J - S J S'` to vanish there, to
the block's `atol`, and a stated [`NoiseCovariance`](@ref) requires
`V - K` and `V + K` to be positive semidefinite, to the block's `atol`
or the covariance's, whichever is larger, relative to the square of the
largest entry of `S`, and a covariance must have finite, Hermitian
entries, to the same tolerance of its own largest entry. Throws
otherwise, naming the block `name` and the frequency `at`.

The covariance a completed noise model states is held to its entries
like any other, since the completion repairs a covariance and does not
stand in for one: what it adds is the negative part of a Hermitian
matrix, which is not defined for data that is neither. It is held to
them over the whole domain the completion reads, the modes of each
ladder padded by the noise model's multiples of the pump frequency (see
`paddedladder`), and not over the modes of the solve alone, since that
is the matrix the negative part is taken of; the completion then reads
data this has passed. Only the positive semidefiniteness a completed
model establishes is left unchecked, so such a block reads its stated
covariance alone and never assembles the scattering matrix it would be
completed against.
"""
checkpumpedblock(block::LinearizedScattering, wmodes::AbstractVector, K::AbstractMatrix{Int}, name, at) =
    checkpumpedblock(block, wmodes, wmodes, K, name, at)
function checkpumpedblock(block::LinearizedScattering, rows::AbstractVector, cols::AbstractVector, K::AbstractMatrix{Int}, name, at)
    if block.noise isa NoiseCovariance && block.noise.completed
        for g in pumpladders(block.wp, rows)
            prows = paddedladder(block, Float64[rows[i] for i in g])
            checkcovarianceentries(block, statedcovariance(block, prows), name,
                lazy"$(at), over the ladder from $(first(prows)) to $(last(prows)) rad/s its completion is formed on")
        end
        return nothing
    end
    S, Kc, V = pumpednoisematrices(block, rows, cols, K)
    scale = max(1.0, maximum(abs, S))^2
    if isnothing(V)
        worst = maximum(abs, Kc)
        tol = block.atol*scale
        worst <= tol || throw(ArgumentError(lazy"the pumped block at $(name) declares noise = Lossless(), but at $(at) the largest entry of J - S J S' over the modes of the solve, with J the signs of the mode frequencies, is $(worst) against a tolerance of $(tol): a device with loss or gain states its noise with NoiseCovariance, and a fit states the noise its commutator requires."))
        return nothing
    end
    checkcovarianceentries(block, V, name, at)
    margin = min(minimum(real.(eigvals(Hermitian(V - Kc)))), minimum(real.(eigvals(Hermitian(V + Kc)))))
    tol = max(block.atol, block.noise.atol)*scale
    margin < -tol && throw(ArgumentError(lazy"the stated covariance of the pumped block at $(name) is less than the commutation relations require at $(at): the smallest eigenvalue of V - K or V + K, with K = J - S J S' and J the signs of the mode frequencies, is $(margin) against a tolerance of $(tol); see NoiseCovariance."))
    return nothing
end

# the entries a covariance must have wherever it is stated, finite and
# Hermitian, to the block's tolerance or the covariance's, whichever is
# larger, of its largest entry: what the completion of a covariance
# assumes of it, and what the channels formed from one assume
function checkcovarianceentries(block::LinearizedScattering, V::AbstractMatrix, name, at)
    all(isfinite, V) || throw(ArgumentError(lazy"the stated covariance of the pumped block at $(name) is not finite at $(at)."))
    skew = maximum(abs, V .- V')
    skew <= max(block.atol, block.noise.atol)*max(1.0, maximum(abs, V)) || throw(ArgumentError(lazy"the stated covariance of the pumped block at $(name) is not Hermitian at $(at): the largest entry of V - V' over those modes is $(skew); a covariance is Hermitian, V_{-k}(nu) = V_k(nu - k wp)'."))
    return nothing
end

# whether a table holds the frequency `nu`, to the roundoff its
# evaluation admits at the edges
function tablecovers(t::TabulatedMatrixProvider, nu::Real)
    f = t.frequencies
    edgetol = 8eps(Float64)*max(abs(f[1]), abs(f[end]))
    return f[1] - edgetol <= nu <= f[end] + edgetol
end

# whether a provider holds data at the frequency `nu`: a table within
# its knots, or everywhere when it declares how it extrapolates, a
# piecewise table within one of its bands, the data being the samples
# and what lies between them, with no declaration made beyond them,
# and a provider of any other kind, a callable, a constant or a
# filter, everywhere
providercovers(p::TabulatedMatrixProvider, nu::Real) = p.extrapolation != :error || holdsdata(p, nu)
providercovers(p::PiecewiseTabulatedProvider, nu::Real) = holdsdata(p, nu)
providercovers(p, nu::Real) = true

# whether a provider holds a sample of its own at `nu`, the knots of a
# table reaching it before any extrapolation, which is where its data can
# be checked against a relation it must obey; a provider which is not
# tabulated states its value everywhere
holdsdata(p::TabulatedMatrixProvider, nu::Real) = tablecovers(p, nu)
holdsdata(p::PiecewiseTabulatedProvider, nu::Real) = any(t -> tablecovers(t, nu), p.tables)
holdsdata(p, nu::Real) = true

# the knots of a tabulated provider, and none for one of any other kind
tableknots(p::TabulatedMatrixProvider) = p.frequencies
tableknots(p::PiecewiseTabulatedProvider) = piecewisefrequencies(p)
tableknots(p) = Float64[]

# a provider at the signed frequencies its data covers, into
# `dest[:, :, i]`, and zero where it does not
function evaluatecovered!(dest::AbstractArray{Complex{Float64},3}, p::AbstractMatrixProvider, ws::AbstractVector)
    n = size(dest, 1)
    fill!(dest, 0)
    inside = [providercovers(p, w) for w in ws]
    any(inside) || return dest
    buf = Array{Complex{Float64},3}(undef, n, n, count(inside))
    evaluateprovider!(buf, p, ws[inside])
    dest[:, :, inside] .= buf
    return dest
end

# the image of the signed frequency `nu` on the conjugate ladder of the
# harmonic `k`: the frequency whose harmonic covariance is the transpose
# of the one at `nu` (see conjugateladder)
ladderimage(k::Int, wp::Real, nu::Real) = -nu - k*wp

# a harmonic covariance of a block's stated noise at the signed
# frequencies `ws`, into `dest[:, :, i]`, zero where neither the
# frequency nor its image on the conjugate ladder is covered: the noise
# of a real wave obeys `V_k(-nu - k wp) = transpose(V_k(nu))` (see
# conjugateladder), so data of one sign states the noise at the other,
# which a covariance built from a solve needs, holding its rows at the
# modes of the solve and not at their conjugates
function evaluatecoveredcovariance!(dest::AbstractArray{Complex{Float64},3}, p::AbstractMatrixProvider,
        k::Int, wp::Real, ws::AbstractVector)
    n = size(dest, 1)
    fill!(dest, 0)
    inside = [providercovers(p, w) for w in ws]
    if any(inside)
        buf = Array{Complex{Float64},3}(undef, n, n, count(inside))
        evaluateprovider!(buf, p, [ws[i] for i in eachindex(ws) if inside[i]])
        dest[:, :, inside] .= buf
    end
    mirrored = [!inside[i] && providercovers(p, ladderimage(k, wp, ws[i])) for i in eachindex(ws)]
    if any(mirrored)
        buf = Array{Complex{Float64},3}(undef, n, n, count(mirrored))
        evaluateprovider!(buf, p, [ladderimage(k, wp, ws[i]) for i in eachindex(ws) if mirrored[i]])
        c = 0
        for i in eachindex(ws)
            mirrored[i] || continue
            c += 1
            dest[:, :, i] .= transpose(view(buf, :, :, c))
        end
    end
    return dest
end

# the harmonic transfer functions of a block, or with `covariance` the
# harmonic covariances of its stated noise, at the signed frequencies
# each harmonic's data covers, rotated by the pump phase as
# `evaluateharmonics!` and `evaluateharmoniccovariances!` rotate them,
# and zero where the data does not reach, which the harmonic map of a
# family never reads. A covariance is read on the conjugate ladder as
# well, the relation carrying it there being the same at either sign of
# the rotation
function evaluatecoveredharmonics!(dest::AbstractArray{Complex{Float64},4}, block::LinearizedScattering,
        ws::AbstractVector; covariance::Bool = false)
    n = block.nports
    buf = Array{Complex{Float64},3}(undef, n, n, length(ws))
    for (j, k) in enumerate(block.harmonics)
        if covariance
            evaluatecoveredcovariance!(buf, block.noise.provider[j], k, block.wp, ws)
        else
            evaluatecovered!(buf, block.providers[j], ws)
        end
        dest[:, :, j, :] .= cis(k*block.phase) .* buf
    end
    return dest
end

# whether the block holds the entry from the input frequency `nu`
# through the signed harmonic `k`, read as the stamp reads it: the
# harmonic `k` at `nu`, or the conjugate of the harmonic `-k` at `-nu`
function entrycovered(block::LinearizedScattering, k::Int, nu::Real)
    j = findfirst(==(abs(k)), block.harmonics)
    isnothing(j) && return false
    return providercovers(block.providers[j], k >= 0 ? nu : -nu)
end

# whether the block's declaration reaches the output mode at `nu`: a
# stated covariance must hold its own row there, `V_0(nu)`, for the
# commutation relations to be checked, at the frequency or at its image
# on the conjugate ladder, `-nu`, which states the same row transposed
# (see conjugateladder); a completed covariance reaches every output the
# harmonics do, being zero where it holds no data and completed to what
# the commutator requires there
function outputcovered(block::LinearizedScattering, nu::Real)
    block.noise isa NoiseCovariance || return true
    block.noise.completed && return true
    p = block.noise.provider[1]
    return providercovers(p, nu) || providercovers(p, ladderimage(0, block.wp, nu))
end

"""
    conjugateladder(block::LinearizedScattering)

Check the stated harmonic covariances of a pumped block against the
conjugate ladder and throw where they disagree. The noise of a real
wave has the covariance between the conjugates of two modes equal to
the conjugate of the covariance between the modes, so the harmonic
covariances obey

    V_k(-nu - k wp) = transpose(V_k(nu)),

which is how the noise at a frequency a covariance holds no data at is
read from the data at its image (see [`NoiseCovariance`](@ref)). Where
a table holds both a knot and its image the two state the same noise
and must agree, to the block's `atol` or the covariance's, whichever is
larger, of the largest entry of the table; a constant covariance is its
own image and so is symmetric, a Hermitian one real. A callable is
what it is at each frequency and is not checked.
"""
function conjugateladder(block::LinearizedScattering)
    block.noise isa NoiseCovariance || return nothing
    n = block.nports
    wp = block.wp
    atol = max(block.atol, block.noise.atol)
    for (j, k) in enumerate(block.harmonics)
        p = block.noise.provider[j]
        if p isa ConstantMatrixProvider
            M = Array{Complex{Float64},3}(undef, n, n, 1)
            evaluateprovider!(M, p, [0.0])
            d = maximum(abs, view(M, :, :, 1) .- transpose(view(M, :, :, 1)))
            d <= atol*max(1.0, maximum(abs, M)) || throw(ArgumentError(lazy"the constant covariance of the harmonic $(k) is not that of a real wave: it differs from its transpose by $(d), and a covariance which is the same at a frequency and at its image on the conjugate ladder is symmetric, V_k(-nu - k wp) = transpose(V_k(nu)). State a covariance which depends on frequency as a table or a callable."))
            continue
        end
        nus = tableknots(p)
        isempty(nus) && continue
        images = [ladderimage(k, wp, nu) for nu in nus]
        both = findall(im -> holdsdata(p, im), images)
        isempty(both) && continue
        A = Array{Complex{Float64},3}(undef, n, n, length(both))
        B = Array{Complex{Float64},3}(undef, n, n, length(both))
        evaluateprovider!(A, p, [nus[i] for i in both])
        evaluateprovider!(B, p, [images[i] for i in both])
        scale = max(1.0, maximum(abs, A))
        for c in eachindex(both)
            d = maximum(abs, view(B, :, :, c) .- transpose(view(A, :, :, c)))
            d <= atol*scale || throw(ArgumentError(lazy"the stated covariance of the harmonic $(k) is not that of a real wave: at $(nus[both[c]]) rad/s and at its image $(images[both[c]]) rad/s on the conjugate ladder the table differs from its transpose by $(d) against the largest entry $(scale), where the two state the same noise, V_k(-nu - k wp) = transpose(V_k(nu)); raise atol to admit a discrepancy of the data."))
        end
    end
    return nothing
end

"""
    pumpedfamily(block::LinearizedScattering, nus; reach = true, keep = false)

The output modes a pumped block's data reaches from the frequencies
`nus`, every input mode which feeds them, and the harmonic coupling
each pair, as `(rows, cols, K)` for [`pumpednoisematrices`](@ref): the
outputs are the signed frequencies a harmonic of the block apart from
one of `nus`, or `nus` themselves with `reach = false`, into which the
block holds an entry from an input a harmonic apart that its data
covers, and at which a stated covariance holds its row, and the inputs
are every such input, which reach one harmonic further out than the
outputs where the data goes so far. Frequencies within roundoff of one
another are one mode, and `K[m, n]` is `typemin(Int)` where the block
holds no entry from `cols[n]` to `rows[m]`. Validating the outputs
against every input which reaches them, rather than one set of modes
against itself, keeps a truncation of the inputs at the edge of the set
from counting as a violation by the data.

With `keep` every frequency is an output whether the block's data
reaches it or not, which is how a completed covariance is formed (see
[`completedcovariance`](@ref)): a frequency the data does not reach is
one the block scatters nothing at, an output of the circuit all the
same, and its row of the commutator is the one it has there.
"""
function pumpedfamily(block::LinearizedScattering, nus; reach::Bool = true, keep::Bool = false)
    kmax = maximum(block.harmonics)
    wp = block.wp
    near(a, b) = abs(a - b) <= 1e-9*(abs(a) + wp)
    function dedupe(ws)
        out = Float64[]
        for w in sort(ws)
            (isempty(out) || !near(w, out[end])) && push!(out, w)
        end
        return out
    end
    reached = reach ? dedupe(Float64[nu + j*wp for nu in nus for j in -kmax:kmax]) : dedupe(Float64.(collect(nus)))
    rows = Float64[]
    entries = Tuple{Int,Float64,Int}[]
    for r in reached
        (keep || outputcovered(block, r)) || continue
        found = false
        for kk in block.harmonics, k in (kk == 0 ? (0,) : (kk, -kk))
            c = r - k*wp
            entrycovered(block, k, c) || continue
            found || push!(rows, r)
            found = true
            push!(entries, (length(rows), c, k))
        end
        (found || !keep) || push!(rows, r)
    end
    cols = dedupe(Float64[c for (_, c, _) in entries])
    K = fill(typemin(Int), length(rows), length(cols))
    for (m, c, k) in entries
        i = searchsortedfirst(cols, c)
        nn = i <= length(cols) && near(cols[i], c) ? i : i - 1
        K[m, nn] = k
    end
    return rows, cols, K
end

# the violation of a pumped block's declaration over the signed output
# modes `rows` fed by the input modes `cols` through `K` (see
# pumpedfamily), or over one set of modes coupled by `K`, relative to
# the square of the largest entry of its multi-mode scattering matrix:
# the largest entry of `J - S J S'` for a declared lossless block, and
# how far below zero the smallest eigenvalue of `V - K` or `V + K`
# falls for one which states its noise; what checkpumpedblock holds to
# the tolerance
pumpedviolation(block::LinearizedScattering, wmodes::AbstractVector, K::AbstractMatrix{Int}) =
    pumpedviolation(block, wmodes, wmodes, K)
function pumpedviolation(block::LinearizedScattering, rows::AbstractVector, cols::AbstractVector, K::AbstractMatrix{Int})
    isempty(rows) && return 0.0
    S, Kc, V = pumpednoisematrices(block, rows, cols, K)
    scale = max(1.0, maximum(abs, S))^2
    isnothing(V) && return maximum(abs, Kc)/scale
    margin = min(minimum(real.(eigvals(Hermitian(V - Kc)))), minimum(real.(eigvals(Hermitian(V + Kc)))))
    return max(0.0, -margin, maximum(abs, V .- V'))/scale
end

"""
    pumpedblocknoisewaves!(noiseoutputwave, block::LinearizedScattering, sb,
        wmodes, phiadj, rowoffset, Nmodes, K)

The noise output waves of the channels of a pumped block which states
its noise, into the rows of `noiseoutputwave` after `rowoffset`, as
`blocknoisewaves!` does for a block which does not convert.

The block's noise is one wave over all its ports and modes at once, of
covariance `V` and commutator `K = J - S J S'` (see
[`pumpednoisematrices`](@ref)), so its channels are the columns of the
factors of `(V + K)/2`, which emit like modes in their vacuum, and of
`(V - K)/2`, which emit like their conjugates; each column spans the
port current rows of every mode. The block has `2 nports` channel slots
of `Nmodes` rows each in the noise scattering matrix, and the columns
are laid out over those rows, a column per row, the first `nports`
slots holding the first kind and the rest the second, so that every row
is one channel; their sign kinds are fixed (see
[`noisechannelsigns`](@ref)). Each row is the contraction of its column
against the adjoint solution over all the block's rows, with the sign
and the square root of each mode's frequency, and is multiplied by the
sign of its own row's mode frequency, which the adjoint route's sign
(see [`adjointnoisesigns!`](@ref)) then undoes, so that after it the
entry is the wave of the channel in the signed frequency convention.
"""
function pumpedblocknoisewaves!(noiseoutputwave::AbstractMatrix,
    block::LinearizedScattering, sb::StampedScatteringBlock,
    wmodes::AbstractVector, phiadj::AbstractMatrix, rowoffset::Integer,
    Nmodes::Integer, K::AbstractMatrix{Int})

    n = block.nports
    N = n*Nmodes
    nrhs = size(phiadj, 2)
    S, Kc, V = pumpednoisematrices(block, wmodes, K)
    # the factors by eigendecomposition: the covariance of a block with a
    # few internal channels is far from full rank, and the pivot rule of
    # the triangular factor cannot tell a pivot of roundoff from a small
    # one, where an eigenvalue can; the covariance was checked against
    # the commutator before the sweep (see checkblocknoisemodels)
    L = psdfactor((V .+ Kc) ./ 2)
    M = psdfactor((V .- Kc) ./ 2)
    @inbounds for kind in 1:2
        F = kind == 1 ? L : M
        for c in 1:N
            # the slot and the row of column c: the first nports slots
            # hold the first kind
            row = rowoffset + ((kind-1)*n)*Nmodes + (c-1) + 1
            mrow = (c-1) % Nmodes + 1
            srow = sign(wmodes[mrow])
            for k in 1:nrhs
                acc = zero(Complex{Float64})
                for p in 1:n, m in 1:Nmodes
                    iszero(wmodes[m]) && continue
                    acc += sign(wmodes[m])*sqrt(abs(wmodes[m]))*F[(p-1)*Nmodes + m, c]*
                        phiadj[sb.auxbase + (p-1)*Nmodes + m, k]
                end
                # the sign of the row's own mode and of the column's, which
                # the adjoint route's signs undo, so that the entry carries
                # the sign of each injection's mode alone
                noiseoutputwave[row, k] = srow*sign(wmodes[(k-1) % Nmodes + 1])*acc
            end
        end
    end
    return noiseoutputwave
end

# the noise output waves of one block, behind a function barrier because
# the block is stored untyped: the formula of `scatteringnoisewaves!`
function blocknoisewaves!(noiseoutputwave::AbstractMatrix,
    block::ScatteringParameters, sb::StampedScatteringBlock,
    wmodes::AbstractVector, phiadj::AbstractMatrix, rowoffset::Integer,
    Nmodes::Integer, work::ScatteringNoiseWorkspace)

    n = block.nports
    nrhs = size(phiadj, 2)
    stated = statednoise(block)
    nchannels = stated ? 2n : n
    if size(work.S) != (n, n, Nmodes)
        work.S = Array{Complex{Float64},3}(undef, n, n, Nmodes)
        resize!(work.L, n*n)
    end
    S = work.S
    L = work.L
    M = work.M
    evaluatescattering!(S, block, wmodes, work.absws)
    if stated
        if size(work.V) != (n, n, Nmodes)
            work.V = Array{Complex{Float64},3}(undef, n, n, Nmodes)
            resize!(work.M, n*n)
        end
        evaluatecovariance!(work.V, block, wmodes, work.absws)
        block.noise.completed && completestated!(work.V, S)
    end
    for m in 1:Nmodes
        if iszero(wmodes[m])
            for c in 1:nchannels
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
        if stated
            statednoisefactors!(L, M, work.V, m, n)
        else
            psdcholesky!(L, 0, n)
        end
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
            stated || continue
            for c in 1:n
                acc = zero(Complex{Float64})
                for p in 1:n
                    acc += M[(c-1)*n + p]*
                        phiadj[sb.auxbase + (p-1)*Nmodes + m, k]
                end
                noiseoutputwave[rowoffset + (n+c-1)*Nmodes + m, k] = kw*acc
            end
        end
    end
    return noiseoutputwave
end

# the stated covariance of an ordinary block at every mode completed to
# the commutation relations of its scattering matrix there (see
# completecovariance), for a block whose covariance is completed
function completestated!(V::AbstractArray{<:Complex,3}, S::AbstractArray{<:Complex,3})
    n = size(S, 1)
    for m in axes(S, 3)
        Sm = view(S, :, :, m)
        K = Matrix{Complex{Float64}}(I, n, n) - Sm*Sm'
        V[:, :, m] .= completecovariance(view(V, :, :, m), K)
    end
    return V
end

"""
    statednoisefactors!(L, M, V, m, n)

The factors of the channels of a block which states its noise: on entry
the lower triangle of `L` holds `K = I - S S'` at mode `m` and `V[:,:,m]`
the stated covariance there; on exit `L` holds the triangular factor of
`(V + K)/2` and `M` that of `(V - K)/2`, both by [`psdcholesky!`](@ref),
which are positive semidefinite for a covariance
[`checkblocknoisemodels`](@ref) has admitted or which is completed.
"""
function statednoisefactors!(L, M, V, m::Integer, n::Integer)
    @inbounds for c in 1:n
        for p in c:n
            k = L[(c-1)*n + p]
            v = V[p, c, m]
            L[(c-1)*n + p] = (v + k)/2
            M[(c-1)*n + p] = (v - k)/2
        end
    end
    psdcholesky!(L, 0, n)
    psdcholesky!(M, 0, n)
    return nothing
end

"""
    checkblocknoisemodels(ssys::Union{Nothing,ScatteringStampSystem}, w,
        wpumpmodes)

Check the stated noise of the blocks of `ssys` at every frequency the
sweep will solve at, every signal frequency at every pump mode, and
throw if one is not met there: a block which states its noise with a
[`NoiseCovariance`](@ref) must state at least what the commutation
relations require (see [`quantumnoisemargin`](@ref)), to the `atol` of
the model. A pumped block is checked by
[`checkpumpedblockmodels`](@ref).

Stored data is checked when a block is constructed, but a callable
cannot be, and a fitted block differs from its data. These are exactly
the frequencies the channels are formed at, so nothing the sweep will
use goes unchecked, and the check is done here rather than in the
sweep so that a refusal reaches the caller before any solve, as a
plain error rather than the failure of a task.
"""
function checkblocknoisemodels(ssys::Union{Nothing,ScatteringStampSystem}, w,
    wpumpmodes)

    isnothing(ssys) && return nothing
    any(sb -> statednoise(sb.block) && !(sb.block isa LinearizedScattering), ssys.blocks) || return nothing
    ws = Float64[w[i] + m for i in eachindex(w) for m in wpumpmodes]
    filter!(!iszero, ws)
    isempty(ws) && return nothing
    work = HybridWorkspace()
    for sb in ssys.blocks
        (statednoise(sb.block) && !(sb.block isa LinearizedScattering)) || continue
        checkoneblockstated(sb.block, sb.name, ws, work)
    end
    return nothing
end

"""
    checkpumpedblockmodels(ssys::Union{Nothing,ScatteringStampSystem}, w,
        wpumpmodes)

Check what every pumped block of `ssys` declares, its losslessness or
its stated covariance, over the modes of the solve at each signal
frequency (see [`checkpumpedblock`](@ref)), and throw where it is not
met. The sweep's modes are the one place the multi-mode matrices of a
pumped block are known, and the family its noise is completed over.
Runs at every solve, whether or not a noise output is asked for, since
a lossless pumped block emits nothing only if it is lossless: a block
is accepted or refused on its declaration and not on the outputs. A
block whose covariance is completed has what it states checked for
finite, Hermitian entries and the relations it meets by construction
left alone.
"""
function checkpumpedblockmodels(ssys::Union{Nothing,ScatteringStampSystem}, w,
    wpumpmodes)

    isnothing(ssys) && return nothing
    wmodes = zeros(Float64, length(wpumpmodes))
    for (bi, sb) in enumerate(ssys.blocks)
        sb.block isa LinearizedScattering || continue
        j = findfirst(==(bi), ssys.pumped)
        for i in eachindex(w)
            wmodes .= w[i] .+ wpumpmodes
            checkpumpedblock(sb.block, wmodes, ssys.pumpedk[j], sb.name, lazy"the signal frequency $(w[i]) rad/s")
        end
    end
    return nothing
end

# behind a function barrier: the block is stored untyped
function checkoneblockstated(block::ScatteringParameters, name, ws, work)
    n = block.nports
    S = Array{Complex{Float64},3}(undef, n, n, length(ws))
    V = Array{Complex{Float64},3}(undef, n, n, length(ws))
    evaluatescattering!(S, block, ws, work.absws)
    evaluatecovariance!(V, block, ws, work.absws)
    atol = block.noise.atol
    for k in eachindex(ws)
        Vk = view(V, :, :, k)
        skew = maximum(abs, Vk .- Vk')
        skew <= atol*max(1.0, maximum(abs, Vk)) || throw(ArgumentError(lazy"the noise covariance of the scattering block at $(name) is not Hermitian at $(ws[k]) rad/s: the largest entry of V - V' is $(skew)."))
        block.noise.completed && continue
        margin = quantumnoisemargin(Vk, view(S, :, :, k))
        if margin < -atol
            throw(ArgumentError(lazy"the noise covariance of the scattering block at $(name) is less than the commutation relations require at $(ws[k]) rad/s: the smallest eigenvalue of V - K or V + K, with K = I - S S', is $(margin). An amplifier of power gain G has to emit at least G - 1 at its output; see NoiseCovariance."))
        end
    end
    return nothing
end

"""
    noisechanneltemperatures(psc::CompiledCircuit,
        noiseportimpedanceindices, noiseplan, ssys, temperature)

The temperature of each row of the noise scattering matrix, in the order
[`noisechannelnames`](@ref) gives them.

`temperature` is the analysis default, which every dissipative element takes
unless it states one of its own. A lumped component states it as
`Resistor(R; temperature = T)` and a [`ScatteringParameters`](@ref) as
`noise = ThermalEquilibrium(T)`, both of which are recorded by
[`compile`](@ref) as it lowers the circuit. Only the typed circuit
format carries them; a netlist of tuples states none and everything in it
takes the default.

The channels of a block which states its noise with a
[`NoiseCovariance`](@ref) are at zero temperature: the covariance it
states is the whole of its noise, and the occupation of one is what
leaves it so.
"""
function noisechanneltemperatures(psc, noiseportimpedanceindices, noiseplan,
    ssys, temperature)

    stated = psc.componenttemperatures
    ts = Float64[get(stated, i, Float64(temperature))
        for i in noiseportimpedanceindices]
    isnothing(noiseplan) && return ts
    # a block states its temperature on its own noise model
    for (e, bi) in enumerate(noiseplan.blockindices)
        sb = ssys.blocks[bi]
        t = sb.block.noise isa ThermalEquilibrium ?
            Float64(sb.block.noise.temperature) :
            statednoise(sb.block) ? 0.0 : Float64(temperature)
        for _ in 1:noiseplan.channelcounts[e]
            push!(ts, t)
        end
    end
    return ts
end

"""
    noisechannelsigns(noiseportimpedanceindices, noiseplan, ssys)

The sign kind of each channel of the noise scattering matrix in the
commutation relations, in the order [`noisechannelnames`](@ref) gives
them, or `nothing` when every channel is of the first kind: `1` for a
channel which emits like a mode in its vacuum, whose rows count with
the sign of their mode frequency; `-1` for a channel of the conjugate
kind, the second half of the channels of a block which states its noise
(see [`ScatteringNoisePlan`](@ref)), whose rows count with the opposite
sign, since what it emits is the conjugate of a mode; and `2` and `-2`
for the channels of a pumped block which states its noise, which span
all its modes at once, so that every row of one counts with the fixed
sign of its kind (see [`pumpedblocknoisewaves!`](@ref)).
"""
function noisechannelsigns(noiseportimpedanceindices, noiseplan, ssys)
    isnothing(noiseplan) && return nothing
    any(e -> statednoise(ssys.blocks[noiseplan.blockindices[e]].block),
        eachindex(noiseplan.blockindices)) || return nothing
    signs = ones(Float64, length(noiseportimpedanceindices))
    for (e, bi) in enumerate(noiseplan.blockindices)
        block = ssys.blocks[bi].block
        n = block.nports
        count = noiseplan.channelcounts[e]
        fixed = block isa LinearizedScattering ? 2.0 : 1.0
        for c in 1:count
            push!(signs, c > n ? -fixed : fixed)
        end
    end
    return signs
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
