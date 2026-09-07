# The physical baths of a transient circuit: one independent equilibrium
# bath per matched port termination and per finite internal resistor, the
# amplitudes of the Norton currents that represent one quadrature pair of
# a bath at one frequency, and the accumulation of covariance and
# commutator from responses. The port owned resistor is the external
# environment, not a second internal source. The noise itself is
# [`transientnoise`](@ref), the tangent or the adjoint of the recorded
# steps driven at the baths.

"""
    TransientNoiseBath

One bath: its name, the rows of the equations its source enters with
its weights there, the port it terminates (zero otherwise), its
resistance (zero for a block channel) and its temperature. A resistor's
source is a Norton current into its first terminal and out of its
second, which for a port's termination is the port's own source
convention, so the current enters the port wave as a port source's
does; a block channel's source is an emitted noise wave, a combination
of the block's ports weighted by an eigenvector of `I - S S'` times the
square root of its eigenvalue, entering the port current rows as the
source `2 eta` of the hybrid equation. [`transientinjection`](@ref)
builds the columns.
"""
struct TransientNoiseBath
    name::String
    # the rows of the equations the bath's source enters and its weights
    # there: a resistor's Norton current into its two nodes, a block
    # channel's emitted wave into the block's port current rows
    rows::Vector{Int}
    weights::Vector{Float64}
    port::Int
    # the resistance of a resistor bath, or zero for a wave channel
    resistance::Float64
    temperature::Float64
end

"""
    TransientNoiseBaths

The baths of a [`TransientProblem`](@ref), from [`transientnoisebaths`](@ref).
"""
struct TransientNoiseBaths
    problem::TransientProblem
    channels::Vector{TransientNoiseBath}
    # the groups of channels whose covariance depends on frequency: the
    # ports of a rational block as channels `first:last` of the list, and
    # the block, whose scattering matrix at a frequency gives the
    # covariance `I - S S'` of its emitted wave
    groups::Vector{@NamedTuple{channels::UnitRange{Int}, block::Int}}
end
TransientNoiseBaths(p::TransientProblem, channels::Vector{TransientNoiseBath}) =
    TransientNoiseBaths(p, channels, @NamedTuple{channels::UnitRange{Int}, block::Int}[])
Base.length(b::TransientNoiseBaths) = length(b.channels)

# the covariance of a group's emitted wave at a frequency in Hz
function groupcovariance(baths::TransientNoiseBaths, group, frequency)
    b = baths.problem.blocks[group.block]
    S = zeros(ComplexF64, length(b.signal), length(b.signal), 1)
    evaluateprovider!(S, b.definition.provider, [2pi*frequency])
    K = Hermitian(I - S[:, :, 1]*S[:, :, 1]')
    # a block active at the frequency has no equilibrium noise; its
    # construction should have refused it, and this is the last check
    worst = minimum(eigvals(K))
    worst >= -1e-6 || throw(ArgumentError(
        lazy"the scattering block at $(b.path) is active at $(frequency) Hz: the minimum eigenvalue of I - S S' is $(worst), so it has no equilibrium noise."))
    return K
end

"""
    transientnoisebaths(problem; temperature = 0.0, porttemperatures = nothing)

The independent equilibrium baths of a circuit in time: every matched,
port owned termination is one external bath, every finite internal
resistor one internal bath, and every lossy scattering block the
independent channels of its emitted noise wave, whose covariance is
`I - S S'` (Bosma's relation, as the linearized solver has it), one
channel per positive eigenvalue; from the compiler's termination
ownership, the bound values, the component temperatures and the
blocks' noise models, `ThermalEquilibrium(T)` stating a block's
temperature, `Passive()` taking the default, `Lossless()` asserting the
block emits nothing, which is checked. `temperature` is the default in
kelvin and `porttemperatures` overrides the external baths in compiled
port order. Every port must own a matched finite termination; an open
resistor adds no bath. The same temperatures and models set the noise
of [`hblinsolve`](@ref), so the two solvers compare.
"""
function transientnoisebaths(p::TransientProblem; temperature = 0.0, porttemperatures = nothing)
    c, vvn = p.circuit, p.matrices.vvn
    t = Float64(temperature)
    isfinite(t) && t >= 0 || throw(ArgumentError("the bath temperature must be finite and nonnegative."))
    isnothing(porttemperatures) || length(porttemperatures) == length(c.ports) || throw(DimensionMismatch(
        "give one external temperature per port."))
    channels = TransientNoiseBath[]
    groups = @NamedTuple{channels::UnitRange{Int}, block::Int}[]
    for (j, port) in enumerate(c.ports)
        isapprox(p.portconductances[j]*p.portimpedances[j], 1; rtol = 1e-12) || throw(ArgumentError(
            lazy"port $(port.number) needs its own matched termination for a bath."))
        temp = isnothing(porttemperatures) ? get(c.componenttemperatures, port.environment, t) :
            Float64(porttemperatures[j])
        push!(channels, TransientNoiseBath("port $(port.number)", noderows(p.portpositive[j], p.portnegative[j])...,
            j, p.portimpedances[j], temp))
    end
    for k in noiseindices(c, vvn)
        c.componenttypes[k] == :R || throw(ArgumentError("the transient baths are real resistors and scattering blocks."))
        r = vvn[k]
        isfinite(r) || continue
        push!(channels, TransientNoiseBath(c.componentnames[k], noderows(c.nodeindices[1, k] - 1, c.nodeindices[2, k] - 1)...,
            0, r, get(c.componenttemperatures, k, t)))
    end
    # A block emits the noise its loss requires, the wave with covariance
    # I - S S' per Bosma, as the linearized solver computes it: the
    # independent channels are the eigenvectors of that covariance scaled
    # by the square roots of its eigenvalues, each entering the block's
    # port current rows as the source 2 eta of the hybrid equation. A
    # lossless block emits nothing, and says so or is checked.
    for b in p.blocks
        noise = b.definition.noise
        noise isa NoiseCovariance && throw(ArgumentError(
            lazy"the scattering block at $(b.path) has an arbitrary noise covariance, which the transient does not support; use Passive, Lossless or ThermalEquilibrium."))
        K = Symmetric(I - b.S*transpose(b.S))
        if noise isa Lossless
            # a declared lossless rational block is validated by its norms
            lossless = size(b.A, 1) > 0 ? losslessnorms(b.definition.provider) : maximum(abs, K) <= 1e-10
            lossless || throw(ArgumentError(
                lazy"the scattering block at $(b.path) declares noise = Lossless(), but I - S S' does not vanish at every frequency; a block which dissipates must carry the noise its loss requires."))
            continue
        end
        temp = noise isa ThermalEquilibrium ? Float64(noise.temperature) : t
        n = length(b.signal)
        if size(b.A, 1) > 0
            # a rational block's covariance depends on frequency: one
            # channel per port, correlated by the group's covariance in
            # the contraction
            first = length(channels) + 1
            for q in 1:n
                push!(channels, TransientNoiseBath(string(b.path, " port ", q), [b.auxbase + q], [2.0], 0, 0.0, temp))
            end
            push!(groups, (channels = first:first + n - 1, block = findfirst(x -> x === b, p.blocks)))
            continue
        end
        values, vectors = eigen(K)
        for c in 1:n
            values[c] > 1e-12 || continue
            push!(channels, TransientNoiseBath(string(b.path, " channel ", c), collect(b.auxbase + 1:b.auxbase + n),
                2sqrt(values[c]) .* vectors[:, c], 0, 0.0, temp))
        end
    end
    all(b -> isfinite(b.temperature) && b.temperature >= 0, channels) || throw(ArgumentError(
        "the bath temperatures must be finite and nonnegative."))
    isempty(channels) && throw(ArgumentError("the circuit has no dissipative element to be a bath."))
    return TransientNoiseBaths(p, channels, groups)
end

# The contraction of the covariance and the commutator for the groups
# whose channels are correlated: the group's emitted wave has the
# covariance `K(f) = I - S S'`, so `H conj(K) H'` is added, with
# `H = R_c - i R_s` the complex response of the group's cosine and sine
# quadratures in the demodulation's phasor convention, the real part
# weighted by the occupation into the covariance and the imaginary part
# into the commutator, as the independent channels are; the group's
# columns are masked out of the independent accumulation, so a nearly
# lossless block's small covariance is never the difference of two
# large ones. That convention is the conjugate of the linearized
# solver's, whose `S` gives `K`, so the wave's covariance reads
# `conj(K)` here; the conjugate only matters when `K` has imaginary
# off-diagonal entries, a scalar or real `K` hides it, and the block
# with a feedthrough of opposite signs at its ports in the tests does
# not. The loss matrix of every group at every frequency is evaluated
# once, since it depends on no condition, and the group's columns of a
# condition's response are gathered and contracted where the response
# lives, on the backend, so nothing of the response is copied to the
# host.
struct GroupCorrection{C, M, H, I, W}
    columns::Vector{Vector{I}}
    losses::Vector{Vector{C}}
    occupations::Vector{Vector{Float64}}
    # the mask of the independent columns, zero on the groups'
    mask::W
    R::M
    H::H
    T::H
    A::H
end

function groupcorrection(baths::TransientNoiseBaths, frequencies, m::Int, backend)
    nb = length(baths)
    columns, losses, occupations = Vector{Vector{Int}}[], Vector{Matrix{ComplexF64}}[], Vector{Float64}[]
    pmax = 0
    mask = ones(2nb*length(frequencies))
    for group in baths.groups
        cols = [[2*((f - 1)*nb + b - 1) + q for b in group.channels for q in 1:2] for f in eachindex(frequencies)]
        foreach(c -> (mask[c] .= 0), cols)
        Ks = [Matrix(conj(groupcovariance(baths, group, frequency))) for frequency in frequencies]
        occ = [thermaloccupation(2pi*frequency, baths.channels[first(group.channels)].temperature)/2 for frequency in frequencies]
        push!(columns, cols); push!(losses, Ks); push!(occupations, occ)
        pmax = max(pmax, length(group.channels))
    end
    allocate = (T, dims...) -> KernelAbstractions.zeros(backend, T, dims...)
    return GroupCorrection([[tobackend(backend, c) for c in cols] for cols in columns], [[tobackend(backend, K) for K in Ks] for Ks in losses],
        occupations, tobackend(backend, mask), allocate(Float64, m, 2pmax), allocate(ComplexF64, m, pmax), allocate(ComplexF64, m, pmax),
        allocate(ComplexF64, m, m))
end

# the groups' columns of a response masked out, for the independent
# accumulation, once the groups have been contracted
maskgroups!(response, gc::GroupCorrection) = (response .*= transpose(gc.mask); response)

function groupcorrection!(covariance, commutator, response, gc::GroupCorrection)
    for g in eachindex(gc.columns), f in eachindex(gc.columns[g])
        K = gc.losses[g][f]
        p = size(K, 1)
        R = view(gc.R, :, 1:2p)
        R .= view(response, :, gc.columns[g][f])
        H = view(gc.H, :, 1:p)
        H .= view(R, :, 1:2:2p) .- im .* view(R, :, 2:2:2p)
        T = view(gc.T, :, 1:p)
        mul!(T, H, K)
        mul!(gc.A, T, H')
        covariance .+= gc.occupations[g][f] .* real.(gc.A)
        commutator .+= imag.(gc.A)
    end
    return nothing
end

# the baths as targets of the tangent and the adjoint: the injection of a
# unit current into the first terminal and out of the second, and the port
# each bath terminates
function transientinjection(p::TransientProblem, baths::TransientNoiseBaths)
    # the baths are the circuit's, shared by every rebinding of its sources
    (baths.problem.circuit === p.circuit && baths.problem.matrices === p.matrices) || throw(ArgumentError(
        "the baths belong to another circuit."))
    rows, cols, vals = Int[], Int[], Float64[]
    for (k, b) in enumerate(baths.channels), (r, w) in zip(b.rows, b.weights)
        push!(rows, r); push!(cols, k); push!(vals, w)
    end
    return sparse(rows, cols, vals, length(p), length(baths))
end
targetports(::TransientProblem, baths::TransientNoiseBaths) = [b.port for b in baths.channels]

"""
    bathamplitude(bath, frequency, weight)

The peak amplitude `2 sqrt(h f df/R)` of the cosine and sine Norton currents
representing one quadrature pair of `bath` at the frequency `f` in Hz
with the quadrature weight `df` in Hz. With the independent quadratures at
variance `thermaloccupation(2pi f, T)/2 = nbar + 1/2` this gives the
bilateral symmetrized current spectral density `h f/R coth(h f/2kT)`,
`2kT/R` classically.
"""
bathamplitude(b::TransientNoiseBath, frequency, weight) = iszero(b.resistance) ?
    sqrt(planck_constant*frequency*weight) : 2sqrt(planck_constant*frequency*weight/b.resistance)

# the rows and weights of a Norton current out of the first node and into
# the second, skipping ground
function noderows(n1::Int, n2::Int)
    rows, weights = Int[], Float64[]
    n1 > 0 && (push!(rows, n1); push!(weights, 1.0))
    n2 > 0 && (push!(rows, n2); push!(weights, -1.0))
    return rows, weights
end

"""
    transientnoiseaccumulate!(covariance, commutator, response, variances)

Add the contribution of the responses `response[:, 2j-1:2j]` of the
quadrature pairs to the covariance, weighted by the pairs' `variances`,
and to the commutator, which the occupation does not weight.
"""
function transientnoiseaccumulate!(covariance, commutator, response, variances)
    size(response, 2) == length(variances) && iseven(length(variances)) || throw(DimensionMismatch(
        "the responses come in quadrature pairs."))
    x = view(response, :, 1:2:size(response, 2))
    p = view(response, :, 2:2:size(response, 2))
    mul!(commutator, x, transpose(p), 1.0, 1.0)
    mul!(commutator, p, transpose(x), -1.0, 1.0)
    weighted = response .* transpose(sqrt.(variances))
    mul!(covariance, weighted, transpose(weighted), 1.0, 1.0)
    return nothing
end

# the conditions grouped by identical initial state, since a stationary
# operator depends on nothing else: the groups as lists of condition
# indices
function initialstategroups(x0s)
    groups = Vector{Int}[]
    keys = Vector{Float64}[]
    for j in axes(x0s, 2)
        x = Array(view(x0s, :, j))
        g = findfirst(==(x), keys)
        isnothing(g) ? (push!(keys, x); push!(groups, [j])) : push!(groups[g], j)
    end
    return groups
end

# The stationary response of the circuit at the initial state to a
# sinusoidal drive of the baths, `[-w^2 C + i w G + L + J'(x_0)] x = b`, in
# the scaled system of the step: the fluctuations already stored in the
# capacitors and inductors before the drive, and their correlation with
# the drive to come, which starting the responses from zero would drop.
# The real and imaginary parts of the complex response are the initial
# flux and rate of the cosine and sine quadratures.
# The stationary operator of the circuit linearized at `x0` at the
# angular frequency `w`, on the flux phasors and, with lines, the phasors
# of the waves leaving the line ports: the step's own equations in
# frequency, `-w^2 C + i w G + L + J'(x0)` on the nodes with the lines'
# forced currents `2 q / sqrt(Z)` from the far ports' waves a delay
# earlier, `q = P a` with `P` the swap times `exp(-i w tau)`, and the
# waves' own equations `a = i w phi0 (E' x) / sqrt(Z) - P a`. Nothing is
# inverted in the lines, so the operator is regular at a line's half wave
# resonances, where its admittance is not.
function stationaryoperator(sys::TransientSystem, x0, w)
    p = sys.problem
    n = length(p)
    C, G, L = hostsparse(sys.C), hostsparse(sys.G), hostsparse(sys.L)
    RJ = hostsparse(sys.RJ)
    phi = RJ*Array(x0)
    F = -w^2 .* C .+ (im*w) .* G .+ L .+ transpose(RJ)*Diagonal(Array(sys.lmolj) .* cos.(phi))*RJ
    blockstates(p) == 0 || (F = F .+ rationalmatrix(p, im*w, sys.Lscale, n))
    isempty(p.lines) && return F
    nl2 = 2length(p.lines)
    E = hostsparse(sys.linescatter)
    Linj = hostsparse(sys.lineinjection)
    z = [line.Z for line in p.lines for _ in 1:2]
    P = sparse([2l - 1 for l in eachindex(p.lines)], [2l for l in eachindex(p.lines)], [cis(-w*line.delay) for line in p.lines], nl2, nl2)
    P = P + transpose(P)
    return [F  (-Linj*Diagonal(2 ./ sqrt.(z))*P); (-(im*w*phi0) .* Diagonal(1 ./ sqrt.(z))*sparse(transpose(E)))  (I + P)]
end

# The stationary response of the circuit at the initial state to a
# sinusoidal drive of the baths, in the scaled system of the step: the
# fluctuations already stored in the capacitors, inductors and lines
# before the drive, and their correlation with the drive to come, which
# starting the responses from zero would drop. The real and imaginary
# parts of the complex response are the initial flux and rate of the
# cosine and sine quadratures, and the samples of the waves leaving the
# line ports over the prehistory the responses read.
# the incident waves of the rational blocks as a linear map on the
# stationary solution, `(nports, n + nl2)` complex: the port voltage
# phasor `i w phi0 (x_signal - x_ref)` and the current `phi0 u / Lscale`
function incidentmap(sys::TransientSystem, w, ncols)
    p = sys.problem
    rows, cols, vals = Int[], Int[], ComplexF64[]
    prow = 0
    for b in p.blocks
        for q in eachindex(b.signal)
            b.signal[q] > 0 && (push!(rows, prow + q); push!(cols, b.signal[q]); push!(vals, im*w*phi0/(2sqrt(b.R[q]))))
            b.ref[q] > 0 && (push!(rows, prow + q); push!(cols, b.ref[q]); push!(vals, -im*w*phi0/(2sqrt(b.R[q]))))
            push!(rows, prow + q); push!(cols, b.auxbase + q); push!(vals, sqrt(b.R[q])*phi0/(2sys.Lscale))
        end
        prow += length(b.signal)
    end
    return sparse(rows, cols, vals, prow, ncols)
end

# the state phasors of the rational blocks from the incident wave
# phasors, `(i w I - A)^(-1) B a` per block
function statephasors(p::TransientProblem, w, a)
    z = zeros(ComplexF64, blockstates(p), size(a, 2))
    prow = 0
    for b in p.blocks
        np, nz = length(b.signal), size(b.A, 1)
        nz > 0 && (z[b.zbase + 1:b.zbase + nz, :] .= (im*w*I - b.A) \ (b.B*a[prow + 1:prow + np, :]))
        prow += np
    end
    return z
end

function stationaryresponses(sys::TransientSystem, x0, injection, frequencies, t0, reference)
    p = sys.problem
    n = length(p)
    nl2 = 2length(p.lines)
    nzs = blockstates(p)
    npre = lineprehistory(p, sys.h)
    # the injection in the scaled system, as the tangent scales it
    injection = (sys.Lscale/phi0) .* injection
    nb = size(injection, 2)
    nf = length(frequencies)
    flux = zeros(n, 2nb*nf)
    rate = zeros(n, 2nb*nf)
    waves = zeros(nl2, npre, 2nb*nf)
    states = zeros(nzs, 2nb*nf)
    tpre = [t0 - (npre - j)*sys.h for j in 1:npre]
    for (f, frequency) in enumerate(frequencies)
        w = 2pi*frequency
        F = lu(stationaryoperator(sys, x0, w))
        Ga = nzs > 0 ? incidentmap(sys, w, n + nl2) : nothing
        for b in 1:nb
            # a unit cosine current `cos(w (t - reference))` at the bath: the
            # complex amplitude `exp(-i w reference)`, and the sine `-i`
            # times it, whose response is `-i` times the cosine's; at time
            # t the response is Re(x exp(i w t)) and its rate
            # Re(i w x exp(i w t))
            rhs = zeros(ComplexF64, n + nl2)
            rhs[1:n] .= cispi(-2frequency*reference) .* Vector(injection[:, b])
            y = F \ rhs
            x = y[1:n]
            z = x .* cispi(2frequency*t0)
            col = 2*((f - 1)*nb + b - 1)
            flux[:, col + 1] .= real.(z)
            rate[:, col + 1] .= real.(im*w .* z)
            flux[:, col + 2] .= real.(-im .* z)
            rate[:, col + 2] .= real.(w .* z)
            for j in 1:npre
                a = y[n + 1:n + nl2] .* cispi(2frequency*tpre[j])
                waves[:, j, col + 1] .= real.(a)
                waves[:, j, col + 2] .= real.(-im .* a)
            end
            if nzs > 0
                zc = statephasors(p, w, Ga*y) .* cispi(2frequency*t0)
                states[:, col + 1] .= real.(zc)
                states[:, col + 2] .= real.(-im .* zc)
            end
        end
    end
    all(isfinite, flux) && all(isfinite, rate) && all(isfinite, waves) && all(isfinite, states) || throw(ArgumentError(
        "the stationary response of a bath is singular: an undamped mode at that frequency needs an explicit initial state."))
    return flux, rate, waves, states
end

# The stationary initial term of the adjoint method: the contraction of
# the adjoint's initial flux and rate with the stationary responses of
# every bath at every frequency, without building those responses. For an
# objective `o` and a bath `b` at `w`, with `z = F^{-1} inj_b exp(i w (t0 - reference))`
# and `F = -w^2 C + i w G + L + J'(x0)` symmetric, the term is
# `Re((lambda_x + i w lambda_v)' z)` for the cosine quadrature and the
# imaginary part for the sine, so one solve of `F` per frequency and
# objective, `y = F^{-1} (lambda_x + i w lambda_v)`, and the sparse
# products `inj' y` give every bath. Returns the terms in the column order
# of the responses, (frequency, bath, cosine/sine).
function stationaryinitialterm(sys::TransientSystem, x0, injection, frequencies, t0, reference, initialflux, initialrate,
        initialwaves = nothing, initialstates = nothing)
    return stationaryinitialterms(sys, reshape(x0, :, 1), injection, frequencies, t0, reference,
        reshape(initialflux, size(initialflux, 1), :, 1), reshape(initialrate, size(initialrate, 1), :, 1),
        isnothing(initialwaves) ? nothing : reshape(initialwaves, size(initialwaves, 1), size(initialwaves, 2), :, 1),
        isnothing(initialstates) ? nothing : reshape(initialstates, size(initialstates, 1), :, 1))[1]
end

# the terms of every condition, `(n, m, N)` initial fluxes and rates, one
# factorization per frequency and group of conditions sharing an initial
# state, solved for the objectives of the whole group together
function stationaryinitialterms(sys::TransientSystem, x0s, injection, frequencies, t0, reference, initialflux, initialrate,
        initialwaves = nothing, initialstates = nothing)
    p = sys.problem
    n = length(p)
    nl2 = 2length(p.lines)
    nzs = blockstates(p)
    npre = lineprehistory(p, sys.h)
    injection = (sys.Lscale/phi0) .* injection
    injectiont = sparse(transpose(injection))
    nb, nf = size(injection, 2), length(frequencies)
    lx, lv = Array(initialflux), Array(initialrate)
    la = nl2 > 0 ? Array(initialwaves) : zeros(0, npre, size(lx, 2), size(lx, 3))
    lz = nzs > 0 ? Array(initialstates) : zeros(0, size(lx, 2), size(lx, 3))
    m, N = size(lx, 2), size(lx, 3)
    terms = [zeros(m, 2nb*nf) for _ in 1:N]
    for group in initialstategroups(x0s)
        LX = reshape(view(lx, :, :, group), :, m*length(group))
        LV = reshape(view(lv, :, :, group), :, m*length(group))
        LA = reshape(view(la, :, :, :, group), nl2, npre, m*length(group))
        LZ = reshape(view(lz, :, :, group), nzs, m*length(group))
        for (f, frequency) in enumerate(frequencies)
            w = 2pi*frequency
            # the adjoint's stationary solve is the transposed system, the
            # same as the system where the operator is symmetric; the
            # prehistory's cotangents enter on the wave rows with the
            # phase of each sample's time relative to the start
            F = lu(stationaryoperator(sys, view(x0s, :, first(group)), w))
            rhs = zeros(ComplexF64, n + nl2, m*length(group))
            rhs[1:n, :] .= LX .+ (im*w) .* LV
            for j in 1:npre
                rhs[n + 1:n + nl2, :] .+= cispi(-2frequency*(npre - j)*sys.h) .* view(LA, :, j, :)
            end
            if nzs > 0
                # the states' cotangents through the transposed state map
                # and the transposed incident wave map onto the solution
                Ga = incidentmap(sys, w, n + nl2)
                abar = zeros(ComplexF64, size(Ga, 1), size(rhs, 2))
                prow = 0
                for b in p.blocks
                    np, nz = length(b.signal), size(b.A, 1)
                    nz > 0 && (abar[prow + 1:prow + np, :] .= transpose(b.B)*(transpose(im*w*I - b.A) \ view(LZ, b.zbase + 1:b.zbase + nz, :)))
                    prow += np
                end
                rhs .+= transpose(Ga)*abar
            end
            Y = (transpose(F) \ rhs)[1:n, :]
            P = (injectiont*Y) .* cispi(2frequency*(t0 - reference))
            all(isfinite, P) || throw(ArgumentError(
                "the stationary response of a bath is singular: an undamped mode at that frequency needs an explicit initial state."))
            for (jj, j) in enumerate(group), b in 1:nb
                col = 2*((f - 1)*nb + b - 1)
                cols = (jj - 1)*m + 1:jj*m
                terms[j][:, col + 1] .= real.(view(P, b, cols))
                terms[j][:, col + 2] .= imag.(view(P, b, cols))
            end
        end
    end
    return terms
end

# the Fourier bins of an input plan the bath frequencies fall on, for the
# incremental gain from the input modes; the bath is then the periodic
# Fourier mode of the record
function noiseinputbins(inputs, measurement, frequencies, weights)
    isnothing(inputs) && return Int[]
    inputs.times == measurement.times || throw(ArgumentError("the input and output modes must share the record."))
    period = length(inputs.times)*inputs.dt
    all(w -> isapprox(w, inv(period); rtol = 1e-10), weights) || throw(ArgumentError(
        "the gain from input modes needs the Fourier bath weights 1/(N dt) of the record."))
    bins = round.(Int, frequencies .* period)
    all(k -> 1 <= k <= length(inputs.frequencies), bins) &&
        all(j -> isapprox(frequencies[j], inputs.frequencies[bins[j]]; rtol = 1e-10), eachindex(bins)) ||
        throw(ArgumentError("the gain from input modes needs bath frequencies on the Fourier bins of the record."))
    all(j -> isapprox(sum(abs2, view(inputs.coefficients, bins, j)), 1; rtol = 1e-10), eachindex(inputs.ports)) ||
        throw(ArgumentError("the bath grid must cover every input mode."))
    return bins
end

function accumulatenoisegain!(gain, response, inputs, bins, frequencies, baths)
    isnothing(inputs) && return nothing
    nb, nf = length(baths), length(frequencies)
    mixing = zeros(size(response, 2), 2length(inputs.ports))
    for f in 1:nf, b in 1:nb, k in eachindex(inputs.ports)
        baths.channels[b].port == inputs.ports[k] || continue
        c = inputs.coefficients[bins[f], k]
        j = (f - 1)*nb + b
        mixing[2j - 1, 2k - 1] = mixing[2j, 2k] = real(c)
        mixing[2j, 2k - 1] = imag(c)
        mixing[2j - 1, 2k] = -imag(c)
    end
    mul!(gain, response, tobackend(KernelAbstractions.get_backend(gain), mixing), 1.0, 1.0)
    return nothing
end

# The memory the noise contraction may take on its backend, in bytes: a
# quarter of the free memory when the reference holds zero, which it
# does unless a test sets it, so that the tiling of the contraction can
# be exercised on a small problem.
const noisememorybudget = Ref(0)
noisebudget(backend) = noisememorybudget[] > 0 ? noisememorybudget[] : freememory(backend) ÷ 4

# The classical equilibrium the noise starts from. The state must be a
# fixed point of the circuit under the drive at that time: the drift
# balances the drive, no inductive flux and no junction phase moves, the
# drive did not change just before, and every junction sits at a phase of
# positive differential inductance, so that the passive linearization the
# stationary response uses is stable. A node with only conductance and
# capacitance may carry a constant voltage, since its flux never enters
# an equation. The prehistory of an arbitrary waveform cannot be proved
# from its values, so that the drive was constant before the record
# remains the user's contract.
function transientstationary(sys::TransientSystem, x, v, t, p::TransientProblem = sys.problem, waves = zeros(2length(p.lines)),
        states = zeros(blockstates(p)))
    G, L, RJ = hostsparse(sys.G), hostsparse(sys.L), hostsparse(sys.RJ)
    lmolj = Array(sys.lmolj)
    xh, vh = Array(x), Array(v)
    # the drive with the lines' arriving waves and the blocks' resting
    # waves, the rest of the equilibrium
    linevalues = lineforcing(p, arrivingwaves(p, Array(waves)))
    resting = blockstates(p) > 0 ? restingwaves(sys, reshape(Array(states), :, 1))[:, 1] : nothing
    b = hostdrivecurrent(sys, t, p, linevalues, resting)
    phi = RJ*xh
    scale = max(norm(b, Inf), 1.0)
    r = G*vh .+ L*xh .+ transpose(RJ)*(lmolj .* sin.(phi)) .- b
    # a rational block's states are at rest under the incident waves of
    # the state: `A z + B a = 0`
    for bl in p.blocks
        nz = size(bl.A, 1)
        nz == 0 && continue
        rate = [(bl.signal[q] > 0 ? vh[bl.signal[q]] : 0.0) - (bl.ref[q] > 0 ? vh[bl.ref[q]] : 0.0) for q in eachindex(bl.signal)]
        current = xh[bl.auxbase + 1:bl.auxbase + length(bl.signal)] .* (phi0/sys.Lscale)
        a = (phi0 .* rate ./ sqrt.(bl.R) .+ sqrt.(bl.R) .* current) ./ 2
        z = Array(states)[bl.zbase + 1:bl.zbase + nz]
        drift = bl.A*z .+ bl.B*a
        norm(drift, Inf) <= 1e-6*max(norm(bl.A*z, Inf), norm(bl.B*a, Inf), 1e-30) + 1e-12*norm(bl.A, Inf)*max(norm(z, Inf), 1.0) ||
            throw(ArgumentError(lazy"the noise needs a classical equilibrium at the start of the record: the states of the scattering block at $(bl.path) are not at rest under the incident waves there; start the solve before the drive."))
    end
    norm(r, Inf) <= 1e-6*scale || throw(ArgumentError(
        "the noise needs a classical equilibrium at the start of the record, where the stationary response of the baths is the circuit's prehistory: the drift does not balance the drive there; start the solve before the drive."))
    motion = max(norm(L*vh, Inf), norm(RJ*vh, Inf))
    motion <= 1e-6*max(norm(vh, Inf), 1.0) || throw(ArgumentError(
        "the noise needs a classical equilibrium at the start of the record: an inductive flux or a junction phase is moving there, so the circuit is not stationary; start the solve before the drive."))
    before = hostdrivecurrent(sys, t - sys.h, p, linevalues, resting)
    norm(before .- b, Inf) <= 1e-6*scale || throw(ArgumentError(
        "the noise needs a constant drive before the start of the record, whose stationary state is the circuit's prehistory; start the solve before the drive changes."))
    all(>(0), cos.(phi)) || throw(ArgumentError(
        "a junction starts beyond a quarter flux quantum, where its differential inductance is negative and the stationary linearization is not a stable prehistory."))
    return nothing
end

"""
    transientnoise(solution, measurement; frequencies, weights,
        baths = transientnoisebaths(solution.problem), method = :adjoint,
        inputs = nothing, commutationrtol = 1e-3, reuse = nothing)

The symmetrized quantum noise of a recorded transient, or of every
condition of a [`TransientBatchSolution`](@ref) on one pass, in the temporal
modes of `measurement`, a [`TransientQuantumPlan`](@ref) on a window of
the recorded times: the physical baths, the port terminations and the internal
resistors, propagated through the linearization about the complete
recorded trajectory, pump and signals together. For a batch the
covariance, the commutator and the gain carry the conditions as the
trailing dimension and the diagnostics are a vector. `frequencies` are the
positive nodes in Hz of a quadrature over the bath spectrum and `weights`
its weights in Hz; for a periodic bath over the record of duration `T`
use `f = k/T` and weights `1/T`. Each bath at each frequency is a pair of
cosine and sine Norton currents of amplitude `2 sqrt(h f df/R)`, whose
independent quadratures have variance `nbar + 1/2`, started from the
stationary response of the circuit to them at the initial state, so
that the fluctuations stored before the record and their correlation
with the forcing are kept; the trajectory must therefore start at a
classical equilibrium under a constant drive, with no inductive flux or
junction phase moving, and the measurement window may begin after the
drive has settled. With an input plan the bath is the periodic Fourier
bath of the measurement window, its tones extended over the whole record
and started stationary, and the gain is the response to those periodic
modes; a probe applied only inside the window is a different quantity,
to be computed by driving that waveform.

`method = :forward` propagates the bath quadratures forward, two tangent
directions per bath and frequency through [`transienttangent`](@ref);
`method = :adjoint` propagates the measured quadratures backward through
[`transientadjoint`](@ref), whose derivatives with respect to the bath
currents at every recorded time are the bath kernels, contracted against
the frequencies as a discrete Fourier sum at the grid and the stage
times the adjoint hands over, so the
frequency count costs sums rather than integrations. Both return the
same `covariance`, `commutator`, `expectedcommutator`, `diagnostics`
(from [`transientquantumdiagnostics`](@ref)), and, with an input plan
`inputs` on the same record, the incremental quadrature `gain` from its
modes. A failed diagnostic is returned as such, not rescaled away; refine
the bath cutoff, the frequency spacing, the record and the step
independently before reading a quantum efficiency with
[`transientquantumefficiency`](@ref).
"""
function transientnoise(sol::Union{TransientSolution,TransientBatchSolution}, measurement::TransientQuantumPlan;
        frequencies, weights, baths = nothing,
        method::Symbol = :adjoint, inputs = nothing, commutationrtol = 1e-3,
        reuse = nothing)
    recordedsolution(sol)
    problems, _, x0s, v0s = batchview(sol)
    N = length(problems)
    p = first(problems)
    baths = isnothing(baths) ? transientnoisebaths(p) : baths
    baths.problem.circuit === p.circuit || throw(ArgumentError("the baths belong to another circuit."))
    method in (:forward, :adjoint) || throw(ArgumentError("method must be :forward or :adjoint."))
    isnothing(inputs) || inputs isa TransientQuantumPlan || throw(ArgumentError("inputs must be a TransientQuantumPlan."))
    fs, ws = Float64.(collect(frequencies)), Float64.(collect(weights))
    !isempty(fs) && length(fs) == length(ws) && all(x -> isfinite(x) && x > 0, fs) &&
        all(x -> isfinite(x) && x > 0, ws) && issorted(fs) && allunique(fs) || throw(ArgumentError(
        "give distinct increasing positive bath frequencies and positive weights in Hz."))
    backend = KernelAbstractions.get_backend(sol.finalflux)
    KernelAbstractions.get_backend(measurement.weights) == backend || throw(ArgumentError(
        "the measurement and the solution must share a backend."))
    # the measurement's record is a window of the recorded times, whose
    # start is the phase reference of the bath quadratures; the baths are
    # driven from the record's start, which must be a classical equilibrium
    nm = length(measurement.times)
    offset = windowoffset(sol, measurement)
    reference = first(measurement.times)
    maximum(fs) < 0.5/measurement.dt || throw(ArgumentError("the bath cutoff must lie below the Nyquist frequency of the record."))
    np = length(p.portimpedances)
    maximum(measurement.ports) <= np || throw(DimensionMismatch("a measurement port is absent from the circuit."))
    bins = noiseinputbins(inputs, measurement, fs, ws)
    fact = transientfactorization(backend)
    sys = transientsystem(reuse, p, sol.dt, sol.method, backend, fact)
    # the classical equilibrium at the start of every condition, whose
    # stationary response to the baths is the circuit's prehistory
    w0s, z0s = batchview(sol)[7], batchview(sol)[8]
    for j in 1:N
        transientstationary(sys, view(x0s, :, j), view(v0s, :, j), first(sol.times), problems[j], view(w0s, :, j), view(z0s, :, j))
    end
    nb, nf, nt = length(baths), length(fs), length(sol.times)
    m = 2length(measurement.ports)
    stageoffsets = stagetimeoffsets(sol)
    injection = transientinjection(p, baths)
    amplitudes = [bathamplitude(b, fs[f], ws[f]) for b in baths.channels, f in 1:nf]
    variances = [thermaloccupation(2pi*fs[f], baths.channels[b].temperature)/2 for f in 1:nf for b in 1:nb for _ in 1:2]
    ngain = isnothing(inputs) ? 0 : 2length(inputs.ports)
    covariance, commutator, gain = zeros(m, m, N), zeros(m, m, N), zeros(m, ngain, N)
    if method == :forward
        # the responses of every quadrature pair in the measured modes of
        # every condition, as the columns of one matrix per condition in
        # the order (frequency, bath, cosine/sine)
        response = zeros(m, 2nb*nf, N)
        correction = isempty(baths.groups) ? nothing : groupcorrection(baths, fs, m, CPU())
        # the bath quadratures as the directions of one tangent over the
        # batch, the currents shared and the stationary response of each
        # condition's initial state, for unit currents scaled by the
        # amplitudes, as its initial perturbation
        scale = [amplitudes[b, f] for f in 1:nf for b in 1:nb for _ in 1:2]
        ndir = 2nb*nf
        # the quadratures at the grid and at the stage times of each step
        currents = zeros(nb, 3, nt, ndir)
        for f in 1:nf, b in 1:nb, s in 1:3
            col = 2*((f - 1)*nb + b - 1)
            ts = sol.times .+ stageoffsets[s] .- reference
            currents[b, s, :, col + 1] .= amplitudes[b, f] .* cospi.(2fs[f] .* ts)
            currents[b, s, :, col + 2] .= amplitudes[b, f] .* sinpi.(2fs[f] .* ts)
        end
        n = length(p)
        npre = lineprehistory(p, sys.h)
        flux0s, rate0s, waves0s = zeros(n, ndir, N), zeros(n, ndir, N), zeros(2length(p.lines), npre, ndir, N)
        states0s = zeros(blockstates(p), ndir, N)
        for j in 1:N
            flux0, rate0, waves0, states0 = stationaryresponses(sys, view(x0s, :, j), injection, fs, first(sol.times), reference)
            flux0s[:, :, j] .= flux0 .* transpose(scale)
            rate0s[:, :, j] .= rate0 .* transpose(scale)
            waves0s[:, :, :, j] .= waves0 .* reshape(scale, 1, 1, :)
            states0s[:, :, j] .= states0 .* transpose(scale)
        end
        # the responses measured as the tangent produces them
        measured, outputsink = measurementsink(measurement, offset, ndir*N, backend)
        transienttangent(sol, currents; targets = baths, initialstate = (flux0s, rate0s, waves0s, states0s), reuse, outputsink)
        response .= reshape(Array(measured), m, ndir, N)
        for j in 1:N
            if !isnothing(correction)
                groupcorrection!(view(covariance, :, :, j), view(commutator, :, :, j), view(response, :, :, j), correction)
                maskgroups!(view(response, :, :, j), correction)
            end
            transientnoiseaccumulate!(view(covariance, :, :, j), view(commutator, :, :, j), view(response, :, :, j), variances)
            accumulatenoisegain!(view(gain, :, :, j), view(response, :, :, j), inputs, bins, fs, baths)
        end
    else
        # the bath kernels: the derivatives of the measured quadratures with
        # respect to a current at each bath and recorded time
        weightsout = zeros(np, nt, m)
        w = Array(measurement.weights)
        for (j, port) in enumerate(measurement.ports), q in 1:2
            weightsout[port, offset:offset + nm - 1, 2(j - 1) + q] .= w[:, 2(j - 1) + q]
        end
        # Each column of the bath kernels, as the adjoint hands it over
        # with the objectives of every condition, at a grid time or at a
        # stage time, is contracted with the cosine and the sine of every
        # frequency at that time into a (bath, objective, condition) by
        # (frequency, quadrature) accumulator on the backend, so no kernel
        # is stored: the columns are buffered over a block of times and
        # contracted by one product, so the accumulator is read and
        # written once per block rather than once per time. The
        # accumulator of all conditions and frequencies at once is
        # `16 nb m nf N` bytes; within a quarter of the backend's free
        # memory it is one, beyond it the conditions are tiled, each tile
        # one adjoint over its conditions, and then the frequencies, each
        # tile one more adjoint over the same conditions, a trade of
        # passes for memory. The covariance, the commutator and the gain
        # are accumulated from each tile on the backend, so the host never
        # holds a kernel or a response.
        budget = noisebudget(backend)
        percondition = 16*nb*m*nf
        ctile = clamp(budget ÷ percondition, 1, N)
        dcov, dcomm = [KernelAbstractions.zeros(backend, Float64, m, m, N) for _ in 1:2]
        dgain = KernelAbstractions.zeros(backend, Float64, m, ngain, N)
        # the frequency tiles outermost, so that the loss matrices of the
        # blocks' groups over a tile are evaluated once for every tile of
        # conditions
        nft = clamp(budget ÷ (16*nb*m*ctile), 1, nf)
        for ftile in [f:min(f + nft - 1, nf) for f in 1:nft:nf]
            correction = isempty(baths.groups) ? nothing : groupcorrection(baths, fs[ftile], m, backend)
            for conditions in [j:min(j + ctile - 1, N) for j in 1:ctile:N]
                sub = length(conditions) == N ? sol : sol[conditions]
                noisetile!(view(dcov, :, :, conditions), view(dcomm, :, :, conditions), view(dgain, :, :, conditions),
                    sub, sys, view(x0s, :, conditions), weightsout, baths, injection, fs, ftile, amplitudes,
                    reference, stageoffsets, inputs, bins, reuse, correction)
            end
        end
        covariance .= Array(dcov)
        commutator .= Array(dcomm)
        gain .= Array(dgain)
    end
    expected = copy(measurement.commutator)
    diagnostics = map(1:N) do j
        transientquantumdiagnostics(covariance[:, :, j], commutator[:, :, j], expected; rtol = commutationrtol)
    end
    if sol isa TransientSolution
        return (; covariance = covariance[:, :, 1], commutator = commutator[:, :, 1], expectedcommutator = expected,
            diagnostics = diagnostics[1], gain = gain[:, :, 1], measurement, inputs, baths, frequencies = fs, weights = ws)
    end
    return (; covariance, commutator, expectedcommutator = expected, diagnostics, gain,
        measurement, inputs, baths, frequencies = fs, weights = ws)
end

# One tile of the noise contraction: the adjoint of the conditions of
# `sub` with the measured quadratures as objectives, its columns
# contracted against the frequencies of the tile at the grid and the
# stage times through a block buffer, the stationary initial terms of
# the conditions through the adjoint's initial flux and rate, and the
# responses formed on the backend, scaled by the amplitudes, and
# accumulated into the covariance, the commutator and the gain.
function noisetile!(covariance, commutator, gain, sub, sys, x0s, weightsout, baths, injection, fs, ftile, amplitudes,
        reference, stageoffsets, inputs, bins, reuse, correction)
    backend = KernelAbstractions.get_backend(covariance)
    nb, m, Nt = length(baths), size(covariance, 1), size(covariance, 3)
    fsl = fs[ftile]
    nfl = length(fsl)
    fsb = tobackend(backend, fsl)
    accumulator = KernelAbstractions.zeros(backend, Float64, nb*m*Nt, 2nfl)
    # the block of columns and their quadratures, contracted when full
    block = 64
    columns = KernelAbstractions.zeros(backend, Float64, nb*m*Nt, block)
    quadratures = KernelAbstractions.zeros(backend, Float64, block, 2nfl)
    filled = Ref(0)
    flush! = () -> begin
        filled[] == 0 && return nothing
        mul!(accumulator, view(columns, :, 1:filled[]), view(quadratures, 1:filled[], :), 1.0, 1.0)
        filled[] = 0
        nothing
    end
    contract! = (t, values) -> begin
        filled[] == block && flush!()
        filled[] += 1
        copyto!(view(columns, :, filled[]), vec(values))
        view(quadratures, filled[], 1:nfl) .= cospi.(2 .* fsb .* t)
        view(quadratures, filled[], nfl + 1:2nfl) .= sinpi.(2 .* fsb .* t)
        nothing
    end
    sink = (k, values) -> contract!(sub.times[k] - reference, values)
    stagesink = (k, i, values) -> contract!(sub.times[k] + stageoffsets[i + 1] - reference, values)
    adjoint = transientadjoint(sub, weightsout; quantity = :outgoing, targets = baths, reuse, sink, stagesink)
    flush!()
    # the stationary initial terms of every condition of the tile
    initials = stationaryinitialterms(sys, x0s, injection, fsl, first(sub.times), reference,
        reshape(adjoint.initialflux, :, m, Nt), reshape(adjoint.initialrate, :, m, Nt),
        isnothing(adjoint.initialwaves) ? nothing : reshape(adjoint.initialwaves, size(adjoint.initialwaves, 1), size(adjoint.initialwaves, 2), m, Nt),
        isnothing(adjoint.initialstates) ? nothing : reshape(adjoint.initialstates, size(adjoint.initialstates, 1), m, Nt))
    acc = reshape(accumulator, nb, m, Nt, 2nfl)
    amp = tobackend(backend, amplitudes[:, ftile])
    variances = tobackend(backend, [thermaloccupation(2pi*fsl[f], baths.channels[b].temperature)/2 for f in 1:nfl for b in 1:nb for _ in 1:2])
    response = KernelAbstractions.zeros(backend, Float64, m, 2, nb, nfl)
    flat = reshape(response, m, 2nb*nfl)
    for j in 1:Nt
        # the (objective, quadrature, bath, frequency) response of the
        # condition from the accumulator's cosine and sine halves and the
        # initial term, times the bath amplitudes
        view(response, :, 1, :, :) .= permutedims(view(acc, :, :, j, 1:nfl), (2, 1, 3))
        view(response, :, 2, :, :) .= permutedims(view(acc, :, :, j, nfl + 1:2nfl), (2, 1, 3))
        flat .+= tobackend(backend, initials[j])
        response .*= reshape(amp, 1, 1, nb, nfl)
        accumulatenoisegain!(view(gain, :, :, j), flat, inputs, isempty(bins) ? bins : bins[ftile], fsl, baths)
        if !isempty(baths.groups)
            groupcorrection!(view(covariance, :, :, j), view(commutator, :, :, j), flat, correction)
            maskgroups!(flat, correction)
        end
        transientnoiseaccumulate!(view(covariance, :, :, j), view(commutator, :, :, j), flat, variances)
    end
    return nothing
end

# the sink of a tangent that measures the outgoing waves of every column
# in the modes of a plan on its window, accumulating the `(2 nports,
# columns)` quadratures on the backend, and that accumulator
function measurementsink(measurement::TransientQuantumPlan, offset, ncolumns, backend)
    m = 2length(measurement.ports)
    nm = length(measurement.times)
    w = Array(measurement.weights)
    measured = KernelAbstractions.zeros(backend, Float64, m, ncolumns)
    sink = (k, voltage, incident, outgoing) -> begin
        offset <= k <= offset + nm - 1 || return nothing
        kk = k - offset + 1
        for (j, port) in enumerate(measurement.ports), q in 1:2
            view(measured, 2(j - 1) + q, :) .+= w[kk, 2(j - 1) + q] .* view(outgoing, port, :)
        end
        nothing
    end
    return measured, sink
end

# the offsets from a recorded time of the three times a staged current
# holds, the grid time and the two stage times of the step from it
stagetimeoffsets(sol) = sol.method isa GaussLegendre ? (0.0, sol.dt*gausscoefficients().c[1], sol.dt*gausscoefficients().c[2]) : (0.0, 0.0, 0.0)

# the index of a plan's window in a solution's times, the plan's record
# being a contiguous window of them
function windowoffset(sol, plan::TransientQuantumPlan)
    nm = length(plan.times)
    offset = findfirst(t -> isapprox(t, first(plan.times); atol = 1e-6*sol.dt), sol.times)
    (!isnothing(offset) && offset + nm - 1 <= length(sol.times) &&
        isapprox(view(sol.times, offset:offset + nm - 1), plan.times; atol = 1e-6*sol.dt)) ||
        throw(ArgumentError("a plan's record must be a window of the recorded times of the solution."))
    return offset
end

"""
    transientgain(solution, measurement, inputs; reuse = nothing)

The quadrature gain from the temporal modes of `inputs` to those of
`measurement`, both [`TransientQuantumPlan`](@ref)s on windows of the
record, for a probe applied inside the input window only: each input
mode's two quadratures are the incident waves of unit `X` and unit `P`
at the input's port, evaluated at the grid and the stage times inside
the window and zero outside it, driven through [`transienttangent`](@ref)
and read in the measurement's modes. Returns
the `(2 nout, 2 nin)` matrix, with the conditions as a trailing
dimension for a [`TransientBatchSolution`](@ref). This is the causal,
pulsed gain, carrying the transients of the probe's own edges; the
`gain` of [`transientnoise`](@ref) with an input plan is the response to
the periodic Fourier mode of the window, extended over the record and
started stationary, which is what a stationary amplifier's harmonic
balance gain is. The two agree as the window grows past the circuit's
memory.
"""
function transientgain(sol::Union{TransientSolution,TransientBatchSolution}, measurement::TransientQuantumPlan,
        inputs::TransientQuantumPlan; reuse = nothing)
    recordedsolution(sol)
    problems, _, _, _ = batchview(sol)
    N = length(problems)
    p = first(problems)
    backend = KernelAbstractions.get_backend(sol.finalflux)
    (KernelAbstractions.get_backend(measurement.weights) == backend && KernelAbstractions.get_backend(inputs.weights) == backend) ||
        throw(ArgumentError("the plans and the solution must share a backend."))
    np = length(p.portimpedances)
    (maximum(measurement.ports) <= np && maximum(inputs.ports) <= np) || throw(DimensionMismatch("a plan's port is absent from the circuit."))
    mo, io = windowoffset(sol, measurement), windowoffset(sol, inputs)
    nm, ni, nt = length(measurement.times), length(inputs.times), length(sol.times)
    nin, m = length(inputs.ports), 2length(measurement.ports)
    # the incident wave of a unit quadrature of each input mode on its
    # window, `Re` and `-Im` of the mode's spectrum with the bin
    # amplitudes `sqrt(h f / T)`, brought to the grid and, with the phase
    # of each bin advanced to the stage time, to the stages of every step
    # inside the window, as the drive currents `2 a / sqrt(Z)` at its port
    T = ni*inputs.dt
    currents = zeros(np, 3, nt, 2nin)
    spectrum = zeros(ComplexF64, ni)
    offsets = stagetimeoffsets(sol)
    for j in 1:nin, s in 1:3
        fill!(spectrum, 0)
        spectrum[2:length(inputs.frequencies) + 1] .= sqrt.(planck_constant .* inputs.frequencies ./ T) .*
            view(inputs.coefficients, :, j) .* cispi.(-2 .* inputs.frequencies .* offsets[s])
        wave = FFTW.fft(spectrum)
        port = inputs.ports[j]
        scale = 2/sqrt(p.portimpedances[port])
        last = s == 1 ? io + ni - 1 : min(io + ni - 1, nt - 1)
        currents[port, s, io:last, 2j - 1] .= scale .* real.(view(wave, 1:last - io + 1))
        currents[port, s, io:last, 2j] .= -scale .* imag.(view(wave, 1:last - io + 1))
    end
    # the modes measured as the tangent produces the waves
    measured, outputsink = measurementsink(measurement, mo, 2nin*N, backend)
    transienttangent(sol, currents; reuse, outputsink)
    gain = reshape(Array(measured), m, 2nin, N)
    return sol isa TransientSolution ? gain[:, :, 1] : gain
end
