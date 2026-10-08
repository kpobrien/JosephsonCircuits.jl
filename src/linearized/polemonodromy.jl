# The period map of a pumped circuit: every mode once, from the map of its
# linearized dynamics over one pump period, which the transient's tangent
# carries along one recorded period of the orbit.

"""
    Monodromy(; steps = 64, nev = 10, backend = CPU())

Find the least damped modes of a periodic operating point with
[`hbstability`](@ref). The method forms a dense period map, computes all
of its multipliers, and returns up to `nev` numerically resolved modes
(**10 by default**). `nev = :all` requests all resolved modes. Each
selected mode has one frequency representative rather than the aliases
returned by a harmonic polynomial search.

The HB orbit initializes [`transientsolve`](@ref), which records one pump
period `T` in `steps` Gauss–Legendre steps. [`transienttangent`](@ref)
then propagates the constrained circuit perturbations: node fluxes and
rates, rational-block states, and outgoing-wave histories of transmission
lines. The map satisfies `δx(T) = M*δx(0)`; its eigenvalues are multipliers
`μ = exp(s*T)`. Algebraic auxiliary currents are constrained, and uniform
flux gauge directions are removed. This is a circuit state/history map,
not a scattering matrix. Pumped blocks convert fully on, independent of
their startup `envelope`, as they do in harmonic balance.

For each selected eigenvector, propagation and Fourier analysis of its
periodic part choose a dominant-frequency representative of `log(μ)/T`
and supply a mode profile. `Nmodulationharmonics = (H,)` sets the width
of that profile only; it does not truncate the map eigenproblem. Require
`2H + 1 <= steps`. Frequencies outside the temporal sampling bandwidth
can be placed at aliases; increase `steps` to resolve them. A line adds
about `2*delay/(T/steps)` history coordinates, and its delay must be at
least one step. A mode with no node or junction content, internal to
rational blocks or to lines, takes its representative from its block
states and the waves leaving its lines' ports.

Refine both the HB orbit's harmonics and `steps`. `rateerrors` estimates
each rate's timestep error from the selected modes propagated through a
period of twice the steps: the refined map projected onto them, which
follows their mixing among themselves, gives each mode's change of rate,
scaled by 16/15 to the coarse map's error at the rule's fourth order. It
is not a certified bound: mixing with modes beyond `nev` is not followed,
and a map too coarse for its error to fall as the fourth power of the
step can have more error than its estimate. Compare solves at increasing
`steps`, and at a larger `nev` where modes crowd, before deciding the
sign of a small rate. The Gauss–Legendre rule is not L-stable: decay much
faster than a step can appear as slow discrete decay. Resolve the
circuit's fast decays as well as its high frequencies. Delay-history
interpolation also adds error.

`backend` selects propagation and the dense eigensolve:

- `CPU()` computes all multipliers through a real Schur decomposition,
  then forms only the selected left/right vectors on the session's BLAS
  threads.
- `CUDABackend()` requires CUDA.jl and CUDSS.jl, and a CUDA runtime
  with cuSOLVER 11.7.1 or newer for its `geev`, which CUDA.jl's wrapper
  checks at the eigensolve. The device computes all right vectors and
  obtains selected left vectors from their basis, each checked against
  its residual; a singular or ill-conditioned basis falls back to the
  CPU's Schur vectors. Balancing, selection, orbit initialization, and
  constraint setup stay on the host.

For map dimension `d`, storage is `O(d^2)` and the dense eigensolve is
`O(d^3)` even for small `nev`. Reducing `nev` saves selected-vector and
profile work. More steps also enlarge the map when it contains delays.

Without a pump this method uses [`DenseSpectrum`](@ref)`()` instead,
with its default 1000-unknown limit. Exact unpumped delays and any
[`LaplaceResponse`](@ref) model require [`ContourIntegral`](@ref).
A frequency-domain block without a supported transient realization
requires an appropriate polynomial or contour method.
"""
struct Monodromy{B<:Backend}
    steps::Int
    nev::Int
    backend::B
end
function Monodromy(; steps::Integer = 64, nev::Union{Integer,Symbol} = 10, backend::Backend = CPU())
    steps >= 1 || throw(ArgumentError("steps must be positive."))
    count = nev === :all ? typemax(Int) : nev isa Integer ? Int(nev) : 0
    count >= 1 || throw(ArgumentError("nev is a positive integer or :all."))
    return Monodromy(Int(steps), count, backend)
end

# The poles of a pumped circuit from its period map (see Monodromy). `sys`
# is the one-harmonic pole system, which checks the circuit and the
# operating point and gives the pump and frequency scale and the junction
# branches; `H` is the half width of the profiles in harmonics; `chunk` is
# how many directions one tangent carries on the host, a few dozen at
# most, whose work stays near the processor's cache, which every direction
# at once overflows with many times the map's memory, and `devicechunk`
# how many on a device, whose kernels take many columns to fill it;
# `tracking` is the share of a mode which the finer map's mode continuing
# it carries, more than half, where the refinement leaves it the mode it
# was; `ruleorder` is the rule's order, by which the error of a rate
# the steps resolve falls as the steps double; `points` is how many samples
# of a line's history resample it at the finer step: the resampling's
# error biases the estimate where it nears the step's own error in a
# rate, which in a weakly damped mode is far below the step's error in
# its frequency, and twice the six of the lines' own read keep it below
# that in a mode the steps resolve.
function polemonodromy(method::Monodromy, psc::CompiledCircuit, circuitdefs, nonlinear,
        sys::HBPoleSystem, H::Int, factorization; chunk::Int = 16, devicechunk::Int = 1024,
        tracking::Float64 = 0.5, ruleorder::Int = 4, points::Int = 12)
    omega, scale = sys.pumpfrequency, sys.scale
    T = 2pi/omega
    N = method.steps
    # the profiles' harmonics are read from the samples of the period, of
    # which each needs its own
    2H + 1 <= N || throw(ArgumentError(
        lazy"the profiles' $(2H + 1) harmonics, Nmodulationharmonics = ($(H),), need as many samples of the period, and Monodromy(steps = $(N)) takes $(N); give steps >= $(2H + 1) or fewer harmonics."))
    # the transient of the operating point, its pumped blocks converting
    # as harmonic balance and the polynomial take them, always on
    p = alwayson(transientproblem(psc, circuitdefs; sources = isnothing(nonlinear) ? TransientSource[] : orbitsources(psc, nonlinear)))
    # a step reads a line's accepted history only
    for line in p.lines
        line.delay >= T/N || throw(ArgumentError(
            lazy"the transmission line at $(line.path) has a delay of $(line.delay) s, shorter than the step of a period of Monodromy(steps = $(N)); give steps >= $(ceil(Int, T/line.delay))."))
    end
    backend = method.backend
    host = backend isa CPU
    width = host ? chunk : devicechunk
    # The record of the period in `steps` steps from the orbit's state at
    # its start, its lines' history sampled at the step, made consistent
    # on the host's system of those steps, which reads the constraints;
    # with the system its solve and every tangent along it share, and
    # their factorization: the host's, or on a device the device's.
    function record(steps)
        shared = TransientReuse()
        stepped = transientsystem(shared, p, T/steps, GaussLegendre(), CPU(), factorization)
        state = isnothing(nonlinear) ? transientstate(p) :
            projectedstate(orbitstate(p, nonlinear; dt = T/steps), stepped, 0.0)
        fact = host ? factorization : transientfactorization(backend)
        recorded = host ? shared : TransientReuse()
        solved = transientsolve(p, (0.0, T); dt = T/steps, method = GaussLegendre(), record = :phases,
            initialstate = state, factorization = fact, backend, reuse = recorded)
        return (; constraints = stepped, system = host ? stepped : transientsystem(recorded, p, T/steps, GaussLegendre(), backend, fact),
            sol = solved, factorization = fact)
    end
    rec = record(N)
    n, nz = length(p), blockstates(p)
    nn, d = p.Nnodal, 2n + nz
    # the history the map holds, and its maps to the tangent's at the
    # step over `refinement`: none without lines, whose maps a circuit
    # without them then does not compile
    nl2 = 2length(p.lines)
    hrows, hcols = isempty(p.lines) ? (Int[], Int[]) : historycolumns(p, T/N, points)
    nw = length(hrows)
    D = d + nw
    maps(refinement) = nw == 0 ? nothing : Base.invokelatest(historymaps, p, hrows, hcols, T/N, refinement, points)
    history = maps(1)
    # The map acts on the dynamic coordinates: the node fluxes and rates,
    # the blocks' states and the lines' history. The auxiliary currents of
    # coupled inductors and of blocks, with their rates, are algebraic,
    # fixed by the others, and nothing reads their initial value. The
    # uniform flux of each inductively floating subnetwork, which no
    # equation reads either, is a direction the map keeps, of multiplier
    # one; the map on the quotient by those directions sets them aside:
    # each subnetwork's first node is not propagated, and its share is
    # taken out of every image.
    floating = [component .- 1 for component in
        calcstaticfluxcomponents(psc.componenttypes, psc.nodeindices, p.matrices.vvn, psc.Nnodes)]
    dynamic = vcat(setdiff(1:nn, first.(floating)), n .+ (1:nn), 2n .+ (1:nz), d .+ (1:nw))
    nq = length(p.portimpedances)
    # Along a linear algebraic direction which no junction and no source
    # touches, the rule keeps the constraint's residual `Zc L x` as it is
    # (see `InvariantRate`), so a starting flux off the constraint would
    # come back as a multiplier of one, a pole at zero of roundoff's sign.
    # The orbit and every mode about it meet the constraint, and the
    # starting fluxes are projected onto it; those directions map to zero.
    function constrained!(W, rec)
        ir = rec.constraints.invariant
        isnothing(ir) && return W
        X = W[1:n, :]
        readinvariantrate!(X, X, ir, invariantratework(ir, CPU(), n, size(W, 2)))
        W[1:n, :] .= X
        return W
    end
    # the tangent of the directions `W` along the record `rec`, its final
    # states stacked on the host, on the record's system, which each
    # thread's reuse shares, the history taken from the coordinates and
    # returned to them through `history` (see `historymaps`), its states
    # handed to `statesink` if one is given; the port responses, which
    # nothing reads, go to a sink which drops them rather than into stored
    # histories
    discard(k, voltage, incident, outgoing) = nothing
    function carried(W, rec, reuse, history, statesink)
        ncol = size(W, 2)
        waves = nw == 0 ? nothing : reshape(history.into*W[d+1:D, :], nl2, history.npre, ncol)
        t = transienttangent(rec.sol, zeros(nq, length(rec.sol.times), ncol);
            initialstate = (W[1:n, :], W[n+1:2n, :], waves, W[2n+1:d, :]), factorization = rec.factorization,
            reuse, outputsink = discard, statesink)
        return vcat(Array(t.finalflux), Array(t.finalrate),
            isnothing(t.finalstates) ? zeros(0, ncol) : Array(t.finalstates),
            nw == 0 ? zeros(0, ncol) : history.from*reshape(Array(t.finalwaves), :, ncol))
    end
    # The images of directions in the dynamic coordinates, `columns(ch)`
    # those of the columns `ch` of `out`, carried along the record `rec`
    # into them with their gauge shares taken out: the directions are
    # independent, and are carried in chunks of `width`, the same whatever
    # the threads, which share them out, each thread on a workspace of its
    # own from chunk to chunk.
    function images!(out, columns, rec, history)
        slices = chunkranges(size(out, 2), width)
        runchunks(batchchunks(backend, length(slices))) do _, range
            reuse = TransientReuse(rec.system, nothing, nothing, nothing, nothing, nothing)
            for ch in view(slices, range)
                E = zeros(D, length(ch))
                E[dynamic, :] .= columns(ch)
                Y = carried(constrained!(E, rec), rec, reuse, history, nothing)
                for component in floating
                    share = Y[first(component), :]
                    for node in component
                        Y[node, :] .-= share
                    end
                end
                out[:, ch] .= view(Y, dynamic, :)
            end
            return nothing
        end
        return out
    end
    # the map, written column by column into the one matrix the eigensolve
    # overwrites
    nd = length(dynamic)
    function unit(ch)
        E = zeros(nd, length(ch))
        for (k, i) in enumerate(ch)
            E[i, k] = 1.0
        end
        return E
    end
    M = images!(zeros(nd, nd), unit, rec, history)
    # Every multiplier, and the right and left vectors of the `nev`
    # largest the eigensolve resolves, with their conditions, on the
    # backend (see `mapeigen!`): a multiplier is resolved where it exceeds
    # its error bound, and the others are the algebraic node rates and the
    # modes damped beyond a period's precision.
    eig = mapeigen!(M, method.nev, backend)
    values, selected, VR, VL, column = eig.values, eig.selected, eig.right, eig.left, eig.column
    wi = imag.(values)
    count = length(selected)
    lambda = log.(values[selected]) ./ T
    vector(X, k, c) = eigenvector(X, k, c, wi)
    W = zeros(ComplexF64, D, count)
    for (j, k) in enumerate(selected)
        W[dynamic, j] .= vector(VR, k, column)
    end
    # each multiplier's error bound relative to it, which bounds the error
    # of the pole's real part by itself over the period
    residuals = eps() .* eig.norm1 ./ eig.conditions ./ abs.(values[selected])
    # the error of the steps in each selected rate, from the selected modes
    # carried through the period in twice the steps (see `rateerrors`):
    # the real and imaginary parts are carried, the columns of `VR`, a
    # complex pair's once, a line's history resampled at the finer step and
    # read back at the coarser one's columns
    errors = count == 0 ? zeros(0) : rateerrors(values, selected, VR, VL, column, eig.scale,
        images!(zeros(nd, size(VR, 2)), ch -> VR[:, ch], record(2N), maps(2)), residuals, T; tracking, order = ruleorder)
    # Each mode over the period, the modes in chunks of at most `chunk`
    # directions, their real and imaginary parts, whose samples hold no
    # more than the map does, the same whatever the threads, which share
    # them out: the Fourier coefficients over the steps of a mode's
    # periodic part exp(-lambda t) x(t) place its pole at the harmonic
    # which carries most of it, and give its profile in the harmonics
    # about that one.
    jb = sys.junctionbranches
    Rj = p.matrices.Rbnm[jb, :]
    nj = length(jb)
    modechunk = max(1, min(count, width ÷ 2, D^2 ÷ max(1, 2N*(nn + nj + nz + nl2))))
    groups = chunkranges(count, modechunk)
    # the profiles' harmonics, in the order the other methods give them
    modes = calcfreqsdft((H,)).modes
    m = length(modes)
    voltage, flux = zeros(ComplexF64, m, nn, count), zeros(ComplexF64, m, nj, count)
    poles, edges, stateedges = zeros(ComplexF64, count), zeros(count), zeros(count)
    # The block states in the units the polynomial gives them, balanced
    # as `polebalanced` balances each one's inputs over the frequency scale
    # against its outputs, every output of a pumped block counted, and
    # scaled by `Lscale`, and the waves leaving the line ports in the
    # voltage they carry, `sqrt(Z) a`, over the frequency scale: a mode's
    # states and waves then compare with its node voltages and junction
    # fluxes as in the polynomial, whatever the scale of a block's
    # realization, and a mode is internal to its blocks and lines where
    # those are the eigensolve's roundoff beside its states and waves over
    # the period (see `INTERNALFRACTION`).
    statescale = ones(nz)
    for b in p.blocks
        isempty(b.A) && continue
        outputs = vcat(b.C, (m.C for m in b.modulations)...)
        statescale[b.zbase .+ (1:size(b.A, 1))] .= p.Lscale ./ statebalance(b.B, outputs, scale)
    end
    wavescale = [sqrt(line.Z) for line in p.lines for _ in 1:2]
    internal = zeros(Bool, count)
    realmode = [iszero(wi[k]) for k in selected]
    # a device's samples, copied to the host
    onhost(a) = host ? a : Array(a)
    runchunks(batchchunks(backend, length(groups))) do _, range
        reuse = TransientReuse(rec.system, nothing, nothing, nothing, nothing, nothing)
        for js in view(groups, range)
            nc = length(js)
            # the real and imaginary parts of each mode, side by side
            Wc = zeros(D, 2nc)
            Wc[:, 1:2:end] .= real.(W[:, js])
            Wc[:, 2:2:end] .= imag.(W[:, js])
            rates, fluxes, states, waves = zeros(nn, 2nc, N), zeros(nj, 2nc, N), zeros(nz, 2nc, N), zeros(nl2, 2nc, N)
            function sink(k, x, v, z, a)
                k <= N || return nothing
                rates[:, :, k] .= onhost(view(v, 1:nn, :))
                fluxes[:, :, k] .= Rj*onhost(view(x, 1:nn, :))
                isnothing(z) || (states[:, :, k] .= onhost(z))
                isnothing(a) || (waves[:, :, k] .= wavescale .* onhost(a))
                return nothing
            end
            carried(constrained!(Wc, rec), rec, reuse, history, sink)
            for (c, j) in enumerate(js)
                periodic = cis.(-imag(lambda[j]) .* (0:N-1) .* (T/N)) .* exp.(-real(lambda[j]) .* (0:N-1) .* (T/N))
                coefficients(a) = size(a, 1) == 0 ? zeros(ComplexF64, 0, N) :
                    FFTW.fft((view(a, :, 2c - 1, :) .+ im .* view(a, :, 2c, :)) .* transpose(periodic), 2) ./ N
                V, J, Z = phi0 .* coefficients(rates), phi0 .* coefficients(fluxes), coefficients(states)
                # the node voltages over the frequency scale and the junction
                # fluxes, as the polynomial weighs them, or the blocks' states
                # in its units and the lines' waves for a mode internal to
                # them, whose profile is zero
                content = vec(sum(abs2, V ./ scale; dims = 1)) .+ vec(sum(abs2, J; dims = 1))
                statecontent = vec(sum(abs2, Z .* statescale; dims = 1)) .+ vec(sum(abs2, coefficients(waves) ./ scale; dims = 1))
                internal[j] = sqrt(sum(content)) <= INTERNALFRACTION*sqrt(sum(statecontent))
                # a real multiplier's mode is real, its content the same at
                # a harmonic and at that harmonic's mirror across zero
                # frequency, and its pole is placed at the nonnegative
                # frequency of the two, whatever the roundoff between them
                weight = internal[j] ? statecontent : content
                harmonic(i) = i - 1 <= N ÷ 2 ? i - 1 : i - 1 - N
                candidates = filter(i -> !realmode[j] || imag(lambda[j]) + harmonic(i)*omega >= 0, 1:N)
                nstar = harmonic(candidates[argmax(view(weight, candidates))])
                poles[j] = lambda[j] + im*nstar*omega
                window = [mod(nstar + only(mode), N) + 1 for mode in modes]
                edge = [abs(only(mode)) == H for mode in modes]
                total = sum(content[window])
                if !internal[j] && total > 0
                    voltage[:, :, j] .= transpose(V[:, window]) ./ sqrt(total)
                    flux[:, :, j] .= transpose(J[:, window]) ./ sqrt(total)
                    edges[j] = H > 0 ? sum(content[window[edge]])/total : 0.0
                end
                statetotal = sum(statecontent[window])
                H > 0 && statetotal > 0 && (stateedges[j] = sum(statecontent[window[edge]])/statetotal)
            end
        end
        return nothing
    end
    order = poleorder(poles)
    searches = NamedTuple[(; steps = N, period = T, multipliers = nd - eig.unresolved, tested = eig.tested,
        gauge = length(floating), history = nw)]
    return HBStabilityResult(poles[order], residuals[order], errors[order], edges[order], voltage[:, :, order], flux[:, :, order],
        modes, String.(psc.nodenames[2:end]), copy(jb), zeros(Int, count), searches, true,
        eig.unresolved, 0, omega, scale, stateedges[order], BitVector(internal[order]), zeros(Int, count), method)
end

# The error of the steps in each of the period map's selected rates (see
# `polemonodromy`): the selected modes `selected` of the multipliers
# `values`, their right and left vectors `VR` and `VL` as `mapeigen!` holds
# them at `column`, both of a complex pair's, with the balancing's
# `scale`, and the images `Y2` of `VR`'s columns through the period in
# twice the steps give the finer map projected onto those columns along
# the left ones, `G = (VL' VR) \ VL' Y2`, real, whose eigenvalues are the
# multipliers there of the modes the columns hold, with their mixing among
# them followed. A mode continues into the eigenvector of `G` which
# carries more than `tracking` of it, by the norms of its components along
# the modes, a complex pair's two from its real and imaginary columns, in
# the balanced coordinates the conditions are taken in, and the change of
# its rate, scaled to the coarser map's error at the rule's `order`, with
# the eigensolve's bound on it, `residuals` over the period `T`, estimates
# its error; a mode which no eigenvector carries so, or more than one
# does, has no estimate.
function rateerrors(values::Vector{ComplexF64}, selected::Vector{Int}, VR::Matrix{Float64}, VL::Matrix{Float64},
        column::Vector{Int}, scale::Vector{Float64}, Y2::Matrix{Float64}, residuals::Vector{Float64}, T::Float64;
        tracking::Float64, order::Int)
    wi = imag.(values)
    refined = eigen((VL'*VR) \ (VL'*Y2))
    # the modes the columns hold, a pair's two, and their lengths
    held = sort!(unique!([q for k in selected for q in (wi[k] > 0 ? (k, k + 1) : wi[k] < 0 ? (k - 1, k) : (k,))]))
    lengths = [norm(eigenvector(VR, k, column, wi) ./ scale) for k in held]
    # the component along the mode `k` of a vector `w` over the columns
    along(w, k) = wi[k] == 0 ? abs(w[column[k]]) : wi[k] > 0 ? abs(w[column[k]] - im*w[column[k + 1]])/2 :
        abs(w[column[k - 1]] + im*w[column[k]])/2
    position = zeros(Int, length(values))
    position[selected] .= eachindex(selected)
    errors, claims = fill(Inf, length(selected)), zeros(Int, length(selected))
    richardson = 2.0^order/(2.0^order - 1)
    for c in eachindex(refined.values)
        w = view(refined.vectors, :, c)
        shares = [along(w, k)*l for (k, l) in zip(held, lengths)]
        i = argmax(shares)
        j = position[held[i]]
        j > 0 && shares[i] > tracking*sum(shares) || continue
        claims[j] += 1
        errors[j] = (richardson*abs(log(abs(refined.values[c])/abs(values[selected[j]]))) + residuals[j])/T
    end
    errors[claims .!= 1] .= Inf
    return errors
end

# The history the period map holds, as the line port and the column at
# the step `h` of each coordinate: of the wave leaving each port, the
# columns from the first which a read of its far port at or after the
# start reaches, at the step or at half of it, and at least `points` of
# them, which `historymaps` resamples a mode's history from at half the
# step, on to the start. A column which no read at the step reaches is
# the mode's own history all the same, of the map's image as of its
# eigenvectors, and adds a zero multiplier alone.
function historycolumns(p::TransientProblem, h::Float64, points::Int)
    npre, npre2 = lineprehistory(p, h), lineprehistory(p, h/2)
    stencil = zeros(6)
    rows, cols = Int[], Int[]
    for (l, line) in enumerate(p.lines)
        before, _ = linestencil!(stencil, -line.delay, linestart(0.0, npre, h), h, npre)
        before2, _ = linestencil!(stencil, -line.delay, linestart(0.0, npre2, h/2), h/2, npre2)
        # the finer step's first column read, at or after this column
        start = max(1, min(before + 1, npre - cld(npre2 - before2 - 1, 2), npre - points + 1))
        for e in 1:2, c in start:npre
            push!(rows, 2(l - 1) + e)
            push!(cols, c)
        end
    end
    return rows, cols
end

# The maps between the period map's history coordinates and the history
# a tangent at the step `h/refinement` takes and returns. The
# coordinates are the columns `cols` of the line ports `rows` at the step
# `h`, each port's contiguous and in order, in the units of the rates,
# `sqrt(Z) a/phi0`, the voltage the wave carries. `into` gives the
# `npre` columns of the history before the start at the finer step,
# `(port, column)` flattened: a coordinate at its own column, and the
# Lagrange interpolation of `points` coordinates between them (see
# `lagrangestencil`). `from` reads the coordinates back from a history at
# the end, at their columns, but for a column the finer step's history
# no longer holds, which no read at either step reaches.
function historymaps(p::TransientProblem, rows::Vector{Int}, cols::Vector{Int}, h::Float64, refinement::Int,
        points::Int)
    nl2 = 2length(p.lines)
    npre, npres = lineprehistory(p, h), lineprehistory(p, h/refinement)
    unit = [phi0/sqrt(p.lines[(r + 1) ÷ 2].Z) for r in rows]
    ii, jj, vv = Int[], Int[], Float64[]
    fi, fj, fv = Int[], Int[], Float64[]
    k1 = 0
    while k1 < length(rows)
        # the port's coordinates, the run of its rows
        k0 = k1 + 1
        r = rows[k0]
        k1 = k0
        while k1 < length(rows) && rows[k1 + 1] == r
            k1 += 1
        end
        ks = k0:k1
        # the finer step's column of the port's first coordinate
        c1 = npres - refinement*(npre - cols[first(ks)])
        for c in max(c1, 1):npres
            j0, weights = lagrangestencil((c - c1)/refinement, length(ks), points)
            for (i, w) in enumerate(weights)
                iszero(w) && continue
                push!(ii, r + (c - 1)*nl2); push!(jj, ks[j0 + i]); push!(vv, w*unit[ks[j0 + i]])
            end
        end
        for k in ks
            c = npres - refinement*(npre - cols[k])
            c >= 1 || continue
            push!(fi, k); push!(fj, r + (c - 1)*nl2); push!(fv, 1/unit[k])
        end
    end
    return (; npre = npres, into = sparse(ii, jj, vv, nl2*npres, length(rows)),
        from = sparse(fi, fj, fv, length(rows), nl2*npres))
end

# The Lagrange weights at `x` of `points` consecutive samples of the `n`
# at `0:n-1`, of all where there are fewer: centered on `x` where the
# samples allow, and shifted at the ends rather than narrowed, since a
# resampling, unlike a line's read at every step, need not be a
# contraction. Returns the sample before the first, and the weights.
function lagrangestencil(x, n::Int, points::Int)
    k = min(points, n)
    j0 = clamp(floor(Int, x) - (k ÷ 2 - 1), 0, n - k)
    return j0, [prod((x - j0 - m)/(i - m) for m in 0:k - 1 if m != i; init = 1.0) for i in 0:k - 1]
end

# The multipliers of the period map `M`, which it overwrites, every one,
# and the right and left vectors of the `nev` largest the eigensolve
# resolves, with their conditions, on `backend` (see `mapspectrum!`), or
# on the host where a device's vectors fail their residuals. `M` is
# balanced as LAPACK's geevx balances it, by a diagonal scaling, and a
# multiplier is resolved where it exceeds its error bound, eps times the
# balanced map's 1-norm `norm1` over its condition, the cosine of the
# angle between its left and right vectors: the largest are tested so in
# turn, and the others, whose vectors are not formed, as a well
# conditioned multiplier is, against eps times `norm1`. `unresolved`
# counts the multipliers which fail, and `tested` those tested against
# their own condition. The vectors are held as LAPACK holds
# them (see `eigenvector`), in the map's coordinates: `right` and `left`
# hold the columns of the chosen ones, both of a complex pair's, at
# `column[k]` for the multiplier `k`; `scale` is the balancing's.
function mapeigen!(M::Matrix{Float64}, nev::Int, backend::Backend)
    ilo, ihi, scale = LAPACK.gebal!('S', M)
    norm1 = opnorm(M, 1)
    chosen = mapchoose(mapspectrum!(M, ilo, ihi, backend)..., norm1, nev)
    # a device's vectors short of their residuals, which leaves the map as
    # it was: the host's in their place
    isnothing(chosen) && (chosen = mapchoose(mapspectrum!(M, ilo, ihi, CPU())..., norm1, nev))
    values, selected, column, right, left = chosen.values, chosen.selected, chosen.column, chosen.right, chosen.left
    wi = imag.(values)
    # the chosen ones' columns alone, in the map's coordinates
    keep = sort!(unique!(Int[c for k in selected for c in (wi[k] > 0 ? (column[k], column[k] + 1) :
        wi[k] < 0 ? (column[k] - 1, column[k]) : (column[k],))]))
    at = zeros(Int, size(right, 2))
    at[keep] .= eachindex(keep)
    column .= [c == 0 ? 0 : at[c] for c in column]
    return (; values, norm1, scale, selected, conditions = chosen.conditions, unresolved = chosen.unresolved,
        tested = chosen.tested, column,
        right = right[:, keep] .* scale, left = left[:, keep] ./ scale)
end

# The choice of `mapeigen!` from the multipliers `values` of the balanced
# map and the function `vectors` giving the vectors of some (see
# `mapspectrum!`), in the balanced coordinates; nothing where `vectors`
# gives none.
function mapchoose(values::Vector{ComplexF64}, vectors, norm1::Float64, nev::Int)
    n = length(values)
    wi = imag.(values)
    zerobound = eps()*norm1
    order = sortperm(abs.(values); rev = true)
    unresolved = sum(z -> abs(z) <= zerobound, values; init = 0)
    right, left, column = zeros(n, 0), zeros(n, 0), zeros(Int, n)
    selected, conditions = Int[], Float64[]
    next, tested = 1, 0
    while length(selected) < nev && next <= n && abs(values[order[next]]) > zerobound
        stop, need = next, min(nev - length(selected), n)
        while stop < min(n, next + need - 1) && abs(values[order[stop + 1]]) > zerobound
            stop += 1
        end
        batch = order[next:stop]
        next = stop + 1
        # the vectors not yet formed: a pair's two share theirs
        new = filter(k -> column[k] == 0, batch)
        if !isempty(new)
            formed = vectors(new)
            isnothing(formed) && return nothing
            R, L, at = formed
            for (k, c) in at
                column[k] = size(right, 2) + c
            end
            right, left = hcat(right, R), hcat(left, L)
        end
        tested += length(batch)
        for k in batch
            v, u = eigenvector(right, k, column, wi), eigenvector(left, k, column, wi)
            s = abs(innerproduct(u, v))/(norm(u)*norm(v))
            abs(values[k]) > zerobound/s ? (push!(selected, k); push!(conditions, s)) : (unresolved += 1)
        end
    end
    return (; values, selected, conditions, unresolved, tested, column, right, left)
end

# The multipliers of the balanced map `M`, which it overwrites, and a
# function giving the right and left vectors of chosen ones, as `vectors`
# below does, on the host: the Schur form, as geevx forms it, and the
# vectors of the chosen ones alone, of the Schur form by back
# substitution and carried back through its vectors, where geevx forms
# every one. A device's (CUDA extension) takes every right vector and
# the left ones from their inverse, and gives nothing where those fail.
function mapspectrum!(M::Matrix{Float64}, ilo::Int, ihi::Int, ::Backend)
    n = size(M, 1)
    _, tau = LAPACK.gehrd!(ilo, ihi, M)
    T = triu(M, -1)
    LAPACK.orghr!(ilo, ihi, M, tau)
    _, Z, values = LAPACK.hseqr!('S', 'V', ilo, ihi, T, M)
    wi = imag.(values)
    # The right and left vectors of the multipliers `ks` and of their
    # pairs' partners, as LAPACK holds them, and the column of each one
    # formed, `k => c`.
    function vectors(ks)
        select = zeros(BlasInt, n)
        select[ks] .= 1
        yl, yr = schurvectors!(select, T, wi)
        at, c = Pair{Int,Int}[], 1
        for k in 1:n
            select[k] == 0 && continue
            push!(at, k => c)
            wi[k] == 0 || push!(at, k + 1 => c + 1)
            c += wi[k] == 0 ? 1 : 2
        end
        return Z*yr, Z*yl, at
    end
    return values, vectors
end

# the vector of the eigenvalue `k` from the columns of `X` at `c[k]`, as
# LAPACK holds them: those of a complex pair are the columns of its real
# and imaginary parts, the second's conjugate the first's
eigenvector(X, k, c, wi) = wi[k] > 0 ? complex.(view(X, :, c[k]), view(X, :, c[k + 1])) :
    wi[k] < 0 ? complex.(view(X, :, c[k - 1]), .-view(X, :, c[k])) : complex.(view(X, :, c[k]))

# The left and right vectors of the eigenvalues `select` marks of the
# quasi-triangular Schur form `T`, by back substitution (LAPACK's dtrevc),
# as LAPACK holds them; `select` is left marking a pair's first.
function schurvectors!(select::Vector{BlasInt}, T::Matrix{Float64}, wi::Vector{Float64})
    n = size(T, 1)
    mm = sum(k -> select[k] == 0 ? 0 : wi[k] == 0 ? 1 : 2, 1:n; init = 0)
    VL, VR = zeros(n, mm), zeros(n, mm)
    m, info, work = Ref{BlasInt}(0), Ref{BlasInt}(0), zeros(3n)
    ccall((LinearAlgebra.BLAS.@blasfunc(dtrevc_), LinearAlgebra.LAPACK.liblapack), Cvoid,
        (Ref{UInt8}, Ref{UInt8}, Ptr{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ptr{Float64},
         Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt}, Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Clong, Clong),
        'B', 'S', select, n, T, max(1, n), VL, max(1, n), VR, max(1, n), mm, m, work, info, 1, 1)
    LinearAlgebra.LAPACK.chklapackerror(info[])
    return VL[:, 1:m[]], VR[:, 1:m[]]
end

# the ranges of `1:n` in chunks of `width`, the last one shorter
chunkranges(n::Int, width::Int) = [j0:min(j0 + width - 1, n) for j0 in 1:width:n]

# `u' v`, a loop the first call of the map infers at once where `dot`
# infers its BLAS call
function innerproduct(u, v)
    s = zero(promote_type(eltype(u), eltype(v)))
    for i in eachindex(u, v)
        s += conj(u[i])*v[i]
    end
    return s
end
