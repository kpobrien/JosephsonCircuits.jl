# Fitting a rational scattering block to sampled scattering data: common
# poles by the relaxed pole relocation of Gustavsen and Semlyen, residues
# and feedthrough by least squares at the settled poles, a real state
# space realization, and passivity enforcement by the smallest residue
# perturbation that brings every singular value back under one. Nothing
# here needs a package beyond linear algebra.

# A pole is real when its imaginary part is below this fraction of its
# magnitude. One shared constant rather than a keyword: the basis, the
# residue map, the real pole matrix, the relocation's pairing and the
# realization must all classify a pole the same way, or a conjugate
# pair splits across representations.
const realpoletolerance = 1e-6

# A singular value of the residue basis below this fraction of the
# largest is a direction the samples do not determine, which the residue
# solve leaves out. One shared constant rather than a keyword: the
# relocation, the pruning and the final solve all judge a fit by the
# residues this solve returns, so they must leave out the same directions.
const fitranktolerance = 1e-12

"""
    RationalScattering(block::ScatteringParameters, npoles;
        frequencies = nothing, iterations = 30, passivity = true, atol = 1e-8,
        margin = 1e-6, rounds = 20, pruneslack = 0.05)

A [`RationalScattering`](@ref) block fitted to the scattering data of
`block` at `npoles` common poles by vector fitting: starting poles
spread over the band are relocated by the relaxed iteration of Gustavsen
and Semlyen, for at most `iterations` rounds, until they settle, and the
residues of every entry and the constant term follow by least squares.
The data is sampled at `frequencies` in Hz, by default a tabulated
block's own. Poles the data does not need drift out of the band or
coalesce and are pruned, so `npoles` is a budget rather than the order
returned; too few poles settle on a poor fit that typically cannot be
made passive, which is an error. The result is a real state space
realization, stable by construction, with as many states per pole as its
residue has rank. The fitted block keeps the reference impedances,
grounding and noise model of `block`. A delay is not a rational
function: model a cable as a [`TransmissionLine`](@ref) of its delay in
cascade with a fit of the data with that delay removed.

With `passivity = true` the fit is perturbed wherever its largest
singular value crosses one, by the smallest change of its residues and
constant that brings each crossing to `1 - margin`, over at most
`rounds` rounds, and the result is validated by the same test as any
rational block, to `atol`. `passivity = false` returns the raw fit,
without the enforcement or its repair of an active constant term, but
still rejects a fit that fails validation. The test decides on a lower
bound of the largest singular value, so a returned fit can have a true
norm above one by up to twice the norm search's relative tolerance;
[`passivityassessment`](@ref) reports that uncertainty, and a solve
which cannot tolerate a block active within it should ask for a `margin`
above it.

If `block` states its zero frequency behavior, through the `dcmodel` it
was built with, the fit meets the statement exactly: the value at zero
is linear in the residues once the poles settle, so it is imposed as an
equality, on the residue fit and through the passivity enforcement. The
statement is the caller's to make and is not checked against the data,
which begins above zero and cannot check it. A statement far below the
band may need poles placed there: on the connector data in the tests, a
stated through or open is refused at twenty poles and accepted once a
few poles sit below the band. A statement of unit norm -- a through, an
open, a short -- pins the norm of any fit meeting it at one, so such a
fit is accepted at one to `atol` rather than contracted, which would
move the statement; where the statement and passivity cannot both hold,
the fit is refused.

`pruneslack` bounds how much the fit error may grow over the whole
pruning, as a fraction of the error before any pole was dropped. It
bounds the pruned fit, not the returned block: the feedthrough repair
and the passivity enforcement come after. `margin` must be well below
the block's own dissipation `1 - sigma`, or the enforcement writes over
the loss it should preserve.

The dissipation `I - S S'` is a difference of nearly equal quantities
when the block is nearly lossless, so a fit error `E` appears in it as
roughly `2E`. A block with dissipation of order `1e-5` needs a fit
accurate to well under that before its noise means anything; compare the
fitted dissipation against the data's before trusting the noise of a
fit.
"""
function RationalScattering(block::ScatteringParameters, npoles::Integer; frequencies = nothing,
        iterations::Integer = 30, passivity::Bool = true, atol::Real = 1e-8,
        margin::Real = 1e-6, rounds::Integer = 20, pruneslack::Real = 0.05)
    npoles >= 1 || throw(ArgumentError("fit at least one pole."))
    (isfinite(margin) && 0 <= margin < 1) || throw(ArgumentError(
        "margin must be finite and in [0, 1): it is how far below one a singular value is put."))
    (isfinite(atol) && atol >= 0) || throw(ArgumentError("atol must be finite and nonnegative."))
    iterations >= 1 || throw(ArgumentError("give at least one relocation iteration."))
    rounds >= 1 || throw(ArgumentError("give at least one round of passivity enforcement."))
    (isfinite(pruneslack) && pruneslack >= 0) || throw(ArgumentError("pruneslack must be finite and nonnegative."))
    fs, S = samplescattering(block, frequencies)
    return fitsampled(block, S, fs, Int(npoles); iterations = Int(iterations),
        passivity = passivity, atol = atol, margin = Float64(margin),
        rounds = Int(rounds), pruneslack = Float64(pruneslack))
end

# The sample frequencies a fit accepts, checked once for both fitters: a
# sample at zero is data, but at least one positive frequency is needed
# to set the frequency scale, and duplicate frequencies are refused
# because they make the divided differences of the Loewner pencil
# singular.
function checkfrequencies(fs::AbstractVector)
    all(f -> isfinite(f) && f >= 0, fs) || throw(ArgumentError(
        "the sample frequencies must be finite and nonnegative."))
    issorted(fs) || throw(ArgumentError("the sample frequencies must be increasing."))
    any(k -> fs[k] == fs[k + 1], 1:length(fs) - 1) && throw(ArgumentError(
        "the sample frequencies must be strictly increasing; two samples share a frequency."))
    any(>(0), fs) || throw(ArgumentError(
        "the fit needs at least one positive frequency; every sample is at zero."))
    return fs
end

# The frequencies in Hz the fit is to match the block at, and the block
# sampled there, split out because the order search fits the same
# samples many times and a computed provider is not cheap to evaluate.
function samplescattering(block::ScatteringParameters, frequencies)
    # A tabulated block is evaluated at the angular frequencies it
    # stores: dividing by 2 pi and multiplying back is off by an ulp,
    # which puts the first sample outside the table's own range. The
    # frequencies in Hz are still returned, since that is what the fit
    # reports in, but nothing is evaluated at them.
    ws = if isnothing(frequencies)
        block.provider isa TabulatedMatrixProvider || throw(ArgumentError(
            "give the frequencies in Hz to sample the block at; only a tabulated block has its own."))
        copy(block.provider.frequencies)
    else
        2pi .* Float64.(collect(frequencies))
    end
    fs = ws ./ (2pi)
    checkfrequencies(fs)
    S = zeros(ComplexF64, block.nports, block.nports, length(ws))
    evaluatescattering!(S, block, ws)
    return fs, S
end

# One fit of already sampled data at a fixed order.
function fitsampled(block::ScatteringParameters, S, fs, npoles::Int; iterations::Int,
        passivity::Bool, atol::Real, margin::Float64, rounds::Int, pruneslack::Float64)
    # a block which states what it does at zero frequency has the fit meet
    # it exactly; one which states nothing leaves it to the extrapolation
    dc = block.dcmodel isa ScatteringLimit ? nothing :
        dcscatteringmatrix(block.dcmodel, block.nports)
    # the constant term is repaired only for a caller who asked for the
    # fit to be made passive, and to that caller's tolerances
    poles, residues, D = vectorfit(S, 2pi .* fs, npoles, iterations;
        pruneslack = pruneslack, dc = dc,
        constanttol = passivity ? atol : nothing, constantmargin = margin)
    A, B, C = realization(poles, residues, block.nports)
    if passivity
        A, B, C, D = enforcepassivity(A, B, C, D, 2pi .* fs; atol = atol,
            margin = margin, rounds = rounds, dc = dc)
    end
    return RationalScattering(A, B, C, D; zref = block.zref,
        grounded = block.grounded, noise = block.noise, atol = atol)
end

# The number of poles the samples can determine, from the numerical
# rank of the Loewner pencil of Mayo and Antoulas: the samples split
# alternately into a left and a right set with the interpolation
# directions cycling the ports, and the rank of `[L Ls]` is the McMillan
# degree of the system the data describes. One small singular value
# decomposition, available before any fitting. The degree bounds the
# pole count from above: a pole takes as many states as its residue has
# rank, so dividing the degree by the port count would undercount, but
# every pole takes at least one state.
function supporteddegree(S::AbstractArray{<:Complex,3}, ws::AbstractVector,
        noisefloor::Float64; maxsamples::Int = 400)
    n, K = size(S, 1), length(ws)
    # The two halves of the pencil carry different units, a divided
    # difference of the response being an inverse frequency and a
    # shifted one dimensionless, so a rank threshold on the stacked pair
    # counts differently as the unit of frequency changes. Normalizing
    # to the geometric centre of the band makes the count unit
    # independent.
    lo = findfirst(>(0), ws)
    isnothing(lo) && return typemax(Int)
    wref = sqrt(ws[lo]*maximum(ws))
    xs = ws ./ wref
    # the pencil is quadratic in the samples, so dense data is thinned:
    # the rank is a property of the system, not of the sampling density
    idx = K <= maxsamples ? collect(1:K) : unique(round.(Int, range(1, K; length = maxsamples)))
    length(idx) >= 4 || return typemax(Int)
    right, left = idx[1:2:end], idx[2:2:end]
    nr, nl = length(right), length(left)
    (nr >= 1 && nl >= 1) || return typemax(Int)
    L, Ls = zeros(ComplexF64, nl, nr), zeros(ComplexF64, nl, nr)
    for i in 1:nl, j in 1:nr
        mu, lam = im*xs[left[i]], im*xs[right[j]]
        di, dj = (i - 1) % n + 1, (j - 1) % n + 1
        a = S[di, dj, left[i]]
        b = S[di, dj, right[j]]
        L[i, j] = (a - b)/(mu - lam)
        Ls[i, j] = (mu*a - lam*b)/(mu - lam)
    end
    sv = svdvals(hcat(L, Ls))
    (isempty(sv) || !(sv[1] > 0)) && return typemax(Int)
    return max(1, count(>(noisefloor*sv[1]), sv))
end

# How far a fit is from its data, as a fraction of the largest
# response: the worst sample, not the mean, since a fit which is
# excellent almost everywhere and wrong at one resonance is not a model
# of the block.
function relativefiterror(fit, S, fs)
    P = fit.provider
    nz = size(P.A, 1)
    worst, scale = 0.0, 0.0
    # one factorization for the whole sweep, not one solve per sample
    rf = resolventfactors(P.A, P.B)
    CZ = P.C*rf.Z
    for (k, f) in enumerate(fs)
        F = transferat(rf, CZ, P.D, 2pi*f, nz)
        worst = max(worst, opnorm(F .- view(S, :, :, k)))
        scale = max(scale, opnorm(view(S, :, :, k)))
    end
    return worst/max(scale, floatmin(Float64))
end

"""
    RationalScattering(block::ScatteringParameters; tol, minpoles = 4,
        maxpoles = nothing, noisefloor = 1e-12, frequencies = nothing,
        iterations = 30, passivity = true, atol = 1e-8, margin = 1e-6,
        rounds = 20, pruneslack = 0.05)

A [`RationalScattering`](@ref) block fitted to the scattering data of
`block` at the fewest poles that meet `tol`, the largest allowed error
over the samples as a fraction of the largest response, measured on the
block that is returned, after any passivity enforcement. The order
matters beyond tidiness because the states of the fit are the states a
transient solve steps.

The search scans one order at a time from `minpoles`, because more
poles do not always fit better: past the order the data supports, the
pole relocation is decided by directions the samples do not determine,
and the fit is then not merely inaccurate but often cannot be made
passive at all. The orders meeting a tolerance are therefore a window
rather than a tail, so there is nothing to bisect on, and an order that
fails to fit counts as one that missed the tolerance. `minpoles` is the
way not to pay for orders a block is known not to need.

Without `maxpoles` the search budgets itself by the degree the samples
determine, the numerical rank of their Loewner pencil with `noisefloor`
as the rank threshold, as a fraction of the largest singular value; it
scans to twice that degree or to as many poles as samples, whichever is
fewer, and stops sooner if four consecutive orders produce no fit at
all, which is the model class running out. The degree is an estimate,
not a bound: it moves with `noisefloor`, which measured data with a
real noise floor wants larger than the default, which suits data good
to nearly full precision; the pencil is built along cycling coordinate
directions from at most four hundred samples, so dynamics weak or
narrow in the directions it does not probe can be missed; and a
constant term contributes to it. `maxpoles` overrides all of this and
is a hard ceiling.

If no order meets `tol`, the closest fit found, its order, the degree
the samples determine, and why any orders failed to fit are reported as
an error, since a block quietly less accurate than asked for is worse
than none: loosen `tol`, raise `maxpoles`, sample the block more
finely, or fit a narrower band.

See the `npoles` method for the meaning of the remaining arguments, and
for the warning about fitting a block with little loss: a tolerance
which looks tight against `S` may still be far too loose against the
dissipation `I - S S'`, where the error appears roughly doubled.
"""
function RationalScattering(block::ScatteringParameters; tol::Real,
        minpoles::Integer = 4, maxpoles = nothing, noisefloor::Real = 1e-12,
        frequencies = nothing, iterations::Integer = 30, passivity::Bool = true,
        atol::Real = 1e-8, margin::Real = 1e-6, rounds::Integer = 20,
        pruneslack::Real = 0.05)
    (isfinite(tol) && tol > 0) || throw(ArgumentError("tol must be finite and positive."))
    minpoles >= 1 || throw(ArgumentError("give minpoles >= 1."))
    (isfinite(noisefloor) && noisefloor > 0) || throw(ArgumentError("noisefloor must be finite and positive."))
    (isfinite(margin) && 0 <= margin < 1) || throw(ArgumentError(
        "margin must be finite and in [0, 1): it is how far below one a singular value is put."))
    (isfinite(atol) && atol >= 0) || throw(ArgumentError("atol must be finite and nonnegative."))
    iterations >= 1 || throw(ArgumentError("give at least one relocation iteration."))
    rounds >= 1 || throw(ArgumentError("give at least one round of passivity enforcement."))
    (isfinite(pruneslack) && pruneslack >= 0) || throw(ArgumentError("pruneslack must be finite and nonnegative."))
    fs, S = samplescattering(block, frequencies)
    # fitting N poles needs N + 1 samples, so a scan from minpoles needs
    # at least that many; refusing here names the samples rather than
    # reporting a scan of no orders
    length(fs) >= Int(minpoles) + 1 || throw(ArgumentError(
        lazy"fitting $(minpoles) poles needs at least $(Int(minpoles) + 1) sample frequencies; there are $(length(fs))."))
    supported = supporteddegree(S, 2pi .* fs, Float64(noisefloor))
    # the ceiling on the search: the degree the samples determine, since
    # a pole the fit keeps takes at least one state and poles beyond the
    # degree are decided by directions the samples do not carry
    ceiling = isnothing(maxpoles) ?
        (supported == typemax(Int) ? 64 : max(supported, Int(minpoles))) : Int(maxpoles)
    minpoles <= ceiling || throw(ArgumentError(lazy"minpoles is $(minpoles) but the search ceiling is $(ceiling); give a smaller minpoles or a larger maxpoles."))
    besterror, bestpoles = Inf, 0
    # why the failed orders failed, keyed by message with the numbers
    # stripped, so the search reports a reason rather than only that no
    # order met the tolerance -- often every order fails the same way
    failures = Dict{String,Vector{Int}}()
    attempt = np -> begin
        fit, err = try
            candidate = fitsampled(block, S, fs, np; iterations = Int(iterations),
                passivity = passivity, atol = atol, margin = Float64(margin),
                rounds = Int(rounds), pruneslack = Float64(pruneslack))
            (candidate, relativefiterror(candidate, S, fs))
        catch e
            e isa ArgumentError || rethrow()
            # the message without the orders and numbers in it, so that
            # the same failure at twenty orders is one line and not twenty
            key = replace(sprint(showerror, e), r"[-+]?[0-9][0-9.e+-]*" => "N")
            push!(get!(failures, key, Int[]), np)
            (nothing, Inf)
        end
        err < besterror && ((besterror, bestpoles) = (err, np))
        return (fit, err)
    end
    # One order at a time from `minpoles` upward, returning the first
    # that meets the tolerance, which is therefore the fewest poles that
    # do. There is nothing to bisect on: more poles do not always fit
    # better, since past the order the data supports the relocation is
    # decided by directions the samples do not determine, so the orders
    # meeting a tolerance are a window and not a tail -- on the
    # connector data in the tests, the error at thirteen poles is
    # 1.2e-3, at fourteen 4.5e-4 and at fifteen 1.3e-3. An order stepped
    # over is an order about which nothing is known, so every order
    # below the answer has to be fitted, and `minpoles` is the way not
    # to pay for orders a block is known not to need.
    #
    # Without `maxpoles` the degree estimate is a budget and not a wall:
    # it can be short of what the block needs, so the scan continues to
    # twice the estimate, and never past the sample count, beyond which
    # a fit interpolates noise. A caller's `maxpoles` is a wall.
    stop = isnothing(maxpoles) ? min(2*ceiling, length(fs)) : ceiling
    np, barren = Int(minpoles), 0
    while np <= stop
        fit, err = attempt(np)
        err <= tol && return fit
        # A run of orders producing no fit at all is the model class
        # running out, and that does not reverse: past the degree the
        # samples determine, the relocation is decided by directions
        # they do not carry and the result cannot be made passive, at a
        # cost that rises with the order. Stopping needs a run rather
        # than one refusal, which can be a numerical accident: on a
        # delay estimated at four poles, seven cannot be made passive
        # and eight fits to 8e-13. (A test on the error improving is the
        # wrong shape here: the orders meeting a tolerance are a window,
        # so orders of no improvement establish nothing about the next.)
        barren = isfinite(err) ? 0 : barren + 1
        barren >= 4 && break
        np += 1
    end
    reached = min(np, stop)
    why = isempty(failures) ? "" :
        " Some orders could not be fitted at all: " *
        join(("$(length(v)) of them ($(first(v)) to $(last(v)) poles) with \"$(k)\""
              for (k, v) in sort(collect(failures); by = x -> -length(x[2]))), "; ") * "."
    throw(ArgumentError(lazy"no fit between $(minpoles) and $(reached) poles met a tolerance of $(tol): every order in that range was fitted, and the closest was $(besterror) at $(bestpoles) poles, while the samples determine a degree of about $(supported).$(why) Loosen tol, raise maxpoles, sample the block more finely, or fit a narrower band."))
end

# The relaxed vector fit: with the poles `a`, the unknowns of every entry
# are its residues and constant, and the shared unknowns the residues
# and the constant of the weight `sigma(s) = d + sum r_p/(s - a_p)` whose
# zeros are the next poles, with the relaxed normalization that the real
# part of `sigma` averages to one over the samples so that the trivial
# solution is excluded; complex poles come in conjugate pairs and enter
# through the real basis, so every unknown is real. Returns the poles,
# the residue matrices `(n, n, npoles)` and the constant `(n, n)`.
function vectorfit(S::AbstractArray{<:Complex,3}, ws::AbstractVector, npoles::Int, iterations::Int;
        pruneslack::Real = 0.05, dc::Union{Nothing,AbstractMatrix} = nothing,
        constanttol::Union{Nothing,Real} = nothing, constantmargin::Real = 1e-6,
        startdamping::Real = 0.01)
    # a real rational function is real at zero frequency, so a complex
    # statement there cannot be met by any fit
    isnothing(dc) || maximum(abs, imag.(dc)) == 0 || throw(ArgumentError(
        "the value stated at zero frequency must be real: a real rational function is real at zero."))
    # The fit is done at frequencies of order one, about the geometric
    # centre of the band, and the poles and residues are scaled back at
    # the end: in rad/s the basis columns and the constant term's column
    # differ by many orders and the value at zero frequency cannot be
    # written down accurately at all, while the scaled fit is the same
    # function, a pole at `p` with residue `r` standing for a pole at
    # `p*wref` with residue `r*wref`. A band which begins at zero has
    # its lowest positive frequency stand in for the lowest, a sample at
    # zero being legitimate data.
    lo = findfirst(>(0), ws)
    isnothing(lo) && throw(ArgumentError("the fit needs at least one positive frequency; every sample is at zero or below."))
    wref = sqrt(ws[lo]*maximum(ws))
    xs = ws ./ wref
    n, K = size(S, 1), length(xs)
    # the band the starting poles are spread over is the positive part of
    # it: a pole started at zero is a pole on top of a sample there, where
    # the basis is infinite
    wmin, wmax = xs[lo], maximum(xs)
    # the starting poles: complex pairs along the band, a real one if odd
    poles = ComplexF64[]
    npairs = npoles ÷ 2
    spread = npairs == 1 ? [(wmin + wmax)/2] : collect(range(wmin, wmax; length = npairs))
    for w in spread
        push!(poles, complex(-startdamping*w, w))
        push!(poles, complex(-startdamping*w, -w))
    end
    isodd(npoles) && push!(poles, complex(-(wmin + wmax)/2, 0.0))
    poles = converge(S, xs, poles, iterations; dc = dc)
    # a relocation which has lost the data leaves poles which are not
    # finite; fail here, naming the fit, rather than in a factorization
    # built from them
    all(isfinite, poles) || throw(ArgumentError(lazy"the pole relocation diverged at $(npoles) poles: the relocated poles are not finite. Fit with fewer poles, or over a narrower band."))
    poles = prunepoles(S, xs, poles, iterations, pruneslack; dc = dc)
    # the poles are settled from the data before the value at zero is
    # stated, so a condition outside the band never moves them
    residues, D = fitresidues(S, xs, poles; dc = dc)
    (all(isfinite, residues) && all(isfinite, D)) || throw(ArgumentError(lazy"the residues at $(length(poles)) poles are not finite: the least squares of the fit is singular. Fit with fewer poles, or over a narrower band."))
    # A constant term above one is unreachable by the enforcement,
    # which perturbs over a band and cannot reach infinite frequency, so
    # it is brought under here, where the residues can be refitted
    # around it -- and only when the caller asked for passivity: samples
    # of `1.5 - 1.2/(s + 1)` have an exact one pole fit whose constant
    # is 1.5, and a raw fit must return it.
    if !isnothing(constanttol)
        D, clamped = passiveconstant(D; tol = Float64(constanttol), margin = Float64(constantmargin))
        if clamped
            residues, D = fitresidues(S, xs, poles; constant = D, dc = dc)
            all(isfinite, residues) || throw(ArgumentError(lazy"the residues at $(length(poles)) poles are not finite once the constant term is made passive: the least squares of the fit is singular. Fit with fewer poles, or over a narrower band."))
        end
    end
    return poles .* wref, residues .* wref, D
end

# the relocation iterated until the poles settle to a part in a billion,
# or until the fit at them is at the roundoff of the least squares, where
# the weight's columns are combinations of the entries' and the
# relocation would be arbitrary
# the scale an error is measured against, in the norm it is measured in
roundoff(S; scale::Real = 1e-12) = scale*maximum(k -> opnorm(view(S, :, :, k)), axes(S, 3))
# how far above one a fit may be and still be scaled under it rather than
# refused: a block scaled from further away than this is mostly loss
const passivityscalelimit = 1e-2
function converge(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64}, iterations::Int;
        stallpatience::Int = 5, settletol::Real = 1e-9,
        dc::Union{Nothing,AbstractMatrix} = nothing)
    # a pole is compared against the band, and a band which begins at zero
    # has its lowest positive frequency stand in
    lo = findfirst(>(0), ws)
    wfloor = isnothing(lo) ? 1.0 : ws[lo]
    order(x) = (real(x), imag(x))
    exact = roundoff(S)
    # The iterate kept is the one measured to fit best, not the one
    # whose poles last moved least: a relocation which is circling can
    # pass its best fit on an iteration where the poles are moving
    # quickly. Every iterate is measured, at the cost of one residue
    # solve each.
    best, besterror = copy(poles), Inf
    bestchange, stalled = Inf, 0
    for _ in 1:iterations
        err = fiterror(S, ws, poles; dc = dc)
        err < besterror && ((best, besterror) = (copy(poles), err))
        # a fit which is already exact settles its poles too
        err <= exact && return best
        newpoles = relocate(S, ws, poles)
        # sorted against sorted, so that a displacement is divided by the
        # size of the pole it belongs to and not of whichever pole held
        # that position before the sort
        previous = sort(poles; by = order)
        change = length(newpoles) == length(poles) ?
            maximum(abs.(sort(newpoles; by = order) .- previous) ./ max.(abs.(previous), wfloor)) : Inf
        poles = newpoles
        if change < bestchange
            bestchange, stalled = change, 0
        else
            stalled += 1
            # the poles are moving as much as they were several steps ago,
            # so they are circling rather than settling
            stalled >= stallpatience && break
        end
        change < settletol && break
    end
    # whatever ended it, the iterate that is returned is the one measured
    # to fit best, which the last one is only if it was measured
    err = fiterror(S, ws, poles; dc = dc)
    return err < besterror ? poles : best
end

# the largest deviation of the fit at the given poles from the data
# over the samples, in the spectral norm, with the residues refitted:
# the measure the relocation and the pruning are judged on
function fiterror(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64};
        dc::Union{Nothing,AbstractMatrix} = nothing)
    residues, D = fitresidues(S, ws, poles; dc = dc)
    n = size(D, 1)
    # accumulated in place, since a sum over a generator of matrices
    # would allocate one per pole per sample, and this runs once per
    # relocation iteration and once per pruning candidate
    F = Matrix{ComplexF64}(undef, n, n)
    err = 0.0
    @inbounds for (k, w) in enumerate(ws)
        for j in 1:n, i in 1:n
            F[i, j] = D[i, j] - S[i, j, k]
        end
        for p in eachindex(poles)
            c = 1/(im*w - poles[p])
            for j in 1:n, i in 1:n
                F[i, j] += residues[i, j, p]*c
            end
        end
        err = max(err, opnorm(F))
    end
    return err
end
# The poles the fit does not need, once the relocation has settled: a
# needless pole drifts out of the band, above it fitting a constant the
# constant term holds, or below it where the data cannot see it,
# carrying next to nothing; or it settles in a cluster with the pole it
# copies, split by less than the data can resolve, the members carrying
# large residues of opposite sign. Needless poles are dropped, least
# contribution first, clusters are replaced by their mean and relocated
# again, and each simplification is kept only if the refitted fit stays
# within `pruneslack` of the error before any pole was dropped -- a
# budget for the whole pruning, not for each deletion.
function prunepoles(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64},
        iterations::Int, pruneslack::Real = 0.05; dc::Union{Nothing,AbstractMatrix} = nothing)
    # the budget is measured from before any pole was dropped
    baseline = fiterror(S, ws, poles; dc = dc)
    limit = (1 + pruneslack)*baseline + roundoff(S)
    acceptable(err, _) = err <= limit
    before = baseline
    while true
        n0 = length(poles)
        poles, before = dropneedless(S, ws, poles, before, acceptable; dc = dc)
        poles, before = mergeclusters(S, ws, poles, iterations, before, acceptable; dc = dc)
        length(poles) == n0 && break
    end
    return poles
end
# a pole and its conjugate are dropped when the fit without them is as
# close, the pole contributing least to the fit tried first
function dropneedless(S, ws, poles, before, acceptable; dc = nothing)
    while true
        residues, _ = fitresidues(S, ws, poles; dc = dc)
        # the norm of the residue does not vary over the samples, so it
        # is taken once per pole rather than once per pole and sample
        contribution(p) = opnorm(view(residues, :, :, p))*
            maximum(1/abs(im*w - poles[p]) for w in ws)
        candidates = sort(filter(p -> imag(poles[p]) >= 0, eachindex(poles)); by = contribution)
        dropped = false
        for p in candidates
            keep = [q for q in eachindex(poles) if q != p && !(imag(poles[p]) != 0 && poles[q] == conj(poles[p]))]
            trial = poles[keep]
            err = fiterror(S, ws, trial; dc = dc)
            if acceptable(err, before)
                # the last pole is needless only if a constant reproduces
                # the data, which is then no rational block
                isempty(trial) && throw(ArgumentError("no pole is needed: a constant reproduces the data, so give it as a ScatteringParameters matrix."))
                poles, before, dropped = trial, err, true
                break
            end
        end
        dropped || return poles, before
    end
end
function mergeclusters(S, ws, poles, iterations, before, acceptable; dc = nothing)
    # the resolution of the data at a pole: the spacing of the samples
    # around its frequency, the first spacing below the band
    resolution(a) = begin
        k = searchsortedfirst(ws, abs(imag(a)))
        k <= 1 ? ws[2] - ws[1] : k > length(ws) ? ws[end] - ws[end - 1] : ws[k] - ws[k - 1]
    end
    # the representatives, real or in the upper half plane, a pole within
    # the resolution of the real axis being real
    reps = ComplexF64[]
    for a in poles
        imag(a) < -resolution(a) && continue
        push!(reps, abs(imag(a)) <= resolution(a) ? complex(real(a), 0.0) : a)
    end
    # the clusters, by the nearest representative within the resolution
    cluster = collect(1:length(reps))
    root(k) = (while cluster[k] != k; k = cluster[k]; end; k)
    for k in eachindex(reps), l in 1:k - 1
        abs(reps[k] - reps[l]) <= resolution(reps[k]) && (cluster[root(k)] = root(l))
    end
    merged = ComplexF64[]
    for r in unique(root.(eachindex(reps)))
        members = reps[root.(eachindex(reps)) .== r]
        a = sum(members)/length(members)
        if imag(a) == 0
            push!(merged, a)
        else
            push!(merged, a)
            push!(merged, conj(a))
        end
    end
    length(merged) == length(poles) && return poles, before
    merged = converge(S, ws, merged, iterations; dc = dc)
    err = fiterror(S, ws, merged; dc = dc)
    return acceptable(err, before) ? (merged, err) : (poles, before)
end

# the real basis of a pole set: columns of the matrix `Phi(s)` such that a
# real coefficient vector `c` gives `sum_p r_p/(s - a_p)` with conjugate
# residues on conjugate poles; and the map back from coefficients to the
# complex residues
function realbasis(ws::AbstractVector, poles::Vector{ComplexF64})
    K, N = length(ws), length(poles)
    Phi = zeros(ComplexF64, K, N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= realpoletolerance*abs(a)
            Phi[:, p] .= 1 ./ (im .* ws .- real(a))
            p += 1
        else
            Phi[:, p] .= 1 ./ (im .* ws .- a) .+ 1 ./ (im .* ws .- conj(a))
            Phi[:, p + 1] .= im ./ (im .* ws .- a) .- im ./ (im .* ws .- conj(a))
            p += 2
        end
    end
    return Phi
end
function complexresidues(c::AbstractVector, poles::Vector{ComplexF64})
    N = length(poles)
    r = zeros(ComplexF64, N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= realpoletolerance*abs(a)
            r[p] = c[p]
            p += 1
        else
            r[p] = complex(c[p], c[p + 1])
            r[p + 1] = conj(r[p])
            p += 2
        end
    end
    return r
end

# one relocation of the poles: the least squares of every entry with the
# shared weight, the weight's residues, and the zeros of the weight
function relocate(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64};
        weightfloor::Real = 1e-8, conjtol::Real = 1e-8)
    n, K, N = size(S, 1), length(ws), length(poles)
    Phi = realbasis(ws, poles)
    ne = n*n
    # The unknowns are every entry's residues and constant plus the
    # shared weight's residues and constant, all real; the equations,
    # every entry's samples with a zero right hand side plus the
    # relaxation. An entry's own unknowns enter only its own rows, so
    # they are eliminated entry by entry: the QR of the entry's block
    # `[A_e B_e]` leaves the trailing block of its triangle as the
    # entry's contribution to the weight's system, and the weight is
    # solved from all of them and the relaxation at once, so nothing the
    # size of every entry's every sample is ever formed (Deschrijver's
    # fast fit).
    K >= N + 1 || throw(ArgumentError(lazy"fitting $(N) poles needs at least $(N + 1) sample frequencies."))
    Ae = [real.(Phi) ones(K); imag.(Phi) zeros(K)]
    # `A_e` is the same matrix for every entry, so its reflectors are
    # computed once and applied to each entry's `B_e`: with
    # `A_e = Q [R1; 0]`, the triangle of `[A_e B_e]` has as its trailing
    # block the triangle of the rows of `Q' B_e` below `R1`, which is a
    # factorization of `N + 1` columns per entry instead of `2N + 2`.
    F = qr(Ae)
    # the weight's columns scaled to their size over every entry
    scale = zeros(N + 1)
    for e in 1:ne
        i, j = (e - 1) % n + 1, (e - 1) ÷ n + 1
        Se = view(S, i, j, :)
        for p in 1:N
            scale[p] += sum(abs2, Se .* view(Phi, :, p))
        end
        scale[N + 1] += sum(abs2, Se)
    end
    scale = max.(sqrt.(scale), floatmin(Float64))
    reduced = zeros(ne*(N + 1) + 1, N + 1)
    rhs = zeros(ne*(N + 1) + 1)
    Be = zeros(2K, N + 1)
    for e in 1:ne
        i, j = (e - 1) % n + 1, (e - 1) ÷ n + 1
        Se = view(S, i, j, :)
        for p in 1:N
            v = Se .* view(Phi, :, p)
            Be[1:K, p] .= -real.(v) ./ scale[p]
            Be[K + 1:2K, p] .= -imag.(v) ./ scale[p]
        end
        Be[1:K, N + 1] .= -real.(Se) ./ scale[N + 1]
        Be[K + 1:2K, N + 1] .= -imag.(Se) ./ scale[N + 1]
        lmul!(F.Q', Be)
        R = qr(Be[N + 2:2K, :]).R
        reduced[(e - 1)*(N + 1) + 1:e*(N + 1), :] .= R
    end
    # the relaxation: the real part of the weight averages to one; scaled
    # to the size of the data so that it weighs as one sample of it
    weight = sqrt(sum(abs2, S)/K)
    for p in 1:N
        reduced[end, p] = weight*sum(real, view(Phi, :, p))/K/scale[p]
    end
    reduced[end, N + 1] = weight/scale[N + 1]
    rhs[end] = weight*1.0
    x = (reduced \ rhs) ./ scale
    csigma = x[1:N]
    dsigma = x[N + 1]
    # a weight whose constant vanished has its zeros at its poles: it is
    # held away from zero, as Gustavsen does
    abs(dsigma) < weightfloor && (dsigma = copysign(weightfloor, dsigma))
    # the zeros of the weight: the eigenvalues of A - b r'/d in the real
    # form of the pole set
    Ar, br = realpolematrix(poles)
    cr = realresiduerow(csigma, poles)
    zeros_ = eigvals(Ar .- br*transpose(cr) ./ dsigma)
    newpoles = ComplexF64[]
    for z in zeros_
        real(z) > 0 && (z = complex(-real(z), imag(z)))
        push!(newpoles, z)
    end
    # conjugate pairs paired, real poles real
    sort!(newpoles; by = x -> (real(x), imag(x)))
    out = ComplexF64[]
    used = falses(length(newpoles))
    for (k, z) in enumerate(newpoles)
        used[k] && continue
        if abs(imag(z)) <= realpoletolerance*abs(z)
            push!(out, complex(real(z), 0.0))
            used[k] = true
        else
            partner = findfirst(l -> !used[l] && l != k && abs(newpoles[l] - conj(z)) < conjtol*abs(z), eachindex(newpoles))
            if isnothing(partner)
                push!(out, complex(real(z), 0.0))
                used[k] = true
            else
                push!(out, complex(real(z), abs(imag(z))))
                push!(out, complex(real(z), -abs(imag(z))))
                used[k] = true
                used[partner] = true
            end
        end
    end
    return out
end

# the real block diagonal matrix of a pole set, its input column of ones
# in the real form, and the coefficient row of the weight's residues in
# the same form, so that `A - b c'` has the weight's zeros as eigenvalues
function realpolematrix(poles::Vector{ComplexF64})
    N = length(poles)
    A = zeros(N, N)
    b = zeros(N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= realpoletolerance*abs(a)
            A[p, p] = real(a)
            b[p] = 1.0
            p += 1
        else
            A[p, p] = real(a); A[p, p + 1] = imag(a)
            A[p + 1, p] = -imag(a); A[p + 1, p + 1] = real(a)
            b[p] = 2.0; b[p + 1] = 0.0
            p += 2
        end
    end
    return A, b
end
function realresiduerow(c::AbstractVector, poles::Vector{ComplexF64})
    N = length(poles)
    row = zeros(N)
    p = 1
    while p <= N
        a = poles[p]
        if imag(a) == 0 || abs(imag(a)) <= realpoletolerance*abs(a)
            row[p] = c[p]
            p += 1
        else
            row[p] = c[p]; row[p + 1] = c[p + 1]
            p += 2
        end
    end
    return row
end

# The least squares of the residue solve, factored once and applied to
# the right hand side of every entry. Directions of the basis whose
# singular value is below `fitranktolerance` of the largest are left out
# of the solution: poles which have coalesced give columns the samples
# cannot tell apart, and a plain least squares answers them with large
# cancelling residues which corrupt the fit error and the value at zero
# frequency. The minimum norm solution leaves those directions at zero.
struct FitLeastSquares
    U::Matrix{Float64}
    s::Vector{Float64}
    V::Matrix{Float64}
    scale::Vector{Float64}
end
function FitLeastSquares(A::AbstractMatrix)
    # the rank is judged with every column at unit norm: a column
    # carries the scale of its basis function, and the smallest singular
    # value of the unscaled basis measures that scale rather than what
    # the samples determine
    scale = [norm(view(A, :, j)) for j in axes(A, 2)]
    for j in eachindex(scale)
        scale[j] > 0 || (scale[j] = 1.0)
    end
    F = svd(A ./ scale')
    kept = isempty(F.S) ? 0 : count(>(fitranktolerance*first(F.S)), F.S)
    return FitLeastSquares(F.U[:, 1:kept], F.S[1:kept], F.V[:, 1:kept], scale)
end
fitleastsquares(F::FitLeastSquares, b::AbstractVector) =
    (F.V*((F.U'*b) ./ F.s)) ./ F.scale

# the residues of every entry and the constant at fixed poles, by least
# squares in the real basis. With `constant` given, that matrix is held
# and only the strictly proper part is fitted, to `S - constant`.
function fitresidues(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64};
        constant::Union{Nothing,AbstractMatrix} = nothing,
        dc::Union{Nothing,AbstractMatrix} = nothing)
    n, K, N = size(S, 1), length(ws), length(poles)
    fixed = !isnothing(constant)
    Phi = realbasis(ws, poles)
    cols = fixed ? N : N + 1
    M = zeros(2K, cols)
    M[1:K, 1:N] .= real.(Phi)
    M[K + 1:2K, 1:N] .= imag.(Phi)
    fixed || (M[1:K, N + 1] .= 1.0)
    residues = zeros(ComplexF64, n, n, N)
    D = fixed ? Matrix{Float64}(constant) : zeros(n, n)
    if isnothing(dc)
        F = FitLeastSquares(M)
        for i in 1:n, j in 1:n
            rhs = fixed ?
                vcat(real.(view(S, i, j, :)) .- D[i, j], imag.(view(S, i, j, :))) :
                vcat(real.(view(S, i, j, :)), imag.(view(S, i, j, :)))
            x = fitleastsquares(F, rhs)
            residues[i, j, :] .= complexresidues(x[1:N], poles)
            fixed || (D[i, j] = x[N + 1])
        end
        return residues, D
    end
    # The value at zero frequency is the same combination of the same
    # basis, read at zero instead of on the axis, and linear in the
    # unknowns once the poles are settled, so stating it is one equality
    # per entry, met exactly: a particular solution, then a least
    # squares over the null space of the equality. A sample near zero
    # would instead ask the fit to take that value at a frequency where
    # it does not, at an excursion far worse than the extrapolation it
    # corrects.
    row = vec(real.(realbasis([0.0], poles)))
    c = fixed ? row : vcat(row, 1.0)
    cc = dot(c, c)
    cc > 0 || throw(ArgumentError("the zero frequency row of the basis vanishes, so the value there cannot be stated."))
    Z = nullspace(reshape(c, 1, :))
    F = FitLeastSquares(M*Z)
    for i in 1:n, j in 1:n
        target = Float64(real(dc[i, j])) - (fixed ? D[i, j] : 0.0)
        x0 = c .* (target/cc)
        rhs = fixed ?
            vcat(real.(view(S, i, j, :)) .- D[i, j], imag.(view(S, i, j, :))) :
            vcat(real.(view(S, i, j, :)), imag.(view(S, i, j, :)))
        x = x0 .+ Z*fitleastsquares(F, rhs .- M*x0)
        residues[i, j, :] .= complexresidues(x[1:N], poles)
        fixed || (D[i, j] = x[N + 1])
    end
    return residues, D
end

"""
    passiveconstant(D; tol = 1e-6, margin = tol)

The constant term `D` with every singular value above `1 + tol` brought
to `1 - margin`, leaving the singular vectors alone. `D` is the model
at infinite frequency, which no perturbation over a band of frequencies
can reach, so an active constant term would leave the passivity
enforcement searching a band which begins at infinity; bringing it
under one before the residues are fitted removes that failure at its
source. `tol` is how far above one a singular value must stand to be
treated as active, and belongs with the tolerance the block is
validated against; `margin` is how far under one it is put, and belongs
with the margin the enforcement aims for. Follows Gustavsen, IEEE
Transactions on Electromagnetic Compatibility 67(3), 2025, section X-A.
"""
function passiveconstant(D::AbstractMatrix; tol::Real = 1e-6, margin::Real = tol)
    F = svd(D)
    # only a constant term genuinely above one is brought down. A block
    # which is lossless at infinite frequency has a unitary `D` by right,
    # a series inductor between ports for one, and moving its singular
    # values would put an error into a fit which was exact.
    any(>(1 + tol), F.S) || return Matrix{Float64}(D), false
    σ = [s > 1 + tol ? 1 - margin : s for s in F.S]
    return F.U*Diagonal(σ)*F.Vt, true
end

# the real state space realization of poles and residue matrices: every
# real pole `a` with residue `R` is `n` states `A = a I, B = I, C = R`, and
# every conjugate pair `a, conj(a)` with residues `R, conj(R)` the `2n`
# states `A = [Re a I  -Im a I; Im a I  Re a I]`, `B = [I; 0]`,
# `C = [2 Re R  -2 Im R]`
function realization(poles::Vector{ComplexF64}, residues::AbstractArray{<:Complex,3}, n::Int;
        ranktol::Real = 1e-12)
    N = length(poles)
    blocks = Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}}[]
    p = 1
    while p <= N
        a = poles[p]
        R = residues[:, :, p]
        # the states of a pole are the rank of its residue, `R = U S V'`
        # split as `C = U sqrt(S)`, `B = sqrt(S) V'`, so that the
        # realization is minimal and balanced; a real pole with a
        # residue of rank `r` takes `r` states, a pair `2r`
        F = svd(R)
        r = count(s -> s > ranktol*max(F.S[1], floatmin(Float64)), F.S)
        Bc = Diagonal(sqrt.(F.S[1:r]))*F.Vt[1:r, :]
        Cc = F.U[:, 1:r]*Diagonal(sqrt.(F.S[1:r]))
        if imag(a) == 0 || abs(imag(a)) <= realpoletolerance*abs(a)
            r > 0 && push!(blocks, (real(a) .* Matrix(1.0I, r, r), real.(Bc), real.(Cc)))
            p += 1
        else
            # the complex state `z' = a z + Bc u`, `y = 2 Re(Cc z)`, as its
            # real and imaginary parts
            α, β = real(a), imag(a)
            r > 0 && push!(blocks, ([α .* Matrix(1.0I, r, r)  -β .* Matrix(1.0I, r, r); β .* Matrix(1.0I, r, r)  α .* Matrix(1.0I, r, r)],
                [real.(Bc); imag.(Bc)], [2 .* real.(Cc)  -2 .* imag.(Cc)]))
            p += 2
        end
    end
    isempty(blocks) && return zeros(0, 0), zeros(0, n), zeros(n, 0)
    nz = sum(size(b[1], 1) for b in blocks)
    A, B, C = zeros(nz, nz), zeros(nz, n), zeros(n, nz)
    offset = 0
    for (Ab, Bb, Cb) in blocks
        k = size(Ab, 1)
        A[offset + 1:offset + k, offset + 1:offset + k] .= Ab
        B[offset + 1:offset + k, :] .= Bb
        C[:, offset + 1:offset + k] .= Cb
        offset += k
    end
    return A, B, C
end

# The passivity enforcement: the bands where the largest singular value
# of the fit stands above one are found by sweeping the fit over
# frequency, the residues and the constant are perturbed by the least
# change, in the norm of the fit over the samples, that brings each
# band's worst point to one less a margin -- a linear constraint on the
# perturbation through the singular vectors -- and the rounds repeat
# until the exact norm test at the end passes, with a bound on their
# number.
#
# The resolvent through a Schur form: `A` and `B` are fixed while the
# enforcement perturbs `C` and `D`, so the factorization is taken once
# and every later evaluation is a triangular solve, `O(nz^2 n)` rather
# than the `O(nz^3)` of a fresh dense solve at each frequency. The
# violation sweep evaluates the fit over a grid in every band and
# dominates the enforcement, so this is where its time goes.
struct ResolventFactors{TZ,TT,TB}
    Z::TZ
    T::TT
    ZtB::TB
end
function resolventfactors(A, B)
    F = schur(complex(Matrix(A)))
    return ResolventFactors(F.Z, F.T, F.Z'*B)
end
# `(i w I - A)^-1 B` through the held Schur factors
resolventat(rf::ResolventFactors, w) =
    rf.Z*(UpperTriangular(im*w*I - rf.T) \ rf.ZtB)

# the fit at one frequency, `D + C (i w I - A)^-1 B`, through the held
# factors; a realization whose poles have run away cannot be evaluated,
# and the failure is named here rather than surfacing from a
# factorization downstream
function transferat(rf::ResolventFactors, CZ, D, w, nstates)
    isfinite(w) || throw(ArgumentError(
        lazy"the fit with $(nstates) states is not passive at infinite frequency: its feedthrough alone has a singular value above one, so no perturbation over the band can make it passive. Fit with fewer poles, or over a narrower band."))
    F = D .+ CZ*(UpperTriangular(im*w*I - rf.T) \ rf.ZtB)
    all(isfinite, F) || throw(ArgumentError(
        lazy"the fit with $(nstates) states cannot be evaluated at $(w) rad/s: its realization is not finite there. Fit with fewer poles, or over a narrower band."))
    return F
end

# The unknowns of the enforcement, gathered by output port: a
# perturbation of `S[i,j]` reads row `i` of `dC` and the single entry
# `dD[i,j]`, so the normal matrix `H` is block diagonal with `n` blocks
# of `nz + n` rather than one dense block of `n nz + n^2`, and the
# blocks are coupled only through the constraints (Gustavsen, IEEE
# Transactions on Electromagnetic Compatibility 67(3), 2025, section
# III-C).
function portblocks(n::Integer, nz::Integer)
    return [vcat([(k - 1)*n + i for k in 1:nz], [n*nz + (i - 1)*n + j for j in 1:n])
        for i in 1:n]
end

# The multipliers of the least change which meets every linearized
# passivity constraint at once: with `H` the normal matrix of the fit
# over the samples and the constraints `Ceq delta <= ceq`, the
# correction is `delta = -H^-1 Ceq' mu` and `mu` minimizes the dual
# `mu' G mu / 2 + ceq' mu` over `mu >= 0`, with `G = Ceq H^-1 Ceq'`, a
# nonnegative least squares problem solved by Lawson and Hanson's
# algorithm: admit the constraint whose multiplier most wants to be
# positive, solve on the working set, and take the longest step toward
# that solution which keeps every multiplier nonnegative, retiring
# those which reach zero. The gradient `G mu + ceq` is `ceq - Ceq
# delta`, so an entry of it below zero is exactly a constraint the
# current correction still violates, and the sweep for one is the
# sweep for the other.
#
# Termination says every multiplier is nonnegative and no constraint
# asks to be admitted, not that the constraints can all be met: `Ceq`
# has far fewer rows than columns here and they are generically
# independent, so they can be, but where they cannot the dual has no
# minimum, the Gram matrix is singular along the direction that shows
# it, and the least norm multipliers satisfy the gradient test without
# meeting the constraints. The caller therefore tests the correction
# against the constraints rather than trusting the flag alone.
function lawsonhanson(G::AbstractMatrix, c::AbstractVector; admittol::Real = 1e-12)
    m = length(c)
    μ = zeros(m)
    free = falses(m)
    # relative to the constraints, not to one: a correction to a block
    # with little loss has them at 1e-10, and a threshold with a unit
    # floor admits nothing and stops immediately
    scale = maximum(abs, c; init = 0.0)
    cap = 4*m + 40
    for _ in 1:cap
        g = G*μ .+ c
        admit, worst = 0, -admittol*scale
        for k in 1:m
            free[k] && continue
            g[k] < worst && ((admit, worst) = (k, g[k]))
        end
        admit == 0 && return μ
        free[admit] = true
        for _ in 1:cap
            P = findall(free)
            zp = gramsolve(G[P, P], view(c, P))
            all(isfinite, zp) || return μ
            if minimum(zp) >= 0
                fill!(μ, 0.0)
                μ[P] .= zp
                break
            end
            # the longest step toward `zp` which keeps every multiplier
            # nonnegative; at least one reaches zero and retires
            α = Inf
            for (q, k) in enumerate(P)
                zp[q] < 0 && (α = min(α, μ[k]/(μ[k] - zp[q])))
            end
            isfinite(α) || return μ
            for (q, k) in enumerate(P)
                μ[k] += α*(zp[q] - μ[k])
            end
            for k in P
                if μ[k] <= 0
                    μ[k] = 0.0
                    free[k] = false
                end
            end
            # a degenerate step which retires the constraint just admitted
            # makes no progress, and repeating it would not either
            free[admit] || return μ
        end
    end
    return μ
end

# The sweep above is finite and exact where the working set is
# independent. Where it is not, it can stop on a working set whose
# equalities contradict each other -- two constraints with parallel
# gradients and different bounds cannot both hold with equality -- and
# the least norm solve then returns a compromise satisfying neither,
# unnoticed because the sweep tests the gradient only of the entries
# it left out: on `min x'x/2` subject to `2x <= -2` and `x <= -1.5`
# the answer is `x = -1.5` with only the second constraint active, and
# the compromise is `x = -1.1`. Coordinate descent on the same dual
# finishes it: each step `mu_k <- max(0, mu_k - g_k/G_kk)` is exact
# along its own coordinate, asks nothing of `G` beyond positive
# semidefiniteness, converges on a dependent working set where a solve
# of the equalities cannot, and finds nothing to do where the sweep is
# already right. It runs to the conditions, not to a count.
function dualactiveset(G::AbstractMatrix, c::AbstractVector; reltol::Real = 1e-10)
    m = length(c)
    μ = lawsonhanson(G, c)
    # nonnegative multipliers, a zero gradient where one is positive, and
    # a nonnegative gradient where one is zero: all of the conditions, not
    # the half the sweep checks
    residual = ν -> begin
        g = G*ν .+ c
        worst = 0.0
        for k in 1:m
            worst = max(worst, ν[k] > 0 ? abs(g[k]) : max(0.0, -g[k]))
        end
        worst
    end
    tol = reltol*max(maximum(abs, c; init = 0.0), floatmin(Float64))
    for _ in 1:200
        residual(μ) <= tol && break
        for k in 1:m
            G[k, k] > 0 || continue
            gk = c[k]
            for j in 1:m
                gk += G[k, j]*μ[j]
            end
            μ[k] = max(0.0, μ[k] - gk/G[k, k])
        end
    end
    return μ, residual(μ) <= tol
end

# `G z = -c` on the working set. Several frequencies, or several singular
# values at one frequency, can ask for the same direction and leave the
# Gram matrix singular; the least norm multipliers are the right answer
# there, and every choice among them gives the same correction.
function gramsolve(Gp::AbstractMatrix, cp::AbstractVector)
    S = Symmetric(Matrix(Gp))
    F = cholesky(S; check = false)
    issuccess(F) && return F \ (-collect(cp))
    return pinv(Matrix(S))*(-collect(cp))
end

# How far to contract a block toward an anchor: the largest step `t`
# for which `measure(t)`, the quantity to bring under `target`, does
# so. Along the path `S0 + t (S - S0)` the response is affine in `t`,
# so its norm is convex and the steps meeting the target are an
# interval containing zero whose end the search finds. The step is
# found on the path itself because no formula stands in for it:
# scaling by `1/M`, with `M` a bound on the norm, answers only the
# `S0 = 0` case, and the triangle bound `t <= (1 - a)/(M - a)`, with
# `a` the anchor's own norm, is sufficient but not necessary -- at
# `a = 1` it says nothing, yet `S = -1.001 + 2.001/(s + 1)` anchored
# at `S(0) = 1` contracts at `t = 2/2.001` to the all pass
# `(1 - s)/(1 + s)`.
function contractionstep(measure, target::Real; tol::Real = 1e-12)
    f1 = measure(1.0) - target
    f1 <= 0 && return 1.0
    f0 = measure(0.0) - target
    f0 <= 0 || return 0.0
    # The measure is convex, so the chord from a point under the target
    # to one over it crosses the target on the near side of the true
    # boundary: every chord step is a valid step, the steps approach the
    # boundary from below, and where the measure is linear, as it is
    # toward an anchor of zero, the first chord lands on it. A chord
    # which makes no progress falls back to bisecting the bracket, so
    # the search cannot stall. `tol` bounds the step gained per
    # iteration and the bracket at which the search stops; the shortfall
    # of the returned step from the boundary is what the caller pays as
    # needless contraction, so it is kept small.
    lo, flo, hi, fhi = 0.0, f0, 1.0, f1
    for _ in 1:60
        t = fhi > flo ? lo - flo*(hi - lo)/(fhi - flo) : (lo + hi)/2
        lo < t < hi || (t = (lo + hi)/2)
        ft = measure(t) - target
        if ft <= 0
            gain = t - lo
            lo, flo = t, ft
            # a step landing exactly on the target is the boundary
            (gain <= tol || flo == 0) && break
        else
            hi, fhi = t, ft
        end
        (hi - lo) <= tol && break
    end
    return lo
end

function enforcepassivity(A, B, C, D, ws; atol = 1e-8, rounds::Int = 20, margin::Real = 1e-6,
        dc::Union{Nothing,AbstractMatrix} = nothing, regularization::Real = 1e-12,
        slackreltol::Real = 1e-8, contractiontol::Real = 1e-12,
        scalelimit::Real = passivityscalelimit, dcchecktol::Real = 1e-9,
        warnfloor::Real = 1e-9, bandpoints::Int = 51)
    n, nz = size(D, 1), size(A, 1)
    level = 1 + atol/2
    # `S(inf) = D`, which no perturbation over a band can reach: a
    # feedthrough above one leaves the norm attained at infinite frequency,
    # where the band search has nothing to work with. Bring it under first.
    if opnorm(D) > 1
        # only a feedthrough which is nearly contractive is scaled under
        # one. Scaling one which is far above it would divide the whole fit
        # by that factor and return a block which is passive because it
        # transmits nothing, which is worse than refusing.
        opnorm(D) <= 1 + scalelimit || throw(ArgumentError(
            lazy"the fit is not passive at infinite frequency: its feedthrough alone has a largest singular value of $(opnorm(D)), which no perturbation over the band can bring under one. Fit with fewer poles, or over a narrower band."))
        # Scaling the whole fit multiplies its value at zero frequency by
        # the same factor. A block which states one is contracted toward
        # that statement instead, which holds it exactly: `S -> S0 +
        # t (S - S0)` is `C -> t C` and `D -> (1 - t) S0 + t D`, and at
        # zero it gives `S0` back.
        anchor0 = isnothing(dc) ? zeros(size(D)) : Float64.(dc)
        t = contractionstep(τ -> opnorm((1 - τ) .* anchor0 .+ τ .* D), 1.0)
        if t > 0
            C, D = t .* C, (1 - t) .* anchor0 .+ t .* D
        elseif opnorm(D) > 1 + atol
            # no step exists, the statement being of unit norm itself, and
            # the excess is more than the tolerance: the two cannot both
            # hold
            throw(ArgumentError(
                lazy"the fit's feedthrough has a largest singular value of $(opnorm(D)) and the value the block states at zero frequency has one of $(opnorm(Float64.(dc))): contracting toward a statement which is itself of unit norm cannot bring anything under one. Fit with fewer poles, over a narrower band, or without stating the value at zero."))
        end
        # an excess within the tolerance and no step to take is a
        # feedthrough which is passive to the tolerance asked for, and is
        # left where it is
    end
    # The perturbation is confined to the changes which leave the value
    # at zero frequency where the residue solve put it. That value is
    # `S(0) = D + C X0` with `X0 = (-A)^-1 B`, so a change preserves it
    # exactly when `dC X0 + dD = 0`; restricted to one output port's
    # unknowns, its row of `C` then its row of `D`, that reads
    # `[X0' I] w = 0`, the same `n` by `nz + n` block for every port,
    # and the perturbation is confined to its nullspace, `nz` free
    # directions per port instead of `nz + n`. Eliminating rather than
    # penalizing is what makes the statement hold exactly: the least
    # squares is solved in the eliminated coordinates, so the answer is
    # the least change which holds it, not the least change which
    # nearly does.
    dcbasis, X0 = nothing, nothing
    if !isnothing(dc)
        X0 = (-A) \ B
        dcbasis = nullspace(hcat(transpose(X0), Matrix{Float64}(I, n, n)))
        size(dcbasis, 2) == nz || throw(ArgumentError(
            lazy"the zero frequency condition leaves $(size(dcbasis, 2)) free directions per port where $(nz) were expected; report this."))
    end
    blocks = portblocks(n, nz)
    blockfactors = nothing
    rf = resolventfactors(A, B)
    for roundindex in 1:rounds
        bands = sampledviolations(rf, C, D, ws; atol = atol)
        # the rounds refine; the guarantee below is what every path ends at
        isempty(bands) && break
        # the worst point of each band; constraining every grid point
        # above the level instead over-constrains the perturbation and
        # keeps the rounds from converging
        CZ = C*rf.Z
        points = Float64[]
        for (w1, w2) in bands
            grid = range(w1, w2; length = bandpoints)
            k = argmax([opnorm(transferat(rf, CZ, D, w, nz)) for w in grid])
            push!(points, grid[k])
        end
        # the perturbation of C and D, `delta S(i w) = delta C (i w I - A)^(-1) B + delta D`,
        # as a real vector; its effect at the samples for the norm and at
        # the points for the constraints
        nunk = n*nz + n*n
        function perturbation(w)
            X = resolventat(rf, w)
            # delta S[i, j] = sum_k delta C[i, k] X[k, j] + delta D[i, j]
            M = zeros(ComplexF64, n*n, nunk)
            for i in 1:n, j in 1:n
                row = (i - 1)*n + j
                for k in 1:nz
                    M[row, (k - 1)*n + i] = X[k, j]
                end
                M[row, n*nz + (i - 1)*n + j] = 1.0
            end
            return M
        end
        # the norm: the fit over the samples, as the normal matrix
        # accumulated sample by sample, never the stacked samples. It reads
        # only `A` and `B`, which the enforcement leaves alone, so it is
        # built once and kept across the rounds.
        if isnothing(blockfactors)
            # Restricted to port `i`'s own unknowns the perturbation
            # reads `[X(w)' I]`, which does not depend on `i`: the
            # ports differ only in which columns they occupy, so one
            # block of `nz + n` is accumulated and one factorization
            # serves them all (Gustavsen, section III-C, equal
            # weighting).
            nb = nz + n
            Hb = zeros(nb, nb)
            G = zeros(ComplexF64, n, nb)
            for w in ws
                X = resolventat(rf, w)
                fill!(G, 0)
                for j in 1:n
                    for k in 1:nz
                        G[j, k] = X[k, j]
                    end
                    G[j, nz + j] = 1.0
                end
                Hb .+= real.(G'*G)
            end
            for i in 1:nb
                Hb[i, i] += regularization
            end
            # in the eliminated coordinates when a value at zero is held
            isnothing(dcbasis) || (Hb = transpose(dcbasis)*Hb*dcbasis)
            blockfactors = cholesky(Symmetric(Hb))
        end
        # The constraints, one per singular value which stands above the
        # level at each point, not only the largest: a perturbation which
        # pulls the largest under can push the next one over, and the
        # rounds are then spent trading them. The relation between a
        # perturbation and a singular value is the first order one,
        # `d sigma_j = Re(u_j' dS v_j)`.
        rows = Vector{Float64}[]
        rhs = Float64[]
        for w in points
            Sw = transferat(rf, CZ, D, w, nz)
            F = svd(Sw)
            M = perturbation(w)
            for j in 1:n
                F.S[j] > level || continue
                u, v = F.U[:, j], F.V[:, j]
                row = zeros(ComplexF64, nunk)
                for i in 1:n, l in 1:n
                    row .+= conj(u[i])*v[l] .* view(M, (i - 1)*n + l, :)
                end
                push!(rows, real.(row))
                push!(rhs, 1 - margin - F.S[j])
            end
        end
        isempty(rows) && break
        Ceq = reduce(vcat, transpose.(rows))
        ceq = rhs
        # The constraints are inequalities, not equalities: forcing
        # each to hold with equality would drag a singular value which
        # only just crosses one the whole way to the target, and with
        # every violating singular value constrained at once the extra
        # equalities would over-determine the perturbation. A
        # multiplier of the wrong sign says exactly that a constraint
        # is met without being pushed, and it is dropped.
        m = length(ceq)
        # `H^-1 Ceq'` and the Gram matrix `Ceq H^-1 Ceq'` of the dual,
        # built once for all of this round's constraints, one output
        # port block at a time, the rows reduced through the zero
        # frequency nullspace where there is one
        blockrows = i -> isnothing(dcbasis) ? Matrix(view(Ceq, :, blocks[i])) :
            Matrix(view(Ceq, :, blocks[i]))*dcbasis
        Cb = [blockrows(i) for i in 1:n]
        G = zeros(m, m)
        Wb = [blockfactors \ Matrix(transpose(Cb[i])) for i in 1:n]
        for i in 1:n
            G .+= Cb[i]*Wb[i]
        end
        μ, solved = dualactiveset(Symmetric(G), ceq)
        δ = zeros(nunk)
        for i in 1:n
            ui = .-(Wb[i]*μ)
            δ[blocks[i]] .= isnothing(dcbasis) ? ui : dcbasis*ui
        end
        # The multipliers are nonnegative by construction; the
        # correction must also meet the constraints it was built from,
        # or the round leaves the fit as it stands for the guarantee
        # below to scale. The test is measured against the size of the
        # numbers involved, since a correction to a block with little
        # loss has constraints of order 1e-9, where an absolute floor
        # would accept anything.
        residuals = Ceq*δ
        slack = maximum(residuals .- ceq; init = -Inf)
        slacktol = slackreltol*max(maximum(abs, ceq; init = 0.0),
                            maximum(abs, residuals; init = 0.0), eps())
        # the elimination is exact, so this is a check on the arithmetic
        isnothing(dcbasis) || (maximum(abs, reshape(δ[1:n*nz], n, nz)*X0 .+
            transpose(reshape(δ[n*nz + 1:end], n, n)); init = 0.0) <= dcchecktol*(1 + maximum(abs, δ))) ||
            throw(ArgumentError("the passivity correction moved the stated value at zero frequency; report this."))
        (solved && slack <= slacktol) || break
        # the unknowns were laid out as delta C[i, k] at (k - 1) n + i and
        # delta D[i, j] at (i - 1) n + j
        C = C .+ reshape(δ[1:n*nz], n, nz)
        D = D .+ transpose(reshape(δ[n*nz + 1:end], n, n))
        # a singular KKT system leaves a correction which is not finite,
        # and the failure would otherwise surface as an error from the
        # next factorization built from `C` and `D`
        (all(isfinite, C) && all(isfinite, D)) || throw(ArgumentError(
            lazy"the passivity enforcement diverged in round $(roundindex) with $(nz) states: its least norm correction is not finite. Fit with fewer poles, or over a narrower band."))
    end
    # The rounds refine; this guarantees. A block whose largest
    # singular value exceeds one is active and an interconnection built
    # on it can grow without bound, so a fit which is nearly passive is
    # contracted the rest of the way rather than refused: toward
    # nothing where the block states no value at zero, toward the
    # statement where it makes one. Only a fit far from passive is
    # refused, since scaling that far would describe a block which is
    # mostly loss. `worst` is the largest singular value the norm
    # search evaluated and `ceiling` the level it established nothing
    # reaches; the decision is taken on `worst`, so what is guaranteed
    # is passivity to within the search's tolerance, `2 rtol`, and not
    # passivity outright. Deciding on `ceiling` would scale every block
    # whose feedthrough is unitary and buy nothing: a peak which
    # exceeds a level by less than roundoff brings its two crossings
    # together into a nearly double eigenvalue that leaves the axis, so
    # the level cannot be brought to one and tightening the tolerance
    # stops working before it gets there.
    worst, _, ceiling = hinfnorm(A, B, C, D)
    isfinite(ceiling) || throw(ArgumentError(
        lazy"the fit's largest singular value could not be bounded after $(rounds) rounds: the level set search did not converge, and the largest value it saw was $(worst). Fit with fewer poles, or over a narrower band."))
    # a statement of unit norm at zero frequency -- a through, an
    # open, a short -- pins the norm of any fit meeting it at one, so
    # nothing can be contracted away and such a fit is accepted at one
    # to the tolerance it is about to be validated against, which is
    # what a lossless block is accepted at in any case; a statement
    # with room under one leaves room, and the contraction below uses
    # it
    worst <= 1 && return A, B, C, D
    worst <= 1 + scalelimit || throw(ArgumentError(
        lazy"the fit could not be made passive in $(rounds) rounds: its largest singular value is still $(worst), too far above one to scale under it without describing a block which is mostly loss. Fit with fewer poles, or over a narrower band."))
    # The worst point moves as the block is scaled and the norm is
    # found by a search, so a single step to a margin guessed in
    # advance can leave the next measurement a shade above one: the
    # scaling is repeated against fresh measurements until the block
    # measures contractive, which takes one further step at most. A
    # block which states a value at zero is contracted toward that
    # statement, `C -> t C` and `D -> (1 - t) S0 + t D`, which leaves
    # the statement exactly where it is; that contraction is weaker,
    # pulling toward a matrix of norm up to one rather than toward
    # nothing, so it too is measured rather than trusted.
    S0 = isnothing(dc) ? zeros(size(D)) : Float64.(dc)
    anchor = opnorm(S0)
    Cbefore, Dbefore = copy(C), copy(D)
    s, measured = 1.0, worst
    for _ in 1:8
        measured <= 1 && break
        previous = measured
        # the largest step along the path which measures contractive,
        # found on the norm itself rather than from a bound on it
        step = contractionstep(1.0; tol = contractiontol) do τ
            # at no step at all the block is the anchor, a constant, whose
            # norm is its own; the pencil of a realization with no output
            # coupling is degenerate and does not need solving
            τ <= 0 && return anchor
            try
                first(hinfnorm(A, B, τ .* C, (1 - τ) .* S0 .+ τ .* D))
            catch e
                # a step the norm cannot be measured at is a step to move
                # away from, not one to accept
                (e isa LinearAlgebra.LAPACKException || e isa SingularException) || rethrow()
                Inf
            end
        end
        step > 0 || break
        C, D = step .* C, (1 - step) .* S0 .+ step .* D
        s *= step
        measured, _ = hinfnorm(A, B, C, D)
        measured >= previous && break
    end
    # The contraction is bounded by the same principle as the scaling of
    # the feedthrough: a fit within `scalelimit` of passive needs a
    # contraction of about `scalelimit` to repair, and the repeated
    # measurement can spend that much again, so a contraction below
    # `1 - 2 scalelimit` describes a block which is mostly the anchor
    # rather than the data -- passive because nothing of the fit remains,
    # which is worse than refusing. The norm search decides each step on
    # a lower bound, so without this floor a step of no size at all can
    # measure contractive against an anchor of unit norm and be accepted.
    s >= 1 - 2*scalelimit || throw(ArgumentError(isnothing(dc) ?
        lazy"making the fit passive took a contraction to $(s) of the data, which describes a block which is mostly loss rather than the data. Fit with fewer poles, or over a narrower band." :
        lazy"making the fit passive took a contraction to $(s) of the data, which describes a block which is mostly the value it states at zero frequency rather than the data. Fit with fewer poles, over a narrower band, or without stating the value at zero."))
    # one is the floor for a fit anchored at a statement of unit norm,
    # not a target it can be brought under; it is accepted there to the
    # tolerance it is validated against
    (measured <= 1 || (anchor > 0 && measured <= 1 + atol)) || throw(ArgumentError(
        isnothing(dc) ?
        lazy"contracting the fit by $(s) left its largest singular value at $(measured), which should not happen; report this." :
        lazy"contracting the fit toward the value it states at zero frequency, by $(s), left its largest singular value at $(measured); the statement has a norm of $(anchor). Fit with fewer poles, over a narrower band, or without stating the value at zero."))
    if worst - 1 > warnfloor
        # report what the contraction moved, measured over the samples:
        # a contraction toward a nonzero anchor does not add loss
        # uniformly, so no scalar formula describes it
        Xs = [(im*w*I - A) \ B for w in ws]
        moved = maximum(k -> opnorm((C - Cbefore)*Xs[k] .+ (D - Dbefore)), eachindex(ws))
        @warn "the rational fit was contracted by $(s) to make it passive, which moved its response by up to $(moved) over the samples; compare the fitted dissipation against the data's before trusting the noise of this block."
    end
    return A, B, C, D
end

# The bands of angular frequency where the largest singular value
# stands above one, found by sweeping the fit on a logarithmic grid.
# The exact pencil of the passivity test is what the guarantee at the
# end of the enforcement rests on, but it is a generalized eigenproblem
# of twice the state dimension and costs far more than the sweep, which
# is one triangular solve per point once the resolvent is factored; the
# sweep's risk is a violation narrower than the grid spacing, which the
# exact norm at the end still catches (Gustavsen, IEEE Transactions on
# Electromagnetic Compatibility 67(3), 2025, section X-B).
function sampledviolations(rf::ResolventFactors, C, D, ws; atol = 1e-8, density::Int = 12,
        pad::Real = 10.0)
    level = 1 + atol/2
    n, nz = size(D, 1), size(rf.T, 1)
    # the logarithmic grid starts from the lowest positive sample, a band
    # which begins at zero having no logarithm; zero itself is on the grid
    # in either case, as its own point
    lo = ws[something(findfirst(>(0), ws), lastindex(ws))]
    hi = last(ws)
    grid = vcat(0.0, exp.(range(log(lo/pad), log(hi*pad); length = density*nz)))
    CZ = C*rf.Z
    above = falses(length(grid))
    for (k, w) in enumerate(grid)
        above[k] = opnorm(transferat(rf, CZ, D, w, nz)) > level
    end
    bands = Tuple{Float64,Float64}[]
    k = 1
    while k <= length(grid)
        if above[k]
            j = k
            while j < length(grid) && above[j + 1]
                j += 1
            end
            # widen by one spacing on each side, so the true crossing lies
            # inside the band the search then refines
            push!(bands, (grid[max(k - 1, 1)], grid[min(j + 1, length(grid))]))
            k = j + 1
        else
            k += 1
        end
    end
    return bands
end

