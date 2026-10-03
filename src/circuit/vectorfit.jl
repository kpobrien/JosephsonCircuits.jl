# Fitting a rational scattering block to sampled scattering data: common
# poles by the relaxed pole relocation of Gustavsen and Semlyen, residues
# and feedthrough by least squares at the settled poles, a real state
# space realization, and passivity enforcement by the least change of the
# residues and the feedthrough that brings every singular value back
# under one. Nothing here needs a package beyond linear algebra.

# A pole is real when its imaginary part is below this fraction of its
# magnitude. One shared constant rather than a keyword: the basis, the
# residue map, the real pole matrix, the relocation's pairing and the
# realization must all classify a pole the same way, or a conjugate
# pair splits across representations.
const realpoletolerance = 1e-6
isrealpole(a) = abs(imag(a)) <= realpoletolerance*abs(a)

# A singular value of the residue basis below this fraction of the
# largest is a direction the samples do not determine, which the residue
# solve leaves out. One shared constant rather than a keyword: the
# relocation, the pruning and the final solve all judge a fit by the
# residues this solve returns, so they must leave out the same directions.
const fitranktolerance = 1e-12

"""
    VectorFitting(; iterations = 30, stallpatience = 5, start = :linear,
        pruneslack = 0.05, ranktol = 1e-12, zerotol = 1e-12)

The numerical parameters of the vector fit of [`RationalScattering`](@ref)
at a number of poles.

- `iterations` is the most rounds of pole relocation. The relocation ends
  sooner when `stallpatience` rounds in a row bring the largest deviation
  over the samples no lower than a thousandth of itself under the least so
  far, which is a fit that has settled, or one circling, or when the fit
  meets its samples to their roundoff. The rounds go on while the fit
  improves, however far its poles move, and the iterate returned is the
  one measured to fit best.
- `start` is where the relocation starts. `:linear` spreads complex
  pairs evenly over the positive band and `:log` evenly over the
  logarithm of its frequency, the lowest positive sample standing in for
  zero, with a real pole at the centre of the band in the same spacing
  when the count is odd, each damped by a hundredth of its frequency.
  `:linlog` spreads about half of the pairs each way, the logarithmic ones
  short of the band's ends, with the real pole at the logarithmic centre,
  as the linlogcmplx start of Gustavsen's VFdriver does, and five poles or
  fewer as `:log`. Where the relocation ends can depend on where it
  starts, and no spacing is the better on every data set. A vector of
  poles in rad/s starts from those, as many as the poles fitted: each
  finite and in the open left half plane, and each complex one given with
  its conjugate. The order search fits many orders from a spacing and
  takes no vector.
- `pruneslack` bounds how much the fit error may grow over the whole
  pruning of needless poles, as a fraction of the error before any pole
  was dropped, in the largest deviation over the samples and in the rms
  deviation over every entry and sample, each against its own value:
  the largest deviation alone can be pinned by one feature of the data,
  a seam or a spike, and leave the rest of the fit free to worsen. It
  bounds the pruned fit, not the returned block: the repair of the
  constant term and the passivity enforcement come after.
- A pole takes as many states as its residue has singular values above
  `ranktol` of the largest, decided on the fit before its passivity
  enforcement: the residue is taken within the space of the singular
  vectors kept, which the enforcement changes it within and the
  realization spans, so that the block returned is the one made
  passive.
- In the fit of a [`LinearizedScattering`](@ref) block, the cosine or
  sine part of a harmonic whose samples all stand below `zerotol` of the
  harmonic's largest is realized as zero, with no state.

The relocation shares its columns among the threads Julia runs with,
where their least squares are large enough for that to pay, each thread
running BLAS on one thread, and gives the same poles on any number of
them; on a single Julia thread BLAS keeps its own threads,
whose products round differently in the last digits. A relocation on
data its poles cannot follow exactly carries such a difference into a
different set of poles, as it does any other change of rounding, so that
a fit is the same on any number of threads only where BLAS runs on one
thread throughout.
"""
struct VectorFitting
    iterations::Int
    stallpatience::Int
    start::Union{Symbol,Vector{ComplexF64}}
    pruneslack::Float64
    ranktol::Float64
    zerotol::Float64
end
function VectorFitting(; iterations::Integer = 30, stallpatience::Integer = 5, start = :linear,
        pruneslack::Real = 0.05, ranktol::Real = 1e-12, zerotol::Real = 1e-12)
    iterations >= 1 || throw(ArgumentError("give at least one relocation iteration."))
    stallpatience >= 1 || throw(ArgumentError("stallpatience must be at least one round."))
    for (name, x) in ((:pruneslack, pruneslack), (:ranktol, ranktol), (:zerotol, zerotol))
        (isfinite(x) && x >= 0) || throw(ArgumentError(lazy"$(name) must be finite and nonnegative."))
    end
    return VectorFitting(Int(iterations), Int(stallpatience), startingpoles(start), pruneslack, ranktol, zerotol)
end
# The start of the relocation, checked: a spacing by name, or poles
# which can start it, the pairs side by side as the real basis reads them,
# each complex pole given with its conjugate to within `conjtol` of its
# size
function startingpoles(start; conjtol::Real = 1e-8)
    if start isa Symbol
        start in (:linear, :log, :linlog) || throw(ArgumentError(
            lazy"start is :linear, :log, :linlog or a vector of starting poles in rad/s, not :$(start)."))
        return start
    end
    start isa AbstractVector{<:Number} || throw(ArgumentError(
        "start is :linear, :log, :linlog or a vector of starting poles in rad/s."))
    given = ComplexF64.(start)
    isempty(given) && throw(ArgumentError("give at least one starting pole."))
    all(a -> isfinite(a) && real(a) < 0, given) || throw(ArgumentError(
        "the starting poles must be finite and in the open left half plane: a pole on or right of the imaginary axis is no stable pole to start from."))
    # a complex pole stands with its conjugate, as every pole of a real
    # rational function does; the pole in the upper half plane is kept
    # with its exact conjugate
    upper = filter(a -> !isrealpole(a) && imag(a) > 0, given)
    lower = filter(a -> !isrealpole(a) && imag(a) < 0, given)
    used = falses(length(lower))
    paired = length(upper) == length(lower) && all(upper) do a
        l = findfirst(l -> !used[l] && abs(lower[l] - conj(a)) < conjtol*abs(a), eachindex(lower))
        isnothing(l) || (used[l] = true)
        return !isnothing(l)
    end
    paired || throw(ArgumentError(
        "the complex starting poles must come in conjugate pairs: give the conjugate of each."))
    poles = ComplexF64[]
    for a in given
        if isrealpole(a)
            push!(poles, complex(real(a), 0.0))
        elseif imag(a) > 0
            push!(poles, a, conj(a))
        end
    end
    return poles
end

# the time, in seconds, a fit's certified sweeps may run where the caller
# sets none: the enforcement's `maxtime` by default (see
# PassivityEnforcement), and the check of a fit made without it (see
# checkrawpassive)
const sweepmaxtime = 300.0

"""
    PassivityEnforcement(; margin = 1e-6, rounds = 60, scalelimit = 1e-2,
        maxtime = 300.0, regularization = 1e-12)

The numerical parameters of the passivity enforcement of a fit by
[`RationalScattering`](@ref).

- A round searches the fit for bands where its dissipation `I - S'S`
  falls below `-atol`, the tolerance the block's samples and constant
  term are held to, its largest singular value above `sqrt(1 + atol)`,
  about one plus half of `atol`, on a logarithmic grid of points per pole
  reaching past the band on either side, and at each complex pole and
  half a width and a width either side of it, a pole found above the
  level adding a band of half widths either side. Where the grid finds
  none, the certified sweep below decides: it establishes the fit
  passive, or finds the bands the grid passed between, as where a
  correction holds the response at one beside the points it constrains
  and a bump narrower than the grid's spacing rises. It constrains the
  worst of a few dozen points across each band, and infinite frequency
  where the constant term stands above the level; the points of the
  earlier rounds, infinite frequency among them, are constrained again
  at their largest singular value whether or not it stands above one, so
  that a correction cannot push a point back over one that an earlier
  one brought under. It changes the residues and the constant term by
  the least amount, in the norm of the fit over the samples, each entry
  weighed by its weight where the fit takes weights (see
  [`RationalScattering`](@ref)), that brings every singular value above
  one there to `1 - margin`, over at most `rounds` rounds. `margin`
  must be well below the block's own dissipation `1 - sigma`, or the
  enforcement writes over the loss it should preserve.
- The certified sweep settles every interval of the frequency axis,
  infinity included: it bounds the fit's dissipation `I - S'S` from below
  across an interval, by its value and derivative at the interval's
  centre and the distances of the poles from the interval, and halves an
  interval its bound does not settle, so that no band passes between its
  points. A sweep which settles every interval under that level
  establishes that the fit is passive to it; otherwise its bands are
  those where the fit stands above the level, or within the roundoff of
  it, and the worst point of each is constrained. The dissipation of a
  lossless fit vanishes, residues and all, and a bound over the whole
  axis from its residues and constant settles such a fit, and its norm,
  before any interval. Elsewhere the intervals it needs grow as the
  margin by which the fit stands under the level shrinks, over the bands
  where it stands near it. Where the fit is lossless in some directions
  at every frequency and lossy in the others, a bound by the
  dissipation's second derivative as well settles an interval the first
  order would halve many times, and they grow as the inverse cube root
  of the margin rather than its inverse square root. Each generation of
  intervals is shared among the threads Julia runs with, where the form
  is large enough for that to pay, each thread running BLAS on one
  thread, and gives the same result on any number of them; on a single
  Julia thread BLAS keeps its own threads, whose products round
  differently in the last digits, so that the result is the same there
  only where BLAS runs on one thread as well. An interval costs as the
  cube of the ports, so the sweeps are bounded in time rather than in
  intervals: the sweeps and norm searches of the enforcement stop at
  their first reading of the clock `maxtime` seconds or more after it
  began, and a fit they have not settled by then is refused. Only they
  are timed, not the grid searches and corrections between them, so that
  the enforcement can run past `maxtime`, and the bound over the whole
  axis, which takes no interval, settles a lossless fit however little
  time is left. A fit settled is the same whatever the time left, and
  `maxtime = Inf` lets the sweeps run to the floating-point resolution.
- The least change is found in the coordinates of the residues at the
  frequency scaled to the band: each residue's change within the space
  of its rows which its states span (`ranktol` of
  [`VectorFitting`](@ref)), so that the block keeps its states, measured
  by the Gram matrix of the pole basis over the samples, which carries no
  unit, and weighted by each entry's weights where the fit takes them, a
  Gram matrix for each entry. Samples which are reciprocal to `atol` have
  the transpose of each constraint as well. Where every residue has full
  rank, and an entry and its transpose have the same weights, the least
  change is then reciprocal, to the accuracy of its solve, which above
  the band, where the samples hardly determine the change, can fall short
  of `atol`; a residue of lower rank changes within the space of its
  rows, which its transpose need not lie in, and the fit can lose its
  reciprocity by as much as the change.
- The normal matrix of that least squares carries `regularization` times
  its mean diagonal entry on its diagonal, which must be positive: it
  keeps the change finite in the directions the samples hardly
  determine, above the band, and a larger one keeps it smaller there at
  the cost of a larger change over the samples.
- A result the certified sweep has not established passive is measured
  by its largest singular value over every frequency, found by the same
  bounds: the response at each interval's centre raises a lower bound,
  and an interval whose bound is under that bound times one and a
  relative tolerance is settled, as is one halving cannot settle, its
  centre within the roundoff of that level, its bound's own width within
  the roundoff, or its half width at the floating-point resolution of its
  centre, whose bound then counts toward a ceiling no frequency reaches,
  so that a tolerance under the roundoff gives a ceiling looser than
  itself. The fit is accepted where that ceiling stands under the level,
  never on the lower bound alone. The tolerance is `1e-8`, or half the
  margin between one and the level where that is finer, so that the
  ceiling of a fit which reaches one, as a fit held at a value of unit
  norm at zero frequency does, stands under the level.
- A fit the rounds leave within `scalelimit` above one, by what it
  reaches, is contracted the rest of the way: toward nothing until its
  ceiling stands under one, by the step the ceiling gives where that
  measures under one and otherwise by a step searched for on the path;
  toward a value stated at zero frequency, whose own norm may be one, by
  such a step until its ceiling stands under the level. One further
  above is refused, since a contraction that large describes a block
  which is mostly loss rather than the data, as is one whose ceiling a
  contraction cannot bring under its target; a contraction which moves
  the response over the samples by more than `atol` is warned of.
"""
struct PassivityEnforcement
    margin::Float64
    rounds::Int
    scalelimit::Float64
    maxtime::Float64
    regularization::Float64
end
function PassivityEnforcement(; margin::Real = 1e-6, rounds::Integer = 60, scalelimit::Real = 1e-2,
        maxtime::Real = sweepmaxtime, regularization::Real = 1e-12)
    (isfinite(margin) && 0 <= margin < 1) || throw(ArgumentError(
        "margin must be finite and in [0, 1): it is how far below one a singular value is put."))
    rounds >= 1 || throw(ArgumentError("give at least one round of passivity enforcement."))
    # a contraction to below `1 - 2 scalelimit` of the data is refused (see
    # enforcepassivity), so a limit of a half or more would refuse nothing
    (isfinite(scalelimit) && 0 <= scalelimit < 0.5) || throw(ArgumentError(
        "scalelimit must be finite and in [0, 0.5): it is how far above one a fit may stand and still be contracted under it."))
    maxtime >= 0 || throw(ArgumentError("maxtime must be nonnegative: it is how long, in seconds, the certified sweeps may run."))
    (isfinite(regularization) && regularization > 0) || throw(ArgumentError(
        "regularization must be finite and positive: the normal matrix of the least change is singular in the directions the samples do not determine."))
    return PassivityEnforcement(margin, Int(rounds), scalelimit, Float64(maxtime), regularization)
end

"""
    RationalScattering(block::ScatteringParameters, npoles;
        frequencies = nothing, tol = 1e-2, atol = 1e-8,
        fitting = VectorFitting(), passivity = PassivityEnforcement(),
        delays = nothing, maxstates = defaultmaxstates(), weights = nothing)

A [`RationalScattering`](@ref) block fitted to the scattering data of
`block` at `npoles` common poles by vector fitting: starting poles,
spread over the band or given, are relocated by the relaxed iteration of
Gustavsen and Semlyen until they settle, and the residues of every entry
and the constant term follow by least squares, with the numerical
parameters of `fitting` (see [`VectorFitting`](@ref)).
The data is sampled at `frequencies` in Hz, by default a tabulated
block's own. Poles the data does not need drift out of the band or
coalesce and are pruned, so `npoles` is a budget rather than the order
returned; too few poles settle on a poor fit. The block returned is
measured against the samples, by its largest deviation from them in the
spectral norm as a fraction of the largest response, as the `tol`
method of the search measures it, and refused where that exceeds `tol`,
since a block quietly far from its data is worse than none; where no fit
with the poles the relocation settles on can come within `tol`, as a
weighted least squares at them bounds, it is refused before the
passivity enforcement, which keeps the poles. The result is a real state
space realization, stable by construction, with as many states per pole
as its residue has rank. The fitted block keeps the reference
impedances, grounding and noise model of `block`.

`weights`, an array of the samples' shape, `(nports, nports, K)` for `K`
frequencies, of positive numbers, weighs each entry at each sample: the
least squares of the relocation and of the residues, the errors the
pruning holds and `tol` bounds, and the least change of the passivity
enforcement all take each entry's deviation times its weight, and the
largest response is that of the samples so weighted. `1 ./ abs.(S)` of
the samples `S` fits every entry to the same relative accuracy, as the
individual element weighting of Gustavsen's VFdriver does, which small
entries such as an isolation or a crosstalk need, and a weight the same
for every entry at a sample weighs the samples. Only the weights' ratios
matter, and weights the same everywhere are no weights. Weighted, the
relocation runs on the entries, each with a least squares of its own,
rather than on their singular components, and so do the residues, so
that a fit costs as the square of the ports rather than as the
components; the entries are shared among Julia's threads.

The realization holds dense matrices of the states, so the memory a fit
takes grows as the square of its states. A fit of more than `maxstates`
states is refused before any of them is formed. By default that is the
most states whose matrices take half the memory, as
`Sys.total_memory()` reports it. Weighted, the passivity enforcement
holds a metric of its own for each entry, or where a residue has less
than full rank for each port, of every unknown of a port's row, and a
correction whose metric needs more than the states leave of that memory
is refused before the metric is formed.

A delay is not a rational function. A stable proper rational function's
impulse response starts at zero time, so a fit of a delay answers
throughout the delay rather than after it, and no pole budget removes
that: the fit follows the phase of `exp(-i w tau)` over the band it is
given and no further, and the turns above the band are what arrive
early. `delays`, one per port in seconds, takes a delay out of the data
before the fit, as a cable's is taken out before its fit: the entry
from port `q` to port `p` is fitted with `exp(i w (tau_p + tau_q))`
taken out, so the block returned is the device with a lossless line of
delay `tau_p` cut off each port, and is put back with a
[`TransmissionLine`](@ref) of that delay in cascade at the port, which
carries it exactly. A stated covariance is a correlation of the waves
the block emits and is rotated by the difference of the delays where
the scatter carries their sum, so delays which are all equal, and
uncorrelated ports, leave it alone; one they do turn is turned at the
frequency it is read at, whatever kind of provider states it. The
delays leave a block passive, a diagonal phase on each side of `S`
holding its singular values, and leave the zero frequency statement of
a `dcmodel` alone.

With `passivity`, a [`PassivityEnforcement`](@ref), the fit is perturbed
wherever its dissipation `I - S'S` falls below `-atol`, its largest
singular value above `sqrt(1 + atol)`, about one plus half of `atol`,
found on a grid over the band and at each resonance of the fit, or by a
sweep of the frequency axis which bounds the fit over every interval, by
the smallest change of its residues and constant that brings each such
point under one by its `margin`, over at most its `rounds` rounds, and
the result is established passive to `atol` by the sweep's bounds, at
every frequency and not only at samples. A peak the rounds leave under
`sqrt(1 + atol)` stands where it is, so a solve which needs a block
further under one asks for a smaller `atol`. A
fit the rounds leave within its `scalelimit` of one is contracted the
rest of the way, which moves its response by about that much, and a
warning says how far where that exceeds `atol`; one further above is
refused, since a contraction that large describes a block which is
mostly loss rather than the data.
`passivity = nothing` returns the raw fit, without the enforcement or
its repair of an active constant term, but still refuses one which the
same sweep finds above `sqrt(1 + atol)` at some frequency, or does not
settle within the enforcement's default `maxtime`, 300 s; enforcing
passivity with a larger `maxtime` gives the sweep longer. A
block which states its noise with a [`NoiseCovariance`](@ref) may be
active, so it is fitted as it is, with neither the enforcement nor the
validation, and `passivity` is moot; the fit is still stable by
construction, and the stated covariance is held to what the fitted
scattering matrix requires wherever a solver evaluates it.

If `block` states its zero frequency behavior, through the `dcmodel` it
was built with, the fit meets the statement exactly: the value at zero
is linear in the residues once the poles settle, so it is imposed as an
equality, on the residue fit and through the passivity enforcement. The
statement is the caller's to make and is not checked against the data,
which begins above zero and cannot check it. A statement far below the
band may need poles placed there, which a start given as a vector of
poles places. A statement of unit norm -- a through, an open, a short --
pins the norm of any fit meeting it at one, so such a fit is accepted at
one to `atol` rather than contracted, which would move the statement;
where the statement and passivity cannot both hold, the fit is refused.

The dissipation `I - S S'` is a difference of nearly equal quantities
when the block is nearly lossless, so a fit error `E` appears in it as
roughly `2E`. A block with dissipation of order `1e-5` needs a fit
accurate to well under that before its noise means anything; compare the
fitted dissipation against the data's before trusting the noise of a
fit.
"""
function RationalScattering(block::ScatteringParameters, npoles::Integer; frequencies = nothing,
        tol::Real = 1e-2, atol::Real = 1e-8, fitting::VectorFitting = VectorFitting(),
        passivity::Union{Nothing,PassivityEnforcement} = PassivityEnforcement(), delays = nothing,
        maxstates::Integer = defaultmaxstates(), weights = nothing)
    npoles >= 1 || throw(ArgumentError("fit at least one pole."))
    (isfinite(tol) && tol >= 0) || throw(ArgumentError("tol must be finite and nonnegative."))
    fs, S, options = fitsetup(block, frequencies, delays; atol, fitting, passivity, maxstates, weights)
    (; fit, err, contraction, states) = fitsampled(block, S, fs, Int(npoles); tol, options...)
    states > maxstates && throw(ArgumentError(lazy"the fit at $(npoles) poles has $(states) states, more than maxstates = $(maxstates): its realization and passivity enforcement would hold dense matrices of that many states. Fit with fewer poles, or raise maxstates where the memory allows."))
    isnothing(fit) && throw(ArgumentError(lazy"no fit at $(npoles) poles comes within the tol of $(tol): every fit with the poles the relocation settled on misses the data by at least $(err) of the largest response over the samples. Fit with more poles or over a narrower band, take a delay out with delays, or raise tol."))
    err <= tol || throw(ArgumentError(lazy"the fit at $(npoles) poles misses the data by $(err) of the largest response over the samples, against the tol of $(tol): fit with more poles or over a narrower band, take a delay out with delays, or raise tol to accept a fit that far from the data."))
    warncontraction(contraction)
    return fit
end

# The options both fitting methods take, checked, and the block sampled
# for them with any delay taken out: the frequencies in Hz, the samples,
# and the options of `fitsampled`, the noise the fitted block carries,
# and the components of the samples the relocation runs on and their
# largest response among them, the latter computed once since the order
# search fits the same samples at every order.
function fitsetup(block::ScatteringParameters, frequencies, delays; atol, fitting, passivity,
        maxstates, weights)
    checkatol(atol)
    maxstates >= 1 || throw(ArgumentError("maxstates must be at least one."))
    taus = fitdelays(delays, block.nports)
    fs, S = samplescattering(block, frequencies; delays = taus)
    W = fitweights(weights, size(S))
    options = (; atol, fitting, passivity, noise = undelaynoise(block.noise, taus), weights = W,
        components = relocationcolumns(S, W), scale = max(largestopnorm(weightedsamples(S, W)), floatmin(Float64)),
        maxstates = Int(maxstates))
    return fs, S, options
end

# The weights a fit takes, checked against the samples' shape `dims`:
# none, or a positive finite weight for each entry at each sample. Every
# measure of a weighted fit is relative, so the weights are taken relative
# to the largest, in their own arithmetic, which a common scale can then
# neither underflow nor overflow; weights the same everywhere weigh every
# entry and sample alike, which is the fit without weights.
function fitweights(weights, dims::NTuple{3,Int})
    isnothing(weights) && return noweights
    weights isa AbstractArray{<:Real,3} && size(weights) == dims || throw(ArgumentError(
        lazy"give weights as an array of the samples' shape, one for each entry at each sample: $(dims)."))
    all(w -> isfinite(w) && w > 0, weights) || throw(ArgumentError("the weights must be finite and positive."))
    largest = maximum(weights)
    all(==(largest), weights) && return noweights
    W = Array{Float64,3}(weights ./ largest)
    all(>(0), W) || throw(ArgumentError(
        "the weights span more than Float64 holds: relative to the largest, some are zero in Float64."))
    return W
end

# The states a fit may take by default. Its realization, and the passivity
# enforcement's normal matrix where the samples carry no weights, hold
# dense matrices of the states: together at their peak, with what the
# collector has yet to free, under `statebytes` bytes per squared state.
# The default keeps them within `share` of the memory, as
# `Sys.total_memory` reports it, which a container's limit constrains; a
# weighted metric is held to what the states leave of the same budget (see
# weightedmetricbytes).
const statebytes = 256
defaultmaxstates(; share::Real = 0.5) = floor(Int, sqrt(share*Sys.total_memory()/statebytes))

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
function samplescattering(block::ScatteringParameters, frequencies; delays = nothing)
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
    undelayscattering!(S, ws, delays)
    return fs, S
end

# One finite nonnegative delay per port, or nothing where none is asked
# for and where every one of them is zero, which is the same block.
function fitdelays(delays, nports::Int)
    isnothing(delays) && return nothing
    taus = Float64.(collect(delays))
    length(taus) == nports && all(t -> isfinite(t) && t >= 0, taus) || throw(ArgumentError(
        lazy"give one finite nonnegative delay per port ($(nports))."))
    return all(iszero, taus) ? nothing : taus
end

# The delays taken out of the samples: a wave from port `q` to port `p`
# travels the line of both, so the block left to fit is the device with
# a lossless line of `tau_p` cut off each port.
function undelayscattering!(S, ws, taus)
    isnothing(taus) && return S
    @inbounds for i in eachindex(ws), q in axes(S, 2), p in axes(S, 1)
        S[p, q, i] *= cis(ws[i]*(taus[p] + taus[q]))
    end
    return S
end

# A stated covariance at the same reference planes. It is a correlation
# of the waves the block emits, each delayed by the line of its own
# port, so it carries the difference of the two delays where the scatter
# carries their sum. The phase is put on the covariance the provider
# returns at each frequency, so it composes with the interpolation and
# the extrapolation the covariance declares and holds between the
# samples as well as on them. Nothing to rotate where that difference is
# zero for every pair the covariance is nonzero on: delays which are all
# equal, and a covariance of uncorrelated ports.
function undelaynoise(noise, taus)
    (isnothing(taus) || !(noise isa NoiseCovariance)) && return noise
    all(==(first(taus)), taus) && return noise
    v = noise.provider
    v isa ConstantMatrixProvider && isdiag(v.A) && return noise
    return NoiseCovariance(RotatedMatrixProvider(v, taus), noise.interpolation,
        noise.extrapolation, noise.atol, noise.completed, noise.padding)
end

# One fit of already sampled data at a fixed order, as the fitted block
# `fit`; its error `err` against the samples (see relativefiterror); the
# `contraction` the passivity enforcement applied to it, `nothing` for
# none, which the caller warns of if it returns the fit; the error `raw`
# of the fit before its enforcement, as a fraction of the largest
# response alike; and its `states`. `components` are the columns the
# relocation runs on (see relocationcolumns) and `scale` the samples'
# largest response, weighted where they are.
# A fit of more than `maxstates` states goes no further than the count:
# no fit is returned, and an error of `Inf`. Nor is one returned where
# the poles admit no fit within a finite `tol`, which is then not
# enforced: the error is a lower bound on that of any fit with them.
function fitsampled(block::ScatteringParameters, S, fs, npoles::Int; tol::Real, atol::Real,
        fitting::VectorFitting, passivity::Union{Nothing,PassivityEnforcement},
        components::Matrix{ComplexF64}, scale::Real, maxstates::Int, noise = block.noise,
        weights::Array{Float64,3} = noweights)
    # a block which states what it does at zero frequency has the fit meet
    # it exactly; one which states nothing leaves it to the extrapolation
    dc = block.dcmodel isa ScatteringLimit ? nodc : dcscatteringmatrix(block.dcmodel, block.nports)
    # the constant term is repaired only for a caller who asked for the
    # fit to be made passive, and to that caller's tolerances; a block
    # which states its noise may be active and is fitted as it is, and
    # neither enforced nor tested
    tested = !(noise isa NoiseCovariance)
    enforce = tested && !isnothing(passivity)
    poles, residues, D = vectorfit(S, 2pi .* fs, npoles, fitting; dc = dc,
        constanttol = enforce ? atol : nothing, constantmargin = enforce ? passivity.margin : 0.0,
        weights = weights, components = components, scale = scale)
    raw = first(residualerrors(S, 2pi .* fs, poles, residues, D; weights = weights))/scale
    # the rank of each residue is decided here, once: the enforcement
    # changes each within its space and the realization spans the same
    # spaces, so the block has these states
    spaces = residuespaces(poles, residues, fitting.ranktol)
    states = statecount(spaces)
    states > maxstates && return (; fit = nothing, err = Inf, contraction = nothing, raw, states)
    contraction = nothing
    if enforce
        # The enforcement changes the residues and the constant and keeps
        # the poles, so the fit it returns comes no closer to the samples
        # than the closest fit with these poles; where even that misses a
        # finite `tol`, the enforcement, the costly part of a fit, is not
        # run. The raw fit is one fit with these poles, so where it meets
        # `tol` no bound on them can exceed it.
        if isfinite(tol) && raw > tol
            bound = isempty(weights) ? deviationfloor(components, 2pi .* fs, poles, size(S, 1), tol*scale) :
                deviationfloor(S, weights, 2pi .* fs, poles, tol*scale)
            bound > tol*scale && return (; fit = nothing, err = bound/scale, contraction, raw, states)
        end
        # samples reciprocal to the block's tolerance have a reciprocal fit,
        # which the correction keeps so far as its spaces allow (see
        # PassivityEnforcement)
        reciprocal = maximum(k -> opnorm(view(S, :, :, k) .- transpose(view(S, :, :, k))),
            axes(S, 3); init = 0.0) <= atol
        residues, D, contraction = enforcepassivity(poles, residues, D, 2pi .* fs, passivity;
            atol = atol, dc = dc, spaces = spaces, reciprocal = reciprocal, weights = weights,
            memory = statebytes*(Float64(maxstates)^2 - Float64(states)^2))
    elseif tested
        checkrawpassive(poles, residues, D, 2pi .* fs, spaces; atol = atol)
    end
    A, B, C = realization(poles, residues, block.nports; spaces = spaces)
    # the fit's passivity is tested on its residues, by the enforcement or
    # checkrawpassive, which the validation of the realization would repeat
    fit = rationalblock(A, B, C, D; zref = block.zref, grounded = block.grounded,
        noise = noise, atol = atol, normtested = true)
    return (; fit, err = relativefiterror(fit, S, fs; scale = scale, weights = weights), contraction, raw, states)
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

# The fewest poles a fit of the samples `S` could meet `tol` with, from
# the singular values of the samples, the norms of their `components`
# (see relocationcomponents). A real common pole fit of `N` poles and a
# constant is a combination of `N + 1` real functions of frequency, so
# the samples of every entry, stacked as real and imaginary parts, are
# fitted with an error of Frobenius norm at least that of their best
# approximation of rank `N + 1`, the root sum of squares of the singular
# values beyond it (Eckart and Young). Over `K` samples of `n` ports the
# largest spectral norm of the error at a sample is at least that over
# `sqrt(K n)`, against `tol` of the largest response, `scale`.
function fewestpoles(components::AbstractMatrix, S::AbstractArray{<:Complex,3}, tol::Real;
        scale::Real = largestopnorm(S))
    n, K = size(S, 1), size(S, 3)
    squares = [sum(abs2, view(components, :, j)) for j in axes(components, 2)]
    tail = reverse!(cumsum(reverse(squares)))
    for N in 0:length(squares)
        beyond = N + 2 <= length(tail) ? tail[N + 2] : 0.0
        sqrt(beyond/(K*n)) <= tol*scale && return N
    end
    return length(squares)
end

# A lower bound on the largest spectral norm of the deviation from the
# samples of any fit with the poles `poles`, whatever its residues and
# constant. For weights `w_k` over the samples summing to one, the
# largest `|E_k|^2` is at least `sum_k w_k |E_k|_F^2 / n`, and that at
# least the least weighted sum of squares over the poles' basis, taken
# here on the samples' `components` (see relocationcomponents), whose
# squares are those of the entries less what their truncation drops.
# Lawson's reweighting, each `w_k` in proportion to `w_k |E_k|_F`, raises
# the bound toward the least largest Frobenius norm over `sqrt(n)`, which
# the largest `|E_k|_F / sqrt(n)` of each step's fit bounds from above;
# the steps end once the bound exceeds `target` or that falls to it,
# which settles all a caller asks, or after `steps`.
function deviationfloor(components::AbstractMatrix{<:Complex}, ws::AbstractVector,
        poles::Vector{ComplexF64}, n::Integer, target::Real; steps::Integer = 8)
    K, r = size(components)
    Phi = realbasis(ws, poles)
    M = [real.(Phi) ones(K); imag.(Phi) zeros(K)]
    X = [real.(components); imag.(components)]
    R = similar(X)
    w, e = fill(1/K, K), zeros(K)
    bound = 0.0
    for _ in 1:steps
        sw = sqrt.(vcat(w, w))
        F = qr!(sw .* M)
        # the weighted samples less their projection on the weighted
        # basis, which the reflectors span with any column it lacks
        R .= sw .* X
        lmul!(F.Q', R)
        R[1:size(M, 2), :] .= 0
        lmul!(F.Q, R)
        bound = max(bound, sqrt(sum(abs2, R)/n))
        bound > target && break
        e .= 0
        @inbounds for j in 1:r, k in 1:K
            e[k] += abs2(R[k, j]) + abs2(R[K + k, j])
        end
        for k in 1:K
            e[k] = w[k] > 0 ? sqrt(e[k]/w[k]) : 0.0
        end
        maximum(e)/sqrt(n) <= target && break
        total = dot(w, e)
        total > 0 || break
        w .*= e ./ total
    end
    return bound
end

# The same bound where the samples carry weights, each entry's deviation
# at each sample weighed by its weight: each entry then has a least
# squares of its own (see fitcoefficients), its rows weighed by its
# weights and by the square roots of the Lawson weights `w_k`, and a
# sample's squared error is the sum over its weighted entries.
function deviationfloor(S::AbstractArray{<:Complex,3}, weights::Array{Float64,3}, ws::AbstractVector,
        poles::Vector{ComplexF64}, target::Real; steps::Integer = 8)
    n, K = size(S, 1), size(S, 3)
    Phi = realbasis(ws, poles)
    M = [real.(Phi) ones(K); imag.(Phi) zeros(K)]
    Sf, Wf = reshape(S, n*n, K), reshape(weights, n*n, K)
    w, e = fill(1/K, K), zeros(K)
    rw, r = zeros(2K), zeros(2K)
    bound = 0.0
    for _ in 1:steps
        e .= 0
        for f in 1:n*n
            @inbounds for k in 1:K
                rw[k] = rw[K + k] = sqrt(w[k])*Wf[f, k]
                r[k], r[K + k] = rw[k]*real(Sf[f, k]), rw[K + k]*imag(Sf[f, k])
            end
            # the entry's weighted samples less their projection on its
            # weighted basis
            F = qr!(rw .* M)
            lmul!(F.Q', r)
            r[1:size(M, 2)] .= 0
            lmul!(F.Q, r)
            @inbounds for k in 1:K
                e[k] += abs2(r[k]) + abs2(r[K + k])
            end
        end
        bound = max(bound, sqrt(sum(e)/n))
        bound > target && break
        for k in 1:K
            e[k] = w[k] > 0 ? sqrt(e[k]/w[k]) : 0.0
        end
        maximum(e)/sqrt(n) <= target && break
        total = dot(w, e)
        total > 0 || break
        w .*= e ./ total
    end
    return bound
end

# How far a fit is from its data, as a fraction of the largest
# response: the worst sample, not the mean, since a fit which is
# excellent almost everywhere and wrong at one resonance is not a model
# of the block. `scale` is the samples' largest response where the caller
# has it.
function relativefiterror(fit, S, fs; scale::Union{Nothing,Real} = nothing,
        weights::Array{Float64,3} = noweights)
    worst, largest = fitdeviation(fit.provider, 2pi .* fs, S; scale = scale, weights = weights)
    return worst/max(largest, floatmin(Float64))
end
# the largest deviation of a provider from the samples `S` at the angular
# frequencies `ws`, in the spectral norm, and the largest sample, `scale`
# where it is given; each entry at each sample weighed by its weight
# where `weights` are given
function fitdeviation(provider, ws, S; scale::Union{Nothing,Real} = nothing,
        weights::Array{Float64,3} = noweights)
    F = similar(S, ComplexF64)
    evaluateprovider!(F, provider, ws)
    F .-= S
    isempty(weights) || (F .*= weights)
    return largestopnorm(F), isnothing(scale) ? largestopnorm(weightedsamples(S, weights)) : scale
end

# The largest spectral norm of the matrices `M[:, :, k]`, with a singular
# value decomposition only where one could matter: a spectral norm is at
# most the Frobenius norm, so the matrices taken in decreasing Frobenius
# norm stop once the next one's is under the largest spectral norm found;
# a matrix the first screen leaves is passed over where cheaper bounds put
# it under that too (boundedby!).
function largestopnorm(M::AbstractArray{<:Number,3})
    fro = [norm(view(M, :, :, k)) for k in axes(M, 3)]
    worst = 0.0
    rows, W = zeros(size(M, 1)), Matrix{float(eltype(M))}(undef, size(M, 2), size(M, 2))
    for k in sortperm(fro; rev = true)
        fro[k] <= worst && break
        boundedby!(W, rows, view(M, :, :, k), worst) && continue
        worst = max(worst, opnorm(view(M, :, :, k)))
    end
    return worst
end

# `sqrt(|A|_1 |A|_inf)`, an upper bound on the spectral norm, from the
# largest column sum and the largest row sum, gathered into `rows`, in one
# pass over `A`
function holderbound(A::AbstractMatrix, rows::Vector{Float64})
    fill!(rows, 0.0)
    columns = 0.0
    @inbounds for j in axes(A, 2)
        c = 0.0
        for i in axes(A, 1)
            a = abs(A[i, j])
            c += a
            rows[i] += a
        end
        columns = max(columns, c)
    end
    return sqrt(columns)*sqrt(maximum(rows))
end

# Whether the largest singular value of `S` exceeds `level`: whether
# `level^2 I - S'S` fails to be positive definite, which a Cholesky
# factorization decides at a fraction of the cost of a singular value
# decomposition. `W` is its scratch.
function exceeds!(W::AbstractMatrix, S::AbstractMatrix, level::Real)
    mul!(W, S', S, -1, 0)
    for i in axes(W, 1)
        W[i, i] += level^2
    end
    return !issuccess(cholesky!(Hermitian(W); check = false))
end

# Whether the spectral norm of `A` is under `level` by a bound cheaper
# than a decomposition: `sqrt(|A|_1 |A|_inf)` (holderbound), then the
# Cholesky test (exceeds!) at `level` less the `n^2 eps` of it by which
# forming `A'A` and factoring it can err. `W` and `rows` are scratch.
function boundedby!(W::AbstractMatrix, rows::Vector{Float64}, A::AbstractMatrix, level::Real)
    holderbound(A, rows) <= level && return true
    return level > 0 && !exceeds!(W, A, level*(1 - 2*length(A)*eps()))
end

# The fraction of itself by which the fit's error before its enforcement
# must fall for the order search to count it as progress (see the search's
# stop on progress)
const searchprogress = 0.1

"""
    RationalScattering(block::ScatteringParameters; tol, minpoles = 4,
        maxpoles = nothing, noisefloor = 1e-12, frequencies = nothing,
        atol = 1e-8, fitting = VectorFitting(),
        passivity = PassivityEnforcement(), delays = nothing,
        maxstates = defaultmaxstates(), weights = nothing)

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
way not to pay for orders a block is known not to need. The samples
themselves rule out the lowest orders of a block with many ports: a fit
of `N` poles and a constant is a combination of `N + 1` real functions
of frequency, so its error is at least the part of the samples beyond
their best approximation of that rank (Eckart and Young), and orders
whose bound exceeds `tol` are not fitted. Nor is an order up to the
degree the samples determine made passive whose poles admit no fit
within `tol`: the enforcement keeps the poles, and a least squares at
them, weighted toward the samples it misses most, bounds the error of
every fit with them from below. Past that degree every order is made
passive, since which of them can be decides where the scan ends.

Without `maxpoles` the search budgets itself by the degree the samples
determine, the numerical rank of their Loewner pencil with `noisefloor`
as the rank threshold, as a fraction of the largest singular value; it
scans to twice that degree or to one pole fewer than the samples,
whichever is fewer, and stops sooner if four consecutive orders past
that degree produce no fit at all, which is the model class running out.
Below the degree an order which produces no fit has too few poles for
the data, and the scan goes on. The degree is an estimate, not a bound:
it moves with `noisefloor`, which measured data with a real noise floor
wants larger than the default, which suits data good to nearly full
precision; the pencil is built along cycling coordinate directions from
at most four hundred samples, so dynamics weak or narrow in the
directions it does not probe can be missed; and a constant term
contributes to it. Data of a high degree, as measured interconnects with
long lines are, can put the end of the scan hundreds of orders away,
each dearer than the last, so the scan also stops once the order has
doubled since the fit's own error, before its passivity enforcement,
last fell by a tenth of itself, counted from the first fit whose
error is under the largest response. An error that has stalled can still
fall again at much higher orders, so this is a budget rather than a
proof, and the refusal says where the scan stopped and why. `maxpoles`
overrides all of this and is a hard ceiling, which the scan reaches
unless an order meets `tol`.

Whatever the ceiling, a fit of more than `maxstates` states ends the
scan, since the memory a fit takes grows as the square of its states
and the states grow with the order (see the `npoles` method).

If no order meets `tol`, the closest fit found and its order, the
nearest of the orders refused before the enforcement and its bound, the
degree the samples determine, where the scan ended if it stopped on its
progress or on `maxstates`, and why any orders failed to fit are
reported as an error, since a block quietly less accurate than asked for
is worse than none: loosen `tol`, raise `maxpoles`, sample the block more
finely, or fit a narrower band. Where no order produced a fit at all,
the error says so, and names why they failed.

A `delays` given here is taken out before the search, so the order it
reports is the order of what is left after the delay, which for a cable
is a small fraction of what the delay itself would cost.

A contraction made to bring a fit under one is warned of for the fit
returned alone: an order the search discards describes nothing the
caller receives.

Every order starts from the spacing `fitting` names; starting poles
given as a vector fix one order, and are refused here.

With `weights` every error is the weighted one (see the `npoles`
method), and the samples' singular values rule out no order, each entry
being fitted by the basis weighed by its own weights; the least squares
bound on the poles of an order still refuses it before its enforcement.

See the `npoles` method for the meaning of the remaining arguments, and
for the warning about fitting a block with little loss: a tolerance
which looks tight against `S` may still be far too loose against the
dissipation `I - S S'`, where the error appears roughly doubled.
"""
function RationalScattering(block::ScatteringParameters; tol::Real,
        minpoles::Integer = 4, maxpoles = nothing, noisefloor::Real = 1e-12,
        frequencies = nothing, atol::Real = 1e-8, fitting::VectorFitting = VectorFitting(),
        passivity::Union{Nothing,PassivityEnforcement} = PassivityEnforcement(), delays = nothing,
        maxstates::Integer = defaultmaxstates(), weights = nothing)
    (isfinite(tol) && tol > 0) || throw(ArgumentError("tol must be finite and positive."))
    minpoles >= 1 || throw(ArgumentError("give minpoles >= 1."))
    (isfinite(noisefloor) && noisefloor > 0) || throw(ArgumentError("noisefloor must be finite and positive."))
    fitting.start isa Symbol || throw(ArgumentError("the search fits many orders and starting poles given as a vector fix one: fit at their order with RationalScattering(block, npoles), or start the search from a spacing, :linear, :log or :linlog."))
    fs, S, options = fitsetup(block, frequencies, delays; atol, fitting, passivity, maxstates, weights)
    # fitting N poles needs N + 1 samples, so a scan from minpoles needs
    # at least that many; refusing here names the samples rather than
    # reporting a scan of no orders
    length(fs) >= Int(minpoles) + 1 || throw(ArgumentError(
        lazy"fitting $(minpoles) poles needs at least $(Int(minpoles) + 1) sample frequencies; there are $(length(fs))."))
    supported = supporteddegree(S, 2pi .* fs, Float64(noisefloor))
    # weighted, each entry is fitted by the basis weighed by its own
    # weights, so no bound on the rank of their fit rules out an order
    floorpoles = isempty(options.weights) ? fewestpoles(options.components, S, tol; scale = options.scale) : 0
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
    # the orders refused before their enforcement, and the least of the
    # bounds that refused them and its order
    ruledout, leastbound, boundpoles = Int[], Inf, 0
    # the fit's own error before its enforcement where it last fell by
    # `searchprogress` of itself, and that order
    gained, gainedpoles = Inf, 0
    attempt = np -> begin
        (; fit, err, contraction, raw, states) = try
            # past the ceiling every order is enforced: the end of the
            # scan there is decided by which orders produce a fit at all
            fitsampled(block, S, fs, np; tol = np <= ceiling ? tol : Inf, options...)
        catch e
            e isa ArgumentError || rethrow()
            # the message without the orders and numbers in it, so that
            # the same failure at twenty orders is one line and not twenty
            key = replace(sprint(showerror, e), r"[-+]?[0-9][0-9.e+-]*" => "N")
            push!(get!(failures, key, Int[]), np)
            (; fit = nothing, err = Inf, contraction = nothing, raw = Inf, states = 0)
        end
        # a fit whose error is the largest response or more has caught
        # nothing of the data, which zero does as well
        if raw < 1 && raw < (1 - searchprogress)*gained
            gained, gainedpoles = raw, np
        end
        if isnothing(fit) && isfinite(err)
            push!(ruledout, np)
            err < leastbound && ((leastbound, boundpoles) = (err, np))
        else
            err < besterror && ((besterror, bestpoles) = (err, np))
        end
        return (fit, err, contraction, states)
    end
    # One order at a time from `minpoles` upward, returning the first
    # that meets the tolerance, which is therefore the fewest poles that
    # do. There is nothing to bisect on: more poles do not always fit
    # better, since past the order the data supports the relocation is
    # decided by directions the samples do not determine, so the orders
    # meeting a tolerance are a window and not a tail. An order stepped
    # over is an order about which nothing is known, so every order
    # below the answer has to be fitted, and `minpoles` is the way not
    # to pay for orders a block is known not to need.
    #
    # Without `maxpoles` the degree estimate is a budget and not a wall:
    # it can be short of what the block needs, so the scan continues to
    # twice the estimate, and never to the sample count, since `N` poles
    # take `N + 1` samples. A caller's `maxpoles` is a wall.
    stop = isnothing(maxpoles) ? min(2*ceiling, length(fs) - 1) : ceiling
    floorpoles <= stop || throw(ArgumentError(lazy"no fit of at most $(stop) poles can meet a tolerance of $(tol): the samples' own singular values put the error of any fit of that many poles above it. Loosen tol, raise maxpoles, or fit a narrower band."))
    np, barren, stalled, overstates = max(Int(minpoles), floorpoles), 0, false, 0
    while np <= stop
        fit, err, contraction, states = attempt(np)
        # only the fit returned warns of its contraction
        err <= tol && (warncontraction(contraction); return fit)
        # the states grow with the order, and the memory as their square
        if states > maxstates
            overstates = states
            break
        end
        # A run of orders producing no fit at all past the degree the
        # samples determine is the model class running out, and that
        # does not reverse: there the relocation is decided by
        # directions the samples do not carry and the result cannot be
        # made passive, at a cost that rises with the order. Below the
        # degree an order producing no fit has too few poles, which the
        # next orders remedy, so only the orders past it count. Stopping
        # needs a run rather than one refusal, which can be a numerical
        # accident, an order which cannot be made passive beside one
        # which fits to the roundoff. (The orders meeting a
        # tolerance are a window, so orders of no improvement establish
        # nothing about the next, and the stop on progress below is a
        # budget, not a proof.)
        barren = (isfinite(err) || np <= ceiling) ? 0 : barren + 1
        barren >= 4 && break
        # Without `maxpoles` the scan also ends on measured progress: once
        # the order has doubled since the fit's own error, before its
        # enforcement, last fell by `searchprogress` of itself. An error that
        # stalls can still fall again at much higher orders, so this is a
        # budget rather than a proof; the refusal says where it stopped,
        # and a `maxpoles` scans on to it.
        if isnothing(maxpoles) && gainedpoles > 0 && np >= 2*gainedpoles
            stalled = true
            break
        end
        np += 1
    end
    reached = min(np, stop)
    why = join(("$(length(v)) of them ($(first(v)) to $(last(v)) poles) with \"$(k)\""
        for (k, v) in sort(collect(failures); by = x -> -length(x[2]))), "; ")
    ended = overstates > 0 ? " The scan ended at $(np) poles, whose fit has $(overstates) states, more than maxstates = $(maxstates); a larger maxstates scans on where the memory allows." :
        stalled ? " The scan ended at $(np) poles, twice the order at which the fit's own error before its passivity enforcement last fell by $(searchprogress) of itself, to $(gained) at $(gainedpoles) poles; a maxpoles scans on to it." : ""
    # with nothing fitted, refused or failed, the first order's states
    # ended the scan
    (isfinite(besterror) || !isempty(ruledout) || !isempty(failures)) || throw(ArgumentError(
        lazy"the fit at $(np) poles, where the scan starts, has $(overstates) states, more than maxstates = $(maxstates). Raise maxstates where the memory allows, loosen tol, or fit a narrower band."))
    (isfinite(besterror) || !isempty(ruledout)) || throw(ArgumentError(lazy"no order between $(minpoles) and $(reached) poles could be fitted at all, while the samples determine a degree of about $(supported): $(why).$(ended) Raise maxpoles, sample the block more finely, or fit a narrower band."))
    failed = isempty(failures) ? "" : " Some orders could not be fitted at all: $(why)."
    closest = isfinite(besterror) ? "the closest was $(besterror) at $(bestpoles) poles" :
        "no order came near enough to be made passive"
    refused = isempty(ruledout) ? "" : " $(length(ruledout)) orders ($(first(ruledout)) to $(last(ruledout)) poles) were refused before their passivity enforcement: no fit with the poles they settled on comes within the tolerance, and the nearest, at $(boundpoles) poles, misses the data by at least $(leastbound)."
    throw(ArgumentError(lazy"no fit between $(minpoles) and $(reached) poles met a tolerance of $(tol): $(closest), while the samples determine a degree of about $(supported).$(ended)$(refused)$(failed) Loosen tol, raise maxpoles, sample the block more finely, or fit a narrower band."))
end

# No value stated at zero frequency, which the fit's functions take as an
# empty matrix rather than `nothing`: a statement is then one type, and
# the functions are compiled once for fits with and without one.
const nodc = zeros(0, 0)
# No weights on the samples, taken the same way.
const noweights = zeros(0, 0, 0)

# The relaxed vector fit: with the poles `a`, the unknowns of every entry
# are its residues and constant, and the shared unknowns the residues
# and the constant of the weight `sigma(s) = d + sum r_p/(s - a_p)` whose
# zeros are the next poles, with the relaxed normalization that the real
# part of `sigma` averages to one over the samples so that the trivial
# solution is excluded; complex poles come in conjugate pairs and enter
# through the real basis, so every unknown is real. Returns the poles,
# the residue matrices `(n, n, npoles)` and the constant `(n, n)`.
# `lastpole` is what becomes of the last pole when a constant reproduces
# the data: `:refuse` throws, since the data is then no rational block,
# `:drop` returns no poles and the constant, and `:keep` keeps the pole,
# for a strictly proper fit which cannot be a constant. `proper` holds
# the constant at zero in every least squares, the relocation's and the
# pruning's as well as the residues', for a fit which vanishes at infinite
# frequency. `weights`, where the samples carry them, weigh each entry at
# each sample in every least squares and every measure of the fit (see
# RationalScattering).
# `components` are those of `S` the relocation runs on (see
# relocationcolumns), and `scale` its largest response, weighted where the
# samples are, which the errors' roundoff is measured against, given by a
# caller which fits the same samples more than once.
function vectorfit(S::AbstractArray{<:Complex,3}, ws::AbstractVector, npoles::Int,
        fitting::VectorFitting; dc::AbstractMatrix = nodc,
        constanttol::Union{Nothing,Real} = nothing, constantmargin::Real = 1e-6,
        proper::Bool = false, lastpole::Symbol = :refuse, weights::Array{Float64,3} = noweights,
        components::Matrix{ComplexF64} = relocationcolumns(S, weights),
        scale::Real = largestopnorm(weightedsamples(S, weights)))
    lastpole in (:refuse, :drop, :keep) || throw(ArgumentError("lastpole is :refuse, :drop or :keep."))
    # a real rational function is real at zero frequency, so a complex
    # statement there cannot be met by any fit
    maximum(abs, imag.(dc); init = 0.0) == 0 || throw(ArgumentError(
        "the value stated at zero frequency must be real: a real rational function is real at zero."))
    dc = Matrix{Float64}(real.(dc))
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
    start = fitting.start
    poles = if start isa Symbol
        spreadpoles(start, wmin, wmax, npoles)
    else
        length(start) == npoles || throw(ArgumentError(lazy"$(length(start)) starting poles are given for a fit at $(npoles) poles; give as many as the poles fitted."))
        start ./ wref
    end
    poles = converge(S, xs, poles, fitting; dc = dc, weights = weights, components = components, scale = scale,
        proper = proper)
    # a relocation which has lost the data leaves poles which are not
    # finite; fail here, naming the fit, rather than in a factorization
    # built from them
    all(isfinite, poles) || throw(ArgumentError(lazy"the pole relocation diverged at $(npoles) poles: the relocated poles are not finite. Fit with fewer poles, or over a narrower band."))
    poles = prunepoles(S, xs, poles, fitting; dc = dc, lastpole = lastpole, weights = weights,
        components = components, scale = scale, proper = proper)
    # the poles are settled from the data before the value at zero is
    # stated, so a condition outside the band never moves them
    residues, D = fitresidues(S, xs, poles; dc = dc, weights = weights,
        constant = proper ? zeros(n, n) : nothing)
    (all(isfinite, residues) && all(isfinite, D)) || throw(ArgumentError(lazy"the residues at $(length(poles)) poles are not finite: the least squares of the fit is singular. Fit with fewer poles, or over a narrower band."))
    # A constant term above one is brought under here, where the residues
    # can be refitted around it over the samples, rather than by the
    # enforcement's corrections, which move it with the residues by the
    # least change over the samples -- and only when the caller asked for
    # passivity: samples of `1.5 - 1.2/(s + 1)` have an exact one pole fit
    # whose constant is 1.5, and a raw fit must return it.
    if !isnothing(constanttol)
        D, clamped = passiveconstant(D; tol = Float64(constanttol), margin = Float64(constantmargin))
        if clamped
            residues, D = fitresidues(S, xs, poles; constant = D, dc = dc, weights = weights)
            all(isfinite, residues) || throw(ArgumentError(lazy"the residues at $(length(poles)) poles are not finite once the constant term is made passive: the least squares of the fit is singular. Fit with fewer poles, or over a narrower band."))
        end
    end
    return poles .* wref, residues .* wref, D
end

# the starting poles spread over the positive band `[wmin, wmax]`:
# complex pairs evenly in frequency (`:linear`) or in its logarithm
# (`:log`), or `cld(npoles - 1, 4)` of them evenly in frequency and
# `npoles ÷ 4` evenly in its logarithm without the band's ends (`:linlog`,
# the linlogcmplx start of VFdriver, which spreads five or fewer as
# `:log`), damped by `damping` of their frequency, and a real pole at the
# centre of the band, in the logarithm but for `:linear`, if the count is
# odd
function spreadpoles(spacing::Symbol, wmin::Real, wmax::Real, npoles::Int; damping::Real = 0.01)
    npairs = npoles ÷ 2
    spacing == :linlog && npoles < 6 && (spacing = :log)
    centre = spacing == :linear ? (wmin + wmax)/2 : sqrt(wmin*wmax)
    spread = npairs == 1 ? [centre] : spacing == :log ?
        exp.(range(log(wmin), log(wmax); length = npairs)) : spacing == :linlog ?
        sort!(vcat(range(wmin, wmax; length = cld(npoles - 1, 4)),
            exp.(range(log(wmin), log(wmax); length = npoles ÷ 4 + 2))[2:end - 1])) :
        collect(range(wmin, wmax; length = npairs))
    poles = ComplexF64[]
    for w in spread
        push!(poles, complex(-damping*w, w))
        push!(poles, complex(-damping*w, -w))
    end
    isodd(npoles) && push!(poles, complex(-centre, 0.0))
    return poles
end

# the scale an error is measured against, in the norm it is measured in,
# as a fraction of the samples' largest response
const fitroundoff = 1e-12
# The samples as the fit's measures read them: each entry at each sample
# times its weight, where the samples carry weights (see
# RationalScattering).
weightedsamples(S::AbstractArray{<:Complex,3}, weights::Array{Float64,3}) = isempty(weights) ? S : weights .* S
# The columns the relocation runs on: the samples' components (see
# relocationcomponents), or where the samples carry weights, which weigh
# each entry's rows by its own, the entries, a column each (see relocate).
relocationcolumns(S::AbstractArray{<:Complex,3}, weights::Array{Float64,3}) =
    isempty(weights) ? relocationcomponents(S) : permutedims(reshape(S, :, size(S, 3)))
# The relocation iterated from `poles`, each iterate measured on the
# samples `S` and relocated on their `components` (see relocationcolumns),
# which relocate them the same, until its error settles or reaches the
# roundoff of the samples, where the weight's columns are combinations of
# the entries' and the relocation would be arbitrary. `scale` is the
# samples' largest response, which a caller that has it gives rather than
# have it computed again.
function converge(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64},
        fitting::VectorFitting; dc::Matrix{Float64} = nodc, weights::Array{Float64,3} = noweights,
        components::Matrix{ComplexF64} = relocationcolumns(S, weights),
        scale::Real = largestopnorm(weightedsamples(S, weights)), proper::Bool = false,
        errorsettletol::Real = 1e-3)
    (; iterations, stallpatience) = fitting
    exact = fitroundoff*scale
    # The iterate kept is the one measured to fit best, not the one
    # whose poles last moved least: a relocation which is circling can
    # pass its best fit on an iteration where the poles are moving
    # quickly. Every iterate is measured, at the cost of one residue
    # solve each.
    best, besterror = copy(poles), Inf
    # The rounds since the best error last fell by errorsettletol of
    # itself, which end the relocation at stallpatience: a fit which has
    # settled while its poles still creep, or one circling. How far the
    # poles move is no measure of it, since on a wide band the poles near
    # its low end move by their own size from round to round while the
    # fit still improves.
    since = 0
    for _ in 1:iterations
        err = fiterror(S, ws, poles; dc = dc, weights = weights, proper = proper)
        since = err < (1 - errorsettletol)*besterror ? 0 : since + 1
        err < besterror && ((best, besterror) = (copy(poles), err))
        # a fit which is already exact settles its poles too
        err <= exact && return best
        since >= stallpatience && return best
        poles = relocate(components, ws, poles; weights = weights, proper = proper)
    end
    # out of iterations, the iterate that is returned is the one measured
    # to fit best, the last one measured as well
    err = fiterror(S, ws, poles; dc = dc, weights = weights, proper = proper)
    return err < besterror ? poles : best
end

# The deviation of the fit at the given poles from the data over the
# samples, with the residues refitted, in two measures from one residue
# solve: the largest over the samples in the spectral norm, which the
# relocation, the pruning and the acceptance are judged on, and the rms
# over every entry and sample, which the pruning holds as well. The
# coefficients are the residue solve's (see fitcoefficients), so what is
# measured is the fit that solve returns: its samples, `M x` for an
# entry's coefficients `x`, are formed `tile` samples at a time as one
# product into a buffer held over the tiles, a sample's entries
# contiguous, and the deviation at the samples that could hold the
# largest spectral norm is formed again from the coefficients. Where the
# samples carry `weights`, each entry's deviation at each sample is
# weighed by its weight, in both measures and in the solve. A `proper` fit
# holds its constant at zero.
function fiterrors(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64};
        dc::Matrix{Float64} = nodc, tile::Int = 64, weights::Array{Float64,3} = noweights, proper::Bool = false)
    n, K = size(S, 1), length(ws)
    X, M = fitcoefficients(S, ws, poles; dc = dc, weights = weights, constant = proper ? zeros(n, n) : nothing)
    Sf = reshape(S, n*n, K)
    weighted = !isempty(weights)
    Wf = weighted ? reshape(weights, n*n, K) : zeros(0, 0)
    squares = zeros(K)
    FR, FI = Matrix{Float64}(undef, n*n, min(tile, K)), Matrix{Float64}(undef, n*n, min(tile, K))
    for lo in 1:tile:K
        ids = lo:min(lo + tile - 1, K)
        R, J = view(FR, :, 1:length(ids)), view(FI, :, 1:length(ids))
        mul!(R, X, transpose(view(M, ids, :)))
        mul!(J, X, transpose(view(M, K .+ ids, :)))
        @inbounds for (j, k) in enumerate(ids)
            acc = 0.0
            for e in 1:n*n
                d = complex(R[e, j], J[e, j]) - Sf[e, k]
                acc += abs2(weighted ? Wf[e, k]*d : d)
            end
            squares[k] = acc
        end
    end
    # The largest spectral norm, the samples taken in decreasing Frobenius
    # norm until the next one's is under it (see largestopnorm). A sample
    # is measured by products of its own, which round differently from
    # its tile's: two evaluations of the same inner products of `N + 1`
    # terms differ by at most `2 gamma |X| |m|` an entry, so a sample's
    # Frobenius norm is bounded by its tile's plus
    # `2 gamma |X|_F |[m_k; m_(K+k)]|`. That is nothing in a well
    # conditioned basis, and in one whose large coefficients cancel it is
    # as large as the error the rounding leaves unresolved, where the
    # screen passes over nothing.
    u = size(X, 2)*eps()
    reach = 2*u/(1 - u)*norm(X)
    # weighted, by the largest weight at the sample
    bounds = [sqrt(squares[k]) + reach*(weighted ? maximum(view(Wf, :, k)) : 1.0)*
        sqrt(sum(abs2, view(M, k, :)) + sum(abs2, view(M, K + k, :))) for k in 1:K]
    worst = 0.0
    E, W = Matrix{ComplexF64}(undef, n, n), Matrix{ComplexF64}(undef, n, n)
    er, ei, rows = zeros(n*n), zeros(n*n), zeros(n)
    for k in sortperm(bounds; rev = true)
        bounds[k] <= worst && break
        mul!(er, X, view(M, k, :))
        mul!(ei, X, view(M, K + k, :))
        @inbounds for e in 1:n*n
            E[e] = complex(er[e], ei[e]) - Sf[e, k]
            weighted && (E[e] *= Wf[e, k])
        end
        boundedby!(W, rows, E, worst) && continue
        worst = max(worst, opnorm(E))
    end
    return worst, sqrt(sum(squares)/length(S))
end

# the largest spectral norm and the rms of the deviation of the residues
# and constant `residues`, `D` at the poles `poles` from the samples,
# weighed by the samples' `weights` where they carry them
function residualerrors(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64},
        residues::AbstractArray{<:Complex,3}, D::AbstractMatrix; weights::Array{Float64,3} = noweights)
    n = size(D, 1)
    # accumulated in place, since a sum over a generator of matrices
    # would allocate one per pole per sample, and this runs once per
    # relocation iteration and once per pruning candidate
    F = Matrix{ComplexF64}(undef, n, n)
    squares = Vector{Float64}(undef, length(ws))
    total = 0.0
    for k in eachindex(ws)
        squares[k] = sum(abs2, fitresidual!(F, S, D, residues, poles, ws, k, weights))
        total += squares[k]
    end
    # the largest spectral norm, the samples taken in decreasing Frobenius
    # norm until the next one's is under it (see largestopnorm)
    worst = 0.0
    rows, W = zeros(n), Matrix{ComplexF64}(undef, n, n)
    for k in sortperm(squares; rev = true)
        sqrt(squares[k]) <= worst && break
        fitresidual!(F, S, D, residues, poles, ws, k, weights)
        boundedby!(W, rows, F, worst) && continue
        worst = max(worst, opnorm(F))
    end
    return worst, sqrt(total/length(S))
end

# `F` the fit less the samples at the `k`th sample, each entry weighed by
# its weight there where `weights` are given
function fitresidual!(F::Matrix{ComplexF64}, S::AbstractArray{<:Complex,3}, D::AbstractMatrix,
        residues::AbstractArray{<:Complex,3}, poles::Vector{ComplexF64}, ws::AbstractVector, k::Int,
        weights::Array{Float64,3} = noweights)
    n, w = size(D, 1), ws[k]
    @inbounds for j in 1:n, i in 1:n
        F[i, j] = D[i, j] - S[i, j, k]
    end
    # a pole pair's residues `r` and `conj(r)` add
    # `Re r (u + v) + Im r i (u - v)`, its two functions formed without
    # their difference (see pairbasis)
    p = 1
    @inbounds while p <= length(poles)
        a = poles[p]
        if isrealpole(a)
            c = 1/(im*w - real(a))
            for j in 1:n, i in 1:n
                F[i, j] += residues[i, j, p]*c
            end
            p += 1
        else
            f1, f2 = pairbasis(a, im*w)
            for j in 1:n, i in 1:n
                r = residues[i, j, p]
                F[i, j] += real(r)*f1 + imag(r)*f2
            end
            p += 2
        end
    end
    isempty(weights) || (F .*= view(weights, :, :, k))
    return F
end
# the largest deviation alone
fiterror(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64};
    dc::Matrix{Float64} = nodc, weights::Array{Float64,3} = noweights, proper::Bool = false) =
    first(fiterrors(S, ws, poles; dc = dc, weights = weights, proper = proper))
# The poles the fit does not need, once the relocation has settled: a
# needless pole drifts out of the band, above it fitting a constant the
# constant term holds, or below it where the data cannot see it,
# carrying next to nothing; or it settles in a cluster with the pole it
# copies, split by less than the data can resolve, the members carrying
# large residues of opposite sign. Needless poles are dropped, least
# contribution first, clusters are replaced by their mean and relocated
# again, and each simplification is kept only if the refitted fit stays
# within `pruneslack` of the fit before any pole was dropped, in its
# largest deviation over the samples and in its rms over every entry and
# sample, each against its own value -- a budget for the whole pruning,
# not for each deletion. The rms is held as well because the largest
# deviation can be pinned by one feature no pole set follows, a seam or
# a spike in measured data, and would let the rest of the fit worsen
# unseen.
function prunepoles(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64},
        fitting::VectorFitting; dc::Matrix{Float64} = nodc, lastpole::Symbol = :refuse,
        weights::Array{Float64,3} = noweights, components::Matrix{ComplexF64} = relocationcolumns(S, weights),
        scale::Real = largestopnorm(weightedsamples(S, weights)), proper::Bool = false)
    # the budget is measured from before any pole was dropped, with the
    # roundoff of the data, which bounds the rms per entry as well as the
    # largest deviation
    worst, rms = fiterrors(S, ws, poles; dc = dc, weights = weights, proper = proper)
    worstlimit = (1 + fitting.pruneslack)*worst + fitroundoff*scale
    rmslimit = (1 + fitting.pruneslack)*rms + fitroundoff*scale
    acceptable(err, _) = err[1] <= worstlimit && err[2] <= rmslimit
    before = (worst, rms)
    # dropneedless returns where no drop is acceptable, so a pass whose
    # merge leaves the poles as the drops left them has no successor to
    # change them
    while true
        poles, before = dropneedless(S, ws, poles, before, acceptable; dc = dc, lastpole = lastpole,
            scale = scale, weights = weights, proper = proper)
        kept = length(poles)
        poles, before = mergeclusters(S, ws, poles, fitting, before, acceptable; dc = dc,
            components = components, scale = scale, weights = weights, proper = proper)
        length(poles) == kept && break
    end
    return poles
end
# a pole and its conjugate are dropped when the fit without them is as
# close, the pole contributing least to the fit tried first
function dropneedless(S, ws, poles, before, acceptable; dc::Matrix{Float64} = nodc, lastpole::Symbol = :refuse,
        scale::Real, weights::Array{Float64,3} = noweights, proper::Bool = false)
    n = size(S, 1)
    while true
        residues, _ = fitresidues(S, ws, poles; dc = dc, weights = weights, constant = proper ? zeros(n, n) : nothing)
        # the norm of the residue does not vary over the samples, so it
        # is taken once per pole rather than once per pole and sample;
        # weighted, it does, and the weighed residue's Frobenius norm at
        # each sample stands in, the order being a guess the errors of
        # the trials settle
        contribution(p) = isempty(weights) ? opnorm(view(residues, :, :, p))*
            maximum(1/abs(im*w - poles[p]) for w in ws) :
            maximum(norm(view(weights, :, :, k) .* view(residues, :, :, p))/abs(im*ws[k] - poles[p])
                for k in eachindex(ws))
        upper = filter(p -> imag(poles[p]) >= 0, eachindex(poles))
        candidates = upper[sortperm([contribution(p) for p in upper])]
        dropped = false
        for p in candidates
            keep = [q for q in eachindex(poles) if q != p && !(imag(poles[p]) != 0 && poles[q] == conj(poles[p]))]
            trial = poles[keep]
            err = fiterrors(S, ws, trial; dc = dc, weights = weights, proper = proper)
            if acceptable(err, before)
                # the last pole is needless where a constant fits the data
                # as closely as the poles do: refused, returned as the
                # constant, or kept as the caller asked
                if isempty(trial)
                    lastpole == :refuse && throw(ArgumentError(nopolereason(first(err), scale)))
                    lastpole == :keep && continue
                end
                poles, before, dropped = trial, err, true
                break
            end
        end
        dropped || return poles, before
    end
end
# Why a fit which keeps no pole is refused: a constant which reproduces
# the data is no rational block, and poles which fit the data no closer
# than a constant does are too few for how fast it turns over the band,
# as a delay of many turns does; `err` is the constant's error, `scale`
# the samples' largest response.
function nopolereason(err, scale)
    err <= fitroundoff*scale && return "no pole is needed: a constant reproduces the data, so give it as a ScatteringParameters matrix."
    return "the poles fit the data no closer than a constant, which misses it by $(err/scale) of its largest response: the data turns over the band faster than the poles can follow, as a delay of many turns does. Fit with more poles, over a narrower band, or take a delay out of the data with delays."
end

function mergeclusters(S, ws, poles, fitting::VectorFitting, before, acceptable;
        dc::Matrix{Float64} = nodc, components::Matrix{ComplexF64}, scale::Real,
        weights::Array{Float64,3} = noweights, proper::Bool = false)
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
        push!(merged, a)
        imag(a) == 0 || push!(merged, conj(a))
    end
    length(merged) == length(poles) && return poles, before
    merged = converge(S, ws, merged, fitting; dc = dc, weights = weights, components = components, scale = scale,
        proper = proper)
    err = fiterrors(S, ws, merged; dc = dc, weights = weights, proper = proper)
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
        if isrealpole(a)
            Phi[:, p] .= 1 ./ (im .* ws .- real(a))
            p += 1
        else
            # `u + v` and `i (u - v)`, `u = 1/(s - a)` and `v` its conjugate
            # pole's, formed without their difference (see pairbasis)
            for k in 1:K
                Phi[k, p], Phi[k, p + 1] = pairbasis(a, im*ws[k])
            end
            p += 2
        end
    end
    return Phi
end

# A pole pair's two functions of the real basis at the point `s`, `u + v`
# and `i (u - v)` with `u = 1/(s - a)` and `v = 1/(s - conj(a))`, formed as
# `2 (s - Re a)/d` and `-2 Im a/d` with `d = (s - a)(s - conj(a))`: the
# difference of `u` and `v` cancels where the pair lies near the real axis,
# and the quotients carry a few units of roundoff of their own values.
function pairbasis(a::ComplexF64, s::Complex)
    d = (s - a)*(s - conj(a))
    return 2*(s - real(a))/d, -2*imag(a)/d
end
function complexresidues(c::AbstractVector, poles::Vector{ComplexF64})
    N = length(poles)
    r = zeros(ComplexF64, N)
    p = 1
    while p <= N
        a = poles[p]
        if isrealpole(a)
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

# The data the relocation runs on. The relaxed weight depends on the
# entries only through quadratic forms in them, each entry taken as the
# real `2K`-vector `[real(S_ij); imag(S_ij)]` of its samples: the
# entry's reduced block, the column scaling and the relaxation row. Any
# set of vectors with the entries' Gram matrix `sum_ij x_ij x_ij'`
# relocates the poles the same, and the singular components of the
# `2K x n^2` matrix of stacked entries are one, at most `2K` of them
# where there are `n^2` entries: component `j` is `U[:, j] sigma_j`
# read back as a complex `K`-vector. Every component whose singular
# value stands above the roundoff of the decomposition, `eps` times its
# larger dimension times the largest singular value, is kept, so the
# relocation is the same to roundoff, at a cost in proportion to the
# components rather than to the square of the port count. Returns the
# components as the columns of a `K x r` matrix.
function relocationcomponents(S::AbstractArray{<:Complex,3})
    n, K = size(S, 1), size(S, 3)
    X = Matrix{Float64}(undef, 2K, n*n)
    @inbounds for j in 1:n, i in 1:n
        e = (j - 1)*n + i
        for k in 1:K
            X[k, e] = real(S[i, j, k])
            X[K + k, e] = imag(S[i, j, k])
        end
    end
    # the left singular vectors alone: the right ones are as many as the
    # entries and are not needed
    U, s, _ = LAPACK.gesvd!('S', 'N', X)
    r = isempty(s) ? 0 : count(>(eps(Float64)*max(2K, n*n)*first(s)), s)
    C = Matrix{ComplexF64}(undef, K, r)
    @inbounds for j in 1:r, k in 1:K
        C[k, j] = complex(U[k, j], U[K + k, j])*s[j]
    end
    return C
end

# `f` applied to `1:m` in chunks shared among the threads Julia runs
# with, each running BLAS on one thread, whose own threads would otherwise
# multiply with Julia's, where each of the `m` pieces takes `work` floating
# point operations that pay for it (see sharedwork); a single thread, or
# pieces smaller than that, take the whole of it and leave BLAS the
# threads it has
function sharedchunks(f, m::Int, work::Real)
    if Threads.nthreads() == 1 || work < sharedwork
        f(1:m)
    else
        withoneblasthread() do
            @sync for chunk in Iterators.partition(1:m, max(cld(m, Threads.nthreads()), 1))
                Threads.@spawn f(chunk)
            end
        end
    end
    return nothing
end

# one relocation of the poles: the least squares of every column of `E`,
# the entries or their components (see relocationcomponents), with the
# shared weight, the weight's residues, and the zeros of the weight.
# `weights`, where the samples carry them (see RationalScattering), weigh
# each entry's rows at each sample, and `E` is then the entries, the
# `(j - 1) n + i`th column the samples of entry `(i, j)`. A `proper` fit's
# columns have no constant of their own.
function relocate(E::AbstractMatrix{<:Complex}, ws::AbstractVector, poles::Vector{ComplexF64};
        weights::Array{Float64,3} = noweights, proper::Bool = false, weightfloor::Real = 1e-8)
    K, ne, N = size(E, 1), size(E, 2), length(poles)
    weighted = !isempty(weights)
    weighted && size(weights, 1)^2 != ne && throw(ArgumentError("weighted samples are relocated on their entries."))
    Wf = weighted ? reshape(weights, ne, K) : zeros(0, 0)
    Phi = realbasis(ws, poles)
    # The unknowns are every column's residues and constant plus the
    # shared weight's residues and constant, all real; the equations,
    # every column's samples with a zero right hand side plus the
    # relaxation. A column's own unknowns enter only its own rows, so
    # they are eliminated column by column: the QR of the column's block
    # `[A_e B_e]` leaves the trailing block of its triangle as the
    # column's contribution to the weight's system, and the weight is
    # solved from all of them and the relaxation at once, so nothing the
    # size of every column's every sample is ever formed (Deschrijver's
    # fast fit).
    K >= N + 1 || throw(ArgumentError(lazy"fitting $(N) poles needs at least $(N + 1) sample frequencies."))
    Ae = proper ? [real.(Phi); imag.(Phi)] : [real.(Phi) ones(K); imag.(Phi) zeros(K)]
    m = size(Ae, 2)
    # `A_e` is the same matrix for every column, so its reflectors are
    # computed once and applied to each column's `B_e`: with
    # `A_e = Q [R1; 0]`, the triangle of `[A_e B_e]` has as its trailing
    # block the triangle of the rows of `Q' B_e` below `R1`, which is a
    # factorization of `N + 1` columns per column of `E` instead of
    # `N + 1 + m`. Weighted, each entry's rows are its own weights times
    # those of `A_e`, and each entry takes reflectors of its own.
    F = weighted ? nothing : qr(Ae)
    # the weight's columns scaled to their size over every column, from
    # each sample's energy over the columns, `sum_e |E[k, e]|^2`, each
    # entry weighted where the samples are
    energy = zeros(K)
    @inbounds for e in 1:ne, k in 1:K
        energy[k] += abs2(weighted ? Wf[e, k]*E[k, e] : E[k, e])
    end
    squares = [[sum(k -> energy[k]*abs2(Phi[k, p]), 1:K) for p in 1:N]; sum(energy)]
    scale = max.(sqrt.(squares), floatmin(Float64))
    reduced = zeros(ne*(N + 1) + 1, N + 1)
    rhs = zeros(ne*(N + 1) + 1)
    # The columns are independent, each writing its own rows, so they are
    # shared among the threads, each with a block of its own. Every
    # variable the closure reads is assigned once, so that none of them is
    # boxed.
    function rows!(cols)
        Be = zeros(2K, N + 1)
        Aw = zeros(weighted ? 2K : 0, m)
        for e in cols
            @inbounds for p in 1:N + 1, k in 1:K
                v = p <= N ? E[k, e]*Phi[k, p] : E[k, e]
                weighted && (v *= Wf[e, k])
                Be[k, p] = -real(v)/scale[p]
                Be[K + k, p] = -imag(v)/scale[p]
            end
            if weighted
                @inbounds for p in 1:m, k in 1:K
                    Aw[k, p] = Wf[e, k]*Ae[k, p]
                    Aw[K + k, p] = Wf[e, k]*Ae[K + k, p]
                end
                lmul!(qr!(Aw).Q', Be)
            else
                lmul!(F.Q', Be)
            end
            T = qr!(view(Be, m + 1:2K, :)).factors
            reduced[(e - 1)*(N + 1) + 1:e*(N + 1), :] .= UpperTriangular(view(T, 1:N + 1, :))
        end
        return nothing
    end
    sharedchunks(rows!, ne, 4K*(N + 1)^2)
    # the relaxation: the real part of the weight averages to one; scaled
    # to the size of the data so that it weighs as one sample of it
    weight = sqrt(squares[N + 1]/K)
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
    zeros_ = eigvals(Ar .- br*transpose(csigma) ./ dsigma)
    newpoles = ComplexF64[]
    for z in zeros_
        real(z) > 0 && (z = complex(-real(z), imag(z)))
        push!(newpoles, z)
    end
    # conjugate pairs paired, real poles real: the eigenvalues of a real
    # matrix come in pairs of exact conjugates, which the reflection into
    # the left half plane keeps
    sort!(newpoles; by = x -> (real(x), imag(x)))
    out = ComplexF64[]
    used = falses(length(newpoles))
    for (k, z) in enumerate(newpoles)
        used[k] && continue
        used[k] = true
        partner = isrealpole(z) ? nothing : findfirst(l -> !used[l] && newpoles[l] == conj(z), eachindex(newpoles))
        if isnothing(partner)
            push!(out, complex(real(z), 0.0))
        else
            push!(out, complex(real(z), abs(imag(z))), complex(real(z), -abs(imag(z))))
            used[partner] = true
        end
    end
    return out
end

# the real block diagonal matrix of a pole set and its input column of
# ones in the real form, so that `A - b c'/d` has the weight's zeros as
# eigenvalues, with `c` the weight's coefficients in the real basis,
# which are its residues in that form
function realpolematrix(poles::Vector{ComplexF64})
    N = length(poles)
    A = zeros(N, N)
    b = zeros(N)
    p = 1
    while p <= N
        a = poles[p]
        if isrealpole(a)
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

# The least squares of the residue solve, factored once and applied to
# every entry's right hand side at once (see fitcoefficients). Directions
# of the basis whose singular value is below `fitranktolerance` of the
# largest are left out of the solution: poles which have coalesced give
# columns the samples cannot tell apart, and a plain least squares
# answers them with large cancelling residues which corrupt the fit
# error and the value at zero frequency. The minimum norm solution leaves
# those directions at zero.
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

# The solution FitLeastSquares gives for one right hand side `b`, from the
# QR of `A`, whose triangle has the singular values of `A`, and the singular
# value decomposition of the triangle, `N + 1` square, rather than of `A`;
# `A` and `b` are overwritten.
function leastsquaressolution!(A::Matrix{Float64}, b::Vector{Float64})
    scale = [norm(view(A, :, j)) for j in axes(A, 2)]
    for j in eachindex(scale)
        scale[j] > 0 || (scale[j] = 1.0)
    end
    A ./= transpose(scale)
    F = qr!(A)
    lmul!(F.Q', b)
    k = min(size(A)...)
    T = svd(F.R)
    kept = isempty(T.S) ? 0 : count(>(fitranktolerance*first(T.S)), T.S)
    return (view(T.V, :, 1:kept)*((transpose(view(T.U, :, 1:kept))*view(b, 1:k)) ./ view(T.S, 1:kept))) ./ scale
end

# the residues of every entry and the constant at fixed poles, by least
# squares in the real basis (see fitcoefficients). With `constant` given,
# that matrix is held and only the strictly proper part is fitted, to
# `S - constant`.
function fitresidues(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64};
        constant::Union{Nothing,AbstractMatrix} = nothing,
        dc::Matrix{Float64} = nodc, weights::Array{Float64,3} = noweights)
    n, N = size(S, 1), length(poles)
    fixed = !isnothing(constant)
    residues = zeros(ComplexF64, n, n, N)
    D = fixed ? Matrix{Float64}(constant) : zeros(n, n)
    X, _ = fitcoefficients(S, ws, poles; constant = constant, dc = dc, weights = weights)
    for j in 1:n, i in 1:n
        e = (j - 1)*n + i
        residues[i, j, :] .= complexresidues(view(X, e, 1:N), poles)
        fixed || (D[i, j] = X[e, N + 1])
    end
    return residues, D
end

# Every entry's coefficients in the real basis at fixed poles, as the rows
# of `X`, and the basis `M`, its functions' real parts at the samples over
# their imaginary parts, whose columns the coefficients weigh: the
# residues' real and imaginary parts, and the constant where it is not
# held. Every entry shares the least squares, so its coefficients on the
# kept left singular vectors, `U' b` for the entry's real samples `b`, are
# one product with the samples read as the complex `n^2 x K` matrix they
# are stored as, and the solutions `V (U' b) ./ s ./ scale` one more (see
# FitLeastSquares); a held constant and a stated value at zero frequency
# enter both as terms of rank one.
function fitcoefficients(S::AbstractArray{<:Complex,3}, ws::AbstractVector, poles::Vector{ComplexF64};
        constant::Union{Nothing,AbstractMatrix} = nothing,
        dc::Matrix{Float64} = nodc, weights::Array{Float64,3} = noweights)
    n, K, N = size(S, 1), length(ws), length(poles)
    fixed = !isnothing(constant)
    Phi = realbasis(ws, poles)
    cols = fixed ? N : N + 1
    M = zeros(2K, cols)
    M[1:K, 1:N] .= real.(Phi)
    M[K + 1:2K, 1:N] .= imag.(Phi)
    fixed || (M[1:K, N + 1] .= 1.0)
    held = fixed ? vec(Matrix{Float64}(constant)) : zeros(0)
    Sf = reshape(S, n*n, K)
    # The value at zero frequency is the same combination of the same
    # basis, read at zero instead of on the axis, and linear in the
    # unknowns once the poles are settled, so stating it is one equality
    # per entry, met exactly: a particular solution, then a least squares
    # over the null space of the equality. A sample near zero would
    # instead ask the fit to take that value at a frequency where it does
    # not, at an excursion far worse than the extrapolation it corrects.
    stated = !isempty(dc)
    if stated
        row = vec(real.(realbasis([0.0], poles)))
        c = fixed ? row : vcat(row, 1.0)
        cc = dot(c, c)
        cc > 0 || throw(ArgumentError("the zero frequency row of the basis vanishes, so the value there cannot be stated."))
        Z = nullspace(reshape(c, 1, :))
        # the particular solution `c t_e/cc` of each entry's target `t_e`,
        # and what its samples `g t_e` take from the coefficients
        g = M*c ./ cc
        t = vec(real.(dc)) .- (fixed ? held : 0.0)
    end
    if !isempty(weights)
        X = stated ? weightedcoefficients(M, Sf, reshape(weights, n*n, K), held, Z, t .* transpose(c ./ cc), t, g) :
            weightedcoefficients(M, Sf, reshape(weights, n*n, K), held, nothing, zeros(n*n, 0), zeros(0), zeros(0))
        return X, M
    end
    if !stated
        F = FitLeastSquares(M)
        return fitsolutions(F, first(leftcoefficients(F, Sf, held))), M
    end
    F = FitLeastSquares(M*Z)
    C, Ut, Ub = leftcoefficients(F, Sf, held)
    C .-= t .* transpose(transpose(Ut)*view(g, 1:K) .+ transpose(Ub)*view(g, K + 1:2K))
    return t .* transpose(c ./ cc) .+ fitsolutions(F, C)*transpose(Z), M
end

# Weighted, each entry has a least squares of its own, the rows of the
# basis `M` and of the entry's samples weighed by its weights, `[w_e; w_e]`,
# the entries shared among the threads (see sharedchunks). A held constant
# `held` comes off the real parts of the samples; a stated value, where
# `Z` is the null space of its equality, leaves each entry's particular
# solution, `particular[e, :]`, whose samples are `t[e] g`, and a least
# squares over `Z`.
function weightedcoefficients(M::Matrix{Float64}, Sf::AbstractMatrix{<:Complex}, Wf::AbstractMatrix{Float64},
        held::Vector{Float64}, Z::Union{Nothing,Matrix{Float64}}, particular::Matrix{Float64},
        t::Vector{Float64}, g::Vector{Float64})
    ne, K = size(Sf)
    B = isnothing(Z) ? M : M*Z
    X = zeros(ne, size(M, 2))
    function solve!(entries)
        b, w, A = zeros(2K), zeros(2K), similar(B)
        for e in entries
            @inbounds for k in 1:K
                w[k] = w[K + k] = Wf[e, k]
                b[k] = real(Sf[e, k]) - (isempty(held) ? 0.0 : held[e])
                b[K + k] = imag(Sf[e, k])
            end
            isnothing(Z) || (b .-= t[e] .* g)
            A .= w .* B
            b .*= w
            x = leastsquaressolution!(A, b)
            isnothing(Z) ? (X[e, :] .= x) : (X[e, :] .= view(particular, e, :) .+ Z*x)
        end
        return nothing
    end
    sharedchunks(solve!, ne, 4K*size(B, 2)^2)
    return X
end

# every entry's coefficients on the kept left singular vectors of `F`,
# `Re(Sf (Ut - i Ub))`, less those of a held constant `held`, empty for
# none, whose samples are `[D_e; 0]`; and `Ut`, `Ub`
function leftcoefficients(F::FitLeastSquares, Sf::AbstractMatrix{<:Complex}, held::Vector{Float64})
    K = size(Sf, 2)
    Ut, Ub = F.U[1:K, :], F.U[K + 1:2K, :]
    C = real.(Sf*complex.(Ut, .-Ub))
    isempty(held) || (C .-= held .* transpose(vec(sum(Ut; dims = 1))))
    return C, Ut, Ub
end
# the solutions of the coefficients `C`, `(C ./ s') V' ./ scale'`
fitsolutions(F::FitLeastSquares, C::Matrix{Float64}) = ((C ./ transpose(F.s))*transpose(F.V)) ./ transpose(F.scale)

"""
    passiveconstant(D; tol = 1e-6, margin = tol)

The constant term `D` with every singular value above `sqrt(1 + tol)`,
where its dissipation `I - D'D` falls below `-tol`, brought to
`1 - margin`, leaving the singular vectors alone. `D` is the model at
infinite frequency. The passivity enforcement changes it with the
residues, by the least change over the samples; bringing it under one
before the residues are fitted instead has the residues fitted to what
is left over the samples, the least change of the fit with its poles.
`tol` is how far below zero the dissipation must fall for a singular
value to be treated as active, and belongs with the tolerance the block
is validated against;
`margin` is how far under one it is put, and belongs with the margin the
enforcement aims for. Follows Gustavsen, IEEE Transactions on
Electromagnetic Compatibility 67(3), 2025, section X-A.
"""
function passiveconstant(D::AbstractMatrix; tol::Real = 1e-6, margin::Real = tol)
    F = svd(D)
    # only a constant term genuinely above one is brought down. A block
    # which is lossless at infinite frequency has a unitary `D` by right,
    # a series inductor between ports for one, and moving its singular
    # values would put an error into a fit which was exact.
    level = passivelevel(tol)
    any(>(level), F.S) || return Matrix{Float64}(D), false
    σ = [s > level ? 1 - margin : s for s in F.S]
    return F.U*Diagonal(σ)*F.Vt, true
end

# the real state space realization of poles and residue matrices, each
# residue `R` taken within the space of its rows `Q` (see residuespaces):
# `R Q = U S V'` split as `C = U sqrt(S)`, `B = sqrt(S) V' Q'`, so that
# `C B = R Q Q'` and the realization is balanced; a real pole `a` with a
# space of dimension `r` is the `r` states `A = a I`, and a conjugate
# pair `a, conj(a)` the `2r` states of the real and imaginary parts of
# `r` complex ones, each state's two adjacent,
# `A = I kron [Re a  -Im a; Im a  Re a]`, so that `A` is block diagonal
# in blocks of one and two and is its own real Schur form
function realization(poles::Vector{ComplexF64}, residues::AbstractArray{<:Complex,3}, n::Int;
        ranktol::Real = 1e-12, spaces = residuespaces(poles, residues, ranktol))
    N = length(poles)
    blocks = Tuple{Matrix{Float64},Matrix{Float64},Matrix{Float64}}[]
    p = 1
    while p <= N
        a, Q = poles[p], spaces[p]
        r = size(Q, 2)
        if isrealpole(a)
            if r > 0
                Qr = real.(Q)
                F = svd(real.(view(residues, :, :, p))*Qr)
                Bc = Diagonal(sqrt.(F.S))*F.Vt*transpose(Qr)
                Cc = F.U*Diagonal(sqrt.(F.S))
                push!(blocks, (real(a) .* Matrix(1.0I, r, r), Bc, Cc))
            end
            p += 1
        else
            if r > 0
                F = svd(view(residues, :, :, p)*Q)
                Bc = Diagonal(sqrt.(F.S))*F.Vt*Q'
                Cc = F.U*Diagonal(sqrt.(F.S))
                # the complex states `z' = a z + Bc u`, `y = 2 Re(Cc z)`, as
                # their real and imaginary parts, a state's two adjacent
                α, β = real(a), imag(a)
                Bp, Cp = zeros(2r, n), zeros(n, 2r)
                Bp[1:2:end, :] .= real.(Bc)
                Bp[2:2:end, :] .= imag.(Bc)
                Cp[:, 1:2:end] .= 2 .* real.(Cc)
                Cp[:, 2:2:end] .= -2 .* imag.(Cc)
                push!(blocks, (kron(Matrix(1.0I, r, r), [α -β; β α]), Bp, Cp))
            end
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

# the numerical rank of a residue whose singular values, largest first,
# are `s`
residuerank(s::AbstractVector{<:Real}, ranktol::Real) = count(>(ranktol*max(first(s), floatmin(Float64))), s)

# The spaces of the rows of the residues which a realization gives states
# to: for each pole, the right singular vectors of its residue's singular
# values above `ranktol` of the largest, real for a real pole, and for a
# pair the first's and their conjugates for the second; a residue of full
# rank spans every row, and its space is given by the identity, which
# depends on no choice of singular vectors. The rank of a fit is decided
# here once: the passivity enforcement changes each residue within its
# space and the realization spans the same, so that the block realized
# is the one the enforcement made passive.
function residuespaces(poles::Vector{ComplexF64}, residues::AbstractArray{<:Complex,3},
        ranktol::Real)
    n = size(residues, 1)
    spaces = Vector{Matrix{ComplexF64}}(undef, length(poles))
    space(F) = (r = residuerank(F.S, ranktol); r == n ? Matrix{ComplexF64}(I, n, n) : F.V[:, 1:r])
    p = 1
    while p <= length(poles)
        if isrealpole(poles[p])
            spaces[p] = space(svd(real.(view(residues, :, :, p))))
            p += 1
        else
            spaces[p] = space(svd(residues[:, :, p]))
            spaces[p + 1] = conj.(spaces[p])
            p += 2
        end
    end
    return spaces
end

# the states of a realization through `spaces`, a pair's two poles a
# space each
statecount(spaces::Vector{Matrix{ComplexF64}}) = sum(Q -> size(Q, 2), spaces; init = 0)

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
# meeting the constraints. The conditions of the optimum, which
# dualactiveset tests in full, hold only where the correction meets every
# constraint.
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
# Gram matrix singular; the least norm multipliers are the answer where
# those constraints agree, and dualactiveset finishes the solve where they
# do not.
function gramsolve(Gp::AbstractMatrix, cp::AbstractVector)
    S = Symmetric(Matrix(Gp))
    F = cholesky(S; check = false)
    issuccess(F) && return F \ (-collect(cp))
    return pinv(Matrix(S))*(-collect(cp))
end

# How far to contract a block toward an anchor: the largest step `t` in
# `[lo, 1]` for which `measure(t)`, the quantity to bring under `target`,
# does so, `flo` and `fhi` being the measure at `lo` and at one where the
# caller has them. Along the path `S0 + t (S - S0)` the response is affine
# in `t`, so its norm is convex and the steps meeting the target are an
# interval containing zero whose end the search finds. The step is
# found on the path itself because no formula stands in for it:
# scaling by `1/M`, with `M` a bound on the norm, answers only the
# `S0 = 0` case, and the triangle bound `t <= (1 - a)/(M - a)`, with
# `a` the anchor's own norm, is sufficient but not necessary -- at
# `a = 1` it says nothing, yet `S = -1.001 + 2.001/(s + 1)` anchored
# at `S(0) = 1` contracts at `t = 2/2.001` to the all pass
# `(1 - s)/(1 + s)`.
function contractionstep(measure, target::Real; lo::Real = 0.0, fhi::Real = measure(1.0),
        flo::Real = measure(lo), tol::Real = 1e-12, steps::Integer = 60)
    fhi <= target && return 1.0
    flo <= target || return Float64(lo)
    # The measure is convex, so the chord from a point under the target
    # to one over it crosses the target on the near side of the true
    # boundary: every chord step is a valid step, and where the measure is
    # linear, as it is toward an anchor of zero, the first chord lands on
    # the boundary. On a curved measure the chord alone creeps up on the
    # boundary while the far end of the bracket stays put, so where the
    # same end is replaced twice in a row the value kept at the other is
    # halved, which moves that end too (the Illinois step); a chord which
    # makes no progress falls back to bisecting the bracket. The search
    # stops where the bracket is within `tol`, or after `steps`
    # evaluations; the shortfall of the returned step from the boundary is
    # what the caller pays as needless contraction, so it is kept small.
    a, b, left, right = flo - target, fhi - target, Float64(lo), 1.0
    kept = 0
    for _ in 1:steps
        t = left - a*(right - left)/(b - a)
        left < t < right || (t = (left + right)/2)
        ft = measure(t) - target
        if ft <= 0
            left, a = t, ft
            kept == -1 && (b /= 2)
            kept = -1
            # a step landing exactly on the target is the boundary
            a == 0 && break
        else
            right, b = t, ft
            kept == 1 && (a /= 2)
            kept = 1
        end
        right - left <= tol && break
    end
    return left
end

# The contraction of a returned fit, which the caller is warned of: a
# fit the order search discards is not, since its contraction describes
# nothing the caller receives.
warncontraction(::Nothing) = nothing
warncontraction(c::NamedTuple) = @warn "the rational fit was contracted by $(c.factor) to make it passive, which moved its response by up to $(c.moved) over the samples; compare the fitted dissipation against the data's before trusting the noise of this block."

# what the rounds of the enforcement did, for a refusal to say
function enforcementrecord(stopped::Symbol, corrections::Int)
    rounds = corrections == 1 ? "one round of correction" : "$(corrections) rounds of correction"
    stopped === :clean && corrections == 0 && return "its search found no violation to correct"
    stopped === :clean && return "$(rounds) left no violation its search could find"
    stopped === :infeasible && return "after $(rounds) the next correction could not meet its own constraints"
    return "$(rounds) did not remove every violation"
end


# A fitted block as its poles and, over the real basis of `realbasis`,
# real coefficient matrices at the frequency scaled by `wref`:
# `S(x) = D + sum_c X[:, :, c] phi_c(x)` with `x = w/wref`. It is the form
# the fit produces, and a frequency costs `O(n^2 N)` in it where the
# Schur path of its realization costs `O(n^2 nz)`, `nz` about `n N`. Its
# residues lie within the `spaces` of their rows (see residuespaces), and
# its realization spans the same. The enforcement changes the
# coefficients within those spaces, and `D`; the poles stay.
struct ResidueForm
    poles::Vector{ComplexF64}
    X::Array{Float64,3}
    D::Matrix{Float64}
    wref::Float64
    spaces::Vector{Matrix{ComplexF64}}
end

# the form of the poles and residues of `vectorfit`, in rad/s, each
# residue taken within its space, `R Q Q'`
function ResidueForm(poles::Vector{ComplexF64}, residues::AbstractArray{<:Complex,3},
        D::AbstractMatrix, wref::Real, spaces::Vector{Matrix{ComplexF64}})
    n, N = size(D, 1), length(poles)
    X = zeros(n, n, N)
    p = 1
    while p <= N
        Q = spaces[p]
        if isrealpole(poles[p])
            R = real.(view(residues, :, :, p))
            size(Q, 2) < n && (R = (R*real.(Q))*transpose(real.(Q)))
            X[:, :, p] .= R ./ wref
            p += 1
        else
            R = residues[:, :, p]
            size(Q, 2) < n && (R = (R*Q)*Q')
            X[:, :, p] .= real.(R) ./ wref
            X[:, :, p + 1] .= imag.(R) ./ wref
            p += 2
        end
    end
    return ResidueForm(poles ./ wref, X, Matrix{Float64}(D), Float64(wref), spaces)
end

# the realization of the form at the scaled frequency
realization(m::ResidueForm) =
    realization(m.poles, scaledresidues(m), size(m.D, 1); spaces = m.spaces)

# the complex residues at the scaled frequency, conjugate pairs completed
function scaledresidues(m::ResidueForm)
    n, N = size(m.D, 1), length(m.poles)
    R = zeros(ComplexF64, n, n, N)
    p = 1
    while p <= N
        if isrealpole(m.poles[p])
            R[:, :, p] .= view(m.X, :, :, p)
            p += 1
        else
            R[:, :, p] .= complex.(view(m.X, :, :, p), view(m.X, :, :, p + 1))
            R[:, :, p + 1] .= conj.(view(R, :, :, p))
            p += 2
        end
    end
    return R
end

# `phi` the real basis of the poles at the scaled frequency `x`, the row
# of `realbasis` there with a pair's two formed without cancellation (see
# pairbasis), and zero at infinite frequency
function basisat!(phi::AbstractVector, poles::Vector{ComplexF64}, x::Real)
    fill!(phi, 0)
    isinf(x) && return phi
    s = im*x
    p = 1
    @inbounds while p <= length(poles)
        a = poles[p]
        if isrealpole(a)
            phi[p] = 1/(s - real(a))
            p += 1
        else
            phi[p], phi[p + 1] = pairbasis(a, s)
            p += 2
        end
    end
    return phi
end

# `F` the block at the scaled frequency `x`, through `phi`
function responseat!(F::AbstractMatrix, m::ResidueForm, x::Real, phi::AbstractVector)
    basisat!(phi, m.poles, x)
    n = size(m.D, 1)
    F .= m.D
    @inbounds for c in eachindex(phi)
        f = phi[c]
        iszero(f) && continue
        for j in 1:n, i in 1:n
            F[i, j] += m.X[i, j, c]*f
        end
    end
    return F
end

# The unknowns of one output port's row of the correction, the same for
# every port: its row of the change of `D`, then pole by pole its row of
# the change of the residue in the orthonormal basis `Q_p` of the
# residue's space (see residuespaces), `dR_p = Z_p Q_p'`, so that a
# residue stays within its space and the realization keeps its states,
# and the coordinates carry no scale of their own. Each block
# `(c, cols, M)` says that the port's unknowns `x[cols]` enter its row as
# `phi_c(x) (M x[cols])'`, `c = 0` being the constant.
function correctioncoordinates(m::ResidueForm)
    n, N = size(m.D, 1), length(m.poles)
    blocks = Tuple{Int,UnitRange{Int},Matrix{Float64}}[(0, 1:n, Matrix{Float64}(I, n, n))]
    k = n
    p = 1
    while p <= N
        Q = m.spaces[p]
        r = size(Q, 2)
        if isrealpole(m.poles[p])
            push!(blocks, (p, k + 1:k + r, real.(Q)))
            k += r
            p += 1
        else
            # a row `z Q'` of the pair's residue, `z = a + i b`, has real
            # part `a QR' + b QI'` and imaginary part `b QR' - a QI'`
            QR, QI = real.(Q), imag.(Q)
            push!(blocks, (p, k + 1:k + 2r, [QR QI]))
            push!(blocks, (p + 1, k + 1:k + 2r, [-QI QR]))
            k += 2r
            p += 2
        end
    end
    return blocks, k
end

# `g` the unknowns' response to the direction `v` at a frequency whose
# basis is `phi`: `g[u] = f_u(x) v` for the row function `f_u` of each
# unknown
function coordinateresponse!(g::AbstractVector, blocks, phi::AbstractVector, v::AbstractVector)
    fill!(g, 0)
    for (c, cols, M) in blocks
        f = c == 0 ? one(eltype(phi)) : phi[c]
        iszero(f) && continue
        mul!(view(g, cols), transpose(M), v, f, true)
    end
    return g
end

# The norm of the correction over the samples in these coordinates, for
# one port's unknowns: `sum_k ||dS_i(x_k)||^2 = x_i' H x_i`, with
# `H = sum_{c, c'} Re(P_cc') M_c' M_c'` and `P` the Gram matrix of the
# basis over the samples, so the samples enter once, through the
# `(N + 1)`-square `P` (with the residues' own coordinates it is the same
# block for every entry: Gustavsen, IEEE Transactions on Electromagnetic
# Compatibility 67(3), 2025, section III-C).
function correctiongram(blocks, nu::Int, poles::Vector{ComplexF64}, xs::AbstractVector)
    Phi = hcat(ones(length(xs)), realbasis(xs, poles))
    P = real.(transpose(Phi)*conj.(Phi))
    H = zeros(nu, nu)
    for (c1, cols1, M1) in blocks, (c2, cols2, M2) in blocks
        w = P[c1 + 1, c2 + 1]
        iszero(w) && continue
        mul!(view(H, cols1, cols2), transpose(M1), M2, w, true)
    end
    return (H .+ transpose(H))./2
end

# The metric of the least change, factored: `L L'` the normal matrix of one
# port's unknowns (see correctiongram), with `regularization` times its mean
# diagonal entry on its diagonal, in the coordinates `Z` leaves free where a
# value at zero frequency is held, `nothing` where none is. That value is
# held exactly: the change there is `sum_c phi_c(0) dX_c + dD`, which for
# one port's row reads `E x_i` with `E` the same for every port, so the
# correction is confined to the nullspace of `E`. Where every residue has
# full rank its space is the identity (see residuespaces) and a port's
# unknowns are the changes of its row's coefficients, `n` to a basis
# function, so the normal matrix is `P kron I` with `P` the Gram matrix of
# the basis over the samples, and `E` is `e' kron I` with `e` the basis at
# zero: the regularization, the nullspace and the factor are `P`'s and
# `e`'s, `N + 1` square, taken `stride = n` unknowns at a time along the
# basis (see alongbasis). Otherwise they are the whole normal matrix's, and
# `stride` is one. `L` holds that one factor, the same for every port,
# where the samples are not weighted. Weighted, each entry `(i, j)` has its
# own Gram matrix of the basis, its samples' squared weights weighing them,
# and each port its own normal matrix: `L` holds a factor for each entry,
# the `(j - 1) n + i`th of the port `i`'s unknowns `j`, `j + n`, ... along
# the basis where every residue has full rank, and otherwise a factor for
# each port (see portsolve).
struct CorrectionMetric
    L::Vector{LowerTriangular{Float64,Matrix{Float64}}}
    Z::Union{Nothing,Matrix{Float64}}
    stride::Int
end

# The bytes the correction metric of a weighted fit takes at its peak, in
# the representation correctionmetric chooses: each entry's Gram matrix of
# the basis, `N + 1` square, and a factor of it for each entry where every
# residue has full rank, and otherwise a factor of the `n` plus the states'
# unknowns of a port's row for each port; and the few square matrices of
# its construction. The full rank metric drops each Gram matrix once it is
# factored, and they are counted all the same, the collector being free to
# hold them until the metric is done. Without weights the metric is one
# factor, within what `statebytes` allows a state (see defaultmaxstates).
function weightedmetricbytes(spaces::Vector{Matrix{ComplexF64}}, n::Int)
    q = length(spaces) + 1.0
    fullrank = all(Q -> size(Q, 2) == n, spaces)
    k, factors = fullrank ? (q, n^2) : (n + statecount(spaces), n)
    return 8*((factors + 3)*k^2 + n^2*q^2)
end

function correctionmetric(m::ResidueForm, blocks, nu::Int, xs::AbstractVector,
        regularization::Real, dc::Matrix{Float64}, weights::Array{Float64,3} = noweights)
    n, N = size(m.D, 1), length(m.poles)
    phi = basisat!(zeros(ComplexF64, N), m.poles, 0.0)
    # a normal matrix regularized and taken into the free coordinates,
    # factored; one which is not positive definite in Float64, where the
    # squares of an entry's weights underflow, or where the regularization
    # lifts too little of what the samples leave undetermined, is refused
    # as a fit which cannot be made passive
    function factored(H, Z)
        H = (H .+ transpose(H))./2
        view(H, diagind(H)) .+= regularization*tr(H)/size(H, 1)
        isnothing(Z) || (H = transpose(Z)*H*Z)
        F = cholesky(Symmetric(H); check = false)
        issuccess(F) || throw(ArgumentError("the fit could not be made passive: the metric of its correction, its change over the samples weighed by the squares of any weights, is not positive definite in floating point. Narrow the range of the weights, raise regularization, or fit with fewer poles."))
        return F.L
    end
    fullrank = all(Q -> Q == I, m.spaces)
    Z = nothing
    if fullrank
        isempty(dc) || (Z = nullspace(permutedims(vcat(1.0, real.(phi)))))
    elseif !isempty(dc)
        E = zeros(n, nu)
        for (c, cols, M) in blocks
            E[:, cols] .+= real(c == 0 ? 1.0 : phi[c]) .* M
        end
        Z = nullspace(E)
        size(Z, 2) == nu - n || throw(ArgumentError(
            lazy"the zero frequency condition leaves $(size(Z, 2)) free directions per port where $(nu - n) were expected; report this."))
    end
    isempty(weights) && !fullrank &&
        return CorrectionMetric([factored(correctiongram(blocks, nu, m.poles, xs), Z)], Z, 1)
    Phi = hcat(ones(length(xs)), realbasis(xs, m.poles))
    isempty(weights) && return CorrectionMetric([factored(real.(transpose(Phi)*conj.(Phi)), Z)], Z, n)
    # each entry's Gram matrix of the basis, its samples weighed by their
    # squared weights, which where every residue has full rank is needed
    # only for that entry's factor and is dropped once factored
    gram(i, j) = real.(transpose(Phi)*(view(weights, i, j, :).^2 .* conj.(Phi)))
    fullrank && return CorrectionMetric([factored(gram(i, j), Z) for j in 1:n for i in 1:n], Z, n)
    grams = [gram(i, j) for i in 1:n, j in 1:n]
    # each port's normal matrix: its row's change over the samples,
    # `sum_j sum_k w_ij(k)^2 |sum_c phi_c(x_k) (M_c x[cols_c])_j|^2`
    factors = map(1:n) do i
        H = zeros(nu, nu)
        for j in 1:n
            P = grams[i, j]
            for (c1, cols1, M1) in blocks, (c2, cols2, M2) in blocks
                w = P[c1 + 1, c2 + 1]
                iszero(w) && continue
                mul!(view(H, cols1, cols2), view(M1, j:j, :)', view(M2, j:j, :), w, true)
            end
        end
        factored(H, Z)
    end
    return CorrectionMetric(factors, Z, 1)
end

# `X`, a port's unknowns down each column, through port `i`'s factor of a
# weighted `metric`, `L_i^-1 X`, or with `adjoint` through its
# transpose's: one factor for the port, or where the factors are the
# entries' one for each of the port's unknowns `j`, `j + n`, ... along
# the basis (see CorrectionMetric)
function portsolve(metric::CorrectionMetric, i::Int, X::AbstractMatrix{Float64}; adjoint::Bool = false)
    (; L, stride) = metric
    stride == 1 && return adjoint ? transpose(L[i]) \ X : L[i] \ X
    Y = similar(X)
    for j in 1:stride
        rows = j:stride:size(X, 1)
        F = L[(j - 1)*stride + i]
        Y[rows, :] .= adjoint ? transpose(F) \ X[rows, :] : F \ X[rows, :]
    end
    return Y
end

# `f` applied to the columns of `G` along the basis: `G`'s rows are a
# port's unknowns, `stride` to a basis function, and `f` acts on the basis
# index of each of them and of each column alike, as `f kron I` would on
# `G` itself
function alongbasis(f, G::AbstractMatrix, stride::Int)
    stride == 1 && return f(G)
    q, k = size(G, 1) ÷ stride, size(G, 2)
    X = reshape(permutedims(reshape(G, stride, q, k), (2, 1, 3)), q, stride*k)
    Y = f(X)
    return reshape(permutedims(reshape(Y, size(Y, 1), stride, k), (2, 1, 3)), stride*size(Y, 1), k)
end

# The bands of scaled frequency where the largest singular value stands
# above `level`, found by sweeping a logarithmic grid of `density` points
# per pole reaching `gridpad` times past the samples, each decided by
# `exceeds!`, and five points at each distinct complex pole, where a
# resonance narrower than the grid's spacing peaks, a pole found above the
# level adding a band of `polespan` half widths either side.
# The grid follows the poles and not the states: a multiport's residues
# give a pole as many states as their rank, but the largest singular
# value varies only as fast as the poles let it. The certified sweep the
# enforcement falls back on settles every frequency (see passivitysweep),
# but where the fit stands near the level its intervals shrink as the
# square root of the margin, and the grid's points are few; Gustavsen
# (IEEE Transactions on Electromagnetic Compatibility 67(3), 2025,
# section X-B) orders a grid and an exact test the same way.
function residueviolations(m::ResidueForm, xs::AbstractVector, level::Real; density::Integer = 12)
    n = size(m.D, 1)
    F, W, phi = zeros(ComplexF64, n, n), zeros(ComplexF64, n, n), zeros(ComplexF64, length(m.poles))
    above(x) = exceeds!(W, responseat!(F, m, x, phi), level)
    lo = xs[something(findfirst(>(0), xs), lastindex(xs))]
    hi = last(xs)
    grid = vcat(0.0, exp.(range(log(lo/gridpad), log(hi*gridpad); length = density*max(length(m.poles), 1))))
    bands = hitbands(grid, [above(x) for x in grid])
    for l in distinctpoles(m.poles)
        imag(l) > 0 || continue
        w0, half = imag(l), max(abs(real(l)), eps())
        any(t -> above(w0 + t*half), (-1.0, -0.5, 0.0, 0.5, 1.0)) || continue
        push!(bands, polebracket(l, polespan))
    end
    return bands
end

# The sweep of a form at a level over every frequency, zero and infinity
# included, which settles each interval of the axis by a bound on the
# form's dissipation across it. On the axis the dissipation
# `Phi = I - S' S` is Hermitian, and the form stands under the level where
# its least eigenvalue is at least `1 - level^2`. It is the rational
# function `I - S(-s)^T S(s)`, whose poles are the form's and their
# negatives (see SweepWork). The dissipation of a lossless form vanishes,
# residues and all, and so does its bound, where a bound on the largest
# singular value, which stands at one, settles an interval only once it is
# narrow. On an interval of centre `c` and half width `h` of the scaled
# frequency `x`, the dissipation at `c + t` is `P0 + t P1` and a remainder:
# `P0 = I - S0' S0` and `P1 = -(S1' S0 + S0' S1)`, `S0` and `S1` the response
# and its derivative at the centre, summed in the form's real basis. With
# `a_j = i c - q_j` for a pole `q_j` of the dissipation and `A_j` its
# residue, its term is exactly
# `A_j/a_j - i t A_j/a_j^2 - t^2 A_j/(a_j^2 (a_j + i t))`, so the remainder is
# at most `h^2 sum_j ||A_j||/(|a_j|^2 d_j)`, `d_j` the distance from `q_j` to
# the interval. The least eigenvalue of the Hermitian matrix `P0 + t P1` is
# concave in `t`, so it is least at an end. Past four times the largest
# pole magnitude the variable is `y = 1/x`, a term being `A_j y/(i - q_j y)`,
# whose remainder's coefficient is `||A_j||/(|b_j|^2 e_j)` with
# `b_j = i - q_j y0` and `e_j` the distance from `i/q_j` to the interval, so
# that the last interval ends at `y = 0`, where the dissipation is
# `I - D' D`. Before any interval the dissipation is bounded over the
# whole axis at once (see everywhereunder), which settles a lossless form
# without one. An interval whose bound, less the roundoff of the response,
# of the products and of the factorization, is at least `1 - level^2` is
# settled by two Cholesky tests (see settled). One whose centre stands
# above the level belongs to a band. Where the centre stands under the
# level by less than the remainder over the interval narrowed `halvings`
# times, so that halving would take as many generations at least, and the
# roundoff leaves the margin for it (see secondguard), the expansion to
# second order is tried (see secondorder!), which settles the intervals
# where the form is lossless in some directions and lossy in others with
# far fewer halvings. One whose centre stands within the roundoff of the
# level, where the remainder over it is within the roundoff too or the
# roundoff leaves the second order no room, or whose bound's own width is
# within the roundoff, halving cannot settle, and it belongs to a band;
# any other is halved.
# The intervals go a generation at a time, the halves of one making the
# next, and a generation's are shared among the threads (see
# eachinterval!). Returns the verdict, `:passive` where every interval was
# settled, `:active` where a centre stood above the level,
# `:indeterminate` where intervals were left unresolved and none above it,
# and `:budget` where the clock, `time_ns()`, reached `deadline` before
# every interval was settled; the bands, merged, in the scaled frequency,
# the interval at which the time ran out for `:budget`; the centre of each
# band where the response stood highest; and the intervals evaluated, none
# where the bound over the whole axis settled the form.
function passivitysweep(m::ResidueForm, level::Real, deadline::Real; halvings::Integer = 3)
    if isempty(m.poles)
        inside = opnorm(m.D) <= level
        return (verdict = inside ? :passive : :active,
            bands = inside ? Tuple{Float64,Float64}[] : [(0.0, Inf)], peaks = inside ? Float64[] : [Inf], evaluated = 1)
    end
    work = SweepWork(m)
    everywhereunder(m, work, level) &&
        return (verdict = :passive, bands = Tuple{Float64,Float64}[], peaks = Float64[], evaluated = 0)
    deficit = (level - 1)*(level + 1)
    n = size(m.D, 1)
    # an interval's outcome, and at a band the largest singular value at
    # its centre
    function classify(w, (a, b, tail))
        c, h = (a + b)/2, (b - a)/2
        M, roundoff, magnitude, e = expand!(w, m, c, h, tail)
        shift, allowance = testshift(w, h, M, roundoff, magnitude, level)
        settled(w, h, shift) && return (:settled, 0.0)
        short = h^2*M/4.0^halvings
        # A centre within the allowance of the level, where the remainder
        # over the interval is within it too, halving cannot settle: the
        # interval belongs to a band, as it does where the centre stands
        # above the level
        if deficit - allowance > 2*secondguard(n, deficit)
            if !positiveby!(w.W, w.P0, w.G, 0.0, deficit - max(short, allowance))
                positiveby!(w.W, w.P0, w.G, 0.0, deficit) || return (:violating, opnorm(w.S0))
                short <= allowance && return (:unresolved, opnorm(w.S0))
                secondorder!(w, m, c, h, tail, deficit, allowance, e, short) && return (:settled, 0.0)
            end
        elseif !positiveby!(w.W, w.P0, w.G, 0.0, deficit - allowance)
            positiveby!(w.W, w.P0, w.G, 0.0, deficit) || return (:violating, opnorm(w.S0))
            return (:unresolved, opnorm(w.S0))
        end
        unsettleable(w, c, h, M, allowance) && return (:unresolved, opnorm(w.S0))
        return (:halved, 0.0)
    end
    intervals = sweepintervals(m.poles)
    works = [work]
    outcomes = Tuple{Symbol,Float64}[]
    found = NTuple{4,Float64}[]
    violating, unresolved, evaluated = false, false, 0
    while !isempty(intervals)
        resize!(outcomes, length(intervals))
        fill!(outcomes, (:pending, 0.0))
        eachinterval!(classify, outcomes, intervals, works, deadline)
        halves = similar(intervals, 0)
        for ((a, b, tail), (kind, value)) in zip(intervals, outcomes)
            c = (a + b)/2
            band = tail ? (1/b, a > 0 ? 1/a : Inf) : (a, b)
            kind === :pending && return (verdict = :budget, bands = [band], peaks = [tail ? 1/c : c],
                evaluated = evaluated + count(o -> first(o) !== :pending, outcomes))
            if kind === :halved
                push!(halves, (a, c, tail), (c, b, tail))
            elseif kind !== :settled
                violating |= kind === :violating
                unresolved |= kind === :unresolved
                # a band, with its centre and the largest singular value there
                push!(found, (band..., tail ? (c > 0 ? 1/c : Inf) : c, value))
            end
        end
        evaluated += length(intervals)
        intervals = halves
    end
    verdict = violating ? :active : unresolved ? :indeterminate : :passive
    bands, peaks = mergebands(found)
    return (; verdict, bands, peaks, evaluated)
end

# The work of a sweep of a form (see passivitysweep): the form's
# coefficients in its real basis, a column of `n^2` for each basis
# function, and their Frobenius norms, which bound the roundoff of the
# response; the form's poles, a real one's exactly real, and at each a
# bound on the norm of the dissipation's residue, which the pole's negative
# shares; the response and its derivative at a centre, with
# `P0 = I - S0' S0` and `G = S1' S0`, the dissipation's derivative being
# `-(G + G')`; and the response's second derivative, and the arrays in
# which the expansion to second order is turned and tested (see
# secondorder!).
struct SweepWork
    coefficients::Matrix{Float64}
    xfrob::Vector{Float64}
    dnorm::Float64
    poles::Vector{ComplexF64}
    rho::Vector{Float64}
    basis::Matrix{Float64}
    terms::Matrix{Float64}
    S0::Matrix{ComplexF64}
    S1::Matrix{ComplexF64}
    P0::Matrix{ComplexF64}
    G::Matrix{ComplexF64}
    W::Matrix{ComplexF64}
    rows::Vector{Float64}
    basis2::Matrix{Float64}
    terms2::Matrix{Float64}
    S2::Matrix{ComplexF64}
    R0::Matrix{ComplexF64}
    R1::Matrix{ComplexF64}
    R2::Matrix{ComplexF64}
    H::Matrix{ComplexF64}
    T::Matrix{ComplexF64}
end
# The residue of the dissipation at `p_k` is `A_k = -S(-p_k)^T R_k`, `R_k`
# the form's, and at `-p_k` it is `-A_k^T`, since the dissipation is its own
# `Phi(-s)^T`, so the two have one norm. `S(-p_k)` is summed in the form's
# real basis at the point `-p_k` (see formatpoint!) to within `(N + 10) eps`
# of its scale in the Frobenius norm, and its product with `R_k` adds
# `(n + 2) eps` of the scale times `||R_k||_F`; twice their sum, added to the
# norm of the residue computed, bounds the norm of the residue.
function SweepWork(m::ResidueForm)
    n, N = size(m.D, 1), length(m.poles)
    poles = ComplexF64[isrealpole(a) ? complex(real(a)) : a for a in m.poles]
    xfrob = [norm(view(m.X, :, :, c)) for c in 1:N]
    rho = zeros(N)
    S, R, A = zeros(ComplexF64, n, n), zeros(ComplexF64, n, n), zeros(ComplexF64, n, n)
    p = 1
    while p <= N
        a = poles[p]
        pair = !isrealpole(a)
        scale = formatpoint!(S, m, poles, xfrob, -a)
        R .= pair ? complex.(view(m.X, :, :, p), view(m.X, :, :, p + 1)) : view(m.X, :, :, p)
        mul!(A, transpose(S), R)
        bound = opnorm(A) + 2*(n + N + 12)*eps()*scale*norm(R)
        rho[p] = bound
        pair && (rho[p + 1] = bound)
        p += pair ? 2 : 1
    end
    square() = zeros(ComplexF64, n, n)
    return SweepWork(reshape(m.X, n*n, N), xfrob, norm(m.D), poles, rho, zeros(N, 4), zeros(n*n, 4),
        square(), square(), square(), square(), square(), zeros(n), zeros(N, 2), zeros(n*n, 2),
        square(), square(), square(), square(), square(), square())
end
# a work for another task: the form's arrays and bounds shared, the
# arrays of a centre its own
function SweepWork(work::SweepWork)
    n, N = size(work.S0, 1), length(work.poles)
    square() = zeros(ComplexF64, n, n)
    return SweepWork(work.coefficients, work.xfrob, work.dnorm, work.poles, work.rho, zeros(N, 4), zeros(n*n, 4),
        square(), square(), square(), square(), square(), zeros(n), zeros(N, 2), zeros(n*n, 2),
        square(), square(), square(), square(), square(), square())
end

# The outcome of each interval, `outcomes[k] = f(work, intervals[k])`, the
# intervals shared among the threads where the form's order pays for it
# (see sharedwork), `blocks` blocks of them a task, taken in turn
# so that a task whose intervals cost more is made up for by the others;
# each task has a work of its own and one BLAS thread. An interval not
# reached before the clock, `time_ns()`, reaches `deadline` keeps the
# outcome it had; the clock is read every `clockevery` intervals, since a
# read costs as much as a small interval's arithmetic. The outcome of an
# interval depends on that interval alone and on the threads BLAS runs
# its products on, so that the outcomes are the same however many tasks
# share them; a single task leaves BLAS the threads it has.
function eachinterval!(f, outcomes::Vector, intervals::Vector{Tuple{Float64,Float64,Bool}},
        works::Vector{SweepWork}, deadline::Real; blocks::Integer = 4, clockevery::Integer = 16)
    K = length(intervals)
    tasks = 8*size(first(works).S0, 1)^3 >= sharedwork ? min(Threads.nthreads(), K) : 1
    if tasks == 1
        w = first(works)
        for k in 1:K
            (k - 1) % clockevery == 0 && time_ns() >= deadline && break
            outcomes[k] = f(w, intervals[k])
        end
        return outcomes
    end
    while length(works) < tasks
        push!(works, SweepWork(first(works)))
    end
    width = cld(K, blocks*tasks)
    next = Threads.Atomic{Int}(1)
    function claim(w::SweepWork)
        while true
            start = Threads.atomic_add!(next, width)
            start > K && return nothing
            for k in start:min(start + width - 1, K)
                (k - start) % clockevery == 0 && time_ns() >= deadline && return nothing
                outcomes[k] = f(w, intervals[k])
            end
        end
    end
    withoneblasthread() do
        @sync for j in 1:tasks
            Threads.@spawn claim(works[j])
        end
    end
    return outcomes
end

# The arithmetic of a piece of work, in floating point operations, from
# which the pieces are shared among the threads: under it a piece's
# products and factorizations are so small that the BLAS library's calls
# from several threads at once cost more than the threads save. An
# interval of the sweep of a form of order `n` takes products of `8 n^3`,
# so a form's intervals are shared from the order 16 (see eachinterval!),
# and a least squares of `m` rows and `c` columns a factorization of
# `2 m c^2` (see sharedchunks).
const sharedwork = 2^15

# `S` the form at the point `s` of the complex scaled frequency, summed in
# its real basis; returns the scale of its roundoff,
# `||D||_F + sum_c ||X_c||_F |f_c|`, `f_c` the basis functions (see
# pairbasis)
function formatpoint!(S::Matrix{ComplexF64}, m::ResidueForm, poles::Vector{ComplexF64},
        xfrob::Vector{Float64}, s::ComplexF64)
    S .= m.D
    scale = norm(m.D)
    p = 1
    while p <= length(poles)
        a = poles[p]
        if isrealpole(a)
            u = 1/(s - a)
            S .+= view(m.X, :, :, p) .* u
            scale += xfrob[p]*abs(u)
            p += 1
        else
            f1, f2 = pairbasis(a, s)
            S .+= view(m.X, :, :, p) .* f1 .+ view(m.X, :, :, p + 1) .* f2
            scale += xfrob[p]*abs(f1) + xfrob[p + 1]*abs(f2)
            p += 2
        end
    end
    return scale
end

# Whether the form stands under `level` at every frequency by a bound of
# its dissipation over the whole axis. Its residues at the form's poles
# `p_k` are `-E_k` and those at `-p_k` their negative transposes (see
# SweepWork), so on the axis the dissipation is `C - F - F'`, with
# `C = I - D^T D` and `F = sum_k E_k/(i x - p_k)`, and `|i x - p_k|` is at least
# `-Re p_k`: its least eigenvalue is at least that of `C` less
# `2 sum_k rho_k/(-Re p_k)`. The form is under the level where `C`, shifted by
# `level^2 - 1` less that sum, is positive definite, which a Cholesky
# factorization decides; the shift is less the roundoff of `C`,
# `2 (n + 2) eps (||D||_F^2 + ||C||_F)`, and that of the factorization of a
# matrix under `I` plus the shift (see testshift). A lossless form, whose
# residues vanish, is settled here without an interval.
function everywhereunder(m::ResidueForm, work::SweepWork, level::Real)
    n, N = size(m.D, 1), length(m.poles)
    fbound = 0.0
    for k in 1:N
        fbound += work.rho[k]/(-real(work.poles[k]))
    end
    C = Matrix{Float64}(I, n, n)
    mul!(C, transpose(m.D), m.D, -1, 1)
    shift = (level - 1)*(level + 1) - 2*fbound*(1 + 4*(N + 2)*eps()) - 2*(n + 2)*eps()*(norm(m.D)^2 + norm(C))
    shift -= (n + 2)^2*eps()*(1 + abs(shift))
    for i in 1:n
        C[i, i] += shift
    end
    return issuccess(cholesky!(Symmetric(C); check = false))
end

# The first intervals of a sweep, as `(a, b, tail)`: the poles' frequencies
# and a few half widths either side, where a resonance changes fastest,
# and eight points a decade from a tenth of the smallest pole magnitude to
# four times the largest; then, in `y = 1/x`, the interval from there to
# infinity.
function sweepintervals(poles::Vector{ComplexF64})
    hi = 4*maximum(abs, poles)
    lo = minimum(abs, poles)/10
    edges = [0.0, hi]
    for l in distinctpoles(poles)
        imag(l) > 0 || continue
        for s in (-3, -1, 0, 1, 3)
            x = imag(l) + s*abs(real(l))
            0 < x < hi && push!(edges, x)
        end
    end
    append!(edges, exp.(range(log(lo), log(hi); length = max(2, ceil(Int, 8*log10(hi/lo))))))
    sort!(unique!(edges))
    stack = Tuple{Float64,Float64,Bool}[(edges[k], edges[k + 1], false) for k in 1:length(edges) - 1]
    push!(stack, (0.0, 1/hi, true))
    return stack
end

# `S0` and `S1` of the work the response and its derivative at the centre
# `c` of an interval of half width `h`, of `x` or in the tail of `y`, and
# from them the upper triangle of `P0` and `G`; returns the remainder's
# coefficient, the roundoff of the dissipation at the ends, a bound on the
# largest eigenvalue there, and the response's roundoff, a few units of it
# for each term summed, which bounds `||dS0|| + h ||dS1||` and which the
# products carry to the dissipation with their own.
function expand!(work::SweepWork, m::ResidueForm, c::Float64, h::Float64, tail::Bool)
    n, N = size(m.D, 1), length(m.poles)
    scale, M = sweepbasis!(work, c, h, tail)
    mul!(work.terms, work.coefficients, work.basis)
    (; S0, S1, P0, G, terms) = work
    @inbounds for j in 1:n, i in 1:n
        k = i + n*(j - 1)
        S0[i, j] = complex(m.D[i, j] + terms[k, 1], terms[k, 2])
        S1[i, j] = complex(terms[k, 3], terms[k, 4])
    end
    BLAS.herk!('U', 'C', -1.0, S0, 0.0, P0)
    for i in 1:n
        P0[i, i] += 1
    end
    mul!(G, S1', S0)
    e = 2*(N + 8)*eps()*(work.dnorm + scale)
    s0, s1 = norm(S0) + e, h*norm(S1) + e
    roundoff = (2*(s0 + s1) + e)*e + 2*(n + 3)*eps()*(s0 + s1)^2
    return M*(1 + 8*(N + 1)*eps()), roundoff, 1 + 2*s0*s1, e
end

# The shift of an interval's test at `level`, `level^2 - 1` less the
# remainder and the roundoff, and the roundoff: that of the dissipation,
# and that of a Cholesky factorization of a Hermitian matrix of order `n`,
# which, where it succeeds, shows the matrix positive definite less
# `(n + 2)^2 eps` of its norm (Higham, Accuracy and Stability of Numerical
# Algorithms, 2nd ed., theorem 10.5), the matrix then positive semidefinite
# with its norm its largest eigenvalue, within `magnitude` and the shift
function testshift(work::SweepWork, h::Float64, M::Float64, roundoff::Float64, magnitude::Float64, level::Real)
    deficit = (level - 1)*(level + 1)
    allowance = roundoff + (size(work.S0, 1) + 2)^2*eps()*(magnitude + abs(deficit) + h^2*M + 2*roundoff)
    return deficit - h^2*M - allowance, allowance
end

# whether the interval is settled: the dissipation's affine part at each
# end, shifted by `shift`, positive definite
settled(work::SweepWork, h::Float64, shift::Float64) =
    positiveby!(work.W, work.P0, work.G, -h, shift) && positiveby!(work.W, work.P0, work.G, h, shift)

# Whether `P0 - t (G + G')`, shifted by `shift`, is positive definite, by a
# Cholesky factorization of its upper triangle in `W`, `P0` Hermitian and
# given by its upper triangle
function positiveby!(W::Matrix{ComplexF64}, P0::Matrix{ComplexF64}, G::Matrix{ComplexF64}, t::Float64,
        shift::Float64)
    n = size(W, 1)
    @inbounds for j in 1:n
        for i in 1:j - 1
            W[i, j] = P0[i, j] - t*(G[i, j] + conj(G[j, i]))
        end
        W[j, j] = real(P0[j, j]) - 2t*real(G[j, j]) + shift
    end
    return issuccess(cholesky!(Hermitian(W, :U); check = false))
end

# whether halving cannot settle the interval: its bound's own width is
# within the roundoff, or the interval is at the floating-point resolution
# of its centre, where its halves would be no narrower
unsettleable(work::SweepWork, c::Float64, h::Float64, M::Float64, allowance::Float64) =
    2h*holderbound(work.G, work.rows) + h^2*M <= 4*allowance || h <= eps(c)

# The form's real basis and its derivative at the centre `c` of an
# interval of half width `h`, of `x` or in the tail of `y = 1/x` (see
# passivitysweep), as the columns of `basis`: the real and imaginary parts
# of the basis, then of its derivative. Returns the scale of the terms,
# `sum ||X_c||_F (|f_c| + h |f_c'|)` over the basis functions, the
# derivative of a pair's first function counting `|u'| + |v'|` in place of
# its own size (see pairterms); and the remainder's coefficient over the
# dissipation's poles, in which the negative of a pole, its mirror across
# the axis, which lies as far from every point of it, counts as the pole
# does.
function sweepbasis!(work::SweepWork, c::Float64, h::Float64, tail::Bool)
    (; basis, poles, rho, xfrob) = work
    scale, M = 0.0, 0.0
    p = 1
    @inbounds while p <= length(poles)
        a = poles[p]
        if isrealpole(a)
            u, du, r = sweepterm(a, c, h, tail)
            basis[p, 1], basis[p, 2], basis[p, 3], basis[p, 4] = real(u), imag(u), real(du), imag(du)
            scale += xfrob[p]*(abs(u) + h*abs(du))
            M += 2*rho[p]*r
            p += 1
        else
            # the pair's functions `u + v` and `i (u - v)` and their
            # derivatives, `v` its conjugate pole's term (see pairterms)
            f1, f2, d1, d2, dsize, r = pairterms(a, c, h, tail)
            basis[p, 1], basis[p, 2], basis[p, 3], basis[p, 4] = real(f1), imag(f1), real(d1), imag(d1)
            basis[p + 1, 1], basis[p + 1, 2], basis[p + 1, 3], basis[p + 1, 4] = real(f2), imag(f2), real(d2), imag(d2)
            scale += xfrob[p]*(abs(f1) + h*dsize) + xfrob[p + 1]*(abs(f2) + h*abs(d2))
            M += 2*rho[p]*r
            p += 2
        end
    end
    return scale, M
end

# a pole's term at the centre, its derivative, and its remainder's
# coefficient over the interval, in `x` or in the tail's `y`
function sweepterm(a::ComplexF64, c::Float64, h::Float64, tail::Bool)
    if tail
        b = im - a*c
        q = im/a
        return c/b, im/b^2, 1/(abs2(b)*sqrt(imag(q)^2 + max(0.0, abs(c - real(q)) - h)^2))
    end
    s = im*c - a
    return 1/s, -im/s^2, 1/(abs2(s)*sqrt(real(a)^2 + max(0.0, abs(c - imag(a)) - h)^2))
end

# A pole pair's two functions at the centre `c` and their derivatives, in
# `x` or in the tail's `y`, formed from the differences of the point and
# the two poles without the difference of the pole's term and its
# conjugate's (see pairbasis), with the pair's remainder coefficient over
# the interval (see sweepterm). In `x`, with `t = i c - Re a` and
# `d = (i c - a)(i c - conj(a))`, the functions are `2t/d`, `-2 Im a/d`,
# `2i ((Im a)^2 - t^2)/d^2` and `4i Im a t/d^2`; in `y`, with
# `e = (i - a c)(i - conj(a) c)`, `2c (i - c Re a)/e`, `-2 Im a c^2/e`,
# `2 (2 c Re a - i (1 + ((Im a)^2 - (Re a)^2) c^2))/e^2` and
# `4 Im a c (1 + i c Re a)/e^2`. The numerator of the first derivative is a
# sum whose terms are as large as `|u'| + |v'|` times the denominator, and
# carries the roundoff of those, the size returned with it; the others
# carry a few units of their own.
function pairterms(a::ComplexF64, c::Float64, h::Float64, tail::Bool)
    ar, ai = real(a), imag(a)
    if tail
        za, zb = im - a*c, im - conj(a)*c
        na, nb = abs2(za), abs2(zb)
        qa, qb = im/a, im/conj(a)
        r = 1/(na*sqrt(imag(qa)^2 + max(0.0, abs(c - real(qa)) - h)^2)) +
            1/(nb*sqrt(imag(qb)^2 + max(0.0, abs(c - real(qb)) - h)^2))
        w = 1/(za*zb)
        return 2c*(im - ar*c)*w, -2ai*c^2*w, 2*(2ar*c - im*(1 + (ai^2 - ar^2)*c^2))*w^2,
            4ai*c*(1 + im*ar*c)*w^2, 1/na + 1/nb, r
    end
    za, zb = im*c - a, im*c - conj(a)
    na, nb = abs2(za), abs2(zb)
    r = 1/(na*sqrt(ar^2 + max(0.0, abs(c - ai) - h)^2)) + 1/(nb*sqrt(ar^2 + max(0.0, abs(c + ai) - h)^2))
    t = im*c - ar
    w = 1/(za*zb)
    return 2t*w, -2ai*w, 2im*(ai^2 - t^2)*w^2, 4im*ai*t*w^2, 1/na + 1/nb, r
end

# `S2` of the work the response's second derivative at the centre `c`, in
# `x` or in the tail's `y`, summed in the form's real basis as `S0` and `S1`
# are; returns its roundoff, a few units of `sum ||X_c||_F |f_c''|` over
# the basis, a function's size being that of the terms it sums. In `x` a
# pole pair's second derivatives are `-4 t (t^2 - 3 (Im a)^2)/d^3` and
# `4 Im a (3 t^2 - (Im a)^2)/d^3`, with `t = i c - Re a` and
# `d = (i c - a)(i c - conj(a))` (see pairterms); in the tail, where no pole
# is near, they are the sum and the difference of the pole's term and its
# conjugate's.
function secondderivative!(work::SweepWork, m::ResidueForm, c::Float64, tail::Bool)
    (; basis2, poles, xfrob, S2, terms2) = work
    n, N = size(S2, 1), length(poles)
    scale = 0.0
    p = 1
    @inbounds while p <= N
        a = poles[p]
        if isrealpole(a)
            u = tail ? 2im*a/(im - a*c)^3 : -2/(im*c - a)^3
            basis2[p, 1], basis2[p, 2] = real(u), imag(u)
            scale += xfrob[p]*abs(u)
            p += 1
        else
            if tail
                u, v = 2im*a/(im - a*c)^3, 2im*conj(a)/(im - conj(a)*c)^3
                f1, f2 = u + v, im*(u - v)
                size1 = size2 = abs(u) + abs(v)
            else
                t, b = im*c - real(a), imag(a)
                d = (im*c - a)*(im*c - conj(a))
                f1, f2 = -4t*(t^2 - 3b^2)/d^3, 4b*(3t^2 - b^2)/d^3
                size1, size2 = 4abs(t)*(abs2(t) + 3b^2)/abs(d)^3, 4abs(b)*(3abs2(t) + b^2)/abs(d)^3
            end
            basis2[p, 1], basis2[p, 2] = real(f1), imag(f1)
            basis2[p + 1, 1], basis2[p + 1, 2] = real(f2), imag(f2)
            scale += xfrob[p]*size1 + xfrob[p + 1]*size2
            p += 2
        end
    end
    mul!(terms2, work.coefficients, basis2)
    @inbounds for j in 1:n, i in 1:n
        k = i + n*(j - 1)
        S2[i, j] = complex(terms2[k, 1], terms2[k, 2])
    end
    return 2*(N + 8)*eps()*scale
end

# The coefficient of the dissipation's remainder past its second order over
# an interval, `|t|^3` times it bounding the remainder: a term
# `A_j/(a_j + i t)` leaves `|t|^3 ||A_j||/(|a_j|^3 d_j)` in `x`, and one
# `A_j y/(i - q_j y)` of the tail `|t|^3 ||A_j|| |q_j|/(|b_j|^3 e_j)` (see
# passivitysweep), a pole's negative counting as the pole does, with the
# roundoff of the sum
function cubicremainder(work::SweepWork, c::Float64, h::Float64, tail::Bool)
    M3 = 0.0
    @inbounds for (q, r) in zip(work.poles, work.rho)
        if tail
            b, y = im - q*c, im/q
            M3 += 2r*abs(q)/(abs(b)^3*sqrt(imag(y)^2 + max(0.0, abs(c - real(y)) - h)^2))
        else
            a = im*c - q
            M3 += 2r/(abs(a)^3*sqrt(real(q)^2 + max(0.0, abs(c - imag(q)) - h)^2))
        end
    end
    return M3*(1 + 8*(length(work.poles) + 1)*eps())
end

# Whether the dissipation's expansion to second order settles an interval
# the first order leaves unsettled, its centre standing within `short` of
# the level (see passivitysweep). At `c + t` the dissipation is
# `Q(t) = P0 + t P1 + t^2 P2` and a remainder of norm at most `|t|^3 M3`
# (see cubicremainder), `P2 = -(S1' S1 + (S2' S0 + S0' S2)/2)` with `S2` the
# response's second derivative, and the interval is settled where
# `Q + gamma` is positive semidefinite over it, `gamma` being `level^2 - 1`
# less the remainder and the roundoff, for which `P0 + gamma` must be
# positive definite. Over the interval the first- and second-order terms
# move no eigenvalue by more than `eta = h ||P1|| + h^2 ||P2||`. A pivoted
# Cholesky factorization of `P0`, stopped where its pivots fall under
# `(2 eta - gamma)/n`, spans the directions it dissipates in, `W`, and
# leaves along the rest, `V`, a block whose eigenvalues stand under
# `2 eta - gamma`; where it leaves none, `P0 + gamma - eta` positive definite
# settles the interval. Otherwise, in the coordinates the reflectors of
# the smaller span give, the blocks of `P_k` are `A_k` along `V`, `B_k` from
# `V` to `W` and `C_k` along `W`. Along `W` the expansion stays above
# `C = C0 + gamma - eta_W`, `eta_W = h ||C1|| + h^2 ||C2||`, and where that is
# positive definite `Q + gamma` is positive semidefinite where the Schur
# complement of that block is,
# `K(t) = A0 + gamma + t A1 + t^2 A2 - B(t)' C^-1 B(t)` with
# `B(t) = B0 + t B1 + t^2 B2`, or `K0 + t K1 + t^2 K2 + t^3 K3 + t^4 K4`,
# whose coefficients are Hermitian products: with `C = U' U` and
# `Y_k = U^-H B_k`, `B_i' C^-1 B_j` is `Y_i' Y_j`.
# Where the form is lossless along `V` the dissipation stays null there:
# `K1` vanishes and `K2` is of the order of `h`, the second-order term
# making up for the coupling to `W` through which the first order's least
# eigenvalue dips toward the ends of the interval, so that the intervals
# such a form needs grow as the inverse cube root of the margin rather
# than its inverse square root. `K0 + t K1 + t^2 K2` stands above
# `K0 - h ||K1|| + t^2 K2`, which is affine in `t^2`, and positive definite
# at both ends of `[0, h^2]` it is so between. The norms are Holder's
# bounds (see holderbound), and `h^3 ||K3|| + h^4 ||K4||` is taken off.
# The roundoff is the first order's `allowance` (see testshift), the
# second-order term's, carried from the response's roundoff `e` (see
# expand!) and the second derivative's, with `||S0||` at most the level
# at a centre under it, and `(n + 2)^2 eps` of the expansion's size for
# the reflectors, the products and the factorizations; the part of it
# which does not shrink with the interval is secondguard's. The centre
# stands within `short` of the level, so that the test can settle the
# interval only where the remainder fits under `short` less the roundoff
# that does not shrink with it. The coefficients are turned, and the
# complement formed, in the work's arrays.
function secondorder!(work::SweepWork, m::ResidueForm, c::Float64, h::Float64, tail::Bool, deficit::Real,
        allowance::Float64, e::Float64, short::Float64)
    (; S0, S1, S2, G, W, P0, R0, R1, R2, H, T, rows) = work
    n = size(S0, 1)
    room = deficit - allowance - secondguard(n, deficit)
    M3 = cubicremainder(work, c, h, tail)
    h^3*M3 < short + room - deficit || return false
    e2 = secondderivative!(work, m, c, tail)
    s0, s1, s2 = sqrt(1 + deficit), h*norm(S1), h^2*norm(S2)
    size2 = 2*s0*s1 + s1^2 + s0*s2
    gamma = room - h^3*M3 - (2*s1 + s2)*e - s0*h^2*e2 - (n + 2)^2*eps()*(size2 + h^3*M3)
    positiveby!(W, P0, G, 0.0, gamma) || return false
    # `P0`, `P1 = -(G + G')` and `P2` in full, in `R0`, `R1` and `R2`
    mul!(R2, S2', S0)
    @inbounds for j in 1:n, i in 1:j
        R0[i, j], R0[j, i] = P0[i, j], conj(P0[i, j])
        x = -(G[i, j] + conj(G[j, i]))
        R1[i, j], R1[j, i] = x, conj(x)
        x = -(R2[i, j] + conj(R2[j, i]))/2
        R2[i, j], R2[j, i] = x, conj(x)
    end
    mul!(R2, S1', S1, -1, 1)
    eta = h*holderbound(R1, rows) + h^2*holderbound(R2, rows)
    # the pivoted factor `U = [U1 U2]`: its rows span `W`, and `V` is spanned
    # by the solutions of `U z = 0`, `z = [-U1^-1 U2; I] y`, in its order; the
    # smaller span's basis goes to the first columns of `H`
    copyto!(W, P0)
    _, piv, q, _ = LAPACK.pstrf!('U', W, (2eta - gamma)/n)
    q == n && return positiveby!(W, P0, G, 0.0, gamma - eta)
    r = n - q
    Y = view(H, :, 1:min(q, r))
    fill!(Y, 0)
    if q <= r
        for k in 1:q, i in k:n
            Y[piv[i], k] = conj(W[k, i])
        end
        w, v = 1:q, q + 1:n
    else
        Z = view(W, 1:q, q + 1:n)
        BLAS.trsm!('L', 'U', 'N', 'N', one(ComplexF64), view(W, 1:q, 1:q), Z)
        for k in 1:r
            Y[piv[q + k], k] = 1
            for i in 1:q
                Y[piv[i], k] = -Z[i, k]
            end
        end
        w, v = r + 1:n, 1:r
    end
    # its QR factorization's reflectors, `I - Y F Y'` the unitary which turns
    # the coefficients
    if q > 0
        F = view(T, 1:size(Y, 2), 1:size(Y, 2))
        LAPACK.geqrt!(Y, F)
        for k in axes(Y, 2)
            Y[1:k - 1, k] .= 0
            Y[k, k] = 1
        end
        for R in (R0, R1, R2)
            turn!(R, Y, F, W)
        end
    end
    etaw = h*holderbound(view(R1, w, w), rows) + h^2*holderbound(view(R2, w, w), rows)
    C = view(R0, w, w)
    for i in axes(C, 1)
        C[i, i] += gamma - etaw
    end
    issuccess(cholesky!(Hermitian(C, :U); check = false)) || return false
    # `B_k` turned into `Y_k = U^-H B_k`, `C = U' U`, and the upper triangles
    # of `K0`, `K1` and `K2` formed in those of the blocks along `V`, those of
    # `-K3` and `-K4` in `W`
    B0, B1, B2 = view(R0, w, v), view(R1, w, v), view(R2, w, v)
    for B in (B0, B1, B2)
        BLAS.trsm!('L', 'U', 'C', 'N', one(ComplexF64), C, B)
    end
    K0, K1, K2 = view(R0, v, v), view(R1, v, v), view(R2, v, v)
    for i in 1:r
        K0[i, i] += gamma
    end
    BLAS.herk!('U', 'C', -1.0, B0, 1.0, K0)
    BLAS.her2k!('U', 'C', -one(ComplexF64), B0, B1, 1.0, K1)
    BLAS.herk!('U', 'C', -1.0, B1, 1.0, K2)
    BLAS.her2k!('U', 'C', -one(ComplexF64), B0, B2, 1.0, K2)
    K = view(W, 1:r, 1:r)
    BLAS.her2k!('U', 'C', one(ComplexF64), B1, B2, 0.0, K)
    cubic = h^3*holderbound(Hermitian(K, :U), rows)
    BLAS.herk!('U', 'C', 1.0, B2, 0.0, K)
    cubic += h^4*holderbound(Hermitian(K, :U), rows)
    # `K0 - delta` and `K0 + h^2 K2 - delta`, the latter in `W`
    delta = h*holderbound(Hermitian(K1, :U), rows) + cubic
    @inbounds for j in 1:r, i in 1:j
        K[i, j] = K0[i, j] + h^2*K2[i, j]
    end
    for i in 1:r
        K0[i, i] -= delta
        K[i, i] -= delta
    end
    return issuccess(cholesky!(Hermitian(K0, :U); check = false)) &&
        issuccess(cholesky!(Hermitian(K, :U); check = false))
end

# `R` turned in place to `Q' R Q`, `Q = I - Y F Y'` the unitary of a QR
# factorization's reflectors `Y` and triangular factor `F` (see
# LAPACK.geqrt!), `X` the products' scratch
function turn!(R::Matrix{ComplexF64}, Y::AbstractMatrix{ComplexF64}, F::AbstractMatrix{ComplexF64},
        X::Matrix{ComplexF64})
    k = size(Y, 2)
    A = view(X, :, 1:k)
    mul!(A, R, Y)
    BLAS.trmm!('R', 'U', 'N', 'N', one(ComplexF64), F, A)
    mul!(R, A, Y', -1, 1)
    B = view(X, 1:k, :)
    mul!(B, Y', R)
    BLAS.trmm!('L', 'U', 'C', 'N', one(ComplexF64), F, B)
    mul!(R, Y, B, -1, 1)
    return R
end

# the roundoff of the second-order test which does not shrink with the
# interval, `(n + 2)^2 eps` of the dissipation's size at a centre under the
# level (see secondorder!); the test is tried only where the margin the
# first order's roundoff leaves exceeds it twice, since past that the
# margin it works in is mostly its own roundoff
secondguard(n::Integer, deficit::Real) = (n + 2)^2*eps()*(2 + 2*abs(deficit))

# The largest singular value of the form over every frequency, by the
# bounds of passivitysweep: the response at each interval's centre raises
# a lower bound, from the constant and the response at the poles'
# frequencies; a form the bound over the whole axis puts under that lower
# bound times `1 + rtol` needs no interval (see everywhereunder); an
# interval whose bound is under the lower bound times `1 + rtol` is
# settled, as is one the expansion to second order settles where the
# centre stands near that level and leaves it the room of its roundoff
# (see passivitysweep), and one halving cannot settle, its centre within
# the roundoff of that level where the remainder over it is too, or its
# bound's own width within the roundoff, whose bound then counts toward
# the ceiling; and any other is halved. The intervals go a generation at
# a time, as in passivitysweep, each measured against the lower bound as
# the generations before left it and its own centre, so that the search
# is the same however many threads share them (see eachinterval!).
# Returns the lower bound, the frequency in rad/s where it was attained,
# and the ceiling no frequency reaches, `Inf` where the clock,
# `time_ns()`, reached `deadline` first.
function residuenorm(m::ResidueForm; rtol::Real, deadline::Real = Inf, halvings::Integer = 3)
    n = size(m.D, 1)
    bound, where = opnorm(m.D), Inf
    isempty(m.poles) && return bound, where, bound
    work = SweepWork(m)
    for l in distinctpoles(m.poles)
        imag(l) > 0 || continue
        expand!(work, m, imag(l), 0.0, false)
        s = opnorm(work.S0)
        s > bound && ((bound, where) = (s, imag(l)))
    end
    everywhereunder(m, work, bound*(1 + rtol)) && return bound, where*m.wref, bound*(1 + rtol)
    intervals = sweepintervals(m.poles)
    works = [work]
    outcomes = Tuple{Symbol,Float64,Float64}[]
    ceiling = 0.0
    while !isempty(intervals)
        resize!(outcomes, length(intervals))
        fill!(outcomes, (:pending, 0.0, 0.0))
        # an interval's outcome, the largest singular value at its centre
        # where that stands above the lower bound `below`, and its bound
        # where halving cannot settle it
        let below = bound
            eachinterval!(outcomes, intervals, works, deadline) do w, (a, b, tail)
                c, h = (a + b)/2, (b - a)/2
                M, roundoff, magnitude, e = expand!(w, m, c, h, tail)
                s = positiveby!(w.W, w.P0, w.G, 0.0, (below - 1)*(below + 1)) ? 0.0 : opnorm(w.S0)
                top = max(below, s)
                shift, allowance = testshift(w, h, M, roundoff, magnitude, top*(1 + rtol))
                settled(w, h, shift) && return (:settled, s, 0.0)
                deficit, short = (top*(1 + rtol) - 1)*(top*(1 + rtol) + 1), h^2*M/4.0^halvings
                # the interval's bound, which counts toward the ceiling where
                # halving cannot settle the interval
                ceilingof() = sqrt(top^2 + 2h*holderbound(w.G, w.rows) + h^2*M + 2*allowance)
                if !positiveby!(w.W, w.P0, w.G, 0.0, deficit - max(short, allowance))
                    # A centre within the allowance of the level, where the
                    # remainder is within it too or the roundoff leaves the
                    # second order no room, halving cannot settle, as in the
                    # sweep; the second order is tried where the centre
                    # leaves it the room its roundoff takes, whatever the
                    # level, which for a form under one stands under one.
                    if !positiveby!(w.W, w.P0, w.G, 0.0, deficit - allowance)
                        (short <= allowance || deficit - allowance <= 2*secondguard(n, deficit)) &&
                            return (:unresolved, s, ceilingof())
                    elseif positiveby!(w.W, w.P0, w.G, 0.0, deficit - allowance - 2*secondguard(n, deficit))
                        secondorder!(w, m, c, h, tail, deficit, allowance, e, short) && return (:settled, s, 0.0)
                    end
                end
                unsettleable(w, c, h, M, allowance) && return (:unresolved, s, ceilingof())
                return (:halved, s, 0.0)
            end
        end
        halves = similar(intervals, 0)
        for ((a, b, tail), (kind, s, value)) in zip(intervals, outcomes)
            kind === :pending && return bound, where*m.wref, Inf
            c = (a + b)/2
            s > bound && ((bound, where) = (s, tail ? (c > 0 ? 1/c : Inf) : c))
            kind === :unresolved && (ceiling = max(ceiling, value))
            kind === :halved && push!(halves, (a, c, tail), (c, b, tail))
        end
        intervals = halves
    end
    return bound, where*m.wref, max(ceiling, bound*(1 + rtol))
end

# The bands `(a, b, centre, value)` merged where they meet or overlap, and
# of each merged band the centre with the largest value
function mergebands(found::Vector{NTuple{4,Float64}})
    bands, peaks, values = Tuple{Float64,Float64}[], Float64[], Float64[]
    for (a, b, x, v) in sort(found)
        if !isempty(bands) && a <= last(bands)[2]
            bands[end] = (last(bands)[1], max(last(bands)[2], b))
            v > values[end] && ((peaks[end], values[end]) = (x, v))
        else
            push!(bands, (a, b)); push!(peaks, x); push!(values, v)
        end
    end
    return bands, peaks
end

# each run of consecutive points of `grid` whose `hits` are true, from the
# point before it to the one after, as a band
function hitbands(grid::AbstractVector, hits::AbstractVector{Bool})
    bands = Tuple{Float64,Float64}[]
    k = 1
    while k <= length(grid)
        if hits[k]
            j = k
            while j < length(grid) && hits[j + 1]
                j += 1
            end
            push!(bands, (grid[max(k - 1, 1)], grid[min(j + 1, length(grid))]))
            k = j + 1
        else
            k += 1
        end
    end
    return bands
end

# The passivity enforcement of a fit given as its poles and residues in
# rad/s and its constant, at the sample angular frequencies `ws`. The
# bands where the largest singular value stands above the level,
# `sqrt(1 + atol)`, are found on a grid, or by the certified sweep where
# the grid finds none (see passivitysweep); the residues and the constant
# are changed by the least amount, in the norm of the fit over the
# samples, that brings the worst point of each band to one less a margin
# -- a linear constraint through the singular vectors -- and the rounds
# repeat, with a bound on their number, until none is found. A certified
# sweep which settles every frequency under the level establishes the
# fit passive to it; otherwise the ceiling of the norm search decides,
# and a fit whose ceiling stands above the level is contracted until it
# stands under one, or toward a value stated at zero frequency under the
# level. The least change is solved in the residues' own coordinates at
# the frequency scaled to the band, where its normal matrix is that of
# the pole basis over the samples and carries no unit.
# Returns the residues, the constant, and the contraction the end
# applied, its factor and how far it moved the response over the
# samples, or `nothing` where it moved it by no more than `atol`; the
# caller warns of the contraction of a fit it returns. `bandpoints` points
# across a band find its worst, a norm search resolves to `normrtol` of
# the norm or to half the margin between one and the level where that is
# finer, and a correction which moves a value held at zero frequency by
# more than `dcchecktol` is an error of the arithmetic.
# The residues are taken within their `spaces` (see residuespaces) from
# the start, and the residues returned lie within them, so that the
# realization through the same spaces is the block made passive.
# `reciprocal` says the samples are, and constrains the transpose of each
# point as well. `memory` is the bytes a weighted correction's metric may
# take (see weightedmetricbytes), which is refused before it is built
# where it needs more.
function enforcepassivity(poles::Vector{ComplexF64}, residues::AbstractArray{<:Complex,3},
        D::AbstractMatrix, ws::AbstractVector, passivity::PassivityEnforcement = PassivityEnforcement();
        atol = 1e-8, dc::Matrix{Float64} = nodc,
        spaces::Vector{Matrix{ComplexF64}} = residuespaces(poles, residues, 1e-12), reciprocal::Bool = false,
        weights::Array{Float64,3} = noweights, memory::Real = Inf, bandpoints::Integer = 51,
        normrtol::Real = 1e-8, dcchecktol::Real = 1e-9)
    (; margin, rounds, scalelimit, maxtime, regularization) = passivity
    # the time by which the certified sweeps and norm searches must be done
    deadline = time_ns() + 1e9*maxtime
    lo = findfirst(>(0), ws)
    isnothing(lo) && throw(ArgumentError("the enforcement needs at least one positive sample frequency."))
    wref = sqrt(ws[lo]*maximum(ws))
    xs = ws ./ wref
    m = ResidueForm(poles, residues, D, wref, spaces)
    X, D = m.X, m.D
    n = size(D, 1)
    # A residue taken within its space loses its part `P` outside it, and
    # with it `-P/a` of the value at zero frequency, `a` its pole; the
    # constant takes that back, so that a value stated there, which the
    # fit meets, is the one held from here on.
    if !isempty(dc)
        for p in eachindex(poles)
            Q = spaces[p]
            size(Q, 2) < n || continue
            R = isrealpole(poles[p]) ? complex.(real.(view(residues, :, :, p))) : residues[:, :, p]
            D .-= real.((R .- (R*Q)*Q') ./ poles[p])
        end
    end
    # The level the rounds hold the fit to: its dissipation `I - S'S` no
    # lower than `-atol` at any frequency, as the block's samples and its
    # constant term are held (see checkblockcontract and checkpassive), so
    # that a fit the rounds leave passive passes the constructor's checks
    level = passivelevel(atol)
    # A value held at zero frequency holds every fit which meets it there,
    # so one above the level leaves no passive fit to find; where a fit
    # holds one, the correction at zero is held at zero, and a constraint
    # there says nothing.
    isempty(dc) || opnorm(dc) <= level || throw(ArgumentError(
        lazy"the value the block states at zero frequency has a largest singular value of $(opnorm(dc)), above sqrt(1 + atol) for atol = $(atol): no fit which meets it is passive to atol. State a passive value, or raise atol."))
    # the norm search resolves to half the margin between one and the
    # level where `normrtol` is coarser, so that the ceiling of a fit which
    # reaches one, as a fit held at a value of unit norm at zero frequency
    # does, stands under the level
    normsearch(form) = residuenorm(form; rtol = min(normrtol, (level - 1)/2), deadline)
    blocks, nu = correctioncoordinates(m)
    nz = nu - n
    # The metric of the least change depends on the poles, the residues'
    # spaces, the samples, the weights and a value held at zero, none of
    # which the rounds change: it is built at the first correction, so that
    # a fit the rounds find passive at once never builds it, and kept.
    metricref = Ref{Union{Nothing,CorrectionMetric}}(nothing)
    phi = zeros(ComplexF64, length(m.poles))
    F = zeros(ComplexF64, n, n)
    g = zeros(ComplexF64, nu)
    corrections, stopped = 0, :rounds
    # The points of every round, infinite frequency among them, are
    # constrained again in the rounds after, linearized afresh, at their
    # largest singular value whether or not it stands above the level: a
    # change which is cheap over the samples moves the response freely
    # beside them, above the band most, and a point constrained only while
    # it stands above the level is pushed back over by the next correction.
    kept = Float64[]
    # the worst point of each band, among its `peaks` where they are
    # given, where it stands above the level, or wherever it stands with
    # `above = false`; a band reaching infinite
    # frequency is sampled over three decades from its start, and at
    # infinity, where the response is `D`. `m` holds `X` and `D`, which the
    # rounds change in place.
    bandgrid(a, b) = isinf(b) ? vcat(exp.(range(log(max(a, eps())), log(1e3*max(a, 1.0)); length = bandpoints - 1)), Inf) :
        range(a, b; length = bandpoints)
    worstpoints(bands, peaks = fill(NaN, length(bands)); above = true) = Float64[x for x in (begin
        grid = isnan(peak) ? bandgrid(a, b) : vcat(bandgrid(a, b), peak)
        grid[argmax([opnorm(responseat!(F, m, x, phi)) for x in grid])]
    end for ((a, b), peak) in zip(bands, peaks)) if !above || opnorm(responseat!(F, m, x, phi)) > level]
    certified = false
    for roundindex in 1:rounds
        bands = residueviolations(m, xs, level)
        found = worstpoints(bands)
        if isempty(found)
            # The grid can pass between the edges of a narrow band; the
            # certified sweep cannot (see passivitysweep). Where it settles
            # every frequency under the level, infinity among them, the fit
            # is passive to the level; otherwise its bands are where the fit
            # stands above it, or within the roundoff of it, and the worst
            # point of each is constrained, wherever it stands.
            sweep = passivitysweep(m, level, deadline)
            sweep.verdict === :budget && throw(ArgumentError(
                lazy"the fit's passivity could not be settled within maxtime = $(maxtime) s, the sweep of the frequency axis reaching the interval from $(first(only(sweep.bands))*wref) to $(last(only(sweep.bands))*wref) rad/s. Raise maxtime, fit with fewer poles, or over a narrower band."))
            if sweep.verdict === :passive
                certified, stopped = true, :clean
                break
            end
            bands = sweep.bands
            found = worstpoints(bands, sweep.peaks; above = false)
        end
        # infinite frequency, where no band reaches: the constant alone
        opnorm(D) > level && push!(found, Inf)
        isempty(dc) || filter!(!iszero, found)
        isempty(found) && (stopped = :clean; break)
        points = unique!(append!(kept, found))
        # one constraint per singular value above the level at each point,
        # and one on the largest whatever it is,
        # `d sigma_j = Re(u_j' dS v_j)`, which for port `i` reads
        # `x_i' Re(conj(u_ij) g_j)` with `g_j` the unknowns' response to
        # `v_j`: the same `g_j` for every port. A reciprocal block has each
        # constraint's transpose as well, so that where every residue has
        # full rank the least change is itself reciprocal; the space of a
        # residue of lower rank need not hold the transpose of its change.
        Gs = Vector{ComplexF64}[]
        As = Vector{ComplexF64}[]
        rhs = Float64[]
        for x in points
            S = responseat!(F, m, x, phi)
            Fs = svd(S)
            for j in 1:n
                (j == 1 || Fs.S[j] > level) || continue
                pairs = reciprocal ? ((Fs.U[:, j], Fs.V[:, j]), (conj.(Fs.V[:, j]), conj.(Fs.U[:, j]))) :
                    ((Fs.U[:, j], Fs.V[:, j]),)
                for (u, v) in pairs
                    coordinateresponse!(g, blocks, phi, v)
                    push!(Gs, copy(g))
                    push!(As, conj.(u))
                    push!(rhs, 1 - margin - Fs.S[j])
                end
            end
        end
        mc = length(rhs)
        # The dual of the least change: with `H = L L'` and `h = L^-1 g`,
        # the Gram matrix of the constraints over the ports is
        # `sum_i Re(a_i h)' Re(a_j h)`, built from the real and imaginary
        # parts of `a` and `h` without forming any constraint over all the
        # ports' unknowns. `h` is solved for every constraint at once, as
        # two real triangular solves: the real factor would otherwise be
        # converted to a complex one for each. The solves, and the
        # projection onto the coordinates a value held at zero leaves free,
        # are taken along the basis where the metric is factored so (see
        # CorrectionMetric).
        if isnothing(metricref[])
            need = isempty(weights) ? 0.0 : weightedmetricbytes(m.spaces, n)
            need <= memory || throw(ArgumentError(lazy"the passivity correction of this weighted fit needs a metric of about $(round(need/1e9; sigdigits = 3)) GB, more than the $(round(max(memory, 0.0)/1e9; sigdigits = 3)) GB that maxstates leaves it once the fit's states are held: each entry has a metric of its own, and a residue of less than full rank gives each of the $(n) ports one of every unknown of its row. Fit with fewer poles, without weights, or raise maxstates where the memory allows."))
            metricref[] = correctionmetric(m, blocks, nu, xs, regularization, dc, weights)
        end
        metric = metricref[]::CorrectionMetric
        (; L, Z, stride) = metric
        Gm = reduce(hcat, Gs)
        isnothing(Z) || (Gm = alongbasis(G -> transpose(Z)*G, Gm, stride))
        Am = reduce(hcat, As)
        AR, AI = real.(Am), imag.(Am)
        if length(L) == 1
            HR = alongbasis(G -> only(L) \ G, real.(Gm), stride)
            HI = alongbasis(G -> only(L) \ G, imag.(Gm), stride)
            Gd = (transpose(AR)*AR) .* (transpose(HR)*HR) .- (transpose(AR)*AI) .* (transpose(HR)*HI) .-
                 (transpose(AI)*AR) .* (transpose(HI)*HR) .+ (transpose(AI)*AI) .* (transpose(HI)*HI)
        else
            # A metric of each port's own: port `i`'s rows of the
            # constraints, `Re(a_i g)`, taken through its factor, `h_i`,
            # and the Gram matrix summed over the ports, `sum_i h_i' h_i`;
            # `h_i` is not held for every port, and the correction takes
            # the constraints' combination through the factor instead.
            GR, GI = real.(Gm), imag.(Gm)
            rowsof(i) = portsolve(metric, i, GR .* transpose(view(AR, i, :)) .- GI .* transpose(view(AI, i, :)))
            Gd = zeros(length(rhs), length(rhs))
            for i in 1:n
                h = rowsof(i)
                mul!(Gd, transpose(h), h, true, true)
            end
        end
        Gd = (Gd .+ transpose(Gd))./2
        μ, solved = dualactiveset(Symmetric(Gd), rhs)
        # multipliers which meet every condition of the dual's optimum
        # give a change which meets its constraints (see dualactiveset)
        solved || (stopped = :infeasible; break)
        # the correction, every port's unknowns a column, and what it does
        # to the constraints it was built from, `-G mu`
        Xc = if length(L) == 1
            -alongbasis(G -> transpose(only(L)) \ G, HR*(μ .* transpose(AR)) .- HI*(μ .* transpose(AI)), stride)
        else
            # port `i`'s change, `-H_i^-1 Re(a_i g) mu`, the constraints
            # combined before the solves, one column where `h_i` has one per
            # constraint
            reduce(hcat, [-portsolve(metric, i, portsolve(metric, i,
                reshape(GR*(μ .* view(AR, i, :)) .- GI*(μ .* view(AI, i, :)), :, 1)); adjoint = true) for i in 1:n])
        end
        isnothing(Z) || (Xc = alongbasis(G -> Z*G, Xc, stride))
        dX, dD = zeros(size(X)), zeros(n, n)
        for (c, cols, M) in blocks
            change = transpose(M*view(Xc, cols, :))
            c == 0 ? (dD .+= change) : (dX[:, :, c] .+= change)
        end
        # the elimination is exact, so this is a check on the arithmetic
        if !isempty(dc)
            basisat!(phi, m.poles, 0.0)
            at0 = copy(dD)
            for c in eachindex(phi)
                at0 .+= real(phi[c]) .* view(dX, :, :, c)
            end
            maximum(abs, at0; init = 0.0) <= dcchecktol*(1 + maximum(abs, Xc)) ||
                throw(ArgumentError("the passivity correction moved the stated value at zero frequency; report this."))
        end
        X .+= dX
        D .+= dD
        corrections += 1
        (all(isfinite, X) && all(isfinite, D)) || throw(ArgumentError(
            lazy"the passivity enforcement diverged in round $(roundindex) with $(nz) states: its least norm correction is not finite. Fit with fewer poles, or over a narrower band."))
    end
    certified && return scaledresidues(m) .* wref, D, nothing
    # The norm search brackets the largest singular value between what the
    # fit reaches and a ceiling no frequency reaches; only the ceiling
    # establishes the fit under the level, and a fit whose ceiling stands
    # above it is contracted by the ceiling, while what it reaches decides
    # whether it is too far above one to be contracted at all.
    worst, peak, ceiling = normsearch(m)
    isfinite(ceiling) || throw(ArgumentError(
        lazy"the fit's largest singular value could not be bounded within maxtime = $(maxtime) s, and the largest value it saw was $(worst), at $(peak) rad/s. Raise maxtime, fit with fewer poles, or over a narrower band."))
    ceiling <= level && return scaledresidues(m) .* wref, D, nothing
    worst <= 1 + scalelimit || throw(ArgumentError(
        lazy"the fit could not be made passive: $(enforcementrecord(stopped, corrections)), and its largest singular value is $(worst), at $(peak) rad/s, more than scalelimit = $(scalelimit) above one, too far to scale under it without describing a block which is mostly loss. Fit with fewer poles, or over a narrower band."))
    # The contraction along the path `S0 + t (S - S0)`, toward nothing or
    # toward the value stated at zero, measured by the ceiling; the time
    # once out, every later search is out of it as well, so the first
    # search along the path that runs out refuses the fit. A contraction
    # toward nothing brings the ceiling under one; one toward a stated
    # value brings it under the level, since the statement's own norm may
    # be one, and no point of the path then has a ceiling under one.
    S0 = isempty(dc) ? zeros(size(D)) : dc
    anchor = opnorm(S0)
    target = anchor > 0 ? level : 1.0
    Xbefore, Dbefore = copy(X), copy(D)
    function alongpath(τ)
        τ <= 0 && return anchor
        c = last(normsearch(ResidueForm(m.poles, τ .* X, (1 - τ) .* S0 .+ τ .* D, wref, m.spaces)))
        isfinite(c) || throw(ArgumentError(
            lazy"the contraction of the fit could not be measured within maxtime = $(maxtime) s. Raise maxtime, fit with fewer poles, or over a narrower band."))
        return c
    end
    # Toward nothing the norm of `t S` is `t` times the norm of `S`, so the
    # step the ceiling gives is tried first; the ceiling measured afresh
    # carries the search's roundoff, and where that leaves it above the
    # target, as with a statement, the step is searched for on the path,
    # from the least a contraction may take, `1 - 2 scalelimit`, which is
    # measured first, so that a fit needing more is refused at once.
    least = 1 - 2*scalelimit
    direct = (1 - 4eps())/ceiling
    s = if anchor == 0 && direct >= least && alongpath(direct) <= target
        direct
    else
        atleast = alongpath(least)
        atleast <= target || throw(ArgumentError(isempty(dc) ?
            lazy"making the fit passive takes a contraction to less than $(least) of the data, which describes a block which is mostly loss rather than the data. Fit with fewer poles, or over a narrower band." :
            lazy"making the fit passive takes a contraction to less than $(least) of the data, which describes a block which is mostly the value it states at zero frequency rather than the data. Fit with fewer poles, over a narrower band, or without stating the value at zero."))
        contractionstep(alongpath, target; lo = least, flo = atleast, fhi = ceiling)
    end
    X .*= s
    D .= (1 - s) .* S0 .+ s .* D
    # what the contraction moved, over the samples, which is warned of
    # where it exceeds the tolerance the block is validated to
    before = ResidueForm(m.poles, Xbefore, Dbefore, wref, m.spaces)
    G = zeros(ComplexF64, n, n)
    moved = maximum(xs) do x
        opnorm(responseat!(F, m, x, phi) .- responseat!(G, before, x, phi))
    end
    moved > atol || return scaledresidues(m) .* wref, D, nothing
    return scaledresidues(m) .* wref, D, (factor = s, moved = moved)
end

# A fit not made passive is refused as any rational block is (see
# checkpassive), where its largest singular value exceeds the level of
# `atol` (see passivelevel) at some frequency: at infinite frequency,
# where it is the feedthrough's, and elsewhere where the sweep of its
# residue form at that level finds it (see passivitysweep). A fit the
# sweep finds under the level by no more than the roundoff, which it
# cannot settle either way, is accepted, as a realization the validation
# cannot settle is; one whose response it evaluates above the level is
# refused however small the excess, so that at an `atol` of zero a
# lossless fit, which stands at the level, is accepted or refused as the
# roundoff of its evaluation falls. The sweep stops at its first reading
# of the clock `maxtime` seconds after it began, the enforcement's
# default, and a fit it has not settled by then is refused, as the
# enforcement refuses one. The poles and residues are in rad/s, at the
# angular frequencies `ws`, each residue taken within its space (see
# residuespaces), as the realization is.
function checkrawpassive(poles::Vector{ComplexF64}, residues::AbstractArray{<:Complex,3}, D::AbstractMatrix,
        ws::AbstractVector, spaces::Vector{Matrix{ComplexF64}}; atol::Real, maxtime::Real = sweepmaxtime)
    margin = passivitymargin(D)
    margin < -atol && throw(ArgumentError(lazy"The rational scattering block is not passive at infinite frequency: the minimum eigenvalue of I - D*D' is $(margin)."))
    lo = findfirst(>(0), ws)
    wref = sqrt(ws[lo]*maximum(ws))
    m = ResidueForm(poles, residues, D, wref, spaces)
    sweep = passivitysweep(m, passivelevel(atol), time_ns() + 1e9*maxtime)
    sweep.verdict === :budget && throw(ArgumentError(
        lazy"the raw fit's passivity could not be settled within $(maxtime) s, the sweep of the frequency axis reaching the interval from $(first(only(sweep.bands))*wref) to $(last(only(sweep.bands))*wref) rad/s. Enforce passivity, whose maxtime sets the sweep's time, or fit with fewer poles or over a narrower band."))
    sweep.verdict === :active || return nothing
    F, phi = zeros(ComplexF64, size(m.D)), zeros(ComplexF64, length(m.poles))
    worst, k = findmax(x -> opnorm(responseat!(F, m, x, phi)), sweep.peaks)
    throw(ArgumentError(lazy"The rational scattering block is not passive: its largest singular value over all frequencies is at least $(worst), at $(sweep.peaks[k]*wref) rad/s, which is above the tolerance $(atol)."))
end

"""
    RationalScattering(block::LinearizedScattering, npoles; frequencies = nothing,
        band = nothing, delays = nothing, tol = 1e-2, noisetol = 5e-3,
        padding = 4, fitting = VectorFitting())

The [`LinearizedScattering`](@ref) block with every harmonic transfer
function fitted to a stable rational realization, which is how the
transient realizes it: `H_0` as an ordinary rational function with its
constant term, or as its constant alone with no state where the data
needs no pole, and each `H_k` for `k > 0` as the pair of real rational
functions of its cosine and sine parts (see
[`ModulatedRationalProvider`](@ref)), strictly proper, since a
conversion vanishes at infinite frequency. Each is fitted at `npoles`
poles by the vector fit of the `ScatteringParameters` method, with the
parameters of `fitting` (see [`VectorFitting`](@ref)), at the
`frequencies` in Hz, by default the magnitudes of the frequencies the
harmonic is tabulated at, within `band = (flo, fhi)` in Hz when given:
a block built from a solve carries every sideband its mode truncation
reached, far above the band a signal occupies, and a harmonic with no
sample in the band is realized as zero; no passivity is enforced, since
a pumped block is lossless as a whole and its parts are not. A
harmonic is read where its data covers a frequency and is zero beyond
its tables, and `H_0`, a real function, is read at whichever sign of a
frequency its data holds and mirrored to the other,
`H_0(-nu) = conj(H_0(nu))`, so a table of one sign fits. A fit
within a band is an approximation within it: a solve evaluates the
block at every sideband its harmonics reach from a signal, where such
a fit only extrapolates, so it serves signals and pulses in the band
and states, through its completed noise, that it is no better outside
it. A lumped device fits over all of its sidebands at a few poles
each; a long line does not, its sidebands being dispersive delay of
many turns of phase. `delays`, one per port in seconds, removes a
delay from the data before the fit, as a cable's is removed before its
fit: the entry from port `q` to port `p` at the harmonic `k` is fitted
with `exp(i (nu + k wp) tau_p + i nu tau_q)` taken out, so the fitted
block is the device with a lossless line of delay `tau_p` cut off each
port, and is put back with a [`TransmissionLine`](@ref) of that delay
in cascade at the port; the delay of a line's signal band is not that
of its sidebands, which keep the difference. The pump phase of
`block` is folded into the fitted functions, and the fitted block keeps
its pump, ports, noise model and envelope; a stated covariance is
rotated by the phase and the delays as the functions are.

The data must meet what the block declares, that it is lossless or the
covariance it states, over the modes its harmonics reach from every
frequency it holds, each output against every input which feeds it
(see [`pumpedfamily`](@ref)), to the block's `atol` and a covariance's,
which the fit checks first and refuses otherwise. The fit itself meets
the declaration no better than its error, so the fitted block does not
declare it: its noise is the covariance the block states, zero for a
lossless one, completed to the commutation relations of the fitted
functions over the ladder of the modes of a solve padded by `padding`
multiples of the pump frequency (see [`NoiseCovariance`](@ref)), so
that it adds, whatever modes a solve keeps, the noise its own
commutator requires, the least a channel with the fitted functions can
add for a lossless device, and its output obeys the commutation
relations exactly. That noise is what the fit costs, and the fit is
refused when it exceeds, in quanta, `noisetol` of the square of the largest entry
over the modes the data reaches from its frequencies and from the
midpoints between them; the block's `atol` stays that of the data. The fit is held to the data as
well, and refused where it misses a sample by more than `tol` of the
largest response, in the spectral norm, as the `ScatteringParameters`
method measures a fit: a stated covariance large enough covers the
commutator of a poor fit at no noise, so the noise a fit adds says
nothing of its accuracy, and the two are held apart. A `band` must leave the
unconverted response a sample. A `dcmodel` the block states is met by
the fit of `H_0` exactly. The harmonic balance solvers evaluate the
fitted block too, so the two describe the same block.
"""
function RationalScattering(block::LinearizedScattering, npoles::Integer; frequencies = nothing,
        band = nothing, delays = nothing, tol::Real = 1e-2, noisetol::Real = 5e-3, padding::Integer = 4,
        fitting::VectorFitting = VectorFitting())
    npoles >= 1 || throw(ArgumentError("fit at least one pole."))
    padding >= 0 || throw(ArgumentError("padding must be nonnegative."))
    (isfinite(tol) && tol >= 0) || throw(ArgumentError("tol must be finite and nonnegative."))
    (isfinite(noisetol) && noisetol >= 0) || throw(ArgumentError("noisetol must be finite and nonnegative."))
    isnothing(band) || (length(band) == 2 && 0 <= band[1] < band[2]) || throw(ArgumentError("band is (flo, fhi) in Hz with 0 <= flo < fhi."))
    taus = something(fitdelays(delays, block.nports), zeros(block.nports))
    n = block.nports
    providers = AbstractMatrixProvider[]
    # the delays taken out: of the output port at the output frequency
    # and of the input port at the input frequency
    undelay!(H, nus, k) = for (i, nu) in enumerate(nus), q in 1:n, pp in 1:n
        H[pp, q, i] *= cis((nu + k*block.wp)*taus[pp] + nu*taus[q])
    end
    # a harmonic at the frequencies its data covers, and zero beyond
    sample(p, nus) = evaluatecovered!(Array{Complex{Float64},3}(undef, n, n, length(nus)), p, nus)
    # the largest deviation of a fitted function from its data over the
    # samples, in the spectral norm, and the largest response, the data
    # as it is fitted, with the pump phase and the delays folded in
    fiterr, datascale = 0.0, 0.0
    function measure!(provider, nus, H)
        worst, scale = fitdeviation(provider, nus, H)
        fiterr, datascale = max(fiterr, worst), max(datascale, scale)
        return nothing
    end
    for (j, k) in enumerate(block.harmonics)
        p = block.providers[j]
        knots = tableknots(p)
        fs = if !isnothing(frequencies)
            Float64.(collect(frequencies))
        elseif !isempty(knots)
            sort!(unique!(filter(>(0), abs.(knots ./ (2pi)))))
        else
            throw(ArgumentError(lazy"give the frequencies in Hz to sample the harmonic $(k) at; only a tabulated harmonic has its own."))
        end
        isnothing(band) || filter!(f -> band[1] <= f <= band[2], fs)
        if isempty(fs)
            # a harmonic with nothing in the band converts nothing there;
            # the unconverted response is the block's reflection and
            # transmission, which the fit cannot leave out
            k == 0 && throw(ArgumentError("the band leaves the unconverted response of the block without a sample; widen it to the frequencies the block is tabulated at."))
            push!(providers, ModulatedRationalProvider(emptyrational(n), emptyrational(n)))
            continue
        end
        checkfrequencies(fs)
        ws = 2pi .* fs
        # the harmonic at the positive and the negative frequencies where
        # its data covers them and zero beyond, the data being the
        # samples and what lies between them, with the block's pump phase
        # folded in
        rot = cis(k*block.phase)
        Hp = sample(p, ws)
        Hp .*= rot
        undelay!(Hp, ws, k)
        Hm = sample(p, -ws)
        Hm .*= rot
        undelay!(Hm, -ws, k)
        if k == 0
            # the unconverted response is a real function, `H_0(-nu) =
            # conj(H_0(nu))`, so a frequency the data holds at one sign
            # alone, an idler's, which a block built from a solve holds
            # at the negative frequency of its mode, or a table of one
            # sign, is fitted from that sign
            for i in eachindex(ws)
                providercovers(p, ws[i]) || (Hp[:, :, i] .= conj.(view(Hm, :, :, i)))
            end
            # a block which states what it does at zero frequency has the
            # fit of its unconverted response meet it exactly; a response
            # which needs no pole is its constant, with no state
            dc = block.dcmodel isa ScatteringLimit ? nodc : dcscatteringmatrix(block.dcmodel, n)
            poles, residues, D = vectorfit(Hp, ws, Int(npoles), fitting; dc = dc, lastpole = :drop)
            A, B, C = realization(poles, residues, n; ranktol = fitting.ranktol)
            push!(providers, RationalScatteringProvider(A, B, C, D))
            measure!(providers[end], ws, Hp)
            continue
        end
        Gc = (Hp .+ conj.(Hm)) ./ 2
        Gs = (Hp .- conj.(Hm)) ./ (2im)
        scale = max(maximum(abs, Hp), maximum(abs, Hm), floatmin(Float64))
        parts = RationalScatteringProvider[]
        for G in (Gc, Gs)
            if maximum(abs, G) <= fitting.zerotol*scale
                push!(parts, emptyrational(n))
                continue
            end
            poles, residues, D = vectorfit(G, ws, Int(npoles), fitting; proper = true, lastpole = :keep)
            A, B, C = realization(poles, residues, n; ranktol = fitting.ranktol)
            push!(parts, RationalScatteringProvider(A, B, C, zeros(n, n)))
        end
        push!(providers, ModulatedRationalProvider(parts[1], parts[2]))
        measure!(providers[end], ws, Hp)
        measure!(providers[end], -ws, Hm)
    end
    # the fit against its data: refused where it misses a sample by more
    # than tol of the largest response
    fiterr <= tol*datascale || throw(ArgumentError(lazy"the fit misses the data by $(fiterr/datascale) of the largest response over the samples, against the tol of $(tol): fit with more poles or over a narrower band, take a delay out, or raise tol to accept a fit that far from the data."))
    # the pump phase and the delays are folded into the fitted functions,
    # so a stated covariance is rotated the same way, `<n(nu + k wp) n(nu)'>`
    # between the ports `p` and `q` by the delays of both emitted waves,
    # each at the frequency its own wave is emitted at; a lossless block
    # states no noise, a covariance of zero
    noise = block.noise
    if noise isa NoiseCovariance
        vps = AbstractMatrixProvider[]
        for (j, k) in enumerate(block.harmonics)
            v = noise.provider[j]
            if block.phase != 0 || any(!iszero, taus)
                v = RotatedMatrixProvider(v, taus; offset = k*block.wp, phase = cis(k*block.phase))
            end
            push!(vps, v)
        end
        noise = NoiseCovariance(vps, noise.interpolation, noise.extrapolation, noise.atol, true, Int(padding))
    else
        noise = NoiseCovariance(AbstractMatrixProvider[ConstantMatrixProvider(zeros(Complex{Float64}, n, n)) for _ in block.harmonics],
            :cubic, :error, 1e-8, true, Int(padding))
    end
    # the fitted block with its covariance completed, over which the
    # noise the completion adds is measured against the covariance as
    # stated
    fitted = LinearizedScattering(block.harmonics, providers, block.wp, 0.0, n, block.zref,
        block.grounded, noise, block.dcmodel, block.envelope, block.atol)
    nus = reduce(vcat, (tableknots(p) for p in block.providers); init = Float64[])
    isempty(nus) && !isnothing(frequencies) && (nus = 2pi .* Float64.(collect(frequencies)))
    nus = sort!(unique!(abs.(nus)))
    filter!(>(0), nus)
    probes = isempty(nus) ? nus : vcat(nus, (nus[1:end - 1] .+ nus[2:end]) ./ 2)
    # the data must meet the declaration itself, at its own frequencies
    # and over the modes it holds there, each output with every input
    # which feeds it, to the tolerance of the block and of a covariance,
    # so that what the fit adds is its own error alone and never a
    # declaration the data violated: every block is checked here,
    # whatever kind its data is, since the declaration is the block's and
    # not its data's
    declared = declaredtolerance(block)
    for nu in nus
        rows, cols, K = pumpedfamily(block, (nu,))
        isempty(rows) && continue
        v = pumpedviolation(block, rows, cols, K)
        v <= declared || throw(ArgumentError(lazy"the block's data does not meet what it declares: over the modes its harmonics reach from $(rows[argmin(abs.(rows))]) rad/s the violation of its losslessness or of the commutation relations of its stated covariance is $(v) of the square of its largest entry, against the $(declared) of its atol and its noise model's; state its noise, or raise atol if that is meant."))
    end
    # the noise the completion adds to the fit, the completed covariance
    # of the fitted functions less the stated one, over the modes the
    # data reaches from its frequencies and from the midpoints between
    # them, since a solve evaluates the fit at frequencies of its own,
    # relative to the square of the largest entry: what the fit costs,
    # refused beyond what is accepted
    added = 0.0
    for nu in probes
        rows, cols, K = pumpedfamily(block, (nu,))
        isempty(rows) && continue
        S, Kc, V = pumpednoisematrices(fitted, rows, cols, K; complete = false)
        extra = pumpednoisematrices(fitted, rows, cols, K)[3] - V
        added = max(added, maximum(real.(eigvals(Hermitian(extra))))/max(1.0, maximum(abs, S))^2)
    end
    added <= noisetol || throw(ArgumentError(lazy"the fit adds noise of $(added) of the square of its largest entry to obey the commutation relations, against the noisetol of $(noisetol) accepted: fit with more poles or over a narrower band, or raise noisetol to accept a block which adds that much."))
    return fitted
end

# the realization of no states and no feedthrough, the fit of a part
# which is zero
emptyrational(n::Int) = RationalScatteringProvider(zeros(0, 0), zeros(0, n), zeros(n, 0), zeros(n, n))
