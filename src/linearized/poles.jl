# Temporal Floquet exponents of a periodic operating point. The nonlinear
# solve supplies only the orbit; the polynomial below is the physical
# variational equation, without its Newton gauge pins or port-wave scales.

# The fraction of a mode's states below which its node voltages and
# junction fluxes are the roundoff an eigensolve leaves in coordinates
# the mode does not reach, which makes it internal to its blocks, or to
# its lines under the period map, whose waves count with the states.
const INTERNALFRACTION = sqrt(eps(Float64))

"""
    HBStabilityResult

The finite poles [`hbstability`](@ref) found, with time dependence
`exp(s*t)`, in decreasing order of their real parts, every array in that
order:

- `poles`: the growth rates `s` in inverse seconds, growing where the
  real part is positive.
- `residuals`: each pole's relative backward error in the scaled
  polynomial, which meets the method's `tol`, or under
  [`Monodromy`](@ref) with a pump a condition-based estimate of multiplier
  roundoff relative to the multiplier. It excludes timestep error,
  which `rateerrors` estimates separately.
- `rateerrors`: for pumped [`Monodromy`](@ref), an estimate of each real
  part's timestep error. The selected modes are propagated through a
  period in twice the steps, their line histories resampled and read
  back in the coarse coordinates, and the refined map projected onto
  their right and left vectors, `G = (VL' VR) \\ VL' M2 VR`, which follows
  their mixing among themselves. Each mode continues into the
  eigenvector of `G` that carries most of it, of multiplier `μ2`; the
  estimate is `(16/15*abs(log(abs(μ2/μ))) + residuals[k])/T`, the coarse
  map's error where it falls as the fourth power of the step. `Inf`
  means no refined mode continues the mode. Mixing with modes beyond
  `nev` is not followed, and a map too coarse for its error to fall so
  can have more error than its estimate. `abs(real(poles[k])) >
  rateerrors[k]` is evidence of a resolved sign, growth where
  `real(poles[k]) > rateerrors[k]` and decay where
  `real(poles[k]) < -rateerrors[k]`, not a certificate: confirm it by
  doubling `steps`, which cuts the estimate about sixteenfold once the
  steps resolve the mode, and where modes crowd by a larger `nev`.
  `nothing` for other methods and without a pump.
- `edgeweights`: the fraction of each mode's squared amplitude, node
  voltage over `frequencyscale` and junction flux together, in the two
  outermost harmonics, zero without a pump. A large one says to retain
  more harmonics in a polynomial method, and a wider profile under the
  period map; a small one bounds no truncation error.
- `nodevoltage`, `junctionflux`: each mode's coefficients, `(mode, node,
  pole)` in volts and `(mode, junction, pole)` in webers, normalized
  together, along `modes`, `nodes` (ground excluded) and
  `junctionbranches` (the junctions' branches in the compiled circuit,
  whose orientation sets the sign). For a nonzero profile,
  `sum(abs2, nodevoltage[:,:,k]/frequencyscale) + sum(abs2, junctionflux[:,:,k])`
  is one. The mode's amplitude and global phase are arbitrary; these
  coefficients do not predict a driven oscillation amplitude. Harmonic
  `n` has angular frequency `imag(s) + n*pumpfrequency`.
- `stateedgeweights`: the same fraction for the states of rational
  blocks, which depends on their realization, and under
  [`Monodromy`](@ref) the waves leaving the lines' ports with them, in
  the voltage they carry over `frequencyscale`.
- `internalonly`: modes whose node/junction content is negligible beside
  their balanced rational-state content or, under pumped Monodromy, the
  waves leaving their lines' ports; their exported profiles are zero and
  their dominant harmonic is selected from those states and waves. Zero
  profiles do not imply that a mode is absent or stable.
- `harmonics`: each mode's dominant harmonic, the `n` whose coefficients
  carry most of its squared amplitude, by the norm of `edgeweights` (by
  its states' and waves' for a mode internal to blocks and lines): the
  mode lives mostly at
  the angular frequency `imag(s) + n*pumpfrequency`, which tells a pole's
  Floquet aliases apart by their content, and one whose dominant harmonic
  is an outermost retained one is shaped by the truncation.
- `shiftindices`: the shift of [`ShiftInvert`](@ref) each pole was found
  near, zero for the other methods.
- `searches`: diagnostics for each search. A shift records requested and
  converged counts and products taken. A contour records quadrature,
  rank and capacity, saturation, changes in poles and moments, and the
  nearest pole's distance to the circle. Pumped Monodromy records steps,
  period, the multipliers not classified as unresolved (`multipliers`),
  how many of the largest were tested against their own condition-based
  bounds (`tested`), gauge directions removed, and line-history columns.
  Multipliers within `eps()*opnorm(M,1)` of zero after balancing are
  classified as unresolved without an individual test and do not count
  in `tested`. With `nev = :all` every multiplier above that threshold is
  tested against its own bound; with a smaller `nev` the untested ones
  are taken against the threshold alone, as well conditioned, so the
  counts can depend on `nev`. The array is empty for DenseSpectrum and
  the unpumped Monodromy fallback.
- `converged`: whether the finite eigensolve/search met its acceptance
  checks. Polynomial/contour candidates must meet `tol`; `rejected`
  counts failed candidates. Pumped Monodromy classifies unresolved
  multipliers separately and can return `converged=true` with `Inf`
  rate errors. This field does not establish discretization accuracy,
  search completeness for the physical circuit, or stability.
- `infinite`: the infinite eigenvalues [`DenseSpectrum`](@ref) set aside,
  those [`ShiftInvert`](@ref)'s searches came upon (see there), or the
  multipliers of the period map within the eigensolve's error of zero,
  eps times the norm of the balanced map, and those of the largest within
  their own, that over their condition: its algebraic directions and the
  modes damped beyond a period's precision. The map tests its largest
  multipliers against their own bounds until `nev` pass, and the others
  against the bound of a well conditioned one alone, so that the count
  is the eigensolve's for every multiplier with `nev = :all` and can be
  smaller with fewer; `nothing` for [`ContourIntegral`](@ref).
- `pumpfrequency`, `frequencyscale`: angular frequencies in rad/s.
- `method`: the method which found the poles.

A physical mode has Floquet aliases `s + im*k*pumpfrequency`. Harmonic
truncation can perturb these copies differently; polynomial methods
return the candidates they find. A pole found near two shifts is returned
twice. The period map returns one representative per selected mode, up to `nev`.
`converged` concerns the eigensolve:
neither it nor a dense spectrum establishes that the circuit has no
unstable mode the truncation in harmonics or in steps misses.
"""
struct HBStabilityResult
    poles::Vector{ComplexF64}
    residuals::Vector{Float64}
    rateerrors::Union{Nothing,Vector{Float64}}
    edgeweights::Vector{Float64}
    nodevoltage::Array{ComplexF64,3}
    junctionflux::Array{ComplexF64,3}
    modes::Vector{Tuple{Int}}
    nodes::Vector{String}
    junctionbranches::Vector{Int}
    shiftindices::Vector{Int}
    searches::Vector{NamedTuple}
    converged::Bool
    infinite::Union{Nothing,Int}
    rejected::Int
    pumpfrequency::Float64
    frequencyscale::Float64
    stateedgeweights::Vector{Float64}
    internalonly::BitVector
    harmonics::Vector{Int}
    method::Any
end

"""
    ShiftInvert(shifts; nev = 6, krylovdim = max(30, 3nev), restarts = 200,
        tol = 1e-9)

The method of [`hbstability`](@ref) which finds the `nev` poles nearest
each of the complex `shifts`, in inverse seconds; `imag(shift) = 2pi*f`
searches near the frequency `f` in Hz. Each shift takes one
factorization of the polynomial there, and ArnoldiMethod's restarted
Schur iteration on the companion pencil shifted and inverted through it,
which never forms the pencil, with at most `krylovdim` basis vectors and
`restarts` restarts, to the tolerance `tol`, which bounds the residual a
pole is accepted with as well. A shift exactly at a pole makes the
polynomial singular there; move it off by a little.

A search for more poles than lie near its shift comes upon the pencil's
infinite eigenvalues, which the inversion maps to zero. Of its
candidates which fail the polynomial, those which outnumber the finite
poles it can still have are `infinite`: the degree of the polynomial's
determinant, which the sums of its rows' and of its columns' degrees
bound, less the poles the search accepted. The others are `rejected`, a
remote finite pole the search did not resolve among them, and the
search has not converged.
"""
struct ShiftInvert
    shifts::Vector{ComplexF64}
    nev::Int
    krylovdim::Int
    restarts::Int
    tol::Float64
end
function ShiftInvert(shifts; nev::Integer = 6, krylovdim::Integer = max(30, 3nev),
        restarts::Integer = 200, tol::Real = 1e-9)
    s = shifts isa Number ? ComplexF64[shifts] : ComplexF64.(collect(shifts))
    (!isempty(s) && all(isfinite, s)) || throw(ArgumentError("shifts must be nonempty and finite."))
    (1 <= nev < krylovdim && restarts >= 0) || throw(ArgumentError(
        "require 1 <= nev < krylovdim, a basis with room to grow, and restarts >= 0."))
    (isfinite(tol) && tol > 0) || throw(ArgumentError("tol must be finite and positive."))
    return ShiftInvert(s, Int(nev), Int(krylovdim), Int(restarts), Float64(tol))
end

"""
    DenseSpectrum(; tol = 1e-9, maxunknowns = 1000)

The method of [`hbstability`](@ref) which finds every finite pole of the
truncated problem, its Floquet aliases included, by the QZ algorithm on
its companion pencil formed densely, twice the unknowns square: the
reference for a small circuit. Its time grows as the cube of the
unknowns, every node's harmonics counted, so a problem of more than
`maxunknowns` is refused before its pencil is formed; [`Monodromy`](@ref)
finds the selected least damped modes of a pumped circuit, and
[`ShiftInvert`](@ref) the poles near given frequencies. A pole is
accepted with a residual under `tol`.
"""
struct DenseSpectrum
    tol::Float64
    maxunknowns::Int
end
function DenseSpectrum(; tol::Real = 1e-9, maxunknowns::Integer = 1000)
    (isfinite(tol) && tol > 0) || throw(ArgumentError("tol must be finite and positive."))
    maxunknowns >= 1 || throw(ArgumentError("maxunknowns must be positive."))
    return DenseSpectrum(Float64(tol), Int(maxunknowns))
end

"""
    ContourIntegral(center, radius; quadrature = 64, moments = 4,
        probes = 12, ranktol = 1e-10, changetol = 1e-6, refinements = 3,
        tol = 1e-9)

The method of [`hbstability`](@ref) which finds the poles inside the
circle of `center` and `radius`, in inverse seconds, by the block moment
integral of Beyn (Linear Algebra Appl. 436, 2012, 3839), the method for
[`LaplaceResponse`](@ref) models and exact [`TransmissionLine`](@ref)
delays, whose operator is not a polynomial; [`Monodromy`](@ref) takes
the lines of a pumped circuit as well. It factors the operator at
`quadrature` points of the circle, refactoring its fixed pattern, and
doubles the points, keeping those it has, up to `refinements` times,
until the poles move by less than `changetol` of the radius and the
moments by less than `changetol` of the mean norm of the solves sampled
on the circle. `probes` right hand sides and `moments` moments bound how
many poles it holds, `ranktol` cuts the moments' rank relative to the
larger of their largest singular value and that mean norm, so that an
empty circle converges, and a pole is accepted with a residual under
`tol`. A saturated rank or a failed check sets
`converged = false`; the checks are numerical, not a certified count of
the poles. Keep the poles away from the circle, and the branch cuts of an
analytic response outside it.
"""
struct ContourIntegral
    center::ComplexF64
    radius::Float64
    quadrature::Int
    moments::Int
    probes::Int
    ranktol::Float64
    changetol::Float64
    refinements::Int
    tol::Float64
end
function ContourIntegral(center::Number, radius::Real; quadrature::Integer = 64,
        moments::Integer = 4, probes::Integer = 12, ranktol::Real = 1e-10,
        changetol::Real = 1e-6, refinements::Integer = 3, tol::Real = 1e-9)
    (isfinite(center) && isfinite(radius) && radius > 0) || throw(ArgumentError(
        "the contour's center must be finite and its radius finite and positive, in inverse seconds."))
    (quadrature >= 8 && moments >= 1 && probes >= 1 && refinements >= 1) || throw(ArgumentError(
        "require quadrature >= 8, moments >= 1, probes >= 1 and refinements >= 1."))
    all(t -> isfinite(t) && 0 < t < 1, (ranktol, changetol)) || throw(ArgumentError(
        "ranktol and changetol must lie strictly between zero and one."))
    (isfinite(tol) && tol > 0) || throw(ArgumentError("tol must be finite and positive."))
    return ContourIntegral(ComplexF64(center), Float64(radius), Int(quadrature), Int(moments),
        Int(probes), Float64(ranktol), Float64(changetol), Int(refinements), Float64(tol))
end

"""
    hbstability(circuit, circuitdefs = Dict(); nonlinear = nothing,
        Nmodulationharmonics = (8,), method = Monodromy(),
        pumpfrequency = nothing, frequencyscale = nothing,
        factorization = KLUfactorization())

The local stability of the periodic operating point `nonlinear`, a
converged [`hbnlsolve`](@ref) solution with one pump frequency, or of
zero junction phase without one: the temporal poles of the circuit
linearized about it, as an [`HBStabilityResult`](@ref), whose real parts
tell whether a small perturbation of the orbit grows. They say nothing of
what larger perturbations do, and settle the orbit's stability only as
far as the truncation, in harmonics or steps, is converged. The circuit
and its definitions must be those of the operating point; commensurate
drives are harmonics of its one pump.

The perturbation is expanded in every harmonic `-H:H` of
`Nmodulationharmonics = (H,)`, whatever the pump's symmetry, over the
junction modulation [`hblinsolve`](@ref) samples from the same operating
point; without an operating point or `pumpfrequency` it has the one
harmonic `0`. `pumpfrequency` gives the fundamental of pumped blocks
without an operating point, and must agree with one. The polynomial
methods truncate the perturbation at these harmonics: converge their
poles in the pump's harmonics and in these, which `edgeweights` and a
solve at a larger `H` tell. The period map, [`Monodromy`](@ref), takes
them as the width of its profiles alone, requiring `2H + 1 <= steps`.
Converge its poles in the pump's harmonics and in its `steps`. Its
`rateerrors` is an estimate from the map at twice the steps, not a
bound; near zero growth, compare full solves at increasing `steps`.

Choose a method for the desired spectral coverage:

- [`Monodromy`](@ref)`()` returns up to 10 least damped resolved modes of
  a pumped circuit, one representative each, including line delays.
  Set `nev` to change that count. Without a pump it uses
  `DenseSpectrum()` with the default 1000-unknown limit.
- [`DenseSpectrum`](@ref) returns accepted finite poles of a small
  harmonic truncation, including aliases.
- [`ShiftInvert`](@ref) searches near specified complex frequencies.
- [`ContourIntegral`](@ref) searches inside a circle. It is required for
  [`LaplaceResponse`](@ref) models and exact delays without a pump.

A shift or contour search can miss instabilities outside its region.
`factorization` factors the polynomial on the host,
[`KLUfactorization`](@ref) or [`LUfactorization`](@ref), refactoring
its fixed pattern, and the transient's steps for the period map on the
host. `frequencyscale`, an angular frequency, scales the polynomial: the
pump frequency by default, or a natural frequency of the unpumped
circuit; it changes no pole.

The variational equation is `Q(s) = K + G (s + D) + C (s + D)^2`, with
`D = diag(im*n*pumpfrequency)` and the junction stiffness along the orbit
in `K`. Floating flux references are replaced by voltages to retain
physical zero poles and RC decay. Ports keep their terminations, and
independent sources are held fixed.

Supported models include real constant lumped elements, mutual inductors,
junctions, nonlinear inductors, constant and [`RationalScattering`](@ref)
blocks with internal states, and [`LinearizedScattering`](@ref) blocks.
Use a polynomial method for constant/rational conversion harmonics, or
[`ContourIntegral`](@ref) for analytic Laplace harmonics. Pumped
[`Monodromy`](@ref) requires a supported time-domain realization.
Tabulated data or a callable of real frequency needs a causal rational
fit or [`LaplaceResponse`](@ref); complex constant lumped values are
refused. Refine a positive pole in harmonics or steps before interpreting
it as an instability of the physical circuit.

# Examples
```jldoctest
julia> c = Circuit([(:r, 1, 0, Resistor(2.0)),
                   (:c, 1, 0, Capacitor(0.5)),
                   (:l, 1, 0, Inductor(2.0))]);

julia> p = hbstability(c);

julia> sort(imag.(p.poles)) ≈ [-sqrt(3)/2, sqrt(3)/2]
true

julia> all(isapprox.(real.(p.poles), -0.5)) && p.converged
true
```
"""
function hbstability(circuit::CompilableCircuit,
        circuitdefs::AbstractDict = Dict{Symbol,Any}(); nonlinear = nothing,
        Nmodulationharmonics = (8,), method = Monodromy(),
        pumpfrequency = nothing, frequencyscale = nothing,
        factorization = KLUfactorization())
    method isa Union{Monodromy,ShiftInvert,DenseSpectrum,ContourIntegral} || throw(ArgumentError(
        "method is Monodromy(), DenseSpectrum(), ShiftInvert(shifts) or ContourIntegral(center, radius)."))
    factorization isa Union{KLUfactorization,LUfactorization} || throw(ArgumentError(
        "hbstability factors on the host: give KLUfactorization() or LUfactorization()."))
    psc = compile(circuit)
    if method isa Monodromy
        Nmodulationharmonics isa Tuple{Int} && Nmodulationharmonics[1] >= 0 ||
            throw(ArgumentError("Nmodulationharmonics must be a tuple (H,) with H >= 0."))
        sys = hbpolesystem(psc, circuitdefs, nonlinear, (0,), frequencyscale; pumpfrequency)
        # the transient carries a line's history, and no Laplace model
        laplacemodels(psc, circuitdefs) && throw(ArgumentError(
            "LaplaceResponse models need method = ContourIntegral(center, radius)."))
        # Without a pump the map of a period is the exponential of the
        # polynomial, whose roots are found directly, and a delay's are
        # not. Each path is reached dynamically, so that a call compiles
        # the one it takes and not the other with it.
        if iszero(sys.pumpfrequency)
            isempty(sys.terms) || throw(ArgumentError(
                "without a pump, exact delays need method = ContourIntegral(center, radius)."))
            return Base.invokelatest(densepoles, sys, DenseSpectrum(), method)::HBStabilityResult
        end
        return Base.invokelatest(polemonodromy, method, psc, circuitdefs, nonlinear, sys,
            only(Nmodulationharmonics), factorization)::HBStabilityResult
    end
    sys = hbpolesystem(psc, circuitdefs, nonlinear, Nmodulationharmonics,
        frequencyscale; pumpfrequency)
    method isa ContourIntegral || isempty(sys.terms) || throw(ArgumentError(
        "exact delays and LaplaceResponse models need method = ContourIntegral(center, radius), and the lines of a pumped circuit take Monodromy() too."))
    return polesolve(method, sys, factorization)
end

# whether a circuit holds a LaplaceResponse model, as a lumped value or a
# block's response, which the transient does not realize; the pole
# system has refused any other callable already
function laplacemodels(psc::CompiledCircuit, circuitdefs)
    _, dynamic = polelumpedvalues(psc, resolvedvalues(psc, circuitdefs))
    isempty(dynamic) || return true
    return any(psc.scatteringblocks) do cb
        d = cb.definition
        any(p -> p isa CallableMatrixProvider, d isa LinearizedScattering ? d.providers : (d.provider,))
    end
end

# every finite pole, from the dense companion pencil
polesolve(method::DenseSpectrum, sys, factorization) = densepoles(sys, method, method)

# the dense spectrum's poles, reported as found by `reported`
function densepoles(sys, method::DenseSpectrum, reported)
    n = size(sys.Q0, 1)
    n <= method.maxunknowns || throw(ArgumentError(reported isa Monodromy ?
        lazy"without a pump the poles are those of the dense spectrum, and its $(n) unknowns exceed the $(method.maxunknowns) it takes by default, its time growing as their cube; give method = DenseSpectrum(maxunknowns = $(n)) for every pole, or ShiftInvert(shifts) for the poles near given frequencies." :
        lazy"the dense spectrum of $(n) unknowns exceeds maxunknowns = $(method.maxunknowns), its time growing as their cube; find selected least damped modes of a pumped circuit with Monodromy(), or the poles near given frequencies with ShiftInvert(shifts)."))
    A, B = polecompanion(sys)
    alpha, beta, _, V = LAPACK.ggev!('N', 'V', A, B)
    finite = .!iszero.(beta)
    # a 0/0 generalized eigenvalue is an undetermined circuit, not an
    # infinite algebraic mode to set aside
    any(iszero.(alpha) .& iszero.(beta)) && throw(ArgumentError(
        "the pole pencil is singular; the circuit has undetermined variables."))
    values = alpha[finite] ./ beta[finite]
    return poleresult(sys, values, V[1:n, finite], zeros(Int, length(values)),
        NamedTuple[], true, count(!, finite), method.tol, reported)
end

# the poles nearest each shift, by the shifted and inverted companion
# pencil, each shift's polynomial refactored on the pattern of the first
function polesolve(method::ShiftInvert, sys, factorization)
    n = size(sys.Q0, 1)
    work = PoleMatrixWorkspace(sys)
    cache = FactorizationCache()
    values, vectors = ComplexF64[], Matrix{ComplexF64}(undef, n, 0)
    shiftindices, searches = Int[], NamedTuple[]
    converged = true
    infinite = 0
    for (k, shift) in enumerate(method.shifts)
        z = shift/sys.scale
        tryfactorize!(cache, factorization, polematrix!(work, sys, z))
        op = PoleShiftInvert(cache.factorization, sys.Q1, sys.Q2, z, zeros(ComplexF64, n), zeros(ComplexF64, n))
        count = min(method.nev, 2n)
        maxdim = min(method.krylovdim, 2n)
        # a restart keeps half the basis, the requested poles at least
        decomp, history = ArnoldiMethod.partialschur(op; nev = count, which = :LM,
            tol = method.tol, mindim = clamp(maxdim ÷ 2, count, maxdim), maxdim,
            restarts = method.restarts)
        mu, V = ArnoldiMethod.partialeigen(decomp)
        # the companion's infinite eigenvalues map to zero: one there
        # exactly is infinite, and the others are candidates, checked
        # against the polynomial (see poleresult)
        finite = findall(x -> isfinite(x) && !iszero(x), mu)
        infinite += sum(x -> isfinite(x) && iszero(x), mu; init = 0)
        append!(values, z .+ inv.(mu[finite]))
        vectors = hcat(vectors, V[1:n, finite])
        append!(shiftindices, fill(k, length(finite)))
        push!(searches, (requested = count, converged = history.nconverged,
            products = history.mvproducts))
        converged &= history.converged
    end
    return poleresult(sys, values, vectors, shiftindices, searches, converged,
        infinite, method.tol, method; finite = finitebound(sys))
end

# An upper bound on the number of finite poles, the degree of the
# polynomial's determinant: the smaller of the sums of its rows' and of
# its columns' degrees, each the highest power whose coefficient has an
# entry there, which leaves out the infinite eigenvalues of every row or
# column without the polynomial's leading terms.
function finitebound(sys)
    n = size(sys.Q0, 1)
    rows, cols = zeros(Int, n), zeros(Int, n)
    for (degree, Q) in ((1, sys.Q1), (2, sys.Q2)), j in 1:n, p in nzrange(Q, j)
        iszero(nonzeros(Q)[p]) && continue
        i = rowvals(Q)[p]
        rows[i] = max(rows[i], degree)
        cols[j] = max(cols[j], degree)
    end
    return min(sum(rows), sum(cols))
end

# Scaled polynomial and maps back to physical perturbations. A common
# coordinate is voltage / scale; the other node coordinates are flux
# relative to that component's reference. Auxiliary currents are Lscale*i.
struct HBPoleSystem
    Q0::SparseMatrixCSC{ComplexF64,Int}
    Q1::SparseMatrixCSC{ComplexF64,Int}
    Q2::SparseMatrixCSC{ComplexF64,Int}
    fluxmap::SparseMatrixCSC{Float64,Int}
    commonmap::SparseMatrixCSC{Float64,Int}
    junctionmap::SparseMatrixCSC{Float64,Int}
    offsets::Vector{Float64}
    modes::Vector{Tuple{Int}}
    nodes::Vector{String}
    junctionbranches::Vector{Int}
    scale::Float64
    pumpfrequency::Float64
    terms::Vector{PoleResponseTerm}
    stateindices::Vector{Int}
end

function hbpolesystem(psc::CompiledCircuit, circuitdefs, nonlinear,
        Nmodulationharmonics, frequencyscale; pumpfrequency = nothing)
    Nmodulationharmonics isa Tuple{Int} && Nmodulationharmonics[1] >= 0 ||
        throw(ArgumentError("Nmodulationharmonics must be a tuple (H,) with H >= 0."))
    values = resolvedvalues(psc, circuitdefs)
    vvn, dynamic = polelumpedvalues(psc, values)
    checkstaticstiffnessvalues(psc.componenttypes, vvn)
    checkisolatedsubnetworks(psc)
    omega = isnothing(pumpfrequency) ? 0.0 : Float64(pumpfrequency)
    (isfinite(omega) && (isnothing(pumpfrequency) || omega > 0)) ||
        throw(ArgumentError("pumpfrequency must be finite and positive."))
    if !isnothing(nonlinear)
        nonlinear isa NonlinearHB || throw(ArgumentError("nonlinear must be a NonlinearHB solution."))
        nonlinear.solverinfo.converged || throw(ArgumentError("the nonlinear operating point has not converged."))
        length(nonlinear.w) == 1 || throw(ArgumentError(
            "hbstability requires a periodic solution with one fundamental frequency; express commensurate drives as harmonics."))
        isnothing(pumpfrequency) || samefrequency(omega, only(nonlinear.w), only(nonlinear.w)) ||
            throw(ArgumentError("pumpfrequency differs from the nonlinear operating point."))
        omega = Float64(only(nonlinear.w))
        isfinite(omega) && omega > 0 || throw(ArgumentError("the pump frequency must be finite and positive."))
    end
    freq = calcfreqsdft(iszero(omega) ? (0,) : Nmodulationharmonics)
    m = length(freq.modes)
    nm = numericmatrices(psc, vvn; Nmodes = m)
    if !isnothing(nonlinear)
        originalnm = isempty(dynamic) ? nm : numericmatrices(psc, values; Nmodes = m)
        (nonlinear.nodes == psc.nodenames && nonlinear.Ljb == originalnm.Ljb &&
            nonlinear.Lb == originalnm.Lb) || throw(ArgumentError(
            "the circuit nodes or inductances differ from the operating point."))
    end
    modulation = linearizedmodulation(psc, nm, freq, nonlinear)
    offsets = omega .* only.(freq.modes)
    # For a static circuit a characteristic frequency from both damping
    # and stiffness avoids baking an arbitrary GHz unit into small examples.
    Cnorm, Gnorm = norm(nm.Cnm, Inf), norm(nm.Gnm, Inf)
    Knorm = norm(nm.invLnm, Inf) + sum(abs, inv.(nonzeros(nm.Ljb)))
    natural = Cnorm > 0 ? max(Gnorm/Cnorm, sqrt(Knorm/Cnorm)) :
        Gnorm > 0 ? Knorm/Gnorm : 1.0
    natural = max(natural, maximum(b -> polefrequency(b.definition), psc.scatteringblocks; init = 0.0))
    scale = isnothing(frequencyscale) ? (omega > 0 ? omega :
        (isfinite(natural) && natural > 0 ? natural : 1.0)) : Float64(frequencyscale)
    isfinite(scale) && scale > 0 || throw(ArgumentError("frequencyscale must be finite and positive."))
    coupled = mnacoupledbranches(nm.Mb)
    nnodal = (psc.Nnodes - 1)*m
    blocks, n = poleblocks(psc, omega, m, nnodal + length(coupled)*m)
    naux = n - nnodal
    n > 0 || throw(ArgumentError("the pole system has no unknowns."))
    # Scale currents into flux units before any row equilibration. This
    # matters particularly for coupled inductors with algebraic rows.
    Lscale = abs(real(calcsolverscale([scale], psc.componenttypes,
        nm.vvn, nm.portimpedances, nm.Lmean)))
    isfinite(Lscale) && Lscale > 0 || (Lscale = 1.0)
    C = mnapad(ComplexF64.(nm.Cnm)*Lscale, naux)
    G = mnapad(ComplexF64.(nm.Gnm)*Lscale, naux)
    K = mnapad(ComplexF64.(nm.invLnm)*Lscale, naux)
    M = calcAmnaind(coupled, nm.Lb, nm.Mb, psc.topology.Rbn,
        m, nnodal, n, Lscale)
    R = hcat(nm.Rbnm, spzeros(eltype(nm.Rbnm), size(nm.Rbnm, 1), naux))
    # The shared constructor stamps the junction derivative. Its Lscale
    # argument is one, so multiply the modulation after construction.
    lsys = HBLinearizedSystem(modulation.Amatrixindices, nm.Ljb, R,
        m, psc.topology.Nbranches, modulation.phimatrix, K, G, C,
        K, G, C, false, M, modulation.wpumpmodes, psc.Nnodes)
    K = copy(lsys.Asparse)
    K.nzval .*= Lscale
    K += lsys.invLnm + M
    stampoleblocks!(K, G, blocks, freq.modes, m, Lscale, scale)

    # K*Tcommon is identically zero: a uniform flux shift inside any
    # component of the L/Lj graph changes no constitutive branch flux.
    # Dividing those columns by (lambda + i*n*omega/scale) replaces the
    # reference flux with its voltage, without a root cancellation test.
    floating = calcstaticfluxcomponents(psc.componenttypes, psc.nodeindices,
        vvn, psc.Nnodes)
    rows, cols, vals = collect(1:n), collect(1:n), ones(Float64, n)
    common = falses(n)
    for component in floating
        ref = first(component) - 1
        for mode in 1:m
            col = (ref-1)*m + mode
            common[col] = true
            for node in component[2:end]
                push!(rows, (node-2)*m + mode)
                push!(cols, col)
                push!(vals, 1.0)
            end
        end
    end
    T = sparse(rows, cols, vals, n, n)
    keep = Float64.(.!common)
    fluxmap = T[1:nnodal, :]*spdiagm(0 => keep)
    commonmap = T[1:nnodal, :]*spdiagm(0 => Float64.(common))
    Ct = C*T*scale^2
    Gt = G*T*scale
    Kt = K*T
    # Enforce the analytic zeros of common flux columns instead of leaving
    # cancellation roundoff from a large junction network in the pencil.
    Kt = Kt*spdiagm(0 => keep)
    d = im .* repeat(offsets ./ scale, n÷m)
    Q0 = Kt + Gt*spdiagm(0 => ifelse.(common, 1.0+0.0im, d)) +
        Ct*spdiagm(0 => ifelse.(common, d, d.^2))
    Q1 = Gt*spdiagm(0 => keep) + Ct*spdiagm(0 => ifelse.(common, 1.0+0.0im, 2 .* d))
    Q2 = Ct*spdiagm(0 => keep)
    fluxrows, commonrows = sparse(transpose(fluxmap)), sparse(transpose(commonmap))
    terms = poleblockterms(blocks, freq.modes, offsets, fluxrows, commonrows, scale, Lscale)
    append!(terms, polelumpedterms(psc, values, dynamic, m, offsets, fluxrows, commonrows, scale, Lscale))
    # Row equilibration does not change the roots or the physical vectors.
    weights = vec(sum(abs, Q0; dims = 2) + sum(abs, Q1; dims = 2) + sum(abs, Q2; dims = 2))
    for term in terms
        localweights = vec(sum(abs, poleterm(term, 1.0+0.0im, scale); dims = 2))
        for (i, row) in enumerate(term.rows)
            weights[row] += localweights[i]
        end
    end
    all(x -> isfinite(x) && x > 0, weights) || throw(ArgumentError(
        "the pole system has a zero or nonfinite equation; check the circuit constraints."))
    E = spdiagm(0 => inv.(weights))
    Q0, Q1, Q2 = E*Q0, E*Q1, E*Q2
    dropzeros!(Q0); dropzeros!(Q1); dropzeros!(Q2)
    for term in terms
        term.left ./= weights[term.rows]
    end
    stateindices = Int[]
    for block in blocks, filter in block.filters
        append!(stateindices, filter.statebase .+ (1:polenstates(filter.provider)*m))
    end
    jb = nm.Ljb.nzind
    jrows = [(b-1)*m + mode for b in jb for mode in 1:m]
    junctionmap = real.(nm.Rbnm[jrows, :])*fluxmap
    return HBPoleSystem(Q0, Q1, Q2, fluxmap, commonmap, junctionmap,
        offsets, freq.modes, String.(psc.nodenames[2:end]),
        copy(jb), scale, omega, terms, stateindices)
end

function polecompanion(sys::HBPoleSystem)
    n = size(sys.Q0, 1)
    A, B = zeros(ComplexF64, 2n, 2n), zeros(ComplexF64, 2n, 2n)
    for i in 1:n
        A[i, n+i] = 1
        B[i, i] = 1
    end
    A[n+1:end, 1:n] = -sys.Q0
    A[n+1:end, n+1:end] = -sys.Q1
    B[n+1:end, n+1:end] = sys.Q2
    return A, B
end

struct PoleShiftInvert{F}
    factor::F
    Q1::SparseMatrixCSC{ComplexF64,Int}
    Q2::SparseMatrixCSC{ComplexF64,Int}
    shift::ComplexF64
    work::Vector{ComplexF64}
    rhs::Vector{ComplexF64}
end

Base.size(op::PoleShiftInvert) = (2length(op.work), 2length(op.work))
Base.size(op::PoleShiftInvert, i::Integer) = size(op)[i]
Base.eltype(::PoleShiftInvert) = ComplexF64

function LinearAlgebra.mul!(y::AbstractVector, op::PoleShiftInvert, x::AbstractVector)
    n = length(op.work)
    u, v = view(x, 1:n), view(x, n+1:2n)
    a, b = view(y, 1:n), view(y, n+1:2n)
    # (A-shift*B)^-1 B [u;v], eliminating b = u + shift*a.
    @. op.work = v + op.shift*u
    mul!(op.rhs, op.Q2, op.work)
    mul!(op.rhs, op.Q1, u, 1, 1)
    ldiv!(op.factor, op.rhs)
    @. a = -op.rhs
    @. b = u + op.shift*a
    return y
end

# the order of poles by decreasing real part, then increasing imaginary
# part, as a result lists them
poleorder(poles) = sortperm([(-real(z), imag(z)) for z in poles])

function poleresult(sys, values, vectors, shiftindices, searches,
        converged, infinite, tol, method; finite = nothing)
    m, nn, nj = length(sys.modes), length(sys.nodes), length(sys.junctionbranches)
    count = length(values)
    voltage = zeros(ComplexF64, m, nn, count)
    flux = zeros(ComplexF64, m, nj, count)
    errors, edges, stateedges = zeros(count), zeros(count), zeros(count)
    internal = falses(count)
    dominant = zeros(Int, count)
    norms = (opnorm(sys.Q0, Inf), opnorm(sys.Q1, Inf), opnorm(sys.Q2, Inf))
    H = maximum(abs(only(mode)) for mode in sys.modes)
    edge = findall(mode -> abs(only(mode)) == H, sys.modes)
    for k in eachindex(values)
        z = values[k]
        q = view(vectors, :, k)
        nq = norm(q)
        if !(isfinite(z) && isfinite(nq) && nq > 0)
            errors[k] = Inf
            continue
        end
        q ./= nq
        r = sys.Q0*q + z*(sys.Q1*q) + z^2*(sys.Q2*q)
        # Allow absolute roundoff in Q0 at a physical zero pole, even if
        # Q0 is identically zero (a capacitor with no discharge path).
        # Keep the degree-one and degree-two weights coefficient-specific:
        # a vanishing leading coefficient must not make enormous spurious
        # finite approximations of algebraic infinity pass this check.
        qmax = norm(q, Inf)
        denom = (maximum(norms) + abs(z)*norms[2] + abs2(z)*norms[3])*qmax
        for term in sys.terms
            T = poleterm(term, z, sys.scale)
            for (j, col) in enumerate(term.cols), (i, row) in enumerate(term.rows)
                r[row] += T[i,j]*q[col]
            end
            denom += opnorm(T, Inf)*qmax
        end
        errors[k] = norm(r, Inf)/max(denom, floatmin(Float64))
        isfinite(errors[k]) && errors[k] <= tol || continue
        f = reshape(sys.fluxmap*q, m, nn)
        v = reshape(sys.commonmap*q, m, nn)
        v .+= (z .+ im .* sys.offsets ./ sys.scale) .* f
        j = reshape(sys.junctionmap*q, m, nj)
        # Normalize the physical mode, with voltage/scale and flux in the
        # same units; retain their relative amplitudes and phases.
        amplitude = sqrt(sum(abs2, v) + sum(abs2, j))
        states = isempty(sys.stateindices) ? zeros(ComplexF64, m, 0) : reshape(q[sys.stateindices], m, :)
        internal[k] = amplitude <= INTERNALFRACTION*norm(states)
        if !isempty(sys.stateindices) && H > 0
            total = sum(abs2, states)
            total > 0 && (stateedges[k] = sum(abs2, states[edge, :])/total)
        end
        # the harmonic which carries most of the mode, by the same norm, or
        # by its states' for a mode internal to its blocks, whose profile is
        # zero
        content = internal[k] ? vec(sum(abs2, states; dims = 2)) : vec(sum(abs2, v; dims = 2)) .+ vec(sum(abs2, j; dims = 2))
        dominant[k] = only(sys.modes[argmax(content)])
        internal[k] && continue
        v ./= amplitude
        j ./= amplitude
        voltage[:, :, k] = sys.scale .* v
        flux[:, :, k] = j
        if H > 0
            edges[k] = sum(abs2, view(v, edge, :)) + sum(abs2, view(j, edge, :))
        end
    end
    accepted = findall(e -> isfinite(e) && e <= tol, errors)
    # A search's eigenvalues are the polynomial's poles, finite or not:
    # of the candidates which fail the polynomial, as many as outnumber
    # the finite poles it can have besides those the search accepted,
    # `finite` at most, are infinite, and the others rejected.
    atinfinity = 0
    if !isnothing(finite)
        for k in unique(shiftindices)
            members = sum(==(k), shiftindices; init = 0)
            found = sum(i -> shiftindices[i] == k, accepted; init = 0)
            atinfinity += max(0, members - found - max(0, finite - found))
        end
    end
    rejected = count - length(accepted) - atinfinity
    isnothing(infinite) || (infinite += atinfinity)
    converged &= iszero(rejected)
    order = accepted[poleorder(view(values, accepted))]
    return HBStabilityResult(sys.scale .* values[order], errors[order], nothing, edges[order],
        voltage[:, :, order], flux[:, :, order], sys.modes, sys.nodes,
        sys.junctionbranches, shiftindices[order], searches, converged,
        infinite, rejected, sys.pumpfrequency, sys.scale, stateedges[order], internal[order],
        dominant[order], method)
end
