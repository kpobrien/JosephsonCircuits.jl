# Units and conventions

## Frequencies

| API or quantity | Unit |
|---|---|
| HB `ws`, `wp`, and `w` | rad/s |
| `FrequencyDependent(w -> ...)` | rad/s, evaluated at the magnitude of the mode frequency |
| Scattering callables, tabulated scattering frequencies, and `RationalScattering` fitting `frequencies` and `band` | rad/s |
| `transientdemodulate`, `transientiqplan`, `transientquantumplan`, and the plans' `frequencies`, `bandwidth3db` and `noisebandwidth` | rad/s |
| `transientnoise` `frequencies`, `cutoff` and quadrature `weights` | rad/s |
| Time, line delays, transient step size | s |

Use `w = 2pi*f` to convert Hz to rad/s. For ordinary component values,
negative-frequency evaluation follows conjugate symmetry. A scattering
provider can instead supply its own signed-frequency data; see
[`ScatteringParameters`](@ref).

## Current amplitudes

HB sources use complex Fourier coefficients. A single nonzero coefficient
`Ip` at angular frequency `wp` represents

```math
I(t) = 2\operatorname{Re}\{I_p e^{i\omega_p t}\}.
```

Its peak current is `2abs(Ip)` and its RMS current is `sqrt(2)*abs(Ip)`.
The zero-frequency coefficient is the DC current itself, without doubling.
`TransientSource` takes the instantaneous real current. Thus the same
real pump is written as `current = Icoeff` in HB and as
`t -> 2Icoeff*cos(wp*t)` in a transient solve.

A port source is a Norton current injected into the port's first
(positive) terminal. With a matched termination of resistance `R`, a
sinusoid of peak current `Ipeak` launches available power `Ipeak^2*R/8`.
A [`CurrentSource`](@ref) drives its current through itself from its
first terminal to its second: it draws the current from the node at its
first terminal and delivers it to the node at its second, the opposite
sense of a port source. A port and a `CurrentSource` written on the same
nodes in the same order therefore drive in opposite senses; see the
[direct current example](recipes/dc.md).

### Check the convention on a resistor

The circuit below is just a 50 Ω port termination. Its voltage Fourier
coefficient is `50Icoeff`; its instantaneous voltage has twice that peak.
The transient initial voltage is supplied because a resistor cannot start
at zero voltage under a nonzero current.

```@example amplitude
using JosephsonCircuits
c = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0))])
f, Icoeff = 1e9, 1e-6
hb = hbnlsolve((2pi*f,), (1,), [(mode = (1,), port = 1, current = Icoeff)],
    c; keyedarrays = false)
Vcoeff = im*2pi*f*JosephsonCircuits.phi0*only(hb.nodeflux)
@assert isapprox(Vcoeff, 50Icoeff; rtol = 1e-10)
drive(I, f) = t -> 2I*cospi(2f*t)
p = transientproblem(c; sources = [TransientSource(1, drive(Icoeff, f))])
initial = transientstate(p; voltage = [2*50Icoeff])
td = transientsolve(p, (0.0, 1/f); dt = 1/(100f), initialstate = initial)
@assert isapprox(td.voltage[1, 1], 2real(Vcoeff); rtol = 1e-10)
(real(Vcoeff), td.voltage[1, 1])
```

## Modes and signed frequencies

A pump mode is an integer tuple `m`, with physical angular frequency
`sum(m .* wp)`. Independent pumps form a multidimensional Fourier grid.
The real-transform representation stores nonnegative indices along the
first grid axis and both signs along the others, with redundant conjugate
modes removed. This is an index-space convention: a stored mode can have
negative physical frequency when there is more than one pump.

A signal mode is an offset from the signal frequency:

```math
\omega_m = \omega_s + \boldsymbol m\cdot\boldsymbol\omega_p.
```

For a single pump and a positive signal below the pump:

| Mode | Signed frequency | Interpretation |
|---|---|---|
| `(0,)` | `ws` | Signal |
| `(-2,)` | `ws - 2wp[1]` | Conjugate idler coordinate in four-wave mixing |
| `(-1,)` | `ws - wp[1]` | Conjugate idler coordinate in three-wave mixing |

A negative idler coordinate represents the conjugate of the physical
positive-frequency tone. The sign matters for phases and commutators.
For commensurate pumps, specify one fundamental frequency and drive its
harmonics through different source mode indices. Independent incommensurate
pumps describe a quasiperiodic state rather than one finite-period orbit.

## Photon gain and power gain

`LinearizedHB.S` is normalized to photon flux. For nonzero frequencies,

```math
G_{\mathrm{photon}} = |S_{oi}|^2,\qquad
G_{\mathrm{power}} = \frac{|\omega_o|}{|\omega_i|}|S_{oi}|^2.
```

Multiply `S` by `sqrt(abs(wout)/abs(win))` for the corresponding power-wave
amplitude. The absolute values are required for signed idler coordinates.
For equal input and output frequencies, the two gains coincide.

Transient `incident` and `outgoing` arrays are instantaneous real power
waves in `sqrt(W)`. A demodulated peak wave amplitude `a` at one tone has
average power `abs2(a)/2`. Do not compare these arrays directly with
photon-normalized HB coefficients.

`NonlinearHB.S` contains output-to-drive ratios for the solved strong-drive
state. With several driven port modes, each column contains the response
to all sources divided by that column's incident wave. Use `LinearizedHB.S`
for the differential scattering response about the operating point.

## Flux, phase, and state arrays

Write physical node flux as `Φ`, in webers, and reduced flux as
`φ = Φ/φ₀`, where `φ₀ = ħ/(2e)` is the reduced flux quantum.

| Quantity | Meaning |
|---|---|
| HB `nodeflux` | Fourier coefficients of reduced node flux |
| HB voltage coefficient at nonzero `w` | `im*w*phi0*nodeflux` |
| HB `dcnodevoltage` | Average node voltage in volts; separate from the static flux |
| `transientstate` input `flux`, `voltage` | Physical node flux in Wb and node voltage in V |
| Transient `voltage` | Port voltages in V |
| Transient `flux`, `rate`, `finalflux`, `finalrate` | Solver-scaled state arrays, including auxiliary coordinates |

Use `transientstate(solution)` to continue a trajectory. It also preserves
block states and line history; copying only the final node arrays loses
that information. For state layout details, see [`TransientState`](@ref).

## Noise normalization and temperature

Noise is in quanta, and the vacuum is half a photon. Every covariance is
symmetrized, `⟨{Δr, Δr†}⟩/2`. A vacuum mode has the covariance `1/2`, and
a thermal mode `nbar + 1/2`, where `nbar` is its occupation, the mean
number of photons, zero in the vacuum.

The noise outputs of harmonic balance, `QE`, `QEideal`, `nbar`, `Vout`,
`Cnoise` and `Snoise`, are defined in
[what the solvers report](portnoise.md#What-the-solvers-report).

Temporal-mode quadratures are ordered `X1, P1, X2, P2, ...`, with
`[X, P] = im` and vacuum covariance `I/2` for orthonormal modes. This is
the same convention: a harmonic balance covariance converted to these
quadratures needs no factor.
`transientnoise(...).covariance` is the **total** measured output covariance,
HB's `Vout` in these quadratures, and `transientnoise(...).addedcovariance`
the noise the circuit **adds**, HB's `Cnoise`. See the
[worked comparison](transientnoise.md#A-passive-two-port).

The quantum-optics functions of quadrature covariances take `hbar`, `1` by
default: [`JosephsonCircuits.is_cptp`](@ref) and `is_cptp_quadrature_pair`
and `_block`, `rand_cptp_quadrature_pair` and `_block`,
`B_from_X_Y_quadrature` and its `_pair` and `_block` forms,
`X_Y_to_symplectic_pair` and `_block`, and `Ymin_from_X` and its forms.
`hbar` fixes the vacuum covariance of the quadratures at `(hbar/2) I`, so
`hbar = 1` is this convention and `hbar = 2` is the shot-noise units of
quantum optics. The ladder functions and the conversions between the bases
take none: a ladder covariance counts the vacuum as `I/2`, and
`ladder_to_quadrature_pair` and `_block` take it to the quadrature
covariance at `hbar = 1`.

These functions are internal to the package, written qualified with its
name and documented in the [appendix](api/internals.md). A function ending
in `_pair` orders the operators of `n` modes in pairs, each mode's two
together, `(x₁, p₁, x₂, p₂, ...)` or `(a₁, a₁†, ...)`; one ending in
`_block` orders them in blocks, `(x₁, ..., xₙ, p₁, ..., pₙ)`.
`pair_to_block`, `block_to_pair` and the `R_` permutations convert
between the two. The symplectic, Bogoliubov and covariance matrices of
`n` modes are `2n × 2n`. The predicates (`is_unitary`,
`is_symplectic_pair`, ...) test a square matrix, dense or sparse, to the
tolerances of `isapprox`; `is_positive_semi_definite` and the `is_cptp`
family take a dense one. The decompositions (`williamson_pair`,
`bloch_messiah_block`, `autonne_takagi`, `polar`, the Iwasawa
decompositions, `symplectic_normal_form_pair`, ...) take a square dense
matrix, a view of one, or a `Symmetric` or `Diagonal` matrix;
`williamson_pair` and `williamson_block` also take a `SparseMatrixCSC`,
whose Cholesky factorization stays sparse. Give the others a sparse
matrix as `Matrix(M)`. The Williamson and Bloch–Messiah decompositions
take real matrices, and `autonne_takagi` a complex or a real symmetric
one.

Each source of noise takes its temperature as follows, in both solvers:

| Contribution | Temperature or noise model |
|---|---|
| External port inputs | The temperature of the port's termination, `Port(n; termination = MatchedTermination(temperature = T))`, zero unless it states one; the analysis `temperature` does not warm it |
| Internal resistor or other supported dissipative element, a resistor across a port included | Its stated temperature; otherwise the analysis default |
| `Passive()` block | Analysis default |
| `ThermalEquilibrium(T)` block | Specified `T`, in kelvin |
| `NoiseCovariance(V)` block | Supplied covariance, independent of the analysis temperature |
| `Lossless()` block | No emitted noise; the declaration is checked where the model permits |

### Temperatures and noise temperatures

A temperature describes a bath; a noise temperature is a noise spectral
density in kelvin. One thermal state has three different "temperatures".
At 5 GHz and a physical 50 mK:

| Quantity | Formula | Value |
|---|---|---|
| Physical temperature | the `T` with `nbar = 1/(exp(ħω/kT) - 1)` | 50 mK, `nbar = 0.0083` |
| Symmetrized noise temperature | `ħω (nbar + 1/2)/k` | 122 mK; the vacuum alone gives `ħω/2k = 120` mK |
| Rayleigh-Jeans temperature | `ħω nbar/k` | 2 mK; zero in the vacuum |

Physical temperatures describe baths, terminations, qubits and
resonators. Symmetrized noise temperatures describe amplifiers and
measurement chains. The Rayleigh-Jeans temperature applies only to data
that excludes the vacuum. The solvers report noise in quanta, and these
functions convert it:

- [`thermaloccupation`](@ref)`(ω, T)` is `nbar` at the physical
  temperature `T`, and [`effectivetemperature`](@ref)`(ω, nbar)` is the
  inverse: the temperature of the bath with that occupation.
- [`noisetemperature`](@ref)`(ω, N)` is `ħ|ω| N/k` for noise `N` in quanta,
  and [`noisequanta`](@ref)`(ω, T)` is the inverse. Pass `nbar + 1/2` for
  the symmetrized noise temperature and `nbar` for the Rayleigh-Jeans
  temperature. `1/(2QE)` gives a chain's system noise temperature,
  referred to the input the QE is taken from, with that input's vacuum
  included.

An input-referred temperature refers to the frequency of the input mode,
which differs from the output's for an idler.

## Derivatives and numerical accuracy

Component sensitivities use a relative perturbation `p -> r*p` at `r=1`:
`dS/dr = p*dS/dp`. `designsensitivities` returns derivatives with respect
to the named design parameters themselves, including dependencies of
multiple component values on the same parameter.

An exact tangent or adjoint differentiates the discrete equations used by
the solver, to its numerical tolerances. It does not remove truncation,
time-step, fitting, or bath-quadrature error. Check those independently
when interpreting a derivative or a noise result.
