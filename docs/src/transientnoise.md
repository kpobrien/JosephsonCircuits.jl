# Quantum noise in time

[`transientnoise`](@ref) propagates small fluctuations about a recorded
classical trajectory. The trajectory can include a pump, strong signals,
depletion, and intermodulation. The calculation is a Gaussian
linearization about that mean, not a simulation of the full nonlinear
quantum state.

The workflow is to record the trajectory, define the output temporal
modes, choose physical baths and their frequency quadrature, and inspect
both the noise and its convergence.

## A passive two-port

Start with a circuit whose answer is known: a passive network at zero
temperature must preserve vacuum when all baths are included. The resistor
below joins two capacitively loaded 50 Ω ports.

```@example noise
using JosephsonCircuits, LinearAlgebra
circuit = Circuit([
    (:p1, 1, 0, Port(1)), (:p2, 2, 0, Port(2)),
    (:c1, 1, 0, Capacitor(0.3e-12)),
    (:c2, 2, 0, Capacitor(0.5e-12)),
    (:loss, 1, 2, Resistor(30.0)),
])
problem = transientproblem(circuit)
n, T, w = 512, 1e-9, 2pi*3e9
solution = transientsolve(problem, (0.0, T*(n - 1)/n);
    dt = T/n, record = :phases)
measurement = transientquantumplan(solution, solution.times, [w, w]; ports = [1, 2])

# This circuit is time invariant and the modes are single Fourier bins.
# Only the bath frequency at that bin contributes to each stationary mode.
noise = transientnoise(solution, measurement;
    frequencies = [w], weights = [2pi/T], inputs = measurement)
@assert noise.diagnostics.passed
@assert isapprox(noise.covariance, measurement.vacuum; rtol = 1e-5, atol = 1e-6)
round.(noise.covariance; digits = 4)
```

The four rows are `X1, P1, X2, P2`, and the expected covariance is `I/2`.
The single-frequency bath is appropriate for this stationary,
bin-aligned example. **A pulse or a pumped trajectory generally requires
many bath frequencies**, including frequencies outside the measurement
band. See [choosing the bath grid](#Choosing-the-bath-grid).

`N` uniform samples define a Fourier interval of duration `N*dt`, even
though the last sample is at `(N-1)*dt`. That is why the endpoint above is
slightly less than `T`: the 3 GHz measurement then lies exactly on a bin
of a 1 ns interval.

### Compare gain and thermal noise with harmonic balance

For this time-invariant circuit, convert each complex HB entry into the
real quadrature block `[real(z) imag(z); -imag(z) real(z)]`. The signs
follow the temporal mode's cosine/sine convention.

```@example noise
quadratures(M) = reduce(vcat, [reduce(hcat,
    [[real(z) imag(z); -imag(z) real(z)] for z in M[j, :]])
    for j in axes(M, 1)])
hb = hblinsolve([w], circuit; keyedarrays = false, returnCnoise = true)
@assert isapprox(noise.gain, quadratures(hb.S[:, :, 1]); rtol = 1e-4, atol = 1e-4)

# Warm only internal loss. The port terminations stay at zero kelvin.
baths = transientnoisebaths(problem; temperature = 0.3)
warm = transientnoise(solution, measurement; frequencies = [w],
    weights = [2pi/T], baths)
hbwarm = hblinsolve([w], circuit; keyedarrays = false,
    temperature = 0.3, returnCnoise = true)
S = hbwarm.S[:, :, 1]
expected = quadratures(S*S'/2 + hbwarm.Cnoise[:, :, 1])
@assert isapprox(warm.covariance, expected; rtol = 1e-4, atol = 1e-5)
@assert isapprox(warm.addedcovariance, quadratures(hbwarm.Cnoise[:, :, 1]);
    rtol = 1e-4, atol = 1e-5)
[b.temperature for b in baths.channels]
```

The port temperatures are zero and the resistor is at 0.3 K. The total
covariance holds both the incident vacuum, half a photon at each input
(`S*S'/2`), and the noise the circuit adds, `Cnoise`, which is the
transient's `addedcovariance`: the part from the internal baths, the
resistors and the blocks, without the port terminations'. Both solvers
count the vacuum as one half, so the quadrature form needs no further
factor.

For a nonlinear amplifier with a smooth pump ramp, distinct finite input
and output windows, and independent timestep/bath refinements, continue
with the [pumped-noise example](recipes/pumped-noise.md).

## Temporal modes

[`transientquantumplan`](@ref) defines photon-normalized measurements of
outgoing power waves. A mode is a unit-norm vector of coefficients over
positive Fourier bins. For `N` samples with spacing `dt`, those bins are at
the angular frequencies `2pi*k/(N*dt)` for `k=1:fld(N-1,2)`; DC and the
self-conjugate Nyquist bin are excluded.

For a canonical bin of angular frequency `ω_k` over duration `T=N*dt`, the
wave is

```math
a_k(t)=\sqrt{\frac{\hbar\omega_k}{T}}
\left[X_k\cos(\omega_k(t-t_0))+P_k\sin(\omega_k(t-t_0))\right].
```

Here `ħ` is the reduced Planck constant. The quadratures obey `[X,P]=im`;
vacuum has variance one half in each quadrature. The coherent photon number
is `(X^2+P^2)/2` for a normalized mode.

A frequency-only constructor selects bin-aligned modes. Supplied sampled
envelopes are projected onto the positive-frequency bins and normalized.
The weighting `1/sqrt(hbar*w)` is applied separately at each bin before
combination, which matters for broadband pulses. With coefficients `c`,
the measured annihilator is `A_j=sum(conj(c[k,j])*a_k)`.

Modes on the same port can overlap. `plan.gram` records their overlap;
`plan.vacuum` and `plan.commutator` include the resulting cross-mode terms.
Use those matrices as the reference, rather than assuming `I/2` and
independent canonical pairs. [`transientquantumvjp!`](@ref) is the
transpose of the measurement. [`transientiq`](@ref) instead measures
classical peak phasors.

## Baths and the initial fluctuations

The classical record must start at an equilibrium compatible with a
constant passive prehistory. The check includes node flux and rate, line
waves, and rational-block states. For a pumped nonlinear circuit, a smooth
ramp from an equilibrium is one way to meet this requirement.

Zero classical initial voltage does not mean zero initial quantum noise.
For each bath frequency, the calculation initializes the stationary
response of the circuit about the initial state. This includes energy
stored in reactive elements and its correlation with subsequent forcing.
Lines carry the fluctuations that entered before the record; rational
blocks carry their stationary internal fluctuations.

Trapezoidal and backward Euler use the stationary response of their own
discrete rules. Gauss–Legendre uses the continuous stationary response,
whose discretization mismatch is fourth order at resolved frequencies.
The initial-noise treatment must therefore be included in time-step and
cutoff convergence checks. See [theory](transienttheory.md#Quantum-noise).

### Bath types and temperatures

[`transientnoisebaths`](@ref) constructs external port baths, internal
resistor baths, and supported scattering-block channels. Each external
port must own a finite matched termination. An infinite resistor is open
and adds no bath.

Ports, internal resistors and scattering blocks take the temperatures
and noise models of the
[temperature table](conventions.md#Noise-normalization-and-temperature),
as in harmonic balance; the temperature passed to `transientnoisebaths`
is that of a resistor which states none.

For a resistor `R` at a positive angular frequency `w`, a quadrature weight
`dw`, both in rad/s, gives cosine and sine Norton-current amplitudes

```math
I_{\mathrm{peak}}=2\sqrt{\frac{\hbar\omega\,d\omega}{2\pi R}}.
```

Their independent quadratures have variance `nbar + 1/2`, with
`nbar = thermaloccupation(w, T)`.
The current's variance over the band is then
`(2*hbar*w/R)*coth(hbar*w/(2kT))*dw/(2pi)`, approaching Johnson's
`4kT/R*dw/(2pi)` at high temperature.

## The calculation

The example's result fields are:

| Field | Meaning |
|---|---|
| `covariance` | Total symmetrized output quadrature covariance, harmonic balance's `Vout` |
| `addedcovariance` | Its part from the internal baths, the noise the circuit adds, harmonic balance's `Cnoise` |
| `commutator` | Commutator propagated from the baths |
| `expectedcommutator` | Algebra defined by the measurement modes |
| `gain` | Stationary quadrature response to `inputs`, when requested |
| `diagnostics` | Commutator and uncertainty-relation checks |

The default `method=:adjoint` propagates the measured quadratures backward
and contracts their bath-response kernels with the frequency quadrature.
`method=:forward` propagates two directions per bath and frequency. It is
useful as a small reference calculation but can be much more expensive.
Both methods contract the same discrete responses.

## Choosing the bath grid

With no custom grid, the bath uses all positive Fourier bins of the
record, at the angular frequencies `2pi*k/T` with spacing and weight
`2pi/T`, below Nyquist. A `cutoff` in rad/s limits that grid.
Alternatively, supply positive `frequencies` and corresponding positive
quadrature `weights`, both in rad/s.

For a pulse, begin with enough frequency coverage to include resonances
and pump-converted contributions to the measured modes. Refine the
frequency spacing to resolve the measurement window's spectrum. A list of
just the applied tones and nominal idlers is not generally a complete
bath: the loaded trajectory can mix a continuous band into the window.

A calculation continuing the setup above, using the default record bins,
can be written as:

```julia
full = transientnoise(solution, measurement)
bounded = transientnoise(solution, measurement; cutoff = 2pi*20e9)
```

These calls illustrate grid selection; agreement must be checked for the
particular trajectory. Smooth measurement envelopes usually suppress
remote spectral leakage more effectively than rectangular windows.

When `inputs` is supplied for stationary gain, the bath must cover its
Fourier coefficients on the matching bins with weights `2pi/T`. Do not use
arbitrary integration weights and interpret that gain as the same input
normalization.

## Pulsed gain or stationary gain

| Calculation | Input being measured |
|---|---|
| `transientgain(solution, measurement, inputs)` | Probe confined to its input window, including the response to its edges |
| `transientnoise(solution, measurement; inputs).gain` | Periodic Fourier mode extended over the record and initialized stationary |

Use `transientgain` for a causal pulse response. An output window ending
before the input begins has zero response. Use the stationary definition
for comparison with settled HB gain. They approach one another when the
windows are long compared with the circuit's memory; for short windows,
their difference is expected.

## Diagnostics and convergence

The diagnostics compare propagated and expected commutators and test the
uncertainty relation

```math
V+\frac{i}{2}\Omega_{\mathrm{expected}}\succeq0.
```

A failed check is reported, not repaired by rescaling the result. Passing
is necessary but does not establish convergence. Refine independently:

1. Time step, with the classical trajectory and derivatives recomputed.
2. Bath cutoff, including relevant circuit resonances and conversion bands.
3. Bath-frequency spacing and quadrature weights.
4. Record length, settling time, and measurement window.
5. Fitting band/order and harmonic coverage for fitted pumped blocks.

Compare covariance and the selected gain, not only the commutator.
Near an amplifier threshold, frequency warping and incomplete settling
can dominate a comparison with HB. The [JPA validation recipe](recipes/transient-wrspice.md)
illustrates the role of settling time.

## Scalar gain and efficiency

[`transientquantumefficiency`](@ref) accepts a phase-insensitive one-mode
channel: a scaled rotation for phase-preserving gain, or a scaled
reflection for phase conjugation, together with isotropic covariance
`v*I`. It reports photon gain `G` and `QE=G/(2v)`, with the package's
ideal-gain reference. The covariance holds every bath in its state, the
input's own port included, so the QE is that of the measurement, as
`hbsolve` reports it.

Phase-sensitive gain and anisotropic covariance need the full matrices;
the scalar reduction rejects them. See [`transientquantumefficiency`](@ref)
for its tolerances and the convention below unity gain.

## Scattering blocks and correlated baths

A passive block emits covariance `(nbar + 1/2)(I-S*S')`; an active block can supply
`NoiseCovariance(V)`. Rational-block covariances depend on frequency and
are contracted directly, retaining small losses without subtracting large
output covariances. Their internal states also enter the initial-noise
response.

A pumped block correlates bath frequencies separated by pump harmonics,
including conjugate partners. The bath grid must include those partners.
The covariance and commutator are assembled on padded pump ladders,
consistent with the HB model. A fitted block's covariance is completed
against its fitted transfer functions. See [pumped-block contracts](scattering.md#Correlated-noise-and-fit-tolerances).

Correlated pumped-block frequencies cannot be tiled independently, so the
required frequency accumulator must fit the memory budget. For ordinary
independent baths, condition and frequency tiling trade extra adjoint
passes for memory. See [noise cost](performance.md#Noise-calculation-cost)
and [implementation notes](implementation.md#Pumped-block-noise).
