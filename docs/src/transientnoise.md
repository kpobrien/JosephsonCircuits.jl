# Quantum noise in time

The noise calculation linearizes about the complete recorded trajectory,
pump, signals, depletion and intermodulation products together, and
propagates the fluctuations of the physical baths, the port terminations
and the resistors, through that time dependent linear system: the
Gaussian approximation of harmonic balance's noise calculation carried
into time, small fluctuations about a classical mean rather than the
full quantum state. It runs on the tangent and the adjoint of the
recorded steps, so it is exact for the discrete trajectory it is given
and needs no interpolation, replay or quadrature grid of its own.

## Temporal modes

[`transientquantumplan`](@ref) defines photon normalized measurements of
the outgoing power waves in sqrt(W). Where [`transientiq`](@ref) returns
classical peak phasors, this plan returns real quadratures in the row
order `X1, P1, X2, P2, ...`, with `[X, P] = im`, vacuum covariance `I/2`,
and coherent photon number `(X^2 + P^2)/2` for a normalized mode.

A record of `N` uniform samples covers the half open interval
`[t0, t0 + N*dt)`; its positive frequency bins are `f = k/(N*dt)` for
`k = 1:fld(N - 1, 2)`, without DC and without the self conjugate Nyquist
bin. With `T = N*dt` the physical wave of one canonical bin is

```math
w_k(t) = \sqrt{\frac{h f_k}{T}}\,[X_k\cos(2\pi f_k(t - t_0)) + P_k\sin(2\pi f_k(t - t_0))].
```

A mode is a column of unit norm coefficients over those bins, and the
measured annihilator is `A_j = sum(conj(c[k, j])*a_k)`, with the
`1/sqrt(h*f)` weighting applied per bin before the bins are combined,
which matters for a pulse spread over gigahertz. The frequency form of
the constructor makes bin aligned monochromatic modes, or projects
sampled envelopes onto the positive bins and normalizes them. Modes on
the same port need not be orthogonal: `plan.gram` is their overlap and
`plan.vacuum` and `plan.commutator` the corresponding cross mode
quadrature matrices. [`transientquantumvjp!`](@ref) is the exact
transpose of the measurement.

## Baths and the initial fluctuations

[`transientnoisebaths`](@ref)`(problem)` builds one independent
equilibrium bath per matched port termination and per finite internal
resistor, from the compiled port ownership, the bound values and the
component temperatures; a typed [`Circuit`](@ref) expresses which
resistor a port owns. For positive frequency `f` and quadrature weight
`df` in Hz, a bath of resistance `R` injects cosine and sine Norton
currents of peak amplitude

```math
I_{\mathrm{peak}} = 2\sqrt{h f\,df/R},
```

whose independent quadratures have variance
`thermaloccupation(2pi*f, T)/2 = nbar + 1/2`, so the bilateral
symmetrized current spectral density is `h*f/R*coth(h*f/(2kT))`, the
classical `2kT/R` at high temperature. A port's bath is injected as a
port source is, into the port, so its current enters the port's wave as
a source's does; an internal resistor's bath is injected between the
resistor's terminals. As targets of [`transienttangent`](@ref) and
[`transientadjoint`](@ref) the baths are what the noise drives.

The record must start at a classical equilibrium with a constant passive
prehistory, and the check reads the whole of the state the solve started
from: the node fluxes and rates, the waves on the lines, which carry a
direct current, and the states of the rational blocks, which must be at
rest under the incident waves there. For each bath and frequency the
stationary response
`[-w^2 C + i w G + L + J'(x_0)] x = b` of the circuit at the initial
state starts the cosine and sine responses, so the noise already stored
in the capacitors and inductors, and its correlation with the forcing to
come, are kept; zero classical initial voltage does not mean zero
quantum noise.

## The calculation

```julia
using JosephsonCircuits

# a solve recording the junction phases, which is all the noise reads
problem = transientproblem(circuit; sources = [TransientSource(1, pump)])
n, T = 20000, 100e-9
solution = transientsolve(problem, (0.0, T*(n - 1)/n); dt = T/n, record = :phases)
measurement = transientquantumplan(solution.times, [5e9]; ports = [2])
inputs = transientquantumplan(solution.times, [5e9]; ports = [1])
df = 1/T
frequencies = collect(df:df:25e9)   # a grid to refine, not a prescription
noise = transientnoise(solution, measurement; frequencies, weights = fill(df, length(frequencies)), inputs)

noise.diagnostics
noise.covariance          # the symmetrized output covariance
noise.commutator          # the propagated bath commutator
noise.expectedcommutator  # the measurement's own mode algebra
noise.gain                # the incremental quadrature response to the inputs
```

`method = :adjoint`, the default, propagates the measured quadratures
backward through [`transientadjoint`](@ref): its derivatives with
respect to the bath currents at every recorded time are the bath
kernels, and each column of them is contracted against the cosine and
sine of every frequency at the grid or stage time the adjoint hands it
over at, a rank one update
of a (bath, objective) by (frequency, quadrature) accumulator on the
backend, so nothing is stored per recorded time and the frequency count
costs sums rather than integrations; the columns are buffered over a
block of times and contracted by one product. That accumulator holds
`16 nb m nf N` bytes for `nb` baths, `m` measured quadratures, `nf`
frequencies and `N` conditions, which for distributed loss and a dense
band is the memory of the method: within a quarter of the backend's free
memory it is one, beyond it the conditions are tiled, each tile one
adjoint over its conditions, and then the frequencies, each tile one
more adjoint over the same conditions, a trade of passes for memory
rather than a bound; the frequency count itself is set by where the
measured modes respond, a few bands around the tones and their idlers
at the resolution of the window, not by the whole record. The
covariance, the commutator and the gain accumulate from each tile on
the backend, and the stationary operator is factorized once per
frequency for all conditions sharing an initial state. The stationary initial term enters
through the adjoint's initial flux and rate with one stationary solve
per frequency and objective, so no bath's initial state is ever built;
the memory of the adjoint method is the record of the solve itself plus
that accumulator. `method = :forward` propagates the bath
quadratures forward, two tangent directions per bath and frequency
through [`transienttangent`](@ref), all directions of every condition
of a batch on each step's factorization; it is the reference the adjoint
is checked against, and the two agree to roundoff since they contract
the same responses.

Two gains are on offer. [`transientgain`](@ref)`(solution, measurement,
inputs)` drives each input mode's unit quadratures as incident waves
inside the input window only, evaluated at the grid and the stage times
of the steps inside it and zero elsewhere, through the tangent, and
reads them in the measurement's modes: the causal, pulsed gain, with
the transients of the probe's own edges, independent of any bath grid,
and nothing of it reaches a window that ends before the probe starts.
The `gain` of
[`transientnoise`](@ref) with an input plan is the response to the
periodic Fourier mode of the window extended over the record and
started stationary, which is what a stationary amplifier's harmonic
balance gain is; it needs the bath on the window's Fourier bins with
weights `1/T`. The two agree as the window grows past the circuit's
memory.

Every dissipative element is a bath: the port terminations, the
internal resistors, and the lossy scattering blocks, whose emitted
noise wave has the covariance `I - S S'` of Bosma's relation, the same
the linearized solver uses, factored into independent channels that
enter the block's port current rows as the source of the hybrid
equation; a lossless block, a circulator or a through, emits nothing.
A [`RationalScattering`](@ref) block's covariance depends on
frequency, so its ports are channels correlated by `I - S S'` at each
bath frequency, contracted directly with that covariance so that a
block lossless to a part in a million is never the difference of two
large terms, and its states enter the stationary response and its
initial term; a lossy rational block leaves the vacuum the vacuum when
cold and has the excess `hblinsolve` gives it when warm. Nothing is
inferred about a rational block's loss: it keeps its channels however
small its loss, and a declared `Lossless()` is validated by the norms
of the block and of its inverse over all frequencies, which can only
refuse it.
A block's temperature is its `ThermalEquilibrium(T)`, a resistor's its
`temperature`, and the rest the analysis default, in both solvers
alike, so a warm attenuator in time has the covariance
[`hblinsolve`](@ref) gives it, and equals its own resistive network at
the same temperature.

A transmission line stores the fluctuations that entered it before the
record began. The stationary response of the circuit to each bath tone
is solved on the flux phasors together with the phasors of the waves
leaving the line ports, nothing inverted in the lines, so it is regular
at a line's half wave resonances; the forward method hands the tangent
those waves over the prehistory the delays reach into, and the adjoint's
initial term contracts its cotangents of that prehistory through the
transposed system. A cold line between mismatched loads then leaves the
vacuum the vacuum, and a pumped amplifier behind a cable has the gain
and the quantum efficiency the linearized solver gives it with the same
line.

The bath quadrature is any set of positive frequency nodes with positive
weights in Hz; a list of driven tones and idlers is not a complete bath
for a pulse, since the loaded trajectory mixes the whole band into a
measured window. The optional input plan returns the incremental
quadrature gain from the same responses, and requires the measurement's
record, bin aligned frequencies with weights `1/T`, and coverage of every
input coefficient.

What convergence looks like on a pulsed, distributed loss case, a two
hundred cell line of series junctions with a loss resistor in every
cell, pumped by a pulse and measured in Hann windows of 2 ns on the
rise, the plateau and the fall of the pulse with the output window
delayed by the line's transit: the step converges at fourth order, the
covariance and the pulsed gain moving sixteen times less from 2.5 ps to
1.25 ps than from 5 ps to 2.5 ps; the bath spacing converges as the
window's spectrum is resolved, at a part in ten thousand by a quarter of
the window's bin; and the cutoff must pass the junction plasma frequency
of the line, near 29 GHz there, where the line responds most and the
pump mixes into the measured band, a band stopping short of it missing
a part in a thousand and one past it changing nothing further. A
rectangular window converges slowly in the cutoff, since its sidelobes
carry noise from far away; a smooth envelope does not. The test suite
runs that case.

The diagnostics check the commutator against the measurement's algebra
and the uncertainty relation `covariance + im*expectedcommutator/2 >= 0`,
and a failure is returned as such rather than rescaled away. Passing is
necessary, not sufficient: refine the bath cutoff, the frequency
spacing, the record and the step independently. On a passive two port
the quadrature gain agrees with [`hblinsolve`](@ref) to a part in ten
thousand at 128 samples of a 3 GHz record and to three parts in a hundred
thousand at 512, converging at second order with the step.

[`transientquantumefficiency`](@ref) reduces a phase preserving single
mode result, a scaled rotation for the gain and an isotropic covariance,
to photon gain `G`, added noise `v/G - 1/2`, `QE = G/(2v)` and the
package's ideal `G/(2G - 1)`; phase sensitive gain or anisotropic noise
is rejected and the full matrices remain. On a pumped Josephson
amplifier measured after 200 ns of settling with the stationary Floquet
frequencies as the bath, the gain and the quantum efficiency converge to
[`hbsolve`](@ref) at the order of the stepping rule: under the
trapezoidal rule, whose frequency warping a resonator's detuning and an
amplifier's bifurcation magnify, the gain is 12% high at 5 ps and half
a percent at 0.16 ps; under [`GaussLegendre`](@ref) it agrees to 1.5e-3
at 5 ps and to 2e-5 at 2.5 ps, with the quantum efficiency to 6e-6, in
about a second of solve and noise together. On a traveling wave amplifier of two hundred series junctions pumped to
2.6 radians of node flux, where the accumulated propagation phase and
the phase matching are the constraint, the periodic gain and the
quantum efficiency agree with [`hbsolve`](@ref) to 4e-4 at 5 ps and
3e-5 at 2.5 ps, again at fourth order. Harmonic balance takes
milliseconds on these stationary problems and stays their reference,
the transient being for pulses and drives too heavily loaded for a
harmonic grid.
