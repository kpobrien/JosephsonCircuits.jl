# Scattering blocks and fitted data

A scattering block represents a microwave network through its port waves.
Use it for measured data, a simulated linear network, or the linearized
response of a pumped device. The model keeps its own reference impedances;
the surrounding circuit determines the loading at its terminals.

## A constant two-port

This 3 dB matched attenuator connects two 50 Ω analysis ports. With the
default `grounded=true`, a block's node list contains one signal node per
port; its reference terminals are ground.

```@example pad
using JosephsonCircuits
attenuation = 10^(-3/20)
pad = ScatteringParameters([0.0 attenuation; attenuation 0.0]; zref = 50.0)
circuit = Circuit([
    (:p1, 1, 0, Port(1)),
    (:pad, 1, 2, pad),
    (:p2, 2, 0, Port(2)),
])
response = hblinsolve([2pi*5e9], circuit; keyedarrays = false)
@assert isapprox(abs2(response.S[2, 1, 1]), 10^(-3/10); rtol = 1e-12)
round(10log10(abs2(response.S[2, 1, 1])); digits = 3)
```

A block with `grounded=false` exposes both terminals of each port. In a
connection group, `(:block, p, t)` addresses terminal `t` of port `p`.
This supports floating or differential ports. See [`ScatteringParameters`](@ref)
for the input forms and terminal rules.

## Data providers and frequency coverage

| Provider | Interpretation |
|---|---|
| Matrix | Constant scattering matrix |
| Callable | Matrix evaluated at angular frequency in rad/s |
| `(frequencies, S)` | Tabulated angular frequencies and a port-by-port-by-frequency array |
| Touchstone path | File data and its native reference impedances |

Tabulated data uses cubic interpolation and, by default, rejects
frequencies outside the supplied band. HB can request many frequencies
beyond the signal band: pump harmonics, idlers, and intermodulation
products. Check coverage before treating measured data as a circuit element.

A rational fit supplies a stable time-domain realization and can be
evaluated outside the measured band. Its extrapolation is a model
assumption. Passivity and a small in-band fit error do not establish that
out-of-band behavior matches the physical device.

## Fit a rational model

This complete example fits synthetic data from a passive low-pass two-port.
Tabulated input frequencies are angular frequencies; the optional
`frequencies` keyword used to sample a callable for fitting is in Hz.

```@example fitting
using JosephsonCircuits, LinearAlgebra
fc = 5e9
sample_f = collect(range(0.0, 20e9; length = 81))
lowpass(w) = (h = inv(1 + im*w/(2pi*fc)); [0 h; h 0])
Sdata = cat((lowpass(2pi*f) for f in sample_f)...; dims = 3)
data = ScatteringParameters((2pi .* sample_f, Sdata); nports = 2, zref = 50.0)
fitted = RationalScattering(data, 4; tol = 1e-4)

network(block) = Circuit([
    (:p1, 1, 0, Port(1)), (:device, 1, 2, block), (:p2, 2, 0, Port(2)),
])
# Check held-out points against the analytic response, not only the samples.
check_f = collect(range(0.25e9, 19.75e9; length = 40))
fit_response = hblinsolve(2pi .* check_f, network(fitted); keyedarrays = false)
error = maximum(abs(fit_response.S[2, 1, k] - lowpass(2pi*f)[2, 1])
    for (k, f) in enumerate(check_f))
@assert error < 1e-4
round(error; sigdigits = 3)
```

Increase the requested order if the fit misses the data. The fitter can
remove redundant poles, but more poles do not compensate for insufficient
frequency coverage or unresolved data. `tol` bounds the accepted fit error;
`atol` controls physical-contract checks. See [`VectorFitting`](@ref) and
[`PassivityEnforcement`](@ref) for the algorithm options.

For a nearly lossless device, inspect `I-S*S'` as well as `S`. Noise depends
on the small difference from unitarity, so a modest scattering error can
be a large relative error in loss.

### Use the fitted block in time

Continuing the fit above:

```@example fitting
rise(t) = t <= 0 ? 0.0 : t >= 0.2e-9 ? 1.0 : sinpi(t/0.4e-9)^2
problem = transientproblem(network(fitted);
    sources = [TransientSource(1, t -> 1e-6*rise(t)*cospi(2*3e9*t))])
solution = transientsolve(problem, (0.0, 1e-9); dt = 2e-12,
    method = GaussLegendre())
@assert all(isfinite, solution.voltage)
nothing # hide
```

Constant real blocks and rational realizations use `GaussLegendre()`.
General complex matrices and arbitrary frequency callables do not define
instantaneous real time-domain models.

### Supply a state-space realization

A rational block has `S(s) = D + C*(s*I-A)^(-1)*B`. For example, a series
inductor between 50 Ω reference ports has this one-state realization:

```@example realization
using JosephsonCircuits, LinearAlgebra
L, Z = 1e-9, 50.0
a = 2Z/L
inductor = RationalScattering(fill(-a, 1, 1), [1.0 -1.0],
    -a .* [1.0; -1.0;;], Matrix(1.0I, 2, 2); zref = Z)
nothing # hide
```

The constructor checks stability and the declared noise/passivity contract.
The passivity calculation searches the largest singular value over
frequency using a level-set method and numerical tolerances; it is not
just a check at the supplied samples. The [theory](transienttheory.md#Fitting-measured-data)
describes this calculation.

## Transmission lines and delay removal

[`TransmissionLine`](@ref) takes a characteristic impedance and length.
Its default phase velocity is the speed of light; use `vp` for another
velocity. A lossless cable with impedance 60 Ω and a 300 ps delay at the
default velocity can be written approximately as:

```julia
line = TransmissionLine(60.0, 0.09)
```

In time, the delay must be at least one integration step. Use at least
four steps per delay to retain the Gauss–Legendre rule's fourth-order
accuracy with the history interpolation.

A long propagation delay is expensive to approximate with poles. The
fitting keyword `delays` removes one reference-plane delay per port before
fitting. Restore those delays as explicit lines in the circuit, using the
appropriate reference impedance and phase velocity. Supplied noise
covariances are shifted to the same reference planes. Delay removal does
not infer which part of an arbitrary measured reflection is a physical
cable; the specified reference-plane shifts are part of the model.

## Zero-frequency behavior

The default `ScatteringLimit()` evaluates the block at zero frequency.
When measured data does not reach DC, or a formula is singular there even
though its limit exists, provide `dcmodel`: `OpenDC()`, `ShortDC()`,
`ThroughDC()`, or `ScatteringDC(S0)`. This controls only the DC rows;
nonzero-frequency evaluation still uses the original provider.

## Noise models

A passive block emits a noise wave with the commutator `K = I-S*S'` and
the covariance `(nbar + 1/2)K`, `K/2` in the vacuum.
`Passive()` uses the analysis temperature, `ThermalEquilibrium(T)` supplies
its own, and `Lossless()` asserts that the block emits no noise. A rational
block retains its noise channels even when its loss is very small.

An active block must state its noise with `NoiseCovariance(V)`. `V` is the
symmetrized added covariance in quanta, the units of `Cnoise`, where the
vacuum counts as one half. It is the noise the block adds, not the noise of
what drives it: the fields entering its ports come from the rest of the
circuit, which counts them where they arise. For an ordinary frequency-preserving block it
must satisfy

```math
V-K/2\succeq0,\qquad V+K/2\succeq0,\qquad K=I-SS^\dagger.
```

For power gain `G` from port 1 to port 2, the minimum corresponding added
noise is `(G-1)/2` in `V[2,2]`. An input-referred added noise of `nadd`
photons is represented by `V[2,2]=G*nadd`; a datasheet noise temperature
`T_N` is `nadd = noisequanta(w, T_N)`. A passive thermal model is
`V=(thermaloccupation(w,T) + 1/2)*K`.

These checks enforce the model's noise contract. They do not establish
its accuracy against a measured device. See [noise conventions](conventions.md#Noise-normalization-and-temperature).

## Pumped devices

[`LinearizedScattering`](@ref) represents a device about a periodic pumped
state. Its harmonic transfer function `H_k(nu)` maps a wave at `nu` to
one at `nu+k*wp`. It retains frequency conversion and the interaction of
the idlers with the surrounding circuit.

The following standalone recipe builds a block from a lossless JPA model.
It is a larger calculation than the small examples above:

```julia
using JosephsonCircuits
jpa = Circuit([
    (:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(1e-12)),
])
wp = (2pi*4.75001e9,)
ws = 2pi .* (4.5:0.01:5.0) .* 1e9
device = hbsolve(ws, wp, [(mode = (1,), port = 1, current = 0.003e-6)],
    (8,), (16,), jpa; returnCnoise = true)
@assert device.nonlinear.solverinfo.converged
block = LinearizedScattering(device.linearized, wp[1];
    noise = NoiseCovariance(device.linearized.Cnoise))
chain = Circuit([(:p1, 1, 0, Port(1)), (:amp, 1, block)])
response = hbsolve(ws, wp, [], (8,), (16,), chain)
```

### Excitation and phase

Include the block's pump frequency in the enclosing HB analysis. When the
block is the only pumped element, it needs no explicit pump source: its
conversion is already in the supplied data. `phase` shifts the pump phase
relative to that data.

The strong-drive pump solve supports coupling between its retained modes.
It rejects a driven calculation when the block also couples a retained
mode to another's conjugate. The zero-drive operating point does not have
that restriction on the subsequent linearized sweep. This restriction
concerns a block's strong-drive embedding, not small-signal conversion.

### Fit for transient use

`RationalScattering(block, npoles; band=...)` fits the harmonic transfer
functions for time-domain use; `band` is in Hz. Choose a band that covers
the signals and conversion products needed by the calculation. A long
traveling-wave device usually also needs delay removal. Outside the
fitted band the filters extrapolate, so convergence with fitting band and
order must be checked against the source data.

An `envelope(t)` scales the conversion terms. It does not model the
underlying pump turning on: the unconverted response, filters, and stated
covariance remain those of the fitted operating point. Simulate the
nonlinear circuit when pump turn-on dynamics are the question.

For noise with a ramped conversion envelope, read the result after the
conversion has been on longer than the block's memory. Without an
envelope, the noise prehistory is taken from the unconverted response;
see [initial fluctuations](transientnoise.md#Baths-and-the-initial-fluctuations).

### Correlated noise and fit tolerances

A lossy pumped block can use `NoiseCovariance(device.linearized.Cnoise)`
from a solve with `returnCnoise=true`. Its covariance contains correlations
between modes, including conjugate idler coordinates. Data supplied at
both a frequency and its conjugate must obey the corresponding transpose
symmetry. The block's `atol`, together with the covariance tolerance,
controls these checks.

A fitted block has a different commutator from the original data. Its noise
is completed to satisfy the commutation relations of the fitted transfer
functions. This completion uses padded frequency ladders so the declared
model does not change just because a sweep retains fewer modes.
`tol` limits scattering fit error and `noisetol` limits the noise added by
the completion. A physical completion does not certify the fitted device's
out-of-band accuracy. Details are in [implementation notes](implementation.md#Pumped-block-noise).
