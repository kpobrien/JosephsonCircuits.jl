# From a pumped JPA to a transient scattering model

Solve a lossy Josephson amplifier, export its small-signal conversion and
noise, fit a causal model, and measure its gain and noise in time. The
checks compare both the fit and the transient result with the original
circuit. This example requires `JosephsonCircuits` and `LinearAlgebra`.

## Solve the original operating point

The 20 kΩ resistor represents internal dissipation. Port and resistor baths
are at their default zero temperature. The pump current is a complex HB
Fourier coefficient; a real sinusoid would have twice this peak amplitude.

```text
 1                  2
 o------[cc]--------o--------+--------+
 |                  |        |        |
[p1]              [jj]     [cj]    [loss]
 |                  |        |        |
 o------------------o--------+--------+
 0
```

```@example fittedamp
using JosephsonCircuits, LinearAlgebra
jpa = Circuit([
    (:p1, 1, 0, Port(1)),
    (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1e-9)),
    (:cj, 2, 0, Capacitor(1e-12)),
    (:loss, 2, 0, Resistor(2e4)),
])
fp = 4.75e9
wp = (2pi*fp,)
sources = [(mode = (1,), port = 1, current = 0.00565e-6)]
ws = 2pi .* collect(range(4e9, 5.5e9; length = 76))
device = hbsolve(ws, wp, sources, (4,), (8,), jpa;
    atol = 1e-14, returnCnoise = true, returnVout = true)
@assert device.nonlinear.solverinfo.converged
nothing # hide
```

`Cnoise` is the internal added covariance used to export the device;
`Vout` also includes noise entering through the external port. Substituting
`Vout` for `Cnoise` would count that port noise twice in the new circuit.

## Export and fit

`LinearizedScattering` retains the pump frequency and sideband conversion.
`RationalScattering` supplies states that can be advanced in time. Its fit
uses the tabulated signal and sideband frequencies by default. Here the
device is lumped and needs no reference-plane delay removal; see the
[scattering guide](../scattering.md) for fitting measured cables.

```@example fittedamp
rise(t) = t <= 0 ? 0.0 : t >= 10e-9 ? 1.0 : (1 - cospi(t/10e-9))/2
exported = LinearizedScattering(device.linearized, wp[1];
    noise = NoiseCovariance(device.linearized.Cnoise), envelope = rise)
fitted = RationalScattering(exported, 10; tol = 1e-3, noisetol = 5e-3)
model = Circuit([(:p1, 1, 0, Port(1)), (:amp, 1, fitted)])
nothing # hide
```

The exported model describes **perturbations about this one pump operating
point**. Its output does not contain the original pump carrier, and it
cannot predict depletion or changes of pump loading. Recompute HB and
refit when changing the pump or its embedding environment.

The envelope gates conversion terms; it does not solve the original
junction's pump turn-on. The static response and noise model still belong
to the exported operating point. It lets this example begin with
unconverted stationary prehistory and wait for a settled measurement.
For physical pump transients, simulate the original nonlinear circuit as
in the [pulsed-noise recipe](pumped-noise.md).

`tol` controls the scattering fit and `noisetol` the noise completion in
quanta. Neither bounds the final circuit error. A pumped amplifier is
active, so passive-fit criteria are not an amplifier stability test;
check the original operating point with [stability analysis](../stability.md).

## Check frequencies outside the fitting grid

Use the same converged pump state to evaluate the original circuit at
three additional signal frequencies. HB evaluates the exported conversion
fully on, without its transient envelope. The fitted circuit needs no
classical pump source because its conversion coefficients already contain
the pump dependence.

```@example fittedamp
checkw = 2pi .* [4.61e9, 4.67e9, 4.73e9]
reference = hblinsolve(checkw, jpa; nonlinear = device.nonlinear,
    Nmodulationharmonics = (4,), returnCnoise = true, returnVout = true)
fitcheck = hbsolve(checkw, wp, [], (4,), (8,), model;
    returnCnoise = true, returnVout = true)
fiterror = maximum(abs.(fitcheck.linearized.S((0,), 1, (0,), 1, :) .-
    reference.S((0,), 1, (0,), 1, :)))
@assert fiterror < 1e-3
fiterror
```

This checks signal reflection inside the sampled band. For a different
measurement, check its conversion channels and frequencies too; a fit
accepted on its training data does not justify extrapolation.

## Measure settled gain and noise in time

There is no mean probe here: tangent/noise propagation extracts the
infinitesimal response about the zero perturbation trajectory. Record
checkpoints for the response calculation, wait 40 ns, and measure a 10 ns
rectangular Fourier window at 4.7 GHz.

```@example fittedamp
f, settle, T = 4.7e9, 40e-9, 10e-9
function measure(dt)
    sol = transientsolve(transientproblem(model), (0.0, settle + T - dt);
        dt, record = :checkpoints)
    times = sol.times[round(Int, settle/dt) + 1:end]
    plan = transientquantumplan(sol, times, [2pi*f])
    # The even pump ladder that mixes into this settled Fourier bin.
    frequencies = sort!(abs.([2pi*(f + 2k*fp) for k in -2:2]))
    transientnoise(sol, plan; frequencies,
        weights = fill(2pi/T, length(frequencies)), inputs = plan,
        commutationrtol = 3e-3)
end
noise = measure(2.5e-12)
@assert noise.diagnostics.passed
nothing # hide
```

The optional `inputs` gain uses periodically extended input modes with
stationary prehistory. Together with a settled output window, this makes
it comparable to monochromatic HB. A probe applied only within a finite
input window belongs in `transientgain` and generally gives another result.

This sparse bath quadrature is specific to the settled Fourier-bin check:
both the carrier and conversion frequencies fall on the window grid, and
only the retained even pump ladder contributes in the periodic limit.
It is not a general noise quadrature for a ramp or a short pulse. Use a
broad bath band and refine its spacing and cutoff for those measurements,
as in the [finite-window example](pumped-noise.md).

```@example fittedamp
original = hblinsolve([2pi*f], jpa; nonlinear = device.nonlinear,
    Nmodulationharmonics = (4,), returnVout = true, returnCnoise = true)
S = original.S((0,), 1, (0,), 1, 1)
V = real(original.Vout((0,), 1, (0,), 1, 1))
# The package's (X,P) convention for this phase-insensitive signal mode.
G = [real(S) imag(S); -imag(S) real(S)]
relative_error(a, b) = norm(a - b)/norm(b)
errors = (gain = relative_error(noise.gain, G),
    covariance = relative_error(noise.covariance, V*Matrix{Float64}(I, 2, 2)))
@assert maximum(values(errors)) < 1e-4
errors
```

Both errors are a few parts per million for these settings. The covariance
includes external port vacuum and the fitted internal noise; vacuum alone
would be `I/2`. This comparison checks the combined fit, sideband
truncation, settling, and time integration at one operating point, rather
than providing a general accuracy guarantee.

```@example fittedamp
coarse = measure(5e-12)
changes = (gain = relative_error(coarse.gain, noise.gain),
    covariance = relative_error(coarse.covariance, noise.covariance))
@assert coarse.diagnostics.passed
@assert maximum(values(changes)) < 1e-4
changes
```

For a new device, also increase the settling time, the HB harmonic limits,
and the fitting sample density/order until the measured observables stop
changing at the required accuracy. Passing the commutator check is
necessary evidence about noise consistency, not an error bound on gain.
