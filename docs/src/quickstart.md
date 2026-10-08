# A first simulation

This example computes the reflection gain of a current-pumped Josephson
parametric amplifier (JPA). It needs only `JosephsonCircuits`; plotting is
optional. For installation, see the [home page](index.md#Installation).

## Build the circuit

The port owns a matched 50 Ω termination. A coupling capacitor connects it
to a junction shunted by a capacitor. Node `0` is ground; node `2` is the
junction node.

```@example quickstart
using JosephsonCircuits

circuit = Circuit([
    (:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1000e-12)),
    (:cj, 2, 0, Capacitor(1000e-15)),
])
nothing # hide
```

`JosephsonJunction(Lj)` takes the small-signal junction inductance in
henries. A port already supplies its termination; adding a resistor
across it adds another physical load. See [ports](circuits.md#Ports-and-sources).

## Solve the pump and sweep the signal

Frequencies passed to harmonic balance are angular frequencies in rad/s.
The pump's `current` is a Fourier coefficient in amperes: a real value
`Icoeff` represents a cosine of peak current `2Icoeff`.

```@example quickstart
fp = 4.75001e9                       # pump frequency, Hz
Icoeff = 0.00565e-6                  # pump Fourier coefficient, A
ws = 2pi .* (4.5:0.01:5.0) .* 1e9   # signal angular frequencies, rad/s
wp = (2pi*fp,)                      # one independent pump
sources = [(mode = (1,), port = 1, current = Icoeff)]

sol = hbsolve(ws, wp, sources, (8,), (16,), circuit)
@assert sol.nonlinear.solverinfo.converged
S11 = sol.linearized.S(outputmode = (0,), outputport = 1,
    inputmode = (0,), inputport = 1, freqindex = :)
gain_dB = 10 .* log10.(abs2.(S11))
@assert all(isfinite, gain_dB) # hide
round(maximum(gain_dB); digits = 1)
```

The gain peaks near the pump frequency. Here the input and output
frequencies are equal, so photon gain and power gain coincide.
`(8,)` sets the retained modulation harmonics for the signal calculation;
`(16,)` sets the retained harmonics of the strong pump. These are numerical
resolution choices, not device parameters.

Convergence of the nonlinear solver checks its residual on this grid.
Before using the result quantitatively, increase the pump and modulation
harmonic limits and compare the gain. Refine the nonlinear evaluation grid
separately if necessary. See [convergence](harmonicbalance.md#Checking-convergence).

To plot the result, install `Plots` and continue in the same session:

```julia
using Plots
plot(ws ./ (2pi*1e9), gain_dB; xlabel = "Signal frequency (GHz)",
    ylabel = "Reflection gain (dB)", label = "JPA")
```

## Apply the same pump in time

A `TransientSource` returns instantaneous current, so the waveform includes
the factor of two. The smooth ramp starts the circuit at rest.

```@example quickstart
ramp(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : (1 - cospi(t/2e-9))/2
pump(t) = 2Icoeff*ramp(t)*cospi(2fp*t)
problem = transientproblem(circuit; sources = [TransientSource(1, pump)])
solution = transientsolve(problem, (0.0, 5e-9); dt = 2.5e-12)
@assert all(isfinite, solution.voltage) # hide
(size(solution.voltage, 1), first(solution.times), last(solution.times))
```

The voltage array has one row per port and one column per saved time.
This short record illustrates turn-on; it is not a settled measurement of
the frequency-domain gain. A transient gain calculation also needs a probe
or a tangent response. Continue with the [transient guide](transient.md).

## Choose the next analysis

| Question | Analysis |
|---|---|
| What is the steady pump response? | `hbnlsolve` |
| What weak signals and idlers does it amplify? | `hblinsolve` about that operating point, or `hbsolve` for both steps |
| Do perturbations about this state grow? | [`hbstability`](stability.md), after checking the operating point |
| How does a pulse propagate or deplete the pump? | `transientsolve` |
| What is the incremental response about a recorded trajectory? | `transienttangent`, `transientadjoint`, or `transientgain` |
| What noise is measured in a time window? | `transientnoise` with a temporal-mode plan |

See [conventions](conventions.md) for current amplitudes, signed modes,
wave normalization, and the units of the returned states.
