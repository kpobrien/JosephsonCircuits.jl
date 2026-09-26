# Direct current

Distinguish static flux from average voltage using a resistor carrying DC. The executable example checks Ohm's law without adding an artificial inductive path to ground.

The harmonic balance state is periodic node flux, and a voltage is its time
derivative, so a mode of frequency `w` has `V = i*w*phi0*phi`. At zero
frequency that is zero, which is right for a capacitor and wrong for a
resistor: a resistor would be an open circuit at DC and a current source
driving one could not develop `I*R`.

The solver carries the missing coordinate, the average node voltage,
separately. A finite inductor or a zero-voltage Josephson junction has zero
average voltage across it, so the average voltage is constant on each
connected group of inductors and junctions, and the direct current problem
reduces to those coordinates alone. They are carried as a small block of
unknowns beside the periodic state, with those equations as their rows, and
only when some direct current is injected: with none every average voltage
is zero, the periodic system is exact as it stands, and a circuit with no
direct current drive pays nothing.

`hbnlsolve` returns the result as `dcnodevoltage`, in volts, indexed like
`nodeflux`: ground is excluded, so the first entry is the first real node.
It is a vector of zeros when the circuit has a zero frequency mode and no
direct current is drawn, and `nothing` only when the analysis has no zero
frequency mode. It is not the same thing as the zero frequency entry of
`nodeflux`, which remains the static periodic flux that sets inductor
currents and junction phases.

```@example dc
using JosephsonCircuits

# a current source into a resistor, with a capacitor which carries no
# direct current
Idc = 1.0e-6
circuit = Circuit(
    [(:p1, 1, 0, Port(1; Z0 = 50.0)),
     (:rl, 1, 0, Resistor(150.0)),
     (:c1, 1, 0, Capacitor(1.0e-12))])

sol = hbnlsolve((2*pi*5e9,), (1,), [(mode = (0,), port = 1, current = Idc)],
    circuit; dc = true, odd = true, keyedarrays = false)

# the port environment in parallel with the load, so V = I*(50 || 150)
expected_voltage = Idc*inv(1/50.0 + 1/150.0)
@assert sol.solverinfo.converged
@assert isapprox(sol.dcnodevoltage[1], expected_voltage; rtol = 1e-10)
(sol.dcnodevoltage[1], expected_voltage)
```

A finite inductor across a resistor sets its average voltage to zero. A
Josephson junction does the same on a zero-voltage branch, provided that
branch exists at the imposed bias. A junction driven into a running-phase
state is outside this periodic-flux DC model. A resistor between distinct
voltage groups can carry current with a voltage difference of `I*R`.

Two things are refused rather than approximated. A component whose
conductance at zero frequency is not finite and real has no direct current
behaviour to use -- a frequency dependent resistance whose limit at DC is
complex or unbounded, for instance -- and is reported with the entry named.
And a direct current with nowhere to go, injected into a group of nodes
which no resistor, inductor or junction connects to anything else, has no
bounded solution and is reported as such before the solve starts.

A scattering block carries direct current according to its zero frequency
limit: by default the limit of its own data (`ScatteringLimit()`), which a
constant or tabulated block reaching zero frequency supplies, and
otherwise the model stated with its `dcmodel` keyword, `OpenDC()`,
`ShortDC()`, `ThroughDC()` or `ScatteringDC(S0)`. A block whose limit is
a short or a through constrains the direct voltages of its ports instead
of conducting between them, which the solver handles by an explicit
direct current block rather than by the elimination above.
