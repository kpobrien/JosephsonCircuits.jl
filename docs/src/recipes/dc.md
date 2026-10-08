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

The port's source drives the resistor `rl`; the second example adds the
current source `i1` across it:

```text
 1
 o-------+--------+--------+
 |       |        |        |
[p1]   [rl]     [c1]     [i1]
 |       |        |        |
 o-------+--------+--------+
 0
```

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

The port source injects its current into node 1, the port's first
terminal. A [`CurrentSource`](@ref) in the circuit drives its current
through itself from its first terminal to its second: it draws the
current from the node at its first terminal and delivers it to the node
at its second, the opposite sense. Written from node 1 to ground, the
same current develops the opposite voltage across the same load, and
beside the port source it cancels it:

```@example dc
reversed = Circuit(
    [(:p1, 1, 0, Port(1; Z0 = 50.0)),
     (:i1, 1, 0, CurrentSource(Idc)),
     (:rl, 1, 0, Resistor(150.0)),
     (:c1, 1, 0, Capacitor(1.0e-12))])
alone = hbnlsolve((2*pi*5e9,), (1,), [], reversed;
    dc = true, odd = true, keyedarrays = false)
both = hbnlsolve((2*pi*5e9,), (1,), [(mode = (0,), port = 1, current = Idc)],
    reversed; dc = true, odd = true, keyedarrays = false)
@assert isapprox(alone.dcnodevoltage[1], -expected_voltage; rtol = 1e-10)
@assert abs(both.dcnodevoltage[1]) < 1e-12*expected_voltage
(alone.dcnodevoltage[1], both.dcnodevoltage[1])
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

## Initialize a biased junction in time

A constant drive present at the first transient sample requires a
consistent initial state. This example biases a junction through a series
resistor. The port termination also carries DC, so the injected current
is larger than the junction current. KCL and the zero-voltage junction
relation determine the initial voltage and flux analytically.

```text
 drive             junction
 o-------[r]-------o--------+
 |                 |        |
[p1]             [jj]     [cj]
 |                 |        |
 o-----------------o--------+
 0
```

```@example biasedstate
using JosephsonCircuits
Lj, R, Z0 = 1e-9, 100.0, 50.0
phi0 = JosephsonCircuits.phi0
Ic = phi0/Lj
Ij = 0.3Ic
Idc = Ij*(R + Z0)/Z0
circuit = Circuit([
    (:p1, "drive", 0, Port(1; Z0)),
    (:r, "drive", "junction", Resistor(R)),
    (:jj, "junction", 0, JosephsonJunction(Lj)),
    (:cj, "junction", 0, Capacitor(1e-12)),
])
compiled = compile(circuit)
problem = transientproblem(compiled; sources = [TransientSource(1, Returns(Idc))])

# Ground is omitted. Use the compiled node order, not netlist position.
names = compiled.nodenames[2:end]
flux = [name == "junction" ? phi0*asin(Ij/Ic) : 0.0 for name in names]
voltage = [name == "drive" ? R*Ij : 0.0 for name in names]
initial = transientstate(problem; flux, voltage)
solution = transientsolve(problem, (0.0, 1e-9); dt = 5e-12,
    initialstate = initial, record = :states)
@assert maximum(abs.(solution.voltage .- R*Ij)) < 1e-12
solution.voltage[:, end]
```

Transient flux is in webers; HB `nodeflux` is divided by `phi0`. The
junction begins at phase `asin(Ij/Ic)` and zero voltage. The drive node
has voltage `R*Ij`; its flux then grows linearly. This is compatible with
HB's separate average-voltage coordinate:

```@example biasedstate
hb = hbnlsolve((2pi*1e9,), (1,), [(mode = (0,), port = 1, current = Idc)],
    compiled; dc = true, atol = 1e-12)
@assert hb.solverinfo.converged
@assert isapprox(hb.dcnodevoltage, voltage; rtol = 1e-9)
hb.dcnodevoltage
```

This is an analytic initialization for this circuit, not a general HB to
transient adapter. For a driven periodic orbit, reconstructing a state also
requires the correct time origin, rates, auxiliary currents, and any line
or rational-block history. Continuing from `transientstate(solution)`
preserves that transient history.

The default zero state fails here because the uncapacitated drive node
must satisfy algebraic KCL immediately. Check that the solver reports the
inconsistent initial condition:

```@example biasedstate
failure = try
    transientsolve(problem, (0.0, 1e-9); dt = 5e-12)
    nothing
catch err
    err
end
@assert failure isa ArgumentError
@assert occursin("initial state", sprint(showerror, failure))
sprint(showerror, failure)
```

Increasing the iteration limit cannot repair this initial condition.
Supply the consistent state above, or start at equilibrium and smoothly
ramp the source to let the circuit approach the intended branch. A ramp
can select another branch in a multistable circuit, so inspect the settled
state before using it.
