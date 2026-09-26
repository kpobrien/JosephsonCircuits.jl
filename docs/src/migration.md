# Migrating circuit definitions

This page describes the typed-circuit interface on the development branch.
Use the version selector when comparing with a registered release; the
changes below are not a promise about a particular future release number.

## Tuple netlists and port terminations

The older `(name, node1, node2, value)` netlist infers component types from
name prefixes. It remains readable with a deprecation warning. New code
should use explicit components:

```julia
using JosephsonCircuits
circuit = Circuit([
    (:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:c1, 1, 0, Capacitor(1e-12)),
    (:jj, 1, 0, JosephsonJunction(1e-9)),
])
```

A legacy port adopts the resistor across it as its termination. When
converting that pair to `Port(1; Z0=50.0)`, remove the adopted resistor.
Keeping it would add a second physical load and an internal noise source.
Retain resistors that model actual device dissipation.

## Frequency-dependent values

Replace expressions using the solver's deprecated `symfreqvar` keyword
with `FrequencyDependent` values. The callable receives angular frequency
in rad/s:

```julia
wc = 2pi*10e9
resistor = Resistor(FrequencyDependent(w -> 50.0*(1 + im*w/wc)))
```

Frequency-dependent leaves also combine with parameter expressions.
`FrequencyDependent(identity)` represents the angular frequency when that
form is more convenient.

## Parameter sweeps and sensitivities

`hbcache` and `designsensitivities` take a typed circuit and parameter
definitions, rather than a function that rebuilds the circuit at each
point. Write changing values as parameters such as `Capacitor(:Cc)`, then
supply `Dict(:Cc => 100e-15)`. See the
[sensitivity example](recipes/sensitivities.md) and
[cache example](performance.md#Reuse-across-a-sweep-of-values).

## Node ordering and topology

Choose node order when compiling a typed circuit:

```julia
compiled = compile(circuit; sorting = :number)
```

Pass the compiled circuit to the solver. The solver-level `sorting`
keyword remains available for legacy tuple netlists.

A compiled circuit carries its own topology. Solver entry points no
longer take a separate circuit graph. `calccircuitgraph` remains available
for inspecting topology and its diagnostics.

## Noise in quanta

Noise covariances are symmetrized and counted in quanta, the vacuum being
half a photon, in harmonic balance, the transient and the network cascade.
For the same physics, `Cnoise`, the data of a [`NoiseCovariance`](@ref),
`calcCnoise`, and the noise correlations of `connectS` and `solveS` are
half their previous values. Halve a covariance written for the previous
convention before passing it to `NoiseCovariance`. `S`, `Snoise`, `CM` and,
with the ports at zero temperature, `QE` are unchanged.
[`thermaloccupation`](@ref) returns the occupation `nbar`. The default
`noisetol` of a [`RationalScattering`](@ref) fit halves with the noise it
bounds, so the same fits are accepted; halve a `noisetol` given explicitly
to keep its meaning.

The linearized solvers return `nbar`, the occupation of the wave leaving
each port mode, by default. A sweep which asked for `S` alone with
`returnQE = false, returnCM = false` adds `returnnbar = false`.

The quadrature functions of the quantum optics module take `hbar`, one by
default, with the vacuum's covariance `(hbar/2)I`; `hbar = 2` gives their
previous values. The ladder functions count the vacuum as `I/2`.
