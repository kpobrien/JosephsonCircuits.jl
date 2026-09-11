# Defining a circuit

A circuit is a list of component instances and the connections between
their terminals. It is written in one of two forms, and both build the same
`Circuit` object.

## The netlist form

Each entry is a tuple of the instance name, the node of every terminal in
order, and the component: `(name, nodes..., component)`, as a line of a
SPICE netlist. Names are symbols or strings, nodes are integers, strings or
symbols, and node `0` (or `"0"`) is ground. Entries which name the same node
share a net, and the nets take the node names, so they can be found again in
the outputs.

```julia
using JosephsonCircuits

R = 50.0
Cc = 100.0e-15
Lj = 1000.0e-12
Cj = 1000.0e-15

# a Josephson parametric amplifier: a port, a coupling capacitor, and a
# junction shunted by a capacitor, with node 0 the ground
circuit = Circuit(
    [(:p1, 1, 0, Port(1; Z0 = R)),
     (:cc, 1, 2, Capacitor(Cc)),
     (:jj, 2, 0, JosephsonJunction(Lj)),
     (:cj, 2, 0, Capacitor(Cj))])
```

The components are typed: `Capacitor`, `Inductor`, `Resistor`,
`JosephsonJunction`, `Port` (with its reference impedance `Z0`),
`MutualInductor`, `NonlinearInductor`, `ScatteringParameters` for a block
described by its scattering matrix, and a `Circuit` with an interface as a
subcircuit. A component value may be a number, a complex number (a
capacitor with dielectric loss), or a `FrequencyDependent` function of the
mode frequency.

An entry lists one node per terminal, so the form is not limited to two
terminal elements. A subcircuit instance lists its pins in the order they
were declared, a scattering parameter block lists the signal terminal of
each port (or both terminals of each port when the block is not grounded),
and a mutual inductor, which couples two inductor branches rather than
nets, names the two inductors in place of nodes:

```julia
(:k1, :l1, :l2, MutualInductor(0.9))
```

## The connection-group form

The same circuit written as a list of named components and a list of
connection groups, each group naming the terminals which share a net. A
terminal is `(instance, number)`; `Ground` may appear in any group, and may
also be declared as a component, `:gnd => Ground()`, and referred to
through its single terminal.

```julia
circuit = Circuit(
    [:p1 => Port(1; Z0 = R),
     :cc => Capacitor(Cc),
     :jj => JosephsonJunction(Lj),
     :cj => Capacitor(Cj),
     :gnd => Ground()],
    [[(:p1, 1), (:cc, 1)],
     [(:cc, 2), (:jj, 1), (:cj, 1)],
     [(:p1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
```

This is the form the netlist form expands to. It is what to use when a
connection is not a node list: the bundled port views of scattering blocks
in pair connections, or nets named explicitly with `Net`. The two forms may
be mixed freely across a hierarchy, since either produces a `Circuit`.

## Subcircuits

A `Circuit` given an interface through the `pins` keyword is a component,
and is instanced in either form like any other. The pins map an interface
pin number to a terminal of an inner component. The examples below build
traveling wave amplifiers from unit cells and a snake amplifier from
hierarchical subcircuits this way; a subcircuit instanced many times, such
as a unit cell, is defined once, and the flattened circuit names its inner
nets by path.

```julia
# one unit cell of a transmission line, exposed through pins 1 and 2
cell(Lj, Cj, Cg) = Circuit(
    [(:jj, 1, 2, JosephsonJunction(Lj)),
     (:cj, 1, 2, Capacitor(Cj)),
     (:cg, 1, 0, Capacitor(Cg))];
    pins = [1 => (:jj, 1), 2 => (:jj, 2)])

# three cells in a chain between two ports
line = Circuit(
    [(:p1, 1, 0, Port(1)),
     (:cell1, 1, 2, cell(1e-9, 50e-15, 40e-15)),
     (:cell2, 2, 3, cell(1e-9, 50e-15, 40e-15)),
     (:cell3, 3, 4, cell(1e-9, 50e-15, 40e-15)),
     (:p2, 4, 0, Port(2))])
```

The older netlist of `(name, node1, node2, value)` tuples with the
component type given by the prefix of the name, `("C1", "1", "0", 1e-12)`,
is still read; see the `Circuit` docstring.


## Nonlinear elements and their current-phase relations

A [`JosephsonJunction`](@ref) is the sinusoidal relation
`I(φ) = (phi0/Lj)*sin(φ)`, and it is what almost every circuit uses. An
element whose relation is something else is a
[`NonlinearInductor`](@ref), written as its small signal inductance and a
relation of unit slope at zero:

```julia
using JosephsonCircuits
# the effective relation of a SNAIL, biased away from its symmetric point:
# the quadratic term is what makes it a three wave mixer
snail = NonlinearInductor(1e-9, PolynomialCPR([1.0, 0.3, -1/6]))
```

A [`PolynomialCPR`](@ref) is given by the coefficients of its expansion,
`f(φ) = c[1]*φ + c[2]*φ^2 + ...`, with `c[1] = 1` so that the `L0` of the
element is the small signal inductance. It is the way to write an element
whose junctions you do not want to wire up: a SNAIL, a SQUID, a Quarton,
a kinetic inductor, or an array of `N` junctions in series, which divides
the phase and so has the relation `N*sin(φ/N)`.

The element is a junction to everything else. It makes the same branch,
enters the same matrices, and is indexed with the junctions, so a circuit
which mixes the two kinds is ordinary. What differs is the relation the
solver evaluates at the junction phases, and its derivative, which is
where `cos` would otherwise stand.

Two things follow from a polynomial not being a sine. It is not bounded,
so nothing warns that an element is past the range its coefficients were
fitted on, and it is not band limited: a term of degree `d` generates
harmonics to `d` times the drive, which the harmonic count has to cover.
Both are yours to judge.

Both solvers evaluate any of these relations. Harmonic balance takes it
in the residual, the Jacobian, the Hessian and the pump modulation of the
linearized system; the transient solver steps it, and its tangent, its
adjoint and the linearization its noise is taken about all read the same
derivative.


## Scattering blocks, transmission lines and fitted data

A [`ScatteringParameters`](@ref) block is a multiport given by its
scattering matrix: a constant matrix, a callable of angular frequency, a
table of frequencies and matrices, or a Touchstone file, at its own
reference impedances, with the default `grounded = true` tying every
reference terminal to ground. The harmonic balance solvers evaluate a
block at every frequency they need; the time domain solver takes a block
with a constant real matrix, an attenuator, a circulator, a through,
a short or an open, and a [`RationalScattering`](@ref) block, a real
state space realization `S(s) = D + C (s I - A)^(-1) B`, which is the
form a fit of measured or simulated data takes.

A tabulated block is interpolated with the cubic spline through each
entry's samples and, by default, refuses a frequency outside its band.
The harmonic balance solvers place mixing products at sums and
differences of the pump harmonics and the signal, which measured data
often does not cover, and the noise of a lossy block reads the
dissipation `I - S S'`, in which an error of the data or its interpolant
appears roughly doubled. Measured data meant for those solvers is
therefore best fitted once with [`RationalScattering`](@ref): the fit
extrapolates as a passive rational function, is passive at every
frequency by construction, and the same block runs in the frequency and
the time domain solvers.

```julia
# a fit of tabulated data at as many poles as it might need; the poles it
# does not need are dropped
data = ScatteringParameters((2pi .* frequencies, S); nports = 2, zref = 50.0)
fitted = RationalScattering(data, 8)

# a hand written realization: a series inductor L between 50 ohm ports
a = 2*50.0/L
inductor = RationalScattering(fill(-a, 1, 1), [1.0 -1.0], -a .* [1.0; -1.0;;], Matrix(1.0I, 2, 2); zref = 50.0)
```

A rational block is validated as stable and passive at construction by
its largest singular value over every frequency, which finds a peak
however narrow, and a fit is made passive where it strays. A
[`TransmissionLine`](@ref) is an ideal lossless line of a characteristic
impedance and a length, with a phase velocity that defaults to the speed
of light; in the frequency domain it is the exact line, in time the
method of characteristics with a history of the waves. A delay is not a
rational function, so a lossy cable is a line in cascade with a fit of
its data with the delay removed.

## Noise models and temperatures

Every dissipative element is a bath at a temperature: a
[`Resistor`](@ref) at its `temperature` keyword or the analysis default,
a port's termination at the port's, and a scattering block at its noise
model, [`Passive`](@ref) taking the analysis default,
[`ThermalEquilibrium`](@ref)`(T)` its own temperature, and
[`Lossless`](@ref) asserting that the block emits nothing, which is
validated where it can be. A lossy block emits the noise wave of
Bosma's relation, of covariance `I - S S'`; nothing is inferred about a
rational block's loss, which keeps its noise channels however small it
is. The same models and temperatures set the noise of the linearized
solver and of the time domain solver, so the two compare on the same
circuit.
