# Defining a circuit

A [`Circuit`](@ref) contains named component instances and their
connections. Use a netlist for a circuit with named nodes, or connection
groups to connect terminals directly. Both forms produce the same circuit
model and can be mixed across a hierarchy.

## The netlist form

Each entry is `(name, nodes..., component)`. The node arguments follow the
component's terminal order. Instance names can be symbols or strings;
node names can also be integers. Node `0`, `"0"`, or `Ground` is ground.
Entries naming the same node are connected.

```@example circuit
using JosephsonCircuits
R, Cc, Lj, Cj = 50.0, 100e-15, 1000e-12, 1000e-15
circuit = Circuit([
    (:p1, 1, 0, Port(1; Z0 = R)),
    (:cc, 1, 2, Capacitor(Cc)),
    (:jj, 2, 0, JosephsonJunction(Lj)),
    (:cj, 2, 0, Capacitor(Cj)),
])
nothing # hide
```

The example is a JPA: a port couples through `cc` to a junction shunted by
`cj`. Component values use SI units. `JosephsonJunction(Lj)` takes its
small-signal inductance; [`IctoLj`](@ref) converts a critical current to
that inductance.

## Ports and sources

`Port(1; Z0=50.0)` identifies port 1 and supplies a matched 50 Ω external
termination across its two terminals. It does not connect either terminal
to ground unless the circuit does so. An extra resistor across the port is
an additional load, with its own internal noise.

Use `termination=nothing` when the port should add no physical loading.
This is a current-source and impedance-probe boundary, with a different
reflection interpretation; see [`Port`](@ref). Transient noise baths
require matched port terminations.

Drive amplitudes belong to the analysis. HB uses Fourier coefficients;
`TransientSource` uses instantaneous current. A port source injects its
current into the port's first (positive) terminal. A
[`CurrentSource`](@ref) component drives its constant current through
itself from its first terminal to its second: it draws the current from
the node at its first terminal and delivers it to the node at its second,
the opposite sense of a port source. A port and a `CurrentSource` written
on the same nodes in the same order therefore drive in opposite senses.
See [current conventions](conventions.md#Current-amplitudes).

## The connection-group form

A connection group lists terminals that share a net. `(instance, number)`
selects a terminal; `Ground` can occur directly in a group or as a named
component. Continuing the setup above:

```@example circuit
connected = Circuit(
    [:p1 => Port(1; Z0 = R), :cc => Capacitor(Cc),
     :jj => JosephsonJunction(Lj), :cj => Capacitor(Cj), :gnd => Ground()],
    [[(:p1, 1), (:cc, 1)],
     [(:cc, 2), (:jj, 1), (:cj, 1)],
     [(:p1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
# Both descriptions give the same reflection coefficient.
a = hblinsolve([2pi*4e9], circuit; keyedarrays = false)
b = hblinsolve([2pi*4e9], connected; keyedarrays = false)
@assert isapprox(a.S, b.S; rtol = 1e-12)
nothing # hide
```

Use this form when connections are easier to express through component
interfaces than through node names. [`Net`](@ref) names a connection
explicitly; [`PortRef`](@ref) and [`PinRef`](@ref) address interfaces.

## Mutual inductors

A mutual-inductor entry names two inductor instances rather than nodes:

```@example mutual
using JosephsonCircuits
coupled = Circuit([
    (:p1, 1, 0, Port(1)),
    (:l1, 1, 0, Inductor(1e-9)),
    (:l2, 2, 0, Inductor(2e-9)),
    (:p2, 2, 0, Port(2)),
    (:k1, :l1, :l2, MutualInductor(0.9)),
])
nothing # hide
```

For positive coupling, currents entering the first terminal of each
inductor add flux to both. Reversing one inductor's terminals, or changing
the sign of the coefficient, reverses the coupling orientation.

## Subcircuits

Give a circuit a `pins` interface to use it as a component. Each interface
pin maps to a terminal of an internal instance. A netlist instance lists
the external nodes in the declared pin order.

```@example subcircuit
using JosephsonCircuits
cell(Lj, Cj, Cg) = Circuit([
    (:jj, 1, 2, JosephsonJunction(Lj)),
    (:cj, 1, 2, Capacitor(Cj)),
    (:cg, 1, 0, Capacitor(Cg)),
]; pins = [1 => (:jj, 1), 2 => (:jj, 2)])

line = Circuit([
    (:p1, 1, 0, Port(1)),
    (:cell1, 1, 2, cell(1e-9, 50e-15, 40e-15)),
    (:cell2, 2, 3, cell(1e-9, 50e-15, 40e-15)),
    (:cell3, 3, 4, cell(1e-9, 50e-15, 40e-15)),
    (:p2, 4, 0, Port(2)),
])
compiled = compile(line)
nothing # hide
```

Repeated instances can share one subcircuit definition. Elaboration
flattens the hierarchy and qualifies internal names by their instance
paths. The [JTWPA](recipes/traveling-wave.md) and
[snake amplifier](recipes/lesa.md) examples use this pattern at larger
scales.

## Component values and design parameters

A component value can be a number or a parameter such as `:Lj`. Supply a
dictionary of parameter values when solving. For expressions involving
several parameters, use `JosephsonCircuits.@params`; see
[design sensitivities](recipes/sensitivities.md).
[`symbolicmatrices`](@ref) gives the circuit's matrices with the
parameters left as expressions:

```@example parameters
using JosephsonCircuits
jpa = Circuit([
    (:p1, 1, 0, Port(1; Z0 = :R)),
    (:cc, 1, 2, Capacitor(:Cc)),
    (:jj, 2, 0, JosephsonJunction(:Lj)),
    (:cj, 2, 0, Capacitor(:Cj)),
])
matrices = symbolicmatrices(jpa)
(capacitance = matrices.Cnm, conductance = matrices.Gnm)
```

Frequency-domain analyses also accept complex values and
[`FrequencyDependent`](@ref) closures. The callable receives nonnegative
angular frequency in rad/s; negative-frequency modes use conjugate
symmetry. For example, a simple frequency-dependent impedance is:

```@example values
using JosephsonCircuits
wc = 2pi*10e9
load = Resistor(FrequencyDependent(w -> 50.0*(1 + im*w/wc)))
nothing # hide
```

Such a value is not automatically a causal time-domain realization. Use
explicit circuit elements or a [rational scattering model](scattering.md)
for transient simulation.

## Nonlinear elements and their current-phase relations

A [`JosephsonJunction`](@ref) uses `I(φ) = (phi0/Lj)*sin(φ)`, where `φ`
is reduced branch flux and `phi0` is the reduced flux quantum.
[`NonlinearInductor`](@ref) accepts a different current-phase relation.
The inductance sets the current scale; the relation has unit slope at zero.

```@example cpr
using JosephsonCircuits
# An effective asymmetric relation with quadratic and cubic terms.
element = NonlinearInductor(1e-9, PolynomialCPR([1.0, 0.3, -1/6]))
nothing # hide
```

[`PolynomialCPR`](@ref) represents `f(φ) = c[1]*φ + c[2]*φ^2 + ...`,
with `c[1]=1`. It can approximate a biased SNAIL, a SQUID, a kinetic
inductor, or a junction array over a specified phase range. For example,
`N` identical series junctions have normalized relation `N*sin(φ/N)`;
expand that relation to obtain a local polynomial approximation.

For a kinetic inductor with differential inductance
`L(I)=L0*(1+I^2/Istar^2)`, the expansion through fifth order is
`[1, 0, -(IL/Istar)^2/3, 0, (IL/Istar)^4/3]`, where `IL=phi0/L0`.

Check two approximations independently:

1. Keep the simulated phase within the range where the fitted relation is
   physically meaningful. Polynomial extrapolation can predict unphysical currents.
2. Refine the Fourier grids or time step. A degree-`d` polynomial can
   expand the Fourier support of its input phase by a factor of `d`;
   a sine also generates higher harmonics. The drive frequency alone does
   not determine the required grid.

Both solvers use the chosen relation and its derivatives for the
nonlinear response, linearization, and sensitivities.

## Supported analyses

| Representation | Harmonic balance | Transient simulation |
|---|---|---|
| Real constant R, C, L, mutual inductors, junctions, polynomial CPRs | Supported | Supported; initial algebraic constraints must hold |
| Complex or frequency-dependent lumped values | Supported | Use a causal circuit or block realization instead |
| Constant real scattering matrix | Supported | `GaussLegendre()` |
| Tabulated, Touchstone, or callable scattering data | Supported within the provider's domain | Fit a rational model first |
| Rational scattering realization | Supported | `GaussLegendre()` |
| Pumped `LinearizedScattering` | Supported subject to its excitation restrictions | Fit with `RationalScattering`; `GaussLegendre()` |
| Ideal `TransmissionLine` | Supported | `GaussLegendre()`; delay at least one step |

Noise calculations have additional restrictions. In particular, noise
outputs for lossy mutually coupled inductors are rejected; an S-only
frequency-domain calculation does not supply a noise model for that loss.
Transient noise requires a passive equilibrium prehistory and supported
baths. See [quantum noise in time](transientnoise.md).

## Storage and compilation

For a generated circuit, use vectors of entries or connection groups.
Component types may differ, so a heterogeneous description is normal.
Compilation gathers values into concrete numerical arrays for the solver;
the entire topology need not be encoded in a Julia tuple type.

The connection-group constructor retains the supplied collections.
Subsequent elaboration observes edits to them. A compiled circuit or solver
cache represents a compiled structure: rebuild it after changing topology.
For repeated value changes, use the [cache interface](performance.md).

The [implementation notes](implementation.md#Circuit-compilation) describe
parsing, hierarchy reuse, and storage. Custom interface keys must obey
Julia's `isequal`/`hash` contract.

### Interface changes

See [migration](migration.md) for code written for v0.5.4: tuple
netlists, symbolic values, node order, and frequency-dependent
expressions.

## Scattering blocks, transmission lines and fitted data

The [scattering-block guide](scattering.md) covers construction, fitting,
delays, pumped devices, and noise contracts with separate examples.

## Noise models and temperatures

Ports, internal losses and scattering blocks take the temperatures and
noise models of the
[temperature table](conventions.md#Noise-normalization-and-temperature);
[noise at the ports](portnoise.md) works through warm ports and a readout
chain.
