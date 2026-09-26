# Transient simulation

[`transientsolve`](@ref) integrates the circuit equations under prescribed
current waveforms. Use it for pump turn-on, pulses, depletion, and drives
whose many independent tones make a harmonic grid impractical. The tangent
and adjoint calculate incremental responses about the recorded trajectory.

## Example

This one-port junction circuit is driven by a smooth pump and a smaller
signal. The source returns instantaneous current in amperes.

```@example transient
using JosephsonCircuits
circuit = Circuit([
    (:P1, 1, 0, Port(1; Z0 = 50.0)),
    (:C1, 1, 0, Capacitor(1e-12)),
    (:Lj1, 1, 0, JosephsonJunction(1e-9)),
])
rise(t) = t <= 0 ? 0.0 : t >= 2e-9 ? 1.0 : sinpi(t/4e-9)^2
drive(t) = 0.12e-6*rise(t)*cospi(2*3e9*t) +
    2e-9*rise(t - 2e-9)*cospi(2*1.1e9*t)
problem = transientproblem(circuit; sources = [TransientSource(1, drive)])
solution = transientsolve(problem, (0.0, 8e-9); dt = 2e-12)
@assert all(isfinite, solution.voltage) # hide
(size(solution.voltage), last(solution.times))
```

Rows follow `problem.circuit.ports`; columns are saved times in seconds.
The main outputs are:

| Field | Meaning |
|---|---|
| `voltage` | Port voltage, V |
| `incident`, `outgoing` | Instantaneous real power waves, `sqrt(W)` |
| `times` | Saved sample times, s |
| `stats` | Solver work and convergence statistics |
| `finalflux`, `finalrate` | Scaled final state arrays; use `transientstate(solution)` to continue |

For a matched termination `R`, a sinusoidal Norton current of peak
amplitude `Ipeak` launches available power `Ipeak^2*R/8`. The same pump's
HB coefficient is half its peak amplitude; see [conventions](conventions.md).

`TransientSource("I1", waveform)` replaces the constant value of a named
`CurrentSource`. Multiple waveforms on the same target sum. Source
functions must return finite real values and be deterministic: response
calculations and checkpoint replay evaluate them again.

## Choosing the rule, the step and the record

### Stepping rule

| Method | Nominal differential order | Practical use |
|---|---:|---|
| `GaussLegendre()` | 4 | Default; low phase error, blocks, lines, and batched conditions |
| `Trapezoidal()` | 2 | Nondamping reference; supports an optional iterative linear solve |
| `BackwardEuler()` | 1 | Strongly damping comparison method |

Gauss–Legendre and trapezoidal integration do not damp a resolved linear
lossless LC oscillation. This does not make large steps accurate: phase
error shifts resonances and can strongly change amplifier gain near a
threshold. Neither method damps unresolved fast modes as an L-stable
method would. Algebraic variables can have lower accuracy than the
nominal differential order; see [theory](transienttheory.md).

### Time step and output sampling

`dt` is the maximum step. A uniform grid is chosen and shortened slightly
if needed to land on the final time. There is no adaptive temporal-error
control. `rtol` and `atol` control the nonlinear residual at each step,
not the discretization error.

Repeat with `dt/2` and compare the quantities you need: amplitudes, phases,
weak products, and pulse energy. Resolve generated harmonics and circuit
resonances, not only the applied drives. Align discontinuities with the
grid where possible; a smooth turn-on generally requires fewer high
frequencies to resolve.

`saveevery` decimates saved samples without changing the integration grid
and always preserves the endpoints. It applies no anti-alias filter.
The saved rate must still resolve anything you demodulate later.

### Recording level

| `record` | Use | Stored trajectory |
|---|---|---|
| `:ports` | Port waveforms and continuation | Default; port traces and endpoint states |
| `:phases` | Tangent, adjoint, and noise about a nonlinear trajectory | Junction phases at every integration step |
| `:states` | State inspection and state-dependent component derivatives | Full flux/rate history as well as phases |
| `:checkpoints` | Response calculations on long records | States at checkpoints, with intermediate steps replayed |

The last three require `saveevery=1`. Checkpoints are supported under
`GaussLegendre()`; `checkpointevery` chooses their spacing. They reduce
stored state history, while port traces and some input/output arrays still
grow with the record length. See [memory and replay](performance.md#Transient-records-and-reuse).

## Supported circuits and the initial state

The [support table](circuits.md#Supported-analyses) distinguishes component
representations. Real constant lumped elements work directly. Constant
real scattering blocks, rational blocks, fitted pumped blocks, and ideal
transmission lines require `GaussLegendre()`. General complex or
frequency-dependent component values need a causal time-domain model.

The default initial state is zero. Use [`transientstate`](@ref) to supply
node fluxes in webers and node voltages in volts, in compiled node order.
It constructs the scaled state and the associated auxiliary coordinates.
The solver checks the initial algebraic constraints; it does not search
for an operating point or project an inconsistent initial state.

For example, a resistor driven by nonzero current at the first sample
needs its initial voltage supplied. A shunt capacitor can instead begin
charging from zero. The [amplitude example](conventions.md#Check-the-convention-on-a-resistor)
shows an explicitly initialized resistor.

To continue a solve, keep the full state:

```@example transient
continued = transientsolve(problem, (8e-9, 9e-9); dt = 2e-12,
    initialstate = transientstate(solution))
@assert all(isfinite, continued.voltage) # hide
nothing # hide
```

With transmission lines this includes the required prehistory; with
rational blocks it includes the internal states. At a different time step,
line history is read through the line interpolation. There is currently
no adapter from an HB operating point to a transient initial state.

## Demodulate a port trace

[`transientdemodulate`](@ref) returns a complex peak amplitude at a
frequency in Hz. Continuing the example above:

```@example transient
window(t) = 4e-9 <= t <= 8e-9 ? sinpi((t - 4e-9)/4e-9)^2 : 0.0
amplitude = transientdemodulate(solution, 1, 1.1e9;
    quantity = :outgoing, window)
@assert isfinite(amplitude) # hide
nothing # hide
```

A smooth window reduces leakage from a strong pump. It does not separate
arbitrarily close tones: their separation also sets the required
observation time. The result above includes the response to both applied
drives; use a tangent calculation to isolate an incremental probe.

## Tangent and adjoint of the recorded steps

The tangent propagates a small change to the drive about the entire
trajectory, including the pump, signal, and depletion. The adjoint gives
the derivative of a selected output objective with respect to the drives.
This example checks their transpose identity:

```@example transient
recorded = transientsolve(problem, (0.0, 2e-9); dt = 2e-12, record = :phases)
deltacurrent = zeros(size(recorded.voltage))
deltacurrent[1, :] .= 1e-9 .* sinpi.(2*1.8e9 .* recorded.times)
response = transienttangent(recorded, deltacurrent)

weights = zeros(size(recorded.outgoing))
weights[1, end] = 1.0   # final outgoing wave at port 1
adj = transientadjoint(recorded, weights; quantity = :outgoing)
forward = sum(weights .* response.outgoing)
reverse = sum(adj.currents .* deltacurrent)
@assert isapprox(forward, reverse; rtol = 1e-8, atol = 1e-15)
(forward, reverse)
```

For an integral objective, include quadrature weights in `weights`.
Use separate real and imaginary objectives for a complex demodulation.
`adj.initialflux` and `adj.initialrate` describe sensitivity to the initial
scaled state; blocks and lines have additional state contributions.

Port numbers or component names select drive targets. A trailing array
dimension carries several tangent directions or adjoint objectives through
the same step factorization. The derivatives are exact for the discretized
equations to the response-solve tolerances; refine `dt` to check physical
accuracy.

### Component sensitivities

[`transientsensitivity`](@ref) differentiates port responses with respect
to relative changes of named `C`, `L`, `R`, and `Lj` values. For `p -> r*p`,
it returns the derivative at `r=1`. `transientadjoint(...; components=...)`
contracts the same derivatives against the objective.

```@example transient
states = transientsolve(problem, (0.0, 2e-9); dt = 2e-12, record = :states)
sensitivity = transientsensitivity(states, ["Lj1", "C1"])
adj_components = transientadjoint(states, weights; quantity = :outgoing,
    components = ["Lj1", "C1"])
@assert isapprox(sum(weights .* sensitivity.outgoing[:, :, 1]),
    adj_components.sensitivity[1]; rtol = 1e-8, atol = 1e-15)
nothing # hide
```

Capacitor, inductor, and resistor derivatives require states or checkpoint
replay; junction derivatives can use the phase record. Changing a port
termination also changes its wave normalization, which is included in the
derivative.

## Many drive conditions as one solve

Rebind the waveforms of a compiled problem to simulate several drive
conditions under `GaussLegendre()`. The topology, driven targets, and
sources left constant must be the same; only the waveforms differ.

```@example transient
base = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
pump(Ipeak) = t -> Ipeak*rise(t)*cospi(2*3e9*t)
problems = [transientproblem(base; sources = [TransientSource(1, pump(a))])
    for a in (0.04e-6, 0.08e-6, 0.12e-6)]
batch = transientsolve(problems, (0.0, 2e-9); dt = 5e-12, record = :phases)
member = batch[2]
size(batch.voltage)   # port, time, condition
```

A member is an ordinary solution for demodulation, responses, or noise.
Batched responses put the condition dimension last. Each condition retains
its own stiffness, factorization, convergence check, and initial state.
See [performance](performance.md) for CPU/GPU tradeoffs and workspace reuse.

## Scattering blocks, transmission lines and fitted data

See [scattering blocks](scattering.md) for complete fitting examples,
reference-plane delays, and pumped-block restrictions. In an initial state,
`linecurrents` sets line DC currents; `waves` and `blockstates` hold the
line and rational-block state. A noise calculation additionally requires
an equilibrium prehistory.

## GPU execution

See [GPU execution](performance.md#GPU-execution). The solution arrays stay
on the backend. `transientdemodulate` downloads the selected port trace;
`Array(solution.outgoing)` explicitly downloads all outgoing traces.

## A Josephson transmission line with ten pulsed signals

The [ten-tone line recipe](recipes/transient-line.md) constructs a small
junction line and measures transmitted tones, an idler, and a third
harmonic in a smooth time window.

## A JPA against WRspice at the signal frequency

The [JPA comparison](recipes/transient-wrspice.md) evaluates the
small-signal reflection by tangent propagation and compares it with HB and
WRspice, both after settling and during pump turn-on.
