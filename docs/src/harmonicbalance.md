# Harmonic balance

Harmonic balance solves for the Fourier coefficients of a driven circuit.
A second, linearized calculation gives its response to weak signals and
their idlers. Start with the [quickstart](quickstart.md) for a complete JPA
example, or use the small setup below while reading this guide.

## Running a solve

```@example hbguide
using JosephsonCircuits
circuit = Circuit([
    (:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1e-9)),
    (:cj, 2, 0, Capacitor(1e-12)),
])
ws = 2pi .* [4.6e9, 4.7e9, 4.8e9]
wp = (2pi*4.75001e9,)
sources = [(mode = (1,), port = 1, current = 0.00565e-6)]
sol = hbsolve(ws, wp, sources, (8,), (16,), circuit)
@assert sol.nonlinear.solverinfo.converged
round.(10 .* log10.(abs2.(sol.linearized.S((0,), 1, (0,), 1, :))); digits = 3)
```

Frequencies are in rad/s. `current` is a complex Fourier coefficient in
amperes; a nonzero real coefficient `Ip` represents `2Ip*cos(w*t)`.
A DC coefficient is not doubled. See [conventions](conventions.md).

| Entry point | Result |
|---|---|
| `hbnlsolve(wp, Npumpharmonics, sources, circuit)` | Strong-drive operating point |
| `hblinsolve(ws, circuit; nonlinear, Nmodulationharmonics)` | Small-signal sweep about an operating point |
| `hbstability(circuit; nonlinear, method)` | Local stability, the temporal poles; see [stability](stability.md) |
| `hblinsolve(ws, circuit)` | Linear response, with junctions linearized at zero phase |
| `hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics, circuit)` | Operating point and linearized sweep |

Parameterized circuits take their definitions after the circuit argument.
Sources are named tuples `(mode, port, current)`. With two independent
pumps, use `(1, 0)` and `(0, 1)` for their fundamentals. For commensurate
pumps, use one fundamental in `wp` and drive its harmonics with different
mode indices.

## Reading the results

`hbsolve` returns `sol.nonlinear` and `sol.linearized`.

### The operating point

| `NonlinearHB` field | Meaning |
|---|---|
| `nodeflux` | Reduced node-flux Fourier coefficients, indexed by `outputmode, node` |
| `dcnodevoltage` | Average node voltages in V when a zero-frequency mode exists; otherwise `nothing` |
| `modes` | Retained pump-mode tuples |
| `solverinfo` | Convergence flag, residual history, and stage diagnostics |
| `S` | Output-wave/incident-drive ratios of this strong-drive solution |
| `sources` | The drive, as `(mode, port, current)`; at a nonzero mode of frequency `f` the current in time is `2 real(current*cis(f*t))` |

`NonlinearHB.S` is not the differential scattering response. With several
driven port modes, each column includes the response to all sources.
Use `sol.linearized.S` for small-signal scattering.

At nonzero angular frequency `w`, the physical voltage coefficient is
`im*w*phi0*nodeflux`. The zero-frequency `nodeflux` is static flux; it
is distinct from `dcnodevoltage`.

### The small-signal response

Keyed arrays allow selection by physical port number and mode tuple:

```@example hbguide
Ssignal = sol.linearized.S(outputmode = (0,), outputport = 1,
    inputmode = (0,), inputport = 1, freqindex = :)
Sidler = sol.linearized.S((-2,), 1, (0,), 1, :)
(size(Ssignal), size(Sidler))
```

Mode `(0,)` is the signal at `ws`; mode `(k,)` is at `ws + k*wp[1]`.
In four-wave mixing, `(-2,)` is the conjugate idler coordinate at
`ws - 2wp[1]`. With three-wave mixing enabled, the analogous idler is
`(-1,)`. See the [mode table](conventions.md#Modes-and-signed-frequencies).

| Output | Keyed axes, in order | Interpretation |
|---|---|---|
| `S`, `QE`, `QEideal` | `outputmode, outputport, inputmode, inputport, freqindex` | Scattering and efficiency for each input/output pair |
| `CM`, `nbar` | `outputmode, outputport, freqindex` | Output commutator; occupation of the wave leaving each port mode |
| `Snoise` | `inputmode, component, outputmode, outputport, freqindex` | Response from internal noise channels |
| `Cnoise`, `Vout` | `outputmode, outputport, conjoutputmode, conjoutputport, freqindex` | Added output-noise covariance; total output covariance |
| `nodeflux`, `voltage` | `outputmode, node, inputmode, inputport, freqindex` | Internal response to a unit input current coefficient |
| `Ssensitivity` | `outputmode, outputport, inputmode, inputport, component, freqindex` | Relative component derivative of `S` |

Use `returnSnoise`, `returnCnoise`, `returnVout`, `returnnodeflux`,
`returnvoltage`, or `returnSsensitivity` to request optional outputs;
`nbar` is returned by default, like `QE` and `CM`. Unrequested array outputs
are empty. With `keyedarrays=false`, modes and ports are flattened with
mode varying fastest; see [`LinearizedHB`](@ref JosephsonCircuits.LinearizedHB) for the plain-array layout.

`abs2(S)` is photon gain. Power gain is
`abs(wout)/abs(win) * abs2(S)` for nonzero frequencies. The two coincide
for reflection at the signal frequency. Frequency-converting power-wave
amplitudes require the square root of that absolute frequency ratio.

## The modes and their truncation

The harmonic limits select unknowns on a finite Fourier grid. The
nonlinearity is evaluated on a separate, usually larger grid.

| Keyword | Role |
|---|---|
| `Npumpharmonics` | Largest retained harmonic index along each pump axis |
| `Nevaluationharmonics` | Nonlinear evaluation grid; twice the retained limits by default |
| `Nmodulationharmonics` | Pump-harmonic offsets retained around the signal |
| `dc` | Include the zero-frequency mode |
| `fourwavemixing`, `threewavemixing` | Select the corresponding parity of pump harmonics and signal offsets |
| `maxpumpintermodorder`, `maxmodulationintermodorder` | Limit intermodulation order by the sum of absolute indices |
| `frequencywindow` | Bounds on the absolute frequency of retained pump modes |

The default evaluation padding avoids aliasing of the leading cubic
products into retained modes. It is not an exact representation of a sine
at arbitrary phase excursion. The [theory page](harmonicbalancetheory.md#The-mode-set-and-the-transforms)
explains the grids and aliasing. The [multi-tone tutorial](recipes/multitone.md)
visualizes the retained lattice and works through a three-tone solve.

## Checking convergence

Check three different questions:

1. **Did the nonlinear solve converge on this grid?** Inspect
   `sol.nonlinear.solverinfo.converged` and its residual history.
2. **Is the Fourier truncation adequate?** Increase pump harmonics,
   modulation harmonics, and evaluation harmonics separately. Compare the
   quantities of interest, including weak conversion products.
3. **Is the noise/mode representation consistent?** Check `CM` against
   `+1` at positive output frequency and `-1` at negative frequency.
   This includes internal noise contributions where present.

A small residual or commutator error does not by itself establish
convergence of gain or noise. Near a threshold, small shifts in the
operating point can produce much larger changes in gain.

## The nonlinear solver

| Method | Use |
|---|---|
| `NewtonKrylov()` | Default; exact real Jacobian-vector products with preconditioned GMRES |
| `Newton()` | Assemble and factorize the exact real Jacobian |
| `Staged()` | Increase drive and retained harmonics in stages when a cold solve is difficult |
| `QuasiNewton()` | Approximate holomorphic Jacobian with Anderson acceleration |

A solve that fails to converge returns its last iterate and records the
failure in `solverinfo`. It may warn about an exhausted work budget, a
line search, or stagnation. Do not interpret that iterate as a converged
operating point.

`Staged()` reuses converged states along a source-continuation path. Its
history can indicate where the path became difficult, but failure is not
proof that no operating point exists. Another initial state, a finer
schedule, or a different method can reach another branch.

The following alternatives continue the setup at the start of this page;
choose one rather than running every method:

```julia
sol = hbsolve(ws, wp, sources, (8,), (16,), circuit; method = Staged())
sol = hbsolve(ws, wp, sources, (8,), (16,), circuit;
    method = NewtonKrylov(preconditioner = MeasuredBand()))
```

`atol` controls the scaled nonlinear residual. It is not a bound on gain
error. Method options belong to the method object, for example
`NewtonKrylov(refresh=Probe())` or
`NewtonKrylov(linesearch=Backtracking(interpolate=false))`.
See [performance](performance.md) for preconditioners and reuse, and
[interoperability](interop.md) for external solvers.

## Recovering from a failed solve

This controlled example gives exact Newton only one iteration, so it
returns an unconverged iterate. It then solves the same requested drive
with source continuation. The warning in the first call is expected.

```@example hbguide
failed = hbnlsolve(wp, (8,), sources, circuit;
    method = Newton(), iterations = 1, atol = 1e-10)
@assert !failed.solverinfo.converged
stage = only(failed.solverinfo.stages)
@assert stage.reason == :iterations
(stage.reason, stage.normresidual, stage.alpha, stage.backtracks)
```

Here the residual decreased and the full step was accepted: the small
iteration budget is the immediate problem. Increasing that budget would
also be reasonable. `Staged()` demonstrates recovery by solving a sequence
of easier drive/grid problems instead:

```@example hbguide
recovered = hbnlsolve(wp, (8,), sources, circuit;
    method = Staged(), atol = 1e-10)
@assert recovered.solverinfo.converged
(recovered.solverinfo.finalresidual, length(recovered.solverinfo.stages))
```

The final convergence flag applies to the requested drive and grid, not
merely to an intermediate continuation stage. A finite `sourcefold` means
the continuation reported a branch ending at that fraction of the target
drive. It is not proof that no other branch reaches the target.

For an ordinary Newton or Newton–Krylov solve, each stage is an
[`IterationInfo`](@ref JosephsonCircuits.IterationInfo). A `Staged()` record
contains its inner solver records; inspect those when diagnosing a stage.

| Evidence | What to try next |
|---|---|
| `reason == :iterations`, residual still falling, reasonable `alpha` | Increase the iteration budget; compare the residual decrease, not just the count |
| `reason == :work`, repeated GMRES stagnation, or large `residualratio` relative to `forcing` in `stage.krylov` | Inspect preconditioner refresh/escalation and operator-product counts; on a small circuit compare with `Newton()` |
| `reason == :linesearch`, repeated backtracks or tiny `alpha` | Check units and the initial state; reduce the parameter/drive step or try `Staged()`; a larger GMRES budget alone need not help |
| `reason == :progress` | Inspect the residual history and continuation stages; the solver has already attempted recovery before declaring a stall |
| `converged == true`, but gain changes on harmonic refinement | Refine pump, modulation, and evaluation grids separately; nonlinear iteration tolerance does not control truncation error |
| A sweep jumps despite convergence | Compare state/response continuity and forward/backward sweeps; inspect cache retry records and check [stability](stability.md) |

These observations guide diagnosis; none uniquely identifies a physical
bifurcation. Strong nonlinearity can make a continuation step fail without
a branch ending. Compare cold and warm starts near a suspected jump. The
[cache workflow](performance.md#Reuse-across-a-sweep-of-values) explains
its automatic cold retry and how that can change the branch being followed.

## The linearized sweep

Reuse a converged operating point when changing the signal frequencies or
modulation truncation:

```@example hbguide
lin = hblinsolve(ws, circuit; nonlinear = sol.nonlinear,
    Nmodulationharmonics = (8,))
@assert isapprox(lin.S, sol.linearized.S; rtol = 1e-10)
nothing # hide
```

Each signal frequency requires a linear system over its signal and idler
modes. The solver reuses the sparse pattern across the sweep.
A mode whose signal-plus-pump frequency is numerically zero is rejected:
photon-wave normalization has no DC limit. Study a limiting response with
small nonzero frequencies instead.

## Noise and quantum efficiency

The noise calculation includes the field each port's termination sends
in and the noise emitted by supported internal losses. A termination is at
zero temperature, sending in the vacuum, unless it states one:
`Port(1; termination = MatchedTermination(temperature = 0.05))`. Internal
component temperatures and block noise models follow the
[temperature table](conventions.md#Noise-normalization-and-temperature);
the analysis `temperature` does not warm the ports.

For a selected input/output pair, `QE` is its photon gain divided by twice
the total noise at the output, with every input in its state, the
selected input's own included: a warm source lowers the efficiency of a
measurement of its signal. `QEideal` is the package's ideal-amplifier
reference at that gain; `QE/QEideal` compares the device with that
reference. `nbar` is the occupation of the wave leaving each port mode:
at an output, the photons the measurement receives; at the port facing a
device, what the circuit sends back toward it.

`Cnoise` contains only the internal added covariance, so a solved device
can be embedded as a block with `NoiseCovariance(Cnoise)`. `Vout` is the
whole output covariance, `S*Diagonal(sigma)*S' + Cnoise`, with `sigma` the
`nbar + 1/2` of every input. `Snoise` describes transfer coefficients and
is independent of temperature; with `channeltemperatures` it gives each
internal channel's share of the noise at an output, a noise budget.
[Noise at the ports](portnoise.md) works these through an input line and
a readout chain.

Continuing the circuit above:

```@example hbguide
noisy = hbsolve(ws, wp, sources, (8,), (16,), circuit;
    temperature = 0.05, returnCnoise = true)
@assert noisy.nonlinear.solverinfo.converged # hide
noisy.linearized.QE((0,), 1, (0,), 1, :) ./
    noisy.linearized.QEideal((0,), 1, (0,), 1, :)
```

This circuit has no internal dissipative component, so changing the
analysis temperature alone does not add thermal noise. To model a warm
internal load, add a resistor or a block with a thermal noise model.
Noise outputs for lossy mutually coupled inductors are not supported;
see [component support](circuits.md#Supported-analyses).

## Sensitivities

`sensitivitynames` selects relative component derivatives,
`dS/dr` for `p -> r*p` at `r=1`. By default the derivative includes the
shift of the pump operating point. Use
`sensitivityoperatingpoint=false` only when the operating point is meant
to remain fixed.

```julia
sensitive = hbsolve(ws, wp, sources, (8,), (16,), circuit;
    sensitivitynames = ["jj", "cc"], returnSsensitivity = true)
dS_dlnLj = sensitive.linearized.Ssensitivity((0,), 1, (0,), 1, "jj", :)
```

[`designsensitivities`](@ref) instead differentiates with respect to the
parameters used to define component values. It combines every component's
contribution by the chain rule. The [worked example](recipes/sensitivities.md)
checks a design derivative against finite differences and converts it to
a gain derivative.

## Direct current and flux pumping

Use `dc=true` and a source with `mode=(0,)` for DC bias. A nonzero
`CurrentSource` component also requires a zero-frequency mode. A flux pump
can drive a bias inductor coupled to the device by a mutual inductor.

Static flux determines inductor currents and junction phases. Average node
voltage is a separate coordinate returned as `dcnodevoltage`. Thus a
resistively grounded circuit can carry DC even without an inductive path
to ground. A gauge fixes an undetermined static flux offset; it does not
remove resistive DC conduction.

A current with no supported DC return path is rejected, as is a component
whose required DC conductance is not finite and real. Scattering blocks
use their stated zero-frequency limit. See the [DC example](recipes/dc.md).

## Execution on a device

See [GPU execution](performance.md#GPU-execution) for dependencies,
backend selection, and host/device boundaries.

## Reuse across a sweep of values

See [cache reuse](performance.md#Reuse-across-a-sweep-of-values) for a
complete parameter-sweep example.
