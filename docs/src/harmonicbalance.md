# Harmonic balance

The frequency domain solvers find the periodic steady state of a circuit
driven by one or more strong tones, the pumps, and then the response of
weak signals through the circuit linearized about that state, with its
noise, quantum efficiency and sensitivities. This page is the usage
guide; the [theory and implementation](harmonicbalancetheory.md) page
explains what the solvers do and why, and the [examples](examples.md)
run them on amplifiers from a single junction to a traveling wave line.

## Running a simulation

Three functions run the analyses. Add a question mark `?` in front of a
function to read its docstring, for example `?hbsolve`.

- `hbnlsolve(wp, Npumpharmonics, sources, circuit)` solves the nonlinear
  circuit driven by the pumps: the harmonic balance solution of the pump
  and its harmonics at one operating point.
- `hblinsolve(ws, circuit; nonlinear = ...)` sweeps weak signals through
  the circuit linearized about that operating point, or through a linear
  circuit when no operating point is given.
- `hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics, circuit)`
  runs both in sequence, which is what the examples below do.

The drive is a vector of sources, each a mode, a port and a current
amplitude: `(mode = (1,), port = 1, current = Ip)` is a current of amplitude
`Ip` at the first pump frequency `wp[1]` applied to port 1, and
`(mode = (0,), port = 2, current = Idc)` a direct current bias on port 2
(with `dc = true`). With two pumps the mode is a pair, `(1, 0)` and
`(0, 1)`. The harmonic counts say how many pump harmonics the nonlinear
solve keeps and how many signal and idler modes the linearized sweep
keeps, one count per pump.

## Reading the results

`hbnlsolve` returns a `NonlinearHB`, the operating point. Its main fields:

- `nodeflux`: the node flux at every retained mode of every node, a keyed
  array with axes `outputmode` and `node` (a plain matrix with
  `keyedarrays = false`). The zero mode is the static flux; the voltage of
  a mode at frequency `w` is `i*w*phi0` times its flux.
- `S`: the scattering parameters at the pump frequencies, which measure
  how much of each pump is reflected and converted.
- `dcnodevoltage`: the average voltage of each node in volts, when the
  analysis has a zero frequency mode. See the direct current example.
- `solverinfo`: whether the solve converged (`solverinfo.converged`) and
  the residual history and step record of every solver stage, for
  diagnosing a solve which did not.
- `modes`: the retained modes as tuples of harmonic indices, which label
  the mode axes of the arrays above.

`hblinsolve` returns a `LinearizedHB`, the small signal response. Its
frequency dependent outputs are indexed by output mode, output port, input
mode, input port and signal frequency, and are keyed arrays: the gain of
the JPA below is read as
`S(outputmode = (0,), outputport = 1, inputmode = (0,), inputport = 1, freqindex = :)`,
or positionally as `S((0,), 1, (0,), 1, :)`. Mode `(0,)` is the signal
itself and `(k,)` the idler offset by `k` pump harmonics, so
`S((-1,), 1, (0,), 1, :)` is the conversion from the signal to the first
idler. The fields:

- `w`: the signal frequencies of the sweep, and `modes` the retained
  signal and idler modes.
- `S`: the scattering parameters in units of photon flux, so that gain in
  dB is `10*log10.(abs2.(S))` and a conversion between frequencies is
  read in photons; multiply by `sqrt(w_out/w_in)` for power.
- `QE` and `QEideal`: the quantum efficiency of each output, and that of
  an ideal amplifier with the same gain, so `QE ./ QEideal` is the
  fraction of the ideal.
- `CM`: the commutation relation of each output, which is `1` when the
  scattering matrix is complete; its deviation from `1` measures modes the
  truncation left out.
- `Snoise` and `Cnoise` (on request): the scattering from the noise
  channels of the dissipative elements to the ports, and the added noise
  covariance at a given temperature.
- `nodeflux` and `voltage` (on request, `returnnodeflux = true` and
  `returnvoltage = true`): the node fluxes and voltages resulting from a
  unit input at each port and mode, for looking inside the circuit.
- `Ssensitivity` (on request): the derivative of `S` with respect to the
  named components or design parameters.

`hbsolve` returns an `HB` holding both, as `nonlinear` and `linearized`;
the examples below read the gain from `sol.linearized.S` and the operating
point from `sol.nonlinear`. Every output which was not requested is an
empty array, and the `return...` keywords of `hbsolve` and `hblinsolve`
choose which are computed.


## The modes and their truncation

A solve retains a finite set of modes, the harmonics and
intermodulation products of the pumps, and the keywords choose which.
Every mode is a tuple of harmonic indices, one per pump, and its
frequency is the dot product with the pump frequencies.

| Keyword | What it selects |
| --- | --- |
| `Npumpharmonics` | the largest harmonic index of each pump kept as an unknown of the nonlinear solve |
| `Nevaluationharmonics` | the grid the nonlinearity is sampled on, twice the retained set by default, which dealiases the cubic products of a junction |
| `Nmodulationharmonics` | the harmonics of each pump kept around the signal in the linearized solve, which sets the signal and idler modes |
| `dc` | the zero frequency mode, for a direct current bias |
| `fourwavemixing`, `threewavemixing` | the odd and the even pump harmonics, the ones four and three wave mixing couple through; in the linearized solve, the even and the odd offsets from the signal |
| `maxpumpintermodorder`, `maxmodulationintermodorder` | a diamond truncation of the multi pump lattice by the absolute sum of the indices |
| `frequencywindow` | a lower and upper bound on the absolute frequency of the retained pump modes |

Convergence in the harmonics is checked by raising `Npumpharmonics` and
`Nmodulationharmonics` until the outputs stop moving; the commutation
relation `CM` of each output measures the modes the truncation left out,
and is `1` when the scattering matrix is complete. Two pumps should be
incommensurate; a commensurate pair is written as one frequency with the
other as a source at the mode index of the ratio, and a product that
lands on zero frequency is refused.

## The nonlinear solver

`method` chooses how the operating point is solved. The default
[`NewtonKrylov`](@ref) is matrix free: Newton steps whose linear
systems GMRES solves with the exact Jacobian-vector product through the
transforms and a preconditioner, [`Automatic`](@ref) by default, which
picks the mode block diagonal for one pump and the full Jacobian in
single precision block factors for two or more when they fit in half the
free memory. A preconditioner that stalls is grown, so the method is
never less robust than a direct solve. [`Newton`](@ref) assembles the
exact real Jacobian and factorizes it, [`QuasiNewton`](@ref) uses the
holomorphic approximation with Anderson acceleration, and
[`Staged`](@ref) is source continuation on a ladder of harmonic grids,
the method for an operating point the others cannot reach from a cold
start, and the one that tells a hard operating point from one that does
not exist.

```julia
sol = hbsolve(ws, wp, sources, (2,), (8,), circuit; method = NewtonKrylov(preconditioner = MeasuredBand()))
sol = hbsolve(ws, wp, sources, (2,), (8,), circuit; method = Staged())
sol.nonlinear.solverinfo.converged
```

`ftol` is the residual tolerance of the scaled system, independent of
the units, and is raised to the rounding floor of the source when that
is larger. A solve that does not converge returns its last iterate with
`solverinfo.converged = false` and warns with the reason it stopped: the
iterations spent, the work budget spent, a line search without decrease,
or a residual history that projects no convergence. Check the flag
before using a result. `NewtonKrylov(precision = Float32)` iterates in
single precision, and `refresh = Probe()` rebuilds the preconditioner
only when a measurement says it pays.

## The linearized sweep

[`hblinsolve`](@ref) sweeps the signal frequencies through the circuit
linearized about the operating point, or through a linear circuit when
none is given. Each signal frequency is one sparse linear system in the
signal and idler modes, factorized and solved for a unit current at
every port and mode; the frequencies are split into `nbatches` batches
over the threads, and on a device a batch is assembled and solved as
one uniform batch. `factorization` is the sparse factorization of the
system, [`KLUfactorization`](@ref) on the host for one pump and the
dense node block [`BlockFactorization`](@ref) for two or more when its
factors fit in memory, in single precision refined to double.

```julia
linear = hblinsolve(2pi*(1:0.01:10)*1e9, circuit)                 # a linear circuit
lin = hblinsolve(ws, circuit; nonlinear = sol.nonlinear, Nmodulationharmonics = (2,))
lin.S((0,), 2, (0,), 1, :)        # signal transmission from port 1 to port 2
lin.S((-1,), 2, (0,), 1, :)       # conversion to the first idler
```

## Noise and quantum efficiency

Every dissipative element is a noise channel: the port terminations,
the resistors, the lossy capacitors and inductors, and the lossy
scattering blocks, whose emitted noise wave has the covariance
`I - S S'` of Bosma's relation. The linearized solve propagates every
channel to the ports by one transposed solve per port, and reports the
quantum efficiency `QE` of each output, the ratio of the signal it
carries to everything it carries, that of an ideal amplifier of the same
gain `QEideal`, and the commutation relation `CM`. `temperature` sets
the occupation of every channel that does not state its own, a resistor
by its `temperature` keyword and a block by its noise model, and
`returnCnoise = true` returns the added noise covariance at the ports;
`Snoise` is the scattering from the channels to the ports and does not
depend on temperature.

```julia
sol = hbsolve(ws, wp, sources, (2,), (8,), circuit; temperature = 0.05, returnCnoise = true)
sol.linearized.QE((0,), 1, (0,), 1, :) ./ sol.linearized.QEideal((0,), 1, (0,), 1, :)
```

## Sensitivities

`sensitivitynames` names the components whose relative perturbation the
scattering parameters are differentiated with respect to, by the
adjoint method, at a fixed operating point or, with
`sensitivityoperatingpoint = true`, including the shift of the pump
operating point through the exact real Jacobian; near the gain peak of a
strongly pumped amplifier the shift is the larger term.
[`designsensitivities`](@ref) differentiates with respect to the
parameters of a circuit builder by the chain rule through the
components, with the exact direction of every dependent value, which is
what a gradient based optimizer wants.

```julia
sol = hbsolve(ws, wp, sources, (2,), (8,), circuit; sensitivitynames = ["Lj", "Cc"], returnSsensitivity = true)
sol.linearized.Ssensitivity
out, dSdp = designsensitivities(make, (Lj = 1e-9, Cc = 100e-15), ws, wp, sources, (2,), (8,))
```

## Direct current and flux pumping

`dc = true` retains the zero frequency mode, whose flux is the static
flux setting the inductor currents and junction phases; a direct
current bias is a source with `mode = (0,)`, and a flux pump a current
source through a mutual inductor. The average node voltages, which the
periodic state alone does not carry, are solved beside it and returned
as `dcnodevoltage`; a resistor is open at direct current in the periodic
state and carries its direct current there. A subnetwork no inductor or
junction connects to ground has a free static flux, fixed by a gauge
row, and a direct current injected into such a subnetwork is refused,
since no periodic solution exists.

## Execution on a device

```julia
using JosephsonCircuits, CUDA, CUDSS
sol = hbsolve(ws, wp, sources, (2,), (8,), circuit; backend = CUDABackend())
```

The nonlinear solve assembles, factorizes and iterates on the device,
its transforms through the device FFT, its preconditioner through cuDSS
or the batched block factorization; the linearized sweep assembles the
system matrices of a batch of frequencies with one kernel and factorizes
and solves them as a uniform batch, falling back to the host for what it
cannot serve, a symbolic frequency variable in the component values.

## Reuse across a sweep of values

A sweep over component values builds the structure once and moves
values: [`hbcache`](@ref) holds what a solve built and
[`hbsolve!`](@ref) reuses it, so a parameter sweep pays the compile, the
symbolic analysis and the plans once. The [other solvers](interop.md)
page shows the problem object every external solver can drive.
