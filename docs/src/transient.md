# The circuit in time

Harmonic balance solves for the amplitudes of a set of tones. When the
drive is a pulse, or carries more tones than a harmonic grid can hold, the
alternative is to integrate the circuit directly in physical time. The
transient solver does that on the same compiled circuit and the same
unknowns as harmonic balance: the node fluxes, with the auxiliary branch
currents of the mutually coupled inductors and the gauge rows of the
floating subnetworks, at one mode. Every source acts on one state, so
pump harmonics, idlers, depletion and intermodulation products need no
additional unknowns, and the state size is that of the circuit and not of
the number of tones. It is a deterministic solver with the exact tangent
and adjoint of its time steps; it does not by itself compute added
thermal or quantum noise.


This page is the usage guide: how to set a circuit up in time, what a
solve returns, how to run many drive conditions at once, how to take the
tangent and the adjoint of a recorded solve, and how to run on a GPU,
with worked examples. The [theory and implementation](transienttheory.md)
page explains what the solver does and why, and the
[quantum noise](transientnoise.md) page how fluctuations are propagated
through a recorded trajectory.

## A first pulse

```julia
using JosephsonCircuits

circuit = [
    ("P1", "1", "0", 1),
    ("R1", "1", "0", 50.0),
    ("C1", "1", "0", 1e-12),
    ("Lj1", "1", "0", 1e-9),
]
rise(t) = t <= 0 ? 0.0 : t >= 15e-9 ? 1.0 : sinpi(t/30e-9)^2
pulse(t) = rise(t - 20e-9)*rise(120e-9 - t)
drive(t) = 0.12e-6*rise(t)*cospi(2*3e9*t) +
    pulse(t)*(2e-9*cospi(2*1.1e9*t) + 2e-9*cospi(2*1.8e9*t))

problem = transientproblem(circuit; sources = [TransientSource(1, drive)])
solution = transientsolve(problem, (0.0, 140e-9); dt = 1e-12)

# port rows follow problem.circuit.ports; time is in seconds
solution.voltage                       # volts
solution.incident, solution.outgoing   # instantaneous waves in sqrt(W)
solution.stats
```

`TransientSource(1, drive)` is an instantaneous Norton current in Amperes,
positive into the positive terminal of port 1, and it adds no termination:
the port's own termination is already part of the compiled circuit. For a
matched port of resistance `R` a sinusoidal peak current `Ip` launches the
available power `Ip^2*R/8`. Unlike a harmonic balance source, the callable
returns the physical waveform, not a Fourier coefficient.

`TransientSource("I1", waveform)` instead replaces the constant value of a
named `CurrentSource`; that component's current flows out of its first
terminal and into its second. Several waveforms on one target add. A source
callable must return finite real values and must be deterministic, because
the tangent and adjoint evaluate it again on the recorded grid.

The circuit may be a typed [`Circuit`](@ref), a compiled circuit, or a
legacy netlist, exactly as for [`hbsolve`](@ref).


## Choosing the rule, the step and the record

[`Trapezoidal`](@ref), the default, applies the trapezoidal rule to the
flux and to its rate, which is Newmark's rule with the averaging
parameters: second order, and free of numerical damping, so a lossless LC
oscillation keeps its energy. [`GaussLegendre`](@ref) is the two stage
Gauss-Legendre collocation, fourth order, A-stable and symplectic, with
the same freedom from damping; on a resonator it is the rule to use, see
below. [`BackwardEuler`](@ref) is first order and strongly damping, a
reference for checking that a result does not depend on the rule. A
trapezoidal step is one implicit equation in the new flux,

The choice of rule is a matter of samples per period. The trapezoidal
rule warps every frequency by `(2 pi f dt)^2/12`, a fifth of a percent
at 42 samples of a 4.75 GHz period, and a resonator's detuning and an
amplifier's bifurcation magnify that: the pumped amplifier of the noise
example has its pump response 123% off harmonic balance at 5 ps and
needs 0.3 ps for half a percent. The Gauss-Legendre rule warps by
`(2 pi f dt)^4/720`: the same amplifier is 5e-4 off at 5 ps and 3e-5 at
2.5 ps, converging as the fourth power, at six microseconds a step on
the CPU where the trapezoidal rule at 5 ps refactorizes every step and
takes twenty five. Neither rule is L-stable: an unresolved fast mode is
not damped, and a sharp edge in a drive is resolved by the grid, not
smoothed by the rule.

`dt` is a maximum step; the uniform grid is shortened slightly to land on
the final time. `rtol` and `atol` control the Newton residual, not the
temporal error, and there is no step control: repeat at `dt/2` and compare
the amplitudes, phases, weak products and pulse energy. Resolve the
highest generated harmonic and the circuit's resonances, not only the
drive. A discontinuous drive needs a grid aligned with it.

By default only the port waveforms and the final state are kept.
`saveevery` decimates the saved samples without changing the integration
grid, always keeping both endpoints; it has no antialiasing filter, so the
saved rate must still resolve whatever is demodulated afterwards.

[`transientdemodulate`](@ref)`(solution, port, frequency; quantity,
window)` integrates a windowed waveform and returns a complex peak
amplitude. A smooth window suppresses leakage from a strong pump; close
products still need enough observation time.


With `linearsolver = GMRES()` each trapezoidal Newton correction is
instead solved matrix free by the package's Krylov solver, with the last factorization
as the preconditioner and a refresh of it only when the iteration count
says the phases have moved away from it, so a whole solve can run on one
factorization. On a one dimensional line on the CPU the direct step is
faster, since KLU refactorizes such a pattern in a tenth of the step;
the iterative step is for patterns with fill, where a refactorization
grows faster than a product does.

A [`TransientReuse`](@ref) passed as `reuse` to a solve, a tangent or an
adjoint carries the scaled system with its Jacobian plan, the
factorization as it was left and the Krylov workspace to the next call
on the same problem at the same step, rule and backend, so a parameter
sweep over the drives, or the tangent and adjoint after a solve, pay the
setup once; it is the counterpart of the reuse between the solves of an
[`hbcache`](@ref).

## Supported circuits and the initial state

Supported are real, constant resistors, capacitors, inductors, mutual
inductors, sinusoidal Josephson junctions, current sources, ports,
[`ScatteringParameters`](@ref) blocks with a constant real matrix,
[`RationalScattering`](@ref) blocks, and ideal
[`TransmissionLine`](@ref)s, the last three under
[`GaussLegendre`](@ref). Frequency dependent or complex values and
other blocks are rejected, since they need a causal realization in
time. An infinite resistance is an open.


[`transientstate`](@ref) builds the initial state, the pair of the scaled
fluxes and their rates, from node fluxes in Weber and node voltages in
Volts in the compiled order; the auxiliary currents follow from the
constitutive equations and the gauge is normalized as harmonic balance
normalizes an initial guess. The default is the zero state. The solver
checks that the state satisfies the algebraic rows at the start, the
nodes without capacitance and the augmentation, and otherwise throws: it
does not look for an operating point or project the state. A resistor
driven by a nonzero current at the first sample needs its initial voltage
supplied; a shunt capacitor can start charging from zero. An adapter
from a harmonic balance operating point to a transient initial state is
future work.


## Scattering blocks, transmission lines and fitted data

A circuit in time may contain constant real scattering blocks, an
attenuator, a circulator, an ideal through, short or open, rational
scattering blocks, which is the form a fit of measured or simulated
scattering data takes, and ideal transmission lines, all under
[`GaussLegendre`](@ref). They are the same components the harmonic
balance solvers take, with the same reference impedances, noise models
and temperatures, so a result can be compared between the two domains
with the same circuit.

```julia
using JosephsonCircuits

# a cable of 60 ohms and 90 ps in front of a two port whose measured
# scattering parameters were fitted at four poles, at 0.3 kelvin
data = ScatteringParameters((2pi .* frequencies, S); nports = 2, zref = 50.0,
    noise = ThermalEquilibrium(0.3))
fitted = RationalScattering(data, 4)
circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:cable, 1, 2, TransmissionLine(60.0, 0.09)),
    (:block, 2, 3, fitted),
    (:cc, 3, 4, Capacitor(100e-15)),
    (:jj, 4, 0, JosephsonJunction(1000e-12)),
    (:cj, 4, 0, Capacitor(1000e-15))])
problem = transientproblem(circuit; sources = [TransientSource(1, pump)])
solution = transientsolve(problem, (0.0, 200e-9); dt = 2.5e-12, method = GaussLegendre())
```

A [`TransmissionLine`](@ref) is given its impedance and length, with a
phase velocity that defaults to the speed of light; its delay must be at
least a step, and at least four steps for the read of its history to
keep the rule's order. A [`RationalScattering`](@ref) block is given as
a real state space realization, validated as stable and passive, or is
fitted from any [`ScatteringParameters`](@ref) block, tabulated,
Touchstone or callable, by `RationalScattering(block, npoles)`; ask for
as many poles as the data might need, since the poles it does not need
are dropped. A delay is not a rational function, so a cable is a line in
cascade with a fit of the data with that delay removed. A block's
`ThermalEquilibrium(T)` sets the temperature of the noise its loss
emits, and a declared `Lossless()` is validated. A
[`transientstate`](@ref) of a circuit with lines takes their direct
currents as `linecurrents`, and holds the waves on the lines and the
states of the blocks as its third and fourth members.

## Many drive conditions as one solve

The typical use of the transient is one circuit under many pumps or
signals. [`transientproblem`](@ref)`(problem; sources)` rebinds the
sources of a compiled problem without compiling again, and
[`transientsolve`](@ref) on a vector of such problems steps them as one
system under [`GaussLegendre`](@ref): the states are matrices with the
conditions as columns, every product and residual takes all conditions
in one call, the junction stiffness and the complex factorization are
one per condition, KLU on the CPU and the uniform cuDSS batch on a
device, and the Newton engine accepts and refreshes per condition. The
problems of a batch drive the same targets and leave the same sources
constant, so that only their waveforms differ, and each condition's
initial state is checked under its own drive.

```julia
base = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
pump(ip) = t -> 2ip*ramp(t)*cospi(2fp*t)
problems = [transientproblem(base; sources = [TransientSource(1, pump(ip))]) for ip in amplitudes]
batch = transientsolve(problems, (0.0, 300e-9); dt = 5e-12, record = :phases)
batch.voltage            # port by time by condition
member = batch[3]        # an ordinary solution: demodulate it, take its adjoint, its noise
```

On a device this is where the throughput lies: a step's launches serve
every condition, so the time per step per condition falls with the
batch until the device is busy. A single problem is a batch of one, on
the same path. The tangent, the adjoint and the noise of a batch run
on the same principle, every condition's directions or objectives in
one pass on one factorization per condition, and return the conditions
as the trailing dimension.

By default a solve keeps the port waves and the first and last states,
which is what the port responses need; `record = :phases` adds the
junction phases at every step, the least the responses and the noise
read, `record = :states` the whole flux and rate history, and
`record = :checkpoints` only the state every `checkpointevery` steps,
from which the responses replay each window of steps as they walk it,
so the phases of a record of any length cost the checkpoints, the
square root of the steps by default, and one window, for one more solve;
what remains per time is the port waves of the record, the weights of an
objective and the currents of a direction, which are per port rather
than per state.


## Tangent and adjoint of the recorded steps

An adjoint given a `sink`, a function of the recorded index and the
final column of the currents, stores no currents at all: each column is
handed over once no remaining step touches it, in decreasing time, which
is how the noise contracts a long record without memory per time.

Record the complete state only when it is needed:

```julia
solution = transientsolve(problem, (0.0, 2e-9); dt = 1e-12, record = :phases)
deltacurrent = zeros(size(solution.voltage))  # A, port rows and time columns
deltacurrent[1, :] .= 1e-9*sinpi.(2*1.8e9*solution.times)
response = transienttangent(solution, deltacurrent)

weights = zeros(size(solution.outgoing))
weights[1, end] = 1.0  # the objective: the final outgoing wave at port 1
adjoint = transientadjoint(solution, weights; quantity = :outgoing)
# sum(weights .* response.outgoing) == sum(adjoint.currents .* deltacurrent)
```

The tangent is linearized about the full recorded trajectory, pump and
signals together, so the loaded junction phases enter every response. The
adjoint is the exact transpose of the steps taken, on the factorization
of the step itself since the step matrix is symmetric, including both
source endpoints of the trapezoidal rule and the direct feedthrough of a
port current into the measured wave; `adjoint.initialflux` and
`adjoint.initialrate` are the sensitivities to what the circuit stored at
the start. For a time integral put the quadrature weights into `weights`,
and use separate real and imaginary objectives for a complex
demodulation. Recording the states costs the state size times the number
of steps.

The tangent and the adjoint take any targets, port numbers or component
names, and a trailing dimension of directions or objectives propagated
together on each step's factorization; the [quantum noise](transientnoise.md)
is built on them, with the baths as targets.

## GPU execution

Every transient solve accepts the package's KernelAbstractions backend
convention:

```julia
using JosephsonCircuits, CUDA, CUDSS
CUDA.allowscalar(false)
solution = transientsolve(problem, (0.0, 140e-9); dt = 1e-12, backend = CUDABackend())
outgoing = Array(solution.outgoing)   # download once for host analysis
```

The default is `backend = CPU()`. The sparse factorization is KLU on the
CPU and the package's [`CUDSSFactorization`](@ref) on an NVIDIA device,
through the same `factorize` and `refactorize!` interface as the harmonic
balance solvers, so another device needs a factorization for it. On a
device the state, the products, the junction term, the Jacobian assembly,
the factorization and the solves all run there, through the package's
device sparse matrices and its assembly kernels; the compiled circuit, the
sample times and the source callables stay on the host, and one vector of
instantaneous drive currents is transferred each step. The Newton control
needs scalar norms, so the loop synchronizes between kernels, and a small
circuit is not expected to run faster on a GPU.

`solution.voltage`, `incident`, `outgoing`, `finalflux`, `finalrate` and
the optional `flux` and `rate` stay on the backend.
[`transientdemodulate`](@ref) downloads only the requested port trace.
[`transienttangent`](@ref) and [`transientadjoint`](@ref) run on the
solution's backend and return arrays there.


## A Josephson transmission line with ten pulsed signals

The line below is a small demonstration circuit, not a calibrated
travelling wave amplifier. A pump at 7.5 GHz turns on smoothly and stays
on while a 100 ns pulse of ten signals between 4 and 6.4 GHz passes. The
solve keeps only the port waveforms and the final state, and the outgoing
signals, the pump's third harmonic and an intermodulation product are read
by demodulating the port 2 wave through a smooth window.

```julia
using JosephsonCircuits

function transientline(cells)
    circuit = Tuple{String,String,String,Float64}[]
    push!(circuit, ("P1", "1", "0", 1.0), ("R1", "1", "0", 50.0))
    for k in 1:cells
        push!(circuit, ("Lj$k", string(k), string(k + 1), 100e-12))
        push!(circuit, ("Cj$k", string(k), string(k + 1), 20e-15))
        push!(circuit, ("Cg$k", string(k), "0", 40e-15))
    end
    last = string(cells + 1)
    push!(circuit, ("Cend", last, "0", 40e-15), ("P2", last, "0", 2.0),
        ("R2", last, "0", 50.0))
    return circuit
end

rise(t, width) = t <= 0 ? 0.0 : t >= width ? 1.0 : sinpi(t/(2width))^2
pulse(t) = rise(t - 20e-9, 15e-9)*rise(120e-9 - t, 15e-9)

cells, ntones = 64, 10
frequencies = collect(range(4e9, 6.4e9; length = ntones))
phases = [pi*j*(j - 1)/ntones for j in 1:ntones]
fp, Ip, Is = 7.5e9, 2e-6, 5e-9
function drive(t)
    pump = Ip*rise(t, 15e-9)*cospi(2fp*t)
    signals = sum(cospi(2frequencies[j]*t + phases[j]/pi) for j in 1:ntones)
    return pump + Is*pulse(t)*signals
end

problem = transientproblem(transientline(cells); sources = [TransientSource(1, drive)])
solution = transientsolve(problem, (0.0, 140e-9); dt = 2e-12)

window(t) = 40e-9 <= t <= 100e-9 ? sinpi((t - 40e-9)/60e-9)^2 : 0.0
for f in frequencies
    a = transientdemodulate(solution, 2, f; window)
    println("$(f/1e9) GHz: $(abs(a)) sqrt(W) at $(angle(a)) rad")
end
idler = transientdemodulate(solution, 2, 2fp - first(frequencies); window)
third = transientdemodulate(solution, 2, 3fp; window)
```

Repeat with half the step and compare the amplitudes: the step controls
the temporal error, and nothing in the solver estimates it (see
[Choosing the rule, the step and the record](@ref)).


## A JPA against WRspice at the signal frequency

The JPA of the first example of the manual, pumped near its resonance,
solved in time and compared with a WRspice transient simulation of the
same circuit at the signal frequency, and with harmonic balance. WRspice
comes with the [XicTools_jll](https://github.com/JuliaBinaryWrappers/XicTools_jll.jl/)
package on x86_64 Linux; elsewhere install it and use
`JosephsonCircuits.wrspice_cmd()` for the executable.

The pump alone is solved once, for 200 ns at eighty steps per pump
period with the `1 - sech(t/10 ns)` rise of the WRspice input. The
signals are then directions of the tangent about that recorded
trajectory: a unit current at each signal frequency, its cosine and its
sine, with the same rise, all on one pass. The reflection at a signal
frequency comes from the last pump period exactly as
`JosephsonCircuits.wrspice_calcS_paramp` reads it: the cosine and sine responses
combined as `cos + i sin` carry the signal at its own frequency and the
idler at the negative of its own, so demodulated at the signal
frequency the idler is periodic in the pump and averages to nothing
over one period, and a Norton current `I` into the 50 ohm port gives
`V = 25 I (1 + S11)`.

```julia
using JosephsonCircuits, XicTools_jll, Plots

circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1000e-12)), (:cj, 2, 0, Capacitor(1000e-15))])
fp, Ip = 4.75001e9, 0.00565e-6
fs = (4.5:0.01:5.0)*1e9

# harmonic balance
jpa = hbsolve(2pi*fs, (2pi*fp,), [(mode = (1,), port = 1, current = Ip)], (8,), (16,), circuit)
hbgain = 10*log10.(abs2.(jpa.linearized.S(outputmode = (0,), outputport = 1,
    inputmode = (0,), inputport = 1, freqindex = :)))

# the pump alone in time, eighty steps per pump period, WRspice's rise
rise(t) = 1 - 2/(exp(t/10e-9) + exp(-t/10e-9))
problem = transientproblem(circuit; sources = [TransientSource(1, t -> 2Ip*rise(t)*cospi(2fp*t))])
steps = 80
dt = 1/(steps*fp)
pump = transientsolve(problem, (0.0, 76000dt); dt, method = GaussLegendre(), record = :phases)

# the signals as directions of the tangent: unit cosine and sine at each frequency
currents = zeros(1, length(pump.times), 2length(fs))
for (k, f) in enumerate(fs)
    currents[1, :, 2k - 1] .= rise.(pump.times) .* cospi.(2f .* pump.times)
    currents[1, :, 2k] .= rise.(pump.times) .* sinpi.(2f .* pump.times)
end
signal = transienttangent(pump, currents)

# the reflection over the last pump period
reflection(v, f, t) = 2*sum(v .* cispi.(-2f .* t))/length(t)/50 - 1
last = length(pump.times) - steps + 1:length(pump.times)
S11 = [reflection(signal.voltage[1, last, 2k - 1] .+ im .* signal.voltage[1, last, 2k], fs[k], pump.times[last])
    for k in eachindex(fs)]
tdgain = 10*log10.(abs2.(S11))

# WRspice: the pump alone, then the pump with a small sine and a small cosine at each frequency
netlist = JosephsonCircuits.exportnetlist(circuit)
input = JosephsonCircuits.wrspice_input_paramp(netlist.netlist, 2pi*fs, 2pi*fp, 2Ip, (0, 1), (0, 1);
    stepsperperiod = steps)
output = JosephsonCircuits.spice_run(input, XicTools_jll.wrspice())
wrgain = 10*log10.(abs2.(JosephsonCircuits.wrspice_calcS_paramp(output, 2pi*fs, netlist.Nnodes;
    stepsperperiod = steps).S11))

plot(fs/1e9, hbgain; label = "harmonic balance", xlabel = "Frequency (GHz)", ylabel = "Gain (dB)")
plot!(fs/1e9, tdgain; label = "JosephsonCircuits.jl transient", seriestype = :scatter)
plot!(fs/1e9, wrgain; label = "WRspice", seriestype = :scatter, marker = :x)
```

On one core the pump solve and the tangent of a hundred and two
directions take about ten seconds each, and the hundred and three
WRspice simulations about twenty five. The three agree:

| Signal (GHz) | Harmonic balance (dB) | Transient (dB) | WRspice (dB) |
|---|---|---|---|
| 4.50 | 0.0027 | 0.0027 | 0.0027 |
| 4.70 | 0.7108 | 0.7102 | 0.7120 |
| 4.73 | 4.1802 | 4.1781 | 4.1874 |
| 4.75 | 13.3023 | 13.1726 | 13.2485 |
| 4.77 | 4.1857 | 4.1843 | 4.1930 |
| 4.80 | 0.7115 | 0.7115 | 0.7130 |
| 4.90 | 0.0193 | 0.0194 | 0.0193 |

The two time domain results are within 0.03 dB of each other and of
harmonic balance across the band, and within 0.08 dB at the gain peak,
where the amplifier sits nearest its threshold and the gain has not
finished settling at 200 ns, as the next comparison shows.

The port response resolved in time at one signal frequency, 4.76 GHz,
is the same demodulation slid along the record, one pump period at a
time, through the rise of the pump. From WRspice it is the difference
between the runs with and without the signal, which the tangent gives
directly:

```julia
k = findfirst(==(4.76e9), fs)
Is = 1e-13
vw = ((output[2k + 1].values["V"][1, :] .- output[1].values["V"][1, :]) .+
    im .* (output[2k].values["V"][1, :] .- output[1].values["V"][1, :])) ./ Is
vt = signal.voltage[1, :, 2k - 1] .+ im .* signal.voltage[1, :, 2k]
js = steps:steps:length(pump.times)
envelope(v) = [10*log10(abs2(reflection(v[j - steps + 1:j], fs[k], pump.times[j - steps + 1:j]))) for j in js]
plot(pump.times[js]*1e9, envelope(vt); label = "JosephsonCircuits.jl transient",
    xlabel = "Time (ns)", ylabel = "Signal gain (dB)")
plot!(pump.times[js]*1e9, envelope(vw); label = "WRspice", linestyle = :dash)
```

| Time (ns) | Transient (dB) | WRspice (dB) |
|---|---|---|
| 4.2 | -1.093 | -1.093 |
| 12.6 | -4.650 | -4.650 |
| 21.1 | -3.255 | -3.255 |
| 42.1 | 0.869 | 0.872 |
| 63.2 | 5.319 | 5.333 |
| 84.2 | 7.734 | 7.757 |
| 109.5 | 8.454 | 8.483 |
| 160.0 | 8.153 | 8.176 |
| 197.9 | 8.181 | 8.205 |

The signal is first absorbed as the pump rises, then amplified as it
passes threshold, overshoots and settles, and the two simulations
follow each other to 0.03 dB at every time.

