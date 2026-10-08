# A Josephson transmission line with ten pulsed signals

The line below is a small demonstration circuit, not a calibrated
travelling wave amplifier. A pump at 7.5 GHz turns on smoothly and stays
on while a 100 ns pulse of ten signals between 4 and 6.4 GHz passes. The
solve keeps only the port waveforms and the final state, and the outgoing
signals, the pump's third harmonic and an intermodulation product are read
by demodulating the port 2 wave through a smooth window.

```text
 1                      2                            65
 o---[Lj1 || Cj1]-------o---[Lj2 || Cj2]--- ... -----o
 |                      |                            |
[P1 || Cg1]           [Cg2]                   [Cend || P2]
 |                      |                            |
 o----------------------o----------------------------o
 0
```

```@example transientline
using JosephsonCircuits

function transientline(cells)
    circuit = Any[(:P1, 1, 0, Port(1; Z0 = 50.0))]
    for k in 1:cells
        push!(circuit, (Symbol(:Lj, k), k, k + 1, JosephsonJunction(100e-12)))
        push!(circuit, (Symbol(:Cj, k), k, k + 1, Capacitor(20e-15)))
        push!(circuit, (Symbol(:Cg, k), k, 0, Capacitor(40e-15)))
    end
    push!(circuit, (:Cend, cells + 1, 0, Capacitor(40e-15)),
        (:P2, cells + 1, 0, Port(2; Z0 = 50.0)))
    return Circuit(circuit)
end

rise(t, width) = t <= 0 ? 0.0 : t >= width ? 1.0 : sinpi(t/(2width))^2
pulse(t) = rise(t - 20e-9, 15e-9)*rise(120e-9 - t, 15e-9)

cells, ntones = 64, 10
ws = 2pi .* collect(range(4e9, 6.4e9; length = ntones))
phases = [pi*j*(j - 1)/ntones for j in 1:ntones]
wp, Ip, Is = 2pi*7.5e9, 2e-6, 5e-9
# the drive as a closure over its values, which the solver calls at every step
function drive(wp, Ip, ws, phases, Is)
    return function (t)
        pump = Ip*rise(t, 15e-9)*cos(wp*t)
        signals = sum(cos(ws[j]*t + phases[j]) for j in eachindex(ws))
        return pump + Is*pulse(t)*signals
    end
end

lineproblem(cells) = transientproblem(transientline(cells);
    sources = [TransientSource(1, drive(wp, Ip, ws, phases, Is))])
nothing # hide
```

The full 64-cell, 140 ns run and its demodulation are:

```julia
problem = lineproblem(cells)
solution = transientsolve(problem, (0.0, 140e-9); dt = 2e-12)

window(t) = 40e-9 <= t <= 100e-9 ? sinpi((t - 40e-9)/60e-9)^2 : 0.0
for w in ws
    a = transientdemodulate(solution, 2, w; window)
    println("$(round(w/(2pi*1e9); digits = 3)) GHz: $(abs(a)) sqrt(W) at $(angle(a)) rad")
end
idler = transientdemodulate(solution, 2, 2wp - first(ws); window)
third = transientdemodulate(solution, 2, 3wp; window)
```

Repeat with half the step and compare the amplitudes: the step controls
the temporal error, and nothing in the solver estimates it (see
[choosing a step](../transient.md#Choosing-the-rule,-the-step-and-the-record)).

## A small executable check

The documentation build uses four cells and runs through the start of the
signal pulse. It checks the same circuit builder, all ten source tones,
the transient solve, and demodulation without the full device's cost.
The short record and smaller circuit are a syntax/numerics check, not a
prediction of the full line's gain or spectral resolution.

```@example transientline
small = transientsolve(lineproblem(4), (0.0, 24e-9); dt = 4e-12)
@assert all(isfinite, small.outgoing)
@assert pulse(23e-9) > 0
small_window(t) = 21e-9 <= t <= 24e-9 ? sinpi((t - 21e-9)/3e-9)^2 : 0.0
measured = transientdemodulate(small, 2, first(ws); window = small_window)
@assert isfinite(measured) && abs(measured) > 0
nothing # hide
```
