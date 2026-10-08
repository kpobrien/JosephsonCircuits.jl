# JPA with a frequency dependent environmental impedance

The current-pumped JPA of the [first example](jpa.md), connected to a
40 Ω source through 5 cm of 50 Ω cable. The JPA no longer sees a matched
environment but the source transformed by the cable, a different complex
impedance at the signal, the idler and every pump harmonic, so the
operating point moves as well as the signal response. Requires
`JosephsonCircuits` and `Plots`.

```text
 src                 1                 2
 o------[line]-------o------[cc]-------o--------+
 |                                     |        |
[p1]                                 [jj]     [cj]
 |                                     |        |
 o-------------------------------------o--------+
 0
```

The source is the port. A [`Port`](@ref) owns a matched source and load of
its reference impedance, so `Port(1; Z0 = 40.0)` is the 40 Ω source, and
its scattering parameters are referenced to it: `abs2(S11)` is the
reflected power over the power the source makes available. The cable is a
[`TransmissionLine`](@ref) between the source and the JPA, so the pump
injected at the source reaches the junction through it, and the solver
produces the mismatch between source, cable and JPA. The source's 40 Ω and
the cable's 50 Ω, 5 cm and phase velocity of `2e8` m/s are assumed; the
JPA and its pump are those of the first example.

```@example environment
using JosephsonCircuits
Cc, Lj, Cj = 100.0e-15, 1000.0e-12, 1000.0e-15
Z0, Rs, len, vp = 50.0, 40.0, 0.05, 2.0e8

ideal = Circuit([
    (:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:cc, 1, 2, Capacitor(Cc)),
    (:jj, 2, 0, JosephsonJunction(Lj)),
    (:cj, 2, 0, Capacitor(Cj))])
cable = Circuit([
    (:p1, "src", 0, Port(1; Z0 = Rs)),
    (:line, "src", 1, TransmissionLine(Z0, len; vp)),
    (:cc, 1, 2, Capacitor(Cc)),
    (:jj, 2, 0, JosephsonJunction(Lj)),
    (:cj, 2, 0, Capacitor(Cj))])
nothing # hide
```

With the pump off the circuit is linear. The JPA's impedance `ZJ`, its
junction an inductance at zero phase, seen through the line,
`Zin = Z0 (ZJ + i Z0 tan(ωℓ/vp))/(Z0 + i ZJ tan(ωℓ/vp))`, reflects
`(Zin - Rs)/(Zin + Rs)` back to the source:

```@example environment
ws = 2pi .* [4.6e9, 4.75e9, 4.9e9]
off = hblinsolve(ws, cable; keyedarrays = false)
ZJ(w) = 1/(im*w*Cc) + 1/(im*w*Cj + 1/(im*w*Lj))
Zin(w) = Z0*(ZJ(w) + im*Z0*tan(w*len/vp))/(Z0 + im*ZJ(w)*tan(w*len/vp))
closed = [(Zin(w) - Rs)/(Zin(w) + Rs) for w in ws]
@assert isapprox(off.S[1, 1, :], closed; atol = 1e-12)
(solver = off.S[1, 1, :], closed = closed)
```

The pumped sweep, with the pump current of the first example:

```@example environment
using Plots
wsweep = 2pi .* (4.5:0.002:5.0) .* 1e9
wp = (2pi*4.75001e9,)
sources = [(mode = (1,), port = 1, current = 0.00565e-6)]
matched = hbsolve(wsweep, wp, sources, (8,), (16,), ideal)
@assert matched.nonlinear.solverinfo.converged
throughcable = hbsolve(wsweep, wp, sources, (8,), (16,), cable)
@assert throughcable.nonlinear.solverinfo.converged
gain(sol) = 10 .* log10.(abs2.(sol.linearized.S((0,), 1, (0,), 1, :)))
plot(wsweep ./ (2pi*1e9), [gain(matched) gain(throughcable)],
    label = ["matched 50 ohm source" "40 ohm source through the cable"],
    xlabel = "Frequency (GHz)", ylabel = "Gain (dB)")
```

The matched JPA peaks at 13.3 dB at the pump frequency; through the cable,
at the same pump current, it peaks at 3.2 dB. At the pump frequency the
JPA sees the source through the cable as about `58 + 9i` Ω instead of
50 Ω, a weaker coupling and a shifted resonance, and the pump, which also
reaches it from a source of lower available power, leaves it far from its
gain point:

```@example environment
Zenv(w) = Z0*(Rs + im*Z0*tan(w*len/vp))/(Z0 + im*Rs*tan(w*len/vp))
(matched = maximum(gain(matched)), throughcable = maximum(gain(throughcable)),
    environment = Zenv(wp[1]))
```

Retune the pump for the environment rather than correcting the readout.
