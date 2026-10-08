# Sensitivity to frequency dependent scattering parameters

Differentiate the reflection phase with respect to a matched line length. The analytic reference is the round-trip propagation phase. Block derivatives describe the data at one design point; rebuilding a block is required when its captured parameters change.

The plotting code requires `Plots` in addition to `JosephsonCircuits`.

```text
 1                  2                  3
 o------[tl]--------o-------[cc]-------o--------+
 |                                     |        |
[p1]                                 [jj]     [c2]
 |                                     |        |
 o-------------------------------------o--------+
 0
```

A [`ScatteringParameters`](@ref) block depends on a design parameter through the analytic derivative its `derivatives` keyword states: here a matched transmission line section in front of the amplifier, with the derivative of the reflection with respect to the line length. A matched lossless line delays the reflection on its way in and out, so the gain does not depend on the length, and the phase of `S11` turns by `-2*w*len/vphase` per fractional change of the length, which the example compares with. A block which states no derivative (measured Touchstone data, say) is parameter independent and costs nothing.

```julia
using JosephsonCircuits
using Plots

vphase = 1.2e8 # the phase velocity of the line in m/s
len = 2.0e-3   # its length in m
# a matched, lossless transmission line: S21 = exp(-im*w*len/vphase)
tlineS(w) = (t = exp(-im*w*len/vphase); [0 t; t 0])
# the analytic derivative of S with respect to the length
dSdlen(w) = (d = -im*w/vphase*exp(-im*w*len/vphase); [0 d; d 0])
line = ScatteringParameters(tlineS; nports = 2, noise = Lossless(),
    derivatives = (len = dSdlen,))
circuit = Circuit(
    [(:p1, 1, 0, Port(1)),
     (:tl, 1, 2, line),
     (:cc, 2, 3, Capacitor(100.0e-15)),
     (:jj, 3, 0, JosephsonJunction(:Lj)),
     (:c2, 3, 0, Capacitor(1000.0e-15))])

p = Dict(:Lj => 1000.0e-12)
ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = (2*pi*4.75001*1e9,)
sources = [(mode=(1,),port=1,current=0.00565e-6)]

# the parameters: the line length the block states a derivative for, and Lj
r = designsensitivities(circuit, p, ws, wp, sources, (8,), (16,))

# the derivative of the phase of S11 with respect to the fractional change
# of the length, len*imag(dS/dlen/S), against the delay of the line there
# and back
S = r.out.linearized.S((0,),1,(0,),1,:)
dphase = len .* imag.(r.dSdp((0,),1,(0,),1,:len,:) ./ S)
delay = -2 .* ws .* len ./ vphase
println("largest difference from -2*w*len/vphase: ",
    round(maximum(abs, dphase .- delay); sigdigits = 2), " rad")
plot(ws/(2*pi*1e9), dphase,
    label="JosephsonCircuits.jl",
    xlabel="Frequency (GHz)",
    ylabel="d arg(S11)/dln(len) (rad)")
plot!(ws/(2*pi*1e9), delay,
    label="-2*w*len/vphase",
    linestyle=:dash)
```

```
largest difference from -2*w*len/vphase: 1.3e-7 rad
```
