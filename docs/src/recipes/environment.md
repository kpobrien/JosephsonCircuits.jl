# JPA with a frequency dependent environmental impedance

Compare the same JPA with a matched environment and a mismatched cable. This example changes the pump loading as well as the signal response; it is not just a phase correction to the output.

The plotting code requires `Plots` in addition to `JosephsonCircuits`.

Any component value can be a function of frequency: `FrequencyDependent`
wraps an arbitrary Julia closure of a positive frequency in radians per
second, and both the nonlinear pump solve and the linearized sweep evaluate
it at the magnitude of their own mode frequencies, taking the complex
conjugate at a negative one, since a physical impedance obeys
`Z(-w) = conj(Z(w))`. Here the JPA of the first example sees its
environment through five centimeters of slightly mismatched cable -- a 50 ohm
line terminated by a 40 ohm source -- so the port resistor becomes the
complex input impedance of that line. The transformed environment reshapes
and lowers the gain, and because the pump harmonics feel it too, the
operating point itself shifts, not just the readout.

```julia
using JosephsonCircuits
using Plots

Cc = 100.0e-15
Lj = 1000.0e-12
Cj = 1000.0e-15

# the input impedance of a mismatched cable between source and amplifier
Z0 = 50.0    # cable characteristic impedance, ohms
ZL = 40.0    # source impedance, ohms
len = 0.05   # cable length, meters
vp = 2.0e8   # cable phase velocity, meters per second
Zenv(w) = Z0*(ZL + im*Z0*tan(w*len/vp))/(Z0 + im*ZL*tan(w*len/vp))

amplifier(Zr) = Circuit(
    [(:p1, 1, 0, Port(1; Z0 = Zr)),
     (:cc, 1, 2, Capacitor(Cc)),
     (:jj, 2, 0, JosephsonJunction(Lj)),
     (:cj, 2, 0, Capacitor(Cj))])

ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = (2*pi*4.75001*1e9,)
sources = [(mode=(1,), port=1, current=0.00565e-6)]

@time ideal = hbsolve(ws, wp, sources, (8,), (16,), amplifier(50.0))
@assert ideal.nonlinear.solverinfo.converged
@time cable = hbsolve(ws, wp, sources, (8,), (16,),
    amplifier(FrequencyDependent(Zenv)))
@assert cable.nonlinear.solverinfo.converged

gain(sol) = 10*log10.(abs2.(sol.linearized.S(outputmode=(0,),
    outputport=1, inputmode=(0,), inputport=1, freqindex=:)))
plot(ws/(2*pi*1e9), [gain(ideal) gain(cable)],
    label=["ideal 50 ohm environment" "through mismatched cable"],
    xlabel="Frequency (GHz)", ylabel="Gain (dB)")
```
