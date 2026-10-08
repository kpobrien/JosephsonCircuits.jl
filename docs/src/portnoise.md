# Noise at the ports

A port's termination is the measurement line: the source and load that
send a field into the port and absorb what leaves it. Its temperature sets
the field it sends in; the circuit's losses and blocks set the noise the
circuit adds. This page states port temperatures and reads two quantities
of a readout chain: the measurement efficiency from the device to room
temperature, and the thermal occupation the chain sends back toward the
device.

## Warm terminations

`Port(n; termination = MatchedTermination(temperature = T))` gives a port a
matched termination at the physical temperature `T` in kelvin. It sends in
a thermal field, the symmetrized noise `nbar + 1/2` in each mode, with
`nbar = thermaloccupation(ω, T)` at the mode's frequency, so a signal and
its idler see different occupations. The
[temperature table](conventions.md#Noise-normalization-and-temperature)
gives the temperature of each source of noise. For a line whose noise is
not that of one temperature, put its attenuators in the circuit at their
own temperatures with the port at the top, as below, or give the port the
effective temperature of the field that reaches it
([`effectivetemperature`](@ref)).

A passive circuit whose ports and losses share one temperature is in
equilibrium: every wave leaving it is thermal.

```@example portnoise
using JosephsonCircuits, LinearAlgebra
w = 2pi*3e9
T = 0.3
warm(n) = Port(n; termination = MatchedTermination(temperature = T))
circuit = Circuit([
    (:p1, 1, 0, warm(1)), (:c1, 1, 0, Capacitor(0.3e-12)),
    (:loss, 1, 2, Resistor(30.0; temperature = T)),
    (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, warm(2)),
])
sol = hblinsolve([w], circuit; keyedarrays = false, returnVout = true)
n = thermaloccupation(w, T)
@assert isapprox(sol.Vout[:, :, 1], (n + 1/2)*I; atol = 1e-12)
(occupation = n, nbar = sol.nbar[:, 1])
```

## What the solvers report

Noise is in quanta, the vacuum being half a photon (see
[conventions](conventions.md#Noise-normalization-and-temperature)).

| Output | Meaning |
|---|---|
| `QE` | photon gain over twice the noise at the output, every input in its state, the signal's own included |
| `QEideal` | the ideal phase-preserving amplifier at the same gain, its inputs in the vacuum |
| `nbar` | the occupation of the wave leaving each port mode |
| `Vout` | the output covariance, `S*Diagonal(sigma)*S' + Cnoise`, with `sigma` the `nbar + 1/2` of every input (`returnVout = true`) |
| `Cnoise` | the noise the circuit adds, never its ports' (`returnCnoise = true`) |
| `Snoise`, `channeltemperatures` | the transfer from each internal noise channel, and the channel's temperature (`returnSnoise = true`) |

`QE` and `nbar` read the same noise at an output, `v = Vout[i,i]`:
`QE = G/(2v)` and `nbar = v - 1/2`. With a warm source the QE counts the
thermal photons arriving with the signal, as a measurement of the signal to
noise ratio does; with the source cold it is the device's own efficiency,
and `QE/QEideal` its distance from the phase-preserving quantum limit.
Every source `s` puts `(nbar_s + 1/2)*abs2(t_s)` into `v`, a port mode
through `S` and an internal channel through `Snoise`: that is the noise
budget of an output.

## The occupation reaching the device

An input line brings the room-temperature load to the device through
attenuators at the fridge's stages. The occupation reaching the device is
the cascade of the stages, each stage's occupation weighted by its
transmission to the device.

```@example portnoise
att(dB, T) = ScatteringParameters([0.0 10^(-dB/20); 10^(-dB/20) 0.0];
    zref = 50.0, noise = ThermalEquilibrium(T))
at(T) = MatchedTermination(temperature = T)
line = Circuit([
    (:device, 1, 0, Port(1; termination = at(0.02))),
    (:a1, 1, 2, att(20, 0.02)), (:a2, 2, 3, att(20, 0.1)),
    (:a3, 3, 4, att(20, 4.0)), (:room, 4, 0, Port(2; termination = at(300.0))),
])
w = 2pi*5e9
sol = hblinsolve([w], line; returnSnoise = true)
nback = sol.nbar((0,), 1, 1)
A = 0.01
cascade = thermaloccupation(w, 300.0)*A^3 + thermaloccupation(w, 4.0)*(1 - A)*A^2 +
    thermaloccupation(w, 0.1)*(1 - A)*A + thermaloccupation(w, 0.02)*(1 - A)
@assert isapprox(nback, cascade; rtol = 1e-10)
# the budget: each source's thermal photons, weighted by its transmission;
# a channel is named for its block and port
blockof(channel) = String(first(split(channel, '/')))
budget = Dict("room" => thermaloccupation(w, 300.0)*abs2(sol.S((0,), 1, (0,), 2, 1)))
for (c, T) in zip(sol.Snoise.component, sol.channeltemperatures)
    photons = thermaloccupation(w, T)*abs2(sol.Snoise((0,), c, (0,), 1, 1))
    budget[blockof(c)] = get(budget, blockof(c), 0.0) + photons
end
@assert isapprox(sum(values(budget)), nback; rtol = 1e-10)
(nbar = nback, temperature = effectivetemperature(w, nback),
    budget = sort(collect(budget); by = last, rev = true))
```

The device sees the field of a bath at the effective temperature, above
the mixing chamber's 20 mK; most of it comes from the 4 K and 0.1 K stages
and from the room-temperature load.

## A readout chain

A JPA read out through a circulator, an isolator at the mixing chamber, a
cable at 4 K and a HEMT stated by its noise temperature, to a room
temperature load. The quantum efficiency from the device's port to the room
temperature port is the measurement efficiency of the chain, and its
budget says what limits it: the idler of the JPA, the losses in front of the
HEMT, and the HEMT.

```@example portnoise
iso = ScatteringParameters([0.0 10^(-20/20); 10^(-0.45/20) 0.0]; zref = 50.0,
    noise = ThermalEquilibrium(0.02))
circulator = ScatteringParameters([0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0];
    zref = 50.0, noise = Lossless())
fs, fp, ip = 4.748e9, 4.75e9, 0.0058e-6
ws = 2pi*fs
GH, TN = 1e4, 2.0
nadd = noisequanta(ws, TN)
# the HEMT's input emits backward like a matched load at 4 K
hemt = ScatteringParameters([0.0 0.0; sqrt(GH) 0.0]; zref = 50.0,
    noise = NoiseCovariance([thermaloccupation(ws, 4.0)+1/2 0.0; 0.0 GH*nadd]))
chain = Circuit([
    (:device, 1, 0, Port(1; termination = at(0.02))),
    (:circ, 1, 2, 3, circulator),
    (:cc, 2, 4, Capacitor(100e-15)), (:jj, 4, 0, JosephsonJunction(1000e-12)),
    (:cj, 4, 0, Capacitor(1000e-15)),
    (:iso, 3, 5, iso), (:cable, 5, 6, att(1, 4.0)), (:hemt, 6, 7, hemt),
    (:room, 7, 0, Port(2; termination = at(300.0))),
])
sol = hbsolve([ws], (2pi*fp,), [(mode = (1,), port = 1, current = ip)], (8,), (16,),
    chain; returnSnoise = true).linearized
QE = sol.QE((0,), 2, (0,), 1, 1)
G = abs2(sol.S((0,), 2, (0,), 1, 1))
# every source's share of the noise at the output, referred to the input:
# each port in each mode, and each block's channels in every mode
wmode(m) = ws + m[1]*2pi*fp
share(n, t) = (n + 1/2)*abs2(t)/G
budget = Dict(string("port ", p, " ", m) =>
    share(thermaloccupation(wmode(m), p == 1 ? 0.02 : 300.0), sol.S((0,), 2, m, p, 1))
    for p in (1, 2) for m in sol.modes)
for (c, T) in zip(sol.Snoise.component, sol.channeltemperatures), m in sol.modes
    photons = share(thermaloccupation(wmode(m), T), sol.Snoise(m, c, (0,), 2, 1))
    budget[blockof(c)] = get(budget, blockof(c), 0.0) + photons
end
@assert isapprox(sum(values(budget)), 1/(2QE); rtol = 1e-10)
# against the Friis cascade: the JPA's idler, then each loss and the HEMT
# referred through the gain before it
ns = thermaloccupation(ws, 0.02)
ni = thermaloccupation(abs(ws - 2*2pi*fp), 0.02)
eiso, ecable = 10^(-0.45/10), 10^(-1/10)
GJ = G/(eiso*ecable*GH)
loss(eta, T) = (1 - eta)/eta*(thermaloccupation(ws, T) + 1/2)
friis = (GJ - 1)/GJ*(ni + 1/2) +
    (loss(eiso, 0.02) + loss(ecable, 4.0)/eiso + nadd/(eiso*ecable))/GJ
@assert isapprox(1/(2QE) - (ns + 1/2), friis; rtol = 1e-3)
(jpagain = 10*log10(GJ), QE = QE, fraction = QE/sol.QEideal((0,), 2, (0,), 1, 1),
    systemtemperature = noisetemperature(ws, 1/(2QE)),
    budget = sort(filter(x -> last(x) > 1e-3, collect(budget)); by = last, rev = true))
```

The budget's entries are photons referred to the device's port. They add
to `1/(2QE)`, the signal's own half photon included; the JPA's idler
contributes about half a photon at this gain, the quantum limit of a phase
preserving amplifier. The system noise temperature is the same noise in
kelvin ([`noisetemperature`](@ref)).

## The noise the chain sends back

`nbar` at the device's port is the thermal occupation the chain sends back
toward the device. Here it is the HEMT's backward emission and the 4 K
cable's, through the 20 dB of the isolator, and the isolator's own at the
mixing chamber.

```@example portnoise
nback = sol.nbar((0,), 1, 1)
gi2 = 10^(-20/10)
back = gi2*thermaloccupation(ws, 4.0) + (1 - gi2)*thermaloccupation(ws, 0.02)
@assert isapprox(nback, back; rtol = 1e-8)
(nbar = nback, temperature = effectivetemperature(ws, nback))
```

A circulator of finite isolation would add the JPA's amplified output here,
and a real HEMT emits backward more than a matched load at its stage
temperature; both are blocks with the noise they state.

## In time

The transient's port baths take the same termination temperatures, so
[`transientnoise`](@ref) measures warm ports as harmonic balance does: its
`covariance` is the output covariance `Vout` in the measured modes and its
`addedcovariance` the noise the circuit adds, `Cnoise`, and
[`transientquantumefficiency`](@ref) returns the same QE of a
phase-insensitive mode. See [quantum noise in time](transientnoise.md).
