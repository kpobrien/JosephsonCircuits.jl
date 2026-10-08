# Examples

These recipes model complete devices. Start with the
[quickstart](quickstart.md) for a small first calculation. Each recipe
lists its dependencies and draws its circuit: `[x]` is the component
named `x` in the code, `||` joins components in parallel, and node `0` is
ground.

| Example | What it calculates |
|---|---|
| [Josephson parametric amplifier (JPA)](recipes/jpa.md) | The reflection gain of a current-pumped JPA, against WRspice |
| [JPA with a frequency dependent environmental impedance](recipes/environment.md) | The same JPA through a cable from a mismatched source, the pump-off reflection against a closed form |
| [Multi-tone Fourier grids](recipes/multitone.md) | The retained modes and evaluation grids of one, two and three pump tones, and their refinement |
| [Double-pumped JPA](recipes/double-pump.md) | The gain between two pump tones, against WRspice |
| [Flux-pumped JPA](recipes/flux-pump.md) | Three-wave mixing in a flux-biased SQUID, against WRspice, and the resonance against the bias |
| [SNAIL parametric amplifier](recipes/snail.md) | Pump-on and pump-off gain near a Kerr-free point, against WRspice |
| [Josephson traveling wave parametric amplifier (JTWPA)](recipes/traveling-wave.md) | Gain, idlers, quantum efficiency and commutator of a resonantly phase-matched line, the pump-off line against the chain matrices of its cells |
| [Floquet JTWPA](recipes/floquet.md) | A tapered line, lossless and with dielectric loss, the pump-off line against the chain matrices of its cells |
| [Impedance-engineered JPA](recipes/lesa.md) | A snake amplifier built from nested subcircuits, signal and idler in power units, the pump-off reflection against a closed form |
| [Design parameter sensitivities](recipes/sensitivities.md) | Total and frozen-pump derivatives against finite differences, as gain derivatives |
| [Sensitivity to frequency dependent scattering parameters](recipes/scattering-sensitivities.md) | The derivative of a reflection phase against the round-trip propagation phase |
| [Direct current](recipes/dc.md) | Static flux and average voltage, a current source's sign, and a biased junction in time |
| [Pumped finite-window noise](recipes/pumped-noise.md) | The gain and quantum noise of a pumped amplifier in a finite time window |
| [Fit a pumped amplifier](recipes/fitted-amplifier.md) | A causal fit of a pumped JPA's conversion and noise, measured in time against the circuit |
| [Pulsed Josephson line](recipes/transient-line.md) | Ten pulsed signals through a Josephson line in time |
| [Transient WRspice comparison](recipes/transient-wrspice.md) | A JPA in time against WRspice |

The [stability guide](stability.md) works through a JPA, mode tracking
across a sweep, a traveling-wave amplifier's spatial mode profile and a
delay circuit, telling an operating point that converged from one that is
dynamically stable. The [transient guide](transient.md) solves a circuit
in time and differentiates its response, the
[quantum noise guide](transientnoise.md) measures the vacuum and thermal
noise of a passive two-port, and the [scattering guide](scattering.md)
fits scattering data and simulates the result.
