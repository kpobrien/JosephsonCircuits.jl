# Examples

These recipes cover complete device models. Start with the
[quickstart](quickstart.md) if you want a small first calculation. Each page
lists its dependencies and what to inspect in the output.

The original section headings below remain as links so existing bookmarks
continue to identify the corresponding example.

## Josephson parametric amplifier (JPA)

[Josephson parametric amplifier (JPA)](recipes/jpa.md). Calculate reflection gain near the resonance of a current-pumped JPA, then compare with WRspice. The signal-frequency sweep uses the small-signal response about the strong pump. Refine both harmonic limits before interpreting the peak gain.

## JPA with a frequency dependent environmental impedance

[JPA with a frequency dependent environmental impedance](recipes/environment.md). Compare the same JPA with a matched environment and a mismatched cable. This example changes the pump loading as well as the signal response; it is not just a phase correction to the output.

## Double-pumped Josephson parametric amplifier (JPA)

[Double-pumped Josephson parametric amplifier (JPA)](recipes/double-pump.md). Drive the JPA with two independent strong tones and inspect the response between them. Each pump adds a Fourier-grid dimension, so refine the harmonic limits while watching memory use.

## Flux-pumped Josephson parametric amplifier (JPA)

[Flux-pumped Josephson parametric amplifier (JPA)](recipes/flux-pump.md). Bias a SQUID through a mutual inductor and apply a pump near twice its resonance. The DC and three-wave-mixing options retain the modes needed by the flux bias and pump. The final sweep maps the pump-off resonance against bias.

## SNAIL Parametric Amplifier

[SNAIL Parametric Amplifier](recipes/snail.md). Compare pump-on and pump-off gain for an explicitly modeled SNAIL. The quadratic nonlinearity supports three-wave mixing. Similar resonance frequencies in the two curves indicate operation near the chosen Kerr-free point.

## Josephson traveling wave parametric amplifier (JTWPA)

[Josephson traveling wave parametric amplifier (JTWPA)](recipes/traveling-wave.md). Build a resonant-phase-matched JTWPA from repeated cells. Inspect forward and reverse gain, conversion to idlers, quantum efficiency, and the commutator error. This is a full device example; reduce the cell count for a quick syntax check, but do not expect the same gain.

## Floquet JTWPA

[Floquet JTWPA](recipes/floquet.md). Taper the unit-cell parameters of a traveling-wave amplifier, then add dielectric loss. The second section reuses `floquetcircuit` from the first. Compare gain and normalized quantum efficiency, and check convergence of the harmonic truncations.

## Floquet JTWPA with dissipation

[Floquet JTWPA with dissipation](recipes/floquet.md#Floquet-JTWPA-with-dissipation). Taper the unit-cell parameters of a traveling-wave amplifier, then add dielectric loss. The second section reuses `floquetcircuit` from the first. Compare gain and normalized quantum efficiency, and check convergence of the harmonic truncations.

## Impedance-engineered JPA

[Impedance-engineered JPA](recipes/lesa.md). Build an impedance-engineered JPA from nested snake and SQUID subcircuits. The signal and conjugate-idler curves are plotted in power units, so the idler conversion includes the absolute frequency ratio.

## Design parameter sensitivities

[Design parameter sensitivities](recipes/sensitivities.md). Differentiate the JPA response with respect to its design parameters. The small executable example first checks a derivative against finite differences; the plotting recipe then expresses gain derivatives in dB per fractional parameter change.

## Sensitivity to frequency dependent scattering parameters

[Sensitivity to frequency dependent scattering parameters](recipes/scattering-sensitivities.md). Differentiate the reflection phase with respect to a matched line length. The analytic reference is the round-trip propagation phase. Block derivatives describe the data at one design point; rebuilding a block is required when its captured parameters change.

## Direct current

[Direct current](recipes/dc.md). Distinguish static flux from average voltage using a resistor carrying DC. The executable example checks Ohm's law without adding an artificial inductive path to ground.

## Transient and noise examples

- [Transient solve and response derivatives](transient.md)
- [A Josephson line with ten pulsed signals](recipes/transient-line.md)
- [A JPA compared with WRspice in time](recipes/transient-wrspice.md)
- [Vacuum and thermal noise of a passive two-port](transientnoise.md)
- [Fit scattering data and simulate the result](scattering.md)
