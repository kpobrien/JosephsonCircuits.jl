```@raw html
---
layout: home
hero:
  name: JosephsonCircuits.jl
  text: Simulate superconducting circuits
  tagline: Harmonic balance and transient simulation of Josephson circuits. Calculate gain, frequency conversion, stability, noise, and sensitivities.
  image:
    light: /logo.svg
    dark: /logo-dark.svg
    alt: JosephsonCircuits.jl
  actions:
    - theme: brand
      text: Get started
      link: /quickstart
    - theme: alt
      text: Examples
      link: /examples
    - theme: alt
      text: View on GitHub
      link: https://github.com/kpobrien/JosephsonCircuits.jl
features:
  - title: Build a circuit
    details: Connect junctions, linear components, subcircuits, scattering models, and transmission lines.
    link: /circuits
  - title: Harmonic balance
    details: Find a driven operating point and calculate its small-signal response, noise, and sensitivities.
    link: /harmonicbalance
  - title: Transient simulation
    details: Simulate pulses and multiple strong tones. Measure the response and quantum noise in time windows.
    link: /transient
  - title: API reference
    details: Constructors, solvers, result types, and advanced interfaces.
    link: /reference
---
```

JosephsonCircuits.jl simulates nonlinear microwave circuits containing
Josephson junctions. It is designed for devices such as Josephson
parametric amplifiers and traveling-wave amplifiers with thousands of
components.

Use **harmonic balance** for a steady driven state and its small-signal
response. Use **transient simulation** for pump turn-on, pulses, and drives
with too many independent tones for a practical harmonic grid. Both
analyses use the same circuit description, subject to the
[component restrictions](circuits.md#Supported-analyses).

The noise calculations propagate small fluctuations about a classical
operating point or trajectory. They describe Gaussian quantum noise;
they do not evolve the full quantum state of a nonlinear circuit.

## Installation

In a Julia environment, install the registered release with:

```julia
using Pkg
Pkg.add("JosephsonCircuits")
```

Then follow the [first calculation](quickstart.md). The plotting examples
also need `Plots`; examples using external solvers or WRspice list their
additional dependencies.

The version selector distinguishes the released manual from development
versions. Development documentation can describe APIs that are not yet in
the registered release. To use the current development branch:

```julia
using Pkg
Pkg.add(name = "JosephsonCircuits", rev = "main")
```

If an example does not match your installation, check `Pkg.status()` and
`VERSION` first. Include those values and a minimal reproducer in a bug
report. See the [migration guide](migration.md) for interface changes.

## A first calculation

The [quickstart](quickstart.md) builds a one-port JPA, solves its pumped
operating point, and reads its reflection gain. It also shows how to apply
the same pump in time, including the Fourier-amplitude convention.

## How the documentation is organized

| Task | Start here |
|---|---|
| Run a first simulation | [Quickstart](quickstart.md) |
| Check units or normalization | [Conventions](conventions.md) |
| Define components and connections | [Circuits](circuits.md) |
| Include measured data or a fitted amplifier | [Scattering blocks](scattering.md) |
| Find a driven operating point and gain | [Harmonic balance](harmonicbalance.md) |
| Test temporal stability of a driven state | [Harmonic-balance stability](stability.md) |
| Warm ports, a readout chain's efficiency, the noise reaching a qubit | [Noise at the ports](portnoise.md) |
| Simulate pulses | [Transient simulation](transient.md) |
| Calculate fluctuations in temporal modes | [Quantum noise in time](transientnoise.md) |
| Reuse work or run on a GPU | [Performance](performance.md) |
| Supply another numerical solver | [Using other solvers](interop.md) |
| Understand the equations and algorithms | [HB theory](harmonicbalancetheory.md), [transient theory](transienttheory.md) |

## Supported calculations

- Nonlinear steady states with one or more strong drive tones.
- Small-signal scattering, frequency conversion, quantum efficiency, and
  sensitivities about a driven state.
- Temporal stability poles about a driven state, with method-specific limits
  and convergence checks in the [stability guide](stability.md).
- Linear circuit responses and symbolic circuit matrices.
- Transient trajectories, their tangent and adjoint responses, and noise
  in selected temporal modes.

The [examples](examples.md) include comparisons with WRspice. Reference 6
below describes comparisons of the frequency-domain solver with ADS and
Fourier analysis of WRspice simulations.

## Contributing

We welcome contributions in the form of issues/bug reports or pull requests. This project uses the [MIT open source license](https://opensource.org/license/MIT). You retain the copyright to any code you contribute.

## References

1. Andrew J. Kerman "Efficient numerical simulation of complex Josephson quantum circuits" [arXiv:2010.14929 (2020)](https://doi.org/10.48550/arXiv.2010.14929) 
2. Ji&#345;&#237; Vlach and Kishore Singhal "Computer Methods for Circuit Analysis and Design" 2nd edition, [Springer New York, NY (1993)](https://link.springer.com/book/9780442011949)
3. Stephen A. Maas "Nonlinear Microwave and RF Circuits" 2nd edition, [Artech House (1997)](https://us.artechhouse.com/Nonlinear-Microwave-and-RF-Circuits-Second-Edition-P1097.aspx)
4. Jos&#233; Carlos Pedro, David E. Root, Jianjun Xu, and Lu&#237;s C&#243;timos Nunes. "Nonlinear Circuit Simulation and Modeling: Fundamentals for Microwave Design" The Cambridge RF and Microwave Engineering Series, [Cambridge University Press (2018)](https://www.cambridge.org/core/books/nonlinear-circuit-simulation-and-modeling/1705F3B449B4313A2BE890599DAC0E38)
5. David E. Root, Jan Verspecht, Jason Horn, and Mihai Marcu. "X-Parameters: Characterization, Modeling, and Design of Nonlinear RF and Microwave Components" The Cambridge RF and microwave engineering series, [Cambridge University Press (2013)](https://www.cambridge.org/sb/academic/subjects/engineering/rf-and-microwave-engineering/x-parameters-characterization-modeling-and-design-nonlinear-rf-and-microwave-components)
6. Kaidong Peng, Rick Poore, Philip Krantz, David E. Root, and Kevin P. O'Brien "X-parameter based design and simulation of Josephson traveling-wave parametric amplifiers for quantum computing applications" [IEEE International Conference on Quantum Computing & Engineering (QCE22) (2022)](http://arxiv.org/abs/2211.05328)

## Philosophy

The motivation for developing this package is to simulate the gain and noise performance of ultra low noise amplifiers for quantum computing applications such as the [Josephson traveling-wave parametric amplifier](https://www.science.org/doi/10.1126/science.aaa8525), which have thousands of linear and nonlinear circuit elements. 

We prioritize speed (including compile time and time to first use), simplicity, and scalability.

## Related packages and software
* [Xyce.jl](https://github.com/JuliaComputing/Xyce.jl) provides a wrapper for [Xyce](https://xyce.sandia.gov/), the open source parallel circuit simulator from Sandia National Laboratories which can perform time domain and harmonic balance method simulations.
* [NgSpice.jl](https://github.com/JuliaComputing/Ngspice.jl) and [LTspice.jl](https://github.com/cstook/LTspice.jl) provide wrappers for [NgSpice](http://ngspice.sourceforge.net/) and [LTspice](https://www.analog.com/en/design-center/design-tools-and-calculators/ltspice-simulator.html), respectively.  
* [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) supports time domain circuit simulations from [scratch](https://mtk.sciml.ai/stable/tutorials/acausal_components) and using their [standard library](https://docs.sciml.ai/ModelingToolkitStandardLibrary/stable/tutorials/rc_circuit/)
* [ACME.jl](https://github.com/HSU-ANT/ACME.jl) simulates electrical circuits in the time domain with an emphasis on audio effect circuits.
* [Cedar EDA](https://cedar-eda.com) is a Julia-based commercial cloud service for circuit simulations.
* [Keysight ADS](https://www.keysight.com/us/en/products/software/pathwave-design-software/pathwave-advanced-design-system.html), [Cadence AWR](https://www.cadence.com/en_US/home/tools/system-analysis/rf-microwave-design/awr-microwave-office.html), [Cadence Spectre RF](https://www.cadence.com/en_US/home/tools/custom-ic-analog-rf-design/circuit-simulation/spectre-rf-option.html), and [Qucs](http://qucs.sourceforge.net/) are capable of time and frequency domain analysis of nonlinear circuits. [WRSPICE](http://wrcad.com/wrspice.html) performs time domain simulations of Josephson junction containing circuits and frequency domain simulations of linear circuits.

## Funding
We gratefully acknowledge funding from the [AWS Center for Quantum Computing](https://aws.amazon.com/blogs/quantum-computing/announcing-the-opening-of-the-aws-center-for-quantum-computing/) and the [MIT Center for Quantum Engineering (CQE)](https://cqe.mit.edu/).
