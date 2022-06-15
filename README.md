
# JosephsonCircuits.jl 

[JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl) is a high-performance frequency domain simulator for nonlinear circuits containing Josephson junctions, capacitors, inductors, mutual inductors, and resistors. [JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl) simulates the frequency domain behavior using a variant [1] of nodal analysis [2] and the harmonic balance method [3-5] with an analytic Jacobian. Noise performance, quantified by quantum efficiency, is efficiently simulated through an adjoint method.

Frequency dependent circuit parameters are supported to model realistic impedance environments or dissipative components. Dissipation can be modeled by (potentially) frequency dependent resistors or capacitors with an imaginary capacitance proportional to the loss tangent.

[JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl) supports the following:
* Nonlinear simulations in which the user defines a circuit, the drive current, frequency, and number of harmonics and the code calculates the node flux or equivalently node voltage at each harmonic.
* Simulations linearized about the nonlinear operating point calculated above. This effectively simulates the small signal response of a time dependent linear circuit and is useful for simulating parametric amplification and frequency conversion in the undepleted (strong) pump limit. Calculation of X parameters [4-5], which are a generalization of scattering parameters which quantifies how waves with different frequencies interact with the circuit and with each other. 
* Linear simulations of linear circuits. Calculation of node fluxes (or node voltages) and scattering parameters.
* Optional calculation of symbolic capacitance and inverse inductance matrices.

# Examples:
Josephson parametric amplifier (a driven nonlinear LC resonator). Add a circuit diagram. 

```julia
using JosephsonCircuits

@syms R Cc Lj Cj
circuit = Array{Tuple{String,String,String,Any},1}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",R))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))

circuitdefs = Dict(
    Lj =>1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    R => 50.0,
)

@time jpa1 = hbsolve(2*pi*(4.5:0.001:5.0)*1e9,
    2*pi*4.75001*1e9,0.00565e-6,8,8,circuit,circuitdefs,
    pumpports=[1],symfreqvar=w);

using Plots
plot(jpa1.signal.w/(2*pi*1e9),
	10*log10.(abs2.(jpa1.signal.S[
	jpa1.signal.signalindex,
	jpa1.signal.signalindex,:])),
    xlabel="Frequency (GHz)",ylabel="Gain (dB)")
```
```
  0.004098 seconds (68.61 k allocations: 5.766 MiB)
```
![alt text](https://qce.mit.edu/plot.png "Title")

We find excellent agreement with the gain obtained by Fourier analysis of a time domain simulation performed by [WRSPICE](http://wrcad.com/wrspice.html) using the following [netlist]().

We can also sweep the pump frequency.

* LC ladder transmission line. That would be nice because i could compare with the analytic solution for the scattering parameters. I should intentionally impedance mismatch it. Add a circuit diagram.

* JTWPA (with dissipation)

* Floquet JTWPA

# Additional Examples:


# Installation:

Until this package is registered, you can install it by starting Julia using the command:
```
JULIA_PKG_USE_CLI_GIT=true julia
```
then type ] followed by (replacing kpobrien with your GitHub username):
```
add git@github.com:kpobrien/JosephsonCircuits.jl.git
```
or once we make the repository public:
```
]add https://github.com/kpobrien/JosephsonCircuits.jl
```

# References:

1. Andrew J. Kerman "Efficient numerical simulation of complex Josephson quantum circuits" [arXiv:2010.14929 (2020)](https://doi.org/10.48550/arXiv.2010.14929) 
2. Ji&#345;&#237; Vlach and Kishore Singhal "Computer Methods for Circuit Analysis and Design" 2nd edition, [Springer New York, NY (1993)](https://link.springer.com/book/9780442011949)
3. Stephen A. Maas "Nonlinear Microwave and RF Circuits" 2nd edition, [Artech House (1997)](https://us.artechhouse.com/Nonlinear-Microwave-and-RF-Circuits-Second-Edition-P1097.aspx)
4. Jos&#233; Carlos Pedro, David E. Root, Jianjun Xu, and Lu&#237;s C&#243;timos Nunes. "Nonlinear Circuit Simulation and Modeling: Fundamentals for Microwave Design" The Cambridge RF and Microwave Engineering Series, [Cambridge University Press (2018)](https://www.cambridge.org/core/books/nonlinear-circuit-simulation-and-modeling/1705F3B449B4313A2BE890599DAC0E38)
5. David E. Root, Jan Verspecht, Jason Horn, and Mihai Marcu. "X-Parameters: Characterization, Modeling, and Design of Nonlinear RF and Microwave Components" The Cambridge RF and microwave engineering series, [Cambridge University Press (2013)](https://www.cambridge.org/sb/academic/subjects/engineering/rf-and-microwave-engineering/x-parameters-characterization-modeling-and-design-nonlinear-rf-and-microwave-components)

# Philosophy:

The motivation for developing this package is to simulate the gain and noise performance of ultra low noise amplifiers for quantum computing applications such as the [Josephson traveling-wave parametric amplifier](https://www.science.org/doi/10.1126/science.aaa8525), which have thousands of linear and nonlinear circuit elements. 

We prioritize speed (including compile time and time to first use), simplicity, and scalability over style. This motivates, for example, circuit input as an array rather than a DSL (domain-specific language).

# Future developments:

* Flux pumping, DC biasing, three wave mixing, fluxoid offsets, and two-tone harmonic balance.
* Design optimization with analytic derivatives or autodiff.
* Defining new nonlinear components by their branch flux response. 
* Iterative solvers for larger scale problems and shooting method for more nonlinear behavior.
* Interoperability with [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl) and [HarmonicBalance.jl](https://github.com/NonlinearOscillations/HarmonicBalance.jl) to perform bifurcation analysis and find unstable operating points.
* Interoperability with [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl) for time domain simulations. 


# Related packages and software:
* [Xyce.jl](https://github.com/JuliaComputing/Xyce.jl) provides a wrapper for [Xyce](https://xyce.sandia.gov/), the open source parallel circuit simulator from Sandia National Laboratories which can perform time domain and harmonic balance method simulations.
* [NgSpice.jl](https://github.com/JuliaComputing/Ngspice.jl) and [LTspice.jl](https://github.com/cstook/LTspice.jl) provide wrappers for [NgSpice](http://ngspice.sourceforge.net/) and [LTspice](https://www.analog.com/en/design-center/design-tools-and-calculators/ltspice-simulator.html), respectively.  
* [ModelingToolkit.jl](https://github.com/SciML/ModelingToolkit.jl) supports time domain circuit simulations from [scratch](https://mtk.sciml.ai/stable/tutorials/acausal_components) and using their [standard library](https://mtkstdlib.sciml.ai/dev/tutorials/rc_circuit)
* [ACME.jl](https://github.com/HSU-ANT/ACME.jl) simulates electrical circuits in the time domain with an emphasis on audio effect circuits.
* [Cedar EDA](https://cedar-eda.com) is a Julia-based commercial cloud service for circuit simulations.
* [Keysight ADS](https://www.keysight.com/us/en/products/software/pathwave-design-software/pathwave-advanced-design-system.html), [Cadence AWR](https://www.awr.com/), [Cadence Spectre RF](https://www.cadence.com/en_US/home/tools/custom-ic-analog-rf-design/circuit-simulation/spectre-rf-option.html), and [Qucs](http://qucs.sourceforge.net/) are capable of time and frequency domain analysis of nonlinear circuits. [WRSPICE](http://wrcad.com/wrspice.html) performs time domain simulations of Josephson junction containing circuits and frequency domain simulations of linear circuits. 

# Funding
We gratefully acknowledge funding from the [AWS Center for Quantum Computing](https://aws.amazon.com/blogs/quantum-computing/announcing-the-opening-of-the-aws-center-for-quantum-computing/) and the [MIT Center for Quantum Engineering (CQE)](https://cqe.mit.edu/).
