# JosephsonCircuits.jl

[![Code coverage](https://codecov.io/gh/kpobrien/JosephsonCircuits.jl/branch/main/graphs/badge.svg)](https://codecov.io/gh/kpobrien/JosephsonCircuits.jl)
[![Build Status](https://github.com/kpobrien/JosephsonCircuits.jl/actions/workflows/CI.yml/badge.svg
)](https://github.com/kpobrien/JosephsonCircuits.jl/actions?query=workflow) [![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/J/JosephsonCircuits.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/J/JosephsonCircuits.html) [![Stable docs](https://img.shields.io/badge/docs-stable-blue.svg)](https://josephsoncircuits.org/stable)
 [![Dev docs](https://img.shields.io/badge/docs-dev-blue.svg)](https://josephsoncircuits.org/dev)

[JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl) is a high-performance frequency domain simulator for nonlinear circuits containing Josephson junctions, capacitors, inductors, mutual inductors, and resistors. [JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl) simulates the frequency domain behavior using a variant [1] of nodal analysis [2] and the harmonic balance method [3-5] with an analytic Jacobian. Noise performance, quantified by quantum efficiency, is efficiently simulated through an adjoint method.

Frequency dependent circuit parameters are supported to model realistic impedance environments or dissipative components. Dissipation can be modeled by capacitors with an imaginary capacitance or frequency dependent resistors. 

[JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl) supports the following:
* Nonlinear simulations in which the user defines a circuit, the drive current, frequency, and number of harmonics and the code calculates the node flux or node voltage at each harmonic.
* Linearized simulations about the nonlinear operating point calculated above. This simulates the small signal response of a periodically time varying linear circuit and is useful for simulating parametric amplification and frequency conversion in the undepleted (strong) pump limit. Calculation of node fluxes (or node voltages) and scattering parameters of the linearized circuit [4-5].
* Linear simulations of linear circuits. Calculation of node fluxes (or node voltages) and scattering parameters.
* Calculation of symbolic capacitance and inverse inductance matrices.

As detailed in [6], we find excellent agreement with [Keysight ADS](https://www.keysight.com/us/en/products/software/pathwave-design-software/pathwave-advanced-design-system.html) simulations and Fourier analysis of time domain simulation performed by [WRSPICE](http://wrcad.com/wrspice.html).

**Warning:** this package is under heavy development and there will be breaking changes. We will keep the examples updated to ease the burden of any breaking changes.

## New Feature: Taylor Expansion Nonlinearities

JosephsonCircuits.jl now supports Taylor expansion nonlinearities (NL elements) in addition to Josephson junctions. This enables modeling of nonlinear inductors with polynomial current-flux relationships, useful for simulating DC-biased RF SQUID TWPAs, KTWPAs, and other nonlinear inductance-based devices.

### Mathematical Model

The NL element models nonlinear inductors of the form:
```
L(φ) = L₀(1 + c₁φ + c₂φ² + c₃φ³ + c₄φ⁴)
```

Where:
- `L₀` is the linear inductance
- `c₁, c₂, c₃, c₄` are the Taylor expansion coefficients  
- `φ` is the flux

### Current-Phase Relationship

The corresponding current-phase relation is derived from `φ₀ dφ/dt = L di/dt` where `φ₀` is the reduced flux quantum:

```
I(φ) = φ₀/L₀ (φ - c₁φ²/2 + (c₁² - c₂)φ³/3 - (c₁³ - 2c₁c₂ + c₃)φ⁴/4 + (c₁⁴ - 3c₁²c₂ + c₂² + 2c₁c₃ - c₄)φ⁵/5)
```

### Usage

#### Basic NL Element Definition
```julia
# Define a circuit with Taylor expansion nonlinearity
circuit = [
    ("P1", "1", "0", "1"),
    ("R1", "1", "0", "50"),
    ("NL1", "1", "2", "poly 1e-9, 0.0, 0.5, 0.0, 0.1"),  # L0=1nH, c2=0.5, c4=0.1
    ("C1", "2", "0", "1e-15"),
    ("P2", "2", "0", "2"),
    ("R2", "2", "0", "50")
]
```

#### Using Symbolic Variables
```julia
# Circuit with symbolic parameters
circuit = [
    ("NL1", "1", "2", "poly L0val, c1val, c2val, c3val, c4val")
]

# Define parameters in dictionary
circuitdefs = Dict(
    "L0val" => 1e-9,    # Base inductance
    "c1val" => 0.0,     # Linear term (usually 0)
    "c2val" => 0.5,     # Quadratic term
    "c3val" => 0.0,     # Cubic term
    "c4val" => 0.1      # Quartic term
)
```

#### Approximating a Josephson Junction with Taylor Expansion
```julia
# Josephson junction circuit
jj_circuit = [("B1", "1", "0", "1e-6")]  # 1 μA critical current

# Equivalent Taylor approximation (sin(φ) ≈ φ - φ³/6)
# For a JJ: L_J = `\phi_0`/(2π*Ic) = 329 pH for Ic = 1 μA
nl_circuit = [("NL1", "1", "0", "poly 329e-12, 0.0, 0.5")]
```

### Technical Summary

- **Component Type**: New `:NL` component type for nonlinear inductors
- **Syntax**: `"poly L0[, c1][, c2][, c3][, c4]"` format with support for symbolic parameters (coefficients are optional, default to 0)
- **Integration**: Extends existing harmonic balance solver through unified FFT machinery
- **Mixed Circuits**: Supports circuits with both Josephson junctions and Taylor expansion elements

For detailed implementation details, see [docs/nl_implementation.md](docs/nl_implementation.md).

# Acknowledgments

Original JosephsonCircuits.jl developed by [Kevin O'Brien](https://github.com/kpobrien). Taylor expansion nonlinearity feature contributed by [Maxime Malnou](https://github.com/MaxMalnou).

# Installation:

To install the latest release of the package, install Julia using [Juliaup](https://github.com/JuliaLang/juliaup), start Julia, and enter the following command:
```
using Pkg
Pkg.add("JosephsonCircuits")
```

To install the development version, start Julia and enter the command:
```
using Pkg
Pkg.add(name="JosephsonCircuits",rev="main")
```

To run the examples below, you will need to install Plots.jl using the command:
```
Pkg.add("Plots")
```

To run the Taylor expansion nonlinearity comparison examples, you will also need to install CairoMakie.jl:
```
Pkg.add("CairoMakie")
```

If you get errors when running the examples, please try installing the latest version of Julia and updating to the latest version of JosephsonCircuits.jl by running:
```
Pkg.update()
```

Then check that you are running the latest version of the package with:
```
Pkg.status()
```

Simulations of the linearized system can be effectively parallelized, so we suggest starting Julia with the number of threads equal to the number of physical cores. This can be done with the command line argument `--threads` or by setting the environmental variable `JULIA_NUM_THREADS`. See the [Julia documentation](https://docs.julialang.org/en/v1/manual/multi-threading) for the more details. Verify you are using the desired number of threads by running:
```
Threads.nthreads()
```
For context, the simulation times reported for the examples below use 16 threads on an AMD Ryzen 9 9950X system running Linux.

The examples can be run in the command line (REPL) after starting Julia or you can run them in a Jupyter notebook with [IJulia](https://github.com/JuliaLang/IJulia.jl) or in Visual Studio Code with the [Julia extension](https://code.visualstudio.com/docs/languages/julia).

# Usage:
Generate a netlist using circuit components including capacitors `C`, inductors `L`, Josephson junctions described by the Josephson inductance `Lj`, nonlinear inductors described by Taylor expansion coefficients `NL`, mutual inductors described by the mutual coupling coefficient `K`, and resistors `R`. See the [SPICE netlist format](https://duckduckgo.com/?q=spice+netlist+format), docstrings, and examples below for usage. Run the harmonic balance analysis using [`hbnlsolve`](https://josephsoncircuits.org/stable/reference/#JosephsonCircuits.hbnlsolve-Union{Tuple{K},%20Tuple{N},%20Tuple{NTuple{N,%20Number},%20Any,%20JosephsonCircuits.Frequencies{N},%20JosephsonCircuits.FourierIndices{N},%20JosephsonCircuits.ParsedSortedCircuit,%20JosephsonCircuits.CircuitGraph,%20JosephsonCircuits.CircuitMatrices}}%20where%20{N,%20K}) to solve a nonlinear system at one operating point, [`hblinsolve`](https://josephsoncircuits.org/dev/reference/#JosephsonCircuits.hblinsolve-Union{Tuple{K},%20Tuple{Any,%20Any,%20Any}}%20where%20K) to solve a linear (or linearized) system at one or more frequencies, or [`hbsolve`](https://josephsoncircuits.org/dev/reference/#JosephsonCircuits.hbsolve-Union{Tuple{K},%20Tuple{M},%20Tuple{N},%20Tuple{Any,%20NTuple{N,%20Number},%20Vector,%20NTuple{M,%20Int64},%20NTuple{N,%20Int64},%20Any,%20Any}}%20where%20{N,%20M,%20K}) to run both analyses. Add a question mark `?` in front of a function to access the docstring. For example, type (don't copy-paste) the following to see the documentation for `hbsolve`:
```
?hbsolve
```

# Examples:
## Josephson parametric amplifier (JPA)
A driven nonlinear LC resonator.

**Note**: Timing results shown below are from the original JosephsonCircuits.jl package and may differ when using the Taylor expansion features.

```julia
using JosephsonCircuits
using Plots

@variables R Cc Lj Cj
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("C1","1","2",Cc),
    ("Lj1","2","0",Lj),
    ("C2","2","0",Cj)]

circuitdefs = Dict(
    Lj =>1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    R => 50.0)

ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = (2*pi*4.75001*1e9,)
Ip = 0.00565e-6
sources = [(mode=(1,),port=1,current=Ip)]
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)

@time jpa = hbsolve(ws, wp, sources, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs)

plot(
    jpa.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpa.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    label="JosephsonCircuits.jl",
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
)
```

```
  0.001817 seconds (12.99 k allocations: 4.361 MiB)
```

![JPA simulation with JosephsonCircuits.jl](https://qce.mit.edu/JosephsonCircuits.jl/jpa.png)


Compare with WRspice. Please note that on Linux you can install the [XicTools_jll](https://github.com/JuliaBinaryWrappers/XicTools_jll.jl/) package which provides WRspice for x86_64. For other operating systems and platforms, you can install WRspice yourself and substitute `XicTools_jll.wrspice()` with `JosephsonCircuits.wrspice_cmd()` which will attempt to provide the path to your WRspice executable. 

```julia
using XicTools_jll

wswrspice=2*pi*(4.5:0.01:5.0)*1e9
n = JosephsonCircuits.exportnetlist(circuit,circuitdefs);
input = JosephsonCircuits.wrspice_input_paramp(n.netlist,wswrspice,wp[1],2*Ip,(0,1),(0,1));

@time output = JosephsonCircuits.spice_run(input,XicTools_jll.wrspice());
S11,S21=JosephsonCircuits.wrspice_calcS_paramp(output,wswrspice,n.Nnodes);

plot!(wswrspice/(2*pi*1e9),10*log10.(abs2.(S11)),
    label="WRspice",
    seriestype=:scatter)

```

```
 12.743245 seconds (32.66 k allocations: 499.263 MiB, 0.41% gc time)
```

![JPA simulation with JosephsonCircuits.jl and WRspice](https://qce.mit.edu/JosephsonCircuits.jl/jpa_WRspice.png)

### JPA with Taylor Expansion Nonlinearities: JJ vs NL Comparison

The following example demonstrates the Taylor expansion nonlinearity feature by comparing a JPA implemented with Josephson junctions versus nonlinear inductors:

```julia
using JosephsonCircuits
using CairoMakie

println("=== Simple JPA Test (No DC Bias) ===")

# Define circuit parameter values as variables
Lj = 1000.0e-12  # 1 nH
Cc = 100.0e-15   # 100 fF
Cj = 1000.0e-15  # 1 pF
R = 50.0

# Circuit components for JJ version - use actual variables
jj_circuit = Tuple{String,String,String,Any}[
    ("P1", "1", "0", 1),
    ("R1", "1", "0", R),     # Use the actual variable value
    ("C1", "1", "2", Cc),    # Use the actual variable value
    ("Lj1", "2", "0", Lj),   # Use the actual variable value
    ("C2", "2", "0", Cj)     # Use the actual variable value
]

# Circuit components for NL version
nl_circuit = Tuple{String,String,String,Any}[
    ("P1", "1", "0", 1),
    ("R1", "1", "0", R),
    ("C1", "1", "2", Cc),
    ("NL1", "2", "0", "poly Lj, 0.0, 0.5"),  # This stays as string
    ("C2", "2", "0", Cj)
]

# Circuit parameters dictionary with symbols
circuitdefs = Dict(
    :Lj => Lj,
    :Cc => Cc,
    :Cj => Cj,
    :R => R
)

# Simulation parameters
ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = (2*pi*4.75001*1e9,)
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)

# Sources
pump_current_jj = 0.00565e-6
pump_current_nl = 0.00565e-6

sources_jj = [(mode=(1,), port=1, current=pump_current_jj)]
sources_nl = [(mode=(1,), port=1, current=pump_current_nl)]

# Run simulations
println("Running JJ version...")
println("JJ pump current: ", pump_current_jj*1e9, " nA")
sol_jj = hbsolve(ws, wp, sources_jj, Nmodulationharmonics, Npumpharmonics, 
                 jj_circuit, circuitdefs, sorting=:name)

println("Running NL version...")
println("NL pump current: ", pump_current_nl*1e9, " nA")
sol_nl = hbsolve(ws, wp, sources_nl, Nmodulationharmonics, Npumpharmonics, 
                 nl_circuit, circuitdefs, sorting=:name)

# Extract results
freq_GHz = ws./(2*pi*1e9)
S11_jj = abs2.(sol_jj.linearized.S(outputmode=(0,), outputport=1, 
                                   inputmode=(0,), inputport=1, freqindex=:))
S11_nl = abs2.(sol_nl.linearized.S(outputmode=(0,), outputport=1, 
                                   inputmode=(0,), inputport=1, freqindex=:))
QE_jj = sol_jj.linearized.QE((0,),1,(0,),1,:) ./ 
        sol_jj.linearized.QEideal((0,),1,(0,),1,:)
QE_nl = sol_nl.linearized.QE((0,),1,(0,),1,:) ./ 
        sol_nl.linearized.QEideal((0,),1,(0,),1,:)

# Create figure with CairoMakie
fig = Figure(size = (600, 600))

# Gain plot
ax1 = Axis(fig[1, 1], 
    xlabel = "Frequency [GHz]",
    ylabel = "S11 Gain [dB]",
    title = "JPA Gain Comparison"
)

lines!(ax1, freq_GHz, 10*log10.(S11_jj), label="JJ", linewidth=2, color=:blue)
lines!(ax1, freq_GHz, 10*log10.(S11_nl), label="NL", linewidth=2, color=:red)
vlines!(ax1, [4.75001], label="pump", linestyle=:dot, color=:purple, alpha=0.5)
axislegend(ax1, position = :lt)

# QE plot
ax2 = Axis(fig[2, 1],
    xlabel = "Frequency [GHz]",
    ylabel = "QE/QE_ideal",
    title = "Quantum Efficiency"
)

lines!(ax2, freq_GHz, QE_jj, label="JJ", linewidth=2, color=:blue)
lines!(ax2, freq_GHz, QE_nl, label="NL", linewidth=2, color=:red)
hlines!(ax2, [1.0], label="Ideal", linestyle=:dash, color=:black, alpha=0.7)
ylims!(ax2, 0, 1.05)
axislegend(ax2, position = :lb)

# Print summary
println("\n=== Performance Summary ===")
println("JJ max gain: $(round(maximum(10*log10.(S11_jj)), digits=1)) dB")
println("NL max gain: $(round(maximum(10*log10.(S11_nl)), digits=1)) dB")
println("Max QE (JJ): $(round(maximum(QE_jj), digits=3))")
println("Max QE (NL): $(round(maximum(QE_nl), digits=3))")
```

Output:
```julia
=== Simple JPA Test (No DC Bias) ===
Running JJ version...
JJ pump current: 5.65 nA
Running NL version...
NL pump current: 5.65 nA

=== Performance Summary ===
JJ max gain: 13.3 dB
NL max gain: 14.3 dB
Max QE (JJ): 1.0
Max QE (NL): 1.0
```

<img src="examples/jpa_comparison.png" width="60%">

## SNAIL Parametric Amplifier
Circuit parameters from [here](https://doi.org/10.1103/PhysRevApplied.10.054020). Notice that the resonance frequency is similar for pump-on and pump-off, indicating it is operating near the Kerr-free point.

<details>

<summary>Code</summary>

```julia
using JosephsonCircuits
using Plots

@variables R Cc Cj Lj Cr Lr Ll Ldc K Lg
alpha = 0.29
Z0 = 50
w0 = 2*pi*8e9
l=10e-3
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    # a very large inductor so the DC node flux of this node isn't floating
    ("L0","1","0",Lg), 
    ("C1","1","2",Cc),
    ("L1","2","3",Lr),
    ("C2","2","0",Cr),
    ("Lj1","3","0",Lj/alpha),
    ("Cj1","3","0",Cj/alpha),
    ("L2","3","4",Ll),
    ("Lj2","4","5",Lj),
    ("Cj2","4","5",Cj),
    ("Lj3","5","6",Lj),
    ("Cj3","5","6",Cj),
    ("Lj4","6","0",Lj),
    ("Cj4","6","0",Cj),
    ("L3","7","0",Ldc), 
    ("K1","L2","L3",K),
    # a port with a very large resistor so we can apply the bias across the port
    ("P2","7","0",2),
    ("R2","7","0",1000.0),
] 

circuitdefs = Dict(
    Lj => 60e-12,
    Cj => 10.0e-15, 
    Lr =>0.4264e-9*1.25,
    Cr => 0.4e-12*1.25,
    Lg => 100.0e-9,
    Cc => 0.048e-12,
    R => 50.0, 
    Ll => 34e-12, 
    K => 0.999, # the inverse inductance matrix for K=1.0 diverges, so set K<1.0
    Ldc => 0.74e-12,
)

# ws = 2*pi*(9.7:0.0001:9.8)*1e9
# ws = 2*pi*(5.0:0.001:11)*1e9
ws = 2*pi*(7.8:0.001:8.2)*1e9
wp = (2*pi*16.00*1e9,)
Ip = 4.4e-6
# Idc = 140.3e-6
Idc = 0.000159
# add the DC bias and pump to port 2
sourcespumpon = [(mode=(0,),port=2,current=Idc),(mode=(1,),port=2,current=Ip)]
sourcespumpoff = [(mode=(0,),port=2,current=Idc),(mode=(1,),port=2,current=0.0)]
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)
@time jpapumpon = hbsolve(ws, wp, sourcespumpon, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs, dc = true, threewavemixing=true,fourwavemixing=true) # enable dc and three wave mixing
@time jpapumpoff = hbsolve(ws, wp, sourcespumpoff, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs, dc = true, threewavemixing=true,fourwavemixing=true) # enable dc and three wave mixing

p1 = plot(
    jpapumpon.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpapumpon.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
    label="pump on",
)

plot!(
    jpapumpoff.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpapumpoff.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    label="pump off",
)

p2 = plot(
    jpapumpon.linearized.w/(2*pi*1e9),
    angle.(
        jpapumpon.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    ),
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
    label="pump on",
)

plot!(
    jpapumpoff.linearized.w/(2*pi*1e9),
    angle.(
        jpapumpoff.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    ),
    label="pump off",
)
plot(p1,p2,layout=(2,1))
```

</details>


```
  0.010345 seconds (16.74 k allocations: 40.025 MiB)
  0.011252 seconds (16.68 k allocations: 39.985 MiB)
```

![SNAIL parametric amplifier simulation with JosephsonCircuits.jl](https://qce.mit.edu/JosephsonCircuits.jl/snail.png)


and compare with WRspice
<details>

<summary>Code</summary>

```julia
using XicTools_jll

# simulate the JPA in WRSPICE
wswrspice=2*pi*(7.8:0.005:8.2)*1e9
n = JosephsonCircuits.exportnetlist(circuit,circuitdefs);
input = JosephsonCircuits.wrspice_input_paramp(n.netlist,wswrspice,[0.0,wp[1]],[Idc,2*Ip],[(0,1)],[(0,7),(0,7)];trise=10e-9,tstop=600e-9);

@time output = JosephsonCircuits.spice_run(input,XicTools_jll.wrspice());
S11,S21=JosephsonCircuits.wrspice_calcS_paramp(output,wswrspice,n.Nnodes);

# plot the output
plot(
    jpapumpon.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpapumpon.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
    label="JosephsonCircuits.jl",
)

plot!(wswrspice/(2*pi*1e9),10*log10.(abs2.(S11)),
    label="WRspice",
    seriestype=:scatter)
```

</details>

```
2067.364975 seconds (149.73 k allocations: 29.873 GiB, 0.01% gc time)
```

![SNAIL parametric amplifier simulation with JosephsonCircuits.jl and WRspice](https://qce.mit.edu/JosephsonCircuits.jl/snail_WRspice.png)

### SNAIL Parametric Amplifier with Taylor Expansion Nonlinearities: JJ vs NL Comparison

The following example demonstrates the Taylor expansion implementation for a SNAIL parametric amplifier. Note that the NL version requires 94% of the original JJ DC bias current to achieve similar performance:

<details>

<summary>Code</summary>

```julia
using JosephsonCircuits
using CairoMakie

println("=== SNAILPA Test Example ===")

# Circuit parameters
α = 0.29
Lj = 60e-12
Lj_large = 60e-12 / α
Cj = 10.0e-15
Cj_large = 10.0e-15 / α
Lr = 0.4264e-9 * 1.25
Cr = 0.4e-12 * 1.25
Lg = 100.0e-9
Cc = 0.048e-12
R = 50.0
Ll = 34e-12
K = 0.999
Ldc = 0.74e-12
Rdc = 1000.0

circuitdefs = Dict(
    :Lj => Lj,
    :Lj_large => Lj_large,
    :Cj => Cj,
    :Cj_large => Cj_large,
    :Lr => Lr,
    :Cr => Cr,
    :Lg => Lg,
    :Cc => Cc,
    :R => R,
    :Ll => Ll,
    :K => K,
    :Ldc => Ldc,
    :Rdc => Rdc
)

# JJ circuit
jj_circuit = Tuple{String,String,String,Any}[
    ("P1", "1", "0", 1),
    ("R1", "1", "0", R),
    ("L0", "1", "0", Lg),
    ("C1", "1", "2", Cc),
    ("L1", "2", "3", Lr),
    ("C2", "2", "0", Cr),
    ("Lj1", "3", "0", Lj_large),
    ("Cj1", "3", "0", Cj_large),
    ("L2", "3", "4", Ll),
    ("Lj2", "4", "5", Lj),
    ("Cj2", "4", "5", Cj),
    ("Lj3", "5", "6", Lj),
    ("Cj3", "5", "6", Cj),
    ("Lj4", "6", "0", Lj),
    ("Cj4", "6", "0", Cj),
    ("L3", "7", "0", Ldc),
    ("K1", "L2", "L3", K),
    ("P2", "7", "0", 2),
    ("R2", "7", "0", Rdc)
]

# NL circuit
nl_circuit = Tuple{String,String,String,Any}[
    ("P1", "1", "0", 1),
    ("R1", "1", "0", R),
    ("L0", "1", "0", Lg),
    ("C1", "1", "2", Cc),
    ("L1", "2", "3", Lr),
    ("C2", "2", "0", Cr),
    ("NL1", "3", "0", "poly Lj_large, 0.0, 0.5"),
    ("Cj1", "3", "0", Cj_large),
    ("L2", "3", "4", Ll),
    ("NL2", "4", "5", "poly Lj, 0.0, 0.5"),
    ("Cj2", "4", "5", Cj),
    ("NL3", "5", "6", "poly Lj, 0.0, 0.5"),
    ("Cj3", "5", "6", Cj),
    ("NL4", "6", "0", "poly Lj, 0.0, 0.5"),
    ("Cj4", "6", "0", Cj),
    ("L3", "7", "0", Ldc),
    ("K1", "L2", "L3", K),
    ("P2", "7", "0", 2),
    ("R2", "7", "0", Rdc)
]

# Simulation parameters
ws = 2*pi*(7.8:0.001:8.2)*1e9
wp = (2*pi*16.0*1e9,)
dc_current_jj = 159e-6
dc_current_nl = 159e-6 * 0.94
pump_current_jj = 4.4e-6
pump_current_nl = 4.4e-6

sources_jj = [
    (mode=(0,), port=2, current=dc_current_jj),
    (mode=(1,), port=2, current=pump_current_jj)
]

sources_nl = [
    (mode=(0,), port=2, current=dc_current_nl),
    (mode=(1,), port=2, current=pump_current_nl)
]

Npumpharmonics = (16,)
Nmodulationharmonics = (8,)


# Run simulations
println("Running JJ version with DC bias...")
sol_jj = hbsolve(ws, wp, sources_jj, Nmodulationharmonics, Npumpharmonics,
                 jj_circuit, circuitdefs, dc=true, threewavemixing=true,
                 fourwavemixing=true, sorting=:name)



println("Running NL version with scaled DC bias...")
sol_nl = hbsolve(ws, wp, sources_nl, Nmodulationharmonics, Npumpharmonics,
                 nl_circuit, circuitdefs, dc=true, threewavemixing=true,
                 fourwavemixing=true, sorting=:name)

# Extract results
freq_GHz = ws./(2*pi*1e9)
S11_jj = abs2.(sol_jj.linearized.S(outputmode=(0,), outputport=1,
                                   inputmode=(0,), inputport=1, freqindex=:))
S11_nl = abs2.(sol_nl.linearized.S(outputmode=(0,), outputport=1,
                                   inputmode=(0,), inputport=1, freqindex=:))

# Create figure
fig = Figure(size = (600, 400))

ax = Axis(fig[1, 1],
    xlabel = "Frequency [GHz]",
    ylabel = "S11 Gain [dB]",
    title = "SNAILPA with DC Bias"
)

lines!(ax, freq_GHz, 10*log10.(S11_jj), label="JJ", linewidth=2, color=:blue)
lines!(ax, freq_GHz, 10*log10.(S11_nl), label="NL", linewidth=2, color=:red)
axislegend(ax)

println("\n=== Performance Summary ===")
println("JJ max gain: $(round(maximum(10*log10.(S11_jj)), digits=1)) dB")
println("NL max gain: $(round(maximum(10*log10.(S11_nl)), digits=1)) dB")
println("JJ DC bias: $(round(dc_current_jj*1e6, digits=1)) μA")
println("NL DC bias: $(round(dc_current_nl*1e6, digits=1)) μA")
println("JJ pump current: $(round(pump_current_jj*1e6, digits=1)) μA")
println("NL pump current: $(round(pump_current_nl*1e6, digits=1)) μA")
```

</details>

Output:
```Julia
=== SNAILPA Test Example ===
Running JJ version with DC bias...
Running NL version with scaled DC bias...

=== Performance Summary ===
JJ max gain: 21.1 dB
NL max gain: 19.6 dB
JJ DC bias: 159.0 μA
NL DC bias: 149.5 μA
JJ pump current: 4.4 μA
NL pump current: 4.4 μA
```

<img src="examples/snailpa_comparison.png" width="60%">
