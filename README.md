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

For detailed mathematical derivations and implementation details, see [docs/nl_implementation.md](docs/nl_implementation.md).

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

See [docs/nl_implementation.md](docs/nl_implementation.md) for complete implementation details.

# Acknowledgments

Original JosephsonCircuits.jl developed by Kevin O'Brien. Taylor expansion nonlinearity feature contributed by Maxime Malnou.

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

### JPA with Taylor Expansion Nonlinearities: JJ vs NL Comparison

The following example demonstrates the Taylor expansion nonlinearity feature by comparing a JPA implemented with Josephson junctions versus nonlinear inductors:

```julia
using JosephsonCircuits
using CairoMakie

# Circuit parameters
Lj = 1000.0e-12  # 1 nH
Cc = 100.0e-15   # 100 fF
Cj = 1000.0e-15  # 1 pF
R = 50.0

# JJ version
jj_circuit = [
    ("P1", "1", "0", 1),
    ("R1", "1", "0", R),
    ("C1", "1", "2", Cc),
    ("Lj1", "2", "0", Lj),
    ("C2", "2", "0", Cj)
]

# NL version - Taylor approximation of Josephson junction
nl_circuit = [
    ("P1", "1", "0", 1),
    ("R1", "1", "0", R),
    ("C1", "1", "2", Cc),
    ("NL1", "2", "0", "poly Lj, 0.0, 0.5"),  # sin(φ) ≈ φ - φ³/6
    ("C2", "2", "0", Cj)
]

circuitdefs = Dict(
    :Lj => Lj,
    :Cc => Cc,
    :Cj => Cj,
    :R => R
)

# Simulation parameters
ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = (2*pi*4.75001*1e9,)
pump_current = 0.00565e-6
sources = [(mode=(1,), port=1, current=pump_current)]
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)

# Run both simulations
sol_jj = hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics, 
                 jj_circuit, circuitdefs)
sol_nl = hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics, 
                 nl_circuit, circuitdefs)

# Extract results
freq_GHz = ws./(2*pi*1e9)
S11_jj = abs2.(sol_jj.linearized.S((0,), 1, (0,), 1, :))
S11_nl = abs2.(sol_nl.linearized.S((0,), 1, (0,), 1, :))

# Plot comparison
fig = Figure(size = (600, 400))
ax = Axis(fig[1, 1], 
    xlabel = "Frequency [GHz]",
    ylabel = "S11 Gain [dB]",
    title = "JPA: Josephson Junction vs Taylor Expansion"
)

lines!(ax, freq_GHz, 10*log10.(S11_jj), label="JJ", linewidth=2)
lines!(ax, freq_GHz, 10*log10.(S11_nl), label="NL (Taylor)", linewidth=2)
axislegend(ax)
```

![JPA Comparison](examples/jpa_comparison.png)