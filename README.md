# JosephsonCircuits.jl (Fork with Taylor Expansion Nonlinearities)

> **Note**: This is a fork of the original [JosephsonCircuits.jl](https://github.com/kpobrien/JosephsonCircuits.jl) by Kevin O'Brien with added support for Taylor expansion nonlinearities (NL elements).

## Fork Features

### ✨ What's New in This Fork

**Taylor Expansion Nonlinearities (NL elements)**: Support for nonlinear inductors with polynomial current-flux relationships, enabling simulation of DC-biased RF SQUID TWPAs, KTWPAs, and other nonlinear inductance-based devices.

### 📝 Implementation Summary

This fork extends JosephsonCircuits.jl to support Taylor expansion nonlinearities in addition to Josephson junctions. This enables modeling of nonlinear inductors of the form:

```
L(φ) = L₀(1 + c₁φ + c₂φ² + c₃φ³ + c₄φ⁴)
```
corresponding to a current-phase relation:
```
I(φ) = φ₀/L₀ (φ - c₁φ²/2 + (c₁² - c₂)φ³/3 - (c₁³ - 2 c₁c₂ + c₃)φ⁴/4 + (c₁⁴ - 3 c₁²c₂ + c₂² + 2c₁c₃ - c₄)φ⁵/5)
```

Where:
- `L₀` is the linear inductance = φ₀/Ic
- `c₁, c₂, c₃, c₄` are the Taylor expansion coefficients
- `φ` is the flux

### 🔧 Technical Implementation

#### New Component Type: NL (Nonlinear Inductor)
- **Syntax**: `("NL1", "node1", "node2", "poly, L0, c1, c2, c3, c4")`
- **Example**: `("NL1", "1", "2", "poly, 329e-12, 0.0, 0.5")` approximates a 329 pH Josephson junction
- Supports symbolic variables defined in `circuitdefs` dictionary

#### Code Infrastructure

**Parser Modifications** (`parseinput.jl`)
- Added `:NL` to `allowedcomponenttypes`
- Parser recognizes "poly" syntax with symbolic variables
- Created `parse_nl_value()` to parse L(φ).
- Created `convert_poly_to_taylor_coeffs` to convert L(φ) to I(φ) coefficients
- Added `PolyNL` struct to hold polynomial NL information
- Added `NonlinearElement` struct (in `hbsolve.jl`) to track different types (josephson, taylor)
- Implemented `identify_nonlinear_elements()` function

**Unified FFT Machinery** (`hbsolve.jl`)
- Both JJ and NL are passed to the FFT machinery when applying the nonlinearity
- Created `apply_nonlinearities!` for generating appropriate functions per element type
- All nonlinearities (Josephson and Taylor) now go through `applynl_mixed!` (defined in `fftutils.jl`)
- Modified `calcfj2!` to handle all nonlinear elements uniformly
- Updated `hbnlsolve` and `hblinsolve` to use `all_nl_branches` instead of just `Ljb`

### 📊 Example Comparisons

The `examples/` folder contains side-by-side comparisons between the original JJ-based examples and equivalent NL element implementations, where Josephson junctions are replaced by their Taylor expansion approximation (expanded to second order, "poly, L0, c1, c2"). These comparisons are also documented in the `julia_wrapper_examples.ipynb` notebook in the TWPA Design Package.

### ✅ Comprehensive Testing

The fork has been tested against the original implementation:

**Basic Functionality**
- ✅ Linear NL elements match regular inductors exactly
- ✅ Taylor approximation of sin(φ) matches JJ for fundamental frequency
- ✅ Example: `"poly, 329e-12, 0.0, 0.5"` approximates a 329 pH Josephson junction
- ✅ No numerical issues even at high currents (120% of Ic)

**JosephsonCircuits Examples Verification**
1. ✅ **JJ-JPA**: Functions identically in both Registered and Forked versions
2. ✅ **NL-JPA**: Taylor approximation gives a similar result than JJ-JPA (slightly more gain) for the same pump parameters
3. ✅ **JJ-JTWPA**: Functions similarly in both Registered and Forked versions
4. ✅ **NL-JTWPA**: Taylor approximation gives a similar result than JJ-JTWPA (slightly more gain) for the same pump parameters
5. ✅ **JJ-SNAIL-PA**: Functions similarly in both versions
6. ✅ **NL-SNAIL-PA**: Taylor approximation needs lower dc bias (~ 0.94 Id) to give a similar result than JJ-SNAIL-PA
7. ✅ **JJ-flux-JTWPA**: Functions similarly in both Registered and Forked versions
8. ✅ **NL-flux-JTWPA**: Taylor approximation needs slightly lower dc bias (0.985 Id) and lower pump amplitude (0.7 Ip) to give a similar result than JJ-flux-JTWPA
9. ✅ **JJ-Floquet-JTWPA with loss**: Functions similarly in both Registered and Forked versions
10. ✅ **NL-Floquet-JTWPA with loss**: Taylor approximation gives a similar result than NL-Floquet-JTWPA (slightly more gain) for the same pump parameters

**Small Signal Regime**
- ✅ `hbnlsolve`: L, JJ, and NL all give identical results in small current regime
- ✅ `hbnlsolve`: JJ and NL give consistent results in high current regime

### 🚀 Future Improvements

The code would benefit from:
- ⭕ True DC analysis (currently not supported)

### 📁 Key Modified Files

- `parseinput.jl`: Added `:NL` to allowed components, created parser for polynomial syntax
- `hbsolve.jl`: Main modifications for unified nonlinear solver
- `calcfj2!`: Core function handling all nonlinearities
- `identify_nonlinear_elements()`: New function to find and categorize nonlinear elements
- `apply_nonlinearities!`: Handles mixed nonlinearity types
- `applynl_mixed!`: FFT machinery for multiple element types

---

## Original JosephsonCircuits.jl Documentation

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

# Installation:

## Installing This Fork

To use this fork with the TWPA Design Package, it will be automatically installed when you first use the julia_wrapper module. No manual installation needed.

For manual installation in Julia:
```julia
using Pkg
Pkg.add(url="https://github.com/MaxMalnou/JosephsonCircuits.jl")
```

### Running the Examples

The examples folder includes a `Project.toml` that specifies all required dependencies (including CairoMakie for plotting). When you run any example, it will automatically activate this environment and install the necessary packages:

```julia
# The examples automatically handle this:
using Pkg
Pkg.activate("path/to/examples")
Pkg.instantiate()  # Installs CairoMakie and other dependencies
```

No manual package installation is needed for running the examples.

## Installing Original Version

To install the latest release of the original package, start Julia and enter:
```julia
using Pkg
Pkg.add("JosephsonCircuits")
```

To install the development version:
```julia
using Pkg
Pkg.add(name="JosephsonCircuits",rev="main")
```

To run the examples below, you will need to install Plots.jl:
```julia
Pkg.add("Plots")
```

If you get errors when running the examples, please try installing the latest version of Julia and updating to the latest version of JosephsonCircuits.jl by running:
```julia
Pkg.update()
```

# Usage Examples

## Using NL Elements (Fork Feature)

### Basic NL Element Definition
```julia
# Define a circuit with Taylor expansion nonlinearity
circuit = [
    ("P1", "1", "0", "1"),
    ("R1", "1", "0", "50"),
    ("NL1", "1", "2", "poly, 1e-9, 0.0, 0.5, 0.0, 0.1"),  # L0=1nH, c2=0.5, c4=0.1
    ("C1", "2", "0", "1e-15"),
    ("P2", "2", "0", "2"),
    ("R2", "2", "0", "50")
]
```

### Using Symbolic Variables
```julia
# Circuit with symbolic parameters
circuit = [
    ("NL1", "1", "2", "poly, L0val, c1val, c2val, c3val, c4val")
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

### Approximating a Josephson Junction with Taylor Expansion
```julia
# Josephson junction circuit
jj_circuit = [("B1", "1", "0", "100e-6")]  # 100 μA critical current

# Equivalent Taylor approximation (sin(φ) ≈ φ - φ³/6)
# For a JJ: L_J = Φ₀/(2π*Ic) = 329 pH for Ic = 1 mA
nl_circuit = [("NL1", "1", "0", "poly, 329e-12, 0.0, 0.5, 0.0, 0.0")]
```

## Complete Example: TWPA Simulation

```julia
using JosephsonCircuits
using Plots

# Generate a netlist for a TWPA (example)
netlist = generate_twpa_netlist()  # Your netlist generation function

# Parse circuit
circuit = JosephsonCircuits.parseinputfile(netlist)

# Setup sources
sources = [(mode=(1,), port=1, current=1e-6)]  # 1 μA pump

# Run harmonic balance
solution = hbsolve(circuit, sources, freq=8e9, Nharmonics=10)

# Extract S-parameters
S = JosephsonCircuits.sparams(solution)
```

# References

[1-6] See original JosephsonCircuits.jl documentation for references

# Fork Maintenance

This fork is maintained by Maxime Malnou for use with the [TWPA Design Package](https://github.com/MaxMalnou/twpa_design).

## Reporting Issues
- **Fork-specific issues** (NL elements, Taylor expansions): Open an issue in this repository
- **General JosephsonCircuits.jl issues**: Report to the [original repository](https://github.com/kpobrien/JosephsonCircuits.jl)

## Acknowledgments

Original JosephsonCircuits.jl developed by Kevin O'Brien. Fork extensions for Taylor expansion nonlinearities developed for TWPA simulations.