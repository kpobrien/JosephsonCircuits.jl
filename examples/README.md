# NL Element Examples

These examples demonstrate the Taylor expansion (NL) nonlinearity functionality in this fork of JosephsonCircuits.jl.

## Running the Examples

Each example automatically sets up its environment and installs required dependencies (including CairoMakie for plotting):

```julia
# From the Julia REPL, navigate to the examples folder and run:
include("01_simple_jpa_example.jl")

# Or run all examples:
include("run_all_examples.jl")
```

The examples will automatically:
- Activate the correct Julia environment
- Install necessary packages (JosephsonCircuits, CairoMakie)
- Run the simulations and generate comparison plots

## Expected Results

Each example compares Josephson Junction (JJ) implementations with Taylor expansion (NL) approximations:

- `01_simple_jpa_example.jl`: JPA with no DC bias - slight gain difference
- `02_jtwpa_example.jl`: 2048-junction TWPA - slight gain difference
- `03_snailpa_example.jl`: SNAIL with DC bias - NL requires 0.94× DC scaling
- `04_flux_jtwpa_example.jl`: Flux-pumped TWPA - some gain difference (DC sensitive)
- `05_floquet_jtwpa_loss_example.jl`: Gaussian-weighted TWPA with dielectric loss - slight gain difference

All examples generate comparison plots saved as PNG files.

## NL Element Syntax

For detailed information about the NL (nonlinear) element syntax and implementation, see the main README in the parent directory.

Quick example - approximating a Josephson junction with Taylor expansion:
```julia
# Josephson junction
("Lj1", "1", "0", "Lj")  

# Equivalent NL element (sin(φ) ≈ φ - φ³/6), corresponding to L(φ) ≈ L₀(φ + φ²/2)
("NL1", "1", "0", "poly Lj, 0.0, 0.5")
```