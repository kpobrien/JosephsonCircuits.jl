# Taylor Expansion Nonlinearity Implementation

## Overview

This document details the implementation of Taylor expansion nonlinearities (NL elements) in JosephsonCircuits.jl, enabling simulation of DC-biased RF SQUID TWPAs, KTWPAs, and other nonlinear inductance-based devices.

## Mathematical Foundation

### Inductance Model

The NL element models nonlinear inductors with a polynomial inductance-flux relationship:

```
L(φ) = L₀(1 + c₁φ + c₂φ² + c₃φ³ + c₄φ⁴)
```

Where:
- `L₀` is the linear inductance
- `c₁, c₂, c₃, c₄` are the Taylor expansion coefficients
- `φ` is the flux across the inductor

### Current-Phase Relationship

The corresponding current-phase relation is derived from `I = ∫(φ/L(φ))dφ`:

```
I(φ) = φ₀/L₀ (φ - c₁φ²/2 + (c₁² - c₂)φ³/3 - (c₁³ - 2c₁c₂ + c₃)φ⁴/4 + (c₁⁴ - 3c₁²c₂ + c₂² + 2c₁c₃ - c₄)φ⁵/5)
```

### Conversion Formulas

The implementation converts polynomial inductance coefficients to Taylor series current coefficients:

- **1st order**: `I₁ = 1.0` (normalized by φ₀/L₀)
- **2nd order**: `I₂ = -c₁/2`
- **3rd order**: `I₃ = (c₁² - c₂)/3`
- **4th order**: `I₄ = (-c₃ + 2c₁c₂ - c₁³)/4`
- **5th order**: `I₅ = (c₁⁴ - 3c₁²c₂ + c₂² + 2c₁c₃ - c₄)/5`

## Component Syntax

### Basic Usage
```julia
("NL1", "node1", "node2", "poly L0, c1, c2, c3, c4")
```

### Examples
```julia
# Fixed coefficients
("NL1", "1", "2", "poly 329e-12, 0.0, 0.5")  # Approximates 1μA JJ

# Symbolic parameters
("NL1", "1", "2", "poly L0val, c1val, c2val")
# With circuitdefs = Dict("L0val" => 329e-12, "c1val" => 0.0, "c2val" => 0.5)
```

## Implementation Details

### Modified Files and Functions

#### 1. `src/parseinput.jl` - Parser Extensions

**New Data Structures**

```julia
struct PolyNL
    L0                    # Base inductance (can be symbolic)
    c_coeffs::Vector     # Polynomial coefficients [c1, c2, c3, c4]
    value::Float64       # Numeric value for arithmetic operations
end
```
- Represents polynomial inductance L(φ)
- Includes arithmetic operator overloading for compatibility
- Provides conversion methods to Float64 and ComplexF64

**Modified Functions**

`parsecircuit()` - Extended to recognize NL components:
```julia
# Added "NL" to recognized component types
allowedcomponents = ["Lj","NL","L","C","K","I","R","P"]

# Special handling for NL components with "poly" syntax
if componenttype == :NL && isa(value, AbstractString) && startswith(value, "poly")
    componentvalues[i] = parse_nl_value(value)
end
```

`extractbranches!()` - Include NL in branch-capable components:
```julia
allowedcomponenttypes = [:Lj,:L,:I,:P,:V,:NL]
```

**New Functions**

`parse_nl_value(value_str::AbstractString)`
- Parses "poly L0, c1, c2, c3, c4" format
- Supports numeric and symbolic coefficients
- Returns PolyNL struct

`convert_poly_to_taylor_coeffs(poly::PolyNL)`
- Converts polynomial L(φ) to Taylor series I(φ)
- Implements mathematical conversion formulas
- Returns Dict with :coeffs, :powers, :inductance

`identify_nonlinear_elements(componenttypes, componentvalues, nodeindices, edge2indexdict, circuitdefs)`
- Creates unified nonlinearity tracking for JJ and NL elements
- Maps components to branch indices
- Returns Dict{Int, NonlinearElement}

#### 2. `src/hbsolve.jl` - Solver Integration

**New Data Structures**

```julia
struct NonlinearElement
    type::Symbol  # :josephson or :taylor
    indices::Vector{Int}  # Branch indices
    params::Dict{Symbol,Any}  # Parameters for the nonlinearity
end
```

**Modified Functions**

`hbsolve()` - Main entry point:
- Added `x0 = nothing` parameter for initial guess
- Calls `identify_nonlinear_elements()` to catalog all nonlinearities
- Passes `nonlinear_elements` through to solver

`hbnlsolve()` - Core nonlinear solver:
- Creates sparse vector including ALL nonlinear elements (not just Josephson)
- Recalculates Lmean based on all nonlinear inductances
- Uses unified `all_nl_branches` instead of just `Ljb`

`calcfj2!()` - Jacobian and function evaluation:
- Extracts all nonlinear branches (Josephson and Taylor)
- Applies nonlinearities through unified FFT machinery
- Handles inductance-specific scaling for each type

**New Functions**

`apply_nonlinearities!(phimatrix, phimatrixtd, nonlinear_elements, mode, irfftplan, rfftplan)`
- Applies appropriate nonlinearity based on element type
- Modes: `:function` (sin/polynomial) or `:jacobian` (cos/derivative)
- Creates function array for mixed nonlinearity evaluation

#### 3. `src/fftutils.jl` - FFT Extensions

**New Functions**

`applynl_mixed!(fd, td, nl_functions, irfftplan, rfftplan)`
- Extension of original `applynl!()` for column-specific functions
- Each column gets its own nonlinear function
- Enables mixed Josephson/Taylor circuits

Process:
1. Transform to time domain via IFFT
2. Apply column-specific nonlinearities
3. Transform back to frequency domain via FFT
4. Apply normalization

`applynl_mixed_verbose!(fd, td, nl_functions, irfftplan, rfftplan)`
- Debug version with extensive logging
- Tracks transformations per column
- Uses flag to prevent multiple debug outputs

#### 4. `src/capindmat.jl` - Matrix Calculations

**Modified Functions**

`numericmatrices()` - Convert NL Dict values to numeric:
```julia
for (i, type) in enumerate(psc.componenttypes)
    if type == :NL && isa(vvn[i], Dict)
        coeffs = vvn[i][:coeffs]
        vvn[i] = phi0 / coeffs[1]  # L0 value
    end
end
```

`calcLjb()` - Include NL elements:
```julia
# Modified to include both Josephson and NL elements
return calcbranchvector(..., [:Lj, :NL], ..., [:Lj, :NL], ...)
```

`calcbranchvector()` - Major refactor:
- Accepts multiple component types
- Handles component filtering to prevent double-counting
- Preserves original data structure (Dict for NL)

`calcLmean()` - Include NL in mean calculation:
- Forces Float64 type when NL present
- Includes NL inductances in Lmean normalization

#### 5. `src/JosephsonCircuits.jl` - Module Exports

**New Exports**
```julia
export NonlinearElement
export identify_nonlinear_elements
export debug_log, get_debug_log, clear_debug_log
export warning_log, get_warning_log, clear_warning_log
```

## Design Decisions

### 1. Unified FFT Machinery
Both Josephson and Taylor nonlinearities use the same FFT-based evaluation pipeline, maintaining Kevin O'Brien's efficient architecture while extending it for heterogeneous circuits.

### 2. Backward Compatibility
All original functionality is preserved. Circuits with only Josephson junctions work identically to the original implementation.

### 3. Normalization Consistency
All inductances (Josephson and Taylor) are normalized by Lmean, preserving the numerical stability of the original solver.

### 4. Sparse Matrix Preservation
The implementation maintains sparsity patterns even with mixed nonlinearities, ensuring scalability to large circuits.

## Testing and Validation

### Test Coverage
- Linear NL elements match regular inductors exactly
- Taylor approximation of sin(φ) matches JJ behavior for small amplitudes
- Mixed circuits with both JJ and NL elements converge properly
- Symbolic parameter evaluation works correctly

### Example Validations
The `examples/` folder contains side-by-side comparisons demonstrating that:
- NL approximations of Josephson junctions produce similar gain curves
- Parameter adjustments (bias current, pump power) can compensate for approximation differences
- Both element types work in complex circuits (SNAIL-PA, flux-pumped TWPA, Floquet TWPA)

## Usage Examples

### Simple NL Circuit
```julia
circuit = [
    ("NL1", "1", "0", "poly 329e-12, 0.0, 0.5"),  # Nonlinear inductor
    ("C1", "1", "0", "1e-15"),  # Shunt capacitor
    ("P1", "1", "0", "1")  # Port
]
```

### Mixed JJ/NL Circuit
```julia
circuit = [
    ("B1", "1", "2", "1e-6"),  # Josephson junction
    ("NL1", "2", "3", "poly 500e-12, 0.0, 0.3"),  # Nonlinear inductor
    ("C1", "3", "0", "2e-15")  # Capacitor
]
```

### Symbolic Parameters
```julia
circuit = [("NL1", "1", "0", "poly Lj*exp(alpha), 0, beta")]
circuitdefs = Dict("Lj" => 300e-12, "alpha" => 0.1, "beta" => 0.3)
```

## Future Enhancements

Potential improvements include:
- True DC analysis support
- Higher-order polynomial terms (beyond 4th order)
- Automatic Taylor expansion order optimization
- Performance optimizations for large mixed circuits