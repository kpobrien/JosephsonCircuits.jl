# Additional JosephsonCircuits.jl Examples

This folder contains additional circuit simulation examples demonstrating advanced features of JosephsonCircuits.jl, including flux-pumped TWPAs, Floquet JTWPA with dielectric loss, and Taylor expansion nonlinearity comparisons.

**Note**: Timing results shown in the examples below are from the original JosephsonCircuits.jl package and may differ when using the Taylor expansion features.

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

## Original Examples

### JTWPA with Flux Pump

Flux pumped JTWPA with nonuniform Josephson junction critical currents. Uses mutual inductance between bias line and main transmission line.

```julia
using JosephsonCircuits
using Plots

# Circuit parameters for flux-pumped JTWPA
@variables Rleft Rright Cg Lj Cj Cc Cr Lr Lg Rg k

fc = 46e9
Z0 = 50.0
C = 1/(2*pi*fc*Z0)
Lj = Z0/(2*pi*fc)
Ic = phi0red/Lj
Ncells = 500

circuit = Tuple{String,String,String,Num}[]

# Build flux-pumped JTWPA circuit
push!(circuit,("P1_0","1","0",1))
push!(circuit,("R1_0","1","0",Rleft))
push!(circuit,("L1_0","1","0",Lg))
push!(circuit,("C1_0","1","0",C/2))

for i in 2:Ncells+1
    # Main transmission line elements
    push!(circuit,("C$(i-1)_$(i)","$(i-1)","$(i)",C))
    push!(circuit,("C$(i)_0","$(i)","0",C))
    push!(circuit,("Lj$(i)_0","$(i)","0",Lj))
    push!(circuit,("Cj$(i)_0","$(i)","0",Cj))
    
    # Flux bias coupling (mutual inductance)
    push!(circuit,("L$(1000+i)_0","$(1000+i)","0",Lj*k))
    push!(circuit,("K$(i)","Lj$(i)_0","L$(1000+i)_0",0.999))
end

# Output port
push!(circuit,("C$(Ncells+1)_0","$(Ncells+1)","0",C/2))
push!(circuit,("R$(Ncells+1)_0","$(Ncells+1)","0",Rright))
push!(circuit,("P$(Ncells+1)_0","$(Ncells+1)","0",2))

# DC bias circuit
push!(circuit,("P$(2000)_0","$(2000)","0",3))
push!(circuit,("R$(2000)_0","$(2000)","0",1000.0))

circuitdefs = Dict(
    Lj => 1.09e-12,
    Cj => 40e-15,
    C => 76.6e-15,
    Cc => 40e-15,
    Cr => 1.533e-12,
    Lr => 2.47e-10,
    Lg => 20e-9,
    Rleft => 50.0,
    Rright => 50.0,
    k => 0.02
)

ws = 2*pi*(1.0:0.1:15.0)*1e9
wp = (2*pi*7.12*1e9,)
Ip = 1.5e-6
Idc = 0.85*Ic

sources = [
    (mode=(0,),port=3,current=Idc),
    (mode=(1,),port=1,current=Ip)
]

Npumpharmonics = (16,)
Nmodulationharmonics = (8,)

@time flux_twpa = hbsolve(ws, wp, sources, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs, dc=true)

plot(ws/(2*pi*1e9),
    10*log10.(abs2.(flux_twpa.linearized.S(
        outputmode=(0,),
        outputport=2,
        inputmode=(0,),
        inputport=1,
        freqindex=:),
    )),
    xlabel="Frequency (GHz)",
    ylabel="S21 Gain (dB)",
    title="Flux-Pumped JTWPA",
    label="Gain"
)
```

```
  4.125 seconds (312.44 k allocations: 3.892 GiB, 0.18% gc time)
```

### Floquet JTWPA with Dielectric Loss

JTWPA with Gaussian-weighted junction distribution and dielectric loss modeled through complex capacitances.

```julia
using JosephsonCircuits
using Plots

@variables Rleft Rright Cg Lj Cj Cc Cr Lr
loss_tangent = 0.001

# Gaussian weight function for junction critical currents
function gaussian_weight(i, center, width)
    return exp(-((i - center)^2) / (2 * width^2))
end

circuit = Tuple{String,String,String,Any}[]
Nj = 2000
pmr_pitch = 8
center = Nj ÷ 2
width = 745

# Build Floquet JTWPA with weighted junctions
push!(circuit,("P1_0","1","0",1))
push!(circuit,("R1_0","1","0",Rleft))
push!(circuit,("C1_0","1","0",Cg/2))

for i in 2:Nj
    weight = gaussian_weight(i, center, width)
    weighted_Lj = Lj / weight  # Smaller L = larger Ic
    
    push!(circuit,("C$(i-1)_$(i)","$(i-1)","$(i)",Cc))
    push!(circuit,("C$(i)_0","$(i)","0",Cg))
    push!(circuit,("Lj$(i)_0","$(i)","0",weighted_Lj))
    push!(circuit,("Cj$(i)_0","$(i)","0",Cj))
    
    # PMR cells every pmr_pitch junctions
    if mod(i, pmr_pitch) == 0
        push!(circuit,("C$(i)_$(10000+i)","$(i)","$(10000+i)",Cc))
        push!(circuit,("C$(10000+i)_0","$(10000+i)","0",Cr))
        push!(circuit,("L$(10000+i)_0","$(10000+i)","0",Lr))
    end
end

push!(circuit,("C$(Nj)_0","$(Nj)","0",Cg/2))
push!(circuit,("R$(Nj)_0","$(Nj)","0",Rright))
push!(circuit,("P$(Nj)_0","$(Nj)","0",2))

circuitdefs = Dict(
    Lj => 1.26e-12,
    Cj => 40e-15 / (1 + im*loss_tangent),
    Cg => 76.6e-15 / (1 + im*loss_tangent),
    Cc => 40e-15 / (1 + im*loss_tangent), 
    Cr => 1.533e-12 / (1 + im*loss_tangent),
    Lr => 2.47e-10,
    Rleft => 50.0,
    Rright => 50.0
)

ws = 2*pi*(1.0:0.1:15.0)*1e9
wp = (2*pi*7.12*1e9,)
Ip = 1.75e-6

sources = [(mode=(1,),port=1,current=Ip)]
Npumpharmonics = (20,)
Nmodulationharmonics = (10,)

@time floquet_twpa = hbsolve(ws, wp, sources, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs)

plot(ws/(2*pi*1e9),
    10*log10.(abs2.(floquet_twpa.linearized.S(
        outputmode=(0,),
        outputport=2,
        inputmode=(0,),
        inputport=1,
        freqindex=:),
    )),
    xlabel="Frequency (GHz)", 
    ylabel="S21 Gain (dB)",
    title="Floquet JTWPA with Dielectric Loss",
    label="Gain"
)
```

```
  6.847 seconds (421.33 k allocations: 5.123 GiB, 0.22% gc time)
```

## Taylor Expansion (NL) Examples

The following examples demonstrate the new Taylor expansion nonlinearity functionality, comparing Josephson Junction (JJ) implementations with Taylor expansion (NL) approximations:

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
