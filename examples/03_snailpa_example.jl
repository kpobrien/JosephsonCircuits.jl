# examples/03_snailpa_example.jl

using Pkg
Pkg.activate(@__DIR__)  # Use the examples environment
Pkg.instantiate()  # Install dependencies if needed

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

save(joinpath(@__DIR__, "snailpa_comparison.png"), fig)
display(fig)

println("\n=== Performance Summary ===")
println("JJ max gain: $(round(maximum(10*log10.(S11_jj)), digits=1)) dB")
println("NL max gain: $(round(maximum(10*log10.(S11_nl)), digits=1)) dB")
println("JJ DC bias: $(round(dc_current_jj*1e6, digits=1)) μA")
println("NL DC bias: $(round(dc_current_nl*1e6, digits=1)) μA")
println("JJ pump current: $(round(pump_current_jj*1e6, digits=1)) μA")
println("NL pump current: $(round(pump_current_nl*1e6, digits=1)) μA")