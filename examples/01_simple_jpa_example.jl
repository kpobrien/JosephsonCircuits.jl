# examples/01_simple_jpa_example.jl
using Pkg
Pkg.activate(@__DIR__)  # Use the examples environment
Pkg.instantiate()  # Install dependencies if needed

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

save(joinpath(@__DIR__, "jpa_comparison.png"), fig)
display(fig)

# Print summary
println("\n=== Performance Summary ===")
println("JJ max gain: $(round(maximum(10*log10.(S11_jj)), digits=1)) dB")
println("NL max gain: $(round(maximum(10*log10.(S11_nl)), digits=1)) dB")
println("Max QE (JJ): $(round(maximum(QE_jj), digits=3))")
println("Max QE (NL): $(round(maximum(QE_nl), digits=3))")