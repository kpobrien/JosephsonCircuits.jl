# examples/02_jtwpa_example.jl

using Pkg
Pkg.activate(@__DIR__)  # Use the examples environment
Pkg.instantiate()  # Install dependencies if needed

using JosephsonCircuits
using CairoMakie

println("=== JTWPA Test Example ===")

# Physical constants
const ϕ₀ = 2.0678338484619295e-15
const ϕ₀_red = ϕ₀ / (2π)

# Circuit parameters
Ic = 3.4e-6  # Critical current
Lj = ϕ₀_red / Ic  # Junction inductance
nr_junctions = 2048
pmr_pitch = 4

# Define all parameter values
Cg = 45.0e-15
Cg_half = 22.5e-15
Cg_minus_Cc = 15.0e-15
Cc = 30.0e-15
Cr = 2.8153e-12
Lr = 1.70e-10
Cj = 55e-15
Rleft = 50.0
Rright = 50.0

circuitdefs = Dict(
    :Lj => Lj,
    :Cg => Cg,
    :Cg_half => Cg_half,
    :Cg_minus_Cc => Cg_minus_Cc,
    :Cc => Cc,
    :Cr => Cr,
    :Lr => Lr,
    :Cj => Cj,
    :Rleft => Rleft,
    :Rright => Rright
)

# Build circuit function
function build_jtwpa_circuit(use_nl::Bool)
    circuit = Tuple{String,String,String,Any}[]
    
    # Port 1
    push!(circuit, ("P1_0", "1", "0", 1))
    push!(circuit, ("R1_0", "1", "0", Rleft))
    push!(circuit, ("C1_0", "1", "0", Cg_half))
    
    # First junction
    if use_nl
        push!(circuit, ("NL1_2", "1", "2", "poly Lj, 0, 0.5"))
    else
        push!(circuit, ("Lj1_2", "1", "2", Lj))
    end
    push!(circuit, ("C1_2", "1", "2", Cj))
    
    # Build middle cells
    node = 2
    for i in 2:nr_junctions-1
        if i % pmr_pitch == pmr_pitch ÷ 2
            # PMR cell
            push!(circuit, ("C$(node)_0", "$node", "0", Cg_minus_Cc))
            if use_nl
                push!(circuit, ("NL$(node)_$(node+2)", "$node", "$(node+2)", "poly Lj, 0, 0.5"))
            else
                push!(circuit, ("Lj$(node)_$(node+2)", "$node", "$(node+2)", Lj))
            end
            push!(circuit, ("C$(node)_$(node+2)", "$node", "$(node+2)", Cj))
            
            # PMR branch
            push!(circuit, ("C$(node)_$(node+1)", "$node", "$(node+1)", Cc))
            push!(circuit, ("C$(node+1)_0", "$(node+1)", "0", Cr))
            push!(circuit, ("L$(node+1)_0", "$(node+1)", "0", Lr))
            node += 2
        else
            # Regular cell
            push!(circuit, ("C$(node)_0", "$node", "0", Cg))
            if use_nl
                push!(circuit, ("NL$(node)_$(node+1)", "$node", "$(node+1)", "poly Lj, 0, 0.5"))
            else
                push!(circuit, ("Lj$(node)_$(node+1)", "$node", "$(node+1)", Lj))
            end
            push!(circuit, ("C$(node)_$(node+1)", "$node", "$(node+1)", Cj))
            node += 1
        end
    end
    
    # Last components
    push!(circuit, ("C$(node)_0", "$node", "0", Cg_half))
    push!(circuit, ("R$(node)_0", "$node", "0", Rright))
    push!(circuit, ("P$(node)_0", "$node", "0", 2))
    
    return circuit
end

# Build both circuits
jj_circuit = build_jtwpa_circuit(false)
nl_circuit = build_jtwpa_circuit(true)

# Simulation parameters
ws = 2*pi*(1.0:0.1:14.0)*1e9
wp = (2*pi*7.12*1e9,)
pump_current_jj = 1.85e-6
pump_current_nl = 1.85e-6
sources_jj = [(mode=(1,), port=1, current=pump_current_jj)]
sources_nl = [(mode=(1,), port=1, current=pump_current_nl)]
Npumpharmonics = (20,)
Nmodulationharmonics = (10,)

# Run simulations

println("Running JJ version ($(length(jj_circuit)) components)...")
sol_jj = hbsolve(ws, wp, sources_jj, Nmodulationharmonics, Npumpharmonics,
                 jj_circuit, circuitdefs, sorting=:name)
                  

println("Running NL version ($(length(nl_circuit)) components)...")
sol_nl = hbsolve(ws, wp, sources_nl, Nmodulationharmonics, Npumpharmonics,
                 nl_circuit, circuitdefs, sorting=:name)

# Extract results
freq_GHz = ws./(2*pi*1e9)
S21_jj = abs2.(sol_jj.linearized.S(outputmode=(0,), outputport=2,
                                   inputmode=(0,), inputport=1, freqindex=:))
S21_nl = abs2.(sol_nl.linearized.S(outputmode=(0,), outputport=2,
                                   inputmode=(0,), inputport=1, freqindex=:))

QE_jj = sol_jj.linearized.QE((0,),2,(0,),1,:) ./
        sol_jj.linearized.QEideal((0,),2,(0,),1,:)
QE_nl = sol_nl.linearized.QE((0,),2,(0,),1,:) ./
        sol_nl.linearized.QEideal((0,),2,(0,),1,:)

# Create figure
fig = Figure(size = (600, 600))

ax1 = Axis(fig[1, 1],
    xlabel = "Frequency [GHz]",
    ylabel = "S21 Gain [dB]",
    title = "JTWPA Gain Comparison"
)

lines!(ax1, freq_GHz, 10*log10.(S21_jj), label="JJ", linewidth=2, color=:blue)
lines!(ax1, freq_GHz, 10*log10.(S21_nl), label="NL", linewidth=2, color=:red)
vlines!(ax1, [7.12], label="pump", linestyle=:dot, color=:purple, alpha=0.5)
axislegend(ax1)

ax2 = Axis(fig[2, 1],
    xlabel = "Frequency [GHz]",
    ylabel = "QE/QE_ideal",
    title = "Quantum Efficiency"
)

lines!(ax2, freq_GHz, QE_jj, label="JJ", linewidth=2, color=:blue)
lines!(ax2, freq_GHz, QE_nl, label="NL", linewidth=2, color=:red)
hlines!(ax2, [1.0], label="Ideal", linestyle=:dash, color=:black, alpha=0.7)
ylims!(ax2, 0, 1.05)
axislegend(ax2)

save(joinpath(@__DIR__, "jtwpa_comparison.png"), fig)
display(fig)

println("\n=== Performance Summary ===")
println("JJ max S21 gain: $(round(maximum(10*log10.(S21_jj)), digits=1)) dB")
println("NL max S21 gain: $(round(maximum(10*log10.(S21_nl)), digits=1)) dB")