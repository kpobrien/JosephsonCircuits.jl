# examples/05_floquet_jtwpa_loss_example.jl

using Pkg
Pkg.activate(@__DIR__)  # Use the examples environment
Pkg.instantiate()  # Install dependencies if needed

using JosephsonCircuits
using CairoMakie

println("=== Floquet JTWPA with Dielectric Loss ===")

# Physical constants
const ϕ₀ = 2.0678338484619295e-15
const ϕ₀_red = ϕ₀ / (2π)

# Circuit parameters
Ic = 1.75e-6
Lj = ϕ₀_red / Ic
loss_tangent = 0.001

# Base parameters (with loss applied to capacitors)
Rleft = 50.0
Rright = 50.0
Cg = 76.6e-15 / (1 + im*loss_tangent)
Cc = 40.0e-15 / (1 + im*loss_tangent)
Cr = 1.533e-12 / (1 + im*loss_tangent)
Lr = 2.47e-10
Cj = 40e-15 / (1 + im*loss_tangent)
c1 = 0.0
c2 = 0.5

circuitdefs = Dict(
    :Rleft => Rleft,
    :Rright => Rright,
    :Lj => Lj,
    :Cg => Cg,
    :Cc => Cc,
    :Cr => Cr,
    :Lr => Lr,
    :Cj => Cj,
    :c1 => c1,
    :c2 => c2
)

# Build parameters
nr_junctions = 2000
pmr_pitch = 8
weight_width = 745

# Weight function - returns numeric value
weight(n, N, w) = exp(-(n - N/2)^2/(w)^2)

# Build circuit function
function build_floquet_circuit(use_nl::Bool)
    circuit = Tuple{String,String,String,Any}[]
    
    # Port on the left
    push!(circuit, ("P1_0", "1", "0", 1))
    push!(circuit, ("R1_0", "1", "0", Rleft))
    
    # First half cap - compute weighted value
    push!(circuit, ("C1_0", "1", "0", Cg/2 * weight(0.5, nr_junctions, weight_width)))
    
    # First junction
    if use_nl
        push!(circuit, ("NL1_2", "1", "2", "poly Lj*$(weight(1, nr_junctions, weight_width)), c1, c2"))
    else
        push!(circuit, ("Lj1_2", "1", "2", Lj * weight(1, nr_junctions, weight_width)))
    end
    push!(circuit, ("C1_2", "1", "2", Cj / weight(1, nr_junctions, weight_width)))
    
    # Build middle cells
    j = 2
    for i in 2:nr_junctions-1
        w_val = weight(i-0.5, nr_junctions, weight_width)
        w_junc = weight(i, nr_junctions, weight_width)
        
        if i % pmr_pitch == pmr_pitch ÷ 2
            # PMR cell
            push!(circuit, ("C$(j)_0", "$j", "0", (Cg-Cc) * w_val))
            
            if use_nl
                push!(circuit, ("NL$(j)_$(j+2)", "$j", "$(j+2)", 
                               "poly Lj*$(w_junc), c1, c2"))
            else
                push!(circuit, ("Lj$(j)_$(j+2)", "$j", "$(j+2)", Lj * w_junc))
            end
            
            push!(circuit, ("C$(j)_$(j+2)", "$j", "$(j+2)", Cj / w_junc))
            
            # PMR components
            push!(circuit, ("C$(j)_$(j+1)", "$j", "$(j+1)", Cc * w_val))
            push!(circuit, ("C$(j+1)_0", "$(j+1)", "0", Cr))
            push!(circuit, ("L$(j+1)_0", "$(j+1)", "0", Lr))
            
            j += 1
        else
            # Regular cell
            push!(circuit, ("C$(j)_0", "$j", "0", Cg * w_val))
            
            if use_nl
                push!(circuit, ("NL$(j)_$(j+1)", "$j", "$(j+1)", 
                               "poly Lj*$(w_junc), c1, c2"))
            else
                push!(circuit, ("Lj$(j)_$(j+1)", "$j", "$(j+1)", Lj * w_junc))
            end
            
            push!(circuit, ("C$(j)_$(j+1)", "$j", "$(j+1)", Cj / w_junc))
        end
        j += 1
    end
    
    # Last components
    w_last = weight(nr_junctions-0.5, nr_junctions, weight_width)
    push!(circuit, ("C$(j)_0", "$j", "0", Cg/2 * w_last))
    push!(circuit, ("R$(j)_0", "$j", "0", Rright))
    push!(circuit, ("P$(j)_0", "$j", "0", 2))
    
    return circuit
end

# Build circuits
jj_circuit = build_floquet_circuit(false)
nl_circuit = build_floquet_circuit(true)

# Simulation parameters
ws = 2*pi*(1.0:0.1:14.0)*1e9
wp = (2*pi*7.9*1e9,)

# Pump current with loss compensation
pump_base_jj = 1.1e-6
pump_base_nl = 1.1e-6

pump_adjusted_jj = pump_base_jj * (1 + 125 * loss_tangent)
pump_adjusted_nl = pump_base_nl * (1 + 125 * loss_tangent)

sources_jj = [(mode=(1,), port=1, current=pump_adjusted_jj)]
sources_nl = [(mode=(1,), port=1, current=pump_adjusted_nl)]

Npumpharmonics = (20,)
Nmodulationharmonics = (10,)

# Run simulations
println("Running JJ version with loss ($(length(jj_circuit)) components)...")
sol_jj = hbsolve(ws, wp, sources_jj, Nmodulationharmonics, Npumpharmonics,
                 jj_circuit, circuitdefs, sorting=:name)


println("Running NL version with loss...")
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
    title = "Floquet JTWPA with Loss (tan δ = $loss_tangent)"
)

lines!(ax1, freq_GHz, 10*log10.(S21_jj), label="JJ", linewidth=2, color=:blue)
lines!(ax1, freq_GHz, 10*log10.(S21_nl), label="NL", linewidth=2, color=:red)
vlines!(ax1, [7.9], label="pump", linestyle=:dot, color=:purple, alpha=0.5)
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

save(joinpath(@__DIR__, "floquet_jtwpa_loss_comparison.png"), fig)
display(fig)

println("\n=== Performance Summary ===")
println("Loss tangent: $loss_tangent")
println("JJ base pump: $(round(pump_base_jj*1e6, digits=2)) μA")
println("JJ adjusted pump: $(round(pump_adjusted_jj*1e6, digits=2)) μA")
println("NL base pump: $(round(pump_base_nl*1e6, digits=2)) μA")
println("NL adjusted pump: $(round(pump_adjusted_nl*1e6, digits=2)) μA")
println("JJ max gain: $(round(maximum(10*log10.(S21_jj)), digits=1)) dB")
println("NL max gain: $(round(maximum(10*log10.(S21_nl)), digits=1)) dB")