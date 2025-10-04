# examples/04_flux_jtwpa_example.jl

using Pkg
Pkg.activate(@__DIR__)  # Use the examples environment
Pkg.instantiate()  # Install dependencies if needed

using JosephsonCircuits
using CairoMakie

println("=== Flux JTWPA Test Example ===")

# Physical constants
const ϕ₀ = 2.0678338484619295e-15
const ϕ₀_red = ϕ₀ / (2π)

# Circuit parameters
fc = 46e9  # Cutoff frequency
Z0 = 50.0
C = 1 / (2π * fc * Z0)
Lj = Z0 / (2π * fc)
Ic = ϕ₀_red / Lj
nr_cells = 500
k = 0.02  # Coupling coefficient
M = k * Lj  # Mutual inductance

# Calculate junction capacitance
Jc_density = 3e6  # A/m²
jj_area = Ic / Jc_density
jj_cap_density = 50e-15 / (1e-6)^2
Cj = jj_cap_density * jj_area

# Define parameter values
kappa = 0.999
Lg = 20.0e-9
Rport = 50.0
C_half = C/2
Lpump = 1.1 * Lj
Cpump = 1.1 * Lj / Z0^2
Cpump_half = Cpump / 2
Lsmall = M^2 / Lj

circuitdefs = Dict(
    :kappa => kappa,
    :Lg => Lg,
    :Rport => Rport,
    :C => C,
    :C_half => C_half,
    :Lj => Lj,
    :Lpump => Lpump,
    :Cpump => Cpump,
    :Cpump_half => Cpump_half,
    :Lsmall => Lsmall,
    :Cj => Cj
)

# Build circuit function
function build_flux_jtwpa_circuit(use_nl::Bool)
    circuit = Tuple{String,String,String,Any}[]
    
    node = 1
    
    # Ports
    push!(circuit, ("P$node", "$node", "0", 1))
    push!(circuit, ("R$node", "$node", "0", Rport))
    
    node_p3 = node + 1
    push!(circuit, ("P$node_p3", "$node_p3", "0", 3))
    push!(circuit, ("R$node_p3", "$node_p3", "0", Rport))
    
    # Build cells
    for i in 1:nr_cells
        # Signal line
        cap_val = i == 1 ? C_half : C
        push!(circuit, ("C$(node)_0", "$node", "0", cap_val))
        
        # SQUID structure
        if use_nl
            push!(circuit, ("NL_a$(node)_$(node+3)", "$node", "$(node+3)", "poly Lj, 0, 0.5"))
            push!(circuit, ("NL_b$(node+2)_$(node+3)", "$(node+2)", "$(node+3)", "poly Lj, 0, 0.5"))
        else
            push!(circuit, ("Lj_a$(node)_$(node+3)", "$node", "$(node+3)", Lj))
            push!(circuit, ("Lj_b$(node+2)_$(node+3)", "$(node+2)", "$(node+3)", Lj))
        end
        
        push!(circuit, ("Cj_a$(node)_$(node+3)", "$node", "$(node+3)", Cj))
        push!(circuit, ("L$(node)_$(node+2)", "$node", "$(node+2)", Lsmall))
        push!(circuit, ("Cj_b$(node+2)_$(node+3)", "$(node+2)", "$(node+3)", Cj))
        
        # Pump line
        push!(circuit, ("L$(node+1)_$(node+4)", "$(node+1)", "$(node+4)", Lpump))
        cap_val = i == 1 ? Cpump_half : Cpump
        push!(circuit, ("C$(node+1)_0", "$(node+1)", "0", cap_val))
        
        # Mutual coupling
        push!(circuit, ("K$node", "L$(node)_$(node+2)", "L$(node+1)_$(node+4)", kappa))
        
        node += 3
    end
    
    # Final components
    push!(circuit, ("C$(node)_0", "$node", "0", C_half))
    push!(circuit, ("P$node", "$node", "0", 2))
    push!(circuit, ("R$node", "$node", "0", Rport))
    
    push!(circuit, ("C$(node+1)_0", "$(node+1)", "0", Cpump_half))
    push!(circuit, ("P$(node+1)", "$(node+1)", "0", 4))
    push!(circuit, ("R$(node+1)", "$(node+1)", "0", Rport))
    push!(circuit, ("L$(node+1)_0", "$(node+1)", "0", Lg))
    
    return circuit
end

# Build circuits
jj_circuit = build_flux_jtwpa_circuit(false)
nl_circuit = build_flux_jtwpa_circuit(true)

# Simulation parameters
ws = 2*pi*(5:0.05:19.5)*1e9
# ws = 2*pi*(6:0.1:12)*1e9
wp = (2*pi*20e9,)

# DC flux bias
dc_flux = ϕ₀ / 3  # Φ₀/3
dc_current_jj = dc_flux / M
dc_current_nl = dc_current_jj * 0.985
# pump amplitude
modulation_jj = 0.06
modulation_nl = 0.06 * 0.7

sources_jj = [
    (mode=(0,), port=3, current=dc_current_jj),
    (mode=(1,), port=3, current=modulation_jj * dc_current_jj)
]

sources_nl = [
    (mode=(0,), port=3, current=dc_current_nl),
    (mode=(1,), port=3, current=modulation_nl * dc_current_nl)
]

Npumpharmonics = (8,)
Nmodulationharmonics = (4,)


# Run simulations
println("Running JJ version ($(length(jj_circuit)) components)...")
sol_jj = hbsolve(ws, wp, sources_jj, Nmodulationharmonics, Npumpharmonics,
                 jj_circuit, circuitdefs, dc=true, threewavemixing=true,
                 fourwavemixing=true, sorting=:name,
                 switchofflinesearchtol=0.0, alphamin=1e-7, iterations=200)



println("Running NL version...")
sol_nl = hbsolve(ws, wp, sources_nl, Nmodulationharmonics, Npumpharmonics,
                 nl_circuit, circuitdefs, dc=true, threewavemixing=true,
                 fourwavemixing=true, sorting=:name,
                 switchofflinesearchtol=0.0, alphamin=1e-7, iterations=1000, ftol=1e-5)

# Extract and plot
freq_GHz = ws./(2*pi*1e9)
S21_jj = abs2.(sol_jj.linearized.S(outputmode=(0,), outputport=2,
                                   inputmode=(0,), inputport=1, freqindex=:))
S21_nl = abs2.(sol_nl.linearized.S(outputmode=(0,), outputport=2,
                                   inputmode=(0,), inputport=1, freqindex=:))

fig = Figure(size = (600, 400))

ax = Axis(fig[1, 1],
    xlabel = "Frequency [GHz]",
    ylabel = "S21 Gain [dB]",
    title = "Flux-pumped JTWPA (Φ₀/3 bias)"
)

lines!(ax, freq_GHz, 10*log10.(S21_jj), label="JJ", linewidth=2, color=:blue)
lines!(ax, freq_GHz, 10*log10.(S21_nl), label="NL", linewidth=2, color=:red)
axislegend(ax)

save(joinpath(@__DIR__, "flux_jtwpa_comparison.png"), fig)
display(fig)

println("\n=== Performance Summary ===")
println("DC flux bias: Φ₀/3")
println("JJ DC current: $(round(dc_current_jj*1e6, digits=1)) μA")
println("NL DC current: $(round(dc_current_nl*1e6, digits=1)) μA")
println("JJ pump current: $(round(modulation_jj * dc_current_jj*1e6, digits=1)) μA")
println("NL pump current: $(round(modulation_nl * dc_current_nl*1e6, digits=1)) μA")
println("JJ max gain: $(round(maximum(10*log10.(S21_jj)), digits=1)) dB")
println("NL max gain: $(round(maximum(10*log10.(S21_nl)), digits=1)) dB")