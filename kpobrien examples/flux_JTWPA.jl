using JosephsonCircuits
using CairoMakie

const magnetic_flux_quantum = 2.0678338484619295e-15
const reduced_magnetic_flux_quantum = magnetic_flux_quantum / (2*pi)

@variables Rport C Cj Lj Lpump Cpump kappa Lg Lsmall

# From Section V. POSSIBLE CIRCUIT DESIGN
cutoff_frequency = 46e9 # [Hz]
transmission_line_impedance = 50.0

capacitance = 1 / (2 * pi * cutoff_frequency * transmission_line_impedance) # equation (51) [H]
junction_inductance = transmission_line_impedance / (2 * pi * cutoff_frequency) # L = L', equation (52) [F]
critical_current = reduced_magnetic_flux_quantum / junction_inductance

critical_current_density = 3e6 # Typical value mentioned in paper [A/m^2]
jj_area = critical_current / critical_current_density # [m^2]
jj_cap_density = 50 * 1e-15 / (1e-6)^2 # Typical [F/m^2]
jj_capacitance = jj_cap_density * jj_area

nr_cells = 500
modulation_parameter = 0.06

coupling = 0.02 # = M / L', equation (8)
mutual_inductance = coupling * junction_inductance
# Add small linear inductance in dc-SQUID loop to couple pump line inductance with
linear_squid_loop_inductance = mutual_inductance^2 / junction_inductance

# in order to reduce critical current of dc squid to critical current of junction
optimal_dc_flux = magnetic_flux_quantum / 3
Idc = optimal_dc_flux / mutual_inductance
Ip = modulation_parameter * Idc

pump_line_inductance = 1.1 * junction_inductance

pump_frequency = 20e9 # [Hz]
frequency_detuning = range(-0.5, 1.5, 500)
signal_frequency = pump_frequency / 2 .* (frequency_detuning .+ 1)

circuit = Tuple{String,String,String,Num}[]
entry = (elem, n1, n2, value) -> push!(circuit, ("$(elem)$(n1)_$(n2)", "$n1", "$n2", value))

function build_circuit()
    node = 1

    node_p1 = node # P1: Start of transmission line
    entry("P", node_p1, 0, 1)
    entry("R", node_p1, 0, Rport)

    node_p3 = node+1 # P3: Start of pump line
    entry("P", node_p3, 0, 3)
    entry("R", node_p3, 0, Rport)

    for cell_index in 1:nr_cells
        if cell_index == 1
            entry("C", node, 0, C/2)
        else
            entry("C", node, 0, C)
        end
        entry("Lj_a", node, node+3, Lj)
        entry("Cj_a", node, node+3, Cj)
        entry("L", node, node+2, Lsmall)
        entry("Lj_b", node+2, node+3, Lj)
        entry("Cj_b", node+2, node+3, Cj)

        entry("L", node+1, node+4, Lpump)
        if cell_index == 1
            entry("C", node+1, 0, Cpump/2)
        else
            entry("C", node+1, 0, Cpump)
        end
        push!(circuit, ("K$(node)", "L$(node)_$(node+2)", "L$(node+1)_$(node+4)", kappa))

        node += 3
    end

    entry("C", node, 0, C/2)
    entry("P", node, 0, 2) # P2: End of transmission line
    entry("R", node, 0, Rport)

    entry("C", node+1, 0, Cpump/2)
    entry("P", node+1, 0, 4) # P4: End of pump line
    entry("R", node+1, 0, Rport)
    entry("L", node+1, 0, Lg)

end

build_circuit()

circuitdefs = Dict(
    kappa => 0.999,
    Lg => 20.0e-9, # inductance to ground, required for solver
    Rport => 50.0,
    C => capacitance,
    Lj => junction_inductance,
    Lpump => pump_line_inductance,
    Cpump => pump_line_inductance / transmission_line_impedance^2,
    Lsmall => linear_squid_loop_inductance,
    Cj => jj_capacitance,
)

ws = 2*pi*signal_frequency
wp = (2*pi*pump_frequency,)

# add the DC bias and pump to port 3
sourcespumpon = [(mode=(0,),port=3,current=Idc),(mode=(1,),port=3,current=Ip)]
Npumpharmonics = (8,)
Nmodulationharmonics = (4,)
@time sol = hbsolve(ws, wp, sourcespumpon, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs;
    dc = true, threewavemixing=true,fourwavemixing=true,
    switchofflinesearchtol=0.0,alphamin=1e-7,iterations=200) # enable dc and three wave mixing

# Create figure with CairoMakie
fig = Figure(resolution = (1200, 800))

# S-parameters plot
ax1 = Axis(fig[1, 1], 
    xlabel = "Signal Frequency (GHz)",
    ylabel = "dB",
    title = "Scattering Parameters",
    limits = (nothing, (-40, 30))
)

freqs_ghz = sol.linearized.w/(2*pi*1e9)

# S21
s21_db = 10*log10.(abs2.(sol.linearized.S(
    outputmode=(0,),
    outputport=2,
    inputmode=(0,),
    inputport=1,
    freqindex=:)
))
lines!(ax1, freqs_ghz, s21_db, label="S21")

# S12
s12_db = 10*log10.(abs2.(sol.linearized.S((0,),1,(0,),2,:)))
lines!(ax1, freqs_ghz, s12_db, label="S12")

# S11
s11_db = 10*log10.(abs2.(sol.linearized.S((0,),1,(0,),1,:)))
lines!(ax1, freqs_ghz, s11_db, label="S11")

# S22
s22_db = 10*log10.(abs2.(sol.linearized.S((0,),2,(0,),2,:)))
lines!(ax1, freqs_ghz, s22_db, label="S22")

axislegend(ax1, position = :rb)

# Quantum efficiency plot
ax2 = Axis(fig[1, 2],
    xlabel = "Signal Frequency (GHz)",
    ylabel = "QE/QE_ideal",
    title = "Quantum efficiency",
    limits = (nothing, (0, 1.05))
)

qe_ratio = sol.linearized.QE((0,),2,(0,),1,:)./sol.linearized.QEideal((0,),2,(0,),1,:)
lines!(ax2, freqs_ghz, qe_ratio)

# All idlers plot
ax3 = Axis(fig[2, 1],
    xlabel = "Signal Frequency (GHz)",
    ylabel = "dB",
    title = "All idlers",
    limits = (nothing, (-40, 30))
)

all_idlers_db = 10*log10.(abs2.(sol.linearized.S(:,2,(0,),1,:)'))
for i in 1:size(all_idlers_db, 2)
    lines!(ax3, freqs_ghz, all_idlers_db[:, i])
end

# Commutation relation error plot
ax4 = Axis(fig[2, 2],
    xlabel = "Signal Frequency (GHz)",
    ylabel = "Commutation relation error",
    title = "Commutation relation error"
)

comm_error = 1 .- sol.linearized.CM((0,),2,:)
lines!(ax4, freqs_ghz, comm_error)

# Save the figure
save("flux_TWPA_results.png", fig)

# Display the figure
fig