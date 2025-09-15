using JosephsonCircuits
using CairoMakie

@variables R Cc Cj Lj Cr Lr Ll Ldc K Lg
alpha = 0.29
Z0 = 50
w0 = 2*pi*8e9
l=10e-3
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    # a very large inductor so the DC node flux of this node isn't floating
    ("L0","1","0",Lg), 
    ("C1","1","2",Cc),
    ("L1","2","3",Lr),
    ("C2","2","0",Cr),
    ("Lj1","3","0",Lj/alpha),
    ("Cj1","3","0",Cj/alpha),
    ("L2","3","4",Ll),
    ("Lj2","4","5",Lj),
    ("Cj2","4","5",Cj),
    ("Lj3","5","6",Lj),
    ("Cj3","5","6",Cj),
    ("Lj4","6","0",Lj),
    ("Cj4","6","0",Cj),
    ("L3","7","0",Ldc), 
    ("K1","L2","L3",K),
    # a port with a very large resistor so we can apply the bias across the port
    ("P2","7","0",2),
    ("R2","7","0",1000.0),
] 

circuitdefs = Dict(
    Lj => 60e-12,
    Cj => 10.0e-15, 
    Lr =>0.4264e-9*1.25,
    Cr => 0.4e-12*1.25,
    Lg => 100.0e-9,
    Cc => 0.048e-12,
    R => 50.0, 
    Ll => 34e-12, 
    K => 0.999, # the inverse inductance matrix for K=1.0 diverges, so set K<1.0
    Ldc => 0.74e-12,
)

# ws = 2*pi*(9.7:0.0001:9.8)*1e9
# ws = 2*pi*(5.0:0.001:11)*1e9
ws = 2*pi*(7.8:0.001:8.2)*1e9
wp = (2*pi*16.00*1e9,)
Ip = 4.4e-6
# Idc = 140.3e-6
Idc = 0.000159
# add the DC bias and pump to port 2
sourcespumpon = [(mode=(0,),port=2,current=Idc),(mode=(1,),port=2,current=Ip)]
sourcespumpoff = [(mode=(0,),port=2,current=Idc),(mode=(1,),port=2,current=0.0)]
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)
@time jpapumpon = hbsolve(ws, wp, sourcespumpon, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs, dc = true, threewavemixing=true,fourwavemixing=true) # enable dc and three wave mixing
@time jpapumpoff = hbsolve(ws, wp, sourcespumpoff, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs, dc = true, threewavemixing=true,fourwavemixing=true) # enable dc and three wave mixing

# Create figure with CairoMakie
fig = Figure(resolution = (800, 800))

# Top plot - Gain (dB)
ax1 = Axis(fig[1, 1],
    xlabel = "Frequency (GHz)",
    ylabel = "Gain (dB)",
    title = "SNAIL Parametric Amplifier Response"
)

# Extract frequency and S-parameters
freqs_ghz = jpapumpon.linearized.w/(2*pi*1e9)

# S11 for pump on
s11_pumpon_db = 10*log10.(abs2.(
    jpapumpon.linearized.S(
        outputmode=(0,),
        outputport=1,
        inputmode=(0,),
        inputport=1,
        freqindex=:
    )
))

# S11 for pump off
s11_pumpoff_db = 10*log10.(abs2.(
    jpapumpoff.linearized.S(
        outputmode=(0,),
        outputport=1,
        inputmode=(0,),
        inputport=1,
        freqindex=:
    )
))

lines!(ax1, freqs_ghz, s11_pumpon_db, label="pump on", linewidth=2)
lines!(ax1, freqs_ghz, s11_pumpoff_db, label="pump off", linewidth=2)
axislegend(ax1)

# Bottom plot - Phase
ax2 = Axis(fig[2, 1],
    xlabel = "Frequency (GHz)",
    ylabel = "Phase (rad)",
    title = "Phase Response"
)

# Phase for pump on
phase_pumpon = angle.(
    jpapumpon.linearized.S(
        outputmode=(0,),
        outputport=1,
        inputmode=(0,),
        inputport=1,
        freqindex=:
    )
)

# Phase for pump off
phase_pumpoff = angle.(
    jpapumpoff.linearized.S(
        outputmode=(0,),
        outputport=1,
        inputmode=(0,),
        inputport=1,
        freqindex=:
    )
)

lines!(ax2, freqs_ghz, phase_pumpon, label="pump on", linewidth=2)
lines!(ax2, freqs_ghz, phase_pumpoff, label="pump off", linewidth=2)
axislegend(ax2)

# Save the figure
save("SNAILPA_results.png", fig)

# Display the figure
fig