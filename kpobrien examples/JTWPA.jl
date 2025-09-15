using JosephsonCircuits
using CairoMakie

@variables Rleft Rright Cg Lj Cj Cc Cr Lr
circuit = Tuple{String,String,String,Num}[]

# port on the input side
push!(circuit,("P$(1)_$(0)","1","0",1))
push!(circuit,("R$(1)_$(0)","1","0",Rleft))
Nj=2048
pmrpitch = 4
#first half cap to ground
push!(circuit,("C$(1)_$(0)","1","0",Cg/2))
#middle caps and jj's
push!(circuit,("Lj$(1)_$(2)","1","2",Lj)) 
push!(circuit,("C$(1)_$(2)","1","2",Cj)) 

j=2
for i = 2:Nj-1
    
    if mod(i,pmrpitch) == pmrpitch÷2

        # make the jj cell with modified capacitance to ground
        push!(circuit,("C$(j)_$(0)","$(j)","$(0)",Cg-Cc))
        push!(circuit,("Lj$(j)_$(j+2)","$(j)","$(j+2)",Lj))

        push!(circuit,("C$(j)_$(j+2)","$(j)","$(j+2)",Cj))
        
        #make the pmr
        push!(circuit,("C$(j)_$(j+1)","$(j)","$(j+1)",Cc))
        push!(circuit,("C$(j+1)_$(0)","$(j+1)","$(0)",Cr))
        push!(circuit,("L$(j+1)_$(0)","$(j+1)","$(0)",Lr))
        
        # increment the index
        global j+=1
    else
        push!(circuit,("C$(j)_$(0)","$(j)","$(0)",Cg))
        push!(circuit,("Lj$(j)_$(j+1)","$(j)","$(j+1)",Lj))
        push!(circuit,("C$(j)_$(j+1)","$(j)","$(j+1)",Cj))
    end
    
    # increment the index
    global j+=1

end

#last jj
push!(circuit,("C$(j)_$(0)","$(j)","$(0)",Cg/2))
push!(circuit,("R$(j)_$(0)","$(j)","$(0)",Rright))
# port on the output side
push!(circuit,("P$(j)_$(0)","$(j)","$(0)",2))

circuitdefs = Dict(
    Lj => IctoLj(3.4e-6),
    Cg => 45.0e-15,
    Cc => 30.0e-15,
    Cr =>  2.8153e-12,
    Lr => 1.70e-10,
    Cj => 55e-15,
    Rleft => 50.0,
    Rright => 50.0,
)

ws=2*pi*(1.0:0.1:14)*1e9
wp=(2*pi*7.12*1e9,)
Ip=1.85e-6

sources = [(mode=(1,),port=1,current=Ip)]
Npumpharmonics = (20,)
Nmodulationharmonics = (10,)

@time rpm = hbsolve(ws, wp, sources, Nmodulationharmonics,
    Npumpharmonics, circuit, circuitdefs)

# Create figure with CairoMakie
fig = Figure(resolution = (1200, 800))

# S-parameters plot
ax1 = Axis(fig[1, 1], 
    xlabel = "Signal Frequency (GHz)",
    ylabel = "dB",
    title = "Scattering Parameters",
    limits = (nothing, (-40, 30))
)

freqs_ghz = ws/(2*pi*1e9)

# S21
s21_db = 10*log10.(abs2.(rpm.linearized.S(
    outputmode=(0,),
    outputport=2,
    inputmode=(0,),
    inputport=1,
    freqindex=:)
))
lines!(ax1, freqs_ghz, s21_db, label="S21")

# S12
s12_db = 10*log10.(abs2.(rpm.linearized.S((0,),1,(0,),2,:)))
lines!(ax1, freqs_ghz, s12_db, label="S12")

# S11
s11_db = 10*log10.(abs2.(rpm.linearized.S((0,),1,(0,),1,:)))
lines!(ax1, freqs_ghz, s11_db, label="S11")

# S22
s22_db = 10*log10.(abs2.(rpm.linearized.S((0,),2,(0,),2,:)))
lines!(ax1, freqs_ghz, s22_db, label="S22")

axislegend(ax1, position = :rb)

# Quantum efficiency plot
ax2 = Axis(fig[1, 2],
    xlabel = "Signal Frequency (GHz)",
    ylabel = "QE/QE_ideal",
    title = "Quantum efficiency",
    limits = (nothing, (0, 1.05))
)

qe_ratio = rpm.linearized.QE((0,),2,(0,),1,:)./rpm.linearized.QEideal((0,),2,(0,),1,:)
lines!(ax2, freqs_ghz, qe_ratio)

# All idlers plot
ax3 = Axis(fig[2, 1],
    xlabel = "Signal Frequency (GHz)",
    ylabel = "dB",
    title = "All idlers",
    limits = (nothing, (-40, 30))
)

all_idlers_db = 10*log10.(abs2.(rpm.linearized.S(:,2,(0,),1,:)'))
for i in 1:size(all_idlers_db, 2)
    lines!(ax3, freqs_ghz, all_idlers_db[:, i])
end

# Commutation relation error plot
ax4 = Axis(fig[2, 2],
    xlabel = "Signal Frequency (GHz)",
    ylabel = "Commutation relation error",
    title = "Commutation relation error"
)

comm_error = 1 .- rpm.linearized.CM((0,),2,:)
lines!(ax4, freqs_ghz, comm_error)

# Save the figure
save("TWPA_results.png", fig)

# Display the figure
fig