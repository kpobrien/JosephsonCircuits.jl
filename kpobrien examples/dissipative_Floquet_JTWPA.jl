results = []
tandeltas = [1.0e-6,1.0e-3, 2.0e-3, 3.0e-3]
for tandelta in tandeltas
    circuitdefs = Dict(
        Rleft => 50,
        Rright => 50,
        Lj => IctoLj(1.75e-6),
        Cg => 76.6e-15/(1+im*tandelta),
        Cc => 40.0e-15/(1+im*tandelta),
        Cr => 1.533e-12/(1+im*tandelta),
        Lr => 2.47e-10,
        Cj => 40e-15,
    )  
    wp=(2*pi*7.9*1e9,)
    ws=2*pi*(1.0:0.1:14)*1e9
    Ip=1.1e-6*(1+125*tandelta)
    sources = [(mode=(1,),port=1,current=Ip)]
    Npumpharmonics = (20,)
    Nmodulationharmonics = (10,)
    @time floquet = hbsolve(ws, wp, sources, Nmodulationharmonics,
        Npumpharmonics, circuit, circuitdefs)
    push!(results,floquet)
end

p1 = plot(title="Gain (S21)")
for i = 1:length(results)
        plot!(ws/(2*pi*1e9),
            10*log10.(abs2.(results[i].linearized.S((0,),2,(0,),1,:))),
            ylim=(-60,30),label="tanδ=$(tandeltas[i])",
            legend=:bottomleft,
            xlabel="Signal Frequency (GHz)",ylabel="dB")
end

p2 = plot(title="Quantum Efficiency")
for i = 1:length(results)
        plot!(ws/(2*pi*1e9),
            results[i].linearized.QE((0,),2,(0,),1,:)./results[i].linearized.QEideal((0,),2,(0,),1,:),
            ylim=(0.6,1.05),legend=false,
            title="Quantum efficiency",
            ylabel="QE/QE_ideal",xlabel="Signal Frequency (GHz)")
end

p3 = plot(title="Reverse Gain (S12)")
for i = 1:length(results)
        plot!(ws/(2*pi*1e9),
            10*log10.(abs2.(results[i].linearized.S((0,),1,(0,),2,:))),
            ylim=(-10,1),legend=false,
            xlabel="Signal Frequency (GHz)",ylabel="dB")
end

p4 = plot(title="Commutation \n relation error")
for i = 1:length(results)
        plot!(ws/(2*pi*1e9),
            1 .- results[i].linearized.CM((0,),2,:),
            legend=false,
            ylabel="Commutation\n relation error",xlabel="Signal Frequency (GHz)")
end

plot(p1, p2, p3,p4,layout = (2, 2))