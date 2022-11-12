
function tdsolvexyce(ws::AbstractArray{Float64,1},wp::AbstractFloat,Ip::AbstractFloat,circuit,circuitdefs,port=1,jj=true)


    c = parsecomponents(circuit,circuitdefs);

    # # calculate the component dictionaries
    # c = consolidatecomponents(c0.typearray,c0.typearrayw,nodearraysorted,
    #     nodearraysortedw,c0.mutualinductorbranches,c0.mutualinductorbranchesw,
    #     c0.valuearray, c0.valuearrayw,c0.uniquenodearray)
    # cdict = c.cdict

    # n = exportnetlist(c.cdict,c.Nnodes,port,jj)

    n = exportnetlist(circuit,circuitdefs,port=port,jj=jj)
    netlist = n.netlist

    #convert from angular frequency
    fp = wp/(2*pi)


    #define the frequencies and currents
    # fs = 6.0e9
    Is = 1e-13

    stepsperperiod = 80
    tstep = 1/(fp*stepsperperiod)
    tstop = 200e-9
    trise = 10e-9

    inputpump = xyce_input_transient(netlist,Ip,fp,0.0,0*Is,0.0,0.0,tstep,tstop,trise);
    input = [inputpump]

    for i = 1:length(ws)
        fs = ws[i]/(2*pi)

        inputsin = xyce_input_transient(netlist,Ip,fp,0.0,Is,fs,0.0,tstep,tstop,trise);
        push!(input,inputsin)

        inputcos = xyce_input_transient(netlist,Ip,fp,0.0,Is,fs,pi/2,tstep,tstop,trise);
        push!(input,inputcos)
    end

    out = spice_run(input,xyce_cmd())
    return out,input

end


function tdsolvewrspice(ws::AbstractArray{Float64,1},wp::AbstractFloat,Ip::AbstractFloat,circuit,circuitdefs,port=1,jj=true)


    c = parsecomponents(circuit,circuitdefs);

    # # calculate the component dictionaries
    # c = consolidatecomponents(c0.typearray,c0.typearrayw,nodearraysorted,
    #     nodearraysortedw,c0.mutualinductorbranches,c0.mutualinductorbranchesw,
    #     c0.valuearray, c0.valuearrayw,c0.uniquenodearray)
    # cdict = c.cdict

    # n = exportnetlist(c.cdict,c.Nnodes,port,jj)

    n = exportnetlist(circuit,circuitdefs,port=port,jj=jj)
    netlist = n.netlist

    #convert from angular frequency
    fp = wp/(2*pi)


    #define the frequencies and currents
    # fs = 6.0e9
    Is = 1e-13

    stepsperperiod = 80
    tstep = 1/(fp*stepsperperiod)
    tstop = 200e-9
    trise = 10e-9

    inputpump = wrspice_input_transient(netlist,Ip,fp,0.0,0*Is,0.0,0.0,tstep,tstop,trise);
    input = [inputpump]

    for i = 1:length(ws)
        fs = ws[i]/(2*pi)

        inputsin = wrspice_input_transient(netlist,Ip,fp,0.0,Is,fs,0.0,tstep,tstop,trise);
        push!(input,inputsin)

        inputcos = wrspice_input_transient(netlist,Ip,fp,0.0,Is,fs,pi/2,tstep,tstop,trise);
        push!(input,inputcos)
    end

    out = spice_run(input,wrspice_cmd())
    return out,input

end


"""

ideally i should make this function automatically find the ports.
and calculate the multi-frequency scattering matrix. similar to how i do it
for the pump.

"""

function calcgxyce(out,wswrspice,circuit,circuitdefs,port=1,jj=true)

    stepsperperiod=20

    c = parsecomponents(circuit,circuitdefs)
    # n = exportnetlist(c.cdict,c.Nnodes,port,jj)
    n = exportnetlist(circuit,circuitdefs,port=port,jj=jj)

    netlist = n.netlist
    Nnodes = length(c.uniquenodearray)-1


    # g = zeros(Float64,length(wswrspice))
    S11 = zeros(Complex{Float64},length(wswrspice))
    S21 = zeros(Complex{Float64},length(wswrspice))

    for i = 1:length(wswrspice)

        vpump = out[1][1].values["voltage"][1:Nnodes-1,end-stepsperperiod+1:end]
        vsinpump = out[1][2*i].values["voltage"][1:Nnodes-1,end-stepsperperiod+1:end]
        vcospump = out[1][2*i+1].values["voltage"][1:Nnodes-1,end-stepsperperiod+1:end]
        vsin = vsinpump .- vpump
        vcos = vcospump .- vpump    

        phasepump = out[1][1].values["voltage"][Nnodes:end,end-stepsperperiod+1:end]
        phasesinpump = out[1][2*i].values["voltage"][Nnodes:end,end-stepsperperiod+1:end]
        phasecospump = out[1][2*i+1].values["voltage"][Nnodes:end,end-stepsperperiod+1:end]
        phasesin = phasesinpump .- phasepump
        phasecos = phasecospump .- phasepump

        t = out[1][1].values["time"][end-stepsperperiod+1:end]

        vsincos = vcos .+ im*vsin
        # vsincos = vcos .- im*vsin

        phasesincos = phasecos .+ im*phasesin

        for k = 1:(size(phasesincos)[1])
            for j = 1:stepsperperiod
                phasesincos[k,j] = phasesincos[k,j]*exp.(-im*wswrspice[i]*t[j])
            end
        end

        for k = 1:Nnodes-1
            for j = 1:stepsperperiod
                vsincos[k,j] = vsincos[k,j]*exp.(-im*wswrspice[i]*t[j])
                # vsincos[k,j] = vsincos[k,j]*exp.(im*wswrspice[i]*t[j])
            end
        end

        # println(size(vsincos))
        # println(size(phasesincos))

        ftvsincos = FFTW.fft(vsincos,[2])/stepsperperiod;
        ftphasesincos = FFTW.fft(phasesincos,[2])/stepsperperiod;



        #calculate S11 and S21 from the wrspice simulation
        Vth=1e-13*50

        V1=ftvsincos[1,1]
        V2=ftvsincos[end,1]
        # phase1=ftphasesincos[1,1]
        # phase2=ftphasesincos[end,1]

        S11[i] =(2 .* V1 .- Vth) ./ Vth
        S21[i] = 2*V2 ./ Vth .* sqrt(50/50)

    end
    return (S11=S11,S21=S21)
end


function calcgwrspice(out,wswrspice,circuit,circuitdefs,port=1,jj=true)

    stepsperperiod=20

    c = parsecomponents(circuit,circuitdefs)
    # n = exportnetlist(c.cdict,c.Nnodes,port,jj)
    n = exportnetlist(circuit,circuitdefs,port=port,jj=jj)

    netlist = n.netlist
    Nnodes = length(c.uniquenodearray)-1


    # g = zeros(Float64,length(wswrspice))
    S11 = zeros(Complex{Float64},length(wswrspice))
    S21 = zeros(Complex{Float64},length(wswrspice))

    for i = 1:length(wswrspice)

        vpump = out[1][1].values["V"][1:Nnodes-1,end-stepsperperiod+1:end]
        vsinpump = out[1][2*i].values["V"][1:Nnodes-1,end-stepsperperiod+1:end]
        vcospump = out[1][2*i+1].values["V"][1:Nnodes-1,end-stepsperperiod+1:end]
        vsin = vsinpump .- vpump
        vcos = vcospump .- vpump    

        phasepump = out[1][1].values["V"][Nnodes:end,end-stepsperperiod+1:end]
        phasesinpump = out[1][2*i].values["V"][Nnodes:end,end-stepsperperiod+1:end]
        phasecospump = out[1][2*i+1].values["V"][Nnodes:end,end-stepsperperiod+1:end]
        phasesin = phasesinpump .- phasepump
        phasecos = phasecospump .- phasepump

        t = out[1][1].values["S"][end-stepsperperiod+1:end]

        vsincos = vcos .+ im*vsin
        # vsincos = vcos .- im*vsin

        phasesincos = phasecos .+ im*phasesin

        for k = 1:(size(phasesincos)[1])
            for j = 1:stepsperperiod
                phasesincos[k,j] = phasesincos[k,j]*exp.(-im*wswrspice[i]*t[j])
            end
        end

        for k = 1:Nnodes-1
            for j = 1:stepsperperiod
                vsincos[k,j] = vsincos[k,j]*exp.(-im*wswrspice[i]*t[j])
                # vsincos[k,j] = vsincos[k,j]*exp.(im*wswrspice[i]*t[j])
            end
        end

        # println(size(vsincos))
        # println(size(phasesincos))

        ftvsincos = FFTW.fft(vsincos,[2])/stepsperperiod;
        ftphasesincos = FFTW.fft(phasesincos,[2])/stepsperperiod;



        #calculate S11 and S21 from the wrspice simulation
        Vth=1e-13*50

        V1=ftvsincos[1,1]
        V2=ftvsincos[end,1]
        # phase1=ftphasesincos[1,1]
        # phase2=ftphasesincos[end,1]

        S11[i] =(2 .* V1 .- Vth) ./ Vth
        S21[i] = 2*V2 ./ Vth .* sqrt(50/50)

    end
    return (S11=S11,S21=S21)
end



function linsolvewrspice(circuit::Vector,circuitdefs::Dict,w,
    port::Int = true,jj::Bool = true)

    # c = parsecircuit(circuit,circuitdefs);

    # n = exportnetlist(c.cdict,c.Nnodes,port,jj)

    n = exportnetlist(circuit,circuitdefs,port=port,jj=jj)

    # netlist = join(n.netlist,"\n")
    netlist = n.netlist

    input = wrspice_input_ac(netlist,w/(2*pi),n.portnodes,n.portcurrent)
    output = spice_run(input,wrspice_cmd())

    return output[2]
    # return input
end


