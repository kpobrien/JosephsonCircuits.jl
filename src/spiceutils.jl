"""
    tdsolvewrspice(ws, wp, Ip, circuit, circuitdefs; port = 1, jj = true,
        stepsperperiod = 80, Is = 1e-13, tstop = 200e-9, trise = 10e-9)

Generate the WRSPICE input files for a time domain domain simulation in which
the signal angular frequency is swept over `ws` with pump angular frequency
`wp` and pump current `Ip`. 
"""
function tdsolvewrspice(ws, wp, Ip, circuit, circuitdefs; port = 1, jj = true,
    stepsperperiod = 80, Is = 1e-13, tstop = 200e-9, trise = 10e-9)

    n = exportnetlist(circuit,circuitdefs,port=port,jj=jj)
    netlist = n.netlist

    #convert from angular frequency
    fp = wp/(2*pi)

    #define the signal frequencies and currents
    # Is = 1e-13

    # define the time step
    tstep = 1/(fp*stepsperperiod)

    # pump only simulation
    inputpump = wrspice_input_transient(netlist,Ip,fp,0.0,0*Is,0.0,0.0,tstep,tstop,trise);
    input = [inputpump]

    for i = 1:length(ws)
        fs = ws[i]/(2*pi)

        inputsin = wrspice_input_transient(netlist,Ip,fp,0.0,Is,fs,0.0,tstep,tstop,trise);
        push!(input,inputsin)

        inputcos = wrspice_input_transient(netlist,Ip,fp,0.0,Is,fs,pi/2,tstep,tstop,trise);
        push!(input,inputcos)
    end

    return input

end

"""
    calcgwrspice(out, wswrspice, circuit, circuitdefs; port=1, jj = true,
        stepsperperiod = 80)


"""
function calcgwrspice(out, wswrspice, circuit, circuitdefs; port=1, jj = true,
    stepsperperiod = 80)

    # generate the netlist
    n = exportnetlist(circuit,circuitdefs,port=port,jj=jj)
    netlist = n.netlist
    Nnodes = n.Nnodes # includes the ground node

    # empty arrays for the scattering parameters
    S11 = zeros(Complex{Float64},length(wswrspice))
    S21 = zeros(Complex{Float64},length(wswrspice))

    for i = 1:length(wswrspice)

        vpump = out[1][1].values["V"][1:Nnodes-1,end-stepsperperiod+1:end]
        vsinpump = out[1][2*i].values["V"][1:Nnodes-1,end-stepsperperiod+1:end]
        vcospump = out[1][2*i+1].values["V"][1:Nnodes-1,end-stepsperperiod+1:end]
        vsin = vsinpump .- vpump
        vcos = vcospump .- vpump

        # the JJ phase is stored as a voltage node after the regular voltage
        # nodes
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

