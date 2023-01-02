"""
    wrspice_input_paramp(netlist, ws, wp, Ip; stepsperperiod = 80, Is = 1e-13,
        tstop = 200e-9, trise = 10e-9)

Generate the WRSPICE input files for a time domain domain simulation in which
the signal angular frequency is swept over `ws` with pump angular frequency
`wp` and pump current `Ip`. 
"""
function wrspice_input_paramp(netlist, ws, wp, Ip; stepsperperiod = 80, Is = 1e-13,
    tstop = 200e-9, trise = 10e-9)

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
    wrspice_calcS_paramp(out, wswrspice, Nnodes, stepsperperiod = 80,
        Is = 1e-13)

This function assume the first node is the input port and the last node is the
output port.

Nnodes is the number of nodes including the ground node.


# Examples
```jldoctest
Nnodes = 2
ws = 2*pi*5e9
wp = 2*pi*6e9
stepsperperiod = 80
t = LinRange(0,2*pi/wp,stepsperperiod)
Is = 1e-13
Vpump = zeros(1,length(t))
Vsignalsin = zeros(1,length(t))
Vsignalcos = zeros(1,length(t))
Vpump[1,:] .= sin.(2*pi*wp*t)
Vsignalsin[1,:] .= sin.(2*pi*wp*t)+Is/50*sin.(2*pi*ws*t)
Vsignalcos[1,:] .= sin.(2*pi*wp*t)+Is/50*cos.(2*pi*ws*t)

out = [[(values=Dict("S"=>t,"V"=>Vpump),),
        (values=Dict("S"=>t,"V"=>Vsignalsin),),
        (values=Dict("S"=>t,"V"=>Vsignalcos),),
        ]];
out[1][1].values["V"];
out[1][2].values["V"];
out[1][3].values["V"];
JosephsonCircuits.wrspice_calcS_paramp(out, 2*pi, Nnodes;
    stepsperperiod = stepsperperiod, Is = Is)

# output
(S11 = ComplexF64[-0.9999710828404902 + 2.7521387922213123e-5im], S21 = ComplexF64[2.8917159509832197e-5 + 2.7521387922213123e-5im])
```
"""
function wrspice_calcS_paramp(out, wswrspice, Nnodes; stepsperperiod = 80,
    Is = 1e-13)

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
        Vth=Is*50

        V1=ftvsincos[1,1]
        V2=ftvsincos[end,1]
        # phase1=ftphasesincos[1,1]
        # phase2=ftphasesincos[end,1]

        S11[i] =(2 .* V1 .- Vth) ./ Vth
        S21[i] = 2*V2 ./ Vth .* sqrt(50/50)

    end
    return (S11=S11,S21=S21)
end

