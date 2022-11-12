
"""

    wrspice_input_transient(netlist,idrive,fdrive,tstep,tstop)

Generate the WRSPICE input file for a transient simulation using circuit 
parameters from the given netlist, the source current and frequency, and the 
time step and stop time. Example usage:

c = twpa_params().c
netlist = twpa_netlist(c,true);
input = wrspice_input_transient(netlist,1e-6,5e9,10e-12,20e-9);
out = wrspice_run(input);

#plot the voltage vs time at the last unit cell node
plot(out[1]*1e9,out[2][c.Nnodes,:]/25)

#plot the voltage vs position at final time
plot(out[2][1:c.Nnodes,end]/25)

#plot the phase across the junctions vs position at the final time 
plot(out[2][c.Nnodes+1:2*c.Nnodes-1,end])

#plot the phase across the junction vs time at the final JJ
plot(out[2][2*c.Nnodes-1,:])

#plot the phase vs position and time
plot(out[2][c.Nnodes+1:2*c.Nnodes-1,end-50:end]',seriestype=:heatmap,
xlabel="position",ylabel="time")
"""
function wrspice_input_transient(netlist,idrive,fdrive,thetadrive,idrive2,fdrive2,thetadrive2,tstep,tstop,trise)

    # Does changing tstart speed up the simulation? No and sometime it slows it 
    # down because the maximum timestep tmax is chosen to be the smaller of 
    # tstep and (tstop-tstart)/50. Not allowing tstart as an input parameter.

    #A.3.6 Circuit 6 on page 389 shows how to use arbitrary functions in sources.
    # trise = 10e-9

    # * isrc 1 0 sin(0 1u 5g)
    # * isrc 1 0 sin(0 $(idrive*1e6)u $(fdrive*1e-9)g)
    # * isrc2 1 0 sin(0 $(idrive2*1e6)u $(fdrive2*1e-9)g)
    # * isrc 1 0 $(idrive)*sin($(2*pi*fdrive)*x)+$(idrive2)*sin($(2*pi*fdrive2)*x)
    # * isrc 1 0 $(idrive*1e6)u*sin($(2*pi*fdrive*1e-9)g*x)+$(idrive2*1e6)u*sin($(2*pi*fdrive2*1e-9)g*x)

    # *isrc 1 0 0.5*(1+erf((x-$trise)*$(1/trise)))*($(idrive*1e6)u*sin($(2*pi*fdrive*1e-9)g*x)+$(idrive2*1e6)u*sin($(2*pi*fdrive2*1e-9)g*x))

    # * exponential
    # * isrc 1 0 exp(0 1 0n $(trise*1e9)n $(tstop*1e9)n $(trise*1e9)n)*sin(0 $(idrive*1e6)u $(fdrive*1e-9)g)
    # * isrc2 1 0 exp(0 1 0n $(trise*1e9)n $(tstop*1e9)n $(trise*1e9)n)*sin(0 $(idrive2*1e6)u $(fdrive2*1e-9)g)


    control="""

    * Current source
    * 1-hyperbolic secant rise
    * isrc 0 1 $(idrive*1e6)u*sin($(2*pi*fdrive*1e-9)g*x+$thetadrive)*(1-2/(exp(x/$trise)+exp(-x/$trise)))
    * isrc2 0 1 $(idrive2*1e6)u*sin($(2*pi*fdrive2*1e-9)g*x+$thetadrive2)*(1-2/(exp(x/$trise)+exp(-x/$trise)))

    * idc 0 1 $(1e-3*1e6)u*(1-2/(exp(x/$trise)+exp(-x/$trise)))

    * isrc 1 0 sin(0 $(idrive*1e6)u $(fdrive*1e-9)g 0.0 $(thetadrive*180/pi))
    * isrc2 1 0 sin(0 $(idrive2*1e6)u $(fdrive2*1e-9)g 0.0 $(thetadrive2*180/pi))


    isrc 1 0 sin(0 $(idrive*1e6)u $(fdrive*1e-9)g)
    isrc2 1 0 sin(0 $(idrive2*1e6)u $(fdrive2*1e-9)g)

    *isrc 0 1 $(idrive*1e6)u*sin($(2*pi*fdrive*1e-9)g*x+$thetadrive)
    *isrc2 0 1 $(idrive2*1e6)u*sin($(2*pi*fdrive2*1e-9)g*x+$thetadrive2)

    * isrc 1 0 sin(0 $(-idrive*1e6)u $(fdrive*1e-9)g 0.0 0.0 $(thetadrive*180/pi))
    * isrc2 1 0 sin(0 $(-idrive2*1e6)u $(fdrive2*1e-9)g 0.0 0.0 $(thetadrive2*180/pi))


    * Set up the transient simulation
    * .tran 5p 10n
    .tran $(tstep*1e12)p $(tstop*1e9)n uic

    * The control block
    .control

    * Maximum size of data to export in kilobytes from 1e3 to 2e9 with 
    * default 2.56e5. This has to come before the run command
    set maxdata = 10e6

    * jjaccel : Causes a faster convergence testing and iteration control algorithm 
    * to be used, rather than the standard more comprehensive algorithm suitable 
    * for all devices
    set jjaccel=1

    * dphimax : The maximum allowed phase change per time step. 
    * decreasing dphimax from the default of pi/5 to a smaller value is critical
    * for matching the accuracy of the harmonic balance method simulations. this
    * increases simulation time by pi/5/(dphimax)
    set dphimax=0.01

    * Run the simulation
    run

    * Binary files are faster to save and load. 
    set filetype=binary

    * Leave filename empty so we can add that as a command line argument.
    * Don't specify any variables so it saves everything.    
    write

    .endc

    """

    input = netlist*control

    return input


end



"""

    xyce_input_transient(netlist,idrive,fdrive,tstep,tstop)

Generate the Xyce input file for a transient simulation using circuit 
parameters from the given netlist, the source current and frequency, and the 
time step and stop time. Example usage:

"""
function xyce_input_transient(netlist,idrive,fdrive,thetadrive,idrive2,fdrive2,thetadrive2,tstep,tstop,trise)

    #A.3.6 Circuit 6 on page 389 shows how to use arbitrary functions in sources.
    # trise = 10e-9

    # * isrc 1 0 sin(0 1u 5g)
    # * isrc 1 0 sin(0 $(idrive*1e6)u $(fdrive*1e-9)g)
    # * isrc2 1 0 sin(0 $(idrive2*1e6)u $(fdrive2*1e-9)g)
    # * isrc 1 0 $(idrive)*sin($(2*pi*fdrive)*x)+$(idrive2)*sin($(2*pi*fdrive2)*x)
    # * isrc 1 0 $(idrive*1e6)u*sin($(2*pi*fdrive*1e-9)g*x)+$(idrive2*1e6)u*sin($(2*pi*fdrive2*1e-9)g*x)

    # *isrc 1 0 0.5*(1+erf((x-$trise)*$(1/trise)))*($(idrive*1e6)u*sin($(2*pi*fdrive*1e-9)g*x)+$(idrive2*1e6)u*sin($(2*pi*fdrive2*1e-9)g*x))

    # * exponential
    # * isrc 1 0 exp(0 1 0n $(trise*1e9)n $(tstop*1e9)n $(trise*1e9)n)*sin(0 $(idrive*1e6)u $(fdrive*1e-9)g)
    # * isrc2 1 0 exp(0 1 0n $(trise*1e9)n $(tstop*1e9)n $(trise*1e9)n)*sin(0 $(idrive2*1e6)u $(fdrive2*1e-9)g)


    control="""

    * Current source
    isrc 1 0 sin(0.0V $(-idrive*1e6)uV $(fdrive*1e-9)g 0.0 0.0 $(thetadrive*180/pi))
    isrc2 1 0 sin(0.0V $(-idrive2*1e6)uV $(fdrive2*1e-9)g 0.0 0.0 $(thetadrive2*180/pi))


    * Set up the transient simulation
    .tran $(tstep*1e12)p  $(tstop*1e9)n 0n $(tstep*1e12)p

    *.options timeint method=trap minord=2 reltol=1e-6 abstol=1e-6

    .end

    """

    input = netlist*control

    return input


end


"""
    wrspice_input_ac(netlist,nsteps,fstart,fstop)

Generate the WRSPICE input file for an AC small signal simulation using circuit
parameters from the given netlist, and the specified frequency range. Example
usage:


c = twpa_params().c
netlist = twpa_netlist(c,true);
input = wrspice_input_ac(netlist,1000,4e9,10e9);
out = wrspice_run(input);

#plot the voltage vs frequency at the last unit cell
plot(out[1]/1e9,10*log10.(abs2.(out[2][c.Nnodes,:]/25)),ylim=(-1,1))

#plot the voltage vs position at the last frequency
plot(real.(out[2][1:c.Nnodes,end]/25))

#plot the phase vs frequency at the last unit cell
plot(out[1]/1e9,real.(out[2][c.Nnodes,:]/25),ylim=(-1,1))

#CANNOT plot the phase vs position at the last frequency because it is not
#calculated in this mode
plot(real.(out[2][c.Nnodes+1:2*c.Nnodes-1,end]/25))

#plot the voltage vs pmr position at the last frequency
plot(real.(out[2][2*c.Nnodes:end,end]/25))

#plot the voltage vs frequency at the first pmr
plot(real.(out[2][2*c.Nnodes,:]/25))


"""

function wrspice_input_ac(netlist::String,freqs::AbstractArray{Float64,1},
    portnodes::Tuple{Int, Int},portcurrent::Complex{Float64})
    if length(freqs) == 1
        return wrspice_input_ac(netlist,1,freqs[1],freqs[1],portnodes,portcurrent)
    else
        return wrspice_input_ac(netlist,length(freqs)-2,freqs[1],freqs[end],portnodes,portcurrent)
    end
end

function wrspice_input_ac(netlist::String,freqs::AbstractRange{Float64},
    portnodes::Tuple{Int, Int},portcurrent::Complex{Float64})

    return wrspice_input_ac(netlist,length(freqs)-2,freqs[1],freqs[end],portnodes,portcurrent)
end

function wrspice_input_ac(netlist::String,freqs::Float64,
    portnodes::Tuple{Int, Int},portcurrent::Complex{Float64})

    return wrspice_input_ac(netlist,1,freqs,freqs,portnodes,portcurrent)
end

function wrspice_input_ac(netlist,nsteps,fstart,fstop,portnodes,portcurrent)

    control="""

    * AC current source with magnitude 1 and phase 0
    isrc $(portnodes[2]-1) $(portnodes[1]-1) ac $(abs(portcurrent)) $(angle(portcurrent))

    * Set up the AC small signal simulation
    .ac lin $(nsteps) $(fstart*1e-9)g $(fstop*1e-9)g

    * The control block
    .control

    * Maximum size of data to export in kilobytes from 1e3 to 2e9 with 
    * default 2.56e5. This has to come before the run command
    set maxdata = 10e6

    * Run the simulation
    run

    * Binary files are faster to save and load. 
    set filetype=binary

    * Leave filename empty so we can add that as a command line argument.
    * Don't specify any variables so it saves everything.    
    write

    .endc

    """

    input = netlist*control

    return input
end



function xyce_input_ac(netlist::String,freqs::AbstractArray{Float64,1},
    portnodes::Tuple{Int, Int},portcurrent::Complex{Float64})
    if length(freqs) == 1
        return xyce_input_ac(netlist,1,freqs[1],freqs[1],portnodes,portcurrent)
    else
        return xyce_input_ac(netlist,length(freqs),freqs[1],freqs[end],portnodes,portcurrent)
    end
end

function xyce_input_ac(netlist::String,freqs::AbstractRange{Float64},
    portnodes::Tuple{Int, Int},portcurrent::Complex{Float64})

    return xyce_input_ac(netlist,length(freqs),freqs[1],freqs[end],portnodes,portcurrent)
end

function xyce_input_ac(netlist::String,freqs::Float64,
    portnodes::Tuple{Int, Int},portcurrent::Complex{Float64})

    return xyce_input_ac(netlist,1,freqs,freqs,portnodes,portcurrent)
end

function xyce_input_ac(netlist,nsteps,fstart,fstop,portnodes,portcurrent)

    control="""

    * AC current source with magnitude 1 and phase 0
    isrc $(portnodes[2]-1) $(portnodes[1]-1) ac $(abs(portcurrent)) $(angle(portcurrent))

    * Set up the AC small signal simulation
    .ac lin $(nsteps) $(fstart*1e-9)g $(fstop*1e-9)g 

    * .options linsol-ac type=klu
    * .print ac

    .end

    """

    input = netlist*control

    return input
end



"""
    wrspice_run(input::String)

Argument is a string containing the input commands for wrspice.  This function 
saves the string to disk, runs wrspice, parses the results with wrsplice_load(), 
then returns those parsed results.

The input should not should have a file name listed after the write command 
in the .control block so that we can specify the raw output file with a command 
line argument.
"""

function wrspice_cmd()
    # Note: This code has been tested on Linux but not macOS or Windows. 
    if Sys.isunix()
        wrspicecmd = "/usr/local/xictools/bin/wrspice"
    elseif Sys.iswindows()
        wrspicecmd = "C:/usr/local/xictools/bin/wrspice.bat"        
    else
        error("Operating system not supported")
    end

    # Check if the wrspice executable exists
    if !islink(wrspicecmd) && !isfile(wrspicecmd)
        error("WRSPICE executable not found at $(wrspicecmd). Please install WRSPICE or change
        the location of the executable (or edit the wrspicecmd line in the code.")
    end

    return wrspicecmd

end

function xyce_cmd()
    # Note: This code has been tested on Linux but not macOS or Windows. 
    if Sys.isunix()
        spicecmd = "/usr/bin/Xyce"
    elseif Sys.iswindows()
        spicecmd = "C:/usr/local/xictools/bin/wrspice.bat"  
        error("Operating system not supported")  
    else
        error("Operating system not supported")
    end

    # Check if the wrspice executable exists
    if !islink(spicecmd) && !isfile(spicecmd)
        error("Xyce executable not found at $(spicecmd). Please install Xyce or change
        the location of the executable (or edit the spicecmd line in the code.")
    end

    return spicecmd

end


"""
    spice_run(input::String,spicecmd::String)

Argument is a string containing the input commands for wrspice.  This function 
saves the string to disk, runs spice, parses the results with wrsplice_load(), 
then returns those parsed results.

The input should not should have a file name listed after the write command 
in the .control block so that we can specify the raw output file with a command 
line argument.
"""
function spice_run(input::String,spicecmd::String)

    # #save a circuit file manually
    # write("julia.cir",input)
    #
    # #load a file as input
    # input = read("julia.cir", String)


    #find the temporary directory
    path = tempdir()

    #generate unique filenames
    inputfilename = joinpath(path,"spice-"* string(UUIDs.uuid1()) * ".cir")
    outputfilename = joinpath(path,"spice-"* string(UUIDs.uuid1()) * ".raw")

    # ##manually load these files for testing
    # inputfilename = "julia.cir"
    # outputfilename = "julia.raw"


    #save the input file
    open(inputfilename, "w") do f
        write(f, input)
    end

    # run the simulation
    # simple version that lets wrspice output go to stdout.
    # run(`$wrspicecmd -b -r $outputfilename $inputfilename`)
    # redirect stdout to the variable reader.
    reader = read(`$spicecmd -b -r $outputfilename $inputfilename`,String)

    # if we need to use this in the future.
    # reader = @async read(`$wrspicecmd -b -r $outputfilename $inputfilename`,String)
    # fetch(reader)

    #parse the output
    output=spice_raw_load(outputfilename)

    #clean up the input and output files
    rm(inputfilename)
    rm(outputfilename)

    return output
end

"""
    spice_run(input::AbstractArray{String,1})

If the input to wrspice_run() is an array of strings, then call 
multiple wrspice processes using in parallel. The number of
parallel processes is decided from Threads.nthreads(). It can
be changed manually. 

"""
function spice_run(input::AbstractArray{String,1},spicecmd::String)
    # Set the number of simultaneous parallel simulations equal to 
    # the number of threads. ntasks can be manually changed here 
    # or by changing the number of threads. Note this is only faster because
    # we are calling an external program which launches a new processes. 
    ntasks = Threads.nthreads()

    return asyncmap((input) -> spice_run(input,spicecmd),input;ntasks=ntasks);
end


function spice_hb_load(filename)

    #open a file handle
    io = open(filename)

    data = Array{Float64}(undef,0)
    header = []

    #loop over the contents of the header
    i=0
    # while true
    # line = readline(io)


    while !eof(io)

        i+=1
        
        line = readline(io)

        # linename = split(line,":",limit=2)[1]

        if line[1:5] == "Index"
            header = split(strip(line),r"\s+")
        elseif line == "End of Xyce(TM) Simulation"
            break
        else
            vals = parse.(Float64,split(strip(line),r"\s+"))
            for val in vals
                push!(data,val)
            end
        end

        # println(line)

    end


    index = data[1:length(header):end]
    f = data[2:length(header):end];

    data1 = zeros(Complex{Float64},(length(header)-2)÷2,length(data) ÷ length(header))

    for i = 1:(length(header)-2)÷2
        data1[i,:] .= data[2+i:length(header):end]
        data1[i,:] .+= im.*data[3+i:length(header):end]
    end    

    return (data=data1,f=f,index=index,header=header)

end

