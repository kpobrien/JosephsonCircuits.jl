
"""
    wrspice_input_transient(netlist,idrive,fdrive,tstep,tstop)

Generate the WRSPICE input file for a transient simulation using circuit 
parameters from the given netlist, the source current and frequency, and the 
time step and stop time. Example usage:

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

# Examples
```jldoctest
julia> println(JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",1e-6,5e9,3.14,1e-3,6e9,6.28,1e-9,100e-9,10e-9))
* SPICE Simulation
* Current source
* 1-hyperbolic secant rise
isrc 0 1 1.0u*sin(31.41592653589793g*x+3.14)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))
isrc2 0 1 1000.0u*sin(37.69911184307752g*x+6.28)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))

* isrc 1 0 sin(0 1.0u 5.0g 0.0 179.90874767107852)
* isrc2 1 0 sin(0 1000.0u 6.0g 0.0 359.81749534215703)

* Set up the transient simulation
* .tran 5p 10n
.tran 1000.0000000000001p 100.0n uic

* The control block
.control
set maxdata=1.0e7
set jjaccel=1
set dphimax=0.01
run
set filetype=binary
write
.endc

```
"""
function wrspice_input_transient(netlist,idrive,fdrive,thetadrive,idrive2,
    fdrive2,thetadrive2,tstep,tstop,trise;maxdata=10e6,jjaccel=1,dphimax=0.01)

    control="""

    * Current source
    * 1-hyperbolic secant rise
    isrc 0 1 $(idrive*1e6)u*sin($(2*pi*fdrive*1e-9)g*x+$thetadrive)*(1-2/(exp(x/$trise)+exp(-x/$trise)))
    isrc2 0 1 $(idrive2*1e6)u*sin($(2*pi*fdrive2*1e-9)g*x+$thetadrive2)*(1-2/(exp(x/$trise)+exp(-x/$trise)))

    * isrc 1 0 sin(0 $(idrive*1e6)u $(fdrive*1e-9)g 0.0 $(thetadrive*180/pi))
    * isrc2 1 0 sin(0 $(idrive2*1e6)u $(fdrive2*1e-9)g 0.0 $(thetadrive2*180/pi))

    * Set up the transient simulation
    * .tran 5p 10n
    .tran $(tstep*1e12)p $(tstop*1e9)n uic

    * The control block
    .control
    set maxdata=$(maxdata)
    set jjaccel=$(jjaccel)
    set dphimax=$(dphimax)
    run
    set filetype=binary
    write
    .endc

    """

    input = netlist*control

    return input


end


"""
    wrspice_input_ac(netlist,nsteps,fstart,fstop)

Generate the WRSPICE input file for an AC small signal simulation using circuit
parameters from the given netlist, and the specified frequency range. Example
usage:

# Examples
```jldoctest
julia> println(JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",100,4e9,5e9,[1,2],1e-6))
* SPICE Simulation
* AC current source with magnitude 1 and phase 0
isrc 1 0 ac 1.0e-6 0.0

* Set up the AC small signal simulation
.ac lin 100 4.0g 5.0g

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

```
```jldoctest
julia> println(JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",(4:0.01:5)*1e9,[1,2],1e-6))
* SPICE Simulation
* AC current source with magnitude 1 and phase 0
isrc 1 0 ac 1.0e-6 0.0

* Set up the AC small signal simulation
.ac lin 99 4.0g 5.0g

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

```
```jldoctest
julia> println(JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",collect((4:0.01:5)*1e9),[1,2],1e-6))
* SPICE Simulation
* AC current source with magnitude 1 and phase 0
isrc 1 0 ac 1.0e-6 0.0

* Set up the AC small signal simulation
.ac lin 99 4.0g 5.0g

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

```
"""
function wrspice_input_ac(netlist::String,freqs::AbstractArray{Float64,1},
    portnodes,portcurrent)
    if length(freqs) == 1
        return wrspice_input_ac(netlist,1,freqs[1],freqs[1],portnodes,portcurrent)
    else
        return wrspice_input_ac(netlist,length(freqs)-2,freqs[1],freqs[end],portnodes,portcurrent)
    end
end

function wrspice_input_ac(netlist::String,freqs::AbstractRange{Float64},
    portnodes,portcurrent)

    return wrspice_input_ac(netlist,length(freqs)-2,freqs[1],freqs[end],portnodes,portcurrent)
end

function wrspice_input_ac(netlist::String,freqs::Float64,
    portnodes,portcurrent)

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


"""
    wrspice_cmd()

This returns the path of the WRSPICE executable.
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
        error("WRSPICE executable not found. Please install WRSPICE or supply a path manually if already installed.")
    end

    return wrspicecmd

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
function spice_run(input::String,spicecmd)

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

If the input to wrspice_run() is an array of strings, then call multiple
processes in parallel. The number of parallel processes is decided from
Threads.nthreads(). It can be changed manually. 

"""
function spice_run(input::AbstractArray{String,1},spicecmd;
    ntasks = Threads.nthreads())
    # Set the number of simultaneous parallel simulations equal to 
    # the number of threads. ntasks can be manually changed here 
    # or by changing the number of threads. Note this is only faster because
    # we are calling an external program which launches a new processes. 

    return asyncmap((input) -> spice_run(input,spicecmd),input;ntasks=ntasks);
end

"""
    spice_hb_load(filename)

Load a Xyce harmonic balance simulation.
"""
function spice_hb_load(filename)

    #open a file handle
    io = open(filename)

    data = Array{Float64}(undef,0)
    header = []

    #loop over the contents of the header
    i=0
    while !eof(io)

        i+=1

        line = readline(io)

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

