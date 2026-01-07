
"""
    wrspice_input_transient(netlist::String, current, frequency, phase, tstep,
        tstop, trise; maxdata = 2e9, jjaccel = 1, dphimax = 0.01,
        filetype = "binary")

Generate the WRSPICE input file for a transient simulation using circuit 
parameters from the given netlist, the source current, frequency phase, and
nodes for the sources, and the time step and stop time. Leave filename empty
so we can add that as a command line argument. Don't specify any variables
so it saves everything.

# Arguments
- `netlist`: String containing the circuit netlist, excluding sources.
- `current`: Vector of current source amplitudes in Ampere.
- `frequency`: Vector of current source frequencies in Hz.
- `phase`: Vector of current source phases in radians.
- `sourcenodes`: Vector of tuples of nodes (src,dst) at which to place the
    current source(s).
- `tstep`: Time step in seconds.
- `tstop`: Time for which to run the simulation in seconds.
- `trise`: The simulation ramps up the current source amplitude with a
    1-sech(t/trise) envelope which reaches 35 percent of the peak in one
    `trise`.

# Keywords
- `maxdata = 2e9`: Maximum size of data to export in kilobytes from 1e3 to
    2e9 with WRspice default 2.56e5. This has to come before the run command.
- `jjaccel = 1`: Causes a faster convergence testing and iteration control
    algorithm to be used, rather than the standard more comprehensive
    algorithm suitable for all devices.
- `dphimax = 0.01`: The maximum allowed phase change per time step. Decreasing
    dphimax from the default of pi/5 to a smaller value is critical for
    matching the accuracy of the harmonic balance method simulations. This
    increases simulation time by pi/5/(dphimax).
- `filetype = "binary" or "ascii"`: Binary files are faster to save and load.

# Examples
```jldoctest
julia> println(JosephsonCircuits.wrspice_input_transient("* SPICE Simulation",[1e-6,1e-3],[5e9,6e9],[3.14,6.28],[(1,0),(1,0)],1e-9,100e-9,10e-9))
* SPICE Simulation
* Current source
* 1-hyperbolic secant rise
isrc1 1 0 1.0u*cos(31.41592653589793g*x+3.14)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))
isrc2 1 0 1000.0u*cos(37.69911184307752g*x+6.28)*(1-2/(exp(x/1.0e-8)+exp(-x/1.0e-8)))
* Set up the transient simulation
* .tran 5p 10n
.tran 1000.0000000000001p 100.0n uic

* The control block
.control
set maxdata=2.0e9
set jjaccel=1
set dphimax=0.01
run
set filetype=binary
write
.endc
```
"""
function wrspice_input_transient(netlist::String, current, frequency, phase,
    sourcenodes, tstep, tstop, trise; maxdata = 2e9, jjaccel = 1,
    dphimax = 0.01, filetype = "binary")

    if length(current) == 1
        if length(sourcenodes) == 2 
            if (eltype(sourcenodes) <: String || eltype(sourcenodes) <: Int)
                sourcenodes = [sourcenodes]
            else
                throw(ArgumentError("Source nodes not strings or integers."))
            end
        end
    end

    if length(current) != length(frequency) || length(current) != length(phase) || length(current) != length(sourcenodes)
        throw(ArgumentError("Input vector lengths not equal."))
    end

    for s in sourcenodes
        if length(s) != 2
            throw(ArgumentError("Two nodes are required per source."))
        end
        if !(eltype(s) <: String || eltype(s) <: Int)
            throw(ArgumentError("Nodes are not an integer or string."))
        end
    end

    control = ""
    control *="""

    * Current source
    * 1-hyperbolic secant rise
    """

    for i in 1:length(current)
        control*="""isrc$(i) $(sourcenodes[i][1]) $(sourcenodes[i][2]) $(current[i]*1e6)u*cos($(2*pi*frequency[i]*1e-9)g*x+$(phase[i]))*(1-2/(exp(x/$trise)+exp(-x/$trise)))\n"""
    end

    control *="""
    * Set up the transient simulation
    * .tran 5p 10n
    .tran $(tstep*1e12)p $(tstop*1e9)n uic

    * The control block
    .control
    set maxdata=$(maxdata)
    set jjaccel=$(jjaccel)
    set dphimax=$(dphimax)
    run
    set filetype=$(filetype)
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
set maxdata=2.0e9

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
set maxdata=2.0e9

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
    portnodes,portcurrent; maxdata = 2e9)
    if length(freqs) == 1
        return wrspice_input_ac(netlist,1,freqs[1],freqs[1],portnodes,portcurrent; maxdata = maxdata)
    else
        return wrspice_input_ac(netlist,length(freqs)-2,freqs[1],freqs[end],portnodes,portcurrent; maxdata = maxdata)
    end
end

function wrspice_input_ac(netlist::String,freqs::AbstractRange{Float64},
    portnodes,portcurrent; maxdata = 2e9)

    return wrspice_input_ac(netlist,length(freqs)-2,freqs[1],freqs[end],portnodes,portcurrent; maxdata = maxdata)
end

function wrspice_input_ac(netlist::String,freqs::Float64,
    portnodes,portcurrent; maxdata = 2e9)

    return wrspice_input_ac(netlist,1,freqs,freqs,portnodes,portcurrent; maxdata = maxdata)
end

function wrspice_input_ac(netlist,nsteps,fstart,fstop,portnodes,portcurrent; maxdata = 2e9)

    control="""

    * AC current source with magnitude 1 and phase 0
    isrc $(portnodes[2]-1) $(portnodes[1]-1) ac $(abs(portcurrent)) $(angle(portcurrent))

    * Set up the AC small signal simulation
    .ac lin $(nsteps) $(fstart*1e-9)g $(fstop*1e-9)g

    * The control block
    .control

    * Maximum size of data to export in kilobytes from 1e3 to 2e9 with
    * default 2.56e5. This has to come before the run command
    set maxdata=$(maxdata)

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

    if Sys.iswindows()
        wrspicecmd = "C:/usr/local/xictools/bin/wrspice.bat"
    # Note: This code has been tested on Linux but not macOS or Windows. 
    else
        wrspicecmd = "/usr/local/xictools/bin/wrspice"
    end

    # Check if the wrspice executable exists
    if !islink(wrspicecmd) && !isfile(wrspicecmd)
        # if the path isn't valid, try checking if we have imported
        # XicTools_jll and use that.
        if :XicTools_jll in names(Main,imported=true) && isdefined(Main.XicTools_jll,:wrspice)
           wrspicecmd =  Main.XicTools_jll.wrspice()
        else
            error("WRSPICE executable not found. Please install WRSPICE, load XicTools_jll, or supply a path manually if installed elsewhere.")
        end
    end

    return wrspicecmd

end

"""
    spice_run(input, spicecmd::String)

Argument is a string or command containing the input commands for wrspice. 
This function saves the string to disk, runs spice, parses the results with
wrsplice_load(), then returns those parsed results.

The input should not should have a file name listed after the write command in
the .control block so that we can specify the raw output file with a command
line argument.

"""
function spice_run(input,spicecmd)

    #find the temporary directory
    path = tempdir()

    #generate unique filenames
    inputfilename = joinpath(path,"spice-"* string(UUIDs.uuid1()) * ".cir")
    outputfilename = joinpath(path,"spice-"* string(UUIDs.uuid1()) * ".raw")

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
    spice_run(input::AbstractVector,spicecmd; ntasks = Threads.nthreads())

If the input to wrspice_run() is an array of strings, then call multiple
processes in parallel. The number of parallel processes is decided from
Threads.nthreads(). It can be changed manually.

"""
function spice_run(inputs::AbstractVector,spicecmd;ntasks::Int = Threads.nthreads())
    # Set the number of simultaneous parallel simulations equal to
    # the number of threads. ntasks can be manually changed here
    # or by changing the number of threads. Note this is only faster because
    # we are calling an external program which launches a new processes.
    tsk(input) = spice_run(input,spicecmd)
    return asyncmap(tsk,inputs;ntasks=ntasks);
end

"""
    spice_hb_load(filename)

Load a Xyce harmonic balance simulation.
"""
function spice_hb_load(filename)

    data = Array{Float64}(undef,0)
    header = []
    i=0

    #open a file handle
    open(filename, "r") do io
        #loop over the contents of the header and data
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

