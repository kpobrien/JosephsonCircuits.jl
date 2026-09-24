
"""
    wrspice_input_transient(netlist::String, current, frequency, phase,
        sourcenodes, tstep, tstop, trise; maxdata = 2e9, jjaccel = 1,
        dphimax = 0.01, filetype = "binary")

Generate the WRSPICE input for a transient simulation of the circuit in
`netlist`, driven by one sinusoidal current source per entry of `current`,
`frequency`, `phase` and `sourcenodes`, with the time step and stop time
given. The output file name is left out of the `write` command so it can
be given on the command line, and no variables are named so that every
node is saved.

# Arguments
- `netlist`: String containing the circuit netlist, excluding sources.
- `current`: Vector of current source amplitudes in Ampere.
- `frequency`: Vector of current source frequencies in Hz.
- `phase`: Vector of current source phases in radians.
- `sourcenodes`: Vector of tuples of nodes `(src, dst)` at which to place the
    current source(s). A source draws its current from `src` and injects it
    into `dst`, as SPICE orients a current source, so `(0, 1)` drives node
    `1` from ground.
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
                throw(ArgumentError(lazy"Source nodes not strings or integers."))
            end
        end
    end

    if length(current) != length(frequency) || length(current) != length(phase) || length(current) != length(sourcenodes)
        throw(ArgumentError(lazy"Input vector lengths not equal."))
    end

    for s in sourcenodes
        if length(s) != 2
            throw(ArgumentError(lazy"Two nodes are required per source."))
        end
        if !(eltype(s) <: String || eltype(s) <: Int)
            throw(ArgumentError(lazy"Nodes are not an integer or string."))
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
    wrspice_input_ac(netlist, nsteps, fstart, fstop, portnodes, portcurrent;
        maxdata = 2e9)
    wrspice_input_ac(netlist, freqs, portnodes, portcurrent; maxdata = 2e9)

Generate the WRSPICE input for an AC small signal simulation of the circuit
in `netlist`, driven by an AC current source of the complex amplitude
`portcurrent` across the node pair `portnodes`, given as one-based node
indices with ground as index 1: the source draws its current from
`portnodes[1]` and injects it into `portnodes[2]`, each decremented to its
SPICE node label, so `[1, 2]` drives the SPICE node `1` from ground, with
the magnitude and the phase of `portcurrent`, the phase written in degrees
as SPICE reads it. The analysis is `.ac lin nsteps fstart fstop`, over
linearly spaced frequencies from `fstart` to `fstop` in Hz. The second
form takes the frequencies as a single number, or as a vector or range of
which only the first and last entries are used, with `length(freqs) - 2`
passed as `nsteps`, which WRSPICE answers with `length(freqs)` points.
`maxdata` is the WRSPICE limit on the size of the data written, in
kilobytes.

# Examples
```jldoctest
julia> println(JosephsonCircuits.wrspice_input_ac("* SPICE Simulation",100,4e9,5e9,[1,2],1e-6))
* SPICE Simulation
* AC current source into the port
isrc 0 1 ac 1.0e-6 0.0

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
* AC current source into the port
isrc 0 1 ac 1.0e-6 0.0

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

    * AC current source into the port
    isrc $(portnodes[1]-1) $(portnodes[2]-1) ac $(abs(portcurrent)) $(rad2deg(angle(portcurrent)))

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

# the command a loaded provider registers, which XicTools_jll's
# extension fills with its wrspice
const wrspicedefaultcmd = Ref{Any}(nothing)

# where WRSPICE installs itself, which `wrspice_cmd` prefers to a loaded
# provider; the tests ask whether an installation is there too.
# Note: This code has been tested on Linux but not macOS or Windows.
wrspicestandardpath() = Sys.iswindows() ?
    "C:/usr/local/xictools/bin/wrspice.bat" : "/usr/local/xictools/bin/wrspice"

"""
    wrspice_cmd()

The command which runs WRSPICE: the executable at WRSPICE's standard
installation path if one is installed there, else the one a loaded
provider registered, which loading the `XicTools_jll` package does on the
platforms its artifact supports. Throws when neither is available;
[`WRspice`](@ref) and [`spice_run`](@ref) take an executable directly
for one installed elsewhere.
"""
function wrspice_cmd()
    wrspicecmd = wrspicestandardpath()
    (islink(wrspicecmd) || isfile(wrspicecmd)) && return wrspicecmd
    isnothing(wrspicedefaultcmd[]) || return wrspicedefaultcmd[]
    error("WRSPICE executable not found. Please install WRSPICE, load XicTools_jll, or supply a path manually if installed elsewhere.")
end

"""
    spice_run(input, spicecmd)

Run WRSPICE in batch mode on the input `input`, a string, with the command
`spicecmd`, the path of the executable or the command [`wrspice_cmd`](@ref)
returns, and return the rawfile the run writes as [`spice_raw_load`](@ref)
reads it. The input and the rawfile are written in a temporary directory,
which is removed when the run ends, however it ends. The `write` command of
the input's control block must not name a file, so that the rawfile is the
one given on the command line.

A run which fails, writes no rawfile, or aborts, when WRSPICE writes the
plot of its constants in place of the analysis, is an error which quotes
what WRSPICE printed.
"""
function spice_run(input, spicecmd)
    return mktempdir() do path
        inputfilename = joinpath(path, "spice.cir")
        outputfilename = joinpath(path, "spice.raw")
        write(inputfilename, input)
        # run in batch mode, keeping what WRSPICE prints for a failure
        printed, errors = IOBuffer(), IOBuffer()
        process = run(pipeline(
            ignorestatus(`$spicecmd -b -r $outputfilename $inputfilename`);
            stdout = printed, stderr = errors))
        failed(why) = error("WRSPICE $(why); it printed:\n" *
            strip(String(take!(printed)) * "\n" * String(take!(errors))))
        success(process) || failed("exited with code $(process.exitcode)")
        isfile(outputfilename) || failed("wrote no rawfile")
        output = spice_raw_load(outputfilename)
        lowercase(output.header.plotname) == "constants" &&
            failed("aborted the analysis and wrote its constants")
        return output
    end
end

"""
    spice_run(inputs::AbstractVector, spicecmd; ntasks::Int = Sys.CPU_THREADS)

Run each input of `inputs` as [`spice_run`](@ref) does, `ntasks` WRSPICE
processes at a time, and return their rawfiles in the order of the
inputs. The processes run in parallel whatever the number of Julia
threads, since the tasks only wait on them, so `ntasks` is the number of
logical processors by default.
"""
function spice_run(inputs::AbstractVector, spicecmd; ntasks::Int = Sys.CPU_THREADS)
    return asyncmap(input -> spice_run(input, spicecmd), inputs; ntasks = ntasks)
end

"""
    spice_hb_load(filename)

Load the frequency domain output of a Xyce harmonic balance simulation, a
`.HB.FD.prn` file. Returns a named tuple with the value of each output
variable at each frequency in `data`, one row per variable in the order
of the header, the name of each row in `variables`, the frequencies in
`f`, the index column in `index`, and the column names in `header`. A
variable printed as the two columns `Re(name)` and `Im(name)` is the row
`name`, its columns paired by their names; any other column, an
expression `{...}` say, is a row of its own, with the column as its real
part.
"""
function spice_hb_load(filename)

    data = Float64[]
    header = SubString{String}[]

    open(filename, "r") do io
        for line in eachline(io)
            s = strip(line)
            isempty(s) && continue
            if startswith(s, "Index")
                append!(header, split(s, r"\s+"))
            elseif s == "End of Xyce(TM) Simulation"
                break
            else
                append!(data, parse.(Float64, split(s, r"\s+")))
            end
        end
    end

    values = reshape(data, length(header), :)
    index = values[1, :]
    f = values[2, :]

    # the columns of each variable, its real and imaginary parts by name
    variables = String[]
    recolumn = Int[]
    imcolumn = Int[]
    for c in 3:length(header)
        m = match(r"^(Re|Im)\((.*)\)$", header[c])
        name = isnothing(m) ? String(header[c]) : String(m.captures[2])
        k = findfirst(==(name), variables)
        if isnothing(k)
            push!(variables, name)
            push!(recolumn, 0)
            push!(imcolumn, 0)
            k = length(variables)
        end
        if !isnothing(m) && m.captures[1] == "Im"
            imcolumn[k] = c
        else
            recolumn[k] = c
        end
    end

    data1 = zeros(Complex{Float64}, length(variables), size(values, 2))
    for k in eachindex(variables)
        recolumn[k] > 0 && (data1[k, :] .+= view(values, recolumn[k], :))
        imcolumn[k] > 0 && (data1[k, :] .+= im .* view(values, imcolumn[k], :))
    end

    return (data=data1, f=f, index=index, header=header, variables=variables)
end
