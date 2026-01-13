
"""
    SpiceRawHeader(title::String, date::String, plotname::String,
        flags::String, nvariables::Int, npoints::Int, command::String,
        option::String)

A simple structure to hold the SPICE raw file header.
"""
struct SpiceRawHeader
    title::String
    date::String
    plotname::String
    flags::String
    nvariables::Int
    npoints::Int
    command::String
    option::String
end

"""
    SpiceRaw(header::SpiceRawHeader, variables::Dict{String, Vector{String}},
        values::Dict{String,T})

A simple structure to hold the SPICE raw file contents including the header,
variables, and values.
"""
struct SpiceRaw{T}
    header::SpiceRawHeader
    variables::Dict{String, Vector{String}}
    values::Dict{String,T}
end

"""
    spice_raw_load(filename)

Parse the binary raw output file from WRSPICE or Xyce. Tested for transient
analysis and frequency domain analysis. The file format is documented in the
[WRSPICE manual](http://www.srware.com/xictools/docs/wrsmanual-4.3.13.pdf)
in Appendix 1, File Formats, A.1 Rawfile Format.

The Xyce rawfile format is very similar and described
[here](https://xyce.sandia.gov/files/xyce/Reference_Guide.pdf#section.8.2).

The function outputs a header, the times/frequencies, the currents, and the
voltages. The voltage and current arrays have dimensions nVoltages by nPoints
and nCurrents by nPoints.

"""
function spice_raw_load(filename)

    #open a file handle
    io = open(filename)

    # the contents of the header
    title = ""
    date = ""
    plotname = ""
    flags = ""
    nvariables = 0
    npoints = 0
    command = ""
    option = ""

    #other variables
    headerlength = 0
    filetype = ""
    ncurrents = 0
    nvoltages = 0
    ntimes = 0
    # separator = " "

    #loop over the contents of the header
    i=0
    # while true
    while !eof(io)

        i+=1
        
        line = readline(io)
        linesplit = split(line,":",limit=2)

        if length(linesplit) == 2
            linename = linesplit[1]
            linevalue = strip(linesplit[2])
        else
            throw(ArgumentError("Line doesn't have the correct format."))
        end

        if linename == "Title"
            title = linevalue
        elseif linename == "Date"
            date = linevalue
        elseif linename == "Plotname"
            plotname = linevalue
        elseif linename == "Flags"
            flags = linevalue
        elseif linename == "No. Variables"
            nvariables =  parse(Int,linevalue)
        elseif linename == "No. Points"
            npoints = parse(Int,linevalue)
        elseif linename == "Command"
            command = linevalue
        elseif linename == "Option"
            option = linevalue
        elseif linename == "Variables"
            headerlength = i
            break
        end
    end

    header = SpiceRawHeader(title, date, plotname, flags, nvariables, npoints,
        command, option)

    #array of strings to contain all of the variable names
    variables =  Dict{String,Vector{String}}()
    indices = Dict{String,Vector{Int}}()

    #loop over the variable names
    i = 0
    # while true
    while !eof(io)

        i+=1
        # read the line, remove leading and trailing whitespace
        line = strip(readline(io))

        # println(line)

        # break if we reach the end of the variables section
        if line == "Binary:" || line == "Values:"
            filetype = line
            break
        end

        # split the line and assign the variable to a dictionary
        splitline = split(line,r"\s+")

        if length(splitline) > 3
            @warn "Variable line has additional parameters which we are ignoring."
        end 
        # variables[i-headerlength] = splitline[3]
        if !haskey(variables,splitline[3])
            variables[splitline[3]] = Vector{String}(undef,0)
            indices[splitline[3]] = Vector{Int}(undef,0)
        end

        # store the variable and the index at which it occurs
        push!(variables[splitline[3]],splitline[2])
        # push!(indices[splitline[3]],i)
        push!(indices[splitline[3]],i)

    end 

    # read the data
    if filetype == "Binary:"

        # make an empty array for the binary data
        if flags == "real"
            sf = Array{Float64}(undef,nvariables,npoints)
            values = Dict{String,Matrix{Float64}}()
        elseif flags =="complex"
            sf = Array{Complex{Float64}}(undef,nvariables,npoints)
            values = Dict{String,Matrix{Complex{Float64}}}()
        else
            throw(ArgumentError("Unknown flag."))
        end

        # read the binary data
        read!(io,sf)

        #close the file
        close(io)
    else
        throw(ArgumentError("This function only handles Binary files not ASCII."))
    end
    
    # sort the labels. voltages such as  "V(1)","V(10)","V(100)"."V(101)"
    sortperms = calcspicesortperms(variables)

    # use the sort permutation to sort the rest of the data
    # values = Dict()
    for (label,sp) in sortperms
        values[label] = sf[indices[label][sp],:]
        variables[label] = variables[label][sp]
    end

    return SpiceRaw(header, variables, values)
end

"""
    calcspicesortperms(variabledict::Dict{String,Vector{String}})

Calculate the sortperms which will sort the variable and node names.
"""
function calcspicesortperms(variabledict::Dict{String,Vector{String}})

    sortperms = Dict{String,Vector{Int}}()

    for (label,variables) in variabledict
        sortvariables = Dict()
        for variable in variables
            key, val = parsespicevariable(variable)
            if !haskey(sortvariables,key)
                sortvariables[key] = typeof(val)[]
            end
            push!(sortvariables[key],val)
        end
        #loop over the sorted outer arrays
        sp = Int[]
        # for (key,val) in sort(collect(sortvariables), by=x->x)
        for (key,val) in sort(collect(sortvariables), by=identity)
            #sort the dictionary
            p = sortperm(val)
            sp = vcat(sp,p .+ length(sp))
        end

        sortperms[label] = sp
    end

    return sortperms
end

"""
    parsespicevariable(variable::String)

Parse a variable name string into the variable name and node number. Will this
work with arbitrary node strings?

# Examples
```jldoctest
julia> JosephsonCircuits.parsespicevariable("V1(5)")
("V1", 5)

julia> JosephsonCircuits.parsespicevariable("V1")
("V", 1)

julia> JosephsonCircuits.parsespicevariable("V-1")
("V", 1)

julia> JosephsonCircuits.parsespicevariable("frequency")
("frequency", "frequency")
```
"""
function parsespicevariable(variable::String)

    s = variable
    m1 = match(r"^\w+",s)
    if isnothing(m1)
        throw(ArgumentError("No match found."))
    end

    m2 = match(r"\d+",s)
    m3 = match(r"\d+",s[m1.offset+length(m1.match):end])

    if !isnothing(m3)
        #if there is a separate number, use that
        key = m1.match
        val = parse(Int,m3.match)
    elseif !isnothing(m2)
        key = m1.match[1:m2.offset-1]
        #otherwise if there is a symbol separated by a number 
        val = parse(Int,m2.match)
    else
        key = m1.match
        val = m1.match
    end

    return key,val
  end