
"""
    h,t,i,v = spice_raw_load(filename)

Parse the binary raw output file from WRSPICE or Xyce. 
Tested for transient analysis and frequency domain analysis.
The file format is documented in the WRSPICE manual in 
Appendix 1, File Formats, A.1 Rawfile Format
http://www.srware.com/xictools/docs/wrsmanual-4.3.13.pdf

The Xyce rawfile format is very similar and described in:
https://xyce.sandia.gov/files/xyce/Reference_Guide.pdf#section.8.2

The function outputs a header, the times/frequencies,
the currents, and the voltages. The voltage and current
arrays have dimensions nVoltages by nPoints and nCurrents by
nPoints. 

"""
function spice_raw_load(filename)

    #open a file handle
    io = open(filename)

    # define a dictionary for the contents of the header
    header = Dict(
        "title" => "",
        "date" => "",
        "plotname" => "",
        "flags" => "",
        "nvariables" => 0,
        "npoints" => 0,
        "command" => "",
        "option" => "")

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

        linename = split(line,":",limit=2)[1]

        if linename == "Title"
            header["title"] = split(line,": ",limit=2)[2]
        elseif linename == "Date"
            header["date"] = split(line,": ",limit=2)[2]
        elseif linename == "Plotname"
            header["plotname"] = split(line,": ",limit=2)[2]
        elseif linename == "Flags"
            header["flags"] = split(line,": ",limit=2)[2]
        elseif linename == "No. Variables"
            header["nvariables"] = parse(Int,split(line,": ",limit=2)[2])
        elseif linename == "No. Points"
            header["npoints"] = parse(Int,split(line,": ",limit=2)[2])
        elseif linename == "Command"
            header["command"] = split(line,": ",limit=2)[2]
        elseif linename == "Option"
            header["option"] = split(line,": ",limit=2)[2]
        elseif linename == "Variables"
            headerlength = i
            break
        end
    end

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
            @warn "Variable line has additional parameters which we are ignoring"
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
        if header["flags"] == "real"
            sf = Array{Float64}(undef,header["nvariables"],header["npoints"])
        elseif header["flags"] =="complex"
            sf = Array{Complex{Float64}}(undef,header["nvariables"],header["npoints"])
        else
            error("Unknown flag")
        end

        # read the binary data
        read!(io,sf)

        #close the file
        close(io)
    else
        error("This function only handles Binary files not ASCII")
    end
    
    # sort the labels. voltages such as  "V(1)","V(10)","V(100)"."V(101)"
    sortperms = calcspicesortperms(variables)

    # use the sort permutation to sort the rest of the data
    values = Dict()
    for (label,sp) in sortperms
        values[label] = sf[indices[label][sp],:]
        variables[label] = variables[label][sp]
    end

    return (variables=variables,values=values)

end

"""
    calcspicesortperms(variabledict::Dict{String,Vector{String}})

# Examples
```jldoctest
julia> JosephsonCircuits.calcspicesortperms(Dict("V" => ["v(1)", "v(2)", "v(3)"], "Hz" => ["frequency"]))
Dict{String, Vector{Int64}} with 2 entries:
  "V"  => [1, 2, 3]
  "Hz" => [1]
```
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
    m2 = match(r"\d+",s)
    m3 = match(r"\d+",s[m1.offset+length(m1.match):end])

    if m3 !== nothing
        #if there is a separate number, use that
        key = m1.match
        val = parse(Int,m3.match)
    elseif m2 !== nothing
        key = m1.match[1:m2.offset-1]
        #otherwise if there is a symbol separated by a number 
        val = parse(Int,m2.match)
    else
        key = m1.match
        val = m1.match
    end

    return key,val
  end