
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

Read the binary rawfile of a WRSPICE or Xyce analysis, transient or
frequency domain. The file format is documented in the
[WRSPICE manual](http://www.srware.com/xictools/docs/wrsmanual-4.3.13.pdf)
in Appendix 1, File Formats, A.1 Rawfile Format; the Xyce rawfile format
is very similar and described
[here](https://xyce.sandia.gov/files/xyce/Reference_Guide.pdf#section.8.2).

Returns a [`SpiceRaw`](@ref): the header, and for each type of variable,
`V` for the node voltages, `S` for the time of a transient or `Hz` for
the frequency of an AC analysis, the names of the variables of that type
in `variables` and their values in `values`, a matrix with one row per
variable and one column per point, real or complex as the header's
flags say. The variables of each type are sorted by
[`calcspicesortperms`](@ref). A rawfile in the ASCII format is refused.
"""
function spice_raw_load(filename)
    header, variables, indices, sf = open(filename) do io

        # the contents of the header
        title = ""
        date = ""
        plotname = ""
        flags = ""
        nvariables = 0
        npoints = 0
        command = ""
        option = ""

        #loop over the contents of the header
        while !eof(io)

            line = readline(io)
            linesplit = split(line,":",limit=2)

            if length(linesplit) == 2
                linename = linesplit[1]
                linevalue = strip(linesplit[2])
            else
                throw(ArgumentError(lazy"Line doesn't have the correct format."))
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
                break
            end
        end

        header = SpiceRawHeader(title, date, plotname, flags, nvariables,
            npoints, command, option)

        # the variable names and the rows they occupy, grouped by type
        variables =  Dict{String,Vector{String}}()
        indices = Dict{String,Vector{Int}}()

        #loop over the variable names
        filetype = ""
        i = 0
        while !eof(io)

            i+=1
            # read the line, remove leading and trailing whitespace
            line = strip(readline(io))

            # break if we reach the end of the variables section
            if line == "Binary:" || line == "Values:"
                filetype = line
                break
            end

            # the index, the name and the type of the variable; the
            # constants an aborted WRSPICE run writes have no type
            splitline = split(line,r"\s+")

            if length(splitline) > 3
                @warn lazy"Variable line has additional parameters which we are ignoring."
            end
            type = length(splitline) >= 3 ? String(splitline[3]) : ""

            # store the variable and the index at which it occurs
            push!(get!(variables, type, String[]), splitline[2])
            push!(get!(indices, type, Int[]), i)
        end

        # read the data
        if filetype == "Binary:"
            if flags == "real"
                sf = Array{Float64}(undef,nvariables,npoints)
            elseif flags =="complex"
                sf = Array{Complex{Float64}}(undef,nvariables,npoints)
            else
                throw(ArgumentError(lazy"Unknown flag."))
            end
            read!(io,sf)
        else
            throw(ArgumentError(lazy"This function only handles Binary files not ASCII."))
        end
        return header, variables, indices, sf
    end

    # sort the labels. voltages such as  "V(1)","V(10)","V(100)"."V(101)"
    sortperms = calcspicesortperms(variables)

    # use the sort permutation to sort the rest of the data
    values = Dict{String,typeof(sf)}()
    for (label,sp) in sortperms
        values[label] = sf[indices[label][sp],:]
        variables[label] = variables[label][sp]
    end

    return SpiceRaw(header, variables, values)
end

"""
    calcspicesortperms(variabledict::Dict{String,Vector{String}})

The permutation which sorts the variables of each type of a rawfile, by
the name and the node [`parsespicevariable`](@ref) reads from each: by
name, and for a name by node, the numbered nodes in numerical order
before the nodes named by words, which keep their order. So a rawfile
may mix nodes named by numbers and by words, and `v(2)` comes before
`v(10)`.

# Examples
```jldoctest
julia> JosephsonCircuits.calcspicesortperms(Dict("V" => ["v(10)", "v(out)", "v(2)", "v(1)"]))["V"]
4-element Vector{Int64}:
 4
 3
 1
 2
```
"""
function calcspicesortperms(variabledict::Dict{String,Vector{String}})
    # numbers numerically first, then everything else by its text
    sortkey(v) = v isa Integer ? (0, Int(v), "") : (1, 0, string(v))
    sortperms = Dict{String,Vector{Int}}()
    for (label, variables) in variabledict
        parsed = map(parsespicevariable, variables)
        sortperms[label] = sortperm(parsed;
            by = p -> (string(first(p)), sortkey(last(p))))
    end
    return sortperms
end

"""
    parsespicevariable(variable::String)

The name and the node of a rawfile variable, which
[`calcspicesortperms`](@ref) sorts by. A number after the leading word of
the variable, `V1(5)` or `V-1`, is the node and the word the name; else a
number within the leading word splits it, `V1` being the name `V` and the
node `1`; and a variable without a number is its own name and node. A
variable which does not start with a word character is refused.

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
        throw(ArgumentError(lazy"No match found."))
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