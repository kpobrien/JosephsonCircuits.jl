"""
    touchstone_load(filename)

Read a file in the Touchstone format. Standards compliant. This is the 1.1 spec:
https://ibis.org/connector/touchstone_spec11.pdf
and the 2.0 spec:
https://ibis.org/touchstone_ver2.0/touchstone_ver2_0.pdf

Outputs un-normalized network parameters with frequency units of Hz. 
"""
function touchstone_load(filename)

    # Don't enforce any particular file extension or try to infer the number of
    # ports from the extension.

    #open a file handle
    io = open(filename)

    return touchstone_parse(io)
end

function touchstone_parse(io)

    # define a dictionary for the contents of the header
    header = Dict()

    #other variables
    readoptionline = true

    #loop over the contents of the header
    while !eof(io)
        
        # read the line
        line = readline(io)

        # discard anything after the comment character !
        m = findfirst('!',line)
        if m !== nothing
            if !haskey(header,:comments)
                header[:comments] = []
            end
            push!(header[:comments],line[m+1:end])
            line = line[1:m-1]
        end

        # parse the line
        if length(line) == 0
            # don't bother parsing the line if it has zero length
            nothing
        elseif (m = match(r"^\[version\]",lowercase(line))) !== nothing
            # parse the [Version] line
            if haskey(header,:version)
                error("Only one [Version] keyword allowed:\n$(line)")
            else
                header[:version]=parse(Float64,line[m.offset+length(m.match):end])
            end

        elseif match(r"^#",lowercase(line)) !== nothing && readoptionline
            # parse the option line, if it is the first option line, otherwise ignore
            # per the spec "additional option lines after the first one shall be ignored"

            # In summary, the option line should read:
            # For 1-port files: # [Hz|kHz|MHz|GHz] [S|Y|Z] [DB|MA|RI] [R n]
            # For 2-port files: # [Hz|kHz|MHz|GHz] [S|Y|Z|G|H] [DB|MA|RI] [R n]
            # For 3-port and beyond files: # [Hz|kHz|MHz|GHz] [S|Y|Z] [DB|MA|RI] [R n]
            # For mixed-mode files: # [Hz|kHz|MHz|GHz] [S|Y|Z] [DB|MA|RI] [R n]

            # All of these parameters are optional. A line containing only # 
            # is a valid option line. I'm not checking for order although I
            # probably should. 
            splitline = split(strip(lowercase(line)),r"\s+")[2:end]
            i=1
            for s in splitline
                if s == "hz"
                    header[:frequencyunit] = 1.0e0
                elseif s == "khz"
                    header[:frequencyunit] = 1.0e3
                elseif s == "mhz"
                    header[:frequencyunit] = 1.0e6
                elseif s == "ghz"
                    header[:frequencyunit] = 1.0e9
                elseif s == "s"
                    header[:parameter] = "S"
                elseif s == "y"
                    header[:parameter] = "Y"
                elseif s == "z"
                    header[:parameter] = "Z"
                elseif s == "g"
                    header[:parameter] = "G"
                elseif s == "h"
                    header[:parameter] = "H"
                elseif s == "db"
                    header[:format] = "DB"
                elseif s == "ma"
                    header[:format] = "MA"
                elseif s == "ri"
                    header[:format] = "RI"
                elseif s == "r"
                    # i should switch this to parsing the next entry
                    header[:R] = parse.(Float64,splitline[i+1])
                    break
                else
                    error("Unknown parameter in options line:\n$(line)")
                end
                i+=1
            end

            # Rules for Version 1.0 Files:
            # For Version 1.0 files, the option line shall precede any data lines 
            # and shall be the first non-comment, non-blank line in the file.   

            # Rules for Version 2.0 Files:
            # For Version 2.0 files, the option line shall follow the [Version] 
            # keyword and argument and precede the [Number of Ports] keyword and 
            # argument.

            if  !haskey(header,:version)
                #this is a version 1.0 or 1.1 file. break the loop and move to the data
                # header["version"] = 1.0
                break
            elseif haskey(header,:version) && !haskey(header,:nports)
                # this is a version 2.0 file and the option line preceeds the 
                # number of ports keyword as required. continue reading the header.
                nothing
            else
                error("Option line position not valid for Touchstone v 1.0 or v 2.0:\n$(line)")
            end


            # set this to false so that any subsequent option lines are ignored
            # during header parsing
            readoptionline = false      

        elseif match(r"^#",lowercase(line)) !== nothing && !readoptionline
            # the spec says a second options line should be ignored, but the 
            # golden parser tsck2 throws an error. 
            error("Error: Invalid secondary options line:\n$(line)")
        elseif (m = match(r"^\[number of ports\]",lowercase(line))) !== nothing
            # parse the [Number of Ports] line
            if haskey(header,:numberofports)
                error("Error: Only one [Number of Ports] keyword allowed:\n$(line)")
            else
                header[:numberofports]=parse(Int64,line[m.offset+length(m.match):end])
            end
        
        elseif (m = match(r"^\[two-port data order\]",lowercase(line))) !== nothing
            # parse the [Two-Port Data Order] line. required if number of ports is 2
            if header[:numberofports] == 2
                if haskey(header,:twoportdataorder)
                    error("Error: Only one [Two-Port Data Order] line allowed:\n$(line)")
                else
                    if strip(line[m.offset+length(m.match):end]) == "12_21"
                        header[:twoportdataorder] = "12_21"
                    elseif strip(line[m.offset+length(m.match):end]) == "21_12"
                        header[:twoportdataorder] = "21_12"
                    else
                        error("Error: Unknown [Two-Port Data Order] parameter:\n$(line)")
                    end
                end
            else
                error("Error: [Two-Port Data Order] is only allowed if [Number of Ports] is 2:\n$(line)")
            end

        elseif (m = match(r"^\[begin information\]",lowercase(line))) !== nothing
            if !haskey(header,:numberofports)
                error("Error: [Number of Ports] must be before [Begin Information]")
            end

            header[:information] = []

            #loop over the contents of the header
            while !eof(io)
                # read the line
                line = readline(io)

                # end the loop once we reach the end information line
                if (m = match(r"^\[end information\]",lowercase(line))) !== nothing
                    break
                end

                # discard anything after the comment symbol !
                m = findfirst('!',line)
                if m !== nothing
                    if !haskey(header,:comments)
                        header[:comments] = []
                    end
                    push!(header[:comments],line[m+1:end])
                    line = line[1:m-1]
                end

                # store the contents of the line in the information array
                if length(line) > 0
                    push!(header[:information],line)
                end

            end

        elseif (m = match(r"^\[number of frequencies\]",lowercase(line))) !== nothing
            # parse the [Number of Frequencies] line
            header[:numberoffrequencies]=parse(Int64,line[m.offset+length(m.match):end])
            if header[:numberoffrequencies] <= 0
                error("Error: Number of frequencies must be an integer greater than zero:\n$(line)")
            end

        elseif (m = match(r"^\[number of noise frequencies\]",lowercase(line))) !== nothing
            # parse the [Number of Noise Frequencies] line
            header[:numberofnoisefrequencies]=parse(Int64,line[m.offset+length(m.match):end])
            if header[:numberofnoisefrequencies] <= 0
                error("Error: Number of noise frequencies must be an integer greater than zero:\n$(line)")
            end
        
        elseif (m = match(r"^\[reference\]",lowercase(line))) !== nothing
            # parse the [Reference] line
            if haskey(header,:reference)
                error("Only one [Reference] keyword allowed:\n$(line)")
            end

            # if there is a number on this line, parse the line
            if match(r"\d",lowercase(line[m.offset+length(m.match):end])) != nothing
                header[:reference]=parse.(Float64,split(strip(line[m.offset+length(m.match):end]),r"\s+"))
            end

            # check if we have all of the references and take more data if we do not.
            if !haskey(header,:reference) || (length(header[:reference]) < header[:numberofports])

                #loop over the contents of the header
                while !eof(io)

                    # read the line
                    line = readline(io)

                    # discard anything after the comment symbol !
                    m = findfirst('!',line)
                    if m !== nothing
                        if !haskey(header,:comments)
                            header[:comments] = []
                        end
                        push!(header[:comments],line[m+1:end])
                        line = line[1:m-1]
                    end

                    if length(line) > 0
                        if haskey(header,:reference)
                            header[:reference] = vcat(
                                header[:reference],
                                parse.(Float64,split(strip(line),r"\s+"))
                            )
                        else
                            header[:reference] = parse.(Float64,split(strip(line),r"\s+"))
                        end
                        if length(header[:reference]) == header[:numberofports]
                            break
                        elseif length(header[:reference]) > header[:numberofports]
                            error("Too many values on [Reference] line:\n$(line)")
                        end
                    end

                end

            end
        
        elseif (m = match(r"^\[matrix format\]",lowercase(line))) !== nothing
            # parse the [Matrix Format] line
            format = strip(lowercase(line[m.offset+length(m.match):end]))
            if format == "full"
                header[:matrixformat] = "Full"
            elseif format == "lower"
                header[:matrixformat] = "Lower"
            elseif format == "upper"
                header[:matrixformat] = "Upper"                
            else
                error("Unknown format:\n$(line)")
            end         
        elseif (m = match(r"^\[mixed-mode order\]",lowercase(line))) !== nothing
            # parse the [Mixed-Mode Order] line
            splitline = split(strip(lowercase(line[m.offset+length(m.match):end])),r"\s+");

            header[:mixedmodeorder] = [] 

            for l in splitline
                splitline2 = parse.(Int64,split(strip(lowercase(l[2:end])),r","));
                
                push!(header[:mixedmodeorder],(uppercase(l[1]),splitline2))
            end

        elseif line == "[Network Data]"
            break
        else
            # println(line)
            error("Unknown line type:\n$(line)")
        end
    end

    # do some checks on the header and assign default parameters
    # if haskey(header,"twoportdataorder")
    # need separate checks for version 1.0 and version 2.0

    #here are some checks
        # In summary, the option line should read:
    # For 1-port files: # [Hz|kHz|MHz|GHz] [S|Y|Z] [DB|MA|RI] [R n]
    # For 2-port files: # [Hz|kHz|MHz|GHz] [S|Y|Z|G|H] [DB|MA|RI] [R n]
    # For 3-port and beyond files: # [Hz|kHz|MHz|GHz] [S|Y|Z] [DB|MA|RI] [R n]
    # For mixed-mode files: # [Hz|kHz|MHz|GHz] [S|Y|Z] [DB|MA|RI] [R n]


    ### checks for both v1 and v2

    # check that the required keywords are present


    ### checks for v1 only, use if else to only select one of these two. 
    #check that v1 files only have the option line.

    ### checks for v2 only

    if haskey(header,:version) 
        if header[:version] != 2.0
            error("Error: Only version 1 (without a version keyword) or 2 are supported.")
        end
    end


    if haskey(header,:version) && header[:version] >= 2.0
        if haskey(header,:numberofports) && header[:numberofports] == 2 && !haskey(header,:twoportdataorder)
            error("[Two-Port Data Order] is required when the number of ports is 2.")
        end
    end

    # a version 1.0 or 1.1 file will only have the option line. 

    # set some default parameters

    # default format is MA (magnitude angle)
    if !haskey(header,:format)
        header[:format] = "MA"
    end

    # default frequency is GHz
    if !haskey(header,:frequencyunit)
        header[:frequencyunit] = 1.0e9
    end

    # default network parameter is S
    if !haskey(header,:parameter)
        header[:parameter] = "S"
    end

    # default reference resistance is 50 ohms. Note: this is overridden by [Reference]
    if !haskey(header,:R)
        header[:R] = 50.0
    end


    # finished reading the header, now read the network and noise data

    networkdata = Array{Float64}(undef,0)
    noisedata = Array{Float64}(undef,0)

    ndatalines = 0
    nvals = 0
    freq = 0.0
    nfreq = 0
    nnoisefreq = 0

    while !eof(io)
        
        # read the line
        line = readline(io)

        # discard anything after the comment symbol !
        if ( m = findfirst("!",line)) !== nothing
            if m[1] == 1
                # i should make a different dictionary for this
                if !haskey(header,:comments)
                    header[:comments] = []
                end
                push!(header[:comments],line)
            end
            line = line[1:m[1]-1]
        end

        if length(line) < 1
            #skip any lines that don't have anything
            nothing
        elseif match(r"^#",lowercase(line)) !== nothing
            # skip any additional option lines. these could occur in v1 files
            # after the first option line. 
            nothing 
        elseif lowercase(strip(line)) == lowercase("[End]")
            break
        elseif lowercase(strip(line)) == lowercase("[Noise Data]")
            ndatalines = 0
            nvals = 0
            freq = 0.0
            nnoisefreq = 0
            break
        else
            # parse the network data lines
            vals = parse.(Float64,split(strip(line),r"\s+"))

            # check if the current line is a new frequency or not
            # println(length(vals))

            # println(vals)

            # for v2 files, we know the number of ports so we can read the correct
            # number of data points even if spread across multiple lines. 
            # for v1 files, we do not know the number of ports, so we have to read
            # we find a different frequency or until we reach the end of the file
            # either will give us the number of ports. 

            if isodd(length(vals))
                # this line is a new frequency

                if vals[1] < freq
                    # if the new frequency is less than the old frequency, then 
                    # we have moved on to noise data. add that value, break the loop
                    # and move to the noise data reading loop. 
                    for val in vals
                        push!(noisedata,val)
                    end
                    # check the length to see if it is consistent with noise data
                    if length(vals) != 5
                        error("Noise data lines must have 5 entries:\n$(line)")
                    end
                    ndatalines = 1
                    freq = vals[1]
                    nnoisefreq = 1
                    nvals=length(vals)
                   break

                else
                    ndatalines +=1

                    if ndatalines > 1
                        # i think i need a new formula here if the matrix type is upper or lower?
                        # if the matrix is a full matrix, it will have nvals = 2*nports^2+1 values (two values per element
                        # plus a frequency). so nports = (nvals-1)/2. if the matrix is an upper or lower
                        # triangular matrix then there are nvals = nports*(nports+1) +1  or
                        # nports = (sqrt(4*nvals-3)-1)/2
                        # 
                        if !haskey(header,:matrixformat) || header[:matrixformat] == "Full"
                            nports = convert(Int64,sqrt((nvals -1)/2))
                        else
                            nports = convert(Int64,(sqrt(4*nvals-3)-1)/2)
                        end

                        if !haskey(header,:numberofports)
                            header[:numberofports] = nports
                        elseif header[:numberofports] != nports
                            error("Number of ports in file and in header not consistent:\n$(line)")
                        end
                    end

                    # if this is the first frequency or the frequency is larger
                    # than the previous one, increment the number of frequencies
                    # counter
                    freq = vals[1]
                    nfreq+=1

                    # add the values to the network data array
                    nvals = length(vals)
                    for val in vals
                        push!(networkdata,val)
                    end
                end

            else
                # this line is a continuation of data for a previous frequency
                ndatalines+=1
                nvals +=length(vals)
                # add the values to the network data array
                for val in vals
                    push!(networkdata,val)
                end
            end
        end
    end

    # check that the number of frequencies are consistent. 
    if haskey(header,:version) && header[:version] == 2.0
        # for v 2.0 files
        if header[:numberoffrequencies] != nfreq
            # println(nfreq)
            # println(header["numberoffrequencies"] )
            error("Number of frequencies from header is not equal to number of frequencies in network data.")
        end
    else
        #for v 1.0 files, we have no other frequency information. 
        header[:numberoffrequencies] = nfreq
    end

    #it's a bit complicated because i might have read in the first line of noise 
    # data in a v1 file. 
    while !eof(io)
        
        # read the line
        line = readline(io)

        # discard anything after the comment symbol !
        if ( m = findfirst("!",line)) !== nothing
            if m[1] == 1
                # i should make a different dictionary for this
                if !haskey(header,:comments)
                    header[:comments] = []
                end
                push!(header[:comments],line)
            end
            line = line[1:m[1]-1]
        end

        if length(line) < 1
            #skip any lines that don't have anything
            nothing
        elseif match(r"^#",lowercase(line)) !== nothing
            # skip any additional option lines. these could occur in v1 files
            # after the first option line. 
            nothing 
        elseif lowercase(strip(line)) == lowercase("[End]")
            break
        elseif lowercase(strip(line)) == lowercase("[Noise Data]")
            error("Only one [Noise Data] keyword allowed.")
        else
            # parse the noise data lines
            vals = parse.(Float64,split(strip(line),r"\s+"))

            if length(vals) != 5
                error("Noise data lines must have 5 entries.")
            end

            # add the values to the noise data array
            nvals = length(vals)
            for val in vals
                push!(noisedata,val)
            end

            if vals[1] < freq
                # if the new frequency is less than the old frequency, then 
                # we have moved on to noise data. add that value, break the loop
                # and move to the noise data reading loop. 
                error("Frequencies descending in noise data")
            else
                freq = vals[1]
                nnoisefreq+=1
            end

        end
    end

    if !haskey(header,:numberofports)

        if !haskey(header,:matrixformat) || header[:matrixformat] == "Full"
            nports = convert(Int64,sqrt((nvals -1)/2))
        else
            nports = convert(Int64,(sqrt(4*nvals-3)-1)/2)
        end

        if !haskey(header,:numberofports)
            header[:numberofports] = nports
            # println(header["numberofports"])
        elseif header[:numberofports] != nports
            error("Number in file and in header not consistent")
        end

    end

    # #check that the length of the network data agrees with the expected number of frequencies
    # println((2*header["numberofports"]*header["numberofports"]+1)*header["numberoffrequencies"])
    # println(length(networkdata))


    #make an empty array for the network data
    N = Array{Complex{Float64}}(undef,header[:numberofports],header[:numberofports],header[:numberoffrequencies])
    f = Array{Float64}(undef,header[:numberoffrequencies])

    # for i = 1:header["numberoffrequencies"]
    #     f[i] = networkdata[(2*header["numberofports"]^2+1)*(i-1)+1]
    # end

    if haskey(header,:twoportdataorder) && haskey(header,:matrixformat)
        indices = matrixindices(header[:numberofports],header[:matrixformat],header[:twoportdataorder])
    elseif haskey(header,:twoportdataorder)
        indices = matrixindices(header[:numberofports],"Full",header[:twoportdataorder])
    elseif haskey(header,:matrixformat)
        indices = matrixindices(header[:numberofports],header[:matrixformat])
    else
        if !haskey(header,:version) && header[:numberofports] == 2
            # 21_12 is the default for v 1 files when there are two ports. 
            indices = matrixindices(header[:numberofports],"Full","21_12")
        else        
            indices = matrixindices(header[:numberofports],"Full")
        end
    end


    if (header[:parameter] == "Z" || header[:parameter] == "Y" || header[:parameter] == "G" || header[:parameter] == "H") && !haskey(header,:version)
        normalization = header[:R]
    elseif (header[:parameter] == "Z" || header[:parameter] == "Y" || header[:parameter] == "G" || header[:parameter] == "H") && header[:version] == 2.0
        normalization = 1
    elseif header[:parameter] == "S"
        normalization = 1
    else
        error("Error: Unknown format or version")
    end

    f .= networkdata[1:(2*length(indices)+1):end]

    for n = 1:length(indices)
        if header[:format] == "DB"
            N[indices[n],:] .= 10 .^((1/20).*networkdata[2*n:(2*length(indices)+1):end].*normalization)
            N[indices[n],:]  .*= exp.((pi/180*im).*networkdata[2*n+1:(2*length(indices)+1):end])
        elseif header[:format] == "MA"
            N[indices[n],:] .= networkdata[2*n:(2*length(indices)+1):end].*normalization
            N[indices[n],:] .*= exp.((pi/180*im).*networkdata[2*n+1:(2*length(indices)+1):end])
        elseif header[:format] == "RI"
            N[indices[n],:] .= networkdata[2*n:(2*length(indices)+1):end].*normalization
            N[indices[n],:] .+= im.*networkdata[2*n+1:(2*length(indices)+1):end].*normalization
        end

        # if the format is upper or lower, the matrix is assumed to be symmetric
        # so fill in the other parts
        if haskey(header,:matrixformat) && header[:matrixformat] != "Full"
            if indices[n][1] != indices[n][2]
                N[indices[n][2],indices[n][1],:] .= N[indices[n][1],indices[n][2],:]
            end
        end
    end


    if !haskey(header,:reference)
        header[:reference] = header[:R].*ones(Float64,header[:numberofports])
    end

    close(io)


    # networkdata = (Symbol(header["parameter"]) = N, Symbol("frequencies") = f)

    key = (Symbol(header[:parameter]), Symbol("f"))
    val = (N, f*header[:frequencyunit])    

    networkdata = (; zip(key, val)...) 

    # is there any reason i shouldn't use symbols everywhere?

    header = (; zip((keys(header)), values(header))...)

    return (header=header,networkdata=networkdata,noisedata=noisedata)
end

"""
    matrixindices(nports,format)

"""
function matrixindices(nports,format)
    return matrixindices(nports,format,"12_21")
end

function matrixindices(nports,format,twoportdataorder)
    indices = []
    printflag = false
    
    if format == "Lower"
        ncol = 1
        ntotal = nports*(1+nports) ÷ 2
    elseif format == "Upper"
        ncol = nports
        ntotal = nports*(1+nports) ÷ 2
    else
        ntotal = nports^2
        ncol = nports
    end    

    i = 1
    j = 1
    nold = 1
    for n = 1:ntotal
        # this is the row index.
        # increment it when we have written ncol
        # data points since the last one. 
        if n - nold == ncol
            i+=1
            nold = n
            if format == "Lower"
                ncol+= 1
            elseif format == "Upper"
                ncol-= 1
            end
        end

        #j = mod(n-1,ncol)+1 # this is the column index
        # this is the column index
        j = n - nold + 1

        if format == "Upper"
            if printflag
                print("$(i)$(j+nrow-ncol) ")
            end
            push!(indices,CartesianIndex(i,j+nrow-ncol))
        else
            if printflag
                print("$(i)$(j) ")
            end
            push!(indices,CartesianIndex(i,j))

        end
        if j == ncol
            if printflag
                print("\n")
            end
            if format == "Upper"
                for k = 1:nrow-ncol+1
                    print("   ")
                end
            end
        end
    end
    
    if twoportdataorder == "21_12" && nports == 2
        for n = 1:length(indices)
            if indices[n][1] != indices[n][2]
                indices[n] = CartesianIndex(indices[n][2],indices[n][1])
            end
        end
    elseif twoportdataorder == "21_12"
        error("Two port data order = 21_12 is only allowed if the number of ports is two")
    elseif twoportdataorder == "12_21"
        nothing
    else
        error("Unknown two port data order string.")
    end
    return indices
end


"""
    touchstone_save(frequencies::AbstractArray,N::AbstractArray,
    filename::String;version=1.0,reference=[50.0,50.0],R = 50.0,format="RI",
    comments=[""],twoportdataorder="12_21",matrixformat="Full",frequencyunit="Hz")

Write a file in the Touchstone format. Standards compliant except does not
support writing noise data. This is the 1.1 spec:
https://ibis.org/connector/touchstone_spec11.pdf
and the 2.0 spec:
https://ibis.org/touchstone_ver2.0/touchstone_ver2_0.pdf

"""
function touchstone_save(frequencies::AbstractArray,N::AbstractArray,
    filename::String;version=1.0,reference=[50.0,50.0],R = 50.0,format="RI",
    comments=[""],twoportdataorder="12_21",matrixformat="Full",frequencyunit="Hz")


    header = Dict()

    header[:version] = version
    header[:reference] = reference
    header[:format] = format
    header[:comments] = comments
    header[:twoportdataorder] = twoportdataorder
    header[:matrixformat] = matrixformat
    header[:frequencyunit] = frequencyunit
    header[:R] = R

    # header[:version] = 2.0
    # header[:reference] = [50.0, 50.0]
    # header[:format] = "RI"
    # header[:comments] = [""]
    # header[:twoportdataorder] = "12_21"
    # header[:matrixformat] = "Full"
    # header[:frequencyunit] = "Hz"


    # set the reference impedance to 50.0 Ohms if no impedance given.
    if !haskey(header,:R)
        header[:R] = 50.0
    end

    # set the frequency unit to GHz if one is not given
    if !haskey(header,:frequencyunit)
        header[:frequencyunit] = "GHz"
    elseif lowercase(header[:frequencyunit]) == lowercase("Hz")
        header[:frequencyunit] = "Hz"
    elseif lowercase(header[:frequencyunit]) == lowercase("kHz")
        header[:frequencyunit] = "kHz"
    elseif lowercase(header[:frequencyunit]) == lowercase("MHz")
        header[:frequencyunit] = "MHz"
    elseif lowercase(header[:frequencyunit]) == lowercase("GHz")
        header[:frequencyunit] = "GHz"
    else
        error("Error: Unknown option for the frequency unit.")
    end

    # set the parameter to S if none is given
    if !haskey(header,:parameter)
        header[:parameter] = "S"
    elseif lowercase(header[:parameter]) == lowercase("S")
        header[:parameter] = "S"
    elseif lowercase(header[:parameter]) == lowercase("Y")
        header[:parameter] = "Y"
    elseif lowercase(header[:parameter]) == lowercase("Z")
        header[:parameter] = "Z"
    elseif lowercase(header[:parameter]) == lowercase("H")
        header[:parameter] = "H"
    elseif lowercase(header[:parameter]) == lowercase("G")
        header[:parameter] = "G"
    else
        error("Error: Unknown option for the parameter.")
    end

    # set the default format to MA if none is given
    if !haskey(header,:format)
        header[:format] = "MA"
    elseif lowercase(header[:format]) == lowercase("DB")
        header[:format] = "DB"
    elseif lowercase(header[:format]) == lowercase("MA")
        header[:format] = "MA"
    elseif lowercase(header[:format]) == lowercase("RI")
        header[:format] = "RI"
    else
        error("Error: Unknown option for the format.")
    end

    # check that the network data doesn't have too many dimensions
    if size(N)[1] == size(N)[2]
        header[:numberofports] = size(N)[1]
    else
        error("Error: Network data arrays are not square.")
    end

    if length(size(N)) == 2
        header[:numberoffrequencies] = 1
    elseif length(size(N)) == 3
        header[:numberoffrequencies] = size(N)[3]
    else
        error("Error: Network data array has too many dimensions.")
    end

    # check that the dimensions of the frequency data are consistent with the
    # inferred number of frequencies from the network parameters
    if length(frequencies) != header[:numberoffrequencies]
        error("Error: Number of frequencies from length of frequencies variable and from network data not consistent.")
    end

    # check the filename. it the file doesn't have a .ts or .sNp extension,
    # then add the .sNp extension, where N is the number of ports. 
    if (m = match(r"\.ts$|\.s\dp$",lowercase(filename))) != nothing
        if header[:version] == 1.0 || header[:version] == 1.1  || !haskey(header,:version)
            if m.match != ".s$(header[:numberofports])p"
                println("Warning: Extension of $(m.match) is not the recommended extension of .s$(header[:numberofports])p for a version 1.0 file with $(header[:numberofports]) ports.")
            end
        else
            if m.match != ".s$(header[:numberofports])p" && m.match != ".ts"
                println("Warning: Extension of $(m.match) is not the recommended extension of .ts or .s$(header[:numberofports])p for a file with $(header[:numberofports]) ports.")
            end
        end

    else
        println("Warning: Adding extension of .s$(header[:numberofports])p")
        filename = filename*".s$(header[:numberofports])p"
    end

    # if reference impedance is given, make sure there is one per port
    if haskey(header,:reference)
        if length(header[:reference]) !=  header[:numberofports]
            error("Error: Number of reference impedances not equal to number of ports.")
        end
    end

    # if we are writing a version 1.0 file
    if header[:version] == 1.0 || header[:version] == 1.1  || !haskey(header,:version)

        # open the file

        open(filename,"w") do io

            # write the comments
            for comment in header[:comments]
                write(io,"!$(comment)\n")
            end

            # write the option line
            write(io,"# $(header[:frequencyunit]) $(header[:parameter]) $(header[:format]) R $(header[:R])\n")

            if header[:numberofports] == 2
                indices = matrixindices(header[:numberofports],"Full","21_12")
            else
                indices = matrixindices(header[:numberofports],"Full")
            end

            # write a comment with the ports
            write(io,"! freq ")
            for j = 1:length(indices)
                if header[:format] == "RI"
                    write(io,"Re$(header[:parameter])$(indices[j][1])$(indices[j][2]) Im$(header[:parameter])$(indices[j][1])$(indices[j][2]) ")
                elseif header[:format] == "MA"
                    write(io,"mag$(header[:parameter])$(indices[j][1])$(indices[j][2]) ang$(header[:parameter])$(indices[j][1])$(indices[j][2]) ")
                elseif header[:format] == "DB"
                    write(io,"logmag$(header[:parameter])$(indices[j][1])$(indices[j][2]) ang$(header[:parameter])$(indices[j][1])$(indices[j][2]) ")
                else
                    error("Error: Unknown format")
                end
            end
            write(io,"\n")

            # normalize Z, Y, G, or H by the reference impedance per the spec. 
            if header[:parameter] == "Z" || header[:parameter] == "Y" || header[:parameter] == "G" || header[:parameter] == "H"
                normalization = header[:R]
            elseif header[:parameter] == "S"
                normalization = 1
            else
                error("Error: Unknown format")
            end

            # write the network data
            for i = 1:length(frequencies)
                write(io,"$(frequencies[i])")
                k=1
                for j = 1:length(indices)

                    # do different things for each format
                    if header[:format] == "RI"
                        s1 = "$(real(N[indices[j],i]/normalization) )"
                        s2 = "$(imag(N[indices[j],i]/normalization) )"

                    elseif header[:format] == "MA"
                        s1 = "$(abs(N[indices[j],i]/normalization) )"
                        s2 = "$(180/pi*angle(N[indices[j],i]/normalization) )"

                    elseif header[:format] == "DB"
                        s1 = "$(10*log10(abs2(N[indices[j],i]/normalization)) )"
                        s2 = "$(180/pi*angle(N[indices[j],i]/normalization) )"
                    else
                        error("Error: Unknown format")
                    end
                    # write the strings
                    write(io," "*s1)
                    write(io," "*s2)

                    # only 4 pairs allowed per line
                    # the convention is to || mod(j,header[:numberofports]) == 0)
                    if  k  == 4  && j != length(indices)
                        k = 1
                        write(io,"\n")
                    elseif header[:numberofports] == 2
                        nothing
                    elseif mod(j,header[:numberofports]) == 0  && j != length(indices)
                        k = 0
                        write(io,"\n")
                    end
                    k+=1
                end
                write(io,"\n")
            end

        end

        # write the noise data. not implemented yet


    elseif header[:version] == 2.0


        # open the file

        open(filename,"w") do io

            # write the comments
            for comment in header[:comments]
                write(io,"!$(comment)\n")
            end

            # write the version keyword
            write(io,"[Version] 2.0\n")

            # write the option line
            write(io,"# $(header[:frequencyunit]) $(header[:parameter]) $(header[:format]) R $(header[:R])\n")

            # write the number of ports
            write(io,"[Number of Ports] $(header[:numberofports])\n")

            # write the two-port data order keyword if there are two ports
            if header[:numberofports] == 2
                write(io,"[Two-Port Data Order] $(header[:twoportdataorder])\n")
            end

            # write the number of frequencies
            write(io,"[Number of Frequencies] $(header[:numberoffrequencies])\n")

            if haskey(header,:numberofnoisefrequencies)
                # write the number of noise frequencies
                write(io,"[Number of Noise Frequencies] $(header[:numberofnoisefrequencies])\n")
            end

            # write the reference keyword
            write(io,"[Reference]")
            for r in header[:reference]
                write(io," $(real(r))")
            end

            write(io,"\n")


            # write the matrix format


            # write the network data keyword
            write(io,"[Network Data]\n")

            # calculate the indices
            if header[:numberofports] == 2
                indices = matrixindices(header[:numberofports],header[:matrixformat],header[:twoportdataorder])
            else
                indices = matrixindices(header[:numberofports],header[:matrixformat])
            end

            # write a comment with the ports
            write(io,"! freq ")
            for j = 1:length(indices)
                if header[:format] == "RI"
                    write(io,"Re$(header[:parameter])$(indices[j][1])$(indices[j][2]) Im$(header[:parameter])$(indices[j][1])$(indices[j][2]) ")
                elseif header[:format] == "MA"
                    write(io,"mag$(header[:parameter])$(indices[j][1])$(indices[j][2]) ang$(header[:parameter])$(indices[j][1])$(indices[j][2]) ")
                elseif header[:format] == "DB"
                    write(io,"logmag$(header[:parameter])$(indices[j][1])$(indices[j][2]) ang$(header[:parameter])$(indices[j][1])$(indices[j][2]) ")
                else
                    error("Error: Unknown format")
                end
            end
            write(io,"\n")

            # normalize Z, Y, G, or H by the reference impedance per the spec. 
            if header[:parameter] == "Z" || header[:parameter] == "Y" || header[:parameter] == "G" || header[:parameter] == "H"
                normalization = header[:R]
            elseif header[:parameter] == "S"
                normalization = 1
            else
                error("Error: Unknown format")
            end

            # write the network data
            for i = 1:length(frequencies)
                write(io,"$(frequencies[i])")
                k=1
                for j = 1:length(indices)

                    # do different things for each format
                    if header[:format] == "RI"
                        s1 = "$(real(N[indices[j],i]/normalization) )"
                        s2 = "$(imag(N[indices[j],i]/normalization) )"

                    elseif header[:format] == "MA"
                        s1 = "$(abs(N[indices[j],i]/normalization) )"
                        s2 = "$(180/pi*angle(N[indices[j],i]/normalization) )"

                    elseif header[:format] == "DB"
                        s1 = "$(10*log10(abs2(N[indices[j],i]/normalization)) )"
                        s2 = "$(180/pi*angle(N[indices[j],i]/normalization) )"
                    else
                        error("Error: Unknown format")
                    end
                    # write the strings
                    write(io," "*s1)
                    write(io," "*s2)

                    # only 4 pairs allowed per line
                    # the convention is to || mod(j,header[:numberofports]) == 0)
                    if  k  == 4  && j != length(indices)
                        k = 1
                        write(io,"\n")
                    elseif header[:numberofports] == 2
                        nothing
                    elseif mod(j,header[:numberofports]) == 0  && j != length(indices)
                        k = 0
                        write(io,"\n")
                    end
                    k+=1
                end
                write(io,"\n")
            end

            write(io,"[End]\n")


        end

        # write the noise data. not implemented yet


    else
        error("Error: unsupported version $(header[:version])")

    end
end


