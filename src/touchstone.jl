"""
    TouchstoneFile(f::Vector{Float64},
        N::Array{Complex{Float64}},
        frequencyunit::String,
        parameter::String,
        format::String,
        R::Float64,
        version::Float64,
        numberofports::Int,
        twoportdataorder::String,
        numberoffrequencies::Int,
        numberofnoisefrequencies::Int,
        reference::Vector{Float64},
        information::Vector{String},
        matrixformat::String,
        mixedmodeorder::Vector{Tuple{Char, Vector{Int}}},
        comments::Vector{String},
        networkdata::Vector{Float64},
        noisedata::Vector{Float64})

A structure to hold the data contained in a Touchstone file. In most cases, 
the user will not generate the struct directly. Instead, they will load a 
Touchstone file with [`touchstone_load`](@ref), parse an IOStream or IOBuffer
with [`touchstone_parse`](@ref), or generate a TouchstoneFile struct with
[`touchstone_file`](@ref).
"""
struct TouchstoneFile
    f::Vector{Float64}
    N::Array{Complex{Float64}}
    frequencyunit::String
    parameter::String
    format::String
    R::Float64
    version::Float64
    numberofports::Int
    twoportdataorder::String
    numberoffrequencies::Int
    numberofnoisefrequencies::Int
    reference::Vector{Float64}
    information::Vector{String}
    matrixformat::String
    mixedmodeorder::Vector{Tuple{Char, Vector{Int}}}
    comments::Vector{String}
    networkdata::Vector{Float64}
    noisedata::Vector{Float64}
end

struct TouchstoneOptionLine
  frequencyunit::String
  parameter::String
  format::String
  R::Float64
end

"""
    touchstone_load(filename)

Read a file in the Touchstone format. Standard compliant.

This is the 1.1 spec:
https://ibis.org/connector/touchstone_spec11.pdf
and the 2.0 spec:
https://ibis.org/touchstone_ver2.0/touchstone_ver2_0.pdf

Outputs un-normalized network parameters with frequency units of Hz. Don't
enforce any particular file extension or try to infer the number of ports from
the extension.
"""
function touchstone_load(filename)

    #open a file handle
    io = open(filename)

    out = touchstone_parse(io)

    close(io)

    return out
end

"""
    touchstone_parse(io::IO)

Parse the Touchstone file described by the IOBuffer or IOStream `io`.

This is the 1.1 spec:
https://ibis.org/connector/touchstone_spec11.pdf
and the 2.0 spec:
https://ibis.org/touchstone_ver2.0/touchstone_ver2_0.pdf

# Examples
```jldoctest
str="!Example 1:\n!1-port S-parameter file, single frequency point\n# MHz S MA R 50\n!freq magS11 angS11\n2.000 0.894 -12.136";
println(str);
JosephsonCircuits.touchstone_parse(IOBuffer(str))

# output
!Example 1:
!1-port S-parameter file, single frequency point
# MHz S MA R 50
!freq magS11 angS11
2.000 0.894 -12.136
JosephsonCircuits.TouchstoneFile([2.0e6], [0.874020294860635 - 0.18794819544685323im;;;], "mhz", "s", "ma", 50.0, 1.0, 1, "12_21", 1, 0, [50.0], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 1:", "1-port S-parameter file, single frequency point", "freq magS11 angS11"], [2.0, 0.894, -12.136], Float64[])
```
```jldoctest
str="!Example 4 (Version 2.0):\n! 4-port S-parameter data\n! Default impedance is overridden by [Reference]\n! Data cannot be represented using 1.0 syntax\n! Note that the [Reference] keyword arguments appear on a separate line\n[Version] 2.0\n# GHz S MA R 50\n[Number of Ports] 4\n[Reference]\n50 75 0.01 0.01\n[Number of Frequencies] 1\n[Network Data]\n5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1\n0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2\n0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3\n0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4";
println(str);
JosephsonCircuits.touchstone_parse(IOBuffer(str))

# output
!Example 4 (Version 2.0):
! 4-port S-parameter data
! Default impedance is overridden by [Reference]
! Data cannot be represented using 1.0 syntax
! Note that the [Reference] keyword arguments appear on a separate line
[Version] 2.0
# GHz S MA R 50
[Number of Ports] 4
[Reference]
50 75 0.01 0.01
[Number of Frequencies] 1
[Network Data]
5.00000 0.60 161.24 0.40 -42.20 0.42 -66.58 0.53 -79.34 !row 1
0.40 -42.20 0.60 161.20 0.53 -79.34 0.42 -66.58 !row 2
0.42 -66.58 0.53 -79.34 0.60 161.24 0.40 -42.20 !row 3
0.53 -79.34 0.42 -66.58 0.40 -42.20 0.60 161.24 !row 4
JosephsonCircuits.TouchstoneFile([5.0e9], [-0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im; 0.29632183851470006 - 0.2686882357291961im -0.5679895560694177 + 0.1933594171383067im 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im; 0.16693665375723588 - 0.38539869438327984im 0.09803970583787712 - 0.5208533537179372im -0.5681244079815996 + 0.1929628385351877im 0.29632183851470006 - 0.2686882357291961im; 0.09803970583787712 - 0.5208533537179372im 0.16693665375723588 - 0.38539869438327984im 0.29632183851470006 - 0.2686882357291961im -0.5681244079815996 + 0.1929628385351877im;;;], "ghz", "s", "ma", 50.0, 2.0, 4, "12_21", 1, 0, [50.0, 75.0, 0.01, 0.01], String[], "Full", Tuple{Char, Vector{Int64}}[], ["Example 4 (Version 2.0):", " 4-port S-parameter data", " Default impedance is overridden by [Reference]", " Data cannot be represented using 1.0 syntax", " Note that the [Reference] keyword arguments appear on a separate line", "row 1", "row 2", "row 3", "row 4"], [5.0, 0.6, 161.24, 0.4, -42.2, 0.42, -66.58, 0.53, -79.34, 0.4  …  0.4, -42.2, 0.53, -79.34, 0.42, -66.58, 0.4, -42.2, 0.6, 161.24], Float64[])
```
"""
function touchstone_parse(io::IO)

    # header
    readoptionline = true
    optionline = TouchstoneOptionLine("Hz","S","MA",50.0)
    comments = String[]
    version = 0.0
    numberofports = 0
    twoportdataorder = ""
    numberoffrequencies = 0
    numberofnoisefrequencies = 0
    reference = Float64[]
    information = String[]
    matrixformat = "Full"
    mixedmodeorder = Tuple{Char, Vector{Int}}[]

    # network and noise data
    networkdata = Float64[]
    noisedata = Float64[]
    normalization = 1.0

    if eof(io)
        error("Input is empty and thus not a valid touchstone file.")
    end

    #loop over the contents of the header
    while !eof(io)

        # discard anything after the comment character !
        line = stripcommentslowercase!(comments,readline(io))

        # parse the line
        if isempty(line)
            # don't bother parsing if the string is empty
            nothing
        
        elseif isversion(line)
            # parse the [Version] line
            if !iszero(version)
                error("Only one [Version] keyword allowed: $(line)")
            else
                version=parseversion(line)
            end

        elseif isoptionline(line) && readoptionline
            optionline = parseoptionline(line)

            # set this to false so that any subsequent option lines are ignored
            # during header parsing
            readoptionline = false

            # if this is a version 1.0 or 1.1 file then the header is done
            # and we break this loop and collect the network data.
            if iszero(version)
                version = 1.0
                break
            end

        elseif isoptionline(line) && !readoptionline
            # the spec says a second options line should be ignored, but the 
            # golden parser tsck2 throws an error. 
            error("Invalid secondary options line: $(line)")
        
        elseif isnumberofports(line)
            # parse the [Number of Ports] line
            if !iszero(numberofports)
                error("Only one [Number of Ports] keyword allowed: $(line)")
            else
                numberofports = parsenumberofports(line)
            end
        
        elseif istwoportdataorder(line)
            # parse the [Two-Port Data Order] line. required if number of ports is 2
            if numberofports == 2
                if !isempty(twoportdataorder)
                    error("Only one [Two-Port Data Order] line allowed: $(line)")
                else
                    twoportdataorder = parsetwoportdataorder(line)
                end
            else
                error("[Two-Port Data Order] is only allowed if [Number of Ports] is 2: $(line)")
            end

        elseif isbegininformation(line)
            if iszero(numberofports)
                error("[Number of Ports] must be before [Begin Information]")
            end
            parseinformation!(information, comments, io)

        elseif isendinformation(line)
            error("The [End Information] line should have been parsed by the parseinformation() function.")
        
        elseif isnumberoffrequencies(line)
            numberoffrequencies = parsenumberoffrequencies(line)
        
        elseif isnumberofnoisefrequencies(line)
            numberofnoisefrequencies = parsenumberofnoisefrequencies(line)
        
        elseif isreference(line)
            if !isempty(reference)
                error("Only one [Reference] keyword allowed: $(line)")
            end
            if iszero(numberofports)
                error("[Number of Ports] must be before [Reference]")
            end
            parsereference!(reference, comments, line, numberofports, io)
        
        elseif ismatrixformat(line)
            matrixformat = parsematrixformat(line)
        elseif ismixedmodeorder(line)
            parsemixedmodeorder!(mixedmodeorder,line)
        elseif line == "[network data]"
            break
        else
            error("Unknown line type: $(line)")
        end
    end


    # parse the network data
    nfreq,nvals =  parsenetworkdata!(networkdata, comments, io)

    if nvals == 0
        error("No network data.")
    end

    # check if the number of ports and frequencies inferred from reading the
    # data are consistent with information from the header, if available.
    if version < 2.0
        # version 1 files have no information on the number of ports in the
        # header and the information in the filename is unreliable so take the
        # number of ports form the network data.
        #  v1 files only have full matrices
        @assert numberofports == 0
        numberofports = convert(Int,sqrt((nvals -1)/2))
    else
        # compare the number of ports in the header with the number of ports
        # inferred from reading the network data. 
        if matrixformat == "Full"
            if numberofports != convert(Int,sqrt((nvals -1)/2))
                error("Number of ports not consistent between header and network data.")
            end
        else
            if numberofports != convert(Int,(sqrt(4*nvals-3)-1)/2)
                error("Number of ports not consistent between header and network data.")
            end
        end
    end


    # check if the number of frequencies inferred from reading the
    # data are consistent with information from the header, if available.
    if version < 2.0
        # version 1 files have no information on the number of frequencies
        # in the header so take that from the noise data.
        @assert numberoffrequencies == 0
        numberoffrequencies = nfreq
    else
        # compare the number of frequencies specified in the header
        # with the number of frequencies in the data file.
        @assert numberoffrequencies == nfreq
    end

    # parse the noise data
    nnoisefreq = parsenoisedata!(noisedata, comments, io)

    # check if the number of noise frequencies inferred from reading the
    # data are consistent with information from the header, if available.
    if version < 2.0
        # version 1 files have no information on the number of noise frequencies
        # in the header so take that from the noise data.
        @assert numberofnoisefrequencies == 0
        numberofnoisefrequencies = nnoisefreq
    else
        # compare the number of noise frequencies specified in the header
        # with the number of noise frequencies in the data file.
        @assert numberofnoisefrequencies == nnoisefreq
    end

    # set the default matrix format if one is not already set. 
    if isempty(matrixformat)
        matrixformat = "Full"
    end

    # calculate the cartesian indices of each element of the network data
    # corresponding to the position in the scattering matrix. 
    if version < 2
        if numberofports == 2
        # a v1 file with two ports always has the 21_12 port order
            twoportdataorder = "21_12"
        else
            # a v1 file with anything other than two ports has the 12_21 order
            twoportdataorder = "12_21"
        end
    else
        if numberofports == 2 && isempty(twoportdataorder)
            error("two-port data order must be defined for a v2 file with two ports.")
        elseif numberofports != 2 && !isempty(twoportdataorder)
            error("two-port data order keyword is not allowed for a v2 file with more than two ports.")
        elseif isempty(twoportdataorder)
            twoportdataorder = "12_21"
        end
    end

    frequencies, N = networkdatatoarray(networkdata, numberofports,
        numberoffrequencies, matrixformat, twoportdataorder,
        optionline.parameter, optionline.frequencyunit, optionline.format,
        optionline.R, version)

    if isempty(reference)
        for i in 1:numberofports
            push!(reference, optionline.R)
        end
    end

    return TouchstoneFile(frequencies, N, optionline.frequencyunit,
        optionline.parameter, optionline.format, optionline.R, version,
        numberofports, twoportdataorder, numberoffrequencies,
        numberofnoisefrequencies, reference, information, matrixformat,
        mixedmodeorder, comments, networkdata, noisedata)
end



"""
    touchstone_save(filename::String,frequencies::AbstractArray,
        N::AbstractArray;version=1.0,reference=[50.0,50.0],R = 50.0,format="RI",
        comments=[""],twoportdataorder="12_21",matrixformat="Full",frequencyunit="Hz")

Write a file in the Touchstone format. Standards compliant except does not
support writing noise data.

This is the 1.1 spec:
https://ibis.org/connector/touchstone_spec11.pdf
and the 2.0 spec:
https://ibis.org/touchstone_ver2.0/touchstone_ver2_0.pdf
"""
function touchstone_save(filename::String,frequencies::AbstractVector,
    N::AbstractArray;version=1.0,reference=[50.0,50.0],R = 50.0,format="RI",
    parameter = "S", comments=[""],twoportdataorder="12_21",matrixformat="Full",
    frequencyunit="Hz")

    # check that the network data doesn't have too many dimensions
    if size(N)[1] == size(N)[2]
        numberofports = size(N)[1]
    else
        error("Network data arrays are not square.")
    end

    # check the filename. it the file doesn't have a .ts or .sNp extension,
    # then add the .sNp extension, where N is the number of ports. 
    if (m = match(r"\.ts$|\.s\dp$",lowercase(filename))) != nothing
        if version == 1.0 || version == 1.1
            if m.match != ".s$(numberofports)p"
                println("Warning: Extension of $(m.match) is not the recommended extension of .s$(numberofports)p for a version 1.0 file with $(numberofports) ports.")
            end
        else
            if m.match != ".s$(numberofports)p" && m.match != ".ts"
                println("Warning: Extension of $(m.match) is not the recommended extension of .ts or .s$(numberofports)p for a file with $(numberofports) ports.")
            end
        end
    else
        println("Warning: Adding extension of .s$(numberofports)p")
        filename = filename*".s$(numberofports)p"
    end

    # open up a file and create an IOStream io
    io = open(filename,"w")

    # write the touchstone file to io.
    touchstone_write(io, frequencies, N; version = version, reference = reference,
        R = R, format = format, parameter = parameter, comments = comments,
        twoportdataorder = twoportdataorder, matrixformat = matrixformat,
        frequencyunit = frequencyunit)

    close(io)

    return nothing
end


"""
    touchstone_write(io::IO,frequencies::AbstractVector,N::AbstractArray;
        version=1.0,reference=[50.0,50.0],R = 50.0,format="RI",parameter = "S",
        comments=[""],twoportdataorder="12_21",matrixformat="Full",frequencyunit="Hz")

Write a Touchstone file to the IOStream or IOBuffer `io`.

# Examples
```jldoctest
julia> io = IOBuffer();JosephsonCircuits.touchstone_write(io, [1.0e9, 2.0e9, 10.0e9], [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im];version=1.0,R=50.0,format="RI",frequencyunit="Hz",comments=["Example 4:","2-port S-parameter file, three frequency points"]);println(String(take!(io)))
!Example 4:
!2-port S-parameter file, three frequency points
# Hz S RI R 50.0
! freq ReS11 ImS11 ReS21 ImS21 ReS12 ImS12 ReS22 ImS22 
1.0e9 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211
2.0e9 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054
1.0e10 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336


julia> io = IOBuffer();JosephsonCircuits.touchstone_write(io, [1.0e9, 2.0e9, 10.0e9], [0.3926 - 0.1211im -0.0003 - 0.0021im; -0.0003 - 0.0021im 0.3926 - 0.1211im;;; 0.3517 - 0.3054im -0.0096 - 0.0298im; -0.0096 - 0.0298im 0.3517 - 0.3054im;;; 0.3419 + 0.3336im -0.0134 + 0.0379im; -0.0134 + 0.0379im 0.3419 + 0.3336im];version=2.0,R=50.0,format="RI",frequencyunit="Hz",comments=["Example 4:","2-port S-parameter file, three frequency points"]);println(String(take!(io)))
!Example 4:
!2-port S-parameter file, three frequency points
[Version] 2.0
# Hz S RI R 50.0
[Number of Ports] 2
[Two-Port Data Order] 12_21
[Number of Frequencies] 3
[Reference] 50.0 50.0
[Network Data]
! freq ReS11 ImS11 ReS12 ImS12 ReS21 ImS21 ReS22 ImS22 
1.0e9 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211
2.0e9 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054
1.0e10 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336
[End]

```
"""
function touchstone_write(io::IO,
    f::Vector{Float64}, N::Array{Complex{Float64},3};
    frequencyunit::String = "GHz",
    parameter::String = "S",
    format::String = "MA",
    R::Float64 = 50.0, 
    version::Float64 = 2.0,
    twoportdataorder::String = "",
    reference::Vector{Float64} = Float64[],
    information::Vector{String} = String[],
    matrixformat::String = "Full",
    mixedmodeorder::Vector{Tuple{Char, Vector{Int}}} = Tuple{Char, Vector{Int}}[],
    comments::Vector{String} = String[],
    noisedata::Vector{Float64} = Float64[],
    )
return touchstone_write(io,touchstone_file(f, N;frequencyunit, parameter, format, R, version,
    twoportdataorder, reference, information, matrixformat, mixedmodeorder,
    comments, noisedata))
end

"""
    touchstone_write(io::IO, ts::TouchstoneFile)

Write a Touchstone file specified by the TouchstoneFile object `ts` to the 
IOStream or IOBuffer `io`.
"""
function touchstone_write(io::IO,ts::TouchstoneFile)

    # write a verion 1.0 file
    if ts.version < 2.0
        # write the comments
        for comment in ts.comments
            write(io,"!$(comment)\n")
        end

        # write the option line
        write(io,"# $(ts.frequencyunit) $(ts.parameter) $(ts.format) R $(ts.R)\n")

        indices = matrixindices(ts.numberofports, ts.matrixformat, ts.twoportdataorder)

        # write a comment with the ports
        write(io,"! freq ")
        fmt = lowercase(ts.format)
        for j = 1:length(indices)
            if fmt == "ri"
                write(io,"Re$(ts.parameter)$(indices[j][1])$(indices[j][2]) Im$(ts.parameter)$(indices[j][1])$(indices[j][2]) ")
            elseif fmt == "ma"
                write(io,"mag$(ts.parameter)$(indices[j][1])$(indices[j][2]) ang$(ts.parameter)$(indices[j][1])$(indices[j][2]) ")
            elseif fmt == "db"
                write(io,"logmag$(ts.parameter)$(indices[j][1])$(indices[j][2]) ang$(ts.parameter)$(indices[j][1])$(indices[j][2]) ")
            else
                error("Error: Unknown format")
            end
        end
        write(io,"\n")


        if !isempty(ts.mixedmodeorder)
            error("Version 1.0 or 1.1 files do not support mixed-mode order (common and differential modes. Save as a Version 2.0 file or delete mixe-mode order data.")
        end

        # normalize Z, Y, G, or H by the reference impedance per the spec.
        p = lowercase(ts.parameter)
        if p == "z" || p == "y" || p == "g" || p == "h"
            normalization = ts.R
        elseif p == "s"
            normalization = 1
        else
            error("Unknown parameter")
        end

        # write the network data
        for i = 1:length(ts.f)
            k=1
            # write the frequencies
            write(io,"$(ts.networkdata[(i-1)*(length(indices)*2+1)+1])")

            for j = 1:length(indices)
                # write the network parameters
                write(io," $(ts.networkdata[(i-1)*(length(indices)*2+1)+2*j])")
                write(io," $(ts.networkdata[(i-1)*(length(indices)*2+1)+2*j+1])")

                # only 4 pairs allowed per line
                if  k  == 4  && j != length(indices)
                    k = 1
                    write(io,"\n")
                elseif ts.numberofports == 2
                    nothing
                elseif mod(j,ts.numberofports) == 0  && j != length(indices)
                    k = 0
                    write(io,"\n")
                end
                k+=1
            end
            write(io,"\n")
        end

        # write the noise data
        if !isempty(ts.noisedata)
            write(io,"[Noise Data]\n")
            for i in 1:length(ts.noisedata)
                write(io,"$(ts.noisedata[i])")
                if mod(i,5) == 0
                    write(io,"\n")
                else
                    write(io," ")
                end
            end
        end

    # write a version 2.0 file
    else
        # write the comments
        for comment in ts.comments
            write(io,"!$(comment)\n")
        end

        # write the version keyword
        write(io,"[Version] $(ts.version)\n")

        # write the option line
        write(io,"# $(ts.frequencyunit) $(ts.parameter) $(ts.format) R $(ts.R)\n")

        indices = matrixindices(ts.numberofports, ts.matrixformat, ts.twoportdataorder)

        # write the number of ports
        write(io,"[Number of Ports] $(ts.numberofports)\n")

        # write the information section
        if !isempty(ts.information)
            write(io,"[Begin Information]")
            for inform in ts.information
                write(io," $(inform)\n")
            end
            write(io,"[End Information]")

        end

        # write the two-port data order keyword if there are two ports
        if ts.numberofports == 2
            write(io,"[Two-Port Data Order] $(ts.twoportdataorder)\n")
        end

        # write the number of frequencies
        write(io,"[Number of Frequencies] $(ts.numberoffrequencies)\n")

        if !iszero(ts.numberofnoisefrequencies)
            # write the number of noise frequencies
            write(io,"[Number of Noise Frequencies] $(ts.numberofnoisefrequencies)\n")
        end

        # write the reference keyword
        write(io,"[Reference]")
        for r in ts.reference
            write(io," $(r)")
        end
        write(io,"\n")

        # write the matrix format
        if ts.matrixformat != "Full"
            write(io,"[Matrix Format] $(ts.matrixformat)\n")
        end

        # write the mixed-mode order information
        if !isempty(ts.mixedmodeorder)
            if length(ts.mixedmodeorder) != ts.numberofports
                error("Number of mixed mode order descriptors must match the number of ports.")
            end
        end


        # write the network data keyword
        write(io,"[Network Data]\n")

        # write a comment with the ports
        write(io,"! freq ")
        fmt = lowercase(ts.format)
        for j = 1:length(indices)
            if fmt == "ri"
                write(io,"Re$(ts.parameter)$(indices[j][1])$(indices[j][2]) Im$(ts.parameter)$(indices[j][1])$(indices[j][2]) ")
            elseif fmt == "ma"
                write(io,"mag$(ts.parameter)$(indices[j][1])$(indices[j][2]) ang$(ts.parameter)$(indices[j][1])$(indices[j][2]) ")
            elseif fmt == "db"
                write(io,"logmag$(ts.parameter)$(indices[j][1])$(indices[j][2]) ang$(ts.parameter)$(indices[j][1])$(indices[j][2]) ")
            else
                error("Error: Unknown format")
            end
        end
        write(io,"\n")

        # write the network data
        for i = 1:length(ts.f)
            k=1
            # write the frequencies
            write(io,"$(ts.networkdata[(i-1)*(length(indices)*2+1)+1])")

            for j = 1:length(indices)
                # write the network parameters
                write(io," $(ts.networkdata[(i-1)*(length(indices)*2+1)+2*j])")
                write(io," $(ts.networkdata[(i-1)*(length(indices)*2+1)+2*j+1])")

                # only 4 pairs allowed per line
                if  k  == 4  && j != length(indices)
                    k = 1
                    write(io,"\n")
                elseif ts.numberofports == 2
                    nothing
                elseif mod(j,ts.numberofports) == 0  && j != length(indices)
                    k = 0
                    write(io,"\n")
                end
                k+=1
            end
            write(io,"\n")
        end

        # write the noise data
        if !isempty(ts.noisedata)
            write(io,"[Noise Data]\n")
            for i in 1:length(ts.noisedata)
                write(io,"$(ts.noisedata[i])")
                if mod(i,5) == 0
                    write(io,"\n")
                else
                    write(io," ")
                end
            end
        end

        write(io,"[End]")

    end

    return nothing
end


"""
    touchstone_file(f::Vector{Float64}, N::Array{Complex{Float64}};
        frequencyunit::String = "GHz",
        parameter::String = "S",
        format::String = "MA",
        R::Float64 = 50.0, 
        version::Float64 = 2.0,
        twoportdataorder::String = "",
        numberofnoisefrequencies::Int = 0,
        reference::Vector{Float64} = Float64[],
        information::Vector{String} = String[],
        matrixformat::String = "Full",
        mixedmodeorder::Vector{Tuple{Char, Vector{Int}}} = Tuple{Char, Vector{Int}}[],
        comments::Vector{String} = String[],
        noisedata::Vector{Float64} = Float64[])

Generate a TouchstoneFile struct from the frequency vector `f` in units of Hz
and the complex network data array `N` where the frequency axis is the last
dimension. All other arguments are optional.
"""
function touchstone_file(f::Vector{Float64}, N::Array{Complex{Float64},3};
    frequencyunit::String = "GHz",
    parameter::String = "S",
    format::String = "MA",
    R::Float64 = 50.0, 
    version::Float64 = 2.0,
    # numberofports::Int = 0,
    twoportdataorder::String = "",
    # numberoffrequencies::Int,
    # numberofnoisefrequencies::Int = 0,
    reference::Vector{Float64} = Float64[],
    information::Vector{String} = String[],
    matrixformat::String = "Full",
    mixedmodeorder::Vector{Tuple{Char, Vector{Int}}} = Tuple{Char, Vector{Int}}[],
    comments::Vector{String} = String[],
    # networkdata::Vector{Float64} = Float64[],
    noisedata::Vector{Float64} = Float64[]
)

# determine the number of frequencies from `f`
numberoffrequencies = length(f)

# determine the number of ports and number of frequencies from the network 
# data `N`. check for consistency.
if numberoffrequencies != size(N,3)
    error("The size of the last axis of network data N must equal the number of frequecies")
end

if size(N,1) != size(N,2)
    error("Network data matrix must be square")
end
numberofports = size(N,1)

# check the version
# should be 1.0, 1.1, or 2.0.
if !(version == 1.0 || version == 1.1 || version == 2.0)
    error("Version must be 1.0, 1.1, or 2.0")
end

# check the frequencyunit
# should be Hz, kHz, MHz, or GHz, not case sensitive.
fu = lowercase(frequencyunit)
if !(fu == "hz" || fu == "khz" || fu == "mhz" || fu == "ghz")
    error("Unknown frequency unit")
end

# check the format
# should be MA, RI, or DB, not case sensitive.
fmt = lowercase(format)
if !(fmt == "ma" || fmt == "ri" || fmt == "db")
    error("Unknown format")
end

# check the parameter
# should be MA, RI, or DB, not case sensitive.
p = lowercase(parameter)
if !(p == "s" || p == "y" || p == "z" || p == "h" || p == "g")
    error("Unknown parameter")
end

# check the impedance `R` and per port impedance `reference`.
# if the per port impedance `reference` is given, it overrides `R`. if 
# `reference is not given, then populate `reference` from `R` and the number
# of ports
if isempty(reference)
    reference = ones(Float64,numberofports)*R
end

if length(reference) != numberofports
    error("Number of per port impedances must equal number of ports")
end

if version < 2.0 && !allequal(reference)
    error("The port impedances are not equal, so we cannot generate a Touchstone file with version < 2.0.")
end

# check the matrixformat
if isempty(matrixformat)
    matrixformat = "Full"
end

if version < 2.0
    if lowercase(matrixformat) != "full"
        error("For Touchstone files with version less than 2.0 only the Full matrix format is required.") 
    end
else
    if !(lowercase(matrixformat) == "full" || lowercase(matrixformat) == "lower" || lowercase(matrixformat) == "upper")
        error("For Touchstone files with version 2.0 the valid formats are Full, Lower, or Upper.") 
    end
end

# check the twoportdataorder
if version < 2
    if numberofports == 2
    # a v1 file with two ports always has the 21_12 port order
        twoportdataorder = "21_12"
    else
        # a v1 file with anything other than two ports has the 12_21 order
        twoportdataorder = "12_21"
    end
else
    if isempty(twoportdataorder)
        twoportdataorder = "12_21"
        # error("two-port data order must be defined for a v2 file with two ports.")
    elseif twoportdataorder == "12_21" || twoportdataorder == "21_12"
        nothing
    else
        error("Unknown two-port data order.")
    end
end

# check the mixedmodeorder

# check the noisedata
numberofnoisefrequencies = 0
if rem(length(noisedata),5) == 0
    numberofnoisefrequencies = length(noisedata) ÷ 5
else
    error("Invalid number of elements in noise data.")
end

# generate the networkdata
networkdata = arraytonetworkdata(f, N, numberofports, numberoffrequencies, 
    matrixformat, twoportdataorder, parameter, frequencyunit, format, R, version)

# return the TouchstoneFile structure
return TouchstoneFile(f, N, frequencyunit,
    parameter, format, R, version,
    numberofports, twoportdataorder, numberoffrequencies,
    numberofnoisefrequencies, reference, information, matrixformat,
    mixedmodeorder, comments, networkdata, noisedata)
end


function arraytonetworkdata(frequencies,N, numberofports, numberoffrequencies, 
    matrixformat, twoportdataorder, parameter, frequencyunit, format, R, version)

    nvals = 0
    # the number of values per frequency
    mf = lowercase(matrixformat)
    if mf == "full"
        nvals = (2*numberofports^2+1)
    elseif mf == "lower" || mf == "upper"
        nvals = ((2*numberofports+1)^2+3)÷4
    else
        error("Unknown matrixformat.")
    end

    indices = matrixindices(numberofports,matrixformat,twoportdataorder)

    # scale the parameters if it's a v1 file and not a scattering parameter. 
    p = lowercase(parameter)
    if p == "s"
        normalization = 1.0
    elseif (p == "z" || p == "y" || p  == "g" || p  == "h") && version < 2.0
        normalization = R
    elseif (p == "z" || p == "y" || p == "g" || p  == "h") && version == 2.0
        normalization = 1.0
    else
        error("Unknown format or version")
    end

    networkdata = Vector{Float64}(undef,nvals*numberoffrequencies)

    # copy over the frequency data
    networkdata[1:(2*length(indices)+1):end] .= frequencies ./ frequencyscale(frequencyunit)

    # and  the network data
    fmt = lowercase(format)
    for n = 1:length(indices)
        if fmt == "db"
            networkdata[2*n:(2*length(indices)+1):end] .= 10*log10.(abs2.(N[indices[n],:]/normalization))
            networkdata[2*n+1:(2*length(indices)+1):end] .= 180/pi*angle.(N[indices[n],:]/normalization)
        elseif fmt == "ma"
            networkdata[2*n:(2*length(indices)+1):end] .= abs.(N[indices[n],:]/normalization)
            networkdata[2*n+1:(2*length(indices)+1):end] .= 180/pi*angle.(N[indices[n],:]/normalization)
        elseif fmt == "ri"
            networkdata[2*n:(2*length(indices)+1):end] .= real.(N[indices[n],:]/normalization)
            networkdata[2*n+1:(2*length(indices)+1):end] .= imag.(N[indices[n],:]/normalization)
        end
    end

  return networkdata
end


function networkdatatoarray(networkdata, numberofports, numberoffrequencies, 
    matrixformat, twoportdataorder, parameter, frequencyunit, format, R, version)

    # #make an empty array for the network data
    N = Array{Complex{Float64}}(undef,numberofports,numberofports,numberoffrequencies)
    frequencies = Array{Float64}(undef,numberoffrequencies)

    indices = matrixindices(numberofports,matrixformat,twoportdataorder)

    # scale the parameters if it's a v1 file and not a scattering parameter. 
    p = lowercase(parameter)
    if p == "s"
        normalization = 1.0
    elseif (p == "z" || p  == "y" || p  == "g" || p  == "h") && version < 2.0
        normalization = R
    elseif (p == "z" || p == "y" || p == "g" || p == "h") && version == 2.0
        normalization = 1.0
    else
        error("Error: Unknown parameter.")
    end

    # copy over the frequency data
    frequencies .= networkdata[1:(2*length(indices)+1):end] .* frequencyscale(frequencyunit)

    # and  the network data
    fmt = lowercase(format)
    for n = 1:length(indices)
        if fmt == "db"
            N[indices[n],:] .= 10 .^((1/20).*networkdata[2*n:(2*length(indices)+1):end].*normalization)
            N[indices[n],:]  .*= exp.((pi/180*im).*networkdata[2*n+1:(2*length(indices)+1):end])
        elseif fmt == "ma"
            N[indices[n],:] .= networkdata[2*n:(2*length(indices)+1):end].*normalization
            N[indices[n],:] .*= exp.((pi/180*im).*networkdata[2*n+1:(2*length(indices)+1):end])
        elseif fmt == "ri"
            N[indices[n],:] .= networkdata[2*n:(2*length(indices)+1):end].*normalization
            N[indices[n],:] .+= im.*networkdata[2*n+1:(2*length(indices)+1):end].*normalization
        end

        # if the format is upper or lower, the matrix is assumed to be symmetric
        # so fill in the other parts
        if lowercase(matrixformat) != "full"
            if indices[n][1] != indices[n][2]
                N[indices[n][2],indices[n][1],:] .= N[indices[n][1],indices[n][2],:]
            end
        end
    end

  return frequencies, N
end


"""
    matrixindices(nports,format)

Return the cartesian indices of the elements of a scattering matrix given the
number of ports `nports` and the format `format` which can be "Full", "Upper",
or "Lower". 

# Examples
```jldoctest
julia> JosephsonCircuits.matrixindices(2,"Full",printflag=true)
11 12 
21 22 
4-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 1)
 CartesianIndex(1, 2)
 CartesianIndex(2, 1)
 CartesianIndex(2, 2)

julia> JosephsonCircuits.matrixindices(2,"Upper",printflag=true)
11 12 
   22 
      3-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 1)
 CartesianIndex(1, 2)
 CartesianIndex(2, 2)

julia> JosephsonCircuits.matrixindices(2,"Lower",printflag=true)
11 
21 22 
3-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 1)
 CartesianIndex(2, 1)
 CartesianIndex(2, 2)
```
"""
function matrixindices(nports,format; printflag = false)
    return matrixindices(nports,format,"12_21", printflag = printflag)
end

"""
    matrixindices(nports,format,twoportdataorder)

Return the cartesian indices of the elements of a scattering matrix given the
number of ports `nports` and the format `format` which can be "Full", "Upper",
or "Lower". The two port data order `twoportdataorder` can be "`12_21`" or "`21_12`"
for 2 ports but must be "`12_21`" for other numbers of ports.

# Examples
```jldoctest
julia> JosephsonCircuits.matrixindices(2,"Full","12_21")
4-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 1)
 CartesianIndex(1, 2)
 CartesianIndex(2, 1)
 CartesianIndex(2, 2)

julia> JosephsonCircuits.matrixindices(2,"Full","21_12")
4-element Vector{CartesianIndex{2}}:
 CartesianIndex(1, 1)
 CartesianIndex(2, 1)
 CartesianIndex(1, 2)
 CartesianIndex(2, 2)
```
"""
function matrixindices(nports,format,twoportdataorder;printflag = false)
    indices = CartesianIndex{2}[]
    
    fmt = lowercase(format)
    if fmt == "lower"
        ncol = 1
        ntotal = nports*(1+nports) ÷ 2
    elseif fmt == "upper"
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
            if fmt == "lower"
                ncol+= 1
            elseif fmt == "upper"
                ncol-= 1
            end
        end

        #j = mod(n-1,ncol)+1 # this is the column index
        # this is the column index
        j = n - nold + 1

        if fmt == "upper"
            if printflag
                print("$(i)$(j+nports-ncol) ")
            end
            push!(indices,CartesianIndex(i,j+nports-ncol))
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
            if fmt == "upper"
                for k = 1:1:nports-ncol+1
                    if printflag
                        print("   ")
                    end
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
    isoptionline(line::String)

Return `true` if the string `line` is the option line of a Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.isoptionline("# MHz Z MA R 75")
true

julia> JosephsonCircuits.isoptionline("[number of ports] 1")
false
```
"""
function isoptionline(line::String)
    return findfirst('#',line) == 1
end

"""
    parseoptionline(line::String)

Return a struct TouchstoneOptionLine which contains the option line of a
Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.parseoptionline("# MHz Z MA R 75")
JosephsonCircuits.TouchstoneOptionLine("MHz", "Z", "MA", 75.0)

julia> JosephsonCircuits.parseoptionline("# MHz H RI R 75")
JosephsonCircuits.TouchstoneOptionLine("MHz", "H", "RI", 75.0)
```
"""
function parseoptionline(line::String)
    # parse the elements of the option line. if nothing is found use the
    # default. an empty option line is perfectly acceptable. 

    indices = findfirst(r"(?i)(hz|khz|mhz|ghz)",line)
    if indices == nothing
        frequencyunit = "ghz"
    else
        frequencyunit = line[indices]
    end

    indices = findfirst(r"(?i)( s | y | z | g | h )",line)
    if indices == nothing
        parameter = "s"
    else
        parameter = strip(line[indices])
    end

    indices = findfirst(r"(?i)(db|ma|ri)",line)
    if indices == nothing
        format = "ma"
    else
        format = line[indices]
    end

    indices = findfirst(r"(?i)( r )",line)
    if indices == nothing
        R = 50.0
    else
        R = parse.(Float64,line[last(findfirst(r"(?i)( r )",line))+1:end])
    end

    return TouchstoneOptionLine(frequencyunit, parameter, format, R)
end

function frequencyscale(frequencyunit::String)
    fu = lowercase(frequencyunit)
    if fu == "hz"
        return 1.0
    elseif fu == "khz"
        return 1.0e3
    elseif fu == "mhz"
        return 1.0e6
    elseif fu == "ghz"
        return 1.0e9
    else
        error("Unknown frequency unit $frequencyunit")
    end
end


"""
    isversion(line::String)

Return `true` if the string `line` is the [version] line of a Touchstone file. 

# Examples
```jldoctest
julia> JosephsonCircuits.isversion("[version] 1.0")
true

julia> JosephsonCircuits.isversion("[number of ports] 1")
false
```
"""
function isversion(line::String)
    indices = findfirst("[version]",line)
    return indices != nothing && first(indices) == 1
end

"""
    parseversion(line::String)

Return the version parsed from the [version] line of a Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.parseversion("[version] 1.0")
1.0
```
"""
function parseversion(line::String)
    return parse(Float64,line[10:end])
end

"""
    isnumberofports(line::String)

Return `true` if the string `line` is the [number of ports] line of a Touchstone
file.

# Examples
```jldoctest
julia> JosephsonCircuits.isnumberofports("[number of ports] 1")
true

julia> JosephsonCircuits.isnumberofports("[version] 1.0")
false
```
"""
function isnumberofports(line::String)
    indices = findfirst("[number of ports]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parsenumberofports(line::String)

Return the number of ports parsed from the [number of ports] line of a Touchstone
file.

# Examples
```jldoctest
julia> JosephsonCircuits.parsenumberofports("[number of ports] 1")
1
```
"""
function parsenumberofports(line::String)
    return parse(Int,line[18:end])
end

"""
    istwoportdataorder(line::String)

Return `true` if the string `line` is the [two-port data order] line of a
Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.istwoportdataorder("[two-port data order] 12_21")
true

julia> JosephsonCircuits.istwoportdataorder("[version] 1.0")
false
```
"""
function istwoportdataorder(line::String)
    indices = findfirst("[two-port data order]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parsetwoportdataorder(line::String)

Return the two-port data order string parsed from the [two-port data order]
line of a Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.parsetwoportdataorder("[two-port data order] 12_21")
"12_21"

julia> JosephsonCircuits.parsetwoportdataorder("[two-port data order] 21_12")
"21_12"
```
"""
function parsetwoportdataorder(line::String)
    twoportdataorder = strip(line[22:end])
    if !(twoportdataorder == "12_21" || twoportdataorder == "21_12")
        error("Unknown [Two-Port Data Order] parameter:\n$(line)")
    end
    return String(twoportdataorder)
end

"""
    isbegininformation(line::String)

Return `true` if the string `line` is the [begin information] line of a
Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.isbegininformation("[begin information]")
true

julia> JosephsonCircuits.isbegininformation("[version] 1.0")
false
```
"""
function isbegininformation(line::String)
    indices = findfirst("[begin information]", line)
    return indices != nothing && first(indices) == 1
end

"""
    isendinformation(line::String)

Return `true` if the string `line` is the [end information] line of a
Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.isendinformation("[end information]")
true

julia> JosephsonCircuits.isendinformation("[version] 1.0")
false
```
"""
function isendinformation(line::String)
    indices = findfirst("[end information]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parseinformation!(information::Vector{String},io::IO)

Append the contents of the information section of a Touchstone file from the
IOBuffer or IOStream `io` to the vector `information`.

# Examples
```jldoctest
information = String[]
comments = String[]
io = IOBuffer("This is an information section.\n[End Information]")
JosephsonCircuits.parseinformation!(information,comments,io)
println(information)

# output
["this is an information section."]
```
"""
function parseinformation!(information::Vector{String},
        comments::Vector{String}, io::IO)
    # loop over and store the information
    while !eof(io)
        line = stripcommentslowercase!(comments,readline(io))

        # end the loop once we reach the end information line
        if isendinformation(line)
            break
        end

        # store the contents of the line in the information array
        if !isempty(line)
            push!(information,line)
        end
    end
    return nothing
end

"""
    isnumberoffrequencies(line::String)

Return `true` if the string `line` is the [number of frequencies] line of a
Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.isnumberoffrequencies("[number of frequencies]")
true

julia> JosephsonCircuits.isnumberoffrequencies("[version] 1.0")
false
```
"""
function isnumberoffrequencies(line::String)
    indices = findfirst("[number of frequencies]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parsenumberoffrequencies(line::String)

Return the number of frequencies parsed from the [number of frequencies] line
of a Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.parsenumberoffrequencies("[number of frequencies] 10")
10
```
"""
function parsenumberoffrequencies(line::String)
    numberoffrequencies = parse(Int,line[24:end])
    if numberoffrequencies < 1
        error("Error: Number of frequencies must be an integer greater than zero:\n$(line)")
    end
    return numberoffrequencies
end

"""
    isnumberofnoisefrequencies(line::String)

Return `true` if the string `line` is the [number of noise frequencies] line
of a Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.isnumberofnoisefrequencies("[number of noise frequencies]")
true

julia> JosephsonCircuits.isnumberofnoisefrequencies("[version] 1.0")
false
```
"""
function isnumberofnoisefrequencies(line::String)
    indices = findfirst("[number of noise frequencies]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parsenumberofnoisefrequencies(line::String)

Return the number of noise frequencies parsed from the [number of noise
frequencies] line of a Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.parsenumberofnoisefrequencies("[number of noise frequencies] 10")
10
```
"""
function parsenumberofnoisefrequencies(line::String)
    numberofnoisefrequencies = parse(Int,line[30:end])
    if numberofnoisefrequencies < 1
        error("Error: Number of noise frequencies must be an integer greater than zero:\n$(line)")
    end
    return numberofnoisefrequencies
end

"""
    isreference(line::String)

Return `true` if the string `line` is the [reference] line of a Touchstone
file.

# Examples
```jldoctest
julia> JosephsonCircuits.isreference("[reference]")
true

julia> JosephsonCircuits.isreference("[version] 1.0")
false
```
"""
function isreference(line::String)
    indices = findfirst("[reference]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parsereference!(reference::Vector{Float64}, comments::Vector{String},
        line::String, numberofports::Int, io::IO)

Append the contents of the [reference] section of a Touchstone file from the
IOBuffer or IOStream `io` to the vector `reference`. The reference impedance
values can be spread across multiple lines. 

# Examples
```jldoctest
io = IOBuffer("[Reference] 50.0 60.0 75.0")
numberofports = 3
comments = String[]
reference = Float64[]
line = JosephsonCircuits.stripcommentslowercase!(comments,readline(io))
JosephsonCircuits.parsereference!(reference, comments, line, numberofports, io)
println(reference)

# output
[50.0, 60.0, 75.0]
```
```jldoctest
io = IOBuffer("[Reference] 50.0 \n60.0 75.0")
numberofports = 3
comments = String[]
reference = Float64[]
line = JosephsonCircuits.stripcommentslowercase!(comments,readline(io))
JosephsonCircuits.parsereference!(reference, comments, line, numberofports, io)
println(reference)

# output
[50.0, 60.0, 75.0]
```
```jldoctest
io = IOBuffer("[Reference] 50.0 \n60.0 75.0\n[Number of Frequencies] 1")
numberofports = 3
comments = String[]
reference = Float64[]
line = JosephsonCircuits.stripcommentslowercase!(comments,readline(io))
JosephsonCircuits.parsereference!(reference, comments, line, numberofports, io)
println(reference)

# output
[50.0, 60.0, 75.0]
```
"""
function parsereference!(reference::Vector{Float64}, comments::Vector{String},
    line::String, numberofports::Int, io::IO)
    # reference can be multiple lines.

    # if there is a number on this line, parse the line
    index = findfirst(r"\d",line)
    if index != nothing
        append!(reference,parse.(Float64,split(strip(line[first(index):end]),r"\s+")))
    end

    # check if we have all of the reference values and take more data
    # if we do not.
    if length(reference) < numberofports
        while !eof(io)
            line = stripcommentslowercase!(comments,readline(io))

            if !isempty(line)
                append!(reference,parse.(Float64,split(strip(line),r"\s+")))
                if length(reference) == numberofports
                    break
                elseif length(reference) > numberofports
                    error("Too many values on [Reference] line: $(line)")
                end
            end
        end
    end

    return nothing
end

"""
    ismatrixformat(line::String)

Return `true` if the string `line` is the [matrix format] line of a Touchstone
file.

# Examples
```jldoctest
julia> JosephsonCircuits.ismatrixformat("[matrix format] full")
true

julia> JosephsonCircuits.ismatrixformat("[version] 1.0")
false
```
"""
function ismatrixformat(line::String)
    indices = findfirst("[matrix format]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parsematrixformat(line::String)

Return the two-port data order string parsed from the [two-port data order]
line of a Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.parsematrixformat("[matrix format] lower")
"Lower"

julia> JosephsonCircuits.parsematrixformat("[matrix format] upper")
"Upper"

julia> JosephsonCircuits.parsematrixformat("[matrix format] full")
"Full"
```
"""
function parsematrixformat(line::String)
    matrixformat = strip(line[16:end])
    if matrixformat == "full"
        return "Full"
    elseif matrixformat == "lower"
        return "Lower"
    elseif matrixformat == "upper"
        return "Upper"
    else
        error("Unknown format:\n$(matrixformat)")
    end
end

"""
    ismixedmodeorder(line::String)

Return `true` if the string `line` is the [mixed-mode order] line of a
Touchstone file.

# Examples
```jldoctest
julia> JosephsonCircuits.ismixedmodeorder("[mixed-mode order] full")
true

julia> JosephsonCircuits.ismixedmodeorder("[version] 1.0")
false
```
"""
function ismixedmodeorder(line::String)
    indices = findfirst("[mixed-mode order]", line)
    return indices != nothing && first(indices) == 1
end

"""
    parsemixedmodeorder!(mixedmodeorder::Vector{Tuple{Char, Vector{Int}}}, line::String)

Append the contents of the [mixed-mode order] line of a Touchstone file
from the to the vector `mixedmodeorder`.

# Examples
```jldoctest
julia> mixedmodeorder = Tuple{Char, Vector{Int}}[];JosephsonCircuits.parsemixedmodeorder!(mixedmodeorder,"[Mixed-Mode Order] D2,3 D6,5 C2,3 C6,5 S4 S1");mixedmodeorder
6-element Vector{Tuple{Char, Vector{Int64}}}:
 ('D', [2, 3])
 ('D', [6, 5])
 ('C', [2, 3])
 ('C', [6, 5])
 ('S', [4])
 ('S', [1])
```
"""
function parsemixedmodeorder!(mixedmodeorder::Vector{Tuple{Char, Vector{Int}}},
        line::String)
    # parse the [Mixed-Mode Order] line
    splitline = split(strip(lowercase(line[19:end])),r"\s+");

    for l in splitline
        splitline2 = parse.(Int,split(strip(lowercase(l[2:end])),r","));
        push!(mixedmodeorder,(uppercase(l[1]),splitline2))
    end
    return nothing
end

"""
    isnetworkdata(line::String)

Return `true` if the string `line` is the [network data] line of a
Touchstone file. 

# Examples
```jldoctest
julia> JosephsonCircuits.isnetworkdata("[network data]")
true

julia> JosephsonCircuits.isnetworkdata("[version] 1.0")
false
```
"""
function isnetworkdata(line::String)
    indices = findfirst("[network data]", line)
    return indices != nothing && first(indices) == 1
end

"""
    isend(line::String)

Return `true` if the string `line` is the [end] line of a
Touchstone file. 

# Examples
```jldoctest
julia> JosephsonCircuits.isend("[end]")
true

julia> JosephsonCircuits.isend("[version] 1.0")
false
```
"""
function isend(line::String)
    indices = findfirst("[end]", line)
    return indices != nothing && first(indices) == 1
end

"""
    isnoisedata(line::String)

Return `true` if the string `line` is the [noise data] line of a
Touchstone file. 

# Examples
```jldoctest
julia> JosephsonCircuits.isnoisedata("[noise data]")
true

julia> JosephsonCircuits.isnoisedata("[version] 1.0")
false
```
"""
function isnoisedata(line::String)
    indices = findfirst("[noise data]", line)
    return indices != nothing && first(indices) == 1
end

"""
    stripcommentslowercase!(comments::Vector{String},line::String)

Append the comment portion of a line `line` of a Touchstone file to the vector
`comments`. Return the line with the comments removed and made lowercase. 

# Examples
```jldoctest
julia> comments=String[];println(JosephsonCircuits.stripcommentslowercase!(comments,"! This is a comment"));println(comments)

[" This is a comment"]

julia> comments=String[];println(JosephsonCircuits.stripcommentslowercase!(comments,"This is a !comment"));println(comments)
this is a 
["comment"]
```
"""
function stripcommentslowercase!(comments::Vector{String},line::String)
    index = findfirst('!',line)
    if index != nothing
        push!(comments,line[index+1:end])
        return lowercase(line[1:index-1])
    else
        return lowercase(line)
    end
end


"""
    parsenetworkdata!(networkdata::Vector{Float64}, comments::Vector{String}, io::IO)

Append the contents of the networkdata section of a Touchstone file from the
IOBuffer or IOStream `io` to the vector `networkdata`.

# Examples
```jldoctest
networkdata = Float64[]
comments = String[]
io = IOBuffer("1.0000 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211\n2.0000 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054\n10.000 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336")
JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
println(networkdata)

# output
[1.0, 0.3926, -0.1211, -0.0003, -0.0021, -0.0003, -0.0021, 0.3926, -0.1211, 2.0, 0.3517, -0.3054, -0.0096, -0.0298, -0.0096, -0.0298, 0.3517, -0.3054, 10.0, 0.3419, 0.3336, -0.0134, 0.0379, -0.0134, 0.0379, 0.3419, 0.3336]
```
```jldoctest
networkdata = Float64[]
comments = String[]
io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n! NOISE PARAMETERS\n4 .7 .64 69 .38\n18 2.7 .46 -33 .40")
JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
println(networkdata)

# output
[2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0]
```
```jldoctest
networkdata = Float64[]
comments = String[]
io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n4 .7 .64 69 19\n18 2.7 .46 -33 20\n[End]")
JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
println(networkdata)

# output
[2.0, 0.95, -26.0, 3.57, 157.0, 0.04, 76.0, 0.66, -14.0, 22.0, 0.6, -144.0, 1.3, 40.0, 0.14, 40.0, 0.56, -85.0]
```
"""
function parsenetworkdata!(networkdata::Vector{Float64},
        comments::Vector{String}, io::IO)

    # finished reading the header, now read the network and noise data

    # parse the network data and the noise data
    ndatalines = 0
    nvals = 0
    freq = 0.0
    nfreq = 0
    nvalsold = 0

    # position of IO stream or buffer
    pos = 0

    #loop over the network data
    while !eof(io)
        # find the current position of the IO stream or buffer. 
        pos = position(io)

        # read a line
        line = stripcommentslowercase!(comments,readline(io))

        if isempty(line)
            #skip any lines that don't have anything
            nothing
        elseif isoptionline(line)
            # skip any additional option lines. these could occur in v1 files
            # after the first option line. 
            nothing 
        elseif isend(line)
            break
        elseif isnoisedata(line)
            # if we find a noise data line, then break the loop and read
            # the noise data in a separate function
            break
        else
            # parse the network data lines
            vals = parse.(Float64,split(strip(line),r"\s+"))

            # check if the current line is a new frequency or not
            # for v2 files, we know the number of ports so we can read the correct
            # number of data points even if spread across multiple lines. 
            # for v1 files, we do not know the number of ports, so we have to read
            # we find a different frequency or until we reach the end of the file
            # either will give us the number of ports. 

            if isodd(length(vals))
                # this line is a new frequency

                if vals[1] < freq
                    # if the new frequency is less than the old frequency, then 
                    # we have moved on to noise data. break the loop
                    # and move to the noise data reading loop. 
                    # return to the previous position. we will read the
                    # noise data in a separate function
                    seek(io,pos)

                    break

                else
                    ndatalines +=1

                    if ndatalines > 1
                        # this should be consistent between frequencies
                        if nvalsold == 0
                            nvalsold = nvals
                        elseif nvalsold != nvals
                            error("Number of ports are not consistent between lines:\n$(line)")
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
    return nfreq,nvals
end

"""
    parsenoisedata!(noisedata::Vector{Float64}, comments::Vector{String}, io::IO)

Append the contents of the networkdata section of a Touchstone file from the
IOBuffer or IOStream `io` to the vector `networkdata`.

# Examples
```jldoctest
networkdata = Float64[]
noisedata = Float64[]
comments = String[]
io = IOBuffer("1.0000 0.3926 -0.1211 -0.0003 -0.0021 -0.0003 -0.0021 0.3926 -0.1211\n2.0000 0.3517 -0.3054 -0.0096 -0.0298 -0.0096 -0.0298 0.3517 -0.3054\n10.000 0.3419 0.3336 -0.0134 0.0379 -0.0134 0.0379 0.3419 0.3336")
JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
JosephsonCircuits.parsenoisedata!(noisedata,comments,io)
println(noisedata)

# output
Float64[]
```
```jldoctest
networkdata = Float64[]
noisedata = Float64[]
comments = String[]
io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n! NOISE PARAMETERS\n4 .7 .64 69 .38\n18 2.7 .46 -33 .40")
JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
JosephsonCircuits.parsenoisedata!(noisedata,comments,io)
println(noisedata)

# output
[4.0, 0.7, 0.64, 69.0, 0.38, 18.0, 2.7, 0.46, -33.0, 0.4]
```
```jldoctest
networkdata = Float64[]
noisedata = Float64[]
comments = String[]
io = IOBuffer("2 .95 -26 3.57 157 .04 76 .66 -14\n22 .60 -144 1.30 40 .14 40 .56 -85\n[Noise Data]\n4 .7 .64 69 19\n18 2.7 .46 -33 20\n[End]")
JosephsonCircuits.parsenetworkdata!(networkdata,comments,io)
JosephsonCircuits.parsenoisedata!(noisedata,comments,io)
println(noisedata)

# output
[4.0, 0.7, 0.64, 69.0, 19.0, 18.0, 2.7, 0.46, -33.0, 20.0]
```
"""
function parsenoisedata!(noisedata, comments, io)

    # parse the network data and the noise data
    # ndatalines = 0
    nvals = 0
    freq = 0.0
    # nfreq = 0
    nvalsold = 0
    nnoisefreq = 0

    #loop over the noise data
    #it's a bit complicated because i might have read in the first line of noise 
    # data in a v1 file. 
    while !eof(io)
        
        line = stripcommentslowercase!(comments,readline(io))

        if isempty(line)
            #skip any lines that don't have anything
            nothing
        elseif isoptionline(line)
            # skip any additional option lines. these could occur in v1 files
            # after the first option line. 
            nothing 
        elseif isend(line)
            break
        elseif isnoisedata(line)
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
    return nnoisefreq
end


