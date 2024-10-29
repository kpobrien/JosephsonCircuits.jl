
"""
    Frequencies(Nharmonics::NTuple{N, Int}, Nw::NTuple{N,Int}, Nt::NTuple{N,Int},
        coords::Vector{CartesianIndex{N}}, modes::Vector{NTuple{N,Int})

A simple structure to hold time and frequency domain information for the
signals. See also [`calcfreqsrdft`](@ref) and [`calcfreqsdft`](@ref).

# Fields
- `Nharmonics::NTuple{N, Int}`: The number of harmonics for each frequency.
- `Nw::NTuple{N,Int}`: The dimensions of the frequency domain signal for a
    single node.
- `Nt::NTuple{N,Int}`: The dimensions of the time domain signal for a single
    node.
- `coords::Vector{CartesianIndex{N}}`: The coordinates of each mixing products.
- `modes::Vector{NTuple{N,Int}}`: The mode indices of each mixing product, eg.
     (0,0), (1,0), (2,1).

# Examples
```jldoctest
Nharmonics = (2,1)
Nw = (3, 3)
Nt = (4, 3)
coords = CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3)]
modes = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, -1), (1, -1), (2, -1)]
JosephsonCircuits.Frequencies(Nharmonics, Nw,Nt,coords,modes)

# output
JosephsonCircuits.Frequencies{2}((2, 1), (3, 3), (4, 3), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3)], [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, -1), (1, -1), (2, -1)])
```
"""
struct Frequencies{N}
    Nharmonics::NTuple{N, Int}
    Nw::NTuple{N, Int}
    Nt::NTuple{N, Int}
    coords::Vector{CartesianIndex{N}}
    modes::Vector{NTuple{N, Int}}
end

"""
    FourierIndices(conjsymdict::Dict{CartesianIndex{N},CartesianIndex{N}},
        vectomatmap::Vector{Int}, conjsourceindices::Vector{Int},
        conjtargetindices::Vector{Int}, hbmatmodes::Matrix{NTuple{N, Int}},
        hbmatindices::Matrix{Int})

A simple structure to hold time and frequency domain information for the
signals, particularly the indices for converting between the node flux vectors
and matrices. See also [`fourierindices`](@ref).
"""
struct FourierIndices{N}
    conjsymdict::Dict{CartesianIndex{N},CartesianIndex{N}}
    vectomatmap::Vector{Int}
    conjsourceindices::Vector{Int}
    conjtargetindices::Vector{Int}
    hbmatmodes::Matrix{NTuple{N, Int}}
    hbmatindices::Matrix{Int}
end

"""
    fourierindices(freq::Frequencies)

Generate the indices used in the RDFT or DFT and inverse RDFT or DFT and
converting between a node flux vector for solving system and the matrices for
the Fourier analysis. See also [`FourierIndices`](@ref), [`Frequencies`](@ref),
[`calcfreqsrdft`](@ref) and [`calcfreqsdft`](@ref).

"""
function fourierindices(freq::Frequencies)

    conjsymdict = conjsym(freq)
    freqindexmap, conjsourceindices, conjtargetindices = calcphiindices(freq,conjsymdict)
    Amatrixmodes, Amatrixindices = hbmatind(freq)

    return FourierIndices(
        conjsymdict,
        freqindexmap,
        conjsourceindices,
        conjtargetindices,
        Amatrixmodes,
        Amatrixindices,
    )
end

"""
    calcfreqsrdft(Nharmonics::NTuple{N,Int})

Calculate the dimensions of the RDFT in the frequency domain
and the time domain given a tuple of the number of harmonics. Eg. 0,w,2w,3w
would be 3 harmonics. Also calculate the possible modes and their coordinates
in the frequency domain RDFT array.

# Arguments
- `Nharmonics`: is a tuple of the number of harmonics to calculate for
    each frequency.

# Returns
- `Frequencies`: A simple structure to hold time and frequency domain
    information for the signal for a single node. See [`Frequencies`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.calcfreqsrdft((1,))
JosephsonCircuits.Frequencies{1}((1,), (2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])

julia> JosephsonCircuits.calcfreqsrdft((2,))
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)])

julia> JosephsonCircuits.calcfreqsrdft((3,))
JosephsonCircuits.Frequencies{1}((3,), (4,), (6,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,), CartesianIndex(4,)], [(0,), (1,), (2,), (3,)])

julia> JosephsonCircuits.calcfreqsrdft((3,3))
JosephsonCircuits.Frequencies{2}((3, 3), (4, 7), (6, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2), CartesianIndex(1, 3), CartesianIndex(2, 3)  …  CartesianIndex(3, 5), CartesianIndex(4, 5), CartesianIndex(1, 6), CartesianIndex(2, 6), CartesianIndex(3, 6), CartesianIndex(4, 6), CartesianIndex(1, 7), CartesianIndex(2, 7), CartesianIndex(3, 7), CartesianIndex(4, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (3, 1), (0, 2), (1, 2)  …  (2, -3), (3, -3), (0, -2), (1, -2), (2, -2), (3, -2), (0, -1), (1, -1), (2, -1), (3, -1)])
```
"""
function calcfreqsrdft(Nharmonics::NTuple{N,Int}) where N

    # double the size of all but the first n because 
    # only the first axis of a multi-dimensional rfft has
    # only positive frequencies.
    Nw=NTuple{N,Int}(ifelse(i == 1, val+1, 2*val+1) for (i,val) in enumerate(Nharmonics))
    
    # choose the number of time points based on the number of fourier
    # coefficients
    # changed to below because i wasn't using enough points when Nmodes=1.
    # the results contained only real values. 
    if Nw[1] == 2
        stepsperperiod = 2*Nw[1]-1
    else
        stepsperperiod = 2*Nw[1]-2
    end

    Nt = NTuple{N,Int}(ifelse(i == 1, stepsperperiod, val) for (i,val) in enumerate(Nw))

    return calcfreqs(Nharmonics, Nw, Nt)
end

"""
    calcfreqsdft(Nharmonics::NTuple{N,Int})

Calculate the dimensions of the DFT in the frequency domain
and the time domain given a tuple of the number of harmonics. Eg. 0,w,2w,3w
would be 3 harmonics. Also calculate the possible modes and their coordinates
in the frequency domain DFT array.

# Arguments
- `Nharmonics`: is a tuple of the number of harmonics to calculate for
    each frequency.

# Returns
- `Frequencies`: A simple structure to hold time and frequency domain
    information for the signal for a single node. See [`Frequencies`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.calcfreqsrdft((1,))
JosephsonCircuits.Frequencies{1}((1,), (2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])

julia> JosephsonCircuits.calcfreqsrdft((2,))
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)])

julia> JosephsonCircuits.calcfreqsrdft((3,))
JosephsonCircuits.Frequencies{1}((3,), (4,), (6,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,), CartesianIndex(4,)], [(0,), (1,), (2,), (3,)])

julia> JosephsonCircuits.calcfreqsrdft((3,3))
JosephsonCircuits.Frequencies{2}((3, 3), (4, 7), (6, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2), CartesianIndex(1, 3), CartesianIndex(2, 3)  …  CartesianIndex(3, 5), CartesianIndex(4, 5), CartesianIndex(1, 6), CartesianIndex(2, 6), CartesianIndex(3, 6), CartesianIndex(4, 6), CartesianIndex(1, 7), CartesianIndex(2, 7), CartesianIndex(3, 7), CartesianIndex(4, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (3, 1), (0, 2), (1, 2)  …  (2, -3), (3, -3), (0, -2), (1, -2), (2, -2), (3, -2), (0, -1), (1, -1), (2, -1), (3, -1)])
```
"""
function calcfreqsdft(Nharmonics::NTuple{N,Int}) where N
    Nw=NTuple{N,Int}(2*val+1 for (i,val) in enumerate(Nharmonics))
    return calcfreqs(Nharmonics, Nw, Nw)
end

"""
    calcfreqs(Nharmonics::NTuple{N,Int}, Nw::NTuple{N,Int}, Nt::NTuple{N,Int}) 

Calculate the dimensions of the DFT or RFDT in the frequency domain
and the time domain given a tuple of the number of harmonics. Eg. 0,w,2w,3w
would be 3 harmonics. Also calculate the possible modes and their coordinates
in the frequency domain RDFT array. See also [`calcfreqsrdft`](@ref)
and [`calcfreqsdft`](@ref).
"""
function calcfreqs(Nharmonics::NTuple{N,Int}, Nw::NTuple{N,Int},
    Nt::NTuple{N,Int}) where N

    # the coordinates of each mixing products
    coords = Array{CartesianIndex{N},1}(undef,prod(Nw))

    # the values of the mixing products in the form of multiples of the 
    # input frequencies
    modes = Vector{NTuple{N,Int}}(undef,prod(Nw))

    # a temporary array for calculating the mixing products
    nvals = zeros(Int, N)

    index = 1
    for i in CartesianIndices(Nw)
        for (ni,nval) in enumerate(i.I)
            if nval <= Nharmonics[ni] + 1
                nvals[ni] = nval-1
            else
                nvals[ni] = -Nw[ni]+nval-1
            end
        end
        coords[index] = i
        modes[index] = NTuple{N,Int}(nvals)
        index+=1
    end

    return Frequencies(Nharmonics, Nw, Nt, coords, modes)
end

"""
    removeconjfreqs(frequencies::Frequencies{N})

Return a new Frequencies struct with the conjugate symmetric terms in the DFT or
RDFT removed.

# Examples
```jldoctest
julia> JosephsonCircuits.removeconjfreqs(JosephsonCircuits.Frequencies{1}((1,), (2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)]))
JosephsonCircuits.Frequencies{1}((1,), (2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])

julia> frequencies = JosephsonCircuits.Frequencies{2}((2,2), (3, 5), (4, 5), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(1, 4), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(1, 5), CartesianIndex(2, 5), CartesianIndex(3, 5)], [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2), (0, -2), (1, -2), (2, -2), (0, -1), (1, -1), (2, -1)]);JosephsonCircuits.removeconjfreqs(frequencies).modes
11-element Vector{Tuple{Int64, Int64}}:
 (0, 0)
 (1, 0)
 (2, 0)
 (0, 1)
 (1, 1)
 (2, 1)
 (0, 2)
 (1, 2)
 (2, 2)
 (1, -2)
 (1, -1)

julia> JosephsonCircuits.removeconjfreqs(JosephsonCircuits.calcfreqsrdft((2,2))).modes
11-element Vector{Tuple{Int64, Int64}}:
 (0, 0)
 (1, 0)
 (2, 0)
 (0, 1)
 (1, 1)
 (2, 1)
 (0, 2)
 (1, 2)
 (2, 2)
 (1, -2)
 (1, -1)
```
"""
function removeconjfreqs(frequencies::Frequencies)
    conjsymdict = conjsym(frequencies)
    return removefreqs(frequencies,collect(values(conjsymdict)))
end

"""
    keepfreqs(frequencies::Frequencies{N},
        keepmodes::AbstractVector{NTuple{N,Int}})

Return a new Frequencies struct with all coordinates and modes except the ones
in keepmodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,2)),[(0,0),(1,0),(0,1),(1,1)])
JosephsonCircuits.Frequencies{2}((2, 2), (3, 5), (4, 5), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2), CartesianIndex(2, 2)], [(0, 0), (1, 0), (0, 1), (1, 1)])

julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,2)),Tuple{Int64,Int64}[])
JosephsonCircuits.Frequencies{2}((2, 2), (3, 5), (4, 5), CartesianIndex{2}[], Tuple{Int64, Int64}[])
```
"""
function keepfreqs(frequencies::Frequencies{N},
    keepmodes::AbstractVector{NTuple{N,Int}}) where N
    Nt = frequencies.Nt
    Nw = frequencies.Nw
    coords = frequencies.coords
    modes = frequencies.modes

    keepmodesdict = Dict{eltype(keepmodes),Int}()
    sizehint!(keepmodesdict,length(keepmodes))

    keepmodessorted = Vector{eltype(modes)}(undef,0)
    sizehint!(keepmodessorted,length(keepmodes))

    keepcoords = Vector{eltype(coords)}(undef,0)
    sizehint!(keepcoords,length(keepmodes))

    for (i,mode) in enumerate(keepmodes)
        keepmodesdict[mode] = i
    end

    for i in eachindex(modes)
        if haskey(keepmodesdict,modes[i])
            push!(keepmodessorted,modes[i])
            push!(keepcoords,coords[i])
        end
    end

    return Frequencies(frequencies.Nharmonics,Nw,Nt,keepcoords,keepmodessorted)
end

"""
    keepfreqs(frequencies::Frequencies{N},
        keepcoords::AbstractVector{CartesianIndex{N}})

Return a new Frequencies struct with all coordinates and modes except the ones
in keepmodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[], Tuple{Int64}[])

julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,)])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,)], [(0,)])
```
"""
function keepfreqs(frequencies::Frequencies{N},
    keepcoords::AbstractVector{CartesianIndex{N}}) where N
    Nt = frequencies.Nt
    Nw = frequencies.Nw
    coords = frequencies.coords
    modes = frequencies.modes

    keepcoordsdict = Dict{eltype(keepcoords),Int}()
    sizehint!(keepcoordsdict,length(keepcoords))

    keepmodes = Vector{eltype(modes)}(undef,0)
    sizehint!(keepmodes,length(keepcoords))

    keepcoordssorted = Vector{eltype(coords)}(undef,0)
    sizehint!(keepcoordssorted,length(keepcoords))

    for (i,coord) in enumerate(keepcoords)
        keepcoordsdict[coord] = i
    end

    for i in eachindex(coords)
        if haskey(keepcoordsdict,coords[i])
            push!(keepmodes,modes[i])
            push!(keepcoordssorted,coords[i])
        end
    end

    return Frequencies(frequencies.Nharmonics,Nw,Nt,keepcoordssorted,keepmodes)
end

"""
    removefreqs(frequencies::Frequencies{N},
        removemodes::AbstractVector{NTuple{N,Int}})

Return a new Frequency struct with the coordinates and modes for the modes in
removemodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[(2,)])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])

julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[(0,),(1,),(2,),(3,)])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[], Tuple{Int64}[])

julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)])
```
"""
function removefreqs(frequencies::Frequencies{N},
    removemodes::AbstractVector{NTuple{N,Int}}) where N
    Nt = frequencies.Nt
    Nw = frequencies.Nw
    coords = frequencies.coords
    modes = frequencies.modes

    # estimate the size of the output
    if length(removemodes) >= length(modes)
        sizeestimate = 0
    else
        sizeestimate = length(modes) - length(removemodes)
    end

    removemodesdict = Dict{eltype(removemodes),Int}()
    sizehint!(removemodesdict,length(removemodes))

    keepmodessorted = Vector{eltype(modes)}(undef,0)
    sizehint!(keepmodessorted,sizeestimate)

    keepcoords = Vector{eltype(coords)}(undef,0)
    sizehint!(keepcoords,sizeestimate)

    for (i,mode) in enumerate(removemodes)
        removemodesdict[mode] = i
    end

    for i in eachindex(modes)
        if !haskey(removemodesdict,modes[i])
            push!(keepmodessorted,modes[i])
            push!(keepcoords,coords[i])
        end
    end

    return Frequencies(frequencies.Nharmonics,Nw,Nt,keepcoords,keepmodessorted)
end

"""
    removefreqs(frequencies::Frequencies{N},
        removecoords::AbstractVector{CartesianIndex{N}})

Return a new Frequency struct with the coordinates and modes for the modes in
removemodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,)])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(3,)], [(1,), (2,)])

julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,),CartesianIndex(2,),CartesianIndex(3,),CartesianIndex(4,)])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[], Tuple{Int64}[])

julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[])
JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)])
```
"""
function removefreqs(frequencies::Frequencies{N},
    removecoords::AbstractVector{CartesianIndex{N}}) where N
    Nt = frequencies.Nt
    Nw = frequencies.Nw
    coords = frequencies.coords
    modes = frequencies.modes

    # estimate the size of the output
    if length(removecoords) >= length(modes)
        sizeestimate = 0
    else
        sizeestimate = length(modes) - length(removecoords)
    end

    removecoordsdict = Dict{eltype(removecoords),Int}()
    sizehint!(removecoordsdict,length(removecoords))

    keepmodes = Vector{eltype(modes)}(undef,0)
    sizehint!(keepmodes,sizeestimate)

    keepcoords = Vector{eltype(coords)}(undef,0)
    sizehint!(keepcoords,sizeestimate)

    for (i,coord) in enumerate(removecoords)
        removecoordsdict[coord] = i
    end

    for i in eachindex(coords)
        if !haskey(removecoordsdict,coords[i])
            push!(keepmodes,modes[i])
            push!(keepcoords,coords[i])
        end
    end

    return Frequencies(frequencies.Nharmonics,Nw,Nt,keepcoords,keepmodes)
end

"""
    truncfreqs(frequencies::Frequencies;
        dc = false, odd = true, even = false, maxintermodorder = Inf)

Return a new Frequencies struct with the coordinates and modes truncated
according to the user specified criteria.

# Examples
```jldoctest
julia> JosephsonCircuits.truncfreqs(JosephsonCircuits.calcfreqsrdft((3,3));maxintermodorder=2).modes
12-element Vector{Tuple{Int64, Int64}}:
 (0, 0)
 (1, 0)
 (2, 0)
 (3, 0)
 (0, 1)
 (1, 1)
 (0, 2)
 (0, 3)
 (0, -3)
 (0, -2)
 (0, -1)
 (1, -1)

julia> JosephsonCircuits.truncfreqs(JosephsonCircuits.calcfreqsrdft((3,3));dc=false,even=false,maxintermodorder=3).modes
10-element Vector{Tuple{Int64, Int64}}:
 (1, 0)
 (3, 0)
 (0, 1)
 (2, 1)
 (1, 2)
 (0, 3)
 (0, -3)
 (1, -2)
 (0, -1)
 (2, -1)

julia> JosephsonCircuits.truncfreqs(JosephsonCircuits.calcfreqsrdft((3,3));maxintermodorder=2)
JosephsonCircuits.Frequencies{2}((3, 3), (4, 7), (6, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(1, 3), CartesianIndex(1, 4), CartesianIndex(1, 5), CartesianIndex(1, 6), CartesianIndex(1, 7), CartesianIndex(2, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (0, 2), (0, 3), (0, -3), (0, -2), (0, -1), (1, -1)])
```
"""
function  truncfreqs(frequencies::Frequencies;
    dc = true, odd = true, even = true, maxintermodorder = Inf)

    coords = frequencies.coords
    modes = frequencies.modes

    # keepmodesdict = Dict{eltype(keepmodes),Int}()
    keepmodes = Vector{eltype(modes)}(undef,0)
    sizehint!(keepmodes,length(modes))

    keepcoords = Vector{eltype(coords)}(undef,0)
    sizehint!(keepcoords,length(modes))

    for (i,nvals) in enumerate(modes)

        # to be returned as a valid frequency, the point has to match the
        # criteria of dc, even, or odd, and either only contain a single frequency
        # or be less than the max intermod order if it contains multiple
        # frequencies. 
        if (
                # test for DC
                (dc && all(==(0),nvals)) ||
                # test for one of the fundamental frequencies
                # sum(abs,nvals) == 1 ||
                # test for even (and not DC)
                (even && mod(sum(abs,nvals),2) == 0 && sum(abs,nvals) > 0) ||
                # test for odd
                (odd && mod(sum(abs,nvals),2) == 1)
            ) && # test for containing only one frequency or less than maxintermodorder
                # (sum( nvals .!== 0) == 1 || sum(abs,nvals) <= maxintermodorder)
                (count(!=(0), nvals) == 1 || sum(abs,nvals) <= maxintermodorder)

            push!(keepcoords,coords[i])
            push!(keepmodes,nvals)
        end
    end

    return Frequencies(frequencies.Nharmonics, frequencies.Nw, frequencies.Nt,
        keepcoords, keepmodes)
end

"""
    calcmodefreqs(w::NTuple{N},modes::Vector{NTuple{N,Int}})

Calculate the frequencies of the modes given a tuple of fundamental frequencies
and a vector of tuples containing the mixing products and harmonics.

# Examples
```jldoctest
julia> @variables wp1 wp2;JosephsonCircuits.calcmodefreqs((wp1, wp2),[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)])
6-element Vector{Num}:
          0
        wp1
       2wp1
        wp2
  wp1 + wp2
 2wp1 + wp2

julia> JosephsonCircuits.calcmodefreqs((1., 1.1),[(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1)])
6-element Vector{Float64}:
 0.0
 1.0
 2.0
 1.1
 2.1
 3.1
```
"""
function calcmodefreqs(w::NTuple{N,Any},modes::Vector{NTuple{N,Int}}) where N
    # generate the frequencies of the modes
    wmodes = Vector{eltype(w)}(undef, length(modes))
    for (i,mode) in enumerate(modes)
        wmodes[i] = dot(w,mode)
    end
    return wmodes
end

"""
    visualizefreqs(w::NTuple{N,Any}, freq::Frequencies{N})

Create a vector or array containing the mixing products for visualization
purposes.

# Examples
```jldoctest
@variables wp1
w = (wp1,)
freq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsrdft((3,)),
        dc=true, odd=true, even=true, maxintermodorder=Inf,
)
JosephsonCircuits.visualizefreqs(w,freq)

# output
4-element Vector{Num}:
    0
  wp1
 2wp1
 3wp1
```
```jldoctest
@variables wp1,wp2
w = (wp1,wp2)
freq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsrdft((3,3)),
        dc=true, odd=true, even=true, maxintermodorder=3,
)
JosephsonCircuits.visualizefreqs(w,freq)

# output
4×7 Matrix{Num}:
    0         wp2        2wp2  3wp2  -3wp2       -2wp2        -wp2
  wp1   wp1 + wp2  wp1 + 2wp2     0      0  wp1 - 2wp2   wp1 - wp2
 2wp1  2wp1 + wp2           0     0      0           0  2wp1 - wp2
 3wp1           0           0     0      0           0           0
```
```jldoctest
w = (1.1,1.2)
freq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsrdft((3,3)),
        dc=true, odd=true, even=true, maxintermodorder=3,
)
JosephsonCircuits.visualizefreqs(w,freq)

# output
4×7 Matrix{Float64}:
 0.0  1.2  2.4  3.6  -3.6  -2.4  -1.2
 1.1  2.3  3.5  0.0   0.0  -1.3  -0.1
 2.2  3.4  0.0  0.0   0.0   0.0   1.0
 3.3  0.0  0.0  0.0   0.0   0.0   0.0
```
"""
function visualizefreqs(w::NTuple{N,Any}, freq::Frequencies{N}) where N
    wmodes = calcmodefreqs(w,freq.modes)

    s = zeros(eltype(wmodes),freq.Nw)
    for (i,coord) in enumerate(freq.coords)
        s[coord] = wmodes[i]
    end

    return s
end

"""
    conjsym(Nw::NTuple{N, Int}, Nt::NTuple{N, Int})

Calculate the conjugate symmetries in the multi-dimensional frequency domain
data.

# Examples
```jldoctest
julia> JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsrdft((2,)))
Dict{CartesianIndex{1}, CartesianIndex{1}}()

julia> JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsdft((2,)))
Dict{CartesianIndex{1}, CartesianIndex{1}} with 2 entries:
  CartesianIndex(2,) => CartesianIndex(5,)
  CartesianIndex(3,) => CartesianIndex(4,)

julia> JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsrdft((2,1)))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 2 entries:
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
  CartesianIndex(3, 2) => CartesianIndex(3, 3)

julia> JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsdft((2,1)))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 7 entries:
  CartesianIndex(2, 3) => CartesianIndex(5, 2)
  CartesianIndex(2, 1) => CartesianIndex(5, 1)
  CartesianIndex(3, 3) => CartesianIndex(4, 2)
  CartesianIndex(3, 1) => CartesianIndex(4, 1)
  CartesianIndex(2, 2) => CartesianIndex(5, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
  CartesianIndex(3, 2) => CartesianIndex(4, 3)

julia> JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsrdft((2,1,1)))
Dict{CartesianIndex{3}, CartesianIndex{3}} with 8 entries:
  CartesianIndex(1, 2, 1) => CartesianIndex(1, 3, 1)
  CartesianIndex(1, 2, 3) => CartesianIndex(1, 3, 2)
  CartesianIndex(1, 2, 2) => CartesianIndex(1, 3, 3)
  CartesianIndex(3, 2, 1) => CartesianIndex(3, 3, 1)
  CartesianIndex(1, 1, 2) => CartesianIndex(1, 1, 3)
  CartesianIndex(3, 2, 3) => CartesianIndex(3, 3, 2)
  CartesianIndex(3, 2, 2) => CartesianIndex(3, 3, 3)
  CartesianIndex(3, 1, 2) => CartesianIndex(3, 1, 3)
```
"""
function conjsym(Nw::NTuple{N, Int}, Nt::NTuple{N, Int}) where N

    conjsymdict = Dict{CartesianIndex{N},CartesianIndex{N}}()

    for k in CartesianIndices(Nw)
        # check that none of the indices are equal but that they are valid
        # indices. 
        if any(k.I .- 1 .!= mod.(Nt .- (k.I .- 1),Nt)) && all(mod.(Nt .- (k.I .- 1),Nt) .< Nw)
            # sort the terms to consistently decide which to call the conjugate
            tmp1 = k.I
            tmp2 = mod.(Nt .- (k.I .- 1),Nt) .+ 1
            if tmp1 < tmp2
                if !haskey(conjsymdict,CartesianIndex(tmp1))
                    conjsymdict[CartesianIndex(tmp1)] = CartesianIndex(tmp2)
                end
            else
                if !haskey(conjsymdict,CartesianIndex(tmp2))
                    conjsymdict[CartesianIndex(tmp2)] = CartesianIndex(tmp1)
                end
            end
        end
    end
    return conjsymdict
end

function conjsym(frequencies::Frequencies{N}) where N
    return conjsym(frequencies.Nw, frequencies.Nt)
end

"""
    printsymmetries(Nw::NTuple{N, Int}, Nt::NTuple{N, Int})

Print the conjugate symmetries in the multi-dimensional DFT or RDFT from the
dimensions of the signal in the frequency domain and the time domain. Negative numbers
indicate that element is the complex conjugate of the corresponding positive
number. A zero indicates that element has no corresponding complex conjugate.

# Examples
```jldoctest
julia> JosephsonCircuits.printsymmetries((3,),(4,))
3-element Vector{Int64}:
 0
 0
 0

julia> JosephsonCircuits.printsymmetries((4,),(4,))
4-element Vector{Int64}:
  0
  1
  0
 -1

julia> JosephsonCircuits.printsymmetries((3,3),(4,3))
3×3 Matrix{Int64}:
 0  1  -1
 0  0   0
 0  2  -2

julia> JosephsonCircuits.printsymmetries((4,3),(4,3))
4×3 Matrix{Int64}:
  0   2  -2
  1   3   5
  0   4  -4
 -1  -5  -3
```
"""
function printsymmetries(Nw::NTuple{N, Int}, Nt::NTuple{N, Int}) where N

    d=conjsym(Nw,Nt)

    z=zeros(Int,Nw)
    i = 1
    for (key,val) in sort(OrderedCollections.OrderedDict(d))
        z[key] = i
        z[val] = -i
        i+=1
    end
    return z
end

"""
    printsymmetries(freq::Frequencies)

See  [`printsymmetries`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.printsymmetries(JosephsonCircuits.calcfreqsrdft((2,)))
3-element Vector{Int64}:
 0
 0
 0

julia> JosephsonCircuits.printsymmetries(JosephsonCircuits.calcfreqsdft((2,)))
5-element Vector{Int64}:
  0
  1
  2
 -2
 -1

julia> JosephsonCircuits.printsymmetries(JosephsonCircuits.calcfreqsrdft((2,2)))
3×5 Matrix{Int64}:
 0  1  3  -3  -1
 0  0  0   0   0
 0  2  4  -4  -2
```
"""
function printsymmetries(freq::Frequencies)
    return printsymmetries(freq.Nw, freq.Nt)
end

"""
    calcindexdict(N::Tuple)

Return a dictionary of Cartesian indices where the Cartesian index is the key
and the index giving the order is the value.

# Examples
```jldoctest
julia> JosephsonCircuits.calcindexdict((2,3))
Dict{CartesianIndex{2}, Int64} with 6 entries:
  CartesianIndex(2, 3) => 6
  CartesianIndex(2, 1) => 2
  CartesianIndex(1, 3) => 5
  CartesianIndex(1, 1) => 1
  CartesianIndex(2, 2) => 4
  CartesianIndex(1, 2) => 3
```
"""
function calcindexdict(N::Tuple)
    d = Dict{CartesianIndex{length(N)},Int}()

    for (i,index) in enumerate(CartesianIndices(N))
        d[index] = i
    end
    return d
end

"""
    calcindexdict(N::Int)

Return a dictionary of Cartesian indices where the Cartesian index is the key
and the index giving the order is the value.

# Examples
```jldoctest
julia> JosephsonCircuits.calcindexdict(3)
Dict{CartesianIndex{1}, Int64} with 3 entries:
  CartesianIndex(2,) => 2
  CartesianIndex(3,) => 3
  CartesianIndex(1,) => 1
```
"""
function calcindexdict(N)
    return calcindexdict(Tuple(N))
end

"""
    calcphiindices(frequencies::Frequencies{N},
        conjsymdict::Dict{CartesianIndex{N},CartesianIndex{N}})

Return the indices which map the elements of the frequency domain vector
to the corresponding elements of the frequency domain array. Also
return the indices `conjsourceindices` whose data should be copied from the
vector to `conjtargetindices` in the array then complex conjugated.

# Arguments
- `Nt`: tuple with dimensions of signal in time domain 
- `dropdict`: dictionary of elements of frequency domain signal to drop where
    the key is the Cartesian index and the value is the value. 

# Returns
- `indexmap`: the indices which map the elements of the frequency domain
    vector elements to the corresponding elements of the frequency domain array
- `conjsourceindices`: data should be copied from here
- `conjtargetindices`: data should be copied to here and conjugated

# Examples
```jldoctest
freq = JosephsonCircuits.Frequencies{2}((4, 3), (5, 7), (8, 7), CartesianIndex{2}[CartesianIndex(2, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(3, 2), CartesianIndex(2, 3), CartesianIndex(1, 4), CartesianIndex(2, 6), CartesianIndex(3, 7)], [(1, 0), (3, 0), (0, 1), (2, 1), (1, 2), (0, 3), (1, -2), (2, -1)])
conjsymdict = Dict{CartesianIndex{2}, CartesianIndex{2}}(CartesianIndex(5, 4) => CartesianIndex(5, 5), CartesianIndex(1, 3) => CartesianIndex(1, 6), CartesianIndex(5, 2) => CartesianIndex(5, 7), CartesianIndex(1, 4) => CartesianIndex(1, 5), CartesianIndex(1, 2) => CartesianIndex(1, 7), CartesianIndex(5, 3) => CartesianIndex(5, 6))
JosephsonCircuits.calcphiindices(freq, conjsymdict)

# output
([2, 4, 6, 8, 12, 16, 27, 33], [6, 16], [31, 21])
```
```jldoctest
freq = JosephsonCircuits.calcfreqsrdft((4,3));
truncfreq = JosephsonCircuits.truncfreqs(freq;dc=false,odd=true,even=false,maxintermodorder=3)
noconjtruncfreq = JosephsonCircuits.removeconjfreqs(truncfreq)
conjsymdict = JosephsonCircuits.conjsym(noconjtruncfreq)
JosephsonCircuits.calcphiindices(noconjtruncfreq,conjsymdict)

# output
([2, 4, 6, 8, 12, 16, 27, 33], [6, 16], [31, 21])
```
"""
function calcphiindices(frequencies::Frequencies{N},
    conjsymdict::Dict{CartesianIndex{N},CartesianIndex{N}}) where N

    modes = frequencies.modes
    coords = frequencies.coords
    Nw = frequencies.Nw
    Nt = frequencies.Nt

    coordsdict = Dict{CartesianIndex{N},Int}()
    sizehint!(coordsdict,length(coords))
    for (i,coord) in enumerate(coords)
        coordsdict[coord] = i
    end

    # empty vector to hold the map between indices in the vector and the
    # matrix
    indexmap = Vector{Int}(undef,length(coords))

    # the index of the element of the N dimensional array in the frequency domain
    # that i should copy to conjtargetindices and take the complex conjugate of. 
    conjsourceindices = Array{Int}(undef,0)

    # the index of the element of the N dimensional array in the frequency domain
    # that i should take the complex conjugate of
    conjtargetindices = Vector{Int}(undef,0)

    # create a dictionary that maps between the CartesianIndex coordinates
    # and the index in the array at which they occur. 
    carttoint = calcindexdict(Nw)

    # generate the index maps to convert between the vector
    # and matrix. Loop over the modes and coords vectors to keep the order
    # between.
    for (i,coord) in enumerate(coords)
        # push!(indexmap,carttoint[coord])
        indexmap[i] = carttoint[coord]
        if haskey(conjsymdict,coord)
            push!(conjsourceindices,carttoint[coord])
            push!(conjtargetindices,carttoint[conjsymdict[coord]])
        end
    end
    return indexmap, conjsourceindices, conjtargetindices
end

"""
    phivectortomatrix!(phivector::AbstractVector,phimatrix::AbstractArray,
        indexmap::Vector{Int},conjsourceindices::Vector{Int},
        conjtargetindices::Vector{Int},Nbranches::Int)

The harmonic balance method requires a vector with all of the conjugate symmetric
terms removed and potentially other terms dropped if specified by the user (
for example, intermodulation products which are not of interest) whereas the
Fourier transform operates on multidimensional arrays with the proper
conjugate symmetries and with dropped terms set to zero. This function converts
a vector to an array with the above properties.

# Examples
```jldoctest
freqindexmap = [2, 4, 6, 8, 12, 16, 27, 33]
conjsourceindices = [16, 6]
conjtargetindices = [21, 31]
Nbranches = 1

phivector = 1im.*Complex.(1:Nbranches*length(freqindexmap));
phimatrix=zeros(Complex{Float64},5,7,1)

JosephsonCircuits.phivectortomatrix!(phivector,
    phimatrix,
    freqindexmap,
    conjsourceindices,
    conjtargetindices,
    Nbranches,
)
phimatrix

# output
5×7×1 Array{ComplexF64, 3}:
[:, :, 1] =
 0.0+0.0im  0.0+3.0im  0.0+0.0im  0.0+6.0im  0.0-6.0im  0.0+0.0im  0.0-3.0im
 0.0+1.0im  0.0+0.0im  0.0+5.0im  0.0+0.0im  0.0+0.0im  0.0+7.0im  0.0+0.0im
 0.0+0.0im  0.0+4.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+8.0im
 0.0+2.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im  0.0+0.0im
```
"""
function phivectortomatrix!(phivector::AbstractVector, phimatrix::AbstractArray,
    indexmap::Vector{Int}, conjsourceindices::Vector{Int},
    conjtargetindices::Vector{Int}, Nbranches::Int)

    if length(indexmap)*Nbranches != length(phivector)
        throw(DimensionMismatch("Unexpected length for phivector"))
    end

    if length(phivector) == 0
        Nvector = 0
    else
        Nvector = length(phivector)÷ Nbranches
    end

    # Nvector = length(phivector)÷ Nbranches
    Nmatrix = prod(size(phimatrix)[1:end-1])

    # fill the matrix with zeros
    fill!(phimatrix,0)

    # copy the mixing products from the vector to the matrix
    for i in 1:Nbranches
        for j in 1:length(indexmap)
            phimatrix[indexmap[j]+(i-1)*Nmatrix] = phivector[j+(i-1)*Nvector]
        end
    end

    # for the complex conjugates, copy and complex conjugate them from source
    # index of the matrx to the destination index of the matrix.
    for i in 1:Nbranches
        for j in 1:length(conjtargetindices)
            phimatrix[conjtargetindices[j]+ (i-1)*Nmatrix] = conj(phimatrix[conjsourceindices[j]+ (i-1)*Nmatrix])
        end
    end
    return nothing
end

"""
    phimatrixtovector!(phivector::Vector, phimatrix::Array,
        indexmap::Vector{Int}, conjsourceindices::Vector{Int},
        conjtargetindices::Vector{Int}, Nbranches::Int)

The harmonic balance method requires a vector with all of the conjugate symmetric
terms removed and potentially other terms dropped if specified by the user (
for example, intermodulation products which are not of interest) whereas the
Fourier transform operates on multidimensional arrays with the proper
conjugate symmetries and with dropped terms set to zero. This function converts
an array to a vector with the above properties.

# Examples
```jldoctest
freqindexmap = [2, 4, 6, 8, 12, 16, 27, 33]
conjsourceindices = [16, 6]
conjtargetindices = [21, 31]
Nbranches = 1

phivector = zeros(Complex{Float64}, Nbranches*length(freqindexmap))
phimatrix = [0.0 + 0.0im 0.0 + 3.0im 0.0 + 0.0im 0.0 + 6.0im 0.0 - 6.0im 0.0 + 0.0im 0.0 - 3.0im; 0.0 + 1.0im 0.0 + 0.0im 0.0 + 5.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 7.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 4.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 8.0im; 0.0 + 2.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;;;]

JosephsonCircuits.phimatrixtovector!(phivector,
    phimatrix,
    freqindexmap,
    conjsourceindices,
    conjtargetindices,
    Nbranches,
)
phivector

# output
8-element Vector{ComplexF64}:
 0.0 + 1.0im
 0.0 + 2.0im
 0.0 + 3.0im
 0.0 + 4.0im
 0.0 + 5.0im
 0.0 + 6.0im
 0.0 + 7.0im
 0.0 + 8.0im
```
"""
function phimatrixtovector!(phivector::AbstractVector, phimatrix::AbstractArray,
    indexmap::Vector{Int}, conjsourceindices::Vector{Int},
    conjtargetindices::Vector{Int}, Nbranches::Int)

    # Nvector = length(phivector)÷ Nbranches

    if length(phivector) == 0
        Nvector = 0
    else
        Nvector = length(phivector)÷ Nbranches
    end

    Nmatrix = prod(size(phimatrix)[1:end-1])

    # fill the vector with zeros
    fill!(phivector,0)

    if length(indexmap)*Nbranches != length(phivector)
        throw(DimensionMismatch("Unexpected length for phivector"))
    end

    for i in 1:Nbranches
        for j in 1:length(indexmap)
            phivector[j+(i-1)*Nvector] = phimatrix[indexmap[j]+(i-1)*Nmatrix]
        end
    end
    return nothing
end

"""
    applynl(am::Array{Complex{Float64}}, f::Function)

Perform the inverse discrete Fourier transform on an array `am` of complex
frequency domain data, apply the function `f` in the time domain, then perform
the discrete Fourier transform to return to the frequency domain. Apply the
Fourier transform on all but the last dimensions. See also [`applynl!`](@ref)
and [`plan_applynl`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.applynl([[0, 0.2+0.0im, 0, 0];;],cos)
4×1 Matrix{ComplexF64}:
   0.9603980498951228 + 0.0im
                  0.0 + 0.0im
 -0.01966852794611884 + 0.0im
                  0.0 + 0.0im

julia> JosephsonCircuits.applynl([[0, 0.2+0.0im];;],cos)
2×1 Matrix{ComplexF64}:
   0.9603980498951228 + 0.0im
 -0.01966852794611884 + 0.0im

julia> JosephsonCircuits.applynl([0.0 + 0.0im 0.45 + 0.0im 0.45 + 0.0im; 0.55 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;;;],sin)
3×3×1 Array{ComplexF64, 3}:
[:, :, 1] =
 -0.0209812+0.0im   0.295151+0.0im   0.295151+0.0im
   0.359826+0.0im  -0.041417+0.0im  -0.041417+0.0im
 0.00788681+0.0im  -0.110947+0.0im  -0.110947+0.0im
```
"""
function applynl(fd::Array{Complex{Float64}}, f)

    td, irfftplan, rfftplan = plan_applynl(fd)
    fdcopy = copy(fd)

    applynl!(fdcopy, td, f, irfftplan, rfftplan)

    return fdcopy
end

"""
    plan_applynl(fd::Array{Complex{T}})

Creates an empty time domain data array and the inverse and forward plans
for the RFFT of an array of frequency domain data. See also [`applynl!`](@ref).

"""
function plan_applynl(fd::Array{Complex{T}}) where T

    #choose the number of time points based on the number of fourier
    #coefficients
    # changed to below because i wasn't using enough points when Nmodes=1.
    # the results contained only real values. 
    sizefd = size(fd)
    if sizefd[1] == 2
        stepsperperiod = 2*sizefd[1]-1
    else
        stepsperperiod = 2*sizefd[1]-2
    end

    # generate the time domain array with the appropriate dimensions
    td = Array{T}(
        undef,
        NTuple{length(sizefd),Int}(ifelse(i == 1, stepsperperiod, val) for (i,val) in enumerate(sizefd)),
    )

    # make the irfft plan
    irfftplan = FFTW.plan_irfft(fd,stepsperperiod,1:length(size(fd))-1; flags=FFTW.ESTIMATE, timelimit=Inf)

    # make the rfft plan
    rfftplan = FFTW.plan_rfft(td,1:length(size(fd))-1; flags=FFTW.ESTIMATE, timelimit=Inf)

    return td, irfftplan, rfftplan
end

"""
    applynl!(fd::Array{Complex{T}}, td::Array{T}, f::Function, irfftplan,
        rfftplan)

Apply the nonlinear function f to the frequency domain data by transforming
to the time domain, applying the function, then transforming back to the
frequency domain, overwriting the contents of fd and td in the process. We
use plans for the forward and reverse RFFT prepared by [`plan_applynl`](@ref).

# Examples
```jldoctest
fd=ones(Complex{Float64},3,2)
td, irfftplan, rfftplan = JosephsonCircuits.plan_applynl(fd)
JosephsonCircuits.applynl!(fd, td, cos, irfftplan, rfftplan)
fd

# output
3×2 Matrix{ComplexF64}:
  0.586589+0.0im   0.586589+0.0im
 -0.413411+0.0im  -0.413411+0.0im
 -0.413411+0.0im  -0.413411+0.0im
```
"""
function applynl!(fd::Array{Complex{T}}, td::Array{T}, f, irfftplan,
    rfftplan) where T

    #transform to the time domain
    mul!(td, irfftplan, fd)

    # normalize the fft
    normalization = prod(size(td)[1:end-1])
    invnormalization = 1/normalization

    #apply the nonlinear function
    for i in eachindex(td)
        td[i] = f(td[i]*normalization)
    end

    # transform to the frequency domain
    mul!(fd, rfftplan, td)

    # normalize
    for i in eachindex(fd)
        fd[i] = fd[i]*invnormalization
    end

    return nothing
end

"""
    hbmatind(truncfrequencies::Frequencies{N})

Returns a matrix describing which indices of the frequency domain matrix
(from the RFFT) to pull out and use in the harmonic balance matrix. A negative
index means we take the complex conjugate of that element. A zero index means
that term is not present, so skip it. The harmonic balance matrix describes
the coupling between different frequency modes.

# Examples
```jldoctest
julia> freq = JosephsonCircuits.calcfreqsrdft((5,));JosephsonCircuits.hbmatind(JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(freq;dc=false,odd=true,even=false,maxintermodorder=2)))[2]
3×3 Matrix{Int64}:
 1  -3  -5
 3   1  -3
 5   3   1

julia> freq = JosephsonCircuits.calcfreqsrdft((3,));JosephsonCircuits.hbmatind(JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(freq;dc=true,odd=true,even=true,maxintermodorder=2)))[2]
4×4 Matrix{Int64}:
 1  -2  -3  -4
 2   1  -2  -3
 3   2   1  -2
 4   3   2   1

julia> freq = JosephsonCircuits.calcfreqsrdft((2,2));JosephsonCircuits.hbmatind(JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(freq;dc=true,odd=true,even=true,maxintermodorder=2)))[1]
7×7 Matrix{Tuple{Int64, Int64}}:
 (0, 0)   (-1, 0)  (-2, 0)   (0, -1)  (-1, -1)  (0, -2)  (-1, 1)
 (1, 0)   (0, 0)   (-1, 0)   (1, -1)  (0, -1)   (1, -2)  (0, 1)
 (2, 0)   (1, 0)   (0, 0)    (2, -1)  (1, -1)   (2, -2)  (1, 1)
 (0, 1)   (-1, 1)  (-2, 1)   (0, 0)   (-1, 0)   (0, -1)  (-1, 2)
 (1, 1)   (0, 1)   (-1, 1)   (1, 0)   (0, 0)    (1, -1)  (0, 2)
 (0, 2)   (-1, 2)  (-2, 2)   (0, 1)   (-1, 1)   (0, 0)   (-1, 3)
 (1, -1)  (0, -1)  (-1, -1)  (1, -2)  (0, -2)   (1, -3)  (0, 0)

julia> freq = JosephsonCircuits.calcfreqsrdft((2,2));JosephsonCircuits.hbmatind(JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(freq;dc=true,odd=true,even=true,maxintermodorder=2)))[2]
7×7 Matrix{Int64}:
  1   -2   -3  13   -5  10  -14
  2    1   -2  14   13  11    4
  3    2    1  15   14  12    5
  4  -14  -15   1   -2  13  -11
  5    4  -14   2    1  14    7
  7  -11  -12   4  -14   1    0
 14   13   -5  11   10   0    1
```
"""
function hbmatind(truncfrequencies::Frequencies{N}) where N
    frequencies = calcfreqs(truncfrequencies.Nharmonics,
        truncfrequencies.Nw, truncfrequencies.Nt)
    return hbmatind(frequencies, truncfrequencies)
end

"""
    hbmatind(frequencies::Frequencies{N},
        truncfrequencies::Frequencies{N})

Returns a matrix describing which indices of the frequency domain matrix
(from the RFFT or FFT) to pull out and use in the harmonic balance matrix.
A negative index means we take the complex conjugate of that element. A zero
index means that term is not present, so skip it. The harmonic balance matrix
describes the coupling between different frequency modes.

# Examples
```jldoctest
pumpfreq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsrdft((4,)))
signalfreq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsdft((4,));
    dc=false,odd=true,even=false,maxintermodorder=2,
)
JosephsonCircuits.hbmatind(pumpfreq, signalfreq)[2]

# output
4×4 Matrix{Int64}:
  1  -3  5   3
  3   1  0   5
 -5   0  1  -3
 -3  -5  3   1
```
"""
function hbmatind(frequencies::Frequencies{N},
    truncfrequencies::Frequencies{N}) where N

    modes = frequencies.modes
    truncmodes = truncfrequencies.modes

    # this is calculating the frequency domain input output relations
    # first calculate this in terms of the modes
    Amatrixmodes = Matrix{NTuple{N,Int}}(undef,length(truncmodes),length(truncmodes))
    for i in 1:length(truncmodes)
        for j in 1:length(truncmodes)
            mode = NTuple{N,Int}(truncmodes[i][k]-truncmodes[j][k] for k in 1:length(truncmodes[i]))
            Amatrixmodes[i,j] = mode
        end
    end

    # now i need to find the keys that are in the rfft matrix and their
    # locations
    modesdict = Dict{eltype(modes),Int}()
    for (i,mode) in enumerate(modes)
        modesdict[mode] = i
    end

    # now we want to know where in the dft matrix we should pull
    # these modes. to do these, we use the un-truncated frequencies struct.
    Amatrixindices = zeros(Int,length(truncmodes),length(truncmodes))
    for (i,mode) in enumerate(Amatrixmodes)
        conjmode = NTuple{N}(-m for m in mode)
        if haskey(modesdict,mode)
            Amatrixindices[i] = modesdict[mode]
        elseif haskey(modesdict,conjmode)
            Amatrixindices[i] = -modesdict[conjmode]
        end
    end

    return Amatrixmodes, Amatrixindices
end
