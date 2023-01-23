
"""
    Frequencies(Nw::NTuple{N,Int}, Nt::NTuple{N,Int},
        coords::Vector{CartesianIndex{N}}, modes::Vector{NTuple{N,Int})

A simple structure to hold time and frequency domain information for the
signals. See also [`calcfreqsrdft`](@ref).

# Fields
- `Nw::NTuple{N,Int}`: The dimensions of the frequency domain signal for a
    single node.
- `Nt::NTuple{N,Int}`: The dimensions of the time domain signal for a single
    node.
- `coords::Vector{CartesianIndex{N}}`: The coordinates of each mixing products.
- `modes::Vector{NTuple{N,Int}}`: The mode indices of each mixing product, eg.
     (0,0), (1,0), (2,1).

# Examples
```jldoctest
Nw = (3, 3)
Nt = (4, 3)
coords = CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3)]
modes = [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, -1), (1, -1), (2, -1)]
JosephsonCircuits.Frequencies(Nw,Nt,coords,modes)

# output
JosephsonCircuits.Frequencies{2}((3, 3), (4, 3), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3)], [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, -1), (1, -1), (2, -1)])
```
"""
struct Frequencies{N}
    Nw::NTuple{N,Int}
    Nt::NTuple{N,Int}
    coords::Vector{CartesianIndex{N}}
    modes::Vector{NTuple{N,Int}}
end

"""
    calcfreqsrdft(Nharmonics::NTuple{N,Int})

Calculate the dimensions of the RFDT in the frequency domain
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
JosephsonCircuits.Frequencies{1}((2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])

julia> JosephsonCircuits.calcfreqsrdft((2,))
JosephsonCircuits.Frequencies{1}((3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)])

julia> JosephsonCircuits.calcfreqsrdft((3,))
JosephsonCircuits.Frequencies{1}((4,), (6,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,), CartesianIndex(4,)], [(0,), (1,), (2,), (3,)])

julia> JosephsonCircuits.calcfreqsrdft((3,3))
JosephsonCircuits.Frequencies{2}((4, 7), (6, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2), CartesianIndex(1, 3), CartesianIndex(2, 3)  …  CartesianIndex(3, 5), CartesianIndex(4, 5), CartesianIndex(1, 6), CartesianIndex(2, 6), CartesianIndex(3, 6), CartesianIndex(4, 6), CartesianIndex(1, 7), CartesianIndex(2, 7), CartesianIndex(3, 7), CartesianIndex(4, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (3, 1), (0, 2), (1, 2)  …  (2, -3), (3, -3), (0, -2), (1, -2), (2, -2), (3, -2), (0, -1), (1, -1), (2, -1), (3, -1)])
```
"""
function calcfreqsrdft(Nharmonics::NTuple{N,Int}) where N
    
    n = Nharmonics .+ 1
    
    # double the size of all but the first n because 
    # only the first axis of a multi-dimensional rfft has
    # only positive frequencies.
    Nw=NTuple{N,Int}(ifelse(i == 1, val, 2*val-1) for (i,val) in enumerate(n))
    
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

    # the coordinates of each mixing products
    coords = Array{CartesianIndex{length(n)},1}(undef,prod(Nw))

    # the values of the mixing products in the form of multiples of the 
    # input frequencies
    modes = Vector{NTuple{length(n),Int}}(undef,prod(Nw))

    # a temporary array for calculating the mixing products
    nvals = zeros(Int,length(n))

    index = 1
    for i in CartesianIndices(Nw)
        for (ni,nval) in enumerate(i.I)
            if nval <= n[ni]
                nvals[ni] = nval-1
            else
                nvals[ni] = -Nw[ni]+nval-1
            end
        end
        coords[index] = i
        modes[index] = NTuple{N,Int}(nvals)
        index+=1
    end

    return Frequencies(Nw,Nt,coords,modes)
end

"""
    removeconjfreqsrdft(frequencies::JosephsonCircuits.Frequencies{N}) where N

Return a new Frequency struct with the conjugate symmetric terms in the RDFT
removed.

# Examples
```jldoctest
julia> JosephsonCircuits.removeconjfreqsrdft(JosephsonCircuits.Frequencies{1}((2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)]))
JosephsonCircuits.Frequencies{1}((2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])

julia> frequencies = JosephsonCircuits.Frequencies{2}((3, 5), (4, 5), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(1, 4), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(1, 5), CartesianIndex(2, 5), CartesianIndex(3, 5)], [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2), (0, -2), (1, -2), (2, -2), (0, -1), (1, -1), (2, -1)]);JosephsonCircuits.removeconjfreqsrdft(frequencies).modes
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

julia> JosephsonCircuits.removeconjfreqsrdft(JosephsonCircuits.calcfreqsrdft((2,2))).modes
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
function removeconjfreqsrdft(frequencies::JosephsonCircuits.Frequencies)
    conjsymdict = conjsymrdft(frequencies.Nt)
    return removefreqs(frequencies,collect(values(conjsymdict)))
end

"""
    keepfreqs(frequencies::JosephsonCircuits.Frequencies{N},
        keepmodes::AbstractVector{NTuple{N,Int}}) where N

Return a new Frequency struct with all coordinates and modes except the ones
in keepmodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,2)),[(0,0),(1,0),(0,1),(1,1)])
JosephsonCircuits.Frequencies{2}((3, 5), (4, 5), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2), CartesianIndex(2, 2)], [(0, 0), (1, 0), (0, 1), (1, 1)])

julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,2)),Tuple{Int64,Int64}[])
JosephsonCircuits.Frequencies{2}((3, 5), (4, 5), CartesianIndex{2}[], Tuple{Int64, Int64}[])
```
"""
function keepfreqs(frequencies::JosephsonCircuits.Frequencies{N},
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

    return Frequencies(Nw,Nt,keepcoords,keepmodessorted)
end


"""
    keepfreqs(frequencies::JosephsonCircuits.Frequencies{N},
        keepcoords::AbstractVector{CartesianIndex{N}}) where N

Return a new Frequency struct with all coordinates and modes except the ones
in keepmodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[])
JosephsonCircuits.Frequencies{1}((3,), (4,), CartesianIndex{1}[], Tuple{Int64}[])

julia> JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,)])
JosephsonCircuits.Frequencies{1}((3,), (4,), CartesianIndex{1}[CartesianIndex(1,)], [(0,)])
```
"""
function keepfreqs(frequencies::JosephsonCircuits.Frequencies{N},
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

    return Frequencies(Nw,Nt,keepcoordssorted,keepmodes)
end


"""
    removefreqs(frequencies::JosephsonCircuits.Frequencies{N},
        removemodes::AbstractVector{NTuple{N,Int}}) where N

Return a new Frequency struct with the coordinates and modes for the modes in
removemodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[(2,)])
JosephsonCircuits.Frequencies{1}((3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])

julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[])
JosephsonCircuits.Frequencies{1}((3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)])
```
"""
function removefreqs(frequencies::JosephsonCircuits.Frequencies{N},
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

    return Frequencies(Nw,Nt,keepcoords,keepmodessorted)
end

"""
    removefreqs(frequencies::JosephsonCircuits.Frequencies{N},
        removecoords::AbstractVector{CartesianIndex{N}}) where N

Return a new Frequency struct with the coordinates and modes for the modes in
removemodes removed.

# Examples
```jldoctest
julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[])
JosephsonCircuits.Frequencies{1}((3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)])

julia> JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,)])
JosephsonCircuits.Frequencies{1}((3,), (4,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(3,)], [(1,), (2,)])
```
"""
function removefreqs(frequencies::JosephsonCircuits.Frequencies{N},
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

    return Frequencies(Nw,Nt,keepcoords,keepmodes)
end

"""
    truncfreqsrdft(frequencies::JosephsonCircuits.Frequencies;
        dc = false, odd = true, even = false, maxintermodorder = Inf)

Return a new Frequency struct with the coordinates and modes truncated
according to the user specified criteria.

# Examples
```jldoctest
julia> JosephsonCircuits.truncfreqsrdft(JosephsonCircuits.calcfreqsrdft((3,3));maxintermodorder=2).modes
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

julia> JosephsonCircuits.truncfreqsrdft(JosephsonCircuits.calcfreqsrdft((3,3));dc=false,even=false,maxintermodorder=3).modes
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

julia> JosephsonCircuits.truncfreqsrdft(JosephsonCircuits.calcfreqsrdft((3,3));maxintermodorder=2)
JosephsonCircuits.Frequencies{2}((4, 7), (6, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(1, 3), CartesianIndex(1, 4), CartesianIndex(1, 5), CartesianIndex(1, 6), CartesianIndex(1, 7), CartesianIndex(2, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (0, 2), (0, 3), (0, -3), (0, -2), (0, -1), (1, -1)])
```
"""
function  truncfreqsrdft(frequencies::JosephsonCircuits.Frequencies;
    dc = true, odd = true, even = true, maxintermodorder = Inf)

    Nt = frequencies.Nt
    Nw = frequencies.Nw
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
                sum(abs,nvals) == 1 ||
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

    return Frequencies(Nw,Nt,keepcoords,keepmodes)
end


"""
    calcfrequencies(Nharmonics::Tuple; maxintermodorder = Inf,
        dc = true, even = true, odd = true)

Given the number of harmonics for each frequency
`Nharmonics`, return the dimensions of the frequency domain representation of
the pumps with harmonics and intermodulation products up to `maxintermodorder`.

# Arguments
- `Nharmonics::Tuple`: is a tuple of the number of harmonics to calculate for
    each frequency.

# Keywords
- `maxintermodorder::Number = Inf`: the maximum intermod order
- `dc::Bool = true`: whether to include the DC term (zero frequency) (0*w)
- `even::Bool = true`: whether to include even order terms (2*w, 4*w, 6*w, etc)
- `odd::Bool = true`: whether to include odd order terms (w, 3*w, 5*w, etc)

# Returns
- `Nw`: tuple with dimensions of signal in frequency domain
- `coords`: vector of Cartesian indices of harmonics and intermods we are keeping
- `values`: vector of values of harmonics and intermods we are keeping (also
    potentially including conjugate symmetric terms)
- `dropcoords`: vector of Cartesian indices for intermods we have dropped
- `dropvalues`: vector of intermod values we have dropped

# Examples
```jldoctest
Nharmonics = (2,3)
maxintermodorder = 2
Nw,coords,values,dropcoords,dropvalues = JosephsonCircuits.calcfrequencies(Nharmonics,
    maxintermodorder=maxintermodorder,dc=false,even=false,odd=true);
println(Nw)
println(coords)
println(values)
println(dropcoords)
println(dropvalues)

# output
(3, 7)
CartesianIndex{2}[CartesianIndex(2, 1), CartesianIndex(1, 2), CartesianIndex(1, 4), CartesianIndex(1, 5), CartesianIndex(1, 7)]
[(1, 0), (0, 1), (0, 3), (0, -3), (0, -1)]
CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(3, 1), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(2, 5), CartesianIndex(3, 5), CartesianIndex(1, 6), CartesianIndex(2, 6), CartesianIndex(3, 6), CartesianIndex(2, 7), CartesianIndex(3, 7)]
[(0, 0), (2, 0), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2), (1, 3), (2, 3), (1, -3), (2, -3), (0, -2), (1, -2), (2, -2), (1, -1), (2, -1)]
```
"""
function calcfrequencies(Nharmonics::Tuple; maxintermodorder::Number = Inf,
    dc::Bool = true, even::Bool = true, odd::Bool = true)

    # if length(w) !== length(Nharmonics)
    #     error("Each frequency must have a number of harmonics.")
    # end
    
    n = Nharmonics .+ 1
    
    # double the size of all but the first n because 
    # only the first axis of a multi-dimensional rfft has
    # only positive frequencies.
    Nw=NTuple{length(n),Int}(ifelse(i == 1, val, 2*val-1) for (i,val) in enumerate(n))
    
    # store the coordinates in an array of CartesianIndex structures
    # not sure if i want to use an array of tuples or cartesianindices
    coords = Array{CartesianIndex{length(n)},1}(undef,0)
    dropcoords = Array{CartesianIndex{length(n)},1}(undef,0)

    # store the values in an array of the same type as w
    values = Vector{NTuple{length(n),Int}}(undef,0)
    dropvalues = Vector{NTuple{length(n),Int}}(undef,0)
    
    nvals = zeros(eltype(n),length(n))

    calcfrequencies!(coords, values, dropcoords, dropvalues, n, nvals, Nw, dc,
        even, odd, maxintermodorder)

    return Nw,coords,values,dropcoords,dropvalues
end

"""
    calcfrequencies!(coords, values, dropcoords, dropvalues, n, nvals,
        Nw, dc, even, odd, maxintermodorder)

See the description for [`calcfrequencies`](@ref).
"""
function calcfrequencies!(coords, values, dropcoords, dropvalues, n, nvals,
    Nw, dc, even, odd, maxintermodorder)

    for i in CartesianIndices(Nw)
        for (ni,nval) in enumerate(i.I)
            if nval <= n[ni]
                nvals[ni] = nval-1
            else
                nvals[ni] = -Nw[ni]+nval-1
            end
        end

        # to be returned as a valid frequency, the point has to match the
        # criteria of dc, even, or odd, and either only contain a single frequency
        # or be less than the max intermod order if it contains multiple
        # frequencies. 
        if (
                # test for DC
                (dc && all(==(0),nvals)) ||
                # test for one of the fundamental frequencies
                sum(abs,nvals) == 1 ||
                # test for even (and not DC)
                (even && mod(sum(abs,nvals),2) == 0 && sum(abs,nvals) > 0) ||
                # test for odd
                (odd && mod(sum(abs,nvals),2) == 1)
            ) && # test for containing only one frequency or less than maxintermodorder
                # (sum( nvals .!== 0) == 1 || sum(abs,nvals) <= maxintermodorder)
                (count(!=(0), nvals) == 1 || sum(abs,nvals) <= maxintermodorder)

            push!(coords,i)
            # push!(values,dot(w,nvals))
            push!(values,Tuple(nvals))
        else
            push!(dropcoords,i)
            # push!(dropvalues,dot(w,nvals))
            push!(dropvalues,Tuple(nvals))

        end
    end

    return nothing
end



"""
    calcfrequencies(w::Tuple, Nharmonics::Tuple; maxintermodorder = Inf,
        dc = true, even = true, odd = true)

Given the frequencies `w` and the number of harmonics for each frequency
`Nharmonics`, return the dimensions of the frequency domain representation of
the pumps with harmonics and intermodulation products up to `maxintermodorder`.

# Arguments
- `w::Tuple`: is a tuple of angular frequencies such as (2*pi*5.0e9, 2*pi*6.0e9)
- `Nharmonics::Tuple`: is a tuple of the number of harmonics to calculate for
    each frequency.

# Keywords
- `maxintermodorder::Number = Inf`: the maximum intermod order
- `dc::Bool = true`: whether to include the DC term (zero frequency) (0*w)
- `even::Bool = true`: whether to include even order terms (2*w, 4*w, 6*w, etc)
- `odd::Bool = true`: whether to include odd order terms (w, 3*w, 5*w, etc)

# Returns
- `Nw`: tuple with dimensions of signal in frequency domain
- `coords`: vector of Cartesian indices of harmonics and intermods we are keeping
- `values`: vector of values of harmonics and intermods we are keeping (also
    potentially including conjugate symmetric terms)
- `dropcoords`: vector of Cartesian indices for intermods we have dropped
- `dropvalues`: vector of intermod values we have dropped

# Examples
```jldoctest
@variables w1,w2
w = (w1,w2)
Nharmonics = (2,3)
maxintermodorder = 2
Nw,coords,values,dropcoords,dropvalues = JosephsonCircuits.calcfrequencies(w,Nharmonics,
    maxintermodorder=maxintermodorder,dc=false,even=false,odd=true);
println(Nw)
println(coords)
println(values)
println(dropcoords)
println(dropvalues)

# output
(3, 7)
CartesianIndex{2}[CartesianIndex(2, 1), CartesianIndex(1, 2), CartesianIndex(1, 4), CartesianIndex(1, 5), CartesianIndex(1, 7)]
Num[w1, w2, 3w2, -3w2, -w2]
CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(3, 1), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(2, 5), CartesianIndex(3, 5), CartesianIndex(1, 6), CartesianIndex(2, 6), CartesianIndex(3, 6), CartesianIndex(2, 7), CartesianIndex(3, 7)]
Num[0, 2w1, w1 + w2, w2 + 2w1, 2w2, w1 + 2w2, 2w1 + 2w2, w1 + 3w2, 2w1 + 3w2, w1 - 3w2, 2w1 - 3w2, -2w2, w1 - 2w2, 2w1 - 2w2, w1 - w2, 2w1 - w2]
```
"""
function calcfrequencies(w::Tuple, Nharmonics::Tuple; maxintermodorder::Number = Inf,
    dc::Bool = true, even::Bool = true, odd::Bool = true)

    if length(w) !== length(Nharmonics)
        error("Each frequency must have a number of harmonics.")
    end
    
    n = Nharmonics .+ 1
    
    # double the size of all but the first n because 
    # only the first axis of a multi-dimensional rfft has
    # only positive frequencies.
    Nw=NTuple{length(n),Int}(ifelse(i == 1, val, 2*val-1) for (i,val) in enumerate(n))
    
    # store the coordinates in an array of CartesianIndex structures
    # not sure if i want to use an array of tuples or cartesianindices
    coords = Array{CartesianIndex{length(n)},1}(undef,0)
    dropcoords = Array{CartesianIndex{length(n)},1}(undef,0)

    # store the values in an array of the same type as w
    values = Array{eltype(w),1}(undef,0)
    dropvalues = Array{eltype(w),1}(undef,0)
    
    nvals = zeros(eltype(n),length(n))

    calcfrequencies!(coords, values, dropcoords, dropvalues, w, n, nvals, Nw, dc,
        even, odd, maxintermodorder)

    return Nw,coords,values,dropcoords,dropvalues
end

"""
    calcfrequencies(w::Number, Nharmonics::Number; maxintermodorder=Inf,
        dc=true, even=true, odd=true)

See the description for [`calcfrequencies`](@ref).

# Examples
```jldoctest
@variables w1
w = w1
Nharmonics = 2
maxintermodorder = 2
Nw,coords,values,dropcoords,dropvalues = JosephsonCircuits.calcfrequencies(w,Nharmonics,
    maxintermodorder=maxintermodorder,dc=false,even=false,odd=true);
println(Nw)
println(coords)
println(values)
println(dropcoords)
println(dropvalues)

# output
(3,)
CartesianIndex{1}[CartesianIndex(2,)]
Num[w1]
CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(3,)]
Num[0, 2w1]
```
"""
function calcfrequencies(w::Number, Nharmonics::Number; maxintermodorder::Number = Inf,
    dc::Bool = true, even::Bool = true, odd::Bool = true)
    return calcfrequencies(Tuple(w), Tuple(Nharmonics),
        maxintermodorder = maxintermodorder, dc = dc, even = even, odd = odd)
end

"""
    calcfrequencies!(coords, values, dropcoords, dropvalues, w, n, nvals,
        Nw, dc, even, odd, maxintermodorder)

See the description for [`calcfrequencies`](@ref).
"""
function calcfrequencies!(coords, values, dropcoords, dropvalues, w, n, nvals,
    Nw, dc, even, odd, maxintermodorder)

    for i in CartesianIndices(Nw)
        for (ni,nval) in enumerate(i.I)
            if nval <= n[ni]
                nvals[ni] = nval-1
            else
                nvals[ni] = -Nw[ni]+nval-1
            end
        end

        # to be returned as a valid frequency, the point has to match the
        # criteria of dc, even, or odd, and either only contain a single frequency
        # or be less than the max intermod order if it contains multiple
        # frequencies. 
        if (
                # test for DC
                (dc && all(==(0),nvals)) ||
                # test for one of the fundamental frequencies
                sum(abs,nvals) == 1 ||
                # test for even (and not DC)
                (even && mod(sum(abs,nvals),2) == 0 && sum(abs,nvals) > 0) ||
                # test for odd
                (odd && mod(sum(abs,nvals),2) == 1)
            ) && # test for containing only one frequency or less than maxintermodorder
                # (sum( nvals .!== 0) == 1 || sum(abs,nvals) <= maxintermodorder)
                (count(!=(0), nvals) == 1 || sum(abs,nvals) <= maxintermodorder)

            push!(coords,i)
            push!(values,dot(w,nvals))
        else
            push!(dropcoords,i)
            push!(dropvalues,dot(w,nvals))
        end
    end

    return nothing
end


"""
    vectortodense(coords::Vector, values::Vector, Nharmonics)

Returns a dense array (vector or matrix) from a vector of coordinates, a
vector of values, and a tuple or single number giving the dimensions of the
array. The inputs can be generated with [`calcfrequencies`](@ref).

# Examples
```jldoctest
@variables w1
Nharmonics = (5);
coords = [CartesianIndex(2,), CartesianIndex(3,), CartesianIndex(4,), CartesianIndex(5,), CartesianIndex(6,)];
values = [w1, 2w1, 3w1, 4w1, 5w1];
JosephsonCircuits.vectortodense(coords,values,Nharmonics)

# output
6-element Vector{Num}:
   0
  w1
 2w1
 3w1
 4w1
 5w1
```
```jldoctest
@variables w1
Nharmonics = (5);
coords = [CartesianIndex(2,), CartesianIndex(4,), CartesianIndex(6,)];
values = [w1, 3w1, 5w1];
JosephsonCircuits.vectortodense(coords,values,Nharmonics)

# output
6-element Vector{Num}:
   0
  w1
   0
 3w1
   0
 5w1
```
```jldoctest
@variables w1, w2
Nharmonics = (3,2);
coords = [CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(1, 3), CartesianIndex(1, 4), CartesianIndex(1, 5), CartesianIndex(2, 5)]
values = [w1, 2w1, 3w1, w2, w1 + w2, 2w2, -2w2, -w2, w1 - w2]
JosephsonCircuits.vectortodense(coords,values,Nharmonics)

# output
4×5 Matrix{Num}:
   0       w2  2w2  -2w2      -w2
  w1  w1 + w2    0     0  w1 - w2
 2w1        0    0     0        0
 3w1        0    0     0        0
```
```jldoctest
@variables w1, w2
Nharmonics = (3,2);
coords = [CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(4, 3), CartesianIndex(1, 4), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(4, 4), CartesianIndex(1, 5), CartesianIndex(2, 5), CartesianIndex(3, 5), CartesianIndex(4, 5)];
values = [w1, 2w1, 3w1, w2, w1 + w2, w2 + 2w1, w2 + 3w1, 2w2, w1 + 2w2, 2w1 + 2w2, 2w2 + 3w1, -2w2, w1 - 2w2, 2w1 - 2w2, 3w1 - 2w2, -w2, w1 - w2, 2w1 - w2, 3w1 - w2];
JosephsonCircuits.vectortodense(coords,values,Nharmonics)

# output
4×5 Matrix{Num}:
   0        w2        2w2       -2w2       -w2
  w1   w1 + w2   w1 + 2w2   w1 - 2w2   w1 - w2
 2w1  w2 + 2w1  2w1 + 2w2  2w1 - 2w2  2w1 - w2
 3w1  w2 + 3w1  2w2 + 3w1  3w1 - 2w2  3w1 - w2
```
"""
function vectortodense(coords::Vector, values::Vector, Nharmonics)

    if length(Nharmonics) != length(eltype(coords))
        throw(DimensionMismatch("Dimensions of coords elements and Nharmonics must be consistent."))
    end

    if length(eltype(coords)) == 1
        s = zeros(eltype(values),Nharmonics[1]+1)
        for (i,coord) in enumerate(coords)
            s[coord] = values[i]
        end
    elseif length(eltype(coords)) == 2
        s = zeros(eltype(values),Nharmonics[1]+1,2*(Nharmonics[2]+1)-1)
        for (i,coord) in enumerate(coords)
            s[coord] = values[i]
        end
        
    else
        throw(DimensionMismatch("Not designed to visualize higher dimensional arrays"))
    end
    return s
end


"""
    conjsymdft(N::Tuple)

Calculate the conjugate symmetries in the multi-dimensional DFT (or FFT).

# Examples
```jldoctest
julia> JosephsonCircuits.conjsymdft(4)
Dict{CartesianIndex{1}, CartesianIndex{1}} with 1 entry:
  CartesianIndex(2,) => CartesianIndex(4,)

julia> JosephsonCircuits.conjsymdft((2,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 2 entries:
  CartesianIndex(2, 2) => CartesianIndex(2, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)

julia> JosephsonCircuits.conjsymdft((3,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 4 entries:
  CartesianIndex(2, 3) => CartesianIndex(3, 2)
  CartesianIndex(2, 1) => CartesianIndex(3, 1)
  CartesianIndex(2, 2) => CartesianIndex(3, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
```
"""
function conjsymdft(Nt)

    Nw=NTuple{length(Nt),Int}(val for (i,val) in enumerate(Nt))

    d = Dict{CartesianIndex{length(Nt)},CartesianIndex{length(Nt)}}()
    for k in CartesianIndices(Nw)
        # check that none of the indices are equal but that they are valid
        # indices. 
        if any(k.I .- 1 .!= mod.(Nt .- (k.I .- 1),Nt)) && all(mod.(Nt .- (k.I .- 1),Nt) .< Nt)
            tmp = sort([k.I,mod.(Nt .- (k.I .- 1),Nt) .+ 1])
            if !haskey(d,CartesianIndex(tmp[1]))
                d[CartesianIndex(tmp[1])] = CartesianIndex(tmp[2])
            end
        end
    end
    return d
end


"""
    conjsymdft(N::AbstractArray)

Calculate the conjugate symmetries in the multi-dimensional DFT (or FFT).

# Examples
```jldoctest
julia> JosephsonCircuits.conjsymdft([1 2;3 4])
Dict{CartesianIndex{2}, CartesianIndex{2}}()
```
"""
function conjsymdft(A::AbstractArray)
    return conjsymdft(size(A))
end

"""
    printdftsymmetries(N)

Print the conjugate symmetries in the multi-dimensional DFT of a real
signal from the dimensions of the signal in the time domain. Negative numbers
indicate that element is the complex conjugate of the corresponding positive
number. A zero indicates that element has no corresponding complex conjugate.

# Examples
```jldoctest
julia> JosephsonCircuits.printdftsymmetries(4)
4-element Vector{Int64}:
  0
  1
  0
 -1

julia> JosephsonCircuits.printdftsymmetries((3,3))
3×3 Matrix{Int64}:
  0   2  -2
  1   3   4
 -1  -4  -3

julia> JosephsonCircuits.printdftsymmetries((4,3))
4×3 Matrix{Int64}:
  0   2  -2
  1   3   5
  0   4  -4
 -1  -5  -3

julia> JosephsonCircuits.printdftsymmetries((3,4))
3×4 Matrix{Int64}:
  0   2   0  -2
  1   3   4   5
 -1  -5  -4  -3
```
"""
function printdftsymmetries(N)

    d=conjsymdft(N)

    z=zeros(Int,N)
    i = 1
    for (key,val) in sort(OrderedDict(d))
        z[key] = i
        z[val] = -i
        i+=1
    end
    return z
end

"""
    printdftsymmetries(A::AbstractArray)

# Examples
```jldoctest
julia> JosephsonCircuits.printdftsymmetries([1 2 3;4 5 6;7 8 9])
3×3 Matrix{Int64}:
  0   2  -2
  1   3   4
 -1  -4  -3
```
"""
function printdftsymmetries(A::AbstractArray)
    return printdftsymmetries(size(A))
end

"""
    conjsymrdft(Nt)

Calculate the conjugate symmetries in the multi-dimensional RDFT (DFT of
a real signal).

# Examples
```jldoctest
julia> JosephsonCircuits.conjsymrdft(4)
Dict{CartesianIndex{1}, CartesianIndex{1}}()

julia> JosephsonCircuits.conjsymdft(4)
Dict{CartesianIndex{1}, CartesianIndex{1}} with 1 entry:
  CartesianIndex(2,) => CartesianIndex(4,)

julia> JosephsonCircuits.conjsymrdft((2,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 2 entries:
  CartesianIndex(2, 2) => CartesianIndex(2, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)

julia> JosephsonCircuits.conjsymrdft((3,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 1 entry:
  CartesianIndex(1, 2) => CartesianIndex(1, 3)

julia> JosephsonCircuits.conjsymrdft((2,3,3))
Dict{CartesianIndex{3}, CartesianIndex{3}} with 8 entries:
  CartesianIndex(1, 2, 1) => CartesianIndex(1, 3, 1)
  CartesianIndex(2, 2, 3) => CartesianIndex(2, 3, 2)
  CartesianIndex(2, 2, 1) => CartesianIndex(2, 3, 1)
  CartesianIndex(1, 2, 2) => CartesianIndex(1, 3, 3)
  CartesianIndex(1, 1, 2) => CartesianIndex(1, 1, 3)
  CartesianIndex(2, 2, 2) => CartesianIndex(2, 3, 3)
  CartesianIndex(2, 1, 2) => CartesianIndex(2, 1, 3)
  CartesianIndex(1, 2, 3) => CartesianIndex(1, 3, 2)
```
"""
function conjsymrdft(Nt)

    Nw=NTuple{length(Nt),Int}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))

    d = Dict{CartesianIndex{length(Nt)},CartesianIndex{length(Nt)}}()

    for k in CartesianIndices(Nw)
        # check that none of the indices are equal but that they are valid
        # indices. 
        if any(k.I .- 1 .!= mod.(Nt .- (k.I .- 1),Nt)) && all(mod.(Nt .- (k.I .- 1),Nt) .< Nw)
            # sort the terms to consistently decide which to call the conjugate
            tmp1 = k.I
            tmp2 = mod.(Nt .- (k.I .- 1),Nt) .+ 1
            if tmp1 < tmp2
                if !haskey(d,CartesianIndex(tmp1))
                    d[CartesianIndex(tmp1)] = CartesianIndex(tmp2)
                end
            else
                if !haskey(d,CartesianIndex(tmp2))
                    d[CartesianIndex(tmp2)] = CartesianIndex(tmp1)
                end
            end
        end
    end
    return d
end

"""
    conjsymrdft(A::AbstractArray)

Calculate the conjugate symmetries in the multi-dimensional RDFT (DFT of
a real signal).

# Examples
```jldoctest
julia> JosephsonCircuits.conjsymrdft([1 2;3 4])
Dict{CartesianIndex{2}, CartesianIndex{2}}()
```
"""
function conjsymrdft(A::AbstractArray)
    return conjsymdft(size(A))
end


"""
    printrdftsymmetries(N::Tuple)

Print the conjugate symmetries in the multi-dimensional RDFT (DFT of a real
signal) from the dimensions of the signal in the time domain. Negative numbers
indicate that element is the complex conjugate of the corresponding positive
number. A zero indicates that element has no corresponding complex conjugate.

# Examples
```jldoctest
julia> JosephsonCircuits.printrdftsymmetries(4)
3-element Vector{Int64}:
 0
 0
 0

julia> JosephsonCircuits.printrdftsymmetries((3,3))
2×3 Matrix{Int64}:
 0  1  -1
 0  0   0

julia> JosephsonCircuits.printrdftsymmetries((4,3))
3×3 Matrix{Int64}:
 0  1  -1
 0  0   0
 0  2  -2

julia> JosephsonCircuits.printrdftsymmetries((3,4))
2×4 Matrix{Int64}:
 0  1  0  -1
 0  0  0   0
```
"""
function printrdftsymmetries(Nt)

    Nw=NTuple{length(Nt),Int}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))

    d=conjsymrdft(Nt)

    z=zeros(Int,Nw)
    i = 1
    for (key,val) in sort(OrderedDict(d))
        z[key] = i
        z[val] = -i
        i+=1
    end
    return z
end

"""
    printrdftsymmetries(A::AbstractArray)

Print the conjugate symmetries in the multi-dimensional RDFT (DFT of a real
signal) from the dimensions of the signal in the time domain. Negative numbers
indicate that element is the complex conjugate of the corresponding positive
number. A zero indicates that element has no corresponding complex conjugate.

# Examples
```jldoctest
julia> JosephsonCircuits.printrdftsymmetries([1 2 3;4 5 6;7 8 9])
2×3 Matrix{Int64}:
 0  1  -1
 0  0   0
```
"""
function printrdftsymmetries(A::AbstractArray)
    return printrdftsymmetries(size(A))
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
    calcfrequencies2(Nt,coords,values)

Return the vector of values with any conjugate symmetric terms removed. Requires
the dimensions of the time domain signal, the vector of coordinates, and the
vector of values.

# Examples
```jldoctest
@variables w1, w2
Nt = (7, 5)
coords = [ CartesianIndex(2, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(3, 2), CartesianIndex(2, 3), CartesianIndex(4, 3), CartesianIndex(2, 4), CartesianIndex(4, 4), CartesianIndex(1, 5), CartesianIndex(3, 5)]
values = [w1, 3w1, w2, w2 + 2w1, w1 + 2w2, 2w2 + 3w1, w1 - 2w2, 3w1 - 2w2, -w2, 2w1 - w2]
JosephsonCircuits.calcfrequencies2(Nt,coords,values)

# output
9-element Vector{Num}:
        w1
       3w1
        w2
  w2 + 2w1
  w1 + 2w2
 2w2 + 3w1
  w1 - 2w2
 3w1 - 2w2
  2w1 - w2
```
"""
function calcfrequencies2(Nt, coords, values)

    indexdict = conjsymrdft(Nt);
    Nw=NTuple{length(Nt),Int}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))
    reverseindexdict = Dict{CartesianIndex{length(Nt)},CartesianIndex{length(Nt)}}()
    for (key,val) in indexdict
        reverseindexdict[val] = key
    end

    values2 = Array{eltype(values),1}(undef,0)

    carttoint = calcindexdict(Nw)
    inttocart = CartesianIndices(Nw)

    for (i,coord) in enumerate(coords)
        if !haskey(reverseindexdict,coord)
            push!(values2,values[i])
        end
    end

    return values2
end

"""
    calcphiindices(Nt, dropdict)

Return the indices which map the elements of the frequency domain vector
elements to the corresponding elements of the frequency domain array. Also
return the indices `conjsourceindices` whose data should be copied from the
vector to `conjtargetindices` in the array then complex conjugated

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
@variables w1, w2
Nt = (7, 5);
dropdict = Dict{CartesianIndex{2}, Num}(CartesianIndex(4, 2) => w2 + 3w1, CartesianIndex(1, 3) => 2w2, CartesianIndex(1, 1) => 0, CartesianIndex(2, 5) => w1 - w2, CartesianIndex(3, 3) => 2w1 + 2w2, CartesianIndex(3, 1) => 2w1, CartesianIndex(2, 2) => w1 + w2, CartesianIndex(1, 4) => -2w2, CartesianIndex(3, 4) => 2w1 - 2w2, CartesianIndex(4, 5) => 3w1 - w2);
indexmap,conjsourceindices,conjtargetindices = JosephsonCircuits.calcphiindices(Nt,dropdict);
println(indexmap)
println(conjsourceindices)
println(conjtargetindices)

# output
[2, 4, 5, 7, 10, 12, 14, 16, 19]
[5]
[17]
```
"""
function calcphiindices(Nt, dropdict)

    indices = conjsymrdft(Nt);
    Nw=NTuple{length(Nt),Int}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))
    conjindices = Dict{CartesianIndex{length(Nt)},CartesianIndex{length(Nt)}}()
    for (key,val) in indices
        conjindices[val] = key
    end

    indexmap = Array{Int,1}(undef,0)

    # the index of the element of the N dimensional array in the frequency domain
    # that i should copy to conjtargetindices and take the complex conjugate of. 
    conjsourceindices = Array{Int,1}(undef,0)

    # the index of the element of the N dimensional array in the frequency domain
    # that i should take the complex conjugate of
    conjtargetindices = Array{Int,1}(undef,0)

    carttoint = calcindexdict(Nw)
    inttocart = CartesianIndices(Nw)

    # generate the index maps to convert between the vector
    # and matrix
    for (j,index) in enumerate(CartesianIndices(Nw))
        # don't include the index in the maps if it shows up
        # in dropindices
        if !haskey(dropdict,index)
            if haskey(conjindices,index)
                push!(conjtargetindices,carttoint[index])
                push!(conjsourceindices,carttoint[conjindices[index]])
            else
                push!(indexmap,j)
            end
        end
    end

    return indexmap, conjsourceindices, conjtargetindices
end

"""
    calcphiindices2(Nt, dropdict)

Return the indices which map the elements of the frequency domain vector
elements to the corresponding elements of the frequency domain array. Also
return the indices `conjsourceindices` whose data should be copied from the
vector to `conjtargetindices` in the array then complex conjugated

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
freq = JosephsonCircuits.Frequencies{2}((5, 7), (8, 7), CartesianIndex{2}[CartesianIndex(2, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(3, 2), CartesianIndex(2, 3), CartesianIndex(1, 4), CartesianIndex(2, 6), CartesianIndex(3, 7)], [(1, 0), (3, 0), (0, 1), (2, 1), (1, 2), (0, 3), (1, -2), (2, -1)])
conjsymdict = Dict{CartesianIndex{2}, CartesianIndex{2}}(CartesianIndex(5, 4) => CartesianIndex(5, 5), CartesianIndex(1, 3) => CartesianIndex(1, 6), CartesianIndex(5, 2) => CartesianIndex(5, 7), CartesianIndex(1, 4) => CartesianIndex(1, 5), CartesianIndex(1, 2) => CartesianIndex(1, 7), CartesianIndex(5, 3) => CartesianIndex(5, 6))
JosephsonCircuits.calcphiindices2(freq, conjsymdict)

# output
([2, 4, 6, 8, 12, 16, 27, 33], [6, 16], [31, 21])
```
```jldoctest
freq = JosephsonCircuits.calcfreqsrdft((4,3));
truncfreq = JosephsonCircuits.truncfreqsrdft(freq;dc=false,odd=true,even=false,maxintermodorder=3)
noconjtruncfreq = JosephsonCircuits.removeconjfreqsrdft(truncfreq)
conjsymdict = JosephsonCircuits.conjsymrdft(noconjtruncfreq.Nt)
JosephsonCircuits.calcphiindices2(noconjtruncfreq,conjsymdict)

# output
([2, 4, 6, 8, 12, 16, 27, 33], [6, 16], [31, 21])
```
"""
function calcphiindices2(frequencies::JosephsonCircuits.Frequencies{N},
        conjsymdict::Dict{CartesianIndex{N},CartesianIndex{N}}) where N

    modes = frequencies.modes
    coords = frequencies.coords
    Nw = frequencies.Nw
    Nt = frequencies.Nt

    coordsdict = Dict{CartesianIndex{N},Int}()
    for (i,coord) in enumerate(coords)
        coordsdict[coord] = i
    end

    # empty vector to hold the map between indices in the vector and the
    # matrix
    indexmap = Vector{Int}(undef,0)

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
        push!(indexmap,carttoint[coord])
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

    Nvector = length(phivector)÷ Nbranches
    Nmatrix = prod(size(phimatrix)[1:end-1])

    # fill the matrix with zeros
    fill!(phimatrix,0)

    for i in 1:Nbranches
        for j in 1:length(indexmap)
            phimatrix[indexmap[j]+(i-1)*Nmatrix] = phivector[j+(i-1)*Nvector]
        end
    end

    for i in 1:Nbranches
        for j in 1:length(conjtargetindices)
            phimatrix[conjtargetindices[j]+ (i-1)*Nmatrix] = conj(phimatrix[conjsourceindices[j]+ (i-1)*Nmatrix])
        end
    end
    return nothing
end

"""
    phimatrixtovector!(phivector::Vector,phimatrix::Array,
        indexmap::Vector{Int},conjsourceindices::Vector{Int},
        conjtargetindices::Vector{Int},Nbranches::Int)

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

    Nvector = length(phivector)÷ Nbranches
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
julia> JosephsonCircuits.applynl([[0, 0.2+0.0im, 0, 0];;],(x)->cos(x))
4×1 Matrix{ComplexF64}:
   0.9603980498951228 + 0.0im
                  0.0 + 0.0im
 -0.01966852794611884 + 0.0im
                  0.0 + 0.0im

julia> JosephsonCircuits.applynl([[0, 0.2+0.0im];;],(x)->cos(x))
2×1 Matrix{ComplexF64}:
   0.9603980498951228 + 0.0im
 -0.01966852794611884 + 0.0im

julia> JosephsonCircuits.applynl([0.0 + 0.0im 0.45 + 0.0im 0.45 + 0.0im; 0.55 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;;;],(x)->sin(x))
3×3×1 Array{ComplexF64, 3}:
[:, :, 1] =
 -0.0209812+0.0im   0.295151+0.0im   0.295151+0.0im
   0.359826+0.0im  -0.041417+0.0im  -0.041417+0.0im
 0.00788681+0.0im  -0.110947+0.0im  -0.110947+0.0im
```
"""
function applynl(fd::Array{Complex{Float64}}, f::Function)

    td, irfftplan, rfftplan = plan_applynl(fd)
    fdcopy = copy(fd)

    applynl!(fdcopy, td, f, irfftplan, rfftplan)

    return fdcopy
end

"""
    plan_applynl(fd)

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
        rfftplan) where T

Apply the nonlinear function f to the frequency domain data by transforming
to the time domain, applying the function, then transforming back to the
frequency domain, overwriting the contents of fd and td in the process. We
use plans for the forward and reverse RFFT prepared by [`plan_applynl`](@ref).

# Examples
```jldoctest
fd=ones(Complex{Float64},3,2)
td, irfftplan, rfftplan = JosephsonCircuits.plan_applynl(fd)
JosephsonCircuits.applynl!(fd, td, (x)->cos(x), irfftplan, rfftplan)
fd

# output
3×2 Matrix{ComplexF64}:
  0.586589+0.0im   0.586589+0.0im
 -0.413411+0.0im  -0.413411+0.0im
 -0.413411+0.0im  -0.413411+0.0im
```
"""
function applynl!(fd::Array{Complex{T}}, td::Array{T}, f::Function, irfftplan,
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
    hbmatind(Nharmonics::Tuple; maxintermodorder::Number = Inf,
        dc::Bool = true, even::Bool = true, odd::Bool = true)

Returns a matrix describing which indices of the frequency domain matrix
(from the RFFT) to pull out and use in the harmonic balance matrix. A negative
index means we take the complex conjugate of that element. A zero index means
that term is not present, so skip it. The harmonic balance matrix describes
the coupling between different frequency modes.

# Examples
```jldoctest
julia> JosephsonCircuits.hbmatind((5,);maxintermodorder=2,dc=false, even=false, odd=true)
3×3 Matrix{Int64}:
 1  -3  -5
 3   1  -3
 5   3   1

julia> JosephsonCircuits.hbmatind((3,);maxintermodorder=2,dc=true, even=true, odd=true)
4×4 Matrix{Int64}:
 1  -2  -3  -4
 2   1  -2  -3
 3   2   1  -2
 4   3   2   1

julia> JosephsonCircuits.hbmatind((2,2);maxintermodorder=2,dc=true, even=true, odd=true)
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
function hbmatind(Nharmonics::Tuple; maxintermodorder::Number = Inf,
    dc::Bool = true, even::Bool = true, odd::Bool = true)


    # this code is ugly but the idea is that we need all of the possible values
    # to make satisfy transitions between elements and to find the right indices
    Nw,coords,valuesall,dropcoords,dropvalues = JosephsonCircuits.calcfrequencies(
        Nharmonics,
        maxintermodorder=Inf,
        dc=true,
        even=true,
        odd=true,
    )

    # now we calculate the reduced set of frequency indices
    Nw,coords,values,dropcoords,dropvalues = JosephsonCircuits.calcfrequencies(
        Nharmonics,
        maxintermodorder=maxintermodorder,
        dc=dc,
        even=even,
        odd=odd,
    )
    Nt=NTuple{length(Nw),Int}(ifelse(i == 1, 2*val-1, val) for (i,val) in enumerate(Nw))
    values2 = JosephsonCircuits.calcfrequencies2(Nt,coords,values)

    # assign the frequencies
    wmodes = values2[:]
    Nmodes = length(wmodes)

    # println(wmodes)
    # this is calculating the frequency domain input output relations
    Amatrixkeys = Matrix{eltype(wmodes)}(undef,length(wmodes),length(wmodes))
    for i in 1:length(wmodes)
        for j in 1:length(wmodes)
            key = Tuple(wmodes[i][k]-wmodes[j][k] for k in 1:length(wmodes[i]))
            Amatrixkeys[i,j] = key
        end
    end

    # println(Amatrixkeys)

    # now i need to find the keys that are in the rfft matrix and their locations
    valuesdict = Dict{eltype(valuesall),Int}()
    for (i,v) in enumerate(valuesall)
        valuesdict[v] = i
    end

    Amatrixindices = zeros(Int,length(wmodes),length(wmodes))
    for i in 1:length(Amatrixkeys)
        key = Amatrixkeys[i]
        conjkey = Tuple(-k for k in key)
        if haskey(valuesdict,key)
            Amatrixindices[i] = valuesdict[key]
        elseif haskey(valuesdict,conjkey)
            Amatrixindices[i] = -valuesdict[conjkey]
        end
    end

    return Amatrixindices
end


"""
    hbmatind2(frequencies::JosephsonCircuits.Frequencies{N},
        truncfrequencies::JosephsonCircuits.Frequencies{N}) where N

Returns a matrix describing which indices of the frequency domain matrix
(from the RFFT) to pull out and use in the harmonic balance matrix. A negative
index means we take the complex conjugate of that element. A zero index means
that term is not present, so skip it. The harmonic balance matrix describes
the coupling between different frequency modes.

# Examples
```jldoctest
julia> freq = JosephsonCircuits.calcfreqsrdft((5,));JosephsonCircuits.hbmatind2(freq,JosephsonCircuits.removeconjfreqsrdft(JosephsonCircuits.truncfreqsrdft(freq;dc=false,odd=true,even=false,maxintermodorder=2)))
3×3 Matrix{Int64}:
 1  -3  -5
 3   1  -3
 5   3   1

julia> freq = JosephsonCircuits.calcfreqsrdft((3,));JosephsonCircuits.hbmatind2(freq,JosephsonCircuits.removeconjfreqsrdft(JosephsonCircuits.truncfreqsrdft(freq;dc=true,odd=true,even=true,maxintermodorder=2)))
4×4 Matrix{Int64}:
 1  -2  -3  -4
 2   1  -2  -3
 3   2   1  -2
 4   3   2   1

julia> freq = JosephsonCircuits.calcfreqsrdft((2,2));JosephsonCircuits.hbmatind2(freq,JosephsonCircuits.removeconjfreqsrdft(JosephsonCircuits.truncfreqsrdft(freq;dc=true,odd=true,even=true,maxintermodorder=2)))
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
function hbmatind2(frequencies::JosephsonCircuits.Frequencies{N},
    truncfrequencies::JosephsonCircuits.Frequencies{N}) where N

    Nw = frequencies.Nw
    Nt = frequencies.Nt
    modes = frequencies.modes
    truncmodes = truncfrequencies.modes

    if length(truncmodes) > length(modes)
        error("Number of truncated frequencies should be not greater than number of frequencies.")
    end

    # this is calculating the frequency domain input output relations
    Amatrixkeys = Matrix{eltype(truncmodes)}(undef,length(truncmodes),length(truncmodes))
    for i in 1:length(truncmodes)
        for j in 1:length(truncmodes)
            key = NTuple{N,Int}(truncmodes[i][k]-truncmodes[j][k] for k in 1:length(truncmodes[i]))
            Amatrixkeys[i,j] = key
        end
    end

    # now i need to find the keys that are in the rfft matrix and their locations
    modesdict = Dict{eltype(modes),Int}()
    for (i,mode) in enumerate(modes)
        modesdict[mode] = i
    end

    Amatrixindices = zeros(Int,length(truncmodes),length(truncmodes))
    for i in 1:length(Amatrixkeys)
        key = Amatrixkeys[i]
        conjkey = Tuple(-k for k in key)
        if haskey(modesdict,key)
            Amatrixindices[i] = modesdict[key]
        elseif haskey(modesdict,conjkey)
            Amatrixindices[i] = -modesdict[conjkey]
        end
    end

    return Amatrixindices
end