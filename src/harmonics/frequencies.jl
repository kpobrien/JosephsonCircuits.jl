
"""
    Frequencies(Nharmonics::NTuple{N, Int}, Nw::NTuple{N,Int}, Nt::NTuple{N,Int},
        coords::Vector{CartesianIndex{N}}, modes::Vector{NTuple{N,Int}})

A simple structure to hold time and frequency domain information for the
signals. See also [`calcfreqsrdft`](@ref) and [`calcfreqsdft`](@ref).

# Fields
- `Nharmonics::NTuple{N, Int}`: The number of harmonics of each frequency
    the grid samples, which sets `Nw` and `Nt`; the modes kept may be fewer.
- `Nw::NTuple{N,Int}`: The dimensions of the frequency domain signal for a
    single node.
- `Nt::NTuple{N,Int}`: The dimensions of the time domain signal for a single
    node.
- `coords::Vector{CartesianIndex{N}}`: The coordinates of each mixing products.
- `modes::Vector{NTuple{N,Int}}`: The mode indices of each mixing product, eg.
     (0,0), (1,0), (2,1).
"""
struct Frequencies{N}
    Nharmonics::NTuple{N, Int}
    Nw::NTuple{N, Int}
    Nt::NTuple{N, Int}
    coords::Vector{CartesianIndex{N}}
    modes::Vector{NTuple{N, Int}}
end

"""
    ModeDifferences{N} <: AbstractMatrix{NTuple{N,Int}}

The mode difference `modes[i] .- modes[j]` of every pair of modes, the
harmonic offset of each mode coupling, computed when an entry is read
rather than stored.
"""
struct ModeDifferences{N} <: AbstractMatrix{NTuple{N,Int}}
    modes::Vector{NTuple{N,Int}}
end

Base.size(D::ModeDifferences) = (length(D.modes), length(D.modes))
Base.@propagate_inbounds Base.getindex(D::ModeDifferences, i::Int,
    j::Int) = D.modes[i] .- D.modes[j]

"""
    FourierIndices(vectomatmap::Vector{Int}, conjsourceindices::Vector{Int},
        conjtargetindices::Vector{Int}, hbmatmodes::ModeDifferences{N},
        hbconjmatindices::Matrix{Int})

A simple structure to hold time and frequency domain information for the
signals, particularly the indices for converting between the node flux vectors
and matrices. The `hbmatmodes` matrix holds the differences of the modes,
which describe the coupling between the modes (the derivative of the residual
with respect to the node fluxes), while the `hbconjmatindices` matrix is built
from the sums of the modes, aliased back onto the sampled grid, and describes
the coupling between the modes and the complex conjugates of the modes (the
derivative of the residual with respect to the complex conjugates of the node
fluxes).

The differences are computed entry by entry when read
([`ModeDifferences`](@ref)): a solve places them on the sampled grid with
[`hbmatindices`](@ref), with aliasing for the exact Jacobian and without it
for the holomorphic one, and selects the preconditioner's bands by them.
See also [`fourierindices`](@ref).
"""
struct FourierIndices{N}
    vectomatmap::Vector{Int}
    conjsourceindices::Vector{Int}
    conjtargetindices::Vector{Int}
    hbmatmodes::ModeDifferences{N}
    hbconjmatindices::Matrix{Int}
end

"""
    fourierindices(freq::Frequencies)

Generate the indices used in the RDFT or DFT and inverse RDFT or DFT and
converting between a node flux vector for solving system and the matrices for
the Fourier analysis. See also [`FourierIndices`](@ref), [`Frequencies`](@ref),
[`calcfreqsrdft`](@ref) and [`calcfreqsdft`](@ref).

"""
function fourierindices(freq::Frequencies)

    freqindexmap, conjsourceindices, conjtargetindices =
        calcphiindices(freq, conjsym(freq))
    # the differences of the modes, computed when an entry is read
    Amatrixmodes = ModeDifferences(freq.modes)
    Amatrixconjindices = hbconjmatind(freq)

    return FourierIndices(
        freqindexmap,
        conjsourceindices,
        conjtargetindices,
        Amatrixmodes,
        Amatrixconjindices,
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
"""
function calcfreqsrdft(Nharmonics::NTuple{N,Int}) where N

    # every dimension but the first holds both signs of frequency; the
    # first dimension of a real transform holds only the positive ones
    Nw=NTuple{N,Int}(ifelse(i == 1, val+1, 2*val+1) for (i,val) in enumerate(Nharmonics))
    
    # the number of time points in each dimension; an odd number in the
    # first dimension so that the highest mode is not the Nyquist mode,
    # whose coefficient a real transform forces to be real
    Nt =  NTuple{N,Int}(ifelse(i == 1, 2*Nw[1]-1, val) for (i,val) in enumerate(Nw))

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
"""
function removeconjfreqs(frequencies::Frequencies)
    conjsymdict = conjsym(frequencies)
    return removefreqs(frequencies,collect(values(conjsymdict)))
end

"""
    removefreqs(frequencies::Frequencies{N},
        removemodes::AbstractVector{NTuple{N,Int}})
    removefreqs(frequencies::Frequencies{N},
        removecoords::AbstractVector{CartesianIndex{N}})

Return a new [`Frequencies`](@ref) without the modes `removemodes`, or
without the modes at the coordinates `removecoords`, keeping the others in
their order.
"""
removefreqs(frequencies::Frequencies{N},
    removemodes::AbstractVector{NTuple{N,Int}}) where N =
    keepunlisted(frequencies, frequencies.modes, Set(removemodes))
removefreqs(frequencies::Frequencies{N},
    removecoords::AbstractVector{CartesianIndex{N}}) where N =
    keepunlisted(frequencies, frequencies.coords, Set(removecoords))

# the frequencies whose key, the mode or the coordinate, is not in `removed`
function keepunlisted(frequencies::Frequencies, keys, removed)
    keep = [!(k in removed) for k in keys]
    return Frequencies(frequencies.Nharmonics, frequencies.Nw, frequencies.Nt,
        frequencies.coords[keep], frequencies.modes[keep])
end

"""
    truncfreqs(frequencies::Frequencies; maxharmonics = frequencies.Nharmonics,
        maxintermodorder = Inf, dc = true, odd = true, even = true,
        w = nothing, frequencywindow = (0, Inf))

Return a new [`Frequencies`](@ref) with the coordinates and modes truncated
to those which satisfy the criteria: the zero frequency mode when `dc`,
modes of odd or even total harmonic order when `odd` or `even`, modes
which are a harmonic of a single tone or whose absolute harmonic indices
sum to at most `maxintermodorder`, and modes whose absolute harmonic
index in each tone is at most `maxharmonics` for that tone.

With the tone frequencies `w` given and a window other than the default
`(0, Inf)`, a mode is also required to lie in the `frequencywindow`:
`wmin <= abs(dot(w, mode)) <= wmax`, in the units of `w`. The zero
frequency mode is governed by `dc` alone. This is a truncation by frequency
rather than by order: incommensurate tones scatter combination frequencies
of high order arbitrarily close to zero, where a floating circuit's linear
response is enormous although nothing excites those modes, and near the
junction plasma frequency at the other end, at the edge of the sampled grid
where their products alias. Such modes carry almost no flux at the
operating point, and their nearly singular blocks are what a block diagonal
preconditioner, and every preconditioner built from it, inverts badly: the
ceiling can decide whether such a preconditioner converges a solve, while
the retained modes barely change.

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
JosephsonCircuits.Frequencies{2}((3, 3), (4, 7), (7, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(1, 3), CartesianIndex(1, 4), CartesianIndex(1, 5), CartesianIndex(1, 6), CartesianIndex(1, 7), CartesianIndex(2, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (0, 2), (0, 3), (0, -3), (0, -2), (0, -1), (1, -1)])

julia> JosephsonCircuits.truncfreqs(JosephsonCircuits.calcfreqsrdft((3,3));dc=false,even=false,maxharmonics=(2,2),maxintermodorder=3).modes
7-element Vector{Tuple{Int64, Int64}}:
 (1, 0)
 (0, 1)
 (2, 1)
 (1, 2)
 (1, -2)
 (0, -1)
 (2, -1)
```
"""
function  truncfreqs(frequencies::Frequencies{N};
    maxharmonics::NTuple{N,Int} = frequencies.Nharmonics,
    maxintermodorder = Inf,
    dc::Bool = true, odd::Bool = true, even::Bool = true,
    w = nothing, frequencywindow = (0, Inf)) where N

    coords = frequencies.coords
    modes = frequencies.modes

    length(frequencywindow) == 2 || throw(ArgumentError(
        "`frequencywindow` is a tuple `(wmin, wmax)`."))
    wmin, wmax = frequencywindow
    0 <= wmin <= wmax || throw(ArgumentError(
        lazy"`frequencywindow` = $(frequencywindow) must satisfy 0 <= wmin <= wmax."))
    windowed = !isnothing(w) && (wmin > 0 || isfinite(wmax))
    if windowed
        length(w) == N || throw(DimensionMismatch(
            lazy"`w` has $(length(w)) tones but the frequencies have $(N)."))
    end
    inwindow(nvals) = !windowed || all(==(0), nvals) ||
        (wmin <= abs(sum(w[k]*nvals[k] for k in 1:N)) <= wmax)

    keepmodes = Vector{eltype(modes)}(undef,0)
    sizehint!(keepmodes,length(modes))

    keepcoords = Vector{eltype(coords)}(undef,0)
    sizehint!(keepcoords,length(modes))

    for (i,nvals) in enumerate(modes)

        # a mode is kept when it matches the dc, even or odd criterion, and
        # is either a harmonic of a single tone or within the intermodulation
        # order, and is within the per tone harmonic bound
        if (
                # test for DC
                (dc && all(==(0),nvals)) ||
                # test for even (and not DC)
                (even && mod(sum(abs,nvals),2) == 0 && sum(abs,nvals) > 0) ||
                # test for odd
                (odd && mod(sum(abs,nvals),2) == 1)
            ) && # a harmonic of one tone, or within the intermodulation order
                (count(!=(0), nvals) == 1 || sum(abs,nvals) <= maxintermodorder) &&
                # and less than the maxharmonics
                all(map(<=,abs.(nvals),maxharmonics)) &&
                # and, with the tone frequencies given, inside the window
                inwindow(nvals)

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
    conjsym(Nw::NTuple{N, Int}, Nt::NTuple{N, Int})

Calculate the conjugate symmetries in the multi-dimensional frequency domain
data.
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
    calcphiindices(frequencies::Frequencies{N},
        conjsymdict::Dict{CartesianIndex{N},CartesianIndex{N}})

Return the indices which map the elements of the frequency domain vector
to the corresponding elements of the frequency domain array. Also
return the indices `conjsourceindices` whose data should be copied from the
vector to `conjtargetindices` in the array then complex conjugated.

# Arguments
- `frequencies`: the retained frequencies, whose coordinates are the
    positions of the vector's elements in the frequency domain array.
- `conjsymdict`: the conjugate symmetric pairs of coordinates of the array
    (see [`conjsym`](@ref)); a retained coordinate with a pair has its
    conjugate written at the paired coordinate.

# Returns
- `indexmap`: the indices which map the elements of the frequency domain
    vector elements to the corresponding elements of the frequency domain array
- `conjsourceindices`: data should be copied from here
- `conjtargetindices`: data should be copied to here and conjugated

# Examples
```jldoctest
freq = JosephsonCircuits.Frequencies{2}((4, 3), (5, 7), (9, 7), CartesianIndex{2}[CartesianIndex(2, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(3, 2), CartesianIndex(2, 3), CartesianIndex(1, 4), CartesianIndex(2, 6), CartesianIndex(3, 7)], [(1, 0), (3, 0), (0, 1), (2, 1), (1, 2), (0, 3), (1, -2), (2, -1)])
conjsymdict = Dict{CartesianIndex{2}, CartesianIndex{2}}(CartesianIndex(1, 2) => CartesianIndex(1, 7), CartesianIndex(1, 3) => CartesianIndex(1, 6), CartesianIndex(1, 4) => CartesianIndex(1, 5))
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

    coords = frequencies.coords
    Nw = frequencies.Nw

    # the position in the frequency domain array of each mode of the vector
    indexmap = Vector{Int}(undef,length(coords))

    # the positions in the frequency domain array which are copied,
    # conjugated, to `conjtargetindices`
    conjsourceindices = Array{Int}(undef,0)

    # the positions which receive the conjugates
    conjtargetindices = Vector{Int}(undef,0)

    # the position of each coordinate in the frequency domain array
    carttoint = LinearIndices(Nw)

    # the vector to matrix map, in the order of `coords`
    for (i,coord) in enumerate(coords)
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
        throw(DimensionMismatch(lazy"Unexpected length for phivector"))
    end

    if length(phivector) == 0
        Nvector = 0
    else
        Nvector = length(phivector)÷ Nbranches
    end

    Nmatrix = prod(size(phimatrix)[1:end-1])

    # fill the matrix with zeros
    fill!(phimatrix,0)

    # copy the mixing products from the vector to the matrix
    for i in 1:Nbranches
        for j in 1:length(indexmap)
            phimatrix[indexmap[j]+(i-1)*Nmatrix] = phivector[j+(i-1)*Nvector]
        end
    end

    # the conjugate modes, copied and conjugated within the matrix
    for i in 1:Nbranches
        for j in 1:length(conjtargetindices)
            phimatrix[conjtargetindices[j]+ (i-1)*Nmatrix] = conj(phimatrix[conjsourceindices[j]+ (i-1)*Nmatrix])
        end
    end
    return nothing
end

"""
    plan_applynl(fd::AbstractArray{Complex{T}}, backend::Backend = CPU())

Creates an empty time domain data array and the inverse and forward plans
for the RFFT of an array of frequency domain data. See also [`applynl!`](@ref).
A system with no junctions has nothing to transform: for an empty `fd` the
two plans are `nothing`, and applying them is the identity.

"""
function plan_applynl(fd::AbstractArray{Complex{T}},
    backend::Backend = CPU()) where T

    sizefd = size(fd)
    stepsperperiod = 2*sizefd[1]-1

    # the time domain array, on the same device as the frequency domain
    # array
    dims = (stepsperperiod, sizefd[2:end]...)
    td = similar(fd, T, dims)

    # A system with no Josephson junctions has an empty batch dimension and
    # nothing to transform. FFTW plans that but cuFFT rejects it, so the
    # plan is `nothing` and applying it is the identity (see `applyifft!`).
    irfftplan, rfftplan = if isempty(fd)
        (nothing, nothing)
    else
        fftplans(fd, td, stepsperperiod, backend)
    end

    return td, irfftplan, rfftplan
end

"""
    fftplans(fd::AbstractArray{Complex{T}}, td::AbstractArray{T},
        stepsperperiod::Int, backend::Backend)

Create the inverse real transform plan from the frequency domain array `fd`
to the time domain array `td`, and the forward plan back, on the given
KernelAbstractions backend. The transform runs over all but the last
dimension, the last being the Josephson junction index. The inverse plan is
the unnormalized backward transform, which gives the time domain samples
in the convention of [`applyifft!`](@ref) with no scaling pass; the
forward plan is unnormalized too, and [`applyfft!`](@ref) scales it.

The `CPU()` method uses FFTW. A device backend supplies its own method, which
is the only thing the residual and the matrix-free products need that the
core package cannot provide without taking on the device dependency: every
other step of those is either a plain array operation or a kernel of
[`NonlinearTermPlan`](@ref). Load the package extension for the device (for
CUDA, `using CUDA`) to get its method.
"""
function fftplans(fd::AbstractArray{Complex{T}}, td::AbstractArray{T},
    stepsperperiod::Int, backend::CPU) where T
    dims = 1:length(size(fd))-1
    irfftplan = FFTW.plan_brfft(fd, stepsperperiod, dims;
        flags = FFTW.ESTIMATE, timelimit = Inf)
    rfftplan = FFTW.plan_rfft(td, dims; flags = FFTW.ESTIMATE, timelimit = Inf)
    return irfftplan, rfftplan
end

function fftplans(fd::AbstractArray, td::AbstractArray, stepsperperiod::Int,
    backend::Backend)
    throw(ArgumentError(lazy"no real transform plans are defined for the backend $(backend). Load the package extension for the device (for CUDA, `using CUDA`) or use the CPU() backend."))
end

"""
    applynl!(fd::AbstractArray{Complex{T}}, td::AbstractArray{T}, f,
        irfftplan, rfftplan)
    applynl!(fd, td, f, ::Nothing, ::Nothing)

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
  0.856732+0.0im   0.856732+0.0im
 -0.143268+0.0im  -0.143268+0.0im
 -0.143268+0.0im  -0.143268+0.0im
```
"""
function applynl!(fd::AbstractArray{Complex{T}}, td::AbstractArray{T}, f, irfftplan,
    rfftplan) where T

    applyifft!(td, fd, irfftplan)
    # broadcasting keeps this device generic
    td .= f.(td)
    applyfft!(fd, td, rfftplan)
    return nothing
end

# with no Josephson junctions there is nothing to transform and no plan; see
# [`plan_applynl`](@ref)
applynl!(fd::AbstractArray{Complex{T}}, ::AbstractArray{T}, f, ::Nothing,
    ::Nothing) where T = fd

"""
    applyrelationnl!(fd, td, work, relations, coefficients, trig, irfftplan,
        rfftplan)

[`applynl!`](@ref) with the per junction relation of `relations` in place
of a single function: `coefficients` are the polynomial coefficients to
evaluate, one row per polynomial junction of the table, and `trig` is what
its sinusoidal junctions take instead. `work` is a time domain array the size of
`td`, which the Horner loop reads while it writes `td`.
"""
function applyrelationnl!(fd::AbstractArray{Complex{T}}, td::AbstractArray{T},
        work::AbstractArray{T}, relations, coefficients, trig::F, irfftplan,
        rfftplan) where {T,F}
    applyifft!(work, fd, irfftplan)
    applyrelationlast!(td, work, coefficients, relations, trig)
    applyfft!(fd, td, rfftplan)
    return nothing
end

"""
    hbmatindices(frequencies::Frequencies{N},
        Amatrixmodes::AbstractMatrix{NTuple{N,Int}}; alias::Bool = false)

The index matrix of the harmonic balance matrix, which describes the
coupling between the modes: for the mode difference of each pair of modes
in `Amatrixmodes` (the [`ModeDifferences`](@ref) of the retained modes),
the element of the frequency domain array of the untruncated grid
`frequencies` (from the RFFT or FFT) that coupling reads. A negative index
means the complex conjugate of that element is taken, and a zero index
that the term is not present, so it is skipped. With `alias = true` a
difference which falls outside the sampled grid is aliased back onto it by
the periodicity of the discrete transform ([`aliasmode`](@ref)) rather than
dropped; the linearized solver uses `alias = false`, which makes the
assembled matrix an explicit truncation. The differences are an argument
so that a solve which needs the indices both with and without aliasing
forms them once.

# Examples
```jldoctest
pumpfreq = JosephsonCircuits.calcfreqsrdft((4,))
signalfreq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsdft((4,));
    dc=false,odd=true,even=false,maxintermodorder=2,
)
JosephsonCircuits.hbmatindices(pumpfreq,
    JosephsonCircuits.ModeDifferences(signalfreq.modes))

# output
4×4 Matrix{Int64}:
  1  -3  5   3
  3   1  0   5
 -5   0  1  -3
 -3  -5  3   1
```
```jldoctest
pumpfreq = JosephsonCircuits.calcfreqsrdft((4,))
signalfreq = JosephsonCircuits.truncfreqs(
    JosephsonCircuits.calcfreqsdft((4,));
    dc=false,odd=true,even=false,maxintermodorder=2,
)
JosephsonCircuits.hbmatindices(pumpfreq,
    JosephsonCircuits.ModeDifferences(signalfreq.modes); alias = true)

# output
4×4 Matrix{Int64}:
  1  -3   5   3
  3   1  -4   5
 -5   4   1  -3
 -3  -5   3   1
```
"""
function hbmatindices(frequencies::Frequencies{N},
    Amatrixmodes::AbstractMatrix{NTuple{N,Int}}; alias::Bool = false) where N
    modesdict = modeindexdict(frequencies)
    Nt = frequencies.Nt
    return [storedmodeindex(modesdict, mode, Nt, alias) for mode in Amatrixmodes]
end

# the position of each mode of the untruncated grid in the frequency domain
# array
function modeindexdict(frequencies::Frequencies)
    modesdict = Dict{eltype(frequencies.modes),Int}()
    for (i, mode) in enumerate(frequencies.modes)
        modesdict[mode] = i
    end
    return modesdict
end

# The signed position of `mode` in the frequency domain array whose modes
# `modesdict` indexes: negative for the complex conjugate of a stored mode
# and zero for a mode which is not stored. With `alias`, a mode which falls
# outside the grid is first taken back onto it by the periodicity of the
# transform.
function storedmodeindex(modesdict, mode::NTuple{N,Int}, Nt::NTuple{N,Int},
    alias::Bool) where N
    if alias
        aliasedmode, conjflag = aliasmode(mode, Nt)
        k = get(modesdict, aliasedmode, 0)
        return conjflag ? -k : k
    end
    k = get(modesdict, mode, 0)
    iszero(k) || return k
    return -get(modesdict, map(-, mode), 0)
end


"""
    aliasmode(mode::NTuple{N,Int}, Nt::NTuple{N,Int})

Alias the mode `mode` back onto the sampled grid with `Nt` time domain
samples along each dimension, returning the canonical stored mode and
whether the stored element is the complex conjugate of the requested
mode. The first dimension is the RDFT dimension of a real signal, which
stores only the modes from 0 to Nt[1] ÷ 2; a mode which aliases onto the
other half of that dimension is stored as the complex conjugate at the
negated mode. The other dimensions are full DFT dimensions whose stored
modes range from -(Nt[d]-1) ÷ 2 - iseven(Nt[d]) to (Nt[d]-1) ÷ 2.

# Examples
```jldoctest
julia> JosephsonCircuits.aliasmode((3,), (8,))
((3,), false)

julia> JosephsonCircuits.aliasmode((5,), (8,))
((3,), true)

julia> JosephsonCircuits.aliasmode((8,), (8,))
((0,), false)

julia> JosephsonCircuits.aliasmode((1, 3), (8, 4))
((1, -1), false)
```
"""
function aliasmode(mode::NTuple{N,Int}, Nt::NTuple{N,Int}) where N
    # wrap each dimension onto the sampled grid 0:Nt-1
    w1 = mod(mode[1], Nt[1])
    # the first dimension is the RDFT dimension of a real signal, which
    # stores only the modes from 0 to Nt[1] ÷ 2. a mode on the other half
    # of that dimension is stored as the complex conjugate at the negated
    # mode.
    conjflag = w1 > Nt[1] ÷ 2
    m = conjflag ? ntuple(d -> -mode[d], Val(N)) : mode
    # represent the full DFT dimensions with the signed modes used by the
    # frequency structures, where the wrapped indices above Nt ÷ 2 are
    # the negative modes.
    s = ntuple(Val(N)) do d
        wd = mod(m[d], Nt[d])
        if d == 1 || 2*wd <= Nt[d] - 1
            wd
        else
            wd - Nt[d]
        end
    end
    return s, conjflag
end

"""
    selfconjmodes(frequencies::Frequencies)

Return a vector of booleans indicating which of the modes of `frequencies`
are their own complex conjugates on the sampled grid, which is the case
when twice the mode aliases to zero along every dimension. The Fourier
coefficients of a real signal at those modes, such as dc, are purely
real.

# Examples
```jldoctest
julia> freq = JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(JosephsonCircuits.calcfreqsrdft((3,));dc=true,odd=true,even=true));JosephsonCircuits.selfconjmodes(freq)
4-element Vector{Bool}:
 1
 0
 0
 0
```
"""
function selfconjmodes(frequencies::Frequencies)
    Nt = frequencies.Nt
    return [all(d -> mod(2*m[d], Nt[d]) == 0, 1:length(Nt))
        for m in frequencies.modes]
end

"""
    hbconjmatind(truncfrequencies::Frequencies{N})

Returns a matrix describing which indices of the frequency domain matrix
(from the RFFT) to pull out and use in the conjugate harmonic balance
matrix, which is built from the sums of the modes, aliased back onto the
sampled grid, and describes the coupling between the modes and the
complex conjugates of the modes (the derivative of the residual with
respect to the complex conjugates of the node fluxes). A negative index
means we take the complex conjugate of that element. A zero index means
that term is not present, so skip it. See also [`hbmatindices`](@ref).

# Examples
```jldoctest
julia> freq = JosephsonCircuits.calcfreqsrdft((3,));JosephsonCircuits.hbconjmatind(JosephsonCircuits.removeconjfreqs(JosephsonCircuits.truncfreqs(freq;dc=true,odd=true,even=true,maxintermodorder=2)))
4×4 Matrix{Int64}:
 1   2   3   4
 2   3   4  -4
 3   4  -4  -3
 4  -4  -3  -2
```
"""
function hbconjmatind(truncfrequencies::Frequencies{N}) where N
    frequencies = calcfreqs(truncfrequencies.Nharmonics,
        truncfrequencies.Nw, truncfrequencies.Nt)
    return hbconjmatind(frequencies, truncfrequencies)
end

"""
    hbconjmatind(frequencies::Frequencies{N},
        truncfrequencies::Frequencies{N})

Returns a matrix describing which indices of the frequency domain matrix
(from the RFFT) to pull out and use in the conjugate harmonic balance
matrix, which is built from the sums of the modes
`truncfrequencies.modes[i] + truncfrequencies.modes[j]`, aliased back
onto the sampled grid described by `frequencies`. A negative index means
we take the complex conjugate of that element. A zero index means that
term is not present, so skip it. See also [`hbmatindices`](@ref).
"""
function hbconjmatind(frequencies::Frequencies{N},
    truncfrequencies::Frequencies{N}) where N

    truncmodes = truncfrequencies.modes
    modesdict = modeindexdict(frequencies)
    Nt = frequencies.Nt
    # the coupling between the modes and the complex conjugates of the
    # modes involves the sums of the modes, aliased back onto the sampled
    # grid of the untruncated frequencies; the sums are read and not kept
    return [storedmodeindex(modesdict, truncmodes[i] .+ truncmodes[j], Nt, true)
        for i in eachindex(truncmodes), j in eachindex(truncmodes)]
end
