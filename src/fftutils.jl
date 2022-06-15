# w is a vector of frequencies
# n is a vector of number of modes
# dimensionality of the output is given by the length of w and n.
function calcfrequencies(w::Number,Nharmonics::Number;maxintermodorder=Inf,
    dc=true,even=true,odd=true)
    return calcfrequencies(Tuple(w),Tuple(Nharmonics),
        maxintermodorder=maxintermodorder,dc=dc,even=even,odd=odd)
end
function calcfrequencies(w::Tuple,Nharmonics::Tuple;maxintermodorder=Inf,
    dc=true,even=true,odd=true)
#     out = zeros(typeof(w1),nmax,mmax)
    if length(w) !== length(Nharmonics)
        error("Each frequency must have a number of harmonics.")
    end
    
    n = Nharmonics .+ 1
    
    # double the size of all but the first n because 
    # only the first axis of a multi-dimensional rfft has
    # only positive frequencies.
    Nw=NTuple{length(n),Int64}(ifelse(i == 1, val, 2*val-1) for (i,val) in enumerate(n))
    
    # store the coordinates in an array of CartesianIndex structures
    # not sure if i want to use an array of tuples or cartesianindices
    coords = Array{CartesianIndex{length(n)},1}(undef,0)
    dropcoords = Array{CartesianIndex{length(n)},1}(undef,0)

#     coords = Array{NTuple{length(n),Int64},1}(undef,0)

    # store the values in an array of the same type as w
    values = Array{eltype(w),1}(undef,0)
    dropvalues = Array{eltype(w),1}(undef,0)
    
    nvals = zeros(eltype(n),length(n))
#     oddonly = true

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
                (dc && all(nvals .== 0)) ||
                # test for even (and not DC)
                (even && mod(sum(abs,nvals),2) == 0 && sum(abs,nvals) > 0) ||
                # test for odd
                (odd && mod(sum(abs,nvals),2) == 1)
            ) && # test for containing only one frequency or less than maxintermodorder
                (sum( nvals .!== 0) == 1 || sum(abs,nvals) <= maxintermodorder)

            push!(coords,i)
            push!(values,dot(w,nvals))
        else
            push!(dropcoords,i)
            push!(dropvalues,dot(w,nvals))
        end
    end
    return Nw,coords,values,dropcoords,dropvalues
end

function vectortodense(coords,values,Nharmonics)
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
        error("Not designed to visualize higher dimensional arrays")
    end
    return s
end


"""
    calcdftsymmetries(N::Tuple)

Calculate the conjugate symmetries in the multi-dimensional DFT (or FFT).

# Examples
```jldoctest
julia> JosephsonCircuits.calcdftsymmetries((2,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 2 entries:
  CartesianIndex(2, 2) => CartesianIndex(2, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
```
```jldoctest
julia> JosephsonCircuits.calcdftsymmetries((3,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 4 entries:
  CartesianIndex(2, 3) => CartesianIndex(3, 2)
  CartesianIndex(2, 1) => CartesianIndex(3, 1)
  CartesianIndex(2, 2) => CartesianIndex(3, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
```
"""
function calcdftsymmetries(N::Tuple)

    d = Dict{CartesianIndex{length(N)},CartesianIndex{length(N)}}()
    for k in CartesianIndices(N)
        # check that none of the indices are equal but that they are valid
        # indices. 
        if any(k.I .- 1 .!= mod.(N .- (k.I .- 1),N)) && all(mod.(N .- (k.I .- 1),N) .< N)
            tmp = sort([k.I,mod.(N .- (k.I .- 1),N) .+ 1])
            if !haskey(d,CartesianIndex(tmp[1]))
                d[CartesianIndex(tmp[1])] = CartesianIndex(tmp[2])
            end
        end
    end
    return d
end

function calcdftsymmetries(A::AbstractArray)
    return calcdftsymmetries(size(A))
end

"""
    printdftsymmetries(N::Tuple)

Calculate the conjugate symmetries in the multi-dimensional DFT (or FFT).

# Examples
```jldoctest
julia> JosephsonCircuits.printdftsymmetries((3,3))
3×3 Matrix{Int64}:
  0   2  -2
  1   3   4
 -1  -4  -3
```
```jldoctest
julia> JosephsonCircuits.printdftsymmetries((4,3))
4×3 Matrix{Int64}:
  0   2  -2
  1   3   5
  0   4  -4
 -1  -5  -3
```
```jldoctest
julia> JosephsonCircuits.printdftsymmetries((3,4))
3×4 Matrix{Int64}:
  0   2   0  -2
  1   3   4   5
 -1  -5  -4  -3
```
"""
function printdftsymmetries(N::Tuple)

    d=calcdftsymmetries(N)

    z=zeros(Int64,N)
    i = 1
    for (key,val) in sort(OrderedDict(d))
        z[key] = i
        z[val] = -i
        i+=1
    end
    return z
end

function printdftsymmetries(A::AbstractArray)
    return printdftsymmetries(size(A))
end

"""
    calcrdftsymmetries(N::Tuple)

Calculate the conjugate symmetries in the multi-dimensional RDFT (DFT of
a real signal).

# Examples
```jldoctest
julia> JosephsonCircuits.calcrdftsymmetries((2,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 2 entries:
  CartesianIndex(2, 2) => CartesianIndex(2, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
```
```jldoctest
julia> JosephsonCircuits.calcrdftsymmetries((3,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 1 entry:
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
```
"""
function calcrdftsymmetries(Nt::Tuple)

    Nw=NTuple{length(Nt),Int64}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))


    d = Dict{CartesianIndex{length(Nt)},CartesianIndex{length(Nt)}}()
    for k in CartesianIndices(Nw)
        # check that none of the indices are equal but that they are valid
        # indices. 
        if any(k.I .- 1 .!= mod.(Nt .- (k.I .- 1),Nt)) && all(mod.(Nt .- (k.I .- 1),Nt) .< Nw)
            tmp = sort([k.I,mod.(Nt .- (k.I .- 1),Nt) .+ 1])
            if !haskey(d,CartesianIndex(tmp[1]))
                d[CartesianIndex(tmp[1])] = CartesianIndex(tmp[2])
            end
        end
    end
    return d
end

function calcrdftsymmetries(A::AbstractArray)
    return calcdftsymmetries(size(A))
end


"""
    printrdftsymmetries(N::Tuple)

Calculate the conjugate symmetries in the multi-dimensional RDFT (DFT of a real
signal).

# Examples
```jldoctest
julia> JosephsonCircuits.printrdftsymmetries((3,3))
2×3 Matrix{Int64}:
 0  1  -1
 0  0   0
```
```jldoctest
julia> JosephsonCircuits.printrdftsymmetries((4,3))
3×3 Matrix{Int64}:
 0  1  -1
 0  0   0
 0  2  -2
```
```jldoctest
julia> JosephsonCircuits.printrdftsymmetries((3,4))
2×4 Matrix{Int64}:
 0  1  0  -1
 0  0  0   0
```
"""
function printrdftsymmetries(Nt::Tuple)

    Nw=NTuple{length(Nt),Int64}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))

    d=calcrdftsymmetries(Nt)

    z=zeros(Int64,Nw)
    i = 1
    for (key,val) in sort(OrderedDict(d))
        z[key] = i
        z[val] = -i
        i+=1
    end
    return z
end

function printrdftsymmetries(A::AbstractArray)
    return printrdftsymmetries(size(A))
end

function calcindexdict(N)
#     d = Dict()
    d = Dict{CartesianIndex{length(N)},Int64}()

    for (i,index) in enumerate(CartesianIndices(N))
        d[index] = i
    end
    return d
end



function calcfrequencies2(Nt,coords,values)

    indexdict = calcrdftsymmetries(Nt);
    N2=NTuple{length(Nt),Int64}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))
    reverseindexdict = Dict{CartesianIndex{length(Nt)},CartesianIndex{length(Nt)}}()
    for (key,val) in indexdict
        reverseindexdict[val] = key
    end

    values2 = Array{eltype(values),1}(undef,0)

    carttoint = calcindexdict(N2)
    inttocart = CartesianIndices(N2)

    for (i,coord) in enumerate(coords)
        if !haskey(reverseindexdict,coord)
            push!(values2,values[i])
        end
    end

    return values2
end



function calcphiindices(N,dropindices)

    indices = calcrdftsymmetries(N);
    N2=NTuple{length(N),Int64}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(N))
    conjindices = Dict{CartesianIndex{length(N)},CartesianIndex{length(N)}}()
    for (key,val) in indices
        conjindices[val] = key
    end

    indexmap = Array{Int64,1}(undef,0)

    # the index of the element of the N dimensional array in the frequency domain
    # that i should copy to conjtargetindices and take the complex conjugate of. 
    conjsourceindices = Array{Int64,1}(undef,0)

    # the index of the element of the N dimensional array in the frequency domain
    # that i should take the complex conjugate of
    conjtargetindices = Array{Int64,1}(undef,0)

    carttoint = calcindexdict(N2)
    inttocart = CartesianIndices(N2)

    # generate the index maps to convert between the vector
    # and matrix
    for (j,index) in enumerate(CartesianIndices(N2))
        # don't include the index in the maps if it shows up
        # in dropindices
        if !haskey(dropindices,index)
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


function phivectortomatrix!(phivector::Vector,phimatrix::Array,
    indexmap::Vector{Int64},conjsourceindices::Vector{Int64},
    conjtargetindices::Vector{Int64},Nbranches::Int64)


    if length(indexmap)*Nbranches != length(phivector)
        error("Unexpected length for phivector")
    end

    Nvector = length(phivector)÷ Nbranches
    Nmatrix = prod(size(phimatrix)[1:end-1])

    # fill the matrix with zeros
    fill!(phimatrix,0)

    for i in 1:Nbranches
    # @inbounds for i in 1:Nbranches
        for j in 1:length(indexmap)
            phimatrix[indexmap[j]+(i-1)*Nmatrix] = phivector[j+(i-1)*Nvector]
        end
    end

    for i in 1:Nbranches
    # @inbounds for i in 1:Nbranches
        for j in 1:length(conjtargetindices)
            phimatrix[conjtargetindices[j]+ (i-1)*Nmatrix] = conj(phimatrix[conjsourceindices[j]+ (i-1)*Nmatrix])
        end
    end
    return nothing
end


function phimatrixtovector!(phivector::Vector,phimatrix::Array,
    indexmap::Vector{Int64},conjsourceindices::Vector{Int64},
    conjtargetindices::Vector{Int64},Nbranches::Int64)

    Nvector = length(phivector)÷ Nbranches
    Nmatrix = prod(size(phimatrix)[1:end-1])

    # fill the vector with zeros
    fill!(phivector,0)


    # @inbounds for i in 1:Nbranches
    for i in 1:Nbranches
        for j in 1:length(indexmap)
            phivector[j+(i-1)*Nvector] = phimatrix[indexmap[j]+(i-1)*Nmatrix]
        end
    end
    return nothing
end
