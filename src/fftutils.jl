
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
                # (dc && all(nvals .== 0)) ||
                (dc && all(==(0),nvals)) ||
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
julia> JosephsonCircuits.calcdftsymmetries(4)
Dict{CartesianIndex{1}, CartesianIndex{1}} with 1 entry:
  CartesianIndex(2,) => CartesianIndex(4,)

julia> JosephsonCircuits.calcdftsymmetries((2,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 2 entries:
  CartesianIndex(2, 2) => CartesianIndex(2, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)

julia> JosephsonCircuits.calcdftsymmetries((3,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 4 entries:
  CartesianIndex(2, 3) => CartesianIndex(3, 2)
  CartesianIndex(2, 1) => CartesianIndex(3, 1)
  CartesianIndex(2, 2) => CartesianIndex(3, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)
```
"""
function calcdftsymmetries(Nt)

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

function calcdftsymmetries(A::AbstractArray)
    return calcdftsymmetries(size(A))
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

    d=calcdftsymmetries(N)

    z=zeros(Int,N)
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
    calcrdftsymmetries(Nt)

Calculate the conjugate symmetries in the multi-dimensional RDFT (DFT of
a real signal).

# Examples
```jldoctest
julia> JosephsonCircuits.calcrdftsymmetries(4)
Dict{CartesianIndex{1}, CartesianIndex{1}}()

julia> JosephsonCircuits.calcdftsymmetries(4)
Dict{CartesianIndex{1}, CartesianIndex{1}} with 1 entry:
  CartesianIndex(2,) => CartesianIndex(4,)

julia> JosephsonCircuits.calcrdftsymmetries((2,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 2 entries:
  CartesianIndex(2, 2) => CartesianIndex(2, 3)
  CartesianIndex(1, 2) => CartesianIndex(1, 3)

julia> JosephsonCircuits.calcrdftsymmetries((3,3))
Dict{CartesianIndex{2}, CartesianIndex{2}} with 1 entry:
  CartesianIndex(1, 2) => CartesianIndex(1, 3)

julia> JosephsonCircuits.calcrdftsymmetries((2,3,3))
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
function calcrdftsymmetries(Nt)

    Nw=NTuple{length(Nt),Int}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))

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
function printrdftsymmetries(Nt)

    Nw=NTuple{length(Nt),Int}(ifelse(i == 1, (val÷2)+1, val) for (i,val) in enumerate(Nt))

    d=calcrdftsymmetries(Nt)

    z=zeros(Int,Nw)
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

"""
    calcindexdict(N::Tuple)

Return a dictionary of Cartesian indices where the Cartesian index is the key
and the index giving the order is the value.

# Examples
```jldoctest
julia> JosephsonCircuits.calcindexdict(3)
Dict{CartesianIndex{1}, Int64} with 3 entries:
  CartesianIndex(2,) => 2
  CartesianIndex(3,) => 3
  CartesianIndex(1,) => 1

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

function calcindexdict(N::Number)
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

    indexdict = calcrdftsymmetries(Nt);
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
return the indices `conjsourceindices` whose data should be copied and complex
conjugated to `conjtargetindices`.

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

    indices = calcrdftsymmetries(Nt);
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
    phivectortomatrix!(phivector::Vector,phimatrix::Array,
        indexmap::Vector{Int},conjsourceindices::Vector{Int},
        conjtargetindices::Vector{Int},Nbranches::Int)

The harmonic balance method requires a vector with all of the conjugate symmetric
terms removed and potentially other terms dropped if specified by the user (
for example, intermodulation products which are not of interest) whereas the
Fourier transform operates on multidimensional arrays with the proper
conjugate symmetries and with dropped terms set to zero. This function converts
a vector to an array with the above properties.

"""
function phivectortomatrix!(phivector::Vector, phimatrix::Array,
    indexmap::Vector{Int}, conjsourceindices::Vector{Int},
    conjtargetindices::Vector{Int}, Nbranches::Int)

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

"""
function phimatrixtovector!(phivector::Vector, phimatrix::Array,
    indexmap::Vector{Int}, conjsourceindices::Vector{Int},
    conjtargetindices::Vector{Int}, Nbranches::Int)

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
