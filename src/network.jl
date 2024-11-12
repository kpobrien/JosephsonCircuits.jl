
"""
    connectS(Sx::AbstractArray,k::Int,l::Int)

Connect ports `k` and `l` on the same `m` port microwave network represented by
the scattering parameter matrix `Sx`, resulting in an `(m-2)` port network, as illustrated below:

Input network:
```
      m |         | l+1    
        |   ...   |         l
        |_________|__________ 
        |         |          |
        |   Sx    |  ...     |
        |  m x m  |          |
    ____|_________|_____ k+1 |
    1   |   ...   |          |
        |         | k        |
      2 |         |__________|
```
Output network:
```
    m-2 |         | l-1     
        |         |         
        |   ...   |         
        |_________|         
        |         |         
        |    S    |  ...    
        |m-2 x m-2|         
    ____|_________|_________
    1   |   ...         k   
        |                   
        |                   
      2 |                   
```
# Arguments
- `Sx::Array`: Array of scattering parameters representing the network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `k::Int`: First port to connect, with one based indexing.
- `l::Int`: Second port to connect, with one based indexing.

# References
Haifang Liao and Wayne Wei-Ming Dai, "Capturing Time-of-flight Delay For
Transient Analysis Based On Scattering Parameter Macromodel,"
IEEE/ACM International Conference on Computer-Aided Design,
San Jose, CA, USA, 1994, pp. 412-417, doi: 10.1109/ICCAD.1994.629836.

R. C. Compton, "Perspectives in microwave circuit analysis," Proceedings of
the 32nd Midwest Symposium on Circuits and Systems,
Champaign, IL, USA, 1989, pp. 716-718 vol.2, doi: 10.1109/MWSCAS.1989.101955.

G. Filipsson, "A New General Computer Algorithm for S-Matrix Calculation of
Interconnected Multiports," 1981 11th European Microwave Conference,
Amsterdam, Netherlands, 1981, pp. 700-704, doi: 10.1109/EUMA.1981.332972.
"""
function connectS(Sx::AbstractArray{T,N},k::Int,l::Int) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeS = NTuple{N}(ifelse(i<=2,size(Sx,i)-2,size(Sx,i)) for i in 1:ndims(Sx))

    # allocate an array of zeros of the same type as Sx
    Sout = similar(Sx,sizeS)

    # remove the self loop
    connectS!(Sout,Sx,k,l)

    return Sout
end

"""
    connectS!(Sout,Sx,k::Int,l::Int)

See [`connectS`](@ref) for description.

"""
function connectS!(Sout,Sx,k::Int,l::Int)

    # validate all of the inputs
    if ndims(Sx) != ndims(Sout)
        throw(DimensionMismatch("`Sout` and `Sx` must have the same number of dimensions."))
    end

    if ndims(Sx) < 2
        throw(DimensionMismatch("`Sout` and `Sx` must have at least two dimensions."))
    end

    if size(Sx,1) != size(Sx,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sx` must be equal."))
    end

    if size(Sout,1) != size(Sout,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."))
    end

    if size(Sx,1) -2 != size(Sout,1)
        throw(DimensionMismatch("Length of first two dimensions must be 2 smaller for `Sout` than `Sx` because we are merging two ports."))
    end

    for i in 3:ndims(Sx)
        if size(Sx,i) != size(Sout,i)
            throw(DimensionMismatch("Non-port axis lengths of `Sx` and `Sout` must be equal."))
        end
    end

    if k > size(Sx,1)
        throw(ArgumentError("Port `k` is larger than number of ports in `Sx`."))
    end

    if l > size(Sx,1)
        throw(ArgumentError("Port `l` is larger than number of ports in `Sx`."))
    end

    if l < 1
        throw(ArgumentError("Port `l` is smaller than one."))
    end

    if k < 1
        throw(ArgumentError("Port `k` is smaller than one."))
    end

    if l == k
        throw(ArgumentError("`k` and `l` cannot be equal because a port cannot be merged with itself."))
    end
  
    # the number of ports in the input matrix
    m = size(Sx,1)

    # make the indices so we can skip k and l
    xindices = zeros(Int,m-2)
    iout = 0
    for i in 1:m
        if i != k && i != l
            # if the current index is neither
            # k nor l, then keep it
            iout+=1
            xindices[iout] = i
        end
    end
  
    # loop over the dimensions of the array greater than 2
    for i in CartesianIndices(axes(Sout)[3:end])
        connectS_inner!(view(Sout,:,:,i),view(Sx,:,:,i),k,l,m,xindices)
    end

    return nothing
end

"""
    connectS_inner!(Sout,Sx,k,l,m,xindices)

See [`connectS`](@ref) for description.

"""
function connectS_inner!(Sout,Sx,k,l,m,xindices)

    # Eq. 16.3
    oneoverdelta = 1/((1 - Sx[l,k])*(1 - Sx[k,l]) - Sx[l,l]*Sx[k,k])

    # generate the scattering parameters
    # by looping over the output matrix indices
    for i in 1:m-2

        # input matrix index
        xi = xindices[i]

        # Eq. 16.1, 16.2
        al = (Sx[l,xi]*Sx[k,k]+Sx[k,xi]*(1 - Sx[l,k]))*oneoverdelta
        ak = (Sx[k,xi]*Sx[l,l]+Sx[l,xi]*(1 - Sx[k,l]))*oneoverdelta

        for j in 1:m-2

            # input matrix index
            xj = xindices[j]

            # Eq. 15
            Sout[j,i] =Sx[xj,xi] + Sx[xj,l]*al + Sx[xj,k]*ak
        end
    end
    return nothing
end


"""
    connectS(Sx::AbstractArray,Sy::AbstractArray,k::Int,l::Int)

Connect port `k` on an `m` port network, represented by the scattering
parameter matrix `Sx`, to port `l` on an `n` port network, represented by the
scattering parameter matrix `Sy`, resulting in a single `(m+n-2)` port
network, as illustrated below:

Input network:
```
      m |        | k+1                       | 2
        |        |                           |
        |   ...  |                     ...   |
        |________|                  _________|________
        |        |                  |        |       1
        |   Sx   |                  |   Sy   |
        |  m x m |                  |  n x n |
    ____|________|__________________|________|
    1   |   ...     k           l   |   ...  |
        |                           |        |
        |                           |        |
      2 |                       l+1 |        | n
```
Output network:
```
    m-1 |        | k      | m+1    
        |        |        |        
        |   ...  |   ...  |        
        |________|________|________
        |                 |     m  
        |        S        |        
        |  m+n-2 x m+n-2  |        
    ____|_________________|        
    1   |   ...  |   ...  |        
        |        |        |        
        |        |        |        
      2 |        |        |  m+n-2 
                m-1+l              
```
# Arguments
- `Sx::Array`: Array of scattering parameters representing the first network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `Sy::Array`: Array of scattering parameters representing the second network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `k::Int`: Port on first network, with one based indexing.
- `l::Int`: Port on second network, with one based indexing.

# References
Haifang Liao and Wayne Wei-Ming Dai, "Capturing Time-of-flight Delay For
Transient Analysis Based On Scattering Parameter Macromodel,"
IEEE/ACM International Conference on Computer-Aided Design,
San Jose, CA, USA, 1994, pp. 412-417, doi: 10.1109/ICCAD.1994.629836.

R. C. Compton, "Perspectives in microwave circuit analysis," Proceedings of
the 32nd Midwest Symposium on Circuits and Systems,
Champaign, IL, USA, 1989, pp. 716-718 vol.2, doi: 10.1109/MWSCAS.1989.101955.

G. Filipsson, "A New General Computer Algorithm for S-Matrix Calculation of
Interconnected Multiports," 1981 11th European Microwave Conference,
Amsterdam, Netherlands, 1981, pp. 700-704, doi: 10.1109/EUMA.1981.332972.

"""
function connectS(Sx::AbstractArray{T,N},Sy::AbstractArray{T,N},k::Int,l::Int) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeSx = size(Sx)
    sizeSy = size(Sy)
    sizeS = NTuple{N}(ifelse(i<=2,sizeSx[i]+sizeSy[i]-2,sizeSx[i]) for i in 1:length(sizeSx))

    # allocate an array of zeros of the same type as Sx
    Sout = similar(Sx,sizeS)

    # connect the networks
    connectS!(Sout,Sx,Sy,k,l)

    return Sout
end

"""
    connectS!(Sout,Sx,Sy,k,l)

See [`connectS`](@ref) for description.

"""
function connectS!(Sout,Sx,Sy,k::Int,l::Int)
  

    # validate all of the inputs
    if ndims(Sx) != ndims(Sy)
        throw(DimensionMismatch("`Sx` and `Sy` must have the same number of dimensions."))
    end

    if ndims(Sx) != ndims(Sout)
        throw(DimensionMismatch("`Sout`, `Sx`, and `Sy` must have the same number of dimensions."))
    end

    if ndims(Sx) < 2
        throw(DimensionMismatch("`Sout`, `Sx`, and `Sy` must have atleast two dimensions."))
    end

    if size(Sx,1) != size(Sx,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sx` must be equal."))
    end

    if size(Sy,1) != size(Sy,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sy` must be equal."))
    end

    if size(Sout,1) != size(Sout,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."))
    end

    if size(Sx,1) + size(Sy,1) - 2 != size(Sout,1)
        throw(DimensionMismatch("First two dimensions of `Sout` must be `m+n-2`."))
    end

    for i in 3:ndims(Sx)
        if size(Sx,i) != size(Sout,i)
            throw(DimensionMismatch("Non-port axis lengths of `Sx`, `Sy, and `Sout` must be equal."))
        end
    end

    if k > size(Sx,1)
        throw(ArgumentError("Port `k` is larger than number of ports in `Sx`."))
    end

    if l > size(Sy,1)
        throw(ArgumentError("Port `l` is larger than number of ports in `Sy`."))
    end

    if l < 1
        throw(ArgumentError("Port `l` is smaller than one."))
    end

    if k < 1
        throw(ArgumentError("Port `k` is smaller than one."))
    end

  
    # the number of ports in the input matrix
    m = size(Sx,1)
    n = size(Sy,1)
    
    # make the indices so we can skip k
    xindices = zeros(Int,m-1)
    for iout in 1:k-1
        xindices[iout] = iout
    end
    for iout in k:m-1
        xindices[iout] = iout+1
    end

    # make the indices so we can skip l
    yindices = zeros(Int,n-1)
    for iout in 1:l-1
        yindices[iout] = iout
    end
    for iout in l:n-1
        yindices[iout] = iout+1
    end

    # loop over the dimensions of the array greater than 2
    for i in CartesianIndices(axes(Sout)[3:end])
        connectS_inner!(view(Sout,:,:,i),view(Sx,:,:,i),view(Sy,:,:,i),k,l,m,
            n,xindices,yindices)
    end

    return nothing
end

"""
    connectS_inner!(Sout,Sx,Sy,k,l,m,n,xindices,yindices)

See [`connectS`](@ref) for description.

"""
function connectS_inner!(Sout,Sx,Sy,k,l,m,n,xindices,yindices)

    # use a separate loop for each
    # quadrant of the output matrix

    # calculate the inverse of the denominator
    oneoverdenom = 1/(1-Sx[k,k]*Sy[l,l])

    # upper left quadrant, i,j in Sx
    for i in 1:m-1
        for j in 1:m-1
            # indices for the input matrices
            xi = xindices[i]
            xj = xindices[j] 

            # Eq. 10.1
            Sout[j,i] = Sx[xj,xi] + Sx[k,xi]*Sy[l,l]*Sx[xj,k]*oneoverdenom
        end
    end

    # upper right  quadrant, i in Sy, j in Sx
    for i in m:m+n-2
        for j in 1:m-1
            # indices for the input matrices
            yi = yindices[i - m + 1]
            xj = xindices[j]

            # Eq. 10.3
            Sout[j,i] = Sx[xj,k]*Sy[l,yi]*oneoverdenom
        end
    end

    # lower left quadrant, i in Sx, j in Sy
    for i in 1:m-1
        for j in m:m+n-2
            # indices for the input matrices
            xi = xindices[i]
            yj = yindices[j - m + 1]

            # Eq. 10.3
            Sout[j,i] = Sx[k,xi]*Sy[yj,l]*oneoverdenom 
        end
    end

    # lower right quadrant, i,j in Sy
    for i in m:m+n-2
        for j in m:m+n-2
            # indices for the input matrices
            yi = yindices[i - m + 1]
            yj = yindices[j - m + 1] 
              
            # Eq. 10.2
            Sout[j,i] = Sy[yj,yi] + Sy[l,yi]*Sx[k,k]*Sy[yj,l]*oneoverdenom
        end
    end
    return nothing
end

"""
    cascadeS(Sa,Sb)

Cascade the scattering parameter matrix `Sa` with the scattering matrix `Sb`
and return the combined scattering matrix.

# Examples
```jldoctest
julia> Sa = [0.0 0.5;0.5 0.0];Sb = [0.0 0.1;0.1 0.0];JosephsonCircuits.cascadeS(Sa,Sb)
2×2 Matrix{Float64}:
 0.0   0.05
 0.05  0.0

julia> Sa=rand(Complex{Float64},2,2);Sb=rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.cascadeS(Sa,Sb),JosephsonCircuits.AtoS(JosephsonCircuits.StoA(Sa)*JosephsonCircuits.StoA(Sb)))
true

julia> Sa=[rand(Complex{Float64},2,2) for i in 1:10];Sb=[rand(Complex{Float64},2,2) for i in 1:10];isapprox(JosephsonCircuits.cascadeS.(Sa,Sb),JosephsonCircuits.AtoS.(JosephsonCircuits.StoA.(Sa).*JosephsonCircuits.StoA.(Sb)))
true
```

# References
D. J. R. Stock and L. J. Kaplan, "A Comment on the Scattering Matrix of 
Cascaded 2n-Ports (Correspondence)," in IRE Transactions on Microwave 
Theory and Techniques, vol. 9, no. 5, pp. 454-454, September 1961, doi:
10.1109/TMTT.1961.1125369 .
"""
function cascadeS(Sa,Sb)

    S = similar(Sa)
    # make a view of T,S and loop
    # make a temporary array 

    # assume the port impedances are all the same for all ports and
    # frequencies. loop over the dimensions of the array greater than 2
    for i in CartesianIndices(axes(S)[3:end])
        cascadeS!(view(S,:,:,i),view(Sa,:,:,i),view(Sb,:,:,i))
    end

    return S

end

"""
    cascadeS!(S, Sa, Sb)

See [`cascadeS`](@ref) for description.

"""
function cascadeS!(S::AbstractMatrix, Sa::AbstractMatrix, Sb::AbstractMatrix)
    range1 = 1:size(Sa,1)÷2
    range2 = size(Sa,1)÷2+1:size(Sa,1)

    # make views of the block matrices
    S11a = view(Sa,range1,range1)
    S12a = view(Sa,range1,range2)
    S21a = view(Sa,range2,range1)
    S22a = view(Sa,range2,range2)

    S11b = view(Sb,range1,range1)
    S12b = view(Sb,range1,range2)
    S21b = view(Sb,range2,range1)
    S22b = view(Sb,range2,range2)

    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    tmp1 = (I-S22a*S11b)\S21a
    tmp2 = (I-S11b*S22a)\S12b
    S11 .= S11a + S12a*S11b*tmp1 # typo in ref: S21 should be S21a
    S12 .= S12a*tmp2 # typo in ref: S11a should be S11a
    S21 .= S21b*tmp1
    S22 .= S22b + S21b*S22a*tmp2

    return nothing
end

# generate the non in-place versions of the network parameter conversion
# functions.

# functions without an impedance argument
for f in [:StoT, :TtoS, :AtoB, :BtoA, :AtoZ, :ZtoA, :AtoY, :YtoA, :BtoY, :YtoB, :BtoZ, :ZtoB, :ZtoY, :YtoZ]
    @eval function ($f)(x::AbstractMatrix)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # evaluate the in place version of the function
        ($(Symbol(String(f)*"!")))(y,x,tmp)
        return y
    end
    
    @eval function ($f)(x::AbstractArray)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # evaluate the in place version of the function
        for i in CartesianIndices(axes(x)[3:end])
            ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp)
        end
        return y
    end
end

# functions using portimpedances
for f in [:ABCDtoS, :StoABCD]
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # evaluate the in place version of the function
        ($(Symbol(String(f)*"!")))(y,x,portimpedances)
        return y
    end
    
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if ndims(portimpedances) == 0
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),portimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),view(portimpedances,:,i))
            end
        end

        return y
    end
end

# functions using sqrtportimpedances
for f in [:StoZ, :YtoS]
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,sqrtportimpedances)
        else
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,Diagonal(sqrtportimpedances))
        end
        return y
    end
    
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if ndims(portimpedances) == 0
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,sqrtportimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,Diagonal(view(sqrtportimpedances,:,i)))
            end
        end

        return y
    end
end


# functions using oneoversqrtportimpedances
for f in [:StoY, :ZtoS]
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # calculate the inverse of the square root of the port impedances
        oneoversqrtportimpedances = 1 ./sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,oneoversqrtportimpedances)
        else
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,Diagonal(oneoversqrtportimpedances))
        end
        return y
    end
    
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # calculate the inverse of the square root of the port impedances
        oneoversqrtportimpedances = 1 ./sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if ndims(portimpedances) == 0
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,oneoversqrtportimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,Diagonal(view(oneoversqrtportimpedances,:,i)))
            end
        end

        return y
    end
end

# functions using sqrtportimpedances split into two
for f in [:StoA, :AtoS, :StoB, :BtoS]
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,sqrtportimpedances,
                sqrtportimpedances)
        else
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,
                Diagonal(sqrtportimpedances[1:size(sqrtportimpedances,1)÷2]),
                Diagonal(sqrtportimpedances[size(sqrtportimpedances,1)÷2+1:size(sqrtportimpedances,1)]))
        end
        return y
    end
    
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if ndims(portimpedances) == 0
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,
                    sqrtportimpedances, sqrtportimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,
                Diagonal(view(sqrtportimpedances,1:size(sqrtportimpedances,1)÷2,i)),
                Diagonal(view(sqrtportimpedances,size(sqrtportimpedances,1)÷2+1:size(sqrtportimpedances,1),i)))
            end
        end

        return y
    end
end


@doc """
    StoT(S)

Convert the scattering parameter matrix `S` to a transmission matrix
`T` and return the result.

# Examples
```jldoctest
julia> S = [0.0 1.0;1.0 0.0];JosephsonCircuits.StoT(S)
2×2 Matrix{Float64}:
  1.0  -0.0
 -0.0   1.0

julia> @variables S11 S12 S21 S22 Z0;S = [S11 S12;S21 S22];JosephsonCircuits.StoT(S)
2×2 Matrix{Num}:
 S12 + (-S11*S22) / S21  S11 / S21
           (-S22) / S21    1 / S21
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoT


"""
    StoT!(T::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix)

See [`StoT`](@ref) for description.

"""
function StoT!(T::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(T,1)÷2
    range2 = size(T,1)÷2+1:size(T,1)

    # tmp = [-I S11; 0 S21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp[range1,range2] .= S[range1, range1]
    tmp[range2,range2] .= S[range2, range1]

    # T = [-S12 0; -S22 I]
    fill!(T,zero(eltype(T)))
    for d in range2
        T[d,d] = 1
    end

    T[range1,range1] .= -S[range1, range2]
    T[range2,range1] .= -S[range2, range2]

    # perform the left division
    # T = inv(tmp)*T = [-I S11; 0 S21] \ [-S12 0; -S22 I]
    T .= tmp \ T

    return nothing
end


@doc """
    StoZ(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to an impedance parameter matrix
`Z` and return the result. Assumes a port impedance of 50 Ohms unless
specified with the `portimpedances` keyword argument.

# Examples
```jldoctest
julia> S = [0.0 0.0;0.0 0.0];JosephsonCircuits.StoZ(S)
2×2 Matrix{Float64}:
 50.0   0.0
  0.0  50.0

julia> S = [0.0 0.999;0.999 0.0];JosephsonCircuits.StoZ(S)
2×2 Matrix{Float64}:
 49975.0  49975.0
 49975.0  49975.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoZ

"""
    StoZ!(Z::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

See [`StoZ`](@ref) for description.

"""
function StoZ!(Z::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

    # tmp = (I - S)
    copy!(tmp,S)
    rmul!(tmp,-1)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end
    
    # Z = (I + S)*sqrt(portimpedances)
    copy!(Z,S)
    for d in 1:size(Z,1)
        Z[d,d] += 1
    end
    rmul!(Z,sqrtportimpedances)

    # perform the left division
    # Z = inv(tmp)*Z = (I - S) \ ((I + S)*sqrt(portimpedances))
    Z .= tmp \ Z

    # left multiply by sqrt(z)
    # compute sqrt(portimpedances)*((I - S) \ ((I + S)*sqrt(portimpedances)))
    lmul!(sqrtportimpedances,Z)
    return nothing
end

@doc """
    StoY(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to an admittance parameter matrix
`Y` and return the result. Assumes a port impedance of 50 Ohms unless
specified with the `portimpedances` keyword argument.

# Examples
```jldoctest
julia> S = [0.0 0.0;0.0 0.0];JosephsonCircuits.StoY(S)
2×2 Matrix{Float64}:
  0.02  -0.0
 -0.0    0.02

julia> S = [0.0 0.999;0.999 0.0];JosephsonCircuits.StoY(S)
2×2 Matrix{Float64}:
  19.99  -19.99
 -19.99   19.99
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoY

"""
    StoY!(Y::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportadmittances)

See [`StoY`](@ref) for description.

"""
function StoY!(Y::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)
    

    # tmp = (I + S)
    copy!(tmp,S)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # Y = (I - S)*oneoversqrtportimpedances
    copy!(Y,S)
    rmul!(Y,-1)
    for d in 1:size(Y,1)
        Y[d,d] += 1
    end
    rmul!(Y,oneoversqrtportimpedances)


    # perform the left division
    # Y = inv(tmp)*Y = (I + S) \ ((I - S)*oneoversqrtportimpedances)
    Y .= tmp \ Y

    # left multiply by sqrt(z)
    # compute oneoversqrtportimpedances*((I - S) \ ((I + S)*oneoversqrtportimpedances))
    lmul!(oneoversqrtportimpedances,Y)
    return nothing
end


@doc """
    StoA(S)

Convert the scattering parameter matrix `S` to the chain (ABCD) matrix `A` and
return the result.

# Examples
```jldoctest
julia> S = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.ZtoA(JosephsonCircuits.StoZ(S)),JosephsonCircuits.StoA(S))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoA

"""
    StoA!(A::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`StoA`](@ref) for description.

"""
function StoA!(A::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    A11 = view(A,range1,range1)
    A12 = view(A,range1,range2)
    A21 = view(A,range2,range1)
    A22 = view(A,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # define the matrices

    # tmp = -[(I-S11)/g1 -(I+S11)*g1; -S21*g1 -S21*g1]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    tmp11 .= -(I - S11)/sqrtportimpedances1
    tmp12 .= (I + S11)*sqrtportimpedances1
    tmp21 .= S21/sqrtportimpedances1
    tmp22 .= S21*sqrtportimpedances1

    # A = [-S12/g2 S12*g2; (I-S22)/g2 (I+S22)*g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    A11 .= -S12/sqrtportimpedances2
    A12 .= S12*sqrtportimpedances2
    A21 .= (I-S22)/sqrtportimpedances2
    A22 .= (I+S22)*sqrtportimpedances2

    # perform the left division

    # A = inv(tmp)*A
    A .= tmp \ A

    return nothing
end


@doc """
    StoB(S)

Convert the scattering parameter matrix `S` to the inverse chain (ABCD) matrix
`B` and return the result. Note that despite the name, the inverse of the chain
matrix is not equal to the inverse chain matrix, inv(A) ≠ B.

# Examples
```jldoctest
julia> S = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.BtoS(JosephsonCircuits.StoB(S)),S)
true

julia> S = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.AtoS(JosephsonCircuits.BtoA(JosephsonCircuits.StoB(S))),S)
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoB

"""
    StoB!(B::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`StoB`](@ref) for description.

"""
function StoB!(B::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    B11 = view(B,range1,range1)
    B12 = view(B,range1,range2)
    B21 = view(B,range2,range1)
    B22 = view(B,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # define the matrices

    # tmp = [S12/g2 S12*g2; -(I-S22)/g2 (I+S22)*g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    tmp11 .= S12/sqrtportimpedances2
    tmp12 .= S12*sqrtportimpedances2
    tmp21 .= -(I-S22)/sqrtportimpedances2
    tmp22 .= (I+S22)*sqrtportimpedances2

    # B = [(I-S11)/g1 (I+S11)*g1; -S21/g1 S21*g1]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    B11 .= (I-S11)/sqrtportimpedances1
    B12 .= (I+S11)*sqrtportimpedances1
    B21 .= -S21/sqrtportimpedances1
    B22 .= S21*sqrtportimpedances1

    # perform the left division

    # B = inv(tmp)*B
    B .= tmp \ B

    return nothing
end


@doc """
    TtoS(T)

Convert the transmission matrix `T` to a scattering parameter matrix `S` and return the result.

# Examples
```jldoctest
julia> T = [1.0 0.0;0.0 1.0];JosephsonCircuits.TtoS(T)
2×2 Matrix{Float64}:
 -0.0   1.0
  1.0  -0.0

julia> @variables T11 T12 T21 T22;T = [T11 T12;T21 T22];JosephsonCircuits.TtoS(T)
2×2 Matrix{Num}:
 T12 / T22  T11 + (-T12*T21) / T22
   1 / T22            (-T21) / T22
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of sign on T11 and T21 terms (suspected typo).
""" TtoS

"""
    TtoS!(S::AbstractMatrix,T::AbstractMatrix,tmp::AbstractMatrix)

See [`TtoS`](@ref) for description.

"""
function TtoS!(S::AbstractMatrix,T::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(T,1)÷2
    range2 = size(T,1)÷2+1:size(T,1)

    # make views of the block matrices
    # S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    # S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    T11 = view(T,range1,range1)
    T12 = view(T,range1,range2)
    T21 = view(T,range2,range1)
    T22 = view(T,range2,range2)

    # tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    # tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # tmp = [-I T12; 0 T22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp12 .= T12
    tmp22 .= T22

    # S = [0 T11; I T21]
    fill!(S,zero(eltype(S)))
    for d in range1
        S[d+size(T,1)÷2,d] = 1
    end

    S12 .= -T11
    S22 .= -T21

    # perform the left division
    # S = inv(tmp)*S
    S .= tmp \ S

    return nothing
end


@doc """
    ZtoS(Z;portimpedances=50.0)

Convert the impedance parameter matrix `Z` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances. Assumes a port impedance of 50 Ohms unless specified with
the `portimpedances` keyword argument.

# Examples
```jldoctest
julia> Z = [50.0 0.0;0.0 50.0];JosephsonCircuits.ZtoS(Z)
2×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0

julia> Z = [0.0 0.0;0.0 0.0];JosephsonCircuits.ZtoS(Z)
2×2 Matrix{Float64}:
 -1.0   0.0
  0.0  -1.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ZtoS

"""
    ZtoS!(S::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

See [`ZtoS`](@ref) for description.

"""
function ZtoS!(S::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)
    
    # compute \tilde{Z} = oneoversqrtportimpedances*Z*oneoversqrtportimpedances in S and tmp
    copy!(S,Z)
    rmul!(S,oneoversqrtportimpedances)
    lmul!(oneoversqrtportimpedances,S)
    copy!(tmp,S)

    # tmp = (\tilde{Z} + I)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # S = (\tilde{Z} - I)
    for d in 1:size(S,1)
        S[d,d] -= 1
    end

    # perform the left division
    # S = inv(tmp)*S = (\tilde{Z} + I) \ (\tilde{Z} - I)
    S .= tmp \ S

    return nothing
end

function ZtoY!(Y, Z, tmp)
    Y .= inv(Z)
    return nothing
end

function YtoZ!(Z, Y, tmp)
    Z .= inv(Y)
    return nothing
end

@doc """
    ZtoA(Z)

Convert the impedance matrix `Z` to the ABCD matrix `A` and return the result.

# Examples
```jldoctest
julia> Z=[50.0 50.0;50.0 50.0];JosephsonCircuits.ZtoA(Z)
2×2 Matrix{Float64}:
 1.0   -0.0
 0.02   1.0

julia> @variables Z11 Z12 Z21 Z22;JosephsonCircuits.AtoZ([Z11 Z12;Z21 Z22])
2×2 Matrix{Num}:
 Z11 / Z21  -Z12 + (Z11*Z22) / Z21
   1 / Z21               Z22 / Z21
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ZtoA

"""
    ZtoA!(A::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)

See [`ZtoA`](@ref) for description.

"""
function ZtoA!(A::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)
    return AtoZ!(A,Z,tmp)
end

@doc """
    ZtoB(Z)

Convert the impedance matrix `Z` to the inverse chain matrix `B` and return
the result.

# Examples
```jldoctest
julia> Z=[50.0 50;50 50];JosephsonCircuits.ZtoB(Z)
2×2 Matrix{Float64}:
 1.0   -0.0
 0.02   1.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ZtoB

"""
    ZtoB!(B::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)

See [`ZtoB`](@ref) for description.

"""
function ZtoB!(B::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(Z,1)÷2
    range2 = size(Z,1)÷2+1:size(Z,1)

    # tmp = [0 Z12; -I Z22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d+size(Z,1)÷2,d] = -1
    end
    tmp[range1,range2] .= Z[range1, range2]
    tmp[range2,range2] .= Z[range2, range2]

    # B = [I Z11; 0 Z21]
    fill!(B,zero(eltype(B)))
    for d in range1
        B[d,d] = 1
    end

    B[range1,range2] .= Z[range1, range1]
    B[range2,range2] .= Z[range2, range1]

    # println(tmp)
    # println(B)
    # perform the left division
    # B = inv(tmp)*B = [0 Z12; -I Z22] \ [I Z11; 0 Z21]
    B .= tmp \ B

    return nothing
end

@doc """
    YtoS(Y;portimpedances=50.0)

Convert the admittance parameter matrix `Y` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances.

# Examples
```jldoctest
julia> Y = [1/50.0 0.0;0.0 1/50.0];JosephsonCircuits.YtoS(Y)
2×2 Matrix{Float64}:
  0.0  -0.0
 -0.0   0.0

julia> Y = [0.0 0.0;0.0 0.0];JosephsonCircuits.YtoS(Y)
2×2 Matrix{Float64}:
  1.0  -0.0
 -0.0   1.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" YtoS

"""
    YtoS!(S::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

See [`YtoS`](@ref) for description.

"""
function YtoS!(S::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)
    
    # compute \tilde{Y} = sqrtportimpedances*Y*sqrtportimpedances in S and tmp
    copy!(S,Y)
    rmul!(S,sqrtportimpedances)
    lmul!(sqrtportimpedances,S)
    copy!(tmp,S)
    rmul!(S,-1)

    # tmp = (\tilde{Y} + I)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # S = (-\tilde{Y} + I)
    for d in 1:size(S,1)
        S[d,d] += 1
    end

    # perform the left division
    # S = inv(tmp)*S = (I + \tilde{Y}) \ (I - \tilde{Y})
    S .= tmp \ S

    return nothing
end

@doc """
    YtoA(Y)

Convert the admittance matrix `Y` to the chain (ABCD) matrix `A` and return
the result.

# Examples
```jldoctest
julia> Y=[1/50 1/50;1/50 1/50];JosephsonCircuits.YtoA(Y)
2×2 Matrix{Float64}:
 -1.0  -50.0
  0.0   -1.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of overall sign on (suspected typo).
""" YtoA

"""
    YtoA!(A::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)

See [`YtoA`](@ref) for description.

"""
function YtoA!(A::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [-Y11 I; Y21 0]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d+size(A,1)÷2] = 1
    end
    tmp[range1,range1] .= -Y[range1, range1]
    tmp[range2,range1] .= Y[range2, range1]

    # A = [Y12 0; -Y22 -I]
    fill!(A,zero(eltype(A)))
    for d in range2
        A[d,d] = -1
    end

    A[range1,range1] .= Y[range1, range2]
    A[range2,range1] .= -Y[range2, range2]

    # perform the left division
    # A = inv(tmp)*A = [-Y11 I; Y21 0] \ [Y12 0; -Y22 -I]
    A .= tmp \ A

    return nothing
end

@doc """
    YtoB(Y)

Convert the admittance matrix `Y` to the inverse chain matrix `B` and return
the result.

# Examples
```jldoctest
julia> Y=[1/50 1/50;1/50 1/50];JosephsonCircuits.YtoB(Y)
2×2 Matrix{Float64}:
 -1.0  -50.0
  0.0   -1.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" YtoB

"""
    YtoB!(B::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)

See [`YtoB`](@ref) for description.

"""
function YtoB!(B::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(Y,1)÷2
    range2 = size(Y,1)÷2+1:size(Y,1)

    # tmp = [-Y12 0; -Y22 I]
    fill!(tmp,zero(eltype(tmp)))
    for d in range2
        tmp[d,d] = 1
    end
    tmp[range1,range1] .= -Y[range1, range2]
    tmp[range2,range1] .= -Y[range2, range2]

    # B = [Y11 I; Y21 0]
    fill!(B,zero(eltype(B)))
    for d in range1
        B[d,d+size(B,1)÷2] = 1
    end

    B[range1,range1] .= Y[range1, range1]
    B[range2,range1] .= Y[range2, range1]

    # perform the left division
    # B = inv(tmp)*B = [-Y12 0; -Y22 I] \ [Y11 I; Y21 0]
    B .= tmp \ B

    return nothing
end


@doc """
    AtoS(A)

Convert the chain (ABCD) matrix `A` to the scattering parameter matrix `S` and
return the result.

# Examples
```jldoctest
julia> S = rand(Complex{Float64},2,2);isapprox(S,JosephsonCircuits.AtoS(JosephsonCircuits.StoA(S)))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoS

"""
    AtoS!(S::AbstractMatrix, A::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`AtoS`](@ref) for description.

"""
function AtoS!(S::AbstractMatrix, A::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    A11 = view(A,range1,range1)
    A12 = view(A,range1,range2)
    A21 = view(A,range2,range1)
    A22 = view(A,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # tmp = [-g1 A11*g2+A12/g2; 1/g1 A21*g2+A22/g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    fill!(tmp,zero(eltype(tmp)))

    for d in range1
        tmp[d,d] = 1
    end

    for d in range1
        tmp[d+size(A,1)÷2,d] = 1
    end

    tmp11 .= -tmp11*sqrtportimpedances1
    tmp12 .= A11*sqrtportimpedances2+A12/sqrtportimpedances2
    tmp21 .= tmp21/sqrtportimpedances1
    tmp22 .= A21*sqrtportimpedances2+A22/sqrtportimpedances2

    # S = [g1 -A11*g2+A12/g2; 1/g1 -A21*g2+A22/g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    fill!(S,zero(eltype(S)))

    for d in range1
        S[d,d] = 1
    end

    for d in range1
        S[d+size(A,1)÷2,d] = 1
    end

    S11 .= S11*sqrtportimpedances1
    S12 .= -A11*sqrtportimpedances2+A12/sqrtportimpedances2
    S21 .=  S21/sqrtportimpedances1
    S22 .= -A21*sqrtportimpedances2+A22/sqrtportimpedances2

    # perform the left division
    # S = inv(tmp)*S
    S .= tmp \ S

    return nothing
end

@doc """
    AtoZ(A)

Convert the ABCD matrix `A` to the impedance matrix `Z` and return the result.

# Examples
```jldoctest
julia> A=[1.0 0.0;1/50 1.0];JosephsonCircuits.AtoZ(A)
2×2 Matrix{Float64}:
 50.0  50.0
 50.0  50.0

julia> @variables A B C D;JosephsonCircuits.AtoZ([A B;C D])
2×2 Matrix{Num}:
 A / C  -B + (A*D) / C
 1 / C           D / C
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoZ

"""
    AtoZ!(Z::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

See [`AtoZ`](@ref) for description.

"""
function AtoZ!(Z::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [-I A11; 0 A21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp[range1,range2] .= A[range1, range1]
    tmp[range2,range2] .= A[range2, range1]

    # Z = [0 A12; I A22]
    fill!(Z,zero(eltype(A)))
    for d in range1
        Z[d+size(A,1)÷2,d] = 1
    end

    Z[range1,range2] .= A[range1, range2]
    Z[range2,range2] .= A[range2, range2]

    # perform the left division
    # Z = inv(tmp)*Z = [-I A11; 0 A21] \ [0 A12; I A22]
    Z .= tmp \ Z

    return nothing

end

@doc """
    AtoY(A)

Convert the chain (ABCD) matrix `A` to the admittance matrix `Y` and return
the result.

# Examples
```jldoctest
julia> A=[-1 -50.0;0 -1];JosephsonCircuits.AtoY(A)
2×2 Matrix{Float64}:
 0.02  0.02
 0.02  0.02
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoY

"""
    AtoY!(Y::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

See [`AtoY`](@ref) for description.

"""
function AtoY!(Y::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [0 A12; I A22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d+size(A,1)÷2,d] = 1
    end
    tmp[range1,range2] .= A[range1, range2]
    tmp[range2,range2] .= A[range2, range2]

    # Y = [-I A11; 0 A21]
    fill!(Y,zero(eltype(A)))
    for d in range1
        Y[d,d] = -1
    end

    Y[range1,range2] .= A[range1, range1]
    Y[range2,range2] .= A[range2, range1]

    # perform the left division
    # Y = inv(tmp)*Z = [0 A12; I A22] \ [-I A11; 0 A21]
    Y .= tmp \ Y

    return nothing

end

@doc """
    AtoB(A)

Convert the chain (ABCD) matrix `A` to the inverse chain matrix `B` and return
the result. Note that despite the name, the inverse of the chain
matrix is not equal to the inverse chain matrix, inv(A) ≠ B.

# Examples
```jldoctest
julia> A=[1.0 0.0;1/50 1.0];JosephsonCircuits.AtoB(A)
2×2 Matrix{Float64}:
 1.0   0.0
 0.02  1.0

julia> S = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.AtoB(JosephsonCircuits.ZtoA(JosephsonCircuits.StoZ(S))),JosephsonCircuits.StoB(S))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoB

"""
    AtoB!(B::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

See [`AtoB`](@ref) for description.

"""
function AtoB!(B::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [A11 -A12; A21 -A22]
    copy!(tmp,A)
    tmp[range1,range2] .*= -1
    tmp[range2,range2] .*= -1

    # B = [I 0; 0 I]
    fill!(B,zero(eltype(A)))
    for d in range1
        B[d,d] = 1
    end
    for d in range2
        B[d,d] = -1
    end

    # perform the left division
    # B = inv(tmp)*B = [A11 -A12; A21 -A22] \ [I 0; 0 I]
    B .= tmp \ B

    return nothing

end

@doc """
    BtoS(B)

Convert the inverse chain (ABCD) matrix `B` to the scattering parameter matrix
`S` and return the result.

# Examples
```jldoctest
julia> S = rand(Complex{Float64},2,2);isapprox(S,JosephsonCircuits.BtoS(JosephsonCircuits.StoB(S)))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of overall sign (suspected typo).
""" BtoS

"""
    BtoS!(S::AbstractMatrix, B::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`BtoS`](@ref) for description.

"""
function BtoS!(S::AbstractMatrix, B::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    B11 = view(B,range1,range1)
    B12 = view(B,range1,range2)
    B21 = view(B,range2,range1)
    B22 = view(B,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # tmp = [B11*g1+B12/g1 -g2; B21*g1+B22/g1 1/g2]
    fill!(tmp,zero(eltype(tmp)))

    for d in range1
        tmp[d,d+size(B,1)÷2] = 1
    end

    for d in range2
        tmp[d,d] = 1
    end

    tmp11 .= B11*sqrtportimpedances1+B12/sqrtportimpedances1
    tmp12 .= -tmp12*sqrtportimpedances2
    tmp21 .= B21*sqrtportimpedances1+B22/sqrtportimpedances1
    tmp22 .= tmp22/sqrtportimpedances2

    # S = [-B11*g1+B12/g1 g2; -B21*g1+B22/g1 1/g2]
    fill!(S,zero(eltype(S)))

    for d in range1
        S[d,d+size(B,1)÷2] = 1
    end

    for d in range2
        S[d,d] = 1
    end

    S11 .= -B11*sqrtportimpedances1+B12/sqrtportimpedances1
    S12 .= S12*sqrtportimpedances2
    S21 .= -B21*sqrtportimpedances1+B22/sqrtportimpedances1
    S22 .= S22/sqrtportimpedances2

    # perform the left division
    # S = inv(tmp)*S
    S .= tmp \ S

    return nothing
end

@doc """
    BtoZ(A)

Convert the inverse chain matrix `B` to the impedance matrix `Z` and return
the result.

# Examples
```jldoctest
julia> B=[1.0 0.0;1/50 1];JosephsonCircuits.BtoZ(B)
2×2 Matrix{Float64}:
 50.0  50.0
 50.0  50.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of sign on B21 and B22 terms (suspected typo).
""" BtoZ

"""
    BtoZ!(Z::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)

See [`BtoZ`](@ref) for description.

"""
function BtoZ!(Z::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # tmp = [B11 -I; B21 0]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d+size(B,1)÷2] = -1
    end
    tmp[range1,range1] .= B[range1, range1]
    tmp[range2,range1] .= -B[range2, range1]

    # Z = [B12 0; B22 -I]
    fill!(Z,zero(eltype(B)))
    for d in range2
        Z[d,d] = -1
    end

    Z[range1,range1] .= B[range1, range2]
    Z[range2,range1] .= -B[range2, range2]

    # perform the left division
    # Y = inv(tmp)*Y = [B12 0; B22 I] \ [B11 -I; B21 0]
    # Z = inv(tmp)*Z = [B11 -I; B21 0] \ [B12 0; B22 -I]
    Z .= tmp \ Z

    return nothing
end

@doc """
    BtoY(A)

Convert the inverse chain matrix `B` to the admittance matrix `Y` and return
the result.

# Examples
```jldoctest
julia> @variables A B C D;JosephsonCircuits.BtoY([A B;C D])
2×2 Matrix{Num}:
          A / B       -1 / B
 C + (-A*D) / B  -((-D) / B)

julia> B=[-1 -50.0;0 -1];JosephsonCircuits.BtoY(B)
2×2 Matrix{Float64}:
 0.02  0.02
 0.02  0.02
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" BtoY

"""
    BtoY!(Y::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)

See [`BtoY`](@ref) for description.

"""
function BtoY!(Y::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # tmp = [B12 0; B22 I]
    fill!(tmp,zero(eltype(tmp)))
    for d in range2
        tmp[d,d] = 1
    end
    tmp[range1,range1] .= B[range1, range2]
    tmp[range2,range1] .= B[range2, range2]

    # Y = [B11 -I; B21 0]
    fill!(Y,zero(eltype(B)))
    for d in range1
        Y[d,d+size(B,1)÷2] = -1
    end

    Y[range1,range1] .= B[range1, range1]
    Y[range2,range1] .= B[range2, range1]

    # perform the left division
    # Y = inv(tmp)*Y = [B12 0; B22 I] \ [B11 -I; B21 0]
    Y .= tmp \ Y

    return nothing
end

@doc """
    BtoA(B)

Convert the inverse chain matrix `B` to the chain (ABCD) matrix `A` and return
the result. Note that despite the name, the inverse of the chain
matrix is not equal to the inverse chain matrix, inv(A) ≠ B.

# Examples
```jldoctest
julia> B=[1.0 0.0;1/50 1.0];JosephsonCircuits.BtoA(B)
2×2 Matrix{Float64}:
 1.0   0.0
 0.02  1.0
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" BtoA


"""
    BtoA!(A::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)

See [`BtoA`](@ref) for description.

"""
function BtoA!(A::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)
    return AtoB!(A,B,tmp)
end

@doc """
    ABCDtoS(ABCD;portimpedances=50.0)

Convert the 2 port chain (ABCD) matrix `ABCD` to the scattering parameter
matrix `S` and return the result. Assumes a port impedance of 50 Ohms unless
specified with the `portimpedances` keyword argument.

# Examples
```jldoctest
julia> A = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.AtoS(A),JosephsonCircuits.ABCDtoS(A))
true

julia> A = rand(Complex{Float64},2,2,10);isapprox(JosephsonCircuits.AtoS(A),JosephsonCircuits.ABCDtoS(A))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ABCDtoS

function ABCDtoS!(S,ABCD,portimpedances)
    return ABCDtoS!(S,ABCD,first(portimpedances),last(portimpedances))
end

function ABCDtoS!(S,A,RS,RL)
    S[1,1] = (A[1,1]*RL+A[1,2]-A[2,1]*RS*RL-A[2,2]*RS)/(A[1,1]*RL+A[1,2]+A[2,1]*RS*RL+A[2,2]*RS)
    S[1,2] = 2*sqrt(RS*RL)*(A[1,1]*A[2,2]-A[1,2]*A[2,1])/(A[1,1]*RL+A[1,2]+A[2,1]*RS*RL+A[2,2]*RS)
    S[2,1] = 2*sqrt(RS*RL)/(A[1,1]*RL+A[1,2]+A[2,1]*RS*RL+A[2,2]*RS)
    S[2,2] = (-A[1,1]*RL+A[1,2]-A[2,1]*RS*RL+A[2,2]*RS)/(A[1,1]*RL+A[1,2]+A[2,1]*RS*RL+A[2,2]*RS)
    return nothing
end

@doc """
    StoABCD(S;portimpedances=50.0))

Convert the scattering parameter matrix `S` to the 2 port chain (ABCD) matrix and
return the result. Assumes a port impedance of 50 Ohms unless specified with the
`portimpedances` keyword argument.

# Examples
```jldoctest
julia> S = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.StoA(S),JosephsonCircuits.StoABCD(S))
true

julia> S = rand(Complex{Float64},2,2,10);isapprox(JosephsonCircuits.StoA(S),JosephsonCircuits.StoABCD(S))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoABCD

function StoABCD!(ABCD,S,portimpedances)
    return StoABCD!(ABCD,S,first(portimpedances),last(portimpedances))
end
function StoABCD!(ABCD,S,RS,RL)
    ABCD[1,1] = sqrt(RS/RL)*((1+S[1,1])*(1-S[2,2])+S[2,1]*S[1,2])/(2*S[2,1])
    ABCD[1,2] = sqrt(RS*RL)*((1+S[1,1])*(1+S[2,2])-S[2,1]*S[1,2])/(2*S[2,1])
    ABCD[2,1] = 1/sqrt(RS*RL)*((1-S[1,1])*(1-S[2,2])-S[2,1]*S[1,2])/(2*S[2,1])
    ABCD[2,2] = sqrt(RL/RS)*((1-S[1,1])*(1+S[2,2])+S[2,1]*S[1,2])/(2*S[2,1])
    return nothing
end


"""
    ABCD_seriesZ(Z)

Return the ABCD matrix for a series impedance `Z`.
```
o---Z1---o
          
          
o--------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_seriesZ(50)
2×2 Matrix{Int64}:
 1  50
 0   1
```
"""
function ABCD_seriesZ(Z)
    return [1 Z;0 1]
end

"""
    ABCD_shuntY(Y)

Return the ABCD matrix for a shunt admittance `Y`.
```
o---------o
     |
     Y
     |
o---------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_shuntY(1/50)
2×2 Matrix{Float64}:
 1.0   0.0
 0.02  1.0
```
"""
function ABCD_shuntY(Y)
    return [1 0;Y 1]
end

"""
    ABCD_PiY(Y1,Y2,Y3)

Return the ABCD matrix for a Pi network of admittances `Y1`, `Y2`, and `Y3`.
```
o----Y3-----o
   |     |   
   Y1    Y2  
   |     |   
o-----------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_PiY(1,2,4)
2×2 Matrix{Float64}:
 1.5  0.25
 3.5  1.25
```
"""
function ABCD_PiY(Y1,Y2,Y3)
    return [1+Y2/Y3 1/Y3;Y1+Y2+Y1*Y2/Y3 1+Y1/Y3]
end

"""
    Y_PiY(Y1,Y2,Y3)

Return the admittance matrix for a Pi network of admittances `Y1`, `Y2`, and
`Y3`.
```
o----Y3-----o
   |     |   
   Y1    Y2  
   |     |   
o-----------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.Y_PiY(1.0,2.0,4.0)
2×2 Matrix{Float64}:
  5.0  -4.0
 -4.0   6.0

julia> Y1=1.0;Y2=2.0;Y3=4.0;isapprox(JosephsonCircuits.YtoA(JosephsonCircuits.Y_PiY(Y1,Y2,Y3)),JosephsonCircuits.ABCD_PiY(Y1,Y2,Y3))
true
```
"""
function Y_PiY(Y1,Y2,Y3)
    return [Y1+Y3 -Y3;-Y3 Y2+Y3]
end

"""
    ABCD_TZ(Z1,Z2,Z3)

Return the ABCD matrix for a T network of impedances `Z1`, `Z2`, and `Z3`.
```
o--Z1-----Z2--o
       |       
      Z3       
       |       
o-------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_TZ(1,2,4)
2×2 Matrix{Float64}:
 1.25  3.5
 0.25  1.5
```
"""
function ABCD_TZ(Z1,Z2,Z3)
    return [1+Z1/Z3 Z1+Z2+Z1*Z2/Z3;1/Z3 1+Z2/Z3]
end

"""
    Z_TZ(Z1,Z2,Z3)

Return the ABCD matrix for a T network of impedances `Z1`, `Z2`, and `Z3`.
```
o--Z1-----Z2--o
       |       
      Z3       
       |       
o-------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_TZ(1,2,4)
2×2 Matrix{Float64}:
 1.25  3.5
 0.25  1.5

julia> Z1=1.0;Z2=2.0;Z3=4.0;isapprox(JosephsonCircuits.ZtoA(JosephsonCircuits.Z_TZ(Z1,Z2,Z3)),JosephsonCircuits.ABCD_TZ(Z1,Z2,Z3))
true
```
"""
function Z_TZ(Z1,Z2,Z3)
    return [Z1+Z3 Z3;Z3 Z2+Z3]
end

"""
    ABCD_tline(theta,Z0)

Return the ABCD matrix for a transmission line described by phase delay
`theta` in radians and characteristic impedance `Z0` in Ohms.
```
   theta, Z0  
o--========--o
              
              
o------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_tline(pi/4,1)
2×2 Matrix{ComplexF64}:
 0.707107+0.0im            0.0+0.707107im
      0.0+0.707107im  0.707107+0.0im
```
"""
function ABCD_tline(theta,Z0)
    return [cos(theta) im*Z0*sin(theta);im/Z0*sin(theta) cos(theta)]
end

"""
    Z_tline(theta, Z0)

Return the impedance matrix for a transmission line described by phase delay
`theta` in radians and characteristic impedance `Z0` in Ohms.
```
   theta, Z0  
o--========--o
              
              
o------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.Z_tline(pi/4,1)
2×2 Matrix{ComplexF64}:
 0.0-1.0im      0.0-1.41421im
 0.0-1.41421im  0.0-1.0im

julia> isapprox(JosephsonCircuits.Z_tline(pi/4,1),JosephsonCircuits.AtoZ(JosephsonCircuits.ABCD_tline(pi/4,1)))
true
```
"""
function Z_tline(theta, Z0)
    return [-im*Z0*cot(theta) -im*Z0*csc(theta);-im*Z0*csc(theta) -im*Z0*cot(theta)]
end

"""
    ABCD_coupled_tline(thetae, thetao, Z0e, Z0o)

Return the ABCD matrix for two coupled transmission lines described by
even and odd mode phase delays `thetae` and `thetao`, and even and odd mode
impedances `Z0e` and `Z0o`.
```
thetae, Z0e
thetao, Z0o

V1, I1 -->  ======== <-- I3, V3
V2, I2 -->  ======== <-- I4, V4

[(V1+V2)/2, (I1+I2)/2] = ABCDe * [(V3+V4)/2, -(I3+I4)/2]
[(V1-V2)/2, (I1-I2)/2] = ABCDo * [(V3-V4)/2, -(I3-I4)/2]

[V1, V2, I1, I2] = ABCD_coupled_tline * [V3, V4, -I3, -I4]
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_coupled_tline(pi/4,pi/4,50,50)
4×4 Matrix{ComplexF64}:
 0.707107+0.0im             0.0+0.0im        …       0.0+0.0im
      0.0+0.0im        0.707107+0.0im                0.0+35.3553im
      0.0+0.0141421im       0.0+0.0im                0.0+0.0im
      0.0+0.0im             0.0+0.0141421im     0.707107+0.0im

julia> isapprox(JosephsonCircuits.AtoZ(JosephsonCircuits.ABCD_coupled_tline(pi/4,pi/4,55,50)),JosephsonCircuits.Z_coupled_tline(pi/4,pi/4,55,50))
true
```
"""
function ABCD_coupled_tline(thetae, thetao, Z0e, Z0o)
    A11 = A22 = A33 = A44 = 1/2*(cos(thetae) + cos(thetao))
    A12 = A21 = A34 = A43 = 1/2*(cos(thetae) - cos(thetao))
    A13 = A24 = im/2*(Z0e*sin(thetae) + Z0o*sin(thetao))
    A31 = A42 = im/2*(sin(thetae)/Z0e + sin(thetao)/Z0o)
    A14 = A23 = im/2*(Z0e*sin(thetae) - Z0o*sin(thetao))
    A41 = A32 = im/2*(sin(thetae)/Z0e - sin(thetao)/Z0o)
    return [A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44]
end

"""
    Z_coupled_tline(thetae, thetao, Z0e, Z0o)

Return the impedance matrix for two coupled transmission lines described by
even and odd mode phase delays `thetae` and `thetao`, and even and odd mode
impedances `Z0e` and `Z0o`.
```
thetae, Z0e
thetao, Z0o

V1, I1 -->  ======== <-- I3, V3
V2, I2 -->  ======== <-- I4, V4

[(V1+V2)/2, (I1+I2)/2] = ABCDe * [(V3+V4)/2, -(I3+I4)/2]
[(V1-V2)/2, (I1-I2)/2] = ABCDo * [(V3-V4)/2, -(I3-I4)/2]

[V1, V2, V3, V4] = Z_coupled_tline * [I1, I2, I3, I4]
```
# Examples
```jldoctest
julia> JosephsonCircuits.Z_coupled_tline(pi/4,pi/4,50,50)
4×4 Matrix{ComplexF64}:
 0.0-50.0im     0.0-0.0im      0.0-70.7107im  0.0-0.0im
 0.0-0.0im      0.0-50.0im     0.0-0.0im      0.0-70.7107im
 0.0-70.7107im  0.0-0.0im      0.0-50.0im     0.0-0.0im
 0.0-0.0im      0.0-70.7107im  0.0-0.0im      0.0-50.0im
```
"""
function Z_coupled_tline(thetae,thetao,Z0e,Z0o)
    Z11 = Z22 = Z33 = Z44 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
    Z12 = Z21 = Z34 = Z43 = -im/2*(Z0e*cot(thetae) - Z0o*cot(thetao)) 
    Z13 = Z31 = Z42 = Z24 = -im/2*(Z0e*csc(thetae) + Z0o*csc(thetao))
    Z14 = Z41 = Z32 = Z23 = -im/2*(Z0e*csc(thetae) - Z0o*csc(thetao))
    return [Z11 Z12 Z13 Z14; Z21 Z22 Z23 Z24; Z31 Z32 Z33 Z34; Z41 Z42 Z43 Z44]
end

"""
    A_coupled_tlines(L,Cmaxwell,omega,l)

Returns the 2mx2m chain (ABCD) matrix for a port number symmetric multi-port
network of m coupled transmission lines described by a symmetric mxm
Maxwell inductance (per unit length) matrix `L`, a symmetric mxm Maxwell
capacitance matrix (per unit length) `Cmaxwell`, an angular frequency `omega`,
and a physical length `l`.
```
V_1, I_1 -->  ======== <-- I_{m+1}, V_{m+1}
V_2, I_2 -->  ======== <-- I_{m+2}, V_{m+2}
          .
          .
          .
V_m, I_m -->  ======== <-- I_n, V_n
              <---l-->

where n=2*m.

[V_1, ...V_m, I_1, ...I_m] = A_coupled_tline * [V_{m+1}, ...V_n, I_{m+1}, ...I_n]
```

# Examples
```jldoctest
Zeven = 51.0
Zodd = 49.0
neven = 1.1
nodd = 1.08
l = 3.5e-3
c = JosephsonCircuits.speed_of_light
omega = 2*pi*5e9

L, C = JosephsonCircuits.even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
A1 = JosephsonCircuits.A_coupled_tlines(L,C,omega,l)
A2 = JosephsonCircuits.ABCD_coupled_tline(neven*omega/c*l,nodd*omega/c*l,Zeven,Zodd)
println(isapprox(A1,A2))

# output
true
```

# References
Paul, Clayton R. Analysis of Multiconductor Transmission Lines, Second
Edition. Wiley, 2008.
"""
function A_coupled_tlines(L,Cmaxwell,omega,l)

    N = size(L,1)

    ZC, TI, TV, theta, U, lambda, S = ZC_basis_coupled_tlines(L,Cmaxwell)
    TIinv = inv(TI)
    TVinv = inv(TV)
    YC = inv(ZC)

    gamma = im*lambda*omega

    phi11 = 1/2*ZC*TI*(exp(gamma*l)+exp(-gamma*l))*TIinv*YC
    phi12 = -1/2*ZC*TI*(exp(-gamma*l)-exp(gamma*l))*TIinv
    phi21 = -1/2*TI*(exp(-gamma*l)-exp(gamma*l))*TIinv*YC
    phi22 = 1/2*TI*(exp(gamma*l)+exp(-gamma*l))*TIinv

    # phi11 = 1/2*TV*(exp(gamma*l)+exp(-gamma*l))*TVinv
    # phi12 = -1/2*TV*(exp(gamma*l)-exp(-gamma*l))*TVinv*ZC
    # phi21 = -1/2*YC*TV*(exp(gamma*l)-exp(-gamma*l))*TVinv
    # phi22 = 1/2*YC*TV*(exp(gamma*l)+exp(-gamma*l))*TVinv*ZC

    A = zeros(Complex{Float64},2*N,2*N)
    A[1:N,1:N] .= phi11
    A[1:N,N+1:end] .= phi12
    A[N+1:end,1:N] .= phi21
    A[N+1:end,N+1:end] .= phi22
    return A
end

"""
    ZC_basis_coupled_tlines(L, Cmaxwell)

    Returns the characteristic impedance matrix `ZC` and eigenbasis for
    current `TI` and voltage `TV` from the inductance per unit length matrix
    `L` and Maxwell capacitance per unit length matrix `Cmaxwell`.

# Arguments
- `L`: inductance per unit length matrix.
- `C`: Maxwell capacitance per unit length matrix.

# Returns
- `ZC`: characteristic impedance matrix.
- `TI`: matrix which transforms mode currents to currents, I = TI*Im. Computed
        from TI = U*theta*S.
- `TV`: matrix which transforms mode voltages to voltages, V = TV*Vm. Computed
        from TV = U*inv(theta)*S.
- `theta`: Diagonal matrix with the square of the eigenvalues of Cmaxwell
        along the diagonals.
- `U`: eigenvectors of Cmaxwell.
- `lambda`: Diagonal matrix with the square of the eigenvalues of
        theta*Ut*L*U*theta along the diagonals. `lambda` is related
        to the propagation constant, gamma, as gamma^2 = -omega^2*lambda^2.
- `S`: eigenvectors of theta*Ut*L*U*theta.

# Examples
```jldoctest
Zeven = 51.0
Zodd = 49.0
neven = 1.1
nodd = 1.08
c = JosephsonCircuits.speed_of_light

L, C = JosephsonCircuits.even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
ZC, TI, TV, theta, U, lambda, S = JosephsonCircuits.ZC_basis_coupled_tlines(L,C)
@show ZC
@show TI
@show TV
@show Matrix(theta)
@show U
@show Matrix(lambda)
@show S
println(isapprox(Zeven,ZC[1,1]+ZC[1,2]))
println(isapprox(Zodd,ZC[1,1]-ZC[1,2]))
println(isapprox(neven,lambda[2,2]*c))
println(isapprox(nodd,lambda[1,1]*c))

# output
ZC = [49.999999999999986 0.9999999999999929; 0.9999999999999929 49.999999999999986]
TI = [-6.063012846509498e-6 -5.997716107132906e-6; 6.063012846509498e-6 -5.997716107132906e-6]
TV = [-82467.25063230767 -83365.06614665617; 82467.25063230767 -83365.06614665617]
Matrix(theta) = [8.48205146197092e-6 0.0; 0.0 8.574394996376037e-6]
U = [-0.7071067811865475 -0.7071067811865475; -0.7071067811865475 0.7071067811865475]
Matrix(lambda) = [3.6024922281400425e-9 0.0; 0.0 3.669205047179673e-9]
S = [0.0 1.0; 1.0 0.0]
true
true
true
true
```

# References
Paul, Clayton R. Analysis of Multiconductor Transmission Lines, Second
Edition. Wiley, 2008.
"""
function ZC_basis_coupled_tlines(L,Cmaxwell)

    # compute U and theta
    vals, vecs = eigen(Cmaxwell)
    U = vecs
    theta = Diagonal(sqrt.(vals))

    # compute S and lambda
    vals, vecs = eigen(theta*U'*L*U*theta)
    lambda = Diagonal(sqrt.(vals))
    S = vecs

    # compute TI and TV
    TI = U*theta*S
    TV = U*inv(theta)*S

    # compute ZC and YC
    ZC = U*inv(theta)*S*lambda*S'*inv(theta)*U'

    return (ZC = ZC, TI = TI, TV = TV, theta = theta, U = U, lambda = lambda,
        S = S)
end


"""
    maxwell_to_mutual(Cmaxwell::AbstractMatrix)

Return the mutual capacitance matrix from the Maxwell capacitance matrix
`Cmaxwell`.

The Maxwell capacitance `Cmaxwell` is the relationship between charge and
voltage on each node, Q = C V or dQi/dVj = C_ij where `C` is the
Maxwell capacitance matrix.

Each element of the mutual capacitance matrix `Cmutual` is the value of a
physical capacitor placed between two nodes in a circuit or between a node
and ground.

# Examples
```jldoctest
julia> @variables C11, C12, C21, C22;C = [C11 C12;C21 C22];JosephsonCircuits.maxwell_to_mutual(C)
2×2 Matrix{Num}:
 C11 + C12       -C12
      -C21  C21 + C22

julia> C = [1.0 -0.1;-0.1 2.0];JosephsonCircuits.maxwell_to_mutual(C)
2×2 Matrix{Float64}:
 0.9  0.1
 0.1  1.9
```
"""
function maxwell_to_mutual(Cmaxwell::AbstractMatrix)
  return Diagonal(Cmaxwell)-Cmaxwell + Diagonal(sum(Cmaxwell,dims=2)[:])
end

"""
    maxwell_to_mutual(Cmutual::AbstractMatrix)

Return the Maxwell capacitance matrix from the mutual capacitance matrix
`Cmutual`.

The Maxwell capacitance `Cmaxwell` is the relationship between charge and
voltage on each node, Q = C V or dQi/dVj = C_ij where `C` is the
Maxwell capacitance matrix.

Each element of the mutual capacitance matrix `Cmutual` is the value of a
physical capacitor placed between two nodes in a circuit or between a node
and ground.

# Examples
```jldoctest
julia> @variables Cg, Cm;C = [Cg Cm;Cm Cg];JosephsonCircuits.mutual_to_maxwell(C)
2×2 Matrix{Num}:
 Cg + Cm      -Cm
     -Cm  Cg + Cm

julia> C = [0.9 0.1;0.1 1.9];JosephsonCircuits.mutual_to_maxwell(C)
2×2 Matrix{Float64}:
  1.0  -0.1
 -0.1   2.0
```
"""
function mutual_to_maxwell(C::AbstractMatrix)
  return maxwell_to_mutual(C)
end

"""
    maxwell_to_even_odd(L, Cmaxwell)

Return the even and odd mode impedances and the even and odd mode indices from
the inductance matrix `L` and the Maxwell capacitance matrix `Cmaxwell`.

# Examples
```jldoctest
@variables C11, C12, L11, L12
C = [C11 C12;C12 C11]
L = [L11 L12;L12 L11]
Zeven, Zodd, neven, nodd = JosephsonCircuits.maxwell_to_even_odd(L,C)
@show Zeven
@show Zodd
@show neven
@show nodd
;

# output
Zeven = sqrt((L11 + L12) / (C11 + C12))
Zodd = sqrt((L11 - L12) / (C11 - C12))
neven = 2.99792458e8sqrt((C11 + C12)*(L11 + L12))
nodd = 2.99792458e8sqrt((C11 - C12)*(L11 - L12))
```
"""
function maxwell_to_even_odd(L, Cmaxwell)

    # consider erroring if not 2x2 matrices with capacitance and inductance
    # to ground equal for both transmission lines.
    # also error if not symmetric
    # (Cmaxwell[1,1] != Cmaxwell[2,2]) || (Cmaxwell[1,2] != Cmaxwell[2,1])
    # size(Cmaxwell) != (2,2)
    # same for inductance matrix

    c = JosephsonCircuits.speed_of_light
    Zeven = sqrt((L[1,1]+L[1,2])/(Cmaxwell[1,1]+Cmaxwell[1,2]))
    Zodd = sqrt((L[1,1]-L[1,2])/(Cmaxwell[1,1]-Cmaxwell[1,2]))
    neven = c*sqrt((L[1,1]+L[1,2])*(Cmaxwell[1,1]+Cmaxwell[1,2]))
    nodd = c*sqrt((L[1,1]-L[1,2])*(Cmaxwell[1,1]-Cmaxwell[1,2]))
    return (Zeven = Zeven, Zodd = Zodd, neven = neven, nodd = nodd)
end

"""
    mutual_to_even_odd(L, Cmutual)

Return the even and odd mode impedances and the even and odd mode indices from
the inductance matrix `L` and the mutual capacitance matrix `Cmutual`.

# Examples
```jldoctest
@variables Cg, Cm, Ls, Lm
C = [Cg Cm; Cm Cg]
L = [Ls Lm; Lm Ls]
Zeven, Zodd, neven, nodd = JosephsonCircuits.mutual_to_even_odd(L,C)
@show Zeven
@show Zodd
@show neven
@show nodd
;

# output
Zeven = sqrt((Lm + Ls) / Cg)
Zodd = sqrt((-Lm + Ls) / (Cg + 2Cm))
neven = 2.99792458e8sqrt(Cg*(Lm + Ls))
nodd = 2.99792458e8sqrt((Cg + 2Cm)*(-Lm + Ls))
```
"""
function mutual_to_even_odd(L, Cmutual)
    return maxwell_to_even_odd(L, maxwell_to_mutual(Cmutual))
end


"""
    even_odd_to_maxwell(Zeven, Zodd, neven, nodd)

Return the inductance matrix and Maxwell capacitance matrix for two coupled
transmission lines with even and odd mode impedances `Zeven`, `Zodd` and
even and odd mode indices `neven`, `nodd`.

# Examples
```jldoctest
L1 = [1.1 0.1;0.1 1.1]
C1 = [2.0 -0.4;-0.4 2.0]
L2, C2 = JosephsonCircuits.even_odd_to_maxwell(JosephsonCircuits.maxwell_to_even_odd(L1,C1)...)
isapprox(L1,L2) && isapprox(C1,C2)

# output
true
```
"""
function even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
    c = JosephsonCircuits.speed_of_light
    L11 = (neven*Zeven+nodd*Zodd)/(2*c)
    L12 = (neven*Zeven-nodd*Zodd)/(2*c)
    C11 = (neven*Zodd+nodd*Zeven)/(2*c*Zeven*Zodd)
    C12 = (neven*Zodd-nodd*Zeven)/(2*c*Zodd*Zeven)
    return (L = [L11 L12;L12 L11], C = [C11 C12;C12 C11])
end

"""
    even_odd_to_mutual(Zeven, Zodd, neven, nodd)

Return the inductance matrix and mutual capacitance matrix for two coupled
transmission lines with even and odd mode impedances `Zeven`, `Zodd` and
even and odd mode indices `neven`, `nodd`.

# Examples
```jldoctest
L1 = [1.1 0.1;0.1 1.1]
C1 = [1.6 0.4;0.4 1.6]
L2, C2 = JosephsonCircuits.even_odd_to_mutual(JosephsonCircuits.mutual_to_even_odd(L1,C1)...)
isapprox(L1,L2) && isapprox(C1,C2)

# output
true
```
"""
function even_odd_to_mutual(Zeven, Zodd, neven, nodd)
    L, Cmaxwell = even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
    return (L = L, C = maxwell_to_mutual(Cmaxwell))
end


"""
    Z_canonical_coupled_line_circuit(i::Int, thetae, thetao, Z0e, Z0o)

Return the impedance matrix for the `i`'th canonical coupled line circuit, as
a function of the even mode phase delay `thetae` in radians, the odd mode
phase delay `thetao` in radians, the even mode characteristic impedance `Z0e`
in Ohms, and the odd mode characteristic impedance `Z0o` in Ohms.
1) low pass
```
   gnd--==========
1--> o--==========--o <--2
```
2) band pass
```
   gnd--==========--o <--2
1--> o--==========--gnd
```
3) band pass
```
        ==========--o <--2
1--> o--==========
```
4) band pass
```
1--> o--==========--gnd
2--> o--==========
```
5) all pass
```
        ==========
1--> o--==========--o  <--2
```
6) all pass
```
   gnd--==========--gnd
1--> o--==========--o  <--2
```
7) all pass
```
1--> o--==========--|
2--> o--==========--|
```
8) all stop
```
   gnd--==========--o  <--2
1--> o--==========
```
9) all stop
```
1--> o--==========--gnd
2--> o--==========--gnd
```
10) all stop
```
1--> o--==========
2--> o--==========
```
# Examples
```jldoctest
julia> @variables θe, θo, Ze, Zo;JosephsonCircuits.Z_canonical_coupled_line_circuits(3,θe,θo,Ze,Zo)
2×2 Matrix{Complex{Num}}:
 -0.5(Ze*cot(θe) + Zo*cot(θo))*im  -0.5(Ze*csc(θe) - Zo*csc(θo))*im
 -0.5(Ze*csc(θe) - Zo*csc(θo))*im  -0.5(Ze*cot(θe) + Zo*cot(θo))*im
```
# References
E. M. T. Jones, "Coupled-Strip-Transmission-Line Filters and Directional
Couplers," in IRE Transactions on Microwave Theory and Techniques, vol. 4,
no. 2, pp. 75-81, April 1956, doi: 10.1109/TMTT.1956.1125022.

Pozar, D. M. Microwave Engineering (4 ed.). John Wiley & Sons (2011)
ISBN 9780470631553.
"""
function Z_canonical_coupled_line_circuits(i::Int, thetae, thetao, Z0e, Z0o)

    if i == 1
        Z11 = -2*im*Z0e*Z0o*cos(thetae)*cos(thetao)/(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao))
        Z22 = im*(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))/(2*(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao)))
        Z12 = Z21 = -im*Z0e*Z0o*(cos(thetae)+cos(thetao))/(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao))
    elseif i == 2
        Z11 = Z22 = -2*im*Z0e*Z0o*(Z0e*cos(thetao)*sin(thetae)+Z0o*cos(thetae)*sin(thetao))/(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
        Z12 = Z21 = 2*im*Z0e*Z0o*(Z0e*sin(thetae)-Z0o*sin(thetao))/(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
    elseif i == 3
        Z11 = Z22 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
        Z12 = Z21 = -im/2*(Z0e*csc(thetae) - Z0o*csc(thetao))
    elseif i == 4
        Z11 = im*(2*Z0e*Z0o-2*Z0e*Z0o*cos(thetae)*cos(thetao)+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))/(2*(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao)))
        Z22 = im*(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))/(2*(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao)))
        Z12 = Z21 = im*(Z0e-Z0o)*(Z0e+Z0o)/(2*(Z0e*cot(thetae)+Z0o*cot(thetao)))
    elseif i == 5
        Z11 = Z22 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
        Z12 = Z21 = -im/2*(Z0e*csc(thetae) + Z0o*csc(thetao))
    elseif i == 6
        Z11 = Z22 = -2*im*Z0e*Z0o*(Z0e*cos(thetao)*sin(thetae)+Z0o*cos(thetae)*sin(thetao))/(2*Z0e*Z0o-2*Z0e*Z0o*cos(thetae)*cos(thetao)+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
        Z12 = Z21 = -2*im*Z0e*Z0o*(Z0e*sin(thetae)+Z0o*sin(thetao))/(2*Z0e*Z0o-2*Z0e*Z0o*cos(thetae)*cos(thetao)+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
    elseif i == 7
        Z11 = Z22 = -im/2*(Z0e*cot(thetae)-Z0o*tan(thetao))
        Z21 = Z12 = -im/2*(Z0e*cot(thetae)+Z0o*tan(thetao))
    elseif i == 8
        Z11 = -2*im*Z0e*Z0o*cot(thetae)*cot(thetao)/(Z0e*cot(thetae)+Z0o*cot(thetao))
        Z22 = -im*(Z0e*cot(thetae) + Z0o*cot(thetao) - Z0e*csc(thetae)-Z0o*csc(thetao))*(Z0e*cot(thetae)+Z0e*csc(thetae)+Z0o*cot(thetao)+Z0o*csc(thetao))/(2*(Z0e*cot(thetae)+Z0o*cot(thetao)))
        Z12 = Z21 = im*Z0e*Z0o*(cos(thetae)-cos(thetao))*csc(thetae)*csc(thetao)/(Z0e*cot(thetae)+Z0o*cot(thetao))
    elseif i == 9
        Z11 = Z22 = im/2*(Z0e*tan(thetae)+Z0o*tan(thetao))
        Z12 = Z21 = im/2*(Z0e*tan(thetae)-Z0o*tan(thetao))
    elseif i == 10
        Z11 = Z22 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
        Z12 = Z21 = -im/2*(Z0e*cot(thetae) - Z0o*cot(thetao))
    else
        throw(ArgumentError("Canonical coupled line circuit number must be 1-10."))
    end

    return [Z11 Z12; Z21 Z22]
end

function canonical_coupled_line_circuits(i::Int, ne, no, Z0e, Z0o)
    c = JosephsonCircuits.speed_of_light

    if i == 3
        L1 = L2 = (ne*Z0e+no*Z0o)/(6*c)
        M = -(ne*Z0e-no*Z0o)/(12*c)
        C1 = C2 = ne/(c*Z0e)
        Cm = (no*Z0e-ne*Z0o)/(2*c*Z0e*Z0o)
        return (L1 = L1, L2 = L2, M = M, C1 = C1, C2 = C2, Cm = Cm)
    elseif i == 8
        L1 = 2*Z0e*Z0o*(no^3*Z0e+ne^3*Z0o)/(3*c*(Z0e*no+Z0o*ne)^2)
        L2 = (ne*Z0e+no*Z0o)/(2*c)
        M = -Z0e*Z0o*(ne^2-no^2)/(2*c*(no*Z0e+ne*Z0o))
        C1 = (no*Z0e+ne*Z0o)/(2*c*Z0e*Z0o)
        return (L1 = L1, L2 = L2, M = M, C1 = C1)
    elseif i == 9
        L1 = L2 = (ne*Z0e+no*Z0o)/(2*c)
        M = (ne*Z0e-no*Z0o)/(2*c)
        return (L1 = L1, L2 = L2, M = M)
    elseif i == 10
        L1 = L2 = (ne*Z0e+no*Z0o)/(6*c)
        M = (ne*Z0e-no*Z0o)/(6*c)
        C1 = C2 = ne/(c*Z0e)
        Cm = (no*Z0e-ne*Z0o)/(2*c*Z0e*Z0o)
        return (L1 = L1, L2 = L2, M = M, C1 = C1, C2 = C2, Cm = Cm)
    else
        throw(ArgumentError("Canonical coupled line circuit number must be 1-10."))
    end
end

"""
    maxwell_combine(n::Int, d::Dict{NTuple{N, Int}, T}) where {N,T<:AbstractMatrix}

Return the Maxwell capacitance matrix for an `n` terminal system from the
Maxwell capacitance matrices for sets of terminals stored in the dictionary
`d`. The dictionary keys are tuples of the terminal numbers for the
capacitance matrices and the values are the capacitance matrices.

# Examples
```jldoctest
julia> @variables C11, C12, C13, C21, C22, C23, C31, C32, C33;JosephsonCircuits.maxwell_combine(3, Dict((1,2)=>[C11 C12;C21 C22],(1,3)=>[C11 C13;C31 C33],(2,3)=>[C22 C23;C32 C33]))
3×3 Matrix{Num}:
 C11  C12  C13
 C21  C22  C23
 C31  C32  C33

julia> @variables C11, C12, C13, C21, C22, C23, C31, C32, C33;JosephsonCircuits.maxwell_combine(3, Dict((1,2,3)=>[C11 C12 C13;C21 C22 C23; C31 C32 C33]))
3×3 Matrix{Num}:
 C11  C12  C13
 C21  C22  C23
 C31  C32  C33
```
"""
function maxwell_combine(n::Int, d::Dict{NTuple{N, Int}, T}) where {N,T<:AbstractMatrix}

    C = zeros(eltype(T), n, n)
    for (key,val) in d
        for i in 1:N
            for j in 1:N
                if iszero(C[key[i],key[j]])
                    C[key[i],key[j]] = val[i,j]
                else
                    C[key[i],key[j]] += val[i,j]
                    C[key[i],key[j]] /= 2
                end
            end
        end
    end
    return C
end
