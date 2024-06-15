
"""
  connectS(Sx::Array{T,N},k::Int,l::Int)

Connect ports `k' and `l` on the same `m` port microwave network represented by
the scattering parameter matrix `Sx`, resulting in an `(m-2)` port scattering
parameter matrix, as illustrated below:

Input network:

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

Output network:

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
function connectS(Sx::Array{T,N},k::Int,l::Int) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeS = NTuple{N}(ifelse(i<=2,size(Sx,i)-2,size(Sx,i)) for i in 1:ndims(Sx))

    # allocate an array of zeros of the same type as Sx
    Sout = zeros(T,sizeS)

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
  connectS(Sx::Array{T,N},Sy::Array{T,N},k::Int,l::Int)

Connect port `k' on an `m` port network, represented by the scattering
parameter matrix `Sx`, to port `l` on an `n` port network, represented by the
scattering parameter matrix `Sy`, resulting in a single `(m-2)` port
network, as illustrated below:

Input network:

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

Output network:

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
function connectS(Sx::Array{T,N},Sy::Array{T,N},k::Int,l::Int) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeSx = size(Sx)
    sizeSy = size(Sy)
    sizeS = NTuple{N}(ifelse(i<=2,sizeSx[i]+sizeSy[i]-2,sizeSx[i]) for i in 1:length(sizeSx))

    # allocate an array of zeros of the same type as Sx
    Sout = zeros(T,sizeS)

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
    StoZ(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to an impedance parameter matrix
`Z` and return the result. 

``Z=\\sqrt{z}(1_{\\!N}-S)^{-1}(1_{\\!N}+S)\\sqrt{z}``

``\\sqrt{z}^{-1}Z=(1_{\\!N}-S)^{-1}(1_{\\!N}+S)\\sqrt{z}``

``\\sqrt{z}^{-1}Z=(1_{\\!N}-S) \\div (1_{\\!N}+S)\\sqrt{z}``

``Z= \\sqrt{z}((1_{\\!N}-S) \\div (1_{\\!N}+S)\\sqrt{z})``

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
"""
function StoZ(S;portimpedances=50.0)

    Z = similar(S)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = zeros(eltype(S),size(S,1),size(S,2))
    
    sqrtportimpedances = sqrt.(portimpedances)

    if ndims(portimpedances) == 0
        # assume the port impedances are all the same for all ports and
        # frequencies. loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoZ!(view(Z,:,:,i),view(S,:,:,i),tmp,sqrtportimpedances)
        end
    else
        # assume the port impedances are given for each port and frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoZ!(view(Z,:,:,i),view(S,:,:,i),tmp,Diagonal(view(sqrtportimpedances,:,i)))
        end
    end
    return Z
end


"""
    StoZ!(Z::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)


"""
function StoZ!(Z::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)
    
    # compute (I + S)*sqrt(portimpedances)
    copy!(Z,S)
    for d in 1:size(Z,1)
        Z[d,d] += 1
    end
    rmul!(Z,sqrtportimpedances)

    # compute (I - S)
    copy!(tmp,S)
    rmul!(tmp,-1)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # factorize the matrix
    A = lu!(tmp)

    # perform the left division
    # compute (I - S) \ ((I + S)*sqrt(portimpedances))
    ldiv!(A,Z)

    # left multiply by sqrt(z)
    # compute sqrt(portimpedances)*((I - S) \ ((I + S)*sqrt(portimpedances)))
    lmul!(sqrtportimpedances,Z)
    return nothing
end


"""
    ZtoS(Z;portimpedances=50.0)

Convert the impedance parameter matrix `Z` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances.

``S=(\\sqrt{y}Z\\sqrt{y}-1_{\\!N})(\\sqrt{y}Z\\sqrt{y}+1_{\\!N})^{-1}``
``S =(\\sqrt{y}Z\\sqrt{y}+1_{\\!N})^{-1}(\\sqrt{y}Z\\sqrt{y}-1_{\\!N})``
``S =(\\sqrt{y}Z\\sqrt{y}+1_{\\!N}) \\div (\\sqrt{y}Z\\sqrt{y}-1_{\\!N})``

where ```\\sqrt{y}=(\\sqrt{z})^{-1}`` where `z` is a diagonal matrix of port impedances. 

First compute ``\tilde{Z} = \\sqrt{y}Z\\sqrt{y}``, then:

``S =(\\tilde{Z}+1_{\\!N}) \\div (\\tilde{Z}-1_{\\!N})``

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
"""
function ZtoS(Z;portimpedances=50.0)

    S = similar(Z)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = zeros(eltype(Z),size(Z,1),size(Z,2))
    
    oneoversqrtportimpedances = 1 ./sqrt.(portimpedances)

    # assume the port impedances are all the same for all ports and frequencies
    if ndims(portimpedances) == 0
        # # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(Z)[3:end])
            ZtoS!(view(S,:,:,i),view(Z,:,:,i),tmp,oneoversqrtportimpedances)
        end
    else
        for i in CartesianIndices(axes(Z)[3:end])
            ZtoS!(view(S,:,:,i),view(Z,:,:,i),tmp,Diagonal(view(oneoversqrtportimpedances,:,i)))
        end
    end
    return S
end


"""
    ZtoS!(S::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

"""
function ZtoS!(S::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)
    
    # compute \tilde{Z} = oneoversqrtportimpedances*Z*oneoversqrtportimpedances in S and tmp
    copy!(S,Z)
    rmul!(S,oneoversqrtportimpedances)
    lmul!(oneoversqrtportimpedances,S)
    copy!(tmp,S)

    # compute (\tilde{Z} - I)
    for d in 1:size(S,1)
        S[d,d] -= 1
    end

    # compute (\tilde{Z} + I)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # factorize the matrix
    A = lu!(tmp)

    # perform the left division
    # compute (\tilde{Z} + I) \ (\tilde{Z} - I)
    ldiv!(A,S)

    return nothing
end


"""
    StoY(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to an admittance parameter matrix
`Y` and return the result. 

``Y=\\sqrt{y}}(I_{N}-S)(I_{N}+S)^{-1}{\\sqrt{y}}``
``={\\sqrt {y}}(I_{N}+S)^{-1}(I_{N}-S){\\sqrt {y}}``

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
"""
function StoY(S;portimpedances=50.0)

    Y = similar(S)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = zeros(eltype(S),size(S,1),size(S,2))
    
    oneoversqrtportimpedances = 1 ./sqrt.(portimpedances)

    if ndims(portimpedances) == 0
        # assume the port impedances are all the same for all ports and
        # frequencies. loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoY!(view(Y,:,:,i),view(S,:,:,i),tmp,oneoversqrtportimpedances)
        end
    else
        # assume the port impedances are given for each port and frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoY!(view(Y,:,:,i),view(S,:,:,i),tmp,Diagonal(view(oneoversqrtportimpedances,:,i)))
        end
    end
    return Y
end


"""
    StoY!(Y::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportadmittances)


"""
function StoY!(Y::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)
    
    # compute (I - S)*oneoversqrtportimpedances
    copy!(Y,S)
    rmul!(Y,-1)
    for d in 1:size(Y,1)
        Y[d,d] += 1
    end
    rmul!(Y,oneoversqrtportimpedances)

    # compute (I + S)
    copy!(tmp,S)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # factorize the matrix
    A = lu!(tmp)

    # perform the left division
    # compute (I + S) \ ((I - S)*oneoversqrtportimpedances)
    ldiv!(A,Y)

    # left multiply by sqrt(z)
    # compute oneoversqrtportimpedances*((I - S) \ ((I + S)*oneoversqrtportimpedances))
    lmul!(oneoversqrtportimpedances,Y)
    return nothing
end


"""
    YtoS(Y;portimpedances=50.0)

Convert the admittance parameter matrix `Y` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances.

``S=(I_{N}-{\\sqrt {z}}Y{\\sqrt {z}})(I_{N}+{\\sqrt {z}}Y{\\sqrt {z}})^{-1}``
``   =(I_{N}+{\\sqrt {z}}Y{\\sqrt {z}})^{-1}(I_{N}-{\\sqrt {z}}Y{\\sqrt{z}}``

where ```\\sqrt{y}=(\\sqrt{z})^{-1}`` where `z` is a diagonal matrix of port impedances. 

First compute ``\tilde{Y} = \\sqrt{z}Y\\sqrt{z}``, then:

``S =(1_{\\!N}+\\tilde{Y}) \\div (1_{\\!N}-\\tilde{Y})``

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
"""
function YtoS(Y;portimpedances=50.0)

    S = similar(Y)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = zeros(eltype(Y),size(Y,1),size(Y,2))
    
    sqrtportimpedances = sqrt.(portimpedances)

    # assume the port impedances are all the same for all ports and frequencies
    if ndims(portimpedances) == 0
        # # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(Y)[3:end])
            YtoS!(view(S,:,:,i),view(Y,:,:,i),tmp,sqrtportimpedances)
        end
    else
        for i in CartesianIndices(axes(Y)[3:end])
            YtoS!(view(S,:,:,i),view(Y,:,:,i),tmp,Diagonal(view(sqrtportimpedances,:,i)))
        end
    end
    return S
end


"""
    YtoS!(S::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

"""
function YtoS!(S::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)
    
    # compute \tilde{Y} = sqrtportimpedances*Y*sqrtportimpedances in S and tmp
    copy!(S,Y)
    rmul!(S,sqrtportimpedances)
    lmul!(sqrtportimpedances,S)
    copy!(tmp,S)
    rmul!(S,-1)

    # compute (-\tilde{Y} + I)
    for d in 1:size(S,1)
        S[d,d] += 1
    end

    # compute (\tilde{Y} + I)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # factorize the matrix
    A = lu!(tmp)

    # perform the left division
    # compute (I + \tilde{Y}) \ (I - \tilde{Y})
    ldiv!(A,S)

    return nothing
end




















