
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
Communications Engineering, Second Edition. Artech House, 2006
"""
function StoT(S)

    T = similar(S)
    # make a view of T,S and loop
    # make a temporary array 
    tmp = similar(S,axes(S)[1:2])
    
    # assume the port impedances are all the same for all ports and
    # frequencies. loop over the dimensions of the array greater than 2
    for i in CartesianIndices(axes(S)[3:end])
        StoT!(view(T,:,:,i),view(S,:,:,i),tmp)
    end

    return T

end

"""
    StoT!(T::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix)

"""
function StoT!(T::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(T,1)÷2
    range2 = size(T,1)÷2+1:size(T,1)

    # compute [-S12 0; -S22 I]
    fill!(T,zero(eltype(T)))
    for d in range2
        T[d,d] = 1
    end

    T[range1,range1] .= -S[range1, range2]
    T[range2,range1] .= -S[range2, range2]

    # compute [-I S11; 0 S21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp[range1,range2] .= S[range1, range1]
    tmp[range2,range2] .= S[range2, range1]

    # perform the left division
    # compute [-I S11; 0 S21] \ [-S12 0; -S22 I]
    T .= tmp \ T

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

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
"""
function StoZ(S;portimpedances=50.0)

    Z = similar(S)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(S,axes(S)[1:2])
    
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

    # perform the left division
    # compute (I - S) \ ((I + S)*sqrt(portimpedances))
    Z .= tmp \ Z

    # left multiply by sqrt(z)
    # compute sqrt(portimpedances)*((I - S) \ ((I + S)*sqrt(portimpedances)))
    lmul!(sqrtportimpedances,Z)
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

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
"""
function StoY(S;portimpedances=50.0)

    Y = similar(S)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(S,axes(S)[1:2])
    
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

    # perform the left division
    # compute (I + S) \ ((I - S)*oneoversqrtportimpedances)
    Y .= tmp \ Y

    # left multiply by sqrt(z)
    # compute oneoversqrtportimpedances*((I - S) \ ((I + S)*oneoversqrtportimpedances))
    lmul!(oneoversqrtportimpedances,Y)
    return nothing
end


"""
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
Communications Engineering, Second Edition. Artech House, 2006
"""
function StoA(S;portimpedances=50.0)

    A = similar(S)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(S,axes(S)[1:2])
    
    sqrtportimpedances = sqrt.(portimpedances)

    if ndims(portimpedances) == 0
        # assume the port impedances are all the same for all ports and
        # frequencies. loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoA!(view(A,:,:,i),view(S,:,:,i),tmp,sqrtportimpedances, sqrtportimpedances)
        end
    else
        # assume the port impedances are given for each port and frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoA!(view(A,:,:,i),view(S,:,:,i),tmp,
                Diagonal(view(sqrtportimpedances,1:size(sqrtportimpedances,1)÷2,i)),
                Diagonal(view(sqrtportimpedances,size(sqrtportimpedances,1)÷2+1:size(sqrtportimpedances,1),i)))
        end
    end
    return A
end


"""
    StoA!(A::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

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
    tmp11 .= -(I - S11)/sqrtportimpedances1
    tmp12 .= (I + S11)*sqrtportimpedances1
    tmp21 .= S21/sqrtportimpedances1
    tmp22 .= S21*sqrtportimpedances1

    A11 .= -S12/sqrtportimpedances2
    A12 .= S12*sqrtportimpedances2
    A21 .= (I-S22)/sqrtportimpedances2
    A22 .= (I+S22)*sqrtportimpedances2

    # perform the left division
    A .= tmp \ A

    return nothing
end


"""
    StoB(S)

Convert the scattering parameter matrix `S` to the inverse chain (ABCD) matrix
`B` and return the result. Note that despite the name, the inverse of the chain
matrix is not equal to the inverse chain matrix, inv(A) ≠ B.

# Examples
```
# below won't work yet
julia> S = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.ZtoB(JosephsonCircuits.StoZ(S)),JosephsonCircuits.StoB(S))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
"""
function StoB(S;portimpedances=50.0)

    B = similar(S)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(S,axes(S)[1:2])
    
    sqrtportimpedances = sqrt.(portimpedances)

    if ndims(portimpedances) == 0
        # assume the port impedances are all the same for all ports and
        # frequencies. loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoB!(view(B,:,:,i),view(S,:,:,i),tmp,sqrtportimpedances, sqrtportimpedances)
        end
    else
        # assume the port impedances are given for each port and frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(S)[3:end])
            StoB!(view(B,:,:,i),view(S,:,:,i),tmp,
                Diagonal(view(sqrtportimpedances,1:size(sqrtportimpedances,1)÷2,i)),
                Diagonal(view(sqrtportimpedances,size(sqrtportimpedances,1)÷2:size(sqrtportimpedances,1),i)))
        end
    end
    return B
end


"""
    StoB!(B::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

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
    tmp11 .= S12./sqrtportimpedances2
    tmp12 .= S12.*sqrtportimpedances2
    tmp21 .= -(I-S22)./sqrtportimpedances2
    tmp22 .= (I+S22).*sqrtportimpedances2

    B11 .= (I-S11)./sqrtportimpedances1
    B12 .= (I+S11).*sqrtportimpedances1
    B21 .= -S21./sqrtportimpedances1
    B22 .= S21.*sqrtportimpedances1

    # perform the left division
    B .= tmp \ B

    return nothing
end


"""
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
"""
function TtoS(T)

    S = similar(T)
    # make a view of T,S and loop
    # make a temporary array 
    tmp = similar(T,axes(T)[1:2])
    
    # assume the port impedances are all the same for all ports and
    # frequencies. loop over the dimensions of the array greater than 2
    for i in CartesianIndices(axes(T)[3:end])
        TtoS!(view(S,:,:,i),view(T,:,:,i),tmp)
    end

    return S

end

"""
    TtoS!(S::AbstractMatrix,T::AbstractMatrix,tmp::AbstractMatrix)

"""
function TtoS!(S::AbstractMatrix,T::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(T,1)÷2
    range2 = size(T,1)÷2+1:size(T,1)

    # compute [0 T11; I T21]
    fill!(S,zero(eltype(S)))
    for d in range1
        S[d+size(T,1)÷2,d] = 1
    end

    S[range1,range2] .= -T[range1, range1]
    S[range2,range2] .= -T[range2, range1]

    # compute [-I T12; 0 T22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp[range1,range2] .= T[range1, range2]
    tmp[range2,range2] .= T[range2, range2]

    # perform the left division
    # compute [-I T12; 0 T22] \ [0 T11; I T21]
    S .= tmp \ S

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

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
"""
function ZtoS(Z;portimpedances=50.0)

    S = similar(Z)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(Z,axes(Z)[1:2])
    
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

    # perform the left division
    # compute (\tilde{Z} + I) \ (\tilde{Z} - I)
    S .= tmp \ S

    return nothing
end

"""
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
Communications Engineering, Second Edition. Artech House, 2006
"""
function ZtoA(Z)
    return AtoZ(Z)
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

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
"""
function YtoS(Y;portimpedances=50.0)

    S = similar(Y)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(Y,axes(Y)[1:2])
    
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

    # perform the left division
    # compute (I + \tilde{Y}) \ (I - \tilde{Y})
    S .= tmp \ S

    return nothing
end


"""
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
Communications Engineering, Second Edition. Artech House, 2006
"""
function AtoS(A;portimpedances=50.0)

    S = similar(A)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(A,axes(A)[1:2])
    
    sqrtportimpedances = sqrt.(portimpedances)

    if ndims(portimpedances) == 0
        # assume the port impedances are all the same for all ports and
        # frequencies. loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(A)[3:end])
            AtoS!(view(S,:,:,i),view(A,:,:,i),tmp,sqrtportimpedances, sqrtportimpedances)
        end
    else
        # assume the port impedances are given for each port and frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(A)[3:end])
            AtoS!(view(S,:,:,i),view(A,:,:,i),tmp,
                Diagonal(view(sqrtportimpedances,1:size(sqrtportimpedances,1)÷2,i)),
                Diagonal(view(sqrtportimpedances,size(sqrtportimpedances,1)÷2+1:size(sqrtportimpedances,1),i)))
        end
    end
    return S
end


"""
    AtoS!(S::AbstractMatrix, A::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

"""
function AtoS!(S::AbstractMatrix, A::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # compute 
    fill!(S,zero(eltype(S)))

    for d in range1
        S[d,d] = 1
    end

    for d in range1
        S[d+size(A,1)÷2,d] = 1
    end

    S[range1,range1] .= S[range1,range1]*sqrtportimpedances1
    S[range1,range2] .= -A[range1, range1]*sqrtportimpedances2+A[range1, range2]/sqrtportimpedances2
    S[range2,range1] .=  S[range2,range1]/sqrtportimpedances1
    S[range2,range2] .= -A[range2, range1]*sqrtportimpedances2+A[range2, range2]/sqrtportimpedances2


    fill!(tmp,zero(eltype(tmp)))

    for d in range1
        tmp[d,d] = 1
    end

    for d in range1
        tmp[d+size(A,1)÷2,d] = 1
    end

    tmp[range1,range1] .= tmp[range1,range1]*sqrtportimpedances1
    tmp[range1,range2] .= -A[range1, range1]*sqrtportimpedances2-A[range1, range2]/sqrtportimpedances2
    tmp[range2,range1] .= -tmp[range2,range1]/sqrtportimpedances1
    tmp[range2,range2] .= -A[range2, range1]*sqrtportimpedances2-A[range2, range2]/sqrtportimpedances2

    S .= -tmp \ S

    return nothing
end

"""
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
Communications Engineering, Second Edition. Artech House, 2006
"""
function AtoZ(A)

    Z = similar(A)

    # make a view of Z,A and loop
    # make a temporary array 
    tmp = similar(A,axes(A)[1:2])
    
    # assume the port impedances are all the same for all ports and
    # frequencies. loop over the dimensions of the array greater than 2
    for i in CartesianIndices(axes(A)[3:end])
        AtoZ!(view(Z,:,:,i),view(A,:,:,i),tmp)
    end

    return Z
end

"""
    AtoZ!(Z::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

"""
function AtoZ!(Z::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # compute [0 A12; I A22]
    fill!(Z,zero(eltype(A)))
    for d in range1
        Z[d+size(A,1)÷2,d] = 1
    end

    Z[range1,range2] .= A[range1, range2]
    Z[range2,range2] .= A[range2, range2]

    # compute [-I A11; 0 A21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp[range1,range2] .= A[range1, range1]
    tmp[range2,range2] .= A[range2, range1]

    # perform the left division
    # compute [-I A11; 0 A21] \ [0 A12; I A22]
    Z .= tmp \ Z

    return nothing

end


"""
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
Communications Engineering, Second Edition. Artech House, 2006
"""
function AtoB(A)

    B = similar(A)

    # make a view of Z,A and loop
    # make a temporary array 
    tmp = similar(A,axes(A)[1:2])
    
    # assume the port impedances are all the same for all ports and
    # frequencies. loop over the dimensions of the array greater than 2
    for i in CartesianIndices(axes(A)[3:end])
        AtoB!(view(B,:,:,i),view(A,:,:,i),tmp)
    end

    return B
end

"""
    AtoB!(B::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

"""
function AtoB!(B::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # compute [I 0; 0 I]
    fill!(B,zero(eltype(A)))
    for d in range1
        B[d,d] = 1
    end
    for d in range2
        B[d,d] = -1
    end

    # compute [A11 -A12; A21 -A22]
    copy!(tmp,A)
    tmp[range1,range2] .*= -1
    tmp[range2,range2] .*= -1

    # perform the left division
    # compute [A11 -A12; A21 -A22] \ [I 0; 0 I]
    B .= tmp \ B

    return nothing

end

"""
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
"""
function BtoS(B;portimpedances=50.0)

    S = similar(B)
    # make a view of Z,S and loop
    # make a temporary array 
    tmp = similar(B,axes(B)[1:2])
    
    sqrtportimpedances = sqrt.(portimpedances)

    if ndims(portimpedances) == 0
        # assume the port impedances are all the same for all ports and
        # frequencies. loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(B)[3:end])
            BtoS!(view(S,:,:,i),view(B,:,:,i),tmp,sqrtportimpedances, sqrtportimpedances)
        end
    else
        # assume the port impedances are given for each port and frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(B)[3:end])
            BtoS!(view(S,:,:,i),view(B,:,:,i),tmp,
                Diagonal(view(sqrtportimpedances,1:size(sqrtportimpedances,1)÷2,i)),
                Diagonal(view(sqrtportimpedances,size(sqrtportimpedances,1)÷2+1:size(sqrtportimpedances,1),i)))
        end
    end
    return S
end


"""
    BtoS!(S::AbstractMatrix, B::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

"""
function BtoS!(S::AbstractMatrix, B::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)


    fill!(tmp,zero(eltype(tmp)))

    for d in range1
        tmp[d,d+size(B,1)÷2] = 1
    end

    for d in range2
        tmp[d,d] = 1
    end

    # tmp[range1,range1] .= tmp[range1,range1]*sqrtportimpedances1
    # tmp[range1,range2] .= -B[range1, range1]*sqrtportimpedances2-B[range1, range2]/sqrtportimpedances2
    # tmp[range2,range1] .= -tmp[range2,range1]/sqrtportimpedances1
    # tmp[range2,range2] .= -B[range2, range1]*sqrtportimpedances2-B[range2, range2]/sqrtportimpedances2

    tmp[range1,range1] .= -B[range1, range1]*sqrtportimpedances1-B[range1, range2]/sqrtportimpedances1
    tmp[range1,range2] .= tmp[range1,range2]*sqrtportimpedances2
    tmp[range2,range1] .= -B[range2, range1]*sqrtportimpedances1-B[range2, range2]/sqrtportimpedances2
    tmp[range2,range2] .= -tmp[range2,range2]/sqrtportimpedances2


    # compute 
    fill!(S,zero(eltype(S)))

    for d in range1
        S[d,d+size(B,1)÷2] = 1
    end

    for d in range2
        S[d,d] = 1
    end

    # S[range1,range1] .= S[range1,range1]*sqrtportimpedances1
    # S[range1,range2] .= -B[range1, range1]*sqrtportimpedances2+B[range1, range2]/sqrtportimpedances2
    # S[range2,range1] .=  S[range2,range1]/sqrtportimpedances1
    # S[range2,range2] .= -B[range2, range1]*sqrtportimpedances2+B[range2, range2]/sqrtportimpedances2

    S[range1,range1] .= -B[range1, range1]*sqrtportimpedances1+B[range1, range2]/sqrtportimpedances1
    S[range1,range2] .= S[range1,range2]*sqrtportimpedances2
    S[range2,range1] .= -B[range2, range1]*sqrtportimpedances1+B[range2, range2]/sqrtportimpedances1
    S[range2,range2] .= S[range2,range2]/sqrtportimpedances2

    # println(tmp)
    # println(S)

    S .= -tmp \ S

    return nothing
end

"""
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
Communications Engineering, Second Edition. Artech House, 2006
"""
function BtoA(B)
    return AtoB(B)
end


"""
    ABCDtoS(ABCD;portimpedances=50.0)

Convert the 2 port chain (ABCD) matrix `ABCD` to the scattering parameter
matrix `S` and return the result.

# Examples
```jldoctest
julia> A = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.AtoS(A),JosephsonCircuits.ABCDtoS(A))
true

julia> A = rand(Complex{Float64},2,2,10);isapprox(JosephsonCircuits.AtoS(A),JosephsonCircuits.ABCDtoS(A))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
"""
function ABCDtoS(ABCD;portimpedances=50.0)

    S = similar(ABCD)
    # make a view of ABCD and loop
    
    # assume the port impedances are all the same for all ports and frequencies
    if ndims(portimpedances) == 0
        # # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(ABCD)[3:end])
            ABCDtoS!(view(S,:,:,i),view(ABCD,:,:,i),portimpedances)
        end
    else
        for i in CartesianIndices(axes(ABCD)[3:end])
            ABCDtoS!(view(S,:,:,i),view(ABCD,:,:,i),view(portimpedances,:,i))
        end
    end
    return S
end

# function ABCDtoS(ABCD;portimpedances=50.0)
#     return ABCDtoS(ABCD,first(portimpedances),last(portimpedances))
# end


"""
    StoABCD(S;portimpedances=50.0))

Convert the scattering parameter matrix `S` to the 2 port chain (ABCD) matrix and
return the result.

# Examples
```jldoctest
julia> S = rand(Complex{Float64},2,2);isapprox(JosephsonCircuits.StoA(S),JosephsonCircuits.StoABCD(S))
true

julia> S = rand(Complex{Float64},2,2,10);isapprox(JosephsonCircuits.StoA(S),JosephsonCircuits.StoABCD(S))
true
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
"""
function StoABCD(S;portimpedances=50.0)

    ABCD = similar(S)
    # make a view of ABCD and loop
    
    # assume the port impedances are all the same for all ports and frequencies
    if ndims(portimpedances) == 0
        # # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(ABCD)[3:end])
            StoABCD!(view(ABCD,:,:,i),view(S,:,:,i),portimpedances)
        end
    else
        for i in CartesianIndices(axes(ABCD)[3:end])
            StoABCD!(view(ABCD,:,:,i),view(S,:,:,i),view(portimpedances,:,i))
        end
    end
    return ABCD
end

# function StoABCD(S;portimpedances=50.0)
#     return StoABCD(S,first(portimpedances),last(portimpedances))
# end


function ABCDtos11(A,B,C,D,RS,RL)
    return (A*RL+B-C*RS*RL-D*RS)/(A*RL+B+C*RS*RL+D*RS)
end

function ABCDtos12(A,B,C,D,RS,RL)
    return 2*sqrt(RS*RL)*(A*D-B*C)/(A*RL+B+C*RS*RL+D*RS)
end

function ABCDtos21(A,B,C,D,RS,RL)
    return 2*sqrt(RS*RL)/(A*RL+B+C*RS*RL+D*RS)
end

function ABCDtos22(A,B,C,D,RS,RL)
    return (-A*RL+B-C*RS*RL+D*RS)/(A*RL+B+C*RS*RL+D*RS)
end

function StoA(S11,S12,S21,S22,RS,RL)
    return sqrt(RS/RL)*((1+S11)*(1-S22)+S21*S12)/(2*S21)
end

function StoB(S11,S12,S21,S22,RS,RL)
    return sqrt(RS*RL)*((1+S11)*(1+S22)-S21*S12)/(2*S21)
end

function StoC(S11,S12,S21,S22,RS,RL)
    return 1/sqrt(RS*RL)*((1-S11)*(1-S22)-S21*S12)/(2*S21)
end

function StoD(S11,S12,S21,S22,RS,RL)
    return sqrt(RL/RS)*((1-S11)*(1+S22)+S21*S12)/(2*S21)
end

function StoABCD!(ABCD,S,RS,RL)
    ABCD[1,1] = StoA(S[1,1],S[1,2],S[2,1],S[2,2],RS,RL)
    ABCD[1,2] = StoB(S[1,1],S[1,2],S[2,1],S[2,2],RS,RL)
    ABCD[2,1] = StoC(S[1,1],S[1,2],S[2,1],S[2,2],RS,RL)
    ABCD[2,2] = StoD(S[1,1],S[1,2],S[2,1],S[2,2],RS,RL)
    return nothing
end

function StoABCD!(ABCD,S,portimpedances)
    return StoABCD!(ABCD,S,first(portimpedances),last(portimpedances))
end

# function StoABCD(S,RS,RL)
#     ABCD = similar(S)
#     StoABCD!(ABCD,S,RS,RL)
#     return ABCD
# end


# function ABCDtoS(ABCD,RS,RL)
#     S = similar(ABCD)
#     ABCDtoS!(S,ABCD,RS,RL)
#     return S
# end

function ABCDtoS!(S,ABCD,portimpedances)
    return ABCDtoS!(S,ABCD,first(portimpedances),last(portimpedances))
end

function ABCDtoS!(S,ABCD,RS,RL)
    S[1,1] = ABCDtos11(ABCD[1,1],ABCD[1,2],ABCD[2,1],ABCD[2,2],RS,RL) 
    S[1,2] = ABCDtos12(ABCD[1,1],ABCD[1,2],ABCD[2,1],ABCD[2,2],RS,RL)
    S[2,1] = ABCDtos21(ABCD[1,1],ABCD[1,2],ABCD[2,1],ABCD[2,2],RS,RL)
    S[2,2] = ABCDtos22(ABCD[1,1],ABCD[1,2],ABCD[2,1],ABCD[2,2],RS,RL)
    return nothing
end




"""
    Z_coupled_tline(betae,betao,Z0e,Z0o,l)

betae, Z0e
betao, Z0o

V3, I3 -->  ======== <-- I4, V4
V1, I1 -->  ======== <-- I2, V2
            <---l-->

[(V1+V3)/2, (I1+I3)/2] = ABCDe * [(V2+V4)/2, -(I2+I4)/2]
[(V1-V3)/2, (I1-I3)/2] = ABCDo * [(V2-V4)/2, -(I2-I4)/2]

[V1, V2, V3, V4] = Z_coupled_tline * [I1, I2, I3, I4]

# Examples
```jldoctest
julia> JosephsonCircuits.Z_coupled_tline(1,1,50,50,pi/4)
4×4 Matrix{ComplexF64}:
 0.0-50.0im     0.0-70.7107im  0.0-0.0im      0.0-0.0im
 0.0-70.7107im  0.0-50.0im     0.0-0.0im      0.0-0.0im
 0.0-0.0im      0.0-0.0im      0.0-50.0im     0.0-70.7107im
 0.0-0.0im      0.0-0.0im      0.0-70.7107im  0.0-50.0im
 ```
"""
function Z_coupled_tline(betae,betao,Z0e,Z0o,l)
    Z11 = Z22 = Z33 = Z44 = -im/2*(Z0e*cot(betae*l) + Z0o*cot(betao*l))
    Z12 = Z21 = Z34 = Z43 = -im/2*(Z0e*csc(betae*l) + Z0o*csc(betao*l))
    Z13 = Z31 = Z42 = Z24 = -im/2*(Z0e*cot(betae*l) - Z0o*cot(betao*l))
    Z14 = Z41 = Z32 = Z23 = -im/2*(Z0e*csc(betae*l) - Z0o*csc(betao*l))
    return [Z11 Z12 Z13 Z14; Z21 Z22 Z23 Z24; Z31 Z32 Z33 Z34; Z41 Z42 Z43 Z44]
end

"""
    ABCD_tline(beta,Z0,l)

   <--l--->
o--========--o
   beta, Z0   
              
o------------o

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_tline(1,1,pi/4)
2×2 Matrix{ComplexF64}:
 0.707107+0.0im            0.0+0.707107im
      0.0+0.707107im  0.707107+0.0im
 ```
"""
function ABCD_tline(beta,Z0,l)
    return [cos(beta*l) im*Z0*sin(beta*l);im/Z0*sin(beta*l) cos(beta*l)]
end


"""
    ABCD_seriesZ(Z)

o---Z1---o
          
          
o--------o

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_seriesZ(2)
2×2 Matrix{Int64}:
 1  2
 0  1
 ```
"""
function ABCD_seriesZ(Z)
    return [1 Z;0 1]
end

"""
    ABCD_shuntY(Y)

---------
    |
    Y1
    |
---------

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_shuntY(2)
2×2 Matrix{Int64}:
 1  0
 2  1
 ```
"""
function ABCD_shuntY(Y)
    return [1 0;Y 1]
end

"""
    ABCD_PiY(Y1,Y2,Y3)

o----Y3-----o
   |     |   
   Y1    Y2  
   |     |   
o-----------o

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
    ABCD_TZ(Z1,Z2,Z3)

o--Z1-----Z2--o
       |       
      Z3       
       |       
o-------------o

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

