
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