
"""
    intraconnectS(Sa::AbstractArray, k::Int, l::Int;
        nbatches::Int = Base.Threads.nthreads())

Connect ports `k` and `l` on the same `m` port microwave network represented
by the scattering parameter matrix `Sa`, resulting in an `(m-2)` port network,
as illustrated below:

Input network:
```
      m |         | l+1    
        |   ...   |         l
        |_________|__________ 
        |         |          |
        |   Sa    |  ...     |
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
- `Sa::Array`: Array of scattering parameters representing the network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `k::Int`: First port to connect, with one based indexing.
- `l::Int`: Second port to connect, with one based indexing.

# References
R. C. Compton and D. B. Rutledge, "Perspectives in Microwave Circuit
Analysis," Proceedings of the 32nd Midwest Symposium on Circuits and Systems,
vol. 2, pp. 716–718, Aug. 1989. doi: 10.1109/MWSCAS.1989.101955
V. A. Monaco and P. Tiberio, "Computer-Aided Analysis of Microwave Circuits,"
in IEEE Transactions on Microwave Theory and Techniques, vol. 22, no. 3, pp.
249-263, Mar. 1974, doi: 10.1109/TMTT.1974.1128208.
"""
function intraconnectS(Sa::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeS = NTuple{N}(ifelse(i<=2,size(Sa,i)-2,size(Sa,i)) for i in 1:ndims(Sa))

    # allocate an array of zeros of the same type as Sa
    Sout = similar(Sa,sizeS)

    # remove the self loop
    intraconnectS!(Sout,Sa,k,l;nbatches = nbatches)

    return Sout
end

"""
    intraconnectS!(Sout, Sa, k::Int, l::Int; nbatches::Int = Base.Threads.nthreads())

See [`intraconnectS`](@ref) for description.

"""
function intraconnectS!(Sout, Sa, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())

    # validate all of the inputs
    if ndims(Sa) != ndims(Sout)
        throw(DimensionMismatch("`Sout` and `Sa` must have the same number of dimensions."))
    end

    if ndims(Sa) < 2
        throw(DimensionMismatch("`Sout` and `Sa` must have at least two dimensions."))
    end

    if size(Sa,1) != size(Sa,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."))
    end

    if size(Sout,1) != size(Sout,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."))
    end

    if size(Sa,1) -2 != size(Sout,1)
        throw(DimensionMismatch("Length of first two dimensions must be 2 smaller for `Sout` than `Sa` because we are merging two ports."))
    end

    for i in 3:ndims(Sa)
        if size(Sa,i) != size(Sout,i)
            throw(DimensionMismatch("Non-port axis lengths of `Sa` and `Sout` must be equal."))
        end
    end

    if k > size(Sa,1)
        throw(ArgumentError("Port `k` is larger than number of ports in `Sa`."))
    end

    if l > size(Sa,1)
        throw(ArgumentError("Port `l` is larger than number of ports in `Sa`."))
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

    # loop over the dimensions of the array greater than 2
    indices = CartesianIndices(axes(Sout)[3:end])
    if  nbatches > 1 && length(indices) > nbatches
        batches = Base.Iterators.partition(1:length(indices),1+(length(indices)-1)÷nbatches)
        Threads.@sync for batch in batches
            Base.Threads.@spawn intraconnectS_inner!(Sout,Sa,k,l,batch)
        end
    else
        intraconnectS_inner!(Sout,Sa,k,l,indices)
    end

    return Sout
end

"""
    intraconnectS_inner!(Sout, Sa, k::Int, l::Int, batch::AbstractArray)

See [`intraconnectS`](@ref) for description.

"""
function intraconnectS_inner!(Sout, Sa, k::Int, l::Int, batch::AbstractArray)

    # the number of ports in the input matrix
    m = size(Sa,1)

    # order the ports as k < l
    k, l = ifelse(k>l,(l,k),(k,l))

    range1 = 1:k-1
    range2 = k+1:l-1
    range3 = l+1:m
    ranges = (range1,range2,range3)

    @inbounds for b in batch
        gammacc_Scc = StaticArrays.SMatrix{2,2}(
            -Sa[k,k,b],
            one(Sa[l,k,b])-Sa[l,k,b],
            one(Sa[k,l,b])-Sa[k,l,b],
            -Sa[l,l,b]
        )
        gammacc_Scc_lu = lu_2x2(gammacc_Scc)
        # gammacc_Scc_lu =  lu(gammacc_Scc)

        # ii and jj are the indices which extend up to m and skip k,l
        # i and j extend up to m-2 and are consecutive 
        for jindex in eachindex(ranges)
            for jj in ranges[jindex]
                j = jj-jindex+1

                Scp = StaticArrays.SVector{2}(Sa[k,jj,b],Sa[l,jj,b])

                # solve the linear system
                ac1jj, ac2jj = ldiv_2x2(gammacc_Scc_lu,Scp)
                # ac1jj, ac2jj = gammacc_Scc_lu \ Scp

                for iindex in eachindex(ranges)
                    for ii in ranges[iindex]
                        i = ii-iindex+1
                        Sout[i,j,b] = Sa[ii,jj,b] + Sa[ii,k,b]*ac1jj + Sa[ii,l,b]*ac2jj
                    end
                end
            end
        end
    end

    return Sout
end

"""
    intraconnectS(Sa::AbstractArray, Ca::AbstractArray, k::Int, l::Int;
        nbatches::Int = Base.Threads.nthreads())

Connect ports `k` and `l` on the same `m` port microwave network represented
by the scattering parameter matrix `Sa`, and noise correlation matrix `Ca` 
resulting in an `(m-2)` port network, as illustrated below:

Input network:
```
      m |         | l+1    
        |   ...   |         l
        |_________|__________ 
        |         |          |
        |   Sa    |  ...     |
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
- `Sa::Array`: Array of scattering parameters representing the network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `Ca::Array`: Array of noise correlation parameters of the same dimensions as
    `Sa`.
- `k::Int`: First port to connect, with one based indexing.
- `l::Int`: Second port to connect, with one based indexing.

# References
S. W. Wedge, "Computer-aided design of low noise microwave circuits," PhD
thesis (1991).
R. C. Compton and D. B. Rutledge, "Perspectives in Microwave Circuit
Analysis," Proceedings of the 32nd Midwest Symposium on Circuits and Systems,
vol. 2, pp. 716–718, Aug. 1989. doi: 10.1109/MWSCAS.1989.101955
V. A. Monaco and P. Tiberio, "Computer-Aided Analysis of Microwave Circuits,"
in IEEE Transactions on Microwave Theory and Techniques, vol. 22, no. 3, pp.
249-263, Mar. 1974, doi: 10.1109/TMTT.1974.1128208.
"""
function intraconnectS(Sa::AbstractArray{T,N}, Ca::AbstractArray{T,N}, k::Int,
    l::Int; nbatches::Int = Base.Threads.nthreads()) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeS = NTuple{N}(ifelse(i<=2,size(Sa,i)-2,size(Sa,i)) for i in 1:ndims(Sa))

    # allocate an array of zeros of the same type as Sa
    Sout = similar(Sa,sizeS)
    Cout = similar(Ca,sizeS)

    # remove the self loop
    intraconnectS!(Sout,Cout,Sa,Ca,k,l;nbatches = nbatches)

    return Sout, Cout
end

"""
    intraconnectS!(Sout, Cout, Sa, Ca, k::Int, l::Int; nbatches::Int = Base.Threads.nthreads())

See [`intraconnectS`](@ref) for description.

"""
function intraconnectS!(Sout, Cout, Sa, Ca, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())

    # validate all of the inputs
    if ndims(Sa) != ndims(Sout)
        throw(DimensionMismatch("`Sout` and `Sa` must have the same number of dimensions."))
    end

    if ndims(Sa) < 2
        throw(DimensionMismatch("`Sout` and `Sa` must have at least two dimensions."))
    end

    if size(Sa,1) != size(Sa,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."))
    end

    if size(Sout,1) != size(Sout,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."))
    end

    if size(Sa,1) -2 != size(Sout,1)
        throw(DimensionMismatch("Length of first two dimensions must be 2 smaller for `Sout` than `Sa` because we are merging two ports."))
    end

    for i in 3:ndims(Sa)
        if size(Sa,i) != size(Sout,i)
            throw(DimensionMismatch("Non-port axis lengths of `Sa` and `Sout` must be equal."))
        end
    end

    if k > size(Sa,1)
        throw(ArgumentError("Port `k` is larger than number of ports in `Sa`."))
    end

    if l > size(Sa,1)
        throw(ArgumentError("Port `l` is larger than number of ports in `Sa`."))
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

    if size(Ca) != size(Sa)
        throw(DimensionMismatch("The size of `Ca` must the same as the size of `Sa`."))
    end

    if size(Cout) != size(Sout)
        throw(DimensionMismatch("The size of `Cout` must the same as the size of `Sout`."))
    end

    # loop over the dimensions of the array greater than 2
    indices = CartesianIndices(axes(Sout)[3:end])
    if  nbatches > 1 && length(indices) > nbatches
        batches = Base.Iterators.partition(1:length(indices),1+(length(indices)-1)÷nbatches)
        Threads.@sync for batch in batches
            Base.Threads.@spawn intraconnectS_inner!(Sout,Cout,Sa,Ca,k,l,batch)
        end
    else
        intraconnectS_inner!(Sout,Cout,Sa,Ca,k,l,indices)
    end

    return Sout, Cout
end

function intraconnectS_inner!(Sout, Cout, Sa, Ca, k::Int, l::Int,
    batch::AbstractArray)

    # the number of ports in the input matrix
    m = size(Sa,1)

    # order the ports as k < l
    k, l = ifelse(k>l,(l,k),(k,l))

    range1 = 1:k-1
    range2 = k+1:l-1
    range3 = l+1:m
    ranges = (range1,range2,range3)

    il_lk_ll_ik = similar(Sa,m-2)
    ik_kl_kk_il = similar(Sa,m-2)

    @inbounds for b in batch
        gammacc_Scc = StaticArrays.SMatrix{2,2}(
            -Sa[k,k,b],
            one(Sa[k,l,b])-Sa[k,l,b],
            one(Sa[l,k,b])-Sa[l,k,b],
            -Sa[l,l,b]
        )
        gammacc_Scc_lu = lu_2x2(gammacc_Scc)
        # gammacc_Scc_lu =  lu(gammacc_Scc)

        # compute the terms we will use in the inner loop
        for iindex in eachindex(ranges)
            for ii in ranges[iindex]
                i = ii-iindex+1

                # solve the linear system to compute il_lk_ll_ik and
                # ik_kl_kk_il where
                # il_lk_ll_ik = (Sa[i,l]*(1-Sa[l,k])+Sa[l,l]*Sa[i,k])/denom
                # ik_kl_kk_il = (Sa[i,k]*(1-Sa[k,l])+Sa[k,k]*Sa[i,l])/denom
                # denom = (1-Sa[k,l])*(1-Sa[l,k]) - Sa[k,k]*Sa[l,l]
                Scp = StaticArrays.SVector{2}(Sa[ii,k,b],Sa[ii,l,b])
                il_lk_ll_ik[i], ik_kl_kk_il[i] = ldiv_2x2(gammacc_Scc_lu,Scp)
                # il_lk_ll_ik[i], ik_kl_kk_il[i] = gammacc_Scc_lu \ Scp

            end
        end

        # ii and jj are the indices which extend up to m and skip k,l
        # i and j extend up to m-2 and are consecutive 
        for jindex in eachindex(ranges)
            for jj in ranges[jindex]
                j = jj-jindex+1

                # solve the linear system to compute jl_lk_ll_jk and
                # jk_kl_kk_jl where
                # jl_lk_ll_jk = (Sa[j,l]*(1-Sa[l,k])+Sa[l,l]*Sa[j,k])/denom
                # jk_kl_kk_jl = (Sa[j,k]*(1-Sa[k,l])+Sa[k,k]*Sa[j,l])/denom
                # denom = (1-Sa[k,l])*(1-Sa[l,k]) - Sa[k,k]*Sa[l,l]
                Scp = StaticArrays.SVector{2}(Sa[jj,k,b],Sa[jj,l,b])
                jl_lk_ll_jk, jk_kl_kk_jl = ldiv_2x2(gammacc_Scc_lu,Scp)
                # jl_lk_ll_jk, jk_kl_kk_jl = gammacc_Scc_lu \ Scp

                for iindex in eachindex(ranges)
                    for ii in ranges[iindex]
                        i = ii-iindex+1

                        # compute the scattering matrix
                        # Wedge thesis Eq. 3.14
                        Sout[i,j,b] = Sa[ii,jj,b] +
                           Sa[l,jj,b]*ik_kl_kk_il[i] +
                           Sa[k,jj,b]*il_lk_ll_ik[i]

                        # compute the noise correlation matrix
                        # Wedge thesis Eq. 3.17
                        Cout[i,j,b] = Ca[ii,jj,b] +
                            Ca[l,k,b]*ik_kl_kk_il[i]*conj(jl_lk_ll_jk) +
                            Ca[k,l,b]*il_lk_ll_ik[i]*conj(jk_kl_kk_jl) +
                            Ca[l,l,b]*ik_kl_kk_il[i]*conj(jk_kl_kk_jl) +
                            Ca[k,k,b]*il_lk_ll_ik[i]*conj(jl_lk_ll_jk) +
                            Ca[l,jj,b]*ik_kl_kk_il[i] +
                            Ca[k,jj,b]*il_lk_ll_ik[i] +
                            Ca[ii,l,b]*conj(jk_kl_kk_jl) +
                            Ca[ii,k,b]*conj(jl_lk_ll_jk)
                    end
                end
            end
        end
    end
    return Sout, Cout
end

"""
    interconnectS(Sa::AbstractArray, Sb::AbstractArray, Ca::AbstractArray,
        Cb::AbstractArray,k::Int, l::Int;
        nbatches::Int = Base.Threads.nthreads())

Connect port `k` on an `m` port network, represented by the scattering
parameter matrix `Sa`, to port `l` on an `n` port network, represented by the
scattering parameter matrix `Sb`, resulting in a single `(m+n-2)` port
network, as illustrated below:

Input network:
```
      m |        | k+1                       | 2
        |        |                           |
        |   ...  |                     ...   |
        |________|                  _________|________
        |        |                  |        |       1
        |   Sa   |                  |   Sb   |
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
- `Sa::Array`: Array of scattering parameters representing the first network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `Sb::Array`: Array of scattering parameters representing the second network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `Ca::Array`: Array of noise correlation parameters of the same dimensions as
    `Sa`.
- `Cb::Array`: Array of noise correlation parameters of the same dimensions as
    `Sb`.
- `k::Int`: Port on first network, with one based indexing.
- `l::Int`: Port on second network, with one based indexing.

# References
S. W. Wedge, "Computer-aided design of low noise microwave circuits," PhD
thesis (1991).
R. C. Compton and D. B. Rutledge, "Perspectives in Microwave Circuit
Analysis," Proceedings of the 32nd Midwest Symposium on Circuits and Systems,
vol. 2, pp. 716–718, Aug. 1989. doi: 10.1109/MWSCAS.1989.101955
V. A. Monaco and P. Tiberio, "Computer-Aided Analysis of Microwave Circuits,"
in IEEE Transactions on Microwave Theory and Techniques, vol. 22, no. 3, pp.
249-263, Mar. 1974, doi: 10.1109/TMTT.1974.1128208.
"""
function interconnectS(Sa::AbstractArray{T,N}, Sb::AbstractArray{T,N},
    Ca::AbstractArray{T,N}, Cb::AbstractArray{T,N}, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads()) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeSa = size(Sa)
    sizeSb = size(Sb)
    sizeS = NTuple{N}(ifelse(i<=2,sizeSa[i]+sizeSb[i]-2,sizeSa[i]) for i in 1:length(sizeSa))

    # allocate an array of zeros of the same type as Sa
    Sout = similar(Sa,sizeS)
    Cout = similar(Ca,sizeS)

    # connect the networks
    interconnectS!(Sout,Cout,Sa,Sb,Ca,Cb,k,l;nbatches = nbatches)

    return Sout, Cout
end


"""
    interconnectS!(Sout, Cout, Sa, Sb, Ca, Cb, k, l)

See [`interconnectS`](@ref) for description.

"""
function interconnectS!(Sout, Cout, Sa, Sb, Ca, Cb, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())

    # validate all of the inputs
    if ndims(Sa) != ndims(Sb)
        throw(DimensionMismatch("`Sa` and `Sb` must have the same number of dimensions."))
    end

    if ndims(Sa) != ndims(Sout)
        throw(DimensionMismatch("`Sout`, `Sa`, and `Sb` must have the same number of dimensions."))
    end

    if ndims(Sa) < 2
        throw(DimensionMismatch("`Sout`, `Sa`, and `Sb` must have atleast two dimensions."))
    end

    if size(Sa,1) != size(Sa,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."))
    end

    if size(Sb,1) != size(Sb,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sb` must be equal."))
    end

    if size(Sout,1) != size(Sout,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."))
    end

    if size(Sa,1) + size(Sb,1) - 2 != size(Sout,1)
        throw(DimensionMismatch("First two dimensions of `Sout` must be `m+n-2`."))
    end

    for i in 3:ndims(Sa)
        if size(Sa,i) != size(Sout,i)
            throw(DimensionMismatch("Non-port axis lengths of `Sa`, `Sb`, and `Sout` must be equal."))
        end
    end

    if k > size(Sa,1)
        throw(ArgumentError("Port `k` is larger than number of ports in `Sa`."))
    end

    if l > size(Sb,1)
        throw(ArgumentError("Port `l` is larger than number of ports in `Sb`."))
    end

    if l < 1
        throw(ArgumentError("Port `l` is smaller than one."))
    end

    if k < 1
        throw(ArgumentError("Port `k` is smaller than one."))
    end

    if size(Ca) != size(Sa)
        throw(DimensionMismatch("The size of `Ca` must the same as the size of `Sa`."))
    end

    if size(Cb) != size(Sb)
        throw(DimensionMismatch("The size of `Cb` must the same as the size of `Sb`."))
    end

    if size(Cout) != size(Sout)
        throw(DimensionMismatch("The size of `Cout` must the same as the size of `Sout`."))
    end

    # loop over the dimensions of the array greater than 2
    indices = CartesianIndices(axes(Sout)[3:end])
    if nbatches > 1 && length(indices) > nbatches
        batches = Base.Iterators.partition(1:length(indices),1+(length(indices)-1)÷nbatches)
        Threads.@sync for batch in batches
            Base.Threads.@spawn interconnectS_inner!(Sout,Cout,Sa,Sb,Ca,Cb,k,l,batch)
        end

    else
        interconnectS_inner!(Sout,Cout,Sa,Sb,Ca,Cb,k,l,indices)
    end

    return Sout
end

"""
    interconnectS_inner!(Sout, Cout, Sa, Sb, Ca, Cb, k::Int, l::Int,
        batch::AbstractArray)

See [`interconnectS`](@ref) for description.

"""
function interconnectS_inner!(Sout, Cout, Sa, Sb, Ca, Cb, k::Int, l::Int,
    batch::AbstractArray)

    # the number of ports in the input matrix
    m = size(Sa,1)
    n = size(Sb,1)

    range1a = 1:k-1
    range1b = k+1:m
    range2a = m+1:m+l-1
    range2b = m+l+1:m+n

    # this indexes across the first part
    ranges1 = (range1a, range1b)

    # this indexes across the second part
    ranges2 = (range2a, range2b)

    a_ik = similar(Sa,m)
    a_ik_b_ll = similar(Sa,m)
    b_il = similar(Sb,n)
    b_il_a_kk = similar(Sb,n)

    # loop over the axes of the scattering parameter matrices after the first
    # two (eg. frequencies).
    @inbounds for b in batch

        gammacc_Scc = StaticArrays.SMatrix{2,2}(
            -Sa[k,k,b],
            one(Sa[k,k,b]),
            one(Sa[k,k,b]),
            -Sb[l,l,b]
        )
        # gammacc_Scc_lu = lu(gammacc_Scc)
        gammacc_Scc_lu = lu_2x2(gammacc_Scc)

        # solve the linear system to compute Sa[i,k]*Sb[l,l]/denom and
        # Sa[i,k]/denom
        for iindex in eachindex(ranges1)
            for ii in ranges1[iindex]
                Scp = StaticArrays.SVector{2}(Sa[ii,k,b],0)
                a_ik_b_ll[ii], a_ik[ii] = ldiv_2x2(gammacc_Scc_lu,Scp)
            end
        end

        # solve the linear system to compute Sb[i,l]*Sa[k,k]/denom and
        # Sb[i,l]/denom
        for iindex in eachindex(ranges2)
            for ii in ranges2[iindex]
                Scp = StaticArrays.SVector{2}(0,Sb[ii-m,l,b])
                b_il[ii-m], b_il_a_kk[ii-m] = ldiv_2x2(gammacc_Scc_lu,Scp)
            end
        end

        # # calculate the denominator, for debugging
        # denom = one(Sa[k,k,b])-Sa[k,k,b]*Sb[l,l,b]

        # ii and jj index across the concatenated Sa and Sb arrays skipping
        # the k'th and l'th elements. Imagine the Sa nxn array in the upper
        # left quadrant and the Sb mxm in the lower right quadrant. These are
        # used to pick out the elements of Sa and Sb. i and j index into the 
        # m+n-2 output array. These indices are contiguous.

        # left side
        for jindex in eachindex(ranges1)

            for jj in ranges1[jindex]
                j = jj-jindex+1

                # upper left quadrant
                for iindex in eachindex(ranges1)
                    for ii in ranges1[iindex]
                        i = ii-iindex+1

                        # Eq. 3.15a from Wedge thesis
                        Sout[i,j,b] = Sa[ii,jj,b] + Sa[k,jj,b]*a_ik_b_ll[ii]
                        # Eq. 3.19a from Wedge thesis
                        Cout[i,j,b] = Ca[ii,jj,b] +
                            Ca[ii,k,b]*conj(a_ik_b_ll[jj]) +
                            Ca[k,jj,b]*a_ik_b_ll[ii] +
                            Cb[l,l,b]*a_ik[ii]*conj(a_ik[jj]) +
                            Ca[k,k,b]*a_ik_b_ll[ii]*conj(a_ik_b_ll[jj])

                        # # a non-optimized version of the above for debugging
                        # Sout[i,j,b] = Sa[ii,jj,b] + Sa[k,jj,b]*Sb[l,l,b]*Sa[ii,k,b]/denom
                        # Cout[i,j,b] = Ca[ii,jj,b] +
                        #     Ca[ii,k,b]*conj(Sa[jj,k,b]*Sb[l,l,b]/denom) +
                        #     Ca[k,jj,b]*Sa[ii,k,b]*Sb[l,l,b]/denom +
                        #     Cb[l,l,b]*Sa[ii,k,b]/denom*conj(Sa[jj,k,b]/denom) +
                        #     Ca[k,k,b]*Sb[l,l,b]*Sa[ii,k,b]/denom*conj(Sb[l,l,b]*Sa[jj,k,b]/denom)
                    end
                end

                # lower left quadrant
                for iindex in eachindex(ranges2)
                    for ii in ranges2[iindex]
                        i = ii-iindex+1-1

                        # Eq. 3.15b from Wedge thesis
                        Sout[i,j,b] = Sa[k,jj,b]*b_il[ii-m]
                        # Eq. 3.19b from Wedge thesis
                        Cout[i,j,b] = Cb[ii-m,l,b]*conj(a_ik[jj]) +
                            Ca[k,jj,b]*b_il[ii-m] +
                            Cb[l,l,b]*b_il_a_kk[ii-m]*conj(a_ik[jj]) +
                            Ca[k,k,b]*b_il[ii-m]*conj(a_ik_b_ll[jj])

                        # # a non-optimized version of the above for debugging
                        # Sout[i,j,b] = Sa[k,jj,b]*Sb[ii-m,l,b]/denom
                        # Cout[i,j,b] = Cb[ii-m,l,b]*conj(Sa[jj,k,b]/denom) +
                        #     Ca[k,jj,b]*Sb[ii-m,l,b]/denom +
                        #     Cb[l,l,b]*Sb[ii-m,l,b]*Sa[k,k,b]/denom*conj(Sa[jj,k,b]/denom) +
                        #     Ca[k,k,b]*Sb[ii-m,l,b]/denom*conj(Sa[jj,k,b]*Sb[l,l,b]/denom)
                    end
                end

            end

        end

        # right side
        for jindex in eachindex(ranges2)

            for jj in ranges2[jindex]
                j = jj-jindex+1-1

                # upper right quadrant
                for iindex in eachindex(ranges1)
                    for ii in ranges1[iindex]
                        i = ii-iindex+1

                        # Eq. 3.15b from Wedge thesis
                        Sout[i,j,b] = Sb[l,jj-m,b]*a_ik[ii]
                        # Eq. 3.19b from Wedge thesis
                        Cout[i,j,b] = Ca[ii,k,b]*conj(b_il[jj-m]) +
                            Cb[l,jj-m,b]*a_ik[ii] +
                            Ca[k,k,b]*a_ik_b_ll[ii]*conj(b_il[jj-m]) +
                            Cb[l,l,b]*a_ik[ii]*conj(b_il_a_kk[jj-m])

                        # # a non-optimized version of the above for debugging
                        # Sout[i,j,b] = Sb[l,jj-m,b]*Sa[ii,k,b]/denom
                        # Cout[i,j,b] = Ca[ii,k,b]*conj(Sb[jj-m,l,b]/denom) +
                        #     Cb[l,jj-m,b]*Sa[ii,k,b]/denom +
                        #     Ca[k,k,b]*Sa[ii,k,b]*Sb[l,l,b]/denom*conj(Sb[jj-m,l,b]/denom) +
                        #     Cb[l,l,b]*Sa[ii,k,b]/denom*conj(Sb[jj-m,l,b]*Sa[k,k,b]/denom)
                    end
                end

                # lower right quadrant
                for iindex in eachindex(ranges2)
                    for ii in ranges2[iindex]
                        i = ii-iindex+1-1

                        # Eq. 3.15a from Wedge thesis
                        Sout[i,j,b] = Sb[ii-m,jj-m,b] + Sb[l,jj-m,b]*b_il_a_kk[ii-m]
                        # Eq. 3.19a from Wedge thesis
                        Cout[i,j,b] = Cb[ii-m,jj-m,b] +
                            Cb[ii-m,l,b]*conj(b_il_a_kk[jj-m]) +
                            Cb[l,jj-m,b]*b_il_a_kk[ii-m] +
                            Ca[k,k,b]*b_il[ii-m]*conj(b_il[jj-m]) +
                            Cb[l,l,b]*b_il_a_kk[ii-m]*conj(b_il_a_kk[jj-m])

                        # # a non-optimized version of the above for debugging
                        # Sout[i,j,b] = Sb[ii-m,jj-m,b] + Sb[l,jj-m,b]*Sa[k,k,b]*Sb[ii-m,l,b]/denom
                        # Cout[i,j,b] = Cb[ii-m,jj-m,b] +
                        #     Cb[ii-m,l,b]*conj(Sb[jj-m,l,b]*Sa[k,k,b]/denom) +
                        #     Cb[l,jj-m,b]*Sb[ii-m,l,b]*Sa[k,k,b]/denom +
                        #     Ca[k,k,b]*Sb[ii-m,l,b]/denom*conj(Sb[jj-m,l,b]/denom) +
                        #     Cb[l,l,b]*Sa[k,k,b]*Sb[ii-m,l,b]/denom*conj(Sa[k,k,b]*Sb[jj-m,l,b]/denom)
                    end
                end

            end

        end

    end

    return Sout, Cout

end

"""
    interconnectS(Sa::AbstractArray, Sb::AbstractArray, k::Int, l::Int;
        nbatches::Int = Base.Threads.nthreads())

Connect port `k` on an `m` port network, represented by the scattering
parameter matrix `Sa`, to port `l` on an `n` port network, represented by the
scattering parameter matrix `Sb`, resulting in a single `(m+n-2)` port
network, as illustrated below:

Input network:
```
      m |        | k+1                       | 2
        |        |                           |
        |   ...  |                     ...   |
        |________|                  _________|________
        |        |                  |        |       1
        |   Sa   |                  |   Sb   |
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
- `Sa::Array`: Array of scattering parameters representing the first network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `Sb::Array`: Array of scattering parameters representing the second network
    with ports along first two dimensions, followed by an arbitrary number
    of other dimensions (eg. frequency).
- `k::Int`: Port on first network, with one based indexing.
- `l::Int`: Port on second network, with one based indexing.

# References
V. A. Monaco and P. Tiberio, "Computer-Aided Analysis of Microwave Circuits,"
in IEEE Transactions on Microwave Theory and Techniques, vol. 22, no. 3, pp.
249-263, Mar. 1974, doi: 10.1109/TMTT.1974.1128208.
"""
function interconnectS(Sa::AbstractArray{T,N}, Sb::AbstractArray{T,N}, k::Int,
    l::Int; nbatches::Int = Base.Threads.nthreads()) where {T,N}

    # make a tuple with the size of the array
    # the first two dimensions are two smaller
    sizeSa = size(Sa)
    sizeSb = size(Sb)
    sizeS = NTuple{N}(ifelse(i<=2,sizeSa[i]+sizeSb[i]-2,sizeSa[i]) for i in 1:length(sizeSa))

    # allocate an array of zeros of the same type as Sa
    Sout = similar(Sa,sizeS)

    # connect the networks
    interconnectS!(Sout,Sa,Sb,k,l;nbatches = nbatches)

    return Sout
end

"""
    interconnectS!(Sout, Sa, Sb, k, l)

See [`interconnectS`](@ref) for description.

"""
function interconnectS!(Sout, Sa, Sb, k::Int, l::Int;
    nbatches::Int = Base.Threads.nthreads())

    # validate all of the inputs
    if ndims(Sa) != ndims(Sb)
        throw(DimensionMismatch("`Sa` and `Sb` must have the same number of dimensions."))
    end

    if ndims(Sa) != ndims(Sout)
        throw(DimensionMismatch("`Sout`, `Sa`, and `Sb` must have the same number of dimensions."))
    end

    if ndims(Sa) < 2
        throw(DimensionMismatch("`Sout`, `Sa`, and `Sb` must have atleast two dimensions."))
    end

    if size(Sa,1) != size(Sa,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sa` must be equal."))
    end

    if size(Sb,1) != size(Sb,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sb` must be equal."))
    end

    if size(Sout,1) != size(Sout,2)
        throw(DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."))
    end

    if size(Sa,1) + size(Sb,1) - 2 != size(Sout,1)
        throw(DimensionMismatch("First two dimensions of `Sout` must be `m+n-2`."))
    end

    for i in 3:ndims(Sa)
        if size(Sa,i) != size(Sout,i)
            throw(DimensionMismatch("Non-port axis lengths of `Sa`, `Sb`, and `Sout` must be equal."))
        end
    end

    if k > size(Sa,1)
        throw(ArgumentError("Port `k` is larger than number of ports in `Sa`."))
    end

    if l > size(Sb,1)
        throw(ArgumentError("Port `l` is larger than number of ports in `Sb`."))
    end

    if l < 1
        throw(ArgumentError("Port `l` is smaller than one."))
    end

    if k < 1
        throw(ArgumentError("Port `k` is smaller than one."))
    end

    # loop over the dimensions of the array greater than 2
    indices = CartesianIndices(axes(Sout)[3:end])
    if nbatches > 1 && length(indices) > nbatches
        batches = Base.Iterators.partition(1:length(indices),1+(length(indices)-1)÷nbatches)
        Threads.@sync for batch in batches
            Base.Threads.@spawn interconnectS_inner!(Sout,Sa,Sb,k,l,batch)
        end

    else
        interconnectS_inner!(Sout,Sa,Sb,k,l,indices)
    end

    return Sout
end

"""
    interconnectS_inner!(Sout,Sa,Sb,k::Int,l::Int,batch::AbstractArray)

See [`interconnectS`](@ref) for description.

"""
function interconnectS_inner!(Sout, Sa, Sb, k::Int, l::Int, batch::AbstractArray)

    # the number of ports in the input matrix
    m = size(Sa,1)
    n = size(Sb,1)

    range1a = 1:k-1
    range1b = k+1:m
    range2a = m+1:m+l-1
    range2b = m+l+1:m+n

    # this indexes across the entire array
    ranges = (range1a, range1b, range2a, range2b)

    # this indexes across the first part
    ranges1 = (range1a, range1b)

    # this indexes across the second part
    ranges2 = (range2a, range2b)

    # loop over the axes of the scattering parameter matrices after the first
    # two (eg. frequencies).
    @inbounds for b in batch

        gammacc_Scc = StaticArrays.SMatrix{2,2}(
            -Sa[k,k,b],
            one(Sa[k,k,b]),
            one(Sa[k,k,b]),
            -Sb[l,l,b]
        )
        # gammacc_Scc_lu = lu(gammacc_Scc)
        gammacc_Scc_lu = lu_2x2(gammacc_Scc)

        # ii and jj are the indices which extend up to m+n and skip k,l
        # i and j extend up to m+n-2 and are consecutive 

        # left side
        for jindex in eachindex(ranges1)

            for jj in ranges1[jindex]
                j = jj-jindex+1
                Scp = StaticArrays.SVector{2}(Sa[k,jj,b],zero(Sa[k,jj,b]))

                # solve the linear system
                # ac1jj, ac2jj = gammacc_Scc_lu \ Scp
                ac1jj, ac2jj = ldiv_2x2(gammacc_Scc_lu,Scp)

                # upper left quadrant
                for iindex in eachindex(ranges1)
                    for ii in ranges1[iindex]
                        i = ii-iindex+1
                        Sout[i,j,b] = Sa[ii,jj,b] + Sa[ii,k,b]*ac1jj
                    end
                end

                # lower left quadrant
                for iindex in eachindex(ranges2)
                    for ii in ranges2[iindex]
                        i = ii-iindex+1-1
                        Sout[i,j,b] = Sb[ii-m,l,b]*ac2jj
                    end
                end

            end

        end

        # right side
        for jindex in eachindex(ranges2)

            for jj in ranges2[jindex]
                j = jj-jindex+1-1
                Scp = StaticArrays.SVector{2}(zero(Sb[l,jj-m,b]),Sb[l,jj-m,b])

                # solve the linear system.
                # ac1jj, ac2jj = gammacc_Scc_lu \ Scp
                ac1jj, ac2jj = ldiv_2x2(gammacc_Scc_lu,Scp)

                # upper right quadrant
                for iindex in eachindex(ranges1)
                    for ii in ranges1[iindex]
                        i = ii-iindex+1
                        Sout[i,j,b] =  Sa[ii,k,b]*ac1jj
                    end
                end

                # lower right quadrant
                for iindex in eachindex(ranges2)
                    for ii in ranges2[iindex]
                        i = ii-iindex+1-1
                        Sout[i,j,b] = Sb[ii-m,jj-m,b] + Sb[ii-m,l,b]*ac2jj
                    end
                end

            end

        end

    end

    return Sout

end

"""
    intraconnectSports(portsa::AbstractVector{Tuple{T,Int}},k::Int,l::Int) where T

Return a vector of tuples of (networkname, portindex) from `portsa` after
ports `k` and `l` have been connected. See [`connectS`](@ref) for more
information.

# Examples
```jldoctest
julia> JosephsonCircuits.intraconnectSports([(:S1,1),(:S1,2),(:S1,3),(:S1,4),(:S1,5)],3,4)
3-element Vector{Tuple{Symbol, Int64}}:
 (:S1, 1)
 (:S1, 2)
 (:S1, 5)
```
"""
function intraconnectSports(portsa::AbstractVector{Tuple{T,Int}}, k::Int,
    l::Int) where T


    if k > length(portsa)
        throw(ArgumentError("Port `k` is larger than number of ports."))
    end

    if l > length(portsa)
        throw(ArgumentError("Port `l` is larger than number of ports."))
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

    # the ports of the network after the connections
    ports = similar(portsa,length(portsa)-2)

    # loop over the ports
    j = 0
    for i in eachindex(ports)
        j += 1
        # skip over the ports that are connected together
        if j == k || j == l
            j += 1
        end

        if j == k || j == l
            j += 1
        end

        ports[i] = portsa[j]
    end

    return ports
end

"""
    interconnectSports(portsa::AbstractVector{Tuple{T,Int}},
        portsb::AbstractVector{Tuple{T,Int}}, k::Int, l::Int) where T

Return a vector of tuples of (networkname, portindex) with `portsa` from the
first network and `portsb` from the second network after ports `k` and `l`
from the first and second networks have been connected. If the first network
has `n` ports and the second network has `m` ports, then the combined network
has `(m+n-2)` ports. See [`connectS`](@ref) for more information.

# Examples
```jldoctest
julia> JosephsonCircuits.interconnectSports([(:S1,1),(:S1,2),(:S1,3),(:S1,4),(:S1,5)],[(:S2,1),(:S2,2),(:S2,3),(:S2,4),(:S2,5)],3,4)
8-element Vector{Tuple{Symbol, Int64}}:
 (:S1, 1)
 (:S1, 2)
 (:S1, 4)
 (:S1, 5)
 (:S2, 1)
 (:S2, 2)
 (:S2, 3)
 (:S2, 5)
```
"""
function interconnectSports(portsa::AbstractVector{Tuple{T,Int}},
        portsb::AbstractVector{Tuple{T,Int}}, k::Int, l::Int) where T


    if k > length(portsa)
        throw(ArgumentError("Port `k` is larger than number of ports in `portsa`."))
    end

    if l > length(portsb)
        throw(ArgumentError("Port `l` is larger than number of ports in `portsb`."))
    end

    if l < 1
        throw(ArgumentError("Port `l` is smaller than one."))
    end

    if k < 1
        throw(ArgumentError("Port `k` is smaller than one."))
    end

    # the ports of the network after the connections
    ports = similar(portsa,length(portsa)+length(portsb)-2)


    # loop over the first network ports, skipping k
    j = 0
    for i in 1:length(portsa)-1
        j += 1
        # skip over the ports that are connected together
        if j == k
            j += 1
        end
        ports[i] = portsa[j]
    end

    # loop over the second network ports, skipping l
    j = 0
    for i in length(portsa):length(portsa)+length(portsb)-2
        j += 1
        # skip over the ports that are connected together
        if j == l
            j += 1
        end
        ports[i] = portsb[j]
    end

    return ports
end


"""
    cascadeS(Sa, Sb)

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
function cascadeS(Sa, Sb)

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
    S12 .= S12a*tmp2 # typo in ref: S11a should be S11b
    S21 .= S21b*tmp1
    S22 .= S22b + S21b*S22a*tmp2

    return nothing
end

"""
    remove_edge!(g,src_node,edge_index)

Remove an edge from graph `g` specified by the source node `src_node`
and the edge index `edge_index` in the forward adjacency list.

# Examples
```jldoctest
julia> g=JosephsonCircuits.Graphs.SimpleDiGraphFromIterator(JosephsonCircuits.tuple2edge([(1,1),(2,1),(2,3)]));JosephsonCircuits.remove_edge!(g,1,1)
1
```
"""
function remove_edge!(g, src_node, edge_index)

    # remove the source
    dst_node = popat!(g.fadjlist[src_node], edge_index)
    # and the destination
    for k in eachindex(g.badjlist[dst_node])
        if g.badjlist[dst_node][k] == src_node
            deleteat!(g.badjlist[dst_node],k)
            break
        end
    end
    return dst_node
end

"""
    move_fedge!(g,src_node,src_node_new,edge_index,fadjlist1,fadjlist2)

Move an edge from graph `g` at source node `src_node` to the new source node
`src_node_new` with the edge index `edge_index` in the forward adjacency list.
Also perform the same operations on the forward adjacency lists
`fadjlist1` and `fadjlist2`.

# Examples
```jldoctest
julia> g=JosephsonCircuits.Graphs.SimpleDiGraphFromIterator(JosephsonCircuits.tuple2edge([(1,1),(2,1),(2,3)]));JosephsonCircuits.move_fedge!(g,1,2,1,deepcopy(g.fadjlist),deepcopy(g.fadjlist));g.fadjlist
3-element Vector{Vector{Int64}}:
 []
 [1, 3, 1]
 []
```
"""
function move_fedge!(g, src_node, src_node_new, edge_index, fadjlist1,
    fadjlist2)

    # remove the source
    dst_node = popat!(g.fadjlist[src_node], edge_index)
    # add the new source
    push!(g.fadjlist[src_node_new], dst_node)

    if !isempty(fadjlist1)
        push!(fadjlist1[src_node_new], popat!(fadjlist1[src_node], edge_index))
    end

    if !isempty(fadjlist2)
        push!(fadjlist2[src_node_new], popat!(fadjlist2[src_node], edge_index))
    end

    # and update the destination
    for k in eachindex(g.badjlist[dst_node])
        if g.badjlist[dst_node][k] == src_node
            deleteat!(g.badjlist[dst_node],k)
            push!(g.badjlist[dst_node],src_node_new)
            break
        end
    end
    return dst_node
end

"""
    move_bedge!(g,dst_node,dst_node_new,edge_index,fadjlist1,fadjlist2)

Move an edge from graph `g` at destination node `dst_node` to the new
destination node `dst_node_new` with the edge index `edge_index` in the
backwards adjacency list. Also perform the same operations on the forward
adjacency lists `fadjlist1` and `fadjlist2`.

# Examples
```jldoctest
julia> g=JosephsonCircuits.Graphs.SimpleDiGraphFromIterator(JosephsonCircuits.tuple2edge([(1,1),(2,1),(2,3)]));JosephsonCircuits.move_bedge!(g,1,2,1,deepcopy(g.fadjlist),deepcopy(g.fadjlist));g.badjlist
3-element Vector{Vector{Int64}}:
 [2]
 [1]
 [2]
```
"""
function move_bedge!(g,dst_node,dst_node_new,edge_index,fadjlist1,fadjlist2)

    # remove the destination
    src_node = popat!(g.badjlist[dst_node], edge_index)
    # add the new destination
    push!(g.badjlist[dst_node_new],src_node)

    # and update the source
    for k in eachindex(g.fadjlist[src_node])
        if g.fadjlist[src_node][k] == dst_node
            deleteat!(g.fadjlist[src_node],k)
            push!(g.fadjlist[src_node], dst_node_new)

            if !isempty(fadjlist1)
                push!(fadjlist1[src_node], popat!(fadjlist1[src_node], k))
            end

            if !isempty(fadjlist2)
                push!(fadjlist2[src_node], popat!(fadjlist2[src_node], k))
            end

            break
        end
    end
    return src_node
end

"""
    move_fedges!(g,src_node,src_node_new,fadjlist1,fadjlist2)

Move the edges from graph `g` at source node `src_node` to the new source node
`src_node_new` in the forward adjacency list. Also perform the same operations
on the forward adjacency lists `fadjlist1` and `fadjlist2`.

# Examples
```jldoctest
julia> g=JosephsonCircuits.Graphs.SimpleDiGraphFromIterator(JosephsonCircuits.tuple2edge([(1,1),(2,1),(2,3)]));JosephsonCircuits.move_fedges!(g,1,2,deepcopy(g.fadjlist),deepcopy(g.fadjlist));g.fadjlist
3-element Vector{Vector{Int64}}:
 []
 [1, 3, 1]
 []
```
"""
function move_fedges!(g,src_node,src_node_new,fadjlist1,fadjlist2)

    for edge_index in reverse(1:length(g.fadjlist[src_node]))
        move_fedge!(g,src_node,src_node_new,edge_index,fadjlist1,fadjlist2)
    end
    return g
end

"""
    move_bedges!(g,dst_node,dst_node_new,fadjlist1,fadjlist2)

Move the edges from graph `g` at destination node `dst_node` to the new
destination node `dst_node_new` in the backwards adjacency list. Also perform
the same operations on the forward adjacency lists `fadjlist1` and `fadjlist2`.

# Examples
```jldoctest
julia> g=JosephsonCircuits.Graphs.SimpleDiGraphFromIterator(JosephsonCircuits.tuple2edge([(1,1),(2,1),(2,3)]));JosephsonCircuits.move_bedges!(g,1,2,deepcopy(g.fadjlist),deepcopy(g.fadjlist));g.badjlist
3-element Vector{Vector{Int64}}:
 []
 [2, 1]
 [2]
```
"""
function move_bedges!(g,dst_node,dst_node_new,fadjlist1,fadjlist2)

    for edge_index in reverse(1:length(g.badjlist[dst_node]))
        move_bedge!(g,dst_node,dst_node_new,edge_index,fadjlist1,fadjlist2)
    end
    return g
end

"""
    move_edges!(g,node,node_new,fadjlist1,fadjlist2)

Move the edges from graph `g` at node `node` to the new node `node_new`. Also
perform the same operations on the forward adjacency lists `fadjlist1` and
`fadjlist2`.

# Examples
```jldoctest
julia> g=JosephsonCircuits.Graphs.SimpleDiGraphFromIterator(JosephsonCircuits.tuple2edge([(1,1),(2,1),(2,3)]));JosephsonCircuits.move_edges!(g,1,2,deepcopy(g.fadjlist),deepcopy(g.fadjlist));g.fadjlist
3-element Vector{Vector{Int64}}:
 []
 [3, 2, 2]
 []
```
"""
function move_edges!(g,node,node_new,fadjlist1,fadjlist2)
    move_fedges!(g,node,node_new,fadjlist1,fadjlist2)
    move_bedges!(g,node,node_new,fadjlist1,fadjlist2)
    return g
end

"""
    add_ports(networks)

Return the vector of networks `networks` with ports added.

# Examples
```jldoctest
julia> networks = [(:S1,[0.0 1.0;1.0 0.0]),(:S2,[0.5 0.5;0.5 0.5])];JosephsonCircuits.add_ports(networks)
2-element Vector{Tuple{Symbol, Matrix{Float64}, Vector{Tuple{Symbol, Int64}}}}:
 (:S1, [0.0 1.0; 1.0 0.0], [(:S1, 1), (:S1, 2)])
 (:S2, [0.5 0.5; 0.5 0.5], [(:S2, 1), (:S2, 2)])

julia> networks = [(:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S1,2)]),(:S2,[0.5 0.5;0.5 0.5],[(:S3,1),(:S3,2)])];JosephsonCircuits.add_ports(networks)
2-element Vector{Tuple{Symbol, Matrix{Float64}, Vector{Tuple{Symbol, Int64}}}}:
 (:S1, [0.0 1.0; 1.0 0.0], [(:S1, 1), (:S1, 2)])
 (:S2, [0.5 0.5; 0.5 0.5], [(:S3, 1), (:S3, 2)])
```
"""
function add_ports(networks::AbstractVector)
    return [(network[1],network[2],get_ports(network)) for network in networks]
end

function add_ports(networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}}) where {T,N}
    return networks
end

"""
    get_ports(network::Tuple{T, N}) where {T,N}

Return the ports for a network `network`. The ports are generated based on
the network name.

# Examples
```jldoctest
julia> JosephsonCircuits.get_ports((:S1,[0.0 1.0;1.0 0.0]))
2-element Vector{Tuple{Symbol, Int64}}:
 (:S1, 1)
 (:S1, 2)
```
"""
function get_ports(network::Tuple{T, N}) where {T,N}
    # a vector of tuples containing the ports for each of the networks in the
    # same order the networks were supplied. eg.
    # [(:S1,1),(:S1,2)]
    return [(network[1],i) for i in 1:size(network[2],1)]
end

"""
    get_ports(network::Tuple{T, N, Vector{Tuple{T, Int}}}) where {T,N}

Return the ports for a network `network`. The ports are already present in the
network.

# Examples
```jldoctest
julia> JosephsonCircuits.get_ports((:S1,[0.0 1.0;1.0 0.0],[(:S1,1),(:S2,3)]))
2-element Vector{Tuple{Symbol, Int64}}:
 (:S1, 1)
 (:S2, 3)
```
"""
function get_ports(network::Tuple{T, N, Vector{Tuple{T, Int}}}) where {T,N}
    return network[3]
end

"""
    find_duplicate_network_names(
        networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}}) where {T,N}

Return a vector of tuples of (networkname, counts) where counts is the number
of times a given network name appears.

"""
function find_duplicate_network_names(
    networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}}) where {T,N}
    # find the duplicates to return a useful error message
    networkdatacounts = Dict{T,Int}()
    # count the number of times each network occurs
    for network in networks
        if haskey(networkdatacounts,network[1])
            networkdatacounts[network[1]] += 1
        else
            networkdatacounts[network[1]] = 1
        end
    end
    # report any networks that occur more than once
    networkdataduplicates = Tuple{T,Int}[]
    for (key,val) in networkdatacounts
        if val > 1
            push!(networkdataduplicates,(key,val))
        end
    end
    return networkdataduplicates
end

"""
    find_duplicate_connections(
        connections::AbstractVector{Tuple{T,T,Int,Int}}) where {T}

Return a vector of tuples of (connection, counts) where counts is the number
of times a given connection appears.
"""
function find_duplicate_connections(
    connections::AbstractVector{Tuple{T,T,Int,Int}}) where {T}
    # find the duplicates to return a useful error message
    connectionscounts = Dict{Tuple{T,Int},Int}()
    # count the number of times each network occurs
    for connection in connections
        key1 = (connection[1],connection[3])
        if haskey(connectionscounts,key1)
            connectionscounts[key1] += 1
        else
            connectionscounts[key1] = 1
        end
        key2 = (connection[2],connection[4])
        if haskey(connectionscounts,key2)
            connectionscounts[key2] += 1
        else
            connectionscounts[key2] = 1
        end
    end
    # report any networks that occur more than once
    connectionsduplicates = Tuple{Tuple{T,Int},Int}[]
    for (key,val) in connectionscounts
        if val > 1
            push!(connectionsduplicates,(key,val))
        end
    end
    return connectionsduplicates
end

"""
    add_modes(connections::AbstractVector{Tuple{T,T,Int,Int}},
        Nmodes::Integer) where {T}

Assume the scattering parameter matrices are multi-mode and `connections`
specifies the connections between physical ports. Return the connections
vector with the existing connections re-numbered and with added connections
between every mode associated with each physical port.

# Examples
```jldoctest
julia> connections = [("S1","S2",3,4),("S1","S2",2,1)];JosephsonCircuits.add_modes(connections,2)
4-element Vector{Tuple{String, String, Int64, Int64}}:
 ("S1", "S2", 5, 7)
 ("S1", "S2", 6, 8)
 ("S1", "S2", 3, 1)
 ("S1", "S2", 4, 2)
```
"""
function add_modes(connections::AbstractVector{Tuple{T,T,Int,Int}},
    Nmodes::Integer) where {T}

    # initialize an empty vector for the connections that is Nmodes times
    # longer than the input
    connections_modes = Vector{Tuple{T,T,Int,Int}}(undef,Nmodes*length(connections))
    
    # loop over the connections and for each connection between two ports, add
    # a connection between each of the same two ports for each mode. adjust
    # the indices of the ports to follow the all modes for the first port,
    # then all the modes for the next port convention.
    for k in eachindex(connections)
        (name1,name2,port1,port2) = connections[k]
        for i in 1:Nmodes
            connections_modes[(k-1)*Nmodes+i] = (name1,name2,(port1-1)*Nmodes+i,(port2-1)*Nmodes+i)
        end
    end
    return connections_modes
end

"""
    add_modes(connections::AbstractVector{Vector{Tuple{T,Int}}},
        Nmodes::Integer) where {T}

Assume the scattering parameter matrices are multi-mode and `connections`
specifies the connections between physical ports. Return the connections
vector with the existing connections re-numbered and with added connections
between every mode associated with each physical port.

# Examples
```jldoctest
julia> connections = [[("S1",3),("S2",4)],[("S1",2),("S2",1)]];JosephsonCircuits.add_modes(connections,2)
4-element Vector{Vector{Tuple{String, Int64}}}:
 [("S1", 5), ("S2", 7)]
 [("S1", 6), ("S2", 8)]
 [("S1", 3), ("S2", 1)]
 [("S1", 4), ("S2", 2)]
```
"""
function add_modes(connections::AbstractVector{Vector{Tuple{T,Int}}},
    Nmodes::Integer) where {T}

    # initialize an empty vector for the connections that is Nmodes times
    # longer than the input
    connections_modes = Vector{Vector{Tuple{T,Int}}}(undef,Nmodes*length(connections))
    
    # loop over the connections and for each connection between two ports, add
    # a connection between each of the same two ports for each mode. adjust
    # the indices of the ports to follow the all modes for the first port,
    # then all the modes for the next port convention.
    for k in eachindex(connections)
        for i in 1:Nmodes
            connections_modes[(k-1)*Nmodes+i] = [(name,(port-1)*Nmodes+i) for (name,port) in connections[k]]
        end
    end
    return connections_modes
end

"""
    add_splitters(networks::AbstractVector{Tuple{T,N}},
        connections::AbstractVector{<:AbstractVector{Tuple{T,Int}}};
        splitter_name_length = 20) where {T,N}

Return the `networks` and `connections` with splitters (ideal lossless
symmetrical reciprocal networks) and connections to the splitters added when
more than two ports intersect. `connections` is also converted from a vector
of vectors of tuples where the tuple contains the network and the port such as
[[(:S1,1),(:S2,1)]] to a vector of tuples where the tuple contains the two
networks and ports being connected [(S1,:S2,1,1)].

# References
S. F. Cao, Y. C. Jiao, and Z. Zhang. "Applications of Generalized Cascade
Scattering Matrix on the Microwave Circuits and Antenna Arrays". International
Journal of Antennas and Propagation Vol. 2015, 759439,
doi:10.1155/2015/759439.
"""
function add_splitters(networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}},
    connections::AbstractVector{Vector{Tuple{T,Int}}};
    small_splitters = true) where {T,N}

    # compute the size and element type of the first scattering matrix and
    # check the rest against this
    sizeS = size(networks[1][2])[3:end]
    typeS = eltype(networks[1][2])
    for network in networks
        if size(network[2],1) != size(network[2],2)
            throw(ArgumentError("The sizes of the first two dimensions ($(size(network[2],1)),$(size(network[2],2))) of the scattering matrix $(network[1]) must be the same."))
        end
        if sizeS != size(network[2])[3:end]
            throw(ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of $(networks[1][1]) is $(size(networks[1][2])) and size of $(network[1]) is $(size(network[2]))."))
        end
        if typeS != eltype(network[2])
            throw(ArgumentError("The element types of the scattering matrices must be the same. Element type of $(networks[1][1]) is $(eltype(networks[1][2])) and element type of $(network[1]) is $(eltype(network[2]))."))
        end
    end

    # copy the networks vector. don't deepcopy so we don't duplicate all of
    # the arrays contained in the vector.
    netflat = copy(networks)

    # define a new vector to store the connections with splitters added.
    conflat = Tuple{T,T,Int,Int}[]

    # Store the scattering parameter matrices of the splitters, so we can
    # reuse the matrices if the same splitter is used twice. The key is the
    # size the scattering matrix and the value is the matrix itself.
    splitters = Dict{Int,N}()

    # loop over the connections, converting to the flattened format and adding
    # splitters where more than two ports are connected.
    # for c in connections
    for k in eachindex(connections)
        c = connections[k]
        # if less than two tuples, not a valid connection
        if length(c) < 2
            throw(ArgumentError("Invalid connection $(c) with only network and port."))
        # if two tuples, this is a single connection
        elseif length(c) == 2
            push!(conflat,(c[1][1],c[2][1],c[1][2],c[2][2]))
        # if more than two tuples, add a splitter
        # and make all of the connections
        else
            # if small_splitters is true, then generate the N port splitter
            # by combining (N-2) 3 port splitters. if small_splitters is
            # false, then make the N port splitter and connect the components
            # to it.
            if small_splitters
                Nsplitters = length(c)-2

                # make all of the splitters
                for i in 1:Nsplitters

                    # make a new name for the splitter
                    id = T(string(UUIDs.uuid1()))

                    # compute the size of the splitter
                    # assume we will always make a 3 port splitter
                    sizeS = NTuple{ndims(netflat[1][2])}(ifelse(j<=2,3,size(netflat[1][2],j)) for j in 1:ndims(netflat[1][2]))

                    # check if we have already made this splitter and make a new
                    # one if not.
                    if !haskey(splitters,2)
                        splitters[2] = S_splitter!(similar(N,sizeS))
                    end

                    # add the splitter to the vector of networks
                    push!(netflat,(id,splitters[2],get_ports((id,splitters[2]))))

                end

                # always make the first two connections between the splitter and
                # the components
                # first port of first splitter to first component
                push!(conflat,(netflat[end-Nsplitters+1][1],c[1][1],1,c[1][2]))

                # second port of first splitter to second component
                push!(conflat,(netflat[end-Nsplitters+1][1],c[2][1],2,c[2][2]))

                # if only one splitter make a connection between the third port
                # of the splitter to the third component
                if Nsplitters == 1
                    push!(conflat,(netflat[end-Nsplitters+1][1],c[3][1],3,c[3][2]))
                # if more than one splitter, then connect the third port of the
                # splitter to the first port of the next splitter
                else
                    push!(conflat,(netflat[end-Nsplitters+1][1],netflat[end-Nsplitters+2][1],3,1))
                end

                # loop through the rest of the splitters
                # add groups of two connections
                for i in 2:Nsplitters
                    # first is always between the second port of the splitter
                    # and a component
                    push!(conflat,(netflat[end-Nsplitters+i][1],c[i+1][1],2,c[i+1][2]))

                    if i == Nsplitters
                        push!(conflat,(netflat[end-Nsplitters+i][1],c[i+2][1],3,c[i+2][2]))
                    else
                        push!(conflat,(netflat[end-Nsplitters+i][1],netflat[end-Nsplitters+i+1][1],3,1))
                    end
                end

            else
                # make a new name
                id = T(string(UUIDs.uuid1()))

                # compute the size of the splitter
                sizeS = NTuple{ndims(netflat[1][2])}(ifelse(j<=2,length(c),size(netflat[1][2],j)) for j in 1:ndims(netflat[1][2]))

                # check if we have already made this splitter and make a new one
                # if not.
                if !haskey(splitters,length(c))
                    splitters[length(c)] = S_splitter!(similar(N,sizeS))
                end

                # add the splitter to the vector of networks
                push!(netflat,(id,splitters[length(c)],get_ports((id,splitters[length(c)]))))
        
                # add the connections to the splitter
                for j in eachindex(c)
                    push!(conflat,(id,c[j][1],j,c[j][2]))
                end
            end
        end
    end

    return netflat, conflat
end

"""
    add_splitters(networks, connections::AbstractVector{Tuple{T,T,Int,Int}};
        kwargs...)) where T

If the connections are already in the correct format, just return them. This
function assumes ports have already been added to `networks`.
"""
function add_splitters(networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}},
    connections::AbstractVector{Tuple{T,T,Int,Int}};
    small_splitters = true) where {T,N}
    return networks,connections
end

"""
    make_connection!(g,fconnectionlist,fweightlist,ports,networkdata,src_node,
    connection_index)

Apply the connection specified by the source node `src_node` and the index
of the connection in the forward adjacency list `connection_index`. Modify the
arguments and return nothing.

"""
function make_connection!(g::Graphs.SimpleGraphs.SimpleDiGraph{Int},
    fconnectionlist::AbstractVector{<:AbstractVector{Tuple{T,T,Int,Int}}},
    fweightlist::AbstractVector{<:AbstractVector{Int}},
    ports::AbstractVector{<:AbstractVector{Tuple{T,Int}}},
    networkdata::AbstractVector{N},
    src_node::Int,connection_index::Int,nbatches::Int,
    userinput::AbstractVector{Bool},
    storage::Dict{Int,N}) where {T,N}

    # remove the edge from the graph
    dst_node = remove_edge!(g,src_node,connection_index)

    # remove and store the connection associated with that edge
    connection = popat!(fconnectionlist[src_node], connection_index)

    # remove the weight for that connection
    weight = popat!(fweightlist[src_node], connection_index)

    # the source and destination ports eg. (:S1,1) and (:S2,1)
    src_port = (connection[1],connection[3])
    dst_port = (connection[2],connection[4])

    # the indices at which these ports are located in the source
    # destination scattering parameter matrix
    src_port_index = findfirst(isequal(src_port),ports[src_node])
    dst_port_index = findfirst(isequal(dst_port),ports[dst_node])

    if isnothing(src_port_index)
        throw(ArgumentError("Source port $(src_port) not found in the ports $(ports[src_node]) of the source node $(src_node)."))
    end
    if isnothing(dst_port_index)
        throw(ArgumentError("Destination port $(dst_port) not found in the ports $(ports[dst_node]) of the destination node $(dst_node)."))
    end
    # println("src_node => dst_node: ",src_node," => ",dst_node)
    # println("src_port => dst_port: ",src_port," => ",dst_port)
    # println("src_port_index => dst_port_index: ",src_port_index," => ",dst_port_index)

    # if src_node == dst_node, then make a self connection
    if src_node == dst_node
        # connect the networks and find the ports of the connected network
        connected_network = intraconnectS(networkdata[src_node],src_port_index,dst_port_index;nbatches=nbatches)
        connected_ports = intraconnectSports(ports[src_node],src_port_index,dst_port_index)
    
        # delete the ports for the dst. update the ports for the src.
        ports[src_node] = connected_ports

        # update the networkdata for the src and replace the src with an empty array.
        networkdata[src_node] = connected_network
    else

        # when making a new array. try to see if that same sizes is in storage.
        # if so, reuse that. if not, make a new one and add to storage.
        # then switch from connectS to connectS!
        d = size(networkdata[src_node],1)+size(networkdata[dst_node],1)-2
        if haskey(storage,d)
            connected_network = pop!(storage,d)
        else
            # make a tuple with the size of the array
            # the first two dimensions are two smaller
            sizeSx = size(networkdata[src_node])
            sizeSy = size(networkdata[dst_node])
            sizeS = NTuple{length(sizeSx)}(ifelse(i<=2,sizeSx[i]+sizeSy[i]-2,sizeSx[i]) for i in 1:length(sizeSx))

            # allocate an array of zeros of the same type as Sx
            connected_network = similar(networkdata[src_node],sizeS)
        end

        # connect the networks and find the ports of the connected network
         interconnectS!(connected_network,networkdata[src_node],networkdata[dst_node],src_port_index,dst_port_index;nbatches=nbatches)
        connected_ports = interconnectSports(ports[src_node],ports[dst_node],src_port_index,dst_port_index)

        # delete the ports for the dst. update the ports for the src.
        ports[src_node] = connected_ports
        ports[dst_node] = []


        # if the nodes are not user input, then store the arrays
        if !userinput[src_node]
            storage[size(networkdata[src_node],1)]= networkdata[src_node]
        end
        if !userinput[dst_node]
            storage[size(networkdata[dst_node],1)]= networkdata[dst_node]
        end

        # update the networkdata for the src and replace the src with an empty array.
        networkdata[src_node] = connected_network
        networkdata[dst_node] = Array{eltype(connected_network)}(undef,ntuple(zero,ndims(connected_network)))
        
        # mark the source and destination nodes as non-user input
        # so we can re-use the arrays if desired
        userinput[src_node] = false
        userinput[dst_node] = false

        # move the edges away from the dst node
        move_edges!(g,dst_node,src_node,fconnectionlist,fweightlist)
    end

    # applying a connection changes the size of the scattering parameter
    # matrices. first update the weights for the connections originating
    # from the src_node. fadjlist[src_node] has the indices of the nodes
    # for each connection originating there, so loop over those, and
    # update the weights
    for k in eachindex(g.fadjlist[src_node])
        # update the weights of fweightlist[src_node][k]
        if src_node == g.fadjlist[src_node][k]
            # self connections always reduce the size so give them zero weight
            fweightlist[src_node][k] = 0
        else
            src_weight = size(networkdata[src_node],1)-1
            dst_weight = size(networkdata[g.fadjlist[src_node][k]],1)-1
            connection_weight = src_weight*dst_weight
            fweightlist[src_node][k] = connection_weight
        end
    end
    # also update the weights for any connection that ends at the
    # src_node. those nodes are found in g.badjlist[src_node]
    for k in g.badjlist[src_node]
        for l in eachindex(g.fadjlist[k])
            if g.fadjlist[k][l] == src_node
                 # update the weights of fweightlist[k][l]
                if src_node == k
                    # self connections always reduce the size so give them zero weight
                    fweightlist[k][l] = 0
                else
                    src_weight = size(networkdata[src_node],1)-1
                    dst_weight = size(networkdata[k],1)-1
                    connection_weight = src_weight*dst_weight
                    fweightlist[k][l] = connection_weight
                end
            end
        end
    end
    return nothing
end

"""
    connectS_initialize(networks::AbstractVector, connections::AbstractVector;
    small_splitters::Bool = true)

Return a directed graph of connections between the networks.

# Examples
```jldoctest
networks = [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5])];
connections = [[("S1",1),("S2",2)]];
JosephsonCircuits.connectS_initialize(networks,connections)

# output
(Graphs.SimpleGraphs.SimpleDiGraph{Int64}(2, [[2], Int64[]], [Int64[], [1]]), [[("S1", "S2", 1, 2)], Tuple{String, String, Int64, Int64}[]], [[1], Int64[]], [[("S1", 1), ("S1", 2)], [("S2", 1), ("S2", 2)]], [[0.0 1.0; 1.0 0.0], [0.5 0.5; 0.5 0.5]])
```
"""
function connectS_initialize(networks::AbstractVector, connections::AbstractVector;
    small_splitters::Bool = true, Nmodes::Integer = 1) 

    networks_ports = add_ports(networks)

    connections_modes = add_modes(connections, Nmodes)

    networks_flat, connnections_flat = add_splitters(networks_ports,
        connections_modes; small_splitters = small_splitters)

    return connectS_initialize(networks_flat, connnections_flat)
end

"""
    connectS_initialize(networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}},
        connections::AbstractVector{Tuple{T,T,Int,Int}}) where {T,N}

Return a directed graph of connections between the networks.

# Examples
```jldoctest
networks = [("S1", [0.0 1.0; 1.0 0.0], [("S1", 1), ("S1", 2)]), ("S2", [0.5 0.5; 0.5 0.5], [("S2", 1), ("S2", 2)])];
connections = [("S1","S2",1,2)];
JosephsonCircuits.connectS_initialize(networks,connections)

# output
(Graphs.SimpleGraphs.SimpleDiGraph{Int64}(2, [[2], Int64[]], [Int64[], [1]]), [[("S1", "S2", 1, 2)], Tuple{String, String, Int64, Int64}[]], [[1], Int64[]], [[("S1", 1), ("S1", 2)], [("S2", 1), ("S2", 2)]], [[0.0 1.0; 1.0 0.0], [0.5 0.5; 0.5 0.5]])
```
"""
function connectS_initialize(networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}},
    connections::AbstractVector{Tuple{T,T,Int,Int}}) where {T,N}

    # network data is associated with each node, so also store those as a
    # vector of matrices where the index is the node index.
    networkdata = [network[2] for network in networks]

    # compute the size and element type of the first scattering matrix and
    # check the rest against this
    sizeS = size(networkdata[1])[3:end]
    typeS = eltype(networkdata[1])
    for i in eachindex(networkdata)
        if size(networkdata[i],1) != size(networkdata[i],2)
            throw(ArgumentError("The sizes of the first two dimensions ($(size(networkdata[i],1)),$(size(networkdata[i],2))) of the scattering matrix $(networks[i][1]) must be the same."))
        end
        if sizeS != size(networkdata[i])[3:end]
            throw(ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of $(networks[1][1]) is $(size(networks[1][2])) and size of $(networks[i][1]) is $(size(networks[i][2]))."))
        end
        if typeS != eltype(networkdata[i])
            throw(ArgumentError("The element types of the scattering matrices must be the same. Element type of $(networks[1][1]) is $(eltype(networks[1][2])) and element type of $(networks[i][1]) is $(eltype(networks[i][2]))."))
        end
    end

    # make a dictionary where the keys are the network names and the values
    # are the node indices.
    networkindices = Dict(network[1]=>i for (i,network) in enumerate(networks))

    # check if there are duplicate networks
    if length(networkindices) != length(networkdata)
        throw(ArgumentError("Duplicate network names detected [(networkname,count)]: $(find_duplicate_network_names(networks))."))
    end

    # check if the port indices are unique
    duplicateconnections = find_duplicate_connections(connections)
    if !isempty(duplicateconnections)
        throw(ArgumentError("Duplicate connections detected [(networkname,port),counts]: $(duplicateconnections)."))
    end

    # portnames are associated with each node, so store those as a vector
    # where the index is the node index
    ports = [network[3] for network in networks]

    # make the port dictionary
    portdict = Dict{Tuple{T,Int},Tuple{Int,Int}}()
    for i in eachindex(networks)
        for (j,port) in enumerate(networks[i][3])
            if haskey(portdict,port)
                throw(ArgumentError("Duplicate port $(port) in network $(networks[i][1])."))
            else
                portdict[port] = (i,j)
            end
        end
    end

    # make the adjacency lists for the connections
    fadjlist = Vector{Vector{Int}}(undef,length(networks))
    badjlist = Vector{Vector{Int}}(undef,length(networks))
    fconnectionlist = Vector{Vector{Tuple{T, T, Int64, Int64}}}(undef,length(networks))
    fweightlist = Vector{Vector{Int}}(undef,length(networks))


    # fill them with empty vectors
    for i in eachindex(fadjlist)
        fadjlist[i] = []
    end
    for i in eachindex(badjlist)
        badjlist[i] = []
    end
    for i in eachindex(fconnectionlist)
        fconnectionlist[i] = []
    end
    for i in eachindex(fweightlist)
        fweightlist[i] = []
    end

    # loop through the connections and populate the adjacency lists
    for (src_name, dst_name, src_port, dst_port) in connections

        src = (src_name,src_port)
        dst = (dst_name,dst_port)

        # check if the source and destination networks exist
        if !haskey(portdict,src)
            throw(ArgumentError("Source (network name, port number) ($(src_name), $(src_port)) not found for connection ($(src_name),$(dst_name),$(src_port),$(dst_port))."))
        end
        if !haskey(portdict,dst)
            throw(ArgumentError("Destination (network name, port number) ($(dst_name), $(dst_port)) not found for connection ($(src_name),$(dst_name),$(src_port),$(dst_port))."))
        end

        src_index, src_port_index = portdict[src]
        dst_index, dst_port_index = portdict[dst]

        # the source node entry points to the destination node 
        push!(fadjlist[src_index],dst_index)

        # the destination node entry points to the source node
        push!(badjlist[dst_index],src_index)

        # only store the source connections
        push!(fconnectionlist[src_index],(src_name, dst_name, src_port, dst_port))

        src_weight = size(networkdata[src_index],1)-1
        dst_weight = size(networkdata[dst_index],1)-1
        connection_weight = src_weight*dst_weight
        if src_index == dst_index
            connection_weight = 0
        end
        push!(fweightlist[src_index],connection_weight)

    end

    # turn the adjacency lists into a graph
    g = Graphs.SimpleGraphs.SimpleDiGraph(length(networkdata), fadjlist,badjlist)

    return (g, fconnectionlist, fweightlist, ports, networkdata)
end

"""
    connectS!(g::Graphs.SimpleGraphs.SimpleDiGraph{Int},
        fconnectionlist::AbstractVector{<:AbstractVector{Tuple{T,T,Int,Int}}},
        fweightlist::AbstractVector{<:AbstractVector{Int}},
        ports::AbstractVector{<:AbstractVector{Tuple{T,Int}}},
        networkdata::AbstractVector{N};
        nbatches::Int = Base.Threads.nthreads()) where {T,N}

Return the non-empty elements of the updated `networkdata` and `ports` after
applying all of the connections in the connection forward adjacency list
`fconnectionlist` to the graph `g`, the forward adjacency weight list
`fweightlist`, the vector of ports `ports`, and the vector of scattering
parameter matrices `networkdata`.

# Examples
```jldoctest
networks = [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5])];
connections = [[("S1",1),("S2",2)]];
init = JosephsonCircuits.connectS_initialize(networks, connections);
JosephsonCircuits.connectS!(init...)

# output
(S = [[0.5 0.5; 0.5 0.5]], ports = [[("S1", 2), ("S2", 1)]])
```
"""
function connectS!(g::Graphs.SimpleGraphs.SimpleDiGraph{Int},
    fconnectionlist::AbstractVector{<:AbstractVector{Tuple{T,T,Int,Int}}},
    fweightlist::AbstractVector{<:AbstractVector{Int}},
    ports::AbstractVector{<:AbstractVector{Tuple{T,Int}}},
    networkdata::AbstractVector{N};
    nbatches::Int = Base.Threads.nthreads()) where {T,N}

    # copy the graph, fconnectionlist, fweightlist, ports, and networkdata
    # so we don't mutate these, so we can apply connectS! multiple times only
    # modifying the input networks.
    g = deepcopy(g)
    fconnectionlist = deepcopy(fconnectionlist)
    fweightlist = deepcopy(fweightlist)
    ports = deepcopy(ports)
    # we don't modify the scattering parameter matrices that make up networkdata so
    # there is no need to deepycopy networkdata.
    networkdata = copy(networkdata)


    userinput = ones(Bool,length(networkdata))
    storage = Dict{Int,N}()
    # find the minimum weight and the second to minimum weight
    # we want unique weights, eg, both shouldn't be the same weight
    minweight = Inf
    secondtominweight = Inf
    for i in eachindex(fweightlist)
        for j in eachindex(fweightlist[i])
            weight = fweightlist[i][j]
            if weight < minweight
                # this is the new minimum weight and the old
                # minimum weight is now the second to minimum weight
                secondtominweight = minweight
                minweight = weight
            elseif weight > minweight && weight < secondtominweight
                # this is the new second minimum weight that is greater
                # than the minweight but less than the secondtominweight
                secondtominweight = weight
            end
        end
    end
    # println("minweight ",minweight)
    while !all(isempty.(fweightlist))
        for i in eachindex(fweightlist)
            j = 1
            n = length(fweightlist[i]) 
            while j <= n
                weight = fweightlist[i][j]
                if weight <= minweight
                    # perform the connection if this is the minimum weight
                    # println("i ",i," j ",j," N ",N," length(fweightlist[i]) ",length(fweightlist[i]))
                    # println(fweightlist[i])
                    make_connection!(g,fconnectionlist,fweightlist,ports,networkdata,i,j,nbatches,userinput,storage)
                    # set j = 1 to start looping through again
                    j = 1
                    n = length(fweightlist[i])
                    # we want to make all of the minimum weight connections
                    minweight = weight
                else
                    if weight < secondtominweight
                        # record the second to minimum weight
                        # which is greater than the minweight and less
                        # than the secondtominweight
                        secondtominweight = weight
                    end
                    j+=1
                end
            end
        end
        minweight = secondtominweight
        secondtominweight = Inf
    end
    return (S=networkdata[map(!isempty,networkdata)],ports=ports[map(!isempty,networkdata)])
end

"""
    connectS(networks, connections; small_splitters::Bool = true,
        nbatches::Int = Base.Threads.nthreads())

Return the network and ports resulting from connecting the networks in
`networks` according to the connections in `connections`. `networks` is a
vector of tuples of the network name and scattering parameter matrix such as
[("network1name",rand(Complex{Float64},2,2),
("network2name",rand(Complex{Float64},2,2)]. `connections` is a vector of
vectors of tuples of networks names and ports such as [[("network1name",1),
("network2name",2)]] where network1 and network2 are the two networks being
connected and 1 and 2 are integers describing the ports to connect.

This function supports connections between more than two ports by
automatically adding splitters.

# Examples
```jldoctest
networks = [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5])];
connections = [[("S1",1),("S2",2)]];
JosephsonCircuits.connectS(networks,connections)

# output
(S = [[0.5 0.5; 0.5 0.5]], ports = [[("S1", 2), ("S2", 1)]])
```
```jldoctest
networks = [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5],[("S3",5),("S3",6)])];
connections = [("S1","S3",1,6)];
JosephsonCircuits.connectS(networks,connections)

# output
(S = [[0.5 0.5; 0.5 0.5]], ports = [[("S1", 2), ("S3", 5)]])
```
"""
function connectS(networks, connections; small_splitters::Bool = true,
    Nmodes::Integer = 1, nbatches::Int = Base.Threads.nthreads())
    init = connectS_initialize(networks,connections;
        small_splitters = small_splitters, Nmodes = Nmodes)
    return connectS!(init...;nbatches = nbatches)
end

"""
    parse_connections_sparse(networks::AbstractVector{Tuple{T,N}},
        connections::AbstractVector{Tuple{T,T,Int,Int}}) where {T,N}

Return the indices of the internal ports `portc_indices`, the external ports
`portp_indices`, the vector of port tuples `ports`, the vector of scattering
parameter data `networkdata`, the connection matrix `gamma`, the sparse matrix
containing indices in networkdata `Sindices`, and an empty sparse matrix of
scattering parameter data `S`. The scattering parameter data consists of the
input networks assembled as a block diagonal matrix.

# References
V. A. Monaco and P. Tiberio, "Computer-Aided Analysis of Microwave Circuits,"
in IEEE Transactions on Microwave Theory and Techniques, vol. 22, no. 3, pp.
249-263, Mar. 1974, doi: 10.1109/TMTT.1974.1128208.
"""
function parse_connections_sparse(networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}},
    connections::AbstractVector{Tuple{T,T,Int,Int}}) where {T,N}

    # network data is associated with each node, so also store those as a
    # vector of matrices where the index is the node index.
    networkdata = [network[2] for network in networks]

    # first dimension of the final scattering matrix
    m = 0
    # number of nonzero elements in S
    M = 0
    # compute the size and element type of the first scattering matrix and
    # check the rest against this
    sizeS = size(networkdata[1])[3:end]
    typeS = eltype(networkdata[1])
    for i in eachindex(networkdata)
        if size(networkdata[i],1) != size(networkdata[i],2)
            throw(ArgumentError("The sizes of the first two dimensions ($(size(networkdata[i],1)),$(size(networkdata[i],2))) of the scattering matrix $(networks[i][1]) must be the same."))
        end
        if sizeS != size(networkdata[i])[3:end]
            throw(ArgumentError("The sizes of the third and higher dimensions of the scattering matrices must be the same. Size of $(networks[1][1]) is $(size(networks[1][2])) and size of $(networks[i][1]) is $(size(networks[i][2]))."))
        end
        if typeS != eltype(networkdata[i])
            throw(ArgumentError("The element types of the scattering matrices must be the same. Element type of $(networks[1][1]) is $(eltype(networks[1][2])) and element type of $(networks[i][1]) is $(eltype(networks[i][2]))."))
        end
        m+=size(networkdata[i],1)
        M+=size(networkdata[i],1)*size(networkdata[i],2)
    end

    # the indices in the sparse matrix formed by placing all of the scattering
    # matrices along the diagonal.
    networkdataindices = zeros(Int,length(networks))
    networkdataindices[1] = 1
    for i in 2:length(networkdataindices)
        networkdataindices[i] = networkdataindices[i-1] + size(networkdata[i-1],1)
    end

    # make a dictionary where the keys are the network names and the values
    # are the node indices.
    networkindices = Dict(network[1]=>i for (i,network) in enumerate(networks))

    if length(networkindices) != length(networkdata)
        throw(ArgumentError("Duplicate network names detected [(networkname,count)]: $(find_duplicate_network_names(networks))."))
    end

    # a vector of tuples containing the ports for each of the networks in the
    # same order the networks were supplied. eg.
    # [(:S1,1),(:S1,2),(:S2,1),(:S2,2)]
    ports = Vector{Tuple{T,Int}}(undef,m)

    # make the port dictionary
    portdict = Dict{Tuple{T,Int},Int}()

    k = 1
    for i in eachindex(networks)
        for port in networks[i][3]
            if haskey(portdict,port)
                throw(ArgumentError("Duplicate port $(port) in network $(networks[i][1])."))
            else
                portdict[port] = k
                ports[k] = port
            end
            k+=1
        end
    end

    # Compute the sparse connection matrix gamma and the port indices as which
    # connections occur which are also called internal ports.
    Igamma = Vector{Int}(undef,2*length(connections))
    Jgamma = Vector{Int}(undef,2*length(connections))
    Vgamma = Vector{Int}(undef,2*length(connections))
    empty!(Igamma)
    empty!(Jgamma)
    empty!(Vgamma)

    portc_indices = Vector{Int}(undef,2*length(connections))
    empty!(portc_indices)

    # loop through the connections and compute the sparse connection matrix
    # gamma
    for (src_name, dst_name, src_port, dst_port) in connections

        src = (src_name,src_port)
        dst = (dst_name,dst_port)

        # check if the source and destination networks exist
        if !haskey(portdict,src)
            throw(ArgumentError("Source (network name, port number) ($(src_name), $(src_port)) not found for connection ($(src_name),$(dst_name),$(src_port),$(dst_port))."))
        end
        if !haskey(portdict,dst)
            throw(ArgumentError("Destination (network name, port number) ($(dst_name), $(dst_port)) not found for connection ($(src_name),$(dst_name),$(src_port),$(dst_port))."))
        end

        src_index = portdict[src]
        dst_index = portdict[dst]

        push!(Igamma,src_index)
        push!(Jgamma,dst_index)
        push!(Vgamma,1)

        push!(Igamma,dst_index)
        push!(Jgamma,src_index)
        push!(Vgamma,1)

        # add to the vector of internal ports
        push!(portc_indices,src_index)
        push!(portc_indices,dst_index)

    end

    # then sort the internal port indices
    sort!(portc_indices)

    # check if the port indices are unique
    if !allunique(portc_indices)
        throw(ArgumentError("Duplicate connections detected [(networkname,port),counts]: $(find_duplicate_connections(connections))."))
    end

    # and generate the external port indices
    portp_indices = collect(1:m)
    deleteat!(portp_indices,portc_indices)

    # make the sparse connection matrix
    gamma = sparse(Igamma,Jgamma,Vgamma,m,m)

    # compute the block diagonal scattering matrix formed by placing all of
    # the input networks along the diagonal. Sindices contains the
    # indices from networkdata from which the components should come from.
    # S has the same sparsity pattern so we can use Sindices to copy over
    # the values from networkdata. this is useful when looping over different
    # frequencies.
    Iindices =  Vector{Int}(undef,M)
    Jindices = Vector{Int}(undef,M)
    Vindices = Vector{Tuple{Int,Int,Int}}(undef,M)
    empty!(Iindices)
    empty!(Jindices)
    empty!(Vindices)

    for i in eachindex(networkdata)
        for c in CartesianIndices(axes(networkdata[i])[1:2])
            push!(Iindices,networkdataindices[i]+c[1]-1)
            push!(Jindices,networkdataindices[i]+c[2]-1)
            push!(Vindices,(i,c[1],c[2]))
        end
    end

    # make the sparse matrix, erroring if there are multiple elements with the
    # same coordinates
    Sindices = sparse(Iindices,Jindices,Vindices,m,m,error)

    # compute an empty block diagonal scattering parameter matrix to be
    # populated later. we can reuse the sparsity structure from Sindices.
    S = SparseMatrixCSC(Sindices.m,Sindices.n,copy(Sindices.colptr),copy(Sindices.rowval),Vector{eltype(N)}(undef,M))

    return portc_indices, portp_indices, ports, networkdata, gamma, Sindices, S
end

function solveS_initialize(networks::AbstractVector,
        connections::AbstractVector; small_splitters::Bool = true,
        factorization = KLUfactorization(), internal_ports::Bool = false,
        Nmodes::Integer = 1, nbatches::Integer = Base.Threads.nthreads())

    networks_ports = add_ports(networks)

    connections_modes = add_modes(connections, Nmodes)

    networks_flat,connnections_flat = add_splitters(networks_ports,
        connections_modes; small_splitters = small_splitters)

    return solveS_initialize(networks_flat, connnections_flat;
        factorization = factorization, internal_ports = internal_ports,
        nbatches = nbatches)
end

function solveS_initialize(networks::AbstractVector{Tuple{T,N,Vector{Tuple{T, Int}}}},
        connections::AbstractVector{Tuple{T,T,Int,Int}};
        factorization = KLUfactorization(), internal_ports::Bool = false,
        nbatches::Integer = Base.Threads.nthreads()) where {T,N}

    portc_indices, portp_indices, ports, networkdata, gamma, Sindices, S = parse_connections_sparse(networks,connections)

    portsp = ports[portp_indices]
    portsc = ports[portc_indices]

    Scp_indices = Sindices[portc_indices,portp_indices]
    Spc_indices = Sindices[portp_indices,portc_indices]
    Spp_indices = Sindices[portp_indices,portp_indices]
    Scc_indices = Sindices[portc_indices,portc_indices]

    Scp = S[portc_indices,portp_indices]
    Spc = S[portp_indices,portc_indices]
    Spp = S[portp_indices,portp_indices]
    Scc = S[portc_indices,portc_indices]

    gammacc = gamma[portc_indices,portc_indices]

    sizeSp = NTuple{ndims(networkdata[1]),Int}(ifelse(i<=2,length(portsp),size(networkdata[1],i)) for i in 1:ndims(networkdata[1]))
    Sp = zeros(eltype(N),sizeSp)

    sizeSc = tuple(length(portsc),length(portsp),[size(networkdata[1],i) for i in 3:ndims(networkdata[1])]...)
    Sc = ifelse(internal_ports,zeros(eltype(N),sizeSc),zeros(eltype(N),0))

    # make gammacc - Scc. Scc can contain zeros (eg. a match at some
    # frequency). we want to retain these structural zeros because the
    # scattering matrix structure must not change.
    gammacc_Scc = spaddkeepzeros(gammacc,-Scc)

    # generate the index maps so we can perform non-allocating sparse matrix
    # addition
    gammacc_indexmap = sparseaddmap(gammacc_Scc,gammacc)
    Scc_indexmap = sparseaddmap(gammacc_Scc,Scc)

    return Sp, Sc, portsp, portsc, gammacc, Spp, Spc, Scp, Scc, Spp_indices,
        Spc_indices, Scp_indices, Scc_indices, gammacc_indexmap, Scc_indexmap,
        networkdata, nbatches, factorization, internal_ports
end


"""
    solveS_update!(Spp, Spc, Scp, Scc, Spp_indices, Spc_indices,
        Scp_indices, Scc_indices, networkdata, i)

Update the sparse matrices `Spp`, `Spc`, `Scp`, and `Scc` using the indices
from `Spp_indices`, `Spc_indices`, `Scp_indices`, and `Scc_indices` which
are indices into `networkdata` with frequency index `i`.
"""
function solveS_update!(Spp, Spc, Scp, Scc, Spp_indices, Spc_indices,
    Scp_indices, Scc_indices, networkdata, i)

    for (j,c) in enumerate(Spp_indices.nzval)
        Spp.nzval[j] = networkdata[c[1]][c[2],c[3],i]
    end

    for (j,c) in enumerate(Spc_indices.nzval)
        Spc.nzval[j] = networkdata[c[1]][c[2],c[3],i]
    end

    for (j,c) in enumerate(Scp_indices.nzval)
        Scp.nzval[j] = networkdata[c[1]][c[2],c[3],i]
    end

    for (j,c) in enumerate(Scc_indices.nzval)
        Scc.nzval[j] = networkdata[c[1]][c[2],c[3],i]
    end

    return nothing
end


function solveS_inner!(Sp,Sc,gammacc,Spp, Spc, Scp, Scc, Spp_indices, Spc_indices,
            Scp_indices, Scc_indices, gammacc_indexmap, Scc_indexmap,
            networkdata, indices, batch, factorization)

    # make a copy of the scattering matrices for each thread
    Spp = copy(Spp)
    Spc = copy(Spc)
    Scp = copy(Scp)
    Scc = copy(Scc)

    # a dense matrix version of Scp
    Scp_dense = zeros(eltype(Sp),size(Scp,1),size(Scp,2))
    ac = similar(Sp,size(Scc,1),size(Spp,1))

    # generate an empty FactorizationCache struct
    cache = FactorizationCache()

    # update the scattering matrices for the first element in the batch
    solveS_update!(Spp, Spc, Scp, Scc, Spp_indices, Spc_indices,
        Scp_indices, Scc_indices, networkdata, indices[batch[1]])

    # make gammacc - Scc. Scc can contain zeros (eg. a match at some
    # frequency). we want to retain these structural zeros because the
    # scattering matrix structure must not change.
    gammacc_Scc = spaddkeepzeros(gammacc,-Scc)

    # loop over the dimensions of the array greater than 2
    for (i,j) in enumerate(batch)
        # only perform the updates for the second or later element in the
        # batch
        if i > 1
            solveS_update!(Spp, Spc, Scp, Scc, Spp_indices, Spc_indices,
                Scp_indices, Scc_indices, networkdata, indices[j])
            fill!(gammacc_Scc,0)
            sparseadd!(gammacc_Scc,1,gammacc,gammacc_indexmap)
            sparseadd!(gammacc_Scc,-1,Scc,Scc_indexmap)
        end

        # copy only the nonzero elements. the rest of the temporary array
        # is always zero
        # Scptmp .= Scp
        for k in 1:length(Scp.colptr)-1
            for l in Scp.colptr[k]:(Scp.colptr[k+1]-1)
                Scp_dense[Scp.rowval[l],k] = Scp.nzval[l]
            end
        end

        # perform a factorization or update the factorization
        tryfactorize!(cache, factorization, gammacc_Scc)

        # solve the linear system
        # Eq. 26
        trysolve!(ac, cache.factorization, Scp_dense)

        # Eq. 28
        Sp[:,:,j] .= Spp + Spc*ac
        if !isempty(Sc)
            # derived from Eqns. 24, 28
            Sc[:,:,j] .= Scp + Scc*ac
        end
    end

    return nothing
end

"""
    solveS!(Sp, Sc, portsp, portsc, gammacc, Spp, Spc, Scp, Scc,
        Spp_indices, Spc_indices, Scp_indices, Scc_indices, networkdata,
        nbatches, factorization, internal_ports)

In-place version of `solveS`. See [`solveS`](@ref) for description. The use-
case for this function is to perform in-place updates of a network connection,
for example, by changing the arrays that are referenced in `networks` then
recomputing the scattering parameters for the connected system.

# Examples
```jldoctest
networks = [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5])];
connections = [[("S1",1),("S2",2)]];
init = JosephsonCircuits.solveS_initialize(networks, connections);
JosephsonCircuits.solveS!(init...)

# output
(S = [0.5 0.5; 0.5 0.5], ports = [("S1", 2), ("S2", 1)], Sinternal = Float64[], portsinternal = [("S1", 1), ("S2", 2)])
```

# References
V. A. Monaco and P. Tiberio, "Computer-Aided Analysis of Microwave Circuits,"
in IEEE Transactions on Microwave Theory and Techniques, vol. 22, no. 3, pp.
249-263, Mar. 1974, doi: 10.1109/TMTT.1974.1128208.
"""
function solveS!(Sp, Sc, portsp, portsc, gammacc, Spp, Spc, Scp, Scc,
    Spp_indices, Spc_indices, Scp_indices, Scc_indices, gammacc_indexmap,
    Scc_indexmap, networkdata, nbatches, factorization, internal_ports)

    # solve the linear system for the specified frequencies. the response for
    # each frequency is independent so it can be done in parallel; however
    # we want to reuse the factorization object and other input arrays. 
    # perform array allocations and factorization "nbatches" times.
    # parallelize using tasks
    indices = CartesianIndices(axes(networkdata[1])[3:end])
    batches = Base.Iterators.partition(1:length(indices),1+(length(indices)-1)÷nbatches)
    Threads.@sync for batch in batches
        Base.Threads.@spawn solveS_inner!(Sp,Sc,gammacc,Spp, Spc, Scp, Scc, Spp_indices, Spc_indices,
            Scp_indices, Scc_indices, gammacc_indexmap, Scc_indexmap,
            networkdata, indices, batch, factorization)
    end

    return (S=Sp, ports=portsp, Sinternal=Sc, portsinternal = portsc)
end

"""
    solveS(networks, connections; small_splitters::Bool = true,
        factorization = KLUfactorization(), internal_ports::Bool = false,
        nbatches::Integer = Base.Threads.nthreads())

Perform the connections between the networks in `networks` specified by the
vector of vectors of tuples `connections`. Return the sparse matrix of
scattering parameters for the external ports `S` and the external ports
`ports`. Also return the internal port scattering parameters `Sinternal` and
the internal ports `portsinternal`.

# Arguments
- `networks`: a vector of tuples of the network name, scattering parameter
    matrix, and optionally the ports  such as
    [("network1name",rand(Complex{Float64},2,2))] or
    [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5])].
- `connections::AbstractVector{<:AbstractVector{Tuple{T,Int}}}`: a vector of
    vectors of tuples of networks names and ports such as [[("S1",1),("S2",2)]]
    or [[("network1name",1),("network2name",2)]] where network1 and network2
    are the two networks being connected and 1 and 2 are integers describing
    the ports to connect.

# Keywords
- `small_splitters::Bool = true`: if true, then generate any N port splitter
    by combining (N-2) 3 port splitters. if false, then make the N port
    splitter and connect the components to it.
- `factorization = KLUfactorization()`: use KLU factorization by default. 
    JosephsonCircuits.LUfactorization() is another good choice. Keyword
    arguments can be passed to the solver as keyword arguments to these
    functions.
- `internal_ports::Bool = false`: return the scattering parameters for the
    internal ports.
- `nbatches::Integer = Base.Threads.nthreads()`: the number of batches to run
    on threads. Defaults to the number of threads with which Julia was
    launched.

# Returns
- `S`: sparse matrix of scattering parameters for the external ports.
- `ports`: the vector of tuples of network name and port number for the
    external ports.
- `Sinternal`: sparse matrix of scattering parameters for the internal ports.
- `portsinternal`: the vector of tuples of network name and port number for the
    internal ports.

# Examples
```jldoctest
networks = [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5])];
connections = [[("S1",1),("S2",2)]];
JosephsonCircuits.solveS(networks,connections;internal_ports=true)

# output
(S = [0.5 0.5; 0.5 0.5], ports = [("S1", 2), ("S2", 1)], Sinternal = [1.0 0.0; 0.5 0.5], portsinternal = [("S1", 1), ("S2", 2)])
```
```jldoctest
networks = [("S1",[0.0 1.0;1.0 0.0]),("S2",[0.5 0.5;0.5 0.5],[("S3",5),("S3",6)])];
connections = [("S1","S3",1,6)];
JosephsonCircuits.solveS(networks,connections)

# output
(S = [0.5 0.5; 0.5 0.5], ports = [("S1", 2), ("S3", 5)], Sinternal = Float64[], portsinternal = [("S1", 1), ("S3", 6)])
```

# References
V. A. Monaco and P. Tiberio, "Computer-Aided Analysis of Microwave Circuits,"
in IEEE Transactions on Microwave Theory and Techniques, vol. 22, no. 3, pp.
249-263, Mar. 1974, doi: 10.1109/TMTT.1974.1128208.
"""
function solveS(networks::AbstractVector, connections::AbstractVector;
    small_splitters::Bool = true, factorization = KLUfactorization(),
    internal_ports::Bool = false, Nmodes::Integer = 1,
    nbatches::Integer = Base.Threads.nthreads())

    init = solveS_initialize(networks,connections;
        small_splitters = small_splitters, factorization = factorization,
        internal_ports = internal_ports, Nmodes = Nmodes, nbatches = nbatches)

    return solveS!(init...)
end
