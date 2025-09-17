
"""
    diagrepeat(A::AbstractVecOrMat, Nmodes::Integer)

Return a matrix with each element of `A` duplicated along the diagonal
`Nmodes` times.

# Examples
```jldoctest
julia> JosephsonCircuits.diagrepeat([1 2;3 4],2)
4×4 Matrix{Int64}:
 1  0  2  0
 0  1  0  2
 3  0  4  0
 0  3  0  4

julia> JosephsonCircuits.diagrepeat([1,2],2)
4-element Vector{Int64}:
 1
 1
 2
 2
```
"""
function diagrepeat(A::AbstractVecOrMat, Nmodes::Integer)
    out = zeros(eltype(A),size(A).*Nmodes)
    diagrepeat!(out,A,Nmodes)
    return out
end


"""
    diagrepeat(A::AbstractArray, Nmodes::Integer)

Return a array with each element of the first two axes of `A` duplicated along
the diagonal `Nmodes` times.

# Examples
```jldoctest
julia> JosephsonCircuits.diagrepeat([1 2;3 4;;;],2)
4×4×1 Array{Int64, 3}:
[:, :, 1] =
 1  0  2  0
 0  1  0  2
 3  0  4  0
 0  3  0  4
```
"""
function diagrepeat(A::AbstractArray, Nmodes::Integer)
    # only scale the first two dimensions
    sizeout = NTuple{ndims(A),Int}(ifelse(i == 1 || i == 2, Nmodes*val, val) for (i,val) in enumerate(size(A)))
    out = zeros(eltype(A),sizeout)
    return diagrepeat!(out,A,Nmodes)
end

"""
    diagrepeat!(out::AbstractVecOrMat, A::AbstractVecOrMat, Nmodes::Integer)

Overwrite `out` with the elements of `A` duplicated `Nmodes` times along
the diagonal.

# Examples
```jldoctest
julia> A = [1 2;3 4];out = zeros(eltype(A),4,4);JosephsonCircuits.diagrepeat!(out,A,2);out
4×4 Matrix{Int64}:
 1  0  2  0
 0  1  0  2
 3  0  4  0
 0  3  0  4
```
"""
function diagrepeat!(out::AbstractVecOrMat, A::AbstractVecOrMat, Nmodes::Integer)

    if size(A).*Nmodes != size(out)
        throw(DimensionMismatch("Sizes not consistent"))
    end

    @inbounds for coord in CartesianIndices(A)
        if !iszero(A[coord])
            for i in 1:Nmodes
                out[CartesianIndex((coord.I .- 1).*Nmodes .+ i)] = A[coord]
            end
        end
    end

    return out
end

function diagrepeat!(out::AbstractArray, A::AbstractArray, Nmodes::Integer)
    # use views to loop over the dimensions of the
    # array higher than 2.
    for i in CartesianIndices(axes(A)[3:end])
        diagrepeat!(view(out,:,:,i),view(A,:,:,i),Nmodes)
    end
    return out
end

"""
    diagrepeat(A::Diagonal, Nmodes::Integer)

Return a diagonal matrix with each element of `A` duplicated along the
diagonal `Nmodes` times.

# Examples
```jldoctest
julia> JosephsonCircuits.diagrepeat(JosephsonCircuits.LinearAlgebra.Diagonal([1,2]),2)
4×4 LinearAlgebra.Diagonal{Int64, Vector{Int64}}:
 1  ⋅  ⋅  ⋅
 ⋅  1  ⋅  ⋅
 ⋅  ⋅  2  ⋅
 ⋅  ⋅  ⋅  2
```
"""
function diagrepeat(A::Diagonal, Nmodes::Integer)
    out = zeros(eltype(A),length(A.diag)*Nmodes)
    diagrepeat!(out,A.diag,Nmodes)
    return Diagonal(out)
end

"""
    diagrepeat(A::SparseMatrixCSC, Nmodes::Integer)

Return a sparse matrix with each element of `A` duplicated along the diagonal 
`Nmodes` times.

# Examples
```jldoctest
julia> JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,1,2,2], [1,2,1,2], [1,2,3,4],2,2),2)
4×4 SparseArrays.SparseMatrixCSC{Int64, Int64} with 8 stored entries:
 1  ⋅  2  ⋅
 ⋅  1  ⋅  2
 3  ⋅  4  ⋅
 ⋅  3  ⋅  4
```
"""
function diagrepeat(A::SparseMatrixCSC, Nmodes::Integer)

    # column pointer has length number of columns + 1
    colptr = Vector{Int}(undef,Nmodes*A.n+1)
    
    # the sum of the number of nonzero elements is an upper bound 
    # for the number of nonzero elements in the sum.
    # set rowval and nzval to be that size then reduce size later
    rowval = Vector{Int}(undef,Nmodes*nnz(A))
    nzval = Vector{eltype(A.nzval)}(undef,Nmodes*nnz(A))

    diagrepeat!(colptr,rowval,nzval,A,Nmodes)

    return SparseMatrixCSC(A.m*Nmodes,A.n*Nmodes,colptr,rowval,nzval)
end

function diagrepeat!(colptr::Vector, rowval::Vector, nzval::Vector,
    A::SparseMatrixCSC, Nmodes::Integer)

    colptr[1] = 1
    # loop over the columns
    @inbounds for i in 2:length(A.colptr)
        # the diagonally repeated elements are additional columns
        # in between the original columns with the elements shifted
        # down.
        for k in 1:Nmodes
            idx = (i-2)*Nmodes+k+1
            colptr[idx] = colptr[(i-2)*Nmodes+k]
            for j in A.colptr[i-1]:(A.colptr[i]-1)
                rowval[colptr[idx]] = (A.rowval[j]-1)*Nmodes+k
                nzval[colptr[idx]] = A.nzval[j]
                colptr[idx] += 1
            end
        end
    end
    return nothing
end

"""
    diagrepeat(A::SparseVector, Nmodes::Integer)

Return a sparse vector with each element of `A` duplicated along the diagonal 
`Nmodes` times.

# Examples
```jldoctest
julia> JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparsevec([1,2],[1,2]),2)
4-element SparseArrays.SparseVector{Int64, Int64} with 4 stored entries:
  [1]  =  1
  [2]  =  1
  [3]  =  2
  [4]  =  2
```
"""
function diagrepeat(A::SparseVector, Nmodes::Integer)

    # define empty vectors for the rows, columns, and values
    nzind = Vector{eltype(A.nzind)}(undef,nnz(A)*Nmodes)
    nzval = Vector{eltype(A.nzval)}(undef,nnz(A)*Nmodes)

    @inbounds for i in 1:length(A.nzind)
        for j in 1:Nmodes
            nzind[(i-1)*Nmodes+j] = (A.nzind[i]-1)*Nmodes+j
            nzval[(i-1)*Nmodes+j] = A.nzval[i]
        end
    end

    return SparseVector(A.n*Nmodes,nzind,nzval)
end

"""
    diagcombine(x::Vector{T}) where T<:AbstractArray

Accept a vector of abstract arrays `x` where each element is a matrix or array
of scattering parameters for one mode. Returns a single matrix or array of the 
multi-mode scattering parameter matrices.

# Examples
```jldoctest
julia> JosephsonCircuits.diagcombine([[111 121;211 221],[112 122;212 222],[113 123;213 223]])
6×6 Matrix{Int64}:
 111    0    0  121    0    0
   0  112    0    0  122    0
   0    0  113    0    0  123
 211    0    0  221    0    0
   0  212    0    0  222    0
   0    0  213    0    0  223

julia> JosephsonCircuits.diagcombine([[111 121;211 221;;;],[112 122;212 222;;;],[113 123;213 223;;;]])
6×6×1 Array{Int64, 3}:
[:, :, 1] =
 111    0    0  121    0    0
   0  112    0    0  122    0
   0    0  113    0    0  123
 211    0    0  221    0    0
   0  212    0    0  222    0
   0    0  213    0    0  223
```
"""
function diagcombine(x::Vector{T}) where T<:AbstractArray
    
    # check if the sizes of the matrices are consistent
    sizex = size(first(x))
    for xi in x
        if size(xi) != sizex
            throw(DimensionMismatch("Sizes are not consistent."))
        end
    end

    # prepare a new tuple with the size of the scattering
    # parameter array.
    # the first two dimensions are multiplied by the length
    # of the 
    outsize = NTuple{length(sizex),Int}(ifelse(i == 1 || i == 2, length(x)*val, val) for (i,val) in enumerate(sizex))

    out = zeros(eltype(eltype(x)),outsize)

    # now assign the elements
    for (i,xi) in enumerate(x)
        diagcombine!(out,xi,i)
    end

    return out
end

function diagcombine!(out::AbstractVecOrMat,A::AbstractVecOrMat,mode_index::Int)
    # should this only operate on matrices?
    # i think this function is really simple
    # should i check that the size is sufficient?
#     Nports = size(A,1)
    Nmodes = size(out,1) ÷ size(A,1)
    
    if mode_index <= 0
        throw(ArgumentError("mode_index, $(mode_index), is less than or equal to zero."))
    end
    
    if mode_index > Nmodes
        throw(ArgumentError("mode_index, $(mode_index), cannot be larger than Nmodes, which is $(Nmodes) from the matrix sizes."))
    end
        
    for coord in CartesianIndices(A)
        out[CartesianIndex((coord.I .- 1).*Nmodes .+ mode_index)] = A[coord]
    end
    return out
end

function diagcombine!(out::AbstractArray, A::AbstractArray, mode_index::Integer)
    # use views to loop over the dimensions of the
    # array higher than 2.
    for i in CartesianIndices(axes(A)[3:end])
        diagcombine!(view(out,:,:,i),view(A,:,:,i),mode_index)
    end
    return out
end


"""
    axis_to_modes(S::AbstractArray, modes_axis::Integer)


# Examples
```jldoctest
julia> JosephsonCircuits.axis_to_modes([111 121;211 221;;; 112 122;212 222;;; 113 123;213 223],3)
6×6 Matrix{Int64}:
 111    0    0  121    0    0
   0  112    0    0  122    0
   0    0  113    0    0  123
 211    0    0  221    0    0
   0  212    0    0  222    0
   0    0  213    0    0  223

julia> JosephsonCircuits.axis_to_modes([111 121;211 221;;;; 112 122;212 222;;;; 113 123;213 223],4)
6×6×1 Array{Int64, 3}:
[:, :, 1] =
 111    0    0  121    0    0
   0  112    0    0  122    0
   0    0  113    0    0  123
 211    0    0  221    0    0
   0  212    0    0  222    0
   0    0  213    0    0  223

julia> JosephsonCircuits.axis_to_modes([111 121;211 221;;;; 112 122;212 222;;;; 113 123;213 223],3)
2×2×3 Array{Int64, 3}:
[:, :, 1] =
 111  121
 211  221

[:, :, 2] =
 112  122
 212  222

[:, :, 3] =
 113  123
 213  223
```
"""
function axis_to_modes(S::AbstractArray, modes_axis::Integer)

    if modes_axis > ndims(S)
        error("modes_axis larger than number of dimensions in input array.")
    end

    if modes_axis < 3
        error("modes_axis must be 3 or more (the first two dimensions are ports.")
    end

    if ndims(S) < 3
        error("The input array needs 3 or more dimensions (two for ports and one for modes).")
    end

    # the size of S with the dimension we will delete removed
    indicesS = NTuple{ndims(S)-1,Int}(ifelse(i < modes_axis, i, i+1) for i in 1:ndims(S)-1)

    # allocate an array with the
    # don't assume the number of input ports are equal to the number
    # of output ports
    sizeout = NTuple{ndims(S)-1,Int}(ifelse(i == 1 || i == 2, size(S,modes_axis)*size(S,i), size(S,i)) for i in indicesS)


    # make sure to allocate a matrix of zeros because some of the
    # elements will be zero
    out = zeros(eltype(S),sizeout)

    # copy over the data
    Nmodes = size(S,modes_axis)
    
    # loop over the dimensions of the array greater than 2
    for c in CartesianIndices(axes(out)[3:end])
        # 
        if modes_axis == 3
            axis_to_modes!(view(out,:,:,c),view(S,:,:,:,c.I[(modes_axis-2):end]...))
         else
            axis_to_modes!(view(out,:,:,c),view(S,:,:,c.I[1:(modes_axis-3)]...,:,c.I[(modes_axis-2):end]...))
        end
    end

    return out
end

function axis_to_modes!(X::AbstractMatrix,S::AbstractArray)

    Nmodes = size(S,3)

    # check that the dimensions are consistent
    if ndims(S) != 3
        error("S parameter array must have 3 dimensions (the first two dimensions are ports and the last is the modes).")
    end

    if size(X,1) != size(S,1)*Nmodes || size(X,2) != size(S,2)*Nmodes
        error("Dimensions of X and S not consistent.")
    end
    # copy over the data
    for k in axes(S,3)
        for j in axes(S,2)
            for i in axes(S,1)
                X[(i-1)*Nmodes+k,(j-1)*Nmodes+k] = S[i,j,k]
            end
        end
    end

    return X
end



"""
    spaddkeepzeros(A::SparseMatrixCSC, B::SparseMatrixCSC)

Add sparse matrices `A` and `B` and return the result, keeping any structural
zeros, unlike the default Julia sparse matrix addition functions. 

# Examples
```jldoctest
julia> A = JosephsonCircuits.SparseArrays.sprand(10,10,0.2); B = JosephsonCircuits.SparseArrays.sprand(10,10,0.2);JosephsonCircuits.spaddkeepzeros(A,B) == A+B
true
```
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,0],2,2);
B = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1],2,2);
JosephsonCircuits.spaddkeepzeros(A,B)

# output
2×2 SparseArrays.SparseMatrixCSC{Int64, Int64} with 3 stored entries:
 2  0
 ⋅  3
```
"""
function spaddkeepzeros(A::SparseMatrixCSC, B::SparseMatrixCSC)

    if !(A.m == B.m && A.n == B.n)
        throw(DimensionMismatch("argument shapes must match"))
    end

    # column pointer has length number of columns + 1
    colptr = Vector{Int}(undef,A.n+1)

    # the sum of the number of nonzero elements is an upper bound 
    # for the number of nonzero elements in the sum.
    # set rowval and nzval to be that size then reduce size later
    rowval = Vector{Int}(undef,nnz(A)+nnz(B))
    nzval = Vector{promote_type(eltype(A.nzval),eltype(B.nzval))}(undef,nnz(A)+nnz(B))
    fill!(nzval,0)

    colptr[1] = 1
    # loop over the columns and combine the row elements
    @inbounds for i in 2:length(A.colptr)
        j = A.colptr[i-1]
        jmax = A.colptr[i]-1
        k = B.colptr[i-1]
        kmax = B.colptr[i]-1
        colptr[i] = colptr[i-1]
        while j <= jmax  || k <= kmax
            if k > kmax
                rowval[colptr[i]] = A.rowval[j]
                nzval[colptr[i]] += A.nzval[j]
                j+=1
            elseif j > jmax
                rowval[colptr[i]] = B.rowval[k]
                nzval[colptr[i]] += B.nzval[k]
                k+=1
            elseif A.rowval[j] < B.rowval[k]
                rowval[colptr[i]] = A.rowval[j]
                nzval[colptr[i]] += A.nzval[j]
                j+=1
            elseif A.rowval[j] > B.rowval[k]
                rowval[colptr[i]] = B.rowval[k]
                nzval[colptr[i]] += B.nzval[k]
                k+=1
            else
                rowval[colptr[i]] = A.rowval[j]
                nzval[colptr[i]] += A.nzval[j] + B.nzval[k]
                j+=1
                k+=1
            end
            colptr[i] += 1
        end
    end
    resize!(rowval,colptr[end]-1)
    resize!(nzval,colptr[end]-1)
    return SparseMatrixCSC(A.m,A.n,colptr,rowval,nzval)
end

"""
    sprandsubset(A::SparseMatrixCSC, p::AbstractFloat, dropzeros = true)

Given a sparse matrix `A`, return a sparse matrix with random values in some
fraction of the non-zero elements with probability p. If `dropzeros = false`,
then the zeros will be retained as structural zeros otherwise they are dropped.

This is used for testing non-allocating sparse matrix addition.

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sprand(2,2,0.5)
B = JosephsonCircuits.sprandsubset(A, 0.1)
length(A.nzval) >= length(B.nzval)

# output
true
```
```jldoctest
A = JosephsonCircuits.SparseArrays.sprand(100,100,0.5)
B = JosephsonCircuits.sprandsubset(A, 0.1)
length(A.nzval) >= length(B.nzval)

# output
true
```
"""
function sprandsubset(A::SparseMatrixCSC, p::AbstractFloat, dropzeros = true)
    B = copy(A)
    for i in 1:nnz(A)
        if rand(1)[1] <= p
            B.nzval[i] = 0
        else
            B.nzval[i] = A.nzval[i]
        end
    end
    if dropzeros
        dropzeros!(B)
    end
    return B
end

"""
    sparseadd!(A::SparseMatrixCSC, As::SparseMatrixCSC, indexmap)

Add sparse matrices `A` and `As` and return the result in `A` without
performing any allocations. This is only possible if the positions of elements
in `As` are a subset of the positions of elements in `A`. The `indexmap` can
be generated with [`sparseaddmap`](@ref).

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
indexmap = JosephsonCircuits.sparseaddmap(A,As)
JosephsonCircuits.sparseadd!(A,As,indexmap)
A

# output
2×2 SparseArrays.SparseMatrixCSC{Int64, Int64} with 3 stored entries:
 4  1
 ⋅  2
```
"""
function sparseadd!(A::SparseMatrixCSC,As::SparseMatrixCSC,indexmap)
    if nnz(A) < nnz(As)
        throw(DimensionMismatch("As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch("The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch("A and As must be the same size."))
    end

    for i in 1:nnz(As)
        A.nzval[indexmap[i]] += As.nzval[i]
    end
    return A
end

"""
    sparseadd!(A::SparseMatrixCSC, c::Number, As::SparseMatrixCSC, indexmap)

Add sparse matrices `A` and `c*As` and return the result in `A`. The sparse
matrix `As` must have nonzero entries only in a subset of the positions in `A`
which have nonzero (structural zeros are ok) entries.

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
indexmap = JosephsonCircuits.sparseaddmap(A,As)
JosephsonCircuits.sparseadd!(A,2,As,indexmap)
A

# output
2×2 SparseArrays.SparseMatrixCSC{Int64, Int64} with 3 stored entries:
 7  5
 ⋅  2
```
"""
function sparseadd!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,indexmap::Vector)

    if nnz(A) < nnz(As)
        throw(DimensionMismatch("As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch("The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch("A and As must be the same size."))
    end

    for i in 1:nnz(As)
        A.nzval[indexmap[i]] += c*As.nzval[i]
    end
    return A
end

"""
    sparseadd!(A::SparseMatrixCSC, c::Number, As::SparseMatrixCSC,
        Ad::Diagonal, indexmap)

Add sparse matrices `A` and `c*As*Ad` and return the result in `A`. The sparse
matrix `As` must have nonzero entries only in a subset of the positions in `A`
which have nonzero (structural zeros are ok) entries.

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
indexmap = JosephsonCircuits.sparseaddmap(A,As)
JosephsonCircuits.sparseadd!(A,2,As,Ad,indexmap)
A

# output
2×2 SparseArrays.SparseMatrixCSC{Int64, Int64} with 3 stored entries:
 7  -19
 ⋅    2
```
"""
function sparseadd!(A::SparseMatrixCSC, c::Number, As::SparseMatrixCSC,
    Ad::Diagonal, indexmap::Vector)

    if nnz(A) < nnz(As)
        throw(DimensionMismatch("As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch("The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch("A and As must be the same size."))
    end

    if size(A) != size(Ad)
        throw(DimensionMismatch("A and Ad must be the same size."))
    end

    for i in 1:length(As.colptr)-1
        for j in As.colptr[i]:(As.colptr[i+1]-1)
            A.nzval[indexmap[j]] += c*Ad[i,i]*As.nzval[j]
        end
    end
    return A
end

"""
    sparseadd!(A::SparseMatrixCSC, c::Number, Ad::Diagonal,
        As::SparseMatrixCSC, indexmap::Vector)

Add sparse matrices `A` and `c*Ad*As` and return the result in `A`. The sparse
matrix `As` must have nonzero entries only in a subset of the positions in `A`
which have nonzero (structural zeros are ok) entries.

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3,4],2,2)
indexmap = JosephsonCircuits.sparseaddmap(A,As)
JosephsonCircuits.sparseadd!(A,2,Ad,As,indexmap)
A

# output
2×2 SparseArrays.SparseMatrixCSC{Int64, Int64} with 3 stored entries:
 7  5
 ⋅  2
```
"""
function sparseadd!(A::SparseMatrixCSC, c::Number, Ad::Diagonal,
    As::SparseMatrixCSC, indexmap::Vector)

    if nnz(A) < nnz(As)
        throw(DimensionMismatch("As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch("The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch("A and As must be the same size."))
    end

    if size(A) != size(Ad)
        throw(DimensionMismatch("A and Ad must be the same size."))
    end

    for i in 1:length(As.colptr)-1
        for j in As.colptr[i]:(As.colptr[i+1]-1)
            A.nzval[indexmap[j]] += c*Ad[As.rowval[j],As.rowval[j]]*As.nzval[j]
        end
    end
    return A
end


# """
#   sparseaddconj!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,
#     Ad::Diagonal,indexmap::Vector,conjflag::Diagonal)

# Perform the operation A+c*As*Ad and return the result in A. Take the complex
# conjugate of As for any column where conjflag is true. 

# The sparse matrix As must have nonzero elements only in a subset of the 
# positions in A which has nonzero lements.

# # Examples
# ```jldoctest
# A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
# Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
# As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
# wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
# indexmap = JosephsonCircuits.sparseaddmap(A,As)
# JosephsonCircuits.sparseaddconj!(A,2,As,Ad,indexmap,wmodesm .< 0)
# A

# # output
# 2×2 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 3 stored entries:
#  7.0-3.0im  -19.0-12.0im
#      ⋅        2.0+1.0im
# ```
# """
# function sparseaddconj!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,
#     Ad::Diagonal,indexmap::Vector,conjflag::Diagonal)

#     if nnz(A) < nnz(As)
#         throw(DimensionMismatch("As cannot have more nonzero elements than A"))
#     end

#     if nnz(As) != length(indexmap)
#         throw(DimensionMismatch("The indexmap must be the same length as As"))
#     end

#     if size(A) != size(As)
#         throw(DimensionMismatch("A and As must be the same size."))
#     end

#     if size(A) != size(Ad)
#         throw(DimensionMismatch("A and Ad must be the same size."))
#     end

#     if size(A) != size(conjflag)
#         throw(DimensionMismatch("A and conjflag must be the same size."))
#     end

#     for i in 1:length(As.colptr)-1
#         for j in As.colptr[i]:(As.colptr[i+1]-1)
#             if conjflag[i,i]
#                 A.nzval[indexmap[j]] += c*Ad[i,i]*conj(As.nzval[j])
#             else
#                 A.nzval[indexmap[j]] += c*Ad[i,i]*As.nzval[j]
#             end
#         end
#     end
#     return nothing
# end

"""
    sparseaddconjsubst!(A::SparseMatrixCSC, c::Number, As::SparseMatrixCSC,
        Ad::Diagonal, indexmap, conjflag::Diagonal, wmodesm::Diagonal,
        freqsubstindices::Vector, symfreqvar)

Perform the operation `A+c*As*Ad` and return the result in `A`. Take the
complex conjugate of `As` for any column where `conjflag = true`.

The sparse matrix `As` must have nonzero elements only in a subset of the 
positions in `A` which has nonzero lements.

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
indexmap = JosephsonCircuits.sparseaddmap(A,As)
freqsubstindices  = JosephsonCircuits.symbolicindices(As)
JosephsonCircuits.sparseaddconjsubst!(A,2,As,Ad,indexmap,wmodesm .< 0,wmodesm,freqsubstindices,nothing)
A

# output
2×2 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 3 stored entries:
 7.0-3.0im  -19.0-12.0im
     ⋅        2.0+1.0im
```
"""
function sparseaddconjsubst!(A::SparseMatrixCSC, c::Number,
    As::SparseMatrixCSC, Ad::Diagonal, indexmap, conjflag::Diagonal,
    wmodesm::Diagonal, freqsubstindices::Vector, symfreqvar)

    if nnz(A) < nnz(As)
        throw(DimensionMismatch("As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch("The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch("A and As must be the same size."))
    end

    if size(A) != size(Ad)
        throw(DimensionMismatch("A and Ad must be the same size."))
    end

    if size(A) != size(conjflag)
        throw(DimensionMismatch("A and conjflag must be the same size."))
    end

    if size(A) != size(wmodesm)
        throw(DimensionMismatch("A and wmodesm must be the same size."))
    end

    # if !(symfreqvar isa Symbolic || symfreqvar isa Num) && length(freqsubstindices) > 0
    #     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
    # end

    k = 1

    # println(length(indexmap))
    # println(length(A.nzval))
    # println(maximum(indexmap))
    # println(1:length(As.colptr)-1)
    # println(length(freqsubstindices))
    # println((freqsubstindices))
    # println(" ")

    for i in 1:length(As.colptr)-1
        for j in As.colptr[i]:(As.colptr[i+1]-1)

            # if length(freqsubstindices) > 0 && j == freqsubstindices[k]
            #     # tmp = substitute(As.nzval[j],Dict(symfreqvar=>wmodesm[i,i]))
            #     # tmp = valuetonumber(As.nzval[j],Dict(symfreqvar=>wmodesm[i,i]))
            #     tmp = valuetonumber(As.nzval[j],symfreqvar=>wmodesm[i,i])

            #     k+=1
            # else
            #     tmp = As.nzval[j]
            # end

            # i don't notice a difference between these two. so i think i 
            # should remove the freqsubstindices functionality. 
            tmp = valuetonumber(As.nzval[j],symfreqvar=>wmodesm[i,i])
            # tmp = As.nzval[j]


            if conjflag[i,i]
                # println(j)
                A.nzval[indexmap[j]] += c*Ad[i,i]*conj(tmp)
            else
                A.nzval[indexmap[j]] += c*Ad[i,i]*tmp
            end
        end
    end
    return A
end

"""
    sparseaddmap(A::SparseMatrixCSC, B::SparseMatrixCSC)

Return a vector of length `nnz(B)` which maps the indices of elements of `B`
in `B.nzval` to the corresponding indices in `A.nzval`. The sparse matrix `B`
must have elements in a subset of the positions in `A` which have nonzero
entries (structural zeros are elements).

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
As = JosephsonCircuits.SparseArrays.sparse([1], [2], [4],2,2)
JosephsonCircuits.sparseaddmap(A,As)

# output
1-element Vector{Int64}:
 2
```
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,-3],2,2)
As = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [4,2],2,2)
JosephsonCircuits.sparseaddmap(A,As)

# output
2-element Vector{Int64}:
 1
 3
```
"""
function sparseaddmap(A::SparseMatrixCSC, B::SparseMatrixCSC)

    if size(A) != size(B)
        throw(DimensionMismatch("A and B must be the same size."))
    end

    kstart=1
    lstart=1
    indexmap = zeros(Int,nnz(B))

    for i in 1:length(B.colptr)-1
        for j in B.colptr[i]:(B.colptr[i+1]-1)
            rowindexB = B.rowval[j]
            colindexB = i
            valB = B.nzval[j]
            indexB = j
            # println("searching for:   ", rowindexB," ",colindexB," ",indexB," ",abs.(valB))
            kstart,lstart,indexA = sparseaddmap_innerloop(A,B,rowindexB,colindexB,kstart,lstart)
            indexmap[indexB] = indexA
        end
    end
    return indexmap
end

function sparseaddmap_innerloop(A::SparseMatrixCSC, B::SparseMatrixCSC,
    rowindexB::Int, colindexB::Int, kstart::Int, lstart::Int)
    @inbounds for k in kstart:length(A.colptr) -1
        # println("k: ",k, " kstart: ",kstart, " length(A.colptr) -1): ",(length(A.colptr) -1))
        tmp = A.colptr[k]:(A.colptr[k+1]-1)
        for l in lstart:length(tmp)
            rowindexA = A.rowval[tmp[l]]
            colindexA = k
            valA = A.nzval[tmp[l]]
            indexA = tmp[l]

            # if we have reached the end of a loop, reset lstart to 1
            if l == length(tmp)
                lstart = 1
            end

            if rowindexA == rowindexB && colindexA == colindexB
                # println("match found:     ",rowindexA," ",colindexA," ",indexA," ",abs(valA))

                if l == length(tmp)
                    # if we have reached the end of a 
                    # column, then reset lstart back to
                    # 1 and increase kstart by 1
                    return k+1,1,indexA
                else
                    # if we haven't reached the end of
                    # a column, increase lstart by 1
                    return k,l+1,indexA
                end
            # else
                # println(l," ",lstart," ",length(tmp))
                # println("match not found: ",rowindexA," ",colindexA," ",indexA," ",abs(valA))
            end
        end
            
    end
    error("Coordinate not found. Are the positions of elements in As a subset of the positions of elements in A?")
end

"""
    conjnegfreq(A, wmodes)

Take the complex conjugate of any element of `A` which would be negative when
multipled from the right by a diagonal matrix consisting of `wmodes`
replicated along the diagonal.

Each axis of `A` should be an integer multiple of the length of `wmodes`.

# Examples
```jldoctest
julia> A = JosephsonCircuits.SparseArrays.sparse([1,2,1,2], [1,1,2,2], [1+1im,1+1im,1+1im,1+1im],2,2);JosephsonCircuits.conjnegfreq(A,[-1,1])
2×2 SparseArrays.SparseMatrixCSC{Complex{Int64}, Int64} with 4 stored entries:
 1-1im  1+1im
 1-1im  1+1im

julia> A = JosephsonCircuits.SparseArrays.sparse([1,2,1,2], [1,1,2,2], [1im,1im,1im,1im],2,2);all(A*JosephsonCircuits.LinearAlgebra.Diagonal([-1,1]) .== JosephsonCircuits.conjnegfreq(A,[-1,1]))
true
```
"""
function conjnegfreq(A::SparseMatrixCSC, wmodes::Vector)
    B = copy(A)
    conjnegfreq!(B,wmodes)
    return B
end

"""
    conjnegfreq!(A, wmodes)

Take the complex conjugate of any element of `A` which would be negative when
multipled from the right by a diagonal matrix consisting of `wmodes`
replicated along the diagonal. Overwrite `A` with the output.

Each axis of `A` should be an integer multiple of the length of `wmodes`.

# Examples
```jldoctest
julia> A = JosephsonCircuits.SparseArrays.sparse([1,2,1,2], [1,1,2,2], [1+1im,1+1im,1+1im,1+1im],2,2);JosephsonCircuits.conjnegfreq!(A,[-1,1]);A
2×2 SparseArrays.SparseMatrixCSC{Complex{Int64}, Int64} with 4 stored entries:
 1-1im  1+1im
 1-1im  1+1im
```
"""
function conjnegfreq!(A::SparseMatrixCSC, wmodes::Vector)

    for i in size(A)
        if i % length(wmodes) != 0
            throw(DimensionMismatch("The dimensions of A must be integer multiples of the length of wmodes."))
        end
    end

    @inbounds for i in 1:length(A.colptr)-1
        for j in A.colptr[i]:(A.colptr[i+1]-1)
          if wmodes[((i-1) % length(wmodes)) + 1] < 0
            A.nzval[j] = conj(A.nzval[j])
          end
        end
      end
    return A
end

"""
    symbolicindices(A)

Return the indices in `A.nzval` where the elements of the matrix `A` are
symbolic variables.

# Examples
```jldoctest
julia> @variables w;A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,1.0,3*w+1]);println(A.nzval);JosephsonCircuits.symbolicindices(A)
Num[w, 1 + 3w, 1.0]
2-element Vector{Int64}:
 1
 2

julia> A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,1.0,2+3im]);JosephsonCircuits.symbolicindices(A)
Int64[]
```
"""
function symbolicindices(A)
    
    indices = Vector{Int}(undef,0)
    
    # @inbounds for i = 1:length(A.colptr)-1
    #     for j in A.colptr[i]:(A.colptr[i+1]-1)
    #         if checkissymbolic(A.nzval[j])
    #             push!(indices,j)
    #         end
    #     end
    # end
    # return indices

    for (i,j) in enumerate(A)
        if checkissymbolic(j)
            push!(indices,i)
        end
    end
    return indices

end

function symbolicindices(A::SparseMatrixCSC)
    return symbolicindices(A.nzval)
end

"""
    checkissymbolic(a)

Check if `a` is a symbolic variable. Define a function to do this because a
different function call is required for `@syms` vs `@variables`.

# Examples
```jldoctest
julia> @variables w;JosephsonCircuits.checkissymbolic(w)
true

julia> JosephsonCircuits.checkissymbolic(1.0)
false
```
"""
function checkissymbolic(a)
    return a isa Symbolic
end

"""
    checkissymbolic(a::Num)

Check if `a` is a symbolic variable.

# Examples
```jldoctest
julia> @variables w;JosephsonCircuits.checkissymbolic(w)
true
```
"""
function checkissymbolic(a::Symbolics.Num)
    return !(Symbolics.value(a) isa Number)
end

"""
    freqsubst(A::SparseMatrixCSC, wmodes::Vector, symfreqvar)

Substitute the frequency dependent elements of `A` using the vector of mode
frequencies `wmodes` and the symbolic frequency variable `symfreqvar`. Returns
a sparse matrix with type `Complex{Float64}`.

# Examples
```jldoctest
@variables w
wmodes = [-1,2];
A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,2*w,3*w],2,2),2);
JosephsonCircuits.freqsubst(A,wmodes,w)

# output
4×4 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 6 stored entries:
 -1.0+0.0im      ⋅      -3.0+0.0im      ⋅    
      ⋅      2.0+0.0im       ⋅      6.0+0.0im
      ⋅          ⋅      -2.0+0.0im      ⋅    
      ⋅          ⋅           ⋅      4.0+0.0im
```
```jldoctest
wmodes = [-1,2];
A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,3],2,2),2);
JosephsonCircuits.freqsubst(A,wmodes,nothing)

# output
4×4 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 6 stored entries:
 1.0+0.0im      ⋅      3.0+0.0im      ⋅    
     ⋅      1.0+0.0im      ⋅      3.0+0.0im
     ⋅          ⋅      2.0+0.0im      ⋅    
     ⋅          ⋅          ⋅      2.0+0.0im
```
"""
function freqsubst(A::SparseMatrixCSC, wmodes::Vector, symfreqvar)

    for i in size(A)
        if i % length(wmodes) != 0
            throw(DimensionMismatch("The dimensions of A must be integer multiples of the length of wmodes."))
        end
    end

    if !isnothing(symfreqvar)
        if !checkissymbolic(symfreqvar)
            error("symfreqvar must be a symbolic variable (or nothing if no symbolic variables)")
        end
    end

    # set the output to be Complex{Float64} since the input type may be something
    # weird like a symbolic type
    nzval = zeros(Complex{Float64},length(A.nzval))

    @inbounds for i in 1:length(A.colptr)-1
        for j in A.colptr[i]:(A.colptr[i+1]-1)
            if checkissymbolic(A.nzval[j])
                if isnothing(symfreqvar)
                    error("Set symfreqvar equal to the symbolic variable representing frequency.")
                else
                    nzval[j] = valuetonumber(A.nzval[j],Dict(symfreqvar=>wmodes[((i-1) % length(wmodes)) + 1]))
                end
            else
                nzval[j] = A.nzval[j]
            end
        end
    end

    return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval,nzval)

end

"""
    spmatmul!(C::SparseMatrixCSC, A::SparseMatrixCSC, B::SparseMatrixCSC,
        xb::Vector{Bool})

Non-allocating sparse matrix multiplication of `A` and `B` when sparsity pattern
of product `C` is known. Based on spmatmul from
[SparseArrays.jl](https://github.com/JuliaSparse/SparseArrays.jl/blob/main/src/linalg.jl).

# Examples
```jldoctest
julia> a = JosephsonCircuits.sprand(100,100,0.1);b = JosephsonCircuits.sprand(100,100,0.1);c = a*b; d = copy(c);xb = fill(false, size(a,1));JosephsonCircuits.spmatmul!(c,a,b,xb);c == d
true
```
"""
function spmatmul!(C::SparseMatrixCSC, A::SparseMatrixCSC, B::SparseMatrixCSC,
    xb::Vector{Bool})

    if size(A,2) != size(B,1)
        throw(DimensionMismatch("Number of columns in A must equal number of rows in B."))
    end

    if size(C,1) != size(A,1)
        throw(DimensionMismatch("Number of rows in C must equal number of rows in A."))
    end

    if size(C,2) != size(B,2)
        throw(DimensionMismatch("Number of columns in C must equal number of columns in B."))
    end

    if length(xb) != size(A,1)
        throw(DimensionMismatch("Length of xb vector must equal number of rows in A."))
    end


    nnzC = nnz(C)
    ip = 1
    fill!(xb,false)
    for i in 1:size(B,2)
        if ip + size(A,1) - 1 > nnzC
            nnzC += max(size(A,1), nnzC>>2)
            resize!(C.rowval, nnzC)
            resize!(C.nzval, nnzC)
        end
        ip = SparseArrays.spcolmul!(C.rowval, C.nzval, xb, i, ip, A, B)
    end
    resize!(C.rowval, ip - 1)
    resize!(C.nzval, ip - 1)
    return C
end

"""
    ldiv_2x2(fact,b)

Solve the linear system A*x = b for x using left division when given `fact` 
which is the LU factorization of `A`.
"""
function ldiv_2x2(fact::Union{LU,StaticArrays.LU},b::AbstractVector)
    p1, p2 = fact.p

    # solve P*L*y = b where L = [l11 0; l21 l22] and P = [1 0; 0 1] if not
    # pivoting and [0 1;1 0] if pivoting
    if p1 == 1 && p2 == 2
        y1 = b[1]/fact.L[1,1]
        y2 = (b[2]-y1*fact.L[2,1])/fact.L[2,2]
    elseif p1 == 2 && p2 == 1
        y1 = b[2]/fact.L[1,1]
        y2 = (b[1]-y1*fact.L[2,1])/fact.L[2,2]
    else
        throw(ArgumentError("Unknown pivot."))
    end

    # solve U*x = y where U = [u11 u12; 0 u22]
    if iszero(fact.U[2,2])
        # if U[2,2] is zero, the matrix is singular and the linear system
        # has potentially no solution or no unique solution. assume the matrix
        # is rank 1 (solution not unique) then check if the solution we find
        # solves the linear system.
        x2 = zero(y2)
    else
        x2 = y2/fact.U[2,2]
    end

    x1 = (y1-x2*fact.U[1,2])/fact.U[1,1]
    x = StaticArrays.SVector{2}(x1, x2)

    if iszero(fact.U[2,2])
        # if the matrix is singular, check that we are returning a valid
        # solution to the linear system
        if p1 == 1 && p2 == 2
            P = StaticArrays.SMatrix{2,2}(1,0,0,1)
        else
            P = StaticArrays.SMatrix{2,2}(0,1,1,0)
        end
        if !(isequal(P*fact.L*fact.U*x,b) || isapprox(P*fact.L*fact.U*x,b))
            throw(ArgumentError("Failed to solve linear system."))
        end
    end

    return x
end

"""
    lu_2x2(A)

Return the LU factorization of a 2 by 2 matrix as a StaticArrays.LU struct.
Perform the LU factorization even if `A` is singular.
"""
function lu_2x2(A::AbstractArray)
    # decide whether or not to pivot
    if iszero(A[2,1]) && iszero(A[1,1])
        # if A[2,1] and A[1,1] are both zero, no point in pivoting
        # the matrix is singular, but still has an LU decomposition.
        u11 = zero(A[1,1])
        u12 = A[1,2]
        l21 = zero(A[2,1])
        u22 = A[2,2] - l21*u12
        p = StaticArrays.SVector{2}(1,2)
    elseif  pivot_rows(A[1,1],A[2,1]) 
        u11 = A[2,1]
        u12 = A[2,2]
        l21 = A[1,1]/u11
        u22 = A[1,2] - l21*u12
        p = StaticArrays.SVector{2}(2,1)
    else
        u11 = A[1,1]
        u12 = A[1,2]
        l21 = A[2,1]/u11
        u22 = A[2,2] - l21*u12
        p = StaticArrays.SVector{2}(1,2)
    end

    L = LinearAlgebra.LowerTriangular(StaticArrays.SMatrix{2,2}(one(l21),l21,zero(l21),one(l21)))
    U = LinearAlgebra.UpperTriangular(StaticArrays.SMatrix{2,2}(u11,zero(u11),u12,u22))
    
    return StaticArrays.LU(L,U,p)
end

"""
    pivot_rows(A11::Union{T,Complex{T}},
    A21::Union{T,Complex{T}}) where {T<:AbstractFloat}

Return true if pivoting during LU decomposition.

# Examples
```jldoctest
julia> JosephsonCircuits.pivot_rows(0.1+0.0im,0.9+0.1im)
true

julia> JosephsonCircuits.pivot_rows(0.9+0.1im,0.1+0.0im)
false
```
"""
function pivot_rows(A11::Union{T,Complex{T}},A21::Union{T,Complex{T}}) where {T<:AbstractFloat}
    # this help stability but doesn't work for Symbolic variables.
    if abs(A11) < abs(A21)
        return true
    else
        return false
    end
end

"""
    pivot_rows(A11,A21)

Return true if pivoting during LU decomposition.

# Examples
```jldoctest
julia> @variables A21;JosephsonCircuits.pivot_rows(0,A21)
true

julia> @variables A11 A21;JosephsonCircuits.pivot_rows(A11,A21)
false
```
"""
function pivot_rows(A11,A21)
    # this works for Symbolic variables, but isn't the best for numerical
    # stability.
    if iszero(A11)
        return true
    else
        return false
    end
end