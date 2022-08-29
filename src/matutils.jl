
"""
    diagrepeat(A::Matrix,counts::Integer)

Return a matrix with each element of A duplicated along the diagonal "counts"
times.

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
function diagrepeat(A::AbstractArray,counts::Integer)
    out = zeros(eltype(A),size(A).*counts)
    diagrepeat!(out,A,counts)
    return out
end

"""
    diagrepeat!(out,A,counts::Integer)

Overwrite "out" with the elements of A duplicated "counts" times along
the diagonal. 
"""
function diagrepeat!(out,A,counts::Integer)

    if size(A).*counts != size(out)
        error("Sizes not consistent")
    end

    @inbounds for coord in CartesianIndices(A)
        if A[coord] != 0
            for i in 1:counts
                out[CartesianIndex((coord.I .- 1).*counts .+ i)] = A[coord]
            end
        end
    end

    return nothing
end

"""
    diagrepeat(A::Diagonal,counts::Integer)

Return a diagonal matrix with each element of A duplicated along the diagonal 
counts times.

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
function diagrepeat(A::Diagonal,counts::Integer)
    out = zeros(eltype(A),length(A.diag)*counts)
    diagrepeat!(out,A.diag,counts)
    return Diagonal(out)
end

"""
    diagrepeat(A::SparseMatrixCSC,counts::Integer)

Return a sparse matrix with each element of A duplicated along the diagonal 
counts times.

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
function diagrepeat(A::SparseMatrixCSC,counts::Integer)

    # define empty vectors for the rows, columns, and values
    I = Vector{eltype(A.rowval)}(undef,nnz(A)*counts)
    J = Vector{eltype(A.rowval)}(undef,nnz(A)*counts)
    V = Vector{eltype(A.nzval)}(undef,nnz(A)*counts)

    @inbounds for i in 1:length(A.colptr)-1
        for j in A.colptr[i]:(A.colptr[i+1]-1)
            for k in 1:counts
                I[(j-1)*counts+k] = (A.rowval[j]-1)*counts+k
                J[(j-1)*counts+k] = (i-1)*counts+k
                V[(j-1)*counts+k] = A.nzval[j]
            end
        end
    end

    return sparse(I,J,V,A.m*counts,A.n*counts)
end


"""
    diagrepeat(A::SparseVector,counts::Integer)

Return a sparse vector with each element of A duplicated along the diagonal 
counts times.

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
function diagrepeat(A::SparseVector,counts::Integer)

    # define empty vectors for the rows, columns, and values
    I = Vector{eltype(A.nzind)}(undef,nnz(A)*counts)
    V = Vector{eltype(A.nzval)}(undef,nnz(A)*counts)

    @inbounds for i in 1:length(A.nzind)
        for j in 1:counts
            I[(i-1)*counts+j] = (A.nzind[i]-1)*counts+j
            V[(i-1)*counts+j] = A.nzval[i]
        end
    end

    return sparsevec(I,V,A.n*counts)
end


"""
    sprandsubset(A::SparseMatrixCSC,p::AbstractFloat)

Given a sparse matrix A, return a sparse matrix with random values in some
fraction of the non-zero elements with probability p. If dropzeros is set to 
false, then the zeros will be retained as structural zeros otherwise they
are dropped. 

This is used for testing non-allocating sparse matrix addition. 
"""
function sprandsubset(A::SparseMatrixCSC,p::AbstractFloat,dropzeros=true)
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
  sparseadd!(A::SparseMatrixCSC,As::SparseMatrixCSC,indexmap)

Add sparse matrices A and As and return the result in A without performing any
allocations. This is only possible if the positions of elements in As are a 
subset of the positions of elements in A. The indexmap can be generated with 
[`sparseaddmap`](@ref).

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
    return nothing
end

"""
  sparseadd!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,indexmap)

Add sparse matrices A and c*As and return the result in A. The sparse matrix 
As must have nonzero entries only in a subset of the positions in A which have 
nonzero (structural zeros are ok) entries.

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
    return nothing
end

"""
  sparseadd!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,Ad::Diagonal,indexmap)

Add sparse matrices A and c*As*Ad and return the result in A. The sparse matrix 
As must have nonzero entries only in a subset of the positions in A which have 
nonzero (structural zeros are ok) entries.

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
function sparseadd!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,Ad::Diagonal,indexmap::Vector)

    if nnz(A) < nnz(As)
        throw(DimensionMismatch("As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch("The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch("A and As must be the same size."))
    end


    for i in 1:length(As.colptr)-1
        for j in As.colptr[i]:(As.colptr[i+1]-1)
            A.nzval[indexmap[j]] += c*Ad[i,i]*As.nzval[j]
        end
    end
    return nothing
end

"""
  sparseadd!(A::SparseMatrixCSC,c::Number,Ad::Diagonal,As::SparseMatrixCSC,indexmap::Vector)

Add sparse matrices A and c*Ad*As and return the result in A. The sparse matrix 
As must have nonzero entries only in a subset of the positions in A which have 
nonzero (structural zeros are ok) entries.

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
function sparseadd!(A::SparseMatrixCSC,c::Number,Ad::Diagonal,As::SparseMatrixCSC,indexmap::Vector)

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
    return nothing
end


"""
  sparseaddconj!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,
    Ad::Diagonal,indexmap::Vector,conjflag::Diagonal)

Perform the operation A+c*As*Ad and return the result in A. Take the complex
conjugate of As for any column where conjflag is true. 

The sparse matrix As must have nonzero elements only in a subset of the 
positions in A which has nonzero lements.

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
indexmap = JosephsonCircuits.sparseaddmap(A,As)
JosephsonCircuits.sparseaddconj!(A,2,As,Ad,indexmap,wmodesm .< 0)
A

# output
2×2 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 3 stored entries:
 7.0-3.0im  -19.0-12.0im
     ⋅        2.0+1.0im
```
"""
function sparseaddconj!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,
    Ad::Diagonal,indexmap::Vector,conjflag::Diagonal)

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

    for i in 1:length(As.colptr)-1
        for j in As.colptr[i]:(As.colptr[i+1]-1)
            if conjflag[i,i]
                A.nzval[indexmap[j]] += c*Ad[i,i]*conj(As.nzval[j])
            else
                A.nzval[indexmap[j]] += c*Ad[i,i]*As.nzval[j]
            end
        end
    end
    return nothing
end


"""
  sparseaddconjsubst!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,
    Ad::Diagonal,indexmap,conjflag::Diagonal,wmodesm::Diagonal,
    freqsubstindices::Vector,symfreqvar)

Perform the operation A+c*As*Ad and return the result in A. Take the complex
conjugate of As for any column where conjflag is true. 

The sparse matrix As must have nonzero elements only in a subset of the 
positions in A which has nonzero lements.

# Examples
```jldoctest
A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1.0+1.0im,2.0+1.0im,-3.0+0.0im],2,2)
Ad = JosephsonCircuits.LinearAlgebra.Diagonal([1,-2])
As = JosephsonCircuits.SparseArrays.sparse([1,1], [1,2], [3.0+2.0im,4.0+3.0im],2,2)
wmodesm = JosephsonCircuits.LinearAlgebra.Diagonal([-1,1])
indexmap = JosephsonCircuits.sparseaddmap(A,As)
JosephsonCircuits.sparseaddconj!(A,2,As,Ad,indexmap,wmodesm .< 0)
A

# output
2×2 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 3 stored entries:
 7.0-3.0im  -19.0-12.0im
     ⋅        2.0+1.0im
```
"""
function sparseaddconjsubst!(A::SparseMatrixCSC,c::Number,As::SparseMatrixCSC,
    Ad::Diagonal,indexmap,conjflag::Diagonal,wmodesm::Diagonal,freqsubstindices::Vector,symfreqvar)

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
    return nothing
end

"""
  sparseaddmap(A::SparseMatrixCSC,B::SparseMatrixCSC)

Return a vector of length nnz(B) which maps the indices of elements
of B in B.nzval to the corresponding indices in A.nzval. The sparse matrix B 
must have elements in a subset of the positions in A which have nonzero entries
(structural zeros are elements).

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
function sparseaddmap(A::SparseMatrixCSC,B::SparseMatrixCSC)

    if size(A) != size(B)
        throw(DimensionMismatch("A and B must be the same size."))
    end

    kstart=1
    lstart=1
    indexmap = zeros(Int64,nnz(B))

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

function sparseaddmap_innerloop(A::SparseMatrixCSC,B::SparseMatrixCSC,
    rowindexB::Int64,colindexB::Int64,kstart::Int64,lstart::Int64)
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
            else
                # println(l," ",lstart," ",length(tmp))
                # println("match not found: ",rowindexA," ",colindexA," ",indexA," ",abs(valA))
            end
        end
            
    end
    error("Coordinate not found. Are the positions of elements in As a subset of the positions of elements in A?")
end

"""

  conjnegfreq(A,wmodes)

Take the complex conjugate of any element of A which would be negative when
multipled from the right by a diagonal matrix consisting of wmodes replicated
along the diagonal.

Each axis of A should be an integer multiple of the length of wmodes.

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
function conjnegfreq(A::SparseMatrixCSC,wmodes::Vector)
    B = copy(A)
    conjnegfreq!(B,wmodes)
    return B
end

"""
  conjnegfreq!(A,wmodes)

Take the complex conjugate of any element of A which would be negative when
multipled from the right by a diagonal matrix consisting of wmodes replicated
along the diagonal. Overwrite A with the output. 

Each axis of A should be an integer multiple of the length of wmodes.

# Examples
```jldoctest
julia> A = JosephsonCircuits.SparseArrays.sparse([1,2,1,2], [1,1,2,2], [1+1im,1+1im,1+1im,1+1im],2,2);JosephsonCircuits.conjnegfreq!(A,[-1,1]);A
2×2 SparseArrays.SparseMatrixCSC{Complex{Int64}, Int64} with 4 stored entries:
 1-1im  1+1im
 1-1im  1+1im
```
"""
function conjnegfreq!(A::SparseMatrixCSC,wmodes::Vector)

    for i in size(A)
        if i % length(wmodes) != 0
            error("The dimensions of A must be integer multiples of the length of wmodes.")
        end
    end

    @inbounds for i in 1:length(A.colptr)-1
        for j in A.colptr[i]:(A.colptr[i+1]-1)
          if wmodes[((i-1) % length(wmodes)) + 1] < 0
            A.nzval[j] = conj(A.nzval[j])
          end
        end
      end
    return nothing
end


"""
    symbolicindices(A)

Return the indices in A.nzval where the elements of the matrix A are symbolic
variables.

# Examples
```jldoctest
julia> @syms w;A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [w,1.0,3*w+1]);println(A.nzval);JosephsonCircuits.symbolicindices(A)
Any[w, 1 + 3w, 1.0]
2-element Vector{Int64}:
 1
 2

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
    
    indices = Vector{Int64}(undef,0)
    
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

Check if a is a symbolic variable. Define a function to do this because a
different function call is required for @syms vs @variables.

# Examples
```jldoctest
julia> @syms w;JosephsonCircuits.checkissymbolic(w)
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

Check if a is a symbolic variable.

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
    freqsubst(A::SparseMatrixCSC,wmodes::Vector,symfreqvar)

Substitute the frequency dependent elements of A using the vector of mode
frequencies wmodes and the symbolic frequency variable symfreqvar. Returns a
sparse matrix with type Complex{Float64}.

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
function freqsubst(A::SparseMatrixCSC,wmodes::Vector,symfreqvar)

    for i in size(A)
        if i % length(wmodes) != 0
            error("The dimensions of A must be integer multiples of the length of wmodes.")
        end
    end

    if !(symfreqvar == nothing || checkissymbolic(symfreqvar))
        error("Error: symfreqvar must be a symbolic variable (or nothing if no symbolic variables)")
    end

    # set the output to be Complex{Float64} since the input type may be something
    # weird like a symbolic type
    nzval = zeros(Complex{Float64},length(A.nzval))

    @inbounds for i in 1:length(A.colptr)-1
        for j in A.colptr[i]:(A.colptr[i+1]-1)
            if checkissymbolic(A.nzval[j] )

                if !checkissymbolic(symfreqvar)
                    error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                end
                nzval[j] = valuetonumber(A.nzval[j],Dict(symfreqvar=>wmodes[((i-1) % length(wmodes)) + 1]))
            else
                nzval[j] = A.nzval[j]
            end
        end
    end

    return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval,nzval)

end