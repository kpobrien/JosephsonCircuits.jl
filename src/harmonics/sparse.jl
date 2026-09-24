
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

Write the elements of `A`, each repeated `Nmodes` times along the
diagonal, into `out`. Only the nonzero elements of `A` are written, so the
other entries of `out` keep what they held, zero for the result to be the
repeated matrix.

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
        throw(DimensionMismatch(lazy"Sizes not consistent"))
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
    
    # each stored entry of `A` is repeated once per mode
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
        throw(DimensionMismatch(lazy"argument shapes must match"))
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
        throw(DimensionMismatch(lazy"As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch(lazy"The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch(lazy"A and As must be the same size."))
    end

    for i in 1:nnz(As)
        A.nzval[indexmap[i]] += c*As.nzval[i]
    end
    return A
end

"""
    sparseaddconjsubst!(A::SparseMatrixCSC, c::Number, As::SparseMatrixCSC,
        indexmap, wmodes::AbstractVector, power::Integer)

Perform `A += c*As*Ad` with `Ad` the implicit diagonal whose entry in column
`i` is the signed mode frequency of that column raised to `power`,
`wmodes[(i-1) % length(wmodes) + 1]^power`, with the mode index fastest over
the nodes and any auxiliary variables. The stored value of `As` is complex
conjugated in every column whose mode frequency is negative, and a
frequency dependent entry is resolved at the mode frequency of its
column. The frequency and the conjugation of a column are computed from its
index rather than read from materialized diagonals, so the assembly loop of
[`hblinsolve`](@ref) allocates nothing at a signal frequency.
"""
function sparseaddconjsubst!(A::SparseMatrixCSC, c::Number,
    As::SparseMatrixCSC, indexmap, wmodes::AbstractVector, power::Integer)

    if nnz(A) < nnz(As)
        throw(DimensionMismatch(lazy"As cannot have more nonzero elements than A"))
    end

    if nnz(As) != length(indexmap)
        throw(DimensionMismatch(lazy"The indexmap must be the same length as As"))
    end

    if size(A) != size(As)
        throw(DimensionMismatch(lazy"A and As must be the same size."))
    end

    Nmodes = length(wmodes)
    for i in 1:length(As.colptr)-1
        wm = wmodes[(i-1) % Nmodes + 1]
        Ad = power == 0 ? one(wm) : power == 1 ? wm : wm^power
        for j in As.colptr[i]:(As.colptr[i+1]-1)
            tmp = substitutefreq(As.nzval[j], wm)
            A.nzval[indexmap[j]] += c*Ad*modevalue(tmp, wm)
        end
    end
    return A
end

"""
    sparseaddmap(A::SparseMatrixCSC, B::SparseMatrixCSC)

Return a vector of length `nnz(B)` which maps the indices of elements of `B`
in `B.nzval` to the corresponding indices in `A.nzval`. The sparse matrix `B`
must have elements in a subset of the positions in `A` which have nonzero
entries (structural zeros are elements); an element of `B` with no position
in `A` is an `ArgumentError`. Each column is one merge of the sorted row
indices of the two matrices, whatever their index type.

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
        throw(DimensionMismatch(lazy"A and B must be the same size."))
    end
    indexmap = Vector{Int}(undef, nnz(B))
    Arows = rowvals(A); Brows = rowvals(B)
    @inbounds for j in 1:size(B, 2)
        k = SparseArrays.getcolptr(A)[j]
        kmax = SparseArrays.getcolptr(A)[j+1] - 1
        for ptr in nzrange(B, j)
            r = Brows[ptr]
            while k <= kmax && Arows[k] < r
                k += 1
            end
            (k <= kmax && Arows[k] == r) || throw(ArgumentError(
                lazy"the entry ($(r), $(j)) has no position in the matrix it is mapped into."))
            indexmap[ptr] = k
        end
    end
    return indexmap
end

"""
    modevalue(v, w)

The value a linear term matrix entry contributes in a column belonging to a
mode with (signed) frequency `w`: the stored value for the non negative
frequency modes and its complex conjugate for the negative frequency modes.
This is the single definition of the negative frequency conjugation
convention, used by [`conjnegfreq!`](@ref) (which bakes it into a matrix),
by [`sparseaddconjsubst!`](@ref) (which applies it during assembly), and
by `sensitivitystampvalue` (which applies it to the sensitivity stamps).
"""
@inline modevalue(v, w) = real(w) < 0 ? conj(v) : v

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
            throw(DimensionMismatch(lazy"The dimensions of A must be integer multiples of the length of wmodes."))
        end
    end

    @inbounds for i in 1:length(A.colptr)-1
        wm = wmodes[((i-1) % length(wmodes)) + 1]
        for j in A.colptr[i]:(A.colptr[i+1]-1)
            A.nzval[j] = modevalue(A.nzval[j], wm)
        end
    end
    return A
end

"""
    freqsubst(A::SparseMatrixCSC, wmodes::Vector)

Resolve the frequency dependent elements of `A` at the vector of mode
frequencies `wmodes`, each at the magnitude of its column's mode frequency
(see [`substitutefreq`](@ref)). Returns a sparse matrix with type
`Complex{Float64}`.

# Examples
```jldoctest
wmodes = [-1,2];
f = JosephsonCircuits.FrequencyDependent;
A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [f(w->w),f(w->2*w),f(w->3*w)],2,2),2);
JosephsonCircuits.freqsubst(A,wmodes)

# output
4×4 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 6 stored entries:
 1.0+0.0im      ⋅      3.0+0.0im      ⋅
     ⋅      2.0+0.0im      ⋅      6.0+0.0im
     ⋅          ⋅      2.0+0.0im      ⋅
     ⋅          ⋅          ⋅      4.0+0.0im
```
```jldoctest
wmodes = [-1,2];
A = JosephsonCircuits.diagrepeat(JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,2,3],2,2),2);
JosephsonCircuits.freqsubst(A,wmodes)

# output
4×4 SparseArrays.SparseMatrixCSC{ComplexF64, Int64} with 6 stored entries:
 1.0+0.0im      ⋅      3.0+0.0im      ⋅    
     ⋅      1.0+0.0im      ⋅      3.0+0.0im
     ⋅          ⋅      2.0+0.0im      ⋅    
     ⋅          ⋅          ⋅      2.0+0.0im
```
"""
function freqsubst(A::SparseMatrixCSC, wmodes::Vector)

    for i in size(A)
        if i % length(wmodes) != 0
            throw(DimensionMismatch(lazy"The dimensions of A must be integer multiples of the length of wmodes."))
        end
    end

    # the output is Complex{Float64} whatever the (possibly symbolic) input
    # element type
    nzval = zeros(Complex{Float64},length(A.nzval))

    @inbounds for i in 1:length(A.colptr)-1
        for j in A.colptr[i]:(A.colptr[i+1]-1)
            if checkissymbolic(A.nzval[j])
                # `substitutefreq` evaluates the frequency dependent
                # provider leaves at the mode frequency
                substituted = substitutefreq(A.nzval[j],
                    wmodes[((i-1) % length(wmodes)) + 1])
                if checkissymbolic(substituted)
                    error(lazy"The matrix contains the symbolic value $(A.nzval[j]). If this represents a frequency dependent component, write it as a FrequencyDependent closure of the frequency. If it contains variables which should have numerical values, add them to the circuit definitions dictionary circuitdefs.")
                end
                nzval[j] = substituted
            else
                nzval[j] = A.nzval[j]
            end
        end
    end

    return SparseMatrixCSC(A.m, A.n, A.colptr, A.rowval,nzval)

end

"""
    nzposition(S::SparseMatrixCSC, i::Integer, j::Integer)

The position of the stored entry `(i, j)` in the value array of `S`, whose
row indices are sorted within a column, found by a binary search of the
column. An entry `S` does not store is an `ArgumentError`.
"""
function nzposition(S::SparseMatrixCSC, i::Integer, j::Integer)
    r = nzrange(S, j)
    k = searchsortedfirst(view(rowvals(S), r), i)
    (k <= length(r) && rowvals(S)[r[k]] == i) ||
        throw(ArgumentError(lazy"($(i), $(j)) is not a stored entry of the pattern."))
    return r[k]
end
