
"""
    testshow(io::IO,S)

Print `S` to the IOStream or IOBuffer `io`. This is used to generate the
inputs for testing purposes. The default `show` function doesn't always
produce an output which can be evaluated as the input. For example, for sparse
vectors.

# Examples
```jldoctest
julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.SparseArrays.sparsevec([1],[2],3))
sparsevec([1], [2], 3)

julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.SparseArrays.sparsevec([],Nothing[],3))
sparsevec(Int64[], Nothing[], 3)

julia> JosephsonCircuits.testshow(IOBuffer(),JosephsonCircuits.AxisKeys.KeyedArray(rand(Int8, 2,10), ([:a, :b], 10:10:100)))
```
"""
function testshow(io::IO,S::JosephsonCircuits.AbstractSparseVector)
    I = S.nzind
    V = S.nzval
    n = S.n
    print(io,"sparsevec(", I, ", ", V, ", ", n, ")")
end

testshow(io::IO,S) = show(io,S)
testshow(io::IO,S::JosephsonCircuits.HB) = showstruct(io,S)
testshow(io::IO,S::JosephsonCircuits.NonlinearHB) = showstruct(io,S)
testshow(io::IO,S::JosephsonCircuits.LinearizedHB) = showstruct(io,S)
testshow(io::IO,S::JosephsonCircuits.CircuitMatrices) = showstruct(io,S)
testshow(io::IO,S::JosephsonCircuits.AxisKeys.KeyedArray) = showstruct(io,S)

"""
    showstruct(io::IO,out)

Recursively print the struct `out` to the IOStream or IOBuffer `io`.

# Examples
```jldoctest
julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.warmupnumericmatrices())
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Nothing[], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [2], Int64[], 1.0e-9, Real[1, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])

julia> JosephsonCircuits.testshow(IOBuffer(),JosephsonCircuits.warmupsyms())
```
"""
function showstruct(io::IO,out)
  tn = typeof(out)
  fn = fieldnames(tn)
  print(io,tn,"(")
  for i in 1:length(fn)-1
    testshow(io,getfield(out,fn[i]))
    print(io,", ")
  end
  testshow(io,getfield(out,fn[end]))
  print(io,")")
end

"""
    comparestruct(x,y)

Compare two structures for testing purposes.

# Examples
```jldoctest
julia> JosephsonCircuits.comparestruct(JosephsonCircuits.warmupnumericmatrices(),JosephsonCircuits.warmupnumericmatrices())
true

julia> JosephsonCircuits.comparestruct(JosephsonCircuits.warmup(),JosephsonCircuits.warmup())
true

julia> JosephsonCircuits.comparestruct(nothing,nothing)
true

julia> JosephsonCircuits.compare(nothing,nothing)
true

julia> cg = JosephsonCircuits.CircuitGraph(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), JosephsonCircuits.SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 2), [(1, 2), (1, 3)], Tuple{Int64, Int64}[], [(1, 2), (1, 3)], Vector{Int64}[], Int64[], JosephsonCircuits.Graphs.SimpleGraphs.SimpleGraph{Int64}(2, [[2, 3], [1], [1]]), 2);JosephsonCircuits.compare(cg,cg)
true
```
"""
function comparestruct(x,y)
  tn = typeof(x)
  fn = fieldnames(tn)
  out = true
  for i in 1:length(fn)
    # println(fn[i])
    fieldx = getfield(x,fn[i])
    fieldy = getfield(y,fn[i])
    # println(compare(fieldx,fieldy))
    out*=compare(fieldx,fieldy)
  end
  return out
end

"""
    comparearray(x::AbstractArray{T},y::AbstractArray{T}) where T

Compare two arrays for testing purposes.

# Examples
```jldoctest
julia> JosephsonCircuits.comparearray([1,2],[1,2,3])
false

julia> JosephsonCircuits.comparearray([1,2],[1,2,])
true
```
"""
function comparearray(x::AbstractArray{T},y::AbstractArray{T}) where T
    if size(x) == size(y)
        z = similar(x)
        for i in eachindex(x)
            z[i] = x[i]-y[i]
        end
        return LinearAlgebra.norm(z) <= 1e-6
    else
        return false
    end
end

compare(x,y)::Bool = isequal(x,y)
compare(x::AbstractArray{Complex{Float64}},y::AbstractArray{Complex{Float64}}) = comparearray(x,y)
compare(x::AbstractArray{Float64},y::AbstractArray{Float64}) = comparearray(x,y)
compare(x::StepRangeLen, y::StepRangeLen) = true

compare(x::Nothing,y::Nothing) = true
compare(x::JosephsonCircuits.AbstractSparseVector,y::JosephsonCircuits.AbstractSparseVector) = compare(x.nzval,y.nzval) && compare(x.nzind,y.nzind)
compare(x::JosephsonCircuits.SparseMatrixCSC{Nothing, Int64},y::JosephsonCircuits.SparseMatrixCSC{Nothing, Int64}) = true
compare(x::JosephsonCircuits.HB,y::JosephsonCircuits.HB) = comparestruct(x,y)
compare(x::JosephsonCircuits.NonlinearHB,y::JosephsonCircuits.NonlinearHB) = comparestruct(x,y)
compare(x::JosephsonCircuits.LinearizedHB,y::JosephsonCircuits.LinearizedHB) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitMatrices,y::JosephsonCircuits.CircuitMatrices) = comparestruct(x,y)
compare(x::JosephsonCircuits.ParsedSortedCircuit,y::JosephsonCircuits.ParsedSortedCircuit) = comparestruct(x,y)
compare(x::JosephsonCircuits.ParsedCircuit,y::JosephsonCircuits.ParsedCircuit) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitGraph,y::JosephsonCircuits.CircuitGraph) = comparestruct(x,y)
compare(x::JosephsonCircuits.Frequencies,y::JosephsonCircuits.Frequencies) = comparestruct(x,y)
compare(x::JosephsonCircuits.NonlinearElement,y::JosephsonCircuits.NonlinearElement) = comparestruct(x,y)
function compare(x::Dict{Int,JosephsonCircuits.NonlinearElement},y::Dict{Int,JosephsonCircuits.NonlinearElement})
    if length(x) != length(y)
        return false
    end
    for (k, v) in x
        if !haskey(y, k) || !compare(v, y[k])
            return false
        end
    end
    return true
end
