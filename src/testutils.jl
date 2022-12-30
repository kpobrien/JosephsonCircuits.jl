

"""
  testshow(io::IO,S::JosephsonCircuits.AbstractSparseVector)

# Examples
```jldoctest
julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.SparseArrays.sparsevec([1],[2],3))
sparsevec([1], [2], 3)

julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.SparseArrays.sparsevec([],Nothing[],3))
sparsevec(Int64[], Nothing[], 3)
```
"""
function testshow(io::IO,S::JosephsonCircuits.AbstractSparseVector)
    I = S.nzind
    V = S.nzval
    n = S.n
    print(io,"sparsevec(", I, ", ", V, ", ", n, ")")
end


"""
  showstruct(io::IO,out)

# Examples
```jldoctest
julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.warmupnumericmatrices())
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Nothing[], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [3], Int64[], 1.0e-9, Real[1, 1.0e-8, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])
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

testshow(io::IO,S) = show(io,S)
testshow(io::IO,S::JosephsonCircuits.HB) = showstruct(io,S)
testshow(io::IO,S::JosephsonCircuits.NonlinearHB) = showstruct(io,S)
testshow(io::IO,S::JosephsonCircuits.LinearHB) = showstruct(io,S)
testshow(io::IO,S::JosephsonCircuits.CircuitMatrices) = showstruct(io,S)
# show(io::IO,S::StepRangeLen) = 

"""
  comparestruct(x,y)

# Examples
```jldoctest
julia> JosephsonCircuits.comparestruct(JosephsonCircuits.warmupnumericmatrices(),JosephsonCircuits.warmupnumericmatrices())
true

julia> JosephsonCircuits.comparestruct(JosephsonCircuits.warmup(),JosephsonCircuits.warmup())
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

compare(x,y) = isequal(x,y)
compare(x::Vector{Complex{Float64}},y::Vector{Complex{Float64}}) = isapprox(x,y)
compare(x::Vector{Float64},y::Vector{Float64}) = isapprox(x,y)

compare(x::Nothing,y::Nothing) = true
compare(x::JosephsonCircuits.AbstractSparseVector,y::JosephsonCircuits.AbstractSparseVector) = compare(x.nzval,y.nzval) && compare(x.nzind,y.nzind)
compare(x::JosephsonCircuits.SparseMatrixCSC{Nothing, Int64},y::JosephsonCircuits.SparseMatrixCSC{Nothing, Int64}) = true
compare(x::JosephsonCircuits.HB,y::JosephsonCircuits.HB) = comparestruct(x,y)
compare(x::JosephsonCircuits.NonlinearHB,y::JosephsonCircuits.NonlinearHB) = comparestruct(x,y)
compare(x::JosephsonCircuits.LinearHB,y::JosephsonCircuits.LinearHB) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitMatrices,y::JosephsonCircuits.CircuitMatrices) = comparestruct(x,y)
compare(x::JosephsonCircuits.ParsedSortedCircuit,y::JosephsonCircuits.ParsedSortedCircuit) = comparestruct(x,y)
compare(x::JosephsonCircuits.ParsedCircuit,y::JosephsonCircuits.ParsedCircuit) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitGraph,y::JosephsonCircuits.CircuitGraph) = comparestruct(x,y)
