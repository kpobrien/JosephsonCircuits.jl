

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

julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.warmupsyms().signal)
JosephsonCircuits.LinearHB([0.9999636971626665 + 0.026577597828732157im 0.024776280908444095 + 0.004460970065009537im; 0.02006351652496576 - 0.01520590853015591im 0.8793137132801259 - 0.4769079135295333im;;; 0.8793041537769276 + 0.4769257467844491im 0.020066345817057213 + 0.015208700168465511im; 0.0247802510228797 - 0.004461160668153754im 0.9999643476849712 - 0.02655684759593593im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993670379479913 0.0006329620520086488; 0.0006329620520086476 0.9993670379479913;;; 0.9993668400044972 0.0006331599955028162; 0.0006331599955028159 0.9993668400044972], [0.9993670379479914 1.0; 1.0 0.9993670379479901;;; 0.9993668400044982 1.0; 1.0 0.9993668400044973], [-0.9999999999999998 -0.9999999999999991; 1.000000000000001 0.9999999999999998], ComplexF64[], ComplexF64[], ComplexF64[], 2, 3, 2, 2, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(-2,), (0,)])

julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.warmupsyms())
JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB((2.9845193040956104e10,), JosephsonCircuits.Frequencies{1}((4,), (5,), (8,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(4,)], [(1,), (3,)]), ComplexF64[-0.01317268153461008 - 0.008620192636063851im, 3.3877483616232297e-6 + 8.323112211967887e-6im, 0.12179774887406367 + 0.07965319542032037im, 2.1979490981718878e-5 + 7.557330878003312e-7im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], [1.0e-9], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), JosephsonCircuits.LinearHB([0.9999636971626665 + 0.026577597828732157im 0.024776280908444095 + 0.004460970065009537im; 0.02006351652496576 - 0.01520590853015591im 0.8793137132801259 - 0.4769079135295333im;;; 0.8793041537769276 + 0.4769257467844491im 0.020066345817057213 + 0.015208700168465511im; 0.0247802510228797 - 0.004461160668153754im 0.9999643476849712 - 0.02655684759593593im], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.9993670379479913 0.0006329620520086488; 0.0006329620520086476 0.9993670379479913;;; 0.9993668400044972 0.0006331599955028162; 0.0006331599955028159 0.9993668400044972], [0.9993670379479914 1.0; 1.0 0.9993670379479901;;; 0.9993668400044982 1.0; 1.0 0.9993668400044973], [-0.9999999999999998 -0.9999999999999991; 1.000000000000001 0.9999999999999998], ComplexF64[], ComplexF64[], ComplexF64[], 2, 3, 2, 2, 2.827433388230814e10:3.141592653589793e9:3.141592653589793e10, [(-2,), (0,)]))
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

compare(x,y) = isequal(x,y)
compare(x::AbstractArray{Complex{Float64}},y::AbstractArray{Complex{Float64}}) = isapprox(x,y,atol=1e-6)
compare(x::AbstractArray{Float64},y::AbstractArray{Float64}) = isapprox(x,y,atol=1e-6)
compare(x::StepRangeLen, y::StepRangeLen) = true

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
compare(x::JosephsonCircuits.Frequencies,y::JosephsonCircuits.Frequencies) = comparestruct(x,y)
