

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

julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.LinearHB([0.7298749531683114 + 0.6835806848774727im;;;], Array{ComplexF64, 3}(undef, 0, 0, 0), [1.0;;;], [1.0;;;], [0.9999999999999993;;], ComplexF64[], ComplexF64[], 1, 3, 2, 1, 3.141592653589793e10))
JosephsonCircuits.LinearHB([0.7298749531683114 + 0.6835806848774727im;;;], Array{ComplexF64, 3}(undef, 0, 0, 0), [1.0;;;], [1.0;;;], [0.9999999999999993;;], ComplexF64[], ComplexF64[], 1, 3, 2, 1, 3.141592653589793e10)

julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB(ComplexF64[-0.01317268153461008 - 0.008620192636063851im, 3.3877483616232297e-6 + 8.323112211967887e-6im, 0.12179774887406367 + 0.07965319542032037im, 2.1979490981718878e-5 + 7.557330878003312e-7im], JosephsonCircuits.SparseArrays.sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), JosephsonCircuits.SparseArrays.sparsevec([2], [1.0e-9], 2), JosephsonCircuits.SparseArrays.sparsevec(Int64[], Nothing[], 2), JosephsonCircuits.SparseArrays.sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), ComplexF64[0.9789325751107439 + 0.0im; 0.12051187816111208 + 0.07881274759778019im; -0.004230418875131352 - 0.00964520752500062im; 0.00019110015334431832 + 0.0im;;], JosephsonCircuits.LinearHB([-4.634426524178756 - 0.9355487549085718im -4.228439402152284 - 1.8637223245828725im; -4.235145340771 + 1.8484330183613003im -4.637776367875941 + 0.9187986952504916im;;;], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.511439992761294 0.4885600072387061; 0.48856000723870646 0.5114399927612936;;;], [0.5114399927612937 0.511988590306596; 0.511988590306596 0.5114399927612937;;;], [-1.0000000000000178; 0.9999999999999822;;], ComplexF64[], ComplexF64[], 2, 3, 2, 2, [2.984513020910303e10])))
JosephsonCircuits.HB(JosephsonCircuits.NonlinearHB(ComplexF64[-0.01317268153461008 - 0.008620192636063851im, 3.3877483616232297e-6 + 8.323112211967887e-6im, 0.12179774887406367 + 0.07965319542032037im, 2.1979490981718878e-5 + 7.557330878003312e-7im], sparse([1, 2, 3, 4], [1, 2, 3, 4], [1, 1, 1, 1], 4, 4), sparsevec([2], [1.0e-9], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([3, 4], [1.0e-9, 1.0e-9], 4), 2, 2, ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im]), ComplexF64[0.9789325751107439 + 0.0im; 0.12051187816111208 + 0.07881274759778019im; -0.004230418875131352 - 0.00964520752500062im; 0.00019110015334431832 + 0.0im;;], JosephsonCircuits.LinearHB([-4.634426524178756 - 0.9355487549085718im -4.228439402152284 - 1.8637223245828725im; -4.235145340771 + 1.8484330183613003im -4.637776367875941 + 0.9187986952504916im;;;], Array{ComplexF64, 3}(undef, 0, 0, 0), [0.511439992761294 0.4885600072387061; 0.48856000723870646 0.5114399927612936;;;], [0.5114399927612937 0.511988590306596; 0.511988590306596 0.5114399927612937;;;], [-1.0000000000000178; 0.9999999999999822;;], ComplexF64[], ComplexF64[], 2, 3, 2, 2, [2.984513020910303e10]))
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

julia> JosephsonCircuits.compares(nothing,nothing)
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
