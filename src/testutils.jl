# Helpers used by the test suite to record solver output as Julia source
# and to compare structures field by field with a tolerance.

"""
    testshow(io::IO,S)

Print `S` to `io` in a form which can be pasted back into a test as Julia
source. The default `show` does not always produce such a form (a sparse
vector, for example), and a parameterized struct would print its full type
parameters, which are implementation details; this prints sparse vectors
as `sparsevec(...)` calls and the solver result structs as their
constructor applied to their fields.

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

Print the struct `out` to `io` as its constructor name (without type
parameters) applied to its fields, each printed with [`testshow`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.NoiseReduction([1.0, 2.0], [3.0, -4.0]))
JosephsonCircuits.NoiseReduction{Vector{Float64}}([1.0, 2.0], [3.0, -4.0])

julia> JosephsonCircuits.testshow(IOBuffer(),JosephsonCircuits.warmupsyms())
```
"""
function showstruct(io::IO,out)
  tn = typeof(out)
  fn = fieldnames(tn)
  # the constructor name without type parameters, so the printed form does
  # not depend on field types
  print(io,Base.typename(tn).wrapper,"(")
  for i in 1:length(fn)-1
    testshow(io,getfield(out,fn[i]))
    print(io,", ")
  end
  testshow(io,getfield(out,fn[end]))
  print(io,")")
end

"""
    comparestruct(x,y)

Compare two structures of the same type field by field with
[`compare`](@ref), which compares floating point arrays with a tolerance
and ignores the solver diagnostics.

# Examples
```jldoctest
julia> JosephsonCircuits.comparestruct(JosephsonCircuits.NoiseReduction([1.0], [2.0]),JosephsonCircuits.NoiseReduction([1.0], [2.0]))
true

julia> JosephsonCircuits.comparestruct(JosephsonCircuits.warmup(),JosephsonCircuits.warmup())
true

julia> JosephsonCircuits.comparestruct(nothing,nothing)
true

julia> JosephsonCircuits.compare(nothing,nothing)
true

julia> t = JosephsonCircuits.CircuitTopology(Dict((1, 2) => 1, (3, 1) => 2, (1, 3) => 2, (2, 1) => 1), JosephsonCircuits.SparseArrays.sparse([1, 2], [1, 2], [1, 1], 2, 2), 2);JosephsonCircuits.compare(t,t)
true
```
"""
function comparestruct(x,y)
  tn = typeof(x)
  fn = fieldnames(tn)
  out = true
  for i in 1:length(fn)
    fieldx = getfield(x,fn[i])
    fieldy = getfield(y,fn[i])
    out*=compare(fieldx,fieldy)
  end
  return out
end

"""
    comparearray(x::AbstractArray{T},y::AbstractArray{T}; rtol = 1e-6) where T

Whether two arrays have the same size and differ in the 2-norm by at most
`rtol` times the larger of their norms, so that the comparison holds the
same digits whatever the scale of the values, and an array compares equal
to zero only when it is zero.

# Examples
```jldoctest
julia> JosephsonCircuits.comparearray([1,2],[1,2,3])
false

julia> JosephsonCircuits.comparearray([1,2],[1,2,])
true

julia> JosephsonCircuits.comparearray([3e-8, 4e-8], [0.0, 0.0])
false

julia> JosephsonCircuits.comparearray([3e-8, 4e-8], [3e-8, 4e-8*(1 + 1e-9)])
true
```
"""
function comparearray(x::AbstractArray{T},y::AbstractArray{T};
        rtol = 1e-6) where T
    size(x) == size(y) || return false
    isempty(x) && return true
    return LinearAlgebra.norm(x[i] - y[i] for i in eachindex(x)) <=
        rtol*max(LinearAlgebra.norm(x), LinearAlgebra.norm(y))
end

"""
    compare(x,y)

Compare two values for the tests: `isequal`, except that floating point
arrays are compared with a tolerance ([`comparearray`](@ref)), the solver
result structures are compared field by field ([`comparestruct`](@ref)),
and the solver diagnostics are ignored.
"""
compare(x,y)::Bool = isequal(x,y)
compare(x::AbstractArray{Complex{Float64}},y::AbstractArray{Complex{Float64}}) = comparearray(x,y)
compare(x::AbstractArray{Float64},y::AbstractArray{Float64}) = comparearray(x,y)

compare(x::JosephsonCircuits.AbstractSparseVector,y::JosephsonCircuits.AbstractSparseVector) = compare(x.nzval,y.nzval) && compare(x.nzind,y.nzind)
compare(x::JosephsonCircuits.HB,y::JosephsonCircuits.HB) = comparestruct(x,y)
compare(x::JosephsonCircuits.NonlinearHB,y::JosephsonCircuits.NonlinearHB) = comparestruct(x,y)
# the solver diagnostics depend on the solver method and the iteration
# counts, so two solutions with different diagnostics still compare equal
compare(x::JosephsonCircuits.SolverInfo,y::JosephsonCircuits.SolverInfo) = true
compare(x::JosephsonCircuits.LinearizedHB,y::JosephsonCircuits.LinearizedHB) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitMatrices,y::JosephsonCircuits.CircuitMatrices) = comparestruct(x,y)
compare(x::JosephsonCircuits.CompiledCircuit,y::JosephsonCircuits.CompiledCircuit) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitTopology,y::JosephsonCircuits.CircuitTopology) = comparestruct(x,y)
compare(x::JosephsonCircuits.Frequencies,y::JosephsonCircuits.Frequencies) = comparestruct(x,y)
compare(x::JosephsonCircuits.PassiveNetwork,y::JosephsonCircuits.PassiveNetwork) = comparestruct(x,y)
