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
julia> JosephsonCircuits.testshow(stdout,JosephsonCircuits.warmupnumericmatrices())
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Nothing[], 2), sparsevec(Int64[], Nothing[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [50.0], [2], Int64[], 1.0e-9, [50.0, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])

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
    fieldx = getfield(x,fn[i])
    fieldy = getfield(y,fn[i])
    out*=compare(fieldx,fieldy)
  end
  return out
end

"""
    comparearray(x::AbstractArray{T},y::AbstractArray{T}) where T

Whether two arrays have the same size and differ by at most `1e-6` in the
2-norm.

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
compare(x::StepRangeLen, y::StepRangeLen) = true

compare(x::Nothing,y::Nothing) = true
compare(x::JosephsonCircuits.AbstractSparseVector,y::JosephsonCircuits.AbstractSparseVector) = compare(x.nzval,y.nzval) && compare(x.nzind,y.nzind)
compare(x::JosephsonCircuits.SparseMatrixCSC{Nothing, Int64},y::JosephsonCircuits.SparseMatrixCSC{Nothing, Int64}) = true
compare(x::JosephsonCircuits.HB,y::JosephsonCircuits.HB) = comparestruct(x,y)
compare(x::JosephsonCircuits.NonlinearHB,y::JosephsonCircuits.NonlinearHB) = comparestruct(x,y)
# the solver diagnostics depend on the solver method and the iteration
# counts, so two solutions with different diagnostics still compare equal
compare(x::JosephsonCircuits.SolverInfo,y::JosephsonCircuits.SolverInfo) = true
compare(x::JosephsonCircuits.LinearizedHB,y::JosephsonCircuits.LinearizedHB) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitMatrices,y::JosephsonCircuits.CircuitMatrices) = comparestruct(x,y)
compare(x::JosephsonCircuits.CompiledCircuit,y::JosephsonCircuits.CompiledCircuit) = comparestruct(x,y)
compare(x::JosephsonCircuits.CircuitGraph,y::JosephsonCircuits.CircuitGraph) = comparestruct(x,y)
compare(x::JosephsonCircuits.Frequencies,y::JosephsonCircuits.Frequencies) = comparestruct(x,y)
compare(x::JosephsonCircuits.PassiveNetwork,y::JosephsonCircuits.PassiveNetwork) = comparestruct(x,y)
compare(x::String,y::String) = isequal(x,y)

"""
    structurejacobian(d, Amatrixindices, Amatrixconjindices, Ljb, Lmean, Rbnm,
        Nmodes, Nbranches, Nfreq, invLnm, Gnm, Cnm, rl, cl)

The sparsity structure of the real Jacobian restricted to the mode
coupling described by `Amatrixindices` and `Amatrixconjindices`, together
with the [`StructureRealJacobianPlan`](@ref) which assembles it, taking
the linear term matrices from `d.sys` so that the assembly is the one the
solver performs. `d` is the named tuple returned by
`hbnlsolve(...; debugJacobian = true)`. Used only by the tests.
"""
function structurejacobian(d, Amatrixindices::Matrix,
    Amatrixconjindices::Matrix, Ljb, Lmean, Rbnm, Nmodes, Nbranches, Nfreq,
    invLnm, Gnm, Cnm, rl, cl)

    P, nodesandsigns = realjacobianstructure(Amatrixindices,
        Amatrixconjindices, Ljb, Rbnm, Nmodes, Nbranches, invLnm, Gnm, Cnm,
        rl, cl)
    plan = planstructurerealjacobian(P, eltype(P), Amatrixindices,
        Amatrixconjindices, Ljb, Lmean, nodesandsigns, d.sys.invLnm,
        d.sys.Gnm, d.sys.Cnm, d.sys.wmodesm, d.sys.wmodes2m, rl, cl, Nmodes,
        Nfreq, CPU(); transposed = false)
    return P, plan
end
