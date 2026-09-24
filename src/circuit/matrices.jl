
"""
    CircuitMatrices(Cnm::SparseMatrixCSC, Gnm::SparseMatrixCSC, Lb::SparseVector
        Lbm::SparseVector, Ljb::SparseVector, Ljbm::SparseVector,
        Mb::SparseMatrixCSC, invLnm::SparseMatrixCSC,
        Rbnm::SparseMatrixCSC{Int, Int}, portindices::Vector{Int},
        portnumbers::Vector{Int}, portimpedances::Vector,
        portenvironmentindices::Vector{Int},
        noiseportimpedanceindices::Vector{Int}, Lmean, vvn)

The matrices of a compiled circuit at a given mode count: the capacitance,
conductance and inverse inductance matrices in the node basis, the
inductance vectors in the branch basis, the mutual inductance matrix, the
incidence matrix, the port data, and the resolved component values. Built
by [`numericmatrices`](@ref) and [`symbolicmatrices`](@ref).

# Fields
- `Cnm`: the capacitance matrix in the node basis with each
    element duplicated along the diagonal Nmodes times.
- `Gnm`: the conductance matrix in the node basis with each
    element duplicated along the diagonal Nmodes times.
- `Lb`: vector of branch linear inductances.
- `Lbm`: vector of branch linear inductances with each element
    duplicated Nmodes times.
- `Ljb`: vector of branch Josephson junction inductances.
- `Ljbm`: vector of branch Josephson junction inductances with
    each element duplicated Nmodes times.
- `Mb`: the mutual inductance matrix in the branch basis.
- `invLnm`: the inverse inductance matrix in the node basis, with each
    element duplicated along the diagonal Nmodes times. It excludes the
    mutually coupled inductor branches, which the solvers represent with
    auxiliary branch current variables (see circuit/mna.jl).
- `Rbnm::SparseMatrixCSC{Int, Int}`: incidence matrix to convert between the
    node and branch bases.
- `portindices::Vector{Int}`: vector of indices at which ports occur.
- `portnumbers::Vector{Int}`: vector of port numbers.
- `portimpedances::Vector`: the reference impedance of each port, ordered by
    port number. This is what the waves are normalized to, and it is defined
    for every port whether or not the port owns an environment.
- `portenvironmentindices::Vector{Int}`: vector of indices at which the port
    owned environments occur, ordered by port number, with zero for a port
    which owns none.
- `noiseportimpedanceindices::Vector{Int}`: the indices of the components
    which add thermal noise, for the noise calculations: the resistors
    other than a port's own termination, and the lossy capacitors and
    inductors.
- `Lmean`: the mean of the linear and Josephson inductances, zero when the
    circuit has none; the solvers replace it by the solver scale of
    [`calcsolverscale`](@ref) under the same name.
- `vvn`: the vector of component values with the definitions substituted.
"""
struct CircuitMatrices{TC,TG,TLb,TLbm,TLj,TLjm,TM,TiL,TLmean}
    Cnm::TC
    Gnm::TG
    Lb::TLb
    Lbm::TLbm
    Ljb::TLj
    Ljbm::TLjm
    Mb::TM
    invLnm::TiL
    Rbnm::SparseMatrixCSC{Int, Int}
    portindices::Vector{Int}
    portnumbers::Vector{Int}
    portimpedances::Vector
    portenvironmentindices::Vector{Int}
    noiseportimpedanceindices::Vector{Int}
    Lmean::TLmean
    vvn::Vector{Any}
end

# the flat value table is stored as `Vector{Any}` whatever it was built as,
# so that the matrices of a circuit have one type per element type of the
# assembled groups rather than one per way the table was typed
function CircuitMatrices(Cnm, Gnm, Lb, Lbm, Ljb, Ljbm, Mb, invLnm,
        Rbnm::SparseMatrixCSC{Int,Int}, portindices::Vector{Int},
        portnumbers::Vector{Int}, portimpedances::Vector,
        portenvironmentindices::Vector{Int},
        noiseportimpedanceindices::Vector{Int}, Lmean, vvn::AbstractVector)
    return CircuitMatrices(Cnm, Gnm, Lb, Lbm, Ljb, Ljbm, Mb, invLnm, Rbnm,
        portindices, portnumbers, portimpedances, portenvironmentindices,
        noiseportimpedanceindices, Lmean, Vector{Any}(vvn))
end

"""
    symbolicmatrices(circuit; Nmodes = 1)
    symbolicmatrices(psc::CompiledCircuit; Nmodes = 1)

The [`CircuitMatrices`](@ref) of a circuit with its component values left
symbolic, so that the capacitance and inverse inductance matrices can be
inspected as expressions. The mutually coupled inductor branches are
excluded from the inverse inductance matrix and represented by auxiliary
branch currents instead (see circuit/mna.jl), so no symbolic linear solve is
needed.

See also  [`CircuitMatrices`](@ref), [`numericmatrices`](@ref),
[`assemblematrices`](@ref), [`orderedports`](@ref),
[`portreferenceimpedances`](@ref), and [`noiseindices`](@ref).

# Examples
```julia
@variables Ipump Rleft Cc Lj Cj
circuit = Circuit(
    [:p1 => Port(1; Z0 = Rleft),
     :i1 => CurrentSource(Ipump),
     :cc => Capacitor(Cc),
     :jj => JosephsonJunction(Lj),
     :cj => Capacitor(Cj),
     :gnd => Ground()],
    [[(:p1, 1), (:i1, 1), (:cc, 1)],
     [(:cc, 2), (:jj, 1), (:cj, 1)],
     [(:p1, 2), (:i1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
JosephsonCircuits.testshow(stdout,symbolicmatrices(circuit))

# output
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}[Cc, -Cc, -Cc, Cc + Cj], 2, 2), sparse([1], [1], SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}[1 / Rleft], 2, 2), sparsevec(Int64[], Float64[], 2), sparsevec(Int64[], Float64[], 2), sparsevec([2], SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}[Lj], 2), sparsevec([2], SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}[Lj], 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], SymbolicUtils.BasicSymbolicImpl.var"typeof(BasicSymbolicImpl)"{SymReal}[Rleft], [2], Int64[], Lj, Any[Rleft, Rleft, Ipump, Cc, Lj, Cj])
```
"""
function symbolicmatrices(circuit::CompilableCircuit; Nmodes::Int = 1)
    return numericmatrices(circuit, Dict(), Nmodes = Nmodes)
end

"""
    numericmatrices(circuit, circuitdefs; Nmodes = 1)
    numericmatrices(psc::CompiledCircuit, circuitdefs;
        Nmodes = 1)
    numericmatrices(psc::CompiledCircuit, vvn; Nmodes = 1)

The [`CircuitMatrices`](@ref) of a circuit with its component values
resolved to numbers with `circuitdefs`, at the mode count `Nmodes`, with
every matrix entry repeated `Nmodes` times along the diagonal. The third
form takes the already resolved values `vvn`, so that a second call at a
different mode count (the signal grid of [`hblinsolve`](@ref) after the
pump grid of [`hbnlsolve`](@ref)) does not resolve them again.

See also [`CircuitMatrices`](@ref), [`numericmatrices`](@ref),
[`assemblematrices`](@ref), [`orderedports`](@ref),
[`portreferenceimpedances`](@ref), and [`noiseindices`](@ref).

# Examples
```jldoctest
circuit = Circuit(
    [:p1 => Port(1; Z0 = :Rleft),
     :i1 => CurrentSource(:Ipump),
     :cc => Capacitor(:Cc),
     :jj => JosephsonJunction(:Lj),
     :cj => Capacitor(:Cj),
     :gnd => Ground()],
    [[(:p1, 1), (:i1, 1), (:cc, 1)],
     [(:cc, 2), (:jj, 1), (:cj, 1)],
     [(:p1, 2), (:i1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
circuitdefs = Dict(:Lj => 1000.0e-12, :Cc => 100.0e-15, :Cj => 1000.0e-15, :Rleft => 50.0, :Ipump => 1.0e-8)
JosephsonCircuits.testshow(stdout,numericmatrices(circuit,circuitdefs))

# output
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Float64[], 2), sparsevec(Int64[], Float64[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [50.0], [2], Int64[], 1.0e-9, Any[50.0, 50.0, 1.0e-8, 1.0e-13, 1.0e-9, 1.0e-12])
```
```jldoctest
circuit = Circuit(
    [:p1 => Port(1; Z0 = :Rleft),
     :i1 => CurrentSource(:Ipump),
     :cc => Capacitor(:Cc),
     :jj => JosephsonJunction(:Lj),
     :cj => Capacitor(:Cj),
     :gnd => Ground()],
    [[(:p1, 1), (:i1, 1), (:cc, 1)],
     [(:cc, 2), (:jj, 1), (:cj, 1)],
     [(:p1, 2), (:i1, 2), (:jj, 2), (:cj, 2), (:gnd, 1)]])
circuitdefs = Dict(:Lj => 1000.0e-12, :Cc => 100.0e-15, :Cj => 1000.0e-15, :Rleft => 50.0, :Ipump => 1.0e-8)
psc = JosephsonCircuits.compile(circuit)
JosephsonCircuits.testshow(stdout,numericmatrices(psc, circuitdefs))

# output
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), sparsevec(Int64[], Float64[], 2), sparsevec(Int64[], Float64[], 2), sparsevec([2], [1.0e-9], 2), sparsevec([2], [1.0e-9], 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse(Int64[], Int64[], Float64[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [50.0], [2], Int64[], 1.0e-9, Any[50.0, 50.0, 1.0e-8, 1.0e-13, 1.0e-9, 1.0e-12])
```
"""
function numericmatrices(circuit::CompilableCircuit,
    circuitdefs::AbstractDict; Nmodes::Int = 1)
    return numericmatrices(compile(circuit), circuitdefs, Nmodes = Nmodes)
end

function numericmatrices(psc::CompiledCircuit,
    circuitdefs::AbstractDict; Nmodes::Int = 1)

    # convert as many values as we can to numerical values using definitions
    # from circuitdefs, in the one dictionary type the solvers take, so that
    # a dictionary of another type is accepted and compiles nothing new
    vvn = componentvaluestonumber(psc.componentvalues,
        definitiontable(circuitdefs))
    return numericmatrices(psc, vvn; Nmodes = Nmodes)
end

# the same, from already resolved component values, so a second call at a
# different mode count (the signal grid of hblinsolve after the pump grid
# of hbnlsolve) does not redo the symbolic value resolution.
function numericmatrices(psc::CompiledCircuit,
    vvn::AbstractVector; Nmodes::Int = 1)
    return assemblematrices(circuitmatrixplan(psc; Nmodes = Nmodes),
        bindvalues(psc, vvn))
end

"""
    calcjunctionrelations(componenttypes::Vector{Symbol},
        nodeindices::Matrix{Int}, junctioncprs::AbstractDict,
        edge2indexdict::Dict, Ljb::SparseVector)

The [`JunctionRelations`](@ref) of the junctions of a circuit, ordered by
the nonzero entries of the branch inductance vector `Ljb`, which is the
order every solver indexes its junction axis by. Returns `nothing` when
every junction is the sinusoidal Josephson one, which is the case the
solvers evaluate as plain `sin` and `cos`.

`junctioncprs` is keyed by the flat component index, as it is on a
[`CompiledCircuit`](@ref); a junction is placed by the branch its two
nodes make, the same way the assembly places its inductance.
"""
function calcjunctionrelations(componenttypes::Vector{Symbol},
        nodeindices::Matrix{Int}, junctioncprs::AbstractDict,
        edge2indexdict::Dict, Ljb::SparseVector)
    isempty(junctioncprs) && return nothing
    # the relation of each branch which holds a junction
    bycpr = Dict{Int,Any}()
    for (i, type) in enumerate(componenttypes)
        type === :Lj || continue
        b = edge2indexdict[(nodeindices[1,i], nodeindices[2,i])]
        bycpr[b] = get(junctioncprs, i, nothing)
    end
    return junctionrelations([get(bycpr, b, nothing) for b in Ljb.nzind])
end

combine_reciprocal_sum(x1,x2) = x1*x2/(x1+x2)
combine_error(x1,x2) = throw(ArgumentError(lazy"Components $(x1) and $(x2) cannot be combined to a single element. Please place the two components between different nodes."))

"""
    branchendpoints(Rbn::SparseMatrixCSC, Nbranches::Int)

The two nodes of every branch in the order the incidence matrix `Rbn`
oriented it, as `(from, to)`: branch `b` leaves `from[b]` and enters
`to[b]`. `Rbn` carries `-1` at a branch's source and `1` at its
destination, and has no column for the ground node, which is node 1, so a
branch with a terminal there keeps the initial value at that end.

One walk over the stored entries names the endpoints of every branch at
once, which is where the mutual couplings read theirs from.
"""
function branchendpoints(Rbn::SparseMatrixCSC, Nbranches::Int)
    from = fill(1, Nbranches)
    to = fill(1, Nbranches)
    rows = rowvals(Rbn)
    vals = nonzeros(Rbn)
    for j in axes(Rbn, 2)
        for p in nzrange(Rbn, j)
            if vals[p] < 0
                from[rows[p]] = j + 1
            else
                to[rows[p]] = j + 1
            end
        end
    end
    return from, to
end
