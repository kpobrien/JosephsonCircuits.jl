
"""
    CircuitMatrices(Cnm::SparseMatrixCSC, Gnm::SparseMatrixCSC, Lb::SparseVector
        Lbm::SparseVector, Ljb::SparseVector, Ljbm::SparseVector,
        Mb::SparseMatrixCSC, invLnm::SparseMatrixCSC,
        Rbnm::SparseMatrixCSC{Int, Int}, portindices::Vector{Int},
        portnumbers::Vector{Int}, portimpedanceindices::Vector{Int}
        noiseportimpedanceindices::Vector{Int}, Lmean, vvn)

A simple structure to hold the circuit matrices including the capacitance
matrix, the conductance matrix, the inductance vectors, the Josephson
inductance vectors, the mutual inductance matrix, the inverse inductance
matrix, the incidence matrix, the dictionary of port and resistor values where
the nodes are the keys and the values are the values, and the mean of the
inductances. See also [`numericmatrices`](@ref) and [`symbolicmatrices`](@ref).

# Fields
- `Cnm::SparseMatrixCSC`: the capacitance matrix in the node basis with each
    element duplicated along the diagonal Nmodes times.
- `Gnm::SparseMatrixCSC`: the conductance matrix in the node basis with each
    element duplicated along the diagonal Nmodes times.
- `Lb::SparseVector`: vector of branch linear inductances.
- `Lbm::SparseVector`: vector of branch linear inductances with each element
    duplicated Nmodes times.
- `Ljb::SparseVector`: vector of branch Josephson junction inductances.
- `Ljbm::SparseVector`: vector of branch Josephson junction inductances with
    each element duplicated Nmodes times.
- `Mb::SparseMatrixCSC`: the mutual inductance matrix in the branch basis with
    each element duplicated along the diagonal Nmodes times.
- `invLnm::SparseMatrixCSC`: the inverse inductance matrix in the node basis
    with each element duplicated along the diagonal Nmodes times.
- `Rbnm::SparseMatrixCSC{Int, Int}`: incidence matrix to convert between the
    node and branch bases.
- `portindices::Vector{Int}`: vector of indices at which ports occur.
- `portnumbers::Vector{Int}`: vector of port numbers.
- `portimpedanceindices::Vector{Int}`: vector of indices at which port
    impedances occur.
- `noiseportimpedanceindices::Vector{Int}`: vector of indices at which
    resistive elements other than port impedances occur, for noise
    calculations.
- `Lmean`: the mean of all of the geometric and Josephson inductances.
- `vvn`: the vector of component values with numbers substituted in.

"""
struct CircuitMatrices
    Cnm::SparseMatrixCSC
    Gnm::SparseMatrixCSC
    Lb::SparseVector
    Lbm::SparseVector
    Ljb::SparseVector
    Ljbm::SparseVector
    Mb::SparseMatrixCSC
    invLnm::SparseMatrixCSC
    Rbnm::SparseMatrixCSC{Int, Int}
    portindices::Vector{Int}
    portnumbers::Vector{Int}
    portimpedanceindices::Vector{Int}
    noiseportimpedanceindices::Vector{Int}
    Lmean
    vvn
end

"""
    symbolicmatrices(circuit; Nmodes = 1, sorting = :number)

Return the symbolic matrices describing the circuit properties.

See also  [`CircuitMatrices`](@ref), [`numericmatrices`](@ref), [`calcCn`](@ref),
[`calcGn`](@ref), [`calcLb`](@ref),[`calcLjb`](@ref), [`calcMb`](@ref),
[`calcinvLn`](@ref), [`calcLmean`](@ref), [`calcportindicesnumbers`](@ref),
[`calcportimpedanceindices`](@ref), and [`calcnoiseportimpedanceindices`](@ref).

# Examples
```jldoctest
@variables Ipump Rleft Cc Lj Cj
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))
println(symbolicmatrices(circuit))

# output
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], SymbolicUtils.BasicSymbolic{Real}[Cc, -Cc, -Cc, Cc + Cj], 2, 2), sparse([1], [1], SymbolicUtils.BasicSymbolic{Real}[1 / Rleft], 2, 2), 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries, 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries,   [2]  =  Lj,   [2]  =  Lj, sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [3], Int64[], Lj, Any[1, Ipump, Rleft, Cc, Lj, Cj])
```
```jldoctest
@variables Ipump Rleft Cc Lj Cj
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))
psc = JosephsonCircuits.parsesortcircuit(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
println(symbolicmatrices(psc, cg))

# output
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], SymbolicUtils.BasicSymbolic{Real}[Cc, -Cc, -Cc, Cc + Cj], 2, 2), sparse([1], [1], SymbolicUtils.BasicSymbolic{Real}[1 / Rleft], 2, 2), 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries, 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries,   [2]  =  Lj,   [2]  =  Lj, sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [3], Int64[], Lj, Any[1, Ipump, Rleft, Cc, Lj, Cj])
```
"""
function symbolicmatrices(circuit; Nmodes = 1, sorting = :number)
    return numericmatrices(circuit, Dict(), Nmodes = Nmodes, sorting = sorting)
end

function symbolicmatrices(psc::ParsedSortedCircuit, cg::CircuitGraph; Nmodes = 1)
    return numericmatrices(psc, cg, Dict(), Nmodes = Nmodes)
end

"""
    numericmatrices(circuit, circuitdefs; Nmodes = 1, sorting = :number)

Return the numeric matrices describing the circuit properties.

See also [`CircuitMatrices`](@ref), [`numericmatrices`](@ref),
[`calcCn`](@ref), [`calcGn`](@ref), [`calcLb`](@ref),[`calcLjb`](@ref),
[`calcMb`](@ref), [`calcinvLn`](@ref), [`calcLmean`](@ref),
[`calcportindicesnumbers`](@ref), [`calcportimpedanceindices`](@ref), and
[`calcnoiseportimpedanceindices`](@ref).

# Examples
```jldoctest
@variables Ipump Rleft Cc Lj Cj
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))
circuitdefs = Dict(Lj =>1000.0e-12,Cc => 100.0e-15,Cj => 1000.0e-15,Rleft => 50.0,Ipump => 1.0e-8)
println(numericmatrices(circuit,circuitdefs))

# output
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries, 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries,   [2]  =  1.0e-9,   [2]  =  1.0e-9, sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [3], Int64[], 1.0e-9, Real[1, 1.0e-8, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])
```
```jldoctest
@variables Ipump Rleft Cc Lj Cj
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("C1","1","2",Cc)) 
push!(circuit,("Lj1","2","0",Lj)) 
push!(circuit,("C2","2","0",Cj))
circuitdefs = Dict(Lj =>1000.0e-12,Cc => 100.0e-15,Cj => 1000.0e-15,Rleft => 50.0,Ipump => 1.0e-8)
psc = JosephsonCircuits.parsesortcircuit(circuit)
cg = JosephsonCircuits.calccircuitgraph(psc)
println(numericmatrices(psc, cg, circuitdefs))

# output
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], [1.0e-13, -1.0e-13, -1.0e-13, 1.1e-12], 2, 2), sparse([1], [1], [0.02], 2, 2), 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries, 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries,   [2]  =  1.0e-9,   [2]  =  1.0e-9, sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [3], Int64[], 1.0e-9, Real[1, 1.0e-8, 50.0, 1.0e-13, 1.0e-9, 1.0e-12])
```
"""
function numericmatrices(circuit, circuitdefs; Nmodes = 1, sorting = :number)

    # parse the circuit
    psc = parsesortcircuit(circuit, sorting = sorting)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)

    return numericmatrices(psc, cg, circuitdefs, Nmodes = Nmodes)
end

function numericmatrices(psc::ParsedSortedCircuit, cg::CircuitGraph,
    circuitdefs; Nmodes = 1)

    # convert as many values as we can to numerical values using definitions
    # from circuitdefs
    vvn = componentvaluestonumber(psc.componentvalues, circuitdefs)
    
    # capacitance matrix
    Cnm = calcCn(psc.componenttypes, psc.nodeindices, vvn, Nmodes, psc.Nnodes)

    # conductance matrix
    Gnm = calcGn(psc.componenttypes, psc.nodeindices, vvn, Nmodes, psc.Nnodes)

    # branch inductance vector with Nmodes = 1
    Lb = calcLb(psc.componenttypes, psc.nodeindices, vvn, cg.edge2indexdict,
        1, cg.Nbranches)

    # branch inductance vector with Nmodes = Nmodes
    Lbm = calcLb(psc.componenttypes, psc.nodeindices, vvn, cg.edge2indexdict,
        Nmodes, cg.Nbranches)

    # branch Josephson inductance vector with Nmodes = 1
    Ljb = calcLjb(psc.componenttypes, psc.nodeindices, vvn, cg.edge2indexdict,
        1, cg.Nbranches)

    # branch Josephson inductance vector with Nmodes = Nmodes
    Ljbm = calcLjb(psc.componenttypes, psc.nodeindices, vvn, cg.edge2indexdict,
        Nmodes, cg.Nbranches)

    # mutual branch inductance matrix
    Mb = calcMb(psc.componenttypes, psc.nodeindices, vvn, psc.componentnamedict,
        psc.mutualinductorbranchnames, cg.edge2indexdict, 1, cg.Nbranches)

    # inverse nodal inductance matrix from branch inductance vector and branch
    # inductance matrix
    invLnm = calcinvLn(Lb, Mb, cg.Rbn, Nmodes)

    # expand the size of the incidence matrix
    Rbnm = diagrepeat(cg.Rbn, Nmodes)

    # calculate Lmean
    Lmean = calcLmean(psc.componenttypes, vvn)

    portindices, portnumbers = calcportindicesnumbers(psc.componenttypes,
        psc.nodeindices, psc.mutualinductorbranchnames, vvn)

    portimpedanceindices = calcportimpedanceindices(psc.componenttypes,
        psc.nodeindices, psc.mutualinductorbranchnames, vvn)

    noiseportimpedanceindices = calcnoiseportimpedanceindices(psc.componenttypes,
        psc.nodeindices, psc.mutualinductorbranchnames, vvn)

    return CircuitMatrices(Cnm, Gnm, Lb, Lbm, Ljb, Ljbm, Mb, invLnm, Rbnm,
        portindices, portnumbers, portimpedanceindices,
        noiseportimpedanceindices, Lmean, vvn)
end

"""
    calcIb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
        componentvalues::Vector, edge2indexdict::Dict, Nmodes, Nbranches)

Calculate the sparse branch current source vector whose length is
`Nbranches*Nmodes`. Note that `nodeindices` is "one indexed" so 1 is the
ground node.

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
componenttypes = [:I,:C,:L,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [1e-9, 0.2, 4e-9, 1e-12]
componentnamedict = Dict{Symbol, Int}(:C2 => 4,:L1 => 3,:I1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ib = JosephsonCircuits.calcIb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 1 stored entry:
  [1]  =  1.0e-9
```
```jldoctest
@variables I1 C1 L1 C2
Nmodes = 1
Nbranches = 2
componenttypes = [:I,:C,:L,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [I1, C1, L1, C2]
componentnamedict = Dict{Symbol, Int}(:C2 => 4,:L1 => 3,:I1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ib = JosephsonCircuits.calcIb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 1 stored entry:
  [1]  =  I1
```
"""
function calcIb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    componentvalues::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(componenttypes, nodeindices, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:I]), edge2indexdict,
        Nmodes, Nbranches, :I)
end

"""
    calcVb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
        componentvalues::Vector, edge2indexdict::Dict, Nmodes, Nbranches)

Calculate the sparse branch voltage source vector whose length is
`Nbranches*Nmodes`. Note that `nodeindices` is "one indexed" so 1 is the
ground node.

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
componenttypes = [:V,:C,:L1,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [1e-9, 0.2, 4e-9, 1e-12]
componentnamedict = Dict{Symbol, Int}(:C2 => 4,:L1 => 3,:V1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Vb = JosephsonCircuits.calcVb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 1 stored entry:
  [1]  =  1.0e-9
```
```jldoctest
@variables V1 C1 L1 C2
Nmodes = 1
Nbranches = 2
componenttypes = [:V,:C,:L,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [V1, C1, L1, C2]
componentnamedict = Dict{Symbol, Int}(:C2 => 4,:L1 => 3,:V1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Vb = JosephsonCircuits.calcVb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 1 stored entry:
  [1]  =  V1
```
"""
function calcVb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    componentvalues::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(componenttypes, nodeindices, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:V]), edge2indexdict,
        Nmodes, Nbranches, :V)
end

"""
    calcLb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
        componentvalues::Vector, edge2indexdict::Dict, Nmodes, Nbranches)

Calculate the sparse branch inductance vector whose length is
`Nbranches*Nmodes`. Note that `nodeindices` is "one indexed" so 1 is the
ground node.

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
componenttypes = [:L,:K,:L,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [1e-9, 0.2, 4e-9, 1e-12]
componentnamedict = Dict{Symbol, Int}(:C2 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Lb = JosephsonCircuits.calcLb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 2 stored entries:
  [1]  =  1.0e-9
  [2]  =  4.0e-9
```
```jldoctest
@variables L1 K1 L2 C1
Nmodes = 1
Nbranches = 2
componenttypes = [:L,:K,:L,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [L1, K1, L2, C1]
componentnamedict = Dict{Symbol, Int}(:C1 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Lb = JosephsonCircuits.calcLb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 2 stored entries:
  [1]  =  L1
  [2]  =  L2
```
"""
function calcLb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    componentvalues::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(componenttypes, nodeindices, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:L,:K]),
        edge2indexdict, Nmodes, Nbranches, :L)
end

"""
    calcLjb(componenttypes, nodeindices, componentvalues, edge2indexdict,
        Nmodes, Nbranches)

Calculate the sparse branch Josephson inductance vector whose length is
`Nbranches*Nmodes`. Note that `nodeindices` is "one indexed" so 1 is the
ground node.

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
componenttypes = [:Lj,:C,:Lj,:C]
nodeindices = [2 3 3 3; 1 2 1 1]
componentvalues = [1e-9, 1e-12, 4e-9, 1e-12]
componentnamedict = Dict{Symbol, Int}(:C2 => 4,:L2 => 3,:L1 => 1,:Cc => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ljb = JosephsonCircuits.calcLjb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 2 stored entries:
  [1]  =  1.0e-9
  [2]  =  4.0e-9
```
```jldoctest
@variables Lj1 K1 Lj2 C1
Nmodes = 1
Nbranches = 2
componenttypes = [:Lj,:K,:Lj,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [Lj1, K1, Lj2, C1]
componentnamedict = Dict{Symbol, Int}(:C1 => 4,:Lj2 => 3,:Lj1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ljb = JosephsonCircuits.calcLjb(componenttypes,nodeindices,componentvalues,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 2 stored entries:
  [1]  =  Lj1
  [2]  =  Lj2
```
"""
function calcLjb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    componentvalues::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(componenttypes, nodeindices, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:Lj]), edge2indexdict,
        Nmodes, Nbranches, :Lj)
end

"""
    calcbranchvector(componenttypes::Vector{Symbol},
        nodeindices::Matrix{Int}, componentvalues::Vector,
        valuecomponenttypes::Vector, edge2indexdict::Dict, Nmodes, Nbranches,
        component::Symbol)

Calculate the sparse branch vector whose length is `Nbranches*Nmodes` for the
given component symbol. Note that `nodeindices` is "one indexed" so 1 is
the ground node.
"""
function calcbranchvector(componenttypes::Vector{Symbol},
    nodeindices::Matrix{Int}, componentvalues::Vector,
    valuecomponenttypes::Vector, edge2indexdict::Dict, Nmodes, Nbranches,
    component::Symbol)

    # calculate the expected number of elements
    Nelements = 0
    for (i,type) in enumerate(componenttypes)
        if type == component
            Nelements += 1
        end
    end

    # define empty vectors for the indices and values
    Ib = Vector{Int}(undef, Nelements)
    Vb = Vector{eltype(valuecomponenttypes)}(undef, Nelements)

    # copy the components over
    j = 1
    for (i,type) in enumerate(componenttypes)
        if type == component
            Ib[j] = edge2indexdict[(nodeindices[1,i],nodeindices[2,i])]
            Vb[j] = convert(eltype(valuecomponenttypes), componentvalues[i])
            j += 1
        end
    end

    # return a sparse vector
    branchvector = sparsevec(Ib,Vb,Nbranches)
    if Nmodes == 1
        return branchvector
    else
        return diagrepeat(branchvector, Nmodes)
    end
end

"""
    calcMb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
        componentvalues::Vector, componentnamedict::Dict,
        mutualinductorbranchnames::Vector, edge2indexdict::Dict, Nmodes,
        Nbranches)

Returns the branch mutual inductance matrix. Note that `nodeindices` is
"one indexed" so 1 is the ground node.

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
componenttypes = [:L,:K,:L,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [1e-9, 0.2, 2e-9, 1e-12]
componentnamedict = Dict{Symbol, Int}(:C2 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
mutualinductorbranchnames = [ :L1, :L2]
Mb = JosephsonCircuits.calcMb(componenttypes,nodeindices,componentvalues,componentnamedict,mutualinductorbranchnames,edge2indexdict,Nmodes,Nbranches)

# output
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
  ⋅           2.82843e-10
 2.82843e-10   ⋅ 
```
```jldoctest
@variables L1 L2 K1 C1
Nmodes = 2
Nbranches = 2
componenttypes = [:L,:K,:L,:C]
nodeindices = [2 0 3 3; 1 0 1 1]
componentvalues = [L1, K1, L2, C1]
componentnamedict = Dict{Symbol, Int}(:C1 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
mutualinductorbranchnames = [ :L1, :L2]
Mb = JosephsonCircuits.calcMb(componenttypes,nodeindices,componentvalues,componentnamedict,mutualinductorbranchnames,edge2indexdict,Nmodes,Nbranches)

# output
4×4 SparseArrays.SparseMatrixCSC{Num, Int64} with 4 stored entries:
              ⋅               ⋅  K1*sqrt(L1*L2)               ⋅
              ⋅               ⋅               ⋅  K1*sqrt(L1*L2)
 K1*sqrt(L1*L2)               ⋅               ⋅               ⋅
              ⋅  K1*sqrt(L1*L2)               ⋅               ⋅
```
"""
function calcMb(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    componentvalues::Vector, componentnamedict::Dict,
    mutualinductorbranchnames::Vector, edge2indexdict::Dict, Nmodes,
    Nbranches)
    return calcMb_inner(componenttypes, nodeindices, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:L,:K]),
        componentnamedict, mutualinductorbranchnames, edge2indexdict, Nmodes,
        Nbranches)
end

function calcMb_inner(componenttypes::Vector{Symbol},
    nodeindices::Matrix{Int}, componentvalues::Vector,
    valuecomponenttypes::Vector, componentnamedict::Dict,
    mutualinductorbranchnames::Vector, edge2indexdict::Dict, Nmodes,
    Nbranches)

    # define empty vectors of zero length for the row indices, column indices,
    # and values
    Ib = Vector{Int}(undef, 0)
    Jb = Vector{Int}(undef, 0)
    Vb = Vector{eltype(valuecomponenttypes)}(undef, 0)

    n = 1
    #loop through componenttypes for mutual inductors
    @inbounds for (i,type) in enumerate(componenttypes)
        # when we find a mutual inductor:
        # -find the value of the mutual inductor in componentvalues[i]
        # -find the names of the two inductors it couples together from
        #   mutualinductorbranchnames[n]
        #  -look up the index of the inductors in
        #     index=componentnamedict[inductorsymbol] for each inductor symbol
        #  -given the index of the inductor, look of the value of the inductor
        #     from componentvalues
        #  -then compute the value of the mutual inductance from the two
        #     inductor values and K

        # then use the index of the inductors to get the nodes from
        # nodeindices. use that as a key in edge2indexdict to look up the
        # branch index then assign those to I, J, V for the sparse array.
        # then do the usual step of expanding that to Nmodes after finishing
        # this loop.

        if type == :K
            # value of K
            K = componentvalues[i]
            # names of inductors
            inductor1name = mutualinductorbranchnames[2*n-1]
            inductor2name = mutualinductorbranchnames[2*n]

            # indices of inductors
            inductor1index = componentnamedict[inductor1name]
            inductor2index = componentnamedict[inductor2name]

            # values of inductors
            inductor1value = componentvalues[inductor1index]
            inductor2value = componentvalues[inductor2index]

            # values of mutual inductance Lm
            Lm = K*sqrt(inductor1value*inductor2value)

            inductor1edge = (nodeindices[1,inductor1index],nodeindices[2,inductor1index])
            inductor2edge = (nodeindices[1,inductor2index],nodeindices[2,inductor2index])

            # add the edges
            push!(Ib,edge2indexdict[inductor1edge])
            push!(Jb,edge2indexdict[inductor2edge])
            pushval!(Vb,Lm,1,false)

            push!(Ib,edge2indexdict[inductor2edge])
            push!(Jb,edge2indexdict[inductor1edge])
            pushval!(Vb,Lm,1,false)

            n+=1
        end

    end

    # if there is only one frequency mode return the sparse matrix generated
    # from these vectors. if multiple modes then duplicate along the diagonal
    # Nmodes times.
    if Nmodes == 1
        return sparse(Ib,Jb,Vb,Nbranches,Nbranches)
    else
        return diagrepeat(sparse(Ib,Jb,Vb,Nbranches,Nbranches),Nmodes)
    end
end

"""
    calcinvLn(Lb::SparseVector, Rbn::SparseMatrixCSC, Nmodes)

Returns the nodal inverse inductance matrix. Accepts the vector of branch
inductances `Lb` and the incidence matrix `Rbn`.

# Examples
```jldoctest
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[1e-9,4e-9])
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
JosephsonCircuits.calcinvLn(Lb,Rbn,Nmodes)

# output
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 1.0e9   ⋅ 
  ⋅     2.5e8
```
```jldoctest
@variables L1 L2
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[L1,L2])
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
JosephsonCircuits.calcinvLn(Lb,Rbn,Nmodes)

# output
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 2 stored entries:
 1 / L1       ⋅
      ⋅  1 / L2
```
```jldoctest
@variables L1 L2
Nmodes = 2
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[L1,L2])
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
JosephsonCircuits.calcinvLn(Lb,Rbn,Nmodes)

# output
4×4 SparseArrays.SparseMatrixCSC{Num, Int64} with 4 stored entries:
 1 / L1       ⋅       ⋅       ⋅
      ⋅  1 / L1       ⋅       ⋅
      ⋅       ⋅  1 / L2       ⋅
      ⋅       ⋅       ⋅  1 / L2
```
```jldoctest
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([],Nothing[])
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
JosephsonCircuits.calcinvLn(Lb,Rbn,Nmodes).nzval

# output
Nothing[]
```
```jldoctest
@syms L1 L2
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[L1,L2])
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
JosephsonCircuits.calcinvLn(Lb,Rbn,Nmodes).nzval

# output
2-element Vector{Any}:
 1 / L1
 1 / L2
```
"""
function calcinvLn(Lb::SparseVector, Rbn::SparseMatrixCSC, Nmodes)
    if nnz(Lb)>0
        s = transpose(Rbn[Lb.nzind,:])*spdiagm(0 => 1 ./Lb.nzval)*Rbn[Lb.nzind,:]
        if Nmodes == 1
            return s
        else
            return diagrepeat(s,Nmodes)
        end 
    else
        return spzeros(eltype(Lb),Nmodes*size(Rbn)[2],Nmodes*size(Rbn)[2])
    end
end

"""
    calcinvLn(Lb::SparseVector, Mb::SparseMatrixCSC,
        Rbn::SparseMatrixCSC, Nmodes)

Returns the nodal inverse inductance matrix. Accepts the vector of branch
inductances `Lb`, the branch mutual inductance matrix `Mb`, and the incidence
matrix `Rbn`.

Using ldiv instead of an inverse: (where the extra div is an escape sequence)
Can solve A x = B with: x = A \\ B or x = invA * B, so we can perform the
inverse here with:
s = RbnT * invL * Rbn or s = RbnT * (L \\ Rbn), the latter of which should be
faster and more numerically stable.

# Examples
```jldoctest
Nmodes = 2
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[1e-9,4e-9])
Mb = JosephsonCircuits.SparseArrays.sparse([2,1], [1,2], [4e-10,4e-10])
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
JosephsonCircuits.calcinvLn(Lb,Mb,Rbn,Nmodes)

# output
4×4 SparseArrays.SparseMatrixCSC{Float64, Int64} with 8 stored entries:
  1.04167e9    ⋅         -1.04167e8    ⋅ 
   ⋅          1.04167e9    ⋅         -1.04167e8
 -1.04167e8    ⋅          2.60417e8    ⋅ 
   ⋅         -1.04167e8    ⋅          2.60417e8
```
```jldoctest
@variables L1 L2 Lm
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[L1,L2]);
Mb = JosephsonCircuits.SparseArrays.sparse([2,1], [1,2], [Lm,Lm]);
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1.0,1.0])
println(JosephsonCircuits.calcinvLn(Lb,Mb,Rbn,Nmodes))

# output
sparse([1, 2, 1, 2], [1, 1, 2, 2], Num[(1.0 + (Lm*(Lm / L1)) / (L2 + (-(Lm^2)) / L1)) / L1, (-(Lm / L1)) / (L2 + (-(Lm^2)) / L1), (-(Lm / (L2 + (-(Lm^2)) / L1))) / L1, 1.0 / (L2 + (-(Lm^2)) / L1)], 2, 2)
```
```jldoctest
@syms L1 L2 Lm
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[L1,L2]);
Mb = JosephsonCircuits.SparseArrays.sparse([2,1], [1,2], [Lm,Lm]);
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1.0,1.0])
println(JosephsonCircuits.calcinvLn(Lb,Mb,Rbn,Nmodes))

# output
sparse([1, 2, 1, 2], [1, 1, 2, 2], Num[(1.0 + (Lm*(Lm / L1)) / (L2 + (-(Lm^2)) / L1)) / L1, (-(Lm / L1)) / (L2 + (-(Lm^2)) / L1), (-(Lm / (L2 + (-(Lm^2)) / L1))) / L1, 1.0 / (L2 + (-(Lm^2)) / L1)], 2, 2)
```
```jldoctest
@variables L1 L2
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[L1,L2]);
Mb = JosephsonCircuits.SparseArrays.sparse([], [], Nothing[]);
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
println(JosephsonCircuits.calcinvLn(Lb,Mb,Rbn,Nmodes))

# output
sparse([1, 2], [1, 2], Num[1 / L1, 1 / L2], 2, 2)
```
"""
function calcinvLn(Lb::SparseVector, Mb::SparseMatrixCSC,
    Rbn::SparseMatrixCSC, Nmodes)
    if nnz(Lb) > 0 &&  nnz(Mb) > 0
        valuetype = promote_type(eltype(Lb),eltype(Mb))
    else 
        valuetype = eltype(Lb)
    end

    return calcinvLn_inner(Lb,Mb,Rbn,Nmodes,Array{valuetype,1}(undef,0))
end

function calcinvLn_inner(Lb::SparseVector, Mb::SparseMatrixCSC,
    Rbn::SparseMatrixCSC, Nmodes, valuecomponenttypes::Vector)

    # if there are no mutual inductors, return the inverse of the diagonal
    # elements
    if nnz(Lb) > 0 &&  nnz(Mb) == 0
        return calcinvLn(Lb,Rbn,Nmodes)
    elseif nnz(Lb) > 0 &&  nnz(Mb) > 0
        # add the mutual inductance matrix to a diagonal matrix made from the
        # inductance vector.
        # we pick out only the indices where there are inductors for
        # efficiency reasons.
        # calculate the symbolic inductance matrix
        if eltype(valuecomponenttypes) <: Symbolic

            # take a subset of the arrays
            Mbs = Mb[Lb.nzind,Lb.nzind]
            Lbs = Lb[Lb.nzind]

            # add together the two sparse arrays

            # define empty vectors for the rows, columns, and values
            In = Vector{Int}(undef, nnz(Mbs)+nnz(Lbs))
            Jn = Vector{Int}(undef, nnz(Mbs)+nnz(Lbs))
            Vn = Vector{eltype(valuecomponenttypes)}(undef, nnz(Mbs)+nnz(Lbs))

            for i = 1:length(Mbs.colptr)-1
                for j in Mbs.colptr[i]:(Mbs.colptr[i+1]-1)
                    In[j] = Mbs.rowval[j]
                    Jn[j] = i
                    Vn[j] = Mbs.nzval[j]
                end
            end

            for i = 1:length(Lbs.nzind)
                In[i+nnz(Mbs)] = Lbs.nzind[i]
                Jn[i+nnz(Mbs)] = Lbs.nzind[i]
                Vn[i+nnz(Mbs)] = Lbs.nzval[i]
            end

            Lsparse = sparse(In,Jn,Vn,Mbs.m,Mbs.n)


            # L = Mb[Lb.nzind,Lb.nzind] + Diagonal(Lb[Lb.nzind])
            # L = Mb[Lb.nzind,Lb.nzind] + spdiagm(0=>Lb[Lb.nzind])

            # println(eltype(valuecomponenttypes)
            # println("type: ",eltype(valuecomponenttypes) <: Symbolic)


            L = Array{Any,2}(undef,size(Lsparse))
            fill!(L,0)

            for i = 1:length(Lsparse.colptr)-1
                for j in Lsparse.colptr[i]:(Lsparse.colptr[i+1]-1)
                    L[Lsparse.rowval[j],i] = Lsparse.nzval[j]
                end
            end

        else
            # calculate the numeric inductance matrix
            L = Mb[Lb.nzind,Lb.nzind] + Diagonal(Lb[Lb.nzind])
        end

        # calculate the inverse inductance matrix
        if eltype(valuecomponenttypes) <: Union{AbstractFloat, Complex}
            # using ldiv with klu. fastest option in most cases
            s = transpose(Rbn[Lb.nzind,:])*sparse(KLU.klu(L) \ Matrix(Rbn[Lb.nzind,:]))
        else
            s = calcsymbolicinvLn(L,Lb,Rbn)
        end 

        if Nmodes == 1
            return s
        else
            return diagrepeat(s,Nmodes)
        end 
    else
        return spzeros(eltype(Lb),Nmodes*size(Rbn)[2],Nmodes*size(Rbn)[2])
    end
end

# Symbolics is needed for calculating the inverse inductance matrix if there
# are mutual inductors. For large systems with many mutual inductors this may
# result in very large equations so is not recommended. 
function calcsymbolicinvLn(L,Lb,Rbn)
    s =  sparse(transpose(Rbn[Lb.nzind,:])*(Symbolics.sym_lu(Matrix(L))\ Matrix(Rbn[Lb.nzind,:])))
    return SparseMatrixCSC(s.m, s.n, s.colptr, s.rowval,s.nzval)
    # return SparseMatrixCSC(s.m, s.n, s.colptr, s.rowval,Symbolics.value.(s.nzval))
end

"""
    calcLmean(componenttypes::Vector{Symbol}, componentvalues::Vector)

Return the mean of the linear and Josephson inductors.

# Examples
```jldoctest
julia> JosephsonCircuits.calcLmean([:R,:L,:C,:Lj],[10,4,5,1])
2.5

julia> @variables R1 L1 C1 Lj1;JosephsonCircuits.calcLmean([:R,:L,:C,:Lj],[R1, L1, C1, Lj1])
(1//2)*(L1 + Lj1)
```
"""
function calcLmean(componenttypes::Vector{Symbol}, componentvalues::Vector)
    return calcLmean_inner(componenttypes, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:Lj, :L]))
end

"""
    calcLmean_inner(componenttypes::Vector, componentvalues::Vector,
        valuecomponenttypes::Vector{Nothing})

Return the mean of the linear and Josephson inductors. Return 0 if the expected
return type is Nothing.

# Examples
```jldoctest
julia> JosephsonCircuits.calcLmean_inner([:R,:C,:C,:P],[10,4,5,1],Nothing[])
0

julia> @variables R1 L1 C1 Lj1;JosephsonCircuits.calcLmean_inner([:R,:C,:C,:C],[R1, L1, C1, Lj1],Nothing[])
0
```
"""
function calcLmean_inner(componenttypes::Vector, componentvalues::Vector,
    valuecomponenttypes::Vector{Nothing})
    return 0
end

"""
    calcLmean_inner(componenttypes::Vector, componentvalues::Vector,
        valuecomponenttypes::Vector)

Return the mean of the linear and Josephson inductors.

# Examples
```jldoctest
julia> JosephsonCircuits.calcLmean_inner([:R,:L,:C,:Lj],[10,4,5,1],Float64[])
2.5

julia> JosephsonCircuits.calcLmean_inner([:R,:C,:C,:C],[10,4,5,1],Float64[])
0.0

julia> @variables R1 L1 C1 Lj1;JosephsonCircuits.calcLmean_inner([:R,:L,:C,:Lj],[R1, L1, C1, Lj1], Num[])
(1//2)*(L1 + Lj1)
```
"""
function calcLmean_inner(componenttypes::Vector, componentvalues::Vector,
    valuecomponenttypes::Vector)

    if length(componenttypes) != length(componentvalues)
        throw(DimensionMismatch("componenttypes and componentvalues should have the same length"))
    end

    # count the number of inductors
    ninductors = 0
    for (i,type) in enumerate(componenttypes)
        if type == :L || type == :Lj
            ninductors += 1
        end
    end

    # it's a litle absurd but we have to copy the inductance values into
    # a new array to perform a type stable mean over some elements of
    # componentvalues
    Vn = Array{eltype(valuecomponenttypes), 1}(undef, ninductors)
    j = 1
    for (i,type) in enumerate(componenttypes)
        if type == :L || type == :Lj
            Vn[j] = convert(eltype(valuecomponenttypes),componentvalues[i])
            j += 1
        end
    end

    # take the mean. mean will return NaN if there are no elements so return
    # zero manually if Vn is empty.
    if isempty(Vn)
        return zero(eltype(valuecomponenttypes))
    else
        return Statistics.mean(Vn)
    end
end

"""
    calcCn(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
        componentvalues::Vector, Nmodes, Nnodes)

Returns the node capacitance matrix from the capacitance values in
`componentvalues` when `componenttypes` has the symbol `:C` with node indices
from `nodeindices`. Other symbols are ignored. Capacitances to ground
become diagonal elements. Capacitance between elements is an off-diagonal
element with a minus sign and is added to the diagonal with a plus sign. The
dimensions of the output are `(Nnodes-1)*Nmodes` by `(Nnodes-1)` times
`Nmodes` where `Nnodes` is the number of nodes including ground and `Nmodes`
is the number of different frequencies. Note that `nodeindices` is
"one indexed" so 1 is the ground node.

# Examples
```jldoctest
julia> JosephsonCircuits.calcCn([:C,:C],[2 3;1 1],[1.0,2.0],1,3)
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 1.0   ⋅ 
  ⋅   2.0

julia> JosephsonCircuits.calcCn([:C,:C,:C],[2 2 3;1 3 1],[1.0,0.1,2.0],1,3)
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 4 stored entries:
  1.1  -0.1
 -0.1   2.1

julia> JosephsonCircuits.calcCn([:C,:C,:C],[2 2 3;1 3 1],[1.0,0.1,2.0],2,3)
4×4 SparseArrays.SparseMatrixCSC{Float64, Int64} with 8 stored entries:
  1.1    ⋅   -0.1    ⋅ 
   ⋅    1.1    ⋅   -0.1
 -0.1    ⋅    2.1    ⋅ 
   ⋅   -0.1    ⋅    2.1

julia> @variables Cg1 Cg2;JosephsonCircuits.calcCn([:C,:C],[2 3;1 1],[Cg1,Cg2],1,3)
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 2 stored entries:
 Cg1    ⋅
   ⋅  Cg2

julia> @variables Cg1 Cc Cg2;JosephsonCircuits.calcCn([:C,:C,:C],[2 2 3;1 3 1],[Cg1, Cc, Cg1],1,3)
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 4 stored entries:
 Cc + Cg1       -Cc
      -Cc  Cc + Cg1

julia> @variables Cg1 Cc Cg2;JosephsonCircuits.calcCn([:C,:C,:C],[2 2 3;1 3 1],[Cg1, Cc, Cg1],2,3)
4×4 SparseArrays.SparseMatrixCSC{Num, Int64} with 8 stored entries:
 Cc + Cg1         ⋅       -Cc         ⋅
        ⋅  Cc + Cg1         ⋅       -Cc
      -Cc         ⋅  Cc + Cg1         ⋅
        ⋅       -Cc         ⋅  Cc + Cg1
```
"""
function calcCn(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    componentvalues::Vector, Nmodes, Nnodes)
    return calcnodematrix(componenttypes, nodeindices, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:C]), Nmodes, Nnodes, :C, false)
end

"""
    calcGn(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
        componentvalues::Vector, Nmodes, Nnodes)

Returns the node conductance matrix from the resistance values in
`componentvalues` when `componenttypes` has the symbol `:R`. The node indices
are taken from `nodeindices`. Conductances to ground are diagonal elements.
Conductance between elements is an off-diagonal element with a minus sign and
is added to the diagonal with a plus sign. The dimensions of the output are
`(Nnodes-1)` times `Nmodes` by `(Nnodes-1)` times `Nmodes`. Note that
`nodeindices` is "one indexed" so 1 is the ground node.

We have to calculate the inverse of the individual components so select a type
that allows that.

# Examples
```jldoctest
julia> JosephsonCircuits.calcGn([:R,:R],[2 3;1 1],[1.0,2.0],1,3)
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
 1.0   ⋅ 
  ⋅   0.5

julia> JosephsonCircuits.calcGn([:R,:R,:R],[2 2 3;1 3 1],[1.0,100.0,2.0],1,3)
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 4 stored entries:
  1.01  -0.01
 -0.01   0.51

julia> JosephsonCircuits.calcGn([:R,:R,:R],[1 3 1;2 2 3],[1.0,100.0,2.0],1,3)
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 4 stored entries:
  1.01  -0.01
 -0.01   0.51

julia> JosephsonCircuits.calcGn([:R,:R,:R],[2 2 3;1 3 1],[1.0,100.0,2.0],2,3)
4×4 SparseArrays.SparseMatrixCSC{Float64, Int64} with 8 stored entries:
  1.01    ⋅    -0.01    ⋅ 
   ⋅     1.01    ⋅    -0.01
 -0.01    ⋅     0.51    ⋅ 
   ⋅    -0.01    ⋅     0.51

julia> @variables Rg1 Rg2;JosephsonCircuits.calcGn([:R,:R],[2 3;1 1],[Rg1,Rg2],1,3)
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 2 stored entries:
 1 / Rg1        ⋅
       ⋅  1 / Rg2

julia> @variables Rg1 Rc Rg2;JosephsonCircuits.calcGn([:R,:R,:R],[2 2 3;1 3 1],[Rg1,Rc,Rg2],1,3)
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 4 stored entries:
 1 / Rc + 1 / Rg1           -1 / Rc
          -1 / Rc  1 / Rc + 1 / Rg2

julia> @variables Rg1 Rc Rg2;JosephsonCircuits.calcGn([:R,:R,:R],[2 2 3;1 3 1],[Rg1,Rc,Rg2],2,3)
4×4 SparseArrays.SparseMatrixCSC{Num, Int64} with 8 stored entries:
 1 / Rc + 1 / Rg1                 ⋅           -1 / Rc                 ⋅
                ⋅  1 / Rc + 1 / Rg1                 ⋅           -1 / Rc
          -1 / Rc                 ⋅  1 / Rc + 1 / Rg2                 ⋅
                ⋅           -1 / Rc                 ⋅  1 / Rc + 1 / Rg2
```
"""
function calcGn(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
    componentvalues::Vector, Nmodes, Nnodes)

    return calcnodematrix(componenttypes, nodeindices, componentvalues,
        calcvaluetype(componenttypes, componentvalues, [:R]), Nmodes, Nnodes,
        :R, true)
end

"""
    calcnodematrix(componenttypes::Vector{Symbol}, nodeindices::Matrix{Int},
        componentvalues::Vector, valuecomponenttypes::Vector, Nmodes, Nnodes,
        component::Symbol, invert::Bool)

Returns either the capacitance or conductance matrix depending on the values
of `component` and `invert`. `:C` and `false` for capacitance and `:R` and
`true` for conductance. The dimensions of the output are `(Nnodes-1)` times
`Nmodes` by `(Nnodes-1)` times `Nmodes`. Note that `nodeindices` is
"one indexed" so 1 is the ground node.
"""
function calcnodematrix(componenttypes::Vector{Symbol},
    nodeindices::Matrix{Int}, componentvalues::Vector,
    valuecomponenttypes::Vector, Nmodes, Nnodes, component::Symbol,
    invert::Bool)

    if length(componenttypes) != size(nodeindices,2) || length(componenttypes) != length(componentvalues)
        throw(DimensionMismatch("componenttypes, nodeindices, and componentvalues should have the same length"))
    end

    if size(nodeindices,1) != 2
        throw(DimensionMismatch("nodeindices should have a first dimension size of 2."))
    end

    # calculate the expected number of elements
    Nelements = 0
    for (i,type) in enumerate(componenttypes)
        if type == component
            if nodeindices[1,i] == 1
                Nelements += 1
            elseif nodeindices[2,i] == 1
                Nelements += 1
            else
                Nelements += 4
            end
        end
    end

    # define empty vectors for the row indices, column indices, and values
    In = Vector{Int}(undef, Nelements)
    Jn = Vector{Int}(undef, Nelements)
    Vn = Vector{eltype(valuecomponenttypes)}(undef, Nelements)

    j=1
    # generate the capacitance or conductance matrix values for Nmodes=1
    for (i,type) in enumerate(componenttypes)
        if type == component

            if nodeindices[1,i] == 1
                # capacitance to ground, add to diagonal
                In[j] = nodeindices[2,i]-1
                Jn[j] = nodeindices[2,i]-1
                # convert is necessary here to avoid allocations
                Vn[j] = convert(eltype(valuecomponenttypes),componentvalues[i])
                if invert
                    Vn[j] = 1/Vn[j]
                end
                j+=1

            elseif nodeindices[2,i] == 1
                # capacitance to ground, add to diagonal
                In[j] = nodeindices[1,i]-1
                Jn[j] = nodeindices[1,i]-1
                # convert is necessary here to avoid allocations
                Vn[j] = convert(eltype(valuecomponenttypes),componentvalues[i])
                if invert
                    Vn[j] = 1/Vn[j]
                end
                j+=1

            else
                # diagonal elements
                In[j] = nodeindices[1,i]-1
                Jn[j] = nodeindices[1,i]-1
                # convert is necessary here to avoid allocations
                Vn[j] = convert(eltype(valuecomponenttypes),componentvalues[i])
                if invert
                    Vn[j] = 1/Vn[j]
                end
                j+=1

                In[j] = nodeindices[2,i]-1
                Jn[j] = nodeindices[2,i]-1
                # convert is necessary here to avoid allocations
                Vn[j] = convert(eltype(valuecomponenttypes),componentvalues[i])
                if invert
                    Vn[j] = 1/Vn[j]
                end
                j+=1

                # off diagonal elements
                In[j] = nodeindices[1,i]-1
                Jn[j] = nodeindices[2,i]-1
                # convert is necessary here to avoid allocations
                Vn[j] = convert(eltype(valuecomponenttypes),componentvalues[i])
                Vn[j] = -Vn[j]
                if invert
                    Vn[j] = 1/Vn[j]
                end
                j+=1

                In[j] = nodeindices[2,i]-1
                Jn[j] = nodeindices[1,i]-1
                # convert is necessary here to avoid allocations
                Vn[j] = convert(eltype(valuecomponenttypes),componentvalues[i])
                Vn[j] = -Vn[j]
                if invert
                    Vn[j] = 1/Vn[j]
                end
                j+=1
            end
        end
    end

    # if there is only one frequency mode return the sparse matrix generated
    # from these vectors. if multiple modes then duplicate along the diagonal
    # Nmodes times.
    if Nmodes == 1
        return sparse(In,Jn,Vn,(Nnodes-1),(Nnodes-1))
    else
        return diagrepeat(sparse(In,Jn,Vn,(Nnodes-1),(Nnodes-1)),Nmodes)
    end
end

"""
    pushval!(V::Vector, val, c, invert::Bool)

Append the value `val` of capacitance or conductance to the vector `V`. Scale
the value by `c`. If `invert = true`, append `c/val` otherwise append `c*val`.

# Examples
```jldoctest
julia> V = Array{Float64, 1}(undef, 0);JosephsonCircuits.pushval!(V,2.0,-1.0,false);V
1-element Vector{Float64}:
 -2.0

julia> V = Array{Float64, 1}(undef, 0);JosephsonCircuits.pushval!(V,2.0,-1.0,true);V
1-element Vector{Float64}:
 -0.5
```
"""
function pushval!(V::Vector, val, c, invert::Bool)
    if invert
        push!(V,c/val)
    else
        push!(V,c*val)
    end
    return nothing
end
