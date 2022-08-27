
"""
    CircuitMatrices((Cnm,Gnm,Lb,Lbm,Ljb,Ljbm,Mb,invLnm,Rbnm,portdict,
        resistordict,Lmean)

A simple structure to hold the circuit matrices including the capacitance
matrix, the conductance matrix, the inductance vectors, the Josephson inductance
vectors, the mutual inductance matrix, the inverse inductance matrix,
the incidence matrix, the dictionary of port and resistor values where the nodes
are the keys and the values are the values, and the mean of the inductances. 
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
    Rbnm::SparseMatrixCSC{Int64, Int64}
    portindices::Vector{Int64}
    portnumbers::Vector{Int64}
    portimpedanceindices::Vector{Int64}
    noiseportimpedanceindices::Vector{Int64}
    Lmean
    vvn
end


"""
    symbolicmatrices(circuit;Nmodes=1,sorting=:number)

Return the symbolic matrices describing the circuit properties. 

See also  [`CircuitMatrices`](@ref), [`numericmatrices`](@ref), [`calcCn`](@ref), 
[`calcGn`](@ref), [`calcLb`](@ref),[`calcLjb`](@ref), [`calcMb`](@ref),
[`calcinvLn`](@ref), [`componentdictionaryP`](@ref),
[`componentdictionaryR`](@ref), [`calcLmean`](@ref), and [`diagrepeat`](@ref).


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
JosephsonCircuits.CircuitMatrices(sparse([1, 2, 1, 2], [1, 1, 2, 2], SymbolicUtils.Symbolic{Real}[Cc, -Cc, -Cc, Cc + Cj], 2, 2), sparse([1], [1], SymbolicUtils.Symbolic{Real}[1 / Rleft], 2, 2), 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries, 2-element SparseArrays.SparseVector{Nothing, Int64} with 0 stored entries,   [2]  =  Lj,   [2]  =  Lj, sparse(Int64[], Int64[], Nothing[], 2, 2), sparse(Int64[], Int64[], Nothing[], 2, 2), sparse([1, 2], [1, 2], [1, 1], 2, 2), [1], [1], [3], Int64[], Lj, Any[1, Ipump, Rleft, Cc, Lj, Cj])
```
"""
function symbolicmatrices(circuit; Nmodes = 1, sorting = :number)
    return numericmatrices(circuit, Dict(), Nmodes = Nmodes, sorting = sorting)
end

function symbolicmatrices(psc::ParsedSortedCircuit, cg::CircuitGraph; Nmodes = 1,
    sorting = :number)
    return numericmatrices(psc, cg, Dict(), Nmodes = Nmodes)
end

"""
    numericmatrices(circuit;Nmodes=1,sorting=:number)

Return the symbolic matrices describing the circuit properties. 

See also  [`CircuitMatrices`](@ref), [`numericmatrices`](@ref), [`calcCn`](@ref), 
[`calcGn`](@ref), [`calcLb`](@ref),[`calcLjb`](@ref), [`calcMb`](@ref),
[`calcinvLn`](@ref), [`componentdictionaryP`](@ref),
[`componentdictionaryR`](@ref), [`calcLmean`](@ref), and [`diagrepeat`](@ref).


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
    vvn = valuevectortonumber(psc.valuevector, circuitdefs)
    
    # capacitance matrix
    Cnm = calcCn(psc.typevector, psc.nodeindexarraysorted, vvn, Nmodes, psc.Nnodes)

    # conductance matrix
    Gnm = calcGn(psc.typevector, psc.nodeindexarraysorted, vvn, Nmodes, psc.Nnodes)

    # branch inductance vector with Nmodes = 1
    Lb = calcLb(psc.typevector, psc.nodeindexarraysorted, vvn, cg.edge2indexdict,
        1, cg.Nbranches)

    # branch inductance vector with Nmodes = Nmodes
    Lbm = calcLb(psc.typevector, psc.nodeindexarraysorted, vvn, cg.edge2indexdict,
        Nmodes, cg.Nbranches)

    # branch Josephson inductance vector with Nmodes = 1
    Ljb = calcLjb(psc.typevector, psc.nodeindexarraysorted, vvn, cg.edge2indexdict,
        1, cg.Nbranches)

    # branch Josephson inductance vector with Nmodes = Nmodes
    Ljbm = calcLjb(psc.typevector, psc.nodeindexarraysorted, vvn, cg.edge2indexdict,
        Nmodes, cg.Nbranches)

    # mutual branch inductance matrix
    Mb = calcMb(psc.typevector, psc.nodeindexarraysorted, vvn, psc.namedict,
        psc.mutualinductorvector, cg.edge2indexdict, 1, cg.Nbranches)

    # inverse nodal inductance matrix from branch inductance vector and branch
    # inductance matrix
    invLnm = calcinvLn(Lb, Mb, cg.Rbn, Nmodes)

    # expand the size of the incidence matrix
    Rbnm = diagrepeat(cg.Rbn, Nmodes)

    # calculate Lmean
    Lmean = calcLmean(psc.typevector, vvn)

    portindices, portnumbers = calcportindicesnumbers(psc.typevector,
        psc.nodeindexarraysorted, psc.mutualinductorvector, vvn)

    portimpedanceindices = calcportimpedanceindices(psc.typevector,
        psc.nodeindexarraysorted, psc.mutualinductorvector, vvn)

    noiseportimpedanceindices = calcnoiseportimpedanceindices(psc.typevector,
        psc.nodeindexarraysorted, psc.mutualinductorvector, vvn)

    return CircuitMatrices(Cnm, Gnm, Lb, Lbm, Ljb, Ljbm, Mb, invLnm, Rbnm,
        portindices, portnumbers, portimpedanceindices,
        noiseportimpedanceindices, Lmean, vvn)
end


"""
    calcIb(typevector, nodeindexarray, valuevector, edge2indexdict, Nmodes, Nbranches)

Calculate the sparse branch current source vector whose length is Nbranches*Nmodes.
Note that the nodeindexarray is "one indexed" so 1 is the ground node. 

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
typevector = [:I,:C,:L,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [1e-9, 0.2, 4e-9, 1e-12]
namedict = Dict{Symbol, Int64}(:C2 => 4,:L1 => 3,:I1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ib = JosephsonCircuits.calcIb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 1 stored entry:
  [1]  =  1.0e-9
```
```jldoctest
@variables I1 C1 L1 C2
Nmodes = 1
Nbranches = 2
typevector = [:I,:C,:L,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [I1, C1, L1, C2]
namedict = Dict{Symbol, Int64}(:C2 => 4,:L1 => 3,:I1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ib = JosephsonCircuits.calcIb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 1 stored entry:
  [1]  =  I1
```
"""
function calcIb(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(typevector, nodeindexarray, valuevector,
        calcvaluetype(typevector, valuevector, [:I]), edge2indexdict, Nmodes,
        Nbranches, :I)
end


"""
    calcVb(typevector, nodeindexarray, valuevector, edge2indexdict, Nmodes, Nbranches)

Calculate the sparse branch voltage source vector whose length is Nbranches*Nmodes.
Note that the nodeindexarray is "one indexed" so 1 is the ground node. 

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
typevector = [:V,:C,:L1,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [1e-9, 0.2, 4e-9, 1e-12]
namedict = Dict{Symbol, Int64}(:C2 => 4,:L1 => 3,:V1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Vb = JosephsonCircuits.calcVb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 1 stored entry:
  [1]  =  1.0e-9
```
```jldoctest
@variables V1 C1 L1 C2
Nmodes = 1
Nbranches = 2
typevector = [:V,:C,:L,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [V1, C1, L1, C2]
namedict = Dict{Symbol, Int64}(:C2 => 4,:L1 => 3,:V1 => 1,:C1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Vb = JosephsonCircuits.calcVb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 1 stored entry:
  [1]  =  V1
```
"""
function calcVb(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(typevector, nodeindexarray, valuevector,
        calcvaluetype(typevector, valuevector, [:V]), edge2indexdict, Nmodes,
        Nbranches, :V)
end

"""
    calcLb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

Calculate the sparse branch inductance vector whose length is Nbranches*Nmodes.
Note that the nodeindexarray is "one indexed" so 1 is the ground node. 

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
typevector = [:L,:K,:L,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [1e-9, 0.2, 4e-9, 1e-12]
namedict = Dict{Symbol, Int64}(:C2 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Lb = JosephsonCircuits.calcLb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 2 stored entries:
  [1]  =  1.0e-9
  [2]  =  4.0e-9
```
```jldoctest
@variables L1 K1 L2 C1
Nmodes = 1
Nbranches = 2
typevector = [:L,:K,:L,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [L1, K1, L2, C1]
namedict = Dict{Symbol, Int64}(:C1 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Lb = JosephsonCircuits.calcLb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 2 stored entries:
  [1]  =  L1
  [2]  =  L2
```
"""
function calcLb(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(typevector, nodeindexarray, valuevector,
        calcvaluetype(typevector, valuevector, [:L,:K]), edge2indexdict,
        Nmodes, Nbranches, :L)
end

"""
    calcLjb(typevector, nodeindexarray, valuevector, edge2indexdict, Nmodes, Nbranches)

Calculate the sparse branch Josephson inductance vector whose length is
Nbranches*Nmodes. Note that the nodeindexarray is "one indexed" so 1 is the
ground node.

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
typevector = [:Lj,:C,:Lj,:C]
nodeindexarray = [2 3 3 3; 1 2 1 1]
valuevector = [1e-9, 1e-12, 4e-9, 1e-12]
namedict = Dict{Symbol, Int64}(:C2 => 4,:L2 => 3,:L1 => 1,:Cc => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ljb = JosephsonCircuits.calcLjb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Float64, Int64} with 2 stored entries:
  [1]  =  1.0e-9
  [2]  =  4.0e-9
```
```jldoctest
@variables Lj1 K1 Lj2 C1
Nmodes = 1
Nbranches = 2
typevector = [:Lj,:K,:Lj,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [Lj1, K1, Lj2, C1]
namedict = Dict{Symbol, Int64}(:C1 => 4,:Lj2 => 3,:Lj1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
Ljb = JosephsonCircuits.calcLjb(typevector,nodeindexarray,valuevector,edge2indexdict,Nmodes,Nbranches)

# output
2-element SparseArrays.SparseVector{Num, Int64} with 2 stored entries:
  [1]  =  Lj1
  [2]  =  Lj2
```
"""
function calcLjb(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, edge2indexdict::Dict, Nmodes, Nbranches)
    return calcbranchvector(typevector, nodeindexarray, valuevector,
        calcvaluetype(typevector, valuevector, [:Lj]), edge2indexdict, Nmodes,
        Nbranches, :Lj)
end

"""
    calcbranchvector(typevector, nodeindexarray, valuevector, edge2indexdict,
        Nmodes, Nbranches)

Calculate the sparse branch vector whose length is Nbranches*Nmodes for the 
given component symbol. Note that the nodeindexarray is "one indexed" so 1 is
the ground node. 
"""
function calcbranchvector(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, valuetypevector::Vector, edge2indexdict::Dict,
    Nmodes, Nbranches, component::Symbol)

    # define empty vectors of zero length for the indices and values
    Ib = Array{Int64, 1}(undef, 0)
    Vb = Array{eltype(valuetypevector), 1}(undef, 0)

    @inbounds for (i,type) in enumerate(typevector)
        if type == component
            push!(Ib,edge2indexdict[(nodeindexarray[1,i],nodeindexarray[2,i])])
            push!(Vb,valuevector[i])
        end
    end
    if Nmodes == 1
        return sparsevec(Ib,Vb,Nbranches)
    else
        #define empty vectors for the row indices and values
        Ibm = Array{eltype(Ib), 1}(undef, length(Ib)*Nmodes)
        Vbm = Array{eltype(valuetypevector), 1}(undef, length(Vb)*Nmodes)

        # repeat the elements Nmodes times
        n = 1
        @inbounds for i in eachindex(Ib)
            for j = 1:Nmodes
                Ibm[n] = (Ib[i]-1)*Nmodes+j
                Vbm[n] = Vb[i]
                n+=1
            end
        end

        # make the sparse array
        return sparsevec(Ibm,Vbm,Nbranches*Nmodes)

    end
end


"""
    calcMb(typevector, nodeindexarray, valuevector, namedict,
        mutualinductorvector, edge2indexdict, Nmodes, Nbranches)

Returns the branch mutual inductance matrix. Note that the
nodeindexarray is "one indexed" so 1 is the ground node.

# Examples
```jldoctest
Nmodes = 1
Nbranches = 2
typevector = [:L,:K,:L,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [1e-9, 0.2, 2e-9, 1e-12]
namedict = Dict{Symbol, Int64}(:C2 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
mutualinductorvector = [ :L1, :L2]
Mb = JosephsonCircuits.calcMb(typevector,nodeindexarray,valuevector,namedict,mutualinductorvector,edge2indexdict,Nmodes,Nbranches)

# output
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 2 stored entries:
  ⋅           2.82843e-10
 2.82843e-10   ⋅ 
```
```jldoctest
@variables L1 L2 K1 C1
Nmodes = 1
Nbranches = 2
typevector = [:L,:K,:L,:C]
nodeindexarray = [2 0 3 3; 1 0 1 1]
valuevector = [L1, K1, L2, C1]
namedict = Dict{Symbol, Int64}(:C1 => 4,:L2 => 3,:L1 => 1,:K1 => 2)
edge2indexdict = Dict{Tuple{Int64, Int64}, Int64}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
mutualinductorvector = [ :L1, :L2]
Mb = JosephsonCircuits.calcMb(typevector,nodeindexarray,valuevector,namedict,mutualinductorvector,edge2indexdict,Nmodes,Nbranches)

# output
2×2 SparseArrays.SparseMatrixCSC{Num, Int64} with 2 stored entries:
              ⋅  K1*sqrt(L1*L2)
 K1*sqrt(L1*L2)               ⋅
```
"""

function calcMb(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, namedict::Dict, mutualinductorvector::Vector,
    edge2indexdict::Dict, Nmodes, Nbranches)
    return calcMb_inner(typevector,nodeindexarray,valuevector,
        calcvaluetype(typevector,valuevector,[:L,:K]),namedict,
        mutualinductorvector,edge2indexdict,Nmodes,Nbranches)
end

function calcMb_inner(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, valuetypevector::Vector, namedict::Dict,
    mutualinductorvector::Vector, edge2indexdict::Dict, Nmodes, Nbranches)

    # define empty vectors of zero length for the row indices, column indices,
    # and values
    Ib = Array{Int64, 1}(undef, 0)
    Jb = Array{Int64, 1}(undef, 0)
    Vb = Array{eltype(valuetypevector), 1}(undef, 0)

    n = 1
    #loop through typevector for mutual inductors
    @inbounds for (i,type) in enumerate(typevector)
        # when we find a mutual inductor:
        #  find the value of the mutual inductor in valuevector[i]
        #  find the names of the two inductors it couples together from mutualinductorvector[n]
        #  look up the index of the inductors in index=namedict[inductorsymbol] for each inductor symbol
        #  given the index of the inductor, look of the value of the inductor from
        #     valuevector
        # then compute the value of the mutual inductance from the two inductor values and K

        # then use the index of the inductors to get the nodes from nodeindexarray
        # use that as a key in edge2indexdict to look up the branch index
        # then assign those to I, J, V for the sparse array. 
        # then do the usual step of expanding that to Nmodes after finishing this loop

        if type == :K
            # value of K
            K = valuevector[i]

            # names of inductors
            inductor1name = mutualinductorvector[2*n-1]
            inductor2name = mutualinductorvector[2*n]

            # indices of inductors
            inductor1index = namedict[inductor1name]
            inductor2index = namedict[inductor2name]

            # values of inductors
            inductor1value = valuevector[inductor1index]
            inductor2value = valuevector[inductor2index]

            # values of mutual inductance Lm
            Lm = K*sqrt(inductor1value*inductor2value)

            inductor1edge = (nodeindexarray[1,inductor1index],nodeindexarray[2,inductor1index])
            inductor2edge = (nodeindexarray[1,inductor2index],nodeindexarray[2,inductor2index])

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

    if Nmodes == 1
        return sparse(Ib,Jb,Vb,Nbranches,Nbranches)
    else
        # define empty vectors for the row indices, column indices, and values
        Ibm = Array{eltype(Ib), 1}(undef, length(Ib)*Nmodes)
        Jbm = Array{eltype(Jb), 1}(undef, length(Jb)*Nmodes)
        Vbm = Array{eltype(valuetypevector), 1}(undef, length(Vb)*Nmodes)

        # repeat the elements along the diagonal Nmodes times
        n = 1
        @inbounds for i in eachindex(Ib)
            for j = 1:Nmodes
                Ibm[n] = (Ib[i]-1)*Nmodes+j
                Jbm[n] = (Jb[i]-1)*Nmodes+j
                Vbm[n] = Vb[i]
                n+=1
            end
        end

        # make the sparse array        
        return sparse(Ibm,Jbm,Vbm,Nbranches*Nmodes,Nbranches*Nmodes)

    end
end


"""
    calcinvLn(Lb, Rbn, Nmodes)

Returns the nodal inverse inductance matrix. Accepts the vector of branch
inductances Lb and the incidence matrix Rbn.

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
"""
function calcinvLn(Lb::SparseVector, Rbn::SparseMatrixCSC, Nmodes)
    if nnz(Lb)>0
        # s = transpose(Rbn[Lb.nzind,:])*sparse(Diagonal(1 ./Lb.nzval))*Rbn[Lb.nzind,:]
        s = transpose(Rbn[Lb.nzind,:])*spdiagm(0 => 1 ./Lb.nzval)*Rbn[Lb.nzind,:]
        if Nmodes > 1
            return diagrepeat(s,Nmodes)
        else
            return s
        end 
    else
        return spzeros(eltype(Lb),Nmodes*size(Rbn)[2],Nmodes*size(Rbn)[2])
    end
end

"""
    calcinvLn(Lb, Mb, Rbn)

Returns the nodal inverse inductance matrix. Accepts the vector of branch
inductances Lb, the branch mutual inductance matrix Mb, and the incidence
matrix Rbn.

Using ldiv instead of an inverse: (where the extra div is an escape sequence)
Can solve A x = B with: x = A \\ B or x = invA * B, so we can perform the
inverse here with:
s = RbnT * invL * Rbn or s = RbnT * (L \\ Rbn), the latter of which should be
faster and more numerically stable.


# Examples
```jldoctest
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[1e-9,4e-9])
Mb = JosephsonCircuits.SparseArrays.sparse([2,1], [1,2], [4e-10,4e-10])
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
JosephsonCircuits.calcinvLn(Lb,Mb,Rbn,Nmodes)

# output
2×2 SparseArrays.SparseMatrixCSC{Float64, Int64} with 4 stored entries:
  1.04167e9  -1.04167e8
 -1.04167e8   2.60417e8
```
```jldoctest
@variables L1 L2 Lm
Nmodes = 1
Lb = JosephsonCircuits.SparseArrays.sparsevec([1,2],[L1,L2]);
Mb = JosephsonCircuits.SparseArrays.sparse([2,1], [1,2], [Lm,Lm]);
Rbn = JosephsonCircuits.SparseArrays.sparse([1,2], [1,2], [1,1])
println(JosephsonCircuits.calcinvLn(Lb,Mb,Rbn,Nmodes))

# output
sparse([1, 2, 1, 2], [1, 1, 2, 2], Num[(1 + (Lm*(Lm / L1)) / (L2 + (-(Lm^2)) / L1)) / L1, (-(Lm / L1)) / (L2 + (-(Lm^2)) / L1), (-(Lm / (L2 + (-(Lm^2)) / L1))) / L1, 1 / (L2 + (-(Lm^2)) / L1)], 2, 2)
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
    Rbn::SparseMatrixCSC, Nmodes, valuetypevector::Vector)
    if nnz(Lb) > 0 &&  nnz(Mb) > 0
            # add the mutual inductance matrix to a diagonal matrix made from the
            # inductance vector. 
            # we pick out only the indices where there are inductors for
            # efficiency reasons. 

        if eltype(valuetypevector) <: Symbolic

            # take a subset of the arrays
            Mbs = Mb[Lb.nzind,Lb.nzind]
            Lbs = Lb[Lb.nzind]

            # add together the two sparse arrays

            # define empty vectors for the rows, columns, and values
            In = Array{Int64, 1}(undef, nnz(Mbs)+nnz(Lbs))
            Jn = Array{Int64, 1}(undef, nnz(Mbs)+nnz(Lbs))
            Vn = Vector{eltype(valuetypevector)}(undef, nnz(Mbs)+nnz(Lbs))

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

            # println(eltype(valuetypevector)
            # println("type: ",eltype(valuetypevector) <: Symbolic)


            L = Array{Any,2}(undef,size(Lsparse))
            fill!(L,0)

            for i = 1:length(Lsparse.colptr)-1
                for j in Lsparse.colptr[i]:(Lsparse.colptr[i+1]-1)
                    L[Lsparse.rowval[j],i] = Lsparse.nzval[j]
                end
            end

        else
            L = Mb[Lb.nzind,Lb.nzind] + Diagonal(Lb[Lb.nzind])
        end

            ## using the inverse directly (bad)
            # s = transpose(Rbn[Lb.nzind,:])*sparse(inv(Matrix(L)))*Rbn[Lb.nzind,:]
            if eltype(valuetypevector) <: Union{AbstractFloat, Complex}
                # using ldiv with klu. fastest option in most cases
                s = transpose(Rbn[Lb.nzind,:])*sparse(KLU.klu(L) \ Matrix(Rbn[Lb.nzind,:]))
                
                # if above is not allowed, use ldiv with a dense matrix. this can be much
                # slower but will work for symbolic matrices. 
                # s = transpose(Rbn[Lb.nzind,:])*sparse(Matrix(L) \ Matrix(Rbn[Lb.nzind,:]))
            else

                # println(L)
                # println(typeof(L))
                # println(typeof(transpose(Rbn[Lb.nzind,:])))
                # println( typeof(Matrix(Rbn[Lb.nzind,:])))
                # s = transpose(Rbn[Lb.nzind,:])*(Matrix(L) \ Matrix(Rbn[Lb.nzind,:]))
                # s = transpose(Rbn[Lb.nzind,:])*(Matrix(L) \ Matrix(Rbn[Lb.nzind,:]))

                try
                    s = calcsymbolicinvLn(L,Lb,Rbn)
                catch e
                    # throw(MethodError("Calculating the inverse inductance matrix with calcsymbolicinvLn(L,Lb,Rbn) when there are mutual inductors requires Symbolics.jl. Run the command:  using Symbolics"))
                    println("Calculating the inverse inductance matrix with calcsymbolicinvLn(L,Lb,Rbn) when there are mutual inductors requires Symbolics.jl. Run the command:  using Symbolics")
                    throw(e)
                end

                # s =  sparse(transpose(Rbn[Lb.nzind,:])*(sym_lu(L)\ Matrix(Rbn[Lb.nzind,:])))
                # s = SparseMatrixCSC(s.m, s.n, s.colptr, s.rowval,value.(s.nzval))
            end 

            if Nmodes > 1
                return diagrepeat(s,Nmodes)
            else
                return s
            end 
    elseif nnz(Lb) > 0
        return calcinvLn(Lb,Rbn,Nmodes)
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
    calcLmean(typevector,valuevector)

Return the mean of the linear and Josephson inductors. 

# Examples
```jldoctest
julia> JosephsonCircuits.calcLmean([:R,:L,:C,:Lj],[10,4,5,1])
2.5

julia> @variables R1 L1 C1 Lj1;JosephsonCircuits.calcLmean([:R,:L,:C,:Lj],[R1, L1, C1, Lj1])
(1//2)*L1 + (1//2)*Lj1
```
"""
function calcLmean(typevector::Vector{Symbol}, valuevector::Vector)
    return calcLmean_inner(typevector, valuevector,
        calcvaluetype(typevector, valuevector, [:Lj, :L]))
end
function calcLmean_inner(typevector::Vector, valuevector::Vector,
    valuetypevector::Vector)

    if length(typevector) !== length(valuevector)
        throw(DimensionMismatch("typevector and valuevector should have the same length"))
    end

    # it's a litle absurd but we have to copy the inductance values into
    # a new array to perform a type stable mean over some elements of
    # valuevector
    Vn = Array{eltype(valuetypevector), 1}(undef, 0)
    @inbounds for (i,type) in enumerate(typevector)
        if type == :L || type == :Lj
            push!(Vn,valuevector[i])
        end
    end

    # take the mean
    # Lmean = zero(eltype(valuetypevector))
    Lmean = 0
    ninductors = 0
    for v in Vn
        ninductors+=1
        Lmean = Lmean + (v-Lmean)/ninductors
    end
    return Lmean
end


"""
    calcCn(typevector, nodeindexarray, valuevector, Nmodes, Nnodes)

Returns the node capacitance matrix from the capacitance values in valuevector
when typevector has the symbol :C with node indices from nodeindexarray.  Other symbols are 
ignored. Capacitances to ground become diagonal elements. Capacitance between
elements is an off-diagonal element with a minus sign and is added to the
diagonal with a  plus sign. The dimensions of the output are (Nnodes-1)*Nmodes
by (Nnodes-1) times Nmodes where Nnodes is the number of nodes including 
ground and Nmodes is the number of different frequencies. Note that the 
nodeindexarray is "one indexed" so 1 is the ground node. 

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
function calcCn(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, Nmodes, Nnodes)
    return calcnodematrix(typevector, nodeindexarray, valuevector,
        calcvaluetype(typevector, valuevector, [:C]), Nmodes, Nnodes, :C, false)
end

"""
    calcGn(typevector, nodeindexarray, valuevector, Nmodes, Nnodes)

Returns the node conductance matrix from the resistance values in valuevector
when typevector has the symbol :R. The node indices are taken from nodeindexarray.
Conductances to ground are diagonal elements. Conductance between elements is
an off-diagonal element with a minus sign and is added to the diagonal with a
plus sign. The dimensions of the output are (Nnodes-1) times Nmodes by
(Nnodes-1) times Nmodes. Note that the  nodeindexarray is "one indexed" so 1
is the ground node.

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
function calcGn(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, Nmodes, Nnodes)

    return calcnodematrix(typevector, nodeindexarray, valuevector,
        calcvaluetype(typevector, valuevector, [:R]), Nmodes, Nnodes, :R, true)
end


"""
    calcnodematrix(typevector, nodeindexarray, valuevector, Nmodes, Nnodes,
        component, invert)

Returns either the capacitance or conductance matrix depending on the values
of component and invert. :C and false for capacitance and :R and true for
conductance. The dimensions of the output are (Nnodes-1) times Nmodes by 
(Nnodes-1) times Nmodes. Note that the  nodeindexarray is "one indexed" so 1
is the ground node. 

"""
function calcnodematrix(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int64},
    valuevector::Vector, valuetypevector::Vector, Nmodes, Nnodes,
    component::Symbol, invert::Bool)

    if length(typevector) !== size(nodeindexarray,2) || length(typevector) !== length(valuevector)
        throw(DimensionMismatch("typevector, nodeindexarray, and valuevector should have the same length"))
    end

    if size(nodeindexarray,1) !== 2
        throw(DimensionMismatch("nodeindexarray should have a first dimension size of 2."))
    end

    # define empty vectors of zero length for the row indices, column indices,
    # and values
    In = Array{Int64, 1}(undef, 0)
    Jn = Array{Int64, 1}(undef, 0)
    Vn = Array{eltype(valuetypevector), 1}(undef, 0)

    # generate the capacitance or conductance matrix values for Nmodes=1
    @inbounds for (i,type) in enumerate(typevector)
        if type == component

            if nodeindexarray[1,i] == 1
                # capacitance to ground, add to diagonal
                push!(In,nodeindexarray[2,i])
                push!(Jn,nodeindexarray[2,i])
                pushval!(Vn,valuevector[i],1,invert)

            elseif nodeindexarray[2,i] == 1
                # capacitance to ground, add to diagonal
                push!(In,nodeindexarray[1,i])
                push!(Jn,nodeindexarray[1,i])
                pushval!(Vn,valuevector[i],1,invert)

            else
                # diagonal elements
                push!(In,nodeindexarray[1,i])
                push!(Jn,nodeindexarray[1,i])
                pushval!(Vn,valuevector[i],1,invert)

                push!(In,nodeindexarray[2,i])
                push!(Jn,nodeindexarray[2,i])
                pushval!(Vn,valuevector[i],1,invert)

                # off diagonal elements
                push!(In,nodeindexarray[1,i])
                push!(Jn,nodeindexarray[2,i])
                pushval!(Vn,valuevector[i],-1,invert)

                push!(In,nodeindexarray[2,i])
                push!(Jn,nodeindexarray[1,i])
                pushval!(Vn,valuevector[i],-1,invert)
            end
        end
    end

    if Nmodes == 1
        return sparse(In.-1,Jn.-1,Vn,(Nnodes-1),(Nnodes-1))
    else
        # define empty vectors for the row indices, column indices, and values
        Inm = Array{eltype(In), 1}(undef, length(In)*Nmodes)
        Jnm = Array{eltype(Jn), 1}(undef, length(Jn)*Nmodes)
        Vnm = Array{eltype(valuetypevector), 1}(undef, length(Vn)*Nmodes)

        # repeat the elements along the diagonal Nmodes times
        n = 1
        @inbounds for i in eachindex(In)
            for j = 1:Nmodes
                Inm[n] = (In[i]-2)*Nmodes+j
                Jnm[n] = (Jn[i]-2)*Nmodes+j
                Vnm[n] = Vn[i]
                n+=1
            end
        end

        # make the sparse array        
        return sparse(Inm,Jnm,Vnm,(Nnodes-1)*Nmodes,(Nnodes-1)*Nmodes)
    end
end


"""
    pushval!(V, val::Number, c::Number, invert::Bool)

Append the value val of capacitance or conductance to the vector V. Scale the
value by c. If invert is true, append c/val otherwise append c*val.

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
