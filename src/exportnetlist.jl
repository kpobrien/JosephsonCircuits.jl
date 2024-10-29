

"""
    export_netlist(filename, circuit, circuitdefs)

Export the netlist in `circuit` to the file with name and path `filename`.
"""
function export_netlist(filename, circuit, circuitdefs)
    open(filename, "w") do io
        export_netlist!(io, circuit, circuitdefs)
    end
    return nothing
end

"""
    export_netlist(filename, circuit)

Export the netlist in `circuit` to the file with name and path `filename`.
"""
function export_netlist(filename, circuit)
    return export_netlist(filename, circuit, Dict())
end

"""
    export_netlist!(io::IO, circuit, circuitdefs)

Export the netlist in `circuit` to the IOBuffer or IOStream `io`.

# Examples
```jldoctest
julia> io = IOBuffer();JosephsonCircuits.export_netlist!(io, [("P","1","0",1),("R","1","0",50.0)],Dict());println(String(take!(io)))
P 1 0 1
R 1 0 50.0
```
"""
function export_netlist!(io::IO, circuit, circuitdefs)
    for c in circuit
        for (i,ci) in enumerate(c)
            if i > 1
                write(io," ")
            end
            write(io,string(JosephsonCircuits.Symbolics.substitute(ci,circuitdefs)))
        end
        write(io,"\n")
    end
end

"""
    import_netlist(filename)

Import the netlist from the IOBuffer or IOStream `io` and return a vector of
tuples representing the circuit.
"""
function import_netlist(filename)
    circuit = Tuple{String,String,String,Num}[]
    open(filename, "r") do io
        import_netlist!(io, circuit)
    end
    return circuit
end

"""
    import_netlist!(io::IO, circuit)

Import the netlist from the IOBuffer or IOStream `io` to the vector of tuples
`circuit`.

# Examples
```jldoctest
julia> io = IOBuffer();circuit1=[("P","1","0",1),("R","1","0",50.0)];JosephsonCircuits.export_netlist!(io,circuit1,Dict());circuit2 = Tuple{String,String,String,Num}[];JosephsonCircuits.import_netlist!(io,circuit2);circuit2
2-element Vector{Tuple{String, String, String, Num}}:
 ("P", "1", "0", 1.0)
 ("R", "1", "0", 50.0)

julia> io = IOBuffer("P 1 0 1\\n R 1 0 50.0");circuit2 = Tuple{String,String,String,Num}[];JosephsonCircuits.import_netlist!(io,circuit2);circuit2
2-element Vector{Tuple{String, String, String, Num}}:
 ("P", "1", "0", 1.0)
 ("R", "1", "0", 50.0)
```
"""
function import_netlist!(io::IO, circuit)
    seekstart(io)
    for line in eachline(io)
        split_line = split(strip(line),r"\s+")
        if length(split_line) != 4
            error("each line should have component name, node1, node2, component value")
        end
        value = try
            parse(Float64,split_line[4])
        catch
            # https://docs.sciml.ai/Symbolics/stable/manual/parsing/
            # JosephsonCircuits.Symbolics.parse_expr_to_symbolic(Meta.parse(split_line[4]),Main)
            JosephsonCircuits.Symbolics.parse_expr_to_symbolic(Meta.parse(split_line[4]),@__MODULE__)
        end
        push!(circuit,(split_line[1],split_line[2],split_line[3],value))
    end
    return nothing
end



# export_netlist("test1.net", circuit,circuitdefs)
# export_netlist("test2.net", circuit,Dict())

# circuit = import_netlist("test2.net")



"""
    sumvalues(type::Symbol, value1, value2)

Sum together two values in different ways depending on the circuit component
type.

# Examples
```jldoctest
julia> JosephsonCircuits.sumvalues(:L, 1.0, 4.0)
0.8

julia> JosephsonCircuits.sumvalues(:Lj, 1.0, 4.0)
0.8

julia> JosephsonCircuits.sumvalues(:C, 1.0, 4.0)
5.0

julia> JosephsonCircuits.sumvalues(:K, 1.0, 4.0)
5.0
```
"""
function sumvalues(type::Symbol, value1, value2)
    if type == :C || type == :K
        return value1+value2
    elseif type == :Lj || type == :L
        return 1/(1/value1+1/value2)
    else
        error("unknown component type in sumvalues")
    end
end

"""
    calcnodes(nodeindex::Int, mutualinductorindex::Int,
        componenttypes::Vector{Symbol}, nodeindexarray::Matrix,
        componentnamedict::Dict, mutualinductorbranchnames::Vector{String})

Calculate the two nodes (or mutual inductor indices) given the index in the
typvector and the component type. For component types where order matters,
such as mutual inductors, the nodes are not sorted. For other component types
where order does not matter, the nodes are sorted. 

# Examples
```jldoctest
@variables R Cc L1 L2 Cj1 Cj2 I1 V1
@variables Ipump Rleft L1 K1 K2 L2 C2 C3
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1))
push!(circuit,("K1","L1","L2",K1))
push!(circuit,("K2","L1","L2",K2))
push!(circuit,("L2","2","0",L2))
push!(circuit,("C2","2","0",C2))
push!(circuit,("C3","2","0",C3))
psc = JosephsonCircuits.parsesortcircuit(circuit)
println(JosephsonCircuits.calcnodes(1,1,psc.componenttypes,psc.nodeindices, psc.componentnamedict,psc.mutualinductorbranchnames))
println(JosephsonCircuits.calcnodes(5,1,psc.componenttypes,psc.nodeindices, psc.componentnamedict,psc.mutualinductorbranchnames))

# output
(1, 2)
(4, 7)
```
"""
function calcnodes(nodeindex::Int, mutualinductorindex::Int,
    componenttypes::Vector{Symbol}, nodeindexarray::Matrix,
    componentnamedict::Dict, mutualinductorbranchnames::Vector{String})

    # calculate the nodes
    if componenttypes[nodeindex] == :K
        # don't sort these because the mutual inductance changes sign
        # if the nodes are changed. this is OK because only values with
        # the same inductor ordering will be summed. 

        # names of inductors
        inductor1name = mutualinductorbranchnames[2*mutualinductorindex-1]
        inductor2name = mutualinductorbranchnames[2*mutualinductorindex]

        # indices of inductors
        return componentnamedict[inductor1name], componentnamedict[inductor2name]

    else
        if nodeindexarray[1,nodeindex] < nodeindexarray[2,nodeindex]
            return nodeindexarray[1,nodeindex], nodeindexarray[2,nodeindex]
        else
            return nodeindexarray[2,nodeindex], nodeindexarray[1,nodeindex]
        end
    end
end

"""
    componentdictionaries(componenttypes::Vector{Symbol},
        nodeindexarray::Matrix{Int}, componentnamedict::Dict,
        mutualinductorbranchnames::Vector)

# Examples
```jldoctest
@variables Ipump Rleft L1 K1 L2 C2 C3
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1))
push!(circuit,("K1","L1","L2",K1))
push!(circuit,("L2","2","0",L2))
push!(circuit,("C2","2","0",C2))
push!(circuit,("C3","2","0",C3))
psc = parsesortcircuit(circuit)
countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.componenttypes,psc.nodeindices,psc.componentnamedict,psc.mutualinductorbranchnames)

println(countdict)
println(indexdict)

# output
Dict((:L, 1, 3) => 1, (:K, 4, 6) => 1, (:R, 1, 2) => 1, (:I, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 2, (:L, 1, 2) => 1)
Dict((:C, 1, 3, 1) => 7, (:I, 1, 2, 1) => 2, (:R, 1, 2, 1) => 3, (:L, 1, 3, 1) => 6, (:C, 1, 3, 2) => 8, (:L, 1, 2, 1) => 4, (:P, 1, 2, 1) => 1, (:K, 4, 6, 1) => 5)
```
```jldoctest
@variables Ipump Rleft L1 K1 K2 L2 C2 C3
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("I1","1","0",Ipump))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1))
push!(circuit,("K1","L1","L2",K1))
push!(circuit,("K2","L1","L2",K2))
push!(circuit,("L2","2","0",L2))
push!(circuit,("C2","2","0",C2))
push!(circuit,("C3","2","0",C3))
psc = parsesortcircuit(circuit)
countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.componenttypes,psc.nodeindices,psc.componentnamedict,psc.mutualinductorbranchnames)

println(countdict)
println(indexdict)

# output
Dict((:L, 1, 3) => 1, (:K, 4, 7) => 2, (:R, 1, 2) => 1, (:I, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 2, (:L, 1, 2) => 1)
Dict((:C, 1, 3, 1) => 8, (:I, 1, 2, 1) => 2, (:R, 1, 2, 1) => 3, (:K, 4, 7, 1) => 5, (:K, 4, 7, 2) => 6, (:L, 1, 2, 1) => 4, (:L, 1, 3, 1) => 7, (:P, 1, 2, 1) => 1, (:C, 1, 3, 2) => 9)
```
"""
function componentdictionaries(componenttypes::Vector{Symbol},
    nodeindexarray::Matrix{Int}, componentnamedict::Dict,
    mutualinductorbranchnames::Vector{String})

    if  length(componenttypes) != size(nodeindexarray,2)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end

    # key = (componenttype,node1,node2), value = counts
    countdict = Dict{Tuple{eltype(componenttypes),eltype(nodeindexarray),eltype(nodeindexarray)},Int}()
    sizehint!(countdict,length(componenttypes))

    # key = (node1,node2,count), value = index in componenttypes
    indexdict = Dict{Tuple{eltype(componenttypes),eltype(nodeindexarray),eltype(nodeindexarray),Int},Int}()
    sizehint!(indexdict,length(componenttypes))

    mutualinductorindex = 0
    for i in eachindex(componenttypes)

        if componenttypes[i] == :K
            mutualinductorindex+=1
        end

        node1, node2 = calcnodes(i, mutualinductorindex, componenttypes,
            nodeindexarray, componentnamedict, mutualinductorbranchnames)

        countkey = (componenttypes[i], node1, node2)
        if haskey(countdict,countkey)
            countdict[countkey] += 1
        else
            countdict[countkey] = 1
        end

        indexkey = (componenttypes[i], node1, node2, countdict[countkey])
        indexdict[indexkey] = i
    end

    return countdict, indexdict
end

"""
    sumbranchvalues!(type::Symbol, node1::Int, node2::Int,componentvalues::Vector,
        countdict, indexdict)

Given a branch and a type, return the sum of all of the values of the same
type and branch. The sum will behave differently depending on the type.

# Examples
```jldoctest
vvn = Real[1, 50.0, 1.0e-13, 2.0e-9, 2.0e-9, 5.0e-13, 5.0e-13, 0.1]
countdict = Dict((:L, 1, 3) => 2, (:R, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 2, (:C, 2, 3) => 1, (:I, 1, 3) => 1)
indexdict = Dict((:C, 2, 3, 1) => 3, (:C, 1, 3, 1) => 6, (:R, 1, 2, 1) => 2, (:L, 1, 3, 1) => 4, (:C, 1, 3, 2) => 7, (:L, 1, 3, 2) => 5, (:P, 1, 2, 1) => 1, (:I, 1, 3, 1) => 8)
println(JosephsonCircuits.sumbranchvalues!(:C, 1, 3, vvn, countdict, indexdict))

# output
(true, 1.0e-12, 6)
```
"""
function sumbranchvalues!(type::Symbol, node1::Int, node2::Int,
    componentvalues::Vector, countdict::Dict, indexdict::Dict)
    countkey = (type, node1, node2)
    countflag = false
    value = zero(eltype(componentvalues))
    index = 0

    if haskey(countdict,countkey)
        counts = countdict[countkey]
        if counts > 0
            countflag = true
            index = indexdict[(type, node1, node2, 1)]
            value = componentvalues[index]

            for count in 2:counts
                index1 = indexdict[(type, node1, node2, count)]
                value = sumvalues(type, value, componentvalues[index1])
            end
            countdict[countkey] = 0
        end
    end

    return countflag, value, index

end

"""
    calcCjIcmean(componenttypes::Vector{Symbol}, nodeindexarray::Matrix{Int},
        componentvalues::Vector, componentnamedict::Dict,
        mutualinductorbranchnames::Vector{String}, countdict::Dict,
        indexdict::Dict)

Calculate the junction properties including the max and min critical currents
and ratios of critical current to junction capacitance. This is necessary in
order to set the junction properties of the JJ model in WRSPICE.

# Examples
```jldoctest
componenttypes = [:P, :R, :C, :Lj, :C, :C, :Lj, :C]
nodeindexarray = [2 2 2 3 3 3 4 4; 1 1 3 1 1 4 1 1]
componentvalues = Real[1, 50.0, 1.0e-13, 1.0e-9, 1.0e-12, 1.0e-13, 1.1e-9, 1.2e-12]
componentnamedict = Dict("R1" => 2, "Cc2" => 6, "Cj2" => 8, "Cj1" => 5, "P1" => 1, "Cc1" => 3, "Lj2" => 7, "Lj1" => 4)
mutualinductorbranchnames = String[]
countdict = Dict((:Lj, 1, 4) => 1, (:C, 3, 4) => 1, (:C, 1, 4) => 1, (:Lj, 1, 3) => 1, (:R, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 1, (:C, 2, 3) => 1)
indexdict = Dict((:C, 2, 3, 1) => 3, (:Lj, 1, 3, 1) => 4, (:C, 1, 3, 1) => 5, (:R, 1, 2, 1) => 2, (:C, 3, 4, 1) => 6, (:P, 1, 2, 1) => 1, (:C, 1, 4, 1) => 8, (:Lj, 1, 4, 1) => 7)
Cj, Icmean = JosephsonCircuits.calcCjIcmean(componenttypes, nodeindexarray,
    componentvalues, componentnamedict,mutualinductorbranchnames, countdict, indexdict)

# output
(3.1100514732000003e-13, 3.1414661345454545e-7)
```
```jldoctest
componenttypes = [:P, :R, :C, :Lj, :C, :C, :Lj, :C]
nodeindexarray = [2 2 2 3 3 3 4 4; 1 1 3 1 1 4 1 1]
componentvalues = Real[1, 50.0, 1.0e-13, 2.0e-9, 1.0e-12, 1.0e-13, 1.1e-9, 1.2e-12]
componentnamedict = Dict("R1" => 2, "Cc2" => 6, "Cj2" => 8, "Cj1" => 5, "P1" => 1, "Cc1" => 3, "Lj2" => 7, "Lj1" => 4)
mutualinductorbranchnames = String[]
countdict = Dict((:Lj, 1, 4) => 1, (:C, 3, 4) => 1, (:C, 1, 4) => 1, (:Lj, 1, 3) => 1, (:R, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 1, (:C, 2, 3) => 1)
indexdict = Dict((:C, 2, 3, 1) => 3, (:Lj, 1, 3, 1) => 4, (:C, 1, 3, 1) => 5, (:R, 1, 2, 1) => 2, (:C, 3, 4, 1) => 6, (:P, 1, 2, 1) => 1, (:C, 1, 4, 1) => 8, (:Lj, 1, 4, 1) => 7)
Cj, Icmean = JosephsonCircuits.calcCjIcmean(componenttypes, nodeindexarray,
    componentvalues, componentnamedict,mutualinductorbranchnames, countdict, indexdict)

# output
(2.2955141825999997e-13, 2.3187011945454544e-7)
```
"""
function calcCjIcmean(componenttypes::Vector{Symbol},
    nodeindexarray::Matrix{Int}, componentvalues::Vector,
    componentnamedict::Dict, mutualinductorbranchnames::Vector{String},
    countdict::Dict, indexdict::Dict)

    # make a copy of these dictionaries so that i don't modify them
    countdictcopy = copy(countdict)
    indexdictcopy = copy(indexdict)

    # first loop to calculate the junction and junction capacitance parameters.
    # in WRSPICE, a JJ needs a capacitor. 
    Icmean = 0
    Icmax = 0
    Icmin = 0
    Cjmean = 0
    CjoIc = 0

    nJJ = 0
    mutualinductorindex = 0
    for i in eachindex(componenttypes)

        if componenttypes[i] == :K
            mutualinductorindex+=1
        end

        node1, node2 = calcnodes(i, mutualinductorindex, componenttypes,
            nodeindexarray, componentnamedict, mutualinductorbranchnames)

        # sum up the values on the branch
        flag, value, index = sumbranchvalues!(componenttypes[i], node1, node2, componentvalues, countdictcopy, indexdictcopy)
        # println(componenttypes[i]," ",flag," ",value)

        if flag == true && componenttypes[i] == :Lj
            nJJ += 1

            capflag, capvalue, capindex = sumbranchvalues!(:C, node1, node2, componentvalues, countdictcopy, indexdictcopy)
            # println(:C," ",capflag," ", capvalue)

            # if !capflag
            #     error("Each Josephson junction needs a capacitor.")
            # end

            Ictmp = real(LjtoIc(value))
            Icmean = Icmean + (Ictmp-Icmean)/nJJ

            CjoIctmp = real(capvalue/Ictmp)

            if nJJ == 1
                CjoIc = CjoIctmp
                Icmin = Ictmp
            end
            
            if Ictmp < Icmin
                Icmin = Ictmp
            end

            if Ictmp > Icmax
                Icmax = Ictmp
            end

            # we can always add a separate junction capacitance so we want to find the minimum
            # Cj / Ic to use in the jj model. 
            if CjoIctmp == 0.0
                error("Cj cannot be zero in the WRSPICE JJ model.")
            elseif CjoIctmp < CjoIc
                CjoIc = CjoIctmp
            end
        end
    end

    # decide on the circuit parameters
    #decide on the JJ parameters and write the junction model.
    # Icmean*CjoIc

    # println(Icmin," ",Icmax," ",Icmean)
    # println(Icmin/Icmean)
    # println(Icmax/Icmean)

    # check if the junction sizes are within the range allowed by WRSPICE
    if Icmin/Icmean < 0.02
        error("Minimum junction too much smaller than average for WRSPICE.")
    end
    if Icmax/Icmean > 50.0
        error("Maximum junction too much larger than average for WRSPICE.")
    end

    # check if the ratio of Cj / Ic is within the range allowed by WRSPICE
    if CjoIc > 0.99e-6
        CjoIc = 0.99e-6
    end

    return CjoIc*Icmean, Icmean
end

"""
    exportnetlist(circuit::Vector,circuitdefs::Dict,port::Int = true,
        jj::Bool = true)

# Examples
```jldoctest
@variables R Cc Lj Cj I
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("C1","1","2",Cc),
    ("Lj1","2","0",Lj),
    ("C2","2","0",Cj)]

circuitdefs = Dict(
    Lj =>1000.0e-12,
    Cc => 100.0e-15,
    Cj => 1000.0e-15,
    R => 50.0)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
R1 1 0 50.0
C1 1 2 100.0f
B1 2 0 3 jjk ics=0.32910597599999997u
C2 2 0 674.18508376f
.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=325.81491624f,force=1,vm=9.9

* SPICE Simulation
R1 1 0 50.0
C1 1 2 100.0f
Lj1 2 0 1000.0000000000001p
C2 2 0 1000.0f
```
```jldoctest
@variables R Cc L1 L2 Cj1 Cj2 I1 V1
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("C1","1","2",Cc),
    ("L1","2","0",L1),
    ("L2","2","0",L2),
    ("C2","2","0",Cj1),
    ("C3","2","0",Cj2),
    ("I1","2","0",I1)]

circuitdefs = Dict(
    L1 =>2000.0e-12,
    L2 =>2000.0e-12,
    Cc => 100.0e-15,
    Cj1 => 500.0e-15,
    Cj2 => 500.0e-15,
    R => 50.0,
    I1 =>0.1)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
R1 1 0 50.0
C1 1 2 100.0f
L1 2 0 1000.0000000000001p
C2 2 0 1000.0f

* SPICE Simulation
R1 1 0 50.0
C1 1 2 100.0f
L1 2 0 1000.0000000000001p
C2 2 0 1000.0f
```
```jldoctest
@variables Rleft L1 K1 L2 C2 C3 Lj1
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1))
push!(circuit,("Lj1","2","0",Lj1))
push!(circuit,("K1","L1","L2",K1))
push!(circuit,("L2","2","0",L2))
push!(circuit,("C2","2","0",C2))
push!(circuit,("C3","2","0",C3))
circuitdefs = Dict(
    Rleft => 50.0,
    L1 => 1000.0e-12,
    Lj1 => 1000.0e-12,
    K1 => 0.1,
    L2 => 1000.0e-12,
    C2 => 1000.0e-15,
    C3 => 1000.0e-15)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
R1 1 0 50.0
L1 1 0 1000.0000000000001p
B1 2 0 3 jjk ics=0.32910597599999997u
C2 2 0 1674.18508376f
K1 L1 L2 0.1
L2 2 0 1000.0000000000001p
.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=325.81491624f,force=1,vm=9.9

* SPICE Simulation
R1 1 0 50.0
L1 1 0 1000.0000000000001p
Lj1 2 0 1000.0000000000001p
K1 L1 L2 0.1
L2 2 0 1000.0000000000001p
C2 2 0 2000.0f
```
```jldoctest
@variables Rleft L1 K1 L2 C2 C3 Lj1
circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
push!(circuit,("P1","1","0",1))
push!(circuit,("R1","1","0",Rleft))
push!(circuit,("L1","1","0",L1))
push!(circuit,("Lj1","2","0",Lj1))
push!(circuit,("K1","L2","L1",K1))
push!(circuit,("L2","2","0",L2))
push!(circuit,("C2","2","0",C2))
push!(circuit,("C3","2","0",C3))
circuitdefs = Dict(
    Rleft => 50.0,
    L1 => 1000.0e-12,
    Lj1 => 1000.0e-12,
    K1 => 0.1,
    L2 => 1000.0e-12,
    C2 => 1000.0e-15,
    C3 => 1000.0e-15)

println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = true).netlist)
println("")
println(JosephsonCircuits.exportnetlist(circuit, circuitdefs;port = 1, jj = false).netlist)

# output
* SPICE Simulation
R1 1 0 50.0
L1 1 0 1000.0000000000001p
B1 2 0 3 jjk ics=0.32910597599999997u
C2 2 0 1674.18508376f
K1 L2 L1 0.1
L2 2 0 1000.0000000000001p
.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=325.81491624f,force=1,vm=9.9

* SPICE Simulation
R1 1 0 50.0
L1 1 0 1000.0000000000001p
Lj1 2 0 1000.0000000000001p
K1 L2 L1 0.1
L2 2 0 1000.0000000000001p
C2 2 0 2000.0f
```
"""
function exportnetlist(circuit::Vector,circuitdefs::Dict;port::Int = 1,
        jj::Bool = true)

    # set these to 1 for now, but i should consider how or whether to handle
    # multi-port devices.
    portnodes = 1
    portcurrent = 1

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting=:number)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)
    
    # convert as many values as we can to numerical values using definitions
    # from circuitdefs
    componentvalues = componentvaluestonumber(psc.componentvalues,circuitdefs)

    countdict, indexdict = componentdictionaries(
        psc.componenttypes,
        psc.nodeindices,
        psc.componentnamedict,
        psc.mutualinductorbranchnames,
        )

    Nnodes = length(psc.nodenames)
    componenttypes = psc.componenttypes
    componentnames = psc.componentnames
    nodeindexarray = psc.nodeindices
    uniquenodevector = psc.nodenames
    mutualinductorbranchnames = psc.mutualinductorbranchnames
    componentnamedict = psc.componentnamedict

    # calculate the junction properties
    Cj, Icmean = calcCjIcmean(componenttypes, nodeindexarray, componentvalues, componentnamedict,
        mutualinductorbranchnames, countdict, indexdict)

    CjoIc = Cj/Icmean

    # define scale factors for prefixes
    # multiply by these scale factors
    femto = 1e15
    pico = 1e12
    nano = 1e9
    micro = 1e6
    giga = 1e-9

    # Set vm, (reference icrit)*rsub, which determines the junction resistance
    # default is 16.5e-3 which is extremely lossy. The allowed range is 8e-3 to
    # 100e-3. Once the force flag is enabled we can increase beyond this limit.
    # To turn off the force flag remove force=1 from the jj model argument
    # setting that force=0 does nothing.
    # http://www.wrcad.com/ftp/pub/jj.va
    vm = 99e-1

    # define an array of strings for the netlist
    netlist =  ["* SPICE Simulation"]

    # write the netlist
    # make a copy of the dictionaries so we don't modify the originals
    # not strictly necessary since we don't use them again after the loop below.
    countdictcopy = copy(countdict)
    indexdictcopy = copy(indexdict)
    nJJ = 0
    mutualinductorindex = 0
    for i in eachindex(componenttypes)

        if componenttypes[i] == :K
            mutualinductorindex+=1
        end

        node1, node2 = calcnodes(i, mutualinductorindex, componenttypes,
            nodeindexarray, componentnamedict, mutualinductorbranchnames)

        # sum up the values on the branch
        flag, value, index = sumbranchvalues!(componenttypes[i], node1, node2, componentvalues, countdictcopy, indexdictcopy)

        if flag == true && componenttypes[i] == :Lj

            Ictmp = real(LjtoIc(value))

            # if jj == true, then write the JJ otherwise write and inductor
            if jj == true
                nJJ += 1
                # push!(netlist,"B$(nJJ) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(Nnodes+nJJ-1) jjk ics=$(real(LjtoIc(value)*micro))u")
                push!(netlist,"B$(componentnames[i][3:end]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(Nnodes+nJJ-1) jjk ics=$(real(LjtoIc(value)*micro))u")
                capflag, capvalue, capindex = sumbranchvalues!(:C, node1, node2, componentvalues, countdictcopy, indexdictcopy)

                # add any additional capacitance
                if real(capvalue) > Ictmp*CjoIc
                    push!(netlist,"$(componentnames[capindex]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(femto*real(capvalue-Ictmp*CjoIc))f")
                end
            else
                push!(netlist,"$(componentnames[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value*pico))p")
            end
        elseif flag == true && componenttypes[i] == :L
            push!(netlist,"$(componentnames[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value*pico))p")
        elseif flag == true && componenttypes[i] == :C
            push!(netlist,"$(componentnames[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value*femto))f")
        elseif flag == true && componenttypes[i] == :K
            push!(netlist,"$(componentnames[i]) $(mutualinductorbranchnames[2*mutualinductorindex-1]) $(mutualinductorbranchnames[2*mutualinductorindex]) $(real(value))")
        elseif flag == true && componenttypes[i] == :R
            push!(netlist,"$(componentnames[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value))")
        end
    end

    # # find the nodes for the selected port
    # # i should also support multiple ports. maybe pass in an array of port indices.
    # portnodes=findall(x->x == port, cdict[:P])
    # portcurrent = 0.0
    # if isempty(portnodes)
    #     error("Port $(port) does not exist in dictionary.")
    # else
    #     portnodes=first(portnodes)
    #     # portcurrent=cdict[:I][portnodes]
    # end

    if jj == true && nJJ > 0
        push!(netlist,".model jjk jj(rtype=0,cct=1,icrit=$(micro*Icmean)u,cap=$(femto*Icmean*real(CjoIc))f,force=1,vm=$(vm)")
    end

    return  (netlist=join(netlist,"\n"),portnodes=portnodes,port=port,portcurrent=portcurrent,Nnodes = Nnodes)
end
