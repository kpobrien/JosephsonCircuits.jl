
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

julia> JosephsonCircuits.sumvalues(:V, 1.0, 4.0)
ERROR: unknown component type in sumvalues
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
    sumbranchvalues!(type::Symbol, node1::Int, node2::Int,valuevector::Vector,
        countdict, indexdict)

Given a branch and a type, return the sum of all of the values of the same
type and branch. The sum will behave differently depending on the type.

# Examples
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

psc = JosephsonCircuits.parsesortcircuit(circuit)
vvn = JosephsonCircuits.valuevectortonumber(psc.valuevector,circuitdefs)
countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.typevector,psc.nodeindexarraysorted,psc.mutualinductorvector,psc.namedict)

println(JosephsonCircuits.sumbranchvalues!(:C, 1, 3, vvn, countdict, indexdict))
println(countdict)
println(indexdict)

# output
(true, 1.0e-12, 6)
Dict((:L, 1, 3) => 2, (:R, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 0, (:C, 2, 3) => 1, (:I, 1, 3) => 1)
Dict((:C, 2, 3, 1) => 3, (:C, 1, 3, 1) => 6, (:R, 1, 2, 1) => 2, (:L, 1, 3, 1) => 4, (:C, 1, 3, 2) => 7, (:L, 1, 3, 2) => 5, (:P, 1, 2, 1) => 1, (:I, 1, 3, 1) => 8)
```
"""
function sumbranchvalues!(type::Symbol, node1::Int, node2::Int,
    valuevector::Vector, countdict, indexdict)
    countkey = (type, node1, node2)
    countflag = false
    value = zero(eltype(valuevector))
    index = 0

    if haskey(countdict,countkey)
        counts = countdict[countkey]
        if counts > 0
            countflag = true
            index = indexdict[(type, node1, node2, 1)]
            value = valuevector[index]

            for count in 2:counts
                index1 = indexdict[(type, node1, node2, count)]
                value = sumvalues(type, value, valuevector[index1])
            end
            countdict[countkey] = 0
        end
    end

    return countflag, value, index

end

"""
    componentdictionaries(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int},
        mutualinductorvector::Vector, namedict::Dict)

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
countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.typevector,psc.nodeindexarraysorted,psc.mutualinductorvector,psc.namedict)

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
countdict, indexdict = JosephsonCircuits.componentdictionaries(psc.typevector,psc.nodeindexarraysorted,psc.mutualinductorvector,psc.namedict)

println(countdict)
println(indexdict)

# output
Dict((:L, 1, 3) => 1, (:K, 4, 7) => 2, (:R, 1, 2) => 1, (:I, 1, 2) => 1, (:P, 1, 2) => 1, (:C, 1, 3) => 2, (:L, 1, 2) => 1)
Dict((:C, 1, 3, 1) => 8, (:I, 1, 2, 1) => 2, (:R, 1, 2, 1) => 3, (:K, 4, 7, 1) => 5, (:K, 4, 7, 2) => 6, (:L, 1, 2, 1) => 4, (:L, 1, 3, 1) => 7, (:P, 1, 2, 1) => 1, (:C, 1, 3, 2) => 9)
```
"""
function componentdictionaries(typevector::Vector{Symbol}, nodeindexarray::Matrix{Int},
    mutualinductorvector::Vector, namedict::Dict)

    if  length(typevector) != size(nodeindexarray,2)
        throw(DimensionMismatch("Input arrays must have the same length"))
    end

    if length(size(nodeindexarray)) != 2
        throw(DimensionMismatch("The nodeindexarray must have two dimensions"))
    end

    if size(nodeindexarray,1) != 2
        throw(DimensionMismatch("The length of the first axis must be 2"))
    end

    # key = (componenttype,node1,node2), value = counts
    # countdict = Dict{Tuple{eltype(typevector),eltype(nodeindexarray),eltype(nodeindexarray)},Int}()
    countdict = Dict{Tuple{eltype(typevector),eltype(nodeindexarray),eltype(nodeindexarray)},Int}()
    sizehint!(countdict,length(typevector))

    # key = (node1,node2,count), value = index in typevector
    # indexdict = Dict{Tuple{eltype(nodeindexarray),eltype(nodeindexarray),Int},Int}()
    indexdict = Dict{Tuple{eltype(typevector),eltype(nodeindexarray),eltype(nodeindexarray),Int},Int}()
    sizehint!(indexdict,length(typevector))

    n = 1
    for i in eachindex(typevector)
        # i think i need to take the indices from mutualinductorvector
        if typevector[i] == :K
            # don't sort these because the mutual inductance changes sign
            # if the nodes are changed. this is OK because only values with
            # the same inductor ordering will be summed. 
            # names of inductors
            inductor1name = mutualinductorvector[2*n-1]
            inductor2name = mutualinductorvector[2*n]

            # indices of inductors
            node1 = namedict[inductor1name]
            node2 = namedict[inductor2name]

            countkey = (typevector[i], node1, node2)

            if haskey(countdict,countkey)
                countdict[countkey] += 1
            else
                countdict[countkey] = 1
            end

            indexkey = (typevector[i], node1, node2, countdict[countkey])
            indexdict[indexkey] = i
            n+=1

        else

            node1 = nodeindexarray[1,i]
            node2 = nodeindexarray[2,i]
            if node1 < node2
                countkey = (typevector[i], node1, node2)
            else
                countkey = (typevector[i], node2, node1)
            end

            if haskey(countdict,countkey)
                countdict[countkey] += 1
            else
                countdict[countkey] = 1
            end

            if node1 < node2
                indexkey = (typevector[i], node1, node2, countdict[countkey])
            else
                indexkey = (typevector[i], node2, node1, countdict[countkey])
            end

            indexdict[indexkey] = i
        end
    end

    return countdict, indexdict
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
C2 2 0 670.8940240000001f
.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9

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
C2 2 0 1670.894024f
K1 L1 L2 0.1
L2 2 0 1000.0000000000001p
.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9

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
C2 2 0 1670.894024f
K1 L2 L1 0.1
L2 2 0 1000.0000000000001p
.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9

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


    portnodes = 1
    portcurrent = 1

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting=:number)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)
    
    # convert as many values as we can to numerical values using definitions
    # from circuitdefs
    valuevector = valuevectortonumber(psc.valuevector,circuitdefs)

    countdict, indexdict = componentdictionaries(
        psc.typevector,
        psc.nodeindexarraysorted,
        psc.mutualinductorvector,
        psc.namedict)

    countdictcopy = copy(countdict)
    indexdictcopy = copy(indexdict)


    # c = consolidatecomponents(psc.typevector,psc.nodeindexarraysorted,
    #     psc.mutualinductorvector,vvn,psc.uniquenodevectorsorted)

    Nnodes = length(psc.uniquenodevectorsorted)
    typevector = psc.typevector
    namevector = psc.namevector
    nodeindexarray = psc.nodeindexarraysorted
    uniquenodevector = psc.uniquenodevectorsorted
    mutualinductorvector = psc.mutualinductorvector
    namedict = psc.namedict

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


    # i think i need to run through the loop multiple times
    # some things need to be pre-calculated. 
    # some questions:
    # if jj=true, are there capacitors for every JJ? what value of specific
    # capacitance should i use for the jjs and what capacitance should i add
    # across the jj?

    # if jj = true, i can loop over typearray, and whenever i find a jj (check
    # jj count too and sum up) , i can check the count of capacitors for that
    # node (with both node orders) to see if there is a capacitor. i can also
    # record the value of that capacitance and inductance to determine the
    # junction capacitance. after i find a copy, i need to set the counts to
    # zero so i don't count it again. i should work on a copy of the dictionary
    # so i can use it again.
    # i should also have a consolidatecomponents argument. if jj is true, i have
    # to do it for capacitors on the jj. i should have an option to do it for
    # other components also. or should i? i think so, it seems to help with 
    # efficiency in wrspice. 
    # oh yeah, i almost forgot, i have to deal with the stupid phase nodes.
    # fortunately the code i already have seems to work well for that. 

    Icmean = 0
    Icmax = 0
    Icmin = 0
    Cjmean = 0
    CjoIc = 0


    # first loop to calculate the junction parameters
    nJJ = 0
    nK = 0
    for i in eachindex(typevector)

        # calculate the nodes
        if typevector[i] == :K
            nK+=1
            # don't sort these because the mutual inductance changes sign
            # if the nodes are changed. this is OK because only values with
            # the same inductor ordering will be summed. 

            # names of inductors
            inductor1name = mutualinductorvector[2*nK-1]
            inductor2name = mutualinductorvector[2*nK]

            # indices of inductors
            node1 = namedict[inductor1name]
            node2 = namedict[inductor2name]

        else
            if nodeindexarray[1,i] < nodeindexarray[2,i]
                node1 = nodeindexarray[1,i]
                node2 = nodeindexarray[2,i]
            else
                node1 = nodeindexarray[2,i]
                node2 = nodeindexarray[1,i]
            end
        end

        # sum up the values on the branch
        flag, value, index = sumbranchvalues!(typevector[i], node1, node2, valuevector, countdict, indexdict)
        # println(typevector[i]," ",flag," ",value)

        if flag == true && typevector[i] == :Lj
            nJJ += 1

            capflag, capvalue, capindex = sumbranchvalues!(:C, node1, node2, valuevector, countdict, indexdict)
            # println(:C," ",capflag," ", capvalue)

            if !capflag
                error("Each Josephson junction needs a capacitor.")
            end

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
    Icmean*CjoIc

    # check if the junction sizes are within the range allowed by WRSPICE
    if Icmin/Icmean < 0.02
        error("Minimum junction too much smaller than average for WRSPICE.")
    end
    if Icmax/Icmean > 50.0
        error("Maximum junction too much larger than average for WRSPICE.")
    end

    # check if the ratio of Cj / Ic is within the range allowed by WRSPICE
    if CjoIc > 1e-6
        CjoIc = 1e-6
    end


    # make a second loop to write the netlist
    nJJ = 0
    nK = 0
    for i in eachindex(typevector)

        # calculate the nodes
        if typevector[i] == :K
            nK+=1
            # don't sort these because the mutual inductance changes sign
            # if the nodes are changed. this is OK because only values with
            # the same inductor ordering will be summed. 

            # names of inductors
            inductor1name = mutualinductorvector[2*nK-1]
            inductor2name = mutualinductorvector[2*nK]

            # indices of inductors
            node1 = namedict[inductor1name]
            node2 = namedict[inductor2name]

        else
            if nodeindexarray[1,i] < nodeindexarray[2,i]
                node1 = nodeindexarray[1,i]
                node2 = nodeindexarray[2,i]
            else
                node1 = nodeindexarray[2,i]
                node2 = nodeindexarray[1,i]
            end
        end

        # sum up the values on the branch
        flag, value, index = sumbranchvalues!(typevector[i], node1, node2, valuevector, countdictcopy, indexdictcopy)
        # println(typevector[i]," ",flag," ",value)

        if flag == true && typevector[i] == :Lj

            Ictmp = real(LjtoIc(value))
            # CjoIctmp = real(capvalue/Ictmp)

            # if jj == true, then write the JJ otherwise write and inductor
            # Ictmp = LjtoIc(real(val))
            if jj == true
                nJJ += 1
                # push!(netlist,"B$(nJJ) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(Nnodes+nJJ-1) jjk ics=$(real(LjtoIc(value)*micro))u")
                push!(netlist,"B$(namevector[i][3:end]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(Nnodes+nJJ-1) jjk ics=$(real(LjtoIc(value)*micro))u")
                capflag, capvalue, capindex = sumbranchvalues!(:C, node1, node2, valuevector, countdictcopy, indexdictcopy)

                # add any additional capacitance
                if real(capvalue) > Ictmp*CjoIc
                    push!(netlist,"$(namevector[capindex]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(femto*real(capvalue-Ictmp*CjoIc))f")
                end
            else
                push!(netlist,"$(namevector[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value*pico))p")
            end
        elseif flag == true && typevector[i] == :L
            push!(netlist,"$(namevector[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value*pico))p")
        elseif flag == true && typevector[i] == :C
            push!(netlist,"$(namevector[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value*femto))f")
        elseif flag == true && typevector[i] == :K
            push!(netlist,"$(namevector[i]) $(mutualinductorvector[2*nK-1]) $(mutualinductorvector[2*nK]) $(real(value))")
        elseif flag == true && typevector[i] == :R
            push!(netlist,"$(namevector[i]) $(uniquenodevector[nodeindexarray[1, i]]) $(uniquenodevector[nodeindexarray[2, i]]) $(real(value))")
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
