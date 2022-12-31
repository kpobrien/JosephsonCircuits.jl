
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
            # if the nodes are changed
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
    consolidatecomponents()

Calculate the component dictionaries where the key is the edge and the value
is the value of the component. 

This function should be retired. Check if it's used in the spice export. 

"""
function consolidatecomponents(typevector,nodeindexarraysorted,
    mutualinductorvector,valuevector,uniquenodevector)

    Nnodes = length(uniquenodevector)

    # nodearraysorted,nodearraysortedw,valuevector)
    componenttypes =unique(typevector)

    # dictionary to hold dictionaries of components
    cdict = Dict()

    # dictionaries to hold arrays of edges for each component type
    edgearray = Dict()

    for componenttype in componenttypes
        # make a different one for mutual inductors, K
        if componenttype == :K
            error("Error: Mutual inductors, :K, not handled yet")
        else
            if componenttype ==:P
                tmp = Dict{Tuple{Int,Int},Int}()
            else
                tmp = Dict{Tuple{Int,Int},Complex{Float64}}()
            end

            for (i,j) in enumerate(findall(x->x==componenttype,typevector))
                key= (nodeindexarraysorted[1,j],nodeindexarraysorted[2,j])
                if haskey(tmp,key)
                    if componenttype == :C
                        # add capacitor in parallel
                        tmp[key]+=valuevector[j]
                    elseif componenttype == :Lj
                        error("No support for multiple JJs on same branch")
                    elseif componenttype == :R || componenttype == :L
                        # add inductor or resistor in parallel
                        tmp[key]=1/(1/valuevector[j]+1/tmp[key])
                    else
                        error("Only one of these per branch.")
                    end
                else
                    tmp[key]=valuevector[j]
                end
            end
            cdict[componenttype]=tmp
        end
    end


    for componenttype in componenttypes
        edgearray[componenttype] = collect(keys(cdict[componenttype]))
    end


    # set some values to zero
    Lmean = 0.0
    ninductors = 0

    # calculate the mean inductance
    for (i,j) in enumerate(findall(x->(x==:Lj) || (x==:L),typevector))
        ninductors+=1
        Lmean = Lmean + (valuevector[j]-Lmean)/ninductors
    end
    if ninductors == 0
        Lmean = 1.0+0.0im
    end

    # return (cdict=cdict,cdictw=cdictw,edgearray=edgearray,Lmean=Lmean,Nnodes=Nnodes)
    return (cdict=cdict,edgearray=edgearray,Lmean=Lmean,Nnodes=Nnodes)
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
B1 2 0 3 jjk ics=0.32910597599999997u
Cj1 2 0 670.8940240000001f
.model jjk jj(rtype=0,cct=1,icrit=0.32910597599999997u,cap=329.105976f,force=1,vm=9.9
C1 1 2 100.0f
R1 1 0 50.0

* SPICE Simulation
Lj1 2 0 1000.0000000000001p
Cj1 2 0 1000.0f
C1 1 2 100.0f
R1 1 0 50.0
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
L1 2 0 1000.0000000000001p
C1 1 2 100.0f
C2 2 0 1000.0f
R1 1 0 50.0

* SPICE Simulation
L1 2 0 1000.0000000000001p
C1 1 2 100.0f
C2 2 0 1000.0f
R1 1 0 50.0
```
"""
function exportnetlist(circuit::Vector,circuitdefs::Dict;port::Int = 1,
        jj::Bool = true)

    # parse and sort the circuit
    psc = parsesortcircuit(circuit, sorting=:number)

    # calculate the circuit graph
    cg = calccircuitgraph(psc)
    
    # convert as many values as we can to numerical values using definitions
    # from circuitdefs
    vvn = valuevectortonumber(psc.valuevector,circuitdefs)

    c = consolidatecomponents(psc.typevector,psc.nodeindexarraysorted,
        psc.mutualinductorvector,vvn,psc.uniquenodevectorsorted)

    Nnodes = length(psc.uniquenodevectorsorted)
    cdict = c.cdict
    typearray = psc.typevector
    namearray = psc.namevector
    nodearray = psc.nodeindexarraysorted
    uniquenodearray = psc.uniquenodevectorsorted
    mutualinductorbranches = psc.mutualinductorvector


    # find the nodes for the selected port
    # i should also support multiple ports. maybe pass in an array of port indices.
    portnodes=findall(x->x == port, cdict[:P])
    portcurrent = 0.0
    if isempty(portnodes)
        error("Port $(port) does not exist in dictionary.")
    else
        portnodes=first(portnodes)
        # portcurrent=cdict[:I][portnodes]
    end

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


    Icmean = 0
    Icmax = 0
    Icmin = 0
    Cjmean = 0

    nJJ=0
    if haskey(cdict,:C)
        Cdict = copy(cdict[:C])
    end
    CjoIc = 0

    inductorlabels = Dict()

    if haskey(cdict,:Lj)

        for (Ljedge,val) in cdict[:Lj]
            nJJ+=1
            Ictmp = LjtoIc(real(cdict[:Lj][Ljedge]))
            Icmean = Icmean + (Ictmp-Icmean)/nJJ
            if !haskey(Cdict,Ljedge)
                error("Error: No junction capacitance found for junction at: $(Ljedge)")
            end

            CjoIctmp = real(pop!(Cdict,Ljedge)/Ictmp)
            # println(nJJ)

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


            # if abs((CjoIctmp - CjoIc)/CjoIc) > 1e-8
            #     error("Cj over critical current ratio not constant. This is required for the WRSPICE JJ model. ")
            # # idea, i can find the minimum CjoIc ratio, set that equal to the junction capacitance, then
            # # add additional capacitance. 
            # else
            #     CjoIc = CjoIctmp
            # end
        end
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
        # if CjoIc > 1e-9
        #     CjoIc = 1e-9
        # end


        ## write the JJs (or linear inductors in place of JJs) ##
        i = 1
        for (key,val) in sort(collect(cdict[:Lj]), by=x->x[1][1])
            Ictmp = LjtoIc(real(val))
            if jj == true
                push!(netlist,"B$(i) $(key[1]-1) $(key[2]-1) $(Nnodes+i-1) jjk ics=$(micro*Ictmp)u")
                # add any additional capacitance
                if abs(real(cdict[:C][key]) - Ictmp*CjoIc) > 1e-18
                    push!(netlist,"Cj$(i) $(key[1]-1) $(key[2]-1) $(femto*real(cdict[:C][key]-Ictmp*CjoIc))f")
                end
            else
                push!(netlist,"Lj$(i) $(key[1]-1) $(key[2]-1) $(pico*real(val))p")
                inductorlabels[(key[1],key[2])] = "Lj$(i)"
                inductorlabels[(key[2],key[1])] = "Lj$(i)"
                push!(netlist,"Cj$(i) $(key[1]-1) $(key[2]-1) $(femto*real(cdict[:C][key]))f")
            end
            i+=1
        end

        if jj == true
            push!(netlist,".model jjk jj(rtype=0,cct=1,icrit=$(micro*Icmean)u,cap=$(femto*Icmean*real(CjoIc))f,force=1,vm=$(vm)")
        end
    end

    if haskey(cdict,:L)
        ## write the linear inductors ##
        i = 1
        for (key,val) in sort(collect(cdict[:L]), by=x->x[1][1])
            push!(netlist,"L$(i) $(key[1]-1) $(key[2]-1) $(pico*real(val))p")
            inductorlabels[(key[1],key[2])] = "L$(i)"
            inductorlabels[(key[2],key[1])] = "L$(i)"
            i+=1
        end
    end

    # if haskey(cdict,:Lm)
    #     ## write the mutual inductors ##
    #     # loop over the mutual inductors 
    #     i = 1
    #     for (key,val) in sort(collect(cdict[:Lm]), by=x->x[1][1])
    #         L1label = inductorlabels[(key[1],key[2])]
    #         L2label = inductorlabels[(key[3],key[4])]
    #         L1 = cdict[:L][(key[1],key[2])]
    #         L2 = cdict[:L][(key[3],key[4])]

    #         K = real(val/sqrt(L1*L2))
    #         push!(netlist,"K$(i) $(L1label) $(L2label) $(K)")
    #         i+=1
    #     end
    # end

    if haskey(cdict,:C)
        ## write the capacitors ##
        i = 1
        for (key,val) in sort(collect(Cdict), by=x->x[1][1])
            push!(netlist,"C$(i) $(key[1]-1) $(key[2]-1) $(femto*real(val))f")
            i+=1
        end
    end

    ## write the resistors ##
    if haskey(cdict,:R)
        i = 1
        for (key,val) in sort(collect(cdict[:R]), by=x->x[1][1])
            push!(netlist,"R$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
            i+=1
        end
    end

    # if haskey(cdict,:V)
    #     ## write the voltage sources ##
    #     i = 1
    #     for (key,val) in sort(collect(cdict[:V]), by=x->x[1][1])
    #         push!(netlist,"V$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
    #         i+=1
    #     end
    # end

    if haskey(cdict,:I)
        ## write the current sources ##
        for (key,val) in sort(collect(cdict[:I]), by=x->x[1][1])
            #push!(netlist,"I$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
            i+=1
        end
    end


    # end


    # nK = 1
    # nJJ = 1
    # for (i,type) in enumerate(typearray)
    #     if type == :K
    #         push!(netlist,"$(namearray[i]) $(mutualinductorbranches[2*nK-1]) $(c0.mutualinductorbranches[2*nK]) $(real(numberarray[i]))")
    #         nK+=1
    #     elseif type == :P
    #         nothing
    #     elseif type == :I
    #         nothing
    #     elseif type == :V
    #         nothing
    #     elseif type == :L
    #         push!(netlist,"$(namearray[i]) $(uniquenodearray[nodearray[2*i-1]]) $(uniquenodearray[nodearray[2*i]]) $(real(numberarray[i]*pico))p")
    #     elseif type == :C
    #         push!(netlist,"$(namearray[i]) $(uniquenodearray[nodearray[2*i-1]]) $(uniquenodearray[nodearray[2*i]]) $(real(numberarray[i]*femto))f")
    #     elseif type == :Lj
    #         # push!(netlist,"$(namearray[i]) $(uniquenodearray[nodearray[2*i-1]]) $(uniquenodearray[nodearray[2*i]]) $(real(numberarray[i]*pico))p")

    #         # Ictmp = LjtoIc(real(val))
    #         if jj == true
    #             push!(netlist,"B$(nJJ) $(uniquenodearray[nodearray[2*i-1]]) $(uniquenodearray[nodearray[2*i]]) $(Nnodes+nJJ-1) jjk ics=$(real(LjtoIc(numberarray[i])*micro))u")
    #             nJJ += 1
    #             # # add any additional capacitance
    #             # if real(cdict[:C][key]) > Ictmp*CjoIc
    #             #     push!(netlist,"Cj$(i) $(key[1]-1) $(key[2]-1) $(femto*real(cdict[:C][key]-Ictmp*CjoIc))f")
    #             # end       
    #         else
    #             # push!(netlist,"Lj$(i) $(key[1]-1) $(key[2]-1) $(pico*real(val))p")
    #             # inductorlabels[(key[1],key[2])] = "Lj$(i)"
    #             # inductorlabels[(key[2],key[1])] = "Lj$(i)"
    #             # push!(netlist,"Cj$(i) $(key[1]-1) $(key[2]-1) $(femto*real(cdict[:C][key]))f")
    #             push!(netlist,"$(namearray[i]) $(uniquenodearray[nodearray[2*i-1]]) $(uniquenodearray[nodearray[2*i]]) $(real(numberarray[i]*pico))p")
    #         end


    #     else            
    #         push!(netlist,"$(namearray[i]) $(uniquenodearray[nodearray[2*i-1]]) $(uniquenodearray[nodearray[2*i]]) $(real(numberarray[i]))")
    #     end
    # end



    # if jj == true && haskey(cdict,:Lj)
    #     push!(netlist,".model jjk jj(rtype=0,cct=1,icrit=$(micro*Icmean)u,cap=$(femto*Icmean*real(CjoIc))f,force=1,vm=$(vm)")
    # end

    return  (netlist=join(netlist,"\n"),portnodes=portnodes,port=port,portcurrent=portcurrent,Nnodes = Nnodes)
end

# function exportnetlist(cdict::Dict{Any, Any},Nnodes::Int,port::Int = 1,
#     jj::Bool = true)

#     # find the nodes for the selected port
#     # i should also support multiple ports. maybe pass in an array of port indices.
#     portnodes=findall(x->x == port, cdict[:P])
#     portcurrent = 0.0
#     if isempty(portnodes)
#         error("Port $(port) does not exist in dictionary.")
#     else
#         portnodes=first(portnodes)
#         portcurrent=cdict[:I][portnodes]
#     end

#     # define scale factors for prefixes
#     # multiply by these scale factors
#     femto = 1e15
#     pico = 1e12
#     nano = 1e9
#     micro = 1e6
#     giga = 1e-9

#     # Set vm, (reference icrit)*rsub, which determines the junction resistance
#     # default is 16.5e-3 which is extremely lossy. The allowed range is 8e-3 to 
#     # 100e-3. Once the force flag is enabled we can increase beyond this limit.
#     # To turn off the force flag remove force=1 from the jj model argument 
#     # setting that force=0 does nothing. 
#     # http://www.wrcad.com/ftp/pub/jj.va
#     vm = 99e-1

#     # define an array of strings for the netlist
#     netlist =  ["* SPICE Simulation"]

#     Icmean = 0
#     Icmax = 0
#     Icmin = 0
#     Cjmean = 0

#     nJJ=0
#     if haskey(cdict,:C)
#         Cdict = copy(cdict[:C])
#     end
#     CjoIc = 0

#     if haskey(cdict,:Lj)

#         for (Ljedge,val) in cdict[:Lj]
#             nJJ+=1
#             Ictmp = LjtoIc(real(cdict[:Lj][Ljedge]))
#             Icmean = Icmean + (Ictmp-Icmean)/nJJ
#             if !haskey(Cdict,Ljedge)
#                 error("Error: No junction capacitance found for junction at: $(Ljedge)")
#             end

#             CjoIctmp = real(pop!(Cdict,Ljedge)/Ictmp)
#             # println(nJJ)

#             if nJJ == 1
#                 CjoIc = CjoIctmp
#                 Icmin = Ictmp
#             end
            
#             if Ictmp < Icmin
#                 Icmin = Ictmp
#             end

#             if Ictmp > Icmax
#                 Icmax = Ictmp
#             end

#             # we can always add a separate junction capacitance so we want to find the minimum
#             # Cj / Ic to use in the jj model. 
#             if CjoIctmp == 0.0
#                 error("Cj cannot be zero in the WRSPICE JJ model.")
#             elseif CjoIctmp < CjoIc
#                 CjoIc = CjoIctmp
#             end


#             # if abs((CjoIctmp - CjoIc)/CjoIc) > 1e-8
#             #     error("Cj over critical current ratio not constant. This is required for the WRSPICE JJ model. ")
#             # # idea, i can find the minimum CjoIc ratio, set that equal to the junction capacitance, then
#             # # add additional capacitance. 
#             # else
#             #     CjoIc = CjoIctmp
#             # end
#         end
#         #decide on the JJ parameters and write the junction model.
#         Icmean*CjoIc

#         # check if the junction sizes are within the range allowed by WRSPICE
#         if Icmin/Icmean < 0.02
#             error("Minimum junction too much smaller than average for WRSPICE.")
#         end
#         if Icmax/Icmean > 50.0
#             error("Maximum junction too much larger than average for WRSPICE.")
#         end

#         # check if the ratio of Cj / Ic is within the range allowed by WRSPICE
#         if CjoIc > 1e-6
#             CjoIc = 1e-6
#         end
#         # if CjoIc > 1e-9
#         #     CjoIc = 1e-9
#         # end


#         ## write the JJs (or linear inductors in place of JJs) ##
#         inductorlabels = Dict()
#         i = 1
#         for (key,val) in sort(collect(cdict[:Lj]), by=x->x[1][1])
#             Ictmp = LjtoIc(real(val))
#             if jj == true
#                 push!(netlist,"B$(i) $(key[1]-1) $(key[2]-1) $(Nnodes+i-1) jjk ics=$(micro*Ictmp)u")
#                 # add any additional capacitance
#                 if real(cdict[:C][key]) > Ictmp*CjoIc
#                     push!(netlist,"Cj$(i) $(key[1]-1) $(key[2]-1) $(femto*real(cdict[:C][key]-Ictmp*CjoIc))f")
#                 end       
#             else
#                 push!(netlist,"Lj$(i) $(key[1]-1) $(key[2]-1) $(pico*real(val))p")
#                 inductorlabels[(key[1],key[2])] = "Lj$(i)"
#                 inductorlabels[(key[2],key[1])] = "Lj$(i)"
#                 push!(netlist,"Cj$(i) $(key[1]-1) $(key[2]-1) $(femto*real(cdict[:C][key]))f")
#             end
#             i+=1
#         end

#         if jj == true
#             push!(netlist,".model jjk jj(rtype=0,cct=1,icrit=$(micro*Icmean)u,cap=$(femto*Icmean*real(CjoIc))f,force=1,vm=$(vm)")
#         end
#     end

#     if haskey(cdict,:L)
#         ## write the linear inductors ##
#         i = 1
#         for (key,val) in sort(collect(cdict[:L]), by=x->x[1][1])
#             push!(netlist,"L$(i) $(key[1]-1) $(key[2]-1) $(pico*real(val))p")
#             inductorlabels[(key[1],key[2])] = "L$(i)"
#             inductorlabels[(key[2],key[1])] = "L$(i)"
#             i+=1
#         end
#     end

#     if haskey(cdict,:Lm)
#         ## write the mutual inductors ##
#         # loop over the mutual inductors 
#         i = 1
#         for (key,val) in sort(collect(cdict[:Lm]), by=x->x[1][1])
#             L1label = inductorlabels[(key[1],key[2])]
#             L2label = inductorlabels[(key[3],key[4])]
#             L1 = cdict[:L][(key[1],key[2])]
#             L2 = cdict[:L][(key[3],key[4])]

#             K = real(val/sqrt(L1*L2))
#             push!(netlist,"K$(i) $(L1label) $(L2label) $(K)")
#             i+=1
#         end
#     end

#     if haskey(cdict,:C)
#         ## write the capacitors ##
#         i = 1
#         for (key,val) in sort(collect(Cdict), by=x->x[1][1])
#             push!(netlist,"C$(i) $(key[1]-1) $(key[2]-1) $(femto*real(val))f")
#             i+=1
#         end
#     end

#     ## write the resistors ##
#     if haskey(cdict,:R)
#         i = 1
#         for (key,val) in sort(collect(cdict[:R]), by=x->x[1][1])
#             push!(netlist,"R$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
#             i+=1
#         end
#     end

#     if haskey(cdict,:V)
#         ## write the voltage sources ##
#         i = 1
#         for (key,val) in sort(collect(cdict[:V]), by=x->x[1][1])
#             push!(netlist,"V$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
#             i+=1
#         end
#     end

#     if haskey(cdict,:I)
#         ## write the current sources ##
#         for (key,val) in sort(collect(cdict[:I]), by=x->x[1][1])
#             #push!(netlist,"I$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
#             i+=1
#         end
#     end

#     return (netlist=join(netlist,"\n"),portnodes=portnodes,port=port,portcurrent=portcurrent)
# end
