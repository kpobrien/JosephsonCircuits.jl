

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
                tmp = Dict{Tuple{Int64,Int64},Int64}()       
            else
                tmp = Dict{Tuple{Int64,Int64},Complex{Float64}}()
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
                        tmp[key]=1/(1/valuevector[j]+1/cdict[componenttype][key])
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
    exportnetlist(circuit::Vector,circuitdefs::Dict,port::Int64 = true,
        jj::Bool = true)

"""
function exportnetlist(circuit::Vector,circuitdefs::Dict;port::Int64 = 1,
        jj::Bool = true)

    # parse the circuit components
    c0 = parsecircuit(circuit)

    numberarray = valuearraytonumber(c0.valuearray,circuitdefs)

    # sort the nodes
    nodearraysorted = sortnodes(c0.uniquenodearray,c0.nodearray,sorting=:number)

    edgearray = extractedges(c0.typearray,nodearraysorted)

    g = calcgraphs(edgearray,length(c0.uniquenodearray))
    
    c = consolidatecomponents(c0.typearray,nodearraysorted,
        c0.mutualinductorbranches,numberarray,c0.uniquenodearray)

    Nnodes = length(c0.uniquenodearray)
    cdict = c.cdict
    typearray = c0.typearray
    namearray = c0.namearray
    nodearray = c0.nodearray
    uniquenodearray = c0.uniquenodearray
    mutualinductorbranches = c0.mutualinductorbranches


    # find the nodes for the selected port
    # i should also support multiple ports. maybe pass in an array of port indices.
    portnodes=findall(x->x == port, cdict[:P])
    portcurrent = 0.0
    if isempty(portnodes)
        error("Port $(port) does not exist in dictionary.")
    else
        portnodes=first(portnodes)
        portcurrent=cdict[:I][portnodes]
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
        inductorlabels = Dict()
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

    if haskey(cdict,:Lm)
        ## write the mutual inductors ##
        # loop over the mutual inductors 
        i = 1
        for (key,val) in sort(collect(cdict[:Lm]), by=x->x[1][1])
            L1label = inductorlabels[(key[1],key[2])]
            L2label = inductorlabels[(key[3],key[4])]
            L1 = cdict[:L][(key[1],key[2])]
            L2 = cdict[:L][(key[3],key[4])]

            K = real(val/sqrt(L1*L2))
            push!(netlist,"K$(i) $(L1label) $(L2label) $(K)")
            i+=1
        end
    end

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

    if haskey(cdict,:V)
        ## write the voltage sources ##
        i = 1
        for (key,val) in sort(collect(cdict[:V]), by=x->x[1][1])
            push!(netlist,"V$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
            i+=1
        end
    end

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

function exportnetlist(cdict::Dict{Any, Any},Nnodes::Int64,port::Int64 = 1,
    jj::Bool = true)

    # find the nodes for the selected port
    # i should also support multiple ports. maybe pass in an array of port indices.
    portnodes=findall(x->x == port, cdict[:P])
    portcurrent = 0.0
    if isempty(portnodes)
        error("Port $(port) does not exist in dictionary.")
    else
        portnodes=first(portnodes)
        portcurrent=cdict[:I][portnodes]
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
        inductorlabels = Dict()
        i = 1
        for (key,val) in sort(collect(cdict[:Lj]), by=x->x[1][1])
            Ictmp = LjtoIc(real(val))
            if jj == true
                push!(netlist,"B$(i) $(key[1]-1) $(key[2]-1) $(Nnodes+i-1) jjk ics=$(micro*Ictmp)u")
                # add any additional capacitance
                if real(cdict[:C][key]) > Ictmp*CjoIc
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

    if haskey(cdict,:Lm)
        ## write the mutual inductors ##
        # loop over the mutual inductors 
        i = 1
        for (key,val) in sort(collect(cdict[:Lm]), by=x->x[1][1])
            L1label = inductorlabels[(key[1],key[2])]
            L2label = inductorlabels[(key[3],key[4])]
            L1 = cdict[:L][(key[1],key[2])]
            L2 = cdict[:L][(key[3],key[4])]

            K = real(val/sqrt(L1*L2))
            push!(netlist,"K$(i) $(L1label) $(L2label) $(K)")
            i+=1
        end
    end

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

    if haskey(cdict,:V)
        ## write the voltage sources ##
        i = 1
        for (key,val) in sort(collect(cdict[:V]), by=x->x[1][1])
            push!(netlist,"V$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
            i+=1
        end
    end

    if haskey(cdict,:I)
        ## write the current sources ##
        for (key,val) in sort(collect(cdict[:I]), by=x->x[1][1])
            #push!(netlist,"I$(i) $(key[1]-1) $(key[2]-1) $(real(val))")
            i+=1
        end
    end

    return (netlist=join(netlist,"\n"),portnodes=portnodes,port=port,portcurrent=portcurrent)
end
