


# """
#   calcphibports(nodevalues,branches,Nmodes)

# Calculate the branch values (variables) at the branches specified in the array 
# of tuples which contain the two nodes that make up the branch. Nmodes gives
# the number of frequency modes.

# # Examples
# ```jldoctest
# julia> JosephsonCircuits.calcbranchvalues([1 2 3;4 5 6],[(2,1)],1)
# 1×3 Matrix{Int64}:
#  1  2  3
# ```
# """
# function calcbranchvalues(nodevalues,branches,Nmodes)

#     # branchvalues = Matrix(eltype(nodevalues),length(branches)*Nmodes,size(nodevalues,2))
#     branchvalues = Matrix{eltype(nodevalues)}(undef,length(branches)*Nmodes,size(nodevalues,2))

#     calcbranchvalues!(branchvalues,nodevalues,branches,Nmodes)

#     return branchvalues

# end


# """
#     calcbranchvalues!(branchvalues,nodevalues,branches,Nmodes)

# Calculate the branch values (variables) at the branches specified in the array 
# of tuples which contain the two nodes that make up the branch. Return the result
# in branchvalues. Nmodes gives the number of frequency modes.

# Note: I think this is using the wrong sign convention.
# In Kerman et al.  
# phib_{branch_{i,j}} = phin_j - phin_i
# In Rasmussen it's end - start.

# I need to find out where else this is swapped because changing the signs here
# messes up the S11 calculation. 

# """
# function calcbranchvalues!(branchvalues,nodevalues,branches,Nmodes)

#     if size(branchvalues,2) != size(nodevalues,2)
#         error("Error: inconsistent second dimension for branchflux and nodeflux.")
#     end

#     if size(branchvalues,1) != length(branches)*Nmodes
#         error("Error: inconsistent dimensions for branchflux and number of ports and modes.")
#     end

#     for (i,key) in enumerate(branches)
#         # calculate the branch flux from the node flux. if there are flux
#         # quanta in the loop this might be incorrect. if any of the keys
#         # are equal to 1, then it is the ground node and that flux will be
#         # zero

#         if key[1] == 1
#             for j in 1:Nmodes
#                 for k in 1:size(nodevalues,2)
#                     branchvalues[(i-1)*Nmodes+j,k] = -nodevalues[(key[2]-2)*Nmodes+j,k]
#                 end
#             end
#         elseif key[2] == 1
#             for j in 1:Nmodes
#                 for k in 1:size(nodevalues,2)
#                     branchvalues[(i-1)*Nmodes+j,k] =  nodevalues[(key[1]-2)*Nmodes+j,k]
#                 end
#             end
#         else
#             for j in 1:Nmodes
#                 for k in 1:size(nodevalues,2)
#                     branchvalues[(i-1)*Nmodes+j,k] =  nodevalues[(key[1]-2)*Nmodes+j,k] 
#                     branchvalues[(i-1)*Nmodes+j,k] -= nodevalues[(key[2]-2)*Nmodes+j,k]
#                 end
#             end
#         end
#     end
#     return nothing
# end

# function calcportcurrent!(input,bnm,resistordict,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     if size(input,1) != length(resistordict)*Nmodes
#         throw(DimensionMismatch("First axis of input matrix must be number of ports times number of modes"))
#     end

#     if size(input,2) != size(bnm,2)
#         throw(DimensionMismatch("Size of second dimension of input matrix must be the same as the size of the second dimension of bnm"))
#     end

#     # loop over the ports
#     for (j,(key,val)) in enumerate(resistordict)
#         # loop over the modes
#         for k = 1:Nmodes
#             # if symbolic, convert to numeric
#             if val isa Symbolic
#                 if !(symfreqvar isa Symbolic)
#                     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
#                 end
#                 resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
#             else
#                 resistance = val
#             end

#             # find the right term in bnm
#             if key[1] == 1
#                 Ip = -bnm[(key[2]-2)*Nmodes+k,:]
#             elseif key[2] == 1
#                 Ip =  bnm[(key[1]-2)*Nmodes+k,:]
#             else
#                 Ip =  bnm[(key[1]-2)*Nmodes+k,:] 
#                 Ip .-= bnm[(key[2]-2)*Nmodes+k,:]
#             end

#             input[(j-1)*Nmodes+k,:] .= 
#                 Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#                 # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))

#             # input[(j-1)*Nmodes+k,(j-1)*Nmodes+k] = 
#             #     Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#             #     # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#         end
#     end

#     return nothing
# end

# # i should get this from the boundary condition matrix.
# function calcinput!(input,Ip,phibports,resistordict,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     # i should have a case for when there is no resistor at the port. 
#     for (j,(key,val)) in enumerate(resistordict)
#         for k = 1:Nmodes

#             # if symbolic, convert to numeric
#             if val isa Symbolic
#                 if !(symfreqvar isa Symbolic)
#                     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
#                 end
#                 resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
#             else
#                 resistance = val
#             end
#             input[(j-1)*Nmodes+k,(j-1)*Nmodes+k] = 
#                 Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#                 # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))

#         end
#     end

#     return nothing
# end


# """

# """
# function calcphibthevenin!(input,bnm,resistordict,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     if size(input,1) != length(resistordict)*Nmodes
#         throw(DimensionMismatch("First axis of input matrix must be number of ports times number of modes"))
#     end

#     if size(input,2) != size(bnm,2)
#         throw(DimensionMismatch("Size of second dimension of input matrix must be the same as the size of the second dimension of bnm"))
#     end

#     # loop over the ports
#     for (j,(key,val)) in enumerate(resistordict)
#         # loop over the modes
#         for k = 1:Nmodes
#             # if symbolic, convert to numeric
#             if val isa Symbolic
#                 if !(symfreqvar isa Symbolic)
#                     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
#                 end
#                 resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
#             else
#                 resistance = val
#             end

#             # should i loop over bnm instead?
#             # then take 1:Nmodes terms for each of those. 
#             # this way i don't need to know where the inputs are. 
#             # 
#             # 

#             # # find the right term in bnm
#             # if key[1] == 1
#             #     Ip = -bnm[(key[2]-2)*Nmodes+k,(j-1)*Nmodes+k]
#             # elseif key[2] == 1
#             #     Ip =  bnm[(key[1]-2)*Nmodes+k,(j-1)*Nmodes+k]
#             # else
#             #     Ip =  bnm[(key[1]-2)*Nmodes+k,(j-1)*Nmodes+k] 
#             #     Ip -= bnm[(key[2]-2)*Nmodes+k,(j-1)*Nmodes+k]
#             # end

#             # find the right term in bnm
#             if key[1] == 1
#                 Ip = -bnm[(key[2]-2)*Nmodes+k,:]
#             elseif key[2] == 1
#                 Ip =  bnm[(key[1]-2)*Nmodes+k,:]
#             else
#                 Ip =  bnm[(key[1]-2)*Nmodes+k,:] 
#                 Ip .-= bnm[(key[2]-2)*Nmodes+k,:]
#             end

#             input[(j-1)*Nmodes+k,:] .= 
#                 Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#                 # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))

#             # input[(j-1)*Nmodes+k,(j-1)*Nmodes+k] = 
#             #     Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#             #     # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#         end
#     end

#     return nothing
# end


# # i should get this from the boundary condition matrix.
# function calcinput!(input,Ip,phibports,portdict,resistordict,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     j=0
#     # i should have a case for when there is no resistor at the port. 
#     for (key,val) in portdict
#         j+=1
#         for k = 1:Nmodes

#             # if symbolic, convert to numeric
#             if resistordict[key] isa Symbolic
#                 if !(symfreqvar isa Symbolic)
#                     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
#                 end
#                 resistance = substitute(resistordict[key],Dict(symfreqvar=>wmodes[k]))
#             else
#                 resistance = resistordict[key]
#             end
#             input[(j-1)*Nmodes+k,(j-1)*Nmodes+k] = 
#                 Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
#                 # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))

#         end
#     end

#     return nothing
# end

# function calcoutput!(output,phibports,resistordict,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     if size(output,2) != size(phibports,2)
#         error("Error: inconsistent dimensions for output and phibports.")
#     end

#     for (j,(key,val)) in enumerate(resistordict)
#         # for k = 1:Nmodes
#         #     output[(j-1)*Nmodes+k,:] .= 
#         #         im*wmodes[k]*phibports[(j-1)*Nmodes+k,:]./sqrt(resistordict[key])/sqrt(abs(wmodes[k]))
#         # end
#         for k = 1:Nmodes
#             # if symbolic, convert to numeric
#             if val isa Symbolic
#                 if !(symfreqvar isa Symbolic)
#                     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
#                 end
#                 resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
#             else
#                 resistance = val
#             end

#             for l = 1:size(output,2)
#                 output[(j-1)*Nmodes+k,l] = 
#                     im*wmodes[k]*phibports[(j-1)*Nmodes+k,l]/sqrt(resistance)/sqrt(abs(wmodes[k]))
#             end
#         end
#     end

#     return nothing
# end

# """


# portimpedance can be for several sets of independent inputs. 
# second dimension is arbitrary
# first dimension is length(resistorvalues)*Nmodes

# """
# function calcportimpedance!(portimpedance,resistorvalues,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     # if size(branchvalues,2) != size(nodevalues,2)
#     #     error("Error: inconsistent second dimension for branchflux and nodeflux.")
#     # end

#     if size(portimpedance,1) != length(resistorvalues)*Nmodes
#         error("Error: inconsistent dimensions for portimpedance and number of resistorvalues and modes.")
#     end


#     for (j,val) in enumerate(resistorvalues)
#         for k = 1:Nmodes
#             # if symbolic, convert to numeric
#             # i might want to put the check for type outside of this loop
#             if val isa Symbolic
#                 if !(symfreqvar isa Symbolic)
#                     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
#                 end
#                 resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
#             else
#                 resistance = val
#             end

#             # take the complex conjugate if frequency is negative
#             if wmodes[k] < 0
#                 resistance = conj(resistance)
#             end

#             for l = 1:size(portimpedance,2)
#                 # current[(j-1)*Nmodes+k,l] *= resistance/sqrt(resistance)/sqrt(abs(wmodes[k]))
#                 portimpedance[(j-1)*Nmodes+k,l] = resistance
#             end
#         end
#     end

#     return nothing
# end


# function scaleby!(portimpedance,wmodes)

#     Nmodes = length(wmodes)
#     Nbranches = div(size(portimpedance,1),Nmodes)

#     if size(portimpedance,1) != Nbranches*Nmodes
#         error("length of first axis of portimpedance must be an integer multiple of Nmodes")
#     end

#     for i = 1:Nbranches
#         for j = 1:Nmodes
#             for k = 1:size(portimpedance,2)
#                 portimpedance[(i-1)*Nmodes+j,k] *= wmodes[j]
#             end
#         end
#     end

#     return nothing
# end


# function calcoutput!(output,phibports,resistorvalues,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     if size(output,2) != size(phibports,2)
#         error("Error: inconsistent dimensions for output and phibports.")
#     end

#     # for (j,(key,val)) in enumerate(resistordict)
#     for (j,val) in enumerate(resistorvalues)
#         # for k = 1:Nmodes
#         #     output[(j-1)*Nmodes+k,:] .= 
#         #         im*wmodes[k]*phibports[(j-1)*Nmodes+k,:]./sqrt(resistordict[key])/sqrt(abs(wmodes[k]))
#         # end
#         for k = 1:Nmodes
#             # if symbolic, convert to numeric
#             if val isa Symbolic
#                 if !(symfreqvar isa Symbolic)
#                     error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
#                 end
#                 resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
#             else
#                 resistance = val
#             end

#             for l = 1:size(output,2)
#                 output[(j-1)*Nmodes+k,l] = 
#                     im*wmodes[k]*phibports[(j-1)*Nmodes+k,l]/sqrt(resistance)/sqrt(abs(wmodes[k]))
#             end
#         end
#     end

#     return nothing
# end


"""
  calcS!

Return the scattering parameters for the system linearized around the strong pump. 
"""


# function calcS!(S,inputwave,outputwave,phin,bnm,inputportbranches,outputportbranches,
#     inputportimpedance,outputportimpedance,wmodes,symfreqvar)
function calcS!(S,inputwave,outputwave,phin,bnm,inputportindices,outputportindices,
    inputportimpedances,outputportimpedances,nodeindexarraysorted,typevector,wmodes,symfreqvar)
    # check the size of S

    # check the size of inputwave

    # check the size of outputwave

    # check the sizes of all of the inputs

    # loop over input branches and modes to define inputwaves
    Ninputports = length(inputportindices)
    Noutputports = length(outputportindices)
    Nsolutions = size(phin,2)
    Nmodes = length(wmodes)

    for i = 1:Ninputports
        for j = 1:Nmodes
            for k = 1:Nsolutions

                # key = inputportbranches[i]
                key = (nodeindexarraysorted[1,inputportindices[i]],nodeindexarraysorted[2,inputportindices[i]])

                # calculate the branch source currents at the ports from the node
                # source current array bnm
                if key[1] == 1
                    sourcecurrent = -bnm[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k]
                else
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k] 
                    sourcecurrent -= bnm[(key[2]-2)*Nmodes+j,k]
                end

                # calculate the branch fluxes at the ports from the node flux array phin
                if key[1] == 1
                    portvoltage = -phin[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k]
                else
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k] 
                    portvoltage -= phin[(key[2]-2)*Nmodes+j,k]
                end

                # scale the branch flux by frequency to get voltage
                portvoltage *= im*wmodes[j]

                # # calculate the port impedance
                # if inputportimpedances[i] isa Symbolic
                #     if !(symfreqvar isa Symbolic)
                #         error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                #     end
                #     portimpedance = substitute(inputportimpedances[i],Dict(symfreqvar=>wmodes[j]))
                # else
                #     portimpedance = inputportimpedances[i]
                # end

                # # portimpedance = inputportimpedances[i]
                # # take the complex conjugate if frequency is negative
                # if wmodes[j] < 0
                #     portimpedance = conj(portimpedance)
                # end
                # println(typevector[inputportindices[i]])
                portimpedance = calcimpedance(inputportimpedances[i],typevector[inputportindices[i]],wmodes[j],symfreqvar)

                # calculate the current flowing through the port
                portcurrent = sourcecurrent - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(real(portimpedance))

                # convert from sqrt(power) to sqrt(photons/second)
                kval *= 1 /sqrt(abs(wmodes[j]))

                # calculate the input and output power waves as defined in (except in
                # units of sqrt(photons/second) instead of sqrt(power)
                # K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE Trans.
                # Micr. Theory and Tech. 13, 194–202 (1965) 
                # doi: 10.1109/TMTT.1965.1125964
                # inputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage + portimpedance * portcurrent)
                # we can simplify the above to:
                inputwave[(i-1)*Nmodes+j,k] = 1/2*kval * portimpedance * sourcecurrent
                # outputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage - conj(portimpedance) * portcurrent)

            end
        end
    end

    # loop over output branches and modes to define outputwaves
    for i = 1:Noutputports
        for j = 1:Nmodes
            for k = 1:Nsolutions

                # key = outputportbranches[i]
                key = (nodeindexarraysorted[1,outputportindices[i]],nodeindexarraysorted[2,outputportindices[i]])

                # calculate the branch source currents at the ports from the node
                # source current array bnm
                if key[1] == 1
                    sourcecurrent = -bnm[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k]
                else
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k] 
                    sourcecurrent -= bnm[(key[2]-2)*Nmodes+j,k]
                end

                # calculate the branch fluxes at the ports from the node flux array phin
                if key[1] == 1
                    portvoltage = -phin[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k]
                else
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k] 
                    portvoltage -= phin[(key[2]-2)*Nmodes+j,k]
                end

                # scale the branch flux by frequency to get voltage
                portvoltage *= im*wmodes[j]

                # # calculate the port impedance
                # if outputportimpedances[i] isa Symbolic
                #     if !(symfreqvar isa Symbolic)
                #         error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                #     end
                #     portimpedance = substitute(outputportimpedances[i],Dict(symfreqvar=>wmodes[j]))
                # else
                #     portimpedance = outputportimpedances[i]
                # end

                # # portimpedance = vvn[outputportindices[i]]
                # # take the complex conjugate if frequency is negative
                # if wmodes[j] < 0
                #     portimpedance = conj(portimpedance)
                # end
                portimpedance = calcimpedance(outputportimpedances[i],typevector[outputportindices[i]],wmodes[j],symfreqvar)


                # calculate the current flowing through the port
                portcurrent = sourcecurrent - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(real(portimpedance))

                # convert from sqrt(power) to sqrt(photons/second)
                kval *= 1 /sqrt(abs(wmodes[j]))

                # calculate the input and output power waves as defined in (except in
                # units of sqrt(photons/second) instead of sqrt(power)
                # K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE Trans.
                # Micr. Theory and Tech. 13, 194–202 (1965) 
                # doi: 10.1109/TMTT.1965.1125964
                # inputwave[(i-1)*Nmodes+j,k] = 1/2*k * (portvoltage + portimpedance * portcurrent)
                # we can simplify the above to:
                # inputwave[(i-1)*Nmodes+j,k] = 1/2*k * portimpedance * sourcecurrent
                outputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage - conj(portimpedance) * portcurrent)
            end
        end
    end
 
    # scattering matrix is defined as outputwave = S * inputwave
    rdiv!(outputwave,inputwave)
    copy!(S,outputwave)

    return nothing
end



# this is a bit of a hack but i ran into issues with complex capacitance when
# the capacitor was at the same branch as a current source. the calcS function
# would use that current source in calculating the output waves, which it should
# not do. 
function calcSnoise!(S,inputwave,outputwave,phin,bnm,inputportindices,outputportindices,
    inputportimpedances,outputportimpedances,nodeindexarraysorted,typevector,wmodes,symfreqvar)
    # check the size of S

    # check the size of inputwave

    # check the size of outputwave

    # check the sizes of all of the inputs

    # loop over input branches and modes to define inputwaves
    Ninputports = length(inputportindices)
    Noutputports = length(outputportindices)
    Nsolutions = size(phin,2)
    Nmodes = length(wmodes)

    for i = 1:Ninputports
        for j = 1:Nmodes
            for k = 1:Nsolutions

                # key = inputportbranches[i]
                key = (nodeindexarraysorted[1,inputportindices[i]],nodeindexarraysorted[2,inputportindices[i]])

                # calculate the branch source currents at the ports from the node
                # source current array bnm
                if key[1] == 1
                    sourcecurrent = -bnm[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k]
                else
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k] 
                    sourcecurrent -= bnm[(key[2]-2)*Nmodes+j,k]
                end

                # calculate the branch fluxes at the ports from the node flux array phin
                if key[1] == 1
                    portvoltage = -phin[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k]
                else
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k] 
                    portvoltage -= phin[(key[2]-2)*Nmodes+j,k]
                end

                # scale the branch flux by frequency to get voltage
                portvoltage *= im*wmodes[j]

                # # calculate the port impedance
                # if inputportimpedances[i] isa Symbolic
                #     if !(symfreqvar isa Symbolic)
                #         error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                #     end
                #     portimpedance = substitute(inputportimpedances[i],Dict(symfreqvar=>wmodes[j]))
                # else
                #     portimpedance = inputportimpedances[i]
                # end

                # # portimpedance = inputportimpedances[i]
                # # take the complex conjugate if frequency is negative
                # if wmodes[j] < 0
                #     portimpedance = conj(portimpedance)
                # end
                # println(typevector[inputportindices[i]])
                portimpedance = calcimpedance(inputportimpedances[i],typevector[inputportindices[i]],wmodes[j],symfreqvar)

                # calculate the current flowing through the port
                portcurrent = sourcecurrent - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(real(portimpedance))

                # convert from sqrt(power) to sqrt(photons/second)
                kval *= 1 /sqrt(abs(wmodes[j]))

                # calculate the input and output power waves as defined in (except in
                # units of sqrt(photons/second) instead of sqrt(power)
                # K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE Trans.
                # Micr. Theory and Tech. 13, 194–202 (1965) 
                # doi: 10.1109/TMTT.1965.1125964
                # inputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage + portimpedance * portcurrent)
                # we can simplify the above to:
                inputwave[(i-1)*Nmodes+j,k] = 1/2*kval * portimpedance * sourcecurrent
                # outputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage - conj(portimpedance) * portcurrent)

            end
        end
    end

    # loop over output branches and modes to define outputwaves
    for i = 1:Noutputports
        for j = 1:Nmodes
            for k = 1:Nsolutions

                # key = outputportbranches[i]
                key = (nodeindexarraysorted[1,outputportindices[i]],nodeindexarraysorted[2,outputportindices[i]])

                # calculate the branch source currents at the ports from the node
                # source current array bnm
                if key[1] == 1
                    sourcecurrent = -bnm[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k]
                else
                    sourcecurrent =  bnm[(key[1]-2)*Nmodes+j,k] 
                    sourcecurrent -= bnm[(key[2]-2)*Nmodes+j,k]
                end

                # calculate the branch fluxes at the ports from the node flux array phin
                if key[1] == 1
                    portvoltage = -phin[(key[2]-2)*Nmodes+j,k]
                elseif key[2] == 1
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k]
                else
                    portvoltage =  phin[(key[1]-2)*Nmodes+j,k] 
                    portvoltage -= phin[(key[2]-2)*Nmodes+j,k]
                end

                # scale the branch flux by frequency to get voltage
                portvoltage *= im*wmodes[j]

                # # calculate the port impedance
                # if outputportimpedances[i] isa Symbolic
                #     if !(symfreqvar isa Symbolic)
                #         error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                #     end
                #     portimpedance = substitute(outputportimpedances[i],Dict(symfreqvar=>wmodes[j]))
                # else
                #     portimpedance = outputportimpedances[i]
                # end

                # # portimpedance = vvn[outputportindices[i]]
                # # take the complex conjugate if frequency is negative
                # if wmodes[j] < 0
                #     portimpedance = conj(portimpedance)
                # end
                portimpedance = calcimpedance(outputportimpedances[i],typevector[outputportindices[i]],wmodes[j],symfreqvar)


                # calculate the current flowing through the port
                # portcurrent = sourcecurrent - portvoltage / portimpedance
                portcurrent = - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(real(portimpedance))

                # convert from sqrt(power) to sqrt(photons/second)
                kval *= 1 /sqrt(abs(wmodes[j]))

                # calculate the input and output power waves as defined in (except in
                # units of sqrt(photons/second) instead of sqrt(power)
                # K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE Trans.
                # Micr. Theory and Tech. 13, 194–202 (1965) 
                # doi: 10.1109/TMTT.1965.1125964
                # inputwave[(i-1)*Nmodes+j,k] = 1/2*k * (portvoltage + portimpedance * portcurrent)
                # we can simplify the above to:
                # inputwave[(i-1)*Nmodes+j,k] = 1/2*k * portimpedance * sourcecurrent
                outputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage - conj(portimpedance) * portcurrent)
            end
        end
    end
 
    # scattering matrix is defined as outputwave = S * inputwave
    rdiv!(outputwave,inputwave)
    copy!(S,outputwave)

    return nothing
end




function calcimpedance(c,type,w,symfreqvar)
    if type == :R
        if w >= 0
            return c
        else
            return conj(c)
        end
    elseif type == :C
        if w >= 0
            return 1/abs(w*imag(c))
            # return -1/(im*w*c)

        else
            return 1/abs(w*imag(c))
            # return 1/imag(w*c)
        end
    else
        error("Unknown component type")
    end
end

function calcimpedance(c::Symbolics.Num,type,w,symfreqvar)
    if type == :R
        if w >= 0
            return Symbolics.substitute(Symbolics.unwrap(c),symfreqvar => w)
        else
            return conj(Symbolics.substitute(Symbolics.unwrap(c),symfreqvar => w))
        end
    elseif type == :C
        if w >= 0
            return 1/(im*w*Symbolics.substitute(Symbolics.unwrap(c),symfreqvar => w))
        else
            return conj(1/(im*w*Symbolics.substitute(Symbolics.unwrap(c),symfreqvar => w)))
        end
    else
        error("Unknown component type")
    end
end

function calcimpedance(c::Symbolics.Symbolic,type,w,symfreqvar)
    if type == :R
        if w >= 0
            return Symbolics.substitute(c,symfreqvar => w)
        else
            return conj(Symbolics.substitute(c,symfreqvar => w))
        end
    elseif type == :C
        if w >= 0
            return 1/(im*w*Symbolics.substitute(c,symfreqvar => w))
        else
            return conj(1/(im*w*Symbolics.substitute(c,symfreqvar => w)))
        end
    else
        error("Unknown component type")
    end
end


"""
    calccm(S,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Sum the abs2 of each element along the horizontal
axis, applying a minus sign if the corresponding frequency is negative. Represents
energy conservation. 

this should equal 3.0e-200
JosephsonCircuits.calccm([1 1e-100 2e-100 1;1 1 1 1],[1 -1])


"""
function calccm(S,w)

    cm = zeros(Float64,size(S,1))

    calccm!(cm,S,w)

    return cm
end

"""
    calccm!(cm,S,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Overwrites cm with output. 

"""
function calccm!(cm,S,w)

    m = length(w)

    if any(mod.(size(S),m) .!=0)
        throw(DimensionMismatch("Dimensions of scattering matrix must be an integer multiple of the number of frequencies"))
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch("First dimension of scattering matrix must equal the length of cm."))        
    end

    # for i = 1:size(S,1)
    #     cm[i] = 0.0
    #     for j = 1:size(S,2)
    #         cm[i] += abs2(S[i,j])*sign(w[(j-1) % m + 1])
    #     end
    # end

    for i = 1:size(S,1)
        c = 0.0
        cm[i] = abs2(S[i,1])*sign(w[(1-1) % m + 1])
        for j = 2:size(S,2)
            t = cm[i] + abs2(S[i,j])*sign(w[(j-1) % m + 1])
            if abs(cm[i]) >= abs2(S[i,j])
                c += (cm[i]-t) + abs2(S[i,j])*sign(w[(j-1) % m + 1])
            else
                c += (abs2(S[i,j])*sign(w[(j-1) % m + 1])-t) + cm[i]
            end
            cm[i] = t
        end
        cm[i]+=c
    end

    return nothing
end

"""
    calccm(S,Snoise,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Sum the abs2 of each element along the horizontal
axis, applying a minus sign if the corresponding frequency is negative. Represents
energy conservation. 

this should equal 6.0e-200

JosephsonCircuits.calccm([1 1e-100 2e-100 1;1 1 1 1],[1 1e-100 2e-100 1;1 1 1 1],[1 -1])


"""
function calccm(S,Snoise,w)

    cm = zeros(Float64,size(S,1))

    calccm!(cm,S,Snoise,w)

    return cm
end

"""
    calccm!(cm,S,Snoise,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Overwrites cm with output. 

"""
function calccm!(cm,S,Snoise,w)

    m = length(w)

    if any(mod.(size(S),m) .!=0)
        throw(DimensionMismatch("Dimensions of scattering matrix must be an integer multiple of the number of frequencies"))
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch("First dimension of scattering matrix must equal the length of cm."))        
    end

    if size(S,1) .!= size(Snoise,1)
        throw(DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
    end


    # for i = 1:size(S,1)
    #     for j = 1:size(S,2)
    #         cm[i] += abs2(S[i,j])*sign(w[(j-1) % m + 1])
    #     end
    #     for j = 1:size(Snoise,2)
    #         cm[i] += abs2(Snoise[i,j])*sign(w[(j-1) % m + 1])
    #     end
    # end

    # use a Kahan, Babushka, Neumaier compensated sum. 
    for i = 1:size(S,1)
        c = 0.0
        cm[i] = abs2(S[i,1])*sign(w[(1-1) % m + 1])
        for j = 2:size(S,2)
            t = cm[i] + abs2(S[i,j])*sign(w[(j-1) % m + 1])
            if abs(cm[i]) >= abs2(S[i,j])
                c += (cm[i]-t) + abs2(S[i,j])*sign(w[(j-1) % m + 1])
            else
                c += (abs2(S[i,j])*sign(w[(j-1) % m + 1])-t) + cm[i]
            end
            cm[i] = t
        end

        for j = 1:size(Snoise,2)
            t = cm[i] + abs2(Snoise[i,j])*sign(w[(j-1) % m + 1])
            if abs(cm[i]) >= abs2(Snoise[i,j])
                c += (cm[i]-t) + abs2(Snoise[i,j])*sign(w[(j-1) % m + 1])
            else
                c += (abs2(Snoise[i,j])*sign(w[(j-1) % m + 1])-t) + cm[i]
            end
            cm[i] = t
        end
        cm[i]+=c
    end

    return nothing
end

"""
    calcqe(S)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. 
"""
function calcqe(S)

    qe = zeros(Float64,size(S))

    calcqe!(qe,S)

    return qe
end

"""
    calcqe!(qe,S)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. Overwrites qe with output. 
"""
function calcqe!(qe,S)

    if any(size(qe) .!= size(S))
        throw(DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."))
    end

    for i = 1:size(S,1)
        denom=0
        for j = 1:size(S,2)
            denom += abs2(S[i,j])
        end
        for j = 1:size(S,2)
            qe[i,j] = abs2(S[i,j]) / denom
        end

        # qe[i,:] = abs2.(S[i,:]) ./ denom
    end

    return nothing
end

"""
    calcqe(S,Snoise)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. 
"""
function calcqe(S,Snoise)

    qe = zeros(Float64,size(S))

    calcqe!(qe,S,Snoise)

    return qe
end

"""
    calcqe!(qe,S,Snoise)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. Overwrites qe with output. 
"""
function calcqe!(qe,S,Snoise)

    if any(size(qe) .!= size(S))
        throw(DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."))
    end

    if size(S,1) .!= size(Snoise,1)
        throw(DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
    end

    for i = 1:size(S,1)
        denom=0
        for j = 1:size(S,2)
            denom += abs2(S[i,j])
        end

        for j = 1:size(Snoise,2)
            denom += abs2(Snoise[i,j])
        end
        for j = 1:size(S,2)
            qe[i,j] = abs2(S[i,j]) / denom
        end
        # qe[i,:] = abs2.(S[i,:]) ./ denom
    end

    return nothing
end

