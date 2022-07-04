


"""
  calcphibports(nodevalues,branches,Nmodes)

Calculate the branch values (variables) at the branches specified in the array 
of tuples which contain the two nodes that make up the branch. Nmodes gives
the number of frequency modes.

# Examples
```jldoctest
julia> JosephsonCircuits.calcbranchvalues([1 2 3;4 5 6],[(2,1)],1)
1×3 Matrix{Int64}:
 1  2  3
```
"""
function calcbranchvalues(nodevalues,branches,Nmodes)

    # branchvalues = Matrix(eltype(nodevalues),length(branches)*Nmodes,size(nodevalues,2))
    branchvalues = Matrix{eltype(nodevalues)}(undef,length(branches)*Nmodes,size(nodevalues,2))

    calcbranchvalues!(branchvalues,nodevalues,branches,Nmodes)

    return branchvalues

end


"""
    calcbranchvalues!(branchvalues,nodevalues,branches,Nmodes)

Calculate the branch values (variables) at the branches specified in the array 
of tuples which contain the two nodes that make up the branch. Return the result
in branchvalues. Nmodes gives the number of frequency modes.

Note: I think this is using the wrong sign convention.
In Kerman et al.  
phib_{branch_{i,j}} = phin_j - phin_i
In Rasmussen it's end - start.

I need to find out where else this is swapped because changing the signs here
messes up the S11 calculation. 

"""
function calcbranchvalues!(branchvalues,nodevalues,branches,Nmodes)

    if size(branchvalues,2) != size(nodevalues,2)
        error("Error: inconsistent second dimension for branchflux and nodeflux.")
    end

    if size(branchvalues,1) != length(branches)*Nmodes
        error("Error: inconsistent dimensions for branchflux and number of ports and modes.")
    end

    for (i,key) in enumerate(branches)
        # calculate the branch flux from the node flux. if there are flux
        # quanta in the loop this might be incorrect. if any of the keys
        # are equal to 1, then it is the ground node and that flux will be
        # zero

        if key[1] == 1
            for j in 1:Nmodes
                for k in 1:size(nodevalues,2)
                    branchvalues[(i-1)*Nmodes+j,k] = -nodevalues[(key[2]-2)*Nmodes+j,k]
                end
            end
        elseif key[2] == 1
            for j in 1:Nmodes
                for k in 1:size(nodevalues,2)
                    branchvalues[(i-1)*Nmodes+j,k] =  nodevalues[(key[1]-2)*Nmodes+j,k]
                end
            end
        else
            for j in 1:Nmodes
                for k in 1:size(nodevalues,2)
                    branchvalues[(i-1)*Nmodes+j,k] =  nodevalues[(key[1]-2)*Nmodes+j,k] 
                    branchvalues[(i-1)*Nmodes+j,k] -= nodevalues[(key[2]-2)*Nmodes+j,k]
                end
            end
        end
    end
    return nothing
end

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


"""

"""
function calcphibthevenin!(input,bnm,resistordict,wmodes,symfreqvar)

    Nmodes = length(wmodes)

    if size(input,1) != length(resistordict)*Nmodes
        throw(DimensionMismatch("First axis of input matrix must be number of ports times number of modes"))
    end

    if size(input,2) != size(bnm,2)
        throw(DimensionMismatch("Size of second dimension of input matrix must be the same as the size of the second dimension of bnm"))
    end

    # loop over the ports
    for (j,(key,val)) in enumerate(resistordict)
        # loop over the modes
        for k = 1:Nmodes
            # if symbolic, convert to numeric
            if val isa Symbolic
                if !(symfreqvar isa Symbolic)
                    error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                end
                resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
            else
                resistance = val
            end

            # should i loop over bnm instead?
            # then take 1:Nmodes terms for each of those. 
            # this way i don't need to know where the inputs are. 
            # 
            # 

            # # find the right term in bnm
            # if key[1] == 1
            #     Ip = -bnm[(key[2]-2)*Nmodes+k,(j-1)*Nmodes+k]
            # elseif key[2] == 1
            #     Ip =  bnm[(key[1]-2)*Nmodes+k,(j-1)*Nmodes+k]
            # else
            #     Ip =  bnm[(key[1]-2)*Nmodes+k,(j-1)*Nmodes+k] 
            #     Ip -= bnm[(key[2]-2)*Nmodes+k,(j-1)*Nmodes+k]
            # end

            # find the right term in bnm
            if key[1] == 1
                Ip = -bnm[(key[2]-2)*Nmodes+k,:]
            elseif key[2] == 1
                Ip =  bnm[(key[1]-2)*Nmodes+k,:]
            else
                Ip =  bnm[(key[1]-2)*Nmodes+k,:] 
                Ip .-= bnm[(key[2]-2)*Nmodes+k,:]
            end

            input[(j-1)*Nmodes+k,:] .= 
                Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
                # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))

            # input[(j-1)*Nmodes+k,(j-1)*Nmodes+k] = 
            #     Ip*resistance/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
            #     # 1*resistance/phi0/sqrt(resistance)/2/sqrt(abs(wmodes[k]))
        end
    end

    return nothing
end


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

"""


portimpedance can be for several sets of independent inputs. 
second dimension is arbitrary
first dimension is length(resistorvalues)*Nmodes

"""
function calcportimpedance!(portimpedance,resistorvalues,wmodes,symfreqvar)

    Nmodes = length(wmodes)

    # if size(branchvalues,2) != size(nodevalues,2)
    #     error("Error: inconsistent second dimension for branchflux and nodeflux.")
    # end

    if size(portimpedance,1) != length(resistorvalues)*Nmodes
        error("Error: inconsistent dimensions for portimpedance and number of resistorvalues and modes.")
    end


    for (j,val) in enumerate(resistorvalues)
        for k = 1:Nmodes
            # if symbolic, convert to numeric
            # i might want to put the check for type outside of this loop
            if val isa Symbolic
                if !(symfreqvar isa Symbolic)
                    error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                end
                resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
            else
                resistance = val
            end

            # take the complex conjugate if frequency is negative
            if wmodes[k] < 0
                resistance = conj(resistance)
            end

            for l = 1:size(portimpedance,2)
                # current[(j-1)*Nmodes+k,l] *= resistance/sqrt(resistance)/sqrt(abs(wmodes[k]))
                portimpedance[(j-1)*Nmodes+k,l] = resistance
            end
        end
    end

    return nothing
end


function scaleby!(portimpedance,wmodes)

    Nmodes = length(wmodes)
    Nbranches = div(size(portimpedance,1),Nmodes)

    if size(portimpedance,1) != Nbranches*Nmodes
        error("length of first axis of portimpedance must be an integer multiple of Nmodes")
    end

    for i = 1:Nbranches
        for j = 1:Nmodes
            for k = 1:size(portimpedance,2)
                portimpedance[(i-1)*Nmodes+j,k] *= wmodes[j]
            end
        end
    end

    return nothing
end

# """


# Convert current from amperes to units of sqrt energy or sqrt(photons/second), 
# up to a factor of Planck's constant which we will not apply because it cancels
# out (and is a very small number).
# """
# function currenttosqrtphotonflux!(current,resistorvalues,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     for (j,val) in enumerate(resistorvalues)
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

#             # take the complex conjugate if frequency is negative
#             if wmodes[k] < 0
#                 resistance = conj(resistance)
#             end

#             for l = 1:size(current,2)
#                 # current[(j-1)*Nmodes+k,l] *= resistance/sqrt(resistance)/sqrt(abs(wmodes[k]))
#                 current[(j-1)*Nmodes+k,l] *= (resistance)/sqrt(real(resistance))/sqrt(abs(wmodes[k]))

#             end
#         end
#     end

#     return nothing
# end


# function currenttosqrtphotonfluxconj!(current,resistorvalues,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     for (j,val) in enumerate(resistorvalues)
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

#             # take the complex conjugate if frequency is negative
#             if wmodes[k] < 0
#                 resistance = conj(resistance)
#             end

#             for l = 1:size(current,2)
#                 # current[(j-1)*Nmodes+k,l] *= resistance/sqrt(resistance)/sqrt(abs(wmodes[k]))
#                 current[(j-1)*Nmodes+k,l] *= conj(resistance)/sqrt(real(resistance))/sqrt(abs(wmodes[k]))

#             end
#         end
#     end

#     return nothing
# end


# """

# Convert flux from units of volt / frequency to units of sqrt energy or sqrt(photons/second), up to a factor
# of Planck's constant which we will not apply because it cancels out (and is a 
# very small number).

# """
# function fluxtosqrtphotonflux!(flux,resistorvalues,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     for (j,val) in enumerate(resistorvalues)
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

#             for l = 1:size(flux,2)
#                 # flux[(j-1)*Nmodes+k,l] *= im*wmodes[k]/sqrt(resistance)/sqrt(abs(wmodes[k]))
#                 flux[(j-1)*Nmodes+k,l] *= im*wmodes[k]/sqrt(real(resistance))/sqrt(abs(wmodes[k]))
#             end
#         end
#     end

#     return nothing
# end


function calcoutput!(output,phibports,resistorvalues,wmodes,symfreqvar)

    Nmodes = length(wmodes)

    if size(output,2) != size(phibports,2)
        error("Error: inconsistent dimensions for output and phibports.")
    end

    # for (j,(key,val)) in enumerate(resistordict)
    for (j,val) in enumerate(resistorvalues)
        # for k = 1:Nmodes
        #     output[(j-1)*Nmodes+k,:] .= 
        #         im*wmodes[k]*phibports[(j-1)*Nmodes+k,:]./sqrt(resistordict[key])/sqrt(abs(wmodes[k]))
        # end
        for k = 1:Nmodes
            # if symbolic, convert to numeric
            if val isa Symbolic
                if !(symfreqvar isa Symbolic)
                    error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                end
                resistance = substitute(val,Dict(symfreqvar=>wmodes[k]))
            else
                resistance = val
            end

            for l = 1:size(output,2)
                output[(j-1)*Nmodes+k,l] = 
                    im*wmodes[k]*phibports[(j-1)*Nmodes+k,l]/sqrt(resistance)/sqrt(abs(wmodes[k]))
            end
        end
    end

    return nothing
end

# function calcoutput!(output,phibports,portdict,resistordict,wmodes,symfreqvar)

#     Nmodes = length(wmodes)

#     if size(output,2) != size(phibports,2)
#         error("Error: inconsistent dimensions for output and phibports.")
#     end

#     j=0
#     # i should have a case for when there is no resistor at the port. 
#     for (key,val) in portdict
#         j+=1
#         # for k = 1:Nmodes
#         #     output[(j-1)*Nmodes+k,:] .= 
#         #         im*wmodes[k]*phibports[(j-1)*Nmodes+k,:]./sqrt(resistordict[key])/sqrt(abs(wmodes[k]))
#         # end
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

#             for l = 1:size(output,2)
#                 output[(j-1)*Nmodes+k,l] = 
#                     im*wmodes[k]*phibports[(j-1)*Nmodes+k,l]/sqrt(resistance)/sqrt(abs(wmodes[k]))
#             end
#         end
#     end

#     return nothing
# end

"""
  calcS(phibports,resistordict,wmodes)

Return the generalized scattering parameters for the system linearized around
the strong pump. 
"""
# function calcS(input,output,)

#   S = zeros(Complex{Float64},Nports*Nmodes,Nports*Nmodes)

#   calcS!(S,phibports,resistordict,wmodes)

#   return S

# end

"""
  calcS!

Return the scattering parameters for the system linearized around the strong pump. 
"""
function calcS!(S,inputwave,outputwave,phin,bnm,inputportbranches,outputportbranches,
    inputportimpedance,outputportimpedance,wmodes,symfreqvar)

    # check the size of S

    # check the size of inputwave

    # check the size of outputwave

    # check the sizes of all of the inputs

    # loop over input branches and modes to define inputwaves
    Ninputports = length(inputportbranches)
    Noutputports = length(outputportbranches)
    Nsolutions = size(phin,2)
    Nmodes = length(wmodes)

    for i = 1:Ninputports
        for j = 1:Nmodes
            for k = 1:Nsolutions

                key = inputportbranches[i]

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

                # calculate the port impedance
                if inputportimpedance[i] isa Symbolic
                    if !(symfreqvar isa Symbolic)
                        error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                    end
                    portimpedance = substitute(inputportimpedance[i],Dict(symfreqvar=>wmodes[j]))
                else
                    portimpedance = inputportimpedance[i]
                end

                # take the complex conjugate if frequency is negative
                if wmodes[j] < 0
                    portimpedance = conj(portimpedance)
                end


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

                key = outputportbranches[i]

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

                # calculate the port impedance
                if outputportimpedance[i] isa Symbolic
                    if !(symfreqvar isa Symbolic)
                        error("Error: Set symfreqvar equal to the symbolic variable representing frequency.")
                    end
                    portimpedance = substitute(outputportimpedance[i],Dict(symfreqvar=>wmodes[j]))
                else
                    portimpedance = outputportimpedance[i]
                end

                # take the complex conjugate if frequency is negative
                if wmodes[j] < 0
                    portimpedance = conj(portimpedance)
                end

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


# """
#     calcqe(S,Sn)

# Calculate the quantum efficiency matrix for a scattering matrix S for the input
# and output ports and an additional scattering matrix Sn representing waves at the
# output port from additional input ports (ex quantum noise from lossy capacitors 
# modelled as infinite transmission lines). S and Sn are in the field ladder 
# operator basis. 
# """
# function calcqe(S::AbstractArray{Complex{Float64},2},
#     Sn::AbstractArray{Complex{Float64},2})

#     qe = zeros(Float64,size(S))

#     calcqe!(S,Sn,qe)

#     return qe
# end

# """
#     calcqe(S,Sn,qe)

# Calculate the quantum efficiency matrix for a scattering matrix S for the input
# and output ports and an additional scattering matrix Sn representing waves at the
# output port from additional input ports (ex quantum noise from lossy capacitors 
# modelled as infinite transmission lines). S and Sn are in the field ladder 
# operator basis. Overwrites qe with output. 
# """
# function calcqe!(S::AbstractArray{Complex{Float64},2},
#     Sn::AbstractArray{Complex{Float64},2},qe::AbstractArray{Float64,2})

#     if size(S)[1] != size(Sn)[1]
#         error("Scattering matrix and noise scattering matrix must be the same height")
#     end

#     for i = 1:size(S)[1]
#         denom:: Float64 = 0.0
#         for j = 1:size(S)[2]
#             denom += abs2.(S[i,j])
#         end
#         for j = 1:size(Sn)[2]
#             denom += abs2.(Sn[i,j])
#         end
#         qe[i,:] = abs2.(S[i,:]) ./ denom
#     end

#     return nothing
# end

# function calcqe!(S::AbstractArray{Complex{Float64},2},
#         Sn1::AbstractArray{Complex{Float64},2},Sn2::AbstractArray{Complex{Float64},2},
#         qe::AbstractArray{Float64,2})

#     if size(S)[1] != size(Sn1)[1]
#         error("Scattering matrix and noise scattering matrix must be the same height")
#     end
#     if size(S)[1] != size(Sn2)[1]
#         error("Scattering matrix and noise scattering matrix must be the same height")
#     end

#     for i = 1:size(S)[1]
#         denom:: Float64 = 0.0
#         for j = 1:size(S)[2]
#             denom += abs2.(S[i,j])
#         end
#         for j = 1:size(Sn1)[2]
#             denom += abs2.(Sn1[i,j])
#         end
#         for j = 1:size(Sn2)[2]
#             denom += abs2.(Sn2[i,j])
#         end
#         qe[i,:] = abs2.(S[i,:]) ./ denom
#     end

#     return nothing
# end




# """
#     SphitoS(S,w,Zin,Zout)

# Convert from a scattering matrix in the node flux basis to the field ladder
# operator basis. 
# """
# function SphitoS(Sphi::AbstractArray{Complex{Float64},2},w::AbstractArray{Float64,1},
#     Zin::AbstractArray{Float64,1},Zout::AbstractArray{Float64,1})

#     S = zeros(Complex{Float64},size(Sphi))

#     copy!(S,Sphi)
#     SphitoS!(S,w,Zin,Zout)

#     return S
# end

# """
#     SphitoS!(Sphi,w,Zin,Zout,S)

# Convert from a scattering matrix in the node flux basis to the field ladder
# operator basis. Overwrites S with output. 
# """
# function SphitoS!(Sphi::AbstractArray{Complex{Float64},2},w::AbstractArray{Float64,1},
#     Zin::AbstractArray{Float64,1},Zout::AbstractArray{Float64,1},
#     S::AbstractArray{Complex{Float64},2})

#     if any(size(S) .!= size(Sphi))
#         error("Dimensions of input and output scattering matrices must be identical")
#     end

#     copy!(S,Sphi)
#     SphitoS!(S,w,Zin,Zout)

#     return nothing
# end

# """
#     SphitoS!(S,w,Zin,Zout)

# Convert from a scattering matrix in the node flux basis to the field ladder
# operator basis. Overwrites S with output. 
# """
# function SphitoS!(Sphi::AbstractArray{Complex{Float64},2},w::AbstractArray{Float64,1},
#     Zin::AbstractArray{Float64,1},Zout::AbstractArray{Float64,1})

#     m = length(w)

#     #if m > 1 && all(x->x >= 0, w)
#     #    error("If m>1 then there should be a negative frequency idler. Are you
#     #    using vector of signal frequencies instead?")
#     #end

#     if any(mod.(size(Sphi),m) .!=0)
#         error("Dimensions of scattering matrix must be an integer multiple of the number of frequencies")
#     end

#     if size(Sphi)[2] != length(Zin) 
#         error("Input impedance length incorrect")
#     end

#     if size(Sphi)[1] != length(Zout)
#         error("Output impedance length incorrect")
#     end


#     for j = 1:size(Sphi)[2]
#         for i = 1:size(Sphi)[1]
#             Sphi[i,j] *= sqrt(abs(w[(i-1) % m + 1]/w[(j-1) % m + 1]))*sqrt(Zin[j]/Zout[i])
#         end
#     end
#     return nothing
# end


# """
#     StoSphi(S,w,Zin,Zout)

# Convert from a scattering matrix in the field ladder operator basis to the node 
# flux basis.
# """
# function StoSphi(S::AbstractArray{Complex{Float64},2},w::AbstractArray{Float64,1},
#     Zin::AbstractArray{Float64,1},Zout::AbstractArray{Float64,1})

#     Sphi = zeros(Complex{Float64},size(S))

#     copy!(Sphi,S)
#     SphitoS!(Sphi,w,Zin,Zout)

#     return Sphi
# end

# """
#     StoSphi!(S,w,Zin,Zout,Sphi)

# Convert from a scattering matrix in the field ladder operator basis to the node 
# flux basis. Overwrites Sphi with output. 
# """
# function StoSphi!(S::AbstractArray{Complex{Float64},2},w::AbstractArray{Float64,1},
#     Zin::AbstractArray{Float64,1},Zout::AbstractArray{Float64,1},
#     Sphi::AbstractArray{Complex{Float64},2})

#     if any(size(Sphi) .!= size(S))
#         error("Dimensions of input and output scattering matrices must be identical")
#     end

#     copy!(Sphi,S)
#     StoSphi!(Sphi,w,Zin,Zout)

#     return nothing
# end

# """
#     StoSphi!(S,w)

# Convert from a scattering matrix in the field ladder operator basis to the node 
# flux basis. Overwrites S with output. 
# """
# function StoSphi!(S::AbstractArray{Complex{Float64},2},w::AbstractArray{Float64,1},
# Zin::AbstractArray{Float64,1},Zout::AbstractArray{Float64,1})

#     m = length(w)

#     if any(mod.(size(S),m) .!=0)
#         error("Dimensions of scattering matrix must be an integer multiple of the number of frequencies")
#     end

#     for j = 1:size(S)[2]
#         for i = 1:size(S)[1]
#             #S[i,j] /= sqrt.(abs.(w[(i-1) % m + 1]/w[(j-1) % m + 1]))
#             S[i,j] /= sqrt(abs(w[(i-1) % m + 1]/w[(j-1) % m + 1]))*sqrt(Zin[j]/Zout[i])
#         end
#     end

#     return nothing
# end

