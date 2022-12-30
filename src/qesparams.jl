
"""
  calcS!

Return the scattering parameters for the system linearized around the strong pump. 
"""
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

    for i in 1:Ninputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

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

                # calculate the port impedance
                portimpedance = calcimpedance(inputportimpedances[i],typevector[inputportindices[i]],wmodes[j],symfreqvar)

                # calculate the current flowing through the port
                portcurrent = sourcecurrent - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(Complex(real(portimpedance)))


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
    for i in 1:Noutputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

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

                # calculate the port impedance
                portimpedance = calcimpedance(outputportimpedances[i],typevector[outputportindices[i]],wmodes[j],symfreqvar)

                # calculate the current flowing through the port
                portcurrent = sourcecurrent - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(Complex(real(portimpedance)))

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
  calcSnoise!

Return the scattering parameters for the system linearized around the strong pump. 

This is a bit of a hack but I ran into issues with complex capacitance when
the capacitor was at the same branch as a current source. the calcS function
would use that current source in calculating the output waves, which it should
not do.
"""
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

    for i in 1:Ninputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

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

                # calculate the port impedance
                portimpedance = calcimpedance(inputportimpedances[i],typevector[inputportindices[i]],wmodes[j],symfreqvar)

                # calculate the current flowing through the port
                portcurrent = sourcecurrent - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(Complex(real(portimpedance)))

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
    for i in 1:Noutputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

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

                # calculate the port impedance
                portimpedance = calcimpedance(outputportimpedances[i],typevector[outputportindices[i]],wmodes[j],symfreqvar)

                # calculate the current flowing through the port
                # portcurrent = sourcecurrent - portvoltage / portimpedance
                portcurrent = - portvoltage / portimpedance

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(Complex(real(portimpedance)))

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
    calcimpedance(c::Union{Integer,T,Complex{T}},type,w,symfreqvar) where {T<:AbstractFloat}

# Examples
```jldoctest
julia> JosephsonCircuits.calcimpedance(30.0,:C,1.0,nothing)
0.0 - 0.03333333333333333im

julia> JosephsonCircuits.calcimpedance(30.0,:R,1.0,nothing)
30.0 + 0.0im

julia> JosephsonCircuits.calcimpedance(30.0,:C,-1.0,nothing)
-0.0 + 0.03333333333333333im

julia> JosephsonCircuits.calcimpedance(30.0,:R,-1.0,nothing)
30.0 + 0.0im
```
"""
# function calcimpedance(c::Union{Integer,T,Complex{T}},type,w,symfreqvar) where {T<:AbstractFloat}
function calcimpedance(c::Union{T,Complex{T}},type,w,symfreqvar) where {T<:Union{AbstractFloat,Integer}}

    if type == :R
        if w >= 0
            return c+0.0im
        else
            return conj(c)+0.0im
        end
    elseif type == :C
        if w >= 0
            return 1/(im*w*c)
        else
            return 1/(im*w*conj(c))
        end
    else
        error("Unknown component type")
    end
end


"""
    calcimpedance(c,type,w,symfreqvar)

# Examples
```jldoctest
julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:R,2.0,w)
60.0 + 0.0im

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:C,2.0,w)
0.0 - 0.008333333333333333im

julia> @variables w;JosephsonCircuits.calcimpedance(30,:R,-2.0,w)
30.0 + 0.0im

julia> @variables w;JosephsonCircuits.calcimpedance(30,:C,-2.0,w)
-0.0 + 0.016666666666666666im
```
"""
function calcimpedance(c,type,w,symfreqvar)
    if type == :R
        if w >= 0
            return valuetonumber(c,symfreqvar => w)+0.0im
        else
            return conj(valuetonumber(c,symfreqvar => w))+0.0im
        end
    elseif type == :C
        if w >= 0
            return 1/(im*w*valuetonumber(c,symfreqvar => w))
        else
            return 1/(im*w*conj(valuetonumber(c,symfreqvar => w)))
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

# Examples
```jldoctest
julia> JosephsonCircuits.calccm(Complex{Float64}[3/5 4/5;4/5 3/5],[1])
2-element Vector{Float64}:
 1.0
 1.0

julia> JosephsonCircuits.calccm([1 1e-100 2e-100 1;1 0 0 1],[1 -1])
2-element Vector{Float64}:
 3.0e-200
 0.0

julia> @variables a b;JosephsonCircuits.calccm([a b; b a],[1 -1])
2-element Vector{Num}:
 abs2(a) - abs2(b)
 abs2(b) - abs2(a)
```
"""
function calccm(S::AbstractArray{T}, w) where {T}
    cm = zeros(T,size(S,1))
    calccm!(cm,S,w)
    return cm
end

function calccm(S::AbstractArray{Complex{T}}, w) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    cm = zeros(T,size(S,1))
    calccm!(cm,S,w)
    return cm
end

"""
    calccm!(cm,S,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Overwrites cm with output. Use a compensated sum
to reduce floating point errors.

# Examples
```jldoctest
julia> JosephsonCircuits.calccm(Complex{Float64}[1 1e-100 2e-100 1;1 1 1 1],Complex{Float64}[1 1e-100 2e-100 1;1 1 1 1],[1 -1])
2-element Vector{Float64}:
 6.0e-200
 0.0

julia> JosephsonCircuits.calccm([1 1e-100 2e-100 1;1 1 1 1],[1 1e-100 2e-100 1;1 1 1 1],[1 -1])
2-element Vector{Float64}:
 6.0e-200
 0.0
```
"""
function calccm!(cm::AbstractArray{T},S,w) where {T<:AbstractFloat}

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch("First dimension of scattering matrix must equal the length of cm."))        
    end

    # use a Kahan, Babushka, Neumaier compensated sum. more cache efficient version
    fill!(cm,zero(T))
    c = zeros(eltype(cm),size(cm))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            t = cm[i] + abs2(S[i,j])*sign(w[(j-1) % m + 1])
            c[i] += ifelse(
                abs(cm[i]) >= abs2(S[i,j]),
                (cm[i]-t) + abs2(S[i,j])*sign(w[(j-1) % m + 1]),
                (abs2(S[i,j])*sign(w[(j-1) % m + 1])-t) + cm[i])
            cm[i]  = t
        end
    end

    @inbounds for i in 1:size(S,1)
        cm[i]+=c[i]
    end


    return nothing
end

"""
    calccm!(cm,S,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Overwrites cm with output. 

"""
function calccm!(cm,S,w)

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch("First dimension of scattering matrix must equal the length of cm."))        
    end

    # more cache friendly version
    fill!(cm,zero(eltype(cm)))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            cm[i] += abs2(S[i,j])*sign(w[(j-1) % m + 1])
        end
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
function calccm(S::AbstractArray{T},Snoise::AbstractArray{T},w) where {T}
    cm = zeros(T,size(S,1))
    calccm!(cm,S,Snoise,w)
    return cm
end

function calccm(S::AbstractArray{Complex{T}},Snoise::AbstractArray{Complex{T}},w) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    cm = zeros(T,size(S,1))
    calccm!(cm,S,Snoise,w)
    return cm
end


"""
    calccm!(cm,S,Snoise,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Overwrites cm with output.  Use a compensated sum
to reduce floating point errors.

"""
function calccm!(cm::AbstractArray{T},S,Snoise,w) where {T<:AbstractFloat}

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    for d in size(Snoise)
        if mod(d, m) != 0
            throw(DimensionMismatch("Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch("First dimension of scattering matrix must equal the length of cm."))        
    end

    if size(S,1) != size(Snoise,1)
        throw(DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
    end

    # use a Kahan, Babushka, Neumaier compensated sum. more cache efficient version
    fill!(cm,zero(T))
    c = zeros(eltype(cm),size(cm))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            t = cm[i] + abs2(S[i,j])*sign(w[(j-1) % m + 1])
            c[i] += ifelse(
                abs(cm[i]) >= abs2(S[i,j]),
                (cm[i]-t) + abs2(S[i,j])*sign(w[(j-1) % m + 1]),
                (abs2(S[i,j])*sign(w[(j-1) % m + 1])-t) + cm[i])
            cm[i]  = t
        end
    end

    @inbounds for j in 1:size(Snoise,2)
        for i in 1:size(Snoise,1)
            t = cm[i] + abs2(Snoise[i,j])*sign(w[(j-1) % m + 1])
            c[i] += ifelse(
                abs(cm[i]) >= abs2(Snoise[i,j]),
                (cm[i]-t) + abs2(Snoise[i,j])*sign(w[(j-1) % m + 1]),
                (abs2(Snoise[i,j])*sign(w[(j-1) % m + 1])-t) + cm[i])
            cm[i]  = t
        end
    end

    @inbounds for i in 1:size(S,1)
        cm[i]+=c[i]
    end

    return nothing
end


"""
    calccm!(cm,S,Snoise,w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Overwrites cm with output.

"""
function calccm!(cm,S,Snoise,w)

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    for d in size(Snoise)
        if mod(d, m) != 0
            throw(DimensionMismatch("Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch("First dimension of scattering matrix must equal the length of cm."))        
    end

    if size(S,1) != size(Snoise,1)
        throw(DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
    end

    # more cache friendly version
    fill!(cm,zero(eltype(cm)))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            cm[i] += abs2(S[i,j])*sign(w[(j-1) % m + 1])
        end
    end

    @inbounds for j in 1:size(Snoise,2)
        for i in 1:size(Snoise,1)
            cm[i] += abs2(Snoise[i,j])*sign(w[(j-1) % m + 1])
        end
    end

    return nothing
end


"""
    calcqe(S)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. 

# Examples
```jldoctest
julia> JosephsonCircuits.calcqe([3/5 4/5;4/5 3/5])
2×2 Matrix{Float64}:
 0.36  0.64
 0.64  0.36

julia> JosephsonCircuits.calcqe(Complex{Float64}[3/5 4/5;4/5 3/5])
2×2 Matrix{Float64}:
 0.36  0.64
 0.64  0.36

julia> @variables a b c d;JosephsonCircuits.calcqe([a b; c d])
2×2 Matrix{Num}:
 abs2(a) / (abs2(a) + abs2(b))  abs2(b) / (abs2(a) + abs2(b))
 abs2(c) / (abs2(c) + abs2(d))  abs2(d) / (abs2(c) + abs2(d))
```
"""
function calcqe(S::AbstractArray{T}) where {T}
    qe = zeros(T,size(S))
    calcqe!(qe,S)
    return qe
end

function calcqe(S::AbstractArray{Complex{T}}) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    qe = zeros(T,size(S))
    calcqe!(qe,S)
    return qe
end

"""
    calcqe!(qe,S)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. Overwrites qe with output. 
"""
function calcqe!(qe,S)

    if size(qe) != size(S)
        throw(DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."))
    end

    # more cache efficient version of QE calculation
    denom = zeros(eltype(qe),size(S,1))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            denom[i] += abs2(S[i,j])
        end
    end

    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            qe[i,j] = abs2(S[i,j]) / denom[i]
        end
    end

    return nothing
end

"""
    calcqe(S,Snoise)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. 

# Examples
```jldoctest
julia> @variables a b c d an bn cn dn;JosephsonCircuits.calcqe([a b; c d],[an bn; cn dn])
2×2 Matrix{Num}:
 abs2(a) / (abs2(a) + abs2(an) + abs2(b) + abs2(bn))  …  abs2(b) / (abs2(a) + abs2(an) + abs2(b) + abs2(bn))
 abs2(c) / (abs2(c) + abs2(cn) + abs2(d) + abs2(dn))     abs2(d) / (abs2(c) + abs2(cn) + abs2(d) + abs2(dn))
```
"""
function calcqe(S::AbstractArray{T},Snoise::AbstractArray{T}) where {T}
    qe = zeros(T,size(S))
    calcqe!(qe,S,Snoise)
    return qe
end

function calcqe(S::AbstractArray{Complex{T}},Snoise::AbstractArray{Complex{T}}) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    qe = zeros(T,size(S))
    calcqe!(qe,S,Snoise)
    return qe
end

"""
    calcqe!(qe,S,Snoise)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. Overwrites qe with output. 

"""
function calcqe!(qe,S,Snoise)

    if size(qe) != size(S)
        throw(DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."))
    end

    if size(S,1) != size(Snoise,1)
        throw(DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
    end

    # more cache efficient version of QE calculation
    denom = zeros(eltype(qe),size(S,1))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            denom[i] += abs2(S[i,j])
        end
    end

    @inbounds for j in 1:size(Snoise,2)
        for i in 1:size(Snoise,1)
            denom[i] += abs2(Snoise[i,j])
        end
    end

    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            qe[i,j] = abs2(S[i,j]) / denom[i]
        end
    end

    return nothing
end


"""
    calcqeideal(S)

Calculate the ideal (best possible) quantum efficiency for each element of a
scattering matrix. 

# Examples
```jldoctest
julia> JosephsonCircuits.calcqeideal([3/5 4/5;4/5 3/5])
2×2 Matrix{Float64}:
 1.0  1.0
 1.0  1.0

julia> JosephsonCircuits.calcqeideal(Complex{Float64}[3/5 4/5;4/5 3/5])
2×2 Matrix{Float64}:
 1.0  1.0
 1.0  1.0

julia> @variables S11 S12 S21 S22;println(JosephsonCircuits.calcqeideal([S11 S12;S21 S22]))
Num[ifelse(abs2(S11) <= 1, 1, 1 / (2 + -1 / abs2(S11))) ifelse(abs2(S12) <= 1, 1, 1 / (2 + -1 / abs2(S12))); ifelse(abs2(S21) <= 1, 1, 1 / (2 + -1 / abs2(S21))) ifelse(abs2(S22) <= 1, 1, 1 / (2 + -1 / abs2(S22)))]
 ```
"""
function calcqeideal(S::AbstractArray{T}) where {T}
    qeideal = zeros(T,size(S))
    calcqeideal!(qeideal,S)
    return qeideal
end

function calcqeideal(S::AbstractArray{Complex{T}}) where {T}
    # quantum efficiency is real so if the type of complex, use this
    # parametric method to define a real matrix.
    qeideal = zeros(T,size(S))
    calcqeideal!(qeideal,S)
    return qeideal
end

function calcqeideal!(qeideal,S)
    if size(qeideal) != size(S)
        throw(DimensionMismatch("Sizes of QE and S matrices must be equal."))
    end
    for i in eachindex(S)
        abs2S = abs2(S[i])
        qeideal[i] = ifelse(abs2S <= 1,one(eltype(qeideal)),1 /(2 - 1 /abs2S))
    end
    return nothing
end