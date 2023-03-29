
"""
    calcinputoutput!(inputwave, outputwave, phin, bnm, inputportindices,
        outputportindices, inputportimpedances, outputportimpedances,
        nodeindices, componenttypes, wmodes, symfreqvar)

Return the input and output waves for the system linearized around the strong pump.

# Examples
```jldoctest
inputwave = JosephsonCircuits.LinearAlgebra.Diagonal(ComplexF64[0])
outputwave = ComplexF64[0;;]
bnm = ComplexF64[1; 0;;]
portimpedanceindices = [3]
portimpedances = ComplexF64[50]
nodeindices = [2 2 2 2 0 3 3; 1 1 1 1 0 1 1]
componenttypes = [:P, :I, :R, :L, :K, :L, :C]
wmodes = [1]
phin = ComplexF64[0;0;;]
symfreqvar = nothing
JosephsonCircuits.calcinputoutput!(inputwave,outputwave,phin,bnm,portimpedanceindices,
    portimpedanceindices,portimpedances,portimpedances,nodeindices,componenttypes,
    wmodes,symfreqvar)
println(inputwave)
println(outputwave)

# output
ComplexF64[3.5355339059327378 + 0.0im;;]
ComplexF64[-3.5355339059327378 + 0.0im;;]
```
```jldoctest
inputwave = JosephsonCircuits.LinearAlgebra.Diagonal(ComplexF64[0])
outputwave = ComplexF64[0;;]
bnm = ComplexF64[1; 0;;]
portimpedanceindices = [3]
portimpedances = ComplexF64[50]
nodeindices = [2 2 2 2 0 3 3; 1 1 1 1 0 1 1]
componenttypes = [:P, :I, :R, :L, :K, :L, :C]
wmodes = [1]
phin = ComplexF64[50/(im*wmodes[1]);0;;]
symfreqvar = nothing
JosephsonCircuits.calcinputoutput!(inputwave,outputwave,phin,bnm,portimpedanceindices,
    portimpedanceindices,portimpedances,portimpedances,nodeindices,componenttypes,
    wmodes,symfreqvar)
println(inputwave)
println(outputwave)

# output
ComplexF64[3.5355339059327378 + 0.0im;;]
ComplexF64[3.5355339059327378 + 0.0im;;]
```
```jldoctest
inputwave = JosephsonCircuits.LinearAlgebra.Diagonal(ComplexF64[0])
outputwave = ComplexF64[0;;]
bnm = ComplexF64[1; 0;;]
portimpedanceindices = [3]
portimpedances = ComplexF64[50]
nodeindices = [1 1 1 1 0 1 1; 2 2 2 2 0 3 3;]
componenttypes = [:P, :I, :R, :L, :K, :L, :C]
wmodes = [1]
phin = ComplexF64[50/(im*wmodes[1]);0;;]
symfreqvar = nothing
JosephsonCircuits.calcinputoutput!(inputwave,outputwave,phin,bnm,portimpedanceindices,
    portimpedanceindices,portimpedances,portimpedances,nodeindices,componenttypes,
    wmodes,symfreqvar)
println(inputwave)
println(outputwave)

# output
ComplexF64[-3.5355339059327378 - 0.0im;;]
ComplexF64[-3.5355339059327378 + 0.0im;;]
```
```jldoctest
inputwave = JosephsonCircuits.LinearAlgebra.Diagonal(ComplexF64[0])
outputwave = ComplexF64[0;;]
bnm = ComplexF64[-1; 1;;]
portimpedanceindices = [2]
portimpedances = ComplexF64[50.0 + 0.0im]
nodeindices = [2 2 2 2 3; 3 3 1 1 1]
componenttypes = [:P, :R, :L, :C, :C]
wmodes = [1]
phin = ComplexF64[0;0;;]
symfreqvar = nothing
JosephsonCircuits.calcinputoutput!(inputwave,outputwave,phin,bnm,portimpedanceindices,
    portimpedanceindices,portimpedances,portimpedances,nodeindices,componenttypes,
    wmodes,symfreqvar)
println(inputwave)
println(outputwave)

# output
ComplexF64[-7.0710678118654755 + 0.0im;;]
ComplexF64[7.0710678118654755 + 0.0im;;]
```
```jldoctest
inputwave = JosephsonCircuits.LinearAlgebra.Diagonal(ComplexF64[0])
outputwave = ComplexF64[0;;]
bnm = ComplexF64[-1; 1;;]
portimpedanceindices = [2]
portimpedances = ComplexF64[50.0 + 0.0im]
nodeindices = [2 2 2 2 3; 3 3 1 1 1]
componenttypes = [:P, :R, :L, :C, :C]
wmodes = [1]
phin = ComplexF64[-50/(im*wmodes[1]);50/(im*wmodes[1]);;]
symfreqvar = nothing
JosephsonCircuits.calcinputoutput!(inputwave,outputwave,phin,bnm,portimpedanceindices,
    portimpedanceindices,portimpedances,portimpedances,nodeindices,componenttypes,
    wmodes,symfreqvar)
println(inputwave)
println(outputwave)

# output
ComplexF64[-7.0710678118654755 + 0.0im;;]
ComplexF64[-7.0710678118654755 + 0.0im;;]
```
"""
function calcinputoutput!(inputwave, outputwave, phin, bnm, inputportindices,
    outputportindices, inputportimpedances, outputportimpedances,
    nodeindices, componenttypes, wmodes, symfreqvar)
    return calcinputoutput_inner!(inputwave, outputwave, phin, bnm, inputportindices,
        outputportindices, inputportimpedances, outputportimpedances,
        nodeindices, componenttypes, wmodes, symfreqvar, false)
end

"""
    calcinputoutputnoise!(S, inputwave, outputwave, phin, bnm, inputportindices,
        outputportindices, inputportimpedances, outputportimpedances,
        nodeindices, componenttypes, wmodes, symfreqvar)

Return the input and output waves for the system linearized around the strong pump. 

This is a bit of a hack but I ran into issues with complex capacitance when
the capacitor was at the same branch as a current source. the calcS function
would use that current source in calculating the output waves, which it should
not do.

# Examples
```jldoctest
inputwave = JosephsonCircuits.LinearAlgebra.Diagonal(ComplexF64[0])
noiseoutputwave = ComplexF64[0;;]
phin = ComplexF64[-2.5000000000007394e-10 - 0.000795774715459398im; 1.983790476804266e-20 + 3.141592641138603e-16im;;]
bnm = ComplexF64[1.0 + 0.0im; 0.0 + 0.0im;;]
portimpedanceindices = [2]
noiseportimpedanceindices = [6]
portimpedances = [50]
noiseportimpedances = [1]
nodeindices = [2 2 2 3 3 3; 1 1 3 1 1 1]
componenttypes = [:P, :R, :C, :Lj, :C, :R]
wmodes = [2*pi*5e9]
symfreqvar = nothing
JosephsonCircuits.calcinputoutputnoise!(inputwave,noiseoutputwave,
    phin,bnm,portimpedanceindices,noiseportimpedanceindices,
    portimpedances,noiseportimpedances,nodeindices,
    componenttypes,wmodes,symfreqvar)
println(inputwave)
println(noiseoutputwave)

# output
ComplexF64[1.994711402007163e-5 + 0.0im;;]
ComplexF64[-5.568327974762547e-11 + 3.516177070001411e-15im;;]
```
"""
function calcinputoutputnoise!(inputwave, outputwave, phin, bnm, inputportindices,
    outputportindices, inputportimpedances, outputportimpedances,
    nodeindices, componenttypes, wmodes, symfreqvar)
    return calcinputoutput_inner!(inputwave, outputwave, phin, bnm, inputportindices,
        outputportindices, inputportimpedances, outputportimpedances,
        nodeindices, componenttypes, wmodes, symfreqvar, true)
end


"""
    calcinputoutput_inner!(inputwave, outputwave, phin, bnm, inputportindices,
        outputportindices, inputportimpedances, outputportimpedances,
        nodeindices, componenttypes, wmodes, symfreqvar, nosource)

Calculate the input and output power waves as defined in (except in
units of sqrt(photons/second) instead of sqrt(power)
K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE Trans.
Micr. Theory and Tech. 13, 194–202 (1965) 
doi: 10.1109/TMTT.1965.1125964
inputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage + portimpedance * portcurrent)
we can simplify the above to:
inputwave[(i-1)*Nmodes+j,k] = 1/2*kval * portimpedance * sourcecurrent
outputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage - conj(portimpedance) * portcurrent)

"""
function calcinputoutput_inner!(inputwave, outputwave, nodeflux, bnm, inputportindices,
    outputportindices, inputportimpedances, outputportimpedances,
    nodeindices, componenttypes, wmodes, symfreqvar, nosource)

    # check the size of inputwave

    # check the size of outputwave

    # check the sizes of all of the inputs

    # loop over input branches and modes to define inputwaves
    Ninputports = length(inputportindices)
    Noutputports = length(outputportindices)
    Nsolutions = size(nodeflux,2)
    Nmodes = length(wmodes)

    for i in 1:Ninputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

                sourcecurrent = calcsourcecurrent(
                    nodeindices[1,inputportindices[i]],
                    nodeindices[2,inputportindices[i]],
                    bnm,Nmodes,j,k)

                # calculate the port impedance
                portimpedance = calcimpedance(
                    inputportimpedances[i],
                    componenttypes[inputportindices[i]],
                    wmodes[j],symfreqvar)

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(Complex(real(portimpedance)))

                # this will give NaN for DC, so set kval=0 in that case
                if wmodes[j] == 0
                    kval = 0
                else
                    kval *= 1 /sqrt(abs(wmodes[j]))
                end

                # calculate the input and output power waves as defined in (except in
                # units of sqrt(photons/second) instead of sqrt(power)
                # K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE Trans.
                # Micr. Theory and Tech. 13, 194–202 (1965) 
                # doi: 10.1109/TMTT.1965.1125964
                inputwave[(i-1)*Nmodes+j,k] = 1/2*kval * portimpedance * sourcecurrent
            end
        end
    end

    # loop over output branches and modes to define outputwaves
    for i in 1:Noutputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

                sourcecurrent = calcsourcecurrent(
                    nodeindices[1,outputportindices[i]],
                    nodeindices[2,outputportindices[i]],
                    bnm,Nmodes,j,k)

                portvoltage = calcportvoltage(
                    nodeindices[1,outputportindices[i]],
                    nodeindices[2,outputportindices[i]],
                    nodeflux,
                    wmodes,
                    Nmodes,j,k)

                # calculate the port impedance
                portimpedance = calcimpedance(
                    outputportimpedances[i],
                    componenttypes[outputportindices[i]],
                    wmodes[j],symfreqvar)

                # calculate the current flowing through the port
                if nosource
                    portcurrent = - portvoltage / portimpedance
                else
                    portcurrent = sourcecurrent - portvoltage / portimpedance
                end

                # calculate the scaling factor for the waves
                kval = 1 / sqrt(Complex(real(portimpedance)))

                # convert from sqrt(power) to sqrt(photons/second)
                # this will give NaN for DC, so set kval=0 in that case
                if wmodes[j] == 0
                    kval = 0
                else
                    kval *= 1 /sqrt(abs(wmodes[j]))
                end

                # calculate the input and output power waves as defined in (except in
                # units of sqrt(photons/second) instead of sqrt(power)
                # K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE Trans.
                # Micr. Theory and Tech. 13, 194–202 (1965) 
                # doi: 10.1109/TMTT.1965.1125964
                outputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage - conj(portimpedance) * portcurrent)
            end
        end
    end
 
    return nothing
end

"""
    calcscatteringmatrix!(S,inputwave::Diagonal,outputwave)

The scattering matrix is defined as outputwave = S * inputwave .

# Examples
```jldoctest
julia> inputwave=JosephsonCircuits.LinearAlgebra.Diagonal([1.0,1.0]);outputwave=[im/sqrt(2) 1/sqrt(2);1/sqrt(2) im/sqrt(2)];S = zeros(Complex{Float64},2,2);JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave);S
2×2 Matrix{ComplexF64}:
      0.0+0.707107im  0.707107+0.0im
 0.707107+0.0im            0.0+0.707107im
```
"""
function calcscatteringmatrix!(S,inputwave::Diagonal,outputwave)
    # copy!(S,outputwave)
    # rdiv!(S,inputwave)
    rdiv!(outputwave,inputwave)
    copy!(S,outputwave)
    
    return nothing
end

"""
    calcscatteringmatrix!(S,inputwave,outputwave)

The scattering matrix is defined as outputwave = S * inputwave .

# Examples
```jldoctest
julia> inputwave=[1.0 0.0;0.0 1.0];outputwave=[im/sqrt(2) 1/sqrt(2);1/sqrt(2) im/sqrt(2)];S = zeros(Complex{Float64},2,2);JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave);S
2×2 Matrix{ComplexF64}:
      0.0+0.707107im  0.707107+0.0im
 0.707107+0.0im            0.0+0.707107im

julia> inputwave = rand(Complex{Float64},2,2);outputwave = rand(Complex{Float64},2,2);S=zeros(Complex{Float64},2,2);JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave);isapprox(S*inputwave,outputwave)
true
```
"""
function calcscatteringmatrix!(S,inputwave,outputwave)
    S .= outputwave / inputwave
    return nothing
end

"""
    calcscatteringmatrix!(S,inputwave::Vector,outputwave::Vector)

The scattering matrix is defined as outputwave = S * inputwave .

# Examples
```jldoctest
julia> inputwave=[1.0,0.0];outputwave=[im/sqrt(2), 1/sqrt(2)];S = zeros(Complex{Float64},2,2);JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave);S
2×2 Matrix{ComplexF64}:
      0.0+0.707107im  0.0+0.0im
 0.707107+0.0im       0.0+0.0im
```
"""
function calcscatteringmatrix!(S,inputwave::Vector,outputwave::Vector)
    if size(S,1) != length(outputwave)
        throw(DimensionMismatch("First dimension of scattering matrix not consistent with first dimensions of outputwave."))
    end
    if size(S,2) != length(inputwave)
        throw(DimensionMismatch("Second dimension of scattering matrix not consistent with first dimension of input wave."))
    end

    fill!(S,0)
    for j in eachindex(inputwave)
        if !iszero(inputwave[j])
            for i in eachindex(outputwave)
                S[i,j] = outputwave[i]/inputwave[j]
            end
        end
    end
    return nothing
end

function calcportvoltage(key1,key2,phin,wmodes,Nmodes,j,k)

    # calculate the branch fluxes at the ports from the node flux array phin
    if key1 == 1
        portvoltage = -phin[(key2-2)*Nmodes+j,k]
    elseif key2 == 1
        portvoltage =  phin[(key1-2)*Nmodes+j,k]
    else
        portvoltage =  phin[(key1-2)*Nmodes+j,k] 
        portvoltage -= phin[(key2-2)*Nmodes+j,k]
    end

    # scale the branch flux by frequency to get voltage
    portvoltage *= im*wmodes[j]

    return portvoltage
end

function calcsourcecurrent(key1,key2,bnm,Nmodes,j,k)

    if key1 == 1
        sourcecurrent = -bnm[(key2-2)*Nmodes+j,k]
    elseif key2 == 1
        sourcecurrent =  bnm[(key1-2)*Nmodes+j,k]
    else
        sourcecurrent =  bnm[(key1-2)*Nmodes+j,k] 
        sourcecurrent -= bnm[(key2-2)*Nmodes+j,k]
    end
    return sourcecurrent
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

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:R,-2.0,w)
-60.0 + 0.0im

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:C,-2.0,w)
0.0 - 0.008333333333333333im
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

julia> JosephsonCircuits.calccm([1 1e-100 2e-100 1;1 0 0 1],[1, -1])
2-element Vector{Float64}:
 3.0e-200
 0.0

julia> @variables a b;JosephsonCircuits.calccm([a b; b a],[1, -1])
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
julia> cm=Float64[0,0];JosephsonCircuits.calccm!(cm,[3/5 4/5;4/5 3/5],[-1,1]);cm
2-element Vector{Float64}:
  0.28000000000000014
 -0.28000000000000014
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

# Examples
```jldoctest
julia> @variables a b;cm=Num[0,0];JosephsonCircuits.calccm!(cm,[a b; b a],[-1,1]);cm
2-element Vector{Num}:
 abs2(b) - abs2(a)
 abs2(a) - abs2(b)
```
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

# Examples
```jldoctest
julia> JosephsonCircuits.calccm([1 1e-100 2e-100 1;1 1 1 1],[1 1e-100 2e-100 1;1 1 1 1],[1, -1])
2-element Vector{Float64}:
 6.0e-200
 0.0

julia> JosephsonCircuits.calccm(Complex{Float64}[1 1e-100 2e-100 1;1 1 1 1],Complex{Float64}[1 1e-100 2e-100 1;1 1 1 1],[1, -1])
2-element Vector{Float64}:
 6.0e-200
 0.0

julia> @variables a b c d an bn cn dn;JosephsonCircuits.calccm([a b; c d],[an bn; cn dn],[1, -1])
2-element Vector{Num}:
 abs2(a) + abs2(an) - abs2(b) - abs2(bn)
 abs2(c) + abs2(cn) - abs2(d) - abs2(dn)
```
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

# Examples
```jldoctest
julia> cm=Float64[0, 0];JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2 3 4;5 6 7 8],[-1,1]);cm
2-element Vector{Float64}:
 13.0
 33.0
```
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

# Examples
```jldoctest
julia> @variables a b c d an bn cn dn;cm = Num[0, 0];JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn],[1, -1]);cm
2-element Vector{Num}:
 abs2(a) + abs2(an) - abs2(b) - abs2(bn)
 abs2(c) + abs2(cn) - abs2(d) - abs2(dn)
```
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
julia> JosephsonCircuits.calcqe([3/5 4/5;4/5 3/5],[0.0 0.0;0.0 0.0])
2×2 Matrix{Float64}:
 0.36  0.64
 0.64  0.36

julia> JosephsonCircuits.calcqe(Complex{Float64}[3/5 4/5;4/5 3/5],Complex{Float64}[0.0 0.0;0.0 0.0])
2×2 Matrix{Float64}:
 0.36  0.64
 0.64  0.36

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

# Examples
```jldoctest
julia> qe=Float64[1 2;3 4];JosephsonCircuits.calcqe!(qe,[1 2;3 4],[1 2 3;4 5 6]);qe
2×2 Matrix{Float64}:
 0.0526316  0.210526
 0.0882353  0.156863
```
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
scattering matrix. See also [`calcqeideal!`](@ref). 

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


"""
    calcqeideal!(qeideal,S)

See [`calcqeideal`](@ref). 

"""
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

"""
    calcinputcurrentoutputvoltage!(inputcurrent, outputvoltage, nodeflux,
        bnm, inputportindices, outputportindices, nodeindices, wmodes)

Calculate the elements of the Z matrix.

# Examples
```jldoctest
inputwave = JosephsonCircuits.LinearAlgebra.Diagonal(ComplexF64[0])
outputwave = ComplexF64[0;;]
bnm = ComplexF64[-1; 1;;]
portimpedanceindices = [2]
portimpedances = ComplexF64[50.0 + 0.0im]
nodeindices = [2 2 2 2 3; 3 3 1 1 1]
componenttypes = [:P, :R, :L, :C, :C]
wmodes = [1]
phin = ComplexF64[-50/(im*wmodes[1]);50/(im*wmodes[1]);;]
symfreqvar = nothing
JosephsonCircuits.calcinputcurrentoutputvoltage!(inputwave,outputwave,phin,bnm,portimpedanceindices,
    portimpedanceindices,nodeindices,wmodes)
println(inputwave)
println(outputwave)

# output
ComplexF64[-2.0 + 0.0im;;]
ComplexF64[-100.0 + 0.0im;;]
```
"""
function calcinputcurrentoutputvoltage!(inputcurrent, outputvoltage, nodeflux,
    bnm, inputportindices, outputportindices, nodeindices, wmodes)

    # check the size of inputwave

    # check the size of outputwave

    # check the sizes of all of the inputs

    # loop over input branches and modes to define inputwaves
    Ninputports = length(inputportindices)
    Noutputports = length(outputportindices)
    Nsolutions = size(nodeflux,2)
    Nmodes = length(wmodes)

    for i in 1:Ninputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

                sourcecurrent = calcsourcecurrent(
                    nodeindices[1,inputportindices[i]],
                    nodeindices[2,inputportindices[i]],
                    bnm,Nmodes,j,k)

                # this will give NaN for DC, so set kval=0 in that case
                # if wmodes[j] == 0
                    # inputcurrent[(i-1)*Nmodes+j,k] = 0
                # else
                inputcurrent[(i-1)*Nmodes+j,k] = sourcecurrent/sqrt(abs(wmodes[j]))
                # end
            end
        end
    end

    # loop over output branches and modes to define outputwaves
    for i in 1:Noutputports
        for j in 1:Nmodes
            for k in 1:Nsolutions

                portvoltage = calcportvoltage(
                    nodeindices[1,outputportindices[i]],
                    nodeindices[2,outputportindices[i]],
                    nodeflux,
                    wmodes,
                    Nmodes,j,k)

                # convert from sqrt(power) to sqrt(photons/second)
                # this will give NaN for DC, so set kval=0 in that case
                # if wmodes[j] == 0
                    # outputvoltage[(i-1)*Nmodes+j,k] = 0
                # else
                outputvoltage[(i-1)*Nmodes+j,k] = portvoltage/sqrt(abs(wmodes[j]))
                # end
            end
        end
    end
 
    return nothing
end