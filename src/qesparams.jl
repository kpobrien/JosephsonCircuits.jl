
"""
    calcinputoutput!(inputwave, outputwave, phin, bnm, inputportindices,
        outputportindices, inputportimpedances, outputportimpedances,
        nodeindices, componenttypes, wmodes, symfreqvar)

Return the input and output waves for the system linearized around the strong
pump.

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
println(outputwave)

# output
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
println(outputwave)

# output
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
println(outputwave)

# output
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
println(outputwave)

# output
ComplexF64[3.5355339059327378 + 0.0im;;]
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
println(outputwave)

# output
ComplexF64[-10.606601717798213 + 0.0im;;]
```
"""
function calcinputoutput!(inputwave, outputwave, phin, bnm, inputportindices,
    outputportindices, inputportimpedances, outputportimpedances,
    nodeindices, componenttypes, wmodes, symfreqvar)
    return calcinputoutput_inner!(inputwave, outputwave, phin, bnm,
        inputportindices, outputportindices, inputportimpedances,
        outputportimpedances, nodeindices, componenttypes, wmodes, symfreqvar,
        false)
end

"""
    calcinputoutputnoise!(S, inputwave, outputwave, phin, bnm,
        inputportindices, outputportindices, inputportimpedances,
        outputportimpedances, nodeindices, componenttypes, wmodes, symfreqvar)

Return the input and output waves for the system linearized around the strong
pump.

This is a bit of a hack but I ran into issues with complex capacitance when
the capacitor was at the same branch as a current source. The calcS function
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
println(noiseoutputwave)

# output
ComplexF64[-5.568327974762547e-11 + 3.516177070001411e-15im;;]
```
"""
function calcinputoutputnoise!(inputwave, outputwave, phin, bnm,
    inputportindices, outputportindices, inputportimpedances,
    outputportimpedances, nodeindices, componenttypes, wmodes, symfreqvar)
    return calcinputoutput_inner!(inputwave, outputwave, phin, bnm,
        inputportindices, outputportindices, inputportimpedances,
        outputportimpedances, nodeindices, componenttypes, wmodes, symfreqvar,
        true)
end

"""
    thermaloccupation(w, temperature)

The factor `2*nbar + 1 = coth(hbar*abs(w)/(2*k*T))` by which a mode at
angular frequency `w` in thermal equilibrium at `temperature` in Kelvin
exceeds its vacuum noise, and `1` at zero temperature or zero frequency.

This is the whole of the temperature dependence of the noise. A dissipative
element adds noise whose covariance is its vacuum covariance times this,
which in the Gaussian channel picture scales the `Y` of the map without
touching its `X`.

Zero frequency returns one rather than diverging: the wave normalization of
a noise channel is already zero there (see [`portwavescale`](@ref)), so the
factor multiplies nothing.

# Examples
```jldoctest
julia> JosephsonCircuits.thermaloccupation(2*pi*5e9, 0.0)
1.0

julia> round(JosephsonCircuits.thermaloccupation(2*pi*5e9, 0.1), digits=5)
1.19962

julia> round(JosephsonCircuits.thermaloccupation(2*pi*5e9, 1.0), digits=4)
8.3746
```
"""
function thermaloccupation(w, temperature)
    (iszero(temperature) || iszero(w)) && return one(Float64)
    x = reduced_planck_constant*abs(w)/(2*boltzmann_constant*temperature)
    return coth(x)
end

"""
    calcnoisecovariance!(Cnoise, Snoise, occupation)

The added noise covariance at the output ports,

    Cnoise[i,i'] = sum_c occupation[c] Snoise[c,i] conj(Snoise[c,i'])

which is the `Y` of the Gaussian channel whose `X` is the scattering matrix:
the map takes an input covariance to `X sigma X' + Y`. `Snoise` describes the
transformation and `occupation` the state of each channel, so this is where
the two meet and where a temperature shows up.

In the normalization of the rest of these outputs a vacuum channel counts as
one, so at zero temperature the diagonal is `sum(abs2, Snoise[:,i])`, which
is exactly the noise term in the denominator of [`calcqe!`](@ref).
"""
function calcnoisecovariance!(Cnoise::AbstractMatrix, Snoise::AbstractMatrix,
    occupation = nothing)
    np = size(Snoise, 2)
    if size(Cnoise) != (np, np)
        throw(DimensionMismatch(lazy"`Cnoise` has size $(size(Cnoise)) but the $(np) output port modes need ($(np), $(np))."))
    end
    fill!(Cnoise, zero(eltype(Cnoise)))
    @inbounds for c in axes(Snoise, 1)
        f = isnothing(occupation) ? 1.0 : occupation[c]
        iszero(f) && continue
        for ip in 1:np
            a = f*Snoise[c, ip]
            for iq in 1:np
                Cnoise[ip, iq] += a*conj(Snoise[c, iq])
            end
        end
    end
    return Cnoise
end

"""
    noiseoccupation!(occupation, temperatures, wmodes, Nmodes)

Fill `occupation` with `2*nbar + 1` for each row of the noise scattering
matrix, whose rows run over the noise channels with the modes innermost.

This is what a channel's vacuum share of the noise is multiplied by to get
the noise it actually carries. It is applied where the noise power is asked
for, by [`calcqe!`](@ref) and by the noise covariance, and never to
`Snoise` itself, which is a scattering matrix and does not depend on
temperature. The commutation relations, which are a statement about the
transformation rather than about the state of anything, therefore do not see
it at all.

`temperatures` is one temperature per noise channel, in the order
[`noisechannelnames`](@ref) gives them. `nothing`, and every temperature
being zero, both give all ones.
"""
noiseoccupation!(occupation::AbstractVector, ::Nothing, wmodes,
    Nmodes::Integer) = fill!(occupation, 1.0)
function noiseoccupation!(occupation::AbstractVector, temperatures, wmodes,
    Nmodes::Integer)
    if length(occupation) != length(temperatures)*Nmodes
        throw(DimensionMismatch(lazy"`occupation` has length $(length(occupation)) but $(length(temperatures)) channels of $(Nmodes) modes need $(length(temperatures)*Nmodes)."))
    end
    @inbounds for c in eachindex(temperatures)
        for m in 1:Nmodes
            occupation[(c-1)*Nmodes + m] =
                thermaloccupation(wmodes[m], temperatures[c])
        end
    end
    return occupation
end

"""
    portwavescale(portimpedance, w)

The scale factor of the Kurokawa power waves at a port with impedance
`portimpedance` and (signed) mode frequency `w`, in units of
sqrt(photons/second) rather than sqrt(power):
`1/sqrt(real(Z))/sqrt(abs(w))`, and zero at zero frequency, where the wave
normalization is singular. This is the single definition used by the
scattering parameter calculation ([`calcinputoutput_inner!`](@ref)) and by
the sensitivity scaling ([`calcsensitivityscaling!`](@ref)), so the two
cannot drift apart.
"""
@inline function portwavescale(portimpedance, w)
    kval = 1/sqrt(Complex(real(portimpedance)))
    if w == 0
        return zero(kval)
    end
    return kval/sqrt(abs(w))
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
.

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
            # the port impedance and the wave scale depend on the port and
            # the mode, not on the drive column, so they are computed once
            # per (port, mode) rather than once per solution.
            portimpedance = calcimpedance(
                inputportimpedances[i],
                componenttypes[inputportindices[i]],
                wmodes[j],symfreqvar)
            kval = portwavescale(portimpedance, wmodes[j])
            for k in 1:Nsolutions

                sourcecurrent = calcsourcecurrent(
                    nodeindices[1,inputportindices[i]],
                    nodeindices[2,inputportindices[i]],
                    bnm,Nmodes,j,k)

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
            # the port impedance and the wave scale depend on the port and
            # the mode, not on the drive column, so they are computed once
            # per (port, mode) rather than once per solution.
            portimpedance = calcimpedance(
                outputportimpedances[i],
                componenttypes[outputportindices[i]],
                wmodes[j],symfreqvar)
            kval = portwavescale(portimpedance, wmodes[j])
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

                # calculate the current flowing through the port
                if nosource
                    portcurrent = - portvoltage / portimpedance
                else
                    portcurrent = sourcecurrent - portvoltage / portimpedance
                end

                # calculate the input and output power waves as defined in
                # (except in units of sqrt(photons/second) instead of
                # sqrt(power)
                # K. Kurokawa, "Power Waves and the Scattering Matrix", IEEE
                # Trans. Micr. Theory and Tech. 13, 194–202 (1965) 
                # doi: 10.1109/TMTT.1965.1125964
                outputwave[(i-1)*Nmodes+j,k] = 1/2*kval * (portvoltage - conj(portimpedance) * portcurrent)
            end
        end
    end
 
    return nothing
end

"""
    calcscatteringmatrix!(S, inputwave::Diagonal, outputwave)

The scattering matrix is defined as `outputwave = S * inputwave`.

# Examples
```jldoctest
julia> inputwave=JosephsonCircuits.LinearAlgebra.Diagonal([1.0,1.0]);outputwave=[im/sqrt(2) 1/sqrt(2);1/sqrt(2) im/sqrt(2)];S = zeros(Complex{Float64},2,2);JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave);S
2×2 Matrix{ComplexF64}:
      0.0+0.707107im  0.707107+0.0im
 0.707107+0.0im            0.0+0.707107im
```
"""
function calcscatteringmatrix!(S, inputwave::Diagonal, outputwave)
    # copy!(S,outputwave)
    # rdiv!(S,inputwave)
    rdiv!(outputwave,inputwave)
    copy!(S,outputwave)
    
    return nothing
end

"""
    calcscatteringmatrix!(S, inputwave, outputwave)

The scattering matrix is defined as `outputwave = S * inputwave`.

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
function calcscatteringmatrix!(S, inputwave, outputwave)
    S .= outputwave / inputwave
    return nothing
end

"""
    calcscatteringmatrix!(S, inputwave::Vector, outputwave::Vector)

The scattering matrix is defined as `outputwave = S * inputwave`.

# Examples
```jldoctest
julia> inputwave=[1.0,0.0];outputwave=[im/sqrt(2), 1/sqrt(2)];S = zeros(Complex{Float64},2,2);JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave);S
2×2 Matrix{ComplexF64}:
      0.0+0.707107im  0.0+0.0im
 0.707107+0.0im       0.0+0.0im
```
"""
function calcscatteringmatrix!(S, inputwave::Vector, outputwave::Vector)
    if size(S,1) != length(outputwave)
        throw(DimensionMismatch(lazy"First dimension of scattering matrix not consistent with first dimensions of outputwave."))
    end
    if size(S,2) != length(inputwave)
        throw(DimensionMismatch(lazy"Second dimension of scattering matrix not consistent with first dimension of input wave."))
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

"""
    calcportvoltage(key1, key2, phin, wmodes, Nmodes, j, k)

"""
function calcportvoltage(key1, key2, phin, wmodes, Nmodes, j, k)

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
    # portvoltage *= im*abs(wmodes[j])
    portvoltage *= im*wmodes[j]

    return portvoltage
end

"""
    calcsourcecurrent(key1, key2, bnm, Nmodes, j, k)

"""
function calcsourcecurrent(key1, key2, bnm, Nmodes, j, k)

    if key1 == 1
        sourcecurrent = -bnm[(key2-2)*Nmodes+j,k]
    elseif key2 == 1
        sourcecurrent =  bnm[(key1-2)*Nmodes+j,k]
    else
        sourcecurrent =  bnm[(key1-2)*Nmodes+j,k] 
        sourcecurrent -= bnm[(key2-2)*Nmodes+j,k]
        sourcecurrent /= 2
    end
    return sourcecurrent
end

# """
#     IMPEDANCE_R, IMPEDANCE_C, IMPEDANCE_L

# The component types [`impedance`](@ref) knows, as integers.

# A kernel cannot carry a `Symbol`, so the shared impedance is written against
# these and [`impedancecode`](@ref) maps a component type to one at the boundary.
# """
const IMPEDANCE_R, IMPEDANCE_C, IMPEDANCE_L = Int32(1), Int32(2), Int32(3)

"""
    impedancecode(type)

The [`impedance`](@ref) code of a component type, or an error for a type which
has no impedance.
"""
function impedancecode(type)
    type === :R && return IMPEDANCE_R
    type === :C && return IMPEDANCE_C
    type === :L && return IMPEDANCE_L
    error(lazy"Unknown component type")
end

"""
    impedance(c, code::Integer, w)

The impedance of a component of value `c` and type `code` at frequency `w`,
conjugating the stored value at a negative frequency.

The single implementation, called by [`calcimpedance`](@ref) on the host and
directly from the kernels which compute power waves on a backend.
"""
@inline function impedance(c, code::Integer, w)
    cc = real(w) >= 0 ? c : conj(c)
    code == IMPEDANCE_R && return cc + 0.0im
    code == IMPEDANCE_C && return 1/(im*w*cc)
    return im*w*cc
end

"""
    calcimpedance(c::Union{Integer,T,Complex{T}}, type, w, symfreqvar,
        ) where {T<:AbstractFloat}

# Examples
```jldoctest
julia> JosephsonCircuits.calcimpedance(30.0,:C,1.0,nothing)
0.0 - 0.03333333333333333im

julia> JosephsonCircuits.calcimpedance(30.0,:L,1.0,nothing)
0.0 + 30.0im

julia> JosephsonCircuits.calcimpedance(30.0,:R,1.0,nothing)
30.0 + 0.0im

julia> JosephsonCircuits.calcimpedance(30.0,:C,-1.0,nothing)
-0.0 + 0.03333333333333333im

julia> JosephsonCircuits.calcimpedance(30.0,:L,-1.0,nothing)
-0.0 - 30.0im

julia> JosephsonCircuits.calcimpedance(30.0,:R,-1.0,nothing)
30.0 + 0.0im

```
"""
function calcimpedance(c::Union{T,Complex{T}}, type, w, symfreqvar,
    ) where {T<:Union{AbstractFloat,Integer}}
    return impedance(c, impedancecode(type), w)
end


"""
    calcimpedance(c, type, w, symfreqvar)

# Examples
```jldoctest
julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:R,2.0,w)
60.0 + 0.0im

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:C,2.0,w)
0.0 - 0.008333333333333333im

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:L,2.0,w)
0.0 + 120.0im

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:R,-2.0,w)
-60.0 + 0.0im

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:C,-2.0,w)
0.0 - 0.008333333333333333im

julia> @variables w;JosephsonCircuits.calcimpedance(30*w,:L,-2.0,w)
0.0 + 120.0im
```
"""
function calcimpedance(c, type, w, symfreqvar)
    # substitutefreq evaluates FrequencyDependent provider leaves at the
    # signed mode frequency whether or not a symbolic frequency variable
    # is in use, and substitutes symfreqvar when one is; on a plain
    # number it is the identity
    v = substitutefreq(c, symfreqvar, w)
    if type == :R
        if w >= 0
            return v+0.0im
        else
            return conj(v)+0.0im
        end
    elseif type == :C
        if w >= 0
            return 1/(im*w*v)
        else
            return 1/(im*w*conj(v))
        end
    elseif type == :L
        if w >= 0
            return (im*w*v)
        else
            return (im*w*conj(v))
        end
    else
        error(lazy"Unknown component type")
    end
end

"""
    calccm(S::AbstractArray, w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Sum the abs2 of each element along the horizontal
axis, applying a minus sign if the corresponding frequency is negative.
Represents energy conservation.

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
```
"""
function calccm(S::AbstractArray{T}, w) where {T}
    cm = zeros(T,size(S,1))
    return calccm!(cm,S,w)
end

function calccm(S::AbstractArray{Complex{T}}, w) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    cm = zeros(T,size(S,1))
    return calccm!(cm,S,w)
end

"""
    calccm!(cm::AbstractArray{T}, S, w) where {T<:AbstractFloat}

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
function calccm!(cm::AbstractArray{T}, S, w) where {T<:AbstractFloat}

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch(lazy"Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch(lazy"First dimension of scattering matrix must equal the length of cm."))        
    end

    # use a Kahan, Babushka, Neumaier compensated sum. more cache efficient
    # version
    fill!(cm,zero(T))
    c = zeros(eltype(cm),size(cm))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            t = cm[i] + abs2(S[i,j])*sign(real(w[(j-1) % m + 1]))
            c[i] += ifelse(
                abs(cm[i]) >= abs2(S[i,j]),
                (cm[i]-t) + abs2(S[i,j])*sign(real(w[(j-1) % m + 1])),
                (abs2(S[i,j])*sign(real(w[(j-1) % m + 1]))-t) + cm[i])
            cm[i]  = t
        end
    end

    @inbounds for i in 1:size(S,1)
        cm[i]+=c[i]
    end

    return cm
end

"""
    calccm!(cm, S, w)

Calculate the bosonic commutation relations for a scattering matrix S in the
field ladder operator basis. Overwrites cm with output.

# # Examples
# ```jldoctest
# julia> @variables a b;cm=Num[0,0];JosephsonCircuits.calccm!(cm,[a b; b a],[-1,1]);cm
# 2-element Vector{Num}:
#  -abs2(a) + abs2(b)
#   abs2(a) - abs2(b)
# ```
"""
function calccm!(cm, S, w)

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch(lazy"Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch(lazy"First dimension of scattering matrix must equal the length of cm."))        
    end

    # more cache friendly version
    fill!(cm,zero(eltype(cm)))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            cm[i] += abs2(S[i,j])*sign(w[(j-1) % m + 1])
        end
    end

    return cm
end

"""
    calccm(S::AbstractArray, Snoise::AbstractArray, w)

Calculate the bosonic commutation relations for a scattering matrix S in the 
field ladder operator basis. Sum the abs2 of each element along the horizontal
axis, applying a minus sign if the corresponding frequency is negative.
Represents energy conservation.

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
```
"""
function calccm(S::AbstractArray{T}, Snoise::AbstractArray{T}, w) where {T}
    cm = zeros(T,size(S,1))
    return calccm!(cm,S,Snoise,w)
end

function calccm(S::AbstractArray{Complex{T}}, Snoise::AbstractArray{Complex{T}},
    w) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    cm = zeros(T,size(S,1))
    return calccm!(cm,S,Snoise,w)
end


"""
    calccm!(cm::AbstractArray{T}, S, Snoise, w) where {T<:AbstractFloat}

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
function calccm!(cm::AbstractArray{T}, S, Snoise, w) where {T<:AbstractFloat}

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch(lazy"Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    for d in size(Snoise)
        if mod(d, m) != 0
            throw(DimensionMismatch(lazy"Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch(lazy"First dimension of scattering matrix must equal the length of cm."))        
    end

    if size(S,1) != size(Snoise,1)
        throw(DimensionMismatch(lazy"First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
    end

    # use a Kahan, Babushka, Neumaier compensated sum. more cache efficient
    # version
    fill!(cm,zero(T))
    c = zeros(eltype(cm),size(cm))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            t = cm[i] + abs2(S[i,j])*sign(real(w[(j-1) % m + 1]))
            c[i] += ifelse(
                abs(cm[i]) >= abs2(S[i,j]),
                (cm[i]-t) + abs2(S[i,j])*sign(real(w[(j-1) % m + 1])),
                (abs2(S[i,j])*sign(real(w[(j-1) % m + 1]))-t) + cm[i])
            cm[i]  = t
        end
    end

    @inbounds for j in 1:size(Snoise,2)
        for i in 1:size(Snoise,1)
            t = cm[i] + abs2(Snoise[i,j])*sign(real(w[(j-1) % m + 1]))
            c[i] += ifelse(
                abs(cm[i]) >= abs2(Snoise[i,j]),
                (cm[i]-t) + abs2(Snoise[i,j])*sign(real(w[(j-1) % m + 1])),
                (abs2(Snoise[i,j])*sign(real(w[(j-1) % m + 1]))-t) + cm[i])
            cm[i]  = t
        end
    end

    @inbounds for i in 1:size(S,1)
        cm[i]+=c[i]
    end

    return cm
end


"""
    calccm!(cm, S, Snoise, w)

Calculate the bosonic commutation relations for a scattering matrix S in the
field ladder operator basis. Overwrites cm with output.

# # Examples
# ```jldoctest
# julia> @variables a b c d an bn cn dn;cm = Num[0, 0];JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn],[1, -1]);cm
# 2-element Vector{Num}:
#  -abs2(bn) + abs2(a) - abs2(b) + abs2(an)
#   abs2(c) - abs2(dn) - abs2(d) + abs2(cn)
# ```
"""
function calccm!(cm, S, Snoise, w)

    m = length(w)

    for d in size(S)
        if mod(d, m) != 0
            throw(DimensionMismatch(lazy"Dimensions of scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    for d in size(Snoise)
        if mod(d, m) != 0
            throw(DimensionMismatch(lazy"Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."))
        end
    end

    if size(S,1) != length(cm)
        throw(DimensionMismatch(lazy"First dimension of scattering matrix must equal the length of cm."))        
    end

    if size(S,1) != size(Snoise,1)
        throw(DimensionMismatch(lazy"First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
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

    return cm
end

"""
    NoiseReduction

The noise scattering matrix reduced to what the quantum efficiency and the
commutation relations read of it.

Both consume `Snoise` only through a sum over the noise index: the quantum
efficiency through `sum(abs2, Snoise[:,i])` weighted by the occupation of
each channel, and the commutation relations through the same sum weighted by
the sign of each noise index's mode frequency instead. The occupation enters
the first and not the second, which is why they are separate sums and why
the commutation relations do not depend on temperature. On a circuit whose loss is spread along the line that is a
reduction of thousands of rows to one number per port mode, so when the noise
scattering parameters are not themselves an output there is no reason to form
the matrix on the host at all.

Passed to [`calcqe!`](@ref) and [`calccm!`](@ref) in place of the matrix.
"""
struct NoiseReduction{V}
    denom::V
    signed::V
end

"""
    calccm!(cm, S, noise::NoiseReduction, w)

The commutation relations from the scattering matrix and the reduction of the
noise scattering matrix, rather than the matrix itself.

`calccm!(cm, S, Snoise, w)` reads `Snoise` only as the sum of `abs2` over each
row weighted by the sign of the mode frequency of the column, so a caller
which has that sum already passes it here. See [`NoiseReduction`](@ref).

The scattering matrix's own contribution is left to `calccm!(cm, S, w)`, so
this shares that method's summation, compensated where the commutation
relations are real, rather than repeating it.
"""
function calccm!(cm, S, noise::NoiseReduction, w)
    return calccmreduced!(cm, S, noise, w)
end

# the same, for the real commutation relations, whose own method would
# otherwise be neither more nor less specific than the one above
function calccm!(cm::AbstractArray{T}, S, noise::NoiseReduction,
    w) where {T<:AbstractFloat}
    return calccmreduced!(cm, S, noise, w)
end

function calccmreduced!(cm, S, noise::NoiseReduction, w)

    if size(S,1) != length(noise.signed)
        throw(DimensionMismatch(lazy"First dimension of the scattering parameter matrix must equal the length of the noise reduction."))
    end

    calccm!(cm, S, w)
    @inbounds for i in eachindex(cm)
        cm[i] += noise.signed[i]
    end

    return cm
end


"""
    calcqe(S::AbstractArray)

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
```
"""
function calcqe(S::AbstractArray{T}) where {T}
    qe = zeros(T,size(S))
    return calcqe!(qe,S)
end

function calcqe(S::AbstractArray{Complex{T}}) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    qe = zeros(T,size(S))
    return calcqe!(qe,S)
end

"""
    calcqe!(qe, S)

Calculate the quantum efficiency matrix for a scattering matrix in the field
ladder operator basis. Overwrites qe with output.
"""
function calcqe!(qe, S)

    if size(qe) != size(S)
        throw(DimensionMismatch(lazy"Dimensions of quantum efficiency and scattering parameter matrices must be equal."))
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

    return qe
end

"""
    calcqe(S::AbstractArray, Snoise::AbstractArray)

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
```
"""
function calcqe(S::AbstractArray{T}, Snoise::AbstractArray{T}) where {T}
    qe = zeros(T,size(S))
    return calcqe!(qe,S,Snoise)
end

function calcqe(S::AbstractArray{Complex{T}},
    Snoise::AbstractArray{Complex{T}}) where {T}
    # commutation relations are real so if the type of complex, use this
    # parametric method to define a real matrix.
    qe = zeros(T,size(S))
    return calcqe!(qe,S,Snoise)
end

"""
    calcqe!(qe, S, Snoise)

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
function calcqe!(qe, S, Snoise, occupation = nothing)

    if size(qe) != size(S)
        throw(DimensionMismatch(lazy"Dimensions of quantum efficiency and scattering parameter matrices must be equal."))
    end

    if size(S,1) != size(Snoise,1)
        throw(DimensionMismatch(lazy"First dimensions of scattering parameter matrice and noise scattering matrix must be equal."))
    end

    # more cache efficient version of QE calculation
    denom = zeros(eltype(qe),size(S,1))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            denom[i] += abs2(S[i,j])
        end
    end

    # each noise channel contributes the noise it actually carries, which is
    # its vacuum share times the occupation of the channel. `Snoise` is a
    # scattering matrix and does not depend on temperature; the occupation
    # is applied here, where the noise power is what is being asked for.
    @inbounds for j in 1:size(Snoise,2)
        f = isnothing(occupation) ? one(eltype(denom)) : occupation[j]
        for i in 1:size(Snoise,1)
            denom[i] += f*abs2(Snoise[i,j])
        end
    end

    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            qe[i,j] = abs2(S[i,j]) / denom[i]
        end
    end

    return qe
end

"""
    calcqe!(qe, S, noise::NoiseReduction)

The quantum efficiency from the scattering matrix and the reduction of the
noise scattering matrix, rather than the matrix itself.

`calcqe!(qe, S, Snoise)` reads `Snoise` only as `sum(abs2, Snoise[i,:])` per
row, so a caller which has that sum already, as one computed on a backend has,
passes it here and never forms the matrix. See [`NoiseReduction`](@ref).
"""
function calcqe!(qe, S, noise::NoiseReduction, occupation)
    # the reduction's `denom` already carries the occupation, applied on the
    # backend where it was summed, so it is not applied again here
    return calcqe!(qe, S, noise)
end

function calcqe!(qe, S, noise::NoiseReduction)

    if size(qe) != size(S)
        throw(DimensionMismatch(lazy"Dimensions of quantum efficiency and scattering parameter matrices must be equal."))
    end

    if size(S,1) != length(noise.denom)
        throw(DimensionMismatch(lazy"First dimension of the scattering parameter matrix must equal the length of the noise reduction."))
    end

    denom = zeros(eltype(qe),size(S,1))
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            denom[i] += abs2(S[i,j])
        end
    end

    @inbounds for i in 1:size(S,1)
        denom[i] += noise.denom[i]
    end

    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            qe[i,j] = abs2(S[i,j]) / denom[i]
        end
    end
    return qe
end


"""
    calcqeideal(S::AbstractArray)

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
    return calcqeideal!(qeideal,S)
end

function calcqeideal(S::AbstractArray{Complex{T}}) where {T}
    # quantum efficiency is real so if the type of complex, use this
    # parametric method to define a real matrix.
    qeideal = zeros(T,size(S))
    return calcqeideal!(qeideal,S)
end


"""
    calcqeideal!(qeideal,S)

See [`calcqeideal`](@ref).

"""
function calcqeideal!(qeideal,S)
    if size(qeideal) != size(S)
        throw(DimensionMismatch(lazy"Sizes of QE and S matrices must be equal."))
    end
    for i in eachindex(S)
        abs2S = abs2(S[i])
        qeideal[i] = ifelse(abs2S <= 1,one(eltype(qeideal)),1 /(2 - 1 /abs2S))
    end
    return qeideal
end

"""
    calcCnoise(S::AbstractMatrix{T}) where {T}

Return the noise wave covariance matrix computed using Bosma's theorem for a
passive linear network with scattering parameter matrix `S`. The network can
be lossy and non-reciprocal.

This function assumes vacuum fluctuations as input to the ports, but could be
extended to allow arbitrary noise temperatures.

# Examples
```jldoctest
julia> JosephsonCircuits.calcCnoise(JosephsonCircuits.S_splitter!(zeros(Complex{Float64},2,2)))
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im

julia> S = rand(Complex{Float64},3,3);isapprox(JosephsonCircuits.calcCnoise(S),[1.0 0 0;0 1 0;0 0 1].-S*S')
true
```
"""
function calcCnoise(S::AbstractMatrix{T}) where {T}
    Cnoise = zeros(T,size(S))
    return calcCnoise!(Cnoise,S)
end

"""
    calcCnoise!(Cnoise, S)

Calculate the noise wave covariance matrix for a scattering matrix in the field
ladder operator basis. Overwrites `Cnoise` with output.
"""
function calcCnoise!(Cnoise::AbstractMatrix, S)

    if size(Cnoise) != size(S)
        throw(DimensionMismatch(lazy"The dimensions of the noise wave covariance and scattering parameter matrices must be equal."))
    end

    if size(S,1) != size(S,2)
        throw(DimensionMismatch(lazy"The scattering parameter and noise wave covariance matrices must be square."))
    end

    # compute C = I - S*S' from Bosma's theorem
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            # if i == j
            #     Cnoise[i,j] = one(eltype(Cnoise))
            # else
            #     Cnoise[i,j] = zero(eltype(Cnoise))
            # end
            Cnoise[i,j] = zero(eltype(Cnoise))
            for k in 1:size(S,2)
                # use abs2 as a cludge to make sure QE is identical for
                # symbolic math with real variables.
                Cnoise[i,j] -= ifelse(i==j,abs2(S[i,k]),S[i,k]*conj(S[j,k]))
                # Cnoise[i,j] -= S[i,k]*conj(S[j,k])
            end
        end
    end

    @inbounds for k in 1:size(S,2)
        Cnoise[k,k] += one(eltype(Cnoise))
    end

    return Cnoise
end

"""
    calcCnoise(S::AbstractArray{T}, Snoise::AbstractArray{T}) where {T}

Calculate the noise wave covariance matrix for a scattering matrix in the
field ladder operator basis.

# Examples
```jldoctest
julia> JosephsonCircuits.calcCnoise([3/5 4/5;4/5 3/5],[0.0 0.0;0.0 0.0])
2×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0

julia> JosephsonCircuits.calcCnoise(Complex{Float64}[3/5 4/5;4/5 3/5],Complex{Float64}[0.0 0.0;0.0 0.0])
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  0.0+0.0im
```
"""
function calcCnoise(S::AbstractArray{T}, Snoise::AbstractArray{T}) where {T}
    Cnoise = zeros(T,size(S))
    return calcCnoise!(Cnoise,S,Snoise)
end

"""
    calcCnoise!(Cnoise, S, Snoise)

Calculate the noise wave covariance matrix for a scattering matrix in the
field ladder operator basis. Overwrites `Cnoise` with output.

# Examples
```jldoctest
julia> C=zeros(Float64,2,2);JosephsonCircuits.calcCnoise!(C,[1 2;3 4],[0.0 0 0;0 0 0]);C
2×2 Matrix{Float64}:
 0.0  0.0
 0.0  0.0
```
"""
function calcCnoise!(Cnoise, S, Snoise)

    if size(Cnoise) != size(S)
        throw(DimensionMismatch(lazy"The dimensions of the noise wave covariance and scattering parameter matrices must be equal."))
    end

    if size(S,1) != size(Snoise,1)
        throw(DimensionMismatch(lazy"The first dimensions of the scattering parameter and noise scattering parameter matrices must be equal."))
    end

    # add the noise covariance from the noise ports to the
    # physical ports
    @inbounds for j in 1:size(S,2)
        for i in 1:size(S,1)
            Cnoise[i,j] = zero(eltype(Cnoise))
            for k in 1:size(Snoise,2)
                # use abs2 as a cludge to make sure QE is identical for
                # symbolic math with real variables.
                Cnoise[i,j] += ifelse(i==j,abs2(Snoise[i,k]),Snoise[i,k]*conj(Snoise[j,k]))
                # Cnoise[i,j] -= Snoise[i,k]*conj(Snoise[j,k])

            end
        end
    end

    return Cnoise
end

"""
    calcqe_S_Cnoise(S::AbstractArray, Cnoise::AbstractArray)

Calculate the noise wave covariance matrix from the scattering parameter
matrix and the noise covariance matrix, both in the field ladder operator
(sqrt photon number) basis.

# Examples
```jldoctest
julia> S = JosephsonCircuits.ABCDtoS(JosephsonCircuits.ABCD_attenuator_T(50,10));isapprox(JosephsonCircuits.calcqe_S_Cnoise(S,JosephsonCircuits.calcCnoise(S)),[0 0.1;0.1 0])
true

julia> S = JosephsonCircuits.ABCDtoS(JosephsonCircuits.ABCD_attenuator_T(50,10).+0im);isapprox(JosephsonCircuits.calcqe_S_Cnoise(S,JosephsonCircuits.calcCnoise(S)),[0 0.1;0.1 0])
true
```
"""
function calcqe_S_Cnoise(S::AbstractArray{T}, Cnoise::AbstractArray{T}) where {T}
    qe = zeros(T,size(S))
    return calcqe_S_Cnoise!(qe,S,Cnoise)
end

function calcqe_S_Cnoise(S::AbstractArray{Complex{T}},
    Cnoise::AbstractArray{Complex{T}}) where {T}
    # QE is real so if the type is complex, use this
    # parametric method to define a real matrix.
    qe = zeros(T,size(S))
    return calcqe_S_Cnoise!(qe,S,Cnoise)
end

"""
    calcqe_S_Cnoise!(qe, S, Cnoise)

Calculate the quantum efficiency matrix from the scattering parameter matrix
and the noise wave covariance matrix, both in the field ladder operator (sqrt
photon number) basis. Overwrites qe with output.

"""
function calcqe_S_Cnoise!(qe, S, Cnoise)

    if size(qe) != size(S)
        throw(DimensionMismatch(lazy"The dimensions of the quantum efficiency and scattering parameter matrices must be equal."))
    end

    if size(S) != size(Cnoise)
        throw(DimensionMismatch(lazy"The dimensions of the noise wave covariance and scattering parameter matrices must be equal."))
    end

    # more cache efficient version of QE calculation
    denom = zeros(eltype(qe),size(Cnoise,1))
    @inbounds for i in 1:size(Cnoise,1)
            denom[i] += Cnoise[i,i]
    end

    # more cache efficient version of QE calculation
    # denom = zeros(eltype(qe),size(S,1))
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

    return qe
end

