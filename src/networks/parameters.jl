# The network parameter conversions, S to Z, Z to S, S to T and the rest.
# Each is a kernel which converts one matrix into another,
# `f!(y, x, tmp, ...)`, defined further down beside its docstring, and the
# per frequency forms, on a matrix or on an array of matrices with one per
# frequency, allocating and in place, come from the driver here. A port
# impedance argument is a number, the same at every port and frequency, a
# vector with one value per port, or a matrix with one row per port and one
# column per frequency. The allocating forms return the element type the
# conversion needs, see [`conversiontype`](@ref), whatever that of the
# input.
#
# The scattering parameters are those of the waves
# `a_k = (V_k + Z_k I_k)/(2 sqrt(Z_k))` and `b_k = (V_k - Z_k I_k)/(2 sqrt(Z_k))`
# at each port `k` of impedance `Z_k`, with the principal square root. For a
# real `Z_k` these are the usual power waves; for a complex one they are
# pseudo-waves, which a load of impedance `Z_k` (not its conjugate)
# matches, and with which a lossless network need not have a unitary
# scattering matrix.

"""
    PortDiagonal(values, rows = Colon())

A port argument handed to a conversion kernel as the diagonal matrix of
its values at the frequency, restricted to the ports `rows`: `values` is a
vector with one value per port, the same at every frequency, or a matrix
with one row per port and one column per frequency. See
[`convertperfrequency!`](@ref).
"""
struct PortDiagonal{A<:AbstractArray,R}
    values::A
    rows::R
end
PortDiagonal(values::AbstractArray) = PortDiagonal(values, Colon())

"""
    atfrequency(a, i)

The port argument `a` at the frequency index `i`: a number is the same
everywhere, a vector holds one value per port for every frequency, a
matrix has one column per frequency, and a [`PortDiagonal`](@ref) is the
diagonal matrix of its ports there.
"""
atfrequency(a::Number, i) = a
atfrequency(a::AbstractVector, i) = a
atfrequency(a::AbstractMatrix, i) = view(a, :, i)
atfrequency(d::PortDiagonal{<:AbstractVector}, i) =
    Diagonal(view(d.values, d.rows))
atfrequency(d::PortDiagonal{<:AbstractMatrix}, i) =
    Diagonal(view(d.values, d.rows, i))

# the frequency axes a port argument carries, `nothing` when it is the same
# at every frequency
frequencyaxes(a::Number) = nothing
frequencyaxes(a::AbstractVector) = nothing
frequencyaxes(a::AbstractMatrix) = axes(a)[2:end]
frequencyaxes(d::PortDiagonal) = frequencyaxes(d.values)

"""
    convertperfrequency!(kernel!, y, x, args...)
    convertperfrequency(kernel!, x, args...)

Convert the matrix `x[:, :, i]` into `y[:, :, i]` at every frequency index
`i` of the trailing dimensions of `x` by the kernel
`kernel!(y_i, x_i, tmp, args_i...)`, with `tmp` a scratch matrix and each
port argument given at that frequency by [`atfrequency`](@ref). The
allocating form converts into an array like `x`, and for a matrix `x`,
one frequency, calls the kernel directly. There is one method per number
of port arguments, none, one or two, so that no argument map is compiled
per kernel and argument type. The type parameter makes the method
specialize on the kernel, which Julia does not do on its own for a
function argument.
"""
function convertperfrequency!(kernel!::F, y::AbstractArray,
        x::AbstractArray) where {F}
    tmp = checkconversion(y, x)
    for i in CartesianIndices(axes(x)[3:end])
        kernel!(view(y, :, :, i), view(x, :, :, i), tmp)
    end
    return y
end

function convertperfrequency!(kernel!::F, y::AbstractArray, x::AbstractArray,
        a) where {F}
    tmp = checkconversion(y, x, a)
    for i in CartesianIndices(axes(x)[3:end])
        kernel!(view(y, :, :, i), view(x, :, :, i), tmp, atfrequency(a, i))
    end
    return y
end

function convertperfrequency!(kernel!::F, y::AbstractArray, x::AbstractArray,
        a, b) where {F}
    tmp = checkconversion(y, x, a, b)
    for i in CartesianIndices(axes(x)[3:end])
        kernel!(view(y, :, :, i), view(x, :, :, i), tmp, atfrequency(a, i),
            atfrequency(b, i))
    end
    return y
end

convertperfrequency(kernel!::F, x::AbstractArray, args...) where {F} =
    convertperfrequency!(kernel!, similar(x, conversiontype(x, args...)), x,
        args...)

# a matrix is one frequency and is converted by the kernel directly: no
# loop, no views and no scratch of its own to compile per kernel and
# argument type, which for the single matrix calls of user code is nearly
# all of what the loop form would compile
function convertperfrequency(kernel!::F, x::AbstractMatrix) where {F}
    T = conversiontype(x)
    y = similar(x, T)
    kernel!(y, x, similar(x, T))
    return y
end
function convertperfrequency(kernel!::F, x::AbstractMatrix, a) where {F}
    checkconversion(x, x, a)
    T = conversiontype(x, a)
    y = similar(x, T)
    kernel!(y, x, similar(x, T), atmatrix(a))
    return y
end
function convertperfrequency(kernel!::F, x::AbstractMatrix, a, b) where {F}
    checkconversion(x, x, a, b)
    T = conversiontype(x, a, b)
    y = similar(x, T)
    kernel!(y, x, similar(x, T), atmatrix(a), atmatrix(b))
    return y
end

# the port argument of a single matrix: as `atfrequency` but with the
# diagonal built from the vector itself, not a view of it
atmatrix(a::Number) = a
atmatrix(a::AbstractVector) = a
atmatrix(d::PortDiagonal{<:AbstractVector}) =
    Diagonal(d.rows isa Colon ? d.values : d.values[d.rows])

"""
    convertcopy(kernel!, x, a)

The two port chain conversions in place on one copy of `x`: the kernel
`kernel!(y_i, a_i)` converts the matrix at every frequency index of the
copy, with the port argument given at that frequency by
[`atfrequency`](@ref); a matrix `x` is converted directly.
"""
function convertcopy(kernel!::F, x::AbstractArray, a) where {F}
    checkconversion(x, x, a)
    y = copyto!(similar(x, conversiontype(x, a)), x)
    for i in CartesianIndices(axes(x)[3:end])
        kernel!(view(y, :, :, i), atfrequency(a, i))
    end
    return y
end
function convertcopy(kernel!::F, x::AbstractMatrix, a) where {F}
    checkconversion(x, x, a)
    y = copyto!(similar(x, conversiontype(x, a)), x)
    kernel!(y, atmatrix(a))
    return y
end

"""
    conversiontype(x, args...)

The element type of the conversion of the array `x` with the port
arguments `args`: the kernels divide, so it is a floating point type (or
whatever type the division of `x`'s element type gives), and it is
complex when `x` or a port argument is.
"""
conversiontype(x::AbstractArray, args...) =
    promote_type(typeof(one(eltype(x)) / one(eltype(x))), map(porttype, args)...)
porttype(a::Number) = typeof(a)
porttype(a::AbstractArray) = eltype(a)
porttype(d::PortDiagonal) = eltype(d.values)

# the checks of the driver, and its scratch matrix, of the element type of
# the output: the output like the input, and every port argument with one
# column per frequency of the input or none
function checkconversion(y, x, args...)
    axes(y) == axes(x) || throw(DimensionMismatch(
        lazy"Sizes of output $(size(y)) and input $(size(x)) must be equal."))
    trailing = axes(x)[3:end]
    for a in args
        fa = frequencyaxes(a)
        (isnothing(fa) || fa == trailing) || throw(ArgumentError(
            lazy"A port argument of size $(size(a isa PortDiagonal ? a.values : a)) does not have one column per frequency of the input of size $(size(x))."))
    end
    return similar(x, eltype(y), axes(x)[1:2])
end

"""
    portdiagonal(a)
    porthalves(a)

A port impedance argument as the kernels take it: a number stays a
number and an array becomes a [`PortDiagonal`](@ref); `porthalves` splits
the ports in two, the first half the input ports of a chain matrix and the
second half its output ports.
"""
portdiagonal(a::Number) = a
portdiagonal(a::AbstractArray) = PortDiagonal(a)
porthalves(a::Number) = (a, a)
function porthalves(a::AbstractArray)
    n = size(a, 1)
    return PortDiagonal(a, 1:n÷2), PortDiagonal(a, n÷2+1:n)
end

# the factor of column `j` of a port argument of a kernel: a number is the
# same for every port, a diagonal holds one per port
portfactor(g::Number, j) = g
portfactor(g::Diagonal, j) = g.diag[j]

# the index ranges of the first and second halves of the ports of a matrix
# whose conversion splits its ports in two: the input and output ports of
# a chain or transmission matrix, which is square with an even number of
# rows
function porthalfranges(A::AbstractMatrix)
    n = size(A, 1)
    if size(A, 2) != n || isodd(n)
        throw(DimensionMismatch(lazy"A matrix whose ports are split into inputs and outputs must be square with an even number of rows, not of size $(size(A))."))
    end
    return 1:n÷2, n÷2+1:n
end

# the chain matrix and the scattering matrix of a two port are 2 by 2
function checktwoport(A::AbstractMatrix)
    if size(A) != (2, 2)
        throw(DimensionMismatch(lazy"The matrix of a two port is 2 by 2, not of size $(size(A))."))
    end
    return nothing
end

# the per frequency forms: allocating, and in place through a copy

StoT(x::AbstractArray) = convertperfrequency(StoT!, x)
TtoS(x::AbstractArray) = convertperfrequency(TtoS!, x)
AtoB(x::AbstractArray) = convertperfrequency(AtoB!, x)
BtoA(x::AbstractArray) = convertperfrequency(BtoA!, x)
AtoZ(x::AbstractArray) = convertperfrequency(AtoZ!, x)
ZtoA(x::AbstractArray) = convertperfrequency(ZtoA!, x)
AtoY(x::AbstractArray) = convertperfrequency(AtoY!, x)
YtoA(x::AbstractArray) = convertperfrequency(YtoA!, x)
BtoY(x::AbstractArray) = convertperfrequency(BtoY!, x)
YtoB(x::AbstractArray) = convertperfrequency(YtoB!, x)
BtoZ(x::AbstractArray) = convertperfrequency(BtoZ!, x)
ZtoB(x::AbstractArray) = convertperfrequency(ZtoB!, x)
ZtoY(x::AbstractArray) = convertperfrequency(ZtoY!, x)
YtoZ(x::AbstractArray) = convertperfrequency(YtoZ!, x)
StoT!(x::AbstractArray) = copy!(x, StoT(x))
TtoS!(x::AbstractArray) = copy!(x, TtoS(x))
AtoB!(x::AbstractArray) = copy!(x, AtoB(x))
BtoA!(x::AbstractArray) = copy!(x, BtoA(x))
AtoZ!(x::AbstractArray) = copy!(x, AtoZ(x))
ZtoA!(x::AbstractArray) = copy!(x, ZtoA(x))
AtoY!(x::AbstractArray) = copy!(x, AtoY(x))
YtoA!(x::AbstractArray) = copy!(x, YtoA(x))
BtoY!(x::AbstractArray) = copy!(x, BtoY(x))
YtoB!(x::AbstractArray) = copy!(x, YtoB(x))
BtoZ!(x::AbstractArray) = copy!(x, BtoZ(x))
ZtoB!(x::AbstractArray) = copy!(x, ZtoB(x))
ZtoY!(x::AbstractArray) = copy!(x, ZtoY(x))
YtoZ!(x::AbstractArray) = copy!(x, YtoZ(x))

# the chain matrix conversions of a two port work on a copy of the input
# with the port impedances as they are
ABCDtoS(x::AbstractArray; portimpedances = 50.0) =
    convertcopy(ABCDtoS!, x, portimpedances)
StoABCD(x::AbstractArray; portimpedances = 50.0) =
    convertcopy(StoABCD!, x, portimpedances)
ABCDtoS!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, ABCDtoS(x; portimpedances))
StoABCD!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, StoABCD(x; portimpedances))

# the scattering conversions take the square roots of the port impedances
StoZ(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(StoZ!, x, portdiagonal(sqrt.(portimpedances)))
YtoS(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(YtoS!, x, portdiagonal(sqrt.(portimpedances)))
StoY(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(StoY!, x, portdiagonal(1 ./ sqrt.(portimpedances)))
ZtoS(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(ZtoS!, x, portdiagonal(1 ./ sqrt.(portimpedances)))
StoZ!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, StoZ(x; portimpedances))
YtoS!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, YtoS(x; portimpedances))
StoY!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, StoY(x; portimpedances))
ZtoS!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, ZtoS(x; portimpedances))

# the chain conversions of scattering parameters take the square roots of
# the port impedances split between the input and the output ports
StoA(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(StoA!, x, porthalves(sqrt.(portimpedances))...)
AtoS(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(AtoS!, x, porthalves(sqrt.(portimpedances))...)
StoB(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(StoB!, x, porthalves(sqrt.(portimpedances))...)
BtoS(x::AbstractArray; portimpedances = 50.0) =
    convertperfrequency(BtoS!, x, porthalves(sqrt.(portimpedances))...)
StoA!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, StoA(x; portimpedances))
AtoS!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, AtoS(x; portimpedances))
StoB!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, StoB(x; portimpedances))
BtoS!(x::AbstractArray; portimpedances = 50.0) =
    copy!(x, BtoS(x; portimpedances))

@doc """
    StoT(S)

Convert the scattering parameter matrix `S` to a transmission matrix
`T` and return the result.

# Examples
```jldoctest
julia> S = Complex{Float64}[0.0 1.0;1.0 0.0];JosephsonCircuits.StoT(S)
2×2 Matrix{ComplexF64}:
 1.0+0.0im  -0.0+0.0im
 0.0-0.0im   1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoT


"""
    StoT!(T::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix)

See [`StoT`](@ref) for description.

"""
function StoT!(T::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix)

    range1, range2 = porthalfranges(S)

    # tmp = [-I S11; 0 S21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -one(eltype(tmp))
    end
    @views tmp[range1,range2] .= S[range1, range1]
    @views tmp[range2,range2] .= S[range2, range1]

    # T = [-S12 0; -S22 I]
    fill!(T,zero(eltype(T)))
    for d in range2
        T[d,d] = one(eltype(T))
    end

    @views T[range1,range1] .= .-S[range1, range2]
    @views T[range2,range1] .= .-S[range2, range2]

    # perform the left division
    # T = inv(tmp)*T = [-I S11; 0 S21] \ [-S12 0; -S22 I]
    ldiv!(lu!(tmp), T)

    return nothing
end


@doc """
    StoZ(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to an impedance parameter matrix
`Z` and return the result. Assumes a port impedance of 50 Ohms unless
specified with the `portimpedances` keyword argument, a scalar, vector, or
matrix of port impedances. A complex port impedance defines pseudo-waves,
as described in [`ZtoS`](@ref).

# Examples
```jldoctest
julia> S = Complex{Float64}[0.0 0.0;0.0 0.0];JosephsonCircuits.StoZ(S)
2×2 Matrix{ComplexF64}:
 50.0+0.0im   0.0+0.0im
  0.0+0.0im  50.0+0.0im

julia> S = Complex{Float64}[0.0 0.999;0.999 0.0];JosephsonCircuits.StoZ(S)
2×2 Matrix{ComplexF64}:
 49975.0+0.0im  49975.0+0.0im
 49975.0+0.0im  49975.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoZ

"""
    StoZ!(Z::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

See [`StoZ`](@ref) for description.

"""
function StoZ!(Z::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

    # tmp = (I - S)
    copy!(tmp,S)
    rmul!(tmp,-1)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end
    
    # Z = (I + S)*sqrt(portimpedances)
    copy!(Z,S)
    for d in 1:size(Z,1)
        Z[d,d] += 1
    end
    rmul!(Z,sqrtportimpedances)

    # perform the left division
    # Z = inv(tmp)*Z = (I - S) \ ((I + S)*sqrt(portimpedances))
    ldiv!(lu!(tmp), Z)

    # left multiply by sqrt(z)
    # compute sqrt(portimpedances)*((I - S) \ ((I + S)*sqrt(portimpedances)))
    lmul!(sqrtportimpedances,Z)
    return nothing
end

@doc """
    StoY(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to an admittance parameter matrix
`Y` and return the result. Assumes a port impedance of 50 Ohms unless
specified with the `portimpedances` keyword argument.

# Examples
```jldoctest
julia> S = Complex{Float64}[0.0 0.999;0.999 0.0];JosephsonCircuits.StoY(S)
2×2 Matrix{ComplexF64}:
  19.99+0.0im  -19.99+0.0im
 -19.99+0.0im   19.99+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoY

"""
    StoY!(Y::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)

In place version of [`StoY`](@ref), writing into `Y` with `tmp` as scratch
and `oneoversqrtportimpedances` the reciprocal square roots of the port
impedances.

"""
function StoY!(Y::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)
    
    # tmp = (I + S)
    copy!(tmp,S)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # Y = (I - S)*oneoversqrtportimpedances
    copy!(Y,S)
    rmul!(Y,-1)
    for d in 1:size(Y,1)
        Y[d,d] += 1
    end
    rmul!(Y,oneoversqrtportimpedances)

    # perform the left division
    # Y = inv(tmp)*Y = (I + S) \ ((I - S)*oneoversqrtportimpedances)
    ldiv!(lu!(tmp), Y)

    # left multiply by 1/sqrt(z)
    # compute oneoversqrtportimpedances*((I + S) \ ((I - S)*oneoversqrtportimpedances))
    lmul!(oneoversqrtportimpedances,Y)
    return nothing
end


@doc """
    StoA(S; portimpedances = 50.0)

Convert the scattering parameter matrix `S` to the chain (ABCD) matrix `A` and
return the result. The first half of the ports are the inputs of the chain
matrix and the second half its outputs, each with its own port impedances:
`portimpedances` is a scalar, a vector with one value per port, or a matrix
with one row per port and one column per frequency, 50 Ohms unless specified.

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoA

"""
    StoA!(A::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`StoA`](@ref) for description.

"""
function StoA!(A::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1, range2 = porthalfranges(S)
    h = length(range1)

    # tmp = [(S11-I)/g1 (I+S11)*g1; S21/g1 S21*g1]
    # A = [-S12/g2 S12*g2; (I-S22)/g2 (I+S22)*g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2, each
    # scaling the columns of its block
    for j in 1:h
        g1 = portfactor(sqrtportimpedances1, j)
        g2 = portfactor(sqrtportimpedances2, j)
        for i in 1:h
            δ = i == j
            S11 = S[i, j]
            S12 = S[i, h+j]
            S21 = S[h+i, j]
            S22 = S[h+i, h+j]
            tmp[i, j] = (S11 - δ)/g1
            tmp[i, h+j] = (S11 + δ)*g1
            tmp[h+i, j] = S21/g1
            tmp[h+i, h+j] = S21*g1
            A[i, j] = -S12/g2
            A[i, h+j] = S12*g2
            A[h+i, j] = (δ - S22)/g2
            A[h+i, h+j] = (S22 + δ)*g2
        end
    end

    # perform the left division
    # A = inv(tmp)*A
    ldiv!(lu!(tmp), A)

    return nothing
end


@doc """
    StoB(S; portimpedances = 50.0)

Convert the scattering parameter matrix `S` to the inverse chain (ABCD) matrix
`B` and return the result. Note that despite the name, the inverse of the
chain matrix is not equal to the inverse chain matrix, inv(A) ≠ B. The first
half of the ports are the inputs of the chain matrix and the second half its
outputs, each with its own port impedances: `portimpedances` is a scalar, a
vector with one value per port, or a matrix with one row per port and one
column per frequency, 50 Ohms unless specified.

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoB

"""
    StoB!(B::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`StoB`](@ref) for description.

"""
function StoB!(B::AbstractMatrix, S::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1, range2 = porthalfranges(S)
    h = length(range1)

    # tmp = [S12/g2 S12*g2; (S22-I)/g2 (I+S22)*g2]
    # B = [(I-S11)/g1 (I+S11)*g1; -S21/g1 S21*g1]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2, each
    # scaling the columns of its block
    for j in 1:h
        g1 = portfactor(sqrtportimpedances1, j)
        g2 = portfactor(sqrtportimpedances2, j)
        for i in 1:h
            δ = i == j
            S11 = S[i, j]
            S12 = S[i, h+j]
            S21 = S[h+i, j]
            S22 = S[h+i, h+j]
            tmp[i, j] = S12/g2
            tmp[i, h+j] = S12*g2
            tmp[h+i, j] = (S22 - δ)/g2
            tmp[h+i, h+j] = (S22 + δ)*g2
            B[i, j] = (δ - S11)/g1
            B[i, h+j] = (S11 + δ)*g1
            B[h+i, j] = -S21/g1
            B[h+i, h+j] = S21*g1
        end
    end

    # perform the left division
    # B = inv(tmp)*B
    ldiv!(lu!(tmp), B)

    return nothing
end


@doc """
    TtoS(T)

Convert the transmission matrix `T` to a scattering parameter matrix `S` and return the result.

# Examples
```jldoctest
julia> T = Complex{Float64}[1.0 0.0;0.0 1.0];JosephsonCircuits.TtoS(T)
2×2 Matrix{ComplexF64}:
 0.0-0.0im  1.0+0.0im
 1.0+0.0im  0.0-0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of sign on T11 and T21 terms (suspected typo).
""" TtoS

"""
    TtoS!(S::AbstractMatrix,T::AbstractMatrix,tmp::AbstractMatrix)

See [`TtoS`](@ref) for description.

"""
function TtoS!(S::AbstractMatrix,T::AbstractMatrix,tmp::AbstractMatrix)

    range1, range2 = porthalfranges(T)

    # tmp = [-I T12; 0 T22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    @views tmp[range1,range2] .= T[range1,range2]
    @views tmp[range2,range2] .= T[range2,range2]

    # S = [0 -T11; I -T21]
    fill!(S,zero(eltype(S)))
    for d in range1
        S[d+length(range1),d] = 1
    end
    @views S[range1,range2] .= .-T[range1,range1]
    @views S[range2,range2] .= .-T[range2,range1]

    # perform the left division
    # S = inv(tmp)*S
    ldiv!(lu!(tmp), S)

    return nothing
end


@doc """
    ZtoS(Z;portimpedances=50.0)

Convert the impedance parameter matrix `Z` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances. Assumes a port impedance of 50 Ohms unless specified with
the `portimpedances` keyword argument. A complex port impedance `Z_k`
defines the pseudo-waves `a_k = (V_k + Z_k I_k)/(2 sqrt(Z_k))` and
`b_k = (V_k - Z_k I_k)/(2 sqrt(Z_k))`, which a load of impedance `Z_k`
matches, rather than power waves.

# Examples
```jldoctest
julia> Z = Complex{Float64}[0.0 0.0;0.0 0.0];JosephsonCircuits.ZtoS(Z)
2×2 Matrix{ComplexF64}:
 -1.0+0.0im   0.0+0.0im
  0.0+0.0im  -1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ZtoS

"""
    ZtoS!(S::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)

In place version of [`ZtoS`](@ref), writing into `S` with `tmp` as scratch
and `oneoversqrtportimpedances` the reciprocal square roots of the port
impedances.

"""
function ZtoS!(S::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix,oneoversqrtportimpedances)
    
    # compute \tilde{Z} = oneoversqrtportimpedances*Z*oneoversqrtportimpedances in S and tmp
    copy!(S,Z)
    rmul!(S,oneoversqrtportimpedances)
    lmul!(oneoversqrtportimpedances,S)
    copy!(tmp,S)

    # tmp = (\tilde{Z} + I)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # S = (\tilde{Z} - I)
    for d in 1:size(S,1)
        S[d,d] -= 1
    end

    # perform the left division
    # S = inv(tmp)*S = (\tilde{Z} + I) \ (\tilde{Z} - I)
    ldiv!(lu!(tmp), S)

    return nothing
end

function ZtoY!(Y, Z, tmp)
    # Y = inv(Z), as Z \ I through the factorization of a copy of Z
    copy!(tmp, Z)
    fill!(Y, zero(eltype(Y)))
    for d in 1:size(Y, 1)
        Y[d, d] = one(eltype(Y))
    end
    ldiv!(lu!(tmp), Y)
    return nothing
end

function YtoZ!(Z, Y, tmp)
    # the same inversion as ZtoY!
    return ZtoY!(Z, Y, tmp)
end

@doc """
    ZtoA(Z)

Convert the impedance matrix `Z` to the ABCD matrix `A` and return the result.

# Examples
```jldoctest
julia> Z = Complex{Float64}[50.0 50.0;50.0 50.0];JosephsonCircuits.ZtoA(Z)
2×2 Matrix{ComplexF64}:
  1.0+0.0im  0.0-0.0im
 0.02+0.0im  1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ZtoA

"""
    ZtoA!(A::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)

See [`ZtoA`](@ref) for description.

"""
function ZtoA!(A::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)
    # the conversion between Z and A is its own inverse
    return AtoZ!(A,Z,tmp)
end

@doc """
    ZtoB(Z)

Convert the impedance matrix `Z` to the inverse chain matrix `B` and return
the result.

# Examples
```jldoctest
julia> Z = Complex{Float64}[50.0 50;50 50];JosephsonCircuits.ZtoB(Z)
2×2 Matrix{ComplexF64}:
  1.0+0.0im  0.0-0.0im
 0.02+0.0im  1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ZtoB

"""
    ZtoB!(B::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)

See [`ZtoB`](@ref) for description.

"""
function ZtoB!(B::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(Z)

    # tmp = [0 Z12; -I Z22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d+length(range1),d] = -1
    end
    @views tmp[range1,range2] .= Z[range1, range2]
    @views tmp[range2,range2] .= Z[range2, range2]

    # B = [I Z11; 0 Z21]
    fill!(B,zero(eltype(B)))
    for d in range1
        B[d,d] = 1
    end

    @views B[range1,range2] .= Z[range1, range1]
    @views B[range2,range2] .= Z[range2, range1]

    # perform the left division
    # B = inv(tmp)*B = [0 Z12; -I Z22] \ [I Z11; 0 Z21]
    ldiv!(lu!(tmp), B)

    return nothing
end

@doc """
    YtoS(Y;portimpedances=50.0)

Convert the admittance parameter matrix `Y` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances. Assumes a port impedance of 50 Ohms unless specified with
the `portimpedances` keyword argument.

# Examples
```jldoctest
julia> Y = Complex{Float64}[1/50.0 0.0;0.0 1/50.0];JosephsonCircuits.YtoS(Y)
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0-0.0im
 0.0-0.0im  0.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" YtoS

"""
    YtoS!(S::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

See [`YtoS`](@ref) for description.

"""
function YtoS!(S::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)
    
    # compute \tilde{Y} = sqrtportimpedances*Y*sqrtportimpedances in S and tmp
    copy!(S,Y)
    rmul!(S,sqrtportimpedances)
    lmul!(sqrtportimpedances,S)
    copy!(tmp,S)
    rmul!(S,-1)

    # tmp = (\tilde{Y} + I)
    for d in 1:size(tmp,1)
        tmp[d,d] += 1
    end

    # S = (-\tilde{Y} + I)
    for d in 1:size(S,1)
        S[d,d] += 1
    end

    # perform the left division
    # S = inv(tmp)*S = (I + \tilde{Y}) \ (I - \tilde{Y})
    ldiv!(lu!(tmp), S)

    return nothing
end

@doc """
    YtoA(Y)

Convert the admittance matrix `Y` to the chain (ABCD) matrix `A` and return
the result.

# Examples
```jldoctest
julia> Y = Complex{Float64}[1/50 1/50;1/50 1/50];JosephsonCircuits.YtoA(Y)
2×2 Matrix{ComplexF64}:
 -1.0+0.0im  -50.0+0.0im
  0.0+0.0im   -1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of overall sign on (suspected typo).
""" YtoA

"""
    YtoA!(A::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)

See [`YtoA`](@ref) for description.

"""
function YtoA!(A::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(Y)

    # tmp = [-Y11 I; Y21 0]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d+length(range1)] = 1
    end
    @views tmp[range1,range1] .= .-Y[range1, range1]
    @views tmp[range2,range1] .= Y[range2, range1]

    # A = [Y12 0; -Y22 -I]
    fill!(A,zero(eltype(A)))
    for d in range2
        A[d,d] = -1
    end

    @views A[range1,range1] .= Y[range1, range2]
    @views A[range2,range1] .= .-Y[range2, range2]

    # perform the left division
    # A = inv(tmp)*A = [-Y11 I; Y21 0] \ [Y12 0; -Y22 -I]
    ldiv!(lu!(tmp), A)

    return nothing
end

@doc """
    YtoB(Y)

Convert the admittance matrix `Y` to the inverse chain matrix `B` and return
the result.

# Examples
```jldoctest
julia> Y = Complex{Float64}[1/50 1/50;1/50 1/50];JosephsonCircuits.YtoB(Y)
2×2 Matrix{ComplexF64}:
 -1.0+0.0im  -50.0+0.0im
  0.0+0.0im   -1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" YtoB

"""
    YtoB!(B::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)

See [`YtoB`](@ref) for description.

"""
function YtoB!(B::AbstractMatrix,Y::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(Y)

    # tmp = [-Y12 0; -Y22 I]
    fill!(tmp,zero(eltype(tmp)))
    for d in range2
        tmp[d,d] = 1
    end
    @views tmp[range1,range1] .= .-Y[range1, range2]
    @views tmp[range2,range1] .= .-Y[range2, range2]

    # B = [Y11 I; Y21 0]
    fill!(B,zero(eltype(B)))
    for d in range1
        B[d,d+length(range1)] = 1
    end

    @views B[range1,range1] .= Y[range1, range1]
    @views B[range2,range1] .= Y[range2, range1]

    # perform the left division
    # B = inv(tmp)*B = [-Y12 0; -Y22 I] \ [Y11 I; Y21 0]
    ldiv!(lu!(tmp), B)

    return nothing
end


@doc """
    AtoS(A; portimpedances = 50.0)

Convert the chain (ABCD) matrix `A` to the scattering parameter matrix `S` and
return the result. The first half of the ports are the inputs of the chain
matrix and the second half its outputs, each with its own port impedances:
`portimpedances` is a scalar, a vector with one value per port, or a matrix
with one row per port and one column per frequency, 50 Ohms unless specified.

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoS

"""
    AtoS!(S::AbstractMatrix, A::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`AtoS`](@ref) for description.

"""
function AtoS!(S::AbstractMatrix, A::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1, range2 = porthalfranges(A)
    h = length(range1)

    # tmp = [-I*g1 A11*g2+A12/g2; I/g1 A21*g2+A22/g2]
    # S = [I*g1 -A11*g2+A12/g2; I/g1 -A21*g2+A22/g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2, each
    # scaling the columns of its block
    for j in 1:h
        g1 = portfactor(sqrtportimpedances1, j)
        g2 = portfactor(sqrtportimpedances2, j)
        for i in 1:h
            δ = i == j
            A11 = A[i, j]
            A12 = A[i, h+j]
            A21 = A[h+i, j]
            A22 = A[h+i, h+j]
            tmp[i, j] = -δ*g1
            tmp[i, h+j] = A11*g2 + A12/g2
            tmp[h+i, j] = δ/g1
            tmp[h+i, h+j] = A21*g2 + A22/g2
            S[i, j] = δ*g1
            S[i, h+j] = -A11*g2 + A12/g2
            S[h+i, j] = δ/g1
            S[h+i, h+j] = -A21*g2 + A22/g2
        end
    end

    # perform the left division
    # S = inv(tmp)*S
    ldiv!(lu!(tmp), S)

    return nothing
end

@doc """
    AtoZ(A)

Convert the ABCD matrix `A` to the impedance matrix `Z` and return the result.

# Examples
```jldoctest
julia> A = Complex{Float64}[1.0 0.0;1/50 1.0];JosephsonCircuits.AtoZ(A)
2×2 Matrix{ComplexF64}:
 50.0+0.0im  50.0+0.0im
 50.0+0.0im  50.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoZ

"""
    AtoZ!(Z::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

See [`AtoZ`](@ref) for description.

"""
function AtoZ!(Z::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(A)

    # tmp = [-I A11; 0 A21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    @views tmp[range1,range2] .= A[range1, range1]
    @views tmp[range2,range2] .= A[range2, range1]

    # Z = [0 A12; I A22]
    fill!(Z,zero(eltype(Z)))
    for d in range1
        Z[d+length(range1),d] = 1
    end

    @views Z[range1,range2] .= A[range1, range2]
    @views Z[range2,range2] .= A[range2, range2]

    # perform the left division
    # Z = inv(tmp)*Z = [-I A11; 0 A21] \ [0 A12; I A22]
    ldiv!(lu!(tmp), Z)

    return nothing

end

@doc """
    AtoY(A)

Convert the chain (ABCD) matrix `A` to the admittance matrix `Y` and return
the result.

# Examples
```jldoctest
julia> A = Complex{Float64}[-1 -50.0;0 -1];JosephsonCircuits.AtoY(A)
2×2 Matrix{ComplexF64}:
 0.02+0.0im  0.02+0.0im
 0.02+0.0im  0.02+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoY

"""
    AtoY!(Y::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

See [`AtoY`](@ref) for description.

"""
function AtoY!(Y::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(A)

    # tmp = [0 A12; I A22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d+length(range1),d] = 1
    end
    @views tmp[range1,range2] .= A[range1, range2]
    @views tmp[range2,range2] .= A[range2, range2]

    # Y = [-I A11; 0 A21]
    fill!(Y,zero(eltype(Y)))
    for d in range1
        Y[d,d] = -1
    end

    @views Y[range1,range2] .= A[range1, range1]
    @views Y[range2,range2] .= A[range2, range1]

    # perform the left division
    # Y = inv(tmp)*Y = [0 A12; I A22] \ [-I A11; 0 A21]
    ldiv!(lu!(tmp), Y)

    return nothing

end

@doc """
    AtoB(A)

Convert the chain (ABCD) matrix `A` to the inverse chain matrix `B` and return
the result. Note that despite the name, the inverse of the chain
matrix is not equal to the inverse chain matrix, inv(A) ≠ B.

# Examples
```jldoctest
julia> A = Complex{Float64}[1.0 0.0;1/50 1.0];JosephsonCircuits.AtoB(A)
2×2 Matrix{ComplexF64}:
  1.0+0.0im  0.0+0.0im
 0.02+0.0im  1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" AtoB

"""
    AtoB!(B::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)

See [`AtoB`](@ref) for description.

"""
function AtoB!(B::AbstractMatrix,A::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(A)

    # tmp = [A11 -A12; A21 -A22]
    copy!(tmp,A)
    @views tmp[range1,range2] .*= -1
    @views tmp[range2,range2] .*= -1

    # B = [I 0; 0 -I]
    fill!(B,zero(eltype(B)))
    for d in range1
        B[d,d] = 1
    end
    for d in range2
        B[d,d] = -1
    end

    # perform the left division
    # B = inv(tmp)*B = [A11 -A12; A21 -A22] \ [I 0; 0 -I]
    ldiv!(lu!(tmp), B)

    return nothing

end

@doc """
    BtoS(B; portimpedances = 50.0)

Convert the inverse chain (ABCD) matrix `B` to the scattering parameter matrix
`S` and return the result. The first half of the ports are the inputs of the
chain matrix and the second half its outputs, each with its own port
impedances: `portimpedances` is a scalar, a vector with one value per port, or
a matrix with one row per port and one column per frequency, 50 Ohms unless
specified.

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of overall sign (suspected typo).
""" BtoS

"""
    BtoS!(S::AbstractMatrix, B::AbstractMatrix, tmp::AbstractMatrix,
        sqrtportimpedances1, sqrtportimpedances2)

See [`BtoS`](@ref) for description.

"""
function BtoS!(S::AbstractMatrix, B::AbstractMatrix, tmp::AbstractMatrix,
    sqrtportimpedances1, sqrtportimpedances2)

    range1, range2 = porthalfranges(B)
    h = length(range1)

    # tmp = [B11*g1+B12/g1 -I*g2; B21*g1+B22/g1 I/g2]
    # S = [-B11*g1+B12/g1 I*g2; -B21*g1+B22/g1 I/g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2, each
    # scaling the columns of its block
    for j in 1:h
        g1 = portfactor(sqrtportimpedances1, j)
        g2 = portfactor(sqrtportimpedances2, j)
        for i in 1:h
            δ = i == j
            B11 = B[i, j]
            B12 = B[i, h+j]
            B21 = B[h+i, j]
            B22 = B[h+i, h+j]
            tmp[i, j] = B11*g1 + B12/g1
            tmp[i, h+j] = -δ*g2
            tmp[h+i, j] = B21*g1 + B22/g1
            tmp[h+i, h+j] = δ/g2
            S[i, j] = -B11*g1 + B12/g1
            S[i, h+j] = δ*g2
            S[h+i, j] = -B21*g1 + B22/g1
            S[h+i, h+j] = δ/g2
        end
    end

    # perform the left division
    # S = inv(tmp)*S
    ldiv!(lu!(tmp), S)

    return nothing
end

@doc """
    BtoZ(B)

Convert the inverse chain matrix `B` to the impedance matrix `Z` and return
the result.

# Examples
```jldoctest
julia> B = Complex{Float64}[1.0 0.0;1/50 1];JosephsonCircuits.BtoZ(B)
2×2 Matrix{ComplexF64}:
 50.0+0.0im  50.0+0.0im
 50.0+0.0im  50.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006
with change of sign on B21 and B22 terms (suspected typo).
""" BtoZ

"""
    BtoZ!(Z::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)

See [`BtoZ`](@ref) for description.

"""
function BtoZ!(Z::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(B)

    # tmp = [B11 -I; -B21 0]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d+length(range1)] = -1
    end
    @views tmp[range1,range1] .= B[range1, range1]
    @views tmp[range2,range1] .= .-B[range2, range1]

    # Z = [B12 0; -B22 -I]
    fill!(Z,zero(eltype(Z)))
    for d in range2
        Z[d,d] = -1
    end

    @views Z[range1,range1] .= B[range1, range2]
    @views Z[range2,range1] .= .-B[range2, range2]

    # perform the left division
    # Z = inv(tmp)*Z = [B11 -I; -B21 0] \ [B12 0; -B22 -I]
    ldiv!(lu!(tmp), Z)

    return nothing
end

@doc """
    BtoY(B)

Convert the inverse chain matrix `B` to the admittance matrix `Y` and return
the result.

# Examples
```jldoctest
julia> B = Complex{Float64}[-1 -50.0;0 -1];JosephsonCircuits.BtoY(B)
2×2 Matrix{ComplexF64}:
 0.02+0.0im  0.02+0.0im
 0.02+0.0im  0.02+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" BtoY

"""
    BtoY!(Y::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)

See [`BtoY`](@ref) for description.

"""
function BtoY!(Y::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)
    
    range1, range2 = porthalfranges(B)

    # tmp = [B12 0; B22 I]
    fill!(tmp,zero(eltype(tmp)))
    for d in range2
        tmp[d,d] = 1
    end
    @views tmp[range1,range1] .= B[range1, range2]
    @views tmp[range2,range1] .= B[range2, range2]

    # Y = [B11 -I; B21 0]
    fill!(Y,zero(eltype(Y)))
    for d in range1
        Y[d,d+length(range1)] = -1
    end

    @views Y[range1,range1] .= B[range1, range1]
    @views Y[range2,range1] .= B[range2, range1]

    # perform the left division
    # Y = inv(tmp)*Y = [B12 0; B22 I] \ [B11 -I; B21 0]
    ldiv!(lu!(tmp), Y)

    return nothing
end

@doc """
    BtoA(B)

Convert the inverse chain matrix `B` to the chain (ABCD) matrix `A` and return
the result. Note that despite the name, the inverse of the chain
matrix is not equal to the inverse chain matrix, inv(A) ≠ B.

# Examples
```jldoctest
julia> B = Complex{Float64}[1.0 0.0;1/50 1.0];JosephsonCircuits.BtoA(B)
2×2 Matrix{ComplexF64}:
  1.0+0.0im  0.0+0.0im
 0.02+0.0im  1.0+0.0im
```

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" BtoA


"""
    BtoA!(A::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)

See [`BtoA`](@ref) for description.

"""
function BtoA!(A::AbstractMatrix,B::AbstractMatrix,tmp::AbstractMatrix)
    # the conversion between A and B is its own inverse
    return AtoB!(A,B,tmp)
end

@doc """
    ABCDtoS(ABCD;portimpedances=50.0)

Convert the 2 port chain (ABCD) matrix `ABCD` to the scattering parameter
matrix `S` and return the result. Assumes a port impedance of 50 Ohms unless
specified with the `portimpedances` keyword argument.

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" ABCDtoS

function ABCDtoS!(ABCD::AbstractMatrix,portimpedances)
    return ABCDtoS!(ABCD,first(portimpedances),last(portimpedances))
end

function ABCDtoS!(A::AbstractMatrix,RS,RL)
    checktwoport(A)
    A11 = A[1,1]
    A12 = A[1,2]
    A21 = A[2,1]
    A22 = A[2,2]

    A[1,1] = (A11*RL+A12-A21*RS*RL-A22*RS)/(A11*RL+A12+A21*RS*RL+A22*RS)
    A[1,2] = 2*sqrt(RS*RL)*(A11*A22-A12*A21)/(A11*RL+A12+A21*RS*RL+A22*RS)
    A[2,1] = 2*sqrt(RS*RL)/(A11*RL+A12+A21*RS*RL+A22*RS)
    A[2,2] = (-A11*RL+A12-A21*RS*RL+A22*RS)/(A11*RL+A12+A21*RS*RL+A22*RS)
    return A
end


@doc """
    StoABCD(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to the 2 port chain (ABCD) matrix and
return the result. Assumes a port impedance of 50 Ohms unless specified with the
`portimpedances` keyword argument.

# References
Russer, Peter. Electromagnetics, Microwave Circuit, And Antenna Design for
Communications Engineering, Second Edition. Artech House, 2006.
""" StoABCD

function StoABCD!(S::AbstractMatrix,portimpedances)
    return StoABCD!(S,first(portimpedances),last(portimpedances))
end

function StoABCD!(S::AbstractMatrix,RS,RL)
    checktwoport(S)
    S11 = S[1,1]
    S12 = S[1,2]
    S21 = S[2,1]
    S22 = S[2,2]

    S[1,1] = sqrt(RS/RL)*((1+S11)*(1-S22)+S21*S12)/(2*S21)
    S[1,2] = sqrt(RS*RL)*((1+S11)*(1+S22)-S21*S12)/(2*S21)
    S[2,1] = 1/sqrt(RS*RL)*((1-S11)*(1-S22)-S21*S12)/(2*S21)
    S[2,2] = sqrt(RL/RS)*((1-S11)*(1+S22)+S21*S12)/(2*S21)
    return S
end
