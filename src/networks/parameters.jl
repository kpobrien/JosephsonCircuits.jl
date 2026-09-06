# The network parameter conversions, S to Z, Z to S, S to T and the rest.
# Each is a kernel which converts one matrix into another,
# `f!(y, x, tmp, ...)`, defined further down beside its docstring, and the
# per frequency forms, on a matrix or on an array of matrices with one per
# frequency, allocating and in place, come from the driver here. A port
# impedance argument is a number, the same at every port and frequency, a
# vector with one value per port, or a matrix with one row per port and one
# column per frequency.

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
one frequency, calls the kernel directly. One method per number of port arguments, none, one or two, rather
than a variadic one: the argument map a variadic driver needs is compiled
once per kernel and argument type, and that outweighed what it saved. The
type parameter makes the method specialize on the kernel, which Julia does
not do on its own for a function argument.
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
    convertperfrequency!(kernel!, similar(x), x, args...)

# a matrix is one frequency and is converted by the kernel directly: no
# loop, no views and no scratch of its own to compile per kernel and
# argument type, which for the single matrix calls of user code is nearly
# all of what the loop form would compile
function convertperfrequency(kernel!::F, x::AbstractMatrix) where {F}
    y = similar(x)
    kernel!(y, x, similar(x))
    return y
end
function convertperfrequency(kernel!::F, x::AbstractMatrix, a) where {F}
    checkconversion(x, x, a)
    y = similar(x)
    kernel!(y, x, similar(x), atmatrix(a))
    return y
end
function convertperfrequency(kernel!::F, x::AbstractMatrix, a, b) where {F}
    checkconversion(x, x, a, b)
    y = similar(x)
    kernel!(y, x, similar(x), atmatrix(a), atmatrix(b))
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
    y = copy(x)
    for i in CartesianIndices(axes(x)[3:end])
        kernel!(view(y, :, :, i), atfrequency(a, i))
    end
    return y
end
function convertcopy(kernel!::F, x::AbstractMatrix, a) where {F}
    checkconversion(x, x, a)
    y = copy(x)
    kernel!(y, atmatrix(a))
    return y
end

# the checks of the driver, and its scratch matrix: the output like the
# input, and every port argument with one column per frequency of the
# input or none
function checkconversion(y, x, args...)
    axes(y) == axes(x) || throw(DimensionMismatch(
        lazy"Sizes of output $(size(y)) and input $(size(x)) must be equal."))
    trailing = axes(x)[3:end]
    for a in args
        fa = frequencyaxes(a)
        (isnothing(fa) || fa == trailing) || throw(ArgumentError(
            lazy"A port argument of size $(size(a isa PortDiagonal ? a.values : a)) does not have one column per frequency of the input of size $(size(x))."))
    end
    return similar(x, axes(x)[1:2])
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
 -0.0-0.0im   1.0-0.0im
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
    
    range1 = 1:size(T,1)÷2
    range2 = size(T,1)÷2+1:size(T,1)

    # tmp = [-I S11; 0 S21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -one(eltype(tmp))
    end
    tmp[range1,range2] .= S[range1, range1]
    tmp[range2,range2] .= S[range2, range1]

    # T = [-S12 0; -S22 I]
    fill!(T,zero(eltype(T)))
    for d in range2
        T[d,d] = one(eltype(T))
    end

    T[range1,range1] .= -S[range1, range2]
    T[range2,range1] .= -S[range2, range2]

    # perform the left division
    # T = inv(tmp)*T = [-I S11; 0 S21] \ [-S12 0; -S22 I]
    T .= tmp \ T

    return nothing
end


@doc """
    StoZ(S;portimpedances=50.0)

Convert the scattering parameter matrix `S` to an impedance parameter matrix
`Z` and return the result. Assumes a port impedance of 50 Ohms unless
specified with the `portimpedances` keyword argument.

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
    Z .= tmp \ Z

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
    Y .= tmp \ Y

    # left multiply by sqrt(z)
    # compute oneoversqrtportimpedances*((I - S) \ ((I + S)*oneoversqrtportimpedances))
    lmul!(oneoversqrtportimpedances,Y)
    return nothing
end


@doc """
    StoA(S)

Convert the scattering parameter matrix `S` to the chain (ABCD) matrix `A` and
return the result.

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

    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    A11 = view(A,range1,range1)
    A12 = view(A,range1,range2)
    A21 = view(A,range2,range1)
    A22 = view(A,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # define the matrices

    # tmp = -[(I-S11)/g1 -(I+S11)*g1; -S21*g1 -S21*g1]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    tmp11 .= -(I - S11)/sqrtportimpedances1
    tmp12 .= (I + S11)*sqrtportimpedances1
    tmp21 .= S21/sqrtportimpedances1
    tmp22 .= S21*sqrtportimpedances1

    # A = [-S12/g2 S12*g2; (I-S22)/g2 (I+S22)*g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    A11 .= -S12/sqrtportimpedances2
    A12 .= S12*sqrtportimpedances2
    A21 .= (I-S22)/sqrtportimpedances2
    A22 .= (I+S22)*sqrtportimpedances2

    # perform the left division

    # A = inv(tmp)*A
    A .= tmp \ A

    return nothing
end


@doc """
    StoB(S)

Convert the scattering parameter matrix `S` to the inverse chain (ABCD) matrix
`B` and return the result. Note that despite the name, the inverse of the chain
matrix is not equal to the inverse chain matrix, inv(A) ≠ B.

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

    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    B11 = view(B,range1,range1)
    B12 = view(B,range1,range2)
    B21 = view(B,range2,range1)
    B22 = view(B,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # define the matrices

    # tmp = [S12/g2 S12*g2; -(I-S22)/g2 (I+S22)*g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    tmp11 .= S12/sqrtportimpedances2
    tmp12 .= S12*sqrtportimpedances2
    tmp21 .= -(I-S22)/sqrtportimpedances2
    tmp22 .= (I+S22)*sqrtportimpedances2

    # B = [(I-S11)/g1 (I+S11)*g1; -S21/g1 S21*g1]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    B11 .= (I-S11)/sqrtportimpedances1
    B12 .= (I+S11)*sqrtportimpedances1
    B21 .= -S21/sqrtportimpedances1
    B22 .= S21*sqrtportimpedances1

    # perform the left division

    # B = inv(tmp)*B
    B .= tmp \ B

    return nothing
end


@doc """
    TtoS(T)

Convert the transmission matrix `T` to a scattering parameter matrix `S` and return the result.

# Examples
```jldoctest
julia> T = Complex{Float64}[1.0 0.0;0.0 1.0];JosephsonCircuits.TtoS(T)
2×2 Matrix{ComplexF64}:
 -0.0-0.0im   1.0+0.0im
  1.0-0.0im  -0.0-0.0im
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
    
    range1 = 1:size(T,1)÷2
    range2 = size(T,1)÷2+1:size(T,1)

    # make views of the block matrices
    # S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    # S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    T11 = view(T,range1,range1)
    T12 = view(T,range1,range2)
    T21 = view(T,range2,range1)
    T22 = view(T,range2,range2)

    # tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    # tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # tmp = [-I T12; 0 T22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp12 .= T12
    tmp22 .= T22

    # S = [0 T11; I T21]
    fill!(S,zero(eltype(S)))
    for d in range1
        S[d+size(T,1)÷2,d] = 1
    end

    S12 .= -T11
    S22 .= -T21

    # perform the left division
    # S = inv(tmp)*S
    S .= tmp \ S

    return nothing
end


@doc """
    ZtoS(Z;portimpedances=50.0)

Convert the impedance parameter matrix `Z` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances. Assumes a port impedance of 50 Ohms unless specified with
the `portimpedances` keyword argument.

# Examples
```jldoctest
julia> Z = Complex{Float64}[0.0 0.0;0.0 0.0];JosephsonCircuits.ZtoS(Z)
2×2 Matrix{ComplexF64}:
 -1.0+0.0im   0.0-0.0im
  0.0-0.0im  -1.0+0.0im
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
    S .= tmp \ S

    return nothing
end

function ZtoY!(Y, Z, tmp)
    Y .= inv(Z)
    return nothing
end

function YtoZ!(Z, Y, tmp)
    Z .= inv(Y)
    return nothing
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
    
    range1 = 1:size(Z,1)÷2
    range2 = size(Z,1)÷2+1:size(Z,1)

    # tmp = [0 Z12; -I Z22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d+size(Z,1)÷2,d] = -1
    end
    tmp[range1,range2] .= Z[range1, range2]
    tmp[range2,range2] .= Z[range2, range2]

    # B = [I Z11; 0 Z21]
    fill!(B,zero(eltype(B)))
    for d in range1
        B[d,d] = 1
    end

    B[range1,range2] .= Z[range1, range1]
    B[range2,range2] .= Z[range2, range1]

    # println(tmp)
    # println(B)
    # perform the left division
    # B = inv(tmp)*B = [0 Z12; -I Z22] \ [I Z11; 0 Z21]
    B .= tmp \ B

    return nothing
end

@doc """
    YtoS(Y;portimpedances=50.0)

Convert the admittance parameter matrix `Y` to a scattering parameter matrix
`S` and return the result. `portimpedances` is a scalar, vector, or matrix of
port impedances.

# Examples
```jldoctest
julia> Y = Complex{Float64}[1/50.0 0.0;0.0 1/50.0];JosephsonCircuits.YtoS(Y)
2×2 Matrix{ComplexF64}:
  0.0-0.0im  -0.0-0.0im
 -0.0-0.0im   0.0-0.0im
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
    S .= tmp \ S

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
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [-Y11 I; Y21 0]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d+size(A,1)÷2] = 1
    end
    tmp[range1,range1] .= -Y[range1, range1]
    tmp[range2,range1] .= Y[range2, range1]

    # A = [Y12 0; -Y22 -I]
    fill!(A,zero(eltype(A)))
    for d in range2
        A[d,d] = -1
    end

    A[range1,range1] .= Y[range1, range2]
    A[range2,range1] .= -Y[range2, range2]

    # perform the left division
    # A = inv(tmp)*A = [-Y11 I; Y21 0] \ [Y12 0; -Y22 -I]
    A .= tmp \ A

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
    
    range1 = 1:size(Y,1)÷2
    range2 = size(Y,1)÷2+1:size(Y,1)

    # tmp = [-Y12 0; -Y22 I]
    fill!(tmp,zero(eltype(tmp)))
    for d in range2
        tmp[d,d] = 1
    end
    tmp[range1,range1] .= -Y[range1, range2]
    tmp[range2,range1] .= -Y[range2, range2]

    # B = [Y11 I; Y21 0]
    fill!(B,zero(eltype(B)))
    for d in range1
        B[d,d+size(B,1)÷2] = 1
    end

    B[range1,range1] .= Y[range1, range1]
    B[range2,range1] .= Y[range2, range1]

    # perform the left division
    # B = inv(tmp)*B = [-Y12 0; -Y22 I] \ [Y11 I; Y21 0]
    B .= tmp \ B

    return nothing
end


@doc """
    AtoS(A)

Convert the chain (ABCD) matrix `A` to the scattering parameter matrix `S` and
return the result.

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

    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    A11 = view(A,range1,range1)
    A12 = view(A,range1,range2)
    A21 = view(A,range2,range1)
    A22 = view(A,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # tmp = [-g1 A11*g2+A12/g2; 1/g1 A21*g2+A22/g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    fill!(tmp,zero(eltype(tmp)))

    for d in range1
        tmp[d,d] = 1
    end

    for d in range1
        tmp[d+size(A,1)÷2,d] = 1
    end

    tmp11 .= -tmp11*sqrtportimpedances1
    tmp12 .= A11*sqrtportimpedances2+A12/sqrtportimpedances2
    tmp21 .= tmp21/sqrtportimpedances1
    tmp22 .= A21*sqrtportimpedances2+A22/sqrtportimpedances2

    # S = [g1 -A11*g2+A12/g2; 1/g1 -A21*g2+A22/g2]
    # where g1 = sqrtportimpedances1 and g2 = sqrtportimpedances2
    fill!(S,zero(eltype(S)))

    for d in range1
        S[d,d] = 1
    end

    for d in range1
        S[d+size(A,1)÷2,d] = 1
    end

    S11 .= S11*sqrtportimpedances1
    S12 .= -A11*sqrtportimpedances2+A12/sqrtportimpedances2
    S21 .=  S21/sqrtportimpedances1
    S22 .= -A21*sqrtportimpedances2+A22/sqrtportimpedances2

    # perform the left division
    # S = inv(tmp)*S
    S .= tmp \ S

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
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [-I A11; 0 A21]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d] = -1
    end
    tmp[range1,range2] .= A[range1, range1]
    tmp[range2,range2] .= A[range2, range1]

    # Z = [0 A12; I A22]
    fill!(Z,zero(eltype(A)))
    for d in range1
        Z[d+size(A,1)÷2,d] = 1
    end

    Z[range1,range2] .= A[range1, range2]
    Z[range2,range2] .= A[range2, range2]

    # perform the left division
    # Z = inv(tmp)*Z = [-I A11; 0 A21] \ [0 A12; I A22]
    Z .= tmp \ Z

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
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [0 A12; I A22]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d+size(A,1)÷2,d] = 1
    end
    tmp[range1,range2] .= A[range1, range2]
    tmp[range2,range2] .= A[range2, range2]

    # Y = [-I A11; 0 A21]
    fill!(Y,zero(eltype(A)))
    for d in range1
        Y[d,d] = -1
    end

    Y[range1,range2] .= A[range1, range1]
    Y[range2,range2] .= A[range2, range1]

    # perform the left division
    # Y = inv(tmp)*Z = [0 A12; I A22] \ [-I A11; 0 A21]
    Y .= tmp \ Y

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
 0.02-0.0im  1.0-0.0im
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
    
    range1 = 1:size(A,1)÷2
    range2 = size(A,1)÷2+1:size(A,1)

    # tmp = [A11 -A12; A21 -A22]
    copy!(tmp,A)
    tmp[range1,range2] .*= -1
    tmp[range2,range2] .*= -1

    # B = [I 0; 0 I]
    fill!(B,zero(eltype(A)))
    for d in range1
        B[d,d] = 1
    end
    for d in range2
        B[d,d] = -1
    end

    # perform the left division
    # B = inv(tmp)*B = [A11 -A12; A21 -A22] \ [I 0; 0 I]
    B .= tmp \ B

    return nothing

end

@doc """
    BtoS(B)

Convert the inverse chain (ABCD) matrix `B` to the scattering parameter matrix
`S` and return the result.

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

    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # make views of the block matrices
    S11 = view(S,range1,range1)
    S12 = view(S,range1,range2)
    S21 = view(S,range2,range1)
    S22 = view(S,range2,range2)

    B11 = view(B,range1,range1)
    B12 = view(B,range1,range2)
    B21 = view(B,range2,range1)
    B22 = view(B,range2,range2)

    tmp11 = view(tmp,range1,range1)
    tmp12 = view(tmp,range1,range2)
    tmp21 = view(tmp,range2,range1)
    tmp22 = view(tmp,range2,range2)

    # tmp = [B11*g1+B12/g1 -g2; B21*g1+B22/g1 1/g2]
    fill!(tmp,zero(eltype(tmp)))

    for d in range1
        tmp[d,d+size(B,1)÷2] = 1
    end

    for d in range2
        tmp[d,d] = 1
    end

    tmp11 .= B11*sqrtportimpedances1+B12/sqrtportimpedances1
    tmp12 .= -tmp12*sqrtportimpedances2
    tmp21 .= B21*sqrtportimpedances1+B22/sqrtportimpedances1
    tmp22 .= tmp22/sqrtportimpedances2

    # S = [-B11*g1+B12/g1 g2; -B21*g1+B22/g1 1/g2]
    fill!(S,zero(eltype(S)))

    for d in range1
        S[d,d+size(B,1)÷2] = 1
    end

    for d in range2
        S[d,d] = 1
    end

    S11 .= -B11*sqrtportimpedances1+B12/sqrtportimpedances1
    S12 .= S12*sqrtportimpedances2
    S21 .= -B21*sqrtportimpedances1+B22/sqrtportimpedances1
    S22 .= S22/sqrtportimpedances2

    # perform the left division
    # S = inv(tmp)*S
    S .= tmp \ S

    return nothing
end

@doc """
    BtoZ(A)

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
    
    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # tmp = [B11 -I; B21 0]
    fill!(tmp,zero(eltype(tmp)))
    for d in range1
        tmp[d,d+size(B,1)÷2] = -1
    end
    tmp[range1,range1] .= B[range1, range1]
    tmp[range2,range1] .= -B[range2, range1]

    # Z = [B12 0; B22 -I]
    fill!(Z,zero(eltype(B)))
    for d in range2
        Z[d,d] = -1
    end

    Z[range1,range1] .= B[range1, range2]
    Z[range2,range1] .= -B[range2, range2]

    # perform the left division
    # Y = inv(tmp)*Y = [B12 0; B22 I] \ [B11 -I; B21 0]
    # Z = inv(tmp)*Z = [B11 -I; B21 0] \ [B12 0; B22 -I]
    Z .= tmp \ Z

    return nothing
end

@doc """
    BtoY(A)

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
    
    range1 = 1:size(B,1)÷2
    range2 = size(B,1)÷2+1:size(B,1)

    # tmp = [B12 0; B22 I]
    fill!(tmp,zero(eltype(tmp)))
    for d in range2
        tmp[d,d] = 1
    end
    tmp[range1,range1] .= B[range1, range2]
    tmp[range2,range1] .= B[range2, range2]

    # Y = [B11 -I; B21 0]
    fill!(Y,zero(eltype(B)))
    for d in range1
        Y[d,d+size(B,1)÷2] = -1
    end

    Y[range1,range1] .= B[range1, range1]
    Y[range2,range1] .= B[range2, range1]

    # perform the left division
    # Y = inv(tmp)*Y = [B12 0; B22 I] \ [B11 -I; B21 0]
    Y .= tmp \ Y

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
 0.02-0.0im  1.0-0.0im
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
    StoABCD(S;portimpedances=50.0))

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
