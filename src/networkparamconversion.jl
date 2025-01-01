# generate the non in-place versions of the network parameter conversion
# functions.

# functions without an impedance argument
for f in [:StoT, :TtoS, :AtoB, :BtoA, :AtoZ, :ZtoA, :AtoY, :YtoA, :BtoY, :YtoB, :BtoZ, :ZtoB, :ZtoY, :YtoZ]
    # non-in-place version for matrix input
    @eval function ($f)(x::AbstractMatrix)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # evaluate the in place version of the function
        ($(Symbol(String(f)*"!")))(y,x,tmp)
        return y
    end

    # non-in-place version for array input
    @eval function ($f)(x::AbstractArray)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # evaluate the in place version of the function
        for i in CartesianIndices(axes(x)[3:end])
            ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp)
        end
        return y
    end

    # in place version for matrix or array input, copy results to input
    @eval function ($(Symbol(String(f)*"!")))(x)
        return copy!(x,($f)(x))
    end

end

# functions using portimpedances
for f in [:ABCDtoS, :StoABCD]
    # non-in-place version for matrix input
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = copy(x)
        # evaluate the in place version of the function
        ($(Symbol(String(f)*"!")))(y,portimpedances)
        return y
    end

    # non-in-place version for array input
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = copy(x)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(y)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),portimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(y)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(portimpedances,:,i))
            end
        end
        return y
    end

    # in place version for matrix or array input, copy results to input
    @eval function ($(Symbol(String(f)*"!")))(x;portimpedances=50.0)
        return copy!(x,($f)(x;portimpedances=portimpedances))
    end
end

# functions using sqrtportimpedances
for f in [:StoZ, :YtoS]
    # non-in-place version for matrix input
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,sqrtportimpedances)
        else
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,Diagonal(sqrtportimpedances))
        end
        return y
    end
    
    # non-in-place version for array input
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,sqrtportimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,Diagonal(view(sqrtportimpedances,:,i)))
            end
        end
        return y
    end

    # in place version for matrix or array input, copy results to input
    @eval function ($(Symbol(String(f)*"!")))(x;portimpedances=50.0)
        return copy!(x,($f)(x;portimpedances=portimpedances))
    end
end


# functions using oneoversqrtportimpedances
for f in [:StoY, :ZtoS]
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # calculate the inverse of the square root of the port impedances
        oneoversqrtportimpedances = 1 ./sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,oneoversqrtportimpedances)
        else
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,Diagonal(oneoversqrtportimpedances))
        end
        return y
    end
    
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # calculate the inverse of the square root of the port impedances
        oneoversqrtportimpedances = 1 ./sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,oneoversqrtportimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,Diagonal(view(oneoversqrtportimpedances,:,i)))
            end
        end
        return y
    end

    # in place version for matrix or array input, copy results to input
    @eval function ($(Symbol(String(f)*"!")))(x;portimpedances=50.0)
        return copy!(x,($f)(x;portimpedances=portimpedances))
    end
end

# functions using sqrtportimpedances split into two
for f in [:StoA, :AtoS, :StoB, :BtoS]
    @eval function ($f)(x::AbstractMatrix;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x)
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,sqrtportimpedances,
                sqrtportimpedances)
        else
            # evaluate the in place version of the function
            ($(Symbol(String(f)*"!")))(y,x,tmp,
                Diagonal(sqrtportimpedances[1:size(sqrtportimpedances,1)÷2]),
                Diagonal(sqrtportimpedances[size(sqrtportimpedances,1)÷2+1:size(sqrtportimpedances,1)]))
        end
        return y
    end
    
    @eval function ($f)(x::AbstractArray;portimpedances=50.0)
        # define the output matrix
        y = similar(x)
        # define a temporary storage array
        tmp = similar(x,axes(x)[1:2])
        # calculate the square root of the port impedances
        sqrtportimpedances = sqrt.(portimpedances)
        # if portimpedances is a scalar, pass that, otherwise pass as a
        # diagonal matrix
        if iszero(ndims(portimpedances))
            # assume the port impedances are all the same for all ports and
            # frequencies. loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,
                    sqrtportimpedances, sqrtportimpedances)
            end
        else
            # assume the port impedances are given for each port and frequency
            # loop over the dimensions of the array greater than 2
            for i in CartesianIndices(axes(x)[3:end])
                ($(Symbol(String(f)*"!")))(view(y,:,:,i),view(x,:,:,i),tmp,
                Diagonal(view(sqrtportimpedances,1:size(sqrtportimpedances,1)÷2,i)),
                Diagonal(view(sqrtportimpedances,size(sqrtportimpedances,1)÷2+1:size(sqrtportimpedances,1),i)))
            end
        end
        return y
    end

    # in place version for matrix or array input, copy results to input
    @eval function ($(Symbol(String(f)*"!")))(x;portimpedances=50.0)
        return copy!(x,($f)(x;portimpedances=portimpedances))
    end
end


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
    StoY!(Y::AbstractMatrix,S::AbstractMatrix,tmp::AbstractMatrix,sqrtportadmittances)

See [`StoY`](@ref) for description.

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
    ZtoS!(S::AbstractMatrix,Z::AbstractMatrix,tmp::AbstractMatrix,sqrtportimpedances)

See [`ZtoS`](@ref) for description.

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
