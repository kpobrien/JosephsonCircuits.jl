# generate the network parameter convenience functions based on the functions below

# one argument functions
for f in [:ABCD_seriesZ, :ABCD_shuntY, :Y_seriesY, :Z_shuntZ]

    # non-in-place version for array x1
    @eval function $f(x1::AbstractArray)
        # define the output matrix
        y = zeros(eltype(x1),(2,2,size(x1)...))
        # evaluate the in place version of the function
        $(Symbol(String(f)*"!"))(y,x1)
        return y
    end

    # in-place version for scalar x1
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray,x1::Number)
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1)
        end
        return y
    end

    # in-place version for abstract array x1
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray,x1::AbstractArray)
        if size(y)[3:end] != size(x1)
            throw(ArgumentError("Sizes of output $(size(y)) and input $(size(x1)) not compatible."))
        end
        # assume a value of x1 is given for each frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1[i])
        end
        return y
    end
end

# two argument functions, complex
for f in [:ABCD_tline, :Z_tline]

    # non-in-place version for array x1, x2
    @eval function $f(x1::AbstractArray, x2::AbstractArray)
        # check the sizes of the inputs
        if size(x1) != size(x2)
            throw(ArgumentError("Sizes of inputs $(size(x1)) and $(size(x2)) must be equal."))
        end
        # define the output matrix
        y = zeros(promote_type(typeof(im*x1[1]),eltype(x2)),(2,2,size(x1)...))
        # evaluate the in place version of the function
        $(Symbol(String(f)*"!"))(y,x1,x2)
        return y
    end

    # non-in-place version for scalar x1, array x2
    @eval function $f(x1::Number, x2::AbstractArray)
        # define the output matrix
        y = zeros(promote_type(typeof(im*x1),eltype(x2)),(2,2,size(x2)...))
        # evaluate the in place version of the function
        $(Symbol(String(f)*"!"))(y,x1,x2)
        return y
    end

    # in-place version for scalar x1, x2
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray, x1::Number,
        x2::Number)
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1,x2)
        end
        return y
    end

    # in-place version for scalar x1, abstract array x2
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray,
        x1::Number, x2::AbstractArray)
        # check the sizes of the inputs
        if size(y)[3:end] != size(x2)
            throw(ArgumentError("Sizes of output $(size(y)) and inputs $(size(x2)) not compatible."))
        end
        # assume a value of x1 is given for each frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1,x2[i])
        end
        return y
    end

    # in-place version for abstract array x1, x2
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray,
        x1::AbstractArray, x2::AbstractArray)
        # check the sizes of the inputs
        if size(x1) != size(x2)
            throw(ArgumentError("Sizes of inputs $(size(x1)) and $(size(x2)) must be equal."))
        end
        if size(y)[3:end] != size(x1)
            throw(ArgumentError("Sizes of output $(size(y)) and inputs $(size(x1)) not compatible."))
        end
        # assume a value of x1 is given for each frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1[i],x2[i])
        end
        return y
    end
end

# three argument functions
for f in [:ABCD_PiY, :Y_PiY, :ABCD_TZ, :Z_TZ]

    # non-in-place version for array x1, x2, x3
    @eval function $f(x1::AbstractArray, x2::AbstractArray, x3::AbstractArray)
        # check the sizes of the inputs
        if size(x1) != size(x2) || size(x1) != size(x3)
            throw(ArgumentError("Sizes of inputs $(size(x1)), $(size(x2)), and $(size(x3)) must be equal."))
        end
        # define the output matrix
        y = zeros(promote_type(eltype(x1),eltype(x2),eltype(x3)),(2,2,size(x1)...))
        # evaluate the in place version of the function
        $(Symbol(String(f)*"!"))(y,x1,x2,x3)
        return y
    end

    # in-place version for scalar x1, x2, x3
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray, x1::Number,
        x2::Number, x3::Number)
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1,x2,x3)
        end
        return y
    end

    # in-place version for abstract array x1, x2, x3
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray,
        x1::AbstractArray, x2::AbstractArray, x3::AbstractArray)
        # check the sizes of the inputs
        if size(x1) != size(x2) || size(x1) != size(x3)
            throw(ArgumentError("Sizes of inputs $(size(x1)), $(size(x2)), and $(size(x3)) must be equal."))
        end
        if size(y)[3:end] != size(x1)
            throw(ArgumentError("Sizes of output $(size(y)) and inputs $(size(x1)) not compatible."))
        end
        # assume a value of x1 is given for each frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1[i],x2[i],x3[i])
        end
        return y
    end
end

# four argument functions, complex
for f in [:ABCD_coupled_tline, :Z_coupled_tline]

    # non-in-place version for array x1, x2, x3, x4
    @eval function $f(x1::AbstractArray, x2::AbstractArray,
        x3::AbstractArray, x4::AbstractArray)
        # check the sizes of the inputs
        if size(x1) != size(x2) || size(x1) != size(x3) || size(x1) != size(x4)
            throw(ArgumentError("Sizes of inputs $(size(x1)), $(size(x2)), $(size(x3)), and $(size(x4)) must be equal."))
        end
        # define the output matrix
        y = zeros(promote_type(typeof(im*x1[1]),typeof(im*x2[1]),eltype(x3),eltype(x4)),(4,4,size(x1)...))
        # evaluate the in place version of the function
        $(Symbol(String(f)*"!"))(y,x1,x2,x3,x4)
        return y
    end

    # non-in-place version for scalar x1, x2, array x3, x4
    @eval function $f(x1::Number, x2::Number, x3::AbstractArray,
        x4::AbstractArray)
        # check the sizes of the inputs
        if size(x3) != size(x4)
            throw(ArgumentError("Sizes of inputs $(size(x3)) and $(size(x4)) must be equal."))
        end
        # define the output matrix
        y = zeros(promote_type(typeof(im*x1),typeof(im*x2),eltype(x3),eltype(x4)),(4,4,size(x3)...))
        # evaluate the in place version of the function
        $(Symbol(String(f)*"!"))(y,x1,x2,x3,x4)
        return y
    end

    # in-place version for scalar x1, x2, x3, x4
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray, x1::Number,
        x2::Number, x3::Number, x4::Number)
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1,x2,x3,x4)
        end
        return y
    end

    # in-place version for scalar x1, x2, abstract array x3, x4
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray,
        x1::Number, x2::Number, x3::AbstractArray, x4::AbstractArray)
        # check the sizes of the inputs
        if size(x3) != size(x4)
            throw(ArgumentError("Sizes of inputs $(size(x3)) and $(size(x4)) must be equal."))
        end
        if size(y)[3:end] != size(x3)
            throw(ArgumentError("Sizes of output $(size(y)) and inputs $(size(x3)) not compatible."))
        end
        # assume a value of x1 is given for each frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1,x2,x3[i],x4[i])
        end
        return y
    end

    # in-place version for abstract array x1, x2, x3, x4
    @eval function $(Symbol(String(f)*"!"))(y::AbstractArray,
        x1::AbstractArray, x2::AbstractArray, x3::AbstractArray,
        x4::AbstractArray)
        # check the sizes of the inputs
        if size(x1) != size(x2) || size(x1) != size(x3) || size(x1) != size(x4)
            throw(ArgumentError("Sizes of inputs $(size(x1)), $(size(x2)), $(size(x3)), and $(size(x4)) must be equal."))
        end
        if size(y)[3:end] != size(x1)
            throw(ArgumentError("Sizes of output $(size(y)) and inputs $(size(x1)) not compatible."))
        end
        # assume a value of x1 is given for each frequency
        # loop over the dimensions of the array greater than 2
        for i in CartesianIndices(axes(y)[3:end])
            $(Symbol(String(f)*"!"))(view(y,:,:,i),x1[i],x2[i],x3[i],x4[i])
        end
        return y
    end
end

"""
    ABCD_seriesZ(Z1)

Return the ABCD matrix for a series impedance `Z1`.
```
o---Z1---o
          
          
o--------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_seriesZ(50)
2×2 Matrix{Int64}:
 1  50
 0   1
```
"""
function ABCD_seriesZ(Z1::Number)
    return [one(Z1) Z1;zero(Z1) one(Z1)]
end

"""
    ABCD_seriesZ!(ABCD,Z1)

In-place version of [`ABCD_seriesZ`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_seriesZ!(zeros(2,2),50)
2×2 Matrix{Float64}:
 1.0  50.0
 0.0   1.0
```
"""
function ABCD_seriesZ!(ABCD::AbstractMatrix,Z1::Number)
    if size(ABCD,1) != 2 || size(ABCD,2) != 2
        throw(ArgumentError("Size of output $(size(ABCD)) must be (2, 2)."))
    end
    ABCD[1,1] = one(eltype(ABCD))
    ABCD[1,2] = Z1
    ABCD[2,1] = zero(eltype(ABCD))
    ABCD[2,2] = one(eltype(ABCD))
    return ABCD
end

"""
    Y_seriesY(Y1)

Return the Y matrix for a series admittance `Y1`.
```
o---Y1---o
          
          
o--------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.Y_seriesY(1/50)
2×2 Matrix{Float64}:
  0.02  -0.02
 -0.02   0.02
```
"""
function Y_seriesY(Y1::Number)
    return [Y1 -Y1;-Y1 Y1]
end

"""
    Y_seriesY!(Y,Y1)

In-place version of [`Y_seriesY`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.Y_seriesY!(zeros(2,2),1/50)
2×2 Matrix{Float64}:
  0.02  -0.02
 -0.02   0.02
```
"""
function Y_seriesY!(Y::AbstractMatrix,Y1::Number)
    if size(Y,1) != 2 || size(Y,2) != 2
        throw(ArgumentError("Size of output $(size(Y)) must be (2, 2)."))
    end
    Y[1,1] = Y1
    Y[1,2] = -Y1
    Y[2,1] = -Y1
    Y[2,2] = Y1
    return Y
end

"""
    ABCD_shuntY(Y1)

Return the ABCD matrix for a shunt admittance `Y1`.
```
o---------o
     |
     Y1
     |
o---------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_shuntY(1/50)
2×2 Matrix{Float64}:
 1.0   0.0
 0.02  1.0
```
"""
function ABCD_shuntY(Y1::Number)
    return [one(Y1) zero(Y1);Y1 one(Y1)]
end

"""
    ABCD_shuntY!(ABCD,Y1)

In-place version of [`ABCD_shuntY`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_shuntY!(zeros(2,2),1/50)
2×2 Matrix{Float64}:
 1.0   0.0
 0.02  1.0
```
"""
function ABCD_shuntY!(ABCD::AbstractMatrix,Y1::Number)
    if size(ABCD,1) != 2 || size(ABCD,2) != 2
        throw(ArgumentError("Size of output $(size(ABCD)) must be (2, 2)."))
    end
    ABCD[1,1] = one(eltype(ABCD))
    ABCD[1,2] = zero(eltype(ABCD))
    ABCD[2,1] = Y1
    ABCD[2,2] = one(eltype(ABCD))
    return ABCD
end

"""
    Z_shuntZ(Z1)

Return the Z matrix for a shunt impedance `Z1`.
```
o---------o
     |
     Z1
     |
o---------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.Z_shuntZ(50)
2×2 Matrix{Int64}:
 50  50
 50  50
```
"""
function Z_shuntZ(Z1::Number)
    return [Z1 Z1;Z1 Z1]
end

"""
    Z_shuntZ!(Z,Z1)

In-place version of [`Z_shuntZ`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.Z_shuntZ!(zeros(2,2),50)
2×2 Matrix{Float64}:
 50.0  50.0
 50.0  50.0
```
"""
function Z_shuntZ!(Z::AbstractMatrix,Z1::Number)
    if size(Z,1) != 2 || size(Z,2) != 2
        throw(ArgumentError("Size of output $(size(Z)) must be (2, 2)."))
    end
    Z[1,1] = Z1
    Z[1,2] = Z1
    Z[2,1] = Z1
    Z[2,2] = Z1
    return Z
end

"""
    ABCD_PiY(Y1,Y2,Y3)

Return the ABCD matrix for a Pi network of admittances `Y1`, `Y2`, and `Y3`.
```
o----Y3-----o
   |     |   
   Y1    Y2  
   |     |   
o-----------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_PiY(1,2,4)
2×2 Matrix{Float64}:
 1.5  0.25
 3.5  1.25
```
"""
function ABCD_PiY(Y1::Number,Y2::Number,Y3::Number)
    return [1+Y2/Y3 1/Y3;Y1+Y2+Y1*Y2/Y3 1+Y1/Y3]
end

"""
    ABCD_PiY!(ABCD,Y1,Y2,Y3)

In-place version of [`ABCD_PiY`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_PiY!(zeros(2,2),1,2,4)
2×2 Matrix{Float64}:
 1.5  0.25
 3.5  1.25
```
"""
function ABCD_PiY!(ABCD::AbstractMatrix,Y1::Number,Y2::Number,Y3::Number)
    if size(ABCD,1) != 2 || size(ABCD,2) != 2
        throw(ArgumentError("Size of output $(size(ABCD)) must be (2, 2)."))
    end
    ABCD[1,1] = 1+Y2/Y3
    ABCD[1,2] = 1/Y3
    ABCD[2,1] = Y1+Y2+Y1*Y2/Y3
    ABCD[2,2] = 1+Y1/Y3
    return ABCD
end

"""
    Y_PiY(Y1,Y2,Y3)

Return the admittance matrix for a Pi network of admittances `Y1`, `Y2`, and
`Y3`.
```
o----Y3-----o
   |     |   
   Y1    Y2  
   |     |   
o-----------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.Y_PiY(1.0,2.0,4.0)
2×2 Matrix{Float64}:
  5.0  -4.0
 -4.0   6.0

julia> Y1=1.0;Y2=2.0;Y3=4.0;isapprox(JosephsonCircuits.YtoA(JosephsonCircuits.Y_PiY(Y1,Y2,Y3)),JosephsonCircuits.ABCD_PiY(Y1,Y2,Y3))
true
```
"""
function Y_PiY(Y1::Number,Y2::Number,Y3::Number)
    return [Y1+Y3 -Y3;-Y3 Y2+Y3]
end

"""
    Y_PiY!(Y,Y1,Y2,Y3)

In-place version of [`Y_PiY`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.Y_PiY!(zeros(2,2),1.0,2.0,4.0)
2×2 Matrix{Float64}:
  5.0  -4.0
 -4.0   6.0
```
"""
function Y_PiY!(Y::AbstractMatrix,Y1::Number,Y2::Number,Y3::Number)
    if size(Y,1) != 2 || size(Y,2) != 2
        throw(ArgumentError("Size of output $(size(Y)) must be (2, 2)."))
    end
    Y[1,1] = Y1+Y3
    Y[1,2] = -Y3
    Y[2,1] = -Y3
    Y[2,2] = Y2+Y3
    return Y
end

"""
    ABCD_TZ(Z1,Z2,Z3)

Return the ABCD matrix for a T network of impedances `Z1`, `Z2`, and `Z3`.
```
o--Z1-----Z2--o
       |       
      Z3       
       |       
o-------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_TZ(1,2,4)
2×2 Matrix{Float64}:
 1.25  3.5
 0.25  1.5
```
"""
function ABCD_TZ(Z1::Number,Z2::Number,Z3::Number)
    return [1+Z1/Z3 Z1+Z2+Z1*Z2/Z3;1/Z3 1+Z2/Z3]
end

"""
    ABCD_TZ!(ABCD,Z1,Z2,Z3)

In-place version of [`ABCD_TZ`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_TZ!(ones(2,2),1,2,4)
2×2 Matrix{Float64}:
 1.25  3.5
 0.25  1.5
```
"""
function ABCD_TZ!(ABCD::AbstractMatrix,Z1::Number,Z2::Number,Z3::Number)
    if size(ABCD,1) != 2 || size(ABCD,2) != 2
        throw(ArgumentError("Size of output $(size(ABCD)) must be (2, 2)."))
    end
    ABCD[1,1] = 1+Z1/Z3
    ABCD[1,2] = Z1+Z2+Z1*Z2/Z3
    ABCD[2,1] = 1/Z3
    ABCD[2,2] = 1+Z2/Z3
    return ABCD
end

"""
    Z_TZ(Z1,Z2,Z3)

Return the ABCD matrix for a T network of impedances `Z1`, `Z2`, and `Z3`.
```
o--Z1-----Z2--o
       |       
      Z3       
       |       
o-------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.Z_TZ(1,2,4)
2×2 Matrix{Int64}:
 5  4
 4  6

julia> Z1=1.0;Z2=2.0;Z3=4.0;isapprox(JosephsonCircuits.ZtoA(JosephsonCircuits.Z_TZ(Z1,Z2,Z3)),JosephsonCircuits.ABCD_TZ(Z1,Z2,Z3))
true
```
"""
function Z_TZ(Z1::Number,Z2::Number,Z3::Number)
    return [Z1+Z3 Z3;Z3 Z2+Z3]
end

"""
    Z_TZ!(Z,Z1,Z2,Z3)

In-place version of [`Z_TZ`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.Z_TZ!(ones(2,2),1,2,4)
2×2 Matrix{Float64}:
 5.0  4.0
 4.0  6.0
```
"""
function Z_TZ!(Z::AbstractMatrix,Z1::Number,Z2::Number,Z3::Number)
    if size(Z,1) != 2 || size(Z,2) != 2
        throw(ArgumentError("Size of output $(size(Z)) must be (2, 2)."))
    end
    Z[1,1] = Z1+Z3
    Z[1,2] = Z3
    Z[2,1] = Z3
    Z[2,2] = Z2+Z3
    return Z
end

"""
    ABCD_tline(Z0, theta)

Return the ABCD matrix for a transmission line described by phase delay
`theta` in radians and characteristic impedance `Z0` in Ohms.
```
   theta, Z0  
o--========--o
              
              
o------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_tline(50, pi/4)
2×2 Matrix{ComplexF64}:
 0.707107+0.0im             0.0+35.3553im
      0.0+0.0141421im  0.707107+0.0im
```
"""
function ABCD_tline(Z0::Number, theta::Number)
    return [cos(theta) im*Z0*sin(theta);im/Z0*sin(theta) cos(theta)]
end

"""
    ABCD_tline!(ABCD, Z0, theta)

In-place version of [`ABCD_tline`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_tline!(ones(Complex{Float64},2,2),50, pi/4)
2×2 Matrix{ComplexF64}:
 0.707107+0.0im             0.0+35.3553im
      0.0+0.0141421im  0.707107+0.0im
```
"""
function ABCD_tline!(ABCD::AbstractMatrix,Z0::Number,theta::Number)
    if size(ABCD,1) != 2 || size(ABCD,2) != 2
        throw(ArgumentError("Size of output $(size(ABCD)) must be (2, 2)."))
    end
    ABCD[1,1] = cos(theta)
    ABCD[1,2] = im*Z0*sin(theta)
    ABCD[2,1] = im/Z0*sin(theta)
    ABCD[2,2] = cos(theta)
    return ABCD
end

"""
    Z_tline(Z0, theta)

Return the impedance matrix for a transmission line described by a
characteristic impedance `Z0` in Ohms and phase delay `theta` in radians 
```
   theta, Z0  
o--========--o
              
              
o------------o
```
# Examples
```jldoctest
julia> JosephsonCircuits.Z_tline(50, pi/4)
2×2 Matrix{ComplexF64}:
 0.0-50.0im     0.0-70.7107im
 0.0-70.7107im  0.0-50.0im

julia> isapprox(JosephsonCircuits.Z_tline(50, pi/4),JosephsonCircuits.AtoZ(JosephsonCircuits.ABCD_tline(50, pi/4)))
true
```
"""
function Z_tline(Z0::Number, theta::Number)
    return [-im*Z0*cot(theta) -im*Z0*csc(theta);-im*Z0*csc(theta) -im*Z0*cot(theta)]
end

"""
    Z_tline!(Z, Z0, theta)

In-place version of [`Z_tline`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.Z_tline!(ones(Complex{Float64},2,2),50, pi/4)
2×2 Matrix{ComplexF64}:
 0.0-50.0im     0.0-70.7107im
 0.0-70.7107im  0.0-50.0im
```
"""
function Z_tline!(Z::AbstractMatrix,Z0::Number,theta::Number)
    if size(Z,1) != 2 || size(Z,2) != 2
        throw(ArgumentError("Size of output $(size(Z)) must be (2, 2)."))
    end
    Z[1,1] = -im*Z0*cot(theta)
    Z[1,2] = -im*Z0*csc(theta)
    Z[2,1] = -im*Z0*csc(theta)
    Z[2,2] = -im*Z0*cot(theta)
    return Z
end

"""
    ABCD_coupled_tline(Z0e, Z0o, thetae, thetao)

Return the ABCD matrix for two coupled transmission lines described by
even and odd mode phase delays `thetae` and `thetao`, and even and odd mode
impedances `Z0e` and `Z0o`.
```
thetae, Z0e
thetao, Z0o

V1, I1 -->  ======== <-- I3, V3
V2, I2 -->  ======== <-- I4, V4

[(V1+V2)/2, (I1+I2)/2] = ABCDe * [(V3+V4)/2, -(I3+I4)/2]
[(V1-V2)/2, (I1-I2)/2] = ABCDo * [(V3-V4)/2, -(I3-I4)/2]

[V1, V2, I1, I2] = ABCD_coupled_tline * [V3, V4, -I3, -I4]
```
# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_coupled_tline(50,50,pi/4,pi/4)
4×4 Matrix{ComplexF64}:
 0.707107+0.0im             0.0+0.0im        …       0.0+0.0im
      0.0+0.0im        0.707107+0.0im                0.0+35.3553im
      0.0+0.0141421im       0.0+0.0im                0.0+0.0im
      0.0+0.0im             0.0+0.0141421im     0.707107+0.0im

julia> isapprox(JosephsonCircuits.AtoZ(JosephsonCircuits.ABCD_coupled_tline(55,50,pi/4,pi/4)),JosephsonCircuits.Z_coupled_tline(55,50,pi/4,pi/4))
true
```
"""
function ABCD_coupled_tline(Z0e::Number, Z0o::Number, thetae::Number,
    thetao::Number)
    A11 = A22 = A33 = A44 = 1/2*(cos(thetae) + cos(thetao))
    A12 = A21 = A34 = A43 = 1/2*(cos(thetae) - cos(thetao))
    A13 = A24 = im/2*(Z0e*sin(thetae) + Z0o*sin(thetao))
    A31 = A42 = im/2*(sin(thetae)/Z0e + sin(thetao)/Z0o)
    A14 = A23 = im/2*(Z0e*sin(thetae) - Z0o*sin(thetao))
    A41 = A32 = im/2*(sin(thetae)/Z0e - sin(thetao)/Z0o)
    return [A11 A12 A13 A14; A21 A22 A23 A24; A31 A32 A33 A34; A41 A42 A43 A44]
end

"""
    ABCD_coupled_tline!(A, Z0e, Z0o, thetae, thetao)

In-place version of [`ABCD_coupled_tline`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.ABCD_coupled_tline!(zeros(Complex{Float64},4,4),50,50,pi/4,pi/4)
4×4 Matrix{ComplexF64}:
 0.707107+0.0im             0.0+0.0im        …       0.0+0.0im
      0.0+0.0im        0.707107+0.0im                0.0+35.3553im
      0.0+0.0141421im       0.0+0.0im                0.0+0.0im
      0.0+0.0im             0.0+0.0141421im     0.707107+0.0im
```
"""
function ABCD_coupled_tline!(A::AbstractMatrix, Z0e::Number, Z0o::Number,
    thetae::Number, thetao::Number)
    if size(A,1) != 4 || size(A,2) != 4
        throw(ArgumentError("Size of output $(size(A)) must be (4, 4)."))
    end
    A[1,1] = A[2,2] = A[3,3] = A[4,4] = 1/2*(cos(thetae) + cos(thetao))
    A[1,2] = A[2,1] = A[3,4] = A[4,3] = 1/2*(cos(thetae) - cos(thetao))
    A[1,3] = A[2,4] = im/2*(Z0e*sin(thetae) + Z0o*sin(thetao))
    A[3,1] = A[4,2] = im/2*(sin(thetae)/Z0e + sin(thetao)/Z0o)
    A[1,4] = A[2,3] = im/2*(Z0e*sin(thetae) - Z0o*sin(thetao))
    A[4,1] = A[3,2] = im/2*(sin(thetae)/Z0e - sin(thetao)/Z0o)
    return A
end

"""
    Z_coupled_tline(Z0e, Z0o, thetae, thetao)

Return the impedance matrix for two coupled transmission lines described even
and odd mode impedances `Z0e` and `Z0o` and by even and odd mode phase delays
`thetae` and `thetao`.
```
thetae, Z0e
thetao, Z0o

V1, I1 -->  ======== <-- I3, V3
V2, I2 -->  ======== <-- I4, V4

[(V1+V2)/2, (I1+I2)/2] = ABCDe * [(V3+V4)/2, -(I3+I4)/2]
[(V1-V2)/2, (I1-I2)/2] = ABCDo * [(V3-V4)/2, -(I3-I4)/2]

[V1, V2, V3, V4] = Z_coupled_tline * [I1, I2, I3, I4]
```
# Examples
```jldoctest
julia> JosephsonCircuits.Z_coupled_tline(50,50,pi/4,pi/4)
4×4 Matrix{ComplexF64}:
 0.0-50.0im     0.0-0.0im      0.0-70.7107im  0.0-0.0im
 0.0-0.0im      0.0-50.0im     0.0-0.0im      0.0-70.7107im
 0.0-70.7107im  0.0-0.0im      0.0-50.0im     0.0-0.0im
 0.0-0.0im      0.0-70.7107im  0.0-0.0im      0.0-50.0im
```
"""
function Z_coupled_tline(Z0e::Number, Z0o::Number, thetae::Number,
    thetao::Number)
    Z11 = Z22 = Z33 = Z44 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
    Z12 = Z21 = Z34 = Z43 = -im/2*(Z0e*cot(thetae) - Z0o*cot(thetao)) 
    Z13 = Z31 = Z42 = Z24 = -im/2*(Z0e*csc(thetae) + Z0o*csc(thetao))
    Z14 = Z41 = Z32 = Z23 = -im/2*(Z0e*csc(thetae) - Z0o*csc(thetao))
    return [Z11 Z12 Z13 Z14; Z21 Z22 Z23 Z24; Z31 Z32 Z33 Z34; Z41 Z42 Z43 Z44]
end

"""
    Z_coupled_tline!(Z, Z0e, Z0o, thetae, thetao)

In-place version of [`Z_coupled_tline`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.Z_coupled_tline!(zeros(Complex{Float64},4,4),50,50,pi/4,pi/4)
4×4 Matrix{ComplexF64}:
 0.0-50.0im     0.0-0.0im      0.0-70.7107im  0.0-0.0im
 0.0-0.0im      0.0-50.0im     0.0-0.0im      0.0-70.7107im
 0.0-70.7107im  0.0-0.0im      0.0-50.0im     0.0-0.0im
 0.0-0.0im      0.0-70.7107im  0.0-0.0im      0.0-50.0im
```
"""
function Z_coupled_tline!(Z::AbstractMatrix, Z0e::Number, Z0o::Number,
    thetae::Number, thetao::Number)
    if size(Z,1) != 4 || size(Z,2) != 4
        throw(ArgumentError("Size of output $(size(Z)) must be (4, 4)."))
    end
    Z[1,1] = Z[2,2] = Z[3,3] = Z[4,4] = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
    Z[1,2] = Z[2,1] = Z[3,4] = Z[4,3] = -im/2*(Z0e*cot(thetae) - Z0o*cot(thetao)) 
    Z[1,3] = Z[3,1] = Z[4,2] = Z[2,4] = -im/2*(Z0e*csc(thetae) + Z0o*csc(thetao))
    Z[1,4] = Z[4,1] = Z[3,2] = Z[2,3] = -im/2*(Z0e*csc(thetae) - Z0o*csc(thetao))
    return Z
end

"""
    A_coupled_tlines(L,Cmaxwell,l,omega)

Returns the 2mx2m chain (ABCD) matrix for a port number symmetric multi-port
network of m coupled transmission lines described by a symmetric mxm
Maxwell inductance (per unit length) matrix `L`, a symmetric mxm Maxwell
capacitance matrix (per unit length) `Cmaxwell`, a physical length `l`, and an
angular frequency `omega`.
```
V_1, I_1 -->  ======== <-- I_{m+1}, V_{m+1}
V_2, I_2 -->  ======== <-- I_{m+2}, V_{m+2}
          .
          .
          .
V_m, I_m -->  ======== <-- I_n, V_n
              <---l-->

where n=2*m.

[V_1, ...V_m, I_1, ...I_m] = A_coupled_tline * [V_{m+1}, ...V_n, I_{m+1}, ...I_n]
```

# Examples
```jldoctest
Zeven = 51.0
Zodd = 49.0
neven = 1.1
nodd = 1.08
l = 3.5e-3
c = JosephsonCircuits.speed_of_light
omega = 2*pi*5e9

L, C = JosephsonCircuits.even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
A1 = JosephsonCircuits.A_coupled_tlines(L,C,l,omega)
A2 = JosephsonCircuits.ABCD_coupled_tline(Zeven,Zodd,neven*omega/c*l,nodd*omega/c*l)
println(isapprox(A1,A2))

# output
true
```

# References
Paul, Clayton R. Analysis of Multiconductor Transmission Lines, Second
Edition. Wiley, 2008.
"""
function A_coupled_tlines(L,Cmaxwell,l,omega)

    N = size(L,1)

    # compute the basis for the coupled lines
    b = ZC_basis_coupled_tlines(L,Cmaxwell)
    TI = b.TI
    TV = b.TV
    ZC = b.ZC
    lambda = b.lambda

    TIinv = inv(TI)
    TVinv = inv(TV)
    YC = inv(ZC)

    # pre-compute these matrix products. they are the same for every frequency
    ZC_TI = ZC*TI
    TIinv_YC = TIinv*YC

    # compute the first gamma*l
    gammal = im*lambda*first(omega)*l

    # compute the first cosh(gamma*l) and sinh(gamma*l). reuse these Diagonal
    # matrices for the rest of the frequencies.
    coshgammal = cosh(gammal)
    sinhgammal = sinh(gammal)

    # allocate a temporary array
    phi_tmp = zeros(Complex{Float64},N,N)

    # define the output matrix
    A = zeros(Complex{Float64},ifelse(length(omega)>1,(2*N,2*N,length(omega)),(2*N,2*N)))

    # loop over the frequencies
    for i in eachindex(omega)
        # update coshgammal and sinhgammal
        if i > 1
            gammal .= im*lambda*omega[i]*l
            coshgammal .= cosh(gammal)
            sinhgammal .= sinh(gammal)
        end

        # make views for the sub-matrices
        phi11 = view(A,1:N,1:N,i)
        phi12 = view(A,1:N,N+1:2*N,i)
        phi21 = view(A,N+1:2*N,1:N,i)
        phi22 = view(A,N+1:2*N,N+1:2*N,i)
        
        # compute the ABCD matrix for the coupled lines
        A_coupled_tlines!(phi11,phi12,phi21,phi22,phi_tmp,TI,TIinv,ZC_TI,TIinv_YC,coshgammal,sinhgammal)
    end
    return A
end


function A_coupled_tlines!(phi11, phi12, phi21, phi22, phi_tmp, TI, TIinv,
    ZC_TI, TIinv_YC, coshgammal, sinhgammal)

    # # The equations
    # phi11 = 1/2*ZC*TI*(exp(gamma*l)+exp(-gamma*l))*TIinv*YC
    # phi12 = -1/2*ZC*TI*(exp(-gamma*l)-exp(gamma*l))*TIinv
    # phi21 = -1/2*TI*(exp(-gamma*l)-exp(gamma*l))*TIinv*YC
    # phi22 = 1/2*TI*(exp(gamma*l)+exp(-gamma*l))*TIinv

    # phi11 = 1/2*TV*(exp(gamma*l)+exp(-gamma*l))*TVinv
    # phi12 = -1/2*TV*(exp(gamma*l)-exp(-gamma*l))*TVinv*ZC
    # phi21 = -1/2*YC*TV*(exp(gamma*l)-exp(-gamma*l))*TVinv
    # phi22 = 1/2*YC*TV*(exp(gamma*l)+exp(-gamma*l))*TVinv*ZC

    mul!(phi_tmp,coshgammal,TIinv_YC)
    mul!(phi11,ZC_TI,phi_tmp)

    mul!(phi_tmp,sinhgammal,TIinv)
    mul!(phi12,ZC_TI,phi_tmp)

    mul!(phi_tmp,sinhgammal,TIinv_YC)
    mul!(phi21,TI,phi_tmp)

    mul!(phi_tmp,coshgammal,TIinv)
    mul!(phi22,TI,phi_tmp)
end

"""
    ZC_basis_coupled_tlines(L, Cmaxwell)

    Returns the characteristic impedance matrix `ZC` and eigenbasis for
    current `TI` and voltage `TV` from the inductance per unit length matrix
    `L` and Maxwell capacitance per unit length matrix `Cmaxwell`.

# Arguments
- `L`: inductance per unit length matrix.
- `C`: Maxwell capacitance per unit length matrix.

# Returns
- `ZC`: characteristic impedance matrix.
- `TI`: matrix which transforms mode currents to currents, I = TI*Im. Computed
        from TI = U*theta*S.
- `TV`: matrix which transforms mode voltages to voltages, V = TV*Vm. Computed
        from TV = U*inv(theta)*S.
- `theta`: Diagonal matrix with the square of the eigenvalues of Cmaxwell
        along the diagonals.
- `U`: eigenvectors of Cmaxwell.
- `lambda`: Diagonal matrix with the square of the eigenvalues of
        theta*Ut*L*U*theta along the diagonals. `lambda` is related
        to the propagation constant, gamma, as gamma^2 = -omega^2*lambda^2.
- `S`: eigenvectors of theta*Ut*L*U*theta.

# Examples
```jldoctest
Zeven = 51.0
Zodd = 49.0
neven = 1.1
nodd = 1.08
c = JosephsonCircuits.speed_of_light

L, C = JosephsonCircuits.even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
b = JosephsonCircuits.ZC_basis_coupled_tlines(L,C)
@show b.ZC
@show b.TI
@show b.TV
@show Matrix(b.theta)
@show b.U
@show Matrix(b.lambda)
@show b.S
println(isapprox(Zeven,b.ZC[1,1]+b.ZC[1,2]))
println(isapprox(Zodd,b.ZC[1,1]-b.ZC[1,2]))
println(isapprox(neven,b.lambda[2,2]*c))
println(isapprox(nodd,b.lambda[1,1]*c))

# output
b.ZC = [49.999999999999986 0.9999999999999929; 0.9999999999999929 49.999999999999986]
b.TI = [-6.063012846509498e-6 -5.997716107132906e-6; 6.063012846509498e-6 -5.997716107132906e-6]
b.TV = [-82467.25063230767 -83365.06614665617; 82467.25063230767 -83365.06614665617]
Matrix(b.theta) = [8.48205146197092e-6 0.0; 0.0 8.574394996376037e-6]
b.U = [-0.7071067811865475 -0.7071067811865475; -0.7071067811865475 0.7071067811865475]
Matrix(b.lambda) = [3.6024922281400425e-9 0.0; 0.0 3.669205047179673e-9]
b.S = [0.0 1.0; 1.0 0.0]
true
true
true
true
```

# References
Paul, Clayton R. Analysis of Multiconductor Transmission Lines, Second
Edition. Wiley, 2008.
"""
function ZC_basis_coupled_tlines(L,Cmaxwell)

    # compute U and theta
    vals, vecs = eigen(Cmaxwell)
    U = vecs
    theta = Diagonal(sqrt.(vals))

    # compute S and lambda
    vals, vecs = eigen(theta*U'*L*U*theta)
    lambda = Diagonal(sqrt.(vals))
    S = vecs

    # compute TI and TV
    TI = U*theta*S
    TV = U*inv(theta)*S

    # compute ZC and YC
    ZC = U*inv(theta)*S*lambda*S'*inv(theta)*U'

    return CoupledLinesBasis(ZC, TI, TV, theta, U, lambda,S)
end

"""
    CoupledLinesBasis(ZC, TI, TV, theta, U, lambda, S)

A simple structure to hold the output of `ZC_basis_coupled_tlines`.
"""
struct CoupledLinesBasis
    ZC
    TI
    TV
    theta
    U
    lambda
    S
end

"""
    maxwell_to_mutual(Cmaxwell::AbstractMatrix)

Return the mutual capacitance matrix from the Maxwell capacitance matrix
`Cmaxwell`.

The Maxwell capacitance `Cmaxwell` is the relationship between charge and
voltage on each node, Q = C V or dQi/dVj = C_ij where `C` is the
Maxwell capacitance matrix.

Each element of the mutual capacitance matrix `Cmutual` is the value of a
physical capacitor placed between two nodes in a circuit or between a node
and ground.

# Examples
```jldoctest
julia> @variables C11, C12, C21, C22;C = [C11 C12;C21 C22];JosephsonCircuits.maxwell_to_mutual(C)
2×2 Matrix{Num}:
 C11 + C12       -C12
      -C21  C21 + C22

julia> C = [1.0 -0.1;-0.1 2.0];JosephsonCircuits.maxwell_to_mutual(C)
2×2 Matrix{Float64}:
 0.9  0.1
 0.1  1.9
```
"""
function maxwell_to_mutual(Cmaxwell::AbstractMatrix)
  return Diagonal(Cmaxwell)-Cmaxwell + Diagonal(sum(Cmaxwell,dims=2)[:])
end

"""
    maxwell_to_mutual(Cmutual::AbstractMatrix)

Return the Maxwell capacitance matrix from the mutual capacitance matrix
`Cmutual`.

The Maxwell capacitance `Cmaxwell` is the relationship between charge and
voltage on each node, Q = C V or dQi/dVj = C_ij where `C` is the
Maxwell capacitance matrix.

Each element of the mutual capacitance matrix `Cmutual` is the value of a
physical capacitor placed between two nodes in a circuit or between a node
and ground.

# Examples
```jldoctest
julia> @variables Cg, Cm;C = [Cg Cm;Cm Cg];JosephsonCircuits.mutual_to_maxwell(C)
2×2 Matrix{Num}:
 Cg + Cm      -Cm
     -Cm  Cg + Cm

julia> C = [0.9 0.1;0.1 1.9];JosephsonCircuits.mutual_to_maxwell(C)
2×2 Matrix{Float64}:
  1.0  -0.1
 -0.1   2.0
```
"""
function mutual_to_maxwell(C::AbstractMatrix)
  return maxwell_to_mutual(C)
end

"""
    maxwell_to_even_odd(L, Cmaxwell)

Return the even and odd mode impedances and the even and odd mode indices from
the inductance matrix `L` and the Maxwell capacitance matrix `Cmaxwell`.

# Examples
```jldoctest
@variables C11, C12, L11, L12
C = [C11 C12;C12 C11]
L = [L11 L12;L12 L11]
Zeven, Zodd, neven, nodd = JosephsonCircuits.maxwell_to_even_odd(L,C)
@show Zeven
@show Zodd
@show neven
@show nodd
;

# output
Zeven = sqrt((L11 + L12) / (C11 + C12))
Zodd = sqrt((L11 - L12) / (C11 - C12))
neven = 2.99792458e8sqrt((C11 + C12)*(L11 + L12))
nodd = 2.99792458e8sqrt((C11 - C12)*(L11 - L12))
```
"""
function maxwell_to_even_odd(L, Cmaxwell)

    # consider erroring if not 2x2 matrices with capacitance and inductance
    # to ground equal for both transmission lines.
    # also error if not symmetric
    # (Cmaxwell[1,1] != Cmaxwell[2,2]) || (Cmaxwell[1,2] != Cmaxwell[2,1])
    # size(Cmaxwell) != (2,2)
    # same for inductance matrix

    c = JosephsonCircuits.speed_of_light
    Zeven = sqrt((L[1,1]+L[1,2])/(Cmaxwell[1,1]+Cmaxwell[1,2]))
    Zodd = sqrt((L[1,1]-L[1,2])/(Cmaxwell[1,1]-Cmaxwell[1,2]))
    neven = c*sqrt((L[1,1]+L[1,2])*(Cmaxwell[1,1]+Cmaxwell[1,2]))
    nodd = c*sqrt((L[1,1]-L[1,2])*(Cmaxwell[1,1]-Cmaxwell[1,2]))
    return (Zeven = Zeven, Zodd = Zodd, neven = neven, nodd = nodd)
end

"""
    mutual_to_even_odd(L, Cmutual)

Return the even and odd mode impedances and the even and odd mode indices from
the inductance matrix `L` and the mutual capacitance matrix `Cmutual`.

# Examples
```jldoctest
@variables Cg, Cm, Ls, Lm
C = [Cg Cm; Cm Cg]
L = [Ls Lm; Lm Ls]
Zeven, Zodd, neven, nodd = JosephsonCircuits.mutual_to_even_odd(L,C)
@show Zeven
@show Zodd
@show neven
@show nodd
;

# output
Zeven = sqrt((Lm + Ls) / Cg)
Zodd = sqrt((-Lm + Ls) / (Cg + 2Cm))
neven = 2.99792458e8sqrt(Cg*(Lm + Ls))
nodd = 2.99792458e8sqrt((Cg + 2Cm)*(-Lm + Ls))
```
"""
function mutual_to_even_odd(L, Cmutual)
    return maxwell_to_even_odd(L, maxwell_to_mutual(Cmutual))
end


"""
    even_odd_to_maxwell(Zeven, Zodd, neven, nodd)

Return the inductance matrix and Maxwell capacitance matrix for two coupled
transmission lines with even and odd mode impedances `Zeven`, `Zodd` and
even and odd mode indices `neven`, `nodd`.

# Examples
```jldoctest
L1 = [1.1 0.1;0.1 1.1]
C1 = [2.0 -0.4;-0.4 2.0]
L2, C2 = JosephsonCircuits.even_odd_to_maxwell(JosephsonCircuits.maxwell_to_even_odd(L1,C1)...)
isapprox(L1,L2) && isapprox(C1,C2)

# output
true
```
"""
function even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
    c = JosephsonCircuits.speed_of_light
    L11 = (neven*Zeven+nodd*Zodd)/(2*c)
    L12 = (neven*Zeven-nodd*Zodd)/(2*c)
    C11 = (neven*Zodd+nodd*Zeven)/(2*c*Zeven*Zodd)
    C12 = (neven*Zodd-nodd*Zeven)/(2*c*Zodd*Zeven)
    return (L = [L11 L12;L12 L11], C = [C11 C12;C12 C11])
end

"""
    even_odd_to_mutual(Zeven, Zodd, neven, nodd)

Return the inductance matrix and mutual capacitance matrix for two coupled
transmission lines with even and odd mode impedances `Zeven`, `Zodd` and
even and odd mode indices `neven`, `nodd`.

# Examples
```jldoctest
L1 = [1.1 0.1;0.1 1.1]
C1 = [1.6 0.4;0.4 1.6]
L2, C2 = JosephsonCircuits.even_odd_to_mutual(JosephsonCircuits.mutual_to_even_odd(L1,C1)...)
isapprox(L1,L2) && isapprox(C1,C2)

# output
true
```
"""
function even_odd_to_mutual(Zeven, Zodd, neven, nodd)
    L, Cmaxwell = even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
    return (L = L, C = maxwell_to_mutual(Cmaxwell))
end


"""
    Z_canonical_coupled_line_circuit(i::Int,  Z0e, Z0o, thetae, thetao)

Return the impedance matrix for the `i`'th canonical coupled line circuit, as
a function of the even mode phase delay `thetae` in radians, the odd mode
phase delay `thetao` in radians, the even mode characteristic impedance `Z0e`
in Ohms, and the odd mode characteristic impedance `Z0o` in Ohms.
1) low pass
```
   gnd--==========
1--> o--==========--o <--2
```
2) band pass
```
   gnd--==========--o <--2
1--> o--==========--gnd
```
3) band pass
```
        ==========--o <--2
1--> o--==========
```
4) band pass
```
1--> o--==========--gnd
2--> o--==========
```
5) all pass
```
        ==========
1--> o--==========--o  <--2
```
6) all pass
```
   gnd--==========--gnd
1--> o--==========--o  <--2
```
7) all pass
```
1--> o--==========--|
2--> o--==========--|
```
8) all stop
```
   gnd--==========--o  <--2
1--> o--==========
```
9) all stop
```
1--> o--==========--gnd
2--> o--==========--gnd
```
10) all stop
```
1--> o--==========
2--> o--==========
```
# Examples
```jldoctest
julia> @variables θe, θo, Ze, Zo;JosephsonCircuits.Z_canonical_coupled_line_circuits(3,Ze,Zo,θe,θo)
2×2 Matrix{Complex{Num}}:
 -0.5(Ze*cot(θe) + Zo*cot(θo))*im  -0.5(Ze*csc(θe) - Zo*csc(θo))*im
 -0.5(Ze*csc(θe) - Zo*csc(θo))*im  -0.5(Ze*cot(θe) + Zo*cot(θo))*im
```
# References
E. M. T. Jones, "Coupled-Strip-Transmission-Line Filters and Directional
Couplers," in IRE Transactions on Microwave Theory and Techniques, vol. 4,
no. 2, pp. 75-81, April 1956, doi: 10.1109/TMTT.1956.1125022.

Pozar, D. M. Microwave Engineering (4 ed.). John Wiley & Sons (2011)
ISBN 9780470631553.
"""
function Z_canonical_coupled_line_circuits(i::Int, Z0e, Z0o, thetae, thetao)

    if i == 1
        Z11 = -2*im*Z0e*Z0o*cos(thetae)*cos(thetao)/(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao))
        Z22 = im*(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))/(2*(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao)))
        Z12 = Z21 = -im*Z0e*Z0o*(cos(thetae)+cos(thetao))/(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao))
    elseif i == 2
        Z11 = Z22 = -2*im*Z0e*Z0o*(Z0e*cos(thetao)*sin(thetae)+Z0o*cos(thetae)*sin(thetao))/(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
        Z12 = Z21 = 2*im*Z0e*Z0o*(Z0e*sin(thetae)-Z0o*sin(thetao))/(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
    elseif i == 3
        Z11 = Z22 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
        Z12 = Z21 = -im/2*(Z0e*csc(thetae) - Z0o*csc(thetao))
    elseif i == 4
        Z11 = im*(2*Z0e*Z0o-2*Z0e*Z0o*cos(thetae)*cos(thetao)+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))/(2*(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao)))
        Z22 = im*(-2*Z0e*Z0o*(1+cos(thetae)*cos(thetao))+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))/(2*(Z0o*cos(thetao)*sin(thetae)+Z0e*cos(thetae)*sin(thetao)))
        Z12 = Z21 = im*(Z0e-Z0o)*(Z0e+Z0o)/(2*(Z0e*cot(thetae)+Z0o*cot(thetao)))
    elseif i == 5
        Z11 = Z22 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
        Z12 = Z21 = -im/2*(Z0e*csc(thetae) + Z0o*csc(thetao))
    elseif i == 6
        Z11 = Z22 = -2*im*Z0e*Z0o*(Z0e*cos(thetao)*sin(thetae)+Z0o*cos(thetae)*sin(thetao))/(2*Z0e*Z0o-2*Z0e*Z0o*cos(thetae)*cos(thetao)+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
        Z12 = Z21 = -2*im*Z0e*Z0o*(Z0e*sin(thetae)+Z0o*sin(thetao))/(2*Z0e*Z0o-2*Z0e*Z0o*cos(thetae)*cos(thetao)+(Z0e^2+Z0o^2)*sin(thetae)*sin(thetao))
    elseif i == 7
        Z11 = Z22 = -im/2*(Z0e*cot(thetae)-Z0o*tan(thetao))
        Z21 = Z12 = -im/2*(Z0e*cot(thetae)+Z0o*tan(thetao))
    elseif i == 8
        Z11 = -2*im*Z0e*Z0o*cot(thetae)*cot(thetao)/(Z0e*cot(thetae)+Z0o*cot(thetao))
        Z22 = -im*(Z0e*cot(thetae) + Z0o*cot(thetao) - Z0e*csc(thetae)-Z0o*csc(thetao))*(Z0e*cot(thetae)+Z0e*csc(thetae)+Z0o*cot(thetao)+Z0o*csc(thetao))/(2*(Z0e*cot(thetae)+Z0o*cot(thetao)))
        Z12 = Z21 = im*Z0e*Z0o*(cos(thetae)-cos(thetao))*csc(thetae)*csc(thetao)/(Z0e*cot(thetae)+Z0o*cot(thetao))
    elseif i == 9
        Z11 = Z22 = im/2*(Z0e*tan(thetae)+Z0o*tan(thetao))
        Z12 = Z21 = im/2*(Z0e*tan(thetae)-Z0o*tan(thetao))
    elseif i == 10
        Z11 = Z22 = -im/2*(Z0e*cot(thetae) + Z0o*cot(thetao))
        Z12 = Z21 = -im/2*(Z0e*cot(thetae) - Z0o*cot(thetao))
    else
        throw(ArgumentError("Canonical coupled line circuit number must be 1-10."))
    end

    return [Z11 Z12; Z21 Z22]
end

function canonical_coupled_line_circuits(i::Int, Z0e, Z0o, ne, no)
    c = JosephsonCircuits.speed_of_light

    if i == 3
        L1 = L2 = (ne*Z0e+no*Z0o)/(6*c)
        M = -(ne*Z0e-no*Z0o)/(12*c)
        C1 = C2 = ne/(c*Z0e)
        Cm = (no*Z0e-ne*Z0o)/(2*c*Z0e*Z0o)
        return (L1 = L1, L2 = L2, M = M, C1 = C1, C2 = C2, Cm = Cm)
    elseif i == 8
        L1 = 2*Z0e*Z0o*(no^3*Z0e+ne^3*Z0o)/(3*c*(Z0e*no+Z0o*ne)^2)
        L2 = (ne*Z0e+no*Z0o)/(2*c)
        M = -Z0e*Z0o*(ne^2-no^2)/(2*c*(no*Z0e+ne*Z0o))
        C1 = (no*Z0e+ne*Z0o)/(2*c*Z0e*Z0o)
        return (L1 = L1, L2 = L2, M = M, C1 = C1)
    elseif i == 9
        L1 = L2 = (ne*Z0e+no*Z0o)/(2*c)
        M = (ne*Z0e-no*Z0o)/(2*c)
        return (L1 = L1, L2 = L2, M = M)
    elseif i == 10
        L1 = L2 = (ne*Z0e+no*Z0o)/(6*c)
        M = (ne*Z0e-no*Z0o)/(6*c)
        C1 = C2 = ne/(c*Z0e)
        Cm = (no*Z0e-ne*Z0o)/(2*c*Z0e*Z0o)
        return (L1 = L1, L2 = L2, M = M, C1 = C1, C2 = C2, Cm = Cm)
    else
        throw(ArgumentError("Canonical coupled line circuit number must be 1-10."))
    end
end

"""
    maxwell_combine(n::Int, d::Dict{NTuple{N, Int}, T}) where {N,T<:AbstractMatrix}

Return the Maxwell capacitance matrix for an `n` terminal system from the
Maxwell capacitance matrices for sets of terminals stored in the dictionary
`d`. The dictionary keys are tuples of the terminal numbers for the
capacitance matrices and the values are the capacitance matrices.

# Examples
```jldoctest
julia> @variables C11, C12, C13, C21, C22, C23, C31, C32, C33;JosephsonCircuits.maxwell_combine(3, Dict((1,2)=>[C11 C12;C21 C22],(1,3)=>[C11 C13;C31 C33],(2,3)=>[C22 C23;C32 C33]))
3×3 Matrix{Num}:
 C11  C12  C13
 C21  C22  C23
 C31  C32  C33

julia> @variables C11, C12, C13, C21, C22, C23, C31, C32, C33;JosephsonCircuits.maxwell_combine(3, Dict((1,2,3)=>[C11 C12 C13;C21 C22 C23; C31 C32 C33]))
3×3 Matrix{Num}:
 C11  C12  C13
 C21  C22  C23
 C31  C32  C33
```
"""
function maxwell_combine(n::Int, d::Dict{NTuple{N, Int}, T}) where {N,T<:AbstractMatrix}

    C = zeros(eltype(T), n, n)
    for (key,val) in d
        for i in 1:N
            for j in 1:N
                if iszero(C[key[i],key[j]])
                    C[key[i],key[j]] = val[i,j]
                else
                    C[key[i],key[j]] += val[i,j]
                    C[key[i],key[j]] /= 2
                end
            end
        end
    end
    return C
end


"""
    S_splitter!(S::AbstractArray)

Return the scattering parameters for a N port ideal lossless symmetrical
reciprocal network. Overwrite `S` with the output.

# Examples
```jldoctest
julia> JosephsonCircuits.S_splitter!(ones(2,2))
2×2 Matrix{Float64}:
 0.0  1.0
 1.0  0.0

julia> JosephsonCircuits.S_splitter!(ones(3,3))
3×3 Matrix{Float64}:
 -0.333333   0.666667   0.666667
  0.666667  -0.333333   0.666667
  0.666667   0.666667  -0.333333
```
"""
function S_splitter!(S::AbstractArray)

    # scattering matrices should be square.
    if size(S,1) != size(S,2)
        throw(ArgumentError("The sizes of the first two dimensions ($(size(S,1)),$(size(S,2))) of the scattering matrix `S` must be the same."))
    end
    # scattering matrices need to be matrices or higher dimensional arrays,
    # not vectors
    if ndims(S) < 2
        throw(ArgumentError("The scattering matrix `S` with size $(size(S)) must have at least two dimensions."))
    end

    # fill with 2/N where N is the size of the scattering matrix.
    fill!(S,2/size(S,1))
  
    # loop over the dimensions of the array greater than 2 and subtract one
    # from the diagonals.
    for k in CartesianIndices(axes(S)[3:end])
        for i in 1:size(S,1)
            S[i,i,k] -= 1
        end
    end

    return S
end


"""
    S_short!(S::AbstractArray)

Return the scattering parameters for a N port ideal short. Overwrite`S` with
the output.

# Examples
```jldoctest
julia> JosephsonCircuits.S_short!(ones(1,1))
1×1 Matrix{Float64}:
 -1.0

julia> JosephsonCircuits.S_short!(ones(2,2))
2×2 Matrix{Float64}:
 -1.0   0.0
  0.0  -1.0
```
"""
function S_short!(S::AbstractArray)

    # scattering matrices should be square.
    if size(S,1) != size(S,2)
        throw(ArgumentError("The sizes of the first two dimensions ($(size(S,1)),$(size(S,2))) of the scattering matrix `S` must be the same."))
    end
    # scattering matrices need to be matrices or higher dimensional arrays,
    # not vectors
    if ndims(S) < 2
        throw(ArgumentError("The scattering matrix `S` with size $(size(S)) must have at least two dimensions."))
    end

    # fill with zeros
    fill!(S,zero(eltype(S)))
  
    # loop over the dimensions of the array greater than 2 and set the diagonals equal to 1
    for k in CartesianIndices(axes(S)[3:end])
        for i in 1:size(S,1)
            S[i,i,k] = -one(eltype(S))
        end
    end

    return S
end

"""
    S_open!(S::AbstractArray)

Return the scattering parameters for a N port ideal open. Overwrite`S` with
the output.

# Examples
```jldoctest
julia> JosephsonCircuits.S_open!(ones(1,1))
1×1 Matrix{Float64}:
 1.0

julia> JosephsonCircuits.S_open!(ones(2,2))
2×2 Matrix{Float64}:
 1.0  0.0
 0.0  1.0
```
"""
function S_open!(S::AbstractArray)

    # scattering matrices should be square.
    if size(S,1) != size(S,2)
        throw(ArgumentError("The sizes of the first two dimensions ($(size(S,1)),$(size(S,2))) of the scattering matrix `S` must be the same."))
    end
    # scattering matrices need to be matrices or higher dimensional arrays,
    # not vectors
    if ndims(S) < 2
        throw(ArgumentError("The scattering matrix `S` with size $(size(S)) must have at least two dimensions."))
    end

    # fill with zeros
    fill!(S,zero(eltype(S)))
  
    # loop over the dimensions of the array greater than 2 and set the diagonals equal to 1
    for k in CartesianIndices(axes(S)[3:end])
        for i in 1:size(S,1)
            S[i,i,k] = one(eltype(S))
        end
    end

    return S
end

"""
    S_match!(S::AbstractArray)

Return the scattering parameters for a N port ideal match. Overwrite`S` with
the output.

# Examples
```jldoctest
julia> JosephsonCircuits.S_match!(ones(1,1))
1×1 Matrix{Float64}:
 0.0

julia> JosephsonCircuits.S_match!(ones(1,1,2))
1×1×2 Array{Float64, 3}:
[:, :, 1] =
 0.0

[:, :, 2] =
 0.0
```
"""
function S_match!(S::AbstractArray)

    # scattering matrices should be square.
    if size(S,1) != size(S,2)
        throw(ArgumentError("The sizes of the first two dimensions ($(size(S,1)),$(size(S,2))) of the scattering matrix `S` must be the same."))
    end
    # scattering matrices need to be matrices or higher dimensional arrays,
    # not vectors
    if ndims(S) < 2
        throw(ArgumentError("The scattering matrix `S` with size $(size(S)) must have at least two dimensions."))
    end

    # fill with zeros
    fill!(S,zero(eltype(S)))

    return S
end
