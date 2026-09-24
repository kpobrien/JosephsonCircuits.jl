# MIT license
# Copyright (c) 2012-2025 DSP.jl contributors
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to
# deal in the Software without restriction, including without limitation the
# rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
# sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

# The functions below are adapted from:
# https://github.com/JuliaDSP/DSP.jl/blob/master/src/unwrap.jl
# keeping the unwrapping along one dimension.

"""
    unwrap!(m; kwargs...)

In-place version of [`unwrap`](@ref).
"""
unwrap!(m::AbstractArray; kwargs...) = unwrap!(m, m; kwargs...)

"""
    unwrap!(y, m; dims = nothing, range = 2pi)

Unwrap `m` storing the result in `y`, see [`unwrap`](@ref). `dims` must be
given for an array of more than one dimension.
"""
function unwrap!(y::AbstractArray{T,N}, m::AbstractArray{T,N}; dims=nothing,
    range=2T(pi), kwargs...) where {T<:Real,N}
    if dims === nothing
        if N != 1
            throw(ArgumentError("`unwrap!`: required keyword parameter dims missing"))
        end
        dims = 1
    end
    if dims isa Integer
        accumulate!(unwrap_kernel(range), y, m; dims)
    else
        throw(ArgumentError("`unwrap!`: Invalid dims specified: $dims"))
    end
    return y
end

unwrap_kernel(range) = (x, y) -> y - round((y - x) / range) * range

"""
    unwrap(m; dims = nothing, range = 2pi)

Assumes `m` to be a sequence of real values, such as phases, that has been
wrapped to be inside the given `range` (centered around zero), and undoes
the wrapping by identifying discontinuities: each value is moved by a
multiple of `range` to within half of `range` of the value before it.
`dims` is the dimension along which to unwrap, required for an array of
more than one dimension; the array is unwrapped along that dimension only.
Complex values are not phases and are refused: unwrap `angle.(z)`.

A common usage is a phase measured over time or over frequency, such as
the phase of a scattering parameter, which `angle` wraps to stay within
(-pi, pi].

# Arguments
- `m::AbstractArray{T, N}`: Array of real values to unwrap.
- `dims=nothing`: Dimension along which to unwrap.
- `range=2pi`: Range of wrapped array.
"""
unwrap(m::AbstractArray; kwargs...) = unwrap!(similar(m), m; kwargs...)
