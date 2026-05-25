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

# The two functions below are from:
# https://github.com/JuliaDSP/DSP.jl/blob/master/src/unwrap.jl

"""
    unwrap!(m; kwargs...)

In-place version of [`unwrap`](@ref).
"""
unwrap!(m::AbstractArray; kwargs...) = unwrap!(m, m; kwargs...)

"""
    unwrap!(y, m; kwargs...)

Unwrap `m` storing the result in `y`, see [`unwrap`](@ref).
"""
function unwrap!(y::AbstractArray{T,N}, m::AbstractArray{T,N}; dims=nothing, range=2T(pi), kwargs...) where {T,N}
    if dims === nothing
        if N != 1
            throw(ArgumentError("`unwrap!`: required keyword parameter dims missing"))
        end
        dims = 1
    end
    if dims isa Integer
        accumulate!(unwrap_kernel(range), y, m; dims)
    # elseif dims == 1:N ## commented out. no need for unwrap_nd! at the moment
    #     unwrap_nd!(y, m; range, kwargs...)
    else
        throw(ArgumentError("`unwrap!`: Invalid dims specified: $dims"))
    end
    return y
end

unwrap_kernel(range) = (x, y) -> y - round((y - x) / range) * range

"""
    unwrap(m; kwargs...)

Assumes `m` to be a sequence of values that has been wrapped to be inside the
given `range` (centered around zero), and undoes the wrapping by identifying
discontinuities. If a single dimension is passed to `dims`, then `m` is assumed
to have wrapping discontinuities only along that dimension. If a range of
dimensions, as in `1:ndims(m)`, is passed to `dims`, then `m` is assumed to have
wrapping discontinuities across all `ndims(m)` dimensions.

A common usage for unwrapping across a singleton dimension is for a phase
measurement over time, such as when
comparing successive frames of a short-time Fourier transform, as
each frame is wrapped to stay within (-pi, pi].

A common usage for unwrapping across multiple dimensions is for a phase
measurement of a scene, such as when retrieving the phase information
of an image, as each pixel is wrapped to stay within (-pi, pi].

# Arguments
- `m::AbstractArray{T, N}`: Array to unwrap.
- `dims=nothing`: Dimensions along which to unwrap. If `dims` is an integer, then
    `unwrap` is called on that dimension. If `dims=1:ndims(m)`, then `m` is unwrapped
    across all dimensions.
- `range=2pi`: Range of wrapped array.
- `circular_dims=(false, ...)`:  When an element of this tuple is `true`, the
    unwrapping process will consider the edges along the corresponding axis
    of the array to be connected.
- `rng=default_rng()`: Unwrapping of arrays with dimension > 1 uses a random
    initialization. A user can pass their own RNG through this argument.
"""
unwrap(m::AbstractArray; kwargs...) = unwrap!(similar(m), m; kwargs...)

#= Algorithm based off of
 M. A. Herráez, D. R. Burton, M. J. Lalor, and M. A. Gdeisat,
 "Fast two-dimensional phase-unwrapping algorithm based on sorting by reliability following a noncontinuous path"
 `Applied Optics, Vol. 41, Issue 35, pp. 7437-7444 (2002) <http://dx.doi.org/10.1364/AO.41.007437>`
 and
 H. Abdul-Rahman, M. Gdeisat, D. Burton, M. Lalor,
 "Fast three-dimensional phase-unwrapping algorithm based on sorting by reliability following a non-continuous path",
 `Proc. SPIE 5856, Optical Measurement Systems for Industrial Inspection IV, 32 (2005) <http://dx.doi.ogr/doi:10.1117/12.611415>`
 Code inspired by Scipy's implementation, which is under BSD license.
=#
