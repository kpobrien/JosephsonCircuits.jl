"""
    JosephsonCircuitsForwardDiffExt

ForwardDiff.jl support, loaded with `using ForwardDiff`: the derivatives of
the values function of `designsensitivities` and `designjacobian` when no
`jacobian` is given.
"""
module JosephsonCircuitsForwardDiffExt

using ForwardDiff
import JosephsonCircuits

# The Jacobian of the values `f(p)`, the names `V`, with respect to the
# real design parameters `p`, a row per name and a column per parameter:
# that of their real and imaginary parts as real functions of the real
# parameters, so that a value a parameter rotates in the complex plane is
# exact.
function valuejacobian(f, p::NamedTuple{P}, V) where {P}
    function parts(x)
        y = f(NamedTuple{P}(Tuple(x)))
        keys(y) == V || throw(ArgumentError(lazy"the values function gives the names $(keys(y)) at one point and $(V) at another; it gives the same names in the same order at every point."))
        t = Tuple(y)
        return [real.(t)...; imag.(t)...]
    end
    R = ForwardDiff.jacobian(parts, Float64[x for x in p])
    n = length(V)
    return complex.(R[1:n, :], R[n+1:2n, :])
end

end # module
