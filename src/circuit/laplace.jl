"""
    LaplaceResponse(f)

An explicitly supplied Laplace-domain response, with time convention `exp(s*t)`.
Calling the wrapper at a real angular frequency `w` returns `f(im*w)`, so it
can be used by harmonic balance as well as [`hbstability`](@ref):

```julia
Resistor(FrequencyDependent(LaplaceResponse(s -> 50*(1 + s/1e10))))
ScatteringParameters(LaplaceResponse(s -> fill(0.2/(1 + s/1e10), 1, 1));
    nports = 1, zref = 50.0)
```

For a lumped element, `f` returns its component value (ohms, farads or
henries), not its admittance. For a scattering provider it returns a matrix;
`form = :inplace` and `:entry` also work, with `f(dest, s)` and `f(p, q, s)`.
The response must be analytic wherever the pole search evaluates it and
must represent the intended causal model. Ordinary real components and
unconverted scattering responses require `f(conj(s)) == conj(f(s))`.
For a converted harmonic of `LinearizedScattering`, the opposite harmonic
is supplied by this reflection instead. No passivity, causality or global
analyticity proof is inferred from this wrapper. Internal states hidden
from a transfer function cannot be recovered; use `RationalScattering`
when their poles are required.

Frequency tables and arbitrary frequency-only callbacks do not specify a
Laplace continuation. Fit tables using `RationalScattering`, or explicitly
supply a model with this wrapper; `hbstability` never fits one silently.

# Examples
```jldoctest
julia> f = LaplaceResponse(s -> 50*(1 + s/1e10));

julia> f(1e10) == 50 + 50im
true

julia> r = Resistor(FrequencyDependent(f));
```
"""
struct LaplaceResponse{F}
    f::F
end

(response::LaplaceResponse)(w::Real) = response.f(im*w)
(response::LaplaceResponse)(dest::AbstractMatrix, w::Real) = response.f(dest, im*w)
(response::LaplaceResponse)(p::Integer, q::Integer, w::Real) = response.f(p, q, im*w)

# Only holomorphic operations are accepted in a component expression.
# In particular real(f(s)), imag(f(s)) and conj(f(s)) are not continuations
# of their frequency-axis expressions.
poleanalytic(x::Number) = isreal(x)
poleanalytic(x::CircuitValues.Constant) = isreal(x.val)
poleanalytic(x::CircuitValues.Provider) = x.f isa LaplaceResponse
poleanalytic(x::CircuitValues.Unary) = x.f in (-, inv, sqrt, exp, log) && poleanalytic(x.a)
poleanalytic(x::CircuitValues.Binary) = x.f in (+, -, *, /, ^) && poleanalytic(x.a) && poleanalytic(x.b)
poleanalytic(x) = false

laplacevalue(x::Number, s) = x
laplacevalue(x::CircuitValues.Constant, s) = x.val
laplacevalue(x::CircuitValues.Provider, s) = x.f.f(s)
laplacevalue(x::CircuitValues.Unary, s) = x.f(laplacevalue(x.a, s))
laplacevalue(x::CircuitValues.Binary, s) = x.f(laplacevalue(x.a, s), laplacevalue(x.b, s))
