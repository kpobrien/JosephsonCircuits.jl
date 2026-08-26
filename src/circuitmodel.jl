
"""
    AbstractComponent

The abstract supertype for component models in the typed circuit
representation. Component data belongs in concrete structs which subtype
`AbstractComponent`; behavior is provided by dispatched methods such as
[`nterminals`](@ref) and [`hasports`](@ref). A [`Circuit`](@ref) with an
[`Interface`](@ref) is also usable as a component.
"""
abstract type AbstractComponent end

"""
    ComponentNotSupportedError(msg)

An exception thrown when a component in the typed circuit representation
parses, validates, and elaborates successfully but is not yet supported by
the numerical solvers.
"""
struct ComponentNotSupportedError <: Exception
    msg::String
end
Base.showerror(io::IO, e::ComponentNotSupportedError) = print(io,
    "ComponentNotSupportedError: ", e.msg)

# === lumped linear components ===

"""
    Inductor(L)

A two terminal linear inductor with inductance `L` in Henries. Terminals are
`1` and `2` with orientation from terminal 1 to terminal 2. The value may be
a number or a symbolic variable.

# Examples
```jldoctest
julia> Inductor(1e-9)
Inductor{Float64}(1.0e-9)
```
"""
struct Inductor{T} <: AbstractComponent
    L::T
end

"""
    Capacitor(C)

A two terminal linear capacitor with capacitance `C` in Farads. Terminals are
`1` and `2`. The value may be a number or a symbolic variable.

# Examples
```jldoctest
julia> Capacitor(100e-15)
Capacitor{Float64}(1.0e-13)
```
"""
struct Capacitor{T} <: AbstractComponent
    C::T
end

"""
    Resistor(R)

A two terminal linear resistor with resistance `R` in Ohms. Terminals are `1`
and `2`. The value may be a number or a symbolic variable.

# Examples
```jldoctest
julia> Resistor(50.0)
Resistor{Float64}(50.0)
```
"""
struct Resistor{T} <: AbstractComponent
    R::T
end

# === sources and analysis ports ===

"""
    CurrentSource(I)

A two terminal current source with current `I` in Amperes flowing from
terminal 1 to terminal 2 through the source. The value is typically a
symbolic variable whose numerical value is supplied through `circuitdefs`
or through the analysis sources.
"""
struct CurrentSource{T} <: AbstractComponent
    I::T
end

"""
    VoltageSource(V)

A two terminal voltage source with voltage `V` in Volts between terminal 1
and terminal 2.
"""
struct VoltageSource{T} <: AbstractComponent
    V::T
end

"""
    Port(number::Integer; Z0 = 50.0)

An analysis port marker with the port number `number` and reference impedance
`Z0` in Ohms. A `Port` identifies an electrical port for excitation and
observation; excitation amplitudes belong to the analysis arguments, not to
the circuit topology. `Z0` is carried as metadata for wave-based analyses.

The legacy solver obtains the port termination resistance from a parallel
[`Resistor`](@ref), exactly as with the legacy `P` component.

# Examples
```jldoctest
julia> Port(1)
Port{Int64}(1, 50.0, 1)
```
"""
struct Port{T} <: AbstractComponent
    number::Int
    Z0::Float64
    value::T
end
Port(number::Integer; Z0 = 50.0) = Port{Int}(Int(number), Float64(Z0), Int(number))

# === mutual inductors ===

"""
    MutualInductor(K, inductor1, inductor2)

A mutual inductor with the dimensionless mutual coupling coefficient `K`
coupling the two [`Inductor`](@ref) instances with identifiers `inductor1`
and `inductor2` in the same circuit level. A mutual inductor couples two
named inductor branches rather than nets, so it appears in the component
list with no entries in the connections:

```julia
:k1 => MutualInductor(0.9, :l1, :l2)
```

The referenced identifiers are resolved within the containing circuit during
elaboration, so each instance of a subcircuit couples its own inductors.
"""
struct MutualInductor{T,I1,I2} <: AbstractComponent
    K::T
    inductor1::I1
    inductor2::I2
end

# === nonlinear inductive elements ===

"""
    PolynomialCPR(coefficients)

A current-phase relation (CPR) specified by the coefficients of its
polynomial expansion `f(φ) = coefficients[1]*φ + coefficients[2]*φ^2 + ...`,
where `φ` is the reduced branch phase. The linear coefficient must equal one
so that the `L0` of the containing [`NonlinearInductor`](@ref) is the small
signal inductance. The object is callable and its analytic derivative is
available through [`cprderivative`](@ref).

This supports specifying the effective nonlinearity of a SNAIL, SQUID,
Quarton, or kinetic inductor directly through its expansion coefficients
without wiring up the underlying junction arrangement.

# Examples
```jldoctest
julia> p = PolynomialCPR([1.0, 0.0, -1/6]); p(0.1)
0.09983333333333333

julia> JosephsonCircuits.cprderivative(p)(0.0)
1.0
```
"""
struct PolynomialCPR{T}
    a::Vector{T}   # f(φ) = a[1] + a[2]*φ + a[3]*φ^2 + ...
    PolynomialCPR{T}(a::Vector{T}) where T = new{T}(a)
end

function PolynomialCPR(coefficients::AbstractVector{T}) where T
    if isempty(coefficients)
        throw(ArgumentError("PolynomialCPR requires at least the linear coefficient."))
    end
    c1 = coefficients[1]
    if c1 isa Real && !isapprox(c1, one(c1); atol = 1e-12)
        throw(ArgumentError(lazy"The linear coefficient of a PolynomialCPR must equal one so that L0 is the small signal inductance; got $(c1). Rescale the coefficients and absorb the scale into L0."))
    end
    a = Vector{T}(undef, length(coefficients)+1)
    a[1] = zero(T)
    for i in eachindex(coefficients)
        a[i+1] = coefficients[i]
    end
    return PolynomialCPR{T}(a)
end

(p::PolynomialCPR)(φ) = evalpoly(φ, p.a)

"""
    PolynomialCPRDerivative(a)

The analytic derivative of a [`PolynomialCPR`](@ref), produced by
[`cprderivative`](@ref). Callable.
"""
struct PolynomialCPRDerivative{T}
    a::Vector{T}
end
(p::PolynomialCPRDerivative)(φ) = evalpoly(φ, p.a)

"""
    cprderivative(f)

Return the derivative of a current-phase relation. Built-in CPRs (`sin`,
[`PolynomialCPR`](@ref)) provide analytic derivatives. There is no automatic
differentiation or finite difference fallback: user supplied callables must
provide their derivative explicitly through the three argument
[`NonlinearInductor`](@ref) constructor.
"""
cprderivative(::typeof(sin)) = cos
function cprderivative(p::PolynomialCPR{T}) where T
    n = length(p.a)
    da = Vector{T}(undef, max(n-1,1))
    if n == 1
        da[1] = zero(T)
    else
        for k in 2:n
            da[k-1] = (k-1)*p.a[k]
        end
    end
    return PolynomialCPRDerivative{T}(da)
end
cprderivative(f) = throw(ArgumentError(lazy"No analytic derivative is known for the current-phase relation $(f). Supply it explicitly with NonlinearInductor(L0, f, df); automatic differentiation and finite differences are deliberately not used."))

"""
    NonlinearInductor(L0, cpr, dcpr)
    NonlinearInductor(L0, cpr)

A two terminal nonlinear inductive element defined by its current-phase
relation: `I(φ) = (phi0/L0)*cpr(φ)` where `φ` is the reduced branch phase,
`L0` is the small signal inductance in Henries, and `cpr` is a callable with
unit slope at zero. `dcpr` is the derivative of `cpr`, required by the
harmonic balance Jacobian; it must be supplied explicitly for user defined
callables, while built-in CPRs (`sin`, [`PolynomialCPR`](@ref)) provide
analytic derivatives through [`cprderivative`](@ref).

This supports specifying the effective nonlinearity of a SNAIL, SQUID, or
Quarton directly, as an alternative to composing the underlying junctions.
See also [`JosephsonJunction`](@ref).
"""
struct NonlinearInductor{T,F,DF} <: AbstractComponent
    L0::T
    cpr::F
    dcpr::DF
end
NonlinearInductor(L0, cpr) = NonlinearInductor(L0, cpr, cprderivative(cpr))

"""
    JosephsonJunction(Lj)
    JosephsonJunction(; Ic)

A Josephson junction with junction inductance `Lj` in Henries, or
equivalently critical current `Ic` in Amperes, and the sinusoidal
current-phase relation `I(φ) = Ic*sin(φ)`. Equal to
`NonlinearInductor(Lj, sin, cos)`.

# Examples
```jldoctest
julia> JosephsonJunction(100e-12) == NonlinearInductor(100e-12, sin, cos)
true

julia> JosephsonJunction(Ic = 3.29105976e-6).L0
1.0e-10
```
"""
JosephsonJunction(Lj) = NonlinearInductor(Lj, sin, cos)
JosephsonJunction(; Ic) = NonlinearInductor(phi0/Ic, sin, cos)

Base.:(==)(a::NonlinearInductor, b::NonlinearInductor) =
    isequal(a.L0, b.L0) && a.cpr == b.cpr && a.dcpr == b.dcpr

"""
    issinusoidal(c::NonlinearInductor)

Return true if the current-phase relation of `c` is the sinusoidal Josephson
relation, in which case the component lowers to the legacy `Lj` component.
"""
issinusoidal(c::NonlinearInductor) = c.cpr === sin && c.dcpr === cos

# === frequency dependent matrix providers ===

"""
    AbstractMatrixProvider

The abstract supertype for providers of frequency dependent matrix data such
as scattering parameters, noise covariance, and Gaussian channel `X` and `Y`
matrices. Providers implement [`evaluateprovider!`](@ref),
[`frequencydomain`](@ref), and [`providersize`](@ref). Large arrays and
interpolation data belong to the provider inside a shared component
definition and are never copied per instance.
"""
abstract type AbstractMatrixProvider end

"""
    ConstantMatrixProvider(A)

A frequency independent matrix provider wrapping the matrix `A`.
"""
struct ConstantMatrixProvider{M<:AbstractMatrix} <: AbstractMatrixProvider
    A::M
end

"""
    CallableMatrixProvider(f, n)

A matrix provider which evaluates the callable `f` at each requested angular
frequency, returning an `n` by `n` matrix.
"""
struct CallableMatrixProvider{F} <: AbstractMatrixProvider
    f::F
    n::Int
    # how `f` is called. `:matrix` returns a fresh n by n matrix, which is
    # the natural way to write a block by hand and is the default. `:inplace`
    # writes into one, `f(dest, w)`, which matters when the same provider is
    # evaluated many thousands of times and the returned matrices dominate
    # the allocation of a sweep. `:entry` returns one scalar,
    # `f(p, q, w) -> S[p,q]`, which is the only form a kernel can call, so it
    # is the one that lets a callable block be evaluated on a backend.
    form::Symbol
end

CallableMatrixProvider(f, n::Int) = CallableMatrixProvider(f, n, :matrix)

"""
    CALLABLE_FORMS

The ways a callable provider may be called. See
[`CallableMatrixProvider`](@ref).
"""
const CALLABLE_FORMS = (:matrix, :inplace, :entry)

"""
    TabulatedMatrixProvider(frequencies, values; interpolation = :linear,
        extrapolation = :error)

A matrix provider for tabulated data. `frequencies` is a strictly increasing
vector of angular frequencies in radians per second and `values` is an array
of dimensions (n, n, length(frequencies)). `interpolation` may be `:linear`.
`extrapolation` may be `:error` (the default), `:constant`, or `:linear`.
With `:error`, any requested frequency outside the tabulated range throws an
error listing the offending frequencies; extrapolation of tabulated data is
deliberately opt-in.
"""
struct TabulatedMatrixProvider{T} <: AbstractMatrixProvider
    frequencies::Vector{Float64}
    values::Array{T,3}
    interpolation::Symbol
    extrapolation::Symbol
end

function TabulatedMatrixProvider(frequencies::AbstractVector,
        values::AbstractArray{T,3};
        interpolation::Symbol = :linear,
        extrapolation::Symbol = :error) where T
    if size(values,1) != size(values,2)
        throw(DimensionMismatch(lazy"Tabulated matrix data must be square along the first two dimensions; got size $(size(values))."))
    end
    if size(values,3) != length(frequencies)
        throw(DimensionMismatch(lazy"The number of frequencies $(length(frequencies)) must match the third dimension of the data $(size(values,3))."))
    end
    if length(frequencies) < 1
        throw(ArgumentError("Tabulated matrix data requires at least one frequency."))
    end
    if !issorted(frequencies; lt = <)
        throw(ArgumentError("Tabulated frequencies must be strictly increasing."))
    end
    if interpolation != :linear
        throw(ArgumentError(lazy"Unknown interpolation $(interpolation). Supported: :linear."))
    end
    if !(extrapolation in (:error, :constant, :linear))
        throw(ArgumentError(lazy"Unknown extrapolation $(extrapolation). Supported: :error, :constant, :linear."))
    end
    return TabulatedMatrixProvider{T}(collect(Float64, frequencies),
        Array{T,3}(values), interpolation, extrapolation)
end

"""
    providersize(p::AbstractMatrixProvider)

Return the matrix dimension n of the n by n matrices produced by the
provider.
"""
providersize(p::ConstantMatrixProvider) = size(p.A,1)
providersize(p::CallableMatrixProvider) = p.n
providersize(p::TabulatedMatrixProvider) = size(p.values,1)

"""
    frequencydomain(p::AbstractMatrixProvider)

Return the tuple (fmin, fmax) of the valid angular frequency range of the
provider, or `nothing` if the provider is valid at any frequency.
"""
frequencydomain(p::ConstantMatrixProvider) = nothing
frequencydomain(p::CallableMatrixProvider) = nothing
frequencydomain(p::TabulatedMatrixProvider) = (p.frequencies[1], p.frequencies[end])

"""
    evaluateprovider!(dest::AbstractArray{T,3}, p::AbstractMatrixProvider,
        ws::AbstractVector)

Evaluate the provider at the angular frequencies `ws`, writing the matrix at
frequency `ws[i]` into `dest[:,:,i]`. This is the batch, in-place evaluation
contract used at analysis time: the caller evaluates each definition once
per frequency grid and caches the result, so instances sharing a definition
share one evaluation.
"""
function evaluateprovider!(dest::AbstractArray{T,3},
        p::ConstantMatrixProvider, ws::AbstractVector) where T
    checkdestsize(dest, providersize(p), length(ws))
    for i in eachindex(ws)
        dest[:,:,i] .= p.A
    end
    return dest
end

function evaluateprovider!(dest::AbstractArray{T,3},
        p::CallableMatrixProvider, ws::AbstractVector) where T
    checkdestsize(dest, providersize(p), length(ws))
    if p.form === :inplace
        for i in eachindex(ws)
            p.f(view(dest,:,:,i), ws[i])
        end
        return dest
    elseif p.form === :entry
        n = p.n
        for i in eachindex(ws)
            for q in 1:n, r in 1:n
                dest[r,q,i] = p.f(r, q, ws[i])
            end
        end
        return dest
    end
    for i in eachindex(ws)
        A = p.f(ws[i])
        if size(A) != (p.n, p.n)
            throw(DimensionMismatch(lazy"The callable provider returned a matrix of size $(size(A)) instead of ($(p.n), $(p.n)) at frequency $(ws[i])."))
        end
        dest[:,:,i] .= A
    end
    return dest
end

function evaluateprovider!(dest::AbstractArray{T,3},
        p::TabulatedMatrixProvider, ws::AbstractVector) where T
    checkdestsize(dest, providersize(p), length(ws))
    f = p.frequencies
    if p.extrapolation == :error
        outofrange = [w for w in ws if w < f[1] || w > f[end]]
        if !isempty(outofrange)
            throw(ArgumentError(lazy"The angular frequencies $(outofrange) are outside the tabulated range [$(f[1]), $(f[end])] rad/s. Extrapolation of tabulated data is opt-in: pass extrapolation = :constant or :linear if extrapolation is intended."))
        end
    end
    for i in eachindex(ws)
        w = ws[i]
        if w <= f[1]
            if p.extrapolation == :constant || length(f) == 1 || w == f[1]
                dest[:,:,i] .= view(p.values,:,:,1)
            else # linear extrapolation from the first segment
                lerpslices!(view(dest,:,:,i), p, 1, 2, w)
            end
        elseif w >= f[end]
            if p.extrapolation == :constant || length(f) == 1 || w == f[end]
                dest[:,:,i] .= view(p.values,:,:,length(f))
            else
                lerpslices!(view(dest,:,:,i), p, length(f)-1, length(f), w)
            end
        else
            j = searchsortedlast(f, w)
            if f[j] == w
                dest[:,:,i] .= view(p.values,:,:,j)
            else
                lerpslices!(view(dest,:,:,i), p, j, j+1, w)
            end
        end
    end
    return dest
end

function lerpslices!(dest, p::TabulatedMatrixProvider, j1::Int, j2::Int, w)
    f1 = p.frequencies[j1]
    f2 = p.frequencies[j2]
    t = (w - f1)/(f2 - f1)
    A = view(p.values,:,:,j1)
    B = view(p.values,:,:,j2)
    @. dest = (1-t)*A + t*B
    return dest
end

function checkdestsize(dest, n::Int, nf::Int)
    if size(dest) != (n, n, nf)
        throw(DimensionMismatch(lazy"The destination array has size $(size(dest)) but ($(n), $(n), $(nf)) is required."))
    end
    return nothing
end

"""
    matrixprovider(x, T; n = nothing, interpolation = :linear,
        extrapolation = :error)

Normalize user input into an [`AbstractMatrixProvider`](@ref) with element
type `T`. Accepts a matrix, a tuple `(frequencies, values)` of tabulated
data, an existing provider, or any callable of angular frequency.
"""
matrixprovider(p::AbstractMatrixProvider, ::Type{T}; kwargs...) where T = p
function matrixprovider(A::AbstractMatrix, ::Type{T}; n = nothing,
        form::Symbol = :matrix, kwargs...) where T
    checknoform(form, "a constant matrix")
    if size(A,1) != size(A,2)
        throw(DimensionMismatch(lazy"The matrix must be square; got size $(size(A))."))
    end
    if !isnothing(n) && size(A,1) != n
        throw(DimensionMismatch(lazy"The matrix has size $(size(A)) but dimension $(n) is required."))
    end
    return ConstantMatrixProvider(Matrix{T}(A))
end
function matrixprovider(t::Tuple{<:AbstractVector,<:AbstractArray{<:Any,3}},
        ::Type{T}; n = nothing, interpolation::Symbol = :linear,
        extrapolation::Symbol = :error, form::Symbol = :matrix) where T
    checknoform(form, "tabulated data")
    frequencies, values = t
    if !isnothing(n) && size(values,1) != n
        throw(DimensionMismatch(lazy"The tabulated data has matrix dimension $(size(values,1)) but dimension $(n) is required."))
    end
    return TabulatedMatrixProvider(frequencies, Array{T,3}(values);
        interpolation = interpolation, extrapolation = extrapolation)
end
"""
    checknoform(form::Symbol, what::AbstractString)

Reject a `form` for a provider which is not a callable. How a function is
called means nothing for data which is already stored.
"""
function checknoform(form::Symbol, what::AbstractString)
    if form !== :matrix
        throw(ArgumentError(lazy"form = $(repr(form)) was given with $(what), but it applies only to a callable provider: it says how that function is called."))
    end
    return nothing
end

function matrixprovider(f, ::Type{T}; n = nothing, form::Symbol = :matrix,
        kwargs...) where T
    if !(form in CALLABLE_FORMS)
        throw(ArgumentError(lazy"Unknown form $(repr(form)). Supported: :matrix (f(w) returns an n by n matrix), :inplace (f(dest, w) writes one), :entry (f(p, q, w) returns S[p,q])."))
    end
    if isnothing(n)
        throw(ArgumentError("The matrix dimension cannot be inferred from a callable provider; pass the dimension explicitly (nports for a ScatteringBlock, nmodes for a GaussianChannel)."))
    end
    return CallableMatrixProvider(f, n, form)
end

# === negative frequency rules ===

"""
    ConjugateSymmetry()

The default negative frequency rule for scattering and covariance providers:
data is evaluated at the absolute value of the requested frequency and
conjugated for negative frequencies, imposing S(-ω) = conj(S(ω)). This is
uniformly safe for tabulated positive-frequency data and for user callables.
"""
struct ConjugateSymmetry end

"""
    Native()

A negative frequency rule declaring that the provider is natively valid on
signed frequencies and may be evaluated directly at negative frequencies.
Opt-in for analytic providers whose formulas satisfy the physical
conjugation identity by construction.
"""
struct Native end

# === scattering block noise models ===

"""
    Passive()

The default noise model for a [`ScatteringBlock`](@ref): the block is a
passive deterministic network with no locally specified noise. Temperature
dependent thermal noise is supplied at analysis time, keeping temperature
out of shared component definitions.
"""
struct Passive end

"""
    ThermalEquilibrium(temperature)

A noise model for a passive [`ScatteringBlock`](@ref) in thermal equilibrium
at the common physical temperature `temperature` in Kelvin. The added noise
covariance is determined from the scattering data and the temperature by the
passive thermal channel relation Y(ω) = ν_T(ω) R(I - S(ω) S(ω)').
"""
struct ThermalEquilibrium{T}
    temperature::T
end

"""
    NoiseCovariance(V; interpolation = :linear, extrapolation = :error)

An arbitrary user specified noise covariance for a
[`ScatteringBlock`](@ref). `V` may be a matrix, a callable of angular
frequency, or a tuple `(frequencies, values)` of tabulated data, following
the same provider forms as the scattering data. Constant and tabulated
covariances are validated for Hermitian symmetry at construction.
"""
struct NoiseCovariance{P}
    provider::P
    interpolation::Symbol
    extrapolation::Symbol
end
function NoiseCovariance(V; interpolation::Symbol = :linear,
        extrapolation::Symbol = :error)
    return NoiseCovariance(V, interpolation, extrapolation)
end

# === scattering block ===

"""
    ScatteringBlock(S; nports = nothing, zref = 50.0, grounded = false,
        noise = Passive(), negative_frequency = ConjugateSymmetry(),
        interpolation = :linear, extrapolation = :error, atol = 1e-8)

A multiport component defined by its scattering parameters. `S` may be:

- a constant matrix;
- a callable of angular frequency returning a matrix (requires `nports`);
- a tuple `(frequencies, values)` of tabulated data with `frequencies` in
  radians per second and `values` of size (nports, nports, nfrequencies);
- a path to a Touchstone file, from which the reference impedance is also
  read.

`zref` is the reference impedance in Ohms, a scalar broadcast to all ports
or a vector with one entry per port. Scattering data is used at its native
reference impedance; no renormalization of the data is ever performed when
stamping, and conversion to analysis reference impedances happens only in
the wave domain at analysis boundaries.

Each port `p` has terminals `1` (signal) and `2` (reference), addressed as
`(:instance, p, t)`; ports may be floating or differential. With
`grounded = true` every reference terminal is automatically tied to
[`Ground`](@ref) and `(:instance, p)` in a connection group addresses the
signal terminal of port `p`; explicitly connecting a reference terminal of a
grounded block is an error.

`noise` is [`Passive`](@ref) (default), [`ThermalEquilibrium`](@ref), or
[`NoiseCovariance`](@ref). `negative_frequency` is
[`ConjugateSymmetry`](@ref) (default) or [`Native`](@ref). Passivity of
constant and tabulated scattering data is validated at construction with
absolute tolerance `atol` unless the noise model is a `NoiseCovariance`,
which permits active blocks.

# Examples
```jldoctest
julia> ScatteringBlock([0 1;1 0]).nports
2
```
"""
struct ScatteringBlock{P,N,NF} <: AbstractComponent
    provider::P
    nports::Int
    zref::Vector{Float64}
    grounded::Bool
    noise::N
    negative_frequency::NF
end

function ScatteringBlock(S; nports = nothing, zref = 50.0,
        grounded::Bool = false, noise = Passive(),
        negative_frequency = ConjugateSymmetry(),
        interpolation::Symbol = :linear, extrapolation::Symbol = :error,
        form::Symbol = :matrix,
        atol::Real = 1e-8)
    if S isa AbstractString
        return touchstonescatteringblock(S; nports = nports, zref = zref,
            grounded = grounded, noise = noise,
            negative_frequency = negative_frequency,
            interpolation = interpolation, extrapolation = extrapolation,
            atol = atol)
    end
    provider = matrixprovider(S, Complex{Float64}; n = nports,
        interpolation = interpolation, extrapolation = extrapolation,
        form = form)
    n = providersize(provider)
    if !isnothing(nports) && n != nports
        throw(DimensionMismatch(lazy"nports = $(nports) does not match the scattering data dimension $(n)."))
    end
    zrefvec = zrefvector(zref, n)
    if !(noise isa NoiseCovariance)
        checkpassive(provider; atol = atol)
    end
    checknoise(noise, n)
    return ScatteringBlock(provider, n, zrefvec, grounded, noise,
        negative_frequency)
end

function zrefvector(zref, n::Int)
    if zref isa AbstractVector
        if length(zref) != n
            throw(DimensionMismatch(lazy"zref has length $(length(zref)) but the block has $(n) ports."))
        end
        z = collect(Float64, zref)
    else
        z = fill(Float64(zref), n)
    end
    for zi in z
        if !(zi > 0) || !isfinite(zi)
            throw(ArgumentError(lazy"Reference impedances must be positive and finite; got $(z)."))
        end
    end
    return z
end

function checknoise(noise::NoiseCovariance, n::Int)
    provider = matrixprovider(noise.provider, Complex{Float64}; n = n,
        interpolation = noise.interpolation,
        extrapolation = noise.extrapolation)
    if providersize(provider) != n
        throw(DimensionMismatch(lazy"The noise covariance dimension $(providersize(provider)) does not match the number of ports $(n)."))
    end
    checkhermitian(provider)
    return nothing
end
checknoise(noise, n::Int) = nothing

function checkhermitian(p::ConstantMatrixProvider; atol = 1e-8)
    if !ishermitiantol(p.A, atol)
        throw(ArgumentError("The noise covariance matrix must be Hermitian."))
    end
    return nothing
end
function checkhermitian(p::TabulatedMatrixProvider; atol = 1e-8)
    for k in axes(p.values,3)
        if !ishermitiantol(view(p.values,:,:,k), atol)
            throw(ArgumentError(lazy"The noise covariance matrix at frequency index $(k) must be Hermitian."))
        end
    end
    return nothing
end
checkhermitian(p::AbstractMatrixProvider; atol = 1e-8) = nothing

function ishermitiantol(A, atol)
    for j in axes(A,2), i in axes(A,1)
        if abs(A[i,j] - conj(A[j,i])) > atol
            return false
        end
    end
    return true
end

"""
    passivitymargin(S::AbstractMatrix)

Return the minimum eigenvalue of I - S S', which is nonnegative for a
passive scattering matrix.
"""
function passivitymargin(S::AbstractMatrix)
    n = size(S,1)
    M = Matrix{Complex{Float64}}(I, n, n) - S*S'
    return minimum(real.(eigvals(Hermitian(M))))
end

function checkpassive(p::ConstantMatrixProvider; atol = 1e-8)
    margin = passivitymargin(p.A)
    if margin < -atol
        throw(ArgumentError(lazy"The scattering matrix is not passive: the minimum eigenvalue of I - S*S' is $(margin). For an active block, supply the noise with NoiseCovariance."))
    end
    return nothing
end
function checkpassive(p::TabulatedMatrixProvider; atol = 1e-8)
    worst = Inf
    worstindex = 0
    for k in axes(p.values,3)
        margin = passivitymargin(view(p.values,:,:,k))
        if margin < worst
            worst = margin
            worstindex = k
        end
    end
    if worst < -atol
        throw(ArgumentError(lazy"The tabulated scattering data is not passive: the minimum eigenvalue of I - S*S' is $(worst) at frequency index $(worstindex) ($(p.frequencies[worstindex]) rad/s). For an active block, supply the noise with NoiseCovariance."))
    end
    return nothing
end
checkpassive(p::AbstractMatrixProvider; atol = 1e-8) = nothing

"""
    evaluatescattering!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringBlock, ws::AbstractVector)

Evaluate the scattering parameters of `block` at the signed angular
frequencies `ws`, applying the negative frequency rule of the block, and
write the matrix at frequency `ws[i]` into `dest[:,:,i]`.
"""
function evaluatescattering!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringBlock, ws::AbstractVector,
        absbuffer::Union{Nothing,Vector{Float64}} = nothing)
    if block.negative_frequency isa Native
        evaluateprovider!(dest, block.provider, ws)
    else
        # `abs.(ws)` per call is one allocation per block per frequency on a
        # line whose every cell is a block; the caller may pass a buffer
        absws = if isnothing(absbuffer)
            abs.(ws)
        else
            resize!(absbuffer, length(ws))
            @inbounds for i in eachindex(ws)
                absbuffer[i] = abs(ws[i])
            end
            absbuffer
        end
        evaluateprovider!(dest, block.provider, absws)
        for i in eachindex(ws)
            if ws[i] < 0
                dv = view(dest,:,:,i)
                dv .= conj.(dv)
            end
        end
    end
    return dest
end

# Touchstone loading. The reference impedance is read from the file option
# line; a conflicting explicit zref is an error because the intent is
# ambiguous between correcting a mislabeled file and requesting
# renormalization.
function touchstonescatteringblock(path::AbstractString; nports, zref,
        grounded, noise, negative_frequency, interpolation, extrapolation,
        atol)
    ts = Touchstone.touchstone_load(path)
    filezref = collect(Float64, ts.reference)
    if !(zref isa Number && zref == 50.0) # user supplied a non-default zref
        z = zref isa AbstractVector ? collect(Float64, zref) :
            fill(Float64(zref), length(filezref))
        if length(z) != length(filezref) || any(!isapprox(z[i], filezref[i])
                for i in eachindex(filezref))
            throw(ArgumentError(lazy"The Touchstone file $(path) declares reference impedances $(filezref) Ohms but zref = $(zref) was supplied. The intent is ambiguous between correcting a mislabeled file and requesting renormalization, so this is an error; omit zref to use the file value."))
        end
    end
    frequencies = 2 .* pi .* collect(Float64, ts.f)
    values = Array{Complex{Float64},3}(ts.N)
    provider = TabulatedMatrixProvider(frequencies, values;
        interpolation = interpolation, extrapolation = extrapolation)
    n = providersize(provider)
    if !isnothing(nports) && n != nports
        throw(DimensionMismatch(lazy"nports = $(nports) does not match the Touchstone data dimension $(n)."))
    end
    if !(noise isa NoiseCovariance)
        checkpassive(provider; atol = atol)
    end
    checknoise(noise, n)
    return ScatteringBlock(provider, n, filezref, grounded, noise,
        negative_frequency)
end

"""
    TransmissionLineProvider(Z0, delay)

The scattering parameter provider of an ideal lossless transmission line of
characteristic impedance `Z0` and one way delay `delay` seconds, referenced
to `Z0`: S11 = S22 = 0 and S21 = S12 = exp(-im*ω*delay). Natively valid on
signed frequencies.
"""
struct TransmissionLineProvider <: AbstractMatrixProvider
    Z0::Float64
    delay::Float64
end
providersize(p::TransmissionLineProvider) = 2
frequencydomain(p::TransmissionLineProvider) = nothing
function evaluateprovider!(dest::AbstractArray{T,3},
        p::TransmissionLineProvider, ws::AbstractVector) where T
    checkdestsize(dest, 2, length(ws))
    for i in eachindex(ws)
        s21 = exp(-im*ws[i]*p.delay)
        dest[1,1,i] = 0
        dest[2,2,i] = 0
        dest[2,1,i] = s21
        dest[1,2,i] = s21
    end
    return dest
end

"""
    TransmissionLine(Z0, len; vp = speed_of_light, grounded = false,
        noise = Passive())

An ideal lossless transmission line of characteristic impedance `Z0` Ohms,
length `len` meters, and phase velocity `vp` meters per second, as a two
port [`ScatteringBlock`](@ref) referenced to `Z0` with an analytic provider
which is exact at any frequency ([`Native`](@ref) signed frequency
evaluation).

# Examples
```jldoctest
julia> TransmissionLine(50.0, 1e-3).nports
2
```
"""
function TransmissionLine(Z0, len; vp = speed_of_light,
        grounded::Bool = false, noise = Passive())
    if !(Z0 > 0) || !(len >= 0) || !(vp > 0)
        throw(ArgumentError("TransmissionLine requires Z0 > 0, len >= 0, and vp > 0."))
    end
    provider = TransmissionLineProvider(Float64(Z0), Float64(len)/Float64(vp))
    return ScatteringBlock(provider, 2, fill(Float64(Z0), 2), grounded,
        noise, Native())
end

# === Gaussian channels ===

"""
    symplecticform(n::Integer)

Return the 2n by 2n symplectic form Ω = [0 I; -I 0] in the real quadrature
ordering (x_1,…,x_n,p_1,…,p_n).
"""
function symplecticform(n::Integer)
    Ω = zeros(Float64, 2n, 2n)
    for i in 1:n
        Ω[i, n+i] = 1.0
        Ω[n+i, i] = -1.0
    end
    return Ω
end

"""
    completepositivitymargin(X::AbstractMatrix, Y::AbstractMatrix)

Return the minimum eigenvalue of Y + (i/2)(Ω - X Ω X'), which is nonnegative
for a completely positive Gaussian channel in the real quadrature
representation with vacuum covariance I/2.
"""
function completepositivitymargin(X::AbstractMatrix, Y::AbstractMatrix)
    n2 = size(X,1)
    if iszero(n2 % 2) == false
        throw(DimensionMismatch("Gaussian channel matrices must have even dimension 2n."))
    end
    Ω = symplecticform(n2 ÷ 2)
    M = Complex{Float64}.(Y) .+ (im/2).*(Ω .- X*Ω*transpose(X))
    return minimum(real.(eigvals(Hermitian(M))))
end

"""
    quadraturetransform(A::AbstractMatrix, B::AbstractMatrix)

Convert the complex Bogoliubov transformation b = A a + B conj(a) to the
real quadrature transformation X in the ordering (x_1,…,x_n,p_1,…,p_n), so
that d_out = X d_in. Returns the 2n by 2n real matrix
X = [Re(A+B) -Im(A-B); Im(A+B) Re(A-B)].

# Examples
```jldoctest
julia> quadraturetransform([0 1;1 0], zeros(2,2)) == [0 1 0 0;1 0 0 0;0 0 0 1;0 0 1 0]
true
```
"""
function quadraturetransform(A::AbstractMatrix, B::AbstractMatrix)
    if size(A) != size(B) || size(A,1) != size(A,2)
        throw(DimensionMismatch(lazy"A and B must be square with equal size; got $(size(A)) and $(size(B))."))
    end
    return [real.(A .+ B) -imag.(A .- B); imag.(A .+ B) real.(A .- B)]
end

"""
    GaussianChannel(X, Y; nmodes = nothing, displacement = nothing,
        grounded = false, interpolation = :linear, extrapolation = :error,
        atol = 1e-8)

An arbitrary Gaussian bosonic channel in the canonical real quadrature
representation: with quadratures ordered (x_1,…,x_n,p_1,…,p_n) and vacuum
covariance I/2, the channel acts as d_out = X d_in + d_0 and
V_out = X V_in X' + Y. `X` and `Y` are 2n by 2n real matrices and may each
be a constant matrix, a callable of angular frequency (requires `nmodes`),
or a tuple `(frequencies, values)` of tabulated data.

Complete positivity, Y + (i/2)(Ω - X Ω X') ⪰ 0, and the symmetry of Y are
validated pointwise at construction for constant and tabulated data with
absolute tolerance `atol`; the worst margin is recorded in the `cp_margin`
field (NaN when validation is deferred for callable providers).

Each mode is a two terminal port like a scattering block port, with the same
`grounded` convenience. The deterministic X participates in linearized
analyses; participation in nonlinear pump solves additionally requires a
coherent mean model, which is an analysis capability outside the input
format.

The complex Bogoliubov form b = A a + B conj(a) may be converted to the
deterministic part with [`quadraturetransform`](@ref).

# Examples
```jldoctest
julia> η = 0.5; abs(GaussianChannel(sqrt(η)*[1 0;0 1], (1-η)/2*[1 0;0 1]; nmodes=1).cp_margin) < 1e-10
true
```
"""
struct GaussianChannel{PX,PY,D} <: AbstractComponent
    X::PX
    Y::PY
    nmodes::Int
    displacement::D
    grounded::Bool
    cp_margin::Float64
end

function GaussianChannel(X, Y; nmodes = nothing, displacement = nothing,
        grounded::Bool = false, interpolation::Symbol = :linear,
        extrapolation::Symbol = :error, atol::Real = 1e-8)
    n2 = isnothing(nmodes) ? nothing : 2*nmodes
    Xp = matrixprovider(X, Float64; n = n2, interpolation = interpolation,
        extrapolation = extrapolation)
    Yp = matrixprovider(Y, Float64; n = isnothing(n2) ? providersize(Xp) : n2,
        interpolation = interpolation, extrapolation = extrapolation)
    if providersize(Xp) != providersize(Yp)
        throw(DimensionMismatch(lazy"X has dimension $(providersize(Xp)) but Y has dimension $(providersize(Yp))."))
    end
    if isodd(providersize(Xp))
        throw(DimensionMismatch(lazy"Gaussian channel matrices must have even dimension 2n in the real quadrature representation; got $(providersize(Xp))."))
    end
    n = providersize(Xp) ÷ 2
    if !isnothing(nmodes) && n != nmodes
        throw(DimensionMismatch(lazy"nmodes = $(nmodes) does not match the matrix dimension 2n = $(providersize(Xp))."))
    end
    margin = gaussianchannelmargin(Xp, Yp, atol)
    return GaussianChannel(Xp, Yp, n, displacement, grounded, margin)
end

function gaussianchannelmargin(Xp, Yp, atol)
    if Xp isa ConstantMatrixProvider && Yp isa ConstantMatrixProvider
        checkchannelpoint(Xp.A, Yp.A, atol, nothing)
        return completepositivitymargin(Xp.A, Yp.A)
    elseif Xp isa TabulatedMatrixProvider && Yp isa TabulatedMatrixProvider
        if Xp.frequencies != Yp.frequencies
            throw(ArgumentError("Tabulated X and Y must share the same frequency grid."))
        end
        worst = Inf
        for k in axes(Xp.values,3)
            checkchannelpoint(view(Xp.values,:,:,k), view(Yp.values,:,:,k),
                atol, k)
            worst = min(worst,
                completepositivitymargin(view(Xp.values,:,:,k),
                    view(Yp.values,:,:,k)))
        end
        return worst
    elseif Xp isa ConstantMatrixProvider && Yp isa TabulatedMatrixProvider
        worst = Inf
        for k in axes(Yp.values,3)
            checkchannelpoint(Xp.A, view(Yp.values,:,:,k), atol, k)
            worst = min(worst,
                completepositivitymargin(Xp.A, view(Yp.values,:,:,k)))
        end
        return worst
    elseif Xp isa TabulatedMatrixProvider && Yp isa ConstantMatrixProvider
        worst = Inf
        for k in axes(Xp.values,3)
            checkchannelpoint(view(Xp.values,:,:,k), Yp.A, atol, k)
            worst = min(worst,
                completepositivitymargin(view(Xp.values,:,:,k), Yp.A))
        end
        return worst
    else
        return NaN # deferred for callable providers
    end
end

function checkchannelpoint(X, Y, atol, k)
    where_ = isnothing(k) ? "" : " at frequency index $(k)"
    for j in axes(Y,2), i in axes(Y,1)
        if abs(Y[i,j] - Y[j,i]) > atol
            throw(ArgumentError(lazy"The added covariance Y must be symmetric$(where_)."))
        end
    end
    margin = completepositivitymargin(X, Y)
    if margin < -atol
        throw(ArgumentError(lazy"The Gaussian channel is not completely positive$(where_): the minimum eigenvalue of Y + (i/2)(Ω - X Ω X') is $(margin)."))
    end
    return nothing
end
