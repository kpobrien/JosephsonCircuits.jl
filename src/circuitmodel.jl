
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
    Inductor(L; temperature = nothing)

A two terminal linear inductor with inductance `L` in Henries. Terminals are
`1` and `2` with orientation from terminal 1 to terminal 2. The value may be
a number or a symbolic variable.

`temperature` is the physical temperature in Kelvin, which sets the noise a
lossy instance adds; a lossless one adds none. `nothing`, the default, takes
the temperature the analysis is run at.

# Examples
```jldoctest
julia> Inductor(1e-9)
Inductor{Float64}(1.0e-9, nothing)

julia> Inductor(1e-9; temperature = 4.0).temperature
4.0
```
"""
struct Inductor{T} <: AbstractComponent
    L::T
    # the physical temperature in Kelvin, or `nothing` to take the one the
    # analysis is run at. Only a dissipative instance adds noise, so this is
    # read only where it does.
    temperature::Union{Nothing,Float64}
end
Inductor(L; temperature = nothing) = Inductor(L, temperature)

"""
    Capacitor(C; temperature = nothing)

A two terminal linear capacitor with capacitance `C` in Farads. Terminals are
`1` and `2`. The value may be a number or a symbolic variable.

`temperature` is the physical temperature in Kelvin, which sets the noise a
lossy instance adds; a lossless one adds none. `nothing`, the default, takes
the temperature the analysis is run at.

# Examples
```jldoctest
julia> Capacitor(100e-15)
Capacitor{Float64}(1.0e-13, nothing)

julia> Capacitor(100e-15; temperature = 4.0).temperature
4.0
```
"""
struct Capacitor{T} <: AbstractComponent
    C::T
    # the physical temperature in Kelvin, or `nothing` to take the one the
    # analysis is run at. Only a dissipative instance adds noise, so this is
    # read only where it does.
    temperature::Union{Nothing,Float64}
end
Capacitor(C; temperature = nothing) = Capacitor(C, temperature)

"""
    Resistor(R; temperature = nothing)

A two terminal linear resistor with resistance `R` in Ohms. Terminals are `1`
and `2`. The value may be a number or a symbolic variable.

`temperature` is the physical temperature in Kelvin, which sets the noise a
dissipative instance adds. `nothing`, the default, takes the temperature the
analysis is run at.

# Examples
```jldoctest
julia> Resistor(50.0)
Resistor{Float64}(50.0, nothing)

julia> Resistor(50.0; temperature = 4.0).temperature
4.0
```
"""
struct Resistor{T} <: AbstractComponent
    R::T
    # the physical temperature in Kelvin, or `nothing` to take the one the
    # analysis is run at. Only a dissipative instance adds noise, so this is
    # read only where it does.
    temperature::Union{Nothing,Float64}
end
Resistor(R; temperature = nothing) = Resistor(R, temperature)

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
    AbstractPortTermination

The supertype of the external environments a [`Port`](@ref) may own. A
termination is the source and load the port sees looking outward, which is
distinct from the reference impedance the port normalizes its waves to, even
where the two are numerically equal.
"""
abstract type AbstractPortTermination end

"""
    MatchedTermination()

A source and load environment matched to the port's reference impedance,
acting across the two port terminals. This is the default: a port owns its
environment, so no resistor should be added in order to terminate it.
"""
struct MatchedTermination <: AbstractPortTermination end

"""
    NoPortTermination()

No port owned environment, written as `termination = nothing`. The port
remains an excitation and observation boundary with its reference impedance
intact, but contributes no physical loading of its own.
"""
struct NoPortTermination <: AbstractPortTermination end

"""
    LegacyTermination(component)

The port owned environment of a legacy netlist: a resistor the netlist
already contains, named by its instance identifier.

Internal. A legacy netlist states a port's impedance by placing a resistor
across it and carries no role marker, so the adapter finds that resistor once
and records which one it is. Everything downstream then reads the port's
environment from the port, exactly as for a native matched port, and nothing
searches for a resistor sharing a port's branch.
"""
struct LegacyTermination{I} <: AbstractPortTermination
    component::I
end

porttermination(t::AbstractPortTermination) = t
porttermination(::Nothing) = NoPortTermination()
porttermination(x) = throw(ArgumentError(lazy"The port termination $(x) is not recognized. Write termination = nothing for an unterminated boundary, or omit the keyword for the default matched environment."))

# A port reference impedance must be a finite positive real number once it is
# numeric. Symbolic and deferred values pass through and are checked after
# binding, so the fallback method accepts anything.
checkportimpedance(Z0, number) = nothing
function checkportimpedance(
        Z0::Union{AbstractFloat,Integer,Rational,
            Complex{<:Union{AbstractFloat,Integer,Rational}}}, number)
    if !isreal(Z0) || !isfinite(Z0) || !(real(Z0) > 0)
        throw(ArgumentError(lazy"Port $(number) has reference impedance $(Z0), which must be finite, real, and positive."))
    end
    return nothing
end

"""
    Port(number::Integer; Z0 = 50.0, termination = MatchedTermination())

An analysis port with the port number `number` and reference impedance `Z0`
in Ohms. A `Port` identifies an electrical port for excitation and
observation; excitation amplitudes belong to the analysis arguments, not to
the circuit topology.

By default the port owns a matched external source and load environment of
impedance `Z0` acting across its two terminals, so a port needs no resistor
to define its impedance. The environment acts between the port terminals and
is never tied to [`Ground`](@ref) on its own, so a differential port behaves
the same way as a ground referenced one.

`termination = nothing` keeps `Z0` for wave normalization but adds no
physical loading, which is the right form when the circuit already contains
the resistor that terminates the port.

Any further [`Resistor`](@ref) across the same terminals is an ordinary
device resistor: it loads the port in parallel with the environment, it
remains a dissipative noise source, and it is never mistaken for the port's
own environment.

# Examples
```jldoctest
julia> Port(1)
Port(1; Z0 = 50.0)

julia> Port(2; Z0 = 1000.0, termination = nothing)
Port(2; Z0 = 1000.0, termination = nothing)
```
"""
struct Port{TZ,TT<:AbstractPortTermination} <: AbstractComponent
    number::Int
    Z0::TZ
    termination::TT
end

# One method and not several: keyword arguments do not participate in
# dispatch, so separate methods for the terminated and unterminated spellings
# would redefine one another rather than coexist, and the surviving one would
# silently answer for both.
function Port(number::Integer; Z0 = 50.0, termination = MatchedTermination())
    checkportimpedance(Z0, number)
    return Port(Int(number), Z0, porttermination(termination))
end

function Base.show(io::IO, p::Port)
    # the default is silent; anything else is stated, because a port which
    # owns no environment and one which owns an existing resistor are
    # different circuits
    print(io, "Port(", p.number, "; Z0 = ", p.Z0)
    p.termination isa NoPortTermination && print(io, ", termination = nothing")
    p.termination isa LegacyTermination &&
        print(io, ", termination = LegacyTermination(",
            repr(p.termination.component), ")")
    print(io, ")")
end

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
        throw(ArgumentError("The matrix dimension cannot be inferred from a callable provider; pass the dimension explicitly (nports for a ScatteringParameters, nmodes for a GaussianChannel)."))
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

The default noise model for a [`ScatteringParameters`](@ref): the block is a
passive network with no locally specified noise, so what it adds is set by
what it absorbs and by the temperature the analysis is run at. A dissipative
block adds noise of covariance `I - S S'` (see
[`ScatteringNoisePlan`](@ref)), scaled by the thermal factor of the
`temperature` or `temperatures` given to the analysis, which default to zero
temperature and so to the vacuum covariance itself.

Saying nothing here is what keeps temperature out of a component definition
which several analyses share. [`ThermalEquilibrium`](@ref) is how a block
states its own temperature instead.
"""
struct Passive end

"""
    Lossless()

A noise model for a [`ScatteringParameters`](@ref) whose scattering matrix is
unitary at every frequency, so that it absorbs nothing and adds no noise.

The default [`Passive`](@ref) gives a block noise channels unless its data
can be shown to be unitary, which is possible for a constant matrix and for a
table and not for a callable, whose values away from any sampled frequency
are unknown. A lossless callable therefore carries channels which are
identically zero: on a five hundred cell line that was measured at a third of
the run time on the host and half of it on a backend, spent computing zeros.
This is how to say so and get them back.

Where the data can be checked it is, at construction, so asserting this of a
constant or tabulated block which is not unitary is an error. For a callable
it cannot be checked and is taken on trust: asserting it of a block which
does absorb omits noise the block should add, and the quantum efficiency and
the commutation relations will be wrong by that much.
"""
struct Lossless end

# === the zero frequency model of a scattering block ===
#
# The direct current behavior of a block is its scattering matrix at zero
# frequency, and the default is to ask the block for it. That is right when
# the block can answer, and there are two ways it cannot. A measured or
# tabulated block may start at gigahertz and have no zero frequency entry at
# all. A closed form one may have a limit which exists but is not evaluable
# there: a series capacitor written `1/(im*w*C)` is an open circuit at
# direct current, and evaluating that expression at zero gives infinity
# rather than the open.
#
# So the model is stated separately from the radio frequency data when the
# two do not agree. Every realizable choice is some real `S(0)`, so one
# constructor carries them all and the common ones are named.

"""
    AbstractDCModel

The supertype of the zero frequency models a [`ScatteringParameters`](@ref)
may declare, written as `dcmodel`.

`ScatteringLimit()` (the default) evaluates the block's own scattering data
at zero frequency. [`OpenDC`](@ref), [`ShortDC`](@ref), [`ThroughDC`](@ref)
and [`ScatteringDC`](@ref) state it instead.
"""
abstract type AbstractDCModel end

"""
    ScatteringLimit()

The default zero frequency model of a [`ScatteringParameters`](@ref): its own
scattering data, evaluated at zero frequency.

The result must be real and finite there. A block whose limit exists but is
not evaluable at zero -- a series capacitance written `1/(im*w*C)`, whose
limit is the open circuit -- has to state that limit with one of the other
models rather than be asked for it.
"""
struct ScatteringLimit <: AbstractDCModel end

"""
    OpenDC()

A [`ScatteringParameters`](@ref) which is an open circuit at zero frequency:
`S(0) = I`, so no direct current flows through any port. This is the
constant of a series capacitance and of anything else which blocks direct
current.
"""
struct OpenDC <: AbstractDCModel end

"""
    ShortDC()

A [`ScatteringParameters`](@ref) each of whose ports is shorted to its own
reference terminal at zero frequency: `S(0) = -I`. The port voltages are
held at zero and the currents are whatever the rest of the circuit sends,
which is the ideal short the direct current rows are written to express.
"""
struct ShortDC <: AbstractDCModel end

"""
    ThroughDC()

A two port [`ScatteringParameters`](@ref) which passes direct current
unchanged: `S(0) = [0 1; 1 0]`, so the port voltages are equal and the
currents are equal and opposite. This is the constant of a series inductance
and of a transmission line.
"""
struct ThroughDC <: AbstractDCModel end

"""
    ScatteringDC(S0)

A [`ScatteringParameters`](@ref) whose zero frequency scattering matrix is
stated as `S0`, which must be real, finite, and of the block's dimension.

Every realizable zero frequency behavior is one of these: a resistor, a
transformer, an attenuator, and the open, short and through the named models
are shorthand for. It is validated for passivity at construction on the same
terms as the block's own data, so an active zero frequency model needs the
same `NoiseCovariance` declaration an active block does.
"""
struct ScatteringDC <: AbstractDCModel
    S0::Matrix{Float64}
end
function ScatteringDC(S0::AbstractMatrix)
    size(S0, 1) == size(S0, 2) ||
        throw(DimensionMismatch(lazy"a zero frequency scattering matrix must be square; got $(size(S0))."))
    all(isfinite, S0) ||
        throw(ArgumentError("a zero frequency scattering matrix must be finite."))
    all(iszero∘imag, S0) ||
        throw(ArgumentError("a zero frequency scattering matrix must be real: at zero frequency there is no phase to carry an imaginary part."))
    return ScatteringDC(Matrix{Float64}(real.(S0)))
end

"""
    dcscatteringmatrix(m::AbstractDCModel, nports)

The zero frequency scattering matrix a stated model describes.
"""
dcscatteringmatrix(::OpenDC, n::Integer) = Matrix{Float64}(I, n, n)
dcscatteringmatrix(::ShortDC, n::Integer) = Matrix{Float64}(-I, n, n)
function dcscatteringmatrix(::ThroughDC, n::Integer)
    n == 2 || throw(ArgumentError(lazy"ThroughDC is a two port model and this block has $(n) ports. Write the zero frequency matrix with ScatteringDC."))
    return [0.0 1.0; 1.0 0.0]
end
function dcscatteringmatrix(m::ScatteringDC, n::Integer)
    size(m.S0, 1) == n || throw(DimensionMismatch(lazy"the zero frequency scattering matrix has dimension $(size(m.S0,1)) but the block has $(n) ports."))
    return m.S0
end

"""
    ThermalEquilibrium(temperature)

A noise model for a passive [`ScatteringParameters`](@ref) in thermal equilibrium
at the physical temperature `temperature` in Kelvin. The added noise
covariance is the vacuum covariance `I - S S'` scaled by
`coth(hbar*w/(2*k*T))`, the factor by which a mode at that temperature
exceeds its vacuum noise.

This states the block's temperature where the block is defined, so it
overrides both the `temperature` and the `temperatures` arguments of the
analysis. `ThermalEquilibrium(0)` is [`Passive`](@ref).
"""
struct ThermalEquilibrium{T}
    temperature::T
end

"""
    NoiseCovariance(V; interpolation = :linear, extrapolation = :error)

An arbitrary user specified noise covariance for a
[`ScatteringParameters`](@ref). `V` may be a matrix, a callable of angular
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
    ScatteringParameters(S; nports = nothing, zref = 50.0, grounded = true,
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

With the default `grounded = true` every reference terminal is
automatically tied to [`Ground`](@ref) and `(:instance, p)` in a connection
group addresses the signal terminal of port `p`; explicitly connecting a
reference terminal of a grounded block is an error. With
`grounded = false` each port `p` has terminals `1` (signal) and `2`
(reference), addressed as `(:instance, p, t)`, and ports may be floating
or differential.

`noise` is [`Passive`](@ref) (default), [`Lossless`](@ref),
[`ThermalEquilibrium`](@ref), or [`NoiseCovariance`](@ref). A dissipative
block adds noise, which the noise scattering parameters, the quantum
efficiency and the commutation relations account for; of the four models they
support `Passive` and `Lossless` only, and a `NoiseCovariance` block does not
lower into a circuit at all. `Lossless` is how a callable which is unitary
says so, since unlike stored data it cannot be checked.
`negative_frequency` is [`ConjugateSymmetry`](@ref) (default) or
[`Native`](@ref). Passivity of constant and tabulated scattering data is
validated at construction with absolute tolerance `atol` unless the noise
model is a `NoiseCovariance`, which permits active blocks.

`dcmodel` states the block's zero frequency behavior when its own data does
not give it: [`OpenDC`](@ref), [`ShortDC`](@ref), [`ThroughDC`](@ref) or
[`ScatteringDC`](@ref), defaulting to [`ScatteringLimit`](@ref), which
evaluates the block at zero. Measured data which starts at gigahertz has no
zero frequency entry, and a closed form whose limit exists may not be
evaluable there -- a series capacitance written `1/(im*w*C)` is an open
circuit at direct current and infinite at zero -- so those state the limit
instead. The model is used only by the direct current rows; the alternating
current path always uses the block's own data.

`derivatives` supplies analytic derivatives of the scattering matrix with
respect to design parameters, for [`designsensitivities`](@ref): a named
tuple keyed by parameter name whose values are accepted in the same forms
as `S` (a matrix, a callable of angular frequency, or tabulated data). A
parameter the block depends on but has no entry for is differentiated by
central finite differences through `S` instead. A derivative is not a
scattering matrix and is never passivity checked.

# Examples
```jldoctest
julia> ScatteringParameters([0 1;1 0]).nports
2
```
"""
struct ScatteringParameters{P,N,NF,D,DM<:AbstractDCModel} <: AbstractComponent
    provider::P
    nports::Int
    zref::Vector{Float64}
    grounded::Bool
    noise::N
    negative_frequency::NF
    # analytic dS/dtheta providers keyed by design parameter name, for
    # [`designsensitivities`](@ref); empty when derivatives come from
    # finite differences through `provider`
    derivatives::D
    # the zero frequency behavior, when the block's own data does not give
    # it; see [`AbstractDCModel`](@ref)
    dcmodel::DM
end

# a block with no analytic derivatives and no stated zero frequency model,
# which is every block constructed before those keywords existed
ScatteringParameters(provider, nports::Int, zref::Vector{Float64},
    grounded::Bool, noise, negative_frequency) =
    ScatteringParameters(provider, nports, zref, grounded, noise,
        negative_frequency, NamedTuple(), ScatteringLimit())
ScatteringParameters(provider, nports::Int, zref::Vector{Float64},
    grounded::Bool, noise, negative_frequency, derivatives) =
    ScatteringParameters(provider, nports, zref, grounded, noise,
        negative_frequency, derivatives, ScatteringLimit())

function ScatteringParameters(S; nports = nothing, zref = 50.0,
        grounded::Bool = true, noise = Passive(),
        negative_frequency = ConjugateSymmetry(),
        interpolation::Symbol = :linear, extrapolation::Symbol = :error,
        form::Symbol = :matrix,
        derivatives::NamedTuple = NamedTuple(),
        dcmodel::AbstractDCModel = ScatteringLimit(),
        atol::Real = 1e-8)
    if S isa AbstractString
        return touchstonescatteringblock(S; nports = nports, zref = zref,
            grounded = grounded, noise = noise, dcmodel = dcmodel,
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
    checklossless(noise, provider)
    checknoise(noise, n)
    dprov = NamedTuple(k => begin
            dp = matrixprovider(v, Complex{Float64}; n = n, form = form)
            providersize(dp) == n || throw(DimensionMismatch(lazy"the derivative for parameter $(k) has dimension $(providersize(dp)) but the block has $(n) ports."))
            dp
        end for (k, v) in pairs(derivatives))
    checkdcmodel(dcmodel, n, noise, atol)
    return ScatteringParameters(provider, n, zrefvec, grounded, noise,
        negative_frequency, dprov, dcmodel)
end

# A stated zero frequency matrix is data like the block's own, so it is
# checked on the same terms and at the same time: the size against the
# block, and passivity unless the block declared itself active.
checkdcmodel(::ScatteringLimit, n::Int, noise, atol) = nothing
function checkdcmodel(m::AbstractDCModel, n::Int, noise, atol)
    S0 = dcscatteringmatrix(m, n)
    if !(noise isa NoiseCovariance)
        sv = maximum(svdvals(S0); init = 0.0)
        sv <= 1 + atol || throw(ArgumentError(lazy"the zero frequency scattering matrix has largest singular value $(sv), so it is active at direct current. An active block has to declare its own noise with NoiseCovariance, as it does at every other frequency."))
    end
    return nothing
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
    unitaritydeviation(S::AbstractMatrix)

The largest absolute entry of `I - S S'`, which is zero for a lossless
(unitary) scattering matrix. Unlike [`passivitymargin`](@ref) this sees
gain as well as loss, so it is the quantity a lossless test compares
against a tolerance.
"""
function unitaritydeviation(S::AbstractMatrix)
    n = size(S,1)
    worst = 0.0
    for q in 1:n
        for p in 1:n
            acc = p == q ? one(Complex{Float64}) : zero(Complex{Float64})
            for l in 1:n
                acc -= S[p,l]*conj(S[q,l])
            end
            worst = max(worst, abs(acc))
        end
    end
    return worst
end

"""
    provablylossless(provider::AbstractMatrixProvider; atol = 1e-10)
    provablylossless(block::ScatteringParameters; atol = 1e-10)

Whether the scattering data can be shown to be unitary at every frequency
from the data alone, which is possible for a constant matrix and for a
table and is not for a callable, whose values away from any sampled
frequency are unknown.

A block which is not provably lossless carries vacuum noise channels
([`ScatteringNoisePlan`](@ref)), and those of a block which is in fact
lossless are identically zero. A `false` here therefore costs work and
never correctness, which is why the fallback is `false`.
"""
provablylossless(p::AbstractMatrixProvider; atol = 1e-10) = false
provablylossless(p::ConstantMatrixProvider; atol = 1e-10) =
    unitaritydeviation(p.A) <= atol
function provablylossless(p::TabulatedMatrixProvider; atol = 1e-10)
    for k in axes(p.values,3)
        unitaritydeviation(view(p.values,:,:,k)) <= atol || return false
    end
    return true
end
provablylossless(b::ScatteringParameters; atol = 1e-10) =
    provablylossless(b.provider; atol = atol)

# `Lossless` is an assertion about every frequency, which a constant matrix
# and a table can be held to and a callable cannot
function checklossless(noise::Lossless, provider)
    if provider isa CallableMatrixProvider
        return nothing
    end
    if !provablylossless(provider)
        throw(ArgumentError("noise = Lossless() says the scattering matrix is unitary at every frequency, but this block's data is not: the largest absolute entry of I - S*S' over it is $(worstunitaritydeviation(provider)). Use the default Passive() noise model, which gives a dissipative block the noise its loss requires."))
    end
    return nothing
end
checklossless(noise, provider) = nothing

# the worst deviation over stored data, for the message above
worstunitaritydeviation(p::ConstantMatrixProvider) = unitaritydeviation(p.A)
function worstunitaritydeviation(p::TabulatedMatrixProvider)
    worst = 0.0
    for k in axes(p.values,3)
        worst = max(worst, unitaritydeviation(view(p.values,:,:,k)))
    end
    return worst
end

"""
    evaluatescattering!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringParameters, ws::AbstractVector)

Evaluate the scattering parameters of `block` at the signed angular
frequencies `ws`, applying the negative frequency rule of the block, and
write the matrix at frequency `ws[i]` into `dest[:,:,i]`.
"""
function evaluatescattering!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringParameters, ws::AbstractVector,
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
        atol, dcmodel::AbstractDCModel = ScatteringLimit())
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
    checkdcmodel(dcmodel, n, noise, atol)
    return ScatteringParameters(provider, n, filezref, grounded, noise,
        negative_frequency, NamedTuple(), dcmodel)
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
    TransmissionLine(Z0, len; vp = speed_of_light, grounded = true,
        noise = Passive())

An ideal lossless transmission line of characteristic impedance `Z0` Ohms,
length `len` meters, and phase velocity `vp` meters per second, as a two
port [`ScatteringParameters`](@ref) referenced to `Z0` with an analytic provider
which is exact at any frequency ([`Native`](@ref) signed frequency
evaluation).

# Examples
```jldoctest
julia> TransmissionLine(50.0, 1e-3).nports
2
```
"""
function TransmissionLine(Z0, len; vp = speed_of_light,
        grounded::Bool = true, noise = Passive())
    if !(Z0 > 0) || !(len >= 0) || !(vp > 0)
        throw(ArgumentError("TransmissionLine requires Z0 > 0, len >= 0, and vp > 0."))
    end
    provider = TransmissionLineProvider(Float64(Z0), Float64(len)/Float64(vp))
    return ScatteringParameters(provider, 2, fill(Float64(Z0), 2), grounded,
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
        grounded = true, interpolation = :linear, extrapolation = :error,
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
        grounded::Bool = true, interpolation::Symbol = :linear,
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
