# The component models of the typed circuit representation, and the
# matrix providers a scattering or Gaussian channel block reads its
# frequency dependent data from.

"""
    AbstractComponent

The supertype of the component models of the typed circuit
representation. A component is a concrete struct holding its data;
behavior is given by methods such as [`nterminals`](@ref) and
[`hasports`](@ref) in circuit/parse.jl. A [`Circuit`](@ref) with an
[`Interface`](@ref) is also a component.
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
    # the physical temperature in Kelvin, or `nothing` for the temperature
    # the analysis is run at; read only when the instance is dissipative
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
    # as for `Inductor`
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
    # as for `Inductor`
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

# normalize the `termination` keyword of `Port`: `nothing` means no
# termination
porttermination(t::AbstractPortTermination) = t
porttermination(::Nothing) = NoPortTermination()
porttermination(x) = throw(ArgumentError(lazy"The port termination $(x) is not recognized. Write termination = nothing for an unterminated boundary, or omit the keyword for the default matched environment."))

# A numeric port reference impedance must be finite, real and positive.
# A symbolic one passes through here and is checked after binding.
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
    Port(number::Integer; Z0 = 50.0,
        termination = JosephsonCircuits.MatchedTermination())

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
the resistor that terminates the port. Such a port is an ideal current
source and an impedance probe: a source on it drives its whole current into
the circuit, and since the incident wave is that of the source current with
nothing to absorb it, the reflection reported at the port is
`S = 2Z/Z0 - 1` with `Z` the impedance across its terminals, so the
impedance at the node is `Z = Z0 (1 + S)/2`.

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

# One keyword method covers every spelling; keyword arguments do not
# participate in dispatch, so separate methods would overwrite each other.
function Port(number::Integer; Z0 = 50.0, termination = MatchedTermination())
    checkportimpedance(Z0, number)
    return Port(Int(number), Z0, porttermination(termination))
end

# The default termination is not printed; any other is, since a port with
# no termination and one adopting an existing resistor are different
# circuits.
function Base.show(io::IO, p::Port)
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

In the netlist form of [`Circuit`](@ref) the entry names the two inductors
in place of nodes, and the component is written `MutualInductor(K)` alone:

```julia
(:k1, :l1, :l2, MutualInductor(0.9))
```
"""
struct MutualInductor{T,I1,I2} <: AbstractComponent
    K::T
    inductor1::I1
    inductor2::I2
end

# the inductors are named by the netlist entry
MutualInductor(K) = MutualInductor(K, nothing, nothing)

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
without wiring up the underlying junction arrangement. An array of `N`
identical junctions in series, for example, divides the phase and so has
the relation `N*sin(φ/N)` with small signal inductance `N*Lj`, whose
expansion is `[1, 0, -1/(6N^2), 0, 1/(120N^4), ...]`.

The coefficients are taken as given. A polynomial is not periodic and not
bounded, so where it leaves the range it was fitted on, and what the
harmonic count must be for the harmonics its degree generates, are the
user's to judge: the solver evaluates what it is given.

# Examples
```jldoctest
julia> p = PolynomialCPR([1.0, 0.0, -1/6]); p(0.1)
0.09983333333333333

julia> JosephsonCircuits.cprderivative(p)(0.0)
1.0
```
"""
struct PolynomialCPR{T}
    # the coefficients in `evalpoly` order, `f(φ) = a[1] + a[2]*φ + ...`,
    # with `a[1] = 0` so that the user's `coefficients[k]` is `a[k+1]`
    a::Vector{T}
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
    cprderivative(cpr)

The derivative of a current-phase relation as a callable. `sin` gives
`cos` and a [`PolynomialCPR`](@ref) gives its analytic derivative. Any
other callable throws an `ArgumentError`: there is no automatic
differentiation or finite difference fallback, so a user defined relation
must supply its derivative through the three argument
[`NonlinearInductor`](@ref) constructor.
"""
cprderivative(::typeof(sin)) = cos
function cprderivative(p::PolynomialCPR{T}) where T
    return PolynomialCPRDerivative{T}(differentiatecoefficients(p.a))
end
# The derivative of a derivative is the second derivative, which the
# Hessian of the harmonic balance system and the derivative of the
# linearized system with respect to the operating point both need. Taking
# it twice more gives zero polynomials rather than an error, so a caller
# does not have to know the degree.
cprderivative(p::PolynomialCPRDerivative{T}) where T =
    PolynomialCPRDerivative{T}(differentiatecoefficients(p.a))

# `evalpoly` order in, `evalpoly` order out. A polynomial of degree zero
# differentiates to the zero polynomial, written with the one coefficient
# `evalpoly` needs rather than as an empty vector.
function differentiatecoefficients(a::Vector{T}) where T
    n = length(a)
    n <= 1 && return T[zero(T)]
    da = Vector{T}(undef, n-1)
    for k in 2:n
        da[k-1] = (k-1)*a[k]
    end
    return da
end
cprderivative(f) = throw(ArgumentError(lazy"No analytic derivative is known for the current-phase relation $(f). Supply it explicitly with NonlinearInductor(L0, f, df); automatic differentiation and finite differences are deliberately not used."))

"""
    JunctionRelations(value, derivative, negsecond, sinusoidal, anysinusoidal)

The current-phase relations of the Josephson junction branches of a
circuit, in the order of the nonzero entries of the branch inductance
vector `Ljb`, which is the order of the junction axis of every time domain
array the solvers hold.

A relation is either the sinusoidal Josephson one or a
[`PolynomialCPR`](@ref). The polynomial coefficients of every junction sit
in the rows of `value`, in `evalpoly` order along the second axis and
padded with zeros to one common degree, so that a single Horner loop
evaluates them all; `derivative` and `negsecond` hold the coefficients of
the first derivative and of the *negative* of the second, which is the
combination the Hessian and the derivative of the linearized system with
respect to the operating point are written in, and `third` those of the
third derivative, which the trilinear form of the problem interface
takes. `sinusoidal` is true for
the junctions whose relation no polynomial represents, whose columns the
evaluation writes over with the trigonometric one; `anysinusoidal` is
whether any is, so that a circuit of polynomials alone skips that pass.

A circuit whose junctions are all sinusoidal, which is every circuit that
does not ask for anything else, has no table at all: the solvers hold
`nothing` and take the plain `sin` and `cos` they always did.
"""
struct JunctionRelations{M,V}
    value::M
    derivative::M
    negsecond::M
    third::M
    sinusoidal::V
    anysinusoidal::Bool
end

"""
    junctionrelations(cprs::AbstractVector)

The [`JunctionRelations`](@ref) of the junctions whose relations are
`cprs`, one entry per junction in the order of the junction axis, each
either `nothing` for the sinusoidal Josephson relation or a
[`PolynomialCPR`](@ref). Returns `nothing` when every entry is `nothing`.
"""
function junctionrelations(cprs::AbstractVector)
    any(!isnothing, cprs) || return nothing
    nj = length(cprs)
    # one common degree, so that the Horner loop is the same length for
    # every junction; a sinusoidal one contributes no terms
    nterms = maximum(c -> isnothing(c) ? 0 : length(c.a), cprs)
    value = zeros(Float64, nj, nterms)
    derivative = zeros(Float64, nj, nterms)
    negsecond = zeros(Float64, nj, nterms)
    third = zeros(Float64, nj, nterms)
    sinusoidal = fill(true, nj)
    for (j, c) in enumerate(cprs)
        isnothing(c) && continue
        sinusoidal[j] = false
        d = cprderivative(c)
        d2 = cprderivative(d)
        d3 = cprderivative(d2)
        copycoefficients!(value, j, c.a)
        copycoefficients!(derivative, j, d.a)
        copycoefficients!(negsecond, j, -1 .* d2.a)
        copycoefficients!(third, j, d3.a)
    end
    return JunctionRelations(value, derivative, negsecond, third, sinusoidal,
        any(sinusoidal))
end

"""
    relationat(r::JunctionRelations, phi)

The current-phase relation of every junction at the branch phases `phi`,
whose first axis is the junction: `sin.(phi)` when every relation is the
Josephson one, and the polynomials by Horner otherwise. Allocating, for
the setup, the diagnostics and the noise, which are not inside a step
loop; [`relationinto!`](@ref) is the in place form the steps take.

`phi` and the table must live on the same backend, so a host path takes a
host table from [`hostrelations`](@ref).
"""
relationat(r::JunctionRelations, phi) = relationinto!(similar(phi), r, phi)

"""
    derivativeat(r::JunctionRelations, phi)

The derivative of the relation of every junction at the branch phases
`phi`, `cos.(phi)` for the Josephson relation: the differential
inductance which every Jacobian and every linearization is built from.
The counterpart of [`relationat`](@ref), with
[`derivativeinto!`](@ref) as its in place form.
"""
derivativeat(r::JunctionRelations, phi) = derivativeinto!(similar(phi), r, phi)

"""
    relationinto!(out, r::JunctionRelations, phi)

[`relationat`](@ref) writing into `out`, which may not alias `phi`. A
circuit whose junctions are all sinusoidal takes the single broadcast of
`sin` it always did.
"""
function relationinto!(out, r::JunctionRelations, phi)
    allsinusoidal(r) && return out .= sin.(phi)
    return applyrelationfirst!(out, phi, r.value, r.sinusoidal,
        r.anysinusoidal, sin)
end

"""
    derivativeinto!(out, r::JunctionRelations, phi)

[`derivativeat`](@ref) writing into `out`, which may not alias `phi`. A
circuit whose junctions are all sinusoidal takes the single broadcast of
`cos` it always did.
"""
function derivativeinto!(out, r::JunctionRelations, phi)
    allsinusoidal(r) && return out .= cos.(phi)
    return applyrelationfirst!(out, phi, r.derivative, r.sinusoidal,
        r.anysinusoidal, cos)
end

"""
    hostrelations(r::JunctionRelations)
    hostrelations(r::JunctionRelations, rows::AbstractVector{Int})

The table on the host, for the host loops which cannot read a device
array, and with `rows` the table of that subset of the junctions, in that
order, which is what a projection onto part of the circuit reads.
"""
hostrelations(r::JunctionRelations) = JunctionRelations(Array(r.value),
    Array(r.derivative), Array(r.negsecond), Array(r.third),
    Array(r.sinusoidal), r.anysinusoidal)

function hostrelations(r::JunctionRelations, rows::AbstractVector{Int})
    allsinusoidal(r) && return hostrelations(r)
    sub(m) = Array(m)[rows, :]
    mask = Array(r.sinusoidal)[rows]
    return JunctionRelations(sub(r.value), sub(r.derivative),
        sub(r.negsecond), sub(r.third), mask, any(mask))
end

"""
    allsinusoidal(r::JunctionRelations)

Whether every junction of `r` is the sinusoidal Josephson one, which is
the empty table: the solvers hold a table always, so that their type does
not depend on what the circuit holds, and an empty one means the plain
`sin` and `cos` of the Josephson relation.
"""
allsinusoidal(r::JunctionRelations) = size(r.value, 2) == 0

"""
    emptyrelations(A::AbstractArray)

The empty [`JunctionRelations`](@ref), of the array types of `A`, which is
what a circuit of sinusoidal junctions alone holds.
"""
emptyrelations(A::AbstractArray) = JunctionRelations(similar(A, 0, 0),
    similar(A, 0, 0), similar(A, 0, 0), similar(A, 0, 0),
    similar(A, Bool, 0), true)

"""
    torelations(r, A::AbstractArray)

The relations `r`, `nothing` or a host [`JunctionRelations`](@ref), in the
array types of `A`: the coefficients take the working precision of `A` and
move to the backend it lives on, so that the solvers' broadcasts stay on
one device and in one precision.
"""
torelations(::Nothing, A::AbstractArray) = emptyrelations(A)
function torelations(r::JunctionRelations, A::AbstractArray)
    T = real(eltype(A))
    move(m) = copyto!(similar(A, T, size(m)...), convert(Matrix{T}, m))
    return JunctionRelations(move(r.value), move(r.derivative),
        move(r.negsecond), move(r.third),
        copyto!(similar(A, Bool, length(r.sinusoidal)), r.sinusoidal),
        r.anysinusoidal)
end

function copycoefficients!(m::Matrix, j::Integer, a::AbstractVector)
    for k in eachindex(a)
        m[j, k] = a[k]
    end
    return m
end

# `out .= f.(src)` where `f` is the relation of each junction, given by the
# rows of `c` as a polynomial in `evalpoly` order, and `trig` where the
# junction is sinusoidal. The junction is the last axis of `src`, as it is
# in every time domain array of harmonic balance, so a coefficient column
# is reshaped to broadcast along it. Everything here is a broadcast of
# whole arrays, which is what keeps it device generic: no kernel, no
# scalar indexing, and one pass per degree rather than one per junction.
function applyrelationlast!(out::AbstractArray, src::AbstractArray,
        c::AbstractMatrix, sinusoidal::AbstractVector, anysinusoidal::Bool, trig)
    shape = ntuple(d -> d == ndims(out) ? size(c, 1) : 1, ndims(out))
    column(k) = reshape(view(c, :, k), shape)
    nterms = size(c, 2)
    out .= column(nterms)
    for k in nterms-1:-1:1
        out .= out .* src .+ column(k)
    end
    if anysinusoidal
        mask = reshape(sinusoidal, shape)
        out .= ifelse.(mask, trig.(src), out)
    end
    return out
end

# the same for the transient, whose branch phases carry the junction on
# the first axis, `(junction,)` or `(junction, condition)`, against which a
# coefficient column broadcasts as it is
function applyrelationfirst!(out::AbstractArray, src::AbstractArray,
        c::AbstractMatrix, sinusoidal::AbstractVector, anysinusoidal::Bool, trig)
    nterms = size(c, 2)
    out .= view(c, :, nterms)
    for k in nterms-1:-1:1
        out .= out .* src .+ view(c, :, k)
    end
    if anysinusoidal
        out .= ifelse.(sinusoidal, trig.(src), out)
    end
    return out
end

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
A quadratic term in the relation is what makes such an element a three
wave mixer, which the Josephson relation, being odd, is not.

The relations the solvers evaluate are the sinusoidal Josephson one,
which is what [`JosephsonJunction`](@ref) writes, and a
[`PolynomialCPR`](@ref); another callable is refused when the circuit is
compiled. The element is a junction to everything else in the solvers: it
makes the same branch, enters the same matrices, and `L0` is its
inductance there, so only the pointwise relation differs. Harmonic
balance evaluates it, in the residual, in the Jacobian, in the Hessian
and in the pump modulation of the linearized system, and the transient
solver steps it, in its residual, its Jacobian, its tangent and adjoint,
and the linearization its noise is taken about.

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

Whether the current-phase relation of `c` is the sinusoidal Josephson
relation `(sin, cos)`, in which case the component compiles to the `:Lj`
type the solvers support.
"""
issinusoidal(c::NonlinearInductor) = c.cpr === sin && c.dcpr === cos

"""
    junctioncpr(c::NonlinearInductor, path)

The relation the solvers evaluate for `c`: `nothing` for the sinusoidal
Josephson relation, which they hold as `sin` and `cos`, or the
[`PolynomialCPR`](@ref) they evaluate by Horner. Any other callable throws,
since a relation the solvers cannot write down is one they cannot
transform.

The coefficients are converted to `Float64` here, so that the table the
solvers build is concrete whatever the user wrote them as.
"""
junctioncpr(c::NonlinearInductor, path = "") =
    issinusoidal(c) ? nothing : junctioncpr(c.cpr, path)
junctioncpr(p::PolynomialCPR, path) =
    PolynomialCPR{Float64}(convert(Vector{Float64}, p.a))
junctioncpr(f, path) = throw(ComponentNotSupportedError(lazy"the NonlinearInductor at $(path) has the current-phase relation $(f), which the solvers do not evaluate; give the sinusoidal Josephson relation or a PolynomialCPR."))

# === frequency dependent matrix providers ===

"""
    AbstractMatrixProvider

The supertype of the sources of frequency dependent matrix data: the
scattering parameters and noise covariance of a
[`ScatteringParameters`](@ref) and the `X` and `Y` matrices of a
[`GaussianChannel`](@ref). A provider implements
[`evaluateprovider!`](@ref) and [`providersize`](@ref). Its data lives in
the shared component definition and is never copied per instance.
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
    CallableMatrixProvider(f, n, form = :matrix)

A matrix provider which evaluates the callable `f` at each requested angular
frequency for an `n` by `n` matrix. `form` says how `f` is called:
`:matrix` (`f(w)` returns a fresh `n` by `n` matrix, the natural way to
write a block by hand and the default), `:inplace` (`f(dest, w)` writes
into one, which matters when the same provider is evaluated many thousands
of times), or `:entry` (`f(p, q, w)` returns `S[p, q]`, the only form a
kernel can call and so the one which lets a callable block be evaluated on
a backend); see [`CALLABLE_FORMS`](@ref).
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
    TabulatedMatrixProvider(frequencies, values; interpolation = :cubic,
        extrapolation = :error)

A matrix provider for tabulated data. `frequencies` is a strictly increasing
vector of angular frequencies in radians per second and `values` is an array
of dimensions (n, n, length(frequencies)).

`interpolation` may be `:cubic` (the default) or `:linear`. Cubic
interpolation evaluates the cubic spline through each entry's samples, which
follows a rotating phase far more closely than the chords between them; a
table too short for a cubic takes the highest order it determines, the
parabola through three samples or the line through two. A spline can
overshoot between samples, so a quantity bounded in the data -- a scattering
entry's magnitude, say -- can exceed the bound a little between them; a
bound guaranteed at every frequency takes a passive fit, see
[`RationalScattering`](@ref).

`extrapolation` may be `:error` (the default), `:constant`, `:linear` or
`:zero`. With `:error`, any requested frequency outside the tabulated
range throws an error listing the offending frequencies; extrapolation of
tabulated data is deliberately opt-in. `:constant` holds the end values
beyond the ends, `:linear` continues from them with the interpolant's end
slopes, and `:zero` is zero beyond them, the data of a response known to
vanish outside the band it was sampled over.
"""
struct TabulatedMatrixProvider{T} <: AbstractMatrixProvider
    frequencies::Vector{Float64}
    values::Array{T,3}
    interpolation::Symbol
    extrapolation::Symbol
    # the spline's second derivatives at the knots and its slopes at the two
    # band edges, which the evaluation and the device kernels read; both
    # empty for `:linear`
    curvatures::Array{T,3}
    endslopes::Array{T,3}
end

function TabulatedMatrixProvider(frequencies::AbstractVector,
        values::AbstractArray{T,3};
        interpolation::Symbol = :cubic,
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
    # `issorted` with `<` admits equal neighbors; strictly increasing means
    # each knot exceeds the last, and a repeated knot would put a zero
    # interval under the interpolation
    if !issorted(frequencies; lt = <=)
        throw(ArgumentError("Tabulated frequencies must be strictly increasing."))
    end
    if !(interpolation in (:linear, :cubic))
        throw(ArgumentError(lazy"Unknown interpolation $(interpolation). Supported: :cubic, :linear."))
    end
    if !(extrapolation in (:error, :constant, :linear, :zero))
        throw(ArgumentError(lazy"Unknown extrapolation $(extrapolation). Supported: :error, :constant, :linear, :zero."))
    end
    f = collect(Float64, frequencies)
    vals = Array{T,3}(values)
    n = size(vals, 1)
    if interpolation == :cubic
        curvatures, endslopes = splinecoefficients(f, vals)
    else
        curvatures = Array{T,3}(undef, n, n, 0)
        endslopes = Array{T,3}(undef, n, n, 0)
    end
    return TabulatedMatrixProvider{T}(f, vals, interpolation, extrapolation,
        curvatures, endslopes)
end

# The cubic spline through each entry's samples, held as the spline's second
# derivatives at the knots and its slopes at the two band edges: a piecewise
# cubic with stated knot values and knot second derivatives is exactly the
# spline, and the two arrays are the shape both the evaluation and a device
# kernel index. The spline comes from FastInterpolations; a table of three
# samples takes its parabola and one of two its line, the highest orders
# they determine.
function splinecoefficients(f::Vector{Float64}, values::Array{T,3}) where T
    n = size(values, 1)
    nf = length(f)
    M = zeros(T, n, n, nf)
    slopes = zeros(T, n, n, 2)
    nf == 1 && return M, slopes
    y = Vector{T}(undef, nf)
    for j in 1:n, i in 1:n
        for k in 1:nf
            y[k] = values[i, j, k]
        end
        if nf == 2
            slopes[i, j, 1] = slopes[i, j, 2] = (y[2] - y[1])/(f[2] - f[1])
        elseif nf == 3
            d1 = (y[2] - y[1])/(f[2] - f[1])
            d2 = (y[3] - y[2])/(f[3] - f[2])
            a = (d2 - d1)/(f[3] - f[1])
            for k in 1:3
                M[i, j, k] = 2a
            end
            slopes[i, j, 1] = d1 - a*(f[2] - f[1])
            slopes[i, j, 2] = d2 + a*(f[3] - f[2])
        else
            itp = FastInterpolations.cubic_interp(f, y)
            for k in 1:nf - 1
                c = FastInterpolations.coeffs(itp, (f[k] + f[k + 1])/2)
                M[i, j, k] = 2*c.p[3]
                if k == 1
                    slopes[i, j, 1] = c.p[2]
                end
                if k == nf - 1
                    h = f[nf] - f[nf - 1]
                    M[i, j, nf] = 2*c.p[3] + 6*c.p[4]*h
                    slopes[i, j, 2] = c.p[2] + 2*c.p[3]*h + 3*c.p[4]*h^2
                end
            end
        end
    end
    return M, slopes
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
    # a frequency within roundoff of an end knot is on it: the edge knots
    # of a table assembled from several sources can differ from the
    # frequencies a solver forms by an ulp, and so can a frequency
    # carried to the conjugate ladder and back, or padded by the pump and
    # back. This is the tolerance the coverage test admits (see
    # `tablecovers`), so what a solver takes the table to hold it can
    # evaluate
    edgetol = 8eps(Float64)*max(abs(f[1]), abs(f[end]))
    if p.extrapolation == :error
        outofrange = [w for w in ws if w < f[1] - edgetol || w > f[end] + edgetol]
        if !isempty(outofrange)
            throw(ArgumentError(lazy"The angular frequencies $(outofrange) are outside the tabulated range [$(f[1]), $(f[end])] rad/s. Extrapolation of tabulated data is opt-in: pass extrapolation = :constant or :linear if extrapolation is intended, or fit the block with RationalScattering, which extrapolates as a passive rational function."))
        end
    end
    for i in eachindex(ws)
        w = ws[i]
        if p.extrapolation == :zero && (w < f[1] - edgetol || w > f[end] + edgetol)
            dest[:,:,i] .= zero(T)
            continue
        end
        # a frequency beyond an end knot by that roundoff alone is
        # evaluated at the knot, whatever the block extrapolates by
        if f[1] - edgetol <= w < f[1]
            w = f[1]
        elseif f[end] < w <= f[end] + edgetol
            w = f[end]
        end
        if w <= f[1]
            if p.extrapolation == :constant || length(f) == 1 || w == f[1]
                dest[:,:,i] .= view(p.values,:,:,1)
            elseif p.interpolation == :linear
                # the first segment continues
                lerpslices!(view(dest,:,:,i), p, 1, 2, w)
            else
                edgeslices!(view(dest,:,:,i), p, 1, w)
            end
        elseif w >= f[end]
            if p.extrapolation == :constant || length(f) == 1 || w == f[end]
                dest[:,:,i] .= view(p.values,:,:,length(f))
            elseif p.interpolation == :linear
                lerpslices!(view(dest,:,:,i), p, length(f)-1, length(f), w)
            else
                edgeslices!(view(dest,:,:,i), p, 2, w)
            end
        else
            j = searchsortedlast(f, w)
            if f[j] == w
                dest[:,:,i] .= view(p.values,:,:,j)
            elseif p.interpolation == :linear
                lerpslices!(view(dest,:,:,i), p, j, j+1, w)
            else
                splineslices!(view(dest,:,:,i), p, j, w)
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

# the spline on the segment from knot `j` to knot `j + 1`, from the knot
# values and curvatures
function splineslices!(dest, p::TabulatedMatrixProvider, j::Int, w)
    f1 = p.frequencies[j]
    h = p.frequencies[j+1] - f1
    t = (w - f1)/h
    u = 1 - t
    cu = h^2/6*(u^3 - u)
    ct = h^2/6*(t^3 - t)
    A = view(p.values,:,:,j)
    B = view(p.values,:,:,j+1)
    MA = view(p.curvatures,:,:,j)
    MB = view(p.curvatures,:,:,j+1)
    @. dest = u*A + t*B + cu*MA + ct*MB
    return dest
end

# the tangent beyond a band edge, `e = 1` below the band and `e = 2` above
function edgeslices!(dest, p::TabulatedMatrixProvider, e::Int, w)
    k = e == 1 ? 1 : length(p.frequencies)
    dw = w - p.frequencies[k]
    A = view(p.values,:,:,k)
    s = view(p.endslopes,:,:,e)
    @. dest = A + s*dw
    return dest
end

"""
    PiecewiseTabulatedProvider(tables::Vector{TabulatedMatrixProvider})

Tabulated data in several disjoint bands, each a
[`TabulatedMatrixProvider`](@ref) over its own frequencies, evaluated by
the band a frequency falls in and zero between and beyond them, the
data of a response known only in the bands it was sampled over. The
harmonic transfer functions of a [`LinearizedScattering`](@ref) block built
from a solve are of this kind: the signal band shifted by every
harmonic of the pump, with nothing in between, across which one spline
would swing wildly.
"""
struct PiecewiseTabulatedProvider{T} <: AbstractMatrixProvider
    tables::Vector{TabulatedMatrixProvider{T}}
end
providersize(p::PiecewiseTabulatedProvider) = providersize(p.tables[1])
function evaluateprovider!(dest::AbstractArray{T,3},
        p::PiecewiseTabulatedProvider, ws::AbstractVector) where T
    checkdestsize(dest, providersize(p), length(ws))
    buf = Array{T,3}(undef, providersize(p), providersize(p), 1)
    for i in eachindex(ws)
        w = ws[i]
        dest[:, :, i] .= zero(T)
        for t in p.tables
            f = t.frequencies
            edgetol = 8eps(Float64)*max(abs(f[1]), abs(f[end]))
            if f[1] - edgetol <= w <= f[end] + edgetol
                evaluateprovider!(buf, t, [clamp(w, f[1], f[end])])
                dest[:, :, i] .= view(buf, :, :, 1)
                break
            end
        end
    end
    return dest
end
# the frequencies of every band of a piecewise table, ascending
piecewisefrequencies(p::PiecewiseTabulatedProvider) = sort!(vcat([t.frequencies for t in p.tables]...))

# tabulated data in disjoint bands: the knots split where a gap exceeds
# `gapfactor` times the median spacing, one table per band, each with
# the interpolation given and zero beyond it
function piecewisetable(nus::Vector{Float64}, values::Array{T,3}; interpolation::Symbol = :cubic,
        gapfactor::Real = 10) where T
    n = length(nus)
    spacings = diff(nus)
    median = isempty(spacings) ? 0.0 : sort(spacings)[(length(spacings) + 1) ÷ 2]
    tables = TabulatedMatrixProvider{T}[]
    start = 1
    for k in 1:n
        if k == n || spacings[k] > gapfactor*median
            push!(tables, TabulatedMatrixProvider(nus[start:k], values[:, :, start:k];
                interpolation = interpolation, extrapolation = :zero))
            start = k + 1
        end
    end
    return PiecewiseTabulatedProvider(tables)
end

function checkdestsize(dest, n::Int, nf::Int)
    if size(dest) != (n, n, nf)
        throw(DimensionMismatch(lazy"The destination array has size $(size(dest)) but ($(n), $(n), $(nf)) is required."))
    end
    return nothing
end

"""
    matrixprovider(x, T; n = nothing, interpolation = :cubic,
        extrapolation = :error, form = :matrix)

Normalize user input into an [`AbstractMatrixProvider`](@ref) with element
type `T`: a matrix becomes a [`ConstantMatrixProvider`](@ref), a tuple
`(frequencies, values)` a [`TabulatedMatrixProvider`](@ref), a callable of
angular frequency a [`CallableMatrixProvider`](@ref) (which requires the
dimension `n` and accepts `form`), and an existing provider is returned as
is. `n`, when given, is checked against the data.
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
        ::Type{T}; n = nothing, interpolation::Symbol = :cubic,
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

A noise model for a [`ScatteringParameters`](@ref) asserting that its
scattering matrix is unitary at every frequency, or for a
[`LinearizedScattering`](@ref) that its multi-mode scattering matrix `S`
satisfies `J - S J S' = 0` with `J` the signs of the mode frequencies,
so that the block absorbs nothing and adds no noise.

With the default [`Passive`](@ref) model a block gets noise channels unless
its data can be shown to be unitary, which is possible for a constant
matrix and a table but not for a callable, whose values away from the
evaluated frequencies are unknown. A lossless callable would therefore
carry noise channels which are identically zero, at a real cost on a
circuit with many blocks; `Lossless()` removes them.

For a constant or tabulated block the assertion is checked at
construction, with a fixed tolerance of `1e-10` rather than the block's
`atol`, and is an error if false. For a callable it is taken on trust,
and asserting it of a block which does absorb omits the noise the block
should add, making the quantum efficiency and the commutation relations
wrong by that much. A pumped block is checked when it is built and
again over the modes of every solve which evaluates it, to the block's
`atol`; a fit of a pumped block is not lossless to better than its
error, and states the noise it needs instead (see
[`NoiseCovariance`](@ref)).
"""
struct Lossless end

# === the zero frequency model of a scattering block ===
#
# The direct current behavior of a block is its scattering matrix at zero
# frequency, and by default the block's own data is evaluated there. That
# fails in two ways: tabulated data may not extend down to zero, and a
# closed form may have a limit which exists but is not evaluable at zero
# (a series capacitor written `1/(im*w*C)` is an open circuit at direct
# current, but the expression is infinite there). For those the zero
# frequency model is stated separately. Every realizable choice is some
# real `S(0)`, so `ScatteringDC` carries the general case and the common
# ones are named.

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
are shorthand for. It is validated for passivity when the containing
[`ScatteringParameters`](@ref) is constructed, with that block's `atol`, on
the same terms as the block's own data, so an active zero frequency model
needs the same `NoiseCovariance` declaration an active block does.
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
    dcscatteringmatrix(m::AbstractDCModel, n)

The `n` by `n` real zero frequency scattering matrix of a stated model:
[`OpenDC`](@ref), [`ShortDC`](@ref), [`ThroughDC`](@ref) or
[`ScatteringDC`](@ref). [`ScatteringLimit`](@ref) states none and has no
method here; it is handled by evaluating the block's own data. Throws for a
[`ThroughDC`](@ref) with `n != 2` or a [`ScatteringDC`](@ref) of the wrong
dimension.
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
overrides the `temperature` argument of the analysis. At zero temperature
it coincides with [`Passive`](@ref).
"""
struct ThermalEquilibrium{T}
    temperature::T
end

"""
    NoiseCovariance(V; interpolation = :cubic, extrapolation = :error,
        atol = 1e-8, completed = false, padding = 4)

The noise a [`ScatteringParameters`](@ref) adds, stated outright rather
than derived from its loss, which is how an active block, an amplifier
given by its scattering parameters, declares its noise. `V` is the
symmetrized covariance of the noise wave the block emits at its ports, in
the units of the rest of the noise outputs, where a vacuum channel counts
as one and a channel at temperature `T` as `coth(hbar*w/(2*k*T))`, that is
`2*nbar + 1`: the same units as `Cnoise`. It may be a matrix, a callable
of angular frequency, or a tuple `(frequencies, values)` of tabulated
data, following the same provider forms as the scattering data, and is
evaluated with the block's negative frequency rule.

The block's output must obey the commutation relations, so with
`K = I - S S'` the added noise has the commutator `K`, and a covariance is
realizable only when `V - K` and `V + K` are both positive semidefinite.
`(V + K)/2` is then the covariance of channels which emit like a mode in
its vacuum and `(V - K)/2` of channels which emit like the conjugate of
one, the idler channels of an amplifier; a passive block in thermal
equilibrium is the case `V = coth(hbar*w/(2*k*T)) K`. A phase insensitive
amplifier of power gain `G` from port 1 to port 2 has `K[2,2] = 1 - G`, so
`V[2,2] >= G - 1`, the noise of a quantum limited amplifier, and one with
an input referred added noise of `nadd` photons has `V[2,2] = 2*G*nadd`;
`V[1,1]` is what it emits backward out of its input, `(2*nbar + 1)*K[1,1]`
for an input matched at its physical temperature.

The condition is checked on the eigenvalues to `atol` at construction
where both the scattering data and `V` are stored, and at every frequency
a solver evaluates the block at. `V` is also validated for Hermitian
symmetry. A block with this model carries no temperature: its noise is
`V`, whatever the analysis temperature.

With `completed = true` the covariance is completed to the commutation
relations rather than held to them: `V` is replaced by
`V + neg(V - K) + neg(V + K)`, `neg` taking the negative part of a
Hermitian matrix, the sum of `-lambda v v'` over its negative
eigenvalues, which makes `V - K` and `V + K` positive semidefinite and
adds nothing where they are. For `V = 0` the addition is `|K|`, the
least total noise a Gaussian channel with the block's map can add, the
`Ymin` of the quantum optics functions in the basis of the modes; for
a stated `V` it is a sufficient addition, the least being a
semidefinite program with no closed form. The block then adds the
noise it states and what more its commutator requires, and its output
obeys the commutation relations exactly, whatever `V` and `S` are. An
ordinary block is completed at each frequency. A pumped block's
covariance spans the modes of a solve at once, and the negative part
of a matrix is not that of its parts, so it is completed over the
ladder of the modes padded by `padding` multiples of the pump
frequency on either side, with every input which feeds it, and
restricted to the modes of the solve (see
[`completedcovariance`](@ref)): the block's noise is then one model
whatever modes a solve keeps, to the precision of the padding, which
converges since a fit's conversion vanishes at high frequency, and the
inputs a solve lacks are traced out in their vacuum. `padding` is that
precision and not a bound on it, so a device whose conversion reaches
far is checked by raising it and comparing the covariance the solve
keeps. Every mode a solve
asks for is completed, including one beyond the sidebands the block's
data holds, which it scatters nothing at and so carries the vacuum its
commutator requires, as its stamp takes it. This is how a fit
of a pumped block states its noise (see [`RationalScattering`](@ref)):
a fit is neither lossless nor consistent with a stated covariance to
better than its error, and the completion turns that error into noise
the block emits, where a tolerance would only excuse it.
"""
struct NoiseCovariance{P}
    provider::P
    interpolation::Symbol
    extrapolation::Symbol
    atol::Float64
    completed::Bool
    # the multiples of the pump frequency a pumped block's ladder is
    # padded by on either side of the modes of a solve before its
    # completion
    padding::Int
end
function NoiseCovariance(V; interpolation::Symbol = :cubic,
        extrapolation::Symbol = :error, atol::Real = 1e-8, completed::Bool = false,
        padding::Integer = 4)
    (isfinite(atol) && atol >= 0) || throw(ArgumentError("atol must be finite and nonnegative."))
    padding >= 0 || throw(ArgumentError("padding must be nonnegative."))
    return NoiseCovariance(V, interpolation, extrapolation, Float64(atol), completed, Int(padding))
end

# === scattering block ===

"""
    ScatteringParameters(S; nports = nothing, zref = nothing, grounded = true,
        noise = Passive(), negative_frequency = ConjugateSymmetry(),
        interpolation = :cubic, extrapolation = :error, form = :matrix,
        derivatives = NamedTuple(), dcmodel = ScatteringLimit(),
        atol = 1e-8)

A multiport component defined by its scattering parameters. `S` may be:

- a constant matrix;
- a callable of angular frequency returning a matrix (requires `nports`);
- a tuple `(frequencies, values)` of tabulated data with `frequencies` in
  radians per second and `values` of size (nports, nports, nfrequencies);
- a path to a Touchstone file, from which the reference impedance is also
  read.

Tabulated data is interpolated with the cubic spline through each entry's
samples (`interpolation = :cubic`, the default; `:linear` takes the chords
between them, which lag a rotating phase) and is never extrapolated unless
asked: `extrapolation` is `:error` by default, with `:constant` and
`:linear` to opt in. The harmonic balance solvers evaluate a block wherever
their mixing products fall, which can be far outside the band the data
covers, so measured data meant for them is better fitted with
[`RationalScattering`](@ref), which extrapolates as a passive rational
function and is passive at every frequency by construction, where any
interpolant of passive samples can stray between and beyond them.

`zref` is the reference impedance in Ohms, a scalar broadcast to all ports
or a vector with one entry per port. For a Touchstone file the reference
impedance is read from the file, and a `zref` which disagrees with it is
an error, since the intent is ambiguous between correcting a mislabeled
file and requesting renormalization. Scattering data is used at its native
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
efficiency and the commutation relations account for. A `NoiseCovariance`
states the added noise outright, which is what an active block, an
amplifier given by its scattering parameters, needs, and is held to the
minimum the commutation relations require of it. `Lossless` is how a
unitary callable says so, since unlike stored data it cannot be checked.
`negative_frequency` is [`ConjugateSymmetry`](@ref) (default) or
[`Native`](@ref). `form` says how a callable `S` is called; see
[`CallableMatrixProvider`](@ref). Passivity of constant and tabulated
scattering data is validated at construction with absolute tolerance
`atol` unless the noise model is a `NoiseCovariance`, which permits
active blocks.

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
scattering matrix and is never passivity checked. `derivatives` and `form`
apply when `S` is given as data or a callable; a Touchstone path ignores
them.

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

# positional constructors without derivatives and without a stated zero
# frequency model
ScatteringParameters(provider, nports::Int, zref::Vector{Float64},
    grounded::Bool, noise, negative_frequency) =
    ScatteringParameters(provider, nports, zref, grounded, noise,
        negative_frequency, NamedTuple(), ScatteringLimit())
ScatteringParameters(provider, nports::Int, zref::Vector{Float64},
    grounded::Bool, noise, negative_frequency, derivatives) =
    ScatteringParameters(provider, nports, zref, grounded, noise,
        negative_frequency, derivatives, ScatteringLimit())

function ScatteringParameters(S; nports = nothing, zref = nothing,
        grounded::Bool = true, noise = Passive(),
        negative_frequency = ConjugateSymmetry(),
        interpolation::Symbol = :cubic, extrapolation::Symbol = :error,
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
    # omitted, the reference impedance is 50 Ohms at every port
    zrefvec = zrefvector(something(zref, 50.0), n)
    if !(noise isa NoiseCovariance)
        checkpassive(provider; atol = atol)
    end
    checklossless(noise, provider)
    noise = preparenoise(noise, provider, n)
    dprov = NamedTuple(k => begin
            dp = matrixprovider(v, Complex{Float64}; n = n, form = form)
            providersize(dp) == n || throw(DimensionMismatch(lazy"the derivative for parameter $(k) has dimension $(providersize(dp)) but the block has $(n) ports."))
            dp
        end for (k, v) in pairs(derivatives))
    checkdcmodel(dcmodel, n, noise, atol)
    return ScatteringParameters(provider, n, zrefvec, grounded, noise,
        negative_frequency, dprov, dcmodel)
end

# A stated zero frequency matrix is checked like the block's own data: its
# size against the block, and passivity unless the block declared itself
# active with a `NoiseCovariance`.
checkdcmodel(::ScatteringLimit, n::Int, noise, atol) = nothing
function checkdcmodel(m::AbstractDCModel, n::Int, noise, atol)
    S0 = dcscatteringmatrix(m, n)
    if !(noise isa NoiseCovariance)
        sv = maximum(svdvals(S0); init = 0.0)
        sv <= 1 + atol || throw(ArgumentError(lazy"the zero frequency scattering matrix has largest singular value $(sv), so it is active at direct current. An active block has to declare its own noise with NoiseCovariance, as it does at every other frequency."))
    end
    return nothing
end

# the reference impedances as a vector of `n` positive finite reals
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

"""
    preparenoise(noise, provider, n)

The noise model of a block as it is stored: a [`NoiseCovariance`](@ref)
with its data as a matrix provider of the block's dimension `n`, validated
for Hermitian symmetry and, where both it and the scattering `provider`
are stored data, for the minimum noise the commutation relations require
(see [`quantumnoisemargin`](@ref)); any other model unchanged.
"""
function preparenoise(noise::NoiseCovariance, provider, n::Int)
    vp = matrixprovider(noise.provider, Complex{Float64}; n = n,
        interpolation = noise.interpolation,
        extrapolation = noise.extrapolation)
    if providersize(vp) != n
        throw(DimensionMismatch(lazy"The noise covariance dimension $(providersize(vp)) does not match the number of ports $(n)."))
    end
    checkhermitian(vp)
    prepared = NoiseCovariance(vp, noise.interpolation, noise.extrapolation,
        noise.atol, noise.completed, noise.padding)
    checkquantum(prepared, provider, n)
    return prepared
end
preparenoise(noise, provider, n::Int) = noise

"""
    quantumnoisemargin(V::AbstractMatrix, S::AbstractMatrix)

The smallest eigenvalue of `V - K` and of `V + K`, with `K = I - S S'`,
which is nonnegative when the noise covariance `V` is one a block with
the scattering matrix `S` can add without violating the commutation
relations: `(V + K)/2` and `(V - K)/2` are then the covariances of its
channels of either kind (see [`NoiseCovariance`](@ref)).
"""
function quantumnoisemargin(V::AbstractMatrix, S::AbstractMatrix)
    n = size(S, 1)
    K = Matrix{Complex{Float64}}(I, n, n) - S*S'
    Vh = Matrix{Complex{Float64}}(V)
    return min(minimum(real.(eigvals(Hermitian(Vh - K)))),
        minimum(real.(eigvals(Hermitian(Vh + K)))))
end

# The frequencies at which stored data of both providers can be
# compared: the knots of whichever is tabulated, within the range of the
# other. A callable is unknown between evaluations, so a pair with one is
# checked only where a solver evaluates it.
function storedfrequencies(a, b)
    ta = a isa TabulatedMatrixProvider
    tb = b isa TabulatedMatrixProvider
    ta || tb || return Float64[0.0]
    fs = Float64[]
    ta && append!(fs, a.frequencies)
    tb && append!(fs, b.frequencies)
    inrange(p, f) = !(p isa TabulatedMatrixProvider) ||
        p.extrapolation != :error || (p.frequencies[1] <= f <= p.frequencies[end])
    return unique!(sort!([f for f in fs if inrange(a, f) && inrange(b, f)]))
end

# whether the check can run on data alone
storedprovider(p) = p isa ConstantMatrixProvider || p isa TabulatedMatrixProvider

function checkquantum(noise::NoiseCovariance, provider, n::Int)
    # a covariance completed to the commutation relations meets them by
    # construction
    noise.completed && return nothing
    storedprovider(provider) && storedprovider(noise.provider) || return nothing
    fs = storedfrequencies(provider, noise.provider)
    isempty(fs) && return nothing
    S = Array{Complex{Float64},3}(undef, n, n, length(fs))
    V = Array{Complex{Float64},3}(undef, n, n, length(fs))
    evaluateprovider!(S, provider, fs)
    evaluateprovider!(V, noise.provider, fs)
    for k in eachindex(fs)
        margin = quantumnoisemargin(view(V, :, :, k), view(S, :, :, k))
        if margin < -noise.atol
            throw(ArgumentError(lazy"The noise covariance is less than the commutation relations require of this block: at $(fs[k]) rad/s the smallest eigenvalue of V - K or V + K, with K = I - S S', is $(margin). An amplifier of power gain G has to emit at least G - 1 at its output; see NoiseCovariance."))
        end
    end
    return nothing
end

# whether the provider's matrices are Hermitian to `atol` at every sample, the
# validation of a `NoiseCovariance`
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

# whether the provider's matrices are passive to `atol` at every sample
# (no singular value above one), the validation of a `ScatteringParameters`
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
    RationalScatteringProvider(A, B, C, D)

A scattering matrix given as the real state space realization
`S(s) = D + C (s I - A)^(-1) B` of a passive rational multiport, the
form a vector fit of measured or simulated scattering data takes and the
one the transient solver realizes in time, `dz/dt = A z + B a`,
`b = C z + D a` on the incident and reflected power waves. Evaluated at
a signed angular frequency by one dense solve, without forming an
inverse. Built by [`RationalScattering`](@ref), which validates it.
"""
struct RationalScatteringProvider <: AbstractMatrixProvider
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
end
providersize(p::RationalScatteringProvider) = size(p.D, 1)
function evaluateprovider!(dest::AbstractArray{T,3},
        p::RationalScatteringProvider, ws::AbstractVector) where T
    n = size(p.D, 1)
    checkdestsize(dest, n, length(ws))
    nz = size(p.A, 1)
    if nz == 0
        for i in eachindex(ws)
            dest[:, :, i] .= p.D
        end
        return dest
    end
    # over many frequencies the resolvent is taken through the Schur form
    # of the state matrix, `A = Q T Q'` with `T` upper triangular, once,
    # each frequency then costing a triangular solve of the states squared
    # rather than a factorization of the states cubed, and as accurately,
    # the Schur vectors being unitary whatever the conditioning of the
    # eigenvectors, which a nearly defective realization, an all pass, has
    # no usable set of
    if T <: Complex && length(ws) > 8
        F = schur(Matrix{Complex{Float64}}(p.A))
        QB = F.Z'*p.B
        CQ = p.C*F.Z
        M = similar(F.T)
        X = similar(QB)
        for i in eachindex(ws)
            M .= .-F.T
            for k in axes(M, 1)
                M[k, k] += im*ws[i]
            end
            X .= QB
            ldiv!(UpperTriangular(M), X)
            dest[:, :, i] .= p.D .+ CQ*X
        end
        return dest
    end
    for i in eachindex(ws)
        dest[:, :, i] .= p.D .+ p.C*((im*ws[i]*I - p.A) \ p.B)
    end
    return dest
end

# The passivity of a rational realization: the bounded real lemma's
# Hamiltonian test where the feedthrough is strictly contractive, exact
# and dependency free, since `S` has a singular value crossing one on the
# imaginary axis if and only if the Hamiltonian matrix has an imaginary
# eigenvalue; and where the feedthrough reaches one, as a lossless block's
# does, a fine sampling over the band of the poles.
function checkpassive(p::RationalScatteringProvider; atol = 1e-8)
    A, B, C, D = p.A, p.B, p.C, p.D
    margin = passivitymargin(D)
    margin < -atol && throw(ArgumentError(lazy"The rational scattering block is not passive at infinite frequency: the minimum eigenvalue of I - D*D' is $(margin)."))
    size(A, 1) == 0 && return nothing
    verdict, worst, level, w = passivityassessment(A, B, C, D; atol = atol)
    verdict === :active && throw(ArgumentError(lazy"The rational scattering block is not passive: its largest singular value over all frequencies is at least $(worst), at $(w) rad/s, which is above the tolerance $(atol)."))
    return nothing
end

# the realization in the frequency unit of the poles, `S(s) = D + C (s/w I - A/w)^(-1) B/w`,
# with each state scaled to balance its input row and output column, the
# similarity `A -> T^(-1) A T`, `B -> T^(-1) B`, `C -> C T`, which leaves
# `S` alone; and the scale
function balancedrealization(A, B, C)
    nz = size(A, 1)
    wscale = max(maximum(abs, eigvals(A)), floatmin(Float64))
    An, Bn = A ./ wscale, B ./ wscale
    t = [(cb = norm(view(C, :, k)); bb = norm(view(Bn, k, :)); cb > 0 && bb > 0 ? sqrt(bb/cb) : cb > 0 ? 1/cb : bb > 0 ? bb : 1.0) for k in 1:nz]
    return (1 ./ t) .* An .* transpose(t), Bn ./ t, C .* transpose(t), wscale
end

"""
    hinfnorm(A, B, C, D; rtol = 1e-8, span = 8.0, refinements = 32,
        pad = 10.0)

The largest singular value of the real rational matrix
`S(s) = D + C (s I - A)^(-1) B` over every frequency, by the level set
iteration of Boyd, Balakrishnan, Bruinsma and Steinbuch: a lower bound
from the feedthrough, samples spanning the poles' frequencies, and a
peak search around each pole is raised by a hair to a level, the
frequencies where a singular value equals the level are the imaginary
eigenvalues of the pencil of [`passivitycrossings`](@ref) for `S` over
the level, and the largest singular value between consecutive ones
raises the bound, until no singular value reaches the level. A peak
however narrow is found, since the pencil finds every crossing of the
level, which a sample can miss.

Returns three values: a lower bound on the norm, the frequency in
rad/s where it was attained, and the level the search established
nothing reaches, `Inf` where termination established no such level.
The first value is only the largest value the search evaluated; what
termination proves is the third. They differ by `2 rtol`, which
matters wherever the answer is compared against one: at the default
tolerance a norm returned as `1 - 1e-9` is consistent with a true norm
of `1 + 1e-8`, so calling a block passive on the first value is
calling it passive on a lower bound.

`rtol` is not worth pushing far below its default. The level is a
bound only so far as the pencil resolves the crossings of it, and a
peak which exceeds a level by less than roundoff brings its two
crossings together into a nearly double eigenvalue which leaves the
imaginary axis and is lost: a tolerance below the square root of `eps`
asks the pencil for a resolution it does not have, and returns a level
which is not a bound.

`span` is how many pole half widths the peak search brackets either
side of each complex pole, `refinements` how many golden section steps
refine each bracket, and `pad` how far past the outermost pole
magnitudes the probe grid extends.
"""
function hinfnorm(A, B, C, D; rtol = 1e-8, span::Real = 8.0, refinements::Int = 32,
        pad::Real = 10.0)
    An, Bn, Cn, wscale = balancedrealization(A, B, C)
    # A dense solve at every evaluation rather than held Schur
    # factors: the search refines around each pole, which is where a
    # triangular solve of `i w I - T` is at its worst, the diagonal
    # entry it divides by nearly zero with nothing to pivot, and the
    # peak of a narrow resonance depends on those digits. The dense
    # solve pivots, and the accuracy buys the answer.
    S = w -> try
        D .+ Cn*((im*w*I - An) \ Bn)
    catch e
        e isa SingularException ? fill(Inf, size(D)) : rethrow()
    end
    λs = eigvals(An)
    # The grid spans the magnitudes of the poles the system has, not
    # fixed decades around one: a peak can lie far from every pole
    # frequency, as `k s/((s + a)(s + b))` peaking at `sqrt(a b)` shows,
    # so the grid has to cover the whole range the poles set.
    mags = [abs(l) for l in λs if abs(l) > 0]
    glo = isempty(mags) ? 1e-2 : minimum(mags)/pad
    ghi = isempty(mags) ? 1e2 : maximum(mags)*pad
    probes = vcat(0.0, [abs(imag(l)) for l in λs if abs(imag(l)) > 0],
                  exp.(range(log(glo), log(ghi); length = max(9, 2*length(λs)))))
    sort!(probes)
    bound, where = opnorm(D), Inf
    for w in probes
        s = opnorm(S(w))
        isfinite(s) || return Inf, w*wscale, Inf
        s > bound && ((bound, where) = (s, w))
    end
    isfinite(bound) || return Inf, where*wscale, Inf
    # A golden section search sharpens the lower bound: a resonance
    # peaks near its pole's frequency but not at it, and for a narrow
    # one the difference exceeds the tolerance being tested against.
    # This only ever raises a lower bound, so it cannot make the answer
    # wrong, and the bound is what decides passivity, since the level
    # cannot be sharpened past the pencil's resolution (see the
    # docstring). The brackets are a few half widths either side of
    # each complex pole, plus the span between the neighbours of the
    # best probe, which covers a peak the grid only straddled -- a real
    # pole has no resonance of its own, and a peak it takes part in
    # need not be near it, so the second bracket is the one that covers
    # it.
    brackets = Tuple{Float64,Float64}[]
    for l in λs
        imag(l) > 0 || continue
        w0, half = imag(l), max(abs(real(l)), eps())
        push!(brackets, (max(w0 - span*half, 0.0), w0 + span*half))
    end
    if isfinite(where)
        j = searchsortedfirst(probes, where)
        push!(brackets, (probes[max(j - 1, 1)], probes[min(j + 1, length(probes))]))
    end
    for (a, b) in brackets
        b > a || continue
        φ = (sqrt(5) - 1)/2
        u, v = b - φ*(b - a), a + φ*(b - a)
        fu, fv = opnorm(S(u)), opnorm(S(v))
        for _ in 1:refinements
            if fu > fv
                b, v, fv = v, u, fu
                u = b - φ*(b - a); fu = opnorm(S(u))
            else
                a, u, fu = u, v, fv
                v = a + φ*(b - a); fv = opnorm(S(v))
            end
        end
        isfinite(fu) && isfinite(fv) || return Inf, a*wscale, Inf
        fu > bound && ((bound, where) = (fu, u))
        fv > bound && ((bound, where) = (fv, v))
    end
    # the level the search ends on, once nothing reaches it
    ceiling = Inf
    for iteration in 1:100
        level = (1 + 2rtol)*bound
        # the pencil at a level above the bound is regular, a singular
        # value equal to the level everywhere being impossible, so its
        # eigenvalues are taken as they come: a nearly lossless block only
        # scales the pencil's determinant by a small constant, and a
        # singular value test on it would mistake that for singularity and
        # miss a peak between the samples
        crossings, _ = pencilcrossings(An, Bn, Cn ./ level, D ./ level; singulartest = false)
        if length(crossings) < 2
            # no singular value reaches the level, anywhere
            ceiling = level
            break
        end
        raised = false
        for k in 1:length(crossings) - 1
            w = (crossings[k] + crossings[k + 1])/2
            s = opnorm(S(w))
            isfinite(s) || return Inf, w*wscale, Inf
            s > bound*(1 + rtol) && ((bound, where, raised) = (s, w, true))
        end
        # a level with crossings is reached somewhere, so midpoints which
        # cannot improve on the bound are a failure of the search rather
        # than a proof, and nothing is established
        raised || break
    end
    return bound, where*wscale, ceiling
end

"""
    passivityassessment(A, B, C, D; atol = 1e-8, rtol = 1e-8)

Whether the real rational block `S(s) = D + C (s I - A)^(-1) B` is
passive to within `atol`, as `(:passive, :active, :indeterminate)`,
together with the lower bound on its largest singular value over all
frequencies, the level the search ended at, and the frequency in rad/s
where the lower bound was attained.

Three answers rather than two, because [`hinfnorm`](@ref) returns two
numbers which straddle the truth and the question can fall between them.
A block is `:active` when even the lower bound exceeds `1 + atol`, which
settles it; `:passive` when the level does not, which also settles it,
the level being what the search's termination establishes. In between
nothing is settled: the block may be passive or may not, and the caller
is told so rather than given whichever bound suits.

The middle case is not rare. A lossless block has a largest singular
value of exactly one, so its level stands at `1 + 2 rtol` and it is
`:indeterminate` at any `atol` below that. Callers which must decide
regardless decide on the lower bound, and this makes that choice
explicit rather than implicit in a two valued answer.
"""
function passivityassessment(A, B, C, D; atol = 1e-8, rtol = 1e-8)
    lower, where, level = hinfnorm(A, B, C, D; rtol = rtol)
    verdict = if !isfinite(level) || !isfinite(lower)
        lower > 1 + atol ? :active : :indeterminate
    elseif lower > 1 + atol
        :active
    elseif level <= 1 + atol
        :passive
    else
        :indeterminate
    end
    return verdict, lower, level, where
end

"""
    passivitycrossings(A, B, C, D)

The frequencies in rad/s at which a singular value of the real rational
scattering matrix `S(s) = D + C (s I - A)^(-1) B` equals one, sorted,
and the scale of the poles: the finite eigenvalues on the imaginary axis
of the pencil of the equations `i w x = A x + B u`, `-i w y = A' y + C' w`,
`w = C x + D u`, `u = B' y + D' w`, which say `S(i w)' S(i w) u = u`,
whose matrices are formed without inverting `I - D' D`, so a feedthrough
on the unit circle is no obstacle. Returns `nothing` for the crossings
when the pencil is singular, which is when a singular value is one at
every frequency, as a lossless block's are. The pencil for `S` over a
level finds the crossings of that level, which is how [`hinfnorm`](@ref)
finds the largest singular value.
"""
function passivitycrossings(A, B, C, D)
    An, Bn, Cn, wscale = balancedrealization(A, B, C)
    crossings, _ = pencilcrossings(An, Bn, Cn, D)
    return isnothing(crossings) ? nothing : crossings .* wscale, wscale
end

# the crossings of one in the units of the balanced realization
function pencilcrossings(An, Bn, Cn, D; singulartest::Bool = true)
    nz, m = size(An, 1), size(D, 1)
    Z = zeros
    H = [An Z(nz, nz) Bn Z(nz, m); Z(nz, nz) transpose(An) Z(nz, m) transpose(Cn);
        Cn Z(m, nz) D -Matrix(1.0I, m, m); Z(m, nz) transpose(Bn) -Matrix(1.0I, m, m) transpose(D)]
    E = Matrix(Diagonal(vcat(ones(nz), -ones(nz), zeros(2m))))
    # a singular pencil, one whose determinant vanishes at every point,
    # has no meaningful eigenvalues: it is told by its rank at two points
    # off the axis
    if singulartest
        for l0 in (complex(0.7, 1.3), complex(-1.1, 0.4))
            sv = svdvals(l0 .* E .- H)
            sv[end] <= 1e-10*sv[1] && return nothing, 1.0
        end
    end
    lambda = eigvals(H, E)
    finite = [l for l in lambda if isfinite(real(l)) && isfinite(imag(l)) && abs(l) <= 1e8]
    crossings = sort!([abs(imag(l)) for l in finite if abs(real(l)) <= 1e-8*(abs(l) + 1)])
    # a crossing and its conjugate are one
    merged = Float64[]
    for w in crossings
        (isempty(merged) || w - last(merged) > 1e-9*(w + 1)) && push!(merged, w)
    end
    return merged, 1.0
end

"""
    RationalScattering(A, B, C, D; zref = 50.0, grounded = true,
        noise = Passive(), atol = 1e-8)

A [`ScatteringParameters`](@ref) block from the real state space
realization `S(s) = D + C (s I - A)^(-1) B` of a passive rational
multiport, with `A` the `nz` by `nz` state matrix, `B` `nz` by `nports`,
`C` `nports` by `nz` and `D` `nports` by `nports`, all real and finite,
`A` stable. The block is validated as passive by the bounded real
lemma's Hamiltonian test, or by sampling where its feedthrough is
lossless, and it is rejected otherwise, unless it states its noise with
a [`NoiseCovariance`](@ref), which is how an active block, an amplifier
given by its scattering parameters, declares it; stability is required
of every realization. It is evaluated by the harmonic balance solvers
at every frequency and realized in time by the transient solver with
its states, so the two describe the same block, and its noise is the
noise of its loss at every frequency by Bosma's relation, or the noise
it states. The realization is what a vector fit of measured or
simulated data delivers; a lossless line is [`TransmissionLine`](@ref)
instead, which needs no states.
"""
function RationalScattering(A, B, C, D; zref = 50.0, grounded::Bool = true, noise = Passive(), atol::Real = 1e-8)
    Am, Bm, Cm, Dm = Matrix{Float64}(A), Matrix{Float64}(B), Matrix{Float64}(C), Matrix{Float64}(D)
    n, nz = size(Dm, 1), size(Am, 1)
    size(Dm) == (n, n) && size(Am) == (nz, nz) && size(Bm) == (nz, n) && size(Cm) == (n, nz) || throw(DimensionMismatch(
        lazy"the realization needs A of size (nz, nz), B (nz, nports), C (nports, nz) and D (nports, nports); got $(size(Am)), $(size(Bm)), $(size(Cm)), $(size(Dm))."))
    all(M -> all(isfinite, M), (Am, Bm, Cm, Dm)) || throw(ArgumentError("the realization must be finite."))
    nz == 0 || maximum(real.(eigvals(Am))) < 0 || throw(ArgumentError(
        lazy"the realization is unstable: the largest real part of an eigenvalue of A is $(maximum(real.(eigvals(Am)))) per second."))
    provider = RationalScatteringProvider(Am, Bm, Cm, Dm)
    noise isa NoiseCovariance || checkpassive(provider; atol = atol)
    checklossless(noise, provider)
    noise = preparenoise(noise, provider, n)
    z = zref isa Number ? fill(Float64(zref), n) : Float64.(collect(zref))
    length(z) == n && all(x -> isfinite(x) && x > 0, z) || throw(ArgumentError("give one positive reference impedance per port."))
    return ScatteringParameters(provider, n, z, grounded, noise, ConjugateSymmetry())
end

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
# A rational block is lossless when every singular value of `S` is one at
# every frequency: the largest at most one, by the largest singular value
# over all frequencies, and the smallest at least one, by the largest
# singular value of `S^(-1)` over all frequencies, which is a rational
# block of its own when the feedthrough is invertible, and without an
# invertible feedthrough the block is not lossless. Both are found by the
# level set iteration, which finds a peak or a notch however narrow. No
# test of the coefficients of `I - S(-s)' S(s)` can: a notch of relative
# width `eps` has coefficients of order `eps^2` there and reaches one
# at its center, so only the values on the axis tell. The tolerance is
# a part in 1e8 of a singular value, a loss no noise sees.
# Nothing is inferred about a rational block by default: it keeps the
# channels of its loss, however small, and only a declaration of
# `Lossless()` is validated, by `losslessnorms`, which can only refuse.
provablylossless(p::RationalScatteringProvider; atol = 1e-8) = false
function losslessnorms(p::RationalScatteringProvider; atol = 1e-8)
    nz = size(p.A, 1)
    unitaritydeviation(p.D) <= atol || return false
    nz == 0 && return true
    hinfnorm(p.A, p.B, p.C, p.D)[1] <= 1 + atol || return false
    Dinv = inv(p.D)
    return hinfnorm(p.A - p.B*Dinv*p.C, p.B*Dinv, -Dinv*p.C, Dinv)[1] <= 1 + atol
end
provablylossless(p::ConstantMatrixProvider; atol = 1e-10) =
    unitaritydeviation(p.A) <= atol
function provablylossless(p::TabulatedMatrixProvider; atol = 1e-10)
    # linear extrapolation leaves the stored range unbounded, so nothing
    # can be proved about it, and zero beyond the band absorbs everything
    (p.extrapolation == :linear || p.extrapolation == :zero) && return false
    return worstunitaritydeviation(p) <= atol
end
provablylossless(b::ScatteringParameters; atol = 1e-10) =
    provablylossless(b.provider; atol = atol)

# `Lossless` asserts unitarity at every frequency, which stored data can be
# checked for and a callable cannot
function checklossless(noise::Lossless, provider)
    if provider isa CallableMatrixProvider
        return nothing
    end
    if provider isa RationalScatteringProvider
        losslessnorms(provider) || throw(ArgumentError("noise = Lossless() says the scattering matrix is unitary at every frequency, but this rational block's largest singular value over all frequencies, or its inverse's, exceeds one by more than a part in 1e8. Use the default Passive() noise model, which gives a dissipative block the noise its loss requires."))
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
# The interpolant between two knots is (1-t)*A + t*B, whose deviation from
# unitarity with A and B unitary is t*(1-t)*(A*B' + B*A' - 2I): largest at
# the midpoint, so the knots and the midpoints together bound the whole
# table. Two unitary knots do not make a unitary interpolant (S = 1 and
# S = -1 interpolate to a perfect absorber halfway between).
function worstunitaritydeviation(p::TabulatedMatrixProvider)
    worst = 0.0
    for k in axes(p.values,3)
        worst = max(worst, unitaritydeviation(view(p.values,:,:,k)))
    end
    for k in 1:size(p.values,3)-1
        mid = (view(p.values,:,:,k) .+ view(p.values,:,:,k+1))./2
        if p.interpolation == :linear
            worst = max(worst, unitaritydeviation(mid))
        else
            # The spline is the chord plus the curvature correction
            # `(h^2/6)*((u^3 - u)*M_k + (t^3 - t)*M_(k+1))`, whose operator
            # norm never exceeds `(h^2/6)*(2/(3*sqrt(3)))*(|M_k| + |M_(k+1)|)`,
            # and a perturbation `E` of a matrix of norm at most one moves
            # `I - S*S'` by at most `2*|E| + |E|^2`; the chord's own deviation
            # peaks at the midpoint, so knots, midpoints and this bound
            # together bound the spline over the whole table.
            h = p.frequencies[k+1] - p.frequencies[k]
            e = 2/(3*sqrt(3))/6*h^2*(opnorm(Matrix(view(p.curvatures,:,:,k))) +
                opnorm(Matrix(view(p.curvatures,:,:,k+1))))
            worst = max(worst, unitaritydeviation(mid) + 2e + e^2)
        end
    end
    return worst
end

"""
    evaluatescattering!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringParameters, ws::AbstractVector,
        absbuffer = nothing)

Evaluate the scattering parameters of `block` at the signed angular
frequencies `ws`, applying the block's negative frequency rule, and write
the matrix at `ws[i]` into `dest[:,:,i]`. With [`ConjugateSymmetry`](@ref)
the provider is evaluated at `abs.(ws)` and conjugated where `ws[i] < 0`;
`absbuffer`, a `Vector{Float64}` the caller may pass, holds the absolute
frequencies and avoids one allocation per call.
"""
function evaluatescattering!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringParameters, ws::AbstractVector,
        absbuffer::Union{Nothing,Vector{Float64}} = nothing)
    if block.negative_frequency isa Native
        evaluateprovider!(dest, block.provider, ws)
    else
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

"""
    evaluatecovariance!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringParameters, ws::AbstractVector,
        absbuffer = nothing)

Evaluate the stated noise covariance of `block`, whose noise model is a
[`NoiseCovariance`](@ref), at the signed angular frequencies `ws` with the
block's negative frequency rule, as [`evaluatescattering!`](@ref) does
the scattering parameters: the covariance of a wave at a negative
frequency is the conjugate of that at the positive one.
"""
function evaluatecovariance!(dest::AbstractArray{Complex{Float64},3},
        block::ScatteringParameters, ws::AbstractVector,
        absbuffer::Union{Nothing,Vector{Float64}} = nothing)
    provider = block.noise.provider
    if block.negative_frequency isa Native
        evaluateprovider!(dest, provider, ws)
    else
        absws = if isnothing(absbuffer)
            abs.(ws)
        else
            resize!(absbuffer, length(ws))
            @inbounds for i in eachindex(ws)
                absbuffer[i] = abs(ws[i])
            end
            absbuffer
        end
        evaluateprovider!(dest, provider, absws)
        for i in eachindex(ws)
            if ws[i] < 0
                dv = view(dest,:,:,i)
                dv .= conj.(dv)
            end
        end
    end
    return dest
end

# A block loaded from a Touchstone file. The reference impedance is read
# from the file's option line; an explicit `zref` which disagrees with it
# is an error, since it is ambiguous between correcting a mislabeled file
# and asking for renormalization.
function touchstonescatteringblock(path::AbstractString; nports, zref,
        grounded, noise, negative_frequency, interpolation, extrapolation,
        atol, dcmodel::AbstractDCModel = ScatteringLimit())
    ts = Touchstone.touchstone_load(path)
    filezref = collect(Float64, ts.reference)
    # any explicit `zref` is checked against the file, an explicit 50 Ohms
    # included; only an omitted one defers to the file
    if !isnothing(zref)
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
    noise = preparenoise(noise, provider, n)
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
port [`ScatteringParameters`](@ref) referenced to `Z0`. Its
[`TransmissionLineProvider`](@ref) is exact at every signed frequency, so
the block uses [`Native`](@ref) negative frequency evaluation. `grounded`
and `noise` are as for `ScatteringParameters`.

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
        preparenoise(noise, provider, 2), Native())
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
        grounded = true, interpolation = :cubic, extrapolation = :error,
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

Each mode is a two terminal port addressed like a port of a
[`ScatteringParameters`](@ref), with the same `grounded` behavior.
`displacement` is the mean displacement `d_0` and is stored but unused. A
`GaussianChannel` is accepted by the circuit representation, but the
harmonic balance solvers do not support it yet and [`compile`](@ref)
throws a [`ComponentNotSupportedError`](@ref) for it.

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
        grounded::Bool = true, interpolation::Symbol = :cubic,
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

# the worst complete positivity margin over the frequencies where both X
# and Y are stored data, validating each point; NaN when either is a
# callable and nothing can be checked
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
        return NaN
    end
end

# check the symmetry of Y and complete positivity at one frequency point
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

# === a pumped scattering block ===

"""
    LinearizedScattering(linearized, wp; ports = nothing, zref = 50.0,
        grounded = true, noise = Lossless(), phase = 0.0,
        interpolation = :cubic, atol = 1e-6, dcmodel = ScatteringLimit(),
        envelope = nothing)
    LinearizedScattering(H, wp; harmonics, nports, zref = 50.0,
        grounded = true, noise = Lossless(), phase = 0.0,
        dcmodel = ScatteringLimit(), envelope = nothing)

The linearized scattering of a pumped device, a linear time-periodic
multiport: a parametric amplifier, converter or isolator in its
periodic steady state, whose small signal response converts between
frequencies separated by harmonics of its pump `wp` (radians per
second). It is described by
harmonic transfer functions `H_k(nu)`, the wave leaving at the absolute
frequency `nu + k*wp` per unit wave incident at `nu`, for the harmonics
`k >= 0` it converts by; a real device has `H_{-k}(nu) = conj(H_k(-nu))`,
which supplies the rest. `H_0` is an ordinary scattering matrix, and
every `H_k` is a function of the signed frequency, so the block is
evaluated natively at negative frequencies.

The harmonic transfer functions act on power waves, as the hybrid stamp
does; the scattering matrix the solvers report is in waves of photons
per second, so a conversion entry of the two differs by the square root
of the ratio of the frequencies, which the first form applies.

The first form builds one from the `linearized` output of
[`hbsolve`](@ref) of the device, with its keyed scattering matrix over
the modes of one pump and its frequencies: the entry from input mode
`n` to output mode `m` at the signal frequency `w` is a sample of
`H_{m-n}` at `w + n*wp`, and the samples of every harmonic from every
mode pair are collected, folded onto `k >= 0`, and tabulated band by
band over the shifted bands with `interpolation`, zero between and
beyond them, since the data says nothing there (see
[`PiecewiseTabulatedProvider`](@ref)). Samples which fall on the same
frequency from different mode pairs must agree to `atol`, relative to
the largest entry, which is what makes the data that of one periodic
steady state, and the block built from the tables is checked against
its declaration over the modes of the solve at every frequency, as
every solve checks it, since the tables hold at one frequency samples
from solves of neighboring signal frequencies whose mode truncations
differ; a device whose mode truncation was too tight fails both, and
`atol` admits the discrepancy of one which was nearly so.
`ports` selects and orders the device's ports which become the block's,
by default all of them, and `zref` gives their reference impedances.
`phase` rotates `H_k` by `exp(im*k*phase)`, the phase of the block's
pump relative to the one the data was computed with, which matters when
other elements of the circuit share the pump.

The second form takes the harmonic transfer functions directly: `H` is a
vector of providers, one per entry of `harmonics` (nonnegative, ascending,
beginning with zero), each a matrix, a callable of the signed angular
frequency, or a tuple `(frequencies, values)` tabulated over signed
frequencies.

The block is stamped by the harmonic balance solvers as a coupling
between the modes of the circuit whose frequencies differ by its
harmonics, so the circuit must be solved with the block's pump: its
`wp` must be a harmonic combination of the circuit's pump frequencies
and the mode set must reach the harmonics the block converts by, which
is how a circuit with no junctions is solved with a pumped block, by
giving `hbsolve` the pump frequency and no source at it. The block is
linear, so it acts on the circuit's own pump harmonics the same way.

`noise = Lossless()` asserts that the device is lossless: its
multi-mode scattering matrix `S` satisfies `J - S J S' = 0` with `J` the
signs of the mode frequencies, which is checked on the block as built
and again over the modes of every solve which evaluates it, in harmonic
balance at every signal frequency and in time at every bath frequency,
and the block then emits no noise. A device with loss states the noise it adds
with `noise = NoiseCovariance(linearized.Cnoise)`, the covariance its
solve reports with `returnCnoise = true`, over the same modes, ports
and frequencies: in the units of `Cnoise`, where a vacuum channel counts
as one, it becomes the harmonic covariances
`V_k(nu) = <n(nu + k wp) n(nu)'>`, sampled like the transfer functions,
and is held to the minimum the commutation relations require,
`V - K` and `V + K` positive semidefinite with `K = J - S J S'`, on the
data and again at every frequency of a sweep. The block then carries
channels of both kinds over all its modes at once, as an active block
does (see [`NoiseCovariance`](@ref)), so it adds the noise its solve
found, correlated across the modes, and its output obeys the
commutation relations. Given by its harmonic transfer functions, the
stated noise is one covariance provider per harmonic. `atol` is the
tolerance of the block's data, relative to the square of the largest
entry of the multi-mode scattering matrix: the consistency of the data
of one periodic steady state, the losslessness declared or the
commutation relations a stated covariance must satisfy, which the
block is checked for on its stored data when it is built, a table at
its knots and a block of constants which does not convert at any one
frequency, and over the modes of every solve, whatever outputs the
solve is asked for, a callable and constants which convert being
checkable only there, and an entry of the data below it is zero where
the pump solve asks whether the block converts the conjugate of a
mode; a covariance's own `atol` counts for its checks as well. A fit of the
block is held to no tolerance: it states the noise its own commutator
requires, a covariance completed to the commutation relations (see
[`NoiseCovariance`](@ref)). `grounded` and
`dcmodel` are as for [`ScatteringParameters`](@ref); the direct current
behavior is that of `H_0` at zero unless stated.

In time the block is realized by [`RationalScattering`](@ref)`(block,
npoles)`, which fits every harmonic transfer function to stable filters:
`H_0` as an ordinary rational block, and each `H_k` as the pair of real
filters of its cosine and sine parts, strictly proper, whose outputs the
transient multiplies by `2 cos(k wp t)` and `-2 sin(k wp t)`. `envelope`
is a callable of the time in seconds which multiplies every harmonic but
the zeroth there, a prescribed gate of the conversion and not a model
of the pump being switched: `H_0`, the filters and a stated covariance
stay those of the pumped device while the conversion is scaled, so the
block satisfies the commutation relations only where the envelope is
one, or zero for a device whose `H_0` is that of the device unpumped.
What the gate is for is the start of a record: ramped from zero, it
leaves the circuit time invariant before the record a noise calculation
needs, as the pumps of the junction circuits are, and the noise is read
once the conversion has been on longer than the block's memory.
`nothing` is a conversion always on, which a noise calculation then
takes as periodic from the start of its record with the fluctuations
before it those of the unconverted response, not of the periodic
device.
"""
struct LinearizedScattering{N,DM<:AbstractDCModel} <: AbstractComponent
    harmonics::Vector{Int}
    providers::Vector{AbstractMatrixProvider}
    wp::Float64
    phase::Float64
    nports::Int
    zref::Vector{Float64}
    grounded::Bool
    noise::N
    dcmodel::DM
    # the gate of the conversion in time, a callable of the time in
    # seconds multiplying every harmonic but the zeroth, or nothing for
    # a block whose conversion is always on
    envelope::Any
    # the tolerance of the block's data, relative to the square of the
    # largest entry of its multi-mode scattering matrix: the consistency
    # of the data, the declaration it must meet, with a covariance's own
    # tolerance, when it is built and at every frequency a solve
    # evaluates it at, and the entries taken as zero where the pump solve
    # asks whether the block converts the conjugate of a mode
    atol::Float64
end

function LinearizedScattering(H::AbstractVector, wp::Real; harmonics::AbstractVector{<:Integer},
        nports::Integer, zref = 50.0, grounded::Bool = true, noise = Lossless(),
        phase::Real = 0.0, dcmodel::AbstractDCModel = ScatteringLimit(),
        envelope = nothing, atol::Real = 1e-6)
    (isfinite(wp) && wp > 0) || throw(ArgumentError("the pump frequency must be positive and finite."))
    (isfinite(atol) && atol >= 0) || throw(ArgumentError("atol must be finite and nonnegative."))
    ks = collect(Int, harmonics)
    (!isempty(ks) && ks[1] == 0 && issorted(ks; lt = <=) && all(>=(0), ks)) || throw(ArgumentError(
        "the harmonics must be nonnegative, strictly increasing, and begin with zero."))
    length(H) == length(ks) || throw(DimensionMismatch(lazy"give one provider per harmonic; $(length(H)) providers for $(length(ks)) harmonics."))
    n = Int(nports)
    providers = AbstractMatrixProvider[]
    for h in H
        p = matrixprovider(h, Complex{Float64}; n = n)
        providersize(p) == n || throw(DimensionMismatch(lazy"a harmonic transfer function has dimension $(providersize(p)) but the block has $(n) ports."))
        push!(providers, p)
    end
    if noise isa NoiseCovariance
        V = noise.provider
        (V isa AbstractVector && length(V) == length(ks)) || throw(ArgumentError(
            "the stated noise of a pumped block given by its harmonic transfer functions is one covariance provider per harmonic, V_k(nu) = <n(nu + k wp) n(nu)'> in the units of Cnoise."))
        vps = AbstractMatrixProvider[]
        for v in V
            vp = matrixprovider(v, Complex{Float64}; n = n)
            providersize(vp) == n || throw(DimensionMismatch(lazy"a harmonic covariance has dimension $(providersize(vp)) but the block has $(n) ports."))
            push!(vps, vp)
        end
        noise = NoiseCovariance(vps, noise.interpolation, noise.extrapolation, noise.atol, noise.completed, noise.padding)
    else
        noise isa Lossless || throw(ArgumentError("a pumped block declares noise = Lossless(), or states its noise with a NoiseCovariance over its harmonics."))
    end
    isfinite(phase) || throw(ArgumentError("the pump phase must be finite."))
    checkdcmodel(dcmodel, n, noise, 1e-8)
    built = LinearizedScattering(ks, providers, Float64(wp), Float64(phase), n,
        zrefvector(zref, n), grounded, noise, dcmodel, envelope, Float64(atol))
    # the stored data is checked where it is stored: a table at its
    # knots, over the modes the harmonics reach from each with every
    # input which feeds them, and a block of constants which does not
    # convert at any one frequency, its data being the same at every;
    # a callable, and constants which convert, whose entries carry the
    # photon factors of the modes, are checked at the solve
    conjugateladder(built)
    stored = vcat(providers, noise isa NoiseCovariance ? noise.provider : AbstractMatrixProvider[])
    nus = Float64[]
    for p in stored
        append!(nus, tableknots(p))
    end
    isempty(nus) && ks == [0] && all(p -> p isa ConstantMatrixProvider, stored) && push!(nus, Float64(wp))
    declared = noise isa NoiseCovariance ? max(Float64(atol), noise.atol) : Float64(atol)
    for nu in unique!(sort!(nus))
        rows, cols, K = pumpedfamily(built, (nu,))
        v = pumpedviolation(built, rows, cols, K)
        v <= declared || throw(ArgumentError(lazy"the block's data does not meet what it declares: over the modes its harmonics reach from $(nu) rad/s the violation of its losslessness or of the commutation relations of its stated covariance is $(v) of the square of its largest entry, against the $(declared) of atol and its noise model's. A device with loss or gain states its noise with NoiseCovariance; raise atol to admit a discrepancy of the data."))
    end
    return built
end

function LinearizedScattering(lin, wp::Real; ports = nothing, zref = 50.0,
        grounded::Bool = true, noise = Lossless(), phase::Real = 0.0,
        interpolation::Symbol = :cubic, atol::Real = 1e-6,
        dcmodel::AbstractDCModel = ScatteringLimit(), envelope = nothing)
    (hasproperty(lin, :S) && hasproperty(lin, :w)) || throw(ArgumentError(
        "give the linearized output of hbsolve, which carries the scattering matrix and its frequencies, or the harmonic transfer functions with their harmonics."))
    S = lin.S
    S isa AxisKeys.KeyedArray || throw(ArgumentError(
        "the scattering matrix must be keyed by mode and port; solve with keyedarrays = true."))
    (isfinite(wp) && wp > 0) || throw(ArgumentError("the pump frequency must be positive and finite."))
    (isfinite(atol) && atol >= 0) || throw(ArgumentError("atol must be finite and nonnegative."))
    noise isa Lossless || noise isa NoiseCovariance || throw(ArgumentError(
        "a pumped block declares noise = Lossless(), or states the covariance its solve reports with noise = NoiseCovariance(linearized.Cnoise)."))
    outmodes = collect(AxisKeys.axiskeys(S, 1))
    outports = collect(AxisKeys.axiskeys(S, 2))
    inmodes = collect(AxisKeys.axiskeys(S, 3))
    inports = collect(AxisKeys.axiskeys(S, 4))
    all(m -> length(m) == 1, outmodes) && all(m -> length(m) == 1, inmodes) || throw(ArgumentError(
        "a pumped block has one pump; the data has the modes of several."))
    outmodes == inmodes || throw(ArgumentError("the output and input modes of the data differ."))
    modes = Int[m[1] for m in outmodes]
    selected = isnothing(ports) ? outports : collect(ports)
    all(p -> p in outports && p in inports, selected) || throw(ArgumentError(
        lazy"the ports $(selected) are not all ports of the data, $(outports)."))
    allunique(selected) || throw(ArgumentError("a port is selected twice."))
    n = length(selected)
    n >= 1 || throw(ArgumentError("select at least one port."))
    po = [findfirst(==(p), outports) for p in selected]
    qi = [findfirst(==(p), inports) for p in selected]
    A = Array(S)
    wv = Float64.(collect(lin.w))
    length(wv) == size(A, 5) || throw(DimensionMismatch("the frequencies do not match the scattering data."))
    nm = length(modes)
    wp = Float64(wp)
    # the samples of each harmonic from every mode pair at every
    # frequency, folded onto nonnegative harmonics by the realness of the
    # device
    samples = Dict{Int,Vector{Tuple{Float64,Matrix{Complex{Float64}}}}}()
    mirrored = Tuple{Float64,Matrix{Complex{Float64}}}[]
    for (a, mo) in enumerate(modes), (b, mi) in enumerate(modes), i in eachindex(wv)
        k = mo - mi
        nu = wv[i] + mi*wp
        # the solver's waves are in photons per second, the stamp's in
        # power: a conversion between frequencies carries the ratio of
        # the photon energies
        power = sqrt(abs(nu + k*wp)/abs(nu))
        M = Matrix{Complex{Float64}}(undef, n, n)
        for q in 1:n, p in 1:n
            M[p, q] = power*A[a, po[p], b, qi[q], i]
        end
        if k >= 0
            push!(get!(samples, k, Tuple{Float64,Matrix{Complex{Float64}}}[]), (nu, M))
        else
            push!(get!(samples, -k, Tuple{Float64,Matrix{Complex{Float64}}}[]), (-nu, conj(M)))
        end
        # the unconverted response is a real function, `H_0(-nu) =
        # conj(H_0(nu))`, so its table holds both signs of every sample,
        # and is evaluated at a mode the data reached only at the other
        # sign, an idler's; a mirrored sample is added only where the data
        # has none of its own
        k == 0 && push!(mirrored, (-nu, conj(M)))
    end
    if haskey(samples, 0)
        direct = sort!([nu for (nu, M) in samples[0]])
        for (nu, M) in mirrored
            i = searchsortedfirst(direct, nu)
            near = (i <= length(direct) && abs(direct[i] - nu) <= 1e-9*(abs(nu) + wp)) ||
                (i > 1 && abs(direct[i - 1] - nu) <= 1e-9*(abs(nu) + wp))
            near || push!(samples[0], (nu, M))
        end
    end
    scale = max(1.0, maximum(abs, A))
    ks = sort!(collect(keys(samples)))
    providers = AbstractMatrixProvider[]
    for k in ks
        list = sort!(samples[k]; by = first)
        nus = Float64[]
        mats = Matrix{Complex{Float64}}[]
        for (nu, M) in list
            if !isempty(nus) && abs(nu - nus[end]) <= 1e-9*(abs(nu) + wp)
                # the same frequency reached from two mode pairs: one
                # periodic steady state gives the same value from both
                d = maximum(abs, M .- mats[end])
                d <= atol*scale || throw(ArgumentError(lazy"the data is not that of one periodic steady state: the samples of the harmonic $(k) at $(nu) rad/s from two mode pairs differ by $(d) against the largest entry $(scale). A tighter mode truncation, or data not from the linearized solve of one pumped device, does this; atol admits a discrepancy."))
                continue
            end
            push!(nus, nu)
            push!(mats, M)
        end
        values = Array{Complex{Float64},3}(undef, n, n, length(nus))
        for j in eachindex(nus)
            values[:, :, j] .= mats[j]
        end
        push!(providers, piecewisetable(nus, values; interpolation = interpolation))
    end
    # the stated noise: the covariance the solve reports, over the same
    # modes, ports and frequencies, as harmonic covariances
    # V_k(nu) = <n(nu + k wp) n(nu)'> for k >= 0, in the units of Cnoise,
    # the entry from mode n to mode m at the signal frequency w being a
    # sample of V_{m-n} at w + n*wp and the Hermitian symmetry supplying
    # the negative harmonics
    stated = if noise isa NoiseCovariance
        C = noise.provider
        C isa AxisKeys.KeyedArray || throw(ArgumentError(
            "the stated noise of a pumped block is the keyed covariance its solve reports, noise = NoiseCovariance(linearized.Cnoise), from hbsolve with returnCnoise = true."))
        (collect(AxisKeys.axiskeys(C, 1)) == outmodes && collect(AxisKeys.axiskeys(C, 2)) == outports &&
            collect(AxisKeys.axiskeys(C, 3)) == inmodes && collect(AxisKeys.axiskeys(C, 4)) == inports &&
            size(C, 5) == length(wv) && collect(AxisKeys.axiskeys(C, 5)) == collect(AxisKeys.axiskeys(S, 5))) || throw(ArgumentError(
            "the stated covariance does not share the modes, ports and frequencies of the scattering matrix; take both from one solve."))
        Ca = Array(C)
        vsamples = Dict{Int,Vector{Tuple{Float64,Matrix{Complex{Float64}}}}}()
        for (a, mo) in enumerate(modes), (b, mi) in enumerate(modes), i in eachindex(wv)
            k = mo - mi
            nu = wv[i] + mi*wp
            M = Matrix{Complex{Float64}}(undef, n, n)
            for q in 1:n, p in 1:n
                M[p, q] = Ca[a, po[p], b, qi[q], i]
            end
            if k >= 0
                push!(get!(vsamples, k, Tuple{Float64,Matrix{Complex{Float64}}}[]), (nu, M))
            else
                push!(get!(vsamples, -k, Tuple{Float64,Matrix{Complex{Float64}}}[]), (nu + k*wp, Matrix(M')))
            end
        end
        vscale = max(1.0, maximum(abs, Ca))
        vproviders = AbstractMatrixProvider[]
        for k in ks
            list = sort!(get(vsamples, k, Tuple{Float64,Matrix{Complex{Float64}}}[]); by = first)
            nus = Float64[]
            mats = Matrix{Complex{Float64}}[]
            for (nu, M) in list
                if !isempty(nus) && abs(nu - nus[end]) <= 1e-9*(abs(nu) + wp)
                    d = maximum(abs, M .- mats[end])
                    d <= noise.atol*vscale || throw(ArgumentError(lazy"the stated covariance is not Hermitian over the modes: the samples of the harmonic $(k) at $(nu) rad/s from two mode pairs differ by $(d) against the largest entry $(vscale)."))
                    continue
                end
                push!(nus, nu)
                push!(mats, M)
            end
            values = Array{Complex{Float64},3}(undef, n, n, length(nus))
            for j in eachindex(nus)
                values[:, :, j] .= mats[j]
            end
            push!(vproviders, piecewisetable(nus, values; interpolation = interpolation))
        end
        NoiseCovariance(vproviders, noise.interpolation, noise.extrapolation, noise.atol, noise.completed, noise.padding)
    else
        noise
    end
    isfinite(phase) || throw(ArgumentError("the pump phase must be finite."))
    checkdcmodel(dcmodel, n, stated, 1e-8)
    built = LinearizedScattering(ks, providers, wp, Float64(phase), n,
        zrefvector(zref, n), grounded, stated, dcmodel, envelope, Float64(atol))
    # the block as every solve evaluates it, from its tables, over the
    # modes of the solve at each of its frequencies: a lossless device
    # has a symplectic multi-mode scattering matrix there, and a stated
    # covariance is held to the minimum the commutation relations
    # require of that matrix; the tables hold, at one frequency, samples
    # from solves of neighboring signal frequencies, whose mode
    # truncations differ, so the block meets its declaration only to
    # the discrepancy of the data, which atol admits and which is found
    # here, when the block is built, rather than at its first solve
    conjugateladder(built)
    declared = stated isa NoiseCovariance ? max(Float64(atol), stated.atol) : Float64(atol)
    sq = Vector{Float64}(undef, nm)
    for i in eachindex(wv)
        sq .= wv[i] .+ modes .* wp
        v = pumpedviolation(built, sq, pumpedharmonics(built, sq))
        v <= declared || throw(ArgumentError(lazy"the block does not meet what it declares: at the frequency index $(i), over the modes of the solve, the violation of its losslessness or of the commutation relations of its stated covariance is $(v) of the square of its largest entry, against the $(declared) of atol and its noise model's. Declare its noise with noise = NoiseCovariance(linearized.Cnoise) from a solve with returnCnoise = true, leave a pump port out with the ports keyword, or raise atol to admit the discrepancy between the solves the tables hold at one frequency."))
    end
    return built
end

"""
    evaluateharmoniccovariances!(dest::AbstractArray{Complex{Float64},4},
        block::LinearizedScattering, ws::AbstractVector)

Evaluate the harmonic covariances of the stated noise of `block`, whose
noise model is a [`NoiseCovariance`](@ref) over its harmonics, at the
signed angular frequencies `ws`: `dest[:, :, j, i]` is
`V_k(ws[i]) = <n(ws[i] + k wp) n(ws[i])'>` for `k = block.harmonics[j]`,
in the units of `Cnoise`, where a vacuum channel counts as one, rotated
by the block's pump phase as the transfer functions are. The negative
harmonics follow from `V_{-k}(nu) = V_k(nu - k wp)'`, and the conjugate
ladder of a frequency from `V_k(-nu - k wp) = transpose(V_k(nu))` (see
[`conjugateladder`](@ref)), which is how the covariance a solve reports,
whose rows are the modes of the solve, states the noise at their
conjugates.

This is the covariance every solve reads, zero at a frequency neither
the data nor its conjugate ladder reaches.
"""
function evaluateharmoniccovariances!(dest::AbstractArray{Complex{Float64},4},
        block::LinearizedScattering, ws::AbstractVector)
    n = block.nports
    nk = length(block.harmonics)
    size(dest) == (n, n, nk, length(ws)) || throw(DimensionMismatch(lazy"the destination has size $(size(dest)) but ($(n), $(n), $(nk), $(length(ws))) is required."))
    return evaluatecoveredharmonics!(dest, block, ws; covariance = true)
end

"""
    evaluateharmonics!(dest::AbstractArray{Complex{Float64},4},
        block::LinearizedScattering, ws::AbstractVector)

Evaluate the harmonic transfer functions of `block` at the signed
angular frequencies `ws`: `dest[:, :, j, i]` is `H_k(ws[i])` for the
harmonic `k = block.harmonics[j]`, rotated by the block's pump phase.
The negative harmonics follow from `H_{-k}(nu) = conj(H_k(-nu))`, which a
caller evaluates at the negated frequencies.
"""
function evaluateharmonics!(dest::AbstractArray{Complex{Float64},4},
        block::LinearizedScattering, ws::AbstractVector)
    n = block.nports
    nk = length(block.harmonics)
    size(dest) == (n, n, nk, length(ws)) || throw(DimensionMismatch(lazy"the destination has size $(size(dest)) but ($(n), $(n), $(nk), $(length(ws))) is required."))
    buf = Array{Complex{Float64},3}(undef, n, n, length(ws))
    for (j, k) in enumerate(block.harmonics)
        evaluateprovider!(buf, block.providers[j], ws)
        rot = cis(k*block.phase)
        @inbounds for i in eachindex(ws), q in 1:n, p in 1:n
            dest[p, q, j, i] = rot*buf[p, q, i]
        end
    end
    return dest
end

# the scattering matrix a pumped block has without conversion, `H_0`,
# which is what its zero frequency rows and its share of an ordinary
# stamp read; natively valid at signed frequencies
function evaluatescattering!(dest::AbstractArray{Complex{Float64},3},
        block::LinearizedScattering, ws::AbstractVector,
        absbuffer::Union{Nothing,Vector{Float64}} = nothing)
    return evaluateprovider!(dest, block.providers[1], ws)
end

"""
    pumpedharmonics(block::LinearizedScattering, offsets::AbstractVector)

The harmonic of `block` which couples each pair of modes of a circuit, as
a matrix over (output mode, input mode) of the multiple of the block's
pump by which the modes' frequency `offsets` differ, or `typemin(Int)`
where they differ by no such multiple or by one the block does not
convert by. Throws if the block converts and no pair of modes is a pump
apart, which is a circuit solved without the block's pump.
"""
function pumpedharmonics(block::LinearizedScattering, offsets::AbstractVector)
    nm = length(offsets)
    K = fill(typemin(Int), nm, nm)
    found = false
    for m in 1:nm, n in 1:nm
        d = Float64(offsets[m]) - Float64(offsets[n])
        k = round(Int, d/block.wp)
        abs(d - k*block.wp) <= 1e-6*block.wp || continue
        abs(k) in block.harmonics || continue
        K[m, n] = k
        k == 0 || (found = true)
    end
    # a block with no positive harmonic does not convert and needs no
    # pair of modes a harmonic apart; its unconverted response is stamped
    # on the diagonal like any other block's
    (found || !any(>(0), block.harmonics)) || throw(ArgumentError(lazy"a pumped block converts by harmonics of its pump at $(block.wp) rad/s, but no two modes of this solve differ by one of them: solve the circuit with the block's pump frequency and modes which reach its harmonics, as hbsolve does when given the pump frequency."))
    return K
end

"""
    ModulatedRationalProvider(cosine::RationalScatteringProvider,
        sine::RationalScatteringProvider)

The harmonic transfer function of a [`LinearizedScattering`](@ref) block for a
harmonic `k > 0` as fitted for the transient: `H_k(nu) = G_c(i nu) +
i G_s(i nu)` with `G_c` and `G_s` real rational functions, the cosine and
sine parts `G_c(nu) = (H_k(nu) + conj(H_k(-nu)))/2` and
`G_s(nu) = (H_k(nu) - conj(H_k(-nu)))/(2i)`, each realized as a stable
state space with no feedthrough, since a conversion vanishes at infinite
frequency. Natively valid at signed frequencies. In time the block
multiplies the output of the cosine filter by `2 cos(k wp t)` and that of
the sine filter by `-2 sin(k wp t)`, which is `2 Re[exp(i k wp t) (h_k * a)]`,
the sum of the harmonics `k` and `-k` of a real system.
"""
struct ModulatedRationalProvider <: AbstractMatrixProvider
    cosine::RationalScatteringProvider
    sine::RationalScatteringProvider
end
providersize(p::ModulatedRationalProvider) = size(p.cosine.D, 1)
function evaluateprovider!(dest::AbstractArray{T,3},
        p::ModulatedRationalProvider, ws::AbstractVector) where T
    n = providersize(p)
    checkdestsize(dest, n, length(ws))
    buf = Array{Complex{Float64},3}(undef, n, n, length(ws))
    evaluateprovider!(dest, p.cosine, ws)
    evaluateprovider!(buf, p.sine, ws)
    dest .+= im .* buf
    return dest
end

# whether every harmonic of a pumped block has a realization in time
realizedintime(block::LinearizedScattering) = all(p -> p isa RationalScatteringProvider ||
    p isa ModulatedRationalProvider, block.providers)
