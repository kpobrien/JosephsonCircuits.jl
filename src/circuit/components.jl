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

`temperature` is the physical temperature in Kelvin, finite and
nonnegative, which sets the noise a lossy instance adds; a lossless one adds
none. `nothing`, the default, takes the temperature the analysis is run at.

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
Inductor(L; temperature = nothing) = Inductor(L, checktemperature(temperature, "an Inductor"))

"""
    Capacitor(C; temperature = nothing)

A two terminal linear capacitor with capacitance `C` in Farads. Terminals are
`1` and `2`. The value may be a number or a symbolic variable.

`temperature` is the physical temperature in Kelvin, finite and
nonnegative, which sets the noise a lossy instance adds; a lossless one adds
none. `nothing`, the default, takes the temperature the analysis is run at.

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
Capacitor(C; temperature = nothing) = Capacitor(C, checktemperature(temperature, "a Capacitor"))

"""
    Resistor(R; temperature = nothing)

A two terminal linear resistor with resistance `R` in Ohms. Terminals are `1`
and `2`. The value may be a number or a symbolic variable.

`temperature` is the physical temperature in Kelvin, finite and
nonnegative, which sets the noise a dissipative instance adds. `nothing`, the
default, takes the temperature the analysis is run at.

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
Resistor(R; temperature = nothing) = Resistor(R, checktemperature(temperature, "a Resistor"))

# A temperature a component, a noise model or an analysis states:
# `nothing`, for the one the analysis is run at, or a finite nonnegative
# number of Kelvin, since a negative one would give a mode a negative
# occupation. `what` says whose it is, for the message.
checktemperature(::Nothing, what::AbstractString) = nothing
function checktemperature(t::Real, what::AbstractString)
    (isfinite(t) && t >= 0) || throw(ArgumentError(lazy"The temperature of $(what) is $(t) K; a temperature must be finite and nonnegative."))
    return Float64(t)
end

# === sources and analysis ports ===

"""
    CurrentSource(I)

A two terminal source of the constant current `I` in Amperes, which it
drives through itself from its first terminal to its second: it draws the
current from the node at its first terminal and delivers it to the node at
its second, the opposite sense of a port source, which injects its current
into the port's first (positive) terminal. The value is a number, or a
parameter or an expression in parameters whose number `circuitdefs`
supplies.

Harmonic balance drives the zero frequency mode with the constant current,
so a solve of a circuit holding a nonzero one retains that mode,
`dc = true`. A transient solve takes the constant as well, unless a
[`TransientSource`](@ref) names the source, whose waveform then replaces
it.
"""
struct CurrentSource{T} <: AbstractComponent
    I::T
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
    MatchedTermination(; temperature = 0.0)

A source and load environment matched to the port's reference impedance,
acting across the two port terminals, at the physical `temperature` in
kelvin, finite and nonnegative. This is the default: a port owns its
environment, so no resistor should be added in order to terminate it.

The environment sends into the port the thermal field of a matched load:
the symmetrized noise `nbar + 1/2` in each mode, `nbar` its
[`thermaloccupation`](@ref) at the mode's frequency, which is the vacuum's
half photon at zero temperature. The temperature is the port's own and is
zero unless it is given here, whatever `temperature` an analysis gives its
dissipative elements: the termination is the measurement line, not the
device. For a line whose noise is not that of one temperature, give the
temperature of the field arriving at the port
([`effectivetemperature`](@ref)), or put the line's attenuators in the
circuit at their own temperatures.

# Examples
```jldoctest
julia> MatchedTermination(temperature = 0.05).temperature
0.05
```
"""
struct MatchedTermination <: AbstractPortTermination
    temperature::Float64
    MatchedTermination(; temperature = 0.0) =
        new(checktemperature(temperature, "a port's matched termination"))
end

# the physical temperature of a port's environment, zero for a port which
# owns none
porttemperature(t::MatchedTermination) = t.temperature
porttemperature(::AbstractPortTermination) = 0.0

"""
    NoPortTermination()

No port owned environment, written as `termination = nothing`. The port
remains an excitation and observation boundary with its reference impedance
intact, but contributes no physical loading of its own.
"""
struct NoPortTermination <: AbstractPortTermination end

# normalize the `termination` keyword of `Port`: `nothing` means no
# termination
porttermination(t::AbstractPortTermination) = t
porttermination(::Nothing) = NoPortTermination()
porttermination(x) = throw(ArgumentError(lazy"The port termination $(x) is not recognized. Write termination = nothing for an unterminated boundary, or omit the keyword for the default matched environment."))

# the resistor of the circuit a termination names as the port's own
# environment, by its instance identifier, or `nothing` for a termination
# which names none; only the deprecated tuple netlist's names one (see
# circuit/legacy.jl)
namedtermination(::AbstractPortTermination) = nothing

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
        termination = MatchedTermination())

An analysis port with the port number `number` and reference impedance `Z0`
in Ohms. A `Port` identifies an electrical port for excitation and
observation; excitation amplitudes belong to the analysis arguments, not to
the circuit topology.

By default the port owns a matched external source and load environment of
impedance `Z0` acting across its two terminals, so a port needs no resistor
to define its impedance. The environment is at zero temperature unless
`termination = MatchedTermination(temperature = T)` warms it (see
[`MatchedTermination`](@ref)); its thermal field then enters the port as
the field of a matched load at `T`. The environment acts between the port terminals and
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
    showtermination(io, p.termination)
    print(io, ")")
end
showtermination(io::IO, t::MatchedTermination) = iszero(t.temperature) ? nothing :
    print(io, ", termination = MatchedTermination(temperature = ", t.temperature, ")")
showtermination(io::IO, ::NoPortTermination) = print(io, ", termination = nothing")

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

The sign of `K` follows the order each inductor's terminals are declared
in: with `K > 0` currents entering the first terminal of each inductor add
flux to both, as with the dots of a SPICE coupling statement at the first
named nodes, and reversing the terminals of one inductor, or the sign of
`K`, opposes them.

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
    PolynomialCPR(coefficients; atol = 1e-12)

A current-phase relation (CPR) specified by the coefficients of its
polynomial expansion `f(φ) = coefficients[1]*φ + coefficients[2]*φ^2 + ...`,
where `φ` is the reduced branch phase. The linear coefficient must equal one,
to `atol`, so that the `L0` of the containing [`NonlinearInductor`](@ref) is
the small signal inductance. The object is callable and its analytic
derivative is available through [`cprderivative`](@ref).

This supports specifying the effective nonlinearity of a SNAIL, SQUID,
Quarton, or kinetic inductor directly through its expansion coefficients
without wiring up the underlying junction arrangement. An array of `N`
identical junctions in series, for example, divides the phase and so has
the relation `N*sin(φ/N)` with small signal inductance `N*Lj`, whose
expansion is `[1, 0, -1/(6N^2), 0, 1/(120N^4), ...]`. A kinetic inductor
whose inductance rises with its current as `L(I) = L0*(1 + I^2/Istar^2)`
carries the flux `L0*(I + I^3/(3*Istar^2))`, which inverted for the
current is the relation `[1, 0, -(IL/Istar)^2/3, 0, (IL/Istar)^4/3, ...]`,
with `IL = phi0/L0` the current scale of the element; the quintic term
matters only as the current approaches `Istar`.

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

function PolynomialCPR(coefficients::AbstractVector{T}; atol::Real = 1e-12) where T
    if isempty(coefficients)
        throw(ArgumentError("PolynomialCPR requires at least the linear coefficient."))
    end
    c1 = coefficients[1]
    if c1 isa Real && !isapprox(c1, one(c1); atol = atol)
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

# two relations with the same coefficients are one relation, so that
# nonlinear inductors written from equal expansions compare and hash as
# equal; the coefficients compare as `isequal` does, which the hash of
# the coefficient vector agrees with
Base.:(==)(p::PolynomialCPR, q::PolynomialCPR) = isequal(p.a, q.a)
Base.hash(p::PolynomialCPR, h::UInt) = hash(p.a, hash(:PolynomialCPR, h))

"""
    PolynomialCPRDerivative(a)

The analytic derivative of a [`PolynomialCPR`](@ref), produced by
[`cprderivative`](@ref). Callable.
"""
struct PolynomialCPRDerivative{T}
    a::Vector{T}
end
(p::PolynomialCPRDerivative)(φ) = evalpoly(φ, p.a)
Base.:(==)(p::PolynomialCPRDerivative, q::PolynomialCPRDerivative) =
    isequal(p.a, q.a)
Base.hash(p::PolynomialCPRDerivative, h::UInt) = hash(p.a, hash(:PolynomialCPRDerivative, h))

"""
    cprderivative(cpr)

The analytic derivative of a [`PolynomialCPR`](@ref), or of one of its
derivatives, as a callable [`PolynomialCPRDerivative`](@ref). The solvers
differentiate the other relation they evaluate, the sinusoidal Josephson
one, as `cos` and `-sin` themselves.
"""
function cprderivative end

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

"""
    JunctionRelations(value, derivative, negsecond, third,
        polynomialjunctions, sinusoidaljunctions)

The current-phase relations of the Josephson junction branches of a
circuit, whose junction axis is the order of the nonzero entries of the
branch inductance vector `Ljb`, the order of the junction axis of every
time domain array the solvers hold.

A relation is either the sinusoidal Josephson one or a
[`PolynomialCPR`](@ref), and the junctions are grouped by kind once:
`polynomialjunctions` and `sinusoidaljunctions` are the ascending
positions on the junction axis of each, so that each relation is
evaluated over its own junctions alone. The polynomial coefficients of
the polynomial junctions sit in the rows of `value`, in the order of
`polynomialjunctions`, in `evalpoly` order along the second axis and
padded with zeros to one common degree, so that a single Horner loop
evaluates them all; `derivative` and `negsecond` hold the coefficients of
the first derivative and of the *negative* of the second, which is the
combination the Hessian and the derivative of the linearized system with
respect to the operating point are written in, and `third` those of the
third derivative, which the trilinear form of the problem interface
takes.

A circuit whose junctions are all sinusoidal, which is every circuit that
does not ask for anything else, has the empty table (see
[`allsinusoidal`](@ref)), and the solvers take the plain `sin` and `cos`
of the Josephson relation.
"""
struct JunctionRelations{M,V}
    value::M
    derivative::M
    negsecond::M
    third::M
    polynomialjunctions::V
    sinusoidaljunctions::V
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
    polynomial = findall(!isnothing, cprs)
    # one common degree, so that the Horner loop is the same length for
    # every junction
    nterms = maximum(j -> length(cprs[j].a), polynomial)
    value = zeros(Float64, length(polynomial), nterms)
    derivative = zeros(Float64, length(polynomial), nterms)
    negsecond = zeros(Float64, length(polynomial), nterms)
    third = zeros(Float64, length(polynomial), nterms)
    for (row, j) in enumerate(polynomial)
        c = cprs[j]
        d = cprderivative(c)
        d2 = cprderivative(d)
        d3 = cprderivative(d2)
        copycoefficients!(value, row, c.a)
        copycoefficients!(derivative, row, d.a)
        copycoefficients!(negsecond, row, -1 .* d2.a)
        copycoefficients!(third, row, d3.a)
    end
    return JunctionRelations(value, derivative, negsecond, third, polynomial,
        findall(isnothing, cprs))
end

"""
    sinusoidalmask(r::JunctionRelations)

Whether each junction of `r`, along the junction axis, has the sinusoidal
Josephson relation, on the host.
"""
function sinusoidalmask(r::JunctionRelations)
    mask = fill(false, length(r.polynomialjunctions) + length(r.sinusoidaljunctions))
    mask[Array(r.sinusoidaljunctions)] .= true
    return mask
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
circuit whose junctions are all sinusoidal takes a single broadcast of
`sin`.
"""
function relationinto!(out, r::JunctionRelations, phi)
    allsinusoidal(r) && return out .= sin.(phi)
    return applyrelationfirst!(out, phi, r.value, r, sin)
end

"""
    derivativeinto!(out, r::JunctionRelations, phi)

[`derivativeat`](@ref) writing into `out`, which may not alias `phi`. A
circuit whose junctions are all sinusoidal takes a single broadcast of
`cos`.
"""
function derivativeinto!(out, r::JunctionRelations, phi)
    allsinusoidal(r) && return out .= cos.(phi)
    return applyrelationfirst!(out, phi, r.derivative, r, cos)
end

"""
    negsecondat(r::JunctionRelations, phi)

The negative of the second derivative of the relation of every junction
at the branch phases `phi`, `sin.(phi)` for the Josephson relation: the
change of the differential inductance with the phase, which the
derivative of a linearization with respect to the phases carries. The
counterpart of [`derivativeat`](@ref).
"""
negsecondat(r::JunctionRelations, phi) = negsecondinto!(similar(phi), r, phi)

function negsecondinto!(out, r::JunctionRelations, phi)
    allsinusoidal(r) && return out .= sin.(phi)
    return applyrelationfirst!(out, phi, r.negsecond, r, sin)
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
    Array(r.polynomialjunctions), Array(r.sinusoidaljunctions))

function hostrelations(r::JunctionRelations, rows::AbstractVector{Int})
    allsinusoidal(r) && return hostrelations(r)
    # the row of the tables each junction's coefficients are in, zero for
    # a sinusoidal junction
    polynomial = Array(r.polynomialjunctions)
    tablerow = zeros(Int, length(polynomial) + length(r.sinusoidaljunctions))
    tablerow[polynomial] .= eachindex(polynomial)
    kept = [tablerow[j] for j in rows if tablerow[j] > 0]
    sub(m) = Array(m)[kept, :]
    return JunctionRelations(sub(r.value), sub(r.derivative), sub(r.negsecond),
        sub(r.third), findall(j -> tablerow[j] > 0, rows), findall(j -> tablerow[j] == 0, rows))
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
    similar(A, Int, 0), similar(A, Int, 0))

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
    index(v) = copyto!(similar(A, Int, length(v)), v)
    return JunctionRelations(move(r.value), move(r.derivative),
        move(r.negsecond), move(r.third), index(r.polynomialjunctions),
        index(r.sinusoidaljunctions))
end

function copycoefficients!(m::Matrix, j::Integer, a::AbstractVector)
    for k in eachindex(a)
        m[j, k] = a[k]
    end
    return m
end

# `out .= f.(src)` where `f` is the relation of each junction of `r`: a
# polynomial, the row of `c` for the junction, in `evalpoly` order, and
# `trig` for a sinusoidal junction, each evaluated over its own junctions
# alone. The junction is the last axis of `src`, as it is in every time
# domain array of harmonic balance, so a coefficient column is reshaped to
# broadcast along it. Everything here is a broadcast of whole arrays or of
# their views onto the junctions of one kind, which is what keeps it
# device generic: no kernel, no scalar indexing, and one pass per degree
# rather than one per junction. The relation `trig` is a type parameter,
# since a function only passed on to a broadcast is not specialized on
# and would be dispatched at run time.
function applyrelationlast!(out::AbstractArray, src::AbstractArray,
        c::AbstractMatrix, r::JunctionRelations, trig::F) where F
    s = r.sinusoidaljunctions
    isempty(s) && return hornerlast!(out, src, c)
    N, p = ndims(out), r.polynomialjunctions
    hornerlast!(selectdim(out, N, p), selectdim(src, N, p), c)
    selectdim(out, N, s) .= trig.(selectdim(src, N, s))
    return out
end
function hornerlast!(out::AbstractArray, src::AbstractArray, c::AbstractMatrix)
    shape = ntuple(d -> d == ndims(out) ? size(c, 1) : 1, ndims(out))
    column(k) = reshape(view(c, :, k), shape)
    nterms = size(c, 2)
    out .= column(nterms)
    for k in nterms-1:-1:1
        out .= out .* src .+ column(k)
    end
    return out
end

# the same for the transient, whose branch phases carry the junction on
# the first axis, `(junction,)` or `(junction, condition)`, against which a
# coefficient column broadcasts as it is
function applyrelationfirst!(out::AbstractArray, src::AbstractArray,
        c::AbstractMatrix, r::JunctionRelations, trig::F) where F
    s = r.sinusoidaljunctions
    isempty(s) && return hornerfirst!(out, src, c)
    p = r.polynomialjunctions
    hornerfirst!(selectdim(out, 1, p), selectdim(src, 1, p), c)
    selectdim(out, 1, s) .= trig.(selectdim(src, 1, s))
    return out
end
function hornerfirst!(out::AbstractArray, src::AbstractArray, c::AbstractMatrix)
    nterms = size(c, 2)
    out .= view(c, :, nterms)
    for k in nterms-1:-1:1
        out .= out .* src .+ view(c, :, k)
    end
    return out
end

"""
    NonlinearInductor(L0, cpr)

A two terminal nonlinear inductive element defined by its current-phase
relation: `I(φ) = (phi0/L0)*cpr(φ)` where `φ` is the reduced branch phase,
`L0` is the small signal inductance in Henries, and `cpr` is a callable with
unit slope at zero. The solvers take its derivatives analytically, those
of a [`PolynomialCPR`](@ref) through [`cprderivative`](@ref).

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
struct NonlinearInductor{T,F} <: AbstractComponent
    L0::T
    cpr::F
end

"""
    JosephsonJunction(Lj)
    JosephsonJunction(; Ic)

A Josephson junction with junction inductance `Lj` in Henries, or
equivalently critical current `Ic` in Amperes, and the sinusoidal
current-phase relation `I(φ) = Ic*sin(φ)`. Equal to
`NonlinearInductor(Lj, sin)`.

# Examples
```jldoctest
julia> JosephsonJunction(100e-12) == NonlinearInductor(100e-12, sin)
true

julia> JosephsonJunction(Ic = 1e-6).L0
3.291059784754534e-10
```
"""
JosephsonJunction(Lj) = NonlinearInductor(Lj, sin)
JosephsonJunction(; Ic) = NonlinearInductor(phi0/Ic, sin)

# equal by their inductance and their relation, and hashed alike, so
# that `isequal`, sets and dictionaries of components agree with `==`
Base.:(==)(a::NonlinearInductor, b::NonlinearInductor) =
    isequal(a.L0, b.L0) && a.cpr == b.cpr
Base.hash(c::NonlinearInductor, h::UInt) =
    hash(c.cpr, hash(c.L0, hash(:NonlinearInductor, h)))

"""
    issinusoidal(c::NonlinearInductor)

Whether the current-phase relation of `c` is the sinusoidal Josephson
relation `sin`, which the solvers evaluate as `sin` and `cos` (see
[`junctioncpr`](@ref)). Every `NonlinearInductor` compiles to the `:Lj`
type, and one with another relation has it recorded beside the table (see
[`CompiledCircuit`](@ref)).
"""
issinusoidal(c::NonlinearInductor) = c.cpr === sin

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
[`ScatteringParameters`](@ref). A provider implements
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
write a block by hand and the default) or `:inplace` (`f(dest, w)` writes
into one, which matters when the same provider is evaluated many thousands
of times); see [`CALLABLE_FORMS`](@ref).
"""
struct CallableMatrixProvider{F} <: AbstractMatrixProvider
    f::F
    n::Int
    # how `f` is called. `:matrix` returns a fresh n by n matrix, which is
    # the natural way to write a block by hand and is the default. `:inplace`
    # writes into one, `f(dest, w)`, which matters when the same provider is
    # evaluated many thousands of times and the returned matrices dominate
    # the allocation of a sweep.
    form::Symbol
end

CallableMatrixProvider(f, n::Int) = CallableMatrixProvider(f, n, :matrix)

"""
    CALLABLE_FORMS

The ways a callable provider may be called. See
[`CallableMatrixProvider`](@ref).
"""
const CALLABLE_FORMS = (:matrix, :inplace)

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
    if !all(isfinite, frequencies)
        throw(ArgumentError("Tabulated frequencies must be finite."))
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

# A frequency within roundoff of an end knot of a table is on the knot:
# the edge knots of a table assembled from several sources can differ
# from the frequencies a solver forms by an ulp, and so can a frequency
# carried to the conjugate ladder and back, or padded by the pump and
# back. The evaluation and the coverage test admit the same roundoff,
# so what a solver takes a table to hold it can evaluate.
edgetolerance(f1::Real, fn::Real) = 8eps(Float64)*max(abs(f1), abs(fn))
edgetolerance(f::AbstractVector) = edgetolerance(f[1], f[end])

# whether a table holds the frequency `nu`, to the roundoff its
# evaluation admits at the edges
function tablecovers(t::TabulatedMatrixProvider, nu::Real)
    f = t.frequencies
    edgetol = edgetolerance(f)
    return f[1] - edgetol <= nu <= f[end] + edgetol
end

function evaluateprovider!(dest::AbstractArray{T,3},
        p::TabulatedMatrixProvider, ws::AbstractVector) where T
    checkdestsize(dest, providersize(p), length(ws))
    if p.extrapolation == :error && !all(w -> tablecovers(p, w), ws)
        f = p.frequencies
        outofrange = [w for w in ws if !tablecovers(p, w)]
        throw(ArgumentError(lazy"The angular frequencies $(outofrange) are outside the tabulated range [$(f[1]), $(f[end])] rad/s. Extrapolation of tabulated data is opt-in: pass extrapolation = :constant, :linear or :zero if extrapolation is intended, or fit the block with RationalScattering, which extrapolates as a passive rational function."))
    end
    for i in eachindex(ws)
        tablevalue!(view(dest, :, :, i), p, ws[i])
    end
    return dest
end

# the table at the frequency `w` into the matrix `d`: the knot, the
# interpolant between knots, or beyond them the extrapolation the table
# declares
@inline function tablevalue!(d::AbstractMatrix{T}, p::TabulatedMatrixProvider, w) where T
    f = p.frequencies
    p.extrapolation == :zero && !tablecovers(p, w) && return fill!(d, zero(T))
    # a frequency beyond an end knot by that roundoff alone is
    # evaluated at the knot, whatever the block extrapolates by
    edgetol = edgetolerance(f)
    if f[1] - edgetol <= w < f[1]
        w = f[1]
    elseif f[end] < w <= f[end] + edgetol
        w = f[end]
    end
    if w <= f[1]
        if p.extrapolation == :constant || length(f) == 1 || w == f[1]
            d .= view(p.values,:,:,1)
        elseif p.interpolation == :linear
            # the first segment continues
            lerpslices!(d, p, 1, 2, w)
        else
            edgeslices!(d, p, 1, w)
        end
    elseif w >= f[end]
        if p.extrapolation == :constant || length(f) == 1 || w == f[end]
            d .= view(p.values,:,:,length(f))
        elseif p.interpolation == :linear
            lerpslices!(d, p, length(f)-1, length(f), w)
        else
            edgeslices!(d, p, 2, w)
        end
    else
        j = searchsortedlast(f, w)
        if f[j] == w
            d .= view(p.values,:,:,j)
        elseif p.interpolation == :linear
            lerpslices!(d, p, j, j+1, w)
        else
            splineslices!(d, p, j, w)
        end
    end
    return d
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
    for i in eachindex(ws)
        w = ws[i]
        d = view(dest, :, :, i)
        j = findfirst(t -> tablecovers(t, w), p.tables)
        isnothing(j) ? fill!(d, zero(T)) : tablevalue!(d, p.tables[j], w)
    end
    return dest
end

# Tabulated data in disjoint bands, one table per band, each with the
# interpolation given and zero beyond it. The ascending knots `nus` are
# samples taken `step` apart at most within a band, so consecutive knots
# further apart than that, beyond the roundoff of shifting them, which
# `shifttol` of their size and the step bounds, are in different bands,
# and the data is interpolated between its samples and never across a
# gap; a step of zero, the samples of one frequency, puts each knot in a
# band of its own.
function piecewisetable(nus::Vector{Float64}, values::Array{T,3}; interpolation::Symbol = :cubic,
        step::Real, shifttol::Real = 1e-9) where T
    n = length(nus)
    tables = TabulatedMatrixProvider{T}[]
    start = 1
    for k in 1:n
        if k == n || nus[k + 1] - nus[k] > step + shifttol*(abs(nus[k + 1]) + step)
            push!(tables, TabulatedMatrixProvider(nus[start:k], values[:, :, start:k];
                interpolation = interpolation, extrapolation = :zero))
            start = k + 1
        end
    end
    return PiecewiseTabulatedProvider(tables)
end

"""
    RotatedMatrixProvider(provider, delays; offset = 0.0, phase = 1)

The covariance of the noise a block emits, stated by `provider`, moved
to other reference planes: the entry from `q` to `p` is multiplied by
`phase` and by `cis((w + offset)*delays[p] - w*delays[q])`, which takes
a lossless line of `delays[p]` seconds off each port, the wave `p` emits
carrying `cis(-w*delay)` along it and the wave it is correlated with the
conjugate. `offset` is how far above the frequency asked for the
emitting side is read: zero for the covariance of an ordinary block, and
`k*wp` for the harmonic covariance `V_k(nu) = <n(nu + k wp) n(nu)'>` of a
pumped block, whose correlated waves are the harmonic `k` of its pump
apart. `phase` has unit modulus, a pump phase `cis(k*phi)` for such a
harmonic.

The two sides carry opposite phases, as a correlation's do, where the
scattering matrix of the same block carries the sum of its delays, so
this is the transformation of a covariance, or of a pumped harmonic
covariance, and a block refuses it as scattering data. The phase is
applied to what `provider` returns at each requested frequency rather
than to data it stores, so it composes with the interpolation and the
extrapolation that provider declares instead of changing what they mean,
and a provider of any kind can carry it. This is the change of reference
plane a fit makes to a stated covariance when a delay is taken out of the
data (see [`RationalScattering`](@ref)).
"""
struct RotatedMatrixProvider{P<:AbstractMatrixProvider} <: AbstractMatrixProvider
    provider::P
    delays::Vector{Float64}
    offset::Float64
    phase::Complex{Float64}
end

function RotatedMatrixProvider(provider::AbstractMatrixProvider, delays::AbstractVector;
        offset::Real = 0.0, phase::Number = 1)
    taus = collect(Float64, delays)
    length(taus) == providersize(provider) || throw(DimensionMismatch(
        lazy"give one delay per port ($(providersize(provider))); got $(length(taus))."))
    all(isfinite, taus) || throw(ArgumentError("the delays must be finite."))
    (isfinite(offset) && isfinite(phase)) || throw(ArgumentError("the offset and the phase must be finite."))
    # a change of reference plane turns the covariance and scales nothing
    abs(abs(phase) - 1) <= 8eps(Float64) || throw(ArgumentError(lazy"the phase of a rotation has unit modulus, since moving the reference planes of a covariance scales nothing; got $(phase)."))
    return RotatedMatrixProvider(provider, taus, Float64(offset), Complex{Float64}(phase))
end

providersize(p::RotatedMatrixProvider) = providersize(p.provider)

# The provider a rotation turns, through any number of them. A rotation
# holds data wherever the provider it turns does, so this says whether a
# check can run on stored data and at which samples; the check reads
# what the rotation returns there.
unrotated(p::RotatedMatrixProvider) = unrotated(p.provider)
unrotated(p) = p

# === which data a provider holds, and where ===

# whether a provider holds data at the frequency `nu`: a table within
# its knots, or everywhere when it declares how it extrapolates, a
# piecewise table within one of its bands, the data being the samples
# and what lies between them, with no declaration made beyond them, a
# rotation wherever the provider it turns does, and a provider of any
# other kind, a callable, a constant or a filter, everywhere
providercovers(p::TabulatedMatrixProvider, nu::Real) = p.extrapolation != :error || holdsdata(p, nu)
providercovers(p::PiecewiseTabulatedProvider, nu::Real) = holdsdata(p, nu)
providercovers(p::RotatedMatrixProvider, nu::Real) = providercovers(p.provider, nu)
providercovers(p, nu::Real) = true

# whether a provider holds a sample of its own at `nu`, the knots of a
# table reaching it before any extrapolation, which is where its data can
# be checked against a relation it must obey, and the same samples
# through a rotation; a provider which is not tabulated states its value
# everywhere
holdsdata(p::TabulatedMatrixProvider, nu::Real) = tablecovers(p, nu)
holdsdata(p::PiecewiseTabulatedProvider, nu::Real) = any(t -> tablecovers(t, nu), p.tables)
holdsdata(p::RotatedMatrixProvider, nu::Real) = holdsdata(p.provider, nu)
holdsdata(p, nu::Real) = true

# the knots of a tabulated provider, ascending, of a piecewise one those
# of every band, of a rotation the knots of what it turns, and none for a
# provider of any other kind
tableknots(p::TabulatedMatrixProvider) = p.frequencies
tableknots(p::PiecewiseTabulatedProvider) = sort!(reduce(vcat, (t.frequencies for t in p.tables)))
tableknots(p::RotatedMatrixProvider) = tableknots(p.provider)
tableknots(p) = Float64[]

function evaluateprovider!(dest::AbstractArray{T,3},
        p::RotatedMatrixProvider, ws::AbstractVector) where T
    n = providersize(p)
    checkdestsize(dest, n, length(ws))
    evaluateprovider!(dest, p.provider, ws)
    # the rotation is the outer product of a phase per emitting port
    # with one per receiving port, so it is two phases a port at a
    # frequency and not one an entry
    u = Vector{Complex{Float64}}(undef, n)
    v = Vector{Complex{Float64}}(undef, n)
    @inbounds for i in eachindex(ws)
        w = ws[i]
        for r in 1:n
            u[r] = p.phase*cis((w + p.offset)*p.delays[r])
            v[r] = cis(-w*p.delays[r])
        end
        for q in 1:n, r in 1:n
            dest[r, q, i] *= u[r]*v[q]
        end
    end
    return dest
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

# The scattering data of a block as a provider. A rotation carries opposite
# phases on its two sides, the transformation of a covariance (see
# RotatedMatrixProvider), where the scattering matrix of a block moved to
# other reference planes carries the sum of the delays on both, so it is
# refused here rather than read as a different block.
function scatteringprovider(S; kwargs...)
    p = matrixprovider(S, Complex{Float64}; kwargs...)
    p isa RotatedMatrixProvider && throw(ArgumentError("a RotatedMatrixProvider moves the reference planes of a covariance and is not scattering data: take a delay out of scattering data before it is stored, as RationalScattering does with its delays keyword, or state it as a TransmissionLine in cascade."))
    return p
end

function matrixprovider(f, ::Type{T}; n = nothing, form::Symbol = :matrix,
        kwargs...) where T
    if !(form in CALLABLE_FORMS)
        throw(ArgumentError(lazy"Unknown form $(repr(form)). Supported: :matrix (f(w) returns an n by n matrix) and :inplace (f(dest, w) writes one)."))
    end
    if isnothing(n)
        throw(ArgumentError("The matrix dimension cannot be inferred from a callable provider; pass the dimension explicitly (nports for a ScatteringParameters)."))
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
block adds noise of symmetrized covariance `(nbar + 1/2)(I - S S')` (see
[`ScatteringNoisePlan`](@ref)), with `nbar` the occupation at the
`temperature` given to the analysis, which defaults to zero temperature and
so to the vacuum's `(I - S S')/2`.

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
its data can be shown to be unitary (see [`provablylossless`](@ref)),
which is possible for a constant matrix and a table but not for a
callable, whose values away from the evaluated frequencies are unknown.
A lossless callable would therefore carry noise channels which are
identically zero, at a real cost on a circuit with many blocks;
`Lossless()` removes them.

The assertion is held to the block's `atol` (see
[`checkblockcontract`](@ref)): when a constant, tabulated or rational
block is built, at every sample of stored data and between and beyond
them, and over every frequency of a realization, and it is an error if
false. A callable is checked by a sweep at every frequency it evaluates
the block at and is otherwise taken on trust, and asserting it of a
block which does absorb omits the noise the block should add, making
the quantum efficiency and the commutation relations wrong by that
much. A pumped block is checked when it is built and again over the
modes of every solve which evaluates it, to the block's `atol`; a fit of
a pumped block is not lossless to better than its error, and states the
noise it needs instead (see [`NoiseCovariance`](@ref)).
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
    # the one constructor, so that a matrix of any element type, the
    # `Matrix{Float64}` it is stored as included, is checked
    function ScatteringDC(S0::AbstractMatrix)
        size(S0, 1) == size(S0, 2) ||
            throw(DimensionMismatch(lazy"a zero frequency scattering matrix must be square; got $(size(S0))."))
        all(isfinite, S0) ||
            throw(ArgumentError("a zero frequency scattering matrix must be finite."))
        all(iszero∘imag, S0) ||
            throw(ArgumentError("a zero frequency scattering matrix must be real: at zero frequency there is no phase to carry an imaginary part."))
        return new(Matrix{Float64}(real.(S0)))
    end
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
at the physical temperature `temperature` in Kelvin, finite and
nonnegative. The added noise
covariance is `(nbar + 1/2)(I - S S')`, the commutator `I - S S'` weighted
by the symmetrized noise of a mode at that temperature, `nbar` its
[`thermaloccupation`](@ref).

This states the block's temperature where the block is defined, so it
overrides the `temperature` argument of the analysis. At zero temperature,
under an analysis at zero temperature, it coincides with
[`Passive`](@ref).
"""
struct ThermalEquilibrium{T}
    temperature::T
    function ThermalEquilibrium(temperature::Real)
        checktemperature(temperature, "a ThermalEquilibrium noise model")
        return new{typeof(temperature)}(temperature)
    end
end

"""
    NoiseCovariance(V; interpolation = :cubic, extrapolation = :error,
        atol = 1e-8, completed = false, padding = 4)

The noise a [`ScatteringParameters`](@ref) adds, stated outright rather
than derived from its loss, which is how an active block, an amplifier
given by its scattering parameters, declares its noise. `V` is the
symmetrized covariance of the noise wave the block emits at its ports, in
quanta, the units of the rest of the noise outputs, where a vacuum channel
counts as `1/2` and a channel at temperature `T` as `nbar + 1/2`: the same
units as `Cnoise`. It may be a matrix, a callable
of angular frequency, or a tuple `(frequencies, values)` of tabulated
data, following the same provider forms as the scattering data, and is
evaluated with the block's negative frequency rule.

The block's output must obey the commutation relations, so with
`K = I - S S'` the added noise has the commutator `K`, and a covariance is
realizable only when `V - K/2` and `V + K/2` are both positive
semidefinite. `V + K/2` then factors into channels which emit like a mode
in its vacuum and `V - K/2` into channels which emit like the conjugate of
one, the idler channels of an amplifier; each channel carries the vacuum's
half photon, so the block adds half their sum, `V`. A passive block in
thermal equilibrium is the case
`V = (nbar + 1/2) K`. A phase insensitive amplifier of power gain `G` from
port 1 to port 2 has `K[2,2] = 1 - G`, so `V[2,2] >= (G - 1)/2`, the noise
of a quantum limited amplifier, and one with an input referred added noise
of `nadd` photons has `V[2,2] = G*nadd` (see [`noisequanta`](@ref) for a
noise temperature); `V[1,1]` is what it emits backward out of its input,
`(nbar + 1/2)*K[1,1]` for an input matched at its physical temperature.

The condition is checked on the eigenvalues to `atol` when the block is
built, at the samples where both the scattering data and `V` are stored,
and at every frequency a solver evaluates the block at. `V` is also held
to Hermitian symmetry, to `atol` of its largest entry, or of one for a
covariance whose entries are below one (see
[`checkblockcontract`](@ref)). A block with this model carries no
temperature: its noise is `V`, whatever the analysis temperature.

With `completed = true` the covariance is completed to the commutation
relations rather than held to them: `V` is replaced by
`V + neg(V - K/2) + neg(V + K/2)`, `neg` taking the negative part of a
Hermitian matrix, the sum of `-lambda v v'` over its negative
eigenvalues, which makes `V - K/2` and `V + K/2` positive semidefinite and
adds nothing where they are. For `V = 0` the addition is `|K|/2`, the
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
commutator requires, as its stamp takes it. This is how a fit states
its noise, of a pumped block and of an ordinary one which states a
covariance (see [`RationalScattering`](@ref)): a fit is neither lossless
nor consistent with a stated covariance to better than its error, and
the completion turns that error into noise the block emits, where a
tolerance would only excuse it.
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
    checkatol(atol)
    padding >= 0 || throw(ArgumentError("padding must be nonnegative."))
    return NoiseCovariance(V, interpolation, extrapolation, Float64(atol), completed, Int(padding))
end

# a tolerance of a block's or a covariance's data
checkatol(atol) = (isfinite(atol) && atol >= 0) || throw(ArgumentError("atol must be finite and nonnegative."))

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
- a path to a Touchstone file of S, Y or Z parameters, from which the
  reference impedances are also read; Y and Z parameters are converted to
  S at those impedances.

Tabulated data is interpolated with the cubic spline through each entry's
samples (`interpolation = :cubic`, the default; `:linear` takes the chords
between them, which lag a rotating phase) and is never extrapolated unless
asked: `extrapolation` is `:error` by default, with `:constant`,
`:linear` and `:zero` to opt in (see [`TabulatedMatrixProvider`](@ref)).
The harmonic balance solvers evaluate a block wherever their mixing
products fall, which can be far outside the band the data covers, so measured data meant for them is better fitted with
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
minimum the commutation relations require of it. `Lossless` says the
block is unitary, which spares a callable, whose data cannot be shown to
be, the noise channels of a loss it does not have.
`negative_frequency` is [`ConjugateSymmetry`](@ref) (default) or
[`Native`](@ref). `form` says how a callable `S` is called; see
[`CallableMatrixProvider`](@ref).

`atol` is the tolerance of the block's data, which the block keeps: it
is passive, the smallest eigenvalue of `I - S S'` no lower than `-atol`,
unless it states its noise with a `NoiseCovariance`, which permits an
active block, and declared `Lossless` it is unitary, no entry of
`I - S S'` larger than `atol` (see [`checkblockcontract`](@ref)). The
data is held to this when the block is built, wherever it is stored: at
the samples of constant and tabulated data, over every frequency of a
[`RationalScattering`](@ref) realization, and between and beyond the
samples of a table declared lossless. A callable is known only where it
is evaluated.

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
as `S` (a matrix, a callable of angular frequency, or tabulated data),
tabulated data interpolated and extrapolated as `S` is and a callable
called in the `form` `S` is. The block holds them as `name => provider`
pairs in the order given, which are no part of its type, so blocks which
differ in their derivatives alone are one type, for which the solvers are
compiled once. A block depends on a design parameter through these
entries alone. A
derivative is not a scattering matrix and is never passivity checked. The
block's scattering matrix and its derivatives describe one design point:
the definitions move the parameters a component value is written in, not
values a block's data or closure has captured, so a block whose data
depends on a parameter must be restated at each point along with its
derivatives.
A Touchstone path is tabulated data: its derivatives are read as a
table's are, and `form`, which says how a callable is called, is refused
with it as with any other data.

# Examples
```jldoctest
julia> ScatteringParameters([0 1;1 0]).nports
2
```
"""
struct ScatteringParameters{P,N,NF,DM<:AbstractDCModel} <: AbstractComponent
    provider::P
    nports::Int
    zref::Vector{Float64}
    grounded::Bool
    noise::N
    negative_frequency::NF
    # the providers of the analytic dS/dtheta by design parameter name, in
    # the order given, for [`designsensitivities`](@ref); empty for a
    # block which depends on no design parameter. The providers are held
    # as `Any`, so that neither the names nor the providers' types are
    # part of the block's type: blocks which differ in their derivatives
    # alone share one, and what is compiled for a block is compiled once
    # for them. A derivative is read into a block of its own (see
    # [`designblockjacobian`](@ref)), whose type is concrete.
    derivatives::Vector{Pair{Symbol,Any}}
    # the zero frequency behavior, when the block's own data does not give
    # it; see [`AbstractDCModel`](@ref)
    dcmodel::DM
    # the tolerance of the block's data: how far from passive, or from
    # unitary where it declares itself lossless, it may be wherever it is
    # read; see [`checkblockcontract`](@ref)
    atol::Float64
end

# positional constructors without derivatives, without a stated zero
# frequency model and at the default tolerance, for a block built without
# the checks: a derivative, or a block of zeros
ScatteringParameters(provider, nports::Int, zref::Vector{Float64},
    grounded::Bool, noise, negative_frequency) =
    ScatteringParameters(provider, nports, zref, grounded, noise,
        negative_frequency, Pair{Symbol,Any}[], ScatteringLimit(), 1e-8)

function ScatteringParameters(S; nports = nothing, zref = nothing,
        grounded::Bool = true, noise = Passive(),
        negative_frequency = ConjugateSymmetry(),
        interpolation::Symbol = :cubic, extrapolation::Symbol = :error,
        form::Symbol = :matrix,
        derivatives::NamedTuple = NamedTuple(),
        dcmodel::AbstractDCModel = ScatteringLimit(),
        atol::Real = 1e-8)
    # the derivatives' names are no part of the block's type, nor of what
    # is compiled to build it
    @nospecialize derivatives
    # a Touchstone file states its reference impedances; omitted otherwise,
    # the reference impedance is 50 Ohms at every port
    provider, zrefs = if S isa AbstractString
        touchstoneprovider(S, zref; interpolation, extrapolation, form)
    else
        (scatteringprovider(S; n = nports, interpolation, extrapolation, form),
            something(zref, 50.0))
    end
    n = providersize(provider)
    if !isnothing(nports) && n != nports
        throw(DimensionMismatch(lazy"nports = $(nports) does not match the scattering data dimension $(n)."))
    end
    return checkedblock(provider, n, zrefvector(zrefs, n), grounded, noise,
        negative_frequency, derivativeproviders(derivatives, n; interpolation, extrapolation, form),
        dcmodel, atol)
end

# The providers of a block's derivatives, of its dimension `n`, by
# parameter name in the order given: one given as data is interpolated and
# extrapolated as the block's data is, and one given as a callable is
# called in the `form` the block's callable is. The named tuple is read by
# its field names, so nothing here is compiled for the names.
function derivativeproviders(derivatives::NamedTuple, n::Int;
        interpolation::Symbol, extrapolation::Symbol, form::Symbol)
    @nospecialize derivatives
    out = Pair{Symbol,Any}[]
    for k in fieldnames(typeof(derivatives))
        v = getfield(derivatives, k)
        dp = matrixprovider(v, Complex{Float64}; n = n,
            interpolation = interpolation, extrapolation = extrapolation,
            form = v isa Union{AbstractMatrix,Tuple,AbstractMatrixProvider} ? :matrix : form)
        providersize(dp) == n || throw(DimensionMismatch(lazy"the derivative for parameter $(k) has dimension $(providersize(dp)) but the block has $(n) ports."))
        push!(out, k => dp)
    end
    return out
end

# A block from its parts with its noise prepared and its data held to its
# contract, whatever built it, from a matrix, a file, a line or a
# realization: at every sample it stores, between and beyond them where
# the data determines it, and in its zero frequency model. The noise
# model and the negative frequency rule are ones the solvers read.
@noinline function checkedblock(provider, n::Int, zref::Vector{Float64}, grounded::Bool,
        noise, negative_frequency, derivatives, dcmodel, atol::Real; normtested::Bool = false)
    checkatol(atol)
    noise isa Union{Passive,Lossless,ThermalEquilibrium,NoiseCovariance} || throw(ArgumentError(
        lazy"the noise model of a scattering block is Passive(), Lossless(), ThermalEquilibrium(T) or NoiseCovariance(V), not $(repr(noise))."))
    negative_frequency isa Union{ConjugateSymmetry,Native} || throw(ArgumentError(
        lazy"the negative frequency rule of a scattering block is ConjugateSymmetry() or Native(), not $(repr(negative_frequency))."))
    block = ScatteringParameters(provider, n, zref, grounded, preparenoise(noise, n),
        negative_frequency, derivatives, dcmodel, Float64(atol))
    checkstoreddata(block)
    checkbeyondsamples(block; normtested)
    checkdcmodel(dcmodel, n, block.noise, block.atol)
    return block
end

# A stated zero frequency matrix is checked like the block's own data: its
# size against the block, and passivity to the block's tolerance, the
# smallest eigenvalue of `I - S S'` no lower than `-atol` (see
# checkblockcontract), unless the block declared itself active with a
# `NoiseCovariance`.
checkdcmodel(::ScatteringLimit, n::Int, noise, atol) = nothing
function checkdcmodel(m::AbstractDCModel, n::Int, noise, atol)
    S0 = dcscatteringmatrix(m, n)
    if !(noise isa NoiseCovariance)
        margin = passivitymargin(S0)
        margin >= -atol || throw(ArgumentError(lazy"the zero frequency scattering matrix is active at direct current: the smallest eigenvalue of I - S S' is $(margin), below -atol for the block's atol of $(atol). An active block has to declare its own noise with NoiseCovariance, as it does at every other frequency."))
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
    preparenoise(noise, n)

The noise model of a block as it is stored: a [`NoiseCovariance`](@ref)
with its data as a matrix provider of the block's dimension `n`, and any
other model unchanged. The data is checked with the block's (see
[`checkblockcontract`](@ref)).
"""
function preparenoise(noise::NoiseCovariance, n::Int)
    vp = matrixprovider(noise.provider, Complex{Float64}; n = n,
        interpolation = noise.interpolation,
        extrapolation = noise.extrapolation)
    if providersize(vp) != n
        throw(DimensionMismatch(lazy"The noise covariance dimension $(providersize(vp)) does not match the number of ports $(n)."))
    end
    return NoiseCovariance(vp, noise.interpolation, noise.extrapolation,
        noise.atol, noise.completed, noise.padding)
end
preparenoise(noise, n::Int) = noise

"""
    commutationmargin(V::AbstractMatrix, K::AbstractMatrix)

The smallest eigenvalue of `V - K/2` and of `V + K/2`, which is
nonnegative when the symmetrized noise covariance `V` meets the commutation
relations of the commutator `K`: `V + K/2` and `V - K/2` then factor into
the channels of either kind (see [`NoiseCovariance`](@ref)).
"""
commutationmargin(V::AbstractMatrix, K::AbstractMatrix) =
    min(minimum(eigvals(Hermitian(V - K/2))), minimum(eigvals(Hermitian(V + K/2))))

"""
    quantumnoisemargin(V::AbstractMatrix, S::AbstractMatrix)

The [`commutationmargin`](@ref) of the noise covariance `V` against
`K = I - S S'`, which is nonnegative when `V` is one a block with the
scattering matrix `S` can add without violating the commutation
relations.
"""
function quantumnoisemargin(V::AbstractMatrix, S::AbstractMatrix)
    n = size(S, 1)
    K = Matrix{Complex{Float64}}(I, n, n) - S*S'
    return commutationmargin(Matrix{Complex{Float64}}(V), K)
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

# The largest singular value at which the dissipation `I - S'S` reaches
# `-atol`, the passivity a block's data is held to (see
# checkblockcontract): the level every test of a realization or a fit
# over frequency holds its largest singular value to.
passivelevel(atol) = sqrt(1 + atol)

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

# The two frequency identities of a pumped block's modes: whether the
# difference `d` of two frequencies is the harmonic `k` of the pump `wp`,
# to `harmonictolerance` of the pump frequency, and whether two
# frequencies on its ladders are one, to `ladderroundoff` of their size
# and the pump frequency's. Every such comparison takes these tolerances,
# the grouping of frequencies into the ladders of a pump (pumpladders)
# among them.
const harmonictolerance = 1e-6
const ladderroundoff = 1e-9
isharmonic(d, k, wp) = abs(d - k*wp) <= harmonictolerance*wp
samefrequency(a, b, wp) = abs(a - b) <= ladderroundoff*(abs(a) + wp)

"""
    checkblockcontract(block::ScatteringParameters, S, V, w; name = nothing)

Hold the data of `block` at the angular frequency `w` to what the block
declares, and throw an `ArgumentError` saying where it does not. `S` is
the block's scattering matrix at `w` and `V` the covariance its
[`NoiseCovariance`](@ref) states there, either `nothing` where it is not
known:

- the data is finite;
- a stated covariance is Hermitian, no entry of `V - V'` larger than the
  covariance's `atol` times its largest entry, or than its `atol` for a
  covariance of entries below one;
- it is at least what the commutation relations require of the block,
  the smallest eigenvalue of `V - K/2` and of `V + K/2` with `K = I - S S'`
  (see [`quantumnoisemargin`](@ref)) no lower than `-atol` with the
  covariance's `atol`, unless it is completed to them;
- a block declared [`Lossless`](@ref) is unitary, no entry of `I - S S'`
  (see [`unitaritydeviation`](@ref)) larger than the block's `atol`;
- any other block which does not state its noise is passive, the
  smallest eigenvalue of `I - S S'` (see [`passivitymargin`](@ref)) no
  lower than `-atol` with the block's `atol`.

`w` is `nothing` for data which is the same at every frequency, and
`name` is the block's instance in a circuit; both are for the message.
The construction of a block holds every sample its data stores to this
contract, and a solver holds the data it evaluates to it, so the same
data is given the same verdict wherever it is read.
"""
Base.@nospecializeinfer function checkblockcontract(@nospecialize(block::ScatteringParameters), S, V, w;
        name = nothing)
    subject = isnothing(name) ? "the scattering block" : "the scattering block at $(name)"
    at = isnothing(w) ? "" : " at $(w) rad/s"
    noise = block.noise
    if !isnothing(V)
        all(isfinite, V) || throw(ArgumentError(lazy"the noise covariance of $(subject) is not finite$(at)."))
        skew = maximum(abs, V .- V'; init = 0.0)
        tol = noise.atol*max(1.0, maximum(abs, V; init = 0.0))
        skew <= tol || throw(ArgumentError(lazy"the noise covariance of $(subject) is not Hermitian$(at): the largest entry of V - V' is $(skew), against $(tol), the covariance's atol of its largest entry."))
    end
    isnothing(S) && return nothing
    all(isfinite, S) || throw(ArgumentError(lazy"the scattering matrix of $(subject) is not finite$(at)."))
    if noise isa NoiseCovariance
        (isnothing(V) || noise.completed) && return nothing
        margin = quantumnoisemargin(V, S)
        margin >= -noise.atol || throw(ArgumentError(lazy"the noise covariance of $(subject) is less than the commutation relations require$(at): the smallest eigenvalue of V - K/2 or V + K/2, with K = I - S S', is $(margin), below the covariance's atol of $(noise.atol). An amplifier of power gain G has to emit at least (G - 1)/2 at its output; see NoiseCovariance."))
    elseif noise isa Lossless
        deviation = unitaritydeviation(S)
        deviation <= block.atol || throw(ArgumentError(lazy"$(subject) declares noise = Lossless(), but$(at) the largest entry of I - S S' is $(deviation), above its atol of $(block.atol). A block which dissipates must carry the noise its loss requires; use the default Passive() noise model."))
    else
        margin = passivitymargin(S)
        margin >= -block.atol || throw(ArgumentError(lazy"$(subject) is not passive$(at): the smallest eigenvalue of I - S S' is $(margin), below its atol of $(block.atol). An active block states its noise with NoiseCovariance."))
    end
    return nothing
end

# whether a provider's values are stored data, which a check can read
# when the block is built: a constant, a table, a piecewise table, or a
# rotation of one; a callable, a realization and a line compute theirs
# where a solver asks for them
isstored(p) = unrotated(p) isa Union{ConstantMatrixProvider,TabulatedMatrixProvider,PiecewiseTabulatedProvider}

# The frequencies at which the stored data of the providers `ps` is read
# together: the knots of every table among them, or, where none is a
# table, zero alone, at which a constant states what it states at every
# frequency. The second value says which, so that a message names a
# frequency only for data which has one.
function samplefrequencies(ps...)
    knots = sort!(unique!(reduce(vcat, (tableknots(p) for p in ps); init = Float64[])))
    isempty(knots) && return [0.0], false
    return knots, true
end

# a provider read at those of the frequencies `fs` at which its data is
# stored and covers them, and which those are: none for a provider whose
# values are computed
function storedsamples(p, fs::AbstractVector)
    inside = [isstored(p) && providercovers(p, w) for w in fs]
    n = isnothing(p) ? 0 : providersize(p)
    values = Array{Complex{Float64},3}(undef, n, n, count(inside))
    any(inside) && evaluateprovider!(values, p, fs[inside])
    return values, inside
end

# The contract of a block at every sample its data stores: every
# frequency at which its scattering data or its stated covariance is
# stored, each read there where it covers the frequency. A provider which
# computes its values is held to the contract where a solver evaluates
# it.
function checkstoreddata(block::ScatteringParameters)
    S = block.provider
    V = block.noise isa NoiseCovariance ? block.noise.provider : nothing
    fs, isdata = samplefrequencies(S, V)
    Ss, inS = storedsamples(S, fs)
    Vs, inV = storedsamples(V, fs)
    # A rotation of a constant covariance turns with the frequency and a
    # constant scattering matrix does not, so the commutation relations
    # between the two hold or fail frequency by frequency and no sample of
    # either decides them: the pair has no stored frequency, and a solver
    # checks it at every frequency it evaluates the block at. The
    # rotation's Hermitian symmetry is the same at every frequency (see
    # RotatedMatrixProvider) and is read at its sample. The scattering
    # matrix of a block which states its noise enters only the relations.
    turned = !isdata && V isa RotatedMatrixProvider
    turned && (inS .= false)
    kS = kV = 0
    for k in eachindex(fs)
        Sk = inS[k] ? view(Ss, :, :, kS += 1) : nothing
        Vk = inV[k] ? view(Vs, :, :, kV += 1) : nothing
        checkblockcontract(block, Sk, Vk, isdata ? fs[k] : nothing)
    end
    return nothing
end

# The resolvent of a realization through the real Schur form of its state
# matrix, `A = Z T Z'` with `Z` orthogonal and `T` quasi-upper-triangular:
# `(i w I - A)^-1 B = Z (i w I - T)^-1 Z' B`. The form is taken once, and a
# frequency then costs a quasi-triangular solve of the states squared
# rather than a factorization of the states cubed, as accurately, the
# Schur vectors being orthogonal whatever the conditioning of the
# eigenvectors, which a nearly defective realization, an all pass, has no
# usable set of. The form stays real, a complex pair of eigenvalues being
# a 2 by 2 block of its diagonal, so the solve multiplies real entries of
# `T` into the complex solution. A state matrix already quasi-upper-
# triangular, as a fitted block's block diagonal one is, is its own form,
# `Z` the identity, and `A`, `B` and `C` are held as they are. `poles[k]`
# is the eigenvalue of `T` at state `k`, a 2 by 2 block's pair laid out on
# its two states, and `first[k]` the first nonzero row of column `k` of
# `T`, where the back substitution's update from that column starts, so
# that a block diagonal `T` is solved in time proportional to the states.
struct ResolventFactors
    T::Matrix{Float64}
    ZtB::Matrix{Float64}
    CZ::Matrix{Float64}
    poles::Vector{ComplexF64}
    first::Vector{Int}
end
function resolventfactors(A, B, C)
    Am, Bm, Cm = convert(Matrix{Float64}, A), convert(Matrix{Float64}, B), convert(Matrix{Float64}, C)
    T, ZtB, CZ = isquasitriangular(Am) ? (Am, Bm, Cm) : (F = schur(Am); (F.T, F.Z'*Bm, Cm*F.Z))
    first = [something(findfirst(!iszero, view(T, 1:k, k)), k) for k in axes(T, 2)]
    return ResolventFactors(T, ZtB, CZ, schureigenvalues(T), first)
end

# whether `A` is a real Schur form as LAPACK standardizes one: zero below
# its first subdiagonal, with no two consecutive entries of that nonzero,
# so that its diagonal blocks are 1 by 1 or 2 by 2, and each 2 by 2 block
# `[α b; c α]` with `b c < 0`, a complex pair `α ± i sqrt(-b c)`. A 2 by 2
# block of any other shape can hold two real eigenvalues, which its
# entries give only through a difference that cancels, so a matrix with
# one is not taken as its own form.
function isquasitriangular(A::AbstractMatrix)
    n = size(A, 1)
    for j in 1:n, i in j + 2:n
        iszero(A[i, j]) || return false
    end
    for k in 2:n - 1
        (iszero(A[k, k - 1]) || iszero(A[k + 1, k])) || return false
    end
    for k in 2:n
        iszero(A[k, k - 1]) && continue
        (A[k - 1, k - 1] == A[k, k] && signbit(A[k - 1, k]) != signbit(A[k, k - 1]) &&
            !iszero(A[k - 1, k])) || return false
    end
    return true
end

# `X <- ((s I - T)^-1 X')'` for the factors' quasi-upper-triangular `T`
# and a complex shift `s`, `X` holding a row per port and a column per
# state: back substitution a block of the diagonal at a time, a 1 by 1
# block a real eigenvalue and a 2 by 2 block, marked by its nonzero
# subdiagonal entry, a complex pair, each solved by the inverse of the
# shifted block for every port at once and then taken out of the columns
# of the states above it, from the first nonzero row `first` of its
# columns of `T`. The inverse of a shifted 2 by 2 block is its adjugate
# over `(s - λ1)(s - λ2)`, divided a factor at a time: expanded, the
# determinant would lose a narrow resonance's width to the square of the
# frequency, and overflow where that square does.
function quasitriangularsolve!(X::AbstractMatrix, rf::ResolventFactors, s::Number)
    T, poles, first = rf.T, rf.poles, rf.first
    k = size(T, 1)
    @inbounds while k >= 1
        if k > 1 && !iszero(T[k, k - 1])
            a11, a12, a21, a22 = s - T[k, k], T[k - 1, k], T[k, k - 1], s - T[k - 1, k - 1]
            f1, f2 = inv(s - poles[k - 1]), inv(s - poles[k])
            for q in axes(X, 1)
                x1, x2 = X[q, k - 1], X[q, k]
                X[q, k - 1], X[q, k] = (a11*x1 + a12*x2)*f1*f2, (a21*x1 + a22*x2)*f1*f2
            end
            for j in min(first[k - 1], first[k]):k - 2
                t1, t2 = T[j, k - 1], T[j, k]
                for q in axes(X, 1)
                    X[q, j] += t1*X[q, k - 1] + t2*X[q, k]
                end
            end
            k -= 2
        else
            f = inv(s - T[k, k])
            for q in axes(X, 1)
                X[q, k] *= f
            end
            for j in first[k]:k - 1
                t = T[j, k]
                for q in axes(X, 1)
                    X[q, j] += t*X[q, k]
                end
            end
            k -= 1
        end
    end
    return X
end

# the eigenvalues `α ± i sqrt(-b c)` of the standardized 2 by 2 block
# `[α b; c α]`: a rotation's `|b|` exactly; otherwise the root of the
# product, rounded twice, wherever the product is a normal number, and
# beyond that the product of the roots, which neither overflows nor
# underflows
function blockpoles(α, b, c)
    p = abs(b)*abs(c)
    β = abs(b) == abs(c) ? abs(b) :
        floatmin(p) <= p <= floatmax(p) ? sqrt(p) : sqrt(abs(b))*sqrt(abs(c))
    return complex(α, β), complex(α, -β)
end

# the eigenvalues of a real Schur form, read off its diagonal blocks
function schureigenvalues(T::AbstractMatrix)
    λs = Complex{Float64}[]
    k = 1
    nz = size(T, 1)
    while k <= nz
        if k < nz && !iszero(T[k + 1, k])
            push!(λs, blockpoles(T[k, k], T[k, k + 1], T[k + 1, k])...)
            k += 2
        else
            push!(λs, T[k, k])
            k += 1
        end
    end
    return λs
end

# The scratch of the solves on one set of factors of a square block, the
# solution a row per port and a column per state, which a caller holds so
# that concurrent callers do not share it, and the real and imaginary
# parts of the solution and of the output, which the output matrix, being
# real, multiplies as two real products
struct ResolventWorkspace
    solution::Matrix{Complex{Float64}}
    re::Matrix{Float64}
    im::Matrix{Float64}
    outre::Matrix{Float64}
    outim::Matrix{Float64}
end
function ResolventWorkspace(rf::ResolventFactors)
    nz, n = size(rf.ZtB)
    return ResolventWorkspace(Matrix{Complex{Float64}}(undef, n, nz), Matrix{Float64}(undef, n, nz),
        Matrix{Float64}(undef, n, nz), Matrix{Float64}(undef, n, n), Matrix{Float64}(undef, n, n))
end

# `work.solution = ((i w I - T)^-1 Z' B)'`
schurresolvent!(work::ResolventWorkspace, rf::ResolventFactors, w) =
    quasitriangularsolve!(transpose!(work.solution, rf.ZtB), rf, im*w)

# `work.outre + i work.outim = C Z (i w I - T)^-1 Z' B`, the states' part
# of the response at `w`
function resolventoutput!(work::ResolventWorkspace, rf::ResolventFactors, w)
    X = schurresolvent!(work, rf, w)
    work.re .= real.(X)
    work.im .= imag.(X)
    mul!(work.outre, rf.CZ, transpose(work.re))
    mul!(work.outim, rf.CZ, transpose(work.im))
    return work
end

function rationaltransfer!(dest, rf::ResolventFactors, D, w, work::ResolventWorkspace)
    resolventoutput!(work, rf, w)
    dest .= complex.(work.outre, work.outim) .+ D
    return dest
end

# The response of a realization whose state matrix is block diagonal in
# the blocks the fitter writes (see realization), a real pole `a` as `a`
# and a complex one `a = α + iβ` as `[α -β; β α]`, each pole a run of
# equal blocks, as many as its residue's rank. Over its columns `c` of
# `C` and rows `b'` of `B`, a run of real blocks is `P/(s - a)` with
# `P = sum c b'`, and a run of complex ones `R/(s - a) + conj(R)/(s - conj(a))`
# with `R = sum (c1 - i c2)(b1 + i b2)'/2`: the fitter's pole-residue form,
# exact, which is evaluated in the fitter's real basis (basisat!), whose
# factors `1/(s - a)` lose nothing however narrow the resonance. A
# frequency then costs the ports squared per pole rather than per state,
# and a sweep is one product of the coefficients `columns`, an `n x n`
# matrix a column (`P`, or `Re R` and `Im R`), with the basis. Only runs
# of two blocks or more are grouped, each column doing the work of as
# many states; a single block, a residue of rank one, costs the ports
# squared a frequency through its column or its states alike, and is left
# to its states, which hold two rows of the ports rather than their
# square. `poles` is laid out as the basis takes them.
struct PoleGroups
    poles::Vector{ComplexF64}
    columns::Matrix{Float64}
end
# the pole groups of a realization and the states left to the Schur path:
# the single blocks, or every state where `A` is not block diagonal in
# the fitter's blocks
function polegroups(A::Matrix{Float64}, B::Matrix{Float64}, C::Matrix{Float64})
    nz, n = size(A, 1), size(C, 1)
    ungrouped = (PoleGroups(ComplexF64[], zeros(n*n, 0)), collect(1:nz))
    # the blocks of the diagonal, real poles and rotations, and nothing
    # off them
    starts, poles = Int[], ComplexF64[]
    k = 1
    while k <= nz
        push!(starts, k)
        if k < nz && !iszero(A[k + 1, k])
            α, β = A[k, k], A[k + 1, k]
            A[k + 1, k + 1] == α && A[k, k + 1] == -β && !isrealpole(complex(α, β)) || return ungrouped
            push!(poles, complex(α, β))
            k += 2
        else
            push!(poles, A[k, k])
            k += 1
        end
    end
    blockof = zeros(Int, nz)
    for (b, k) in enumerate(starts)
        blockof[k:(b < length(starts) ? starts[b + 1] - 1 : nz)] .= b
    end
    for j in 1:nz, i in 1:nz
        blockof[i] == blockof[j] || iszero(A[i, j]) || return ungrouped
    end
    # the runs of equal poles: a single block's states are left out, and a
    # longer run takes a column per real pole and two per complex
    runs, states = UnitRange{Int}[], Int[]
    b = 1
    while b <= length(starts)
        e = b
        while e < length(starts) && poles[e + 1] == poles[b]
            e += 1
        end
        e > b ? push!(runs, b:e) : append!(states, starts[b]:starts[b] + (isreal(poles[b]) ? 0 : 1))
        b = e + 1
    end
    columns = sum(r -> isreal(poles[first(r)]) ? 1 : 2, runs; init = 0)
    X = zeros(n*n, columns)
    layout = ComplexF64[]
    col = 0
    for r in runs
        a, first_ = poles[first(r)], starts[r]
        if isreal(a)
            X[:, col + 1] .= vec(C[:, first_]*B[first_, :])
            push!(layout, a)
            col += 1
        else
            C1, C2, B1, B2 = C[:, first_], C[:, first_ .+ 1], B[first_, :], B[first_ .+ 1, :]
            X[:, col + 1] .= vec(C1*B1 .+ C2*B2) ./ 2
            X[:, col + 2] .= vec(C1*B2 .- C2*B1) ./ 2
            push!(layout, a, conj(a))
            col += 2
        end
    end
    return PoleGroups(layout, X), states
end

"""
    RationalScatteringProvider(A, B, C, D)

A scattering matrix given as the real state space realization
`S(s) = D + C (s I - A)^(-1) B` of a passive rational multiport, the
form a vector fit of measured or simulated scattering data takes and the
one the transient solver realizes in time, `dz/dt = A z + B a`,
`b = C z + D a` on the incident and reflected power waves. The provider
holds a copy of the realization and takes the real Schur form of `A`
when it is built, `A` itself where it is already quasi-triangular in
LAPACK's standard form, so a signed angular frequency is evaluated by a
quasi-triangular solve on those factors, without forming an inverse;
the realization is fixed from then on. A realization whose `A` is block
diagonal in real poles and rotations `[α -β; β α]`, each pole a run of
equal blocks as many as its residue's rank, as a fitted block's is, is
evaluated pole by pole instead, in the fit's pole-residue form, at the
cost of the ports squared per pole rather than per state, for each pole
whose residue has rank above one; the Schur form is then taken of the
states of the others alone. Built by
[`RationalScattering`](@ref), which validates it.
"""
struct RationalScatteringProvider <: AbstractMatrixProvider
    A::Matrix{Float64}
    B::Matrix{Float64}
    C::Matrix{Float64}
    D::Matrix{Float64}
    # what every evaluation reads, either of which may be empty: the runs
    # of equal poles grouped (see PoleGroups), and the real Schur form of
    # the other states (see ResolventFactors)
    groups::PoleGroups
    states::ResolventFactors
    function RationalScatteringProvider(A::AbstractMatrix, B::AbstractMatrix,
            C::AbstractMatrix, D::AbstractMatrix)
        Am, Bm, Cm, Dm = Matrix{Float64}(A), Matrix{Float64}(B), Matrix{Float64}(C), Matrix{Float64}(D)
        groups, rest = polegroups(Am, Bm, Cm)
        states = length(rest) == size(Am, 1) ? resolventfactors(Am, Bm, Cm) :
            resolventfactors(Am[rest, rest], Bm[rest, :], Cm[:, rest])
        return new(Am, Bm, Cm, Dm, groups, states)
    end
end

# `dest[:, :, i] += c (D + G(i ws[i]))`, `G` the grouped poles' part of
# the response, or `=` where `add` is false; `tile` frequencies to a
# product
function addgrouped!(dest::AbstractArray{<:Any,3}, g::PoleGroups, D::Matrix{Float64},
        ws::AbstractVector, c::Number; tile::Int = 64, add::Bool = true)
    n, N, L = size(D, 1), length(g.poles), min(tile, length(ws))
    phi, basis, out = zeros(ComplexF64, N), zeros(N, 2L), zeros(n*n, 2L)
    for lo in 1:tile:length(ws)
        ids = lo:min(lo + tile - 1, length(ws))
        m = length(ids)
        # the basis's real parts in the first `m` columns and its
        # imaginary parts in the next, so that one product reads the
        # coefficients once for both
        for (j, i) in enumerate(ids)
            basisat!(phi, g.poles, ws[i])
            basis[:, j] .= real.(phi)
            basis[:, m + j] .= imag.(phi)
        end
        mul!(view(out, :, 1:2m), g.columns, view(basis, :, 1:2m))
        @inbounds for (j, i) in enumerate(ids), q in 1:n, r in 1:n
            e = (q - 1)*n + r
            value = c*(D[r, q] + complex(out[e, j], out[e, m + j]))
            dest[r, q, i] = add ? dest[r, q, i] + value : value
        end
    end
    return dest
end
providersize(p::RationalScatteringProvider) = size(p.D, 1)
function evaluateprovider!(dest::AbstractArray{T,3},
        p::RationalScatteringProvider, ws::AbstractVector) where T
    checkdestsize(dest, providersize(p), length(ws))
    return addrational!(dest, p, ws, 1; add = false)
end

# `dest[:, :, i] += c S(i ws[i])`, or `=` where `add` is false, for the
# realization of `p`: its pole groups and its other states through their
# Schur factors, each where it has poles, the first with the feedthrough
# and as `add` says, the second added to it
function addrational!(dest::AbstractArray{<:Any,3}, p::RationalScatteringProvider,
        ws::AbstractVector, c::Number; add::Bool = true)
    g, rf, n = p.groups, p.states, size(p.D, 1)
    grouped = !isempty(g.poles)
    grouped && addgrouped!(dest, g, p.D, ws, c; add)
    grouped && isempty(rf.poles) && return dest
    d, add = grouped ? (zero(c), true) : (c, add)
    work = ResolventWorkspace(rf)
    for i in eachindex(ws)
        resolventoutput!(work, rf, ws[i])
        @inbounds for q in 1:n, r in 1:n
            value = d*p.D[r, q] + c*complex(work.outre[r, q], work.outim[r, q])
            dest[r, q, i] = add ? dest[r, q, i] + value : value
        end
    end
    return dest
end

# The passivity of a rational realization over every frequency: its
# feedthrough, which is its value at infinite frequency, and then the
# verdict of passivityassessment at the level of `atol` (see
# passivelevel), whose crossing test finds a peak however narrow where a
# sample can miss one. A block it cannot settle either way is accepted
# (see RationalScattering). `normtested` leaves the test out,
# for a fit whose passivity the fitter has tested on its residues (see
# fitsampled).
function checkpassive(p::RationalScatteringProvider; atol = 1e-8, normtested::Bool = false)
    A, B, C, D = p.A, p.B, p.C, p.D
    margin = passivitymargin(D)
    margin < -atol && throw(ArgumentError(lazy"The rational scattering block is not passive at infinite frequency: the minimum eigenvalue of I - D*D' is $(margin)."))
    (size(A, 1) == 0 || normtested) && return nothing
    (; verdict, lower, frequency) = passivityassessment(A, B, C, D; atol = atol)
    verdict === :active && throw(ArgumentError(lazy"The rational scattering block is not passive: its largest singular value over all frequencies is at least $(lower), at $(frequency) rad/s, which is above the tolerance $(atol)."))
    return nothing
end

# the realization in the frequency unit of the poles, `S(s) = D + C (s/w I - A/w)^(-1) B/w`,
# with its states scaled by the similarity `A -> T^(-1) A T`, `B -> T^(-1) B`,
# `C -> C T`, which leaves `S` alone; and the scale. The scales balance the
# system matrix `[A B; C 0]` by the iteration of LAPACK's gebal over the
# states, as SLICOT's TB01ID does for the norm of AB13DD: each state's row,
# its input and its couplings from the other states, against its column,
# its output and its couplings to them, the diagonal aside, so that no
# coupling is scaled out of proportion to the rest of the matrix, whose
# Schur form would then lose a narrow resonance. A state's scale is
# doubled or halved while that brings the two norms together, and kept
# where it takes their sum under `improvement` of what it was; the sweeps
# over the states end when one changes none, after `sweeps` at most. A
# state with no column, which no output and no other state sees, or with
# no row, which no input and no other state drives, has no part in `S`
# and no scale which balances it: the rest of it is cleared. Each scale
# is a power of two, which scales exactly, so that a state matrix in real
# Schur form, as a fit's is, stays in it to the bit and is not factored
# again. A step's scale and the norms it moves are held in the range
# gebal holds them in, so that no entry overflows and no row or column
# of a state with a part in `S` underflows to nothing, which would clear
# it; a norm of finite entries past floatmax is taken at floatmax, which
# sets a first step, and the next sweep balances the rest.
function balancedrealization(A, B, C; improvement::Real = 0.95, sweeps::Integer = 100)
    nz = size(A, 1)
    wscale = max(maximum(abs, isquasitriangular(A) ? schureigenvalues(A) : eigvals(A); init = 0.0), floatmin(Float64))
    An, Bn, Cn = A ./ wscale, B ./ wscale, Matrix{Float64}(C)
    tiny = 2floatmin(Float64)/eps(Float64)
    for _ in 1:sweeps
        changed = false
        for k in 1:nz
            c = hypot(norm(view(An, 1:k - 1, k)), norm(view(An, k + 1:nz, k)), norm(view(Cn, :, k)))
            r = hypot(norm(view(An, k, 1:k - 1)), norm(view(An, k, k + 1:nz)), norm(view(Bn, k, :)))
            c, r = min(c, floatmax(Float64)), min(r, floatmax(Float64))
            if iszero(c) || iszero(r)
                d = An[k, k]
                if !iszero(r)
                    view(An, k, :) .= 0
                    view(Bn, k, :) .= 0
                    changed = true
                elseif !iszero(c)
                    view(An, :, k) .= 0
                    view(Cn, :, k) .= 0
                    changed = true
                end
                An[k, k] = d
                continue
            end
            f, s = 1.0, c + r
            while c < r/2 && max(f, c) < 1/tiny && r > tiny
                f, c, r = 2f, 2c, r/2
            end
            while c/2 >= r && r < 1/tiny && min(f, c) > tiny
                f, c, r = f/2, c/2, 2r
            end
            c + r < improvement*s || continue
            changed = true
            d = An[k, k]
            view(An, k, :) ./= f
            view(An, :, k) .*= f
            An[k, k] = d
            view(Bn, k, :) ./= f
            view(Cn, :, k) .*= f
        end
        changed || break
    end
    return An, Bn, Cn, wscale
end

# The poles once each. An n-port's fit repeats every pole once per port
# among the eigenvalues of its state matrix, and a probe or a bracket
# repeated is work repeated; eigenvalues within the square root of eps of
# each other, relative to their size, are the same pole.
function distinctpoles(poles; tol::Real = sqrt(eps()))
    out = ComplexF64[]
    for l in sort(complex.(poles); by = x -> (imag(x), real(x)))
        (isempty(out) || abs(l - last(out)) > tol*abs(l)) && push!(out, l)
    end
    return out
end

# The searches for where a rational block peaks, passivityassessment's and
# the grid of the passivity enforcement, bracket each complex pole over
# `polespan` of its half widths either side, and spread their grids over
# the poles' magnitudes to `gridpad` times past the outermost either side.
const polespan = 8.0
const gridpad = 10.0
# the golden section steps which refine each of passivityassessment's
# brackets
const peakrefinements = 32
# how far under a level the feedthrough must stand for the crossings of the
# level to be taken from the Hamiltonian matrix rather than the pencil (see
# unitcrossings)
const hamiltonianmargin = 1e-3

# The span of `span` half widths either side of the complex pole `l`,
# where its resonance peaks: a resonance narrower than the spacing of a
# grid spread over the band peaks between its points, and is found by
# looking where its pole is
function polebracket(l, span::Real)
    w0, half = imag(l), max(abs(real(l)), eps())
    return max(w0 - span*half, 0.0), w0 + span*half
end

"""
    passivityassessment(A, B, C, D; atol = 1e-8)

Whether the real rational block `S(s) = D + C (s I - A)^(-1) B` is
passive to within `atol`, its dissipation `I - S'S` no lower than
`-atol` at any frequency, the tolerance a block's data is held to: its
largest singular value at most `sqrt(1 + atol)`, the level. The answer
is the named tuple `(; verdict, lower, upper, frequency)`: `verdict` one
of `:passive`, `:active` and `:indeterminate`, `lower` the largest
singular value found, at `frequency` in rad/s, and `upper` a level no
singular value reaches, `sqrt(1 + atol)` where the block is passive and
`Inf` where none was established.

The largest singular value is sought where a peak can stand: at the
feedthrough, at zero, at the frequency of each pole and on a grid
spanning the poles' magnitudes, and by a golden section search either
side of each complex pole and between the neighbours of the best sample,
since a narrow resonance peaks near its pole but not at it. A peak can
be narrower than any search finds, so the frequencies where a singular
value equals the level are found as well, as the imaginary eigenvalues
of the Hamiltonian matrix of `S` over that level where the feedthrough
stands well under it, and otherwise of a pencil which inverts nothing,
so that a feedthrough at the level is no obstacle (Boyd and
Balakrishnan), and the largest singular value between consecutive
crossings is measured. A block is `:active` where a value found exceeds
the level, which settles it, and `:passive` where none does and the
level has no crossing, which settles it too.

Otherwise it is `:indeterminate`. The crossing test resolves a pair of
crossings only so far as the roundoff lets their eigenvalues stand apart
on the axis: a peak above the level by less than that brings its two
crossings together into a nearly double eigenvalue which leaves the
axis, and a narrow resonance under the level can leave eigenvalues near
enough to the axis to be taken for crossings no value between them
reaches. The caller is told so rather than given a side. The frequency
axis alone is assessed: the stability of `A`, which passivity needs as
well, is the caller's to check.
"""
function passivityassessment(A, B, C, D; atol::Real = 1e-8)
    level = passivelevel(atol)
    An, Bn, Cn, wscale = balancedrealization(A, B, C)
    # the largest singular value at `w` in the unit of the balanced
    # realization, through its Schur form taken once (see
    # ResolventFactors), and `Inf` where the realization is not finite, at
    # a pole on the axis
    rf = resolventfactors(An, Bn, Cn)
    work, F = ResolventWorkspace(rf), similar(D, Complex{Float64})
    sigma = w -> (rationaltransfer!(F, rf, D, w, work); all(isfinite, F) ? opnorm(F) : Inf)
    worst, where = opnorm(D), Inf
    # The grid spans the magnitudes of the poles the system has, not
    # fixed decades around one: a peak can lie far from every pole
    # frequency, as `k s/((s + a)(s + b))` peaking at `sqrt(a b)` shows.
    # It has two points per distinct pole, nine at least, and the probes
    # at the poles' frequencies and the brackets below are taken once per
    # distinct pole as well: a multiport's states repeat a pole once per
    # rank of its residue, which adds no peak to find.
    distinct = distinctpoles(rf.poles)
    mags = [abs(l) for l in rf.poles if abs(l) > 0]
    glo = isempty(mags) ? 1e-2 : minimum(mags)/gridpad
    ghi = isempty(mags) ? 1e2 : maximum(mags)*gridpad
    probes = sort!(vcat(0.0, [abs(imag(l)) for l in distinct if abs(imag(l)) > 0],
        exp.(range(log(glo), log(ghi); length = max(9, 2*length(distinct))))))
    for w in probes
        s = sigma(w)
        s > worst && ((worst, where) = (s, w))
    end
    # A golden section search around each complex pole, and between the
    # neighbours of the best probe, which covers a peak the grid only
    # straddled: a real pole has no resonance of its own, and a peak it
    # takes part in need not be near it.
    brackets = Tuple{Float64,Float64}[polebracket(l, polespan) for l in distinct if imag(l) > 0]
    if isfinite(where)
        j = searchsortedfirst(probes, where)
        push!(brackets, (probes[max(j - 1, 1)], probes[min(j + 1, length(probes))]))
    end
    for (a, b) in brackets
        b > a || continue
        φ = (sqrt(5) - 1)/2
        u, v = b - φ*(b - a), a + φ*(b - a)
        fu, fv = sigma(u), sigma(v)
        for _ in 1:peakrefinements
            if fu > fv
                b, v, fv = v, u, fu
                u = b - φ*(b - a); fu = sigma(u)
            else
                a, u, fu = u, v, fv
                v = a + φ*(b - a); fv = sigma(v)
            end
        end
        fu > worst && ((worst, where) = (fu, u))
        fv > worst && ((worst, where) = (fv, v))
    end
    worst > level && return (; verdict = :active, lower = worst, upper = Inf, frequency = where*wscale)
    # The pencil at a level above every value found is regular, a
    # singular value equal to the level everywhere being impossible.
    # Between consecutive crossings as many singular values stand above
    # the level at every frequency, and below the first and past the last
    # none do, the value at zero and the feedthrough having been measured.
    found = unitcrossings(An, Bn, Cn ./ level, D ./ level)
    for k in 1:length(found) - 1
        w = (found[k] + found[k + 1])/2
        s = sigma(w)
        s > worst && ((worst, where) = (s, w))
    end
    worst > level && return (; verdict = :active, lower = worst, upper = Inf, frequency = where*wscale)
    passive = isempty(found) && !isnan(worst)
    return (; verdict = passive ? :passive : :indeterminate, lower = worst, upper = passive ? level : Inf,
        frequency = where*wscale)
end

# The frequencies, in the unit of the balanced realization, at which a
# singular value of `S(s) = D + Cn (s I - An)^(-1) Bn` equals one, sorted:
# the imaginary eigenvalues of the Hamiltonian matrix of `S` where the
# feedthrough stands at least `hamiltonianmargin` under one, and otherwise
# the finite eigenvalues on the imaginary axis of the pencil of the
# equations `i w x = A x + B u`, `-i w y = A' y + C' w`, `w = C x + D u`,
# `u = B' y + D' w`, which say `S(i w)' S(i w) u = u`. The Hamiltonian
# matrix inverts `I - D' D` and `I - D D'`: it is a standard eigenproblem
# of twice the states, several times faster than the pencil's generalized
# one of twice the states and ports and accurate while those inverses are
# well conditioned, but as the feedthrough nears one it moves the
# eigenvalues at the crossings off the axis and loses them. The pencil's
# matrices are formed without inverting anything, so a feedthrough on the
# unit circle is no obstacle to it. The crossings of `S` over a level are
# the crossings of that level.
function unitcrossings(An, Bn, Cn, D)
    1 - opnorm(D) >= hamiltonianmargin && return hamiltoniancrossings(An, Bn, Cn, D)
    return pencilcrossings(An, Bn, Cn, D)
end

# the crossings of one by the imaginary eigenvalues of the Hamiltonian
# matrix of the bounded real lemma, for a feedthrough under one
function hamiltoniancrossings(An, Bn, Cn, D)
    R = cholesky(Symmetric(I - transpose(D)*D))
    Q = cholesky(Symmetric(I - D*transpose(D)))
    F = An + Bn*(R \ (transpose(D)*Cn))
    H = [F Bn*(R \ transpose(Bn)); -transpose(Cn)*(Q \ Cn) -transpose(F)]
    return axiscrossings(eigvals(H))
end

# the crossings of one by the pencil
function pencilcrossings(An, Bn, Cn, D)
    nz, m = size(An, 1), size(D, 1)
    Z = zeros
    H = [An Z(nz, nz) Bn Z(nz, m); Z(nz, nz) transpose(An) Z(nz, m) transpose(Cn);
        Cn Z(m, nz) D -Matrix(1.0I, m, m); Z(m, nz) transpose(Bn) -Matrix(1.0I, m, m) transpose(D)]
    E = Matrix(Diagonal(vcat(ones(nz), -ones(nz), zeros(2m))))
    # the eigenvalues by LAPACK's ggev, whose QZ is unblocked: the blocked
    # QZ of ggev3, which eigvals(H, E) calls, writes past the end of its
    # eigenvalue arrays at some orders of this pencil
    alphar, alphai, beta = LAPACK.ggev!('N', 'N', H, E)
    return axiscrossings(complex.(alphar, alphai) ./ beta)
end

# the frequencies of the finite eigenvalues on the imaginary axis, in the
# unit of the balanced realization: an eigenvalue of magnitude up to
# `largest` whose real part is within `axistol` of its magnitude plus one,
# a crossing and its conjugate being one, as are two crossings within
# `mergetol` of their frequency plus one
function axiscrossings(lambda; largest::Real = 1e8, axistol::Real = 1e-8, mergetol::Real = 1e-9)
    finite = [l for l in lambda if isfinite(real(l)) && isfinite(imag(l)) && abs(l) <= largest]
    crossings = sort!([abs(imag(l)) for l in finite if abs(real(l)) <= axistol*(abs(l) + 1)])
    merged = Float64[]
    for w in crossings
        (isempty(merged) || w - last(merged) > mergetol*(w + 1)) && push!(merged, w)
    end
    return merged
end

"""
    RationalScattering(A, B, C, D; zref = 50.0, grounded = true,
        noise = Passive(), atol = 1e-8)

A [`ScatteringParameters`](@ref) block from the real state space
realization `S(s) = D + C (s I - A)^(-1) B` of a passive rational
multiport, with `A` the `nz` by `nz` state matrix, `B` `nz` by `nports`,
`C` `nports` by `nz` and `D` `nports` by `nports`, all real and finite,
`A` stable. The block is validated as passive to `atol` over every
frequency, its dissipation `I - S'S` no lower than `-atol` as its
samples would be, by a search for its peaks and the crossing test of
the level `sqrt(1 + atol)` (see [`passivityassessment`](@ref)), and it
is rejected where a singular value is found above that level, unless it
states its noise with a [`NoiseCovariance`](@ref), which is how an
active block, an amplifier given by its scattering parameters, declares
it; declared [`Lossless`](@ref), it is held to singular values whose
squares stand within `atol` of one at every frequency, by the same test
on it and on its inverse. A block the test cannot settle either way is
accepted: the crossing test cannot tell a peak above the level by less
than its roundoff from one just under it, and a lossless block, within
`atol` of the level at every frequency, often leaves it undecided.
Stability is required of every realization. It is evaluated by the harmonic balance solvers
at every frequency and realized in time by the transient solver with
its states, so the two describe the same block, and its noise is the
noise of its loss at every frequency by Bosma's relation, or the noise
it states. The realization is what a vector fit of measured or
simulated data delivers; a lossless line is [`TransmissionLine`](@ref)
instead, which needs no states.
"""
RationalScattering(A, B, C, D; zref = 50.0, grounded::Bool = true, noise = Passive(), atol::Real = 1e-8) =
    rationalblock(A, B, C, D; zref, grounded, noise, atol)

# The block of a realization, validated as the constructor above says.
# `normtested` is for a fit whose passivity the fitter has tested on its
# residues (see fitsampled), which the validation would repeat; its
# feedthrough is still checked.
function rationalblock(A, B, C, D; zref, grounded::Bool, noise, atol::Real,
        normtested::Bool = false)
    Am, Bm, Cm, Dm = Matrix{Float64}(A), Matrix{Float64}(B), Matrix{Float64}(C), Matrix{Float64}(D)
    n, nz = size(Dm, 1), size(Am, 1)
    size(Dm) == (n, n) && size(Am) == (nz, nz) && size(Bm) == (nz, n) && size(Cm) == (n, nz) || throw(DimensionMismatch(
        lazy"the realization needs A of size (nz, nz), B (nz, nports), C (nports, nz) and D (nports, nports); got $(size(Am)), $(size(Bm)), $(size(Cm)), $(size(Dm))."))
    all(M -> all(isfinite, M), (Am, Bm, Cm, Dm)) || throw(ArgumentError("the realization must be finite."))
    provider = RationalScatteringProvider(Am, Bm, Cm, Dm)
    # the eigenvalues are the poles the provider holds, in its groups and
    # its Schur factors
    abscissa = maximum(real, [provider.groups.poles; provider.states.poles]; init = -Inf)
    abscissa < 0 || throw(ArgumentError(
        lazy"the realization is unstable: the largest real part of an eigenvalue of A is $(abscissa) per second."))
    return checkedblock(provider, n, zrefvector(zref, n), grounded, noise,
        ConjugateSymmetry(), Pair{Symbol,Any}[], ScatteringLimit(), atol; normtested)
end

"""
    unitaritybound(provider)

An upper bound on how far the data of `provider` is from unitary at any
frequency, as the largest entry of `I - S S'`, from what it stores: its
samples (see [`unitaritydeviation`](@ref)), and between the knots of a
table the bound of its interpolant. `Inf` where nothing stored bounds
it: a callable or a realization, whose values are computed, a table
extrapolated linearly, which is unbounded beyond its knots, and data
which is zero somewhere, a table extrapolated by zero or a piecewise
table, which absorbs everything there. A lossless line is unitary by
construction. A rotation of unit phases on each side keeps the size of
every entry of `I - S S'`, so it is bounded as what it turns is.
"""
function unitaritybound(provider::AbstractMatrixProvider)
    p = unrotated(provider)
    isstored(p) || return Inf
    p isa PiecewiseTabulatedProvider && return Inf
    p isa TabulatedMatrixProvider && p.extrapolation in (:linear, :zero) && return Inf
    fs, _ = samplefrequencies(p)
    S, _ = storedsamples(p, fs)
    worst = maximum(k -> unitaritydeviation(view(S, :, :, k)), axes(S, 3); init = 0.0)
    p isa TabulatedMatrixProvider && (worst = max(worst, interpolantdeviation(p)))
    return worst
end

# How far a table's interpolant can be from unitary between its knots.
# Two unitary knots do not make a unitary interpolant: `S = 1` and
# `S = -1` interpolate to a perfect absorber halfway between. On the
# segment between two knots the interpolant is a polynomial in the
# fraction `t` of the segment, of degree one for the chord and three for
# the spline, and `t` is real, so every entry of `I - S S'` is a
# polynomial of twice that degree, `m`. A polynomial of degree `m` is
# bounded on the segment by its largest value at the `m + 1` Chebyshev
# points times their Lebesgue constant, which is at most
# `(2/pi) log(m + 1) + 1` (Rivlin), 2.24 for the spline: the bound is
# rigorous, and within that factor of the deviation it bounds, which for
# a spline through unitary samples falls as the fourth power of their
# spacing.
function interpolantdeviation(p::TabulatedMatrixProvider)
    m = p.interpolation == :linear ? 2 : 6
    points = [(1 - cospi((2j + 1)/(2m + 2)))/2 for j in 0:m]
    lebesgue = 2/pi*log(m + 1) + 1
    f = p.frequencies
    n = size(p.values, 1)
    S = Matrix{eltype(p.values)}(undef, n, n)
    worst = 0.0
    for k in 1:length(f) - 1, t in points
        w = f[k] + t*(f[k + 1] - f[k])
        p.interpolation == :linear ? lerpslices!(S, p, k, k + 1, w) : splineslices!(S, p, k, w)
        worst = max(worst, unitaritydeviation(S))
    end
    return lebesgue*worst
end

"""
    provablylossless(block::ScatteringParameters)

Whether the scattering data can be shown to be unitary at every frequency to
the block's `atol` from the data alone (see [`unitaritybound`](@ref)): a
constant and a table can be, a lossless line is by construction, and a
callable, whose values away from any sampled frequency are unknown, and a
realization, whose norms are costly to bound, are not.

A block which is not provably lossless carries noise channels
([`ScatteringNoisePlan`](@ref)), and those of a block which is in fact
lossless are identically zero, so a `false` costs work. A `true` leaves
out channels whose commutator `I - S S'` is within the block's `atol` of
zero, the tolerance a block declared [`Lossless`](@ref) is held to.
"""
Base.@nospecializeinfer provablylossless(@nospecialize(block::ScatteringParameters)) =
    unitaritybound(block.provider) <= block.atol

# A rational block is lossless when every singular value of `S` is one at
# every frequency: the largest at most one, which is its passivity, and
# the smallest at least one, which is the largest singular value of
# `S^(-1)` at most one, a rational block of its own when the feedthrough
# is invertible; without an invertible feedthrough the block is not
# lossless. Declared lossless to `atol`, every squared singular value
# stands within `atol` of one, as a scalar sample's does: the largest
# singular value at most `sqrt(1 + atol)`, the level of `atol` (see
# passivelevel), to which checkpassive holds every block that does not
# state its noise, or the fitter a fit (see rationalblock), and the
# smallest at least `sqrt(1 - atol)`, the inverse's largest at most
# `1/sqrt(1 - atol)`, the level of `atol/(1 - atol)`; an `atol` of one or
# more bounds nothing from below. This holds the smallest, with the
# feedthrough's unitarity, and refuses the inverse where
# passivityassessment finds a value above its level, whose crossing test
# finds a notch however narrow. No test of the coefficients of
# `I - S(-s)' S(s)` can: a notch of relative width `eps` has coefficients
# of order `eps^2` there and reaches one at its center, so only the
# values on the axis tell. Nothing is inferred about a rational block by
# default: it keeps the channels of its loss, however small, and only a
# declaration of `Lossless()` is validated, which this can only refuse.
function losslessnorms(p::RationalScatteringProvider; atol = 1e-8)
    unitaritydeviation(p.D) <= atol || return false
    (size(p.A, 1) == 0 || atol >= 1) && return true
    inverse = atol/(1 - atol)
    Dinv = inv(p.D)
    return passivityassessment(p.A - p.B*Dinv*p.C, p.B*Dinv, -Dinv*p.C, Dinv;
        atol = inverse).lower <= passivelevel(inverse)
end

# The contract of a block between and beyond its samples, where its data
# determines it: a realization over every frequency, passive unless it
# states its noise and unitary where it declares itself lossless, and a
# table declared lossless between its knots and beyond them. A callable
# is known only where it is evaluated.
function checkbeyondsamples(block::ScatteringParameters; normtested::Bool = false)
    p, noise, atol = block.provider, block.noise, block.atol
    if p isa RationalScatteringProvider
        noise isa NoiseCovariance || checkpassive(p; atol, normtested)
        (noise isa Lossless && !losslessnorms(p; atol = atol)) && throw(ArgumentError(lazy"noise = Lossless() says the scattering matrix is unitary at every frequency, but this rational block's smallest singular value falls below sqrt(1 - atol) on the frequency axis, infinite frequency included, for its atol of $(atol). Use the default Passive() noise model, which gives a dissipative block the noise its loss requires."))
    elseif noise isa Lossless && isstored(p)
        bound = unitaritybound(p)
        bound <= atol || throw(ArgumentError(lazy"noise = Lossless() says the scattering matrix is unitary at every frequency, but this block's data is unitary only at its samples: $(beyondsamples(p, bound, atol)). Use the default Passive() noise model, which gives a dissipative block the noise its loss requires."))
    end
    return nothing
end

# where the data of a table stops being unitary, for the message above
function beyondsamples(p, bound, atol)
    p isa PiecewiseTabulatedProvider && return "between and beyond its bands it is zero, which absorbs everything"
    p.extrapolation == :zero && return "beyond its knots it is zero, which absorbs everything"
    p.extrapolation == :linear && return "beyond its knots it is continued linearly, which bounds nothing"
    return "between its knots its interpolant may depart from unitary by up to $(bound) in an entry of I - S S', above its atol of $(atol); a table sampled more finely is bounded more tightly"
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
    return evaluatesignedprovider!(dest, block.provider, block.negative_frequency,
        ws, absbuffer)
end

function evaluatesignedprovider!(dest, provider, rule, ws, absbuffer)
    if rule isa Native
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
    return evaluatesignedprovider!(dest, block.noise.provider,
        block.negative_frequency, ws, absbuffer)
end

# The scattering data of a Touchstone file as a table, and the reference
# impedances the file's option line states; an explicit `zref` which
# disagrees with them is an error, since it is ambiguous between correcting
# a mislabeled file and asking for renormalization.
function touchstoneprovider(path::AbstractString, zref; interpolation::Symbol,
        extrapolation::Symbol, form::Symbol)
    # a file is tabulated data, which no form calls
    checknoform(form, "a Touchstone file")
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
    provider = TabulatedMatrixProvider(frequencies, touchstonescattering(ts, path);
        interpolation = interpolation, extrapolation = extrapolation)
    return provider, filezref
end

# The network data of a Touchstone file as scattering parameters at the
# file's reference impedances. Y and Z parameters are converted from
# siemens and ohms with YtoS and ZtoS. A version 2 file states them so; a
# version 1 file states them normalized to its R, `z = Z/R` and `y = Y R`,
# and the loader multiplies every entry of such a file by R, which gives
# back an impedance but leaves an admittance at `Y R^2`. In dB it
# multiplies the decibels by R instead, which leaves nothing the data can
# be read back from, so such a file is refused. Hybrid H and G parameters
# mix impedances, admittances and ratios and are refused.
function touchstonescattering(ts, path)
    parameter = lowercase(ts.parameter)
    N = Array{Complex{Float64},3}(ts.N)
    parameter == "s" && return N
    (parameter in ("z", "y") && ts.version < 2 && lowercase(ts.format) == "db") && throw(ArgumentError(
        lazy"The Touchstone file $(path) holds version 1 $(uppercase(ts.parameter)) parameters in dB, which the loader scales wrongly; write the data as RI or MA, or as a version 2 file."))
    z = collect(Float64, ts.reference)
    parameter == "z" && return ZtoS(N; portimpedances = z)
    parameter == "y" && return YtoS(ts.version < 2 ? N ./ ts.R^2 : N; portimpedances = z)
    throw(ArgumentError(lazy"The Touchstone file $(path) holds $(uppercase(ts.parameter)) parameters, which mix impedances, admittances and ratios; a scattering block reads S, Y or Z parameters, so convert the data to one of those."))
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
# a lossless line is unitary at every frequency, so it may state so
unitaritybound(::TransmissionLineProvider) = 0.0
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
    return checkedblock(provider, 2, fill(Float64(Z0), 2), grounded, noise,
        Native(), Pair{Symbol,Any}[], ScatteringLimit(), 1e-8)
end

# === a pumped scattering block ===

"""
    LinearizedScattering(linearized, wp; ports = nothing, zref = 50.0,
        grounded = true, noise = Lossless(), phase = 0.0,
        interpolation = :cubic, atol = 1e-6, dcmodel = ScatteringLimit(),
        envelope = nothing)
    LinearizedScattering(H, wp; harmonics, nports, zref = 50.0,
        grounded = true, noise = Lossless(), phase = 0.0,
        dcmodel = ScatteringLimit(), envelope = nothing, atol = 1e-6)

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
[`PiecewiseTabulatedProvider`](@ref)): a band is samples no further
apart than the solve's frequencies are, and a solve at one frequency
gives each sample a band of its own. Samples which fall on the same
frequency from different mode pairs must agree to `atol`, relative to
the largest entry, which is what makes the data that of one periodic
steady state, and the block built from the tables is checked against
its declaration over the modes of the solve at every frequency, as
every solve checks it, since the tables hold at one frequency samples
from solves of neighboring signal frequencies whose mode truncations
differ; a device whose mode truncation was too tight fails both, and
`atol` admits the discrepancy of one which was nearly so.
`ports` selects and orders the device's ports which become the block's,
by default all of them, and `zref` gives their reference impedances. A
port left out takes its termination into the block, so a stated
covariance gains the noise that termination sends in, scattered to the
ports the block keeps, at the temperature the solve records for it
(`porttemperatures`).
`phase` rotates `H_k` by `exp(im*k*phase)`, the phase of the block's
pump relative to the one the data was computed with, which matters when
other elements of the circuit share the pump.

The second form takes the harmonic transfer functions directly: `H` is a
vector of providers, one per entry of `harmonics` (nonnegative, ascending,
beginning with zero), each a matrix, a callable of the signed angular
frequency, or a tuple `(frequencies, values)` tabulated over signed
frequencies. `H_0` is the response of a real system,
`H_0(-nu) = conj(H_0(nu))`: a constant one must be real, and a table
holding both signs of a frequency is checked for it there, to `atol`.

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
as `1/2`, it becomes the harmonic covariances
`V_k(nu) = <n(nu + k wp) n(nu)'>`, sampled like the transfer functions,
and is held to the minimum the commutation relations require,
`V - K/2` and `V + K/2` positive semidefinite with `K = J - S J S'`, on the
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
    checkpumpoptions(wp, atol, phase)
    ks = collect(Int, harmonics)
    (!isempty(ks) && ks[1] == 0 && issorted(ks; lt = <=) && all(>=(0), ks)) || throw(ArgumentError(
        "the harmonics must be nonnegative, strictly increasing, and begin with zero."))
    length(H) == length(ks) || throw(DimensionMismatch(lazy"give one provider per harmonic; $(length(H)) providers for $(length(ks)) harmonics."))
    n = Int(nports)
    providers = AbstractMatrixProvider[]
    for h in H
        p = scatteringprovider(h; n = n)
        providersize(p) == n || throw(DimensionMismatch(lazy"a harmonic transfer function has dimension $(providersize(p)) but the block has $(n) ports."))
        push!(providers, p)
    end
    checkrealresponse(providers[1], Float64(atol))
    if noise isa NoiseCovariance
        V = noise.provider
        (V isa AbstractVector && length(V) == length(ks)) || throw(ArgumentError(
            "the stated noise of a pumped block given by its harmonic transfer functions is one covariance provider per harmonic, V_k(nu) = <n(nu + k wp) n(nu)'> in the units of Cnoise."))
        vps = AbstractMatrixProvider[]
        for v in V
            vp = matrixprovider(v, Complex{Float64}; n = n,
                interpolation = noise.interpolation, extrapolation = noise.extrapolation)
            providersize(vp) == n || throw(DimensionMismatch(lazy"a harmonic covariance has dimension $(providersize(vp)) but the block has $(n) ports."))
            push!(vps, vp)
        end
        noise = NoiseCovariance(vps, noise.interpolation, noise.extrapolation, noise.atol, noise.completed, noise.padding)
    else
        noise isa Lossless || throw(ArgumentError("a pumped block declares noise = Lossless(), or states its noise with a NoiseCovariance over its harmonics."))
    end
    checkdcmodel(dcmodel, n, noise, Float64(atol))
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
    isempty(nus) && ks == [0] && all(p -> unrotated(p) isa ConstantMatrixProvider, stored) && push!(nus, Float64(wp))
    declared = declaredtolerance(built)
    for nu in unique!(sort!(nus))
        rows, cols, K = pumpedfamily(built, (nu,))
        v = pumpedviolation(built, rows, cols, K)
        v <= declared || throw(ArgumentError(lazy"the block's data does not meet what it declares: over the modes its harmonics reach from $(nu) rad/s the violation of its losslessness or of the commutation relations of its stated covariance is $(v) of the square of its largest entry, against the $(declared) of atol and its noise model's. A device with loss or gain states its noise with NoiseCovariance; raise atol to admit a discrepancy of the data."))
    end
    return built
end

# the pump frequency, the tolerance and the pump phase of a pumped block
function checkpumpoptions(wp, atol, phase)
    (isfinite(wp) && wp > 0) || throw(ArgumentError("the pump frequency must be positive and finite."))
    checkatol(atol)
    isfinite(phase) || throw(ArgumentError("the pump phase must be finite."))
    return nothing
end

# the tolerance a pumped block's data is held to against what it
# declares: its own, and a stated covariance's where that is looser
declaredtolerance(block::LinearizedScattering) =
    block.noise isa NoiseCovariance ? max(block.atol, block.noise.atol) : block.atol

function LinearizedScattering(lin, wp::Real; ports = nothing, zref = 50.0,
        grounded::Bool = true, noise = Lossless(), phase::Real = 0.0,
        interpolation::Symbol = :cubic, atol::Real = 1e-6,
        dcmodel::AbstractDCModel = ScatteringLimit(), envelope = nothing)
    (hasproperty(lin, :S) && hasproperty(lin, :w)) || throw(ArgumentError(
        "give the linearized output of hbsolve, which carries the scattering matrix and its frequencies, or the harmonic transfer functions with their harmonics."))
    S = lin.S
    S isa AxisKeys.KeyedArray || throw(ArgumentError(
        "the scattering matrix must be keyed by mode and port; solve with keyedarrays = true."))
    checkpumpoptions(wp, atol, phase)
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
            near = (i <= length(direct) && samefrequency(nu, direct[i], wp)) ||
                (i > 1 && samefrequency(nu, direct[i - 1], wp))
            near || push!(samples[0], (nu, M))
        end
    end
    scale = max(1.0, maximum(abs, A))
    # the samples of one band are the signal frequencies, shifted: no two
    # of them further apart than the widest step between those
    step = length(wv) > 1 ? maximum(diff(sort(wv))) : 0.0
    ks = sort!(collect(keys(samples)))
    providers = harmonictables(samples, ks, n, wp, step, interpolation, atol*scale,
        (k, nu, d) -> "the data is not that of one periodic steady state: the samples of the harmonic $(k) at $(nu) rad/s from two mode pairs differ by $(d) against the largest entry $(scale). A tighter mode truncation, or data not from the linearized solve of one pumped device, does this; atol admits a discrepancy.")
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
        # the covariance is tabulated band by band as the transfer
        # functions are, zero between and beyond its bands, where the
        # solve says nothing, so it extrapolates by no other rule
        noise.extrapolation in (:error, :zero) || throw(ArgumentError(
            lazy"the covariance of a solve is zero between and beyond its bands, as the transfer functions are, so it cannot be extrapolated by $(repr(noise.extrapolation)); leave extrapolation at its default."))
        (collect(AxisKeys.axiskeys(C, 1)) == outmodes && collect(AxisKeys.axiskeys(C, 2)) == outports &&
            collect(AxisKeys.axiskeys(C, 3)) == inmodes && collect(AxisKeys.axiskeys(C, 4)) == inports &&
            size(C, 5) == length(wv) && collect(AxisKeys.axiskeys(C, 5)) == collect(AxisKeys.axiskeys(S, 5))) || throw(ArgumentError(
            "the stated covariance does not share the modes, ports and frequencies of the scattering matrix; take both from one solve."))
        Ca = Array(C)
        # the ports the block leaves out become part of it, their
        # terminations with them: what the field each sends in scatters to
        # the ports the block keeps is noise the block emits, at the
        # temperature of that termination (the vacuum where the solve
        # records none)
        temps = hasproperty(lin, :porttemperatures) ? lin.porttemperatures : zeros(length(inports))
        dropped = [(r, Float64(temps[r])) for r in eachindex(inports) if !(inports[r] in selected)]
        vsamples = Dict{Int,Vector{Tuple{Float64,Matrix{Complex{Float64}}}}}()
        for (a, mo) in enumerate(modes), (b, mi) in enumerate(modes), i in eachindex(wv)
            k = mo - mi
            nu = wv[i] + mi*wp
            M = Matrix{Complex{Float64}}(undef, n, n)
            for q in 1:n, p in 1:n
                M[p, q] = Ca[a, po[p], b, qi[q], i]
            end
            for (r, T) in dropped, (c, mc) in enumerate(modes)
                d = thermalnoise(wv[i] + mc*wp, T)
                for q in 1:n, p in 1:n
                    M[p, q] += A[a, po[p], c, r, i]*d*conj(A[b, po[q], c, r, i])
                end
            end
            if k >= 0
                push!(get!(vsamples, k, Tuple{Float64,Matrix{Complex{Float64}}}[]), (nu, M))
            else
                push!(get!(vsamples, -k, Tuple{Float64,Matrix{Complex{Float64}}}[]), (nu + k*wp, Matrix(M')))
            end
        end
        vscale = max(1.0, maximum(abs, Ca))
        vproviders = harmonictables(vsamples, ks, n, wp, step, noise.interpolation, noise.atol*vscale,
            (k, nu, d) -> "the stated covariance is not Hermitian over the modes: the samples of the harmonic $(k) at $(nu) rad/s from two mode pairs differ by $(d) against the largest entry $(vscale).")
        NoiseCovariance(vproviders, noise.interpolation, :zero, noise.atol, noise.completed, noise.padding)
    else
        noise
    end
    checkdcmodel(dcmodel, n, stated, Float64(atol))
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
    declared = declaredtolerance(built)
    sq = Vector{Float64}(undef, nm)
    for i in eachindex(wv)
        sq .= wv[i] .+ modes .* wp
        v = pumpedviolation(built, sq, pumpedharmonics(built, sq))
        v <= declared || throw(ArgumentError(lazy"the block does not meet what it declares: at the frequency index $(i), over the modes of the solve, the violation of its losslessness or of the commutation relations of its stated covariance is $(v) of the square of its largest entry, against the $(declared) of atol and its noise model's. Declare its noise with noise = NoiseCovariance(linearized.Cnoise) from a solve with returnCnoise = true, leave a pump port out with the ports keyword, or raise atol to admit the discrepancy between the solves the tables hold at one frequency."))
    end
    return built
end

# The unconverted response of a pumped block is that of a real system,
# `H_0(-nu) = conj(H_0(nu))`, which the block evaluates as given at signed
# frequencies: a constant meets it only if it is real, and a table
# holding frequencies of both signs is held to it at each of its knots
# whose negative it also holds, to `atol` of its largest entry, or of one
# for entries below one. A callable is taken as given.
function checkrealresponse(p::AbstractMatrixProvider, atol::Float64)
    if p isa ConstantMatrixProvider
        d = maximum(x -> abs(imag(x)), p.A; init = 0.0)
        d <= atol*max(1.0, maximum(abs, p.A; init = 0.0)) || throw(ArgumentError(lazy"a constant unconverted response H_0 is that of a real system, H_0(-nu) = conj(H_0(nu)), only if it is real, and this one has an imaginary part of $(d); state a response which depends on frequency as a table or a callable."))
        return nothing
    end
    nus = filter(nu -> holdsdata(p, nu) && holdsdata(p, -nu), unique!(abs.(tableknots(p))))
    isempty(nus) && return nothing
    n = providersize(p)
    A = Array{Complex{Float64},3}(undef, n, n, length(nus))
    B = similar(A)
    evaluateprovider!(A, p, nus)
    evaluateprovider!(B, p, -nus)
    scale = max(1.0, maximum(abs, A))
    for i in eachindex(nus)
        d = maximum(abs, view(B, :, :, i) .- conj.(view(A, :, :, i)))
        d <= atol*scale || throw(ArgumentError(lazy"the unconverted response H_0 is that of a real system, H_0(-nu) = conj(H_0(nu)), but at $(nus[i]) rad/s the table differs from that by $(d) against its largest entry $(scale); raise atol to admit a discrepancy of the data."))
    end
    return nothing
end

# The samples `(nu, M)` of each harmonic in `ks`, from every mode pair of
# a solve, as a table band by band (see piecewisetable): the samples at
# one frequency from several mode pairs are one sample, which one
# periodic steady state gives alike, so they must agree to `tol`, and
# `mismatch(k, nu, d)` says where they do not.
function harmonictables(samples, ks, n::Int, wp, step, interpolation::Symbol, tol, mismatch)
    providers = AbstractMatrixProvider[]
    for k in ks
        list = sort!(samples[k]; by = first)
        nus = Float64[]
        mats = Matrix{Complex{Float64}}[]
        for (nu, M) in list
            if !isempty(nus) && samefrequency(nu, nus[end], wp)
                d = maximum(abs, M .- mats[end])
                d <= tol || throw(ArgumentError(mismatch(k, nu, d)))
                continue
            end
            push!(nus, nu)
            push!(mats, M)
        end
        values = Array{Complex{Float64},3}(undef, n, n, length(nus))
        for j in eachindex(nus)
            values[:, :, j] .= mats[j]
        end
        push!(providers, piecewisetable(nus, values; interpolation = interpolation, step = step))
    end
    return providers
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
        isharmonic(d, k, block.wp) || continue
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
    checkdestsize(dest, providersize(p), length(ws))
    addrational!(dest, p.cosine, ws, 1; add = false)
    return addrational!(dest, p.sine, ws, im)
end

# whether every harmonic of a pumped block has a realization in time
realizedintime(block::LinearizedScattering) = all(p -> p isa RationalScatteringProvider ||
    p isa ModulatedRationalProvider, block.providers)
