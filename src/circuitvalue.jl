module CircuitValues
export @params

abstract type CircuitValue end
struct Parameter <: CircuitValue; name::Symbol; end
struct Constant  <: CircuitValue; val::ComplexF64; end
struct Unary{F}  <: CircuitValue; f::F; a::CircuitValue; end
struct Binary{F} <: CircuitValue; f::F; a::CircuitValue; b::CircuitValue; end
# A frequency-dependent component value: an opaque callable evaluated at
# each signed mode frequency by `evalproviders`. The user's frequency law
# is arbitrary Julia inside the closure; the closed operator set of this
# module only ever combines whole component values (series/parallel stamp
# arithmetic), which is bounded by the physics.
struct Provider  <: CircuitValue; f::Any; end

Base.:(==)(a::Parameter,b::Parameter)=a.name==b.name
Base.:(==)(a::Constant,b::Constant)=a.val==b.val
Base.:(==)(a::Unary,b::Unary)=a.f===b.f && a.a==b.a
Base.:(==)(a::Binary,b::Binary)=a.f===b.f && a.a==b.a && a.b==b.b
Base.:(==)(a::Provider,b::Provider)=a.f===b.f
Base.:(==)(::CircuitValue,::CircuitValue)=false
Base.hash(p::Parameter,h::UInt)=hash(p.name,hash(:P,h))
Base.hash(c::Constant,h::UInt)=hash(c.val,hash(:C,h))
Base.hash(u::Unary,h::UInt)=hash(u.a,hash(u.f,hash(:U,h)))
Base.hash(b::Binary,h::UInt)=hash(b.b,hash(b.a,hash(b.f,hash(:B,h))))
Base.hash(p::Provider,h::UInt)=hash(objectid(p.f),hash(:F,h))

tocv(x::CircuitValue)=x; tocv(x::Number)=Constant(x)

# Smart constructors. Constant folding is not cosmetic here: a component
# value assembled from a parameterized netlist accumulates a great deal of
# `0*x`, `1*x` and `x+0` on its way through the stamp arithmetic, and an
# expression which folds as it is built stays the size of the expression the
# user wrote. It also makes structurally equal values compare equal, which
# is what `substituteparams` and the value type machinery rely on.
_isc(x) = x isa Constant
_z(x) = _isc(x) && iszero(x.val)
_o(x) = _isc(x) && isone(x.val)

mk(::typeof(+),a,b) = _z(a) ? b : _z(b) ? a :
    (_isc(a)&&_isc(b)) ? Constant(a.val+b.val) : Binary(+,a,b)
mk(::typeof(-),a,b) = _z(b) ? a :
    (_isc(a)&&_isc(b)) ? Constant(a.val-b.val) : _z(a) ? mk(-,b) : Binary(-,a,b)
mk(::typeof(*),a,b) = (_z(a)||_z(b)) ? Constant(0) : _o(a) ? b : _o(b) ? a :
    (_isc(a)&&_isc(b)) ? Constant(a.val*b.val) : Binary(*,a,b)
mk(::typeof(/),a,b) = _z(a) ? Constant(0) : _o(b) ? a :
    (_isc(a)&&_isc(b)) ? Constant(a.val/b.val) : Binary(/,a,b)
mk(::typeof(^),a,b) = _o(b) ? a : _z(b) ? Constant(1) :
    (_isc(a)&&_isc(b)) ? Constant(a.val^b.val) : Binary(^,a,b)
mk(f,a) = _isc(a) ? Constant(f(a.val)) : Unary(f,a)
mk(::typeof(-),a) = _isc(a) ? Constant(-a.val) : Unary(-,a)

for op in (:+,:-,:*,:/,:^)
    @eval begin
        Base.$op(a::CircuitValue,b::CircuitValue)=mk($op,a,b)
        Base.$op(a::CircuitValue,b::Number)=mk($op,a,tocv(b))
        Base.$op(a::Number,b::CircuitValue)=mk($op,tocv(a),b)
    end
end
for op in (:-,:inv,:sqrt,:exp,:log,:conj,:real,:imag)
    @eval Base.$op(a::CircuitValue)=mk($op,a)
end

Base.zero(::Type{<:CircuitValue})=Constant(0); Base.one(::Type{<:CircuitValue})=Constant(1)

macro params(names...)
    ex=Expr(:block)
    for n in names; push!(ex.args,:($(esc(n))=$(Parameter)($(QuoteNode(n))))); end
    push!(ex.args,Expr(:tuple,(esc(n) for n in names)...)); ex
end

# Promotion. The node types are parametric in the operator, so
# `promote_type(Binary{typeof(*)}, Binary{typeof(+)})` is the UnionAll
# `Binary`, which is not a `DataType` and breaks `calcvaluetype`, whose
# type store is keyed by `DataType`. Collapsing every promotion to the
# abstract `CircuitValue` -- which is a `DataType` -- keeps the value type
# machinery working, the way a single concrete `Num` does for Symbolics.
Base.promote_rule(::Type{<:CircuitValue}, ::Type{<:CircuitValue}) = CircuitValue
Base.promote_rule(::Type{<:CircuitValue}, ::Type{<:Number}) = CircuitValue
Base.convert(::Type{CircuitValue}, x::Number) = Constant(x)
Base.convert(::Type{CircuitValue}, x::CircuitValue) = x

# Scalar semantics. A CircuitValue stands for one component value, so it
# must broadcast as a scalar the way a number does; without this,
# `substitutefreq.(vvn, symfreqvar, w)` tries to iterate the frequency
# variable.
Base.length(::CircuitValue) = 1
Base.size(::CircuitValue) = ()
Base.ndims(::Type{<:CircuitValue}) = 0
Base.iterate(x::CircuitValue) = (x, nothing)
Base.iterate(::CircuitValue, ::Nothing) = nothing
Base.broadcastable(x::CircuitValue) = Ref(x)
Base.isequal(a::CircuitValue, b::CircuitValue) = a == b

# Readable printing. Component values appear verbatim in error messages
# about undefined variables, so the default struct show turns a one line
# expression into an unreadable wall of nested type names.
Base.show(io::IO, p::Parameter) = print(io, p.name)
Base.show(io::IO, p::Provider) = print(io, "FrequencyDependent(", p.f, ")")
function Base.show(io::IO, c::Constant)
    v = c.val
    print(io, iszero(imag(v)) ? real(v) : v)
end
Base.show(io::IO, u::Unary) = print(io, nameof(u.f), "(", u.a, ")")
function Base.show(io::IO, b::Binary)
    op = nameof(b.f)
    if op in (:+, :-, :*, :/, :^)
        print(io, "(", b.a, " ", op, " ", b.b, ")")
    else
        print(io, op, "(", b.a, ", ", b.b, ")")
    end
end

parameters(e)=(s=Set{Symbol}(); _p!(s,e); s)
_p!(s,p::Parameter)=(push!(s,p.name);s); _p!(s,::Constant)=s
_p!(s,::Provider)=s
_p!(s,u::Unary)=_p!(s,u.a); _p!(s,b::Binary)=(_p!(s,b.a);_p!(s,b.b);s)

#     substituteparams(expr, d)
#
# Replace the parameters named in `d` by their values, leaving the rest free.
#
# Because the smart constructors fold constants, an expression whose leaves
# all resolve collapses to a `Constant`, and one which still depends on a
# free parameter -- the symbolic frequency variable, resolved later per mode
# by `freqsubst` -- comes back as a `CircuitValue`. That is what lets a
# frequency dependent component value survive `valuetonumber` and be
# evaluated once per mode afterwards.
substituteparams(c::Constant, d) = c
substituteparams(p::Provider, d) = p
substituteparams(q::Parameter, d) =
    haskey(d, q.name) ? Constant(d[q.name]) : q
substituteparams(u::Unary, d) = mk(u.f, substituteparams(u.a, d))
substituteparams(b::Binary, d) =
    mk(b.f, substituteparams(b.a, d), substituteparams(b.b, d))

evaluate(c::Constant,d)=c.val; evaluate(q::Parameter,d)=d[q.name]
evaluate(::Provider,d)=error("a frequency dependent value cannot be evaluated without a frequency; it is resolved per mode by freqsubst")

#     evalproviders(expr, w)
#
# Replace every frequency-dependent provider leaf by its value at the
# signed frequency `w`. The constant-folding constructors collapse the
# result, so an expression whose only unresolved leaves were providers
# comes back a `Constant`.
evalproviders(c::Constant, w) = c
evalproviders(q::Parameter, w) = q
evalproviders(p::Provider, w) = Constant(ComplexF64(p.f(w)))
evalproviders(u::Unary, w) = mk(u.f, evalproviders(u.a, w))
evalproviders(b::Binary, w) = mk(b.f, evalproviders(b.a, w), evalproviders(b.b, w))

evaluate(u::Unary,d)=u.f(evaluate(u.a,d)); evaluate(b::Binary,d)=b.f(evaluate(b.a,d),evaluate(b.b,d))

# ---------------------------------------------------------------------
# Parsing a component value from a string, replacing
# Symbolics.parse_expr_to_symbolic.
# ---------------------------------------------------------------------
module Parsing
import ..CircuitValues: CircuitValue, Constant, Parameter
const CONSTS = Dict{Symbol,Any}(:im => im, :pi => pi)
const BINOPS = Dict{Symbol,Function}(:+ => +, :- => -, :* => *, :/ => /, :^ => ^)
const UNOPS  = Dict{Symbol,Function}(:- => -, :inv => inv, :sqrt => sqrt,
                                     :exp => exp, :log => log, :conj => conj,
                                     :real => real, :imag => imag)
fromexpr(x::Number) = Constant(x)
function fromexpr(s::Symbol)
    haskey(CONSTS, s) && return Constant(CONSTS[s])
    return Parameter(s)
end
function fromexpr(e::Expr)
    e.head === :call || error("unsupported expression head $(e.head)")
    op = e.args[1]; args = map(fromexpr, e.args[2:end])
    length(args) == 1 && haskey(UNOPS, op) && return UNOPS[op](args[1])
    haskey(BINOPS, op) && return reduce(BINOPS[op], args)
    error("unsupported operator $(op) in a component value")
end
end

using .Parsing: fromexpr
const _fromexpr = fromexpr
end # module CircuitValues

using .CircuitValues
const CircuitValue = CircuitValues.CircuitValue

"""
    FrequencyDependent(f)

A frequency-dependent component value: `f` is called with each signed
mode frequency (radians per second) and returns the component value
there. The function may be arbitrary Julia -- closures over parameters,
special functions, table interpolation:

```julia
R0 = 50.0; wc = 2*pi*10e9
("R1", "1", "0", FrequencyDependent(w -> R0*(1 + im*w/wc)))
```

Negative-frequency convention: `f` receives signed frequencies, exactly
as a symbolic frequency variable did; a law defined only for positive
frequencies should implement its own conjugate rule inside the closure.

The value may be combined with numbers and other component values with
`+ - * / ^ inv sqrt exp log conj`; for anything richer, put the whole
expression inside the closure, where all of Julia is available.
"""
struct FrequencyDependent{F}
    f::F
end

# lower to a provider leaf so the whole expression machinery applies
CircuitValues.tocv(x::FrequencyDependent) = CircuitValues.Provider(x.f)
Base.convert(::Type{CircuitValues.CircuitValue}, x::FrequencyDependent) =
    CircuitValues.Provider(x.f)
Base.promote_rule(::Type{<:FrequencyDependent}, ::Type{<:Number}) = CircuitValue
Base.promote_rule(::Type{<:FrequencyDependent}, ::Type{<:CircuitValues.CircuitValue}) = CircuitValue
for op in (:+, :-, :*, :/, :^)
    @eval begin
        Base.$op(a::FrequencyDependent, b::FrequencyDependent) =
            CircuitValues.mk($op, CircuitValues.tocv(a), CircuitValues.tocv(b))
        Base.$op(a::FrequencyDependent, b::Union{Number,CircuitValues.CircuitValue}) =
            CircuitValues.mk($op, CircuitValues.tocv(a), CircuitValues.tocv(b))
        Base.$op(a::Union{Number,CircuitValues.CircuitValue}, b::FrequencyDependent) =
            CircuitValues.mk($op, CircuitValues.tocv(a), CircuitValues.tocv(b))
    end
end
for op in (:-, :inv, :sqrt, :exp, :log, :conj, :real, :imag)
    @eval Base.$op(a::FrequencyDependent) =
        CircuitValues.mk($op, CircuitValues.tocv(a))
end

# === resolving a written value to a number ===
#
# A component value is written as a number, a name to be looked up in
# `circuitdefs`, an expression in the types above, or a callable of
# frequency. This is where each becomes a number.
"""
    componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)

Convert the array of component values to numbers, if defined in `circuitdefs`.
This function is not type stable by design because we want the output array
to use a concrete type if all of the values are evaluated to numbers.

# Examples
```jldoctest
julia> JosephsonCircuits.componentvaluestonumber([:Lj1,:Lj2],Dict(:Lj1=>1e-12,:Lj2=>2e-12))
2-element Vector{Float64}:
 1.0e-12
 2.0e-12

julia> JosephsonCircuits.@params Lj1 Lj2;JosephsonCircuits.componentvaluestonumber([Lj1,Lj1+Lj2],Dict(Lj1=>1e-12,Lj2=>2e-12))
2-element Vector{Float64}:
 1.0e-12
 3.0e-12
```
"""
function componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)
    # A comprehension, not `map` over the values paired with a repeated
    # dictionary. The two agree on every value and on the element type they
    # produce -- which is the point of doing this at all, since the output
    # should be concretely typed when every value resolves to a number -- but
    # `map` over an iterator of unknown length cannot preallocate and widens
    # the result as it goes. Over 8192 values that is 1229 us against 7 us.
    return [valuetonumber(value,circuitdefs) for value in componentvalues]
end

"""
    valuetonumber(value::Symbol,circuitdefs)

If the component value is a symbol, assume it is a dictionary key.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(:Lj1,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0e-12
```
"""
function valuetonumber(value::Symbol,circuitdefs)
    return circuitdefs[value]
end

"""
    valuetonumber(value::String,circuitdefs)

If the component value is a string, assume it is a dictionary key.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber("Lj1",Dict("Lj1"=>1e-12,"Lj2"=>2e-12))
1.0e-12
```
"""
function valuetonumber(value::String,circuitdefs)
    return circuitdefs[value]
end

# the circuit definitions arrive either as a dictionary or, from the per
# mode frequency substitution in `sparseaddconjsubst!`, as a single
# `symfreqvar => w` pair. `Symbolics.substitute` accepts both, so this must
# too.
_definitionpairs(d::AbstractDict) = pairs(d)
_definitionpairs(d::Pair) = (d,)
_definitionpairs(d) = d

"""
    valuetonumber(value::FrequencyDependent, circuitdefs)

A frequency dependent value has no number to substitute for: it lowers to a
provider leaf which survives the definitions and is evaluated at each mode
frequency later, by [`freqsubst`](@ref).
"""
function valuetonumber(value::FrequencyDependent, circuitdefs)
    # a frequency dependent value lowers to a provider leaf and survives to
    # be resolved per mode by freqsubst
    return CircuitValues.Provider(value.f)
end

"""
    valuetonumber(value::CircuitValue, circuitdefs)

Substitute the definitions in `circuitdefs` into a parameterized component
value. `circuitdefs` may be keyed by `Symbol` or by the parameter itself.

A value which is fully defined comes back as a plain number, real when its
imaginary part is zero. One which is not stays an expression, which is how a
component value depending on the symbolic frequency variable reaches
[`freqsubst`](@ref) with that variable still free.
"""
function valuetonumber(value::CircuitValue, circuitdefs)
    d = Dict{Symbol,ComplexF64}()
    for (k,v) in _definitionpairs(circuitdefs)
        key = k isa CircuitValues.Parameter ? k.name : Symbol(k)
        d[key] = ComplexF64(v)
    end
    # substitute what is defined and leave the rest free. A component value
    # which still depends on the symbolic frequency variable comes back as a
    # CircuitValue and is resolved per mode later by `freqsubst`.
    out = CircuitValues.substituteparams(value, d)
    out isa CircuitValues.Constant || return out
    return iszero(imag(out.val)) ? real(out.val) : out.val
end

# the `Num` and `SymbolicT` methods live in the Symbolics extension, which
# documents them there: a docstring here would attach to whatever follows it

"""
    valuetonumber(value, circuitdefs)

If the component value `value` is a number (or a type we haven't considered,
return it as is.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(1.0,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0
```
"""
function valuetonumber(value, circuitdefs)
    return value
end
