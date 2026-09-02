# A small expression type for parameterized component values.
#
# A component value may be written as an expression in named parameters
# (`Lj/2`, `1/(im*w*C)`) which is resolved to a number once the parameters
# are defined, possibly in two steps: the circuit definitions first, and
# the mode frequency later. `CircuitValue` is the tree such an expression is
# stored as. It is what the Symbolics extension lowers a `Num` to and what
# a parameterized netlist file expression parses to, so the numeric path of
# the package never needs Symbolics itself.
#
# The operator set is closed and small on purpose: `+ - * / ^` and the
# unary `- inv sqrt exp log conj real imag`. Anything richer belongs inside
# a `FrequencyDependent` closure, where all of Julia is available.
module CircuitValues
export @params

abstract type CircuitValue end
struct Parameter <: CircuitValue; name::Symbol; end
struct Constant  <: CircuitValue; val::ComplexF64; end
struct Unary{F}  <: CircuitValue; f::F; a::CircuitValue; end
struct Binary{F} <: CircuitValue; f::F; a::CircuitValue; b::CircuitValue; end
# A frequency dependent leaf: an opaque callable of the signed mode
# frequency, evaluated by `evalproviders`. The frequency law itself is
# arbitrary Julia inside the closure; the operators of this module only
# combine whole component values with each other.
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

# lift a number to a leaf; a value is returned unchanged
tocv(x::CircuitValue)=x; tocv(x::Number)=Constant(x)

# Constructors which fold constants as the tree is built: `0*x` becomes
# `0`, `x + 0` becomes `x`, and an operation on two constants is evaluated.
# This keeps an expression the size the user wrote it, even after the stamp
# arithmetic has combined it with many zeros and ones, and it is what lets
# `substituteparams` collapse a fully defined expression to a `Constant`.
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

# Promotion. The node types are parametric in their operator, so the
# promotion of `Binary{typeof(*)}` with `Binary{typeof(+)}` would be the
# UnionAll `Binary`, which is not a `DataType`. `calcvaluetype` in
# capindmat.jl keys its type table by `DataType`, so every promotion is
# collapsed to the abstract `CircuitValue` instead, which is one.
Base.promote_rule(::Type{<:CircuitValue}, ::Type{<:CircuitValue}) = CircuitValue
Base.promote_rule(::Type{<:CircuitValue}, ::Type{<:Number}) = CircuitValue
Base.convert(::Type{CircuitValue}, x::Number) = Constant(x)
Base.convert(::Type{CircuitValue}, x::CircuitValue) = x

# Scalar semantics. A `CircuitValue` stands for one component value, so it
# must broadcast as a scalar the way a number does. Without this a
# broadcast such as `substitutefreq.(vvn, symfreqvar, w)` would try to
# iterate the frequency variable.
Base.length(::CircuitValue) = 1
Base.size(::CircuitValue) = ()
Base.ndims(::Type{<:CircuitValue}) = 0
Base.iterate(x::CircuitValue) = (x, nothing)
Base.iterate(::CircuitValue, ::Nothing) = nothing
Base.broadcastable(x::CircuitValue) = Ref(x)
Base.isequal(a::CircuitValue, b::CircuitValue) = a == b

# Printing as an expression rather than as nested structs. Component values
# appear verbatim in error messages about undefined parameters, where the
# default `show` of the tree would be unreadable.
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

# the set of parameter names an expression depends on
parameters(e)=(s=Set{Symbol}(); _p!(s,e); s)
_p!(s,p::Parameter)=(push!(s,p.name);s); _p!(s,::Constant)=s
_p!(s,::Provider)=s
_p!(s,u::Unary)=_p!(s,u.a); _p!(s,b::Binary)=(_p!(s,b.a);_p!(s,b.b);s)

#     substituteparams(expr, d)
#
# Replace the parameters named in the dictionary `d` by their values and
# leave the rest free. Because the constructors fold constants, an
# expression whose parameters are all defined collapses to a `Constant`,
# while one which still depends on an undefined parameter (typically the
# symbolic frequency variable, which `freqsubst` in matutils.jl resolves
# once per mode) comes back as a tree.
substituteparams(c::Constant, d) = c
substituteparams(p::Provider, d) = p
substituteparams(q::Parameter, d) =
    haskey(d, q.name) ? Constant(d[q.name]) : q
substituteparams(u::Unary, d) = mk(u.f, substituteparams(u.a, d))
substituteparams(b::Binary, d) =
    mk(b.f, substituteparams(b.a, d), substituteparams(b.b, d))

# evaluate an expression to a number, with every parameter looked up in `d`
evaluate(c::Constant,d)=c.val; evaluate(q::Parameter,d)=d[q.name]
evaluate(::Provider,d)=error("a frequency dependent value cannot be evaluated without a frequency; it is resolved per mode by freqsubst")

#     evalproviders(expr, w)
#
# Replace every `Provider` leaf by its value at the signed frequency `w`.
# The constant folding constructors collapse the result, so an expression
# whose only unresolved leaves were providers comes back a `Constant`.
evalproviders(c::Constant, w) = c
evalproviders(q::Parameter, w) = q
evalproviders(p::Provider, w) = Constant(ComplexF64(p.f(w)))
evalproviders(u::Unary, w) = mk(u.f, evalproviders(u.a, w))
evalproviders(b::Binary, w) = mk(b.f, evalproviders(b.a, w), evalproviders(b.b, w))

evaluate(u::Unary,d)=u.f(evaluate(u.a,d)); evaluate(b::Binary,d)=b.f(evaluate(b.a,d),evaluate(b.b,d))

# === parsing a component value from an expression ===
#
# A value written in a netlist file arrives as a parsed Julia `Expr`. This
# converts it to a `CircuitValue`, accepting only the closed operator set
# above plus the constants `im` and `pi`, so no Symbolics parser is needed.
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

A frequency dependent component value. `f` is called with each signed mode
frequency in radians per second and returns the component value at that
frequency. The function may be arbitrary Julia: a closure over other
parameters, a special function, an interpolation of tabulated data.

```julia
R0 = 50.0; wc = 2*pi*10e9
("R1", "1", "0", FrequencyDependent(w -> R0*(1 + im*w/wc)))
```

`f` receives signed frequencies, as a symbolic frequency variable does, so
a law defined only for positive frequencies should apply its own conjugate
rule for negative ones inside the closure.

The value may be combined with numbers and other component values using
`+ - * / ^` and the unary `- inv sqrt exp log conj real imag`. For anything
richer, put the whole expression inside the closure.
"""
struct FrequencyDependent{F}
    f::F
end

# A `FrequencyDependent` lowers to a `Provider` leaf so that the expression
# arithmetic of `CircuitValues` applies to it.
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
# A component value is written as a number, as a symbol or string looked up
# in `circuitdefs`, as a `CircuitValue` expression, or as a callable of
# frequency. `valuetonumber` turns each into a number, or into an expression
# which still depends on the frequency.

"""
    componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)

Resolve each component value in `componentvalues` with [`valuetonumber`](@ref)
and return the results as a vector. The element type of the result is
inferred from the resolved values, so it is `Vector{Float64}` when every
value resolves to a real number; this is deliberately not type stable.

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
    # A comprehension over a vector of known length preallocates its result;
    # `map` over `zip(values, Iterators.repeated(dict))` cannot, and widens
    # the result element by element, which is far slower on a large circuit.
    return [valuetonumber(value,circuitdefs) for value in componentvalues]
end

"""
    valuetonumber(value::Symbol,circuitdefs)

A symbol is a key of `circuitdefs`; return the value stored under it.

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

A string is a key of `circuitdefs`; return the value stored under it.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber("Lj1",Dict("Lj1"=>1e-12,"Lj2"=>2e-12))
1.0e-12
```
"""
function valuetonumber(value::String,circuitdefs)
    return circuitdefs[value]
end

# The definitions arrive either as a dictionary or, from the per mode
# frequency substitution in `sparseaddconjsubst!` (matutils.jl), as a
# single `symfreqvar => w` pair. `Symbolics.substitute` accepts both, and so
# does this.
_definitionpairs(d::AbstractDict) = pairs(d)
_definitionpairs(d::Pair) = (d,)
_definitionpairs(d) = d

"""
    valuetonumber(value::FrequencyDependent, circuitdefs)

A frequency dependent value has no number to resolve to yet. It is lowered
to a `Provider` leaf, which passes through the definitions unchanged and is
evaluated at each mode frequency later by [`freqsubst`](@ref).
"""
function valuetonumber(value::FrequencyDependent, circuitdefs)
    return CircuitValues.Provider(value.f)
end

"""
    valuetonumber(value::CircuitValue, circuitdefs)

Substitute the definitions in `circuitdefs`, which may be keyed by `Symbol`
or by the parameter objects themselves, into a parameterized component
value.

A fully defined value comes back as a plain number, real when its imaginary
part is zero. A value which still depends on an undefined parameter comes
back as an expression; this is how a value depending on the symbolic
frequency variable reaches [`freqsubst`](@ref) with that variable free.
"""
function valuetonumber(value::CircuitValue, circuitdefs)
    d = Dict{Symbol,ComplexF64}()
    for (k,v) in _definitionpairs(circuitdefs)
        key = k isa CircuitValues.Parameter ? k.name : Symbol(k)
        d[key] = ComplexF64(v)
    end
    out = CircuitValues.substituteparams(value, d)
    out isa CircuitValues.Constant || return out
    return iszero(imag(out.val)) ? real(out.val) : out.val
end

# The methods for Symbolics `Num` and `BasicSymbolic` values are defined in
# ext/JosephsonCircuitsSymbolicsExt.jl.

"""
    valuetonumber(value, circuitdefs)

A number, or any other type not handled by a more specific method, is
returned unchanged.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(1.0,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0
```
"""
function valuetonumber(value, circuitdefs)
    return value
end
