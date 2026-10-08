# A small expression type for parameterized component values.
#
# A component value may be written as an expression in named parameters
# (`Lj/2`, `1/(im*w*C)`) which is resolved to a number once the parameters
# are defined, possibly in two steps: the circuit definitions first, and
# the mode frequency later. `CircuitValue` is the tree such an expression is
# stored as: what the parameters of `@params` build, and what the
# deprecated netlist file reader of circuit/legacy.jl parses an expression
# to, so the numeric path of the package never needs Symbolics itself. A
# Symbolics `Num` is not lowered to one; the extension substitutes the
# definitions into it.
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
# A frequency dependent leaf: an opaque callable of a positive frequency,
# evaluated by `evalproviders`. The frequency law itself is
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

"""
    @params name1 name2 ...

Declare circuit parameters: each name becomes a `CircuitValues.Parameter`
bound to that name, which component values may be written in terms of
(`Lj/2`, `1/(Cc*w)`), and which `circuitdefs` supplies a number for at solve
time, keyed by the parameter or by its symbol. Returns the tuple of the
parameters. This is the dependency free counterpart of Symbolics'
`@variables`.
"""
macro params(names...)
    ex=Expr(:block)
    for n in names; push!(ex.args,:($(esc(n))=$(Parameter)($(QuoteNode(n))))); end
    push!(ex.args,Expr(:tuple,(esc(n) for n in names)...)); ex
end

# Promotion. The node types are parametric in their operator, so the
# promotion of `Binary{typeof(*)}` with `Binary{typeof(+)}` would be the
# UnionAll `Binary`, which is not a `DataType`. Every promotion is collapsed
# to the abstract `CircuitValue` instead, which is one, so that a group of
# values (`grouptype` in circuit/bind.jl) has a `DataType` element type.
Base.promote_rule(::Type{<:CircuitValue}, ::Type{<:CircuitValue}) = CircuitValue
Base.promote_rule(::Type{<:CircuitValue}, ::Type{<:Number}) = CircuitValue
Base.convert(::Type{CircuitValue}, x::Number) = Constant(x)
Base.convert(::Type{CircuitValue}, x::CircuitValue) = x

# Scalar semantics. A `CircuitValue` stands for one component value, so it
# must broadcast as a scalar the way a number does. Without this a
# broadcast such as `substitutefreq.(vvn, w)` would try to iterate the
# value.
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
# every binary operation of the closed set, `+ - * / ^`, is written infix
Base.show(io::IO, b::Binary) = print(io, "(", b.a, " ", nameof(b.f), " ", b.b, ")")

# the set of parameter names an expression depends on
parameters(e)=(s=Set{Symbol}(); _p!(s,e); s)
_p!(s,p::Parameter)=(push!(s,p.name);s); _p!(s,::Constant)=s
_p!(s,::Provider)=s
_p!(s,u::Unary)=_p!(s,u.a); _p!(s,b::Binary)=(_p!(s,b.a);_p!(s,b.b);s)

# whether an expression depends on the frequency, through a `Provider` leaf
hasprovider(::Provider) = true
hasprovider(u::Unary) = hasprovider(u.a)
hasprovider(b::Binary) = hasprovider(b.a) || hasprovider(b.b)
hasprovider(_) = false

#     substituteparams(expr, d)
#
# Replace the parameters named in the dictionary `d` by their values, a
# number or a value of this module, and leave the rest free. Because the
# constructors fold constants, an expression whose parameters are all
# defined collapses to a `Constant`, while one which still depends on an
# undefined parameter comes back as a tree. A frequency dependent leaf is
# a `Provider` rather than a parameter and passes through, for `freqsubst`
# in harmonics/sparse.jl to resolve once per mode.
substituteparams(c::Constant, d) = c
substituteparams(p::Provider, d) = p
substituteparams(q::Parameter, d) =
    haskey(d, q.name) ? tocv(d[q.name]) : q
substituteparams(u::Unary, d) = mk(u.f, substituteparams(u.a, d))
substituteparams(b::Binary, d) =
    mk(b.f, substituteparams(b.a, d), substituteparams(b.b, d))

#     evalproviders(expr, w)
#
# Replace every `Provider` leaf by its value at the frequency `w`.
# The constant folding constructors collapse the result, so an expression
# whose only unresolved leaves were providers comes back a `Constant`.
evalproviders(c::Constant, w) = c
evalproviders(q::Parameter, w) = q
evalproviders(p::Provider, w) = Constant(ComplexF64(p.f(w)))
evalproviders(u::Unary, w) = mk(u.f, evalproviders(u.a, w))
evalproviders(b::Binary, w) = mk(b.f, evalproviders(b.a, w), evalproviders(b.b, w))

#     derivative(expr, name)
#
# The derivative of an expression with respect to the parameter `name`, as
# an expression, for the design sensitivities. The constructors fold the
# constants, so the derivative of an expression which does not depend on
# the parameter collapses to `Constant(0)`. A frequency dependent leaf is
# a closure of the frequency alone and depends on no parameter. A design
# parameter is real, so `conj`, `real` and `imag` commute with the
# derivative.
derivative(p::Parameter, name::Symbol) = Constant(p.name === name ? 1 : 0)
derivative(::Constant, name::Symbol) = Constant(0)
derivative(::Provider, name::Symbol) = Constant(0)
function derivative(u::Unary, name::Symbol)
    a = u.a
    da = derivative(a, name)
    f = u.f
    f === (-) && return mk(-, da)
    f === inv && return mk(-, mk(/, da, mk(*, a, a)))
    f === sqrt && return mk(/, da, mk(*, Constant(2), mk(sqrt, a)))
    f === exp && return mk(*, mk(exp, a), da)
    f === log && return mk(/, da, a)
    (f === conj || f === real || f === imag) && return mk(f, da)
    throw(ArgumentError(lazy"no derivative for $(f)."))
end
function derivative(b::Binary, name::Symbol)
    a, c = b.a, b.b
    da, dc = derivative(a, name), derivative(c, name)
    f = b.f
    f === (+) && return mk(+, da, dc)
    f === (-) && return mk(-, da, dc)
    f === (*) && return mk(+, mk(*, da, c), mk(*, a, dc))
    f === (/) && return mk(/, mk(-, mk(*, da, c), mk(*, a, dc)), mk(*, c, c))
    if f === (^)
        # c a^(c-1) da, and a^c log(a) dc when the exponent depends on the
        # parameter
        t = mk(*, mk(*, c, mk(^, a, mk(-, c, Constant(1)))), da)
        _z(dc) && return t
        return mk(+, t, mk(*, mk(*, mk(^, a, c), mk(log, a)), dc))
    end
    throw(ArgumentError(lazy"no derivative for $(f)."))
end

end # module CircuitValues

using .CircuitValues

"""
    CircuitValue

The package's own expression type for a component value written in terms
of parameters: a parameter, a constant, a unary or binary operation on
them, or a frequency dependent provider. Built by arithmetic on the
parameters of [`@params`](@ref) and resolved to a number by
[`valuetonumber`](@ref) with the definitions in `circuitdefs`.
"""
const CircuitValue = CircuitValues.CircuitValue

"""
    FrequencyDependent(f)

A frequency dependent component value. `f` is called with a nonnegative
frequency in radians per second, zero at a direct current mode, and returns
the component value at that frequency. The function may be arbitrary Julia:
a closure over other parameters, a special function, an interpolation of
tabulated data.

```julia
R0 = 50.0; wc = 2*pi*10e9
Resistor(FrequencyDependent(w -> R0*(1 + im*w/wc)))
```

The value describes the element at positive frequencies. A mode of
negative frequency, such as an idler, takes the complex conjugate of the
value at the magnitude of its frequency, `conj(f(abs(w)))`, which is the
value of a real element there, as the data of a scattering block is
extended by [`ConjugateSymmetry`](@ref). `FrequencyDependent(identity)` is
the magnitude of the frequency, which may be written into an expression
like any other value.

The value is a [`CircuitValue`](@ref) whose leaf is the closure, so it
combines with numbers and other component values using the operators of
that type, `+ - * / ^` and the unary `- inv sqrt exp log conj real imag`.
For anything richer, put the whole expression inside the closure.
"""
FrequencyDependent(f) = CircuitValues.Provider(f)

# === resolving a written value to a number ===
#
# A component value is written as a number, as a symbol or string naming a
# parameter, as a `CircuitValue` expression, or as a symbolic value of the
# Symbolics extension. A parameter is defined by a number, by a value
# written in other parameters, to any depth, or by a frequency dependent
# value. The definitions are resolved once, each name to what its
# definition comes to (`resolvedefinitions`), and every value substitutes
# that one table, whether it names a parameter or is an expression in
# several. The design derivatives follow the same definitions by the chain
# rule (`DesignChain`).

"""
    componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)

Resolve each component value in `componentvalues` at the definitions
`circuitdefs`, resolved once for the whole table
([`resolvedefinitions`](@ref)), and return the results as a `Vector{Any}`:
the table mixes real and
complex values (a port's entry is its reference impedance), symbolic values
and frequency dependent providers,
and the groups the assembly reads are typed when they are gathered from it
(see `grouptype`), so the table itself has one type for every circuit.

# Examples
```jldoctest
julia> JosephsonCircuits.componentvaluestonumber([:Lj1,:Lj2],Dict(:Lj1=>1e-12,:Lj2=>2e-12))
2-element Vector{Any}:
 1.0e-12
 2.0e-12

julia> JosephsonCircuits.@params Lj1 Lj2;JosephsonCircuits.componentvaluestonumber([Lj1,Lj1+Lj2],Dict(Lj1=>1e-12,Lj2=>2e-12))
2-element Vector{Any}:
 1.0e-12
 3.0e-12
```
"""
function componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)
    # A comprehension over a vector of known length preallocates its result,
    # where `map` over `zip(values, Iterators.repeated(dict))` widens the
    # result element by element. A table of floating point numbers is its
    # own resolution; for any other (one holding a name, an expression, a
    # symbolic value) the definitions are resolved once for the whole table
    # rather than once per value.
    all(v -> v isa AbstractFloat || v isa Complex, componentvalues) &&
        return Any[valuetonumber(value,circuitdefs) for value in componentvalues]
    resolved = resolvedefinitions(circuitdefs)
    return Any[resolvevalue(value, resolved) for value in componentvalues]
end

"""
    ResolvedDefinitions

The definitions of the parameters, each resolved once by
[`resolvedefinitions`](@ref): `values` holds what each name comes to, a
number as it was given or what a definition written in other parameters
resolves to, and `expressions` the same as a [`CircuitValue`](@ref)
wherever it is a number or one, for an expression to substitute. Mutable
so that it is a heap object: each value of a table reaches its method by
a dynamic call, which would copy an immutable one to the heap.
"""
mutable struct ResolvedDefinitions
    values::Dict{Symbol,Any}
    expressions::Dict{Symbol,CircuitValue}
end

"""
    resolvedefinitions(circuitdefs)

Resolve the definitions of the parameters once, for every value written in
them: a definition which is a number is its own value, and one written in
other parameters, as an expression or a symbolic value, has their values
substituted, resolved first, so a parameter may be defined in terms of
others to any depth (`:Lj => Lj0/2` beside `:Lj0 => 2e-9`). A frequency
dependent leaf stays for [`substitutefreq`](@ref) to evaluate, a parameter
the definitions do not give stays free, for the check of the values to
name, and any other definition is its own value. A parameter defined in
terms of itself, directly or through others, is refused. Returns a
[`ResolvedDefinitions`](@ref); the keys may be anything
[`definitionsbyname`](@ref) reads.
"""
function resolvedefinitions(circuitdefs)
    byname = definitionsbyname(circuitdefs)
    resolved = ResolvedDefinitions(Dict{Symbol,Any}(),
        Dict{Symbol,CircuitValue}())
    sizehint!(resolved.values, length(byname))
    sizehint!(resolved.expressions, length(byname))
    # the names whose definitions are being resolved, innermost last
    pending = Symbol[]
    for name in keys(byname)
        resolvename!(resolved, pending, byname, name)
    end
    return resolved
end

# Resolve the definition of `name`, those of the parameters it is written
# in first.
function resolvename!(resolved::ResolvedDefinitions, pending::Vector{Symbol},
        byname::Dict{Symbol,Any}, name::Symbol)
    haskey(resolved.values, name) && return nothing
    v = byname[name]
    if checkissymbolic(v)
        k = findfirst(==(name), pending)
        isnothing(k) || throwcircular(vcat(pending[k:end], name))
        push!(pending, name)
        for n in valuenames(v)
            haskey(byname, n) && resolvename!(resolved, pending, byname, n)
        end
        pop!(pending)
        v = resolvevalue(v, resolved)
    end
    resolved.values[name] = v
    if v isa CircuitValue
        resolved.expressions[name] = v
    elseif isnumeric(v)
        resolved.expressions[name] = CircuitValues.Constant(ComplexF64(v))
    end
    return nothing
end

function throwcircular(chain::Vector{Symbol})
    names = join(chain, " -> ")
    throw(ArgumentError(lazy"The parameter $(first(chain)) is defined in terms of itself: $(names)."))
end

# whether a value is a number rather than one which still depends on a
# parameter or on the frequency
isnumeric(v) = v isa Number && !checkissymbolic(v)

# the names of the parameters a value is written in, as it is written: the
# name a symbol or a string gives, the parameters of an expression, and
# the variables of a symbolic value of the Symbolics extension
valuenames(v::Union{Symbol,AbstractString}) = (definitionname(v),)
valuenames(v::CircuitValue) = CircuitValues.parameters(v)
valuenames(v) = (definitionname(x) for x in circuitvariables(v))

"""
    resolvevalue(value, resolved::ResolvedDefinitions)

A value at the resolved definitions `resolved`: a name, written as a
symbol or a string, is what its definition came to, and an expression
substitutes the values of its parameters, so a name means the same
whether it is a whole value or part of one. A parameter the definitions do
not give stays free. Anything else is resolved by its
[`valuetonumber`](@ref) method; the Symbolics extension adds the method of
the four argument form, `resolvevalue(value, byname, resolved,
circuitdefs)` with `byname` the resolved values, for its own values.
"""
resolvevalue(value, resolved::ResolvedDefinitions) =
    resolvevalue(value, resolved.values, resolved, resolved)
function resolvevalue(value::Union{Symbol,AbstractString}, byname, resolved,
        circuitdefs)
    name = definitionname(value)
    return haskey(byname, name) ? byname[name] : CircuitValues.Parameter(name)
end
resolvevalue(value::CircuitValue, byname, resolved, circuitdefs) =
    valuetonumber(value, resolved)
resolvevalue(value, byname, resolved, circuitdefs) =
    valuetonumber(value, circuitdefs)

# A definition the Symbolics extension does not substitute into a
# symbolic value, one which is not a number, is the value of a bare
# variable it defines, resolved as an expression is.
resolvedefinition(v::CircuitValue, d) = valuetonumber(v, d)
resolvedefinition(v, d) = v

"""
    valuetonumber(value::Symbol,circuitdefs)

A symbol names a parameter; return the value `circuitdefs` defines it as,
under whichever key names it: its symbol, its string or its parameter
object (see [`definitionsbyname`](@ref)). A definition written in other
parameters is resolved with their definitions, to any depth (see
[`resolvedefinitions`](@ref)). A name `circuitdefs` does not define comes
back as the parameter, which the check of the values reports naming its
component.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(:Lj1,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0e-12

julia> JosephsonCircuits.valuetonumber(:Lj1,Dict("Lj1"=>1e-12))
1.0e-12
```
"""
valuetonumber(value::Symbol,circuitdefs) =
    resolvevalue(value, resolvedefinitions(circuitdefs))

"""
    valuetonumber(value::String,circuitdefs)

A string names a parameter, as a symbol does; see
[`valuetonumber(::Symbol, ::Any)`](@ref).

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber("Lj1",Dict("Lj1"=>1e-12,"Lj2"=>2e-12))
1.0e-12
```
"""
valuetonumber(value::String,circuitdefs) =
    valuetonumber(Symbol(value), circuitdefs)

"""
    valuetonumber(value::CircuitValue, circuitdefs)

Substitute the definitions in `circuitdefs`, which may be keyed by `Symbol`,
by `String`, or by the parameter objects themselves, into a parameterized
component value, each resolved first (see [`resolvedefinitions`](@ref)),
so a parameter defined in terms of others means here what it means as a
whole value.

A fully defined value comes back as a plain number, real when its imaginary
part is zero. A value which still depends on an undefined parameter comes
back as an expression, which the check of the values reports naming its
component.
"""
valuetonumber(value::CircuitValue, circuitdefs) =
    valuetonumber(value, resolvedefinitions(circuitdefs))

# the resolved form: the definitions resolved once, for every value
function valuetonumber(value::CircuitValue, resolved::ResolvedDefinitions)
    out = CircuitValues.substituteparams(value, resolved.expressions)
    out isa CircuitValues.Constant || return out
    return iszero(imag(out.val)) ? real(out.val) : out.val
end

"""
    normalizedefinitions(byname::Dict{Symbol,Any})

The definitions a parameterized value substitutes, from those keyed by
parameter name (see [`definitionsbyname`](@ref)), as a dictionary of
complex numbers. An entry whose value is not a number (a symbolic
expression included) cannot be substituted into an expression and is left
out, so that a definition which serves another purpose (a component whose
value is itself a key, for instance) does not stop every parameterized
value from resolving; a parameter left undefined comes back unresolved
from [`valuetonumber`](@ref).
"""
function normalizedefinitions(byname::Dict{Symbol,Any})
    d = Dict{Symbol,ComplexF64}()
    for (name, v) in byname
        (v isa Number && !checkissymbolic(v)) || continue
        d[name] = ComplexF64(v)
    end
    return d
end

"""
    definitionsbyname(circuitdefs)

The definitions keyed by the name of the parameter each key names (see
`definitionname`), with their values as given; a key which names
no parameter is left out. A parameter may be defined under its symbol,
its string or its parameter object, and one defined under two of them with
different values is refused rather than resolved by the order of the
dictionary. A table already keyed by name is returned as it is, and
[`ResolvedDefinitions`](@ref) as the values they resolved to.
"""
function definitionsbyname(circuitdefs)
    byname = Dict{Symbol,Any}()
    # a dictionary, or any collection of pairs
    for (k, v) in circuitdefs
        name = definitionname(k)
        isnothing(name) && continue
        if haskey(byname, name) && !isequal(byname[name], v)
            throw(ArgumentError(lazy"The parameter $(name) is defined twice, as $(byname[name]) and as $(v)."))
        end
        byname[name] = v
    end
    return byname
end
definitionsbyname(byname::Dict{Symbol,Any}) = byname
definitionsbyname(resolved::ResolvedDefinitions) = resolved.values

"""
    definitiontable(circuitdefs)

The component definitions as a `Dict{Any,Any}`, whatever the key and value
types of the dictionary given.
"""
definitiontable(circuitdefs::Dict{Any,Any}) = circuitdefs
definitiontable(circuitdefs::AbstractDict) = Dict{Any,Any}(circuitdefs)

# The name of the design parameter a definition key names, and nothing for
# a key which names no parameter: a parameter may be defined under its
# parameter object, under its symbol or under its string, and the solver
# cache and the design sensitivities both ask which parameter a key is for.
# The Symbolics extension adds the method for a `Num`.
definitionname(k::CircuitValues.Parameter) = k.name
definitionname(k::Symbol) = k
definitionname(k::AbstractString) = Symbol(k)
definitionname(k) = nothing

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

# === resolving and checking the values of a circuit ===

"""
    symbolicindices(A)

Return the indices in `A.nzval` where the elements of the matrix `A` are
symbolic variables.

# Examples
```jldoctest
julia> A = JosephsonCircuits.SparseArrays.sparse([1,2,1], [1,2,2], [1,1.0,2+3im]);JosephsonCircuits.symbolicindices(A)
Int64[]
```
"""
function symbolicindices(A)

    indices = Vector{Int}(undef,0)

    for (i,j) in enumerate(A)
        if checkissymbolic(j)
            push!(indices,i)
        end
    end
    return indices

end

function symbolicindices(A::SparseMatrixCSC)
    return symbolicindices(A.nzval)
end

"""
    checkissymbolic(a)

Check if `a` is a symbolic variable. Define a function to do this because
the test depends on which representation the value came from: the core
answer for `CircuitValue`, which a frequency dependent closure is a leaf
of, and the Symbolics extension adds the methods for its own wrappers.

# Examples
```jldoctest
julia> JosephsonCircuits.@params w;JosephsonCircuits.checkissymbolic(w)
true

julia> JosephsonCircuits.checkissymbolic(1.0)
false
```
"""
function checkissymbolic(a)
    return a isa CircuitValue
end

"""
    circuitvariables(a)

The free parameters of a component value. Returns an empty collection for
a numeric value and the `CircuitValues.Parameter`s of a
[`CircuitValue`](@ref) (not `Symbol`s). The Symbolics extension adds a
method for `Num`.
"""
circuitvariables(a) = Symbol[]

"""
    substitutefreq(value, w)

Resolve a component value at the mode frequency `w`: the identity for a
plain number, and for a [`CircuitValue`](@ref) the evaluation of its
frequency dependent leaves (see [`FrequencyDependent`](@ref)) at the
magnitude of `w`, followed by the constant folding of the expression around
them. A value states the element at positive frequencies, and a mode of
negative frequency takes its conjugate where the value is placed (see
[`modevalue`](@ref)), which is the value a real element has there. A value
which does not resolve to a number is returned as it is, for the caller to
diagnose. The Symbolics extension adds the `Num` method.
"""
substitutefreq(value, w) = value
function substitutefreq(value::CircuitValue, w)
    v = CircuitValues.evalproviders(value, abs(w))
    return v isa CircuitValues.Constant ?
        (iszero(imag(v.val)) ? real(v.val) : v.val) : v
end

# the parameters as `Parameter` objects rather than bare symbols, so that
# they print and compare as the user wrote them
circuitvariables(a::CircuitValue) =
    [CircuitValues.Parameter(n) for n in sort!(collect(CircuitValues.parameters(a)))]

"""
    checkcomponentvaluesdefined(componentnames::Vector, vvn::Vector)

Check that no circuit component value still depends on a free parameter.
One which does indicates a parameter which was not assigned a numerical
value in the circuit definitions dictionary `circuitdefs`, and an
informative `ArgumentError` is thrown naming the components and the
undefined parameters. A frequency dependent value carries a closure of
the frequency rather than a parameter and is resolved per frequency by
[`freqsubst`](@ref), so it passes. Called by [`hbnlsolve`](@ref) and
[`hblinsolve`](@ref) before any computation, so a forgotten entry in
`circuitdefs` fails immediately with the actual cause instead of a
downstream error about a symbolic value in a matrix.

# Examples
```jldoctest
julia> JosephsonCircuits.checkcomponentvaluesdefined(["P1","R1"], Any[1, FrequencyDependent(w -> 1/(w*50.0))])

julia> JosephsonCircuits.@params R2;try JosephsonCircuits.checkcomponentvaluesdefined(["P1","R1"], Any[1, R2]) catch e; occursin("R1 has the value", sprint(showerror, e)) end
true
```
"""
function checkcomponentvaluesdefined(componentnames::Vector, vvn::Vector)
    messages = String[]
    for i in eachindex(vvn)
        if checkissymbolic(vvn[i])
            undefined = circuitvariables(vvn[i])
            if !isempty(undefined)
                push!(messages, string("The component ", componentnames[i],
                    " has the value ", vvn[i],
                    ", which contains the symbolic variables [",
                    join(string.(undefined), ", "), "] that were not "*
                    "assigned numerical values."))
            end
        end
    end
    if !isempty(messages)
        throw(ArgumentError(join(messages, " ")*" Add the missing "*
            "variables to the circuit definitions dictionary "*
            "circuitdefs. If a variable represents the frequency of a "*
            "frequency dependent component, write the value as a "*
            "FrequencyDependent closure of the frequency instead."))
    end
    return nothing
end

# the `Num` method is defined in the Symbolics extension
