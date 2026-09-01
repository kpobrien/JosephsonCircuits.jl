module CircuitValues
using RuntimeGeneratedFunctions
RuntimeGeneratedFunctions.init(@__MODULE__)
export @params

abstract type CircuitValue end
struct Parameter <: CircuitValue; name::Symbol; end
struct Constant  <: CircuitValue; val::ComplexF64; end
struct Hole      <: CircuitValue; i::Int; end
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
Base.:(==)(a::Hole,b::Hole)=a.i==b.i
Base.:(==)(a::Unary,b::Unary)=a.f===b.f && a.a==b.a
Base.:(==)(a::Binary,b::Binary)=a.f===b.f && a.a==b.a && a.b==b.b
Base.:(==)(a::Provider,b::Provider)=a.f===b.f
Base.:(==)(::CircuitValue,::CircuitValue)=false
Base.hash(p::Parameter,h::UInt)=hash(p.name,hash(:P,h))
Base.hash(c::Constant,h::UInt)=hash(c.val,hash(:C,h))
Base.hash(x::Hole,h::UInt)=hash(x.i,hash(:H,h))
Base.hash(u::Unary,h::UInt)=hash(u.a,hash(u.f,hash(:U,h)))
Base.hash(b::Binary,h::UInt)=hash(b.b,hash(b.a,hash(b.f,hash(:B,h))))
Base.hash(p::Provider,h::UInt)=hash(objectid(p.f),hash(:F,h))

tocv(x::CircuitValue)=x; tocv(x::Number)=Constant(x)

# Smart constructors. Constant folding is not cosmetic here: symbolic
# differentiation generates a great deal of `0*x` and `1*x`, and without
# folding the derivative expressions are larger than the primal and the
# shape dedup stops working. Folding at construction keeps both compact.
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

# ---------- symbolic differentiation ----------
# The operator set is tiny and closed, so the exact derivative of every
# component value with respect to every design parameter is available
# directly. It is a CircuitValue like any other, so it compiles and
# evaluates through the same machinery as the primal, and it is exact for
# complex component values, which is where forward mode dual numbers give
# trouble: a Complex{Dual} cannot be written into the real output array
# ForwardDiff.jacobian allocates.
derivative(x, ::Symbol) = Constant(0)
derivative(c::Constant, ::Symbol) = Constant(0)
derivative(h::Hole, ::Symbol) = Constant(0)
derivative(q::Parameter, nm::Symbol) = Constant(q.name === nm ? 1 : 0)
derivative(::Provider, ::Symbol) = Constant(0)
function derivative(u::Unary, nm::Symbol)
    da = derivative(u.a, nm)
    _z(da) && return Constant(0)
    u.f === (-)    && return mk(-, da)
    u.f === inv    && return mk(/, mk(-, da), mk(^, u.a, Constant(2)))
    u.f === sqrt   && return mk(/, da, mk(*, Constant(2), mk(sqrt, u.a)))
    u.f === exp    && return mk(*, mk(exp, u.a), da)
    u.f === log    && return mk(/, da, u.a)
    u.f === conj   && return mk(conj, da)
    # valid because the design parameters are real, so d/dp commutes
    # with taking real and imaginary parts
    u.f === real   && return mk(real, da)
    u.f === imag   && return mk(imag, da)
    error("no derivative rule for \$(u.f)")
end
function derivative(b::Binary, nm::Symbol)
    da = derivative(b.a, nm); db = derivative(b.b, nm)
    b.f === (+) && return mk(+, da, db)
    b.f === (-) && return mk(-, da, db)
    b.f === (*) && return mk(+, mk(*, da, b.b), mk(*, b.a, db))
    b.f === (/) && return mk(/, mk(-, mk(*, da, b.b), mk(*, b.a, db)),
                             mk(^, b.b, Constant(2)))
    if b.f === (^)
        _z(db) && return mk(*, mk(*, b.b, mk(^, b.a, mk(-, b.b, Constant(1)))), da)
        return mk(*, mk(^, b.a, b.b),
                  mk(+, mk(*, db, mk(log, b.a)), mk(/, mk(*, b.b, da), b.a)))
    end
    error("no derivative rule for \$(b.f)")
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
Base.show(io::IO, h::Hole) = print(io, "#", h.i)
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
_p!(s,p::Parameter)=(push!(s,p.name);s); _p!(s,::Constant)=s; _p!(s,::Hole)=s
_p!(s,::Provider)=s
_p!(s,u::Unary)=_p!(s,u.a); _p!(s,b::Binary)=(_p!(s,b.a);_p!(s,b.b);s)

shapeof(e)=(cs=ComplexF64[]; (_sh(e,cs),cs))
_sh(p::Parameter,cs)=p
_sh(c::Constant,cs)=(push!(cs,c.val); Hole(length(cs)))
_sh(u::Unary,cs)=Unary(u.f,_sh(u.a,cs))
_sh(b::Binary,cs)=Binary(b.f,_sh(b.a,cs),_sh(b.b,cs))

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
substituteparams(h::Hole, d) = h
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
evalproviders(h::Hole, w) = h
evalproviders(q::Parameter, w) = q
evalproviders(p::Provider, w) = Constant(ComplexF64(p.f(w)))
evalproviders(u::Unary, w) = mk(u.f, evalproviders(u.a, w))
evalproviders(b::Binary, w) = mk(b.f, evalproviders(b.a, w), evalproviders(b.b, w))

evaluate(u::Unary,d)=u.f(evaluate(u.a,d)); evaluate(b::Binary,d)=b.f(evaluate(b.a,d),evaluate(b.b,d))

# ---------- Expr codegen ----------
# A literal is emitted as its real part when it is exactly real, so the
# generated code stays in the parameter's own number type and survives
# ForwardDiff.Dual without being forced to Complex.
_lit(v::ComplexF64) = iszero(imag(v)) ? real(v) : v
_emit(p::Parameter, idx) = :(p[$(idx[p.name])])
_emit(c::Constant, idx)  = _lit(c.val)
_emit(h::Hole, idx)      = :(c[o+$(h.i)])
_emit(u::Unary, idx)     = Expr(:call, u.f, _emit(u.a, idx))
_emit(b::Binary, idx)    = Expr(:call, b.f, _emit(b.a, idx), _emit(b.b, idx))

function buildshapefn(shape, idx, nh)
    body = _emit(shape, idx)
    ex = :((out, p, c, tgt) -> begin
        @inbounds for j in eachindex(tgt)
            o = (j-1)*$nh
            out[tgt[j]] = $body
        end
        nothing
    end)
    return @RuntimeGeneratedFunction(ex)
end

# Real valued literals are stored in a real array. A complex array would
# force every expression that touches a literal into Complex, which costs
# nothing numerically but blocks ForwardDiff: a Complex{Dual} cannot be
# written into the Vector{Dual} that ForwardDiff.jacobian allocates. Since
# the point of this map is to be differentiable with respect to the design
# parameters, the narrowing is load bearing.
_narrow(v::Vector{ComplexF64}) =
    all(iszero, imag.(v)) ? Vector{Float64}(real.(v)) : v

#     ParameterMap(values, names)
#
# A compiled map from a vector of design parameter values to the numeric
# component values of a circuit.
#
# Component values which do not depend on a design parameter are resolved
# once, at construction. The rest are grouped by the *shape* of their
# expression -- the expression with every numeric literal replaced by a hole
# -- so a Floquet weighted line, whose thousands of component values differ
# only in a per cell weight, compiles to a handful of shapes rather than
# thousands of separate expressions. One function is generated per shape and
# evaluated in a tight loop over its literals.
#
# Apply with [`bind!`](@ref).
struct ParameterMap
    constant::AbstractVector
    fns::Vector{Any}
    consts::Vector{Any}
    targets::Vector{Vector{Int}}
    names::Vector{Symbol}
    nshapes::Int
end

function ParameterMap(values::AbstractVector, names::Vector{Symbol})
    idx = Dict(n=>i for (i,n) in enumerate(names)); nameset=Set(names)
    n=length(values); constant=zeros(ComplexF64,n)
    groups=Dict{CircuitValue,Vector{Int}}(); lits=Dict{CircuitValue,Vector{ComplexF64}}()
    nh=Dict{CircuitValue,Int}()
    for i in 1:n
        v=values[i]
        if v isa CircuitValue
            if isempty(intersect(parameters(v),nameset))
                constant[i]=evaluate(v,Dict{Symbol,ComplexF64}())
            else
                sh,cs=shapeof(v)
                push!(get!(Vector{Int},groups,sh),i)
                append!(get!(Vector{ComplexF64},lits,sh),cs); nh[sh]=length(cs)
            end
        else
            constant[i]=ComplexF64(v)
        end
    end
    fns=Any[]; tg=Vector{Int}[]; cs=Any[]
    for (sh,t) in groups
        push!(fns,buildshapefn(sh,idx,nh[sh])); push!(tg,t)
        push!(cs, iszero(nh[sh]) ? Float64[] : _narrow(lits[sh]))
    end
    ParameterMap(_narrow(constant),fns,cs,tg,names,length(fns))
end

#     JacobianMap(values, names)
#
# The exact Jacobian of the component values with respect to the design
# parameters, as one compiled [`ParameterMap`](@ref) per parameter. Costs
# one bind per parameter to evaluate, is allocation free, needs no AD
# package, and is exact for complex component values.
struct JacobianMap
    maps::Vector{ParameterMap}
    names::Vector{Symbol}
    ncomponents::Int
    col::Vector{ComplexF64}   # scratch, so no SubArray is allocated
end

function JacobianMap(values::AbstractVector, names::Vector{Symbol})
    maps = [ParameterMap([derivative(v, nm) for v in values], names)
            for nm in names]
    return JacobianMap(maps, names, length(values),
                       Vector{ComplexF64}(undef, length(values)))
end

function jacobian!(J::AbstractMatrix, jm::JacobianMap, p::AbstractVector)
    col = jm.col
    @inbounds for j in eachindex(jm.names)
        bind!(col, jm.maps[j], p)
        for i in eachindex(col)
            J[i, j] = col[i]
        end
    end
    return J
end

#     bind!(out, m::ParameterMap, p)
#
# Write the component values for the design parameter vector `p` into `out`.
function bind!(out::AbstractVector, m::ParameterMap, p::AbstractVector)
    copyto!(out, m.constant)
    @inbounds for g in eachindex(m.fns)
        m.fns[g](out, p, m.consts[g], m.targets[g])
    end
    return out
end

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

julia> @variables Lj1 Lj2;JosephsonCircuits.componentvaluestonumber([Lj1,Lj1+Lj2],Dict(Lj1=>1e-12,Lj2=>2e-12))
2-element Vector{Float64}:
 1.0e-12
 3.0e-12
```
```jldoctest
# define a frequency dependent impedance function
Zfun(w,R) = ifelse(w>10,R,100*R);
# create symbolic variables including a two argument function
@variables w R
@register_symbolic Zfun(w,R)
# substitute in numerical values and functions for everything but w
out=JosephsonCircuits.componentvaluestonumber([R,Zfun(w,R)],Dict(R=>50));
println(out)
# evaluate with w = 2
println(Symbolics.value.(Symbolics.substitute.(out,(Dict(w=>2),);fold=Val(true))))
# evaluate with w = 11
println(Symbolics.value.(Symbolics.substitute.(out,(Dict(w=>11),);fold=Val(true))))

# output
Any[50, Zfun(w, 50)]
[50, 5000]
[50, 50]
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
    valuetonumber(value::Symbolics.Num,circuitdefs)

If the component value is Symbolics.Num, then try substituting in the definition
from `circuitdefs`.

# Examples
```jldoctest
julia> @variables Lj1;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
3.0e-12

julia> @variables Lj1 Lj2;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
4.0e-12
```
"""
function valuetonumber(value::FrequencyDependent, circuitdefs)
    # a frequency dependent value lowers to a provider leaf and survives to
    # be resolved per mode by freqsubst
    return CircuitValues.Provider(value.f)
end

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

# the `Num` methods live in the Symbolics extension

# """
#     valuetonumber(value::Complex{Symbolics.Num},circuitdefs)

# If the component value `value` is Complex{Symbolics.Num}, then try substituting in the
# definition from `circuitdefs`. This function is currently broken as of 
# Symbolics v7.1.1

# # Examples
# ```jldoctest
# julia> @variables Lj1::Complex;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
# 3.0e-12

# julia> @variables Lj1::Complex Lj2::Complex;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
# ComplexTerm(real(Lj2) + real(Lj1) + im*(imag(Lj2) + imag(Lj1)))
# ```
# """
# function valuetonumber(value::Complex{Symbolics.Num},circuitdefs)
#     return Symbolics.unwrap(Symbolics.substitute(value,circuitdefs))
# end

"""
    valuetonumber(value::Symbolics.SymbolicT, circuitdefs)

If the component value `value` has a type Symbolics.SymbolicT, then try
substituting in the definition from `circuitdefs`.

# Examples
```jldoctest
julia> @syms Lj1;JosephsonCircuits.valuetonumber(Lj1,Dict(Lj1=>3.0e-12))
3.0e-12

julia> @syms Lj1 Lj2;JosephsonCircuits.valuetonumber(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12))
4.0e-12
```
"""
# the `SymbolicT` method lives in the Symbolics extension

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
