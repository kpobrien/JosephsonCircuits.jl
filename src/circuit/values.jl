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

A sweep evaluates each distinct function, as `===` tells them apart, once
at each mode frequency however many components share it, so `f` depends
on the frequency alone.

The value describes the element at positive frequencies. A mode of
negative frequency, such as an idler, takes the complex conjugate of the
value at the magnitude of its frequency, `conj(f(abs(w)))`, which is the
value of a real element there, as the data of a scattering block is
extended by [`ConjugateSymmetry`](@ref).

The value is the whole function: it does not combine with numbers or
other values, so a value derived from others is computed inside the
closure.
"""
struct FrequencyDependent{F}
    f::F
end

# === resolving a written value to a number ===
#
# A component value is written as a number, as a frequency dependent
# value, as a name (a symbol or a string) the definitions give, or as a
# symbolic value of the Symbolics extension. A definition is a number or a
# frequency dependent value; a value derived from others is computed in
# Julia and defined as the number it comes to.

"""
    componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)

Resolve each component value in `componentvalues` at the definitions
`circuitdefs` with [`valuetonumber`](@ref), the definitions gathered by
name once for the whole table ([`definitionsbyname`](@ref)), and return
the results as a `Vector{Any}`: the table mixes real and complex values (a
port's entry is its reference impedance), symbolic values and frequency
dependent values, and the groups the assembly reads are typed when they
are gathered from it (see `grouptype`), so the table itself has one type
for every circuit.

# Examples
```jldoctest
julia> JosephsonCircuits.componentvaluestonumber([:Lj1,:Lj2],Dict(:Lj1=>1e-12,:Lj2=>2e-12))
2-element Vector{Any}:
 1.0e-12
 2.0e-12
```
"""
function componentvaluestonumber(componentvalues::Vector,circuitdefs::Dict)
    byname = definitionsbyname(circuitdefs)
    return Any[valuetonumber(value, byname) for value in componentvalues]
end

"""
    valuetonumber(value, circuitdefs)

The component value `value` at the definitions `circuitdefs`. A name,
written as a symbol or a string, is the definition of that name under
whichever key gives it (see [`definitionsbyname`](@ref)), and a name the
definitions do not give is returned as it is, for the check of the values
to report naming its component ([`checkcomponentvaluesdefined`](@ref)).
Any other value, a number or a frequency dependent value, is returned
unchanged. The Symbolics extension adds the method for its symbolic
values, which substitutes the definitions of their variables.

# Examples
```jldoctest
julia> JosephsonCircuits.valuetonumber(:Lj1,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0e-12

julia> JosephsonCircuits.valuetonumber("Lj1",Dict(:Lj1=>1e-12))
1.0e-12

julia> JosephsonCircuits.valuetonumber(1.0,Dict(:Lj1=>1e-12,:Lj2=>2e-12))
1.0
```
"""
valuetonumber(value, circuitdefs) = value
valuetonumber(value::Union{Symbol,AbstractString}, circuitdefs) =
    get(definitionsbyname(circuitdefs), definitionname(value), value)

"""
    definitionsbyname(circuitdefs)

The definitions keyed by the name each key gives (see `definitionname`),
with their values as given; a key which gives no name is left out. A name
may be defined under its symbol, its string or, with the Symbolics
extension, its variable, and one defined under two of them with different
values is refused rather than resolved by the order of the dictionary. A
table already keyed by name is returned as it is.
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

"""
    definitiontable(circuitdefs)

The component definitions as a `Dict{Any,Any}`, whatever the key and value
types of the dictionary given.
"""
definitiontable(circuitdefs::Dict{Any,Any}) = circuitdefs
definitiontable(circuitdefs::AbstractDict) = Dict{Any,Any}(circuitdefs)

# The name a definition key or a value written as a name gives, and
# nothing for a key which gives none: a symbol is its own name and a
# string its symbol. The solver cache and the design sensitivities ask
# which name a key or a value is for. The Symbolics extension adds the
# method for its variables.
definitionname(k::Symbol) = k
definitionname(k::AbstractString) = Symbol(k)
definitionname(k) = nothing

# === checking the values of a circuit ===

"""
    checkissymbolic(a)

Whether `a` is a symbolic value of the Symbolics extension, which adds the
methods, rather than a number: without the extension no value is.

# Examples
```jldoctest
julia> JosephsonCircuits.checkissymbolic(1.0)
false
```
"""
checkissymbolic(a) = false

"""
    circuitvariables(a)

The variables of a symbolic value of the Symbolics extension, which
defines the methods, for the check of the values to name.
"""
function circuitvariables end

"""
    substitutefreq(value, w)

Resolve a component value at the mode frequency `w`: the identity for a
plain number, and for a [`FrequencyDependent`](@ref) value its function at
the magnitude of `w`. A value states the element at positive frequencies,
and a mode of negative frequency takes its conjugate where the value is
placed (see [`modevalue`](@ref)), which is the value a real element has
there. A value which does not resolve to a number is returned as it is,
for the caller to diagnose. The Symbolics extension adds the `Num` method.
"""
substitutefreq(value, w) = value
substitutefreq(value::FrequencyDependent, w) = value.f(abs(w))

"""
    checkcomponentvaluesdefined(componentnames::Vector, vvn::Vector)

Check that every circuit component value resolved: a value which is still
a name, one the circuit definitions `circuitdefs` do not give, or a
symbolic value with a variable they do not define, is reported in an
`ArgumentError` naming its component and the names it depends on. A
frequency dependent value is a closure of the frequency, resolved per
frequency by [`freqsubst`](@ref), and passes. Called by [`hbnlsolve`](@ref)
and [`hblinsolve`](@ref) before any computation, and by the binding of the
values to the circuit when a name is left, so a forgotten definition fails
immediately with the actual cause.

# Examples
```jldoctest
julia> JosephsonCircuits.checkcomponentvaluesdefined(["P1","R1"], Any[1, FrequencyDependent(w -> 1/(w*50.0))])

julia> try JosephsonCircuits.checkcomponentvaluesdefined(["P1","R1"], Any[1, :R2]) catch e; occursin("R1 has the value", sprint(showerror, e)) end
true
```
"""
function checkcomponentvaluesdefined(componentnames::Vector, vvn::Vector)
    messages = String[]
    for i in eachindex(vvn)
        v = vvn[i]
        # a value of the types a resolved table holds, told without a
        # dynamic call
        v isa Union{Float64,ComplexF64,FrequencyDependent} && continue
        undefined = undefinednames(v)
        isempty(undefined) && continue
        push!(messages, string("The component ", componentnames[i],
            " has the value ", v, ", which depends on [",
            join(string.(undefined), ", "), "], not defined as numbers."))
    end
    if !isempty(messages)
        throw(ArgumentError(join(messages, " ")*" Add the missing "*
            "definitions to the circuit definitions dictionary circuitdefs. "*
            "A value derived from others is computed in Julia and defined "*
            "as the number it comes to, and a value which depends on the "*
            "frequency is written as a FrequencyDependent closure of the "*
            "frequency."))
    end
    return nothing
end

# the names a value still depends on: itself when it is a name the
# definitions did not give, and the variables of a symbolic value
undefinednames(v::Union{Symbol,AbstractString}) = (v,)
undefinednames(v) = checkissymbolic(v) ? circuitvariables(v) : ()
