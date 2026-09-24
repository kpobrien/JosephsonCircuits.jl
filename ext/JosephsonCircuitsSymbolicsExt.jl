"""
    JosephsonCircuitsSymbolicsExt

Symbolics.jl support, loaded with `using Symbolics`.

The numeric path of the package works on `CircuitValue`, its own
dependency free expression type, so Symbolics is needed only to accept
component values written as `Num` or `BasicSymbolic` expressions, by
substituting the definitions into them directly.
"""
module JosephsonCircuitsSymbolicsExt

using Symbolics, SymbolicUtils
import JosephsonCircuits
const JC = JosephsonCircuits
const CV = JosephsonCircuits.CircuitValues

const SymAny = Union{Num,SymbolicUtils.BasicSymbolic}


# === methods of the core's value handling functions for symbolic values ===

# Substitution keeps Symbolics' partial substitution semantics rather than
# lowering to a `CircuitValue` first: a value may still contain a free
# variable afterwards (the symbolic frequency variable, resolved per mode
# by `freqsubst`), and a full evaluation would turn that into a KeyError.
JC.valuetonumber(v::Num, circuitdefs) =
    Symbolics.value(Symbolics.substitute(v, circuitdefs; fold=Val(true)))
JC.valuetonumber(v::SymbolicUtils.BasicSymbolic, circuitdefs) =
    Symbolics.value(Symbolics.substitute(v, circuitdefs; fold=Val(true)))

# the port number of a deprecated tuple netlist entry (circuit/legacy.jl)
JC.unwrapvalue(v::Num) = Symbolics.value(v)
# A symbolic value carries no frequency of its own: it is either already a
# number, in which case unwrapping it makes that visible (SymbolicUtils
# keeps folded constants wrapped, and `checkissymbolic` on the wrapper
# would reject them), or it still depends on a parameter the definitions
# did not give and stays symbolic for the caller to diagnose.
JC.substitutefreq(v::SymAny, w) = Symbolics.value(v)
JC.substitutedefs(v::SymAny, circuitdefs) =
    Symbolics.substitute(v, circuitdefs)

# Both `Num` and the unwrapped `BasicSymbolic` must be covered: an
# unwrapped symbolic value which is missed here is treated as numeric and
# fails much later, converting to ComplexF64.
JC.checkissymbolic(a::Num) = !(Symbolics.value(a) isa Number)
JC.checkissymbolic(a::SymbolicUtils.BasicSymbolic) = true
JC.circuitvariables(a::SymAny) = Symbolics.get_variables(a)

# the design sensitivities: a `Num` key of the definitions names the
# parameter, and the derivative of a `Num` value with respect to the
# parameter of that name is Symbolics' own, evaluated at the definitions
JC.definitionname(k::SymAny) = Symbolics.tosymbol(k; escape = false)
function JC.designderivative(v::SymAny, name::Symbol, definitions)
    # the variables are iterated rather than indexed: `get_variables` gives
    # an ordered set in recent versions of Symbolics and a vector in older
    # ones, and only iteration is common to the two
    for variable in Symbolics.get_variables(v)
        Symbolics.tosymbol(variable; escape = false) === name || continue
        d = JC.valuetonumber(Symbolics.derivative(v, variable), definitions)
        d isa Number || throw(ArgumentError(lazy"the derivative of the value $(v) with respect to $(name) is $(d) at these definitions, which is not a number; design sensitivities need every parameter defined."))
        return ComplexF64(d)
    end
    return zero(ComplexF64)
end

end # module
