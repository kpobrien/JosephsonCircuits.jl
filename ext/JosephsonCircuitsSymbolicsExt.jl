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

JC.unwrapvalue(v::Num) = Symbolics.value(v)
# Fold and unwrap so that a fully resolved value comes back as a plain
# number (SymbolicUtils keeps folded constants wrapped otherwise, and
# `checkissymbolic` on the wrapper would reject them); a partially resolved
# value stays symbolic for the caller to diagnose.
JC.substitutefreq(v::SymAny, symfreqvar, w) =
    Symbolics.value(Symbolics.substitute(v, symfreqvar => w; fold=Val(true)))
JC.substitutedefs(v::SymAny, circuitdefs) =
    Symbolics.substitute(v, circuitdefs)

# Both `Num` and the unwrapped `BasicSymbolic` must be covered: an
# unwrapped symbolic value which is missed here is treated as numeric and
# fails much later, converting to ComplexF64.
JC.checkissymbolic(a::Num) = !(Symbolics.value(a) isa Number)
JC.checkissymbolic(a::SymbolicUtils.BasicSymbolic) = true
JC.circuitvariables(a::SymAny) = Symbolics.get_variables(a)

end # module
