"""
    JosephsonCircuitsSymbolicsExt

Symbolics.jl support, loaded with `using Symbolics`: component values
written as `Num` or `BasicSymbolic` expressions, with the definitions of
their variables substituted; the matrices of a circuit as expressions
(`symbolicmatrices`); and the expressions of a netlist file
(`import_netlist`).
"""
module JosephsonCircuitsSymbolicsExt

using Symbolics, SymbolicUtils
import JosephsonCircuits
const JC = JosephsonCircuits

const SymAny = Union{Num,SymbolicUtils.BasicSymbolic}


# === methods of the core's value handling functions for symbolic values ===

# A symbolic value at the definitions. A bare variable is a name, and is
# its definition whatever that is, as a symbol is; an expression has the
# definitions of its variables which are numbers substituted. A variable
# the definitions do not give stays in the value, for the check of the
# values to name. Each variable is looked up under whichever key names it,
# its own, its symbol or its string.
function JC.valuetonumber(v::SymAny, circuitdefs)
    x = Symbolics.value(v)
    x isa SymbolicUtils.BasicSymbolic || return x
    byname = JC.definitionsbyname(circuitdefs)
    if SymbolicUtils.issym(x)
        def = get(byname, JC.definitionname(x), x)
        return def isa SymAny ? Symbolics.value(def) : def
    end
    d = Dict{Any,Any}()
    for variable in Symbolics.get_variables(x)
        def = get(byname, JC.definitionname(variable), nothing)
        def isa Number && (d[variable] = def)
    end
    return Symbolics.value(Symbolics.substitute(x, d; fold = Val(true)))
end

# the port number of a deprecated tuple netlist entry (circuit/legacy.jl)
JC.unwrapvalue(v::Num) = Symbolics.value(v)
# A symbolic value carries no frequency of its own: it is either already a
# number, in which case unwrapping it makes that visible (SymbolicUtils
# keeps folded constants wrapped, and `checkissymbolic` on the wrapper
# would reject them), or it still depends on a parameter the definitions
# did not give and stays symbolic for the caller to diagnose.
JC.substitutefreq(v::SymAny, w) = Symbolics.value(v)

# Both `Num` and the unwrapped `BasicSymbolic` must be covered: an
# unwrapped symbolic value which is missed here is treated as numeric and
# fails much later, converting to ComplexF64.
JC.checkissymbolic(a::Num) = !(Symbolics.value(a) isa Number)
JC.checkissymbolic(a::SymbolicUtils.BasicSymbolic) = true
JC.circuitvariables(a::SymAny) = Symbolics.get_variables(a)

# a `Num` variable names the parameter of its symbol, as a key of the
# definitions and as a component value
JC.definitionname(k::SymAny) = Symbolics.tosymbol(k; escape = false)

# the variable of a name, for `symbolicmatrices`
variable(name::Symbol) = Symbolics.variable(name)

# A value of a netlist file written as an expression, as Symbolics reads
# it (`import_netlist`): the operators and the functions and constants
# named here are evaluated, and any other name is a variable, so a netlist
# runs no code of its own.
const netlistnames = Dict{Symbol,Any}(:im => im, :pi => pi, :inv => inv,
    :conj => conj, :real => real, :imag => imag)
netlistvalue(e) = Symbolics.parse_expr_to_symbolic(e, netlistnames)

end # module
