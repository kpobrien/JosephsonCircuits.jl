"""
    JosephsonCircuitsSymbolicsExt

Symbolics.jl support, loaded with `using Symbolics`.

The numeric path of the package works on `CircuitValue`, its own
dependency free expression type, so Symbolics is needed only to:

  - accept component values written as `Num` or `BasicSymbolic`
    expressions, by lowering them to `CircuitValue` or substituting into
    them directly;
  - compute the symbolic inverse inductance matrix of a circuit with
    mutual inductors, which needs a symbolic linear solve
    (`Symbolics.sym_lu`).
"""
module JosephsonCircuitsSymbolicsExt

using Symbolics, SymbolicUtils
import JosephsonCircuits
const JC = JosephsonCircuits
const CV = JosephsonCircuits.CircuitValues
import LinearAlgebra, SparseArrays
using SparseArrays: sparse, SparseMatrixCSC

const SymAny = Union{Num,SymbolicUtils.BasicSymbolic}


# === Num -> CircuitValue ===

"""
    tocircuitvalue(x)

Convert a Symbolics expression to a `CircuitValue`. Only the operators
`+ - * / ^ inv sqrt exp log conj` are accepted; any other operator is an
error rather than a silent approximation.
"""
tocircuitvalue(x::Number) = CV.Constant(x)
tocircuitvalue(x::CV.CircuitValue) = x
tocircuitvalue(n::Num) = tocircuitvalue(Symbolics.value(n))

# Symbolics represents a complex expression as a `Complex` whose parts are
# symbolic. That is a `Number`, so without this method it would be caught
# by the `Number` method above and fail to convert.
tocircuitvalue(z::Complex{<:Union{Num,SymbolicUtils.BasicSymbolic}}) =
    tocircuitvalue(real(z)) + im*tocircuitvalue(imag(z))

function tocircuitvalue(x)
    x isa Number && return CV.Constant(x)
    # a leaf is anything which is not a call; `issym` does not recognize
    # every SymbolicUtils leaf type
    if SymbolicUtils.iscall(x)
        op   = SymbolicUtils.operation(x)
        args = map(tocircuitvalue, SymbolicUtils.arguments(x))
        op === (+) && return reduce(+, args)
        op === (*) && return reduce(*, args)
        op === (/) && return args[1]/args[2]
        op === (^) && return args[1]^args[2]
        op === (-) && return length(args) == 1 ? -args[1] : args[1]-args[2]
        for f in (inv, sqrt, exp, log, conj)
            op === f && return f(args[1])
        end
        error("no CircuitValue equivalent for the operator $(op); the "*
              "closed operator set is + - * / ^ inv sqrt exp log conj")
    end
    # a leaf: a variable, or a numeric literal which survived wrapping
    xv = try Symbolics.value(x) catch; x end
    xv isa Number && return CV.Constant(xv)
    str = string(xv)
    n = tryparse(Float64, str)
    isnothing(n) || return CV.Constant(n)
    return CV.Parameter(Symbol(str))
end

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
JC.issymbolicvaluetype(::Type{T}) where {T<:SymAny} = true

# The symbolic inverse inductance matrix `transpose(Rbn) * (L \ Rbn)`,
# restricted to the inductive branches, which needs a symbolic linear solve.
function JC.calcsymbolicinvLn(L::AbstractMatrix, Lb, Rbn)
    # `sym_lu` is defined for `Num`, but the matrix takes its element type
    # from the component values, which may be unwrapped `BasicSymbolic`
    # (from `@syms`) rather than `Num` (from `@variables`), so it is
    # converted here.
    s = sparse(transpose(Rbn[Lb.nzind,:]) *
        (Symbolics.sym_lu(Num.(Matrix(L))) \ Matrix(Rbn[Lb.nzind,:])))
    return SparseMatrixCSC(s.m, s.n, s.colptr, s.rowval, s.nzval)
end

end # module
