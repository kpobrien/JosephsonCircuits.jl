"""
    JosephsonCircuitsSymbolicsExt

Symbolics.jl support for JosephsonCircuits.jl, loaded when Symbolics is
available.

The numeric path of the package uses `CircuitValue`, a dependency-free
expression type, so Symbolics is needed only for:

  - accepting netlists whose component values are `Num`, which is what
    every script written against earlier versions does,
  - the symbolic capacitance and inverse inductance matrices, the only
    place a computer algebra system is genuinely required (a symbolic
    linear solve when mutual inductors are present).
"""
module JosephsonCircuitsSymbolicsExt

using Symbolics, SymbolicUtils
import JosephsonCircuits
const JC = JosephsonCircuits
const CV = JosephsonCircuits.CircuitValues
import LinearAlgebra, SparseArrays
using SparseArrays: sparse, SparseMatrixCSC

const SymAny = Union{Num,SymbolicUtils.BasicSymbolic}


# ---------------------------------------------------------------------
# Num -> CircuitValue
# ---------------------------------------------------------------------
"""
    tocircuitvalue(x)

Convert a Symbolics expression to a `CircuitValue`. Only the closed
operator set of the circuit value type is accepted; anything else is a
hard error rather than a silent approximation.
"""
tocircuitvalue(x::Number) = CV.Constant(x)
tocircuitvalue(x::CV.CircuitValue) = x
tocircuitvalue(n::Num) = tocircuitvalue(Symbolics.value(n))

# Symbolics represents a complex expression as a Complex whose parts are
# themselves symbolic, which is a Number and would otherwise be caught by
# the method above and fail to convert.
tocircuitvalue(z::Complex{<:Union{Num,SymbolicUtils.BasicSymbolic}}) =
    tocircuitvalue(real(z)) + im*tocircuitvalue(imag(z))

function tocircuitvalue(x)
    x isa Number && return CV.Constant(x)
    # "not a call" is a stable definition of a leaf; `issym` does not
    # recognize the current SymbolicUtils leaf types.
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

# ---------------------------------------------------------------------
# the hooks the core leaves open
# ---------------------------------------------------------------------
# Keep Symbolics' own partial-substitution semantics rather than routing
# through CircuitValue evaluation. A value may legitimately still contain a
# free variable after substitution -- the symbolic frequency variable is
# resolved later by `freqsubst`, not here -- and forcing full evaluation
# turns that into a KeyError. `CircuitValue` has no symbolic frequency
# support by design, so this is the one place the two paths genuinely
# differ.
JC.valuetonumber(v::Num, circuitdefs) =
    Symbolics.value(Symbolics.substitute(v, circuitdefs; fold=Val(true)))
JC.valuetonumber(v::SymbolicUtils.BasicSymbolic, circuitdefs) =
    Symbolics.value(Symbolics.substitute(v, circuitdefs; fold=Val(true)))

JC.unwrapvalue(v::Num) = Symbolics.value(v)
# Fold and unwrap like `valuetonumber`: a fully resolved entry must come
# back as a plain Julia number (the new SymbolicUtils keeps folded
# constants wrapped, and an unconditional `checkissymbolic` on the wrapper
# would reject them); a partially resolved entry stays symbolic for the
# caller to diagnose.
JC.substitutefreq(v::SymAny, symfreqvar, w) =
    Symbolics.value(Symbolics.substitute(v, symfreqvar => w; fold=Val(true)))
JC.substitutedefs(v::SymAny, circuitdefs) =
    Symbolics.substitute(v, circuitdefs)

# Both wrappers must be covered. The core previously tested against
# `Symbolics.SymbolicT`, which catches the unwrapped `BasicSymbolic` as
# well as `Num`; an extension that only handles `Num` silently lets
# unwrapped symbolic values be treated as numeric, and they fail much
# later trying to convert to ComplexF64.
JC.checkissymbolic(a::Num) = !(Symbolics.value(a) isa Number)
JC.checkissymbolic(a::SymbolicUtils.BasicSymbolic) = true
JC.circuitvariables(a::SymAny) = Symbolics.get_variables(a)
JC.issymbolicvaluetype(::Type{T}) where {T<:SymAny} = true

# the symbolic inverse inductance matrix, the one genuine CAS use
function JC.calcsymbolicinvLn(L::AbstractMatrix, Lb, Rbn)
    # `sym_lu` is defined for `Num`, but the caller's matrix takes its
    # element type from the component values, which may be unwrapped
    # `BasicSymbolic` (from `@syms`) rather than `Num` (from `@variables`).
    # The normalization belongs here, in the extension that knows about
    # `Num`, rather than in the core, which must stay free of Symbolics
    # types.
    s = sparse(transpose(Rbn[Lb.nzind,:]) *
        (Symbolics.sym_lu(Num.(Matrix(L))) \ Matrix(Rbn[Lb.nzind,:])))
    return SparseMatrixCSC(s.m, s.n, s.colptr, s.rowval, s.nzval)
end

end # module
