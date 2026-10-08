"""
    JosephsonCircuitsKrylovExt

Krylov.jl as the linear solver of the Newton step of `nlsolvekrylov!`, and
of the transient's step solved iteratively, selected with
`linearsolver = JosephsonCircuits.KrylovJL(:gmres)` (or any other Krylov.jl
solver name). Only the linear solve changes: the forcing term, the line
search, the preconditioner escalation and the stagnation handling are
those of `nlsolvekrylov!`, and the package's preconditioner is passed to
Krylov.jl as its right preconditioner `N`. A solve is given the iteration
budget of the package's own GMRES, and a method which can restart
(`:gmres`, `:fgmres`, `:fom`) restarts at its restart length. The Krylov.jl
workspace is made at the first solve of a system and kept in the
package's workspace for the others. Deflation harvesting is
unavailable here: it reads the package's own Arnoldi workspace, which
Krylov.jl does not expose, so the per cycle callback `oncycle` is accepted
and ignored.
"""
module JosephsonCircuitsKrylovExt

using Krylov, LinearAlgebra
import JosephsonCircuits
const JC = JosephsonCircuits

# `nlsolvekrylov!` hands over the Jacobian as an operator supporting
# `mul!`, which Krylov.jl accepts directly, and the preconditioner as an
# `AbstractPreconditioner`, which Krylov.jl applies through `mul!`. The
# transient's step hands over a bare closure `Mop!(z, r)`, which is
# wrapped into that interface.
struct MopWrap{F} <: JC.AbstractPreconditioner; Mop!::F; end
JC.applypreconditioner!(z, m::MopWrap, r) = (m.Mop!(z, r); z)
aspreconditioner(M::JC.AbstractPreconditioner) = M
aspreconditioner(M) = MopWrap(M)

# An operator of dimension `n` applied through `mul!` and reporting the
# element type `T` of the iteration. Krylov.jl compares the element type of
# the operator with that of the vectors and warns when they differ, and
# the package's operators and preconditioners report double precision
# whatever the precision of the iteration they serve.
struct TypedOperator{T,A}
    op::A
    n::Int
end
TypedOperator{T}(op, n::Integer) where {T} = TypedOperator{T,typeof(op)}(op, Int(n))
Base.size(A::TypedOperator) = (A.n, A.n)
Base.size(A::TypedOperator, i::Integer) = A.n
Base.eltype(::TypedOperator{T}) where {T} = T
LinearAlgebra.mul!(y::AbstractVector, A::TypedOperator, x::AbstractVector) =
    mul!(y, A.op, x)

# the methods of Krylov.jl which restart after `memory` iterations when
# asked to, and otherwise keep every basis vector until `itmax`
const RESTARTABLE = (:gmres, :fgmres, :fom)

# the keywords which size a Krylov.jl workspace rather than set a solve
const WORKSPACEKEYWORDS = (:memory, :window)

function JC.hblinearsolve!(ls::JC.KrylovJL, deltax, jvp, F, ws, Mop!;
        rtol, atol, maxrestarts, oncycle = nothing)
    n = length(F)
    T = real(eltype(F))
    A = TypedOperator{T}(jvp, n)
    # the budget of the package's own GMRES: `maxrestarts` cycles of the
    # workspace's restart length, and for a method which restarts, the same
    # cycles, so its basis is that length rather than every step taken
    m = size(ws.H, 2)
    itmax = m*maxrestarts
    restarts = ls.method in RESTARTABLE ? (; restart = true, memory = m) : (;)
    # Krylov.jl defaults `atol` to `sqrt(eps())`, which is wrong inside a
    # Newton loop: the right hand side is the residual being driven to
    # zero, and an absolute floor would eventually accept every solve
    # without doing anything. `nlsolvekrylov!` passes a tenth of its own
    # `atol`.
    # Krylov.jl records the residual history only when asked, and the
    # residual a solve reports is the last of that history, so it is asked
    # for: without it, a solve would report its starting residual. The
    # solver's own keywords win over all of these, and the tolerances are
    # taken in the precision of the iteration, which Krylov.jl requires.
    kw = merge((; atol = atol, history = true), restarts, ls.kwargs)
    kw = merge(kw, (; atol = T(kw.atol)))
    wnames = filter(in(keys(kw)), WORKSPACEKEYWORDS)
    workspace = krylovworkspace!(ws, ls.method, A, F, kw[wnames])
    skw = Base.structdiff(kw, NamedTuple{wnames})
    if isnothing(Mop!)
        krylov_solve!(workspace, A, F; rtol = T(rtol), itmax = itmax, skw...)
    else
        krylov_solve!(workspace, A, F;
            N = TypedOperator{T}(aspreconditioner(Mop!), n), rtol = T(rtol),
            itmax = itmax, skw...)
    end
    copyto!(deltax, Krylov.solution(workspace))
    st = Krylov.statistics(workspace)
    # the record `nlsolvekrylov!` expects; Krylov.jl has no notion of
    # restart cycles, so one is reported, and no product count, so the
    # diagnostics count zero Jacobian products for a Krylov.jl solve
    return (converged = st.solved, residual = isempty(st.residuals) ?
                norm(F) : last(st.residuals),
            iterations = st.niter, cycles = 1,
            reason = st.solved ? :converged : :notconverged,
            precondtime = NaN)
end

# The Krylov.jl workspace of `method` for the system of `F`, sized by the
# keywords `wkw`, kept in the package's workspace `ws` between the solves
# of one system. Krylov.jl allocates a restarting method's basis in full
# when the workspace is made, `memory` vectors, so a workspace made per
# solve would allocate the whole basis at every Newton step, whatever the
# solve took.
function krylovworkspace!(ws, method::Symbol, A, F, wkw::NamedTuple)
    key = (method, wkw, typeof(F))
    held = ws.external
    held isa Tuple && first(held) == key && return last(held)
    workspace = krylov_workspace(Val(method), A, F; wkw...)
    ws.external = (key, workspace)
    return workspace
end

end # module
