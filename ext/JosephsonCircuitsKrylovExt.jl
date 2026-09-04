"""
    JosephsonCircuitsKrylovExt

Krylov.jl as the linear solver of the Newton step of `nlsolvekrylov!`,
selected with `linearsolver = KrylovJL(:gmres)` (or any other Krylov.jl
solver name). Only the linear solve changes: the forcing term, the line
search, the preconditioner escalation and the stagnation handling are those
of `nlsolvekrylov!`, and the package's preconditioner is passed to Krylov.jl
as its right preconditioner `N`.
"""
module JosephsonCircuitsKrylovExt

using Krylov, LinearAlgebra
import JosephsonCircuits
const JC = JosephsonCircuits

# `nlsolvekrylov!` hands over the Jacobian as an operator supporting
# `mul!`, which Krylov.jl accepts directly, and the preconditioner as an
# `AbstractPreconditioner`, which Krylov.jl applies through `mul!`. A bare
# closure `Mop!(z, r)`, the older spelling, is wrapped into that interface.
# Either way the preconditioner is applied standalone here: a form which
# fuses its application with the Jacobian product inside the package's own
# GMRES (`preconditionedproduct!`) falls back to its unfused, exact form.
struct MopWrap{F} <: JC.AbstractPreconditioner; Mop!::F; end
JC.applypreconditioner!(z, m::MopWrap, r) = (m.Mop!(z, r); z)
aspreconditioner(M::JC.AbstractPreconditioner) = M
aspreconditioner(M) = MopWrap(M)

function JC.hblinearsolve!(ls::JC.KrylovJL, deltax, jvp, F, ws, Mop!;
        rtol, atol, maxrestarts)
    n = length(F)
    A = jvp
    solver = getfield(Krylov, ls.method)
    # the package's own GMRES restarts up to `maxrestarts` times over a
    # Krylov space of `size(ws.H, 2)`, so the comparable iteration budget
    # is their product
    itmax = max(size(ws.H, 2)*(maxrestarts+1), 10)
    # Krylov.jl defaults `atol` to `sqrt(eps())`, which is wrong inside a
    # Newton loop: the right hand side is the residual being driven to
    # zero, and an absolute floor would eventually accept every solve
    # without doing anything. `nlsolvekrylov!` passes `atol = ftol/10`; an
    # explicit `atol` in the solver's own keywords wins.
    kw = haskey(ls.kwargs, :atol) ? ls.kwargs : merge((; atol = atol), ls.kwargs)
    x, st = if isnothing(Mop!)
        solver(A, F; rtol = rtol, itmax = itmax, kw...)
    else
        solver(A, F; N = JC.SizedPreconditioner(aspreconditioner(Mop!), n),
            rtol = rtol, itmax = itmax, kw...)
    end
    copyto!(deltax, x)
    # the record `nlsolvekrylov!` expects; Krylov.jl has no notion of
    # restart cycles, so one is reported
    return (converged = st.solved, residual = isempty(st.residuals) ?
                norm(F) : last(st.residuals),
            iterations = st.niter, cycles = 1,
            reason = st.solved ? :converged : :notconverged,
            precondtime = NaN)
end

end # module
