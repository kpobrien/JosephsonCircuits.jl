"""
    JosephsonCircuitsKrylovExt

Krylov.jl as the linear solver for the Newton step, selected with
`linearsolver = KrylovJL(:gmres)`.

Only the linear solve changes. The forcing term, the line search, the
preconditioner escalation and the stagnation handling of `nlsolvekrylov!`
are untouched, and the mode coupling preconditioner is passed through
unchanged because it is applied by `ldiv!`, which is what Krylov.jl's `N`
argument consumes.
"""
module JosephsonCircuitsKrylovExt

using Krylov, LinearAlgebra
import JosephsonCircuits
const JC = JosephsonCircuits

# The seam now receives a `mul!`-able operator, so there is nothing to
# wrap: `nlsolvekrylov!` normalizes the Jacobian with `asoperator` and the
# preconditioner is applied through the package's `AbstractPreconditioner`
# interface. Only the preconditioner still arrives as a closure.
struct MopWrap{F} <: JC.AbstractPreconditioner; Mop!::F; end
JC.applypreconditioner!(z, m::MopWrap, r) = (m.Mop!(z, r); z)

function JC.hblinearsolve!(ls::JC.KrylovJL, deltax, jvp, F, ws, Mop!;
        rtol, atol, maxrestarts)
    n = length(F)
    A = jvp
    solver = getfield(Krylov, ls.method)
    # the internal solver restarts up to maxrestarts times over a Krylov
    # space of size(ws.H,2), so the comparable iteration budget is their
    # product
    itmax = max(size(ws.H, 2)*(maxrestarts+1), 10)
    # Krylov.jl defaults `atol` to sqrt(eps()), which is right for a
    # standalone linear solve and wrong inside a Newton loop: the right hand
    # side is the residual being driven to zero, so an absolute floor
    # eventually swallows it whole and every solve returns success having
    # done nothing. The seam is handed `atol = ftol/10` by
    # `nlsolvekrylov!`; an explicit `atol` in the solver's own keywords
    # still wins.
    kw = haskey(ls.kwargs, :atol) ? ls.kwargs : merge((; atol = atol), ls.kwargs)
    x, st = if isnothing(Mop!)
        solver(A, F; rtol = rtol, itmax = itmax, kw...)
    else
        solver(A, F; N = JC.SizedPreconditioner(MopWrap(Mop!), n),
            rtol = rtol, itmax = itmax, kw...)
    end
    copyto!(deltax, x)
    # the output shape `nlsolvekrylov!` records: `cycles` has no analogue,
    # so it reports one
    return (converged = st.solved, residual = isempty(st.residuals) ?
                norm(F) : last(st.residuals),
            iterations = st.niter, cycles = 1,
            reason = st.solved ? :converged : :notconverged,
            precondtime = NaN)
end

end # module
