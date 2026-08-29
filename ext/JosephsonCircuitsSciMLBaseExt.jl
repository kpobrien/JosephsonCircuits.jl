"""
    JosephsonCircuitsSciMLBaseExt

The harmonic balance system as a `SciMLBase.NonlinearProblem`, so
NonlinearSolve.jl can drive it.

The value here is not speed: `nlsolvekrylov!` already has a good
Jacobian-free Newton-Krylov with a domain specific preconditioner. It is
access to globalization strategies this package does not implement --
trust region, pseudo-transient continuation, Levenberg-Marquardt, and the
polyalgorithm fallback chain -- which matter for getting onto a strongly
pumped branch from a cold start.

The problem is posed in the equivalent real representation, because the
harmonic balance residual is not complex differentiable.
"""
module JosephsonCircuitsSciMLBaseExt

using SciMLBase, SparseArrays
import JosephsonCircuits
const JC = JosephsonCircuits

"""
    SciMLBase.NonlinearProblem(prob::HBNonlinearProblem; u0, sparsejac)

A `NonlinearProblem` over the real representation.

The exact matrix-free Jacobian-vector product is supplied as `jvp`, so a
Krylov `linsolve` never assembles or differentiates anything: one product
is two transforms plus the linear term. `jac` and `jac_prototype` are
supplied when the problem carries an assembled Jacobian, for the direct
factorization path.
"""
function SciMLBase.NonlinearProblem(prob::JC.HBNonlinearProblem;
        u0 = prob.u0, sparsejac::Bool = true)
    f! = (F, u, p) -> (JC.hbresidual!(F, prob, u); nothing)
    jvp! = (Jv, v, u, p) -> (JC.hbjvp!(Jv, prob, u, v); nothing)
    nf = if sparsejac && !isnothing(prob.jacobian)
        jac! = (J, u, p) -> (JC.hbjacobian!(J, prob, u); nothing)
        NonlinearFunction{true}(f!; jvp = jvp!, jac = jac!,
            jac_prototype = copy(prob.jacobian),
            resid_prototype = similar(prob.u0))
    else
        NonlinearFunction{true}(f!; jvp = jvp!,
            resid_prototype = similar(prob.u0))
    end
    return NonlinearProblem(nf, collect(u0), SciMLBase.NullParameters())
end

end # module
