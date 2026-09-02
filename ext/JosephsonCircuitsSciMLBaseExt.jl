"""
    JosephsonCircuitsSciMLBaseExt

The harmonic balance system as a `SciMLBase.NonlinearProblem`, so that
NonlinearSolve.jl can solve it. This gives access to globalization
strategies the package does not implement (trust region, pseudo-transient
continuation, Levenberg-Marquardt, polyalgorithms).

The problem is posed in the equivalent real representation, because the
harmonic balance residual is not complex differentiable.
"""
module JosephsonCircuitsSciMLBaseExt

using SciMLBase, SparseArrays
import JosephsonCircuits
const JC = JosephsonCircuits

"""
    SciMLBase.NonlinearProblem(prob::HBNonlinearProblem; u0 = prob.u0,
        sparsejac = true)

A `NonlinearProblem` over the real representation of `prob`, starting
from `u0`. The exact matrix-free Jacobian-vector product is supplied as
`jvp`, so a Krylov linear solver never assembles or differentiates
anything. When `sparsejac = true` and `prob` carries an assembled
Jacobian, `jac` and `jac_prototype` are supplied as well for the direct
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
