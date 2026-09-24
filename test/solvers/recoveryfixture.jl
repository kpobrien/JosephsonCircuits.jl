# A solver which really solves, with failures injected by attempt number.
#
# The recovery tests need an attempt which fails at a chosen point and one
# which converges, both on the same circuit, so that what a cache or a
# staged schedule does after a failure can be asserted without relying on
# a circuit that happens to be hard to solve on one platform and easy on
# another. `failat` names the attempts, counted from one, which return a
# recognizable state and report failure; every other attempt runs Newton's
# method on the package's own residual and Jacobian. The returned starts
# and solutions record what each attempt was handed and what it gave back.
function recovery_solver(; failat = Int[])
    starts = Vector{Float64}[]
    solutions = Vector{Float64}[]
    method = JosephsonCircuits.ExternalSolver() do prob, u0
        push!(starts, copy(u0))
        u = copy(u0)
        if length(starts) in failat
            fill!(u, 7.0)
            push!(solutions, copy(u))
            return u, false
        end
        F = similar(u)
        J = copy(prob.jacobian)
        for _ in 1:40
            JosephsonCircuits.hbresidual!(F, prob, u)
            if norm(F) < 1e-12
                push!(solutions, copy(u))
                return u, true
            end
            JosephsonCircuits.hbjacobian!(J, prob, u)
            u .-= J\F
        end
        error("the recovery fixture's Newton solve did not converge")
    end
    return (; method, starts, solutions)
end
