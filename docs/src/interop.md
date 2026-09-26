# Using other solvers

The harmonic-balance problem exposes residuals, derivatives, and linear
operators for external numerical solvers. The core interface requires no
extension. `KrylovJL` is enabled by loading Krylov.jl, and
`SciMLBase.NonlinearProblem(prob)` by loading SciMLBase.

This page first constructs a complete problem. The optional examples below
continue that setup and identify the extra packages they need.

## The problem object

```@example problem
using JosephsonCircuits, LinearAlgebra
circuit = Circuit([
    (:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(1e-12)),
])
wp, Nharmonics = (2pi*4.75001e9,), (8,)
sources = [(mode = (1,), port = 1, current = 0.002e-6)]
prob = hbnonlinearproblem(wp, Nharmonics, sources, circuit;
    assemblejacobian = false)
u = zeros(length(prob))
F = similar(u)
hbresidual!(F, prob, u)
@assert all(isfinite, F)
nothing # hide
```

`assemblejacobian=false` skips the full real-Jacobian assembly plan.
This avoids its storage on a long circuit or a multi-tone grid. Omit it
when the external method needs an assembled Jacobian.

Everything a solver asks for is a method on `prob`:

| Interface | Operation |
|---|---|
| `hbresidual!(F, prob, u)` | `F(u)` |
| `hbjvp!(Jv, prob, u, v)` | `J(u)*v`, matrix free |
| `hbvjp!(Jtw, prob, u, w)` | `transpose(J(u))*w`, matrix free |
| `hbjacobian!(J, prob, u)` | the assembled real Jacobian |
| `hbd2F!` / `hbd3F!` | exact second and third directional derivatives |
| `hbdFdp!(out, prob)` | the derivative with respect to the drive |
| `JacobianOperator(prob, u)` | the Jacobian as a `mul!`-able operator |
| `preconditioner(prob, u)` | the mode coupling preconditioner |

The state uses the **equivalent real representation**. The residual
depends on complex coefficients and their conjugates, so exact derivatives
use this real form. It stores the real and imaginary parts of each complex
mode amplitude side by side, and only the real part of a self-conjugate
mode; `JosephsonCircuits.real_to_complex(u, prob.modelayout.isreal)` and
`JosephsonCircuits.complex_to_real(x, prob.modelayout.isreal)` convert.
A DC drive augments this state with average node voltages. Check
`JosephsonCircuits.isaugmented(prob)` before applying the complex/real
conversion helpers to a problem state.

## A linear solver

`JacobianOperator` implements `size`, `eltype`, `mul!`, `adjoint` and
`transpose`; the preconditioner implements `ldiv!`, `mul!` and `\`. A
`JacobianOperator` freezes its evaluation point at construction. Its
`mul!` reuses the transformed junction derivative. Construct a new operator
after moving `u`. The following example uses Krylov.jl, which must be
installed separately:

```julia
using Krylov
J = JacobianOperator(prob, u)
P = preconditioner(prob, u)
dx, stats = Krylov.gmres(J, -F; N = P, rtol = 1e-6, atol = 0.0)
```

!!! warning "Set the absolute tolerance to zero inside a Newton loop"
    `rtol` is relative to the norm of the right hand side, and in a Newton
    loop that right hand side is the residual being driven to zero. Any
    absolute floor eventually exceeds it, at which point the linear solver
    correctly reports success for a system it never touched, the Newton
    step is zero, and the iteration stagnates with nothing reporting a
    failure.

    Krylov.jl defaults to `atol = sqrt(eps())`, about `1.5e-8`, which is
    sensible standalone but lies above the package's default tolerance of
    `1e-8`: once the residual falls below it, every linear solve returns a
    zero step, and the Newton iteration stalls short of converging.

To use an external Krylov solver for the Newton step of this package's own
solver, rather than writing the loop yourself:

```julia
using Krylov
hbnlsolve(wp, Nharmonics, sources, circuit;
    method = NewtonKrylov(linearsolver = KrylovJL(:fgmres)))
```

`GMRES()` is the default. Only the linear solve changes: the
forcing term, line search, preconditioner escalation and stagnation
handling are untouched.

## A nonlinear solver

Pass a solver object as `method`. `NewtonKrylov()`, `Newton()` and
`QuasiNewton()` are the built-ins; `ExternalSolver(f)` takes a root finder
of your own, which receives the problem and the initial value and returns
`(u, converged)`.

This minimal Newton loop uses the `Krylov` and `LinearAlgebra` imports
above. It has no line search or continuation, so it is an interface
example rather than a robust replacement for the built-in methods.

```julia
mysolver = ExternalSolver() do prob, u0
    u = copy(u0); F = similar(u)
    hbresidual!(F, prob, u)
    # built once, and refactorized at each new point on its structure
    P = preconditioner(prob, u)
    for k in 1:40
        norm(F) <= prob.atol && return (u, true)
        J = JacobianOperator(prob, u)
        d, st = Krylov.gmres(J, -F; N = P, rtol = 1e-10, atol = 0.0)
        st.solved || return (u, false)
        u .+= d
        hbresidual!(F, prob, u)
        JosephsonCircuits.updatepreconditioner!(P, u)
    end
    return (u, norm(F) <= prob.atol)
end

hbnlsolve(wp, Nharmonics, sources, circuit; method = mysolver)
```

`prob.atol` is the tolerance of the solve, which the root is held to
whatever the solver reports.

To use NonlinearSolve.jl, install and load `SciMLBase` and
`NonlinearSolve`. Rebuild the same problem with an assembled Jacobian plan
for a method that uses it:

```julia
using SciMLBase, NonlinearSolve
assembled = hbnonlinearproblem(wp, Nharmonics, sources, circuit)
external_problem = SciMLBase.NonlinearProblem(assembled)
external_solution = NonlinearSolve.solve(external_problem)
```

The adapter provides a Jacobian-vector product and an assembled Jacobian.
Check the external solver's termination status and the residual before
using its root.

## Continuation

`setdrive!(prob, s)` scales the drive in place. If `b0` is the unscaled
source vector, the residual has source term `-s*b0`, so `dF/ds=-b0`.
Changing `s` retains the structural plans, though a preconditioner may
need updating as the state moves.

Increase `s` from zero and carry the converged state to each new drive.
This continuation can reach an operating point that a cold solve misses,
but does not guarantee reaching every branch. The following BifurcationKit
setup also requires Accessors for `@optic`:

```julia
using BifurcationKit, Accessors
branch_residual(u, p) = (Fv = similar(u); drivenresidual!(Fv, prob, u, p.s); Fv)
Jmf(u, p)  = (setdrive!(prob, p.s); dx -> hbjvp!(similar(dx), prob, u, dx))
Jadj(u, p) = (setdrive!(prob, p.s); dw -> hbvjp!(similar(dw), prob, u, dw))

bp = BifurcationProblem(branch_residual, zeros(length(prob)), (s = 0.0,), (@optic _.s);
    J = Jmf, Jᵗ = Jadj,
    d2F = (u,p,v,w)   -> hbd2F!(similar(u), prob, u, v, w),
    d3F = (u,p,v,w,z) -> hbd3F!(similar(u), prob, u, v, w, z))
```

!!! tip "Preconditioning pseudo-arclength continuation"
    Use `BorderingBLS(solver = ls)`, not `MatrixFreeBLS(ls)`. Pseudo
    arclength solves an `(n+1)` by `(n+1)` bordered system; `MatrixFreeBLS`
    attacks that directly, so an `n` by `n` preconditioner is a
    `DimensionMismatch`. `BorderingBLS` decomposes it into two `n`
    dimensional solves, where the preconditioner fits.

    On a large circuit, include a compatible preconditioner when comparing
    continuation methods; an unpreconditioned failure can reflect the
    linear solve rather than the nonlinear branch.

!!! warning "The HB Jacobian does not determine physical stability"
    `∂F/∂u` differentiates an algebraic residual, not the time evolution.
    Its eigenvalues alone do not give growth or decay rates of physical
    perturbations. Do not interpret a continuation library's dynamical
    bifurcation labels as classifications of the pumped circuit.

A singular Jacobian can mark a branch degeneracy, but identifying a fold
requires the relevant continuation conditions. Failure to converge is even
weaker evidence: it can reflect resolution, initialization, or numerical
work limits. `Staged()` reports its search history and possible fold
brackets without proving nonexistence.

Physical stability requires an analysis of perturbation dynamics, such as
Floquet multipliers for a periodic orbit. A singular small-signal operator
can indicate a resonant or threshold condition, but a finite real-frequency
sweep of `hblinsolve` is not a complete stability certificate.

## Sensitivities

[`designsensitivities`](@ref) returns derivatives with respect to named
design parameters, combining all dependent component values by the chain
rule. See the complete [parameter-sensitivity recipe](recipes/sensitivities.md).

For a scattering block, supply derivatives through the `derivatives`
keyword of [`ScatteringParameters`](@ref). These describe the block at
one design point; changing parameters in a circuit dictionary does not
update numbers captured by a block's closure. Rebuild the block's data
and derivative at a new point. The
[line-length example](recipes/scattering-sensitivities.md) shows this interface.

## Testing

The interoperability tests live in their own environment so the main test
suite carries no resolve or precompile cost for Krylov.jl or SciMLBase:

```
julia test/interop/runtests.jl
```
