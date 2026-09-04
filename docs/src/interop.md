# Using other solvers

The nonlinear system this package solves is available as an object, so a
solver it has never heard of can drive it. Nothing on this page needs a
package extension: the interface is `mul!`, `ldiv!` and a handful of
in-place functions.

## The problem object

```julia
prob = hbnonlinearproblem(wp, Nharmonics, sources, circuit, circuitdefs)
```

builds the harmonic balance system without solving it. Pass
`assemblejacobian = false` for a matrix-free solver, which skips building
the real Jacobian plan. On a multi-tone problem that plan is the largest
object in the solve; on a 1024 cell travelling wave amplifier with twenty
pump harmonics, building without it took 0.03 s against 0.27 s with, and
avoided allocating a matrix of 1.26 million entries.

Everything a solver asks for is a method on `prob`:

| | |
|---|---|
| `hbresidual!(F, prob, u)` | `F(u)` |
| `hbjvp!(Jv, prob, u, v)` | `J(u)*v`, matrix free |
| `hbvjp!(Jtw, prob, u, w)` | `transpose(J(u))*w`, matrix free |
| `hbjacobian!(J, prob, u)` | the assembled real Jacobian |
| `hbd2F!` / `hbd3F!` | exact second and third directional derivatives |
| `hbdFdp!(out, prob)` | the derivative with respect to the drive |
| `JacobianOperator(prob, u)` | the Jacobian as a `mul!`-able operator |
| `preconditioner(prob, u)` | the mode coupling preconditioner |

The state is the **equivalent real representation**. The harmonic balance
residual is not complex differentiable, so the implicit function theorem
does not hold with the holomorphic Jacobian; anything relying on a Jacobian
needs the real form. `tocomplex` and `toreal` convert.

## A linear solver

`JacobianOperator` implements `size`, `eltype`, `mul!`, `adjoint` and
`transpose`; the preconditioner implements `ldiv!`, `mul!` and `\`. A
`JacobianOperator` freezes its evaluation point at construction -- one
forward transform -- and its `mul!` pays only the product, so inside a
Krylov loop nothing is recomputed per iteration; construct a new operator
after moving `u` (construction is the point update). Between
them that covers Krylov.jl, IterativeSolvers.jl, KrylovKit.jl,
LinearSolve.jl and LinearMaps.jl.

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
    sensible standalone but causes issues here. On a JPA: with the default, 40
    Newton iterations and a final residual of `3.2e-10`, never converged;
    with `atol = 0.0`, 7 iterations and `7.4e-17`.

To use an external Krylov solver for the Newton step of this package's own
solver, rather than writing the loop yourself:

```julia
hbsolve(ws, wp, sources, Nmod, Npump, circuit, circuitdefs;
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

```julia
mysolver = ExternalSolver() do prob, u0
    u = copy(u0); F = similar(u)
    hbresidual!(F, prob, u)
    for k in 1:40
        J = JacobianOperator(prob, u)
        P = preconditioner(prob, u)
        d, st = Krylov.gmres(J, -F; N = P, rtol = 1e-10, atol = 0.0)
        st.solved || return (u, false)
        u .+= d
        hbresidual!(F, prob, u)
    end
    return (u, norm(F) < tol)
end

hbnlsolve(wp, Nharmonics, sources, circuit, circuitdefs; method = mysolver)
```

With NonlinearSolve.jl, `SciMLBase.NonlinearProblem(prob)` builds a
`NonlinearFunction` carrying the matrix-free product as `jvp` and the
assembled Jacobian as `jac`.

## Continuation

`setdrive!(prob, s)` scales the drive in place. The residual is
`B(sin(A*u)) + K*u - b` and the drive enters only through `b`, so this is
the one parameter that can be varied without touching sparsity, plans or
the preconditioner, and `dF/ds = -b` is exact.

Stepping `s` from zero and carrying the converged state forward walks onto
the driven branch, which is the reliable way to reach an operating point a
cold solve cannot find. With BifurcationKit:

```julia
F(u, p) = (Fv = similar(u); drivenresidual!(Fv, prob, u, p.s); Fv)
Jmf(u, p)  = (setdrive!(prob, p.s); dx -> hbjvp!(similar(dx), prob, u, dx))
Jadj(u, p) = (setdrive!(prob, p.s); dw -> hbvjp!(similar(dw), prob, u, dw))

bp = BifurcationProblem(F, zeros(length(prob)), (s = 0.0,), (@optic _.s);
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

    This matters at scale. On a 2560 unknown travelling wave amplifier with
    no Jacobian assembled, preconditioned matrix-free continuation reaches
    full drive in twelve steps; the same run without a preconditioner fails
    to compute even the initial tangent.

!!! danger "Eigenvalues of the harmonic balance Jacobian are not stability"
    `∂F/∂u` is the derivative of an *algebraic* residual, not a linearized
    flow, so its eigenvalues have no direct physical meaning.

    Folds are real and useful: `∂F/∂u` becoming singular is exactly the
    turning point where the solution branch loses existence, which is the
    bistability and hysteresis of a driven parametric amplifier.

    Anything a continuation library labels a Hopf bifurcation here is a
    numerical artifact. The physical instability of a pumped operating
    point is a Neimark-Sacker bifurcation of the underlying periodic orbit,
    a Floquet question which the linearized system of [`hblinsolve`](@ref)
    answers in a different form: it appears as a signal frequency at which
    the linearized system matrix becomes singular, which is the parametric
    oscillation threshold.

## Optimization

`hbobjective` is the closure a gradient based optimizer calls: a vector of
physical design parameters to a gain vector and its exact derivative.

```julia
Ic, adens, Ccoup = @params Ic adens Ccoup
circuit = [("P1","1","0",1), ("R1","1","0",50.0),
           ("C1","1","2", Ccoup),
           ("Lj1","2","0", phi0/Ic),
           ("C2","2","0", Ic*adens)]

obj = hbobjective(ws, wp, sources, Nmod, Npump, circuit,
                  [Ic, adens, Ccoup], Dict())
objectivevalue(obj, p)      # gain in dB at each frequency
objectivejacobian(obj, p)   # its derivative, one column per parameter
```

The gradient comes from the adjoint, so its cost grows far more slowly with
the number of parameters than differencing does.

The objective memoizes on `p`, because an interior point solver evaluates
the value and the Jacobian at the same point in sequence and would
otherwise pay two harmonic balance solves per iteration.

A trial point which does not converge raises `HBConvergenceFailure` rather
than returning a state that looks like a solution.

!!! note "Wiring this into JuMP"
    Use `MOI.VectorNonlinearOracle`, not `@operator`. Operators take scalar
    arguments and return scalars, so a gain target at `m` frequencies would
    need `m` operators over the same expensive solve; the oracle registers
    the whole vector as one constraint block with a shared Jacobian, which
    is what one solve plus one adjoint produces.

    Scale the design variables to order one. Ipopt relaxes variable bounds
    slightly before solving, by an amount which does not shrink with the
    bounds themselves, so bounds of order `1e-13` are not enforced in any
    useful sense: optimizing a coupling capacitance directly in farads
    produced a converged solution with a *negative* capacitance, well
    outside the bounds given. Optimizing in units of the starting point,
    with the Jacobian chain ruled accordingly, fixes it.

    Signal a failed evaluation by filling the output with `NaN`. JuMP does
    not catch exceptions thrown from an oracle callback, so throwing aborts
    the whole optimization; `NaN` is the signal Ipopt understands and it
    cuts the trial step.

    Set `hessian_approximation = "limited-memory"`. The adjoint gives first
    derivatives only.

## Sensitivities

`sensitivityparameters` names physical parameters rather than components.
Every component whose value depends on a parameter contributes to that
parameter's derivative, and the result has one slot per parameter.

```julia
sol = hblinsolve(ws, circuit, circuitdefs;
    sensitivityparameters = [Ic, Cg], returnSsensitivity = true)
sol.Ssensitivity     # (out, in, parameter, frequency)
```

A scattering block has no component value to differentiate, so its
contribution arrives as an explicit `dS/dθ` keyed by `(component, parameter)`
in `blockjacobians`, and lands in the same slot as the lumped components of
that parameter.

## Testing

The interoperability tests live in their own environment so the main test
suite carries no resolve or precompile cost for Krylov.jl or SciMLBase:

```
julia test/interop/runtests.jl
```
