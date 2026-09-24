# Using other solvers

The nonlinear system this package solves is available as an object, so a
solver it has never heard of can drive it. The problem object needs no
package extension: its interface is `mul!`, `ldiv!` and a handful of
in-place functions. Two of the calls below do: [`KrylovJL`](@ref) needs
Krylov.jl loaded, and `SciMLBase.NonlinearProblem(prob)` SciMLBase.

## The problem object

```julia
prob = hbnonlinearproblem(wp, Nharmonics, sources, circuit, circuitdefs)
```

builds the harmonic balance system without solving it. Pass
`assemblejacobian = false` for a matrix-free solver, which skips building
the real Jacobian plan. That plan holds the structure of the whole real
Jacobian, every pair of coupled modes at every junction, which on a long
line or a multi-tone grid is the largest object of the solve.

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
needs the real form. It holds the real and imaginary parts of each complex
mode amplitude side by side, and only the real part of a self conjugate
mode; `JosephsonCircuits.real_to_complex(u, prob.modelayout.isreal)` and
`JosephsonCircuits.complex_to_real(x, prob.modelayout.isreal)` convert.
When the circuit injects direct current the unknowns are the canonical
state instead, which carries the average node voltages as well, and
`JosephsonCircuits.isaugmented(prob)` says so.

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
    sensible standalone but lies above the package's default tolerance of
    `1e-8`: once the residual falls below it, every linear solve returns a
    zero step, and the Newton iteration stalls short of converging.

To use an external Krylov solver for the Newton step of this package's own
solver, rather than writing the loop yourself:

```julia
using Krylov
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

hbnlsolve(wp, Nharmonics, sources, circuit, circuitdefs; method = mysolver)
```

`prob.atol` is the tolerance of the solve, which the root is held to
whatever the solver reports.

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

    This matters at scale: on a travelling wave amplifier with no Jacobian
    assembled, preconditioned matrix-free continuation reaches full drive,
    and without a preconditioner it fails to compute even the initial
    tangent.

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

## Sensitivities

[`designsensitivities`](@ref) differentiates the scattering parameters
with respect to the design parameters a circuit's values are written in
terms of, rather than with respect to its components. Every component
whose value depends on a parameter contributes to that parameter's
derivative, with the exact derivative of its value, and the result has
one slot per parameter.

```julia
out, dSdp = designsensitivities(circuit, circuitdefs, ws, wp, sources,
    Nmodulationharmonics, Npumpharmonics; parameters = [:Ic, :Cg])
dSdp    # (outputmode, outputport, inputmode, inputport, parameter, freqindex)
```

A scattering block has no component value to differentiate, so it states
its derivative with respect to a parameter itself, through the
`derivatives` keyword of [`ScatteringParameters`](@ref), and its
contribution lands in the same slot as the lumped components of that
parameter.

## Testing

The interoperability tests live in their own environment so the main test
suite carries no resolve or precompile cost for Krylov.jl or SciMLBase:

```
julia test/interop/runtests.jl
```
