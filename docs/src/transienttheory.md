# Theory and implementation of the time domain solver

The transient solver integrates the same circuit equations used by
harmonic balance. Its tangent and adjoint differentiate the recorded
discrete trajectory. The [usage guide](transient.md) covers simulation
choices; the [noise guide](transientnoise.md) explains temporal modes and
physical baths.

## The equations

For physical node flux `Φ` in webers, the lumped circuit equation is

```math
C\ddot\Phi+G\dot\Phi+\Lambda\Phi+
R_J^T\left[I_c\odot\sin(R_J\Phi/\varphi_0)\right]=I_{\mathrm{src}}(t).
```

`C`, `G`, and `Λ` are the physical capacitance, conductance, and inductive
stiffness matrices. `R_J` maps node flux to junction branch flux;
`φ₀=ħ/(2e)` and `Ic=φ₀/Lj`. Other current-phase relations replace the sine.
This is the same notation as the [HB derivation](harmonicbalancetheory.md#The-equations).

The implementation uses reduced flux `φ=Φ/φ₀`, auxiliary branch currents,
and scaled rows. Let `x` denote that augmented state and write its scaled
lumped equations as

```math
\widetilde C\ddot x+\widetilde G\dot x+\widetilde Kx+j(x)=b(t).
```

The tildes distinguish these matrices from physical nodal matrices.
The row scaling is equivalent to multiplying physical current balance by
`Lscale/φ₀`. `Lscale` is a problem-specific inductance scale derived from
the circuit and kept fixed across steps and continuations. Auxiliary
currents use scale `Lscale*i/φ₀`.

Mutually coupled inductors retain branch currents, avoiding inversion of a
nearly singular inductance matrix. Gauge rows fix undetermined static flux
offsets. Scattering port currents and rational-block or line states extend
the equations as described below.

### Differential and algebraic directions

The capacitance matrix can be singular. The solver does not invert it or
add artificial shunts to make it invertible. Instead, it identifies the
directions constrained by the algebraic equations.

| Direction | Local equation | Example |
|---|---|---|
| Nonzero capacitance | Differential flux and voltage dynamics | A shunt LC resonator |
| No capacitance, with conductance | The equation determines a rate | A resistively terminated node without shunt capacitance |
| Neither capacitance nor conductance | The equation constrains flux | An inductive or junction direction driven by current |

Capacitor-free islands and auxiliary rows must be considered together.
A block's port-current equations can constrain voltages even when a graph
of its terminals alone would suggest otherwise. The solver derives the
constraint directions from the rate equations, rather than guessing their
rank from connectivity.

The initial state must satisfy these constraints. For a flux constraint,
the initial rate must also satisfy its differentiated equation. Full
nominal order applies to differential variables; algebraic output rates
can have lower order and need their own convergence checks.

## The stepping rules

The solver supports trapezoidal integration, two-stage Gauss–Legendre
collocation, and backward Euler. Their nominal differential orders are
two, four, and one. Trapezoidal and Gauss–Legendre preserve linear
lossless oscillations without numerical damping; backward Euler damps
them. Symplecticity of Gauss collocation applies to the unconstrained
Hamiltonian problem, not automatically to every projected circuit DAE.

### Solve for increments

A trapezoidal step solves for `d=x_{n+1}-x_n`:

```math
\left[\frac{4}{h^2}\widetilde C+
\frac{2}{h}\widetilde G+\widetilde K\right]d
+j(x_n+d)=r_n.
```

The right-hand side contains the sources, the previous rate, and the old
stiffness and nonlinear-current contributions. Formulating the step in
increments avoids a large term proportional to `C*x_n` when flux
accumulates without contributing to capacitor current. The residual
tolerance then scales with the currents being balanced.

Newton reuses a factorization while corrections reduce the residual
sufficiently. It refreshes when the contraction deteriorates. A linear
circuit can reuse one factorization; nonlinear circuits refresh according
to their trajectory. A failed step throws rather than silently continuing.

### Gauss–Legendre stages

The two-stage tableau and nodes are

```math
A=\begin{pmatrix}
1/4&1/4-\sqrt3/6\\
1/4+\sqrt3/6&1/4
\end{pmatrix},\qquad
c=\begin{pmatrix}1/2-\sqrt3/6\\1/2+\sqrt3/6\end{pmatrix}.
```

With step size `h`, stage increments `d_i`, and initial rate `v_n`, the
scaled lumped stage equations are

```math
\sum_l\frac{(A^{-2})_{il}}{h^2}\widetilde C d_l+
\sum_l\frac{(A^{-1})_{il}}{h}\widetilde G d_l+
\widetilde K(x_n+d_i)+j(x_n+d_i)
=b(t_n+c_i h)+\frac{(A^{-1}\mathbf1)_i}{h}\widetilde C v_n.
```

The **inverse** tableau has eigenvalues `μ=3±i√3`; `A` itself has
eigenvalues `1/4±i√3/12`. Freezing the junction stiffness at the mean of
the two stages decouples the approximate Newton operator into a complex
system and its conjugate:

```math
(\mu/h)^2\widetilde C+(\mu/h)\widetilde G+
\widetilde K+J_*.
```

Only one complex factorization is required. The iteration still evaluates
the true stage residual, so the frozen stiffness is an iteration device,
not a replacement for the nonlinear equations. It is refreshed when
corrections no longer contract adequately.

The endpoint is obtained from the collocation polynomial, then projected
where constraints require it. The next step extrapolates the polynomial
for its predictor. Tolerances are checked per equation and per condition,
including a roundoff floor based on the magnitudes of the balanced terms.
A small row is not judged solely against the largest drive elsewhere.

### Source values at stages

A callable drive is evaluated at the stage times. A tangent current given
only on the recorded grid is interpolated with a cubic Lagrange stencil
when enough points exist, or a line on very short records. Inputs supplied
at both grid and stage times use those values directly. Pulsed gain and
noise use this form to preserve probe support and sinusoidal phase.

### The projection of the endpoint

A collocation endpoint need not exactly satisfy a nonlinear algebraic
constraint. The solver projects it to the roundoff of the terms in that
constraint. Trapezoidal and backward Euler endpoints are also projected
where their Newton tolerance leaves a constraint residual.

Let `Z_R` span constrained state directions and `Z_L^T` select the
corresponding constraint equations. A block can make the operator
unsymmetric, so these are not interchangeable. The correction `Z_R*α`
solves

```math
Z_L^T[\widetilde K(x+Z_R\alpha)+j(x+Z_R\alpha)-b]=0,
```

with the small Newton matrix
`Z_L^T*(K̃+J')*Z_R`. Linear invariants that the steps already preserve need
no repeated nonlinear correction. Index-one rates and block currents are
read from the rate equations on their determined range.

For other algebraic directions, reported rates come from differentiated
constraints. The source derivative is evaluated by a small central
difference. Where a block couples its moving internal state into the
constraint, Gauss–Legendre reads the projected-direction rate from the
cubic through the initial state, two stages, and endpoint. Its rate
accuracy is third order in that case. Tangent and adjoint calculations
differentiate the projection and rate-reading operations as well.

## Scattering blocks

A constant real scattering block uses the hybrid relation

```math
(I-S)R^{-1/2}v-(I+S)R^{1/2}i=0.
```

Its port currents are auxiliary unknowns in the nodal balances. This
avoids converting an open, short, or through into a singular admittance.
The block is evaluated at its own reference impedances; shared nodes
supply the loading from the surrounding circuit. The resulting operator
can be unsymmetric, so adjoints use transposed solves.

### Rational internal states

For `S(s)=D+C_b*(sI-A_b)^(-1)*B_b`, denote internal states by `z` and the
incident waves at the two stages by `U`. Their stage values satisfy

```math
Z=Z_z z_n+Z_u U,\qquad
Z_z=(I-hA\otimes A_b)^{-1}(\mathbf1\otimes I),\qquad
Z_u=(I-hA\otimes A_b)^{-1}(hA\otimes B_b).
```

The outgoing stages are `(I₂⊗D)U+(I₂⊗C_b)Z`. The Gauss weights are one
half, giving

```math
z_{n+1}=E_z z_n+E_u U,\qquad
E_z=I+\frac h2(A_b,A_b)Z_z,\qquad
E_u=\frac h2[(A_b,A_b)Z_u+(B_b,B_b)].
```

These are algebraic eliminations of the discrete block stages, not a
separate time integrator. In the complex stage basis the frozen operator
contains the block evaluated at `s=μ/h`. The block state is retained in
continuations, tangent responses, and adjoint initial-state derivatives.

For a fitted pumped block, harmonic filters are followed by periodic
modulation. The within-step modulation is handled as a correction to the
unconverted stage operator. Its algebra and transpose are included in the
responses. See [pumped blocks](scattering.md#Pumped-devices).

## Transmission lines

An ideal line uses traveling waves and an explicit delay. At a port, with
voltage `v`, line impedance `Z`, and incoming wave `q`,

```math
i=v/Z-2q/\sqrt Z,\qquad p=v/\sqrt Z-q.
```

The outgoing wave `p` becomes the incoming wave at the opposite port after
the propagation delay. This method of characteristics avoids a nodal
admittance singularity at a half-wave resonance: resonant round trips
emerge through repeated delayed propagation.

### History interpolation and accuracy

A delay must be at least one step so all history read by a step has
already been accepted. The implementation uses centered Lagrange
interpolation with as many neighboring samples as are available, up to
three on each side: linear, cubic, or quintic.

The chosen centered stencils do not amplify a Fourier component on a
round trip. A one-sided high-order stencil does not share that property
and can destabilize a mismatched lossless line. At least four steps per
delay allow the stencil needed to retain fourth-order accuracy. Shorter
delays use lower-order interpolation.

Interpolation spreads an edge over its stencil, so a small numerical
precursor can occur before the physical arrival. Refine the step when
arrival timing or weak precursors matter. The physical line is lossless;
the finite-step interpolation still has numerical dispersion and error.

### State and responses

The history is stored in a ring sized by the longest delay and stencil
reach. Checkpoints retain the history needed to replay each interval.
`transientstate` can initialize DC line currents and corresponding waves;
continuing a solution preserves the full recent history.

The tangent reads perturbed waves with the same interpolation. The adjoint
scatters through the transpose of those reads and returns sensitivity to
the prehistory. Thus line memory is part of both the state and the noise
initialization, not just a delay applied to the final port trace.

## Fitting measured data

The fitter uses common poles across scattering entries. Relaxed vector
fitting relocates them through a shared least-squares problem, then fits
residues and the direct term. Per-entry QR reductions eliminate local
unknowns before the shared solve, avoiding a dense matrix spanning every
entry and sample. See [Gustavsen and Semlyen, Gustavsen, and Deschrijver et al.](numerical-references.md#Rational-fitting).

Pole pruning and merging are accepted only when the refitted response
stays within the allowed error. A state-space realization uses the ranks
of the residue matrices. Explicit delay removal can reduce the rational
order needed for a long cable or traveling-wave device.

The search for the fewest poles meeting a tolerance starts where the
samples allow: a fit of `N` poles and a constant combines `N + 1` real
functions of frequency, so its error is at least what the samples leave
beyond their best approximation of that rank. An order is made passive
only where its poles could meet the tolerance: the passivity enforcement
keeps the poles, and a least squares at them, reweighted toward the
samples it misses most as in Lawson's algorithm, bounds the error of
every fit with them from below. The scan stops on measured progress,
once the order has doubled since the fit's own error last fell by a set
fraction of itself, which bounds the work on data of a high degree. It
also stops at a fit with too many states for memory to hold the dense
matrices of its realization and its passivity enforcement. See
[Eckart and Young, and Rice and Usow](numerical-references.md#Rational-fitting).

### Passivity and fit accuracy

For a passive model, the largest singular value of `S(iw)` must not exceed
one. Enforcement perturbs residues and the direct term at violation
points, found on a frequency grid and by a sweep of the whole frequency
axis which bounds the fit's dissipation `I-S'S` from below over each
interval, from its value and derivative at the interval's centre and the
distances of the poles, so that no band passes between samples. The
sweep establishes a fit passive to its tolerance, and tests a fit which
is not enforced as well.

A realization supplied directly has no residues to bound, and is tested
by a search for the peaks of its largest singular value and by the
crossings of the level `sqrt(1 + atol)`, at which its dissipation `I-S'S`
reaches `-atol`, found through the Hamiltonian matrix where the direct
term stands well under the level, and otherwise through a matrix pencil,
which needs no inverse of `I-D'D`. A peak within the roundoff of the
level can leave the test undecided, and the block is then accepted; see
[`passivityassessment`](@ref JosephsonCircuits.passivityassessment) for its verdicts and
[`PassivityEnforcement`](@ref) for margins. A declared lossless model is
checked using both the model and its inverse; small apparent loss is not
a reason to discard noise channels.

Active blocks with an explicit noise covariance follow that noise
contract rather than passive enforcement. In all cases, stable poles and
physical noise do not establish accuracy outside the fitted band. Validate
scattering, small loss, and relevant time responses against the data or an
independent circuit model.

## The tangent and the adjoint

The tangent differentiates the stage equations about the recorded
junction phases, including block states, line interpolation, endpoint
projection, and output feedthrough. Iteration on the step's frozen
factorization solves the differentiated equations to the response tolerance.

The adjoint applies the transposes of the same operations in reverse time.
Its initial-state derivatives include flux, rate, block states, and line
prehistory. On an unsymmetric block circuit, using the original rather
than transposed operator would be incorrect.

Checkpoint replay reconstructs each interval from its saved state,
predictor, block states, and line history. Factorizations restart at the
same boundaries in the original and replayed solve. A replay that fails
to reach the next checkpoint to roundoff is an error. This controls state
history storage; it does not eliminate every array proportional to the
record length.

## Quantum noise

Bath kernels are adjoint derivatives of measured quadratures with respect
to physical noise sources. Contracting them with sinusoidal bath
quadratures gives covariance and commutator contributions. A stationary
initial term includes fluctuations already stored at the start.

For trapezoidal and backward Euler, the stationary operator uses discrete
rates `i*(2/h)*tan(w*h/2)` and `(1-exp(-i*w*h))/h` respectively.
Gauss–Legendre uses the continuous stationary response. At resolved
frequencies its rate mismatch is fourth order, but near Nyquist the
continuous and discrete responses can differ substantially. This is why
cutoff, grid spacing, and time step require separate convergence checks.

Rational and line states participate in the stationary solve. Correlated
block noise is contracted as a covariance, including complex off-diagonal
terms with the measurement's phasor convention. Pumped blocks require
joint frequency ladders; see [implementation notes](implementation.md#Pumped-block-noise).

## Execution on a device

See [performance](performance.md#GPU-execution) for the user-facing
choices and [implementation notes](implementation.md#Device-execution)
for kernels, synchronization, and host-assisted projection.

## Validation

The suites compare analytic RC/LC/junction limits, time-step orders,
stationary scattering and noise against HB, explicit networks against
blocks, finite-difference and tangent derivatives, and tangent/adjoint dot
products. WRspice comparisons are documented in the
[JPA recipe](recipes/transient-wrspice.md). GPU tests compare the device
paths with their CPU counterparts.

See [numerical references](numerical-references.md) for Gauss collocation,
vector fitting, norm calculations, and traveling-wave formulations.
