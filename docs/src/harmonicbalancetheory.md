# Theory and implementation of harmonic balance

This page derives the harmonic-balance equations and explains numerical
choices that affect convergence. The [usage guide](harmonicbalance.md)
covers calls and outputs; [implementation notes](implementation.md) map
the algorithms to source files.

## The equations

Let `Φ` denote physical node flux in webers, so `V=Φ̇`. For linear RLC
elements and sinusoidal Josephson junctions, nodal current balance is

```math
C\ddot\Phi+G\dot\Phi+\Lambda\Phi+
R_J^T\left[I_c\odot\sin(R_J\Phi/\varphi_0)\right]=I_{\mathrm{src}}(t).
```

| Symbol | Meaning |
|---|---|
| `C` | Nodal capacitance matrix, F |
| `G` | Nodal conductance matrix, S |
| `Λ` | Nodal inductive stiffness, inverse henries |
| `R_J` | Incidence from node fluxes to junction branch fluxes |
| `Ic` | Vector of junction critical currents, A |
| `φ₀` | Reduced flux quantum `ħ/(2e)`, Wb |

`⊙` denotes elementwise multiplication. Ground is excluded from the nodal
unknowns. For another current-phase relation, replace the sine by the
element's normalized relation.

Define reduced flux `φ=Φ/φ₀`. Dividing by `φ₀` gives

```math
C\ddot\phi+G\dot\phi+\Lambda\phi+
R_J^T\left[L_j^{-1}\odot\sin(R_J\phi)\right]
=I_{\mathrm{src}}(t)/\varphi_0.
```

Here `Lj⁻¹` is the vector of inverse junction inductances, since
`Ic=φ₀/Lj`. The solver further multiplies the rows by `Lscale`. In HB,
`Lscale=Z0/w0`, using representative geometric-mean port impedance and
drive frequency. The scaled residual is dimensionless. Its tolerance is
not a direct bound on gain error.

### Modified nodal analysis and DC

The displayed equation is the physical starting point. The actual system
includes auxiliary variables where elimination would be singular or
poorly conditioned:

- Mutually coupled inductors retain branch currents. Their constitutive
  rows remain bounded as coupling approaches unity, whereas the eliminated
  inverse-inductance matrix can diverge.
- Scattering blocks retain port currents and hybrid constitutive rows.
- Inductive subnetworks with undetermined static flux receive gauge rows.

A periodic flux has zero mean derivative. DC transport therefore needs
average node voltages in addition to periodic flux coefficients. Resistive
conductances and zero-frequency block models determine those voltages and
currents. This affine augmentation is part of the canonical problem state.
A gauge fixes a static flux offset; it does not prohibit resistive DC
conduction. See the [DC example](recipes/dc.md).

## The mode set and the transforms

For independent pump frequencies `ωp`, expand reduced flux on an integer
lattice:

```math
\phi(t)=\sum_{\boldsymbol m}\phi_{\boldsymbol m}
 e^{i(\boldsymbol m\cdot\boldsymbol\omega_p)t},\qquad
\phi_{-\boldsymbol m}=\overline{\phi_{\boldsymbol m}}.
```

One commensurate fundamental gives a periodic state. Independent
incommensurate pumps give a quasiperiodic state, represented as a periodic
function of the pump phase coordinates.

A real transform along the first grid axis and complex transforms along
the others store nonnegative first-axis indices and both signs on the
remaining axes, with redundant conjugate modes removed. A stored
multi-tone mode can still have negative physical frequency: its frequency
is the dot product of the index tuple with `ωp`.

Retained modes obey harmonic, parity, intermodulation-order, and optional
absolute-frequency limits. Commensurate drives should be supplied as
harmonics of one fundamental. A distinct nonzero mode at zero physical
frequency would duplicate the DC coordinate and is rejected.

### Retained modes versus evaluation grid

The linear term is diagonal in harmonic index, with
`K(w)=-w^2*C+im*w*G+Λ`. A junction nonlinearity couples modes. Its branch
flux is gathered onto the Fourier grid, transformed to phase coordinates,
evaluated pointwise, transformed back, and restricted to retained modes.

The first transform dimension uses `Nt=2Nw-1`, so its highest coefficient
is not a self-conjugate Nyquist coefficient. The default evaluation
harmonic limit is twice the retained limit.

Padding controls aliasing, not the number of unknowns. A cubic product of
a phase bandlimited to order `N` reaches order `3N`. A grid of roughly
`3N` samples can fold those products back into the retained band; a grid
of `4N+1` samples keeps that cubic aliasing out. Higher powers in a sine or
polynomial can require more resolution. Refine retained and evaluation
grids independently.

## The residual and its derivatives

Suppressing the affine DC augmentation, write the residual as

```math
F(x)=B\,g(Ax)+Kx-b.
```

`A` gathers reduced junction phases onto the evaluation grid. `B` applies
the forward transform, current scaling, and incidence map back to nodal
rows. `g` is the current-phase relation, a sine for a Josephson junction.
The exact directional derivative is

```math
J(x)v=B\left[g'(Ax)\odot Av\right]+Kv.
```

A matrix-free product uses two transforms and the linear product. Junction
derivatives at the current state can be cached across Krylov iterations.
Second and third directional derivatives use the corresponding derivatives
of `g` and products of directions.

### Why the Jacobian is real

The retained complex coefficients represent a real waveform. Applying the
nonlinearity couples each coefficient to conjugates of others, so the
residual is not holomorphic. Exact Newton uses real and imaginary parts as
real unknowns, with only a real unknown for a self-conjugate mode.

The assembled real Jacobian and matrix-free product represent the same
derivative. Junction blocks use Fourier coefficients of `g'(φ(t))` at
mode differences. `QuasiNewton()` instead retains a holomorphic
approximation; it is not the exact Newton Jacobian.

## The nonlinear solver

Newton–Krylov solves `J*d=-F` by restarted, right-preconditioned GMRES.
Its forcing tolerance adapts to nonlinear progress, and an Armijo
backtracking line search accepts the step. Work budgets and stagnation
checks stop unsuccessful solves. `solverinfo` records the outcome.

### Preconditioning

| Choice | Retained coupling |
|---|---|
| `BlockDiagonal()` | Independent factorization for each mode |
| `HarmonicBand` | Selected harmonic offsets |
| `MeasuredBand()` | Band selected from the Fourier content of the junction derivative |
| `Clusters` | Mode groups merged according to coupling strength |
| `FullJacobian()` | All real-Jacobian couplings |

`Automatic()` uses pump count and predicted factor memory to choose a
strategy. Sparse or dense node-block factorizations solve the retained
systems. The block method uses graph supernodes without first constructing
a scalar sparse matrix. Mixed-precision factors can be refined in double
precision. See [performance](performance.md) for current defaults.

If GMRES stalls, the solver can expand the coupling or increase factor
precision within its memory budget. Mode-local criteria cannot detect
every nearly singular direction of the coupled operator. Optional Floquet
deflation retains difficult correction directions estimated from the
Arnoldi residual image; it is off by default.

By default, the preconditioner is rebuilt before each Newton step.
`Probe()` compares the cost of extra Arnoldi iterations with a fresh
factorization and may reuse the old one. Structural plans and symbolic
factorization can be reused across compatible value sweeps.

### Source continuation

`Staged()` combines source continuation with a ladder of retained grids.
It solves inexpensive intermediate points loosely, transfers states by
mode tuple, and uses the requested tolerance for the final full-drive
solve. The nonlinear evaluation grid stays fixed across stages.

A failed drive increment can be reduced; a stalled coarse grid can be
grown; a transferred state that fails on a finer grid can retreat in
drive. Different truncations can have different branches and convergence
regions.

A stall after a converged point on the finest grid is recorded as a
possible fold bracket. This describes the search, not a proof of
nonexistence or physical instability. The method can also attempt a
direct full-drive solve from the last converged point. Check stage history
and `solverinfo.converged`. See [continuation](interop.md#Continuation).

## The linearized sweep

About a pumped operating point, a weak signal at `ωs` couples to modes at
`ωm=ωs+m⋅ωp`. Before augmentation and signed-frequency representation
choices, its operator consists of junction modulation and the linear
circuit terms:

```math
A(\omega_s)=A_J+\Lambda+i\,\mathrm{diag}(\omega_m)G
-\mathrm{diag}(\omega_m^2)C.
```

`A_J` is assembled from the Fourier coefficients of the junction derivative
at the pump solution. They produce signal-to-idler coupling. The sparse
pattern is analyzed once; each signal frequency updates and refactorizes
its values.

Sources excite individual port modes. For a real reference impedance,
port-wave normalization includes `1/sqrt(abs(w)*Z)`, so scattering
coefficients measure photon flux. There is no such wave normalization at
zero frequency. See [conventions](conventions.md#Photon-gain-and-power-gain).

### Noise and adjoint solves

An adjoint driven at an observed port mode gives its response to a source
anywhere in the circuit. Contracting with the noise-source stamps avoids
a full solve per internal bath. The transpose must preserve the signed
frequency convention, including conjugate-idler commutator signs.

Each internal noise channel carries its symmetrized noise `nbar + 1/2`,
half a photon in the vacuum, and so does each port mode, at the
temperature of its termination. `QE` compares a selected signal
contribution with the total output noise, every input in its state, and
`nbar` is that noise less the vacuum's half photon. `CM` includes signed input and internal-bath
contributions and is `+1` or `-1` for a complete physical output mode.
`Cnoise` contains added noise, not total output covariance.

## Scattering blocks in the frequency domain

At its reference impedance matrix `R`, a frequency-preserving block uses

```math
(I-S)R^{-1/2}v-(I+S)R^{1/2}i=0.
```

Port currents enter nodal balances as auxiliary unknowns. No inverse of
`I±S` is required, so ideal opens, shorts, and throughs need no singular
admittance conversion. Analysis-impedance conversion happens at the
port-wave boundary.

DC rows use the stated zero-frequency limit. Passive loss emits the
covariance `(nbar + 1/2)(I-S*S')`, of commutator `I-S*S'`; an active block
states compatible noise. Pumped blocks couple
harmonic modes as described in the [scattering guide](scattering.md#Pumped-devices).

## Sensitivities

A parameter can change both the linearized circuit and the pump state.
At fixed pump, differentiate component stamps and contract with forward
and adjoint solutions. The pump contribution follows

```math
J\frac{dx}{dp}=-\frac{\partial F}{\partial p},
```

with the exact real Jacobian at a converged root. Near singularity, the
sensitivity can be large; its accuracy depends on the operating point and
linear solves.

The implementation chooses forward contractions per component or reverse
contractions per output functional according to their relative counts.
Design derivatives combine component-value derivatives by the chain rule.
See the [verified example](recipes/sensitivities.md).

## Validation

Tests compare matrix-free products with assembled Jacobians, auxiliary
current formulations with their Schur complements, noise sources with
equivalent ports, blocks with explicit circuits, and sensitivities with
finite differences. Device examples compare with WRspice and the published
ADS comparisons listed on the [home page](index.md#References).

See [numerical references](numerical-references.md) for algorithm sources
and [implementation notes](implementation.md) for source-code organization.
