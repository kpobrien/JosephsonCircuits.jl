# Theory and implementation of harmonic balance

Harmonic balance finds a periodic steady state as the coefficients of
its Fourier series: the circuit's equations, with the state written as a
sum over a finite set of modes, become a nonlinear algebraic system in
those coefficients, whose residual is evaluated by transforming the
junction fluxes to the time domain, applying the sine there, and
transforming back. This page describes the equations, the mode set, the
residual and its derivatives, the nonlinear solver, the linearized
sweep with its noise and sensitivities, the scattering blocks in the
frequency domain, and what runs on a device. The [usage page](harmonicbalance.md)
shows the calls.

The code lives in `src/circuit/` for the compiled circuit, its matrices
and the modified nodal analysis augmentation; `src/harmonics/` for the
mode set, the transforms, the residual and the Jacobian plans;
`src/solvers/` for the nonlinear solvers, their preconditioners and
factorizations; and `src/linearized/` for the sweep, its outputs, the
noise, the sensitivities and the device paths.

## The equations

The unknowns are the node fluxes `phi`, the integrals of the node
voltages, so that a capacitor, a resistor and an inductor to ground
contribute `C phi''`, `G phi'` and `phi/L` to a node's current balance
and a junction the current `I_c sin(phi_b)` of its branch flux in units
of the reduced flux quantum `phi0`. The equations are the Kirchhoff
current laws of every node but ground,

```math
C\ddot\phi + G\dot\phi + L^{-1}\phi + R_b^T I_c \sin(R_b \phi/\varphi_0) = I(t),
```

with `C`, `G` and `L^{-1}` the capacitance, conductance and inverse
inductance matrices of [`numericmatrices`](@ref), `R_b` the incidence
matrix of the junction branches, and `I(t)` the currents the sources
inject. For a periodic state with modes at the frequencies `w_m` the
linear part is diagonal in the modes: the mode `m` of the state sees
`K(w_m) = -w_m^2 C + i w_m G + L^{-1}`, and the junction term couples
every mode to every other.

The formulation is a modified nodal analysis. A mutually coupled
inductor keeps its branch current as an auxiliary unknown with its
constitutive equation `R_b phi - (L_b/Lscale) u = 0`, whose entries stay
bounded as the coupling coefficient approaches one, where the nodal
inverse inductance diverges as `1/(1 - k^2)`; a scattering block keeps
its port currents as auxiliary unknowns with the hybrid equation of the
block; and a subnetwork no inductor or junction connects to ground,
whose static flux the equations do not determine, has it fixed to a
reference by one gauge row per subnetwork and zero frequency mode. The
rows are scaled by `Lscale/phi0` with `Lscale = Z0/w0`, the geometric
mean port impedance over the geometric mean drive frequency, so that the
entries are of order one for a circuit driven near its characteristic
impedance, the auxiliary currents are comparable to the fluxes, and the
residual tolerance means the same in any unit system.

The zero frequency mode of the periodic state is the static flux, and a
voltage is its derivative, so the periodic state carries no average
voltage: a resistor is open there and a block sees no current. With a
direct current bias the average node voltages are carried as unknowns
beside the periodic state, with the resistors' conductances as their
equations and their currents coupled into the nodal rows, and the
transport rows of the zero frequency blocks beside them; that
augmentation is affine and rides on the harmonic system in a canonical
layout, so a solver written against the problem object does not see
whether the circuit has direct current in it (`src/harmonics/directcurrent.jl`).

## The mode set and the transforms

With `N` pumps the modes are tuples of `N` harmonic indices on a lattice,
the frequency of a mode the dot product of its indices with the pump
frequencies. The retained set keeps the indices up to `Nharmonics` per
pump, the odd or the even sums the mixing process couples through, the
zero mode on request, the modes within an intermodulation order, and the
modes within a frequency window; a pair of pumps whose products reach
zero frequency is refused, since such a mode would duplicate the direct
current coordinate. The state is real, so only the positive frequency
half of the lattice is stored: the first dimension holds the harmonics
of the first pump from zero up, and the others both signs.

The nonlinearity is evaluated on a time domain grid by the real inverse
transform of the junction fluxes, the sine applied at every point, and
the forward transform back, `Nt = 2 Nw - 1` points in the first
dimension so that the highest mode is not the Nyquist point, whose
coefficient a real transform forces real. The grid is
`Nevaluationharmonics`, twice the retained set by default: a grid of `M`
harmonics folds a product of order `p` back to `p - (2M + 1)`, products
of two retained modes reach twice the retained set and fold out of it on
a grid half again as large, but the leading nonlinearity of a junction is
cubic and its products reach three times the set, which the three halves
grid folds onto the tone itself. Twice the set dealiases the cubic
products; only the transforms grow, the unknowns are the same. Measured
on a 64 junction line with three tones retaining `(6, 4, 4)` against a
grid four times the set, the unpadded grid is off by 1.5e-4 of the
strongest mode in the fourth order modes, the three halves grid by
1.8e-7, and twice by 2.7e-11.

## The residual and its derivatives

The residual of the harmonic balance system is

```math
F(x) = B\,\sin(A x) + K x - b,
```

with `A` the map from the unknowns to the junction fluxes on the time
grid, the incidence matrix, the gather onto the junction branches and
the inverse transform folded together, `B` the map back, the forward
transform, the `Lscale/L_j` scaling and the transposed incidence, `K`
the linear term with the mode frequencies substituted, and `b` the
sources. Both maps are precomputed index plans, one work item per output
slot reading only index maps, with no scatter and no atomics, so they run
unchanged on a device, and one plan serves every entry point: the
residual applies the sine between them, the Jacobian-vector product
`J v = B (cos(A x) .* A v) + K v` the cosine, and the Hessian-vector
product the negative sine. A product costs two transforms and the linear
term, with the time domain cosine cached across products at one point.

The residual is not complex differentiable, since the sine of a real
flux couples a mode to the conjugates of the others, so the solvers
work in the equivalent real representation, the real and imaginary
parts of every mode as real unknowns, where the Jacobian is exact and
the implicit function theorem applies. The exact real Jacobian, when
assembled, is the incidence triple product of one dense mode block per
junction, from the Fourier coefficients of `cos(phi(t))` on the aliased
mode differences, deposited at the four node pairs of each junction;
the assembly reads that table backwards from each output entry, so it
stores no gather. The holomorphic Jacobian of the quasi-Newton method
keeps the truncated differences and is an approximation.

## The nonlinear solver

The default is Jacobian-free Newton-Krylov: each Newton step solves
`J d = -F` by restarted GMRES over the matrix-free product, right
preconditioned, to the Eisenstat-Walker forcing tolerance clamped to
`[1e-10, 0.9]`, and takes the step through an Armijo line search with
halving and safeguards. The solve ends promptly when it cannot succeed
and says why: the Newton steps spent, the Arnoldi steps beyond the work
budget of one restart length per step, a line search that finds no
decrease after a rebuilt preconditioner, or a residual history that
projects no convergence within the budget after a recovery of exact
Newton steps. A stall outside the Newton basin is the continuation
problem the staged method exists for.

The preconditioner is the Jacobian with its mode coupling restricted to
a selected set and reduced to the mode diagonal elsewhere, factorized.
The block diagonal is one small factorization per mode; a harmonic band
keeps the couplings within a number of offsets per pump, the restriction
the Toeplitz structure of the junction term asks for, and a measured
band sets its width from the Fourier coefficients of the cosine;
clusters of modes are merged in decreasing coupling strength until what
is left between them is contractive; and the full Jacobian is an exact
solve. Any of them is factorized either by the backend's sparse solver
or by the block factorization, which eliminates the circuit graph as
supernodes of dense mode blocks ordered by KLU's analysis of the node
graph, assembled straight from the Fourier coefficients with no sparse
matrix formed, in single precision refined to double when asked; on
three tones it is the fastest measured. [`Automatic`](@ref) chooses by
the number of pumps and the memory the factors would take, sized from
the symbolic analysis. A preconditioner that leaves the preconditioned
operator nearly singular in a few directions stalls GMRES, and no per
mode criterion predicts which directions those are, so the rescue is
escalation: a band grows by one offset per pump and any other set to
the full Jacobian, within the memory the grown factors are predicted to
take, and a single precision full set escalates to double. This fires
once or twice on a strongly pumped line and not at all otherwise. The
Floquet deflation instead measures those directions, harvesting the
singular directions of the residual image from the Arnoldi basis as
physical correction vectors carried across the Newton path, which needs
no sparse factorization at all; it is off by default, having measured a
net loss against escalation on long lines.

The preconditioner is rebuilt before every Newton step by default, a
reproducible path; the probe policy measures the one-step reduction of
the stale preconditioner against the fresh one and skips a rebuild when
the predicted extra Arnoldi steps are cheaper. Across a sweep the
structure, the coupling set and the symbolic factorization are reused,
with the values rebound.

The staged method is source continuation on a ladder of retained
harmonic sets. Near a critical drive the Newton basin is small and the
iterations many, so they are spent where they are cheap: the drive is
climbed from a first fraction in warm started steps with a small set of
harmonics retained as unknowns, each larger set warm started from the
last by matching mode tuples, while the nonlinearity is always evaluated
on the full transform grid so that every stage sees the same aliasing.
The schedule adapts in both directions, since each truncation has its
own solvability boundary: a stalled drive step is halved, a stall at the
minimum step grows the retained set at the current drive, and a point
that fails to reconverge after growth retreats the drive on the new set.
Interior points converge only loosely under a small budget, and the one
expensive solve, the finest set at full drive, starts inside the basin.
A stall from a point converged on the finest set is reported as
bracketing a fold, the end of the solution branch, between the last
converged drive fraction and the stalled one, after one more attempt at
full drive with the caller's method; no path throws, and the walk is
recorded in `solverinfo.stages`.

## The linearized sweep

About a converged operating point the circuit is a linear system varying
periodically in time, and a weak signal at `w_s` scatters into the
idlers at `w_s + k w_p`. The signal modes are the signal offset by the
retained pump harmonics, and for each signal frequency the system matrix
is

```math
A(\omega_s) = A_{L_j} + L^{-1} + i\,\mathrm{diag}(\omega_m) G - \mathrm{diag}(\omega_m^2) C,
```

with `w_m` the frequency of each signal mode, the entries of the
negative frequency modes conjugated, and `A_{L_j}` the modulation term,
the incidence triple product of the Fourier coefficients of
`cos(phi(t))` of the pump at the mode differences, which is what couples
the signal to its idlers. The sparsity pattern is that of the nonlinear
Jacobian and is analyzed once; each frequency refactorizes the values.
A unit current source at every port and mode is the right hand side,
and the scattering parameters are the ratio of the outgoing to the
incoming power waves at the ports, in units of photon flux: the wave of
a mode is scaled by `1/sqrt(|w| Z)`, so a conversion between frequencies
is read in photons, and the zero frequency mode has no wave.

The noise and the sensitivities read the transposed system at the same
factorization. By the adjoint identity the response at an output port to
a source anywhere in the circuit is that source contracted against the
transposed solution driven at the port, one solve per port rather than
one per noise channel, so the noise scattering matrix from every
dissipative element to the ports costs the ports' solves; the sign the
transposed route owes between a positive and a negative frequency mode
is restored, so that a resistor's channel is exactly a port of its
impedance in vacuum. The quantum efficiency of an output is the signal
power it carries over the total, the scattered inputs plus the noise
channels weighted by their occupation `2 nbar + 1` at their temperature;
the commutation relation is the signed sum of the scattered powers, `+1`
for a positive frequency output when the scattering matrix is complete;
and the added noise covariance is the occupation weighted product of the
noise scattering matrix with itself, the `Y` of the Gaussian channel
whose `X` is `S`.

The frequencies are split into batches over the threads, each with its
own workspace; on a device a batch of frequencies, which shares one
pattern, is assembled by one kernel and factorized and solved as a
uniform batch, by cuDSS for a sparse factorization and by the batched
block factorization for dense node blocks.

## Scattering blocks in the frequency domain

A scattering block enters as the hybrid equation of each port,
`(I - S(w)) R^{-1/2} v - (I + S(w)) R^{1/2} i = 0` at the block's own
reference impedances, with the port currents as auxiliary unknowns
coupled into the Kirchhoff rows, so nothing of `S` is inverted and an
ideal through, short or open is exact. The data is used at its native
reference impedances; the conversion to the analysis impedances happens
in the wave domain at the ports. At zero frequency a block states its
limit, an open, a short, a through, or its own evaluation. A lossy block
emits a noise wave of covariance `I - S S'`, Bosma's relation, one
channel per port, correlated by that covariance, at its own temperature;
a block declared lossless emits nothing, and the declaration is
validated where the data allows. Passivity of constant and tabulated
data is checked at construction, and a rational block's by its largest
singular value over all frequencies.

## Sensitivities

The derivative of the scattering parameters with respect to a relative
perturbation of a component value has two parts. At a fixed operating
point each component's stamp, its own contribution to the linear term
applied to the solution, is contracted against the forward and the
transposed solutions, one product per component and frequency. The shift
of the operating point follows from the implicit function theorem on the
residual in the real representation, `dx/dr = -J^{-1} dF/dr` with the
exact real Jacobian assembled and factorized once at the converged
point; its effect on the scattering parameters is contracted either
forward, one product per component, or in reverse, the output
functionals pushed through the transposed pump Jacobian once per output
port and mode, which is chosen when the components outnumber the
outputs. A design parameter of a circuit builder carries the exact
direction of every component value it touches, merged into one
contraction, which is exact for a parameter that rotates a complex value.

## Validation

The solvers are checked against each other and against references: the
matrix-free Jacobian-vector product against the assembled real Jacobian
to machine precision, the auxiliary current formulation against the
nodal one as a Schur complement identity, the noise of a resistor
against a port of its impedance, the scattering blocks against their
explicit networks, the sensitivities against finite differences, and the
amplifiers against Keysight ADS and against Fourier analysis of WRspice
in the reference below; the time domain solver's stationary limits are
checked against these solvers in turn.
