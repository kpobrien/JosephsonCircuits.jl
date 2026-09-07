# Theory and implementation of the time domain solver

The time domain solver integrates the same compiled circuit, on the
same node flux unknowns, that the harmonic balance solvers use, so that
a circuit is written once and a result can be compared between the two
domains. It is written for the way the package is used: many pump or
signal conditions of one circuit as one solve, on a CPU or on a GPU,
with the exact tangent and adjoint of the steps taken, and the quantum
noise about the recorded trajectory. This page describes the equations,
the stepping rules, the treatment of the algebraic directions, the
scattering blocks, lines and fitted data, the responses, and what runs
where; the [usage page](transient.md) shows the calls, and the
[noise page](transientnoise.md) the fluctuations.

The code lives in `src/transient/`: `system.jl` compiles a
[`transientproblem`](@ref) from a circuit, classifying its directions;
`solve.jl` builds the scaled matrices and the trapezoidal and backward
Euler steps; `gauss.jl` the Gauss-Legendre rule, its projection of the
endpoint and the stage algebra of the rational blocks; `batch.jl` the
stepping of a batch of conditions, the rational operators, the line
histories, the tangent and the adjoint; `sensitivity.jl` the argument
handling of the responses; `noise.jl` the baths and the contraction;
`quantum.jl` and `iq.jl` the temporal mode and I/Q measurements. The
vector fitter is `src/circuit/vectorfit.jl`, next to the components it
produces.

## The equations

The unknowns are the node fluxes `phi` in units of the reduced flux
quantum, the same node flux basis as harmonic balance, and the equations
are

```math
C\ddot\phi + G\dot\phi + L^{-1}\phi + I_c\sin\phi_b = I(t),
```

with the capacitance, conductance and inverse inductance matrices of
[`numericmatrices`](@ref) at one mode, the junction term on the junction
branches through the incidence matrix, and the modified nodal analysis
augmentation of [`hbnlsolve`](@ref): a mutually coupled inductor keeps its
branch current as an auxiliary unknown with its constitutive equation,
which keeps the system well posed as the coupling approaches unity, and a
subnetwork no element connects to ground has its free flux offset fixed
by a gauge row. The rows are scaled by `Lscale/phi0` with `Lscale` the
mean inductance of the circuit, fixed for the problem so that a state,
auxiliary currents included, means the same at every step and the final
state of one solve starts another at any step. The capacitance matrix
may be singular: the solver treats the system as the differential
algebraic system it is, without inverting `C` or adding artificial
conductances, and checks at the start that the state satisfies the
algebraic equations along every direction without capacitance to ground,
a node no capacitor touches, a capacitive island, a coupled inductor or
gauge row, and along every direction without conductance either, a
subnetwork no capacitor or resistor connects to ground, that the rate
satisfies the differentiated equation.

Three classes of direction follow, and the rules treat them differently.
Where the capacitance is nonsingular the equations are differential and
each rule has its full order. Along a direction without capacitance
that a resistor reaches, the equation is algebraic in the rate, index
one in the flux and rate: the flux keeps the rule's order and the rate
along it converges at second order under either rule. Along a direction
without conductance either, the equation constrains the flux alone and
the rate is its derivative, index two. Which directions without
capacitance are which is read off the equations themselves: along the
capacitor free subnetworks and the scattering blocks' port currents the
equations are linear in the rates and the currents, and the null space
of that rate system gives the flux directions no equation's rate
determines, as columns, and the combinations of the equations, node
rows and block rows, that constrain them, in which the block currents
cancel. So an open block port adds no conductance and leaves a
junction's node constrained, a short or a through joins the nodes it
ties, a resistive port grounds one, and the two sides of a through
between a driven junction and an inductor are one direction, as the
equations say rather than as a graph of the ports would guess. A
constraint no source drives, a coupled inductor's row, a gauge row or
an inductor between two nodes, is an invariant of both rules and holds
exactly, while a
junction on such a direction or a source driving it is not. There the
trapezoidal rule, which solves the constraint at every endpoint, keeps
second order from a consistent start, and the Gauss-Legendre rule, whose
endpoint is the collocation quadratic extrapolated, would not converge in
the rate at all, so it projects each endpoint onto the constraint and
takes the rate along the direction from the derivative of the cubic
through the state, the two stages and the projected endpoint, third
order; the tangent and the adjoint differentiate the projection with the
step. A junction across an unterminated port with no capacitance is
the smallest such circuit, and its flux is then the constraint's own
solution to the Newton tolerance.


## The stepping rules

Three rules step the equations: [`Trapezoidal`](@ref), the trapezoidal
rule on the flux and on its rate, which is Newmark's rule with the
averaging parameters, second order and free of numerical damping;
[`GaussLegendre`](@ref), the two stage Gauss-Legendre collocation,
fourth order, A-stable and symplectic; and [`BackwardEuler`](@ref),
first order and strongly damping, a reference for checking that a
result does not depend on the rule. A trapezoidal step is one implicit
equation in the new flux,

```math
[\alpha C + \beta G + L^{-1}]\,\phi_{n+1} + I_c\sin\phi_{b,n+1} = r_n,
```

solved by Newton from the previous increment, with the Jacobian
assembled by the real Jacobian plan of harmonic balance at one mode,
whose pattern and symbolic analysis are fixed for the whole solve. A
factorization is kept across steps while Newton converges in one
correction with it, and refreshed at the current iterate when a second
correction is needed, so a linear circuit factorizes once and a junction
driven weakly nearly so.

A Gauss-Legendre step solves the two stage equations of the collocation
together, on the stage increments. The stage matrix of the tableau has
the conjugate eigenvalues `3 ± i sqrt(3)`, so with one junction
stiffness frozen at the mean of the two stages the stages decouple into
the complex system `(mu/h)^2 C + (mu/h) G + L + J*` and its conjugate,
of which one is factorized: the real part is assembled by the same real
Jacobian plan and the imaginary part is constant on its pattern. The
Newton iteration on the true stage residual is then a simplified one,
converging linearly at a rate set by how far the two stage stiffnesses
sit from their mean, about two or three corrections per step, and the
factorization is refreshed only when the contraction of the residual
says the frozen operator has drifted, so a pumped amplifier runs its
whole trajectory on one complex factorization. The tangent and the
adjoint differentiate the full stage equations, with the two stage
stiffnesses the solve recorded as `phases`, and the projection of the
endpoint where the circuit has an algebraic direction to project, with
the endpoint phases of its junctions recorded as `endphases`, and solve
them exactly by iteration on the same complex factorization; a current
on the recorded grid is read at the stage times through a cubic Lagrange
stencil, so a smooth current keeps the fourth order, and a current given
at the grid and the stage times of each step is read as it is, which is
how the pulsed gain and the noise drive their probes and baths, so a
pulse keeps its support and a tone its exact phase.

The step of the Gauss-Legendre rule, in the scaled unknowns, solves for
the two stage increments `d_1, d_2` of the flux from the state `x_n`
with rate `v_n`,

```math
\sum_l \frac{(A^{-2})_{il}}{h^2} C d_l + \sum_l \frac{(A^{-1})_{il}}{h} G d_l + L (x_n + d_i) + J(x_n + d_i)
= b(t_n + c_i h) + \frac{(A^{-1} \mathbf 1)_i}{h} C v_n,
```

with `A` the tableau, `c` its nodes and `J` the junction current, and
takes the endpoint as `x_{n+1} = x_n + e_x \cdot d` and
`v_{n+1} = v_n + e_v \cdot d / h` with the weights of the collocation
polynomial. The next step's stages start on the collocation polynomial
extrapolated. The frozen operator of the Newton iteration is the complex
matrix `(mu/h)^2 C + (mu/h) G + L + J*` with `mu` the eigenvalue of the
tableau's inverse and `J*` the mean stiffness, factorized once and
reused while the first correction of a step contracts the residual by
more than a quarter. The residual is evaluated per condition and per
direction against its own right hand side down to a roundoff floor set
by the magnitudes of the terms it sums, so a weak direction is not left
at the tolerance of a strong one.

The Newton engine of a batch accepts per condition: each column has its
own tolerance from the size of its right hand side, a rejected
correction is retried from its base point after a fresh factorization,
and a step that does not converge on every condition throws rather than
proceed. The factorization of a batch is one per condition, KLU on the
CPU and the uniform batch of cuDSS on a device, on one pattern with one
symbolic analysis; a circuit with scattering blocks has an unsymmetric
operator, and its adjoint uses the transposed factorization of the same
matrix.

### The projection of the endpoint

Along an algebraic direction the Gauss-Legendre endpoint, the
collocation quadratic extrapolated, satisfies the constraint only to the
order of the extrapolation, and the rate along it does not converge at
all, so every endpoint is projected. The directions `Z` and the
constraints `Z'`, the right and left null vectors of the rate system,
are kept apart, since a block makes the operator unsymmetric: Newton on
the coefficients `alpha` of `Z` solves `Z' (L (x + Z alpha) + J(x + Z alpha) - b) = 0`
per condition with the small Jacobians `Z' (L + J'(x)) Z` on the host,
from two products on the backend, to the step's tolerance; a linear
constraint is met by one correction, and a direction no junction,
drive, line or block touches is left alone, since the stages keep a
linear constraint with a constant right hand side exactly. The rate
along the direction is then the derivative at the endpoint of the cubic
through the state, the two stages and the projected endpoint, all of
which satisfy the constraint, extracted along each direction and put
back through `Z`. Last the index one unknowns, the rates along the
capacitor free islands and the block port currents, are read from their
equations at the endpoint through the pseudoinverse of the rate system
on its range, which determines what the equations determine and leaves
the algebraic directions to the projection. The tangent and the adjoint
differentiate the projection with the step, the adjoint through the
transposes of the same small solves with the constraints' rows and the
directions' columns in their places.

## Scattering blocks

A constant real block, an attenuator, a circulator, an ideal through,
short or open, is realized as the linearized solver stamps it: its port
currents are auxiliary unknowns after the coupled inductor currents,
each port carries the hybrid constitutive equation
`(I - S) R^(-1/2) v - (I + S) R^(1/2) i = 0` at the block's own reference
impedances, and the currents enter the node equations. Nothing is
inverted, so a short or an open is exact, and a block between ports of
different reference impedances needs no renormalization, since the
mismatch is the shared node. The rows are unsymmetric, which the adjoint
follows with transposed solves, on a device by a second factorization of
the transposed operator built when an adjoint first asks for it. An
attenuator equals its resistive network to roundoff, and a pumped
amplifier behind a 3 dB pad has the gain the linearized solver gives it.

```julia
g = 10^(-3/20)
pad = ScatteringParameters([0.0 g; g 0.0]; zref = 50.0)
circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:pad, 1, 2, pad),
    (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)), (:cj, 3, 0, Capacitor(1000e-15))])
problem = transientproblem(circuit; sources = [TransientSource(1, t -> 2Ip*rise(t)*cospi(2fp*t))])
solution = transientsolve(problem, (0.0, 200e-9); dt = 2.6e-12, method = GaussLegendre())
```

A [`RationalScattering`](@ref) block is a passive multiport given as the
real state space realization `S(s) = D + C (s I - A)^(-1) B`, the form
a vector fit of measured or simulated scattering data takes, validated
as stable and passive at construction by its largest singular value
over every frequency, found by the level set iteration: a lower bound
from samples is raised to a level, the frequencies where a singular
value equals the level are the imaginary eigenvalues of a pencil formed
without inverting `I - D'D`, and the largest singular value between
consecutive ones raises the bound until none reaches the level, which
finds a peak however narrow, where a sample or a test at the crossings
of one would miss it. Nothing is inferred about a rational block's
loss: it keeps the channels of its loss, however small, in both
solvers, and only a declaration of `Lossless()` is validated, by the
same iteration on the block and on its inverse, which can only refuse
it; no test of the coefficients of `I - S(-s)'S(s)` could stand in,
since a narrow notch has coefficients of the square of its width and a
loss of one at its center. The harmonic balance solvers
evaluate it at every frequency; in time its states step with the
circuit: the two stage values of the states are linear in the state at
the start and the incident waves at the stages, so the reflected waves
at the stages are an exact algebra on the stage unknowns, constant
matrices factorized once per step size, and in the complex stage basis
that algebra is the block's own scattering matrix at the stage
frequency `mu/h`, which the frozen operator carries exactly. A series
inductor written as a one state block equals the explicit inductor to
roundoff at every step size. The block's states are the fourth member
of a [`transientstate`](@ref), at rest under the port voltages by
default, and `record = :states` keeps them as `blockstates`. The
tangent and the adjoint carry the perturbation of the states and its
cotangent through the same algebra and its transpose, a perturbation of
the states being the fourth member of the tangent's initial state and
the adjoint returning theirs as `initialstates`; on the inductor block
both equal the explicit inductor's to roundoff.


The stage algebra of a rational block is exact. With the tableau `A`
and the block's realization, the two stage values of the states are
`Z = Z_z z + Z_u U` with `Z_z = (I - h A \otimes A_b)^{-1} (\mathbf 1 \otimes I)`
and `Z_u = (I - h A \otimes A_b)^{-1} (h A \otimes B_b)`, the reflected
waves at the stages `C_b Z`, and the state after the step
`z' = E_z z + E_u U` with `E_z = I + (h/2)(A_b, A_b) Z_z` and
`E_u = (h/2)((A_b, A_b) Z_u + (B_b, B_b))`, the Gauss weights being
halves. The incident waves `U` at the stages are linear in the stage
rates and currents, so the whole coupling of every block of a batch is
grouped into sparse operators on the stage stacked unknowns
`[stage 1; stage 2]`: `M_d = S C Z_u W_d G` on the increments,
`M_x = S C Z_u W_x G` on the stage values, `M_s = S C Z_z` on the
states, and `E_d, E_x, E_z` for the states' update, with `G` the gather
of the rates across the ports and the port currents, `W` the incident
waves and `S` the scatter onto the blocks' rows; their transposes serve
the adjoint, and `(S_1 C)'` carries the cotangent of the resting waves
the endpoint reading sees to the states. All of them are built once per
step size and live on the backend, so the reflected waves of every
block at both stages are a few products with no per block or per
condition work. In the complex stage basis the algebra is the block's
own scattering matrix at the stage frequency `mu/h`, which the frozen
operator carries exactly through its entries on the Jacobian's pattern.

## Transmission lines

An ideal line is the method of characteristics: at each port the current
into the line is `v/Z - 2q/sqrt(Z)`, a conductance at the line's own
impedance plus a current from the wave that entered the far port a delay
earlier, and the wave leaving each port, `v/sqrt(Z) - q`, is kept at
every step as the history the far port reads a delay later; a mismatch
to the circuit is the shared node. Nothing resonates in the equations,
since the round trips that make a line's admittance singular at its half
wave resonances in frequency are the recursion through the history, so a
line between mismatched loads has the scattering parameters of the
linearized solver at those resonances too. The delay must be at least a
step, so that every wave a step reads is accepted history, and a read
interpolates the endpoint history with the centered Lagrange stencil of
as many samples as the history holds on both sides of the query, up to
three on each side: the quintic where the delay leaves three accepted
samples past the query, the cubic with two, the linear with one. Every
centered stencil is a contraction, its magnitude at most one at every
frequency, so a wave gains nothing on a round trip through any delay
and a passive line stays passive however short; a one sided stencil is
not, and a mismatched line of 1.2 steps' delay read through one grew
without bound. The quintic keeps the rule's fourth order, so a delay
of at least four steps keeps the accuracy, and a shorter delay reads
at lower order. A centered read spreads a pulse's edge over the stencil,
so a wave shows a small precursor before its arrival, of the read's
error. The history before the start is the initial waves, constant, so
the solve and its responses share one linear map from the start: the
tangent along a drive is the linear circuit's own solve under it. The
history lives in a ring as long as the longest delay and the stencil's
reach, so its memory does not grow with the record; the record of the
waves leaving every port is the `linewaves` of the solution, and with
checkpoints the history before each checkpoint is kept in their `waves`
instead, from which a replay starts. A line's state at the start is the
wave leaving each port, from the port voltages and the direct current
the line carries, which [`transientstate`](@ref) takes as `linecurrents`.
The tangent and the adjoint carry their own rings of the perturbation,
read at the stages and the endpoints through the same stencils and
scattered back through their transpose; a perturbation of the
prehistory is the third member of the tangent's initial state, and the
adjoint returns the cotangent of the prehistory as `initialwaves`.


## Fitting measured data

Measured or simulated data become such a block by fitting. Given a
tabulated or Touchstone block, or any block sampled at given
frequencies, `RationalScattering(block, npoles)` fits a real state
space realization at a set of common poles by vector fitting, makes it
passive where the fit strays above unit singular value, and keeps the
data's reference impedances, grounding and noise model, so a fit of a
block at a temperature emits the noise of its loss at that temperature
in time as in the frequency domain. The fit eliminates each entry's own
unknowns from its own samples before solving for the shared poles, so
its storage grows with the entries and the samples rather than their
product, though every entry still visits every sample, and its
realization takes one state per rank of each residue. The poles the
data does not need
drift out of the band or settle on one another and are dropped, so ask
for as many as the data might need; too few settle on a poor fit. A
delay is not a rational function, so a cable is a
[`TransmissionLine`](@ref) of its delay in cascade with a fit of the
data with that delay removed.

```julia
data = ScatteringParameters((2pi .* frequencies, S); nports = 2, zref = 50.0,
    noise = ThermalEquilibrium(0.05))
fitted = RationalScattering(data, 8)
```


The fit is the relaxed vector fitting of Gustavsen and Semlyen: a set of
common poles, started as complex pairs spread over the band, is
relocated by the zeros of a weight function whose residues and constant
are fitted together with every entry's residues and constant, under the
relaxed normalization that the real part of the weight averages to one
over the samples; without the weight's constant as an unknown the exact
poles are not a fixed point, and only fits with a spare pole converge.
An entry's own unknowns enter only its own rows, so each entry's block
is reduced by its QR to the trailing triangle of the weight's columns
and the weight is solved from all the reductions at once, as
Deschrijver's fast fit does, so nothing the size of every entry's every
sample is formed. The iteration stops when the poles settle to a part in
a billion or the fit is at roundoff, where the weight's columns are
combinations of the entries' and the relocation would be arbitrary. The
poles the data does not need drift out of the band or settle in
clusters split by less than the data resolves; each candidate drop or
merge is kept only if the residues refitted reproduce the data as
closely as before, within a factor of two above the roundoff floor. The
realization takes one state per rank of each residue. Passivity is
enforced where the largest singular value exceeds one by the least norm
perturbation of the residues and the constant at the worst points of
the violation bands, and the result is validated as any rational block
is: by its largest singular value over every frequency, found by the
level set iteration of Boyd, Balakrishnan, Bruinsma and Steinbuch on a
pencil formed without inverting `I - D'D`, which finds a peak however
narrow. Nothing is inferred about a block's losslessness; a declared
`Lossless()` is validated by the same iteration on the block and on its
inverse.

## The tangent and the adjoint

The tangent linearizes the recorded steps about the recorded junction
phases: on each step it solves the linearized stage equations with the
same frozen operator, iterating to the exact solution, reads a current
on the recorded grid at the stage times through a cubic Lagrange
stencil or takes it at the stages when it is given there, carries the
perturbation of the block states and of the line histories, and
differentiates the projection of the endpoint. The adjoint is the exact
transpose of those steps in reverse: the multipliers of the stage
equations on the transposed factorization, the cotangent of the
projection through the transposed small solves, the cotangents of the
states through the transposed grouped operators, and of the line waves
through the scatter of the stencils onto the samples each read touched,
handed back as currents on the grid and at the stages, and the
sensitivities to the initial flux, rate, waves and states. A tangent or
adjoint of a batch runs every condition's directions or objectives in
one pass on one factorization per condition. Under a record of
checkpoints the phases of a window of steps are replayed from the
checkpoint's state, stage predictor, block states and line history
before it, and the window is walked forward by the tangent or backward
by the adjoint, so the memory of a record of any length is bounded by
the checkpoints and one window; a replay that does not reach the next
checkpoint is an error, not a warning.

The histories of the lines are rings as long as the longest delay and
the stencil's reach in the solve, the tangent and the adjoint alike: the
adjoint clears a column's slot once its step has consumed it, since a
read at a step touches only columns before that step's endpoint, and
what the ring holds at the end is the cotangent of the prehistory. An
adjoint given a sink stores no currents at all: each column is handed
over once no remaining step touches it.

## Quantum noise

The [noise page](transientnoise.md) gives the physics. In the
implementation every bath is a target of the adjoint, and the noise of
a window of measured modes is one adjoint with the measured quadratures
as objectives: each column of the bath kernels is contracted against
the cosine and sine of every bath frequency at the time the adjoint
hands it over, through a block buffer and one product on the backend,
into an accumulator of `16 nb m nf N` bytes that is tiled over
conditions and frequencies to a quarter of the backend's free memory.
The stationary initial term of every bath enters through the adjoint's
initial flux, rate, waves and states with one stationary solve per
frequency and objective, on the joint system of the flux phasors, the
waves leaving the line ports and the block states, so no bath's initial
state is ever built. A rational block's ports are correlated channels
contracted directly with their covariance `I - S S'` on the backend,
the loss matrices evaluated once per frequency tile; the demodulation's
phasor convention is the conjugate of the linearized solver's, so the
covariance reads `conj(K)` there, which only a loss matrix with
imaginary off-diagonal entries can tell.

## Execution on a device

On a device the state, the products, the junction term, the Jacobian
assembly, the factorization and the solves of every step run there
through the package's device sparse matrices and kernels, and so do the
rational operators, the line histories with their gather and scatter
kernels, the records, the checkpoints and the noise contraction; the
compiled circuit, the sample times, the source callables and the
stencils of the line reads stay on the host, and one vector of drive
currents and one small table of stencils are transferred each step. The
endpoint projection is host assisted: it downloads the projected
junctions' phases, the line forcing and the resting waves of the blocks,
solves the small systems per condition on the host and uploads the
corrections, which is a handful of small transfers per step. The Newton
control needs scalar norms, so a step synchronizes a few times, and a
small circuit is not faster on a GPU; the throughput lies in the batch,
whose launches serve every condition.

Blocks and lines cost a batch little beyond the circuit itself. The
stage algebra of every rational block is grouped into a few sparse
operators on the stage stacked unknowns of the whole batch, built once
per step size and kept on the backend, so the reflected waves of all
the blocks at both stages and their states after the step are a few
products with no per-block or per-condition work; the transposes serve
the adjoint the same way. The history of the line waves lives on the
backend too, and a read is one gather over the ports and the conditions
through the stencils set on the host, the same for every condition, its
transpose the scatter the adjoint uses. On a small amplifier behind a
six state fitted block or a line, a condition in a batch of 512 on the
device costs within a half of one alone.

```julia
circuit = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:cable, 1, 2, TransmissionLine(60.0, 0.09)),
    (:cc, 2, 3, Capacitor(100e-15)), (:jj, 3, 0, JosephsonJunction(1000e-12)), (:cj, 3, 0, Capacitor(1000e-15))])
```

## Validation

The test suites compare the solver against references rather than
against itself: the trapezoidal and Gauss-Legendre rules against
analytic RC, LC and junction solutions and against each other at their
orders; the stationary limit of a pumped amplifier, its gain and quantum
efficiency, against [`hbsolve`](@ref) to the order of the rule; the
scattering of lines and blocks against [`hblinsolve`](@ref) at their
resonances; the noise of every bath type against the linearized solver's
covariance, cold and warm; the tangent against a driven linear solve and
the adjoint against the tangent by the dot product; a fitted block
against the explicit circuit it was fitted from; and the pump response of
a JPA against WRspice at the signal frequency. The GPU suite compares
every device path against the host to a part in 1e8.
