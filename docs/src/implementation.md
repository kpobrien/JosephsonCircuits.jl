# Implementation notes

These notes describe numerical storage and execution. For constructing a
circuit or choosing an analysis, start with the [user guides](index.md#How-the-documentation-is-organized).

## Source map

| Directory | Responsibility |
|---|---|
| `src/circuit/` | Component models, parsing, elaboration, binding, matrices, MNA, and vector fitting |
| `src/harmonics/` | Mode layout, transforms, residuals, exact derivatives, and DC augmentation |
| `src/solvers/` | Nonlinear methods, continuation, preconditioners, and factorization |
| `src/linearized/` | Signal sweeps, temporal poles, output normalization, noise, and sensitivities |
| `src/transient/` | Time stepping, constraints, responses, temporal measurements, and baths |
| `src/networks/` | Network parameters and their conversions, closed-form networks, the connection of scattering networks with their noise, and quantum optics |
| `src/spice/` | SPICE netlist export, the WRspice runs of the comparisons and of the `WRspice` transient method, and SPICE output readers |
| `ext/` | Extensions loaded with CUDA.jl, CUDSS.jl, Krylov.jl, SciMLBase.jl, Symbolics.jl and XicTools_jll |

The corresponding tests generally mirror these directories. Cross-domain
and external-reference comparisons live in the top-level test files.

## Circuit compilation

Parsing resolves instance and terminal names to integer indices.
Elaboration flattens subcircuits. Compilation builds component tables and
topology, and value binding gathers concrete numerical arrays grouped by
component kind. This keeps flexible input syntax out of numerical kernels.

A circuit description can contain heterogeneous component types. Vectors
avoid making every circuit size a different tuple type. Wire numbers,
net numbers, and instance offsets are runtime data, not type parameters.
Flat arrays with offsets hold terminal connectivity.

Within one elaboration, repeated subcircuits share parsed definitions and
interface-key indexes. Independent elaborations rebuild those indexes,
so they observe edits to retained input collections. Custom keys follow
Julia's `isequal`/`hash` contract.

## Harmonic operators and factorizations

The residual and directional derivatives share precomputed gather and
transform plans. Device kernels assign work by output slot, reading index
maps without scatter or atomics. The junction derivative is cached for
multiple products at one Newton point.

Exact real-Jacobian assembly reads mode-difference coefficients of the
junction derivative. Node-block factorization groups the circuit graph
into supernodes ordered from KLU's symbolic analysis, then assembles dense
mode blocks directly. Its storage depends on graph fill as well as the
number of modes; mode count alone is not a memory estimate.

Linearized frequency batches own their workspaces. Pump sensitivities use
the exact real Jacobian at the converged point. Forward contractions are
organized by parameter; reverse contractions by output functional. The
implementation chooses between them according to their counts.

## Transient steps and replay

`system.jl` constructs the scaled problem and classifies constrained
directions. `constraints.jl` projects endpoints and reads algebraic rates.
`solve.jl` contains the Newton engine and trapezoidal/backward-Euler
integration. `gauss.jl` defines the collocation coefficients and block
stage algebra; `batch.jl` steps conditions and implements their responses.

The Gauss–Legendre frozen operator uses one complex factorization per
condition. Acceptance and refresh are per condition, not governed by the
largest residual of the whole batch. An unsymmetric block operator needs
transposed solves for the adjoint.

Checkpoints save state, stage predictor, rational states, and line
prehistory. The original and replayed trajectories restart factorization
at the same checkpoint boundaries. A replay is checked against the next
saved state to roundoff. Source callables must be deterministic for this
contract to hold.

## Grouped block and line operations

Rational-block stage elimination is assembled into sparse operators on
stage-stacked circuit unknowns. They map stage increments, state values,
and prior block states into reflected waves and updated block states.
Their transposes serve the adjoint. The operators are built once per step
size and apply across every condition in a batch.

A line-history ring spans the longest delay and interpolation stencil.
The tangent reads perturbations through the same stencil; the adjoint
scatters their cotangents into the ring. Once no future reverse step can
touch a column, it is cleared or handed to a streaming sink. The remaining
ring at the end represents sensitivity to the prehistory.

An adjoint `sink` consumes current columns as they become final, in reverse
time order. This is useful for contractions that do not need to retain a
full current-gradient history, including the noise calculation.

## Pumped-block noise

An ordinary passive block has commutator `K=I-S*S'` and thermal covariance
`(nbar + 1/2)K`. An active block's `NoiseCovariance(V)` is checked using
both `V-K/2` and `V+K/2`. Its channels include contributions with both
commutator signs.

A pumped block instead couples frequencies on ladders separated by its
pump. Conjugate partners introduce covariance terms relating frequencies
whose sum is a pump harmonic. A covariance supplied at a conjugate
frequency must satisfy the transpose symmetry
`V_k(-nu-k*wp)=transpose(V_k(nu))`.

A fitted pumped block's transfer functions have their own commutator.
Noise completion is performed on padded mode ladders against that
commutator. Padding keeps the effective noise model consistent when a
calculation requests fewer modes. Every requested output is completed;
a mode beyond the stored sidebands still carries the noise required by
its effective commutator. `noisetol` limits the completion cost separately
from the scattering fit's `tol`.

Frequency grouping avoids building one dense matrix across unrelated
ladders. At fixed ladder width, cost scales with the number of ladders,
plus sorting/grouping. Within a ladder, covariance and partner storage
grow quadratically with width, and dense completion grows cubically.
If every bath frequency belongs to one ladder, that wide ladder remains
a large coupled problem.

## Noise contraction

`quantum.jl` defines temporal-mode measurements; `noise.jl` constructs baths
and contracts their responses. A block buffer of adjoint current columns
is multiplied by sinusoidal bath quadratures on the backend. Conditions
and uncorrelated frequency sets can be tiled to control the accumulator
memory described in [performance](performance.md#Noise-calculation-cost).

The stationary initial term solves for node-flux phasors, outgoing line
waves, and rational-block states together. It contracts with adjoint
initial-state derivatives, avoiding construction of every bath's initial
state separately. Conditions sharing an initial state can share its
stationary factorization.

The temporal cosine/sine convention conjugates the corresponding complex
HB covariance representation. Complex off-diagonal loss entries therefore
matter; tests using only real diagonal losses cannot establish that this
conversion is correct.

## Device execution

Device-resident work includes state products, junction evaluations,
Jacobian assembly, factorization, rational operators, line histories,
checkpoint arrays, and noise contractions. Host work includes source
callables, time grids, line-stencil preparation, and nonlinear control.

Endpoint projection is host-assisted: selected phase, forcing, and block
information is transferred for the small constraint solves, then
corrections return to the device. The Newton loop also reads per-condition
norms and masks. These synchronization costs explain why batching can
matter more than accelerating an isolated small circuit.

On the CPU, batches use KLU factors per condition. Device batches use
uniform cuDSS solves over a shared pattern. An adjoint of an unsymmetric
operator can require a separate transposed factorization on the device.

## Documentation workflow

From the repository root, prepare the documentation environment and build:

```sh
julia --project=docs -e 'using Pkg; Pkg.develop(path="."); Pkg.instantiate()'
julia --project=docs docs/make.jl
```

Named `@example` blocks execute during the build and can share state within
a page. Keep small reference assertions with the example they validate.
Plain `julia` fences are display-only; use them for optional dependencies,
large device runs, or explicitly identified continuations of a setup.
The amplifier recipes share their displayed circuit builders with small
executable checks; their full-size runs remain optional. These checks
cover solver convergence and finite responses, with commutator checks for
the amplifier models, and compare the pump-off circuit with a closed form
where the recipe has no other independent reference. They do not
reproduce the published full-device gain curves.
The package's doctest-only test job does not replace the full docs build.

For a Documenter HTML check without the VitePress/Node rendering step, run
`julia --project=docs docs/check.jl`. Both entry points use the same page
list and execute the same tutorial blocks. The check entry point does not
deploy. Run it when editing links, equations, or executable examples.

Public docstrings are grouped in `docs/src/api/`; `docs/api.jl` defines
their groups once for both the public pages and the internal appendix.
Add a new exported name to the appropriate group. The network library and
the SPICE helpers are grouped by the file that documents them, every
public name of the file, so a name declared `public` in
`src/JosephsonCircuits.jl` appears there and the files' helpers stay in
the appendix. The documentation build stops on a name of a public page
that is neither exported nor public, and on a documented public name that
no page shows.

`docs/live.jl` runs the VitePress development server and rebuilds when
source pages change. A package docstring edit needs a Julia restart or
Revise. Keep optional plotting and external-solver dependencies out of the
small documentation examples unless their environment explicitly provides
them.
