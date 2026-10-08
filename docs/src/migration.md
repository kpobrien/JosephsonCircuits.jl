# Migrating from v0.5.4

This page lists what changes for code written for v0.5.4, the last
registered release, and what to write instead. A deprecated form still
runs and prints a warning naming its replacement; a removed one fails.

## At a glance

| v0.5.4 | This release | What to change |
|---|---|---|
| `@variables`, `@syms`, `@register_symbolic`, `Num` and `Symbolics` exported | Symbolics is an optional extension | Load Symbolics, or name values by symbols; see [symbolic values](#Symbolic-values) |
| Netlists of `(name, node1, node2, value)` tuples, a resistor across each port | [`Circuit`](@ref) of typed components; a [`Port`](@ref) owns its termination | Deprecated; see [netlists](#Tuple-netlists-and-port-terminations) |
| `parsecircuit`, `parsesortcircuit` | Removed | [`compile`](@ref) |
| `calccircuitgraph` | Removed | `compile(circuit).topology` |
| `calcqe`, `calccm`, `calcqe_S_Cnoise`, `calcCnoise(S, Snoise)` | Removed | The `QE`, `CM` and `Cnoise` of a solution (`returnQE`, `returnCM`, `returnCnoise`) |
| Nodes ordered as numbers | A typed circuit orders node names as strings | `compile(circuit; sorting = :number)`; see [node order](#Node-order) |
| `symfreqvar` | [`FrequencyDependent`](@ref) | Deprecated |
| `ftol`, `switchofflinesearchtol`, `alphamin`, `maxharmonics`, `maxpumpharmonics` | `atol`; the others are ignored | Deprecated; see [keywords](#Solver-keywords) |
| `returnZ` and its three relatives, the `Z` outputs | Removed | Convert `S` with `StoZ` |
| Newton's method with a KLU factorization | [`NewtonKrylov`](@ref) | `method = Newton()` for the direct method |
| Noise covariances with the vacuum at one | Symmetrized, the vacuum at one half | Halve covariances written for v0.5.4; see [noise](#Noise-in-quanta) |
| `wrspice_input_transient` and `wrspice_input_ac` in Hz | rad/s | Multiply by `2pi`; see [SPICE](#SPICE-inputs) |
| `phi0` and `Phi0` of CODATA 2014 | From the exact SI `h` and `e` | Results move by about 1e-8 of themselves; see [constants](#Physical-constants) |
| `autonne_takagi` of a real matrix returns `(Λ, M)`, `Λ` increasing | `(Λ, W)`, `Λ` decreasing | See [decompositions](#Quantum-optics-decompositions) |

## Symbolic values

v0.5.4 exported `@variables`, `@syms`, `@register_symbolic`, `Num` and
the module `Symbolics`, and its examples began with
`@variables R Cc Lj Cj`. The package no longer depends on Symbolics, so
after `using JosephsonCircuits` alone `@variables` is undefined. Load
Symbolics beside the package, which supports it as an extension:

```julia
using JosephsonCircuits, Symbolics
@variables Lj Cc Cj
```

or write the values without it. A symbol names a parameter, given a
value by the definitions passed to the analysis, and
`JosephsonCircuits.@params` declares names to combine in expressions:

```@example migration
using JosephsonCircuits
circuit = Circuit([
    (:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:cc, 1, 2, Capacitor(:Cc)),
    (:jj, 2, 0, JosephsonJunction(:Lj)),
    (:cj, 2, 0, Capacitor(:Cj)),
])
defs = Dict(:Lj => 1000e-12, :Cc => 100e-15, :Cj => 1000e-15)
ws = 2pi .* [4.6e9, 4.7e9, 4.8e9]
sol = hblinsolve(ws, circuit, defs)
nothing # hide
```

[`symbolicmatrices`](@ref) of a circuit whose values are symbols carries
the symbols in its matrices, with or without Symbolics.

## Tuple netlists and port terminations

A v0.5.4 netlist, such as the JPA of its README,

```julia
@variables R Cc Lj Cj
circuit = [
    ("P1","1","0",1),
    ("R1","1","0",R),
    ("C1","1","2",Cc),
    ("Lj1","2","0",Lj),
    ("C2","2","0",Cj)]
circuitdefs = Dict(Lj => 1000.0e-12, Cc => 100.0e-15, Cj => 1000.0e-15, R => 50.0)
```

still runs, with a deprecation warning, in `hbsolve`, `hbnlsolve`,
`hblinsolve`, `numericmatrices`, `symbolicmatrices` and `exportnetlist`,
and `Circuit(circuit, circuitdefs)` converts it. Its component types
follow the name prefixes, and a port's reference impedance is the
resistor across it, which the port adopts as its termination. The same
circuit in typed components is the circuit of the
[previous section](#Symbolic-values), or with numbers:

```@example migration
jpa = Circuit([
    (:p1, 1, 0, Port(1; Z0 = 50.0)),
    (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1000e-12)),
    (:cj, 2, 0, Capacitor(1000e-15)),
])
nothing # hide
```

`Port(1; Z0 = 50.0)` owns a matched 50 Ω termination, so the resistor
goes: kept, it would be a second load across the port, with noise of its
own. Keep the resistors that model dissipation in the device. A tuple
netlist file, read and written by `import_netlist` and `export_netlist`,
is deprecated with the format.

A node with no inductive path to ground needs no large inductor to fix
its static flux, as v0.5.4's flux-pumped example had: a gauge fixes it,
and a resistor carries direct current; see
[direct current](recipes/dc.md).

## Compiling a circuit

`parsecircuit` and `parsesortcircuit` are removed. [`compile`](@ref)
returns a [`CompiledCircuit`](@ref JosephsonCircuits.CompiledCircuit) with the tables
`parsesortcircuit` returned, `nodenames`, `nodeindices`,
`componentnames`, `componenttypes`, `componentvalues` and
`componentnamedict`. The solvers, [`numericmatrices`](@ref) and
[`symbolicmatrices`](@ref) take the circuit or its compiled form. A
compiled circuit carries its topology, the branches and their oriented
incidence matrix, as `compile(circuit).topology`, so no solver takes a
separate circuit graph, and `calccircuitgraph`, with the spanning tree
and the loops it added, is removed.

## Node order

A tuple netlist ordered its nodes as numbers, and still does. `compile`
orders a typed circuit's node names as strings by default,
`sorting = :name`, so integer node names past 9 come in the order `1`,
`10`, `11`, `2`, ..., and code that reads `nodeflux` or `voltage` by
position reads other nodes. Read them by node name through the keyed
arrays, or compile with the numeric order and pass the compiled circuit
to the solver:

```@example migration
chain = Circuit(vcat([(:p1, 1, 0, Port(1))],
    [(Symbol(:c, i), i, i + 1, Capacitor(1e-12)) for i in 1:10],
    [(:p2, 11, 0, Port(2))]))
compiled = compile(chain; sorting = :number)
@assert compiled.nodenames == string.(0:11)
sol = hblinsolve(ws, compiled)
nothing # hide
```

The solvers take `sorting` only for a tuple netlist.

## Frequency-dependent values

The `symfreqvar` keyword, which named a symbolic variable as the
frequency in the component values, is deprecated. Write the value as a
[`FrequencyDependent`](@ref) closure of the angular frequency in rad/s:

```@example migration
wc = 2pi*10e9
resistor = Resistor(FrequencyDependent(w -> 50.0*(1 + im*w/wc)))
nothing # hide
```

A closure is evaluated at the magnitude of a mode's frequency, and a
negative frequency takes the complex conjugate, so
`FrequencyDependent(identity)` is the magnitude of the angular frequency.
A frequency-dependent value combines with parameter expressions.

## Solvers

The pump is solved by [`NewtonKrylov`](@ref) by default: Newton's method
whose steps GMRES solves with a preconditioner that `Automatic()` chooses.
v0.5.4 solved it by Newton's method with a KLU factorization of the
Jacobian, which `method = Newton()` selects. In v0.5.4 the `factorization`
keyword of `hbsolve` set the factorization of the pump solve as well as
the linearized sweep's; the nonlinear solve's is now an option of its
method, `Newton(factorization = KLUfactorization())`. The `factorization`
keyword of `hbnlsolve` is deprecated: beside the default method it warns
and solves with `Newton(factorization = ...)`, as v0.5.4 did, and beside
another method it is refused. That of `hbsolve` and `hblinsolve` sets the
sweep's alone. They refuse `QRfactorization()` with an
`ArgumentError`: the sweep solves the transposed system on the factors at
each signal frequency for the noise, quantum efficiency, commutation
relations, sensitivities and adjoint outputs, which the sparse QR
factorization does not provide. Leave the default, or use
`KLUfactorization()` or `LUfactorization()`. Check
`sol.nonlinear.solverinfo.converged` before using a result: a solve that
does not converge returns its last iterate and warns.

v0.5.4 sampled the nonlinearity on the grid of the retained harmonics.
The grid is set separately by `Nevaluationharmonics`, twice the retained
harmonics by default, so that the leading cubic products do not alias
into the retained modes; results differ from v0.5.4's by that aliasing.
`Nevaluationharmonics = Npumpharmonics` samples on v0.5.4's grid.

### Solver keywords

| v0.5.4 keyword | Now |
|---|---|
| `ftol` | `atol`, the absolute residual tolerance; `ftol` warns and is read as `atol` |
| `switchofflinesearchtol`, `alphamin` | Ignored with a warning; the line search is an option of the method, `Backtracking` |
| `maxharmonics` of `hbnlsolve`, `maxpumpharmonics` of `hbsolve` | Ignored with a warning; `Nharmonics` and `Npumpharmonics` are the retained harmonics and `Nevaluationharmonics` the grid |
| `symfreqvar` | Deprecated; [`FrequencyDependent`](@ref) |
| `sorting` | A keyword of [`compile`](@ref); the solvers take it only for a tuple netlist |
| `returnZ`, `returnZadjoint`, `returnZsensitivity`, `returnZsensitivityadjoint` | Removed with a warning; convert `S` with `JosephsonCircuits.StoZ` |
| `factorization` of `hbnlsolve` | Deprecated: an option of the method; beside the default method it warns and solves with `Newton(factorization = ...)` |
| `factorization` of `hbsolve` and `hblinsolve` | The linearized sweep's alone, by default `BlockFactorization()` at two or more tones when its factors fit and the backend's sparse factorization otherwise; `QRfactorization()` is refused |

### Outputs

[`LinearizedHB`](@ref JosephsonCircuits.LinearizedHB) no longer has the
fields `Z`, `Zadjoint`, `Zsensitivity` and `Zsensitivityadjoint`, and
`portimpedanceindices` is `portimpedances`, the port impedances
themselves. It adds `Cnoise`, `Vout` and `nbar`; `nbar` is returned by
default, like `QE` and `CM`, so a sweep which skipped the noise
calculation with `returnQE = false, returnCM = false` adds
`returnnbar = false`.
[`NonlinearHB`](@ref JosephsonCircuits.NonlinearHB) adds `solverinfo`,
`dcnodevoltage` and `sources`.
The matrices [`numericmatrices`](@ref) and [`symbolicmatrices`](@ref)
return no longer hold `Lbm`, the branch inductances repeated per mode, or
`noiseportimpedanceindices`; the noise channels of a circuit at its
values are `JosephsonCircuits.noiseindices(compile(circuit), matrices.vvn)`.

## Noise in quanta

Noise covariances are symmetrized and counted in quanta, the vacuum being
half a photon, in harmonic balance, the transient and the network
cascade. For the same physics, `Cnoise`, the data of a
[`NoiseCovariance`](@ref), `calcCnoise`, which was `I - S*S'` in v0.5.4,
and the noise covariances `connectS` and `solveS` return are half their
previous values. Halve a covariance written for v0.5.4 before passing it
to `NoiseCovariance`, or to `connectS` and `solveS` as the covariance of a
network, `(name, S, C)`: unhalved, it would count that network's noise
twice against the `(I - S*S')/2` the cascade computes for a network given
without one. `S`, `Snoise`, `CM` and, with the ports at zero temperature,
their default, `QE` are unchanged. [`thermaloccupation`](@ref) returns
the occupation `nbar`.

The quadrature functions of the quantum optics module take `hbar`, one by
default, with the vacuum's covariance `(hbar/2)I`; `hbar = 2` gives their
v0.5.4 values. The ladder functions count the vacuum as `I/2`.

## Physical constants

`phi0` and `Phi0` follow from the exact SI values of the Planck constant
and the elementary charge, `elementary_charge`. `phi0` is larger than its
v0.5.4 value by 7.5e-9 of itself, and `Phi0` by 2.2e-10. Results that
convert between a junction's inductance and its critical current, or
count a flux in flux quanta, move by about the same amount; a strongly
nonlinear operating point can amplify the change.

## Quantum optics decompositions

`autonne_takagi` returns the named tuple `(Λ, W)` for a real and a complex
matrix alike, with `Λ` in decreasing order. For a real matrix it returned
`(Λ, M)`, the unitary under the name `M`, with `Λ` in increasing order;
read the unitary as `W`, and reverse any code that relied on the order.

## SPICE inputs

The frequencies `wrspice_input_transient` and `wrspice_input_ac` take,
the sources' frequencies of the first and the sweep of the second
(`ws`, or `wstart` and `wstop`), are angular frequencies in rad/s, like
every frequency the package takes; in v0.5.4 they were in Hz. Multiply a
v0.5.4 argument by `2pi`. `wrspice_input_ac` refuses a vector of
frequencies that is not uniformly spaced, which v0.5.4 replaced by its
first and last frequencies and its length. `wrspice_input_paramp` and
`wrspice_calcS_paramp` took angular frequencies already.
`spice_hb_load`, a reader of Xyce harmonic balance output, is deprecated:
the package neither writes nor runs Xyce. It still reads the file, with a
warning.

## Network functions

- `halmos_dilation` treats a singular value above one by up to the
  square root of the machine epsilon as one by default, as `isapprox`
  does, where v0.5.4's tolerance was below the rounding of a lossless
  cascade.
- `is_unitary` and `is_orthogonal` return `false` for a matrix which is
  not square.
- `connectS` and `solveS` throw the exception of the failed frequency
  batch itself, whatever the number of batches, where several batches
  wrapped it in a `CompositeException`.
- `A_coupled_tlines` given a vector of frequencies returns an array with a
  frequency axis whatever the vector's length.
- `connectS` takes the arrays `solveS` takes, `KeyedArray`s included, and
  a covariance of an array type other than the networks'.

## Other deprecations and removals

| v0.5.4 | Now |
|---|---|
| `hbsolve(ws, wp, Ip, Nsignalmodes, Npumpmodes, circuit, circuitdefs)`, deprecated in v0.5.4 | Still deprecated; `hbsolve(ws, (wp,), sources, (Nmodulationharmonics,), (Npumpharmonics,), circuit, circuitdefs)` |
| `connectS(Sa, k, l)`, `connectS(Sa, Sb, k, l)` and their in-place forms, deprecated in v0.5.4 | Still deprecated; `intraconnectS` and `interconnectS` |
| `solveS!` without the fill reducing ordering | Deprecated; call it with the whole tuple `solveS_initialize` returns |
| `export_netlist`, `import_netlist` | Deprecated with the tuple netlist |
| `printsymmetries`, `visualizefreqs`, `sprandsubset`, `hbmatind` | Removed, with the unexported functions of v0.5.4's parser and matrix builders |

The index matrix that `hbmatind(frequencies, truncfrequencies)[2]`
returned is
`JosephsonCircuits.hbmatindices(frequencies, JosephsonCircuits.ModeDifferences(truncfrequencies.modes))`.
OrderedCollections, which `printsymmetries` used, is no longer a direct
dependency of the package.

## From development versions after v0.5.4

Code written against the development branch between v0.5.4 and this
release passed some frequencies in Hz which are now, like every
frequency the package takes, angular frequencies in rad/s. Multiply such
an argument by `2pi`, and read a returned frequency the same way. These
functions are not in v0.5.4 and have no deprecation, so a forgotten
`2pi` gives a result at the wrong frequency rather than an error.

| Function | Frequencies now in rad/s |
|---|---|
| [`RationalScattering`](@ref) of a `ScatteringParameters` block, given `npoles` or `tol` | `frequencies`, the samples it fits |
| [`RationalScattering`](@ref) of a pumped `LinearizedScattering` block | `frequencies` and `band = (wlo, whi)` |
| [`transientdemodulate`](@ref), [`transientiqplan`](@ref) | the carriers, and a plan's `frequencies`, `bandwidth3db` and `noisebandwidth` |
| [`transientquantumplan`](@ref) | the modes' frequencies, and the plan's `frequencies` |
| [`transientnoise`](@ref) | the bath frequencies, the quadrature weights and `cutoff`, and the frequencies and weights it returns |
