# API reference

The entries below are grouped by task. For units and array conventions,
see [Conventions](conventions.md); for a runnable calculation, see the
[quickstart](quickstart.md). Qualified interfaces such as
`JosephsonCircuits.@params` and `JosephsonCircuits.reset!` are listed with
the related exported functions.

```@contents
Pages = ["reference.md"]
Depth = 2
```

## Circuit construction

Define topology and bind component values. Start with the [circuit guide](circuits.md).

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["Circuit", "Interface", "Instance", "Ground", "Net", "PortRef", "PinRef", "compile", "elaborate", "ElaboratedCircuit", "CompiledCircuit", "symbolicmatrices", "numericmatrices", "calccircuitgraph", "IctoLj", "LjtoIc", "@params", "ComponentNotSupportedError"]))
```

## Components and scattering models

Component values, port models, rational fits, and their noise and DC assumptions.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["Inductor", "Capacitor", "Resistor", "CurrentSource", "VoltageSource", "Port", "MatchedTermination", "MutualInductor", "JosephsonJunction", "NonlinearInductor", "PolynomialCPR", "FrequencyDependent", "ScatteringParameters", "GaussianChannel", "TransmissionLine", "RationalScattering", "VectorFitting", "PassivityEnforcement", "LinearizedScattering", "Passive", "Lossless", "ScatteringLimit", "OpenDC", "ShortDC", "ThroughDC", "ScatteringDC", "ThermalEquilibrium", "NoiseCovariance", "ConjugateSymmetry", "Native", "passivityassessment"]))
```

## Harmonic-balance analyses

Solve operating points and signal responses, reuse setup, and differentiate outputs.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["hbsolve", "hbnlsolve", "hblinsolve", "NonlinearHB", "LinearizedHB", "HB", "designsensitivities", "designjacobian", "hbcache", "hbsolve!", "HBCache", "HBReuse", "reset!"]))
```

## Nonlinear and linear solver options

Methods, preconditioners, refresh policies, and factorization choices.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["NewtonKrylov", "Newton", "QuasiNewton", "ExternalSolver", "GMRES", "KrylovJL", "Staged", "BlockDiagonal", "FullJacobian", "HarmonicBand", "MeasuredBand", "Clusters", "CoupledModes", "CouplingMask", "Automatic", "Floquet", "Always", "Probe", "Never", "Backtracking", "KLUfactorization", "LUfactorization", "QRfactorization", "CUDSSFactorization", "BlockFactorization"]))
```

## External solver interface

Residuals and derivatives for [solver integration](interop.md).

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["HBNonlinearProblem", "hbnonlinearproblem", "JacobianOperator", "preconditioner", "hbresidual!", "hbjvp!", "hbvjp!", "hbjacobian!", "hbd2F!", "hbd3F!", "hbdFdp!", "jacobianprototype", "setdrive!", "drivenresidual!"]))
```

## Transient analyses and responses

Time-domain problems, integration rules, saved states, and response derivatives.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["TransientSource", "TransientState", "transientproblem", "transientstate", "transientsolve", "transientsensitivity", "transientdemodulate", "transienttangent", "transientadjoint", "transientinjection", "Trapezoidal", "GaussLegendre", "BackwardEuler", "WRspice", "TransientReuse", "TransientSolution", "TransientBatchSolution", "TransientStepError"]))
```

## Temporal measurements and noise

IQ measurements, canonical temporal modes, baths, noise diagnostics, and the conversions between occupations, temperatures and noise temperatures.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["transientiqplan", "transientiq!", "transientiq", "transientiqvjp!", "transientquantumplan", "transientquantum", "transientquantum!", "transientquantumvjp!", "transientnoisebaths", "transientnoise", "transientgain", "transientquantumdiagnostics", "transientquantumefficiency", "thermaloccupation", "effectivetemperature", "noisetemperature", "noisequanta"]))
```

## Network operations and constants

Connect scattering networks, change quadrature representation, and inspect physical constants.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["connectS", "solveS", "quadraturetransform", "phi0", "Phi0", "speed_of_light", "planck_constant", "reduced_planck_constant", "boltzmann_constant"]))
```

## Additional interfaces and internals

The remaining docstrings include lower-level utilities and internal types.
They are retained here for source readers and existing links; their
presence does not imply the same compatibility guarantees as a public
interface. See [Implementation notes](implementation.md) for a source map.

```@autodocs
Modules = [JosephsonCircuits, JosephsonCircuits.CircuitValues]
Filter = obj -> !any(name -> isdefined(JosephsonCircuits, name) && getfield(JosephsonCircuits, name) === obj, Symbol.(["Circuit", "Interface", "Instance", "Ground", "Net", "PortRef", "PinRef", "compile", "elaborate", "ElaboratedCircuit", "CompiledCircuit", "symbolicmatrices", "numericmatrices", "calccircuitgraph", "IctoLj", "LjtoIc", "@params", "ComponentNotSupportedError", "Inductor", "Capacitor", "Resistor", "CurrentSource", "VoltageSource", "Port", "MatchedTermination", "MutualInductor", "JosephsonJunction", "NonlinearInductor", "PolynomialCPR", "FrequencyDependent", "ScatteringParameters", "GaussianChannel", "TransmissionLine", "RationalScattering", "VectorFitting", "PassivityEnforcement", "LinearizedScattering", "Passive", "Lossless", "ScatteringLimit", "OpenDC", "ShortDC", "ThroughDC", "ScatteringDC", "ThermalEquilibrium", "NoiseCovariance", "ConjugateSymmetry", "Native", "passivityassessment", "hbsolve", "hbnlsolve", "hblinsolve", "NonlinearHB", "LinearizedHB", "HB", "designsensitivities", "designjacobian", "hbcache", "hbsolve!", "HBCache", "HBReuse", "reset!", "NewtonKrylov", "Newton", "QuasiNewton", "ExternalSolver", "GMRES", "KrylovJL", "Staged", "BlockDiagonal", "FullJacobian", "HarmonicBand", "MeasuredBand", "Clusters", "CoupledModes", "CouplingMask", "Automatic", "Floquet", "Always", "Probe", "Never", "Backtracking", "KLUfactorization", "LUfactorization", "QRfactorization", "CUDSSFactorization", "BlockFactorization", "HBNonlinearProblem", "hbnonlinearproblem", "JacobianOperator", "preconditioner", "hbresidual!", "hbjvp!", "hbvjp!", "hbjacobian!", "hbd2F!", "hbd3F!", "hbdFdp!", "jacobianprototype", "setdrive!", "drivenresidual!", "TransientSource", "TransientState", "transientproblem", "transientstate", "transientsolve", "transientsensitivity", "transientdemodulate", "transienttangent", "transientadjoint", "transientinjection", "Trapezoidal", "GaussLegendre", "BackwardEuler", "WRspice", "TransientReuse", "TransientSolution", "TransientBatchSolution", "TransientStepError", "transientiqplan", "transientiq!", "transientiq", "transientiqvjp!", "transientquantumplan", "transientquantum", "transientquantum!", "transientquantumvjp!", "transientnoisebaths", "transientnoise", "transientgain", "transientquantumdiagnostics", "transientquantumefficiency", "thermaloccupation", "effectivetemperature", "noisetemperature", "noisequanta", "connectS", "solveS", "quadraturetransform", "phi0", "Phi0", "speed_of_light", "planck_constant", "reduced_planck_constant", "boltzmann_constant"]))
```
