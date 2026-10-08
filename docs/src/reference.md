# API reference

Choose a topic for its public interfaces. For units and array conventions,
see [Conventions](conventions.md); for a runnable calculation, start with
the [quickstart](quickstart.md). Qualified interfaces such as
`JosephsonCircuits.@params` and `JosephsonCircuits.reset!` are listed with
the related exported functions.

- [Circuit construction](api/circuits.md)
- [Components and scattering models](api/components.md)
- [Harmonic-balance analyses](api/harmonicbalance.md)
- [Nonlinear and linear solver options](api/solvers.md)
- [External solver interface](api/interop.md)
- [Transient analyses and responses](api/transient.md)
- [Temporal measurements and noise](api/noise.md)
- [Network operations and constants](api/networks.md)

The [developer appendix](api/internals.md) contains the remaining internal
and lower-level interfaces, separate from the public reference.

## Existing reference links

Older bookmarks into this page still land at the corresponding entry
below. Follow its link to the full documentation on the new topic page.
These links preserve method-specific anchors as well as function names.
New entries belong on the topic pages; this compatibility index is a
snapshot of the reference before it was split.

## Circuit construction

[Open the full reference for this topic](api/circuits.md).

```@raw html
<a id="JosephsonCircuits.Ground"></a>
```

[`JosephsonCircuits.Ground`](api/circuits.md#JosephsonCircuits.Ground)

```@raw html
<a id="JosephsonCircuits.Circuit"></a>
```

[`JosephsonCircuits.Circuit`](api/circuits.md#JosephsonCircuits.Circuit)

```@raw html
<a id="JosephsonCircuits.Circuit-Tuple{AbstractVector, AbstractDict}"></a>
```

[`JosephsonCircuits.Circuit`](api/circuits.md#JosephsonCircuits.Circuit-Tuple%7BAbstractVector%2C%20AbstractDict%7D)

```@raw html
<a id="JosephsonCircuits.Circuit-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.Circuit`](api/circuits.md#JosephsonCircuits.Circuit-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.CompiledCircuit"></a>
```

[`JosephsonCircuits.CompiledCircuit`](api/circuits.md#JosephsonCircuits.CompiledCircuit)

```@raw html
<a id="JosephsonCircuits.ComponentNotSupportedError"></a>
```

[`JosephsonCircuits.ComponentNotSupportedError`](api/circuits.md#JosephsonCircuits.ComponentNotSupportedError)

```@raw html
<a id="JosephsonCircuits.ElaboratedCircuit"></a>
```

[`JosephsonCircuits.ElaboratedCircuit`](api/circuits.md#JosephsonCircuits.ElaboratedCircuit)

```@raw html
<a id="JosephsonCircuits.Instance"></a>
```

[`JosephsonCircuits.Instance`](api/circuits.md#JosephsonCircuits.Instance)

```@raw html
<a id="JosephsonCircuits.Interface"></a>
```

[`JosephsonCircuits.Interface`](api/circuits.md#JosephsonCircuits.Interface)

```@raw html
<a id="JosephsonCircuits.Net"></a>
```

[`JosephsonCircuits.Net`](api/circuits.md#JosephsonCircuits.Net)

```@raw html
<a id="JosephsonCircuits.PinRef"></a>
```

[`JosephsonCircuits.PinRef`](api/circuits.md#JosephsonCircuits.PinRef)

```@raw html
<a id="JosephsonCircuits.PortRef"></a>
```

[`JosephsonCircuits.PortRef`](api/circuits.md#JosephsonCircuits.PortRef)

```@raw html
<a id="JosephsonCircuits.IctoLj-Tuple{Any}"></a>
```

[`JosephsonCircuits.IctoLj`](api/circuits.md#JosephsonCircuits.IctoLj-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.LjtoIc-Tuple{Any}"></a>
```

[`JosephsonCircuits.LjtoIc`](api/circuits.md#JosephsonCircuits.LjtoIc-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.calccircuitgraph-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.calccircuitgraph`](api/circuits.md#JosephsonCircuits.calccircuitgraph-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.compile-Tuple{ElaboratedCircuit}"></a>
```

[`JosephsonCircuits.compile`](api/circuits.md#JosephsonCircuits.compile-Tuple%7BElaboratedCircuit%7D)

```@raw html
<a id="JosephsonCircuits.elaborate-Tuple{Circuit}"></a>
```

[`JosephsonCircuits.elaborate`](api/circuits.md#JosephsonCircuits.elaborate-Tuple%7BCircuit%7D)

```@raw html
<a id="JosephsonCircuits.numericmatrices-Tuple{Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, AbstractDict}"></a>
```

[`JosephsonCircuits.numericmatrices`](api/circuits.md#JosephsonCircuits.numericmatrices-Tuple%7BUnion%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20AbstractDict%7D)

```@raw html
<a id="JosephsonCircuits.symbolicmatrices-Tuple{Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}}"></a>
```

[`JosephsonCircuits.symbolicmatrices`](api/circuits.md#JosephsonCircuits.symbolicmatrices-Tuple%7BUnion%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%7D)

```@raw html
<a id="JosephsonCircuits.CircuitValues.@params-Tuple"></a>
```

[`JosephsonCircuits.CircuitValues.@params`](api/circuits.md#JosephsonCircuits.CircuitValues.%40params-Tuple)

## Components and scattering models

[Open the full reference for this topic](api/components.md).

```@raw html
<a id="JosephsonCircuits.Capacitor"></a>
```

[`JosephsonCircuits.Capacitor`](api/components.md#JosephsonCircuits.Capacitor)

```@raw html
<a id="JosephsonCircuits.ConjugateSymmetry"></a>
```

[`JosephsonCircuits.ConjugateSymmetry`](api/components.md#JosephsonCircuits.ConjugateSymmetry)

```@raw html
<a id="JosephsonCircuits.CurrentSource"></a>
```

[`JosephsonCircuits.CurrentSource`](api/components.md#JosephsonCircuits.CurrentSource)

```@raw html
<a id="JosephsonCircuits.GaussianChannel"></a>
```

[`JosephsonCircuits.GaussianChannel`](api/components.md#JosephsonCircuits.GaussianChannel)

```@raw html
<a id="JosephsonCircuits.Inductor"></a>
```

[`JosephsonCircuits.Inductor`](api/components.md#JosephsonCircuits.Inductor)

```@raw html
<a id="JosephsonCircuits.LaplaceResponse"></a>
```

[`JosephsonCircuits.LaplaceResponse`](api/components.md#JosephsonCircuits.LaplaceResponse)

```@raw html
<a id="JosephsonCircuits.LinearizedScattering"></a>
```

[`JosephsonCircuits.LinearizedScattering`](api/components.md#JosephsonCircuits.LinearizedScattering)

```@raw html
<a id="JosephsonCircuits.Lossless"></a>
```

[`JosephsonCircuits.Lossless`](api/components.md#JosephsonCircuits.Lossless)

```@raw html
<a id="JosephsonCircuits.MatchedTermination"></a>
```

[`JosephsonCircuits.MatchedTermination`](api/components.md#JosephsonCircuits.MatchedTermination)

```@raw html
<a id="JosephsonCircuits.MutualInductor"></a>
```

[`JosephsonCircuits.MutualInductor`](api/components.md#JosephsonCircuits.MutualInductor)

```@raw html
<a id="JosephsonCircuits.Native"></a>
```

[`JosephsonCircuits.Native`](api/components.md#JosephsonCircuits.Native)

```@raw html
<a id="JosephsonCircuits.NoiseCovariance"></a>
```

[`JosephsonCircuits.NoiseCovariance`](api/components.md#JosephsonCircuits.NoiseCovariance)

```@raw html
<a id="JosephsonCircuits.NonlinearInductor"></a>
```

[`JosephsonCircuits.NonlinearInductor`](api/components.md#JosephsonCircuits.NonlinearInductor)

```@raw html
<a id="JosephsonCircuits.OpenDC"></a>
```

[`JosephsonCircuits.OpenDC`](api/components.md#JosephsonCircuits.OpenDC)

```@raw html
<a id="JosephsonCircuits.Passive"></a>
```

[`JosephsonCircuits.Passive`](api/components.md#JosephsonCircuits.Passive)

```@raw html
<a id="JosephsonCircuits.PassivityEnforcement"></a>
```

[`JosephsonCircuits.PassivityEnforcement`](api/components.md#JosephsonCircuits.PassivityEnforcement)

```@raw html
<a id="JosephsonCircuits.PolynomialCPR"></a>
```

[`JosephsonCircuits.PolynomialCPR`](api/components.md#JosephsonCircuits.PolynomialCPR)

```@raw html
<a id="JosephsonCircuits.Port"></a>
```

[`JosephsonCircuits.Port`](api/components.md#JosephsonCircuits.Port)

```@raw html
<a id="JosephsonCircuits.Resistor"></a>
```

[`JosephsonCircuits.Resistor`](api/components.md#JosephsonCircuits.Resistor)

```@raw html
<a id="JosephsonCircuits.ScatteringDC"></a>
```

[`JosephsonCircuits.ScatteringDC`](api/components.md#JosephsonCircuits.ScatteringDC)

```@raw html
<a id="JosephsonCircuits.ScatteringLimit"></a>
```

[`JosephsonCircuits.ScatteringLimit`](api/components.md#JosephsonCircuits.ScatteringLimit)

```@raw html
<a id="JosephsonCircuits.ScatteringParameters"></a>
```

[`JosephsonCircuits.ScatteringParameters`](api/components.md#JosephsonCircuits.ScatteringParameters)

```@raw html
<a id="JosephsonCircuits.ShortDC"></a>
```

[`JosephsonCircuits.ShortDC`](api/components.md#JosephsonCircuits.ShortDC)

```@raw html
<a id="JosephsonCircuits.ThermalEquilibrium"></a>
```

[`JosephsonCircuits.ThermalEquilibrium`](api/components.md#JosephsonCircuits.ThermalEquilibrium)

```@raw html
<a id="JosephsonCircuits.ThroughDC"></a>
```

[`JosephsonCircuits.ThroughDC`](api/components.md#JosephsonCircuits.ThroughDC)

```@raw html
<a id="JosephsonCircuits.VectorFitting"></a>
```

[`JosephsonCircuits.VectorFitting`](api/components.md#JosephsonCircuits.VectorFitting)

```@raw html
<a id="JosephsonCircuits.VoltageSource"></a>
```

[`JosephsonCircuits.VoltageSource`](api/components.md#JosephsonCircuits.VoltageSource)

```@raw html
<a id="JosephsonCircuits.FrequencyDependent-Tuple{Any}"></a>
```

[`JosephsonCircuits.FrequencyDependent`](api/components.md#JosephsonCircuits.FrequencyDependent-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.JosephsonJunction-Tuple{Any}"></a>
```

[`JosephsonCircuits.JosephsonJunction`](api/components.md#JosephsonCircuits.JosephsonJunction-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.RationalScattering-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.RationalScattering`](api/components.md#JosephsonCircuits.RationalScattering-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.RationalScattering-Tuple{LinearizedScattering, Integer}"></a>
```

[`JosephsonCircuits.RationalScattering`](api/components.md#JosephsonCircuits.RationalScattering-Tuple%7BLinearizedScattering%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.RationalScattering-Tuple{ScatteringParameters, Integer}"></a>
```

[`JosephsonCircuits.RationalScattering`](api/components.md#JosephsonCircuits.RationalScattering-Tuple%7BScatteringParameters%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.RationalScattering-Tuple{ScatteringParameters}"></a>
```

[`JosephsonCircuits.RationalScattering`](api/components.md#JosephsonCircuits.RationalScattering-Tuple%7BScatteringParameters%7D)

```@raw html
<a id="JosephsonCircuits.TransmissionLine-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.TransmissionLine`](api/components.md#JosephsonCircuits.TransmissionLine-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.passivityassessment-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.passivityassessment`](api/components.md#JosephsonCircuits.passivityassessment-NTuple%7B4%2C%20Any%7D)

## Harmonic-balance analyses

[Open the full reference for this topic](api/harmonicbalance.md).

```@raw html
<a id="JosephsonCircuits.ContourIntegral"></a>
```

[`JosephsonCircuits.ContourIntegral`](api/harmonicbalance.md#JosephsonCircuits.ContourIntegral)

```@raw html
<a id="JosephsonCircuits.DenseSpectrum"></a>
```

[`JosephsonCircuits.DenseSpectrum`](api/harmonicbalance.md#JosephsonCircuits.DenseSpectrum)

```@raw html
<a id="JosephsonCircuits.HB"></a>
```

[`JosephsonCircuits.HB`](api/harmonicbalance.md#JosephsonCircuits.HB)

```@raw html
<a id="JosephsonCircuits.HBCache"></a>
```

[`JosephsonCircuits.HBCache`](api/harmonicbalance.md#JosephsonCircuits.HBCache)

```@raw html
<a id="JosephsonCircuits.HBReuse"></a>
```

[`JosephsonCircuits.HBReuse`](api/harmonicbalance.md#JosephsonCircuits.HBReuse)

```@raw html
<a id="JosephsonCircuits.HBStabilityResult"></a>
```

[`JosephsonCircuits.HBStabilityResult`](api/harmonicbalance.md#JosephsonCircuits.HBStabilityResult)

```@raw html
<a id="JosephsonCircuits.LinearizedHB"></a>
```

[`JosephsonCircuits.LinearizedHB`](api/harmonicbalance.md#JosephsonCircuits.LinearizedHB)

```@raw html
<a id="JosephsonCircuits.Monodromy"></a>
```

[`JosephsonCircuits.Monodromy`](api/harmonicbalance.md#JosephsonCircuits.Monodromy)

```@raw html
<a id="JosephsonCircuits.NonlinearHB"></a>
```

[`JosephsonCircuits.NonlinearHB`](api/harmonicbalance.md#JosephsonCircuits.NonlinearHB)

```@raw html
<a id="JosephsonCircuits.ShiftInvert"></a>
```

[`JosephsonCircuits.ShiftInvert`](api/harmonicbalance.md#JosephsonCircuits.ShiftInvert)

```@raw html
<a id="JosephsonCircuits.designjacobian-Tuple{Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, AbstractDict}"></a>
```

[`JosephsonCircuits.designjacobian`](api/harmonicbalance.md#JosephsonCircuits.designjacobian-Tuple%7BUnion%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20AbstractDict%7D)

```@raw html
<a id="JosephsonCircuits.designsensitivities-Tuple{Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, AbstractDict, Vararg{Any, 5}}"></a>
```

[`JosephsonCircuits.designsensitivities`](api/harmonicbalance.md#JosephsonCircuits.designsensitivities-Tuple%7BUnion%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20AbstractDict%2C%20Vararg%7BAny%2C%205%7D%7D)

```@raw html
<a id="JosephsonCircuits.hbcache-Union{Tuple{N}, Tuple{NTuple{N, Number}, NTuple{N, Int64}, Any, Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}}, Tuple{NTuple{N, Number}, NTuple{N, Int64}, Any, Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, AbstractDict}} where N"></a>
```

[`JosephsonCircuits.hbcache`](api/harmonicbalance.md#JosephsonCircuits.hbcache-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Number%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Any%2C%20Union%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%7D%2C%20Tuple%7BNTuple%7BN%2C%20Number%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Any%2C%20Union%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20AbstractDict%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hblinsolve"></a>
```

[`JosephsonCircuits.hblinsolve`](api/harmonicbalance.md#JosephsonCircuits.hblinsolve)

```@raw html
<a id="JosephsonCircuits.hblinsolve-Tuple{Any, JosephsonCircuits.CompiledCircuit, Union{AbstractDict, AbstractVector}, JosephsonCircuits.Frequencies}"></a>
```

[`JosephsonCircuits.hblinsolve`](api/harmonicbalance.md#JosephsonCircuits.hblinsolve-Tuple%7BAny%2C%20JosephsonCircuits.CompiledCircuit%2C%20Union%7BAbstractDict%2C%20AbstractVector%7D%2C%20JosephsonCircuits.Frequencies%7D)

```@raw html
<a id="JosephsonCircuits.hbnlsolve-Union{Tuple{N}, Tuple{NTuple{N, Float64}, NTuple{N, Int64}, Array{@NamedTuple{mode::NTuple{N, Int64}, port::Int64, current::ComplexF64}, 1}, JosephsonCircuits.CompiledCircuit, Dict{Any, Any}}} where N"></a>
```

[`JosephsonCircuits.hbnlsolve`](api/harmonicbalance.md#JosephsonCircuits.hbnlsolve-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Float64%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Array%7B%40NamedTuple%7Bmode%3A%3ANTuple%7BN%2C%20Int64%7D%2C%20port%3A%3AInt64%2C%20current%3A%3AComplexF64%7D%2C%201%7D%2C%20JosephsonCircuits.CompiledCircuit%2C%20Dict%7BAny%2C%20Any%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hbnlsolve-Union{Tuple{N}, Tuple{NTuple{N, Number}, Any, JosephsonCircuits.Frequencies{N}, JosephsonCircuits.FourierIndices{N}, JosephsonCircuits.CompiledCircuit, JosephsonCircuits.CircuitMatrices}} where N"></a>
```

[`JosephsonCircuits.hbnlsolve`](api/harmonicbalance.md#JosephsonCircuits.hbnlsolve-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Number%7D%2C%20Any%2C%20JosephsonCircuits.Frequencies%7BN%7D%2C%20JosephsonCircuits.FourierIndices%7BN%7D%2C%20JosephsonCircuits.CompiledCircuit%2C%20JosephsonCircuits.CircuitMatrices%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hbnlsolve-Union{Tuple{N}, Tuple{NTuple{N, Number}, NTuple{N, Int64}, Any, Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}}, Tuple{NTuple{N, Number}, NTuple{N, Int64}, Any, Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, AbstractDict}} where N"></a>
```

[`JosephsonCircuits.hbnlsolve`](api/harmonicbalance.md#JosephsonCircuits.hbnlsolve-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Number%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Any%2C%20Union%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%7D%2C%20Tuple%7BNTuple%7BN%2C%20Number%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Any%2C%20Union%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20AbstractDict%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hbsolve!-Tuple{HBCache, NamedTuple}"></a>
```

[`JosephsonCircuits.hbsolve!`](api/harmonicbalance.md#JosephsonCircuits.hbsolve!-Tuple%7BHBCache%2C%20NamedTuple%7D)

```@raw html
<a id="JosephsonCircuits.hbsolve-Union{Tuple{M}, Tuple{N}, Tuple{Any, NTuple{N, Number}, Any, NTuple{M, Int64}, NTuple{N, Int64}, Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}}, Tuple{Any, NTuple{N, Number}, Any, NTuple{M, Int64}, NTuple{N, Int64}, Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, AbstractDict}} where {N, M}"></a>
```

[`JosephsonCircuits.hbsolve`](api/harmonicbalance.md#JosephsonCircuits.hbsolve-Union%7BTuple%7BM%7D%2C%20Tuple%7BN%7D%2C%20Tuple%7BAny%2C%20NTuple%7BN%2C%20Number%7D%2C%20Any%2C%20NTuple%7BM%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Union%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%7D%2C%20Tuple%7BAny%2C%20NTuple%7BN%2C%20Number%7D%2C%20Any%2C%20NTuple%7BM%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Union%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20AbstractDict%7D%7D%20where%20%7BN%2C%20M%7D)

```@raw html
<a id="JosephsonCircuits.hbsolve-Union{Tuple{M}, Tuple{N}, Tuple{Vector{Float64}, NTuple{N, Float64}, Array{@NamedTuple{mode::NTuple{N, Int64}, port::Int64, current::ComplexF64}, 1}, NTuple{M, Int64}, NTuple{N, Int64}, JosephsonCircuits.CompiledCircuit, Dict{Any, Any}}} where {N, M}"></a>
```

[`JosephsonCircuits.hbsolve`](api/harmonicbalance.md#JosephsonCircuits.hbsolve-Union%7BTuple%7BM%7D%2C%20Tuple%7BN%7D%2C%20Tuple%7BVector%7BFloat64%7D%2C%20NTuple%7BN%2C%20Float64%7D%2C%20Array%7B%40NamedTuple%7Bmode%3A%3ANTuple%7BN%2C%20Int64%7D%2C%20port%3A%3AInt64%2C%20current%3A%3AComplexF64%7D%2C%201%7D%2C%20NTuple%7BM%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20JosephsonCircuits.CompiledCircuit%2C%20Dict%7BAny%2C%20Any%7D%7D%7D%20where%20%7BN%2C%20M%7D)

```@raw html
<a id="JosephsonCircuits.hbstability"></a>
```

[`JosephsonCircuits.hbstability`](api/harmonicbalance.md#JosephsonCircuits.hbstability)

```@raw html
<a id="JosephsonCircuits.matchpoles-Tuple{HBStabilityResult, HBStabilityResult}"></a>
```

[`JosephsonCircuits.matchpoles`](api/harmonicbalance.md#JosephsonCircuits.matchpoles-Tuple%7BHBStabilityResult%2C%20HBStabilityResult%7D)

```@raw html
<a id="JosephsonCircuits.reset!-Tuple{HBCache}"></a>
```

[`JosephsonCircuits.reset!`](api/harmonicbalance.md#JosephsonCircuits.reset!-Tuple%7BHBCache%7D)

## Nonlinear and linear solver options

[Open the full reference for this topic](api/solvers.md).

```@raw html
<a id="JosephsonCircuits.Always"></a>
```

[`JosephsonCircuits.Always`](api/solvers.md#JosephsonCircuits.Always)

```@raw html
<a id="JosephsonCircuits.Automatic"></a>
```

[`JosephsonCircuits.Automatic`](api/solvers.md#JosephsonCircuits.Automatic)

```@raw html
<a id="JosephsonCircuits.Backtracking"></a>
```

[`JosephsonCircuits.Backtracking`](api/solvers.md#JosephsonCircuits.Backtracking)

```@raw html
<a id="JosephsonCircuits.BlockDiagonal"></a>
```

[`JosephsonCircuits.BlockDiagonal`](api/solvers.md#JosephsonCircuits.BlockDiagonal)

```@raw html
<a id="JosephsonCircuits.BlockFactorization"></a>
```

[`JosephsonCircuits.BlockFactorization`](api/solvers.md#JosephsonCircuits.BlockFactorization)

```@raw html
<a id="JosephsonCircuits.CUDSSFactorization"></a>
```

[`JosephsonCircuits.CUDSSFactorization`](api/solvers.md#JosephsonCircuits.CUDSSFactorization)

```@raw html
<a id="JosephsonCircuits.Clusters"></a>
```

[`JosephsonCircuits.Clusters`](api/solvers.md#JosephsonCircuits.Clusters)

```@raw html
<a id="JosephsonCircuits.CoupledModes"></a>
```

[`JosephsonCircuits.CoupledModes`](api/solvers.md#JosephsonCircuits.CoupledModes)

```@raw html
<a id="JosephsonCircuits.CouplingMask"></a>
```

[`JosephsonCircuits.CouplingMask`](api/solvers.md#JosephsonCircuits.CouplingMask)

```@raw html
<a id="JosephsonCircuits.ExternalSolver"></a>
```

[`JosephsonCircuits.ExternalSolver`](api/solvers.md#JosephsonCircuits.ExternalSolver)

```@raw html
<a id="JosephsonCircuits.Floquet"></a>
```

[`JosephsonCircuits.Floquet`](api/solvers.md#JosephsonCircuits.Floquet)

```@raw html
<a id="JosephsonCircuits.FullJacobian"></a>
```

[`JosephsonCircuits.FullJacobian`](api/solvers.md#JosephsonCircuits.FullJacobian)

```@raw html
<a id="JosephsonCircuits.GMRES"></a>
```

[`JosephsonCircuits.GMRES`](api/solvers.md#JosephsonCircuits.GMRES)

```@raw html
<a id="JosephsonCircuits.HarmonicBand"></a>
```

[`JosephsonCircuits.HarmonicBand`](api/solvers.md#JosephsonCircuits.HarmonicBand)

```@raw html
<a id="JosephsonCircuits.KLUfactorization"></a>
```

[`JosephsonCircuits.KLUfactorization`](api/solvers.md#JosephsonCircuits.KLUfactorization)

```@raw html
<a id="JosephsonCircuits.KrylovJL"></a>
```

[`JosephsonCircuits.KrylovJL`](api/solvers.md#JosephsonCircuits.KrylovJL)

```@raw html
<a id="JosephsonCircuits.LUfactorization"></a>
```

[`JosephsonCircuits.LUfactorization`](api/solvers.md#JosephsonCircuits.LUfactorization)

```@raw html
<a id="JosephsonCircuits.MeasuredBand"></a>
```

[`JosephsonCircuits.MeasuredBand`](api/solvers.md#JosephsonCircuits.MeasuredBand)

```@raw html
<a id="JosephsonCircuits.Never"></a>
```

[`JosephsonCircuits.Never`](api/solvers.md#JosephsonCircuits.Never)

```@raw html
<a id="JosephsonCircuits.Newton"></a>
```

[`JosephsonCircuits.Newton`](api/solvers.md#JosephsonCircuits.Newton)

```@raw html
<a id="JosephsonCircuits.NewtonKrylov"></a>
```

[`JosephsonCircuits.NewtonKrylov`](api/solvers.md#JosephsonCircuits.NewtonKrylov)

```@raw html
<a id="JosephsonCircuits.Probe"></a>
```

[`JosephsonCircuits.Probe`](api/solvers.md#JosephsonCircuits.Probe)

```@raw html
<a id="JosephsonCircuits.QRfactorization"></a>
```

[`JosephsonCircuits.QRfactorization`](api/solvers.md#JosephsonCircuits.QRfactorization)

```@raw html
<a id="JosephsonCircuits.QuasiNewton"></a>
```

[`JosephsonCircuits.QuasiNewton`](api/solvers.md#JosephsonCircuits.QuasiNewton)

```@raw html
<a id="JosephsonCircuits.Staged"></a>
```

[`JosephsonCircuits.Staged`](api/solvers.md#JosephsonCircuits.Staged)

## External solver interface

[Open the full reference for this topic](api/interop.md).

```@raw html
<a id="JosephsonCircuits.HBNonlinearProblem"></a>
```

[`JosephsonCircuits.HBNonlinearProblem`](api/interop.md#JosephsonCircuits.HBNonlinearProblem)

```@raw html
<a id="JosephsonCircuits.JacobianOperator"></a>
```

[`JosephsonCircuits.JacobianOperator`](api/interop.md#JosephsonCircuits.JacobianOperator)

```@raw html
<a id="JosephsonCircuits.drivenresidual!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.HBNonlinearProblem, AbstractVector{&lt;:Real}, Real}"></a>
```

[`JosephsonCircuits.drivenresidual!`](api/interop.md#JosephsonCircuits.drivenresidual!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20Real%7D)

```@raw html
<a id="JosephsonCircuits.hbd2F!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.HBNonlinearProblem, AbstractVector{&lt;:Real}, AbstractVector{&lt;:Real}, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.hbd2F!`](api/interop.md#JosephsonCircuits.hbd2F!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.hbd3F!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.HBNonlinearProblem, AbstractVector{&lt;:Real}, AbstractVector{&lt;:Real}, AbstractVector{&lt;:Real}, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.hbd3F!`](api/interop.md#JosephsonCircuits.hbd3F!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.hbdFdp!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.HBNonlinearProblem}"></a>
```

[`JosephsonCircuits.hbdFdp!`](api/interop.md#JosephsonCircuits.hbdFdp!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.HBNonlinearProblem%7D)

```@raw html
<a id="JosephsonCircuits.hbjacobian!-Tuple{Any, JosephsonCircuits.HBNonlinearProblem, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.hbjacobian!`](api/interop.md#JosephsonCircuits.hbjacobian!-Tuple%7BAny%2C%20JosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.hbjvp!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.HBNonlinearProblem, AbstractVector{&lt;:Real}, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.hbjvp!`](api/interop.md#JosephsonCircuits.hbjvp!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.hbnonlinearproblem"></a>
```

[`JosephsonCircuits.hbnonlinearproblem`](api/interop.md#JosephsonCircuits.hbnonlinearproblem)

```@raw html
<a id="JosephsonCircuits.hbresidual!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.HBNonlinearProblem, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.hbresidual!`](api/interop.md#JosephsonCircuits.hbresidual!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.hbvjp!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.HBNonlinearProblem, AbstractVector{&lt;:Real}, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.hbvjp!`](api/interop.md#JosephsonCircuits.hbvjp!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7B%3C%3AReal%7D%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.jacobianprototype-Tuple{JosephsonCircuits.HBNonlinearProblem}"></a>
```

[`JosephsonCircuits.jacobianprototype`](api/interop.md#JosephsonCircuits.jacobianprototype-Tuple%7BJosephsonCircuits.HBNonlinearProblem%7D)

```@raw html
<a id="JosephsonCircuits.preconditioner-Tuple{JosephsonCircuits.HBNonlinearProblem, AbstractVector}"></a>
```

[`JosephsonCircuits.preconditioner`](api/interop.md#JosephsonCircuits.preconditioner-Tuple%7BJosephsonCircuits.HBNonlinearProblem%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.setdrive!-Tuple{JosephsonCircuits.HBNonlinearProblem, Real}"></a>
```

[`JosephsonCircuits.setdrive!`](api/interop.md#JosephsonCircuits.setdrive!-Tuple%7BJosephsonCircuits.HBNonlinearProblem%2C%20Real%7D)

## Transient analyses and responses

[Open the full reference for this topic](api/transient.md).

```@raw html
<a id="JosephsonCircuits.BackwardEuler"></a>
```

[`JosephsonCircuits.BackwardEuler`](api/transient.md#JosephsonCircuits.BackwardEuler)

```@raw html
<a id="JosephsonCircuits.GaussLegendre"></a>
```

[`JosephsonCircuits.GaussLegendre`](api/transient.md#JosephsonCircuits.GaussLegendre)

```@raw html
<a id="JosephsonCircuits.TransientBatchSolution"></a>
```

[`JosephsonCircuits.TransientBatchSolution`](api/transient.md#JosephsonCircuits.TransientBatchSolution)

```@raw html
<a id="JosephsonCircuits.TransientReuse"></a>
```

[`JosephsonCircuits.TransientReuse`](api/transient.md#JosephsonCircuits.TransientReuse)

```@raw html
<a id="JosephsonCircuits.TransientSolution"></a>
```

[`JosephsonCircuits.TransientSolution`](api/transient.md#JosephsonCircuits.TransientSolution)

```@raw html
<a id="JosephsonCircuits.TransientSource"></a>
```

[`JosephsonCircuits.TransientSource`](api/transient.md#JosephsonCircuits.TransientSource)

```@raw html
<a id="JosephsonCircuits.TransientState"></a>
```

[`JosephsonCircuits.TransientState`](api/transient.md#JosephsonCircuits.TransientState)

```@raw html
<a id="JosephsonCircuits.TransientStepError"></a>
```

[`JosephsonCircuits.TransientStepError`](api/transient.md#JosephsonCircuits.TransientStepError)

```@raw html
<a id="JosephsonCircuits.Trapezoidal"></a>
```

[`JosephsonCircuits.Trapezoidal`](api/transient.md#JosephsonCircuits.Trapezoidal)

```@raw html
<a id="JosephsonCircuits.WRspice"></a>
```

[`JosephsonCircuits.WRspice`](api/transient.md#JosephsonCircuits.WRspice)

```@raw html
<a id="JosephsonCircuits.transientadjoint-Tuple{TransientBatchSolution, AbstractArray{&lt;:Real}}"></a>
```

[`JosephsonCircuits.transientadjoint`](api/transient.md#JosephsonCircuits.transientadjoint-Tuple%7BTransientBatchSolution%2C%20AbstractArray%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.transientadjoint-Tuple{TransientSolution, AbstractArray{&lt;:Real}}"></a>
```

[`JosephsonCircuits.transientadjoint`](api/transient.md#JosephsonCircuits.transientadjoint-Tuple%7BTransientSolution%2C%20AbstractArray%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.transientdemodulate-Tuple{TransientSolution, Integer, Real}"></a>
```

[`JosephsonCircuits.transientdemodulate`](api/transient.md#JosephsonCircuits.transientdemodulate-Tuple%7BTransientSolution%2C%20Integer%2C%20Real%7D)

```@raw html
<a id="JosephsonCircuits.transientinjection-Tuple{JosephsonCircuits.TransientProblem, Any}"></a>
```

[`JosephsonCircuits.transientinjection`](api/transient.md#JosephsonCircuits.transientinjection-Tuple%7BJosephsonCircuits.TransientProblem%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientproblem"></a>
```

[`JosephsonCircuits.transientproblem`](api/transient.md#JosephsonCircuits.transientproblem)

```@raw html
<a id="JosephsonCircuits.transientproblem-Tuple{JosephsonCircuits.TransientProblem}"></a>
```

[`JosephsonCircuits.transientproblem`](api/transient.md#JosephsonCircuits.transientproblem-Tuple%7BJosephsonCircuits.TransientProblem%7D)

```@raw html
<a id="JosephsonCircuits.transientsensitivity-Tuple{Union{TransientBatchSolution, TransientSolution}, Any}"></a>
```

[`JosephsonCircuits.transientsensitivity`](api/transient.md#JosephsonCircuits.transientsensitivity-Tuple%7BUnion%7BTransientBatchSolution%2C%20TransientSolution%7D%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientsolve-Tuple{AbstractVector{JosephsonCircuits.TransientProblem}, Any}"></a>
```

[`JosephsonCircuits.transientsolve`](api/transient.md#JosephsonCircuits.transientsolve-Tuple%7BAbstractVector%7BJosephsonCircuits.TransientProblem%7D%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientsolve-Tuple{JosephsonCircuits.TransientProblem, Any}"></a>
```

[`JosephsonCircuits.transientsolve`](api/transient.md#JosephsonCircuits.transientsolve-Tuple%7BJosephsonCircuits.TransientProblem%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientstate-Tuple{JosephsonCircuits.TransientProblem}"></a>
```

[`JosephsonCircuits.transientstate`](api/transient.md#JosephsonCircuits.transientstate-Tuple%7BJosephsonCircuits.TransientProblem%7D)

```@raw html
<a id="JosephsonCircuits.transienttangent-Tuple{TransientBatchSolution, Union{Nothing, AbstractArray{&lt;:Real}}}"></a>
```

[`JosephsonCircuits.transienttangent`](api/transient.md#JosephsonCircuits.transienttangent-Tuple%7BTransientBatchSolution%2C%20Union%7BNothing%2C%20AbstractArray%7B%3C%3AReal%7D%7D%7D)

```@raw html
<a id="JosephsonCircuits.transienttangent-Tuple{TransientSolution, Union{Nothing, AbstractArray{&lt;:Real}}}"></a>
```

[`JosephsonCircuits.transienttangent`](api/transient.md#JosephsonCircuits.transienttangent-Tuple%7BTransientSolution%2C%20Union%7BNothing%2C%20AbstractArray%7B%3C%3AReal%7D%7D%7D)

## Temporal measurements and noise

[Open the full reference for this topic](api/noise.md).

```@raw html
<a id="JosephsonCircuits.effectivetemperature-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.effectivetemperature`](api/noise.md#JosephsonCircuits.effectivetemperature-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.noisequanta-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.noisequanta`](api/noise.md#JosephsonCircuits.noisequanta-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.noisetemperature-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.noisetemperature`](api/noise.md#JosephsonCircuits.noisetemperature-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.thermaloccupation-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.thermaloccupation`](api/noise.md#JosephsonCircuits.thermaloccupation-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientgain-Tuple{TransientSolution, JosephsonCircuits.TransientQuantumPlan, JosephsonCircuits.TransientQuantumPlan}"></a>
```

[`JosephsonCircuits.transientgain`](api/noise.md#JosephsonCircuits.transientgain-Tuple%7BTransientSolution%2C%20JosephsonCircuits.TransientQuantumPlan%2C%20JosephsonCircuits.TransientQuantumPlan%7D)

```@raw html
<a id="JosephsonCircuits.transientiq!-Tuple{Any, JosephsonCircuits.TransientIQPlan, Any}"></a>
```

[`JosephsonCircuits.transientiq!`](api/noise.md#JosephsonCircuits.transientiq!-Tuple%7BAny%2C%20JosephsonCircuits.TransientIQPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientiq-Tuple{JosephsonCircuits.TransientIQPlan, Any}"></a>
```

[`JosephsonCircuits.transientiq`](api/noise.md#JosephsonCircuits.transientiq-Tuple%7BJosephsonCircuits.TransientIQPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientiqplan-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.transientiqplan`](api/noise.md#JosephsonCircuits.transientiqplan-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientiqvjp!-Tuple{Any, JosephsonCircuits.TransientIQPlan, Any}"></a>
```

[`JosephsonCircuits.transientiqvjp!`](api/noise.md#JosephsonCircuits.transientiqvjp!-Tuple%7BAny%2C%20JosephsonCircuits.TransientIQPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientnoise-Tuple{TransientSolution, JosephsonCircuits.TransientQuantumPlan}"></a>
```

[`JosephsonCircuits.transientnoise`](api/noise.md#JosephsonCircuits.transientnoise-Tuple%7BTransientSolution%2C%20JosephsonCircuits.TransientQuantumPlan%7D)

```@raw html
<a id="JosephsonCircuits.transientnoisebaths-Tuple{JosephsonCircuits.TransientProblem}"></a>
```

[`JosephsonCircuits.transientnoisebaths`](api/noise.md#JosephsonCircuits.transientnoisebaths-Tuple%7BJosephsonCircuits.TransientProblem%7D)

```@raw html
<a id="JosephsonCircuits.transientquantum!-Tuple{Any, JosephsonCircuits.TransientQuantumPlan, Any}"></a>
```

[`JosephsonCircuits.transientquantum!`](api/noise.md#JosephsonCircuits.transientquantum!-Tuple%7BAny%2C%20JosephsonCircuits.TransientQuantumPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientquantum-Tuple{JosephsonCircuits.TransientQuantumPlan, Any}"></a>
```

[`JosephsonCircuits.transientquantum`](api/noise.md#JosephsonCircuits.transientquantum-Tuple%7BJosephsonCircuits.TransientQuantumPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientquantumdiagnostics-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.transientquantumdiagnostics`](api/noise.md#JosephsonCircuits.transientquantumdiagnostics-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientquantumefficiency-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.transientquantumefficiency`](api/noise.md#JosephsonCircuits.transientquantumefficiency-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientquantumplan-Tuple{Any, Any, AbstractMatrix}"></a>
```

[`JosephsonCircuits.transientquantumplan`](api/noise.md#JosephsonCircuits.transientquantumplan-Tuple%7BAny%2C%20Any%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.transientquantumvjp!-Tuple{Any, JosephsonCircuits.TransientQuantumPlan, Any}"></a>
```

[`JosephsonCircuits.transientquantumvjp!`](api/noise.md#JosephsonCircuits.transientquantumvjp!-Tuple%7BAny%2C%20JosephsonCircuits.TransientQuantumPlan%2C%20Any%7D)

## Network operations and constants

[Open the full reference for this topic](api/networks.md).

```@raw html
<a id="JosephsonCircuits.Phi0"></a>
```

[`JosephsonCircuits.Phi0`](api/networks.md#JosephsonCircuits.Phi0)

```@raw html
<a id="JosephsonCircuits.boltzmann_constant"></a>
```

[`JosephsonCircuits.boltzmann_constant`](api/networks.md#JosephsonCircuits.boltzmann_constant)

```@raw html
<a id="JosephsonCircuits.phi0"></a>
```

[`JosephsonCircuits.phi0`](api/networks.md#JosephsonCircuits.phi0)

```@raw html
<a id="JosephsonCircuits.planck_constant"></a>
```

[`JosephsonCircuits.planck_constant`](api/networks.md#JosephsonCircuits.planck_constant)

```@raw html
<a id="JosephsonCircuits.reduced_planck_constant"></a>
```

[`JosephsonCircuits.reduced_planck_constant`](api/networks.md#JosephsonCircuits.reduced_planck_constant)

```@raw html
<a id="JosephsonCircuits.speed_of_light"></a>
```

[`JosephsonCircuits.speed_of_light`](api/networks.md#JosephsonCircuits.speed_of_light)

```@raw html
<a id="JosephsonCircuits.connectS-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.connectS`](api/networks.md#JosephsonCircuits.connectS-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.quadraturetransform-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.quadraturetransform`](api/networks.md#JosephsonCircuits.quadraturetransform-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.solveS-Tuple{AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.solveS`](api/networks.md#JosephsonCircuits.solveS-Tuple%7BAbstractVector%2C%20AbstractVector%7D)

## Additional interfaces and internals

[Open the full reference for this topic](api/internals.md).

```@raw html
<a id="JosephsonCircuits.JosephsonCircuits"></a>
```

[`JosephsonCircuits.JosephsonCircuits`](api/internals.md#JosephsonCircuits.JosephsonCircuits)

```@raw html
<a id="JosephsonCircuits.CALLABLE_FORMS"></a>
```

[`JosephsonCircuits.CALLABLE_FORMS`](api/internals.md#JosephsonCircuits.CALLABLE_FORMS)

```@raw html
<a id="JosephsonCircuits.FCONJ"></a>
```

[`JosephsonCircuits.FCONJ`](api/internals.md#JosephsonCircuits.FCONJ)

```@raw html
<a id="JosephsonCircuits.FORWARDSENSITIVITYSTAMPBYTES"></a>
```

[`JosephsonCircuits.FORWARDSENSITIVITYSTAMPBYTES`](api/internals.md#JosephsonCircuits.FORWARDSENSITIVITYSTAMPBYTES)

```@raw html
<a id="JosephsonCircuits.FWIDE"></a>
```

[`JosephsonCircuits.FWIDE`](api/internals.md#JosephsonCircuits.FWIDE)

```@raw html
<a id="JosephsonCircuits.IMPEDANCE_C"></a>
```

[`JosephsonCircuits.IMPEDANCE_C`](api/internals.md#JosephsonCircuits.IMPEDANCE_C)

```@raw html
<a id="JosephsonCircuits.IMPEDANCE_L"></a>
```

[`JosephsonCircuits.IMPEDANCE_L`](api/internals.md#JosephsonCircuits.IMPEDANCE_L)

```@raw html
<a id="JosephsonCircuits.IMPEDANCE_R"></a>
```

[`JosephsonCircuits.IMPEDANCE_R`](api/internals.md#JosephsonCircuits.IMPEDANCE_R)

```@raw html
<a id="JosephsonCircuits.REVERSESENSITIVITYCHUNKBYTES"></a>
```

[`JosephsonCircuits.REVERSESENSITIVITYCHUNKBYTES`](api/internals.md#JosephsonCircuits.REVERSESENSITIVITYCHUNKBYTES)

```@raw html
<a id="JosephsonCircuits.AbstractComponent"></a>
```

[`JosephsonCircuits.AbstractComponent`](api/internals.md#JosephsonCircuits.AbstractComponent)

```@raw html
<a id="JosephsonCircuits.AbstractDCModel"></a>
```

[`JosephsonCircuits.AbstractDCModel`](api/internals.md#JosephsonCircuits.AbstractDCModel)

```@raw html
<a id="JosephsonCircuits.AbstractFactorization"></a>
```

[`JosephsonCircuits.AbstractFactorization`](api/internals.md#JosephsonCircuits.AbstractFactorization)

```@raw html
<a id="JosephsonCircuits.AbstractHBNonlinearSolver"></a>
```

[`JosephsonCircuits.AbstractHBNonlinearSolver`](api/internals.md#JosephsonCircuits.AbstractHBNonlinearSolver)

```@raw html
<a id="JosephsonCircuits.AbstractMatrixProvider"></a>
```

[`JosephsonCircuits.AbstractMatrixProvider`](api/internals.md#JosephsonCircuits.AbstractMatrixProvider)

```@raw html
<a id="JosephsonCircuits.AbstractModeCoupling"></a>
```

[`JosephsonCircuits.AbstractModeCoupling`](api/internals.md#JosephsonCircuits.AbstractModeCoupling)

```@raw html
<a id="JosephsonCircuits.AbstractPortTermination"></a>
```

[`JosephsonCircuits.AbstractPortTermination`](api/internals.md#JosephsonCircuits.AbstractPortTermination)

```@raw html
<a id="JosephsonCircuits.AbstractPreconditioner"></a>
```

[`JosephsonCircuits.AbstractPreconditioner`](api/internals.md#JosephsonCircuits.AbstractPreconditioner)

```@raw html
<a id="JosephsonCircuits.AbstractPreconditionerSpec"></a>
```

[`JosephsonCircuits.AbstractPreconditionerSpec`](api/internals.md#JosephsonCircuits.AbstractPreconditionerSpec)

```@raw html
<a id="JosephsonCircuits.AbstractStageInfo"></a>
```

[`JosephsonCircuits.AbstractStageInfo`](api/internals.md#JosephsonCircuits.AbstractStageInfo)

```@raw html
<a id="JosephsonCircuits.AbstractTransientIntegrator"></a>
```

[`JosephsonCircuits.AbstractTransientIntegrator`](api/internals.md#JosephsonCircuits.AbstractTransientIntegrator)

```@raw html
<a id="JosephsonCircuits.AbstractWrappedPreconditioner"></a>
```

[`JosephsonCircuits.AbstractWrappedPreconditioner`](api/internals.md#JosephsonCircuits.AbstractWrappedPreconditioner)

```@raw html
<a id="JosephsonCircuits.AndersonState"></a>
```

[`JosephsonCircuits.AndersonState`](api/internals.md#JosephsonCircuits.AndersonState)

```@raw html
<a id="JosephsonCircuits.BathLadder"></a>
```

[`JosephsonCircuits.BathLadder`](api/internals.md#JosephsonCircuits.BathLadder)

```@raw html
<a id="JosephsonCircuits.BlockJacobian"></a>
```

[`JosephsonCircuits.BlockJacobian`](api/internals.md#JosephsonCircuits.BlockJacobian)

```@raw html
<a id="JosephsonCircuits.BlockLU"></a>
```

[`JosephsonCircuits.BlockLU`](api/internals.md#JosephsonCircuits.BlockLU)

```@raw html
<a id="JosephsonCircuits.BlockModulation"></a>
```

[`JosephsonCircuits.BlockModulation`](api/internals.md#JosephsonCircuits.BlockModulation)

```@raw html
<a id="JosephsonCircuits.BlockStructure"></a>
```

[`JosephsonCircuits.BlockStructure`](api/internals.md#JosephsonCircuits.BlockStructure)

```@raw html
<a id="JosephsonCircuits.BoundCircuit"></a>
```

[`JosephsonCircuits.BoundCircuit`](api/internals.md#JosephsonCircuits.BoundCircuit)

```@raw html
<a id="JosephsonCircuits.BranchStampPlan"></a>
```

[`JosephsonCircuits.BranchStampPlan`](api/internals.md#JosephsonCircuits.BranchStampPlan)

```@raw html
<a id="JosephsonCircuits.CallableMatrixProvider"></a>
```

[`JosephsonCircuits.CallableMatrixProvider`](api/internals.md#JosephsonCircuits.CallableMatrixProvider)

```@raw html
<a id="JosephsonCircuits.CanonicalJacobianPlan"></a>
```

[`JosephsonCircuits.CanonicalJacobianPlan`](api/internals.md#JosephsonCircuits.CanonicalJacobianPlan)

```@raw html
<a id="JosephsonCircuits.CanonicalPreconditioner"></a>
```

[`JosephsonCircuits.CanonicalPreconditioner`](api/internals.md#JosephsonCircuits.CanonicalPreconditioner)

```@raw html
<a id="JosephsonCircuits.CanonicalWork"></a>
```

[`JosephsonCircuits.CanonicalWork`](api/internals.md#JosephsonCircuits.CanonicalWork)

```@raw html
<a id="JosephsonCircuits.CircuitGraph"></a>
```

[`JosephsonCircuits.CircuitGraph`](api/internals.md#JosephsonCircuits.CircuitGraph)

```@raw html
<a id="JosephsonCircuits.CircuitMatrices"></a>
```

[`JosephsonCircuits.CircuitMatrices`](api/internals.md#JosephsonCircuits.CircuitMatrices)

```@raw html
<a id="JosephsonCircuits.CircuitMatrixPlan"></a>
```

[`JosephsonCircuits.CircuitMatrixPlan`](api/internals.md#JosephsonCircuits.CircuitMatrixPlan)

```@raw html
<a id="JosephsonCircuits.CircuitMatrixWorkspace"></a>
```

[`JosephsonCircuits.CircuitMatrixWorkspace`](api/internals.md#JosephsonCircuits.CircuitMatrixWorkspace)

```@raw html
<a id="JosephsonCircuits.CircuitTopology"></a>
```

[`JosephsonCircuits.CircuitTopology`](api/internals.md#JosephsonCircuits.CircuitTopology)

```@raw html
<a id="JosephsonCircuits.CircuitValue"></a>
```

[`JosephsonCircuits.CircuitValue`](api/internals.md#JosephsonCircuits.CircuitValue)

```@raw html
<a id="JosephsonCircuits.ClusterBlocks"></a>
```

[`JosephsonCircuits.ClusterBlocks`](api/internals.md#JosephsonCircuits.ClusterBlocks)

```@raw html
<a id="JosephsonCircuits.ClusterProbe"></a>
```

[`JosephsonCircuits.ClusterProbe`](api/internals.md#JosephsonCircuits.ClusterProbe)

```@raw html
<a id="JosephsonCircuits.CompiledPort"></a>
```

[`JosephsonCircuits.CompiledPort`](api/internals.md#JosephsonCircuits.CompiledPort)

```@raw html
<a id="JosephsonCircuits.CompiledScatteringBlock"></a>
```

[`JosephsonCircuits.CompiledScatteringBlock`](api/internals.md#JosephsonCircuits.CompiledScatteringBlock)

```@raw html
<a id="JosephsonCircuits.ComponentPerturbation"></a>
```

[`JosephsonCircuits.ComponentPerturbation`](api/internals.md#JosephsonCircuits.ComponentPerturbation)

```@raw html
<a id="JosephsonCircuits.CompositeLayout"></a>
```

[`JosephsonCircuits.CompositeLayout`](api/internals.md#JosephsonCircuits.CompositeLayout)

```@raw html
<a id="JosephsonCircuits.ConstantMatrixProvider"></a>
```

[`JosephsonCircuits.ConstantMatrixProvider`](api/internals.md#JosephsonCircuits.ConstantMatrixProvider)

```@raw html
<a id="JosephsonCircuits.CoupledLinesBasis"></a>
```

[`JosephsonCircuits.CoupledLinesBasis`](api/internals.md#JosephsonCircuits.CoupledLinesBasis)

```@raw html
<a id="JosephsonCircuits.DCAugmentation"></a>
```

[`JosephsonCircuits.DCAugmentation`](api/internals.md#JosephsonCircuits.DCAugmentation)

```@raw html
<a id="JosephsonCircuits.DCBlockDescriptor"></a>
```

[`JosephsonCircuits.DCBlockDescriptor`](api/internals.md#JosephsonCircuits.DCBlockDescriptor)

```@raw html
<a id="JosephsonCircuits.DCBlockRows"></a>
```

[`JosephsonCircuits.DCBlockRows`](api/internals.md#JosephsonCircuits.DCBlockRows)

```@raw html
<a id="JosephsonCircuits.DCConductancePlan"></a>
```

[`JosephsonCircuits.DCConductancePlan`](api/internals.md#JosephsonCircuits.DCConductancePlan)

```@raw html
<a id="JosephsonCircuits.DCConductanceSolution"></a>
```

[`JosephsonCircuits.DCConductanceSolution`](api/internals.md#JosephsonCircuits.DCConductanceSolution)

```@raw html
<a id="JosephsonCircuits.DCFactorization"></a>
```

[`JosephsonCircuits.DCFactorization`](api/internals.md#JosephsonCircuits.DCFactorization)

```@raw html
<a id="JosephsonCircuits.DCOperatingPoint"></a>
```

[`JosephsonCircuits.DCOperatingPoint`](api/internals.md#JosephsonCircuits.DCOperatingPoint)

```@raw html
<a id="JosephsonCircuits.DCPinning"></a>
```

[`JosephsonCircuits.DCPinning`](api/internals.md#JosephsonCircuits.DCPinning)

```@raw html
<a id="JosephsonCircuits.DCUpdate"></a>
```

[`JosephsonCircuits.DCUpdate`](api/internals.md#JosephsonCircuits.DCUpdate)

```@raw html
<a id="JosephsonCircuits.DeviceBlockNoisePlan"></a>
```

[`JosephsonCircuits.DeviceBlockNoisePlan`](api/internals.md#JosephsonCircuits.DeviceBlockNoisePlan)

```@raw html
<a id="JosephsonCircuits.DeviceNoisePlan"></a>
```

[`JosephsonCircuits.DeviceNoisePlan`](api/internals.md#JosephsonCircuits.DeviceNoisePlan)

```@raw html
<a id="JosephsonCircuits.DeviceProviders"></a>
```

[`JosephsonCircuits.DeviceProviders`](api/internals.md#JosephsonCircuits.DeviceProviders)

```@raw html
<a id="JosephsonCircuits.DeviceScatteringStamps"></a>
```

[`JosephsonCircuits.DeviceScatteringStamps`](api/internals.md#JosephsonCircuits.DeviceScatteringStamps)

```@raw html
<a id="JosephsonCircuits.DeviceSparsePattern"></a>
```

[`JosephsonCircuits.DeviceSparsePattern`](api/internals.md#JosephsonCircuits.DeviceSparsePattern)

```@raw html
<a id="JosephsonCircuits.DeviceSweep"></a>
```

[`JosephsonCircuits.DeviceSweep`](api/internals.md#JosephsonCircuits.DeviceSweep)

```@raw html
<a id="JosephsonCircuits.DeviceValuedSparseMatrix"></a>
```

[`JosephsonCircuits.DeviceValuedSparseMatrix`](api/internals.md#JosephsonCircuits.DeviceValuedSparseMatrix)

```@raw html
<a id="JosephsonCircuits.ErasedFunction"></a>
```

[`JosephsonCircuits.ErasedFunction`](api/internals.md#JosephsonCircuits.ErasedFunction)

```@raw html
<a id="JosephsonCircuits.ErasedPreconditioner"></a>
```

[`JosephsonCircuits.ErasedPreconditioner`](api/internals.md#JosephsonCircuits.ErasedPreconditioner)

```@raw html
<a id="JosephsonCircuits.FactorizationCache"></a>
```

[`JosephsonCircuits.FactorizationCache`](api/internals.md#JosephsonCircuits.FactorizationCache)

```@raw html
<a id="JosephsonCircuits.FillOrdering"></a>
```

[`JosephsonCircuits.FillOrdering`](api/internals.md#JosephsonCircuits.FillOrdering)

```@raw html
<a id="JosephsonCircuits.FloquetPreconditioner"></a>
```

[`JosephsonCircuits.FloquetPreconditioner`](api/internals.md#JosephsonCircuits.FloquetPreconditioner)

```@raw html
<a id="JosephsonCircuits.FloquetState"></a>
```

[`JosephsonCircuits.FloquetState`](api/internals.md#JosephsonCircuits.FloquetState)

```@raw html
<a id="JosephsonCircuits.FourierIndices"></a>
```

[`JosephsonCircuits.FourierIndices`](api/internals.md#JosephsonCircuits.FourierIndices)

```@raw html
<a id="JosephsonCircuits.Frequencies"></a>
```

[`JosephsonCircuits.Frequencies`](api/internals.md#JosephsonCircuits.Frequencies)

```@raw html
<a id="JosephsonCircuits.FrequencySweepPlan"></a>
```

[`JosephsonCircuits.FrequencySweepPlan`](api/internals.md#JosephsonCircuits.FrequencySweepPlan)

```@raw html
<a id="JosephsonCircuits.FunctionOperator"></a>
```

[`JosephsonCircuits.FunctionOperator`](api/internals.md#JosephsonCircuits.FunctionOperator)

```@raw html
<a id="JosephsonCircuits.GMRESWorkspace"></a>
```

[`JosephsonCircuits.GMRESWorkspace`](api/internals.md#JosephsonCircuits.GMRESWorkspace)

```@raw html
<a id="JosephsonCircuits.GMRESWorkspace-Union{Tuple{T}, Tuple{AbstractVector{T}, Integer}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.GMRESWorkspace`](api/internals.md#JosephsonCircuits.GMRESWorkspace-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractVector%7BT%7D%2C%20Integer%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.GroundType"></a>
```

[`JosephsonCircuits.GroundType`](api/internals.md#JosephsonCircuits.GroundType)

```@raw html
<a id="JosephsonCircuits.HBLinearizedSystem"></a>
```

[`JosephsonCircuits.HBLinearizedSystem`](api/internals.md#JosephsonCircuits.HBLinearizedSystem)

```@raw html
<a id="JosephsonCircuits.HBLinearizedSystem-Tuple{Matrix, SparseVector, SparseMatrixCSC, Integer, Integer, Array, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, Any, Any, Any, Bool, SparseMatrixCSC, Any, Integer}"></a>
```

[`JosephsonCircuits.HBLinearizedSystem`](api/internals.md#JosephsonCircuits.HBLinearizedSystem-Tuple%7BMatrix%2C%20SparseVector%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Array%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Any%2C%20Any%2C%20Any%2C%20Bool%2C%20SparseMatrixCSC%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.HBOperatingPoint"></a>
```

[`JosephsonCircuits.HBOperatingPoint`](api/internals.md#JosephsonCircuits.HBOperatingPoint)

```@raw html
<a id="JosephsonCircuits.HBSystem"></a>
```

[`JosephsonCircuits.HBSystem`](api/internals.md#JosephsonCircuits.HBSystem)

```@raw html
<a id="JosephsonCircuits.HybridWorkspace"></a>
```

[`JosephsonCircuits.HybridWorkspace`](api/internals.md#JosephsonCircuits.HybridWorkspace)

```@raw html
<a id="JosephsonCircuits.InverseInductancePlan"></a>
```

[`JosephsonCircuits.InverseInductancePlan`](api/internals.md#JosephsonCircuits.InverseInductancePlan)

```@raw html
<a id="JosephsonCircuits.IterationInfo"></a>
```

[`JosephsonCircuits.IterationInfo`](api/internals.md#JosephsonCircuits.IterationInfo)

```@raw html
<a id="JosephsonCircuits.JunctionRelations"></a>
```

[`JosephsonCircuits.JunctionRelations`](api/internals.md#JosephsonCircuits.JunctionRelations)

```@raw html
<a id="JosephsonCircuits.JunctionStructure"></a>
```

[`JosephsonCircuits.JunctionStructure`](api/internals.md#JosephsonCircuits.JunctionStructure)

```@raw html
<a id="JosephsonCircuits.KrylovSolveInfo"></a>
```

[`JosephsonCircuits.KrylovSolveInfo`](api/internals.md#JosephsonCircuits.KrylovSolveInfo)

```@raw html
<a id="JosephsonCircuits.KrylovVectors"></a>
```

[`JosephsonCircuits.KrylovVectors`](api/internals.md#JosephsonCircuits.KrylovVectors)

```@raw html
<a id="JosephsonCircuits.LegacyTermination"></a>
```

[`JosephsonCircuits.LegacyTermination`](api/internals.md#JosephsonCircuits.LegacyTermination)

```@raw html
<a id="JosephsonCircuits.LinearTermGather"></a>
```

[`JosephsonCircuits.LinearTermGather`](api/internals.md#JosephsonCircuits.LinearTermGather)

```@raw html
<a id="JosephsonCircuits.LinearizedArrays"></a>
```

[`JosephsonCircuits.LinearizedArrays`](api/internals.md#JosephsonCircuits.LinearizedArrays)

```@raw html
<a id="JosephsonCircuits.LinearizedWorkspace"></a>
```

[`JosephsonCircuits.LinearizedWorkspace`](api/internals.md#JosephsonCircuits.LinearizedWorkspace)

```@raw html
<a id="JosephsonCircuits.LinearizedWorkspace-Tuple{JosephsonCircuits.LinearizedArrays, Any, Any, Integer, Integer, Integer, Integer, Any}"></a>
```

[`JosephsonCircuits.LinearizedWorkspace`](api/internals.md#JosephsonCircuits.LinearizedWorkspace-Tuple%7BJosephsonCircuits.LinearizedArrays%2C%20Any%2C%20Any%2C%20Integer%2C%20Integer%2C%20Integer%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.ModeCouplingPreconditioner"></a>
```

[`JosephsonCircuits.ModeCouplingPreconditioner`](api/internals.md#JosephsonCircuits.ModeCouplingPreconditioner)

```@raw html
<a id="JosephsonCircuits.ModeCouplingPreconditioner-Tuple{Any, Matrix, Matrix, SparseVector, Any, SparseMatrixCSC, Integer, Integer, Integer, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, JosephsonCircuits.ModeLayout}"></a>
```

[`JosephsonCircuits.ModeCouplingPreconditioner`](api/internals.md#JosephsonCircuits.ModeCouplingPreconditioner-Tuple%7BAny%2C%20Matrix%2C%20Matrix%2C%20SparseVector%2C%20Any%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Integer%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20JosephsonCircuits.ModeLayout%7D)

```@raw html
<a id="JosephsonCircuits.ModeDifferences"></a>
```

[`JosephsonCircuits.ModeDifferences`](api/internals.md#JosephsonCircuits.ModeDifferences)

```@raw html
<a id="JosephsonCircuits.ModeIndices"></a>
```

[`JosephsonCircuits.ModeIndices`](api/internals.md#JosephsonCircuits.ModeIndices)

```@raw html
<a id="JosephsonCircuits.ModeLayout"></a>
```

[`JosephsonCircuits.ModeLayout`](api/internals.md#JosephsonCircuits.ModeLayout)

```@raw html
<a id="JosephsonCircuits.ModulatedRationalProvider"></a>
```

[`JosephsonCircuits.ModulatedRationalProvider`](api/internals.md#JosephsonCircuits.ModulatedRationalProvider)

```@raw html
<a id="JosephsonCircuits.MutualStampPlan"></a>
```

[`JosephsonCircuits.MutualStampPlan`](api/internals.md#JosephsonCircuits.MutualStampPlan)

```@raw html
<a id="JosephsonCircuits.NewtonTrace"></a>
```

[`JosephsonCircuits.NewtonTrace`](api/internals.md#JosephsonCircuits.NewtonTrace)

```@raw html
<a id="JosephsonCircuits.NoPortTermination"></a>
```

[`JosephsonCircuits.NoPortTermination`](api/internals.md#JosephsonCircuits.NoPortTermination)

```@raw html
<a id="JosephsonCircuits.NodalStampPlan"></a>
```

[`JosephsonCircuits.NodalStampPlan`](api/internals.md#JosephsonCircuits.NodalStampPlan)

```@raw html
<a id="JosephsonCircuits.NoiseReduction"></a>
```

[`JosephsonCircuits.NoiseReduction`](api/internals.md#JosephsonCircuits.NoiseReduction)

```@raw html
<a id="JosephsonCircuits.NonlinearTermPlan"></a>
```

[`JosephsonCircuits.NonlinearTermPlan`](api/internals.md#JosephsonCircuits.NonlinearTermPlan)

```@raw html
<a id="JosephsonCircuits.NonlinearTermTransposePlan"></a>
```

[`JosephsonCircuits.NonlinearTermTransposePlan`](api/internals.md#JosephsonCircuits.NonlinearTermTransposePlan)

```@raw html
<a id="JosephsonCircuits.PaddedLinearTerm"></a>
```

[`JosephsonCircuits.PaddedLinearTerm`](api/internals.md#JosephsonCircuits.PaddedLinearTerm)

```@raw html
<a id="JosephsonCircuits.PaddedLinearTerm-NTuple{13, Any}"></a>
```

[`JosephsonCircuits.PaddedLinearTerm`](api/internals.md#JosephsonCircuits.PaddedLinearTerm-NTuple%7B13%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.PairLadder"></a>
```

[`JosephsonCircuits.PairLadder`](api/internals.md#JosephsonCircuits.PairLadder)

```@raw html
<a id="JosephsonCircuits.ParsedLevel"></a>
```

[`JosephsonCircuits.ParsedLevel`](api/internals.md#JosephsonCircuits.ParsedLevel)

```@raw html
<a id="JosephsonCircuits.PassiveNetwork"></a>
```

[`JosephsonCircuits.PassiveNetwork`](api/internals.md#JosephsonCircuits.PassiveNetwork)

```@raw html
<a id="JosephsonCircuits.PerturbationEntries"></a>
```

[`JosephsonCircuits.PerturbationEntries`](api/internals.md#JosephsonCircuits.PerturbationEntries)

```@raw html
<a id="JosephsonCircuits.PiecewiseTabulatedProvider"></a>
```

[`JosephsonCircuits.PiecewiseTabulatedProvider`](api/internals.md#JosephsonCircuits.PiecewiseTabulatedProvider)

```@raw html
<a id="JosephsonCircuits.PolynomialCPRDerivative"></a>
```

[`JosephsonCircuits.PolynomialCPRDerivative`](api/internals.md#JosephsonCircuits.PolynomialCPRDerivative)

```@raw html
<a id="JosephsonCircuits.PortDiagonal"></a>
```

[`JosephsonCircuits.PortDiagonal`](api/internals.md#JosephsonCircuits.PortDiagonal)

```@raw html
<a id="JosephsonCircuits.PreconditionerPlan"></a>
```

[`JosephsonCircuits.PreconditionerPlan`](api/internals.md#JosephsonCircuits.PreconditionerPlan)

```@raw html
<a id="JosephsonCircuits.RationalScatteringProvider"></a>
```

[`JosephsonCircuits.RationalScatteringProvider`](api/internals.md#JosephsonCircuits.RationalScatteringProvider)

```@raw html
<a id="JosephsonCircuits.ReverseSensitivity"></a>
```

[`JosephsonCircuits.ReverseSensitivity`](api/internals.md#JosephsonCircuits.ReverseSensitivity)

```@raw html
<a id="JosephsonCircuits.ReverseSensitivityBuffers"></a>
```

[`JosephsonCircuits.ReverseSensitivityBuffers`](api/internals.md#JosephsonCircuits.ReverseSensitivityBuffers)

```@raw html
<a id="JosephsonCircuits.RotatedMatrixProvider"></a>
```

[`JosephsonCircuits.RotatedMatrixProvider`](api/internals.md#JosephsonCircuits.RotatedMatrixProvider)

```@raw html
<a id="JosephsonCircuits.ScatteringNoisePlan"></a>
```

[`JosephsonCircuits.ScatteringNoisePlan`](api/internals.md#JosephsonCircuits.ScatteringNoisePlan)

```@raw html
<a id="JosephsonCircuits.ScatteringNoiseWorkspace"></a>
```

[`JosephsonCircuits.ScatteringNoiseWorkspace`](api/internals.md#JosephsonCircuits.ScatteringNoiseWorkspace)

```@raw html
<a id="JosephsonCircuits.ScatteringStampSystem"></a>
```

[`JosephsonCircuits.ScatteringStampSystem`](api/internals.md#JosephsonCircuits.ScatteringStampSystem)

```@raw html
<a id="JosephsonCircuits.ScatteringWorkspace"></a>
```

[`JosephsonCircuits.ScatteringWorkspace`](api/internals.md#JosephsonCircuits.ScatteringWorkspace)

```@raw html
<a id="JosephsonCircuits.SensitivityStamp"></a>
```

[`JosephsonCircuits.SensitivityStamp`](api/internals.md#JosephsonCircuits.SensitivityStamp)

```@raw html
<a id="JosephsonCircuits.SizedPreconditioner"></a>
```

[`JosephsonCircuits.SizedPreconditioner`](api/internals.md#JosephsonCircuits.SizedPreconditioner)

```@raw html
<a id="JosephsonCircuits.SolverInfo"></a>
```

[`JosephsonCircuits.SolverInfo`](api/internals.md#JosephsonCircuits.SolverInfo)

```@raw html
<a id="JosephsonCircuits.SourceTuple"></a>
```

[`JosephsonCircuits.SourceTuple`](api/internals.md#JosephsonCircuits.SourceTuple)

```@raw html
<a id="JosephsonCircuits.SparseBlockFactorization"></a>
```

[`JosephsonCircuits.SparseBlockFactorization`](api/internals.md#JosephsonCircuits.SparseBlockFactorization)

```@raw html
<a id="JosephsonCircuits.SpiceRaw"></a>
```

[`JosephsonCircuits.SpiceRaw`](api/internals.md#JosephsonCircuits.SpiceRaw)

```@raw html
<a id="JosephsonCircuits.SpiceRawHeader"></a>
```

[`JosephsonCircuits.SpiceRawHeader`](api/internals.md#JosephsonCircuits.SpiceRawHeader)

```@raw html
<a id="JosephsonCircuits.StageCorrection"></a>
```

[`JosephsonCircuits.StageCorrection`](api/internals.md#JosephsonCircuits.StageCorrection)

```@raw html
<a id="JosephsonCircuits.StagePlan"></a>
```

[`JosephsonCircuits.StagePlan`](api/internals.md#JosephsonCircuits.StagePlan)

```@raw html
<a id="JosephsonCircuits.StagedStageInfo"></a>
```

[`JosephsonCircuits.StagedStageInfo`](api/internals.md#JosephsonCircuits.StagedStageInfo)

```@raw html
<a id="JosephsonCircuits.StampedScatteringBlock"></a>
```

[`JosephsonCircuits.StampedScatteringBlock`](api/internals.md#JosephsonCircuits.StampedScatteringBlock)

```@raw html
<a id="JosephsonCircuits.StructureComplexJacobianPlan"></a>
```

[`JosephsonCircuits.StructureComplexJacobianPlan`](api/internals.md#JosephsonCircuits.StructureComplexJacobianPlan)

```@raw html
<a id="JosephsonCircuits.StructureComplexJosephsonPlan"></a>
```

[`JosephsonCircuits.StructureComplexJosephsonPlan`](api/internals.md#JosephsonCircuits.StructureComplexJosephsonPlan)

```@raw html
<a id="JosephsonCircuits.StructureRealJacobianPlan"></a>
```

[`JosephsonCircuits.StructureRealJacobianPlan`](api/internals.md#JosephsonCircuits.StructureRealJacobianPlan)

```@raw html
<a id="JosephsonCircuits.SweepEquilibration"></a>
```

[`JosephsonCircuits.SweepEquilibration`](api/internals.md#JosephsonCircuits.SweepEquilibration)

```@raw html
<a id="JosephsonCircuits.TabulatedMatrixProvider"></a>
```

[`JosephsonCircuits.TabulatedMatrixProvider`](api/internals.md#JosephsonCircuits.TabulatedMatrixProvider)

```@raw html
<a id="JosephsonCircuits.TransientBlock"></a>
```

[`JosephsonCircuits.TransientBlock`](api/internals.md#JosephsonCircuits.TransientBlock)

```@raw html
<a id="JosephsonCircuits.TransientIQPlan"></a>
```

[`JosephsonCircuits.TransientIQPlan`](api/internals.md#JosephsonCircuits.TransientIQPlan)

```@raw html
<a id="JosephsonCircuits.TransientLine"></a>
```

[`JosephsonCircuits.TransientLine`](api/internals.md#JosephsonCircuits.TransientLine)

```@raw html
<a id="JosephsonCircuits.TransientNoiseBath"></a>
```

[`JosephsonCircuits.TransientNoiseBath`](api/internals.md#JosephsonCircuits.TransientNoiseBath)

```@raw html
<a id="JosephsonCircuits.TransientNoiseBaths"></a>
```

[`JosephsonCircuits.TransientNoiseBaths`](api/internals.md#JosephsonCircuits.TransientNoiseBaths)

```@raw html
<a id="JosephsonCircuits.TransientProblem"></a>
```

[`JosephsonCircuits.TransientProblem`](api/internals.md#JosephsonCircuits.TransientProblem)

```@raw html
<a id="JosephsonCircuits.TransientQuantumPlan"></a>
```

[`JosephsonCircuits.TransientQuantumPlan`](api/internals.md#JosephsonCircuits.TransientQuantumPlan)

```@raw html
<a id="JosephsonCircuits.TransientSystem"></a>
```

[`JosephsonCircuits.TransientSystem`](api/internals.md#JosephsonCircuits.TransientSystem)

```@raw html
<a id="JosephsonCircuits.TransmissionLineProvider"></a>
```

[`JosephsonCircuits.TransmissionLineProvider`](api/internals.md#JosephsonCircuits.TransmissionLineProvider)

```@raw html
<a id="JosephsonCircuits.TransportRows"></a>
```

[`JosephsonCircuits.TransportRows`](api/internals.md#JosephsonCircuits.TransportRows)

```@raw html
<a id="JosephsonCircuits.ValueMaps"></a>
```

[`JosephsonCircuits.ValueMaps`](api/internals.md#JosephsonCircuits.ValueMaps)

```@raw html
<a id="JosephsonCircuits.WorkerBlockSensitivity"></a>
```

[`JosephsonCircuits.WorkerBlockSensitivity`](api/internals.md#JosephsonCircuits.WorkerBlockSensitivity)

```@raw html
<a id="JosephsonCircuits.ABCD_PiY!-Tuple{AbstractMatrix, Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_PiY!`](api/internals.md#JosephsonCircuits.ABCD_PiY!-Tuple%7BAbstractMatrix%2C%20Number%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_PiY-Tuple{Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_PiY`](api/internals.md#JosephsonCircuits.ABCD_PiY-Tuple%7BNumber%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_TZ!-Tuple{AbstractMatrix, Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_TZ!`](api/internals.md#JosephsonCircuits.ABCD_TZ!-Tuple%7BAbstractMatrix%2C%20Number%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_TZ-Tuple{Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_TZ`](api/internals.md#JosephsonCircuits.ABCD_TZ-Tuple%7BNumber%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_Pi!-Tuple{Any, Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_Pi!`](api/internals.md#JosephsonCircuits.ABCD_attenuator_Pi!-Tuple%7BAny%2C%20Number%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_Pi!-Tuple{Any, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_Pi!`](api/internals.md#JosephsonCircuits.ABCD_attenuator_Pi!-Tuple%7BAny%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_Pi-Tuple{Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_Pi`](api/internals.md#JosephsonCircuits.ABCD_attenuator_Pi-Tuple%7BNumber%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_Pi-Tuple{Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_Pi`](api/internals.md#JosephsonCircuits.ABCD_attenuator_Pi-Tuple%7BNumber%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_T!-Tuple{Any, Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_T!`](api/internals.md#JosephsonCircuits.ABCD_attenuator_T!-Tuple%7BAny%2C%20Number%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_T!-Tuple{Any, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_T!`](api/internals.md#JosephsonCircuits.ABCD_attenuator_T!-Tuple%7BAny%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_T-Tuple{Number, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_T`](api/internals.md#JosephsonCircuits.ABCD_attenuator_T-Tuple%7BNumber%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_attenuator_T-Tuple{Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_attenuator_T`](api/internals.md#JosephsonCircuits.ABCD_attenuator_T-Tuple%7BNumber%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_coupled_tline!-Tuple{AbstractMatrix, Vararg{Number, 4}}"></a>
```

[`JosephsonCircuits.ABCD_coupled_tline!`](api/internals.md#JosephsonCircuits.ABCD_coupled_tline!-Tuple%7BAbstractMatrix%2C%20Vararg%7BNumber%2C%204%7D%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_coupled_tline-NTuple{4, Number}"></a>
```

[`JosephsonCircuits.ABCD_coupled_tline`](api/internals.md#JosephsonCircuits.ABCD_coupled_tline-NTuple%7B4%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_seriesZ!-Tuple{AbstractMatrix, Number}"></a>
```

[`JosephsonCircuits.ABCD_seriesZ!`](api/internals.md#JosephsonCircuits.ABCD_seriesZ!-Tuple%7BAbstractMatrix%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_seriesZ-Tuple{Number}"></a>
```

[`JosephsonCircuits.ABCD_seriesZ`](api/internals.md#JosephsonCircuits.ABCD_seriesZ-Tuple%7BNumber%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_shuntY!-Tuple{AbstractMatrix, Number}"></a>
```

[`JosephsonCircuits.ABCD_shuntY!`](api/internals.md#JosephsonCircuits.ABCD_shuntY!-Tuple%7BAbstractMatrix%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_shuntY-Tuple{Number}"></a>
```

[`JosephsonCircuits.ABCD_shuntY`](api/internals.md#JosephsonCircuits.ABCD_shuntY-Tuple%7BNumber%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_tline!-Tuple{AbstractMatrix, Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_tline!`](api/internals.md#JosephsonCircuits.ABCD_tline!-Tuple%7BAbstractMatrix%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCD_tline-Tuple{Number, Number}"></a>
```

[`JosephsonCircuits.ABCD_tline`](api/internals.md#JosephsonCircuits.ABCD_tline-Tuple%7BNumber%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ABCDtoS"></a>
```

[`JosephsonCircuits.ABCDtoS`](api/internals.md#JosephsonCircuits.ABCDtoS)

```@raw html
<a id="JosephsonCircuits.A_B_to_symplectic_pair-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.A_B_to_symplectic_pair`](api/internals.md#JosephsonCircuits.A_B_to_symplectic_pair-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.A_coupled_tlines-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.A_coupled_tlines`](api/internals.md#JosephsonCircuits.A_coupled_tlines-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.AtoB"></a>
```

[`JosephsonCircuits.AtoB`](api/internals.md#JosephsonCircuits.AtoB)

```@raw html
<a id="JosephsonCircuits.AtoB!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.AtoB!`](api/internals.md#JosephsonCircuits.AtoB!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.AtoS"></a>
```

[`JosephsonCircuits.AtoS`](api/internals.md#JosephsonCircuits.AtoS)

```@raw html
<a id="JosephsonCircuits.AtoS!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any, Any}"></a>
```

[`JosephsonCircuits.AtoS!`](api/internals.md#JosephsonCircuits.AtoS!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.AtoY"></a>
```

[`JosephsonCircuits.AtoY`](api/internals.md#JosephsonCircuits.AtoY)

```@raw html
<a id="JosephsonCircuits.AtoY!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.AtoY!`](api/internals.md#JosephsonCircuits.AtoY!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.AtoZ"></a>
```

[`JosephsonCircuits.AtoZ`](api/internals.md#JosephsonCircuits.AtoZ)

```@raw html
<a id="JosephsonCircuits.AtoZ!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.AtoZ!`](api/internals.md#JosephsonCircuits.AtoZ!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.B_from_X_Y_quadrature-Tuple{AbstractMatrix, AbstractMatrix{&lt;:Real}, AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.B_from_X_Y_quadrature`](api/internals.md#JosephsonCircuits.B_from_X_Y_quadrature-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7B%3C%3AReal%7D%2C%20AbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.B_from_X_Y_quadrature_block-Tuple{AbstractMatrix{&lt;:Real}, AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.B_from_X_Y_quadrature_block`](api/internals.md#JosephsonCircuits.B_from_X_Y_quadrature_block-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%2C%20AbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.B_from_X_Y_quadrature_pair-Tuple{AbstractMatrix{&lt;:Real}, AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.B_from_X_Y_quadrature_pair`](api/internals.md#JosephsonCircuits.B_from_X_Y_quadrature_pair-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%2C%20AbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.BtoA"></a>
```

[`JosephsonCircuits.BtoA`](api/internals.md#JosephsonCircuits.BtoA)

```@raw html
<a id="JosephsonCircuits.BtoA!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.BtoA!`](api/internals.md#JosephsonCircuits.BtoA!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.BtoS"></a>
```

[`JosephsonCircuits.BtoS`](api/internals.md#JosephsonCircuits.BtoS)

```@raw html
<a id="JosephsonCircuits.BtoS!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any, Any}"></a>
```

[`JosephsonCircuits.BtoS!`](api/internals.md#JosephsonCircuits.BtoS!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.BtoY"></a>
```

[`JosephsonCircuits.BtoY`](api/internals.md#JosephsonCircuits.BtoY)

```@raw html
<a id="JosephsonCircuits.BtoY!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.BtoY!`](api/internals.md#JosephsonCircuits.BtoY!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.BtoZ"></a>
```

[`JosephsonCircuits.BtoZ`](api/internals.md#JosephsonCircuits.BtoZ)

```@raw html
<a id="JosephsonCircuits.BtoZ!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.BtoZ!`](api/internals.md#JosephsonCircuits.BtoZ!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.CMtokeyed-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.CMtokeyed`](api/internals.md#JosephsonCircuits.CMtokeyed-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Cnoisetokeyed-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.Cnoisetokeyed`](api/internals.md#JosephsonCircuits.Cnoisetokeyed-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.R_block_to_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.R_block_to_pair`](api/internals.md#JosephsonCircuits.R_block_to_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.R_ladder_to_quadrature_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.R_ladder_to_quadrature_block`](api/internals.md#JosephsonCircuits.R_ladder_to_quadrature_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.R_ladder_to_quadrature_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.R_ladder_to_quadrature_pair`](api/internals.md#JosephsonCircuits.R_ladder_to_quadrature_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.R_pair_to_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.R_pair_to_block`](api/internals.md#JosephsonCircuits.R_pair_to_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.R_quadrature_to_ladder_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.R_quadrature_to_ladder_block`](api/internals.md#JosephsonCircuits.R_quadrature_to_ladder_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.R_quadrature_to_ladder_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.R_quadrature_to_ladder_pair`](api/internals.md#JosephsonCircuits.R_quadrature_to_ladder_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.S_circulator_clockwise!-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.S_circulator_clockwise!`](api/internals.md#JosephsonCircuits.S_circulator_clockwise!-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.S_circulator_clockwise-Tuple{}"></a>
```

[`JosephsonCircuits.S_circulator_clockwise`](api/internals.md#JosephsonCircuits.S_circulator_clockwise-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.S_circulator_counterclockwise!-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.S_circulator_counterclockwise!`](api/internals.md#JosephsonCircuits.S_circulator_counterclockwise!-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.S_circulator_counterclockwise-Tuple{}"></a>
```

[`JosephsonCircuits.S_circulator_counterclockwise`](api/internals.md#JosephsonCircuits.S_circulator_counterclockwise-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.S_directional_coupler!-Tuple{AbstractMatrix, Vararg{Number, 4}}"></a>
```

[`JosephsonCircuits.S_directional_coupler!`](api/internals.md#JosephsonCircuits.S_directional_coupler!-Tuple%7BAbstractMatrix%2C%20Vararg%7BNumber%2C%204%7D%7D)

```@raw html
<a id="JosephsonCircuits.S_directional_coupler-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.S_directional_coupler`](api/internals.md#JosephsonCircuits.S_directional_coupler-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.S_directional_coupler_antisymmetric!-Tuple{Any, Number}"></a>
```

[`JosephsonCircuits.S_directional_coupler_antisymmetric!`](api/internals.md#JosephsonCircuits.S_directional_coupler_antisymmetric!-Tuple%7BAny%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.S_directional_coupler_antisymmetric-Tuple{Number}"></a>
```

[`JosephsonCircuits.S_directional_coupler_antisymmetric`](api/internals.md#JosephsonCircuits.S_directional_coupler_antisymmetric-Tuple%7BNumber%7D)

```@raw html
<a id="JosephsonCircuits.S_directional_coupler_symmetric!-Tuple{Any, Number}"></a>
```

[`JosephsonCircuits.S_directional_coupler_symmetric!`](api/internals.md#JosephsonCircuits.S_directional_coupler_symmetric!-Tuple%7BAny%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.S_directional_coupler_symmetric-Tuple{Number}"></a>
```

[`JosephsonCircuits.S_directional_coupler_symmetric`](api/internals.md#JosephsonCircuits.S_directional_coupler_symmetric-Tuple%7BNumber%7D)

```@raw html
<a id="JosephsonCircuits.S_hybrid_coupler_antisymmetric!-Tuple{Any}"></a>
```

[`JosephsonCircuits.S_hybrid_coupler_antisymmetric!`](api/internals.md#JosephsonCircuits.S_hybrid_coupler_antisymmetric!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.S_hybrid_coupler_antisymmetric-Tuple{}"></a>
```

[`JosephsonCircuits.S_hybrid_coupler_antisymmetric`](api/internals.md#JosephsonCircuits.S_hybrid_coupler_antisymmetric-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.S_hybrid_coupler_symmetric!-Tuple{Any}"></a>
```

[`JosephsonCircuits.S_hybrid_coupler_symmetric!`](api/internals.md#JosephsonCircuits.S_hybrid_coupler_symmetric!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.S_hybrid_coupler_symmetric-Tuple{}"></a>
```

[`JosephsonCircuits.S_hybrid_coupler_symmetric`](api/internals.md#JosephsonCircuits.S_hybrid_coupler_symmetric-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.S_match!-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.S_match!`](api/internals.md#JosephsonCircuits.S_match!-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.S_open!-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.S_open!`](api/internals.md#JosephsonCircuits.S_open!-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.S_short!-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.S_short!`](api/internals.md#JosephsonCircuits.S_short!-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.S_splitter!-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.S_splitter!`](api/internals.md#JosephsonCircuits.S_splitter!-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.Snoisetokeyed-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.Snoisetokeyed`](api/internals.md#JosephsonCircuits.Snoisetokeyed-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Ssensitivitytokeyed-NTuple{7, Any}"></a>
```

[`JosephsonCircuits.Ssensitivitytokeyed`](api/internals.md#JosephsonCircuits.Ssensitivitytokeyed-NTuple%7B7%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.StoA"></a>
```

[`JosephsonCircuits.StoA`](api/internals.md#JosephsonCircuits.StoA)

```@raw html
<a id="JosephsonCircuits.StoA!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any, Any}"></a>
```

[`JosephsonCircuits.StoA!`](api/internals.md#JosephsonCircuits.StoA!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.StoABCD"></a>
```

[`JosephsonCircuits.StoABCD`](api/internals.md#JosephsonCircuits.StoABCD)

```@raw html
<a id="JosephsonCircuits.StoB"></a>
```

[`JosephsonCircuits.StoB`](api/internals.md#JosephsonCircuits.StoB)

```@raw html
<a id="JosephsonCircuits.StoB!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any, Any}"></a>
```

[`JosephsonCircuits.StoB!`](api/internals.md#JosephsonCircuits.StoB!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.StoT"></a>
```

[`JosephsonCircuits.StoT`](api/internals.md#JosephsonCircuits.StoT)

```@raw html
<a id="JosephsonCircuits.StoT!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.StoT!`](api/internals.md#JosephsonCircuits.StoT!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.StoY"></a>
```

[`JosephsonCircuits.StoY`](api/internals.md#JosephsonCircuits.StoY)

```@raw html
<a id="JosephsonCircuits.StoY!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.StoY!`](api/internals.md#JosephsonCircuits.StoY!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.StoZ"></a>
```

[`JosephsonCircuits.StoZ`](api/internals.md#JosephsonCircuits.StoZ)

```@raw html
<a id="JosephsonCircuits.StoZ!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.StoZ!`](api/internals.md#JosephsonCircuits.StoZ!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Stokeyed-NTuple{5, Any}"></a>
```

[`JosephsonCircuits.Stokeyed`](api/internals.md#JosephsonCircuits.Stokeyed-NTuple%7B5%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Stokeyed-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.Stokeyed`](api/internals.md#JosephsonCircuits.Stokeyed-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.TtoS"></a>
```

[`JosephsonCircuits.TtoS`](api/internals.md#JosephsonCircuits.TtoS)

```@raw html
<a id="JosephsonCircuits.TtoS!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.TtoS!`](api/internals.md#JosephsonCircuits.TtoS!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.X_Y_to_bogoliubov_block-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.X_Y_to_bogoliubov_block`](api/internals.md#JosephsonCircuits.X_Y_to_bogoliubov_block-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.X_Y_to_bogoliubov_pair-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.X_Y_to_bogoliubov_pair`](api/internals.md#JosephsonCircuits.X_Y_to_bogoliubov_pair-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.X_Y_to_symplectic_block-Tuple{AbstractMatrix{&lt;:Real}, AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.X_Y_to_symplectic_block`](api/internals.md#JosephsonCircuits.X_Y_to_symplectic_block-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%2C%20AbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.X_Y_to_symplectic_pair-Tuple{AbstractMatrix{&lt;:Real}, AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.X_Y_to_symplectic_pair`](api/internals.md#JosephsonCircuits.X_Y_to_symplectic_pair-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%2C%20AbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.Y_C-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.Y_C`](api/internals.md#JosephsonCircuits.Y_C-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Y_PiY!-Tuple{AbstractMatrix, Number, Number, Number}"></a>
```

[`JosephsonCircuits.Y_PiY!`](api/internals.md#JosephsonCircuits.Y_PiY!-Tuple%7BAbstractMatrix%2C%20Number%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Y_PiY-Tuple{Number, Number, Number}"></a>
```

[`JosephsonCircuits.Y_PiY`](api/internals.md#JosephsonCircuits.Y_PiY-Tuple%7BNumber%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Y_invL-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.Y_invL`](api/internals.md#JosephsonCircuits.Y_invL-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Y_seriesY!-Tuple{AbstractMatrix, Number}"></a>
```

[`JosephsonCircuits.Y_seriesY!`](api/internals.md#JosephsonCircuits.Y_seriesY!-Tuple%7BAbstractMatrix%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Y_seriesY-Tuple{Number}"></a>
```

[`JosephsonCircuits.Y_seriesY`](api/internals.md#JosephsonCircuits.Y_seriesY-Tuple%7BNumber%7D)

```@raw html
<a id="JosephsonCircuits.Ymin_from_X-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.Ymin_from_X`](api/internals.md#JosephsonCircuits.Ymin_from_X-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Ymin_from_X_quadrature_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.Ymin_from_X_quadrature_block`](api/internals.md#JosephsonCircuits.Ymin_from_X_quadrature_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.Ymin_from_X_quadrature_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.Ymin_from_X_quadrature_pair`](api/internals.md#JosephsonCircuits.Ymin_from_X_quadrature_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.YtoA"></a>
```

[`JosephsonCircuits.YtoA`](api/internals.md#JosephsonCircuits.YtoA)

```@raw html
<a id="JosephsonCircuits.YtoA!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.YtoA!`](api/internals.md#JosephsonCircuits.YtoA!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.YtoB"></a>
```

[`JosephsonCircuits.YtoB`](api/internals.md#JosephsonCircuits.YtoB)

```@raw html
<a id="JosephsonCircuits.YtoB!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.YtoB!`](api/internals.md#JosephsonCircuits.YtoB!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.YtoS"></a>
```

[`JosephsonCircuits.YtoS`](api/internals.md#JosephsonCircuits.YtoS)

```@raw html
<a id="JosephsonCircuits.YtoS!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.YtoS!`](api/internals.md#JosephsonCircuits.YtoS!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.ZC_basis_coupled_tlines-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.ZC_basis_coupled_tlines`](api/internals.md#JosephsonCircuits.ZC_basis_coupled_tlines-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Z_L-Tuple{AbstractMatrix, Number}"></a>
```

[`JosephsonCircuits.Z_L`](api/internals.md#JosephsonCircuits.Z_L-Tuple%7BAbstractMatrix%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Z_TZ!-Tuple{AbstractMatrix, Number, Number, Number}"></a>
```

[`JosephsonCircuits.Z_TZ!`](api/internals.md#JosephsonCircuits.Z_TZ!-Tuple%7BAbstractMatrix%2C%20Number%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Z_TZ-Tuple{Number, Number, Number}"></a>
```

[`JosephsonCircuits.Z_TZ`](api/internals.md#JosephsonCircuits.Z_TZ-Tuple%7BNumber%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Z_canonical_coupled_line_circuits-Tuple{Int64, Vararg{Any, 4}}"></a>
```

[`JosephsonCircuits.Z_canonical_coupled_line_circuits`](api/internals.md#JosephsonCircuits.Z_canonical_coupled_line_circuits-Tuple%7BInt64%2C%20Vararg%7BAny%2C%204%7D%7D)

```@raw html
<a id="JosephsonCircuits.Z_coupled_tline!-Tuple{AbstractMatrix, Vararg{Number, 4}}"></a>
```

[`JosephsonCircuits.Z_coupled_tline!`](api/internals.md#JosephsonCircuits.Z_coupled_tline!-Tuple%7BAbstractMatrix%2C%20Vararg%7BNumber%2C%204%7D%7D)

```@raw html
<a id="JosephsonCircuits.Z_coupled_tline-NTuple{4, Number}"></a>
```

[`JosephsonCircuits.Z_coupled_tline`](api/internals.md#JosephsonCircuits.Z_coupled_tline-NTuple%7B4%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Z_invC-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.Z_invC`](api/internals.md#JosephsonCircuits.Z_invC-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.Z_shuntZ!-Tuple{AbstractMatrix, Number}"></a>
```

[`JosephsonCircuits.Z_shuntZ!`](api/internals.md#JosephsonCircuits.Z_shuntZ!-Tuple%7BAbstractMatrix%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Z_shuntZ-Tuple{Number}"></a>
```

[`JosephsonCircuits.Z_shuntZ`](api/internals.md#JosephsonCircuits.Z_shuntZ-Tuple%7BNumber%7D)

```@raw html
<a id="JosephsonCircuits.Z_tline!-Tuple{AbstractMatrix, Number, Number}"></a>
```

[`JosephsonCircuits.Z_tline!`](api/internals.md#JosephsonCircuits.Z_tline!-Tuple%7BAbstractMatrix%2C%20Number%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.Z_tline-Tuple{Number, Number}"></a>
```

[`JosephsonCircuits.Z_tline`](api/internals.md#JosephsonCircuits.Z_tline-Tuple%7BNumber%2C%20Number%7D)

```@raw html
<a id="JosephsonCircuits.ZtoA"></a>
```

[`JosephsonCircuits.ZtoA`](api/internals.md#JosephsonCircuits.ZtoA)

```@raw html
<a id="JosephsonCircuits.ZtoA!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.ZtoA!`](api/internals.md#JosephsonCircuits.ZtoA!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.ZtoB"></a>
```

[`JosephsonCircuits.ZtoB`](api/internals.md#JosephsonCircuits.ZtoB)

```@raw html
<a id="JosephsonCircuits.ZtoB!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.ZtoB!`](api/internals.md#JosephsonCircuits.ZtoB!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.ZtoS"></a>
```

[`JosephsonCircuits.ZtoS`](api/internals.md#JosephsonCircuits.ZtoS)

```@raw html
<a id="JosephsonCircuits.ZtoS!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.ZtoS!`](api/internals.md#JosephsonCircuits.ZtoS!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits._harvestfloquet!-Union{Tuple{T}, Tuple{TJ}, Tuple{TI}, Tuple{JosephsonCircuits.FloquetPreconditioner{TI, TJ, T, TM, TV} where {TM&lt;:AbstractMatrix{T}, TV&lt;:AbstractVector{T}}, AbstractMatrix, AbstractMatrix}} where {TI, TJ, T}"></a>
```

[`JosephsonCircuits._harvestfloquet!`](api/internals.md#JosephsonCircuits._harvestfloquet!-Union%7BTuple%7BT%7D%2C%20Tuple%7BTJ%7D%2C%20Tuple%7BTI%7D%2C%20Tuple%7BJosephsonCircuits.FloquetPreconditioner%7BTI%2C%20TJ%2C%20T%2C%20TM%2C%20TV%7D%20where%20%7BTM%3C%3AAbstractMatrix%7BT%7D%2C%20TV%3C%3AAbstractVector%7BT%7D%7D%2C%20AbstractMatrix%2C%20AbstractMatrix%7D%7D%20where%20%7BTI%2C%20TJ%2C%20T%7D)

```@raw html
<a id="JosephsonCircuits._rebuildfloquet!-Union{Tuple{JosephsonCircuits.FloquetPreconditioner{TI, TJ, T, TM, TV} where {TM&lt;:AbstractMatrix{T}, TV&lt;:AbstractVector{T}}}, Tuple{T}, Tuple{TJ}, Tuple{TI}} where {TI, TJ, T}"></a>
```

[`JosephsonCircuits._rebuildfloquet!`](api/internals.md#JosephsonCircuits._rebuildfloquet!-Union%7BTuple%7BJosephsonCircuits.FloquetPreconditioner%7BTI%2C%20TJ%2C%20T%2C%20TM%2C%20TV%7D%20where%20%7BTM%3C%3AAbstractMatrix%7BT%7D%2C%20TV%3C%3AAbstractVector%7BT%7D%7D%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BTJ%7D%2C%20Tuple%7BTI%7D%7D%20where%20%7BTI%2C%20TJ%2C%20T%7D)

```@raw html
<a id="JosephsonCircuits.activemoderows"></a>
```

[`JosephsonCircuits.activemoderows`](api/internals.md#JosephsonCircuits.activemoderows)

```@raw html
<a id="JosephsonCircuits.add_modes-Union{Tuple{T}, Tuple{AbstractArray{Array{Tuple{T, Int64}, 1}, 1}, Integer}} where T"></a>
```

[`JosephsonCircuits.add_modes`](api/internals.md#JosephsonCircuits.add_modes-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BArray%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%2C%201%7D%2C%20Integer%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.add_modes-Union{Tuple{T}, Tuple{AbstractArray{Tuple{T, T, Int64, Int64}, 1}, Integer}} where T"></a>
```

[`JosephsonCircuits.add_modes`](api/internals.md#JosephsonCircuits.add_modes-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BTuple%7BT%2C%20T%2C%20Int64%2C%20Int64%7D%2C%201%7D%2C%20Integer%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.add_splitters-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{JosephsonCircuits.PassiveNetwork{T, N}, 1}, AbstractArray{Array{Tuple{T, Int64}, 1}, 1}}} where {T, N}"></a>
```

[`JosephsonCircuits.add_splitters`](api/internals.md#JosephsonCircuits.add_splitters-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BJosephsonCircuits.PassiveNetwork%7BT%2C%20N%7D%2C%201%7D%2C%20AbstractArray%7BArray%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%2C%201%7D%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.add_splitters-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{JosephsonCircuits.PassiveNetwork{T, N}, 1}, AbstractArray{Tuple{T, T, Int64, Int64}, 1}}} where {T, N}"></a>
```

[`JosephsonCircuits.add_splitters`](api/internals.md#JosephsonCircuits.add_splitters-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BJosephsonCircuits.PassiveNetwork%7BT%2C%20N%7D%2C%201%7D%2C%20AbstractArray%7BTuple%7BT%2C%20T%2C%20Int64%2C%20Int64%7D%2C%201%7D%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.addblockdc!-Tuple{AbstractVector, JosephsonCircuits.DCBlockRows, AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.addblockdc!`](api/internals.md#JosephsonCircuits.addblockdc!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.DCBlockRows%2C%20AbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.addblocktransport!-Tuple{AbstractVector, JosephsonCircuits.DCBlockRows, AbstractVector}"></a>
```

[`JosephsonCircuits.addblocktransport!`](api/internals.md#JosephsonCircuits.addblocktransport!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.DCBlockRows%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.addconstantsources!-NTuple{9, Any}"></a>
```

[`JosephsonCircuits.addconstantsources!`](api/internals.md#JosephsonCircuits.addconstantsources!-NTuple%7B9%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.addjosephsonterm!-Tuple{AbstractVector, JosephsonCircuits.StructureComplexJosephsonPlan, Any}"></a>
```

[`JosephsonCircuits.addjosephsonterm!`](api/internals.md#JosephsonCircuits.addjosephsonterm!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.StructureComplexJosephsonPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.addsources!-NTuple{11, Any}"></a>
```

[`JosephsonCircuits.addsources!`](api/internals.md#JosephsonCircuits.addsources!-NTuple%7B11%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.adjointdevice-Tuple{JosephsonCircuits.DeviceSweep, Integer}"></a>
```

[`JosephsonCircuits.adjointdevice`](api/internals.md#JosephsonCircuits.adjointdevice-Tuple%7BJosephsonCircuits.DeviceSweep%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.adjointnoisesigns!-Tuple{AbstractMatrix, Any, Integer}"></a>
```

[`JosephsonCircuits.adjointnoisesigns!`](api/internals.md#JosephsonCircuits.adjointnoisesigns!-Tuple%7BAbstractMatrix%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.adjointsolution!-Tuple{Any, JosephsonCircuits.DeviceSweep, Integer}"></a>
```

[`JosephsonCircuits.adjointsolution!`](api/internals.md#JosephsonCircuits.adjointsolution!-Tuple%7BAny%2C%20JosephsonCircuits.DeviceSweep%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.aliasmode-Union{Tuple{N}, Tuple{NTuple{N, Int64}, NTuple{N, Int64}}} where N"></a>
```

[`JosephsonCircuits.aliasmode`](api/internals.md#JosephsonCircuits.aliasmode-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.allsinusoidal-Tuple{JosephsonCircuits.JunctionRelations}"></a>
```

[`JosephsonCircuits.allsinusoidal`](api/internals.md#JosephsonCircuits.allsinusoidal-Tuple%7BJosephsonCircuits.JunctionRelations%7D)

```@raw html
<a id="JosephsonCircuits.amalgamate-Tuple{Any, Any, Any, Integer}"></a>
```

[`JosephsonCircuits.amalgamate`](api/internals.md#JosephsonCircuits.amalgamate-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.andersoncorrection!-Union{Tuple{T}, Tuple{JosephsonCircuits.AndersonState{T}, AbstractVector}} where T"></a>
```

[`JosephsonCircuits.andersoncorrection!`](api/internals.md#JosephsonCircuits.andersoncorrection!-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.AndersonState%7BT%7D%2C%20AbstractVector%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.andersonhistory!-Tuple{JosephsonCircuits.AndersonState, AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.andersonhistory!`](api/internals.md#JosephsonCircuits.andersonhistory!-Tuple%7BJosephsonCircuits.AndersonState%2C%20AbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.andersonrestart!-Tuple{JosephsonCircuits.AndersonState}"></a>
```

[`JosephsonCircuits.andersonrestart!`](api/internals.md#JosephsonCircuits.andersonrestart!-Tuple%7BJosephsonCircuits.AndersonState%7D)

```@raw html
<a id="JosephsonCircuits.applybackwardjosephsontranspose!-Tuple{AbstractArray, JosephsonCircuits.NonlinearTermTransposePlan, JosephsonCircuits.NonlinearTermPlan, AbstractVector}"></a>
```

[`JosephsonCircuits.applybackwardjosephsontranspose!`](api/internals.md#JosephsonCircuits.applybackwardjosephsontranspose!-Tuple%7BAbstractArray%2C%20JosephsonCircuits.NonlinearTermTransposePlan%2C%20JosephsonCircuits.NonlinearTermPlan%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.applybackwardterm!-Tuple{AbstractVector{&lt;:Real}, JosephsonCircuits.NonlinearTermPlan, AbstractArray, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.applybackwardterm!`](api/internals.md#JosephsonCircuits.applybackwardterm!-Tuple%7BAbstractVector%7B%3C%3AReal%7D%2C%20JosephsonCircuits.NonlinearTermPlan%2C%20AbstractArray%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.applydcconductance-Tuple{AbstractVector, JosephsonCircuits.DCConductancePlan, JosephsonCircuits.DCConductanceSolution, Integer}"></a>
```

[`JosephsonCircuits.applydcconductance`](api/internals.md#JosephsonCircuits.applydcconductance-Tuple%7BAbstractVector%2C%20JosephsonCircuits.DCConductancePlan%2C%20JosephsonCircuits.DCConductanceSolution%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.applydcsolve!-Tuple{AbstractVector, AbstractVector, JosephsonCircuits.DCFactorization}"></a>
```

[`JosephsonCircuits.applydcsolve!`](api/internals.md#JosephsonCircuits.applydcsolve!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20JosephsonCircuits.DCFactorization%7D)

```@raw html
<a id="JosephsonCircuits.applydcupdate!-Tuple{AbstractVector, AbstractVector, JosephsonCircuits.DCUpdate}"></a>
```

[`JosephsonCircuits.applydcupdate!`](api/internals.md#JosephsonCircuits.applydcupdate!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20JosephsonCircuits.DCUpdate%7D)

```@raw html
<a id="JosephsonCircuits.applyfft!-Union{Tuple{T}, Tuple{AbstractArray{Complex{T}}, AbstractArray{T}, Any}} where T"></a>
```

[`JosephsonCircuits.applyfft!`](api/internals.md#JosephsonCircuits.applyfft!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BComplex%7BT%7D%7D%2C%20AbstractArray%7BT%7D%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.applyffttranspose!-Union{Tuple{T}, Tuple{Array{Complex{T}}, Array{Complex{T}}, Array{Complex{T}}, Any}} where T"></a>
```

[`JosephsonCircuits.applyffttranspose!`](api/internals.md#JosephsonCircuits.applyffttranspose!-Union%7BTuple%7BT%7D%2C%20Tuple%7BArray%7BComplex%7BT%7D%7D%2C%20Array%7BComplex%7BT%7D%7D%2C%20Array%7BComplex%7BT%7D%7D%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.applyforwardterm!-Tuple{AbstractArray, JosephsonCircuits.NonlinearTermPlan, AbstractVector{&lt;:Real}}"></a>
```

[`JosephsonCircuits.applyforwardterm!`](api/internals.md#JosephsonCircuits.applyforwardterm!-Tuple%7BAbstractArray%2C%20JosephsonCircuits.NonlinearTermPlan%2C%20AbstractVector%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.applyforwardtranspose!-Tuple{AbstractVector, JosephsonCircuits.NonlinearTermTransposePlan, AbstractArray, AbstractVector}"></a>
```

[`JosephsonCircuits.applyforwardtranspose!`](api/internals.md#JosephsonCircuits.applyforwardtranspose!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.NonlinearTermTransposePlan%2C%20AbstractArray%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.applyifft!-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractArray{Complex{T}}, Any}} where T"></a>
```

[`JosephsonCircuits.applyifft!`](api/internals.md#JosephsonCircuits.applyifft!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%7D%2C%20AbstractArray%7BComplex%7BT%7D%7D%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.applynl!-Union{Tuple{T}, Tuple{AbstractArray{Complex{T}}, AbstractArray{T}, Any, Any, Any}} where T"></a>
```

[`JosephsonCircuits.applynl!`](api/internals.md#JosephsonCircuits.applynl!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BComplex%7BT%7D%7D%2C%20AbstractArray%7BT%7D%2C%20Any%2C%20Any%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.applypreconditioner!"></a>
```

[`JosephsonCircuits.applypreconditioner!`](api/internals.md#JosephsonCircuits.applypreconditioner!)

```@raw html
<a id="JosephsonCircuits.applyrealtocomplex!-Tuple{AbstractVector, JosephsonCircuits.NonlinearTermPlan, AbstractVector}"></a>
```

[`JosephsonCircuits.applyrealtocomplex!`](api/internals.md#JosephsonCircuits.applyrealtocomplex!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.NonlinearTermPlan%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.applyrelationnl!-Union{Tuple{F}, Tuple{T}, Tuple{AbstractArray{Complex{T}}, AbstractArray{T}, AbstractArray{T}, Any, Any, F, Any, Any}} where {T, F}"></a>
```

[`JosephsonCircuits.applyrelationnl!`](api/internals.md#JosephsonCircuits.applyrelationnl!-Union%7BTuple%7BF%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BComplex%7BT%7D%7D%2C%20AbstractArray%7BT%7D%2C%20AbstractArray%7BT%7D%2C%20Any%2C%20Any%2C%20F%2C%20Any%2C%20Any%7D%7D%20where%20%7BT%2C%20F%7D)

```@raw html
<a id="JosephsonCircuits.applyscatteringstamps!-Tuple{AbstractMatrix, JosephsonCircuits.DeviceScatteringStamps}"></a>
```

[`JosephsonCircuits.applyscatteringstamps!`](api/internals.md#JosephsonCircuits.applyscatteringstamps!-Tuple%7BAbstractMatrix%2C%20JosephsonCircuits.DeviceScatteringStamps%7D)

```@raw html
<a id="JosephsonCircuits.asoperator-Tuple{Function, Integer}"></a>
```

[`JosephsonCircuits.asoperator`](api/internals.md#JosephsonCircuits.asoperator-Tuple%7BFunction%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.assembleblocks!-Union{Tuple{T}, Tuple{JosephsonCircuits.ClusterBlocks{T}, JosephsonCircuits.BlockStructure, Any}} where T"></a>
```

[`JosephsonCircuits.assembleblocks!`](api/internals.md#JosephsonCircuits.assembleblocks!-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.ClusterBlocks%7BT%7D%2C%20JosephsonCircuits.BlockStructure%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.assemblebranch!-Union{Tuple{F}, Tuple{Vector, Any, JosephsonCircuits.BranchStampPlan, Any, F}} where F"></a>
```

[`JosephsonCircuits.assemblebranch!`](api/internals.md#JosephsonCircuits.assemblebranch!-Union%7BTuple%7BF%7D%2C%20Tuple%7BVector%2C%20Any%2C%20JosephsonCircuits.BranchStampPlan%2C%20Any%2C%20F%7D%7D%20where%20F)

```@raw html
<a id="JosephsonCircuits.assemblebranch-Union{Tuple{F}, Tuple{T}, Tuple{Type{T}, JosephsonCircuits.BranchStampPlan, Any, F, Integer}} where {T, F}"></a>
```

[`JosephsonCircuits.assemblebranch`](api/internals.md#JosephsonCircuits.assemblebranch-Union%7BTuple%7BF%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20JosephsonCircuits.BranchStampPlan%2C%20Any%2C%20F%2C%20Integer%7D%7D%20where%20%7BT%2C%20F%7D)

```@raw html
<a id="JosephsonCircuits.assemblecomplexjacobian!-Tuple{AbstractVector, JosephsonCircuits.StructureComplexJacobianPlan, Any}"></a>
```

[`JosephsonCircuits.assemblecomplexjacobian!`](api/internals.md#JosephsonCircuits.assemblecomplexjacobian!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.StructureComplexJacobianPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.assembleinvinductance-Union{Tuple{T}, Tuple{Type{T}, JosephsonCircuits.InverseInductancePlan, SparseVector, Integer}} where T"></a>
```

[`JosephsonCircuits.assembleinvinductance`](api/internals.md#JosephsonCircuits.assembleinvinductance-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20JosephsonCircuits.InverseInductancePlan%2C%20SparseVector%2C%20Integer%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.assemblematrices!"></a>
```

[`JosephsonCircuits.assemblematrices!`](api/internals.md#JosephsonCircuits.assemblematrices!)

```@raw html
<a id="JosephsonCircuits.assemblematrices-Tuple{JosephsonCircuits.CircuitMatrixPlan, JosephsonCircuits.BoundCircuit}"></a>
```

[`JosephsonCircuits.assemblematrices`](api/internals.md#JosephsonCircuits.assemblematrices-Tuple%7BJosephsonCircuits.CircuitMatrixPlan%2C%20JosephsonCircuits.BoundCircuit%7D)

```@raw html
<a id="JosephsonCircuits.assemblenodal!-Tuple{Vector, Vector{Bool}, JosephsonCircuits.NodalStampPlan, Any}"></a>
```

[`JosephsonCircuits.assemblenodal!`](api/internals.md#JosephsonCircuits.assemblenodal!-Tuple%7BVector%2C%20Vector%7BBool%7D%2C%20JosephsonCircuits.NodalStampPlan%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.assemblenodal-Union{Tuple{T}, Tuple{Type{T}, JosephsonCircuits.NodalStampPlan, Any, Integer}} where T"></a>
```

[`JosephsonCircuits.assemblenodal`](api/internals.md#JosephsonCircuits.assemblenodal-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20JosephsonCircuits.NodalStampPlan%2C%20Any%2C%20Integer%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.assemblerealjacobian!-Tuple{AbstractVector, JosephsonCircuits.StructureRealJacobianPlan, AbstractArray}"></a>
```

[`JosephsonCircuits.assemblerealjacobian!`](api/internals.md#JosephsonCircuits.assemblerealjacobian!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.StructureRealJacobianPlan%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.assemblescattering!"></a>
```

[`JosephsonCircuits.assemblescattering!`](api/internals.md#JosephsonCircuits.assemblescattering!)

```@raw html
<a id="JosephsonCircuits.assemblesweep!-Tuple{AbstractMatrix, JosephsonCircuits.FrequencySweepPlan, AbstractVector}"></a>
```

[`JosephsonCircuits.assemblesweep!`](api/internals.md#JosephsonCircuits.assemblesweep!-Tuple%7BAbstractMatrix%2C%20JosephsonCircuits.FrequencySweepPlan%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.assemblesystemmatrix!-Tuple{SparseMatrixCSC, JosephsonCircuits.HBLinearizedSystem, AbstractVector}"></a>
```

[`JosephsonCircuits.assemblesystemmatrix!`](api/internals.md#JosephsonCircuits.assemblesystemmatrix!-Tuple%7BSparseMatrixCSC%2C%20JosephsonCircuits.HBLinearizedSystem%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.atfrequency-Tuple{Number, Any}"></a>
```

[`JosephsonCircuits.atfrequency`](api/internals.md#JosephsonCircuits.atfrequency-Tuple%7BNumber%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.autonne_takagi-Tuple{AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.autonne_takagi`](api/internals.md#JosephsonCircuits.autonne_takagi-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.autonne_takagi-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.autonne_takagi`](api/internals.md#JosephsonCircuits.autonne_takagi-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.auxcurrentscale-Tuple{Any}"></a>
```

[`JosephsonCircuits.auxcurrentscale`](api/internals.md#JosephsonCircuits.auxcurrentscale-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.backtracking_linesearch!-Tuple{Any, AbstractVector, AbstractVector, AbstractVector, AbstractVector, Real, Real}"></a>
```

[`JosephsonCircuits.backtracking_linesearch!`](api/internals.md#JosephsonCircuits.backtracking_linesearch!-Tuple%7BAny%2C%20AbstractVector%2C%20AbstractVector%2C%20AbstractVector%2C%20AbstractVector%2C%20Real%2C%20Real%7D)

```@raw html
<a id="JosephsonCircuits.backwardjosephsontransposekernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.backwardjosephsontransposekernel!`](api/internals.md#JosephsonCircuits.backwardjosephsontransposekernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.backwardtermkernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.backwardtermkernel!`](api/internals.md#JosephsonCircuits.backwardtermkernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.backwardtermkernelcomplex!-Tuple{Any}"></a>
```

[`JosephsonCircuits.backwardtermkernelcomplex!`](api/internals.md#JosephsonCircuits.backwardtermkernelcomplex!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.batchedinverse!-Union{Tuple{T}, Tuple{AbstractArray{T, 3}, AbstractArray{T, 3}, AbstractArray{T, 3}, KernelAbstractions.CPU}} where T"></a>
```

[`JosephsonCircuits.batchedinverse!`](api/internals.md#JosephsonCircuits.batchedinverse!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%203%7D%2C%20AbstractArray%7BT%2C%203%7D%2C%20AbstractArray%7BT%2C%203%7D%2C%20KernelAbstractions.CPU%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.batchedmul!-Union{Tuple{T}, Tuple{AbstractArray{T, 3}, AbstractArray{T, 3}, AbstractArray{T, 3}, Any, Any, Bool, Bool, KernelAbstractions.CPU}} where T"></a>
```

[`JosephsonCircuits.batchedmul!`](api/internals.md#JosephsonCircuits.batchedmul!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%203%7D%2C%20AbstractArray%7BT%2C%203%7D%2C%20AbstractArray%7BT%2C%203%7D%2C%20Any%2C%20Any%2C%20Bool%2C%20Bool%2C%20KernelAbstractions.CPU%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.bathamplitude-Tuple{JosephsonCircuits.TransientNoiseBath, Any, Any}"></a>
```

[`JosephsonCircuits.bathamplitude`](api/internals.md#JosephsonCircuits.bathamplitude-Tuple%7BJosephsonCircuits.TransientNoiseBath%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.bathfamily-Tuple{LinearizedScattering, Any}"></a>
```

[`JosephsonCircuits.bathfamily`](api/internals.md#JosephsonCircuits.bathfamily-Tuple%7BLinearizedScattering%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.bindvalues-Tuple{JosephsonCircuits.CompiledCircuit, Any}"></a>
```

[`JosephsonCircuits.bindvalues`](api/internals.md#JosephsonCircuits.bindvalues-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.bloch_messiah_block-Tuple{AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.bloch_messiah_block`](api/internals.md#JosephsonCircuits.bloch_messiah_block-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.block_to_pair-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.block_to_pair`](api/internals.md#JosephsonCircuits.block_to_pair-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.block_to_pair-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.block_to_pair`](api/internals.md#JosephsonCircuits.block_to_pair-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.block_to_pair2-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.block_to_pair2`](api/internals.md#JosephsonCircuits.block_to_pair2-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.block_to_pair2-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.block_to_pair2`](api/internals.md#JosephsonCircuits.block_to_pair2-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.block_to_pair_perm-Tuple{Integer}"></a>
```

[`JosephsonCircuits.block_to_pair_perm`](api/internals.md#JosephsonCircuits.block_to_pair_perm-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.blockfactorbytes-Union{Tuple{T}, Tuple{Type{T}, AbstractMatrix{Bool}, Any, Any, Integer, JosephsonCircuits.ModeLayout}} where T"></a>
```

[`JosephsonCircuits.blockfactorbytes`](api/internals.md#JosephsonCircuits.blockfactorbytes-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20AbstractMatrix%7BBool%7D%2C%20Any%2C%20Any%2C%20Integer%2C%20JosephsonCircuits.ModeLayout%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.blocklu!-Union{Tuple{T}, Tuple{JosephsonCircuits.BlockLU{T}, Any}} where T"></a>
```

[`JosephsonCircuits.blocklu!`](api/internals.md#JosephsonCircuits.blocklu!-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.BlockLU%7BT%7D%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.blocklu-Union{Tuple{T}, Tuple{Type{T}, Any, Any}} where T"></a>
```

[`JosephsonCircuits.blocklu`](api/internals.md#JosephsonCircuits.blocklu-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20Any%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.blocknodegraph-Tuple{SparseMatrixCSC, Integer}"></a>
```

[`JosephsonCircuits.blocknodegraph`](api/internals.md#JosephsonCircuits.blocknodegraph-Tuple%7BSparseMatrixCSC%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.blocknoisecontractkernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.blocknoisecontractkernel!`](api/internals.md#JosephsonCircuits.blocknoisecontractkernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.blocknoisefactorkernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.blocknoisefactorkernel!`](api/internals.md#JosephsonCircuits.blocknoisefactorkernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.blockresidual!-Tuple{AbstractArray{&lt;:Any, 3}, JosephsonCircuits.SparseBlockFactorization, AbstractArray{&lt;:Any, 3}, AbstractMatrix}"></a>
```

[`JosephsonCircuits.blockresidual!`](api/internals.md#JosephsonCircuits.blockresidual!-Tuple%7BAbstractArray%7B%3C%3AAny%2C%203%7D%2C%20JosephsonCircuits.SparseBlockFactorization%2C%20AbstractArray%7B%3C%3AAny%2C%203%7D%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.blocksensitivitystamp-Tuple{Any, Any, Integer}"></a>
```

[`JosephsonCircuits.blocksensitivitystamp`](api/internals.md#JosephsonCircuits.blocksensitivitystamp-Tuple%7BAny%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.blocksolve!-Union{Tuple{T}, Tuple{AbstractArray{&lt;:Any, 3}, JosephsonCircuits.SparseBlockFactorization{T}, AbstractArray}} where T"></a>
```

[`JosephsonCircuits.blocksolve!`](api/internals.md#JosephsonCircuits.blocksolve!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7B%3C%3AAny%2C%203%7D%2C%20JosephsonCircuits.SparseBlockFactorization%7BT%7D%2C%20AbstractArray%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.blockstampvals!-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.blockstampvals!`](api/internals.md#JosephsonCircuits.blockstampvals!-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.blockstructure-Union{Tuple{T}, Tuple{Type{T}, Any, Matrix, Matrix, AbstractMatrix{Bool}, SparseMatrixCSC, Integer, Integer, Integer, JosephsonCircuits.ModeLayout, Any}} where T"></a>
```

[`JosephsonCircuits.blockstructure`](api/internals.md#JosephsonCircuits.blockstructure-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20Any%2C%20Matrix%2C%20Matrix%2C%20AbstractMatrix%7BBool%7D%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Integer%2C%20JosephsonCircuits.ModeLayout%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.blocksystembytes-Union{Tuple{T}, Tuple{Type{T}, Any}} where T"></a>
```

[`JosephsonCircuits.blocksystembytes`](api/internals.md#JosephsonCircuits.blocksystembytes-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.branchendpoints-Tuple{SparseMatrixCSC, Int64}"></a>
```

[`JosephsonCircuits.branchendpoints`](api/internals.md#JosephsonCircuits.branchendpoints-Tuple%7BSparseMatrixCSC%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.branchnodesandsigns-Tuple{SparseMatrixCSC, Integer, Integer}"></a>
```

[`JosephsonCircuits.branchnodesandsigns`](api/internals.md#JosephsonCircuits.branchnodesandsigns-Tuple%7BSparseMatrixCSC%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.branchstampplan-Tuple{JosephsonCircuits.CompiledCircuit, Vector{Int64}, Dict, Int64}"></a>
```

[`JosephsonCircuits.branchstampplan`](api/internals.md#JosephsonCircuits.branchstampplan-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Vector%7BInt64%7D%2C%20Dict%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.buildcoupling-Tuple{JosephsonCircuits.PreconditionerPlan, JosephsonCircuits.AbstractModeCoupling, Any, Any}"></a>
```

[`JosephsonCircuits.buildcoupling`](api/internals.md#JosephsonCircuits.buildcoupling-Tuple%7BJosephsonCircuits.PreconditionerPlan%2C%20JosephsonCircuits.AbstractModeCoupling%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcAmna-Tuple{Vector{Int64}, Int64}"></a>
```

[`JosephsonCircuits.calcAmna`](api/internals.md#JosephsonCircuits.calcAmna-Tuple%7BVector%7BInt64%7D%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.calcAmnaind-Tuple{Vector{Int64}, SparseVector, SparseMatrixCSC, SparseMatrixCSC, Int64, Int64, Int64, Any}"></a>
```

[`JosephsonCircuits.calcAmnaind`](api/internals.md#JosephsonCircuits.calcAmnaind-Tuple%7BVector%7BInt64%7D%2C%20SparseVector%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Int64%2C%20Int64%2C%20Int64%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcCjIcmean-Tuple{AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.calcCjIcmean`](api/internals.md#JosephsonCircuits.calcCjIcmean-Tuple%7BAbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.calcCnoise!-Tuple{AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.calcCnoise!`](api/internals.md#JosephsonCircuits.calcCnoise!-Tuple%7BAbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcCnoise!-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.calcCnoise!`](api/internals.md#JosephsonCircuits.calcCnoise!-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcCnoise-Union{Tuple{AbstractMatrix{T}}, Tuple{T}} where T"></a>
```

[`JosephsonCircuits.calcCnoise`](api/internals.md#JosephsonCircuits.calcCnoise-Union%7BTuple%7BAbstractMatrix%7BT%7D%7D%2C%20Tuple%7BT%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.calcCnoise-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractArray{T}}} where T"></a>
```

[`JosephsonCircuits.calcCnoise`](api/internals.md#JosephsonCircuits.calcCnoise-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%7D%2C%20AbstractArray%7BT%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.calcSsensitivity!-NTuple{13, Any}"></a>
```

[`JosephsonCircuits.calcSsensitivity!`](api/internals.md#JosephsonCircuits.calcSsensitivity!-NTuple%7B13%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcSsensitivityreverse!-Tuple{Any, JosephsonCircuits.ReverseSensitivity, Any, Any, Any, Any, Any, Any, JosephsonCircuits.ReverseSensitivityBuffers}"></a>
```

[`JosephsonCircuits.calcSsensitivityreverse!`](api/internals.md#JosephsonCircuits.calcSsensitivityreverse!-Tuple%7BAny%2C%20JosephsonCircuits.ReverseSensitivity%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20JosephsonCircuits.ReverseSensitivityBuffers%7D)

```@raw html
<a id="JosephsonCircuits.calc_noise_covariances-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.calc_noise_covariances`](api/internals.md#JosephsonCircuits.calc_noise_covariances-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.calcblockresidualsensitivity-Tuple{JosephsonCircuits.HBOperatingPoint, JosephsonCircuits.CompiledCircuit, AbstractVector}"></a>
```

[`JosephsonCircuits.calcblockresidualsensitivity`](api/internals.md#JosephsonCircuits.calcblockresidualsensitivity-Tuple%7BJosephsonCircuits.HBOperatingPoint%2C%20JosephsonCircuits.CompiledCircuit%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.calcbranchtimedomainmap-Tuple{Any, Integer, Integer}"></a>
```

[`JosephsonCircuits.calcbranchtimedomainmap`](api/internals.md#JosephsonCircuits.calcbranchtimedomainmap-Tuple%7BAny%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.calccm!"></a>
```

[`JosephsonCircuits.calccm!`](api/internals.md#JosephsonCircuits.calccm!)

```@raw html
<a id="JosephsonCircuits.calcdcgaugeindices-Tuple{Vector{Vector{Int64}}, Vector, Int64}"></a>
```

[`JosephsonCircuits.calcdcgaugeindices`](api/internals.md#JosephsonCircuits.calcdcgaugeindices-Tuple%7BVector%7BVector%7BInt64%7D%7D%2C%20Vector%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.calcfreqs-Union{Tuple{N}, Tuple{NTuple{N, Int64}, NTuple{N, Int64}, NTuple{N, Int64}}} where N"></a>
```

[`JosephsonCircuits.calcfreqs`](api/internals.md#JosephsonCircuits.calcfreqs-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.calcfreqsdft-Union{Tuple{NTuple{N, Int64}}, Tuple{N}} where N"></a>
```

[`JosephsonCircuits.calcfreqsdft`](api/internals.md#JosephsonCircuits.calcfreqsdft-Union%7BTuple%7BNTuple%7BN%2C%20Int64%7D%7D%2C%20Tuple%7BN%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.calcfreqsrdft-Union{Tuple{NTuple{N, Int64}}, Tuple{N}} where N"></a>
```

[`JosephsonCircuits.calcfreqsrdft`](api/internals.md#JosephsonCircuits.calcfreqsrdft-Union%7BTuple%7BNTuple%7BN%2C%20Int64%7D%7D%2C%20Tuple%7BN%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.calcgraphs-Tuple{Vector{Tuple{Int64, Int64}}, Int64}"></a>
```

[`JosephsonCircuits.calcgraphs`](api/internals.md#JosephsonCircuits.calcgraphs-Tuple%7BVector%7BTuple%7BInt64%2C%20Int64%7D%7D%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.calcimpedance-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.calcimpedance`](api/internals.md#JosephsonCircuits.calcimpedance-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcinputwaves!-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.calcinputwaves!`](api/internals.md#JosephsonCircuits.calcinputwaves!-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcjunctionrelations-Tuple{Vector{Symbol}, Matrix{Int64}, AbstractDict, Dict, SparseVector}"></a>
```

[`JosephsonCircuits.calcjunctionrelations`](api/internals.md#JosephsonCircuits.calcjunctionrelations-Tuple%7BVector%7BSymbol%7D%2C%20Matrix%7BInt64%7D%2C%20AbstractDict%2C%20Dict%2C%20SparseVector%7D)

```@raw html
<a id="JosephsonCircuits.calcmodefreqs-Union{Tuple{N}, Tuple{NTuple{N, Any}, Array{NTuple{N, Int64}, 1}}} where N"></a>
```

[`JosephsonCircuits.calcmodefreqs`](api/internals.md#JosephsonCircuits.calcmodefreqs-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Any%7D%2C%20Array%7BNTuple%7BN%2C%20Int64%7D%2C%201%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.calcnodefluxsensitivity-Tuple{JosephsonCircuits.HBOperatingPoint, AbstractMatrix}"></a>
```

[`JosephsonCircuits.calcnodefluxsensitivity`](api/internals.md#JosephsonCircuits.calcnodefluxsensitivity-Tuple%7BJosephsonCircuits.HBOperatingPoint%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.calcnodesorting-Tuple{Vector{String}}"></a>
```

[`JosephsonCircuits.calcnodesorting`](api/internals.md#JosephsonCircuits.calcnodesorting-Tuple%7BVector%7BString%7D%7D)

```@raw html
<a id="JosephsonCircuits.calcnoisecovariance!"></a>
```

[`JosephsonCircuits.calcnoisecovariance!`](api/internals.md#JosephsonCircuits.calcnoisecovariance!)

```@raw html
<a id="JosephsonCircuits.calcoperatingpointstamps-Tuple{JosephsonCircuits.HBOperatingPoint, Any, Any}"></a>
```

[`JosephsonCircuits.calcoperatingpointstamps`](api/internals.md#JosephsonCircuits.calcoperatingpointstamps-Tuple%7BJosephsonCircuits.HBOperatingPoint%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcoutputwaves!-NTuple{8, Any}"></a>
```

[`JosephsonCircuits.calcoutputwaves!`](api/internals.md#JosephsonCircuits.calcoutputwaves!-NTuple%7B8%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcphiindices-Union{Tuple{N}, Tuple{JosephsonCircuits.Frequencies{N}, Dict{CartesianIndex{N}, CartesianIndex{N}}}} where N"></a>
```

[`JosephsonCircuits.calcphiindices`](api/internals.md#JosephsonCircuits.calcphiindices-Union%7BTuple%7BN%7D%2C%20Tuple%7BJosephsonCircuits.Frequencies%7BN%7D%2C%20Dict%7BCartesianIndex%7BN%7D%2C%20CartesianIndex%7BN%7D%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.calcportvoltage-NTuple{7, Any}"></a>
```

[`JosephsonCircuits.calcportvoltage`](api/internals.md#JosephsonCircuits.calcportvoltage-NTuple%7B7%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcqe!"></a>
```

[`JosephsonCircuits.calcqe!`](api/internals.md#JosephsonCircuits.calcqe!)

```@raw html
<a id="JosephsonCircuits.calcqe_S_Cnoise!-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.calcqe_S_Cnoise!`](api/internals.md#JosephsonCircuits.calcqe_S_Cnoise!-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcqe_S_Cnoise-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractArray{T}}} where T"></a>
```

[`JosephsonCircuits.calcqe_S_Cnoise`](api/internals.md#JosephsonCircuits.calcqe_S_Cnoise-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%7D%2C%20AbstractArray%7BT%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.calcqeideal!-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.calcqeideal!`](api/internals.md#JosephsonCircuits.calcqeideal!-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcqeideal-Union{Tuple{AbstractArray{T}}, Tuple{T}} where T"></a>
```

[`JosephsonCircuits.calcqeideal`](api/internals.md#JosephsonCircuits.calcqeideal-Union%7BTuple%7BAbstractArray%7BT%7D%7D%2C%20Tuple%7BT%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.calcresidualsensitivity"></a>
```

[`JosephsonCircuits.calcresidualsensitivity`](api/internals.md#JosephsonCircuits.calcresidualsensitivity)

```@raw html
<a id="JosephsonCircuits.calcscatteringmatrix!-Tuple{Any, AbstractVector, AbstractMatrix}"></a>
```

[`JosephsonCircuits.calcscatteringmatrix!`](api/internals.md#JosephsonCircuits.calcscatteringmatrix!-Tuple%7BAny%2C%20AbstractVector%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.calcscatteringmatrix!-Tuple{Any, Vector, Vector}"></a>
```

[`JosephsonCircuits.calcscatteringmatrix!`](api/internals.md#JosephsonCircuits.calcscatteringmatrix!-Tuple%7BAny%2C%20Vector%2C%20Vector%7D)

```@raw html
<a id="JosephsonCircuits.calcsensitivityscaling!-NTuple{9, Any}"></a>
```

[`JosephsonCircuits.calcsensitivityscaling!`](api/internals.md#JosephsonCircuits.calcsensitivityscaling!-NTuple%7B9%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcsensitivitystamps-Tuple{Any, JosephsonCircuits.CompiledCircuit, JosephsonCircuits.CircuitMatrices, Vararg{Any, 4}}"></a>
```

[`JosephsonCircuits.calcsensitivitystamps`](api/internals.md#JosephsonCircuits.calcsensitivitystamps-Tuple%7BAny%2C%20JosephsonCircuits.CompiledCircuit%2C%20JosephsonCircuits.CircuitMatrices%2C%20Vararg%7BAny%2C%204%7D%7D)

```@raw html
<a id="JosephsonCircuits.calcsolverscale-Tuple{Any, Vector{Symbol}, Vector, Vector, Any}"></a>
```

[`JosephsonCircuits.calcsolverscale`](api/internals.md#JosephsonCircuits.calcsolverscale-Tuple%7BAny%2C%20Vector%7BSymbol%7D%2C%20Vector%2C%20Vector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcsourcecurrent-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.calcsourcecurrent`](api/internals.md#JosephsonCircuits.calcsourcecurrent-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcsources-NTuple{10, Any}"></a>
```

[`JosephsonCircuits.calcsources`](api/internals.md#JosephsonCircuits.calcsources-NTuple%7B10%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.calcspicesortperms-Tuple{Dict{String, Vector{String}}}"></a>
```

[`JosephsonCircuits.calcspicesortperms`](api/internals.md#JosephsonCircuits.calcspicesortperms-Tuple%7BDict%7BString%2C%20Vector%7BString%7D%7D%7D)

```@raw html
<a id="JosephsonCircuits.calcstaticfluxcomponents-Tuple{Vector{Symbol}, Matrix{Int64}, Vector, Int64}"></a>
```

[`JosephsonCircuits.calcstaticfluxcomponents`](api/internals.md#JosephsonCircuits.calcstaticfluxcomponents-Tuple%7BVector%7BSymbol%7D%2C%20Matrix%7BInt64%7D%2C%20Vector%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.candeviceevaluate-Tuple{Any}"></a>
```

[`JosephsonCircuits.candeviceevaluate`](api/internals.md#JosephsonCircuits.candeviceevaluate-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.candidatecount-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.candidatecount`](api/internals.md#JosephsonCircuits.candidatecount-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.canonical_coupled_line_circuits-Tuple{Int64, Vararg{Any, 4}}"></a>
```

[`JosephsonCircuits.canonical_coupled_line_circuits`](api/internals.md#JosephsonCircuits.canonical_coupled_line_circuits-Tuple%7BInt64%2C%20Vararg%7BAny%2C%204%7D%7D)

```@raw html
<a id="JosephsonCircuits.canonicaldim-Tuple{JosephsonCircuits.CompositeLayout}"></a>
```

[`JosephsonCircuits.canonicaldim`](api/internals.md#JosephsonCircuits.canonicaldim-Tuple%7BJosephsonCircuits.CompositeLayout%7D)

```@raw html
<a id="JosephsonCircuits.canonicalfj-Tuple{Any, JosephsonCircuits.CanonicalWork, Any, JosephsonCircuits.CanonicalJacobianPlan}"></a>
```

[`JosephsonCircuits.canonicalfj`](api/internals.md#JosephsonCircuits.canonicalfj-Tuple%7BAny%2C%20JosephsonCircuits.CanonicalWork%2C%20Any%2C%20JosephsonCircuits.CanonicalJacobianPlan%7D)

```@raw html
<a id="JosephsonCircuits.canonicaljacobian!-Tuple{JosephsonCircuits.CanonicalJacobianPlan, SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.canonicaljacobian!`](api/internals.md#JosephsonCircuits.canonicaljacobian!-Tuple%7BJosephsonCircuits.CanonicalJacobianPlan%2C%20SparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.canonicaljacobianplan-Tuple{SparseMatrixCSC, JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.canonicaljacobianplan`](api/internals.md#JosephsonCircuits.canonicaljacobianplan-Tuple%7BSparseMatrixCSC%2C%20JosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.canonicaljvp-Tuple{Any, JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.canonicaljvp`](api/internals.md#JosephsonCircuits.canonicaljvp-Tuple%7BAny%2C%20JosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.canonicalresidual-Tuple{Any, JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.canonicalresidual`](api/internals.md#JosephsonCircuits.canonicalresidual-Tuple%7BAny%2C%20JosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.canonicalresidual-Tuple{SparseMatrixCSC, Vector{Float64}, Integer}"></a>
```

[`JosephsonCircuits.canonicalresidual`](api/internals.md#JosephsonCircuits.canonicalresidual-Tuple%7BSparseMatrixCSC%2C%20Vector%7BFloat64%7D%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.canonicalresult!-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.canonicalresult!`](api/internals.md#JosephsonCircuits.canonicalresult!-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.cansweepondevice-Tuple{JosephsonCircuits.HBLinearizedSystem}"></a>
```

[`JosephsonCircuits.cansweepondevice`](api/internals.md#JosephsonCircuits.cansweepondevice-Tuple%7BJosephsonCircuits.HBLinearizedSystem%7D)

```@raw html
<a id="JosephsonCircuits.cascadeS!-Tuple{AbstractMatrix, AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.cascadeS!`](api/internals.md#JosephsonCircuits.cascadeS!-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.cascadeS-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.cascadeS`](api/internals.md#JosephsonCircuits.cascadeS-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkblockcontract-Tuple{ScatteringParameters, Any, Any, Any}"></a>
```

[`JosephsonCircuits.checkblockcontract`](api/internals.md#JosephsonCircuits.checkblockcontract-Tuple%7BScatteringParameters%2C%20Any%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkblockdeclarations-Tuple{Union{Nothing, JosephsonCircuits.ScatteringStampSystem}, Any, Any}"></a>
```

[`JosephsonCircuits.checkblockdeclarations`](api/internals.md#JosephsonCircuits.checkblockdeclarations-Tuple%7BUnion%7BNothing%2C%20JosephsonCircuits.ScatteringStampSystem%7D%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkblockpair-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.checkblockpair`](api/internals.md#JosephsonCircuits.checkblockpair-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkcachekwargs-Tuple{Any}"></a>
```

[`JosephsonCircuits.checkcachekwargs`](api/internals.md#JosephsonCircuits.checkcachekwargs-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.checkcomponenttypes-Tuple{Vector{String}}"></a>
```

[`JosephsonCircuits.checkcomponenttypes`](api/internals.md#JosephsonCircuits.checkcomponenttypes-Tuple%7BVector%7BString%7D%7D)

```@raw html
<a id="JosephsonCircuits.checkcomponentvaluesdefined-Tuple{Vector, Vector}"></a>
```

[`JosephsonCircuits.checkcomponentvaluesdefined`](api/internals.md#JosephsonCircuits.checkcomponentvaluesdefined-Tuple%7BVector%2C%20Vector%7D)

```@raw html
<a id="JosephsonCircuits.checkcoupledloss-Tuple{JosephsonCircuits.CompiledCircuit, Any}"></a>
```

[`JosephsonCircuits.checkcoupledloss`](api/internals.md#JosephsonCircuits.checkcoupledloss-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkissymbolic-Tuple{Any}"></a>
```

[`JosephsonCircuits.checkissymbolic`](api/internals.md#JosephsonCircuits.checkissymbolic-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.checkjunctioncurrents-Tuple{JosephsonCircuits.NonlinearHB, JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.checkjunctioncurrents`](api/internals.md#JosephsonCircuits.checkjunctioncurrents-Tuple%7BJosephsonCircuits.NonlinearHB%2C%20JosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.checkjunctiondc"></a>
```

[`JosephsonCircuits.checkjunctiondc`](api/internals.md#JosephsonCircuits.checkjunctiondc)

```@raw html
<a id="JosephsonCircuits.checknoform-Tuple{Symbol, AbstractString}"></a>
```

[`JosephsonCircuits.checknoform`](api/internals.md#JosephsonCircuits.checknoform-Tuple%7BSymbol%2C%20AbstractString%7D)

```@raw html
<a id="JosephsonCircuits.checkpumpconjugates-Tuple{Vector{JosephsonCircuits.StampedScatteringBlock}, AbstractVector}"></a>
```

[`JosephsonCircuits.checkpumpconjugates`](api/internals.md#JosephsonCircuits.checkpumpconjugates-Tuple%7BVector%7BJosephsonCircuits.StampedScatteringBlock%7D%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.checkpumpedblock-Tuple{LinearizedScattering, AbstractVector, AbstractMatrix{Int64}, Any, Any}"></a>
```

[`JosephsonCircuits.checkpumpedblock`](api/internals.md#JosephsonCircuits.checkpumpedblock-Tuple%7BLinearizedScattering%2C%20AbstractVector%2C%20AbstractMatrix%7BInt64%7D%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkpumpedblockmodels-Tuple{Union{Nothing, JosephsonCircuits.ScatteringStampSystem}, Any, Any}"></a>
```

[`JosephsonCircuits.checkpumpedblockmodels`](api/internals.md#JosephsonCircuits.checkpumpedblockmodels-Tuple%7BUnion%7BNothing%2C%20JosephsonCircuits.ScatteringStampSystem%7D%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkpumpedblocks-Tuple{JosephsonCircuits.TransientProblem, Any}"></a>
```

[`JosephsonCircuits.checkpumpedblocks`](api/internals.md#JosephsonCircuits.checkpumpedblocks-Tuple%7BJosephsonCircuits.TransientProblem%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.checkstaticstiffnessvalues-Tuple{Vector{Symbol}, Vector}"></a>
```

[`JosephsonCircuits.checkstaticstiffnessvalues`](api/internals.md#JosephsonCircuits.checkstaticstiffnessvalues-Tuple%7BVector%7BSymbol%7D%2C%20Vector%7D)

```@raw html
<a id="JosephsonCircuits.checksweepoptions-Tuple{JosephsonCircuits.CompiledCircuit, Any, Any, Any, Any, Any, Any, Any, Any, Bool}"></a>
```

[`JosephsonCircuits.checksweepoptions`](api/internals.md#JosephsonCircuits.checksweepoptions-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Bool%7D)

```@raw html
<a id="JosephsonCircuits.circuitmatrixplan-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.circuitmatrixplan`](api/internals.md#JosephsonCircuits.circuitmatrixplan-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.circuitnodegraph-Tuple{Any, Any, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, Integer, Integer}"></a>
```

[`JosephsonCircuits.circuitnodegraph`](api/internals.md#JosephsonCircuits.circuitnodegraph-Tuple%7BAny%2C%20Any%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.circuitorder-Tuple{Any, SparseMatrixCSC, Integer, Integer, JosephsonCircuits.ModeLayout}"></a>
```

[`JosephsonCircuits.circuitorder`](api/internals.md#JosephsonCircuits.circuitorder-Tuple%7BAny%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20JosephsonCircuits.ModeLayout%7D)

```@raw html
<a id="JosephsonCircuits.circuittopology-Tuple{Vector{Symbol}, Matrix{Int64}, Int64}"></a>
```

[`JosephsonCircuits.circuittopology`](api/internals.md#JosephsonCircuits.circuittopology-Tuple%7BVector%7BSymbol%7D%2C%20Matrix%7BInt64%7D%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.circuitvariables-Tuple{Any}"></a>
```

[`JosephsonCircuits.circuitvariables`](api/internals.md#JosephsonCircuits.circuitvariables-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.clusterblocks-Union{Tuple{T}, Tuple{Type{T}, Any, Any, Any, Integer, JosephsonCircuits.ModeLayout, Any}} where T"></a>
```

[`JosephsonCircuits.clusterblocks`](api/internals.md#JosephsonCircuits.clusterblocks-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20Any%2C%20Any%2C%20Any%2C%20Integer%2C%20JosephsonCircuits.ModeLayout%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.clustersolve!-Union{Tuple{T}, Tuple{AbstractVector, JosephsonCircuits.ClusterBlocks{T}, AbstractVector, Any}} where T"></a>
```

[`JosephsonCircuits.clustersolve!`](api/internals.md#JosephsonCircuits.clustersolve!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractVector%2C%20JosephsonCircuits.ClusterBlocks%7BT%7D%2C%20AbstractVector%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.clustersymbolic-Tuple{Any, Any, Any, Integer, JosephsonCircuits.ModeLayout}"></a>
```

[`JosephsonCircuits.clustersymbolic`](api/internals.md#JosephsonCircuits.clustersymbolic-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Integer%2C%20JosephsonCircuits.ModeLayout%7D)

```@raw html
<a id="JosephsonCircuits.columnindices-Tuple{JosephsonCircuits.DeviceValuedSparseMatrix{var&quot;#s375&quot;, var&quot;#s374&quot;, V} where {var&quot;#s375&quot;, var&quot;#s374&quot;&lt;:SparseMatrixCSC, V&lt;:AbstractVector{var&quot;#s375&quot;}}}"></a>
```

[`JosephsonCircuits.columnindices`](api/internals.md#JosephsonCircuits.columnindices-Tuple%7BJosephsonCircuits.DeviceValuedSparseMatrix%7Bvar%22%23s375%22%2C%20var%22%23s374%22%2C%20V%7D%20where%20%7Bvar%22%23s375%22%2C%20var%22%23s374%22%3C%3ASparseMatrixCSC%2C%20V%3C%3AAbstractVector%7Bvar%22%23s375%22%7D%7D%7D)

```@raw html
<a id="JosephsonCircuits.commutationmargin-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.commutationmargin`](api/internals.md#JosephsonCircuits.commutationmargin-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.compare-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.compare`](api/internals.md#JosephsonCircuits.compare-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.comparearray-Union{Tuple{T}, Tuple{AbstractArray{T}, AbstractArray{T}}} where T"></a>
```

[`JosephsonCircuits.comparearray`](api/internals.md#JosephsonCircuits.comparearray-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%7D%2C%20AbstractArray%7BT%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.comparestruct-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.comparestruct`](api/internals.md#JosephsonCircuits.comparestruct-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.completecovariance-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.completecovariance`](api/internals.md#JosephsonCircuits.completecovariance-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.completedcovariance-Tuple{LinearizedScattering, AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.completedcovariance`](api/internals.md#JosephsonCircuits.completedcovariance-Tuple%7BLinearizedScattering%2C%20AbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.completepositivitymargin-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.completepositivitymargin`](api/internals.md#JosephsonCircuits.completepositivitymargin-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.complex_to_real!-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractArray{Complex{T}, 1}, AbstractVector{Bool}}} where T&lt;:Real"></a>
```

[`JosephsonCircuits.complex_to_real!`](api/internals.md#JosephsonCircuits.complex_to_real!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractVector%7BT%7D%2C%20AbstractArray%7BComplex%7BT%7D%2C%201%7D%2C%20AbstractVector%7BBool%7D%7D%7D%20where%20T%3C%3AReal)

```@raw html
<a id="JosephsonCircuits.complex_to_real-Union{Tuple{Tj}, Tuple{Ti}, Tuple{T}, Tuple{SparseMatrixCSC{Complex{T}, Ti}, JosephsonCircuits.ModeLayout, JosephsonCircuits.ModeLayout}, Tuple{SparseMatrixCSC{Complex{T}, Ti}, JosephsonCircuits.ModeLayout, JosephsonCircuits.ModeLayout, Type{Tj}}} where {T&lt;:Real, Ti, Tj&lt;:Integer}"></a>
```

[`JosephsonCircuits.complex_to_real`](api/internals.md#JosephsonCircuits.complex_to_real-Union%7BTuple%7BTj%7D%2C%20Tuple%7BTi%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BSparseMatrixCSC%7BComplex%7BT%7D%2C%20Ti%7D%2C%20JosephsonCircuits.ModeLayout%2C%20JosephsonCircuits.ModeLayout%7D%2C%20Tuple%7BSparseMatrixCSC%7BComplex%7BT%7D%2C%20Ti%7D%2C%20JosephsonCircuits.ModeLayout%2C%20JosephsonCircuits.ModeLayout%2C%20Type%7BTj%7D%7D%7D%20where%20%7BT%3C%3AReal%2C%20Ti%2C%20Tj%3C%3AInteger%7D)

```@raw html
<a id="JosephsonCircuits.complexdim-Tuple{Integer, AbstractVector{Bool}}"></a>
```

[`JosephsonCircuits.complexdim`](api/internals.md#JosephsonCircuits.complexdim-Tuple%7BInteger%2C%20AbstractVector%7BBool%7D%7D)

```@raw html
<a id="JosephsonCircuits.complextorealkernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.complextorealkernel!`](api/internals.md#JosephsonCircuits.complextorealkernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.componentindex-Tuple{JosephsonCircuits.CompiledCircuit, Any}"></a>
```

[`JosephsonCircuits.componentindex`](api/internals.md#JosephsonCircuits.componentindex-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.componentlookups-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.componentlookups`](api/internals.md#JosephsonCircuits.componentlookups-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.componentnports-Tuple{ScatteringParameters}"></a>
```

[`JosephsonCircuits.componentnports`](api/internals.md#JosephsonCircuits.componentnports-Tuple%7BScatteringParameters%7D)

```@raw html
<a id="JosephsonCircuits.componentstamp-Tuple{Integer, JosephsonCircuits.CompiledCircuit, JosephsonCircuits.CircuitMatrices, Any, Integer}"></a>
```

[`JosephsonCircuits.componentstamp`](api/internals.md#JosephsonCircuits.componentstamp-Tuple%7BInteger%2C%20JosephsonCircuits.CompiledCircuit%2C%20JosephsonCircuits.CircuitMatrices%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.componentvalues-Tuple{HBCache, NamedTuple}"></a>
```

[`JosephsonCircuits.componentvalues`](api/internals.md#JosephsonCircuits.componentvalues-Tuple%7BHBCache%2C%20NamedTuple%7D)

```@raw html
<a id="JosephsonCircuits.componentvaluestonumber-Tuple{Vector, Dict}"></a>
```

[`JosephsonCircuits.componentvaluestonumber`](api/internals.md#JosephsonCircuits.componentvaluestonumber-Tuple%7BVector%2C%20Dict%7D)

```@raw html
<a id="JosephsonCircuits.compositelayout-Tuple{JosephsonCircuits.ModeLayout, AbstractVector{&lt;:Tuple}}"></a>
```

[`JosephsonCircuits.compositelayout`](api/internals.md#JosephsonCircuits.compositelayout-Tuple%7BJosephsonCircuits.ModeLayout%2C%20AbstractVector%7B%3C%3ATuple%7D%7D)

```@raw html
<a id="JosephsonCircuits.compositelayout-Tuple{JosephsonCircuits.ModeLayout, AbstractVector{Bool}}"></a>
```

[`JosephsonCircuits.compositelayout`](api/internals.md#JosephsonCircuits.compositelayout-Tuple%7BJosephsonCircuits.ModeLayout%2C%20AbstractVector%7BBool%7D%7D)

```@raw html
<a id="JosephsonCircuits.conjnegfreq!-Tuple{SparseMatrixCSC, Vector}"></a>
```

[`JosephsonCircuits.conjnegfreq!`](api/internals.md#JosephsonCircuits.conjnegfreq!-Tuple%7BSparseMatrixCSC%2C%20Vector%7D)

```@raw html
<a id="JosephsonCircuits.conjsym-Union{Tuple{N}, Tuple{NTuple{N, Int64}, NTuple{N, Int64}}} where N"></a>
```

[`JosephsonCircuits.conjsym`](api/internals.md#JosephsonCircuits.conjsym-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.conjugateladder-Tuple{LinearizedScattering}"></a>
```

[`JosephsonCircuits.conjugateladder`](api/internals.md#JosephsonCircuits.conjugateladder-Tuple%7BLinearizedScattering%7D)

```@raw html
<a id="JosephsonCircuits.conjugatemultiplicity-Tuple{AbstractArray, AbstractArray}"></a>
```

[`JosephsonCircuits.conjugatemultiplicity`](api/internals.md#JosephsonCircuits.conjugatemultiplicity-Tuple%7BAbstractArray%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.connectS!-Union{Tuple{N}, Tuple{T}, Tuple{Graphs.SimpleGraphs.SimpleDiGraph{Int64}, AbstractVector{&lt;:AbstractArray{Tuple{T, T, Int64, Int64}, 1}}, AbstractVector{&lt;:AbstractVector{Int64}}, AbstractVector{&lt;:AbstractArray{Tuple{T, Int64}, 1}}, AbstractVector{N}, AbstractVector{N}}} where {T, N}"></a>
```

[`JosephsonCircuits.connectS!`](api/internals.md#JosephsonCircuits.connectS!-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BGraphs.SimpleGraphs.SimpleDiGraph%7BInt64%7D%2C%20AbstractVector%7B%3C%3AAbstractArray%7BTuple%7BT%2C%20T%2C%20Int64%2C%20Int64%7D%2C%201%7D%7D%2C%20AbstractVector%7B%3C%3AAbstractVector%7BInt64%7D%7D%2C%20AbstractVector%7B%3C%3AAbstractArray%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%7D%2C%20AbstractVector%7BN%7D%2C%20AbstractVector%7BN%7D%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.connectS_initialize-Tuple{AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.connectS_initialize`](api/internals.md#JosephsonCircuits.connectS_initialize-Tuple%7BAbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.connectS_initialize-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{JosephsonCircuits.PassiveNetwork{T, N}, 1}, AbstractArray{Tuple{T, T, Int64, Int64}, 1}}} where {T, N}"></a>
```

[`JosephsonCircuits.connectS_initialize`](api/internals.md#JosephsonCircuits.connectS_initialize-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BJosephsonCircuits.PassiveNetwork%7BT%2C%20N%7D%2C%201%7D%2C%20AbstractArray%7BTuple%7BT%2C%20T%2C%20Int64%2C%20Int64%7D%2C%201%7D%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.conversiontype-Tuple{AbstractArray, Vararg{Any}}"></a>
```

[`JosephsonCircuits.conversiontype`](api/internals.md#JosephsonCircuits.conversiontype-Tuple%7BAbstractArray%2C%20Vararg%7BAny%7D%7D)

```@raw html
<a id="JosephsonCircuits.convertcopy-Union{Tuple{F}, Tuple{F, AbstractArray, Any}} where F"></a>
```

[`JosephsonCircuits.convertcopy`](api/internals.md#JosephsonCircuits.convertcopy-Union%7BTuple%7BF%7D%2C%20Tuple%7BF%2C%20AbstractArray%2C%20Any%7D%7D%20where%20F)

```@raw html
<a id="JosephsonCircuits.convertperfrequency!-Union{Tuple{F}, Tuple{F, AbstractArray, AbstractArray}} where F"></a>
```

[`JosephsonCircuits.convertperfrequency!`](api/internals.md#JosephsonCircuits.convertperfrequency!-Union%7BTuple%7BF%7D%2C%20Tuple%7BF%2C%20AbstractArray%2C%20AbstractArray%7D%7D%20where%20F)

```@raw html
<a id="JosephsonCircuits.cosdirectionalderivative!-Tuple{Array, JosephsonCircuits.HBSystem, AbstractVector{&lt;:Complex}}"></a>
```

[`JosephsonCircuits.cosdirectionalderivative!`](api/internals.md#JosephsonCircuits.cosdirectionalderivative!-Tuple%7BArray%2C%20JosephsonCircuits.HBSystem%2C%20AbstractVector%7B%3C%3AComplex%7D%7D)

```@raw html
<a id="JosephsonCircuits.cosphibandwidths-Tuple{Any, Matrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.cosphibandwidths`](api/internals.md#JosephsonCircuits.cosphibandwidths-Tuple%7BAny%2C%20Matrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.cosphimatrix-Tuple{JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.cosphimatrix`](api/internals.md#JosephsonCircuits.cosphimatrix-Tuple%7BJosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.countscatteringports-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.countscatteringports`](api/internals.md#JosephsonCircuits.countscatteringports-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.coupledinductornames-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.coupledinductornames`](api/internals.md#JosephsonCircuits.coupledinductornames-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.coupling_to_even_odd-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.coupling_to_even_odd`](api/internals.md#JosephsonCircuits.coupling_to_even_odd-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.couplingbytes-Tuple{JosephsonCircuits.PreconditionerPlan, JosephsonCircuits.AbstractModeCoupling, Any, Any}"></a>
```

[`JosephsonCircuits.couplingbytes`](api/internals.md#JosephsonCircuits.couplingbytes-Tuple%7BJosephsonCircuits.PreconditionerPlan%2C%20JosephsonCircuits.AbstractModeCoupling%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.couplingmask-Tuple{BlockDiagonal, Integer, Any}"></a>
```

[`JosephsonCircuits.couplingmask`](api/internals.md#JosephsonCircuits.couplingmask-Tuple%7BBlockDiagonal%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.cprderivative"></a>
```

[`JosephsonCircuits.cprderivative`](api/internals.md#JosephsonCircuits.cprderivative)

```@raw html
<a id="JosephsonCircuits.cscvaluepermutation-Tuple{SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.cscvaluepermutation`](api/internals.md#JosephsonCircuits.cscvaluepermutation-Tuple%7BSparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.cubic_trial_step-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.cubic_trial_step`](api/internals.md#JosephsonCircuits.cubic_trial_step-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.dcblockdescriptor-Tuple{JosephsonCircuits.StampedScatteringBlock}"></a>
```

[`JosephsonCircuits.dcblockdescriptor`](api/internals.md#JosephsonCircuits.dcblockdescriptor-Tuple%7BJosephsonCircuits.StampedScatteringBlock%7D)

```@raw html
<a id="JosephsonCircuits.dcblockrows-Tuple{AbstractVector, Vector{Int64}, Integer, Integer, Integer, Any}"></a>
```

[`JosephsonCircuits.dcblockrows`](api/internals.md#JosephsonCircuits.dcblockrows-Tuple%7BAbstractVector%2C%20Vector%7BInt64%7D%2C%20Integer%2C%20Integer%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.dcconductanceplan-Tuple{Vector{Vector{Int64}}, SparseMatrixCSC, AbstractVector, Integer, Integer}"></a>
```

[`JosephsonCircuits.dcconductanceplan`](api/internals.md#JosephsonCircuits.dcconductanceplan-Tuple%7BVector%7BVector%7BInt64%7D%7D%2C%20SparseMatrixCSC%2C%20AbstractVector%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.dccoupling-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dccoupling`](api/internals.md#JosephsonCircuits.dccoupling-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dcinjected-Tuple{JosephsonCircuits.DCConductancePlan, AbstractVector, Integer}"></a>
```

[`JosephsonCircuits.dcinjected`](api/internals.md#JosephsonCircuits.dcinjected-Tuple%7BJosephsonCircuits.DCConductancePlan%2C%20AbstractVector%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.dckeep-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dckeep`](api/internals.md#JosephsonCircuits.dckeep-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dclimit-Tuple{JosephsonCircuits.StampedScatteringBlock, Integer, Real}"></a>
```

[`JosephsonCircuits.dclimit`](api/internals.md#JosephsonCircuits.dclimit-Tuple%7BJosephsonCircuits.StampedScatteringBlock%2C%20Integer%2C%20Real%7D)

```@raw html
<a id="JosephsonCircuits.dcpinning-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dcpinning`](api/internals.md#JosephsonCircuits.dcpinning-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dcresidualsensitivity"></a>
```

[`JosephsonCircuits.dcresidualsensitivity`](api/internals.md#JosephsonCircuits.dcresidualsensitivity)

```@raw html
<a id="JosephsonCircuits.dcscatteringmatrix-Tuple{OpenDC, Integer}"></a>
```

[`JosephsonCircuits.dcscatteringmatrix`](api/internals.md#JosephsonCircuits.dcscatteringmatrix-Tuple%7BOpenDC%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.dcsolutionfrom-Tuple{JosephsonCircuits.DCConductancePlan, AbstractVector}"></a>
```

[`JosephsonCircuits.dcsolutionfrom`](api/internals.md#JosephsonCircuits.dcsolutionfrom-Tuple%7BJosephsonCircuits.DCConductancePlan%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.dcsourcecurrent-Tuple{JosephsonCircuits.DCConductancePlan, AbstractVector, Integer}"></a>
```

[`JosephsonCircuits.dcsourcecurrent`](api/internals.md#JosephsonCircuits.dcsourcecurrent-Tuple%7BJosephsonCircuits.DCConductancePlan%2C%20AbstractVector%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.dcsubsystem-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dcsubsystem`](api/internals.md#JosephsonCircuits.dcsubsystem-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dcsubsystemindices-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dcsubsystemindices`](api/internals.md#JosephsonCircuits.dcsubsystemindices-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dcsubsystemlocal-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dcsubsystemlocal`](api/internals.md#JosephsonCircuits.dcsubsystemlocal-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dcsubsystemrhs-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dcsubsystemrhs`](api/internals.md#JosephsonCircuits.dcsubsystemrhs-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dcupdate-Tuple{JosephsonCircuits.CanonicalWork}"></a>
```

[`JosephsonCircuits.dcupdate`](api/internals.md#JosephsonCircuits.dcupdate-Tuple%7BJosephsonCircuits.CanonicalWork%7D)

```@raw html
<a id="JosephsonCircuits.dcvoltagesensitivity-Tuple{JosephsonCircuits.HBOperatingPoint, AbstractMatrix}"></a>
```

[`JosephsonCircuits.dcvoltagesensitivity`](api/internals.md#JosephsonCircuits.dcvoltagesensitivity-Tuple%7BJosephsonCircuits.HBOperatingPoint%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.defaultgridladder-Union{Tuple{NTuple{N, Int64}}, Tuple{N}} where N"></a>
```

[`JosephsonCircuits.defaultgridladder`](api/internals.md#JosephsonCircuits.defaultgridladder-Union%7BTuple%7BNTuple%7BN%2C%20Int64%7D%7D%2C%20Tuple%7BN%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.definitionsbyname-Tuple{Dict{Symbol, Any}}"></a>
```

[`JosephsonCircuits.definitionsbyname`](api/internals.md#JosephsonCircuits.definitionsbyname-Tuple%7BDict%7BSymbol%2C%20Any%7D%7D)

```@raw html
<a id="JosephsonCircuits.definitiontable-Tuple{Dict{Any, Any}}"></a>
```

[`JosephsonCircuits.definitiontable`](api/internals.md#JosephsonCircuits.definitiontable-Tuple%7BDict%7BAny%2C%20Any%7D%7D)

```@raw html
<a id="JosephsonCircuits.deflationproducts-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.deflationproducts`](api/internals.md#JosephsonCircuits.deflationproducts-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.deflationrebuilds-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.deflationrebuilds`](api/internals.md#JosephsonCircuits.deflationrebuilds-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.deflationsize-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.deflationsize`](api/internals.md#JosephsonCircuits.deflationsize-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.deprecatedsolverkeywords-Tuple{Symbol, Any}"></a>
```

[`JosephsonCircuits.deprecatedsolverkeywords`](api/internals.md#JosephsonCircuits.deprecatedsolverkeywords-Tuple%7BSymbol%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.derivativeat-Tuple{JosephsonCircuits.JunctionRelations, Any}"></a>
```

[`JosephsonCircuits.derivativeat`](api/internals.md#JosephsonCircuits.derivativeat-Tuple%7BJosephsonCircuits.JunctionRelations%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.derivativeinto!-Tuple{Any, JosephsonCircuits.JunctionRelations, Any}"></a>
```

[`JosephsonCircuits.derivativeinto!`](api/internals.md#JosephsonCircuits.derivativeinto!-Tuple%7BAny%2C%20JosephsonCircuits.JunctionRelations%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.designblockjacobian-Tuple{Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, Any}"></a>
```

[`JosephsonCircuits.designblockjacobian`](api/internals.md#JosephsonCircuits.designblockjacobian-Tuple%7BUnion%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.designderivative-Tuple{Number, Symbol, Any}"></a>
```

[`JosephsonCircuits.designderivative`](api/internals.md#JosephsonCircuits.designderivative-Tuple%7BNumber%2C%20Symbol%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.devicecomplexjacobianpattern-Tuple{Integer, Integer, Vararg{Any, 4}}"></a>
```

[`JosephsonCircuits.devicecomplexjacobianpattern`](api/internals.md#JosephsonCircuits.devicecomplexjacobianpattern-Tuple%7BInteger%2C%20Integer%2C%20Vararg%7BAny%2C%204%7D%7D)

```@raw html
<a id="JosephsonCircuits.deviceexpandrealpattern-Tuple{Any, Any, JosephsonCircuits.ModeLayout, Integer, Any}"></a>
```

[`JosephsonCircuits.deviceexpandrealpattern`](api/internals.md#JosephsonCircuits.deviceexpandrealpattern-Tuple%7BAny%2C%20Any%2C%20JosephsonCircuits.ModeLayout%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.devicenoise"></a>
```

[`JosephsonCircuits.devicenoise`](api/internals.md#JosephsonCircuits.devicenoise)

```@raw html
<a id="JosephsonCircuits.devicesolutions"></a>
```

[`JosephsonCircuits.devicesolutions`](api/internals.md#JosephsonCircuits.devicesolutions)

```@raw html
<a id="JosephsonCircuits.diagrepeat-Tuple{SparseMatrixCSC, Integer}"></a>
```

[`JosephsonCircuits.diagrepeat`](api/internals.md#JosephsonCircuits.diagrepeat-Tuple%7BSparseMatrixCSC%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.diagrepeat-Tuple{SparseVector, Integer}"></a>
```

[`JosephsonCircuits.diagrepeat`](api/internals.md#JosephsonCircuits.diagrepeat-Tuple%7BSparseVector%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.direct_sum-Tuple{Any}"></a>
```

[`JosephsonCircuits.direct_sum`](api/internals.md#JosephsonCircuits.direct_sum-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.directcoefficients-Tuple{JosephsonCircuits.TransientSystem, Symbol}"></a>
```

[`JosephsonCircuits.directcoefficients`](api/internals.md#JosephsonCircuits.directcoefficients-Tuple%7BJosephsonCircuits.TransientSystem%2C%20Symbol%7D)

```@raw html
<a id="JosephsonCircuits.dualsearch!-Tuple{Any, AbstractVector, AbstractVector, AbstractVector, AbstractVector, Real, Real, Real, AbstractVector, Real, AbstractVector, AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.dualsearch!`](api/internals.md#JosephsonCircuits.dualsearch!-Tuple%7BAny%2C%20AbstractVector%2C%20AbstractVector%2C%20AbstractVector%2C%20AbstractVector%2C%20Real%2C%20Real%2C%20Real%2C%20AbstractVector%2C%20Real%2C%20AbstractVector%2C%20AbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.eliminationtree-Tuple{Any, AbstractVector{&lt;:Integer}}"></a>
```

[`JosephsonCircuits.eliminationtree`](api/internals.md#JosephsonCircuits.eliminationtree-Tuple%7BAny%2C%20AbstractVector%7B%3C%3AInteger%7D%7D)

```@raw html
<a id="JosephsonCircuits.emptyrelations-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.emptyrelations`](api/internals.md#JosephsonCircuits.emptyrelations-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.ensurecolumns!-Tuple{JosephsonCircuits.GMRESWorkspace, Integer}"></a>
```

[`JosephsonCircuits.ensurecolumns!`](api/internals.md#JosephsonCircuits.ensurecolumns!-Tuple%7BJosephsonCircuits.GMRESWorkspace%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.equilibrate!-Tuple{AbstractMatrix, Any, Any, JosephsonCircuits.SweepEquilibration, Any}"></a>
```

[`JosephsonCircuits.equilibrate!`](api/internals.md#JosephsonCircuits.equilibrate!-Tuple%7BAbstractMatrix%2C%20Any%2C%20Any%2C%20JosephsonCircuits.SweepEquilibration%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.erased-Tuple{Function}"></a>
```

[`JosephsonCircuits.erased`](api/internals.md#JosephsonCircuits.erased-Tuple%7BFunction%7D)

```@raw html
<a id="JosephsonCircuits.escalatepreconditioner!-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.escalatepreconditioner!`](api/internals.md#JosephsonCircuits.escalatepreconditioner!-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.escalatepreconditioner!-Tuple{JosephsonCircuits.ModeCouplingPreconditioner}"></a>
```

[`JosephsonCircuits.escalatepreconditioner!`](api/internals.md#JosephsonCircuits.escalatepreconditioner!-Tuple%7BJosephsonCircuits.ModeCouplingPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.evaluatecovariance!"></a>
```

[`JosephsonCircuits.evaluatecovariance!`](api/internals.md#JosephsonCircuits.evaluatecovariance!)

```@raw html
<a id="JosephsonCircuits.evaluateharmonics!-Tuple{AbstractArray{ComplexF64, 4}, LinearizedScattering, AbstractVector}"></a>
```

[`JosephsonCircuits.evaluateharmonics!`](api/internals.md#JosephsonCircuits.evaluateharmonics!-Tuple%7BAbstractArray%7BComplexF64%2C%204%7D%2C%20LinearizedScattering%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.evaluatehybrid!-Tuple{AbstractArray{ComplexF64, 3}, AbstractArray{ComplexF64, 3}, ScatteringParameters, AbstractVector, JosephsonCircuits.HybridWorkspace}"></a>
```

[`JosephsonCircuits.evaluatehybrid!`](api/internals.md#JosephsonCircuits.evaluatehybrid!-Tuple%7BAbstractArray%7BComplexF64%2C%203%7D%2C%20AbstractArray%7BComplexF64%2C%203%7D%2C%20ScatteringParameters%2C%20AbstractVector%2C%20JosephsonCircuits.HybridWorkspace%7D)

```@raw html
<a id="JosephsonCircuits.evaluatehybridpumped!-Tuple{AbstractArray{ComplexF64, 4}, AbstractArray{ComplexF64, 4}, LinearizedScattering, AbstractVector, AbstractMatrix{Int64}, JosephsonCircuits.ScatteringWorkspace}"></a>
```

[`JosephsonCircuits.evaluatehybridpumped!`](api/internals.md#JosephsonCircuits.evaluatehybridpumped!-Tuple%7BAbstractArray%7BComplexF64%2C%204%7D%2C%20AbstractArray%7BComplexF64%2C%204%7D%2C%20LinearizedScattering%2C%20AbstractVector%2C%20AbstractMatrix%7BInt64%7D%2C%20JosephsonCircuits.ScatteringWorkspace%7D)

```@raw html
<a id="JosephsonCircuits.evaluateprovider!-Union{Tuple{T}, Tuple{AbstractArray{T, 3}, JosephsonCircuits.ConstantMatrixProvider, AbstractVector}} where T"></a>
```

[`JosephsonCircuits.evaluateprovider!`](api/internals.md#JosephsonCircuits.evaluateprovider!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%203%7D%2C%20JosephsonCircuits.ConstantMatrixProvider%2C%20AbstractVector%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.evaluatescattering!"></a>
```

[`JosephsonCircuits.evaluatescattering!`](api/internals.md#JosephsonCircuits.evaluatescattering!)

```@raw html
<a id="JosephsonCircuits.even_odd_to_coupling-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.even_odd_to_coupling`](api/internals.md#JosephsonCircuits.even_odd_to_coupling-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.even_odd_to_maxwell-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.even_odd_to_maxwell`](api/internals.md#JosephsonCircuits.even_odd_to_maxwell-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.even_odd_to_mutual-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.even_odd_to_mutual`](api/internals.md#JosephsonCircuits.even_odd_to_mutual-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.export_netlist!-Tuple{IO, AbstractVector, Dict}"></a>
```

[`JosephsonCircuits.export_netlist!`](api/internals.md#JosephsonCircuits.export_netlist!-Tuple%7BIO%2C%20AbstractVector%2C%20Dict%7D)

```@raw html
<a id="JosephsonCircuits.export_netlist-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.export_netlist`](api/internals.md#JosephsonCircuits.export_netlist-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.export_netlist-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.export_netlist`](api/internals.md#JosephsonCircuits.export_netlist-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.exportnetlist-Tuple{Union{JosephsonCircuits.CompiledCircuit, ElaboratedCircuit, Circuit}, Dict}"></a>
```

[`JosephsonCircuits.exportnetlist`](api/internals.md#JosephsonCircuits.exportnetlist-Tuple%7BUnion%7BJosephsonCircuits.CompiledCircuit%2C%20ElaboratedCircuit%2C%20Circuit%7D%2C%20Dict%7D)

```@raw html
<a id="JosephsonCircuits.extractbranches-Tuple{Vector{Symbol}, Matrix{Int64}}"></a>
```

[`JosephsonCircuits.extractbranches`](api/internals.md#JosephsonCircuits.extractbranches-Tuple%7BVector%7BSymbol%7D%2C%20Matrix%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.factorizationprecision-Tuple{JosephsonCircuits.AbstractFactorization}"></a>
```

[`JosephsonCircuits.factorizationprecision`](api/internals.md#JosephsonCircuits.factorizationprecision-Tuple%7BJosephsonCircuits.AbstractFactorization%7D)

```@raw html
<a id="JosephsonCircuits.factorize-Tuple{BlockFactorization, SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.factorize`](api/internals.md#JosephsonCircuits.factorize-Tuple%7BBlockFactorization%2C%20SparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.fftplans-Union{Tuple{T}, Tuple{AbstractArray{Complex{T}}, AbstractArray{T}, Int64, KernelAbstractions.CPU}} where T"></a>
```

[`JosephsonCircuits.fftplans`](api/internals.md#JosephsonCircuits.fftplans-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BComplex%7BT%7D%7D%2C%20AbstractArray%7BT%7D%2C%20Int64%2C%20KernelAbstractions.CPU%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.fillandfactorize!-Union{Tuple{T}, Tuple{JosephsonCircuits.SparseBlockFactorization{T}, AbstractMatrix}} where T"></a>
```

[`JosephsonCircuits.fillandfactorize!`](api/internals.md#JosephsonCircuits.fillandfactorize!-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.SparseBlockFactorization%7BT%7D%2C%20AbstractMatrix%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.fillordering-Tuple{JosephsonCircuits.AbstractFactorization, Any}"></a>
```

[`JosephsonCircuits.fillordering`](api/internals.md#JosephsonCircuits.fillordering-Tuple%7BJosephsonCircuits.AbstractFactorization%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.find_duplicate_connections-Union{Tuple{AbstractArray{Tuple{T, T, Int64, Int64}, 1}}, Tuple{T}} where T"></a>
```

[`JosephsonCircuits.find_duplicate_connections`](api/internals.md#JosephsonCircuits.find_duplicate_connections-Union%7BTuple%7BAbstractArray%7BTuple%7BT%2C%20T%2C%20Int64%2C%20Int64%7D%2C%201%7D%7D%2C%20Tuple%7BT%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.find_duplicate_network_names-Union{Tuple{AbstractArray{JosephsonCircuits.PassiveNetwork{T, N}, 1}}, Tuple{N}, Tuple{T}} where {T, N}"></a>
```

[`JosephsonCircuits.find_duplicate_network_names`](api/internals.md#JosephsonCircuits.find_duplicate_network_names-Union%7BTuple%7BAbstractArray%7BJosephsonCircuits.PassiveNetwork%7BT%2C%20N%7D%2C%201%7D%7D%2C%20Tuple%7BN%7D%2C%20Tuple%7BT%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.findgroundnodeindex-Tuple{Vector{String}}"></a>
```

[`JosephsonCircuits.findgroundnodeindex`](api/internals.md#JosephsonCircuits.findgroundnodeindex-Tuple%7BVector%7BString%7D%7D)

```@raw html
<a id="JosephsonCircuits.foreachbatch-Union{Tuple{F}, Tuple{F, Any, Integer}} where F"></a>
```

[`JosephsonCircuits.foreachbatch`](api/internals.md#JosephsonCircuits.foreachbatch-Union%7BTuple%7BF%7D%2C%20Tuple%7BF%2C%20Any%2C%20Integer%7D%7D%20where%20F)

```@raw html
<a id="JosephsonCircuits.forwardsolution!-Tuple{Any, JosephsonCircuits.DeviceSweep, Integer}"></a>
```

[`JosephsonCircuits.forwardsolution!`](api/internals.md#JosephsonCircuits.forwardsolution!-Tuple%7BAny%2C%20JosephsonCircuits.DeviceSweep%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.forwardtermkernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.forwardtermkernel!`](api/internals.md#JosephsonCircuits.forwardtermkernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.forwardtermkernelcomplex!-Tuple{Any}"></a>
```

[`JosephsonCircuits.forwardtermkernelcomplex!`](api/internals.md#JosephsonCircuits.forwardtermkernelcomplex!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.forwardtransposekernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.forwardtransposekernel!`](api/internals.md#JosephsonCircuits.forwardtransposekernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.fourierindices-Tuple{JosephsonCircuits.Frequencies}"></a>
```

[`JosephsonCircuits.fourierindices`](api/internals.md#JosephsonCircuits.fourierindices-Tuple%7BJosephsonCircuits.Frequencies%7D)

```@raw html
<a id="JosephsonCircuits.freememory-Tuple{KernelAbstractions.CPU}"></a>
```

[`JosephsonCircuits.freememory`](api/internals.md#JosephsonCircuits.freememory-Tuple%7BKernelAbstractions.CPU%7D)

```@raw html
<a id="JosephsonCircuits.freqsubst-Tuple{SparseMatrixCSC, Vector}"></a>
```

[`JosephsonCircuits.freqsubst`](api/internals.md#JosephsonCircuits.freqsubst-Tuple%7BSparseMatrixCSC%2C%20Vector%7D)

```@raw html
<a id="JosephsonCircuits.frequencydependentcircuit-Tuple{JosephsonCircuits.CompiledCircuit, Any, Any, Symbol}"></a>
```

[`JosephsonCircuits.frequencydependentcircuit`](api/internals.md#JosephsonCircuits.frequencydependentcircuit-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Any%2C%20Any%2C%20Symbol%7D)

```@raw html
<a id="JosephsonCircuits.gathercanonical!-Tuple{AbstractVector, AbstractVector, JosephsonCircuits.CompositeLayout}"></a>
```

[`JosephsonCircuits.gathercanonical!`](api/internals.md#JosephsonCircuits.gathercanonical!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20JosephsonCircuits.CompositeLayout%7D)

```@raw html
<a id="JosephsonCircuits.gatherportrows!-Tuple{AbstractArray{&lt;:Any, 3}, AbstractArray{&lt;:Any, 3}, AbstractVector, Any}"></a>
```

[`JosephsonCircuits.gatherportrows!`](api/internals.md#JosephsonCircuits.gatherportrows!-Tuple%7BAbstractArray%7B%3C%3AAny%2C%203%7D%2C%20AbstractArray%7B%3C%3AAny%2C%203%7D%2C%20AbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.gathervalues!-Tuple{AbstractArray, AbstractVector, AbstractArray}"></a>
```

[`JosephsonCircuits.gathervalues!`](api/internals.md#JosephsonCircuits.gathervalues!-Tuple%7BAbstractArray%2C%20AbstractVector%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.get_ports-Union{Tuple{Tuple{T, N, Array{Tuple{T, Int64}, 1}}}, Tuple{N}, Tuple{T}} where {T, N}"></a>
```

[`JosephsonCircuits.get_ports`](api/internals.md#JosephsonCircuits.get_ports-Union%7BTuple%7BTuple%7BT%2C%20N%2C%20Array%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%7D%7D%2C%20Tuple%7BN%7D%2C%20Tuple%7BT%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.get_ports-Union{Tuple{Union{Tuple{T, N}, Tuple{T, N, N}}}, Tuple{N}, Tuple{T}} where {T, N}"></a>
```

[`JosephsonCircuits.get_ports`](api/internals.md#JosephsonCircuits.get_ports-Union%7BTuple%7BUnion%7BTuple%7BT%2C%20N%7D%2C%20Tuple%7BT%2C%20N%2C%20N%7D%7D%7D%2C%20Tuple%7BN%7D%2C%20Tuple%7BT%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.gmres!-Union{Tuple{T}, Tuple{AbstractVector{T}, Any, AbstractVector{T}, JosephsonCircuits.GMRESWorkspace{T, TV, TM} where {TV&lt;:AbstractVector{T}, TM&lt;:AbstractMatrix{T}}}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.gmres!`](api/internals.md#JosephsonCircuits.gmres!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractVector%7BT%7D%2C%20Any%2C%20AbstractVector%7BT%7D%2C%20JosephsonCircuits.GMRESWorkspace%7BT%2C%20TV%2C%20TM%7D%20where%20%7BTV%3C%3AAbstractVector%7BT%7D%2C%20TM%3C%3AAbstractMatrix%7BT%7D%7D%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.gmres_applyrotations!-Union{Tuple{T}, Tuple{AbstractMatrix{T}, AbstractVector{T}, AbstractVector{T}, AbstractVector{T}, Integer}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.gmres_applyrotations!`](api/internals.md#JosephsonCircuits.gmres_applyrotations!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractMatrix%7BT%7D%2C%20AbstractVector%7BT%7D%2C%20AbstractVector%7BT%7D%2C%20AbstractVector%7BT%7D%2C%20Integer%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.gmres_correction!-Union{Tuple{T}, Tuple{AbstractVector{T}, JosephsonCircuits.GMRESWorkspace{T, TV, TM} where {TV&lt;:AbstractVector{T}, TM&lt;:AbstractMatrix{T}}, Integer, Any}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.gmres_correction!`](api/internals.md#JosephsonCircuits.gmres_correction!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractVector%7BT%7D%2C%20JosephsonCircuits.GMRESWorkspace%7BT%2C%20TV%2C%20TM%7D%20where%20%7BTV%3C%3AAbstractVector%7BT%7D%2C%20TM%3C%3AAbstractMatrix%7BT%7D%7D%2C%20Integer%2C%20Any%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.gmres_givens-Union{Tuple{T}, Tuple{T, T}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.gmres_givens`](api/internals.md#JosephsonCircuits.gmres_givens-Union%7BTuple%7BT%7D%2C%20Tuple%7BT%2C%20T%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.gmres_orthogonalize!-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractMatrix{T}, AbstractMatrix{T}, AbstractVector{T}, AbstractVector{T}, Integer}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.gmres_orthogonalize!`](api/internals.md#JosephsonCircuits.gmres_orthogonalize!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractVector%7BT%7D%2C%20AbstractMatrix%7BT%7D%2C%20AbstractMatrix%7BT%7D%2C%20AbstractVector%7BT%7D%2C%20AbstractVector%7BT%7D%2C%20Integer%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.groupdestinations-Tuple{AbstractVector{Int64}}"></a>
```

[`JosephsonCircuits.groupdestinations`](api/internals.md#JosephsonCircuits.groupdestinations-Tuple%7BAbstractVector%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.halmos_dilation-Tuple{Any}"></a>
```

[`JosephsonCircuits.halmos_dilation`](api/internals.md#JosephsonCircuits.halmos_dilation-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.harmonicritznearzero-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Integer}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.harmonicritznearzero`](api/internals.md#JosephsonCircuits.harmonicritznearzero-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractMatrix%7BT%7D%2C%20Integer%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.harvest!-Tuple{JosephsonCircuits.AbstractPreconditioner, JosephsonCircuits.GMRESWorkspace, NamedTuple}"></a>
```

[`JosephsonCircuits.harvest!`](api/internals.md#JosephsonCircuits.harvest!-Tuple%7BJosephsonCircuits.AbstractPreconditioner%2C%20JosephsonCircuits.GMRESWorkspace%2C%20NamedTuple%7D)

```@raw html
<a id="JosephsonCircuits.harvestcycle!-Tuple{JosephsonCircuits.AbstractPreconditioner, JosephsonCircuits.GMRESWorkspace, Integer}"></a>
```

[`JosephsonCircuits.harvestcycle!`](api/internals.md#JosephsonCircuits.harvestcycle!-Tuple%7BJosephsonCircuits.AbstractPreconditioner%2C%20JosephsonCircuits.GMRESWorkspace%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.harvestdimension-Tuple{JosephsonCircuits.GMRESWorkspace, NamedTuple}"></a>
```

[`JosephsonCircuits.harvestdimension`](api/internals.md#JosephsonCircuits.harvestdimension-Tuple%7BJosephsonCircuits.GMRESWorkspace%2C%20NamedTuple%7D)

```@raw html
<a id="JosephsonCircuits.hasports-Tuple{Union{GaussianChannel, LinearizedScattering, ScatteringParameters}}"></a>
```

[`JosephsonCircuits.hasports`](api/internals.md#JosephsonCircuits.hasports-Tuple%7BUnion%7BGaussianChannel%2C%20LinearizedScattering%2C%20ScatteringParameters%7D%7D)

```@raw html
<a id="JosephsonCircuits.hasrealbackward-Tuple{JosephsonCircuits.NonlinearTermPlan}"></a>
```

[`JosephsonCircuits.hasrealbackward`](api/internals.md#JosephsonCircuits.hasrealbackward-Tuple%7BJosephsonCircuits.NonlinearTermPlan%7D)

```@raw html
<a id="JosephsonCircuits.hbconjmatind-Union{Tuple{JosephsonCircuits.Frequencies{N}}, Tuple{N}} where N"></a>
```

[`JosephsonCircuits.hbconjmatind`](api/internals.md#JosephsonCircuits.hbconjmatind-Union%7BTuple%7BJosephsonCircuits.Frequencies%7BN%7D%7D%2C%20Tuple%7BN%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hbconjmatind-Union{Tuple{N}, Tuple{JosephsonCircuits.Frequencies{N}, JosephsonCircuits.Frequencies{N}}} where N"></a>
```

[`JosephsonCircuits.hbconjmatind`](api/internals.md#JosephsonCircuits.hbconjmatind-Union%7BTuple%7BN%7D%2C%20Tuple%7BJosephsonCircuits.Frequencies%7BN%7D%2C%20JosephsonCircuits.Frequencies%7BN%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hblinearsolve!-Tuple{GMRES, Vararg{Any, 5}}"></a>
```

[`JosephsonCircuits.hblinearsolve!`](api/internals.md#JosephsonCircuits.hblinearsolve!-Tuple%7BGMRES%2C%20Vararg%7BAny%2C%205%7D%7D)

```@raw html
<a id="JosephsonCircuits.hblinsolve_inner!-Tuple{JosephsonCircuits.LinearizedWorkspace, JosephsonCircuits.LinearizedArrays, Vararg{Any, 14}}"></a>
```

[`JosephsonCircuits.hblinsolve_inner!`](api/internals.md#JosephsonCircuits.hblinsolve_inner!-Tuple%7BJosephsonCircuits.LinearizedWorkspace%2C%20JosephsonCircuits.LinearizedArrays%2C%20Vararg%7BAny%2C%2014%7D%7D)

```@raw html
<a id="JosephsonCircuits.hbmatind-Union{Tuple{JosephsonCircuits.Frequencies{N}}, Tuple{N}} where N"></a>
```

[`JosephsonCircuits.hbmatind`](api/internals.md#JosephsonCircuits.hbmatind-Union%7BTuple%7BJosephsonCircuits.Frequencies%7BN%7D%7D%2C%20Tuple%7BN%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hbmatind-Union{Tuple{N}, Tuple{JosephsonCircuits.Frequencies{N}, JosephsonCircuits.Frequencies{N}}} where N"></a>
```

[`JosephsonCircuits.hbmatind`](api/internals.md#JosephsonCircuits.hbmatind-Union%7BTuple%7BN%7D%2C%20Tuple%7BJosephsonCircuits.Frequencies%7BN%7D%2C%20JosephsonCircuits.Frequencies%7BN%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hbmatindices-Union{Tuple{N}, Tuple{JosephsonCircuits.Frequencies{N}, AbstractArray{NTuple{N, Int64}, 2}}} where N"></a>
```

[`JosephsonCircuits.hbmatindices`](api/internals.md#JosephsonCircuits.hbmatindices-Union%7BTuple%7BN%7D%2C%20Tuple%7BJosephsonCircuits.Frequencies%7BN%7D%2C%20AbstractArray%7BNTuple%7BN%2C%20Int64%7D%2C%202%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.hessianvectorproduct!-Tuple{AbstractVector, JosephsonCircuits.HBSystem, AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.hessianvectorproduct!`](api/internals.md#JosephsonCircuits.hessianvectorproduct!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.HBSystem%2C%20AbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.hostrelations-Tuple{JosephsonCircuits.JunctionRelations}"></a>
```

[`JosephsonCircuits.hostrelations`](api/internals.md#JosephsonCircuits.hostrelations-Tuple%7BJosephsonCircuits.JunctionRelations%7D)

```@raw html
<a id="JosephsonCircuits.hostsparse-Tuple{SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.hostsparse`](api/internals.md#JosephsonCircuits.hostsparse-Tuple%7BSparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.hostsystem-Tuple{JosephsonCircuits.HBSystem, Integer, Any, Any}"></a>
```

[`JosephsonCircuits.hostsystem`](api/internals.md#JosephsonCircuits.hostsystem-Tuple%7BJosephsonCircuits.HBSystem%2C%20Integer%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.impedance-Tuple{Any, Integer, Any}"></a>
```

[`JosephsonCircuits.impedance`](api/internals.md#JosephsonCircuits.impedance-Tuple%7BAny%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.impedancecode-Tuple{Any}"></a>
```

[`JosephsonCircuits.impedancecode`](api/internals.md#JosephsonCircuits.impedancecode-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.import_netlist!-Tuple{IO, AbstractVector}"></a>
```

[`JosephsonCircuits.import_netlist!`](api/internals.md#JosephsonCircuits.import_netlist!-Tuple%7BIO%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.import_netlist-Tuple{Any}"></a>
```

[`JosephsonCircuits.import_netlist`](api/internals.md#JosephsonCircuits.import_netlist-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.indefinite_hermitian_form_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.indefinite_hermitian_form_block`](api/internals.md#JosephsonCircuits.indefinite_hermitian_form_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.indefinite_hermitian_form_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.indefinite_hermitian_form_pair`](api/internals.md#JosephsonCircuits.indefinite_hermitian_form_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.initialblockcurrents!-Tuple{AbstractVector, SparseMatrixCSC, Integer, Integer}"></a>
```

[`JosephsonCircuits.initialblockcurrents!`](api/internals.md#JosephsonCircuits.initialblockcurrents!-Tuple%7BAbstractVector%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.initialguess-Tuple{Nothing}"></a>
```

[`JosephsonCircuits.initialguess`](api/internals.md#JosephsonCircuits.initialguess-Tuple%7BNothing%7D)

```@raw html
<a id="JosephsonCircuits.innerpreconditioner"></a>
```

[`JosephsonCircuits.innerpreconditioner`](api/internals.md#JosephsonCircuits.innerpreconditioner)

```@raw html
<a id="JosephsonCircuits.instancedefinition-Tuple{ElaboratedCircuit, Integer}"></a>
```

[`JosephsonCircuits.instancedefinition`](api/internals.md#JosephsonCircuits.instancedefinition-Tuple%7BElaboratedCircuit%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.instanceterminals-Tuple{ElaboratedCircuit, Integer}"></a>
```

[`JosephsonCircuits.instanceterminals`](api/internals.md#JosephsonCircuits.instanceterminals-Tuple%7BElaboratedCircuit%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.interconnectS!-Tuple{Any, Any, Any, Any, Any, Any, Int64, Int64}"></a>
```

[`JosephsonCircuits.interconnectS!`](api/internals.md#JosephsonCircuits.interconnectS!-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Int64%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.interconnectS!-Tuple{Any, Any, Any, Int64, Int64}"></a>
```

[`JosephsonCircuits.interconnectS!`](api/internals.md#JosephsonCircuits.interconnectS!-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Int64%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.interconnectS-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{T, N}, AbstractArray{T, N}, AbstractArray{T, N}, AbstractArray{T, N}, Int64, Int64}} where {T, N}"></a>
```

[`JosephsonCircuits.interconnectS`](api/internals.md#JosephsonCircuits.interconnectS-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%20N%7D%2C%20AbstractArray%7BT%2C%20N%7D%2C%20AbstractArray%7BT%2C%20N%7D%2C%20AbstractArray%7BT%2C%20N%7D%2C%20Int64%2C%20Int64%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.interconnectS-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{T, N}, AbstractArray{T, N}, Int64, Int64}} where {T, N}"></a>
```

[`JosephsonCircuits.interconnectS`](api/internals.md#JosephsonCircuits.interconnectS-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%20N%7D%2C%20AbstractArray%7BT%2C%20N%7D%2C%20Int64%2C%20Int64%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.interconnectS_inner!-Tuple{Any, Any, Any, Any, Any, Any, Int64, Int64, AbstractArray}"></a>
```

[`JosephsonCircuits.interconnectS_inner!`](api/internals.md#JosephsonCircuits.interconnectS_inner!-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Int64%2C%20Int64%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.interconnectS_inner!-Tuple{Any, Any, Any, Int64, Int64, AbstractArray}"></a>
```

[`JosephsonCircuits.interconnectS_inner!`](api/internals.md#JosephsonCircuits.interconnectS_inner!-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Int64%2C%20Int64%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.interconnectSports-Union{Tuple{T}, Tuple{AbstractArray{Tuple{T, Int64}, 1}, AbstractArray{Tuple{T, Int64}, 1}, Int64, Int64}} where T"></a>
```

[`JosephsonCircuits.interconnectSports`](api/internals.md#JosephsonCircuits.interconnectSports-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%2C%20AbstractArray%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%2C%20Int64%2C%20Int64%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.internalpart-Tuple{AbstractVector, JosephsonCircuits.CompositeLayout}"></a>
```

[`JosephsonCircuits.internalpart`](api/internals.md#JosephsonCircuits.internalpart-Tuple%7BAbstractVector%2C%20JosephsonCircuits.CompositeLayout%7D)

```@raw html
<a id="JosephsonCircuits.interpolate_scattering-Tuple{AbstractVector, AbstractArray, AbstractArray}"></a>
```

[`JosephsonCircuits.interpolate_scattering`](api/internals.md#JosephsonCircuits.interpolate_scattering-Tuple%7BAbstractVector%2C%20AbstractArray%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.intraconnectS!-Tuple{Any, Any, Any, Any, Int64, Int64}"></a>
```

[`JosephsonCircuits.intraconnectS!`](api/internals.md#JosephsonCircuits.intraconnectS!-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Any%2C%20Int64%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.intraconnectS!-Tuple{Any, Any, Int64, Int64}"></a>
```

[`JosephsonCircuits.intraconnectS!`](api/internals.md#JosephsonCircuits.intraconnectS!-Tuple%7BAny%2C%20Any%2C%20Int64%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.intraconnectS-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{T, N}, AbstractArray{T, N}, Int64, Int64}} where {T, N}"></a>
```

[`JosephsonCircuits.intraconnectS`](api/internals.md#JosephsonCircuits.intraconnectS-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%20N%7D%2C%20AbstractArray%7BT%2C%20N%7D%2C%20Int64%2C%20Int64%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.intraconnectS-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{T, N}, Int64, Int64}} where {T, N}"></a>
```

[`JosephsonCircuits.intraconnectS`](api/internals.md#JosephsonCircuits.intraconnectS-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%20N%7D%2C%20Int64%2C%20Int64%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.intraconnectS_inner!-Tuple{Any, Any, Int64, Int64, AbstractArray}"></a>
```

[`JosephsonCircuits.intraconnectS_inner!`](api/internals.md#JosephsonCircuits.intraconnectS_inner!-Tuple%7BAny%2C%20Any%2C%20Int64%2C%20Int64%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.intraconnectSports-Union{Tuple{T}, Tuple{AbstractArray{Tuple{T, Int64}, 1}, Int64, Int64}} where T"></a>
```

[`JosephsonCircuits.intraconnectSports`](api/internals.md#JosephsonCircuits.intraconnectSports-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%2C%20Int64%2C%20Int64%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.inv_bogoliubov_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.inv_bogoliubov_block`](api/internals.md#JosephsonCircuits.inv_bogoliubov_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.inv_bogoliubov_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.inv_bogoliubov_pair`](api/internals.md#JosephsonCircuits.inv_bogoliubov_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.inv_symplectic_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.inv_symplectic_block`](api/internals.md#JosephsonCircuits.inv_symplectic_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.inv_symplectic_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.inv_symplectic_pair`](api/internals.md#JosephsonCircuits.inv_symplectic_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.inverseinductanceplan-Tuple{JosephsonCircuits.CircuitTopology, Vector{Int64}, Any}"></a>
```

[`JosephsonCircuits.inverseinductanceplan`](api/internals.md#JosephsonCircuits.inverseinductanceplan-Tuple%7BJosephsonCircuits.CircuitTopology%2C%20Vector%7BInt64%7D%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_bogoliubov_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_bogoliubov_block`](api/internals.md#JosephsonCircuits.is_bogoliubov_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_bogoliubov_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_bogoliubov_pair`](api/internals.md#JosephsonCircuits.is_bogoliubov_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_conjugate_symplectic_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_conjugate_symplectic_block`](api/internals.md#JosephsonCircuits.is_conjugate_symplectic_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_conjugate_symplectic_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_conjugate_symplectic_pair`](api/internals.md#JosephsonCircuits.is_conjugate_symplectic_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_cptp-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.is_cptp`](api/internals.md#JosephsonCircuits.is_cptp-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_cptp_ladder_block-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.is_cptp_ladder_block`](api/internals.md#JosephsonCircuits.is_cptp_ladder_block-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_cptp_ladder_pair-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.is_cptp_ladder_pair`](api/internals.md#JosephsonCircuits.is_cptp_ladder_pair-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_cptp_quadrature_block-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.is_cptp_quadrature_block`](api/internals.md#JosephsonCircuits.is_cptp_quadrature_block-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_cptp_quadrature_pair-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.is_cptp_quadrature_pair`](api/internals.md#JosephsonCircuits.is_cptp_quadrature_pair-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_orthogonal-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_orthogonal`](api/internals.md#JosephsonCircuits.is_orthogonal-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_orthogonal_bogoliubov_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_orthogonal_bogoliubov_block`](api/internals.md#JosephsonCircuits.is_orthogonal_bogoliubov_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_orthogonal_bogoliubov_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_orthogonal_bogoliubov_pair`](api/internals.md#JosephsonCircuits.is_orthogonal_bogoliubov_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_orthogonal_symplectic_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_orthogonal_symplectic_block`](api/internals.md#JosephsonCircuits.is_orthogonal_symplectic_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_orthogonal_symplectic_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_orthogonal_symplectic_pair`](api/internals.md#JosephsonCircuits.is_orthogonal_symplectic_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_positive_definite-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_positive_definite`](api/internals.md#JosephsonCircuits.is_positive_definite-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_positive_definite_symplectic_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_positive_definite_symplectic_block`](api/internals.md#JosephsonCircuits.is_positive_definite_symplectic_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_positive_definite_symplectic_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_positive_definite_symplectic_pair`](api/internals.md#JosephsonCircuits.is_positive_definite_symplectic_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_positive_semi_definite-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_positive_semi_definite`](api/internals.md#JosephsonCircuits.is_positive_semi_definite-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_pseudo_unitary-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.is_pseudo_unitary`](api/internals.md#JosephsonCircuits.is_pseudo_unitary-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_pseudo_unitary_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_pseudo_unitary_block`](api/internals.md#JosephsonCircuits.is_pseudo_unitary_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_pseudo_unitary_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_pseudo_unitary_pair`](api/internals.md#JosephsonCircuits.is_pseudo_unitary_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_symplectic-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.is_symplectic`](api/internals.md#JosephsonCircuits.is_symplectic-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.is_symplectic_block-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_symplectic_block`](api/internals.md#JosephsonCircuits.is_symplectic_block-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_symplectic_pair-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_symplectic_pair`](api/internals.md#JosephsonCircuits.is_symplectic_pair-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.is_unitary-Tuple{Any}"></a>
```

[`JosephsonCircuits.is_unitary`](api/internals.md#JosephsonCircuits.is_unitary-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.isaugmented-Tuple{JosephsonCircuits.HBNonlinearProblem}"></a>
```

[`JosephsonCircuits.isaugmented`](api/internals.md#JosephsonCircuits.isaugmented-Tuple%7BJosephsonCircuits.HBNonlinearProblem%7D)

```@raw html
<a id="JosephsonCircuits.isexactpreconditioner-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.isexactpreconditioner`](api/internals.md#JosephsonCircuits.isexactpreconditioner-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.isgrounded-Tuple{Union{GaussianChannel, LinearizedScattering, ScatteringParameters}}"></a>
```

[`JosephsonCircuits.isgrounded`](api/internals.md#JosephsonCircuits.isgrounded-Tuple%7BUnion%7BGaussianChannel%2C%20LinearizedScattering%2C%20ScatteringParameters%7D%7D)

```@raw html
<a id="JosephsonCircuits.ismnaresistance-Tuple{Any}"></a>
```

[`JosephsonCircuits.ismnaresistance`](api/internals.md#JosephsonCircuits.ismnaresistance-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.isnumericallyzero-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.isnumericallyzero`](api/internals.md#JosephsonCircuits.isnumericallyzero-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.isolatedsubnetworks-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.isolatedsubnetworks`](api/internals.md#JosephsonCircuits.isolatedsubnetworks-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.issinusoidal-Tuple{NonlinearInductor}"></a>
```

[`JosephsonCircuits.issinusoidal`](api/internals.md#JosephsonCircuits.issinusoidal-Tuple%7BNonlinearInductor%7D)

```@raw html
<a id="JosephsonCircuits.iwasawa_block-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.iwasawa_block`](api/internals.md#JosephsonCircuits.iwasawa_block-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.jacobian!-Tuple{JosephsonCircuits.DeviceValuedSparseMatrix, JosephsonCircuits.StructureRealJacobianPlan, JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.jacobian!`](api/internals.md#JosephsonCircuits.jacobian!-Tuple%7BJosephsonCircuits.DeviceValuedSparseMatrix%2C%20JosephsonCircuits.StructureRealJacobianPlan%2C%20JosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.jacobian!-Tuple{SparseMatrixCSC{&lt;:Complex}, JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.jacobian!`](api/internals.md#JosephsonCircuits.jacobian!-Tuple%7BSparseMatrixCSC%7B%3C%3AComplex%7D%2C%20JosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.jacobian!-Tuple{SparseMatrixCSC{&lt;:Real}, JosephsonCircuits.StructureRealJacobianPlan, JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.jacobian!`](api/internals.md#JosephsonCircuits.jacobian!-Tuple%7BSparseMatrixCSC%7B%3C%3AReal%7D%2C%20JosephsonCircuits.StructureRealJacobianPlan%2C%20JosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.jacobianvectorproduct!-Tuple{AbstractVector, JosephsonCircuits.HBSystem, AbstractVector}"></a>
```

[`JosephsonCircuits.jacobianvectorproduct!`](api/internals.md#JosephsonCircuits.jacobianvectorproduct!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.HBSystem%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.jjnodeadjacency-Tuple{SparseVector, Any, Integer}"></a>
```

[`JosephsonCircuits.jjnodeadjacency`](api/internals.md#JosephsonCircuits.jjnodeadjacency-Tuple%7BSparseVector%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.josephsonadjoint!-Tuple{Any, Any, JosephsonCircuits.StructureComplexJosephsonPlan, AbstractVector}"></a>
```

[`JosephsonCircuits.josephsonadjoint!`](api/internals.md#JosephsonCircuits.josephsonadjoint!-Tuple%7BAny%2C%20Any%2C%20JosephsonCircuits.StructureComplexJosephsonPlan%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.junctioncpr"></a>
```

[`JosephsonCircuits.junctioncpr`](api/internals.md#JosephsonCircuits.junctioncpr)

```@raw html
<a id="JosephsonCircuits.junctionpairtable-Union{Tuple{T}, Tuple{Ti}, Tuple{Type{Ti}, Type{T}, SparseVector, Any, Integer}} where {Ti&lt;:Integer, T&lt;:Real}"></a>
```

[`JosephsonCircuits.junctionpairtable`](api/internals.md#JosephsonCircuits.junctionpairtable-Union%7BTuple%7BT%7D%2C%20Tuple%7BTi%7D%2C%20Tuple%7BType%7BTi%7D%2C%20Type%7BT%7D%2C%20SparseVector%2C%20Any%2C%20Integer%7D%7D%20where%20%7BTi%3C%3AInteger%2C%20T%3C%3AReal%7D)

```@raw html
<a id="JosephsonCircuits.junctionrelations-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.junctionrelations`](api/internals.md#JosephsonCircuits.junctionrelations-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.junctionstructure-Union{Tuple{T}, Tuple{Type{T}, Matrix, Matrix, SparseVector, Any, SparseMatrixCSC, Integer, Integer, Integer, Any}} where T&lt;:Real"></a>
```

[`JosephsonCircuits.junctionstructure`](api/internals.md#JosephsonCircuits.junctionstructure-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20Matrix%2C%20Matrix%2C%20SparseVector%2C%20Any%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Integer%2C%20Any%7D%7D%20where%20T%3C%3AReal)

```@raw html
<a id="JosephsonCircuits.klunodeorder-Tuple{AbstractVector{&lt;:AbstractVector{&lt;:Integer}}}"></a>
```

[`JosephsonCircuits.klunodeorder`](api/internals.md#JosephsonCircuits.klunodeorder-Tuple%7BAbstractVector%7B%3C%3AAbstractVector%7B%3C%3AInteger%7D%7D%7D)

```@raw html
<a id="JosephsonCircuits.kluordered-Union{Tuple{SparseMatrixCSC{Tv, Ti}}, Tuple{Ti}, Tuple{Tv}, Tuple{SparseMatrixCSC{Tv, Ti}, Any}} where {Tv, Ti}"></a>
```

[`JosephsonCircuits.kluordered`](api/internals.md#JosephsonCircuits.kluordered-Union%7BTuple%7BSparseMatrixCSC%7BTv%2C%20Ti%7D%7D%2C%20Tuple%7BTi%7D%2C%20Tuple%7BTv%7D%2C%20Tuple%7BSparseMatrixCSC%7BTv%2C%20Ti%7D%2C%20Any%7D%7D%20where%20%7BTv%2C%20Ti%7D)

```@raw html
<a id="JosephsonCircuits.klupivotgrowth-Union{Tuple{KLU.KLUFactorization{Tv}}, Tuple{Tv}} where Tv"></a>
```

[`JosephsonCircuits.klupivotgrowth`](api/internals.md#JosephsonCircuits.klupivotgrowth-Union%7BTuple%7BKLU.KLUFactorization%7BTv%7D%7D%2C%20Tuple%7BTv%7D%7D%20where%20Tv)

```@raw html
<a id="JosephsonCircuits.klurefactor!-Tuple{KLU.KLUFactorization, SparseMatrixCSC, Real}"></a>
```

[`JosephsonCircuits.klurefactor!`](api/internals.md#JosephsonCircuits.klurefactor!-Tuple%7BKLU.KLUFactorization%2C%20SparseMatrixCSC%2C%20Real%7D)

```@raw html
<a id="JosephsonCircuits.ladder_to_quadrature_block-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.ladder_to_quadrature_block`](api/internals.md#JosephsonCircuits.ladder_to_quadrature_block-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.ladder_to_quadrature_block-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.ladder_to_quadrature_block`](api/internals.md#JosephsonCircuits.ladder_to_quadrature_block-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.ladder_to_quadrature_pair-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.ladder_to_quadrature_pair`](api/internals.md#JosephsonCircuits.ladder_to_quadrature_pair-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.ladder_to_quadrature_pair-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.ladder_to_quadrature_pair`](api/internals.md#JosephsonCircuits.ladder_to_quadrature_pair-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.ladder_to_scattering_block-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Any}} where T&lt;:Union{AbstractFloat, Complex{&lt;:AbstractFloat}}"></a>
```

[`JosephsonCircuits.ladder_to_scattering_block`](api/internals.md#JosephsonCircuits.ladder_to_scattering_block-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractMatrix%7BT%7D%2C%20Any%7D%7D%20where%20T%3C%3AUnion%7BAbstractFloat%2C%20Complex%7B%3C%3AAbstractFloat%7D%7D)

```@raw html
<a id="JosephsonCircuits.ladder_to_scattering_pair-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Any}} where T&lt;:Union{AbstractFloat, Complex{&lt;:AbstractFloat}}"></a>
```

[`JosephsonCircuits.ladder_to_scattering_pair`](api/internals.md#JosephsonCircuits.ladder_to_scattering_pair-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractMatrix%7BT%7D%2C%20Any%7D%7D%20where%20T%3C%3AUnion%7BAbstractFloat%2C%20Complex%7B%3C%3AAbstractFloat%7D%7D)

```@raw html
<a id="JosephsonCircuits.ldiv_2x2-Tuple{Union{LU, StaticArrays.LU}, AbstractVector}"></a>
```

[`JosephsonCircuits.ldiv_2x2`](api/internals.md#JosephsonCircuits.ldiv_2x2-Tuple%7BUnion%7BLU%2C%20StaticArrays.LU%7D%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.linearentry-NTuple{15, Any}"></a>
```

[`JosephsonCircuits.linearentry`](api/internals.md#JosephsonCircuits.linearentry-NTuple%7B15%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.lineargather-Union{Tuple{Ti}, Tuple{Type{Ti}, Any, Any, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, Diagonal, Diagonal, Any, Any, Bool, Any}} where Ti"></a>
```

[`JosephsonCircuits.lineargather`](api/internals.md#JosephsonCircuits.lineargather-Union%7BTuple%7BTi%7D%2C%20Tuple%7BType%7BTi%7D%2C%20Any%2C%20Any%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Diagonal%2C%20Diagonal%2C%20Any%2C%20Any%2C%20Bool%2C%20Any%7D%7D%20where%20Ti)

```@raw html
<a id="JosephsonCircuits.linearizedfactorization-Tuple{SparseMatrixCSC, Integer, Integer, Any}"></a>
```

[`JosephsonCircuits.linearizedfactorization`](api/internals.md#JosephsonCircuits.linearizedfactorization-Tuple%7BSparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.linearizedoutputs-Tuple{}"></a>
```

[`JosephsonCircuits.linearizedoutputs`](api/internals.md#JosephsonCircuits.linearizedoutputs-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.linearizedsensitivity-Tuple{}"></a>
```

[`JosephsonCircuits.linearizedsensitivity`](api/internals.md#JosephsonCircuits.linearizedsensitivity-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.linearizedsetup-Tuple{Vector{Float64}, JosephsonCircuits.CompiledCircuit, Vector{Any}, JosephsonCircuits.Frequencies, Any, Any, Any, Any, Vector{String}, Vector{Tuple{String, Int64, ComplexF64}}, Vector{Tuple{String, Int64, Any}}}"></a>
```

[`JosephsonCircuits.linearizedsetup`](api/internals.md#JosephsonCircuits.linearizedsetup-Tuple%7BVector%7BFloat64%7D%2C%20JosephsonCircuits.CompiledCircuit%2C%20Vector%7BAny%7D%2C%20JosephsonCircuits.Frequencies%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Vector%7BString%7D%2C%20Vector%7BTuple%7BString%2C%20Int64%2C%20ComplexF64%7D%7D%2C%20Vector%7BTuple%7BString%2C%20Int64%2C%20Any%7D%7D%7D)

```@raw html
<a id="JosephsonCircuits.linearizedsweep!-Tuple{}"></a>
```

[`JosephsonCircuits.linearizedsweep!`](api/internals.md#JosephsonCircuits.linearizedsweep!-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.linearterm!-Tuple{SparseMatrixCSC, JosephsonCircuits.ValueMaps, Any, Any, Any, Diagonal, Diagonal}"></a>
```

[`JosephsonCircuits.linearterm!`](api/internals.md#JosephsonCircuits.linearterm!-Tuple%7BSparseMatrixCSC%2C%20JosephsonCircuits.ValueMaps%2C%20Any%2C%20Any%2C%20Any%2C%20Diagonal%2C%20Diagonal%7D)

```@raw html
<a id="JosephsonCircuits.linearterm-Union{Tuple{T}, Tuple{Any, Any, Any, Diagonal, Diagonal}, Tuple{Any, Any, Any, Diagonal, Diagonal, Type{T}}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.linearterm`](api/internals.md#JosephsonCircuits.linearterm-Union%7BTuple%7BT%7D%2C%20Tuple%7BAny%2C%20Any%2C%20Any%2C%20Diagonal%2C%20Diagonal%7D%2C%20Tuple%7BAny%2C%20Any%2C%20Any%2C%20Diagonal%2C%20Diagonal%2C%20Type%7BT%7D%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.linesearchevaluate!-Tuple{Any, AbstractVector, AbstractVector, AbstractVector, Real, AbstractVector, Real, Union{Nothing, AbstractVector}}"></a>
```

[`JosephsonCircuits.linesearchevaluate!`](api/internals.md#JosephsonCircuits.linesearchevaluate!-Tuple%7BAny%2C%20AbstractVector%2C%20AbstractVector%2C%20AbstractVector%2C%20Real%2C%20AbstractVector%2C%20Real%2C%20Union%7BNothing%2C%20AbstractVector%7D%7D)

```@raw html
<a id="JosephsonCircuits.linesearchtrialpoint!-Tuple{AbstractVector, AbstractVector, Real, AbstractVector, Real, Union{Nothing, AbstractVector}}"></a>
```

[`JosephsonCircuits.linesearchtrialpoint!`](api/internals.md#JosephsonCircuits.linesearchtrialpoint!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20Real%2C%20AbstractVector%2C%20Real%2C%20Union%7BNothing%2C%20AbstractVector%7D%7D)

```@raw html
<a id="JosephsonCircuits.lu_2x2-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.lu_2x2`](api/internals.md#JosephsonCircuits.lu_2x2-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.make_connection!-Union{Tuple{N}, Tuple{T}, Tuple{Graphs.SimpleGraphs.SimpleDiGraph{Int64}, AbstractVector{&lt;:AbstractArray{Tuple{T, T, Int64, Int64}, 1}}, AbstractVector{&lt;:AbstractVector{Int64}}, AbstractVector{&lt;:AbstractArray{Tuple{T, Int64}, 1}}, AbstractVector{N}, AbstractVector{N}, Int64, Int64, Int64, AbstractVector{Bool}, Dict{Int64, N}, Dict{Int64, N}, Bool}} where {T, N}"></a>
```

[`JosephsonCircuits.make_connection!`](api/internals.md#JosephsonCircuits.make_connection!-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BGraphs.SimpleGraphs.SimpleDiGraph%7BInt64%7D%2C%20AbstractVector%7B%3C%3AAbstractArray%7BTuple%7BT%2C%20T%2C%20Int64%2C%20Int64%7D%2C%201%7D%7D%2C%20AbstractVector%7B%3C%3AAbstractVector%7BInt64%7D%7D%2C%20AbstractVector%7B%3C%3AAbstractArray%7BTuple%7BT%2C%20Int64%7D%2C%201%7D%7D%2C%20AbstractVector%7BN%7D%2C%20AbstractVector%7BN%7D%2C%20Int64%2C%20Int64%2C%20Int64%2C%20AbstractVector%7BBool%7D%2C%20Dict%7BInt64%2C%20N%7D%2C%20Dict%7BInt64%2C%20N%7D%2C%20Bool%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.matrixprovider-Union{Tuple{T}, Tuple{JosephsonCircuits.AbstractMatrixProvider, Type{T}}} where T"></a>
```

[`JosephsonCircuits.matrixprovider`](api/internals.md#JosephsonCircuits.matrixprovider-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.AbstractMatrixProvider%2C%20Type%7BT%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.maxwell_combine-Tuple{Int64, AbstractDict{&lt;:Tuple{Vararg{Int64}}, &lt;:AbstractMatrix}}"></a>
```

[`JosephsonCircuits.maxwell_combine`](api/internals.md#JosephsonCircuits.maxwell_combine-Tuple%7BInt64%2C%20AbstractDict%7B%3C%3ATuple%7BVararg%7BInt64%7D%7D%2C%20%3C%3AAbstractMatrix%7D%7D)

```@raw html
<a id="JosephsonCircuits.maxwell_to_even_odd-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.maxwell_to_even_odd`](api/internals.md#JosephsonCircuits.maxwell_to_even_odd-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.maxwell_to_mutual-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.maxwell_to_mutual`](api/internals.md#JosephsonCircuits.maxwell_to_mutual-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.mergestamps-Tuple{AbstractVector{JosephsonCircuits.SensitivityStamp}, Any}"></a>
```

[`JosephsonCircuits.mergestamps`](api/internals.md#JosephsonCircuits.mergestamps-Tuple%7BAbstractVector%7BJosephsonCircuits.SensitivityStamp%7D%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.merit-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.merit`](api/internals.md#JosephsonCircuits.merit-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.meritslope!-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.meritslope!`](api/internals.md#JosephsonCircuits.meritslope!-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.mnacoupledbranches-Tuple{SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.mnacoupledbranches`](api/internals.md#JosephsonCircuits.mnacoupledbranches-Tuple%7BSparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.mnagaugenormalize!-Tuple{AbstractVector, Vector{Vector{Int64}}, Vector, Int64}"></a>
```

[`JosephsonCircuits.mnagaugenormalize!`](api/internals.md#JosephsonCircuits.mnagaugenormalize!-Tuple%7BAbstractVector%2C%20Vector%7BVector%7BInt64%7D%7D%2C%20Vector%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.mnainitialauxind!-Tuple{AbstractVector, Vector{Int64}, SparseVector, SparseMatrixCSC, SparseMatrixCSC, Int64, Int64, Any}"></a>
```

[`JosephsonCircuits.mnainitialauxind!`](api/internals.md#JosephsonCircuits.mnainitialauxind!-Tuple%7BAbstractVector%2C%20Vector%7BInt64%7D%2C%20SparseVector%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Int64%2C%20Int64%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.mnapad-Tuple{SparseMatrixCSC, Int64}"></a>
```

[`JosephsonCircuits.mnapad`](api/internals.md#JosephsonCircuits.mnapad-Tuple%7BSparseMatrixCSC%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.mnaresistance-Tuple{Real}"></a>
```

[`JosephsonCircuits.mnaresistance`](api/internals.md#JosephsonCircuits.mnaresistance-Tuple%7BReal%7D)

```@raw html
<a id="JosephsonCircuits.mnaungaugedkcl-Tuple{AbstractVector, AbstractVector, Vector{Int64}, Int64}"></a>
```

[`JosephsonCircuits.mnaungaugedkcl`](api/internals.md#JosephsonCircuits.mnaungaugedkcl-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20Vector%7BInt64%7D%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.mnavalidatekcl-Tuple{AbstractVector, AbstractVector, Vector{Int64}, Int64, AbstractVector, Any}"></a>
```

[`JosephsonCircuits.mnavalidatekcl`](api/internals.md#JosephsonCircuits.mnavalidatekcl-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20Vector%7BInt64%7D%2C%20Int64%2C%20AbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.modebandmask-Tuple{AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.modebandmask`](api/internals.md#JosephsonCircuits.modebandmask-Tuple%7BAbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.modeclusters-Tuple{AbstractMatrix{Bool}}"></a>
```

[`JosephsonCircuits.modeclusters`](api/internals.md#JosephsonCircuits.modeclusters-Tuple%7BAbstractMatrix%7BBool%7D%7D)

```@raw html
<a id="JosephsonCircuits.modecouplingmask-Tuple{Integer, Any}"></a>
```

[`JosephsonCircuits.modecouplingmask`](api/internals.md#JosephsonCircuits.modecouplingmask-Tuple%7BInteger%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.modes_ports_to_ports_modes_block-Tuple{AbstractMatrix, Int64}"></a>
```

[`JosephsonCircuits.modes_ports_to_ports_modes_block`](api/internals.md#JosephsonCircuits.modes_ports_to_ports_modes_block-Tuple%7BAbstractMatrix%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.modes_ports_to_ports_modes_pair-Tuple{AbstractMatrix, Int64}"></a>
```

[`JosephsonCircuits.modes_ports_to_ports_modes_pair`](api/internals.md#JosephsonCircuits.modes_ports_to_ports_modes_pair-Tuple%7BAbstractMatrix%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.modes_ports_to_ports_modes_perm-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.modes_ports_to_ports_modes_perm`](api/internals.md#JosephsonCircuits.modes_ports_to_ports_modes_perm-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.modes_ports_to_ports_modes_scattering-Tuple{AbstractMatrix, Int64}"></a>
```

[`JosephsonCircuits.modes_ports_to_ports_modes_scattering`](api/internals.md#JosephsonCircuits.modes_ports_to_ports_modes_scattering-Tuple%7BAbstractMatrix%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.modeslotindex-Tuple{JosephsonCircuits.ModeLayout}"></a>
```

[`JosephsonCircuits.modeslotindex`](api/internals.md#JosephsonCircuits.modeslotindex-Tuple%7BJosephsonCircuits.ModeLayout%7D)

```@raw html
<a id="JosephsonCircuits.modevalue-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.modevalue`](api/internals.md#JosephsonCircuits.modevalue-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.move_bedge!-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.move_bedge!`](api/internals.md#JosephsonCircuits.move_bedge!-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.move_bedges!-NTuple{5, Any}"></a>
```

[`JosephsonCircuits.move_bedges!`](api/internals.md#JosephsonCircuits.move_bedges!-NTuple%7B5%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.move_edges!-NTuple{5, Any}"></a>
```

[`JosephsonCircuits.move_edges!`](api/internals.md#JosephsonCircuits.move_edges!-NTuple%7B5%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.move_fedge!-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.move_fedge!`](api/internals.md#JosephsonCircuits.move_fedge!-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.move_fedges!-NTuple{5, Any}"></a>
```

[`JosephsonCircuits.move_fedges!`](api/internals.md#JosephsonCircuits.move_fedges!-NTuple%7B5%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.mutual_to_even_odd-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.mutual_to_even_odd`](api/internals.md#JosephsonCircuits.mutual_to_even_odd-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.mutual_to_maxwell-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.mutual_to_maxwell`](api/internals.md#JosephsonCircuits.mutual_to_maxwell-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.mutualstampplan-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.mutualstampplan`](api/internals.md#JosephsonCircuits.mutualstampplan-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.ncomponents-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.ncomponents`](api/internals.md#JosephsonCircuits.ncomponents-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.needsadjointsolve"></a>
```

[`JosephsonCircuits.needsadjointsolve`](api/internals.md#JosephsonCircuits.needsadjointsolve)

```@raw html
<a id="JosephsonCircuits.negsecondat-Tuple{JosephsonCircuits.JunctionRelations, Any}"></a>
```

[`JosephsonCircuits.negsecondat`](api/internals.md#JosephsonCircuits.negsecondat-Tuple%7BJosephsonCircuits.JunctionRelations%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.ninstances-Tuple{ElaboratedCircuit}"></a>
```

[`JosephsonCircuits.ninstances`](api/internals.md#JosephsonCircuits.ninstances-Tuple%7BElaboratedCircuit%7D)

```@raw html
<a id="JosephsonCircuits.nlsolve!-Union{Tuple{T}, Tuple{Function, AbstractVector{T}, AbstractArray{T}, AbstractVector{T}}} where T"></a>
```

[`JosephsonCircuits.nlsolve!`](api/internals.md#JosephsonCircuits.nlsolve!-Union%7BTuple%7BT%7D%2C%20Tuple%7BFunction%2C%20AbstractVector%7BT%7D%2C%20AbstractArray%7BT%7D%2C%20AbstractVector%7BT%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.nlsolvekrylov!-Union{Tuple{T}, Tuple{Function, Any, AbstractVector{T}, AbstractVector{T}, JosephsonCircuits.AbstractPreconditioner}, Tuple{Function, Any, AbstractVector{T}, AbstractVector{T}, JosephsonCircuits.AbstractPreconditioner, NewtonKrylov}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.nlsolvekrylov!`](api/internals.md#JosephsonCircuits.nlsolvekrylov!-Union%7BTuple%7BT%7D%2C%20Tuple%7BFunction%2C%20Any%2C%20AbstractVector%7BT%7D%2C%20AbstractVector%7BT%7D%2C%20JosephsonCircuits.AbstractPreconditioner%7D%2C%20Tuple%7BFunction%2C%20Any%2C%20AbstractVector%7BT%7D%2C%20AbstractVector%7BT%7D%2C%20JosephsonCircuits.AbstractPreconditioner%2C%20NewtonKrylov%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.nnets-Tuple{ElaboratedCircuit}"></a>
```

[`JosephsonCircuits.nnets`](api/internals.md#JosephsonCircuits.nnets-Tuple%7BElaboratedCircuit%7D)

```@raw html
<a id="JosephsonCircuits.nodalstampplan-Tuple{JosephsonCircuits.CompiledCircuit, Vector{Int64}, Int64}"></a>
```

[`JosephsonCircuits.nodalstampplan`](api/internals.md#JosephsonCircuits.nodalstampplan-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Vector%7BInt64%7D%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.nodecomponents-Tuple{Int64, Any}"></a>
```

[`JosephsonCircuits.nodecomponents`](api/internals.md#JosephsonCircuits.nodecomponents-Tuple%7BInt64%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.noderenumbering-Tuple{Vector{Int64}}"></a>
```

[`JosephsonCircuits.noderenumbering`](api/internals.md#JosephsonCircuits.noderenumbering-Tuple%7BVector%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.nodevariabletokeyed-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.nodevariabletokeyed`](api/internals.md#JosephsonCircuits.nodevariabletokeyed-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.nodevariabletokeyed-Tuple{AbstractVector, Any}"></a>
```

[`JosephsonCircuits.nodevariabletokeyed`](api/internals.md#JosephsonCircuits.nodevariabletokeyed-Tuple%7BAbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.nodevariabletokeyed-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.nodevariabletokeyed`](api/internals.md#JosephsonCircuits.nodevariabletokeyed-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.noisechannelnames-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.noisechannelnames`](api/internals.md#JosephsonCircuits.noisechannelnames-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.noisechannelsigns-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.noisechannelsigns`](api/internals.md#JosephsonCircuits.noisechannelsigns-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.noisechanneltemperatures-NTuple{5, Any}"></a>
```

[`JosephsonCircuits.noisechanneltemperatures`](api/internals.md#JosephsonCircuits.noisechanneltemperatures-NTuple%7B5%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.noisecovariance!-Tuple{Any, Integer, Integer, Any, Integer}"></a>
```

[`JosephsonCircuits.noisecovariance!`](api/internals.md#JosephsonCircuits.noisecovariance!-Tuple%7BAny%2C%20Integer%2C%20Integer%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.noiseindices"></a>
```

[`JosephsonCircuits.noiseindices`](api/internals.md#JosephsonCircuits.noiseindices)

```@raw html
<a id="JosephsonCircuits.noiseoutputwavekernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.noiseoutputwavekernel!`](api/internals.md#JosephsonCircuits.noiseoutputwavekernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.noisereduction!"></a>
```

[`JosephsonCircuits.noisereduction!`](api/internals.md#JosephsonCircuits.noisereduction!)

```@raw html
<a id="JosephsonCircuits.noisewavescale-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.noisewavescale`](api/internals.md#JosephsonCircuits.noisewavescale-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.nonlinearmatrices-Tuple{JosephsonCircuits.CircuitMatrices, Any, Any, Any}"></a>
```

[`JosephsonCircuits.nonlinearmatrices`](api/internals.md#JosephsonCircuits.nonlinearmatrices-Tuple%7BJosephsonCircuits.CircuitMatrices%2C%20Any%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.nonlinearoutputs-Tuple{}"></a>
```

[`JosephsonCircuits.nonlinearoutputs`](api/internals.md#JosephsonCircuits.nonlinearoutputs-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.nonlinearsetup-Union{Tuple{N}, Tuple{NTuple{N, Float64}, Array{@NamedTuple{mode::NTuple{N, Int64}, port::Int64, current::ComplexF64}, 1}, JosephsonCircuits.Frequencies{N}, JosephsonCircuits.FourierIndices{N}, JosephsonCircuits.CompiledCircuit, Any, Vector{ComplexF64}, Any, Type{&lt;:AbstractFloat}, Union{Nothing, HBReuse}}} where N"></a>
```

[`JosephsonCircuits.nonlinearsetup`](api/internals.md#JosephsonCircuits.nonlinearsetup-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Float64%7D%2C%20Array%7B%40NamedTuple%7Bmode%3A%3ANTuple%7BN%2C%20Int64%7D%2C%20port%3A%3AInt64%2C%20current%3A%3AComplexF64%7D%2C%201%7D%2C%20JosephsonCircuits.Frequencies%7BN%7D%2C%20JosephsonCircuits.FourierIndices%7BN%7D%2C%20JosephsonCircuits.CompiledCircuit%2C%20Any%2C%20Vector%7BComplexF64%7D%2C%20Any%2C%20Type%7B%3C%3AAbstractFloat%7D%2C%20Union%7BNothing%2C%20HBReuse%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.norm2-Tuple{AbstractVector{Float64}}"></a>
```

[`JosephsonCircuits.norm2`](api/internals.md#JosephsonCircuits.norm2-Tuple%7BAbstractVector%7BFloat64%7D%7D)

```@raw html
<a id="JosephsonCircuits.normalizedefinitions-Tuple{Any}"></a>
```

[`JosephsonCircuits.normalizedefinitions`](api/internals.md#JosephsonCircuits.normalizedefinitions-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.nterminals-Tuple{JosephsonCircuits.GroundType}"></a>
```

[`JosephsonCircuits.nterminals`](api/internals.md#JosephsonCircuits.nterminals-Tuple%7BJosephsonCircuits.GroundType%7D)

```@raw html
<a id="JosephsonCircuits.numericvalues-Tuple{JosephsonCircuits.CompiledCircuit, Any}"></a>
```

[`JosephsonCircuits.numericvalues`](api/internals.md#JosephsonCircuits.numericvalues-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.nvoltages-Tuple{JosephsonCircuits.TransportRows}"></a>
```

[`JosephsonCircuits.nvoltages`](api/internals.md#JosephsonCircuits.nvoltages-Tuple%7BJosephsonCircuits.TransportRows%7D)

```@raw html
<a id="JosephsonCircuits.nwindow-Tuple{JosephsonCircuits.CompositeLayout}"></a>
```

[`JosephsonCircuits.nwindow`](api/internals.md#JosephsonCircuits.nwindow-Tuple%7BJosephsonCircuits.CompositeLayout%7D)

```@raw html
<a id="JosephsonCircuits.nzposition-Tuple{SparseMatrixCSC, Integer, Integer}"></a>
```

[`JosephsonCircuits.nzposition`](api/internals.md#JosephsonCircuits.nzposition-Tuple%7BSparseMatrixCSC%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.optimum_eigenvalue_angle-Tuple{Any}"></a>
```

[`JosephsonCircuits.optimum_eigenvalue_angle`](api/internals.md#JosephsonCircuits.optimum_eigenvalue_angle-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.orderedports-Tuple{JosephsonCircuits.CompiledCircuit}"></a>
```

[`JosephsonCircuits.orderedports`](api/internals.md#JosephsonCircuits.orderedports-Tuple%7BJosephsonCircuits.CompiledCircuit%7D)

```@raw html
<a id="JosephsonCircuits.outputnoise!"></a>
```

[`JosephsonCircuits.outputnoise!`](api/internals.md#JosephsonCircuits.outputnoise!)

```@raw html
<a id="JosephsonCircuits.pair_to_block-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.pair_to_block`](api/internals.md#JosephsonCircuits.pair_to_block-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.pair_to_block-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.pair_to_block`](api/internals.md#JosephsonCircuits.pair_to_block-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.pair_to_block2-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.pair_to_block2`](api/internals.md#JosephsonCircuits.pair_to_block2-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.pair_to_block2-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.pair_to_block2`](api/internals.md#JosephsonCircuits.pair_to_block2-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.pair_to_block_perm-Tuple{Integer}"></a>
```

[`JosephsonCircuits.pair_to_block_perm`](api/internals.md#JosephsonCircuits.pair_to_block_perm-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.parametergrouping-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.parametergrouping`](api/internals.md#JosephsonCircuits.parametergrouping-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.parse_connections_sparse-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{JosephsonCircuits.PassiveNetwork{T, N}, 1}, AbstractArray{Tuple{T, T, Int64, Int64}, 1}}} where {T, N}"></a>
```

[`JosephsonCircuits.parse_connections_sparse`](api/internals.md#JosephsonCircuits.parse_connections_sparse-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BJosephsonCircuits.PassiveNetwork%7BT%2C%20N%7D%2C%201%7D%2C%20AbstractArray%7BTuple%7BT%2C%20T%2C%20Int64%2C%20Int64%7D%2C%201%7D%7D%7D%20where%20%7BT%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.parsecircuitlevel"></a>
```

[`JosephsonCircuits.parsecircuitlevel`](api/internals.md#JosephsonCircuits.parsecircuitlevel)

```@raw html
<a id="JosephsonCircuits.parsecomponenttype-Tuple{String, Vector{String}}"></a>
```

[`JosephsonCircuits.parsecomponenttype`](api/internals.md#JosephsonCircuits.parsecomponenttype-Tuple%7BString%2C%20Vector%7BString%7D%7D)

```@raw html
<a id="JosephsonCircuits.parsecomponentvalue-Tuple{AbstractString}"></a>
```

[`JosephsonCircuits.parsecomponentvalue`](api/internals.md#JosephsonCircuits.parsecomponentvalue-Tuple%7BAbstractString%7D)

```@raw html
<a id="JosephsonCircuits.parsespicevariable-Tuple{String}"></a>
```

[`JosephsonCircuits.parsespicevariable`](api/internals.md#JosephsonCircuits.parsespicevariable-Tuple%7BString%7D)

```@raw html
<a id="JosephsonCircuits.passiveconstant-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.passiveconstant`](api/internals.md#JosephsonCircuits.passiveconstant-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.passivitymargin-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.passivitymargin`](api/internals.md#JosephsonCircuits.passivitymargin-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.perfrequency!-Union{Tuple{N}, Tuple{F}, Tuple{F, AbstractArray, Vararg{Any, N}}} where {F, N}"></a>
```

[`JosephsonCircuits.perfrequency!`](api/internals.md#JosephsonCircuits.perfrequency!-Union%7BTuple%7BN%7D%2C%20Tuple%7BF%7D%2C%20Tuple%7BF%2C%20AbstractArray%2C%20Vararg%7BAny%2C%20N%7D%7D%7D%20where%20%7BF%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.perfrequency-Union{Tuple{N}, Tuple{G}, Tuple{F}, Tuple{F, G, Vararg{Any, N}}} where {F, G, N}"></a>
```

[`JosephsonCircuits.perfrequency`](api/internals.md#JosephsonCircuits.perfrequency-Union%7BTuple%7BN%7D%2C%20Tuple%7BG%7D%2C%20Tuple%7BF%7D%2C%20Tuple%7BF%2C%20G%2C%20Vararg%7BAny%2C%20N%7D%7D%7D%20where%20%7BF%2C%20G%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.perronbound-Tuple{AbstractMatrix, Any, Any}"></a>
```

[`JosephsonCircuits.perronbound`](api/internals.md#JosephsonCircuits.perronbound-Tuple%7BAbstractMatrix%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.phivectortomatrix!-Tuple{AbstractVector, AbstractArray, Vector{Int64}, Vector{Int64}, Vector{Int64}, Int64}"></a>
```

[`JosephsonCircuits.phivectortomatrix!`](api/internals.md#JosephsonCircuits.phivectortomatrix!-Tuple%7BAbstractVector%2C%20AbstractArray%2C%20Vector%7BInt64%7D%2C%20Vector%7BInt64%7D%2C%20Vector%7BInt64%7D%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.pivot_rows-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.pivot_rows`](api/internals.md#JosephsonCircuits.pivot_rows-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.pivot_rows-Union{Tuple{T}, Tuple{Union{Complex{T}, T}, Union{Complex{T}, T}}} where T&lt;:AbstractFloat"></a>
```

[`JosephsonCircuits.pivot_rows`](api/internals.md#JosephsonCircuits.pivot_rows-Union%7BTuple%7BT%7D%2C%20Tuple%7BUnion%7BComplex%7BT%7D%2C%20T%7D%2C%20Union%7BComplex%7BT%7D%2C%20T%7D%7D%7D%20where%20T%3C%3AAbstractFloat)

```@raw html
<a id="JosephsonCircuits.plan_applyffttranspose-Union{Tuple{Array{T}}, Tuple{T}} where T"></a>
```

[`JosephsonCircuits.plan_applyffttranspose`](api/internals.md#JosephsonCircuits.plan_applyffttranspose-Union%7BTuple%7BArray%7BT%7D%7D%2C%20Tuple%7BT%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.plan_applynl-Union{Tuple{AbstractArray{Complex{T}}}, Tuple{T}, Tuple{AbstractArray{Complex{T}}, KernelAbstractions.Backend}} where T"></a>
```

[`JosephsonCircuits.plan_applynl`](api/internals.md#JosephsonCircuits.plan_applynl-Union%7BTuple%7BAbstractArray%7BComplex%7BT%7D%7D%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BComplex%7BT%7D%7D%2C%20KernelAbstractions.Backend%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.plancomplexjacobian-Tuple{Matrix, SparseVector, Any, SparseMatrixCSC, Integer, Integer, Integer, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.plancomplexjacobian`](api/internals.md#JosephsonCircuits.plancomplexjacobian-Tuple%7BMatrix%2C%20SparseVector%2C%20Any%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Integer%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.plandeviceblocknoise-Tuple{Any, Nothing, Any, Any}"></a>
```

[`JosephsonCircuits.plandeviceblocknoise`](api/internals.md#JosephsonCircuits.plandeviceblocknoise-Tuple%7BAny%2C%20Nothing%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.plandevicenoise-Tuple{Any, Any, Any, Any, Integer, Any}"></a>
```

[`JosephsonCircuits.plandevicenoise`](api/internals.md#JosephsonCircuits.plandevicenoise-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Any%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.plandeviceproviders-Tuple{Any, Integer, Any, Any, Real}"></a>
```

[`JosephsonCircuits.plandeviceproviders`](api/internals.md#JosephsonCircuits.plandeviceproviders-Tuple%7BAny%2C%20Integer%2C%20Any%2C%20Any%2C%20Real%7D)

```@raw html
<a id="JosephsonCircuits.plandevicescattering-Tuple{Any, AbstractVector{Int64}, Integer, Integer, Any, Integer}"></a>
```

[`JosephsonCircuits.plandevicescattering`](api/internals.md#JosephsonCircuits.plandevicescattering-Tuple%7BAny%2C%20AbstractVector%7BInt64%7D%2C%20Integer%2C%20Integer%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.planequilibration-Tuple{Any, AbstractVector{&lt;:Integer}, Integer, Integer, Any}"></a>
```

[`JosephsonCircuits.planequilibration`](api/internals.md#JosephsonCircuits.planequilibration-Tuple%7BAny%2C%20AbstractVector%7B%3C%3AInteger%7D%2C%20Integer%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.planfrequencysweep-Tuple{JosephsonCircuits.HBLinearizedSystem, Any}"></a>
```

[`JosephsonCircuits.planfrequencysweep`](api/internals.md#JosephsonCircuits.planfrequencysweep-Tuple%7BJosephsonCircuits.HBLinearizedSystem%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.plannonlinearterm"></a>
```

[`JosephsonCircuits.plannonlinearterm`](api/internals.md#JosephsonCircuits.plannonlinearterm)

```@raw html
<a id="JosephsonCircuits.plannonlineartermtranspose-Tuple{JosephsonCircuits.NonlinearTermPlan, Any, AbstractArray, AbstractArray}"></a>
```

[`JosephsonCircuits.plannonlineartermtranspose`](api/internals.md#JosephsonCircuits.plannonlineartermtranspose-Tuple%7BJosephsonCircuits.NonlinearTermPlan%2C%20Any%2C%20AbstractArray%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.planscatteringnoise-Tuple{Nothing}"></a>
```

[`JosephsonCircuits.planscatteringnoise`](api/internals.md#JosephsonCircuits.planscatteringnoise-Tuple%7BNothing%7D)

```@raw html
<a id="JosephsonCircuits.planstructurecomplexjacobian-Union{Tuple{T}, Tuple{SparseMatrixCSC, JosephsonCircuits.JunctionStructure{T}, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, Diagonal, Diagonal, Any}} where T&lt;:Real"></a>
```

[`JosephsonCircuits.planstructurecomplexjacobian`](api/internals.md#JosephsonCircuits.planstructurecomplexjacobian-Union%7BTuple%7BT%7D%2C%20Tuple%7BSparseMatrixCSC%2C%20JosephsonCircuits.JunctionStructure%7BT%7D%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Diagonal%2C%20Diagonal%2C%20Any%7D%7D%20where%20T%3C%3AReal)

```@raw html
<a id="JosephsonCircuits.planstructurecomplexjosephson-Tuple{SparseMatrixCSC, JosephsonCircuits.JunctionStructure, Any}"></a>
```

[`JosephsonCircuits.planstructurecomplexjosephson`](api/internals.md#JosephsonCircuits.planstructurecomplexjosephson-Tuple%7BSparseMatrixCSC%2C%20JosephsonCircuits.JunctionStructure%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.planstructurerealjacobian-Union{Tuple{T}, Tuple{Any, Type{T}, JosephsonCircuits.JunctionStructure{T}, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, Diagonal, Diagonal, JosephsonCircuits.ModeLayout, Any}} where T&lt;:Real"></a>
```

[`JosephsonCircuits.planstructurerealjacobian`](api/internals.md#JosephsonCircuits.planstructurerealjacobian-Union%7BTuple%7BT%7D%2C%20Tuple%7BAny%2C%20Type%7BT%7D%2C%20JosephsonCircuits.JunctionStructure%7BT%7D%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Diagonal%2C%20Diagonal%2C%20JosephsonCircuits.ModeLayout%2C%20Any%7D%7D%20where%20T%3C%3AReal)

```@raw html
<a id="JosephsonCircuits.pointmoved!-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.pointmoved!`](api/internals.md#JosephsonCircuits.pointmoved!-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.pointsystem-Tuple{JosephsonCircuits.HBOperatingPoint}"></a>
```

[`JosephsonCircuits.pointsystem`](api/internals.md#JosephsonCircuits.pointsystem-Tuple%7BJosephsonCircuits.HBOperatingPoint%7D)

```@raw html
<a id="JosephsonCircuits.polar-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.polar`](api/internals.md#JosephsonCircuits.polar-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.portdiagonal-Tuple{Number}"></a>
```

[`JosephsonCircuits.portdiagonal`](api/internals.md#JosephsonCircuits.portdiagonal-Tuple%7BNumber%7D)

```@raw html
<a id="JosephsonCircuits.portreferenceimpedances-Tuple{Vector{JosephsonCircuits.CompiledPort}, Any}"></a>
```

[`JosephsonCircuits.portreferenceimpedances`](api/internals.md#JosephsonCircuits.portreferenceimpedances-Tuple%7BVector%7BJosephsonCircuits.CompiledPort%7D%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.ports_modes_to_modes_ports_block-Tuple{AbstractMatrix, Int64}"></a>
```

[`JosephsonCircuits.ports_modes_to_modes_ports_block`](api/internals.md#JosephsonCircuits.ports_modes_to_modes_ports_block-Tuple%7BAbstractMatrix%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.ports_modes_to_modes_ports_pair-Tuple{AbstractMatrix, Int64}"></a>
```

[`JosephsonCircuits.ports_modes_to_modes_ports_pair`](api/internals.md#JosephsonCircuits.ports_modes_to_modes_ports_pair-Tuple%7BAbstractMatrix%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.ports_modes_to_modes_ports_perm-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.ports_modes_to_modes_ports_perm`](api/internals.md#JosephsonCircuits.ports_modes_to_modes_ports_perm-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.ports_modes_to_modes_ports_scattering-Tuple{AbstractMatrix, Int64}"></a>
```

[`JosephsonCircuits.ports_modes_to_modes_ports_scattering`](api/internals.md#JosephsonCircuits.ports_modes_to_modes_ports_scattering-Tuple%7BAbstractMatrix%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.portsolutionrows-Tuple{Any, Any, Integer}"></a>
```

[`JosephsonCircuits.portsolutionrows`](api/internals.md#JosephsonCircuits.portsolutionrows-Tuple%7BAny%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.portsourcecurrents-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.portsourcecurrents`](api/internals.md#JosephsonCircuits.portsourcecurrents-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.portwavescale-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.portwavescale`](api/internals.md#JosephsonCircuits.portwavescale-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.pre_iwasawa_block-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.pre_iwasawa_block`](api/internals.md#JosephsonCircuits.pre_iwasawa_block-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.preconditionedproduct!-Tuple{AbstractVector, AbstractVector, Any, Any, AbstractVector}"></a>
```

[`JosephsonCircuits.preconditionedproduct!`](api/internals.md#JosephsonCircuits.preconditionedproduct!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20Any%2C%20Any%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.preparenoise-Tuple{NoiseCovariance, Int64}"></a>
```

[`JosephsonCircuits.preparenoise`](api/internals.md#JosephsonCircuits.preparenoise-Tuple%7BNoiseCovariance%2C%20Int64%7D)

```@raw html
<a id="JosephsonCircuits.printsymmetries-Tuple{JosephsonCircuits.Frequencies}"></a>
```

[`JosephsonCircuits.printsymmetries`](api/internals.md#JosephsonCircuits.printsymmetries-Tuple%7BJosephsonCircuits.Frequencies%7D)

```@raw html
<a id="JosephsonCircuits.printsymmetries-Union{Tuple{N}, Tuple{NTuple{N, Int64}, NTuple{N, Int64}}} where N"></a>
```

[`JosephsonCircuits.printsymmetries`](api/internals.md#JosephsonCircuits.printsymmetries-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.probecouplings!-Tuple{JosephsonCircuits.ClusterProbe, JosephsonCircuits.HBSystem, Integer}"></a>
```

[`JosephsonCircuits.probecouplings!`](api/internals.md#JosephsonCircuits.probecouplings!-Tuple%7BJosephsonCircuits.ClusterProbe%2C%20JosephsonCircuits.HBSystem%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.provablylossless-Tuple{ScatteringParameters}"></a>
```

[`JosephsonCircuits.provablylossless`](api/internals.md#JosephsonCircuits.provablylossless-Tuple%7BScatteringParameters%7D)

```@raw html
<a id="JosephsonCircuits.providersize-Tuple{JosephsonCircuits.ConstantMatrixProvider}"></a>
```

[`JosephsonCircuits.providersize`](api/internals.md#JosephsonCircuits.providersize-Tuple%7BJosephsonCircuits.ConstantMatrixProvider%7D)

```@raw html
<a id="JosephsonCircuits.psdcholesky!-Tuple{Any, Integer, Integer}"></a>
```

[`JosephsonCircuits.psdcholesky!`](api/internals.md#JosephsonCircuits.psdcholesky!-Tuple%7BAny%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.pumpedblocknoisewaves!-Tuple{AbstractMatrix, LinearizedScattering, JosephsonCircuits.StampedScatteringBlock, AbstractVector, AbstractMatrix, Integer, Integer, AbstractMatrix{Int64}}"></a>
```

[`JosephsonCircuits.pumpedblocknoisewaves!`](api/internals.md#JosephsonCircuits.pumpedblocknoisewaves!-Tuple%7BAbstractMatrix%2C%20LinearizedScattering%2C%20JosephsonCircuits.StampedScatteringBlock%2C%20AbstractVector%2C%20AbstractMatrix%2C%20Integer%2C%20Integer%2C%20AbstractMatrix%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.pumpedfamily-Tuple{LinearizedScattering, Any}"></a>
```

[`JosephsonCircuits.pumpedfamily`](api/internals.md#JosephsonCircuits.pumpedfamily-Tuple%7BLinearizedScattering%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.pumpedharmonics-Tuple{LinearizedScattering, AbstractVector}"></a>
```

[`JosephsonCircuits.pumpedharmonics`](api/internals.md#JosephsonCircuits.pumpedharmonics-Tuple%7BLinearizedScattering%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.pumpednoisematrices-Tuple{LinearizedScattering, AbstractVector, AbstractMatrix{Int64}}"></a>
```

[`JosephsonCircuits.pumpednoisematrices`](api/internals.md#JosephsonCircuits.pumpednoisematrices-Tuple%7BLinearizedScattering%2C%20AbstractVector%2C%20AbstractMatrix%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.pumpmodeset-Union{Tuple{N}, Tuple{NTuple{N, Real}, NTuple{N, Int64}, NTuple{N, Int64}}} where N"></a>
```

[`JosephsonCircuits.pumpmodeset`](api/internals.md#JosephsonCircuits.pumpmodeset-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Real%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20NTuple%7BN%2C%20Int64%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.pumpsourcecurrents-NTuple{5, Any}"></a>
```

[`JosephsonCircuits.pumpsourcecurrents`](api/internals.md#JosephsonCircuits.pumpsourcecurrents-NTuple%7B5%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.quadratic_trial_step-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.quadratic_trial_step`](api/internals.md#JosephsonCircuits.quadratic_trial_step-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.quadrature_to_ladder_block-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.quadrature_to_ladder_block`](api/internals.md#JosephsonCircuits.quadrature_to_ladder_block-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.quadrature_to_ladder_block-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.quadrature_to_ladder_block`](api/internals.md#JosephsonCircuits.quadrature_to_ladder_block-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.quadrature_to_ladder_pair-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.quadrature_to_ladder_pair`](api/internals.md#JosephsonCircuits.quadrature_to_ladder_pair-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.quadrature_to_ladder_pair-Tuple{AbstractVector}"></a>
```

[`JosephsonCircuits.quadrature_to_ladder_pair`](api/internals.md#JosephsonCircuits.quadrature_to_ladder_pair-Tuple%7BAbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.quadrature_to_scattering_block-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Any}} where T&lt;:Union{AbstractFloat, Complex{&lt;:AbstractFloat}}"></a>
```

[`JosephsonCircuits.quadrature_to_scattering_block`](api/internals.md#JosephsonCircuits.quadrature_to_scattering_block-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractMatrix%7BT%7D%2C%20Any%7D%7D%20where%20T%3C%3AUnion%7BAbstractFloat%2C%20Complex%7B%3C%3AAbstractFloat%7D%7D)

```@raw html
<a id="JosephsonCircuits.quadrature_to_scattering_pair-Union{Tuple{T}, Tuple{AbstractMatrix{T}, Any}} where T&lt;:Union{AbstractFloat, Complex{&lt;:AbstractFloat}}"></a>
```

[`JosephsonCircuits.quadrature_to_scattering_pair`](api/internals.md#JosephsonCircuits.quadrature_to_scattering_pair-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractMatrix%7BT%7D%2C%20Any%7D%7D%20where%20T%3C%3AUnion%7BAbstractFloat%2C%20Complex%7B%3C%3AAbstractFloat%7D%7D)

```@raw html
<a id="JosephsonCircuits.quantumnoisemargin-Tuple{AbstractMatrix, AbstractMatrix}"></a>
```

[`JosephsonCircuits.quantumnoisemargin`](api/internals.md#JosephsonCircuits.quantumnoisemargin-Tuple%7BAbstractMatrix%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.rand_bogoliubov_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_bogoliubov_block`](api/internals.md#JosephsonCircuits.rand_bogoliubov_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_bogoliubov_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_bogoliubov_pair`](api/internals.md#JosephsonCircuits.rand_bogoliubov_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_conjugate_symplectic_block-Tuple{Any, Integer}"></a>
```

[`JosephsonCircuits.rand_conjugate_symplectic_block`](api/internals.md#JosephsonCircuits.rand_conjugate_symplectic_block-Tuple%7BAny%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_cptp_ladder_block-Tuple{Any, Integer}"></a>
```

[`JosephsonCircuits.rand_cptp_ladder_block`](api/internals.md#JosephsonCircuits.rand_cptp_ladder_block-Tuple%7BAny%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_cptp_ladder_pair-Tuple{Any, Integer}"></a>
```

[`JosephsonCircuits.rand_cptp_ladder_pair`](api/internals.md#JosephsonCircuits.rand_cptp_ladder_pair-Tuple%7BAny%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_cptp_quadrature_block-Tuple{Any, Integer}"></a>
```

[`JosephsonCircuits.rand_cptp_quadrature_block`](api/internals.md#JosephsonCircuits.rand_cptp_quadrature_block-Tuple%7BAny%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_cptp_quadrature_pair-Tuple{Any, Integer}"></a>
```

[`JosephsonCircuits.rand_cptp_quadrature_pair`](api/internals.md#JosephsonCircuits.rand_cptp_quadrature_pair-Tuple%7BAny%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_orthogonal_bogoliubov_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_orthogonal_bogoliubov_block`](api/internals.md#JosephsonCircuits.rand_orthogonal_bogoliubov_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_orthogonal_bogoliubov_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_orthogonal_bogoliubov_pair`](api/internals.md#JosephsonCircuits.rand_orthogonal_bogoliubov_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_orthogonal_symplectic_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_orthogonal_symplectic_block`](api/internals.md#JosephsonCircuits.rand_orthogonal_symplectic_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_orthogonal_symplectic_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_orthogonal_symplectic_pair`](api/internals.md#JosephsonCircuits.rand_orthogonal_symplectic_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_positive_definite_symplectic_block-Tuple{Any, Integer}"></a>
```

[`JosephsonCircuits.rand_positive_definite_symplectic_block`](api/internals.md#JosephsonCircuits.rand_positive_definite_symplectic_block-Tuple%7BAny%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_positive_definite_symplectic_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_positive_definite_symplectic_block`](api/internals.md#JosephsonCircuits.rand_positive_definite_symplectic_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_positive_definite_symplectic_pair-Tuple{Any, Integer}"></a>
```

[`JosephsonCircuits.rand_positive_definite_symplectic_pair`](api/internals.md#JosephsonCircuits.rand_positive_definite_symplectic_pair-Tuple%7BAny%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_positive_definite_symplectic_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_positive_definite_symplectic_pair`](api/internals.md#JosephsonCircuits.rand_positive_definite_symplectic_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_positive_semi_definite-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.rand_positive_semi_definite`](api/internals.md#JosephsonCircuits.rand_positive_semi_definite-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.rand_pseudo_unitary_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_pseudo_unitary_block`](api/internals.md#JosephsonCircuits.rand_pseudo_unitary_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_pseudo_unitary_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_pseudo_unitary_pair`](api/internals.md#JosephsonCircuits.rand_pseudo_unitary_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_symplectic_block-Tuple{DataType, Integer}"></a>
```

[`JosephsonCircuits.rand_symplectic_block`](api/internals.md#JosephsonCircuits.rand_symplectic_block-Tuple%7BDataType%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.rand_symplectic_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_symplectic_block`](api/internals.md#JosephsonCircuits.rand_symplectic_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.rand_symplectic_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.rand_symplectic_pair`](api/internals.md#JosephsonCircuits.rand_symplectic_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.real_to_complex!-Union{Tuple{T}, Tuple{AbstractArray{Complex{T}, 1}, AbstractVector{T}, AbstractVector{Bool}}} where T&lt;:Real"></a>
```

[`JosephsonCircuits.real_to_complex!`](api/internals.md#JosephsonCircuits.real_to_complex!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BComplex%7BT%7D%2C%201%7D%2C%20AbstractVector%7BT%7D%2C%20AbstractVector%7BBool%7D%7D%7D%20where%20T%3C%3AReal)

```@raw html
<a id="JosephsonCircuits.realblockterm-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.realblockterm`](api/internals.md#JosephsonCircuits.realblockterm-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.realdim-Tuple{Integer, AbstractVector{Bool}}"></a>
```

[`JosephsonCircuits.realdim`](api/internals.md#JosephsonCircuits.realdim-Tuple%7BInteger%2C%20AbstractVector%7BBool%7D%7D)

```@raw html
<a id="JosephsonCircuits.realjacobiancolumnitem!-NTuple{16, Any}"></a>
```

[`JosephsonCircuits.realjacobiancolumnitem!`](api/internals.md#JosephsonCircuits.realjacobiancolumnitem!-NTuple%7B16%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.realjacobianstructure-Union{Tuple{T}, Tuple{Matrix, Matrix, SparseVector, SparseMatrixCSC, Integer, Integer, Any, Any, Any, JosephsonCircuits.ModeLayout}, Tuple{Matrix, Matrix, SparseVector, SparseMatrixCSC, Integer, Integer, Any, Any, Any, JosephsonCircuits.ModeLayout, Type{T}}} where T&lt;:Real"></a>
```

[`JosephsonCircuits.realjacobianstructure`](api/internals.md#JosephsonCircuits.realjacobianstructure-Union%7BTuple%7BT%7D%2C%20Tuple%7BMatrix%2C%20Matrix%2C%20SparseVector%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Any%2C%20Any%2C%20Any%2C%20JosephsonCircuits.ModeLayout%7D%2C%20Tuple%7BMatrix%2C%20Matrix%2C%20SparseVector%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20Any%2C%20Any%2C%20Any%2C%20JosephsonCircuits.ModeLayout%2C%20Type%7BT%7D%7D%7D%20where%20T%3C%3AReal)

```@raw html
<a id="JosephsonCircuits.realjosephsonentry-Union{Tuple{T}, Tuple{Type{T}, Vararg{Any, 11}}} where T"></a>
```

[`JosephsonCircuits.realjosephsonentry`](api/internals.md#JosephsonCircuits.realjosephsonentry-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20Vararg%7BAny%2C%2011%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.realslot-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.realslot`](api/internals.md#JosephsonCircuits.realslot-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.realstructureentry-Union{Tuple{T}, Tuple{Type{T}, Vararg{Any, 14}}} where T"></a>
```

[`JosephsonCircuits.realstructureentry`](api/internals.md#JosephsonCircuits.realstructureentry-Union%7BTuple%7BT%7D%2C%20Tuple%7BType%7BT%7D%2C%20Vararg%7BAny%2C%2014%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.realtocomplexkernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.realtocomplexkernel!`](api/internals.md#JosephsonCircuits.realtocomplexkernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.rebind!-Tuple{JosephsonCircuits.HBSystem, Vararg{Any, 7}}"></a>
```

[`JosephsonCircuits.rebind!`](api/internals.md#JosephsonCircuits.rebind!-Tuple%7BJosephsonCircuits.HBSystem%2C%20Vararg%7BAny%2C%207%7D%7D)

```@raw html
<a id="JosephsonCircuits.rebind!-Tuple{JosephsonCircuits.ModeCouplingPreconditioner, JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.rebind!`](api/internals.md#JosephsonCircuits.rebind!-Tuple%7BJosephsonCircuits.ModeCouplingPreconditioner%2C%20JosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.refactorize!-Tuple{JosephsonCircuits.ModeCouplingPreconditioner}"></a>
```

[`JosephsonCircuits.refactorize!`](api/internals.md#JosephsonCircuits.refactorize!-Tuple%7BJosephsonCircuits.ModeCouplingPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.refill!-Tuple{JosephsonCircuits.PaddedLinearTerm, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, Any, Any, Any}"></a>
```

[`JosephsonCircuits.refill!`](api/internals.md#JosephsonCircuits.refill!-Tuple%7BJosephsonCircuits.PaddedLinearTerm%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20Any%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.refinedsolve!-Tuple{AbstractArray{&lt;:Any, 3}, JosephsonCircuits.SparseBlockFactorization, AbstractMatrix}"></a>
```

[`JosephsonCircuits.refinedsolve!`](api/internals.md#JosephsonCircuits.refinedsolve!-Tuple%7BAbstractArray%7B%3C%3AAny%2C%203%7D%2C%20JosephsonCircuits.SparseBlockFactorization%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.refreshblockstamps!-Tuple{JosephsonCircuits.WorkerBlockSensitivity, Any}"></a>
```

[`JosephsonCircuits.refreshblockstamps!`](api/internals.md#JosephsonCircuits.refreshblockstamps!-Tuple%7BJosephsonCircuits.WorkerBlockSensitivity%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.refreshlinear!-Tuple{Any, JosephsonCircuits.LinearTermGather, Any, Any, Any, Diagonal, Diagonal, Any}"></a>
```

[`JosephsonCircuits.refreshlinear!`](api/internals.md#JosephsonCircuits.refreshlinear!-Tuple%7BAny%2C%20JosephsonCircuits.LinearTermGather%2C%20Any%2C%20Any%2C%20Any%2C%20Diagonal%2C%20Diagonal%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.refreshvalues!-Tuple{JosephsonCircuits.StructureComplexJacobianPlan, Any, Any, Any, Any, Any, SparseVector, Any}"></a>
```

[`JosephsonCircuits.refreshvalues!`](api/internals.md#JosephsonCircuits.refreshvalues!-Tuple%7BJosephsonCircuits.StructureComplexJacobianPlan%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20SparseVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.refreshvalues!-Tuple{JosephsonCircuits.StructureRealJacobianPlan, Any, Any, Any, Any, Any, SparseVector, Any}"></a>
```

[`JosephsonCircuits.refreshvalues!`](api/internals.md#JosephsonCircuits.refreshvalues!-Tuple%7BJosephsonCircuits.StructureRealJacobianPlan%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20Any%2C%20SparseVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.refreshvalues!-Union{Tuple{T}, Tuple{JosephsonCircuits.BlockStructure{T}, JosephsonCircuits.HBSystem}} where T"></a>
```

[`JosephsonCircuits.refreshvalues!`](api/internals.md#JosephsonCircuits.refreshvalues!-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.BlockStructure%7BT%7D%2C%20JosephsonCircuits.HBSystem%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.refreshvalues!-Union{Tuple{T}, Tuple{JosephsonCircuits.JunctionStructure{T}, SparseVector, Any}} where T"></a>
```

[`JosephsonCircuits.refreshvalues!`](api/internals.md#JosephsonCircuits.refreshvalues!-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.JunctionStructure%7BT%7D%2C%20SparseVector%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.refreshvalues!-Union{Tuple{T}, Tuple{Ti}, Tuple{JosephsonCircuits.NonlinearTermPlan{Ti, T}, JosephsonCircuits.ValueMaps, SparseMatrixCSC, SparseVector, Any}} where {Ti, T}"></a>
```

[`JosephsonCircuits.refreshvalues!`](api/internals.md#JosephsonCircuits.refreshvalues!-Union%7BTuple%7BT%7D%2C%20Tuple%7BTi%7D%2C%20Tuple%7BJosephsonCircuits.NonlinearTermPlan%7BTi%2C%20T%7D%2C%20JosephsonCircuits.ValueMaps%2C%20SparseMatrixCSC%2C%20SparseVector%2C%20Any%7D%7D%20where%20%7BTi%2C%20T%7D)

```@raw html
<a id="JosephsonCircuits.relationat-Tuple{JosephsonCircuits.JunctionRelations, Any}"></a>
```

[`JosephsonCircuits.relationat`](api/internals.md#JosephsonCircuits.relationat-Tuple%7BJosephsonCircuits.JunctionRelations%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.relationinto!-Tuple{Any, JosephsonCircuits.JunctionRelations, Any}"></a>
```

[`JosephsonCircuits.relationinto!`](api/internals.md#JosephsonCircuits.relationinto!-Tuple%7BAny%2C%20JosephsonCircuits.JunctionRelations%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.remove_edge!-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.remove_edge!`](api/internals.md#JosephsonCircuits.remove_edge!-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.removeconjfreqs-Tuple{JosephsonCircuits.Frequencies}"></a>
```

[`JosephsonCircuits.removeconjfreqs`](api/internals.md#JosephsonCircuits.removeconjfreqs-Tuple%7BJosephsonCircuits.Frequencies%7D)

```@raw html
<a id="JosephsonCircuits.removefreqs-Union{Tuple{N}, Tuple{JosephsonCircuits.Frequencies{N}, AbstractArray{NTuple{N, Int64}, 1}}} where N"></a>
```

[`JosephsonCircuits.removefreqs`](api/internals.md#JosephsonCircuits.removefreqs-Union%7BTuple%7BN%7D%2C%20Tuple%7BJosephsonCircuits.Frequencies%7BN%7D%2C%20AbstractArray%7BNTuple%7BN%2C%20Int64%7D%2C%201%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.reparameterize-Tuple{JosephsonCircuits.SensitivityStamp, Number, Integer}"></a>
```

[`JosephsonCircuits.reparameterize`](api/internals.md#JosephsonCircuits.reparameterize-Tuple%7BJosephsonCircuits.SensitivityStamp%2C%20Number%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.residual!-Tuple{AbstractVector, JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.residual!`](api/internals.md#JosephsonCircuits.residual!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.residualstalled"></a>
```

[`JosephsonCircuits.residualstalled`](api/internals.md#JosephsonCircuits.residualstalled)

```@raw html
<a id="JosephsonCircuits.resolveautomatic-Tuple{Any, SparseMatrixCSC, Integer, Integer, JosephsonCircuits.ModeLayout, Any, Any}"></a>
```

[`JosephsonCircuits.resolveautomatic`](api/internals.md#JosephsonCircuits.resolveautomatic-Tuple%7BAny%2C%20SparseMatrixCSC%2C%20Integer%2C%20Integer%2C%20JosephsonCircuits.ModeLayout%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.restrictmodecoupling-Tuple{Matrix, AbstractMatrix{Bool}}"></a>
```

[`JosephsonCircuits.restrictmodecoupling`](api/internals.md#JosephsonCircuits.restrictmodecoupling-Tuple%7BMatrix%2C%20AbstractMatrix%7BBool%7D%7D)

```@raw html
<a id="JosephsonCircuits.rootedtree-Tuple{Any}"></a>
```

[`JosephsonCircuits.rootedtree`](api/internals.md#JosephsonCircuits.rootedtree-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.rowpointer-Tuple{JosephsonCircuits.DeviceValuedSparseMatrix{var&quot;#s375&quot;, var&quot;#s374&quot;, V} where {var&quot;#s375&quot;, var&quot;#s374&quot;&lt;:SparseMatrixCSC, V&lt;:AbstractVector{var&quot;#s375&quot;}}}"></a>
```

[`JosephsonCircuits.rowpointer`](api/internals.md#JosephsonCircuits.rowpointer-Tuple%7BJosephsonCircuits.DeviceValuedSparseMatrix%7Bvar%22%23s375%22%2C%20var%22%23s374%22%2C%20V%7D%20where%20%7Bvar%22%23s375%22%2C%20var%22%23s374%22%3C%3ASparseMatrixCSC%2C%20V%3C%3AAbstractVector%7Bvar%22%23s375%22%7D%7D%7D)

```@raw html
<a id="JosephsonCircuits.scalecolumns-Tuple{SparseMatrixCSC, AbstractVector}"></a>
```

[`JosephsonCircuits.scalecolumns`](api/internals.md#JosephsonCircuits.scalecolumns-Tuple%7BSparseMatrixCSC%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.scattercanonical!-Tuple{AbstractVector, AbstractVector, JosephsonCircuits.CompositeLayout}"></a>
```

[`JosephsonCircuits.scattercanonical!`](api/internals.md#JosephsonCircuits.scattercanonical!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20JosephsonCircuits.CompositeLayout%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_block_perm-Tuple{Vector{Int64}}"></a>
```

[`JosephsonCircuits.scattering_to_block_perm`](api/internals.md#JosephsonCircuits.scattering_to_block_perm-Tuple%7BVector%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_ladder_block-Tuple{AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.scattering_to_ladder_block`](api/internals.md#JosephsonCircuits.scattering_to_ladder_block-Tuple%7BAbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_ladder_block-Tuple{AbstractVector, Any}"></a>
```

[`JosephsonCircuits.scattering_to_ladder_block`](api/internals.md#JosephsonCircuits.scattering_to_ladder_block-Tuple%7BAbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_ladder_pair-Tuple{AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.scattering_to_ladder_pair`](api/internals.md#JosephsonCircuits.scattering_to_ladder_pair-Tuple%7BAbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_ladder_pair-Tuple{AbstractVector, Any}"></a>
```

[`JosephsonCircuits.scattering_to_ladder_pair`](api/internals.md#JosephsonCircuits.scattering_to_ladder_pair-Tuple%7BAbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_pair_perm-Tuple{Vector{Int64}}"></a>
```

[`JosephsonCircuits.scattering_to_pair_perm`](api/internals.md#JosephsonCircuits.scattering_to_pair_perm-Tuple%7BVector%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_quadrature_block-Tuple{AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.scattering_to_quadrature_block`](api/internals.md#JosephsonCircuits.scattering_to_quadrature_block-Tuple%7BAbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_quadrature_block-Tuple{AbstractVector, Any}"></a>
```

[`JosephsonCircuits.scattering_to_quadrature_block`](api/internals.md#JosephsonCircuits.scattering_to_quadrature_block-Tuple%7BAbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_quadrature_pair-Tuple{AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.scattering_to_quadrature_pair`](api/internals.md#JosephsonCircuits.scattering_to_quadrature_pair-Tuple%7BAbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scattering_to_quadrature_pair-Tuple{AbstractVector, Any}"></a>
```

[`JosephsonCircuits.scattering_to_quadrature_pair`](api/internals.md#JosephsonCircuits.scattering_to_quadrature_pair-Tuple%7BAbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scatteringblockindex-Tuple{JosephsonCircuits.CompiledCircuit, Any}"></a>
```

[`JosephsonCircuits.scatteringblockindex`](api/internals.md#JosephsonCircuits.scatteringblockindex-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.scatteringlinearterm-Tuple{JosephsonCircuits.ScatteringStampSystem, AbstractVector}"></a>
```

[`JosephsonCircuits.scatteringlinearterm`](api/internals.md#JosephsonCircuits.scatteringlinearterm-Tuple%7BJosephsonCircuits.ScatteringStampSystem%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.scatteringnoisenames-Tuple{JosephsonCircuits.ScatteringNoisePlan, JosephsonCircuits.ScatteringStampSystem}"></a>
```

[`JosephsonCircuits.scatteringnoisenames`](api/internals.md#JosephsonCircuits.scatteringnoisenames-Tuple%7BJosephsonCircuits.ScatteringNoisePlan%2C%20JosephsonCircuits.ScatteringStampSystem%7D)

```@raw html
<a id="JosephsonCircuits.scatteringnoisewaves!"></a>
```

[`JosephsonCircuits.scatteringnoisewaves!`](api/internals.md#JosephsonCircuits.scatteringnoisewaves!)

```@raw html
<a id="JosephsonCircuits.scatteringstampsystem-Tuple{Vector{JosephsonCircuits.CompiledScatteringBlock}, Integer}"></a>
```

[`JosephsonCircuits.scatteringstampsystem`](api/internals.md#JosephsonCircuits.scatteringstampsystem-Tuple%7BVector%7BJosephsonCircuits.CompiledScatteringBlock%7D%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.scatteringstampsystem-Tuple{Vector{JosephsonCircuits.StampedScatteringBlock}, Integer, Integer, Real}"></a>
```

[`JosephsonCircuits.scatteringstampsystem`](api/internals.md#JosephsonCircuits.scatteringstampsystem-Tuple%7BVector%7BJosephsonCircuits.StampedScatteringBlock%7D%2C%20Integer%2C%20Integer%2C%20Real%7D)

```@raw html
<a id="JosephsonCircuits.scatteringvalues!-Tuple{AbstractVector, JosephsonCircuits.ScatteringStampSystem, AbstractVector, JosephsonCircuits.ScatteringWorkspace}"></a>
```

[`JosephsonCircuits.scatteringvalues!`](api/internals.md#JosephsonCircuits.scatteringvalues!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.ScatteringStampSystem%2C%20AbstractVector%2C%20JosephsonCircuits.ScatteringWorkspace%7D)

```@raw html
<a id="JosephsonCircuits.scattervalues!-Tuple{AbstractVector, AbstractVector, AbstractArray}"></a>
```

[`JosephsonCircuits.scattervalues!`](api/internals.md#JosephsonCircuits.scattervalues!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.seeddeflation!-Tuple{JosephsonCircuits.AbstractPreconditioner, AbstractMatrix}"></a>
```

[`JosephsonCircuits.seeddeflation!`](api/internals.md#JosephsonCircuits.seeddeflation!-Tuple%7BJosephsonCircuits.AbstractPreconditioner%2C%20AbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.seedordering!-Tuple{JosephsonCircuits.FactorizationCache, SparseMatrixCSC, Any}"></a>
```

[`JosephsonCircuits.seedordering!`](api/internals.md#JosephsonCircuits.seedordering!-Tuple%7BJosephsonCircuits.FactorizationCache%2C%20SparseMatrixCSC%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.segmentbydest!-Tuple{AbstractVector, AbstractVector, AbstractVector, Any}"></a>
```

[`JosephsonCircuits.segmentbydest!`](api/internals.md#JosephsonCircuits.segmentbydest!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20AbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.selfconjmodes-Tuple{JosephsonCircuits.Frequencies}"></a>
```

[`JosephsonCircuits.selfconjmodes`](api/internals.md#JosephsonCircuits.selfconjmodes-Tuple%7BJosephsonCircuits.Frequencies%7D)

```@raw html
<a id="JosephsonCircuits.sensitivitydim-Tuple{JosephsonCircuits.HBOperatingPoint}"></a>
```

[`JosephsonCircuits.sensitivitydim`](api/internals.md#JosephsonCircuits.sensitivitydim-Tuple%7BJosephsonCircuits.HBOperatingPoint%7D)

```@raw html
<a id="JosephsonCircuits.sensitivityjacobian-Tuple{JosephsonCircuits.HBOperatingPoint}"></a>
```

[`JosephsonCircuits.sensitivityjacobian`](api/internals.md#JosephsonCircuits.sensitivityjacobian-Tuple%7BJosephsonCircuits.HBOperatingPoint%7D)

```@raw html
<a id="JosephsonCircuits.sensitivitypairtable-Tuple{Any}"></a>
```

[`JosephsonCircuits.sensitivitypairtable`](api/internals.md#JosephsonCircuits.sensitivitypairtable-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.sensitivitystampvalue-Tuple{JosephsonCircuits.SensitivityStamp, Integer, Any, Any}"></a>
```

[`JosephsonCircuits.sensitivitystampvalue`](api/internals.md#JosephsonCircuits.sensitivitystampvalue-Tuple%7BJosephsonCircuits.SensitivityStamp%2C%20Integer%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.setfactorization-Tuple{BlockDiagonal, Any}"></a>
```

[`JosephsonCircuits.setfactorization`](api/internals.md#JosephsonCircuits.setfactorization-Tuple%7BBlockDiagonal%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.setpoint!-Tuple{JosephsonCircuits.HBSystem, AbstractVector{&lt;:Complex}}"></a>
```

[`JosephsonCircuits.setpoint!`](api/internals.md#JosephsonCircuits.setpoint!-Tuple%7BJosephsonCircuits.HBSystem%2C%20AbstractVector%7B%3C%3AComplex%7D%7D)

```@raw html
<a id="JosephsonCircuits.setscatteringindexmap!-Tuple{JosephsonCircuits.ScatteringStampSystem, SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.setscatteringindexmap!`](api/internals.md#JosephsonCircuits.setscatteringindexmap!-Tuple%7BJosephsonCircuits.ScatteringStampSystem%2C%20SparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.shortednets-Tuple{ElaboratedCircuit}"></a>
```

[`JosephsonCircuits.shortednets`](api/internals.md#JosephsonCircuits.shortednets-Tuple%7BElaboratedCircuit%7D)

```@raw html
<a id="JosephsonCircuits.showstruct-Tuple{IO, Any}"></a>
```

[`JosephsonCircuits.showstruct`](api/internals.md#JosephsonCircuits.showstruct-Tuple%7BIO%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.sinusoidalmask-Tuple{JosephsonCircuits.JunctionRelations}"></a>
```

[`JosephsonCircuits.sinusoidalmask`](api/internals.md#JosephsonCircuits.sinusoidalmask-Tuple%7BJosephsonCircuits.JunctionRelations%7D)

```@raw html
<a id="JosephsonCircuits.snapscattering-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.snapscattering`](api/internals.md#JosephsonCircuits.snapscattering-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.solveS!-NTuple{24, Any}"></a>
```

[`JosephsonCircuits.solveS!`](api/internals.md#JosephsonCircuits.solveS!-NTuple%7B24%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.solveS_initialize-Tuple{AbstractVector, AbstractVector}"></a>
```

[`JosephsonCircuits.solveS_initialize`](api/internals.md#JosephsonCircuits.solveS_initialize-Tuple%7BAbstractVector%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.solveS_update!-NTuple{10, Any}"></a>
```

[`JosephsonCircuits.solveS_update!`](api/internals.md#JosephsonCircuits.solveS_update!-NTuple%7B10%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.solvebatch!-Tuple{JosephsonCircuits.DeviceSweep, Integer}"></a>
```

[`JosephsonCircuits.solvebatch!`](api/internals.md#JosephsonCircuits.solvebatch!-Tuple%7BJosephsonCircuits.DeviceSweep%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.solveonbackend!-Tuple{Function, AbstractVector, Any, AbstractVector, Any}"></a>
```

[`JosephsonCircuits.solveonbackend!`](api/internals.md#JosephsonCircuits.solveonbackend!-Tuple%7BFunction%2C%20AbstractVector%2C%20Any%2C%20AbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.solvequasinewton!-Tuple{QuasiNewton}"></a>
```

[`JosephsonCircuits.solvequasinewton!`](api/internals.md#JosephsonCircuits.solvequasinewton!-Tuple%7BQuasiNewton%7D)

```@raw html
<a id="JosephsonCircuits.solverkwargs-Tuple{Union{Nothing, JosephsonCircuits.AbstractFactorization}}"></a>
```

[`JosephsonCircuits.solverkwargs`](api/internals.md#JosephsonCircuits.solverkwargs-Tuple%7BUnion%7BNothing%2C%20JosephsonCircuits.AbstractFactorization%7D%7D)

```@raw html
<a id="JosephsonCircuits.solverprecision-Tuple{NewtonKrylov}"></a>
```

[`JosephsonCircuits.solverprecision`](api/internals.md#JosephsonCircuits.solverprecision-Tuple%7BNewtonKrylov%7D)

```@raw html
<a id="JosephsonCircuits.sortnodes-Tuple{Vector{String}, Vector{Int64}, Vector{Int64}}"></a>
```

[`JosephsonCircuits.sortnodes`](api/internals.md#JosephsonCircuits.sortnodes-Tuple%7BVector%7BString%7D%2C%20Vector%7BInt64%7D%2C%20Vector%7BInt64%7D%7D)

```@raw html
<a id="JosephsonCircuits.sourcetable-Union{Tuple{N}, Tuple{Any, NTuple{N, Number}}} where N"></a>
```

[`JosephsonCircuits.sourcetable`](api/internals.md#JosephsonCircuits.sourcetable-Union%7BTuple%7BN%7D%2C%20Tuple%7BAny%2C%20NTuple%7BN%2C%20Number%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.spaddkeepzeros-Tuple{SparseMatrixCSC, SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.spaddkeepzeros`](api/internals.md#JosephsonCircuits.spaddkeepzeros-Tuple%7BSparseMatrixCSC%2C%20SparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.sparseadd!-Tuple{SparseMatrixCSC, Number, SparseMatrixCSC, Vector}"></a>
```

[`JosephsonCircuits.sparseadd!`](api/internals.md#JosephsonCircuits.sparseadd!-Tuple%7BSparseMatrixCSC%2C%20Number%2C%20SparseMatrixCSC%2C%20Vector%7D)

```@raw html
<a id="JosephsonCircuits.sparseaddconjsubst!-Tuple{SparseMatrixCSC, Number, SparseMatrixCSC, Any, AbstractVector, Integer}"></a>
```

[`JosephsonCircuits.sparseaddconjsubst!`](api/internals.md#JosephsonCircuits.sparseaddconjsubst!-Tuple%7BSparseMatrixCSC%2C%20Number%2C%20SparseMatrixCSC%2C%20Any%2C%20AbstractVector%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.sparseaddmap-Tuple{SparseMatrixCSC, SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.sparseaddmap`](api/internals.md#JosephsonCircuits.sparseaddmap-Tuple%7BSparseMatrixCSC%2C%20SparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.sparsefactorbytes-Union{Tuple{T}, Tuple{SparseMatrixCSC, Type{T}}, Tuple{SparseMatrixCSC, Type{T}, Union{Nothing, JosephsonCircuits.FillOrdering}}} where T"></a>
```

[`JosephsonCircuits.sparsefactorbytes`](api/internals.md#JosephsonCircuits.sparsefactorbytes-Union%7BTuple%7BT%7D%2C%20Tuple%7BSparseMatrixCSC%2C%20Type%7BT%7D%7D%2C%20Tuple%7BSparseMatrixCSC%2C%20Type%7BT%7D%2C%20Union%7BNothing%2C%20JosephsonCircuits.FillOrdering%7D%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.spectralclusters-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.spectralclusters`](api/internals.md#JosephsonCircuits.spectralclusters-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.spice_hb_load-Tuple{Any}"></a>
```

[`JosephsonCircuits.spice_hb_load`](api/internals.md#JosephsonCircuits.spice_hb_load-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.spice_raw_load-Tuple{Any}"></a>
```

[`JosephsonCircuits.spice_raw_load`](api/internals.md#JosephsonCircuits.spice_raw_load-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.spice_run-Tuple{AbstractVector, Any}"></a>
```

[`JosephsonCircuits.spice_run`](api/internals.md#JosephsonCircuits.spice_run-Tuple%7BAbstractVector%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.spice_run-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.spice_run`](api/internals.md#JosephsonCircuits.spice_run-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.spicebranches-Tuple{JosephsonCircuits.CompiledCircuit, AbstractVector}"></a>
```

[`JosephsonCircuits.spicebranches`](api/internals.md#JosephsonCircuits.spicebranches-Tuple%7BJosephsonCircuits.CompiledCircuit%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.splitfrequencydependent-Tuple{SparseMatrixCSC}"></a>
```

[`JosephsonCircuits.splitfrequencydependent`](api/internals.md#JosephsonCircuits.splitfrequencydependent-Tuple%7BSparseMatrixCSC%7D)

```@raw html
<a id="JosephsonCircuits.sprandsubset"></a>
```

[`JosephsonCircuits.sprandsubset`](api/internals.md#JosephsonCircuits.sprandsubset)

```@raw html
<a id="JosephsonCircuits.stagedeviceproviders!-Tuple{AbstractMatrix, JosephsonCircuits.DeviceProviders, Any, Integer, Integer}"></a>
```

[`JosephsonCircuits.stagedeviceproviders!`](api/internals.md#JosephsonCircuits.stagedeviceproviders!-Tuple%7BAbstractMatrix%2C%20JosephsonCircuits.DeviceProviders%2C%20Any%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.stagedhbnlsolve-Union{Tuple{N}, Tuple{Staged, NTuple{N, Float64}, NTuple{N, Int64}, Array{@NamedTuple{mode::NTuple{N, Int64}, port::Int64, current::ComplexF64}, 1}, JosephsonCircuits.CompiledCircuit, Dict{Any, Any}}} where N"></a>
```

[`JosephsonCircuits.stagedhbnlsolve`](api/internals.md#JosephsonCircuits.stagedhbnlsolve-Union%7BTuple%7BN%7D%2C%20Tuple%7BStaged%2C%20NTuple%7BN%2C%20Float64%7D%2C%20NTuple%7BN%2C%20Int64%7D%2C%20Array%7B%40NamedTuple%7Bmode%3A%3ANTuple%7BN%2C%20Int64%7D%2C%20port%3A%3AInt64%2C%20current%3A%3AComplexF64%7D%2C%201%7D%2C%20JosephsonCircuits.CompiledCircuit%2C%20Dict%7BAny%2C%20Any%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.stagescatteringstamps!-Tuple{JosephsonCircuits.DeviceScatteringStamps, Any, Integer, Integer, Any}"></a>
```

[`JosephsonCircuits.stagescatteringstamps!`](api/internals.md#JosephsonCircuits.stagescatteringstamps!-Tuple%7BJosephsonCircuits.DeviceScatteringStamps%2C%20Any%2C%20Integer%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.stageweights!-Tuple{JosephsonCircuits.RationalWork, JosephsonCircuits.TransientSystem, Any, Any}"></a>
```

[`JosephsonCircuits.stageweights!`](api/internals.md#JosephsonCircuits.stageweights!-Tuple%7BJosephsonCircuits.RationalWork%2C%20JosephsonCircuits.TransientSystem%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.stalled!-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.stalled!`](api/internals.md#JosephsonCircuits.stalled!-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.stallmessage-Tuple{Symbol}"></a>
```

[`JosephsonCircuits.stallmessage`](api/internals.md#JosephsonCircuits.stallmessage-Tuple%7BSymbol%7D)

```@raw html
<a id="JosephsonCircuits.statednoise-Tuple{ScatteringParameters}"></a>
```

[`JosephsonCircuits.statednoise`](api/internals.md#JosephsonCircuits.statednoise-Tuple%7BScatteringParameters%7D)

```@raw html
<a id="JosephsonCircuits.statednoisefactors!-Tuple{Any, Any, Any, Integer, Integer}"></a>
```

[`JosephsonCircuits.statednoisefactors!`](api/internals.md#JosephsonCircuits.statednoisefactors!-Tuple%7BAny%2C%20Any%2C%20Any%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.structureassemblerowkernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.structureassemblerowkernel!`](api/internals.md#JosephsonCircuits.structureassemblerowkernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.structureassemblykernel!-Tuple{Any}"></a>
```

[`JosephsonCircuits.structureassemblykernel!`](api/internals.md#JosephsonCircuits.structureassemblykernel!-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.structurejacobian-Tuple{Any, Matrix, Matrix, Vararg{Any, 10}}"></a>
```

[`JosephsonCircuits.structurejacobian`](api/internals.md#JosephsonCircuits.structurejacobian-Tuple%7BAny%2C%20Matrix%2C%20Matrix%2C%20Vararg%7BAny%2C%2010%7D%7D)

```@raw html
<a id="JosephsonCircuits.substitute!-Union{Tuple{T}, Tuple{AbstractArray{&lt;:Any, 3}, JosephsonCircuits.BlockLU{T}, AbstractArray{&lt;:Any, 3}, AbstractArray{&lt;:Any, 3}, Any}} where T"></a>
```

[`JosephsonCircuits.substitute!`](api/internals.md#JosephsonCircuits.substitute!-Union%7BTuple%7BT%7D%2C%20Tuple%7BAbstractArray%7B%3C%3AAny%2C%203%7D%2C%20JosephsonCircuits.BlockLU%7BT%7D%2C%20AbstractArray%7B%3C%3AAny%2C%203%7D%2C%20AbstractArray%7B%3C%3AAny%2C%203%7D%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.substitutefreq-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.substitutefreq`](api/internals.md#JosephsonCircuits.substitutefreq-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.sumvalues-Tuple{Symbol, Any, Any}"></a>
```

[`JosephsonCircuits.sumvalues`](api/internals.md#JosephsonCircuits.sumvalues-Tuple%7BSymbol%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.supportsrecycling-Tuple{JosephsonCircuits.AbstractHBLinearSolver}"></a>
```

[`JosephsonCircuits.supportsrecycling`](api/internals.md#JosephsonCircuits.supportsrecycling-Tuple%7BJosephsonCircuits.AbstractHBLinearSolver%7D)

```@raw html
<a id="JosephsonCircuits.sweepdestinations-Tuple{SparseMatrixCSC, AbstractVector, Bool}"></a>
```

[`JosephsonCircuits.sweepdestinations`](api/internals.md#JosephsonCircuits.sweepdestinations-Tuple%7BSparseMatrixCSC%2C%20AbstractVector%2C%20Bool%7D)

```@raw html
<a id="JosephsonCircuits.sweepfrequencies-Tuple{Any}"></a>
```

[`JosephsonCircuits.sweepfrequencies`](api/internals.md#JosephsonCircuits.sweepfrequencies-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.sweepoutputs!-NTuple{16, Any}"></a>
```

[`JosephsonCircuits.sweepoutputs!`](api/internals.md#JosephsonCircuits.sweepoutputs!-NTuple%7B16%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.symbolicfill-Tuple{SparseMatrixCSC, AbstractVector{&lt;:Integer}}"></a>
```

[`JosephsonCircuits.symbolicfill`](api/internals.md#JosephsonCircuits.symbolicfill-Tuple%7BSparseMatrixCSC%2C%20AbstractVector%7B%3C%3AInteger%7D%7D)

```@raw html
<a id="JosephsonCircuits.symbolicindices-Tuple{Any}"></a>
```

[`JosephsonCircuits.symbolicindices`](api/internals.md#JosephsonCircuits.symbolicindices-Tuple%7BAny%7D)

```@raw html
<a id="JosephsonCircuits.symplectic_form_block-Tuple{Integer}"></a>
```

[`JosephsonCircuits.symplectic_form_block`](api/internals.md#JosephsonCircuits.symplectic_form_block-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.symplectic_form_pair-Tuple{Integer}"></a>
```

[`JosephsonCircuits.symplectic_form_pair`](api/internals.md#JosephsonCircuits.symplectic_form_pair-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.symplectic_normal_form_pair-Tuple{AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.symplectic_normal_form_pair`](api/internals.md#JosephsonCircuits.symplectic_normal_form_pair-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.symplecticform-Tuple{Integer}"></a>
```

[`JosephsonCircuits.symplecticform`](api/internals.md#JosephsonCircuits.symplecticform-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.targetstampsystems-Tuple{Any, Integer, Any}"></a>
```

[`JosephsonCircuits.targetstampsystems`](api/internals.md#JosephsonCircuits.targetstampsystems-Tuple%7BAny%2C%20Integer%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.testshow-Tuple{IO, AbstractSparseVector}"></a>
```

[`JosephsonCircuits.testshow`](api/internals.md#JosephsonCircuits.testshow-Tuple%7BIO%2C%20AbstractSparseVector%7D)

```@raw html
<a id="JosephsonCircuits.thermalnoise!-Tuple{AbstractVector, Nothing, Any, Integer}"></a>
```

[`JosephsonCircuits.thermalnoise!`](api/internals.md#JosephsonCircuits.thermalnoise!-Tuple%7BAbstractVector%2C%20Nothing%2C%20Any%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.tobackend-Tuple{KernelAbstractions.Backend, AbstractArray}"></a>
```

[`JosephsonCircuits.tobackend`](api/internals.md#JosephsonCircuits.tobackend-Tuple%7BKernelAbstractions.Backend%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.tohost-Tuple{Array}"></a>
```

[`JosephsonCircuits.tohost`](api/internals.md#JosephsonCircuits.tohost-Tuple%7BArray%7D)

```@raw html
<a id="JosephsonCircuits.tonefrequencies-Union{Tuple{NTuple{N, Number}}, Tuple{N}} where N"></a>
```

[`JosephsonCircuits.tonefrequencies`](api/internals.md#JosephsonCircuits.tonefrequencies-Union%7BTuple%7BNTuple%7BN%2C%20Number%7D%7D%2C%20Tuple%7BN%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.torelations-Tuple{Nothing, AbstractArray}"></a>
```

[`JosephsonCircuits.torelations`](api/internals.md#JosephsonCircuits.torelations-Tuple%7BNothing%2C%20AbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.tracerestart!-Tuple{JosephsonCircuits.NewtonTrace, Any}"></a>
```

[`JosephsonCircuits.tracerestart!`](api/internals.md#JosephsonCircuits.tracerestart!-Tuple%7BJosephsonCircuits.NewtonTrace%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.tracestalled-Tuple{JosephsonCircuits.NewtonTrace, Integer}"></a>
```

[`JosephsonCircuits.tracestalled`](api/internals.md#JosephsonCircuits.tracestalled-Tuple%7BJosephsonCircuits.NewtonTrace%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.tracestart!-Union{Tuple{T}, Tuple{JosephsonCircuits.NewtonTrace{T}, Any, Any, Any}} where T"></a>
```

[`JosephsonCircuits.tracestart!`](api/internals.md#JosephsonCircuits.tracestart!-Union%7BTuple%7BT%7D%2C%20Tuple%7BJosephsonCircuits.NewtonTrace%7BT%7D%2C%20Any%2C%20Any%2C%20Any%7D%7D%20where%20T)

```@raw html
<a id="JosephsonCircuits.tracestep!-Tuple{JosephsonCircuits.NewtonTrace, Any, Bool}"></a>
```

[`JosephsonCircuits.tracestep!`](api/internals.md#JosephsonCircuits.tracestep!-Tuple%7BJosephsonCircuits.NewtonTrace%2C%20Any%2C%20Bool%7D)

```@raw html
<a id="JosephsonCircuits.tracetrial!"></a>
```

[`JosephsonCircuits.tracetrial!`](api/internals.md#JosephsonCircuits.tracetrial!)

```@raw html
<a id="JosephsonCircuits.transientnoiseaccumulate!-NTuple{4, Any}"></a>
```

[`JosephsonCircuits.transientnoiseaccumulate!`](api/internals.md#JosephsonCircuits.transientnoiseaccumulate!-NTuple%7B4%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transientsystem-Tuple{JosephsonCircuits.TransientProblem, Real, JosephsonCircuits.AbstractTransientIntegrator, KernelAbstractions.Backend, JosephsonCircuits.AbstractFactorization}"></a>
```

[`JosephsonCircuits.transientsystem`](api/internals.md#JosephsonCircuits.transientsystem-Tuple%7BJosephsonCircuits.TransientProblem%2C%20Real%2C%20JosephsonCircuits.AbstractTransientIntegrator%2C%20KernelAbstractions.Backend%2C%20JosephsonCircuits.AbstractFactorization%7D)

```@raw html
<a id="JosephsonCircuits.transportcurrent!-Tuple{AbstractVector, JosephsonCircuits.TransportRows, AbstractVector}"></a>
```

[`JosephsonCircuits.transportcurrent!`](api/internals.md#JosephsonCircuits.transportcurrent!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.TransportRows%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.transportresidual!-Tuple{AbstractVector, JosephsonCircuits.TransportRows, AbstractVector}"></a>
```

[`JosephsonCircuits.transportresidual!`](api/internals.md#JosephsonCircuits.transportresidual!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.TransportRows%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.transportrows-Tuple{JosephsonCircuits.DCConductancePlan, AbstractVector, Integer}"></a>
```

[`JosephsonCircuits.transportrows`](api/internals.md#JosephsonCircuits.transportrows-Tuple%7BJosephsonCircuits.DCConductancePlan%2C%20AbstractVector%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.transposedestinations-Tuple{JosephsonCircuits.DeviceScatteringStamps, AbstractVector{Int64}, Any}"></a>
```

[`JosephsonCircuits.transposedestinations`](api/internals.md#JosephsonCircuits.transposedestinations-Tuple%7BJosephsonCircuits.DeviceScatteringStamps%2C%20AbstractVector%7BInt64%7D%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.transposepattern-Union{Tuple{Ti}, Tuple{JosephsonCircuits.DeviceSparsePattern{Ti, V} where V&lt;:AbstractVector{Ti}, Any}} where Ti"></a>
```

[`JosephsonCircuits.transposepattern`](api/internals.md#JosephsonCircuits.transposepattern-Union%7BTuple%7BTi%7D%2C%20Tuple%7BJosephsonCircuits.DeviceSparsePattern%7BTi%2C%20V%7D%20where%20V%3C%3AAbstractVector%7BTi%7D%2C%20Any%7D%7D%20where%20Ti)

```@raw html
<a id="JosephsonCircuits.treepath-Tuple{Vector{Int64}, Vector{Int64}, Integer, Integer}"></a>
```

[`JosephsonCircuits.treepath`](api/internals.md#JosephsonCircuits.treepath-Tuple%7BVector%7BInt64%7D%2C%20Vector%7BInt64%7D%2C%20Integer%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.truncfreqs-Union{Tuple{JosephsonCircuits.Frequencies{N}}, Tuple{N}} where N"></a>
```

[`JosephsonCircuits.truncfreqs`](api/internals.md#JosephsonCircuits.truncfreqs-Union%7BTuple%7BJosephsonCircuits.Frequencies%7BN%7D%7D%2C%20Tuple%7BN%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.tryfactorize!-Tuple{JosephsonCircuits.FactorizationCache, JosephsonCircuits.AbstractFactorization, Any}"></a>
```

[`JosephsonCircuits.tryfactorize!`](api/internals.md#JosephsonCircuits.tryfactorize!-Tuple%7BJosephsonCircuits.FactorizationCache%2C%20JosephsonCircuits.AbstractFactorization%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.trysolve!-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.trysolve!`](api/internals.md#JosephsonCircuits.trysolve!-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.trysolvetranspose!-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.trysolvetranspose!`](api/internals.md#JosephsonCircuits.trysolvetranspose!-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.tuple2edge-Tuple{Vector{Tuple{Int64, Int64}}}"></a>
```

[`JosephsonCircuits.tuple2edge`](api/internals.md#JosephsonCircuits.tuple2edge-Tuple%7BVector%7BTuple%7BInt64%2C%20Int64%7D%7D%7D)

```@raw html
<a id="JosephsonCircuits.uniformbatchlimit-Tuple{Integer}"></a>
```

[`JosephsonCircuits.uniformbatchlimit`](api/internals.md#JosephsonCircuits.uniformbatchlimit-Tuple%7BInteger%7D)

```@raw html
<a id="JosephsonCircuits.unitaritybound-Tuple{JosephsonCircuits.AbstractMatrixProvider}"></a>
```

[`JosephsonCircuits.unitaritybound`](api/internals.md#JosephsonCircuits.unitaritybound-Tuple%7BJosephsonCircuits.AbstractMatrixProvider%7D)

```@raw html
<a id="JosephsonCircuits.unitaritydeviation-Tuple{AbstractMatrix}"></a>
```

[`JosephsonCircuits.unitaritydeviation`](api/internals.md#JosephsonCircuits.unitaritydeviation-Tuple%7BAbstractMatrix%7D)

```@raw html
<a id="JosephsonCircuits.unscalesolution!-Tuple{AbstractArray{&lt;:Any, 3}, JosephsonCircuits.SweepEquilibration, Any}"></a>
```

[`JosephsonCircuits.unscalesolution!`](api/internals.md#JosephsonCircuits.unscalesolution!-Tuple%7BAbstractArray%7B%3C%3AAny%2C%203%7D%2C%20JosephsonCircuits.SweepEquilibration%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.unwrap!-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.unwrap!`](api/internals.md#JosephsonCircuits.unwrap!-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.unwrap!-Union{Tuple{N}, Tuple{T}, Tuple{AbstractArray{T, N}, AbstractArray{T, N}}} where {T&lt;:Real, N}"></a>
```

[`JosephsonCircuits.unwrap!`](api/internals.md#JosephsonCircuits.unwrap!-Union%7BTuple%7BN%7D%2C%20Tuple%7BT%7D%2C%20Tuple%7BAbstractArray%7BT%2C%20N%7D%2C%20AbstractArray%7BT%2C%20N%7D%7D%7D%20where%20%7BT%3C%3AReal%2C%20N%7D)

```@raw html
<a id="JosephsonCircuits.unwrap-Tuple{AbstractArray}"></a>
```

[`JosephsonCircuits.unwrap`](api/internals.md#JosephsonCircuits.unwrap-Tuple%7BAbstractArray%7D)

```@raw html
<a id="JosephsonCircuits.updatepreconditioner!"></a>
```

[`JosephsonCircuits.updatepreconditioner!`](api/internals.md#JosephsonCircuits.updatepreconditioner!)

```@raw html
<a id="JosephsonCircuits.updatepreconditioner!-Tuple{JosephsonCircuits.ModeCouplingPreconditioner, AbstractVector}"></a>
```

[`JosephsonCircuits.updatepreconditioner!`](api/internals.md#JosephsonCircuits.updatepreconditioner!-Tuple%7BJosephsonCircuits.ModeCouplingPreconditioner%2C%20AbstractVector%7D)

```@raw html
<a id="JosephsonCircuits.usescycleharvest-Tuple{JosephsonCircuits.AbstractPreconditioner}"></a>
```

[`JosephsonCircuits.usescycleharvest`](api/internals.md#JosephsonCircuits.usescycleharvest-Tuple%7BJosephsonCircuits.AbstractPreconditioner%7D)

```@raw html
<a id="JosephsonCircuits.valuemaps-Tuple{JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.valuemaps`](api/internals.md#JosephsonCircuits.valuemaps-Tuple%7BJosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.valuemaps-Union{Tuple{T}, Tuple{Ti}, Tuple{JosephsonCircuits.NonlinearTermPlan{Ti, T}, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, SparseMatrixCSC, SparseVector, JosephsonCircuits.ModeLayout, Vector{Int64}}} where {Ti, T}"></a>
```

[`JosephsonCircuits.valuemaps`](api/internals.md#JosephsonCircuits.valuemaps-Union%7BTuple%7BT%7D%2C%20Tuple%7BTi%7D%2C%20Tuple%7BJosephsonCircuits.NonlinearTermPlan%7BTi%2C%20T%7D%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseMatrixCSC%2C%20SparseVector%2C%20JosephsonCircuits.ModeLayout%2C%20Vector%7BInt64%7D%7D%7D%20where%20%7BTi%2C%20T%7D)

```@raw html
<a id="JosephsonCircuits.valuetonumber-Tuple{Any, Any}"></a>
```

[`JosephsonCircuits.valuetonumber`](api/internals.md#JosephsonCircuits.valuetonumber-Tuple%7BAny%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.valuetonumber-Tuple{JosephsonCircuits.CircuitValues.CircuitValue, Any}"></a>
```

[`JosephsonCircuits.valuetonumber`](api/internals.md#JosephsonCircuits.valuetonumber-Tuple%7BJosephsonCircuits.CircuitValues.CircuitValue%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.valuetonumber-Tuple{String, Any}"></a>
```

[`JosephsonCircuits.valuetonumber`](api/internals.md#JosephsonCircuits.valuetonumber-Tuple%7BString%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.valuetonumber-Tuple{Symbol, Any}"></a>
```

[`JosephsonCircuits.valuetonumber`](api/internals.md#JosephsonCircuits.valuetonumber-Tuple%7BSymbol%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.visualizefreqs-Union{Tuple{N}, Tuple{NTuple{N, Any}, JosephsonCircuits.Frequencies{N}}} where N"></a>
```

[`JosephsonCircuits.visualizefreqs`](api/internals.md#JosephsonCircuits.visualizefreqs-Union%7BTuple%7BN%7D%2C%20Tuple%7BNTuple%7BN%2C%20Any%7D%2C%20JosephsonCircuits.Frequencies%7BN%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.warnduplicatematchedload-NTuple{5, Any}"></a>
```

[`JosephsonCircuits.warnduplicatematchedload`](api/internals.md#JosephsonCircuits.warnduplicatematchedload-NTuple%7B5%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.weightedrowpower!-Tuple{AbstractVector, AbstractVector, AbstractMatrix, Any}"></a>
```

[`JosephsonCircuits.weightedrowpower!`](api/internals.md#JosephsonCircuits.weightedrowpower!-Tuple%7BAbstractVector%2C%20AbstractVector%2C%20AbstractMatrix%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.williamson_block-Tuple{AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.williamson_block`](api/internals.md#JosephsonCircuits.williamson_block-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.williamson_pair-Tuple{AbstractMatrix{&lt;:Real}}"></a>
```

[`JosephsonCircuits.williamson_pair`](api/internals.md#JosephsonCircuits.williamson_pair-Tuple%7BAbstractMatrix%7B%3C%3AReal%7D%7D)

```@raw html
<a id="JosephsonCircuits.windowindex-Tuple{JosephsonCircuits.CompositeLayout, Integer}"></a>
```

[`JosephsonCircuits.windowindex`](api/internals.md#JosephsonCircuits.windowindex-Tuple%7BJosephsonCircuits.CompositeLayout%2C%20Integer%7D)

```@raw html
<a id="JosephsonCircuits.windowindices-Tuple{JosephsonCircuits.CompositeLayout}"></a>
```

[`JosephsonCircuits.windowindices`](api/internals.md#JosephsonCircuits.windowindices-Tuple%7BJosephsonCircuits.CompositeLayout%7D)

```@raw html
<a id="JosephsonCircuits.with-Tuple{Union{JosephsonCircuits.IterationInfo, JosephsonCircuits.KrylovSolveInfo}}"></a>
```

[`JosephsonCircuits.with`](api/internals.md#JosephsonCircuits.with-Tuple%7BUnion%7BJosephsonCircuits.IterationInfo%2C%20JosephsonCircuits.KrylovSolveInfo%7D%7D)

```@raw html
<a id="JosephsonCircuits.withescalation-Tuple{NewtonKrylov, Bool}"></a>
```

[`JosephsonCircuits.withescalation`](api/internals.md#JosephsonCircuits.withescalation-Tuple%7BNewtonKrylov%2C%20Bool%7D)

```@raw html
<a id="JosephsonCircuits.withfactorization-Tuple{JosephsonCircuits.AbstractModeCoupling, Any}"></a>
```

[`JosephsonCircuits.withfactorization`](api/internals.md#JosephsonCircuits.withfactorization-Tuple%7BJosephsonCircuits.AbstractModeCoupling%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.withfactors-Tuple{JosephsonCircuits.DeviceBlockNoisePlan}"></a>
```

[`JosephsonCircuits.withfactors`](api/internals.md#JosephsonCircuits.withfactors-Tuple%7BJosephsonCircuits.DeviceBlockNoisePlan%7D)

```@raw html
<a id="JosephsonCircuits.withprecision-Tuple{JosephsonCircuits.AbstractFactorization, Type{&lt;:AbstractFloat}}"></a>
```

[`JosephsonCircuits.withprecision`](api/internals.md#JosephsonCircuits.withprecision-Tuple%7BJosephsonCircuits.AbstractFactorization%2C%20Type%7B%3C%3AAbstractFloat%7D%7D)

```@raw html
<a id="JosephsonCircuits.wmatrix!-Union{Tuple{N}, Tuple{AbstractMatrix, AbstractVector, NTuple{N, Number}, AbstractArray{NTuple{N, Int64}, 1}}} where N"></a>
```

[`JosephsonCircuits.wmatrix!`](api/internals.md#JosephsonCircuits.wmatrix!-Union%7BTuple%7BN%7D%2C%20Tuple%7BAbstractMatrix%2C%20AbstractVector%2C%20NTuple%7BN%2C%20Number%7D%2C%20AbstractArray%7BNTuple%7BN%2C%20Int64%7D%2C%201%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.wmatrix-Union{Tuple{N}, Tuple{AbstractVector, NTuple{N, Number}, AbstractArray{NTuple{N, Int64}, 1}}} where N"></a>
```

[`JosephsonCircuits.wmatrix`](api/internals.md#JosephsonCircuits.wmatrix-Union%7BTuple%7BN%7D%2C%20Tuple%7BAbstractVector%2C%20NTuple%7BN%2C%20Number%7D%2C%20AbstractArray%7BNTuple%7BN%2C%20Int64%7D%2C%201%7D%7D%7D%20where%20N)

```@raw html
<a id="JosephsonCircuits.workspacetwin-Tuple{JosephsonCircuits.HBSystem}"></a>
```

[`JosephsonCircuits.workspacetwin`](api/internals.md#JosephsonCircuits.workspacetwin-Tuple%7BJosephsonCircuits.HBSystem%7D)

```@raw html
<a id="JosephsonCircuits.wrspice_calcS_paramp-Tuple{Any, Any, Any}"></a>
```

[`JosephsonCircuits.wrspice_calcS_paramp`](api/internals.md#JosephsonCircuits.wrspice_calcS_paramp-Tuple%7BAny%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.wrspice_cmd-Tuple{}"></a>
```

[`JosephsonCircuits.wrspice_cmd`](api/internals.md#JosephsonCircuits.wrspice_cmd-Tuple%7B%7D)

```@raw html
<a id="JosephsonCircuits.wrspice_input_ac-Tuple{String, AbstractVector{Float64}, Any, Any}"></a>
```

[`JosephsonCircuits.wrspice_input_ac`](api/internals.md#JosephsonCircuits.wrspice_input_ac-Tuple%7BString%2C%20AbstractVector%7BFloat64%7D%2C%20Any%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.wrspice_input_paramp-NTuple{6, Any}"></a>
```

[`JosephsonCircuits.wrspice_input_paramp`](api/internals.md#JosephsonCircuits.wrspice_input_paramp-NTuple%7B6%2C%20Any%7D)

```@raw html
<a id="JosephsonCircuits.wrspice_input_transient-Tuple{String, Vararg{Any, 7}}"></a>
```

[`JosephsonCircuits.wrspice_input_transient`](api/internals.md#JosephsonCircuits.wrspice_input_transient-Tuple%7BString%2C%20Vararg%7BAny%2C%207%7D%7D)

```@raw html
<a id="LinearAlgebra.dot-Tuple{AbstractVector{&lt;:Union{Float32, Float64, ComplexF64, ComplexF32}}, JosephsonCircuits.DeviceValuedSparseMatrix, AbstractVector}"></a>
```

[`LinearAlgebra.dot`](api/internals.md#LinearAlgebra.dot-Tuple%7BAbstractVector%7B%3C%3AUnion%7BFloat32%2C%20Float64%2C%20ComplexF64%2C%20ComplexF32%7D%7D%2C%20JosephsonCircuits.DeviceValuedSparseMatrix%2C%20AbstractVector%7D)

```@raw html
<a id="LinearAlgebra.ldiv!-Tuple{AbstractVector, JosephsonCircuits.AbstractPreconditioner, AbstractVector}"></a>
```

[`LinearAlgebra.ldiv!`](api/internals.md#LinearAlgebra.ldiv!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.AbstractPreconditioner%2C%20AbstractVector%7D)

```@raw html
<a id="LinearAlgebra.mul!-Tuple{AbstractVector, JosephsonCircuits.AbstractPreconditioner, AbstractVector}"></a>
```

[`LinearAlgebra.mul!`](api/internals.md#LinearAlgebra.mul!-Tuple%7BAbstractVector%2C%20JosephsonCircuits.AbstractPreconditioner%2C%20AbstractVector%7D)

```@raw html
<a id="LinearAlgebra.mul!-Tuple{AbstractVector{&lt;:Union{Float32, Float64, ComplexF64, ComplexF32}}, JosephsonCircuits.DeviceValuedSparseMatrix, AbstractVector}"></a>
```

[`LinearAlgebra.mul!`](api/internals.md#LinearAlgebra.mul!-Tuple%7BAbstractVector%7B%3C%3AUnion%7BFloat32%2C%20Float64%2C%20ComplexF64%2C%20ComplexF32%7D%7D%2C%20JosephsonCircuits.DeviceValuedSparseMatrix%2C%20AbstractVector%7D)
