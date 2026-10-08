# Network parameters

Convert a network's scattering parameters to and from its impedance,
admittance, chain (ABCD) and transfer parameters, at the reference
impedances of your choice. A conversion takes one matrix, or an array of
matrices with the frequency along the dimensions after the first two; an
in-place form, ending in `!`, writes into an array of the same shape.

See the [API overview](../reference.md) for other topics,
[closed-form networks](networkmodels.md) for networks to convert, and
[connections](connections.md) to join them.

## Convert and renormalize

A line of 40 Ω, 1 cm long, between 50 Ω ports, the measurement system.
Its chain matrix gives its scattering parameters at 50 Ω, and its
impedance matrix carries them to 40 Ω ports, where the line is matched.
At 50 Ω the reflection is that of a line section between two equal media,
`Γ(1 - exp(-2iθ))/(1 - Γ² exp(-2iθ))` with `Γ = (40 - 50)/(40 + 50)` and
`θ` its electrical length; at 6 GHz the line is half a wavelength long
and reflects nothing.

```@example networkparameters
using JosephsonCircuits
w = 2pi .* [4e9, 5e9, 6e9]
len, vp = 0.01, 1.2e8
θ = w .* len ./ vp
# one chain matrix per frequency, 2×2×3
abcd = JosephsonCircuits.ABCD_tline(40.0, θ)
S = JosephsonCircuits.ABCDtoS(abcd; portimpedances = 50.0)
Z = JosephsonCircuits.StoZ(S; portimpedances = 50.0)
S40 = JosephsonCircuits.ZtoS(Z; portimpedances = 40.0)
Γ = (40 - 50)/(40 + 50)
closed = Γ .* (1 .- exp.(-2im .* θ)) ./ (1 .- Γ^2 .* exp.(-2im .* θ))
@assert isapprox(S[1, 1, :], closed; atol = 1e-12)
@assert maximum(abs, S40[1, 1, :]) < 1e-12
(S11 = S[1, 1, :], closed = closed, S11at40 = S40[1, 1, :])
```

A complex reference impedance defines pseudo-waves; see
[`ZtoS`](@ref JosephsonCircuits.ZtoS).

```@autodocs
Modules = [JosephsonCircuits]
Filter = obj -> Main.DocumentationAPI.in_group(obj, :networkparameters)
```
