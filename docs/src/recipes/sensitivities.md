# Design parameter sensitivities

Differentiate the JPA response with respect to its design parameters. The small executable example first checks a derivative against finite differences; the plotting recipe then expresses gain derivatives in dB per fractional parameter change.

The plotting code requires `Plots` in addition to `JosephsonCircuits`.

The first check differentiates a capacitor on a port; the JPA of the
later sections has its junction inductance and capacitances as design
parameters:

```text
 1                  2
 o------[cc]--------o--------+
 |                  |        |
[p1]              [jj]     [cj]
 |                  |        |
 o------------------o--------+
 0
```

The adjoint calculation differentiates scattering parameters with respect to design parameters without re-solving the nonlinear problem for each perturbation. Every component value which depends on a parameter contributes through the chain rule with its exact derivative, including derived values such as `Capacitor(Cj/4)`.

## Check a design derivative

This small check uses a linear RC load, for which a central difference is
inexpensive. It exercises the same parameter interface as the JPA below.

```@example derivative
using JosephsonCircuits
c = Circuit([(:p1, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(:C))])
C = 1e-12
ws, wp = [2pi*3e9], (2pi*5e9,)
sources = [(mode = (1,), port = 1, current = 0.0)]
r = designsensitivities(c, Dict(:C => C), ws, wp, sources, (0,), (1,))
analytic = r.dSdp((0,), 1, (0,), 1, :C, 1)
step = 1e-5*C
response(value) = hblinsolve(ws, c, Dict(:C => value); keyedarrays = false).S[1, 1, 1]
finite_difference = (response(C + step) - response(C - step))/(2step)
@assert isapprox(analytic, finite_difference; rtol = 1e-6)
abs((analytic - finite_difference)/finite_difference)
```

## Derived definitions

A parameter may be defined in terms of others, and the derivative follows
the definition by the chain rule. The forms a definition takes:

| definition | example | derivative |
|---|---|---|
| a number, real or complex | `C => 1e-12` | one with respect to itself; a default parameter |
| a value written in other parameters, to any depth | `:Cc => C/10` | the chain rule through the definition: `Capacitor(:Cc)` has the derivative `1/10` with respect to `C` |
| a frequency dependent value | `:Cd => FrequencyDependent(w -> ...)` | none: no parameter moves it |

The keys may be symbols, strings, the parameters of `@params` or, with
Symbolics loaded, `Num`s. A component value is a number, a name, an
expression in names, or a frequency dependent value, and its derivative
is exact: the expression is differentiated symbolically and evaluated at
the definitions. A derived parameter is not among the defaults; selected
by name, it is differentiated as replacing its definition by its value
would move it, the parameters it is written in held. A component whose
value is frequency dependent cannot be moved by a selected parameter; one
which no selected parameter moves takes no part.

Here the coupling capacitance is defined as a tenth of the resonator's, so
`C` moves both, and the derivative with respect to `C` is checked against
a central difference which moves `C` alone in the definitions.

```@example derived
using JosephsonCircuits
JosephsonCircuits.@params C
c = Circuit([(:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(:Cc)),
    (:cr, 2, 0, Capacitor(C)), (:lr, 2, 0, Inductor(1e-9))])
definitions = Dict(:Cc => C/10, C => 1e-12)
names, values, J = designjacobian(c, definitions)
@assert J[findfirst(==("cc"), names), 1] ≈ 0.1
ws, wp = [2pi*4.9e9], (2pi*5e9,)
sources = [(mode = (1,), port = 1, current = 0.0)]
r = designsensitivities(c, definitions, ws, wp, sources, (0,), (1,))
analytic = r.dSdp((0,), 1, (0,), 1, :C, 1)
step = 1e-6*1e-12
response(value) = hblinsolve(ws, c, merge(definitions, Dict(C => value));
    keyedarrays = false).S[1, 1, 1]
finite_difference = (response(1e-12 + step) - response(1e-12 - step))/(2step)
@assert isapprox(analytic, finite_difference; rtol = 1e-6)
(names .=> real.(J[:, 1]), abs((analytic - finite_difference)/finite_difference))
```

## Total and frozen-pump derivatives

In a pumped circuit, a component change also moves the nonlinear operating
point. `designsensitivities` includes that motion. A finite-difference
check must therefore solve the pump again at both perturbed values.
Here a single design parameter is also a single junction's inductance, so
we can compare with the component sensitivity at fixed pump state.

```@example pumpderivative
using JosephsonCircuits, LinearAlgebra
c = Circuit([
    (:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(:Lj)), (:cj, 2, 0, Capacitor(1e-12)),
])
L = 1e-9
wp = (2pi*4.75001e9,)
ws = 2pi .* [4.69e9, 4.71e9]
sources = [(mode = (1,), port = 1, current = 0.005e-6)]
result = designsensitivities(c, Dict(:Lj => L), ws, wp, sources,
    (4,), (8,); atol = 1e-12)
total = collect(result.dSdp((0,), 1, (0,), 1, :Lj, :))
response(value) = hbsolve(ws, wp, sources, (4,), (8,), c,
    Dict(:Lj => value); atol = 1e-12)
step = 1e-5*L
plus, minus = response(L + step), response(L - step)
@assert all(r -> r.nonlinear.solverinfo.converged, (result.out, plus, minus))
finite_difference = (plus.linearized.S((0,), 1, (0,), 1, :) .-
    minus.linearized.S((0,), 1, (0,), 1, :))/(2step)
total_error = norm(total - finite_difference)/norm(finite_difference)
@assert total_error < 1e-4
total_error
```

To isolate the direct change at **fixed periodic flux**, request
`sensitivityoperatingpoint = false` from `hbsolve`. This switch belongs to
`hbsolve`, not `hblinsolve`. Component sensitivities are relative to the
component value, so divide the junction's result by `L` to compare with
the absolute design derivative `dS/dLj` (units H⁻¹).

```@example pumpderivative
fixed = hbsolve(ws, wp, sources, (4,), (8,), c, Dict(:Lj => L);
    atol = 1e-12, sensitivitynames = ["jj"], returnSsensitivity = true,
    sensitivityoperatingpoint = false)
frozen = collect(fixed.linearized.Ssensitivity((0,), 1, (0,), 1, "jj", :))/L
# Reuse the exact same nonlinear state on both sides of this difference.
fixedresponse(value) = hblinsolve(ws, c, Dict(:Lj => value);
    nonlinear = result.out.nonlinear, Nmodulationharmonics = (4,))
fixed_difference = (fixedresponse(L + step).S((0,), 1, (0,), 1, :) .-
    fixedresponse(L - step).S((0,), 1, (0,), 1, :))/(2step)
frozen_error = norm(frozen - fixed_difference)/norm(fixed_difference)
@assert frozen_error < 1e-4
relative_difference = norm(total - frozen)/norm(total)
@assert relative_difference > 0.1
(total_error, frozen_error, relative_difference)
```

The total and frozen derivatives differ by about 57% in norm at this
operating point; they answer different questions. A component tolerance
or design change with a fixed external pump drive generally needs the
total derivative. The frozen derivative holds the already solved junction
flux waveform fixed even though the perturbed circuit would not produce
that waveform under the same drive.

When checking another design, vary the finite-difference step and tighten
the nonlinear tolerance until the comparison stabilizes. Both perturbed
solves must follow the same operating-point branch. Close to a bifurcation,
branch switching or an ill-conditioned pump Jacobian can invalidate a
naive finite-difference comparison.

## Gain derivatives of a pumped JPA

This independent plotting recipe requires `Plots`. The displayed
normalization multiplies each absolute design derivative by its parameter
value so the curves describe fractional changes.

```julia
using JosephsonCircuits
using Plots

Lj, Cc, Cj = JosephsonCircuits.@params Lj Cc Cj
jpa = Circuit(
    [(:p1, 1, 0, Port(1)),
     (:cc, 1, 2, Capacitor(Cc)),
     (:jj, 2, 0, JosephsonJunction(Lj)),
     (:cj, 2, 0, Capacitor(Cj))])

p = Dict(Lj => 1000.0e-12, Cc => 100.0e-15, Cj => 1000.0e-15)
ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = 2*pi*4.75001*1e9
sources = [(mode=(1,),port=1,current=0.00565e-6)]

r = designsensitivities(jpa, p, ws, (wp,), sources, (8,), (16,))

# the derivative of the gain in dB with respect to each parameter,
# dG/dp = (20/log(10))*real(conj(S)*dS/dp)/abs2(S), scaled by the
# parameter value so the three curves share units: dB of gain per
# fractional change of the parameter
S = r.out.linearized.S((0,),1,(0,),1,:)
plot(ws/(2*pi*1e9),
    [p[q].*(20/log(10)).*
     real.(conj.(S).*r.dSdp((0,),1,(0,),1,q.name,:))./abs2.(S)
     for q in (Lj, Cc, Cj)],
    label=["Lj" "Cc" "Cj"],
    xlabel="Frequency (GHz)",
    ylabel="dG/dln(p) (dB)")
```
