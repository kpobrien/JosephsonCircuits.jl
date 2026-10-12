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

The adjoint calculation differentiates scattering parameters with respect to design parameters without re-solving the nonlinear problem for each perturbation. The components are written as names, and a Julia function of the design parameters gives the value of each name; every component which depends on a parameter contributes through the chain rule, with the derivative of that function, which ForwardDiff takes when it is loaded beside the package (`using JosephsonCircuits, ForwardDiff`), or which `jacobian = p -> J` states.

## Check a design derivative

This small check uses a linear RC load, for which a central difference is
inexpensive. It exercises the same parameter interface as the JPA below.

```@example derivative
using JosephsonCircuits, ForwardDiff
c = Circuit([(:p1, 1, 0, Port(1)), (:c1, 1, 0, Capacitor(:C))])
C = 1e-12
ws, wp = [2pi*3e9], (2pi*5e9,)
sources = [(mode = (1,), port = 1, current = 0.0)]
values(p) = (C = p.C,)
r = designsensitivities(c, values, (C = C,), ws, wp, sources, (0,), (1,))
analytic = r.dSdp((0,), 1, (0,), 1, :C, 1)
step = 1e-5*C
response(value) = hblinsolve(ws, c, Dict(:C => value); keyedarrays = false).S[1, 1, 1]
finite_difference = (response(C + step) - response(C - step))/(2step)
@assert isapprox(analytic, finite_difference; rtol = 1e-6)
abs((analytic - finite_difference)/finite_difference)
```

## Derived values

The function of the design parameters returns the value of every name the
circuit's components are written in, as a `NamedTuple`, and a value
derived from the parameters is computed there, in Julia. Its derivatives
come from ForwardDiff, which differentiates the real and imaginary parts
of each value, so a parameter which rotates a value in the complex plane,
a loss tangent, is exact; `jacobian = p -> J` states them instead, a row
per value and a column per parameter. The function returns numbers: a
frequency dependent value, which no parameter may move, is written into
the circuit as the component's value. A component written as a number
takes no part, and every key of the parameters is a design parameter.

Here the coupling capacitance is a tenth of the resonator's, so `C` moves
both, and the derivative with respect to `C` is checked against a central
difference which moves `C`.

```@example derived
using JosephsonCircuits, ForwardDiff
c = Circuit([(:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(:Cc)),
    (:cr, 2, 0, Capacitor(:Cr)), (:lr, 2, 0, Inductor(1e-9))])
values(p) = (Cc = p.C/10, Cr = p.C)
p = (C = 1e-12,)
names, v, J = designjacobian(c, values, p)
@assert J[findfirst(==("cc"), names), 1] ≈ 0.1
ws, wp = [2pi*4.9e9], (2pi*5e9,)
sources = [(mode = (1,), port = 1, current = 0.0)]
r = designsensitivities(c, values, p, ws, wp, sources, (0,), (1,))
analytic = r.dSdp((0,), 1, (0,), 1, :C, 1)
step = 1e-6*1e-12
response(value) = hblinsolve(ws, c, Dict(pairs(values((C = value,))));
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
using JosephsonCircuits, ForwardDiff, LinearAlgebra
c = Circuit([
    (:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(:Lj)), (:cj, 2, 0, Capacitor(1e-12)),
])
L = 1e-9
wp = (2pi*4.75001e9,)
ws = 2pi .* [4.69e9, 4.71e9]
sources = [(mode = (1,), port = 1, current = 0.005e-6)]
result = designsensitivities(c, p -> (Lj = p.Lj,), (Lj = L,), ws, wp,
    sources, (4,), (8,); atol = 1e-12)
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

This independent plotting recipe requires `Plots` and ForwardDiff. The displayed
normalization multiplies each absolute design derivative by its parameter
value so the curves describe fractional changes.

```julia
using JosephsonCircuits, ForwardDiff
using Plots

jpa = Circuit(
    [(:p1, 1, 0, Port(1)),
     (:cc, 1, 2, Capacitor(:Cc)),
     (:jj, 2, 0, JosephsonJunction(:Lj)),
     (:cj, 2, 0, Capacitor(:Cj))])

# the design parameters are the values themselves
p = (Lj = 1000.0e-12, Cc = 100.0e-15, Cj = 1000.0e-15)
ws = 2*pi*(4.5:0.001:5.0)*1e9
wp = 2*pi*4.75001*1e9
sources = [(mode=(1,),port=1,current=0.00565e-6)]

r = designsensitivities(jpa, identity, p, ws, (wp,), sources, (8,), (16,))

# the derivative of the gain in dB with respect to each parameter,
# dG/dp = (20/log(10))*real(conj(S)*dS/dp)/abs2(S), scaled by the
# parameter value so the three curves share units: dB of gain per
# fractional change of the parameter
S = r.out.linearized.S((0,),1,(0,),1,:)
plot(ws/(2*pi*1e9),
    [p[q].*(20/log(10)).*
     real.(conj.(S).*r.dSdp((0,),1,(0,),1,q,:))./abs2.(S)
     for q in keys(p)],
    label=["Lj" "Cc" "Cj"],
    xlabel="Frequency (GHz)",
    ylabel="dG/dln(p) (dB)")
```
