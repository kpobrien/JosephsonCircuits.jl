# Design parameter sensitivities

Differentiate the JPA response with respect to its design parameters. The small executable example first checks a derivative against finite differences; the plotting recipe then expresses gain derivatives in dB per fractional parameter change.

The plotting code requires `Plots` in addition to `JosephsonCircuits`.

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

@time r = designsensitivities(jpa, p, ws, (wp,), sources, (8,), (16,))

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
