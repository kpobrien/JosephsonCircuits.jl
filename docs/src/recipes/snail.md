# SNAIL Parametric Amplifier

Compare pump-on and pump-off gain for an explicitly modeled SNAIL. The quadratic nonlinearity supports three-wave mixing. Similar resonance frequencies in the two curves indicate operation near the chosen Kerr-free point.

Requires `JosephsonCircuits` and `Plots`. The optional comparison uses
WRspice through `XicTools_jll`, or a local WRspice installation.

Figures and timings come from the original reference run using 16 threads
on an AMD Ryzen 9 9950X under Linux. Rerun the code for your package version
and numerical settings; see [benchmarking](../performance.md#Measuring-performance).

The circuit parameters are from [Frattini et al. (2018)](https://doi.org/10.1103/PhysRevApplied.10.054020).

```julia
using JosephsonCircuits
using Plots

R = 50.0
Cc = 0.048e-12
Cj = 10.0e-15
Lj = 60e-12
Cr = 0.4e-12*1.25
Lr = 0.4264e-9*1.25
Ll = 34e-12
Ldc = 0.74e-12
K = 0.999 # the coupling of the bias inductor ldc to the loop inductor ll

alpha = 0.29
Z0 = 50
w0 = 2*pi*8e9
l=10e-3
circuit = Circuit(
    [(:p1, 1, 0, Port(1; Z0 = R)),
     (:cc, 1, 2, Capacitor(Cc)),
     (:lr, 2, 3, Inductor(Lr)), (:cr, 2, 0, Capacitor(Cr)),
     # the small junction of the SNAIL, across the three large ones
     (:jj1, 3, 0, JosephsonJunction(Lj/alpha)),
     (:cj1, 3, 0, Capacitor(Cj/alpha)),
     (:ll, 3, 4, Inductor(Ll)),
     (:jj2, 4, 5, JosephsonJunction(Lj)), (:cj2, 4, 5, Capacitor(Cj)),
     (:jj3, 5, 6, JosephsonJunction(Lj)), (:cj3, 5, 6, Capacitor(Cj)),
     (:jj4, 6, 0, JosephsonJunction(Lj)), (:cj4, 6, 0, Capacitor(Cj)),
     # the bias inductor, mutually coupled to the loop inductor ll
     (:ldc, 7, 0, Inductor(Ldc)),
     (:k1, :ll, :ldc, MutualInductor(K)),
     # a high impedance port, so the bias may be applied across it
     (:p2, 7, 0, Port(2; Z0 = 1000.0))])
ws = 2*pi*(7.8:0.001:8.2)*1e9
wp = (2*pi*16.00*1e9,)
Ip = 4.4e-6
Idc = 0.000159
# add the DC bias and pump to port 2
sourcespumpon = [(mode=(0,),port=2,current=Idc),(mode=(1,),port=2,current=Ip)]
sourcespumpoff = [(mode=(0,),port=2,current=Idc),(mode=(1,),port=2,current=0.0)]
Npumpharmonics = (16,)
Nmodulationharmonics = (8,)
@time jpapumpon = hbsolve(ws, wp, sourcespumpon, Nmodulationharmonics,
    Npumpharmonics, circuit, dc = true, threewavemixing=true,fourwavemixing=true) # enable dc and three wave mixing
@assert jpapumpon.nonlinear.solverinfo.converged
@time jpapumpoff = hbsolve(ws, wp, sourcespumpoff, Nmodulationharmonics,
    Npumpharmonics, circuit, dc = true, threewavemixing=true,fourwavemixing=true) # enable dc and three wave mixing
@assert jpapumpoff.nonlinear.solverinfo.converged

p1 = plot(
    jpapumpon.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpapumpon.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
    label="pump on",
)

plot!(
    jpapumpoff.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpapumpoff.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    label="pump off",
)

p2 = plot(
    jpapumpon.linearized.w/(2*pi*1e9),
    angle.(
        jpapumpon.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    ),
    xlabel="Frequency (GHz)",
    ylabel="Phase (rad)",
    label="pump on",
)

plot!(
    jpapumpoff.linearized.w/(2*pi*1e9),
    angle.(
        jpapumpoff.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    ),
    label="pump off",
)
plot(p1,p2,layout=(2,1))
```

```
  0.010345 seconds (16.74 k allocations: 40.025 MiB)
  0.011252 seconds (16.68 k allocations: 39.985 MiB)
```

![SNAIL parametric amplifier simulation with JosephsonCircuits.jl](../assets/examples/snail.png)

## Compare with WRspice

```julia
using XicTools_jll

# simulate the JPA in WRSPICE
wswrspice=2*pi*(7.8:0.005:8.2)*1e9
n = JosephsonCircuits.exportnetlist(circuit);
input = JosephsonCircuits.wrspice_input_paramp(n.netlist,wswrspice,[0.0,wp[1]],[Idc,2*Ip],[(0,1)],[(0,7),(0,7)];trise=10e-9,tstop=600e-9);

@time output = JosephsonCircuits.spice_run(input,XicTools_jll.wrspice());
S11,S21=JosephsonCircuits.wrspice_calcS_paramp(output,wswrspice,n.Nnodes);

# plot the output
plot(
    jpapumpon.linearized.w/(2*pi*1e9),
    10*log10.(abs2.(
        jpapumpon.linearized.S(
            outputmode=(0,),
            outputport=1,
            inputmode=(0,),
            inputport=1,
            freqindex=:
        ),
    )),
    xlabel="Frequency (GHz)",
    ylabel="Gain (dB)",
    label="JosephsonCircuits.jl",
)

plot!(wswrspice/(2*pi*1e9),10*log10.(abs2.(S11)),
    label="WRspice",
    seriestype=:scatter)
```

```
2067.364975 seconds (149.73 k allocations: 29.873 GiB, 0.01% gc time)
```

![SNAIL parametric amplifier simulation with JosephsonCircuits.jl and WRspice](../assets/examples/snail_WRspice.png)
