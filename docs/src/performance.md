# Performance and reuse

Measure the analysis you intend to run: a cold solve, a sweep of values,
a signal-frequency sweep, and a batch of transient drives have different
costs. Check numerical convergence before comparing their timings.

## Reuse across a sweep of values

A parameterized circuit lets a cache retain topology, transform plans,
and symbolic factorization while component values change.

```@example cache
using JosephsonCircuits
circuit = Circuit([
    (:P1, 1, 0, Port(1)), (:C1, 1, 2, Capacitor(:Cc)),
    (:Lj1, 2, 0, JosephsonJunction(:Lj)), (:C2, 2, 0, Capacitor(1e-12)),
])
cache = hbcache((2pi*4.75001e9,), (8,),
    [(mode = (1,), port = 1, current = 0.002e-6)],
    circuit, Dict(:Lj => 1e-9, :Cc => 100e-15))
for Lj in (0.95e-9, 1.0e-9, 1.05e-9)
    result = hbsolve!(cache, (Lj = Lj,))
    @assert cache.converged
end
nothing # hide
```

The example reuses the nonlinear operating-point setup. See [`hbcache`](@ref)
for signal-sweep options. A cache is mutable workspace: do not share it
between simultaneous solves. Rebuild after changing topology. For
repeated signal sweeps at one operating point, call `hblinsolve` with the
existing nonlinear result.

## Harmonic-balance solvers

The default `Automatic()` preconditioner uses the full sparse Jacobian for
one pump. With multiple pumps it prefers full single-precision node-block
factors when their estimated memory fits within half the available memory;
otherwise it uses a measured harmonic band. Stalled iterations can expand
the coupling or increase precision within the memory budget.

For a large multi-tone grid, reducing retained coupling lowers factor
storage but can increase Krylov work. Compare total solve time, convergence,
and memory. A cheap preconditioner that needs many iterations can be slower
than a larger factorization.

`NewtonKrylov(refresh=Probe())` can reuse a preconditioner when measured
iteration costs justify it. The default rebuilds before every Newton step.
`NewtonKrylov(precision=Float32)` changes iteration precision; it is
separate from the precision chosen for an approximate factorization.

The linearized sweep can parallelize independent frequencies. Start Julia
with a deliberate thread budget, for example `julia --threads=4`, and
check `Threads.nthreads()`. The best count depends on circuit size, memory
bandwidth, and other concurrent work. See the solver's `nbatches` option
when tuning the frequency batches.

## Transient records and reuse

Use `record=:ports` when only port waveforms and the final state are needed.
Use phases or states for response calculations, or checkpoints under
`GaussLegendre()` for long records. Checkpoint replay trades additional
integration for lower stored state history. Port traces, objective weights,
and some drive-direction arrays can still grow with time.

A [`TransientReuse`](@ref) passed through the `reuse` keyword retains the
scaled system, factorization workspace, and Krylov buffers between
compatible calls. The problem, step, rule, and backend must remain
compatible; see its reference entry for the contract.

For many drive conditions on one topology, use a vector of rebound
problems. Batched Gauss–Legendre steps evaluate the conditions together
while retaining separate factorizations and convergence checks. This
amortizes kernel launches on a GPU and can matter more than accelerating
one small trajectory.

Under `Trapezoidal()`, `linearsolver=GMRES()` uses matrix-free corrections
with a reused factorization as preconditioner. Direct factorization is
often efficient for a one-dimensional sparse line. Iteration is more
attractive when factor fill makes refactorization expensive.

## Noise calculation cost

In harmonic balance, the quantum efficiency, `nbar`, `CM` and `Cnoise` need
the noise of every channel at each frequency: a transposed solve per
frequency and a reduction over the channels. `nbar` and warm port
terminations reuse the output noise the quantum efficiency computes. `Vout`
also needs `Cnoise`, a product over every channel, and adds a product the
size of the scattering matrix. A sweep which needs only `S` skips the noise
with `returnQE = false, returnCM = false, returnnbar = false`.

In a transient, the default adjoint method scales with the measured
quadratures rather than propagating every bath direction through the
trajectory. It still
needs a stationary solve at each bath frequency and contractions over the
frequency grid.

For `nb` baths, `m` measured quadratures, `nf` frequencies, and `N`
conditions, its complex response accumulator takes approximately
`16*nb*m*nf*N` bytes. Ordinary independent baths can tile conditions and
frequencies to fit a fraction of available backend memory, at the cost of
additional adjoint passes.

Pumped-block noise couples frequency partners. Those correlated
frequencies must be contracted together, so frequency tiling cannot remove
that memory requirement. Choose the bath grid by convergence of the
measured modes, including relevant conversion bands and resonances;
reducing it merely to fit memory can change the answer.

## GPU execution

Load CUDA and CUDSS before selecting the CUDA backend. This fragment
continues a circuit/problem setup from the HB or transient guide:

```julia
using JosephsonCircuits, CUDA, CUDSS
CUDA.allowscalar(false)
sol = hbsolve(ws, wp, sources, (8,), (16,), circuit; backend = CUDABackend())
solution = transientsolve(problem, (0.0, 8e-9); dt = 2e-12,
    backend = CUDABackend())
outgoing = Array(solution.outgoing)  # download for host analysis
```

HB transforms, products, and supported factorizations run on the device.
Some work, such as frequency-dependent value evaluation and certain
scattering-parameter sensitivity stamps, uses host paths. A transient
keeps state, assembly, block operators, and line histories on the device,
while source callables and the Newton control remain on the host.

Each transient step synchronizes for convergence decisions, and algebraic
projection can require small host solves. A single small circuit may be
slower on a GPU. Batched conditions provide more parallel work per launch.
Measure transfers and setup as part of the workload you actually use.

## Measuring performance

Report Julia and package versions, hardware, threads, circuit size,
harmonic limits or time step, requested outputs, and whether setup and
compilation are included. Distinguish first-call time from repeated calls
with the same argument types and from cache reuse.

The displayed timings retained in older amplifier examples came from
16 threads on an AMD Ryzen 9 9950X under Linux. They are historical context,
not a benchmark of the current solver revision. The scripts under
`benchmark/` provide reproducible workloads; include their parameters when
reporting a new measurement.
