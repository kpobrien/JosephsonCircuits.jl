# Performance and reuse

Measure the analysis you intend to run: a cold solve, a sweep of values,
a signal-frequency sweep, and a batch of transient drives have different
costs. Check numerical convergence before comparing their timings.

## Reuse across a sweep of values

A parameterized circuit lets a cache retain topology, transform plans,
and symbolic factorization while component values change.

```@example cache
using JosephsonCircuits, LinearAlgebra
circuit = Circuit([
    (:P1, 1, 0, Port(1)), (:C1, 1, 2, Capacitor(:Cc)),
    (:Lj1, 2, 0, JosephsonJunction(:Lj)), (:C2, 2, 0, Capacitor(1e-12)),
])
cache = hbcache((2pi*4.75001e9,), (8,),
    [(mode = (1,), port = 1, current = 0.002e-6)],
    circuit, Dict(:Lj => 1e-9, :Cc => 100e-15))
inductances = (0.95e-9, 1.0e-9, 1.05e-9)
function solve_sweep(cache, inductances)
    map(inductances) do Lj
        result = hbsolve!(cache, (Lj = Lj,))
        @assert result.solverinfo.converged
        result
    end
end
operatingpoints = solve_sweep(cache, inductances)
nothing # hide
```

`hbsolve!` returns a nonlinear operating point. For a signal sweep, pass
that result to `hblinsolve` with the **same component values**:

```@example cache
defs = Dict(:Lj => last(inductances), :Cc => 100e-15)
response = hblinsolve(2pi .* [4.6e9, 4.7e9, 4.8e9], circuit, defs;
    nonlinear = last(operatingpoints), Nmodulationharmonics = (8,))
@assert all(isfinite, response.S)
nothing # hide
```

A cache is mutable workspace: use one per simultaneous solve and rebuild
after changing topology or the retained grid. It reuses the nonlinear
setup and starts from the last converged operating point. If a warm solve
fails, `hbsolve!` retries once from a cold start. Its result contains the
warm attempt's stages followed by the cold attempt's stages. A successful
retry can land on a different branch; convergence alone does not establish
continuity of a sweep. If both attempts fail, the cache retains its last
successful starting point for the next call.

Compare cold and warm answers at a suspect parameter value, and repeat a
sweep in reverse. This weakly driven example reaches the same state:

```@example cache
warm = hbsolve!(cache, (Lj = last(inductances),))
cold = hbsolve!(cache, (Lj = last(inductances),); warmstart = false)
@assert warm.solverinfo.converged && cold.solverinfo.converged
@assert isapprox(warm.nodeflux, cold.nodeflux; rtol = 1e-6, atol = 1e-9)
[(stage.converged, stage.iterations) for stage in cold.solverinfo.stages]
```

`warmstart=false` makes one cold attempt. `JosephsonCircuits.reset!(cache)`
also discards the stored operating point and learned preconditioner state.
Neither operation guarantees selection of a particular physical branch.
For a jump, record the state, response, and stage history, reduce the
parameter spacing, and check [stability](stability.md).

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

### A reproducible CPU measurement

Run the following in a fresh Julia process to include first-use
compilation in `first_call`. In the documentation build, earlier examples
may already have compiled these paths. Use the same Julia version,
thread counts, circuit, tolerances, output flags, and harmonic limits for
comparisons. These are measurements on your machine, not published speed
claims.

```@example timing
using JosephsonCircuits, LinearAlgebra
BLAS.set_num_threads(1)
benchmark_circuit = Circuit([
    (:p1, 1, 0, Port(1)), (:cc, 1, 2, Capacitor(100e-15)),
    (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(1e-12)),
])
function benchmark_solve(; backend = JosephsonCircuits.CPU())
    result = hbsolve(2pi .* [4.6e9, 4.7e9, 4.8e9], (2pi*4.75001e9,),
        [(mode = (1,), port = 1, current = 0.002e-6)],
        (4,), (8,), benchmark_circuit; backend)
    @assert result.nonlinear.solverinfo.converged
    return result
end
first_call = @timed benchmark_solve()
warm_calls = [@timed(benchmark_solve()) for _ in 1:3]
(; julia = VERSION, julia_threads = Threads.nthreads(),
    blas_threads = BLAS.get_num_threads(), first_seconds = first_call.time,
    warm_seconds = [r.time for r in warm_calls],
    warm_allocated_bytes = [r.bytes for r in warm_calls])
```

“Warm” here means compiled code; each call still constructs a fresh solve.
To measure cache reuse, time the complete same-parameter-path sweep with
`hbcache`/`hbsolve!`, including the same cold or warm initialization policy
in every trial. Repeating an already-converged parameter value can do
almost no nonlinear work and is not representative of a design sweep.
Julia threads parallelize independent work; BLAS threads parallelize
supported dense kernels. Start with one BLAS thread when using many Julia
threads, then measure alternatives in separate processes to avoid
oversubscribing the CPU. `@timed.bytes` is allocated host memory over the
call, not peak resident memory or device memory.

The following continues the cache example at the top of this page. Reset
before each timed trial so that each sweep begins cold and then warm-starts
between the same three values. Cache construction is excluded; state and
preconditioner setup within the first solve are included.

```@example cache
cache_trials = map(1:3) do _
    JosephsonCircuits.reset!(cache)
    @timed solve_sweep(cache, inductances)
end
[(seconds = trial.time, allocated_bytes = trial.bytes) for trial in cache_trials]
```

### Compare with a GPU

In the same setup, with compatible CUDA/CUDSS installations and a supported
GPU, load the extensions before warming the device path:

```julia
using CUDA, CUDSS
CUDA.allowscalar(false)
backend = CUDABackend()
gpu_warmup = benchmark_solve(; backend)
CUDA.synchronize()
gpu_seconds = [@elapsed(begin
    CUDA.synchronize()
    result = benchmark_solve(; backend)
    CUDA.synchronize()
end) for _ in 1:3]

# Report download cost separately if downstream work needs host arrays.
download_seconds = @elapsed begin
    host_S = Array(gpu_warmup.linearized.S)
    CUDA.synchronize()
end
@assert isapprox(host_S, Array(first_call.value.linearized.S); rtol = 1e-6, atol = 1e-8)
(; gpu_seconds, download_seconds)
```

The timed device calls include setup and any transfers inside `hbsolve`;
they exclude the final explicit download. Report an end-to-end time too
when the application downloads every result. Synchronization is necessary
because launches may otherwise return before work completes. This tiny
circuit checks the workflow; it is unlikely to demonstrate GPU throughput.
For a useful comparison use a representative circuit/batch that fits both
backends, and check numerical agreement at the same accuracy.

### Estimate storage before scaling up

For retained caps `H` in `d` tone dimensions, the full rectangular lattice
has `prod(2 .* H .+ 1)` points. Parity, conjugate redundancy, and order or
frequency cuts reduce the unknown count; use `length(nonlinear.modes)` for
the actual count. Evaluation caps `E` require a real transform grid with
`prod(2 .* E .+ 1)` samples per evaluated nonlinear waveform. For example,
`H=(4,4,4)` with default `E=(8,8,8)` has 729 lattice points and 4913
evaluation samples. One real double-precision evaluation buffer for 1000
junctions alone is about 39 MB; transforms need several buffers, and
factorizations and Krylov vectors are additional.

For a transient with `nt` saved times, one real trace of `n` variables takes
`8*n*nt` bytes. `record=:ports` retains several port traces;
`record=:states` additionally stores internal histories. At 1000 variables
and 100,000 times, one history is 800 MB. Account separately for final
states, line histories, rational states, tangent/adjoint work, and the
[noise accumulator](#Noise-calculation-cost). `Base.summarysize(result)`
can check a host result's retained storage; it does not capture temporary
peak usage or reliably account for device allocations.
