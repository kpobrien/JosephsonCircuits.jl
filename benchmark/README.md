# The benchmarks

Each script writes one tab separated row per case to standard output, and
takes the label of the revision as its argument:

```sh
julia --startup-file=no --project=. benchmark/frontend.jl main
julia --startup-file=no --project=. benchmark/circuit-preparation.jl main
```

- `frontend.jl`: constructing, parsing, elaborating and compiling circuits
  written in each input form, and a hierarchy of subcircuits.
- `frontend-followup.jl`: the same measurements for the interface-heavy
  inputs, and it includes `frontend.jl`, so it reports both.
- `circuit-preparation.jl`: compiling a long ladder, planning its matrices,
  assembling and refilling them, and evaluating a rational scattering
  block.
- `sensitivities.jl`: the residual derivatives of a pumped junction chain
  and the sensitivity sweeps from its operating point, for two components,
  where the setup of each call is most of the cost, and for every one.

`--smoke` runs every case at a small size with few samples. Nothing it
prints is a measurement; it is what `benchmark/smoke.jl` and the continuous
integration job use to check that the scripts still run.

Compare two revisions by running the same script in a fresh process for
each, with one Julia thread and one BLAS thread, from the same depot, and
take several samples: the single-shot figures are noisy, the allocation
counts are not. A revision's precompilation, load time and first solve are
measured from a clean depot, in a process which does nothing else.
