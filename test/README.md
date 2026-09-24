# The test suites

`Pkg.test()` runs `test/runtests.jl`, which distributes the files of this
directory over worker processes. Three suites are outside it, because each
needs an environment or a process layout of its own: the extension suite
and the serial/threaded parity suite each have a job in
`.github/workflows/CI.yml`, and the GPU suite is run by hand, since it
needs a CUDA device that the hosted runners do not have. A green CI run
therefore says nothing about the GPU backend.

## The package suite

```sh
julia --project=. -e 'using Pkg; Pkg.test()'
```

`JULIA_TEST_SEED` fixes the seed of the random circuits, `JULIA_TEST_WORKERS`
the number of worker processes. A file which is in neither the job list nor
the fixture list fails the "the job list covers the test tree" job, so a new
test file has to be named in `testjobs()`.

## The extensions

```sh
julia --project=test/interop -e 'using Pkg; Pkg.develop(path = "."); Pkg.instantiate()'
julia --project=test/interop test/interop/runtests.jl
```

Krylov.jl, SciMLBase and Symbolics live in `test/interop`'s own environment,
so the package suite pays no resolve or precompile cost for them. The suite
loads the package before its triggers and does not repeat itself in the
other order: these extensions define methods and types and nothing else, so
the order decides when their methods reach the dispatch tables and not what
they do there.

## One thread against four

```sh
julia --project=. test/threaded/compare.jl
```

Starts a one-thread and a four-thread process, runs the batched linearized
sweep, concurrent solve caches over one compiled circuit, and the batched
transient with its tangent, adjoint, noise and gain in each, and checks that
the two processes wrote the same numbers.

## The GPU

```sh
julia --project=test/gpu test/gpu/runtests.jl
```

Needs a CUDA device. It compares every device result against the host
result of the same solve.

## The benchmarks

`benchmark/smoke.jl` runs every benchmark script at a small size and
measures nothing; it is the check that the scripts still run against the
package's current interfaces. `benchmark/README.md` describes the
measurements themselves.
