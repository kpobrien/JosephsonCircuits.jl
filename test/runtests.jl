# The suite runs its files on worker processes, one file per job. Nearly
# all of the suite's time is the compilation each file triggers, and that
# is done in whichever process runs the file, so the files run in parallel
# and the wall time is that of the largest file plus what the files share.
# The workers load the same precompiled image the master built; nothing is
# precompiled twice. A worker is started with this process's own flags
# (`Base.julia_cmd()` carries the code coverage and allocation tracking
# options `Pkg.test` sets) plus bounds checking when this process has it,
# so a worker tests exactly what the master would have, and writes its own
# coverage files, which the coverage tools merge.
#
# As many workers run as this process has threads, each with one thread
# of its own, so that the suite uses the budget its caller set rather
# than the whole machine: `Pkg.test` gives the test process the threads
# named by `JULIA_NUM_THREADS` or its own `-t`, and a shared machine asks
# for few. `JULIA_TEST_WORKERS` overrides that count; `0` runs every job
# in this process in sequence, so that a session kept alive with Revise
# and TestEnv compiles once.

using Test
using Distributed
using Random
using JosephsonCircuits

const TESTDIR = @__DIR__

# The jobs all run from one seed, drawn afresh each run so that the
# suite keeps exploring new draws, and printed so that a run can be
# repeated: the RNG Test prints beneath a failing summary is this
# process's, which runs none of the jobs itself. `JULIA_TEST_SEED` gives
# the seed instead.
const SEED = haskey(ENV, "JULIA_TEST_SEED") ?
    parse(UInt64, ENV["JULIA_TEST_SEED"]) : rand(RandomDevice(), UInt64)

# every `.jl` under `test/` which the job list does not name and which is
# not a fixture the jobs include themselves
function uncoveredtests(listed)
    fixtures = ("runtests.jl", "testcircuits.jl", "docstringcheck.jl",
        "harmonics/layoutreference.jl")
    found = String[]
    for (root, _, names) in walkdir(TESTDIR)
        rel = relpath(root, TESTDIR)
        first(splitpath(rel)) in ("gpu", "interop") && continue
        for n in names
            endswith(n, ".jl") || continue
            push!(found, replace(rel == "." ? n : joinpath(rel, n), '\\' => '/'))
        end
    end
    return sort(setdiff(found, listed, fixtures))
end

# the jobs: a name and the code which runs it, the files heaviest first so
# the worker which draws the last job is not left with a large one
function testjobs()
    jobs = Pair{String,String}[]
    file(f) = f => "include($(repr(joinpath(TESTDIR, f))))"
    # Aqua and the doctests are not run on nightly
    if !occursin("DEV", string(VERSION))
        push!(jobs, "Doctests (Documenter.jl)" => """
            using Documenter, Logging
            DocMeta.setdocmeta!(JosephsonCircuits, :DocTestSetup,
                :(using JosephsonCircuits); recursive = true)
            # Documenter narrates each of its steps at info level, which
            # says nothing about the doctests; a failing doctest is logged
            # as an error and still comes through
            with_logger(ConsoleLogger(stderr, Logging.Warn)) do
                makedocs(remotes = nothing,
                    root = joinpath(dirname(pathof(JosephsonCircuits)), "..", "docs"),
                    modules = [JosephsonCircuits], doctest = :only,
                    sitename = "JosephsonCircuits",
                    format = Documenter.HTML(edit_link = nothing, disable_git = true))
            end
            """)
    end
    # The test tree mirrors src: the tests for `src/<folder>/<file>.jl`
    # are in `test/<folder>/<file>.jl`. A test which compares functions
    # from different files, or drives the whole package, stays at the
    # top level. The files are heaviest first so that the worker which
    # draws the last job is not left with a large one; the order is
    # measured, and worth remeasuring when a file grows.
    files = ("hbsolve.jl", "transient/solve.jl",
            "harmonics/directcurrent.jl", "transient/noise.jl",
            "networks/quantumoptics.jl", "transientpumped.jl",
            "linearized/scatteringblocks.jl", "transient/system.jl",
            "crosscheck.jl", "solvers/modecoupling.jl", "solvers/problem.jl",
            "networks/parameters.jl", "circuit/mna.jl",
            "circuit/vectorfit.jl", "solvers/floquetdeflation.jl",
            "circuit/parse.jl", "nonlinearinductor.jl",
            "linearized/devicesweep.jl", "harmonics/layout.jl",
            "linearized/designsensitivities.jl", "solvers/cache.jl",
            "networks/connections.jl", "solvers/gmres.jl",
            "harmonics/assembly.jl", "circuit/components.jl",
            "solvers/staged.jl", "transient/quantum.jl",
            "harmonics/complexjacobian.jl", "harmonics/system.jl",
            "networks/networks.jl", "solvers/newtonkrylov.jl",
            "wrspicecrosscheck.jl", "harmonics/nonlinearterm.jl",
            "circuit/values.jl", "spice/transient.jl",
            "harmonics/frequencies.jl", "transient/iq.jl", "deprecated.jl",
            "spice/export.jl", "circuit/legacy.jl", "linearized/outputs.jl",
            "circuit/bind.jl", "circuit/graph.jl", "harmonics/sparse.jl",
            "solvers/factorizations.jl", "JosephsonCircuits.jl",
            "spice/raw.jl", "spice/wrapper.jl", "solvers/newton.jl",
            "docstringchecktests.jl", "spice/utils.jl", "networks/unwrap.jl",
            "circuit/matrices.jl", "solvers/solverinfo.jl", "testutils.jl",
            "solvers/linesearch.jl")
    for f in files
        push!(jobs, file(f))
    end
    # a file in neither the list above nor `fixtures` is a file nobody
    # runs; the jobs of `gpu` and `interop` have their own runners
    push!(jobs, "the job list covers the test tree" =>
        "@test $(repr(uncoveredtests(files))) == String[]")

    if !occursin("DEV", string(VERSION))
        push!(jobs, "Code quality (Aqua.jl)" => """
            using Aqua
            Aqua.test_all(JosephsonCircuits; ambiguities = true,
                persistent_tasks = false)
            """)
    end
    return jobs
end

# what a worker is started with: this process's own flags, the test
# environment, bounds checking when this process has it, one thread, and a
# share of the memory
function workerflags(n)
    flags = filter(Base.julia_cmd().exec[2:end]) do f
        !(startswith(f, "--threads") || startswith(f, "-t") ||
          startswith(f, "--project") || startswith(f, "--check-bounds"))
    end
    push!(flags, "--project=$(Base.active_project())")
    Base.JLOptions().check_bounds == 1 && push!(flags, "--check-bounds=yes")
    push!(flags, "--startup-file=no", "--threads=1")
    push!(flags, "--heap-size-hint=$(round(Int, Sys.total_memory()/2^20/(n + 1)))M")
    return Cmd(flags)
end

jobs = testjobs()
nworkers = min(parse(Int, get(ENV, "JULIA_TEST_WORKERS",
    string(Threads.nthreads()))), length(jobs))

println("running $(length(jobs)) jobs on $(nworkers) workers from seed $(repr(SEED))")

if nworkers == 0
    include(joinpath(TESTDIR, "testcircuits.jl"))
    @testset verbose = true "JosephsonCircuits" begin
        for (name, code) in jobs
            @testset "$name" begin
                Random.seed!(SEED)
                include_string(Main, code)
            end
        end
    end
else
    addprocs(nworkers; exeflags = workerflags(nworkers))
    # on every process, the master included: the shared circuits, one BLAS
    # thread so the workers do not oversubscribe the machine
    @everywhere begin
        using Test, Random, JosephsonCircuits
        using LinearAlgebra: BLAS
        BLAS.set_num_threads(1)
        include(joinpath($TESTDIR, "testcircuits.jl"))
    end
    # a job on a worker: its testset runs inside an enclosing one, so it
    # records there rather than reporting itself as the outermost testset,
    # and comes back to the master, which nests it in the suite's.
    #
    # The enclosing testset keeps nothing and reports nothing: were it an
    # ordinary testset it would be the outermost on the worker, and would
    # print a summary of the job at its end, which the suite's own summary
    # on the master then repeats, and would throw when the job failed. A
    # failure or an error still prints where it happens, from the job's
    # testset, along with whatever the code under test logs, and travels
    # to the master in the returned testset, which is what the summary and
    # the exit status are made of. Only the documented interface of Test
    # is used, `record` and `finish` on an `AbstractTestSet`, since its
    # internals moved between Julia versions.
    @everywhere struct SilentTestSet <: Test.AbstractTestSet
        description::String
        SilentTestSet(description; kwargs...) = new(String(description))
    end
    @everywhere Test.record(::SilentTestSet, ::Any) = nothing
    @everywhere Test.finish(ts::SilentTestSet) = ts
    # The job's testset says which type it is, since a testset takes the
    # type of the one enclosing it when it is not told: it would be silent
    # too, and would swallow the job's failures rather than report them.
    # The name is bound here because `@testset` takes the type as a plain
    # name and not as `Test.DefaultTestSet` on the long term support
    # release.
    @everywhere const JobTestSet = Test.DefaultTestSet
    @everywhere function runtestjob(name::String, code::String, seed::UInt64)
        # a testset takes the task's seed as it stands, so seeding before
        # the job's outermost one seeds every testset under it
        Random.seed!(seed)
        job = Ref{Any}(nothing)
        @testset SilentTestSet "worker" begin
            job[] = @testset JobTestSet "$name" begin
                include_string(Main, code)
            end
        end
        return job[]
    end
    # each worker draws the next job until none is left
    queue = copy(jobs)
    results = Dict{String,Any}()
    @sync for w in workers()
        @async while !isempty(queue)
            name, code = popfirst!(queue)
            results[name] = try
                remotecall_fetch(runtestjob, w, name, code, SEED)
            catch e
                e
            end
        end
    end
    rmprocs(workers())
    # the suite's testset, with each job's testset nested in it in the
    # order of the job list, so the summary is the usual table and the
    # exit status the usual one; a worker which failed outright is an error
    @testset verbose = true "JosephsonCircuits" begin
        for (name, _) in jobs
            r = results[name]
            if r isa Test.AbstractTestSet
                Test.record(Test.get_testset(), r)
            else
                @testset "$name" begin
                    error("the worker running $name failed: $(sprint(showerror, r))")
                end
            end
        end
    end
end
