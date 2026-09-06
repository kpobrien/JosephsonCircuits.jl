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
# `JULIA_TEST_WORKERS` sets the worker count; `0` runs every job in this
# process in sequence, which is what a session kept alive with Revise and
# TestEnv wants, since the compilation then happens once per session.

using Test
using Distributed
using JosephsonCircuits

const TESTDIR = @__DIR__

# the jobs: a name and the code which runs it, the files heaviest first so
# the worker which draws the last job is not left with a large one
function testjobs()
    jobs = Pair{String,String}[]
    file(f) = f => "include($(repr(joinpath(TESTDIR, f))))"
    # Aqua and the doctests are not run on nightly
    if !occursin("DEV", string(VERSION))
        push!(jobs, "Doctests (Documenter.jl)" => """
            using Documenter
            DocMeta.setdocmeta!(JosephsonCircuits, :DocTestSetup,
                :(using JosephsonCircuits); recursive = true)
            makedocs(remotes = nothing,
                root = joinpath(dirname(pathof(JosephsonCircuits)), "..", "docs"),
                modules = [JosephsonCircuits], doctest = :only,
                sitename = "JosephsonCircuits",
                format = Documenter.HTML(edit_link = nothing, disable_git = true))
            """)
    end
    for f in ("hbsolve.jl", "directcurrent.jl", "quantumoptics.jl",
            "networkparamconversion.jl", "scatteringblocks.jl", "crosscheck.jl",
            "problem.jl", "modecoupling.jl", "builders.jl", "matrices.jl",
            "exportnetlist.jl", "frequencies.jl", "graph.jl",
            "JosephsonCircuits.jl", "networks.jl", "networkconnection.jl",
            "sparse.jl", "newton.jl", "newtonkrylov.jl", "floquetdeflation.jl",
            "staged.jl", "components.jl", "parse.jl", "bind.jl", "canonical.jl",
            "legacy.jl", "outputs.jl", "layout.jl", "complexjacobian.jl",
            "assembly.jl", "devicesweep.jl", "nonlinearterm.jl", "system.jl",
            "mna.jl", "spiceraw.jl", "spiceutils.jl", "spicewrapper.jl",
            "testutils.jl", "docstringchecktests.jl", "unwrap.jl",
            "deprecated.jl")
        push!(jobs, file(f))
    end
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
    string(Sys.CPU_THREADS))), length(jobs))

if nworkers == 0
    include(joinpath(TESTDIR, "testcircuits.jl"))
    @testset verbose = true "JosephsonCircuits" begin
        for (name, code) in jobs
            @testset "$name" begin
                include_string(Main, code)
            end
        end
    end
else
    addprocs(nworkers; exeflags = workerflags(nworkers))
    # on every process, the master included: the shared circuits, one BLAS
    # thread so the workers do not oversubscribe the machine
    @everywhere begin
        using Test, JosephsonCircuits
        using LinearAlgebra: BLAS
        BLAS.set_num_threads(1)
        include(joinpath($TESTDIR, "testcircuits.jl"))
    end
    # a job on a worker: its testset runs inside an enclosing one, so it
    # records there rather than reporting itself as the outermost testset,
    # and comes back to the master, which nests it in the suite's. The
    # enclosing testset is the outermost on the worker and throws at its
    # end when the job failed; that is caught, since the failures travel
    # in the returned testset. Only the documented interface of Test is
    # used, since its internals moved between Julia versions.
    @everywhere function runtestjob(name::String, code::String)
        job = Ref{Any}(nothing)
        try
            @testset "worker" begin
                job[] = @testset "$name" begin
                    include_string(Main, code)
                end
            end
        catch e
            e isa Test.TestSetException || rethrow()
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
                remotecall_fetch(runtestjob, w, name, code)
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
