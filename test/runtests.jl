using Aqua
using Documenter
using Test
using JosephsonCircuits

include("testcircuits.jl")


@testset verbose = true "JosephsonCircuits" begin

    # don't run Aqua and Doctests on nightly
    if !occursin("DEV", string(VERSION))
        @testset verbose=true "Code quality (Aqua.jl)" begin
            using Aqua
            Aqua.test_all(JosephsonCircuits; ambiguities = true, persistent_tasks=false)
        end

        @testset verbose = true "Doctests (Documenter.jl)" begin
            using Documenter
            # Load the symbolic packages into Main as well as the doctest
            # sandbox: array element types print module-qualified when the
            # type is not visible from Main, so without this `Num[...]`
            # doctests render as `Symbolics.Num[...]` and fail.
            using Symbolics, SymbolicUtils

            DocMeta.setdocmeta!(JosephsonCircuits, 
                :DocTestSetup,
                :(using JosephsonCircuits; using Symbolics; using SymbolicUtils);
                recursive=true)
            makedocs(
                remotes = nothing,
                root = joinpath(dirname(pathof(JosephsonCircuits)), "..", "docs"),
                modules=[JosephsonCircuits],
                doctest = :only,
                sitename="JosephsonCircuits",
                format = Documenter.HTML(edit_link = nothing, disable_git = true),
                )
        end
    end

    include("capindmat.jl")

    include("exportnetlist.jl")

    include("fftutils.jl")

    include("graphproc.jl")

    include("hbsolve.jl")

    include("JosephsonCircuits.jl")

    include("networkparamconversion.jl")

    include("networks.jl")

    include("networkconnection.jl")

    include("matutils.jl")

    include("nlsolve.jl")

    include("krylov.jl")

    include("batchedblocks.jl")

    include("modepreconditioner.jl")

    include("hbnonlinearproblem.jl")

    include("builders.jl")


    include("stagedsolve.jl")

    include("circuitvalue.jl")

    include("circuitmodel.jl")

    include("parseinput.jl")

    include("circuitbind.jl")

    include("legacyadapter.jl")

    include("scatteringstamp.jl")

    include("qesparams.jl")

    include("realcomplexconv.jl")

    include("complexjacobian.jl")

    include("realjacobian.jl")
    include("structureassembly.jl")
    include("devicelinsolve.jl")

    include("nonlinearterm.jl")
    include("hbsystem.jl")

    include("mna.jl")
    include("dcconductance.jl")

    include("spiceraw.jl")

    include("spiceutils.jl")

    include("spicewrapper.jl")

    include("quantumoptics.jl")

    include("testutils.jl")

    include("unwrap.jl")

    include("deprecated.jl")
end
