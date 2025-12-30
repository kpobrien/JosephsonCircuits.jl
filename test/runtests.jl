using Aqua
using Documenter
using Test
using JosephsonCircuits


@testset verbose = true "JosephsonCircuits" begin

    # don't run Aqua and Doctests on nightly
    if !occursin("DEV", string(VERSION))
        @testset verbose=true "Code quality (Aqua.jl)" begin
            using Aqua
            Aqua.test_all(JosephsonCircuits; ambiguities = false, persistent_tasks=false)
        end

        @testset verbose = true "Doctests (Documenter.jl)" begin
            using Documenter

            DocMeta.setdocmeta!(JosephsonCircuits, 
                :DocTestSetup,
                :(using JosephsonCircuits);
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

    include("hbsolveold.jl")

    include("JosephsonCircuits.jl")

    include("networkparamconversion.jl")

    include("networks.jl")

    include("networkconnection.jl")

    include("matutils.jl")

    include("nlsolve.jl")

    include("parseinput.jl")

    include("qesparams.jl")

    include("spiceraw.jl")

    include("spiceutils.jl")

    include("spicewrapper.jl")

    include("testutils.jl")

    include("deprecated.jl")
end
