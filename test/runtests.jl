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

    include("matrices.jl")

    include("exportnetlist.jl")

    include("frequencies.jl")

    include("graph.jl")

    include("hbsolve.jl")

    include("JosephsonCircuits.jl")

    include("networkparamconversion.jl")

    include("networks.jl")

    include("networkconnection.jl")

    include("sparse.jl")

    include("newton.jl")

    include("newtonkrylov.jl")

    include("floquetdeflation.jl")

    include("modecoupling.jl")

    include("problem.jl")

    include("builders.jl")


    include("staged.jl")


    include("components.jl")

    include("parse.jl")

    include("bind.jl")

    include("canonical.jl")

    include("legacy.jl")

    include("scatteringblocks.jl")

    include("outputs.jl")

    include("layout.jl")

    include("complexjacobian.jl")

    include("realjacobian.jl")
    include("assembly.jl")
    include("devicesweep.jl")

    include("nonlinearterm.jl")
    include("system.jl")

    include("mna.jl")
    include("directcurrent.jl")

    include("spiceraw.jl")

    include("spiceutils.jl")

    include("spicewrapper.jl")

    include("quantumoptics.jl")

    include("testutils.jl")
    include("docstringchecktests.jl")

    include("unwrap.jl")

    include("deprecated.jl")
end
