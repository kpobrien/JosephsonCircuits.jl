using Aqua
using Documenter
using Test
using JosephsonCircuits



@testset verbose = true "Code quality (Aqua.jl)" begin
    Aqua.test_all(JosephsonCircuits; ambiguities = false)
end

@testset verbose = true "Doctests (Documenter.jl)" begin
    DocMeta.setdocmeta!(JosephsonCircuits, 
        :DocTestSetup,
        :(using JosephsonCircuits);
        recursive=true)
    makedocs(
        root = joinpath(dirname(pathof(JosephsonCircuits)), "..", "docs"),
        modules=[JosephsonCircuits],
        doctest = :only,
        sitename="JosephsonCircuits",
        format = Documenter.HTML(edit_link = nothing, disable_git = true),
        )
    # doctest(JosephsonCircuits, manual = false)
end

@testset verbose = true "JosephsonCircuits" begin
    # @info ""
    include("capindmat.jl")

    include("exportnetlist.jl")

    include("fftutils.jl")

    include("graphproc.jl")

    include("hbsolve.jl")

    include("hbsolve2.jl")

    include("JosephsonCircuits.jl")

    include("matutils.jl")

    include("nlsolve.jl")

    include("parseinput.jl")

    include("qesparams.jl")

    include("spiceraw.jl")

    include("spiceutils.jl")

    include("spicewrapper.jl")

    include("testutils.jl")

    include("touchstone.jl")

end