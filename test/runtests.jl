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
    #makedocs(sitename="My Documentation",modules=[QCE])
    doctest(JosephsonCircuits, manual = false)
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

    include("parseinput.jl")

    include("qesparams.jl")

    include("spiceraw.jl")

    include("spiceutils.jl")

    include("spicewrapper.jl")

    include("testutils.jl")

    include("touchstone.jl")

end