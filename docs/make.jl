using Test, Documenter, JosephsonCircuits

# DocMeta.setdocmeta!(JosephsonCircuits, :DocTestSetup, :(using JosephsonCircuits); recursive=true)

# makedocs(modules = [JosephsonCircuits], sitename="JosephsonCircuits.jl")

# DocMeta.setdocmeta!(JosephsonCircuits, 
#     :DocTestSetup,
#     :(using JosephsonCircuits);
#     recursive=true)

makedocs(
    root = joinpath(dirname(pathof(JosephsonCircuits)), "..", "docs"),
    modules=[JosephsonCircuits],
    doctest = false,
    sitename="JosephsonCircuits",
    format = Documenter.HTML(edit_link = nothing, disable_git = true,size_threshold_ignore = ["reference.md"]),
    )

deploydocs(
   repo = "github.com/kpobrien/JosephsonCircuits.jl.git",
)
