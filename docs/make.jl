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
    format = Documenter.HTML(edit_link = nothing, disable_git = true,size_threshold_ignore = ["reference.md", "examples.md"]),
    pages = [
        "Home" => "index.md",
        "Circuits" => "circuits.md",
        "Harmonic balance" => ["Usage" => "harmonicbalance.md", "Theory and implementation" => "harmonicbalancetheory.md", "Examples" => "examples.md", "Using other solvers" => "interop.md"],
        "Time domain" => ["Usage" => "transient.md", "Theory and implementation" => "transienttheory.md", "Quantum noise in time" => "transientnoise.md"],
        "Reference" => "reference.md",
    ],
    )

deploydocs(
   repo = "github.com/kpobrien/JosephsonCircuits.jl.git",
)
