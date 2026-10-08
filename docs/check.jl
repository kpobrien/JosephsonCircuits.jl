# Execute the tutorial examples and check Documenter references without Node.
# Package docstring doctests run separately in the package test suite.
using Documenter, JosephsonCircuits
# Load the plotting backend before Documenter evaluates isolated examples.
using Plots
include("pages.jl")
include("api.jl")

makedocs(
    root = @__DIR__,
    build = "build-check",
    modules = [JosephsonCircuits],
    authors = "Kevin P. O'Brien and contributors",
    sitename = "JosephsonCircuits.jl",
    doctest = false,
    remotes = nothing,
    # The developer appendix is longer. Public API pages retain the normal
    # size limit.
    format = Documenter.HTML(disable_git = true, edit_link = nothing,
        repolink = "https://github.com/kpobrien/JosephsonCircuits.jl",
        size_threshold_ignore = ["api/internals.md"]),
    pages = docpages,
)
