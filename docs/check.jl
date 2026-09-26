# Execute the tutorial examples and check Documenter references without Node.
# Package docstring doctests run separately in the package test suite.
using Documenter, JosephsonCircuits
include("pages.jl")

makedocs(
    root = @__DIR__,
    build = "build-check",
    modules = [JosephsonCircuits],
    authors = "Kevin P. O'Brien and contributors",
    sitename = "JosephsonCircuits.jl",
    doctest = false,
    remotes = nothing,
    # Keep the legacy reference anchors on one page; other pages retain the
    # normal HTML size limit. VitePress is the production renderer.
    format = Documenter.HTML(disable_git = true, edit_link = nothing,
        repolink = "https://github.com/kpobrien/JosephsonCircuits.jl",
        size_threshold_ignore = ["reference.md"]),
    pages = docpages,
)
