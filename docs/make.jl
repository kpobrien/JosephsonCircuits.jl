using Documenter, DocumenterVitepress, JosephsonCircuits
# Load the plotting backend before Documenter evaluates isolated examples.
using Plots

include("readme.jl")
include("pages.jl")
include("api.jl")

# The doctests run in the test suite, not here.
makedocs(
    root = @__DIR__,
    modules = [JosephsonCircuits],
    authors = "Kevin P. O'Brien and contributors",
    sitename = "JosephsonCircuits.jl",
    doctest = false,
    format = DocumenterVitepress.MarkdownVitepress(
        repo = "github.com/kpobrien/JosephsonCircuits.jl",
        devbranch = "main",
        devurl = "dev",
        # the site is served from its own domain, so every version is built
        # for the root of that domain rather than for a repository subpath
        deploy_url = "https://josephsoncircuits.org",
        description = "Frequency and time domain simulation of superconducting circuits with Josephson junctions",
        # with DOCS_LIVE set only the markdown is rendered, for a running
        # VitePress development server to pick up
        build_vitepress = !haskey(ENV, "DOCS_LIVE"),
    ),
    pages = docpages,
)

# Deployment is opt-in; local builds and live previews only render the site.
if get(ENV, "DOCS_DEPLOY", "false") == "true" && !haskey(ENV, "DOCS_LIVE")
    DocumenterVitepress.deploydocs(
        repo = "github.com/kpobrien/JosephsonCircuits.jl",
        target = joinpath(@__DIR__, "build"),
        branch = "gh-pages",
        devbranch = "main",
        cname = "josephsoncircuits.org",
    )
end
