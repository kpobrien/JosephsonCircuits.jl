using Documenter, DocumenterVitepress, JosephsonCircuits

# The doctests run in the test suite, not here.
makedocs(
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
    pages = [
        "Home" => "index.md",
        "Circuits" => "circuits.md",
        "Harmonic balance" => [
            "Usage" => "harmonicbalance.md",
            "Theory and implementation" => "harmonicbalancetheory.md",
            "Examples" => "examples.md",
            "Using other solvers" => "interop.md",
        ],
        "Time domain" => [
            "Usage" => "transient.md",
            "Theory and implementation" => "transienttheory.md",
            "Quantum noise in time" => "transientnoise.md",
        ],
        "Reference" => "reference.md",
    ],
)

DocumenterVitepress.deploydocs(
    repo = "github.com/kpobrien/JosephsonCircuits.jl",
    target = joinpath(@__DIR__, "build"),
    branch = "gh-pages",
    devbranch = "main",
    cname = "josephsoncircuits.org",
)
