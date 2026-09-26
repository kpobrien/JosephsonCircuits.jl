# Run docs/make.jl once to install the Node packages, then:
#     julia --project=docs docs/live.jl
# A package docstring edit needs Revise loaded first, or a Julia restart.
using DocumenterVitepress

ENV["DOCS_LIVE"] = "1"
include("make.jl")

# Poll the source tree so edits in nested recipe/asset directories and new
# directories are included. Build output is outside this tree.
function source_snapshot()
    Dict(joinpath(dir, file) => (stat(joinpath(dir, file)).mtime,
        stat(joinpath(dir, file)).size)
        for (dir, _, files) in walkdir(joinpath(@__DIR__, "src")) for file in files)
end

function live_preview()
    server = cd(@__DIR__) do
        run(`$(DocumenterVitepress.node()) node_modules/vitepress/bin/vitepress.js dev build/.documenter`;
            wait = false)
    end
    try
        previous = source_snapshot()
        while true
            sleep(0.5)
            current = source_snapshot()
            current == previous && continue
            previous = current
            try
                include("make.jl")
            catch e
                showerror(stderr, e)
                println(stderr)
            end
        end
    finally
        kill(server)
    end
end

live_preview()
