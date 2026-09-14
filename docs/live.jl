# A live preview of the documentation. The markdown is rendered by
# make.jl, the VitePress development server serves it and reloads the
# browser on every change, and the render repeats whenever a page under
# docs/src changes. A docstring edit needs Revise loaded first, or a
# restart, since the package itself is not reloaded.
#
# Run make.jl once so that the Node packages are installed, then from the
# docs folder, in an environment with the package, Documenter and
# DocumenterVitepress:
#
#     julia --project=. live.jl
#
using DocumenterVitepress, FileWatching

ENV["DOCS_LIVE"] = "1"
include("make.jl")

server = run(`$(DocumenterVitepress.node()) node_modules/vitepress/bin/vitepress.js dev build/.documenter`;
    wait = false)
try
    while true
        watch_folder("src")
        sleep(0.5)
        try
            include("make.jl")
        catch e
            showerror(stderr, e)
            println(stderr)
        end
    end
finally
    unwatch_folder("src")
    kill(server)
end
