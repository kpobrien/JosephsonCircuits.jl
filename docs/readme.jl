# Execute the actual README quickstart, without installation commands or a
# second copy of its circuit. Called by both documentation build entry points;
# also runnable alone with `julia --project=. docs/readme.jl`.
let
    path = joinpath(@__DIR__, "..", "README.md")
    pattern = r"(?s)<!-- readme-example:start -->\s*```julia\n(.*?)\n```\s*<!-- readme-example:end -->"
    blocks = collect(eachmatch(pattern, read(path, String)))
    length(blocks) == 1 || error("README must contain exactly one marked Julia quickstart")
    Base.include_string(Module(:ReadmeExample), only(blocks).captures[1], path)
end
