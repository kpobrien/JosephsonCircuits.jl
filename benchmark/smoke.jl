# Every benchmark script, in a fresh process, at a small size. Nothing is
# measured: a script which no longer runs against the package's current
# interfaces is the failure this catches.
#
#     julia --project=. benchmark/smoke.jl
const PROJECT = Base.active_project()

for file in ("frontend-followup.jl", "circuit-preparation.jl")
    withenv("OPENBLAS_NUM_THREADS" => "1") do
        run(`$(Base.julia_cmd()) --startup-file=no --project=$(PROJECT)
            --threads=1 $(joinpath(@__DIR__, file)) smoke --smoke`)
    end
end
