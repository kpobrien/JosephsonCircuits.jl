# One thread against four: the same solves run in two fresh processes, one
# started with a single thread and one with four, and every number they
# record must agree. The solver's own tests run in whichever process the
# suite was started in, so nothing else in the tree compares the two.
#
#     julia --project=. test/threaded/compare.jl
using Serialization
using Test

# the subprocesses run in whichever environment this process was started
# in, so the suite can be run from the package's own project or from a
# development environment which has it
const PROJECT = Base.active_project()

mktempdir() do dir
    for nthreads in (1, 4)
        cmd = Cmd([first(Base.julia_cmd().exec), "--startup-file=no",
            "--check-bounds=yes", "--project=$(PROJECT)",
            "--threads=$(nthreads)", joinpath(@__DIR__, "runtests.jl"),
            joinpath(dir, "$(nthreads).bin"), string(nthreads)])
        withenv("OPENBLAS_NUM_THREADS" => "1", "JULIA_TEST_WORKERS" => "0") do
            run(cmd)
        end
    end
    one = deserialize(joinpath(dir, "1.bin"))
    four = deserialize(joinpath(dir, "4.bin"))
    @testset "one thread against four" begin
        for field in keys(one)
            if field === :fluxes
                @test all(isapprox(a, b; rtol = 1e-8)
                    for (a, b) in zip(one[field], four[field]))
            else
                @test isapprox(one[field], four[field];
                    rtol = 1e-8, atol = 1e-15)
            end
        end
    end
end
