# The GPU test suite. Deliberately NOT part of `Pkg.test`: CUDA and CUDSS
# live in this directory's own environment so the main test suite carries
# no resolve or precompile cost for them. Run it directly on a machine
# with a CUDA device:
#
#     julia test/gpu/runtests.jl
#
# The script activates its own environment, develops the package into it
# on first run, and instantiates; afterwards it is a plain test run. It
# errors immediately if CUDA is not functional rather than silently
# testing nothing.
import Pkg
Pkg.activate(@__DIR__)
let pkgpath = normpath(joinpath(@__DIR__, "..", ".."))
    # Develop the package on first run, and RE-develop if the manifest has
    # drifted to a registered version instead of this checkout (having the
    # name in the project is not enough).
    manifest = joinpath(@__DIR__, "Manifest.toml")
    devved = isfile(manifest) && let
        entry = get(get(Pkg.TOML.parsefile(manifest), "deps", Dict()),
            "JosephsonCircuits", nothing)
        entry !== nothing && haskey(entry[1], "path") &&
            rstrip(normpath(joinpath(@__DIR__, entry[1]["path"])), '/') ==
            rstrip(pkgpath, '/')
    end
    devved || Pkg.develop(Pkg.PackageSpec(path = pkgpath))
    Pkg.instantiate()
end

using JosephsonCircuits, CUDA, CUDSS, Test, LinearAlgebra
CUDA.functional() || error("CUDA is not functional on this machine; the GPU test suite needs a working CUDA device.")

isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))

# Device/host parity of everything the CUDA and CUDSS extensions cover:
# the same solves on CUDABackend() and CPU() must agree to solver
# tolerance. The circuits are small so the suite is dominated by one-time
# GPU compilation, not by the solves.
@testset verbose = true "GPU device/host parity" begin
    circuit, defs = testchaincircuit()
    w1 = 2*pi*5.0e9; w2 = 2*pi*1.19e9
    src1 = [(mode=(1,), port=1, current=2.0e-6)]
    src2 = [(mode=(1,0), port=1, current=1.0e-6),
            (mode=(0,1), port=1, current=0.5e-6)]

    agree(a, b; rtol = 1e-6) = norm(vec(Array(a)) - vec(Array(b))) <=
        rtol*norm(vec(Array(a)))

    @testset "hbnlsolve newtonkrylov, single tone" begin
        ra = hbnlsolve((w1,), (8,), src1, circuit, defs;
            method = :newtonkrylov)
        rb = hbnlsolve((w1,), (8,), src1, circuit, defs;
            method = :newtonkrylov, backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
    end

    @testset "hbnlsolve newtonkrylov, two tone" begin
        ra = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = :newtonkrylov)
        rb = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = :newtonkrylov,
            backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
    end

    @testset "recycled deflation on the device" begin
        ra = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = :newtonkrylov)
        rb = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = :newtonkrylov,
            krylovrecycle = 8, backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
    end

    @testset "method = :staged on the device" begin
        ra = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = :newton)
        rb = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = :staged,
            backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
    end

    @testset "hbsolve pipeline (nonlinear + linearized sweep)" begin
        ws = 2*pi*(4.55:0.3:5.5)*1e9
        ra = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs)
        rb = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs;
            backend = CUDABackend())
        @test agree(ra.linearized.S, rb.linearized.S)
    end
end
