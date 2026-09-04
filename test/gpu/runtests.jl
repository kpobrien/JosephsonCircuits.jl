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
            method = NewtonKrylov())
        rb = hbnlsolve((w1,), (8,), src1, circuit, defs;
            method = NewtonKrylov(), backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
    end

    @testset "hbnlsolve newtonkrylov, two tone" begin
        ra = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = NewtonKrylov())
        rb = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = NewtonKrylov(),
            backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
    end

    @testset "recycled deflation on the device" begin
        ra = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = NewtonKrylov())
        for form in (:adef1, :adef2, :floquet)
            pre = form === :floquet ? Floquet(size = 8, harvest = 2) :
                Recycling(size = 8, harvest = 2, form = form)
            rb = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
                dc = true, odd = true, even = true, backend = CUDABackend(),
                method = NewtonKrylov(preconditioner = pre, escalate = false))
            @test rb.solverinfo.converged
            @test agree(ra.S, rb.S)
            kr = rb.solverinfo.stages[1].krylov
            # the pair was built and applied on the device
            @test any(k -> k.deflationsize > 0, kr)
            @test kr[end].deflationrebuilds > 0
            @test kr[end].deflationproducts > 0
        end
    end

    @testset "method = Staged() on the device" begin
        ra = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = Newton())
        rb = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = Staged(),
            backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
    end

    @testset "direct current on the device" begin
        # the explicit direct current block: the average voltages appended
        # to the state, the window gathered by index on the device, and the
        # subsystem solved there. The chain is one floating static flux
        # component with a resistor to ground at each end, so a direct
        # current into port 1 develops a voltage across the pair.
        srcdc = [(mode=(1,), port=1, current=2.0e-6),
                 (mode=(0,), port=1, current=1.0e-7)]
        ra = hbnlsolve((w1,), (8,), srcdc, circuit, defs;
            dc = true, odd = true, even = true, method = NewtonKrylov())
        rb = hbnlsolve((w1,), (8,), srcdc, circuit, defs;
            dc = true, odd = true, even = true, method = NewtonKrylov(),
            backend = CUDABackend())
        @test rb.solverinfo.converged
        @test agree(ra.S, rb.S)
        @test agree(ra.dcnodevoltage, rb.dcnodevoltage; rtol = 1e-8)
        @test maximum(abs, ra.dcnodevoltage) > 0

        # and a scattering block whose zero frequency current is one of the
        # subsystem's unknowns
        R = 100.0
        blk = ScatteringParameters(
            w -> JosephsonCircuits.ABCDtoS(
                JosephsonCircuits.ABCD_seriesZ(10.0 + 0im));
            nports = 2, grounded = true, noise = Lossless())
        cblk = Circuit(
            [:p1 => Port(1; Z0 = R), :x => blk, :jj => JosephsonJunction(500e-12),
             :r2 => Resistor(R), :c1 => Capacitor(1e-12), :c2 => Capacitor(1e-12)],
            [[(:p1,1),(:x,1),(:c1,1)], [(:x,2),(:r2,1),(:jj,1),(:c2,1)],
             [(:p1,2),(:r2,2),(:jj,2),(:c1,2),(:c2,2), Ground]])
        srcb = [(mode=(1,), port=1, current=1.0e-6),
                (mode=(0,), port=1, current=1.0e-7)]
        ba = hbnlsolve((w1,), (4,), srcb, cblk, Dict{Any,Any}();
            dc = true, odd = true, even = true, method = NewtonKrylov())
        bb = hbnlsolve((w1,), (4,), srcb, cblk, Dict{Any,Any}();
            dc = true, odd = true, even = true, method = NewtonKrylov(),
            backend = CUDABackend())
        @test bb.solverinfo.converged
        @test agree(ba.S, bb.S)
        @test agree(ba.dcnodevoltage, bb.dcnodevoltage; rtol = 1e-8)
    end

    @testset "a cached sweep on the device" begin
        # the system, the preconditioner and the Krylov vectors are built on
        # the device once and rebound to each point; every point agrees
        # with a fresh device solve and with the host
        make(; Lj, Cg) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0),
            ("Lj1","1","2",Lj), ("C1","1","0",Cg),
            ("Lj2","2","3",Lj), ("C2","2","0",Cg),
            ("Lj3","3","4",Lj), ("C3","3","0",Cg),
            ("C4","4","0",Cg), ("R2","4","0",50.0)]
        p0 = (Lj = 100e-12, Cg = 40e-15)
        cache = hbcache((w1,), (8,), src1, make, p0;
            backend = CUDABackend(), ftol = 1e-10)
        hbsolve!(cache, p0)
        @test cache.converged
        pm = cache.reuse.sys.phimatrix
        for (Lj, Cg) in ((105e-12, 40e-15), (95e-12, 44e-15), (110e-12, 38e-15))
            s = hbsolve!(cache, (Lj = Lj, Cg = Cg))
            @test cache.converged
            @test cache.reuse.sys.phimatrix === pm
            fresh = hbnlsolve((w1,), (8,), src1, make(; Lj = Lj, Cg = Cg);
                backend = CUDABackend(), ftol = 1e-10, keyedarrays = false)
            host = hbnlsolve((w1,), (8,), src1, make(; Lj = Lj, Cg = Cg);
                ftol = 1e-10, keyedarrays = false)
            @test agree(s.nodeflux, fresh.nodeflux; rtol = 1e-7)
            @test agree(s.nodeflux, host.nodeflux; rtol = 1e-7)
        end
    end

    @testset "hbsolve pipeline (nonlinear + linearized sweep)" begin
        ws = 2*pi*(4.55:0.3:5.5)*1e9
        ra = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs)
        rb = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs;
            backend = CUDABackend())
        @test agree(ra.linearized.S, rb.linearized.S)
    end
    @testset "block factorization on the device" begin
        # the dense node-block factorization of the full Jacobian and of a
        # cluster mask, in double and single precision, against the host
        ra = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
            dc = true, odd = true, even = true, method = NewtonKrylov())
        Nm = hbnlsolve((w1,w2), (8,4), src2, circuit, defs; dc = true,
            odd = true, even = true, returnsystem = true).Nmodes
        mask = Matrix{Bool}(I, Nm, Nm)
        for (a, b) in ((1, 2), (2, 3), (1, 3), (4, 5))
            mask[a, b] = mask[b, a] = true
        end
        for (pre, exact) in ((FullJacobian(factorization = BlockFactorization()), true),
                (FullJacobian(factorization = BlockFactorization(; precision = Float32)), false),
                (CouplingMask(mask; factorization = BlockFactorization()), false),
                (Clusters(factorization = BlockFactorization()), false),
                (Clusters(factorization = CUDSSFactorization()), false),
                (Automatic(), false))
            rb = hbnlsolve((w1,w2), (8,4), src2, circuit, defs;
                dc = true, odd = true, even = true, backend = CUDABackend(),
                method = NewtonKrylov(preconditioner = pre))
            @test rb.solverinfo.converged
            @test agree(ra.S, rb.S)
            kr = rb.solverinfo.stages[1].krylov
            # an exact solve in double takes one Arnoldi step per Newton step
            exact && @test all(k -> k.iterations <= 2, kr)
        end
        @test JosephsonCircuits.freememory(CUDABackend()) > 0
    end


end
