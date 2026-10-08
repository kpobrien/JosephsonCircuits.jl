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
using JosephsonCircuits: Clusters, CouplingMask, Floquet, FullJacobian
CUDA.functional() || error("CUDA is not functional on this machine; the GPU test suite needs a working CUDA device.")

isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))
# the I/Q and quantum test functions, taking the backend; the constant keeps
# their files from running their CPU suites here
const TRANSIENTBACKENDTESTS = true
include(joinpath(@__DIR__, "..", "transient", "iq.jl"))
include(joinpath(@__DIR__, "..", "transient", "quantum.jl"))

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

    @testset "a circuit without junctions on the device" begin
        # the junction kernels have nothing to run over: a linear circuit's
        # pump solve and sweep on the device, against the host
        resonator, rdefs = testresonatorcircuit()
        resonatorS(b) = hbsolve([0.9, 1.1], (2.3,), src1, (2,), (4,), resonator,
            rdefs; keyedarrays = false, backend = b).linearized.S
        @test agree(resonatorS(CUDABackend()), resonatorS(JosephsonCircuits.CPU());
            rtol = 1e-10)
    end

    @testset "a polynomial current-phase relation" begin
        # the polynomial relations are evaluated by Horner in whole array
        # broadcasts along the junction axis, and in a circuit which
        # mixes the two kinds each relation over the view of its own
        # junctions: both paths run on the device as they do on the host
        L0 = 1e-9
        w = 2*pi*4.75e9
        src = [(mode = (1,), port = 1, current = 4*0.00565e-6)]
        taylor(n) = [k % 2 == 1 ? (-1.0)^((k-1)÷2)/prod(1.0:k) : 0.0 for k in 1:n]
        poly = NonlinearInductor(L0, PolynomialCPR([1.0, 0.25, -1/6, 0.0, 1/120]))
        jpa = Circuit([("p1","1","0",Port(1)), ("c1","1","2",Capacitor(100e-15)),
            ("lj","2","0",poly), ("c2","2","0",Capacitor(1e-12))])
        # one junction of each kind in one circuit, which takes the mixed
        # path
        mixed = Circuit([("p1","1","0",Port(1)), ("c1","1","2",Capacitor(100e-15)),
            ("jj","2","3",JosephsonJunction(L0)),
            ("lj","3","0",NonlinearInductor(L0, PolynomialCPR(taylor(9)))),
            ("c2","2","0",Capacitor(1e-12))])
        for c in (jpa, mixed)
            ra = hbnlsolve((w,), (8,), src, c, Dict{Symbol,Float64}();
                method = NewtonKrylov(), keyedarrays = false)
            rb = hbnlsolve((w,), (8,), src, c, Dict{Symbol,Float64}();
                method = NewtonKrylov(), backend = CUDABackend(),
                keyedarrays = false)
            @test rb.solverinfo.converged
            @test agree(ra.nodeflux, rb.nodeflux; rtol = 1e-10)
        end
    end

    @testset "the operating point of a device solve" begin
        # an operating point is a host object: a device solve hands it a
        # host twin of its system, relations included, on which its
        # residual vanishes at the point, with the Jacobian the host solve
        # assembles there; the sensitivities through it are the host's
        L0 = 1e-9
        w = 2*pi*4.75e9
        src = [(mode = (1,), port = 1, current = 4*0.00565e-6)]
        poly = NonlinearInductor(L0, PolynomialCPR([1.0, 0.25, -1/6, 0.0, 1/120]))
        jpa = Circuit([("p1","1","0",Port(1)), ("c1","1","2",Capacitor(100e-15)),
            ("lj","2","0",poly), ("c2","2","0",Capacitor(1e-12))])
        opat(; kw...) = hbnlsolve((w,), (8,), src, jpa, Dict{Symbol,Float64}();
            keyedarrays = false, returnoperatingpoint = true, atol = 1e-12,
            kw...).operatingpoint
        oa = opat()
        ob = opat(backend = CUDABackend())
        @test ob.sys.phimatrix isa Array
        @test agree(oa.jacobian, ob.jacobian; rtol = 1e-8)
        xr = JosephsonCircuits.complex_to_real(ob.x, ob.modelayout.isreal)
        JosephsonCircuits.setpoint!(ob.sys, xr)
        F = JosephsonCircuits.residual!(similar(xr), ob.sys)
        @test norm(F) < 1e-10
        ws = 2*pi*(4.55:0.3:5.5)*1e9
        # both contraction orders: the reverse one solves the transposed
        # pump Jacobian, which a device sweep factorizes on the host
        kw = (; returnSsensitivity = true, sensitivitynames = ["Lj1", "C2"],
            keyedarrays = false)
        sa = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs; kw...)
        for mode in (:forward, :reverse)
            sb = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs;
                backend = CUDABackend(), sensitivitymode = mode, kw...)
            @test agree(sa.linearized.Ssensitivity, sb.linearized.Ssensitivity)
        end
        # a single precision device pump is differentiated on a double
        # precision host system, in either order, to what single precision
        # holds of the point
        single(m) = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs;
            backend = CUDABackend(), method = NewtonKrylov(precision = Float32),
            sensitivitymode = m, kw...).linearized.Ssensitivity
        forward = single(:forward)
        @test agree(single(:reverse), forward; rtol = 1e-10)
        @test agree(sa.linearized.Ssensitivity, forward; rtol = 1e-2)
    end

    @testset "block parameters with a device pump" begin
        # a scattering block which states its derivative feeding a JPA: the
        # pump solved on the device, the sweep on the host, which takes a
        # block sensitivity, the block's parameter in either order against
        # the host
        Z0 = 50.0
        seriesS(theta) = w -> (z = 1/(im*w*theta*Z0); [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
        dS(theta) = w -> (z = 1/(im*w*theta*Z0); d = -2*z/(theta*(z+2)^2); [d -d; -d d])
        blk = ScatteringParameters(seriesS(100e-15); nports = 2, grounded = false, derivatives = (theta = dS(100e-15),))
        bjpa = Circuit([:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :cc => blk,
                :jj => JosephsonJunction(:Lj), :c2 => Capacitor(1000.0e-15)],
            [((:p1, 1), (:r1, 1), (:cc, 1, 1)), ((:cc, 2, 1), (:jj, 1), (:c2, 1)),
             ((:cc, 1, 2), (:cc, 2, 2), (:jj, 2), (:c2, 2), (:r1, 2), (:p1, 2), Ground)])
        bdefs = Dict(:Lj => 1000.0e-12)
        bws, bwp = 2pi .* [4.5e9, 4.6e9], (2pi*4.75001e9,)
        bsrc = [(mode = (1,), port = 1, current = 0.00565e-6)]
        host = designsensitivities(bjpa, bdefs, bws, bwp, bsrc, (2,), (8,); sensitivitymode = :forward)
        for mode in (:forward, :reverse)
            device = designsensitivities(bjpa, bdefs, bws, bwp, bsrc, (2,), (8,); sensitivitymode = mode,
                backend = CUDABackend())
            @test agree(host.dSdp, device.dSdp; rtol = 1e-8)
        end
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
        let pre = Floquet(size = 8, harvest = 2)
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

    @testset "the direct methods' factorization on the device" begin
        # with no factorization given, the direct methods factorize where
        # the solve runs
        ra = hbnlsolve((w1,), (8,), src1, circuit, defs; method = Newton())
        for m in (Newton(), QuasiNewton())
            rb = hbnlsolve((w1,), (8,), src1, circuit, defs; method = m,
                backend = CUDABackend())
            @test rb.solverinfo.converged
            @test agree(ra.S, rb.S)
        end
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

        # five of the blocks in a chain, whose subsystem is larger than the
        # one work item solve takes and is solved on the host between a
        # gather and a scatter on the device
        comps = Pair{Symbol,Any}[:p1 => Port(1; Z0 = R)]
        for k in 1:5
            push!(comps, Symbol(:x, k) => blk, Symbol(:c, k) => Capacitor(1e-12))
        end
        push!(comps, :jj => JosephsonJunction(500e-12), :r2 => Resistor(R),
            :cend => Capacitor(1e-12))
        nets = Any[[(:p1,1), (:x1,1), (:c1,1)]]
        for k in 1:4
            push!(nets, [(Symbol(:x, k),2), (Symbol(:x, k + 1),1), (Symbol(:c, k + 1),1)])
        end
        push!(nets, [(:x5,2), (:r2,1), (:jj,1), (:cend,1)])
        push!(nets, vcat(Any[(:p1,2), (:r2,2), (:jj,2), (:cend,2)],
            [(Symbol(:c, k),2) for k in 1:5], [Ground]))
        chain = Circuit(comps, nets)
        ca = hbnlsolve((w1,), (4,), srcb, chain, Dict{Any,Any}();
            dc = true, odd = true, even = true, method = NewtonKrylov())
        cb = hbnlsolve((w1,), (4,), srcb, chain, Dict{Any,Any}();
            dc = true, odd = true, even = true, method = NewtonKrylov(),
            backend = CUDABackend())
        @test cb.solverinfo.converged
        @test agree(ca.S, cb.S)
        @test agree(ca.dcnodevoltage, cb.dcnodevoltage; rtol = 1e-8)
    end

    @testset "a cached sweep on the device" begin
        # the system, the preconditioner and the Krylov vectors are built on
        # the device once and rebound to each point; every point agrees
        # with a fresh device solve and with the host
        chain3 = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:Lj1, 1, 2, JosephsonJunction(:Lj)), (:C1, 1, 0, Capacitor(:Cg)),
            (:Lj2, 2, 3, JosephsonJunction(:Lj)), (:C2, 2, 0, Capacitor(:Cg)),
            (:Lj3, 3, 4, JosephsonJunction(:Lj)), (:C3, 3, 0, Capacitor(:Cg)),
            (:C4, 4, 0, Capacitor(:Cg)), (:R2, 4, 0, Resistor(50.0))])
        p0 = (Lj = 100e-12, Cg = 40e-15)
        cache = hbcache((w1,), (8,), src1, chain3, Dict(pairs(p0));
            backend = CUDABackend(), atol = 1e-10)
        hbsolve!(cache, p0)
        @test cache.converged
        pm = cache.reuse.sys.phimatrix
        for (Lj, Cg) in ((105e-12, 40e-15), (95e-12, 44e-15), (110e-12, 38e-15))
            s = hbsolve!(cache, (Lj = Lj, Cg = Cg))
            @test cache.converged
            @test cache.reuse.sys.phimatrix === pm
            fresh = hbnlsolve((w1,), (8,), src1, chain3, Dict(:Lj => Lj, :Cg => Cg);
                backend = CUDABackend(), atol = 1e-10, keyedarrays = false)
            host = hbnlsolve((w1,), (8,), src1, chain3, Dict(:Lj => Lj, :Cg => Cg);
                atol = 1e-10, keyedarrays = false)
            @test agree(s.nodeflux, fresh.nodeflux; rtol = 1e-7)
            @test agree(s.nodeflux, host.nodeflux; rtol = 1e-7)
        end
        # the direct methods rebind the system and refactorize their
        # Jacobian on the device too
        for m in (Newton(), QuasiNewton())
            c = hbcache((w1,), (8,), src1, chain3, Dict(pairs(p0));
                backend = CUDABackend(), atol = 1e-10, method = m)
            hbsolve!(c, p0)
            held = c.reuse.sys.phimatrix
            s = hbsolve!(c, (Lj = 105e-12, Cg = 40e-15))
            host = hbnlsolve((w1,), (8,), src1, chain3,
                Dict(:Lj => 105e-12, :Cg => 40e-15); atol = 1e-10,
                keyedarrays = false, method = m)
            @test c.converged
            @test c.reuse.sys.phimatrix === held
            @test agree(s.nodeflux, host.nodeflux; rtol = 1e-7)
        end
        # the problem interface and an external solver, which it is handed
        # to, live on the host, and refuse a device backend
        @test_throws ArgumentError hbnonlinearproblem((w1,), (8,), src1,
            circuit, defs; backend = CUDABackend())
        @test_throws ArgumentError hbnlsolve((w1,), (8,), src1, circuit, defs;
            backend = CUDABackend(),
            method = ExternalSolver((prob, u0) -> (u0, false)))
    end

    @testset "hbsolve pipeline (nonlinear + linearized sweep)" begin
        ws = 2*pi*(4.55:0.3:5.5)*1e9
        ra = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs)
        rb = hbsolve(ws, (w1,), src1, (2,), (8,), circuit, defs;
            backend = CUDABackend())
        @test agree(ra.linearized.S, rb.linearized.S)
        # the block factorization on the device: both directions from one
        # factorization per frequency, every output, exact and refined
        # single precision factors
        kw = (; dc = true, threewavemixing = true, fourwavemixing = true,
            returnSnoise = true, keyedarrays = false)
        wl = 2*pi*collect(range(4.41e9, 5.57e9, length = 4))
        rc = hbsolve(wl, (w1,w2), src2, (2,2), (8,4), circuit, defs;
            factorization = KLUfactorization(), kw...)
        for f in (BlockFactorization(), BlockFactorization(precision = Float32))
            rd = hbsolve(wl, (w1,w2), src2, (2,2), (8,4), circuit, defs;
                backend = CUDABackend(), factorization = f, kw...)
            for name in (:S, :Snoise, :QE, :CM)
                @test agree(getfield(rc.linearized, name),
                    getfield(rd.linearized, name))
            end
        end
        # the automatic choice on the device: the block factorization for
        # two tones, agreeing with KLU on the host
        auto = hbsolve(wl, (w1,w2), src2, (2,2), (8,4), circuit, defs;
            backend = CUDABackend(), kw...)
        @test agree(rc.linearized.S, auto.linearized.S)
        @test agree(rc.linearized.QE, auto.linearized.QE)
        Asp = JosephsonCircuits.SparseArrays.sparse(
            ComplexF64[1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1])
        @test JosephsonCircuits.linearizedfactorization(Asp, 2, 2,
            CUDABackend()) isa BlockFactorization
        @test JosephsonCircuits.linearizedfactorization(Asp, 2, 2,
            CUDABackend(); budget = 0) isa CUDSSFactorization
        @test JosephsonCircuits.linearizedfactorization(Asp, 2, 1,
            CUDABackend()) isa CUDSSFactorization
        # single precision solutions on the device
        rs = hbsolve(wl, (w1,w2), src2, (2,2), (8,4), circuit, defs;
            backend = CUDABackend(),
            factorization = BlockFactorization(precision = Float32, refine = false),
            kw...)
        for (name, tol) in ((:S, 2e-3), (:Snoise, 1e-3), (:QE, 1e-3), (:CM, 1e-4))
            @test isapprox(getfield(rc.linearized, name),
                getfield(rs.linearized, name); rtol = tol)
        end
    end
    @testset "a finished sweep leaves its memory to the next" begin
        JC = JosephsonCircuits
        # what CUDA.jl's pool holds without an array in it counts as free
        # memory, since an allocation takes it first
        GC.gc(); CUDA.reclaim()
        free = JC.freememory(CUDABackend())
        x = CUDA.zeros(UInt8, 2^30)
        CUDA.unsafe_free!(x)
        @test abs(JC.freememory(CUDABackend()) - free) < 2^26
        # a sweep returns the arrays of its batches to the pool when it
        # ends, rather than leaving them to the collector, so that the next
        # sweep or solve sizes itself against them: after it the pool holds
        # in arrays less than a tenth of what its batch held, the block
        # factors of every system at two tones, and at one the values,
        # solutions and right-hand sides of both directions of a cuDSS
        # batch (cuDSS allocates its factors itself)
        chain, chaindefs = testchaincircuit(20)
        ws = 2*pi*collect(range(4.4e9, 5.6e9, length = 400))
        for (wp, src, Nmod) in (((w1, w2), src2, (3, 3)), ((w1,), src1, (6,)))
            Npump = length(wp) == 2 ? (8, 4) : (8,)
            nl = hbnlsolve(wp, Npump, src, chain, chaindefs;
                keyedarrays = false, backend = CUDABackend())
            sweep(w; kw...) = hblinsolve(w, chain, chaindefs; nonlinear = nl,
                Nmodulationharmonics = Nmod, keyedarrays = false, kw...)
            d = sweep(ws[1:1]; debuglsys = true)
            A = d.lsys.Asparse
            held = if length(wp) == 2
                length(ws)*JC.blocksystembytes(ComplexF64,
                    JC.blocksymbolic(A, d.Nmodes))
            else
                nrhs = size(d.bnm, 2)
                2*JC.uniformbatchlimit(nrhs)*16*(JC.SparseArrays.nnz(A) +
                    2*size(A, 1)*nrhs)
            end
            sweep(ws; backend = CUDABackend())
            GC.gc(); CUDA.reclaim()
            used = CUDA.used_memory()
            sweep(ws; backend = CUDABackend())
            @test CUDA.used_memory() - used < held ÷ 10
        end
    end
    @testset "a lossy coupled inductor on the device" begin
        # a complex coupled inductance is conjugated at the negative
        # frequency modes, the idlers of the sweep, on the device as on the
        # host
        c = Circuit(Any[(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:L1, 1, 3, Inductor(1e-9*(1 - 0.05im))),
            (:C1, 3, 2, Capacitor(100e-15)),
            (:Lj1, 2, 0, JosephsonJunction(1000e-12)),
            (:C2, 2, 0, Capacitor(1000e-15)),
            (:L2, 4, 0, Inductor(1e-9)), (:C4, 4, 0, Capacitor(100e-15)),
            (:K1, :L1, :L2, MutualInductor(0.3))])
        ws = 2*pi*(4.5:0.25:5.0)*1e9
        sweep(; kw...) = hbsolve(ws, (2*pi*4.75001e9,),
            [(mode = (1,), port = 1, current = 0.00565e-6)], (2,), (8,), c;
            keyedarrays = false, returnQE = false, returnCM = false,
            returnnbar = false,
            kw...).linearized
        ra = sweep()
        rb = sweep(backend = CUDABackend())
        @test agree(ra.S, rb.S)
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
            if exact
                @test all(k -> k.iterations <= 2, kr)
            end
        end
        # a singular diagonal block in one system of a batch throws, naming
        # that system, once the batch is factorized: two uncoupled node
        # blocks, the first singular in the second of three systems
        SA = JosephsonCircuits.SparseArrays
        A = SA.sparse(ComplexF64[2 1 0 0; 1 2 0 0; 0 0 2 1; 0 0 1 2])
        vals = repeat(SA.nonzeros(A), 1, 3)
        vals[1:4, 2] .= 1
        F = JosephsonCircuits.blockanalysis(BlockFactorization(), A;
            blocksize = 2, backend = CUDABackend(), nb = 3)
        err = try
            JosephsonCircuits.fillandfactorize!(F, CuArray(vals))
            nothing
        catch e
            e
        end
        @test err isa SingularException && err.info == 2
        @test JosephsonCircuits.freememory(CUDABackend()) > 0
        # a node block singular in a nonsingular matrix (test/hbsolve.jl): the
        # block factorization the device sweep chose itself gives way to
        # cuDSS, and one given explicitly throws
        resonator, rdefs = testresonatorcircuit()
        resonatorS(; kw...) = hbsolve([1.0], (2.3, 3.7), src2, (1,1), (1,1), resonator,
            rdefs; keyedarrays = false, kw...).linearized.S
        @test agree(resonatorS(backend = CUDABackend()),
            resonatorS(factorization = KLUfactorization()); rtol = 1e-10)
        @test_throws SingularException resonatorS(backend = CUDABackend(),
            factorization = BlockFactorization())
    end

    @testset "single precision cuDSS factors as the preconditioner" begin
        # A preconditioner's factors only have to make the Krylov solve
        # converge, not carry the accuracy of the answer, so they may be
        # held in single precision while the iteration stays double;
        # `applypreconditioner!` converts the residual down and the
        # correction back. cuDSS honours a precision, as the block
        # factorization does: KLU and UMFPACK are compiled for double and
        # promote a single precision matrix.
        ra = hbnlsolve((w1,), (8,), src1, circuit, defs; method = NewtonKrylov())
        its = Int[]
        for prec in (nothing, Float32)
            f = CUDSSFactorization(; precision = prec)
            @test JosephsonCircuits.factorizationprecision(f) === prec
            rb = hbnlsolve((w1,), (8,), src1, circuit, defs;
                backend = CUDABackend(),
                method = NewtonKrylov(preconditioner = FullJacobian(f)))
            @test rb.solverinfo.converged
            # the preconditioner does not enter the answer, only the path
            @test agree(ra.S, rb.S)
            push!(its, sum(st.iterations for st in rb.solverinfo.stages))
        end
        # and it does not cost the solve any Newton steps
        @test its[2] <= its[1] + 2

        # single precision factors of a double precision iteration are a
        # preconditioner, not an exact solve, and their escalation
        # refactorizes the same coupling set in the iteration's precision,
        # after which the preconditioner solve is exact
        d = hbnlsolve((w1,), (8,), src1, circuit, defs;
            backend = CUDABackend(), debugJacobian = true)
        pc = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Rbnm, d.Nmodes,
            d.Nbranches, d.Nfreq, d.modelayout; Amatrixmodes = d.Amatrixmodes,
            spec = FullJacobian(CUDSSFactorization(precision = Float32)))
        @test !JosephsonCircuits.isexactpreconditioner(pc)
        x = CuArray(0.3*randn(length(d.xr)))
        r = CuArray(randn(length(d.xr)))
        z = similar(r); Jz = similar(r)
        # the relative residual of the preconditioner solve against the
        # exact matrix-free product at the same point
        function pcresidual()
            JosephsonCircuits.updatepreconditioner!(pc, x)
            JosephsonCircuits.applypreconditioner!(z, pc, r)
            JosephsonCircuits.jacobianvectorproduct!(Jz, d.sys, z)
            return norm(Jz - r)/norm(r)
        end
        single = pcresidual()
        @test single < 1e-3
        @test JosephsonCircuits.escalatepreconditioner!(pc)
        @test JosephsonCircuits.isexactpreconditioner(pc)
        @test pc.factorization.precision === Float64
        double = pcresidual()
        @test double < 1e-9
        @test double < single/100
    end

    # the circuit in time on the device: the solve through the device
    # sparse products, the assembly plan and cuDSS, its tangent and its
    # adjoint, against the host
    @testset "a polynomial relation in the transient on the device" begin
        # the transient reads the same table: its residual and Jacobian on
        # the device, its tangent and adjoint, and a batch, on a circuit
        # holding one polynomial junction and one Josephson one
        taylor(n) = [k % 2 == 1 ? (-1.0)^((k-1)÷2)/prod(1.0:k) : 0.0 for k in 1:n]
        circuit = Circuit([
            ("p1", "1", "0", Port(1)), ("c1", "1", "2", Capacitor(100e-15)),
            ("jj", "2", "3", JosephsonJunction(1e-9)),
            ("lj", "3", "0", NonlinearInductor(1e-9, PolynomialCPR([1.0, 0.25, -1/6]))),
            ("c2", "2", "0", Capacitor(1e-12))])
        drive(t) = t <= 0 ? 0.0 : 0.3e-6*sinpi(2*3e9*t)
        prob = transientproblem(circuit; sources = [TransientSource(1, drive)])
        args = ((0.0, 0.5e-9),)
        host = transientsolve(prob, args...; dt = 2e-12, record = :phases,
            method = GaussLegendre(), rtol = 1e-12)
        device = transientsolve(prob, args...; dt = 2e-12, record = :phases,
            method = GaussLegendre(), rtol = 1e-12, backend = CUDABackend())
        @test agree(host.outgoing, device.outgoing; rtol = 1e-8)
        currents = [1e-8*sinpi(2*1.1e9*t) for _ in 1:1, t in host.times]
        @test agree(transienttangent(host, currents).outgoing,
            transienttangent(device, currents).outgoing; rtol = 1e-8)
        weights = [cospi(2*1.7e9*t) for _ in 1:1, t in host.times]
        @test agree(transientadjoint(host, weights).currents,
            transientadjoint(device, weights).currents; rtol = 1e-8)
        # a batch, whose stage Jacobian averages the derivative of the
        # relation over the two stages
        half = transientproblem(prob;
            sources = [TransientSource(1, t -> drive(t)/2)])
        bh = transientsolve([prob, half], args...; dt = 2e-12, record = :phases,
            method = GaussLegendre(), rtol = 1e-12)
        bd = transientsolve([prob, half], args...; dt = 2e-12, record = :phases,
            method = GaussLegendre(), rtol = 1e-12, backend = CUDABackend())
        @test agree(bh.outgoing, bd.outgoing; rtol = 1e-8)
    end

    @testset "transient on the device" begin
        circuit = Circuit([
            ("p1", "1", "0", Port(1)), ("p2", "2", "1", Port(2; Z0 = 75.0)),
            ("c1", "1", "0", Capacitor(1e-12)), ("c2", "2", "0", Capacitor(1e-12)),
            ("jj", "2", "1", JosephsonJunction(1e-9)), ("l", "2", "0", Inductor(2e-9))])
        drive(t) = t <= 0 ? 0.0 : 0.1e-6*sinpi(2*3e9*t)
        prob = transientproblem(circuit; sources = [TransientSource(1, drive), TransientSource(2, 2e-8)])
        host = transientsolve(prob, (0.0, 0.5e-9); dt = 2e-12, record = :states, rtol = 1e-12)
        device = transientsolve(prob, (0.0, 0.5e-9); dt = 2e-12, record = :states, rtol = 1e-12,
            backend = CUDABackend())
        for q in (:voltage, :incident, :outgoing, :flux, :rate)
            @test getproperty(device, q) isa CuArray
            @test Array(getproperty(device, q)) ≈ getproperty(host, q) rtol=1e-8 atol=1e-14
        end
        currents = [1e-8*sinpi(2*1.1e9*t + p) for p in 1:2, t in host.times]
        th = transienttangent(host, currents)
        td = transienttangent(device, currents)
        @test Array(td.outgoing) ≈ th.outgoing rtol=1e-8 atol=1e-14
        weights = [cospi(2*1.7e9*t + p) for p in 1:2, t in host.times]
        ah = transientadjoint(host, weights)
        ad = transientadjoint(device, weights)
        @test Array(ad.currents) ≈ ah.currents rtol=1e-8 atol=1e-14
        @test Array(ad.initialflux) ≈ ah.initialflux rtol=1e-8 atol=1e-14
        @test transientdemodulate(device, 2, 2pi*3e9) ≈ transientdemodulate(host, 2, 2pi*3e9) rtol=1e-8
        # the sensitivity to the component values, under both rules, from
        # the states and from checkpoints, and the adjoint's
        names = ["c1", "jj", "l"]
        for (method, record) in ((Trapezoidal(), :states), (GaussLegendre(), :states), (GaussLegendre(), :checkpoints))
            sh = transientsolve(prob, (0.0, 0.5e-9); dt = 2e-12, method, record, rtol = 1e-12)
            sd = transientsolve(prob, (0.0, 0.5e-9); dt = 2e-12, method, record, rtol = 1e-12,
                backend = CUDABackend())
            @test Array(transientsensitivity(sd, names).outgoing) ≈ transientsensitivity(sh, names).outgoing rtol=1e-7 atol=1e-14
            @test Array(transientadjoint(sd, weights; components = names).sensitivity) ≈
                transientadjoint(sh, weights; components = names).sensitivity rtol=1e-7
        end
        # the noise on a device solution against the host: a record that
        # starts at equilibrium, and bath tones on its Fourier bins
        quiet = transientproblem(circuit; sources = [TransientSource(1, drive)])
        qh = transientsolve(quiet, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12)
        qd = transientsolve(quiet, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12,
            backend = CUDABackend())
        dw = 2pi/(length(qh.times)*2e-12)
        plan = transientquantumplan(qh, qh.times, [2dw, 2dw]; ports = [1, 2])
        dplan = transientquantumplan(qh, qh.times, [2dw, 2dw]; ports = [1, 2], backend = CUDABackend())
        nh = transientnoise(qh, plan; frequencies = [2dw, 3dw], weights = fill(dw, 2))
        nd = transientnoise(qd, dplan; frequencies = [2dw, 3dw], weights = fill(dw, 2))
        @test nd.covariance ≈ nh.covariance rtol=1e-8
        @test nd.commutator ≈ nh.commutator rtol=1e-8
        # the forward method contracts the device responses with the device
        # measurement, and agrees with the adjoint on the loaded trajectory
        fd = transientnoise(qd, dplan; frequencies = [2dw, 3dw], weights = fill(dw, 2), method = :forward)
        @test fd.covariance ≈ nd.covariance rtol=1e-8
        @test fd.commutator ≈ nd.commutator rtol=1e-8
        # the Gauss-Legendre rule on the device: the complex stage
        # factorization through cuDSS, the solve, the responses and the
        # noise against the host
        gh = transientsolve(prob, (0.0, 0.5e-9); dt = 2e-12, record = :states, rtol = 1e-12, method = GaussLegendre())
        gd = transientsolve(prob, (0.0, 0.5e-9); dt = 2e-12, record = :states, rtol = 1e-12,
            method = GaussLegendre(), backend = CUDABackend())
        @test gd.phases isa CuArray
        for q in (:voltage, :outgoing, :flux, :rate, :phases)
            @test Array(getproperty(gd, q)) ≈ getproperty(gh, q) rtol=1e-8 atol=1e-14
        end
        @test Array(transienttangent(gd, currents).outgoing) ≈ transienttangent(gh, currents).outgoing rtol=1e-8 atol=1e-14
        gah, gad = transientadjoint(gh, weights), transientadjoint(gd, weights)
        @test Array(gad.currents) ≈ gah.currents rtol=1e-8 atol=1e-14
        @test Array(gad.initialflux) ≈ gah.initialflux rtol=1e-8 atol=1e-14
        gqh = transientsolve(quiet, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre())
        gqd = transientsolve(quiet, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12,
            method = GaussLegendre(), backend = CUDABackend())
        gnh = transientnoise(gqh, plan; frequencies = [2dw, 3dw], weights = fill(dw, 2))
        gnd = transientnoise(gqd, dplan; frequencies = [2dw, 3dw], weights = fill(dw, 2))
        @test gnd.covariance ≈ gnh.covariance rtol=1e-8
        @test gnd.commutator ≈ gnh.commutator rtol=1e-8
        # a batch of conditions on the uniform cuDSS batch, and one larger
        # than a chunk, against the host
        base = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0), TransientSource(2, 2e-8)])
        make(a) = transientproblem(base; sources = [TransientSource(1, let a = a; t -> t <= 0 ? 0.0 : a*sinpi(2*3e9*t); end), TransientSource(2, 2e-8)])
        problems = make.([0.05e-6, 0.1e-6, 0.2e-6])
        bh = transientsolve(problems, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12)
        bd = transientsolve(problems, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, backend = CUDABackend())
        @test bd.voltage isa CuArray
        @test Array(bd.voltage) ≈ bh.voltage rtol=1e-8 atol=1e-14
        @test Array(bd.phases) ≈ bh.phases rtol=1e-8 atol=1e-14
        @test Array(transientadjoint(bd[2], weights).currents) ≈ transientadjoint(bh[2], weights).currents rtol=1e-8 atol=1e-14
        # the noise of a batch on the device against the host
        quietbase = transientproblem(circuit; sources = [TransientSource(1, t -> 0.0)])
        qmake(a) = transientproblem(quietbase; sources = [TransientSource(1, let a = a; t -> t <= 0 ? 0.0 : a*sinpi(2*3e9*t); end)])
        qb = qmake.([0.05e-6, 0.1e-6])
        nbh = transientnoise(transientsolve(qb, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12), plan;
            frequencies = [2dw, 3dw], weights = fill(dw, 2))
        nbd = transientnoise(transientsolve(qb, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, backend = CUDABackend()), dplan;
            frequencies = [2dw, 3dw], weights = fill(dw, 2))
        @test nbd.covariance ≈ nbh.covariance rtol=1e-8
        @test nbd.commutator ≈ nbh.commutator rtol=1e-8
        nfd = transientnoise(transientsolve(qb, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, backend = CUDABackend()), dplan;
            frequencies = [2dw, 3dw], weights = fill(dw, 2), method = :forward)
        @test nfd.covariance ≈ nbh.covariance rtol=1e-8
        # a warm attenuator's channels on the device against the host
        ac = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)),
            (:att, 1, 2, ScatteringParameters([0.0 0.6; 0.6 0.0]; zref = 50.0, noise = ThermalEquilibrium(0.3))),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        ap = transientproblem(ac; sources = [TransientSource(1, t -> t <= 0 ? 0.0 : 0.05e-6*sinpi(2*3e9*t))])
        aplan = transientquantumplan(qh, qh.times, [2dw, 3dw]; ports = [1, 2])
        adplan = transientquantumplan(qh, qh.times, [2dw, 3dw]; ports = [1, 2], backend = CUDABackend())
        anh = transientnoise(transientsolve(ap, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre()), aplan;
            frequencies = [2dw, 3dw], weights = fill(dw, 2))
        and = transientnoise(transientsolve(ap, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend()), adplan;
            frequencies = [2dw, 3dw], weights = fill(dw, 2))
        @test and.covariance ≈ anh.covariance rtol=1e-8
        @test and.addedcovariance ≈ anh.addedcovariance rtol=1e-8
        @test and.commutator ≈ anh.commutator rtol=1e-8
        # a projected algebraic direction on the device against the host:
        # the junction to an inductor node without capacitance
        pc = Circuit([("p1", "1", "0", Port(1)), ("c1", "1", "0", Capacitor(1e-12)), ("jj", "1", "2", JosephsonJunction(1e-9)),
            ("l2", "2", "0", Inductor(2e-9)), ("p2", "2", "0", Port(2; termination = nothing))])
        pp = transientproblem(pc; sources = [TransientSource(1, t -> t <= 0 ? 0.0 : 0.3e-6*sinpi(2*3e9*t)), TransientSource(2, t -> 0.0)])
        ph = transientsolve(pp, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre())
        pd = transientsolve(pp, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend())
        @test Array(pd.voltage) ≈ ph.voltage rtol=1e-8 atol=1e-14
        @test Array(pd.endphases) ≈ ph.endphases rtol=1e-8
        pcurrents = [1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k) for k in 1:2, t in ph.times]
        pweights = [cospi(2*1.7e9*t + k) for k in 1:2, t in ph.times]
        @test Array(transienttangent(pd, pcurrents).outgoing) ≈ transienttangent(ph, pcurrents).outgoing rtol=1e-8
        @test Array(transientadjoint(pd, pweights).currents) ≈ transientadjoint(ph, pweights).currents rtol=1e-8
        # a lattice of junctions without capacitance, every node algebraic
        # and every endpoint projected, each condition on its own kept
        # factorization of the projection: a batch whose strongest drive
        # refreshes its factorization as the junctions' stiffness moves,
        # the tangent, the adjoint and both noise methods through the
        # projection, and the trapezoidal rule's endpoint and forward noise,
        # against the host; and forty steps of the device's stepper on one
        # factorization of its projection
        xsite(i, j) = "n$(i)_$(j)"
        xlattice = Any[("P1", xsite(1, 1), "0", Port(1; Z0 = 50.0)), ("P2", xsite(4, 4), "0", Port(2; Z0 = 50.0))]
        for i in 1:4, j in 1:4
            i < 4 && push!(xlattice, ("Jv$(i)_$(j)", xsite(i, j), xsite(i + 1, j), JosephsonJunction(100e-12)))
            j < 4 && push!(xlattice, ("Jh$(i)_$(j)", xsite(i, j), xsite(i, j + 1), JosephsonJunction(100e-12)))
        end
        xfree = transientproblem(Circuit(xlattice); sources = [TransientSource(1, t -> 0.0)])
        xmake(a) = transientproblem(xfree; sources = [TransientSource(1, let a = a; t -> t <= 0 ? 0.0 : a*sinpi(1e10*t); end)])
        xbatch = xmake.([2e-6, 6e-6, 12e-6])
        xh = transientsolve(xbatch, (0.0, 100e-12); dt = 1e-12, record = :phases)
        xd = transientsolve(xbatch, (0.0, 100e-12); dt = 1e-12, record = :phases, backend = CUDABackend())
        @test Array(xd.voltage) ≈ xh.voltage rtol=1e-8 atol=1e-14
        @test Array(xd.phases) ≈ xh.phases rtol=1e-8 atol=1e-14
        xtimes = xh[1].times
        xcurrents = [1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k) for k in 1:2, t in xtimes]
        xweights = [cospi(2*1.7e9*t + k) for k in 1:2, t in xtimes]
        @test Array(transienttangent(xd[3], xcurrents).outgoing) ≈ transienttangent(xh[3], xcurrents).outgoing rtol=1e-8 atol=1e-14
        @test Array(transientadjoint(xd[3], xweights).currents) ≈ transientadjoint(xh[3], xweights).currents rtol=1e-8 atol=1e-14
        xdw = 2pi/(length(xtimes)*1e-12)
        xplan = transientquantumplan(xh[1], xtimes, [2xdw, 2xdw]; ports = [1, 2])
        xdplan = transientquantumplan(xh[1], xtimes, [2xdw, 2xdw]; ports = [1, 2], backend = CUDABackend())
        xnh = transientnoise(xh, xplan; frequencies = [2xdw, 3xdw], weights = fill(xdw, 2))
        xnd = transientnoise(xd, xdplan; frequencies = [2xdw, 3xdw], weights = fill(xdw, 2))
        xnf = transientnoise(xd, xdplan; frequencies = [2xdw, 3xdw], weights = fill(xdw, 2), method = :forward)
        @test xnd.covariance ≈ xnh.covariance rtol=1e-8
        @test xnf.covariance ≈ xnh.covariance rtol=1e-8
        @test xnf.commutator ≈ xnh.commutator rtol=1e-8
        xth = transientsolve(xbatch[3], (0.0, 100e-12); dt = 1e-12, method = Trapezoidal(), record = :phases)
        xtd = transientsolve(xbatch[3], (0.0, 100e-12); dt = 1e-12, method = Trapezoidal(), record = :phases, backend = CUDABackend())
        @test Array(xtd.voltage) ≈ xth.voltage rtol=1e-8 atol=1e-14
        @test transientnoise(xtd, xdplan; frequencies = [2xdw, 3xdw], weights = fill(xdw, 2), method = :forward).covariance ≈
            transientnoise(xth, xplan; frequencies = [2xdw, 3xdw], weights = fill(xdw, 2)).covariance rtol=1e-8
        xsys = JosephsonCircuits.transientsystem(xbatch[1], 1e-12, GaussLegendre(), CUDABackend(),
            JosephsonCircuits.steppingfactorization(nothing, CUDABackend()))
        xst = JosephsonCircuits.gaussstepper(xsys, [xbatch[1]], 1e-9, 1e-10, 15, JosephsonCircuits.gaussbatchfactor(xsys, 1))
        JosephsonCircuits.setstate!(xst, zeros(length(xbatch[1]), 1), zeros(length(xbatch[1]), 1), nothing)
        for k in 1:40
            JosephsonCircuits.advance!(xst, (k - 1)*1e-12, k*1e-12, k)
        end
        @test xst.pw.factorizations[1] <= 4
        # a constant scattering block on the device against the host: the
        # attenuator and the circulator, whose transposed solves are a
        # second cuDSS factorization
        Sc = [0.0 0.0 1.0; 1.0 0.0 0.0; 0.0 1.0 0.0]
        kc = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)),
            (:att, 1, 4, ScatteringParameters([0.0 0.7; 0.7 0.0]; zref = [50.0, 75.0])),
            (:circ, 4, 2, 3, ScatteringParameters(Sc; zref = 75.0)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(1e-12)), (:p2, 2, 0, Port(2; Z0 = 75.0)),
            (:c3, 3, 0, Capacitor(0.2e-12)), (:p3, 3, 0, Port(3; Z0 = 75.0))])
        kp = transientproblem(kc; sources = [TransientSource(1, t -> t <= 0 ? 0.0 : 1e-6*sinpi(t/1e-9)^2*sinpi(2*4e9*t)),
            TransientSource(2, t -> 0.0), TransientSource(3, t -> 0.0)])
        kh = transientsolve(kp, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre())
        kd = transientsolve(kp, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend())
        @test Array(kd.outgoing) ≈ kh.outgoing rtol=1e-8 atol=1e-14
        bcurrents = [1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k) for k in 1:3, t in kh.times]
        bweights = [cospi(2*1.7e9*t + k) for k in 1:3, t in kh.times]
        @test Array(transienttangent(kd, bcurrents).outgoing) ≈ transienttangent(kh, bcurrents).outgoing rtol=1e-8
        @test Array(transientadjoint(kd, bweights).currents) ≈ transientadjoint(kh, bweights).currents rtol=1e-8
        # a transmission line on the device against the host: the history
        # read on the host, the forcing and the endpoint reading on the
        # device
        lc = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.2e-12)), (:line, 1, 2, TransmissionLine(60.0, 0.09)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.4e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        lp = transientproblem(lc; sources = [TransientSource(1, t -> t <= 0 ? 0.0 : 0.2e-6*sinpi(2*3e9*t))])
        lh = transientsolve(lp, (0.0, 2e-9); dt = 2e-12, record = :states, rtol = 1e-12, method = GaussLegendre())
        ld = transientsolve(lp, (0.0, 2e-9); dt = 2e-12, record = :states, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend())
        @test Array(ld.outgoing) ≈ lh.outgoing rtol=1e-8 atol=1e-14
        @test Array(ld.linewaves) ≈ lh.linewaves rtol=1e-8 atol=1e-14
        # the noise of a line before a pumped junction on the device
        lnh = transientnoise(transientsolve(lp, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre()), plan;
            frequencies = [2dw, 3dw], weights = fill(dw, 2))
        lnd = transientnoise(transientsolve(lp, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend()), dplan;
            frequencies = [2dw, 3dw], weights = fill(dw, 2))
        @test lnd.covariance ≈ lnh.covariance rtol=1e-8
        @test lnd.commutator ≈ lnh.commutator rtol=1e-8
        # a rational block on the device against the host
        ra = 2*50.0/2e-9
        rblock = RationalScattering(fill(-ra, 1, 1), reshape([1.0, -1.0], 1, 2), reshape(-ra .* [1.0, -1.0], 2, 1), Matrix(1.0I, 2, 2); zref = 50.0)
        rc = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)), (:blk, 1, 2, rblock),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        rp = transientproblem(rc; sources = [TransientSource(1, t -> t <= 0 ? 0.0 : 0.3e-6*sinpi(t/1e-9)^2*sinpi(2*3e9*t))])
        rh = transientsolve(rp, (0.0, 1e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre())
        rd = transientsolve(rp, (0.0, 1e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend())
        @test Array(rd.outgoing) ≈ rh.outgoing rtol=1e-8 atol=1e-14
        rph = transientsolve(rp, (0.0, 0.5e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre(), record = :phases)
        rpd = transientsolve(rp, (0.0, 0.5e-9); dt = 2e-12, rtol = 1e-12, method = GaussLegendre(), record = :phases, backend = CUDABackend())
        rcur = [1e-8*sinpi(2*1.1e9*t)^2*cospi(2*0.7e9*t + k) for k in 1:2, t in rph.times]
        rwts = [cospi(2*1.7e9*t + k) for k in 1:2, t in rph.times]
        @test Array(transienttangent(rpd, rcur).outgoing) ≈ transienttangent(rph, rcur).outgoing rtol=1e-8
        @test Array(transientadjoint(rpd, rwts).currents) ≈ transientadjoint(rph, rwts).currents rtol=1e-8
        rnh = transientnoise(rph, plan; frequencies = [2dw, 3dw], weights = fill(dw, 2))
        rnd = transientnoise(rpd, dplan; frequencies = [2dw, 3dw], weights = fill(dw, 2))
        @test rnd.covariance ≈ rnh.covariance rtol=1e-8
        # a checkpointed record replayed on the device against the host
        ch = transientsolve(problems, (0.0, 0.5e-9); dt = 2e-12, record = :checkpoints, checkpointevery = 50, rtol = 1e-12)
        cd_ = transientsolve(problems, (0.0, 0.5e-9); dt = 2e-12, record = :checkpoints, checkpointevery = 50, rtol = 1e-12, backend = CUDABackend())
        @test cd_.checkpoints.flux isa CuArray
        @test Array(transientadjoint(cd_, weights).currents) ≈ transientadjoint(ch, weights).currents rtol=1e-8 atol=1e-14
        @test Array(transientadjoint(cd_, weights).currents) ≈ transientadjoint(bh, weights).currents rtol=1e-8 atol=1e-14
        big = make.(range(0.05e-6, 0.2e-6; length = 130))
        bigh = transientsolve(big, (0.0, 0.2e-9); dt = 2e-12, rtol = 1e-12)
        bigd = transientsolve(big, (0.0, 0.2e-9); dt = 2e-12, rtol = 1e-12, backend = CUDABackend())
        @test Array(bigd.voltage) ≈ bigh.voltage rtol=1e-8 atol=1e-14
        # a condition refreshes its own factorization on the device, where
        # the uniform batch refactorizes a whole chunk: a junction biased
        # hard enough that its frozen operator drifts beside one that
        # never does, against the host, the tangent and the adjoint included
        pcircuit = Circuit([("P1", "1", "0", Port(1; Z0 = 50.0)), ("Lj1", "1", "0", JosephsonJunction(1e-9)),
            ("C1", "1", "0", Capacitor(0.2e-12))])
        pramp(t) = t <= 0 ? 0.0 : t >= 0.5e-9 ? 1.0 : (1 - cospi(t/0.5e-9))/2
        pbase = transientproblem(pcircuit; sources = [TransientSource(1, t -> 0.0)])
        pair = [transientproblem(pbase; sources = [TransientSource(1, let ib = ib; t -> ib*pramp(t) + 1e-9*sinpi(2*3e9*t); end)])
            for ib in (0.2e-6, 0.6e-6)]
        ph = transientsolve(pair, (0.0, 333*15e-12); dt = 15e-12, record = :phases)
        pdv = transientsolve(pair, (0.0, 333*15e-12); dt = 15e-12, record = :phases, backend = CUDABackend())
        @test pdv.stats.factorizations > 10
        @test Array(pdv.voltage) ≈ ph.voltage rtol=1e-8 atol=1e-14
        pweights = [cospi(2*3.1e9*t) for _ in 1:1, t in ph.times]
        @test Array(transientadjoint(pdv, pweights).currents) ≈ transientadjoint(ph, pweights).currents rtol=1e-8 atol=1e-14
        pcurrents = [1e-9*sinpi(2*2.9e9*t) for _ in 1:1, t in ph.times]
        @test Array(transienttangent(pdv, pcurrents).outgoing) ≈ transienttangent(ph, pcurrents).outgoing rtol=1e-8 atol=1e-14
        # a tangent's currents given on the device into a floating branch
        # a port reads, two sources of one waveform into it and out of it,
        # judged on the host form the steps take, against the host (review
        # of bundle 8, finding 2)
        flc = Circuit([("C1", "1", "2", Capacitor(1e-12)), ("P1", "1", "2", Port(1; Z0 = 50.0)),
            ("I1", "0", "1", CurrentSource(0.0)), ("I2", "2", "0", CurrentSource(0.0))])
        flwave(t) = 1e-6*sinpi(2e9*t)
        flp = transientproblem(flc; sources = [TransientSource(:I1, flwave), TransientSource(:I2, flwave)])
        for method in (GaussLegendre(), Trapezoidal())
            flh = transientsolve(flp, (0.0, 0.2e-9); dt = 1e-12, method, record = :phases)
            fld = transientsolve(flp, (0.0, 0.2e-9); dt = 1e-12, method, record = :phases, backend = CUDABackend())
            fl = 1e-7 .* cospi.(3e9 .* flh.times)
            @test Array(transienttangent(fld, CuArray(permutedims([fl fl])); targets = ["I1", "I2"]).voltage) ≈
                transienttangent(flh, permutedims([fl fl]); targets = ["I1", "I2"]).voltage rtol=1e-8 atol=1e-14
            @test_throws ArgumentError transienttangent(fld, CuArray(permutedims([fl zero(fl)])); targets = ["I1", "I2"])
        end
        # the column maxima the Newton control of a batch reads, by one
        # kernel on the device, against the host's loop and the reductions
        # of CUDA.jl on the same arrays: more rows than a workgroup has
        # items, a row far over its tolerance, and junction currents with
        # rows and without
        JC = JosephsonCircuits
        n, N = 600, 3
        res, tol, delta, X, rhs = (randn(n, N, 2) for _ in 1:5)
        tol .= abs.(tol) .+ 0.1
        res[7, 2, 2] = 50.0
        x = randn(n, N)
        for nj in (5, 0)
            jcur = randn(nj, N, 2)
            host = JC.columnmaxima!(zeros(JC.COLUMNMAXIMA, N), nothing, res, tol, delta, x, X, jcur, rhs, JC.CPU())
            d = CuArray.((res, tol, delta, x, X, jcur, rhs))
            dev = JC.columnmaxima!(zeros(JC.COLUMNMAXIMA, N), CUDA.zeros(Float64, JC.COLUMNMAXIMA, N), d..., CUDABackend())
            @test dev == host
            colmax(a) = vec(Array(maximum(abs, a; dims = (1, 3), init = 0.0)))
            @test dev == permutedims(hcat(vec(Array(mapreduce((r, t) -> abs(r)/t, max, d[1], d[2]; dims = (1, 3), init = 0.0))),
                vec(Array(mapreduce((r, t) -> abs(r) > t ? abs(r) : 0.0, max, d[1], d[2]; dims = (1, 3), init = 0.0))),
                colmax(d[3]), vec(Array(maximum(abs, d[4]; dims = 1, init = 0.0))), colmax(d[6]), colmax(d[5]), colmax(d[7])))
        end
    end

    @testset "the batched cuDSS policy reaches both of its callers" begin
        # `_cudss_sweep` sets the package's pivot and refinement policy
        # before applying the caller's own settings, and both the frequency
        # sweep and the transient batch make their solvers through it. The
        # policy can only be checked by reading it back off the solver,
        # since cuDSS applies it inside the factorization.
        SA = JosephsonCircuits.SparseArrays
        At = SA.sparse(transpose(SA.sparse([1,2,3,1], [1,2,3,3],
            ComplexF64[2,3,4,1], 3, 3)))
        rowptr = CuVector{Int32}(SA.getcolptr(At))
        colind = CuVector{Int32}(SA.rowvals(At))
        nzval = CuMatrix{ComplexF64}(repeat(SA.nonzeros(At), 1, 2))
        X = CUDA.zeros(ComplexF64, 3, 1, 2)
        B = CUDA.ones(ComplexF64, 3, 1, 2)
        S = JosephsonCircuits._cudss_sweep(rowptr, colind, nzval, X, B)
        @test CUDSS.cudss_get(S.solver, "pivot_epsilon") == 1e-8
        @test CUDSS.cudss_get(S.solver, "ir_n_steps") == 2
        S2 = JosephsonCircuits._cudss_sweep(rowptr, colind, nzval, X, B;
            pivot_epsilon = 1e-3, ir_n_steps = 5)
        @test CUDSS.cudss_get(S2.solver, "pivot_epsilon") == 1e-3
        @test CUDSS.cudss_get(S2.solver, "ir_n_steps") == 5
        # the batch solves, of two systems and of one, whose descriptors
        # are bound as the single system's matrices, against the host
        xref = Matrix(transpose(At)) \ ones(ComplexF64, 3)
        JosephsonCircuits._cudss_sweepsolve!(S)
        @test Array(X)[:, 1, 1] ≈ xref && Array(X)[:, 1, 2] ≈ xref
        X1 = CUDA.zeros(ComplexF64, 3, 1, 1)
        S1 = JosephsonCircuits._cudss_sweep(rowptr, colind, nzval[:, 1:1], X1,
            CUDA.ones(ComplexF64, 3, 1, 1))
        JosephsonCircuits._cudss_sweepsolve!(S1)
        @test Array(X1)[:, 1, 1] ≈ xref
        Y1 = CUDA.zeros(ComplexF64, 3, 1, 1)
        JosephsonCircuits._cudss_sweepapply!(S1, Y1, CUDA.ones(ComplexF64, 3, 1, 1))
        CUDA.synchronize()
        @test Array(Y1)[:, 1, 1] ≈ xref

        # the transient batch reaches the same call, on a circuit whose
        # scattering block carries the auxiliary port current rows the
        # policy is there for
        Sm = ComplexF64[0 1; 1 0]
        circuit = Circuit([
            ("p1", "1", "0", Port(1)),
            ("b", "1", "2", ScatteringParameters(Sm; zref = 50.0,
                noise = Lossless())),
            ("c", "2", "0", Capacitor(1e-12)),
            ("jj", "2", "0", JosephsonJunction(1e-9))])
        drive(t) = t <= 0 ? 0.0 : 0.1e-6*sinpi(2*3e9*t)
        prob = transientproblem(circuit;
            sources = [TransientSource(1, drive)])
        function onthedevice(; kw...)
            reuse = JosephsonCircuits.TransientReuse()
            sol = transientsolve(prob, (0.0, 0.2e-9); dt = 2e-12,
                method = GaussLegendre(), rtol = 1e-12,
                backend = CUDABackend(), reuse = reuse, kw...)
            return reuse.factor[1].factors[1].solver, sol
        end
        # the step takes the policy's pivot epsilon, which it does not set,
        # and keeps its own refinement, which it does: the forwarding
        # working in both directions at once
        tsolver, tdev = onthedevice()
        @test CUDSS.cudss_get(tsolver, "pivot_epsilon") == 1e-8
        @test CUDSS.cudss_get(tsolver, "ir_n_steps") == 0
        tsolver2, tdev2 = onthedevice(; factorization =
            CUDSSFactorization(pivot_epsilon = 1e-3, ir_n_steps = 5))
        @test CUDSS.cudss_get(tsolver2, "pivot_epsilon") == 1e-3
        @test CUDSS.cudss_get(tsolver2, "ir_n_steps") == 5
        # and the block's rows come back on the host's answer either way
        thost = transientsolve(prob, (0.0, 0.2e-9); dt = 2e-12,
            method = GaussLegendre(), rtol = 1e-12)
        @test Array(tdev.outgoing) ≈ thost.outgoing rtol=1e-8 atol=1e-14
        @test Array(tdev2.outgoing) ≈ thost.outgoing rtol=1e-8 atol=1e-14
    end

    @testset "cuDSS factors of a complex host matrix" begin
        # the Gauss-Legendre step matrix of a host transient is complex:
        # cuDSS factorizes it from the host, and solves the columns of host
        # matrices, which cross to the device through a host vector
        A = JosephsonCircuits.SparseArrays.sparse(ComplexF64[4 1im 0; 1 4 1; 0 -1im 4])
        F = JosephsonCircuits.factorize(CUDSSFactorization(), A)
        B = ComplexF64[1 2im; 3 4; 5im 6]
        X = zeros(ComplexF64, 3, 2)
        JosephsonCircuits.matrixsolve!(view(X, :, 1:2), F, view(B, :, 1:2))
        @test X ≈ Matrix(A) \ B
    end

    @testset "the linearized sweep through cuDSS with a scattering block" begin
        # The batched cuDSS path on a linearized system carrying a
        # scattering block's auxiliary port current rows: both directions,
        # more frequencies than one batch holds so the last batch is short,
        # the outputs against the host, the solutions against the equations
        # assembled at each frequency, and the factorization's options read
        # back off the solver each direction was made with.
        JC = JosephsonCircuits
        Z0 = 50.0
        blockcircuit = Circuit([
            (:p1, 1, 0, Port(1; Z0 = Z0)),
            (:b, 1, 2, ScatteringParameters(ComplexF64[0 1; 1 0]; zref = Z0,
                noise = Lossless())),
            (:cc, 2, 3, Capacitor(100e-15)),
            (:jj, 3, 0, JosephsonJunction(1000e-12)),
            (:cj, 3, 0, Capacitor(1000e-15))])
        wp = (2*pi*5e9,)
        src = [(mode = (1,), port = 1, current = 0.8e-6)]
        ws = 2*pi*collect(range(4.0e9, 6.0e9,
            length = JC.uniformbatchlimit(1) + 2))
        kw = (; returnSnoise = true, keyedarrays = false)
        rh = hbsolve(ws, wp, src, (8,), (8,), blockcircuit; kw...)
        rd = hbsolve(ws, wp, src, (8,), (8,), blockcircuit;
            backend = CUDABackend(), factorization = CUDSSFactorization(),
            kw...)
        @test rd.nonlinear.solverinfo.converged
        for name in (:S, :Snoise, :QE, :CM)
            @test agree(getfield(rh.linearized, name),
                getfield(rd.linearized, name); rtol = 1e-8)
        end
        # the through block is lossless and emits nothing, so its noise
        # outputs are empty; an attenuator at the default passive model
        # emits, and its noise scattering parameters and added noise
        # covariance are compared on something
        @test isempty(rh.linearized.Snoise)
        # the port warm, so that its noise enters the quantum efficiency,
        # the occupations and the output covariance, which the device
        # sweep forms from the covariance its reduction hands back
        lossycircuit = Circuit([
            (:p1, 1, 0, Port(1; Z0 = Z0,
                termination = MatchedTermination(temperature = 0.05))),
            (:b, 1, 2, ScatteringParameters(ComplexF64[0 0.5; 0.5 0];
                zref = Z0)),
            (:cc, 2, 3, Capacitor(100e-15)),
            (:jj, 3, 0, JosephsonJunction(1000e-12)),
            (:cj, 3, 0, Capacitor(1000e-15))])
        kwn = (; returnSnoise = true, returnCnoise = true, returnVout = true,
            keyedarrays = false)
        lh = hbsolve(ws, wp, src, (8,), (8,), lossycircuit; kwn...)
        ld = hbsolve(ws, wp, src, (8,), (8,), lossycircuit;
            backend = CUDABackend(), factorization = CUDSSFactorization(),
            kwn...)
        @test ld.nonlinear.solverinfo.converged
        @test !isempty(lh.linearized.Snoise) && norm(lh.linearized.Snoise) > 0
        @test norm(lh.linearized.Cnoise) > 0
        for name in (:S, :Snoise, :Cnoise, :Vout, :QE, :CM, :nbar)
            @test agree(getfield(lh.linearized, name),
                getfield(ld.linearized, name); rtol = 1e-8)
        end
        # the output covariance without the added one returned
        lv = hbsolve(ws, wp, src, (8,), (8,), lossycircuit;
            backend = CUDABackend(), factorization = CUDSSFactorization(),
            keyedarrays = false, returnVout = true)
        @test agree(lh.linearized.Vout, lv.linearized.Vout; rtol = 1e-8)

        # the sweep itself, batch by batch
        nl = hbnlsolve(wp, (8,), src, blockcircuit; keyedarrays = false)
        psc = JC.compile(blockcircuit)
        sf = JC.removeconjfreqs(JC.truncfreqs(JC.calcfreqsrdft((8,));
            dc = true, odd = true, even = false, maxintermodorder = Inf))
        d = JC.hblinsolve(ws, psc, Dict{Any,Any}(), sf; nonlinear = nl,
            debuglsys = true)
        lsys = d.lsys
        b = Matrix{ComplexF64}(d.bnm)
        spec = (full = true, rows = Int[])
        ds = JC.devicesolutions(lsys, b, ws, CUDABackend(), spec, spec;
            factorization = CUDSSFactorization(pivot_epsilon = 1e-6,
                ir_n_steps = 3))
        @test isnothing(ds.blocks)
        @test ds.nb < length(ws)
        A = copy(lsys.Asparse)
        X = similar(b)
        for lo in 1:ds.nb:length(ws)
            JC.solvebatch!(ds, lo)
            for i in lo:min(lo + ds.nb - 1, length(ws))
                JC.assemblesystemmatrix!(A, lsys, ws[i] .+ d.wpumpmodes)
                At = JC.SparseArrays.sparse(transpose(A))
                # the residual of each direction against the assembled
                # equations, and the solution against the host's
                JC.forwardsolution!(X, ds, i)
                @test norm(A*X - b) <=
                    1e-10*(opnorm(A, 1)*norm(X) + norm(b))
                @test X ≈ A \ b rtol = 1e-8
                JC.adjointsolution!(X, ds, i)
                @test norm(At*X - b) <=
                    1e-10*(opnorm(At, 1)*norm(X) + norm(b))
                @test X ≈ At \ b rtol = 1e-8
            end
        end
        for slot in 1:2
            solver = ds.sweeps[slot].solver
            @test CUDSS.cudss_get(solver, "pivot_epsilon") == 1e-6
            @test CUDSS.cudss_get(solver, "ir_n_steps") == 3
        end

        # the batch is held to the memory budget at cuDSS's estimate of a
        # system, which holds at least the factors cuDSS reports for the
        # pattern; a budget below one system solves one at a time
        bytes = JC.cudsssystembytes(lsys.Asparse, ComplexF64, size(b, 2),
            CUDABackend())
        probe = CudssSolver(CUDA.CUSPARSE.CuSparseMatrixCSR(A), "G", 'F')
        cudss("analysis", probe, CUDA.zeros(ComplexF64, size(A, 1)),
            CUDA.zeros(ComplexF64, size(A, 1)))
        @test bytes >= CUDSS.cudss_get(probe, "lu_nnz")*sizeof(ComplexF64)
        single = JC.devicesolutions(lsys, b, ws, CUDABackend(), spec, spec;
            factorization = CUDSSFactorization(), budget = 1)
        @test single.nb == 1
        JC.solvebatch!(single, 2)
        JC.assemblesystemmatrix!(A, lsys, ws[2] .+ d.wpumpmodes)
        JC.forwardsolution!(X, single, 2)
        @test X ≈ A \ b rtol = 1e-8
    end

    @testset "a sweep a frequency dependent value keeps on the host" begin
        # the stored values change with the frequency, so the sweep runs on
        # host threads with a host factorization whatever the backend, and
        # a device factorization asked for is refused
        lin = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:Z, 1, 2, Resistor(FrequencyDependent(w -> 5.0 + im*w*1e-9))),
            (:C, 2, 0, Capacitor(1e-12))])
        wsl = 2*pi*(4.5:0.1:5.0)*1e9
        la = hblinsolve(wsl, lin; keyedarrays = false)
        lb = hblinsolve(wsl, lin; keyedarrays = false, backend = CUDABackend())
        @test agree(la.S, lb.S; rtol = 1e-12)
        @test_throws ArgumentError hblinsolve(wsl, lin; backend = CUDABackend(),
            factorization = CUDSSFactorization())
    end

    @testset "tabulated blocks beyond their band" begin
        # a lossy block tabulated over 4-6 GHz in front of a pumped JPA:
        # the sidebands leave the band, where the block is zero, holds its
        # end values or continues linearly, its stamps and its noise
        # channels, on the device as on the host
        fs = 2*pi*collect(range(4e9, 6e9; length = 21))
        Stab = cat([ComplexF64[0.1 0.6-0.001k; 0.6-0.001k 0.1]
            for k in eachindex(fs)]...; dims = 3)
        cz(x) = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:x, 1, 2, x),
            (:C1, 2, 3, Capacitor(100e-15)), (:Lj1, 3, 0, JosephsonJunction(1000e-12)),
            (:C2, 3, 0, Capacitor(1000e-15))])
        wsz = 2*pi*(4.5:0.1:5.0)*1e9
        wpz = (2*pi*4.75001e9,)
        srcz = [(mode = (1,), port = 1, current = 0.00565e-6/0.6)]
        run(x; kw...) = hbsolve(wsz, wpz, srcz, (8,), (8,), cz(x);
            keyedarrays = false, kw...).linearized
        for extrapolation in (:zero, :constant, :linear)
            tab = ScatteringParameters((fs, Stab); zref = 50.0,
                extrapolation = extrapolation)
            za = run(tab)
            zb = run(tab; backend = CUDABackend())
            for name in (:S, :QE, :CM)
                @test agree(getfield(za, name), getfield(zb, name); rtol = 1e-10)
            end
        end
        # a table whose data must not be extrapolated is refused where a
        # sweep reads it beyond its band: a lossless one, which its
        # construction holds to its declaration everywhere, reaches the
        # device's own check of the range, and is solved within its band
        lossless = ScatteringParameters((fs, repeat(ComplexF64[0 1; 1 0], 1, 1,
            length(fs))); zref = 50.0, noise = Lossless())
        through = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:x, 1, 2, lossless),
            (:P2, 2, 0, Port(2; Z0 = 50.0))])
        @test_throws ArgumentError hblinsolve(2*pi*[3.0e9, 5.0e9], through;
            backend = CUDABackend())
        @test agree(hblinsolve(2*pi*[4.5e9, 5.0e9], through; keyedarrays = false).S,
            hblinsolve(2*pi*[4.5e9, 5.0e9], through; keyedarrays = false,
                backend = CUDABackend()).S; rtol = 1e-12)
        # a block whose ports have unequal reference impedances, whose
        # coefficients take the impedance of the port each column is of
        unequal = ScatteringParameters(ComplexF64[0.1 0.5im; 0.6 0.3];
            zref = [25.0, 100.0])
        ua = run(unequal)
        ub = run(unequal; backend = CUDABackend())
        for name in (:S, :QE, :CM)
            @test agree(getfield(ua, name), getfield(ub, name); rtol = 1e-10)
        end
    end

    # the period map recorded and carried on the device: the poles, their
    # profiles and their rate errors against the host's
    @testset "the period map on the device" begin
        junction = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)), (:j, 1, 0, JosephsonJunction(1.0)),
            (:c, 1, 0, Capacitor(1.0))])
        pump = hbnlsolve((0.8,), (12,), [(mode = (1,), port = 1, current = 0.04*JosephsonCircuits.phi0)], junction;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        mapped(backend) = hbstability(junction; nonlinear = pump, Nmodulationharmonics = (2,),
            method = Monodromy(; nev = :all, backend))
        host, device = mapped(JosephsonCircuits.CPU()), mapped(CUDABackend())
        # a mode's profile is its vector's, whose phase each eigensolver
        # sets its own way: compared up to it, mode by mode
        alike(a, b; rtol) = all(axes(a, 3)) do j
            x, y = vec(a[:, :, j]), vec(b[:, :, j])
            phase = dot(y, x)
            agree(x, y .* (phase/abs(phase)); rtol)
        end
        @test length(device.poles) == length(host.poles) == 2
        @test agree(host.poles, device.poles; rtol = 1e-12)
        @test alike(host.nodevoltage, device.nodevoltage; rtol = 1e-10)
        @test agree(host.rateerrors, device.rateerrors; rtol = 1e-3)
        # the JPA's two modes, real multipliers whose content is the same at
        # the pump frequency and at its mirror: each placed alike
        jpa = Circuit([(:p, 1, 0, Port(1; Z0 = 50.0)), (:cc, 1, 2, Capacitor(100e-15)),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:cj, 2, 0, Capacitor(1e-12))])
        jpump = hbnlsolve((2pi*4.75001e9,), (16,), [(mode = (1,), port = 1, current = 0.00565e-6)], jpa)
        jmapped(backend) = hbstability(jpa; nonlinear = jpump, Nmodulationharmonics = (2,),
            method = Monodromy(; nev = :all, backend))
        jhost, jdevice = jmapped(JosephsonCircuits.CPU()), jmapped(CUDABackend())
        @test all(s -> imag(s) >= 0, jhost.poles)
        @test agree(jhost.poles, jdevice.poles; rtol = 1e-12)
        @test alike(jhost.nodevoltage, jdevice.nodevoltage; rtol = 1e-10)
        # a junction behind a transmission line, whose history the map
        # carries as each line port's run of coordinates
        behind = Circuit([(:p, 1, 0, Port(1; Z0 = 20.0)), (:line, 1, 2, TransmissionLine(10.0, 1.0; vp = 1.0)),
            (:j, 2, 0, JosephsonJunction(1.0)), (:c, 2, 0, Capacitor(1.0))])
        bpump = hbnlsolve((0.8,), (12,), [(mode = (1,), port = 1, current = 0.04*JosephsonCircuits.phi0)], behind;
            method = Newton(), atol = 1e-12, keyedarrays = false)
        bmapped(backend) = hbstability(behind; nonlinear = bpump, Nmodulationharmonics = (4,),
            method = Monodromy(; nev = 2, steps = 64, backend))
        bhost, bdevice = bmapped(JosephsonCircuits.CPU()), bmapped(CUDABackend())
        @test length(bdevice.poles) == length(bhost.poles) == 2
        @test agree(bhost.poles, bdevice.poles; rtol = 1e-10)
        # an isolated multiplier beside a defective block, whose right
        # vectors are singular together, and the same turned into every
        # coordinate by a reflection, whose inverse holds no left vector to
        # working accuracy: the device's eigensolve gives the host's
        jordan(m) = (J = diagm(vcat(1.01, fill(0.2, m - 1))); foreach(k -> J[k, k + 1] = 1.0, 2:m - 1); J)
        w = collect(1.0:8.0)
        reflect = I - 2*w*w'/(w'*w)
        for M in (jordan(24), reflect*jordan(8)*reflect)
            h = JosephsonCircuits.mapeigen!(copy(M), 1, JosephsonCircuits.CPU())
            d = JosephsonCircuits.mapeigen!(copy(M), 1, CUDABackend())
            @test h.values[h.selected] ≈ d.values[d.selected] ≈ [1.01] && d.conditions ≈ h.conditions
            hv, dv = vec(h.left), vec(d.left)
            @test abs(dot(hv, dv)) ≈ norm(hv)*norm(dv) rtol = 1e-12
        end
        # the device's residual check at its ends: a left vector's residual
        # decides between the device's vectors and the host's, every vector
        # passing a tolerance of infinity, and none of this map's, whose
        # residuals are nonzero, a tolerance of zero
        B = [1/(i + j) + (i == j)*i for i in 1:8, j in 1:8]
        ilo, ihi, _ = LAPACK.gebal!('S', B)
        for (tolerance, kept) in ((Inf, true), (0.0, false))
            multipliers, vectors = JosephsonCircuits.mapspectrum!(copy(B), ilo, ihi, CUDABackend(); tolerance)
            @test isnothing(vectors([argmax(abs.(multipliers))])) == !kept
        end
    end

    # the I/Q and quantum measurements on the device
    testtransientiq(CUDABackend())
    testtransientquantum(CUDABackend())
end
