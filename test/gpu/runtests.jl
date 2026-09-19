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

    @testset "a polynomial current-phase relation" begin
        # the polynomial relations are evaluated by Horner in whole array
        # broadcasts along the junction axis, and the sinusoidal columns of
        # a circuit which mixes the two kinds are written over afterwards:
        # both paths run on the device as they do on the host
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
            backend = CUDABackend(), atol = 1e-10)
        hbsolve!(cache, p0)
        @test cache.converged
        pm = cache.reuse.sys.phimatrix
        for (Lj, Cg) in ((105e-12, 40e-15), (95e-12, 44e-15), (110e-12, 38e-15))
            s = hbsolve!(cache, (Lj = Lj, Cg = Cg))
            @test cache.converged
            @test cache.reuse.sys.phimatrix === pm
            fresh = hbnlsolve((w1,), (8,), src1, make(; Lj = Lj, Cg = Cg);
                backend = CUDABackend(), atol = 1e-10, keyedarrays = false)
            host = hbnlsolve((w1,), (8,), src1, make(; Lj = Lj, Cg = Cg);
                atol = 1e-10, keyedarrays = false)
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
        # the block factorization on the device: both directions from one
        # factorization per frequency, every output, exact and refined
        # single precision factors
        kw = (; dc = true, threewavemixing = true, fourwavemixing = true,
            returnSnoise = true, keyedarrays = false)
        wl = 2*pi*collect(range(4.41e9, 5.57e9, length = 4))
        rc = hbsolve(wl, (w1,w2), src2, (2,2), (8,4), circuit, defs; kw...)
        for f in (BlockFactorization(), BlockFactorization(precision = Float32))
            rd = hbsolve(wl, (w1,w2), src2, (2,2), (8,4), circuit, defs;
                backend = CUDABackend(), factorization = f, kw...)
            for name in (:S, :Snoise, :QE, :CM)
                @test agree(getfield(rc.linearized, name),
                    getfield(rd.linearized, name))
            end
        end
        # the automatic choice on the device: the block factorization for
        # two tones, agreeing with the host
        auto = hbsolve(wl, (w1,w2), src2, (2,2), (8,4), circuit, defs;
            backend = CUDABackend(), kw...)
        @test agree(rc.linearized.S, auto.linearized.S)
        @test agree(rc.linearized.QE, auto.linearized.QE)
        Asp = JosephsonCircuits.SparseArrays.sparse(
            ComplexF64[1 1 0 0; 1 1 1 0; 0 1 1 1; 0 0 1 1])
        @test JosephsonCircuits.linearizedfactorization(Asp, 2, 2,
            CUDABackend()) isa BlockFactorization
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
        @test JosephsonCircuits.freememory(CUDABackend()) > 0
    end

    @testset "single precision cuDSS factors as the preconditioner" begin
        # A preconditioner's factors only have to make the Krylov solve
        # converge, not carry the accuracy of the answer, so they may be
        # held in single precision while the iteration stays double;
        # `applypreconditioner!` converts the residual down and the
        # correction back. cuDSS is the only factorization here which
        # honours a precision: KLU and UMFPACK are compiled for double and
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
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, d.Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; Amatrixmodes = d.Amatrixmodes,
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
        @test transientdemodulate(device, 2, 3e9) ≈ transientdemodulate(host, 2, 3e9) rtol=1e-8
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
        df = 1/(length(qh.times)*2e-12)
        plan = transientquantumplan(qh, qh.times, [2df, 2df]; ports = [1, 2])
        dplan = transientquantumplan(qh, qh.times, [2df, 2df]; ports = [1, 2], backend = CUDABackend())
        nh = transientnoise(qh, plan; frequencies = [2df, 3df], weights = fill(df, 2))
        nd = transientnoise(qd, dplan; frequencies = [2df, 3df], weights = fill(df, 2))
        @test nd.covariance ≈ nh.covariance rtol=1e-8
        @test nd.commutator ≈ nh.commutator rtol=1e-8
        # the forward method contracts the device responses with the device
        # measurement, and agrees with the adjoint on the loaded trajectory
        fd = transientnoise(qd, dplan; frequencies = [2df, 3df], weights = fill(df, 2), method = :forward)
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
        gnh = transientnoise(gqh, plan; frequencies = [2df, 3df], weights = fill(df, 2))
        gnd = transientnoise(gqd, dplan; frequencies = [2df, 3df], weights = fill(df, 2))
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
            frequencies = [2df, 3df], weights = fill(df, 2))
        nbd = transientnoise(transientsolve(qb, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, backend = CUDABackend()), dplan;
            frequencies = [2df, 3df], weights = fill(df, 2))
        @test nbd.covariance ≈ nbh.covariance rtol=1e-8
        @test nbd.commutator ≈ nbh.commutator rtol=1e-8
        nfd = transientnoise(transientsolve(qb, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, backend = CUDABackend()), dplan;
            frequencies = [2df, 3df], weights = fill(df, 2), method = :forward)
        @test nfd.covariance ≈ nbh.covariance rtol=1e-8
        # a warm attenuator's channels on the device against the host
        ac = Circuit([(:p1, 1, 0, Port(1; Z0 = 50.0)), (:c1, 1, 0, Capacitor(0.3e-12)),
            (:att, 1, 2, ScatteringParameters([0.0 0.6; 0.6 0.0]; zref = 50.0, noise = ThermalEquilibrium(0.3))),
            (:jj, 2, 0, JosephsonJunction(1e-9)), (:c2, 2, 0, Capacitor(0.5e-12)), (:p2, 2, 0, Port(2; Z0 = 50.0))])
        ap = transientproblem(ac; sources = [TransientSource(1, t -> t <= 0 ? 0.0 : 0.05e-6*sinpi(2*3e9*t))])
        aplan = transientquantumplan(qh, qh.times, [2df, 3df]; ports = [1, 2])
        adplan = transientquantumplan(qh, qh.times, [2df, 3df]; ports = [1, 2], backend = CUDABackend())
        anh = transientnoise(transientsolve(ap, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre()), aplan;
            frequencies = [2df, 3df], weights = fill(df, 2))
        and = transientnoise(transientsolve(ap, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend()), adplan;
            frequencies = [2df, 3df], weights = fill(df, 2))
        @test and.covariance ≈ anh.covariance rtol=1e-8
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
            frequencies = [2df, 3df], weights = fill(df, 2))
        lnd = transientnoise(transientsolve(lp, (0.0, 0.5e-9); dt = 2e-12, record = :phases, rtol = 1e-12, method = GaussLegendre(), backend = CUDABackend()), dplan;
            frequencies = [2df, 3df], weights = fill(df, 2))
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
        rnh = transientnoise(rph, plan; frequencies = [2df, 3df], weights = fill(df, 2))
        rnd = transientnoise(rpd, dplan; frequencies = [2df, 3df], weights = fill(df, 2))
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
        lossycircuit = Circuit([
            (:p1, 1, 0, Port(1; Z0 = Z0)),
            (:b, 1, 2, ScatteringParameters(ComplexF64[0 0.5; 0.5 0];
                zref = Z0)),
            (:cc, 2, 3, Capacitor(100e-15)),
            (:jj, 3, 0, JosephsonJunction(1000e-12)),
            (:cj, 3, 0, Capacitor(1000e-15))])
        kwn = (; returnSnoise = true, returnCnoise = true,
            keyedarrays = false)
        lh = hbsolve(ws, wp, src, (8,), (8,), lossycircuit; kwn...)
        ld = hbsolve(ws, wp, src, (8,), (8,), lossycircuit;
            backend = CUDABackend(), factorization = CUDSSFactorization(),
            kwn...)
        @test ld.nonlinear.solverinfo.converged
        @test !isempty(lh.linearized.Snoise) && norm(lh.linearized.Snoise) > 0
        @test norm(lh.linearized.Cnoise) > 0
        for name in (:S, :Snoise, :Cnoise, :QE, :CM)
            @test agree(getfield(lh.linearized, name),
                getfield(ld.linearized, name); rtol = 1e-8)
        end

        # the sweep itself, batch by batch
        nl = hbnlsolve(wp, (8,), src, blockcircuit; keyedarrays = false)
        psc = JC.compile(blockcircuit)
        cg = JC.calccircuitgraph(psc)
        sf = JC.removeconjfreqs(JC.truncfreqs(JC.calcfreqsrdft((8,));
            dc = true, odd = true, even = false, maxintermodorder = Inf))
        d = JC.hblinsolve(ws, psc, cg, Dict{Any,Any}(), sf; nonlinear = nl,
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
    end

    # the I/Q and quantum measurements on the device
    testtransientiq(CUDABackend())
    testtransientquantum(CUDABackend())
end
