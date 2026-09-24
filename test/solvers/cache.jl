using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test
isdefined(Main, :recovery_solver) || include("recoveryfixture.jl")

# The solve cache: what it reuses between solves, and the exactness of
# the refill of the linear term it keeps.
@testset verbose=true "hbcache" begin

    wp = (2*pi*4.75001*1e9,)
    src = [(mode=(1,), port=1, current=0.00565e-6)]
    ws = 2*pi*(4.5:0.1:5.0)*1e9
    # a JPA whose junction and coupling capacitance are parameters
    JosephsonCircuits.@params Lj Cc
    circuit = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(Cc)),
        (:Lj1, 2, 0, JosephsonJunction(Lj)), (:C2, 2, 0, Capacitor(1000e-15))])

    @testset "hbcache on a typed circuit" begin
        # a point names the parameters it moves and the rest keep their
        # definition, each solve equals the fresh solve at those values,
        # and the compiled circuit is taken as the circuit is
        defs = Dict(Lj => 1000e-12, Cc => 100e-15)
        cache = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12)
        nl = hbsolve!(cache, (Lj = 1000e-12,); warmstart = false)
        @test cache.converged
        ref = hbnlsolve(wp, (8,), src, circuit, defs; atol = 1e-12, keyedarrays = false)
        @test isapprox(vec(collect(nl.nodeflux)), vec(collect(ref.nodeflux)); rtol = 1e-8)
        moved = hbsolve!(cache, (Lj = 900e-12, Cc = 120e-15))
        @test cache.converged && cache.nsolves == 2
        refmoved = hbnlsolve(wp, (8,), src, circuit, Dict(Lj => 900e-12, Cc => 120e-15); atol = 1e-12, keyedarrays = false)
        @test isapprox(vec(collect(moved.nodeflux)), vec(collect(refmoved.nodeflux)); rtol = 1e-8)
        # the definitions may be keyed by the parameters' symbols, and the
        # compiled circuit given instead of the circuit
        symbolic = hbcache(wp, (8,), src, compile(circuit), Dict(:Lj => 1000e-12, :Cc => 100e-15); atol = 1e-12)
        @test isapprox(vec(collect(hbsolve!(symbolic, (Lj = 1000e-12,)).nodeflux)), vec(collect(ref.nodeflux)); rtol = 1e-8)
        # a parameter the definitions lack is refused when the cache is
        # built, and a value the solve refuses, an infinite junction
        # inductance, when it is solved
        @test_throws ArgumentError hbcache(wp, (8,), src, circuit, Dict(Lj => 1000e-12); atol = 1e-12)
        @test_throws ArgumentError hbsolve!(cache, (Lj = Inf,))
        # the cache evaluates the values at each point, so a frequency
        # dependent one, which has no value until the mode frequency is
        # known, is refused when the cache is built
        fd = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 2, Capacitor(FrequencyDependent(w -> 100e-15))),
            (:Lj1, 2, 0, JosephsonJunction(Lj)),
            (:C2, 2, 0, Capacitor(1000e-15))])
        @test_throws ArgumentError hbcache(wp, (8,), src, fd,
            Dict(Lj => 1000e-12); atol = 1e-12)
    end

    @testset "hbcache reuses the compiled circuit and warm starts" begin
        defs = Dict(:Lj => 1000.0e-12, :Cc => 100.0e-15)
        # the definitions at a point
        at(; p...) = merge(defs, Dict(pairs(p)))
        p = (Lj = 1000.0e-12, Cc = 100.0e-15)
        cache = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12)
        nl = hbsolve!(cache, p; warmstart = false)
        @test cache.converged
        @test cache.nsolves == 1
        ref = hbnlsolve(wp, (8,), src, circuit, at(; p...); atol = 1e-12,
            keyedarrays = false)
        @test isapprox(vec(collect(nl.nodeflux)),
            vec(collect(ref.nodeflux)); rtol = 1e-8)

        # a warm start from the stored point lands in the same place, in
        # fewer iterations
        coldits = sum(st.iterations for st in nl.solverinfo.stages)
        matrixwork = cache.matrixworkspace
        nw = hbsolve!(cache, p)
        @test cache.matrixworkspace === matrixwork
        @test cache.converged
        @test sum(st.iterations for st in nw.solverinfo.stages) <= coldits
        @test isapprox(vec(collect(nw.nodeflux)),
            vec(collect(nl.nodeflux)); rtol = 1e-8)

        # a nearby parameter re-solves against the reference
        p2 = (Lj = 1050.0e-12, Cc = 100.0e-15)
        n2 = hbsolve!(cache, p2)
        @test cache.converged
        ref2 = hbnlsolve(wp, (8,), src, circuit, at(; p2...); atol = 1e-12,
            keyedarrays = false)
        @test isapprox(vec(collect(n2.nodeflux)),
            vec(collect(ref2.nodeflux)); rtol = 1e-6)
        # the waves are normalized to the port's reference impedance
        @test cache.nm.portimpedances == [50.0]
        @test isapprox(Array(n2.S), Array(ref2.S); rtol = 1e-6)

        # the system, the preconditioner and the Krylov vectors were built
        # once and rebound: the same objects carry every point
        pm = cache.reuse.sys.phimatrix
        pc = cache.reuse.preconditioner
        kv = cache.reuse.krylov[]
        n3 = hbsolve!(cache, (Lj = 975.0e-12, Cc = 110.0e-15))
        @test cache.converged
        @test cache.reuse.sys.phimatrix === pm
        @test cache.reuse.preconditioner === pc
        @test cache.reuse.krylov[] === kv
        ref3 = hbnlsolve(wp, (8,), src, circuit, at(; Lj = 975.0e-12,
            Cc = 110.0e-15); atol = 1e-12, keyedarrays = false)
        @test isapprox(vec(collect(n3.nodeflux)),
            vec(collect(ref3.nodeflux)); rtol = 1e-8)
        @test isapprox(Array(n3.S), Array(ref3.S); rtol = 1e-8)

        # a coupling which is zero at the first point keeps its place in the
        # rebound system, so the next point solves the circuit with it
        zc = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12)
        hbsolve!(zc, (Lj = 1000.0e-12, Cc = 0.0); warmstart = false)
        nz = hbsolve!(zc, p)
        @test zc.converged
        @test isapprox(vec(collect(nz.nodeflux)), vec(collect(nl.nodeflux));
            rtol = 1e-8)

        # the assembled Jacobian of a rebound system is that of the new
        # values, not of the first point: the operating point of a cached
        # solve matches a cold one, and a direct Newton cache converges
        # like a cold solve
        opc = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12,
            returnoperatingpoint = true)
        o1 = hbsolve!(opc, p; warmstart = false).operatingpoint
        o2 = hbsolve!(opc, p2).operatingpoint
        r2 = hbnlsolve(wp, (8,), src, circuit, at(; p2...); atol = 1e-12,
            keyedarrays = false, returnoperatingpoint = true).operatingpoint
        @test isapprox(o2.jacobian, r2.jacobian; rtol = 1e-8)
        @test norm(o2.jacobian - r2.jacobian) < 1e-8*norm(r2.jacobian)
        # and it owns its system, which the next solve does not rebind: the
        # first point's operating point is still a cold solve's there, down
        # to the residual derivative the sensitivities read from it
        r1 = hbnlsolve(wp, (8,), src, circuit, at(; p...); atol = 1e-12,
            keyedarrays = false, returnoperatingpoint = true).operatingpoint
        @test isapprox(o1.jacobian, r1.jacobian; rtol = 1e-8)
        nm1 = numericmatrices(opc.compiled, at(; p...); Nmodes = opc.Nmodes)
        idx = [JosephsonCircuits.componentindex(opc.compiled, n) for n in ("Lj1", "C2")]
        @test isapprox(
            Matrix(JosephsonCircuits.calcresidualsensitivity(o1, opc.compiled, nm1, idx)),
            Matrix(JosephsonCircuits.calcresidualsensitivity(r1, opc.compiled, nm1, idx));
            rtol = 1e-8)
        nc = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12, method = Newton())
        hbsolve!(nc, p; warmstart = false)
        nn = hbsolve!(nc, p2; warmstart = false)
        rn = hbnlsolve(wp, (8,), src, circuit, at(; p2...); atol = 1e-12,
            keyedarrays = false, method = Newton())
        @test isapprox(vec(collect(nn.nodeflux)),
            vec(collect(rn.nodeflux)); rtol = 1e-10)
        @test sum(st.iterations for st in nn.solverinfo.stages) ==
            sum(st.iterations for st in rn.solverinfo.stages)

        # a reset discards the stored state
        JosephsonCircuits.reset!(cache)
        @test isnothing(cache.x)
        @test !cache.converged

        # the keywords the cache manages, an unsupported method and a
        # keyword the compiled solve does not take are refused at
        # construction, and a solver keyword still reaches the solve
        for bad in ((x0 = zeros(2),), (keyedarrays = true,),
                (reuse = nothing,), (method = Staged(),), (nosuchkeyword = 1,),
                (maxharmonics = (8,),), (returnsystem = true,),
                (debugJacobian = true,))
            @test_throws ArgumentError hbcache(wp, (8,), src, circuit, defs; bad...)
        end
        # what the cache does anyway is accepted
        @test hbcache(wp, (8,), src, circuit, defs; keyedarrays = false) isa JosephsonCircuits.HBCache
        loose = hbcache(wp, (8,), src, circuit, defs; atol = 1e-2, iterations = 2)
        hbsolve!(loose, p; warmstart = false)
        @test loose.nsolves == 1
        @test sum(st.iterations for st in
            hbsolve!(loose, p; warmstart = false).solverinfo.stages) <= 2
    end

    @testset "failure, rejection, reset, and retained results" begin
        defs = Dict(:Lj => 1e-9, :Cc => 100e-15)
        controlled = recovery_solver(failat = [2])
        cache = hbcache(wp, (4,), src, circuit, defs;
            method = controlled.method, warnnotconverged = false)
        first = hbsolve!(cache, (;))
        # the junction vectors too: the cache refills its matrices at every
        # point, failed or refused, and the result keeps its own
        kept(s) = (copy(s.nodeflux), copy(s.S), copy(s.Ljb), copy(s.Ljbm))
        saved = kept(first)
        matrixwork, capacitance = cache.matrixworkspace, cache.nm.Cnm
        @test cache.converged
        failed = hbsolve!(cache, (Lj = 1.02e-9,))
        @test !cache.converged && !failed.solverinfo.converged
        @test controlled.starts[2] ≈ controlled.solutions[1]
        recovered = hbsolve!(cache, (Lj = 0.98e-9,))
        freshsolver = recovery_solver()
        freshcache = hbcache(wp, (4,), src, circuit, defs;
            method = freshsolver.method)
        fresh = hbsolve!(freshcache, (Lj = 0.98e-9,))
        @test cache.converged && cache.nsolves == 3
        # the point after a failure starts from the last converged point,
        # not from the failed state nor cold
        @test controlled.starts[3] ≈ controlled.solutions[1]
        @test controlled.starts[3] != freshsolver.starts[1]
        @test recovered.nodeflux ≈ fresh.nodeflux rtol = 1e-10
        @test cache.matrixworkspace === matrixwork &&
            cache.nm.Cnm === capacitance
        @test kept(first) == saved
        # a point the solve refuses leaves the cache as it was, so the next
        # point still solves
        retained = copy(cache.x)
        @test_throws ArgumentError hbsolve!(cache, (Lj = Inf,))
        @test cache.x == retained && cache.nsolves == 3
        @test hbsolve!(cache, (Lj = 0.98e-9,)).nodeflux ≈ fresh.nodeflux
        JosephsonCircuits.reset!(cache)
        @test isnothing(cache.x) && !cache.converged
        hbsolve!(cache, (Lj = 0.98e-9,))
        @test controlled.starts[end] == freshsolver.starts[1]
        @test kept(first) == saved
    end

    @testset "values which change a group's element type or behavior" begin
        # a capacitance turned lossy, a resistance at infinity and a mutual
        # coupling at one are points like any other, which the direct solve
        # takes: the matrices are assembled anew where a group's element
        # type changes, and refilled otherwise
        lossy = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(:Cc)),
            (:Lj1, 2, 0, JosephsonJunction(:Lj)), (:C2, 2, 0, Capacitor(1000e-15)),
            (:R1, 2, 0, Resistor(:R))])
        defs = Dict(:Lj => 1e-9, :Cc => 100e-15, :R => 1e4)
        cache = hbcache(wp, (8,), src, lossy, defs; atol = 1e-12)
        for p in ((Cc = 100e-15,), (Cc = 100e-15*(1 - im/100),), (R = Inf,),
                (Cc = 110e-15,))
            sol = hbsolve!(cache, p)
            fresh = hbnlsolve(wp, (8,), src, lossy, merge(defs, Dict(pairs(p)));
                atol = 1e-12, keyedarrays = false)
            @test cache.converged
            @test isapprox(vec(collect(sol.nodeflux)),
                vec(collect(fresh.nodeflux)); rtol = 1e-8)
        end
        transformer = Circuit([:p1 => Port(1), :c1 => Capacitor(1e-13),
            :l1 => Inductor(1e-9), :l2 => Inductor(1e-9),
            :jj => JosephsonJunction(1e-9), :cj => Capacitor(1e-12),
            :k => MutualInductor(:K, :l1, :l2)],
            [[(:p1, 1), (:c1, 1)], [(:c1, 2), (:l1, 1)], [(:l2, 1), (:jj, 1), (:cj, 1)],
             [(:p1, 2), (:l1, 2), (:l2, 2), (:jj, 2), (:cj, 2), Ground]])
        tsrc = [(mode = (1,), port = 1, current = 1e-8)]
        cache = hbcache((2*pi*5e9,), (4,), tsrc, transformer, Dict(:K => 0.9);
            atol = 1e-12)
        for K in (0.9, 1.0)
            sol = hbsolve!(cache, (K = K,))
            fresh = hbnlsolve((2*pi*5e9,), (4,), tsrc, transformer, Dict(:K => K);
                atol = 1e-12, keyedarrays = false)
            @test cache.converged
            @test isapprox(vec(collect(sol.nodeflux)),
                vec(collect(fresh.nodeflux)); rtol = 1e-8)
        end
    end

    @testset "a scattering block under a swept port impedance" begin
        # the port impedance moves the solver scale, which the block's rows
        # of the reused linear term carry: every point of a cached sweep is
        # the fresh solve's, and a reuse is refused at another pump, whose
        # mode frequencies it holds
        capS(C, Z0) = w -> fill((1 - im*w*C*Z0)/(1 + im*w*C*Z0), 1, 1)
        shunt = ScatteringParameters(capS(1000e-15, 50.0); nports = 1,
            grounded = true)
        jpa = Circuit([:p1 => Port(1; Z0 = :Rp), :cc => Capacitor(100e-15),
            :jj => JosephsonJunction(:Lj), :c2 => shunt],
            [((:p1, 1), (:cc, 1)), ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:p1, 2), Ground)])
        defs = Dict(:Rp => 50.0, :Lj => 1000e-12)
        cache = hbcache(wp, (8,), src, jpa, defs; atol = 1e-12)
        for p in ((Rp = 50.0,), (Rp = 49.0,), (Rp = 48.0, Lj = 1010e-12),
                (Rp = 50.0,))
            sol = hbsolve!(cache, p)
            fresh = hbnlsolve(wp, (8,), src, jpa, merge(defs, Dict(pairs(p)));
                atol = 1e-12, keyedarrays = false)
            @test cache.converged
            @test isapprox(vec(collect(sol.nodeflux)),
                vec(collect(fresh.nodeflux)); rtol = 1e-8)
        end
        @test_throws ArgumentError hbnlsolve((2*pi*4.6e9,), cache.sources,
            cache.frequencies, cache.indices, cache.compiled, cache.nm;
            reuse = cache.reuse, keyedarrays = false)
        # a start from node fluxes sets the block's port currents from
        # them, as the cache's warm starts do: a restart from a converged
        # point is converged, rather than starting at the block's response
        cold = hbnlsolve(wp, (8,), src, jpa, defs; atol = 1e-12,
            keyedarrays = false)
        restart = hbnlsolve(wp, (8,), src, jpa, defs; atol = 1e-12,
            keyedarrays = false, x0 = cold.nodeflux)
        @test restart.solverinfo.initialresidual < 1e-12
        @test restart.solverinfo.stages[end].iterations == 0
    end

    @testset "a reset drops what the preconditioner grew" begin
        # a preconditioner escalated at one point stays escalated for the
        # next, which is the point of keeping it; a reset starts the next
        # solve from the preconditioner its method asks for, so it solves
        # step for step as a new cache does
        defs = Dict(:Lj => 1000e-12, :Cc => 100e-15)
        method = NewtonKrylov(preconditioner = BlockDiagonal())
        cache = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12, method)
        hbsolve!(cache, (Lj = 1000e-12,))
        # as a hard point would have
        @test JosephsonCircuits.escalatepreconditioner!(cache.reuse.preconditioner)
        JosephsonCircuits.reset!(cache)
        a = hbsolve!(cache, (Lj = 1000e-12,))
        b = hbsolve!(hbcache(wp, (8,), src, circuit, defs; atol = 1e-12, method),
            (Lj = 1000e-12,))
        steps(s) = [k.iterations for k in s.solverinfo.stages[end].krylov]
        @test steps(a) == steps(b)
        @test a.nodeflux == b.nodeflux
    end

    @testset "Newton-Krylov workspace survives an exhausted solve" begin
        defs = Dict(:Lj => 1e-9, :Cc => 100e-15)
        prepared = hbcache(wp, (4,), src, circuit, defs)
        reuse = prepared.reuse
        opts = (; reuse, keyedarrays = false, atol = 1e-12,
            warnnotconverged = false)
        attempt(d; kw...) = hbnlsolve(prepared.w, prepared.sources,
            prepared.frequencies, prepared.indices, prepared.compiled,
            numericmatrices(prepared.compiled, d; Nmodes = prepared.Nmodes);
            opts..., kw...)
        first = attempt(defs)
        @test first.solverinfo.converged
        # a system handed back would be rebound under its holder
        @test_throws ArgumentError attempt(defs; returnsystem = true)
        @test_throws ArgumentError attempt(defs; debugJacobian = true)
        # the reused system is rebound to each new point's junctions, and
        # the first result keeps its own
        saved = (copy(first.nodeflux), copy(first.S), copy(first.Ljb),
            copy(first.Ljbm))
        workspace = (reuse.sys.phimatrix, reuse.preconditioner, reuse.krylov[])
        moved = Dict(:Lj => 1.02e-9, :Cc => 110e-15)
        # a solve given no Newton iterations at a nonzero drive cannot
        # converge, and rebinds the reuse machinery at the moved values on
        # the way
        failed = attempt(moved; iterations = 0)
        @test !failed.solverinfo.converged
        recovered = attempt(moved)
        fresh = hbnlsolve(wp, (4,), src, circuit, moved;
            keyedarrays = false, atol = 1e-12)
        @test recovered.solverinfo.converged
        @test recovered.nodeflux ≈ fresh.nodeflux rtol = 1e-10
        @test reuse.sys.phimatrix === workspace[1]
        @test reuse.preconditioner === workspace[2]
        @test reuse.krylov[] === workspace[3]
        @test (first.nodeflux, first.S, first.Ljb, first.Ljbm) == saved
    end

    @testset "a solution owns its values under either method" begin
        # a cached solve refills the matrices and, under Newton-Krylov,
        # rebinds the system, so a solution sharing either would move with
        # the next point: its branch vectors, and the junction vectors and
        # the Jacobian of its operating point, are its own
        shunted = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:C1, 1, 2, Capacitor(100e-15)),
            (:Lj1, 2, 0, JosephsonJunction(:Lj)),
            (:C2, 2, 0, Capacitor(1000e-15)), (:L1, 1, 0, Inductor(:L))])
        owned(s) = map(copy, (s.Ljb, s.Lb, s.Ljbm, s.operatingpoint.sys.Ljb,
            s.operatingpoint.sys.Ljbm, s.operatingpoint.jacobian))
        for method in (NewtonKrylov(), Newton())
            c = hbcache(wp, (4,), src, shunted, Dict(:Lj => 1e-9, :L => 10e-9);
                method = method, atol = 1e-12, returnoperatingpoint = true)
            first = hbsolve!(c, (;); warmstart = false)
            saved = owned(first)
            second = hbsolve!(c, (Lj = 1.2e-9, L = 12e-9))
            @test c.converged
            @test second.Ljb != first.Ljb && second.Lb != first.Lb
            @test owned(first) == saved
        end
    end

    @testset "the definitions at a point" begin
        # a parameter defined under its parameter object, its symbol and its
        # string moves under every one of them, a name the definitions lack
        # is refused, and the base definitions stay
        JosephsonCircuits.@params La
        d = Dict{Any,Any}(La => 1.0, :La => 1.0, "La" => 1.0, :Cb => 2.0)
        index = JosephsonCircuits.definitionkeys(d)
        at = JosephsonCircuits.definitionsat(d, index, (La = 3.0, Cb = 4.0))
        @test at[La] == at[:La] == at["La"] == 3.0
        @test at[:Cb] == 4.0
        @test d[La] == 1.0 && d[:Cb] == 2.0
        @test_throws ArgumentError JosephsonCircuits.definitionsat(d, index,
            (La = 3.0, Ln = 5.0))
        # the cache refuses a misspelt parameter, and holds its own copy of
        # the definitions, which the caller's dictionary does not move
        defs = Dict{Any,Any}(:Lj => 1000e-12, :Cc => 100e-15)
        cache = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12)
        a = hbsolve!(cache, (Lj = 1000e-12,))
        @test_throws ArgumentError hbsolve!(cache, (Ljj = 900e-12,))
        @test cache.nsolves == 1
        defs[:Cc] = 150e-15
        b = hbsolve!(cache, (Lj = 1000e-12,); warmstart = false)
        @test isapprox(vec(collect(b.nodeflux)), vec(collect(a.nodeflux));
            rtol = 1e-10)
        # a point moving every parameter of a large set
        n = 4096
        names = ntuple(i -> Symbol(:p, i), n)
        big = Dict{Any,Any}(k => 1.0 for k in names)
        moved = JosephsonCircuits.definitionsat(big, JosephsonCircuits.definitionkeys(big),
            NamedTuple{names}(ntuple(i -> 2.0, n)))
        @test all(moved[k] == 2.0 for k in names)
        # and a cache whose every value is its own parameter, moved at once
        Nj = 8
        ladder = Any[(:P1, 1, 0, Port(1; Z0 = 50.0))]
        for i in 1:Nj
            push!(ladder, (Symbol(:Lj, i), i, i + 1, JosephsonJunction(Symbol(:Lj, i))))
            push!(ladder, (Symbol(:Cg, i), i + 1, 0, Capacitor(Symbol(:Cg, i))))
        end
        push!(ladder, (:P2, Nj + 1, 0, Port(2; Z0 = 50.0)))
        base = Dict{Symbol,Float64}()
        for i in 1:Nj
            base[Symbol(:Lj, i)] = 1e-9; base[Symbol(:Cg, i)] = 100e-15
        end
        pt = NamedTuple(k => 1.1v for (k, v) in base)
        cl = hbcache((2pi*6e9,), (4,), [(mode=(1,), port=1, current=1e-7)],
            Circuit(ladder), base)
        sol = hbsolve!(cl, pt)
        ref = hbnlsolve((2pi*6e9,), (4,), [(mode=(1,), port=1, current=1e-7)],
            Circuit(ladder), Dict(pairs(pt)); keyedarrays = false)
        @test cl.converged
        @test isapprox(vec(collect(sol.nodeflux)), vec(collect(ref.nodeflux)); rtol = 1e-8)
    end

    @testset "the padded linear term is refilled exactly" begin
        # The padding, the augmentation and their union are built once per
        # cache and refilled with values; what the refill produces must be
        # entry for entry what a fresh assembly at the same values produces,
        # including the coupled inductor rows of the augmentation, which are
        # the only value dependent entries it has, and the drive.
        same(a::SparseMatrixCSC, b::SparseMatrixCSC) = size(a) == size(b) &&
            SparseArrays.getcolptr(a) == SparseArrays.getcolptr(b) &&
            rowvals(a) == rowvals(b) && nonzeros(a) == nonzeros(b)
        function refilled(c, w, Nh, sources, defs, p2; kw...)
            cache = hbcache(w, Nh, sources, c, defs; kw...)
            hbsolve!(cache, (;))
            lin = cache.reuse.linear
            hbsolve!(cache, p2)
            fresh = JosephsonCircuits.hbnlsolve(w, Nh, sources, c,
                merge(defs, Dict(pairs(p2))); returnsystem = true,
                keyedarrays = false, kw...)
            return cache.reuse.linear === lin && !isnothing(lin) &&
                same(lin.invLnm, fresh.invLnm) && same(lin.Gnm, fresh.Gnm) &&
                same(lin.Cnm, fresh.Cnm) && same(lin.Rbnm, fresh.Rbnm) &&
                lin.bnm == fresh.sys.bnm &&
                nonzeros(cache.reuse.sys.Knm) == nonzeros(fresh.sys.Knm)
        end
        function chain(Nj)
            c = Any[(:P1, 1, 0, Port(1; Z0 = 50.0))]
            for i in 1:Nj
                push!(c, (Symbol(:Lj, i), i, i + 1, JosephsonJunction(:Lj)))
                push!(c, (Symbol(:Cj, i), i, i + 1, Capacitor(55e-15)))
                push!(c, (Symbol(:Cg, i), i + 1, 0, Capacitor(:Cg)))
            end
            push!(c, (:P2, Nj + 1, 0, Port(2; Z0 = 50.0)))
            return Circuit(c)
        end
        mutual = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)),
            (:L1, 1, 2, Inductor(:L1)), (:L2, 2, 3, Inductor(2e-9)),
            (:L3, 3, 0, Inductor(4e-9)), (:K1, :L1, :L2, MutualInductor(:k)),
            (:Lj1, 1, 0, JosephsonJunction(1e-9)), (:C1, 1, 0, Capacitor(1e-12))])
        Lj0 = IctoLj(3.4e-6)
        @test refilled(chain(8), (2*pi*7e9,), (6,),
            [(mode=(1,), port=1, current=0.5e-6)],
            Dict(:Lj => Lj0, :Cg => 45e-15), (Lj = 1.1*Lj0, Cg = 50e-15))
        @test refilled(mutual, (2*pi*5e9,), (4,),
            [(mode=(1,), port=1, current=1e-7),
             (mode=(0,), port=1, current=2e-8)],
            Dict(:L1 => 1e-9, :k => 0.5), (L1 = 1.5e-9, k = 0.7);
            dc = true, odd = true, even = true)
        # a scattering block restamped at the scale of a new port impedance
        blocked = Circuit([:p1 => Port(1; Z0 = :Rp), :cc => Capacitor(100e-15),
            :jj => JosephsonJunction(1e-9),
            :c2 => ScatteringParameters(w -> fill((1 - im*w*5e-11)/(1 + im*w*5e-11), 1, 1);
                nports = 1, grounded = true)],
            [((:p1, 1), (:cc, 1)), ((:cc, 2), (:jj, 1), (:c2, 1)),
             ((:jj, 2), (:p1, 2), Ground)])
        @test refilled(blocked, (2*pi*4.75e9,), (6,),
            [(mode=(1,), port=1, current=1e-8)], Dict(:Rp => 50.0), (Rp = 30.0,))
    end
end
