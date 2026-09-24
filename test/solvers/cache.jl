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
        # built, a value crossing a structural boundary when it is solved
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

        # the assembled Jacobian of a rebound system is that of the new
        # values, not of the first point: the operating point of a cached
        # solve matches a cold one, and a direct Newton cache converges
        # like a cold solve
        opc = hbcache(wp, (8,), src, circuit, defs; atol = 1e-12,
            returnoperatingpoint = true)
        hbsolve!(opc, p; warmstart = false)
        o2 = hbsolve!(opc, p2).operatingpoint
        r2 = hbnlsolve(wp, (8,), src, circuit, at(; p2...); atol = 1e-12,
            keyedarrays = false, returnoperatingpoint = true).operatingpoint
        @test isapprox(o2.jacobian, r2.jacobian; rtol = 1e-8)
        @test norm(o2.jacobian - r2.jacobian) < 1e-8*norm(r2.jacobian)
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
                (maxharmonics = (8,),))
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
        saved = (copy(first.nodeflux), copy(first.S))
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
        @test controlled.starts[3] == freshsolver.starts[1]
        @test recovered.nodeflux ≈ fresh.nodeflux rtol = 1e-10
        @test cache.matrixworkspace === matrixwork &&
            cache.nm.Cnm === capacitance
        @test first.nodeflux == saved[1] && first.S == saved[2]
        # a point refused before the assembly leaves the cache as it was,
        # so the next point still solves
        retained = copy(cache.x)
        @test_throws ArgumentError hbsolve!(cache, (Lj = Inf,))
        @test_throws ArgumentError hbsolve!(cache,
            (Cc = 100e-15*(1 - im/100),))
        @test cache.x == retained && cache.nsolves == 3
        @test hbsolve!(cache, (Lj = 0.98e-9,)).nodeflux ≈ fresh.nodeflux
        JosephsonCircuits.reset!(cache)
        @test isnothing(cache.x) && !cache.converged
        hbsolve!(cache, (Lj = 0.98e-9,))
        @test controlled.starts[end] == freshsolver.starts[1]
        @test first.nodeflux == saved[1] && first.S == saved[2]
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
        saved = (copy(first.nodeflux), copy(first.S))
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
        @test first.nodeflux == saved[1] && first.S == saved[2]
    end

    @testset "the definitions at a point" begin
        # a parameter defined under its parameter object, its symbol and its
        # string moves under every one of them, a parameter the definitions
        # lack is added under its symbol, and the base definitions stay
        JosephsonCircuits.@params La
        d = Dict{Any,Any}(La => 1.0, :La => 1.0, "La" => 1.0, :Cb => 2.0)
        at = JosephsonCircuits.definitionsat(d, JosephsonCircuits.definitionkeys(d),
            (La = 3.0, Cb = 4.0, Ln = 5.0))
        @test at[La] == at[:La] == at["La"] == 3.0
        @test at[:Cb] == 4.0 && at[:Ln] == 5.0
        @test d[La] == 1.0 && d[:Cb] == 2.0
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
    end
end
