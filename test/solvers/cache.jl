using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))

# The solve cache: what it reuses between solves, and the exactness of
# the refill of the linear term it keeps.
@testset verbose=true "hbcache" begin

    wp = (2*pi*4.75001*1e9,)
    src = [(mode=(1,), port=1, current=0.00565e-6)]
    ws = 2*pi*(4.5:0.1:5.0)*1e9

    @testset "hbcache reuses the parse and warm starts" begin
        make(; Lj, Cc) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",Cc),
            ("Lj1","2","0",Lj), ("C2","2","0",1000e-15)]
        p = (Lj = 1000.0e-12, Cc = 100.0e-15)
        cache = hbcache(wp, (8,), src, make, p; ftol = 1e-12)
        nl = hbsolve!(cache, p; warmstart = false)
        @test cache.converged
        @test cache.nsolves == 1
        ref = hbnlsolve(wp, (8,), src, make(; p...); ftol = 1e-12,
            keyedarrays = false)
        @test isapprox(vec(collect(nl.nodeflux)),
            vec(collect(ref.nodeflux)); rtol = 1e-8)

        # a warm start from the stored point lands in the same place, in
        # fewer iterations
        coldits = sum(st.iterations for st in nl.solverinfo.stages)
        nw = hbsolve!(cache, p)
        @test cache.converged
        @test sum(st.iterations for st in nw.solverinfo.stages) <= coldits
        @test isapprox(vec(collect(nw.nodeflux)),
            vec(collect(nl.nodeflux)); rtol = 1e-8)

        # a nearby parameter re-solves against the reference
        p2 = (Lj = 1050.0e-12, Cc = 100.0e-15)
        n2 = hbsolve!(cache, p2)
        @test cache.converged
        ref2 = hbnlsolve(wp, (8,), src, make(; p2...); ftol = 1e-12,
            keyedarrays = false)
        @test isapprox(vec(collect(n2.nodeflux)),
            vec(collect(ref2.nodeflux)); rtol = 1e-6)
        # the port's reference impedance is the resistor across it, as the
        # legacy netlist states it, not the port number the netlist writes
        # in the value slot; the waves are normalized to it
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
        ref3 = hbnlsolve(wp, (8,), src, make(; Lj = 975.0e-12,
            Cc = 110.0e-15); ftol = 1e-12, keyedarrays = false)
        @test isapprox(vec(collect(n3.nodeflux)),
            vec(collect(ref3.nodeflux)); rtol = 1e-8)
        @test isapprox(Array(n3.S), Array(ref3.S); rtol = 1e-8)

        # the assembled Jacobian of a rebound system is that of the new
        # values, not of the first point: the operating point of a cached
        # solve matches a cold one, and a direct Newton cache converges
        # like a cold solve
        opc = hbcache(wp, (8,), src, make, p; ftol = 1e-12,
            returnoperatingpoint = true)
        hbsolve!(opc, p; warmstart = false)
        o2 = hbsolve!(opc, p2).operatingpoint
        r2 = hbnlsolve(wp, (8,), src, make(; p2...); ftol = 1e-12,
            keyedarrays = false, returnoperatingpoint = true).operatingpoint
        @test isapprox(o2.jacobian, r2.jacobian; rtol = 1e-8)
        @test norm(o2.jacobian - r2.jacobian) < 1e-8*norm(r2.jacobian)
        nc = hbcache(wp, (8,), src, make, p; ftol = 1e-12, method = Newton())
        hbsolve!(nc, p; warmstart = false)
        nn = hbsolve!(nc, p2; warmstart = false)
        rn = hbnlsolve(wp, (8,), src, make(; p2...); ftol = 1e-12,
            keyedarrays = false, method = Newton())
        @test isapprox(vec(collect(nn.nodeflux)),
            vec(collect(rn.nodeflux)); rtol = 1e-10)
        @test sum(st.iterations for st in nn.solverinfo.stages) ==
            sum(st.iterations for st in rn.solverinfo.stages)

        # a reset discards the stored state
        JosephsonCircuits.reset!(cache)
        @test isnothing(cache.x)
        @test !cache.converged

        # a topology change is an error, not a wrong answer
        make2(; Lj, Cc) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0),
            ("Lj1","2","0",Lj)]
        cache.builder = make2
        @test_throws ArgumentError hbsolve!(cache, p)
        # and so is a builder which returns one name twice, which would
        # fill one slot twice and leave another stale
        make3(; Lj, Cc) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",Cc),
            ("Lj1","2","0",Lj), ("C1","2","0",1000e-15)]
        cache.builder = make3
        @test_throws ArgumentError hbsolve!(cache, p)

        # the keywords the cache manages, an unsupported method and a
        # keyword the compiled solve does not take are refused at
        # construction, and a solver keyword still reaches the solve
        for bad in ((x0 = zeros(2),), (keyedarrays = true,),
                (reuse = nothing,), (method = Staged(),), (nosuchkeyword = 1,),
                (maxharmonics = (8,),))
            @test_throws ArgumentError hbcache(wp, (8,), src, make, p; bad...)
        end
        # what the cache does anyway is accepted
        @test hbcache(wp, (8,), src, make, p; keyedarrays = false) isa JosephsonCircuits.HBCache
        loose = hbcache(wp, (8,), src, make, p; ftol = 1e-2, iterations = 2)
        hbsolve!(loose, p; warmstart = false)
        @test loose.nsolves == 1
        @test sum(st.iterations for st in
            hbsolve!(loose, p; warmstart = false).solverinfo.stages) <= 2
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
        function refilled(builder, w, Nh, sources, p1, p2; kw...)
            cache = hbcache(w, Nh, sources, builder, p1; sorting = :name,
                kw...)
            hbsolve!(cache, p1)
            lin = cache.reuse.linear
            hbsolve!(cache, p2)
            fresh = JosephsonCircuits.hbnlsolve(w, Nh, sources,
                builder(; p2...), Dict{Symbol,Number}(); returnsystem = true,
                keyedarrays = false, sorting = :name, kw...)
            return cache.reuse.linear === lin && !isnothing(lin) &&
                same(lin.invLnm, fresh.invLnm) && same(lin.Gnm, fresh.Gnm) &&
                same(lin.Cnm, fresh.Cnm) && same(lin.Rbnm, fresh.Rbnm) &&
                lin.bnm == fresh.sys.bnm &&
                nonzeros(cache.reuse.sys.Knm) == nonzeros(fresh.sys.Knm)
        end
        function chain(; Lj, Cg = 45e-15, Nj = 8)
            c = Tuple{String,String,String,Any}[("P1","1","0",1),
                ("R1","1","0",50.0)]
            for i in 1:Nj
                push!(c, ("Lj$i","$i","$(i+1)",Lj))
                push!(c, ("Cj$i","$i","$(i+1)",55e-15))
                push!(c, ("Cg$i","$(i+1)","0",Cg))
            end
            push!(c, ("P2","$(Nj+1)","0",2)); push!(c, ("R2","$(Nj+1)","0",50.0))
            return c
        end
        mutual(; L1, k) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0), ("L1","1","2",L1),
            ("L2","2","3",2e-9), ("L3","3","0",4e-9), ("K1","L1","L2",k),
            ("Lj1","1","0",1e-9), ("C1","1","0",1e-12)]
        Lj0 = IctoLj(3.4e-6)
        @test refilled(chain, (2*pi*7e9,), (6,),
            [(mode=(1,), port=1, current=0.5e-6)],
            (Lj = Lj0,), (Lj = 1.1*Lj0, Cg = 50e-15))
        @test refilled(mutual, (2*pi*5e9,), (4,),
            [(mode=(1,), port=1, current=1e-7),
             (mode=(0,), port=1, current=2e-8)],
            (L1 = 1e-9, k = 0.5), (L1 = 1.5e-9, k = 0.7);
            dc = true, odd = true, even = true)
    end
end
