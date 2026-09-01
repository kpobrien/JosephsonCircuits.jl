using JosephsonCircuits
using LinearAlgebra
using Test

isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "testcircuits.jl"))

@testset verbose=true "builders, providers and design sensitivities" begin

    wp = (2*pi*4.75001*1e9,)
    src = [(mode=(1,), port=1, current=0.00565e-6)]
    ws = 2*pi*(4.5:0.1:5.0)*1e9

    @testset "FrequencyDependent provider" begin
        # a lossy frequency dependent resistor with an arbitrary law (trig
        # inside the closure -- the generality a closed expression set
        # cannot offer)
        law = w -> 50.0*(1 + 0.1*sin(w/2e10)) + im*abs(w)*1e-9
        circuit = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", FrequencyDependent(law)),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        out = hbnlsolve(wp, (8,), src, circuit;
            keyedarrays = false)
        @test out.solverinfo.converged
        o2 = hbsolve(ws, wp, src, (2,), (8,), circuit)
        @test all(isfinite, Array(o2.linearized.S))

        # a CONSTANT law must agree exactly with a plain numeric value
        circa = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", FrequencyDependent(w -> 50.0)),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        circb = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", 50.0),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        Sa = hbsolve(ws, wp, src, (2,), (8,), circa).linearized.S
        Sb = hbsolve(ws, wp, src, (2,), (8,), circb).linearized.S
        @test isapprox(Array(Sa), Array(Sb), rtol = 1e-12)

        # providers combine with numbers through the value arithmetic
        circc = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0", 2*FrequencyDependent(w -> 25.0)),
            ("C1","1","2",100e-15), ("Lj1","2","0",1000e-12),
            ("C2","2","0",1000e-15)]
        Sc = hbsolve(ws, wp, src, (2,), (8,), circc).linearized.S
        @test isapprox(Array(Sc), Array(Sb), rtol = 1e-12)

        # the linear-only entry point also accepts a circuit alone
        lin = hblinsolve(ws, circb; Nmodulationharmonics = (2,))
        @test all(isfinite, Array(lin.S))
    end

    @testset "designjacobian" begin
        builder(; L, C) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",C),
            ("Lj1","2","0",L), ("C2","2","0",C/4)]
        p = (L = 1000.0e-12, C = 400.0e-15)
        names, v0, J = JosephsonCircuits.designjacobian(builder, p)
        # the port is not among them: its slot in the flat table holds the
        # port number, which is a label rather than a quantity to difference
        @test names == ["R1","C1","Lj1","C2"]
        # exact rows: dC1/dC = 1, dC2/dC = 1/4, dLj1/dL = 1, others zero
        iC1 = findfirst(==("C1"), names); iC2 = findfirst(==("C2"), names)
        iL = findfirst(==("Lj1"), names); iR = findfirst(==("R1"), names)
        jL = 1; jC = 2
        @test isapprox(J[iC1, jC], 1.0; rtol = 1e-6)
        @test isapprox(J[iC2, jC], 0.25; rtol = 1e-6)
        @test isapprox(J[iL, jL], 1.0; rtol = 1e-6)
        @test iszero(J[iR, jL]) && iszero(J[iR, jC])

        # a topology change under perturbation is an error, not a garbage
        # derivative
        badbuilder(; L, C) = L > 995.0e-12 ?
            builder(; L, C) :
            Tuple{String,String,String,Any}[("P1","1","0",1),
                ("R1","1","0",50.0), ("Lj1","2","0",L)]
        @test_throws ArgumentError JosephsonCircuits.designjacobian(
            badbuilder, (L = 1000.0e-12, C = 400.0e-15); delta = 1e-2)
    end

    @testset "designsensitivities matches full-solve finite differences" begin
        make_jpa(; Lj, Cc, Cj) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",Cc),
            ("Lj1","2","0",Lj), ("C2","2","0",Cj),
            ("C3","2","0",Cj/4)]     # derived value: chain-rule factor 1/4
        p = (Lj = 1000.0e-12, Cc = 100.0e-15, Cj = 800.0e-15)
        r = designsensitivities(make_jpa, p, ws, wp, src, (2,), (8,))
        @test size(r.dSdp, 5) == 3

        for q in (:Cj, :Lj)
            h = 1e-6*getproperty(p, q)
            pp = merge(p, NamedTuple{(q,)}((getproperty(p, q) + h,)))
            pm = merge(p, NamedTuple{(q,)}((getproperty(p, q) - h,)))
            Sp = hbsolve(ws, wp, src, (2,), (8,),
                make_jpa(; pp...)).linearized.S
            Sm = hbsolve(ws, wp, src, (2,), (8,),
                make_jpa(; pm...)).linearized.S
            fd = (Array(Sp) .- Array(Sm))./(2h)
            mine = Array(r.dSdp(parameter = q))
            @test norm(vec(mine) .- vec(fd))/norm(vec(fd)) < 1e-4
        end

        # a parameter no component depends on is rejected
        make_jpa2(; Lj, Cc, Cj, unused = 0.0) = make_jpa(; Lj, Cc, Cj)
        @test_throws ArgumentError designsensitivities(make_jpa2,
            merge(p, (unused = 1.0,)), ws, wp, src, (2,), (8,);
            parameters = (:unused,))
    end

    @testset "rotating directions and shared parameters are exact" begin
        # a design parameter which ROTATES a complex component value (a
        # loss tangent) has a direction which is not parallel to the value,
        # so a relative (logarithmic) derivative cannot represent it: the
        # solver must carry the direction dv/dp itself. Also here: one
        # parameter (Ic) touching two components of different kinds (the
        # junction inductance and its capacitance), which exercises the
        # per-parameter stamp merging and the mixed-kind grouping.
        phi0 = 3.29105976e-16
        make(; Ic, adens, t) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","2",100e-15/(1 + im*t)),
            ("Lj1","2","0",phi0/Ic),
            ("C2","2","0",Ic*adens)]
        p = (Ic = phi0/1000e-12, adens = 1000e-15/(phi0/1000e-12), t = 2e-3)
        r = designsensitivities(make, p, ws, wp, src, (2,), (8,))
        for q in (:Ic, :adens, :t)
            h = 1e-6*abs(getproperty(p, q))
            pp = merge(p, NamedTuple{(q,)}((getproperty(p, q) + h,)))
            pm = merge(p, NamedTuple{(q,)}((getproperty(p, q) - h,)))
            Sp = hbsolve(ws, wp, src, (2,), (8,),
                make(; pp...)).linearized.S
            Sm = hbsolve(ws, wp, src, (2,), (8,),
                make(; pm...)).linearized.S
            fd = vec((Array(Sp) .- Array(Sm))./(2h))
            mine = vec(Array(r.dSdp(parameter = q)))
            @test norm(mine .- fd)/norm(fd) < 1e-4
        end
    end

    @testset "reverse contraction agrees with forward" begin
        # the reverse contraction accumulates the operating point shift into
        # the same design parameter slots as the forward one
        make(; Lj, Cc) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",Cc),
            ("Lj1","2","0",Lj), ("C2","2","0",800e-15),
            ("C3","2","3",100e-15), ("Lj2","3","0",2*Lj),
            ("C4","3","0",800e-15)]
        p = (Lj = 1000.0e-12, Cc = 100.0e-15)
        rf = designsensitivities(make, p, ws, wp, src, (2,), (8,);
            sensitivitymode = :forward)
        rr = designsensitivities(make, p, ws, wp, src, (2,), (8,);
            sensitivitymode = :reverse)
        a = vec(Array(rf.dSdp)); b = vec(Array(rr.dSdp))
        @test norm(a .- b)/norm(a) < 1e-8
    end

    @testset "port impedance design parameter" begin
        make(; R, Cc) = Tuple{String,String,String,Any}[
            ("P1","1","0",1), ("R1","1","0",R), ("C1","1","2",Cc),
            ("Lj1","2","0",1000e-12), ("C2","2","0",1000e-15)]
        p = (R = 50.0, Cc = 100.0e-15)
        r = designsensitivities(make, p, ws, wp, src, (2,), (8,))
        h = 1e-6*p.R
        Sp = hbsolve(ws, wp, src, (2,), (8,),
            make(R = p.R + h, Cc = p.Cc)).linearized.S
        Sm = hbsolve(ws, wp, src, (2,), (8,),
            make(R = p.R - h, Cc = p.Cc)).linearized.S
        fd = vec((Array(Sp) .- Array(Sm))./(2h))
        mine = vec(Array(r.dSdp(parameter = :R)))
        @test norm(mine .- fd)/norm(fd) < 1e-4
    end

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
    end

    @testset "scattering block design parameters" begin
        # a theta-dependent series-impedance two port feeding a JPA: the
        # block parameter and a lumped parameter share the output axis,
        # and the block moves the pump operating point
        Z0 = 50.0
        function makeblk(; theta, Lj, analytic = false)
            seriesS(w) = (z = 1/(im*w*theta*Z0);
                [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
            dS(w) = (z = 1/(im*w*theta*Z0); d = -2*z/(theta*(z+2)^2);
                [d -d; -d d])
            blk = analytic ?
                ScatteringParameters(seriesS; nports = 2, grounded = false,
                    derivatives = (theta = dS,)) :
                ScatteringParameters(seriesS; nports = 2, grounded = false)
            return Circuit(
                [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :cc => blk,
                 :jj => JosephsonJunction(Lj),
                 :c2 => Capacitor(1000.0e-15)],
                [((:p1, 1), (:r1, 1), (:cc, 1, 1)),
                 ((:cc, 2, 1), (:jj, 1), (:c2, 1)),
                 ((:cc, 1, 2), (:cc, 2, 2), (:jj, 2), (:c2, 2), (:r1, 2),
                  (:p1, 2), Ground)])
        end
        p = (theta = 100.0e-15, Lj = 1000.0e-12)
        r = designsensitivities((; theta, Lj) -> makeblk(; theta, Lj),
            p, ws, wp, src, (2,), (8,))
        for q in (:theta, :Lj)
            h = 1e-6*abs(getproperty(p, q))
            pp = merge(p, NamedTuple{(q,)}((getproperty(p, q) + h,)))
            pm = merge(p, NamedTuple{(q,)}((getproperty(p, q) - h,)))
            Sp = hbsolve(ws, wp, src, (2,), (8,),
                makeblk(; pp...)).linearized.S
            Sm = hbsolve(ws, wp, src, (2,), (8,),
                makeblk(; pm...)).linearized.S
            fd = vec((Array(Sp) .- Array(Sm))./(2h))
            mine = vec(Array(r.dSdp(parameter = q)))
            @test norm(mine .- fd)/norm(fd) < 1e-4
        end

        # the analytic escape hatch agrees with the finite difference
        # provider to far below the finite difference error
        ra = designsensitivities(
            (; theta, Lj) -> makeblk(; theta, Lj, analytic = true),
            p, ws, wp, src, (2,), (8,))
        a = vec(Array(ra.dSdp(parameter = :theta)))
        f2 = vec(Array(r.dSdp(parameter = :theta)))
        @test norm(a .- f2)/norm(a) < 1e-8

        # a hoisted block is the identical object across builder calls, so
        # it is parameter independent and costs nothing
        fixedS(w) = (z = 1/(im*w*100e-15*Z0);
            [z/(z+2) 2/(z+2); 2/(z+2) z/(z+2)])
        hoisted = ScatteringParameters(fixedS; nports = 2, grounded = false)
        makehoist(; Lj) = Circuit(
            [:p1 => Port(1; termination = nothing), :r1 => Resistor(50.0), :cc => hoisted,
             :jj => JosephsonJunction(Lj), :c2 => Capacitor(1000.0e-15)],
            [((:p1, 1), (:r1, 1), (:cc, 1, 1)),
             ((:cc, 2, 1), (:jj, 1), (:c2, 1)),
             ((:cc, 1, 2), (:cc, 2, 2), (:jj, 2), (:c2, 2), (:r1, 2),
              (:p1, 2), Ground)])
        @test isempty(JosephsonCircuits.designblockjacobian(makehoist,
            (Lj = 1000.0e-12,)))
        rh = designsensitivities(makehoist, (Lj = 1000.0e-12,),
            ws, wp, src, (2,), (8,))
        @test all(isfinite, Array(rh.dSdp))

        # the reverse contraction rejects block parameters
        @test_throws ArgumentError designsensitivities(
            (; theta, Lj) -> makeblk(; theta, Lj), p, ws, wp, src,
            (2,), (8,); sensitivitymode = :reverse)
    end
end
