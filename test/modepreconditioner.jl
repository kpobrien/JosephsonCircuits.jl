using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Test

@testset verbose=true "modepreconditioner" begin

    circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
    push!(circuit,("P1","1","0",1))
    push!(circuit,("R1","1","0",:Rleft))
    push!(circuit,("C1","1","2",:Cc))
    push!(circuit,("Lj1","2","0",:Lj))
    push!(circuit,("C2","2","0",:Cj))
    circuitdefs = Dict{Symbol,Complex{Float64}}(
        :Lj => 1000e-12, :Cc => 100.0e-15, :Cj => 1000e-15, :Rleft => 50.0)
    wp = 2*pi*5e9
    sources = [(mode=(1,),port=1,current=1.0e-6)]

    modeslot(layout) = Int[(Int(layout.inv[j]) - 1) % layout.nmodes + 1
        for j in 1:layout.rdim]

    @testset "modecouplingmask" begin
        @test JosephsonCircuits.modecouplingmask(3, Int[]) == Matrix(I, 3, 3)
        @test all(JosephsonCircuits.modecouplingmask(3, 1:3))
        # a retained set keeps its whole columns plus the diagonal, which is
        # the block lower triangular Gauss-Seidel pattern
        keep = JosephsonCircuits.modecouplingmask(4, [2, 3])
        @test keep == Bool[1 1 1 0; 0 1 1 0; 0 1 1 0; 0 1 1 1]
        @test !keep[1, 4]
        @test keep[4, 2]
        @test_throws ArgumentError JosephsonCircuits.modecouplingmask(3, [4])
        @test_throws ArgumentError JosephsonCircuits.modecouplingmask(0, Int[])
    end

    @testset "restrictmodecoupling" begin
        A = [1 -3 5; 3 1 -3; 5 3 1]
        @test JosephsonCircuits.restrictmodecoupling(A,
            JosephsonCircuits.modecouplingmask(3, [1])) == [1 0 0; 3 1 0; 5 0 1]
        @test JosephsonCircuits.restrictmodecoupling(A,
            JosephsonCircuits.modecouplingmask(3, 1:3)) == A
        @test_throws DimensionMismatch JosephsonCircuits.restrictmodecoupling(
            A, trues(2, 2))
    end

    @testset "the restricted assembly equals the masked full Jacobian" begin
        d = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; debugJacobian = true)
        Nmodes = d.Nmodes
        layout = d.modelayout
        ms = modeslot(layout)

        for S in (Int[], [1], [1, 3], collect(1:Nmodes))
            keep = JosephsonCircuits.modecouplingmask(Nmodes, S)
            P, plan = JosephsonCircuits.structurejacobian(d,
                JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixindicesaliased, keep),
                JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixconjindices, keep),
                d.Ljb, d.Lmean, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq,
                d.invLnm, d.Gnm, d.Cnm, layout, layout)

            x = 0.3*randn(length(d.xr))
            d.fjreal(nothing, d.Jr, x)
            JosephsonCircuits.setpoint!(d.sys, x)
            JosephsonCircuits.jacobian!(P, plan, d.sys)

            Jref = copy(d.Jr)
            rows = rowvals(Jref)
            vals = nonzeros(Jref)
            for j in axes(Jref, 2), k in nzrange(Jref, j)
                keep[ms[rows[k]], ms[j]] || (vals[k] = 0.0)
            end
            @test Matrix(P) == Matrix(Jref)
            @test nnz(P) <= nnz(d.Jr)
        end
    end

    @testset "block lower triangular structure and escalation" begin
        d = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; debugJacobian = true)
        Nmodes = d.Nmodes
        pc = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lmean,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; couplingmodes = [1, 2])
        @test pc.couplingmodes == [1, 2]

        JosephsonCircuits.updatepreconditioner!(pc, 0.3*randn(length(d.xr)))
        ms = modeslot(d.modelayout)
        retained = [m in pc.couplingmodes for m in 1:Nmodes]
        rows = rowvals(pc.P)
        vals = nonzeros(pc.P)
        # no stored nonzero couples a shell column into a retained row
        for j in axes(pc.P, 2), k in nzrange(pc.P, j)
            if !retained[ms[j]] && ms[rows[k]] != ms[j]
                @test iszero(vals[k])
            end
        end

        # the block diagonal really is block diagonal
        pc0 = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lmean,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout)
        @test isempty(pc0.couplingmodes)
        JosephsonCircuits.updatepreconditioner!(pc0, 0.3*randn(length(d.xr)))
        rows0 = rowvals(pc0.P)
        for j in axes(pc0.P, 2), k in nzrange(pc0.P, j)
            @test ms[rows0[k]] == ms[j]
        end
        @test nnz(pc0.P) < nnz(pc.P)

        # escalation goes straight to the full Jacobian, once
        @test JosephsonCircuits.escalatepreconditioner!(pc0)
        @test pc0.couplingmodes == collect(1:Nmodes)
        @test pc0.escalations == 1
        @test !JosephsonCircuits.escalatepreconditioner!(pc0)

        # and the escalated preconditioner is the exact Jacobian
        x = 0.3*randn(length(d.xr))
        JosephsonCircuits.updatepreconditioner!(pc0, x)
        d.fjreal(nothing, d.Jr, x)
        @test Matrix(pc0.P) == Matrix(d.Jr)
        r = randn(length(d.xr))
        z = similar(r)
        JosephsonCircuits.applypreconditioner!(z, pc0, r)
        @test d.Jr*z ≈ r rtol=1e-8

        @test_throws ArgumentError JosephsonCircuits.ModeCouplingPreconditioner(
            d.sys, d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb,
            d.Lmean, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm,
            d.Cnm, d.modelayout; couplingmodes = :bogus)
    end

    @testset "hbnlsolve newtonkrylov agrees with newton" begin
        on = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = :newton, keyedarrays = false)
        @test on.solverinfo.converged
        @test isempty(on.solverinfo.stages[1].krylov)

        for kw in ((;), (; krylovcouplingmodes = :all),
                   (; krylovrecycle = 20), (; krylovcouplingmodes = [1, 3]),
                   (; krylovrecycle = 20, krylovdeflationform = :adef2))
            ok = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
                circuitdefs; method = :newtonkrylov, keyedarrays = false, kw...)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-12*maximum(abs, on.nodeflux))
            st = ok.solverinfo.stages[1]
            @test length(st.krylov) >= st.iterations
        end
        @test_throws ArgumentError JosephsonCircuits.hbnlsolve((wp,), (8,),
            sources, circuit, circuitdefs; method = :newtonkrylov,
            keyedarrays = false, krylovrecycle = 4,
            krylovdeflationform = :bogus)
    end

    @testset "deflation forms on a strongly driven chain" begin
        # a junction chain driven hard enough that the block diagonal alone
        # needs help, solved with escalation disabled so the recycled
        # subspace is what has to carry the solve, in both forms and with
        # the base refreshed eagerly or frozen across the Newton path
        chain = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(chain, ("P1","1","0",1)); push!(chain, ("R1","1","0",:R))
        Ncell = 12
        for i in 1:Ncell
            push!(chain, ("Lj$(i)","$(i)","$(i+1)",:Lj))
            push!(chain, ("C$(i)","$(i)","0",:Cg))
        end
        push!(chain, ("C$(Ncell+1)","$(Ncell+1)","0",:Cg))
        push!(chain, ("R2","$(Ncell+1)","0",:R))
        chaindefs = Dict{Symbol,Complex{Float64}}(
            :Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
        wc = 2*pi*8e9
        # about 0.97 of the critical current at the port; the junction
        # phases reach ~1.4 rad and the block diagonal alone needs several
        # Arnoldi steps per solve, which is what the harvest feeds on
        chainsources = [(mode=(1,),port=1,current=3.2e-6)]
        on = JosephsonCircuits.hbnlsolve((wc,), (8,), chainsources, chain,
            chaindefs; method = :newton, keyedarrays = false)
        @test on.solverinfo.converged
        phimax = maximum(abs, on.nodeflux)
        noescalation = (; krylovescalate = typemax(Int))
        frozen = (; krylovescalate = typemax(Int),
            krylovrefreshiterations = typemax(Int), krylovrefreshrate = 1.0)
        for form in (:adef1, :adef2), (name, kk) in
                (("eager", noescalation), ("frozen", frozen))
            ok = JosephsonCircuits.hbnlsolve((wc,), (8,), chainsources, chain,
                chaindefs; method = :newtonkrylov, keyedarrays = false,
                krylovrecycle = 12, krylovharvest = 4,
                krylovdeflationform = form, krylovkwargs = kk)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-12*phimax)
            st = ok.solverinfo.stages[1]
            # no escalation was allowed, so the base stayed block diagonal
            @test !any(k -> k.escalated, st.krylov)
            # a subspace was harvested, and it was built into the
            # preconditioner: under a frozen base that only happens through
            # the lazy refresh, which is what the count checks
            @test any(k -> k.deflationsize > 0, st.krylov)
            @test st.krylov[end].deflationrebuilds > 0
        end
    end

    @testset "a circuit with a dc mode (real self-conjugate modes)" begin
        # the real representation collapses self-conjugate modes to a single
        # slot, which the layout handling and the restricted assembly must
        # both survive
        dccircuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
        push!(dccircuit,("P1","1","0",1))
        push!(dccircuit,("R1","1","0",:Rleft))
        push!(dccircuit,("L1","1","0",:Lm))
        push!(dccircuit,("K1","L1","L2",:K1))
        push!(dccircuit,("C1","1","2",:Cc))
        push!(dccircuit,("L2","2","3",:Lm))
        push!(dccircuit,("Lj3","3","0",:Lj))
        push!(dccircuit,("Lj4","2","0",:Lj))
        push!(dccircuit,("C2","2","0",:Cj))
        dcdefs = Dict{Symbol,Complex{Float64}}(
            :Lj => 2000e-12, :Lm => 10e-12, :Cc => 200.0e-15,
            :Cj => 900e-15, :Rleft => 50.0, :Rright => 50.0, :K1 => 0.9)
        dcsources = [(mode=(0,),port=1,current=50e-5),
                     (mode=(1,),port=1,current=0.0001e-6)]
        dcw = 2*pi*5e9

        on = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources, dccircuit,
            dcdefs; dc = true, method = :newton, keyedarrays = false)
        @test on.solverinfo.converged

        d = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources, dccircuit,
            dcdefs; dc = true, debugJacobian = true)
        @test any(d.modelayout.isreal)
        Nmodes = d.Nmodes
        ms = modeslot(d.modelayout)
        for S in (Int[], [1], collect(1:Nmodes))
            keep = JosephsonCircuits.modecouplingmask(Nmodes, S)
            P, plan = JosephsonCircuits.structurejacobian(d,
                JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixindicesaliased, keep),
                JosephsonCircuits.restrictmodecoupling(
                    d.Amatrixconjindices, keep),
                d.Ljb, d.Lmean, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq,
                d.invLnm, d.Gnm, d.Cnm, d.modelayout, d.modelayout)
            x = 0.3*randn(length(d.xr))
            d.fjreal(nothing, d.Jr, x)
            JosephsonCircuits.setpoint!(d.sys, x)
            JosephsonCircuits.jacobian!(P, plan, d.sys)
            Jref = copy(d.Jr)
            rows = rowvals(Jref)
            vals = nonzeros(Jref)
            for j in axes(Jref, 2), k in nzrange(Jref, j)
                keep[ms[rows[k]], ms[j]] || (vals[k] = 0.0)
            end
            @test Matrix(P) == Matrix(Jref)
        end

        for kw in ((;), (; krylovrecycle = 10),
                   (; krylovrecycle = 10, krylovdeflationform = :adef2))
            ok = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources,
                dccircuit, dcdefs; dc = true, method = :newtonkrylov,
                keyedarrays = false, kw...)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-9*maximum(abs, on.nodeflux))
        end
        # With direct current injected the state is canonical and the
        # recycler wraps the canonical preconditioner; the harvest must
        # build a subspace there (it used to sit inside the wrapper, where
        # no harvest reached it and recycling was inert on every such
        # circuit). One vector per solve is enough to show it, with
        # escalation off so the deflation is what is used.
        for form in (:adef1, :adef2)
            ok = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources,
                dccircuit, dcdefs; dc = true, method = :newtonkrylov,
                keyedarrays = false, krylovrecycle = 6, krylovharvest = 1,
                krylovdeflationform = form,
                krylovkwargs = (; krylovescalate = typemax(Int)))
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-9*maximum(abs, on.nodeflux))
            st = ok.solverinfo.stages[1]
            @test any(k -> k.deflationsize > 0, st.krylov)
            @test st.krylov[end].deflationrebuilds > 0
        end
    end

    @testset "the deflation subspace is inherited across a cached sweep" begin
        # a parameter sweep through hbcache rebinds the system and the
        # preconditioner rather than rebuilding them; the recycled subspace
        # of the previous point is what the next solve should start from
        builder(; Lj) = [("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","2",100.0e-15), ("Lj1","2","0",Lj), ("C2","2","0",1000e-15)]
        for form in (:adef1, :adef2)
            cache = JosephsonCircuits.hbcache((wp,), (8,), sources, builder,
                (; Lj = 1000e-12); krylovrecycle = 6, krylovharvest = 2,
                krylovdeflationform = form,
                krylovkwargs = (; krylovescalate = typemax(Int)))
            first = JosephsonCircuits.hbsolve!(cache, (; Lj = 1000e-12))
            @test first.solverinfo.converged
            k1 = first.solverinfo.stages[1].krylov
            # a cold start has nothing to deflate at its first solve
            @test k1[1].deflationsize == 0
            @test any(k -> k.deflationsize > 0, k1)
            second = JosephsonCircuits.hbsolve!(cache, (; Lj = 1010e-12))
            @test second.solverinfo.converged
            k2 = second.solverinfo.stages[1].krylov
            # the next point starts from the inherited subspace, rebuilt
            # against the rebound base
            @test k2[1].deflationsize > 0
            @test k2[1].deflationrebuilds >= 1
            @test cache.reuse.recycling isa JosephsonCircuits.RecyclingState
            # and the answer is the answer
            on = JosephsonCircuits.hbnlsolve((wp,), (8,), sources,
                builder(; Lj = 1010e-12), Dict{Symbol,Complex{Float64}}();
                method = :newton, keyedarrays = false)
            @test isapprox(second.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-12*maximum(abs, on.nodeflux))
        end
    end

    @testset "a failed solve does not seed the next point" begin
        # the reuse object commits the candidates of a converged solve only:
        # a solve cut off after one Newton step leaves the previous state in
        # place, and the next converged solve starts from that state
        builder(; Lj) = [("P1","1","0",1), ("R1","1","0",50.0),
            ("C1","1","2",100.0e-15), ("Lj1","2","0",Lj), ("C2","2","0",1000e-15)]
        cache = JosephsonCircuits.hbcache((wp,), (8,), sources, builder,
            (; Lj = 1000e-12); krylovrecycle = 6, krylovharvest = 2,
            krylovkwargs = (; krylovescalate = typemax(Int)))
        first = JosephsonCircuits.hbsolve!(cache, (; Lj = 1000e-12))
        @test first.solverinfo.converged
        committed = cache.reuse.recycling
        @test committed isa JosephsonCircuits.RecyclingState
        @test size(committed.U, 2) > 0
        Xcommitted = copy(committed.U)
        # the solve `hbsolve!` makes, cut off after one Newton step
        failed = JosephsonCircuits.hbnlsolve(cache.w, cache.sources,
            cache.frequencies, cache.indices, cache.compiled, cache.cg,
            cache.nm; keyedarrays = false, reuse = cache.reuse,
            iterations = 1, cache.kwargs...)
        @test !failed.solverinfo.converged
        @test cache.reuse.recycling === committed
        @test cache.reuse.recycling.U == Xcommitted
        again = JosephsonCircuits.hbsolve!(cache, (; Lj = 1005e-12))
        @test again.solverinfo.converged
        @test again.solverinfo.stages[1].krylov[1].deflationsize > 0
        @test cache.reuse.recycling !== committed
    end

    @testset "two tones with a direct current block, both forms" begin
        # the recycler wraps the canonical preconditioner, so its candidates
        # are corrections of the whole canonical state; both forms must
        # reach the direct solve's answer with the deflation active
        w1 = 2*pi*5e9; w2 = 2*pi*1.19e9
        src2 = [(mode=(1,0),port=1,current=3.0e-7),
                (mode=(0,1),port=1,current=2.4e-7)]
        on = JosephsonCircuits.hbnlsolve((w1,w2), (6,3), src2, circuit,
            circuitdefs; dc = true, odd = true, even = true,
            method = :newton, keyedarrays = false)
        @test on.solverinfo.converged
        for form in (:adef1, :adef2)
            ok = JosephsonCircuits.hbnlsolve((w1,w2), (6,3), src2, circuit,
                circuitdefs; dc = true, odd = true, even = true,
                method = :newtonkrylov, keyedarrays = false,
                krylovrecycle = 8, krylovharvest = 2,
                krylovdeflationform = form,
                krylovkwargs = (; krylovescalate = typemax(Int)))
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-9*maximum(abs, on.nodeflux))
            kr = ok.solverinfo.stages[1].krylov
            @test any(k -> k.deflationsize > 0, kr)
            @test kr[end].deflationrebuilds > 0
            @test all(k -> k.products >= k.iterations + k.cycles, kr)
        end
    end

    @testset "the recycling options travel through hbsolve and NewtonKrylov" begin
        on = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = :newton, keyedarrays = false)
        # a solver object carrying preconditioner options, which are not
        # options of the inner Newton-Krylov loop and used to be forwarded
        # to it regardless
        ok = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = JosephsonCircuits.NewtonKrylov(
                krylovrecycle = 6, krylovharvest = 2,
                krylovdeflationform = :adef2, krylovcouplingmodes = :none,
                krylovrestart = 50), keyedarrays = false)
        @test ok.solverinfo.converged
        @test isapprox(ok.nodeflux, on.nodeflux;
            rtol = 1e-6, atol = 1e-12*maximum(abs, on.nodeflux))
        @test any(k -> k.deflationsize > 0, ok.solverinfo.stages[1].krylov)
        # and through hbsolve, which used to drop them
        hs = JosephsonCircuits.hbsolve(2*pi*4.5e9, (wp,), sources, (1,), (8,),
            circuit, circuitdefs; krylovrecycle = 6, krylovharvest = 2,
            krylovdeflationform = :adef2, keyedarrays = false)
        @test hs.nonlinear.solverinfo.converged
        @test any(k -> k.deflationsize > 0,
            hs.nonlinear.solverinfo.stages[1].krylov)
    end

    @testset "KrylovSolveInfo diagnostics" begin
        o = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = :newtonkrylov, keyedarrays = false)
        st = o.solverinfo.stages[1]
        @test st.converged
        @test !isempty(st.krylov)
        @test all(k -> k isa JosephsonCircuits.KrylovSolveInfo, st.krylov)
        # one record per linear solve, at least one per Newton step, and the
        # repeats are exactly the retries and rescues
        @test length(st.krylov) >= st.iterations
        @test issorted([k.iteration for k in st.krylov])
        @test all(k -> k.role in (:step, :retry, :rescue), st.krylov)
        for i in 2:length(st.krylov)
            if st.krylov[i].iteration == st.krylov[i-1].iteration
                @test st.krylov[i].role != :step
            end
        end
        @test JosephsonCircuits.krylovtotaliterations(st.krylov) ==
            sum(k.iterations for k in st.krylov)
        @test JosephsonCircuits.krylovtotaliterations(
            JosephsonCircuits.KrylovSolveInfo[]) == 0
        for k in st.krylov
            @test k.iterations >= 0
            @test k.cycles >= 1
            @test k.reason in (:converged, :breakdown, :stagnation,
                :iterationlimit)
            @test isfinite(k.normF) && k.normF >= 0
            @test 0 < k.forcing <= 1
            @test k.time >= 0
        end
        stepped = filter(k -> isfinite(k.slope), st.krylov)
        @test !isempty(stepped)
        @test all(k -> k.slope < 0, stepped)
        @test all(k -> 0 < k.alpha <= 1, stepped)
        @test issorted([k.time for k in st.krylov])

        on = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = :newton, keyedarrays = false)
        @test isempty(on.solverinfo.stages[1].krylov)

        str = sprint(show, MIME("text/plain"), st.krylov[1])
        @test occursin("KrylovSolveInfo", str)
        @test occursin("eta=", str)
    end

end
