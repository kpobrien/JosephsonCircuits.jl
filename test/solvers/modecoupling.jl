using JosephsonCircuits
using LinearAlgebra
using SparseArrays
using Random

isdefined(Main, :testchaincircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))
using Test

@testset verbose=true "the mode coupling preconditioner" begin

    circuit = Any[]
    push!(circuit,("P1", "1", "0", Port(1; Z0 = :Rleft)))
    push!(circuit,("C1", "1", "2", Capacitor(:Cc)))
    push!(circuit,("Lj1", "2", "0", JosephsonJunction(:Lj)))
    push!(circuit,("C2", "2", "0", Capacitor(:Cj)))
    circuit = Circuit(circuit)
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
                d.Ljb, d.Lscale, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq,
                d.invLnm, d.Gnm, d.Cnm, layout)

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
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; spec = JosephsonCircuits.CoupledModes([1, 2]))
        @test pc.coupling.indices == [1, 2]

        JosephsonCircuits.updatepreconditioner!(pc, 0.3*randn(length(d.xr)))
        ms = modeslot(d.modelayout)
        retained = [m in pc.coupling.indices for m in 1:Nmodes]
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
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout)
        @test pc0.coupling isa JosephsonCircuits.BlockDiagonal
        JosephsonCircuits.updatepreconditioner!(pc0, 0.3*randn(length(d.xr)))
        rows0 = rowvals(pc0.P)
        for j in axes(pc0.P, 2), k in nzrange(pc0.P, j)
            @test ms[rows0[k]] == ms[j]
        end
        @test nnz(pc0.P) < nnz(pc.P)

        # escalation goes straight to the full Jacobian, once
        @test JosephsonCircuits.escalatepreconditioner!(pc0)
        @test pc0.coupling isa JosephsonCircuits.FullJacobian
        @test pc0.escalations == 1
        @test !JosephsonCircuits.escalatepreconditioner!(pc0)

        # a band grows at every escalation until it is the full Jacobian,
        # also on this grid of odd harmonics, whose offsets are all even
        pb = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; spec = HarmonicBand((0,)),
            Amatrixmodes = d.Amatrixmodes)
        stored = [nnz(pb.P)]
        while JosephsonCircuits.escalatepreconditioner!(pb)
            push!(stored, nnz(pb.P))
        end
        @test all(>(0), diff(stored))
        @test pb.coupling isa JosephsonCircuits.FullJacobian

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
            d.Lscale, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm,
            d.Cnm, d.modelayout; spec = JosephsonCircuits.CoupledModes([0]))
    end

    @testset "hbnlsolve newtonkrylov agrees with newton" begin
        on = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = Newton(), keyedarrays = false)
        @test on.solverinfo.converged
        @test isempty(on.solverinfo.stages[1].krylov)

        for m in (NewtonKrylov(), NewtonKrylov(preconditioner = FullJacobian()),
                NewtonKrylov(preconditioner = Floquet(size = 20)),
                NewtonKrylov(preconditioner = CoupledModes([1, 3])))
            ok = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
                circuitdefs; method = m, keyedarrays = false)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-12*maximum(abs, on.nodeflux))
            st = ok.solverinfo.stages[1]
            @test length(st.krylov) >= st.iterations
        end
        @test_throws ArgumentError Floquet(size = 4, harvest = 0, ritz = 0)
    end

    @testset "deflation forms on a strongly driven chain" begin
        # a junction chain driven hard enough that the block diagonal alone
        # needs help, solved with escalation disabled so the recycled
        # subspace is what has to carry the solve, in both forms and with
        # the base refreshed eagerly or frozen across the Newton path
        chain = Any[]
        push!(chain, ("P1", "1", "0", Port(1; Z0 = :R)))
        Ncell = 12
        for i in 1:Ncell
            push!(chain, ("Lj$(i)", "$(i)", "$(i+1)", JosephsonJunction(:Lj)))
            push!(chain, ("C$(i)", "$(i)", "0", Capacitor(:Cg)))
        end
        push!(chain, ("C$(Ncell+1)", "$(Ncell+1)", "0", Capacitor(:Cg)))
        push!(chain, ("R2", "$(Ncell+1)", "0", Resistor(:R)))
        chain = Circuit(chain)
        chaindefs = Dict{Symbol,Complex{Float64}}(
            :Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
        wc = 2*pi*8e9
        # about 0.97 of the critical current at the port; the junction
        # phases reach ~1.4 rad and the block diagonal alone needs several
        # Arnoldi steps per solve, which is what the harvest feeds on
        chainsources = [(mode=(1,),port=1,current=3.2e-6)]
        on = JosephsonCircuits.hbnlsolve((wc,), (8,), chainsources, chain,
            chaindefs; method = Newton(), keyedarrays = false)
        @test on.solverinfo.converged
        phimax = maximum(abs, on.nodeflux)
        for (name, refresh) in (("eager", Always()), ("frozen", Never()))
            ok = JosephsonCircuits.hbnlsolve((wc,), (8,), chainsources, chain,
                chaindefs; keyedarrays = false, method = NewtonKrylov(
                    preconditioner = Floquet(size = 12, harvest = 4),
                    refresh = refresh, escalate = false))
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
        dccircuit = Any[]
        push!(dccircuit,("P1", "1", "0", Port(1; Z0 = :Rleft)))
        push!(dccircuit,("L1", "1", "0", Inductor(:Lm)))
        push!(dccircuit,("K1", "L1", "L2", MutualInductor(:K1)))
        push!(dccircuit,("C1", "1", "2", Capacitor(:Cc)))
        push!(dccircuit,("L2", "2", "3", Inductor(:Lm)))
        push!(dccircuit,("Lj3", "3", "0", JosephsonJunction(:Lj)))
        push!(dccircuit,("Lj4", "2", "0", JosephsonJunction(:Lj)))
        push!(dccircuit,("C2", "2", "0", Capacitor(:Cj)))
        dccircuit = Circuit(dccircuit)
        dcdefs = Dict{Symbol,Complex{Float64}}(
            :Lj => 2000e-12, :Lm => 10e-12, :Cc => 200.0e-15,
            :Cj => 900e-15, :Rleft => 50.0, :Rright => 50.0, :K1 => 0.9)
        dcsources = [(mode=(0,),port=1,current=50e-5),
                     (mode=(1,),port=1,current=0.0001e-6)]
        dcw = 2*pi*5e9

        on = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources, dccircuit,
            dcdefs; dc = true, method = Newton(), keyedarrays = false)
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
                d.Ljb, d.Lscale, d.Rbnm, Nmodes, d.Nbranches, d.Nfreq,
                d.invLnm, d.Gnm, d.Cnm, d.modelayout)
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

        for m in (NewtonKrylov(), NewtonKrylov(preconditioner = Floquet(size = 10)))
            ok = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources,
                dccircuit, dcdefs; dc = true, method = m,
                keyedarrays = false)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-9*maximum(abs, on.nodeflux))
        end
        # With direct current injected the state is canonical and the
        # recycler wraps the canonical preconditioner; the harvest must
        # build a subspace there. One vector per solve is enough to show
        # it, with escalation off so the deflation is what is used.
        # The block diagonal is exact on this circuit (one Arnoldi step per
        # solve), so the benefit filter would rightly leave every harvested
        # direction inactive; it is switched off here so that the active
        # size shows the harvest reached the wrapper.
        let
            ok = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources,
                dccircuit, dcdefs; dc = true, keyedarrays = false,
                method = NewtonKrylov(preconditioner = Floquet(size = 6,
                    harvest = 1, benefittol = 0.0), escalate = false))
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-9*maximum(abs, on.nodeflux))
            st = ok.solverinfo.stages[1]
            @test any(k -> k.deflationsize > 0, st.krylov)
            @test st.krylov[end].deflationrebuilds > 0
            @test st.krylov[end].deflationproducts > 0
        end
    end

    @testset "the deflation subspace is inherited across a cached sweep" begin
        # a parameter sweep through hbcache rebinds the system and the
        # preconditioner rather than rebuilding them; the recycled subspace
        # of the previous point is what the next solve should start from
        jpa = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(100.0e-15)),
            (:Lj1, 2, 0, JosephsonJunction(:Lj)), (:C2, 2, 0, Capacitor(1000e-15))])
        let
            cache = JosephsonCircuits.hbcache((wp,), (8,), sources, jpa,
                Dict(:Lj => 1000e-12); method = NewtonKrylov(preconditioner =
                    Floquet(size = 6, harvest = 2), escalate = false))
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
            @test cache.reuse.recycling isa JosephsonCircuits.FloquetState
            # and the answer is the answer
            on = JosephsonCircuits.hbnlsolve((wp,), (8,), sources,
                jpa, Dict(:Lj => 1010e-12);
                method = Newton(), keyedarrays = false)
            @test isapprox(second.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-12*maximum(abs, on.nodeflux))
        end
    end

    @testset "a failed solve does not seed the next point" begin
        # the reuse object commits the candidates of a converged solve only:
        # a solve cut off after one Newton step leaves the previous state in
        # place, and the next converged solve starts from that state
        jpa = Circuit([(:P1, 1, 0, Port(1; Z0 = 50.0)), (:C1, 1, 2, Capacitor(100.0e-15)),
            (:Lj1, 2, 0, JosephsonJunction(:Lj)), (:C2, 2, 0, Capacitor(1000e-15))])
        cache = JosephsonCircuits.hbcache((wp,), (8,), sources, jpa,
            Dict(:Lj => 1000e-12); method = NewtonKrylov(preconditioner =
                Floquet(size = 6, harvest = 2), escalate = false))
        first = JosephsonCircuits.hbsolve!(cache, (; Lj = 1000e-12))
        @test first.solverinfo.converged
        committed = cache.reuse.recycling
        @test committed isa JosephsonCircuits.FloquetState
        @test size(committed.X, 2) > 0
        Xcommitted = copy(committed.X)
        # the solve `hbsolve!` makes, cut off after one Newton step, which
        # it says
        failed = @test_logs (:warn,) match_mode=:any JosephsonCircuits.hbnlsolve(
            cache.w, cache.sources,
            cache.frequencies, cache.indices, cache.compiled,
            cache.nm; keyedarrays = false, reuse = cache.reuse,
            iterations = 1, cache.kwargs...)
        @test !failed.solverinfo.converged
        @test failed.solverinfo.stages[end].reason == :iterations
        @test cache.reuse.recycling === committed
        @test cache.reuse.recycling.X == Xcommitted
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
            method = Newton(), keyedarrays = false)
        @test on.solverinfo.converged
        let
            ok = JosephsonCircuits.hbnlsolve((w1,w2), (6,3), src2, circuit,
                circuitdefs; dc = true, odd = true, even = true,
                keyedarrays = false, method = NewtonKrylov(preconditioner =
                    Floquet(size = 8, harvest = 2), escalate = false))
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
            circuitdefs; method = Newton(), keyedarrays = false)
        # a solver object carrying preconditioner options, which are not
        # options of the inner Newton-Krylov loop and must not be forwarded
        # to it
        ok = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = JosephsonCircuits.NewtonKrylov(
                preconditioner = Floquet(size = 6, harvest = 2),
                linearsolver = GMRES(restart = 50)), keyedarrays = false)
        @test ok.solverinfo.converged
        @test isapprox(ok.nodeflux, on.nodeflux;
            rtol = 1e-6, atol = 1e-12*maximum(abs, on.nodeflux))
        @test any(k -> k.deflationsize > 0, ok.solverinfo.stages[1].krylov)
        # and through hbsolve
        hs = JosephsonCircuits.hbsolve(2*pi*4.5e9, (wp,), sources, (1,), (8,),
            circuit, circuitdefs; method = NewtonKrylov(preconditioner =
                Floquet(size = 6, harvest = 2)), keyedarrays = false)
        @test hs.nonlinear.solverinfo.converged
        @test any(k -> k.deflationsize > 0,
            hs.nonlinear.solverinfo.stages[1].krylov)
    end

    @testset "KrylovSolveInfo diagnostics" begin
        o = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
            circuitdefs; method = NewtonKrylov(), keyedarrays = false)
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
            circuitdefs; method = Newton(), keyedarrays = false)
        @test isempty(on.solverinfo.stages[1].krylov)

        str = sprint(show, MIME("text/plain"), st.krylov[1])
        @test occursin("KrylovSolveInfo", str)
        @test occursin("eta=", str)
    end

    # the two tone chain with a direct current mode, debugged once for the
    # three testsets which build preconditioners on it
    circuit2, defs2 = testchaincircuit()
    wpb = 2*pi*4.75e9; wsb = 2*pi*5.0e9
    srcb = [(mode=(1,0), port=1, current=0.3e-6),
            (mode=(0,1), port=1, current=0.1e-6)]
    chain2 = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb, circuit2,
        defs2; debugJacobian = true, dc = true, odd = true, even = true,
        keyedarrays = false)

    @testset "block factorization over the circuit graph" begin
        d = chain2
        Nmodes = d.Nmodes
        n = length(d.xr)
        mk(spec; kw...) = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; spec = spec, Amatrixmodes = d.Amatrixmodes, kw...)
        x = 0.3*randn(Random.default_rng(), n)
        d.fjreal(nothing, d.Jr, x)
        r = randn(Random.default_rng(), n)
        z = similar(r)
        # The block factorization inverts its diagonal blocks explicitly
        # and does not pivot across supernodes, so its accuracy is set by
        # the conditioning of those blocks rather than of the matrix; the
        # residual is measured against that.
        blockcond(pre) = maximum(cond(Array(view(Dp, :, :, b)))
            for c in pre.P.clusters for Dp in c.lu.D for b in axes(Dp, 3))
        residualbound(pre, zz, T = Float64) =
            100*eps(T)*blockcond(pre)*norm(d.Jr)*norm(zz)

        # the symbolic pieces on a small chain: KLU's order is a
        # permutation with no fill, and along the chain the elimination
        # tree is a path which amalgamation merges into chains of the target
        adj = [[2], [1, 3], [2, 4], [3, 5], [4]]
        order = JosephsonCircuits.klunodeorder(adj)
        @test sort(order) == 1:5
        parent, post = JosephsonCircuits.eliminationtree(adj, order)
        @test count(==(0), parent) == 1
        @test sort(post) == 1:5
        parent, post = JosephsonCircuits.eliminationtree(adj, 1:5)
        @test parent == [2, 3, 4, 5, 0]
        @test post == [1, 2, 3, 4, 5]
        @test JosephsonCircuits.amalgamate(parent, post, fill(2, 5), 6) ==
            [[1, 2, 3], [4, 5]]
        @test JosephsonCircuits.amalgamate(parent, post, fill(2, 5), 100) ==
            [[1, 2, 3, 4, 5]]
        # a mask's clusters and singletons; the clusters are closed, which
        # is what the block factorization keeps
        mask = Matrix{Bool}(I, Nmodes, Nmodes)
        for (a, b) in ((1, 2), (2, 3), (1, 3), (4, 5))
            mask[a, b] = mask[b, a] = true
        end
        @test JosephsonCircuits.modeclusters(mask) == [[1, 2, 3], [4, 5]]
        @test JosephsonCircuits.singletonmodes(mask) == collect(6:Nmodes)

        # the full coupling set is an exact solve of the Jacobian
        pb = mk(FullJacobian(factorization = BlockFactorization()))
        @test pb.P isa JosephsonCircuits.BlockStructure
        @test length(pb.P.clusters) == 1
        @test isnothing(pb.P.singletons)
        JosephsonCircuits.updatepreconditioner!(pb, x)
        JosephsonCircuits.applypreconditioner!(z, pb, r)
        @test norm(d.Jr*z - r) <= residualbound(pb, z)
        @test JosephsonCircuits.isexactpreconditioner(pb)
        # on the host the substitutions run as loops and allocate nothing;
        # Julia 1.11 still allocates the view wrappers each `mul!` of a
        # slice is handed, which later releases elide
        if VERSION >= v"1.12"
            @test (@allocated JosephsonCircuits.applypreconditioner!(z, pb, r)) == 0
        end
        # in single precision it is a preconditioner
        p32 = mk(FullJacobian(factorization = BlockFactorization(;
            precision = Float32)))
        JosephsonCircuits.updatepreconditioner!(p32, x)
        JosephsonCircuits.applypreconditioner!(z, p32, r)
        # single precision leaves a fixed fraction of the right hand side,
        # or what the conditioning of its blocks allows when that is more
        @test norm(d.Jr*z - r) <=
            max(5e-2*norm(r), residualbound(p32, z, Float32))
        @test eltype(p32.P.clusters[1].lu.D[1]) == Float32

        # a mask made of clusters: the block solve equals the sparse
        # factorization of the same mask
        pm = mk(CouplingMask(mask))
        JosephsonCircuits.updatepreconditioner!(pm, x)
        zm = similar(r)
        JosephsonCircuits.applypreconditioner!(zm, pm, r)
        pc = mk(CouplingMask(mask; factorization = BlockFactorization()))
        @test [c.modes for c in pc.P.clusters] == [[1, 2, 3], [4, 5]]
        @test !isnothing(pc.P.singletons)
        JosephsonCircuits.updatepreconditioner!(pc, x)
        JosephsonCircuits.applypreconditioner!(z, pc, r)
        @test norm(z - zm) <= 100*eps()*blockcond(pc)*norm(zm)
        # a band is factorized on its closure, here everything
        pa = mk(HarmonicBand(1; factorization = BlockFactorization()))
        @test length(pa.P.clusters) == 1
        @test isnothing(pa.P.singletons)
        JosephsonCircuits.updatepreconditioner!(pa, x)
        JosephsonCircuits.applypreconditioner!(z, pa, r)
        @test norm(d.Jr*z - r) <= residualbound(pa, z)
        # escalation from a mask rebuilds the structure on the full set
        @test JosephsonCircuits.escalatepreconditioner!(pc)
        @test length(pc.P.clusters) == 1
        JosephsonCircuits.updatepreconditioner!(pc, x)
        JosephsonCircuits.applypreconditioner!(z, pc, r)
        @test norm(d.Jr*z - r) <= residualbound(pc, z)
        # a rebound structure gives the same solve
        JosephsonCircuits.rebind!(pb, d.sys)
        JosephsonCircuits.updatepreconditioner!(pb, x)
        JosephsonCircuits.applypreconditioner!(z, pb, r)
        @test norm(d.Jr*z - r) <= residualbound(pb, z)
        # and rebound to a system of other component values, the sparse
        # and the block factorizations of the full set invert that
        # system's Jacobian: the constant values the rebind leaves to the
        # refactorization are refreshed there
        d2 = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb, circuit2,
            merge(defs2, Dict(:Lj => 120e-12 + 0im, :Cg => 45e-15 + 0im));
            debugJacobian = true, dc = true, odd = true, even = true,
            keyedarrays = false)
        d2.fjreal(nothing, d2.Jr, x)
        for pr in (mk(FullJacobian()), pb)
            JosephsonCircuits.updatepreconditioner!(pr, x)
            JosephsonCircuits.rebind!(pr, d2.sys)
            JosephsonCircuits.updatepreconditioner!(pr, x)
            JosephsonCircuits.applypreconditioner!(z, pr, r)
            bound = pr.P isa JosephsonCircuits.BlockStructure ?
                100*eps()*blockcond(pr)*norm(d2.Jr)*norm(z) : 1e-10*norm(r)
            @test norm(d2.Jr*z - r) <= bound
        end

        # end to end, against Newton, with the rebuild decided by the count
        # rule and by the probe
        on = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb, circuit2,
            defs2; method = Newton(), dc = true, odd = true, even = true,
            keyedarrays = false)
        for m in (NewtonKrylov(preconditioner = FullJacobian(factorization = BlockFactorization())),
                NewtonKrylov(preconditioner = FullJacobian(factorization = BlockFactorization(; precision = Float32)), refresh = Probe()),
                NewtonKrylov(preconditioner = CouplingMask(mask; factorization = BlockFactorization())),
                NewtonKrylov(preconditioner = MeasuredBand(factorization = BlockFactorization()), refresh = Probe()))
            ok = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb,
                circuit2, defs2; method = m, dc = true, odd = true,
                even = true, keyedarrays = false)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux; rtol = 1e-6,
                atol = 1e-12*maximum(abs, on.nodeflux))
            st = ok.solverinfo.stages[1]
            # every step records whether the preconditioner was rebuilt
            @test count(k -> k.refreshed, st.krylov) >= 1
        end
        @test_throws TypeError NewtonKrylov(refresh = :bogus)
    end


    @testset "clusters measured from the operator" begin
        # the rule on a toy strength matrix: two strong pairs and weak
        # couplings elsewhere give two clusters and a contractive remainder
        W = fill(0.01, 6, 6)
        for i in 1:6; W[i, i] = 0; end
        W[1, 2] = W[2, 1] = 5.0
        W[3, 4] = W[4, 3] = 4.0
        cids, r0, r1, taken = JosephsonCircuits.spectralclusters(W)
        @test r0 > 1
        @test r1 < 1
        @test cids[1] == cids[2] && cids[3] == cids[4]
        @test cids[1] != cids[3] && cids[5] != cids[6]
        @test Set(taken) == Set([(1, 2), (3, 4)])
        m = [cids[i] == cids[j] for i in eachindex(cids), j in eachindex(cids)]
        @test m == m'
        @test JosephsonCircuits.modeclusters(m) == [[1, 2], [3, 4]]
        # nothing to merge when the couplings are already contractive
        cids0, _, _, taken0 = JosephsonCircuits.spectralclusters(0.01*ones(4, 4))
        @test length(unique(cids0)) == 4
        @test isempty(taken0)
        # a coupling strong one way and weak the other, whose power
        # iteration alternates: the radius is the one the eigenvalues give,
        # and the rule takes the coupling
        _, r0p, _, takenp = JosephsonCircuits.spectralclusters([0.0 4.0; 0.5 0.0])
        @test r0p ≈ sqrt(2) rtol = 1e-4
        @test takenp == [(1, 2)]
        # twenty strong pairs in a weak background: the rule stops at the
        # first merge whose remainder contracts, as the eigenvalues of the
        # remainder with those twenty merges and with one fewer show
        Nq = 200
        Wq = (0.2/Nq) .* rand(MersenneTwister(11), Nq, Nq)
        for p in 1:20
            Wq[2p-1, 2p] = Wq[2p, 2p-1] = 1.2
        end
        for i in 1:Nq; Wq[i, i] = 0; end
        cq, _, r1q, takenq = JosephsonCircuits.spectralclusters(Wq)
        between(ids) = [ids[i] == ids[j] ? 0.0 : Wq[i, j] for i in 1:Nq, j in 1:Nq]
        @test length(takenq) == 20
        @test r1q < 1
        @test maximum(abs, eigvals(between(cq))) < 1
        ids19 = collect(1:Nq)
        for (i, j) in takenq[1:19]
            ids19[j] = ids19[i]
        end
        @test maximum(abs, eigvals(between(ids19))) >= 1

        # the probe on a system: strengths are nonnegative with a zero
        # diagonal, and `stalled!` makes the next update remeasure
        d = chain2
        pc = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, d.Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; spec = Clusters(factorization = BlockFactorization()),
            Amatrixmodes = d.Amatrixmodes)
        pr = pc.clusterprobe
        @test pr.probes == 0
        x = 0.3*randn(length(d.xr))
        JosephsonCircuits.updatepreconditioner!(pc, x)
        @test pr.probes == 1
        @test all(>=(0), pr.W)
        @test all(iszero, pr.W[i, i] for i in 1:d.Nmodes)
        @test pc.coupling isa JosephsonCircuits.CouplingMask
        JosephsonCircuits.updatepreconditioner!(pc, x)
        @test pr.probes == 1
        JosephsonCircuits.stalled!(pc)
        JosephsonCircuits.updatepreconditioner!(pc, x)
        @test pr.probes == 2
        # the wrappers forward the notification
        pcw = JosephsonCircuits.SizedPreconditioner(pc, length(d.xr))
        JosephsonCircuits.stalled!(pcw)
        @test pr.reprobe
        @test_throws DimensionMismatch JosephsonCircuits.ModeCouplingPreconditioner(
            d.sys, d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb,
            d.Lscale, d.Rbnm, d.Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm,
            d.Cnm, d.modelayout; spec = CouplingMask(falses(2, 2)))

        # end to end, with a sparse and with the block factorization
        on = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb, circuit2,
            defs2; method = Newton(), dc = true, odd = true, even = true,
            keyedarrays = false)
        for f in (JosephsonCircuits.KLUfactorization(),
                JosephsonCircuits.BlockFactorization(),
                JosephsonCircuits.BlockFactorization(; precision = Float32))
            ok = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb,
                circuit2, defs2; dc = true, odd = true, even = true,
                keyedarrays = false, method = NewtonKrylov(preconditioner =
                    Clusters(factorization = f)))
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux; rtol = 1e-6,
                atol = 1e-12*maximum(abs, on.nodeflux))
            # the probe's products are among those the records count, one
            # per mode at every probe
            @test ok.solverinfo.stages[1].krylov[end].deflationproducts >= d.Nmodes
        end
    end


    @testset "Automatic picks by the tones and the memory" begin
        d = chain2
        sys = d.sys; Nmodes = d.Nmodes; layout = d.modelayout
        nnodes = layout.dim ÷ Nmodes
        ns = JosephsonCircuits.branchnodesandsigns(d.Rbnm, Nmodes, d.Nbranches)
        pairptr, pairrow, _, _ = JosephsonCircuits.junctionpairtable(Int32,
            Float32, sys.Ljb, ns, nnodes)
        adj = JosephsonCircuits.circuitnodegraph(pairptr, pairrow, sys.invLnm,
            sys.Gnm, sys.Cnm, Nmodes, nnodes)
        order = JosephsonCircuits.klunodeorder(adj)
        keep = JosephsonCircuits.modecouplingmask(Nmodes, 1:Nmodes)
        # the predictor counts every array the cluster allocates, exactly,
        # without allocating any of them
        pred = JosephsonCircuits.blockfactorbytes(Float32, keep, adj, order,
            Nmodes, layout)
        C = JosephsonCircuits.clusterblocks(Float32, 1:Nmodes, adj, order,
            Nmodes, layout, JosephsonCircuits.CPU())
        lu = C.lu
        held = (sum(length, lu.D) + sum(length, lu.L) + sum(length, lu.U) +
                sum(length, lu.Dinv))*sizeof(Float32) +
            sum(sizeof, values(lu.scratch); init = 0) +
            sizeof(C.z) + sizeof(C.w) + sizeof(C.tmp)
        @test pred == held
        @test JosephsonCircuits.blockfactorbytes(Float64, keep, adj, order,
            Nmodes, layout) == 2*pred
        @test JosephsonCircuits.freememory(JosephsonCircuits.CPU()) > 0
        @test_throws ArgumentError JosephsonCircuits.freememory(nothing)
        # one tone: the full Jacobian with the backend's sparse
        # factorization, whatever the memory margin; two tones: the full
        # block factors in single precision when they fit, the measured
        # band when they do not
        res(modes; kw...) = JosephsonCircuits.resolveautomatic(sys, d.Rbnm,
            Nmodes, d.Nbranches, layout, modes, JosephsonCircuits.CPU(); kw...)
        @test res(nothing) isa FullJacobian
        @test res(nothing).factorization === nothing    # the backend's default
        @test res(nothing; budget = 0) isa FullJacobian
        @test res(d.Amatrixmodes) isa FullJacobian &&
            res(d.Amatrixmodes).factorization isa BlockFactorization
        @test res(d.Amatrixmodes).factorization.precision == Float32
        @test res(d.Amatrixmodes; budget = pred) isa FullJacobian
        @test res(d.Amatrixmodes; budget = pred - 1) isa MeasuredBand
        # the factors may be held in another precision than the
        # iteration: only cuDSS can, since KLU and UMFPACK are compiled for
        # double and would promote, and a `BlockFactorization` sizes its
        # own dense blocks rather than the matrix it is handed
        @test JosephsonCircuits.factorizationprecision(KLUfactorization()) === nothing
        @test JosephsonCircuits.factorizationprecision(LUfactorization()) === nothing
        @test JosephsonCircuits.factorizationprecision(
            BlockFactorization(precision = Float32)) === nothing
        @test CUDSSFactorization().precision === nothing
        @test CUDSSFactorization(precision = Float32).precision === Float32
        @test JosephsonCircuits.factorizationprecision(
            CUDSSFactorization(precision = Float32)) === Float32
        # and a preconditioner whose factors are single precision still
        # inverts a double precision residual, through the conversion in
        # `applypreconditioner!`
        pc32 = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            layout; spec = FullJacobian(), precision = Float32,
            Amatrixmodes = d.Amatrixmodes)
        x32 = 0.3*randn(length(d.xr))
        JosephsonCircuits.updatepreconditioner!(pc32, x32)
        d.fjreal(nothing, d.Jr, x32)
        @test eltype(pc32.P) === Float32
        r32 = randn(length(d.xr))
        z32 = similar(r32)
        JosephsonCircuits.applypreconditioner!(z32, pc32, r32)
        @test eltype(z32) === Float64          # the iteration keeps its own
        @test d.Jr*z32 ≈ r32 rtol=1e-4         # single precision accuracy

        # an Automatic carries no factorization of its own: the member it
        # resolves to takes the backend's default
        @test JosephsonCircuits.withfactorization(Automatic(),
            KLUfactorization()) === Automatic()
        # and the solve through it agrees with the assembled Newton solve
        on = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb, circuit2,
            defs2; dc = true, odd = true, even = true, method = Newton())
        ok = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb, circuit2,
            defs2; dc = true, odd = true, even = true,
            method = NewtonKrylov(preconditioner = Automatic()))
        @test ok.solverinfo.converged
        @test isapprox(on.S, ok.S; rtol = 1e-6)
        pc = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            layout; spec = Automatic(), Amatrixmodes = d.Amatrixmodes)
        @test pc.coupling isa FullJacobian
        # the full set in single precision is full but not exact: its
        # escalation is the double precision factorization
        @test JosephsonCircuits.isfullcoupling(pc)
        @test !JosephsonCircuits.isexactpreconditioner(pc)
        one = JosephsonCircuits.hbnlsolve((wpb,), (8,),
            [(mode=(1,), port=1, current=0.3e-6)], circuit2, defs2;
            method = NewtonKrylov(preconditioner = Automatic()))
        @test one.solverinfo.converged
        # the memory prediction of any coupling set, from the structure
        # alone, and the escalation budget it is held to
        @test JosephsonCircuits.couplingbytes(pc, FullJacobian(BlockFactorization(precision = Float32))) == pred
        @test JosephsonCircuits.couplingbytes(pc, FullJacobian()) > JosephsonCircuits.couplingbytes(pc, BlockDiagonal()) > 0
        pcb = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            layout; spec = BlockDiagonal(), Amatrixmodes = d.Amatrixmodes)
        @test pcb.budget === nothing
        pcb.budget = 0
        @test !JosephsonCircuits.escalatepreconditioner!(pcb)
        @test pcb.coupling isa BlockDiagonal
        @test pcb.escalations == 0
        pcb.budget = JosephsonCircuits.couplingbytes(pcb, FullJacobian())
        @test JosephsonCircuits.escalatepreconditioner!(pcb)
        @test pcb.coupling isa FullJacobian
        @test pcb.escalations == 1
        # a measured band and a cluster mask grow at an update only within
        # the same budget: with room for the set they start from and no
        # more, a strongly driven point leaves them as they are, where it
        # grows them without the budget
        xg = randn(MersenneTwister(5), length(d.xr))
        for spec in (MeasuredBand(), Clusters())
            pg = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
                d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
                d.Rbnm, Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
                layout; spec = spec, Amatrixmodes = d.Amatrixmodes)
            start = pg.coupling
            pg.budget = JosephsonCircuits.couplingbytes(pg, start)
            JosephsonCircuits.updatepreconditioner!(pg, xg)
            @test pg.coupling === start
            zg = similar(xg)
            JosephsonCircuits.applypreconditioner!(zg, pg, xg)
            @test all(isfinite, zg)
            pg.budget = nothing
            spec isa Clusters && JosephsonCircuits.stalled!(pg)
            JosephsonCircuits.updatepreconditioner!(pg, xg)
            @test pg.coupling !== start
        end

        # the sparse factors are sized under the ordering the factorization
        # takes: on the band of a two tone chain that is nested dissection,
        # whose factors AMD's analysis overestimated twofold
        chain40 = Any[("P1", "1", "0", Port(1; Z0 = :R))]
        for i in 1:40
            push!(chain40, ("Lj$(i)", "$(i)", "$(i+1)", JosephsonJunction(:Lj)),
                ("C$(i)", "$(i)", "0", Capacitor(:Cg)))
        end
        push!(chain40, ("C41", "41", "0", Capacitor(:Cg)), ("R2", "41", "0", Resistor(:R)))
        d40 = JosephsonCircuits.hbnlsolve((2*pi*7e9, 2*pi*7.3e9), (6, 6),
            [(mode = (1, 0), port = 1, current = 1.0e-6),
             (mode = (0, 1), port = 1, current = 1.0e-6)], Circuit(chain40), defs2;
            method = Newton(), iterations = 0, debugJacobian = true,
            keyedarrays = false)
        pband = JosephsonCircuits.ModeCouplingPreconditioner(d40.sys,
            d40.Amatrixindicesaliased, d40.Amatrixconjindices, d40.Ljb,
            d40.Lscale, d40.Rbnm, d40.Nmodes, d40.Nbranches, d40.Nfreq,
            d40.invLnm, d40.Gnm, d40.Cnm, d40.modelayout;
            spec = HarmonicBand((1, 1)), Amatrixmodes = d40.Amatrixmodes)
        predicted = JosephsonCircuits.couplingbytes(pband, pband.coupling)
        JosephsonCircuits.updatepreconditioner!(pband, 0.3*randn(length(d40.xr)))
        built = nnz(pband.cache.factorization)*(sizeof(Float64) + sizeof(Int))
        @test 0.8 < predicted/built < 1.25
    end

    @testset "factors held in less precision than the iteration" begin
        # a factorization holding single precision factors of a double
        # precision iteration is a preconditioner, not an exact solve,
        # whether it is the block factorization or cuDSS: the
        # classification is by the precision of the factors against the
        # iteration's, and the escalation is the same coupling set in the
        # iteration's precision, with the factorization's settings kept.
        # The cuDSS factorization itself needs a device; the
        # classification, the sizing and the rebuild of the matrix do not.
        d = chain2
        mk(spec; kw...) = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, d.Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; spec = spec, Amatrixmodes = d.Amatrixmodes, kw...)
        @test JosephsonCircuits.iterationprecision(d.sys) === Float64
        f32 = CUDSSFactorization(precision = Float32, ir_n_steps = 0)
        pc = mk(FullJacobian(f32))
        @test JosephsonCircuits.factorprecision(pc) === Float32
        @test eltype(pc.P) === Float32
        @test JosephsonCircuits.isfullcoupling(pc)
        @test JosephsonCircuits.reducedprecision(pc)
        @test !JosephsonCircuits.isexactpreconditioner(pc)
        # the promotion is held to the budget like any escalation
        pc.budget = 0
        @test !JosephsonCircuits.escalatepreconditioner!(pc)
        @test pc.factorization === f32
        pc.budget = nothing
        @test JosephsonCircuits.escalatepreconditioner!(pc)
        @test pc.escalations == 1
        @test pc.factorization isa CUDSSFactorization
        @test pc.factorization.precision === Float64
        @test pc.factorization.kwargs == (; ir_n_steps = 0)
        @test pc.coupling.factorization === pc.factorization
        @test pc.plan.precision === Float64
        @test JosephsonCircuits.factorprecision(pc) === Float64
        @test eltype(pc.P) === Float64
        @test JosephsonCircuits.isexactpreconditioner(pc)
        @test !JosephsonCircuits.escalatepreconditioner!(pc)
        # the preconditioner's own precision keyword is the same case for
        # a factorization which has none of its own
        pk = mk(FullJacobian(); precision = Float32)
        @test eltype(pk.P) === Float32
        @test !JosephsonCircuits.isexactpreconditioner(pk)
        # KLU factorizes the single precision matrix in double precision,
        # and solves a double precision residual in it without converting
        xk = 0.3*randn(length(d.xr))
        JosephsonCircuits.updatepreconditioner!(pk, xk)
        rk = randn(length(d.xr)); zk = similar(rk)
        JosephsonCircuits.applypreconditioner!(zk, pk, rk)
        @test norm(Matrix{Float64}(pk.P)*zk - rk) <= 1e-9*norm(rk)
        # the application allocates nothing beyond the solve of its
        # factorization: no conversion and no scratch vector. The solve is
        # the reference rather than zero because KLU.jl's allocates the
        # reference it passes its settings in on some Julia releases.
        fk = pk.cache.factorization
        JosephsonCircuits.trysolve!(zk, fk, rk)
        @test (@allocated JosephsonCircuits.applypreconditioner!(zk, pk, rk)) ==
            (@allocated JosephsonCircuits.trysolve!(zk, fk, rk))
        @test JosephsonCircuits.escalatepreconditioner!(pk)
        @test pk.factorization isa KLUfactorization
        @test eltype(pk.P) === Float64
        @test JosephsonCircuits.isexactpreconditioner(pk)
        # single precision block factors of the full set promote to
        # double, and the promotion is sized as the double precision
        # factors it builds
        pb = mk(FullJacobian(BlockFactorization(precision = Float32)))
        @test !JosephsonCircuits.isexactpreconditioner(pb)
        pb.budget = JosephsonCircuits.couplingbytes(pb,
            FullJacobian(BlockFactorization(precision = Float64))) - 1
        @test !JosephsonCircuits.escalatepreconditioner!(pb)
        pb.budget += 1
        @test JosephsonCircuits.escalatepreconditioner!(pb)
        @test pb.factorization.precision === Float64
        @test pb.coupling.factorization === pb.factorization
        @test JosephsonCircuits.isexactpreconditioner(pb)
        # a direct solve factorizes in the precision of its iteration and
        # refuses a factorization asking for another
        @test_throws ArgumentError Newton(factorization = f32)
        @test_throws ArgumentError QuasiNewton(factorization = f32)
        @test_throws "precision" JosephsonCircuits.factorize(f32,
            sparse(1.0I, 2, 2))
        @test JosephsonCircuits.withprecision(KLUfactorization(), Float64) isa
            KLUfactorization
        # and factors in the iteration's own precision are exact, whatever
        # it is
        s = JosephsonCircuits.hbnlsolve((wpb, wsb), (4, 2), srcb, circuit2,
            defs2; debugJacobian = true, dc = true, odd = true, even = true,
            keyedarrays = false, method = NewtonKrylov(precision = Float32))
        @test JosephsonCircuits.iterationprecision(s.sys) === Float32
        for spec in (FullJacobian(CUDSSFactorization(precision = Float32)),
                     FullJacobian(BlockFactorization(precision = Float32)))
            p = JosephsonCircuits.ModeCouplingPreconditioner(s.sys,
                s.Amatrixindicesaliased, s.Amatrixconjindices, s.Ljb, s.Lscale,
                s.Rbnm, s.Nmodes, s.Nbranches, s.Nfreq, s.invLnm, s.Gnm,
                s.Cnm, s.modelayout; spec = spec, Amatrixmodes = s.Amatrixmodes)
            @test JosephsonCircuits.factorprecision(p) === Float32
            @test !JosephsonCircuits.reducedprecision(p)
            @test JosephsonCircuits.isexactpreconditioner(p)
            @test !JosephsonCircuits.escalatepreconditioner!(p)
        end
    end

    @testset "a singular supernode falls back to the sparse factorization" begin
        # a flux biased SQUID with a weak second tone at its resonance: the
        # node between the loop inductor and a junction has no stiffness of
        # its own at that frequency (its junction and capacitor resonate and
        # the loop inductor's stiffness lives in the promoted branch
        # current), so the supernode holding it alone has a singular
        # diagonal block. Double precision survived on roundoff, single
        # precision does not; the preconditioner must switch to the sparse
        # factorization and converge rather than throw
        case = crosscheckcases()[7]
        ws = case.ws[2]
        srcs = [(mode = (s.mode..., 0), port = s.port, current = s.current)
            for s in case.sources]
        push!(srcs, (mode = (0, 1), port = case.signalport, current = case.Is))
        # (at six harmonics; the singular supernode is a property of the
        # mode grid, and the table's own count is not pinned to it)
        sol = hbnlsolve((case.wp..., ws), (6, 6), srcs,
            case.circuit, case.defs; nonlinearkw(case.kw)..., atol = 1e-12,
            method = NewtonKrylov(preconditioner =
                FullJacobian(BlockFactorization(precision = Float32))))
        @test sol.solverinfo.converged
        ref = hbnlsolve((case.wp..., ws), (6, 6), srcs,
            case.circuit, case.defs; nonlinearkw(case.kw)..., atol = 1e-12,
            method = Newton())
        @test maximum(abs, Array(sol.S) .- Array(ref.S)) < 1e-9
        # the fallback keeps to the budget an escalation keeps to: with no
        # room for the sparse factors of the full set, it takes the mode
        # block diagonal
        d = hbnlsolve((case.wp..., ws), (6, 6), srcs, case.circuit,
            case.defs; nonlinearkw(case.kw)..., method = Newton(),
            iterations = 0, debugJacobian = true)
        pb = JosephsonCircuits.ModeCouplingPreconditioner(d.sys,
            d.Amatrixindicesaliased, d.Amatrixconjindices, d.Ljb, d.Lscale,
            d.Rbnm, d.Nmodes, d.Nbranches, d.Nfreq, d.invLnm, d.Gnm, d.Cnm,
            d.modelayout; Amatrixmodes = d.Amatrixmodes,
            spec = FullJacobian(BlockFactorization(precision = Float32)))
        pb.budget = 0
        JosephsonCircuits.updatepreconditioner!(pb, d.xr)
        @test pb.factorization isa KLUfactorization
        @test pb.coupling isa BlockDiagonal
        @test pb.fallbacks == 1

        # at four harmonics the same single precision factors are not
        # singular but poor: the Krylov solves stagnate, and the escalation
        # the driver asks for is the double precision factorization of the
        # same full set, after which the solve converges
        sol4 = hbnlsolve((case.wp..., ws), (4, 4), srcs, case.circuit,
            case.defs; nonlinearkw(case.kw)..., atol = 1e-12,
            method = NewtonKrylov(preconditioner =
                FullJacobian(BlockFactorization(precision = Float32))))
        @test sol4.solverinfo.converged
        @test count(k -> k.escalated, sol4.solverinfo.stages[end].krylov) >= 1
        ref4 = hbnlsolve((case.wp..., ws), (4, 4), srcs, case.circuit,
            case.defs; nonlinearkw(case.kw)..., atol = 1e-12, method = Newton())
        @test maximum(abs, Array(sol4.S) .- Array(ref4.S)) < 1e-9
    end

end
