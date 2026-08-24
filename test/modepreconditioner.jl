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
            P, plan = JosephsonCircuits.planrealjacobian(
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
                   (; krylovrecycle = 20), (; krylovcouplingmodes = [1, 3]))
            ok = JosephsonCircuits.hbnlsolve((wp,), (8,), sources, circuit,
                circuitdefs; method = :newtonkrylov, keyedarrays = false, kw...)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-12*maximum(abs, on.nodeflux))
            st = ok.solverinfo.stages[1]
            @test length(st.krylov) >= st.iterations
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
            P, plan = JosephsonCircuits.planrealjacobian(
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

        for kw in ((;), (; krylovrecycle = 10))
            ok = JosephsonCircuits.hbnlsolve((dcw,), (2,), dcsources,
                dccircuit, dcdefs; dc = true, method = :newtonkrylov,
                keyedarrays = false, kw...)
            @test ok.solverinfo.converged
            @test isapprox(ok.nodeflux, on.nodeflux;
                rtol = 1e-6, atol = 1e-9*maximum(abs, on.nodeflux))
        end
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
