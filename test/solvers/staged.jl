isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "..", "testcircuits.jl"))
using JosephsonCircuits, Test, LinearAlgebra
isdefined(Main, :recovery_solver) || include("recoveryfixture.jl")

@testset "method = Staged()" begin
    circuit, defs = testchaincircuit()
    w1 = 2*pi*5.0e9; w2 = 2*pi*1.19e9
    src = [(mode=(1,0), port=1, current=1.0e-6),
           (mode=(0,1), port=1, current=0.5e-6)]

    @testset "grid ladder" begin
        @test JosephsonCircuits.defaultgridladder((8,4)) ==
            [(2,2), (4,2), (8,4)]
        @test JosephsonCircuits.defaultgridladder((2,2)) == [(2,2)]
        @test JosephsonCircuits.defaultgridladder((20,)) ==
            [(2,), (3,), (5,), (10,), (20,)]
        # a tone with fewer than two harmonics keeps its own
        @test JosephsonCircuits.defaultgridladder((8,1)) ==
            [(2,1), (4,1), (8,1)]
    end

    @testset "hbsolve integration" begin
        # offset from the 5 GHz pump so no signal + pump mode lands at
        # zero total frequency
        ws = 2*pi*(4.55:0.3:5.5)*1e9
        ra = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs;
            dc = true, threewavemixing = true, fourwavemixing = true)
        rb = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs;
            dc = true, threewavemixing = true, fourwavemixing = true,
            method = Staged())
        a = vec(Array(ra.linearized.S)); b = vec(Array(rb.linearized.S))
        @test norm(a - b)/norm(a) < 1e-6
    end

    @testset "stage records in solverinfo" begin
        r = hbnlsolve((w1,w2), (8,4), src, circuit, defs;
            dc = true, odd = true, even = true, method = Staged())
        st = r.solverinfo.stages
        @test !isempty(st)
        @test all(x -> x isa JosephsonCircuits.StagedStageInfo, st)
        # every attempt appears; the last is the accepted full-drive solve
        # on the finest grid
        @test st[end].action === :final
        @test st[end].accepted && st[end].converged
        @test st[end].grid == (8, 4)
        @test st[end].starget == 1.0
        # the walk starts on the coarsest grid of the default ladder
        @test st[1].grid == JosephsonCircuits.defaultgridladder((8, 4))[1]
        # accepted advances carry the drive monotonically upward
        adv = filter(x -> x.accepted && x.action === :advance, st)
        @test all(x -> x.starget > x.sfrom, adv)
        # each record carries its inner solver stages and a wall time
        @test all(x -> !isempty(x.inner) && x.seconds > 0, st)
        @test all(x -> x.inner[end] isa JosephsonCircuits.IterationInfo, st)
        @test r.solverinfo.converged
        @test r.solverinfo.finalresidual == st[end].finalresidual
    end

    @testset "rejected steps and failed grid growth recover" begin
        c, d = testjpacircuitnumeric()
        pump = (2*pi*4.75001e9,)
        drive = [(mode = (1,), port = 1, current = 2e-9)]
        # half the drive converges, the whole drive fails, three quarters
        # and then the whole drive converge, and the first attempt at the
        # whole drive on the larger grid fails: the schedule retreats on
        # that grid and comes back to the whole drive there
        controlled = recovery_solver(failat = [2, 5])
        r = hbnlsolve(pump, (4,), drive, c, d; keyedarrays = false,
            method = Staged(inner = controlled.method, grids = [(2,), (4,)],
                s0 = 0.5, smin = 0.1))
        stages = r.solverinfo.stages
        @test r.solverinfo.converged
        @test findall(s -> !s.accepted, stages) == [2, 5]
        @test controlled.starts[2] ≈ controlled.solutions[1]
        @test controlled.starts[3] ≈ controlled.solutions[1]
        @test controlled.starts[6] == controlled.starts[5]
        @test stages[6].action == :grow && stages[6].starget < 1
        @test stages[end].grid == (4,) && stages[end].starget == 1
        @test stages[end].accepted && stages[end].action == :final
        @test length(controlled.starts[5]) > length(controlled.starts[4])
        @test isnan(r.solverinfo.sourcefold)
        fresh = hbnlsolve(pump, (4,), drive, c, d; method = Newton(),
            keyedarrays = false, atol = 1e-12)
        @test r.nodeflux ≈ fresh.nodeflux rtol = 1e-8
    end

    @testset "a carried point retreats on the grid it failed on" begin
        # the whole drive converges on the coarsest grid at once, and the
        # point carried to the middle grid fails there: the schedule
        # retreats and climbs back on the middle grid before growing past
        # it, and reaches the direct solve's point
        c, d = testjpacircuitnumeric()
        pump = (2*pi*4.75001e9,)
        drive = [(mode = (1,), port = 1, current = 2e-9)]
        controlled = recovery_solver(failat = [2])
        r = hbnlsolve(pump, (8,), drive, c, d; keyedarrays = false,
            method = Staged(inner = controlled.method,
                grids = [(2,), (4,), (8,)], s0 = 1.0, smin = 0.1))
        @test [st.grid for st in r.solverinfo.stages] ==
            [(2,), (4,), (4,), (4,), (8,)]
        @test r.solverinfo.converged
        fresh = hbnlsolve(pump, (8,), drive, c, d; method = Newton(),
            keyedarrays = false, atol = 1e-12)
        @test r.nodeflux ≈ fresh.nodeflux rtol = 1e-8
        # every solve counts against `maxattempts`, the retreat and the
        # evaluation of the point returned included: four are the coarse
        # solve, the failed growth, one retreat and that evaluation, and the
        # walk stops before climbing back
        budget = recovery_solver(failat = [2])
        rb = @test_logs (:warn,) match_mode = :any hbnlsolve(pump, (8,),
            drive, c, d; keyedarrays = false, method = Staged(
                inner = budget.method, grids = [(2,), (4,), (8,)], s0 = 1.0,
                smin = 0.1, maxattempts = 4))
        @test length(budget.starts) == 4
        @test length(rb.solverinfo.stages) == 4
        @test !rb.solverinfo.converged
        @test rb.solverinfo.stages[end].action === :evaluate
        # a first grid which converges at no drive is reported as that, and
        # not as a fold at zero drive
        stuck = recovery_solver(failat = collect(1:20))
        q = @test_logs (:warn,) hbnlsolve(pump, (4,), drive, c, d;
            keyedarrays = false, method = Staged(inner = stuck.method,
                grids = [(4,)], s0 = 0.5, smin = 0.1))
        @test !q.solverinfo.converged
        @test isnan(q.solverinfo.sourcefold)
    end

    @testset "a netlist source is scaled with the drive and checked once" begin
        # a junction biased near its critical current by a current source
        # of the netlist, and pumped weakly: the continuation scales the
        # bias with the pump, as `setdrive!` scales a problem's drive, so a
        # stage at a fraction of the drive is the direct solve at that
        # fraction of both; and
        # the warning of a junction near its critical current comes once,
        # from the outcome, and not from every stage
        Lj = 1000e-12
        Ic = LjtoIc(Lj)
        biased(Ib) = Circuit([:p1 => Port(1; Z0 = 50.0), :l => Inductor(1e-9),
            :jj => JosephsonJunction(Lj), :c2 => Capacitor(1000e-15),
            :ib => CurrentSource(Ib)],
            [((:p1, 1), (:l, 1)), ((:l, 2), (:jj, 1), (:c2, 1), (:ib, 2)),
             ((:jj, 2), (:c2, 2), (:p1, 2), (:ib, 1), Ground)])
        wb = (2*pi*4.0e9,)
        pump(I) = [(mode = (1,), port = 1, current = I)]
        kw = (; dc = true, odd = true, even = true, keyedarrays = false)
        r = @test_logs (:warn,) hbnlsolve(wb, (8,), pump(1e-9),
            biased(0.995*Ic); method = Staged(s0 = 0.25), kw...)
        @test r.solverinfo.converged
        # one stage at a quarter of the drive, after which two solves leave
        # none short of full drive beside the evaluation of its point, which
        # the schedule returns
        quarter = @test_logs (:warn,) hbnlsolve(wb, (8,), pump(1e-9),
            biased(0.995*Ic); method = Staged(grids = [(8,)], s0 = 0.25,
                interioratol = 1e-12, maxattempts = 2), kw...)
        ref = hbnlsolve(wb, (8,), pump(0.25e-9), biased(0.25*0.995*Ic);
            method = Newton(), atol = 1e-12, kw...)
        @test quarter.nodeflux ≈ ref.nodeflux rtol = 1e-8
    end

    @testset "guards" begin
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = Staged(grids = [(2,2), (4,2)]))
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = Staged(inner = Staged()))
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = Staged(s0 = 0.0))
        # the continuation warm starts its own stages and solves a system of
        # its own at each, so a start, and the system of one solve handed
        # back, are refused
        for kw in ((; x0 = zeros(ComplexF64, 4)), (; returnsystem = true),
                (; debugJacobian = true))
            @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src,
                circuit, defs; dc = true, odd = true, even = true,
                method = Staged(), kw...)
        end
    end

    @testset "a spent schedule returns not converged" begin
        # two solves are one stage on the coarsest grid and the evaluation
        # of its point, which the stage keeps in reserve
        r = @test_logs (:warn,) match_mode=:any hbnlsolve(
            (w1,w2), (8,4), src, circuit, defs; dc = true, odd = true,
            even = true, method = Staged(maxattempts = 2))
        @test !r.solverinfo.converged
        st = r.solverinfo.stages
        @test length(st) == 2
        @test isnan(r.solverinfo.sourcefold)
        # what it returns is at the modes asked for, not at the coarse grid
        # of its last stage, with the operating point which a sensitivity
        # through it reads
        full = JosephsonCircuits.pumpmodeset((w1,w2), (8,4), (16,8);
            dc = true, odd = true, even = true)[1].modes
        @test r.modes == full
        @test st[1].grid != (8,4)
        @test st[end].action === :evaluate && st[end].grid == (8,4)
        o = hbnlsolve((w1,w2), (8,4), src, circuit, defs; dc = true,
            odd = true, even = true, method = Staged(maxattempts = 2),
            returnoperatingpoint = true, warnnotconverged = false)
        @test !o.solverinfo.converged
        @test o.modes == full
        @test o.operatingpoint.Nmodes == length(full)
        # the flag which silences a solve silences the schedule too; with
        # one solve no stage fits beside the evaluation, which takes the
        # cold start
        q = @test_logs min_level=Base.CoreLogging.Warn hbnlsolve(
            (w1,w2), (8,4), src, circuit, defs; dc = true, odd = true,
            even = true, method = Staged(maxattempts = 1), warnnotconverged = false)
        @test !q.solverinfo.converged
        @test only(q.solverinfo.stages).action === :evaluate
    end

end
