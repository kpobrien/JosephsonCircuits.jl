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

    @testset "guards" begin
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = Staged(grids = [(2,2), (4,2)]))
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = Staged(inner = Staged()))
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = Staged(s0 = 0.0))
    end

    @testset "a spent schedule returns not converged" begin
        r = @test_logs (:warn, r"maxattempts") match_mode=:any hbnlsolve(
            (w1,w2), (8,4), src, circuit, defs; dc = true, odd = true,
            even = true, method = Staged(maxattempts = 1))
        @test !r.solverinfo.converged
        @test length(r.solverinfo.stages) == 1
        @test isnan(r.solverinfo.sourcefold)
        # the flag which silences a solve silences the schedule too
        q = @test_logs min_level=Base.CoreLogging.Warn hbnlsolve(
            (w1,w2), (8,4), src, circuit, defs; dc = true, odd = true,
            even = true, method = Staged(maxattempts = 1), warnnotconverged = false)
        @test !q.solverinfo.converged
    end

end
