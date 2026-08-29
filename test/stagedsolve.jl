isdefined(Main, :testjpacircuit) || include(joinpath(@__DIR__, "testcircuits.jl"))
using JosephsonCircuits, Test, LinearAlgebra

@testset "method = :staged" begin
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

    @testset "matches the direct solve" begin
        ra = hbnlsolve((w1,w2), (8,4), src, circuit, defs;
            dc = true, odd = true, even = true, method = :newtonkrylov)
        rb = hbnlsolve((w1,w2), (8,4), src, circuit, defs;
            dc = true, odd = true, even = true, method = :staged)
        @test rb.solverinfo.converged
        a = vec(Array(ra.S)); b = vec(Array(rb.S))
        @test norm(a - b)/norm(a) < 1e-6
    end

    @testset "hbsolve integration" begin
        # offset from the 5 GHz pump so no signal + pump mode lands at
        # zero total frequency
        ws = 2*pi*(4.55:0.3:5.5)*1e9
        ra = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs;
            dc = true, threewavemixing = true, fourwavemixing = true)
        rb = hbsolve(ws, (w1,w2), src, (2,2), (8,4), circuit, defs;
            dc = true, threewavemixing = true, fourwavemixing = true,
            method = :staged)
        a = vec(Array(ra.linearized.S)); b = vec(Array(rb.linearized.S))
        @test norm(a - b)/norm(a) < 1e-6
    end

    @testset "stage records in solverinfo" begin
        r = hbnlsolve((w1,w2), (8,4), src, circuit, defs;
            dc = true, odd = true, even = true, method = :staged)
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

    @testset "guards" begin
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = :staged,
            stagedkwargs = (; grids = [(2,2), (4,2)]))
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = :staged,
            stagedkwargs = (; innermethod = :staged))
        @test_throws ArgumentError hbnlsolve((w1,w2), (8,4), src, circuit,
            defs; dc = true, odd = true, even = true, method = :staged,
            stagedkwargs = (; s0 = 0.0))
    end
end
