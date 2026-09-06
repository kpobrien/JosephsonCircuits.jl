using JosephsonCircuits
using LinearAlgebra
using Test

# The deprecated forms warn and give the same numbers as the forms which
# replace them; the wording of the warnings is not pinned.
@testset verbose=true "deprecated" begin

    @testset "connectS is intraconnectS and interconnectS" begin
        Sa = rand(Complex{Float64},3,3)
        Sb = rand(Complex{Float64},3,3)
        Sout1 = zeros(Complex{Float64},1,1)
        Sout2 = zeros(Complex{Float64},4,4)
        a = @test_logs (:warn,) JosephsonCircuits.connectS(Sa,1,2)
        @test a == JosephsonCircuits.intraconnectS(Sa,1,2)
        b = @test_logs (:warn,) JosephsonCircuits.connectS(Sa,Sb,1,2)
        @test b == JosephsonCircuits.interconnectS(Sa,Sb,1,2)
        @test_logs (:warn,) JosephsonCircuits.connectS!(Sout1,Sa,1,2)
        @test Sout1 == JosephsonCircuits.intraconnectS(Sa,1,2)
        @test_logs (:warn,) JosephsonCircuits.connectS!(Sout2,Sa,Sb,1,2)
        @test Sout2 == JosephsonCircuits.interconnectS(Sa,Sb,1,2)
    end

    @testset "the deprecated solver forms and keywords" begin
        circuit = Array{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}},1}(undef,0)
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        circuitdefs = Dict{Symbol,Complex{Float64}}(
            :Lj =>1000.0e-12, :Cc => 100.0e-15, :Cj => 1000.0e-15,
            :Rleft => 50.0)
        ws = 2*pi*[4.5e9]
        wp = (2*pi*4.75001*1e9,)
        Ip = 0.00565e-6
        sources = [(mode=(1,),port=1,current=Ip)]
        ref = hbsolve(ws, wp, sources, (2,), (2,), circuit, circuitdefs;
            keyedarrays = false)
        same(sol) = isapprox(Array(sol.linearized.S), ref.linearized.S;
            rtol = 1e-10) && isapprox(Array(sol.nonlinear.nodeflux),
            ref.nonlinear.nodeflux; rtol = 1e-10)

        # the single pump frequency, integer harmonic count form of hbsolve:
        # its pump count Npumpmodes is the tuple form's (2*Npumpmodes,), so
        # the operating points agree; its signal mode set is the legacy
        # solver's own, so of the linearized outputs the signal to signal
        # entry is compared, with one signal mode either way
        old = @test_logs (:warn,) match_mode = :any hbsolve(ws, wp[1], Ip, 1, 2,
            circuit, circuitdefs, pumpports = [1], keyedarrays = true)
        oldref = hbsolve(ws, wp, sources, (1,), (4,), circuit, circuitdefs)
        @test isapprox(Array(old.nonlinear.nodeflux(outputmode = (1,))),
            Array(oldref.nonlinear.nodeflux(outputmode = (1,))); rtol = 1e-10)
        @test isapprox(old.linearized.S((0,), 1, (0,), 1, 1),
            oldref.linearized.S((0,), 1, (0,), 1, 1); rtol = 1e-6)

        nlref = hbnlsolve(wp, (2,), sources, circuit, circuitdefs;
            keyedarrays = false)
        for kw in ((switchofflinesearchtol = 1,), (alphamin = 0.1,),
                (maxharmonics = (2,),))
            sol = @test_logs (:warn,) match_mode = :any hbnlsolve(wp, (2,),
                sources, circuit, circuitdefs; keyedarrays = false, kw...)
            @test isapprox(sol.nodeflux, nlref.nodeflux; rtol = 1e-10)
        end
        sol = @test_logs (:warn,) match_mode = :any hbsolve(ws, wp, sources,
            (2,), (2,), circuit, circuitdefs; keyedarrays = false,
            maxpumpharmonics = (2,))
        @test same(sol)

        linref = hblinsolve(ws, circuit, circuitdefs; keyedarrays = false)
        for kw in ((returnZ = true,), (returnZadjoint = true,),
                (returnZsensitivity = true,), (returnZsensitivityadjoint = true,))
            lin = @test_logs (:warn,) match_mode = :any hblinsolve(ws, circuit,
                circuitdefs; keyedarrays = false, kw...)
            @test lin.S == linref.S
        end
    end

end
