using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "deprecated" begin

    @testset "connectS deprecation warnings" begin
        Sa = rand(Complex{Float64},3,3)
        Sb = rand(Complex{Float64},3,3)
        Sout1 = zeros(Complex{Float64},1,1)
        Sout2 = zeros(Complex{Float64},4,4)
        @test_logs((:warn, lazy"connectS(Sa::AbstractArray, k::Int, l::Int)` is deprecated, use `intraconnectS(Sa, k, l)` instead."),JosephsonCircuits.connectS(Sa,1,2));
        @test_logs((:warn, lazy"connectS(Sa::AbstractArray, Sb::AbstractArray, k::Int, l::Int)` is deprecated, use `interconnectS(Sa, Sb, k, l)` instead."),JosephsonCircuits.connectS(Sa,Sb,1,2));
        @test_logs((:warn, lazy"connectS!(Sout, Sa, k::Int, l::Int)` is deprecated, use `intraconnectS!(Sout, Sa, k, l)` instead."),JosephsonCircuits.connectS!(Sout1,Sa,1,2));
        @test_logs((:warn, lazy"connectS!(Sout, Sa, Sb, k::Int, l::Int)` is deprecated, use `interconnectS!(Sout, Sa, Sb, k, l)` instead."),JosephsonCircuits.connectS!(Sout2,Sa,Sb,1,2));
    end

    @testset "hbsolve deprecation warnings" begin
        # define the circuit components
        circuit = Array{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}},1}(undef,0)

        # port on the left side
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)) 
        push!(circuit,("Lj1","2","0",:Lj)) 
        push!(circuit,("C2","2","0",:Cj))

        circuitdefs = Dict{Symbol,Complex{Float64}}(
            :Lj =>1000.0e-12,
            :Cc => 100.0e-15,
            :Cj => 1000.0e-15,
            :Rleft => 50.0,
        )
        @test_logs((:warn,"""Calls the new harmonic balance solvers, [`hbnlsolve`](@ref) and\n[`hblinsolve`](@ref), which work for an arbitrary number of modes and ports),\nusing an identical syntax to the legacy harmonic balance solver, which only\nsupported four wave mixing processes involving single strong tone and an\narbitrary number of tone in the linearized solver. This function is\nprimarily for testing the new solvers and is now deprecated. Please switch\nto the new syntax.\n    """),
            hbsolve(2*pi*(4.5:0.5:5.0)*1e9,2*pi*4.75001*1e9,0.00565e-6,2,2,circuit,circuitdefs,pumpports=[1]),
        )
    end

    @testset "hbnlsolve switchofflinesearchtol and alphamin kwarg deprecation" begin
        # define the circuit components
        circuit = Array{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}},1}(undef,0)

        # port on the left side
        push!(circuit,("P1","1","0",1))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc)) 
        push!(circuit,("Lj1","2","0",:Lj)) 
        push!(circuit,("C2","2","0",:Cj))

        circuitdefs = Dict{Symbol,Complex{Float64}}(
            :Lj =>1000.0e-12,
            :Cc => 100.0e-15,
            :Cj => 1000.0e-15,
            :Rleft => 50.0,
        )

        ws = 2*pi*(4.5:0.5:5.0)*1e9
        wp = (2*pi*4.75001*1e9,)
        sources = [(mode=(1,),port=1,current=0.00565e-6)]
        Nmodulationharmonics = (2,)
        Npumpharmonics = (4,)

        @test_logs((:warn,lazy"The `switchofflinesearchtol` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions."),
            JosephsonCircuits.hbnlsolve(wp, Npumpharmonics, sources, circuit,
                circuitdefs;switchofflinesearchtol = 1),
        )
        @test_logs((:warn,lazy"The `alphamin` kwarg is deprecated and no longer used (and no longer necessary). Please remove it to avoid errors in future versions."),
            JosephsonCircuits.hbnlsolve(wp, Npumpharmonics, sources, circuit, circuitdefs;alphamin = 0.1),
        )
        @test_logs((:warn,lazy"The `maxharmonics` kwarg is deprecated and no longer used. `Nharmonics` is the retained set of modes and `Nevaluationharmonics` the grid on which the nonlinearity is sampled. Please remove it to avoid errors in future versions."),
            JosephsonCircuits.hbnlsolve(wp, Npumpharmonics, sources, circuit, circuitdefs;maxharmonics = (2,)),
        )
        @test_logs((:warn,lazy"The `maxpumpharmonics` kwarg is deprecated and no longer used. `Npumpharmonics` is the retained set of pump modes and `Nevaluationharmonics` the grid on which the nonlinearity is sampled. Please remove it to avoid errors in future versions."),
            JosephsonCircuits.hbsolve(ws, wp, sources, Nmodulationharmonics, Npumpharmonics, circuit, circuitdefs;maxpumpharmonics = (2,)),
        )

        @test_logs((:warn,lazy"The `returnZ`, `returnZadjoint`, `returnZsensitivity`, and `returnZsensitivityadjoint` kwargs have been removed. Please compute them from scattering parameters matrices."),
            JosephsonCircuits.hblinsolve(ws, circuit, circuitdefs;returnZ = true),
        )
        @test_logs((:warn,lazy"The `returnZ`, `returnZadjoint`, `returnZsensitivity`, and `returnZsensitivityadjoint` kwargs have been removed. Please compute them from scattering parameters matrices."),
            JosephsonCircuits.hblinsolve(ws, circuit, circuitdefs;returnZadjoint = true),
        )
        @test_logs((:warn,lazy"The `returnZ`, `returnZadjoint`, `returnZsensitivity`, and `returnZsensitivityadjoint` kwargs have been removed. Please compute them from scattering parameters matrices."),
            JosephsonCircuits.hblinsolve(ws, circuit, circuitdefs;returnZsensitivity = true),
        )
        @test_logs((:warn,lazy"The `returnZ`, `returnZadjoint`, `returnZsensitivity`, and `returnZsensitivityadjoint` kwargs have been removed. Please compute them from scattering parameters matrices."),
            JosephsonCircuits.hblinsolve(ws, circuit, circuitdefs;returnZsensitivityadjoint = true),
        )
    end

end