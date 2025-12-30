using JosephsonCircuits
using LinearAlgebra
using Test

@testset verbose=true "deprecated" begin

    @testset "connectS deprecation warnings" begin
        Sa = rand(Complex{Float64},3,3)
        Sb = rand(Complex{Float64},3,3)
        Sout1 = zeros(Complex{Float64},1,1)
        Sout2 = zeros(Complex{Float64},4,4)
        @test_logs((:warn, "connectS(Sa::AbstractArray, k::Int, l::Int)` is deprecated, use `intraconnectS(Sa, k, l)` instead."),JosephsonCircuits.connectS(Sa,1,2));
        @test_logs((:warn, "connectS(Sa::AbstractArray, Sb::AbstractArray, k::Int, l::Int)` is deprecated, use `interconnectS(Sa, Sb, k, l)` instead."),JosephsonCircuits.connectS(Sa,Sb,1,2));
        @test_logs((:warn, "connectS!(Sout, Sa, k::Int, l::Int)` is deprecated, use `intraconnectS!(Sout, Sa, k, l)` instead."),JosephsonCircuits.connectS!(Sout1,Sa,1,2));
        @test_logs((:warn, "connectS!(Sout, Sa, Sb, k::Int, l::Int)` is deprecated, use `interconnectS!(Sout, Sa, Sb, k, l)` instead."),JosephsonCircuits.connectS!(Sout2,Sa,Sb,1,2));
    end

end