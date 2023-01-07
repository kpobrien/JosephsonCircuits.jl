using JosephsonCircuits
using Test

@testset verbose=true "qesparams" begin

    @testset "calcimpedance" begin
        @test_throws(
            ErrorException("Unknown component type"),
            JosephsonCircuits.calcimpedance(30.0,:D,-1.0,nothing))
    end

    @testset "calcimpedance" begin
        @variables w
        @test_throws(
            ErrorException("Unknown component type"),
            JosephsonCircuits.calcimpedance(30*w,:D,-2.0,w))
    end

    @testset "calccm!" begin
        cm=Float64[0,0]
        @test_throws(
            DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.calccm!(cm,[3/5 4/5;4/5 3/5],[-1,1,2]))
    end

    @testset "calccm!" begin
        cm=Float64[0,0]
        @test_throws(
            DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
            JosephsonCircuits.calccm!(cm,[3/5 4/5;4/5 3/5;0 0;0 0],[-1,1]))
    end

    @testset "calccm!" begin
        @variables a b
        cm=Num[0,0]
        @test_throws(
            DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.calccm!(cm,[a b; b a],[-1,1,2]))
    end

    @testset "calccm!" begin
         @variables a b
         cm=Num[0,0]
        @test_throws(
            DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
            JosephsonCircuits.calccm!(cm,[a b; b a; 0 0; 0 0],[-1,1]))
    end

    @testset "calccm!" begin
        cm=Float64[0, 0]
        @test_throws(
            DimensionMismatch("Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2 3;5 6 7],[-1,1]))
    end

    @testset "calccm!" begin
        cm=Float64[0, 0]
        @test_throws(
            DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."),
            JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2; 3 4; 5 6; 7 8],[-1,1]))
    end

    @testset "calccm!" begin
        cm=Float64[0, 0]
        @test_throws(
            DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2 3 4;5 6 7 8],[-1,1,2]))
    end

    @testset "calccm!" begin
        cm=Float64[0, 0, 0]
        @test_throws(
            DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
            JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2 3 4;5 6 7 8],[-1,1]))
    end

    @testset "calccm!" begin
         @variables a b c d an bn cn dn
         cm = Num[0, 0]
        @test_throws(
            DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn],[1, -1, 2]))
    end

    @testset "calccm!" begin
        @variables a b c d an bn cn dn;cm = Num[0, 0]
        @test_throws(
            DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."),
            JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn; 0 0; 0 0],[1, -1]))
    end

    @testset "calccm!" begin
        @variables a b c d an bn cn dn
        cm = Num[0, 0, 0]
        @test_throws(
            DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
            JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn],[1, -1]))
    end

    @testset "calccm!" begin
        @variables a b c d an bn cn dn
        cm = Num[0, 0]
        @test_throws(
            DimensionMismatch("Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn 0; cn dn 0],[1, -1]))
    end

    @testset "calcqe!" begin
        @test_throws(
            DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcqe!([1 2;3 4],[1 2 3;4 5 6]))
    end

    @testset "calcqe!" begin
        @test_throws(
            DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcqe!([1 2;3 4],[1 2 3;4 5 6],[1 2;3 4]))
    end

    @testset "calcqe!" begin
        @test_throws(
            DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."),
            JosephsonCircuits.calcqe!(Float64[1 2;3 4],[1 2;3 4],[1 2;3 4;5 6]))
    end

    @testset "calcqeideal!" begin
        @test_throws(
            DimensionMismatch("Sizes of QE and S matrices must be equal."),
            JosephsonCircuits.calcqeideal!([1 2;3 4],[1 2 3;4 5 6]))
    end
end