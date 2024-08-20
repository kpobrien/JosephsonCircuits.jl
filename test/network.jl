using JosephsonCircuits
using Test

@testset verbose=true "network" begin

    # one network
    @testset "connectS" begin
        Sx = Float64[1 2;3 4]
        @test_throws(
            ArgumentError("Port `k` is smaller than one."),
            JosephsonCircuits.connectS(Sx,0,1)
        )

        @test_throws(
            ArgumentError("Port `l` is smaller than one."),
            JosephsonCircuits.connectS(Sx,1,0)
        )

        @test_throws(
            ArgumentError("Port `k` is larger than number of ports in `Sx`."),
            JosephsonCircuits.connectS(Sx,3,1)
        )

        @test_throws(
            ArgumentError("Port `l` is larger than number of ports in `Sx`."),
            JosephsonCircuits.connectS(Sx,1,3)
        )

        @test_throws(
            ArgumentError("`k` and `l` cannot be equal because a port cannot be merged with itself."),
            JosephsonCircuits.connectS(Sx,1,1)
        )
    end

    @testset "connectS" begin
        Sx = Float64[1 2 3;4 5 6]
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sx` must be equal."),
            JosephsonCircuits.connectS(Sx,1,2)
        )
    end

    # one network, in place
    @testset "connectS!" begin
        Sx = rand(Complex{Float64},3,3)
        Sout = zeros(Complex{Float64},2,2)
        @test_throws(
            DimensionMismatch("Length of first two dimensions must be 2 smaller for `Sout` than `Sx` because we are merging two ports."),
            JosephsonCircuits.connectS!(Sout,Sx,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},2,2,1)
        @test_throws(
            DimensionMismatch("`Sout` and `Sx` must have the same number of dimensions."),
            JosephsonCircuits.connectS!(Sout,Sx,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},4)
        Sout = zeros(Complex{Float64},2)
        @test_throws(
            DimensionMismatch("`Sout` and `Sx` must have at least two dimensions."),
            JosephsonCircuits.connectS!(Sout,Sx,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},2,3)
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sx,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},3,3,3)
        Sout = zeros(Complex{Float64},1,1,4)
        @test_throws(
            DimensionMismatch("Non-port axis lengths of `Sx` and `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sx,1,2)
        )
    end

    # two networks
    @testset "connectS" begin
        Sx = Float64[1 2;3 4]
        Sy = Float64[5 6;7 8]
        @test_throws(
            ArgumentError("Port `k` is smaller than one."),
            JosephsonCircuits.connectS(Sx,Sy,0,1)
        )

        @test_throws(
            ArgumentError("Port `l` is smaller than one."),
            JosephsonCircuits.connectS(Sx,Sy,1,0)
        )

        @test_throws(
            ArgumentError("Port `k` is larger than number of ports in `Sx`."),
            JosephsonCircuits.connectS(Sx,Sy,3,1)
        )

        @test_throws(
            ArgumentError("Port `l` is larger than number of ports in `Sy`."),
            JosephsonCircuits.connectS(Sx,Sy,1,3)
        )

    end

    @testset "connectS" begin
        Sx = Float64[1 2 3;4 5 6]
        Sy = Float64[5 6;7 8]
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sx` must be equal."),
            JosephsonCircuits.connectS(Sx,Sy,1,2)
        )
    end

    @testset "connectS" begin
        Sx = Float64[1 2;4 5]
        Sy = Float64[5 6 7;8 9 10]
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sy` must be equal."),
            JosephsonCircuits.connectS(Sx,Sy,1,2)
        )
    end

    # two networks, in place
    @testset "connectS!" begin
        Sx = rand(Complex{Float64},3,3)
        Sy = rand(Complex{Float64},3,3)
        Sout = zeros(Complex{Float64},2,2)
        @test_throws(
            DimensionMismatch("First two dimensions of `Sout` must be `m+n-2`."),
            JosephsonCircuits.connectS!(Sout,Sx,Sy,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},3,3)
        Sy = rand(Complex{Float64},3,3,1)
        Sout = zeros(Complex{Float64},4,4)
        @test_throws(
            DimensionMismatch("`Sx` and `Sy` must have the same number of dimensions."),
            JosephsonCircuits.connectS!(Sout,Sx,Sy,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},4,4)
        Sy = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},6,6,1)
        @test_throws(
            DimensionMismatch("`Sout`, `Sx`, and `Sy` must have the same number of dimensions."),
            JosephsonCircuits.connectS!(Sout,Sx,Sy,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},4)
        Sy = rand(Complex{Float64},4)
        Sout = zeros(Complex{Float64},6)
        @test_throws(
            DimensionMismatch("`Sout`, `Sx`, and `Sy` must have atleast two dimensions."),
            JosephsonCircuits.connectS!(Sout,Sx,Sy,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},4,4)
        Sy = rand(Complex{Float64},4,4)
        Sout = zeros(Complex{Float64},6,7)
        @test_throws(
            DimensionMismatch("Lengths of first two dimensions of `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sx,Sy,1,2)
        )
    end

    @testset "connectS!" begin
        Sx = rand(Complex{Float64},3,3,3)
        Sy = rand(Complex{Float64},3,3,3)
        Sout = zeros(Complex{Float64},4,4,4)
        @test_throws(
            DimensionMismatch("Non-port axis lengths of `Sx`, `Sy, and `Sout` must be equal."),
            JosephsonCircuits.connectS!(Sout,Sx,Sy,1,2)
        )
    end

    # check for consistency between the one and two network functions
    @testset "connectS consistency" begin
        Sx = rand(Complex{Float64},3,3,3)
        Sy = rand(Complex{Float64},3,3,3)
        Sboth = zeros(Complex{Float64},6,6,3)
        Sboth[1:size(Sx,1),1:size(Sx,1),:,:] .= Sx
        Sboth[size(Sx,1)+1:end,size(Sx,1)+1:end,:,:] .= Sy
        Sout1 = JosephsonCircuits.connectS(Sboth,1,2+size(Sx,1))
        Sout2 = JosephsonCircuits.connectS(Sx,Sy,1,2)
        @test isapprox(Sout1,Sout2)
    end

    @testset "connectS consistency" begin
        Sx = rand(Complex{Float64},3,3)
        Sy = rand(Complex{Float64},3,3)
        Sboth = zeros(Complex{Float64},6,6)
        Sboth[1:size(Sx,1),1:size(Sx,1),:,:] .= Sx
        Sboth[size(Sx,1)+1:end,size(Sx,1)+1:end,:,:] .= Sy
        Sout1 = JosephsonCircuits.connectS(Sboth,2,1+size(Sx,1))
        Sout2 = JosephsonCircuits.connectS(Sx,Sy,2,1)
        @test isapprox(Sout1,Sout2)
    end

    @testset "StoT TtoS consistency" begin
        S0 = rand(Complex{Float64},4,4,1);
        S1 = JosephsonCircuits.TtoS(JosephsonCircuits.StoT(S0))
        @test isapprox(S0,S1)
    end

    @testset "StoZ ZtoS consistency" begin
        Z0 = rand(Complex{Float64},4,4,1);
        portimpedances = rand(Complex{Float64},4,1)
        S0 = JosephsonCircuits.ZtoS(Z0,portimpedances=portimpedances)
        Z1 = JosephsonCircuits.StoZ(S0,portimpedances=portimpedances)
        @test isapprox(Z0,Z1)
    end

    @testset "StoZ ZtoS consistency" begin
        Z0 = rand(Complex{Float64},4,4,1);
        portimpedances = rand(Complex{Float64},4,1)
        S0 = JosephsonCircuits.ZtoS(Z0,portimpedances=portimpedances)
        Z1 = JosephsonCircuits.StoZ(S0,portimpedances=portimpedances)
        @test isapprox(Z0,Z1)
    end

    @testset "StoZ ZtoS consistency" begin
        Z0 = rand(Complex{Float64},4,4);
        portimpedances = rand(Complex{Float64},4)
        S0 = JosephsonCircuits.ZtoS(Z0,portimpedances=portimpedances)
        Z1 = JosephsonCircuits.StoZ(S0,portimpedances=portimpedances)
        @test isapprox(Z0,Z1)
    end

    @testset "StoZ ZtoS consistency" begin
        Z0 = rand(Complex{Float64},2,2,100);
        portimpedances = rand(Complex{Float64},2,100)
        S0 = JosephsonCircuits.ZtoS(Z0,portimpedances=portimpedances)
        Z1 = JosephsonCircuits.StoZ(S0,portimpedances=portimpedances)
        @test isapprox(Z0,Z1)
    end

    @testset "StoZ ZtoS consistency" begin
        Z0 = rand(Complex{Float64},2,2,100);
        portimpedances = 50.0
        S0 = JosephsonCircuits.ZtoS(Z0,portimpedances=portimpedances)
        Z1 = JosephsonCircuits.StoZ(S0,portimpedances=portimpedances)
        @test isapprox(Z0,Z1)
    end

    @testset "StoY YtoS consistency" begin
        Y0 = rand(Complex{Float64},4,4,1);
        portimpedances = rand(Complex{Float64},4,1)
        S0 = JosephsonCircuits.YtoS(Y0,portimpedances=portimpedances)
        Y1 = JosephsonCircuits.StoY(S0,portimpedances=portimpedances)
        @test isapprox(Y0,Y1)
    end

    @testset "StoY YtoS consistency" begin
        Y0 = rand(Complex{Float64},4,4);
        portimpedances = rand(Complex{Float64},4)
        S0 = JosephsonCircuits.YtoS(Y0,portimpedances=portimpedances)
        Y1 = JosephsonCircuits.StoY(S0,portimpedances=portimpedances)
        @test isapprox(Y0,Y1)
    end

    @testset "StoY YtoS consistency" begin
        Y0 = rand(Complex{Float64},2,2,100);
        portimpedances = rand(Complex{Float64},2,100)
        S0 = JosephsonCircuits.YtoS(Y0,portimpedances=portimpedances)
        Y1 = JosephsonCircuits.StoY(S0,portimpedances=portimpedances)
        @test isapprox(Y0,Y1)
    end

    @testset "StoY YtoS consistency" begin
        Y0 = rand(Complex{Float64},2,2,100);
        portimpedances = 50.0
        S0 = JosephsonCircuits.YtoS(Y0,portimpedances=portimpedances)
        Y1 = JosephsonCircuits.StoY(S0,portimpedances=portimpedances)
        @test isapprox(Y0,Y1)
    end

    @testset "StoA AtoS consistency" begin
        S0 = rand(Complex{Float64},4,4);
        portimpedances = rand(Complex{Float64})
        A0 = JosephsonCircuits.StoA(S0,portimpedances=portimpedances)
        S1 = JosephsonCircuits.AtoS(A0,portimpedances=portimpedances)
        @test isapprox(S0,S1)
    end

    @testset "StoA AtoS consistency" begin
        S0 = rand(Complex{Float64},4,4,10);
        portimpedances = rand(Complex{Float64},4,10)
        A0 = JosephsonCircuits.StoA(S0,portimpedances=portimpedances)
        S1 = JosephsonCircuits.AtoS(A0,portimpedances=portimpedances)
        @test isapprox(S0,S1)
    end

    @testset "StoA AtoS consistency" begin
        S0 = rand(Complex{Float64},4,4,10);
        portimpedances = rand(Complex{Float64})
        A0 = JosephsonCircuits.StoA(S0,portimpedances=portimpedances)
        S1 = JosephsonCircuits.AtoS(A0,portimpedances=portimpedances)
        @test isapprox(S0,S1)
    end

    @testset "StoB BtoS consistency" begin
        S0 = rand(Complex{Float64},4,4);
        portimpedances = rand(Complex{Float64})
        B0 = JosephsonCircuits.StoA(S0,portimpedances=portimpedances)
        S1 = JosephsonCircuits.AtoS(B0,portimpedances=portimpedances)
        @test isapprox(S0,S1)
    end

    @testset "StoB BtoS consistency" begin
        S0 = rand(Complex{Float64},4,4,10);
        portimpedances = rand(Complex{Float64},4,10)
        B0 = JosephsonCircuits.StoA(S0,portimpedances=portimpedances)
        S1 = JosephsonCircuits.AtoS(B0,portimpedances=portimpedances)
        @test isapprox(S0,S1)
    end

    @testset "StoB BtoS consistency" begin
        S0 = rand(Complex{Float64},4,4,10);
        portimpedances = rand(Complex{Float64})
        B0 = JosephsonCircuits.StoA(S0,portimpedances=portimpedances)
        S1 = JosephsonCircuits.AtoS(B0,portimpedances=portimpedances)
        @test isapprox(S0,S1)
    end

end