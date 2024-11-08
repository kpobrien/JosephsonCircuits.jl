using JosephsonCircuits
using LinearAlgebra
using Test
import StaticArrays

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

    @testset "connectS consistency StaticArrays" begin
        Sx = rand(Complex{Float64},3,3,3)
        Sy = rand(Complex{Float64},3,3,3)
        Sboth = zeros(Complex{Float64},6,6,3)
        Sboth[1:size(Sx,1),1:size(Sx,1),:,:] .= Sx
        Sboth[size(Sx,1)+1:end,size(Sx,1)+1:end,:,:] .= Sy
        Sboth = [StaticArrays.MMatrix{size(Sboth,1),size(Sboth,2)}(Sboth[:,:,i]) for i=1:size(Sboth,3)]
        Sx = [StaticArrays.MMatrix{size(Sx,1),size(Sx,2)}(Sx[:,:,i]) for i=1:size(Sx,3)]
        Sy = [StaticArrays.MMatrix{size(Sy,1),size(Sy,2)}(Sy[:,:,i]) for i=1:size(Sy,3)]
        Sout1 = JosephsonCircuits.connectS.(Sboth,1,2+size(Sx,1))
        Sout2 = JosephsonCircuits.connectS.(Sx,Sy,1,2)
        @test isapprox(Sout1,Sout2)
    end

    @testset "connectS consistency StaticArrays" begin
        Sx = rand(Complex{Float64},3,3)
        Sy = rand(Complex{Float64},3,3)
        Sboth = zeros(Complex{Float64},6,6)
        Sboth[1:size(Sx,1),1:size(Sx,1),:,:] .= Sx
        Sboth[size(Sx,1)+1:end,size(Sx,1)+1:end,:,:] .= Sy
        # convert to MMatrix
        Sboth = StaticArrays.MMatrix{size(Sboth,1),size(Sboth,2)}(Sboth)
        Sx = StaticArrays.MMatrix{size(Sx,1),size(Sx,2)}(Sx)
        Sy = StaticArrays.MMatrix{size(Sy,1),size(Sy,2)}(Sy)
        Sout1 = JosephsonCircuits.connectS(Sboth,2,1+size(Sx,1))
        Sout2 = JosephsonCircuits.connectS(Sx,Sy,2,1)
        @test isapprox(Sout1,Sout2)
    end


    @testset "StoZ, StoY, StoA, StoB, StoABCD consistency" begin
        # the different functions we want to test
        for f in [
                (JosephsonCircuits.ZtoS,JosephsonCircuits.StoZ),
                (JosephsonCircuits.YtoS,JosephsonCircuits.StoY),
                (JosephsonCircuits.AtoS,JosephsonCircuits.StoA),
                (JosephsonCircuits.BtoS,JosephsonCircuits.StoB),
                (JosephsonCircuits.ABCDtoS,JosephsonCircuits.StoABCD),
            ]
            # single matrix input
            for portimpedances in [
                    rand(Complex{Float64}), rand(Complex{Float64},2),
                    (StaticArrays.@MVector rand(Complex{Float64},2))
                ]
                for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
                    arg2 = f[1](arg1,portimpedances=portimpedances)
                    arg3 = f[2](arg2,portimpedances=portimpedances)
                    @test isapprox(arg1,arg3)
                end
            end
            # array input
            for portimpedances in [rand(Complex{Float64}), rand(Complex{Float64},2,10)]
                for arg1 in [rand(Complex{Float64},2,2,10)]
                    arg2 = f[1](arg1,portimpedances=portimpedances)
                    arg3 = f[2](arg2,portimpedances=portimpedances)
                    @test isapprox(arg1,arg3)
                end
            end
            # vector of matrices
            for portimpedances in [rand(Complex{Float64}), rand(Complex{Float64},2), (StaticArrays.@MVector rand(Complex{Float64},2))]
                for arg1 in [
                        [rand(Complex{Float64},2,2) for i in 1:10],
                        [(StaticArrays.@MMatrix rand(Complex{Float64},2,2)) for i in 1:10],
                    ]
                    arg2 = [f[1](arg1[i],portimpedances=portimpedances) for i in 1:10]
                    arg3 = [f[2](arg2[i],portimpedances=portimpedances) for i in 1:10]
                    @test isapprox(arg1,arg3)
                end
            end
        end
    end

    @testset "StoT, AtoB, ZtoA, YtoA, YtoB, ZtoB, ZtoY consistency" begin
        # the different functions we want to test
        for f in [
                (JosephsonCircuits.StoT,JosephsonCircuits.TtoS),
                (JosephsonCircuits.AtoB,JosephsonCircuits.BtoA),
                (JosephsonCircuits.ZtoA,JosephsonCircuits.AtoZ),
                (JosephsonCircuits.YtoA,JosephsonCircuits.AtoY),
                (JosephsonCircuits.YtoB,JosephsonCircuits.BtoY),
                (JosephsonCircuits.ZtoB,JosephsonCircuits.BtoZ),
                (JosephsonCircuits.ZtoY,JosephsonCircuits.YtoZ),
            ]
            # single matrix input
            for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
                arg2 = f[1](arg1)
                arg3 = f[2](arg2)
                @test isapprox(arg1,arg3)
            end
            # array input
            for arg1 in [rand(Complex{Float64},2,2,10)]
                arg2 = f[1](arg1)
                arg3 = f[2](arg2)
                @test isapprox(arg1,arg3)
            end
            # vector of matrices
            for arg1 in [
                    [rand(Complex{Float64},2,2) for i in 1:10],
                    [(StaticArrays.@MMatrix rand(Complex{Float64},2,2)) for i in 1:10],
                ]
                arg2 = [f[1](arg1[i]) for i in 1:10]
                arg3 = [f[2](arg2[i]) for i in 1:10]
                @test isapprox(arg1,arg3)
            end
        end
    end

    @testset "ZC_basis_coupled_lines" begin
        Cg = 1.6670474181399462e-10
        Cm = 9.320861870486729e-12
        Ls = 4.167618545349866e-7
        Lm = 2.330215467621682e-8
        C = JosephsonCircuits.LinearAlgebra.Symmetric([Cg -Cm;-Cm Cg])
        L = JosephsonCircuits.LinearAlgebra.Symmetric([Ls Lm;Lm Ls])
        ZC, TI, TV, theta, U, lambda, S = JosephsonCircuits.ZC_basis_coupled_tlines(L,C)

        @test isapprox(U'*C*U,theta^2)
        @test isapprox(S'*(theta*U'*L*U*theta)*S,lambda^2)
        @test isapprox(inv(TI)*C*L*TI,lambda^2)
        @test isapprox(TV'*TI,I)
    end

    @testset "Z_canonical_coupled_line_circuits" begin
        Zeven = 52.0
        Zodd = 48.0
        neven = 2.2
        nodd = 2.1
        l = 2*3.0e-3
        c = 2.998e8
        omega = 2*pi*5e9

        function S_canonical_coupled_line_circuits(i::Int, thetae, thetao, Ze, Zo)

            # define an open
            Sopen = ones(Complex{Float64},1,1)

            # and a short
            Sshort = -ones(Complex{Float64},1,1)
            
            # and a coupled transmission line
            S = JosephsonCircuits.ZtoS(JosephsonCircuits.Z_coupled_tline(thetae,thetao,Zeven,Zodd));

            if i == 1
                S = JosephsonCircuits.connectS(S,Sopen,3,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 2
                S = JosephsonCircuits.connectS(S,Sshort,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 3
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sopen,1,1)
            elseif i == 4
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,3,1)
            elseif i == 5
                S = JosephsonCircuits.connectS(S,Sopen,3,1)
                S = JosephsonCircuits.connectS(S,Sopen,1,1)
            elseif i == 6
                S = JosephsonCircuits.connectS(S,Sshort,3,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 7
                S = JosephsonCircuits.connectS(S,4,3)
            elseif i == 8
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,1,1)
            elseif i == 9
                S = JosephsonCircuits.connectS(S,Sshort,4,1)
                S = JosephsonCircuits.connectS(S,Sshort,3,1)
            elseif i == 10
                S = JosephsonCircuits.connectS(S,Sopen,4,1)
                S = JosephsonCircuits.connectS(S,Sopen,3,1)        
            else
                throw(ArgumentError("Canonical coupled line circuit number must be 1-10."))
            end
                
            return S
        end

        for i in 1:10
            S1 = S_canonical_coupled_line_circuits(i, neven*omega/c*l,nodd*omega/c*l,Zeven,Zodd)
            S2 = JosephsonCircuits.ZtoS(JosephsonCircuits.Z_canonical_coupled_line_circuits(i, neven*omega/c*l, nodd*omega/c*l, Zeven, Zodd));
            @test isapprox(S1,S2)
        end

    end

end