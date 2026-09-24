using JosephsonCircuits
using LinearAlgebra
using Test
import StaticArrays

@testset verbose=true "network parameter conversion" begin

    @testset "StoZ, StoY, StoA, StoB, StoABCD consistency" begin
        # the different functions we want to test
        for f in [
                (JosephsonCircuits.ZtoS,JosephsonCircuits.StoZ,JosephsonCircuits.ZtoS!,JosephsonCircuits.StoZ!),
                (JosephsonCircuits.YtoS,JosephsonCircuits.StoY,JosephsonCircuits.YtoS!,JosephsonCircuits.StoY!),
                (JosephsonCircuits.AtoS,JosephsonCircuits.StoA,JosephsonCircuits.AtoS!,JosephsonCircuits.StoA!),
                (JosephsonCircuits.BtoS,JosephsonCircuits.StoB,JosephsonCircuits.BtoS!,JosephsonCircuits.StoB!),
                (JosephsonCircuits.ABCDtoS,JosephsonCircuits.StoABCD,JosephsonCircuits.ABCDtoS!,JosephsonCircuits.StoABCD!),
            ]
            # single matrix input
            for portimpedances in [
                    rand(Complex{Float64}), rand(Complex{Float64},2)
                ]
                for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
                    arg2 = f[1](arg1,portimpedances=portimpedances)
                    arg3 = f[2](arg2,portimpedances=portimpedances)
                    @test isapprox(arg1,arg3)
                    arg4 = copy(arg1)
                    @test isapprox(arg1,f[4](f[3](arg4,portimpedances=portimpedances),portimpedances=portimpedances))
                end
            end
            # array input
            for portimpedances in [rand(Complex{Float64}), rand(Complex{Float64},2,10)]
                for arg1 in [rand(Complex{Float64},2,2,10)]
                    arg2 = f[1](arg1,portimpedances=portimpedances)
                    arg3 = f[2](arg2,portimpedances=portimpedances)
                    @test isapprox(arg1,arg3)
                    arg4 = copy(arg1)
                    @test isapprox(arg1,f[4](f[3](arg4,portimpedances=portimpedances),portimpedances=portimpedances))
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
                (JosephsonCircuits.StoT,JosephsonCircuits.TtoS,JosephsonCircuits.StoT!,JosephsonCircuits.TtoS!),
                (JosephsonCircuits.AtoB,JosephsonCircuits.BtoA,JosephsonCircuits.AtoB!,JosephsonCircuits.BtoA!),
                (JosephsonCircuits.ZtoA,JosephsonCircuits.AtoZ,JosephsonCircuits.ZtoA!,JosephsonCircuits.AtoZ!),
                (JosephsonCircuits.YtoA,JosephsonCircuits.AtoY,JosephsonCircuits.YtoA!,JosephsonCircuits.AtoY!),
                (JosephsonCircuits.YtoB,JosephsonCircuits.BtoY,JosephsonCircuits.YtoB!,JosephsonCircuits.BtoY!),
                (JosephsonCircuits.ZtoB,JosephsonCircuits.BtoZ,JosephsonCircuits.ZtoB!,JosephsonCircuits.BtoZ!),
                (JosephsonCircuits.ZtoY,JosephsonCircuits.YtoZ,JosephsonCircuits.ZtoY!,JosephsonCircuits.YtoZ!),
            ]
            # single matrix input
            for arg1 in [rand(Complex{Float64},2,2), (StaticArrays.@MMatrix rand(Complex{Float64},2,2))]
                arg2 = f[1](arg1)
                arg3 = f[2](arg2)
                @test isapprox(arg1,arg3)
                arg4 = copy(arg1)
                @test isapprox(arg1,f[4](f[3](arg4)))
            end
            # array input
            for arg1 in [rand(Complex{Float64},2,2,10)]
                arg2 = f[1](arg1)
                arg3 = f[2](arg2)
                @test isapprox(arg1,arg3)
                arg4 = copy(arg1)
                @test isapprox(arg1,f[4](f[3](arg4)))
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

    @testset "Different types of conversions" begin

        S = rand(Complex{Float64},2,2)

        @test isapprox(JosephsonCircuits.ZtoA(JosephsonCircuits.StoZ(S)),JosephsonCircuits.StoA(S))

        @test isapprox(JosephsonCircuits.BtoS(JosephsonCircuits.StoB(S)),S)

        @test isapprox(JosephsonCircuits.AtoS(JosephsonCircuits.BtoA(JosephsonCircuits.StoB(S))),S)

        @test isapprox(S,JosephsonCircuits.AtoS(JosephsonCircuits.StoA(S)))

        @test isapprox(JosephsonCircuits.AtoB(JosephsonCircuits.ZtoA(JosephsonCircuits.StoZ(S))),JosephsonCircuits.StoB(S))

        @test isapprox(S,JosephsonCircuits.BtoS(JosephsonCircuits.StoB(S)))

        A = rand(Complex{Float64},2,2)

        @test isapprox(JosephsonCircuits.AtoS(A),JosephsonCircuits.ABCDtoS(A))

        @test isapprox(JosephsonCircuits.StoA(S),JosephsonCircuits.StoABCD(S))

    end

    @testset "sizes the conversions take" begin
        # the two port conversions take 2 by 2 matrices, and the conversions
        # which split the ports into inputs and outputs square matrices of
        # even size
        @test_throws DimensionMismatch JosephsonCircuits.ABCDtoS(rand(Complex{Float64}, 4, 4))
        @test_throws DimensionMismatch JosephsonCircuits.StoABCD(rand(Complex{Float64}, 4, 4))
        @test_throws DimensionMismatch JosephsonCircuits.ZtoA(rand(Complex{Float64}, 3, 3))
        @test_throws DimensionMismatch JosephsonCircuits.StoT(rand(Complex{Float64}, 3, 3))
    end

    @testset "element types of the conversions" begin
        # integer input, whose conversion is not integral: a 100 Ohm load
        # on each 50 Ohm port reflects 1/3, and a 50 Ohm series impedance
        # between two reflects 1/3 and transmits 2/3
        @test isapprox(JosephsonCircuits.ZtoS([100 0; 0 100]), [1/3 0; 0 1/3])
        @test isapprox(JosephsonCircuits.ZtoS(reshape([100, 200], 1, 1, 2)),
            reshape([1/3, 3/5], 1, 1, 2))
        @test isapprox(JosephsonCircuits.ABCDtoS([1 50; 0 1]), [1/3 2/3; 2/3 1/3])
        @test isapprox(JosephsonCircuits.AtoS([1 50; 0 1]), [1/3 2/3; 2/3 1/3])

        # real input with complex port impedances converts as complex input
        for z in (50.0 + 10.0im, [50.0 + 10.0im, 30.0 - 5.0im])
            Z = [60.0 10.0; 10.0 60.0]
            @test isapprox(JosephsonCircuits.ZtoS(Z; portimpedances = z),
                JosephsonCircuits.ZtoS(complex(Z); portimpedances = z))
            A = [1.0 50.0; 0.0 1.0]
            @test isapprox(JosephsonCircuits.ABCDtoS(A; portimpedances = z),
                JosephsonCircuits.ABCDtoS(complex(A); portimpedances = z))
            @test isapprox(JosephsonCircuits.AtoS(A; portimpedances = z),
                JosephsonCircuits.AtoS(complex(A); portimpedances = z))
            @test isapprox(JosephsonCircuits.StoZ(cat(Z/100, Z/200; dims = 3); portimpedances = z),
                JosephsonCircuits.StoZ(complex(cat(Z/100, Z/200; dims = 3)); portimpedances = z))
        end
    end


end