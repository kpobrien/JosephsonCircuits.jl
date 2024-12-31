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
                    rand(Complex{Float64}), rand(Complex{Float64},2),
                    (StaticArrays.@MVector rand(Complex{Float64},2))
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

end