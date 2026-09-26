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

    @testset "absolute references" begin
        # a T network, series Z1 and Z2 with Z3 shunt between them, whose
        # impedance, admittance and chain matrices are known in closed form
        Z1, Z2, Z3 = 10.0 + 5.0im, 20.0 - 3.0im, 30.0 + 7.0im
        Z = [Z1+Z3 Z3; Z3 Z2+Z3]
        Y = [Z2+Z3 -Z3; -Z3 Z1+Z3] / (Z1*Z2 + Z1*Z3 + Z2*Z3)
        # [V1; I1] = A*[V2; -I2], with the currents into the ports
        A = [1+Z1/Z3 Z1+Z2+Z1*Z2/Z3; 1/Z3 1+Z2/Z3]
        # [V2; I2] = B*[V1; -I1]
        D = Diagonal([1, -1])
        B = D * inv(A) * D
        for (f, x, y) in (
                (JosephsonCircuits.ZtoY, Z, Y), (JosephsonCircuits.YtoZ, Y, Z),
                (JosephsonCircuits.ZtoA, Z, A), (JosephsonCircuits.AtoZ, A, Z),
                (JosephsonCircuits.YtoA, Y, A), (JosephsonCircuits.AtoY, A, Y),
                (JosephsonCircuits.ZtoB, Z, B), (JosephsonCircuits.BtoZ, B, Z),
                (JosephsonCircuits.YtoB, Y, B), (JosephsonCircuits.BtoY, B, Y),
                (JosephsonCircuits.AtoB, A, B), (JosephsonCircuits.BtoA, B, A))
            @test isapprox(f(x), y)
        end

        # the pseudo-waves a = (V + Zr*I)/(2*sqrt(Zr)) and
        # b = (V - Zr*I)/(2*sqrt(Zr)) at each port give
        # S = inv(G)*(Z - Zr)*inv(Z + Zr)*G with G = Diagonal(sqrt.(zr)),
        # for real, unequal, complex and capacitive references
        for zr in ([50.0, 50.0], [50.0, 75.0], [50.0 + 10.0im, 30.0 - 5.0im],
                [-50.0im, -50.0im], [-30.0 + 10.0im, 60.0 + 5.0im])
            Zr = Diagonal(zr)
            G = Diagonal(sqrt.(zr))
            S = inv(G) * (Z - Zr) * inv(Z + Zr) * G
            for (f, x, y) in (
                    (JosephsonCircuits.ZtoS, Z, S), (JosephsonCircuits.StoZ, S, Z),
                    (JosephsonCircuits.YtoS, Y, S), (JosephsonCircuits.StoY, S, Y),
                    (JosephsonCircuits.AtoS, A, S), (JosephsonCircuits.StoA, S, A),
                    (JosephsonCircuits.BtoS, B, S), (JosephsonCircuits.StoB, S, B),
                    (JosephsonCircuits.ABCDtoS, A, S), (JosephsonCircuits.StoABCD, S, A))
                @test isapprox(f(x; portimpedances = zr), y)
            end
        end

        # [b1; a1] = T*[a2; b2]
        S = JosephsonCircuits.ZtoS(Z)
        T = JosephsonCircuits.StoT(S)
        a2, b2 = 0.3 + 0.1im, -0.2 + 0.4im
        a1 = (b2 - S[2, 2]*a2) / S[2, 1]
        b1 = S[1, 1]*a1 + S[1, 2]*a2
        @test isapprox(T * [a2; b2], [b1; a1])
        @test isapprox(JosephsonCircuits.TtoS(T), S)

        # a two port which transmits nothing has no chain matrix
        @test_throws SingularException JosephsonCircuits.StoA([0.5 0.0; 0.0 0.5])
        @test_throws SingularException JosephsonCircuits.StoABCD([0.5 0.0; 0.0 0.5])
    end

    @testset "sizes the conversions take" begin
        # the two port conversions take 2 by 2 matrices, and the conversions
        # which split the ports into inputs and outputs square matrices of
        # even size
        @test_throws DimensionMismatch JosephsonCircuits.ABCDtoS(rand(Complex{Float64}, 4, 4))
        @test_throws DimensionMismatch JosephsonCircuits.StoABCD(rand(Complex{Float64}, 4, 4))
        @test_throws DimensionMismatch JosephsonCircuits.ZtoA(rand(Complex{Float64}, 3, 3))
        @test_throws DimensionMismatch JosephsonCircuits.StoT(rand(Complex{Float64}, 3, 3))
        # one port impedance per port, or for a two port one or two
        @test_throws DimensionMismatch JosephsonCircuits.StoA(
            rand(Complex{Float64}, 4, 4); portimpedances = [1.0, 2, 3, 4, 5, 6])
        @test_throws DimensionMismatch JosephsonCircuits.StoZ(
            rand(Complex{Float64}, 2, 2); portimpedances = [1.0, 2, 3])
        @test_throws DimensionMismatch JosephsonCircuits.StoZ(
            rand(Complex{Float64}, 2, 2, 3); portimpedances = rand(3, 3))
        @test_throws DimensionMismatch JosephsonCircuits.ABCDtoS(
            rand(Complex{Float64}, 2, 2); portimpedances = [1.0, 2, 3])
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