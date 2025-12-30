using JosephsonCircuits
using LinearAlgebra
using Test
import StaticArrays

@testset verbose=true "networks" begin

    # one argument functions
    for fname in [:ABCD_seriesZ, :ABCD_shuntY, :Y_seriesY, :Z_shuntZ]
        f! = @eval JosephsonCircuits.$(Symbol(String(fname)*"!"))
        f = @eval JosephsonCircuits.$fname

        # check errors
        @testset "$(f), $(f!) errors" begin
            @test_throws(
                ArgumentError("Size of output (3, 2) must be (2, 2)."),
                f!(zeros(Complex{Float64},3,2),1.0)
            )
            @test_throws(
                ArgumentError("Sizes of output (2, 2) and input (1,) not compatible."),
                f!(zeros(Complex{Float64},2,2),ones(1))
            )
        end

        # check consistency for different input types
        @testset "$(f), $(f!) consistency" begin
            x1 = rand(Complex{Float64})
            # out-of-place, scalar input
            y1 = f(x1)
            # out-of-place, vector input
            @test isapprox(
                y1,
                f([x1]),
            )
            # in-place matrix, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2),x1),
            )
            # in-place array, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2,1),x1),
            )
            # in-place array, vector input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2,1),[x1]),
            )
        end
    end

    # two argument functions, complex
    for fname in [:ABCD_tline, :Z_tline]
        f! = @eval JosephsonCircuits.$(Symbol(String(fname)*"!"))
        f = @eval JosephsonCircuits.$fname

        # check errors
        @testset "$(f), $(f!) errors" begin
            @test_throws(
                ArgumentError("Sizes of inputs (2,) and (1,) must be equal."),
                f(ones(2),ones(1))
            )
            @test_throws(
                ArgumentError("Size of output (3, 2) must be (2, 2)."),
                f!(zeros(Complex{Float64},3,2),1.0,1.0)
            )
            @test_throws(
                ArgumentError("Sizes of output (2, 2) and inputs (1,) not compatible."),
                f!(zeros(Complex{Float64},2,2),ones(1),ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of output (2, 2) and inputs (1,) not compatible."),
                f!(zeros(Complex{Float64},2,2),1.0,ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of inputs (2,) and (1,) must be equal."),
                f!(zeros(Complex{Float64},2,2,1),ones(2),ones(1))
            )
        end

        # check consistency for different input types
        @testset "$(f), $(f!) consistency" begin
            x1 = rand(Complex{Float64})
            x2 = rand(Complex{Float64})
            # out-of-place, scalar input
            y1 = f(x1,x2)
            # out-of-place, vector input
            @test isapprox(
                y1,
                f([x1],[x2]),
            )
            # out-of-place, scalar and vector input
            @test isapprox(
                y1,
                f(x1,[x2]),
            )
            # in-place matrix, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2),x1,x2),
            )
            # in-place array, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2,1),x1,x2),
            )
            # in-place array, vector input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2,1),[x1],[x2]),
            )
            # in-place array, scalar and vector input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2,1),x1,[x2]),
            )
        end
    end

    # three argument functions
    for fname in [:ABCD_PiY, :Y_PiY, :ABCD_TZ, :Z_TZ]
        f! = @eval JosephsonCircuits.$(Symbol(String(fname)*"!"))
        f = @eval JosephsonCircuits.$fname

        # check errors
        @testset "$(f), $(f!) errors" begin
            @test_throws(
                ArgumentError("Size of output (3, 2) must be (2, 2)."),
                f!(zeros(Complex{Float64},3,2),1.0,1.0,1.0)
            )
            @test_throws(
                ArgumentError("Sizes of inputs (2,), (1,), and (1,) must be equal."),
                f(ones(2),ones(1),ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of inputs (2,), (1,), and (1,) must be equal."),
                f!(zeros(Complex{Float64},2,2),ones(2),ones(1),ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of output (3, 2) and inputs (1,) not compatible."),
                f!(zeros(Complex{Float64},3,2),ones(1),ones(1),ones(1))
            )
        end

        # check consistency for different input types
        @testset "$(f), $(f!) consistency" begin
            x1 = rand(Complex{Float64})
            x2 = rand(Complex{Float64})
            x3 = rand(Complex{Float64})
            # out-of-place, scalar input
            y1 = f(x1,x2,x3)
            # out-of-place, vector input
            @test isapprox(
                y1,
                f([x1],[x2],[x3]),
            )
            # in-place matrix, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2),x1,x2,x3),
            )
            # in-place array, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2,1),x1,x2,x3),
            )
            # in-place array, vector input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},2,2,1),[x1],[x2],[x3]),
            )
        end
    end

    # four argument functions, complex
    for fname in [:ABCD_coupled_tline, :Z_coupled_tline]
        f! = @eval JosephsonCircuits.$(Symbol(String(fname)*"!"))
        f = @eval JosephsonCircuits.$fname

        # check errors
        @testset "$(f), $(f!) errors" begin
            @test_throws(
                ArgumentError("Size of output (4, 3) must be (4, 4)."),
                f!(zeros(Complex{Float64},4,3),1.0,1.0,1.0,1.0)
            )
            @test_throws(
                ArgumentError("Sizes of inputs (2,), (1,), (1,), and (1,) must be equal."),
                f(ones(2),ones(1),ones(1),ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of inputs (2,) and (1,) must be equal."),
                f(1.0,1.0,ones(2),ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of inputs (2,) and (1,) must be equal."),
                f!(zeros(Complex{Float64},4,4),1.0,1.0,ones(2),ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of output (4, 4) and inputs (2,) not compatible."),
                f!(zeros(Complex{Float64},4,4),1.0,1.0,ones(2),ones(2))
            )
            @test_throws(
                ArgumentError("Sizes of inputs (1,), (1,), (2,), and (1,) must be equal."),
                f!(zeros(Complex{Float64},4,4),ones(1),ones(1),ones(2),ones(1))
            )
            @test_throws(
                ArgumentError("Sizes of output (4, 4) and inputs (2,) not compatible."),
                f!(zeros(Complex{Float64},4,4),ones(2),ones(2),ones(2),ones(2))
            )
        end

        # check consistency for different input types
        @testset "$(f), $(f!) consistency" begin
            x1 = rand(Complex{Float64})
            x2 = rand(Complex{Float64})
            x3 = rand(Complex{Float64})
            x4 = rand(Complex{Float64})
            # out-of-place, scalar input
            y1 = f(x1,x2,x3,x4)
            # out-of-place, vector input
            @test isapprox(
                y1,
                f([x1],[x2],[x3],[x4]),
            )
            # out-of-place, scalar and vector input
            @test isapprox(
                y1,
                f(x1,x2,[x3],[x4]),
            )
            # in-place matrix, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},4,4),x1,x2,x3,x4),
            )
            # in-place array, scalar input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},4,4,1),x1,x2,x3,x4),
            )
            # in-place array, vector input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},4,4,1),[x1],[x2],[x3],[x4]),
            )
            # in-place array, scalar and vector input
            @test isapprox(
                y1,
                f!(zeros(Complex{Float64},4,4,1),x1,x2,[x3],[x4]),
            )
        end
    end


    for f in [JosephsonCircuits.S_short!,JosephsonCircuits.S_open!,JosephsonCircuits.S_match!,JosephsonCircuits.S_splitter!]

        @testset "$(f) errors" begin
            @test_throws(
                ArgumentError("The sizes of the first two dimensions (3,2) of the scattering matrix `S` must be the same."),
                f(ones(3,2))
            )
            @test_throws(
                ArgumentError("The scattering matrix `S` with size (1,) must have at least two dimensions."),
                f(ones(1))
            )
        end
    end

    @testset "consistency checks" begin

        x1 = rand(Complex{Float64})
        x2 = rand(Complex{Float64})
        x3 = rand(Complex{Float64})
        x4 = rand(Complex{Float64})

        @test isapprox(
            JosephsonCircuits.ABCD_seriesZ(x1),
            JosephsonCircuits.YtoA(JosephsonCircuits.Y_seriesY(1/x1)),
        )

        @test isapprox(
            JosephsonCircuits.ABCD_shuntY(1/x1),
            JosephsonCircuits.ZtoA(JosephsonCircuits.Z_shuntZ(x1)),
        )

        @test isapprox(
            JosephsonCircuits.ABCD_tline(x1,x2),
            JosephsonCircuits.ZtoA(JosephsonCircuits.Z_tline(x1,x2)),
        )

        @test isapprox(
            JosephsonCircuits.ABCD_PiY(x1,x2,x3),
            JosephsonCircuits.YtoA(JosephsonCircuits.Y_PiY(x1,x2,x3)),
        )

        @test isapprox(
            JosephsonCircuits.ABCD_TZ(x1,x2,x3),
            JosephsonCircuits.ZtoA(JosephsonCircuits.Z_TZ(x1,x2,x3)),
        )

        @test isapprox(
            JosephsonCircuits.ABCD_coupled_tline(x1,x2,x3,x4),
            JosephsonCircuits.ZtoA(JosephsonCircuits.Z_coupled_tline(x1,x2,x3,x4)),
        )

    end

    @testset "Z_invC, Y_C, Z_L, Y_invL, Z_L!, Z_invC!, Y_C!, Y_invL!" begin

        C = rand(Complex{Float64},2,2)
        w=rand()
        @test isapprox(JosephsonCircuits.Z_invC(inv(C),w),inv(JosephsonCircuits.Y_C(C,w)))

        L = rand(Complex{Float64},2,2)
        w=rand()
        @test isapprox(JosephsonCircuits.Z_L(L,w),inv(JosephsonCircuits.Y_invL(inv(L),w)))

        L = rand(Complex{Float64},2,2)
        w=rand()
        Z = zeros(Complex{Float64},2,2)
        @test isapprox(JosephsonCircuits.Z_L(L,w),JosephsonCircuits.Z_L!(Z,L,w))

        invC = rand(Complex{Float64},2,2)
        w=rand()
        Z = zeros(Complex{Float64},2,2)
        @test isapprox(JosephsonCircuits.Z_invC(invC,w),JosephsonCircuits.Z_invC!(Z,invC,w))

        C = rand(Complex{Float64},2,2)
        w=rand()
        Y = zeros(Complex{Float64},2,2)
        @test isapprox(JosephsonCircuits.Y_C(C,w),JosephsonCircuits.Y_C!(Y,C,w))

        invL = rand(Complex{Float64},2,2)
        w=rand()
        Y = zeros(Complex{Float64},2,2)
        @test isapprox(JosephsonCircuits.Y_invL(invL,w),JosephsonCircuits.Y_invL!(Y,invL,w))

    end

    @testset "Z_L!, Z_invC!, Y_C!, Y_invL! errors" begin

        A = rand(Complex{Float64},2,2)
        w=rand()
        B = zeros(Complex{Float64},4,4)
        @test_throws(
            ArgumentError("The size of the input (4, 4) must equal the size of the output (2, 2)."),
            JosephsonCircuits.Z_L!(A,B,w),
        )
        @test_throws(
            ArgumentError("The size of the input (4, 4) must equal the size of the output (2, 2)."),
            JosephsonCircuits.Z_invC!(A,B,w),
        )
        @test_throws(
            ArgumentError("The size of the input (4, 4) must equal the size of the output (2, 2)."),
            JosephsonCircuits.Y_C!(A,B,w),
        )
        @test_throws(
            ArgumentError("The size of the input (4, 4) must equal the size of the output (2, 2)."),
            JosephsonCircuits.Y_invL!(A,B,w),
        )
    end

    @testset "ZC_basis_coupled_lines" begin
        Cg = 1.6670474181399462e-10
        Cm = 9.320861870486729e-12
        Ls = 4.167618545349866e-7
        Lm = 2.330215467621682e-8
        C = JosephsonCircuits.LinearAlgebra.Symmetric([Cg -Cm;-Cm Cg])
        L = JosephsonCircuits.LinearAlgebra.Symmetric([Ls Lm;Lm Ls])
        b = JosephsonCircuits.ZC_basis_coupled_tlines(L,C)
        ZC = b.ZC
        TI = b.TI
        TV = b.TV
        theta = b.theta
        U = b.U
        lambda = b.lambda
        S = b.S
        @test isapprox(U'*C*U,theta^2)
        @test isapprox(S'*(theta*U'*L*U*theta)*S,lambda^2)
        @test isapprox(inv(TI)*C*L*TI,lambda^2)
        @test isapprox(TV'*TI,I)
    end

    @testset "A_coupled_tlines" begin

        Zeven = 51.0
        Zodd = 49.0
        neven = 1.1
        nodd = 1.08
        l = 3.5e-3
        c = JosephsonCircuits.speed_of_light
        omega = 2*pi*(1:10)*1e9

        L, C = JosephsonCircuits.even_odd_to_maxwell(Zeven, Zodd, neven, nodd)
        A1 = JosephsonCircuits.A_coupled_tlines(L,C,l,omega)
        A2 = stack(JosephsonCircuits.ABCD_coupled_tline.(Zeven,Zodd,neven*omega/c*l,nodd*omega/c*l))
        @test isapprox(A1,A2)

    end

    @testset "Z_canonical_coupled_line_circuits" begin
        Zeven = 52.0
        Zodd = 48.0
        neven = 2.2
        nodd = 2.1
        l = 2*3.0e-3
        c = 2.998e8
        omega = 2*pi*5e9

        function S_canonical_coupled_line_circuits(i::Int, Ze, Zo, thetae, thetao)

            # define an open
            Sopen = ones(Complex{Float64},1,1)

            # and a short
            Sshort = -ones(Complex{Float64},1,1)
            
            # and a coupled transmission line
            S = JosephsonCircuits.ZtoS(JosephsonCircuits.Z_coupled_tline(Zeven,Zodd,thetae,thetao));

            if i == 1
                S = JosephsonCircuits.interconnectS(S,Sopen,3,1)
                S = JosephsonCircuits.interconnectS(S,Sshort,1,1)
            elseif i == 2
                S = JosephsonCircuits.interconnectS(S,Sshort,4,1)
                S = JosephsonCircuits.interconnectS(S,Sshort,1,1)
            elseif i == 3
                S = JosephsonCircuits.interconnectS(S,Sopen,4,1)
                S = JosephsonCircuits.interconnectS(S,Sopen,1,1)
            elseif i == 4
                S = JosephsonCircuits.interconnectS(S,Sopen,4,1)
                S = JosephsonCircuits.interconnectS(S,Sshort,3,1)
            elseif i == 5
                S = JosephsonCircuits.interconnectS(S,Sopen,3,1)
                S = JosephsonCircuits.interconnectS(S,Sopen,1,1)
            elseif i == 6
                S = JosephsonCircuits.interconnectS(S,Sshort,3,1)
                S = JosephsonCircuits.interconnectS(S,Sshort,1,1)
            elseif i == 7
                S = JosephsonCircuits.intraconnectS(S,4,3)
            elseif i == 8
                S = JosephsonCircuits.interconnectS(S,Sopen,4,1)
                S = JosephsonCircuits.interconnectS(S,Sshort,1,1)
            elseif i == 9
                S = JosephsonCircuits.interconnectS(S,Sshort,4,1)
                S = JosephsonCircuits.interconnectS(S,Sshort,3,1)
            elseif i == 10
                S = JosephsonCircuits.interconnectS(S,Sopen,4,1)
                S = JosephsonCircuits.interconnectS(S,Sopen,3,1)
            else
                throw(ArgumentError("Canonical coupled line circuit number must be 1-10."))
            end
                
            return S
        end

        for i in 1:10
            S1 = S_canonical_coupled_line_circuits(i, Zeven, Zodd, neven*omega/c*l, nodd*omega/c*l)
            S2 = JosephsonCircuits.ZtoS(JosephsonCircuits.Z_canonical_coupled_line_circuits(i, Zeven, Zodd, neven*omega/c*l, nodd*omega/c*l));
            @test isapprox(S1,S2)
        end

        @test_throws(
            ArgumentError("Canonical coupled line circuit number must be 1-10."),
            JosephsonCircuits.Z_canonical_coupled_line_circuits(11, Zeven, Zodd, neven*omega/c*l, nodd*omega/c*l),
        )

    end

    @testset "canonical_coupled_line_circuits" begin

        Zeven = 52.0
        Zodd = 48.0
        neven = 2.2
        nodd = 2.1

        # 3
        outa = (L1 = 1.1963832214440388e-7, L2 = 1.1963832214440388e-7, M = -3.780393078912389e-9, C1 = 1.4112327104537203e-10, C2 = 1.4112327104537203e-10, Cm = 2.4055103019097465e-12)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(3, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])

        # 8
        outa = (L1 = 1.1941849319292907e-7, L2 = 3.5891496643321166e-7, M = -8.333511920257754e-9, C1 = 1.4352878134728178e-10)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(8, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])

        # 9
        outa = (L1 = 3.5891496643321166e-7, L2 = 3.5891496643321166e-7, M = 2.268235847347433e-8)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(9, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])

        # 10
        outa = (L1 = 1.1963832214440388e-7, L2 = 1.1963832214440388e-7, M = 7.560786157824777e-9, C1 = 1.4112327104537203e-10, C2 = 1.4112327104537203e-10, Cm = 2.4055103019097465e-12)
        outb = JosephsonCircuits.canonical_coupled_line_circuits(10, Zeven, Zodd, neven, nodd)
        @test all([isapprox(outai,outbi) for (outai,outbi) in zip(outa,outb)])


        @test_throws(
            ArgumentError("Canonical coupled line circuit number must be 1-10."),
            JosephsonCircuits.canonical_coupled_line_circuits(11, Zeven, Zodd, neven, nodd),
        )

    end

    @testset "directional couplers" begin

        Z0 = 50.0
        couplingdB = 20

        # compute the scattering parameters for an ideal directional coupler
        S = JosephsonCircuits.S_directional_coupler_symmetric(couplingdB)
        @test isapprox(-10*log10(abs2(S[3,1])),couplingdB)

        # compute the even and odd mode characteristic impedance for a
        # directional coupler, then compare this with the scattering
        # parameters from coupled transmission lines.
        Zeven, Zodd = JosephsonCircuits.coupling_to_even_odd(couplingdB,Z0)
        theta = pi/2
        S = JosephsonCircuits.AtoS(JosephsonCircuits.ABCD_coupled_tline(Zeven,Zodd,theta,theta));
        @test isapprox(-10*log10(abs2(S[2,1])),couplingdB)


        # test that the in-place and out-of-place directional coupler,
        # hybrid coupler, and circulator functions give the same results
        S1 = JosephsonCircuits.S_directional_coupler_symmetric(couplingdB)
        S2 = similar(S1)
        JosephsonCircuits.S_directional_coupler_symmetric!(S2,couplingdB)
        @test isequal(S1,S2)

        S1 = JosephsonCircuits.S_directional_coupler_antisymmetric(couplingdB)
        S2 = similar(S1)
        JosephsonCircuits.S_directional_coupler_antisymmetric!(S2,couplingdB)
        @test isequal(S1,S2)

        S1 = JosephsonCircuits.S_hybrid_coupler_symmetric()
        S2 = similar(S1)
        JosephsonCircuits.S_hybrid_coupler_symmetric!(S2)
        @test isequal(S1,S2)

        S1 = JosephsonCircuits.S_hybrid_coupler_antisymmetric()
        S2 = similar(S1)
        JosephsonCircuits.S_hybrid_coupler_antisymmetric!(S2)
        @test isequal(S1,S2)

        # test for array input
        S1 = JosephsonCircuits.S_directional_coupler_symmetric(couplingdB)
        S2 = similar(S1,4,4,10)
        JosephsonCircuits.S_directional_coupler_symmetric!(S2,couplingdB)
        for i in 1:size(S2,3)
            @test isequal(S1,S2[:,:,i])
        end

        @test_throws(
            ArgumentError("Size of output (2, 2) must be (4, 4)."),
            JosephsonCircuits.S_directional_coupler!(zeros(Complex{Float64},2,2),1,2,3,4),
        )

    end

    @testset "attenuators" begin

        @test_throws(
            ArgumentError("Attenuation of 0 dB is below the minimum attenuation of 3.7653971144370946 dB for a passive circuit given the source and load impedances of 50 and 60 Ohms."),
            JosephsonCircuits.Z_attenuator_inner(50,60,0),
        )
    end


    @testset "circulators" begin

        S1 = JosephsonCircuits.S_circulator_clockwise()
        S2 = similar(S1)
        S3 = zeros(Complex{Float64},3,3,1)
        JosephsonCircuits.S_circulator_clockwise!(S2)
        JosephsonCircuits.S_circulator_clockwise!(S3)
        @test isequal(S1,S2)
        @test isequal(S1,S3[:,:,1])

        @test_throws(
            ArgumentError("Size of output (4, 4) must be (3, 3)."),
            JosephsonCircuits.S_circulator_clockwise!(zeros(Complex{Float64},4,4)),
        )

        S1 = JosephsonCircuits.S_circulator_counterclockwise()
        S2 = similar(S1)
        S3 = zeros(Complex{Float64},3,3,1)
        JosephsonCircuits.S_circulator_counterclockwise!(S2)
        JosephsonCircuits.S_circulator_counterclockwise!(S3)
        @test isequal(S1,S2)
        @test isequal(S1,S3[:,:,1])

        @test_throws(
            ArgumentError("Size of output (4, 4) must be (3, 3)."),
            JosephsonCircuits.S_circulator_counterclockwise!(zeros(Complex{Float64},4,4)),
        )

    end

end
