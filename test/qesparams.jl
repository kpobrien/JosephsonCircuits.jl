using JosephsonCircuits
using Test

@testset verbose=true "qesparams" begin

    @testset "calcimpedance errors" begin

        begin
            inputwave=[1.0, 0.0]
            outputwave=[im/sqrt(2), 1/sqrt(2), 0]
            S = zeros(Complex{Float64},2,2)
            @test_throws(
                DimensionMismatch("First dimension of scattering matrix not consistent with first dimensions of outputwave."),
                JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave))
        end

        begin
            inputwave=[1.0, 0.0, 0.0]
            outputwave=[im/sqrt(2), 1/sqrt(2)]
            S = zeros(Complex{Float64},2,2)
            @test_throws(
                DimensionMismatch("Second dimension of scattering matrix not consistent with first dimension of input wave."),
                JosephsonCircuits.calcscatteringmatrix!(S,inputwave,outputwave))
        end

        begin
            @test_throws(
                ErrorException("Unknown component type"),
                JosephsonCircuits.calcimpedance(30.0,:D,-1.0,nothing))
        end

        begin
            @variables w
            @test_throws(
                ErrorException("Unknown component type"),
                JosephsonCircuits.calcimpedance(30*w,:D,-2.0,w))
        end
    end

    @testset "calccm! errors" begin
        begin
            cm=Float64[0,0]
            @test_throws(
                DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
                JosephsonCircuits.calccm!(cm,[3/5 4/5;4/5 3/5],[-1,1,2]))
        end

        begin
            cm=Float64[0,0]
            @test_throws(
                DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
                JosephsonCircuits.calccm!(cm,[3/5 4/5;4/5 3/5;0 0;0 0],[-1,1]))
        end

        begin
            @variables a b
            cm=Num[0,0]
            @test_throws(
                DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
                JosephsonCircuits.calccm!(cm,[a b; b a],[-1,1,2]))
        end

        begin
            @variables a b
            cm=Num[0,0]
            @test_throws(
                DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
                JosephsonCircuits.calccm!(cm,[a b; b a; 0 0; 0 0],[-1,1]))
        end

        begin
            cm=Float64[0, 0]
            @test_throws(
                DimensionMismatch("Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."),
                JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2 3;5 6 7],[-1,1]))
        end

        begin
            cm=Float64[0, 0]
            @test_throws(
                DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."),
                JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2; 3 4; 5 6; 7 8],[-1,1]))
        end

        begin
            cm=Float64[0, 0]
            @test_throws(
                DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
                JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2 3 4;5 6 7 8],[-1,1,2]))
        end

        begin
            cm=Float64[0, 0, 0]
            @test_throws(
                DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
                JosephsonCircuits.calccm!(cm,[1 2;3 4],[1 2 3 4;5 6 7 8],[-1,1]))
        end

        begin
            @variables a b c d an bn cn dn
            cm = Num[0, 0]
            @test_throws(
                DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
                JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn],[1, -1, 2]))
        end

        begin
            @variables a b c d an bn cn dn;cm = Num[0, 0]
            @test_throws(
                DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."),
                JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn; 0 0; 0 0],[1, -1]))
        end

        begin
            @variables a b c d an bn cn dn
            cm = Num[0, 0, 0]
            @test_throws(
                DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
                JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn; cn dn],[1, -1]))
        end

        begin
            @variables a b c d an bn cn dn
            cm = Num[0, 0]
            @test_throws(
                DimensionMismatch("Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."),
                JosephsonCircuits.calccm!(cm,Num[a b; c d],[an bn 0; cn dn 0],[1, -1]))
        end
    end

    @testset "calcqe! errors" begin
        @test_throws(
            DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcqe!([1 2;3 4],[1 2 3;4 5 6]))

        @test_throws(
            DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcqe!([1 2;3 4],[1 2 3;4 5 6],[1 2;3 4]))

        @test_throws(
            DimensionMismatch("First dimensions of scattering parameter matrice and noise scattering matrix must be equal."),
            JosephsonCircuits.calcqe!(Float64[1 2;3 4],[1 2;3 4],[1 2;3 4;5 6]))
    end

    @testset "calcqeideal!" begin
        @test_throws(
            DimensionMismatch("Sizes of QE and S matrices must be equal."),
            JosephsonCircuits.calcqeideal!([1 2;3 4],[1 2 3;4 5 6]))
    end

    @testset "calcdZdroZ2" begin
        @test_throws(
            ArgumentError("Unknown component."),
            JosephsonCircuits.calcdZdroZ2([1],[:K], [2.0], [1.0],nothing))
    end

    @testset "noise wave covariance matrice and QE" begin

        # symbolic
        N = 3
        for i in 1:N
            indices = collect(1:N)
            popat!(indices,i)

            # generate the `S` matrices. assume `S` is the scattering
            # parameter matrix for a lossless network.
            @variables S11 S12 S13 S21 S22 S23 S31 S32 S33
            S = [S11 S12 S13; S21 S22 S23; S31 S32 S33]

            # pick one port and imagine that it is a resistor with
            # resistance equal to the port impedance. Snoise represents noise emerging
            # from the resistor and propagating to the other ports.
            Snoise = transpose(S[indices,i])

            # generate the noise wave covariance matrices `C`
            # C1 will be zero for a passive network and C2 will be non-zero since we
            # replaced the port with a resistor.
            # C1 = JosephsonCircuits.calcCnoise(S)
            C = JosephsonCircuits.calcCnoise(S[indices,indices],transpose(Snoise))

            # test that the QE's are equal for the original network and the
            # reduced network, with the scattering parameter based QE calculation
            QE1 = JosephsonCircuits.calcqe(S)[indices,indices]
            QE2 = JosephsonCircuits.calcqe(S[indices,indices],transpose(Snoise))
            @test isequal(QE1,QE2)

            # test the QE computed from the covariance matrix is the same
            QE3 = JosephsonCircuits.calcqe_S_Cnoise(S[indices,indices],C)
            @test isequal(QE1,QE3)
        end

        # numeric
        N = 3
        for i in 1:N
            indices = collect(1:N)
            popat!(indices,i)

            # generate the `S` matrices. assume `S` is the scattering
            # parameter matrix for a lossless network.
            S = zeros(Complex{Float64},N,N)

            # pick one port and imagine that it is a resistor with
            # resistance equal to the port impedance. Snoise represents noise emerging
            # from the resistor and propagating to the other ports.
            Snoise = transpose(S[indices,i])

            # generate the noise wave covariance matrices `C`
            # C1 will be zero for a passive network and C2 will be non-zero since we
            # replaced the port with a resistor.
            # C1 = JosephsonCircuits.calcCnoise(S)
            C = JosephsonCircuits.calcCnoise(S[indices,indices],transpose(Snoise))

            # test that the QE's are equal for the original network and the
            # reduced network, with the scattering parameter based QE calculation
            QE1 = JosephsonCircuits.calcqe(S)[indices,indices]
            QE2 = JosephsonCircuits.calcqe(S[indices,indices],transpose(Snoise))
            @test isequal(QE1,QE2)

            # test the QE computed from the covariance matrix is the same
            QE3 = JosephsonCircuits.calcqe_S_Cnoise(S[indices,indices],C)
            @test isequal(QE1,QE3)
        end

    end

end
