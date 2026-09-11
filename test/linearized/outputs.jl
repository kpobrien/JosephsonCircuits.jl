using JosephsonCircuits, LinearAlgebra, Random
using Test

@testset verbose=true "the linearized outputs" begin

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
            JosephsonCircuits.@params w
            @test_throws(
                ErrorException("Unknown component type"),
                JosephsonCircuits.calcimpedance(30*w,:D,-2.0,w))
        end
    end

    @testset "calccm! errors" begin
        cm=Float64[0,0]
        @test_throws(
            DimensionMismatch("Dimensions of scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.calccm!(cm,[3/5 4/5;4/5 3/5],[-1,1,2]))
        @test_throws(
            DimensionMismatch("First dimension of scattering matrix must equal the length of cm."),
            JosephsonCircuits.calccm!(cm,[3/5 4/5;4/5 3/5;0 0;0 0],[-1,1]))
        @test_throws(
            DimensionMismatch("Dimensions of noise scattering matrix must be integer multiples of the number of frequencies."),
            JosephsonCircuits.noisereduction([1 2;3 4;5 6],[-1,1]))
        @test_throws DimensionMismatch JosephsonCircuits.noisereduction!(
            JosephsonCircuits.NoiseReduction(zeros(3), zeros(3)), [1 2;3 4], [-1,1])
        @test_throws DimensionMismatch JosephsonCircuits.noisereduction(
            [1 2;3 4], [-1,1], [1.0, 1.0, 1.0])
        @test_throws(
            DimensionMismatch("First dimension of the scattering parameter matrix must equal the length of the noise reduction."),
            JosephsonCircuits.calccm!(cm,[1 2;3 4],[-1,1],
                JosephsonCircuits.noisereduction([1 2 3;4 5 6],[-1,1])))
        @test_throws DimensionMismatch JosephsonCircuits.weightedrowpower!(
            zeros(3), zeros(3), [1 2;3 4], nothing)
    end

    @testset "calccm! and calcqe!" begin
        # a high gain amplifier's row is a cancellation between the signal
        # and idler terms, which the compensated sum resolves to the unit
        # commutator
        G = 1e8
        S = [sqrt(G) sqrt(G-1); sqrt(G-1) sqrt(G)]
        w = [1, -1]
        @test JosephsonCircuits.calccm(S, w) ≈ [1.0, -1.0] atol = 1e-6
        @test JosephsonCircuits.calccm(S, w) == JosephsonCircuits.calccm(S .+ 0im, w)

        # the noise reduction against the explicit sums, with and without an
        # occupation, and against the diagonal of the noise covariance
        rng = Random.default_rng()
        m = 3
        Snoise = randn(rng, ComplexF64, 5*m, 2*m)
        w = randn(rng, m)
        occ = 1 .+ rand(rng, 5*m)
        n0 = JosephsonCircuits.noisereduction(Snoise, w)
        n1 = JosephsonCircuits.noisereduction(Snoise, w, occ)
        @test n0.denom ≈ vec(sum(abs2, Snoise; dims = 1))
        @test n1.denom ≈ vec(sum(occ .* abs2.(Snoise); dims = 1))
        signs = [sign(w[(c-1) % m + 1]) for c in 1:5*m]
        @test n0.signed ≈ vec(sum(signs .* abs2.(Snoise); dims = 1))
        @test n1.signed == n0.signed
        C = zeros(ComplexF64, 2*m, 2*m)
        JosephsonCircuits.calcnoisecovariance!(C, Snoise, occ)
        @test n1.denom ≈ real(diag(C))

        # the three quantum efficiencies: from the reduction, from the
        # covariance, and from the explicit formula
        S = randn(rng, ComplexF64, 2*m, 2*m)
        qe = JosephsonCircuits.calcqe(S, n1)
        @test qe ≈ JosephsonCircuits.calcqe_S_Cnoise(S, C)
        @test qe ≈ abs2.(S) ./ (vec(sum(abs2, S; dims = 2)) .+ n1.denom)
        @test JosephsonCircuits.calcqe(S) ≈ abs2.(S) ./ vec(sum(abs2, S; dims = 2))
        cm = JosephsonCircuits.calccm(S, w, n0)
        colsigns = [sign(w[(j-1) % m + 1]) for j in 1:2*m]
        @test cm ≈ vec(sum(abs2.(S) .* transpose(colsigns); dims = 2)) .+ n0.signed
        # the scratch is the caller's
        qe2 = similar(qe); cm2 = similar(cm)
        JosephsonCircuits.calcqe!(qe2, S, n1; denom = zeros(2*m), comp = zeros(2*m))
        JosephsonCircuits.calccm!(cm2, S, w, n0; comp = zeros(2*m))
        @test qe2 == qe && cm2 == cm
    end

    @testset "calcqe! errors" begin
        @test_throws(
            DimensionMismatch("Dimensions of quantum efficiency and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcqe!([1 2;3 4],[1 2 3;4 5 6]))
        @test_throws(
            DimensionMismatch("First dimension of the scattering parameter matrix must equal the length of the noise reduction."),
            JosephsonCircuits.calcqe!(Float64[1 2;3 4],[1 2;3 4],
                JosephsonCircuits.noisereduction([1 2 3;4 5 6],[1])))
    end

    @testset "calcqeideal!" begin
        @test_throws(
            DimensionMismatch("Sizes of QE and S matrices must be equal."),
            JosephsonCircuits.calcqeideal!([1 2;3 4],[1 2 3;4 5 6]))
    end

    @testset "calcCnoise! errors" begin
        @test_throws(
            DimensionMismatch("The dimensions of the noise wave covariance and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcCnoise!([1 2;3 4],[1 2 3;4 5 6]))

        @test_throws(
            DimensionMismatch("The dimensions of the noise wave covariance and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcCnoise!([1 2;3 4],[1 2 3;4 5 6],[1 2;3 4]))

        @test_throws(
            DimensionMismatch("The first dimensions of the scattering parameter and noise scattering parameter matrices must be equal."),
            JosephsonCircuits.calcCnoise!([1 2;3 4],[1 2;3 4],[1 2;3 4;5 6]))

    end

    @testset "calcqe_S_Cnoise!(qe, S, Cnoise) errors" begin

        @test_throws(
            DimensionMismatch("The dimensions of the quantum efficiency and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcqe_S_Cnoise!([1 2;3 4],[1 2 3;4 5 6],[1 2;3 4]))

        @test_throws(
            DimensionMismatch("The dimensions of the noise wave covariance and scattering parameter matrices must be equal."),
            JosephsonCircuits.calcqe_S_Cnoise!([1 2;3 4],[1 2;3 4],[1 2;3 4;5 6]))
    
    end

    @testset "noise wave covariance matrice and QE" begin

        # numeric
        N = 3
        for i in 1:N
            indices = collect(1:N)
            popat!(indices,i)

            # generate the `S` matrices: a random unitary, the scattering
            # parameter matrix of a lossless network
            S = Matrix(qr(randn(Random.default_rng(), Complex{Float64}, N, N)).Q)

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
            QE2 = JosephsonCircuits.calcqe(S[indices,indices],
                JosephsonCircuits.noisereduction(Snoise, [1]))
            @test isapprox(QE1,QE2; rtol = 1e-12)

            # test the QE computed from the covariance matrix is the same
            QE3 = JosephsonCircuits.calcqe_S_Cnoise(S[indices,indices],C)
            @test isapprox(QE1,QE3; rtol = 1e-12)
        end

    end

end
