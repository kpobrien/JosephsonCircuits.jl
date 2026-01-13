using JosephsonCircuits
using Test
using SpecialFunctions

@testset verbose=true "fftutils" begin

    @testset "calcfreqsrdft" begin

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsrdft((1,)),
            JosephsonCircuits.Frequencies{1}((1,), (2,), (3,), CartesianIndex{1}[CartesianIndex(1), CartesianIndex(2)], [(0,), (1,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsrdft((2,)),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1), CartesianIndex(2), CartesianIndex(3)], [(0,), (1,), (2,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsrdft((3,)),
            JosephsonCircuits.Frequencies{1}((3,), (4,), (6,), CartesianIndex{1}[CartesianIndex(1), CartesianIndex(2), CartesianIndex(3), CartesianIndex(4)], [(0,), (1,), (2,), (3,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsrdft((3,3)),
            JosephsonCircuits.Frequencies{2}((3, 3), (4, 7), (6, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(4, 3), CartesianIndex(1, 4), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(4, 4), CartesianIndex(1, 5), CartesianIndex(2, 5), CartesianIndex(3, 5), CartesianIndex(4, 5), CartesianIndex(1, 6), CartesianIndex(2, 6), CartesianIndex(3, 6), CartesianIndex(4, 6), CartesianIndex(1, 7), CartesianIndex(2, 7), CartesianIndex(3, 7), CartesianIndex(4, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (0, 1), (1, 1), (2, 1), (3, 1), (0, 2), (1, 2), (2, 2), (3, 2), (0, 3), (1, 3), (2, 3), (3, 3), (0, -3), (1, -3), (2, -3), (3, -3), (0, -2), (1, -2), (2, -2), (3, -2), (0, -1), (1, -1), (2, -1), (3, -1)]),
        )
    end

    @testset "calcfreqsdft" begin


        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsdft((1,)),
            JosephsonCircuits.Frequencies{1}((1,), (3,), (3,), CartesianIndex{1}[CartesianIndex(1), CartesianIndex(2), CartesianIndex(3)], [(0,), (1,), (-1,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsdft((2,)),
            JosephsonCircuits.Frequencies{1}((2,), (5,), (5,), CartesianIndex{1}[CartesianIndex(1), CartesianIndex(2), CartesianIndex(3), CartesianIndex(4), CartesianIndex(5)], [(0,), (1,), (2,), (-2,), (-1,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsdft((3,)),
            JosephsonCircuits.Frequencies{1}((3,), (7,), (7,), CartesianIndex{1}[CartesianIndex(1), CartesianIndex(2), CartesianIndex(3), CartesianIndex(4), CartesianIndex(5), CartesianIndex(6), CartesianIndex(7)], [(0,), (1,), (2,), (3,), (-3,), (-2,), (-1,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.calcfreqsdft((3,3)),
            JosephsonCircuits.Frequencies{2}((3, 3), (7, 7), (7, 7), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(4, 1), CartesianIndex(5, 1), CartesianIndex(6, 1), CartesianIndex(7, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(4, 2), CartesianIndex(5, 2), CartesianIndex(6, 2), CartesianIndex(7, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(4, 3), CartesianIndex(5, 3), CartesianIndex(6, 3), CartesianIndex(7, 3), CartesianIndex(1, 4), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(4, 4), CartesianIndex(5, 4), CartesianIndex(6, 4), CartesianIndex(7, 4), CartesianIndex(1, 5), CartesianIndex(2, 5), CartesianIndex(3, 5), CartesianIndex(4, 5), CartesianIndex(5, 5), CartesianIndex(6, 5), CartesianIndex(7, 5), CartesianIndex(1, 6), CartesianIndex(2, 6), CartesianIndex(3, 6), CartesianIndex(4, 6), CartesianIndex(5, 6), CartesianIndex(6, 6), CartesianIndex(7, 6), CartesianIndex(1, 7), CartesianIndex(2, 7), CartesianIndex(3, 7), CartesianIndex(4, 7), CartesianIndex(5, 7), CartesianIndex(6, 7), CartesianIndex(7, 7)], [(0, 0), (1, 0), (2, 0), (3, 0), (-3, 0), (-2, 0), (-1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (-3, 1), (-2, 1), (-1, 1), (0, 2), (1, 2), (2, 2), (3, 2), (-3, 2), (-2, 2), (-1, 2), (0, 3), (1, 3), (2, 3), (3, 3), (-3, 3), (-2, 3), (-1, 3), (0, -3), (1, -3), (2, -3), (3, -3), (-3, -3), (-2, -3), (-1, -3), (0, -2), (1, -2), (2, -2), (3, -2), (-3, -2), (-2, -2), (-1, -2), (0, -1), (1, -1), (2, -1), (3, -1), (-3, -1), (-2, -1), (-1, -1)]),
        )
    end


    @testset "removeconjfreqs" begin

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.removeconjfreqs(JosephsonCircuits.Frequencies{1}((1,), (2,), (3,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)])),
            JosephsonCircuits.Frequencies{1}((1,), (2,), (3,), CartesianIndex{1}[CartesianIndex(1), CartesianIndex(2)], [(0,), (1,)]),
        )

        frequencies = JosephsonCircuits.Frequencies{2}((2,2), (3, 5), (4, 5), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(3, 1), CartesianIndex(1, 2), CartesianIndex(2, 2), CartesianIndex(3, 2), CartesianIndex(1, 3), CartesianIndex(2, 3), CartesianIndex(3, 3), CartesianIndex(1, 4), CartesianIndex(2, 4), CartesianIndex(3, 4), CartesianIndex(1, 5), CartesianIndex(2, 5), CartesianIndex(3, 5)], [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2), (0, -2), (1, -2), (2, -2), (0, -1), (1, -1), (2, -1)]);
        @test isequal(
            JosephsonCircuits.removeconjfreqs(frequencies).modes,
            [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2), (1, -2), (1, -1)],
        )

        @test isequal(
            JosephsonCircuits.removeconjfreqs(JosephsonCircuits.calcfreqsrdft((2,2))).modes,
            [(0, 0), (1, 0), (2, 0), (0, 1), (1, 1), (2, 1), (0, 2), (1, 2), (2, 2), (1, -2), (1, -1)],
        )

    end

    @testset "keepfreqs" begin

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,2)),[(0,0),(1,0),(0,1),(1,1)]),
            JosephsonCircuits.Frequencies{2}((2, 2), (3, 5), (4, 5), CartesianIndex{2}[CartesianIndex(1, 1), CartesianIndex(2, 1), CartesianIndex(1, 2), CartesianIndex(2, 2)], [(0, 0), (1, 0), (0, 1), (1, 1)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,2)),Tuple{Int64,Int64}[]),
            JosephsonCircuits.Frequencies{2}((2, 2), (3, 5), (4, 5), CartesianIndex{2}[], Tuple{Int64, Int64}[]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[], Tuple{Int64}[]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.keepfreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,)]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1)], [(0,)]),
        )

    end

    @testset "removefreqs" begin

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[(2,)]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,)], [(0,), (1,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[(0,),(1,),(2,),(3,)]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[], Tuple{Int64}[]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),Tuple{Int64}[]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,)]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(2,), CartesianIndex(3,)], [(1,), (2,)]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[CartesianIndex(1,),CartesianIndex(2,),CartesianIndex(3,),CartesianIndex(4,)]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[], Tuple{Int64}[]),
        )

        @test JosephsonCircuits.comparestruct(
            JosephsonCircuits.removefreqs(JosephsonCircuits.calcfreqsrdft((2,)),CartesianIndex{1}[]),
            JosephsonCircuits.Frequencies{1}((2,), (3,), (4,), CartesianIndex{1}[CartesianIndex(1,), CartesianIndex(2,), CartesianIndex(3,)], [(0,), (1,), (2,)]),
        )

    end

    @testset "conjsym" begin

        @test isequal(
            JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsrdft((2,))),
            Dict{CartesianIndex{1}, CartesianIndex{1}}(),
        )

        @test isequal(
            JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsdft((2,))),
            Dict{CartesianIndex{1}, CartesianIndex{1}}(CartesianIndex(3) => CartesianIndex(4), CartesianIndex(2) => CartesianIndex(5)),
        )

        @test isequal(
            JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsrdft((2,1))),
            Dict{CartesianIndex{2}, CartesianIndex{2}}(CartesianIndex(3, 2) => CartesianIndex(3, 3), CartesianIndex(1, 2) => CartesianIndex(1, 3)),
        )

        @test isequal(
            JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsdft((2,1))),
            Dict{CartesianIndex{2}, CartesianIndex{2}}(CartesianIndex(3, 2) => CartesianIndex(4, 3), CartesianIndex(2, 1) => CartesianIndex(5, 1), CartesianIndex(1, 2) => CartesianIndex(1, 3), CartesianIndex(3, 1) => CartesianIndex(4, 1), CartesianIndex(2, 2) => CartesianIndex(5, 3), CartesianIndex(3, 3) => CartesianIndex(4, 2), CartesianIndex(2, 3) => CartesianIndex(5, 2)),
        )

        @test isequal(
            JosephsonCircuits.conjsym(JosephsonCircuits.calcfreqsrdft((2,1,1))),
            Dict{CartesianIndex{3}, CartesianIndex{3}}(CartesianIndex(1, 1, 2) => CartesianIndex(1, 1, 3), CartesianIndex(3, 1, 2) => CartesianIndex(3, 1, 3), CartesianIndex(1, 2, 3) => CartesianIndex(1, 3, 2), CartesianIndex(1, 2, 1) => CartesianIndex(1, 3, 1), CartesianIndex(3, 2, 3) => CartesianIndex(3, 3, 2), CartesianIndex(3, 2, 2) => CartesianIndex(3, 3, 3), CartesianIndex(3, 2, 1) => CartesianIndex(3, 3, 1), CartesianIndex(1, 2, 2) => CartesianIndex(1, 3, 3)),
        )
    end


    @testset "calcindexdict" begin

        @test isequal(
            JosephsonCircuits.calcindexdict(3),
            Dict{CartesianIndex{1}, Int64}(CartesianIndex(3) => 3, CartesianIndex(2) => 2, CartesianIndex(1) => 1),
        )

        @test isequal(
            JosephsonCircuits.calcindexdict((2,3)),
            Dict{CartesianIndex{2}, Int64}(CartesianIndex(1, 1) => 1, CartesianIndex(2, 3) => 6, CartesianIndex(2, 1) => 2, CartesianIndex(1, 2) => 3, CartesianIndex(2, 2) => 4, CartesianIndex(1, 3) => 5),
        )

    end

    @testset "applynl: cos(z*cos(theta))" begin
        # test against Jacobi-Anger expansion for cos(z*cos(theta))
        z = 0.9
        theta = rand()
        mmax = 8
        N = 16
        tolerance = 1e-15

        am = zeros(Complex{Float64},N,1)
        am[2] = z/2*exp(im*(theta)) # z*cos(theta)

        @test isapprox(
            JosephsonCircuits.applynl(am,(x)->cos(x))[1:2:2*mmax],
            # Jacobi-Anger expansion for cos(z*cos(theta))
            [(-1)^m*SpecialFunctions.besselj(2*m,z)*exp(1im*2*m*theta) for m=0:mmax-1],
            atol = tolerance)

        @test isapprox(
            JosephsonCircuits.applynl(am,(x)->cos(x))[2:2:2*mmax],
            zeros(Complex{Float64},mmax),
            atol = tolerance)
    end

    @testset "applynl: sin(z*sin(theta))" begin
        # test against Jacobi-Anger expansion for sin(z*sin(theta))
        z = 0.9
        theta = rand()
        mmax = 8
        N = 16
        tolerance = 1e-15

        am = zeros(Complex{Float64},N,1)
        am[2] = z/2*exp(im*(theta-pi/2)) # z*sin(theta)

        @test isapprox(
            JosephsonCircuits.applynl(am,(x)->sin(x))[2:2:2*mmax],
            # Jacobi-Anger expansion for sin(z*sin(theta))
            [SpecialFunctions.besselj(2*m-1,z)*-im*exp(1im*(2*m-1)*theta) for m=1:mmax],
            atol = tolerance)

        @test isapprox(
            JosephsonCircuits.applynl(am,(x)->sin(x))[1:2:2*mmax],
            zeros(Complex{Float64},mmax),
            atol = tolerance)
    end

    @testset "applynl: x-x^3/6 for z1*cos(theta1)+z2*cos(theta2)" begin
        # test series expansion of sin(z1*cos(theta1)+z2*cos(theta2))
        z1 = 1.1
        z2 = 0.9
        theta1 = rand()
        theta2 = rand()
        N = 8
        tolerance = 1e-15

        am = zeros(Complex{Float64},N,N,1)
        am[2,1,1] = z1/2*exp(im*theta1)
        am[1,2,1] = z2/2*exp(im*theta2)
        am[1,N,1] = z2/2*exp(-im*theta2)

        am2 = JosephsonCircuits.applynl(am,(x)->x-x^3/6)[:,:,1]

        d = Dict(
            CartesianIndex(2,1) => (z1-(3*z1^3/4+3*z1*z2^2/2)/6)*exp(im*theta1)/2, 
            CartesianIndex(1,2) => (z2-(3*z2^3/4+3*z2*z1^2/2)/6)*exp(im*theta2)/2,
            CartesianIndex(4,1) => -z1^3/4/6*exp(im*3*theta1)/2,
            CartesianIndex(1,4) => -z2^3/4/6*exp(im*3*theta2)/2,
            CartesianIndex(2,3) => -3/4*z1*z2^2/6*exp(im*(theta1+2*theta2))/2,
            CartesianIndex(3,2) => -3/4*z1^2*z2/6*exp(im*(2*theta1+theta2))/2,
            CartesianIndex(3,N) => -3/4*z1^2*z2/6*exp(im*(2*theta1-theta2))/2,
            CartesianIndex(2,N-1) => -3/4*z1*z2^2/6*exp(im*(theta1-2*theta2))/2,
            CartesianIndex(1,N) => (z2-(3*z2^3/4+3*z2*z1^2/2)/6)*exp(-im*theta2)/2,
            CartesianIndex(1,N-2) => -z2^3/4/6*exp(-im*3*theta2)/2,
        )
        for i in CartesianIndices(am2)
            if haskey(d,i)
                @test isapprox(am2[i], d[i], atol = tolerance)
            else
                @test isapprox(am2[i], 0, atol = tolerance)
            end
        end
    end

    @testset "applynl: 1-x^2/2 for z1*cos(theta1)+z2*cos(theta2)" begin
        # test series expansion of cos(z1*cos(theta1)+z2*cos(theta2))
        z1 = 1.1
        z2 = 0.9
        theta1 = rand()
        theta2 = rand()
        N = 8
        tolerance = 1e-15

        am = zeros(Complex{Float64},N,N,1)
        am[2,1,1] = z1/2*exp(im*theta1)
        am[1,2,1] = z2/2*exp(im*theta2)
        am[1,N,1] = z2/2*exp(-im*theta2)

        am2 = JosephsonCircuits.applynl(am,(x)->1-x^2/2)[:,:,1]

        d = Dict(
            CartesianIndex(1, 1) => 1-(z1^2/2+z2^2/2)/2,
            CartesianIndex(1, 3) => -((z2^2/2)/2)*exp(im*2*theta2)/2,
            CartesianIndex(3, 1) => -((z1^2/2)/2)*exp(im*2*theta1)/2,
            CartesianIndex(2, 2) => -(z1*z2/2)*exp(im*(theta1+theta2))/2,
            CartesianIndex(2, N) => -(z1*z2/2)*exp(im*(theta1-theta2))/2,
            CartesianIndex(1, N-1) => (-(z2^2/2)/2)*exp(-im*2*theta2)/2,
        )
        for i in CartesianIndices(am2)
            if haskey(d,i)
                @test isapprox(am2[i], d[i], atol = tolerance)
            else
                @test isapprox(am2[i], 0, atol = tolerance)
            end
        end
    end

    @testset "phivectortomatrix! and phimatrixtovector!" begin

        @variables w1,w2
        w = (w1,w2)

        # test whether vector -> matrix -> vector gives the same result
        Nharmonics = (4,3)
        maxintermodorder = 3

        freq = JosephsonCircuits.calcfreqsrdft(Nharmonics);
        truncfreq = JosephsonCircuits.truncfreqs(freq;dc=false,odd=true,even=false,maxintermodorder=maxintermodorder)
        noconjtruncfreq = JosephsonCircuits.removeconjfreqs(truncfreq)
        conjsymdict = JosephsonCircuits.conjsym(noconjtruncfreq)
        freqindexmap, conjsourceindices, conjtargetindices = JosephsonCircuits.calcphiindices(noconjtruncfreq,conjsymdict)
        modes = noconjtruncfreq.modes
        Nt = noconjtruncfreq.Nt
        Nw = noconjtruncfreq.Nw

        Nbranches = 2
        Nfrequencies = length(freqindexmap)

        phivector = rand(Complex{Float64},Nbranches*Nfrequencies);
        phivector1 = rand(Complex{Float64},Nbranches*Nfrequencies);

        phimatrix = zeros(Complex{Float64},(Nw...,Nbranches));

        JosephsonCircuits.phivectortomatrix!(phivector,
            phimatrix,
            freqindexmap,
            conjsourceindices,
            conjtargetindices,
            Nbranches,
        )

        JosephsonCircuits.phimatrixtovector!(phivector1,
            phimatrix,
            freqindexmap,
            conjsourceindices,
            conjtargetindices,
            Nbranches,
        )

        # isapprox(phivector,phivector1)
        @test all(phivector .== phivector1)

    end

    @testset "phimatrixtovector!"  begin

        begin
            freqindexmap = [2, 4, 6, 8, 12, 16, 27, 33]
            conjsourceindices = [16, 6]
            conjtargetindices = [21, 31]
            Nbranches = 1

            phivector = zeros(Complex{Float64}, Nbranches*length(freqindexmap)-1)
            phimatrix = [0.0 + 0.0im 0.0 + 3.0im 0.0 + 0.0im 0.0 + 6.0im 0.0 - 6.0im 0.0 + 0.0im 0.0 - 3.0im; 0.0 + 1.0im 0.0 + 0.0im 0.0 + 5.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 7.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 4.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 8.0im; 0.0 + 2.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;;;]

            @test_throws(
                DimensionMismatch("Unexpected length for phivector"),
                JosephsonCircuits.phimatrixtovector!(phivector,
                    phimatrix,
                    freqindexmap,
                    conjsourceindices,
                    conjtargetindices,
                    Nbranches,
                )
            )
        end

        begin
            freqindexmap = [2, 4, 6, 8, 12, 16, 27, 33]
            conjsourceindices = [16, 6]
            conjtargetindices = [21, 31]
            Nbranches = 1

            phimatrix = [0.0 + 0.0im 0.0 + 3.0im 0.0 + 0.0im 0.0 + 6.0im 0.0 - 6.0im 0.0 + 0.0im 0.0 - 3.0im; 0.0 + 1.0im 0.0 + 0.0im 0.0 + 5.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 7.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 4.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 8.0im; 0.0 + 2.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;;;]
            phivector = 1im.*Complex.(1:Nbranches*length(freqindexmap));

            phimatrix1 = similar(phimatrix)
            phivector1 = similar(phivector)

            JosephsonCircuits.phimatrixtovector!(phivector1,
                phimatrix,
                freqindexmap,
                conjsourceindices,
                conjtargetindices,
                Nbranches,
            )

            @test all(phivector .== phivector1)
        end
    end



    @testset "phivectortomatrix!"  begin
        begin
            freqindexmap = [2, 4, 6, 8, 12, 16, 27, 33]
            conjsourceindices = [16, 6]
            conjtargetindices = [21, 31]
            Nbranches = 1

            phimatrix = [0.0 + 0.0im 0.0 + 3.0im 0.0 + 0.0im 0.0 + 6.0im 0.0 - 6.0im 0.0 + 0.0im 0.0 - 3.0im; 0.0 + 1.0im 0.0 + 0.0im 0.0 + 5.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 7.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 4.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 8.0im; 0.0 + 2.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im;;;]
            phivector = 1im.*Complex.(1:Nbranches*length(freqindexmap));

            phimatrix1 = similar(phimatrix)
            phivector1 = similar(phivector)

            JosephsonCircuits.phivectortomatrix!(phivector,
                phimatrix1,
                freqindexmap,
                conjsourceindices,
                conjtargetindices,
                Nbranches,
            )

            @test all(phimatrix .== phimatrix1)
        end

        begin
            freqindexmap = [2, 4, 6, 8, 12, 16, 27, 33]
            conjsourceindices = [16, 6]
            conjtargetindices = [21, 31]
            Nbranches = 1

            phivector = 1im.*Complex.(1:(Nbranches*length(freqindexmap)-1));
            phimatrix=zeros(Complex{Float64},5,7,1)

            @test_throws(
                DimensionMismatch("Unexpected length for phivector"),
                JosephsonCircuits.phivectortomatrix!(phivector,
                    phimatrix,
                    freqindexmap,
                    conjsourceindices,
                    conjtargetindices,
                    Nbranches,
                )
            )
        end
    end

end
