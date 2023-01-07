using JosephsonCircuits
using Test
using SpecialFunctions

@testset verbose=true "fftutils" begin

    @testset "calcfrequencies" begin

        @variables w1,w2
        w = (w1,w2)
        Nharmonics = (2,)
        maxintermodorder = 2
        @test_throws(
            ErrorException("Each frequency must have a number of harmonics."),
            JosephsonCircuits.calcfrequencies(w,Nharmonics,
                maxintermodorder=maxintermodorder,dc=false,even=false,odd=true)
        )
    end

    @testset "vectortodense" begin
        @test_throws(
            DimensionMismatch("Dimensions of coords elements and Nharmonics must be consistent."),
            JosephsonCircuits.vectortodense([CartesianIndex(1,1)],[1],(1,))
        )
    end

    @testset "vectortodense" begin
        @test_throws(
            DimensionMismatch("Not designed to visualize higher dimensional arrays"),
            JosephsonCircuits.vectortodense([CartesianIndex(1,1,1)],[1],(1,1,1))
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

        Nw,coords,values,dropcoords,dropvalues = JosephsonCircuits.calcfrequencies(w,Nharmonics,
            maxintermodorder=maxintermodorder,dc=false,even=false,odd=true)
        Nt=NTuple{length(Nw),Int64}(ifelse(i == 1, 2*val-1, val) for (i,val) in enumerate(Nw))

        dropdict = Dict(dropcoords .=> dropvalues)

        freqindexmap,conjsourceindices,conjtargetindices = JosephsonCircuits.calcphiindices(Nt,dropdict)

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

    @testset "phimatrixtovector!"  begin

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



    @testset "phivectortomatrix!"  begin

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

    @testset "phivectortomatrix!"  begin

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
