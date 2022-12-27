using JosephsonCircuits
using Test
using SpecialFunctions

@testset verbose = true "fftutils" begin

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

end
