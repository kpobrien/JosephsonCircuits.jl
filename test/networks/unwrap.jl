using JosephsonCircuits


@testset verbose=true "unwrap" begin

    @testset "unwrap!" begin

        theta=0:0.01:6*pi
        @test isapprox(theta,JosephsonCircuits.unwrap(angle.(exp.(im*theta))))

        # complex values are not phases
        @test_throws(
            ArgumentError,
            JosephsonCircuits.unwrap(exp.(im*theta)),
        )
        @test_throws(
            ArgumentError,
            JosephsonCircuits.unwrap!(exp.(im*theta)),
        )

        # integer values unwrap as floating point numbers
        @test JosephsonCircuits.unwrap([0, 7, 1]) ≈ [0, 7 - 2pi, 1]

        # a keyword it does not take, such as `dim` for `dims`, is refused
        @test_throws MethodError JosephsonCircuits.unwrap([0.0, 7.0]; dim = 1)

        @test isapprox(
            JosephsonCircuits.unwrap!(zeros(10,10);dims=1),
            zeros(10,10),
        )

        @test_throws(
            ArgumentError(lazy"`unwrap!`: required keyword parameter dims missing"),
            JosephsonCircuits.unwrap!(zeros(10,10),zeros(10,10)),
        )

        @test_throws(
            ArgumentError(lazy"`unwrap!`: Invalid dims specified: a"),
            JosephsonCircuits.unwrap!(zeros(10),zeros(10);dims='a'),
        )

    end

end