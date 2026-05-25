using JosephsonCircuits


@testset verbose=true "unwrap" begin

    @testset "unwrap!" begin

        theta=0:0.01:6*pi
        @test isapprox(theta,JosephsonCircuits.unwrap(angle.(exp.(im*theta))))

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