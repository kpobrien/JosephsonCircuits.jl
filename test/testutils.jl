using JosephsonCircuits
using Test

@testset verbose=true "testutils" begin

    @testset "compare" begin

        @test JosephsonCircuits.compare(
            JosephsonCircuits.PassiveNetwork("S1",zeros(2,2),zeros(2,2),[("S1",1),("S1",2)]),
            JosephsonCircuits.PassiveNetwork("S1",zeros(2,2),zeros(2,2),[("S1",1),("S1",2)]),
        )

    end
end