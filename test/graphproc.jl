using JosephsonCircuits
using Test

@testset verbose=true "graphproc" begin
    @testset "storeuniqueloops!" begin
        lvarray = Vector{Int}[];
        @test_throws(
            "There should only be one loop associated with each closure branch.",
            JosephsonCircuits.storeuniqueloops!(lvarray,[1, 2, 3],[[1,2,3],[4,5,6]])
        )
    end

    
end