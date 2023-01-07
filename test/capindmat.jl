using JosephsonCircuits
using Test

@testset verbose=true "capindmat" begin
    @testset "calcLmean_inner" begin
        @test_throws(
            "DimensionMismatch: typevector and valuevector should have the same length",
            JosephsonCircuits.calcLmean_inner([:L,:C,:Lj],[10,4,5,1],Float64[])
        )
    end

    @testset "calcnodematrix" begin
        @test_throws(
            "DimensionMismatch: nodeindexarray should have a first dimension size of 2.",
            JosephsonCircuits.calcnodematrix(
                [:R,:R],[2 3;1 1;0 0],[1.0,2.0],Float64[],1,3,:R,false)
        )
        @test_throws(
            "DimensionMismatch: typevector, nodeindexarray, and valuevector should have the same length",
            JosephsonCircuits.calcnodematrix([:R],[2 3;1 1],[1.0,2.0],
                Float64[],1,3,:R,false)
        )
    end

end