using JosephsonCircuits
using Test

@testset verbose=true "capindmat" begin

    @testset "calcMb JJ as first inductor" begin
        Nmodes = 2
        Nbranches = 2
        componenttypes = [:Lj,:K,:L,:C]
        nodeindices = [2 0 3 3; 1 0 1 1]
        componentvalues = [1.0e-9, 0.1, 4.0e-9, 2.0e-12]
        componentnamedict = Dict{Symbol, Int}(:C1 => 4,:L2 => 3,:Lj1 => 1,:K1 => 2)
        edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
        mutualinductorbranchnames = [ :Lj1, :L2]

        @test_throws(
            ArgumentError("Mutual coupling coefficient K must couple two inductors. Lj1 is not an inductor."),
            JosephsonCircuits.calcMb(componenttypes,nodeindices,componentvalues,componentnamedict,mutualinductorbranchnames,edge2indexdict,Nmodes,Nbranches)
        )
    end

    @testset "calcMb JJ as second inductor" begin
        Nmodes = 2
        Nbranches = 2
        componenttypes = [:L,:K,:Lj,:C]
        nodeindices = [2 0 3 3; 1 0 1 1]
        componentvalues = [1.0e-9, 0.1, 4.0e-9, 2.0e-12]
        componentnamedict = Dict{Symbol, Int}(:C1 => 4,:Lj2 => 3,:L1 => 1,:K1 => 2)
        edge2indexdict = Dict{Tuple{Int, Int}, Int}((1, 2) => 1,(3, 1) => 2,(1, 3) => 2,(2, 1) => 1)
        mutualinductorbranchnames = [ :L1, :Lj2]

        @test_throws(
            ArgumentError("Mutual coupling coefficient K must couple two inductors. Lj2 is not an inductor."),
            JosephsonCircuits.calcMb(componenttypes,nodeindices,componentvalues,componentnamedict,mutualinductorbranchnames,edge2indexdict,Nmodes,Nbranches)
        )
    end

    @testset "calcLmean_inner" begin
        @test_throws(
            DimensionMismatch("componenttypes and componentvalues should have the same length"),
            JosephsonCircuits.calcLmean_inner([:L,:C,:Lj],[10,4,5,1],Float64[])
        )
    end

    @testset "calcnodematrix" begin
        @test_throws(
            DimensionMismatch("nodeindices should have a first dimension size of 2."),
            JosephsonCircuits.calcnodematrix(
                [:R,:R],[2 3;1 1;0 0],[1.0,2.0],Float64[],1,3,:R,false)
        )
        @test_throws(
            DimensionMismatch("componenttypes, nodeindices, and componentvalues should have the same length"),
            JosephsonCircuits.calcnodematrix([:R],[2 3;1 1],[1.0,2.0],
                Float64[],1,3,:R,false)
        )
    end

    @testset "combine" begin

        a = rand()
        b = rand()
        @test(JosephsonCircuits.combine_sum(a,b) == a+b)
        @test(JosephsonCircuits.combine_reciprocal_sum(a,b) == a*b/(a+b))
        @test_throws(
            ArgumentError("Components 1 and 2 cannot be combined to a single element. Please place the two components between different nodes."),
            JosephsonCircuits.combine_error(1,2),
        )
    end

    # the element type the matrix builders assemble in
    @testset "calcvaluetype" begin
        @test_throws(
            DimensionMismatch("componenttypes and componentvalues should have the same length"),
            JosephsonCircuits.calcvaluetype(
                [:C,:R],
                [1,2,3],
                [:R]
            )
        )
    end
end
