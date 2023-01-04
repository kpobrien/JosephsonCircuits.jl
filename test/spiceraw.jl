using JosephsonCircuits
using Test

@testset verbose=true "spiceraw" begin

    @testset "spice_raw_load" begin

        filepath = joinpath(dirname(Base.source_path()),"spiceraw","test01.raw")

        out1 = JosephsonCircuits.spice_raw_load(filepath)

        out2 = (variables = Dict("V" => ["v(1)", "v(2)", "v(3)"], "Hz" => ["frequency"]), values = Dict{Any, Any}("V" => ComplexF64[48.87562301047733 - 7.413126995337487im 49.97131616467212 + 1.1949290155299537im 49.02611690128596 - 6.90980805243651im; -10.116167243319213 + 1.534380793728424im 57.578470543293086 + 1.3775359827006193im 12.368446655904192 - 1.743197747303436im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im], "Hz" => ComplexF64[4.0e9 + 0.0im 5.0e9 + 0.0im 6.0e9 + 0.0im]))

        @test out1.variables == out2.variables
        @test out1.values == out2.values
    end

    @testset "spice_raw_load option line" begin

        filepath = joinpath(dirname(Base.source_path()),"spiceraw","test02.raw")

        out1 = JosephsonCircuits.spice_raw_load(filepath)

        # out2 = (variables = Dict("V" => ["v(1)", "v(2)", "v(3)"], "Hz" => ["frequency"]), values = Dict{Any, Any}("V" => ComplexF64[48.87562301047733 - 7.413126995337487im 49.97131616467212 + 1.1949290155299537im 49.02611690128596 - 6.90980805243651im; -10.116167243319213 + 1.534380793728424im 57.578470543293086 + 1.3775359827006193im 12.368446655904192 - 1.743197747303436im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im], "Hz" => ComplexF64[4.0e9 + 0.0im 5.0e9 + 0.0im 6.0e9 + 0.0im]))

        # @test out1.variables == out2.variables
        # @test out1.values == out2.values
    end

    @testset "spice_raw_load extra variables" begin

        filepath = joinpath(dirname(Base.source_path()),"spiceraw","test03.raw")

        out1 = JosephsonCircuits.spice_raw_load(filepath)

        out2 = (variables = Dict("V" => ["v(1)", "v(2)", "v(3)"], "Hz" => ["frequency"]), values = Dict{Any, Any}("V" => ComplexF64[48.87562301047733 - 7.413126995337487im 49.97131616467212 + 1.1949290155299537im 49.02611690128596 - 6.90980805243651im; -10.116167243319213 + 1.534380793728424im 57.578470543293086 + 1.3775359827006193im 12.368446655904192 - 1.743197747303436im; 0.0 + 0.0im 0.0 + 0.0im 0.0 + 0.0im], "Hz" => ComplexF64[4.0e9 + 0.0im 5.0e9 + 0.0im 6.0e9 + 0.0im]))

        @test out1.variables == out2.variables
        @test out1.values == out2.values
    end

    @testset "spice_raw_load ascii" begin

        filepath = joinpath(dirname(Base.source_path()),"spiceraw","test04.raw")

        @test_throws "This function only handles Binary files not ASCII" JosephsonCircuits.spice_raw_load(filepath)

    end

    @testset "spice_raw_load flags" begin

        filepath = joinpath(dirname(Base.source_path()),"spiceraw","test05.raw")

        @test_throws "Unknown flag" JosephsonCircuits.spice_raw_load(filepath)

    end

end