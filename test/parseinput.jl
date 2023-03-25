using JosephsonCircuits
using Test

@testset verbose=true "parseinput" begin

    @testset "parsecircuit" begin

        circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
        push!(circuit,("P1","1","0",1))
        push!(circuit,("I1","1","0",:Ipump))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("L1","1","0",:L1))
        push!(circuit,("L1","1","0",:L1))
        push!(circuit,("K1","L1","L2",:K1))
        push!(circuit,("L2","2","0",:L2))
        push!(circuit,("C2","2","0",:C2))
        @test_throws(
            ArgumentError("""Name "L1" on line 5 is not unique."""),
            parsecircuit(circuit)
        )

    end

    @testset "parsecomponenttype" begin
        @test_throws(
            ArgumentError("parsecomponenttype() currently only works for two letter components"),
            JosephsonCircuits.parsecomponenttype("BAD1",["Lj","BAD","L","C","K","I","R","P"])
        )
    end

    @testset "parsecomponenttype" begin
        @test_throws(
            ArgumentError("No matching component found in allowedcomponents."),
            JosephsonCircuits.parsecomponenttype("B1",["Lj","L","C","K","I","R","P"])
        )
    end

    @testset "checkcomponenttypes" begin
        @test_throws(
            ArgumentError("Allowed components parsing check has failed for Lj. This can happen if a two letter long component comes after a one letter component. Please reorder allowedcomponents."),
            JosephsonCircuits.checkcomponenttypes(["L","Lj","C","K","I","R","P"])
        )
    end

    @testset "extractbranches" begin
        @test_throws(
            DimensionMismatch("componenttypes must have the same length as the number of node indices"),
            JosephsonCircuits.extractbranches(
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3; 1 1 1 3 1]
            )
        )
    end

    @testset "extractbranches" begin
        @test_throws(
            DimensionMismatch("the length of the first axis must be 2"),
            JosephsonCircuits.extractbranches(
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3 3; 1 1 1 3 1 1; 0 0 0 0 0 0],
            )
        )
    end

    @testset "extractbranches!" begin
        @test_throws(
            DimensionMismatch("branchvector should be length zero"),
            JosephsonCircuits.extractbranches!(
                [1],
                [:P,:I,:R,:C,:Lj,:C],
                [2 2 2 2 3 3; 1 1 1 3 1 1],
            )
        )
    end

    @testset "calcnodesorting" begin
        @test_throws(
            ArgumentError("Unknown sorting type."),
            JosephsonCircuits.calcnodesorting(["30","11","0","2"];sorting=:test)
        )
    end

    @testset "calcnodesorting" begin
        @test_throws(
            ArgumentError("No ground node found in netlist."),
            JosephsonCircuits.calcnodesorting(["30","11","1","2"];sorting=:none)
        )
    end

    @testset "calcnodesorting" begin
        @test_throws(
            ArgumentError("No ground node found in netlist."),
            JosephsonCircuits.calcnodesorting(String[];sorting=:none)
        )
    end

    @testset "calcnodesorting" begin
        @test_throws(
            ArgumentError("Failed to parse the nodes as integers. Try setting the keyword argument `sorting=:name` or `sorting=:none`."),
            JosephsonCircuits.calcnodesorting(["30","11","0","a"];sorting=:number)
        )
    end

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

    @testset "calcnoiseportimpedanceindices" begin
        @test_throws(
            DimensionMismatch("Input arrays must have the same length"),
            JosephsonCircuits.calcnoiseportimpedanceindices(
                [:R,:C,:Lj],
                [2 2 3 3; 1 3 1 1],
                [],
                [50,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcnoiseportimpedanceindices" begin
        @test_throws(
            DimensionMismatch("The length of the first axis must be 2"),
            JosephsonCircuits.calcnoiseportimpedanceindices(
                [:R,:C,:Lj,:C],
                [2 2 3 3; 1 3 1 1; 0 0 0 0],
                [],
                [50,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcportindicesnumbers" begin
        @test_throws(
            DimensionMismatch("The length of the first axis must be 2"),
            JosephsonCircuits.calcportindicesnumbers(
                [:P,:R,:C,:Lj,:C],
                [2 2 2 3 3; 1 1 3 1 1; 0 0 0 0 0],
                [],
                [1,50,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcportindicesnumbers" begin
        @test_throws(
            DimensionMismatch("Input arrays must have the same length"),
            JosephsonCircuits.calcportindicesnumbers(
                [:P,:R,:C,:Lj,:C],
                [2 2 2 3; 1 1 3 1],
                [],
                [1,50,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcportindicesnumbers" begin
        @test_throws(
            ArgumentError("Only one port allowed per branch."),
            JosephsonCircuits.calcportindicesnumbers(
                [:P,:P,:C,:Lj,:C],
                [2 2 2 3 3; 1 1 3 1 1],
                [],
                [1,2,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcportindicesnumbers" begin
        @test_throws(
            ArgumentError("Duplicate ports are not allowed."),
            JosephsonCircuits.calcportindicesnumbers(
                [:P,:R,:C,:Lj,:P],
                [2 2 2 3 3; 1 1 3 1 1],
                [],
                [1,50,5e-15,1e-12,1]
            )
        )
    end

    @testset "calcportimpedanceindices" begin
        @test_throws(
            DimensionMismatch("The length of the first axis must be 2"),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:C,:Lj,:C],
                [2 2 2 3 3; 1 1 3 1 1; 0 0 0 0 0],
                [],
                [1,50,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcportimpedanceindices" begin
        @test_throws(
            DimensionMismatch("Input arrays must have the same length"),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:C,:Lj,:C],
                [2 2 2 3; 1 1 3 1],
                [],
                [1,50,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcportimpedanceindices" begin
        @test_throws(
            ArgumentError("Only one port allowed per branch."),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:P,:C,:Lj,:C],
                [2 2 2 3 3; 1 1 3 1 1],
                [],
                [1,2,5e-15,1e-12,30e-15]
            )
        )
    end

    @testset "calcportimpedanceindices" begin
        @test_throws(
            ArgumentError("Duplicate ports are not allowed."),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:C,:Lj,:P],
                [2 2 2 3 3; 1 1 3 1 1],
                [],
                [1,50,5e-15,1e-12,1]
            )
        )
    end

    @testset "calcportimpedanceindices" begin
        @test_throws(
            ArgumentError("Only one resistor allowed per port."),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:R,:Lj,:C],
                [2 2 2 2 3; 1 1 1 3 1],
                [],
                [1,50.0,50.0,1e-12,30e-15]
            )
        )
    end

end