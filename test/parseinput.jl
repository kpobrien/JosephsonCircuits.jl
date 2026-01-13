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

        @variables Ipump Rleft L1 L2 C2
        Kfun(L) = sin(L);@register_symbolic Kfun(L1)
        circuit = Vector{Tuple{String,String,String,Num}}(undef,0)
        push!(circuit,("P1","1","0",1))
        push!(circuit,("I1","1","0",Ipump))
        push!(circuit,("R1","1","0",Rleft))
        push!(circuit,("L1","1","0",L1))
        push!(circuit,("K1","L1","L2",Kfun(L1)))
        push!(circuit,("L2","2","0",L2))
        push!(circuit,("C2","2","0",C2))
        @test JosephsonCircuits.comparestruct(parsecircuit(circuit),JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2], ["1", "0", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Num[1, Ipump, Rleft, L1, Kfun(L1), L2, C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3))


        circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
        push!(circuit,("P1","1","0",1))
        push!(circuit,("I1","1","0",:Ipump))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        @test JosephsonCircuits.comparestruct(parsecircuit(circuit),JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 3, 3, 2, 3, 2], ["1", "0", "2"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Union{Int64, Symbol, ComplexF64}[1, :Ipump, :Rleft, :Cc, :Lj, :Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3))

        circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
        push!(circuit,("P1","One","0",1))
        push!(circuit,("I1","One","0",:Ipump))
        push!(circuit,("R1","One","0",:Rleft))
        push!(circuit,("C1","One","Two",:Cc))
        push!(circuit,("Lj1","Two","0",:Lj))
        push!(circuit,("C2","Two","0",:Cj))
        @test JosephsonCircuits.comparestruct(parsecircuit(circuit),JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 3, 3, 2, 3, 2], ["One", "0", "Two"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Union{Int64, Symbol, ComplexF64}[1, :Ipump, :Rleft, :Cc, :Lj, :Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3))


        circuit = []
        push!(circuit,("P1","1","0",1))
        push!(circuit,("I1","1","0",:Ipump))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("C1","1","2",:Cc))
        push!(circuit,("Lj1","2","0",:Lj))
        push!(circuit,("C2","2","0",:Cj))
        @test JosephsonCircuits.comparestruct(parsecircuit(circuit),JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 3, 3, 2, 3, 2], ["1", "0", "2"], String[], ["P1", "I1", "R1", "C1", "Lj1", "C2"], [:P, :I, :R, :C, :Lj, :C], Any[1, :Ipump, :Rleft, :Cc, :Lj, :Cj], Dict("I1" => 2, "C1" => 4, "C2" => 6, "R1" => 3, "P1" => 1, "Lj1" => 5), 3))


        circuit = Vector{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}}}(undef,0)
        push!(circuit,("P1","1","0",1))
        push!(circuit,("I1","1","0",:Ipump))
        push!(circuit,("R1","1","0",:Rleft))
        push!(circuit,("L1","1","0",:L1))
        push!(circuit,("K1","L1","L2",:K1))
        push!(circuit,("L2","2","0",:L2))
        push!(circuit,("C2","2","0",:C2))
        @test JosephsonCircuits.comparestruct(parsecircuit(circuit),JosephsonCircuits.ParsedCircuit([1, 2, 1, 2, 1, 2, 1, 2, 0, 0, 3, 2, 3, 2], ["1", "0", "2"], ["L1", "L2"], ["P1", "I1", "R1", "L1", "K1", "L2", "C2"], [:P, :I, :R, :L, :K, :L, :C], Union{Int64, Symbol, ComplexF64}[1, :Ipump, :Rleft, :L1, :K1, :L2, :C2], Dict("L1" => 4, "I1" => 2, "L2" => 6, "C2" => 7, "R1" => 3, "P1" => 1, "K1" => 5), 3))


    end

    @testset "parsecomponenttype" begin
        @test_throws(
            ArgumentError("parsecomponenttype() currently only works for two letter components"),
            JosephsonCircuits.parsecomponenttype("BAD1",["Lj","BAD","L","C","K","I","R","P"])
        )

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

        @test_throws(
            ArgumentError("No ground node found in netlist."),
            JosephsonCircuits.calcnodesorting(["30","11","1","2"];sorting=:none)
        )

        @test_throws(
            ArgumentError("No ground node found in netlist."),
            JosephsonCircuits.calcnodesorting(String[];sorting=:none)
        )

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

        @test_throws(
            DimensionMismatch("Input arrays must have the same length"),
            JosephsonCircuits.calcportindicesnumbers(
                [:P,:R,:C,:Lj,:C],
                [2 2 2 3; 1 1 3 1],
                [],
                [1,50,5e-15,1e-12,30e-15]
            )
        )

        @test_throws(
            ArgumentError("Only one port allowed per branch."),
            JosephsonCircuits.calcportindicesnumbers(
                [:P,:P,:C,:Lj,:C],
                [2 2 2 3 3; 1 1 3 1 1],
                [],
                [1,2,5e-15,1e-12,30e-15]
            )
        )

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

        @test_throws(
            DimensionMismatch("Input arrays must have the same length"),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:C,:Lj,:C],
                [2 2 2 3; 1 1 3 1],
                [],
                [1,50,5e-15,1e-12,30e-15]
            )
        )

        @test_throws(
            ArgumentError("Only one port allowed per branch."),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:P,:C,:Lj,:C],
                [2 2 2 3 3; 1 1 3 1 1],
                [],
                [1,2,5e-15,1e-12,30e-15]
            )
        )

        @test_throws(
            ArgumentError("Duplicate ports are not allowed."),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:C,:Lj,:P],
                [2 2 2 3 3; 1 1 3 1 1],
                [],
                [1,50,5e-15,1e-12,1]
            )
        )

        @test_throws(
            ArgumentError("Only one resistor allowed per port."),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:R,:Lj,:C],
                [2 2 2 2 3; 1 1 1 3 1],
                [],
                [1,50.0,50.0,1e-12,30e-15]
            )
        )

        @test_throws(
            ArgumentError("Ports without resistors detected. Each port must have a resistor to define the impedance."),
            JosephsonCircuits.calcportimpedanceindices(
                [:P,:R,:C,:Lj,:C,:P,:C],
                [2 2 2 3 3 3 3; 1 1 3 1 1 1 1],
                [],
                [2,50,5e-15,1e-12,30e-15,1,50.0]
            )
        )

    end

    @testset "Symbolics.jl and complex numbers" begin
    # this functionality is currently broken. this test will error once it is
    # fixed and then we can revise valuetonumber

        @variables Lj1
        @test(JosephsonCircuits.Symbolics.evaluate(Lj1,Dict(Lj1=>3.0e-12))==3.0e-12)

        @variables Lj1::Complex
        @test(JosephsonCircuits.Symbolics.evaluate(Lj1,Dict(Lj1=>3.0e-12))==3.0e-12,broken=true)

        @variables Lj1 Lj2
        @test(JosephsonCircuits.Symbolics.evaluate(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12)) == 4.0e-12)

        @variables Lj1::Complex Lj2::Complex
        @test(JosephsonCircuits.Symbolics.evaluate(Lj1+Lj2,Dict(Lj1=>3.0e-12,Lj2=>1.0e-12)) == 4.0e-12,broken=true)
    end
end