# Circuits shared across test files, so a change to a canonical test
# device is made in one place. Files include this defensively
# (`isdefined(Main, ...) || include(...)`) so each test file still runs
# standalone.

# a JPA: one port, one junction, a coupling capacitor -- symbolic values
# with a definitions dictionary, the standard user-facing form
function testjpacircuit()
    circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
    push!(circuit,("P1","1","0",1))
    push!(circuit,("R1","1","0",:Rleft))
    push!(circuit,("C1","1","2",:Cc))
    push!(circuit,("Lj1","2","0",:Lj))
    push!(circuit,("C2","2","0",:Cj))
    circuitdefs = Dict(
        :Lj => 1000.0e-12,
        :Cc => 100.0e-15,
        :Cj => 1000.0e-15,
        :Rleft => 50.0,
    )
    return circuit, circuitdefs
end

# the same JPA with numeric literal values and an empty definitions
# dictionary, which exercises the numeric-value parsing path
function testjpacircuitnumeric()
    circuit = [("P1","1","0",1), ("R1","1","0",50.0), ("C1","1","2",100e-15),
               ("Lj1","2","0",1000e-12), ("C2","2","0",1000e-15)]
    return circuit, Dict{Any,Any}()
end

# a four junction transmission line chain with a port at each end
function testchaincircuit()
    circuit = Tuple{String,String,String,Union{Complex{Float64},Symbol,Int64}}[]
    push!(circuit, ("P1","1","0",1)); push!(circuit, ("R1","1","0",:R))
    for i in 1:4
        push!(circuit, ("Lj$(i)","$(i)","$(i+1)",:Lj))
        push!(circuit, ("C$(i)","$(i)","0",:Cg))
    end
    push!(circuit, ("C5","5","0",:Cg)); push!(circuit, ("R2","5","0",:R))
    circuitdefs = Dict{Symbol,Complex{Float64}}(
        :Lj => 100e-12, :Cg => 40e-15, :R => 50.0)
    return circuit, circuitdefs
end
