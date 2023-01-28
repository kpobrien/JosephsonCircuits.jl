__precompile__(true)

module JosephsonCircuits

import Graphs
import FFTW
import KLU
import UUIDs
import Symbolics: Symbolic, Sym

using Symbolics
using OrderedCollections
using LinearAlgebra
using SparseArrays
using Printf
using SnoopPrecompile

# define the zero for symbolic numbers so that we can view the sparse arrays
# Base.zero(::Type{Symbolic{Number}}) = 0
# Base.zero(::Type{Any}) = 0
# Base.zero(::Type{Nothing}) = 0


# define the reduced flux quantum in Weber, H*A
# hbar/(2*charge of electron)
const phi0 = 3.29105976e-16

# The file structure below mimics the typical analysis flow. An input file
# is first parsed, then the incidence matrix is calculated, then the 
# capacitance, inverse inductance, and other matrices are calculated, then
# then circuit is solved using the harmonic balance method, and postprocessed
# to determine the scattering parameters and quantum efficiency. 
include("parseinput.jl")
include("graphproc.jl")
include("capindmat.jl")
include("matutils.jl")
include("nlsolve.jl")
include("hbsolve.jl")
include("fftutils.jl")
include("qesparams.jl")

# Experimental multi-tone harmonic balance code.
include("hbsolve2.jl")

# These are for exporting SPICE netlists and running simulations in
# WRSPICE or Xyce. 
include("exportnetlist.jl")
include("spiceutils.jl")
include("spicewrapper.jl")
include("spiceraw.jl")

# This is for reading and writing files in the touchstone format. 
include("touchstone.jl")

# These are functions for perofrming tests
include("testutils.jl")


"""
    LjtoIc(Lj)

Convert the junction inductance to critical current in SI base units. 

# Examples
```jldoctest
julia> LjtoIc(100e-12)
3.29105976e-6
```
"""
function LjtoIc(Lj)
    #reduced flux quantum in Weber, H*A
    # hbar/(2*charge of electron)
    # phi0 = 3.29105976e-16
    return phi0./Lj
end

"""
    IctoLj(Ic)

Convert the junction critical current to inductance in SI base units. 

# Examples
```jldoctest
julia> IctoLj(3.29105976e-6)
1.0e-10
```
"""
function IctoLj(Ic)
    # the formula is the same for the two conversions
    return LjtoIc(Ic)
end

function warmup()

    # define the circuit components
    circuit = Array{Tuple{String,String,String,Union{Complex{Float64}, Symbol,Int}},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
    push!(circuit,("I1","1","0",:Ipump))
    push!(circuit,("R1","1","0",:Rleft))

    push!(circuit,("C1","1","2",:Cc)) 
    push!(circuit,("Lj1","2","0",:Lj)) 
    push!(circuit,("C2","2","0",:Cj))

    circuitdefs = Dict{Symbol,Complex{Float64}}(
        :Lj =>1000.0e-12,
        :Cc => 100.0e-15,
        :Cj => 1000.0e-15,
        :Rleft => 50.0,
        :Ipump => 1.0e-8,
    )

    return hbsolve(2*pi*(4.5:0.5:5.0)*1e9,2*pi*4.75001*1e9,0.00565e-6,2,2,circuit,circuitdefs,pumpports=[1],solver=:klu);

end

function warmupsyms()

    @variables Rleft Cc Lj Cj w L1
    circuit = Tuple{String,String,String,Num}[]
    push!(circuit,("P1","1","0",1))
    push!(circuit,("R1","1","0",Rleft))
    push!(circuit,("C1","1","2",Cc)) 
    push!(circuit,("Lj1","2","0",Lj)) 
    push!(circuit,("C2","2","0",Cj))
    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        Rleft => 50.0,
    )

    return hbsolve(2*pi*(4.5:0.5:5.0)*1e9, 2*pi*4.75001*1e9, 0.00565e-6, 2, 2, circuit, circuitdefs, pumpports=[1]);

end

function warmupsyms2()

    @variables Rleft Cc Lj Cj w L1
    circuit = Tuple{String,String,String,Num}[]
    push!(circuit,("P1","1","0",1))
    push!(circuit,("R1","1","0",Rleft))
    push!(circuit,("C1","1","2",Cc)) 
    push!(circuit,("Lj1","2","0",Lj)) 
    push!(circuit,("C2","2","0",Cj))
    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        Rleft => 50.0,
    )

    return hbsolve2(2*pi*(4.5:0.5:5.0)*1e9, 2*pi*4.75001*1e9, 0.00565e-6, 2, 2, circuit, circuitdefs, pumpports=[1]);

end

function warmupparse()

    @variables Ipump Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
    push!(circuit,("I1","1","0",Ipump))
    # push!(circuit,("R1","1","0",Rleft*w/(w+1)))
    push!(circuit,("R1","1","0",Rleft))

    push!(circuit,("C1","1","2",Cc)) 
    push!(circuit,("Lj1","2","0",Lj)) 
    push!(circuit,("C2","2","0",Cj))

    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        Rleft => 50.0,
        Ipump => 1.0e-8,
    )

    return parsecircuit(circuit)
end

function warmupparsesort()

    @variables Ipump Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
    push!(circuit,("I1","1","0",Ipump))
    # push!(circuit,("R1","1","0",Rleft*w/(w+1)))
    push!(circuit,("R1","1","0",Rleft))

    push!(circuit,("C1","1","2",Cc)) 
    push!(circuit,("Lj1","2","0",Lj)) 
    push!(circuit,("C2","2","0",Cj))

    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        Rleft => 50.0,
        Ipump => 1.0e-8,
    )

    return parsesortcircuit(circuit)
end

function warmupnumericmatrices()

    @variables Ipump Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
    push!(circuit,("I1","1","0",Ipump))
    # push!(circuit,("R1","1","0",Rleft*w/(w+1)))
    push!(circuit,("R1","1","0",Rleft))

    push!(circuit,("C1","1","2",Cc)) 
    push!(circuit,("Lj1","2","0",Lj)) 
    push!(circuit,("C2","2","0",Cj))

    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        Rleft => 50.0,
        Ipump => 1.0e-8,
    )

    return numericmatrices(circuit,circuitdefs)
end

function warmuphblinsolve()

    @variables Ipump Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
    push!(circuit,("I1","1","0",Ipump))
    # push!(circuit,("R1","1","0",Rleft*w/(w+1)))
    push!(circuit,("R1","1","0",Rleft))

    push!(circuit,("C1","1","2",Cc)) 
    push!(circuit,("Lj1","2","0",Lj)) 
    push!(circuit,("C2","2","0",Cj))

    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        Rleft => 50.0,
        Ipump => 1.0e-8,
    )

    return hblinsolve(2*pi*(4.5:0.1:5.0)*1e9,circuit,circuitdefs,Nmodes=1)
end

function warmupvvn()

    @variables Ipump Rleft Cc Lj Cj w L1
    # define the circuit components
    circuit = Array{Tuple{String,String,String,Num},1}(undef,0)

    # port on the left side
    push!(circuit,("P1","1","0",1))
    push!(circuit,("I1","1","0",Ipump))
    # push!(circuit,("R1","1","0",Rleft*w/(w+1)))
    push!(circuit,("R1","1","0",Rleft))

    push!(circuit,("C1","1","2",Cc)) 
    push!(circuit,("Lj1","2","0",Lj)) 
    push!(circuit,("C2","2","0",Cj))

    circuitdefs = Dict(
        Lj =>1000.0e-12,
        Cc => 100.0e-15,
        Cj => 1000.0e-15,
        Rleft => 50.0,
        Ipump => 1.0e-8,
    )

    psc = parsesortcircuit(circuit,sorting=:number)


    return valuevectortonumber(psc.valuevector,circuitdefs)
end


export @syms, hbnlsolve, hblinsolve, hbsolve, hbsolve2, hbnlsolve2, parsecircuit,
    parsesortcircuit, calccircuitgraph, symbolicmatrices, numericmatrices,
    LjtoIc, IctoLj, @variables, @register_symbolic, Num, Symbolics


# the below precompile directives are to help the compiler perform type inference
# during the precompilation stage (when the package is installed) instead of
# when it is loaded or the functions are run. this is a helpful guide:
# https://timholy.github.io/SnoopCompile.jl/stable/snoopi_deep_analysis/#inferrability
# and the basic commands to look at the inference triggers
# julia> using SnoopCompile, JosephsonCircuits
# julia> tinf = @snoopi_deep JosephsonCircuits.warmupsyms();
# julia> itrigs = inference_triggers(tinf)

@precompile_all_calls begin
    warmup()
    warmupsyms()
    warmupsyms2()
end

#end module
end
